function [Y_opt, B_opt, z_opt, a_opt,Zprod_opt, F, out] = solveSDPRelaxation(D, L, depthCounts, regul)
% maxcut_treeZ_sdp  SDP relaxation with UNKNOWN z (depths unknown), optional depth counts
% Rescaled version: uses Yt^{(k)} = 2^k * Y^{(k)} to improve conditioning.
%
% INPUT:
%   D           : n x n symmetric matrix
%   L           : max depth (depths are 1..L)
%   depthCounts : L x 1 vector (# leaves at depth d), or [] if not fixed
%   regul       : [] or struct with fields:
%                   regul.a   
%                   regul.lambda (scalar)
%                 Optional: regul.eps0
%
% OUTPUT:
%   Y_opt : cell, Y_opt{k+1} = Y^{(k)} for k=0..L
%   B_opt : cell, B_opt{k+1} = Y^{(k)} ./ Zprod
%   z_opt : n x 1
%   a_opt : n x L
%   F     : objective value
%   out   : YALMIP output

    n = size(D,1);
    if ~isempty(depthCounts), depthCounts = depthCounts(:); end
    depths = (1:L)';

    Dwork = D;

    % ---- beta
    beta = zeros(L+1,1);
    beta(1) = 2;
    if L >= 1, beta(2) = 2; end
    for k = 2:L
        beta(k+1) = 4^k - 4^(k-1);
    end

    %--------------------------------------------------------------
    % YALMIP variables
    %--------------------------------------------------------------
    a = sdpvar(n, L, 'full');      % depth assignment
    z = sdpvar(n, 1);              % z_i = 2^{-depth(i)} (relaxed)

    % Scaled matrices: Yt{k+1} = 2^k * Y^{(k)}, k=0..L
    Yt = cell(L+1,1);
    for k = 0:L
        Yt{k+1} = sdpvar(n,n,'symmetric');
    end

    % Global cap matrix Zprod ~ z z'
    Zprod = sdpvar(n,n,'symmetric');

    one = ones(n,1);
    Constraints = [];

    %--------------------------------------------------------------
    % a constraints
    %--------------------------------------------------------------
    Constraints = [Constraints, 0 <= a <= 1];
    Constraints = [Constraints, sum(a,2) == 1];
    if ~isempty(depthCounts)
        Constraints = [Constraints, sum(a,1)' == depthCounts];
    end

    % z = sum_d 2^{-d} a_{i,d}
    w = 2.^(-depths);
    Constraints = [Constraints, z == a*w];
    Constraints = [Constraints, sum(z) == 1];
    Constraints = [Constraints, z >= 2^(-L), z <= 1/2];

    % Per-leaf square proxy
    w2 = 2.^(-2*depths);

    %--------------------------------------------------------------
    % Level constraints using Yt
    %--------------------------------------------------------------
    for k = 0:L
        Ytk = Yt{k+1};
        activeDepths = max(1,k):L;

        wAct  = zeros(L,1);
        w2Act = zeros(L,1);
        if ~isempty(activeDepths)
            wAct(activeDepths)  = 2.^(-activeDepths)';
            w2Act(activeDepths) = 2.^(-2*activeDepths)';
        end
        zAct_k = a*wAct;      % n×1
        sAct_k = a*w2Act;     % n×1

        % PSD + entrywise nonnegativity
        Constraints = [Constraints, Ytk >= 0];
        Constraints = [Constraints, Ytk(:) >= 0];

        % diag(Y^{(k)}) = sAct_k  ==>  diag(Yt) = 2^k * sAct_k
        Constraints = [Constraints, diag(Ytk) == 2^k * sAct_k];

        % balance: Y^{(k)}*1 = 2^{-k} zAct_k  ==>  Yt*1 = zAct_k
        Constraints = [Constraints, Ytk*one == zAct_k];

        % 2x2 PSD Cauchy–Schwarz tightening
        mask = zeros(L,1);
        mask(max(1,k):L) = 1;
        qAct_k = a * mask;

        for i = 1:n
            Constraints = [Constraints, [qAct_k(i)  zAct_k(i); ...
                                         zAct_k(i)  sAct_k(i)] >= 0];
        end
    end

    % Nestedness: Y^{(k+1)} <= Y^{(k)}  ==>  Yt^{(k+1)} <= 2 * Yt^{(k)}
    for k = 0:(L-1)
        Constraints = [Constraints, Yt{k+2}(:) <= 2 * Yt{k+1}(:)];
    end

    % Root tightening: [1 z'; z Y^{(0)}] PSD, and Y^{(0)} = Yt^{(0)}
    Constraints = [Constraints, [1 z'; z Yt{1}] >= 0];

    %--------------------------------------------------------------
    % Zprod constraints
    %--------------------------------------------------------------
    Constraints = [Constraints, Zprod >= 0];
    Constraints = [Constraints, Zprod(:) >= 0];
    Constraints = [Constraints, Zprod*one == z];
    Constraints = [Constraints, diag(Zprod) == a*w2];

    % Shor lift
    Constraints = [Constraints, [1 z'; z Zprod] >= 0];

    % Cap off-diagonal: Y^{(k)}(off) <= Zprod(off)  ==>  Yt^{(k)}(off) <= 2^k * Zprod(off)
    Ioff = ~eye(n);
    for k = 0:L
        Constraints = [Constraints, Yt{k+1}(Ioff) <= (2^k) * Zprod(Ioff)];
    end

    % Constraints = [Constraints, Zprod <= z*ones(1,n)];     % Z_ij <= z_i
    % Constraints = [Constraints, Zprod <= ones(n,1)*z'];    % Z_ij <= z_j

    %--------------------------------------------------------------
    % McCormick envelope for Zprod_ij with bounds ell <= z_i <= u
    %--------------------------------------------------------------
    if ~isempty(depthCounts)
        first = find(depthCounts>0, 1, 'first');
        last  = find(depthCounts>0, 1, 'last');
        ell = 2^(-last);
        u   = 2^(-first);
    else
        ell = 2^(-L);
        u   = 1/2;
    end

    I = triu(true(n),1);

    Zi = Zprod(I);
    zi = kron(ones(n,1), z); zi = zi(I(:));
    zj = kron(z, ones(n,1)); zj = zj(I(:));

    Constraints = [Constraints, Zi <= u*zi + ell*zj - u*ell];
    Constraints = [Constraints, Zi <= ell*zi + u*zj - u*ell];
    Constraints = [Constraints, Zi >= ell*zi + ell*zj - ell^2];
    Constraints = [Constraints, Zi >= u*zi + u*zj - u^2];

    %--------------------------------------------------------------
    % t-based tree-feasibility inequalities (optional tightening)
    %--------------------------------------------------------------
    t = sum(a,1)';  % L×1
    for kk = 1:L
        lhs = 0;
        for d = 1:kk
            lhs = lhs + 2^(kk-d) * t(d);
        end
        Constraints = [Constraints, lhs <= 2^kk];
    end

    for kk = 1:L
        if n > 2^kk
            Constraints = [Constraints, sum(t(1:kk)) <= 2^kk - 1];
        end
    end

% %-------------------------------------------------------------
    % Taxa form singletons after they become leaves
    Ioff = ~eye(n);

    for d = 1:L
        Ytd = Yt{d+1};          % level k=d
        M   = 2^(-d);

        ai = a(:,d);            % n×1
        Aj = repmat(ai', n, 1); % row-wise a_j
        Ai = repmat(ai , 1, n); % col-wise a_i

        Constraints = [Constraints, Ytd(Ioff) <= M * (1 - Ai(Ioff))];
        Constraints = [Constraints, Ytd(Ioff) <= M * (1 - Aj(Ioff))];
    end


    % % ==============================================================
    % % Quasi-transitivity for Y
    % for k = 0:L
    %     Yk = Yt{k+1};
    %     s  = 2^k;
    % 
    %     alpha_k = min(1, ell * 2^k);  
    % 
    %     for ell3 = 1:n
    %         for i = 1:n
    %             if i==ell3, continue; end
    %             for j = 1:n
    %                 if j==i || j==ell3, continue; end
    %                 Constraints = [Constraints, ...
    %                     Yk(i,j) >= alpha_k*(Yk(i,ell3) + Yk(ell3,j)) - s*Zprod(i,j)];
    %             end
    %         end
    %     end
    % end

    % % ---- Heuristic transitivity (vectorized)
    % theta = 0.01;    
    % one  = ones(n,1);
    % Ioff = ~eye(n);    
    % mask = cell(n,1);
    % for k = 1:n
    %     mk = Ioff;
    %     mk(k,:) = false;
    %     mk(:,k) = false;
    %     mask{k} = mk;
    % end    
    % for ell = 0:L
    %     Y = Yt{ell+1};   
    %     for k = 1:n
    %         RHS1 = theta*( Y(:,k)*one' );
    %         RHS2 = theta*( one*Y(k,:) );
    %         Constraints = [Constraints, Y(mask{k}) >= RHS1(mask{k})];
    %         Constraints = [Constraints, Y(mask{k}) >= RHS2(mask{k})];
    %     end
    % end


    %--------------------------------------------------------------
    % Objective: trace(D*Y^{(k)}) = trace(D*Yt^{(k)}) / 2^k
    %--------------------------------------------------------------
    Objective = 0;
    for k = 0:L
        Objective = Objective + (beta(k+1)/2^k) * trace(Dwork * Yt{k+1});
    end

    %--------------------------------------------------------------
    % Optional linearized entropy regularization (integrality heuristic)
    %--------------------------------------------------------------
    if ~isempty(regul)
        if isfield(regul,'eps0') && ~isempty(regul.eps0)
            eps0 = regul.eps0;
        else
            eps0 = 1e-6;   % stable default
        end

        a_ref = max(regul.a, eps0);   
        W = -(1 + log(a_ref));            

        % row-centering (per leaf) is the correct stabilization
        W = W - mean(W,2);

        % optional clipping
        W = max(min(W, 10), -10);

        Objective = Objective - regul.lambda * sum(sum(W .* a));
    end

    %--------------------------------------------------------------
    % Solve
    %--------------------------------------------------------------
    opts = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);
    out = optimize(Constraints, Objective, opts);

    if out.problem ~= 0
        warning('YALMIP/MOSEK reported problem code %d: %s', out.problem, out.info);
    end

    F = value(Objective);

    %--------------------------------------------------------------
    % Extract solution (unscale Yt -> Y)
    %--------------------------------------------------------------
    z_opt = value(z);
    a_opt = value(a);

    Y_opt = cell(L+1,1);
    B_opt = cell(L+1,1);
    Zprod_opt = cell(L+1,1);

    % Zden = max(value(Zprod), 1e-12);
    Zdiag = diag(value(Zprod));             
    Zden  = sqrt(Zdiag * Zdiag');           
    for k = 0:L
        Y_opt{k+1} = value(Yt{k+1}) / 2^k;
        B_opt{k+1} = Y_opt{k+1} ./ Zden;
        Zprod_opt{k+1} = value(Zprod);
    end
end