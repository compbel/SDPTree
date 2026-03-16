function out = solveSDP(D,algParam)

    n = size(D,1);
    min_n_curr = 5;
    
    s = zeros(2*n-2,1);
    t = zeros(2*n-2,1);
    nodeCurr = 1:n;
    leavesCurr = (1:n)';
    Dcurr = D;
    e = 1;
    n_curr = n;
    L_curr = algParam.maxDepth;
    for i = 1:(n-1)
        i
        if n_curr > min_n_curr 
            L_curr = ceil(n_curr/2);

            beta = zeros(L_curr+1,1);
            beta(1) = 2;
            if L_curr >= 1 
                beta(2) = 2; 
            end
            for k = 2:L_curr
                beta(k+1) = 4^k - 4^(k-1);
            end
    
            [Y_opt, B_opt, z_opt, a_opt, Zprod_opt, F_lb, out] = solveSDPRelaxation(Dcurr, L_curr,[],[]);
            nCher = zeros(n_curr,n_curr);
            Btot = zeros(n_curr,n_curr);
            K = zeros(n_curr,n_curr);
            AMtot = zeros(n_curr,n_curr);
            for j = 2:length(B_opt)
                B = abs(B_opt{j});
                Y = Y_opt{j};
                % Z_prod = Zprod_opt{j};
                Z_prod = Zprod_opt;
                Bnorm = B;
                scale = max(Bnorm(~eye(n_curr)),[],'all');
                Bnorm(~eye(size(Bnorm,1))) = Bnorm(~eye(n_curr)) / scale;
                Btot = Btot + (B - diag(diag(B)));
                K = K + beta(j)*Y;
                try
                     [thr,nClass,prob] = isBimodMatr(B,0.00001,4);
                catch ME
                    crashDir = 'crash_dumps';
                    ts = getCurrentTask(); wid = 0; 
                    if ~isempty(ts), wid = ts.ID; end
                    fn = fullfile(crashDir, sprintf('crash_w%d_%s.mat', wid, datestr(now,'yyyymmdd_HHMMSSFFF')));
                    save(fn); 
                    save(fn,'ME','-append'); 
                    rethrow(ME);
                end
                if nClass == 1
                    break;
                end
                 
                candCher = nchoosek(1:n_curr,2);
        
                idx = sub2ind(size(B), candCher(:,1), candCher(:,2)); 
                [mval,mind] = max(B(idx));
                r = candCher(mind,1); c = candCher(mind,2);  
                nCher(r,c) = nCher(r,c) + 1;
            end

            d = diag(K);
            Twin = pdist2(K,K,'cityblock') - abs(d*ones(1,n_curr) - K') ...
            - abs(K - ones(n_curr,1)*d');
            Twin(1:n_curr+1:end) = Inf;
        
            nCher(1:n_curr+1:end) = -Inf;
            Badj = (AMtot>0).*Btot;
            G1 = graph(Badj);
    
            if strcmp(algParam.heurRuleCh,'sep')
                [v1,v2] = findedge(G1); 
                idx = find(degree(G1,v1)==1 & degree(G1,v2)==1);                 
                if ~isempty(idx) 
                    [~,k] = max(G1.Edges.Weight(idx)); 
                    r = v1(idx(k));
                    c = v2(idx(k));
                else
                    [r,c] = find(nCher == max(nCher(:)));                       
                end
            end

            if strcmp(algParam.heurRuleCh,'prof')
                [r,c] = find(Twin == min(Twin(:))); 
            end
            leavesCurrPrev = leavesCurr;       
            mergeCher;
        else
            [~, T] = runNJ(Dcurr,algParam);
            root = find(indegree(T)==0);
            ord = flip(dfsearch(T,root));
            nodeCurr = [nodeCurr zeros(1,n_curr-1)];
            int_num = max(nodeCurr)+1;
            for k = 1:length(ord)
                child = successors(T,ord(k));
                if isempty(child)
                    continue;
                end
                nodeCurr(ord(k)) = int_num;
                s(e) = nodeCurr(ord(k));  s(e+1) = nodeCurr(ord(k));
                t(e) = nodeCurr(child(1)); t(e+1) = nodeCurr(child(2));
                int_num = int_num+1;
                e = e+2;
            end
            break;
        end
    
    end
    
    T_SDP = digraph(s,t);
    T_SDP_unr = rooted2unrooted(T_SDP);
    
    F_SDP = calcObjMatrix(T_SDP,D,algParam.unrooted); 
    out.T = T_SDP;
    out.F = F_SDP;
    out.status = 1;
    
    function mergeCher
        u = min([r(1),c(1)]);
        v = max([r(1),c(1)]);
        s(e) = n+i; s(e+1) = n+i;
        t(e) = nodeCurr(u); t(e+1) = nodeCurr(v);
        nodeCurr(u) = n+i;
        nodeCurr(v) = [];
        leavesCurr(v) = [];
        Dcurr(u,:) = (Dcurr(u,:)+Dcurr(v,:))/2;
        Dcurr(:,u) = (Dcurr(:,u)+Dcurr(:,v))/2;
        Dcurr(:,v) = []; Dcurr(v,:) = [];
        Dcurr(u,u) = 0;
        e = e+2;
        n_curr = n_curr-1;
    end

end

