function out = solveSDP(D,algParam)

    n = size(D,1);
    min_n_curr = 5;
    nJoin = algParam.nMergeCher;
    
    s = zeros(2*n-2,1);
    t = zeros(2*n-2,1);
    nodeCurr = 1:n;
    leavesCurr = (1:n)';
    Dcurr = D;
    e = 1;
    n_curr = n;
    L_curr = algParam.maxDepth;
    i = 1;
    while i <= n-1
        i
        if n_curr > min_n_curr 
            if strcmp(algParam.maxDepthType,'lin')
                L_curr = ceil(n_curr/2);
            end
            if strcmp(algParam.maxDepthType,'log')
                L_curr = min(ceil(2*log2(n_curr)),ceil(n_curr/2));
            end

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
            for j = 2:length(B_opt)
                B = abs(B_opt{j});
                Y = Y_opt{j};
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
                if (nClass == 1) && (j >= ceil(log2(n_curr)))
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
    
            if strcmp(algParam.heurRuleCh,'sep')
                % [r,c] = find(nCher == max(nCher(:)));
                nCher = nCher + nCher';
                match = greedyMatching(nCher, nJoin, 'maximize');
            end

            if strcmp(algParam.heurRuleCh,'prof')
                % [r,c] = find(Twin == min(Twin(:))); 
                match = greedyMatching(Twin, nJoin, 'minimize');
            end

            match = sort(match,2);
            for me = 1:size(match,1)
                u = match(me,1);
                v = match(me,2);
                s(e) = n+i; s(e+1) = n+i;
                t(e) = nodeCurr(u); t(e+1) = nodeCurr(v);
                nodeCurr(u) = n+i;
                e = e+2;
                n_curr = n_curr-1;
                i = i+1;
                % G = digraph(s(1:(e-1)),t(1:(e-1)));
                % figure; plot(G,'Layout','layered');
                if i == n-1
                    break;
                end
            end
            nodeCurr(match(:,2)) = [];
            leavesCurr(match(:,2)) = [];
            Dcurr(match(:,1),:) = (Dcurr(match(:,1),:)+Dcurr(match(:,2),:))/2;
            Dcurr(:,match(:,1)) = (Dcurr(:,match(:,1))+Dcurr(:,match(:,2)))/2;
            Dcurr(:,match(:,2)) = []; 
            Dcurr(match(:,2),:) = [];
            Dcurr = Dcurr - diag(diag(Dcurr));
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
    

