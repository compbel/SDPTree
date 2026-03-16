function leafSeq = simulateEvolTree(tree,L,p01,p10)
    % p01 % P(0->1) per edge
    % p10 % P(1->0) per edge
    Lmax = 5*L; % allow auto-increase to guarantee uniqueness
    % ---------------------
    
    N = numnodes(tree);
    
    nLeaves = sum(outdegree(tree)==0); 
    ord = toposort(tree);
    
    while true
        seq = false(N,L); 
        for u = ord(:).'
            for v = successors(tree,u).'
                pu = seq(u,:); r = rand(1,L); ch = pu;
                z = ~pu; o = pu;
                ch(z) = r(z) < p01;       % parent 0 -> 1 with p01
                ch(o) = ~(r(o) < p10);    % parent 1 -> 0 with p10
                seq(v,:) = ch;
            end
        end
        leafSeq = double(seq(1:nLeaves,:));
        if size(unique(leafSeq,'rows'),1) == nLeaves 
            break; 
        end
        L = L + 1;  
        if L > Lmax 
            error('Could not get unique leaves by Lmax=%d. Increase Lmax or change p01/p10.', Lmax); 
        end
    end
