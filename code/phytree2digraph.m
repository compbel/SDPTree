function [G, rootIdx] = phytree2digraph(tr,leafList)
% Convert a phytree to a digraph (parent -> child).

    nL  = get(tr,'NumLeaves');
    nB  = get(tr,'NumBranches');
    N   = nL + nB;

    ptr = get(tr,'Pointers');    
    dst = get(tr,'Distances');  

    parent = zeros(N,1);
    for b = 1:nB
        p  = nL + b;      
        c1 = ptr(b,1);
        c2 = ptr(b,2);
        parent(c1) = p;
        parent(c2) = p;
    end

    rootIdx = find(parent == 0, 1, 'first');

    children = (1:N).';
    keep     = parent ~= 0;
    s = parent(keep);
    t = children(keep);
    w = dst(t);               

    G = digraph(s, t, w, N);

    leafNames   = string(get(tr,'LeafNames'));
    branchNames = string(get(tr,'BranchNames'));
    G.Nodes.Name = [leafNames(:); branchNames(:)];
    types = strings(N,1); types(1:nL) = "leaf"; types(nL+1:end) = "internal";
    G.Nodes.Type = categorical(types);

    S = string(G.Nodes.Name);
    leaves = find(outdegree(G)==0);                             
    [tf,loc] = ismember(string(leafList), S(leaves));          
    idxL = leaves(loc(tf));                                     
    idxI = setdiff(1:numnodes(G), leaves, 'stable');               
    idx = [idxL; idxI'];
    G = reordernodes(G, idx);
end