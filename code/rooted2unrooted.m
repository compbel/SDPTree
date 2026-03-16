function Tunr = rooted2unrooted(T)
    root = find(indegree(T)==0);
    ch = successors(T,root);
    if ismember('Weight', T.Edges.Properties.VariableNames)
        T = addedge(T,ch(1),ch(2),1);
    else
        T = addedge(T,ch(1),ch(2));
    end
    A = adjacency(T)>0;
    idx = true(size(A,1),1); 
    idx(root)=false;
    A = A(idx, idx);
    if ismember('Name', T.Nodes.Properties.VariableNames)
        names = T.Nodes.Name(idx);
        Tunr  = graph(A | A.', names); 
    else
        Tunr  = graph(A | A.'); 
    end
