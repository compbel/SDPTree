function [T,tree] = simulateYuleTree(nLeaves)
    tree(1).parent = 0;
    tree(1).left   = 0;
    tree(1).right  = 0;
    
    leaves = 1;  
    nextNode = 2;   
    
    while numel(leaves) < nLeaves
        idx  = randi(numel(leaves));
        leaf = leaves(idx);
        
        tree(leaf).left  = nextNode;
        tree(leaf).right = nextNode + 1;
        
        tree(nextNode).parent = leaf;
        tree(nextNode).left   = 0;
        tree(nextNode).right  = 0;
        
        tree(nextNode+1).parent = leaf;
        tree(nextNode+1).left   = 0;
        tree(nextNode+1).right  = 0;
        
        leaves(idx) = [];
        leaves = [leaves, nextNode, nextNode+1];
        
        nextNode = nextNode + 2;
    end

    tree = computeSubtreeSizes(tree, 1);
    
    tree = reorderChildren(tree, 1);
    
    tree = assignLevelOrderLabels(tree);
    tree = reconcileTreeLabels(tree);
    tree = embedZ2(tree);
    T = tree2digraph(tree); 

    root = find(indegree(T)==0);
    leaves = find(outdegree(T)==0);
    intern = find(outdegree(T)>0);
    df = flip(dfsearch(T,root));
    ind = ismember(df,intern);
    leaves = leaves(randperm(nLeaves));
    perm = [leaves; df(ind)];
    T = reordernodes(T, perm);
    
end


