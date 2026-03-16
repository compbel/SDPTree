function tree = reconcileTreeLabels(tree)
    nNodes = numel(tree);
    
    [~, sortedIdx] = sort([tree.label]);  

    mapping = zeros(1, nNodes);
    mapping(sortedIdx) = 1:nNodes;
    
    tree = tree(sortedIdx);
    
    for i = 1:nNodes
        if tree(i).parent > 0
            tree(i).parent = mapping(tree(i).parent);
        end
        if tree(i).left > 0
            tree(i).left = mapping(tree(i).left);
        end
        if tree(i).right > 0
            tree(i).right = mapping(tree(i).right);
        end
    end
end