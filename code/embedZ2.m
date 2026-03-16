



function tree = embedZ2(tree)
    nNodes = numel(tree);
    tree(1).labelz2 = 1;
    for i = 1:nNodes
        if tree(i).left + tree(i).right > 0
            tree(tree(i).left).labelz2 = 2*tree(i).labelz2;
            tree(tree(i).right).labelz2 = 2*tree(i).labelz2+1;
        end
    end
end