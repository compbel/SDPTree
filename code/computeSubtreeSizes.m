function tree = computeSubtreeSizes(tree, node)
    if tree(node).left == 0 && tree(node).right == 0
        % Leaf node.
        tree(node).subtreeSize = 1;
    else
        % Recursively compute sizes for left and right children.
        tree = computeSubtreeSizes(tree, tree(node).left);
        tree = computeSubtreeSizes(tree, tree(node).right);
        leftSize  = tree(tree(node).left).subtreeSize;
        rightSize = tree(tree(node).right).subtreeSize;
        tree(node).subtreeSize = leftSize + rightSize;
    end
end