function tree = reorderChildren(tree, node)
% REORDERCHILDREN Ensures that at every internal node, the child with more
% descendants becomes the left child.
    if tree(node).left ~= 0 && tree(node).right ~= 0
        left  = tree(node).left;
        right = tree(node).right;
        leftSize  = tree(left).subtreeSize;
        rightSize = tree(right).subtreeSize;
        
        % Swap children if the right child's subtree is larger.
        if rightSize > leftSize
            tree(node).left = right;
            tree(node).right = left;
        end
        
        % Recursively reorder for the children.
        tree = reorderChildren(tree, tree(node).left);
        tree = reorderChildren(tree, tree(node).right);
    end
end