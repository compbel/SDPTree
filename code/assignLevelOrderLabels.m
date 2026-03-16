

function tree = assignLevelOrderLabels(tree)
% ASSIGNLEVELORDERLABELS Assigns new labels to the nodes using a breadth-first 
% (level order) traversal. This ensures that within each level the labels increase 
% from left to right.
    label = 1;
    queue = 1;  % Start with the root node.
    
    while ~isempty(queue)
        current = queue(1);  % Dequeue the first element.
        queue(1) = [];       % Remove it from the queue.
        
        % Assign the next available label.
        tree(current).label = label;
        label = label + 1;
        
        % Enqueue the children (if any) in left-to-right order.
        if tree(current).left ~= 0
            queue(end+1) = tree(current).left;
        end
        if tree(current).right ~= 0
            queue(end+1) = tree(current).right;
        end
    end
end
