function G = tree2digraph(tree)
%Convert the tree structure to a digraph object.

    nNodes = numel(tree);
    sources = [];
    targets = [];
    for i = 1:nNodes
        if tree(i).left ~= 0
            sources(end+1) = tree(i).label;
            targets(end+1) = tree(tree(i).left).label;
        end
        if tree(i).right ~= 0
            sources(end+1) = tree(i).label;
            targets(end+1) = tree(tree(i).right).label;
        end
    end
    
    G = digraph(sources, targets);    
    nodeLabels = arrayfun(@(x) x.label, tree);
    G.Nodes = table(nodeLabels(:), 'VariableNames', {'NewLabel'});
end