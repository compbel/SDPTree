function [T,tree] = simulateCoalesTree(n)

AM = zeros(2*n-1,2*n-1);
roots = 1:n;
for t = (n+1):(2*n-1)
    coal_lin = randsample(roots, 2);
    AM(t,coal_lin(1)) = 1;
    AM(t,coal_lin(2)) = 1;
    roots(roots==coal_lin(1) | roots==coal_lin(2)) = [];
    roots(end+1) = t;
end

nNodes = size(AM,1);
perm = nNodes:-1:1;
AM = AM(perm,perm);

tree = repmat(struct('parent', 0, 'left', 0, 'right', 0), nNodes, 1);

for i = 1:nNodes
    p = find(AM(:, i), 1);
    if ~isempty(p)
        tree(i).parent = p;
    end
    children = find(AM(i, :));
    if numel(children) == 2
        tree(i).left = children(1);
        tree(i).right = children(2);
    end
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
leaves = leaves(randperm(n));
perm = [leaves; df(ind)];
T = reordernodes(T, perm);
