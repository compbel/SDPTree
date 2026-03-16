function [D,r] = leafDist(tree)
% treeUndir = graph(tree.Edges.EndNodes(:,1), tree.Edges.EndNodes(:,2));
treeUndir = graph(tree.Edges, tree.Nodes);
% figure
% plot(treeUndir,'Layout','layered');
leaves = find(degree(treeUndir)==1);
D_all = distances(treeUndir);
D = D_all(leaves,leaves);

r = [];
if isa(tree,'digraph')
    root = find(indegree(tree)==0);
    r = D_all(:,root);
    r = r(leaves);
end

