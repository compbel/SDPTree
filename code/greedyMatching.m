function edges = greedyMatching(W, k, obj)
%GREEDYMATCHINGGENERAL Greedy matching in a general undirected graph.
% W: NxN weight matrix. 0, Inf, NaN indicate no edge.
% k: desired number of matched edges (optional). If empty, take as many as possible.
% maximize: true for maximum-weight, false for minimum-weight (default true).
%
% edges: m-by-2 array [u v] of chosen edges (1-based indices), m <= k

n = size(W,1);

% Candidate edges: i<j, finite, nonzero
mask = triu(true(n),1) & isfinite(W) & (W ~= 0);

[i,j] = find(mask);
w = W(mask);

% Sort edges by weight
if strcmp(obj,'maximize')
    [w, ord] = sort(w, 'descend');
else
    [w, ord] = sort(w, 'ascend');
end
i = i(ord); j = j(ord);

% Greedy select edges without shared vertices
used = false(n,1);
mmax = min(k, floor(n/2));
edges = zeros(mmax,2);

cnt = 0;
for t = 1:numel(w)
    u = i(t); v = j(t);
    if ~used(u) && ~used(v)
        cnt = cnt + 1;
        edges(cnt,:) = [u v];
        used(u) = true; used(v) = true;
        if cnt == mmax, break; end
    end
end

edges = edges(1:cnt,:);
end