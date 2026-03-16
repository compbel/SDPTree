clear;

addpath(genpath('/home/skumsp/YALMIP'));  % path to YALMIP          
addpath('/home/skumsp/mosek/10.2/toolbox/r2017a','-end');  % path to MOSEK
fmeFolder = './fastme';  % folder with FastME 

n = 10; m = 100;
D = squareform(pdist(rand(n,m))); % Eucledian distance matrix for randomly generated points

maxDepth = n-1;            % maximal depth of the tree

% heuristics to identify leaves merged at agglomerative rounding step 
% 'sep': select leaves unseparable at maximum number of layers
% 'prof': select leaves whose distance profiles to the remaining leaves are most similar

mergeType = 'sep';         
% mergeType = 'prof';      

out = SDPtree(D,maxDepth,mergeType,fmeFolder);

% out: solution structure. out.T - tree, out.F - value of the objective function