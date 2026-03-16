# SDPTree
A Semidefinite Programming Relaxation-based Method for Balanced Minimum Evolution Phylogeny Inference

## Running SDPTree as a MATLAB script 

### **Prerequisites** 
#### 1. MATLAB**  
#### 2. MOSEK (https://www.mosek.com/) 
#### 3. YALMIP (https://yalmip.github.io/) 
#### 4. FastME (http://www.atgc-montpellier.fr/fastme/)

To run SDPTree algorithm, call the function SDPtree with the required parameters.

```out = SDPtree(D,maxDepth,mergeType,fmeFolder)```

### Input: 
* ``D`` : dissimilarity matrix
* ``maxDepth`` : upper bound on the depth of the rooted tree (in the simplest case maxDepth = n-1, where n = size(D,1))
* ``mergeType`` : heuristics to identify leaves merged at agglomerative rounding step
    Possible values:
     - `'sep'` : select leaves unseparable at maximum number of layers
     - `'prof'` : select leaves whose distance profiles to the remaining leaves are most similar
* ``fmeFolder`` : location of fastME executable

### Output: 
* ``out`` : structure with fields
     - `out.T` : tree (as a digraph object)
     - `'out.F'` : value of the objective function

### Example
```
addpath(genpath('/home/skumsp/YALMIP'));  % path to YALMIP          
addpath('/home/skumsp/mosek/10.2/toolbox/r2017a','-end');  % path to MOSEK
fmeFolder = './fastme';  % folder with FastME 

n = 10; m = 100;
D = squareform(pdist(rand(n,m))); % Eucledian distance matrix for randomly generated points

maxDepth = n-1;
mergeType = 'sep';         
% mergeType = 'prof';      

out = SDPtree(D,maxDepth,mergeType,fmeFolder);
```
