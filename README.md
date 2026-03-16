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
