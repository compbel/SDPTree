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
