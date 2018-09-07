## classic MDS algorithm

In CMDSSolver.hpp, CMDSSolver class can be used to get MDS projection using classic MDS. hello.cpp is a demo. CMDSSolver depends on Eigen3.

```
CMDSSolver cMDSSolver(matrixData, nE); // MatrixXd matrixData(nE, nE), int nE
double** projection = cMDSSolver.getProjection();
```

## Failure

I tried to implement nonmetric MDS. However, I failed. Thus, the nMDSSolver class in CMDSSolver.hpp does not work.  
There are two places that may be problematic:
1. the implementation of Pool-adjacent violators algorithm;
2. the parameter settings of initial configuration and iteration step;
3. my understanding of isotonic regression.  
The corresponding code that cannot work has been moved into another branch, namely, failure.