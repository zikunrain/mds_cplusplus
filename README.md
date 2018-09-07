## classic MDS algorithm

In CMDSSolver.hpp, CMDSSolver class can be used to get MDS projection using classic MDS. hello.cpp is a demo. CMDSSolver depends on Eigen3.

```
CMDSSolver cMDSSolver(matrixData, nE); // MatrixXd matrixData(nE, nE), int nE
double** projection = cMDSSolver.getProjection();
```

## Failure

I tried to implement nonmetric MDS. However, I failed.
Thus, the nMDSSolver class does not work.
The corresponding code that cannot work has been moved into the branch of failure.