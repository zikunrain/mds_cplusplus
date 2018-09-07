## classic MDS algorithm

In CMDSSolver.hpp, CMDSSolver class can be used to get MDS projection using classic MDS. hello.cpp is a demo. CMDSSolver depends on Eigen3.
'''CMDSSolver cMDSSolver(matrixData, nE); // MatrixXd matrixData(nE, nE), int nE
  double** projection = cMDSSolver.getProjection();'''