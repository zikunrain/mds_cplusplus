#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
using namespace Eigen;
using namespace std;



void Eig() {
  // typedef Matrix<int, 2, 2> MatrixMDS;
  // MatrixMDS M;
  // M << 1, 2 ,3, 6;
  // cout << M << endl;
  
  // Matrix3d A;
  // A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  
  MatrixXd A(3, 3);
  int i = 1;
  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      A(r, c) = i;
      i++;
    }
  }
  cout << "Here is a 3x3 matrix, A:" << endl << A << endl << endl;
  EigenSolver<Matrix3d> es(A);

  Matrix3d D = es.pseudoEigenvalueMatrix();
  Matrix3d V = es.pseudoEigenvectors();
  cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
  cout << "The pseudo-eigenvector matrix V is:" << endl << V << endl;
  cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;
}
int main() {
  Eig();
  return 0;
}