#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "DistanceMatrix.hpp"
#include "entry.hpp"

#ifndef MATRIXGENERATOR
#define MATRIXGENERATOR
#include "MatrixGenerator.hpp"
#endif

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
// #include <eigen3/Eigen/Diagonal>

#define DBL_MIN -1.79769313486231570E+308

using namespace Eigen;
using namespace std;

template <class Type>
Type stringToNum(const string& str) {
  istringstream iss(str);
  Type num;
  iss >> num;
  return num;
}

vector<Entry> readcsv() {
  ifstream fin("iris.csv");
  string line;
  vector<Entry> entries;
  // Entry entry[];
  while (getline(fin, line)) {
    cout << line << endl;
    istringstream sin(line);
    vector<double> flower;
    int cAttrs = 0;
    string field;
    while (getline(sin, field, ',')) {
      cAttrs++;
      double data;
      data = stringToNum<double>(field);
      flower.push_back(data);
    }
    flower.pop_back();
    Entry entry(flower, cAttrs);
    entries.push_back(entry);
  }
  return entries;
}

double** generateJMatrix(int size) {
  MatrixGenerator generator;

  double** I = generator.generate(size, size);
  for (int r = 0; r < size; r++) {
    for (int c = 0; c < size; c++) {
      if (c == r) {
        I[r][c] = 1.0;
      } else {
        I[r][c] = 0.0;
      }
    }
  }

  double** ones = generator.generate(size, size);
  for (int r = 0; r < size; r++) {
    for (int c = 0; c < size; c++) {
      ones[r][c] = 1.0/size;
    }
  }
  
  double** J = generator.generate(size, size);
  for (int r = 0; r < size; r++) {
    for (int c = 0; c < size; c++) {
      J[r][c] = I[r][c] - ones[r][c];
    }
  }

  return J;
}


int main() {
  int nE;
  cout << "program started" << endl;
  vector<Entry> entries= readcsv();
  nE = entries.size();

  cout << nE << endl;
  
  // initialize matrix
  MatrixGenerator generator;
  double **matrixData = generator.generate(nE, nE);
  for (int r = 0; r < nE; r++) {
    for (int c = 0; c < nE; c++) {
      // the cell of r row and c column is the distance between entries[r] and entries[c]
      matrixData[r][c] = entries[r].getDistance(entries[c]);
    }
  }

  DMatrix matrix(matrixData, nE);
  matrix.setSquareM();
  double** JMatrix = generateJMatrix(nE);
  double** B = matrix.getB(JMatrix);
  

  MatrixXd B2(nE, nE);
  int i = 1;
  for (int r = 0; r < nE; r++) {
    for (int c = 0; c < nE; c++) {
      B2(r, c) = B[r][c];
    }
  }
  cout << "Here is a matrix, B:" << endl << B2 << endl << endl;
  EigenSolver<MatrixXd> es(B2);
  MatrixXd D = es.pseudoEigenvalueMatrix();
  MatrixXd V = es.pseudoEigenvectors();
  cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;

  int maxEigenValueIndex = 0, secondMaxEigenValueIndex = 0;
  double maxEigenValue = D(0, 0), secondMaxEigenValue = DBL_MIN;
  for (int i = 0; i < nE; i++) {
    if (D(i, i) > maxEigenValue) {
      secondMaxEigenValue = maxEigenValue;
      secondMaxEigenValueIndex = maxEigenValueIndex;
      maxEigenValue = D(i, i);
      maxEigenValueIndex = i;
    } else if (D(i, i) > secondMaxEigenValue) {
      secondMaxEigenValue = D(i, i);
      secondMaxEigenValueIndex = i;
    }
  }
  // cout << "The top-2 eigenvalues are:" << endl << maxEigenValue << ' ' << secondMaxEigenValue << endl;
  // cout << "Corresponding index:" << endl << maxEigenValueIndex << ' ' << secondMaxEigenValueIndex << endl;
  // cout << "The pseudo-eigenvector matrix V is:" << endl << V << endl;
  // cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;

  MatrixXd maxVector = V.col(maxEigenValueIndex);
  MatrixXd secondMaxVector = V.col(secondMaxEigenValueIndex);
  // cout << "max vector" << endl << maxVector << endl;
  // cout << "second max vector" << endl << secondMaxVector << endl << endl;


  MatrixXd eigenVectorMatrix(maxVector.rows(), maxVector.cols() + secondMaxVector.cols());
  eigenVectorMatrix << maxVector, secondMaxVector;
  cout << eigenVectorMatrix << endl;

  // Vector2f diag(maxEigenValue, secondMaxEigenValue);
  // // DiagonalMatrix<double, 2> diagMatrix(maxEigenValue, secondMaxEigenValue);
  // cout << diag.asDiagonal() << endl;
  // Vector3f v(1,2,3);
  // cout << v.asDiagonal();
  Eigen::Matrix< double, 2, 1> diagVector;
  diagVector << maxEigenValue, secondMaxEigenValue;
  Eigen::Matrix< double, 2, 2> diagM = diagVector.array().sqrt().matrix().asDiagonal();
  cout << diagM << endl;
  cout << eigenVectorMatrix * diagM << endl;


  ofstream outFile;
	outFile.open("data.csv", ios::out); // 打开模式可省略
  for (int r = 0; r < eigenVectorMatrix.rows(); r++) {
    outFile << eigenVectorMatrix(r, 0) << ',' << eigenVectorMatrix(r, 1) << endl;
  }
	// outFile << "name" << ',' << "age" << ',' << "hobby" << endl;
	// outFile << "Mike" << ',' << 18 << ',' << "paiting" << endl;
	// outFile << "Tom" << ',' << 25 << ',' << "football" << endl;
	// outFile << "Jack" << ',' << 21 << ',' << "music" << endl;
	outFile.close();


  return 0;
}
