#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "entry.hpp"
#include "CMDSSolver.hpp"

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
    // cout << line << endl;
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

int main() {
  int nE;
  // cout << "program started" << endl;
  vector<Entry> entries= readcsv();
  nE = entries.size();

  typedef Matrix<double, Dynamic, Dynamic> MatrixDbDY;
  MatrixXd matrixData(nE, nE);
  for (int r = 0; r < nE; r++) {
    for (int c = 0; c < nE; c++) {
      matrixData(r, c) = entries[r].getDistance(entries[c]);
    }
  }

  CMDSSolver cMDSSolver(matrixData, nE);
  double** projection = cMDSSolver.getProjection();
  
  // MatrixDbDY matrixDataSquare = matrixData.cwiseProduct(matrixData);
  // MatrixDbDY ones =  MatrixXd::Ones(nE, nE) / nE;
  // MatrixDbDY identity = MatrixXd::Identity(nE, nE);
  // MatrixDbDY J = identity - ones;
  // MatrixDbDY B = (J / (-2.0)) * matrixDataSquare * J; 
  

  // EigenSolver<MatrixXd> es(B);
  // MatrixXd D = es.pseudoEigenvalueMatrix();
  // MatrixXd V = es.pseudoEigenvectors();
  // cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
  // cout << "The pseudo-eigenvectors matrix D is:" << endl << V << endl;

  // int maxEigenValueIndex = 0, secondMaxEigenValueIndex = 0;
  // double maxEigenValue = D(0, 0), secondMaxEigenValue = DBL_MIN;
  // for (int i = 0; i < nE; i++) {
  //   if (D(i, i) > maxEigenValue) {
  //     secondMaxEigenValue = maxEigenValue;
  //     secondMaxEigenValueIndex = maxEigenValueIndex;
  //     maxEigenValue = D(i, i);
  //     maxEigenValueIndex = i;
  //   } else if (D(i, i) > secondMaxEigenValue) {
  //     secondMaxEigenValue = D(i, i);
  //     secondMaxEigenValueIndex = i;
  //   }
  // }
  // // cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;

  // MatrixXd maxVector = V.col(maxEigenValueIndex);
  // MatrixXd secondMaxVector = V.col(secondMaxEigenValueIndex);
  // // cout << "max vector" << endl << maxVector << endl;
  // // cout << "second max vector" << endl << secondMaxVector << endl << endl;


  // MatrixXd eigenVectorMatrix(maxVector.rows(), maxVector.cols() + secondMaxVector.cols());
  // eigenVectorMatrix << maxVector, secondMaxVector;
  // cout << eigenVectorMatrix << endl;

  // Eigen::Matrix< double, 2, 1> diagVector;
  // diagVector << maxEigenValue, secondMaxEigenValue;
  // Eigen::Matrix< double, 2, 2> diagM = diagVector.array().sqrt().matrix().asDiagonal();
  // cout << diagM << endl;
  // cout << eigenVectorMatrix * diagM << endl;


  ofstream outFile;
	outFile.open("data.csv", ios::out); // 打开模式可省略
  for (int r = 0; r < nE; r++) {
    outFile << projection[r][0] << ',' << projection[r][1] << endl;
  }
	outFile.close();


  return 0;
}
