#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <time.h>

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
  ifstream fin("data/iris.csv");
  string line;
  vector<Entry> entries;
  while (getline(fin, line)) {
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
  clock_t time1, time2, time3, time4;
  time1 = clock();
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

  time2 = clock();
  CMDSSolver cMDSSolver(matrixData, nE);
  double** projection = cMDSSolver.getProjection();

  time3 = clock();

  ofstream outFile;
	outFile.open("data/data.csv", ios::out); // 打开模式可省略
  for (int r = 0; r < nE; r++) {
    outFile << projection[r][0] << ',' << projection[r][1] << endl;
  }
	outFile.close();

  time4 = clock();
  cout << double(time2 - time1)/ CLOCKS_PER_SEC << "s" << endl;
  cout << double(time3 - time2)/ CLOCKS_PER_SEC << "s" << endl;
  cout << double(time4 - time3)/ CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
