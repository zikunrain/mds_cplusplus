#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
// #include <random>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <time.h>

#define DBL_MIN -1.79769313486231570E+308

using namespace Eigen;
using namespace std;

typedef Matrix<double, Dynamic, Dynamic> MatrixDbDY;

class CMDSSolver {
  public:
    CMDSSolver(MatrixXd matrixData, int nE);
    double** getProjection( void );

  private:
    double** projection;
};

CMDSSolver::CMDSSolver(MatrixXd matrixData, int nE) {
  cout << "CMDSSolver has been created" << endl;

  MatrixDbDY matrixDataSquare = matrixData.cwiseProduct(matrixData);
  MatrixDbDY ones =  MatrixXd::Ones(nE, nE) / nE;
  MatrixDbDY identity = MatrixXd::Identity(nE, nE);
  MatrixDbDY J = identity - ones;
  MatrixDbDY B = (J / (-2.0)) * matrixDataSquare * J; 
  

  EigenSolver<MatrixXd> es(B);
  MatrixXd D = es.pseudoEigenvalueMatrix();
  MatrixXd V = es.pseudoEigenvectors();


  int maxEigenValueIndex = 0, secondMaxEigenValueIndex = 0;
  double maxEigenValue = DBL_MIN, secondMaxEigenValue = DBL_MIN;
  for (int i = 0; i < nE; i++) {
    if (D(i, i) > maxEigenValue) {
      secondMaxEigenValue = maxEigenValue;
      secondMaxEigenValueIndex = maxEigenValueIndex;
      maxEigenValue = D(i, i);
      maxEigenValueIndex = i;
      // cout << D(i, i) << endl;
    } else if (D(i, i) > secondMaxEigenValue) {
      secondMaxEigenValue = D(i, i);
      secondMaxEigenValueIndex = i;
    }
  }
  // cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;

  MatrixXd maxVector = V.col(maxEigenValueIndex);
  MatrixXd secondMaxVector = V.col(secondMaxEigenValueIndex);
  MatrixXd eigenVectorMatrix(maxVector.rows(), maxVector.cols() + secondMaxVector.cols());
  eigenVectorMatrix << maxVector, secondMaxVector;

  cout << maxEigenValueIndex << endl;
  cout << secondMaxEigenValueIndex << endl;
  // cout << D << endl;

  Eigen::Matrix< double, 2, 1> diagVector;
  diagVector << maxEigenValue, secondMaxEigenValue;
  Eigen::Matrix< double, 2, 2> diagM = diagVector.array().sqrt().matrix().asDiagonal();

  MatrixXd result = eigenVectorMatrix * diagM;
  projection = new double* [nE];
  for (int r = 0; r < nE; r++) {
    projection[r] = new double [2];
    projection[r][0] = result(r, 0);
    projection[r][1] = result(r, 1);
  }
}

double** CMDSSolver::getProjection() {
  return projection;
}


class NMDSSolver {
  public:
    NMDSSolver(MatrixXd matrixData, int nE);
    double** getProjection( void );
    MatrixXd getTempMatrix( void );
    bool hasKey(map <double, string> mymap, double key);
  
  private:
    double** projection;
    MatrixXd theltaMatrix;
};

NMDSSolver::NMDSSolver(MatrixXd matrixData, int nE) {
  srand(time(0));
  cout << "NMDSSolver has been created" << endl;
  // cout << "matrixData: " << endl << matrixData << endl;
  map <double, string> dbMapString;
  for (int r = 1; r < nE; r++) {
    for (int c = 0; c < r; c++) {
      double value = matrixData(r,c);
      string rowStr = to_string(r);
      string colStr = to_string(c);
      ostringstream rowColStr;
      rowColStr << rowStr << ',' << colStr;
      if (hasKey(dbMapString, value)) {
        string rowColStrPrev = dbMapString[value];
        ostringstream rowColStrCur;
        rowColStrCur << rowColStrPrev << ";" << rowColStr.str();
        // cout << "duplicate key: " << rowColStrCur.str() << endl;
        dbMapString[value] = rowColStrCur.str();
      } else {
        dbMapString[value] = rowColStr.str();
      }
    }
  }

  MatrixDbDY nMatrixData = MatrixXd::Zero(nE, nE);
  vector<pair<int , int>> orderedRowCol;

  map <double, string>::iterator iter;
  iter = dbMapString.begin();
  int order = 1;
  while(iter != dbMapString.end()) {
    // cout << "checking dbmapstring: " << iter->first << '-' << iter->second << endl;
    string rowColStrs = iter->second;
    istringstream in(rowColStrs);
    string rowColStr;
    int curOrder = order;
    while (getline(in, rowColStr, ';')) {
      // string to char
      char rowColChar[rowColStr.length()];
      int ichar, istr;
      int row, col;
      for (ichar = 0, istr = 0; istr < rowColStr.length(); istr++, ichar++) {
        rowColChar[ichar] = rowColStr[istr];
        if (rowColChar[ichar] == ',') { // comma
          rowColChar[ichar] = '\0';
          row = atoi(rowColChar);
          ichar = -1;
        }
      }
      rowColChar[ichar] = '\0';
      col = atoi(rowColChar);

      nMatrixData(row, col) = curOrder;
      nMatrixData(col, row) = curOrder;

      pair<int, int> p(row, col);
      orderedRowCol.push_back(p);
      order++;
    }
    iter++;
  }
  theltaMatrix = nMatrixData;

  // assign random configuration
  MatrixDbDY conf(nE, 2);
  for (int r = 0; r < nE; r++) {
    conf(r, 0) = rand()%100/(double)101;
    conf(r, 1) = rand()%100/(double)101;
  }

  
  MatrixDbDY D = MatrixXd::Zero(nE, nE);
  MatrixDbDY DHat = MatrixXd::Zero(nE, nE);
  double stress;

  for (int t = 0; t < 100; t++) {
    // calculate d
    double drc;
    double sumDrc = 0.0, sumDrcSquare = 0.0;
    for (int r = 1; r < nE; r++) {
      for (int c = 0; c < r; c++) {
        drc = sqrt((conf(r, 0) - conf(c, 0)) * (conf(r, 0) - conf(c, 0)) + (conf(r, 1) - conf(c, 1)) * (conf(r, 1) - conf(c, 1)));
        double drcValue;
        if (drc < 0.0001) {
          drcValue = 0.01;
        } else {
          drcValue = drc;
        }
        D(r, c) = drcValue;
        sumDrc += drcValue;
        sumDrcSquare += drcValue * drcValue;
      }
    }

    double avgDrc = sumDrc * 2 / (nE * (nE - 1));

    // calculate dhat
    for (int i = 1; i < orderedRowCol.size(); i++) {
      int row = orderedRowCol[i].first;
      int col = orderedRowCol[i].second;
      int prevRow = orderedRowCol[i - 1].first;
      int prevCol = orderedRowCol[i - 1].second;

      if (D(prevRow, prevCol) > D(row, col)) {
        double midean = (D(prevRow, prevCol) + D(row, col)) / 2.0;
        DHat(prevRow, prevCol) = midean;
        DHat(row, col) = midean;
      } else {
        DHat(row, col) = D(row, col);
      }
    }

    double sumDiffSquare = 0.0;
    for (int r = 1; r < nE; r++) {
      for (int c = 0; c < r; c++) {
        sumDiffSquare += (DHat(r, c) - D(r, c)) * (DHat(r, c) - D(r, c));
      }
    }
    
    stress = sqrt(sumDiffSquare / sumDrcSquare);
    cout << " stress : " << stress << endl;

    // if stress is larger than ...
    for (int r = 0; r < nE; r++) {
      double item0 = 0.0;
      for (int c = 0; c < nE; c++) {
        if (c != r) {
          if (c > r) {
            item0 += (1 - (DHat(c, r)/D(c, r))) * (conf(c, 0) - conf(r, 0));
          } else {
            item0 += (1 - (DHat(r, c)/D(r, c))) * (conf(c, 0) - conf(r, 0));
          }
        }
      }
      conf(r, 0) = conf(r, 0) + item0 * 0.4 / (nE - 1);

      double item1 = 0.0;
      for (int c = 0; c < nE; c++) {
        if (c != r) {
          if (c > r) {
            item1 += (1 - (DHat(c, r)/D(c, r))) * (conf(c, 1) - conf(r, 1));
          } else {
            item1 += (1 - (DHat(r, c)/D(r, c))) * (conf(c, 1) - conf(r, 1));
          }
        }
      }
      conf(r, 1) = conf(r, 1) + item1 * 0.4 / (nE - 1);
    }
  }

  cout << conf << endl << endl;
}

MatrixXd NMDSSolver::getTempMatrix() {
  return theltaMatrix;
}


bool NMDSSolver::hasKey(map <double, string> mymap, double key) {
  map <double, string>::iterator it = mymap.find(key);  
  if(it == mymap.end()){    
    return false;    
  } else {    
    return true;    
  }    
}
