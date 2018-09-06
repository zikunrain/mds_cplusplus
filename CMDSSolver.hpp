#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

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
  
  private:
    double** projection;
}

NMDSSolver::NMDSSolver(MatrixXd matrixData, int nE) {
  cout << "NMDSSolver has been created" << endl;
  map <double, string> dbMapString;
  for (int r = 0; r < nE; r++) {
    for (int c = 0; c <= r; c++) {
      double value = matrixData(r,c);
      string rowStr = to_string(r);
      string colStr = to_string(c);
      ostringstream rowColStr;
      rowColStr << rowStr << ',' << colStr;
      if (hasKey(dbMapString, value)) {
        string rowColStrPrev = dbMapString[value];
        ostringstream rowColStrCur;
        rowColStrCur << rowColStrPrev << ";" << rowColStr;
        dbMapString.insert(pair<double, string>(value, rowColStrCur.str()));
      } else {
        dbMapString.insert(pair<double, string>(value, rowColStr.str()));
      }
    }
  }

  MatrixDbDY nMatrixData = MatrixXd::Zero(nE, nE);

  map <double, string>::iterator iter;
  iter = dbMapString.begin();
  int order = 1;
  while(iter != dbMapString.end()) {
    string rowColStrs = iter->second;
    cout << iter->first << endl;
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
          cout << rowColChar << endl;
          ichar = -1;
        }
      }
      rowColChar[ichar] = '\0';
      cout << rowColChar << endl;
      col = atoi(rowColChar);
      cout << row << " "<< col << " order: " << order << endl;

      nMatrixData(row, col) = order;
      nMatrixData(col, row) = order;
      order++;
    }

    // cout << iter->first << " : " << iter->second << endl;
    iter++;
  }


}

bool hasKey(map <double, string> mymap, double key) {
  map <double, string>::iterator it = mymap.find(key);  
  if(it == mymap.end()){    
    return false;    
  } else {    
    return true;    
  }    
}
