#include <map>
#include <string>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;
typedef Matrix<double, Dynamic, Dynamic> MatrixDbDY;

bool hasKey(map <double, string> mymap, double key) {
  map <double, string>::iterator it = mymap.find(key);  
  if(it == mymap.end()){    
    return false;    
  } else {    
    return true;    
  }    
}

int main() {
  map <double, string> dbMapString;

  dbMapString.insert(pair<double, string>(2.6, "1,0"));
  dbMapString.insert(pair<double, string>(3.6, "2,1;3,0"));
  dbMapString.insert(pair<double, string>(4.6, "2,0"));
  dbMapString.insert(pair<double, string>(1.6, "3,1"));
  dbMapString.insert(pair<double, string>(0.6, "3,2"));

  int nE = 4;


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
  cout << nMatrixData << endl;


  return 0;
}


