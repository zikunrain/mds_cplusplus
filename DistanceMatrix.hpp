#include <iostream>

#ifndef MATRIXGENERATOR
#define MATRIXGENERATOR
#include "MatrixGenerator.hpp"
#endif

using namespace std;

class DMatrix {
  public:
    DMatrix(double** matri, int size);
    int getSize( void );
    double** getB( double** J );
    void setSquareM( void );
    // double** multiply( double** );

  private:
    double** m;
    int size;
};

void printFormat( double** m, int size) {
  cout << "print matrix in format" << endl;
  for (int r = 0; r < size; r++) {
    for (int c = 0; c < size; c++) {
      cout << m[r][c] << " ";
    }
    cout << endl;
  }
}

int DMatrix::getSize( void ) {
  return size;
}

DMatrix::DMatrix(double** matrix, int n) {
  cout << "size of matrix " << n << endl;
  size = n;
  m = matrix;
}

void DMatrix::setSquareM( void ) {
  for (int r = 0; r < size; r++) {
    for (int c = 0; c < size; c++) {
      m[r][c] = m[r][c] * m[r][c];
    }
  }
}

double** matrixMultiply( double** a, double** b, int size) {
  MatrixGenerator generator;
  double** result = generator.generate(size, size);

  for (int r = 0; r < size; r++) {
    for (int c = 0; c < size; c++) {
      // r, c => result1[r][c]
      double sum = 0.0;
      for (int index = 0; index < size; index++) {
        sum += a[r][index] * b[index][c];
      }
      result[r][c] = sum;
    }
  }
  return result;
}

void matrixScale(double** m, int size, double s) {
  for (int r = 0; r < size; r++) {
    for (int c = 0; c < size; c++) {
      m[r][c] = m[r][c] * s;
    }
  }
}

double** DMatrix::getB( double** J ) {
  double** result1 = matrixMultiply(J, m, size);
  double** result2 = matrixMultiply(result1, J, size);
  matrixScale(result2, size, -0.5);
  return result2;
}