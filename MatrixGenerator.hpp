#include <iostream>

using namespace std;

class MatrixGenerator {
  public:
    double** generate(int r, int c);
    MatrixGenerator();
};

MatrixGenerator::MatrixGenerator(void) {
  cout << "MatrixGenerator is created" << endl;
}

double** MatrixGenerator::generate(int r, int c) {
  double** m = new double* [r];
  for (int i = 0; i < r; i++) {
    m[i] = new double[c];
  }
  return m;
}