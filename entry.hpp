#include <iostream>
#include <cstring>
#include <cmath>
#include <vector>

using namespace std;

class Entry {
  public:
    Entry(vector<double> vAttributes, int nDimension);
    double getDistance(Entry entry);
    void show( void );

  private:
    vector<double> attrs;
    int nD;
};

Entry::Entry(vector<double> vAttributes, int nDimension) {
  cout << "Object has been created" << endl;
  nD = nDimension;
  attrs.assign(vAttributes.cbegin(), vAttributes.cend());
  // double* attributes = new double[vAttributes.size()];
  // int i = 0;
  // for (i = 0; i < n; i++) {
  //   attributes[i] = vAttributes.at(i);
  // }
  // memcpy(attrs, attributes, sizeof(attributes));
}

double Entry::getDistance(Entry entry) {
  double deltaSum = 0.0, d = 0.0;
  for (int i = 0; i < nD; i++) {
    deltaSum += (entry.attrs[i] - attrs[i]) * (entry.attrs[i] - attrs[i]);
  }
  d = sqrt(deltaSum);
  // cout << d << endl;
  return d;
}

void Entry::show() {
  for (int i = 0; i < attrs.size(); i++) {
    cout << attrs[i] << ' ';
  }
  cout << endl;
}