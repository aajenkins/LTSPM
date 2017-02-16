/**
* @Author: Jenkins Alec <alec>
* @Date:   2017-02-12T14:02:42-08:00
* @Project: LTSPM analysis
* @Last modified by:   alec
* @Last modified time: 2017-02-12T14:11:35-08:00
*/


#include <cmath>
#include <iostream>

using namespace std;

int main() {
  const double pi = 3.14159265358979323846;
  double phi1 = (3.0/5)*pi;
  double phi0 = 0;
  double d_phi = remainder((phi1-phi0),(2*pi));

  cout << "d_phi =  " << d_phi << endl;
}
