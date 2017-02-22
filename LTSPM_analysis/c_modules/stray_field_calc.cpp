/**
* @Author: Jenkins Alec <alec>
* @Date:   2017-02-10T16:37:10-08:00
* @Project: LTSPM analysis
* @Last modified by:   alec
* @Last modified time: 2017-02-16T15:37:27-08:00
*/

#include <iostream>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <stdlib.h>

using namespace std;

const int DLEN = 500;
const double res = 5.0e-7;
const double height = 40.0e-7;
const double pi = 3.14159265358979323846;
const double Ms = 1.0e3;
const double thickness = 1.0e-7;

double get_hz(double x, double y, double mx [DLEN][DLEN], double my [DLEN][DLEN], double mz [DLEN][DLEN]);

int main()
{
  std::ifstream mxs("/Users/alec/UCSB/cofeb_analysis_data/stray_field_test/mnrx_test.dat");
  std::istream_iterator<double> mxstart(mxs), mxend;
  std::vector<double> mxdatalist(mxstart, mxend);
  mxs.close();
  std::ifstream mys("/Users/alec/UCSB/cofeb_analysis_data/stray_field_test/mnry_test.dat");
  std::istream_iterator<double> mystart(mys), myend;
  std::vector<double> mydatalist(mystart, myend);
  mys.close();
  std::ifstream mzs("/Users/alec/UCSB/cofeb_analysis_data/stray_field_test/mnrz_test.dat");
  std::istream_iterator<double> mzstart(mzs), mzend;
  std::vector<double> mzdatalist(mzstart, mzend);
  mzs.close();

  int listsize = DLEN*DLEN;
  double mx [DLEN][DLEN];
  double my [DLEN][DLEN];
  double mz [DLEN][DLEN];
  double hz [DLEN][DLEN];

  for (int i=0; i<listsize; i++) {
    mx[(int) (i/DLEN)][i%DLEN] = mxdatalist[i];
    my[(int) (i/DLEN)][i%DLEN] = mydatalist[i];
    mz[(int) (i/DLEN)][i%DLEN] = mzdatalist[i];
  }
  for (int y=0; y<DLEN; y++) {
    for (int x=0; x<DLEN; x++) {
      hz[y][x] = get_hz(x, y, mx, my, mz);
    }
    std::cout << y << "\n";
  }

  ofstream hz_file;
  hz_file.open("/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/stray_field_test_c.txt");
  for (int y=0; y<DLEN; y++) {
    for (int x=0; x<DLEN-1; x++) {
      hz_file << hz[y][x] << ", ";
    }
    hz_file << hz[y][DLEN-1] << "\n";
  }
  hz_file.close();

}

double get_hz(double xp, double yp, double mx [DLEN][DLEN], double my [DLEN][DLEN], double mz [DLEN][DLEN]) {
  double hz = 0;
  double z = height/res;

  for (int y = 0; y<DLEN; y++) {
    for (int x = 0; x<DLEN; x++) {
      hz += (1/(4*pi)) * (pow(1/res,3)) * Ms * thickness * res * res
      * ( 3*z*(mx[y][x]*(xp-x) + my[y][x]*(yp-y) + mz[y][x]*z)/(pow(pow((xp-x),2) + pow((yp-y),2) + pow(z,2), 5/2))
        + mz[y][x]/(pow(pow((xp-x),2) + pow((yp-y),2) + pow(z,2), 3/2)) );
    }
  }

  return hz;
}
