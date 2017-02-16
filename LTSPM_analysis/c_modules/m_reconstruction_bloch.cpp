/**
* @Author: Jenkins Alec <alec>
* @Date:   2017-02-10T16:37:10-08:00
* @Project: LTSPM analysis
* @Last modified by:   alec
* @Last modified time: 2017-02-13T21:25:23-08:00
*/

#include <iostream>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <ctime>
#include <stdlib.h>

using namespace std;

const int DLEN = 100;
const double pi = 3.14159265358979323846;
const int max_num_steps = 5000;
const double cutoff_residual = 1.0e-5;
const double max_angle_initial = 0.05*pi/180;
const double max_angle_decay = 2500;
const int x_center = 47;
const int y_center = 59;
const double angle_offset = -0.5+(pi/2);
const double momentum = 0.95;

double lagrangian(double mz[DLEN][DLEN], double phi[DLEN][DLEN], int x, int y);
double dphi_lagrangian(double mz[DLEN][DLEN], double phi[DLEN][DLEN], int x, int y);
double dgradphi_lagrangian_x(double mz[DLEN][DLEN], double phi[DLEN][DLEN], int x, int y);
double dgradphi_lagrangian_y(double mz[DLEN][DLEN], double phi[DLEN][DLEN], int x, int y);
double delta_phi_C(double mz[DLEN][DLEN], double phi[DLEN][DLEN], int x, int y);
double get_max(double array[DLEN][DLEN]);
double get_min(double array[DLEN][DLEN]);
double get_residual(double array[DLEN][DLEN]);

int main()
{
  int start_s=clock();
  std::ifstream is("/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/mzdata_norm.txt");
  std::istream_iterator<double> start(is), end;
  std::vector<double> mzdatalist(start, end);
  is.close();
  int listsize = DLEN*DLEN;
  double mz [DLEN][DLEN];
  double phi [DLEN][DLEN];
  double phi_seed [DLEN][DLEN];
  double phi_gradient [DLEN][DLEN];
  double phi_gradient_last [DLEN][DLEN];

  srand (clock());

  for (int i=0; i<listsize; i++) {
    mz[(int) (i/DLEN)][i%DLEN] = 0.9999999*mzdatalist[i];
  }

  for (int x=0; x<DLEN; x++) {
    for (int y=0; y<DLEN; y++) {
      phi_seed[y][x] = angle_offset + atan2((y-y_center),(x-x_center));
      // phi_seed[y][x] = 0; ( (2*(((double) rand())/RAND_MAX) - 1)*(pi/5) )
      phi[y][x] = phi_seed[y][x];
      phi_gradient_last[y][x] = 0;
    }
  }

  double step_size = max_angle_initial;
  double max_abs_gradient = cutoff_residual+1;
  double max_gradient = 0;
  double min_gradient = 0;
  double total_residual = cutoff_residual+1;
  int i = 0;

  while (total_residual > cutoff_residual and i < max_num_steps) {
    for (int y = 2; y<DLEN-2; y++) {
      for (int x = 2; x<DLEN-2; x++) {
        phi_gradient[y][x] = delta_phi_C(mz, phi, x, y);
      }
    }
    max_gradient = get_max(phi_gradient);
		min_gradient = get_min(phi_gradient);
    if (abs(min_gradient) >  max_gradient) {
      max_abs_gradient = abs(min_gradient);
    } else {
      max_abs_gradient = max_gradient;
    }
		step_size = (max_angle_initial*exp(-i/max_angle_decay))/max_abs_gradient;
		total_residual = get_residual(phi_gradient);

    for (int y = 2; y<DLEN-2; y++) {
      for (int x = 2; x<DLEN-2; x++) {
        phi[y][x] = phi[y][x] - step_size*(phi_gradient[y][x] + momentum*phi_gradient_last[y][x]);
        phi_gradient_last[y][x] = phi_gradient[y][x];
      }
    }
    if (i%20 == 0) {
      std::cout << i << ", " << total_residual << ", " << max_abs_gradient << "\n";
    }
    ++i;
  }

  std::cout << "max_abs_gradient = " << max_abs_gradient << "\n";
  std::cout << "total_residual = " << total_residual << "\n";

  ofstream phi_file;
  phi_file.open("/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/phi_c_1760.txt");
  for (int y=0; y<DLEN; y++) {
    for (int x=0; x<DLEN-1; x++) {
      phi_file << phi[y][x] << ", ";
    }
    phi_file << phi[y][DLEN-1] << "\n";
  }
  phi_file.close();

  ofstream phi_seed_file;
  phi_seed_file.open("/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/phi_seed_c_1760.txt");
  for (int y=0; y<DLEN; y++) {
    for (int x=0; x<DLEN-1; x++) {
      phi_seed_file << phi_seed[y][x] << ", ";
    }
    phi_seed_file << phi[y][DLEN-1] << "\n";
  }
  phi_seed_file.close();

  int stop_s=clock();
  cout << "time: " << ((stop_s-start_s)/double(CLOCKS_PER_SEC)) << endl;

}

double lagrangian(double mz[DLEN][DLEN], double phi[DLEN][DLEN], int x, int y) {
  double mr = sqrt(1 - pow(mz[y][x],2));
  double dxphi = remainder((phi[y][x+1]-phi[y][x-1]),(2*pi))/2;
	double dyphi = remainder((phi[y+1][x]-phi[y-1][x]),(2*pi))/2;
	double dxmz = (mz[y][x+1]-mz[y][x-1])/2;
	double dymz = (mz[y+1][x]-mz[y-1][x])/2;

	double l = ( -(mz[y][x]/mr)*((dxmz)*cos(phi[y][x]) + (dymz)*sin(phi[y][x]))
		+ mr*(-dxphi*sin(phi[y][x]) + dyphi*cos(phi[y][x])) );

	return l;
}

double dphi_lagrangian(double mz[DLEN][DLEN], double phi[DLEN][DLEN], int x, int y) {
  double mr = sqrt(1 - pow(mz[y][x],2));
	double dxphi = remainder((phi[y][x+1]-phi[y][x-1]),(2*pi))/2;
	double dyphi = remainder((phi[y+1][x]-phi[y-1][x]),(2*pi))/2;
	double dxmz = (mz[y][x+1]-mz[y][x-1])/2;
	double dymz = (mz[y+1][x]-mz[y-1][x])/2;

	double dphi_l = ( -(mz[y][x]/mr)*(-(dxmz)*sin(phi[y][x]) + (dymz)*cos(phi[y][x]))
		+ mr*(-dxphi*cos(phi[y][x]) - dyphi*sin(phi[y][x])) );

	return dphi_l;
}

double dgradphi_lagrangian_x(double mz[DLEN][DLEN], double phi[DLEN][DLEN], int x, int y) {
  double mr = sqrt(1 - pow(mz[y][x],2));
  double dgradphi_l_x = -mr * sin(phi[y][x]);
	return dgradphi_l_x;
}

double dgradphi_lagrangian_y(double mz[DLEN][DLEN], double phi[DLEN][DLEN], int x, int y) {
  double mr = sqrt(1 - pow(mz[y][x],2));
  double dgradphi_l_y = mr * cos(phi[y][x]);
	return dgradphi_l_y;
}

double delta_phi_C(double mz[DLEN][DLEN], double phi[DLEN][DLEN], int x, int y) {
  double l = lagrangian(mz, phi, x, y);
	double dphi_l = dphi_lagrangian(mz, phi, x, y);

	double grad_ldgradphi_l =  ( ((lagrangian(mz, phi, x+1, y) * dgradphi_lagrangian_x(mz, phi, x+1, y))
						-(lagrangian(mz, phi, x-1, y) * dgradphi_lagrangian_x(mz, phi, x-1, y)))/2
						+((lagrangian(mz, phi, x, y+1) * dgradphi_lagrangian_y(mz, phi, x, y+1))
						-(lagrangian(mz, phi, x, y-1) * dgradphi_lagrangian_y(mz, phi, x, y-1)))/2 );

	return ( l * dphi_l - grad_ldgradphi_l );
}

double get_max(double array[DLEN][DLEN]) {
  double max = 0;
  for (int x = 2; x<DLEN-2; x++) {
    for (int y = 2; y<DLEN-2; y++) {
      if (array[y][x] > max) {
        max = array[y][x];
      }
    }
  }
  return max;
}

double get_min(double array[DLEN][DLEN]) {
  double min = 0;
  for (int x = 2; x<DLEN-2; x++) {
    for (int y = 2; y<DLEN-2; y++) {
      if (array[y][x] < min) {
        min = array[y][x];
      }
    }
  }
  return min;
}

double get_residual(double array[DLEN][DLEN]) {
  double res = 0;
  for (int x = 2; x<DLEN-2; x++) {
    for (int y = 2; y<DLEN-2; y++) {
        res += pow(array[y][x],2);
    }
  }
  return res;
}
