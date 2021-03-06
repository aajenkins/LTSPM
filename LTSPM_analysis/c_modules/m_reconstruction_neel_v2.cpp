/**
* @Author: Jenkins Alec <alec>
* @Date:   2017-02-10T16:37:10-08:00
* @Project: LTSPM analysis
* @Last modified by:   alec
* @Last modified time: 2017-02-20T15:02:02-08:00
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

//seed constants
const int x_center = 47;
const int y_center = 59;
const double angle_offset = -0.5+(pi/2);

//material constants
const double t = 1.0/60.0;
const double h = 72.3/60.0;

// numerical differentiation steps
const double dphi = 0.001*pi/180;
const double dmz = 1.0e-5;
const double dgradmz = 1.0e-6;

// gradient search parameters
const double momentum = 0.95;
const double cutoff_residual = 1.0e-5;
const double max_phi_step_initial = 0.05*pi/180;
const double max_mz_step_initial = 1.0e-4;
const int max_num_steps = 10;
const double max_decay_rate = 2500;

double get_delta_phi_C(double mz [DLEN][DLEN], double phi [DLEN][DLEN], double V_potential [DLEN][DLEN], double alpha_z [3*DLEN][3*DLEN], double alpha_xy [3*DLEN][3*DLEN], int x, int y);
double get_delta_mz_C(double mz [DLEN][DLEN], double phi [DLEN][DLEN], double V_potential [DLEN][DLEN], double alpha_z [3*DLEN][3*DLEN], double alpha_xy [3*DLEN][3*DLEN], int x, int y);
double get_max(double array[DLEN][DLEN]);
double get_min(double array[DLEN][DLEN]);
double get_residual(double array[DLEN][DLEN]);
double get_Fx_xy(double mz_xy, double dx_mz_xy, double phi_xy, double alpha_z [3*DLEN][3*DLEN], double alpha_xy [3*DLEN][3*DLEN], int x, int y);
double get_Fy_xy(double mz_xy, double dy_mz_xy, double phi_xy, double alpha_z [3*DLEN][3*DLEN], double alpha_xy [3*DLEN][3*DLEN], int x, int y);
double get_alpha_z(double x, double y);
double get_alpha_xy(double x, double y);

int main()
{
  // import mzeff as seed and solution to v potential as defined by Dovzhenko
  int start_s=clock();
  std::ifstream ismz("/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/mzdata_norm.txt");
  std::istream_iterator<double> start_mz(ismz), end_mz;
  std::vector<double> mzdatalist(start_mz, end_mz);
  ismz.close();

  std::ifstream isV("/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/Vdata.txt");
  std::istream_iterator<double> start_V(isV), end_V;
  std::vector<double> V_potential_datalist(start_V, end_V);
  isV.close();

  int listsize = DLEN*DLEN;

  // arrays of magnetization parameters
  double V_potential [DLEN][DLEN];
  double dx_V_potential [DLEN][DLEN];
  double dy_V_potential [DLEN][DLEN];
  double Fx_last [DLEN][DLEN];
  double Fy_last [DLEN][DLEN];
  double mz [DLEN][DLEN];
  double mz_seed [DLEN][DLEN];
  double phi [DLEN][DLEN];
  double phi_seed [DLEN][DLEN];
  double phi_gradient [DLEN][DLEN];
  double phi_gradient_last [DLEN][DLEN];
  double mz_gradient [DLEN][DLEN];
  double mz_gradient_last [DLEN][DLEN];
  double alpha_z [3*DLEN][3*DLEN];
  double alpha_xy [3*DLEN][3*DLEN];

  srand (clock());

  // convert imported data to 2d
  for (int i=0; i<listsize; i++) {
    mz_seed[(int) (i/DLEN)][i%DLEN] = 0.9999999*mzdatalist[i];
    V_potential[(int) (i/DLEN)][i%DLEN] = V_potential_datalist[i];
  }

  for (int x=0; x<DLEN; x++) {
    for (int y=0; y<DLEN; y++) {
      phi_seed[y][x] = angle_offset + atan2((y-y_center),(x-x_center));
      // phi_seed[y][x] = 0; ( (2*(((double) rand())/RAND_MAX) - 1)*(pi/5) )
      phi[y][x] = phi_seed[y][x];
      phi_gradient_last[y][x] = 0;
      mz[y][x] = mz_seed[y][x];
      mz_gradient_last[y][x] = 0;
    }
  }
  for (int x=0; x<3*DLEN; x++) {
    for (int y=0; y<3*DLEN; y++) {
      alpha_z[y][x] = get_alpha_z((double)x-DLEN, (double)y-DLEN);
      alpha_xy[y][x] = get_alpha_xy((double)x-DLEN, (double)y-DLEN);
    }
  }

  // calculate V derivative matrices and initial Fx,y arrays
  double phi_xy = phi[y][x];
  double mz_xy = mz[y][x];
  double mr_xy = sqrt(1 - pow(mz_xy,2));
  double dx_mz_xy = (mz[y][x+1]-mz[y][x-1])/2;
  double dy_mz_xy = (mz[y+1][x]-mz[y-1][x])/2;
  for (int x=0; x<DLEN; x++) {
    for (int y=0; y<DLEN; y++) {
      phi_xy = phi[y][x];
      mz_xy = mz[y][x];
      mr_xy = sqrt(1 - pow(mz_xy,2));
      dx_mz_xy = (mz[y][x+1]-mz[y][x-1])/2;
      dy_mz_xy = (mz[y+1][x]-mz[y-1][x])/2;
      dx_V_potential[y][x] = (V_potential[yp][xp+1]-V_potential[yp][xp-1])/2;
      dy_V_potential[y][x] = (V_potential[yp+1][xp]-V_potential[yp-1][xp])/2;
      Fx_last[y][x] = get_Fx
    }
  }

  double step_size_phi = max_phi_step_initial;
  double step_size_mz = max_mz_step_initial;
  double max_abs_phi_gradient = cutoff_residual+1;
  double max_abs_mz_gradient = cutoff_residual+1;
  double max_phi_gradient = 0;
  double min_phi_gradient = 0;
  double max_mz_gradient = 0;
  double min_mz_gradient = 0;
  double total_residual_phi = cutoff_residual+1;
  double total_residual_mz = cutoff_residual+1;
  int i = 0;

  while (i < max_num_steps) {
    for (int y = 2; y<DLEN-2; y++) {
      for (int x = 2; x<DLEN-2; x++) {
        phi_gradient[y][x] = get_delta_phi_C(mz, phi, V_potential, alpha_z, alpha_xy, x, y);
        mz_gradient[y][x] = get_delta_mz_C(mz, phi, V_potential, alpha_z, alpha_xy, x, y);
      }
    }

    max_phi_gradient = get_max(phi_gradient);
		min_phi_gradient = get_min(phi_gradient);
    max_mz_gradient = get_max(mz_gradient);
		min_mz_gradient = get_min(mz_gradient);
    if (abs(min_phi_gradient) >  max_phi_gradient) {
      max_abs_phi_gradient = abs(min_phi_gradient);
    } else {
      max_abs_phi_gradient = max_phi_gradient;
    }
    if (abs(min_mz_gradient) >  max_mz_gradient) {
      max_abs_mz_gradient = abs(min_mz_gradient);
    } else {
      max_abs_mz_gradient = max_mz_gradient;
    }

		step_size_phi = (max_phi_step_initial*exp(-i/max_decay_rate))/max_abs_phi_gradient;
    step_size_mz = (max_mz_step_initial*exp(-i/max_decay_rate))/max_abs_mz_gradient;
		total_residual_phi = get_residual(phi_gradient);
    total_residual_mz = get_residual(mz_gradient);

    for (int y = 2; y<DLEN-2; y++) {
      for (int x = 2; x<DLEN-2; x++) {
        phi[y][x] = phi[y][x] - step_size_phi*(phi_gradient[y][x] + momentum*phi_gradient_last[y][x]);
        phi_gradient_last[y][x] = phi_gradient[y][x];
        mz[y][x] = mz[y][x] - step_size_mz*(mz_gradient[y][x] + momentum*mz_gradient_last[y][x]);
        mz_gradient_last[y][x] = mz_gradient[y][x];
      }
    }
    if (i%1 == 0) {
      std::cout << i << ", phi res = " << total_residual_phi << ", mz res = " << total_residual_mz << "\n";
    }
    ++i;
  }

  std::cout << "total_residual_phi = " << total_residual_phi << "\n";
  std::cout << "total_residual_mz = " << total_residual_mz << "\n";

  ofstream phi_file;
  phi_file.open("/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/phi_c_1760_neel.txt");
  for (int y=0; y<DLEN; y++) {
    for (int x=0; x<DLEN-1; x++) {
      phi_file << phi[y][x] << ", ";
    }
    phi_file << phi[y][DLEN-1] << "\n";
  }
  phi_file.close();

  ofstream phi_seed_file;
  phi_seed_file.open("/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/phi_seed_c_1760_neel.txt");
  for (int y=0; y<DLEN; y++) {
    for (int x=0; x<DLEN-1; x++) {
      phi_seed_file << phi_seed[y][x] << ", ";
    }
    phi_seed_file << phi[y][DLEN-1] << "\n";
  }
  phi_seed_file.close();

  ofstream mz_file;
  mz_file.open("/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/mz_c_1760_neel.txt");
  for (int y=0; y<DLEN; y++) {
    for (int x=0; x<DLEN-1; x++) {
      mz_file << mz[y][x] << ", ";
    }
    mz_file << mz[y][DLEN-1] << "\n";
  }
  mz_file.close();

  ofstream mz_seed_file;
  mz_seed_file.open("/Users/alec/UCSB/cofeb_analysis_data/m_reconstruction/mz_seed_c_1760_neel.txt");
  for (int y=0; y<DLEN; y++) {
    for (int x=0; x<DLEN-1; x++) {
      mz_seed_file << mz_seed[y][x] << ", ";
    }
    mz_seed_file << mz[y][DLEN-1] << "\n";
  }
  mz_seed_file.close();

  int stop_s=clock();
  cout << "time: " << ((stop_s-start_s)/double(CLOCKS_PER_SEC)) << endl;

}

double get_Fx_xy(double mz_xy, double dx_mz_xy, double phi_xy, double alpha_z [3*DLEN][3*DLEN], double alpha_xy [3*DLEN][3*DLEN], int x, int y) {
  double Fx_xy = 0;
  double mr = sqrt(1 - pow(mz_xy,2));
  for (int yp=0; yp<0; yp++) {
      for (int xp=-DLEN; xp<0; xp++) {
        Fx += dx_mz_xy*alpha_z[y-yp+DLEN][x-xp+DLEN] + alpha_xy[y-yp+DLEN][x-xp+DLEN]*mr*cos(phi_xy);
      }
  }
  return Fx_xy;
}

double get_Fy_xy(double mz_xy, double dy_mz_xy, double phi_xy, double alpha_z [3*DLEN][3*DLEN], double alpha_xy [3*DLEN][3*DLEN], int x, int y) {
  double Fy = 0;
  double mr = sqrt(1 - pow(mz_xy,2));
  for (int yp=0; yp<DLEN; yp++) {
      for (int xp=0; xp<DLEN; xp++) {
        Fy += dy_mz_xy*alpha_z[y-yp+DLEN][x-xp+DLEN] + alpha_xy[y-yp+DLEN][x-xp+DLEN]*mr*sin(phi_xy);
      }
  }
  return Fy;
}

double get_alpha_z(double x, double y) {
  return (1/(2*pi))*t/(sqrt(pow(h,2)+pow(x,2)+pow(y,2)));
}

double get_alpha_xy(double x, double y) {
  return (1/(2*pi))*h*t/(pow(pow(h,2)+pow(x,2)+pow(y,2),3/2));
}

double get_delta_phi_C(double mz [DLEN][DLEN], double phi [DLEN][DLEN], double V_potential [DLEN][DLEN], double alpha_z [3*DLEN][3*DLEN], double alpha_xy [3*DLEN][3*DLEN], int x, int y) {
  double phi_xy = phi[y][x];
  double mz_xy = mz[y][x];
  double mr_xy = sqrt(1 - pow(mz_xy,2));
  double dx_mz_xy = (mz[y][x+1]-mz[y][x-1])/2;
  double dy_mz_xy = (mz[y+1][x]-mz[y-1][x])/2;

  double phi_xpyp = 0;
  double mz_xpyp = 0;
  double dx_mz_xpyp = 0;
  double dy_mz_xpyp = 0;
  double Fx = 0;
  double Fy = 0;
  double dx_V = 0;
  double dy_V = 0;

  double delta_phi_C_x = 0;
  double delta_phi_C_y = 0;

  for (int yp=0; yp<DLEN; yp++) {
    for (int xp=0; xp<DLEN; xp++) {
      phi_xpyp = phi[yp][xp];
      mz_xpyp = mz[yp][xp];
      dx_mz_xpyp = (mz[yp][xp+1]-mz[yp][xp-1])/2;
      dy_mz_xpyp = (mz[yp+1][xp]-mz[yp-1][xp])/2;
      Fx = get_Fx(mz_xpyp, dx_mz_xpyp, phi_xpyp, alpha_z, alpha_xy, xp, yp);
      Fy = get_Fy(mz_xpyp, dy_mz_xpyp, phi_xpyp, alpha_z, alpha_xy, xp, yp);
      dx_V = (V_potential[yp][xp+1]-V_potential[yp][xp-1])/2;
      dy_V = (V_potential[yp+1][xp]-V_potential[yp-1][xp])/2;

      delta_phi_C_x += Fx+dx_V;
      delta_phi_C_y += Fy+dy_V;
    }
  }

  double delta_phi_C = -mr_xy*sin(phi_xy)*delta_phi_C_x
                        +mr_xy*cos(phi_xy)*delta_phi_C_y;
	return delta_phi_C;
}

double get_delta_mz_C(double mz [DLEN][DLEN], double phi [DLEN][DLEN], double V_potential [DLEN][DLEN], double alpha_z [3*DLEN][3*DLEN], double alpha_xy [3*DLEN][3*DLEN], int x, int y) {
  double mz_xy = mz[y][x];
  double mr_xy = sqrt(1 - pow(mz_xy,2));
  double dx_mz_xy = (mz[y][x+1]-mz[y][x-1])/2;
  double dy_mz_xy = (mz[y+1][x]-mz[y-1][x])/2;

  double phi_xpyp = 0;
  double mz_xpyp = 0;
  double dx_mz_xpyp = 0;
  double dy_mz_xpyp = 0;
  double Fx = 0;
  double Fy = 0;
  double dx_V = 0;
  double dy_V = 0;

  double delta_mz_C_x = 0;
  double delta_mz_C_y = 0;

  for (int yp=0; yp<DLEN; yp++) {
    for (int xp=0; xp<DLEN; xp++) {
      phi_xpyp = phi[yp][xp];
      mz_xpyp = mz[yp][xp];
      dx_mz_xpyp = (mz[yp][xp+1]-mz[yp][xp-1])/2;
      dy_mz_xpyp = (mz[yp+1][xp]-mz[yp-1][xp])/2;
      Fx = get_Fx(mz_xpyp, dx_mz_xpyp, phi_xpyp, alpha_z, alpha_xy, xp, yp);
      Fy = get_Fy(mz_xpyp, dy_mz_xpyp, phi_xpyp, alpha_z, alpha_xy, xp, yp);

      delta_mz_C_x += Fx+dx_V;
      delta_mz_C_y += Fy+dy_V;
    }
  }

  double delta_mz_C = -((mz_xy/mr_xy)*delta_mz_C_x
                        +(mz_xy/mr_xy)*delta_mz_C_y);
	return delta_mz_C;
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

double get_residual(double array [DLEN][DLEN]) {
  double res = 0;
  for (int x = 2; x<DLEN-2; x++) {
    for (int y = 2; y<DLEN-2; y++) {
        res += pow(array[y][x],2);
    }
  }
  return res;
}
