/**
* @Author: Jenkins Alec <alec>
* @Date:   2017-07-09T09:46:14-07:00
* @Project: LTSPM analysis
 * @Last modified by:   alec
 * @Last modified time: 2017-07-30T10:41:32-07:00
*/

#include <stdio.h>
#include <math.h>

const double pi = 3.14159265358979323846;

int mod(int a, int b);

void dw_dist(const double * xgrid, const double * ygrid, const double * xdw, const double * ydw,
   const double * sdw, int gridsize, int dwsize, double * mindists) {

    double r;
    double r1;
    double r2b;
    double r2;
    double mindist;
    double tempdist;
    int minarg;
    int r1arg;

    for (int i = 0; i < gridsize*gridsize; ++i) {
      mindist = sqrt( pow(xgrid[i]-xdw[0],2) + pow(ygrid[i]-ydw[0],2) );
      minarg = 0;
      for (int m = 1; m < dwsize; ++m) {
        tempdist = sqrt( pow(xgrid[i]-xdw[m],2) + pow(ygrid[i]-ydw[m],2) );
        if (tempdist < mindist) {
          mindist = tempdist;
          minarg = m;
        }
      }
      r1 = mindist;
      r1arg = minarg;
      r2 = sqrt( pow(xgrid[i]-xdw[mod(minarg+1, dwsize)],2) + pow(ygrid[i]-ydw[mod(minarg+1, dwsize)],2) );

      r2b = sqrt( pow(xgrid[i]-xdw[mod(minarg-1, dwsize)],2) + pow(ygrid[i]-ydw[mod(minarg-1, dwsize)],2) );
      if (r2b <= r2) {
        r2 = mindist;
        r1 = r2b;
        r1arg = mod(minarg-1, dwsize);
        // printf("%d\n", minarg);
      }

      r = r2 * sqrt( 1 - pow( (pow(sdw[r1arg],2) + pow(r2,2) - pow(r1,2)) / (2*sdw[r1arg]*r2) ,2) );

      // if (sqrt(pow(xgrid[i],2) + pow(ygrid[i],2)) < sqrt(pow(xdw[minarg],2) + pow(ydw[minarg],2))) {
      //   r = -r;
      //   //mindist = -mindist;
      // }
      mindists[i] = r;
    }


}

int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}
