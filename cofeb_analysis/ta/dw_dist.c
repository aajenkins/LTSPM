/**
* @Author: Jenkins Alec <alec>
* @Date:   2017-07-09T09:46:14-07:00
* @Project: LTSPM analysis
* @Last modified by:   alec
* @Last modified time: 2017-07-10T13:07:53-07:00
*/

#include <stdio.h>
#include <math.h>

void dw_dist(const void * xgridv, const void * ygridv, const void * xdwv, const void * ydwv,
  const void * theta_gradv, int gridsize, int dwsize, void * mindistv, void * theta_gridv) {

    const double * xgrid = (double *) xgridv;
    const double * ygrid = (double *) ygridv;
    const double * xdw = (double *) xdwv;
    const double * ydw = (double *) ydwv;
    const double * theta_grad = (double *) theta_gradv;

    double * mindists = (double *) mindistv;
    double * theta_grid = (double *) theta_gridv;

    double mindist;
    double tempdist;
    int minarg;
    double theta = 0;

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
      if (sqrt(pow(xgrid[i],2) + pow(ygrid[i],2)) < sqrt(pow(xdw[minarg],2) + pow(ydw[minarg],2))) {
        mindist = -mindist;
      }
      mindists[i] = mindist;
      theta_grid[i] = theta_grad[minarg];
    }

}
