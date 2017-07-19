/**
* @Author: Jenkins Alec <alec>
* @Date:   2017-07-09T09:46:14-07:00
* @Project: LTSPM analysis
* @Last modified by:   alec
* @Last modified time: 2017-07-09T12:32:54-07:00
*/

#include <stdio.h>
#include <math.h>

void dw_dist(const void * xgridv, const void * ygridv, const void * xdwv, const void * ydwv,
  int gridsize, int dwsize, void * mindistv, void * minargv) {

    const double * xgrid = (double *) xgridv;
    const double * ygrid = (double *) ygridv;
    const double * xdw = (double *) xdwv;
    const double * ydw = (double *) ydwv;

    double * mindists = (double *) mindistv;
    double * minargs = (double *) minargv;

    double mindist;
    double tempdist;
    int minarg;

    for (int i = 0; i < gridsize * gridsize; ++i) {
      mindist = sqrt( pow(xgrid[i]-xdw[0],2) + pow(ygrid[i]-ydw[0],2) );
      minarg = 0;
      for (int m = 1; m < dwsize; ++m) {
        tempdist = sqrt( pow(xgrid[i]-xdw[0],2) + pow(ygrid[i]-ydw[0],2) );
        if (tempdist < mindist) {
          mindist = tempdist;
          minarg = m;
        }
      }
      mindists[i] = mindist;
      minargs[i] = minarg;
    }

}
