/*###########################################################################
#
#    Isolation Forest -- Anomaly detection using binary trees
#    Copyright (C) 2009 Fei Tony Liu
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
##########################################################################*/

#include <R.h>
#include "rt.h"

void rtDepth(double *x, int *xrow, int *xcol, int *samSize, int *nrnodes, int *ntree, int *hlim, int *nodeStatus, int *lDaughter, int *rDaughter, int  *splitAtt, double *splitPoint, double *ulim, double *llim, int *nSam, int *appRange, double *outF, double *pathLength, int *isoby, double *isoat)
{
  int i,j,idx ;
  double *totDepth, scaleF;
  totDepth = (double *) Calloc(*xrow, double);
  zeroInt(isoby, *xrow * *xcol);
  zeroDouble(totDepth, *xrow);
  zeroDouble(outF, *xrow);
  zeroDouble(pathLength, *xrow);
  zeroDouble(isoat, *xrow * *xcol);

    // travel the ntree number of trees
    for (j = 0; j < *ntree; ++j)  {
      idx = j * *nrnodes;

      trvTree(x, *xrow, *xcol, *hlim, nodeStatus +idx, lDaughter + idx,
             rDaughter +idx, splitAtt+idx, splitPoint+idx,
             ulim+idx, llim+idx, nSam+idx, *appRange, totDepth, isoby, isoat);
    }
    scaleF = avgPL(*samSize);
    for (i=0; i<*xrow ; ++i)
    {
      outF[i] = 1/pow(2,totDepth[i]/(scaleF * (double) *ntree));
      pathLength[i] =  totDepth[i]/ *ntree;
      for (j=0; j<*xcol ; ++j)
          if (isoby[j* *xrow +i] > 0)  // divide the isoat matrix by the number of times of the attribute is being used
                isoat[j * *xrow +i] = isoat[j * *xrow +i] /isoby[j* *xrow +i];
    }
  Free(totDepth);
}

void trvTree(double *x, int xrow, int xcol, int hlim,
             int *nodeStatus, int *lDaughter,
             int *rDaughter, int  *splitAtt, double *splitPoint,
             double *ulim, double *llim, int *nSam, int appRange, double *totDepth, int *isoby, double *isoat)
{
  int i, k, d,p,j, sAtt=0;
  double sVal, lastSplitPoint;
  bool inRange;

  for ( i = 0 ;i < xrow; ++i)  {
    k = 0; d = 0;p=0;
    while (nodeStatus[k]==NODE_INTERIOR && d<hlim) {
        //Does the terminal node applies to the sample at hand?
       // 1. Find the splitting attribute
       sAtt = splitAtt[k] - 1; // - 1 for the indexing system in R

       // 2. Is the value of sAtt within the upper and lower limits of node k?
       sVal = x[i + sAtt * xrow];
  //     isoby[sAtt * xrow + i] ++;
       if (appRange>0)
       {
             inRange=true;
             if ((sVal<llim[k])||(ulim[k]<sVal)) inRange = false;
             p += (inRange) ? 1 : 0;

             if (!inRange)
             {
                isoby[sAtt * xrow + i] ++;
                if  (sVal<llim[k])  isoat[sAtt * xrow + i] += llim[k];
                else  isoat[sAtt * xrow + i] += ulim[k];
             }
       }
       else
       
         p++;
         //isoat[sAtt * xrow + i] += splitPoint[k];

       d++;
       lastSplitPoint = splitPoint[k];  // capture the last splitpoint before k is reassigned.

       k = (sVal<splitPoint[k]) ? lDaughter[k] - 1 : rDaughter[k] - 1;
    }
    totDepth[i] += avgPL(nSam[k]) + (double) p;
    isoby[sAtt * xrow + i]++;
    isoat[sAtt * xrow + i]+= lastSplitPoint;
  }
}

