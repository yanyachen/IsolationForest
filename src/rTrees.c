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

void rTrees(double *x, int *xrow, int *xcol, int *nrnodes, int *ntree, int *hlim, int *rowSamp, int *nRowSamp, int *colSamp, int *nColSamp, int *nmin, double *rFactor, double *colWeight, int *nodeStatus, int *lDaughter, int *rDaughter, int  *splitAtt, double *splitPoint, double *ulim, double *llim, int *nSam, int *ntreeSize)
{
  /********************************************************************
  input:
  x = training samples, enumerate by xrow, xcol. x[row + col * xrow] will return x(row,col)
  xrow = number of cases in x
  xcol = number of variable in x
  nrnodes = number of nodes allocated for each tree
  ntree = number of trees to construct
  hlim = height limit for building trees
  rowSamp = flag for doing row sampling
  nRowSamp = row sampling size
  colSamp = flag for doing col sampling
  nColSamp = col sampling size
  nmin = minimum leaf size
  rFactor = randomization factor
  colWeight = attribute weight
  nodeStatus = or nodes type, see rt.h
  lDaughter = references to left daughters
  rDaughter = references to right daughters
  splitAtt = attribute that was splitted
  splitPoint = the split point where a node splits
  ulim = upper limit for a node
  llim = lower limit for a node
  nSam = number of samples that a nodes has
  ntreeSize = the size (the number of nodes) of each tree
  ********************************************************************/
  int i,n,ktmp,  idx,tmpRowSamp = (*rowSamp) ? *nRowSamp : *xrow;
  int *AttPool = (int *) Calloc(*xcol, int);
  unsigned long int *xref = (unsigned long int *) Calloc(*xrow, unsigned long int);
  bool *AttSkip = (bool *)Calloc(*xcol, bool);
  /* xref is the reference to x, xref helps to save time on swapping cases*/
  for (i=0; i< *xrow; i++) xref[i] = (unsigned long int)i;

  GetRNGstate(); // required for using random number generator, call this before using any random number, e.g. shuffle()
  /* initialize the AttPool array */
  for (i = 0;i <*xcol; ++i) AttPool[i] = i;
  /* prepare AttPool and AttSkip arrays */
  if (*colSamp){
    weightedAttPool(colWeight,AttPool, *nColSamp, *xcol);
    for (i = 0; i<*xcol; ++i) AttSkip[i] = true;
    for (i = 0; i<*nColSamp; ++i) AttSkip[AttPool[i]] = false;
  }
  else {
    for (i = 0; i<*xcol; ++i) AttSkip[i] = false;
    *nColSamp = *xcol;
  }

  for ( i=0 ;i < *ntree; ++i)  {
      idx = i * *nrnodes;  //tree offset
      /* Random sampling */
      if (*rowSamp) {
          /* randomly shuffle the first nRowSamp samples*/

          for (n=0; n<*nRowSamp ; ++n)
          {
            ktmp= (int) (ceil(*xrow * unif_rand())-1); /* only shuffle the data in the sample area*/
              swapUnsignedLongInt(n,ktmp,xref);
          }
      }
    /*build a tree using idx for random samples*/
    rTree(x, xref, *xrow, *xcol, tmpRowSamp, *nrnodes, *hlim, *nmin, *rFactor, AttSkip, AttPool, *nColSamp, nodeStatus+idx, lDaughter+idx, rDaughter+idx, splitAtt+idx, splitPoint+idx, ulim+idx, llim+idx, nSam+idx, ntreeSize+i);
  }

  PutRNGstate(); /* required for using random number generator */
  Free(AttPool);
  Free(xref);
  Free(AttSkip);
}



/* build a single tree, tree matrices are offsetted (nodeStatus, lDaughter, rDaughter,*/
/* splitAtt, splitPoint, ulim, llim, nSam. */
void rTree(double *x, unsigned long int *xref, int xrow, int xcol,  int nsample,int nrnodes, int hlim, int nmin, double rFactor, bool *AttSkip, int *AttPool, int nColSamp, int *nodeStatus, int *lDaughter, int *rDaughter, int *splitAtt, double *splitPoint, double *ulim, double *llim, int *nSam, int *ntreeSize)
{
  int k,nCur, ndStart, ndEnd,ndEndl, *nodeStart, *nodeDepth, *nodeSplitBy;
  double rNum=0;
  /* Pointers:
     k - node to split
     nCur - node that is available to be added
     nodeStart - start of data in x for a nodes
     nSam - number of samples for a nodes
     nodeDepth - node depth, root node has zero depth value
     ndStart - start of data in x for working node k
     ndEnd - end of data in x for working node k
     nodeSplitBy - the reference no of the parent node*/

  nodeStart = (int *) Calloc(nrnodes, int);
  nodeDepth = (int *) Calloc(nrnodes, int);
  nodeSplitBy = (int *) Calloc(nrnodes, int);
  /* initialize some arrays for building tree */
  zeroInt(nodeStatus,nrnodes);
  zeroInt(nodeStart,nrnodes);
  zeroInt(nSam,nrnodes);
  zeroInt(nodeDepth,nrnodes);
  zeroInt(nodeSplitBy, nrnodes);  //-1 means split by no one

  nCur=0;
  nodeStart[0] = 0;
  nodeDepth[0] = 0;
  nSam[0] = nsample;
  nodeStatus[0] = NODE_TOSPLIT;
  nodeSplitBy[0] = -1;

  /* main loop for building tree */
  for (k=0; k<nrnodes-2;++k){
    if (k>nCur || nCur >=nrnodes - 2 ) break;
    /* skip if node is not to be split */
    if (nodeStatus[k] != NODE_TOSPLIT) continue;
    /* Depth control */
    if (nodeDepth[k]>=hlim) {
      nodeStatus[k] = NODE_TERMINAL;
      continue;
    }
    /* initialize for the next call of randomSplit */
    ndStart = nodeStart[k];
    ndEnd = ndStart + nSam[k] - 1;
    rNum =  unif_rand();

    if (rNum <= rFactor)
        randomSplit(x, xref, xrow, xcol, ndStart, ndEnd, splitAtt+k, splitPoint+k,ulim+k,llim+k, &ndEndl, AttPool, nColSamp);
    else
        SDGainSplit(x, xref, xrow, xcol, ndStart, ndEnd, splitAtt+k, splitPoint+k,ulim+k,llim+k, &ndEndl, AttSkip);


    if (splitAtt[k]==-1)
    {
      nodeStatus[k] = NODE_TERMINAL;
      continue;
    }
    nodeStatus[k] = NODE_INTERIOR;
    nodeSplitBy[nCur + 1] = nodeSplitBy[nCur + 2] = k;

    /* leftnode = ncur +1, rightnode = ncur+2*/
    nSam[nCur + 1] = ndEndl - ndStart + 1;
    nSam[nCur + 2] = ndEnd - ndEndl;
    nodeStart[nCur + 1] = ndStart;
    nodeStart[nCur + 2] = ndEndl + 1;
    /* update nodeDepth */
    nodeDepth[nCur + 1] = nodeDepth[nCur + 2] = nodeDepth[k] + 1 ;


    nodeStatus[nCur + 1] = nodeStatus[nCur + 2] = ((nodeDepth[k]+1) >=hlim) ? NODE_TERMINAL : NODE_TOSPLIT;

    if (nSam[nCur + 1] <=nmin )
            nodeStatus[nCur + 1] = NODE_TERMINAL;
    if (nSam[nCur + 2] <=nmin )
            nodeStatus[nCur + 2] = NODE_TERMINAL;


    lDaughter[k] = nCur + 1 + 1; /* + 1 for R indexing */
    rDaughter[k] = nCur + 2 + 1; /* + 1 for R indexing */

    nCur += 2;

  }
  *ntreeSize=k; //capturing the tree size
  Free(nodeSplitBy);
  Free(nodeDepth);
  Free(nodeStart);
}



void SDGainSplit (double *x, unsigned long int *xref, int xrow, int xcol, int ndStart, int ndEnd, int *splitAtt, double *splitPoint, double *ulim, double *llim, int *ndEndl, bool *AttSkip)
{
   //double attMax, attMin, attSplitPoint, attMaxGain=0, maxGain,  baseSd, leftSd, rightSd,  baseSum, baseSqSum, leftSum, leftSqSum,  rightSum, rightSqSum,Gain,xvalue;
   double attMax, attMin, attSplitPoint, attMaxGain=0, maxGain,  baseSd, leftSd, rightSd,  Gain, base;
   long double baseOldM, baseOldQ, leftOldM, leftOldQ,  rightOldM, rightOldQ,baseNewM, baseNewQ, leftNewM, leftNewQ,  rightNewM, rightNewQ,xvalue;
   int n, sAtt, maxGainN;
   bool FirstRun = true;


   *splitAtt = -1;
   if ((ndEnd - ndStart + 1)<3) return; // require three points to split
   for ( sAtt = 0; sAtt< xcol ; ++sAtt)
   {
     if (AttSkip[sAtt]) continue;
     Quicksort( ndStart,ndEnd,sAtt, x, xref, xrow, xcol);
     attMax = x[(int)xref[ndEnd] + sAtt * xrow];
     attMin = x[(int)xref[ndStart] + sAtt * xrow];
     if(attMax==attMin) continue; //no point to step down

     for ( n =ndStart ; n <= ndEnd ; n++)
     {
       xvalue = (long double)x[(int)xref[n]+ sAtt * xrow];
       if (n == ndStart)
       {
         baseOldM =  xvalue;
         baseOldQ = 0;
       }
       else
       {
         baseNewM = baseOldM + (xvalue - baseOldM) / (long double) (n - ndStart +1);
         baseNewQ = baseOldQ + (xvalue - baseOldM) * (xvalue - baseNewM);

         baseOldM = baseNewM;
         baseOldQ = baseNewQ;
       }
     }
     baseSd = sqrt((double)(baseNewQ/(long double)(ndEnd - ndStart+1)));

     rightOldM = baseNewM;
     rightOldQ = baseNewQ;
     for ( n =ndStart ; n < ndEnd ; n++)
     {
      xvalue = x[(int)xref[n]+ sAtt * xrow];

      if (n == ndStart)
      {
        leftOldM = xvalue;
        leftOldQ = 0;
        leftSd = 0;
      }
      else
      {
        leftNewM = leftOldM + (xvalue - leftOldM) / (long double) (n - ndStart + 1);
        leftNewQ = leftOldQ + (xvalue - leftOldM) * (xvalue - leftNewM);

        leftOldM = leftNewM;
        leftOldQ = leftNewQ;

        leftSd = sqrt((double)(leftNewQ/(long double)(n - ndStart +1)));
      }

      if ( (n+1) == ndEnd)
      {
         rightSd = 0;
      }
      else
      {
        rightNewM = (rightOldM * (long double)(ndEnd - n + 1) - xvalue) / (long double)(ndEnd - n);
        rightNewQ =  rightOldQ - (xvalue - rightOldM) * (xvalue - rightNewM);
        rightSd = sqrt((double)(rightNewQ/(long double)(ndEnd - n )));

        rightOldM = rightNewM;
        rightOldQ = rightNewQ;
      }
      base = 2;
      //if (leftSd==0 || rightSd==0) base =1;    // if one of the Sd is zero then base of the average is 1

      Gain = (baseSd - (leftSd+rightSd)/base)/baseSd;
       //Rprintf("sAtt %d,baseSd %f, leftSd %f, rightSd %f, rightSum %f, rightSqSum %f, Gain %f\n",sAtt, baseSd, leftSd, rightSd, rightSum, rightSqSum, Gain);

       if (Gain >= maxGain || n == ndStart)
       {
         maxGainN=n;
         maxGain=Gain;
         attSplitPoint = (x[(int)xref[n+1]+ sAtt * xrow] - x[(int)xref[n]+ sAtt * xrow])/2 + x[(int)xref[n]+ sAtt * xrow];
       }
     }
     if (maxGain > attMaxGain||FirstRun)
     {
        attMaxGain = maxGain;
        *splitAtt = sAtt + 1;
        *splitPoint = attSplitPoint;
        *ndEndl = maxGainN;
        FirstRun=false;
        *ulim = attMax + attSplitPoint - attMin;
        *llim = attMin - (attMax - attSplitPoint);
     }
   }

   if (*splitAtt>-1)
   //quicksort here to restore the ordering of attribute splitAtt
  // Quicksort( ndStart,ndEnd,(*splitAtt-1), x, xref, xrow, xcol);
   *ndEndl = SortAndLocate(ndStart,ndEnd,(*splitAtt-1), x, xref, xrow, xcol, *splitPoint);
}

void randomSplit(double *x, unsigned long int *xref, int xrow, int xcol, int ndStart, int ndEnd, int *splitAtt, double *splitPoint, double *ulim, double *llim, int *ndEndl, int *AttPool, int nColSamp)
{
   double xMax , xMin, *attMax, *attMin, *attSplitPoint;
//   int n, sAtt, xHigh, xLow, sAttIdx, temp;
   int  sAtt, sAttIdx, i,j , *localAttPool;
   //Rprintf("Random Split\n");
   localAttPool = (int *)Calloc(nColSamp, int);
   attMax= (double *)Calloc(xcol, double);
   attMin= (double *)Calloc(xcol, double);
   attSplitPoint= (double *)Calloc(xcol, double);

   for (i = 0 ; i< nColSamp; i++) localAttPool[i]=AttPool[i];  // copy a local version of AttPool
   for (i = ndStart; i<=ndEnd;i++)  // find out all the min and max for each attribute
     for (j = 0; j<xcol; j++)
     {
         if((i==ndStart)||(x[(int)xref[i]+j*xrow] > attMax[j])) attMax[j] = x[(int)xref[i]+j*xrow];
         if((i==ndStart)||(x[(int)xref[i]+j*xrow] < attMin[j])) attMin[j] = x[(int)xref[i]+j*xrow];
     }
   for (j=0; j<xcol; j++)
   {
      attSplitPoint[j] = (ceil(unif_rand() * RAND_MAX) /RAND_MAX) * (attMax[j] - attMin[j]) + attMin[j];
      //   (ceil(unif_rand() * RAND_MAX) /RAND_MAX) changes the range to (0,1] instead of [0,1] in unif_rand()
  //    ulim[j] = attMax[j] + attSplitPoint[j] - attMin[j];
  //    llim[j] = attMin[j] - (attMax[j] - attSplitPoint[j]);
   }


   *splitAtt=-1;
   while (*splitAtt==-1 && nColSamp >0)
   {
     /* This while loop exhaust other attributes to split,
     if a randomly selected attribute cannot be split, then stop splitting.*/
     //sAtt = (int) (unif_rand() * xcol); // Randomly select an attribute
     sAttIdx = (int) (ceil(unif_rand() * (double) nColSamp) - 1);
     sAtt = localAttPool[sAttIdx];
     // if selected attribute is the same as what is being splitted by then select another attribute
     /* sort the cases between ndStart and ndEnd using splitAtt */
     /* find the min and max of splitAtt values as well */

     xMax = attMax[sAtt];
     xMin = attMin[sAtt];
     if (xMax==xMin) {
       // swapping AttPool[sAttIdx] and sAttIdx[nColSamp - 1]
       swapint(sAttIdx,nColSamp - 1, localAttPool);
/*       temp = localAttPool[sAttIdx];
       localAttPool[sAttIdx] = localAttPool[nColSamp - 1];
       localAttPool[nColSamp - 1] = temp;                 */
       nColSamp--; // reducing the pool size
       continue;
     }

     *splitPoint = attSplitPoint[sAtt];

     *ndEndl = SortAndLocate(ndStart,ndEnd,sAtt, x, xref, xrow, xcol, *splitPoint);
     // Rprintf("n2 = %d, splitpoint: %f, x[n]: %f, x[n+1]: %f\n",n, *splitPoint,x[n+sAtt * xrow], x[n+1+sAtt * xrow] );
     *splitAtt = sAtt + 1;  /* + 1 for indexing in R */
   }
   
   if (*splitAtt > -1)
   {
       *ulim = attMax[*splitAtt -1] + attSplitPoint[*splitAtt -1] - attMin[*splitAtt -1];
       *llim = attMin[*splitAtt -1] - (attMax[*splitAtt -1] - attSplitPoint[*splitAtt -1]);
   }

   Free(attSplitPoint);
   Free(attMin);
   Free(attMax);
   Free(localAttPool);


}

