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

unsigned char pack(int nBits, unsigned char *bits) {
    register int i = nBits, temp;
    //Rprintf("pack %d, ", nBits);
    register unsigned char pack=0;
    while (--i >= 0) { temp = (int) bits[i]; pack += (unsigned char) temp << i;}
    return(pack);
}

void unpack(unsigned char pack, unsigned char *bits) {
/* pack is a 1-byte unsigned char.  The sub. returns icat, an integer array of
   zeroes and ones corresponding to the coefficients in the binary expansion 
   of pack. */   
    int i;
    for (i = 0; pack != 0; pack >>= 1, ++i) bits[i] = pack & 1;
}

void packAnode(unsigned char *battValues, int sumnlevels, unsigned char *attValues)
{
  int i;
  for (i=0; i<(int)ceil((double) sumnlevels/8.0); i++)
  {
    if (i==(int)(floor((double) sumnlevels/8.0)))
    {
      if((sumnlevels%8)>0)
         battValues[i]=pack(sumnlevels%8,attValues+i*8);
      else
         battValues[i]=pack(8,attValues+i*8);
    }
    else
      battValues[i]=pack(8,attValues+i*8);
  }
}

void unpackAnode(unsigned char *battValues, int sumnlevels, unsigned char *attValues)
{
    int i;
    for (i=0; i<(int)ceil((double)sumnlevels/8.0); i++)
        unpack(battValues[i], attValues+i*8);
}


void shuffle(int *intArr, int size)
{
    int temp, rNum, last;

    for (last = size; last > 1; last--)
    {
           rNum = (int) floor(unif_rand() * (double)last);
           temp = intArr[rNum];
           intArr[rNum] = intArr[last - 1];
           intArr[last - 1] = temp;
    }
}// end shuffle( )

void weightedShuffle(double *colWeight, int *AttPool, int nColSamp, int xcol)
{ // This function takes in the column weight and Shuffle the attribute according random numbers.
  double tempAttScore, *AttScore = (double *)Calloc(xcol, double) ;
  int i,j, tempj;
  for (i=0;i<xcol;++i)
      AttScore[i] = colWeight[i] * unif_rand();
  for (i=0;i<nColSamp;++i){
    tempAttScore = AttScore[i];
    tempj = i;
    for (j = i+1; j<xcol; ++j)
    {
      if(AttScore[j] > tempAttScore) {
        tempAttScore = AttScore[j];
        tempj = j;
      }
    }
    AttScore[tempj] = AttScore[i];
    AttPool[tempj] = AttPool[i];
    
    AttScore[i] = tempAttScore;
    AttPool[i] = AttPool[tempj];
  }
  Free(AttScore);
}


void weightedAttPool(double *colWeight, int *AttPool, int nColSamp, int xcol)
{ // This function takes in the column weight and selects nColSamp attributes according to column weight.
  double tempColWeight ;
  int i,j, tempj;

  for (i=0;i<nColSamp;++i){
    tempColWeight = colWeight[i];
    tempj = i;
    for (j = i+1; j<xcol; ++j)
    {
      if(colWeight[j] > tempColWeight) {
        tempColWeight = colWeight[j];
        tempj = j;
      }
    }
    colWeight[tempj] = colWeight[i];
    AttPool[tempj] = AttPool[i];

    colWeight[i] = tempColWeight;
    AttPool[i] = AttPool[tempj];
  }
}


void zeroInt(int *x, int length) {
    memset(x, 0, length * sizeof(int));
}

void zeroUnsignedChar(unsigned char *x, int length) {
    memset(x, 0, length * sizeof(unsigned char));
}


void zeroDouble(double *x, int length) {
    memset(x, 0, length * sizeof(double));
}

void Quicksort (int Fp,int Lp,int Att,double *x, unsigned long int *xref, int xrow, int xcol ) //adopted from C4.5
{

    register int Lower, Middle;
    register double Thresh;
    register int i;

    if (Fp < Lp)
    {
        Thresh = x[(int)xref[Lp] + Att * xrow];//CVal (Item[Lp], Att)

        /* Isolate all items with values <= threshold */

        Middle = Fp;

        for (i = Fp; i < Lp; i++)
        {
            if (x[(int)xref[i] + Att * xrow] <= Thresh)
            {
                if (i != Middle)
//                    swapCases (Middle, i,x,xrow,xcol);
                    swapUnsignedLongInt(Middle, i, xref);
                Middle++;
            }
        }

        /* Extract all values equal to the threshold */

        Lower = Middle - 1;

        for (i = Lower; i >= Fp; i--)
        {
            if (x[(int)xref[i] + Att* xrow] == Thresh)
            {
                if (i != Lower)
                 //   swapCases (Lower,i,x,xrow,xcol);
                    swapUnsignedLongInt(Lower,i,xref);
                Lower--;
            }
        }

        /* Sort the lower values */

        Quicksort (Fp, Lower, Att, x,xref,xrow,xcol);

        /* Position the middle element */
         swapUnsignedLongInt(Middle, Lp,xref);
        /* Sort the higher values */

        Quicksort (Middle + 1, Lp, Att, x,xref,xrow,xcol);
    }
}

int SortAndLocate(int Fp, int Lp, int Att, double *x, unsigned long int *xref, int xrow, int xcol, double splitPoint)
{
  register int Middle, i;
  Middle = Fp;
  for (i = Fp; i<=Lp; i++)
  {
     if ( x[(int)xref[i]+Att*xrow] <= splitPoint)
     {
       if (i != Middle) swapUnsignedLongInt(Middle, i,xref);
       //swapCases(Middle, i,x, xrow, xcol);
       Middle++;
     }
  }
  return Middle-1 ;
}


int dSortAndLocate(int Fp, int Lp, int Att, double *x, unsigned long int *xref, int xrow, int xcol, unsigned char *licat)
{
  register int Middle, i;
  Middle = Fp;
  for (i = Fp; i<=Lp; i++)
  {
     if ( licat[(int)(x[(int)xref[i]+Att*xrow]-1.0)] >0)
     {
       if (i != Middle) swapUnsignedLongInt(Middle, i,xref);
       Middle++;
     }
  }
  return Middle-1 ;
}


double rtMean(int Fp,int Lp,int Att,double *x, unsigned long int *xref, int xrow, int xcol)
{
    register int i;
    register double tempMean=0.0;

    for (i = Fp; i <= Lp; i++)
    {
       tempMean += x[(int)xref[i] + Att* xrow];
    }
    return tempMean/(double)(Lp - Fp + 1);
}

double rtSd(int Fp,int Lp,int Att,double *x, unsigned long int *xref, int xrow, int xcol, double tempMean)
{
    register int i;
    double tempSd =0.0;

    for (i = Fp; i <= Lp; i++)
    {
       tempSd += pow((x[(int)xref[i] + Att* xrow] - tempMean), 2.0);
    }
    return tempSd/(double)(Lp - Fp + 1);
}

double rtStdMoment(int Fp,int Lp,int Att,double *x, unsigned long int *xref, int xrow, int xcol, int order, double tempMean, double tempSd)
{
    register int i;
    double tempSum = 0.0;

    if (tempSd <=0)
       return NA_REAL;
    else
    {
      for (i = Fp; i <= Lp; i++)
      {
         tempSum += pow(fabs(x[(int)xref[i] + Att* xrow] - tempMean), (double) order);
         //tempSum += fabs(x[i + Att* xrow] - tempMean);
      }
     return (tempSum / (double)(Lp - Fp + 1))/ pow(tempSd,(double) order) ;
    }
}


double meanIn (int n, double xi, double Oldmean)
{return (Oldmean * (double) n + xi)/(double)(n+1);}

double meanEx (int n, double xi, double Oldmean)
{return (Oldmean * (double) n - xi)/(double)(n-1);}


int locate(double *x, unsigned long int *xref, int xrow, int sAtt, int xHigh, int xLow, double searchPoint)
{
     /* return the instance number of x before the searchPoint */
     int n;
     // there are faster way to find ndEndl, to keep splitting in the middle of ndStart and ndEnd.
     // for (n = ndStart; n<=ndEnd; ++n)
     //   if (x[n + sAtt * xrow]<*splitPoint) *ndEndl = n;
     //   Rprintf("n1 = %d\n",*ndEndl);
     //faster way to find ndEndl, to keep splitting in the middle of xHigh and xLow.
     n = (xHigh + xLow)/2;
     while (!((x[(int)xref[n]+sAtt * xrow] < searchPoint) && (x[(int)xref[n+1]+sAtt * xrow] >= searchPoint)))
     {
       if ( x[(int)xref[n] + sAtt * xrow] < searchPoint)
        xLow = n;
       else
        xHigh = n;
       n = (xHigh + xLow)/2;
       if (xHigh == xLow)
        break;
     }
     return n;
}


double sfunction(double v)
{
  double t = v*2.0 - 1.0;
  return t / (1.0+ pow(t,2.0)) + 0.5;
}


void swapUnsignedLongInt (int a, int b, unsigned long int *x)
{
  register unsigned long int hold;
  hold = x[a];
  x[a] = x[b];
  x[b] = hold;
}


void swapint (int a, int b, int *x)
{
  register int hold;
  hold = x[a];
  x[a] = x[b];
  x[b] = hold;
}

void swapDouble (int a, int b,double *x)
{
       register double hold;
       hold = x[a];
       x[a] = x[b];
       x[b] = hold;
}

void swapCases (int a, int b, double *x, int n, int m)
{
  /* n is the number of rows*/
  /* m is the number of cols*/
  register int i;
  for (i=0;i < m;++i)
    swapDouble( a + i * n, b + i * n, x);
}


double Gaussian(double v, double a, double b, double c)
{
       return  a*exp(-pow(v-b,2)/(2*pow(c,2))); // please see http://en.wikipedia.org/wiki/Gaussian_function
}
