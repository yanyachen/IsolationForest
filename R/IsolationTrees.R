###########################################################################
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
##########################################################################

IsolationTrees<-function(x, ntree=10, hlim=as.integer(ceiling(log2(nrow(x)))), rowSamp=F, nRowSamp=nrow(x), nmin = 1, rFactor =1, colSamp = F, nColSamp = ncol(x), colWeight = c(rep(1,ncol(x))))
# x is a data frame on which the random trees are built.
# ntree is the number of trees to build.
# hlim is the height limit for trees, its default value is log2(nrow(x))
#rowSamp is a logical swith to perform random instance sub-sampling
#nRowSamp is the instance sub-sampling size; it must be less than or equal to the training sample size
# colSamp is alogical switch to perform random attribute sub-sampling}
# nColSamp is the attribute sub-sampling size; it must be less than or equal to the number of attributes
{
     #Trees are implemented as a set of matrices, each row is a tree
     #lDaughter and rDaughter are column pointers to find the respective left and right daughters
     #nodeStatus, see rt.h for nodeStatus
     #splitAtt indicates which attribute is splitted
     #splitPoint store the split point values
     #ulim is the upper limit of splitAtt value to which the node applies
     #llim is the lower limit of splitAtt value to which the node applies
     #nSam is the number of samples that the node had.
     colWeight<-colWeight
     hlim<-hlim
     nRowSamp<-nRowSamp
     nColSamp<-nColSamp
     rowSamp<-rowSamp
     colSamp<-colSamp
     ntree<-ntree

     if (rowSamp && (nRowSamp>= nrow(x))) stop("nRowSamp must be smaller than nrow(x)")
     if (colSamp && (nColSamp>= ncol(x))) stop("nColSamp must be smaller than ncol(x)")

     if (hlim<2)
     {hlim<-2}

     nrnodes <- min(if(rowSamp) nRowSamp else nrow(x), 2^hlim) * 2 - 1

     rtout<-NULL
     ptime<-NULL
     colisfactor<-FALSE
     drtout<-NULL
     dptime<-NULL
     allxcols<-colnames(x)
     if (is.data.frame(x))
     {
        colisfactor<-unlist(lapply(x, is.factor))
        colnlevels<-unlist(lapply(x, nlevels))
     }

     if (any(!colisfactor))
     {
         ncolWeight <- colWeight[!colisfactor]
         xcols<-colnames(x)
         x<-data.matrix(x)
         nrowx<-nrow(x)
         ncolx<-ncol(x)
         storage.mode(x) <- "double"
         ptime<-system.time(rtout<- .C(    "rTrees",
                          x = x,
                          xrow = as.integer(nrowx),
                          xcol = as.integer(ncolx),
                          nrnodes = as.integer(nrnodes),
                          ntree = as.integer(ntree),
                          hlim = as.integer(hlim),
                          rowSamp = as.integer(rowSamp),
                          nRowSamp = as.integer(nRowSamp),
                          colSamp = as.integer(colSamp),
                          nColSamp = as.integer(nColSamp),
                          nmin = as.integer(nmin),
                          rFactor = as.double(rFactor),
                          colWeight = as.double(ncolWeight),
                          nodeStatus = matrix(integer(nrnodes*ntree), ncol=ntree),
                          lDaughter = matrix(integer(nrnodes*ntree), ncol=ntree),
                          rDaughter = matrix(integer(nrnodes*ntree), ncol=ntree),
                          splitAtt = matrix(integer(nrnodes*ntree), ncol=ntree),
                          splitPoint = matrix(double(nrnodes*ntree), ncol=ntree),
                          ulim = matrix(double(nrnodes*ntree), ncol = ntree),
                          llim = matrix(double(nrnodes*ntree), ncol = ntree),
                          nSam = matrix(integer(nrnodes*ntree), ncol=ntree),
                          ntreeSize = integer(ntree),
                          DUP=FALSE,
                          PACKAGE = "IsolationForest")[2:22])
     }
     out<-list(xcols=allxcols,
          trees=rtout,
          ptime=ptime,
          colisfactor=colisfactor,
          colnlevels=colnlevels,
          ntree = ntree,
          hlim = hlim)
     class(out)<-"iForest"
     return(out)
}