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

#' @export
AnomalyScore<-function(x, forest, ntree = forest$ntree, hlim = forest$hlim, appRange=F)
{
 if (!inherits(forest,"iForest"))
    stop("forest is not a iForest object")


    rdout<-NULL

    if (!is.null(forest$trees))
    {
          samSize<- if (forest$trees$rowSamp) forest$trees$nRowSamp
                else forest$trees$xrow

          nrowx<-nrow(x)
          ncolx<-forest$trees$xcol
          xn<-x[,!forest$colisfactor]
          xn<-data.matrix(xn)
          storage.mode(xn) <- "double"

          ptime<-system.time( rdout<- .C( "rtDepth",
                       x=xn,
                       as.integer(nrowx),
                       as.integer(ncolx),
                       as.integer(samSize),
                       as.integer(forest$trees$nrnodes),
                       as.integer(ntree),
                       as.integer(hlim),
                       as.integer(forest$trees$nodeStatus),
                       as.integer(forest$trees$lDaughter),
                       as.integer(forest$trees$rDaughter),
                       as.integer(forest$trees$splitAtt),
                       as.double(forest$trees$splitPoint),
                       as.double(forest$trees$ulim),
                       as.double(forest$trees$llim),
                       as.integer(forest$trees$nSam),
                       as.integer(appRange),
                       outF = double(nrowx),
                       pathLength = double(nrowx),
                       isoby = matrix(integer(nrowx*ncolx), ncol=ncolx),
                       isoat = matrix(double(nrowx*ncolx), ncol=ncolx),
                       DUP=FALSE,
                       PACKAGE = "IsolationForest")[17:20])
          rdout$ptime<-ptime
          rm(xn)
    }

rdout

}
