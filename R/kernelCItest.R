#' Kernel Conditional Independence test
#'
#' Test to check the (conditional) dependence between two variables x and y given a set of variables S, using independence criteria. The kernelCItest() function, uses Distance Covariance or Hilbert-Schmidt Independence Criterion to test for the (conditional) independence between random variables, with an interface that can easily by used in \code{\link{skeleton}}, \code{\link{pc}} or \code{\link{kpc}}.
#' @param x,y,S It is tested, whether x and y are conditionally independent given the subset S of the remaining nodes. x, y, S all are integers, corresponding to variable or node numbers.
#' @param verbose a logical parameter, if TRUE, detailed output is provided.
#' @param suffStat a list of parameters consisting of data, ic.method, p, index, sig, numCol, numCluster, eps, paral
#' @param data numeric matrix witch collumns representing variables and rows representing samples
#' @param ic.method Method for the (conditional) independence test: Distance Covariance (permutation or gamma test), HSIC (permutation or gamma test) or HSIC cluster
#' @param p Number of permutations for Distance Covariance, HSIC permutation and HSIC cluster tests. Default is Distance Covariance
#' @param index Number in (0,2] the power of the distance in the Distance Covariance
#' @param sig Gaussian kernel width for HSIC tests. Default is 1
#' @param numCol Number of columns used in the incomplete Cholesky decomposition. Default is 50
#' @param numCluster Number of clusters for kPC clust algorithm
#' @param eps Normalization parameter for kPC clust. Default is 0.1
#' @param paral Number of cores to use for parallel calculations.
#'
#'
#' @importFrom energy dcov.test
#' @import methods
#'
#' @references G. Szekely, M. Rizzo and N. Bakirov (2007). Measuring and Testing Dependence by Correlation of Distances. The Annals of Statistics 2007, Vol. 35, No. 6, 2769-2794.
#' @references A. Gretton et al. (2005). Kernel Methods for Measuring Independence. JMLR 6 (2005) 2075-2129.
#' @references R. Tillman, A. Gretton and P. Spirtes (2009). Nonlinear directed acyclic structure learning with weakly additive noise model. NIPS 22, Vancouver.
#'
#' @export
#' @return kernelCItest() returns the p-value of the test.
#' @author Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk}) and Nina Ines Bertille Desgranges
#' @examples
#' set.seed(10)
#' library(pcalg)
#' z <- 10*runif(300)
#' w <- 10*runif(300)
#' x <- sin(z) + runif(300)
#' y <- cos(z) + runif(300)
#'
#' data <- cbind(x,y,z,w)
#'
#' #conditionally independent
#' test1a <- kernelCItest(x=1,y=2,S=c(3),suffStat = list(data=data,ic.method="dcc.gamma"))
#' test2a <- kernelCItest(x=1,y=2,S=c(3),suffStat = list(data=data,ic.method="dcc.perm"))
#' test3a <- kernelCItest(x=1,y=2,S=c(3),suffStat = list(data=data,ic.method="hsic.gamma"))
#' test4a <- kernelCItest(x=1,y=2,S=c(3),suffStat = list(data=data,ic.method="hsic.perm"))
#' test5a <- kernelCItest(x=1,y=2,S=c(3),suffStat = list(data=data,ic.method="hsic.clust"))
#' test6a <- gaussCItest( x=1,y=2,S=c(3),suffStat = list(C=cor(data),n=4))
#'
#' test1a
#' test2a
#' test3a
#' test4a
#' test5a
#' test6a
#'
#' #dependent
#' test1b <- kernelCItest(x=1,y=2,S=c(4),suffStat = list(data=data,ic.method="dcc.gamma"))
#' test2b <- kernelCItest(x=1,y=2,S=c(4),suffStat = list(data=data,ic.method="dcc.perm"))
#' test3b <- kernelCItest(x=1,y=2,S=c(4),suffStat = list(data=data,ic.method="hsic.gamma"))
#' test4b <- kernelCItest(x=1,y=2,S=c(4),suffStat = list(data=data,ic.method="hsic.perm"))
#' test5b <- kernelCItest(x=1,y=2,S=c(4),suffStat = list(data=data,ic.method="hsic.clust"))
#' test6b <- gaussCItest( x=1,y=2,S=c(4),suffStat = list(C=cor(data),n=4))
#'
#' test1b
#' test2b
#' test3b
#' test4b
#' test5b
#' test6b

kernelCItest <- function (  x,
                            y,
                            S=NULL,
                            suffStat,
                            verbose=FALSE,
                            data,
                            ic.method=NULL,
                            p=NULL,
                            index=NULL,
                            sig=NULL,
                            numCol=NULL,
                            numCluster=NULL,
                            eps=NULL,
                            paral=NULL
){

  x <- as.numeric(x)
  y <- as.numeric(y)

  if(is.null(suffStat$ic.method)) {suffStat$ic.method <- ic.method}
  if(is.null(suffStat$p)) {suffStat$p <- p}
  if(is.null(suffStat$index)) {suffStat$index <- index}
  if(is.null(suffStat$sig)) {suffStat$sig <- sig}
  if(is.null(suffStat$numCol)) {suffStat$numCol <- numCol}
  if(is.null(suffStat$numCluster)) {suffStat$numCluster <- numCluster}
  if(is.null(suffStat$eps)) {suffStat$eps <- eps}
  if(is.null(suffStat$paral)) {suffStat$paral <- paral}


  if (is.null(suffStat$data)){stop("Need to provide data")}
  n <- nrow(suffStat$data)

  if (is.null(suffStat$ic.method)){suffStat$ic.method<-"dcc.perm"; if(verbose){cat("Independence criterion method was not provided, using the default setting for ic.method distance covariance \n")}}

  if (is.null(suffStat$index) & suffStat$ic.method %in% c("dcc.perm","dcc.gamma") ){suffStat$index<-1; if(verbose){cat("The index for the distance covariance test was not provided, using the default setting index=",1,"\n")}}
  if (is.null(suffStat$p) & suffStat$ic.method %in% c("dcc.perm","hsic.perm","hsic.clust") ){suffStat$p<-100; if(verbose){cat("The number of permutations for the permutation test was not provided, using the default setting p=",100,"\n")}}
  if (is.null(suffStat$sig) & suffStat$ic.method %in% c("hsic.gamma","hsic.perm","hsic.clust") ){suffStat$sig<-1; if(verbose){cat("The kernel width sigma for HSIC not provided, using the default setting sig=",1,"\n")}}
  if (is.null(suffStat$numCol) & suffStat$ic.method %in% c("hsic.gamma","hsic.perm","hsic.clust","dcc.gamma") ){suffStat$numCol<-floor(n/10); if(verbose){cat("The number of columns for the Incomplete Cholesky decomposition was not provided, using the default setting numCol=",floor(n/10),"\n")}}
  if (is.null(suffStat$eps) & suffStat$ic.method %in% c("hsic.gamma","hsic.perm","hsic.clust") ){suffStat$eps<-0.1; if(verbose){cat("The normalizing constanst eps for HSIC cluster was not provided, using the default setting eps=",0.1,"\n")}}
  if (is.null(suffStat$numCluster) & suffStat$ic.method %in% c("hsic.clust") ){suffStat$numCluster<-floor(n/10); if(verbose){cat("The number of clusters for the HSIC cluster was not provided, using the default setting numCluster=",floor(n/10),"\n")}}
  if (is.null(suffStat$paral) & suffStat$ic.method %in% c("hsic.clust") ){suffStat$paral<-1; if(verbose){cat("The number of cores to use for parallel computation was not provided, using the default setting paral=",1,"\n")}}

  if (length(S)==0L){
    if (suffStat$ic.method=='dcc.perm'){pval = dcov.test(x=suffStat$data[,x], y=suffStat$data[,y], index = suffStat$index, R = suffStat$p)$p.value}
    else if (suffStat$ic.method=='dcc.gamma'){pval = dcov.gamma(x=suffStat$data[,x], y=suffStat$data[,y], index = suffStat$index, numCol=suffStat$numCol)$p.value}
    else if (suffStat$ic.method=='hsic.perm'){pval = hsic.perm(x=suffStat$data[,x], y=suffStat$data[,y], sig = suffStat$sig, p=suffStat$p, numCol=suffStat$numCol)$p.value}
    else if ((suffStat$ic.method=='hsic.gamma') || (suffStat$ic.method=='hsic.clust')){pval = hsic.gamma(x=suffStat$data[,x], y=suffStat$data[,y], sig=suffStat$sig, numCol=suffStat$numCol)$p.value}
  }
  else{
    if (suffStat$ic.method=='hsic.clust'){pval <- hsic.clust(x=suffStat$data[,x], y=suffStat$data[,y], z=suffStat$data[,S], sig=suffStat$sig, p=suffStat$p, numCluster=suffStat$numCluster, numCol=suffStat$numCol, eps=suffStat$eps, paral=suffStat$paral)$p.value}
    else {
      residuals <- regrXonS(suffStat$data[,c(x,y)], suffStat$data[,S])
      resx <- residuals[,1]
      resy <- residuals[,2]
      if (suffStat$ic.method=='dcc.perm'){pval = dcov.test(resx, resy, index = suffStat$index, R = suffStat$p)$p.value}
      else if (suffStat$ic.method=='dcc.gamma'){pval = dcov.gamma(x=resx, y=resy, index = suffStat$index, numCol=suffStat$numCol)$p.value}
      else if (suffStat$ic.method=='hsic.perm'){pval = hsic.perm(x=resx, y=resy, sig = suffStat$sig, p=suffStat$p, numCol=suffStat$numCol)$p.value}
      else if (suffStat$ic.method=='hsic.gamma'){pval = hsic.gamma(x=resx, y=resy, sig=suffStat$sig, numCol=suffStat$numCol)$p.value}
    }
  }
  return(pval)
}

