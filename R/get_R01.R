
#' Estimate the 0-1 Correlation matrix by multiple testing
#'@param X gene by individual matrix for estimating correlations
#'@param alpha FDR level for claiming correlations
#'@return a sparse matrix, 1 indicates correlated.
#'@import Rfast
#'@import Matrix
#'@export
#'@details
#'This function implements the one sample Large-Scale Correlation Tests With Normal Approximation
#'proposed by Cai and Liu(2016)

get_R01 = function(X,alpha=0.1){

  p = nrow(X)
  n = ncol(X)

  X = apply(X,2,function(z){z/sum(z)*1e6})


  X.center = scale(t(X),center=TRUE,scale=FALSE)
  S = cova(X.center,center = TRUE)

  # calc S2
  S2 = 0
  for(k in 1:n){
    S2 = S2+(tcrossprod(X.center[k,]))^2
  }

  Tmat = S*(n-1)/sqrt(S2+(2-n)*S^2)
  P = 2*(1-pnorm(abs(Tmat[lower.tri(Tmat)])))
  P.order = sort(P,decreasing = TRUE)



  nt = length(P.order)
  P.adj = P.order*(p^2-p)/2/(nt:1)

  for(t in 1:nt){
    if(P.adj[t]<=alpha){
      break
    }
  }

  bp = 2*(1-pnorm(sqrt(4*log(p)-2*log(log(p)))))
  if(P.order[t]<bp){
    thresh = 2*(1-pnorm(sqrt(4*log(p))))
  }else{
    thresh = P.order[t]
  }

  P.rej = c()
  P.rej[P>thresh] = 0
  P.rej[P<=thresh] = 1
  P.mat = Matrix(0,nrow=p,ncol=p,sparse = T)
  P.mat[lower.tri(P.mat)] = P.rej
  cor.idx = which(P.mat!=0,arr.ind = T)
  cor.idx = rbind(cor.idx,cbind(cor.idx[,2],cor.idx[,1]))


  if(length(cor.idx)==0){
    cor.idx = NULL
  }

  R01 = sparseMatrix(i=cor.idx[,1],j=cor.idx[,2],dims=c(G,G))
  diag(R01) = 1

  return(R01)

}









