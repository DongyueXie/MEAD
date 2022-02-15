#'@title Function performs estimation and inference.
#'@param y a vector of bulk sample
#'@param X reference matrix, estimated relative expression, gene by cell
#'@param Vg variance matrix of X, either G by K or G by K^2.
#'@param w gene weights; if null, w = 1
#'@param hc.type 'hc3', 'hc2', 'hc0', 'cv' or 'cv_indep'
#'@param nfold if using cv, the number of fold to be used. If folds not provided, will run k-medoids on 11' - R01 to determine folds.
#'@param folds if using cv, the folds to be used. folds looks like 111222333.
#'@param R01 the 0-1 correlation matrix
#'@param calc_var whether calculate the covariance matrices.
#'@param use.QP whether use quadratic programming to impose nonnegative constraint
#'@return a list including: 1. p_hat: estimated cell type proportions; 2. p_hat_se: their standard errors; 3. p_group_diff_hat: if groups provided, estimated group mean difference
#'@export

MEAD_est = function(y,X,Vg,
                    w=NULL,
                    hc.type='hc3',
                    centeringXY = FALSE,
                    nfold=10,
                    folds = NULL,
                    R01=NULL,
                    groups=NULL,
                    calc_var=TRUE,
                    use.QP = FALSE
){

  G = nrow(X)
  K = ncol(X)

  # nb is the number of bulk sample
  nb = ncol(y)
  if(is.null(nb)){
    nb=1
  }

  # weights

  if(is.null(w)){
    w = 1
  }
  if(length(w)==1){
    w = rep(w,nrow(X))
  }
  w = w/sum(w)*G

  Xw = X*sqrt(w)
  yw = cbind(y*sqrt(w))
  if(centeringXY){
    yw = apply(yw,2,scale,scale=FALSE)
    Xw = apply(Xw,2,scale,scale=FALSE)
  }
  Vw = Vg*w
  A = crossprod(Xw)
  d_Vg = ncol(Vg)
  if(d_Vg==K^2){
    V = matrix(c(colSums(Vw)),ncol = K)
  }else if(d_Vg==K){
    V = diag(c(colSums(Vw)))
  }else{
    stop('check dimension of Vg')
  }

  #browser()
  # if(verbose){
  #   message("estimating proportions")
  # }
  Q = (A-V)

  Qinv = solve(Q)


  if(use.QP){
    beta_tilde_hat = matrix(nrow=K,ncol=nb)
    for(b in 1:nb){
      QP.out = quadprog::solve.QP(Q,t(Xw)%*%yw[,b],diag(K),rep(0,K))
      beta_tilde_hat[,b] = QP.out$solution
    }

  }else{
    beta_tilde_hat = pmax(Qinv%*%t(Xw)%*%yw,0)
  }


  if(calc_var){

    Q_inv = kronecker(diag(nb),Qinv)
    h = rowSums((Xw%*%Qinv)*Xw)


    Sigma = get_SIGMA(y=yw,
                      X=Xw,
                      beta=beta_tilde_hat,
                      V=Vw,
                      h=h,
                      nb=nb,
                      G=G,
                      K=K,
                      hc.type=hc.type,
                      nfold=nfold,
                      folds=folds,
                      R01=R01)
    covb = Q_inv%*%Sigma%*%Q_inv


    # delta method
    # covb is a (nb*K)*(nb*K) cov matrix
    # formulate Jacobian matrix

    # if(verbose){
    #   message("performing delta method")
    # }

    beta_tilde_hat = cbind(beta_tilde_hat)
    rownames(beta_tilde_hat) = colnames(X)
    J = matrix(0,nrow=(nb*K),ncol=(nb*K))

    for(i in 1:nb){
      J[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = J_sum2one(beta_tilde_hat[,i],K)
    }

    asyV = (J)%*%covb%*%t(J)

    p_hat = apply(beta_tilde_hat,2,function(z){z/sum(z)})
    p_hat_se = sqrt(diag(asyV))
    p_hat_se = matrix(p_hat_se,ncol=nb)
    rownames(p_hat) = colnames(X)
    colnames(p_hat) = colnames(y)
    rownames(p_hat_se) = colnames(X)
    colnames(p_hat_se) = colnames(y)

    if(!is.null(groups)){
      two_group_res = get_two_group_diff(groups,p_hat)
    }else{
      two_group_res = NULL
    }

    # if(verbose){
    #   message("done")
    # }

    return(list(p_hat=p_hat,
                p_hat_se=p_hat_se,
                p_group_diff_hat=two_group_res
    ))

  }else{

    p_hat = apply(beta_tilde_hat,2,function(z){z/sum(z)})
    if(!is.null(groups)){
      two_group_res = get_two_group_diff(groups,p_hat)
    }else{
      two_group_res = NULL
    }
    return(list(p_hat=p_hat,
                p_group_diff_hat=two_group_res
    ))

  }



}


#'@title Jacobian matrix of sum-to-1 scale function
#'@param b beta_hat
#'@param K number of cell types
#'@return Jacobian matrix
J_sum2one = function(b,K){
  J = - (b)%*%t(rep(1,K))
  diag(J) = diag(J) + sum(b)
  J = J/sum(b)^2
  J
}


#'@title return the estimated group mean difference
get_two_group_diff = function(groups,p_hat){
  nb = length(groups)
  if(is.factor(groups)){
    group_name = levels(groups)
  }else{
    group_name = levels(as.factor(groups))
  }

  group1_idx = which(groups==group_name[1])
  a = c()
  a[group1_idx] = 1/length(group1_idx)

  group2_idx = which(groups==group_name[2])
  a[group2_idx] = -1/length(group2_idx)

  diff_group = rowMeans(p_hat[,group1_idx,drop=F]) - rowMeans(p_hat[,group2_idx,drop=F])
  return(diff_group)
}

#' Calculate the Sigma in variance
#'@importFrom cluster pam
get_SIGMA = function(y,X,beta,V,h,nb,G,K,hc.type,nfold,folds,R01){

  res = y - X%*%beta
  h = pmax(pmin(abs(h),1-1/G),0)
  if(hc.type == 'hc0'){
    res.hc = res
  }else if(hc.type == 'hc2'){
    res.hc = res/sqrt(1-h)
  }else if(hc.type == 'hc3'){
    res.hc = res/(1-h)
  }else if(hc.type%in%c('cv','cv_indep')){

    if(is.null(folds)){
      if(is.null(R01)){
        stop('R01 should be provided to get folds when using cv based method')
      }else{
        folds = pam(as.dist(1-R01),nfold,pamonce=5)$clustering
      }
    }

    if(hc.type=='cv'){
      res.hc = get_cv_res(y,X,V,folds=folds)
    }else if(hc.type=='cv_indep'){
      res.hc = get_cv_res_indep(y,X,V,folds=folds,R01=R01)
    }

  }else{
    stop('Invalid hc.type')
  }



  d_V = ncol(V)
  # formulate score matrix, G by nk
  score_mat = matrix(nrow=G,ncol=K*nb)
  for(i in 1:nb){

    if(d_V==K^2){
      #Vbi = t(apply(V,1,function(z){v = matrix(z,ncol=K);v%*%beta[,i]}))
      Vbi = V%*%matrix(c(rep(c(beta[,i],rep(0,K^2)),K-1),beta[,i]),ncol=K)
    }
    if(d_V==K){
      Vbi = (V)%*%diag(c(beta[,i]))
    }
    score_i = (c(res.hc[,i])*X+Vbi)
    score_mat[,((i-1)*K+1):(i*K)] = score_i

  }

  # if(verbose){
  #   message("calculating covariance matrix")
  # }

  if(!is.null(R01)){
    Sigma = crossprod(score_mat,R01)%*%score_mat
  }else{
    Sigma = crossprod(score_mat)
  }


  Sigma = (Sigma+t(Sigma))/2
  return(Sigma)


}




#'@title calculate cv residuals
get_cv_res = function(y,X,V,folds=NULL){
  n = nrow(y)
  n_bulk = ncol(y)
  K = ncol(X)

  d_Vg = ncol(V)

  #browser()


  res = matrix(nrow=n,ncol=n_bulk)
  if(is.null(folds)){
    stop('Folds should be provided')
  }
  nfold = length(table(folds))

  # for the genes in each fold, estimate the beta hat use the genes in the rest folds
  for(f in 1:nfold){
    idx = which(folds==f)
    X.temp = X[-idx,]
    y.temp = y[-idx,]
    if(d_Vg==K^2){
      V.temp = matrix(c(colSums(V[-idx,])),ncol = K)
    }else if(d_Vg==K){
      V.temp = diag(c(colSums(V[-idx,])))
    }else{
      stop('check dimension of V')
    }

    bhat = pmax(solve(crossprod(X.temp)-V.temp)%*%t(X.temp)%*%y.temp,0)
    res[idx,] = y[idx,] - X[idx,]%*%bhat
  }

  res

}


#'@title calculate cv+indep gene residuals
get_cv_res_indep = function(y,X,V,folds=NULL,R01=NULL){
  n = nrow(y)
  n_bulk = ncol(y)
  K = ncol(X)

  d_Vg = ncol(V)

  #browser()


  res = matrix(nrow=n,ncol=n_bulk)
  if(is.null(folds)){
    stop('Folds should be provided')
  }
  if(is.null(R01)){
    stop('R01 should be provided when using cv_indep')
  }

  nfold = length(table(folds))



  for(f in 1:nfold){
    idx = which(folds==f)

    temp = colSums(R01[idx, ])
    #index = ((1:n)[-idx])[which(temp == 0)]
    index = which(temp==0)
    #print(length(index))
    if(length(index)<(n/2)){
      index = ((1:n)[-idx])[(order(temp[-idx],decreasing = F))[1:((n-length(idx))/2)]]
    }

    X.temp = X[index,]
    y.temp = y[index,]
    if(d_Vg==K^2){
      V.temp = matrix(c(colSums(V[index,])),ncol = K)
    }else if(d_Vg==K){
      V.temp = diag(c(colSums(V[index,])))
    }else{
      stop('check dimension of V')
    }

    bhat = pmax(solve(crossprod(X.temp)-V.temp)%*%t(X.temp)%*%y.temp,0)
    res[idx,] = y[idx,] - X[idx,]%*%bhat
  }

  res

}
