


#' Main function performing MEAD
#'@param bulk bulk samples to be deconvolved. A count matrix, whose rownames are gene names and colnames are individual names.
#'@param ref reference dataset, a SingleCellExperiment object, with phenotype: cell_type, individual
#'@param cell_types cell types to be used in the deconvolution. If NULL, will use the cell types in the ref dataset.
#'@param preprocessing_control a list of control parameters passed to MEAD_preprocessing
#'@param estimation_control a list of control parameters passed to MEAD_est
#'@param R01 A 0-1 correlation matrix, from get_R01() function.
#'@return a list of:
#'  *p_hat: estimated cell type proportions
#'  *p_hat_se: standard errors of p_hat
#'  *p_group_diff_hat: estimated two group means different if groups labels are provided.(levels(as.factor(groups)))
#'@export
#'@import SingleCellExperiment

MEAD = function(bulk,
                ref,
                cell_types=NULL,
                preprocessing_control = list(),
                estimation_control = list(),
                R01 = NULL){

  pre_control = preprocessing_control_default()
  pre_control = modifyList(pre_control,preprocessing_control,keep.null = TRUE)
  datax = MEAD_preprocessing(bulk,
                             ref,
                             cell_types=cell_types,
                             marker_gene = pre_control$marker_gene,
                             gene_thresh=pre_control$gene_thresh,
                             max_count_quantile_celltype=pre_control$max_count_quantile_celltype,
                             max_count_quantile_indi = pre_control$max_count_quantile_indi,
                             filter.gene=pre_control$filter.gene,
                             R01=R01)

  ref_mat = MEAD_getX(datax$ref,
                      datax$cell_types,
                      datax$individuals)

  est_control = est_control_default()
  est_control = modifyList(est_control,estimation_control,keep.null = TRUE)
  res = MEAD_est(datax$bulk,
                ref_mat$X,
                ref_mat$V,
                ref_mat$w,
                R01=datax$R01,
                hc.type=est_control$hc.type,
                centeringXY = est_control$centeringXY,
                nfold=est_control$nfold,
                folds = est_control$folds,
                groups=est_control$groups,
                calc_var=est_control$calc_var,
                use.QP = est_control$use.QP)

  return(res)

}


#' Default MEAD_preprocessing input parameters
preprocessing_control_default = function(){
  list(marker_gene = NULL,
       gene_thresh=0.05,
       max_count_quantile_celltype=0.99,
       max_count_quantile_indi = 0.99,
       filter.gene=TRUE)
}

#' Default MEAD_estimation input parameters
est_control_default = function(){
  list(hc.type='hc3',
       centeringXY = FALSE,
       nfold=10,
       folds = NULL,
       groups=NULL,
       calc_var=TRUE,
       use.QP = FALSE)
}

#' Get the confidence interval of cell type proportions
#'@param mead.out output from running MEAD
#'@param alpha (1-alpha) confidence interval
#'@return a list of two matrices: ci_l and ci_r, corresponding to lower and upper bound of the confidence intervals.
#'@export
get_ci = function(mead.out,alpha = 0.05){
  ci_l = mead.out$p_hat - qnorm(1-alpha/2)*mead.out$p_hat_se
  ci_r = mead.out$p_hat + qnorm(1-alpha/2)*mead.out$p_hat_se
  ci_l[is.na(ci_l)] = 0
  ci_r[is.na(ci_r)] = 1
  ci_l = pmax(ci_l,0)
  ci_r = pmin(ci_r,1)
  return(list(ci_l = t(ci_l),ci_r=t(ci_r)))
}







