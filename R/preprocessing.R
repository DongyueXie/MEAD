
#' Prepare reference matrix and its variance for deconvolution
#'@param ref returned by MEAD_preprocessing
#'@param cell_types cell types to be considered
#'@param individuals subject in reference dataset to use
#'@return a list of:
#'  *X: reference matrix, G by K
#'  *V: variance of reference matrix, G by K^2
#'  *w: a vector of gene weights, length G
#'@export
#'@importFrom vashr vash
MEAD_getX = function(ref,cell_types,individuals){

  K = length(cell_types)
  NI = length(individuals)
  G = nrow(ref)
  X_array = array(dim = c(G,K,NI))

  for(i in 1:NI){
    refi = ref[,which(ref$individual==individuals[i])]
    for(k in 1:K){
      ctk = cell_types[k]
      X_array[,k,i] = rowMeans(apply(counts(refi)[,refi$cell_type==ctk,drop=FALSE],2,function(z){z/sum(z)}))
    }
  }

  X = apply(X_array,c(1,2),mean,na.rm=TRUE)
  V = t(apply(X_array,c(1),function(z){(cov(t(z),use = 'complete.obs'))}))/NI

  rownames(X) = rownames(ref)
  colnames(X) = cell_types
  rownames(V) = rownames(ref)


  fit.vash = vash(sqrt(rowSums(V)),df=NI-1)
  w = 1/(fit.vash$sd.post)^2
  w = w/sum(w)
  names(w) = rownames(ref)

  return(list(X=X,V=V,w=w))
}





#' Preprocess bulk and reference datasets
#'@param bulk bulk samples to be deconvolved. A count matrix, whose rownames are gene names and colnames are individual names.
#'@param ref reference dataset, a SingleCellExperiment object, with phenotype: cell_type, individual
#'@param cell_types cell types to be used in the deconvolution. If NULL, will use the cell types in the ref dataset.
#'@param gene_thresh a fraction. Remove genes that have expression in less than (gene_thresh * total numberof cells) cells.
#'@param marker_gene marker genes to use if any. Default is NULL.
#'@param max_count_quantile_celltype a fraction. Remove union of genes that expressed more than max_count_quantile_celltype in each cell type.
#'@param max_count_quantile_indi a fraction. Remove union of genes that expressed more than max_count_quantile_indi in each individual.
#'@param verbose whether print progress
#'@return a list of:
#'  *bulk: bulk samples
#'  *ref: ref samples
#'  *cell_types: cell_types to be considered
#'  *individuals: individuals in the reference dataset
#'@export
#'@import SingleCellExperiment
MEAD_preprocessing = function(bulk,ref,
                              marker_gene = NULL,
                              cell_types=NULL,
                              gene_thresh=0.05,
                              max_count_quantile_celltype=0.99,
                              max_count_quantile_indi = 0.99,
                              filter.gene=TRUE,
                              R01 = NULL){

  # remove cells without expression
  rm.cell = which(colSums(counts(ref))==0)
  if(length(rm.cell)!=0){
    ref = ref[,-rm.cell]
  }

  # what cell types to use
  if(is.null(cell_types)){
    if(is.factor(ref$cell_type)){
      cell_types = levels(ref$cell_type)
    }else{
      cell_types = levels(as.factor(ref$cell_type))
    }
  }
  K = length(cell_types)



  # pick up cells that match cell types we are interested in
  ref = ref[,which(ref$cell_type %in% cell_types)]

  ## drop levels of cell_type, and individuals
  if(is.factor(ref$cell_type)){
    ref$cell_type = droplevels(ref$cell_type)
  }
  if(is.factor(ref$individual)){
    ref$individual = droplevels(ref$individual)
  }

  # individuals to use
  if(is.factor(ref$individual)){
    individuals = levels(ref$individual)
  }else{
    individuals = levels(as.factor(ref$individual))
  }

  # filter out genes

  if(filter.gene){
    if(gene_thresh<1){
      gene_thresh = round(gene_thresh*ncol(ref))
    }
    rm.gene = which(rowSums(counts(ref)!=0)<gene_thresh)

    # remove union of genes that expressed more than max_count_quantile_celltype in each cell type

    if(!is.null(max_count_quantile_celltype)){

      rm.gene.high = c()
      for(k in 1:K){

        cell_k_idx = which(ref$cell_type==cell_types[k])
        if(length(cell_k_idx)!=0){
          gene_counts = rowSums(counts(ref[,cell_k_idx]))
          rm.gene.high = c(rm.gene.high,which(gene_counts>quantile(gene_counts,max_count_quantile_celltype)))
        }

      }

      rm.gene = unique(c(rm.gene,rm.gene.high))

    }

    if(!is.null(max_count_quantile_indi)){

      rm.gene.indi = c()

      for(j in 1:length(individuals)){

        indi_j_idx = which(ref$individual==individuals[j])
        if(length(indi_j_idx)!=0){
          gene_counts = rowSums(counts(ref[,indi_j_idx]))
          rm.gene.indi = c(rm.gene.indi,which(gene_counts>quantile(gene_counts,max_count_quantile_indi)))
        }

      }

      rm.gene = unique(c(rm.gene,rm.gene.indi))

    }

    # find ref genes to use

    if(length(rm.gene)!=0){
      gene_ref = rownames(ref)[-rm.gene]
    }
  }else{
    gene_ref = rownames(ref)
    message('Filtering ref genes recommended')
  }

  if(!is.null(R01)){
    if(is.null(rownames(R01))){
      stop('gene names must be provided for R01 matrix')
    }
    gene_R01 = rownames(R01)
    genes = intersect(intersect(rownames(bulk),gene_ref),gene_R01)
  }else{
    genes = intersect(rownames(bulk),gene_ref)
  }


  if(!is.null(marker_gene)){
    genes = intersect(genes,marker_gene)
  }
  if(length(genes)==0){
    stop('No common genes found. Check gene names.')
  }



  return(list(bulk = bulk[match(genes,rownames(bulk)),],
              ref = ref[match(genes,rownames(ref)),],
              cell_types = cell_types,
              individuals = individuals,
              R01 = R01[match(genes,rownames(R01)),match(genes,rownames(R01))]))

}

