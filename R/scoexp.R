#' Radial basis function kernel
#'
#' @param dis_mat Distance matrix
#' @param sigm Width of rbfk
#' @param zero_diag
#'
#' @return
#' @export
#'
#' @examples rbfk(dis_mat, sigm, zero_diag=F)
rbfk <- function (dis_mat, sigm, zero_diag=T) {
  rbfk_out <- exp(-1*(dis_mat^2)/(2*sigm^2))
  if (zero_diag) diag(rbfk_out) <- 0
  return(rbfk_out)
}

#' Title
#'
#' @param X Expression matrix, n X p
#' @param W Weight matrix, n X n
#' @param method Correlation method, pearson or spearman
#' @param na_zero Na to zero
#'
#' @return Weighted correlation matrix, p X p
#' @export
#'
#' @examples wcor(X=expr_test, W=rbfk_out, method='spearman')
wcor <- function(X, W, method=c('pearson', 'spearman')[1], na_zero=T) {
  if (method=='spearman') X <- apply(X, 2, rank)
  X <- scale(X)
  X[is.nan(X)] <- NA
  W_cov_temp <- (t(X) %*% W %*% X)
  W_diag_mat <- sqrt(diag(W_cov_temp) %*% t(diag(W_cov_temp)))
  cor_mat <- W_cov_temp/W_diag_mat
  if (na_zero) cor_mat[which(is.na(cor_mat), arr.ind=T)] <- 0
  cor_mat
}

#' Remove low correlation genes until reaching threshold
#'
#' @param cor_mat
#' @param ave_cor_cut
#' @param min_n
#' @param max_n
#' @param na_diag
#'
#' @return
#' @export
#'
#' @examples cor_remove(cor_mat, ave_cor_cut = 0.5, min_n=10, max_n=100, na_diag=F)
cor_remove <- function (cor_mat, ave_cor_cut = 0.5, min_n=5, max_n=100, na_diag=F) {
  if (na_diag) {
    diag(cor_mat) <- NA
  }
  cor_mat_temp <- cor_mat
  cor_ave <- rowMeans(cor_mat_temp, na.rm = T)
  cor_max <- max(cor_mat_temp, na.rm = T)
  if ((cor_max < ave_cor_cut)|(nrow(cor_mat_temp)<min_n)) {
    cor_mat_temp <- matrix(NA, 1, 1)
  } else {
    cor_min <- min(cor_ave)
    cor_min_idx <- which.min(cor_ave)
    idx <- 1
    while ((cor_min < ave_cor_cut & idx <= (nrow(cor_mat)-min_n))|((nrow(cor_mat)-idx))>=max_n) {
      cor_mat_temp <- cor_mat_temp[-cor_min_idx, -cor_min_idx]
      cor_ave <- rowMeans(cor_mat_temp, na.rm = T)
      cor_min <- min(cor_ave)
      cor_min_idx <- which.min(cor_ave)
      idx <- idx + 1
    }
  }
  return (cor_mat_temp)
}

#' A wrapper function of ConensusClusterPlus with some default parameters
#'
#' @param d
#' @param maxK
#' @param reps
#' @param distance
#' @param verbose
#' @param ...
#'
#' @return
#' @export
#'
#' @examples cc_wrapper(d, maxK=8, reps=20, distance='spearman', verbose=T)
cc_wrapper <- function(d, maxK=8, reps=20, distance='spearman', verbose=F, plot='png', ...) {
  cc_output <- ConsensusClusterPlus::ConsensusClusterPlus(d, maxK=maxK, reps=reps, distance=distance, verbose=verbose, plot=plot, ...)
  return(cc_output)
}

#' Consensus gene module selection
#'
#' @param cc_res
#' @param cor_mat
#' @param k
#' @param avg_con_min
#' @param avg_cor_min
#' @param min_gen
#' @param max_gen
#'
#' @return
#' @export
#'
#' @examples cc_gene_k(cc_res, cor_mat, k=8, avg_con_min=.5, avg_cor_min=.5, min_gen=20, max_gen=100)
cc_gene_k <- function(cc_res, cor_mat, k=8, avg_con_min=.5, avg_cor_min=.5, min_gen=10, max_gen=100) {
  cc_res_k <- cc_res[[k]]
  con_k <- cc_res_k$consensusMatrix %>% magrittr::set_rownames(names(cc_res_k$consensusClass)) %>% magrittr::set_colnames(names(cc_res_k$consensusClass))
  res <- list()
  for (k_i in 1:k) {
    cat(k_i, ' ')
    gene_k <- names(cc_res_k$consensusClass)[cc_res_k$consensusClass==k_i]
    cor_mat_temp <- cor_mat[gene_k, gene_k, drop=F]
    cor_mat_temp <- cor_remove(cor_mat_temp, max_n=1000, min_n=1, ave_cor_cut=avg_cor_min)
    gene_k <- rownames(cor_mat_temp)
    if (is.null(gene_k)) gene_k <- NA
    con_mat_temp <- con_k[gene_k, gene_k, drop=F]
    con_mat_temp <- cor_remove(con_mat_temp, max_n=max_gen, min_n=1, ave_cor_cut=avg_con_min)
    gene_k <- rownames(con_mat_temp)
    if (is.null(gene_k)) gene_k <- NA
    res[[k_i]] <- gene_k
  }
  names(res) <- paste0('k', 1:k)
  res <- res[sapply(res, function(x) length(x)>=min_gen)]
  res
}

#' Title
#'
#' @param sim_mat similarity matrix
#' @param powerVector power vector
#' @param minClusterSize
#' @param ...
#'
#' @return
#' @export
#'
#' @examples wgcna_wrapper(sim_mat, powerVector=c(1:20), minClusterSize=50)
wgcna_wrapper <- function(sim_mat, powerVector=c(1:20), minClusterSize=50, ...) {
  res <- list()
  ## map to [0, 1] and symmetric ##
  sim_mat <- sim_mat - min(sim_mat, na.rm = T)
  sim_mat <- sim_mat/max(sim_mat, na.rm = T)
  sim_mat <- (sim_mat + t(sim_mat))/2
  sft_thresh <- WGCNA::pickSoftThreshold.fromSimilarity(sim_mat, powerVector=powerVector)
  pw <- sft_thresh$powerEstimate
  res$sft_thresh <- sft_thresh
  disTom <- 1 - (WGCNA::adjacency.fromSimilarity(sim_mat, power=pw) %>% WGCNA::TOMsimilarity(.))
  geneTree <- fastcluster::hclust(as.dist(disTom), method = 'average')
  dynamicMods_cut <- dynamicTreeCut::cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize=minClusterSize)
  gene_vec <- dynamicMods_cut %>% purrr::set_names(rownames(sim_mat))
  gene_vec <- gene_vec[gene_vec!=0]
  res$gv <- gene_vec
  return(res)
}


#' Title
#'
#' @param wgcna_res
#' @param cor_mat
#' @param avg_cor_min
#' @param min_gen
#' @param max_gen
#'
#' @return
#' @export
#'
#' @examples wgcna_gene_k(wgcna_res, cor_mat, avg_cor_min=.5, min_gen=10, max_gen=100)
wgcna_gene_k <- function(wgcna_res, cor_mat, avg_cor_min=.5, min_gen=10, max_gen=100) {
  res <- list()
  gene_vec <- wgcna_res$gv
  for (k_i in sort(unique(gene_vec))) {
    cat(k_i, ' ')
    gene_k <- names(gene_vec)[gene_vec==k_i]
    cor_mat_temp <- cor_mat[gene_k, gene_k, drop=F]
    cor_mat_temp <- cor_remove(cor_mat_temp, max_n=max_gen, min_n=1, ave_cor_cut=avg_cor_min)
    gene_k <- rownames(cor_mat_temp)
    if (is.null(gene_k)) gene_k <- NA
    res[[k_i]] <- gene_k
  }
  names(res) <- paste0('k', sort(unique(gene_vec)))
  res <- res[sapply(res, function(x) length(x)>=min_gen)]
  res
}


#' Title
#'
#' @param schart_inp SChart input on cell of interests
#' @param sigm
#' @param assay
#' @param gene_select
#' @param zero_cutoff
#' @param cor_method
#' @param approach Which approach to use? consensus clustering (cc) or weighted correlation network analysis (wgcna)
#' @param maxK
#' @param k
#' @param avg_con_min
#' @param avg_cor_min
#' @param min_gen
#' @param max_gen
#' @param keep_cc If TRUE, keep the cc model
#' @param keep_wgcna If TRUE, keep the wgcna model
#' @param keep_kern
#' @param keep_wcor
#' @param ...
#'
#' @return
#' @export
#'
#' @examples scoexp(schart_inp, sigm=NULL, assay='RNA', gene_select=NULL, zero_cutoff=5, cor_method='spearman', approach=c('cc', 'wgcna')[1], maxK=8, k=8, avg_con_min=.5, avg_cor_min=.5, min_gen=20, max_gen=100, keep_cc=T, keep_wgcna=T, keep_kern=T, keep_wcor=T)
scoexp <- function(schart_inp, sigm=NULL, assay='RNA', gene_select=NULL, zero_cutoff=5, cor_method='spearman', approach=c('cc', 'wgcna')[1], maxK=8, k=8, avg_con_min=.5, avg_cor_min=.5, min_gen=20, max_gen=100, keep_cc=T, keep_wgcna=T, keep_kern=T, keep_wcor=T, ...) {
  if (!all(c('coord_x', 'coord_y') %in% colnames(schart_inp@meta.data))) stop('coord_x and coord_y not detected in the metadata')
  if (is.null(sigm)) sigm <- schart_inp@images[[1]]@scale.factors$spot_dis
  if (is.null(gene_select)) {
    cat('gene filtering...\n')
    feature_nz <- apply(schart_inp[[assay]]@data, 1, function(x) mean(x!=0)*100)
    features <- names(feature_nz)[feature_nz > zero_cutoff]
    cat(length(features), 'features after filtering...\n')
  } else if (length(gene_select) > 1) {
    features <- intersect(gene_select, rownames(schart_inp[[assay]]@data))
    if (length(features)==0) stop('No genes in gene_select detected')
  }
  schart_inp <- Seurat::ScaleData(schart_inp, features=features)
  res <- list(gs=c(), cc=c(), rbfk=c(), wcor=c())
  dist_mat <- dist(schart_inp@meta.data[, c('coord_x', 'coord_y')]) %>% as.matrix
  kern_mat <- rbfk(dist_mat, sigm=sigm, zero_diag=F)
  expr_mat <- t(as.matrix(schart_inp[[assay]]@scale.data))
  cat('Calculating spatial-weighted cross-correlation...\n')
  wcor_mat <- wcor(X=expr_mat, W=kern_mat, method=cor_method)

  if (approach=='cc') {
    wcor_dis <- as.dist(1-wcor_mat)
    cat('Consensus clustering...\n')
    cc_res <- cc_wrapper(d=wcor_dis, maxK=maxK, ...)
    cons_mat <- cc_res[[k]]$consensusMatrix %>% data.frame
    colnames(cons_mat) <- rownames(cons_mat) <- rownames(wcor_mat)
    cat('Gene module detecting...\n')
    K_gl <- cc_gene_k(cc_res=cc_res, cor_mat=wcor_mat, k=k, avg_con_min=avg_con_min, avg_cor_min=avg_cor_min, min_gen=min_gen, max_gen=max_gen)
    res$gs <- K_gl
    if (keep_cc) res$cc <- cc_res
    if (keep_kern) res$rbfk <- kern_mat
    if (keep_wcor) res$wcor <- wcor_mat
  } else if (approach=='wgcna') {
    cat('WGCNA...\n')
    wgcna_res <- wgcna_wrapper(sim_mat=wcor_mat, minClusterSize=min_gen, ...)
    K_gl <- wgcna_gene_k(wgcna_res, cor_mat=wcor_mat, avg_cor_min=avg_cor_min, min_gen=min_gen, max_gen=max_gen)
    res$gs <- K_gl
    if (keep_wgcna) res$wgcna <- wgcna_res
    if (keep_kern) res$rbfk <- kern_mat
    if (keep_wcor) res$wcor <- wcor_mat
  }
  return(res)
}

