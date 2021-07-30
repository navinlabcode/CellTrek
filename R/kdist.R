#' Title
#'
#' @param inp_df inp_df must contain cell_names, 'coord_x', 'coord_y' columns
#' @param ref Reference groups
#' @param ref_type 'all' or 'each'
#' @param que Query groups
#' @param k The number of nearest neighbors
#' @param new_name New name for kdist
#' @param keep_nn Keep Nearest Neighboor id matrix?
#'
#' @return A list of 1. kdist data frame and 2. a list of knn id matrix
#' @export
#'
#' @examples kdist_out <- kdist(inp_df=test_df, ref=c('A', 'B', 'C'), ref_type='each', que=unique(test_df$cell_names), k=10, keep_nn=T)
kdist <- function(inp_df, ref=NULL, ref_type='all', que=NULL, k=10, new_name='kdist', keep_nn=F) {
  ## Check ##
  if (!all(c('cell_names', 'coord_x', 'coord_y') %in% colnames(inp_df))) {stop("Input data must contain 'cell_names', 'coord_x', 'coord_y' columns")}
  if (any(!c(ref %in% inp_df$cell_names, que %in% inp_df$cell_names))) {stop('Reference or query group not in cell_names')}

  que_dat <- inp_df[inp_df$cell_names %in% que, c('coord_x', 'coord_y')]
  kNN_dist_df <- data.frame(matrix(NA, nrow=nrow(que_dat), ncol=0))
  kNN_nn_list <- list()
  if (ref_type=='each' & length(ref)>0) {
    for (i in 1:length(ref)) {
      ref_i <- ref[i]
      ref_dat <- inp_df[inp_df$cell_names %in% ref_i, c('coord_x', 'coord_y')]
      ## error when k >= nrow(ref_dat), reset k ##
      if (nrow(ref_dat) <= k) {k <- (nrow(ref_dat)-1)}

      kNN_res <- dbscan::kNN(x=ref_dat, k=k, query=que_dat)
      kNN_dist <- apply(kNN_res$dist, 1, mean)
      kNN_dist_df <- base::cbind(kNN_dist_df, kNN_dist)
      if (keep_nn) {
        nn_mat <- kNN_res$id
        nn_mat[] <- rownames(ref_dat)[c(nn_mat)]
        kNN_nn_list[[i]] <- nn_mat
      } else {
        kNN_nn_list[[i]] <- matrix(NA, 0, 0)
      }
    }
    colnames(kNN_dist_df) <- paste0(ref, '_kdist')
    names(kNN_nn_list) <- paste0(ref, '_ref')
  } else {
    ref_dat <- inp_df[inp_df$cell_names %in% ref, c('coord_x', 'coord_y')]
    kNN_res <- dbscan::kNN(x=ref_dat, k=k, query=que_dat)
    kNN_dist_df <- data.frame(apply(kNN_res$dist, 1, mean)) %>% magrittr::set_colnames(new_name)
    if (keep_nn) {
      nn_mat <- kNN_res$id
      nn_mat[] <- rownames(ref_dat)[c(nn_mat)]
      kNN_nn_list$nn_ref <- nn_mat
    } else {
      kNN_nn_list[[1]] <- matrix(NA, 0, 0)
    }
  }
  output <- list(kdist_df=kNN_dist_df, knn_list=kNN_nn_list)
  return(output)
}

#' Title
#'
#' @param schart_inp SChart seurat input
#' @param grp_col Column name in meta data for reference and query groups
#' @param ref Reference groups
#' @param ref_type 'all' or 'each'
#' @param que Query groups
#' @param k The number of nearest neighbors
#' @param new_name New name for kdist
#' @param keep_nn Keep Nearest Neighboor id matrix?
#'
#' @return SChart seurat output
#' @export
#'
#' @examples schart_test <- run_kdist(schart_inp=schart_inp, grp_col='cell_type', ref=c('A', 'B'), ref_type='each', que=unique(test_df$cell_names), k=10)
run_kdist <- function(schart_inp, grp_col='cell_type', ref=NULL, ref_type='all', que=NULL, k=10, new_name='kdist', keep_nn=F) {
  ## Check ##
  if (any(!c(grp_col, 'coord_x', 'coord_y') %in% colnames(schart_inp@meta.data))) {stop("schart_inp metadata must contain grp_col, 'coord_x', 'coord_y' columns")}

  inp_df <- schart_inp@meta.data %>% dplyr::select(cell_names=dplyr::one_of(grp_col), coord_x, coord_y)
  output <- kdist(inp_df=inp_df, ref=ref, ref_type=ref_type, que = que, k=k, new_name=new_name, keep_nn=keep_nn)
  schart_out <- Seurat::AddMetaData(schart_inp, metadata=output$kdist_df)
  return (schart_out)
}
