#' Title
#'
#' @param data_inp Cell-coordinates matrix
#' @param col_cell Column name of cell type, cell type names must be syntactically valid
#'
#' @return dummy data frame with cells on columns
#'
#' @examples cell_dummy_df <- as_dummy_df(celltrek_df, col_cell='Cell')
as_dummy_df <- function(data_inp, col_cell='cell_names') {
  cell_names <- sort(as.character(unique(data_inp[, col_cell])))
  data_inp$id_add <- rownames(data_inp)
  data_inp$col_add <- data_inp[, col_cell]
  data_out <- data_inp %>% dplyr::select(id_add, col_add) %>% dplyr::mutate(value=1) %>%
    tidyr::pivot_wider(id_cols=id_add, names_from=col_add, values_fill=list(value=0)) %>%
    data.frame(check.names=F)
  data_out <- data_out[, c('id_add', cell_names)]
  ## Reorder cell names ##
  data_out <- dplyr::left_join(data_inp, data_out, by='id_add') %>%
    set_rownames(.$id_add) %>% dplyr::select(-c(id_add, col_add))
  return(data_out)
}

#' Title
#'
#' @param x1
#' @param x2
#'
#' @return
#'
#' @examples
euc_dist <- function(x1, x2) {sqrt(sum((x1 - x2) ^ 2))}

#' Title
#'
#' @param data Cell type dummy table
#' @param coord Coordinates data with X and Y columns
#' @param min_num For cells number < min_num, add 1 till min_num
#' @param h Bandwidths for x and y directions, more details in kde2d function in MASS package
#' @param n Number of grid points in each directions, more details in kde2d function in MASS package
#' @param tot_norm Normalization by total?
#' @param Xlim The limits of X-axis
#' @param Ylim The limits of Y-axis
#'
#' @return Kenerl density on grid
#'
#' @examples kern_res <- sp_grid_kern_bin(data=cell_dummy_df, coord=coord_df, h=150, n=25, tot_norm=TRUE, Xlim=range(coord_df$X), Ylim=range(coor_df$Y))
sp_grid_kern_bin <- function(data, coord, min_num=15, h, n=25, tot_norm=TRUE, Xlim=range(coord[, 1]), Ylim=range(coord[, 2])) {

  ## For cells with less than min_num, randomly add 1 ##
  data_ <- data
  if (min_num<0) min_num <- 1
  col_min <- which(colSums(data) < min_num)
  if (length(col_min) > 0) {
    for (i in 1:length(col_min)) {
      set.seed(i)
      data_[, col_min[i]][sample(which(data_[, col_min[i]]==0), min_num-sum(data_[, col_min[i]]))] <- 1
    }
  }


  # data_out_zero <- matrix(0, n^2, length(col_rm)) %>% data.frame %>% set_colnames(names(col_rm))
  # if (length(col_rm)>0) {
  #   data_ <- data[, -col_rm]
  # } else {data_ <- data}

  colnames(coord) <- paste0('X', 1:length(colnames(coord)))
  data_merge <- merge(coord, data_, by='row.names')

  data_out_list <- data_merge %>% dplyr::select(colnames(data_)) %>%
    apply(2, function(x) {
      k2d_temp <- MASS::kde2d(data_merge$X1[which(x==1)], data_merge$X2[which(x==1)], lims=c(Xlim, Ylim), h, n) %>%
        extract2('z') %>% reshape2::melt %>% data.frame %>% set_colnames(c('grid_1', 'grid_2', 'Val'))
      if (tot_norm) {k2d_temp$Val <- k2d_temp$Val/sum(k2d_temp$Val)}
      k2d_temp
    })

  ## Rename the Val column ##
  data_out_list <- lapply(1:length(data_out_list), function(i) {
    colnames(data_out_list[[i]])[3] <- names(data_out_list)[i]
    data_out_list[[i]]
  })

  data_out <- Reduce(function(x, y) merge(x, y, by=c('grid_1', 'grid_2'), suffixes=c(colnames(x)[3], names(y))), data_out_list)
  data_out <- data_out[, c('grid_1', 'grid_2', colnames(data))]
  ## Keep the same cell order with input data ##
  data_out
}


#' Title
#'
#' @param X Density matrix
#' @param eps Small value added
#'
#' @return KL-divergence matrix
#'
#' @examples KL_res <- KL_(X=X, eps=1e-20)
KL_ <- function(X, eps=1e-20) {
  ## Need Suppress print-outs ##
  X_ <- X
  X_[X_<eps] <- eps
  X_ <- scale(X_, center=F, scale=colSums(X_)) %>% data.frame
  KL_D <- suppressMessages(philentropy::KL(t(X_), unit='log'))
  rownames(KL_D) <- colnames(KL_D) <- colnames(X_)
  KL_D[is.infinite(KL_D)] <- 650
  return(KL_D)
}


#' Title
#'
#' @param dummy_df Cell type dummy df
#' @param coord_df Coordinates df
#' @param min_num For cells number < min_num, add 1 till min_num
#' @param Xlim The limits of X-axis
#' @param Ylim The limits of Y-axis
#' @param boot_n Number of bootstrapping iterations
#' @param prop Subsample proportion
#' @param h Bandwidths for x and y directions, more details in kde2d function in MASS package
#' @param n Number of grid points in each directions, more details in kde2d function in MASS package
#' @param tot_norm Normalization by total?
#' @param eps Small value when calculate KL-divergence
#' @param replace should sampling be with replacement?
#'
#' @return A list of 1. Bootstrap KL-divergence; 2. MST consensus matrix
#'
#' @examples boot_mst_res <- KL_boot_mst(dummy_df=cell_type_dummy, coord_df=range(coord_df$X), Xlim=range(coord_df$Y), boot_n=100, prop=0.8, h=150, n=25, tot_norm=T, eps=1e-20)
KL_boot_mst <- function(dummy_df, coord_df, min_num=15, Xlim=range(coord_df[, 1]), Ylim=range(coord_df[, 2]), boot_n=100, prop=0.8, replace=T, h=150, n=25, tot_norm=T, eps=1e-20) {
  n_smp <- nrow(dummy_df)
  n_typ <- ncol(dummy_df)

  dis_boot_array <- array(NA, dim=c(n_typ, n_typ, boot_n),
                          dimnames=list(colnames(dummy_df), colnames(dummy_df), paste0('X_', c(1:boot_n))))
  mst_cons <- matrix(0, n_typ, n_typ) %>% data.frame %>% magrittr::set_rownames(colnames(dummy_df)) %>% magrittr::set_colnames(colnames(dummy_df))
  for (i in 1:boot_n) {
    cat(i, ' ')
    set.seed(i)
    ## Bootstrap ##
    idx <- sample(1:n_smp, round(n_smp*prop), replace=replace)
    data_boot <- dummy_df[idx, ]
    coord_boot <- coord_df[idx, ]
    k2d_boot <- sp_grid_kern_bin(data=data_boot, coord=coord_boot, min_num=min_num, Xlim=Xlim, Ylim=Ylim, h=h, n=n, tot_norm=tot_norm)
    dis_boot <- k2d_boot[, -c(1, 2)] %>% KL_(., eps=eps)
    # dis_boot <- k2d_boot[, -c(1, 2)] %>% KL_(., eps=1e-20)
    ##
    # dis_boot <- 1-k2d_boot[, -c(1, 2)] %>% cor(method='spearman')
    ##
    dis_boot_array[, , i] <- dis_boot
    ## MST ##
    mst_boot <- igraph::graph.adjacency(dis_boot, weighted=T, mode='upper') %>% igraph::mst(., algorithm='prim') %>%
      igraph::as_adjacency_matrix() %>% as.data.frame(make.names=F)
    mst_boot <- mst_boot[colnames(dummy_df), ]
    mst_boot <- mst_boot[, colnames(dummy_df)]
    mst_cons <- mst_cons + mst_boot/boot_n
  }
  output <- list(boot_array=dis_boot_array, mst_cons=mst_cons)
  return(output)
}

#' Title Build delaunay triangulation network
#'
#' @param coord_df Coordinates df
#' @param return_name Return edge list with names?
#' @param dist_cutoff Remove edges with distance >= dist_cutoff
#'
#' @return Delaunay triangulation network edge list
#' @export
#'
#' @examples test <- build_delaunayn(coord_df, return_name=T)
build_delaunayn <- function(coord_df, return_name=T, dist_cutoff=NULL) {
  delaunayn_res <- geometry::delaunayn(coord_df)
  delaunayn_el <- rbind(delaunayn_res[, c(1, 2)], delaunayn_res[, c(1, 3)], delaunayn_res[, c(2, 3)],
                        delaunayn_res[, c(2, 1)], delaunayn_res[, c(3, 1)], delaunayn_res[, c(3, 2)]) %>%
    unique %>% data.frame
  if (!is.null(dist_cutoff)) {
    dist_df <- as.matrix(dist(coord_df))
    delaunayn_el_dist <- dist_df[as.matrix(delaunayn_el)]
    dist_idx <- which(delaunayn_el_dist < dist_cutoff)
    delaunayn_el <- delaunayn_el[dist_idx, ]
  }
  if (return_name) {
    delaunayn_el$X1 <- rownames(coord_df)[delaunayn_el$X1]
    delaunayn_el$X2 <- rownames(coord_df)[delaunayn_el$X2]
  }
  return(delaunayn_el)
}

#' Title
#'
#' @param el_inp Edge list of Delaunay triangulation network
#' @param meta_data Meta data, must contain column of cell type
#' @param col_cell Column name of cell type, cell type names must be syntactically valid
#'
#' @return A list of 1.Network-based cell co-occurrence; 2.Cell co-occurrence log odds ratio
#'
#' @examples test <- edge_odds(el_inp, meta_data, col_cell='cell_names')
edge_odds <- function(el_inp, meta_data, col_cell='cell_names') {
  cell_el <- data.frame(X1=meta_data[el_inp[, 1], col_cell], X2=meta_data[el_inp[, 2], col_cell]) %>%
    dplyr::group_by(X1, X2) %>% dplyr::tally() %>% data.frame
  cell_adj_mat <- cell_el %>%
    tidyr::pivot_wider(id_cols=X1, names_from=X2, values_from=n, values_fill=list(n=0)) %>%
    data.frame(check.names = F) %>% magrittr::set_rownames(.$X1) %>% dplyr::select(-1)
  cell_adj_mat <- cell_adj_mat[, rownames(cell_adj_mat)]

  ## Calculate odds ratio ##
  cell_coo <- cell_adj_mat
  diag(cell_coo) <- 0
  ##
  cell_coo <- cell_coo + 1
  ##
  cell_row <- -sweep(cell_coo, 1, rowSums(cell_coo))
  cell_col <- -sweep(cell_coo, 2, colSums(cell_coo))
  cell_non <- -((cell_row + cell_col)-sum(cell_coo)+cell_coo)
  cell_odd <- (cell_coo*cell_non)/(cell_row*cell_col)
  cell_odd[cell_odd<0] <- 0
  cell_odd_log <- log(cell_odd)
  ## Test ##
  # cell_A='L2/3 IT'
  # cell_B='L4'
  # cell_test <- cell_el
  # cell_test <- cell_test[-which(cell_test$X1==cell_test$X2), ]
  # cell_test$X1[cell_test$X1!=cell_A] <- 'Non_A'
  # cell_test$X2[cell_test$X2!=cell_B] <- 'Non_B'
  # cell_test_count <- cell_test %>% group_by(X1, X2) %>% dplyr::summarise(n=sum(n)) %>%
  #   tidyr::pivot_wider(id_cols=X1, names_from=X2, values_from = n) %>% data.frame %>% set_rownames(.$X1) %>% dplyr::select(-1)
  # fisher.test(cell_test_count)
  ###
  return(list(cell_edgeCnt=cell_adj_mat, cell_logOR=cell_odd_log))
}

#' Title
#'
#' @param meta_df Meta data, must contain cell type column
#' @param coord_df Coordinates data, must be the same order of meta_df
#' @param col_cell Column name of cell type, cell type names must be syntactically valid
#' @param boot_n Number of bootstrapping iterations
#' @param prop Subsample proportion
#' @param dist_cutoff Remove edges with distance >= dist_cutoff
#' @param replace should sampling be with replacement?
#'
#' @return A list of 1. Bootstrap logOR distance; 2. MST consensus matrix
#'
#' @examples test <- DT_boot_mst(meta_df, coord_df, col_cell='cell_names', boot_n=100, prop=.8)
DT_boot_mst <- function(meta_df, coord_df, col_cell='cell_names', boot_n=100, prop=.8, replace=T, dist_cutoff=NULL) {
  cell_names <- unique(meta_df[, col_cell]) %>% sort
  n_smp <- nrow(meta_df)
  n_typ <- length(cell_names)
  dis_boot_array <- array(NA, dim=c(n_typ, n_typ, boot_n), dimnames=list(cell_names, cell_names, paste0('X_', c(1:boot_n))))
  mst_cons <- matrix(0, n_typ, n_typ) %>% data.frame %>% magrittr::set_rownames(cell_names) %>% magrittr::set_colnames(cell_names)
  for (i in 1:boot_n) {
    set.seed(i)
    cat(i, ' ')
    idx <- sort(sample(1:n_smp, round(n_smp*prop), replace=replace))
    data_boot <- meta_df[idx, ]
    coord_boot <- coord_df[idx, ]
    DT_el <- build_delaunayn(coord_boot, return_name = F, dist_cutoff=dist_cutoff)
    DT_OR <- edge_odds(DT_el, meta_data=data_boot, col_cell=col_cell)
    DT_OR_melt <- DT_OR$cell_logOR %>% tibble::rownames_to_column(var='Cell_A') %>%
      tidyr::pivot_longer(cols=-Cell_A, names_to='Cell_B', values_to='logOR') %>% data.frame
    dis_boot_temp <- dis_boot_array[, , i]
    ## Convert logOR to Distance ##
    dis_boot_temp[as.matrix(DT_OR_melt[, c(1, 2)])] <- max(DT_OR_melt$logOR) - DT_OR_melt$logOR + 1
    dis_boot_temp[is.na(dis_boot_temp)] <- max(dis_boot_temp)
    dis_boot_array[, , i] <- dis_boot_temp
    ## MST ##
    mst_boot <- igraph::graph.adjacency(dis_boot_temp, weighted=T, mode='upper') %>% igraph::mst(., algorithm='prim') %>%
      igraph::as_adjacency_matrix() %>% as.data.frame
    mst_boot <- mst_boot[cell_names, ]
    mst_boot <- mst_boot[, cell_names]
    mst_cons <- mst_cons + mst_boot/boot_n
  }
  output <- list(boot_array=dis_boot_array, mst_cons=mst_cons)
  return(output)
}

#' Title
#'
#' @param meta_df Meta data, must contain cell type column
#' @param coord_df Coordinates data, must be the same order of meta_df
#' @param col_cell Column name of cell type, cell type names must be syntactically valid
#' @param boot_n Number of bootstrapping iterations
#' @param k Number of NNs
#' @param prop Subsample proportion
#' @param replace should sampling be with replacement?
#'
#' @return A list of 1. Bootstrap logOR distance; 2. MST consensus matrix
#'
#' @examples test <- KD_boot_mst(meta_df, coord_df, col_cell='cell_names', boot_n=100, prop=.8, eps=1e-5)
KD_boot_mst <- function(meta_df, coord_df, col_cell='cell_names', boot_n=100, prop=.8, replace=T, k=10) {
  cell_names <- unique(meta_df[, col_cell]) %>% sort
  n_smp <- nrow(meta_df)
  n_typ <- length(cell_names)
  dis_boot_array <- array(NA, dim=c(n_typ, n_typ, boot_n), dimnames=list(cell_names, cell_names, paste0('X_', c(1:boot_n))))
  mst_cons <- matrix(0, n_typ, n_typ) %>% data.frame %>% magrittr::set_rownames(cell_names) %>% magrittr::set_colnames(cell_names)
  for (i in 1:boot_n) {
    set.seed(i)
    cat(i, ' ')
    idx <- sort(sample(1:n_smp, round(n_smp*prop), replace=replace))
    data_boot <- meta_df[idx, ]
    coord_boot <- coord_df[idx, ]
    inp_df <- data.frame(data_boot, coord_boot)
    kdist_i <- kdist(inp_df=inp_df, ref=unique(inp_df[, col_cell]), ref_type='each', que=unique(inp_df[, col_cell]), k=k, keep_nn=F)
    kdist_mat <- data.frame(kdist_i$kdist_df, cell_names=inp_df[, col_cell]) %>% dplyr::group_by(cell_names) %>% dplyr::summarise_all(mean) %>%
      data.frame %>% magrittr::set_rownames(.$cell_names) %>% dplyr::select(-1)
    colnames(kdist_mat) <- gsub('_kdist', '', colnames(kdist_mat))
    kdist_mat <- kdist_mat[cell_names, cell_names]
    dis_boot_temp <- log(as.matrix(kdist_mat))
    dis_boot_array[, , i] <- dis_boot_temp

    ## MST ##
    mst_boot <- igraph::graph.adjacency(dis_boot_temp, weighted=T, mode='directed') %>% igraph::mst(., algorithm='prim') %>%
      igraph::as_adjacency_matrix() %>% as.data.frame
    mst_boot <- mst_boot[cell_names, ]
    mst_boot <- mst_boot[, cell_names]
    mst_cons <- mst_cons + mst_boot/boot_n
  }
  output <- list(boot_array=dis_boot_array, mst_cons=mst_cons)
  return(output)
}

#' Title Calculated K-distance between query cells and reference cells based on their spatial coordinates
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

#' Title Run K-distance with CellTrek(Seurat) object and add a metadata column
#'
#' @param celltrek_inp SChart seurat input
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
#' @examples celltrek_test <- run_kdist(celltrek_inp=celltrek_inp, grp_col='cell_type', ref=c('A', 'B'), ref_type='each', que=unique(test_df$cell_names), k=10)
run_kdist <- function(celltrek_inp, grp_col='cell_type', ref=NULL, ref_type='all', que=NULL, k=10, new_name='kdist', keep_nn=F) {
  ## Check ##
  if (any(!c(grp_col, 'coord_x', 'coord_y') %in% colnames(celltrek_inp@meta.data))) {stop("celltrek_inp metadata must contain grp_col, 'coord_x', 'coord_y' columns")}

  inp_df <- celltrek_inp@meta.data %>% dplyr::select(cell_names=dplyr::one_of(grp_col), coord_x, coord_y)
  output <- kdist(inp_df=inp_df, ref=ref, ref_type=ref_type, que = que, k=k, new_name=new_name, keep_nn=keep_nn)
  celltrek_out <- Seurat::AddMetaData(celltrek_inp, metadata=output$kdist_df)
  return (celltrek_out)
}

#' Title SColoc module
#'
#' @param celltrek_inp CellTrek input
#' @param col_cell Column name of cell type, cell type names must be syntactically valid
#' @param h Bandwidths for x and y directions, for KL, more details in kde2d function in MASS package
#' @param n Number of grid points in each directions, for KL, more details in kde2d function in MASS package
#' @param ... See in boot_mst
#' @param use_method Use density-based Kullback Leibler divergence or Delaunay Triangulation network
#'
#' @return A list of 1.Bootstrap distance; 2.MST consensus matrix
#' @export
#'
#' @examples cell_scoloc <- scoloc(celltrek_inp, col_cell='cell_names', use_method='KL', h=140, n=25)
scoloc <- function(celltrek_inp, col_cell='cell_names', use_method=c('KL', 'DT', 'KD')[2], h=celltrek_inp@images[[1]]@scale.factors$spot_dis, n=25, boot_n=20, ...) {
  celltrek_temp <- celltrek_inp
  celltrek_temp$cell_names <- make.names(as.character(celltrek_temp@meta.data[, col_cell]))
  Idents(celltrek_temp) <- celltrek_temp$cell_names
  if (use_method=='KL') {
    cell_dummy_df <- as_dummy_df(celltrek_temp@meta.data %>% dplyr::select(id_new, cell_names, coord_x, coord_y), col_cell='cell_names')
    cell_mst_con <- KL_boot_mst(dummy_df=cell_dummy_df[, c(5:ncol(cell_dummy_df))], coord_df=cell_dummy_df[, c(3, 4)], h=h, n=n, boot_n=boot_n, ...)
  }
  if (use_method=='DT') {
    cell_mst_con <- DT_boot_mst(meta_df = celltrek_temp@meta.data, coord_df = celltrek_temp@meta.data[, c('coord_x', 'coord_y')], boot_n=boot_n, ...)
  }
  if (use_method=='KD') {
    cell_mst_con <- KD_boot_mst(meta_df = celltrek_temp@meta.data, coord_df = celltrek_temp@meta.data[, c('coord_x', 'coord_y')], boot_n=boot_n, ...)
  }
  return(cell_mst_con)
}
