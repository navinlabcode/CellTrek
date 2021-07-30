#' Title
#'
#' @param data_inp Cell-coordinates matrix
#' @param col_cell Column name of cell type, cell type names must be syntactically valid
#'
#' @return dummy data frame with cells on columns
#' @export
#'
#' @examples cell_dummy_df <- as_dummy_df(schart_df, col_cell='Cell')
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
#' @export
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
#' @export
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
        extract2('z') %>% melt %>% data.frame %>% set_colnames(c('grid_1', 'grid_2', 'Val'))
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
#' @export
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
#' @export
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

#' Title
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

## Check Delaunayn ##
# coord_df <- data.frame(x=c(rnorm(100), rnorm(100, 10)), y=c(rnorm(100), rnorm(100, 10)))
# coord_df <- data.frame(x=HBCA40_schart$coord_x, y=HBCA40_schart$coord_y)
# test <- SChartPack::build_delaunayn(coord_df, return_name=F, dist_cutoff=100)
# coord_df_segment <- data.frame(x=coord_df$x[test$X1], y=coord_df$y[test$X1], xend=coord_df$x[test$X2], yend=coord_df$y[test$X2])
# ggplot() + geom_point(data=coord_df, aes(x, y)) + geom_segment(data=coord_df_segment, aes(x=x, y=y, xend=xend, yend=yend))

#' Title
#'
#' @param el_inp Edge list of Delaunay triangulation network
#' @param meta_data Meta data, must contain column of cell type
#' @param col_cell Column name of cell type, cell type names must be syntactically valid
#'
#' @return A list of 1.Network-based cell co-occurrence; 2.Cell co-occurrence log odds ratio
#' @export
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
#' @export
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
#' @export
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



#' Title
#'
#' @param schart_inp SChart input
#' @param col_cell Column name of cell type, cell type names must be syntactically valid
#' @param h Bandwidths for x and y directions, for KL, more details in kde2d function in MASS package
#' @param n Number of grid points in each directions, for KL, more details in kde2d function in MASS package
#' @param ... See in boot_mst
#' @param use_method Use density-based Kullback Leibler divergence or Delaunay Triangulation network
#' @param cell_min Minimum number of cells in each cell type to be calculated
#'
#' @return A list of 1.Bootstrap distance; 2.MST consensus matrix
#' @export
#'
#' @examples cell_scoloc <- scoloc(schart_inp, col_cell='cell_names', use_method='KL', h=140, n=25)
scoloc <- function(schart_inp, col_cell='cell_names', cell_min=15, use_method=c('KL', 'DT', 'KD')[2], h=schart_inp@images[[1]]@scale.factors$spot_dis, n=25, boot_n=20, ...) {
  schart_temp <- schart_inp
  schart_temp$cell_names <- make.names(as.character(schart_temp@meta.data[, col_cell]))
  Idents(schart_temp) <- schart_temp$cell_names
  cell_typ_kept <- names(table(schart_temp@meta.data$cell_names))[table(schart_temp@meta.data$cell_names)>cell_min]
  schart_temp <- suppressWarnings(Seurat::SubsetData(schart_temp, ident.use=cell_typ_kept))
  if (use_method=='KL') {
    cell_dummy_df <- as_dummy_df(schart_temp@meta.data %>% dplyr::select(id_new, cell_names, coord_x, coord_y), col_cell='cell_names')
    cell_mst_con <- KL_boot_mst(dummy_df=cell_dummy_df[, c(5:ncol(cell_dummy_df))], coord_df=cell_dummy_df[, c(3, 4)], h=h, n=n, boot_n=boot_n, ...)
  }
  if (use_method=='DT') {
    cell_mst_con <- DT_boot_mst(meta_df = schart_temp@meta.data, coord_df = schart_temp@meta.data[, c('coord_x', 'coord_y')], boot_n=boot_n, ...)
  }
  if (use_method=='KD') {
    cell_mst_con <- KD_boot_mst(meta_df = schart_temp@meta.data, coord_df = schart_temp@meta.data[, c('coord_x', 'coord_y')], boot_n=boot_n, ...)
  }
  return(cell_mst_con)
}
