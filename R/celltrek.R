#' ST and SC data integration using Seurat transfer anchors
#'
#' @param st_data Seurat ST data object
#' @param sc_data Seurat SC data object
#' @param st_assay ST assay
#' @param sc_assay SC assay
#' @param nfeatures number of features for integration
#' @param cell_names cell cluster/type column name in SC meta data
#' @param coord_xy coordinates column names in ST images slot
#' @param gene_kept selected genes to be kept during integration
#' @param norm normalization method: LogNormalize/SCTransform
#'
#' @return Seurat object of SC-ST co-embeddings
#' @export
#'
#' @examples
#' st_sc_traint <- traint(st_data=brain_st, sc_data=brain_sc, st_assay='Spatial', sc_assay='scint', nfeatures=2000, cell_names='cell_names', coord_xy=c('imagerow', 'imagecol'), gene_kept=NULL)
traint <- function (st_data, sc_data, st_assay='Spatial', sc_assay='scint', norm='LogNormalize', nfeatures=2000,
                    cell_names='cell_names', coord_xy=c('imagerow', 'imagecol'), gene_kept=NULL, ...) {

  st_data$id <- names(st_data$orig.ident)
  sc_data$id <- names(sc_data$orig.ident)
  sc_data$cell_names <- make.names(sc_data@meta.data[, cell_names])
  st_data$type <- 'st'
  sc_data$type <- 'sc'
  st_data$coord_x <- st_data@images[[1]]@coordinates[, coord_xy[1]]
  st_data$coord_y <- st_data@images[[1]]@coordinates[, coord_xy[2]]
  DefaultAssay(st_data) <- st_assay
  DefaultAssay(sc_data) <- sc_assay

  cat('Finding transfer anchors... \n')
  st_idx <- st_data$id
  sc_idx <- sc_data$id

  ## Integration features ##
  sc_st_list <- list(st_data=st_data, sc_data=sc_data)
  sc_st_features <- SelectIntegrationFeatures(sc_st_list, nfeatures=nfeatures)
  if (!is.null(gene_kept)) {
    sc_st_features <- union(sc_st_features, gene_kept)
  }

  sc_st_features <- sc_st_features[(sc_st_features %in% rownames(st_data[[st_assay]]@data)) &
                                     (sc_st_features %in% rownames(sc_data[[sc_assay]]@data))]
  cat('Using', length(sc_st_features), 'features for integration... \n')
  ###

  sc_st_anchors <- Seurat::FindTransferAnchors(reference = sc_data, query = st_data,
                                               reference.assay = sc_assay, query.assay = st_assay,
                                               normalization.method = norm, features = sc_st_features, reduction = 'cca', ...)

  cat('Data transfering... \n')
  st_data_trans <- TransferData(anchorset = sc_st_anchors,
                                refdata = GetAssayData(sc_data, assay = sc_assay, slot='data')[sc_st_features, ],
                                weight.reduction = 'cca')
  st_data@assays$transfer <- st_data_trans

  cat('Creating new Seurat object... \n')
  sc_st_meta <- dplyr::bind_rows(st_data@meta.data, sc_data@meta.data)
  counts_temp <- cbind(data.frame(st_data[['transfer']]@data), data.frame(sc_data[[sc_assay]]@data[sc_st_features, ] %>% data.frame))
  rownames(sc_st_meta) <- make.names(sc_st_meta$id)
  colnames(counts_temp) <- make.names(sc_st_meta$id)
  sc_st_int <- CreateSeuratObject(counts = counts_temp, assay = 'traint', meta.data = sc_st_meta)
  sc_st_int[['traint']]@data <- sc_st_int[['traint']]@counts
  sc_st_int[['traint']]@counts <- matrix(NA, nrow = 0, ncol = 0)

  cat('Scaling -> PCA -> UMAP... \n')
  sc_st_int <- ScaleData(sc_st_int, features = sc_st_features) %>%
    RunPCA(features = sc_st_features)
  sc_st_int <- RunUMAP(sc_st_int, dims = 1:30)
  sc_st_int@images <- st_data@images
  sc_st_int@images[[1]]@coordinates <- data.frame(imagerow=sc_st_int@meta.data$coord_x,
                                                  imagecol=sc_st_int@meta.data$coord_y) %>%
    set_rownames(rownames(sc_st_int@meta.data))

  return (sc_st_int)
}

#' Title Mannually repelling celltrek cells
#'
#' @param celltrek_inp CellTrek(Seurat) object
#' @param repel_r Repelling radius
#' @param repel_iter Repelling iterations
#'
#' @return CellTrek(Seurat) object
#' @export
#'
#' @examples
celltrek_repel <- function(celltrek_inp, repel_r=5, repel_iter=10) {
  celltrek_dr_raw <- Embeddings(celltrek_inp, 'celltrek_raw')
  celltrek_dr <- Embeddings(celltrek_inp, 'celltrek')
  repel_input <- data.frame(celltrek_dr_raw, repel_r=repel_r)

  ## Add noise ##
  theta <- runif(nrow(celltrek_dr), 0, 2*pi)
  alpha <- sqrt(runif(nrow(celltrek_dr), 0, 1))
  repel_input[, 1] <- repel_input[, 1] + sin(theta)*alpha*repel_r
  repel_input[, 2] <- repel_input[, 2] + cos(theta)*alpha*repel_r

  ## Repelling ##
  cat('Repelling points...\n')
  celltrek_repel <- packcircles::circleRepelLayout(repel_input, sizetype='radius', maxiter=repel_iter)
  celltrek_dr[, 1] <- celltrek_repel$layout$x
  celltrek_dr[, 2] <- celltrek_repel$layout$y
  celltrek_out <- celltrek_inp
  celltrek_out@reductions$celltrek@cell.embeddings <- celltrek_dr
  return(celltrek_out)
}

#' Title Calculate the RF-distance between sc and st
#'
#' @param st_sc_int Seurat traint object
#' @param int_assay Name of integration assay
#' @param reduction Dimension reduction method used, usually pca
#' @param intp If TRUE, do interpolation
#' @param intp_pnt Interpolation point number
#' @param intp_lin If TRUE, use linear interpolation
#' @param nPCs Number of PCs used for CellTrek
#' @param ntree Number of trees in random forest
#' @param keep_model If TRUE, return the trained random forest model
#'
#' @return A list of 1. celltrek_distance matrix; 2. trained random forest model (optional)
#'
#' @examples dist_test <- celltrek_dist(st_sc_int=st_sc_int, int_assay='traint', reduction='pca', intp = T, intp_pnt=10000, intp_lin=F, nPCs=30, ntree=1000, keep_model=T)
celltrek_dist <- function (st_sc_int, int_assay='traint', reduction='pca', intp = T, intp_pnt=10000, intp_lin=F, nPCs=30, ntree=1000, keep_model=T) {
  DefaultAssay(st_sc_int) <- int_assay
  kNN_dist <- dbscan::kNN(na.omit(st_sc_int@meta.data[, c('coord_x', 'coord_y')]), k=6)$dist
  spot_dis <- median(kNN_dist) %>% round
  cat('Distance between spots is:', spot_dis, '\n')

  st_sc_int$id <- names(st_sc_int$orig.ident)
  st_idx <- st_sc_int$id[st_sc_int$type=='st']
  sc_idx <- st_sc_int$id[st_sc_int$type=='sc']
  meta_df <- data.frame(st_sc_int@meta.data)

  st_sc_int_pca <- st_sc_int@reductions[[reduction]]@cell.embeddings[, 1:nPCs] %>% data.frame %>%
    mutate(id=st_sc_int$id, type=st_sc_int$type, class=st_sc_int$cell_names,
           coord_x=st_sc_int$coord_x, coord_y=st_sc_int$coord_y)
  st_pca <- st_sc_int_pca %>% dplyr::filter(type=='st') %>% dplyr::select(-c(id:class))

  ## Interpolation ##
  ## Uniform sampling ##
  if (intp) {
    cat ('Interpolating...\n')
    spot_ratio <- intp_pnt/nrow(st_pca)
    st_intp_df <- apply(st_pca[, c('coord_x', 'coord_y')], 1, function(row_x) {
      runif_test <- runif(1)
      if (runif_test < spot_ratio%%1) {
        theta <- runif(ceiling(spot_ratio), 0, 2*pi)
        alpha <- sqrt(runif(ceiling(spot_ratio), 0, 1))
        coord_x <- row_x[1] + (spot_dis/2)*sin(theta)*alpha
        coord_y <- row_x[2] + (spot_dis/2)*cos(theta)*alpha
      } else {
        theta <- runif(floor(spot_ratio), 0, 2*pi)
        alpha <- sqrt(runif(floor(spot_ratio), 0, 1))
        coord_x <- row_x[1] + (spot_dis/2)*sin(theta)*alpha
        coord_y <- row_x[2] + (spot_dis/2)*cos(theta)*alpha
      }
      data.frame(coord_x, coord_y)
    }) %>% Reduce(rbind, .)

    st_intp_df <- apply(st_pca[, 1:nPCs], 2, function(col_x) {
      akima::interpp(x=st_pca$coord_x, y=st_pca$coord_y, z=col_x,
                     linear=intp_lin, xo=st_intp_df$coord_x, yo=st_intp_df$coord_y) %>%
        magrittr::extract2('z')
    }) %>% data.frame(., id='X', type='st_intp', st_intp_df) %>% na.omit
    st_intp_df$id <- make.names(st_intp_df$id, unique = T)
    st_sc_int_pca <- bind_rows(st_sc_int_pca, st_intp_df)
  }

  cat('Random Forest training... \n')
  ## Training on ST ##
  data_train <- st_sc_int_pca %>% dplyr::filter(type=='st') %>% dplyr::select(-c(id:class))
  rf_train <- randomForestSRC::rfsrc(Multivar(coord_x, coord_y) ~ ., data_train, block.size=5, ntree=ntree)

  cat('Random Forest prediction...  \n')
  ## Testing on all ##
  data_test <- st_sc_int_pca
  rf_pred <- randomForestSRC::predict.rfsrc(rf_train, newdata=data_test[, c(1:nPCs)], distance='all')

  cat('Making distance matrix... \n')
  rf_pred_dist <- rf_pred$distance[data_test$type=='sc', data_test$type!='sc'] %>%
    set_rownames(data_test$id[data_test$type=='sc']) %>% set_colnames(data_test$id[data_test$type!='sc'])

  output <- list()
  output$spot_d <- spot_dis
  output$celltrek_dist <- rf_pred_dist
  output$coord_df <- st_sc_int_pca[, c('id', 'type', 'coord_x', 'coord_y')] %>%
    dplyr::filter(type!='sc') %>% magrittr::set_rownames(.$id) %>% dplyr::select(-id)
  if (keep_model) {
    output$model <- rf_train
  }
  return (output)
}

#' Title
#'
#' @param dist_mat Distance matrix of sc-st (sc in rows and st in columns)
#' @param coord_df Coordinates data frame of st (must contain coord_x, coord_y columns, barcode rownames)
#' @param dist_cut Distance cutoff
#' @param top_spot Maximum number of spots that one cell can be charted
#' @param spot_n Maximum number of cells that one spot can contain
#' @param repel_r Repelling radius
#' @param repel_iter Repelling iterations
#'
#' @return SC coordinates
#'
#' @examples chart_test <- function (dist_mat, coord_df, dist_cut=500, top_spot=10, spot_n=10, r=50)
celltrek_chart <- function (dist_mat, coord_df, dist_cut=500, top_spot=10, spot_n=10, repel_r=5, repel_iter=10) {
  cat('Making graph... \n')
  dist_mat[dist_mat>dist_cut] <- NA
  dist_mat_dt <- data.table::data.table(Var1=rownames(dist_mat), dist_mat)
  dist_edge_list <- data.table::melt(dist_mat_dt, id=1, na.rm=T)
  colnames(dist_edge_list) <- c('Var1', 'Var2', 'value')
  dist_edge_list$val_rsc <- scales::rescale(dist_edge_list$value, to=c(0, repel_r))
  dist_edge_list$Var1 %<>% as.character
  dist_edge_list$Var2 %<>% as.character
  dist_edge_list$Var1_type <- 'sc'
  dist_edge_list$Var2_type <- 'non-sc'

  cat('Pruning graph...\n')
  dist_edge_list_sub <- dplyr::inner_join(dist_edge_list %>% group_by(Var1) %>% top_n(n=top_spot, wt=-value),
                                          dist_edge_list %>% group_by(Var2) %>% top_n(n=spot_n, wt=-value)) %>% data.frame

  cat('Spatial Charting SC data...\n')
  sc_coord <- sc_coord_raw <- data.frame(id_raw=dist_edge_list_sub$Var1, id_new=make.names(dist_edge_list_sub$Var1, unique = T))
  sc_coord$coord_x <- sc_coord_raw$coord_x <- coord_df$coord_x[match(dist_edge_list_sub$Var2, rownames(coord_df))]
  sc_coord$coord_y <- sc_coord_raw$coord_y <- coord_df$coord_y[match(dist_edge_list_sub$Var2, rownames(coord_df))]
  ## Add noise ##
  theta <- runif(nrow(dist_edge_list_sub), 0, 2*pi)
  alpha <- sqrt(runif(nrow(dist_edge_list_sub), 0, 1))
  sc_coord$coord_x <- sc_coord$coord_x + dist_edge_list_sub$val_rsc*sin(theta)*alpha
  sc_coord$coord_y <- sc_coord$coord_y + dist_edge_list_sub$val_rsc*cos(theta)*alpha
  ## Point repelling ##
  cat('Repelling points...\n')
  sc_repel_input <- data.frame(sc_coord[, c('coord_x', 'coord_y')], repel_r=repel_r)
  sc_repel <- packcircles::circleRepelLayout(sc_repel_input, sizetype='radius', maxiter=repel_iter)
  sc_coord$coord_x <- sc_repel$layout$x
  sc_coord$coord_y <- sc_repel$layout$y
  return(list(sc_coord_raw, sc_coord))
}

#' Title
#'
#' @param dist_mat Distance matrix of sc-st (sc in rows and st in columns)
#' @param coord_df Coordinates data frame of st (must contain two columns as coord_x, coord_y and rownames as barcodes)
#' @param dist_cut Distance cutoff
#' @param top_spot Maximum number of spots that one cell can be charted
#' @param spot_n Maximum number of cells that one spot can contain
#' @param sc_data SC data
#' @param sc_assay SC assay
#' @param repel_r Repelling radius
#' @param repel_iter Repelling iterations
#' @param st_data ST data, optional
#'
#' @return A list of 1.Seurat object
#' @export
#'
#' @examples celltrek_res <- celltrek_from_dist(dist_mat, coord_df, dist_cut, top_spot=10, spot_n=10, r=NULL, sc_data, sc_assay='RNA')
celltrek_from_dist <- function (dist_mat, coord_df, dist_cut, top_spot=10, spot_n=10, repel_r=5, repel_iter=10, sc_data, sc_assay='RNA', st_data=NULL) {
  colnames(coord_df) <- c('coord_x', 'coord_y')
  spot_dis <- median(unlist(dbscan::kNN(coord_df[, c('coord_x', 'coord_y')], k=4)$dist))
  if (is.null(repel_r)) {repel_r=spot_dis/4}
  sc_coord_list <- celltrek_chart(dist_mat=dist_mat, coord_df=coord_df, dist_cut=dist_cut, top_spot=top_spot, spot_n=spot_n, repel_r=repel_r, repel_iter=repel_iter)
  sc_coord_raw <- sc_coord_list[[1]]
  sc_coord <- sc_coord_list[[2]]
  sc_out <- CreateSeuratObject(counts=sc_data[[sc_assay]]@data[, sc_coord$id_raw] %>% set_colnames(sc_coord$id_new),
                               project='celltrek', assay=sc_assay,
                               meta.data=sc_data@meta.data[sc_coord$id_raw, ] %>%
                                 dplyr::rename(id_raw=id) %>%
                                 mutate(id_new=sc_coord$id_new) %>%
                                 set_rownames(sc_coord$id_new))
  sc_out@meta.data <- dplyr::left_join(sc_out@meta.data, sc_coord) %>% data.frame %>% set_rownames(sc_out$id_new)
  sc_out[[sc_assay]]@data <- sc_out[[sc_assay]]@counts
  sc_out[[sc_assay]]@counts <- matrix(nrow = 0, ncol = 0)
  sc_coord_raw_df <- CreateDimReducObject(embeddings=sc_coord_raw %>%
                                            dplyr::mutate(coord1=coord_y, coord2=max(coord_x)+min(coord_x)-coord_x) %>%
                                            dplyr::select(c(coord1, coord2)) %>% set_rownames(sc_coord_raw$id_new) %>% as.matrix,
                                          assay=sc_assay, key='celltrek_raw')
  sc_coord_dr <- CreateDimReducObject(embeddings=sc_coord %>%
                                        dplyr::mutate(coord1=coord_y, coord2=max(coord_x)+min(coord_x)-coord_x) %>%
                                        dplyr::select(c(coord1, coord2)) %>% set_rownames(sc_coord$id_new) %>% as.matrix,
                                      assay=sc_assay, key='celltrek')
  sc_pca_dr <- CreateDimReducObject(embeddings=sc_data@reductions$pca@cell.embeddings[sc_coord$id_raw, ] %>%
                                      set_rownames(sc_coord$id_new) %>% as.matrix, assay=sc_assay, key='pca')
  sc_umap_dr <- CreateDimReducObject(embeddings=sc_data@reductions$umap@cell.embeddings[sc_coord$id_raw, ] %>%
                                       set_rownames(sc_coord$id_new) %>% as.matrix, assay=sc_assay, key='umap')
  sc_out@reductions$celltrek <- sc_coord_dr
  sc_out@reductions$celltrek_raw <- sc_coord_raw_df
  sc_out@reductions$pca <- sc_pca_dr
  sc_out@reductions$umap <- sc_umap_dr
  if (!is.null(st_data)) {
    sc_out@images <- st_data@images
    sc_out@images[[1]]@assay <- DefaultAssay(sc_out)
    sc_out@images[[1]]@coordinates <- data.frame(imagerow=sc_coord$coord_x, imagecol=sc_coord$coord_y) %>% set_rownames(sc_coord$id_new)
    sc_out@images[[1]]@scale.factors$spot_dis <- spot_dis
  }
  output <- list(celltrek=sc_out)
  return(output)
}

#' Title
#'
#' @param st_sc_int Seurat traint object
#' @param int_assay Integration assay ('traint')
#' @param sc_data SC data, optional
#' @param sc_assay SC assay
#' @param reduction Dimension reduction method, usually 'pca'
#' @param intp If True, do interpolation
#' @param intp_pnt Number of interpolation points
#' @param intp_lin If Ture, do linear interpolation
#' @param nPCs Number of PCs
#' @param ntree Number of Trees
#' @param dist_thresh Distance threshold
#' @param top_spot Maximum number of spots that one cell can be charted
#' @param spot_n Maximum number of cells that one spot can contain
#' @param keep_model If TRUE, return the trained random forest model
#' @param repel_r Repelling radius
#' @param repel_iter Repelling iterations
#' @param ...
#'
#' @return Seurat object
#' @export
#'
#' @examples celltrek_res <- celltrek(st_sc_int, int_assay='traint', sc_data=NULL, sc_assay='RNA', reduction='pca', intp=T, intp_pnt=10000, intp_lin=F, nPCs=30, ntree=1000, dist_thresh=.4, top_spot=10, spot_n=10, r=NULL, keep_model=F, ...)
celltrek <- function (st_sc_int, int_assay='traint', sc_data=NULL, sc_assay='RNA', reduction='pca',
                    intp=T, intp_pnt=10000, intp_lin=F, nPCs=30,
                    ntree=1000, dist_thresh=.4, top_spot=10, spot_n=10, repel_r=5, repel_iter=10, keep_model=F, ...) {
  dist_res <- celltrek_dist(st_sc_int=st_sc_int, int_assay=int_assay, reduction=reduction, intp=intp, intp_pnt=intp_pnt, intp_lin=intp_lin, nPCs=nPCs, ntree=ntree, keep_model=T)
  spot_dis_intp <- median(unlist(dbscan::kNN(dist_res$coord_df[, c('coord_x', 'coord_y')], k=4)$dist))
  if (is.null(repel_r)) {repel_r=spot_dis_intp/4}
  sc_coord_list <- celltrek_chart(dist_mat=dist_res$celltrek_dist, coord_df=dist_res$coord_df, dist_cut=ntree*dist_thresh, top_spot=top_spot, spot_n=spot_n, repel_r=repel_r, repel_iter=repel_iter)
  sc_coord_raw <- sc_coord_list[[1]]
  sc_coord <- sc_coord_list[[2]]
  cat('Creating Seurat Object... \n')
  if (!is.null(sc_data)) {
    cat('sc data...')
    sc_data$id <- Seurat::Cells(sc_data)
    sc_out <- CreateSeuratObject(counts=sc_data[[sc_assay]]@data[, sc_coord$id_raw] %>% set_colnames(sc_coord$id_new),
                                 project='celltrek', assay=sc_assay,
                                 meta.data=sc_data@meta.data[sc_coord$id_raw, ] %>%
                                   dplyr::rename(id_raw=id) %>%
                                   mutate(id_new=sc_coord$id_new) %>%
                                   set_rownames(sc_coord$id_new))
    sc_out@meta.data <- dplyr::left_join(sc_out@meta.data, sc_coord) %>% data.frame %>% set_rownames(sc_out$id_new)

    sc_out[[sc_assay]]@data <- sc_out[[sc_assay]]@counts
    sc_out[[sc_assay]]@counts <- matrix(nrow = 0, ncol = 0)
    sc_coord_raw_df <- CreateDimReducObject(embeddings=sc_coord_raw %>%
                                              dplyr::mutate(coord1=coord_y, coord2=max(coord_x)+min(coord_x)-coord_x) %>%
                                              dplyr::select(c(coord1, coord2)) %>% set_rownames(sc_coord_raw$id_new) %>% as.matrix,
                                            assay=sc_assay, key='celltrek_raw')
    sc_coord_dr <- CreateDimReducObject(embeddings=sc_coord %>%
                                          dplyr::mutate(coord1=coord_y, coord2=max(coord_x)+min(coord_x)-coord_x) %>%
                                          dplyr::select(c(coord1, coord2)) %>% set_rownames(sc_coord$id_new) %>% as.matrix,
                                        assay=sc_assay, key='celltrek')
    sc_pca_dr <- CreateDimReducObject(embeddings=sc_data@reductions$pca@cell.embeddings[sc_coord$id_raw, ] %>%
                                        set_rownames(sc_coord$id_new) %>% as.matrix, assay=sc_assay, key='pca')
    sc_umap_dr <- CreateDimReducObject(embeddings=sc_data@reductions$umap@cell.embeddings[sc_coord$id_raw, ] %>%
                                         set_rownames(sc_coord$id_new) %>% as.matrix, assay=sc_assay, key='umap')
    sc_out@reductions$celltrek <- sc_coord_dr
    sc_out@reductions$celltrek_raw <- sc_coord_raw_df
    sc_out@reductions$pca <- sc_pca_dr
    sc_out@reductions$umap <- sc_umap_dr
  } else {
    cat('no sc data...')
    sc_out <- CreateSeuratObject(counts=st_sc_int[[int_assay]]@data[, sc_coord$id_raw] %>%
                                   set_colnames(sc_coord$id_new),
                                 project='celltrek', assay=int_assay,
                                 meta.data=st_sc_int@meta.data[sc_coord$id_raw, ] %>%
                                   dplyr::rename(id_raw=id) %>%
                                   mutate(id_new=sc_coord$id_new) %>%
                                   set_rownames(sc_coord$id_new))
    sc_out$coord_x <- sc_coord$coord_x[match(sc_coord$id_new, sc_out$id_new)]
    sc_out$coord_y <- sc_coord$coord_y[match(sc_coord$id_new, sc_out$id_new)]

    sc_out[[int_assay]]@counts <- matrix(nrow = 0, ncol = 0)
    sc_out[[int_assay]]@scale.data <- st_sc_int[[int_assay]]@scale.data[, sc_coord$id_raw] %>% set_colnames(sc_coord$id_new)
    sc_coord_raw_df <- CreateDimReducObject(embeddings=sc_coord_raw %>%
                                              dplyr::mutate(coord1=coord_y, coord2=max(coord_x)+min(coord_x)-coord_x) %>%
                                              dplyr::select(c(coord1, coord2)) %>% set_rownames(sc_coord_raw$id_new) %>% as.matrix,
                                            assay=sc_assay, key='celltrek_raw')
    sc_coord_dr <- CreateDimReducObject(embeddings = sc_coord %>%
                                          dplyr::mutate(coord1=coord_y, coord2=max(coord_x)+min(coord_x)-coord_x) %>%
                                          dplyr::select(c(coord1, coord2)) %>%
                                          set_rownames(sc_coord$id_new) %>%
                                          as.matrix,
                                        assay=int_assay, key='celltrek')
    sc_pca_dr <- CreateDimReducObject(embeddings = st_sc_int@reductions$pca@cell.embeddings[sc_coord$id_raw, ] %>%
                                        set_rownames(sc_coord$id_new) %>% as.matrix, assay=int_assay, key='pca')
    sc_umap_dr <- CreateDimReducObject(embeddings = st_sc_int@reductions$umap@cell.embeddings[sc_coord$id_raw, ] %>%
                                         set_rownames(sc_coord$id_new) %>% as.matrix, assay=int_assay, key='umap')

    sc_out@reductions$celltrek <- sc_coord_dr
    sc_out@reductions$celltrek_raw <- sc_coord_raw_df
    sc_out@reductions$pca <- sc_pca_dr
    sc_out@reductions$umap <- sc_umap_dr
  }
  sc_out@images <- st_sc_int@images
  sc_out@images[[1]]@assay <- DefaultAssay(sc_out)
  sc_out@images[[1]]@coordinates <- data.frame(imagerow=sc_coord$coord_x, imagecol=sc_coord$coord_y) %>% set_rownames(sc_coord$id_new)
  sc_out@images[[1]]@scale.factors$spot_dis <- dist_res$spot_d
  sc_out@images[[1]]@scale.factors$spot_dis_intp <- spot_dis_intp

  output <- list(celltrek=sc_out)
  if (keep_model) {
    output[[length(output)+1]] <- dist_res$model
    names(output)[length(output)] <- 'model'
  }
  return(output)
}
