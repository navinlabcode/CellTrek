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
