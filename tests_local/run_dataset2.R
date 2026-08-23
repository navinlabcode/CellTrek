## Second paired dataset test on Seurat v5:
##   ST  = stxBrain anterior1 (10x Visium, updates to Assay5 + VisiumV2 image)
##   scRNA = allen_cortex reference (subclass labels)
## This exercises the VisiumV2 (v5 image class) code path, which the first
## dataset (brain_st_cortex, VisiumV1) did not.
suppressPackageStartupMessages({
  library(CellTrek); library(Seurat); library(dplyr); library(magrittr)
  library(ggplot2)
})
set.seed(1)
cat('Seurat', as.character(packageVersion('Seurat')), '\n')
dir.create('tests_local/figs', showWarnings = FALSE)

## ---- ST: stxBrain anterior1 ----
e <- new.env(); load('data2/anterior1.rda', envir = e)
st <- UpdateSeuratObject(get('anterior1', envir = e))
cat('ST anterior1:', ncol(st), 'spots | image class:', class(st@images[[1]]), '\n')
st <- NormalizeData(st, verbose = FALSE) %>% FindVariableFeatures(nfeatures = 2000, verbose = FALSE)
st <- RenameCells(st, new.names = make.names(Cells(st)))

## ---- scRNA: allen_cortex, subsample + preprocess ----
sc <- readRDS('data2/allen_cortex.rds')
sc$cell_type <- as.character(sc$subclass)
set.seed(42)
md <- sc@meta.data
keep <- unlist(lapply(split(rownames(md), md$cell_type), function(x) sample(x, min(length(x), 70))))
sc <- subset(sc, cells = keep)
sc <- NormalizeData(sc, verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>% RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE)
sc <- RenameCells(sc, new.names = make.names(Cells(sc)))
cat('scRNA allen:', ncol(sc), 'cells |', length(unique(sc$cell_type)), 'cell types\n')

## ---- Step 1: traint ----
cat('\n==== traint (VisiumV2 ST) ====\n')
tr <- CellTrek::traint(st_data = st, sc_data = sc, st_assay = 'Spatial',
                       sc_assay = 'RNA', cell_names = 'cell_type')
cat('traint:', ncol(tr), 'cells | coord_x set:',
    all(!is.na(tr$coord_x[tr$type == 'st'])), '\n')

## ---- Step 2: celltrek ----
cat('\n==== celltrek ====\n')
res <- CellTrek::celltrek(st_sc_int = tr, int_assay = 'traint', sc_data = sc,
                          sc_assay = 'RNA', reduction = 'pca', intp = TRUE,
                          intp_pnt = 5000, intp_lin = FALSE, nPCs = 30, ntree = 1000,
                          dist_thresh = 0.55, top_spot = 5, spot_n = 5,
                          repel_r = 20, repel_iter = 20, keep_model = FALSE)$celltrek
cat('celltrek:', ncol(res), 'charted cells | coords in meta:',
    all(c('coord_x','coord_y') %in% colnames(res@meta.data)), '\n')
saveRDS(res, 'tests_local/dataset2_celltrek_cache.rds')

## ---- Fig: spatial charting ----
df <- res@meta.data
df$cell_type <- factor(df$cell_type, levels = sort(unique(df$cell_type)))
p <- ggplot(df, aes(coord_y, max(coord_x)+min(coord_x)-coord_x, color = cell_type)) +
  geom_point(size = 0.6, alpha = 0.85) + coord_fixed() + theme_void() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = paste0('Dataset 2 (stxBrain anterior + allen_cortex): ',
                      nrow(df), ' cells charted'), color = 'subclass')
ggsave('tests_local/figs/4_dataset2_spatial.png', p, width = 9, height = 7, dpi = 130)
cat('spatial fig saved\n')

## ---- Step 3: scoloc (DT, no external dep) ----
cat('\n==== scoloc (DT) ====\n')
sg <- CellTrek::scoloc(res, col_cell = 'cell_type', use_method = 'DT', boot_n = 10)
cat('scoloc mst_cons:', paste(dim(sg$mst_cons), collapse = 'x'), '\n')

cat('\n==== DATASET 2 ALL STAGES COMPLETED ====\n')
