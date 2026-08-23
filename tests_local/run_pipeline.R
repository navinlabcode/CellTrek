## Local functional test of CellTrek on Seurat v5
## Runs the README workflow on the mouse brain data (subsampled for speed).
suppressPackageStartupMessages({
  library(CellTrek)
  library(Seurat)
  library(dplyr)
  library(magrittr)
})
set.seed(1)
cat('Seurat', as.character(packageVersion('Seurat')), '\n')

data_dir <- 'data'
brain_st <- readRDS(file.path(data_dir, 'brain_st_cortex.rds'))
brain_sc <- readRDS(file.path(data_dir, 'brain_sc.rds'))

## syntactically valid names
brain_st <- RenameCells(brain_st, new.names = make.names(Cells(brain_st)))
brain_sc <- RenameCells(brain_sc, new.names = make.names(Cells(brain_sc)))

## subsample SC for speed, keeping every cell type
set.seed(42)
md <- brain_sc@meta.data
keep <- unlist(lapply(split(rownames(md), md$cell_type), function(x) sample(x, min(length(x), 60))))
brain_sc <- subset(brain_sc, cells = keep)
cat('SC cells after subsample:', ncol(brain_sc), '| cell types:', length(unique(brain_sc$cell_type)), '\n')
cat('ST spots:', ncol(brain_st), '\n')

cache <- 'tests_local/brain_celltrek_cache.rds'
if (file.exists(cache)) {
  cat('\n==== loading cached celltrek object ====\n')
  brain_celltrek <- readRDS(cache)
} else {
## ---- Step 1: traint (co-embedding) ----
cat('\n==== traint ====\n')
brain_traint <- CellTrek::traint(st_data = brain_st, sc_data = brain_sc,
                                 sc_assay = 'RNA', cell_names = 'cell_type')
cat('traint object:', ncol(brain_traint), 'cells |',
    'reductions:', paste(names(brain_traint@reductions), collapse=','), '\n')
stopifnot('pca' %in% names(brain_traint@reductions))

## ---- Step 2: celltrek (charting) ----
cat('\n==== celltrek ====\n')
brain_celltrek <- CellTrek::celltrek(st_sc_int = brain_traint, int_assay = 'traint',
                                     sc_data = brain_sc, sc_assay = 'RNA',
                                     reduction = 'pca', intp = TRUE, intp_pnt = 5000,
                                     intp_lin = FALSE, nPCs = 30, ntree = 1000,
                                     dist_thresh = 0.55, top_spot = 5, spot_n = 5,
                                     repel_r = 20, repel_iter = 20, keep_model = TRUE)$celltrek
saveRDS(brain_celltrek, cache)
}
cat('celltrek object:', ncol(brain_celltrek), 'charted cells\n')
cat('coord range x:', paste(round(range(brain_celltrek$coord_x)), collapse='-'),
    '| y:', paste(round(range(brain_celltrek$coord_y)), collapse='-'), '\n')
stopifnot(all(c('coord_x','coord_y') %in% colnames(brain_celltrek@meta.data)))
stopifnot('celltrek' %in% names(brain_celltrek@reductions))

## ---- Step 3: SColoc (KL) ----
cat('\n==== scoloc (KL) ====\n')
glut_cell <- c('L2/3 IT','L4','L5 IT','L5 PT','NP','L6 IT','L6 CT','L6b')
glut_cell <- glut_cell[glut_cell %in% brain_celltrek$cell_type]
brain_glut <- subset(brain_celltrek, subset = cell_type %in% glut_cell)
cat('glut cells:', ncol(brain_glut), '\n')
brain_sgraph_KL <- CellTrek::scoloc(brain_glut, col_cell = 'cell_type',
                                    use_method = 'KL', eps = 1e-50, boot_n = 10)
cat('\nscoloc mst_cons dim:', paste(dim(brain_sgraph_KL$mst_cons), collapse='x'), '\n')

## ---- Step 4: SCoexp (cc) on L5 IT ----
cat('\n==== scoexp (cc) ====\n')
brain_l5 <- subset(brain_celltrek, subset = cell_type == 'L5 IT')
cat('L5 IT cells:', ncol(brain_l5), '\n')
if (ncol(brain_l5) >= 30) {
  brain_l5 <- FindVariableFeatures(brain_l5, verbose = FALSE)
  feats <- VariableFeatures(brain_l5)
  feats <- head(feats, 300)
  res_cc <- CellTrek::scoexp(celltrek_inp = brain_l5, assay = 'RNA', approach = 'cc',
                             gene_select = feats, sigm = 140, avg_cor_min = .3,
                             zero_cutoff = 3, min_gen = 10, max_gen = 200, maxK = 6, k = 6, reps = 10)
  cat('scoexp modules found:', length(res_cc$gs), '\n')
  cat('module sizes:', paste(sapply(res_cc$gs, length), collapse=','), '\n')
} else {
  cat('skip scoexp: too few L5 IT cells in subsample\n')
}

cat('\n==== ALL STAGES COMPLETED ====\n')
