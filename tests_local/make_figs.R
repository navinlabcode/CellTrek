suppressPackageStartupMessages({
  library(CellTrek); library(Seurat); library(dplyr)
  library(ggplot2); library(pheatmap); library(viridis)
})
set.seed(1)
dir.create('tests_local/figs', showWarnings = FALSE)
ct <- readRDS('tests_local/brain_celltrek_cache.rds')

## ---- Fig 1: spatial charting (charted single cells on tissue coords) ----
df <- ct@meta.data
df$cell_type <- factor(df$cell_type, levels = sort(unique(df$cell_type)))
p1 <- ggplot(df, aes(coord_y, max(coord_x)+min(coord_x)-coord_x, color = cell_type)) +
  geom_point(size = 0.7, alpha = 0.85) +
  coord_fixed() + theme_void() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = paste0('CellTrek charting: ', nrow(df), ' single cells mapped to tissue'),
       color = 'cell type')
ggsave('tests_local/figs/1_celltrek_spatial.png', p1, width = 9, height = 7, dpi = 130)
cat('fig1 done\n')

## ---- Fig 2: SColoc (KL) consensus co-localization matrix ----
glut <- c('L2/3 IT','L4','L5 IT','L5 PT','NP','L6 IT','L6 CT','L6b')
g <- subset(ct, subset = cell_type %in% glut)
sg <- CellTrek::scoloc(g, col_cell = 'cell_type', use_method = 'KL', eps = 1e-50, boot_n = 20)
mst <- as.matrix(sg$mst_cons)
labs <- glut; names(labs) <- make.names(glut)
rownames(mst) <- colnames(mst) <- labs[rownames(mst)]
mst_sym <- mst + t(mst)
png('tests_local/figs/2_scoloc_KL_mst.png', width = 900, height = 800, res = 130)
pheatmap(mst_sym, cluster_rows = TRUE, cluster_cols = TRUE,
         color = viridis(100), display_numbers = TRUE, number_format = '%.2f',
         main = 'SColoc (KL): MST consensus co-localization\n(glutamatergic neurons, 20 bootstraps)')
dev.off()
cat('fig2 done\n')

## ---- Fig 3: SCoexp (cc) spatial co-expression heatmap on L5 IT ----
l5 <- subset(ct, subset = cell_type == 'L5 IT')
l5 <- FindVariableFeatures(l5, verbose = FALSE)
vst_df <- HVFInfo(l5) %>% data.frame %>% mutate(id = rownames(.))
nz <- apply(as.matrix(GetAssayData(l5, assay='RNA', layer='data')), 1, function(x) mean(x!=0)*100)
hz <- names(nz)[nz < 20]
mt <- grep('^Mt-', rownames(l5), value = TRUE)
rp <- grep('^Rpl|^Rps', rownames(l5), value = TRUE)
vst_df <- vst_df %>% filter(!(id %in% c(mt, rp, hz))) %>% arrange(-variance.standardized)
feats <- head(na.omit(vst_df$id), 300)
res <- CellTrek::scoexp(celltrek_inp = l5, assay = 'RNA', approach = 'cc',
                        gene_select = feats, sigm = 140, avg_cor_min = .3,
                        zero_cutoff = 3, min_gen = 15, max_gen = 200, maxK = 6, k = 6, reps = 20)
k_df <- do.call(rbind, lapply(seq_along(res$gs), function(i)
  data.frame(gene = res$gs[[i]], Module = paste0('K', i))))
rownames(k_df) <- k_df$gene; k_df$gene <- NULL
gord <- rownames(k_df)
png('tests_local/figs/3_scoexp_heatmap.png', width = 850, height = 800, res = 130)
pheatmap(res$wcor[gord, gord], clustering_method = 'ward.D2',
         annotation_row = k_df, show_rownames = FALSE, show_colnames = FALSE,
         treeheight_row = 12, treeheight_col = 12, fontsize = 8,
         color = viridis(10), main = 'SCoexp: L5 IT spatial co-expression modules')
dev.off()
cat('fig3 done | modules:', length(res$gs), '| sizes:', paste(sapply(res$gs, length), collapse=','), '\n')
cat('ALL FIGS DONE\n')
