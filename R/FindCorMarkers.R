#' Title
#'
#' @param srt_inp Seurat input
#' @param assay 
#' @param features Features that genes correlated with
#' @param method Correlation method
#'
#' @return
#' @export
#'
#' @examples CorMarkers <- FindCorMarkers(srt_inp=test_seurat, assay='RNA', features=c('k_dist', 'sig_score'), method='spearman')
FindCorMarkers <- function(srt_inp, assay='RNA', features, method='spearman') {
  exp_dat <- as.matrix(srt_inp[[assay]]@data)
  output_df <- data.frame(p_val=numeric(), cor=numeric(), p_val_adj=numeric(), feature=character(), gene=character())
  if (all(features %in% colnames(srt_inp@meta.data)))
  for (i in 1:length(features)) {
    feature_i <- features[i]
    cat('Calculating feature:', feature_i, '... \n')
    feature_dat <- srt_inp@meta.data[, feature_i]
    cor_test <- apply(exp_dat, 1, function(row_x) {
      cor_test_temp=cor.test(row_x, feature_dat, na.action='na.omit', method=method)
      c(p_val=as.numeric(cor_test_temp$p.value), cor=as.numeric(cor_test_temp$estimate))
    }) %>% t %>% data.frame
    cor_test$p_val_adj <- p.adjust(cor_test$p_val, method='bonferroni')
    cor_test$feature <- feature_i
    cor_test$gene <- gsub('(.*)\\..*', '\\1', rownames(exp_dat))
    output_df <- rbind(output_df, cor_test)
  }
  output_df %<>% na.omit
  output_df
}
