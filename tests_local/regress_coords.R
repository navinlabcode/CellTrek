## Regression check after the coordinate-scale fix: traint must resolve coords for
## both VisiumV1 (brain_st_cortex, via @coordinates) and VisiumV2 (stxBrain, via GTC x/y).
suppressPackageStartupMessages({library(CellTrek); library(Seurat); library(dplyr)})
set.seed(1)
mini_sc <- function(sc, per=25){ md<-sc@meta.data; k<-unlist(lapply(split(rownames(md), md$cell_type), function(x) sample(x, min(length(x),per)))); sc[,k] }

## ---- VisiumV1 ----
st1 <- readRDS('data/brain_st_cortex.rds'); sc1 <- readRDS('data/brain_sc.rds')
st1 <- RenameCells(st1, new.names=make.names(Cells(st1))); sc1 <- RenameCells(sc1, new.names=make.names(Cells(sc1)))
st1 <- st1[, sample(Cells(st1), 300)]; sc1 <- mini_sc(sc1, 25)
t1 <- CellTrek::traint(st1, sc1, st_assay='Spatial', sc_assay='RNA', cell_names='cell_type')
c1 <- na.omit(t1@meta.data[t1$type=='st', c('coord_x','coord_y')])
cat('VisiumV1: image', class(st1@images[[1]]), '| coord range', paste(round(range(c1$coord_x)),collapse='-'),
    '| spot_dis', round(median(dbscan::kNN(c1,k=6)$dist)), '| st coords set', all(!is.na(t1$coord_x[t1$type=='st'])), '\n')

## ---- VisiumV2 ----
e <- new.env(); load('data2/anterior1.rda', envir=e); st2 <- UpdateSeuratObject(get('anterior1', envir=e))
st2 <- NormalizeData(st2, verbose=FALSE) %>% FindVariableFeatures(nfeatures=2000, verbose=FALSE)
st2 <- RenameCells(st2, new.names=make.names(Cells(st2)))
sc2 <- readRDS('data2/allen_cortex.rds'); sc2$cell_type <- as.character(sc2$subclass)
sc2 <- mini_sc(sc2, 25)
sc2 <- NormalizeData(sc2,verbose=FALSE) %>% FindVariableFeatures(nfeatures=2000,verbose=FALSE) %>% ScaleData(verbose=FALSE) %>% RunPCA(npcs=30,verbose=FALSE) %>% RunUMAP(dims=1:30,verbose=FALSE)
sc2 <- RenameCells(sc2, new.names=make.names(Cells(sc2)))
t2 <- CellTrek::traint(st2, sc2, st_assay='Spatial', sc_assay='RNA', cell_names='cell_type')
c2 <- na.omit(t2@meta.data[t2$type=='st', c('coord_x','coord_y')])
cat('VisiumV2: image', class(st2@images[[1]]), '| coord range', paste(round(range(c2$coord_x)),collapse='-'),
    '| spot_dis', round(median(dbscan::kNN(c2,k=6)$dist)), '| st coords set', all(!is.na(t2$coord_x[t2$type=='st'])), '\n')
cat('== REGRESS OK ==\n')
