## Lightweight validation of the new traint(coord_df=) parameter on a normal,
## fully-operable object (brain_st_cortex, VisiumV1). Small subsample for speed.
suppressPackageStartupMessages({library(CellTrek); library(Seurat); library(dplyr)})
set.seed(1)
st <- readRDS('data/brain_st_cortex.rds')
sc <- readRDS('data/brain_sc.rds')
st <- RenameCells(st, new.names = make.names(Cells(st)))
sc <- RenameCells(sc, new.names = make.names(Cells(sc)))
## subsample for speed
st <- st[, sample(Cells(st), 250)]
set.seed(2); md <- sc@meta.data
keep <- unlist(lapply(split(rownames(md), md$cell_type), function(x) sample(x, min(length(x), 25))))
sc <- sc[, keep]

## coord_df with NON-standard column names (col1 -> coord_x, col2 -> coord_y)
tc <- GetTissueCoordinates(st)
my_coord <- data.frame(rowpos = tc[, 'imagerow'], colpos = tc[, 'imagecol'])
rownames(my_coord) <- Cells(st)

cat('== Test 1: traint WITH coord_df (non-standard col names) ==\n')
tr <- CellTrek::traint(st, sc, st_assay = 'Spatial', sc_assay = 'RNA',
                       cell_names = 'cell_type', coord_df = my_coord)
st_cells <- tr$id[tr$type == 'st']
ok <- isTRUE(all.equal(unname(tr$coord_x[tr$type=='st']),
                       unname(my_coord[st_cells, 1]), tolerance = 1e-6))
cat('coord_x matches coord_df col1:', ok, '\n')
cat('traint cells:', ncol(tr), '\n')

cat('\n== Test 2: bad rownames -> expect graceful error ==\n')
bad <- my_coord; rownames(bad) <- paste0('X', seq_len(nrow(bad)))
r2 <- tryCatch({ CellTrek::traint(st, sc, sc_assay='RNA', cell_names='cell_type', coord_df=bad); 'NO ERROR (BAD)'},
               error=function(e) paste0('ERRORED as expected: ', conditionMessage(e)))
cat(r2, '\n')

cat('\n== Test 3: auto-detect still works WITHOUT coord_df ==\n')
tr2 <- CellTrek::traint(st, sc, st_assay='Spatial', sc_assay='RNA', cell_names='cell_type')
cat('auto-detect coord_x set:', all(!is.na(tr2$coord_x[tr2$type=='st'])), '\n')

cat('\n== COORD_DF PARAM TEST DONE ==\n')
