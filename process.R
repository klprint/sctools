#!/usr/bin/env Rscript
library(reticulate)
use_condaenv("env_scanpy")
py_config()

library(tidyverse)
library(mclust)
library(Matrix)
library(hdf5r)


get.10x.data <- function(path, sample.name, V = 2, mtx_name = "matrix.mtx", genes_name = "genes.tsv", barcodes_name = "barcodes.tsv"){

  if(V == 3){
    cat("Reding V3 data\n")
    mtx_name <- "matrix.mtx.gz"
    genes_name <- "features.tsv.gz"
    barcodes_name <- "barcodes.tsv.gz"
  }

  barcodes <- read.table(file.path(path, barcodes_name), comment.char = "%")
  genes <- read.table(file.path(path, genes_name), comment.char = "%", sep = "\t")
  mtx <- read.table(file.path(path, mtx_name), comment.char = "%", sep = " ")[-1,]

  A <- Matrix::sparseMatrix(i = mtx[,1], j = mtx[,2], x = mtx[,3])

  row.names(A) <- genes[1:range(mtx[,1])[2],1]
  colnames(A) <- barcodes[1:range(mtx[,2])[2],1]

  colnames(A) <- paste0(sample.name, "_", colnames(A))

  return(A)
}

generate.frac.intron.meta <- function(ex, ft){
  keep.cells <- intersect(colnames(ft), colnames(ex))

  cat("Number of cells in both sets:",length(keep.cells), "\n")
  ex <- ex[,keep.cells]
  ft <- ft[,keep.cells]

  ex.sum <- as.data.frame(Matrix::colSums(ex))
  colnames(ex.sum) = "exonic"

  ft.sum <- as.data.frame(Matrix::colSums(ft))
  colnames(ft.sum) <- "fullTranscript"

  cells.meta <- merge(ft.sum, ex.sum, by = 0)
  colnames(cells.meta)[1] <- "CellID"

  cells.meta <- dplyr::filter(cells.meta, fullTranscript > 0)
  cells.meta <- cells.meta %>% arrange(desc(fullTranscript)) %>%
    mutate(CellRank = 1:nrow(cells.meta)) %>%
    mutate(fracIntronic = 1 - (exonic / fullTranscript))

  return(cells.meta)
}

scrublet <- function(umi, expected_doublet_rate = .1){

  scr <- reticulate::import("scrublet")

  scrub <- scr$Scrublet(as.matrix(t(umi)), expected_doublet_rate = .1)
  doublet_scores <- scrub$scrub_doublets()
  # hist(doublet_scores[[1]], breaks = 100)
  # abline(v = scrub$threshold_)

  # print(
  #   ggplot2::ggplot(NULL, aes(x = doublet_scores[[1]])) +
  #     ggplot2::geom_histogram() +
  #     ggplot2::scale_y_log10() +
  #     ggplot2::geom_vline(xintercept = scrub$threshold_)
  # )

  scrub.meta <- doublet_scores[[1]]
  return(scrub.meta)
}

make.analysis <- function(ex.path, ft.path, sample.name, V=2,
                          empty.cut = 50,
                          max.cellrank = 1e4,
                          man.select = list(cutFracIntron = NULL,
                                            maxCellRank = NULL),
                          min.cells = 25, min.genes=500){
  cat("Reading exon\n")
  ex <- get.10x.data(ex.path, sample.name, V=V)
  cat("Reading pre-mRNA\n")
  ft <- get.10x.data(ft.path, sample.name, V=V)
  empty.drops <- ft

  cat("Exon: ", ncol(ex), " ", nrow(ex), "\n")
  cat("Full Transcript: ", ncol(ft), " ", nrow(ft), "\n")

  cat("Selecting good barcodes\n")
  meta <- generate.frac.intron.meta(ex, ft)

  meta <- meta %>% filter(CellRank <= max.cellrank)

  meta <- meta %>% add_column(Classification = as.character(Mclust(data.frame(.$CellRank, .$fracIntronic), G = 2)$classification))


  cl1.mean <- mean(meta$fracIntronic[meta$Classification == 1])
  cl2.mean <- mean(meta$fracIntronic[meta$Classification == 2])

  cl1.sd <- sd(meta$fracIntronic[meta$Classification == 1])
  cl2.sd <- sd(meta$fracIntronic[meta$Classification == 2])

  if(is.null(man.select$cutFracIntron)){
    cut.fracIntron <- max(cl1.mean, cl2.mean) - 2.5 * ifelse(cl1.mean > cl2.mean, cl1.sd, cl2.sd)
  }else{
    cut.fracIntron <- man.select$cutFracIntron
  }

  if(is.null(man.select$maxCellRank)){
    maxCellRank <- max.cellrank
  }else{
    maxCellRank <- man.select$maxCellRank
  }

  keep.cells <- meta %>% filter(fracIntronic >= cut.fracIntron, CellRank <= maxCellRank) %>% pull(CellID) %>% as.character()

  ft <- ft[,keep.cells]
  ex <- ex[,keep.cells]


  cat("Removing doublets\n")
  ft.scrub <- scrublet(ft)
  ex.scrub <- scrublet(ex)

  ft.scrub.cut <- sort(ft.scrub, decreasing = T, na.last = T)[floor(length(ft.scrub) * .1)]
  ex.scrub.cut <- sort(ex.scrub, decreasing = T, na.last = T)[floor(length(ex.scrub) * .1)]

  # print(ft.scrub.cut)
  # print(length(ft.scrub))
  # print(head(ft.scrub))
  # print(sum(ft.scrub <= ft.scrub.cut))
  # print(dim(ft))

  ft <- ft[,colnames(ft)[ft.scrub <= ft.scrub.cut]]
  ex <- ex[,colnames(ex)[ex.scrub <= ex.scrub.cut]]


  ft <- ft[rowSums(ft>0) >= min.cells,
           colSums(ft>0) >= min.genes]
  ex <- ex[,colnames(ex) %in% colnames(ft)]

  empty.drops <- empty.drops[rownames(ft),
                             colSums(empty.drops>0) < empty.cut & colSums(empty.drops>0) != 0]

  return(list(FT=ft, EXON=ex, META=meta,
              empty = empty.drops))
}




remove.background <- function(umi, empty.drops, min.cells = 25, min.genes=500){
  stopifnot(sum(rownames(umi) == rownames(empty.drops)) == nrow(umi))

  keep.cells <- colnames(umi)[colSums(umi>0) >= min.genes]
  keep.genes <- rownames(umi)[rowSums(umi>0) >= min.cells]

  umi <- umi[keep.genes,]
  empty.drops <- empty.drops[keep.genes,]

  empty.pseudobulk <- matrix(rowSums(empty.drops), ncol=1)
  print(head(empty.pseudobulk))
  rownames(empty.pseudobulk) <- rownames(empty.drops)

  cat("Normalizing \n")
  exprs <- apply(umi, 2, function(x) sqrt(x / sum(x)))
  exprs.drops <- apply(empty.pseudobulk, 2, function(x) sqrt(x / sum(x)))
  exprs.drops <- as.matrix(exprs.drops, ncol=1)

  rownames(exprs.drops) <- rownames(empty.pseudobulk)

  cat("Calculating dot product \n")
  dp <- t(exprs) %*% exprs.drops
  print(head(dp))

  cat("Linear regression\n")
  lmout <- lm(as.matrix(t(exprs)) ~ as.matrix(dp))$residuals
  return(t(lmout))
}


find.hvg <- function(exprs, is.umi = T, min.genes=500, min.cells=25){

  if(is.umi){
    exprs <- exprs[rowSums(exprs>0)>=min.cells,
                   colSums(exprs>0)>=min.genes]

    exprs <- exprs / colSums(exprs)
  }else{
    exprs <- (10^exprs)-1
  }


  means <- rowMeans(exprs)
  vars <- apply(exprs, 1, var)
  disp <- vars/means

  cutoff <- mean(disp,na.rm=T) + sd(disp,na.rm = T)

  plot(means, disp, log="xy",pch=".")
  abline(h=cutoff)

  return(rownames(exprs)[disp >= cutoff])

}

run.umap <- function(exprs, hvg = NULL, n_dims=2){
  umap <- import("umap")

  if(!is.null(hvg)){
    exprs <- exprs[hvg,]
  }

  #u <- umap::umap(t(exprs), method="umap-learn", metric="correlation", min_dist=.3, n_neighbors=30, n_components=n_dims)
  u <- umap$UMAP(metric = "correlation", min_dist = .3, n_neighbors = 30L, n_components = as.integer(n_dims))$fit_transform(t(exprs))
  rownames(u) <- colnames(exprs)
  print(head(u))
  return(u)
}

union.merge <- function(mergelist, all.rn = NULL){
  if(is.null(all.rn)){
    cat("Getting all genes\n")
    all.rn <- NULL
    for(item in mergelist){
      all.rn <- union(rownames(item), all.rn)
    }
  }

  cat("Merging\n")

  mergelist <- lapply(mergelist, function(x){
    add.x <- all.rn[!(all.rn %in% rownames(x))]

    if(length(add.x) < 1){
      return(x[all.rn,])
    }else{
      add.mtx.x <- matrix(0, nrow=length(add.x), ncol=ncol(x))
      rownames(add.mtx.x) <- add.x
      o <- rbind(x, add.mtx.x)
      return(o[all.rn, ])
    }


  })

  out <- do.call(cbind, mergelist)

  return(out)
}

full.analysis <- function(ex.path, ft.path, sample.name, V = 2){
  sample.list <- make.analysis(ex.path, ft.path, sample.name, V = V)
  exprs <- remove.background(sample.list$FT, sample.list$empty)
  hvg <- find.hvg(sample.list$FT[rownames(exprs), ])

  return(list(umi = sample.list, exprs=exprs, hvg=hvg))
}

plot.umap <- function(umap.coords){
  factors <- do.call("c", lapply(rownames(umap.coords), function(x) strsplit(x, "_")[[1]][1]))

  return(
    as_tibble(umap.coords) %>%
      ggplot(aes(x=V1, y=V2)) +
      geom_point(aes(color=factors), size=.25)
  )
}

#### MAIN ####
args = commandArgs(trailingOnly=TRUE)

# 1. Exon output path
# 2. Full transcript output path
# 3. Sample name
# 4. 10x Chromium version
# 5. Output main folder
print(args)
if( length( args ) < 3  ){
  stop("At least three arguments need to be provided", call. = F)
}else if( length( args ) == 3  ){
  args[4] <- "2"
}

args[4] <- as.numeric(args[4])
args[5] <- "bg_removed"

out.dir <- file.path(args[5], args[3])
dir.create(out.dir, recursive = T, showWarnings = F)

pdf(file.path(out.dir, "plots.pdf"))
dfl <- full.analysis(args[1], args[2], args[3], args[4])

cat("Running 2D UMAP\n")
umap2d.dfl <- run.umap(dfl$exprs[dfl$hvg, ])

cat("Running 3D UMAP\n")
umap3d.dfl <- run.umap(dfl$exprs[dfl$hvg, ], n_dims = 3)

plot.umap(umap2d.dfl)

cat("Plotting statistics\n")
ggplot(
  NULL,
  aes(x = colSums(dfl$umi$FT))
) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  xlab("log(nUMI)") +
  ggtitle(paste0(args[3], " nUMI FT"))

ggplot(
  NULL,
  aes(x = colSums(dfl$umi$EXON))
) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  xlab("log(nUMI)") +
  ggtitle(paste0(args[3], " nUMI exonic"))

ggplot(
  NULL,
  aes(x = colSums(dfl$umi$FT>0))
) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  xlab("log(nGene)") +
  ggtitle(paste0(args[3], " nGene FT"))

ggplot(
  NULL,
  aes(x = colSums(dfl$umi$EXON>0))
) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  xlab("log(nGene)") +
  ggtitle(paste0(args[3], " nGene exonic"))

dev.off()

# cat("Writing UMAPS\n")
umap2d <- as.data.frame(umap2d.dfl)
colnames(umap2d) <- c("UMAP1", "UMAP2")
umap2d <- umap2d %>% rownames_to_column("CellID")
# write.table(umap2d, file=file.path(out.dir, "2D_UMAP.csv"), quote=F, sep=",", col.names = T, row.names = F)
umap3d <- as.data.frame(umap3d.dfl)
colnames(umap3d) <- c("UMAP1", "UMAP2", "UMAP3")
umap3d <- umap3d %>% rownames_to_column("CellID")
# write.table(umap3d, file=file.path(out.dir, "3D_UMAP.csv"), quote=F, sep=",", col.names = T, row.names = F)
#
# cat("Writing the output files\n")
# mat.outdir <- file.path(out.dir, "matrices")
# dir.create(mat.outdir, showWarnings = F, recursive = T)
#
# cat(" Writing full transcript\n")
# ft.mtx <- gzfile(file.path(mat.outdir, "full_transcript.csv.gz"))
# write.table(as.matrix(dfl$umi$FT), ft.mtx, sep=",", quote = F)
#
# cat(" Writing exon only\n")
# ex.mtx <- gzfile(file.path(mat.outdir, "exon.csv.gz"))
# write.table(as.matrix(dfl$umi$EXON), ex.mtx, sep=",", quote = F)
#
# cat(" Writing HVG list\n")
# write.table(data.frame(dfl$hvg), file=file.path(out.dir, "HVG.csv"), quote = F, sep=",", row.names = F, col.names = F)
#

cat("Writing data")

h5file <- H5File$new(file.path(out.dir, "data.h5"), mode="w")

# Creating file structure
h5file$create_group("umi")
h5file[["umi"]]$create_group("ft")
h5file[["umi"]]$create_group("exon")

h5file$create_group("umap")

h5file[["umi"]][["ft"]][["counts"]] <- as.matrix(dfl$umi$FT)
h5file[["umi"]][["ft"]][["genes"]] <- rownames(dfl$umi$FT)
h5file[["umi"]][["ft"]][["cells"]] <- colnames(dfl$umi$FT)
h5file[["umi"]][["ft"]][["sf"]] <- colSums(dfl$umi$FT)
h5file[["umi"]][["ft"]][["hvg"]]<- dfl$hvg

h5file[["umi"]][["exon"]][["counts"]] <- as.matrix(dfl$umi$EXON)
h5file[["umi"]][["exon"]][["genes"]] <- rownames(dfl$umi$EXON)
h5file[["umi"]][["exon"]][["cells"]] <- colnames(dfl$umi$EXON)
h5file[["umi"]][["exon"]][["sf"]] <- colSums(dfl$umi$EXON)

h5file[["umap"]][["2d"]] <- umap2d
h5file[["umap"]][["3d"]] <- umap3d

h5file$create_group("exprs")
h5file[["exprs"]][["mat"]] <- as.matrix(dfl$exprs)
h5file[["exprs"]][["genes"]] <- rownames(dfl$exprs)
h5file[["exprs"]][["cells"]] <- colnames(dfl$exprs)


h5file$close_all()
write.table(umap2d, file=file.path(out.dir, "2D_UMAP.csv"), quote=F, sep=",", col.names = T, row.names = F)
write.table(umap3d, file=file.path(out.dir, "3D_UMAP.csv"), quote=F, sep=",", col.names = T, row.names = F)

