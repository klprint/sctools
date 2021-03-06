#!/usr/bin/env Rscript
library(reticulate)
use_condaenv("env_scanpy")
py_config()

library(tidyverse)
library(mclust)
library(Matrix)
library(hdf5r)

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

plot.umap <- function(umap.coords){
  factors <- do.call("c", lapply(rownames(umap.coords), function(x) strsplit(x, "_")[[1]][1]))

  return(
    as_tibble(umap.coords) %>%
      ggplot(aes(x=V1, y=V2)) +
      geom_point(aes(color=factors), size=.25)
  )
}

remove.background <- function(exprs, empty.drops){
  stopifnot(sum(rownames(exprs) == rownames(empty.drops)) == nrow(exprs))

  exprs.drops <- matrix(sqrt(empty.drops/colSums(empty.drops)), ncol=1)
  rownames(exprs.drops) <- rownames(empty.drops)
  cat("Calculating dot product \n")
  dp <- t(exprs) %*% exprs.drops
  print(head(dp))

  cat("Linear regression\n")
  lmout <- lm(as.matrix(t(exprs)) ~ as.matrix(dp))$residuals
  return(t(lmout))
}

#### MAIN ####
args = commandArgs(trailingOnly=TRUE)

# 1. Comma separated list of "process.R" output folders
if( length( args ) < 1  ){
  stop("At least one arguments need to be provided", call. = F)
}

sample.folders <- strsplit(args[1], ",")[[1]]
outputfolder <- args[2]

if(is.na(outputfolder)){
  outputfolder <- "bgr_merged"
}

dir.create(outputfolder)
cat("Reading the data\n")

hvgs <- NULL
exprs.list <- list()
umi.list <- list()
i <- 1
for(s in sample.folders){
  cat("Reading ", s, "\n")
  h5file <- H5File$new(file.path(s, "data.h5"), mode = "r")

  hvgs <- all.genes <- h5file[["umi"]][["ft"]][["genes"]][]
  cids <- h5file[["umi"]][["ft"]][["cells"]][]

  missing.genes <- hvgs[!(hvgs %in% all.genes)]
  present.genes <- hvgs[hvgs %in% all.genes]

  cat("  empty droplets\n")
  empty.genes <- h5file[["umi"]][["empty"]][["genes"]][]
  empty.drops <- matrix(h5file[["umi"]][["empty"]][["counts"]][empty.genes %in% present.genes], ncol=1)
  rownames(empty.drops) <- empty.genes[empty.genes %in% present.genes]

  if(length(missing.genes) > 0){
    filling.mat <- matrix(0, nrow = length(missing.genes),
                          ncol=length(cids))
    rownames(filling.mat) <- missing.genes
  }

  cat("  UMI matrix\n")
  tmp <- as.matrix(h5file[["umi"]][["ft"]][["counts"]][
    all.genes %in% present.genes,
    ])
  rownames(tmp) <- all.genes[all.genes %in% present.genes]

  if(length(missing.genes) > 0){
    cat("Filling genes (",length(missing.genes),")\n")
    tmp <- rbind(tmp, filling.mat)
    empty.drops <- matrix(c(empty.drops[,1], rep(0, length(missing.genes))), ncol=1)
    rownames(empty.drops) <- c(present.genes, missing.genes)
  }

  tmp <- tmp[hvgs, ]
  empty.drops <- matrix(empty.drops[hvgs, ],ncol=1)
  rownames(empty.drops) <- hvgs
  cat("------------------\n")
  print(head(tmp[,1:5]))
  print(head(empty.drops))
  cat("------------------\n")

  colnames(tmp) <- h5file[["umi"]][["ft"]][["cells"]][]
  sf <- h5file[["umi"]][["ft"]][["sf"]][]

  h5file$close_all()

  cat("Normalizing\n")

  umi.list[[i]] <- tmp

  tmp <- sqrt(tmp/sf)

  exprs.list[[i]] <- tmp


  i <- i+1
}

cat("Number of genes read from files: ", length(hvgs), "\n")


cat("Merging the matrices\n")
exprs <- union.merge(exprs.list, hvgs)
umi <- union.merge(umi.list, hvgs)
cat("Number of cells: ", ncol(exprs), "\n")

cat("Writing h5 file\n")
h5file <- H5File$new(file.path(outputfolder, "Merged_Matrix.h5"), mode="w")
h5file[["umi"]] <- umi
h5file[["exprs"]] <- exprs
h5file[["cells"]] <- colnames(umi)
h5file[["genes"]] <- rownames(umi)
h5file$close_all()

gzout <- gzfile(file.path(outputfolder, "Merged_Matrix.csv.gz"))
write.table(exprs, gzout, sep=",", quote = F)
