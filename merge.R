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
i <- 1
for(s in sample.folders){
  cat("Reading ", s, "\n")
  h5file <- H5File$new(file.path(s, "data.h5"), mode = "r")
  hvgs <- union(hvgs, h5file[["umi"]][["ft"]][["hvg"]][])
  all.genes <- h5file[["exprs"]][["genes"]][]

  tmp <- as.matrix(h5file[["exprs"]][["mat"]][
    all.genes %in% hvgs,
  ])

  colnames(tmp) <- h5file[["exprs"]][["cells"]][]

  h5file$close_all()

  rownames(tmp) <- all.genes[
    all.genes %in% hvgs
  ]

  exprs.list[[i]] <- tmp
  print(exprs.list[[i]][1:5,1:5])

  i <- i+1
}

cat("Number of highly variable genes read from files: ", length(hvgs), "\n")

cat("Merging the matrices\n")
exprs <- union.merge(exprs.list, hvgs)


cat("Running UMAP\n")
cat("\t 2D\n")
umap2d <- run.umap(exprs)
cat("\t 3D\n")
umap3d <- run.umap(exprs, n_dims=3)
umap3d <- as.data.frame(umap3d) %>% rownames_to_column("CellID")
colnames(umap3d) <- c("CellID", "UMAP1", "UMAP2", "UMAP3")


umap3d <- umap3d %>% add_column(Batch = do.call("c", lapply(umap3d$CellID, function(x) strsplit(x, "_")[[1]][2])))

write.table(umap3d, file=file.path(outputfolder, "UMAP3d.csv"), sep=",", quote = F, row.names = F)
