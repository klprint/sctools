#!/usr/bin/env Rscript
library(reticulate)
use_condaenv("env_scanpy")
py_config()

library(tidyverse)
library(mclust)
library(Matrix)

run.umap <- function(exprs, hvg = NULL, n_dims=2){
  if(!is.null(hvg)){
    exprs <- exprs[hvg,]
  }

  u <- umap::umap(t(exprs), method="umap-learn", metric="correlation", min_dist=.3, n_neighbors=30, n_dimensions=n_dims)
  return(u$layout)
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

cat("Reding the data\n")

hvgs <- NULL
for(s in sample.folders){
  hvgs <- union(hvgs, read.table(file.path(s, "HVG.csv"))[,1])
}

print(length(hvgs))
