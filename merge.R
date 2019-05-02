#! /usr/bin/env Rscript
library(reticulate)
use_condaenv("env_scanpy")
py_config()

library(optparse)
library(hdf5r)
library(Matrix)
library(tidyverse)

optList <- list(
  make_option(
    c("-i", "--input"),
    type="character",
    help="Comma separated input list",
    default=NULL,
    metavar="character"
  ),

  make_option(
    c("-o", "--output"),
    type="character",
    help="Output folder path",
    default=NULL,
    metavar = "character"
  ),

  make_option(
    c("-n", "--name"),
    type="character",
    help = "Output file name",
    default="data",
    metavar = "character"
  ),

  make_option(
    c("-t", "--threads"),
    type="numeric",
    help="Number of cores",
    default=3,
    metavar = "integer"
  )
)

opt_parser <- OptionParser(option_list = optList)
opt = parse_args(opt_parser)


#####################
##### FUNCTIONS #####
#####################
rem.batch <- function(mtx, l, ncores=3){
  out <- pbmcapply::pbmclapply(unique(l), function(x){
    tmp.mtx <- as.matrix(mtx[,l == x])
    rs <- rowSums(mtx)
    names(rs) <- rownames(mtx)

    tmp.mtx <- sqrt(t(
      t(tmp.mtx) / colSums(tmp.mtx)
    ))

    rs <- sqrt(rs / sum(rs))

    lm(tmp.mtx~rs)$residuals
  }, mc.cores = ncores)

  out <- do.call(cbind,
                 out)

  return(
    out
  )
}

dim.reduce <- function(mtx, n = 100){
  cat("Running PCA\n")
  pca <- irlba::prcomp_irlba(mtx, n = n)
  rownames(pca$x) <- rownames(mtx)
  return(pca$x)
}

get.synbulk <- function(mtx, l){
  synbulk <- lapply(unique(l), function(lsub){
    tmp <- mtx[,l == lsub]
    rowSums(tmp)
  })

  tmp <- do.call(cbind, synbulk)
  colnames(tmp) <- unique(l)
  return(tmp)
}

#########################
if(is.null(opt$input)){
  stop("Please specify input using the '-i' flag.")
}

if(is.null(opt$output)){
  stop("Please specify an output folder using the '-o' flag.")
}

infiles <- opt$input
infiles <- strsplit(infiles, ",")[[1]]

use_genes <- NULL
use_cells <- NULL

meta.data <- list()
meta.data$batches <- do.call("c",lapply(infiles, function(x) rev(strsplit(x, "/")[[1]])[2]))
names(meta.data$batches) <- infiles

cat("Reading cell and gene annotations\n")
for(f in infiles){
  cat("\t file:", f, "\n")
  h5 <- H5File$new(f, mode="r")
  use_genes <- unique(union(use_genes, h5[["umi/ft/genes"]][]))
  use_cells <- c(use_cells, h5[["umi/ft/cells"]][])
  meta.data$hvg <- union(meta.data$hvg, h5[["umi/ft/hvg"]][])


  meta.data$batch <- c(meta.data[["batch"]],
                       rep(meta.data$batches[f], length(h5[["umi/ft/cells"]][])))
  h5$close_all()
}

names(meta.data$batch) <- NULL

mtx <- Matrix(0, ncol = length(use_cells),
              nrow = length(use_genes),
              dimnames = list(use_genes, use_cells),
              sparse = T)

cat("\nReading matrix of the following dimensions:\n")
print(dim(mtx))

##### Reading the data #####

cat("\nReading the matrices\n")
# for(f in infiles){
#   cat("\t file:", f, "\n")
#   h5 <- H5File$new(f, mode="r")
#   tmp_genes <- h5[["umi/ft/genes"]][]
#   tmp_cells <- h5[["umi/ft/cells"]][]
#   tmp_mtx <- Matrix(h5[["umi/ft/counts"]][,], sparse=T,
#                     dimnames = list(tmp_genes, tmp_cells))
#
#   h5$close_all()
#
#   mtx[tmp_genes, tmp_cells] <- tmp_mtx
# }

tmp <- pbmcapply::pbmclapply(infiles, function(f){
  h5 <- H5File$new(f, mode="r")
  tmp_genes <- h5[["umi/ft/genes"]][]
  tmp_cells <- h5[["umi/ft/cells"]][]
  tmp_mtx <- Matrix(h5[["umi/ft/counts"]][,], sparse=T,
                    dimnames = list(tmp_genes, tmp_cells))

  h5$close_all()

  return(tmp_mtx)
}, mc.cores = opt$threads)
names(tmp) <- infiles

cat("\tParsing\n")
for(f in names(tmp)){
  cat("\t", f, "\n")
  mtx[rownames(tmp[[f]]),
      colnames(tmp[[f]])] <- tmp[[f]]
}

rm(tmp)

##### Writing UMI matrix #####
dir.create(opt$output, showWarnings = F)
cat("Writing sparse matrix rds file\n\n")

saveRDS(mtx, file=file.path(opt$output, paste0(opt$name, ".rds")))


##### Stat plots #####
pdf(file.path(
  opt$output,
  "stat_plots.pdf"),
  width=7, height=4
)

ggplot() +
  geom_density(aes(x=colSums(mtx),
                   fill=meta.data$batch,
                   group=meta.data$batch),
               alpha=.5) +
  scale_fill_discrete(name="Run") +
  xlab("nUMI")

ggplot() +
  geom_density(aes(x=colSums(mtx>0),
                   fill=meta.data$batch,
                   group=meta.data$batch),
               alpha=.5) +
  scale_fill_discrete(name="Run") +
  xlab("nGenes")

ggplot(NULL,
       aes(x = colSums(mtx),
           y = colSums(mtx > 0),
           color = meta.data$batch)) +
  geom_point(
    size=.1
  ) +
  geom_rug(alpha=.1) +
  xlab("nUMI") +
  ylab("nGenes") +
  scale_color_discrete(name="Run") +
  scale_x_log10() +
  scale_y_log10()

pheatmap::pheatmap(
  cor(
    get.synbulk(mtx, meta.data$batch)
  )
)

dev.off()


##### Processing data #####
cat("Normalizing\n")
br.mtx <- t(rem.batch(mtx[rownames(mtx) %in% meta.data$hvg, ],
                      meta.data$batch, opt$threads)) # cells x genes

if(ncol(br.mtx) > 400){
  br.mtx <- dim.reduce(br.mtx, n=400)
}


cat("Running UMAP: ")
cat("2D ")
meta.data$umap2d <- umap::umap(br.mtx, method="umap-learn", metric="correlation",
                               min_dist=.3, n_neighbors=30)
cat("Clusterable (15D) ")
meta.data$umap15d <- umap::umap(br.mtx, method="umap-learn", metric="correlation",
                                min_dist=0.0, n_neighbors=30, n_components=15)

cat("3D\n")
meta.data$umap3d <- umap::umap(br.mtx, method="umap-learn", metric="correlation",
                               min_dist=0.0, n_neighbors=30, n_components=3)

cat("Writing meta data\n")
saveRDS(meta.data, file=file.path(opt$output,
                                  paste0(opt$name,"_metadata.rds")))
