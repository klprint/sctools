

#### Functions ####
get.10x.data <- function(path, gzipped = "yes", mtx_name = "matrix.mtx", genes_name = "features.tsv", barcodes_name = "barcodes.tsv"){

  if(gzipped == "yes"){
    mtx_name <- paste0(mtx_name,".gz")
    genes_name <- paste0(genes_name, ".gz")
    barcodes_name <- paste0(barcodes_name, ".gz")
  }

  barcodes <- read.table(file.path(path, barcodes_name), comment.char = "%")
  genes <- read.table(file.path(path, genes_name), comment.char = "%", sep = "\t")
  mtx <- read.table(file.path(path, mtx_name), comment.char = "%", sep = " ")[-1,]

  A <- Matrix::sparseMatrix(i = mtx[,1], j = mtx[,2], x = mtx[,3])

  row.names(A) <- genes[1:range(mtx[,1])[2],1]
  colnames(A) <- barcodes[1:range(mtx[,2])[2],1]

  return(A)
}


generate.frac.intron.meta <- function(exon, fullTranscript){
  # if(!is.matrix(fullTranscript) & !is.matrix(exon)){
  #   ft <- fullTranscript@raw.data
  #   ex <- exon@raw.data
  # }

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


process.frac.intron.meta <- function(meta, max.cellrank = 10000,
                                     cutFracIntron = NULL, cutCellRank = NULL){
  meta <- meta %>% filter(CellRank <= max.cellrank)
  meta <- meta %>% add_column(Classification = as.character(Mclust(data.frame(.$CellRank, .$fracIntronic), G = 2)$classification))

  cl1.mean <- mean(meta$fracIntronic[meta$Classification == 1])
  cl2.mean <- mean(meta$fracIntronic[meta$Classification == 2])

  cl1.sd <- sd(meta$fracIntronic[meta$Classification == 1])
  cl2.sd <- sd(meta$fracIntronic[meta$Classification == 2])

  if(is.null(cutFracIntron)){
    cut.fracIntron <- max(cl1.mean, cl2.mean) - 2.5 * ifelse(cl1.mean > cl2.mean, cl1.sd, cl2.sd)
  }else{
    cut.fracIntron <- cutFracIntron
  }

  if(is.null(cutCellRank)){
    cut.maxCellRank <- max.cellrank
  }else{
    cut.maxCellRank <- cutCellRank
  }

  keep.cells <- meta %>% filter(fracIntronic >= cut.fracIntron,
                                CellRank <= cut.maxCellRank) %>%
    pull(CellID) %>% as.character()


  return(list(metadata = meta, good.cells = keep.cells))
}


scrublet <- function(mat, expected_doublet_rate = .1){
  print(dim(mat))

  scr <- reticulate::import("scrublet")

  scrub <- scr$Scrublet(as.matrix(t(mat)), expected_doublet_rate = .1)
  doublet_scores <- scrub$scrub_doublets()
  # hist(doublet_scores[[1]], breaks = 100)
  # abline(v = scrub$threshold_)

  print(
    ggplot2::ggplot(NULL, aes(x = doublet_scores[[1]])) +
      ggplot2::geom_histogram() +
      ggplot2::scale_y_log10() +
      ggplot2::geom_vline(xintercept = scrub$threshold_)
  )

  db.score <- doublet_scores[[1]]
  return(db.score)
}


get.gene.names <- function(genes, ensembl.mart){
  mart <- useMart("ensembl", ensembl.mart)

  gene.attrs <- c("ensembl_gene_id", "external_gene_name")

  mart.res <- getBM(gene.attrs, filters = "ensembl_gene_id", values = genes, mart = mart)

  return(mart.res)
}



#### MAIN ####

library(optparse)

option_list = list(
  make_option(c("-e", "--exon"),
              type="character",
              default=NULL,
              help="Exon Cellranger output",
              metavar="character"),
  make_option(c("-p", "--premrna"),
              type="character",
              default=NULL,
              help="pre-mRNA Cellranger output",
              metavar="character"),
  make_option(c("-o", "--outdir"),
              type="character",
              default=NULL,
              help="output folder",
              metavar="character"),
  make_option(c("-n", "--samplename"),
              type="character",
              default="scExperiment",
              help="Name of the sample",
              metavar="character"),
  make_option(c("-m", "--mart"),
              type="character",
              default="mmusculus_gene_ensembl",
              help="ENSEMBL mart string",
              metavar="character"),
  make_option(c("-z", "--gzipped"),
              action = "store_true",
              default=F,
              help="Add if the Cellranger outputfiles are gzipped",
              metavar="bool"),
  make_option(c("-c", "--mincells"),
              type="integer",
              default=10,
              help="Minimum number of cells needed to keep the gene",
              metavar="int"),
  make_option(c("-g", "--mingenes"),
              type="integer",
              default=200,
              help="Minimum number of genes needed to keep the cell",
              metavar="int"),
  make_option(c(NULL, "--doubletcut"),
              type="double",
              default=.1,
              help="SCRUBLET hardcut. Removes the top fraction of doublet scoring cells",
              metavar="double"),
  make_option(c(NULL, "--conda"),
              type="character",
              default="env_scanpy",
              help="Conda environment for SCRUBLET",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$exon) | is.null(opt$premrna) | is.null(opt$outdir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

library(reticulate)
library(tidyverse)
library(Seurat)
library(mclust)
library(Matrix)
library(biomaRt)


ex.path <- opt$exon
ft.path <- opt$premrna

outdir <- opt$outdir
sample.name <- opt$samplename
ensembl.mart <- opt$mart
gzipped <- ifelse(opt$gzipped, "yes", "no")
min.genes <- opt$mingenes
min.cells <- opt$mincells
doublet.hardcut <- opt$doubletcut
conda_env <- opt$conda

use_condaenv(conda_env)
py_config()

dir.create(outdir, showWarnings = F)

cat("Reading the data\n")
cat("Reading exonic UMI\n")
ex <- get.10x.data(ex.path, gzipped = gzipped)
cat("Reading pre-mRNA UMI\n")
ft <- get.10x.data(ft.path, gzipped = gzipped)

cat("Selecting good droplets\n")
meta <- generate.frac.intron.meta(ex, ft)
meta <- process.frac.intron.meta(meta)

plot(meta$metadata$CellRank,
     meta$metadata$fracIntronic,
     col = meta$metadata$Classification)


# Samples cells according to genes
mat.ex <- ex[,colnames(ex) %in% meta$good.cells]
mat.ft <- ft[,colnames(ft) %in% meta$good.cells]

mat.ex <- mat.ex[rowSums(mat.ex>0) >= min.cells,
                 colSums(mat.ex>0) >= min.genes]

mat.ft <- mat.ft[rowSums(mat.ft>0) >= min.cells,
                 colSums(mat.ft>0) >= min.genes]

db.score <- scrublet(mat.ft)
not.db <- db.score < quantile(db.score, 1-doublet.hardcut)

mat.ex <- mat.ex[,colnames(mat.ex) %in% colnames(mat.ft)[not.db]]
mat.ft <- mat.ft[,which(not.db)]

colnames(mat.ex) <- paste0(sample.name, "_", colnames(mat.ex))
colnames(mat.ft) <- paste0(sample.name, "_", colnames(mat.ft))

gz.ex <- gzfile(file.path(outdir, "UMI_exon.csv.gz"))
gz.ft <- gzfile(file.path(outdir, "UMI_pre_mRNA.csv.gz"))

cat("Translating ENSEMBL IDs\n")
mart.ex <- get.gene.names(rownames(mat.ex), ensembl.mart)
mart.ft <- get.gene.names(rownames(mat.ft), ensembl.mart)

gnames.ex <- mart.ex$external_gene_name[match(rownames(mat.ex), mart.ex$ensembl_gene_id)]
gnames.ex[is.na(gnames.ex)] <- paste0("NA_", 1:sum(is.na(gnames.ex)))

gnames.ft <- mart.ft$external_gene_name[match(rownames(mat.ft), mart.ft$ensembl_gene_id)]
gnames.ft[is.na(gnames.ft)] <- paste0("NA_", 1:sum(is.na(gnames.ft)))

cat("Prepping output\n")
mat.ex <- as.data.frame(as.matrix(mat.ex)) %>% rownames_to_column("ENSEMBL_ID") %>% add_column("GENE_SYMBOL" = gnames.ex, .after = "ENSEMBL_ID")
mat.ft <- as.data.frame(as.matrix(mat.ft)) %>% rownames_to_column("ENSEMBL_ID") %>% add_column("GENE_SYMBOL" = gnames.ft, .after = "ENSEMBL_ID")

cat(paste0("Saving output to ", outdir, "\n"))
write.table(mat.ex, gz.ex, quote = F, row.names = F, sep=" ")
write.table(mat.ft, gz.ft, quote = F, row.names = F, sep=" ")
