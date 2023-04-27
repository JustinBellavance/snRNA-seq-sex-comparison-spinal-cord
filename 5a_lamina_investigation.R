library(Seurat)
library(SingleCellExperiment)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
library(tidyverse)

#specific lamina study
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#seurat_small <- readRDS("../results/results_small/query.rds")
seurat <- readRDS("../Results/results_small/query.rds")

seurat <- SetIdent(seurat, value = seurat@meta.data$predicted.id)

unique(Idents(seurat) )

#subset SDH neurons.
seurat_SDH <- subset(seurat, idents = c("Excit-01", "Excit-02", "Excit-03", "Excit-06", "Excit-08", "Excit-10", "Excit-12", "Excit-14", "Excit-16", "Excit-18", "Excit-19", "Inhib-04", "Inhib-05", "Inhib-09", "Inhib-13"))
seurat_SDH_inclusive <- subset(seurat, idents = c("Excit-01", "Excit-02", "Excit-03", "Excit-04", "Excit-05", "Excit-06", "Excit-08", "Excit-09","Excit-10", "Excit-12", "Excit-13", "Excit-14", "Excit-15", "Excit-16", "Excit-18", "Excit-19", "Inhib-01", "Inhib-02", "Inhib-03", "Inhib-04", "Inhib-05", "Inhib-06", "Inhib-07", "Inhib-09", "Inhib-10", "Inhib-11", "Inhib-12", "Inhib-13"))

seurat_SDH_inhib <- subset(seurat, idents = c("Inhib-01", "Inhib-02", "Inhib-03", "Inhib-04", "Inhib-05", "Inhib-06", "Inhib-07", "Inhib-09", "Inhib-10", "Inhib-11", "Inhib-12", "Inhib-13"))
seurat_SDH_excit <- subset(seurat, idents = c("Excit-01", "Excit-02", "Excit-03", "Excit-04", "Excit-05", "Excit-06", "Excit-08", "Excit-09","Excit-10", "Excit-12", "Excit-13", "Excit-14", "Excit-15", "Excit-16", "Excit-18", "Excit-19"))

seurat_inhib <- subset(seurat, idents = c("Inhib-01", "Inhib-02", "Inhib-03", "Inhib-04", "Inhib-05", "Inhib-06", "Inhib-07", "Inhib-08", "Inhib-09", "Inhib-10", 
                                          "Inhib-11", "Inhib-12", "Inhib-13", "Inhib-14", "Inhib-15", "Inhib-16", "Inhib-17", "Inhib-18", "Inhib-19", "Inhib-20",
                                          "Inhib-22", "Inhib-23", "Inhib-24", "Inhib-25", "Inhib-27"))

#excitatory neurons
seurat_excit <- subset(seurat, idents = c("Excit-01", "Excit-02", "Excit-03", "Excit-04", "Excit-05", "Excit-06", "Excit-08", "Excit-09", "Excit-10", 
                                          "Excit-11", "Excit-12", "Excit-13", "Excit-14", "Excit-15", "Excit-16", "Excit-18", "Excit-19", "Excit-20",
                                          "Excit-21", "Excit-22", "Excit-23", "Excit-24", "Excit-25", "Excit-26", "Excit-29", "Excit-30",
                                          "Excit-31", "Excit-32", "Excit-33", "Excit-34", "Excit-36", "Excit-38"))

runDESeq <- function(seurat, folderName){
  
  #now run a custom DESeq2 analysis. 1. all cells 2. SDH neurons. 3. DDH neurons. 4. SDH + DDH neurons
  counts <- seurat@assays$RNA@counts
  
  # Set up metadata as desired for aggregation and DE analysis
  metadata <- seurat@meta.data
  
  metadata$sample <- as.factor(metadata$sample)
  
  #make sce
  sce_all <- SingleCellExperiment(assays = list(counts = counts), 
                                  colData = metadata)
  
  groups_all <- colData(sce_all)[, c("sample")]
  
  # Aggregate across cluster-sample groups
  pb_all <- aggregate.Matrix(t(counts(sce_all)), 
                             groupings = groups_all, fun = "sum")
  
  pb_real <- t(pb_all)
  
  sids <- purrr::set_names(levels(sce_all$sample))
  
  m <- match(sids, sce_all$sample)
  
  n_cells <- as.numeric(table(sce_all$sample))
  
  ei <- data.frame(colData(sce_all)[m, ], 
                   n_cells, row.names = NULL) %>% 
    select(-"predicted.id")
  ei #sample 13 - 18 have suspiciously high cells 
  
  metadata_sample <- ei
  
  dds_all <- DESeqDataSetFromMatrix(pb_real,
                                    colData = metadata_sample, 
                                    design = ~ study + sex)
  
  
  dds_wald <- DESeq(dds_all, test="Wald")
  
  rldsva <- rlog(dds_wald)
  
  res <- results(dds_wald)
  
  dir <- paste0("DESeq2/", folderName)
  dir.create(dir)
  
  # Create a tibble for LRT results
  res <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
  # Save all results
  write.csv(res,
            paste0(dir, "all_cells_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  # Subset to return genes with padj < 0.05
  sigWald_genes <- res %>% 
    filter(padj < 0.05)
  
  # Save sig results
  write.csv(sigWald_genes,
            paste0(dir, "all_cells_sig_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
}

runDESeq(seurat_SDH_inclusive, "SDH_all/")
runDESeq(seurat_SDH, "SDH_all2/")
runDESeq(seurat_inhib, "Neurons_inhib/")
runDESeq(seurat_excit, "Neurons_excit/")
runDESeq(seurat_SDH_inhib, "Neurons_inhib/")
runDESeq(seurat_SDH_excit, "Neurons_excit/")
