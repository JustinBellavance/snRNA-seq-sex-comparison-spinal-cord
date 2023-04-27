library(Seurat)
library(tidyverse)
library(Matrix.utils)
library(dplyr)
library(Matrix)
library(purrr)
library(tibble)
library(SingleCellExperiment)
library(apeglm)
library(DESeq2)
library(DEGreport)
library(magrittr)
library(sva)

#to be ran after updated SeqSeek Pipeline. With updated 'query.rds' object.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

seurat_small <- readRDS("../results/results_small/query.rds")

#----------
#consider clustering cell types together
seurat_small <- readRDS("../results/results_small/query.rds")

seurat_small@meta.data$predicted.id <- as.factor(seurat_small@meta.data$predicted.id)
seurat_small@meta.data$coarse_type <- as.factor(seurat_small@meta.data$predicted.id)

levels(seurat_small@meta.data$coarse_type)
levels(seurat_small@meta.data$coarse_type) <- c("Astrocytes", "Astrocytes", "Neurons", "Doublets", "Endothelial", 
                                                "Ependymal", "Neurons", "Neurons", "Neurons", "Neurons",
                                                "Neurons", "Neurons", "Neurons", "Neurons", "Neurons", 
                                                "Neurons", "Neurons", "Neurons", "Neurons", "Neurons", 
                                                "Neurons", "Neurons", "Neurons", "Neurons", "Neurons", 
                                                "Neurons", "Neurons", "Neurons", "Neurons", "Neurons",
                                                "Neurons", "Neurons", "Neurons", "Neurons", "Neurons", 
                                                "Neurons", "Neurons", "Neurons", "Neurons", "Neurons", 
                                                "Neurons", "Neurons", "Neurons", "Neurons", "Neurons", 
                                                "Neurons", "Neurons", "Neurons", "Neurons", "Neurons", 
                                                "Neurons", "Neurons", "Neurons", "Neurons", "Neurons", 
                                                "Neurons", "Neurons", "Neurons", "Neurons", "Neurons",
                                                "Neurons", "Neurons", "Neurons", "Junk", "Meninges",
                                                "Meninges", "Microglia", "Neurons", "Neurons", "Oligodendrocytes",
                                                "Oligodendrocytes", "Oligodendrocytes", "Oligodendrocytes", "Oligodendrocytes", "Pericytes",
                                                "Peripheral Glia.BC Derivative", "Neurons", "Schwann")

#add to metadata
seurat_small <- SetIdent(seurat_small, value = seurat_small@meta.data$coarse_type)

#36013 cells
seurat_small <- subset(seurat_small, idents = c("Junk", "Doublets"), invert = TRUE)

seurat <- seurat_small
rm(seurat_small)

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
  select(-"coarse_type")
ei #sample 13 - 18 have suspiciously high cells 

metadata_sample <- ei

dds_all <- DESeqDataSetFromMatrix(pb_real,
                                  colData = metadata_sample, 
                                  design = ~ study + sex)

#dat  <- counts(dds_all, normalized = FALSE)
#idx  <- rowMeans(dat) > 1 #try to remove this see how that works..
#dat  <- dat[idx, ]
#mod  <- model.matrix(~ study + group_id, colData(dds_all))
#mod0 <- model.matrix(~   1, colData(dds_all))

#No sva gives more significant results??
#run to find how many significant surrogate variables.
#svseq <- svaseq(dat, mod, mod0)

#found 4
##ddssva <- dds_clean
#ddssva$SV1 <- svseq$sv[,1]
#ddssva$SV2 <- svseq$sv[,2]
#ddssva$SV3 <- svseq$sv[,3]
#ddssva$SV4 <- svseq$sv[,4]
#design(ddssva) <- ~ SV1 + SV2 + SV3 + SV4 + study + sex

dds_wald <- DESeq(dds_all, test="Wald")

rldsva <- rlog(dds_wald)

res <- results(dds_wald)

dir <- "DESeq2/results_small/"

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
sigLRT_genes <- res %>% 
  filter(padj < 0.05)

# Save sig results
write.csv(sigLRT_genes,
          paste0(dir, "all_cells_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
