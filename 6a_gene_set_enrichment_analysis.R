#BiocManager::install("fgsea")
library(fgsea)
library(dplyr)
library(stats)
library(ggplot2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data <- read.csv("../R/DESeq2/SDH_all/all_SDH_cells_sig_genes_nopseudo.csv")
GO_file <- "../R/GSEA_data/msigdb_v2023.1.Mm_files_to_download_locally/msigdb_v2023.1.Mm_GMTs/m5.go.cc.v2023.1.Mm.symbols.gmt"

myGO = fgsea::gmtPathways(GO_file)

gene_list <- data$log2FoldChange
names(gene_list) <- data$ï..gene
head(gene_list)
tail(gene_list)
gene_list <- gene_list[!is.na(gene_list)]

# consider removing outliers, Eif2s3y, Uty, Kdm5d, ddx3y, Xist, tsix

if ( any( duplicated(names(gene_list)) )  ) {
  warning("Duplicates in gene names")
  gene_list = gene_list[!duplicated(names(gene_list))]
}
if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
  warning("Gene list not sorted")
  gene_list = sort(gene_list, decreasing = TRUE)
}

head(gene_list)
tail(gene_list)
gene_list <- gene_list[!(names(gene_list) %in% c("Eif2s3y", "Uty", "Kdm5d", "Ddx3y", "Xist", "Tsix")) ]

#barplot(sort(ranks, decreasing= T))
fgRes <- fgsea::fgseaMultilevel(pathways = myGO,
                      stats = gene_list,
                      minSize=15, ## minimum gene set size
                      maxSize=400 ## maximum gene set size
                      ) %>%
  as.data.frame() %>% 
  dplyr::filter(padj < !!0.05) %>% 
  arrange(desc(NES))

message(paste("Number of signficant gene sets =", nrow(fgRes)))

message("Collapsing Pathways -----")
concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                    pathways = myGO,
                                    stats = gene_list)
fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
message(paste("Number of gene sets after collapsing =", nrow(fgRes)))

fgRes$Enrichment = ifelse(fgRes$NES > 0, "Males", "Females")
filtRes = rbind(head(fgRes, n = 10),
                tail(fgRes, n = 10 ))

total_up = sum(fgRes$Enrichment == "Males")
total_down = sum(fgRes$Enrichment == "Females")
header = paste0("Top 10 (Total pathways: Male=", total_up,", Female=",    total_down, ")")

colos = setNames(c("#47A3FC", "#FF69B4"),
                 c("Males", "Females"))

g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
  geom_point( aes(fill = Enrichment, size = size), shape=21) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, color="black"),
    axis.text.y = element_text(size = 12, color="black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size=14, face="bold")
  ) +
  scale_fill_manual(values = colos ) +
  scale_size_continuous(range = c(2,10)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_y_continuous(limits = c(-3.5, 3.5)) +
  #xlim(-3, 3) +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=header)
  

output = list("Results" = fgRes, "Plot" = g1)
output

#sink("GSEA_data/results/sdhneurons_go_cc.txt")
#print(output)
#sink()

