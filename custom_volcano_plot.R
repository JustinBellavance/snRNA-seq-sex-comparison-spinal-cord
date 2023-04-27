#custom volcano plot
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data <- read.csv("data/Neuron_wald_all_genes_nopseudo.csv")

logfc_cols <- "log2FoldChange"
pval_cols <- "padj"

data$log2FoldChange[is.na(data$log2FoldChange)] <- 0
data$pvalue[is.na(data$pvalue)] <- 1
data$padj[is.na(data$padj)] <- 1

data$is_male <- data$log2FoldChange
data$is_male[data$log2FoldChange > 0] <- "TRUE"
data$is_male[data$log2FoldChange < 0] <- "FALSE"

data$is_male[data$padj > 0.05] <- "Not Significant"

data$is_male <- as.factor(data$is_male)

levels(data$is_male)
levels(data$is_male) <- c("Female","Not Significant", "Male")
data$is_male <- factor(data$is_male, levels = c("Female", "Male", "Not Significant"))


# volcano plot in reactive function 
data = data
logfc_col = data$log2FoldChange
pvalue_col = data$padj
gene_col = "gene"
pvalue_thresh = 0.05
logfc_thresh = 0
color_by_de = TRUE
show_logfc_thresh = FALSE
show_pvalue_thresh = FALSE
highlight_genes = NA
x_label = "Log2FC"
y_label = "-log10(padj)"
legend_title = "Sex"
xlim = 5
ylim = 40

# convert pval to -log10(pval)
data$pval <- -log10(data$padj)

# build base of plot
volcano <- ggplot(data, aes(x = log2FoldChange, y = pval))
  
# if show_devec is true color by DE genes
volcano <- volcano +
  geom_point(alpha = .6, aes(colour = is_male)) +
  scale_color_manual(values = c("#FF69B4", "#47A3FC", "#808080" ))
  
volcano <- volcano  + xlim(-2.9, 2.9) +  ylim(0, 45)
  
# add finishing touches to plot
volcanoPlot <- volcano +
  labs(x = x_label, y = y_label, color = legend_title) +
  geom_vline(xintercept=-1,linetype=2) + geom_vline(xintercept=1,linetype=2) +
  theme_classic(base_size = 12)+
  theme(
    axis.text.x = element_text(color="black"),
    axis.text.y = element_text(color="black"),
    axis.ticks = element_line(color = "black"),
    axis.text=element_text(size=14),
    axis.title=element_text(size=14,face="bold")
  )
  
# display plot
volcanoPlot
