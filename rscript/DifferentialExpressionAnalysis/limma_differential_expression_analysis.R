# Differential Expression Analysis: limma
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")
# BiocManager::install("dplyr")
# BiocManager::install("ggrepel")
# BiocManager::install("pheatmap")

library(limma)
library(dplyr)
setwd("/home/lijiawei/Bioinfor/Bioinfor-Analysis-Private/rscript/DifferentialExpressionAnalysis")

# 1\ Read data.
data_counts <- read.table("counts_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
head(data_counts)

# 2\ Get the label list.
data_list <- c(rep("Treat", 66), rep("CK", 148))
data_list <- factor(data_list, levels = c("CK", "Treat"), ordered = FALSE)
data_list <- model.matrix(~ factor(data_list) + 0)
colnames(data_list) <- c("CK", "Treat")
head(data_list)

# 3\ Different Expression Analysis
fit <- limma::lmFit(data_counts, data_list)
contrasts <- limma::makeContrasts(CK - Treat, levels = data_list)
result <- limma::contrasts.fit(fit, contrasts)
result <- limma::eBayes(result)
result <- limma::topTable(result, n = Inf, adjust = "fdr")
head(result)

# 4\ Output the data.
diff_all <- na.omit(result)
# write.csv(diff_all, "limma_out.csv")

foldchange <- 0.5
padj <- 0.05
diff_up_down <- diff_all[(diff_all$P.Value < padj & (diff_all$logFC > foldchange | diff_all$logFC < (- foldchange))), ]
write.csv(diff_up_down, "diff_up_down.csv")

diff_up <- diff_all[(diff_all$P.Value < padj & (diff_all$logFC > foldchange)), ]
# write.csv(diff_up, "diff_up.csv")

diff_down <- diff_all[(diff_all$P.Value < padj & (diff_all$logFC < -foldchange)), ]
# write.csv(diff_down, "diff_down.csv")

# --------------------- Plot ---------------------

# 1\ Volcano Plot.
library(ggplot2)
library(ggrepel)

pdf("volcano.pdf", width = 7, height = 6.5)
significant <- ifelse((diff_all$P.Value < 0.05 & abs(diff_all$logFC) > 0.5), ifelse(diff_all$logFC > 0.5, "Up", "Down"), "Not")

ggplot(diff_all, aes(logFC, -log10(P.Value))) +
geom_point(aes(col = significant)) +
scale_color_manual(values = c("#0072B6", "grey", "#BC3C28")) +
labs(title = " ") +
geom_vline(xintercept = c(-0.5, 0.5), colour = "black", linetype = "dashed") +
geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed") +
theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +
labs(x = "log2(FoldChange)", y = "-log10(Pvalue)") +
theme(axis.text = element_text(size = 13), axis.title = element_text(size = 13)) +
theme_bw()
dev.off()

# 2\ Heatmap Plot.
library(pheatmap)

DEG_id <- unique(rownames(diff_up_down))
DEG_exp <- data_counts[DEG_id, ]
head(DEG_id)

hmexp <- na.omit(DEG_exp)
annotation_col <- data.frame(Group = factor(c(rep("Treat", 66), rep("CK", 148))))
rownames(annotation_col) <- colnames(hmexp)

pdf(file = "heatmap.pdf", height = 8, width = 12)
p <- pheatmap(hmexp,
              annotation_col = annotation_col,
              color = colorRampPalette(c("blue", "white", "red"))(50),
              cluster_cols = FALSE,
              show_rownames = FALSE,
              show_colnames = FALSE,
              scale = "row", ## none, row, column
              fontsize = 12,
              fontsize_row = 12,
              fontsize_col = 6,
              border = FALSE)
print(p)
dev.off()
