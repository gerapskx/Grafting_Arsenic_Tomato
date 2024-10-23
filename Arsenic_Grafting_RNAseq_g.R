library(readxl)
library(DESeq2)
library(tidyverse)
library(airway)
library(ggplot2)

Arsenic <- read.csv("Arsenic_Grafting_Counts_RNAseq.csv", row.names = 1)

data <- read.csv("Arsenic_Grafting_ColData.csv", row.names = 1)

head(Arsenic)

head(data)

all(colnames(Arsenic) %in% rownames(data))

all(colnames(Arsenic) == rownames(data))

dds <-DESeqDataSetFromMatrix(countData = round(Arsenic),
                       colData = data,
                       design = ~Group)
dds

dds2 <- DESeq(dds)

##############ArsenicGrafted vs UngraftedArsenic#########

U_As_VS_G_As <- results(dds2, contrast = c("Group", "Grafted-As", "Ungrafted-As"))

U_As_VS_G_Asfiltered <- subset(U_As_VS_G_As, padj < 0.05)

U_As_VS_G_Asfiltered <- subset(U_As_VS_G_Asfiltered, log2FoldChange > 2 | log2FoldChange < -2)

summary(U_As_VS_G_Asfiltered)

write.csv(U_As_VS_G_Asfiltered, file="U_As_VS_G_As.csv")

##############UngraftedControl vs GraftedArsenic#########

U_C_VS_G_As <- results(dds2, contrast = c("Group", "Grafted-As", "Ungrafted-Control"))

U_C_VS_G_Asfiltered <- subset(U_C_VS_G_As, padj < 0.05)

U_C_VS_G_Asfiltered <- subset(U_C_VS_G_Asfiltered, log2FoldChange > 2 | log2FoldChange < -2)

summary(U_C_VS_G_Asfiltered)

write.csv(U_C_VS_G_Asfiltered, file="U_C_VS_G_As.csv")

##############UngraftedControl vs GraftedControl#########

U_C_VS_G_c <- results(dds2, contrast = c("Group", "Grafted-Control", "Ungrafted-Control"))

U_C_VS_G_cfiltered <- subset(U_C_VS_G_c, padj < 0.05)

U_C_VS_G_cfiltered <- subset(U_C_VS_G_cfiltered, log2FoldChange > 2 | log2FoldChange < -2)

summary(U_C_VS_G_cfiltered)

write.csv(U_C_VS_G_cfiltered, file="U_C_VS_G_c.csv")

##############UngraftedControl vs UngraftedAs#########

U_C_VS_U_As <- results(dds2, contrast = c("Group", "Ungrafted-As", "Ungrafted-Control"))

U_C_VS_U_Asfiltered <- subset(U_C_VS_U_As, padj < 0.05)

U_C_VS_U_Asfiltered <- subset(U_C_VS_U_Asfiltered, log2FoldChange > 2 | log2FoldChange < -2)

summary(U_C_VS_U_Asfiltered)

write.csv(U_C_VS_U_Asfiltered, file="U_C_VS_U_As.csv")


#################################Plotting#############################################################################


###############################################################pca plot#

vsd <- vst(dds2)
pca_resAs <- prcomp(t(assay(vsd)))
pca_varAs <- pca_resAs$sdev^2 / sum(pca_resAs$sdev^2)
pca_var_percentAs <- pca_varAs * 100
pca_var_percentAs

pcaplot <- plotPCA(vsd, intgroup = "Group", ntop = 500, returnData = TRUE)
pcaplot

selected_groupsAs <- c("Grafted-As", "Grafted-Control","Ungrafted-As", "Ungrafted-Control")
vsdAs_filtered <- vsd[, colData(vsd)$Group %in% selected_groupsAs]
pca_resAs <- prcomp(t(assay(vsdAs_filtered)))
pca_varAs <- pca_resAs$sdev^2 / sum(pca_resAs$sdev^2)
pca_var_percentAs <- pca_varAs * 100
pca_var_percentAs

filtered_pcaplotAs <- pcaplot[pcaplot$Group %in% c("Grafted-As", "Grafted-Control","Ungrafted-As", "Ungrafted-Control"), ]
group_colorsAs <- c("Grafted-As" = "lightcoral", "Grafted-Control" = "darkseagreen3", "Ungrafted-As" = "slategray2", "Ungrafted-Control" = "cadetblue2")
filtered_pcaplotAs

pcaplotenhancedAs <- ggplot(filtered_pcaplotAs, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, ) +
  theme_bw() + scale_fill_manual(values = group_colorsAs) + scale_color_manual(values=group_colorsAs) +
  theme(panel.grid.major = element_line(), panel.grid.minor = element_blank()) +
  labs(x = "PC1: 62% variance", y = "PC2: 12% variance") +
  geom_point(shape = 1,size = 4,colour = "gray20", stroke=1) +
  ggforce::geom_mark_ellipse(aes(fill = Group), alpha = 0.4, expand = 0.006)
pcaplotenhancedAs

###############################################################volcanos

#UngraftedAs_vs_GraftedAs#####################################################################################################

Uas_Gas <- as.data.frame(U_As_VS_G_As)
Uas_Gas <- drop_na(Uas_Gas)
Uas_Gas$GeneID <- row.names(Uas_Gas)

Uas_Gaslabel <- c("Solyc11g020290.3", 
                     "Solyc05g041780.1", 
                     "Solyc03g026270.3",
                     "Solyc07g052180.1",
                     "Solyc03g121900.1",
                     "Solyc02g071000.1")

U_As_VS_G_As_volcano <- ggplot(Uas_Gas, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(
      aes(color = ifelse(abs(log2FoldChange) > 2 & -log10(pvalue) > 1.63, 
                         ifelse(log2FoldChange < 0, "lightcoral", "skyblue"), 
                         "grey")),
      size = 3
  ) + geom_point(shape = 1,size = 3.5,colour = "gray20") +
  scale_color_manual(
    values = c("grey", "skyblue", "lightcoral"),
    labels = c("Not significant", "Down-regulated (2578)", "Up-regulated (2073)")
  ) +
  geom_hline(yintercept = -log10(0.023865501), linetype = "dashed", color = "black") +  # Threshold line
  geom_vline(xintercept = 2, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "black") +
  labs(
    x = "Log2FoldChange",
    y = "-log10(p-value)",
    color = "Regulation"
  ) +
  theme_minimal() + geom_label_repel(data = Uas_Gas[Uas_Gas$GeneID %in% Uas_Gaslabel,],
                                    aes(label = GeneID), size = 3, nudge_x = 0.1)

U_As_VS_G_As_volcano


#UngraftedAs_vs_UngraftedControl#####################################################################################################


Uas_Uc <- as.data.frame(U_C_VS_U_As)
Uas_Uc <- drop_na(Uas_Uc)
Uas_Uc$GeneID <- row.names(Uas_Uc)

Uas_Uclabel <- c("Solyc01g101195.1", 
                  "Solyc12g011350.3", 
                  "Solyc04g077990.3",
                  "Solyc04g064880.4",
                  "Solyc05g041780.1",
                  "Solyc11g020290.3")

U_C_VS_U_As_volcano <- ggplot(Uas_Uc, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(
    aes(color = ifelse(abs(log2FoldChange) > 2 & -log10(pvalue) > 1.89, 
                       ifelse(log2FoldChange < 0, "lightcoral", "skyblue"), 
                       "grey")),
    size = 3
  ) + geom_point(shape = 1,size = 3.5,colour = "gray20") +
  scale_color_manual(
    values = c("grey", "skyblue", "lightcoral"),
    labels = c("Not significant", "Down-regulated (822)", "Up-regulated (255)")
  ) +
  geom_hline(yintercept = -log10(0.012916501), linetype = "dashed", color = "black") +  # Threshold line
  geom_vline(xintercept = 2, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "black") +
  labs(
    x = "Log2FoldChange",
    y = "-log10(p-value)",
    color = "Regulation"
  ) +
  theme_minimal() + geom_label_repel(data = Uas_Uc[Uas_Uc$GeneID %in% Uas_Uclabel,],
                                     aes(label = GeneID), size = 3, nudge_x = 0.1)

U_C_VS_U_As_volcano

#Grafted Control vs Ungrafted Control#####################################################################################################

UG_Uc <- as.data.frame(U_C_VS_G_c)
UG_Uc <- drop_na(UG_Uc)
UG_Uc$GeneID <- row.names(UG_Uc)

UG_Uclabel <- c("Solyc06g051680.1", 
                 "Solyc04g078880.3", 
                 "Solyc03g124110.2",
                 "Solyc09g008760.1",
                 "Solyc06g060640.1",
                 "Solyc02g071000.1")

U_C_VS__c_volcano <- ggplot(UG_Uc, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(
    aes(color = ifelse(abs(log2FoldChange) > 2 & -log10(pvalue) > 1.76, 
                       ifelse(log2FoldChange < 0, "lightcoral", "skyblue"), 
                       "grey")),
    size = 3
  ) + geom_point(shape = 1,size = 3.5,colour = "gray20") +
  scale_color_manual(
    values = c("grey", "skyblue", "lightcoral"),
    labels = c("Not significant", "Down-regulated (2067)", "Up-regulated (836)")
  ) +
  geom_hline(yintercept = -log10(0.017479639), linetype = "dashed", color = "black") +  # Threshold line
  geom_vline(xintercept = 2, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "black") +
  labs(
    x = "Log2FoldChange",
    y = "-log10(p-value)",
    color = "Regulation"
  ) +
  theme_minimal() + geom_label_repel(data = UG_Uc[UG_Uc$GeneID %in% UG_Uclabel,],
                                     aes(label = GeneID), size = 3, nudge_x = 0.1)

U_C_VS_G_c_volcano


#Grafted Arsenic vs Ungrafted Control ######################################################################################################

GA_Uc <- as.data.frame(U_C_VS_G_As)
GA_Uc <- drop_na(GA_Uc)
GA_Uc$GeneID <- row.names(GA_Uc)

GA_Uclabel <- c("Solyc06g051680.1", 
                "Solyc03g026270.3", 
                "Solyc06g074420.1",
                "Solyc08g067940.3",
                "Solyc09g062970.1",
                "Solyc02g071000.1")

U_C_VS_G_As_volcano <- ggplot(GA_Uc, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(
    aes(color = ifelse(abs(log2FoldChange) > 2 & -log10(pvalue) > 1.76, 
                       ifelse(log2FoldChange < 0, "lightcoral", "skyblue"), 
                       "grey")),
    size = 3
  ) + geom_point(shape = 1,size = 3.5,colour = "gray20") +
  scale_color_manual(
    values = c("grey", "skyblue", "lightcoral"),
    labels = c("Not significant", "Down-regulated (1735)", "Up-regulated (1107)")
  ) +
  geom_hline(yintercept = -log10(0.017479639), linetype = "dashed", color = "black") +  # Threshold line
  geom_vline(xintercept = 2, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "black") +
  labs(
    x = "Log2FoldChange",
    y = "-log10(p-value)",
    color = "Regulation"
  ) +
  theme_minimal() + geom_label_repel(data = GA_Uc[GA_Uc$GeneID %in% GA_Uclabel,],
                                     aes(label = GeneID), size = 3, nudge_x = 0.1)

U_C_VS_G_As_volcano


#heatmap######################


vsd <- varianceStabilizingTransformation(dds)
vsd

wpn_vsd <- getVarianceStabilizedData(dds2)

rv_wpn <- rowVars(wpn_vsd)

summary(rv_wpn)

q50_wpn <- quantile( rowVars(wpn_vsd), .50)
q50_wpn
q75_wpn <- quantile( rowVars(wpn_vsd), .75)
q75_wpn
q90_wpn <- quantile( rowVars(wpn_vsd), .90)
q90_wpn

q95_wpn <- quantile( rowVars(wpn_vsd), .95)
q95_wpn

q995_wpn <- quantile( rowVars(wpn_vsd), .995)
q995_wpn

q999_wpn <- quantile( rowVars(wpn_vsd), .999)
q999_wpn
wpn

expr_normalized <- wpn_vsd[ rv_wpn > q999_wpn, ]

expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    x = "Sample",
    y = "normalized expression"
  )


write.csv(expr_normalized, "Heatmap_995_normalized.csv")
input_mat = t(expr_normalized)
input_mat[1:5,1:10]


rownames(input_mat)

data

annotation_col <- as.data.frame(data)
annotation_col

write.csv(wpn_vsd, file = "tomergenormalized.csv")
wpn_vsd

mergedwpn_vsd <-read.csv("tomergenormalized.csv", row.names = 1)

mergedMATRIX <- as.matrix(mergedwpn_vsd)

rv_wpnmerged <- rowVars(mergedMATRIX)

summary(rv_wpnmerged)

q999_merged <- quantile( rowVars(mergedMATRIX), .9985)

q999_merged

expr_normalizedmerged <- mergedMATRIX[ rv_wpnmerged > q999_merged, ]

write.csv(expr_normalizedmerged, file = "heatmap.csv")
data_merged <- data.frame(row.names = c("GAs", "GC", "UAs", "UC"),
  Group = c("Grafted-As", "Grafted-Control", "Ungrafted-As", "Ungrafted-Control")
)



heatmap <- pheatmap(expr_normalizedmerged,
         cluster_rows = TRUE,    # Cluster rows
         cluster_cols = TRUE,    # Cluster columns
         scale = "row",
         # Scale each row (gene) to have mean = 0 and standard deviation = 1
         show_rownames = TRUE,   # Display row names (gene names)
         show_colnames = TRUE,   # Display column names (sample names)
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
          fontsize_col = 12, fontsize_row = 10, 
         clustering_method = "ward.D",
         annotation_col = data_merged)

#fpkm

library(GenomicFeatures)
gtf_file <- "ITAG4.1_gene_models.gff"
txdb <- makeTxDbFromGFF(gtf_file, format="gff")
txdb
genes(txdb)


exons_by_gene <- exonsBy(txdb, by="gene")
exons_by_gene
print(length(exons_by_gene))
print(head(exons_by_gene))

first_gene_exons <- exons_by_gene[[1]]
print(first_gene_exons)
reduced_first_gene_exons <- reduce(first_gene_exons)

gene_lengths <- sum(width(reduce(exons_by_gene)))

feature_lengths <- gene_lengths[rownames(dds)]
mcols(dds)$basepairs / 1000 <- feature_lengths
head(mcols(dds)$basepairs)
feature_lengths <- mcols(dds)$basepairs / 1000
feature_lengths

fpkm_values <- fpkm(dds)
head(fpkm_values)

write.csv(fpkm_values, file="ArsenicFPKM.csv")

#tpm

counts <- counts(dds)
feature_lengths <- mcols(dds)$basepairs / 1000
feature_lengths

na_or_zero_lengths <- is.na(feature_lengths) | feature_lengths == 0
feature_lengths[na_or_zero_lengths] <- 1e-6
feature_lengths
rpk <- counts / feature_lengths
rpk
rpk_sum <- colSums(rpk, na.rm = TRUE)
rpk_sum
tpm <- sweep(rpk, 2, rpk_sum, "/") * 1e6

tpm

log_tpm <- log10(tpm + 1)
log_tpm
write.csv(tpm, file="ArsenicTPM.csv")
