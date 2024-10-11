# Data manipulation
library(tidyverse)
library(openxlsx)
library(data.table)

# Differential expression
library(DESeq2)

# Visualization
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

# Functional enrichment
library(biomaRt)
library(clusterProfiler)
library(ReactomePA)

# Organism specific package
library(org.Hs.eg.db)

#SVA
if (!requireNamespace("sva", quietly = TRUE)) {
  install.packages("sva")
}
library(sva)

# Select only relevant columns for further analysis (rows = genes, columns = samples) and place gene names as row names
library(tibble)
counts <- npc_del_counts[, c(1, 7:ncol(npc_del_counts))] %>% 
  column_to_rownames("Geneid")
dim(counts)
head(counts)
names(counts)

# Format sample names - adding names reflecting real names of samples and their replicates
names(counts) <- basename(names(counts))
colnames(counts) <- c("control_line_p2_1", 
                      "control_line_p2_2",
                      "control_line_p2_3",
                      "DS22q008_p1_1",
                      "DS22q008_p1_2",
                      "DS22q008_p2_1",
                      "DS22q008_p2_2",
                      "DS22q008_p2_3",
                      "DS22q008_p3_1",
                      "DS22q008_p3_2",
                      "DS22q036_p1_1",
                      "DS22q008_p3_3",
                      "DS22q009_p1_1",
                      "DS22q009_p2_1",
                      "DS22q009_p3_1",
                      "DS22q009_p3_2",
                      "DS22q011_p1_1",
                      "DS22q011_p2_1",
                      "DS22q011_p2_2",
                      "DS22q011_p3_1",
                      "DS22q011_p3_2",
                      "DS22q036_p1_2",
                      "DS22q035_p1_1",
                      "DS22q035_p1_2",
                      "DS22q035_p1_3",
                      "DS22q035_p2_1",
                      "DS22q035_p2_2",
                      "DS22q035_p2_3",
                      "DS22q035_p3_1",
                      "DS22q008_p1_3",
                      "DS22q009_p1_2",
                      "DS22q009_p1_3",
                      "DS22q036_p1_3",
                      "DS22q009_p2_2",
                      "DS22q009_p2_3",
                      "DS22q009_p3_3",
                      "DS22q009_p4_1",
                      "DS22q009_p4_2",
                      "DS22q009_p4_3",
                      "DS22q011_p1_2",
                      "DS22q011_p1_3",
                      "DS22q011_p2_3",
                      "DS22q011_p3_3",
                      "DS22q036_p2_1",
                      "DS22q035_p3_2",
                      "DS22q035_p3_3",
                      "DS22q037_p1_1",
                      "DS22q037_p1_2",
                      "DS22q037_p1_3",
                      "DS22q036_p2_2",
                      "DS22q037_p2_1",
                      "DS22q037_p2_2",
                      "DS22q037_p2_3",
                      "DS22q010_p1_1",
                      "DS22q010_p1_2",
                      "DS22q010_p1_3",
                      "DS22q036_p2_3",
                      "control_line_p1_1",
                      "control_line_p1_2",
                      "control_line_p1_3")


# Examine again file
dim(counts)  
head(counts)  
names(counts) 

####################### NPC BATCH FINAL #######################


metadata <- pheno_batch

metadata <- data.frame(
  sample_ids = metadata$replicates_ID,
  groups = metadata$group,
  batch = metadata$batch,
  technical_replicates = metadata$replicates_ID,
  biological_replicates = metadata$biological_rep_ID)

rownames(metadata) <- metadata$sample_ids
metadata$groups <- as.factor(metadata$groups)
metadata$batch <- as.factor(metadata$batch)
metadata$sample_ids <- NULL

npc_counts_batch <- as.matrix(counts)
storage.mode(npc_counts_batch) <- 'integer'

# Ensure count data and metadata match in sample order
all_equal <- all(metadata$technical_replicates == names(npc_counts_batch))   


# Batch correction using ComBat-Seq
adjusted_counts <- ComBat_seq(counts = npc_counts_batch, batch = metadata$batch)

length(adjusted_counts)
dim(adjusted_counts)
class(adjusted_counts)


# Differential expression analysis using DESeq2 *sada u count stavljam adjusted_counts koji imaju uveden batch
dds_batch <- DESeqDataSetFromMatrix(countData = adjusted_counts,
                                  colData = metadata,
                                  design = ~groups)  # Adjusting for batch in DESeq2


dds_collapsed <- collapseReplicates(dds_batch, dds_batch$biological_replicates, dds_batch$technical_replicates)
dds_collapsed$runsCollapsed
colData(dds_collapsed)

#Factor levels
metadata$groups <- factor(metadata$groups, levels = c("control", "del"))


# Remove genes with very low counts (filter genes to contain at least 10 reads accross samples)
keep <- rowSums(counts(dds_collapsed)) >= 10
dds_collapsed <- dds_collapsed[keep,]



# Make rlog transformed counts
rld <- rlog(dds_collapsed, blind = T)


# Check PCA plot
plotPCA(rld, intgroup = c("groups"))

# OPTIONAL: Create customized PCA plot
pcaData <- plotPCA(rld, intgroup = c("groups"), returnData = T)
percentVar <- round(100 * attr(pcaData, "percentVar"))

x_range <- range(pcaData$PC1)
y_range <- range(pcaData$PC2)
x_padding <- (x_range[2] - x_range[1]) * 0.1  # 10% padding
y_padding <- (y_range[2] - y_range[1]) * 0.1



###ovo je pca sa nazivima uzoraka
ggplot(pcaData, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 5) +
  geom_label_repel(
    size = 3.5,
    max.overlaps = Inf,  
    box.padding = 0.5,    
    point.padding = 0.5, 
    segment.size = 0.2     ) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_x_continuous(limits = c(x_range[1] - x_padding, x_range[2] + x_padding)) +
  scale_y_continuous(limits = c(y_range[1] - y_padding, y_range[2] + y_padding)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    text = element_text(size = 14),
    legend.text = element_text(size = 14, hjust = 0.5)) +
  scale_color_brewer(palette = "Dark2")


# Make hierarchical clustering heatmap (sample-to-sample distances)
sampleDists <- as.matrix(dist(t(assay(rld))))
pheatmap(sampleDists, col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))


#Run Deseq with local fitting
dds_collapsed <- DESeq(dds_collapsed, fitType = "local", betaPrior = F)
resultsNames(dds_collapsed)

# Evaluate dispersion estimates
plotDispEsts(dds_collapsed)

# Obtain results
res_dds_collapsed_batch <- results(dds_collapsed) 
print(summary(res_dds_collapsed_batch))

#I want to extract all degs noted by deseq, up and down separately
res <- res_dds_collapsed_batch[order(res_dds_collapsed_batch$padj), ] #sortiranje
sig_genes <- subset(res, padj < 0.1) #izvlacenje DEGova po p < 0.1, jer to je treshold koji vidim kad izvucem rezultate
upregulated_genes <- subset(sig_genes, log2FoldChange > 0)  #izdvajanje up gena
#dodajem simbole i entrezID
upregulated_genes$symbol <- mapIds(org.Hs.eg.db,
                                   keys = row.names(upregulated_genes), 
                                   column = "SYMBOL",
                                   keytype = "ENSEMBL",
                                   multiVals = "first")
upregulated_genes$entrez <- mapIds(org.Hs.eg.db,
                                   keys = row.names(upregulated_genes), 
                                   column = "ENTREZID",
                                   keytype = "ENSEMBL",
                                   multiVals = "first")
up_genes <- as.data.frame(upregulated_genes)
num_upregulated <- nrow(upregulated_genes) #vidim da je to isti broj kao kad pokrenem res1


upregulated_genes_df <- rownames_to_column(up_genes, var = "ENSG_ID") #podesavam da su nazivi redova koji sadrze
#ENSG brojeve zapravo budu kolona jer
#inace se ne cuvaju u excellu

output_file <- "C:/Users/ive_s/Documents/UP_911_all_NPCs_batch.xlsx"
write.xlsx(upregulated_genes_df, file = output_file, rowNames = TRUE)

output_file <- "C:/Users/ive_s/Documents/UP_911_all_NPCs_batch.csv"
write.csv(upregulated_genes_df, file = output_file, row.names = TRUE)


#isto radim i za down gene kao i gore za up
downregulated_genes <- subset(sig_genes, log2FoldChange < 0)
#dodajem simbole i entrezID
downregulated_genes$symbol <- mapIds(org.Hs.eg.db,
                                     keys = row.names(downregulated_genes), 
                                     column = "SYMBOL",
                                     keytype = "ENSEMBL",
                                     multiVals = "first")
downregulated_genes$entrez <- mapIds(org.Hs.eg.db,
                                     keys = row.names(downregulated_genes), 
                                     column = "ENTREZID",
                                     keytype = "ENSEMBL",
                                     multiVals = "first")
down_genes <- as.data.frame(downregulated_genes)
num_downregulated <- nrow(downregulated_genes)
down_genes_df <- rownames_to_column(down_genes, var = "ENSG_ID" )

output_file <- "C:/Users/ive_s/Documents/DOWN_670_all_NPCs_batch.xlsx"
write.xlsx(down_genes_df, file = output_file, rowNames = TRUE)

output_file <- "C:/Users/ive_s/Documents/DOWN_670_all_NPCs_batch.csv"
write.csv(upregulated_genes_df, file = output_file, row.names = TRUE)



# Apply shrinkage to generate more accurate log2 fold change estimates
res_shrunken_batch <- lfcShrink(dds_collapsed, coef = "groups_del_vs_control", type = "apeglm")
print(summary(res_shrunken_batch))

# Evaluate the effect of LFC shrinkage with MA plot 
plotMA(res_dds_collapsed_batch, ylim = c(-2, 2), name = "a")
plotMA(res_shrunken_batch, ylim = c(-2, 2))

# Visualize top 100 genes 
n_select_top <- 100 
select <- order(rowMeans(counts(dds_collapsed, normalized = T)),
                decreasing = T)[1:n_select_top] 
topVarGenes <- head(order(rowVars(assay(rld)), 
                          decreasing = T), n_select_top) 


heatmap_genes <- pheatmap(assay(rld)[topVarGenes, ],
                          cluster_rows = T, 
                          show_rownames = F,
                          show_colnames = T,
                          cluster_cols = T)

# Add gene symbols and ENTREZ IDs
res_shrunken_batch$symbol <- mapIds(org.Hs.eg.db,
                                    keys = row.names(res_shrunken_batch), 
                                    column = "SYMBOL",
                                    keytype = "ENSEMBL",
                                    multiVals = "first")
res_shrunken_batch$entrez <- mapIds(org.Hs.eg.db,
                                    keys = row.names(res_shrunken_batch), 
                                    column = "ENTREZID",
                                    keytype = "ENSEMBL",
                                    multiVals = "first")
head(res_shrunken_batch)

# Define P-adjusted significance cutoff and log2 fold change cutoff to filter out significant DEGs
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Create volcano plot 
npc_degs_batch <- EnhancedVolcano(res_shrunken_batch,
                                  subtitle = NULL, caption = NULL, 
                                  lab = res_shrunken_batch$symbol,
                                  pCutoff = padj_cutoff,
                                  FCcutoff = lfc_cutoff,
                                  x = 'log2FoldChange',
                                  y = 'pvalue',
                                  title = "NPCs DEGs batch")
ggsave("npc_degs_batch_Volcano_plot.png", plot = npc_degs_batch, dpi = 300, units = "in", width = 14, height = 8)

# Extract only significant DEGs based on the above defined cutoffs
res_shrunken_sgnf_batch <- res_shrunken_batch %>%
  as.data.frame() %>%
  arrange(padj) %>%
  filter(padj <= padj_cutoff) %>% 
  filter(abs(log2FoldChange) >= lfc_cutoff) %>%
  rownames_to_column("Geneid")

dim(res_shrunken_sgnf_batch)

# How many genes are remained to be significant? How many are up- or down-regulated?
table(res_shrunken_sgnf_batch$log2FoldChange < 0)
table(abs(res_shrunken_batch$log2FoldChange) >= lfc_cutoff, res_shrunken_batch$padj <= padj_cutoff)


# Export results

write.xlsx(res_shrunken_sgnf_batch %>% as.data.frame() %>% rownames_to_column("Geneid"), 
           file = "C:/Users/ive_s/Desktop/TRANSCRIPTOME 22Q11.2DS/DESeq2_results_padj0.05_lfc1_all_NPCs_batch.xlsx")
write.xlsx(res_shrunken_sgnf_batch, "DESeq2_results_padj0.05_lfc1_all_NPCs_batch.xlsx")


#####HEATMAP FOR DEGS#########
degs_npc_batch <- rownames(res_shrunken_batch)[which(res_shrunken_batch$padj < padj_cutoff & abs(res_shrunken_batch$log2FoldChange) > lfc_cutoff)]
dds_degs_batch <- dds_collapsed[degs_npc_batch, ]
dim(dds_degs_batch)
head(dds_degs_batch)
# Apply variance stabilizing transformation
rlog_counts_npc_batch <- rlog(dds_degs_batch, blind = FALSE)
# Extract the transformed count matrix
rlog_mat <- assay(rlog_counts_npc_batch)

heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Plot the heatmap without sample annotations
heatmap_degs <- pheatmap(
  rlog_mat,
  color = heatmap_colors,
  cluster_rows = TRUE,    # Cluster rows to group genes with similar patterns
  cluster_cols = TRUE,    # Cluster columns to group similar samples
  scale = "row",          # Scale expression values by row (genes)
  show_rownames = FALSE,  # Optionally show row names
  show_colnames = TRUE,   # Show column names (sample IDs)
  main = "Significant DEGs"                        
)

ggsave("heatmap_of_sgnf_DEGs_all_NPCs_batch.png", plot = heatmap_degs, dpi = 300, units = "in", width = 14, height = 8)


