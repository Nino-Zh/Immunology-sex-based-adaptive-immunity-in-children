library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)

demographics_file_path <- "C:/Users/ninoz_5pamwj8/OneDrive/Desktop/Immunology/assignment_1_tcell_development/data/metadata/demographics_data.csv"
demographics <- read.csv(demographics_file_path)

rna_seq_path <- "C:/Users/ninoz_5pamwj8/OneDrive/Desktop/Immunology/assignment_1_tcell_development/data/raw/GSE220094_merged_proteincoding_counts.csv.gz"
rna_seq_counts <- read.csv(gzfile(rna_seq_path)) %>%
  rename_with(~str_replace(., "^p006\\.", "HDL042.")) 


# List of target donors (2 males, 2 females)
target_donors <- c("HDL040", "HDL042", "HDL073", "HDL078", "HDL046", "HDL090", "HDL081")


target_donors_all <- c("HDL040", "HDL042", "HDL073", "HDL078", "HDL046", "HDL090",
                       "HDL047", "HDL069", "HDL071", "HDL074", "HDL081")

rna_seq_demographics <- demographics %>%
  filter(Donor.ID %in% target_donors_all)

selected_demographics <- demographics %>%
  filter(Donor.ID %in% target_donors)

# Get all samples from target donors (any tissue)
selected_samples <- colnames(rna_seq_counts)[-1] %>% 
  .[str_extract(., "^HDL[0-9]+") %in% target_donors]


# Create metadata
metadata <- tibble(Sample = selected_samples) %>%
  mutate(
    Donor.ID = str_extract(Sample, "^HDL[0-9]+"),
    Sex = case_when(
      Donor.ID %in% c("HDL040", "HDL073", "HDL090", "HDL081") ~ "Male",
      Donor.ID %in% c("HDL042", "HDL078", "HDL046") ~ "Female"
    )
  ) %>%
  filter(!is.na(Sex))

# Subset RNA Seq data for the selected donor only
count_matrix <- rna_seq_counts %>%
  select(gene_name, all_of(selected_samples)) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()


# DESeq2 analysis for all tissues ------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ Sex
)


dds <- DESeq(dds)
res <- results(dds, contrast = c("Sex", "Female", "Male"))

# Significance calculations
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  mutate(
    padj = ifelse(padj < 1e-10, 1e-10, padj),
    pvalue = ifelse(pvalue < 1e-10, 1e-10, pvalue),
    Sex.Bias = ifelse(log2FoldChange > 0, 'Female', 'Male'),
    log2FoldChange = ifelse(abs(log2FoldChange) > 10, NA, log2FoldChange)
  ) %>%
  filter(!is.na(log2FoldChange))
  

# Get top 10 upregulated genes in females (log2FC > 0)
top_female_genes <- res_df %>%
  filter(Sex.Bias == "Female") %>%
  arrange(desc(log2FoldChange)) %>%
  head(10) %>%
  pull(gene)  # Extract only gene names

# Get top 10 upregulated genes in males (log2FC < 0)
top_male_genes <- res_df %>%
  filter(Sex.Bias == "Male") %>%
  arrange(log2FoldChange) %>%
  head(10) %>%
  pull(gene)  # Extract only gene names

# Save to CSV files without headers
#write_lines(top_female_genes, "top_10_female_genes.csv")
#write_lines(top_male_genes, "top_10_male_genes.csv")



# Volcano Plots for all tissues  -------

# Raw p value
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Sex.Bias), alpha = 0.5, size = 2) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), linetype = 'dashed', color = 'black') + 
  geom_vline(xintercept = c(-5, 5), linetype = 'dashed', color = 'black') +
  scale_color_manual(values = c('Male' = 'green', 'Female' = 'orange')) +
  labs(
    x = 'Log2(Female/Male)',
    y = '-Log10(p-value)',
    color = 'Higher in'
  ) +
  geom_text_repel(
    data = res_df %>% filter(padj < 0.01 & abs(log2FoldChange) > 5),
    aes(label = gene),
    size = 3)+
  theme_minimal()


# Adjusted p value
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Sex.Bias), alpha = 0.3, size = 2) +
  geom_hline(yintercept = -log10(0.01), linetype = 'dashed', color = 'black') + 
  geom_vline(xintercept = c(-5, 5), linetype = 'dashed', color = 'black') +
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  scale_color_manual(values = c('Male' = 'green', 'Female' = 'orange')) +
  labs(
    x = 'Log2(Female/Male)',
    y = '-Log10 (Adjusted p-value)',
    color = 'Higher in'
  ) + 
  geom_text_repel(
    data = res_df %>% filter(padj < 0.01 & abs(log2FoldChange) > 5),
    aes(label = gene),
    size = 3)+
  theme_minimal()


# Heatmap for all tissues ------
# Get top 20 significant genes
top_genes <- res_df %>%
  filter(padj < 0.01) %>%
  arrange(padj) %>%
  head(20) %>%
  pull(gene)

# Convert top_genes to a data frame
top_genes_df <- data.frame(gene = top_genes)

# Write the data frame to a CSV file
#write.csv(top_genes_df, file = "top_genes.csv", row.names = FALSE)

# Normalize counts
vsd <- vst(dds, blind = FALSE)
normalized_counts <- assay(vsd)[top_genes, ]
colnames(normalized_counts) <- gsub(".TEM.*counts", "", colnames(normalized_counts))
# Annotation (ensure Sample is retained)
annotation_col <- metadata %>%
  select(Sample, Sex) %>%
  column_to_rownames("Sample")

sex_colors <- c(Female = "orange", Male = "green")

# Improve heatmap visualization
pheatmap(
  normalized_counts,
  annotation_col   = annotation_col,
  annotation_colors = list(Sex = sex_colors),
  scale           = "row",
  color           = colorRampPalette(c("blue", "white", "red"))(50),
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  fontsize_row    = 8,      # Decrease row label font size
  fontsize_col    = 8,     # Decrease column label font size
  angle_col       = 90,     # Rotate column labels
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,   # Set to FALSE if you want to keep sample order as-is
  cellwidth       = 9,     # Adjust if you need more space per column
  cellheight      = NA      # Adjust if you need more space per row
)




# PCA for under 2 ------
# 1. Remove older children (confirm these are correct donor IDs)
donors_to_remove <- c("HDL046", "HDL081", "HDL090")
filtered_donors <- target_donors_all[!target_donors_all %in% donors_to_remove]

# 2. Create clean metadata (using EXACT column names from your demographics)
metadata <- demographics %>%
  filter(Donor.ID %in% filtered_donors) %>%
  select(Donor.ID, Sex, Race.Ethnicity) %>%
  mutate(
    Race = case_when(
      Race.Ethnicity == "African American" ~ "African American",
      Race.Ethnicity == "White" ~ "White",
      TRUE ~ NA_character_
    )
  ) %>%
  right_join(
    tibble(Sample = selected_samples) %>%
      mutate(Donor.ID = str_extract(Sample, "^HDL[0-9]+")),
    by = "Donor.ID"
  ) %>%
  filter(!is.na(Sex), !is.na(Race)) %>%  # Only keep samples with clear sex/race
  distinct(Sample, .keep_all = TRUE)     # Ensure no duplicate samples

# 3. Create UNSUPERVISED DESeq object (design = ~1)
count_matrix_filtered <- rna_seq_counts %>%
  select(gene_name, all_of(metadata$Sample)) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

dds_unsupervised <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = metadata,
  design = ~ 1  # Critical for unsupervised analysis
)

# 4. PCA Analysis
vsd <- vst(dds_unsupervised, blind = TRUE)  # blind=TRUE for unsupervised
pcaData <- plotPCA(vsd, intgroup = c("Sex", "Race"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# 5. Visualization with clear race/sex encoding
ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(
    aes(color = Sex, shape = Race),
    size = 4,
    alpha = 0.8
  ) +
  scale_color_manual(values = c(Male = "#00FF00", Female = "#FFA500")) +
  scale_shape_manual(values = c("African American" = 17, "White" = 15)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Sex (Color) and Race (Shape)\nFiltered Cohort (≤2 years)") +
  theme_minimal() +

  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )






# PCA for CD4 & CD8 -------

# Modified code for cell type-specific analysis
# 1. Add cell type information to metadata
metadata <- metadata %>%
  mutate(
    CellType = case_when(
      str_detect(Sample, "\\.CD4\\.") ~ "CD4",
      str_detect(Sample, "\\.CD8\\.") ~ "CD8",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(CellType))  # Remove samples without clear cell type

# 2. Split into CD4 and CD8 subsets
cd4_meta <- metadata %>% filter(CellType == "CD4")
cd8_meta <- metadata %>% filter(CellType == "CD8")

# Function to perform cell-type specific PCA
perform_celltype_pca <- function(metadata_subset, cell_type) {
  # Subset count matrix
  count_subset <- rna_seq_counts %>%
    select(gene_name, all_of(metadata_subset$Sample)) %>%
    column_to_rownames("gene_name") %>%
    as.matrix()
  
  # Create DESeq object
  dds <- DESeqDataSetFromMatrix(
    countData = count_subset,
    colData = metadata_subset,
    design = ~ 1
  )
  
  # VST transformation
  vsd <- vst(dds, blind = TRUE)
  
  # PCA analysis
  pca_data <- plotPCA(vsd, intgroup = c("Sex", "Race"), returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  
  # Visualization
  ggplot(pca_data, aes(PC1, PC2)) +
    geom_point(
      aes(color = Sex, shape = Race),
      size = 4,
      alpha = 0.8
    ) +
    scale_color_manual(values = c(Male = "#00FF00", Female = "#FFA500")) +
    scale_shape_manual(values = c("African American" = 17, "White" = 15)) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle(paste("PCA:", cell_type, "T Cells\n(Children ≤2 years)")) +
    theme_minimal() 
}

# Perform analysis for both cell types
cd4_plot <- perform_celltype_pca(cd4_meta, "CD4")
cd8_plot <- perform_celltype_pca(cd8_meta, "CD8")

# Display plots
cd4_plot
cd8_plot