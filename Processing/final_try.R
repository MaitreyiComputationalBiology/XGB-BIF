library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

#Genetic data curation

#Query the genetic data (BRCA)
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
#Download it
GDCdownload(query)

#store it in a data variable
data <- GDCprepare(query, directory = "/Users/priyaltripathi/Multicluster Analysis/Preprocess/TCGA/GDCdata")

# Extract clinical data
clinical_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)
GDCdownload(clinical_query)
clinical_data <- GDCprepare(clinical_query)

# Convert to data frame for analysis
expr_data <- as.data.frame(assay(data))
sample_info <- colData(data) %>% as.data.frame()

# Filter tumor and normal samples
sample_info$sample_type <- as.character(sample_info$sample_type)
expr_data <- expr_data[, sample_info$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]
sample_info <- sample_info[sample_info$sample_type %in% c("Primary Tumor", "Solid Tissue Normal"), ]

# Map Ensembl gene IDs to gene symbols
gene_ids <- rownames(expr_data)
gene_ids <- gsub("\\..*", "", rownames(expr_data))  # Remove everything after the dot
gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

expr_data$Gene_Symbol <- gene_symbols
expr_data <- expr_data[!duplicated(expr_data$Gene_Symbol) & !is.na(expr_data$Gene_Symbol), ]

# Set unique gene symbols as row names
rownames(expr_data) <- expr_data$Gene_Symbol

# Remove the extra column
expr_data$Gene_Symbol <- NULL
sum(duplicated(rownames(expr_data)))

# Replace row names with gene symbols and remove NAs
rownames(expr_data) <- gene_symbols
expr_data <- expr_data[!is.na(rownames(expr_data)), ]

head(rownames(expr_data))
head(keys(org.Hs.eg.db, keytype = "ENSEMBL"))


#Cleaning data
sum(is.na(expr_data))  # Total count of NA values
expr_data <- na.omit(expr_data)

brca_genes <- t(expr_data)
brca_genes <- as.data.frame(brca_genes)
brca_genes <- brca_genes[, colSums(brca_genes != 0) > 0]

# Calculate variance for each gene (col)
gene_variance <- apply(brca_genes, 2, var)

# Set a threshold (e.g., remove genes with variance in the bottom 10%)
variance_threshold <- quantile(gene_variance, 0.10)

# Filter genes with variance above the threshold
brca_genes <- brca_genes[, gene_variance > variance_threshold]

# Count number of samples
tumor_count <- sum(sample_info$sample_type == "Primary Tumor")
normal_count <- sum(sample_info$sample_type == "Solid Tissue Normal")

cat("Processed", tumor_count, "tumor samples and", normal_count, "normal samples.\n")

brca_genes$classes <- sample_info$sample_type
brca_genes_log2 <- brca_genes %>%
  mutate(across(where(is.numeric), ~ log2(. + 1)))

table(brca_genes_log2$classes)

# Save the results
write.csv(brca_genes_log2, "BRCA_gene_expression.csv", row.names = TRUE)
write.csv(sample_info, "BRCA_sample_info.csv", row.names = FALSE)
write.csv(clinical_data, "BRCA_clinical_data.csv", row.names = FALSE)

cat("Data extraction and filtering complete!\n")








#-----------------
#Genetic data curation
#Query the genetic data
query <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
#Download it
GDCdownload(query)

#store it in a data variable
data <- GDCprepare(query)

# Extract clinical data
clinical_query <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)
GDCdownload(clinical_query)
clinical_data <- GDCprepare(clinical_query)

# Convert to data frame for analysis
expr_data <- as.data.frame(assay(data))
sample_info <- colData(data) %>% as.data.frame()

# Filter tumor and normal samples
sample_info$sample_type <- as.character(sample_info$sample_type)
expr_data <- expr_data[, sample_info$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]
sample_info <- sample_info[sample_info$sample_type %in% c("Primary Tumor", "Solid Tissue Normal"), ]

# Map Ensembl gene IDs to gene symbols
gene_ids <- rownames(expr_data)
gene_ids <- gsub("\\..*", "", rownames(expr_data))  # Remove everything after the dot
gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

expr_data$Gene_Symbol <- gene_symbols
expr_data <- expr_data[!duplicated(expr_data$Gene_Symbol) & !is.na(expr_data$Gene_Symbol), ]

# Set unique gene symbols as row names
rownames(expr_data) <- expr_data$Gene_Symbol

# Remove the extra column
expr_data$Gene_Symbol <- NULL
sum(duplicated(rownames(expr_data)))

# Replace row names with gene symbols and remove NAs
rownames(expr_data) <- gene_symbols
expr_data <- expr_data[!is.na(rownames(expr_data)), ]

head(rownames(expr_data))
head(keys(org.Hs.eg.db, keytype = "ENSEMBL"))


#Cleaning data
sum(is.na(expr_data))  # Total count of NA values
sum(is.na(expr_data))  # Total count of NA values
expr_data <- na.omit(expr_data)

lung_genes <- t(expr_data)
lung_genes <- as.data.frame(lung_genes)
lung_genes <- lung_genes[, colSums(lung_genes != 0) > 0]

# Calculate variance for each gene (col)
gene_variance <- apply(lung_genes, 2, var)

# Set a threshold (e.g., remove genes with variance in the bottom 10%)
variance_threshold <- quantile(gene_variance, 0.10)

# Filter genes with variance above the threshold
lung_genes <- lung_genes[, gene_variance > variance_threshold]

# Count number of samples
tumor_count <- sum(sample_info$sample_type == "Primary Tumor")
normal_count <- sum(sample_info$sample_type == "Solid Tissue Normal")

cat("Processed", tumor_count, "tumor samples and", normal_count, "normal samples.\n")

lung_genes$classes <- sample_info$sample_type
lung_genes <- lung_genes %>%
  mutate(across(where(is.numeric), ~ log2(. + 1)))
table(lung_genes$classes)


# Save the results

cat("Data extraction and filtering complete!\n")

write.csv(lung_genes, "Lung_gene_expression.csv", row.names = TRUE)


df <- read.csv("/Users/priyaltripathi/SIP(2023-24)/gastric cancer/trial_dataframe.csv")


