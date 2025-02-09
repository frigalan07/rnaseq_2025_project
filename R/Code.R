## Code

################################################################################
## Author: Frida Galan Hernandez
## Email: fridagh@lcg.unam.mx
## Date: 2025-02-09

## Description: This code implements the workflow for analyzing data from the
## paper "Regional Analysis of the Brain Transcriptome in Mice Bred for High
## and Low Methamphetamine Consumption" by Hitzemann R, Iancu OD, et al. (2019),
## to identify differentially expressed genes in mice bred for high and low
## methamphetamine consumption.

################################################################################

# ==============================================================================
# Section 1: Data loading and pre-processing
# ==============================================================================

## Load the libraries
library(recount3)
library(iSEE)
library("edgeR")
library("ggplot2")
library("limma")
library("pheatmap")

## Explore available mouse datasets in recount3
mouse_projects <- available_projects("mouse")

## Get the project of interest (in this case SRP193734 )
project_info <- subset(
  mouse_projects,
  project == "SRP193734" & project_type == "data_sources"
)

## Create a RangedSummarizedExperiment (RSE) object
rse_gene_SRP193734 <- create_rse(project_info)

## Convert raw counts to read counts to normalize for gene length and sequencing depth
assay(rse_gene_SRP193734, "counts") <- compute_read_counts(rse_gene_SRP193734)

## Expand the attributes of the samples
rse_gene_SRP193734 <- expand_sra_attributes(rse_gene_SRP193734)

colData(rse_gene_SRP193734)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP193734)))
]

## Inspect the attributes of the samples
rse_gene_SRP193734$sra.sample_attributes

## Get the data in the columns
names(colData(rse_gene_SRP193734))

## Casting the data of interest (for a better manipulation)
rse_gene_SRP193734$sra_attribute.selected_line <- factor(rse_gene_SRP193734$sra_attribute.selected_line)
rse_gene_SRP193734$sra_attribute.tissue <- factor(rse_gene_SRP193734$sra_attribute.tissue)

## Check more information about our data
summary(rse_gene_SRP193734$sra_attribute.selected_line)
summary(rse_gene_SRP193734$sra_attribute.tissue)

## Resume of our variable of interest
summary(as.data.frame(colData(rse_gene_SRP193734)[
  , grepl("^sra_attribute.*(selected_line|tissue)", colnames(colData(rse_gene_SRP193734)))
]))


# ==============================================================================
# Section 2: Quality control and data cleaning
# ==============================================================================

## Save a copy of the data before the quality control
rse_gene_SRP193734_unfiltred <- rse_gene_SRP193734

## Calculate the proportion of assigned genes
rse_gene_SRP193734$assigned_gene_prop <- rse_gene_SRP193734$recount_qc.gene_fc_count_all.assigned /
  rse_gene_SRP193734$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP193734$assigned_gene_prop)


## Check if there is a difference between the groups
with(colData(rse_gene_SRP193734), aggregate(assigned_gene_prop,
                                            by = list(sra_attribute.selected_line, sra_attribute.tissue),
                                            FUN = summary)
     )

## Save only the data that pass the quality control
hist(rse_gene_SRP193734$assigned_gene_prop)
table(rse_gene_SRP193734$assigned_gene_prop < 0.3)
### ----------------------
### FALSE
### 143
### ----------------------

rse_gene_SRP193734 <- rse_gene_SRP193734[, rse_gene_SRP193734$assigned_gene_prop > 0.3]

## Calculate the expression levels of the genes ----- using the counts --------
gene_means <- rowMeans(assay(rse_gene_SRP193734, "counts"))
summary(gene_means)

## Delete genes with low expression levels (under 0.1)
rse_gene_SRP193734 <- rse_gene_SRP193734[gene_means > 0.1, ]

## Defining the finals dimension
dim(rse_gene_SRP193734)

## Percentage of genes that we keep
round(nrow(rse_gene_SRP193734) / nrow(rse_gene_SRP193734_unfiltred) * 100, 2)


# ==============================================================================
# Section 3: Normalization of the data
# ==============================================================================

## Normalize the data
dge <- DGEList(
  counts = assay(rse_gene_SRP193734, "counts"),
  genes = rowData(rse_gene_SRP193734)
)

dge <- calcNormFactors(dge)

## Visualize the data after normalization
ggplot(as.data.frame(colData(rse_gene_SRP193734)),
       aes(x = sra_attribute.selected_line, y = assigned_gene_prop, fill = sra_attribute.tissue)) +
  geom_boxplot() +
  labs(x = "Mouse selected line", y = "Assigned gene proportion", fill = "Cerebral tissue") +
  theme_minimal()

# ==============================================================================
# Section 4: Differential expression analysis
# ==============================================================================

## Implement the statistical model
mod <- model.matrix(~ sra_attribute.selected_line + sra_attribute.tissue,
                    data = colData(rse_gene_SRP193734))
colnames(mod)

## Visualize the matrix
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = colData(rse_gene_SRP193734),
  designFormula = ~ sra_attribute.selected_line + sra_attribute.tissue,
  textSizeFitted = 4
)

## Aply the voom transformation (to stabilize the variance)
vGene <- voom(dge, mod, plot = TRUE)

## Fit the linear model for each gene and perform the empirical Bayes moderation
eb_result <- eBayes(lmFit(vGene))

## Extract the results, sorting the genes without any particular order
de_results <- topTable(
  eb_result,
  coef = 2,
  number = nrow(rse_gene_SRP193734),
  sort.by = "none"
)

## Inspect the dimensions and the first rows of the results
dim(de_results)
head(de_results)

## Count the genes with adjusted p-value < 0.05 (significant genes)
table(de_results$adj.P.Val < 0.05)

## Visualization of the changes in the differential expression
plotMA(eb_result, coef = 2)

## Volcano plot highlighting the top 3 genes more significant
volcanoplot(eb_result, coef = 2, highlight = 3, names = de_results$gene_name)
de_results[de_results$gene_name %in% c("Gm9855", "Arhgef15", "Gm1829"), ]


# ==============================================================================
# Section 5: Visualization of the results
# ==============================================================================

## Extract the values of the genes of interest
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

## Create a table with the information of the samples
df <- as.data.frame(colData(rse_gene_SRP193734))

# Subset the columns of interest
df <- df[, c("sra_attribute.selected_line", "sra_attribute.tissue")]
colnames(df) <- c("Selected line", "Tissue")

## Create a heatmap with the expression values of the top 50 genes
pheatmap(exprs_heatmap,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = df
         )


# ==============================================================================
# REPRODUCE THIS CODE
# ==============================================================================
options(width = 120)
sessioninfo::session_info()
