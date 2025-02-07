## Code

## Author: Frida Galan Hernandez
## Email: fridagh@lcg.unam.mx
## Date: 2021-09-30

## Description: This code implements the workflow for analyzing data from the paper "Regional Analysis of the Brain
## Transcriptome in Mice Bred for High and Low Methamphetamine Consumption" by Hitzemann R, Iancu OD, et al. (2019),
## to identify differentially expressed genes in mice bred for high and low methamphetamine consumption.

## The code is divided into four main sections:
## 1) Data loading and exploration
## 2) Data cleaning
## 3) Data analysis
## 4) Data visualization.


# Section 1: Data loading and exploration

## Load the libraries
library(recount3)
library(iSEE)
library("edgeR")
library("ggplot2")
library("limma")

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


## CONTROL OF THE DATA QUALITY

rse_gene_SRP193734$assigned_gene_prop <- rse_gene_SRP193734$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP193734$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP193734$assigned_gene_prop)
with(colData(rse_gene_SRP193734), plot(assigned_gene_prop, sra_attribute.selected_line))
with(colData(rse_gene_SRP193734), plot(assigned_gene_prop, sra_attribute.tissue))

## Save the original data just in case
rse_gene_SRP193734_unfiltred <- rse_gene_SRP193734

## Check if there is a difference between the groups
with(colData(rse_gene_SRP193734), tapply(assigned_gene_prop, sra_attribute.selected_line, sra_attribute.tissue, summary))

## Save only the data that pass the quality control
hist(rse_gene_SRP193734$assigned_gene_prop)
table(rse_gene_SRP193734$assigned_gene_prop < 0.3)
### ----------------------
### FALSE
### 47
### ----------------------

rse_gene_SRP193734 <- rse_gene_SRP193734[, rse_gene_SRP193734$assigned_gene_prop > 0.3]

## Calculate the expression levels of the genes ----- using the counts --------

gene_means <- rowMeans(assay(rse_gene_SRP193734, "counts"))
summary(gene_means)

## Delete genes
rse_gene_SRP193734 <- rse_gene_SRP193734[gene_means > 0.1, ]

## Defining the finals dimension
dim(rse_gene_SRP193734)

## Percentage of genes that we keep
round(nrow(rse_gene_SRP193734) / nrow(rse_gene_SRP193734_unfiltred) * 100, 2)


## NORMALIZE MY DATA

dge <- DGEList(
  counts = assay(rse_gene_SRP193734, "counts"),
  genes = rowData(rse_gene_SRP193734)
)

dge <- calcNormFactors(dge)


## Explore the data
ggplot(
  as.data.frame(colData(rse_gene_SRP193734)),
  aes(y = assigned_gene_prop, x = sra_attribute.selected_line, fill = sra_attribute.selected_line)
) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned gene proportion") +
  xlab("Selected line") +
  ggtitle("Proportion of Assigned Reads by Selected Line") +
  scale_fill_manual(values = c("High Drinker" = "red", "Low Drinker" = "blue"))

## In case of error, we can use the following code to fix it
dev.off()

## Implement the statistical model
mod <- model.matrix(~ sra_attribute.selected_line + sra_attribute.tissue,
                    data = colData(rse_gene_SRP193734))
colnames(mod)

## Analysis of the data

vGene <- voom(dge, mod, plot = TRUE)

eb_result <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_result,
  coef = 2,
  number = nrow(rse_gene_SRP193734),
  sort.by = "none"
)
dim(de_results)
head(de_results)

## Differential expression analysis
table(de_results$adj.P.Val < 0.05)

## Visualization of the results
plotMA(eb_result, coef = 2)

## Volcano plot
volcanoplot(eb_result, coef = 2, highlight = 3, names = de_results$gene_name)
#de_results[de_results$gene_name %in% ]...


