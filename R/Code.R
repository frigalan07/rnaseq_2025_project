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

## Inspect the atrributes of the samples
rse_gene_SRP193734$sra.sample_attributes

## [1] "selected line;;High Drinker|Sex;;male|source_name;;High Drinker_Nucleus Accumbens|tissue;;Dissected Tissue (Brain) - Nucleus Accumbens"
## [2] "selected line;;High Drinker|Sex;;male|source_name;;High Drinker_Nucleus Accumbens|tissue;;Dissected Tissue (Brain) - Nucleus Accumbens"
## [3] "selected line;;High Drinker|Sex;;male|source_name;;High Drinker_Nucleus Accumbens|tissue;;Dissected Tissue (Brain) - Nucleus Accumbens"
## [4] "selected line;;High Drinker|Sex;;male|source_name;;High Drinker_Nucleus Accumbens|tissue;;Dissected Tissue (Brain) - Nucleus Accumbens"
## [5] "selected line;;High Drinker|Sex;;male|source_name;;High Drinker_Nucleus Accumbens|tissue;;Dissected Tissue (Brain) - Nucleus Accumbens"
## ...

## Get the data in the columns
names(colData(rse_gene_SRP193734))

## SELECTING OUR DATA
## Casting the data of interest (for a better manipulation)
rse_gene_SRP193734$sra_attribute.selected_line <- factor(rse_gene_SRP193734$sra_attribute.selected_line)

## Check more information about our data
levels(rse_gene_SRP193734$sra_attribute.selected_line)
summary(rse_gene_SRP193734$sra_attribute.selected_line)

## As a matter of interest, we will filter the data to only include samples from the selected section of the brain (Nucleus Accumbens)
rse_gene_SRP193734 <- rse_gene_SRP193734[, rse_gene_SRP193734$sra_attribute.tissue == "Dissected Tissue (Brain) - Nucleus Accumbens"]
## Casting of the tissue variable
colData(rse_gene_SRP193734)$sra_attribute.tissue <- factor(colData(rse_gene_SRP193734)$sra_attribute.tissue)
## Visualize our actual data
rse_gene_SRP193734$sra.sample_attributes

## Resume of our variable of interest
summary(as.data.frame(colData(rse_gene_SRP193734)[
  , grepl("^sra_attribute.*(selected_line|tissu)", colnames(colData(rse_gene_SRP193734)))
]))


