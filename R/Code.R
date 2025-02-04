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

## Inspect the attributes of the samples
rse_gene_SRP193734 <- expand_sra_attributes(rse_gene_SRP193734)

colData(rse_gene_SRP193734)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP193734)))
]

## Use iSEE to visualize the data
iSEE::iSEE(rse_gene_SRP193734)


