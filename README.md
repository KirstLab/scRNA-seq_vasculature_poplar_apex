<p align="center">

  <p align="center">
    <h3 align="center"> Scripts to reproduce the scRNA-seq analysis of the vascular development in poplar apex</h3>
  </p>
</p>


<!-- ABOUT THE PROJECT -->
## Introduction

This repository contains the source code necessary to reproduce the results described in the manuscript "Single-cell analysis of the shoot apex vascular system differentiation in Populus."

The clustering of the scRNA-seq data was performed on Asc-Seurat v2.1. The clustered datasets used as input for the analysis can be downloaded from FigShare (10.6084/m9.figshare.20321787).

## Data analysis

Before integrating the datasets, it is possible to show that the two datasets share cell types by comparing the expression of orthologous genes.

One plot will be generated for each row of the csv file if the listed genes are detected in the dataset. 

### Generating plots comparing gene expression of orthologous genes in the datasets of poplar and arabidopsis.

**Inputs:**

* vascular_arabidopsis.rds (available on FigShare)
* vascular_poplar.rds (available on FigShare)
* NEW_list_of_genes_for_the_heatmap_final_vasculature.csv
* Mappings_Populus_1to1_Arabidopsis_oct_28.txt

```sh

snakemake -c1 -p comparing_species

```

### Extracting the raw data from the clustered dataset contained in am RDS files

This rule generates the raw counts from datasets contained in an rds file. For example, after clustering the data in Asc-Seurat, one could save the rds file and then use the code below to generate the raw counts necessary for the integration.

**Inputs:**

* vascular_arabidopsis.rds (available on FigShare)
* vascular_poplar.rds (available on FigShare)

**Outputs:**

* data/species1/
* data/species2/

```sh
snakemake -c1 -p generates_raw_counts_from_rds
```

### Using the orthologs mapping from orthofinder to rename Poplar genes

**Inputs:**

* Mappings_Populus_1to1_Arabidopsis_oct_28.txt
* data/species2/genes.tsv

**Outputs:**

* data/species2/barcodes.tsv.gz
* data/species2/features.tsv.gz
* data/species2/matrix.mtx.gz

```sh
snakemake -c1 -p converting_gene_IDs
```

### Integrating the datasets on Asc-Seurat

On Asc-Seurat, perform the integration using the file "configuration_file_for_integration_analysis.csv" as input and the default parameters.

### analysis of the integrated dataset

After obtaining the integrated dataset as an RDS file, which can also be downloaded from FigShare.

**NOTE:** This step requires the installation of dynverse, as shown at https://dynverse.org/users/1-installation/. In addition, Docker must be installed and running since dynverse requires it.

**Inputs:**

* CLUSTERED_vasculature_int_arab_poplar_oct_28.rds (available on FigShare)

**Outputs:**

* Arabidopsis_Phloem_TI.rds
* Poplar_Phloem_TI.rds
* slingshot_TI_Arabidopsis_Phloem.rds
* slingshot_TI_Poplar_Phloem.rds
* TI_scatter_Arabidopsis_Phloem.svg
* TI_dendrogram_Arabidopsis_Phloem.svg
* TI_scatter_Poplar_Phloem.svg
* TI_dendrogram_Poplar_Phloem.svg

```sh
snakemake -c1 -p trajectory_inference_phloem
```

**Inputs:**

* CLUSTERED_vasculature_int_arab_poplar_oct_28.rds (available on FigShare)

**Outputs:**

* Arabidopsis_Xylem_TI.rds
* Poplar_Xylem_TI.rds
* slingshot_TI_Arabidopsis_Xylem.rds
* slingshot_TI_Poplar_Xylem.rds
* TI_scatter_Arabidopsis_Xylem.svg
* TI_dendrogram_Arabidopsis_Xylem.svg
* TI_scatter_Poplar_Xylem.svg
* TI_dendrogram_Poplar_Xylem.svg

```sh
snakemake -c1 -p trajectory_inference_xylem
```

### Differentially expressed genes within the trajectory

Use TRADEseq (https://www.bioconductor.org/packages/release/bioc/html/tradeSeq.html) to identify differentially expressed genes within the trajectories.

**Inputs:**

* CLUSTERED_vasculature_int_arab_poplar_oct_28.rds (available on FigShare)

**Outputs:**

* TRADEseq_Phloem_Arabidopsis.rds
* TRADEseq_Phloem_Arabidopsis.csv
* TRADEseq_Phloem_Poplar.rds
* TRADEseq_Phloem_Poplar.csv

```sh
snakemake -c 2 -p tradeseq_phloem
```

**Inputs:**

* CLUSTERED_vasculature_int_arab_poplar_oct_28.rds  (available on FigShare)

**Outputs:**

* TRADEseq_Xylem_Arabidopsis.rds
* TRADEseq_Xylem_Arabidopsis.csv
* TRADEseq_Xylem_Poplar.rds
* TRADEseq_Xylem_Poplar.csv

```sh
snakemake -c 2 -p tradeseq_xylem
```

Rename poplar genes that were not in the list of orthologs 1-to-1 using orthologs listed on Phytozome.

**Inputs:**

* Mappings_Populus_1to1_Arabidopsis_oct_28.txt
* Ptrichocarpa_533_v4.1.annotation_info.csv
* TRADEseq_Xylem_Poplar.csv
* TRADEseq_Phloem_Poplar.csv

**Outputs:**

* TRADEseq_poplar_xylem_using_phytozome.csv
* TRADEseq_poplar_phloem_using_phytozome.csv

```sh
snakemake -c 2 -p use_phytozome_orthologs
```

Generates the heatmaps comparing the expression of orthologs among species

**Inputs:**

* TRADEseq_Phloem_Arabidopsis.csv, Arabidopsis_Phloem_TI.rds
* TRADEseq_Phloem_Poplar.csv
* Poplar_Phloem_TI.rds
* Ptrichocarpa_533_v4.1.annotation_info.csv
* TRADEseq_Xylem_Arabidopsis.csv
* Arabidopsis_Xylem_TI.rds
* TRADEseq_Xylem_Poplar.csv
* Poplar_Xylem_TI.rds

**Outputs:**

* 30_BINs_list_of_genes_positive_correlation_phloem.csv
* 30_BINs_list_of_genes_positive_correlation_xylem.csv
* 30_BINs_list_of_genes_UNcorrelation_phloem.csv
* 30_BINs_list_of_genes_UNcorrelation_xylem.csv
* 30_BINs_Phloem_all_correlations.csv
* 30_BINs_Phloem_correlated_genes_heatmaps.svg
* 30_BINs_Phloem_UNcorrelated_genes_heatmaps.svg
* 30_BINs_Xylem_all_correlations.csv
* 30_BINs_Xylem_correlated_genes_heatmaps.svg
* 30_BINs_Xylem_UNcorrelated_genes_heatmaps.svg
* poplar_to_arabi_using_phytozome_list.csv
* tradeseq_sig_uniquely_in_poplar_Phloem_COUNTS.csv
* tradeseq_sig_uniquely_in_poplar_Phloem.csv
* tradeseq_sig_uniquely_in_poplar_Xylem_COUNTS.csv
* tradeseq_sig_uniquely_in_poplar_Xylem.csv

```sh
snakemake -c 2 -p generating_heatmaps
```

**Inputs:**

* Arabidopsis_Phloem_TI.rds
* Arabidopsis_Xylem_TI.rds
* Poplar_Phloem_TI.rds
* Poplar_Xylem_TI.rds
* TRADEseq_Phloem_Arabidopsis.csv
* TRADEseq_Phloem_Poplar.csv
* TRADEseq_Xylem_Arabidopsis.csv
* TRADEseq_Xylem_Poplar.csv

**Outputs:**

* All_cells_phloem_correlated_genes_heatmaps_same_dimention.svg
* All_cells_phloem_correlated_genes_heatmaps_corrected_proportion.svg
* All_cells_phloem_UNcorrelated_genes_heatmaps_same_dimention.svg
* All_cells_phloem_UNcorrelated_genes_heatmaps_corrected_proportion.svg
* All_cells_xylem_correlated_genes_heatmaps_same_dimention.svg
* All_cells_xylem_correlated_genes_heatmaps_corrected_proportion.svg
* All_cells_xylem_UNcorrelated_genes_heatmaps_same_dimention.svg
* All_cells_xylem_UNcorrelated_genes_heatmaps_corrected_proportion.svg

```sh
snakemake -c 2 -p generate_heatmap_without_bins
```