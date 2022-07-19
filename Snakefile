configfile: "config.yaml"

# This rule generates plots comparing the expression of gene markers in both samples. If markers correborates the presence of the cell types of interest in both species, it is possible to do the integration of the datasets.
rule comparing_species:
    input:
        expand("{file}", file = config["rds_file_spec1"]),
        expand("{file}", file = config["rds_file_spec2"]),
        expand("{file}", file = config["known_markers_for_vis"]),
        expand("{file}", file = config["1to1_orthologs"])
    params:
        out_dir1 = "markers_expression",
    output:
    shell:
        """
        mkdir -p {params.out_dir1}

        Rscript comparing_markers_among_samples.R --outdir1 {params.out_dir1} {input[0]} {input[1]} {input[2]} {input[3]}
        """

## This rule generates the raw counts from a dataset contained in a rds file. For example, after clutering and selecting the clusters of interest in Asc-Seurat, one could save the rds file, than generate the raw counts that are necessary for the integration.
rule generates_raw_counts_from_rds:
    input:
        expand("{file}", file = config["rds_file_spec1"]),
        expand("{file}", file = config["rds_file_spec2"]),
    params:
        out_dir1 = config["raw_couns_location"]["species1"],
        out_dir2 = config["raw_couns_location"]["species2"],
    output:
        expand("{test}{file}", test=config["raw_couns_location"]["species1"],
        file = config["drop_util_outs"]),
        expand("{test}{file}", test=config["raw_couns_location"]["species2"],
        file = config["drop_util_outs"]),
    shell:
        """
        rm -r data/
        mkdir -p data/

        Rscript generating_raw_counts_from_rds.R --outdir1 {params.out_dir1} --outdir2 {params.out_dir2} {input[0]} {input[1]}
        """

## To integrate the samples, both datasets must have a set of genes in common. Since this pipeline uses different species, it uses orthology information to identify the orthologs between the species. Then, it is necessary to rename the genes in one of the datasets before performing the integration. Here, we rename in the species 2, since species one is a model with more information available.

rule converting_gene_IDs:
    input:
        expand("{file}", file = config["1to1_orthologs"]),
        expand("{test}{file}", test=config["raw_couns_location"]["species2"],
        file = config["drop_util_outs"])[1]
    params:
        out_dir2 = config["raw_couns_location"]["species2"],
    output:
        expand("data/species2/{file}", file = config["inp_for_seurat"])
    shell:
        """
        Rscript converting_poplar_ids_to_arabidopsis.R --out {params.out_dir2}features.tsv {input[0]} {input[1]}

        rm {params.out_dir2}genes.tsv
        gzip {params.out_dir2}*

        """

### -- In Asc-Seurat, perform the integration of the two datasets (as described in supp. figure X). Add the name of the integrated rds in the config.yaml

# Trajectory inference - This will separate the integrated dataset by species and select the cluster of the lineage of interest
## This steps requires the installation of dynverse. Also, Docker must be running.
rule trajectory_inference_phloem:
    input:
        expand("{file}", file = config["integrated_loc"]) # Generated using Asc-Seurat
    params:
        gc = config["groups_of_clusters"]["g1"], # Phloem
        g_name = config["groups_of_clusters"]["g1_name"],
        start = config["groups_of_clusters"]["start_cluster_g1"],
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        spc_name2 = config["groups_of_clusters"]["spc_name2"]
    output:
        # clusters rds
        expand("{spc_name1}_{g_name}_TI.rds",
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        g_name = config["groups_of_clusters"]["g1_name"]),
        expand("{spc_name2}_{g_name}_TI.rds",
        spc_name2 = config["groups_of_clusters"]["spc_name2"],
        g_name = config["groups_of_clusters"]["g1_name"]),
        # TI rds
        expand("slingshot_TI_{spc_name1}_{g_name}.rds",
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        g_name = config["groups_of_clusters"]["g1_name"]),

        expand("slingshot_TI_{spc_name2}_{g_name}.rds",
        spc_name2 = config["groups_of_clusters"]["spc_name2"],
        g_name = config["groups_of_clusters"]["g1_name"]),
        # Svgs of the trajectory
        expand("{TIplots}_{spc_name1}_{g_name}.svg",
        TIplots = config["TIplots"],
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        g_name = config["groups_of_clusters"]["g1_name"]),
        expand("{TIplots}_{spc_name2}_{g_name}.svg",
        TIplots = config["TIplots"],
        spc_name2 = config["groups_of_clusters"]["spc_name2"],
        g_name = config["groups_of_clusters"]["g1_name"])
    shell:
        """
        Rscript trajectory_inference_analysis_and_comparison.R --group_of_clusters {params.gc} --group_name {params.g_name} --starting_cluster {params.start} --species_name1 "Arabidopsis" --species_name2 "Poplar" {input[0]}
        """

rule trajectory_inference_xylem:
    input:
        expand("{file}", file = config["integrated_loc"]) # Generated using Asc-Seurat
    params:
        gc = config["groups_of_clusters"]["g2"], # Phloem
        g_name = config["groups_of_clusters"]["g2_name"],
        start = config["groups_of_clusters"]["start_cluster_g2"],
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        spc_name2 = config["groups_of_clusters"]["spc_name2"]
    output:
        # clusters rds
        expand("{spc_name1}_{g_name}_TI.rds",
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        g_name = config["groups_of_clusters"]["g2_name"]),
        expand("{spc_name2}_{g_name}_TI.rds",
        spc_name2 = config["groups_of_clusters"]["spc_name2"],
        g_name = config["groups_of_clusters"]["g2_name"]),
        # TI rds
        expand("slingshot_TI_{spc_name1}_{g_name}.rds",
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        g_name = config["groups_of_clusters"]["g2_name"]),

        expand("slingshot_TI_{spc_name2}_{g_name}.rds",
        spc_name2 = config["groups_of_clusters"]["spc_name2"],
        g_name = config["groups_of_clusters"]["g2_name"]),
        # Svgs of the trajectory
        expand("{TIplots}_{spc_name1}_{g_name}.svg",
        TIplots = config["TIplots"],
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        g_name = config["groups_of_clusters"]["g2_name"]),
        expand("{TIplots}_{spc_name2}_{g_name}.svg",
        TIplots = config["TIplots"],
        spc_name2 = config["groups_of_clusters"]["spc_name2"],
        g_name = config["groups_of_clusters"]["g2_name"])
    shell:
        """
        Rscript trajectory_inference_analysis_and_comparison.R --group_of_clusters {params.gc} --group_name {params.g_name} --starting_cluster {params.start} --species_name1 "Arabidopsis" --species_name2 "Poplar" {input[0]}
        """
## Use -c 4 to allow the usage of four cores.
rule tradeseq_phloem:
    input:
        expand("{file}", file = config["integrated_loc"]) # Generated using Asc-Seurat
    params:
        gc = config["groups_of_clusters"]["g1"], # Phloem
        g_name = config["groups_of_clusters"]["g1_name"],
        start = config["groups_of_clusters"]["start_cluster_g1"],
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        spc_name2 = config["groups_of_clusters"]["spc_name2"]
    output:
        expand("TRADEseq_{g_name}_{spc_name1}.rds",
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        g_name = config["groups_of_clusters"]["g1_name"]),
        expand("TRADEseq_{g_name}_{spc_name1}.csv",
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        g_name = config["groups_of_clusters"]["g1_name"]),

        expand("TRADEseq_{g_name}_{spc_name2}.rds",
        spc_name2 = config["groups_of_clusters"]["spc_name2"],
        g_name = config["groups_of_clusters"]["g1_name"]),
        expand("TRADEseq_{g_name}_{spc_name2}.csv",
        spc_name2 = config["groups_of_clusters"]["spc_name2"],
        g_name = config["groups_of_clusters"]["g1_name"]),
    shell:
        """
        Rscript TRADESEQ.R --group_of_clusters {params.gc} --group_name {params.g_name} --starting_cluster {params.start} --species_name1 "Arabidopsis" --species_name2 "Poplar" {input[0]}
        """

## Use -c to allow the usage of multiple cores.
rule tradeseq_xylem:
    input:
        expand("{file}", file = config["integrated_loc"]) # Generated using Asc-Seurat
    params:
        gc = config["groups_of_clusters"]["g2"], # Phloem
        g_name = config["groups_of_clusters"]["g2_name"],
        start = config["groups_of_clusters"]["start_cluster_g2"],
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        spc_name2 = config["groups_of_clusters"]["spc_name2"]
    output:
        expand("TRADEseq_{g_name}_{spc_name1}.rds",
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        g_name = config["groups_of_clusters"]["g2_name"]),
        expand("TRADEseq_{g_name}_{spc_name1}.csv",
        spc_name1 = config["groups_of_clusters"]["spc_name1"],
        g_name = config["groups_of_clusters"]["g2_name"]),

        expand("TRADEseq_{g_name}_{spc_name2}.rds",
        spc_name2 = config["groups_of_clusters"]["spc_name2"],
        g_name = config["groups_of_clusters"]["g2_name"]),
        expand("TRADEseq_{g_name}_{spc_name2}.csv",
        spc_name2 = config["groups_of_clusters"]["spc_name2"],
        g_name = config["groups_of_clusters"]["g2_name"]),
    shell:
        """
        Rscript TRADESEQ.R --group_of_clusters {params.gc} --group_name {params.g_name} --starting_cluster {params.start} --species_name1 "Arabidopsis" --species_name2 "Poplar" {input[0]}
        """

## Rename poplar genes that were not in the list of orthologs 1-to-1 using orthologs listed on Phytozome.

rule use_phytozome_orthologs:
    input:
        "Mappings_Populus_1to1_Arabidopsis_oct_28.txt",
        "Ptrichocarpa_533_v4.1.annotation_info.csv",
        "TRADEseq_Xylem_Poplar.csv",
        "TRADEseq_Phloem_Poplar.csv"
    output:
        "TRADEseq_poplar_xylem_using_phytozome.csv",
        "TRADEseq_poplar_phloem_using_phytozome.csv"
    shell:
        """
        Rscript renaming_using_phytozome_ids.R
        """

## This rule shows the heatmaps (grouping cells in bins)
rule generating_heatmaps:
    input:
        "TRADEseq_Phloem_Arabidopsis.csv",
        "Arabidopsis_Phloem_TI.rds",
        "TRADEseq_Phloem_Poplar.csv",
        "Poplar_Phloem_TI.rds",
        "Ptrichocarpa_533_v4.1.annotation_info.csv",
        "TRADEseq_Xylem_Arabidopsis.csv",
        "Arabidopsis_Xylem_TI.rds",
        "TRADEseq_Xylem_Poplar.csv",
        "Poplar_Xylem_TI.rds"
    output:
        "30_BINs_list_of_genes_positive_correlation_phloem.csv",
        "30_BINs_list_of_genes_positive_correlation_xylem.csv",
        "30_BINs_list_of_genes_UNcorrelation_phloem.csv",
        "30_BINs_list_of_genes_UNcorrelation_xylem.csv",
        "30_BINs_Phloem_all_correlations.csv",
        "30_BINs_Phloem_correlated_genes_heatmaps.svg",
        "30_BINs_Phloem_UNcorrelated_genes_heatmaps.svg",
        "30_BINs_Xylem_all_correlations.csv",
        "30_BINs_Xylem_correlated_genes_heatmaps.svg",
        "30_BINs_Xylem_UNcorrelated_genes_heatmaps.svg",
        "poplar_to_arabi_using_phytozome_list.csv",
        "tradeseq_sig_uniquely_in_poplar_Phloem_COUNTS.csv",
        "tradeseq_sig_uniquely_in_poplar_Phloem.csv",
        "tradeseq_sig_uniquely_in_poplar_Xylem_COUNTS.csv",
        "tradeseq_sig_uniquely_in_poplar_Xylem.csv"
    shell:
        """
        Rscript generating_heatmaps.R
        """

## This rule shows the heatmaps (without grouping cells in bins)
rule generate_heatmap_without_bins:
    input:
        "Arabidopsis_Phloem_TI.rds",
        "Arabidopsis_Xylem_TI.rds",
        "Poplar_Phloem_TI.rds",
        "Poplar_Xylem_TI.rds",
        "TRADEseq_Phloem_Arabidopsis.csv",
        "TRADEseq_Phloem_Poplar.csv",
        "TRADEseq_Xylem_Arabidopsis.csv",
        "TRADEseq_Xylem_Poplar.csv"
    output:
        "All_cells_phloem_all_genes_sig_in_both_datasets.svg",
        "All_cells_phloem_genes_sig_ONLY_IN_POPLAR.svg",
        "All_cells_xylem_all_genes_sig_in_both_datasets.svg",
        "All_cells_xylem_genes_sig_ONLY_IN_POPLAR.svg"
    shell:
        """
        Rscript generate_heatmap_without_bins2.R
        """
