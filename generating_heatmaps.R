require(Seurat)
require(tidyverse)
require(dyno)
require(grDevices)
require(slingshot)
require(RColorBrewer)
require(ComplexHeatmap)
require(grDevices)
require(circlize)
require(SingleCellExperiment)
options(scipen = 999)

set.seed(1407)

## Help function that calculates the average of the expression in each bin, accordingly with the order of the cells in the trajectory.
test_function <- function (tradeseq_values_db1 = tradeseq_values,
                           rds_data_db1 = rds_data,
                           pvalue_value_db1 = pvalue_value,
                           N_GROUPS = N_GROUPS,
                           starting_cluster = starting_cluster,
                           bar_plot_prefix = bar_plot_prefix) {
    
    list_of_genes_db1 <- vroom::vroom(tradeseq_values_db1) %>%
        dplyr::filter(fdr < pvalue_value_db1) %>%
        dplyr::filter(!is.na(genes)) %>%
        dplyr::select(genes)
    
    # RDS_files
    sc_data_traj_inf <- readRDS(rds_data_db1)
    
    # Genes with couns equal 0 in all cells are from the other species that was integrated.
    COUNTS <- as.matrix(sc_data_traj_inf@assays$RNA@counts)
    COUNTS <- COUNTS[rowSums(COUNTS) > 0, ]
    
    sce <- SingleCellExperiment(assays = List(counts = as.matrix(COUNTS) ) )
    reducedDims(sce) <- SimpleList(PCA = sc_data_traj_inf@reductions$pca@cell.embeddings)
    colData(sce)$seurat_clusters <- sc_data_traj_inf$seurat_clusters
    
    sds <- slingshot::slingshot(
        sce,
        clusterLabels = "seurat_clusters",
        start.clus = starting_cluster,
        reducedDim = 'PCA',
        shrink = 1L,
        reweight = TRUE,
        reassign = TRUE,
        maxit = 10L,
        smoother = "smooth.spline",
    )
    
    ## I am using the new version of slingshot, that is a bit different of the one in the app. So, I plot here the trajectory to show that it looks similar to the one in the app.
    
    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    plotcol <- colors[cut(sds$slingPseudotime_1, breaks=100)]
    plot(reducedDims(sds)$PCA, col = plotcol, pch=16, asp = 1)
    lines(SlingshotDataSet(sds), lwd=2, col='black')
    
    pst.ord <- order(sds$slingPseudotime_1, na.last = NA)
    heatdata <- assays(sds)$counts
    heatdata <- heatdata[rownames(heatdata) %in% list_of_genes_db1$genes, pst.ord]
    
    cells_order_in_traject <- data_frame(cells = colnames(heatdata))
    # By this point, I have the cells in the correct order in the columns. Now, we need to generate the bins by deviding the cells in the number of selected clusters.
    
    ## Lognormalize the data
    heatdata <- log1p(heatdata)
    
    colnames_heatmap <- c( 
        rep( 1:(N_GROUPS-1),
             each =  floor( ncol(heatdata) / (N_GROUPS) ) ),
        
        rep(N_GROUPS, ncol(heatdata) - ( (N_GROUPS-1) * floor(ncol(heatdata)/(N_GROUPS) ) ) )
    )
    
    print(table(colnames_heatmap))
    
    cells_per_bin <- table(colnames_heatmap)
    
    svg( paste0(bar_plot_prefix,"n_of_cells_per_bin.svg"),
         width = 12,
         height = 10,
         pointsize = 12)
    bar <-  barplot(cells_per_bin,
                    ylim = c(0, max(cells_per_bin) + 10)) 
    dev.off()
    
    colnames(heatdata) <- colnames_heatmap
    
    ## Average the expression in each group
    averaged_heatdata <- data.frame(row.names = rownames(heatdata))
    for(c in unique(colnames_heatmap) ) {
        
        sub <- heatdata[, colnames(heatdata) == c]
        sub2 <- as.data.frame(rowMeans(sub))
        
        averaged_heatdata <- cbind(averaged_heatdata, sub2)
    }
    
    averaged_heatdata <- as.matrix(averaged_heatdata)
    
    heatclus <- sds$seurat_clusters[pst.ord]
    
    cells_order_in_traject_meta <- sc_data_traj_inf@meta.data
    cells_order_in_traject_meta$cells <- rownames(cells_order_in_traject_meta)
    cells_order_in_traject_meta <- cells_order_in_traject_meta[, c(7, 8)]
    
    cells_order_in_traject_cluster <- plyr::join(cells_order_in_traject,
                                                 cells_order_in_traject_meta,
                                                 by = "cells" )
    
    return(list(bin_average = averaged_heatdata,
                cells_order = cells_order_in_traject_cluster))
    
}

######################################################
##################### Comparison #####################
######################################################
prefix = "30_BINs_"
N_GROUPS = 30
PVALUE_THRES = 0.01
starting_cluster = 8

#####################
## Phloem analysis
#####################
arabi_phloem_bins <- test_function(tradeseq_values_db1 = "TRADEseq_Phloem_Arabidopsis.csv",
                                   rds_data_db1 = "Arabidopsis_Phloem_TI.rds",
                                   pvalue_value_db1 = PVALUE_THRES,
                                   N_GROUPS = N_GROUPS,
                                   starting_cluster = starting_cluster,
                                   bar_plot_prefix = paste0(prefix, "_", "Arabidopsis_Phloem"))

arabi_phloem_bins_avg <- data.frame(arabi_phloem_bins[[1]])
colnames(arabi_phloem_bins_avg)[1:ncol(arabi_phloem_bins_avg)] <- paste0("bin_arab", seq(1, N_GROUPS))
arabi_phloem_bins_avg$features <- rownames(arabi_phloem_bins_avg)

print(paste0(nrow(arabi_phloem_bins_avg), " genes in arabidopsis - Phloem"))

pop_phloem_bins <- test_function(tradeseq_values_db1 = "TRADEseq_Phloem_Poplar.csv",
                                 rds_data_db1 = "Poplar_Phloem_TI.rds",
                                 pvalue_value_db1 = PVALUE_THRES,
                                 N_GROUPS = N_GROUPS,
                                 starting_cluster = starting_cluster,
                                 bar_plot_prefix = paste0(prefix, "_", "Poplar_Phloem"))

pop_phloem_bins_avg <- data.frame(pop_phloem_bins[[1]])
colnames(pop_phloem_bins_avg)[1:ncol(pop_phloem_bins_avg)] <- paste0("bin_pop", seq(1, N_GROUPS))
pop_phloem_bins_avg$features <- rownames(pop_phloem_bins_avg)

print(paste0(nrow(pop_phloem_bins_avg), " genes in poplar - Phloem"))

# Identify common genes in both species
# Reads the Phytozome dataset
orthofinder <-  vroom::vroom("Mappings_Populus_1to1_Arabidopsis_oct_28.txt") %>%
    dplyr::select(populus_trichocarpa, arabidopsis_thaliana) %>%
    dplyr::rename( feature_id = populus_trichocarpa ) %>%
    dplyr::rename(Arabidopsis = arabidopsis_thaliana) %>%
    dplyr::mutate(source = "Orthofinder")

Phyto <- vroom::vroom("Ptrichocarpa_533_v4.1.annotation_info.csv", 
                      col_names = T) %>%
    dplyr::rename( feature_id = Gene ) %>%
    dplyr::mutate( combine = paste0(feature_id, "_", Arabidopsis) ) %>%
    dplyr::filter( !duplicated(combine) ) %>%
    dplyr::filter( !is.na(Arabidopsis) ) %>%
    dplyr::mutate( Arabidopsis2 = make.unique(Arabidopsis) ) %>%
    dplyr::select("feature_id", "Arabidopsis2") %>%
    dplyr::rename(Arabidopsis = Arabidopsis2) %>%
    dplyr::filter( ! feature_id %in% orthofinder$feature_id ) %>%
    dplyr::mutate(source = "Phytozome")

Phyto_ortho_combined <- rbind(orthofinder, Phyto) %>%
    dplyr::mutate(Arabidopsis = make.unique(Arabidopsis))

write.table(Phyto_ortho_combined,
            "poplar_to_arabi_using_phytozome_list.csv",
            col.names = T,
            row.names = F,
            sep = ",",
            quote = F)

Phyto_ortho_combined <- Phyto_ortho_combined %>%
    dplyr::select(- source)

to_rename <- pop_phloem_bins_avg %>%
    dplyr::filter( grepl(features,
                         pattern = "Potri") ) %>%
    dplyr::rename(feature_id = features)

to_keep <- pop_phloem_bins_avg %>%
    dplyr::filter( !grepl(features,
                          pattern = "Potri") ) %>%
    dplyr::rename(feature_id = features)

# This will rename any poplar gene that wasnt in the Orthofinder list, but is in the phytozome list.
to_rename_ids <- as.data.frame(to_rename[, ncol(to_rename)])
colnames(to_rename_ids) <- "feature_id"

## Get the names of the genes
to_rename2 <- merge(to_rename_ids, Phyto_ortho_combined, by = "feature_id",  all.x = T)

## Select only the genes that have ID on arabidopsis
to_rename2 <- to_rename2[!is.na(to_rename2$Arabidopsis), ] %>%
    dplyr::mutate( combine = paste0(feature_id, "_", Arabidopsis) ) %>%
    dplyr::filter(!duplicated(combine))

to_rename2 <- merge(to_rename, to_rename2 , by = "feature_id", all.x = T) %>%
    dplyr::rename(feature_id2 = feature_id) %>%
    dplyr::select(- c( "combine") ) %>%
    dplyr::rename(feature_id = Arabidopsis)

no_match_genes <- to_rename2 %>%
    dplyr::filter(is.na(feature_id)) %>%
    dplyr::select(-feature_id) %>%
    dplyr::rename(features = feature_id2) %>%
    dplyr::relocate(features, .after = last_col())

to_rename2 <- to_rename2 %>%
    dplyr::select(-feature_id2) %>%
    dplyr::filter(!is.na(feature_id)) 

rownames(to_rename2) <- to_rename2$feature_id

pop_phloem_bins_avg2 <- dplyr::bind_rows(to_keep,
                                         to_rename2) %>%
    dplyr::rename(features = feature_id)

rownames(pop_phloem_bins_avg2) <- pop_phloem_bins_avg2$features

## Duplicates the arabidopsis entrances of genes which more than one copy are present in the poplar dataset, so we can join the heatmap.

dup_features <- pop_phloem_bins_avg2 %>%
    dplyr::select(features) %>%
    dplyr::filter(grepl(features, pattern = "\\.")) %>%
    dplyr::mutate(features2 = features) %>%
    tidyr::separate(features, into = c("features", "number") ) %>%
    dplyr::select(features, features2)

dup_features_arab <- merge(dup_features, arabi_phloem_bins_avg,  by = "features") %>%
    dplyr::select(-features) %>%
    dplyr::rename(features = features2)

arabi_phloem_bins_avg2 <- dplyr::bind_rows(arabi_phloem_bins_avg,
                                           dup_features_arab)

arabi_phloem_bins_avg2$features <- make.unique(arabi_phloem_bins_avg2$features)
rownames(arabi_phloem_bins_avg2) <- arabi_phloem_bins_avg2$features

combined_avg_bins <- merge(arabi_phloem_bins_avg2, pop_phloem_bins_avg2, by = "features")

combined_avg_bins_names <- combined_avg_bins%>%
    tidyr::separate(features, into = c("name", "number")) %>%
    dplyr::select("name")

print( paste0("There are ", nrow(combined_avg_bins), " common genes in the phloem sample") )

unique_poplar <- pop_phloem_bins_avg2 %>%
    dplyr::filter(!features %in% combined_avg_bins$features) %>%
    dplyr::mutate(features2 = features) %>%
    tidyr::separate(features2, into = c("name", "number")) %>%
    dplyr::filter(!name %in% combined_avg_bins_names$name) %>%
    dplyr::select(-name, - number)

unique_poplar <- unique(rbind(unique_poplar, no_match_genes))

write.table(unique_poplar,
            "tradeseq_sig_uniquely_in_poplar_Phloem_COUNTS.csv",
            col.names = T,
            row.names = F,
            sep = ",",
            quote = F)

write.table(unique_poplar$features,
            "tradeseq_sig_uniquely_in_poplar_Phloem.csv",
            col.names = T,
            row.names = F,
            sep = ",",
            quote = F)
## Generates the correlation for each gene separately. besides the gene name, the first 30 columns are the Arabidopsis data and the other 30 are the poplar.
gene_cor_df_comp <- data.frame()
for( i in 1:nrow(combined_avg_bins) ) {
    
    gene_cor <- cor.test(as.numeric(combined_avg_bins[i, 2:(N_GROUPS + 1)] ),
                         as.numeric(combined_avg_bins[i, (N_GROUPS + 2):ncol(combined_avg_bins)]),
                         method = "pearson" )
    gene_cor_df <- data.frame(gene = combined_avg_bins[i,1],
                              corr = gene_cor$estimate,
                              p_value = gene_cor$p.value)
    
    gene_cor_df_comp <- rbind(gene_cor_df_comp,
                              gene_cor_df)
    
}

## Multiple test (FDR) correction of the p-values
gene_cor_df_comp$fdr <- format( p.adjust(p = gene_cor_df_comp$p_value,
                                         method = "fdr",
                                         n = length(gene_cor_df_comp$p_value) ),
                                scientific = F)

gene_cor_df_comp_sig <- gene_cor_df_comp %>% 
    dplyr::filter(fdr < PVALUE_THRES) %>%
    dplyr::filter(corr > 0.5)

gene_cor_df_comp_diff <- gene_cor_df_comp[gene_cor_df_comp$fdr >= PVALUE_THRES, ]

write.csv(gene_cor_df_comp, paste0(prefix, "Phloem_all_correlations.csv"), row.names = F)

# Heatmap generation

## @@ Correlated genes

## Selecting the expression values of the genes that are correlated to generate the heatmap
corr_genes_arab <- arabi_phloem_bins_avg2 %>%
    dplyr::filter(features %in% gene_cor_df_comp_sig$gene) %>%
    dplyr::arrange(features) %>%
    dplyr::select(-features)

corr_genes_pop <- pop_phloem_bins_avg2 %>%
    dplyr::filter(features %in% gene_cor_df_comp_sig$gene) %>%
    dplyr::arrange(features) %>%
    dplyr::select(-features)

### @@@ >>>> dynutils::scale_quantile Cut off outer quantiles and rescale to a [0, 1] range
corr_genes_arab_scaled <- t(dynutils::scale_quantile(t(corr_genes_arab)))
corr_genes_pop_scaled <- t(dynutils::scale_quantile(t(corr_genes_pop)))

col_fun = colorRamp2( breaks = c(0, 0.5, 1),
                      c("#ffffcc", "#fd8d3c", "#800026") )

heatmap_corr_arab <- ComplexHeatmap::Heatmap(corr_genes_arab_scaled,
                                             border = TRUE,
                                             #rect_gp = grid::gpar(col = "white", lwd = 0.01),
                                             name = "Expression",
                                             cluster_rows = T,
                                             cluster_columns = F,
                                             col = col_fun,
                                             show_row_names = T,
                                             show_column_names = F,
                                             show_row_dend = F,
                                             use_raster = F,
                                             row_names_gp = grid::gpar(fontsize = 10),
                                             column_title = "Arabidopsis")

heatmap_corr_pop <- ComplexHeatmap::Heatmap(corr_genes_pop_scaled,
                                            border = TRUE,
                                            #rect_gp = grid::gpar(col = "white", lwd = 0.01),
                                            name = "Expression",
                                            cluster_rows = T,
                                            cluster_columns = F,
                                            show_row_names = T,
                                            show_column_names = F,
                                            show_row_dend = F,
                                            use_raster = F,
                                            col = col_fun,
                                            row_names_gp = grid::gpar(fontsize = 10),
                                            column_title = "Populus")

svg( paste0(prefix,"Phloem_correlated_genes_heatmaps.svg"),
     width = 10,
     height = 30,
     pointsize = 12)
draw( heatmap_corr_arab + heatmap_corr_pop )
dev.off()

write.csv(gene_cor_df_comp_sig,
          paste0(prefix,"list_of_genes_positive_correlation_phloem.csv"),
          row.names = F)

write.csv(gene_cor_df_comp_diff,
          paste0(prefix,"list_of_genes_UNcorrelation_phloem.csv"),
          row.names = F)

### Genes that does not correlate
## Selecting the expression values of the genes that are correlated to generate the heatmap
uncorr_genes_arab <- arabi_phloem_bins_avg2 %>%
    dplyr::filter(features %in% gene_cor_df_comp_diff$gene) %>%
    dplyr::arrange(features) %>%
    dplyr::select(-features)

uncorr_genes_pop <- pop_phloem_bins_avg2 %>%
    dplyr::filter(features %in% gene_cor_df_comp_diff$gene)%>%
    dplyr::arrange(features) %>%
    dplyr::select(-features)

uncorr_genes_arab_scaled <- t(scale(t(uncorr_genes_arab)))
uncorr_genes_pop_scaled <- t(scale(t(uncorr_genes_pop)))

heatmap_uncorr_arab <- ComplexHeatmap::Heatmap(uncorr_genes_arab_scaled,
                                               border = TRUE,
                                               #rect_gp = grid::gpar(col = "white", lwd = 0.01),
                                               name = "Expression",
                                               cluster_rows = T,
                                               cluster_columns = F,
                                               show_row_names = T,
                                               show_column_names = F,
                                               show_row_dend = F,
                                               use_raster = F,
                                               col = col_fun,
                                               row_names_gp = grid::gpar(fontsize = 10),
                                               column_title = "Arabidopsis")

heatmap_uncorr_pop <- ComplexHeatmap::Heatmap(uncorr_genes_pop_scaled,
                                              border = TRUE,
                                              #rect_gp = grid::gpar(col = "white", lwd = 0.01),
                                              name = "Expression",
                                              cluster_rows = T,
                                              cluster_columns = F,
                                              show_row_names = T,
                                              show_column_names = F,
                                              show_row_dend = F,
                                              use_raster = F,
                                              col = col_fun,
                                              row_names_gp = grid::gpar(fontsize = 10),
                                              column_title = "Populus")

svg( paste0(prefix,"Phloem_UNcorrelated_genes_heatmaps.svg"), width = 10, height = 30, pointsize = 12)
draw( heatmap_uncorr_arab + heatmap_uncorr_pop )
dev.off()

#####################
## Xylem analysis
#####################

arabi_xylem_bins <- test_function(tradeseq_values_db1 = "TRADEseq_Xylem_Arabidopsis.csv",
                                   rds_data_db1 = "Arabidopsis_Xylem_TI.rds",
                                   pvalue_value_db1 = PVALUE_THRES,
                                   N_GROUPS = N_GROUPS,
                                   starting_cluster = starting_cluster,
                                   bar_plot_prefix = paste0(prefix, "_", "Arabidopsis_Xylem"))

arabi_xylem_bins_avg <- data.frame(arabi_xylem_bins[[1]])
colnames(arabi_xylem_bins_avg)[1:ncol(arabi_xylem_bins_avg)] <- paste0("bin_arab", seq(1, N_GROUPS))
arabi_xylem_bins_avg$features <- rownames(arabi_xylem_bins_avg)

print(paste0(nrow(arabi_xylem_bins_avg), " genes in arabidopsis - Xylem"))

pop_xylem_bins <- test_function(tradeseq_values_db1 = "TRADEseq_Xylem_Poplar.csv",
                                 rds_data_db1 = "Poplar_Xylem_TI.rds",
                                 pvalue_value_db1 = PVALUE_THRES,
                                 N_GROUPS = N_GROUPS,
                                 starting_cluster = starting_cluster,
                                 bar_plot_prefix = paste0(prefix, "_", "Poplar_Xylem"))

pop_xylem_bins_avg <- data.frame(pop_xylem_bins[[1]])
colnames(pop_xylem_bins_avg)[1:ncol(pop_xylem_bins_avg)] <- paste0("bin_pop", seq(1, N_GROUPS))
pop_xylem_bins_avg$features <- rownames(pop_xylem_bins_avg)

print(paste0(nrow(pop_xylem_bins_avg), " genes in poplar - Xylem"))

# Identify common genes in both species
# Reads the Phytozome dataset
orthofinder <-  vroom::vroom("Mappings_Populus_1to1_Arabidopsis_oct_28.txt") %>%
    dplyr::select(populus_trichocarpa, arabidopsis_thaliana) %>%
    dplyr::rename( feature_id = populus_trichocarpa ) %>%
    dplyr::rename(Arabidopsis = arabidopsis_thaliana) %>%
    dplyr::mutate(source = "Orthofinder")

Phyto <- vroom::vroom("Ptrichocarpa_533_v4.1.annotation_info.csv", 
                      col_names = T) %>%
    dplyr::rename( feature_id = Gene ) %>%
    dplyr::mutate( combine = paste0(feature_id, "_", Arabidopsis) ) %>%
    dplyr::filter( !duplicated(combine) ) %>%
    dplyr::filter( !is.na(Arabidopsis) ) %>%
    dplyr::mutate( Arabidopsis2 = make.unique(Arabidopsis) ) %>%
    dplyr::select("feature_id", "Arabidopsis2") %>%
    dplyr::rename(Arabidopsis = Arabidopsis2) %>%
    dplyr::filter( ! feature_id %in% orthofinder$feature_id ) %>%
    dplyr::mutate(source = "Phytozome")

Phyto_ortho_combined <- rbind(orthofinder, Phyto) %>%
    dplyr::mutate(Arabidopsis = make.unique(Arabidopsis))

write.table(Phyto_ortho_combined,
            "poplar_to_arabi_using_phytozome_list.csv",
            col.names = T,
            row.names = F,
            sep = ",",
            quote = F)

Phyto_ortho_combined <- Phyto_ortho_combined %>%
    dplyr::select(- source)

to_rename <- pop_xylem_bins_avg %>%
    dplyr::filter( grepl(features,
                         pattern = "Potri") ) %>%
    dplyr::rename(feature_id = features)

to_keep <- pop_xylem_bins_avg %>%
    dplyr::filter( !grepl(features,
                          pattern = "Potri") ) %>%
    dplyr::rename(feature_id = features)

# This will rename any poplar gene that wasnt in the Orthofinder list, but is in the phytozome list.
to_rename_ids <- as.data.frame(to_rename[, ncol(to_rename)])
colnames(to_rename_ids) <- "feature_id"

## Get the names of the genes
to_rename2 <- merge(to_rename_ids, Phyto_ortho_combined, by = "feature_id",  all.x = T)

## Select only the genes that have ID on arabidopsis
to_rename2 <- to_rename2[!is.na(to_rename2$Arabidopsis), ] %>%
    dplyr::mutate( combine = paste0(feature_id, "_", Arabidopsis) ) %>%
    dplyr::filter(!duplicated(combine))

to_rename2 <- merge(to_rename, to_rename2 , by = "feature_id", all.x = T) %>%
    dplyr::rename(feature_id2 = feature_id) %>%
    dplyr::select(- c( "combine") ) %>%
    dplyr::rename(feature_id = Arabidopsis)

no_match_genes <- to_rename2 %>%
    dplyr::filter(is.na(feature_id)) %>%
    dplyr::select(-feature_id) %>%
    dplyr::rename(features = feature_id2) %>%
    dplyr::relocate(features, .after = last_col())

to_rename2 <- to_rename2 %>%
    dplyr::select(-feature_id2) %>%
    dplyr::filter(!is.na(feature_id)) 

rownames(to_rename2) <- to_rename2$feature_id

pop_xylem_bins_avg2 <- dplyr::bind_rows(to_keep,
                                         to_rename2) %>%
    dplyr::rename(features = feature_id)

rownames(pop_xylem_bins_avg2) <- pop_xylem_bins_avg2$features

## Duplicates the arabidopsis entrances of genes which more than one copy are present in the poplar dataset, so we can join the heatmap.

dup_features <- pop_xylem_bins_avg2 %>%
    dplyr::select(features) %>%
    dplyr::filter(grepl(features, pattern = "\\.")) %>%
    dplyr::mutate(features2 = features) %>%
    tidyr::separate(features, into = c("features", "number") ) %>%
    dplyr::select(features, features2)

dup_features_arab <- merge(dup_features, arabi_xylem_bins_avg,  by = "features") %>%
    dplyr::select(-features) %>%
    dplyr::rename(features = features2)

arabi_xylem_bins_avg2 <- dplyr::bind_rows(arabi_xylem_bins_avg,
                                           dup_features_arab)

arabi_xylem_bins_avg2$features <- make.unique(arabi_xylem_bins_avg2$features)
rownames(arabi_xylem_bins_avg2) <- arabi_xylem_bins_avg2$features

combined_avg_bins <- merge(arabi_xylem_bins_avg2, pop_xylem_bins_avg2, by = "features")

combined_avg_bins_names <- combined_avg_bins%>%
    tidyr::separate(features, into = c("name", "number")) %>%
    dplyr::select("name")

print( paste0("There are ", nrow(combined_avg_bins), " common genes in the xylem sample") )

unique_poplar <- pop_xylem_bins_avg2 %>%
    dplyr::filter(!features %in% combined_avg_bins$features) %>%
    dplyr::mutate(features2 = features) %>%
    tidyr::separate(features2, into = c("name", "number")) %>%
    dplyr::filter(!name %in% combined_avg_bins_names$name) %>%
    dplyr::select(-name, - number)

unique_poplar <- unique(rbind(unique_poplar, no_match_genes))

write.table(unique_poplar,
            "tradeseq_sig_uniquely_in_poplar_Xylem_COUNTS.csv",
            col.names = T,
            row.names = F,
            sep = ",",
            quote = F)

write.table(unique_poplar$features,
            "tradeseq_sig_uniquely_in_poplar_Xylem.csv",
            col.names = T,
            row.names = F,
            sep = ",",
            quote = F)

## Generates the correlation for each gene separately. besides the gene name, the first 30 columns are the Arabidopsis data and the other 30 are the poplar.
gene_cor_df_comp <- data.frame()
for( i in 1:nrow(combined_avg_bins) ) {
    
    gene_cor <- cor.test(as.numeric(combined_avg_bins[i, 2:(N_GROUPS + 1)] ),
                         as.numeric(combined_avg_bins[i, (N_GROUPS + 2):ncol(combined_avg_bins)]),
                         method = "pearson" )
    gene_cor_df <- data.frame(gene = combined_avg_bins[i,1],
                              corr = gene_cor$estimate,
                              p_value = gene_cor$p.value)
    
    gene_cor_df_comp <- rbind(gene_cor_df_comp,
                              gene_cor_df)
    
}

## Multiple test (FDR) correction of the p-values
gene_cor_df_comp$fdr <- format( p.adjust(p = gene_cor_df_comp$p_value,
                                         method = "fdr",
                                         n = length(gene_cor_df_comp$p_value) ),
                                scientific = F)

gene_cor_df_comp_sig <- gene_cor_df_comp %>% 
    dplyr::filter(fdr < PVALUE_THRES) %>%
    dplyr::filter(corr > 0.5)

gene_cor_df_comp_diff <- gene_cor_df_comp[gene_cor_df_comp$fdr >= PVALUE_THRES, ]

write.csv(gene_cor_df_comp, paste0(prefix, "Xylem_all_correlations.csv"), row.names = F)

# Heatmap generation

## @@ Correlated genes

## Selecting the expression values of the genes that are correlated to generate the heatmap
corr_genes_arab <- arabi_xylem_bins_avg2 %>%
    dplyr::filter(features %in% gene_cor_df_comp_sig$gene) %>%
    dplyr::arrange(features) %>%
    dplyr::select(-features)

corr_genes_pop <- pop_xylem_bins_avg2 %>%
    dplyr::filter(features %in% gene_cor_df_comp_sig$gene) %>%
    dplyr::arrange(features) %>%
    dplyr::select(-features)

### @@@ >>>> dynutils::scale_quantile Cut off outer quantiles and rescale to a [0, 1] range
corr_genes_arab_scaled <- t(dynutils::scale_quantile(t(corr_genes_arab)))
corr_genes_pop_scaled <- t(dynutils::scale_quantile(t(corr_genes_pop)))

col_fun = colorRamp2( breaks = c(0, 0.5, 1),
                      c("#ffffcc", "#fd8d3c", "#800026") )

heatmap_corr_arab <- ComplexHeatmap::Heatmap(corr_genes_arab_scaled,
                                             border = TRUE,
                                             #rect_gp = grid::gpar(col = "white", lwd = 0.01),
                                             name = "Expression",
                                             cluster_rows = T,
                                             cluster_columns = F,
                                             col = col_fun,
                                             show_row_names = T,
                                             show_column_names = F,
                                             show_row_dend = F,
                                             use_raster = F,
                                             row_names_gp = grid::gpar(fontsize = 10),
                                             column_title = "Arabidopsis")

heatmap_corr_pop <- ComplexHeatmap::Heatmap(corr_genes_pop_scaled,
                                            border = TRUE,
                                            #rect_gp = grid::gpar(col = "white", lwd = 0.01),
                                            name = "Expression",
                                            cluster_rows = T,
                                            cluster_columns = F,
                                            show_row_names = T,
                                            show_column_names = F,
                                            show_row_dend = F,
                                            use_raster = F,
                                            col = col_fun,
                                            row_names_gp = grid::gpar(fontsize = 10),
                                            column_title = "Populus")

svg( paste0(prefix,"Xylem_correlated_genes_heatmaps.svg"),
     width = 10,
     height = 30,
     pointsize = 12)
draw( heatmap_corr_arab + heatmap_corr_pop )
dev.off()

write.csv(gene_cor_df_comp_sig,
          paste0(prefix,"list_of_genes_positive_correlation_xylem.csv"),
          row.names = F)

write.csv(gene_cor_df_comp_diff,
          paste0(prefix,"list_of_genes_UNcorrelation_xylem.csv"),
          row.names = F)

### Genes that does not correlate
## Selecting the expression values of the genes that are correlated to generate the heatmap
uncorr_genes_arab <- arabi_xylem_bins_avg2 %>%
    dplyr::filter(features %in% gene_cor_df_comp_diff$gene) %>%
    dplyr::arrange(features) %>%
    dplyr::select(-features)

uncorr_genes_pop <- pop_xylem_bins_avg2 %>%
    dplyr::filter(features %in% gene_cor_df_comp_diff$gene)%>%
    dplyr::arrange(features) %>%
    dplyr::select(-features)

uncorr_genes_arab_scaled <- t(scale(t(uncorr_genes_arab)))
uncorr_genes_pop_scaled <- t(scale(t(uncorr_genes_pop)))

heatmap_uncorr_arab <- ComplexHeatmap::Heatmap(uncorr_genes_arab_scaled,
                                               border = TRUE,
                                               #rect_gp = grid::gpar(col = "white", lwd = 0.01),
                                               name = "Expression",
                                               cluster_rows = T,
                                               cluster_columns = F,
                                               show_row_names = T,
                                               show_column_names = F,
                                               show_row_dend = F,
                                               use_raster = F,
                                               col = col_fun,
                                               row_names_gp = grid::gpar(fontsize = 10),
                                               column_title = "Arabidopsis")

heatmap_uncorr_pop <- ComplexHeatmap::Heatmap(uncorr_genes_pop_scaled,
                                              border = TRUE,
                                              #rect_gp = grid::gpar(col = "white", lwd = 0.01),
                                              name = "Expression",
                                              cluster_rows = T,
                                              cluster_columns = F,
                                              show_row_names = T,
                                              show_column_names = F,
                                              show_row_dend = F,
                                              use_raster = F,
                                              col = col_fun,
                                              row_names_gp = grid::gpar(fontsize = 10),
                                              column_title = "Populus")

svg( paste0(prefix,"Xylem_UNcorrelated_genes_heatmaps.svg"), width = 10, height = 30, pointsize = 12)
draw( heatmap_uncorr_arab + heatmap_uncorr_pop )
dev.off()
