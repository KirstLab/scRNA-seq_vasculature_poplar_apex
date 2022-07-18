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

scale_data_funct <- function(tradeseq_values_db = tradeseq_values_db,
                             rds_data_db = rds_data_db,
                             fdr_treshold = fdr_treshold,
                             starting_cluster = starting_cluster) {
    
    list_of_genes_db1 <- read.table(tradeseq_values_db,
                                    sep = ",",
                                    header = T) %>%
        dplyr::filter(fdr < fdr_treshold)
    
    sc_data_traj_inf <- readRDS(rds_data_db)
    
    # Not necessary to change the code from here, but to adjust the size of the images at the very end.
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
    
    # colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    # plotcol <- colors[cut(sds$slingPseudotime_1, breaks=100)]
    # plot(reducedDims(sds)$PCA, col = plotcol, pch=16, asp = 1)
    # lines(SlingshotDataSet(sds), lwd=2, col='black')
    
    pst.ord <- order(sds$slingPseudotime_1, na.last = NA)
    heatdata <- assays(sds)$counts
    heatdata <- heatdata[rownames(heatdata) %in% list_of_genes_db1$gene, pst.ord]
    
    cells_order_in_traject <- data.frame(cells = colnames(heatdata))
    
    heatdata <- log1p(heatdata)
    
    #corr_genes_arab_scaled <- t(scale(t(heatdata)))
    corr_genes_arab_scaled <- dynutils::scale_quantile(heatdata) 
    
}

scale_data_funct_poplar <- function(tradeseq_values_db = tradeseq_values_db,
                                    rds_data_db = rds_data_db,
                                    fdr_treshold = fdr_treshold,
                                    starting_cluster = starting_cluster) {
    
    list_of_genes_db1 <- read.table(tradeseq_values_db,
                                    sep = ",",
                                    header = T) %>%
        dplyr::filter(fdr < fdr_treshold)
    
    ids_mapping <- vroom::vroom("poplar_to_arabi_using_phytozome_list.csv",
                                delim = ",") %>%
        dplyr::rename(genes = "feature_id")
    
    list_of_genes_db2 <- merge(list_of_genes_db1,
                               ids_mapping,
                               by = "genes",
                               all.x = T)
    
    list_of_genes_db2$feature_id2 <- ifelse(is.na(list_of_genes_db2$Arabidopsis),
                                            list_of_genes_db2$genes,
                                            list_of_genes_db2$Arabidopsis)
    
    sc_data_traj_inf <- readRDS(rds_data_db)
    
    # Not necessary to change the code from here, but to adjust the size of the images at the very end.
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
    
    # colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    # plotcol <- colors[cut(sds$slingPseudotime_1, breaks=100)]
    # plot(reducedDims(sds)$PCA, col = plotcol, pch=16, asp = 1)
    # lines(SlingshotDataSet(sds), lwd=2, col='black')
    
    pst.ord <- order(sds$slingPseudotime_1, na.last = NA)
    heatdata <- assays(sds)$counts
    
    heatdata <- heatdata[rownames(heatdata) %in% list_of_genes_db2$genes, pst.ord]
    
    heatdata <- as.data.frame(heatdata)
    heatdata$genes <- rownames(heatdata)
    
    # Names still the poplar IDs, so we replace with the gene names that contains the arabidopsis IDs
    
    list_of_genes_db2 <- list_of_genes_db2 %>%
        dplyr::select(c("genes", "feature_id2") )
    
    heatdata2 <- merge(heatdata, list_of_genes_db2, by = "genes")%>%
        dplyr::select(-genes) %>%
        tibble::column_to_rownames("feature_id2")
    
    cells_order_in_traject <- data.frame(cells = colnames(heatdata))
    
    heatdata2 <- log1p(heatdata2)
    
    #corr_genes_arab_scaled <- t(scale(t(heatdata)))
    corr_genes_arab_scaled <- dynutils::scale_quantile(heatdata2) 
    
}

fdr_treshold = 0.01
starting_cluster = 8

### Phloem
###################
## Correlated genes
###################

corr_genes_arab_scaled <- scale_data_funct(
    tradeseq_values_db = "TRADEseq_Phloem_Arabidopsis.csv",
    rds_data_db = "Arabidopsis_Phloem_TI.rds",
    fdr_treshold = fdr_treshold,
    starting_cluster = starting_cluster)

corr_genes_pop_scaled <- scale_data_funct_poplar(
    tradeseq_values_db = "TRADEseq_Phloem_Poplar.csv",
    rds_data_db = "Poplar_Phloem_TI.rds",
    fdr_treshold = fdr_treshold,
    starting_cluster = starting_cluster)

## This section duplicates the entrances of the arabidopsis genes that have more than one orthologs. This is necessary, becase in the heatmaps, a 1-to-1 relashionship is required.

corr_genes_arab_scaled2 <- data.frame(corr_genes_arab_scaled) %>%
    mutate(features = rownames(.))

dup_features <- corr_genes_pop_scaled %>% 
    dplyr::mutate(features = rownames(corr_genes_pop_scaled)) %>%
    dplyr::select(features) %>%
    dplyr::mutate(features2 = features) %>%
    tidyr::separate(features, into = c("features", "number") ) %>%
    dplyr::select(features, features2)

dup_features_arab <- merge(dup_features, corr_genes_arab_scaled2,  by = "features") %>%
    dplyr::select(-features) %>%
    dplyr::rename(features = features2) %>%
    dplyr::filter(grepl(features, pattern = "\\."))

#dup_features_arab <- merge(dup_features, corr_genes_arab_scaled2,  by = "features")

corr_genes_arab_scaled2 <- dplyr::bind_rows(corr_genes_arab_scaled2,
                                            dup_features_arab)

rownames(corr_genes_arab_scaled2) <- corr_genes_arab_scaled2$features

corr_genes_arab_scaled_final <- corr_genes_arab_scaled2[rownames(corr_genes_arab_scaled2) %in% rownames(corr_genes_pop_scaled), ]

corr_genes_pop_scaled_final <- corr_genes_pop_scaled[rownames(corr_genes_pop_scaled) %in% rownames(corr_genes_arab_scaled2), ]

corr_genes_pop_scaled_unique <- corr_genes_pop_scaled[!rownames(corr_genes_pop_scaled) %in% rownames(corr_genes_arab_scaled2), ]

## Sort both datasets, so the heatmap will work propely.
corr_genes_arab_scaled_final <-  corr_genes_arab_scaled_final %>%
    dplyr::arrange(rownames(.)) %>%
    dplyr::select(-features)

corr_genes_pop_scaled_final <- corr_genes_pop_scaled_final %>%
    dplyr::arrange(rownames(.))

col_fun = colorRamp2( breaks = c(0, 0.5, 1),
                      c("#ffffcc", "#fd8d3c", "#800026") )

arabidopsis_ht <- ComplexHeatmap::Heatmap(corr_genes_arab_scaled_final,
                                          border = TRUE,
                                          name = "Expression",
                                          cluster_rows = T,
                                          width = 50,
                                          cluster_columns = F,
                                          col = col_fun,
                                          show_row_names = T,
                                          show_column_names = F,
                                          show_row_dend = F,
                                          use_raster = F,
                                          row_names_gp = grid::gpar(fontsize = 10),
                                          column_title = "Arabidopsis")

poplar_ht <- ComplexHeatmap::Heatmap(corr_genes_pop_scaled_final,
                                     border = TRUE,
                                     #rect_gp = grid::gpar(col = "white", lwd = 0.01),
                                     name = "Expression",
                                     cluster_rows = T,
                                     cluster_columns = F,
                                     col = col_fun,
                                     width = 50,
                                     show_row_names = T,
                                     show_column_names = F,
                                     show_row_dend = F,
                                     use_raster = F,
                                     row_names_gp = grid::gpar(fontsize = 10),
                                     column_title = "Poplar")

poplar_ht_unique <- ComplexHeatmap::Heatmap(corr_genes_pop_scaled_unique,
                                            border = TRUE,
                                            #rect_gp = grid::gpar(col = "white", lwd = 0.01),
                                            name = "Expression",
                                            cluster_rows = T,
                                            cluster_columns = F,
                                            col = col_fun,
                                            width = 50,
                                            show_row_names = T,
                                            show_column_names = F,
                                            show_row_dend = F,
                                            use_raster = F,
                                            row_names_gp = grid::gpar(fontsize = 10),
                                            column_title = "Poplar")

svg( "All_cells_phloem_all_genes_sig_in_both_datasets.svg",
     width = 30,
     height = 60,
     pointsize = 12)
draw(arabidopsis_ht + poplar_ht)
dev.off()

svg("All_cells_phloem_genes_sig_ONLY_IN_POPLAR.svg",
     width = 15,
     height = 60,
     pointsize = 12)
draw(poplar_ht_unique)
dev.off()

### Xylem

###################
## Correlated genes
###################

corr_genes_arab_scaled <- scale_data_funct(
    tradeseq_values_db = "TRADEseq_Xylem_Arabidopsis.csv",
    rds_data_db = "Arabidopsis_Xylem_TI.rds",
    fdr_treshold = fdr_treshold,
    starting_cluster = starting_cluster)

corr_genes_pop_scaled <- scale_data_funct_poplar(
    tradeseq_values_db = "TRADEseq_Xylem_Poplar.csv",
    rds_data_db = "Poplar_Xylem_TI.rds",
    fdr_treshold = fdr_treshold,
    starting_cluster = starting_cluster)

## This section duplicates the entrances of the arabidopsis genes that have more than one orthologs. This is necessary, becase in the heatmaps, a 1-to-1 relashionship is required.

corr_genes_arab_scaled2 <- data.frame(corr_genes_arab_scaled) %>%
    mutate(features = rownames(.))

dup_features <- corr_genes_pop_scaled %>% 
    dplyr::mutate(features = rownames(corr_genes_pop_scaled)) %>%
    dplyr::select(features) %>%
    dplyr::mutate(features2 = features) %>%
    tidyr::separate(features, into = c("features", "number") ) %>%
    dplyr::select(features, features2)

dup_features_arab <- merge(dup_features, corr_genes_arab_scaled2,  by = "features") %>%
    dplyr::select(-features) %>%
    dplyr::rename(features = features2) %>%
    dplyr::filter(grepl(features, pattern = "\\."))

corr_genes_arab_scaled2 <- dplyr::bind_rows(corr_genes_arab_scaled2,
                                            dup_features_arab)

rownames(corr_genes_arab_scaled2) <- corr_genes_arab_scaled2$features

corr_genes_arab_scaled_final <- corr_genes_arab_scaled2[corr_genes_arab_scaled2$features %in% rownames(corr_genes_pop_scaled), ]

corr_genes_pop_scaled_final <- corr_genes_pop_scaled[rownames(corr_genes_pop_scaled) %in% rownames(corr_genes_arab_scaled2), ]

corr_genes_pop_scaled_unique <- corr_genes_pop_scaled[!rownames(corr_genes_pop_scaled) %in% rownames(corr_genes_arab_scaled2), ]

## Sort both datasets, so the heatmap will work propely.
corr_genes_arab_scaled_final <-  corr_genes_arab_scaled_final %>%
    dplyr::arrange(rownames(.)) %>%
    dplyr::select(-features)

corr_genes_pop_scaled_final <- corr_genes_pop_scaled_final %>%
    dplyr::arrange(rownames(.))

col_fun = colorRamp2( breaks = c(0, 0.5, 1),
                      c("#ffffcc", "#fd8d3c", "#800026") )

arabidopsis_ht <- ComplexHeatmap::Heatmap(corr_genes_arab_scaled_final,
                                          border = TRUE,
                                          name = "Expression",
                                          cluster_rows = T,
                                          width = 50,
                                          cluster_columns = F,
                                          col = col_fun,
                                          show_row_names = T,
                                          show_column_names = F,
                                          show_row_dend = F,
                                          use_raster = F,
                                          row_names_gp = grid::gpar(fontsize = 10),
                                          column_title = "Arabidopsis")

poplar_ht <- ComplexHeatmap::Heatmap(corr_genes_pop_scaled_final,
                                     border = TRUE,
                                     #rect_gp = grid::gpar(col = "white", lwd = 0.01),
                                     name = "Expression",
                                     cluster_rows = T,
                                     cluster_columns = F,
                                     col = col_fun,
                                     width = 50,
                                     show_row_names = T,
                                     show_column_names = F,
                                     show_row_dend = F,
                                     use_raster = F,
                                     row_names_gp = grid::gpar(fontsize = 10),
                                     column_title = "Poplar")

poplar_ht_unique <- ComplexHeatmap::Heatmap(corr_genes_pop_scaled_unique,
                                            border = TRUE,
                                            #rect_gp = grid::gpar(col = "white", lwd = 0.01),
                                            name = "Expression",
                                            cluster_rows = T,
                                            cluster_columns = F,
                                            col = col_fun,
                                            width = 50,
                                            show_row_names = T,
                                            show_column_names = F,
                                            show_row_dend = F,
                                            use_raster = F,
                                            row_names_gp = grid::gpar(fontsize = 10),
                                            column_title = "Poplar")

svg("All_cells_xylem_all_genes_sig_in_both_datasets.svg",
     width = 30,
     height = 30,
     pointsize = 12)
draw(arabidopsis_ht + poplar_ht)
dev.off()

svg("All_cells_xylem_genes_sig_ONLY_IN_POPLAR.svg",
     width = 15,
     height = 30,
     pointsize = 12)
draw(poplar_ht_unique)
dev.off()

