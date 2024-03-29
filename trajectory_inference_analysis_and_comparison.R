"Usage:
    trajectory_inference_analysis_and_comparison.R (--group_of_clusters <gc>) (--group_name <gn>) (--starting_cluster <sc>) (--species_name1 <sn1>) (--species_name2 <sn2>) <input1>
    trajectory_inference_analysis_and_comparison.R -h | --help
Options:
    -h --help   show this screen.
    <input1>    RDS file containing the integrated and clustered dataset (generated by Asc-Seurat)
    --group_of_clusters  TEST    
    --group_name    test
    --starting_cluster  TEST
" -> doc

suppressMessages(require(docopt))
suppressMessages(require(Seurat))
suppressMessages(require(tidyverse))
suppressMessages(require(dyno))
suppressMessages(require(grDevices))
suppressMessages(require(slingshot))
suppressMessages(require(RColorBrewer))
suppressMessages(require(ComplexHeatmap))
suppressMessages(require(grDevices))
suppressMessages(require(SingleCellExperiment))
suppressMessages(require(dynplot))
suppressMessages(require(dynwrap))
suppressMessages(require(dynfeature))
suppressMessages(require(ggthemes))
set.seed(1407)

opts <- docopt(doc)
save.image("trajectory_inference_analysis_and_comparison.rda")
#load("trajectory_inference_analysis_and_comparison.rda")

# Loads the dataset
sc_data <- readRDS(opts$`<input1>`)
#DimPlot(sc_data, label = T)

groups <- strsplit(opts$gc, split = ",")
names(groups) <- opts$gn
STARTING_CLUSTER = opts$`<sc>`

#################################
### Generating the Phloem RDS ###
#################################

for (g in 1:length(groups) ) {
    
    to_filter <- base::subset( sc_data,
                               idents = groups[[g]] )
    
    to_filter_ch <- to_filter@meta.data
    to_filter_ch <- base::rownames(to_filter_ch)
    
    sc_data_traj_inf <- base::subset(sc_data,
                                     cells = to_filter_ch)
    
    DefaultAssay(sc_data_traj_inf) <- "RNA"
    
    ### Filtering by the species
    to_filter_spec1 <- sc_data_traj_inf@meta.data %>%
        dplyr::filter(treat == opts$sn1) %>%
        rownames()
    
    to_filter_spec2 <- sc_data_traj_inf@meta.data %>%
        dplyr::filter(treat == opts$sn2) %>%
        rownames()
    
    sc_data_traj_spec1 <- base::subset(sc_data_traj_inf,
                                       cells = to_filter_spec1)
    sc_data_traj_spec2 <- base::subset(sc_data_traj_inf,
                                       cells = to_filter_spec2)
    
    species_data <- list(spec1 = sc_data_traj_spec1,
                         spec2 = sc_data_traj_spec2)
    
    names(species_data) <- c(opts$sn1, opts$sn2)
    
    saveRDS(sc_data_traj_spec1,
            file = paste0(opts$sn1,
                          "_",
                          names(groups)[g],
                          "_TI.rds"),
            compress = T)
    
    saveRDS(sc_data_traj_spec2,
            file = paste0(opts$sn2,
                          "_",
                          names(groups)[g],
                          "_TI.rds"),
            compress = T)
    
    rm(sc_data_traj_spec1, sc_data_traj_spec2)
    
    for (t in 1:length(species_data) ) {
        
        sing_cell_data <- species_data[[t]]
        
        COUNTS <- as.matrix(sing_cell_data@assays$RNA@counts)
        COUNTS <- COUNTS[rowSums(COUNTS) > 0, ]
        
        DATA_counts <- as.matrix(sing_cell_data@assays$RNA@data)
        DATA_counts <- DATA_counts[rowSums(DATA_counts) > 0, ]
        
        print( paste0("In the sample of ",
                      paste0(names(species_data)[t], ", "),
                      paste0(names(groups)[g], ", "),
                      "a total of ",
                      nrow(DATA_counts),
                      " genes and ",
                      ncol(DATA_counts),
                      " cells were used.") )
        
        expressed_genes <- rownames(DATA_counts)
        
        object_expression <- Matrix::t(as(DATA_counts, 'sparseMatrix'))
        object_counts <-  Matrix::t(as(COUNTS, 'sparseMatrix'))
        
        dataset_inf <- wrap_expression(
            
            counts = object_counts,
            expression =  object_expression
            
        )
        
        end_cells <- NULL
        
        sc_meta <-  sing_cell_data[[]] %>%
            mutate(cells = rownames(.))
        
        sc_meta_cluster <- data.frame(cell_id = sc_meta$cells,
                                      group_id = sc_meta$seurat_clusters)
        
        dimred <- sing_cell_data@reductions$pca@cell.embeddings
        
        start_cells <- sc_meta[sc_meta$seurat_clusters == STARTING_CLUSTER, ]
        start_cells <- as.character(start_cells$cells)
        
        # fetch newest version of the method
        method_id <- paste0("dynverse/ti_", "slingshot", ":latest")
        methods_selected <- create_ti_method_container(method_id)
        
        dataset <- add_prior_information(
            
            dataset_inf,
            start_id = start_cells,
            end_id = end_cells,
            groups_id = sc_meta_cluster,
            dimred = dimred
        )
        
        sds <- infer_trajectory(dataset,
                                methods_selected(),
                                verbose = T,
                                give_priors = c("start_id",
                                                #"end_id",
                                                "groups_id",
                                                "dimred"))
        
        ## >>> Parei aqui!        ## Ggsave those
        saveRDS(sds,
                paste0("slingshot_TI_",
                       names(species_data)[t],
                       "_",
                       names(groups)[g],
                       ".rds"),
                compress = T)
        
        p <- dynplot::plot_dimred(sds,
                                  grouping = sc_meta_cluster,
                                  color_density = "grouping") +
            ggtitle("Cell grouping")
        
        ggsave(p, 
               filename = paste0("TI_scatter_",
                                 names(species_data)[t],
                                 "_",
                                 names(groups)[g],
                                 ".svg"),
               width = 10,
               height = 8)
        
        p <- dynplot::plot_dendro(sds,
                                  grouping = sc_meta_cluster) +
            ggtitle("Trajectory")
        
        ggsave(p, 
               filename = paste0("TI_dendrogram_",
                                 names(species_data)[t],
                                 "_",
                                 names(groups)[g],
                                 ".svg"),
               width = 10,
               height = 8)
        
        # p <- dynplot::plot_graph(sds,
        #                          grouping = sc_meta_cluster,
        #                          expression_source = dataset) +
        #     ggtitle("Trajectory represented as a graph")
        # 
        # ggsave(p, 
        #        filename = paste0("TI_graph_",
        #                          names(species_data)[t],
        #                          "_",
        #                          names(groups)[g],
        #                          ".svg"),
        #        width = 10,
        #        height = 8)
    }
}
