"Usage:
    comparing_markers_among_samples.R (--outdir1 <O1>) <input1> <input2> <input3> <input4>
    comparing_markers_among_samples.R -h | --help
Options:
    -h --help   show this screen.
    <input1>    INPUT1
    <input2>    INPUT2
    <input3>    INPUT3
    <input4>    INPUT4
    --outdir1  Folder to output the expression plots - markers list

" -> doc

suppressMessages(require(patchwork))
suppressMessages(require(Seurat))
suppressMessages(require(vroom))
suppressMessages(require(tidyverse))
suppressMessages(require(docopt))
set.seed(1407)

opts <- docopt(doc)

# helper function
my_feature_plot <- function(sc_data = sc_data,
                            gene_name = gene_name) {
    
    assay_id <- "RNA"
    
    minimal <- min(sc_data[[assay_id]]@data[gene_name, ])
    maximal <- max(sc_data[[assay_id]]@data[gene_name, ])
    
    Seurat::FeaturePlot(sc_data,
                        features = gene_name,
                        reduction = "umap",pt.size = 0.4) +
        scale_colour_gradient2(limits=c(minimal, maximal),
                               midpoint = maximal / 2,
                               low = "#ffffcc",
                               mid = "#fd8d3c",
                               high = "#800026"
        ) 
    
}

## Loading the datasets
vasc_arab <- readRDS(opts$`<input1>`)
vasc_poplar <- readRDS(opts$`<input2>`)

markers_of_vasculature <- vroom(opts$`<input3>`,
                                col_names = T) %>%
    rename(populus_trichocarpa = gene)

markers_of_vasculature <- markers_of_vasculature[markers_of_vasculature$Arabidopsis %in% rownames(vasc_arab), ]

for( i in 1:nrow(markers_of_vasculature)) {
    
    gene_name_arab <- markers_of_vasculature[i, 4]
    gene_name_pop <- markers_of_vasculature[i, 1]
    
    cluster_n <- markers_of_vasculature[i, 2]
    
    A <- my_feature_plot(sc_data = vasc_arab,
                         gene_name = as.character(gene_name_arab) )
    
    B <- my_feature_plot(sc_data = vasc_poplar,
                         gene_name = as.character(gene_name_pop) )
    
    plot <- B | A
    plot <- plot +
        plot_annotation( BoldTitle(as.character(markers_of_vasculature[i, 3]) ) )
    
    ggsave(paste0(opts$O1,
                  "/expression_of_the_gene_",
                  gene_name_arab, "_vs_", gene_name_pop,
                  "_marker_of_the_cluster_",
                  cluster_n, 
                  ".png"),
           plot,
           height = 5,
           width = 10)
    
    ggsave(paste0(opts$O1,
                  "/expression_of_the_gene_",
                  gene_name_arab, "_vs_", gene_name_pop,
                  "_marker_of_the_cluster_",
                  cluster_n, 
                  ".svg"),
           plot,
           height = 5,
           width = 10)
    
}