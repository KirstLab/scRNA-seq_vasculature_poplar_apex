suppressMessages( require(vroom) )
suppressMessages( require(tidyverse) )
set.seed(1407)

###################################################
## Identify common genes in both species
# Reads the Phytozome dataset
orthofinder <-  vroom::vroom("Mappings_Populus_1to1_Arabidopsis_oct_28.txt") %>%
    dplyr::select(populus_trichocarpa, arabidopsis_thaliana) %>%
    dplyr::rename( feature_id = populus_trichocarpa ) %>%
    dplyr::rename(Arabidopsis = arabidopsis_thaliana)

Phyto <- vroom::vroom("Ptrichocarpa_533_v4.1.annotation_info.csv", 
                      col_names = T) %>%
    dplyr::rename( feature_id = Gene ) %>%
    dplyr::mutate( combine = paste0(feature_id, "_", Arabidopsis) ) %>%
    dplyr::filter( !duplicated(combine) ) %>%
    dplyr::filter( !is.na(Arabidopsis) ) %>%
    dplyr::mutate( Arabidopsis2 = make.unique(Arabidopsis) ) %>%
    dplyr::select("feature_id", "Arabidopsis2") %>%
    dplyr::rename(Arabidopsis = Arabidopsis2) %>%
    dplyr::filter( ! feature_id %in% orthofinder$feature_id ) 

Phyto_ortho_combined <- rbind(orthofinder, Phyto) %>%
    dplyr::mutate(Arabidopsis = make.unique(Arabidopsis))

## Rename the Xylem genes
poplar_xylem_pop <- vroom::vroom(file = "TRADEseq_Xylem_Poplar.csv",
                                 col_names = T,
                                 delim = ",",
                                 show_col_types = FALSE) %>%
    dplyr::filter( grepl(genes,
                         pattern = "Potri") ) %>%
    dplyr::rename(feature_id = genes)

poplar_xylem_pop2 <- merge(poplar_xylem_pop, Phyto_ortho_combined, by = "feature_id")

poplar_xylem_remaining <- vroom::vroom(file = "TRADEseq_Xylem_Poplar.csv",
                                       col_names = T,
                                       delim = ",",
                                       show_col_types = FALSE) %>%
    dplyr::rename(feature_id = genes) %>%
    dplyr::filter( !feature_id %in% poplar_xylem_pop2$feature_id )

poplar_xylem_pop <- poplar_xylem_pop2 %>%
    dplyr::select(Arabidopsis, waldStat, df, pvalue, meanLogFC, fdr) %>%
    dplyr::rename(feature_id = Arabidopsis)

new_poplar <- rbind(poplar_xylem_pop, poplar_xylem_remaining)

write.table(new_poplar,
            "TRADEseq_poplar_xylem_using_phytozome.csv", 
            col.names = T,
            row.names = F,
            sep = ",",
            quote = F)

## Rename the Phloem genes
poplar_phloem_pop <- vroom::vroom(file = "TRADEseq_Phloem_Poplar.csv",
                                  col_names = T,
                                  delim = ",",
                                  show_col_types = FALSE) %>%
    dplyr::filter( grepl(genes,
                         pattern = "Potri") ) %>%
    dplyr::rename(feature_id = genes)

poplar_phloem_pop2 <- merge(poplar_phloem_pop, Phyto_ortho_combined, by = "feature_id")

poplar_phloem_remaining <- vroom::vroom(file = "TRADEseq_Phloem_Poplar.csv",
                                        col_names = T,
                                        delim = ",",
                                        show_col_types = FALSE) %>%
    dplyr::rename(feature_id = genes) %>%
    dplyr::filter( !feature_id %in% poplar_phloem_pop2$feature_id )

poplar_phloem_pop <- poplar_phloem_pop2 %>%
    dplyr::select(Arabidopsis, waldStat, df, pvalue, meanLogFC, fdr) %>%
    dplyr::rename(feature_id = Arabidopsis)

new_poplar <- rbind(poplar_phloem_pop, poplar_phloem_remaining)

write.table(new_poplar,
            "TRADEseq_poplar_phloem_using_phytozome.csv", 
            col.names = T,
            row.names = F,
            sep = ",",
            quote = F)

