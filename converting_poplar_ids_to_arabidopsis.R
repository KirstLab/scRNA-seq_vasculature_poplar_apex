"Usage:
    converting_poplar_ids_to_arabidopsis.R (--out <O1>) <input1> <input2>
    converting_poplar_ids_to_arabidopsis.R -h | --help
Options:
    -h --help   show this screen.
    <input1>    INPUT1
    <input2>    INPUT2
    --out  File containing the new renamed features

" -> doc

suppressMessages(require(vroom))
suppressMessages(require(tidyverse))
suppressMessages(require(docopt))
set.seed(1407)

opts <- docopt(doc)
save.image("converting_poplar_ids_to_arabidopsis.rda")
#load("converting_poplar_ids_to_arabidopsis.rda")

##########################################
#####   Using a list of orthologs    #####
##########################################

ortho <- vroom::vroom(opts$`<input1>`,
                      col_names = T) %>%
    dplyr::select(arabidopsis_thaliana, populus_trichocarpa)

## Renaming the Poplar genes according with their orthologs.
renaming_features <- vroom::vroom(opts$`<input2>`,
                           col_names = F)

colnames(renaming_features)[1] <- colnames(ortho)[2]

renaming_features2 <- merge(renaming_features, ortho, by = "populus_trichocarpa", all.x = T)

renaming_features2$Arabidopsis <- ifelse( is.na(renaming_features2$arabidopsis_thaliana), 
                                          renaming_features2$populus_trichocarpa,
                                          renaming_features2$arabidopsis_thaliana)

renaming_features2 <- renaming_features2[, c(4,4)]

write.table(renaming_features2,
            opts$O1,
            col.names = F,
            row.names = F,
            sep = "\t",
            quote = F)