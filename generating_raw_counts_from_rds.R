"Usage:
    generating_raw_counts_from_rds.R (--outdir1 <O1>) (--outdir2 <O2>) <input1> <input2>
    generating_raw_counts_from_rds.R -h | --help
Options:
    -h --help   show this screen.
    <input1>    INPUT1
    <input2>    INPUT2
    --outdir1  Folder to write the raw counts of species 1
    --outdir2  Folder to write the raw counts of species 2
    
" -> doc

suppressMessages(require(DropletUtils))
suppressMessages(require(docopt))
set.seed(1407)

opts <- docopt(doc)
save.image("test.rda")

vasc_arab <- readRDS(opts$`<input1>`)
arabidopsis <- vasc_arab@assays$RNA@counts

DropletUtils::write10xCounts(arabidopsis,
                             path = opts$O1)

vasc_poplar <- readRDS(opts$`<input2>`)
poplar <- vasc_poplar@assays$RNA@counts

DropletUtils::write10xCounts(poplar,
                             path = opts$O2)
