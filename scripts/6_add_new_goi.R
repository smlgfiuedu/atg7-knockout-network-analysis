# script to add new goi to the analysis
source("./scripts/src/utils.R")

new.goi <- as.character(read.table("./data/GSE65174/new_goi.txt", sep = ","))

hs.mm.pairs <- getHumanMouseGenes(new.goi)
mm.goi <- hs.mm.pairs$mgi

new.goi[!new.goi %in% hs.mm.pairs$hgnc]

# check if they were already included in the analysis
banjo.input <- read.table("./analysis/banjo/input/mouse_liver_microarray_data.txt", sep = "\t", header = TRUE)

mm.goi %in% colnames(banjo.input) # not included

# They were not included so they will have to be added to the matrix that is input to banjo
original.goi <- read.table("./data/GSE65174/genelist.txt", sep = "\n")[[1]]
updated.goi <- union(original.goi, mm.goi)
write.table(updated.goi, file = "./data/GSE65174/updated_genes_of_interest.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
