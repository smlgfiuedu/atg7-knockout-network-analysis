# prep microarray data for input to banjo
library(openxlsx)
source("./scripts/src/utils.R")

de_genes <- read.xlsx("./analysis/differential_expression/output/de_genes_mouse_liver_network_nov17.xlsx", rowNames = TRUE, colNames = TRUE)
mouse_liver <- read.xlsx("./data/GSE65174/mouse_liver_subset.xlsx")
metadata <- read.table("./data/GSE65174/mouse_liver_metadata.txt", header = TRUE, sep = ",")

# extract differentially expressed genes from the mouse liver microarray data --------
deg <- which(mouse_liver$GENE_SYMBOL %in% de_genes$GENE_SYMBOL)
mouse_deg <- mouse_liver[deg,]
all(mouse_deg$GENE_SYMBOL %in% de_genes$GENE_SYMBOL) #TRUE


microarray_data <- t(as.matrix(mouse_deg[,-(1:2)]))
colnames(microarray_data) <- mouse_deg$GENE_SYMBOL
microarray_data[1:4,1:4]

microarray_data <- scale(microarray_data)
microarray_data[1:4,1:4]




# discretize microarray values --------------------------------------------
microarray_discrete <- apply(microarray_data, MARGIN = 2, discretizeGeneExpr)

#check
microarray_data[1:5,1:5]
microarray_discrete[1:5,1:5]

microarray_final <- as.data.frame(microarray_discrete)

all(rownames(microarray_final)==metadata$SampleName) #check order

metadata$Condition <- factor(metadata$Treatment, levels = c("Atg7flox_flox_mice_liver", "Atg7flox_flox_Alb_Cre_mice_liver"),
                             labels = c("wild_type", "knockout"))

metadata$dummyCondition <- as.numeric(metadata$Condition) - 1
microarray_final$Experiment <- metadata$dummyCondition
microarray_final <- microarray_final[,c(ncol(microarray_final), 1:ncol(microarray_final)-1)]

write.table(microarray_final, file = "./analysis/banjo/input/nov17_mouse_liver_microarray_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#check
ma_data <- read.table("./analysis/banjo/input/nov17_mouse_liver_microarray_data.txt", sep = "\t", header = TRUE)
