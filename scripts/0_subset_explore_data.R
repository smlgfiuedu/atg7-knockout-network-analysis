library(openxlsx)
library(stringr)

mouse_liver <- read.xlsx("./data/GSE65174_series_flox-Alb-Cre mice-cleaned.xlsx")


# format samples list to use as metadata -------------------------------------------------
samples_of_interest <- read.table("./data/samples_of_interest.txt", sep = ",")
colnames(samples_of_interest) <- c("SampleName", "Treatment", "SampleNumber" )

samples_of_interest$Treatment <- str_replace_all(samples_of_interest$Treatment, "; ", "_")
samples_of_interest$Treatment <- str_replace_all(samples_of_interest$Treatment, " ", "_")
samples_of_interest$Treatment <- str_replace_all(samples_of_interest$Treatment, "/", "_")
samples_of_interest$Treatment <- str_replace_all(samples_of_interest$Treatment, "-", "_")



# subset data -------------------------------------------------------------
mouse_liver$GENE_SYMBOL <- str_trim(mouse_liver$GENE_SYMBOL)
select_samples <- c(1:2, which(colnames(mouse_liver) %in% samples_of_interest$SampleName))
mouse_liver <- mouse_liver[, select_samples]

write.xlsx(mouse_liver, file = "./data/mouse_liver_subset.xlsx")
write.table(samples_of_interest, file = "./data/mouse_liver_metadata.txt", sep = ",", quote = FALSE, row.names = FALSE)


# format gene list file ---------------------------------------------------
gene_list <- read.table("./data/genes_of_interest.txt", sep = ",")
gene_list <- str_trim(gene_list)
gene_list <- str_to_sentence(gene_list)
write.table(gene_list, file = "./data/genelist.txt", sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
