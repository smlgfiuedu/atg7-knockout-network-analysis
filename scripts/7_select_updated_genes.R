# select a subset of significant genes to run banjo on
library(openxlsx)
library(stringr)

# read in data ------------------------------------------------------------
load("./analysis/differential_expression/output/contrast_fit_wt_vs_ko.rda")
gene_exp_table <- read.xlsx("./analysis/differential_expression/output/differentially_expressed_mouse_liver_genes_atg7_ko.xlsx",
                            rowNames = TRUE)

genelist <- read.table("./data/GSE65174/updated_genes_of_interest.txt")[[1]]


# extract top genes and genes of interest ------------------------------------
deg_intersect <- intersect(gene_exp_table$GENE_SYMBOL,genelist)
genes <- which(gene_exp_table$GENE_SYMBOL %in% deg_intersect)

genes_of_interest <- gene_exp_table[genes,]

# # check smallest lfc for threshold
# sig_lfc <- gene_exp_table[which(gene_exp_table$adj.P.Val <=.05),]
# 
# min(abs(sig_lfc$logFC[1:200]))
# min(abs(genes_of_interest$logFC))
# 
# rm(sig_lfc)


top_l2fc <- topTable(contrast_fit, coef = 1, number = Inf, adjust.method = "BH", p.value = .05, sort.by = "logFC", lfc = 35)



# check how many genes of interest are present in the top 200
top_intersect <- intersect(top_l2fc$GENE_SYMBOL, genelist)


# merge genes of interest with top 200 genes ------------------------------
already_included <- which(genes_of_interest$GENE_SYMBOL %in% top_intersect)
genes_of_interest_unique <- genes_of_interest[-already_included,]

final_genelist <- rbind(top_l2fc, genes_of_interest_unique)

# look at some stats for the final list
summary(abs(final_genelist$logFC))
summary(final_genelist$adj.P.Val)


# check what goi were included
genelist %in% final_genelist$GENE_SYMBOL
genelist[!genelist %in% final_genelist$GENE_SYMBOL]

genelist[!genelist %in% gene_exp_table$GENE_SYMBOL] # not found in the data
raw.data <- read.xlsx("./data/GSE65174/GSE65174_series_flox-Alb-Cre mice-cleaned.xlsx")
genelist[!genelist %in% trimWhiteSpace(raw.data$GENE_SYMBOL)] # not found in the raw data
# write to file -----------------------------------------------------------

write.xlsx(final_genelist, file = "./analysis/differential_expression/output/updated_goi_de_genes_mouse_liver_network.xlsx", rowNames = TRUE, colNames = TRUE)
