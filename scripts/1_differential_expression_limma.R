#######################################################################################################################################
# Differential expression analysis between Atg7 WT and KO samples in mouse liver microarray data (GSE65174) via limma (user guide pg42)
# "Atg7flox_flox_mice_liver" -> Atg7 wild type
# "Atg7flox_flox_Alb_Cre_mice_liver" -> Atg7 knockout
#######################################################################################################################################
library(openxlsx)
library(stringr)
library(limma)

outputFile <- "differentially_expressed_mouse_liver_genes_atg7_ko"
# read in data ------------------------------------------------------------
mouse_liver <- read.xlsx("./analysis/differential_expression/input/mouse_liver_subset.xlsx", rowNames = TRUE)
mouse_metadata <- read.table("./analysis/differential_expression/input/mouse_liver_metadata.txt", header = TRUE, sep = ",")

mouse_liver$GENE_SYMBOL <- str_trim(mouse_liver$GENE_SYMBOL)

# fit model ---------------------------------------------------------------
group <- factor(mouse_metadata$Treatment, levels = c("Atg7flox_flox_mice_liver", "Atg7flox_flox_Alb_Cre_mice_liver"))
design <- model.matrix(~ 0+group)
colnames(design) <- levels(group)


fit <- lmFit(mouse_liver, design)

contrast.matrix <- makeContrasts(Atg7WT_vs_Atg7KO = Atg7flox_flox_mice_liver-Atg7flox_flox_Alb_Cre_mice_liver, levels = design)
contrast_fit <- contrasts.fit(fit, contrast.matrix)
contrast_fit <- eBayes(contrast_fit)




# summarize results -------------------------------------------------------
full_gene_exp_table <- topTable(contrast_fit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "logFC")
results <- decideTests(contrast_fit)
summary(results)
vennDiagram(results)

save(contrast_fit, file = "./analysis/differential_expression/output/contrast_fit_wt_vs_ko.rda")
write.xlsx(full_gene_exp_table, file = paste0("./analysis/differential_expression/output/", outputFile, ".xlsx"), rowNames = TRUE, colNames = TRUE)
save.image(file = "./analysis/differential_expression/limma_deg_workspace_07262023.RData")


