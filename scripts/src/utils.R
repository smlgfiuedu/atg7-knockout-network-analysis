
# discretize a gene expression vector based on it's Z-score --------
# Coded 0 if 1 > z > -1
#       1 if -1 >= z
#       2 if 1 >= z
discretizeGeneExpr <- function(gene_vector){
  
  is_low <- gene_vector < -1
  is_high <- gene_vector > 1
  is_none <- !is_low & !is_high
  
  gene_vector[is_none] <- 0
  gene_vector[is_low] <- 1
  gene_vector[is_high] <- 2
  
  return(gene_vector)
}



# get human-mouse gene orthologs ------------------------------------------
getHumanMouseGenes <- function(genes){
  mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  # separate human and mouse 
  mouse <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[2]]
  human <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[1]]
  
  # remove some columns
  mouse <- mouse[,c(1,4)]
  human <- human[,c(1,4)]
  
  # merge the 2 dataset  (note that the human list is longer than the mouse one)
  mh_genes <- merge.data.frame(mouse,human,by = "DB.Class.Key",all.y = TRUE) 
  colnames(mh_genes)[2:3] <- c("mgi", "hgnc")
  mh_genes <- mh_genes[!duplicated(mh_genes$hgnc),]
  mh_genes <- mh_genes[!duplicated(mh_genes$mgi),]
  mh_genes <- mh_genes[which(mh_genes$hgnc %in% genes),]
  return(mh_genes)
}