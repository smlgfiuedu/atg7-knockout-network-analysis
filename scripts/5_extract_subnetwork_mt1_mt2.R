library(tidyverse)
library(bnlearn)
source("./scripts/dot_to_str.R")

dotfile <- "./results/raw/mouse.liver.32h.ma.tg.highlight-interestgenes.2023.08.02.17.23.44.txt"
out_prefix <- "./results/raw/mouse_liver_subset_mt1_mt2"

genelist <- read.table("./data/genelist.txt")[[1]]
mliver_data <- read.table("./data/mouse_liver_microarray_data.txt", sep = "\t", header = TRUE)


# recode array data -------------------------------------------------------
colnames(mliver_data)[1] <- "Atg7KO"
mliver_data$Atg7KO <- factor(mliver_data$Atg7KO, labels = c("Atg7WT", "Atg7KO"))
relabel_gene_vector <- function(gene_vector){
  
  is_none <- gene_vector == 0
  is_low <- gene_vector == 1
  is_high <- gene_vector == 2
  
  
  gene_vector[is_none] <- "No Change"
  gene_vector[is_low] <- "Low"
  gene_vector[is_high] <- "High"
  
  gene_vector <- factor(gene_vector, levels = c("No Change", "Low", "High"))
  return(gene_vector)
  
}

for (l in 2:ncol(mliver_data)){
  mliver_data[,l] <- relabel_gene_vector(mliver_data[,l])
}

# convert dot file to model string and bn object --------------------------
str <- TransDottoString(dotfile)
net <- model2network(str)


#relabel experiment to specify condition & transcripts that begin with numbers so R will stop crying and match with dataframe
exp_node <- which(names(net$nodes)=="Experiment")
number_nodes <- grep("^[0-9]", names(net$nodes))
h2_node <- grep("-", names(net$nodes))

new_labels <- names(net$nodes)
new_labels[exp_node] <- "Atg7KO"
new_labels[number_nodes] <- paste0("X", names(net$nodes)[number_nodes])
new_labels[h2_node] <- gsub("-", ".", names(net$nodes)[h2_node])
net <- rename.nodes(net, new_labels)
names(net$nodes)
rm(new_labels, dotfile, str)


# fit entire model --------------------------------------------------------
ordered <- match(names(net$nodes), colnames(mliver_data))
mliver_data <- mliver_data[ ,ordered]
all(names(net$nodes)==colnames(mliver_data))

full.fit <- bn.fit(net, data = mliver_data, method = "bayes", iss = 0.1)
full.fit.mle <- bn.fit(net, data = mliver_data, method = "mle", replace.unidentifiable = TRUE)

full.fit.mle$Atg7KO$prob
score(net, data = mliver_data)
# extract neighbors for markov blanket --------------------------------------------------
atg_mblanket <- mb(net, node = "Atg7KO")

# add MT1 and MT2 to the network:
extract_neighbors <- function(network, node_name){
  node <- network$nodes[[which(names(network$nodes)==node_name)]]
  nbr <- node$nbr
}

mt1_nbr <- extract_neighbors(net, "Mt1")
mt2_nbr <- extract_neighbors(net, "Mt2")
mt_neighbors <- union(mt1_nbr, mt2_nbr)
# parents through great grandparents --------------------------------------
tier1 <- sapply(atg_mblanket, FUN = function(x) parents(net, x)) %>% unlist() %>% unique()
tier2 <- sapply(tier1, FUN = function(x) parents(net, x)) %>% unlist() %>% unique() 
tier3 <- sapply(tier2, FUN = function(x) parents(net, x)) %>% unlist() %>% unique() 

# check which genes of interest are within 3 generations of the mk --------
first_gen_goi <- intersect(genelist, tier1)
second_gen_goi <- intersect(genelist, tier2)
third_gen_goi <- intersect(genelist, tier3)

parent_goi_to_include <- unique(c(first_gen_goi, second_gen_goi, third_gen_goi)) #4 


# for second gen, find which children are parents of mkvs
second_gen_children <- sapply(second_gen_goi, FUN = function(x) children(net, x)) %>% unlist() %>% unique()
gen2c_link <- intersect(second_gen_children, tier1)

# for the third gen, find which children are the parents of  the mkv parents
third_gen_children <- sapply(third_gen_goi, FUN = function(x) children(net, x)) %>% unlist() %>% unique()
gen3c_link <- intersect(third_gen_children, tier2)


parent_goi_link <- union(gen2c_link, gen3c_link)
rm(gen2c_link, gen3c_link)



# children through great grandchildren ------------------------------------
tier1 <- sapply(atg_mblanket, FUN = function(x) children(net, x)) %>% unlist() %>% unique()
tier2 <- sapply(tier1, FUN = function(x) children(net, x)) %>% unlist() %>% unique() 
tier3 <- sapply(tier2, FUN = function(x) children(net, x)) %>% unlist() %>% unique() 

first_gen_goi <- intersect(genelist, tier1)
second_gen_goi <- intersect(genelist, tier2)
third_gen_goi <- intersect(genelist, tier3)

child_goi_to_include <- unique(c(first_gen_goi, second_gen_goi, third_gen_goi)) #8

# for second gen, find which parents are children of mkvs
second_gen_parents <- sapply(second_gen_goi, FUN = function(x) parents(net, x)) %>% unlist() %>% unique()
gen2p_link <- intersect(second_gen_parents, tier1)


# for the third gen, find which parents are the children of  the mkv children
third_gen_parents <- sapply(third_gen_goi, FUN = function(x) parents(net, x)) %>% unlist() %>% unique()
gen3p_link <- intersect(third_gen_parents, tier2)

children_goi_link <- union(gen2p_link, gen3p_link)
rm(gen2p_link, gen3p_link)


# list all genes to include in the subgraph -------------------------------
all_genes <- list(
  exp_parents = net$nodes$Atg7KO$parents, #blue
  exp_children = net$nodes$Atg7KO$children, #red
  parents_children = net$nodes[[which(names(net$nodes)==net$nodes$Atg7KO$children)]]$parents, #purple
  genes_interest = c(union(parent_goi_to_include, child_goi_to_include), "Mt1", "Mt2"), #green
  markov_blanket = atg_mblanket,
  connect_genes = union(parent_goi_link, children_goi_link),
  mt_nbr = mt_neighbors)


genes <- unlist(all_genes) 
names(genes) <- gsub(names(genes), pattern = "[0-9]", replacement = "")

geneset <- data.frame(gene = unique(genes), #41; 57 with Mt1/2 and neighbors
                      atg7ko_parents = 0,
                      atg7ko_children = 0,
                      parents_children = 0,
                      genes_interest = 0,
                      markov_blanket = 0)  
for(g in 1:nrow(geneset)){
  types <- names(genes)[which(genes==geneset$gene[g])]
  
  if("exp_parents" %in% types){
    geneset[g, "atg7ko_parents"] <- 1
  }
  
  if("exp_children" %in% types){
    geneset[g, "atg7ko_children"] <- 1
  }
  
  if("parents_children" %in% types){
    geneset[g, "parents_children"] <- 1
  }
  
  if("genes_interest" %in% types){
    geneset[g, "genes_interest"] <- 1
  }
  
  if("markov_blanket" %in% types){
    geneset[g, "markov_blanket"] <- 1
  }
  
}

# calculate some numbers
length(geneset$gene[which(geneset$genes_interest == 1)]) #12 genes of interest

# assign colors to the nodes  ---------------------------------------------
geneset$color <- "white"

for(g in 1:nrow(geneset)){
  gene <- geneset$gene[g]
  if(gene == "Atg7KO"){
    geneset[g,]$color <- "#ffffb3" #pale yellow
    next
  }
  
  if(geneset[g,]$markov_blanket == 1){
    geneset[g,]$color <- "#ffffb3" #pale yellow
  }
  
  if(geneset[g,]$genes_interest == 1){
    geneset[g,]$color <- "#99d8c9" #light green
  }
  
  if(geneset[g,]$parents_children == 1){
    geneset[g,]$color <- "#bcbddc" #lavender
  }
  
  if(geneset[g,]$atg7ko_children ==1){
    geneset[g,]$color <- "#fb6a4a" #salmon
  }
  
  if(geneset[g,]$atg7ko_parents ==1){
    geneset[g,]$color <- "#9ecae1" #light blue
  }
  
  
  
}

# extract final subgraph --------------------------------------------------
final_subgraph <- subgraph(net, nodes = geneset$gene)
write.dot(file = paste0(out_prefix, ".dot.txt"), graph = final_subgraph)
write.table(geneset, file = paste0(out_prefix, "_color_table.txt"), sep = ",")



# generate file for genie -------------------------------------------------
mliver_subset <- mliver_data[, c(which(colnames(mliver_data) %in% geneset$gene))]
write.csv(mliver_subset, file = paste0(out_prefix, "_data.csv"), row.names = FALSE, quote = FALSE)
# mliver.fit.bayes <- bn.fit(final_subgraph, data = mliver_subset, method = "bayes", iss = 0.1)
# write.net(file = paste0(out_prefix, "_bayes_iss=0.1.net"), fitted = mliver.fit.bayes)

mliver.fit.bayes <- bn.fit(final_subgraph, data = mliver_subset, method = "bayes", iss = 1)
write.net(file = paste0(out_prefix, "_bayes_iss=1.net"), fitted = mliver.fit.bayes)
# rm(list=ls())