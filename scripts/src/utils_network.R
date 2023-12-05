require("bnlearn")
require("tidyverse")
# create structure from dot file for bnlearn package ----------------------
# adapted from Zhengua Gong -- FIU
.processFile = function(filepath) {
  rel = vector()
  ren = vector()
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    if ( grepl("\\[|\\]", line) ) {
      ren = append(ren, line)
      next
    }
    
    if ( grepl("->", line) ) {
      rel = append(rel, line)
      next
    }
    
    else {
      next
    }
    
  }
  
  close(con)
  
  resu = list("name" = ren, "str" = rel)
  return (resu)
}
TransDottoString = function(dotfile){
  node = vector()
  nn = vector()
  nc = vector()
  re = .processFile(dotfile)
  
  na = re$name
  for(j in 1:length(na)){
    nalist <- na[j] %>%
      gsub("\\s", "", .) %>%
      gsub(";", "", .) %>%
      strsplit(., "\\[|label=|\"|\\]") 
    nn = append(nn, nalist[[1]][1])
    nc = append(nc, nalist[[1]][4])
  }
  
  rs = re$str
  for (i in 1:length(rs)){
    stlist<- rs[i] %>%
      gsub("\\s", "", .) %>%
      gsub(";", "", .) %>%
      strsplit(., "->") 
    node = append(node, stlist[[1]])
  }
  
  #   dag = empty.graph(nodes = unique(node))
  #   arc.set <- matrix(node, byrow = TRUE, ncol = 2, dimnames = list(NULL, c("from", "to")))
  #   arcs(dag) <- arc.set
  #   S = modelstring(dag)
  node2 = as.character(factor(node, levels = nn, labels = nc))
  dag2 = empty.graph(nodes = unique(nc))
  arc.set2 <- matrix(node2, byrow = TRUE, ncol = 2, dimnames = list(NULL, c("from", "to")))
  arcs(dag2) <- arc.set2
  S2 = modelstring(dag2)
  
  return(S2)
}

# decode a gene expression vector into a character factor -----------------
.relabelGeneVector <- function(gene_vector){
  
  is_none <- gene_vector == 0
  is_low <- gene_vector == 1
  is_high <- gene_vector == 2
  
  
  gene_vector[is_none] <- "No Change"
  gene_vector[is_low] <- "Low"
  gene_vector[is_high] <- "High"
  
  gene_vector <- factor(gene_vector, levels = c("No Change", "Low", "High"))
  return(gene_vector)
  
}


# fit full model using bnlean ---------------------------------------------
.fitFullModelBN <- function(net, network.data, method = "mle", iss = 0.1){
  fit.message <- paste("Fitting model using method =", method, "and iss = ", iss)
  message(fit.message)
  
  ordered <- match(names(net$nodes), colnames(network.data))
  if(any(is.na(ordered))){
    stop("node names do not match matrix columns - check naming conventions or if there are any missing columns/nodes")
  }
  network.data <- network.data[ ,ordered]
  all(names(net$nodes)==colnames(network.data))
  
  full.fit <- bn.fit(net, data = network.data, method = method, iss = iss)
  
  bde.message <- paste("Model score:", score(net, data = network.data))
  message(bde.message)
  
  return(full.fit)
}

# finds and returns a list of nodes within 3 lineages of the node of interest  --------
.findCloseNodes <- function(net, markov.blanket, genelist){
  message("Finding nodes within 3 generation of the target node")
  # parents through great grandparents --------------------------------------
  tier1 <- sapply(markov.blanket, FUN = function(x) parents(net, x)) %>% unlist() %>% unique()
  tier2 <- sapply(tier1, FUN = function(x) parents(net, x)) %>% unlist() %>% unique() 
  tier3 <- sapply(tier2, FUN = function(x) parents(net, x)) %>% unlist() %>% unique() 
  
  # check which genes of interest are within 3 generations of the mk --------
  first_gen_goi <- intersect(genelist, tier1)
  second_gen_goi <- intersect(genelist, tier2)
  third_gen_goi <- intersect(genelist, tier3)
  
  parent_goi_to_include <- unique(c(first_gen_goi, second_gen_goi, third_gen_goi)) 
  
  
  # for second gen, find which children are parents of mkvs
  second_gen_children <- sapply(second_gen_goi, FUN = function(x) children(net, x)) %>% unlist() %>% unique()
  gen2c_link <- intersect(second_gen_children, tier1)
  
  # for the third gen, find which children are the parents of  the mkv parents
  third_gen_children <- sapply(third_gen_goi, FUN = function(x) children(net, x)) %>% unlist() %>% unique()
  gen3c_link <- intersect(third_gen_children, tier2)
  
  
  parent_goi_link <- union(gen2c_link, gen3c_link)
  
  
  # children through great grandchildren ------------------------------------
  tier1 <- sapply(markov.blanket, FUN = function(x) children(net, x)) %>% unlist() %>% unique()
  tier2 <- sapply(tier1, FUN = function(x) children(net, x)) %>% unlist() %>% unique() 
  tier3 <- sapply(tier2, FUN = function(x) children(net, x)) %>% unlist() %>% unique() 
  
  first_gen_goi <- intersect(genelist, tier1)
  second_gen_goi <- intersect(genelist, tier2)
  third_gen_goi <- intersect(genelist, tier3)
  
  child_goi_to_include <- unique(c(first_gen_goi, second_gen_goi, third_gen_goi)) 
  
  # for second gen, find which parents are children of mkvs
  second_gen_parents <- sapply(second_gen_goi, FUN = function(x) parents(net, x)) %>% unlist() %>% unique()
  gen2p_link <- intersect(second_gen_parents, tier1)
  
  
  # for the third gen, find which parents are the children of  the mkv children
  third_gen_parents <- sapply(third_gen_goi, FUN = function(x) parents(net, x)) %>% unlist() %>% unique()
  gen3p_link <- intersect(third_gen_parents, tier2)
  
  children_goi_link <- union(gen2p_link, gen3p_link)
  
  generations <- list(parent_goi_to_include = parent_goi_to_include,
                      child_goi_to_include = child_goi_to_include,
                      parent_goi_link = parent_goi_link,
                      children_goi_link = children_goi_link)
  return(generations)
  
}


# create genelist dataframe -----------------------------------------------
.createGeneList <- function(net, target.node, markov.blanket, connected.genes){
  message("Creating gene list")
  pc <- which(names(net$nodes)==net$nodes[[target.node]]$children)
  if(length(pc)>0){
    parents_children = net$nodes[[pc]]$parents
  }else{
    parents_children = NA
  }
  
  all.genes <- list(exp = target.node,
                    exp_parents = net$nodes[[target.node]]$parents, #blue
                    exp_children = net$nodes[[target.node]]$children, #red
                    parents_children = parents_children, #purple
                    genes_interest = union(connected.genes$parent_goi_to_include, connected.genes$child_goi_to_include), #green
                    markov_blanket = markov.blanket,
                    connect_genes = union(connected.genes$parent_goi_link, connected.genes$children_goi_link))
  
  genes <- unlist(all.genes) 
  names(genes) <- gsub(names(genes), pattern = "[0-9]", replacement = "")
  
  geneset <- data.frame(gene = unique(genes),
                        parents = 0,
                        children = 0,
                        parents_children = 0,
                        genes_interest = 0,
                        markov_blanket = 0)
  
  na <- which(is.na(geneset))
  if(length(na)>0){
    geneset <- geneset[-na,]
  }
  
  for(g in 1:nrow(geneset)){
    types <- names(genes)[which(genes==geneset$gene[g])]
    
    if("exp_parents" %in% types){
      geneset[g, "parents"] <- 1
    }
    
    if("exp_children" %in% types){
      geneset[g, "children"] <- 1
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
  
  n <- length(geneset$gene[which(geneset$genes_interest == 1)])
  goi.message <- paste(n, "genes of interest are included in the subnetwork")
  message(goi.message)
  return(geneset)
}


# assign colors for use in dot graphic ------------------------------------
.assignColors <- function(geneset, target.node){
  geneset$color <- "white"
  
  for(g in 1:nrow(geneset)){
    gene <- geneset$gene[g]
    if(gene == target.node){
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
    
    if(geneset[g,]$children ==1){
      geneset[g,]$color <- "#fb6a4a" #salmon
    }
    
    if(geneset[g,]$parents ==1){
      geneset[g,]$color <- "#9ecae1" #light blue
    }
  }
  return(geneset)
}



# add color-coded nodes to dot file ---------------------------------------
.annotateNet <- function(filepathI, filepathO, nodes = list(), target.node) {
  
  genes <- nodes$genes
  fr <- file(filepathI, open="rt") #open file connection to read
  fw <- file(filepathO, open="wt") #open file connection to write 
  
  
  lines <- readLines(fr)
  
  edge_start <- "edge "
  graph_start <- "digraph {"
  graph_end <- "}"
  
  for (line in lines){
    # skip lines that direct graph structure
    isLabelLine <- !(grepl(edge_start, line, fixed = TRUE) | grepl(graph_start, line, fixed = TRUE) | grepl(graph_end, line, fixed = TRUE))
    if (isLabelLine){
      gene_interest <- which(sapply(genes, FUN = function(x) grepl(paste0("\\b", x, "\\b"), line)))
      color <- nodes$colors[gene_interest]
      
      line <- paste0("\"", genes[gene_interest], "\" ", "[label=\"", genes[gene_interest], "\", style=\"filled\", fillcolor=\"", color, "\", color=\"black\"];")
      # check if the node is the condition node
      if(grepl(target.node, line, fixed = TRUE)){
        line <- gsub("];", replacement = ", fontname=\"Calibri-BoldItalic\"];", x = line)
        writeLines(line, fw)
      }else{
        writeLines(line, fw)
      }
    }else{
      writeLines(line,fw) #write the rest of the graph
    }
  }
  close(fr);close(fw); #close connections
}

