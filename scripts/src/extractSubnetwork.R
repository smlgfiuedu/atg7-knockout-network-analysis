source("./scripts/src/utils_network.R")
require("bnlearn")
require("tidyverse")

# main extract subnetwork function ----------------------------------------
# returns the location of the final product dotfile
extractSubnetwork <- function(network.data, net, target.node, method = "mle", iss = .1, assign.colors = FALSE, generate.net.file = FALSE, file.out.prefix  = getwd(), genelist = list(), useConnectGenes = TRUE)
{
  # recodify and factor column values
  for (l in 2:ncol(network.data)){
    network.data[,l] <- .relabelGeneVector(network.data[,l])
  }
  
  # fit entire model
  full.fit <- .fitFullModelBN(net, network.data)

  if(useConnectGenes){
    # markov blanket
    markov.blanket <- mb(net, node = target.node)
    
    # extract closely related nodes (within 3 lineages)
    connected.genes <- .findCloseNodes(net, markov.blanket, genelist)
    
    
    geneset <- .createGeneList(net, target.node, markov.blanket, connected.genes)    
  }else{
    geneset <- data.frame(gene = genelist[genelist %in% names(net$nodes)],
                          parents = 0,
                          children = 0,
                          parents_children = 0,
                          genes_interest = 0,
                          markov_blanket = 0)
  
  }

  
  # extract final subgraph
  message("Extracting final subgraph and writing to dot file")
  subset.dotfile <- paste0(file.out.prefix, ".dot.txt")
  final.subgraph <- subgraph(net, nodes = geneset$gene)
  write.dot(file = subset.dotfile, graph = final.subgraph)
  
  #generate genie file
  if(generate.net.file){
    message("generating .net file for GeNIe")
    network.subset <- network.data[, c(which(colnames(network.data) %in% geneset$gene))]
    
    subgraph.fit <- bn.fit(final.subgraph, data = mliver_subset, method = method, iss = iss)
    write.net(file = paste0(file.out.prefix, "_", method, "_iss=", iss, ".net"), fitted = subgraph.fit)
    
  }
  
  
  annotated.out <- paste0(file.out.prefix, "color.annotated.dot.txt")
  if(assign.colors == TRUE){
    # assign colors to the nodes for use in dot graph
    message("Assigning colors for dot graph")
    geneset <- .assignColors(geneset, target.node)
    
    # run code for rewriting dot file to include colors
    message("Writing annotated dot file")
    .annotateNet(filepathI = subset.dotfile, filepathO = annotated.out, 
                nodes = list(genes = geneset$gene,
                             colors = geneset$color),
                target.node = target.node)
    return(annotated.out)
  }else{
    return(subset.dotfile)
  }
}