# relabel experiment to specify condition & transcripts that begin with numbers so R will stop crying and match with dataframe
fixDotNames <- function(net){
  exp_node <- which(names(net$nodes)=="Experiment")
  number_nodes <- grep("^[0-9]", names(net$nodes))
  h2_node <- grep("-", names(net$nodes))
  
  new_labels <- names(net$nodes)
  new_labels[exp_node] <- "Atg7KO"
  new_labels[number_nodes] <- paste0("X", names(net$nodes)[number_nodes])
  new_labels[h2_node] <- gsub("-", ".", names(net$nodes)[h2_node])
  net <- rename.nodes(net, new_labels)
  
  return(net)
  
}