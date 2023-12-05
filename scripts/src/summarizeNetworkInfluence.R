# functions to summarize the network based on a specified set of genes
library(bnlearn)

.makeNodePairs <- function(nodes){
  pairs.table <- t(combn(nodes,2))
  colnames(pairs.table) = c("node1", "node2")
  return(pairs.table)
}


summarizeNetworkInfluence <- function(network, subset){
  node.pairs <- .makeNodePairs(subset)
  subset.graph <- empty.graph(subset)
  
  
  for(p in 1:nrow(node.pairs)){
    A <- node.pairs[p,1]
    B <- node.pairs[p,2]

    isPath <- path.exists(network, from = A, to = B)
    if(isPath){
      cat("Existing path between", A, "and", B, "\n")
      
      # check direction of connection between A and B
      AtoB <- A %in% ancestors(network, B)
      BtoA <- B %in% ancestors(network, A)  
      if(AtoB){
        cat("direction:", A, "->", B, "\n")
        subset.graph <- set.arc(subset.graph, from = A, to = B)
      }else if(BtoA){
        cat("direction:", A, "<-", B, "\n")
        subset.graph <- set.arc(subset.graph, from = B, to = A)
      }
    }
  }
  return(subset.graph)
}