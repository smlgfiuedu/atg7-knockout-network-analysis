data(learning.test)
cpdag = pc.stable(learning.test)
dag = set.arc(cpdag, "A", "B")
graphviz.plot(dag)

network <- dag
subset <- c("A", "F", "D", "E")

# using subgraph() to extract the subnetwork, the connection A->E is lost because B is not a part of the subset, and bridges the connection: A -> B -> E
check <- subgraph(network, subset)
graphviz.plot(check)


# Need a way to demonstrate the connection A -> E without the presence of B


makeNodePairs <- function(nodes){
  pairs.table <- t(combn(nodes,2))
  colnames(pairs.table) = c("node1", "node2")
  # pairs.table <- cbind(pairs.table, pair=paste(pairs.table[,1], "_", pairs.table[,2], sep=''))
  return(pairs.table)
}


node.pairs <- makeNodePairs(subset)
subset.graph <- empty.graph(subset)
graphviz.plot(subset.graph) #demonstrate no connections

for(p in 1:nrow(node.pairs)){
  A <- node.pairs[p,1]
  B <- node.pairs[p,2]
  cat("Node 1:", A, "\t")
  cat("Node 2:", B, "\n")
  
  
  isPath <- path.exists(network, from = A, to = B)
  if(isPath){
    cat("Existing path between", A, "and", B, "\n")
    # print(incoming.arcs(network, A))
    # print(incoming.arcs(network, B))
    
    # check direction of connection between A and B
    AtoB <- A %in% ancestors(network, B)
    BtoA <- B %in% ancestors(network, A)  
    if(AtoB){
      cat(A, "->", B, "\n")
      subset.graph <- set.arc(subset.graph, from = A, to = B)
    }else if(BtoA){
        cat(A, "<-", B, "\n")
      subset.graph <- set.arc(subset.graph, from = B, to = A)
      }
    }
}
graphviz.plot(subset.graph)

intersect(descendants(network, A), ancestors(network, B)) # find path between A and B A->B
intersect(ancestors(network, A), descendants(network, B)) # find path between A and B when B -> A
network <- set.arc(network, from = "D", to = "E")
graphviz.plot(network)



# test function -----------------------------------------------------------
rm(list=ls())
source("./scripts/src/summarizeNetworkInfluence.R")

data(learning.test)
cpdag = pc.stable(learning.test)
dag = set.arc(cpdag, "A", "B")
graphviz.plot(dag)

subgraph <- summarizeNetworkInfluence(network = dag, subset = c("A", "F", "D", "E"))
graphviz.plot(subgraph)




# test function on atg7 data ----------------------------------------------
rm(list=ls())
source("./scripts/src/summarizeNetworkInfluence.R")
source("./scripts/src/extractSubnetwork.R")
source("./scripts/src/utils_atg7.R")

dotfile <- "./analysis/banjo/output/nov17.mouse.liver.32h.ma.top.graph.2023.11.17.11.01.31.txt"
datafile <- "./analysis/banjo/input/nov17_mouse_liver_microarray_data.txt"
genefile <- "./data/GSE65174/select_genes_with_condition.txt"

genelist <- read.table(genefile)[[1]]
mliver.data <- read.table(datafile, sep = "\t", header = TRUE)

# preprocess for input to general functions
# add more descriptive experiment column name
colnames(mliver.data)[1] <- "Atg7KO"
mliver.data$Atg7KO <- factor(mliver.data$Atg7KO, labels = c("Atg7WT", "Atg7KO"))

# convert dot file to model string and bn object
str <- TransDottoString(dotfile)
net <- model2network(str)

# fix naming conventions
net <- fixDotNames(net)

## SUMMARIZE STEP
subgraph <- summarizeNetworkInfluence(network = net, subset = genelist)

dot <- "test.dot"
out.prefix <- "test.summary"
write.dot(file = dot, subgraph)
dot_command <- paste0("dot -Tsvg ", dot, " -o ",  paste0(out.prefix, ".subnetwork.svg"))
system(dot_command)
cat("Final subnetwork graphic can be found in", paste0(out.prefix, ".subnetwork.svg"))
