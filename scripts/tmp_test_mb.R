source("./scripts/src/extractSubnetwork.R")
source("./scripts/src/utils_atg7.R")

datafile <- "./analysis/banjo/input/nov17_mouse_liver_microarray_data.txt"
genefile <- "./data/GSE65174/select_genes_no_condition.txt"
dotfile <- "./analysis/banjo/output/nov17.mouse.liver.32h.ma.top.graph.2023.11.17.11.01.31.txt"
out.prefix <- "./results/raw/updated_goi/TEST_select_genes_no_condition_mouse_liver_subset"

genelist <- read.table(genefile)[[1]]
mliver.data <- read.table(datafile, sep = "\t", header = TRUE)

# preprocess for input to general functions -------------------------------------------------------
# add more descriptive experiment column name
colnames(mliver.data)[1] <- "Atg7KO"
mliver.data$Atg7KO <- factor(mliver.data$Atg7KO, labels = c("Atg7WT", "Atg7KO"))

# convert dot file to model string and bn object
str <- TransDottoString(dotfile)
net <- model2network(str)

# fix naming conventions
net <- fixDotNames(net)


## TEST
network.data <- mliver.data
# recodify and factor column values
for (l in 2:ncol(network.data)){
  network.data[,l] <- .relabelGeneVector(network.data[,l])
}

# fit entire model
full.fit <- .fitFullModelBN(net, network.data)


data(learning.test)
cpdag = pc.stable(learning.test)
dag = set.arc(cpdag, "A", "B")
graphviz.plot(dag)

ancestors(dag, "E")
tn.connections <- list()

tn.connections[["ancestors"]] <- ancestors(net, "Atg7KO")
tn.connections[["descendants"]] <- descendants(net, "Atg7KO")
