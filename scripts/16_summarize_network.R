source("./scripts/src/summarizeNetworkInfluence.R")
source("./scripts/src/extractSubnetwork.R")
source("./scripts/src/utils_atg7.R")

dotfile <- "./analysis/banjo/output/nov17.mouse.liver.32h.ma.top.graph.2023.11.17.11.01.31.txt"
datafile <- "./analysis/banjo/input/nov17_mouse_liver_microarray_data.txt"
genefile <- "./data/GSE65174/select_genes_no_condition.txt"

out.prefix <- "./results/raw/summary_network_select_genes_no_condition"

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

dot <- paste0(out.prefix, ".dot")
write.dot(file = dot, subgraph)
dot_command <- paste0("dot -Tsvg ", dot, " -o ",  paste0(out.prefix, ".subnetwork.svg"))
system(dot_command)
cat("Final subnetwork graphic can be found in", paste0(out.prefix, ".subnetwork.svg"))
