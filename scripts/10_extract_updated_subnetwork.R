source("./scripts/src/extractSubnetwork.R")
source("./scripts/src/utils_atg7.R")

datafile <- "./analysis/banjo/input/updated_goi_mouse_liver_microarray_data.txt"
genefile <- "./data/GSE65174/updated_genes_of_interest.txt"
dotfile <- "./results/raw/updated_goi/updated.goi.mouse.liver.32h.ma.top.graph.2023.10.30.15.18.17.txt"
out.prefix <- "./results/raw/updated_goi/updated_goi_mouse_liver_subset"

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

# extract subnetwork ------------------------------------------------------
dot <- extractSubnetwork(network.data = mliver.data, net = net, target.node = "Atg7KO", 
                         file.out.prefix = out.prefix, assign.colors = TRUE)


# write dot graphic -------------------------------------------------------
dot_command <- paste0("dot -Tsvg ", dot, " -o ",  paste0(out.prefix, ".subnetwork.svg"))
system(dot_command)
cat("Final subnetwork graphic can be found in", paste0(out.prefix, ".subnetwork.svg"))
