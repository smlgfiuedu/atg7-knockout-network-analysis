source("./scripts/src/extractSubnetwork.R")
source("./scripts/src/utils_atg7.R")

dotfile <- "./results/raw/mouse.liver.32h.ma.tg.highlight-interestgenes.2023.08.02.17.23.44.txt"
out.prefix <- "TEST_mouse_liver_subset"

genelist <- read.table("./data/GSE65174/genelist.txt")[[1]]
mliver_data <- read.table("./data/GSE65174/mouse_liver_microarray_data.txt", sep = "\t", header = TRUE)

# preprocess for input to general functions -------------------------------------------------------
# add more descriptive experiment column name
colnames(mliver_data)[1] <- "Atg7KO"
mliver_data$Atg7KO <- factor(mliver_data$Atg7KO, labels = c("Atg7WT", "Atg7KO"))

# convert dot file to model string and bn object
str <- TransDottoString(dotfile)
net <- model2network(str)

# fix naming conventions
net <- fixDotNames(net)

# extract subnetwork ------------------------------------------------------
dot <- extractSubnetwork(network.data = mliver_data, net = net, target.node = "Atg7KO", file.out.prefix = out.prefix, assign.colors = TRUE)


# write dot graphic -------------------------------------------------------
dot_command <- paste0("dot -Tsvg ", dot, " -o ",  paste0(out.prefix, ".subnetwork.svg"))
system(dot_command)



