library(tidyverse)

dotfile <- "./results/raw/updated_goi/updated_goi_mouse_liver_subset.dot.txt"
colors <- "./results/raw/updated_goi/updated_goi_mouse_liver_subset_color_table.txt"
out <- "./results/raw/updated_goi/updated_goi_mouse_liver_annotated.txt"
graph <- "./results/figures/updated_goi/updated_goi_mouse_liver_subnetwork.svg"

annotateNet <- function(filepathI, filepathO, nodes = list(), condition_node) {
  
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
      if(grepl(condition_node, line, fixed = TRUE)){
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


var_colors <- read.table(colors, sep = ",")
annotateNet(filepathI = dotfile, filepathO = out, 
            nodes = list(genes = var_colors$gene,
                         colors = var_colors$color),
            condition_node = "Atg7KO")


dot_command <- paste0("dot -Tsvg ", out, " -o ",  graph)
system(dot_command)
