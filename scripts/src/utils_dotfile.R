# function to highlight the experimental variable of interest in a dot file --------
highlightVariable <- function(filepathI, filepathO, interest_variable) {
  require("tidyverse")  
  fr <- file(filepathI, open="rt") #open file connection to read
  fw <- file(filepathO, open="wt") #open file connection to write 
  
  lines <- readLines(fr)
  
  for (line in lines){
    # skip lines that direct graph structure
    if (grepl(interest_variable, line, fixed = TRUE)){
      line <- gsub("];", replacement = ", style=filled, fillcolor=\"yellow\"];", x = line)
      writeLines(line,fw)
    }else{
      writeLines(line,fw)        
    }
  }
  close(fr);close(fw) #close connections
}