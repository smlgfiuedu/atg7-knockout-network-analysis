source("./scripts/src/utils_dotfile.R")


# original atg7 knockout --------------------------------------------------
highlightVariable(filepathT = "./analysis/banjo/output/mouse.liver.1h.ma.top.graph.2023.08.02.17.34.22.txt",
              filepathO = "./results/data/mouse.liver.1h.ma.tg.highlight.2023.08.02.17.34.22.txt",
              interest_variable = "Experiment")


highlightVariable(filepathT = "./analysis/banjo/output/mouse.liver.4h.ma.top.graph.2023.08.02.17.27.01.txt",
                  filepathO = "./results/data/mouse.liver.4h.ma.tg.highlight.2023.08.02.17.27.01.txt",
                  interest_variable = "Experiment")

highlightVariable(filepathT = "./analysis/banjo/output/mouse.liver.16h.ma.top.graph.2023.08.02.17.29.02.txt",
                  filepathO = "./results/data/mouse.liver.16h.ma.tg.highlight.2023.08.02.17.29.02.txt",
                  interest_variable = "Experiment")

highlightVariable(filepathT = "./analysis/banjo/output/mouse.liver.32h.ma.top.graph.2023.08.02.17.23.44.txt",
                  filepathO = "./results/data/mouse.liver.32h.ma.tg.highlight.2023.08.02.17.23.44.txt",
                  interest_variable = "Experiment")



# updated atg7 knockout ---------------------------------------------------
highlighted.output <- "./results/raw/updated_goi/updated.goi.mouse.liver.32h.ma.top.graph.2023.10.30.15.18.17.txt"

highlightVariable(filepathT = "./analysis/banjo/output/updated.goi.mouse.liver.32h.ma.top.graph.2023.10.30.15.18.17.txt",
                  filepathO = highlighted.output,
                  interest_variable = "Experiment")

dot_command <- paste0("dot -Tsvg ", highlighted.output, " -o ",  "./results/figures/updated_goi/updated.goi.atg7ko.mouseliver.subnetwork.svg")
system(dot_command)
