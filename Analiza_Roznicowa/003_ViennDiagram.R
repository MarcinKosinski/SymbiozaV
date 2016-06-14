# vienn diagram

file1
file2
file3


library(VennDiagram)
library(dplyr)
venn.diagram(list(file1,
                  file2,
                  file3 
                  ),
             imagetype = "png",
             col = "transparent",
             filename = "ViennDiagram.png",
             height = 2000,
             width = 2000,
             fill = 1:3)
             
