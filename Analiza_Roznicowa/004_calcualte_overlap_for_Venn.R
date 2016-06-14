#

library(VennDiagram)
library(dplyr)
calculate.overlap(list(file1,
                       file2,
                       file3
                      )) -> overlap

overlap



x <- list(file1,
          file2,
          file3)

get.venn.partitions(x) -> y
get.venn.partitions(x, force.unique = FALSE)

