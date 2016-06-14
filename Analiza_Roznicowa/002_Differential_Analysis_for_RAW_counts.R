#
# 2.
# - Differential Expression analysis using edgR (0,01%FDR): 

# 2a)  Differential Expression analysis using edgR (0,01%FDR): see annotation file

y <- DGEList(counts = dataMatrix, group = groups)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y, dispersion = "common")
resultMatrix <- topTags(et, n=Inf)$table

write.csv2(resultMatrix,
           file ="wynik.csv",
           quote = FALSE,
           row.names = TRUE)
