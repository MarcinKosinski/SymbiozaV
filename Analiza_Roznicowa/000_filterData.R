# 0. 
# - Wyczyścić dane: 
# Filer the data of the low variability , 
# by removing every column that has 0 in each sample,
# removing every row that has 0 in 70% in samples,
# proceeding with Median Absolute Deviation (MAD) analysis for each row,
# (removing 10% bottom retroelements)


# 0a) Filer the data of the low variability retroelements, 
# by removing every column that has 0 in each sample,

which(colSums(BRCA.rnaseq[,-1]) == 0)


# 0b) removing every row that has 0 in 70% in samples,

(rowSums(BRCA.rnaseq[,-1]) == 0) -> occure0_inROW

ROW2remove  <- which(occure0_inROW/nrow(BRCA.rnaseq) >= 0.7)

BRCA.rnaseq[-ROW2remove, ] -> BRCA.rnaseq.2


# 0c) proceeding with Median Absolute Deviation (MAD) analysis for each row,

mad_BRCA.rnaseq.2 <- apply(BRCA.rnaseq.2[,-1], 1, mad)
quantile(mad_BRCA.rnaseq.2, probs = seq(0,1,0.01)) 

# (removing 10% bottom retroelements)

BRCA.rnaseq.2[which(mad_BRCA.rnaseq.2 > 0), ] -> BRCA.rnaseq.3

write.csv(BRCA.rnaseq.3,
          file = "mad_BRCA.rnaseq.2.csv",
          quote = FALSE,
          row.names = FALSE)
