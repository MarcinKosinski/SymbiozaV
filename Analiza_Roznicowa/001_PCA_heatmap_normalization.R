# 1) plots
# 1.
# - a) PCA plot
# - b) Tabela znormalizowanych wartości (edgeR - normalizacja)
# - c) heatmapa na znormalizowanych wartościach (heatmap'a - heatmap package)

# 1 a)

# dodac grupy
library(RTCGA)
library(dplyr)
BRCA.rnaseq %>%
   mutate(grupa = ifelse(`MDM2|4193` >= median(`MDM2|4193`), "H", "L")) %>%
   select(grupa, bcr_patient_barcode, everything()) -> BRCA.rnaseq.groups
pcaTCGA(BRCA.rnaseq.groups[, -2], group.names = "grupa")

BRCA.rnaseq %>%
   mutate(grupa = cut(`MDM2|4193`,
                      breaks = quantile(BRCA.rnaseq$`MDM2|4193`, seq(0,1,0.25)),
                      include.lowest = TRUE, dig =5)) %>%
   select(grupa, bcr_patient_barcode, everything()) -> BRCA.rnaseq.groups.2
pcaTCGA(BRCA.rnaseq.groups.2[, -2], group.names = "grupa")



expressionsTCGA(OV.rnaseq, BRCA.rnaseq) -> OV_BRCA.rnaseq
OV_BRCA.rnaseq %>%
   mutate(grupa = cut(`MDM2|4193`,
                      breaks = quantile(OV_BRCA.rnaseq$`MDM2|4193`, seq(0,1,0.5)),
                      include.lowest = TRUE, dig =5)) %>%
   mutate(grupa = paste0(grupa, dataset)) %>%
   select(grupa, bcr_patient_barcode, dataset, everything()) -> BRCA.rnaseq.groups.2
pcaTCGA(BRCA.rnaseq.groups.2[, -c(2:3)], group.names = "grupa")




# 1 b) DZIALA NA RAW COUNTS
# - Tabela znormalizowanych wartości (edgeR - normalizacja)
library(edgeR)

BRCA.rnaseq.3_t <- t(BRCA.rnaseq.3)

rownames(PCBC_raw_counts) <- BRCA.rnaseq.3[, 1]
BRCA.rnaseq.4 <- BRCA.rnaseq.3[, -1]
y <- DGEList(BRCA.rnaseq.4)
y <- calcNormFactors(BRCA.rnaseq.4)
BRCA.rnaseq.4_normalize <- cpm(y, normalized.lib.sizes=FALSE)


# 1 c) heatmapa

library(dplyr)
library(pheatmap)
library(RColorBrewer)


# normalnie - featury w kolumnach

# annotacje wierszy i kolumn
# annotation_col <- as.data.frame(grupa)
# rownames(annotation_col) <- rownames(heatmap_data)
# colnames(annotation_col) <- "Category"
# Category <- 
# names(Category) <-




pheatmap(log1p(dane),
         annotation_col = annotation_col,
         cluster_row = F,
         show_colnames = T,
         filename = "heatmap.png",
         fontsize_row = 6, cellheight = 6, 
         fontsize=9, cluster_col=F,
         annotation_colors = list(Category = Category),
         color = rev(brewer.pal(6, "RdYlGn")))

