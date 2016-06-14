library(RTCGA.rnaseq)
library(dplyr)



# pwartosc <- numeric(length(1:20533))
# mean_Yes <- numeric(length(3:20533))
# mean_No <- numeric(length(3:20533))
# median_Yes <- numeric(length(3:20533))
# median_No <- numeric(length(3:20533))
# 
# for(i in 3:20533){
#    kruskal.test(x = d2analyze[, i], g = as.factor(d2analyze$status))$p.value -> pwartosc[i-2]
#    mean(d2analyze[as.factor(d2analyze$status) == 'Yes', i]) -> mean_Yes[i-2]
#    mean(d2analyze[as.factor(d2analyze$status) == 'No', i]) -> mean_No[i-2]
#    median(d2analyze[as.factor(d2analyze$status) == 'Yes', i]) -> median_Yes[i-2]
#    median(d2analyze[as.factor(d2analyze$status) == 'No', i]) -> median_No[i-2]
#    
#    cat(i, "\r")
# }
# 
# head(names(d2analyze))
# tail(names(d2analyze))
# 
# 
# names(d2analyze)[3:20533][which(!is.na(pwartosc))] -> l
# 
# 
# mean_Yes[!is.na(pwartosc)] -> mean_Yes
# mean_No[!is.na(pwartosc)] -> mean_No
# median_No[!is.na(pwartosc)] -> median_No
# median_Yes[!is.na(pwartosc)] -> median_Yes
# 
# na.omit(pwartosc) -> pwartosc
# p.adjust(pwartosc, method = "BH") -> adj_p
# 
# data.frame(gen = l,
#            adj_p = adj_p,
#            mean_Yes = mean_Yes,
#            mean_No = mean_No, 
#            median_No = median_No,
#            median_Yes = median_Yes) -> results
# 
# results %>%
#    arrange(adj_p) -> res_arranged
# 
# 
# write.table(res_arranged, 
#             file = "nazwa_pliku.txt",
#             sep = "\t",
#             quote = FALSE,
#             col.names = TRUE,
#             row.names = FALSE)
