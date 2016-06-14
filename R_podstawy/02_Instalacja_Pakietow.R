#' #  Instalacja dodatkowych pakietów
#' 
#' Czasem niezbędne jest wykorzystanie funkcji i funkcjonalności z pakietów (bibliotek),
#' które nie są dystrybuowane z podstawową wersją R (base R). Aby zainstalować dodatkowe
#' pakiety z najpopularniejszego repozytorium pakietów CRAN (https://cran.r-project.org/)
#' poleceniem jak poniżej
#' 
## ---- eval=FALSE---------------------------------------------------------
## install.packages(c("dplyr", "knitr"))
## # pobiera pakietow o nazwie 'dplyr' i 'knitr'

#' 
#' Można też pobrać developerskie wersje pakietów, o ile istnieją, bezpośrednio z repozytorium na githubie
#' 
## ---- eval=FALSE---------------------------------------------------------
## # wczytanie potrzebnej biblioteki devtools zawierajacej funkcje install_github
## library(devtools) # jezeli nie posiadamy pakietu devtools, nalezy go wpierw zainstalowac z CRAN
## install_github('hadley/dplyr')
## install_github('yihui/knitr')

#' 
#' 
#' Istnieją także pakiety do R o zastosowaniach biologicznych i bioinformatycznych, które przyjęło się publikować
#' i pobierać z repozytorium pakietów BioConductor [https://www.bioconductor.org/](https://www.bioconductor.org/).
#' Przykładowa instalacja pakietu
## ---- eval=FALSE---------------------------------------------------------
## source("http://bioconductor.org/biocLite.R")
## biocLite(c("GenomicFeatures", "AnnotationDbi"))

#' 
#' 
