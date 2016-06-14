# 01 Podstawy R
Marcin Kosiński  
Ostatnia data modyfikacji `r Sys.Date()`  

#  Instalacja dodatkowych pakietów

Czasem niezbędne jest wykorzystanie funkcji i funkcjonalności z pakietów (bibliotek),
które nie są dystrybuowane z podstawową wersją R (base R). Aby zainstalować dodatkowe
pakiety z najpopularniejszego repozytorium pakietów CRAN (https://cran.r-project.org/)
poleceniem jak poniżej


```r
install.packages(c("dplyr", "knitr"))
# pobiera pakietow o nazwie 'dplyr' i 'knitr'
```

Można też pobrać developerskie wersje pakietów, o ile istnieją, bezpośrednio z repozytorium na githubie


```r
# wczytanie potrzebnej biblioteki devtools zawierajacej funkcje install_github
library(devtools) # jezeli nie posiadamy pakietu devtools, nalezy go wpierw zainstalowac z CRAN
install_github('hadley/dplyr')
install_github('yihui/knitr')
```


Istnieją także pakiety do R o zastosowaniach biologicznych i bioinformatycznych, które przyjęło się publikować
i pobierać z repozytorium pakietów BioConductor [https://www.bioconductor.org/](https://www.bioconductor.org/).
Przykładowa instalacja pakietu

```r
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicFeatures", "AnnotationDbi"))
```





# Materiały oparte o

- [Przewodnik po pakiecie R, P. Biecek](http://biecek.pl/R/)
- [Analiza danych z programem R, P. Biecek](http://biecek.pl/Analiza.Danych/)
- [Programowanie w języku R, M. Gągolewski](http://rksiazka.rexamine.com/)
- Darmowy Kurs MOOC Analizy i Przetwarzania Danych w R - [Pogromcy Danych](http://pogromcydanych.icm.edu.pl/)
- [In-depth introduction to machine learning in 15 hours of expert videos](http://www.r-bloggers.com/in-depth-introduction-to-machine-learning-in-15-hours-of-expert-videos/)
- [Microsoft Launches Its First Free Online R Course on edX](http://www.r-bloggers.com/microsoft-launches-its-first-free-online-r-course-on-edx/)
