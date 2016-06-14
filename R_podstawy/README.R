#' ---
#' title: "01 Podstawy R"
#' author: "Marcin Kosiński"
#' date: "Ostatnia data modyfikacji `r Sys.Date()`"
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 2
#'     theme: readable
#'     keep_md: true
#'     number_sections: true
#' ---
#' 
#' 
## ---- echo=FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(comment="", message=FALSE, warning = FALSE, 
               tidy.opts=list(keep.blank.line=TRUE, width.cutoff=150),
							 options(width=150))

#' 
#' 
## ----child='01_Podstawy.Rmd'---------------------------------------------

#' # Podstawy R
#' 
#' ## Podstawowe typy zmiennych
#' 
#' Nie ma skalarów, są wektory
## ------------------------------------------------------------------------
# sprawdzanie typów i klas obiektów
is.vector(5)
class(5)
typeof(5)

#' 
#' ### Typy atomowe
#' 
#' R należy do klasy języków funkcyjnych - każda wykonywana "operacja"jest de facto sprowadzana do wywołania pewnej funkcji. I tak wyrażenie
## ------------------------------------------------------------------------
# przypisanie do zmiennej u liczby 5
u <- 5

#' tak naprawdę intepretowane jest przez parser jako:
## ------------------------------------------------------------------------
'<-'('u',5)

#' co oznacza wywołaj funkcję `'<-'()` z argumentami `'u'` i `5` 
#' 
#' 
#' #### Wartości logiczne
#' 
## ------------------------------------------------------------------------
TRUE
FALSE
# nie trzeba wczesniej deklarowac typow
x <- TRUE
x
y <- c(TRUE, FALSE, TRUE)
y
z <- c(T, F, c(F, F)) # funkcja 'c' tworzy wektory
z
TRUE == FALSE
y == x 
!x # negacja
y | x # alternatywa
y & x # koniunkcja

#' 
#' #### Wektory liczbowe
#' 
## ------------------------------------------------------------------------
# nie trzeba wczesniej deklarowac typow
w1 <- 1:10
w1
# wybieranie elementów
# numerowanie od 1
w1[3:5] # od 3 do 5 elementu
w1[3:5] <- c(30,40,50) # nadpisanie od 3 do 5 elementu
w1[1:5]
w1[c(1,2,3,4,5)] # równoważne odwołanie
w1[c(1, 8)]# pierwszy i osmy element
w2 <- c(4, 2, 4)
class(w2)
typeof(w2)
# sequence/sekwencja
(d  <- seq(from=0, to=1, by=0.25))
## wyniki niepoprawnych odwolan
w1[0]
w1[15]
# dynamiczne rozszerzenie wektorow
w1[12] <- 100
w1
w1[-c(1:5)] # bez pierwszych pięciu

#' 
#' Operacje arytmetyczne
## ------------------------------------------------------------------------
2+4
x <- c(3, 9)
y <- c(5, 10)
x+y
z <- c(1, 2, 3, 4)
z+y # zawijanie wektorow
z*5 # wektoryzacja operacji
x > y

#' 
## ------------------------------------------------------------------------
# INTY
# nie trzeba wczesniej deklarowac typow
x <- c(4L, 2L)
class(x)
typeof(x)
is.integer(x)
is.numeric(x)

#' 
## ------------------------------------------------------------------------
# zespolone
1+3i + 2+6i

#' 
#' 
#' 
#' #### Wektory napisów
#' 
## ------------------------------------------------------------------------
typeof('Hello World!')
class("Hello world!")
length("Hello world!")
length(c("Hello", "world!"))
Hello <- c("Hello", "world!")
cat(Hello)
print(Hello)
str(Hello)

#' 
#' 
#' #### Braki danych
#' 
#' NA - Not available
## ------------------------------------------------------------------------
is.na(NA)
is.na(c(1,NA,3))
c(1,2,NA)*4
# nieskonczonosc
1/0
# nie-liczna
0/0
sqrt(-1)
log(0)

# typ pusty
NULL
is.null(NULL)
is.null(c())# brak informacji o typie

#' 
#' 
#' 
#' 
#' 
#' ## Listy
#' 
#' Ciąg złożony z elementów o dowolnych typach (a więc już niekoniecznie tych samych jak w przypadku wektorów atomowych).
#' 
## ------------------------------------------------------------------------
lista <- list(c(1:10), c("Hello", "World!"), T)
lista
lista[[2]]
lista[2]
lista[[1]][2:4] # odwolanie rekurencyjne
lista2 <- list(lista, (1:5)*2)
lista2

#' 
#' **Lista** jako abstrakcyjna strukturadanych może być niektórym znana z wykładów z algorytmiki. Taka
#' struktura umożliwia szybkie dodawanie nowych elementów(koszt stały, czyli *O*(1)), nieco więcej czasu (pesymistycznie *O*(n)), wymagając od operacji dostępu do konkretnych elementów.
#' 
#' Lista, jako że jest reprezntowana w tapmięci komputera w postaci *tablic*, nie ma powyższej własności. W takiej implementacji nacisk jest położony na szybkie wykonywanie najczęstszej operacji - dostęp do elementów listy ma tutaj stały koszt (*O*(1)), a rzadziej wykorzystywane powiększanie/skracanie listy jest bardziej czasochłonne (*O*(n)).
#' 
#' Osoby posiadające doświadczenie w programowaniu w języku C powinny więc postrzegać reprezentantów typu `list`  jako tablice zawierające wskaźniki do dowolnych obiektów, czyli `void*`.
#' 
#' 
#' ## Najpopularniejsze typy złożone
#' 
#' ### Macierze
#' 
## ------------------------------------------------------------------------
(x <- matrix(1:6, nrow = 2, ncol = 3))
matrix(c("a","b","c", "d"), ncol = 4) # regula zawijania
#  dolaczanie kolumn lub wierszy
cbind(x, c(7, 8)) # column bind
rbind(x, c(10, 11, 12)) # row bind
# odwolania do elementow
x[1, 1]
x[1, 1:3]
x[, 1:2] # wszystkie wiersze oraz 1sza i 2ga kolumna
# x[wiersze, kolumny]
x^2

#' 
#' ### Czynniki - factor
#' 
## ------------------------------------------------------------------------
(f <- factor(c("Warszawa", "Kraków", "Warszawa", "Gdańsk")))
# wektor napisow ze slownikiem
levels(f)
# umozliwia szybka podmiane poziomow
levels(f) <- c("Miasto1", "Miasto2", "Miasto3")
f

#' 
#' 
#' ### Ramki danych
#' 
#' Ramki danych (`data.frame`) to obiekty przypominające zbiory danych znane być może niektórym z programów Microsoft Excel, SPSS czy Statistica. Zazwyczaj wiersze odpowiadają o serwacjom, np. badanym pacjentom, zaś kolumny odpowiadają informacją na temat wartości różnych zmiennych. Nie jest to reguła.
#' 
#' Ramki danych są reprezentowane w R przez listy zawierające wektory atomowe o tej samej długości. Każdy element tej szczególnej listy odpowiada kolumnie ramki danych.
#' 
#' Przykład
#' 
## ------------------------------------------------------------------------
dane <- data.frame(
	plec=c("M", "K", "M"),
	pali=c(T,F,F),
	wiek=c(26,25,22),
	stringsAsFactors = FALSE # czy napisy przechowywać jako factory
)
dane
dim(dane)
colnames(dane)
# odwolania
dane[, 2]
dane[1, 1]
dane[3, 1:2]

#' 
#' 
#' ## Zwektoryzowane Operacje
#' 
#' Większość funkcji i operacji w R jest zwektoryzowana - działa dla każdego elementu osobno bądź zwraca globalny wynik dla całego wektora.
## ------------------------------------------------------------------------
sum(1:10)
cumsum(1:10)
(1:10)^2
matrix(1:10, ncol = 2, nrow = 5) -> macierz
macierz
macierz+macierz

#' 
#' ## (Pre)-Deklaracja zmiennych
#' 
#' Dla wygody pisania, nie zawsze potrzebna jest deklaracja zmiennych,
#' jest to czasem związane z tym, ze R wolniej dziala (gdyż dynamicznie rozszerzając obiekty kopiuje je od nowa i dopisuje nowy element)
## ------------------------------------------------------------------------
x <- numeric(10)
x
x[2] <- 14
x[2]
x[16] <- 5
x
st <- character(0)
ints <- integer(125)
logiczne <- logical(length = length(x))
sum(logiczne)

#' 
#' 
#' 
#' ## Funkcje
#' 
#' 
## ------------------------------------------------------------------------
przykladowa_funkcja_pomnoz_i_podnies_do_potegi <- function(arg1, arg2, power){
	stopifnot(is.numeric(arg1))
	stopifnot(is.numeric(arg2))
	stopifnot(is.numeric(power))
	# stopifnot( is.numeirc(arg1) | is.numeric(arg2) | is.numeric(power) )
	return((arg1*arg2)^power)
}
przykladowa_funkcja_pomnoz_i_podnies_do_potegi(2,3,4)
przykladowa_funkcja_pomnoz_i_podnies_do_potegi(1:4,3:5,2)

#' 
#' 
#' ## Instrukcje warunkowe
#' 
## ------------------------------------------------------------------------
wypisz_badz_podnies_do_kwadratu  <- function(element){
	if(class(element)=="character"){
		cat(element)
	}else{
		if(class(element)=="numeric"){
			element^2
		}
	}
}
# niepotrzebne jest pisanie return
# ani deklarowanie typow parametrow
# ani deklarowanie typow wyjscia z funkcji
wypisz_badz_podnies_do_kwadratu("Hello")
wypisz_badz_podnies_do_kwadratu(4)

#' 
#' 
#'  

#' 
## ----child='02_Instalacja_Pakietow.Rmd'----------------------------------

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

#' 
#' 
## ----child='03_Dodatkowe_Materialy_R.Rmd'--------------------------------

#' # Dodatkowe Materiały
#' 
#' - [Przewodnik po pakiecie R, P. Biecek](http://biecek.pl/R/)
#' - [Analiza danych z programem R, P. Biecek](http://biecek.pl/Analiza.Danych/)
#' - [Programowanie w języku R, M. Gągolewski](http://rksiazka.rexamine.com/)
#' - Darmowy Kurs MOOC Analizy i Przetwarzania Danych w R - [Pogromcy Danych](http://pogromcydanych.icm.edu.pl/)
#' - [In-depth introduction to machine learning in 15 hours of expert videos](http://www.r-bloggers.com/in-depth-introduction-to-machine-learning-in-15-hours-of-expert-videos/)
#' - [Microsoft Launches Its First Free Online R Course on edX](http://www.r-bloggers.com/microsoft-launches-its-first-free-online-r-course-on-edx/)

