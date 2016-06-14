#' ---
#' title: "03 dplyr"
#' subtitle: "Grammar of Data Manipulation"
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
## ----child='01_Zajawka.Rmd'----------------------------------------------

#' # Dlaczego warto znać dplyr'a?
#' 
#' [Data Wrangling with dplyr and tidyr](https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)
#' 
#' Praktycznie każda baza danych różni się listą zaimplementowanych funkcjonalności czy agregatów.
#' 
#' Jeżeli pracujemy z jedną bazą danych to może nam to nie doskiwerać, ale jeżeli trzeba jednocześnie korzystać z MSSQL'a, MySQLa i Postgresa?
#' 
#' Wielu problemów można sobie oszczędzić używając pośrednika do komunikacji z bazą danych takiego jak [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html). Pozwala on do pewnego stopnia na pracę z danymi bez zastanawiania się gdzie te dane aktualnie są i jaką funkcją w aktualnej bazie danych liczy się średnią.
#' 
## ------------------------------------------------------------------------
library(RSQLite)
library(dplyr)

#' 
#' Połączenie z bazą danych
## ------------------------------------------------------------------------
getwd() # zwraca obecną lokalizację ścieżki roboczej
articles_conn_dplyr <- src_sqlite(path = "../02_RSQLite/articles.db")

#' 
#' Obiekt klasy `tbl_sql` można utworzyć podając nazwę tabeli lub pisząc zapytanie typu SELECT 
## ------------------------------------------------------------------------
articles_and_wids_tbl <- tbl(articles_conn_dplyr, "articles_and_wids")
articles_all <- tbl(articles_conn_dplyr, sql(
														 "SELECT *
														 FROM articles_and_wids arts
														 JOIN categories_and_entries cat
														 ON arts.id = cat.id
														 JOIN names_and_idCCs names
														 ON cat.category = names.id_cc "))

#' 
#' Na takim obiekcie można pracować wykorzystując funkcje z pakietu `dplyr`
#' 
## ------------------------------------------------------------------------
# ile jest różnych artykułów przypisanych do danej kategorii
articles_all %>%
	group_by(name) %>%
	summarise(count=n()) %>%
	arrange(desc(count)) -> aggr1
aggr1

#' 
## ------------------------------------------------------------------------
# ile kategorii ma jeden artykuł
articles_all %>%
	group_by(wid) %>%
	summarise(do_ilu_przypisanych_kategorii=n()) %>%
	group_by(do_ilu_przypisanych_kategorii) %>%
	summarise(ile_artykulow=n()) %>%
	arrange(do_ilu_przypisanych_kategorii)-> aggr2
aggr2

#' 
#' Bardzo pożyteczną funkcjonalnością tego pakietu jest odwołanie się do zapytania,
#' którego normalnie trzeba byłoby użyć aby uzyskać taki sam wynik jak przy użyciu funkcji
#' z pakietu R.
#' 
## ------------------------------------------------------------------------
aggr1$query
aggr2$query

#' 
#' *Wyjaśnienie funkcji i składni znajduje się w następnym rozdziale.*
#' 
#' 

#' 
## ----child='02_Operator_Pipe.Rmd'----------------------------------------

#' # Operator Pipe/Then/magrittr [`%>%`](https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html)
#' 
#' 
#' Ten operator pochodzi z pakietu `magrittr` (cytując z jego dokumentacji: **to be pronounced with a sophisticated french accent**) i jest dostępny po włączeniu pakietu `dplyr`.
#' 
#' Jak działa ten operator?
#' 
#' Przekazuje lewą stronę operatora jako pierwszy argument prawej strony tego operatora.
#' 
#' Instrukcja `a %>% f(b)` jest równoważna instrukcji `f(a, b)`.
#' 
#' *Przykład.*
#' 
#' Instrukcja
#' 
## ------------------------------------------------------------------------
summary(dplyr::select(top_n(iris,4), contains("Petal") ) )
# ktora moze byc zapisana 
summary(
	dplyr::select(
		top_n(iris,
					4),
		contains("Petal")
		)
	)
# badz ewentualnie
iris_top4 <- top_n(iris,4)
iris_top4_selected <- dplyr::select(iris_top4, contains("Petal"))
summary(iris_top4_selected)

#' 
#' jest równoważna
#' 
## ------------------------------------------------------------------------
iris %>%
	top_n(4) %>% #wybranie pierwszych 4 obserwacji, dziala i przy grupowaniu
	dplyr::select(contains("Petal")) %>% #wybiera kolumny zawierające w nazwie "Petal"
	summary #zwraca statystyki agregacyjne
	


#' 
#' która jest zdecydowanie czytelniejsza i nie wymaga niepotrzebnej deklaracji zmiennych po drodze.

#' 
## ----child='03_Podstawowe_Funkcjonalnosci_dplyr.Rmd'---------------------

#' # Podstawowe funkcje pakietu dplyr
#' 
#' ## Filtorwanie wierszy
#' 
## ------------------------------------------------------------------------
filter(articles_all,
				 name == "Koszykówka",
				 wid > 481) -> fun1
fun1$query

#' 
#' 
#' ## Dodawanie nowej kolumny
#' 
## ------------------------------------------------------------------------
mutate(articles_all,
			 wid_times_entry = wid*entry) -> fun2
# liczba kolumn
ncol(articles_all)
ncol(fun2) 

#' 
#' ## Wybór zmiennych
#' 
## ------------------------------------------------------------------------
select(articles_all,
			 wid, zajawka) -> fun3
ncol(fun3)

select(articles_all,
			 contains("a")) -> fun4
ncol(fun4)

select(articles_all,
			 -wid, -tresc) -> fun5
ncol(fun5)


#' 
#' 
#' ## Grupowanie i Podsumowania
#' 
## ------------------------------------------------------------------------
# pakiet do pracy z napisami
library(stringi)

# statystyki agregacyjne slow w artykule w podziale na kategorie
articles_all %>%
	as.data.frame  %>% # wczytuje dane do R, tak by mozna bylo
	# uzyc funkcji `stri_count_words` z wewnetrznej biblioteki
	mutate(liczba_slow = stri_count_words(tresc))  %>%
	group_by(name) %>%
	summarise(min_slow = min(liczba_slow),
						sr_slow = mean(liczba_slow) %>% round(2),
						max_slow = max(liczba_slow),
						ile_artykulow = n()) -> aggr3
aggr3

#' 
#' ## Sortowanie wierszy
#' 
## ------------------------------------------------------------------------
aggr3 %>%
	arrange(desc(ile_artykulow))

aggr3 %>%
	arrange(ile_artykulow)

arrange(aggr3, ile_artykulow)

#' 
#' ## Łączenie tabel
#' 
#' 
## ------------------------------------------------------------------------
getwd()
articles_conn <- dbConnect( SQLite(), dbname = "../02_RSQLite/articles.db" )
dbListTables(articles_conn)
articles_and_wids_tbl <- tbl(articles_conn_dplyr, "articles_and_wids")
categories_and_entries_tbl <- tbl(articles_conn_dplyr,"categories_and_entries")
categories_and_entries_tbl %>%
	rename(wid = entry) %>%
	left_join(y = articles_and_wids_tbl,
						by = "wid") -> joined_tbls

joined_tbls %>%
	select(category,wid,tytul) %>%
	head



#' 
## ----child='04_Dodatkowe_Materialy_dplyr.Rmd'----------------------------

#' # Dodatkowe materiały
#' 
#' - Darmowy kurs MOOC [Pogromcy Danych](http://pogromcydanych.icm.edu.pl/)
#' - [Data Wrangling with dplyr and tidyr](https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)
#' - [Data frames](https://cran.r-project.org/web/packages/dplyr/vignettes/data_frames.html)
#' - [Databases](https://cran.r-project.org/web/packages/dplyr/vignettes/databases.html)
#' - [Hybrid evaluation](https://cran.r-project.org/web/packages/dplyr/vignettes/hybrid-evaluation.html)
#' - [Introduction to dplyr](https://cran.r-project.org/web/packages/dplyr/vignettes/introduction.html)
#' - [Adding a new SQL backend](https://cran.r-project.org/web/packages/dplyr/vignettes/new-sql-backend.html)
#' - [Non-standard evaluation](https://cran.r-project.org/web/packages/dplyr/vignettes/nse.html)
#' - [Two-table verbs](https://cran.r-project.org/web/packages/dplyr/vignettes/two-table.html)
#' - [Window functions and grouped mutate/filter](https://cran.r-project.org/web/packages/dplyr/vignettes/window-functions.html)

