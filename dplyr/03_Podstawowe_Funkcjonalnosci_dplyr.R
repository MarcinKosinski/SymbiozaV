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


