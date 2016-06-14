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
