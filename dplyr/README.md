# 03 dplyr
Marcin Kosiński  
Ostatnia data modyfikacji `r Sys.Date()`  






# Dlaczego warto znać dplyr'a?

[Data Wrangling with dplyr and tidyr](https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)

Praktycznie każda baza danych różni się listą zaimplementowanych funkcjonalności czy agregatów.

Jeżeli pracujemy z jedną bazą danych to może nam to nie doskiwerać, ale jeżeli trzeba jednocześnie korzystać z MSSQL'a, MySQLa i Postgresa?

Wielu problemów można sobie oszczędzić używając pośrednika do komunikacji z bazą danych takiego jak [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html). Pozwala on do pewnego stopnia na pracę z danymi bez zastanawiania się gdzie te dane aktualnie są i jaką funkcją w aktualnej bazie danych liczy się średnią.


```r
library(RSQLite)
library(dplyr)
```

Połączenie z bazą danych

```r
getwd() # zwraca obecną lokalizację ścieżki roboczej
```

```
[1] "/home/mkosinski/codepot/codepot-workshop-2015/03_dplyr"
```

```r
articles_conn_dplyr <- src_sqlite(path = "../02_RSQLite/articles.db")
```

Obiekt klasy `tbl_sql` można utworzyć podając nazwę tabeli lub pisząc zapytanie typu SELECT 

```r
articles_and_wids_tbl <- tbl(articles_conn_dplyr, "articles_and_wids")
articles_all <- tbl(articles_conn_dplyr, sql(
														 "SELECT *
														 FROM articles_and_wids arts
														 JOIN categories_and_entries cat
														 ON arts.id = cat.id
														 JOIN names_and_idCCs names
														 ON cat.category = names.id_cc "))
```

Na takim obiekcie można pracować wykorzystując funkcje z pakietu `dplyr`


```r
# ile jest różnych artykułów przypisanych do danej kategorii
articles_all %>%
	group_by(name) %>%
	summarise(count=n()) %>%
	arrange(desc(count)) -> aggr1
aggr1
```

```
Source: sqlite 3.8.6 [../02_RSQLite/articles.db]
From: <derived table> [?? x 2]
Arrange: desc(count) 

                name count
1        Piłka nożna 45505
2  Biznes i ekonomia 21270
3            Kobieta  5688
4          Siatkówka  5687
5         Koszykówka  5157
6          Lifestyle  4464
7             Pogoda  3913
8            Kultura  3772
..               ...   ...
```


```r
# ile kategorii ma jeden artykuł
articles_all %>%
	group_by(wid) %>%
	summarise(do_ilu_przypisanych_kategorii=n()) %>%
	group_by(do_ilu_przypisanych_kategorii) %>%
	summarise(ile_artykulow=n()) %>%
	arrange(do_ilu_przypisanych_kategorii)-> aggr2
aggr2
```

```
Source: sqlite 3.8.6 [../02_RSQLite/articles.db]
From: <derived table> [?? x 2]
Arrange: do_ilu_przypisanych_kategorii 

   do_ilu_przypisanych_kategorii ile_artykulow
1                              1         64464
2                              2         15451
3                              3            30
..                           ...           ...
```

Bardzo pożyteczną funkcjonalnością tego pakietu jest odwołanie się do zapytania,
którego normalnie trzeba byłoby użyć aby uzyskać taki sam wynik jak przy użyciu funkcji
z pakietu R.


```r
aggr1$query
```

```
<Query> SELECT "name", "count"
FROM (SELECT "name", COUNT() AS "count"
FROM (SELECT *
														 FROM articles_and_wids arts
														 JOIN categories_and_entries cat
														 ON arts.id = cat.id
														 JOIN names_and_idCCs names
														 ON cat.category = names.id_cc ) AS "_W1"
GROUP BY "name") AS "_W2"
ORDER BY "count" DESC
<SQLiteConnection>
```

```r
aggr2$query
```

```
<Query> SELECT "do_ilu_przypisanych_kategorii", "ile_artykulow"
FROM (SELECT "do_ilu_przypisanych_kategorii", COUNT() AS "ile_artykulow"
FROM (SELECT "wid", COUNT() AS "do_ilu_przypisanych_kategorii"
FROM (SELECT *
														 FROM articles_and_wids arts
														 JOIN categories_and_entries cat
														 ON arts.id = cat.id
														 JOIN names_and_idCCs names
														 ON cat.category = names.id_cc ) AS "_W1"
GROUP BY "wid") AS "_W3"
GROUP BY "do_ilu_przypisanych_kategorii") AS "_W4"
ORDER BY "do_ilu_przypisanych_kategorii"
<SQLiteConnection>
```

*Wyjaśnienie funkcji i składni znajduje się w następnym rozdziale.*




# Operator Pipe/Then/magrittr [`%>%`](https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html)


Ten operator pochodzi z pakietu `magrittr` (cytując z jego dokumentacji: **to be pronounced with a sophisticated french accent**) i jest dostępny po włączeniu pakietu `dplyr`.

Jak działa ten operator?

Przekazuje lewą stronę operatora jako pierwszy argument prawej strony tego operatora.

Instrukcja `a %>% f(b)` jest równoważna instrukcji `f(a, b)`.

*Przykład.*

Instrukcja


```r
summary(dplyr::select(top_n(iris,4), contains("Petal") ) )
```

```
  Petal.Length    Petal.Width   
 Min.   :4.500   Min.   :1.400  
 1st Qu.:5.100   1st Qu.:1.800  
 Median :5.550   Median :2.000  
 Mean   :5.552   Mean   :2.026  
 3rd Qu.:5.875   3rd Qu.:2.300  
 Max.   :6.900   Max.   :2.500  
```

```r
# ktora moze byc zapisana 
summary(
	dplyr::select(
		top_n(iris,
					4),
		contains("Petal")
		)
	)
```

```
  Petal.Length    Petal.Width   
 Min.   :4.500   Min.   :1.400  
 1st Qu.:5.100   1st Qu.:1.800  
 Median :5.550   Median :2.000  
 Mean   :5.552   Mean   :2.026  
 3rd Qu.:5.875   3rd Qu.:2.300  
 Max.   :6.900   Max.   :2.500  
```

```r
# badz ewentualnie
iris_top4 <- top_n(iris,4)
iris_top4_selected <- dplyr::select(iris_top4, contains("Petal"))
summary(iris_top4_selected)
```

```
  Petal.Length    Petal.Width   
 Min.   :4.500   Min.   :1.400  
 1st Qu.:5.100   1st Qu.:1.800  
 Median :5.550   Median :2.000  
 Mean   :5.552   Mean   :2.026  
 3rd Qu.:5.875   3rd Qu.:2.300  
 Max.   :6.900   Max.   :2.500  
```

jest równoważna


```r
iris %>%
	top_n(4) %>% #wybranie pierwszych 4 obserwacji, dziala i przy grupowaniu
	dplyr::select(contains("Petal")) %>% #wybiera kolumny zawierające w nazwie "Petal"
	summary #zwraca statystyki agregacyjne
```

```
  Petal.Length    Petal.Width   
 Min.   :4.500   Min.   :1.400  
 1st Qu.:5.100   1st Qu.:1.800  
 Median :5.550   Median :2.000  
 Mean   :5.552   Mean   :2.026  
 3rd Qu.:5.875   3rd Qu.:2.300  
 Max.   :6.900   Max.   :2.500  
```

która jest zdecydowanie czytelniejsza i nie wymaga niepotrzebnej deklaracji zmiennych po drodze.


# Podstawowe funkcje pakietu dplyr

## Filtorwanie wierszy


```r
filter(articles_all,
				 name == "Koszykówka",
				 wid > 481) -> fun1
fun1$query
```

```
<Query> SELECT "wid", "tytul", "zajawka", "tresc", "id", "category", "entry", "id:1", "name", "id_cc"
FROM (SELECT *
														 FROM articles_and_wids arts
														 JOIN categories_and_entries cat
														 ON arts.id = cat.id
														 JOIN names_and_idCCs names
														 ON cat.category = names.id_cc ) AS "_W1"
WHERE "name" = 'Koszykówka' AND "wid" > 481.0
<SQLiteConnection>
```


## Dodawanie nowej kolumny


```r
mutate(articles_all,
			 wid_times_entry = wid*entry) -> fun2
# liczba kolumn
ncol(articles_all)
```

```
[1] 10
```

```r
ncol(fun2) 
```

```
[1] 11
```

## Wybór zmiennych


```r
select(articles_all,
			 wid, zajawka) -> fun3
ncol(fun3)
```

```
[1] 2
```

```r
select(articles_all,
			 contains("a")) -> fun4
ncol(fun4)
```

```
[1] 3
```

```r
select(articles_all,
			 -wid, -tresc) -> fun5
ncol(fun5)
```

```
[1] 8
```


## Grupowanie i Podsumowania


```r
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
```

```
Source: local data frame [8 x 5]

               name min_slow sr_slow max_slow ile_artykulow
1 Biznes i ekonomia        1  277.61     3103         21270
2           Kobieta        0  202.31    81195          5688
3        Koszykówka        0  171.17     1572          5157
4           Kultura        0  262.64     5393          3772
5         Lifestyle        0  117.96     2322          4464
6       Piłka nożna        0  145.15     2275         45505
7            Pogoda        0  177.81     2313          3913
8         Siatkówka        0  172.45     1113          5687
```

## Sortowanie wierszy


```r
aggr3 %>%
	arrange(desc(ile_artykulow))
```

```
Source: local data frame [8 x 5]

               name min_slow sr_slow max_slow ile_artykulow
1       Piłka nożna        0  145.15     2275         45505
2 Biznes i ekonomia        1  277.61     3103         21270
3           Kobieta        0  202.31    81195          5688
4         Siatkówka        0  172.45     1113          5687
5        Koszykówka        0  171.17     1572          5157
6         Lifestyle        0  117.96     2322          4464
7            Pogoda        0  177.81     2313          3913
8           Kultura        0  262.64     5393          3772
```

```r
aggr3 %>%
	arrange(ile_artykulow)
```

```
Source: local data frame [8 x 5]

               name min_slow sr_slow max_slow ile_artykulow
1           Kultura        0  262.64     5393          3772
2            Pogoda        0  177.81     2313          3913
3         Lifestyle        0  117.96     2322          4464
4        Koszykówka        0  171.17     1572          5157
5         Siatkówka        0  172.45     1113          5687
6           Kobieta        0  202.31    81195          5688
7 Biznes i ekonomia        1  277.61     3103         21270
8       Piłka nożna        0  145.15     2275         45505
```

```r
arrange(aggr3, ile_artykulow)
```

```
Source: local data frame [8 x 5]

               name min_slow sr_slow max_slow ile_artykulow
1           Kultura        0  262.64     5393          3772
2            Pogoda        0  177.81     2313          3913
3         Lifestyle        0  117.96     2322          4464
4        Koszykówka        0  171.17     1572          5157
5         Siatkówka        0  172.45     1113          5687
6           Kobieta        0  202.31    81195          5688
7 Biznes i ekonomia        1  277.61     3103         21270
8       Piłka nożna        0  145.15     2275         45505
```

## Łączenie tabel



```r
getwd()
```

```
[1] "/home/mkosinski/codepot/codepot-workshop-2015/03_dplyr"
```

```r
articles_conn <- dbConnect( SQLite(), dbname = "../02_RSQLite/articles.db" )
dbListTables(articles_conn)
```

```
[1] "articles_and_wids"      "categories_and_entries" "names_and_idCCs"       
```

```r
articles_and_wids_tbl <- tbl(articles_conn_dplyr, "articles_and_wids")
categories_and_entries_tbl <- tbl(articles_conn_dplyr,"categories_and_entries")
categories_and_entries_tbl %>%
	rename(wid = entry) %>%
	left_join(y = articles_and_wids_tbl,
						by = "wid") -> joined_tbls

joined_tbls %>%
	select(category,wid,tytul) %>%
	head
```

```
  category      wid                                                                tytul
1     1716  8299352                       orange ekstraklasa górnik łęczna korona kielce
2     1716  9719551                                                     kara van bommela
3     1716 17517126        barcelona bayern marc andre ter stegen trochę niesprawiedliwy
4     1716 17517159        barcelona bayern hiszpański media messi zagrać kwadrans finał
5     1716 17505499 ilkay gundogan bayern borussia wolałalby uniknąć zablokować transfer
6     1716 17505499 ilkay gundogan bayern borussia wolałalby uniknąć zablokować transfer
```


# Dodatkowe materiały

- Darmowy kurs MOOC [Pogromcy Danych](http://pogromcydanych.icm.edu.pl/)
- [Data Wrangling with dplyr and tidyr](https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)
- [Data frames](https://cran.r-project.org/web/packages/dplyr/vignettes/data_frames.html)
- [Databases](https://cran.r-project.org/web/packages/dplyr/vignettes/databases.html)
- [Hybrid evaluation](https://cran.r-project.org/web/packages/dplyr/vignettes/hybrid-evaluation.html)
- [Introduction to dplyr](https://cran.r-project.org/web/packages/dplyr/vignettes/introduction.html)
- [Adding a new SQL backend](https://cran.r-project.org/web/packages/dplyr/vignettes/new-sql-backend.html)
- [Non-standard evaluation](https://cran.r-project.org/web/packages/dplyr/vignettes/nse.html)
- [Two-table verbs](https://cran.r-project.org/web/packages/dplyr/vignettes/two-table.html)
- [Window functions and grouped mutate/filter](https://cran.r-project.org/web/packages/dplyr/vignettes/window-functions.html)
