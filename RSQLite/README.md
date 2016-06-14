# RSQLite
Marcin Kosiński  
Ostatnia data modyfikacji `r Sys.Date()`  





# Wczytanie danych z pliku .csv

Do wczytywania danych wykorzystamy pakiet [`data.table`](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.pdf), który udostępnia funkcję
`fread` wczytującą dane z plików szybciej niż podstawowe funkcje z pakietu `base`.

```r
library(data.table)
untar('../100_dane/articles.csv.tar.gz', exdir = '../100_dane/')
articles <- fread("../100_dane/articles.csv",
									data.table = FALSE) # wczytany obiekt będzię klasy data.frame
```

Nazwy kolumn w ramce danych `articles` otrzymuje się poleceniem

```r
names(articles)
```

```
[1] "wid"      "tytul"    "zajawka"  "tresc"    "category" "entry"    "name"     "id_cc"   
```

Dodajmy kolumnę z identyfikatorem wiersza

```r
articles <- cbind(articles, id = 1:nrow(articles))
names(articles)
```

```
[1] "wid"      "tytul"    "zajawka"  "tresc"    "category" "entry"    "name"     "id_cc"    "id"      
```


# Tworzenie bazy danych

## Stworzenie bazy

Aby stworzyć tabele w planowanej bazie danych, stwórzmy mniejsze ramki danych
wybierając tylko niektóre kolumny


```r
articles_and_wids <- articles[, c("wid", "tytul", "zajawka", "tresc", "id")]
categories_and_entries <- articles[, c("category", "entry", "id")]
names_and_idCCs <- articles[, c("name", "id_cc")]
```

Dla tabel zawierających jedynie klucze, możemy wybrać tylko połączenia unikalne

```r
dim(names_and_idCCs)
```

```
[1] 95456     2
```

```r
names_and_idCCs_unique <- unique(names_and_idCCs) # zwraca unikalne wiersze
dim(names_and_idCCs_unique)
```

```
[1] 56  2
```

Połączenie z bazą danych oraz pierwotne stworzenie pliku z bazą danych


```
[1] FALSE
```

```r
library(RSQLite)
articles_conn <- dbConnect( SQLite(), dbname = "articles.db" ) 
# SQLite() - Class SQLiteDriver with constructor SQLite.
```

## Zapis do bazy

Wpisanie nowych tabel do bazy danych z gotowych ramek danych


```r
dbWriteTable( articles_conn, "articles_and_wids", articles_and_wids,
							overwrite = TRUE, row.names = FALSE )
```

```
[1] TRUE
```

```r
dbWriteTable( articles_conn, "categories_and_entries", categories_and_entries,
							overwrite = TRUE, row.names = FALSE )
```

```
[1] TRUE
```

```r
dbWriteTable( articles_conn, "names_and_idCCs", names_and_idCCs_unique,
							overwrite = TRUE, row.names = FALSE )
```

```
[1] TRUE
```

Ewentualnie: stworzenie pustych tabel



```r
dbSendQuery(articles_conn, "
            CREATE TABLE articles_and_wids_PK
            (id INTEGER PRIMARY KEY,
						wid INTEGER,
            tytul VARCHAR(50),
            zajawka VARCHAR(5000),
            tresc VARCHAR(5000000),
						UNIQUE (wid,tytul,zajawka,tresc) ON CONFLICT IGNORE)")
```

```
<SQLiteResult>
```

Wtedy `INSERT` do pustej tabeli można wykonać np. tak


```r
insert <- paste0("INSERT INTO articles_and_wids_PK (id,wid,tytul,zajawka,tresc) VALUES (",
									articles_and_wids[1, 1], ",'",
									articles_and_wids[1, 2], "','",
									articles_and_wids[1, 3], "','",
								  articles_and_wids[1, 4], "','",
									articles_and_wids[1, 5], "')")
substr(insert, start=1, stop = 200)
```

```
[1] "INSERT INTO articles_and_wids_PK (id,wid,tytul,zajawka,tresc) VALUES (8299352,'orange ekstraklasa gĂłrnik Ĺ‚Ä™czna korona kielce','obydwa zespĂłĹ‚ spotkanie duĹĽy znaczenie kolporter korona kielce wal"
```

```r
dbSendQuery(articles_conn, insert)
```

```
<SQLiteResult>
```



# Podstawowe funkcje pakietu [RSQLite](https://cran.r-project.org/web/packages/RSQLite/index.html)

Dla gotowej bazy danych można sprawdzić jakie tabele znajduą się w bazie


```r
dbListTables(articles_conn)
```

```
[1] "articles_and_wids"      "articles_and_wids_PK"   "categories_and_entries" "names_and_idCCs"       
```

lub wyświetlić kolumny w określonej tabeli

```r
dbListFields(articles_conn, "articles_and_wids")
```

```
[1] "wid"     "tytul"   "zajawka" "tresc"   "id"     
```


Można też łatwo wczytać całą tabelę z bazy


```r
dbReadTable(articles_conn, "names_and_idCCs") -> names_and_idCCs_df
head(names_and_idCCs_df)
```

```
           name   id_cc
1 PiĹ‚ka noĹĽna    1716
2 PiĹ‚ka noĹĽna   20488
3       Kultura 1010979
4 PiĹ‚ka noĹĽna 1025721
5       Kultura  113194
6 PiĹ‚ka noĹĽna   17431
```


## Wykonanie zapytania

Do wyciągania danych niewygodne jest używanie funkcji `dbSendQuery`

```r
# ?dbSendQuery - Execute a statement on a given database connection. It does not extracts any records — for that you need to use the function dbFetch, and then you must call dbClearResult when you finish fetching the records you need.
### Example
# res <- dbSendQuery(con, "SELECT * FROM mtcars WHERE cyl = 4;")
# dbFetch(res)
# dbClearResult(res)
```

Przy pomocy tej funkcji najlepiej tworzyć tabele i je usuwać i wykonywać INSERTY


```r
dbSendQuery(articles_conn, "DROP TABLE articles_and_wids_PK")
```

```
<SQLiteResult>
```

Dlatego używać będziemy funkcji `dbGetQuery`

```r
# Send query, retrieve results and then clear result set.
dbGetQuery(articles_conn, statement = "
						SELECT arts.wid, names.name
						FROM(
							SELECT wid, id 
							FROM articles_and_wids a
							LIMIT 5
						) arts
						JOIN categories_and_entries cat
						ON arts.id = cat.id
						JOIN names_and_idCCs names
						ON cat.category = names.id_cc
						" )
```

```
       wid          name
1  8299352 PiĹ‚ka noĹĽna
2  9719551 PiĹ‚ka noĹĽna
3 17517126 PiĹ‚ka noĹĽna
4 17517159 PiĹ‚ka noĹĽna
5 17505499 PiĹ‚ka noĹĽna
```

```r
# dbGetQuery comes with a default implementation that calls dbSendQuery, then if dbHasCompleted is TRUE, it uses fetch to return the results. on.exit is used to ensure the result set is always freed by dbClearResult. Subclasses should override this method only if they provide some sort of performance optimisation.
```

Używając funkcji `dbDisconnect()` można się z bazą danych rozłączyć

```r
# wazne jest by po sobie sprzątać na wypadek gdyby dane z pamięci nie zostały zapisane do pliku
dbDisconnect(articles_conn)
```

```
[1] TRUE
```




# Dodatkowe materiały

- Skrypt z najciekawszymi funkcjami RSQLita: [http://faculty.washington.edu/kenrice/sisg-adv/exampleSQLite.R](http://faculty.washington.edu/kenrice/sisg-adv/exampleSQLite.R)
- [http://www.r-bloggers.com/r-and-sqlite-part-1/](http://www.r-bloggers.com/r-and-sqlite-part-1/)
- [Using SQLite with R](http://rstudio-pubs-static.s3.amazonaws.com/8753_a57d3950027541a590c9b40a045accbf.html#1)


