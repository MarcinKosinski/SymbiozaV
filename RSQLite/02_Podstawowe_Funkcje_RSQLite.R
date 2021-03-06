#' # Podstawowe funkcje pakietu [RSQLite](https://cran.r-project.org/web/packages/RSQLite/index.html)
#' 
#' Dla gotowej bazy danych można sprawdzić jakie tabele znajduą się w bazie
#' 
## ------------------------------------------------------------------------
dbListTables(articles_conn)

#' 
#' lub wyświetlić kolumny w określonej tabeli
## ------------------------------------------------------------------------
dbListFields(articles_conn, "articles_and_wids")

#' 
#' 
#' Można też łatwo wczytać całą tabelę z bazy
#' 
## ------------------------------------------------------------------------
dbReadTable(articles_conn, "names_and_idCCs") -> names_and_idCCs_df
head(names_and_idCCs_df)

#' 
#' 
#' ## Wykonanie zapytania
#' 
#' Do wyciągania danych niewygodne jest używanie funkcji `dbSendQuery`
## ------------------------------------------------------------------------
# ?dbSendQuery - Execute a statement on a given database connection. It does not extracts any records — for that you need to use the function dbFetch, and then you must call dbClearResult when you finish fetching the records you need.
### Example
# res <- dbSendQuery(con, "SELECT * FROM mtcars WHERE cyl = 4;")
# dbFetch(res)
# dbClearResult(res)

#' 
#' Przy pomocy tej funkcji najlepiej tworzyć tabele i je usuwać i wykonywać INSERTY
#' 
## ------------------------------------------------------------------------
dbSendQuery(articles_conn, "DROP TABLE articles_and_wids_PK")

#' 
#' Dlatego używać będziemy funkcji `dbGetQuery`
## ------------------------------------------------------------------------
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
# dbGetQuery comes with a default implementation that calls dbSendQuery, then if dbHasCompleted is TRUE, it uses fetch to return the results. on.exit is used to ensure the result set is always freed by dbClearResult. Subclasses should override this method only if they provide some sort of performance optimisation.


#' 
#' Używając funkcji `dbDisconnect()` można się z bazą danych rozłączyć
## ------------------------------------------------------------------------
# wazne jest by po sobie sprzątać na wypadek gdyby dane z pamięci nie zostały zapisane do pliku
dbDisconnect(articles_conn)

#' 
#' 
