checkNullTable <- function(toCheck) {
  if (is.null(toCheck))
    return(data.frame("Nothing" = ""))
  else
    return(toCheck)
}