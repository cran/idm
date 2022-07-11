tab.disjonctif <- function (tab) 
{
  tab <- data.frame(tab, stringsAsFactors = TRUE)
  modalite.disjonctif <- function(i) {
    moda <- as.factor(tab[, i])
    n <- length(moda)
    x <- matrix(0L, n, nlevels(moda))
    x[(1:n) + n * (unclass(moda) - 1L)] <- 1L
    return(x)
  }
  if (ncol(tab) == 1) {
    res <- modalite.disjonctif(1)
    dimnames(res) <- list(attributes(tab)$row.names, levels(tab[, 
                                                                1]))
  }
  else {
    variable <- rep(attributes(tab)$names, lapply(tab, nlevels))
    listModa <- unlist(lapply(tab, levels))
    wlistModa <- which((listModa) %in% c("y", "n", "Y", "N"))
    if (!is.null(wlistModa)) 
      listModa[wlistModa] <- paste(variable[wlistModa], 
                                   listModa[wlistModa], sep = ".")
    numlistModa <- which(unlist(lapply(listModa, is.numeric)))
    if (!is.null(numlistModa)) 
      listModa[numlistModa] <- paste(variable[numlistModa], 
                                     listModa[numlistModa], sep = ".")
    res <- lapply(1:ncol(tab), modalite.disjonctif)
    res <- as.matrix(data.frame(res, check.names = FALSE, stringsAsFactors = TRUE))
    dimnames(res) <- list(attributes(tab)$row.names, listModa)
  }
  return(res)
}