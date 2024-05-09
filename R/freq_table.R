#' Frequency table of large dataframe
#'
#' Compute a summarized matrix of the contingency table of the dataframe.
#'
#' @param data A matrix or a dataframe where each column is a variable.
#' @param ordered A logical indicating whether to order the output matrix
#' based on the decreasing ordering of the cell frequencies.
#'
#' @return A matrix with all the possible combinations of the variables
#' and a column with the frequency of each combination.
#'
#' @examples
#'
#' set.seed(1)
#' x <- sample(c(0,1), 100, replace = T)
#' y <- sample(c(0,1), 100, replace = T)
#' z <- sample(c(0,1), 100, replace = T)
#' df <- cbind(x, y, z)
#'
#' ftab <- freq_table(df)
#'
#' @keywords internal
#'
freq_table = function(data, ordered=TRUE){

  if(is.matrix(data)){
    data = as.data.frame(data)
  }

  # variables
  name.var = colnames(data)

  # aggregate
  fTable = data %>%
    dplyr::group_by_all() %>%
    dplyr::tally() %>%
    dplyr::ungroup()

  # create all the combinations
  combo = rep(list(c(0,1)),length(name.var))
  names(combo) = name.var
  combo = expand.grid(combo)

  # frequency table
  tabfreq = merge(combo, fTable, all=T)
  tabfreq$n[is.na(tabfreq$n)] = 0

  # order the matrix in decreasing order
  if(ordered==TRUE){
    tabfreq = tabfreq[order(tabfreq$n, decreasing = T),]
  }

  return(tabfreq)
}
