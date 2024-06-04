#' Odds ratio for two dichotomous variables
#'
#' @param x A vector.
#' @param y A vector.
#'
#' @return A number. First it creates a 2x2 table and then it computes the odds ratio.
#'
#' @examples
#'
#' set.seed(1)
#' x <- sample(c(0,1), 100, replace = T)
#' y <- sample(c(0,1), 100, replace = T)
#'
#' or <- OR_tab(x,y)
#'
#' @keywords internal
OR_tab = function(x, y){
  tab = prop.table(table(x, y))
  num = (tab[1,1]*tab[2,2])
  den = (tab[1,2]*tab[2,1])
  if(abs(den)<0.00001){
    warning('OR undefined')
  }
  OR  = num/den
  return(OR)
}
