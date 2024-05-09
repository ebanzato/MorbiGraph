#' Yule's Q
#'
#' Compute the Yule coefficient of association.
#'
#' The coefficient is computed as
#' Q = (OR-1)/(OR+1),
#' where OR is the odds ratio.
#'
#' @param x A symmetric matrix of odds ratios or a dataframe with variables
#'             for which to compute the coefficients.
#' @param OR logical indicating if the input is a matrix of odds ratios (TRUE) or a dataframe (FALSE).
#'
#' @return A symmetric matrix with the Q coefficients computed from the OR or the dataframe.
#'
#' Values of Q can vary between -1 and 1.
#' −1 reflects total negative association,
#' +1 reflects perfect positive association,
#' and 0 reflects no association.
#'
#' @examples
#'
#' # single value
#' qYule(1.5)
#'
#' # matrix of odds ratios
#' OR <- matrix(c(1, 3.5, 0.1,
#'                3.5, 1, 1.3,
#'                0.1, 1.3, 1), ncol=3)
#' qYule(OR)
#'
#' @export
#'
qYule = function(x, OR=TRUE){

  if(isTRUE(OR)){
    if(isSymmetric(x)){
      Q = (x-1)/(x+1)
      return(Q)
    }else{
      print('The OR matrix is not symmetric.')
    }
  }else{
    or = OR(x)
    Q = (x-1)/(x+1)
    return(Q)
  }
}
