#' Yule's Q
#'
#' Compute the Yule coefficient of association.
#'
#' The coefficient is computed as
#' Q = (OR-1)/(OR+1),
#' where OR is the odds ratio.
#'
#' @param OR A vector or a symmetric matrix of odds ratios.
#'
#' @return A vector or a matrix with Q coefficients computed from the OR.
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
qYule = function(OR){
  if(is.matrix(OR)){
    if(isSymmetric(OR)){
      Q = (OR-1)/(OR+1)
      return(Q)
    }else{
      print('The OR matrix is not symmetric.')
    }
  }else{
  Q = (OR-1)/(OR+1)
  return(Q)
  }
}
