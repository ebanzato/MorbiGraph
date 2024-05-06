#' Formula
#'
#' Define the model formula from an adjacency matrix.
#'
#' @param adjm an adjacency matrix.
#'
#' @return A formula for the model.
#'
#' @keywords internal
mFormula = function(adjm){

  adjmatrix = adjm
  adjmatrix[upper.tri(adjmatrix)] = 0

  paste0('n ~ ',paste0(c(colnames(adjmatrix),
                         sapply(which(colSums(adjmatrix)>0), function(j) {
                           n.col = colnames(adjmatrix)[j]
                           n.row = names(adjmatrix[,j])[adjmatrix[,j]==1]
                           paste0(sapply(n.row, function(x) paste0(n.col,':',x)),collapse='+')
                         })),collapse = '+'))
}
