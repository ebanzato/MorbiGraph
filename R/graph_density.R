#' Density index for graphs
#'
#' Compute the density index for a graph.
#'
#' The density index is computed as the ratio between the number of estimated edges
#' and the total number of possible edges.
#' The total number of possible edges for a graph with p nodes is p(p-1)/2.
#'
#' @param graph An igraph object or an adjacency matrix.
#' It can also be a weighted matrix, which will be treated as an adjacency matrix.
#'
#' @return a number indicating the density of the graph.
#' Values of density can vary between 0 and 1, where 0 means that the graph has no edges,
#' while a coefficient equal to 1 represents a complete graph.
#'
#' @examples
#'
#' library(igraph)
#' g <- igraph::make_graph(c(1,2,2,3,3,1,3,4,4,5,5,6,6,3,6,4),directed=F)
#'
#' d <- graph_density(g)
#'
#' @export
graph_density = function(graph) {

  if(is.matrix(graph) | is_igraph(graph)){
    if(is.matrix(graph)){
      p = ncol(graph)
      e_estim = sum(ifelse(graph>0,1,0)[upper.tri(ifelse(graph>0,1,0))])
      e_total = (p*(p-1)/2)
      out = e_estim/e_total
      return(out)
    }

    if(is_igraph(graph)){
      graph = as.matrix(as_adjacency_matrix(graph))
      p = ncol(graph)
      e_estim = sum(ifelse(graph>0,1,0)[upper.tri(ifelse(graph>0,1,0))])
      e_total = (p*(p-1)/2)
      out = e_estim/e_total
      return(out)
    }

  }else{
    stop('graph is not an igraph object or an adjacency matrix.')
  }

}
