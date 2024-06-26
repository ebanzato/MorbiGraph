% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_density.R
\name{graph_density}
\alias{graph_density}
\title{Density index for graphs}
\usage{
graph_density(graph)
}
\arguments{
\item{graph}{An igraph object or an adjacency matrix.
It can also be a weighted matrix, which will be treated as an adjacency matrix.}
}
\value{
a number indicating the density of the graph.
Values of density can vary between 0 and 1, where 0 means that the graph has no edges,
while a coefficient equal to 1 represents a complete graph.
}
\description{
Compute the density index for a graph.
}
\details{
The density index is computed as the ratio between the number of estimated edges
and the total number of possible edges.
The total number of possible edges for a graph with p nodes is p(p-1)/2.
}
\examples{

g <- igraph::make_graph(c(1,2,2,3,3,1,3,4,4,5,5,6,6,3,6,4),directed=F)

d <- graph_density(g)

}
