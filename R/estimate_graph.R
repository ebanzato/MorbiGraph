#' Graph estimate
#'
#' @param data a matrix or a dataframe containing all the variables. All columns must be named.
#' @param v.corr A string with the names of the variables to be considered in the analysis but not in the network.
#'               By design, the associations between them will be zero.
#' @param IC information criteria to use: "AIC", "BIC", "EBIC". EBIC is the default for IsingFit.
#' @param rule rule to apply to the final edges selection, can be one of "AND" or "OR.
#'             "AND" considers the intersection of the edges, "OR" considers the union of edges.
#' @param method estimation method to use: "stepwise", "bestsubset", "IsingFit"
#' @param direction the direction of stepwise search, can be one of "backward" or "forward". Only if method='stepwise'.
#' @param conf.int confidence level required. Default is NULL and only the matrix with the p-values is returned.
#' @param OR logical. If TRUE the function returns an adjacency matrix with odds-ratios, if FALSE with log(OR).
#' @param edge.table logical. If TRUE the function returns the graph also in dataframe mode.
#'
#' @return The function returns the estimated weighted adjacency matrix of the graph and a matrix with the associated p-values.
#'         If conf.int is not null, a matrix with all the associations and the corresponding confidence intervals will be returned.
#'
#' @examples
#'
#' set.seed(1)
#' S <- solve(matrix(c(2,0.9,0,0.9,2,0.9,0,0.9,2), ncol=3))
#' df <- as.data.frame(ifelse(mvtnorm::rmvnorm(1000, sigma=S)>0,1,0))
#' colnames(df) <- paste0('X',1:3)
#'
#' g <- estimate_graph(df, method='stepwise',direction='backward')
#'
#' @export
estimate_graph = function(data, v.conf=NULL, method='IsingFit', direction=NULL, IC='EBIC', rule='AND', OR=TRUE, conf.int=NULL, edge.table=FALSE){

  # learn the structure
  if(method=='stepwise'){
    adjm = node_stepwise(data, v.conf=v.conf, direction=direction, IC=IC, rule=rule)
  }
  if(method=='bestsubset'){
    adjm = node_bestsubset(data, v.conf=v.conf, IC=IC, rule=rule)
  }
  if(method=='IsingFit'){
    AND = ifelse(rule=='AND', TRUE, FALSE)
    ifit = IsingFit::IsingFit(data, family='binomial', AND=AND, gamma=0.25, plot=FALSE)
    adjm = ifelse(ifit$weiadj == 0, 0, 1)

    if(!is.null(v.conf)){
      v.net = colnames(adjm)[!colnames(adjm) %in% v.conf]
      adjm = adjm[c(v.net,v.conf),c(v.net,v.conf)]
      adjm[v.conf,v.conf] = 0
    }
  }

  if(sum(adjm)==0){
    warning('The estimated graph has no connections')
    return(adjm)
  }else{

    # Estimate the parameters
    fTable = freq_table(data)
    mformula = mFormula(adjm)
    mod = glm(mformula, data=fTable, family='poisson')

    p = ncol(data)

    # List of weighted edges
    e.coef = coef(mod)[-c(1:(p+1))]
    e.pval = summary(mod)$coef[-c(1:(p+1)),4]
    e.df = stringr::str_split_fixed(names(e.coef), pattern=':', n=2)
    e.list = data.frame(e.df, e.coef, exp(e.coef), e.pval)
    rownames(e.list) = NULL; colnames(e.list) = c('v1', 'v2', 'coef', 'OR','pvalue')

    # Conf int
    if(!is.null(conf.int)){
      e.confint = confint(mod, parm=names(e.coef), level=conf.int)
      e.list.ci = data.frame(e.df, e.coef, exp(e.coef), e.confint); rownames(e.list.ci) = NULL
      colnames(e.list.ci) = c('v1', 'v2', 'coef', 'OR', 'lower', 'upper')
    }

    # Weighted adj matrix
    w.adj = igraph::as_adjacency_matrix(igraph::graph_from_data_frame(e.list[,-4], directed=FALSE), sparse=FALSE, attr='coef')
    p.adj = igraph::as_adjacency_matrix(igraph::graph_from_data_frame(e.list[,-3], directed=FALSE), sparse=FALSE, attr='pvalue')

    if(isTRUE(OR)){
      w.adj = exp(w.adj)
    }

    # OUTPUT
    if(isFALSE(edge.table) & is.null(conf.int)){
      return(list('w.adj'=w.adj,'p.adj'=p.adj))
    }
    if(isTRUE(edge.table)){
      if(is.null(conf.int)){
        return(list('w.adj'=w.adj,'p.adj'=p.adj, 'e.table'=e.list))
      }else{
        return(list('w.adj'=w.adj,'p.adj'=p.adj, 'e.table'=e.confint))
      }
    }
  }

}
