#' Node-wise best subset selection
#'
#' Structure learning of the graph achieved through best subset selection.
#'
#' @param data a matrix or a dataframe containing all the variables. All columns must be named.
#' @param v.conf A string with the names of the confounding variables. By design, the associations between them will be zero.
#' @param IC information criteria to use: "AIC", "BIC", "EBIC".
#' @param rule Can be "AND" or "OR" to indicate whether the AND-rule or the OR-rule should be used to define the edges in the network.
#'
#' @return The function returns the estimated adjacency matrix of the graph.
#'
#' @examples
#' set.seed(1)
#' df <- matrix(sample(0:1,100,replace=T),ncol=5)
#' colnames(df) <- paste0('X',1:5)
#' node_bestsubset(df)
#'
#' @export
#'
node_bestsubset = function(data, v.conf=NULL, IC='BICg', rule='AND'){

  if(!IC %in% c("AIC", "BIC", "EBIC")){
    stop('\'IC\' should be one of "AIC", "EBIC"')
  }

  if(is.matrix(data)){
    if(is.null(colnames(data))){
      stop('variables names are missing.')
    }
    data = as.data.frame(data)
  }

  p = ncol(data)

  if(!is.null(v.conf)){
    dataM = as.data.frame(data[,!colnames(data) %in% v.conf])
    dataC = as.data.frame(data[, colnames(data) %in% v.conf])
    colnames(dataC) = colnames(data)[colnames(data) %in% v.conf]
    v.net = colnames(dataM)
    data = cbind(dataM,dataC)
  }else{
    v.net = colnames(data)
  }

  adj1 = matrix(FALSE, ncol=p, nrow=p)
  colnames(adj1) = rownames(adj1) = colnames(data)

  ## Xi|X_i regressions
  # if v.corr is not null, do not consider those variables
  for(varY in 1:length(v.net)){
    print(paste0(varY,'/',length(v.net)))

    # order the variables in Xy
    varXY = c(colnames(data)[-varY], colnames(data)[varY])

    m.selec = bestglm::bestglm(data[,varXY], family=binomial, IC=IC)

    v.select = as.matrix(m.selec$BestModels[1,-ncol(m.selec$BestModels)],nrow=1)
    adj1[varY, -varY] = v.select
  }

  # Select edges with the and/or rule
  mat_and_or = adj1 + t(adj1)

  # The and/or rule is applied only to the X matrix, without considering v.corr

  if(rule=='AND'){
    adjmatrix = ifelse(mat_and_or[1:length(v.net),1:length(v.net)] == 2, 1, 0)
  }
  if(rule=='OR'){
    adjmatrix = ifelse(mat_and_or[1:length(v.net),1:length(v.net)] > 0, 1, 0)
  }

  # out
  mat_and_or[1:length(v.net),1:length(v.net)] = adjmatrix
  return(mat_and_or)
}
