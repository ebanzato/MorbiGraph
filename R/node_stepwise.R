#' Node-wise stepwise selection
#'
#' Structure learning of the graph achieved through stepwise selection.
#'
#' @param data a matrix or a dataframe containing all the variables. All columns must be named.
#' @param v.conf A string with the names of the confounding variables. By design, the associations between them will be zero.
#' @param direction the mode of stepwise search, can be one of "backward" or "forward".
#' @param IC information criteria to use: "AIC", "BIC", "EBIC".
#' @param rule Can be "AND" or "OR" to indicate whether the AND-rule or the OR-rule should be used to define the edges in the network.
#'
#' @return The function returns the estimated adjacency matrix of the graph.
#'
#' @examples
#'
#' set.seed(1)
#' df <- matrix(sample(0:1,100,replace=T),ncol=5)
#' colnames(df) <- paste0('X',1:5)
#'
#' net <- node_stepwise(df)
#'
#' @export
#'
node_stepwise = function(data, v.conf=NULL, direction='backward', IC='EBIC', rule='AND'){

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

  # select the penalty
  if(IC=='AIC'){
    k=2               # 2 k - 2 logLik
  }
  if(IC=='BIC'){
    k=log(nrow(data)) # k log(n) - 2 logLik
  }
  if(IC=='EBIC'){
    k=log(nrow(data)) + 2*0.25*log(p-1)
    # k (log(n) + 2 gamma log(p)) - 2 logLik, p: number covariates
  }
  if(!IC %in% c('AIC','BIC','EBIC')){
    stop('\'IC\' should be one of "AIC", "BIC", "EBIC"')
  }

  adj1 = matrix(FALSE, ncol=p, nrow=p)
  colnames(adj1) = rownames(adj1) = colnames(data)

  ## Xi|X_i regressions
  # if v.corr is not null, do not consider those variables
  for(varY in 1:length(v.net)){
    print(paste0(varY,'/',length(v.net)))

    full = paste0(v.net[varY],' ~ .')
    empt = paste0(v.net[varY],' ~ 1')

    m.full = glm(full, data=data, family=binomial)
    m.empt = glm(empt, data=data, family=binomial)

    scope = list(lower=m.empt, upper=m.full)

    if(direction=='backward'){
      m.selec = MASS::stepAIC(m.full, direction=direction, k=k, scope=scope, trace=0)
    }
    if(direction=='forward'){
      m.selec = MASS::stepAIC(m.empt, direction=direction, k=k, scope=scope, trace=0)
    }
    if(!direction %in% c('backward','forward')){
      stop('\'direction\' should be one of "backward", "forward"')
    }

    v.select = names(coef(m.selec)[!names(coef(m.selec)) %in% c('(Intercept)')])
    adj1[varY, -varY] = colnames(data)[-varY] %in% v.select
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
