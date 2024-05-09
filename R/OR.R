#' Odds ratio of dichotomous variables
#'
#' @param data A matrix where each column is a variable.
#' @param test a logical indicating whether to perform chi-squared test for contingency tables.
#'
#' @return A symmetric matrix of odds ratios.
#' For each pair of variables, it creates a 2x2 contingency table and then it computes the odds ratio.
#' If \code{test=TRUE} it returns the pvalues for the test computed for each 2x2 table.
#'
#' @examples
#'
#' set.seed(1)
#' x <- sample(c(0,1), 100, replace = T)
#' y <- sample(c(0,1), 100, replace = T)
#' z <- sample(c(0,1), 100, replace = T)
#' df <- cbind(x, y, z)
#'
#' or <- OR(df)
#'
#' @export
#'
OR = function(data, test=FALSE){
  out = matrix(0, ncol=ncol(data), nrow=ncol(data))
  colnames(out) = rownames(out) = colnames(data)
  if(test==F){
    for(j in 1:(ncol(data)-1)){
      for(i in (j+1):ncol(data)) {
        or = OR_tab(data[,j], data[,i])
        out[j,i] = or
      }
    }
    out = out + t(out)
    diag(out) = 1
    return(out)
  }else{
    out_p = out
    for(j in 1:(ncol(data)-1)){
      for(i in (j+1):ncol(data)) {
        out[j,i] = OR_tab(data[,j], data[,i])
        out_p[j,i] = chisq.test(table(data[,j], data[,i]))$p.value
      }
    }
    out = out + t(out)
    diag(out) = 1
    return(list('OR'=out,'pvalue'=out_p))
  }
}

