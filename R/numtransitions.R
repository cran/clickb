#'function to count number of transitions between states in sequential data (I need it to compute posterior probabilities)
#'
#' @param subdata is a subset of original data referred to one subject
#' @keywords internal
numtransitions<-function(subdata){
split.data.prov<- subdata
split.data.from <- split.data.prov[-length(split.data.prov)]
split.data.to <- split.data.prov[-1]
N.trans <- table(split.data.from, split.data.to)
  return(N.trans)
}