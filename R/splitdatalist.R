#' Sequences as list
#'
#' @param data are a dataframe where each row identify the element of a sequence by time
#' @keywords internal
splitdatalist<-function(data){
  data <- data[order(data$id),]
  data.list <- split(data, data$id)
  data.list <- lapply(data.list, function(x) data.frame(x, obs.num = 1:nrow(x)))
  return(data.list)
}