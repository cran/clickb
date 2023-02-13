#' Marginal likelihood for each Markov model
#'
#' @param A is a transition matrix for the Markov chain
#' @param pi is the vector of initial probabilities for the Markov chain
#' @param data.contr is a subset of the dataset related to one subject
#' @keywords internal
marloglik.contr <- function(A, pi, data.contr) {     
  N.t <- nrow(data.contr)   
  step.t <- array(dim = c(N.t, 1))
  step.t[1,1] <- pi[data.contr$y[1]]   

  for (j in 2:N.t) {
    step.t[j,1] <- step.t[j-1,1]*A[data.contr$y[j-1],data.contr$y[j]]
  }
 
  marloglik.contr <- log(step.t[N.t,1])
  
  return(marloglik.contr)
}