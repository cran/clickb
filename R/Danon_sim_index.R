#' Compute Danon similarity between original and current partition
#'
#' @param data are sequential data whose sequences are assigned to clusters by the algorithm
#' @param label.clus are the original cluster label related to the simulated data
#' @param fir.obs.check is the vector of positions in data for the first elements of each sequence
#' @keywords internal
Danon_sim_index <- function(data, label.clus, fir.obs.check) {
data.sim <-data
Conf.mat <- table(as.factor(data.sim$clus[fir.obs.check]), as.factor(label.clus)  )  # confusion matrix
N_i<-rowSums(Conf.mat) 
N_j<-colSums(Conf.mat)
tot_seq<- sum(Conf.mat)
denomin<-sum( N_i*log(N_i/tot_seq) ) +sum( N_j*log(N_j/tot_seq) )
N_i_N_j<-tcrossprod(as.matrix(N_i),as.matrix(N_j))  
Sim.index <- ( -2*sum( (Conf.mat+0.01)*log(((Conf.mat+0.01)*tot_seq)/N_i_N_j) ) )/denomin 
results <- list(Sim.index, Conf.mat)
names(results) <- c("Sim.index", "Conf.mat")
return(results)
}