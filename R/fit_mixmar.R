#' Bayesian estimation for mixture of Markov models with fixed number of cluster
#'
#' @param data are dataset containing a variable for categorical observations and an ID variable to identify subjects. No missing values allowed
#' @param burn is the number of burn-in iterations
#' @param iter is the number of iterations
#' @param num.cluster is the number of clusters
#' @param states is the number of states of the Markov chain
#' @param A.ini is the list of initial values for the transition matrices
#' @param pi.ini is the list of initial values for the initial probability vectors
#' @param prior.ini is the hyperparameter for the prior distribution of initial probability vector (Dirichlet)
#' @param prior.transrow is the hyperparameter for the prior distribution of transition matrices rows (Dirichlet)
#' @param prior.mixcoef is the hyperparameter for the prior distribution of mixture coefficients (Dirichlet)
#' @param ini.constr is the vector of zeros and ones indicating constrains in initial probabilities (zero: not possible to start from that state)
#' @param trans.constr is the matrix of zeros and ones indicating constrains in transition probabilities (zero: transition is not allowed)
#' @export
fit_mixmar = function(data, iter = 1000, burn = 500, num.cluster, states, A.ini = NULL, 
                             pi.ini = NULL, prior.ini = 1, prior.transrow = 1, prior.mixcoef = 1, 
                             ini.constr = NULL, trans.constr = NULL){     
id <- NULL
if( any( is.na( data ) ) ) stop( "'data' must contain no missing values" ) 
if( iter <= burn )       stop( "'N.iter' must be higher than 'burnin'" )	

data.list<-splitdatalist(data)     # sequences as list

# First steps (to compute Danon similarity index)
first.obs.check<-c(1,cumsum(as.vector(table(data$id)))+1)   # first step for each subject
first.obs.check<-first.obs.check[-length(first.obs.check)]

#  Constrains. Save positions of zeros as impossible links between states and impossible first states 
zeros<-list()
zeros_0 <- which(ini.constr==0) 
for(jj in 1:states)
       zeros[[jj]]  <-  which(trans.constr[jj,]==0) 

# initial values for transition probabilities
A.s<-list()
if(!(is.null(A.ini))){
 for(cl in 1:num.cluster)
       A.s[[cl]]<-A.ini
}else{                    # initial values have to respect constraints
 trans.p.fix <- 1/states*matrix(rep(1, length = (states^2)), nrow = states, byrow = TRUE, dimnames = list(as.character(c(1:states))) )
 for(st in 1:states)
   trans.p.fix[st,zeros[[st]]] <-0
 for(cl in 1:num.cluster)
       A.s[[cl]]<-trans.p.fix/rowSums(trans.p.fix) 
 }
   
A.list<-list()   # list of transition matrices
A.list[[1]]<-A.s

# initial values for initial probabilities
pi.list<-list()
pi.list[[1]] <- list() 
if(!(is.null(pi.ini))){
for(cl in 1:num.cluster)
       pi.list[[1]][[cl]] <- pi.ini
}else{                                          # initial values have to respect constraints
  pi.list.fix <- rep(1/states, length = states)    
  pi.list.fix[zeros_0]<-0 
  for(cl in 1:num.cluster)
       pi.list[[1]][[cl]] <- pi.list.fix/sum(pi.list.fix)
}

clus.label<-array(0,c(burn+iter,length(unique(data$id))))  # array keeping track of cluster labels

# Assign first cluster labels at random 
clus.label[1,] <- sample(1:num.cluster, length(unique(data$id)),prob=rep(1/num.cluster, length=num.cluster), replace=TRUE)


mix.coeff<-list()    # mixture coefficients for each iteration
mix.coeff[[1]]<-rep(1/num.cluster, num.cluster) 


# save the Danon similarity index and Rand index
Sim.index <- Sim.index.r <- vector(mode="numeric", length=iter)
Sim.index[1]<- Sim.index.r[1]<- 0

# Save new values
prov.A.new <- list()
prov.pi.new <- list()

# start sampler
for (i in 2:(burn+iter)){

  A.save.list<-list()
  pi.save.list<-list()

  for (cl in 1:num.cluster){    # for each cluster

      # parameters related to cluster cl
    subj.cl<-which(clus.label[i-1,]==cl)  # subjects in cluster cl
    A.cl<-A.list[[i-1]][[cl]]        
    pi.cl<-pi.list[[i-1]][[cl]]      
    
    data.list.cl <- data.list[subj.cl]  # current cluster cl sequences
    data.list2 <- lapply(data.list.cl, function(x) x[,-c(1,3,4)])   # remove other variables (ID)
    
    # Update
    
    ###update initial 
    # cumulate sequence lengths for each subject
    data.cl <- subset(data, id %in% subj.cl) 
    first.obs<-c(1,cumsum(as.vector(table(data.cl$id)))+1)     # initial step position
    first.obs<-first.obs[1:(length(first.obs)-1)] 

    pi.cl.prov<-rep(0, states)
    pi.cl.prov[-zeros_0]<- MCMCpack::rdirichlet(1,1+prior.ini+table(data.cl$y[first.obs])[-zeros_0])  
    pi.new <- pi.cl.prov
    
    ###update transition
    N.trans <- Reduce('+', lapply(data.list2, numtransitions))
    A.new <- matrix(0, nrow=states, ncol=states)
    for(k in 1:states){     
     A.prov <- rep(0, states)
     A.prov[-zeros[[k]]] <- MCMCpack::rdirichlet(1,1+prior.transrow+ N.trans[k,-zeros[[k]]] )
     A.new[k,]<- A.prov
    }
    

    # save current values
    prov.A.new[[cl]]  <- A.new
    prov.pi.new[[cl]]  <- pi.new

  }
    
  
  mix.coeff.vec<-mix.coeff[[i-1]]    # previous iteration mixture coefficients
    
  ###  Save parameters 
  A.save.list <- prov.A.new
  pi.save.list <- prov.pi.new
  
  ###update cluster membership
          # marginal loglik 
      marginloglik <- parallel::mclapply(1:length(data.list), function (x) { lik<-NULL
      for (jjj in 1:num.cluster) {
           lik<-c(lik, marloglik.contr(A=A.save.list[[jjj]], pi=pi.save.list[[jjj]], data.contr=data.list[[x]]))
           } 
      return(lik)})
      
     # assign subjects in new clusters       
  post.prob<-do.call(rbind,parallel::mclapply(marginloglik, function(x) exp(x-max(x))))
  clus.label[i,]<-apply(post.prob,1,function(t){
    sample(1:num.cluster,1,prob = t)})
    
  
  ###update mixture coefficients
  mix.coeff[[i]]<-MCMCpack::rdirichlet(1,table(clus.label[i,])+1+prior.mixcoef)  
  
    #####update parameters associated with each cluster
  A.up.list<-list()
  pi.up.list<-list()
  
    # update parameters for each cluster after subjects are reassigned according to posterior probabilities
  for (jj in 1:num.cluster){
    cl.indx<-which(clus.label[i,]==jj)
    ###update initial 
    data.cl <- subset(data, id %in% cl.indx) 
    first.obs<-c(1,cumsum(as.vector(table(data.cl$id)))+1)  
    first.obs<-first.obs[1:(length(first.obs)-1)] 

    pi.cl.prov<-rep(0, states)
    pi.cl.prov[-zeros_0]<- MCMCpack::rdirichlet(1,1+prior.ini+table(data.cl$y[first.obs])[-zeros_0]) 
    pi.up.list[[jj]] <- pi.cl.prov 
    
    ###update transition 
    data.list.cl <- data.list[cl.indx]  # sequences in new clusters
    data.list2 <- lapply(data.list.cl, function(x) x[,-c(1,3,4)])   
    N.trans <- Reduce('+', lapply(data.list2, numtransitions))
    A.post <- matrix(,nrow=states,ncol=states, dimnames = list(as.character(c(1:states))) )
    for(k in 1:states){    
     A.prov2 <- rep(0, states)
     A.prov2[-zeros[[k]]] <- MCMCpack::rdirichlet(1,1+prior.transrow+ N.trans[k,-zeros[[k]]] )
     A.post[k,] <- A.prov2
    }
     
    A.up.list[[jj]] <- A.post
  }  
  pi.list[[i]]<-pi.up.list        
  A.list[[i]]<-A.up.list
  
# Compute Danon similarity between original and current partition
Danon_res <- Danon_sim_index(data, clus.label[i-1,], first.obs.check)
Sim.index[i] <- Danon_res$Sim.index
# Compute Rand similarity between original and current partition
Sim.index.r[i] <- mclust::adjustedRandIndex(as.factor(data$clus[first.obs.check]), as.factor(clus.label[i-1,])  ) 
# Save the last Confusion matrix
Conf.mat <- Danon_res$Conf.mat


}

plot(c(1:(i-1)), Sim.index[1:(i-1)], type = "l", lty = 1, xlab="iteration", ylab="Similarity index" )

results <- list(pi.list, A.list, Sim.index, Sim.index.r, Conf.mat)
names(results) <- c("pi.list", "A.list", "Sim.index.Danon", "Sim.index.Rand", "Conf.mat")

return(results)

}
