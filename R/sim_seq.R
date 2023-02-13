#' Simulate data
#'
#' @param M is the number of components
#' @param K is the number of Markov model states 
#' @param ini.prob is a list of initial probability vectors for each component
#' @param trans.prob is a list of transition matrices for each component
#' @param clust.size is a list of components' sizes
#' @param T.range is a vector of two elements: minimum and maximum sequence length
#' @export
#' @importFrom stats rmultinom
sim_seq = function( M, K, ini.prob, trans.prob, clust.size, T.range)
{
data.sim<-NULL
num.cluster = M
states = K
ini = ini.prob
A = trans.prob
N.totsim = unlist(clust.size)
# cluster proportions
clus.prop <- N.totsim/(sum(N.totsim))
# cumulative sum of cluster sizes
cumsum.size <- cumsum(c(0,N.totsim))

# simulate sequences in each cluster
 for(cl in 1:num.cluster){
  for(m in 1:N.totsim[cl]){   # subjects
   num.obs<-rdweibt(1, T.range, .48, .55)  # number of observations for each sequence. truncated Discrete Weibull
   y <- rep(0,num.obs)
   y[1]<-sample(1:states,1,prob=ini[[cl]])   
   for (j in 2:num.obs){
     # multinomial draws
     mChoices = t(apply(t(as.matrix(A[[cl]][y[j-1],])), 1, rmultinom, n = 1, size = 1))
     y[j] = apply(mChoices, 1, function(x) which(x==1))
   }  
   id<-rep(m + cumsum.size[cl], num.obs)      # subject ID
   clus<-rep(cl,num.obs)    # cluster ID
   sub<-cbind(id,y,clus)
   data.sim<-as.data.frame(rbind(data.sim,sub)) 
  } 
 } 
data.sim$y <- as.factor(data.sim$y)
return(data.sim) 
 
} 

