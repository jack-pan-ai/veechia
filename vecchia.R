library(fields)
library(ggplot2)
library(corrplot)
library(GpGp)


full_prob_GP = function(cov, x, mu=0){
  k = length(x)
  prob = -0.5 * (log(det(cov)) + (x -mu) %*% solve(cov) %*% (x -mu) + k * log(2 * pi))
  return(prob)
}

conditional_prob_GP = function(cov_matrix, num_, n_partition_each_block, x){
  length1 = n_partition_each_block
  length2 = num_ - n_partition_each_block
  x1 = x[(length2+1):num_]
  x2 = x[1:length2]
  # Sigma partition
  sigma1 = cov_matrix[(length2 + 1):(length2 + length1), (length2 + 1): (length2 + length1)]
  sigma2 = cov_matrix[1:length2, 1:length2]
  sigma12 = cov_matrix[ (length2 + 1): (length2 + length1), 1:length2]
  sigma21 = t(sigma12)
  # new Sigma for x1|x2
  sigma_conditioned = sigma1 - sigma12 %*% solve(sigma2) %*% sigma21
  mu_conditioned = c(sigma12 %*% solve(sigma2) %*% x2)
  # prob calculation 
  prob_ = full_prob_GP(sigma_conditioned, x1, mu_conditioned)
  
  return (prob_)
}


vecchia_loglik = function(covparams, y, locs, n_partition_each_block = 4, covfun_name='matern'){
  cov = matern15_isotropic(covparams, locs)
  # Block partition, it is the "np" where n is # locations; p << 1
  nn = length(y)
  partition = matrix(1:(nn), ncol = n_partition_each_block, byrow = TRUE)
  num_composite = nn / n_partition_each_block
  ## initialization
  prob = full_prob_GP(cov[1:4, 1:4], y[partition[1, ]])
  ## forloop for vecchia (fake parallel)
  for (i in 2:num_composite){
    num_ = i * n_partition_each_block
    prob = prob + conditional_prob_GP(cov_matrix = cov[1: num_, 1: num_], num_ = num_, 
                                      n_partition_each_block = n_partition_each_block,
                                      x=y) 
  }
  return(prob)
}

# setting, n*n matrix 
n = 12
locs = as.matrix(expand.grid((1:n)/n, (1:n)/n))
# models parameters
covparams = c(1, 0.1689, 0)
cov = matern15_isotropic(covparams, locs)

# simulation 
set.seed(42)
xx = rnorm(n*n)
u = chol(cov)
yy = c(t(u) %*% xx)

# vecchia log-likelihood computed (ours)
vecchia_loglik(covparams, yy, locs)
# vecchia log-likelihood (gpgp)
vecchia_meanzero_loglik(covparams, "matern15_isotropic", yy, locs, find_ordered_nn(locs, n*n-1))
# full prob 
full_prob_GP(cov, yy)

  