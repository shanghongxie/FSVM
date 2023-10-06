### Generate simulated data: binary and continuous outcomes

### Binary

####################################################################################
##   Scenario 1 (Settings 1 - 3): y was generated from a function of FPC scores
####################################################################################

genfSVC.PCA <- function(N, K, bfun, lambdas, grids, eigenfunction, noise.sigma){
  
  # N: sample size
  # K: number of FPCs
  # bfun: boundary function
  # lambdas: eigenvalues
  # grids: time grid
  # eigenfunction: eigenfunction
  # noise.sigma: the variance of noise in functional data
  
  ### calculate fpca scores on the training set (\sqrt(lambda_i)*Zi)
  trainscore = t(diag(sqrt(lambdas))%*%matrix(rnorm(K*N),K,N))
  trainscore = data.frame(trainscore)
  discrete_data = t(eigenfunction)%*%t(trainscore) ## Ntime*N
  
  
  ### add noise for each subject
  library(MASS)
  
  if (noise.sigma > 0){
    noise = mvrnorm(N,rep(0,length(grids)),diag(rep(noise.sigma,length(grids))))
    
    
    discrete_data = discrete_data+t(noise)
  }
  
  
  # ## create boundary function
  boundary=bfun(trainscore[,1],trainscore[,2])
  
  ## add noise on y
  y = boundary + rnorm(N)
  prob = exp(y)/(1+exp(y))
  class = ifelse(prob>0.5, 1,-1)
  
  return (list("classlabel"=class,"boundary"=boundary, "PCscore"=trainscore,"discrete_data"=discrete_data, "y"=y, "prob"=prob))
  
}


#########################################################################################
##   Scenario 2 (Setting 4): y was generated from a scalar-on-function logistic regression
#########################################################################################


genflg <- function(N,time, breaks, beta, alpha, noise.sigma){
  
  library(fda)
  basis <- create.bspline.basis(c(0,1),breaks=breaks)
  
  nbasis=basis$nbasis
  coefs=matrix(rnorm(N*nbasis),nbasis,N)
  
  ### functional object
  library(fda)
  temfd=fd(coefs,basis)
  
  
  ## generate discrete data
  discrete_data = eval.fd(time,temfd)  ## each row is a time point, each column is a subject
  
  # ## add noise for each subject on functional data 
  library(MASS)
  noise = mvrnorm(N,rep(0,length(time)),diag(rep(noise.sigma,length(time))))
  discrete_data = discrete_data+t(noise)
  
  inner = inprod(basis,basis)
  
  #### Create class label
  betat = fd(beta,basis)
  ## calculate intergal of x_i(t)*\beta(t)dt
  A = t(coefs)
  
  c = A%*%inner%*%beta
  
  ## add a noise term on y
  c = c+rnorm(length(c))
  
  prob = exp(alpha+c)
  prob = prob/(1+prob)
  
  class = ifelse(prob>0.5,1,-1)
  
  
  return (list("classlabel"=class,"coefs"=coefs,"temfd"=temfd,"discrete_data"=discrete_data, "betat"=betat))
  
}

#########################################################################################
##   Scenario 3 (Setting 5): two groups with their own mean curve functions
#########################################################################################

gendifmean <- function(N, time, meanf1, meanf2, sd){
  
  library(fda.usc)
  
  ## group 1
  discrete_data1 = meanf1(time)+matrix(rnorm(N/2*length(time), sd=sd),length(time),N/2)
  m1 = mean(discrete_data1)
  discrete_data1 = discrete_data1-m1
  
  ## group 2
  discrete_data2 = meanf2(time)+matrix(rnorm(N/2*length(time), sd=sd),length(time),N/2)
  m2 = mean(discrete_data2)
  discrete_data2 = discrete_data2-m2
  
  discrete_data = cbind(discrete_data1,discrete_data2)
  
  ## class label (first half: class 1; second half: class -1)
  class = c(rep(1,N/2),rep(-1,N/2))
  
  
  return (list("classlabel"=class,"meanf1"=meanf1,"meanf2"=meanf2,"discrete_data"=discrete_data))
}


#########################################################################################
##   Scenario 4 (Setting 6): two groups with the same mean curve functions but different
##.                          covariance functions
#########################################################################################

gendifcov <- function(N, K, meanf, lambdas, grids, eigenfunction, noise.sigma){
  
  ### calculate fpca scores on the training set (\sqrt(lambda_i)*Zi)
  trainscore1 = t(diag(sqrt(2*lambdas))%*%matrix(rnorm(K*(N/2)),K,N/2))
  trainscore1 = data.frame(trainscore1)
  
  trainscore2 = t(diag(sqrt(lambdas))%*%matrix(rnorm(K*(N/2)),K,N/2))
  trainscore2 = data.frame(trainscore2)
  
  # colnames(trainscore1) = colnames(trainscore2) = c("PC1","PC2","PC3","PC4","PC5")
  
  discrete_data1 = meanf(grids) + t(eigenfunction)%*%t(trainscore1) ## Ntime*N
  discrete_data2 = meanf(grids) + t(eigenfunction)%*%t(trainscore2) ## Ntime*N
  
  discrete_data = cbind(discrete_data1, discrete_data2)
  
  
  ### add noise for each subject on functional data
  library(MASS)
  
  if (noise.sigma>0){
    noise = mvrnorm(N, rep(0,length(grids)), diag(rep(noise.sigma,length(grids))))
    discrete_data = discrete_data + t(noise)
  }
  
  ## class label (first half: class 1; second half: class -1)
  y = c(rep(1,N/2),rep(-1,N/2))
  
  
  return (list("y" = y, "discrete_data" = discrete_data))
  
}


#########################################################################################
##   Scenario 5 (Setting 7): y generated from a functional linear discriminant model
#########################################################################################

genflda <- function(N, grids, nbasis, norder, beta1, beta2,alpha, noise.sigma){
  
  library(fda)
  basis <- create.bspline.basis(c(0,1), nbasis = nbasis, norder = norder)
  coefs = matrix(rnorm(N*nbasis),nbasis,N)
  
  
  ### functional object
  library(fda)
  temfd = fd(coefs,basis)
  
  
  ## generate discrete data
  discrete_data = eval.fd(grids,temfd)  ## each row is a time point, each column is a subject
  
  ## add noise for each subject on functional data
  library(MASS)
  noise = mvrnorm(N,rep(0,length(grids)),diag(rep(noise.sigma,length(grids))))
  discrete_data = discrete_data+t(noise)
  
  inner=inprod(basis,basis)
  
  #### Create class label
  betat1 = fd(beta1,basis)
  betat2 = fd(beta2,basis)
  
  ## calculate intergal of x_i(t)*\beta(t)dt
  A = t(coefs)
  
  c1 = A%*%inner%*%beta1
  c2 = A%*%inner%*%beta2
  
  ## add a noise term
  noise = rnorm(length(c1))
  c1 = c1 + noise
  c2 = c2 + noise
  
  class = ifelse(c1>c2,1,-1)
  
  
  return (list("classlabel"=class,"coefs"=coefs,"temfd"=temfd,"discrete_data"=discrete_data, "betat1"=betat1, "betat2"=betat2))
  
}


## Continuous outcomes


####################################################################################
##   Scenario 1 (Settings 1 - 3): y was generated from a function of FPC scores
####################################################################################

genfSVR.PCA <- function(N, K, bfun,lambdas,grids,eigenfunction,noise.sigma){
  
  ### calculate fpca scores on the training set (\sqrt(lambda_i)*Zi)
  trainscore=t(diag(sqrt(lambdas))%*%matrix(rnorm(K*N),K,N))
  trainscore=data.frame(trainscore)
  
  discrete_data=t(eigenfunction)%*%t(trainscore) ## Ntime*N
  
  
  ## add noise for each subject on functional data
  library(MASS)
  
  if (noise.sigma>0){
    noise=mvrnorm(N,rep(0,length(grids)),diag(rep(noise.sigma,length(grids))))
    discrete_data=discrete_data+t(noise)
  }

  # ## create y function
  y=bfun(trainscore[,1],trainscore[,2])
 
  
  return (list("y"=y, "PCscore" = trainscore, "discrete_data"=discrete_data))
 
}


#########################################################################################
##   Scenario 2 (Setting 4): y was generated from a scalar-on-function linear regression
#########################################################################################

genflr <- function(N,time, breaks, beta, alpha, noise.sigma){
  
  library(fda)
  basis <- create.bspline.basis(c(0,1),breaks=breaks)
  
  nbasis=basis$nbasis
  coefs=matrix(rnorm(N*nbasis),nbasis,N)
  
  ### functional object
  library(fda)
  temfd=fd(coefs,basis)
  
  
  ## generate discrete data
  discrete_data = eval.fd(time,temfd)  ## each row is a time point, each column is a subject
  
  # ## add noise for each subject on functional data 
  library(MASS)
  noise = mvrnorm(N,rep(0,length(time)),diag(rep(noise.sigma,length(time))))
  discrete_data = discrete_data+t(noise)
  
  inner = inprod(basis,basis)
  
  #### generate y
  betat = fd(beta,basis)
  ## calculate intergal of x_i(t)*\beta(t)dt
  A = t(coefs)
  
  c = A%*%inner%*%beta
  
  ## add a noise term on y
  c = c+rnorm(length(c))
  y = alpha+c
  
  return (list("y"=y,"coefs"=coefs,"temfd"=temfd,"discrete_data"=discrete_data, "betat"=betat))
  
}

#########################################################################################
##   Scenario 3 (Setting 5): two groups with their own mean curve functions
#########################################################################################

gendifmean_cont <- function(N, time, meanf1, meanf2, sd){
  
  library(fda.usc)
  
  discrete_data1 = meanf1(time) + matrix(rnorm(N/2*length(time), sd=sd),length(time),N/2)
  
  m1 = mean(discrete_data1)
  
  discrete_data1 = discrete_data1-m1
  
  discrete_data2 = meanf2(time) + matrix(rnorm(N/2*length(time), sd=sd),length(time),N/2)
  
  m2 = mean(discrete_data2)
  
  discrete_data2 = discrete_data2-m2
  
  discrete_data = cbind(discrete_data1,discrete_data2)
  
  ## class label (first half: class 1; second half: class -1)
  class=c(rep(1,N/2),rep(-1,N/2))
  
  # add an error
  y = class + rnorm(N)
  
  return (list("y"=y,"meanf1"=meanf1,"meanf2"=meanf2,"discrete_data"=discrete_data))
}

#########################################################################################
##   Scenario 4 (Setting 6): two groups with the same mean curve functions but different
##.                          covariance functions
#########################################################################################

gendifcov_cont <- function(N, K, meanf, lambdas, grids, eigenfunction, noise.sigma){
  
  ### calculate fpca scores on the training set (\sqrt(lambda_i)*Zi)
  trainscore1 = t(diag(sqrt(2*lambdas))%*%matrix(rnorm(K*(N/2)),K,N/2))
  trainscore1 = data.frame(trainscore1)
  
  trainscore2 = t(diag(sqrt(lambdas))%*%matrix(rnorm(K*(N/2)),K,N/2))
  trainscore2 = data.frame(trainscore2)
  
  
  discrete_data1 = meanf(grids) + t(eigenfunction)%*%t(trainscore1) ## Ntime*N
  discrete_data2 = meanf(grids) + t(eigenfunction)%*%t(trainscore2) ## Ntime*N
  
  discrete_data = cbind(discrete_data1, discrete_data2)
  
  
  ### add noise for each subject on functional data
  library(MASS)
  
  if (noise.sigma>0){
    noise = mvrnorm(N, rep(0,length(grids)), diag(rep(noise.sigma,length(grids))))
    discrete_data = discrete_data + t(noise)
  }
  
  ## y label (first half: y 1; second half: y -1)
  class = c(rep(1,N/2),rep(-1,N/2))
  
  # add an error
  y = class + rnorm(N)
  
  return (list("y" = y, "discrete_data" = discrete_data))
  
}
