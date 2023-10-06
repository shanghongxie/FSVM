## Generate simulated data for binary outcomes

####################################################################################
##   Scenario 1 (Settings 1 - 3): y was generated from a function of FPC scores
####################################################################################
genfSVR.PCA <- function(N, K, bfun,lambdas,grids,eigenfunction,noise.sigma){
  
  # N: sample size
  # K: number of FPCs
  # bfun: boundary function
  # lambdas: eigenvalues
  # grids: time grid
  # eigenfunction: eigenfunction
  # noise.sigma: the variance of noise in functional data
  
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
  y = bfun(trainscore[,1],trainscore[,2])
  
  
  return (list("y"=y, "PCscore" = trainscore, "discrete_data"=discrete_data))
  
}


Ntrain=100
Ntime=50
Ntest=10000
nSim=100


grids=seq(0.0001,1,length.out=Ntime)

## Legendre polynomials

eigen1=(35*(2*grids-1)^4-30*(2*grids-1)^2+3)/8
eigen2=(63*(2*grids-1)^5-70*(2*grids-1)^3+15*(2*grids-1))/8
eigen3=2*grids-1
eigen4=(3*(2*grids-1)^2-1)/2
eigen5=(5*(2*grids-1)^3-3*(2*grids-1))/2

library(pracma)
eigen1=eigen1/sqrt(trapz(grids,eigen1*eigen1))
eigen2=eigen2/sqrt(trapz(grids,eigen2*eigen2))
eigen3=eigen3/sqrt(trapz(grids,eigen3*eigen3))
eigen4=eigen4/sqrt(trapz(grids,eigen4*eigen4))
eigen5=eigen5/sqrt(trapz(grids,eigen5*eigen5))

eigenfunction=rbind(eigen1,eigen2,eigen3,eigen4,eigen5)


## generate training data
data_fsvm <- function(Ntrain, K, bfun, lambdas, grids, eigenfunction,noise.sigma){
  
  train = genfSVR.PCA(Ntrain,K, bfun,lambdas,grids,eigenfunction,noise.sigma)
  trainy = train$y
  PCscores = train$PCscore
  discrete_data = train$discrete_data
  
  return(list("train" = train, "trainy" = trainy, "discrete_data" = discrete_data, "PCscores" = PCscores))
}


sims = 1:3
noise.sigmas = c(1, 10)
lambdas=c(1,0.6,0.3,0.1,0.05)
K = 5


for (sim in sims){
  
  cat (sim)
  
  if (sim == 1){
    ## setting 1: linear
    bfun <- function(p1,p2){
      y = p1+0.3*p2
    }
    
  }
  
  if (sim == 2){
    ## setting 2: interaction and quadratic
    bfun <- function(p1,p2){
      y = p1*p2 + p1^2-0.5
    }
    
  }
  
  if (sim == 3){
    ## setting 3: quadratics
    bfun <- function(p1,p2){
      y = p1^2 + p2^2-0.9
    }
    
  }
  
  for (l in 1:length(noise.sigmas)){
    noise.sigma = noise.sigmas[l]
    nSim =  100
    
    ## generate training data
    set.seed(100)
    dataSim=list()
    for (iSim in 1:nSim){
      cat(iSim)
      dataSim[[iSim]] = data_fsvm(Ntrain, K, bfun, lambdas, grids,eigenfunction,noise.sigma)
    }
    
    
    ## generate test data
    Ntest = 10000
    set.seed(2)
    test = genfSVR.PCA(Ntest, K, bfun,lambdas,grids, eigenfunction,noise.sigma)
    test_discrete = test$discrete_data
    testy = test$y
    testPCscore = test$PCscore
    yprop = sum(testy == 1)/Ntest
    
    testboundary=test$boundary
    
  }  ## end of sigma
  
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




data_flr <- function(Ntrain, time, breaks, beta, alpha, noise.sigma){
  
  train = genflr(Ntrain,time, breaks, beta, alpha, noise.sigma)
  trainy = train$y
  discrete_data = train$discrete_data
  
  return(list("train"=train,"trainy"=trainy,"discrete_data"=discrete_data))
}




alpha=0
Ntrain=100
Ntime=50
Ntest = 10000
time=seq(0.05,1,length.out=Ntime)
breaks=seq(0,1,by=0.5)
beta=c(0.8,-1,0.5,1.5,0.5)

library(fda)
basis <- create.bspline.basis(c(0,1),breaks=breaks)
alpha = 0
noise.sigmas = c(1, 10)

nSim=100
for (l in 1:length(noise.sigmas)){
  set.seed(100)
  dataSim=list()
  for (iSim in 1:nSim){
    cat(iSim)
    dataSim[[iSim]] = data_flr(Ntrain, time, breaks, beta, alpha, noise.sigmas[l])
  }
  
  
  test = genflr(Ntest, time, breaks, beta, alpha, noise.sigmas[l])
  testy = test$y
  test_discrete = test$discrete_data

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


meanf1<-function(t){
  cos(t*10-pi/4)+0.5*sin(t*8-pi/4)
}

meanf2 <- function(t) {
  cos(t*8-pi/4)+0.6*sin(t*10-pi/4)
}

data_difmean_cont <- function(Ntrain,time,meanf1, meanf2, sd){
  
  train = gendifmean_cont(Ntrain,time, meanf1, meanf2, sd)
  trainy = train$y
  discrete_data=train$discrete_data
  
  return(list("train"=train,"trainy"=trainy,"discrete_data"=discrete_data))
  
}

sds = c(1, sqrt(10))

nSim=100

for (l in 1:length(sds)){
  set.seed(100)
  dataSim=list()
  for (iSim in 1:nSim){
    cat(iSim)
    dataSim[[iSim]] = data_difmean_cont(Ntrain, time,meanf1, meanf2, sds[l])
  }
  
  Ntest=10000
  test = gendifmean_cont(Ntest,time, meanf1, meanf2, sds[l])
  testy = test$y
  test_discrete = test$discrete_data
  
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


data_difcov_cont <- function(Ntrain, K, meanf, lambdas, grids,eigenfunction,noise.sigma){
  
  train = gendifcov_cont(Ntrain, K,  meanf, lambdas, grids, eigenfunction, noise.sigma)
  
  trainy = train$y
  
  discrete_data = train$discrete_data
  
  return(list("train" = train, "trainy" = trainy, "discrete_data" = discrete_data))
  
}

## define mean function
meanf <-function(t){
  cos(t*10-pi/4)+0.5*sin(t*8-pi/4)
}

Ntrain = 100
Ntime = 50
Ntest = 10000
time=seq(0.0001,1,length.out=Ntime)

eigen1=(35*(2*time-1)^4-30*(2*time-1)^2+3)/8
eigen2=(63*(2*time-1)^5-70*(2*time-1)^3+15*(2*time-1))/8
eigen3=2*time-1
eigen4=(3*(2*time-1)^2-1)/2
eigen5=(5*(2*time-1)^3-3*(2*time-1))/2


library(pracma)
eigen1=eigen1/sqrt(trapz(time,eigen1*eigen1))
eigen2=eigen2/sqrt(trapz(time,eigen2*eigen2))
eigen3=eigen3/sqrt(trapz(time,eigen3*eigen3))
eigen4=eigen4/sqrt(trapz(time,eigen4*eigen4))
eigen5=eigen5/sqrt(trapz(time,eigen5*eigen5))

eigenfunction = rbind(eigen1,eigen2,eigen3,eigen4,eigen5)


lambdas=c(1,0.6,0.3,0.1,0.05)

noise.sigmas = c(1, 10)

grids = time
K = 5


for (l in 1:length(noise.sigmas)){
  
  noise.sigma = noise.sigmas[l]
  nSim = 100
  set.seed(100)
  dataSim=list()
  for (iSim in 1:nSim){
    cat(iSim)
    dataSim[[iSim]] = data_difcov_cont(Ntrain, K, meanf, lambdas, grids, eigenfunction, noise.sigma)
  }
  
  
  Ntest=10000
  set.seed(2)
  test = gendifcov_cont(Ntest, K, meanf, lambdas, grids, eigenfunction, noise.sigma)
  test_discrete = test$discrete_data
  testy = test$y
  
}


