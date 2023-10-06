### Binary example

source("genData.R")
source("FSVC.R")

### generate data
N = 100
K = 5
Ntime = 50
grids = seq(0.0001,1,length.out=Ntime)

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

bfun <- function(p1,p2){
  class = p1^2 + p2^2-0.9
}

lambdas=c(1,0.6,0.3,0.1,0.05)
noise.sigma = 10
data = genfSVC.PCA(N, K, bfun, lambdas, grids, eigenfunction,noise.sigma)
x = t(data$discrete_data)
y = data$classlabel

## Perform FSVC
smoothers = c(0.5,1,5,10)
Cs = seq(0.01, 1, length.out = 5)
npc = 5
Ks = 1:5

fit = FSVC(x, y, kernel = "rbfdot",  Ks, smoothers, Cs = Cs, npc = 5, knots = 35, fold = 5, fit = TRUE)  
opt.k = fit@optk
opt.s = fit@opts
opt.c = fit@optc

## predicted class labels
predclass = fit@predclass
accuracy = fit@accuracy

## Predict on test data set
Ntest = 10000
testdata = genfSVC.PCA(Ntest, K, bfun, lambdas, grids, eigenfunction,noise.sigma)
newdata = t(testdata$discrete_data)
testclass = testdata$classlabel

predtest = predict(fit, newdata = newdata)
accuracy = sum(predtest == testclass)/Ntest




### Continuous example

source("genData.R")
source("FSVR.R")

### generate data
N = 100
K = 5
Ntime = 50
grids = seq(0.0001,1,length.out=Ntime)

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

bfun <- function(p1,p2){
  y = p1^2 + p2^2-0.9
}

lambdas=c(1,0.6,0.3,0.1,0.05)
noise.sigma = 10
data = genfSVR.PCA(N, K, bfun, lambdas, grids, eigenfunction,noise.sigma)
x = t(data$discrete_data)
y = data$y

## Perform FSVR
smoothers=c(0.5,1,5,10)
npc = 5
Ks = 1:5
nus = c(0.1, 0.3, 0.5, 0.8)
Cs = seq(0.01, 1, length.out = 5)

fit = FSVR(x, y, kernel = "rbfdot", Ks, smoothers, Cs = Cs, nus = nus, npc = 5, knots = 35, fold = 5, fit = TRUE)  
opt.k = fit@optk
opt.s = fit@opts
opt.c = fit@optc
opt.nu = fit@optnu

## predicted class labels
predy = fit@predy
error = sqrt(fit@error)

## Predict on test data set
Ntest = 10000
lambdas=c(1,0.6,0.3,0.1,0.05)
testdata = genfSVR.PCA(Ntest, K, bfun, lambdas, grids, eigenfunction,noise.sigma)
newdata = t(testdata$discrete_data)
testy = testdata$y

predtest = predict(fit, newdata = newdata)
error = sqrt(sum((predtest-testy)^2)/Ntest)
