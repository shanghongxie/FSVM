
# R package for FSVM (Functional Support Vector Machine)

<img src="https://img.shields.io/badge/Study%20Status-Results%20Available-yellow.svg" alt="Study Status: Results Available"> 

Methods to perform support vector machine on functional data for both classification and regression problems.


- Title: **Functional Support Vector Machine**

- Authors: **Shanghong Xie<sup>a,b</sup> (shanghongxie@gmail.com), and R. Todd Ogden<sup>b</sup>**

- Affiliations:
   + 1. **School of Statistics, Southwestern University of Finance and Economics, Chengdu, China**
   + 2. **Department of Biostatistics, Mailman School of Public Health, Columbia University, New York, USA**
  



## Setup Requirements
- R


## Code Instructions

- The code for the proposed methodology is included in **FSVM** folder. Please download all the files in the folder to implement the method.
  + The main function for functional support vector classification (FSVC) is **FSVC.R**.
  + The main function for functional support vector regression (FSVR) is **FSVR.R**.

 
- **Examples** folder contains examples.
   + **genData.R**: functions to generate simulated data for all scenarios
   + **SimData_binary.R**: generate data for binary outcomes in the simulation studies
   + **SimData_continuous.R**: generate data for continuous outcomes in the simulation studies
   + **example.R**: examples to implement the method for both binary and continuous outcomes
   + **Sim_binary.R**: run simulations of binary outcomes. Need first use SimData_binary.R to generate data
   + **Sim_continous.R**: run simulations of continuous outcomes. Need first use SimData_continuous.R to generate data

### Main functions: 
### FSVC
#### Usage
FSVC(x = NULL, y = NULL, kernel = "rbfdot", Ks = NULL, smoothers = NULL, Cs = 1, npc = 5, knots = 35, fold = 5, fit = TRUE)

#### Arguments
+ `x`: functional data matrix, an N*Ntime matrix
+ `y`: class label
+ `kernel`: the kernel function. rbfdot: Gaussian kernel, vanilladot: linear kernel
+ `Ks`: the grid of number of FPCs
+ `smoothers`: the grid of smoothing parameter in FPCA
+ `Cs`: the grid of regularization parameter C in SVC
+ `npc`: the maximum number of FPCs
+ `knots`: number of knots to use or the vectors of knots in fpca.face function; defaults to 35
+ `fold`: number of folds for cross-validation
+ `fit`: whether to predict class labels, default to TRUE

#### Value
An S4 object of class "FSVC" containing the fitted model, Accessor functions can be used to access the slots of the object (see examples) which include:
+ `opts`: optimal smoothing parameter
+ `optc`: optimal C parameter
+ `optk`: optimal K parameter
+ `sig2est`: estimated variance from FPCA
+ `score`: estimated FPC scores
+  `eigenfest`: estimated eigenfunctions
+  `evalest`: estimated eigenvalues
+  `muest`: estimated mean function
+  `npc`: the maximum number of FPCs
+  `svm.fit`: "ksvm" object from ksvm function
+  `predclass`: predicted class labels
+  `accuracy`: classification accuracy on fitted data if fit = TRUE

#### Examples
```
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
```

### FSVR
#### Usage
FSVR(x = NULL, y = NULL, kernel = "rbfdot",  Ks = NULL, smoothers = NULL, Cs = 1, nus = NULL, npc = 5, knots = 35, fold = 5, fit = TRUE)

#### Arguments
+ `x`: functional data matrix, an N*Ntime matrix
+ `y`: class label
+ `kernel`: the kernel function. rbfdot: Gaussian kernel, vanilladot: linear kernel
+ `Ks`: the grid of number of FPCs
+ `smoothers`: the grid of smoothing parameter in FPCA
+ `Cs`: the grid of regularization parameter C in SVR
+ `nus`: the grid of tuning parameter nu in SVR
+ `npc`: the maximum number of FPCs
+ `knots`: number of knots to use or the vectors of knots in fpca.face function; defaults to 35
+ `fold`: number of folds for cross-validation
+ `fit`: whether to predict y, default to TRUE

#### Value
An S4 object of class "FSVR" containing the fitted model, Accessor functions can be used to access the slots of the object (see examples) which include:
+ `optS`: optimal smoothing parameter
+ `optc`: optimal C parameter
+ `optk`: optimal K parameter
+ `optnu`: optimal nu parameter
+ `sig2est`: estimated variance from FPCA
+ `score`: estimated FPC scores
+  `eigenfest`: estimated eigenfunctions
+  `evalest`: estimated eigenvalues
+  `muest`: estimated mean function
+  `npc`: the maximum number of FPCs
+  `svm.fit`: "ksvm" object from ksvm function
+  `predy`: predicted y
+  `error`: mean squared error on fitted data if fit = TRUE

#### Examples
```
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
```

### The arguments of other functions are described within R files.
 

