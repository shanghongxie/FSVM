
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

### FSVR
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
+ `fit`: whether to predict class labels, default to TRUE

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


### The arguments of other functions are described within R files.
 

