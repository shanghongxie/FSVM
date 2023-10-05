
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

- The code for the proposed methodology is included in **ICATemporalNetwork** folder. Please download all the files in the folder to implement the method.
  + The main function for the method is **ICATemporalNet.R** and **ICATemporalNetBoot.R** which allows bootstraps.
  + To use **ICATemporalNet.R**, it requires **estICA.R** which estimates the non-Gaussian signals and then removes them from raw data, and **temporalNet.R** which estimates the temporal network. 


 
- **Examples** folder contains examples.
   + **genData.R**: generate simulated data
   + **example.R**: an example to implement the method
   + **Sim_Scenario1.R**: simulations of Scenario 1
   + **Sim_Scenario2.R**: simulations of Scenario 2

### Main function: ICATemporalNet
#### Arguments
+ `Yts`: input data, the user must supply a list of Yts, where each element is a N*K data matrix at time t. N is sample size, K is the number of nodes.
+ `N`: sample size
+ `Ntime`: total number of time points
+ `ncomp`:  maximum number of independent components to be chosen
+  `Ta`: use t<=Ta time points to estimate temporal network A
+  `Tc`: ues t>Tc time points to estimate contemporaneous network Gamma

#### Value
+ `estIC`: results from non-Gaussian estimation step. output from estICA.R
+ `estRts`: R(t),residuals after removing non-Gaussian signals
+ `estS`: independent components S
+ `estUts`: non-Gaussian signals U(t)
+ `estWt`: weight matrix w(t)
+  `nIC`: number of independent components
+  `A`: temporal network
+  `Gamma`: contemporaneous network
+  `Omega`: covariance matrix of e(t), inverse of Gamma

### The arguments of other functions are described within R files.
 

