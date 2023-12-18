## Download the data provided by paper "Using Near-infrared reflectance spectroscopy (NIRS) to predict glucobrassicin concentrations in cabbage and brussels sprout leaf tissue" Ilse E. Renner1 and Vincent A Fritz

library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(scales)
library(ggrepel)

## ## 92 sample, 141 time points
data = read.csv("NIRS.csv")
dataXall = data[,4:144] # exclude response

data_b <- function(dataX, data, b){
  
  index_train = sample(1:nrow(dataX),size = 50,replace=F)
  discrete_data_all = t(dataXall[index_train,])  ## Ntime * Nsample
  test_discrete_all = t(dataXall[-index_train,])
  
  # Fresh weight GBS
  trainGBS = data[index_train,"GBS"]
  testGBS = data[-index_train,"GBS"]
  
  # Dry weight GBS
  trainDWgbs = data[index_train,"DWgbs"]
  testDWgbs = data[-index_train,"DWgbs"]
  
  return(list("trainGBS"=trainGBS, "testGBS"=testGBS, "trainDWgbs" = trainDWgbs, "testDWgbs" = testDWgbs, "discrete_data_all"=discrete_data_all,"test_discrete_all"=test_discrete_all))
}



datab=list()
nboot=100

for (b in 1:nboot){
  cat(b)
  datab[[b]]=data_b(dataXall, data, b)
}




##############################
####      FSVR            ####
##############################

source("FSVR.R")


smoothers=c(0.5,1,5,10)
npc = 10
Ks = 1:10
nus = c(0.1, 0.3, 0.5, 0.8)
Cs = seq(0.01, 1, length.out = 5)

kernel = "rbfdot"
# kernel="vanilladot"


rmses = NULL
opt.Ks = NULL
opt.ss = NULL
opt.Cs = NULL
opt.nus = NULL


nboot = 100
for (b in 1:nboot){
  
  cat(b)
  output =  datab[[b]]
  x = t(output$discrete_data_all) ## N*Ntime
  y = output$trainDWgbs
  # y = output$trainGBS
  testx = output$test_discrete_all
  testy = output$testDWgbs
  # testy = output$testGBS
  
  fit = FSVR(x, y, kernel = kernel, Ks, smoothers, Cs = Cs, nus = nus, npc = npc, knots = 35, fold = 5, fit = TRUE)   
  opt.k = fit@optk
  opt.s = fit@opts
  opt.c = fit@optc
  opt.nu = fit@optnu
  
  opt.Ks = c(opt.Ks, opt.k)
  opt.ss = c(opt.ss, opt.s)
  opt.Cs = c(opt.Cs, opt.c)
  opt.nus = c(opt.nus, opt.nu)
  
  predtest = predict(fit, t(testx))
  
  Ntest = ncol(testx)
  error = sqrt(sum((predtest-testy)^2)/Ntest)
  rmses = c(rmses, error)
  
}

boxplot(rmses)

##############################
###      Multi SVR         ###
##############################


multiSVR <- function(data, yname = "GBS", fold=5, Cs=1, nus, iSim, kernel="rbfdot"){
 
  output =  data[[iSim]]
  traindata = t(output$discrete_data_all)
  
  test_discrete = output$test_discrete_all  ## Ntime *N
  
  if (yname == "GBS"){
    trainy = output$trainGBS
    testy = output$testGBS
  }
  
  if (yname == "DWgbs"){
    trainy = output$trainDWgbs
    testy = output$testDWgbs
  }
  
  
  library(kernlab)
  
  ### do cross-validation, tuning parameters: C, smoothing parameter lambda, number of components
  ### split data 
  m = length(trainy)
  vgr <- split(sample(1:m,m),1:fold)
  
  
  cerrors = matrix(0, fold,length(nus)*length(Cs))
  
  for (i in 1:fold) {
    cind <-  unsplit(vgr[-i],factor(rep((1:fold)[-i],unlist(lapply(vgr[-i],length)))))
    
    sr = 0
    indexsr = NULL
    for (j in 1:length(nus)){
      for (s in 1:length(Cs)){
        sr = sr+1
        svmfit.cv <-ksvm(trainy[cind]~as.matrix(traindata[cind,]), cross=0, nu = nus[j], C = Cs[s], kpar="automatic", kernel=kernel , type = "nu-svr")
        
        testest.cv <- predict(svmfit.cv,as.matrix(traindata[vgr[[i]],]))
        cerrors[i, sr] = cerrors[i, sr] + (sum((testest.cv-trainy[vgr[[i]]])^2))/length(vgr[[i]])
        
        indexsr = rbind(indexsr, c(j,s,sr))
      }
      
    } ## end of j
  } ## end of fold
  
  
  cerror = apply(cerrors,2,mean)
  index = which.min(cerror)[1]
  index.nu = indexsr[index, 1]
  index.c = indexsr[index,2]
  opt.nu = nus[index.nu]
  opt.c = Cs[index.c]
  
  svmfit.mult <- ksvm(trainy~traindata, cross=0, nu = opt.nu, C=opt.c, kpar="automatic", kernel=kernel, type = "nu-svr")
  
  testest.mult <- predict(svmfit.mult,t(test_discrete))
  
  Ntest = ncol(test_discrete)
  error = sqrt(sum((testest.mult-testy)^2)/Ntest)
  
  return(list("error" = error, "opt.c" = opt.c, "opt.nu" = opt.nu))
  
}



Cs = seq(0.01, 1, length.out = 5)
nus = c(0.1, 0.3, 0.5, 0.8)
kernel = "rbfdot"
# kernel="vanilladot"

yname = "DWgbs"
# yname = "GBS"

multi=list()
multi.rmses=NULL

multi.optCs=NULL
multi.nus = NULL
for (b in 1:nboot) {
  cat(b)
  multi[[b]] = multiSVR(datab, yname = yname, fold=5, Cs=Cs, nus = nus, b, kernel=kernel)
  output = multi[[b]]
  
  multi.rmses = c(multi.rmses, output$error)
  multi.optCs = c(multi.optCs, output$opt.c)
  multi.nus = c(multi.nus, output$opt.nu)
}

boxplot(multi.rmses)



################################################
###      Functional linear regression        ###
################################################

freg.fit <- function(data, yname = "GBS", b){  
  library(refund)
  
  output =  data[[b]]
  discrete_data = output$discrete_data_all  ## time by N
  test_discrete = output$test_discrete_all ## time by N
  
  if (yname == "GBS"){
    trainy = output$trainGBS
    testy = output$testGBS
  }
  
  if (yname == "DWgbs"){
    trainy = output$trainDWgbs
    testy = output$testDWgbs
  }
  
  
  X  = t(discrete_data)
  
  freg.model <- try(pfr(trainy ~ lf(X)))
  if (is(freg.model, 'try-error')) freg.model <- try(pfr(trainy ~ lf(X), method = "ML"))
  summary(freg.model)
  
  testf = t(test_discrete)
  
  testest.freg <- predict(freg.model, newdata = list(X = testf), type = 'response')
  error = sqrt(sum((testest.freg-testy)^2)/length(testy))
  
  return(list("error" = error))
  
}


# yname = "GBS"
yname = "DWgbs"

freg=list()
freg.rmses = NULL
nboot = 100
for (b in 1:nboot) {
  cat(b)
  
  freg[[b]]=freg.fit(datab,yname, b)
  freg.output = freg[[b]]
  freg.rmses = c(freg.rmses, freg.output$error)
  
}

boxplot(freg.rmses)

