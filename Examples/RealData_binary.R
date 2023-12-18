
### Download EEG dataset from UCI machine learning repository https://archive.ics.uci.edu/dataset/121/eeg+database
### use the full data set and unzip each subject's folder.


##############################
###      Generate data    ####
##############################

library(ggplot2)
library(dplyr)
library(gridExtra)
library(scales)


list = list.files(path = ".", pattern=NULL, all.files=FALSE,
                  full.names=FALSE)

case_list = c(list[1:65], list[110:121])
control_list = c(list[66:109], list[122])

length(case_list)  ## 77
length(control_list)  ## 45


econame = "P1"

# econame = "P5"

# Case group data
case_data_S2diff = NULL

library(stringr)

seed = 100
set.seed(seed)
for (j in 1:length(case_list)){
  
  cat(j)
  
  filenames = list.files(path = paste(case_list[j],"/", sep = ""), pattern=NULL, all.files=FALSE,
                         full.names=FALSE)
  
  pls = NULL
  for (i in 1:length(filenames)){
    
    my = file(description = paste0(case_list[j],"/", filenames[i]) , open="r", blocking = TRUE)
    pl = readLines(my, n = 4)
    pl = pl[4]
    pl = str_split(pl, ",")[[1]][1]
    
    pls[i] = pl
  }
  
  index_s1 = which(pls == "# S1 obj ")
  index_s2nomatch = which(pls == "# S2 nomatch")
  index_s2match = which(pls == "# S2 match ")
  
  ### aggregate trials with S1 condition (single object); trials with S2 no match; trials with S2 match.
  
  s2match_trail = index_s2match
  s2match_random = sample(s2match_trail, 1)
  
  
  eeg_S2match = read.table(paste0(case_list[j],"/", filenames[s2match_random]))
  
  s2nomatch_trail = index_s2nomatch
  s2nomatch_random = sample(s2nomatch_trail, 1)
  
  eeg_S2nomatch = read.table(paste0(case_list[j],"/", filenames[s2nomatch_random]))
  
  colnames(eeg_S2match) = colnames(eeg_S2nomatch) = c("trail", "position", "Hz", "value")
  
  
  ### S2 match
  eeg_j_S2match = eeg_S2match%>% filter(position == econame)
  eeg_j_wide_S2match = reshape(data = eeg_j_S2match, idvar = "position", v.names = "value", timevar = "Hz", direction = "wide")
  
  eeg_j_wide_S2match = data.frame(ID = j, eeg_j_wide_S2match)
  
  eeg_j_wide_S2match = eeg_j_wide_S2match[,-c(1:3)]
  
  
  ### S2 no match
  eeg_j_S2nomatch = eeg_S2nomatch%>% filter(position == econame)
  eeg_j_wide_S2nomatch = reshape(data = eeg_j_S2nomatch, idvar = "position", v.names = "value", timevar = "Hz", direction = "wide")
  
  eeg_j_wide_S2nomatch = data.frame(ID = j, eeg_j_wide_S2nomatch)
  
  eeg_j_wide_S2nomatch = eeg_j_wide_S2nomatch[,-c(1:3)]
  
  eeg_j_wide_S2diff = -(eeg_j_wide_S2match-eeg_j_wide_S2nomatch)
  
  case_data_S2diff = rbind(case_data_S2diff, eeg_j_wide_S2diff)

}

case_data_S2diff = t(case_data_S2diff) ## Ntime*N


# Control group data
control_data_S2diff = NULL
for (j in 1:length(control_list)){
  
  cat(j)
  
  filenames = list.files(path = paste(control_list[j],"/", sep = ""), pattern=NULL, all.files=FALSE,
                         full.names=FALSE)
  
  pls = NULL
  for (i in 1:length(filenames)){
    
    my = file(description = paste0(control_list[j],"/", filenames[i]) , open="r", blocking = TRUE)
    pl = readLines(my, n = 4)
    pl = pl[4]
    pl = str_split(pl, ",")[[1]][1]
    
    pls[i] = pl
  }
  
  index_s1 = which(pls == "# S1 obj ")
  index_s2nomatch = which(pls == "# S2 nomatch")
  index_s2match = which(pls == "# S2 match ")
  
  ### aggregate trials with S1 condition (single object); trials with S2 no match; trials with S2 match.
  s2match_trail = index_s2match
  s2match_random = sample(s2match_trail, 1)
  
  eeg_S2match = read.table(paste0(control_list[j],"/", filenames[s2match_random]))
  
  s2nomatch_trail = index_s2nomatch
  s2nomatch_random = sample(s2nomatch_trail, 1)
  
  eeg_S2nomatch = read.table(paste0(control_list[j],"/", filenames[s2nomatch_random]))
  
  colnames(eeg_S2match) = colnames(eeg_S2nomatch) = c("trail", "position", "Hz", "value")
  
  
  ### S2 match
  eeg_j_S2match = eeg_S2match%>% filter(position == econame)
  eeg_j_wide_S2match = reshape(data = eeg_j_S2match, idvar = "position", v.names = "value", timevar = "Hz", direction = "wide")
  
  eeg_j_wide_S2match = data.frame(ID = j, eeg_j_wide_S2match)
  
  eeg_j_wide_S2match = eeg_j_wide_S2match[,-c(1:3)]
  
  ### S2 no match
  eeg_j_S2nomatch = eeg_S2nomatch%>% filter(position == econame)
  eeg_j_wide_S2nomatch = reshape(data = eeg_j_S2nomatch, idvar = "position", v.names = "value", timevar = "Hz", direction = "wide")
  
  eeg_j_wide_S2nomatch = data.frame(ID = j, eeg_j_wide_S2nomatch)
  
  eeg_j_wide_S2nomatch = eeg_j_wide_S2nomatch[,-c(1:3)]
  
  eeg_j_wide_S2diff = -(eeg_j_wide_S2match - eeg_j_wide_S2nomatch)
  
  control_data_S2diff = rbind(control_data_S2diff, eeg_j_wide_S2diff)
  
}


control_data_S2diff = t(control_data_S2diff) ## Ntime*N

case_label = rep(1,ncol(case_data_S2diff))
control_label = rep(0, ncol(control_data_S2diff))

label = c(rep(1,ncol(case_data_S2diff)), rep(0, ncol(control_data_S2diff)))


### generate training and test sets
generate_data <- function(case_data_S2diff,  control_data_S2diff, case_label,  control_label, ratio = 0.8){
  index_case = sample(1:ncol(case_data_S2diff), round(ncol(case_data_S2diff)*0.8))
  
  discrete_data_case_S2diff = case_data_S2diff[,index_case]
  
  train_case_class = case_label[index_case]
  
  test_discrete_case_S2diff = case_data_S2diff[,-index_case]
  
  test_case_class = case_label[-index_case]
  
  
  index_control = sample(1:ncol(control_data_S2diff), round(ncol(control_data_S2diff)*0.8))
  discrete_data_control_S2diff = control_data_S2diff[,index_control]
  
  
  train_control_class = control_label[index_control]
  
  test_discrete_control_S2diff = control_data_S2diff[,-index_control]
  
  test_control_class = control_label[-index_control]
  
  
  discrete_data_S2diff = cbind(discrete_data_case_S2diff, discrete_data_control_S2diff)
  
  test_discrete_S2diff = cbind(test_discrete_case_S2diff, test_discrete_control_S2diff)
  
  trainclass = c(train_case_class, train_control_class)
  testclass = c(test_case_class, test_control_class)
  
  return(list("trainclass" = trainclass, "discrete_data_S2diff" = discrete_data_S2diff, "testclass" = testclass, 
              "test_discrete_S2diff" = test_discrete_S2diff))
  
}

### 100 repeats
datab=list()
nboot=100
set.seed(100)
for (b in 1:nboot){
  cat(b)
  set.seed(b)
  datab[[b]] = generate_data(case_data_S2diff, control_data_S2diff, case_label,  control_label, ratio = 0.8)
}



##############################
###      FSVC            ####
##############################

source("FSVC.R")

smoothers = c(0.5,1,5,10)
Cs = seq(0.01, 1, length.out = 5)
npc = 10
Ks = 1:10
accuracys = NULL
opt.Ks = NULL
opt.ss = NULL
opt.Cs = NULL
AUCs = NULL


nboot = 100

set.seed(100)

kernel = "rbfdot"
# kernel="vanilladot"

library(pracma)
for (b in 1:nboot){
  
  cat(b)
  output =  datab[[b]]
  x = t(output$discrete_data_S2diff) ## N*Ntime
  y = output$trainclass
  
  test_discrete_S2diff = output$test_discrete_S2diff
  testclass = output$testclass
  
  fit = FSVC(x, y, kernel = kernel, Ks, smoothers, Cs = Cs,  npc = 10, knots = 15, fold = 5, fit = TRUE)  
  opt.k = fit@optk
  opt.s = fit@opts
  opt.c = fit@optc
  
  opt.Ks = c(opt.Ks, opt.k)
  opt.ss = c(opt.ss, opt.s)
  opt.Cs = c(opt.Cs, opt.c)
  
  predresult= predict(fit, t(test_discrete_S2diff), type = "probabilities")
  predtest = (predresult[,2]>0.5)*1
  
  Ntest = ncol(test_discrete_S2diff)
  accuracy = sum(predtest == testclass)/Ntest
  accuracys = c(accuracys, accuracy)
  
  ### calculate AUC
  TPs=NULL
  FPs=NULL
  TNs=NULL
  FNs=NULL
  
  cutvalue = seq(0.1, 0.9, length.out = 5)
  for (cut in cutvalue){
    predcut = (predresult[,2] > cut)*1
    tp = sum((predcut == 1) & (testclass == 1))/Ntest
    fp = sum((predcut == 1) & (testclass == 0))/Ntest
    tn = sum((predcut == 0) & (testclass == 0))/Ntest
    fn = sum((predcut == 0) & (testclass == 1))/Ntest
    
    TPs=c(TPs,tp)
    
    FPs=c(FPs,fp)
    
    TNs=c(TNs,tn)
    
    FNs=c(FNs,fn)
    
  }
  TPRs=TPs/(TPs+FNs)
  
  FPRs=FPs/(TNs+FPs)
  
  fpr=FPRs
  fpr=sort(fpr)
  fpr=c(0,fpr,1)
  
  tpr=sort(TPRs)
  tpr=c(0,tpr,1)
  
  auc=trapz(fpr,tpr)
  
  
  AUCs=c(AUCs,auc)
  
}

median(accuracys)
mean(accuracys)
boxplot(accuracys)

median(AUCs)
mean(AUCs)
boxplot(AUCs)



##############################
###      Multi SVC         ###
##############################


## Multivaraite approach: Treat each observed time point as a variable 
multi.CV <- function(data, fold=5, Cs=1, iSim, kernel="rbfdot"){
  
  output =  data[[iSim]]
  
  traindata_S2diff = t(output$discrete_data_S2diff) ## N*Ntime
  
  traindata = cbind(traindata_S2diff)
  
  test_discrete_S2diff = output$test_discrete_S2diff
  
  
  test_discrete = rbind(test_discrete_S2diff)
  
  
  Ntest = ncol(test_discrete_S2diff)
  
  trainclass = output$trainclass
  testclass = output$testclass
  
  library(kernlab)
  
  ### do cross-validation, tuning parameters: C, smoothing parameter lambda, number of components
  ### split data 
  m = length(trainclass)
  vgr <- split(sample(1:m,m),1:fold)
  
  # sr = 1
  cerrors = rep(0,length(Cs))
  
  
  for (c in 1:length(Cs)){
    for (i in 1:fold) {
      
      cind <-  unsplit(vgr[-i],factor(rep((1:fold)[-i],unlist(lapply(vgr[-i],length)))))
      
      svmfit.cv <-ksvm(as.factor(trainclass[cind])~as.matrix(traindata[cind,]),cross=0,C=Cs[c],kpar="automatic", kernel=kernel)
      
      
      testest.cv <- predict(svmfit.cv,as.matrix(traindata[vgr[[i]],]))
      cerrors[c] = cerrors[c] + (sum(testest.cv==trainclass[vgr[[i]]]))/length(vgr[[i]])/fold
    } ## end of fold
  } ## end of C
  
  opt.c = Cs[which.max(cerrors)]
  
  svmfit.mult <-ksvm(as.factor(trainclass)~traindata,cross=0,C=opt.c,kpar="automatic",kernel=kernel, prob.model = TRUE)
  
  testest.mult <- predict(svmfit.mult,t(test_discrete), type = "probabilities")
  
  predclass.mult = (testest.mult[,2]>0.5)*1
  
  result = table(predclass.mult,testclass)
  
  accuracy = sum(predclass.mult==testclass)/Ntest
  
  return(list("result" = result, "accuracy" = accuracy, "opt.c" = opt.c, "testest.mult" = testest.mult))
  
}




kernel="rbfdot"

# kernel="vanilladot"

set.seed(100)
multi.accuracy = NULL
multi.classification = matrix(0,2,2)
opt.Cs = NULL
multi = list()
multi.AUCs = NULL

Cs = seq(0.01, 1, length.out = 5)
for (b in 1:nboot) {
  cat(b)
 
  multi[[b]] = multi.CV(datab, fold=5, Cs=Cs, b, kernel= kernel)
  multi.output = multi[[b]]
  multi.accuracy = c(multi.accuracy, multi.output$accuracy)
  opt.Cs = c(opt.Cs, multi.output$opt.c)
  
  testclass = datab[[b]]$testclass
  Ntest = length(testclass)
  
  ### calculate AUC
  
  TPs=NULL
  FPs=NULL
  TNs=NULL
  FNs=NULL
  
  predtest = multi[[b]]$testest.mult[,2]
  cutvalue = seq(0.1, 0.9, length.out = 5)
  for (cut in cutvalue){
    predcut = (predtest > cut)*1
    tp = sum((predcut == 1) & (testclass == 1))/Ntest
    fp = sum((predcut == 1) & (testclass == 0))/Ntest
    tn = sum((predcut == 0) & (testclass == 0))/Ntest
    fn = sum((predcut == 0) & (testclass == 1))/Ntest
    
    TPs=c(TPs,tp)
    
    FPs=c(FPs,fp)
    
    TNs=c(TNs,tn)
    
    FNs=c(FNs,fn)
    
  }
  TPRs=TPs/(TPs+FNs)
  
  FPRs=FPs/(TNs+FPs)
  
  fpr=FPRs
  fpr=sort(fpr)
  fpr=c(0,fpr,1)
  
  tpr=sort(TPRs)
  tpr=c(0,tpr,1)
  
  auc=trapz(fpr,tpr)
  
  
  multi.AUCs = c(multi.AUCs,auc)
  
}

median(multi.accuracy)
mean(multi.accuracy)
boxplot(multi.accuracy)

median(multi.AUCs)
mean(multi.AUCs)
boxplot(multi.AUCs)

##################################################
###     Functional logistic regression         ###
##################################################

flogit.fit <- function(data,testclass, iSim){  
  library(refund)
  
  output =  data[[iSim]]
  trainclass = output$trainclass
  trainclass[trainclass==-1] = 0
  trainclass = as.numeric(trainclass)
  discrete_data = output$discrete_data_S2diff  ## time by N
  
  testclass = output$testclass
  test_discrete = output$test_discrete_S2diff
  
  Ntest = ncol(test_discrete)
  
  X = t(discrete_data) # N by time 
  
  flogit.model <- try(pfr(trainclass ~ lf(X), family = binomial))
  if (is(flogit.model, 'try-error')) flogit.model <- try(pfr(trainclass ~ lf(X), family = binomial, method = "ML"))

  testclass[testclass==-1]=0
 
  testf = t(test_discrete)
  
  testest.flogit <- predict(flogit.model, newdata = list(X = testf), type = 'response')
  
  predclass = 1*(testest.flogit > 0.5)
  
  predclass = as.numeric(predclass)

  result = table(predclass,testclass)
  
  accuracy = sum(predclass==testclass)/Ntest
  
  return(list("result" = result, "accuracy" = accuracy, "testest.flogit" = testest.flogit))
  
}

nboot = 100

flogit=list()
flogit.accuracy = NULL
flogit.classification = matrix(0,2,2)
logit.AUCs = NULL

for (b in 1:nboot) {
  cat(b)
  
  flogit[[b]]=flogit.fit(datab, testclass,b)
  flogit.output = flogit[[b]]
  flogit.accuracy = c(flogit.accuracy, flogit.output$accuracy)
  
  testclass = datab[[b]]$testclass
  Ntest = length(testclass)
  
  ### calculate AUC
  
  TPs=NULL
  FPs=NULL
  TNs=NULL
  FNs=NULL
  
  predtest = flogit[[b]]$testest.flogit
  cutvalue = seq(0.1, 0.9, length.out = 5)

  for (cut in cutvalue){
    predcut = (predtest > cut)*1
    tp = sum((predcut == 1) & (testclass == 1))/Ntest
    fp = sum((predcut == 1) & (testclass == 0))/Ntest
    tn = sum((predcut == 0) & (testclass == 0))/Ntest
    fn = sum((predcut == 0) & (testclass == 1))/Ntest
    
    TPs=c(TPs,tp)
    
    FPs=c(FPs,fp)
    
    TNs=c(TNs,tn)
    
    FNs=c(FNs,fn)
    
  }
  TPRs=TPs/(TPs+FNs)
  
  FPRs=FPs/(TNs+FPs)
  
  fpr=FPRs
  fpr=sort(fpr)
  fpr=c(0,fpr,1)
  
  tpr=sort(TPRs)
  tpr=c(0,tpr,1)
  
  auc=trapz(fpr,tpr)
  
  
  logit.AUCs = c(logit.AUCs,auc)
  
}
flogit.classification=flogit.classification/nboot

median(flogit.accuracy)
boxplot(flogit.accuracy)

mean(logit.AUCs)
median(logit.AUCs)
boxplot(logit.AUCs)


##################################################
###  Functional linear discriminant analysis   ###
##################################################

estPCscore <- function(Y, npc, mu, sigma2, evalues, efunctions){

  Y.pred = Y
  D = NCOL(Y)
  I = NROW(Y)
  I.pred = NROW(Y.pred)
  
  D.inv = diag(1/evalues, nrow = npc, ncol = npc)
  Z = efunctions
  
  scores = matrix(NA, nrow = I.pred, ncol = npc)
  
  Y.tilde = Y.pred - matrix(mu, I.pred, D, byrow = TRUE)
  
  for (i.subj in 1:I.pred) {
    obs.points = which(!is.na(Y.pred[i.subj, ]))
    
    Zcur = matrix(Z[obs.points, ], nrow = length(obs.points), 
                  ncol = dim(Z)[2])
    ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv)
    scores[i.subj, ] = ZtZ_sD.inv %*% t(Zcur) %*% (Y.tilde[i.subj, 
                                                           obs.points])
  }
  
  return(scores)
}

flda.fit <- function(data,iSim){  
  
  library(refund)
  
  output =  data[[iSim]]
  trainclass = output$trainclass
  trainclass = as.vector(trainclass)
  trainclass[trainclass==-1] = 0
  discrete_data = output$discrete_data  ## time by N
  
  testclass = output$testclass
  
  test_discrete_S2diff = output$test_discrete_S2diff
  
  test_discrete = test_discrete_S2diff
  Ntest = ncol(test_discrete)
  
  
  library(fda)
  pca.fit <- fpca.face(Y = t(discrete_data),  simul = F, var = TRUE)
  
  sig2est.cv = pca.fit$sigma2
  score.cv = pca.fit$scores
  eigenfest.cv = pca.fit$efunctions
  evalest.cv = pca.fit$evalues
  muest.cv = pca.fit$mu
  
  npc = ncol(score.cv)
  score.cv = data.frame(score.cv)
  
  model <- lda(trainclass ~ ., data = score.cv)
  
  testclass[testclass==-1]=0
  
  testf = t(test_discrete)
  
  # calculate PC scores on testing in CV
  estscore.cv =  estPCscore(testf, npc = npc, mu=muest.cv, sigma2=sig2est.cv, evalues=evalest.cv, eigenfest.cv)
  
  estscore.cv = data.frame(estscore.cv)
  
  testest.flda = predict(model, estscore.cv)
  
  pred_class = testest.flda$class
  
  posterior = testest.flda$posterior[,2]
  
  result = table(pred_class,testclass)
  
  accuracy = sum(pred_class==testclass)/Ntest
  
  return(list("result" = result, "accuracy" = accuracy, "posterior"= posterior))
}

library(pracma)
flda = list()
flda.accuracy = NULL
flda.classification = matrix(0,2,2)
flda.AUCs = NULL


for (b in 1:nboot) {
  cat(b)
  
  flda[[b]]=flda.fit(datab,b)
  flda.output = flda[[b]]
  flda.accuracy = c(flda.accuracy, flda.output$accuracy)
  
  testclass = datab[[b]]$testclass
  Ntest = length(testclass)
  
  ### calculate AUC
  
  TPs=NULL
  FPs=NULL
  TNs=NULL
  FNs=NULL
  
  predtest = flda[[b]]$posterior
  cutvalue = seq(0.1, 0.9, length.out = 5)
  for (cut in cutvalue){
    predcut = (predtest > cut)*1
    tp = sum((predcut == 1) & (testclass == 1))/Ntest
    fp = sum((predcut == 1) & (testclass == 0))/Ntest
    tn = sum((predcut == 0) & (testclass == 0))/Ntest
    fn = sum((predcut == 0) & (testclass == 1))/Ntest
    
    TPs=c(TPs,tp)
    
    FPs=c(FPs,fp)
    
    TNs=c(TNs,tn)
    
    FNs=c(FNs,fn)
    
  }
  TPRs=TPs/(TPs+FNs)
  
  FPRs=FPs/(TNs+FPs)
  
  fpr=FPRs
  fpr=sort(fpr)
  fpr=c(0,fpr,1)
  
  tpr=sort(TPRs)
  tpr=c(0,tpr,1)
  
  auc=trapz(fpr,tpr)
  
  
  flda.AUCs = c(flda.AUCs,auc)
}

mean(flda.accuracy)
boxplot(flda.accuracy)

mean(flda.AUCs)
median(flda.AUCs)
boxplot(flda.AUCs)
