## Estimate FPC scores
estPCscore <- function(Y, npc, mu, sigma2, evalues, efunctions){
  
  # Y: functional data matrix, an N*Ntime matrix
  # npc: number of maximum FPCs
  # mu: mean function
  # sigma2: variance
  # evalues: eigenvalues
  # efunctions: eigenfunctions
  
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


SVM <- function(score, class,  K, fold=5, C=1:10, kpar,kernel = "rbfdot"){  
  # score: estimated FPC scores from FPCA
  # class: training class labels
  # K: number of FPCs
  # fold: number of folds for cross-validation
  # C: regularization parameter in SVM
  # kpar: the list of hyper-parameters (kernel parameters) in ksvm function
  # kernel: the kernel function in ksvm function
  
  library(kernlab)
  
  if (kernel == "vanilladot"){
    svmfit <-ksvm(as.factor(class)~as.matrix(score[,1:K]),cross = fold, C = C,kernel ="vanilladot", prob.model = TRUE)
  }else{
    svmfit <-ksvm(as.factor(class)~as.matrix(score[,1:K]),cross = fold, C = C, kpar = kpar, prob.model = TRUE)
  }
  return(list("svmfit"=svmfit))
  
}

## Perform functional support vector classification
setClass("FSVC", slots = list(optla = "numeric", optc = "numeric",
         optk = "numeric", sig2est = "numeric",
         score = "matrix", eigenfest = "matrix",
         evalest = "vector", muest = "vector",  npc = "numeric", svm.fit = "ANY", predclass = "ANY", accuracy = "numeric")
)

setGeneric("FSVC", function(x, ...) standardGeneric("FSVC"))
setMethod("FSVC",signature(x = "matrix"),
FSVC <- function(x, y, kernel = "rbfdot", Cs = 1, Ks, lambdas, npc = 5, knots = 35, fold = 10, fit = TRUE){
  
  # x: functional data matrix, an N*Ntime matrix
  # y: class label
  # kernel: the kernel function
  # Cs: the grid of regularization parameter C in SVM
  # Ks: the grid of number of FPCs
  # lambdas: the grid of smoothing parameter in FPCA
  # npc: the maximum number of FPCs
  # knots: number of knots to use or the vectors of knots in fpca.face function; defaults to 35
  # fold: number of folds for cross-validation
  # fit: whether to predict class labels.
  
  library(refund)
  
  ret <- new("FSVC")
  
  ### do cross-validation, tuning parameters: C, smoothing parameter lambda, number of components
  ### split data 
  m = length(y)
  vgr <- split(sample(1:m,m),1:fold)
  
  sr = 1
  accuracys = matrix(0,length(lambdas),length(Cs))
  indexk = NULL
  indexsr = NULL
  for (la in 1:length(lambdas)){
    for (c in 1:length(Cs)){
      accuracy = rep(0, length(Ks))
      
      for (i in 1:fold) {
        
        cind <-  unsplit(vgr[-i],factor(rep((1:fold)[-i],unlist(lapply(vgr[-i],length)))))
        
        pca.fit <- fpca.face(Y = x[cind,], var = TRUE, simul = F, npc = npc, lambda = lambdas[la], knots = knots)
        
        sig2est.cv = pca.fit$sigma2
        score.cv = pca.fit$scores
        eigenfest.cv = pca.fit$efunctions
        evalest.cv = pca.fit$evalues
        muest.cv = pca.fit$mu
        
        ## calculate PC scores on testing in CV
        estscore.cv =  estPCscore(x[vgr[[i]],], npc = npc, mu = muest.cv, sigma2 = sig2est.cv, evalues = evalest.cv, eigenfest.cv)
        
        for (k in 1:length(Ks)){
          
          fSVM.fit<-SVM(score.cv, y[cind], Ks[k], fold = 0, C=Cs[c], kpar="automatic",kernel= kernel)
          
          svm.fit <- fSVM.fit$svmfit
          
          pred.cv <- predict(svm.fit,as.matrix(estscore.cv[,1:Ks[k]]))
          
          ## classification accuray
          accuracy[k] = accuracy[k] + (sum(pred.cv==y[vgr[[i]]]))/length(vgr[[i]])/fold
          
        } ## end of k
        
      }## end of fold
      accuracys[la,c] = max(accuracy)
      indexk = c(indexk, which.max(accuracy))
      
      indexsr = rbind(indexsr, c(la, c, sr))
      
      sr = sr + 1
      
    } ## end of c
  } ## end of lambda
  
  indexs = which(accuracys == max(accuracys), arr.ind = TRUE)
  
  ## optimal tuning parameters
  opt.la = lambdas[indexs[1,1]]; opt.c = Cs[indexs[1,2]]; opt.k = Ks[indexk[which(indexsr[,1]==indexs[1,1] & indexsr[,2]==indexs[1,2])]]
  
  ret@optla <- opt.la
  ret@optc <- opt.c
  ret@optk <- opt.k
  
  
  ## refit on data
  pca.fit <- fpca.face(Y=x, var = TRUE, simul = F, npc = npc, lambda = opt.la, knots = knots)
  
  sig2est = pca.fit$sigma2
  score = pca.fit$scores
  eigenfest = pca.fit$efunctions
  evalest = pca.fit$evalues
  muest = pca.fit$mu
  
  ret@sig2est <- sig2est
  ret@score <- score
  ret@eigenfest <- eigenfest
  ret@evalest <- evalest
  ret@muest <- muest
  ret@npc <- npc
  
  ## calculate FPC scores on  data
  estscore =  estPCscore(x, npc = npc, mu = muest, sigma2 = sig2est, evalues = evalest, eigenfest)
  
  fSVM.fit<-SVM(score, y, opt.k, fold = 0, C = opt.c, kpar = "automatic", kernel = kernel)
  
  svm.fit <- fSVM.fit$svmfit
  
  ret@svm.fit <- svm.fit
  
  ### Perform SVM on data 
  if (fit) {
    predclass = predict(ret, x) 
    accuracy = sum(predclass == y)/m
    
    ret@predclass = predclass
    ret@accuracy = accuracy
    } else NULL
  
  return(ret)
  
})

  

setMethod("predict", signature(object = "FSVC"),
      function (object, newdata, type = "response")
          {
            type <- match.arg(type,c("response","probabilities"))
            
            newdata = as.matrix(newdata)
            newnrows <- nrow(newdata)
            newncols <- ncol(newdata)
            
            estscore =  estPCscore(newdata, npc = object@npc, mu = object@muest, sigma2 = object@sig2est, evalues = object@evalest, efunctions = object@eigenfest)
            
            svm.fit = object@svm.fit
            
            #### Perform SVM on data 
            opt.k = object@optk
            predclass <- predict(svm.fit, as.matrix(estscore[,1:opt.k]))
            
            return(predclass)
          }
)
            

