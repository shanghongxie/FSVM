### Perform functional support vector regression for continuous outcomes

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


SVM.regression <- function(score, y,  K, fold=5, C=1:10, nu = 0.2, kpar, kernel="rbfdot"){  
  # score: estimted FPC scores
  # y: class label
  # K: number of FPCs
  # fold: number of folds for cross-validation
  # C: regularization parameter C
  # nu: tuning parameter nu which controls the fraction of errors
  # kpar: the list of hyper-parameters (kernel parameters) in ksvm function
  # kernel: the kernel function in ksvm function. rbfdot: Gaussian kernel, vanilladot: linear kernel
  
  library(kernlab)
  
  if (kernel=="vanilladot"){
    svmfit <- ksvm(y~as.matrix(score[,1:K]),cross=fold, C=C, kernel="vanilladot", nu = nu, type = "nu-svr")
  }else{
    svmfit <- ksvm(x = as.matrix(score[,1:K]), y = y, C=C, cross=fold, kpar=kpar,  nu = nu, type = "nu-svr")
  }
  
  return(list("svmfit"=svmfit))
  
}

## Perform functional support vector classification
setClass("FSVR", slots = list(opts = "numeric", optc = "numeric",
                              optk = "numeric", optnu = "numeric",  sig2est = "numeric",
                              score = "matrix", eigenfest = "matrix",
                              evalest = "vector", muest = "vector",  npc = "numeric", svm.fit = "ANY", predy = "ANY", error = "numeric")
)

setGeneric("FSVR", function(x, ...) standardGeneric("FSVR"))
setMethod("FSVR",signature(x = "matrix"),
          FSVR <- function(x = NULL, y = NULL, kernel = "rbfdot",  Ks = NULL, smoothers = NULL, Cs = 1, nus = NULL, npc = 5, knots = 35, fold = 5, fit = TRUE){
            
            # x: functional data matrix, an N*Ntime matrix
            # y: class label
            # kernel: the kernel function.  rbfdot: Gaussian kernel, vanilladot: linear kernel
            # Ks: the grid of number of FPCs
            # smoothers: the grid of smoothing parameter in FPCA
            # Cs: the grid of regularization parameter C in SVR
            # nus: the grid of tuning parameter nu in SVR
            # npc: the maximum number of FPCs
            # knots: number of knots to use or the vectors of knots in fpca.face function; defaults to 35
            # fold: number of folds for cross-validation
            # fit: whether to predict class labels.
            
            library(refund)
            
            ret <- new("FSVR")
            
            ### do cross-validation, tuning parameters: C, nu, smoothing parameter smoother, number of components
            ### split data 
            m = length(y)
            vgr <- split(sample(1:m,m),1:fold)
            
            
            lerror = NULL
            for (s in 1:length(smoothers)){
              
              indexk = NULL
              
              cerror = matrix(0, fold, length(Ks)*length(nus)*length(Cs))
              
              for (i in 1:fold) {
                
                cind <-  unsplit(vgr[-i],factor(rep((1:fold)[-i],unlist(lapply(vgr[-i],length)))))
                
                pca.fit <- fpca.face(Y=x[cind,], var = TRUE, simul = F, npc = npc, lambda = smoothers[s], knots = knots)
                
                sig2est.cv = pca.fit$sigma2
                score.cv = pca.fit$scores
                eigenfest.cv = pca.fit$efunctions
                evalest.cv = pca.fit$evalues
                muest.cv = pca.fit$mu
                
                ## calculate PC scores on cross-validation test in CV
                estscore.cv =  estPCscore(x[vgr[[i]],], npc=npc, mu=muest.cv, sigma2=sig2est.cv, evalues=evalest.cv, eigenfest.cv)
                
                sr = 0
                indexsr = NULL
                for (k in 1:length(Ks)){
                  for (j in 1:length(nus)){
                    for (c in 1:length(Cs)){
                      
                      sr = sr+1
                      
                      fSVM.fit<-SVM.regression(score.cv, y[cind], Ks[k], fold = 0, C = Cs[c], nu = nus[j], kpar="automatic",kernel= kernel)
                      
                      svm.fit <- fSVM.fit$svmfit
                      
                      pred.cv <- predict(svm.fit,as.matrix(estscore.cv[,1:Ks[k]]))
                      
                      cerror[i, sr] = cerror[i, sr] + sum((pred.cv-y[vgr[[i]]])^2)/length(vgr[[i]])
                      
                      indexsr = rbind(indexsr, c(k,j, c, sr))
                    }  ## end of C
                    
                  } ## end of j, nu
                  
                } ## end of k
                
              } ## end of fold
              
              lerror = cbind(lerror, apply(cerror,2,mean))  ## average of CV error for each smoother (each column)
            }  ## end of smoother
            
            
            index.smoother = which(lerror == min(lerror), arr.ind = TRUE)[1,2]
            index = which(lerror == min(lerror), arr.ind = TRUE)[1,1]
            index.k = indexsr[index,1]
            index.nu = indexsr[index,2]
            index.c = indexsr[index,3]
            
            opt.s = smoothers[index.smoother]; opt.k = Ks[index.k]; opt.nu = nus[index.nu]; opt.c = Cs[index.c]
            
            ret@opts <- opt.s
            ret@optc <- opt.c
            ret@optk <- opt.k
            ret@optnu <- opt.nu
            
            
            ## refit on training data
            pca.fit <- fpca.face(Y=x, var = TRUE, simul = F, npc = npc, lambda = opt.s, knots = knots)
            
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
            
            ## calculate PC scores on test data
            estscore =  estPCscore(x, npc=npc, mu=muest, sigma2=sig2est, evalues=evalest, eigenfest)
            
            
            fSVM.fit<-SVM.regression(score, y, opt.k, fold = 0, C = opt.c, nu = opt.nu, kpar="automatic",kernel= kernel)
            
            svm.fit <- fSVM.fit$svmfit
            
            ret@svm.fit <- svm.fit
            
            ### Perform SVM on data 
            if (fit) {
              predy = predict(ret, x) 
              error = sum((predy-y)^2)/m
              
              ret@predy = predy
              ret@error = error
            } else NULL
            
            return(ret)
            
          })



setMethod("predict", signature(object = "FSVR"),
          function (object, newdata)
          {
            
            if (is.vector(newdata)) t(t(newdata)) else as.matrix(newdata)
            
            newnrows <- nrow(newdata)
            newncols <- ncol(newdata)
            
            estscore =  estPCscore(newdata, npc = object@npc, mu = object@muest, sigma2 = object@sig2est, evalues = object@evalest, efunctions = object@eigenfest)
            
            svm.fit = object@svm.fit
            
            #### Perform SVM on data 
            opt.k = object@optk
            predy <- predict(svm.fit, as.matrix(estscore[,1:opt.k]))
            
            return(predy)

          }
)


