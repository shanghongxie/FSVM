### Simulation studies for binary outcomes

source("FSVC.R")

smoothers = c(0.5,1,5,10)
Cs = seq(0.01, 1, length.out = 5)
npc = 5
Ks = 1:5
accuracys = NULL
opt.Ks = NULL
opt.las = NULL
opt.Cs = NULL

nSim = 100
for (iSim in 1:nSim){
  
  cat(iSim)
  output =  dataSim[[iSim]]
  x = t(output$discrete_data) ## N*Ntime
  y = output$trainclass
  
  fit = FSVC(x, y, kernel = "rbfdot", Ks, smoothers, Cs = Cs,  npc = 5, knots = 35, fold = 5, fit = TRUE)  
  opt.k = fit@optk
  opt.la = fit@optla
  opt.c = fit@optc
  
  opt.Ks = c(opt.Ks, opt.k)
  opt.las = c(opt.las, opt.la)
  opt.Cs = c(opt.Cs, opt.c)
  
  predtest = predict(fit, t(test_discrete))
  
  Ntest = ncol(test_discrete)
  accuracy = sum(predtest == testclass)/Ntest
  accuracys = c(accuracys, accuracy)
  
}

