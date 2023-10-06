### Simulation studies for continuous outcomes


source("FSVR.R")


smoothers=c(0.5,1,5,10)
npc = 5
Ks = 1:5
nus = c(0.1, 0.3, 0.5, 0.8)
Cs = seq(0.01, 1, length.out = 5)

rmses = NULL
opt.Ks = NULL
opt.ss = NULL
opt.Cs = NULL
opt.nus = NULL

nSim = 100
for (iSim in 1:nSim){
  
  cat(iSim)
  output =  dataSim[[iSim]]
  x = t(output$discrete_data) ## N*Ntime
  y = output$trainy
  
  fit = FSVR(x, y, kernel = "rbfdot", Ks, smoothers, Cs = Cs, nus = nus, npc = 5, knots = 35, fold = 5, fit = TRUE)   
  opt.k = fit@optk
  opt.s = fit@opts
  opt.c = fit@optc
  opt.nu = fit@optnu
  
  opt.Ks = c(opt.Ks, opt.k)
  opt.ss = c(opt.ss, opt.s)
  opt.Cs = c(opt.Cs, opt.c)
  opt.nus = c(opt.nus, opt.nu)
  
  predtest = predict(fit, t(test_discrete))
  
  Ntest = ncol(test_discrete)
  error = sqrt(sum((predtest-testy)^2)/Ntest)
  rmses = c(rmses, error)
  
}

