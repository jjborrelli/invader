# source()

t0 <- Sys.time()

n1 <- lapply(1:500, function(x) n.vals(60, C = .1))
n2m <- lapply(n1, n.mat)       
xiL <- lapply(n2m, function(x){(((10^2)^TrophInd(x)$TL)/100)^-.25*.314})
riL <- lapply(n2m, function(x){(colSums(x) == 0)})
pars <- lapply(1:500, function(x){list(K = 1, x.i = xiL[[x]], yij = 8, eij = .85, xpar = .2, 
                                        B.o = 0.5, r.i = riL[[x]], A = n2m[[x]], G.i = Gi,
                                        FR = Fbd)})

state1 <- lapply(1:500, function(x) runif(60, .5, 1))


#filepath.data <- "D:/jjborrelli/invadr/dynDATA/"
#filepath.data <- "D:/jjborrelli/invadr/dynDATA2/"
#filepath.data <- "D:/jjborrelli/invadr/dynDATAq1/"
filepath.data <- "D:/jjborrelli/invadr/dynDATAcp2/"


cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("n2m", "pars", "state1", "filepath.data"))
registerDoSNOW(cl)
system.time(
  randos <- foreach(i = 1:500, .packages = c("rend", "deSolve")) %dopar% {
    
    dyn <- ode(y = state1[[i]], times = 1:4000, func = CRmod, parms = pars[[i]], 
               events = list(func = ext1, time = 1:4000))
    
    write.csv(dyn, file = paste0(filepath.data, "web-", i, ".csv", collapse = ""))
    return(dyn)
  }
)
stopCluster(cl)

t1 <- Sys.time()
t1-t0

idyn[[1]][1:2, 1:62]


filepath.data <- "D:/jjborrelli/invadr/dynDATA-inv/"

saveRDS(idyn, file = "D:/jjborrelli/invadr/dynDATA-inv/idyn1.Rds")

