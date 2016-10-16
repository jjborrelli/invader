library(parallel)
library(doSNOW)

#filepath.data <- "D:/jjborrelli/invadr/dynDATA/"
#filepath.data <- "D:/jjborrelli/invadr/dynDATA2/"
#filepath.data <- "D:/jjborrelli/invadr/dynDATAq1-inv/"
filepath.data <- "D:/jjborrelli/invadr/dynDATAcp2-inv/"


t1 <- Sys.time()

connecteds <- sapply(1:500, function(x){is.connected(graph.adjacency(n2m[[x]][randos[[x]][4000,-1] > 0,randos[[x]][4000,-1] > 0]))})
  
  
randos2 <- randos[connecteds]

allinvdrs <- list()
for(i in 1:100){  
  
  ##### get invaders
  invdr <- matrix(nrow = 100, ncol = 3)
  for(z in 1:100){
    cond <- FALSE
    while(!cond){
      invdr[z,] <- n.vals(1, .1)
      imat <- n.mat(rbind(n1[connecteds][[i]][which(randos2[[i]][4000,-1] !=0),], invdr[z,]))
      
      cond <- is.connected(graph.adjacency(imat))
    } 
  }
  
  
  allinvdrs[[i]] <- invdr 
  print(i)
}
Sys.time() - t1  
saveRDS(object = allinvdrs, file = paste0(filepath.data, "allinvdrscp2.rds"))  
  
t2 <- Sys.time()
for(i in 1:100){
  invmat <- lapply(1:100, function(x) n.mat(rbind(n1[connecteds][[i]], allinvdrs[[i]][x,])))

  ###### get parameters
  tind <- lapply(invmat, TrophInd)
  state2 <- c(randos2[[i]][4000,-1], .1)
  allpar <- lapply(1:100, function(x){
    xi <- (((10^2)^tind[[x]]$TL)/100)^-.25*0.314
    ri <- c(riL[connecteds][[i]],(sum(invmat[[x]][,ncol(invmat[[x]])]) == 0))
    par <- list(K = 1, x.i = xi, yij = 8, eij = .85, xpar = 1, 
                B.o = 0.5, r.i = ri, A = invmat[[x]], G.i = Gi, FR = Fij)
    return(par)
  })
  
  
  cl <- makeCluster(detectCores()-1)
  clusterExport(cl, c("allpar", "state2", "ext1"))
  registerDoSNOW(cl)
  
  invdyn <- foreach(x = 1:100,.packages = c("deSolve", "rend")) %dopar% {
    
    
    idyn <- ode(y = state2, times = 1:2000, func = CRmod, parms = allpar[[x]], 
                     events = list(func = ext1, time = 1:2000))
    return(idyn)
  }
  
  stopCluster(cl)
  
  saveRDS(invdyn, file = paste0(filepath.data, "inv-", i, ".rds", collapse = ""))
  
  cat(" ", i, ":", Sys.time() - t2, "---")
}

Sys.time()-t2


test <- readRDS(paste0(filepath.data, "inv-2.rds", collapse = ""))
sum(sapply(test, function(x) x[2000, 62]) > 0)

