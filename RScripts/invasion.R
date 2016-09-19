library(igraph)
library(NetIndices)
library(plyr)
library(ggplot2)
library(reshape2)
library(parallel)
library(doSNOW)
library(rend)
library(rootSolve)


n.vals <- function(S,C){
  niche<-runif(S,0,1)
  r<-rbeta(S,1,((1/(2*C))-1))*niche
  
  ci <- runif(S, r/2, niche)
  
  return(cbind(niche = niche, ci = ci, r = r))
}

n.mat <- function(nv){
  S <- nrow(nv)
  new.mat<-matrix(0,nrow=S,ncol=S)
  
  for(i in 1:S){
    
    for(j in 1:S){
      if(nv[j,1]>(nv[i,2]-(.5*nv[i,3])) && nv[j,1]<(nv[i,2]+.5*nv[i,3])){
        new.mat[j,i]<-1
      }
    }
  }
  
  #new.mat<-new.mat[order(apply(new.mat,2,sum)),order(apply(new.mat,2,sum))]
  return(new.mat)
}

motifs2 <- function(g){
  mo <- motifs(g)
  mo2 <- mo[c(5, 8, 12, 3, 7, 14, 9, 10, 6, 13, 16, 15, 11)]  
  names(mo2) <- c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8")
  
  return(mo2)
}

curve_ball<-function(m){
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:(5*l_hp)){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
    
  }
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}

curving <- function(adjmat, n){
  mot <- motif_counter(list(graph.adjacency(adjmat)))
  newmat <- adjmat
  
  for(i in 1:n){
    newmat <- curve_ball(newmat)
    m <- motif_counter(list(graph.adjacency(newmat)))
    mot <- rbind(mot, m)
  }
  return(mot[-1,])
}


nd <- function(gl) nrow(gl) - nrow(unique(aaply(gl, 1, sort)))
# determines the number of double links

dblcan.curve <- function(mat, iter){
  mot <- motif_counter(list(graph.adjacency(mat)))
  
  el <- get.edgelist(graph.adjacency(mat))
  Ne <- nrow(el)
  dbl <- nd(el)
  can <- sum(diag(mat))
  
  for(i in 1:iter){
    ed = TRUE
    dub = TRUE
    ca = TRUE
    while(ed || dub || ca){
      mat2 <- curve_ball(mat)
      
      el2 <- get.edgelist(graph.adjacency(mat2))
      el3 <- unique(el2)
      Ne2 <- nrow(el3)
      dbl2 <- nd(el3)
      can2 <- sum(diag(mat2))
      
      ed <- Ne != Ne2
      dub <- dbl != dbl2
      ca <- can != can2
      
    }
    mat <- mat2
    
    mot <- rbind(mot, motif_counter(list(graph.adjacency(mat))))
  }
  return(M = mot[-1,])
}


N <- sample(50:250, 100)
conn <- runif(100, .1, .3)






n1 <- lapply(1:1000, function(x) n.vals(60, C = .1))
#randos.t3 <- lapply(1:2000, function(x) matrix(NA, ncol = 60, nrow = 1000))

filepath.sink <- "D:/jjborrelli/invadr/"
# Use the filepath.data variable to set the location for 
# storing the dataframes of motif frequencies for each web
filepath.data <- "D:/jjborrelli/invadr/dynDATA/"

cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("n1", "n.mat", "filepath.sink", "filepath.data"))
registerDoSNOW(cl)
system.time(
randos <- foreach(i = 1:1000) %dopar% {
  library(rend)
  library(NetIndices)
  sink(file = paste0(filepath.sink, "web-", i, ".txt"))
  mat1 <- n.mat(n1[[i]])
  xi <- (((10^2)^TrophInd(mat1)$TL)/100)^-.25*.314
  dyn <- CRsimulator(mat1, t = 1:4000, yij = 8, eij = .85, x.i = xi)
  print(dyn)
  sink()
  write.csv(dyn, file = paste0(filepath.data, "web-", i, ".csv"))
  return(dyn)
}
)
stopCluster(cl)




spec <- sapply(lapply(randos, function(x) x[1000,-1]), function(y) sum(y > 0))
sum(spec >=50)
hist(spec)

nm4k <- lapply(1:1000, function(x) n.mat(n1[[x]][which(randos[[x]][4000,-1] != 0),]))
nval4k <- lapply(1:1000, function(x) n1[[x]][which(randos[[x]][4000,-1] != 0),])

con4k <- sapply(nm4k, function(x) sum(x)/(nrow(x)*(nrow(x) -1)))
hist(sqrt(con4k*spec))

eq.comm <- lapply(randos, function(x) which(x[4000,-1] != 0))
eq.comm[[1]]

eq.abund <- lapply(randos, function(x) x[4000, -1][which(x[4000,-1] != 0)])
eq.abund[[1]]


invdrs <- n.vals(100, .1)
inmat1 <- n.mat(rbind(n1[[1]], invdrs[1,]))
inmat2 <- n.mat(rbind(n1[[1]][eq.comm[[1]],], invdrs[1,]))
which(colSums(inmat1) == 0)

is.connected(graph.adjacency(inmat2))
testcon <- lapply(1:1000, function(x){
  myconn <- c()
  for(i in 1:100){
    inmat2 <- n.mat(rbind(n1[[x]][eq.comm[[x]],], invdrs[i,]))
    myconn[i] <- is.connected(graph.adjacency(inmat2))
  }
  return(myconn)
})


woi <- which(sapply(testcon, sum) == 100)


filepath.sink <- "D:/jjborrelli/invadr/"
# Use the filepath.data variable to set the location for 
# storing the dataframes of motif frequencies for each web
filepath.data <- "D:/jjborrelli/invadr/dynDATA-inv/"

cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("n1", "n.mat", "woi", "eq.comm", "eq.abund", "invdrs", "filepath.sink", "filepath.data"))
registerDoSNOW(cl)
system.time(
  invs <- foreach(i = 99:length(woi)) %dopar% {
    library(rend)
    library(NetIndices)
    sink(file = paste0(filepath.sink, "web-", i, ".txt"))
    
    dyn <- list()
    for(x in 1:100){
      mat1 <- n.mat(rbind(n1[[i]][eq.comm[[i]],], invdrs[x,]))
      mat2 <- n.mat(rbind(n1[[i]], invdrs[x,]))
      xi <- (((10^2)^TrophInd(mat1)$TL)/100)^-.25*.314
      ri <- (colSums(mat2) == 0)[c(eq.comm[[i]], 61)]
      dyn[[x]] <- CRsimulator(mat1, states = c(eq.abund[[i]], .1), r = ri,t = 1:2000, 
                              yij = 8, eij = .85, x.i = xi)
    }

    print(dyn)
    sink()
    save(dyn, file = paste0(filepath.data, "webINV-", i, ".Rdata"))
    #return(dyn)
  }
)
stopCluster(cl)



mat2 <- lapply(1:100, function(x) n.mat(rbind(n1[[1]], invdrs[x,])))
mat1 <- lapply(1:100, function(x) n.mat(rbind(n1[[1]][eq.comm[[1]],], invdrs[x,])))
xi <- lapply(mat1, function(x) (((10^2)^TrophInd(x)$TL)/100)^-.25*0.314)

s.t <- Sys.time()
i=1
dyn <- list()
for(x in 1:100){
  ri <- (colSums(mat2[[x]]) == 0)[c(eq.comm[[1]], 61)]
  dyn[[x]] <- CRsimulator(mat1[[x]], states = c(eq.abund[[1]], .1), r = ri,t = 1:2000, 
                          yij = 8, eij = .85, x.i = xi[[x]])
  print(x)
}
e.t <- Sys.time()
e.t-s.t
#1.95 hrs


s.t2 <- Sys.time()
cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("eq.comm", "eq.abund", "invdrs", "mat1", "mat2", "xi"))
registerDoSNOW(cl)

invs <- foreach(x = 1:100,.packages = c("rend")) %dopar% {

  ri <- (colSums(mat2[[x]]) == 0)[c(eq.comm[[1]], 61)]
  dyn[[x]] <- CRsimulator(mat1[[x]], states = c(eq.abund[[1]], .1), r = ri,t = 1:2000, 
                          yij = 8, eij = .85, x.i = xi[[x]])
  
  return(dyn)
}

stopCluster(cl)

e.t2
e.t2-s.t2


par <- list(K = 1, x.i = xi, yij = 8, eij = .85, xpar = .2, 
            B.o = 0.5, r.i = ri, A = mat1[[x]], G.i = Gi, FR = Fij)
state1 <- c(eq.abund[[1]], .1)

out <- ode(y = state1, times = 1:2000, func = CRmod, parms = par, events = list(func = goExtinct, time = 1:2000))




  
  
s.tB <- Sys.time()
i=1
#dyn <- list()
for(x in 1:100){
  ri <- (colSums(mat2[[x]]) == 0)[c(eq.comm[[1]], 61)]
  par <- list(K = 1, x.i = xi[[x]], yij = 8, eij = .85, xpar = .2, 
              B.o = 0.5, r.i = ri, A = mat1[[x]], G.i = Gi, FR = Fij)
  state1 <- c(eq.abund[[1]], .1)
  
  dyn[[x]] <- ode(y = state1, times = 1:2000, func = CRmod, parms = par,
                  events = list(func = goExtinct, time = 1:2000))
  print(x)
}
e.tB <- Sys.time()
e.tB-s.tB


s.t2B <- Sys.time()
cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("eq.comm", "eq.abund", "invdrs", "mat1", "mat2", "xi"))
registerDoSNOW(cl)

invs <- foreach(x = 1:100,.packages = c("rend")) %dopar% {
  
  ri <- (colSums(mat2[[x]]) == 0)[c(eq.comm[[1]], 61)]
  par <- list(K = 1, x.i = xi[[x]], yij = 8, eij = .85, xpar = .2, 
              B.o = 0.5, r.i = ri, A = mat1[[x]], G.i = Gi, FR = Fij)
  state1 <- c(eq.abund[[1]], .1)
  
  dyn <- ode(y = state1, times = 1:2000, func = CRmod, parms = par,
                  events = list(func = goExtinct, time = 1:2000))
  
  return(dyn)
}

stopCluster(cl)

e.t2B
e.t2B-s.t2B
