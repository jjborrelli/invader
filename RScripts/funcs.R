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

conversion <- function(tm){
  for(i in 1:nrow(tm)){
    for(j in 1:ncol(tm)){
      if(tm[i,j] == 1){tm[j,i] <- 1}
    }
  }
  return(tm)
}

n.mod <- function(nm1){
  #nm1 <- n.mat(nv)
  nm2 <- conversion(nm1)
  nc1 <- netcarto(nm2)
  
  return(nc1)
}

nodeprops <- function(nm1, mod1){
  m1 <- mod1[[1]][order(as.numeric(as.character(mod1[[1]]$name)), decreasing = F),]
  tind <- TrophInd(nm1)
  g1 <- graph.adjacency(nm1)
  outd <- degree(g1, mode = "out")
  ind <- degree(g1, mode = "in")
  
  res <- cbind(name = NA, module = NA, connectivity = NA, participation = NA, role = NA, tind, outdeg = outd, indeg = ind)
  res[sort(as.numeric(as.character(mod1[[1]]$name))),1] <- as.numeric(as.character(m1[,1]))
  res[sort(as.numeric(as.character(mod1[[1]]$name))),2] <- m1[,2]
  res[sort(as.numeric(as.character(mod1[[1]]$name))),3] <- m1[,3]
  res[sort(as.numeric(as.character(mod1[[1]]$name))),4] <- m1[,4]
  res[sort(as.numeric(as.character(mod1[[1]]$name))),5] <- m1[,5]
  
  return(res)
}

# compute S, C, L/S
webprops <- function(nm1){
  #nm1 <- n.mat(nv)
  S <- nrow(nm1)
  C <- sum(nm1)/(S*(S-1))
  LperS <- sum(nm1)/S
  g <- graph.adjacency(nm1)
  
  apl <- average.path.length(g)
  troph <- TrophInd(nm1)
  meanTL <- mean(troph$TL)
  
  gen <- mean(colSums(nm1))
  vul <- mean(rowSums(nm1))
  
  clust <- transitivity(graph.adjacency(nm1))
  
  return(c(S = S, C = C, LS = LperS, APL = apl, mTL = meanTL, mGen = gen, mVul = vul, CC = clust))
}