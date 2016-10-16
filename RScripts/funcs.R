library(igraph)
library(NetIndices)
library(rnetcarto)
library(plyr)
library(ggplot2)
library(reshape2)
library(parallel)
library(doSNOW)
library(rend)
library(rootSolve)

ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-10] <- 0 
    return(c(states))
  })
}

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

motif_counter <- function(graph.lists){
  require(igraph)
  
  if(!is.list(graph.lists)){
    stop("The input should be a list of graph objects")
  }
  
  triad.count <- lapply(graph.lists, triad.census)
  triad.matrix <- matrix(unlist(triad.count), nrow = length(graph.lists), ncol = 16, byrow = T)
  colnames(triad.matrix) <- c("empty", "single", "mutual", "s5", "s4", "s1", "d4",
                              "d3", "s2", "s3","d8", "d2", "d1", "d5", "d7", "d6")
  
  triad.df <- as.data.frame(triad.matrix)
  
  motif.data.frame <- data.frame(s1 = triad.df$s1, s2 = triad.df$s2, s3 = triad.df$s3, s4 = triad.df$s4, 
                                 s5 = triad.df$s5, d1 = triad.df$d1, d2 = triad.df$d2, d3 = triad.df$d3, d4 = triad.df$d4,
                                 d5 = triad.df$d5, d6 = triad.df$d6, d7 = triad.df$d7, d8 = triad.df$d8)
  
  return(motif.data.frame)
}

getmot <- function(web){
  com <- combn(length(V(web)), 3)
  adj <- get.adjacency(web, sparse = F)
  
  tbt <- lapply(1:ncol(com), function(x){adj[com[,x], com[,x]]})
  conn <- sapply(tbt, function(x) is.connected(graph.adjacency(x)))
  mots <- motif_counter(lapply(tbt[conn], graph.adjacency))
  
  m1 <- melt(lapply(1:sum(conn), function(x) com[,x]))
  m2 <- mots[m1$L1,]
  m3 <- aggregate(m2, list(m1$value), sum)
  
  
  return(m3)
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