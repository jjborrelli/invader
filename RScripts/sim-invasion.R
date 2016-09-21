t0 <- Sys.time()

library(igraph)
library(rnetcarto)
library(NetIndices)
library(deSolve)
library(rend)



conversion <- function(tm){
  for(i in 1:nrow(tm)){
    for(j in 1:ncol(tm)){
      if(tm[i,j] == 1){tm[j,i] <- 1}
    }
  }
  return(tm)
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

# compute modularity and node participation
n.mod <- function(nm1){
  #nm1 <- n.mat(nv)
  nm2 <- conversion(nm1)
  nc1 <- netcarto(nm2)
  
  return(nc1)
}

nodeprops <- function(nm1, mod1){
  m1 <- mod1[order(as.numeric(rownames(mod1)), decreasing = F),]
  tind <- TrophInd(nm1)
  g1 <- graph.adjacency(nm1)
  outd <- degree(g1, mode = "out")
  ind <- degree(g1, mode = "in")
  
  return(cbind(m1, tind, outdeg = outd, indeg = ind))
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

ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-10] <- 0 
    return(c(states))
  })
}

t1 <- Sys.time()

test1 <- n.vals(60, .15)
testmat <- n.mat(test1)
mod1 <- n.mod(testmat)
mod2 <- mod1[[1]]
wp1 <- c(webprops(testmat), Mod = mod1[[2]])

xi <- (((10^2)^TrophInd(testmat)$TL)/100)^-.25*0.314
ri <- (colSums(testmat) == 0)
par <- list(K = 1, x.i = xi, yij = 8, eij = .85, xpar = .2, 
            B.o = 0.5, r.i = ri, A = testmat, G.i = Gi, FR = Fij)
state1 <- runif(60, .5, 1)

dyn <- ode(y = state1, times = 1:4000, func = CRmod, parms = par, 
           events = list(func = ext1, time = 1:4000))



###
### Equilibrium 

eqcom <- which(dyn[4000, -1] > 0)
eqmat <- testmat[eqcom, eqcom]
mod.eq <- n.mod(eqmat)
wp.eq <- c(webprops(eqmat), Mod = mod.eq[[2]])
wp.eq -wp1


###
### Invasion

invdr <- n.vals(100, .15)

invmat <- list()
idyn <- lapply(1:100, function(x) matrix(0, nrow = 2000, ncol = 61))
m <- list()
wp.inv <- matrix(0, nrow = 100, ncol = 9)
colnames(wp.inv) <- c("N", "C", "LS", "APL", "mTL", "mGen", "mVul", "CC", "Mod")

invmat <- lapply(1:100, function(x) n.mat(rbind(test1, invdr[x,])))
tind <- lapply(invmat, TrophInd)
state2 <- c(dyn[4000,-1], .1)
allpar <- lapply(1:100, function(x){
  xi <- (((10^2)^tind[[x]]$TL)/100)^-.25*0.314
  ri <- (colSums(invmat[[x]]) == 0)
  par <- list(K = 1, x.i = xi, yij = 8, eij = .85, xpar = .2, 
              B.o = 0.5, r.i = ri, A = invmat[[x]], G.i = Gi, FR = Fij)
  return(par)
})


#for(i in 1:100){
  #is.connected(graph.adjacency(invmat))
  
  #imat2 <- n.mat(rbind(test1[eqcom,], invdr[i,]))
  #is.connected(graph.adjacency(invmat))
  
  #modinv <- n.mod(imat2)
  #m[[i]] <- modinv[[1]][order(as.numeric(rownames(modinv[[1]])), decreasing = F),]
  #m[[i]]$names <- c(eqcom, 61)

  #m[[i]] <- cbind(m, tind[c(eqcom, 61),])
  
  #wp.inv[i,] <- c(webprops(imat2), Mod = modinv[[2]])
  
  ####
  #### Post-invasion dynamics
for(i in 1:100){  
  idyn[[i]] <- ode(y = state2, times = 1:2000, func = CRmod, parms = allpar[[i]], 
              events = list(func = ext1, time = 1:2000))
  print(i)
}

matplot(sapply(idyn, function(x) x[,62]), typ = "l", lwd = 2)
sum(sapply(idyn, function(x) x[2000,62]) == 0)


t2 <- Sys.time()

t2 - t1
t2 - t0
