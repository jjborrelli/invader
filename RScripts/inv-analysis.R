source("./invader/RScripts/funcs.R")

filepath.inv <- "D:/jjborrelli/invadr/dynDATA-inv/"
filepath.data <- "D:/jjborrelli/invadr/dynDATA/"

dyna.init <- readRDS(file = paste0(filepath.data, "dyn-webs-conn.rds"))
niches <- readRDS(file = paste0(filepath.data, "webs-all.rds"))
conn.webs <- readRDS(file = paste0(filepath.data, "connecteds.rds"))
invdrs.all <- readRDS(file = paste0(filepath.data, "allinvdrs.rds"))


n2m <- lapply(niches, n.mat)

n2m2 <- n2m[conn.webs][-1]

wp.init <- t(sapply(n2m2, webprops))
wp.init
hist(sqrt(wp.init[,"C"]*wp.init[,"S"]))

wp.final <- t(sapply(lapply(1:330, function(x) n2m2[[x]][tail(dyna.init[[x]], 1)[-1] > 0,tail(dyna.init[[x]], 1)[-1] > 0]), webprops))


mod.init <- lapply(n2m2, n.mod)
mod.final <- lapply(lapply(1:330, function(x) n2m2[[x]][tail(dyna.init[[x]], 1)[-1] > 0,tail(dyna.init[[x]], 1)[-1] > 0]), n.mod)

ifile <- grep(pattern = "inv", list.files(filepath.inv))

dim(wp.init)

num.invs <- c()
for(i in 1:330){
  idy <- readRDS(file = paste0(filepath.inv, list.files(filepath.inv)[ifile[i]]))
  num.invs[i] <- sum(sapply(idy, function(x) tail(x[2000,-1], 1) > 0))
}

d.cv <- sapply(dyna.init, function(x) sd(x[3990:4000, -1])/mean(x[3990:4000, -1]))[1:330]
hist(d.cv)

complx <- sqrt(wp.final[,c("C")]*wp.final[,c("S")])

d.cv <- sapply(dyna.init, function(x) (apply(x[3995:4000, -1], 2, function(y) sd(y)/mean(y))[which(x[4000,-1] > 0)]))[1:330]

summary(glm(cbind(num.invs, 100-num.invs)~d.cv, family = "binomial"))
summary(glm(cbind(num.invs, 100-num.invs)~wp.final[,"C"]+d.cv, family = "binomial"))
summary(glm(cbind(num.invs, 100-num.invs)~wp.final[,c("S","C")], family = "binomial"))
summary(glm(cbind(num.invs, 100-num.invs)~complx, family = "binomial"))

nm2f <- lapply(1:330, function(x) n2m2[[x]][tail(dyna.init[[x]], 1)[-1] > 0,tail(dyna.init[[x]], 1)[-1] > 0])

web.init <- lapply(n2m2, graph.adjacency)
web.final <- lapply(lapply(1:330, function(x) n2m2[[x]][tail(dyna.init[[x]], 1)[-1] > 0,tail(dyna.init[[x]], 1)[-1] > 0]), graph.adjacency)

mot.init <- lapply(web.init, getmot)
mot.final <- lapply(web.final, getmot)

init.mod <- lapply(n2m2, n.mod)
final.mod <- lapply(nm2f, n.mod)

np.init <- lapply(1:330, function(x) nodeprops(n2m2[[x]], init.mod[[x]]))
np.final <- lapply(1:330, function(x) nodeprops(nm2f[[x]], final.mod[[x]]))

motsinit <- do.call(rbind, mot.final)
persists <- unlist(lapply(dyna.init[1:330], function(x) x[4000, -1] > 0))

np.test <- lapply(np.init, function(x) cbind(x$role, x$TL, x$outdeg, x$indeg))
np.test2 <- do.call(rbind, np.test)
fit2 <- glm(persists~np.test2[,3:4], family = "binomial")
summary(fit2)

fit1 <- glm(persists~motsinit$s1+motsinit$s4+motsinit$s5+ rowSums(motsinit[,7:14]), family = "binomial")
summary(fit1)

t1 <- Sys.time()

i = 1

idyn <- readRDS(file = paste0(filepath.inv, list.files(filepath.inv)[ifile[i]]))
inv.init <- lapply(1:100, function(x){
  n.mat(rbind(niches[conn.webs][[i]][tail(dyna.init[[i]], 1)[-1] > 0,], invdrs.all[[i]][x,]))
  })
inv.final <- lapply(1:100, function(x){
  n.mat(rbind(niches[conn.webs][[i]], invdrs.all[[i]][x,])[tail(idyn[[x]], 1)[-1] > 0,])
  })
  
mod.init <- lapply(inv.init, n.mod)
mod.final <- lapply(inv.final, n.mod)
  
wp.inv.in <- cbind(t(sapply(inv.init, webprops)), Mod = sapply(mod.init, "[[", 2)) 
wp.inv.fi <- cbind(t(sapply(inv.final, webprops)), Mod = sapply(mod.final, "[[", 2))
  
np.inv.in <- lapply(1:100, function(x) nodeprops(inv.init[[x]], mod.init[[x]]))
np.inv.fi <- lapply(1:100, function(x) nodeprops(inv.final[[x]], mod.final[[x]]))

isuccess <- as.numeric(sapply(idyn, function(x) tail(x[2000, -1],1) > 0))


Sys.time() - t1