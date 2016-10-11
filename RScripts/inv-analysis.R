filepath.inv <- "D:/jjborrelli/invadr/dynDATA-inv/"
filepath.data <- "D:/jjborrelli/invadr/dynDATA/"

dyna.init <- readRDS(file = paste0(filepath.data, "dyn-webs-conn.rds"))
niches <- readRDS(file = paste0(filepath.data, "webs-all.rds"))
conn.webs <- readRDS(file = paste0(filepath.data, "connecteds.rds"))
invdrs.all <- readRDS(file = paste0(filepath.data, "allinvdrs.rds"))


n2m <- n.mat(niches)

n2m2 <- n2m[connecteds][-1]

wp.init <- t(sapply(n2m2, webprops))
wp.init

grep(pattern = "inv",list.files(filepath.inv))
