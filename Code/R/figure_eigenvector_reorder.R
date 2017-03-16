
setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("E:/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

isSVD = 0
m = 5
dataName = "desikan"
set.seed(12345)


source("function_collection.R")
tmpList = read_data(dataName, DA=F, newGraph=F)
A_all = tmpList[[1]]
n = tmpList[[2]]
M = tmpList[[3]]
rm(tmpList)

# dVec = 1:n
# nD = length(dVec)

A_sum = add(A_all)

source("getElbows.R")
source("USVT.R")

sampleVec = sample.int(M, m)
A_bar = add(A_all[sampleVec])/m

A_bar_diag_aug = diag_aug(A_bar)

# ZG
nElbow = 3
evalVec = ase(A_bar, ceiling(n*3/5), isSVD)[[1]]
dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]

A.ase = ase(A_bar_diag_aug, dZG, isSVD)
d = dZG

xHat <- A.ase[[3]] %*% diag(sqrt(A.ase[[1]]))


hc <- hclust(dist(xHat))
# plot(hc)
nv <- order(hc$order)

require(ggplot2)
# xHat <- xHat[nv,]
df <- data.frame(value=c(xHat),
                 order=rep(sapply(1:n, function(i) {nv[i]}), times=d),
                 v=rep(sapply(1:n, function(i) {if (i < 10) {return(paste0("0", i))} else {return(paste0(i))}}), times=d),
                 d=rep(sapply(1:d, function(i) {if (i < 10) {return(paste0("0", i))} else {return(paste0(i))}}), each=n))
gg <- ggplot(df, aes(x=d, y=reorder(v,order))) + 
  geom_tile(aes(fill = value), colour = "white") + 
  # scale_fill_gradient(low = "white", high = "steelblue")
  scale_fill_gradient(low = "white", high = "grey10") + 
  xlab("dimension") + ylab("vertex")

ggsave("../../Draft/eigenvector_reorder.pdf",
       plot=gg+theme(text=element_text(size=10,family="Times")),
       # plot=gg+theme(text=element_text(size=10,family="CM Roman")),
       width=3.5,height=6.5)