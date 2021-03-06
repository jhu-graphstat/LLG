
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
# PHat <- regularize(xHat %*% t(xHat))

# sum(xHat[, 12]^2)
# 
# require(hisemi)
# theta <- 2*pi/8
# # R12 <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol=2)
# # R <- directSum(R12, diag(d-2))
# R <- directSum(cos(theta), diag(d-2), cos(theta))
# R[d, 1] <- sin(theta)
# R[1, d] <- -sin(theta)

# yHat <- xHat %*% R
# min(yHat)
# max(yHat)
# min(yHat[, 1:2])
# max(yHat[, 1:2])

# xHat <- xHat %*% R

# sum(xHat[, 12]^2)

require(ggplot2)

df <- data.frame(value=c(xHat),
                 v=rep(sapply(1:n, function(i) {if (i < 10) {return(paste0("0", i))} else {return(paste0(i))}}), times=d),
                 d=rep(sapply(1:d, function(i) {if (i < 10) {return(paste0("0", i))} else {return(paste0(i))}}), each=n))
gg <- ggplot(df, aes(d, v)) + 
  geom_tile(aes(fill = value), colour = "white") + 
  # scale_fill_gradient(low = "white", high = "steelblue")
  scale_fill_gradient(low = "white", high = "grey10") + 
  xlab("dimension") + ylab("vertex")

ggsave("../../Draft/eigenvector.pdf",
       plot=gg+theme(text=element_text(size=10,family="Times")),
       # plot=gg+theme(text=element_text(size=10,family="CM Roman")),
       width=3.5,height=6.5)


df <- data.frame(x=xHat[,1], y=xHat[,2])
gg <- ggplot(df, aes(x, y)) + 
  geom_point() + 
  scale_fill_gradient(low = "white", high = "grey10") + 
  xlab("1st dimension") + ylab("2nd dimension")

ggsave("../../Draft/eigenvector_scatter.pdf",
       plot=gg+theme(text=element_text(size=10,family="Times")),
       # plot=gg+theme(text=element_text(size=10,family="CM Roman")),
       width=3.5,height=3)




# xHat
# s = ""
# for (j in 1:2) {
#   for (i in 1:n) {
#     s = paste0(s,",",xHat[i,j])
#   }
# }
# s = substr(s,2,nchar(s))
# write(s,file="../../Result/eigenvector.csv")
# 
# 
# 
# xHat
# s = ""
# for (j in 2) {
#   for (i in 1:n) {
#     s = paste0(s,",",xHat[i,j])
#   }
# }
# s = substr(s,2,nchar(s))
# write(s,file="../../Result/eigenvector_dim2.csv")




hc <- hclust(dist(xHat))
plot(hc)


