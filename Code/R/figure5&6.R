rm(list = ls())
# setwd("E:/GitHub/LLG/Code/R")
# setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

# dataName = "CPAC200"
dataName = "desikan"
# dataName = "JHU"
# dataName = "slab907"
# dataName = "slab1068"
# dataName = "Talairach"

set.seed(12345)

m = 5
isSVD = 0

require(Matrix)
library(latex2exp)
source("function_collection.R")
tmpList = read_data(dataName, DA=F)
A_all = tmpList[[1]]
n = tmpList[[2]]
M = tmpList[[3]]
rm(tmpList)


library(lattice)
new.palette=colorRampPalette(c("white","black"),space="rgb")

myAt <- seq(0, 1, length.out=20)
myCkey <- list(at=myAt)


add <- function(x) Reduce("+", x)
P = add(A_all)/M

pdf("../../Draft/P_desikan.pdf", family="Times", width=4, height=3.5)
# pdf("../../Draft/P_desikan.pdf", family="CM Roman", width=4, height=3.5)
image(Matrix(P),main=list(label=TeX('$P$ for Desikan')),sub="",
      xlab=list(cex=0),ylab=list(cex=0),
      scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
      at=myAt, lwd=0)
dev.off()

sampleVec = sample.int(M, m)
A_bar = add(A_all[sampleVec])/m
pdf("../../Draft/Abar_desikan_m5.pdf", family="Times", width=4, height=3.5)
# pdf("../../Draft/Abar_desikan_m5.pdf", family="CM Roman", width=4, height=3.5)
image(Matrix(A_bar),main=list(label=TeX('$\\bar{A}$ for Desikan with M=5')),
      sub="",xlab=list(cex=0),ylab=list(cex=0),
      scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
      at=myAt, lwd=0)
dev.off()

####### Estimate dimensions ######
source("getElbows.R")
nElb = 3
dMax = ceiling(n*3/5)
evalVec = ase(A_bar, dMax, isSVD)[[1]]
dHat = getElbows(evalVec, n=nElb, plot=F)[[nElb]]

####### Calculate Phat ######
A.ase = ase(diag_aug(A_bar), dHat, isSVD)
if (dHat == 1) {
  Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
} else {
  Ahat <- A.ase[[3]][,1:dHat] %*% diag(A.ase[[1]][1:dHat]) %*% t(A.ase[[2]][,1:dHat])
}
P_hat = regularize(Ahat)

# image(Matrix(P_hat),main=list(label=TeX('$\\hat{P}$ for Desikan with M=5')),
#             sub="",xlab=list(cex=0),ylab=list(cex=0),
#             scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),lwd=0),
#       split=c(x=3,y=1,nx=3,ny=1))
pdf("../../Draft/Phat_desikan_m5.pdf", family="Times", width=4.5, height=3.5)
# pdf("../../Draft/Phat_desikan_m5.pdf", family="CM Roman", width=4.5, height=3.5)
levelplot(P_hat[1:n,n:1],col.regions=new.palette(20),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX('$\\hat{P}$ for Desikan with M=5')),
          at=myAt, colorkey=myCkey, lwd=0)
dev.off()
print(dHat)



####### Plot the difference ######
Diff_A_bar = abs(A_bar - P)
Diff_P_hat = abs(P_hat - P)

# valLow = 0.35
valLow = 0.4

# nv = (Diff_A_bar<valLow)
# nv[upper.tri(nv,diag=T)] = FALSE
# Diff_A_bar[nv] = 0
# new.palette=colorRampPalette(c("white","black"),space="rgb") 
# print(levelplot(Diff_A_bar[1:n,n:1],col.regions=new.palette(20),xlab=list(cex=0),
#                 ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
#                 main=list(label=TeX('$|\\bar{A} - P|$ for desikan with M=5'),cex=2),
#                 colorkey=FALSE),split=c(x=2,y=1,nx=3,ny=1),more=TRUE)
# 
# nv = (Diff_P_hat<valLow)
# nv[upper.tri(nv,diag=T)] = FALSE
# Diff_P_hat[nv] = 0
# new.palette=colorRampPalette(c("white","black"),space="rgb") 
# print(levelplot(Diff_P_hat[1:n,n:1],col.regions=new.palette(20),xlab=list(cex=0),
#                 ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
#                 main=list(label=TeX('$|\\hat{P} - P|$ for desikan with M=5'),cex=2),
#                 colorkey=FALSE),split=c(x=3,y=1,nx=3,ny=1))


nv = (Diff_A_bar<valLow)
nv[upper.tri(nv,diag=T)] = FALSE
Diff_A_bar[nv] = 0
pdf("../../Draft/Diff2_desikan_m5.pdf", family="Times", width=4, height=3.5)
# pdf("../../Draft/Diff2_desikan_m5.pdf", family="CM Roman", width=4, height=3.5)
levelplot(Diff_A_bar[1:n,n:1],col.regions=new.palette(20),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX('$|\\bar{A} - P|$ for Desikan with M=5')),
          colorkey=FALSE)
dev.off()

nv = (Diff_P_hat<valLow)
nv[upper.tri(nv,diag=T)] = FALSE
Diff_P_hat[nv] = 0
<<<<<<< HEAD
pdf("../../Draft/Diff3_desikan_m5.pdf", family="Times", width=4.5, height=3.5)
=======
  pdf("../../Draft/Diff3_desikan_m5.pdf", family="CM Roman", width=4.5, height=3.5)
>>>>>>> TangRunze/master
levelplot(Diff_P_hat[1:n,n:1],col.regions=new.palette(20),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX('$|\\hat{P} - P|$ for Desikan with M=5')),
          colorkey=list(labels=list()))
dev.off()

###### Write to files, for MATLAB use ######

Diff_between = abs(A_bar - P_hat)
Diff_A_bar = abs(A_bar - P)
Diff_P_hat = abs(P_hat - P)

# valLow = sort(Diff_P_hat, decreasing=T)[0.01*n^2]
# nv = (Diff_P_hat<valLow)
# s = ""
# for (i in 1:(n-1)) {
#   for (j in (i+1):n) {
#     if (nv[i,j]==F) {
#       s = paste0(s,",",i,",",j,",",P_hat[i,j]-P[i,j])
#     }
#   }
# }
# s = substr(s,2,nchar(s))
# write(s,file="../../Result/Edge_Diff_Phat_desikan.csv")

valLow = sort(Diff_A_bar-Diff_P_hat, decreasing=T)[50]
nv = ((Diff_A_bar-Diff_P_hat)<valLow)
s = ""
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    if (nv[i,j]==F) {
      s = paste0(s,",",i,",",j,",",P_hat[i,j]-P[i,j])
    }
  }
}
s = substr(s,2,nchar(s))
write(s,file="../../Result/Edge_Diff_between_desikan.csv")




# rowSumDiffBetween = rowSums(Diff_A_bar-Diff_P_hat)

# valLow = sort(rowSumDiffBetween, decreasing=T)[5]
# nv = (rowSumDiffBetween<valLow)
# s = ""
# for (i in 1:n) {
#   if (nv[i]==F) {
#     s = paste0(s,",",i)
#   }
# }
# s = substr(s,2,nchar(s))
# write(s,file="../../Result/Vertex_Diff_Between_desikan.csv")

# # rowSumDiffPhat = rowSums(Diff_P_hat)
# # 
# # valLow = sort(rowSumDiffPhat, decreasing=T)[5]
# # nv = (rowSumDiffPhat<valLow)
# # s = ""
# # for (i in 1:n) {
# #   if (nv[i]==F) {
# #     s = paste0(s,",",i)
# #   }
# # }
# # s = substr(s,2,nchar(s))
# # write(s,file="../../Result/Vertex_Diff_Phat_desikan.csv")

# label_tex <- function (labels, multi_line = TRUE) 
# {
#     labels <- label_value(labels, multi_line = multi_line)
#     if (multi_line) {
#         lapply(unname(labels), lapply, function(values) {
#             c(TeX(string = as.character(values)))
#         })
#     }
#     else {
#         lapply(labels, function(values) {
#             values <- paste0("list(", values, ")")
#             lapply(values, function(expr) c(TeX(string = expr)))
#         })
#     }
# }

label_tex <- function (labels, multi_line = TRUE) 
{
  labels <- label_value(labels, multi_line = multi_line)
  if (multi_line) {
    lapply(unname(labels), lapply, function(values) {
      c(TeX(string = as.character(values)))
    })
  }
  else {
    lapply(labels, function(values) {
      values <- paste0("list(", values, ")")
      lapply(values, function(expr) c(TeX(string = expr)))
    })
  }
}


<<<<<<< HEAD
# ut <- upper.tri(P)
# df_ex <- data.frame(P=P[ut],A_bar=A_bar[ut],P_hat=P_hat[ut])
# df_ex$id <- 1:sum(ut)
# df_ex$better <- with(df_ex,abs(P_hat-P)<abs(A_bar-P))
# df_ex$P_bin <- floor(df_ex$P*20)/20

# df_ex <- df_ex %>% mutate(P_hat_se=(P_hat-P)^2,A_bar_se=(A_bar-P)^2,re=A_bar_se/P_hat_se) 
# df_ex$P_hat_se_fitted=loess(P_hat_se~P,df_ex)$fitted
# df_ex$A_bar_se_fitted=loess(A_bar_se~P,df_ex)$fitted
# df_ex$re_fitted=df_ex$A_bar_se_fitted/df_ex$P_hat_se_fitted


# df_ex %>% group_by(P_bin) %>% summarize(conditional_re=mean(A_bar_se)/mean(P_hat_se),n=n()) %>%
#   group_by(P_bin) %>%
#   mutate(x=min(P_bin+0.025,1),xmin=P_bin, xmax=min(P_bin+.05,1)) %>%
#   ggplot(aes(x=x,xmin=xmin,xmax=xmax,y=conditional_re,ymin=conditional_re-.01,ymax=conditional_re+.01,label=n,size=n))+
#   geom_point()+
#   geom_rect()+
#   geom_text(aes(y=conditional_re+.1),size=4)+
#   theme(panel.grid.major = element_line(colour="grey95"),
#         panel.grid.minor = element_blank())+
#   theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
#   # theme(legend.text=element_text(size=20,face="bold"))+
#   theme(legend.position="bottom")+
#   scale_x_continuous(TeX("P_{ij}"),breaks=0:5/5)+
#   scale_y_continuous("conditional relative error")+
#   guides(size=FALSE)

# ggsave("../../Draft/difference_vs_truth_compare.pdf", plot=gg+theme(text=element_text(size=10,family="Times")),
#     width=6.5,height=3.5)

# ggplot(df_ex,aes(x=P,y=re_fitted))+geom_line()


# df_ex <- data.frame(P=P[ut],A_bar=A_bar[ut],P_hat=P_hat[ut])
# df_ex$id <- 1:sum(ut)
# # df_ex$better <- with(df_ex,abs(P_hat-P)<abs(A_bar-P))
# df_ex_melt <- melt(df_ex,id.vars=c("P","id"),variable.name="estimator")
# df_ex_melt$est_tex <- ifelse(df_ex_melt$estimator=="A_bar","\\bar{A}","\\hat{P}")
# gg <- ggplot(df_ex_melt, aes(x=P, linetype=estimator, y=(value-P)^2))+
#   # geom_point(alpha=.1)+
#   geom_smooth(method="gaussian")+
#   theme(panel.grid.major = element_line(colour="grey95"),
#         panel.grid.minor = element_blank())+
#   theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
#   # theme(legend.text=element_text(size=20,face="bold"))+
#   theme(legend.position="bottom")+
#   scale_x_continuous(TeX("P_{ij}"),breaks=0:5/5)+
#   scale_y_continuous("estimator difference from P")

# print(gg)

# ggsave("../../Draft/difference_vs_truth_compare.pdf", plot=gg+theme(text=element_text(size=10,family="Times")),
#     width=6.5,height=3.5)
=======
  ut <- upper.tri(P)
df_ex <- data.frame(P=P[ut],A_bar=A_bar[ut],P_hat=P_hat[ut])
df_ex$id <- 1:sum(ut)
df_ex$better <- with(df_ex,abs(P_hat-P)<abs(A_bar-P))
df_ex$P_bin <- floor(df_ex$P*20)/20

df_ex <- df_ex %>% mutate(P_hat_se=(P_hat-P)^2,A_bar_se=(A_bar-P)^2,re=A_bar_se/P_hat_se) 
df_ex$P_hat_se_fitted=loess(P_hat_se~P,df_ex)$fitted
df_ex$A_bar_se_fitted=loess(A_bar_se~P,df_ex)$fitted
df_ex$re_fitted=df_ex$A_bar_se_fitted/df_ex$P_hat_se_fitted


df_ex %>% group_by(P_bin) %>% summarize(conditional_re=mean(A_bar_se)/mean(P_hat_se),n=n()) %>%
  group_by(P_bin) %>%
  mutate(x=min(P_bin+0.025,1),xmin=P_bin, xmax=min(P_bin+.05,1)) %>%
  ggplot(aes(x=x,xmin=xmin,xmax=xmax,y=conditional_re,ymin=conditional_re-.01,ymax=conditional_re+.01,label=n,size=n))+
  geom_point()+
  geom_rect()+
  geom_text(aes(y=conditional_re+.1),size=4)+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  # theme(legend.text=element_text(size=20,face="bold"))+
  theme(legend.position="bottom")+
  scale_x_continuous(TeX("P_{ij}"),breaks=0:5/5)+
  scale_y_continuous("conditional relative error")+
  guides(size=FALSE)

ggsave("../../Draft/difference_vs_truth_compare.pdf", plot=gg+theme(text=element_text(size=10,family="CM Roman")),
       width=6.5,height=3.5)

ggplot(df_ex,aes(x=P,y=re_fitted))+geom_line()


df_ex <- data.frame(P=P[ut],A_bar=A_bar[ut],P_hat=P_hat[ut])
df_ex$id <- 1:sum(ut)
# df_ex$better <- with(df_ex,abs(P_hat-P)<abs(A_bar-P))
df_ex_melt <- melt(df_ex,id.vars=c("P","id"),variable.name="estimator")
df_ex_melt$est_tex <- ifelse(df_ex_melt$estimator=="A_bar","\\bar{A}","\\hat{P}")
gg <- ggplot(df_ex_melt, aes(x=P, linetype=estimator, y=(value-P)^2))+
  # geom_point(alpha=.1)+
  geom_smooth(method="gaussian")+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  # theme(legend.text=element_text(size=20,face="bold"))+
  theme(legend.position="bottom")+
  scale_x_continuous(TeX("P_{ij}"),breaks=0:5/5)+
  scale_y_continuous("estimator difference from P")

print(gg)

ggsave("../../Draft/difference_vs_truth_compare.pdf", plot=gg+theme(text=element_text(size=10,family="CM Roman")),
       width=6.5,height=3.5)

