rm(list = ls())

setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("E:/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

# dataName = "CPAC200"
dataName = "desikan"
# dataName = "JHU"

source("function_collection.R")
require(ggplot2)
eigenResult <- list()
pp_scree <- list()


tmpList = read_data(dataName, DA=F, newGraph=F)
A_all = tmpList[[1]]
n = tmpList[[2]]
M = tmpList[[3]]
rm(tmpList)
Abar = add(A_all)/M
AbarDiagAug = diag_aug(Abar)
eigenResult = eigen(AbarDiagAug)$values

eigenResult = 1 - cumsum(sort(abs(eigenResult), decreasing = T))/sum(abs(eigenResult))


yMax = max(eigenResult)
yMin = min(eigenResult)
df <- data.frame(eval=eigenResult, k=1:n)
label_y <- with(df, .75*yMax+.25*yMin)

pp_scree <- ggplot(df,aes(x=k,y=eval))+
  geom_line()+
  scale_linetype_manual(name="",values=c("longdash","dotted","dotdash"))+
  xlab("order in algebraic") + ylab("ratio")+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.position="none")+
  ggtitle("Desikan")

