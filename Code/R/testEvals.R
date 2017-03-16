rm(list = ls())

set.seed(12345)

setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("E:/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

# dataName = "CPAC200"
# dataName = "desikan"
# dataName = "JHU"

dataNameVec = c("JHU", "desikan", "CPAC200")

source("function_collection.R")
require(ggplot2)
library(dplyr)
pp = list()

mVec = c(1, 5, 10, 454)

for (iData in 1:length(dataNameVec)) {
  
  dataName = dataNameVec[iData]
  
  for (iM in 1:length(mVec)) {
    m = mVec[iM]
    
    tmpList = read_data(dataName, DA=F, newGraph=F)
    A_all = tmpList[[1]]
    n = tmpList[[2]]
    M = tmpList[[3]]
    rm(tmpList)
    
    nv = sample(length(A_all), m)
    A_all = A_all[nv]
    
    if (iM == 1) {
      eigenResult = array(0, dim=c(length(mVec), n))
    }
    
    Abar = add(A_all)/m
    AbarDiagAug = diag_aug(Abar)
    eigenResult[iM, ] = eigen(AbarDiagAug)$values
  }

  error_by_dim_df <- data.frame(eval=c(eigenResult), k=rep(1:n, each=length(mVec)),
                                m=rep(mVec,n)) %>%
    mutate(m=factor(paste0("M=",m),sapply(mVec, function(m) {paste0("M=", m)})))
  
  label_y <- with(error_by_dim_df, .75*max(eval)+.25*min(eval))
  
  gg <- ggplot(error_by_dim_df,aes(x=k,y=eval))+
    facet_wrap(~m, nrow=1)+
    scale_linetype_manual(name="",values=c(1))+
    scale_shape_manual(name="",values=c(-1))+
    #   geom_point(dim_selection_df,aes(shape=which))+
    geom_line()+
    # geom_linerange(aes(ymin=lci,ymax=uci),linetype=1,alpha=.5,size=.5)+
    #   geom_vline(data=dim_selection_df,
    #              aes(xintercept=value,color=which,linetype=variable))+
    #   scale_linetype_manual(name="",values=c(1,2,3,4))+
    #   geom_text(data=dim_selection_df %>% filter(variable=="mean"),
    #             aes(x=value+n/30,y=label_y,linetype=variable,label=which,color=which),angle=90)+
    #   scale_color_discrete(guide=FALSE)+
    xlab("")+ylab("Eigenvalue")+
    # theme(strip.text.x = element_text(size=20,face="bold"))+
    # theme(axis.text=element_text(size=15),
    #       axis.title=element_text(size=20,face="bold"))+
    theme(panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor = element_blank())+
    theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
    theme(legend.position="none")+
    ggtitle(paste0(dataName, ", N=", n))
  
  pp[[iData]]=gg
  
  # ggsave(paste0("../../Draft/screeplot_", dataName, "_m_", m, ".pdf"),
  #        # p+theme(text=element_text(size=10,family="Times")),
  #        pp[[iData]]+theme(text=element_text(size=10,family="CM Roman")),
  #        width=2, height=2)
}

multiplot(plotlist = pp, cols = 1)
