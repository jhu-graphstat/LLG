library(plyr)
library(dplyr)
library(ggplot2)

load("../../Data/data_CPAC200.RData")

nmc <- 100

ng <- length(A_all)
nv <- nrow(A_all[[1]])

mrange <- c(1,5,10, 50, 200)

df_res <- data.frame()

Abar <- Reduce("+", A_all)/ng
eval_Abar <- eigen(Abar)$val

for(m in mrange){
	# A_sample <- A_all[sample(ng, 2*m)]
	for( mc in 1:nmc){
		Abar_m <- Reduce("+", A_all[sample(ng,m)])/m
		eval <- eigen(Abar_m)$val
		df_res <- rbind(df_res,
			data.frame(mc=mc,m=m, d=1:nv, eval=eval,eval_Abar=eval_Abar))
	}
}


df_sum <- df_res %>% group_by(m,d) %>%
	summarize(eval_m=mean(eval),eval_err=eval_m-mean(eval_Abar),
		eval_nm=mean(eval)/mean(eval_Abar), eval_s=sd(eval)) %>%
	mutate(dg=round(d/nv*6)+1)

ggplot(df_sum, 
	aes(x=d, y=eval_err,color=factor(m)))+
	geom_line()+scale_x_log10()

df_sum %>% filter(m>2) %>% 
	ggplot(aes(x=d, y=abs(eval_m/eval_s),color=factor(m)))+
	geom_line()+geom_hline(yintercept=10)


df_sum %>% filter(d!=11) %>%
	ggplot(aes(x=factor(m), fill=abs(eval_nm),y=d))+
	geom_raster()

df_sum %>% ggplot(x=eval_m) + geom_histogram()