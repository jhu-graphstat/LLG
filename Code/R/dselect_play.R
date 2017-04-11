library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

load("../../Data/data_CPAC200.RData")

approx <- function(a_eig, d, scale=1){
	with(a_eig,vectors[,d] %*% diag(values[d]*scale) %*% t(vectors[,d]))
}

nmc <- 100

ng <- length(A_all)
nv <- nrow(A_all[[1]])

mrange <- c(1,2,5,10, 20, 50, 100,200)

df_res <- data.frame()

Abar <- Reduce("+", A_all)/ng
eAbar <- eigen(Abar)
eval_Abar <- eigen(Abar)$val

df_best <- data.frame()


for(mc in 1:nmc){
cat(mc,"\r")
for(m in mrange){
	A_sample <- A_all[sample(ng, m)]

	Abar_sample <- Reduce("+", A_sample)/m
	eAbar_sample <- eigen(Abar_sample)

	best_eval <- diag(t(eAbar_sample$vectors) %*% Abar %*% eAbar_sample$vectors)

	df_best <- rbind(df_best,
		data.frame(mc=mc, m=m,d=1:nv,eval=best_eval,true_eval=eAbar_sample$values,abar_eval=eval_Abar))

	usvt_eig <- eAbar_sample$values > 0.7 * sqrt(nv/ng)
	est_usvt_eig <- approx(eAbar_sample, usvt_eig)
	usvt_svd <- abs(eAbar_sample$values) > 0.7 * sqrt(nv/ng)
	est_usvt_svd <- approx(eAbar_sample, usvt_svd)
	est_best <- eAbar_sample$vectors %*% diag(best_eval) %*% t(eAbar_sample$vectors)

	df_res <- rbind(df_res,
		data.frame(mc=mc, m=m, 
			which=c("USVT Eig","USVT SVD", "Best","Abar"), 
			mse=c(  mean((est_usvt_eig-Abar)^2),
					mean((est_usvt_svd-Abar)^2),
					mean((est_best-Abar)^2),
					mean((Abar_sample-Abar)^2)),
			d=c(sum(usvt_eig),sum(usvt_svd), nv, nv)
			))
}
}

df_res %>% group_by(m,which) %>% summarize(mse=mean(mse)) %>% 
	group_by(m) %>% mutate(re = mse/mse[2]) %>% 
	filter(which!="USVT Eig") %>%
	ggplot(aes(x=m,y=re,color=which))+geom_line()+
	scale_x_log10()+
	scale_y_log10()+annotation_logticks()

df_best %>% ggplot(aes(x=true_eval,y=eval,color=d))+
	geom_point()+facet_wrap(~m)

df_best %>% group_by(mc,m) %>% mutate(diff=c(true_eval[1:(nv-1)]-true_eval[2:nv],0)) %>%
	filter(d==199) %>%
	ggplot(aes(x=diff,y=eval-true_eval,color=log(d)))+
	geom_point(alpha=.1)+facet_wrap(~m,scale="free")

df_best %>% group_by(m,d) %>% 
	summarize(mr=mean(eval)/mean(true_eval),sr=mad(eval/true_eval),eval=mean(eval)) %>%
	filter(m>5) %>% 
	ggplot(aes(x=d,y=mr,color=factor(m),weight=abs(eval)))+
		geom_point()+ylim(c(-2,2))

df_best %>% 
	filter(m>5) %>% 
	ggplot(aes(x=d,y=eval/true_eval,color=factor(m),group=factor(m),weight=abs(eval)))+
		geom_smooth()

df_best %>% melt(id.vars=.(mc,m,d),variable.name="which") %>%
	group_by(m,d,which) %>% summarize(meanev=mean(value)) %>%
	filter(which!= "abar_eval") %>%
	ggplot(aes(x=d,y=meanev, color=which))+geom_line()+facet_wrap(~m)