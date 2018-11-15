#|**********************************************************************;
#  This script reproduces the results in 
#  * Project: EXACT HYPOTHESIS TESTING FOR SHRINKAGE BASED GGM (Bernal et al).
#  * Author            : Victor Bernal*, Rainer Bischoff, Victor Guryev, Marco Grzegorczyk, Peter Horvatovich
#  * Date created      : 2018-11-12
#***********************************
# Description
# The aim is to show that the reconstruction of (shrinkage based) GGMs is biased
# unless the shrinkage parameter is included in the inference. 
# To illustrate thi we will highlight that with the new shrunk density 
# (i) p values are correctly distributed i.e. uniform [0 1] under the null hypothesis: no partial correlation
# (ii) expected false positives rate are accurate controlled
# (iii) the positive predicted values is superior
# As gold standard we will use the expensive Monte Carlo p value estimations 
#**********************************
# Prerequisite R library (CRAN)
# requires R libraries "GeneNet" , "stats4", "ggplot2", and "Hmisc".
#************************************************************************
# Parameters       
# * p = Number of variables (e.g. genes)   
# * n = Number of samples  
# * number= Montecarlo iterations
# * rep = Times to repeat the simulation (for fixed p,n)
# 
# * etaA = proportion of TP
# * alpha = significance level
# * eta0 = 1-etaA
#************************************************************************
#  * Revision History  : none
#  **********************************************************************
#  * Details
#  The shrunk probability density is presented here in [] 
#  For empirical p value with Monte Carlo see  [] 
#  The standard probability density is presented in [] 
#  The simulation of data, as well as estimation with the optimal shrinkage (lambda) is done with
#  [Schäfer,J. and Strimmer,K. (2005a),  GeneNet 1.2.13. CRAN.]

# 1) Here we simulate different samples sizes (n[times]) and repeat it (rep).
# We will study the (true and false) positives, positive predictive value PPV vs the sample size

rm(list=ls())

library(GeneNet)
library(stats4)
library(Hmisc)
library(ggplot2)

set.seed(123)
source("shrinkagefunctions.R") # Check that shrinkagefunctions.R in the same directory

#***************************
# Initialize parameters
p<-100  # num of genes
n<-c(1:15)*10  # num of samples
number<-5 # Montecarlo iterations
rep<- 5 # repeat each simulation (for fixed p,n)
etaA=0.03  # proportion of  TP
alpha<-0.01 # significance level
eta0<-1-etaA

sim.pcor<-ggm.simulate.pcor(p, etaA) #simulate partial corr
temp<-sm2vec(sim.pcor)
TP<- which(temp!=0) # true positives
etaA*p*(p-1)/2 # expected number of TP

## Initialize variables
all<-matrix(Inf,3,rep) 
All<-matrix(Inf,3,length(n)) 
se.All<-matrix(Inf,3,length(n))
correct<-matrix(Inf,3,rep) 
incorrect<-matrix(Inf,3,rep) 
ppv<-matrix(Inf,3,rep) 
tp<-matrix(Inf,3,length(n)) 
FP<-matrix(Inf,3,length(n)) 
PPV<-matrix(Inf,3,length(n))  
se.TP<-matrix(Inf,3,length(n)) 
se.FP<-matrix(Inf,3,length(n)) 
se.PPV<-matrix(Inf,3,length(n))  
L<-matrix(Inf,2,length(n))
all.lambda<-matrix(Inf,1,rep) 
rownames(L)<-c("lambda","std lambda")
BH.correct<-matrix(Inf,3,rep) 
BH.incorrect<-matrix(Inf,3,rep) 
BH.ppv<-matrix(Inf,3,rep) 
BH.PPV<-matrix(Inf,3,length(n)) 
BH.se.PPV<-matrix(Inf,3,length(n)) 
BF.correct<-matrix(Inf,3,rep) 
BF.incorrect<-matrix(Inf,3,rep) 
BF.ppv<-matrix(Inf,3,rep) 
BF.PPV<-matrix(Inf,3,length(n)) 
BF.se.PPV<-matrix(Inf,3,length(n)) 
#***************************************


for (times in 1:length(n)) { 
  
  for (again in 1:rep){      
    
    sim.data<- ggm.simulate.data( n[times] , sim.pcor)
    # Uncomment if you want to fix lambda by hand
    #GGM <- pcor.shrink(sim.data,lambda ,verbose=TRUE) 
    GGM <- pcor.shrink(sim.data, verbose=TRUE) # OPTIMAL
    lambda<-attr(GGM, "lambda") 
    
    
    if (lambda<0.99){ # exclude the case of complete shrinkage (lambda=1)
      
      # kappa(GGM, exact=TRUE, norm='2' ) # condition number
      r<-sm2vec(GGM)
      
      
      ##### P-VALUES
      
      # GeneNet  P values 
      # current method is Empirical null fitting (ENF)
      # the test must be in node's order. Sort by node 2 and by node 1
      ENF.test<-network.test.edges(GGM, fdr=TRUE, plot=FALSE)
      ENF.test<-ENF.test[order(ENF.test$node1,ENF.test$node2),]
      BH.ENF.test<-p.adjust(ENF.test$pval ,  method="BH")
      BF.ENF.test<-p.adjust(ENF.test$pval ,  method="bonferroni")
      ### Shrunk MLE  P values
      pval.shrunk<-p.shrunk(r,p,n[times],lambda)
      BH.pval.shrunk<-p.adjust(pval.shrunk ,  method="BH")
      BF.pval.shrunk<-p.adjust(pval.shrunk ,  method="bonferroni")
      ## MONTECARLO P values
      p.monte<-p.montecarlo(r,number,p,n[times],lambda)
      BH.p.monte<-p.adjust(p.monte ,  method="BH")
      BF.p.monte<-p.adjust(p.monte ,  method="bonferroni")
      #################################################
      
      # amount of significant
      all[1,again]<-sum(ENF.test$pval<=alpha)
      all[2,again]<-sum(p.monte<=alpha)
      all[3,again]<-sum( pval.shrunk<=alpha)
      
      # true positives
      correct[1,again]<-sum(which(ENF.test$pval<=alpha) %in% TP )
      correct[2,again]<-sum(which(p.monte<=alpha) %in% TP )
      correct[3,again]<-sum(which(pval.shrunk<=alpha) %in% TP )
      
      ## false positives
      incorrect[1,again]<- sum(!which(ENF.test$pval<=alpha) %in% TP )
      incorrect[2,again]<- sum(!which(p.monte<=alpha) %in% TP )
      incorrect[3,again]<- sum(!which(pval.shrunk<=alpha) %in% TP )
      
      ## ppv
      ppv[1,again]<-correct[1,again]/(correct[1,again]+incorrect[1,again])
      ppv[2,again]<-correct[2,again]/(correct[2,again]+incorrect[2,again])
      ppv[3,again]<-correct[3,again]/(correct[3,again]+incorrect[3,again])
      
      ##############################################################################    
      ## ppv BH
      # true positives
      BH.correct[1,again]<-sum(which(BH.ENF.test<=alpha) %in% TP )
      BH.correct[2,again]<-sum(which(BH.p.monte<=alpha) %in% TP )
      BH.correct[3,again]<-sum(which(BH.pval.shrunk<=alpha) %in% TP )
      
      ## false positives
      BH.incorrect[1,again]<-sum(!which(BH.ENF.test<=alpha) %in% TP )
      BH.incorrect[2,again]<- sum(!which(BH.p.monte<=alpha) %in% TP )
      BH.incorrect[3,again]<- sum(!which(BH.pval.shrunk<=alpha) %in% TP )
      
      BH.ppv[1,again]<-BH.correct[1,again]/(BH.correct[1,again]+BH.incorrect[1,again])
      BH.ppv[2,again]<-BH.correct[2,again]/(BH.correct[2,again]+BH.incorrect[2,again])
      BH.ppv[3,again]<-BH.correct[3,again]/(BH.correct[3,again]+BH.incorrect[3,again])
      ##############################################################################    
      ## ppv BF
      # true positives
      BF.correct[1,again]<-sum(which(BF.ENF.test<=alpha) %in% TP )
      BF.correct[2,again]<-sum(which(BF.p.monte<=alpha) %in% TP )
      BF.correct[3,again]<-sum(which(BF.pval.shrunk<=alpha) %in% TP )
      
      ## false positives
      BF.incorrect[1,again]<-sum(!which(BF.ENF.test<=alpha) %in% TP )
      BF.incorrect[2,again]<- sum(!which(BF.p.monte<=alpha) %in% TP )
      BF.incorrect[3,again]<- sum(!which(BF.pval.shrunk<=alpha) %in% TP )
      
      BF.ppv[1,again]<-BF.correct[1,again]/(BF.correct[1,again]+BF.incorrect[1,again])
      BF.ppv[2,again]<-BF.correct[2,again]/(BF.correct[2,again]+BF.incorrect[2,again])
      BF.ppv[3,again]<-BF.correct[3,again]/(BF.correct[3,again]+BF.incorrect[3,again])
      ################################################################################     
      ## lambda
      all.lambda[1,again]<-lambda
      
      
    }#if 
    else(again=again-1)
  }
  
  ## Compute statistics from the simulations 
  # positives
  
  All[1,times]<-mean(all[1,])
  se.All[1,times]<-sd(all[1,])
  All[2,times]<-mean(all[2,])
  se.All[2,times]<-sd(all[2,])
  All[3,times]<-mean(all[3,])
  se.All[3,times]<-sd(all[3,])
  
  # tp
  tp[1,times]<-mean(correct[1,])
  se.TP[1,times]<-sd(correct[1,])
  tp[2,times]<-mean(correct[2,])
  se.TP[2,times]<-sd(correct[2,])
  tp[3,times]<-mean(correct[3,])
  se.TP[3,times]<-sd(correct[3,])
  
  #fp
  FP[1, times]<-mean(incorrect[1,])
  se.FP[1,times]<-sd(incorrect[1,])
  FP[2,times]<-mean(incorrect[2,])
  se.FP[2,times]<-sd(incorrect[2,])
  FP[3,times]<-mean(incorrect[3,])
  se.FP[3,times]<-sd(incorrect[3,])
  
  #ppv
  PPV[1, times]<-mean(ppv[1,])
  se.PPV[1,times]<- sd(ppv[1,])
  PPV[2,times]<-mean(ppv[2,])
  se.PPV[2,times]<-sd(ppv[2,])
  PPV[3,times]<-mean(ppv[3,])
  se.PPV[3,times]<-sd(ppv[3,])
  
  # shrinkage value
  L[1,times]<-mean(all.lambda)
  L[2,times]<-sd(all.lambda)
  
  #ppv BH
  BH.PPV[1, times]<-mean(BH.ppv[1,])
  BH.se.PPV[1,times]<- sd(BH.ppv[1,])
  BH.PPV[2,times]<-mean(BH.ppv[2,])
  BH.se.PPV[2,times]<-sd(BH.ppv[2,])
  BH.PPV[3,times]<-mean(BH.ppv[3,])
  BH.se.PPV[3,times]<-sd(BH.ppv[3,])
  
  #ppv BF
  BF.PPV[1, times]<-mean(BF.ppv[1,])
  BF.se.PPV[1,times]<- sd(BF.ppv[1,])
  BF.PPV[2,times]<-mean(BF.ppv[2,])
  BF.se.PPV[2,times]<-sd(BF.ppv[2,])
  BF.PPV[3,times]<-mean(BF.ppv[3,])
  BF.se.PPV[3,times]<-sd(BF.ppv[3,])
  
}


save.image(paste(etaA,"PvsN.R"))

########## figures 4 , 6 ###################
# make a folder for the results
cwd <- getwd()# CURRENT dir
newdir <- paste0(cwd,paste0("/Figure4&6",alpha))
dir.create(newdirBH)
setwd(newdir) 
#### Positives vs sample size
Method<-c(rep("ENF",length(n)),rep("MC",length(n)),rep("Shrunk MLE",length(n)))

tgc<-data.frame( Method,rep(n,3),c( All[1,], All[2,], All[3,]),c(se.All[1,],se.All[2,],se.All[3,]) )
colnames(tgc)<-c("Method","n","Positives","se")

p2<- ggplot(tgc, aes(x=n, y= Positives, colour=Method)) + 
  geom_errorbar(aes(ymin= Positives-2*se, ymax= Positives+2*se), width=0.1,colour="gray") +
  geom_line(aes(linetype=Method),alpha = 1,size=1) + #, alpha = 0.7
  geom_point(aes(shape=Method), size=3,stroke = 0.2) +
  scale_color_manual(values=c("gray13","gray21","black"))  +
  scale_linetype_manual(values=c("dotted", "solid", "dashed"))+
  scale_shape_manual( values=c(0,6,20))+
  scale_fill_manual(values=c("black","white","black"))


win.metafile(filename=paste0(etaA,paste0(alpha,"Positives.emf")), width = 7, height = 7, pointsize = 20)
#tiff(filename=paste0(etaA,paste0(alpha,"Positives.tiff")),units = "mm",res=1200, width = 86, height = 86)
p2+ theme_bw(base_size = 14)+ scale_x_continuous(breaks=n)+
  theme(panel.grid.minor = element_blank(),panel.grid.major  = element_blank())+ #legend.justification=c(1,0), legend.position="right"
  geom_hline(yintercept = alpha*p*(p-1)*0.5, size=0.25) 
dev.off()

####### TP vs sample size
Method<-c(rep("ENF",length(n)),rep("MC",length(n)),rep("Shrunk MLE",length(n)))

tgc<-data.frame( Method,rep(n,3),c( tp[1,], tp[2,], tp[3,]),c(se.TP[1,],se.TP[2,],se.TP[3,]) )
colnames(tgc)<-c("Method","n","tp","se")

p2<- ggplot(tgc, aes(x=n, y= tp, colour=Method)) + 
  geom_errorbar(aes(ymin= tp-2*se, ymax= tp+2*se), width=0.1,colour="gray") +
  geom_line(aes(linetype=Method),alpha = 1,size=1) + #, alpha = 0.7
  geom_point(aes(shape=Method), size=3,stroke = 0.2) +
  scale_color_manual(values=c("gray13","gray21","black"))  +
  scale_linetype_manual(values=c("dotted", "solid", "dashed"))+
  scale_shape_manual( values=c(0,6,20))+
  scale_fill_manual(values=c("black","white","black"))


win.metafile(filename=paste0(etaA,paste0(alpha,"TP.emf")), width = 7, height = 7, pointsize = 20)
#tiff(filename=paste0(etaA,paste0(alpha,"TP.tiff")),units = "mm",res=1200, width = 86, height = 86)
p2+ theme_bw(base_size =  14)+ scale_x_continuous(breaks=n)+
  theme(panel.grid.minor = element_blank(),panel.grid.major  = element_blank()) #legend.justification=c(1,0), legend.position="right"
dev.off()

######## PPV vs sample size
tgc<-data.frame( Method,rep(n,3),c(PPV[1,],PPV[2,],PPV[3,]),c(se.PPV[1,],se.PPV[2,],se.PPV[3,]) )
colnames(tgc)<-c("Method","n","ppv","se")

p2<-ggplot(tgc, aes(x=n, y=ppv, colour=Method)) + 
  geom_errorbar(aes(ymin= ppv-2*se, ymax= ppv+2*se), width=0.1,colour="gray") +
  geom_line(aes(linetype=Method),alpha = 1,size=1) + #, alpha = 0.7
  geom_point(aes(shape=Method), size=3,stroke = 0.2) +
  scale_color_manual(values=c("gray13","gray21","black"))  +
  scale_linetype_manual(values=c("dotted", "solid", "dashed"))+
  scale_shape_manual( values=c(0,6,20))+
  scale_fill_manual(values=c("black","white","black"))

win.metafile(filename=paste0(etaA,paste0(alpha,"PPV.emf")), width = 7, height = 7, pointsize = 20)
#tiff(filename=paste0(etaA,paste0(alpha,"PPV.tiff")),units = "mm",res=1200, width = 86, height = 86)
p2+ theme_bw(base_size = 14)+ scale_x_continuous(breaks=n)+
  theme(panel.grid.minor = element_blank(),panel.grid.major  = element_blank()) #legend.justification=c(1,0), legend.position="right"
dev.off()



###### PPV BH vs sample size
tgc<-data.frame( Method,rep(n,3),c(BH.PPV[1,],BH.PPV[2,],BH.PPV[3,]),c(BH.se.PPV[1,],BH.se.PPV[2,],BH.se.PPV[3,]) )
colnames(tgc)<-c("Method","n","ppv","se")

p2<-ggplot(tgc, aes(x=n, y=ppv, colour=Method)) +
  geom_errorbar(aes(ymin= ppv-2*se, ymax= ppv+2*se), width=0.1,colour="gray") +
  geom_line(aes(linetype=Method),alpha = 1,size=1) + #, alpha = 0.7
  geom_point(aes(shape=Method), size=3,stroke = 0.2) +
  scale_color_manual(values=c("gray13","gray21","black"))  +
  scale_linetype_manual(values=c("dotted", "solid", "dashed"))+
  scale_shape_manual( values=c(0,6,20))+
  scale_fill_manual(values=c("black","white","black"))

win.metafile(filename=paste0(etaA,paste0(alpha,"BHPPV.emf")), width = 7, height = 7, pointsize = 20)
#tiff(filename=paste0(etaA,paste0(alpha,"BHPPV.tiff")),units = "mm",res=1200, width = 86, height = 86)
p2+ theme_bw(base_size = 14)+ scale_x_continuous(breaks=n)+
  theme(panel.grid.minor = element_blank(),panel.grid.major  = element_blank()) #legend.justification=c(1,0), legend.position="right"
dev.off()


#### PPV BF vs sample size
tgc<-data.frame( Method,rep(n,3),c(BF.PPV[1,],BF.PPV[2,],BF.PPV[3,]),c(BF.se.PPV[1,],BF.se.PPV[2,],BF.se.PPV[3,]) )
colnames(tgc)<-c("Method","n","ppv","se")

p2<-ggplot(tgc, aes(x=n, y=ppv, colour=Method)) +
  geom_errorbar(aes(ymin= ppv-2*se, ymax= ppv+2*se), width=0.1,colour="gray") +
  geom_line(aes(linetype=Method),alpha = 1,size=1) + #, alpha = 0.7
  geom_point(aes(shape=Method), size=3,stroke = 0.2) +
  scale_color_manual(values=c("gray13","gray21","black"))  +
  scale_linetype_manual(values=c("dotted", "solid", "dashed"))+
  scale_shape_manual( values=c(0,6,20))+
  scale_fill_manual(values=c("black","white","black"))


win.metafile(filename=paste0(etaA,paste0(alpha,"BFPPV.emf")), width = 7, height = 7, pointsize = 20)
#tiff(filename=paste0(etaA,paste0(alpha,"BFPPV.tiff")),units = "mm",res=1200, width = 86, height = 86)
p2+ theme_bw(base_size = 14)+ scale_x_continuous(breaks=n)+
  theme(panel.grid.minor = element_blank(),panel.grid.major  = element_blank()) #legend.justification=c(1,0), legend.position="right"
dev.off()



####  False positives vs sample size
tgc<-data.frame( Method,rep(n,3),c(FP[1,],FP[2,],FP[3,]),c(se.FP[1,],se.FP[2,],se.FP[3,]) )
colnames(tgc)<-c("Method","n","FP","se")

p2<-ggplot(tgc, aes(x=n, y=FP, colour= Method)) + 
  geom_errorbar(aes(ymin= FP-2*se, ymax= FP+2*se), width=0.1,colour="gray") +
  geom_line(aes(linetype=Method),alpha = 1,size=1) + #, alpha = 0.7
  geom_point(aes(shape=Method), size=3,stroke = 0.2) +
  scale_color_manual(values=c("gray13","gray21","black"))  +
  scale_linetype_manual(values=c("dotted", "solid", "dashed"))+
  scale_shape_manual( values=c(0,6,20))+
  scale_fill_manual(values=c("black","white","black"))

win.metafile(filename=paste0(etaA,paste0(alpha,"FP.emf")), width = 7, height = 7, pointsize = 20)
#tiff(filename=paste0(etaA,paste0(alpha,"FP.tiff")),units = "mm",res=1200, width = 86, height = 86)
p2+ theme_bw(base_size = 14)+ scale_x_continuous(breaks=n)+
  theme(panel.grid.minor = element_blank(),panel.grid.major  = element_blank()) #legend.justification=c(1,0), legend.position="right"
dev.off()

setwd(cwd)

#****************************************************************
# 2) Compare densities and histograms
#   We will study the probability density under the null hypothesis (no partial correlation)
#   by using the standard and the novel shrunk density
#************************************************************************
# Parameters       
# * p = Number of variables (e.g. genes)
# * n = Number of samples  
# * number= Montecarlo iterations
# * rep = the simulation is repeated (for fixed p,n)
# 
# * etaA = proportion of TP
# * alpha = significance level
# * eta0 = 1-etaA
#************************************************************************

rm(list=ls())

library(GeneNet)
library(stats4)
library(Hmisc)
library(ggplot2)


set.seed(123)
source("shrinkagefunctions.R") # Check that shrinkagefunctions.R in the same directory


#***************************
# Initialize parameters
p<-100 # genes' number
etaA=0.00 #proportion true positive correlations 
n<-c(20) # sample size
alpha=0.05
number<-5 # MC iterations

## Simulated pcorr and data
sim.pcor<-ggm.simulate.pcor(p, etaA)
sim.data<- ggm.simulate.data( n , sim.pcor)
temp=sim.pcor[lower.tri(sim.pcor, diag = FALSE) ]
TP=which(temp!=0 )
TN=which(temp==0 )
GGM <- pcor.shrink(sim.data, verbose=TRUE) # OPTIMAL
lambda<-attr(GGM, "lambda") 
kappa(GGM, exact=TRUE, norm='2' )
r<-sm2vec(GGM)

# parameters of the mixture distribution used to compute p-values etc.
# c <- fdrtool(sm2vec(GGM), statistic="correlation")
# c$param


## GENENET methods (empirical null fitting + std density)
## WATCH OUT:  the test must be in node order. Sort by node 2 and by node 1
ENF.test<-network.test.edges(GGM, fdr=TRUE, plot=FALSE)
ENF.test<-ENF.test[order(ENF.test$node1,ENF.test$node2),]

#######################################
# MLE (shrunk density)
pval.shrunk<-p.shrunk(r,p,n,lambda)
# MLE (standard/std density)
# pval.std<-p.shrunk(r,p,n,0)
pval.standard <- p.standard ( r, p,n)

# MC
p.monte<-p.montecarlo(r,number,p,n,lambda)

# Estimate the critical value (cv.monte) at alpha with Monte Carlo
r.monte<-matrix(NA, length(r),number)
for (i in 1:number){
  r.data<-ggm.simulate.data(n,diag(p))
  r.monte.GGM<-ggm.estimate.pcor(r.data,lambda=lambda )
  r.monte[,i]<-sm2vec(r.monte.GGM)
}
cv.values<-matrix(NA,1,ncol(r.monte))
for (i in 1:ncol(r.monte)){
  # 2 tails
  cv<-apply( -abs(r.monte),2, function(x) quantile(x,  probs = c((alpha/2)) ) )  
}
cv.monte<-mean(cv)   
cv.sd<-sd(cv)

########## figure 3 ###################
# make a folder for the results
cwd <- getwd()# CURRENT dir
newdir <- paste0(cwd,paste0("/Figure3",alpha))
dir.create(newdirBH)
setwd(newdir)

win.metafile(filename=paste0(etaA,paste0(alpha,"HISTenf.wmf")), width = 10, height = 10, pointsize = 20)
#tiff(filename=paste0(etaA,paste0(n,"HIST.tiff")),units = "mm",res=1200, width = 86, height = 86)
hist(p.monte,20, col=rgb(1,1,1,0.5),main = c(),xlab="p-values",ylim=c(0,max(320)),cex.axis=1,cex.lab=1)
hist(pval.shrunk,20, col=rgb(0.7,0.7,0.7,1),add=T)
hist(ENF.test$pval,20, col=rgb(0.3,0.3,0.3,1),add=T)
legend("bottomright", c("MC", "Shrunk MLE", "ENF"), fill=c("white","gray","gray21"), cex = 1, bg="white")#,bty="n"
dev.off()

win.metafile(filename=paste0(etaA,paste0(alpha,"HISTstd.wmf")), width = 10, height = 10, pointsize = 20)
#tiff(filename=paste0(etaA,paste0(n,"HISTstd.tiff")),units = "mm",res=1200, width = 86, height = 86)
hist(p.monte,20, col=rgb(1,1,1,0.5),main = c(),xlab="p-values",ylim=c(0,max(320)),cex.axis=1,cex.lab=1)
hist(pval.shrunk,20, col=rgb(0.7,0.7,0.7,1),add=T)
hist(pval.standard,20, col=rgb(0.3,0.3,0.3,1),add=T)
legend("bottomright", c("MC", "Shrunk MLE", "Standard MLE"), fill=c("white","gray","gray21"), cex = 1, bg="white")#,bty="n"
dev.off()

###### Estimate degrees of freedom  k
#### Shrunk
sim.pcor.MC<-ggm.simulate.pcor(p, 0)
sim.data.MC<- ggm.simulate.data( n , sim.pcor.MC)
GGM.MC <- pcor.shrink(sim.data.MC,lambda , verbose=TRUE) 
kappa(GGM.MC, exact=TRUE, norm='2' )
r.MC<-sm2vec(GGM.MC)
nlogL.shrunk <- function(k) {
  
  density.shrunk <- function(r.MC) {(  ((1-lambda)^2-r.MC^2) ^((k-3)*0.5)  )/( beta(0.5, 0.5*(k-1))*(1-lambda)*(1-lambda)^(k-3) )}
  
  f<-density.shrunk(r.MC)
  -sum(log(f))
}
### neg log Likelihood
k.fit.shrunk <-  mle(nlogL.shrunk, start = list(k = 100), method = "L-BFGS-B", lower = c(20),
                     upper = c(Inf))
summary(k.fit.shrunk)


######## Standard
nlogL.std <- function(k) {
  density.std <- function(r.MC) {(  ((1)^2-r.MC^2) ^((k-3)*0.5)  )/( beta(0.5, 0.5*(k-1)))}
  
  f<-density.std(r.MC)
  
  -sum(log(f))
}

k.fit.std <-  mle(nlogL.std, start = list(k = 50), method = "L-BFGS-B", lower = c(20),
                  upper = c(Inf))
summary(k.fit.std)

######################################
## Plot the difference of the densities (Standard and Shrunk ) with Maximum Likelihood (MLE). 
density.shrunk <- function(x) {(  ((1-lambda)^2-x^2) ^(( k.fit.shrunk@coef[1]-3)*0.5)  )/( beta(0.5, 0.5*( k.fit.shrunk@coef[1]-1))*(1-lambda)*(1-lambda)^( k.fit.shrunk@coef[1]-3) )}
density.std <- function(x) {(  ((1)^2-x^2) ^((k.fit.std@coef[1]-3)*0.5)  )/( beta(0.5, 0.5*(k.fit.std@coef[1])-1))}
x<-c(-100:100)*(1-lambda)/100

win.metafile(filename=paste0(n,"densities.wmf"),  width = 10, height = 10, pointsize = 20)
#tiff(filename=paste0(n,"densities.tiff"),units = "mm",res=1200, width = 86, height = 86)
plot(x, density.std(x),type="l",col="gray",pch=20,lwd = 1,xlab="pcorr", ylab="Prob density",cex=1, cex.lab=1, cex.axis=1) 
lines(x, density.shrunk(x),lwd = 1,col="black", lty="dashed") 
abline(v =0)
abline(h =0)

legend("topright",
       legend = c("Standard density MLE", "Shrunk density MLE"), 
       col = c("gray", 
               "black"), 
       pch = c("-","-"),
       bty = "n", 
       pt.cex = 1, 
       cex = 1, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))
abline(v=c(cv.monte,-cv.monte), col=c("gray","gray"), lty=c(2,2), lwd=c(2,2))

dev.off()



win.metafile(filename=paste0(n,"diff_densities.wmf"),  width = 10, height = 10, pointsize = 20)
#tiff(filename=paste0(n,"diff_densities.tiff"),units = "mm",res=1200, width = 86, height = 86)
plot(x, density.std(x) - density.shrunk(x),type="l",pch=20,lwd = 2.5, ylab="Prob density",xlab="pcorr",cex=1, cex.lab=1, cex.axis=1 ) 
abline(v =0)
abline(h =0)
legend("bottomright", 
       legend = c("pdf difference", paste("critical value at",alpha)), 
       col = c("black", 
               "black"), 
       pch = c("-","."),
       bg = "white", 
       pt.cex = c(1,4), 
       cex = 1, 
       text.col = "black", 
       horiz = F , 
       inset = c(0, 0.1))
abline(v=c(cv.monte,-cv.monte), col=c("black","black"),lty=c("dotted", "dotted") , lwd=c(2,2))
dev.off()

###########################################
## Plot the difference of the densities (Standard and Shrunk ) with ENF and MLE.
# # parameters of the mixture distribution used to compute p-values etc.
c <- fdrtool(sm2vec(GGM), statistic="correlation")
c$param
data.frame( k.fit.shrunk@coef[1], k.fit.std@coef[1], c$param[5]) 

density.shrunk <- function(x) {(  ((1-lambda)^2-x^2) ^(( k.fit.shrunk@coef[1]-3)*0.5)  )/( beta(0.5, 0.5*( k.fit.shrunk@coef[1]-1))*(1-lambda)*(1-lambda)^( k.fit.shrunk@coef[1]-3) )}
density.std <- function(x) {(  ((1)^2-x^2) ^((c$param[5]-3)*0.5)  )/( beta(0.5, 0.5*(c$param[5])-1))}
x<-c(-100:100)*(1-lambda)/100


win.metafile(filename=paste0(n,"densities_ENF.emf"),  width = 10, height = 10, pointsize = 20)
plot(x, density.std(x),type="l",col="grey", pch=20,lwd = 1, xlab="pcorr", ylab="Prob density")
lines(x, density.shrunk(x),lwd = 1,col="black", lty="dashed")
abline(v =0)
abline(h =0)
legend("topright",
       legend = c("Standard density ENF", "Shrunk density MLE"),
       col = c("grey","black"),
       pch = c("-","-"),
       bty = "n",
       pt.cex = 1,
       cex = 1,
       text.col = "black",
       horiz = F ,
       inset = c(0.1, 0.1))
abline(v=c(cv.monte,-cv.monte), col=c("gray","gray"), lty=c(0.75,0.75), lwd=c(0.75,0.75))
dev.off()


win.metafile(filename=paste0(n,"diff_densities_ENF.emf"),  width = 10, height = 10, pointsize = 20)
#tiff(filename=paste0(n,"diff_densities.tiff"),units = "mm",res=1200, width = 86, height = 86)
plot(x, density.std(x) - density.shrunk(x),type="l",col="black",pch=20,lwd = 2,xlab="pcorr",cex=1, cex.lab=1, cex.axis=1, ylab="Prob density" ) 
abline(v =0)
abline(h =0)
legend("bottomright", 
       legend = c("pdf difference", paste("critical value at",alpha)), 
       col = c("black", 
               "black"), 
       pch = c("-","."),
       bg = "white", 
       pt.cex = c(1,4), 
       cex = 1, 
       text.col = "black", 
       horiz = F , 
       inset = c(0, 0.1))
abline(v=c(cv.monte,-cv.monte), col=c("black","black"),lty=c("dotted", "dotted"), lwd=c(2,2))

dev.off()

setwd(cwd)
#**************************************************************
##  Heatmap of False Positive Rates under the null Hypothesis
#   We will compare the False Positive Rates under the null hypothesis (no partial correlation)
#   with respect to Monte Carlo (the gold standard)
#**************************************************************
rm(list = ls())

library(GeneNet)
library(reshape)
library(ggplot2)
library(stats4)
library(Hmisc)

set.seed(123)
source("shrinkagefunctions.R") # Check that shrinkagefunctions.R in the same directory



#***************************
# Initialize parameters
p<-c( 10,20,40,60,80,100,200,300,400 ) # num of genes   
n<-c(3,5,10,15,20,30,40,70,100 ) # num of samples  

etaA=0.00  # proportion of  TP
alpha<-0.05 # significance level
lambda<-0.3
number= 1 # MC iterations
all.lambda<-matrix(NA,length(n),length(p)) 
diff.mc.shrunk<-matrix(NA,length(n),length(p))
diff.mc.ENF<-matrix(NA,length(n),length(p))

for (times in 1:length(p)){
  
  for (s.size in 1:length(n)){
    
    # ENF    
    true.pcor<-ggm.simulate.pcor(p[times], etaA)
    positives<-which(sm2vec(true.pcor)!=0)
    # same ggm, several simulation
    null.data<-ggm.simulate.data(n[s.size],true.pcor)
    null.data2<-ggm.simulate.data(n[s.size],true.pcor)
    null.data3<-ggm.simulate.data(n[s.size],true.pcor)
    null.data4<-ggm.simulate.data(n[s.size],true.pcor)
    null.data5<-ggm.simulate.data(n[s.size],true.pcor)
    null.data6<-ggm.simulate.data(n[s.size],true.pcor)
    null.data7<-ggm.simulate.data(n[s.size],true.pcor)
    null.data8<-ggm.simulate.data(n[s.size],true.pcor)
    null.data9<-ggm.simulate.data(n[s.size],true.pcor)
    null.data10<-ggm.simulate.data(n[s.size],true.pcor)
    
    GGM<-pcor.shrink(  null.data,lambda  )
    GGM2<-pcor.shrink(  null.data2,lambda  )
    GGM3<-pcor.shrink(  null.data3,lambda  )
    GGM4<-pcor.shrink(  null.data4,lambda  )
    GGM5<-pcor.shrink(  null.data5,lambda  )
    GGM6<-pcor.shrink(  null.data6,lambda  )
    GGM7<-pcor.shrink(  null.data7,lambda  )
    GGM8<-pcor.shrink(  null.data8,lambda  )
    GGM9<-pcor.shrink(  null.data9,lambda  )
    GGM10<-pcor.shrink(  null.data10,lambda  )
    
    r<-sm2vec(GGM)
    r2<-sm2vec(GGM2)
    r3<-sm2vec(GGM3)
    r4<-sm2vec(GGM4)
    r5<-sm2vec(GGM5)
    r6<-sm2vec(GGM6)
    r7<-sm2vec(GGM7)
    r8<-sm2vec(GGM8)
    r9<-sm2vec(GGM9)
    r10<-sm2vec(GGM10)
    #lambda=attr(GGM,"lambda")
    
    test.ENF <- network.test.edges(GGM,plot=FALSE)
    test.ENF<-test.ENF[order(test.ENF$node1,test.ENF$node2),]
    # test.ENF$pval
    
    test.ENF2 <- network.test.edges(GGM2,plot=FALSE)
    test.ENF2<-test.ENF2[order(test.ENF2$node1,test.ENF2$node2),]
    # test.ENF2$pval
    
    test.ENF3 <- network.test.edges(GGM3,plot=FALSE)
    test.ENF3<-test.ENF3[order(test.ENF3$node1,test.ENF3$node2),]
    # test.ENF3$pval
    
    test.ENF4 <- network.test.edges(GGM4,plot=FALSE)
    test.ENF4<-test.ENF4[order(test.ENF4$node1,test.ENF4$node2),]
    # test.ENF4$pval
    
    test.ENF5 <- network.test.edges(GGM5,plot=FALSE)
    test.ENF5<-test.ENF5[order(test.ENF5$node1,test.ENF5$node2),]
    # test.ENF5$pval
    test.ENF6 <- network.test.edges(GGM6,plot=FALSE)
    test.ENF6<-test.ENF5[order(test.ENF6$node1,test.ENF6$node2),]
    
    test.ENF7 <- network.test.edges(GGM7,plot=FALSE)
    test.ENF7<-test.ENF5[order(test.ENF7$node1,test.ENF7$node2),]
    
    test.ENF8 <- network.test.edges(GGM8,plot=FALSE)
    test.ENF8<-test.ENF5[order(test.ENF8$node1,test.ENF8$node2),]
    
    test.ENF9 <- network.test.edges(GGM9,plot=FALSE)
    test.ENF9<-test.ENF5[order(test.ENF9$node1,test.ENF9$node2),]
    
    test.ENF10 <- network.test.edges(GGM10,plot=FALSE)
    test.ENF10<-test.ENF5[order(test.ENF10$node1,test.ENF10$node2),]
    
    
    #Shrunk
    ###  simulate null for shrunk density
    cat("ShrunkMLE")
    
    # ### Shrunk MLE
    pval.shrunk<-p.shrunk(r,p[times],n[s.size],lambda)
    pval.shrunk2<-p.shrunk(r2,p[times],n[s.size],lambda)
    pval.shrunk3<-p.shrunk(r3,p[times],n[s.size],lambda)
    pval.shrunk4<-p.shrunk(r4,p[times],n[s.size],lambda)
    pval.shrunk5<-p.shrunk(r5,p[times],n[s.size],lambda)
    pval.shrunk6<-p.shrunk(r6,p[times],n[s.size],lambda)
    pval.shrunk7<-p.shrunk(r7,p[times],n[s.size],lambda)
    pval.shrunk8<-p.shrunk(r8,p[times],n[s.size],lambda)
    pval.shrunk9<-p.shrunk(r9,p[times],n[s.size],lambda)
    pval.shrunk10<-p.shrunk(r10,p[times],n[s.size],lambda)
    #######
    # MC
    #######
    cat("MC")
    p.monte<-p.montecarlo(r,number,p[times],n[s.size],lambda)
    p.monte2<-p.montecarlo(r2,number,p[times],n[s.size],lambda)
    p.monte3<-p.montecarlo(r3,number,p[times],n[s.size],lambda)
    p.monte4<-p.montecarlo(r4,number,p[times],n[s.size],lambda)
    p.monte5<-p.montecarlo(r5,number,p[times],n[s.size],lambda)
    p.monte6<-p.montecarlo(r6,number,p[times],n[s.size],lambda)
    p.monte7<-p.montecarlo(r7,number,p[times],n[s.size],lambda)
    p.monte8<-p.montecarlo(r8,number,p[times],n[s.size],lambda)
    p.monte9<-p.montecarlo(r9,number,p[times],n[s.size],lambda)
    p.monte10<-p.montecarlo(r10,number,p[times],n[s.size],lambda)
    ###############
    # Monte Carlo - other methods
    # rows=n, col=p
    
    all.lambda[s.size,times]<-lambda
    c<-p[times]*(p[times]-1)*0.5  
    
    
    # false positives    
    diff.mc.shrunk[ s.size,times]<-(1/c)* mean( c(  sum( ! which( p.monte<=alpha)%in%positives) - sum( ! which( pval.shrunk<=alpha)%in%positives),
                                                    sum( ! which( p.monte2<=alpha)%in%positives) - sum( ! which( pval.shrunk2<=alpha)%in%positives),
                                                    sum( ! which( p.monte3<=alpha)%in%positives) - sum( ! which( pval.shrunk3<=alpha)%in%positives),
                                                    sum( ! which( p.monte4<=alpha)%in%positives) - sum( ! which( pval.shrunk4<=alpha)%in%positives),
                                                    sum( ! which( p.monte5<=alpha)%in%positives) - sum( ! which( pval.shrunk5<=alpha)%in%positives),
                                                    sum( ! which( p.monte6<=alpha)%in%positives) - sum( ! which( pval.shrunk6<=alpha)%in%positives),
                                                    sum( ! which( p.monte7<=alpha)%in%positives) - sum( ! which( pval.shrunk7<=alpha)%in%positives),
                                                    sum( ! which( p.monte8<=alpha)%in%positives) - sum( ! which( pval.shrunk8<=alpha)%in%positives),
                                                    sum( ! which( p.monte9<=alpha)%in%positives) - sum( ! which( pval.shrunk9<=alpha)%in%positives),
                                                    sum( ! which( p.monte10<=alpha)%in%positives) - sum( ! which( pval.shrunk10<=alpha)%in%positives)), na.rm = TRUE)  
    
    
    diff.mc.ENF[ s.size,times]<- (1/c)* mean( c(  sum( ! which( p.monte<=alpha)%in%positives) - sum( ! which( test.ENF$pval<=alpha)%in%positives),
                                                  sum( ! which( p.monte2<=alpha)%in%positives) - sum( ! which( test.ENF2$pval<=alpha)%in%positives),
                                                  sum( ! which( p.monte3<=alpha)%in%positives) - sum( ! which( test.ENF3$pval<=alpha)%in%positives),
                                                  sum( ! which( p.monte4<=alpha)%in%positives) - sum( ! which( test.ENF4$pval<=alpha)%in%positives),
                                                  sum( ! which( p.monte5<=alpha)%in%positives) - sum( ! which( test.ENF5$pval<=alpha)%in%positives),
                                                  sum( ! which( p.monte6<=alpha)%in%positives) - sum( ! which( test.ENF6$pval<=alpha)%in%positives),
                                                  sum( ! which( p.monte7<=alpha)%in%positives) - sum( ! which( test.ENF7$pval<=alpha)%in%positives),
                                                  sum( ! which( p.monte8<=alpha)%in%positives) - sum( ! which( test.ENF8$pval<=alpha)%in%positives),
                                                  sum( ! which( p.monte9<=alpha)%in%positives) - sum( ! which( test.ENF9$pval<=alpha)%in%positives),
                                                  sum( ! which( p.monte10<=alpha)%in%positives) - sum( ! which( test.ENF10$pval<=alpha)%in%positives)), na.rm = TRUE)  
    
  }
  
}

save.image(paste(etaA,"heatmap_pn.R"))

########## figure 5 ###################
# make a folder for the results
cwd <- getwd()# CURRENT dir
newdir <- paste0(cwd,paste0("/Figure5",alpha))
dir.create(newdirBH)
setwd(newdir)
# limit for colours
min.colour<-min(diff.mc.ENF,diff.mc.shrunk,na.rm = TRUE)
max.colour<-max(diff.mc.ENF,diff.mc.shrunk,na.rm = TRUE)
num.col<-3

## MC - ENF 

# cols= samples size, rows are p (number of genes)
num_genes <- paste(rep("p", length(p)), p, sep="=") # rows
sample_size <- paste(rep("n", length(n)), n, sep="=") # columns
mc.ENF <- data.frame(genes = num_genes, 
                     diff.mc.ENF,nrow = length(p), ncol =length(n) )

#mc.ENF$genes <- factor(mc.ENF$genes, levels = sort(unique(mc.ENF$genes)))

names(mc.ENF )[2:(length(n)+1)] <- sample_size 
mc.ENF<-mc.ENF[,-c(ncol(mc.ENF)-1,ncol(mc.ENF))]
df_heatmap.gnet <- melt(mc.ENF , id.vars = "genes")
names(df_heatmap.gnet)[2:3] <- c("sample", "P")
head(df_heatmap.gnet)

# lock in factor level 
df_heatmap.gnet$genes <- factor(df_heatmap.gnet$genes, levels = unique(mc.ENF$genes))


win.metafile(paste(etaA,"heatmap_ENF.wmf"), width = 7, height = 7, pointsize = 20)
#tiff(filename=paste0(etaA,"heatmap_ENF.tiff"),units = "mm",res=1200, width = 86, height = 86)
ggplot(df_heatmap.gnet, aes(x=sample, y=genes  )) +
  geom_tile(aes(fill = P), color = "white") +
  # scale_fill_gradientn(colours = terrain.colors(10),limits=c(min.colour,max.colour))+
  # scale_fill_gradientn(colours =rainbow(num.col) , limits=c(min.colour,max.colour))+
  scale_fill_gradient2(low = "black", mid = "white", high = "black" ,  midpoint = 0,limits=c(min.colour,max.colour))+
  
  geom_text(aes(label = round(df_heatmap.gnet$P, 2)),size=4) +
  # scale_colour_gradient2(low = "red", mid = "white",
  #                        high = "steelblue", midpoint = 0, space = "Lab",
  #                        na.value = "grey50", guide = "colourbar")+
  ylab("") +
  xlab("") +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1,size=12),
        axis.text.y = element_text(size=12) )+
  labs(fill = "FPR ENF- FPR MC") # +

dev.off()

# MC- Shrunk
num_genes <- paste(rep("p", length(p)), p, sep="=") # rows
sample_size <- paste(rep("n", length(n)), n, sep="=") # columns


mc.Shrunk <- data.frame(genes = num_genes, 
                        diff.mc.shrunk,nrow = length(p), ncol =length(n) )

names(mc.Shrunk )[2:(length(n)+1)] <- sample_size 

mc.Shrunk<-mc.Shrunk[,-c(ncol(mc.Shrunk)-1,ncol(mc.Shrunk))]

df_heatmap <- melt(mc.Shrunk , id.vars = "genes")
names(df_heatmap)[2:3] <- c("sample", "P")
head(df_heatmap)

# lock in factor level 
df_heatmap$genes <- factor(df_heatmap$genes, levels = unique(mc.Shrunk$genes))

win.metafile(paste(etaA,"heatmap_Shrunk.wmf"), width = 7, height = 7, pointsize = 20)
#tiff(filename=paste0(etaA,"heatmap_Shrunk.tiff"),units = "mm",res=1200, width = 86, height = 86)
ggplot(df_heatmap, aes(x=sample, y=genes  )) +
  geom_tile(aes(fill = P), color = "white") +
  geom_text(aes(label = round(df_heatmap$P, 2)),size=4) +
  # scale_fill_gradientn(colours = terrain.colors(10),limits=c(min.colour,max.colour))+
  # scale_fill_gradientn(colours =rainbow(num.col) , limits=c(min.colour,max.colour))+
  scale_fill_gradient2(low = "black", mid = "white", high = "black" ,  midpoint = 0,limits=c(min.colour,max.colour))+
  
  ylab("") +
  xlab("") +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1,size=12),
        axis.text.y = element_text(size=12) )+
  labs(fill = "FPR Shrunk- FPR MC")
dev.off()

setwd(cwd)


