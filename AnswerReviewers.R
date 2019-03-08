#|**********************************************************************;
#
#                     Answer to the reviewers
#                 Manuscript ID: BIOINF-2018-2426
#
#  * Project: EXACT HYPOTHESIS TESTING FOR SHRINKAGE BASED GGM (Bernal et al).
#  * Author            : Victor Bernal*, Rainer Bischoff, Victor Guryev, Marco Grzegorczyk, Peter Horvatovich
#  * Date created      : 28-Jan-2019
#  * Submitted on: 15-Nov-2018
#  * Review on: 18-Jan-2019
#***********************************
# Description:
# This script reproduces the computations necessary to backup our
# answer to the review.
#***************************************
# Prerequisite R library (CRAN)
# requires R libraries "GeneNet" , "stats4", "ggplot2".
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
#  For computation of empirical p-value via Monte Carlo see  [Martinez, W. L., & Martinez, A. R. (2007). Computational statistics handbook with MATLAB. Chapman and Hall/CRC.] 
#  The standard probability density is presented in [Fisher,R.A. (1924) The distribution of the partial correlation coefficient. Metron, 3,329-332.] 
#  The simulation of data, as well as estimation with the optimal shrinkage (lambda) is done with
#  [Schäfer,J. and Strimmer,K. (2005a),  GeneNet 1.2.13. CRAN.]
#
#************************************************************************;



#***********************************************************************************
# Answer to Reviewer 1 comment 2
# We generated an averaged Histogram and tested H0: the first bin's height 
# are the same as the gold standard MC at 5% level. 
# number of simulations: rep
#***********************************************************************************
rm(list=ls())

library(GeneNet)
library(stats4)
library(ggplot2)


set.seed(123)
source("H:/Semester7(15032017)/RUG/Shrinkage/shrinkagefunctions.R")
#setwd("H:/Semester7(15032017)/RUG/Shrinkage")

#***************************
# Initialize parameters

p <- 100 # genes' number
etaA <- 0.00 #proportion true positive correlations 
n <- c(20) # sample size
alpha <- 0.05
number <- 15 # MC iterations
breaks <- c(0 ,c(1:20)/20)
length(breaks)

## Repeat the simulation (rep) to take the average 
rep<- 25
counts.shrunk <- matrix(Inf, length(breaks)-1, rep)
counts.ENF <- matrix(Inf, length(breaks)-1, rep)
counts.monte <- matrix(Inf, length(breaks)-1, rep)
counts.standard <- matrix(Inf, length(breaks)-1, rep)
sim.pcor<-ggm.simulate.pcor(p, etaA)


for (i in 1:rep){
  
  sim.data<- ggm.simulate.data( n , sim.pcor)
  temp <- sim.pcor[lower.tri(sim.pcor, diag = FALSE) ]
  TP <- which(temp!=0 )
  TN <- which(temp==0 )
  GGM <- pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
  lambda <- attr(GGM, "lambda") 
  
    while (lambda == 1) {
      sim.data <- ggm.simulate.data( n , sim.pcor)
      temp <- sim.pcor[lower.tri(sim.pcor, diag = FALSE) ]
      TP <- which(temp!=0 )
      TN <- which(temp==0 )
      GGM <- pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
      lambda <- attr(GGM, "lambda") 
      print("TRY AGAIN")
    } 
  
    r<-sm2vec(GGM)
    
    #######################################
    # MLE (shrunk density [this work])
    pval.shrunk <- p.shrunk(r,p,n,lambda)
    
    pval.standard <- p.standard ( r, p,n)

    ## GENENET methods (empirical null fitting + std density [Fisher 1921])
    ## WATCH OUT:  the test must be in node order. Sort by node 2 and by node 1
    ENF.test<-network.test.edges(GGM, fdr=TRUE, plot=FALSE, verbose =FALSE)
    ENF.test<-ENF.test[order(ENF.test$node1,ENF.test$node2),]
  
    # Monte Carlo (The gold standard)
    p.monte <- p.montecarlo(r,number,p,n,lambda)
    
    
    
    counts.shrunk[, i] <- hist(pval.shrunk, breaks = breaks,  plot = FALSE)$counts
    counts.standard[, i] <- hist(pval.standard, breaks = breaks ,  plot = FALSE)$counts
    
    counts.ENF[, i] <- hist(ENF.test$pval, breaks = breaks,  plot = FALSE)$counts
    counts.monte[, i] <- hist(p.monte, breaks = breaks ,  plot = FALSE)$counts
    
     
}


# #****** Averaged Histogram
# histgrm <- c(rep("Shrunk MLE", length(breaks)-1 ),rep("ENF", length(breaks)-1), rep("MC", length(breaks)-1))
# 
# 
# 
# tgc<-data.frame( histgrm,
#                  c( rowMeans(counts.shrunk, na.rm = FALSE), 
#                     rowMeans(counts.ENF, na.rm = FALSE),
#                     rowMeans(counts.monte, na.rm = FALSE)),
#                  
#                  c( apply(counts.shrunk,1,sd),
#                     0* apply(counts.ENF,1,sd),
#                     0*apply(counts.monte,1,sd) )/sqrt(n) ,
#                  
#                  rep(hist(p.monte, length(breaks)-1)$mids,3)
#                     
#                  )
# colnames(tgc)<-c("Method", "Average","se","pvalue")
# tgc$Method <-factor(tgc$Method)
# 
# 
# ph <- ggplot(tgc, aes(x = pvalue , y = Average, fill = Method)) +
#   geom_bar(stat="identity", position = "identity", color="black" )+
#   geom_errorbar(aes(ymin= Average-2*se, ymax= Average+2*se, width=0.025 )) +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())+
#   theme_bw(base_size = 14)+ scale_x_continuous(breaks = breaks)+
#   theme(panel.grid.minor = element_blank(), panel.grid.major  = element_blank())+ theme(legend.position=c(0.8,0.2)) 
# 
# win.metafile(filename=paste0(n,paste0(p,"aveHist.emf")), width = 8, height = 8, pointsize = 20)
# ph
# dev.off()

library(Hmisc)

# ENF
win.metafile(filename=paste0(n,paste0(p,"aveHist.emf")), width = 8, height = 8, pointsize = 20)

mp<- barplot(rowMeans(counts.monte, na.rm = FALSE) ,length(breaks)-1, col = rgb(1,1,1,0.5), main = c(),xlab="p-values", 
     ylim=c(0,max(counts.monte)) ) #,cex.axis = 0.75,cex.lab = 0.75
barplot(rowMeans(counts.shrunk , na.rm = FALSE),length(breaks)-1, col = rgb(0.7,0.7,0.7,1), add=T)
barplot(rowMeans(counts.ENF, na.rm = FALSE),length(breaks)-1, col = rgb(0.35,0.35,0.35,1), add=T)
    y <- rowMeans(counts.ENF, na.rm = FALSE)
    delta <- apply(counts.ENF,1,sd)/sqrt(n)
    
    errbar( mp , y, y + 2*delta, y - 2*delta , col="black", add=T , lwd = 3 , errbar.col = "black")
legend("bottomright", c("MC", "Shrunk MLE", "ENF"), fill=c(rgb(1,1,1,0.5),rgb(0.7,0.7,0.7,1),rgb(0.35,0.35,0.35,1)),
       cex = 0.85, bg="white")#,bty="n"
dev.off()


# standard
win.metafile(filename=paste0(n,paste0(p,"aveHist2.emf")), width = 8, height = 8, pointsize = 20)

mp<- barplot(rowMeans(counts.monte, na.rm = FALSE) ,length(breaks)-1, col = rgb(1,1,1,0.5), main = c(),xlab="p-values", 
             ylim=c(0,max(counts.monte)) ) #,cex.axis = 0.75,cex.lab = 0.75
barplot(rowMeans(counts.shrunk , na.rm = FALSE),length(breaks)-1, col = rgb(0.7,0.7,0.7,1), add=T)
barplot(rowMeans(counts.standard, na.rm = FALSE),length(breaks)-1, col = rgb(0.35,0.35,0.35,1), add=T)
    y <- rowMeans(counts.standard, na.rm = FALSE)
    delta <- apply(counts.standard,1,sd)/sqrt(n)

    errbar( mp , y, y + 2*delta, y - 2*delta , col="black", add=T , lwd = 3 , errbar.col = "black")
    legend("bottomright", c("MC", "Shrunk MLE", "Standard MLE"), fill=c(rgb(1,1,1,0.5),rgb(0.7,0.7,0.7,1),rgb(0.35,0.35,0.35,1)),
           cex = 0.85, bg="white")#,bty="n"
dev.off()
#*******************************************************
# Test the difference of means on the first bin
#******************************************************
# One tailed at level 5%  
#alternative = "greater" is the alternative that x has a larger mean than y.

t.test(counts.shrunk[1,], counts.monte[1,], alternative = "less", var.equal = FALSE)
t.test(counts.ENF[1,], counts.monte[1,], alternative = "less", var.equal = FALSE)
    




#***********************************************************************************
# Answer to Reviewer 1 comment 2
# (Averaged) False Positives number with respect to the ratio p/n.
# number of data sets: rep*length(p)
#***********************************************************************************
rm(list=ls())

library(GeneNet)
library(stats4)
library(ggplot2)


set.seed(123)
source("H:/Semester7(15032017)/RUG/Shrinkage/shrinkagefunctions.R")
#setwd("H:/Semester7(15032017)/RUG/Shrinkage")

#***************************
# Initialize parameters
p <- c(50, 100, 150, 200 ) # num of genes
n <- c(20)  # num of samples
number <- 15 # Montecarlo iterations
rep <- 25 # repeat each simulation (for fixed p,n)
etaA <- 0.01  # proportion of  TP
alpha <- 0.05 # significance level
eta0 <- 1-etaA


## Initialize variables
FP <- matrix(Inf,3,length(p)) 
se.FP <- matrix(Inf,3,length(p))
incorrect <- matrix(Inf,3,rep) 
L <- matrix(Inf,2,length(p))
all.lambda <- matrix(Inf,1,rep) 
#***************************************


for (times in 1:length(p)) { 
  
  for (again in 1:rep){      
    
    sim.pcor<-ggm.simulate.pcor(p[times], etaA) #simulate partial corr
    temp<-sm2vec(sim.pcor)
    TP<- which(temp!=0) # true positives
    sim.data<- ggm.simulate.data( n , sim.pcor)
    GGM <- pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
    lambda<-attr(GGM, "lambda") 
    
    
    while (lambda == 1) {
      sim.data<- ggm.simulate.data( n , sim.pcor)
      temp <- sim.pcor[lower.tri(sim.pcor, diag = FALSE) ]
      TP <- which(temp!=0 )
      TN <- which(temp==0 )
      GGM <- pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
      lambda <- attr(GGM, "lambda") 
      print("TRY AGAIN")
     }       
    
     r<-sm2vec(GGM)
      
      
      ##### P-VALUES
      #######################################
      # MLE (shrunk density [this work])
      pval.shrunk <- p.shrunk(r,p[times] ,n,lambda)
      # MLE (standard/std density [Fisher 1921])
      #pval.standard <- p.shrunk(r, p, n, 0)
      
      ## GENENET methods (empirical null fitting + std density [Fisher 1921])
      ## WATCH OUT:  the test must be in node order. Sort by node 2 and by node 1
      ENF.test<-network.test.edges(GGM, fdr=TRUE, plot=FALSE, verbose =FALSE)
      ENF.test<-ENF.test[order(ENF.test$node1,ENF.test$node2),]
      
      
      # Monte Carlo
      p.monte <- p.montecarlo(r,number,p[times],n,lambda)
      
      ## false positives
      incorrect[1,again]<- sum(!which(ENF.test$pval <= alpha) %in% TP )#/length(ENF.test$pval)
      incorrect[2,again]<- sum(!which(p.monte <= alpha) %in% TP )#/length(p.monte)
      incorrect[3,again]<- sum(!which(pval.shrunk <= alpha) %in% TP )#/length(pval.shrunk)
      
      ## lambda
      all.lambda[1,again]<-lambda
      

  #fp
  FP[1, times] <- mean(incorrect[1,])
  se.FP[1,times] <- sd(incorrect[1,])/sqrt(rep)
  FP[2,times] <- mean(incorrect[2,])
  se.FP[2,times] <- sd(incorrect[2,])/sqrt(rep)
  FP[3,times] <- mean(incorrect[3,])
  se.FP[3,times] <- sd(incorrect[3,])/sqrt(rep)
  
  # shrinkage value
  L[1,times] <- mean(all.lambda)

  
}
}

#### Positives vs sample size
Method<-c(rep("ENF",length(p)),rep("MC",length(p)),rep("Shrunk MLE",length(p)))

tgc<-data.frame( Method, rep(p/n,3),c( FP[1,], FP[2,], FP[3,]),c(se.FP[1,],se.FP[2,],se.FP[3,]) )
colnames(tgc)<-c("Method","p2n","FP","se")

p2<- ggplot(tgc, aes(x=p2n, y= FP, colour=Method)) + 
  geom_errorbar(aes(ymin= FP-2*se, ymax= FP+2*se), width=0.2,colour="gray23") +
  geom_line(aes(linetype=Method),alpha = 1, size=1) + #, alpha = 0.7
  geom_point(aes(shape=Method), size=3,stroke = 0.2) +
  scale_color_manual(values=c("gray13","gray21","black"))  +
  scale_linetype_manual(values=c("dotted", "solid", "dashed"))+
  scale_shape_manual( values=c(0,6,20))+
  scale_fill_manual(values=c("black","white","black"))


win.metafile(filename=paste0(etaA,paste0(alpha,"FPvsp2n.emf")), width = 7, height = 7, pointsize = 20)
  p2+ theme_bw(base_size = 18)+ scale_x_continuous(breaks = tgc$p2n )+
    theme(panel.grid.minor = element_blank(),panel.grid.major  = element_blank())+
    xlab("p/n") + 
    ylab("False Positives")
dev.off()




#***********************************************************************************
# Answer to reviewer 3 comment 1
# q-q plot of p-values vs the uniform distribution
# number of data sets: length(n) 
#***********************************************************************************
rm(list=ls())

library(GeneNet)
library(stats4)
library(ggplot2)


set.seed(123)
source("H:/Semester7(15032017)/RUG/Shrinkage/shrinkagefunctions.R")


#***************************
# Initialize parameters
p <- 100# genes' number
etaA <- 0.00 #proportion true positive correlations 
n <- c(15) # sample size
alpha <- 0.05
number <- 15 # MC iterations


sim.pcor<-ggm.simulate.pcor(p, etaA)

    sim.data <- ggm.simulate.data( n , sim.pcor)
    temp <- sim.pcor[lower.tri(sim.pcor, diag = FALSE) ]
    TP <- which(temp!=0 )
    TN <- which(temp==0 )
    GGM <- pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
    lambda <- attr(GGM, "lambda") 
    
    while (lambda == 1) {
      sim.data <- ggm.simulate.data( n , sim.pcor)
      temp <- sim.pcor[lower.tri(sim.pcor, diag = FALSE) ]
      TP <- which(temp!=0 )
      TN <- which(temp==0 )
      GGM <- pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
      lambda <- attr(GGM, "lambda") 
      print("TRY AGAIN")
    } 
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
    #pval.standard<- ENF.test$pval
    
    # MC
    #p.monte<-p.montecarlo(r,number,p,n,lambda)

      #Generate data from Exp(1) distribution
      DATA <- sort(pval.shrunk) #rexp(40,1);
      DATA.ENF <- sort(ENF.test$pval)
      
      #QQ plot of data against Exp(1) distribution
      N         <- length(DATA)
      prob  <-  (1:N )/(N +1) 

      QUANTILES <- qunif(prob, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE) 
      PLOTDATA <- data.frame( Method = c( rep("ShrunkMLE",length(DATA)), rep("ENF",length(DATA))),
                                          pvals= c(DATA, DATA.ENF), Theoretical = c(QUANTILES, QUANTILES))

      QQPLOT<- ggplot(data = PLOTDATA, aes(x = Theoretical, y = pvals, color = Method)) +
        geom_point(size = 0.5, shape = 1 ) +
        geom_abline(intercept = 0, slope = 1, linetype = 'dashed', size=1) + 
        xlab('Theoretical quantiles U[0,1]') + 
        ylab('Empirical quantiles')+ theme_bw(base_size = 18)+
        scale_x_continuous(breaks = seq(0, 1, by = 0.05)) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.05))+
        scale_color_manual(values=c("gray", "gray21", "#56B4E9"))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = c(0.8,0.2), panel.grid.minor = element_blank(),panel.grid.major  = element_blank())+
        
        guides(colour = guide_legend(override.aes = list(size = 5, shape = 16)))

      
win.metafile(filename=paste0(n,paste0(p,"qqplot.emf")), width = 7, height = 7, pointsize = 20)
QQPLOT
dev.off()  

#***********************************************************************************
# Answer to Reviewer 3 comment 2
# Type-I error rate as function of decreasing sample size
# number of data sets: times*length(n) 
#***********************************************************************************
rm(list=ls())

library(GeneNet)
library(stats4)
library(ggplot2)


set.seed(321)
source("H:/Semester7(15032017)/RUG/Shrinkage/shrinkagefunctions.R")
setwd("H:/Semester7(15032017)/RUG/Shrinkage")

#***************************
# Initialize parameters
p <- 100 # genes' number
etaA <- 0.00 #proportion true positive correlations 
n <- c(1:15)*10 # sample size
alpha <- 0.05
number <- 15 # MC iterations
times <- 5

## Simulated 25 pcorr and data
FPR.shrunk <- matrix(Inf, length(n), times)
FPR.ENF <- matrix(Inf,  length(n), times)
FPR.monte <- matrix(Inf, length(n), times)

sim.pcor <- ggm.simulate.pcor(p, etaA)
i <- 1
t <- 1

for (i in 1:length(n)){
  
  for (t in 1:times){ 
    
    sim.data <- ggm.simulate.data( n[i] , sim.pcor)
    temp <- sim.pcor[lower.tri(sim.pcor, diag = FALSE) ]
    TP <- which(temp!=0 )
    TN <- which(temp==0 )
    GGM <- pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
    lambda <- attr(GGM, "lambda") 
      
    while (lambda == 1) {
      sim.data <- ggm.simulate.data( n[i] , sim.pcor)
      temp <- sim.pcor[lower.tri(sim.pcor, diag = FALSE) ]
      TP <- which(temp!=0 )
      TN <- which(temp==0 )
      GGM <- pcor.shrink(sim.data, verbose=FALSE) # OPTIMAL
      lambda <- attr(GGM, "lambda") 
      print("TRY AGAIN")
      } 
          
          r<-sm2vec(GGM)
          
          ## GENENET methods (empirical null fitting + std density)
          ## WATCH OUT:  the test must be in node order. Sort by node 2 and by node 1
          ENF.test<-network.test.edges(GGM, fdr=TRUE, plot=FALSE, verbose =FALSE)
          ENF.test<-ENF.test[order(ENF.test$node1,ENF.test$node2),]
          
          #######################################
          # MLE (shrunk density)
          pval.shrunk<-p.shrunk(r, p, n[i], lambda)
          # MLE (standard/std density)
          # pval.std<-p.shrunk(r,p,n,0)
          
          # ENF.test$pval
          
          # MC
          p.monte<-p.montecarlo(r, number, p, n[i], lambda)
          #********************************************************* 
          FPR.shrunk[i,t]<-sum(!which(pval.shrunk<=alpha) %in% TP )/length(pval.shrunk)
          FPR.ENF[i,t] <- sum(!which(ENF.test$pval<=alpha) %in% TP )/length(ENF.test$pval)
          FPR.monte[i, t]<-sum(!which(p.monte<=alpha) %in% TP )/length(p.monte)
      
          
      
  }

}

#***************************************
pl3 <- c(rep("Shrunk MLE", length(n) ),rep("ENF", length(n)), rep("MC", length(n)))
tgc3<-data.frame( pl3, c(rep(n,3)), 
                 c( rowMeans(FPR.shrunk, na.rm = FALSE), 
                    rowMeans(FPR.ENF, na.rm = FALSE),
                    rowMeans(FPR.monte, na.rm = FALSE)),
                 
                 c( apply(FPR.shrunk,1,sd),
                    apply(FPR.ENF,1,sd),
                    apply(FPR.monte,1,sd) )/sqrt(n) 
)
colnames(tgc3)<-c("Method","n","FPR","se")

ph <- ggplot(tgc3, aes(x = n, y = FPR, fill = Method)) + 
  geom_line(aes(linetype = Method), size= 1.25)+  ylim(0, max(tgc3$FPR+2*tgc3$se))+ 
  geom_point(aes(shape = Method), size=4 , stroke = 0.2)+
  scale_color_manual(values=c("gray13","gray21","black"))  +
  scale_linetype_manual(values=c("dotted", "solid", "dashed"))+
  scale_shape_manual( values=c(0,6,20))+
  scale_fill_manual(values=c("black","white","black")) 

win.metafile(filename=paste0(etaA,paste0(alpha,"type1vsN.emf")), width = 7, height = 7, pointsize = 20)
ph + geom_errorbar(aes(ymin= FPR-2*se, ymax= FPR+2*se, width=0.75 ))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())+
  theme_bw(base_size = 18)+ scale_x_continuous(breaks = n)+
  theme(panel.grid.minor = element_blank(), panel.grid.major  = element_blank(), legend.position = c(0.8, 0.2)) + #legend.justification=c(1,0), legend.position="right"
  geom_hline(yintercept = alpha, size=0.35)
dev.off()
#***************************************

save.image(paste(p,"type1_n20.R"))

