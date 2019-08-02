#' Quantifying evolution and selection on complex traits 
#'
#' G-hat: R function to estimate G-hat from allele frequency and effect size data.
#'
#' @param effects Vector of allele effects.
#' @param change Vector of changes in allele frequency (could be positive, negative or zero).
#' @param method "vanilla" (assumes complete linkage equilibrium between markers), "trim" (excludes markers to approximate linkage equilibrium some of the extreme values) or "scale" (scales results to reflect underlying levels of linkage LD)
#' @param perms Number of permutations to run.
#' @param plot "Ghat", "Cor", or "Both", Should a plot of the Ghat or correlation test be returned?
#' @param blockSize How large should blocks for trimming be? Only required if method = "trim".
#' @param num_eff The effective number of independent markers, to be used only in conjunction with the “scale” method, above (see “ld_decay” function or use help (?ld_decay).
#'
#' @return Ghat Ghat-value
#' @return Cor  Correlation between alleles frequencies and their effects
#' @return p.val two-sided P-value of Evidence of selection
#' @return plot relationship between estimated allelic effects at individual SNPs and the change in allele frequency over generations
#'
#' @importFrom graphics abline close.screen hist legend lines par plot screen split.screen
#' @importFrom stats cor dnorm lm median pnorm sd
#' @importFrom utils head tail
#' @importFrom parallel mclapply
#' @importFrom rrBLUP mixed.solve
#'
#' @export
#'
#' @examples
#' #Example-1 Both SNP effects and change in allele frequency are known
#' maize		<- Maize_wqs[[1]]
#' result.adf	<- Ghat(effects =maize[,1], change=maize[,2], method="scale",
#'                      perms=1000, plot="Ghat", num_eff=54.74819)
#' mtext(paste("WQS ADF test for selection, pval = ", round(result.adf$p.val,4)))
#' message (c(result.adf$Ghat , result.adf$Cor , result.adf$p.va))
#'
#'
#'\donttest{
#' #Example-2 Both SNP effects and change in allele frequency are known
#' ##################################################################
#' ## step 1: #run rrBLUP and estimating allels effects            ##
#' ##################################################################
#' 
#' library(Ghat)
#' library(parallel)
#' library(rrBLUP)
#' phe                 <- Maize_wqs[[2]]
#' map                 <- Maize_wqs[[3]]
#' gen                 <- Maize_wqs[[4]]
#' phe                 <-phe[which(is.na(phe[,2])==FALSE),]
#' gen                 <-gen[which(is.na(phe[,2])==FALSE),]
#' result              <- mixed.solve(phe[,2],
#'                                    Z= as.matrix(gen[,2:ncol(gen)]),
#'                                    X= model.matrix(phe[,2]~phe[,3]),
#'                                    K=NULL, SE=FALSE, return.Hinv=FALSE,
#'                                    method="ML")
#'                                    
#' ##################################################################
#' ## step 2: is to calculate the allele frequency at Cycle 1 and 3##
#' ##################################################################
#' CycleIndicator      <- as.numeric(unlist(strsplit(gen$X,
#'                        split="_C")) [seq(2,2*nrow(gen),2)])
#' Cycle1              <- gen[which(CycleIndicator == 1),]
#' Cycle3              <- gen[which(CycleIndicator == 3),]
#' CycleList           <- list(Cycle1,Cycle3)
#' frequencies         <- matrix(nrow=ncol(gen)-1,ncol=2)
#' for(i in 1:2){
#'   frequencies[,i]   <- colMeans(CycleList[[i]][,-1],na.rm=TRUE)/2
#' }
#' frequencies         <- as.data.frame(frequencies)
#' names(frequencies)  <- c("Cycle1","Cycle3")
#' change<-frequencies$Cycle3-frequencies$Cycle1
#' 
#' ################################################################
#' ## step 3: Calculate LD Decay                                   ##
#' ################################################################
#' ld                  <- ld_decay (gen=gen, map=map,
#'                                  max_win_snp=2000, max.chr=10,
#'                                  cores=1, max_r2=0.03)
#'
#' ##################################################################
#' ## step 4: Calculate Ghat                                       ##
#' ##################################################################
#' Ghat.adf    <- Ghat(effects=result$u, change=change, method = "scale",
#'                     perms=1000,plot="Ghat", num_eff = 54.74819)
#'
#' message (paste("Ghat=" , Ghat.adf$Ghat,
#'             "Cor="  , Ghat.adf$Cor ,
#'             "P-val=", Ghat.adf$p.va, sep = " "))
#'}
#'

Ghat <- function( effects = effects, change=change, method = "scale",perms=1000,
                  plot="Both",blockSize=1000, num_eff = NULL){
  sd_eff    <- NA
  if(method == "vanilla"){
    effects.mat <- as.matrix(effects)
    Ghat <- sum(change*effects,na.rm=TRUE)
    Cor <- cor(change,effects,use="pairwise.complete.obs")
    Ghat_perm <- c()
    Cor_perm <- c()
    message("Doing Permutation Test")
    for(test in 1:perms){
      Ghat_perm[test] <- sum(change*effects[sample(length(effects),replace=F)],na.rm=TRUE)
      Cor_perm[test] <- cor(change,effects[sample(length(effects),replace=F)],use="pairwise.complete.obs")
    }
    if(plot=="Ghat"){
      hist(Ghat_perm,xlim=c(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm))),col="grey",main="Ghat Test for selection",xlab="Null Distribution",freq=F)
      abline(v=Ghat,lwd=3,col="darkblue")
    }
    if(plot=="Cor"){
      hist(Cor_perm,xlim=c(min(c(Cor,Cor_perm)),max(c(Cor,Cor_perm))),col="grey",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    if(plot=="Both"){
      opar <- par(mfrow=c(1,2))
      on.exit(par(opar))
      hist(Ghat_perm,xlim=c(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm))),col="grey",main="Ghat Test for selection",xlab="Null Distribution")
      abline(v=Ghat,lwd=3,col="darkblue")
      hist(Cor_perm,xlim=c(min(c(Cor,Cor_perm)),max(c(Cor,Cor_perm))),col="grey",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    if(Ghat > 0)  p.val <- length(which(Ghat_perm >= Ghat | Ghat_perm <= -1*Ghat))/length(Ghat_perm)
    if(Ghat < 0)  p.val <- length(which(Ghat_perm <= Ghat | Ghat_perm >= -1*Ghat))/length(Ghat_perm)
  }
  if(method == "trim"){
    start <- sample(length(change)/blockSize,1)
    index<-seq(start,length(change),by=blockSize)
    change<-change[index]
    effects.mat <- as.matrix(effects)
    Ghat <- sum(change*effects,na.rm=TRUE)
    Cor <- cor(change,effects,use="pairwise.complete.obs")
    Ghat_perm <- c()
    Cor_perm <- c()
    message("Doing Permutation Test")
    for(test in 1:perms){
      Ghat_perm[test] <- sum(change*effects[sample(length(effects),replace=F)],na.rm=TRUE)
      Cor_perm[test] <- cor(change,effects[sample(length(effects),replace=F)],use="pairwise.complete.obs")
    }
    if(plot=="Ghat"){
      hist(Ghat_perm,xlim=c(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm))),col="grey",main="Ghat Test for selection",xlab="Null Distribution")
      abline(v=Ghat,lwd=3,col="darkblue")
    }
    if(plot=="Cor"){
      hist(Cor_perm,xlim=c(min(c(Cor,Cor_perm)),max(c(Cor,Cor_perm))),col="grey",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    if(plot=="Both"){
      opar <- par(mfrow=c(1,2))
      on.exit(par(opar))
      hist(Ghat_perm,xlim=c(min(c(Ghat,Ghat_perm)),max(c(Ghat,Ghat_perm))),col="grey",main="Ghat Test for selection",xlab="Null Distribution")
      abline(v=Ghat,lwd=3,col="darkblue")
      hist(Cor_perm,xlim=c(min(c(Cor,Cor_perm)),max(c(Cor,Cor_perm))),col="grey",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    if(Ghat > 0)  p.val <- length(which(Ghat_perm >= Ghat | Ghat_perm <= -1*Ghat))/length(Ghat_perm)
    if(Ghat < 0)  p.val <- length(which(Ghat_perm <= Ghat | Ghat_perm >= -1*Ghat))/length(Ghat_perm)
  }
  if(method == "scale"){
    effects.mat <- as.matrix(effects)
    Ghat <- sum(change*effects,na.rm=TRUE)
    Cor <- cor(change,effects,use="pairwise.complete.obs")
    Ghat_perm <- c()
    Cor_perm <- c()
    message("Doing Permutation Test")
    for(test in 1:perms){
      Ghat_perm[test] <- sum(change*effects[sample(length(effects),replace=F)],na.rm=TRUE)
      Cor_perm[test] <- cor(change,effects[sample(length(effects),replace=F)],use="pairwise.complete.obs")
    }
    scale = sqrt(length(effects))/sqrt(num_eff)
    sd_eff = sd(Ghat_perm)*scale
    sd_eff_cor = sd(Cor_perm)*scale
    if(plot=="Ghat"){
      left <- min(mean(Ghat_perm)-4*sd_eff,Ghat)
      right <- max(mean(Ghat_perm)+4*sd_eff,Ghat)
      plot(seq(left,right,length.out=1000),dnorm(seq(left,right,length.out=1000),mean=mean(Ghat_perm),sd=sd_eff),col="red",lwd=3,type="l",ylab="Density",xlab="Null Distribution",main="Scaled Ghat test for selection")
      abline(v=Ghat,lwd=3,col="darkblue")
    }

    if(plot=="Cor"){
      left <- min(mean(Cor_perm)-4*sd_eff_cor,Ghat)
      right <- max(mean(Cor_perm)+4*sd_eff_cor,Ghat)
      plot(seq(left,right,length.out=1000),dnorm(seq(left,right,length.out=1000),mean=mean(Cor_perm),sd=sd_eff_cor),col="red",lwd=3,type="l",ylab="Density",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }

    if(plot=="Both"){
      opar <- par(mfrow=c(1,2))
      on.exit(par(opar))
      left <- min(mean(Ghat_perm)-4*sd_eff,Ghat)
      right <- max(mean(Ghat_perm)+4*sd_eff,Ghat)
      plot(seq(left,right,length.out=1000),dnorm(seq(left,right,length.out=1000),mean=mean(Ghat_perm),sd=sd_eff),col="red",lwd=3,type="l",ylab="Density",main="Ghat Test for selection",xlab="Null Distribution")
      abline(v=Ghat,lwd=3,col="darkblue")
      left <- min(mean(Cor_perm)-4*sd_eff_cor,Cor)
      right <- max(mean(Cor_perm)+4*sd_eff_cor,Cor)
      plot(seq(left,right,length.out=1000),dnorm(seq(left,right,length.out=1000),mean=mean(Cor_perm),sd=sd_eff_cor),col="red",lwd=3,type="l",ylab="Density",main="Cor Test for selection",xlab="Null Distribution")
      abline(v=Cor,lwd=3,col="darkblue")
    }
    if(Ghat > mean(Ghat_perm))  p.val <- pnorm(Ghat,mean=mean(Ghat_perm),sd=sd_eff,lower.tail=F)*2
    if(Ghat < mean(Ghat_perm))  p.val <- pnorm(Ghat,mean=mean(Ghat_perm),sd=sd_eff,lower.tail=T)*2
  }
  return(list(p.val=p.val, Ghat=Ghat, Ghat_perm=Ghat_perm, sd_eff=sd_eff, Cor=Cor, Cor_perm=Cor_perm, effects=effects))
}
"Ghat"