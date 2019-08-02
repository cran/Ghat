#' Evaluation of Linkage Disequilibrium Decay
#'
#' ld_decay: R function for calculating the effective number of independent markers
#'
#' @param gen Matrix of genotype data. Individuals in rows, genotypes (0, 1, 2) in columns.
#' @param map Dataframe inculding the name for each marker with a corresponding chromosome number and physical position.
#' @param max_win_snp The maximum number of markers in each window. Sets the maximum number of markers allowed per window within a chromosome before estimating the LD. Default is 2000.
#' @param max.chr Chromosomes above this number will be excluded from the analysis.
#' @param cores Numer of cores for using parallelized calculation, Default is 1 for windows machine.
#' @param max_r2 the threshold of r^2 to calculate the effective number of independent markers.
#'
#' @return cor: Correlation matrix
#' @return ch_eff_nmark: The Number of independent marker per chromosome
#' @return eff_nmark: The effective number of independent markers
#'
#' @export
#'
#' @examples
#' \donttest{
#' library("parallel")
#' gen         <- Maize_wqs[[4]]
#' map         <- Maize_wqs[[3]]
#' Res_ld <- ld_decay (gen=gen, map=map, max_win_snp=2000,
#'                     max.chr=10, cores=1, max_r2=0.03)
#' }                     
#'                     

ld_decay <- function( gen= gen, map=map, max_win_snp=2000,
                      max.chr=max.chr, cores=1, max_r2=max_r2){



  map         <- map[order(map[,2],map[,3])==1:nrow(map),]
  remMap      <- which(map$Chromosome == 0)
  remGeno     <- remMap+1
  gen         <-gen[,-remGeno]
  map         <-map[-remMap,]

  cor1        <- function(x1,x2){
    x1        <- as.matrix(x1)
    x2        <- as.matrix(x2)
    sx1       <- t(x1)-colMeans(x1)
    sx2       <- t(x2)-colMeans(x2)
    cc        <- (rowSums(sx1*sx2))/(sqrt(rowSums(sx1^2)*rowSums(sx2^2)))
    return(cc)
  }


  alld1       <- matrix(0,nrow=max(map$Chromosome),ncol=max_win_snp)
  amap        <- map
  a           <- gen[,2:ncol(gen)]
  rownames(a) <- gen$X

  for(rr in 1:max.chr){
    message(paste("Chromosome",rr))
    af        <- a[,amap[,2]==rr,drop=FALSE]
    if(ncol(af)==0) next
    message(dim(af))
    ind       <- 1:ncol(af)

    fun       <- function(i){
      if(i<=ncol(af)){ld1 <- mean(cor1(af[,head(ind,-i)],af[,tail(ind,-i)])^2,na.rm=TRUE)}
      else{ld1<- 0}
      return(ld1)
    }
    ld1       <- unlist(mclapply(1:2000, FUN=fun,  mc.preschedule = TRUE, mc.cores = cores))
    alld1[rr,]<- ld1
  }


  laeng       <- unlist(lapply(by(amap[,3],amap[,2],diff),FUN=median))

  f           <- function(x){return(min(which(x<max_r2)))}
  dd          <- apply(alld1,1,f)
  eff_nmark   <- sum(table(map$Chromosome)/dd)

  ##

  ltt         <- rep(1:5,each=8)   ## to set the line types

  sc          <- rbind(c(0,0.8,0,1),c(0.75,1,0,1),c(0.4,0.75,0.5,1))
  split.screen(sc)

  screen(1)
  plot(alld1[1,],ylab="r^2", xlab="Distance in #SNPs", type="l",
       xlim=c(0,500),ylim=c(0,0.5),col="white",bty="n",
       main=paste( "Number of SNPs in LD at Threshold of <", max_r2), adj = 0)
  for(i in 1:nrow(alld1)){
    lines(alld1[i,], lwd=2, col=i,lty=ltt[i])
  }
  abline(h=0.03,col="red",lwd=2)

  screen(2)
  opar <- par(mai=c(0,0,0,0))
  on.exit(par(opar))
  plot(1:1000,1:1000, col="white",axes=FALSE,xlab="", ylab="",bty="n")

  legend(1,900, legend=c("Chr: #SNPs",
                         paste(c(1:10),dd,sep="   :   ")),col=0:nrow(alld1),
         lty=c(0,ltt[1:nrow(alld1)]),lwd=2,bty="n",cex=0.9)

  close.screen(all.screens=TRUE)
###
return (list(cor=alld1, ch_eff_nmark=dd, eff_nmark=eff_nmark))
###
  }
 


"ld_decay"