# <img src="https://github.com/Medhat86/Ghat/blob/master/GOE_Logo_Quer_IPC_Farbe_RGB.png" />


# Ghat 

Functions are provided for quantifying evolution and selection on complex traits.
The package implements effective handling and analysis algorithms scaled for
genome-wide data and calculates a composite statistic, denoted Ghat, which is used
to test for selection on a trait. The package provides a number of simple examples
for handling and analysing the genome data and visualising the output and results.
Beissinger et al., (2018) <doi:10.1534/genetics.118.300857>


## Authors

### Medhat Mahmoud*, ‡,§,1, Mila Tost*,1, ‡, Ngoc-Thuy Ha †, Henner Simianer †, ‡, and Timothy Beissinger*, ‡,1 

#### *Department of Crop Science, University of Goettingen, Goettingen, Germany 

#### †Department of Animal Sciences, University of Goettingen, Goettingen, Germany 

#### ‡Center for Integrated Breeding Research, University of Goettingen, Goettingen, Germany

#### §Creator
#### 1Maintainer


## Installation

### Latest development build

To install the latest development snapshot (see latest changes below), type the following commands into the R console:

```r
library(devtools)
devtools::install_github("Medhat86/Ghat")
```

### Official, stable release

To install the latest stable release from CRAN, type the following command into the R console:

```r
install.packages("Ghat")
```

## Documentation and examples


### Example-1 Both SNP effects and change in allele frequency are known


```r
maize		<- Maize_wqs[[1]]
result.adf	<- Ghat(effects =maize[,1], change=maize[,2], method="scale",
                     perms=1000, plot="Ghat", num_eff=54.74819)
mtext(paste("WQS ADF test for selection, pval = ", round(result.adf$p.val,4)))
message (c(result.adf$Ghat , result.adf$Cor , result.adf$p.va))
```

### Example-2 Both SNP effects and change in allele frequency are known
#### Step 1: Run rrBLUP and estimating allels effects

```r
library(Ghat)
library(parallel)
library(rrBLUP)
phe                 <- Maize_wqs[[2]]
map                 <- Maize_wqs[[3]]
gen                 <- Maize_wqs[[4]]
phe                 <-phe[which(is.na(phe[,2])==FALSE),]
gen                 <-gen[which(is.na(phe[,2])==FALSE),]
result              <- mixed.solve(phe[,2],
                                   Z= as.matrix(gen[,2:ncol(gen)]),
                                   X= model.matrix(phe[,2]~phe[,3]),
                                   K=NULL, SE=FALSE, return.Hinv=FALSE,
                                   method="ML")
``` 
                       
#### Step 2: Is to calculate the allele frequency at Cycle 1 and 3 

```r  
CycleIndicator      <- as.numeric(unlist(strsplit(gen$X,
                       split="_C")) [seq(2,2*nrow(gen),2)])
Cycle1              <- gen[which(CycleIndicator == 1),]
Cycle3              <- gen[which(CycleIndicator == 3),]
CycleList           <- list(Cycle1,Cycle3)
frequencies         <- matrix(nrow=ncol(gen)-1,ncol=2)
for(i in 1:2){
  frequencies[,i]   <- colMeans(CycleList[[i]][,-1],na.rm=TRUE)/2
}
frequencies         <- as.data.frame(frequencies)
names(frequencies)  <- c("Cycle1","Cycle3")
change<-frequencies$Cycle3-frequencies$Cycle1
```

#### Step 3: Calculate LD Decay 

```r
ld                  <- ld_decay (gen=gen, map=map,
                                 max_win_snp=2000, max.chr=10,
                                 cores=1, max_r2=0.03)
```

#### Step 4: Calculate Ghat

```r
Ghat.adf    <- Ghat(effects=result$u, change=change, method = "scale",
                    perms=1000,plot="Ghat", num_eff = 54.74819)

message (paste("Ghat=" , Ghat.adf$Ghat,
            "Cor="  , Ghat.adf$Cor ,
            "P-val=", Ghat.adf$p.va, sep = " "))
```

Please visit [https://github.com/Medhat86/Ghat](https://CRAN.R-project.org/package=Ghat) for documentation and examples.

## Citation

- In case you want / have to cite my package, please use `citation('Ghat')` (https://CRAN.R-project.org/package=Ghat) for citation information. 

- Since core functionality of package depends on the [Ghat-Method](https://cran.r-project.org/package=ggplot2), consider citing this refernce as well: "Tim Beissinger, Jochen Kruppa, David Cavero, Ngoc-Thuy Ha, Malena Erbe, Henner Simianer, A Simple Test Identifies Selection on Complex Traits, Genetics, Volume 209, Issue 1, 1 May 2018, Pages 321–333, https://doi.org/10.1534/genetics.118.300857"
