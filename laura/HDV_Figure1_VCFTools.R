
findCommonVars <- function(varPop1, varPop2) {
    match <- intersect(varPop1, varPop2)
    rowPop1 <- which(varPop1 %in% match)
    rowPop2 <- which(varPop2 %in% match)
    print(length(match))
    print(length(rowPop1))
    print(length(rowPop2))
    return(list(rowPop1  = rowPop1, rowPop2 = rowPop2))
}

includeBiallelicOnly <- function(varPop) {
    dupVariant <- varPop[which(duplicated(varPop) == TRUE)]
    if (length(dupVariant) > 0) {
        biallelic <- varPop[-which(varPop %in% dupVariant)]
        return(biallelic)
    }
}

yMaxLim <- function(x, base)
{ 
  base * round(x/base) 
}

extractFrequencies <- function(frqFileColumn)
{
    colSplit <- strsplit(frqFileColumn, split=":")
    frq <- as.numeric(sapply(colSplit, "[[", 2))
    return(frq)
}


frqFileExtension <- "_dp6_anc_f_dbsnp_snpeff.daf.frq"
delta = 0.4
popList <- c('Mali', 'Zambia') # Here you need to include the names of all your populations.
HDVmat <- matrix(0, nrow=length(popList), ncol=length(popList),dimnames = list(popList,popList))
for (pop1 in 1:(length(popList)-1)) {
  for (pop2 in (pop1+1):length(popList)) {
    print(paste("starting", popList[pop1], "-", popList[pop2]))
    deltaVector <- c()
    ### We first calculate the delta values (differences in derived allele frequencies)
    HDVcounter=0 ## This will tells us whether we have HDVs at all, to write an output log file
    out <- c()
    for (chr in 22:22) {
      pop1FrqFile <- read.delim2(paste(popList[pop1], "_phased_chr", chr, frqFileExtension, sep=""), skip=1, header=T, stringsAsFactors=FALSE)
      pop2FrqFile <- read.delim2(paste(popList[pop2], "_phased_chr", chr, frqFileExtension, sep=""), skip=1, header=T, stringsAsFactors=FALSE)
      varPop1 <- as.vector(pop1FrqFile[,2])
      varPop2 <- as.vector(pop2FrqFile[,2])
      biallelicPop1 <- includeBiallelicOnly(varPop1)
      biallelicPop2 <- includeBiallelicOnly(varPop2)
      pop1FrqFileBiallelic <- pop1FrqFile[which(varPop1 %in% biallelicPop1),]
      pop2FrqFileBiallelic <- pop2FrqFile[which(varPop2 %in% biallelicPop2),]
      if (setequal(biallelicPop1, biallelicPop2) == FALSE) {
          selectRows = findCommonVars(biallelicPop1, biallelicPop2)
          pop1FrqDf = pop1FrqFileBiallelic[selectRows[[1]],]
          pop2FrqDf = pop2FrqFileBiallelic[selectRows[[2]],]
      } else {
          pop1FrqDf = pop1FrqFileBiallelic
          pop2FrqDf = pop2FrqFileBiallelic
      }
      frq1 <- extractFrequencies(as.vector(pop1FrqDf[,5]))
      frq2 <- extractFrequencies(as.vector(pop2FrqDf[,5]))
      deltaVector <- c(deltaVector, frq2-frq1)
      
      ## Now we also create a file with the HDV information for that pairwise population
      HDVpos <- which(abs(deltaVector) >= 0.4)
      if (length(HDVpos) > 0) {
        HDVcounter <- HDVcounter+length(HDVpos)
        out <- rbind(out, cbind(pop1FrqDf[HDVpos,c(1,2,5)],pop2FrqDf[HDVpos,5]))
      }
      print(paste("end chr",chr))
    } #end chr

    if (HDVcounter > 0) { 
      write.table(out, paste("HDVList_",popList[pop1],"-",popList[pop2],"Delta60.txt",sep=""),quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
      ## Now we populate the matrix with pairwise pop HDVs    
      HDVmat[pop1,pop2] <- HDVcounter
    } else {
      ## Now we populate the matrix with pairwise pop HDVs    
      HDVmat[pop1,pop2] <- 0
    }
    ## Now we group the delta values in 0.05 frequency bins
    binVector <- seq(-1,1,by=0.1)
    deltaVectorNoZero <- deltaVector[which(deltaVector != 0)]
    frqList <- round(deltaVectorNoZero, digits=1)
    # Now we save v so that we are able to know in the future what are the exact number of variants for each ferquency bin
    write.table(rbind(binVector[-which(binVector==0)],v), paste(popList[pop1],"-",popList[pop2], "_delta60Distribution.txt",sep=""), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
    
    pdf(file = paste("Figure1_",popList[pop1],popList[pop2],"_Delta", delta, "Distribution.pdf",sep=""),width = 16, height=7)
    layout(matrix(c(0,1,0,2,1,3,2,1,3,2,1,3,0,1,0), nrow=5, ncol=3, byrow = TRUE), widths=c(2,8,2), height=1, respect = FALSE)
    hist(frqList, col="firebrick", breaks=binVector, xlab = paste("Delta =", delta, sep=" "), ylab="Number of variants", main=paste(popList[pop2],"vs.",popList[pop1], sep=" "), cex.axis=1.5 ) # Negative values are variants more frequent in Zambia.
    HDV1 <- deltaVector[which(deltaVector <= -0.4)]
    HDV2 <- deltaVector[which(deltaVector >= 0.4)]
    if (HDVcounter > 0) {
      if ( HDVcounter < 10) {
        yax_max <- 10
      } else {
        yax_max <- max(length(HDV1), length(HDV2))
      }
      par(mar=c(5.1,2.1,2,1))
      hist(round(HDV1, digits = 2), col = rep("blue4",9), breaks = seq(-1, -delta, by=0.05), ylim=c(0, yax_max), xlab = "", ylab="", main=paste(popList[pop2], "=", length(HDV1), sep = " ") )
      par(mar=c(5.1,0,2,1.1))
      hist(round(HDV2, digits=1), col=rep("gold",9), breaks = seq(delta,1,by=0.05), ylim=c(0, yax_max), xlab="", ylab="", main=paste(popList[pop1], "=", length(HDV2), sep="") )
    } 
    dev.off()

  } # end pop1
  
} # end pop2
      

      
    
    
    
    
    
