
findCommonVars <- function(varCol1, varCol2) {
    match <- intersect(varCol1, varCol2)
    rowPop1 <- which(varCol1 %in% match)
    rowPop2 <- which(varCol2 %in% match)
    print(length(match))
    print(length(rowPop1))
    print(length(rowPop2))
    return(list(rowPop1  = rowPop1, rowPop2 = rowPop2))
}

includeBiallelicOnly <- function(varCol) {
    dupVariant <- varCol[which(duplicated(varCol) == TRUE)]
    if (length(dupVariant) > 0) {
        removeDup <- which(varCol %in% dupVariant)
        return(removeDup)
    }
}

yMaxLim <- function(x, base)
{ 
  base * round(x/base) 
}


frqFileExtension <- "_dp6_anc_f_dbsnp_snpeff.daf.frq"
delta = 0.6
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
      pop1FrqFile <- read.delim2(paste(popList[pop1], "_phased_chr", chr, frqFileExtension, sep=""), skip=1, header=T)
      pop2FrqFile <- read.delim2(paste(popList[pop2], "_phased_chr", chr, frqFileExtension, sep=""), skip=1, header=T)
      varCol1 <- as.vector(pop1FrqFile[,2])
      varCol2 <- as.vector(pop2FrqFile[,2])
      duplicatedRow1 <- includeBiallelicOnly(varCol1)
      duplicatedRow2 <- includeBiallelicOnly(varCol2)
      if (length(duplicatedRow1) > 0) {
          pop1FrqFile <- pop1FrqFile[-duplicatedRow1,]
          varCol1 <- varCol1[-duplicatedRow1]
      }
      if (length(duplicatedRow2) > 0) {
         pop2FrqFile <- pop2FrqFile[-duplicatedRow2,]
         varCol2 <- varCol2[-duplicatedRow2]
      }
      if (setequal(varCol1, varCol2) == FALSE) {
          selectRows = findCommonVars(varCol1, varCol2)
          pop1FrqDf = pop1FrqFile[selectRows[[1]],]
          pop2FrqDf = pop2FrqFile[selectRows[[2]],]
      } else {
          pop1FrqDf = pop1FrqFile
          pop2FrqDf = pop2FrqFile
      }
      frqCol1 <- as.vector(pop1FrqDf[,5])
      frqCol1Split <- strsplit(frqCol1, split=":")
      frq1 <- as.numeric(sapply(frqCol1Split, "[[", 2))
      frqCol2 <- as.vector(pop2FrqDf[,5])
      frqCol2Split <- strsplit(frqCol2, split=":")
      frq2 <- as.numeric(sapply(frqCol2Split, "[[", 2))
      deltaVector <- c(deltaVector, frq2-frq1)
      
      ## Now we also create a file with the HDV information for that pairwise population
      HDVpos <- which(abs(deltaVector) >= 0.6)
      if (length(HDVpos) > 0) {
        HDVcounter <- HDVcounter+length(HDVpos)
        out <- rbind(out, cbind(pop1FrqFile[HDVpos,c(1,2,5)],pop2FrqFile[HDVpos,5]))
      }
      print(paste("end chr",chr))
    } #end chr

    if (HDVcounter > 0) { 
      write.table(out, paste("HDVList_",popList[pop1],"-",popList[pop2],"Delta60.txt",sep=""),quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
      ## Now we populate the matrix with pairwise pop HDVs    
      HDVmat[pop1,pop2] <- sum(HDV1,HDV2,na.rm = T)
    } else {
      ## Now we populate the matrix with pairwise pop HDVs    
      HDVmat[pop1,pop2] <- 0
    }
    ## Now we group the delta values in 0.05 frequency bins
    frqBin <- c( (-80:80)/(40*2) )
    binVector <- seq(-1,1,by=0.0500)
    deltaVectorGrouped <- deltaVector[which(deltaVector != 0)]
    numHitsBinTable <- table(deltaVectorGrouped)
    numHitsBin = as.vector(numHitsBinTable) #DELTA_vector #as.vector(t(DeltaCounts[r,]))
    frqList <- as.numeric(names(numHitsBinTable))
    ## We need to round the numbers so that we don't have to deal with the nightmare of equals and decimals in R
    frqList <- round(frqList, digits=2)
    binVector <- round(binVector, digits=2)
    v <- c()
    for (bin in 1:length(binVector)) {
      if (binVector[bin] < 0) {
        v <- c(v, sum(numHitsBin[which(frqList >= binVector[bin] & frqList < binVector[bin+1] ) ]) )
      }  
      if (binVector[bin] > 0) {
        v <- c(v, sum(numHitsBin[which(frqList > binVector[bin-1] & frqList <= binVector[bin] ) ]) )
      }
    }
    # Now we save v so that we are able to know in the future what are the exact number of variants for each ferquency bin
    write.table(rbind(binVector[-which(binVector==0)],v), paste(popList[pop1],"-",popList[pop2], "_delta60Distribution.txt",sep=""), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
    # Now we plot. We generate a barplot and two zoomed in version focusing on HDVs
    t <- barplot(v, plot=FALSE)
    pos <- c(seq(1,20, by=2),seq(22,length(t),by=2))
    lab = seq(-1,1,by=0.1)
    lab = lab[-which(lab==0)]
    
    pdf(file = paste("Figure1_",popList[pop1],popList[pop2],"_Delta60Distribution.pdf",sep=""),width = 16, height=7)
    layout(matrix(c(0,1,0,2,1,3,2,1,3,2,1,3,0,1,0), nrow=5, ncol=3, byrow = TRUE), widths=c(2,8,2), height=1, respect = FALSE)
    barplot(rev(v), col=c(rep("blue4",20),rep("gold",20)), main=paste(pop1,pop2,"Delta distribution"), cex.axis=1.4, cex.main=2 )
    axis(1, at=t[pos], labels = lab, cex.axis=1.4)
    
    HDV1 <- v[which(binVector <= -0.6)]
    HDV2 <- v[which(binVector >= 0.6)]
    
    if (sum(HDV1, HDV2, na.rm = T) > 0) {
      if ( max(HDV1,HDV2,na.rm = T) < 10) {
        yax_max <- 10
      } else {
        yax_max <- yLimMax(max(HDV1,HDV2,na.rm = T),10^(nchar(max(HDV1,HDV2,na.rm = T))-1) )+ (10^(nchar(max(HDV1,HDV2,na.rm = T))-1)/2)
      }
      par(mar=c(5.1,1.1,0,1))
      barplot(HDV1, col=rep("blue4",9), names.arg = seq(-1,-0.60,by=0.05), axes=FALSE, cex.names=1.5, ylim=c(0,yax_max) )
      axis(4, at=seq(from=0, to=yax_max,length.out = 11), labels = , cex.axis=1.3)
      yText = yax_max - (10^(nchar(max(HDV1,HDV2, na.rm = T))-1)/2)
      if (sum(rev(v[which(bin >= 0.6)-1])) != 0 ) {
        text(2.75,yText,paste(pop1,"=", formatC( sum(HDV1, na.rm = T), format="d",big.mark=","  ) ), cex=1.3)
      }
      par(mar=c(5.1,0,0,1.1))
      barplot(HDV2, col=rep("gold",9), names.arg = seq(0.60,1,by=0.05), cex.names=1.5, axes=FALSE, ylim=c(0,yax_max) )
      axis(2, at=seq(from=0, to=yax_max,length.out = 11), labels = , cex.axis=1.3)
      #yText = max(rev(v[which(bin <= -0.6)]))/1.2
      if (sum(rev(v[which(bin <= -0.6)])) != 0 ) {
        text(8,yText,paste(pop2,"=",formatC( sum(HDV2,na.rm = T),format="d",big.mark = "," ) ), cex=1.3)
      }
    } 
    dev.off()

  } # end pop1
  
} # end pop2
      

      
    
    
    
    
    
