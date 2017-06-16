## Generate SFS

setwd("~/HGDPgenomes/")
pop <- c("San", "Mbuti", "Mozabite", "Pathan", "Cambodian", "Yakut", "Maya")
col <-  c(colors()[c(466)], "blue", "turquoise", "green4", "yellow", "orange", "red")

for (popID in 1:length(pop)) {
    allCounts <- c()
    for (chr in 1:22) {
        frqFName = "NameofFrqFile"
        frqFile <- read.delim2(frqFName, skip=1, stringsAsFactors = FALSE, header=FALSE)
        derFrqCol <- as.vector(frqFile[,6])
        derFrq <- as.numeric(sapply(strsplit(derFrqCol, ":"), "[[", 2))    
        allCounts <- c(allCounts, derFrq) ## SFS
    } ## end chr loop
    hist(round(as.numeric(derFrq), digits=2), xlab="Frequency", ylab="Number of variants", main=paste("SFS", pop, sep=" "), col="red")

} ## end loop pop
