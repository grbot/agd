library(data.table)
library(karyoploteR)
library(biomaRt)
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
POP="BOT"
CHR=19
DFIMPACT = fread(sprintf("%s/%s_%d_Effect_Impact.csv",POP,POP,CHR), sep = "\t", header = TRUE)
DFSDT = fread(sprintf("%s/%s_%d_Singletons_Doubletons_Tripletons.csv",POP,POP,CHR), sep = "\t", header = TRUE)
DFFREQ = fread(sprintf("%s/%s_Eagle.baylor.%d_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss_dp6_ancOnly.daf.frq",POP,POP,CHR), sep = "\t",skip=1, header = FALSE)
DFFREQ = DFFREQ[,c("V1","V2","V6")]
colnames(DFFREQ)<-c("CHR","START","FREQ")
DFFREQ$END = DFFREQ$START
colnames(DFSDT)<-c("CHR","POS","SDT","TotalSingletons","PropSingletons","TotalDoubletons","PropDoubletons","TotalTripletons","PropTripletons","TotalVariants")
DFK = as.data.frame(DFIMPACT$POS)
DFK$CHR = DFSDT$CHR
DFK$HIGH = apply(DFIMPACT, 1, function(x) sum(x=="HIGH"))
DFK$TYPE = DFSDT$SDT
DFK$START = DFIMPACT$POS
DFK$END = DFK$START+1
DFK <-DFK[DFK$HIGH>0,]
colnames(DFK)<-c("END","CHR","HIGH","TYPE","START")
DFFREQS=DFFREQ[DFFREQ$START %in% DFK$START,]
DFK$FREQ=DFFREQS$FREQ
freqs=as.array(unlist(strsplit(DFFREQS$FREQ,":")))
freqs= as.double(freqs[(1:length(freqs))%%2==0])
DFK$FREQ=freqs
maxfreq = max(DFK$FREQ)
rm(DFIMPACT)
rm(DFSDT)
rm(DFFREQ)
rm(DFFREQS)
rm(freqs)
DFK = DFK[,c("CHR","START","END","HIGH","TYPE","FREQ")]
chr.regions = paste(CHR,DFK$START,sep=':')
chr.regions = paste(chr.regions,DFK$START,sep=":")

all.results=data.frame() 
#filterlist = list(chr.regions,"protein_coding") 
for (i in 1:length(chr.regions)){
  print(chr.regions[i])
  filterlist=list(chr.regions[i],"protein_coding")

  results=getBM(attributes = c("hgnc_symbol","entrezgene", "chromosome_name", "start_position", "end_position"), 
               filters = c("chromosomal_region","biotype"), 
               values = filterlist, 
               mart = ensembl)
  all.results=rbind(all.results,results)
}

write.csv(all.results, file="all_results.csv",row.names=F)
all.results=read.csv(file="all_results.csv",row.names=F,header=TRUE)
all.results=all.results[!is.na(all.results$entrezgene),]
all.results=all.results[!is.na(all.results$hgnc_symbol),]
all.results=all.results[!duplicated(all.results$entrezgene),]
all.results$hgnc_symbol=gsub('+[[:digit:]]+$', '', all.results$hgnc_symbol)
results=all.results[!duplicated(all.results$hgnc_symbol),]
results=results[order(results$start_position),]
if(dim(results)[1]>0){
	genes <- toGRanges(results[,c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol')])
	seqlevelsStyle(genes) <- "UCSC"
}
DFK$DBFREQ = 100*DFK$FREQ
DFK=DFK[order(DFK$START),]
high.rare.variants.data <- toGRanges(DFK)
names(mcols(high.rare.variants.data)) <- c("HIGH", "TYPE","FREQ","DBFREQ")
seqlevelsStyle(high.rare.variants.data) <- "UCSC"
colortype<-function(type){
	if(type=="S"){
		"red"
	}else{
		if(type=="D"){
		"blue"
		}else{
		"green"
		}
	}
}
cextype<-function(type){
	if(type=="S"){
		1
	}else{
		if(type=="D"){
		1.5
		}else{
		2
		}
	}
}
setEPS()
postscript(sprintf("POP_%s_CHR_%d_HIGH_EFFECTS_FREQ.eps",POP,CHR),horizontal=FALSE, paper="special",height=10,width=10)
kp <- plotKaryotype(plot.type=2,genome = "hg19", chromosomes=paste("chr",CHR,sep=''))
kpAddBaseNumbers(kp)
kpAddCytobandLabels(kp,srt=90, col="blue")

#kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, line.color = "green",r0=0.55, r1=0.75)
kpDataBackground(kp, data.panel = 1, r0=0, r1=1)
kpPoints(kp, data=high.rare.variants.data, y=high.rare.variants.data$HIGH,col=unlist(lapply(high.rare.variants.data$TYPE,colortype)),data.panel = 1,r0=0.1, r1=1,cex=DFK$DBFREQ/max(DFK$DBFREQ))
kpArrows(kp, data=high.rare.variants.data, x0=high.rare.variants.data$START, x1=high.rare.variants.data$START, y0=0, y1=1,cex=DFK$DBFREQ/max(DFK$DBFREQ),col=unlist(lapply(high.rare.variants.data$TYPE,colortype)),r0=0.1, r1=1,data.panel=1,length=0.2)
kpText(kp, data=high.rare.variants.data, x=high.rare.variants.data$START, y=0, labels=high.rare.variants.data$TYPE,col=unlist(lapply(high.rare.variants.data$TYPE,colortype)),cex=0.7,r0=0, r1=0.1,data.panel=1)
kpText(kp, data=high.rare.variants.data, x=0, y=0.5, labels="Variants HEI", col="blue", srt=90,cex=1,data.panel=1)
#kpPlotMarkers(kp, data=high.rare.variants.data, labels=high.rare.variants.data$TYPE, line.color = "blue",r0=0, r1=0.3,cex=0.5,data.panel = 1)
kpDataBackground(kp, data.panel = 2, r0=0, r1=1)
kpAxis(kp, col="black", ymin=0, ymax=max(DFK$DBFREQ), r0=1, r1=0.5, data.panel = 2,tick.pos = c(0,max(DFK$DBFREQ)/4,max(DFK$DBFREQ)/2,3*max(DFK$DBFREQ)/4,max(DFK$DBFREQ)))
kpText(kp, data=high.rare.variants.data, x=0, y=0.75, labels="% (Percent)", col="red", srt=90,cex=0.8,data.panel=2)
kpText(kp, data=high.rare.variants.data, x=0, y=0.25, labels="Protein Coding", col="blue", srt=90,cex=0.8,data.panel=2)
kpLines(kp, data=high.rare.variants.data,x=high.rare.variants.data$START, y=DFK$DBFREQ/max(DFK$DBFREQ), r0=1, r1=0.5,data.panel = 2,col="red")
kpPoints(kp, data=high.rare.variants.data,x=high.rare.variants.data$START, y=DFK$DBFREQ/max(DFK$DBFREQ), r0=1, r1=0.5,cex=DFK$DBFREQ/max(DFK$DBFREQ),col=unlist(lapply(high.rare.variants.data$TYPE,colortype)),data.panel = 2)
if(dim(results)[1]>0){
	kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, line.color = "blue",r0=0, r1=0.4,data.panel = 2,cex=0.4)
}

dev.off()

























