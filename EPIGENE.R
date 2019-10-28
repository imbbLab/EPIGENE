##################################################################
# get script path
##################################################################

args <- commandArgs(trailingOnly = TRUE)

genDat = strsplit(grep('--genome', args, value = TRUE), split = '=')[[1]][[2]]
seqDat = strsplit(grep('--data', args, value = TRUE), split = '=')[[1]][[2]]
fileLoc = strsplit(grep('--filePath', args, value = TRUE), split = '=')[[1]][[2]]


#####################################################################################
#load all the relevant packages
#####################################################################################

pkgLoad("ggplot2")
pkgLoad("GenomicFeatures")
pkgLoad("BSgenome.Hsapiens.UCSC.hg19")
pkgLoad("Biostrings")
pkgLoad("GenomicRanges")
pkgLoad("EPI.genefinder")
pkgLoad("normr")
pkgLoad("bamsignals")

#####################################################################################
#genome binning
#####################################################################################

chr = as.character(c(1:22,"X","Y"))
chrBSVal = paste0("chr",chr)
seqLengths = unlist(lapply(chrBSVal,function(x){
  tmp.genome = BSgenome.Hsapiens.UCSC.hg19
  where <- which(tmp.genome@seqinfo@seqnames == x)
  current.chr <- tmp.genome[[where]]
  current.chr@length
}))

chrBins <- lapply(chrBSVal, seq.check, given.seq = strrep("N",25),tmp.genome=BSgenome.Hsapiens.UCSC.hg19)
binInfo = binning(chrBins,200,seqLengths)
bins = binInfo[[1]]
seBins = bins
seqlevels(seBins) = paste0("chr",seqlevels(seBins))
gContigs = binInfo[[2]]

#####################################################################################
#obtain the read counts for genomic bins
#####################################################################################

featureInfo = read.table(seqDat,header = T,sep = "\t",col.names = c("features","location","sequencing"))

features = as.character(featureInfo$features)
trFPath = as.character(featureInfo$location)
seqD = as.character(featureInfo$sequencing)

counts = mclapply(1:length(features), function(x) {
  if(seqD[x]== "PE")
    bamCount(trFPath[x], bins, mapqual = 30, filteredFlag = 1024, paired.end = "midpoint")
  else
    bamCount(trFPath[x], seBins, mapqual = 30, filteredFlag = 1024, paired.end = "ignore")
}, mc.cores = 10)


names(counts) = features

#####################################################################################
#obtain binarized enrichments
#####################################################################################

runNormR = list(
  H3K27ac = c("H3K27ac", "HM.input"),
  H3K27me3 = c("H3K27me3", "HM.input"),
  H3K36me3 = c("H3K36me3", "HM.input"),
  H3K4me1 = c("H3K4me1", "HM.input"),
  H3K4me3 = c("H3K4me3", "HM.input"),
  H3K9me3 = c("H3K9me3", "HM.input"),
  Pol2 = c("Pol2", "Pol2.input")
)

norm = mclapply(runNormR, function(run) {
  enrichR(treatment = counts[[run[1]]], control = counts[[run[2]]], genome = bins,binFilter = "zero")
}, mc.cores = 10)

clzz = sapply(norm, function(x) { cl = normr::getClasses(x, fdr = 0.2); cl[!is.finite(cl)] = 0; cl})

#####################################################################################
#predicting transcription units
#####################################################################################

train = EPI.genefinder(clzz[, 1:6], emissionProb, NULL, t(transitionMat), seqLenList, updateEmission,
                       t(updateTransition), maxiter = 500, nthreads = 24)
labels = c("tss", "firstExon", "firstIntron", "iExon", "iIntron", "lastExon", "tts", "rtss",
           "rfirstExon", "rfirstIntron", "riExon", "riIntron", "rlastExon", "rtts", "bg1", "bg2", "bg3")

mapData = data.frame("bin"= queryHits(gMap),"contig"=subjectHits(gMap))
transFwd <- mapTssTts(tssFwd,ttsFwd,mapData)
transFwd$Strand<-"+"
transRev <- mapTssTts(tssRev,ttsRev,mapData)
transRev$Strand<-"-"

tmpData <- rbind(transFwd,transRev)
rownames(tmpData)=c(1:nrow(tmpData))

#####################################################################################
#obtain the genomic co-ordinates of transcription units
#####################################################################################

startList <- start(bins)
endList <- end(bins)
chrList <- seqnames(bins)

startDat <- startList[tmpData$TSS]
endDat <- endList[tmpData$TTS]
chrDat <- chrList[tmpData$TSS]

gPredictions <- GRanges(chrDat,IRanges(start = startDat,end = endDat),strand = tmpData$Strand)

#####################################################################################
#visualise estimated parameters
#####################################################################################

emission = data.frame(train$epar)
labelsY = c("tss", "firstExon", "firstIntron", "iExon", "iIntron", "lastExon",
            "tts", "rtss", "rfirstExon", "rfirstIntron", "riExon", "riIntron", "rlastExon", "rtts", "bg1", "bg2", "bg3")
labelsX = colnames(clzz)[1:6]
colnames(emission) = labelsY
emission[,"Histone Marks"] = labelsX
plotEM = reshape2::melt(emission)
colnames(plotEM) = c("HistoneMarks","States","Probability")
myPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
ggplot(plotEM,aes(HistoneMarks,States))+geom_tile(aes(fill = Probability),colour = "black")+scale_fill_gradient(low = "white",high = "darkgreen")+
  theme_bw()+labs(y = "States")+
  theme(axis.text.x = element_text(size = 10,angle = 45, hjust = 1),axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10),legend.title = element_text(size = 10))
