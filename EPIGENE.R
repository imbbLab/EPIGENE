##################################################################
# get script path
##################################################################

args <- commandArgs(trailingOnly = TRUE)

if(("--h" %in% args) & (length(as.character(args)) == 1)){
  cat("Usage: EPIGENE.R --genome [Genome build] --data [location of input information in tab seperated text format] --filePath [location of EPIGENE directory] --cores [number of cores] --h [help]\n")
  q(save = "no", status = 0, runLast = TRUE)
}
if(length(args)>1){
  genDat = strsplit(grep('--genome', args, value = TRUE), split = '=')[[1]][[2]]
  seqDat = strsplit(grep('--data', args, value = TRUE), split = '=')[[1]][[2]]
  fileLoc = strsplit(grep('--filePath', args, value = TRUE), split = '=')[[1]][[2]]
  nCores = strsplit(grep('--cores', args, value = TRUE), split = '=')[[1]][[2]]
}

#####################################################################################
#load/install all the relevant packages
#####################################################################################

source(paste0(fileLoc,"BIN/","functions.R"))

pkgLoad("ggplot2")
pkgLoad("GenomicFeatures")
pkgLoad("Biostrings")
pkgLoad("GenomicRanges")
pkgLoad("EPI.genefinder")
pkgLoad("normr")
pkgLoad("bamsignals")
pkgLoad("RColorBrewer")

if (genDat == "mm10") pkg <- "BSgenome.Mmusculus.UCSC.mm10"
if (genDat == "mm9") pkg  <- "BSgenome.Mmusculus.UCSC.mm9"
if (genDat == "hg19") pkg <- "BSgenome.Hsapiens.UCSC.hg19"
if (genDat == "hg38") pkg <- "BSgenome.Hsapiens.UCSC.hg38"

pkgLoad(pkg)
assign("txdb", eval(parse(text = pkg)))

#####################################################################################
#genome binning
#####################################################################################
cat("Binning the genome............")
if((genDat == "hg19")|(genDat == "hg38")) chr <- as.character(c(1:22,"X","Y"))
if((genDat == "mm9")|(genDat == "mm10")) chr <- as.character(c(1:19,"X","Y"))

chrBSVal = paste0("chr",chr)
seqLengths = unlist(lapply(chrBSVal,function(x){
  tmp.genome = txdb
  where <- which(tmp.genome@seqinfo@seqnames == x)
  current.chr <- tmp.genome[[where]]
  current.chr@length
}))

suppressWarnings(chrBins <- lapply(chrBSVal, seq.check, given.seq = strrep("N",25),tmp.genome=txdb))
suppressWarnings(binInfo <- binning(chrBins,200,seqLengths))
bins = binInfo[[1]]
seBins = bins
seqlevels(seBins) = paste0("chr",seqlevels(seBins))
gContigs = binInfo[[2]]

#####################################################################################
#obtain the read counts for genomic bins
#####################################################################################
cat("Obtain the read counts in genomic bins............")
featureInfo = read.table(seqDat,sep = "\t",col.names = c("features","location","sequencing"))

features = as.character(featureInfo$features)
trFPath = as.character(featureInfo$location)
seqD = as.character(featureInfo$sequencing)

counts = mclapply(1:length(features), function(x) {
  if(seqD[x]== "PE")
    bamCount(trFPath[x], bins, mapqual = 30, filteredFlag = 1024, paired.end = "midpoint")
  else
    bamCount(trFPath[x], seBins, mapqual = 30, filteredFlag = 1024, paired.end = "ignore")
}, mc.cores = nCores)


names(counts) = features

#####################################################################################
#obtain binarized enrichments
#####################################################################################
cat("Binarizing the enrichment values of histone modifications............")
runNormR = list(
  H3K27ac = c("H3K27ac", "HM.input"),
  H3K27me3 = c("H3K27me3", "HM.input"),
  H3K36me3 = c("H3K36me3", "HM.input"),
  H3K4me1 = c("H3K4me1", "HM.input"),
  H3K4me3 = c("H3K4me3", "HM.input"),
  H3K9me3 = c("H3K9me3", "HM.input")
)

norm = mclapply(runNormR, function(run) {
  enrichR(treatment = counts[[run[1]]], control = counts[[run[2]]], genome = bins,binFilter = "zero")
}, mc.cores = nCores)

clzz = sapply(norm, function(x) { cl = normr::getClasses(x, fdr = 0.2); cl[!is.finite(cl)] = 0; cl})

#####################################################################################
#create updateEmission and updateTransition
#####################################################################################
cat("Preparing the input parameters............")
emissionProb = readRDS(paste0(fileLoc,"/DATA/","emissionProb.rds"))
transitionMat = readRDS(paste0(fileLoc,"/DATA/","transitionMat.rds"))
updateEmission = c(rep(FALSE, 14), rep(TRUE, 3))
updateTransition = matrix(FALSE, ncol = 17, nrow = 17)
updateTransition[1, 15:17] = TRUE
updateTransition[15:17, 7] = TRUE
updateTransition[14, 15:17] = TRUE
updateTransition[15:17, 8] = TRUE

updateTransition[15:17, 15:17] = TRUE

updateTransition[14,7] = TRUE
updateTransition[1,7] = TRUE
updateTransition[14,8] = TRUE
updateTransition[1,8] = TRUE

#####################################################################################
#obtain length of contigs
#####################################################################################

gMap <- findOverlaps(bins,gContigs,type = "within")

qHits <- queryHits(gMap)
sHits <- subjectHits(gMap)
mapContigs <- data.frame(qHits,sHits)
colnames(mapContigs) <- c("QueryHits","SubjectHits")

seqLenList <- c()

for(i in 1:length(gContigs)){
  seqLenList <- c(seqLenList,nrow(mapContigs[which(mapContigs$SubjectHits==i),]))
}

names(seqLenList) <- c(1:length(gContigs))

#####################################################################################
#predicting transcription units
#####################################################################################
cat("Predicting active transcription units............")
train = EPI.genefinder(clzz[, 1:6], emissionProb, NULL, t(transitionMat), seqLenList, updateEmission,
                       t(updateTransition), maxiter = 500, nthreads = 24)
labels = c("tss", "firstExon", "firstIntron", "iExon", "iIntron", "lastExon", "tts", "rtss",
           "rfirstExon", "rfirstIntron", "riExon", "riIntron", "rlastExon", "rtts", "bg1", "bg2", "bg3")

tssFwd <- which(train$viterbi$vpath == 1)
ttsFwd <- which(train$viterbi$vpath == 7)
tssRev <- which(train$viterbi$vpath == 14)
ttsRev <- which(train$viterbi$vpath == 8)

mapData = data.frame("bin"= queryHits(gMap),"contig"=subjectHits(gMap))
transFwd <- data.frame(mapTssTts(tssFwd,ttsFwd,mapData),strand = rep("+",nrow(mapTssTts(tssFwd,ttsFwd,mapData))))

transRev <- data.frame(mapTssTts(tssRev,ttsRev,mapData),strand = rep("+",nrow(mapTssTts(tssRev,ttsRev,mapData))))

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
export.bed(gPredictions,con = paste0(fileLoc,"TU_predictions.bed"))

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
estParam = ggplot(plotEM,aes(HistoneMarks,States))+geom_tile(aes(fill = Probability),colour = "black")+scale_fill_gradient(low = "white",high = "darkgreen")+
  theme_bw()+labs(y = "States")+
  theme(axis.text.x = element_text(size = 10,angle = 45, hjust = 1),axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10),legend.title = element_text(size = 10))
ggsave(estParam,filename = paste0(fileLoc,"estimated_params.pdf"))

cat("Done............")
