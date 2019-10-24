#####################################################################################
#check for package installation
#####################################################################################
pkgStatus <- function(pkg){
  if (!pkg %in% installed.packages()) {
    cat(paste0("package ",pkg, " not found. Install now."))
    return(F)
  }else{
    return(T)	
  }
} 

#####################################################################################
#check for package loading
#####################################################################################
pkgLoad <- function(pkg) {
  if(pkgStatus(pkg) == F){
    installPkg(pkg)
  }else{
    cat(paste0("Loading package ", pkg))
    suppressMessages(library(pkg, character.only = TRUE))
    done()
  }
}

#####################################################################################
#install package
#####################################################################################
installPkg <- function(pkg) { 
  
  if(! isTRUE(pkg %in% .packages(all.available=T))) { 
    
    x <- tryCatch({eval(parse(text = sprintf("install.packages(\"%s\", dependencies = T)", pkg)))})
    
    if(is.null(x) & pkg == "EPI.genefinder") {
      print(pkg)
      x <- tryCatch({devtools::install_github("imbeLab/EPI.genefinder")})
    }
    
    if(is.null(x) & pkg == "kfoots") {
      print(pkg)
      x <- tryCatch({devtools::install_github("lamortenera/kfoots")})
    }
    
    if(is.null(x)) {
      x <- tryCatch({source("http://bioconductor.org/biocLite.R")
        eval(parse(text = sprintf("biocLite(\"%s\")", pkg)))})
    }
    
    if(is.null(x)) {
      cat(paste0("Unable to install package: ",pkg,".\n"));
      q();
    }else{
      eval(parse(text = sprintf("require(\"%s\")", pkg)))
    }
  }
}

#####################################################################################
#remove long stretches of N from genome
#####################################################################################
seq.check <- function(given.seq, chr, tmp.genome){
  # subset BSgenome to chromosome of interest:
  where <- which(tmp.genome@seqinfo@seqnames == chr)
  current.chr <- tmp.genome[[where]]
  # get all coordinates that match the pattern:
  tmp.match <- matchPattern(given.seq, getSeq(tmp.genome, chr, start = 1, end = current.chr@length))
  # transform to GRanges and merge adjacent ranges:
  tmp.granges <- GenomicRanges::reduce( GRanges(seqnames = chr, ranges = ranges(tmp.match)))
  return(tmp.granges)
}

#####################################################################################
#genome binning
#####################################################################################
binning <- function(chrBins,widthVal){
  gChrList <- do.call("c",chrBins)
  seqlevels(gChrList) <- sub('chr','',seqlevels(gChrList))
  gHuman <- GRanges(c(1:22,"X","Y"),IRanges(rep(1,23),seqLengths))
  gBinHuman <- unlist(tile(gHuman,width = widthVal))
  nBins <- findOverlaps(gBinHuman,gChrList,type = "within")
  gInvalid <- gBinHuman[queryHits(nBins)]
  gValid <- gBinHuman[setdiff(c(1:length(gBinHuman)),queryHits(nBins))]
  gContigs <-reduce(gValid)
  binData = GRangesList(bins = gValid, contigs = gContigs)
  return(binData)
}

#####################################################################################
#output messages
#####################################################################################
startMsg <- function(m) {cat(paste0("\n--- ",m," ---\n\n"))}
endMsg <- function() {cat("\n\t>>> All done!\n")}
done <- function() {cat("...done\n")}

#####################################################################################
#get counts from bam file
#####################################################################################
getBamCount <- function(trFPath,seqD,bins){
  bf = BamFile(trFPath)
  si = seqinfo(bf)
  # fix chromosome prefix
  print(seqlevels(si)[1])
  gr = bins
  seqlevels(gr) = gsub("chr","",seqlevels(gr))
  if (grepl("chr",seqlevels(si)[1])) {
    seqlevels(gr) = paste0("chr", seqlevels(gr))
  }
  
  gr = gr[seqnames(gr) %in% seqnames(si)]
  cnt = NULL
  if(seqD == "SE"){
    cnt = bamCount(trFPath, gr, mapqual = 30, filteredFlag = 1024, paired.end = "ignore")
  }else if(seqD == "PE"){
    cnt = bamCount(trFPath, gr, mapqual = 30, filteredFlag = 1024, paired.end = "midpoint")
  }
  return(cnt)
}

#####################################################################################
#generate coverage list
#####################################################################################
makeTable <- function(txID){
  txStart = tssBT[[txID]]
  txEnd = ttsBT[[txID]]
  nExon = length(ebt[[txID]])
  binID = txStart:txEnd
  nBins = length(binID)
  cov = data.frame()
  if (nExon != 1 & nBins > 2){
    txStr = txStrand[txID]
    
    res = data.frame(bins = binID, tss = integer(nBins), firstExon = integer(nBins), firstIntron = integer(nBins), iExon = integer(nBins), iIntron = integer(nBins), lastExon = integer(nBins), tts = integer(nBins), row.names = binID)
    
    res$tss[1] = 200
    res$tts[nBins] = 200
    firstExon = 1;
    lastExon = nExon
    
    exonBins = as.character(binsExonBT[[txID]])
    exonRank = ExonBinsBT[[txID]]
    exonRank = exonRank - min(exonRank) + 1
    exonBinsCov = binsExonCovBT[[txID]]
    
    ## first exon
    res[exonBins[exonRank == 1], "firstExon"] = exonBinsCov[exonRank == 1]
    
    
    ## internal exons
    res[exonBins[exonRank > 1 & exonRank < lastExon], "iExon"] = exonBinsCov[exonRank > 1 & exonRank < lastExon]
    
    
    ## last exon
    res[exonBins[exonRank == lastExon], "lastExon"] = exonBinsCov[exonRank == lastExon]
    intronBins = as.character(binsIntronBT[[txID]])
    intronRank = IntronBinsBT[[txID]]
    intronRank = intronRank - min(intronRank)
    intronBinsCov = binsIntronCovBT[[txID]]
    
    if (txStr == "+"){
      intronRank = intronRank + 1;
    } else{
      intronRank = max(intronRank) - intronRank + 1
    }
    nIntron = max(intronRank)
    
    ## first intron
    res[intronBins[intronRank == 1], "firstIntron"] = intronBinsCov[intronRank == 1]
    
    ## internal introns
    ## might be two of them
    intronBins = intronBins[intronRank != 1]
    intronBinsCov = intronBinsCov[intronRank != 1]
    intronBinsCov = tapply(intronBinsCov, intronBins, sum)
    res[names(intronBinsCov), "iIntron"] = intronBinsCov
    cov = res
  } else{
    cov = NULL
  }
  return(cov)
}

#####################################################################################
#obtain the bin information of transcription units
#####################################################################################

mapTssTts <- function(tssStrand,ttsStrand,mapData){
  bin = append(tssStrand,ttsStrand)
  contig = mapData[bin,"contig"]
  df = data.frame("bin" = bin,"contig" = contig)
  tssStrand = data.frame(tssStrand)
  colnames(tssStrand) = "bin"
  tssStrand["site"] = "TSS"
  ttsStrand = data.frame(ttsStrand)
  colnames(ttsStrand) = "bin"
  ttsStrand["site"] = "TTS"
  data = rbind(tssStrand,ttsStrand)
  data = merge(data,df,by="bin")
  final = data[data$contig %in% data[duplicated(data$contig),"contig"],]
  tss = final[which(final$site == "TSS"),]
  tts = final[which(final$site == "TTS"),]
  mapping = data.frame()
  for(i in 1:nrow(tss)){
    mapping[i,"TSS"]=tss[i,"bin"]
    ttsData = tts[which((tts$contig == tss[i,"contig"])&(tts$bin > tss[i,"bin"])),"bin"]
    if(length(ttsData)>0){
      mapping[i,"TTS"]=ttsData[1]
    }
    else{
      mapping[i,"TTS"]=0
    }
  }
  
  mapping = mapping[which(mapping$TTS != 0),]
  return(mapping)
}