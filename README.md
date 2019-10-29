# EPIGENE
EPIgenomic GENE prediction in R.
EPIGENE uses a multivariate HMM to predict active transcription units.
## Contact
sahua@staff.uni-marburg.de
## Installation
##### R PACKAGES

All required R packages are installed at the first run.\
Nothing needs to be installed manually.

> Overview of all packages:
> 
> bamsignals\
> kfoots\
> BSgenome.Mmusculus.UCSC.mm9 # when using genome mm9\
> BSgenome.Mmusculus.UCSC.mm10 # when using genome mm10\
> BSgenome.Hsapiens.UCSC.hg19 # when using genome hg19\
> BSgenome.Hsapiens.UCSC.hg38 # when using genome hg38\
> EPI.genefinder\
> GenomicRanges\
> ggplot2\
> GenomicFeatures\
> GenomicRanges\
> normr\
> Biostrings\
> RColorBrewer\
> GenomicAlignments
## Run EPIGENE

#### Get the help message

Navigate to 'EPIGENE' directory and type 'Rscript EPIGENE.R' or 'Rscript EPIGENE.R -h' to see EPIGENE usage:

>Usage: EPIGENE.R --genome [Genome build] --data [location of input information in tab seperated text format] --filePath [location of EPIGENE directory] --h [help]


#### Input parameters and data preparation

####### --genome

Currently EPIGENE supports hg19,hg38,mm9 and mm10.

####### --data
The only preparation that has to be done is to create a tab delimited info file that lists
the location of all ChIP-seq experiments in bam file format. The --data parameter requires the location of 

