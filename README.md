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

##### Get the help message

Navigate to 'EPIGENE' directory and type 'Rscript EPIGENE.R' or 'Rscript EPIGENE.R -h' to see EPIGENE usage:

>Usage: EPIGENE.R --genome [Genome build] --data [location of input information in tab seperated text format] --filePath [location of EPIGENE directory] --h [help]


#### Input parameters and data preparation

##### --genome

Currently EPIGENE supports hg19,hg38,mm9 and mm10.

##### --data

The only preparation that has to be done is to create a tab delimited info file that lists
the location of all ChIP-seq experiments in bam file format.\ All the bam files should be indexed, the .bai file should be present in the same location as the bam file.\

EPIGENE requires the following histone modifications for active transcription unit prediction:\
H3K27ac, H3K4me3, H3K4me1, H3K36m3, H3K27me3, H3K9me3.\
The info file should contain the location and sequencing details of these histone modifications. The required columns in info file are: "features", 'location" and "sequencing_info".

features        : list of histone modifcations
location        : location of the alignments in bam format, e.g.: '/project/projs-sahu/EPIGENE/DATA/H3K27ac.bam'\
sequencing_info : SE (for single-end sequencing) or PR (for paired-end sequencing)\

The  --data parameter requires the info file (specify the exact location of file if this file is not located in EPIGENE directory). An example info file (fileDetails_K562.txt) can be found in EXAMPLE folder.

##### --filePath

Location of EPIGENE folder.

#### Example run:
Rscript EPIGENE.R --genome hg19 --data EPIGENE/EXAMPLE/fileDetails_K562.txt --filePath /project/projs-sahu/EPIGENE

##### Output

> Bed file containing genomic co-ordinates of active transcription units: 'TU_predictions.bed'
> Heatmap of estimated parameters: 'estimated_params.pdf'

