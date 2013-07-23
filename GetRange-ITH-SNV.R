##################################################
##  Script to get sequence (which can be SNP masked)    
##  around SNV or indels to design primers for   		
##	 Targeted resequencing on the MiSeq 				
##  Aparicio Lab WSOP 2013-001 developed by			
##	 Dr Damian Yap , Research Associate				
##	 dyap@bccrc.ca  Version 2.0 (Jul 2013) 			
##################################################


# These commands must be specifed in order for this script to work
# source("http://www.bioconductor.org/biocLite.R"); source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); library('BSgenome.Hsapiens.UCSC.hg19')

library(Biostrings)
library("IRanges")
library("GenomicRanges")
library(Rsamtools)

#################################################
# Directory structure - uncomment for first running of script
Project="ITH"
# homebase="/Users/dyap/Documents/Breast Cancer" # MOMAC14
homebase="/home/dyap/Projects" # beast

setwd(homebase)
# system('mkdir ITH')
setwd(paste(homebase,Project,sep="/"))
# system('mkdir primer3')
# system('mkdir positions')
# system('mkdir Annotate')
getwd()

#######################################
# Save input files under $homebase/positions#
#######################################

# Check if the hg19 library is install, if not install it and load it
genome="BSgenome.Hsapiens.UCSC.hg19"
library('BSgenome.Hsapiens.UCSC.hg19')
library(SNPlocs.Hsapiens.dbSNP.20120608)

# Load the latest available dnSNP for the version of bioconductor if not installed 
len <-  length(available.SNPs())
dbSNP <- available.SNPs()[len]
SNP <-   installed.SNPs()
# If dbSNP is not currently installed load the latest version from source
# if (!dbSNP %in%  SNP) source("http://www.bioconductor.org/biocLite.R"); biocLite(dbSNP)

# Inject the SNPs into the hg19 ref
SNP_Hsapiens <- injectSNPs(Hsapiens, dbSNP)

##############################################
######            User defined variables               ######
# Directory and file references
basedir=paste(homebase,Project,sep="/")
sourcedir=paste(basedir,"positions", sep="/")
# outdir=paste(basedir,"positions", sep="/")
p3dir=paste(basedir,"primer3", sep="/")
outpath=paste(basedir,"Annotate", sep="/")

######################
# These are the input files
#type="indel"
type="SNV"

file="ITH_pos.txt"

p3file=paste(type,"design.txt",sep="_")

annofile = paste(outpath, paste(type, "Annotate.csv", sep="_") ,sep="/")

outfile=paste(p3dir,p3file,sep="/")
input=paste(sourcedir,file,sep="/")

# offsets (sequences on either side of SNV,indel for design space)
offset=200
WToffset=5

# Select the appropriate Genome (mask) - for reference only
#BSg="Hsapiens" # normal-reference
#BSg="SNP_Hsapiens" # SNP-hard masked genome

##############################################

indf <- read.table(file=input,  stringsAsFactors = FALSE, header=TRUE)

                  
outdf <- data.frame(ID = rep("", nrow(indf)),
		                 Chr = rep("", nrow(indf)),
                     Start = rep(0, nrow(indf)),
                     End = rep(0, nrow(indf)),
                     Indel = rep("", nrow(indf)),
                     Context = rep("", nrow(indf)),
                     Design = rep("", nrow(indf)),
                     stringsAsFactors = FALSE)
                     
                     
for (ri in seq(nrow(indf))) {

  # for format  2:139318499-139318499
   chr <- paste("chr", strsplit(indf[ri,1], split=":")[[1]][1], sep="")
   start <- as.numeric(strsplit(strsplit(indf[ri,1], split=":")[[1]][2], split = "-")[[1]][1])
   end <- as.numeric(strsplit(strsplit(indf[ri,1], split=":")[[1]][2], split = "-")[[1]][2])

  if ( start == end ) id = paste(chr,start,sep="_")
  if ( start != end ) id = paste(paste(chr,start,sep="_"), end, sep="-")

# The indel or SNV has to match exactly so we get the reference sequences from hg19
 idseq <- as.character(getSeq(Hsapiens,chr,start,end))
 
 # If the index <5 or it is an SNV then we need 5 bp on either side to match (for visualization only)
 if (nchar(idseq) < 5) cxt <- as.character(getSeq(Hsapiens,chr,start-WToffset,end+WToffset)) else cxt <- idseq
 
 # This design space comprises of upstream and downstream annotated sequence and the exact reference seq
 # for hg19 for the SNV position and/or indel (This is important for matching)
  dseq <- paste(paste(as.character(getSeq(SNP_Hsapiens,chr,start-offset,start-1)),idseq,sep=""),
              as.character(getSeq(SNP_Hsapiens,chr,end+1,end+offset)),sep="")
              
  outdf$ID[ri] <- id          
  outdf$Chr[ri] <- chr
  outdf$Start[ri] <- start
  outdf$End[ri] <- end
  outdf$Indel[ri] <- idseq
  outdf$Context[ri] <- cxt
  outdf$Design[ri] <- dseq
  }

# Output file ID, chr, start, end, indel sequence, context for matching, design space seq for primer3 design
write.csv(outdf, file = outfile)


######################

# For annotation files

andf1 <- data.frame(Chr1 = rep("", nrow(indf)),
                     Pos1 = rep(0, nrow(indf)),
                     Pos2 = rep(0, nrow(indf)),
                     WT1 = rep("", nrow(indf)),
                     SNV1 = rep("", nrow(indf)),
                     stringsAsFactors = FALSE)
        
                     
for (ri in seq(nrow(indf))) {
# For GetSeq we need the Chr prefix for chromosome but not for ANNOVAR  
# Assume indels are <40bp on the same chromosome

  # for format  2:139318499-139318499
   chr <- paste("chr", strsplit(indf[ri,1], split=":")[[1]][1], sep="")
   chrom <- strsplit(indf[ri,1], split=":")[[1]][1]
   pos1 <- as.numeric(strsplit(strsplit(indf[ri,1], split=":")[[1]][2], split = "-")[[1]][1])
   pos2 <- as.numeric(strsplit(strsplit(indf[ri,1], split=":")[[1]][2], split = "-")[[1]][2])
   
   if ( pos1 == pos2 ) id = paste(chr,pos1,sep="_")
   if ( pos1 != pos2 ) id = paste(paste(chr,pos2,sep="_"), pos2, sep="-")
   
   # It indel is >50 bp, the middle is chosen to check annotation
   diff=pos2-pos1
   if ( diff >= 50 )  pos=pos1+(pos2-pos1)/2 else pos=pos1 
   
   wt <- as.character(getSeq(SNP_Hsapiens,chr,pos,pos))
    
# Fake the SNV to be just the complement of WT position (as SNV allele is not known)

if (wt=="A") snv1 <- "T"
if (wt=="C") snv1 <- "G"
if (wt=="G") snv1 <- "C"
if (wt=="T") snv1 <- "A"
  
  andf1$Chr1[ri] <- chrom
  andf1$Pos1[ri] <- pos
  andf1$Pos2[ri] <- pos
  andf1$WT1[ri] <- wt
  andf1$SNV1[ri] <-snv1


  }
  
  
# Format for ANNOVAR  <15 43762161 43762161 T C>
write.csv(andf1, file = annofile )

