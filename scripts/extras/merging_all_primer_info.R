setwd('/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/SNPpanel/Panel2_200initialpairs/')

# load dependencies
library(stringr)
library(vcfR)
library(adegenet)

# adding template info to primer pairs
pairs <- read.csv('multiplex_wormhole/4_OptimizedSets/Run02_150Loci_allLoci_dimers.csv',
                  colClasses=c("PrimerPairID"="character"))
templates <- read.csv('multiplex_wormhole/0_Inputs/CoastalMartenTemplates_MAF10_CENSOR_Martesmartes_trimmed.csv')
names(templates)
names(pairs)

pairs$SEQUENCE_ID <- unlist(lapply(pairs$PrimerPairID, function(x) strsplit(x, "\\.")[[1]][1]))

merged <- merge(pairs, templates, by="SEQUENCE_ID")
head(merged)

# add primer info
primers <- read.csv('multiplex_wormhole/4_OptimizedSets/Run02_150Loci_allLoci_primers.csv')
design <- read.csv('multiplex_wormhole/2_FilteredPrimers/SpecificityCheckTemplates_passed.csv')
names(primers) <- c("PrimerID", "Sequence")
primers$SEQUENCE_ID <- unlist(lapply(primers$PrimerID, function(x) strsplit(as.character(x), "\\.")[[1]][1]))
primers$Direction <- unlist(lapply(primers$PrimerID, function(x) strsplit(as.character(x), "\\.")[[1]][3]))
primers$PrimerPairID <- unlist(lapply(primers$PrimerID, function(x) {
  split <- strsplit(as.character(x), "\\.")[[1]]
  return (paste(split[1],split[2],sep="."))
}))
primers$Sequence <- str_replace_all(primers$Sequence, "tcgtcggcagcgtcagatgtgtataagagacag", "")
primers$Sequence <- str_replace_all(primers$Sequence, "gtctcgtgggctcggagatgtgtataagagacag", "")
head(primers)
merged_primers <- merge(primers, design, by="PrimerID")
merged_primers <- merged_primers[,c("PrimerID","PrimerPairID","SEQUENCE_ID", "Direction.x", "Sequence.x", "Length","AnnealingTempC","PropBound","AmpliconSize")]
names(merged_primers)[5] <- "PrimerSeq"
names(merged_primers)[6] <- "PrimerLen"

merged$FW <- NA
merged$FWtemp <- NA
merged$FWlen <- NA
merged$FWpropbound <- NA
merged$REV <- NA
merged$REVtemp <- NA
merged$REVlen <- NA
merged$REVpropbound <- NA
for (pair in merged$PrimerPairID){
  primersub = merged_primers[merged_primers$PrimerPairID==pair,]
  if (nrow(primersub)>0){
    merged$FW[merged$PrimerPairID==pair] <- primersub$PrimerSeq[which(primersub$Direction.x=="FW")]
    merged$FWtemp[merged$PrimerPairID==pair] <- primersub$AnnealingTempC[which(primersub$Direction.x=="FW")]
    merged$FWpropbound[merged$PrimerPairID==pair] <- primersub$PropBound[which(primersub$Direction.x=="FW")]
    merged$REV[merged$PrimerPairID==pair] <- primersub$PrimerSeq[which(primersub$Direction=="REV")]
    merged$REVtemp[merged$PrimerPairID==pair] <- primersub$AnnealingTempC[which(primersub$Direction.x=="REV")]
    merged$REVpropbound[merged$PrimerPairID==pair] <- primersub$PropBound[which(primersub$Direction.x=="REV")]
  }#if
}#pair
merged$FWlen <- nchar(merged$FW)
merged$REVlen <- nchar(merged$REV)

head(merged)

# calculate amplicon info
reverse_complement <- function(seq){
  # reverse
  rev <- rev(strsplit(seq, NULL)[[1]])
  # complement
  comp <- sapply(rev, function(base){
    switch(base, "a" = "t", "c" = "g", "g" = "c", "t" = "a")
    }#function
  )#sapply
  return(paste(comp,collapse=""))
}#reverse_complement

merged$Amplicon <- NA
merged$Probe <- NA
merged$FWstart <- NA
merged$FWend <- NA
merged$REVstart <- NA
merged$REVend <- NA
for (i in 1:nrow(merged)){
  merged$FWstart[i] <- str_locate(merged$SEQUENCE_TEMPLATE[i], merged$FW[i])[1]
  merged$FWend[i] <- str_locate(merged$SEQUENCE_TEMPLATE[i], merged$FW[i])[2]
  merged$REVstart[i] <- str_locate(merged$SEQUENCE_TEMPLATE[i], reverse_complement(merged$REV[i]))[1]
  merged$REVend[i] <- str_locate(merged$SEQUENCE_TEMPLATE[i], reverse_complement(merged$REV[i]))[2]
  merged$Amplicon[i] <- substr(merged$SEQUENCE_TEMPLATE[i], merged$FWstart[i], merged$REVend[i])
  merged$Probe[i] <- substr(merged$SEQUENCE_TEMPLATE[i], merged$FWstart[i], merged$REVstart[i]-1)
}#for

merged$AmpSize <- nchar(merged$Amplicon)


# calculate allele frequencies
vcf <- read.vcfR('multiplex_wormhole/0_Inputs/CoastalMartens.maf10.Mar2024.recode.vcf')
fixes <- as.data.frame(getFIX(vcf))

genind <- vcfR2genind(vcf)
maf <- minorAllele(genind)

merged$POS <- NA
merged$MAF <- NA
merged$Major <- NA
merged$Minor <- NA
for (row in 1:nrow(merged)){
  locus = merged$SEQUENCE_ID[row]
  snp = fixes$POS[which(fixes$CHROM==locus)]
  pos <- paste(locus, snp, sep=":")
  merged$POS[row] <- paste(snp, collapse=", ")
  major <- fixes$REF[which(fixes$ID%in%pos)]
  minor <- fixes$ALT[which(fixes$ID%in%pos)]
  af <- maf[rownames(maf)%in%pos]
  if (length(af)>1){
    if(length(unique(af))>1){
      merged$MAF[row] <- paste(af, collapse=", ")
      merged$Major[row] <- paste(major, collapse=", ")
      merged$Minor[row] <- paste(minor, collapse=", ")
    }else{
      merged$MAF[row] <- unique(af)
      merged$Major[row] <- str_flatten(major)
      merged$Minor[row] <- str_flatten(minor)
    }#ifelse
  }else{
    merged$MAF[row] <- af
    merged$Major[row] <- major
    merged$Minor[row] <- minor
  }#ifelse
}#for

# save CSV
View(merged)
write.csv(merged, 'multiplex_wormhole/4_OptimizedSets/FullPrimerInfo_Run12_150allloci.csv')

# check distributions
hist(as.numeric(merged$MAF), main="MAF", xlab="MAF")
hist(c(merged$FWtemp, merged$REVtemp), main="Annealing Temp", xlab="Temp C")
hist(merged$AmpSize, main="Amplicon Size", xlab="BP")
hist(c(merged$FWpropbound, merged$REVpropbound), main="Primer Binding", xlab="Proportion Bound")
hist(c(merged$FWlen, merged$REVlen), main="Primer Length", xlab="BP")


#### make GTseq primerInfo files ####
# grab out fields
primerinfo <- merged[,c("PrimerPairID", "FW", "REV")]
names(primerinfo)[1] <- "Locus"

# capitalize primer sequences
primerinfo$FW <- toupper(primerinfo$FW)
primerinfo$REV <- toupper(primerinfo$REV)

# add prefix to primer pair names (if desired)
primerinfo$PrimerPairID <- paste("MACA", primerinfo$PrimerPairID, sep='_')

# check & export
head(primerinfo)
write.csv(primerinfo, 'GTseq_PrimerInfo_150plex_newPrimers.csv', row.names=FALSE)


#### make GTseq locusInfo file ####
# function to extract probe seqs
extract_probes <- function(row){ #row = merged[index,]
  # check how many SNPs are in locus
  poss <- as.numeric(strsplit(row$POS, ",")[[1]])

  # if there is more than 1 SNP - choose which one to use
  if (length(pos>1)){
    mafs <- as.numeric(strsplit(row$MAF, ",")[[1]])
    # if maf is same- just choose the first
    if (length(mafs==1)){
      index <- 1
    # otherwise, use the SNP with highest MAF
    }else{
     index <- which(mafs == max(mafs))
    }#ifelse

  # otherwise, just use the one SNP
  }else{
    index <- 1
  }#ifelse
  
  # NOTE: This is the position in the template, so we'll have to adjust by where the FW primer starts
  pos <- poss[index]
  pos <- pos - row$FWstart
  
  # find 6 bp before and after SNP
  start <- pos - 6
  end <- pos + 6
  
  # check that 6 bp before SNP isn't in FW primer - adjust if so
  if (start <= row$FWlen){
    startdiff <- row$FWlen - start
    start <- row$FWlen + 1
    # adjust REV based on this diff so that probe length is ~13
    end <- end + startdiff
  }#if
  
  # check that 6 bp before SNP isn't in REV primer - adjust if so
  REVstart <- (row$AmpSize-row$REVlen+1) # this is where the REV primer starts in the amplicon, not be confused with where it starts in the template
  if (end >= REVstart){
    enddiff <- end - REVstart
    end <- REVstart - 1
  }#if
  
  # split amplicon into component bases
  split <- strsplit(row$Amplicon, "")[[1]]
  
  # extract major and minor alleles
  major <- strsplit(row$Major, ",")[[1]][index]
  minor <- strsplit(row$Minor, ",")[[1]][index]
  
  # create probe1
  probe1 <- split
  probe1[pos] <- major # substitute allele1 in position
  probe1 <- probe1[start:end] # extract subset of amplicon
  probe1 <- str_flatten(probe1) # put bases back into string
  probe1 <- toupper(probe1) # capitalize all
  
  # create probe2
  probe2 <- split
  probe2[pos] <- minor # substitute allele1 in position
  probe2 <- probe2[start:end] # extract subset of amplicon
  probe2 <- str_flatten(probe2) # put bases back into string
  probe2 <- toupper(probe2) # capitalize all
  
  # account for other SNPs that may occur in probe....
  snps <- poss[-index]
  i <- 1
  while(i<=length(snps)){
    # adjust position based on amplicon instead of template
    snp <- snps[i]
    snp <- snp - row$FWstart
    # check if snp is in probe region
    if (snp %in% start:end){
      # find position of snp in probe
      snppos <- which(start:end==snp)
      # grab alleles for SNP
      major <- strsplit(row$Major,",")[[1]][i]
      minor <- strsplit(row$Minor,",")[[1]][i]
      # convert to GTseq syntax
      alleles <- str_flatten(c("[", major, minor, "]"))
      # replace position in probes with GTseq SNP
      probe1 <- strsplit(probe1,"")[[1]]
      probe1[snppos] <- alleles
      probe1 <- str_flatten(probe1)
      probe2 <- strsplit(probe2,"")[[1]]
      probe2[snppos] <- alleles
      probe2 <- str_flatten(probe2)
    }#if
    #update iterator
    i=i+1
  }#while

  # return probes
  return(c(probe1,probe2))
}#function extract_probes


# set up output dataframe
locusinfo <- data.frame(LocusID=merged$PrimerPairID, Probe1=NA, Probe2=NA, FW_primer=toupper(merged$FW))

# grab probe info for all rows
for (row in 1:nrow(merged)){
  locusinfo[row, c("Probe1","Probe2")] <- extract_probes(merged[row,])
}#for

# add prefix to primer pair names (if desired)
locusinfo$LocusID <- paste("MACA", locusinfo$LocusID, sep='_')

# check & export
head(locusinfo)
write.csv(primerinfo, 'GTseq_LocusInfo_150plex_newPrimers.csv', row.names=FALSE)
