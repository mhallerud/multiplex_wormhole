# set up env
library(openPrimeR)
library(vcfR)
library(stringr)


#### READ IN SEQUENCES
setwd('/Users/maggiehallerud/Marten_Primer_Design/Plate1_First55Pairs_Sep2023/')
# load template sequences
marten_seq <- read_templates('MAF30_repBaseFiltered_01AUG2023/0_Inputs/13A_MAF30.repBaseFiltered.random200.fa') #input fasta file
#View(as.data.frame(marten_seq$Sequence))

# load VCF
marten_vcf <- read.vcfR('populations.snps.vcf')
marten_fix <- as.data.frame(marten_vcf@fix) #this file holds SNP positions 
#names(marten_fix) <- c("CHROM","POS","ID","MAJOR","MINOR","UNK","UNK","ALLELE_FREQ")

# fix locus IDs
marten_fix$CHROM <- unlist(lapply(marten_fix$ID, function(X) strsplit(X,':')[[1]][1]))

# reset allowed primer binding locations- check VCF for position of SNPs on each sequence
# go through one at a time & fix...
for (i in 1:nrow(marten_seq)){
  id <- str_split(marten_seq$ID[i], '_')[[1]][2]
  snps <- marten_fix$POS[marten_fix$CHROM==as.numeric(id)]
  snps <- as.numeric(snps)
  min_snp <- min(snps)
  max_snp <- max(snps)
  if (min_snp==Inf) min_snp <- marten_seq$Sequence_Length[i]
  if (max_snp==-Inf) max_snp <- 1
  marten_seq$Allowed_End_fw[i] <- min_snp-1
  marten_seq$Allowed_Start_rev[i] <- max_snp+1
  marten_seq$Allowed_End_rev[i] <- marten_seq$Sequence_Length[i]
}#for i

# check all are populated
marten_seq$Allowed_Start_fw
marten_seq$Allowed_End_fw
marten_seq$Allowed_Start_rev
marten_seq$Allowed_End_rev
#marten_seq$Allowed_Start_fw <- 1

# adjust allowed binding regions
for (i in 1:nrow(marten_seq)){
  marten_seq$Allowed_fw[i] <- substr(marten_seq$Sequence[i], marten_seq$Allowed_Start_fw[i], marten_seq$Allowed_End_fw[i])
  marten_seq$Allowed_rev[i] <- substr(marten_seq$Sequence[i], marten_seq$Allowed_Start_rev[i], marten_seq$Allowed_End_rev[i])
}#for i

# check allowed binding regions
marten_seq$Allowed_fw
marten_seq$Allowed_rev

# double check that microhaps are properly represented
counts=data.frame(table(marten_fix$CHROM))
microhaps=counts[counts$Freq>1,]
marten_fix[marten_fix$CHROM==microhaps$Var1[1],]
marten_seq[marten_seq$ID==">CLocus_1",]

# export microhaplotypes
nrow(microhaps)
microhaps$Var1 <- paste0('>CLocus_', microhaps$Var1)
microhaps <- microhaps[order(microhaps$Freq),]
microhaps_seq <- marten_seq[which(marten_seq$ID %in% microhaps$Var1)]
microhaps_seq <- microhaps_seq[which(nchar(microhaps_seq$Allowed_rev)>=18 | nchar(microhap_seq$Allowed_fw)>=18),]
nrow(microhaps)
write.csv(microhaps_seq$ID, 'MartenTemplates_MAF30-repBaseFilter_01Aug2023_microhaplotypes.csv', row.names=FALSE)

# remove any that don't have enough binding space
marten_seq <- marten_seq[which(nchar(marten_seq$Allowed_fw)>=18 | nchar(marten_seq$Allowed_rev)>=18),]

# check that there's only one line per locus
nrow(marten_seq)
length(unique(marten_seq$Header))


# double check that everything looks OK...
View(marten_seq)

# export IDs, templates, targets as CSV for use with primer3
target_start <- marten_seq$Allowed_End_fw + 1
target_end <- marten_seq$Allowed_Start_rev
target_len <- target_end - target_start
targets <- paste0(as.character(target_start), ",", as.character(target_len))
marten_csv <- data.frame(SEQUENCE_ID=str_replace_all(marten_seq$ID,'>',''),
                         SEQUENCE_TEMPLATE=marten_seq$Sequence,
                         SEQUENCE_TARGET=targets)
write.csv(marten_csv, 'MartenTemplates_MAF30-repBaseFilter-random300.csv', row.names=FALSE)                         

