library(vcfR)
library(stringr)
library(adegenet)
library(seqinr)

setwd('/Users/maggiehallerud/Marten_Primer_Design/MAF30_repBaseFiltered_01AUG2023/Random200/')

# load in full VCF and loci list
# load VCF
marten_vcf <- read.vcfR('../0_Inputs/13A_MAF30.recode.vcf')
marten_fix <- as.data.frame(marten_vcf@fix) #this file holds SNP positions 
marten_gt <- as.data.frame(marten_vcf@gt) #this file holds genotypes per sample- rows match marten fix
#names(marten_fix) <- c("CHROM","POS","ID","MAJOR","MINOR","UNK","UNK","ALLELE_FREQ")
# fix locus IDs
marten_fix$CHROM <- unlist(lapply(marten_fix$ID, function(X) strsplit(X,':')[[1]][1]))

# load template sequences
#marten_seq <- read_templates('../0_Inputs/13A_MAF30.recode.vcf') #input fasta file
#View(as.data.frame(marten_seq$Sequence))
loci_set <- read.csv('2_OptimizedSets/MAF30-repBaseFiltered-Random200-02Aug2023_optimizedSet16.csv', header=FALSE)
names(loci_set) <- c('Locus','Dimers')
loci <- unlist(lapply(loci_set$Locus, function(X) str_split(X, '_')[[1]][2]))

# subset VCF for SNPs within these loci
snp_subset_indx <- marten_fix$CHROM %in% loci
panel_results <- marten_vcf[snp_subset_indx]

# extract genotypes per individual
panel_genind <- vcfR2genind(panel_results)
View(panel_genind@tab)

# extract allele freqs- vcf AFs
#allele_freqs <- unlist(lapply(panel_fix$INFO, function(X) as.numeric(str_split(str_split(X,";")[[1]][2],"=")[[1]][2])))
freqs <- as.data.frame(makefreq(panel_genind))

# calculate total allele freqs per locus
allele_freqs <- colSums(freqs,na.rm=T)/nrow(freqs)
hist(allele_freqs)
sum(allele_freqs>=0.3 & allele_freqs<=0.7)/2

# allele freqs- coastal marten only
rownames(freqs)
coastal_freqs <- colSums(freqs[11:29,],na.rm=T)/19
#hist(coastal_freqs)
sum(coastal_freqs>=0.3 & coastal_freqs<=0.7) /2#because diploid

# Dunes freqs
dunes_freqs <- colSums(freqs[c(11:12,19:29),],na.rm=T)/13
#hist(dunes_freqs)
sum(dunes_freqs>=0.3 & dunes_freqs<=0.7)/2

# soOR freqs
soOR_freqs <- colSums(freqs[c(13:18),],na.rm=T)/6
#hist(soOR_freqs)
sum(soOR_freqs>=0.3 & soOR_freqs<=0.7)/2

# BMT freqs
bmt_freqs <- colSums(freqs[c(4:10),],na.rm=T)/7
#hist(bmt_freqs)
sum(bmt_freqs>=0.3 & bmt_freqs<=0.7)/2

# extract these primers and sequences
primers <- read.csv('MAF30-repBaseFiltered-Random200-02Aug2023_specificityCheck_passed.csv')
primers$PairID <- unlist(lapply(1:nrow(primers), function(X) paste(primers$LocusID[X],primers$PrimerPair[X], sep='_')))
primers_sub <- primers[primers$PairID%in%loci_set$Locus,]

fasta <- read_templates('../0_Inputs/13A_MAF30.repBaseFiltered.random200.fa')
fasta$CHROM <- unlist(lapply(fasta$Header, function(X) str_split(X, '_')[[1]][2]))
#amplicons <- fasta[fasta$CHROM %in% loci]
#write.table(out_IDT_order, file='TargetLoci_set22_03Aug2023.fa', row.names=FALSE, quote=FALSE)

out_IDT_order=data.frame(matrix(NA,nrow=0,ncol=4))
for (i in 1:nrow(primers_sub)){
  id <- primers_sub$PrimerID[i]
  seq <- toupper(primers_sub$Sequence[i])
  amt <- '25nm'
  std <- 'STD' #standard desalting
  out_IDT_order <- rbind(out_IDT_order, 
                         c(id, seq, amt, std))
}#
write.csv(out_IDT_order, 'IDT_Primer_Order_03Aug2023_Marten50loci.csv',row.names=FALSE,col.names=FALSE)
