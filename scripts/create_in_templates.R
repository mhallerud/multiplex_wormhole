# load dependencies
library(openPrimeR)
library(vcfR)

# set input variables
fasta <- "MAF30All.CENSORed.fa"
invcf <- "MAF30AllGrayFoxPops.Filtered.recode.vcf"
type <- "refbased"
#type <- "denovo"
minflankbp <- 18 # minimum size of flanking regions (i.e., min primer binding site size)
maxSNP <- 3 # maximum number of SNPs allowed in a template
## TODO: Consider adjusting this to a SNP density, e.g. for variable template lengths


#### READ IN SEQUENCES & SNPs ####
# load template sequences
templates <- openPrimeR::read_templates(fasta) #input fasta file
templates$ID <- gsub(">","", templates$ID)
#View(as.data.frame(templates$Sequence))

# load VCF
vcf <- vcfR::read.vcfR(invcf)
fix <- as.data.frame(vcf@fix) #this file holds SNP positions 



#### NOTE: Before proceeding, make sure that your locus names match with the "CHROM" field in the VCF ####
# For de novo Stacks output, the CHROM can be set based on the SNP ID field:
if(type=="denovo"){
  fix$CHROM <- unlist(lapply(fix$ID, function(X) strsplit(X,':')[[1]][1]))
}#if

# For reference-based output with loci in the format CHROM:startbp-endbp:
if(type=="refbased"){
  # separate the fasta info
  lociinfo <- data.frame(ID=templates$ID)
  lociinfo$CHROM <- unlist(lapply(lociinfo$ID, function(X) strsplit(X,':')[[1]][1]))
  lociinfo$CHROM <- gsub(">","",lociinfo$CHROM)
  lociinfo$RangeBP <- unlist(lapply(lociinfo$ID, function(X) strsplit(X,':')[[1]][2]))
  lociinfo$StartBP <- as.numeric(lapply(lociinfo$RangeBP, function(X) strsplit(X,'-')[[1]][1]))
  lociinfo$EndBP <- as.numeric(lapply(lociinfo$RangeBP, function(X) strsplit(X,'-')[[1]][2]))
  #head(lociinfo)
  # for each site in the VCF, find the associated locus based on the chromosome and position
  newfixes <- data.frame(CHROM=c(),POS=c(),ID=c(), REF=c(), ALT=c())
  for (i in 1:nrow(fix)){
    sub <- lociinfo[lociinfo$CHROM==fix$CHROM[i],] # subset to the chromosome level
    locus <- sub[which(sub$StartBP<fix$POS[i] & sub$EndBP>fix$POS[i]),] # find locus based on SNP position
    # copy over SNP info & create a new row for each template the SNP occurs in
    # (this approach adds one row for each locus the SNP occurs in, in case of overlapping loci)
    new <- fix[i,1:5]
    new <- as.data.frame(matrix(new, nrow=nrow(locus), ncol=5, byrow=T))
    names(new) <- c("CHROM","POS","ID","REF","ALT")
    new$POS <- as.numeric(new$POS)
    # reassign chrom
    new$CHROM <- gsub(">","",locus$ID)
    # adjust SNP position based on position within locus (instead of chrom)
    new$POS <- new$POS-locus$StartBP
    # add to DF
    newfixes <- rbind(newfixes, new) # add to new SNP df
  }#fix
  #head(newfixes)#check
  fix <- newfixes
}#if
## NOTE: The adjusted SNP position will be flanking regions + 1 
## (e.g., if 100-bp flanking regions were extracted with bedtools, POS=101)



#### IDENTIFY PRIMER BINDING LOCATIONS ####
# reset allowed primer binding locations for each template based on SNP position(s)
for (i in 1:nrow(templates)){
  id <- str_split(templates$ID[i], '>')[[1]][2] # grab locusID
  snps <- as.numeric(fix$POS[fix$CHROM==id])
  min_snp <- min(snps)
  max_snp <- max(snps)
  if (min_snp==Inf) min_snp <- templates$Sequence_Length[i]
  if (max_snp==-Inf) max_snp <- 1
  templates$Allowed_End_fw[i] <- min_snp-1
  templates$Allowed_Start_rev[i] <- max_snp+1
  templates$Allowed_End_rev[i] <- templates$Sequence_Length[i]
}#for i

# check all are populated correctly
templates$Allowed_Start_fw #these should all be 1
templates$Allowed_End_fw # these should be 1 bp before the first SNP in the locus
templates$Allowed_Start_rev # these should be 1 bp after the last SNP in the locus
templates$Allowed_End_rev # these should be the last bp in the locus

# extract allowed primer binding regions
for (i in 1:nrow(templates)){
  templates$Allowed_fw[i] <- substr(templates$Sequence[i], templates$Allowed_Start_fw[i], templates$Allowed_End_fw[i])
  templates$Allowed_rev[i] <- substr(templates$Sequence[i], templates$Allowed_Start_rev[i], templates$Allowed_End_rev[i])
}#for i

# check allowed binding regions
head(templates$Allowed_fw)
head(templates$Allowed_rev)



#### FILTERING TEMPLATES ####
## 1. remove any templates with insufficient binding space (due to SNPs at start/end of sequence)
templates$ID[-which(nchar(templates$Allowed_fw)>=minflankbp & nchar(templates$Allowed_rev)>=minflankbp)]
templates <- templates[which(nchar(templates$Allowed_fw)>=minflankbp & nchar(templates$Allowed_rev)>=minflankbp),]
nrow(templates)

## 2. remove any templates with more than 3 SNPs (these are likely paralogs, mRNA transcript variants, etc.)
names(counts)[counts>maxSNP]
templates <- templates[!templates$ID %in% paste0(">",names(counts)[counts>maxSNP]),]

# double check that there's only one line per locus
nrow(templates)
length(unique(templates$Header))


## 3. remove loci with high levels of correlation with other SNPs
## NOTE: This can also be done via LD-pruning, or ex post facto on optimized panel
# extract genotypes
#gt <- extract.gt(vcf, element="GT")
#for (i in 1:nrow(gt)){
#  gt[i,gt[i,]=="0/0"] <- 0
#  gt[i,gt[i,]=="0/1"] <- 1
#  gt[i,gt[i,]=="1/1"] <- 2
#}#for
#gt <- as.data.frame(gt)
#gt <- apply(gt, MARGIN=1, function(X) {as.numeric(X)})

# calc correlations, sum high correlations per SNP
#cors <- cor(gt, use="pairwise.complete")
#cors <- as.data.frame(cors)
#high_cors <- apply(cors, MARGIN=1, function(X) as.numeric(X>0.6))
##sum(high_cors,na.rm=T)-nrow(high_cors)
#totalcorrs <- rowSums(high_cors, na.rm=T)

# grab loci IDs w/ > 400 high corrs (correlated with ~10% of SNPs)
#sum(totalcorrs>400)
#highcorIDs <- rownames(cors)[totalcorrs>400]
#highcor_loci <- unlist(lapply(highcorIDs, function(X){ strsplit(X, ":")[[1]][1]}))

# remove
#templates <- templates[which(!templates$ID %in% paste0(">",highcor_loci)),]



#### EXPORT TEMPLATES ####
# double check that everything looks OK...
View(templates)

# extract microhaplotypes
tbl <- table(fix$CHROM,fix$POS) #rows=CHROM,cols=POS
counts <- apply(tbl, 1, function(X) sum(X>0)) #count SNPs (this accounts for possible duplicates at same POS)
microhaps_seq <- templates[which(templates$ID %in% paste0(">",microhaps$CHROM)),]
nrow(microhaps_seq)
#microhaps_seq <- microhaps_seq[which(nchar(microhaps_seq$Allowed_rev)>=18 | nchar(microhap_seq$Allowed_fw)>=18),]
#write.csv(microhaps_seq$ID, 'microhaplotypes.csv', row.names=FALSE)
# check that primer binding regions are correct for microhaps
microhaps <- fix[fix$CHROM%in%names(counts)[counts>1],]
microhaps[microhaps$CHROM==microhaps$CHROM[1],]
templates[templates$ID==paste0(">",microhaps$CHROM[1]),]

# export IDs, templates, and targets as CSV in primer3 format
target_len <- templates$Allowed_Start_rev - templates$Allowed_End_fw + 1
targets <- paste0(as.character(templates$Allowed_End_fw+1), ",", as.character(target_len))
templates_csv <- data.frame(SEQUENCE_ID=gsub(">","",templates$ID,'>',''),
                            SEQUENCE_TEMPLATE=templates$Sequence,
                            SEQUENCE_TARGET=targets)
write.csv(templates_csv, "GrayFox_primer3Templates.csv", row.names=FALSE)                         
