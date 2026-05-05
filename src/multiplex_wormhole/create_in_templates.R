# load dependencies
#library(openPrimeR)
library(vcfR)

# set input variables
setwd("/Users/maggiehallerud/Desktop/Gray_Fox_SNPs/Input_SNPs/")
fasta <- "MAF30All.CENSORed.fa"
invcf <- "MAF30AllGrayFoxPops.Filtered.recode.vcf"
type <- "refbased"
#type <- "denovo"
prefix <- "" # prefix of FASTA headers that needs to be removed for locusIDs to match CHROM in VCF
minflankbp <- 18 # minimum size of flanking regions (i.e., min primer binding site size)
maxSNP <- 3 # maximum number of SNPs allowed in a template
## TODO: Consider adjusting this to a SNP density, e.g. for variable template lengths


#### READ IN SEQUENCES & SNPs ####
# load template sequences as openPrimerR format
#templates <- openPrimeR::read_templates(fasta) #input fasta file
faraw <- read.table(fasta)$V1
templates <- data.frame(ID=faraw[startsWith(faraw,">")],
                        Header=faraw[startsWith(faraw,">")],
                        Sequence_Length=NA,
                        Allowed_Start_fw=1,
                        Allowed_End_fw=minflankbp,
                        Allowed_Start_rev=NA,
                        Allowed_End_rev=NA,
                        Sequence=faraw[!startsWith(faraw,">")])
templates$ID <- gsub(">","", templates$ID)
templates$ID <- gsub(prefix,"", templates$ID)
# set initial primer binding regions based on minimum flanking region
templates$Sequence_Length <- nchar(templates$Sequence)
templates$Allowed_Start_rev <- templates$Sequence_Length - minflankbp
templates$Allowed_End_rev <- templates$Sequence_Length
View(templates)#check

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
  head(lociinfo)
  
  # CHECK!!
  sum(!lociinfo$CHROM %in% fix$CHROM)
  
  # for each site in the VCF, find the associated locus based on the chromosome and position
  newfixes <- data.frame(CHROM=c(),POS=c(),ID=c(), REF=c(), ALT=c())
  excluded <- newfixes
  for (i in 1:nrow(fix)){
    sub <- lociinfo[lociinfo$CHROM==fix$CHROM[i],] # subset to the chromosome level
    if(nrow(sub)==0){
      warning(paste("WARNING: chromosome ID ",fix$CHROM[i]," not found in templates"))
    }#if
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
    # keep track of SNPs not found in templates
    if(nrow(new)==0){
      excluded <- rbind(excluded, new)
    }#if
  }#fix
  #head(newfixes)#check
  if(nrow(excluded)>0){
    print(paste("WARNING!",nrow(excluded),"SNPs","out of", nrow(fix), "were excluded. Check 'excluded' object for details."))
  }#if
  fix <- newfixes
}#if
## NOTE: The adjusted SNP position will be flanking regions + 1 
## (e.g., if 100-bp flanking regions were extracted with bedtools, POS=101)




#### IDENTIFY PRIMER BINDING LOCATIONS ####
# reset allowed primer binding locations for each template based on SNP position(s)
nosnps <- c()
for (i in 1:nrow(templates)){
  snps <- as.numeric(fix$POS[fix$CHROM==templates$ID[i]])
  if(length(snps)==0){
    nosnps <- templates$ID[i]
  }#if
  min_snp <- min(snps)
  max_snp <- max(snps)
  # remove SNP from targets if too close to start/end of sequence
  while (min_snp<templates$Allowed_End_fw[i]){
    templates$Allowed_Start_fw[i] <- min_snp+1
    snps <- snps[-which(snps==min_snp)]
    ifelse(length(snps)>0,  
           min_snp <- min(snps),
           min_snp <- Inf)
  }#while
  while (max_snp>templates$Allowed_Start_rev[i]){
    templates$Allowed_End_rev[i] <- max_snp-1
    snps <- snps[-which(snps==max_snp)]
    ifelse(length(snps)>0,
           max_snp <- max(snps),
           max_snp <- -Inf)
  }#while
  # define primer binding regions based on remaining snps
  if (!is.infinite(min_snp)) templates$Allowed_End_fw[i] <- min_snp-2
  if (!is.infinite(max_snp)) templates$Allowed_Start_rev[i] <- max_snp+2
}#for i

# raise warning if any sequences contained 0 SNPs
if(length(nosnps)>0){
  print(paste("WARNING:",length(snps), "FASTA sequences contained 0 SNPs in the VCF. See the 'nosnps' object for sequence IDs."))
}#if


# extract allowed primer binding regions
templates$Allowed_fw <- NA
templates$Allowed_rev <- NA
for (i in 1:nrow(templates)){
  templates$Allowed_fw[i] <- substr(templates$Sequence[i], templates$Allowed_Start_fw[i], templates$Allowed_End_fw[i])
  templates$Allowed_rev[i] <- substr(templates$Sequence[i], templates$Allowed_Start_rev[i], templates$Allowed_End_rev[i])
}#for i

# check allowed binding regions
sum(!nchar(templates$Allowed_fw)==(templates$Allowed_End_fw-templates$Allowed_Start_fw+1))#should be 0
# check all are populated correctly
View(templates)



# check whether all templates have SNPs and vice versa
unique(templates$ID[!templates$ID %in% fix$CHROM])#should be 0


#### FILTERING TEMPLATES ####
## 1. remove any templates with insufficient binding space (due to SNPs at start/end of sequence)
remove <- templates$ID[-which(nchar(templates$Allowed_fw)>=minflankbp & nchar(templates$Allowed_rev)>=minflankbp)]
remove
if (length(remove)>0){
  templates <- templates[which(nchar(templates$Allowed_fw)>=minflankbp & nchar(templates$Allowed_rev)>=minflankbp),]
}#if
nrow(templates)

## 2. remove any templates with more than 3 SNPs (these are likely paralogs, mRNA transcript variants, etc.)
counts <- table(fix$CHROM)
names(counts)[counts>maxSNP]
if (length(counts)>0){
  templates <- templates[!templates$ID %in% paste0(">",names(counts)[counts>maxSNP]),]
}#iff

# double check that there's only one line per locus
nrow(templates)==length(unique(templates$Header))#should be True


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
microhaps <- names(counts[counts>1])
microhaps_seq <- templates[which(templates$ID %in% microhaps),]
nrow(microhaps_seq)
# check that primer binding regions are correct for microhaps
#allowed_end_fw should be min POS-2, allowed_start_rev should be max POS+2
microhapvcf <- fix[fix$CHROM%in%names(counts)[counts>1],]
microhapvcf[microhapvcf$CHROM==microhapvcf$CHROM[1],]
templates[templates$ID==paste0(microhapvcf$CHROM[1]),]
# save microhaps separatel
target_len <- nchar(substr(microhaps_seq$Sequence, microhaps_seq$Allowed_End_fw, microhaps_seq$Allowed_Start_rev))
sum(target_len<1)
targets <- paste0(as.character(microhaps_seq$Allowed_End_fw+1), ",", as.character(target_len))
write.csv(data.frame(SEQUENCE_ID=microhaps_seq$ID,
                     SEQUENCE_TEMPLATE=microhaps_seq$Sequence,
                     SEQUENCE_TARGET=targets), 'Microhaplotype_templates.csv', row.names=FALSE)

# export IDs, templates, and targets as CSV in primer3 format
target_len <- templates$Allowed_Start_rev - templates$Allowed_End_fw + 1
sum(target_len<1)#check!
targets <- paste0(as.character(templates$Allowed_End_fw+1), ",", as.character(target_len))
templates_csv <- data.frame(SEQUENCE_ID=templates$ID,
                            SEQUENCE_TEMPLATE=templates$Sequence,
                            SEQUENCE_TARGET=targets)
write.csv(templates_csv, "SingletonSNP_Templates.csv", row.names=FALSE)
