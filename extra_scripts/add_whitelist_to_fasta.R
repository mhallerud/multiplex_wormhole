# set up env
library(seqinr)
library(stringr)


#### READ IN DATA
#setwd('/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/SNP Panel/Panel2_200initialpairs/')
# read in GTseq locusinfo - NOTE: I'm also including failed primer pairs in this so that the same loci aren't attempted
setwd('/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/SNP Panel/Panel2_200initialpairs/')
gtseq <- read.csv('Marten_LocusInfo_whitelist.csv') #input fasta file
templates <- read.csv("CoastalMartenTemplates_MAF10_CENSOR_Martesmartes_trimmed.csv")
fasta <- read_templates("SpecificityCheckTemplates_passed.fa")

# ensure that sequences are in same case
gtseq$ProbeSeq1 <- toupper(gtseq$ProbeSeq1)
gtseq$ProbeSeq2 <- toupper(gtseq$ProbeSeq2)
templates$SEQUENCE_TEMPLATE <- toupper(templates$SEQUENCE_TEMPLATE)

# check for each whitelist locus in the fasta- if present, remove template locus
all_matches <- c()
for (i in 1:nrow(gtseq)){
  probe1_matches <- which(grepl(gtseq$ProbeSeq1[i], templates$SEQUENCE_TEMPLATE, perl=TRUE))
  probe2_matches <- which(grepl(gtseq$ProbeSeq2[i], templates$SEQUENCE_TEMPLATE, perl=TRUE))
  if(length(probe1_matches)==0 & length(probe2_matches)==0){
    print(paste(gtseq$Locus.Name[i],": No matches"))
  }#if
  matches = unique(probe1_matches, probe2_matches)
  if (length(matches)>0){
    all_matches <- append(all_matches, matches)
  }#if
}#for

# check # of matches
length(all_matches)
length(unique(all_matches))

# return list of matching template names
ids <- paste0(">",templates$SEQUENCE_ID[unique(all_matches)])
fasta$locusID <- unlist(lapply(fasta$ID, function(X) str_split(X,"_")[[1]][1]))
fasta_sub <- fasta[-which(fasta$locusID %in% ids),]
nrow(fasta_sub)
nrow(fasta)

# add to fasta formatted list
fasta_sub$ID <- str_replace_all(fasta_sub$ID, ">", ">CLocus_")#rename names (not always necessary)
out <- c()
for (i in 1:nrow(fasta_sub)){
  out <- append(out, fasta_sub$ID[i])
  out <- append(out, fasta_sub$Sequence[i])
}#for

# add whitelist primers
whitelist <- read.csv("/Users/maggiehallerud/Marten_Primer_Design/Plate1_First55Pairs_Sep2023/IDT_Primer_Order_03Aug2023_Marten50loci.csv")
whitelist$loci <- unlist(lapply(whitelist$Name, function(X) str_split(X,"_")[[1]][2]))
failed <- c("131237", "322181", "322183", "402022", "402027", 
            "237309", "34762", "38234","465138","86842")
whitelist <- whitelist[-which(whitelist$loci%in%failed),]

whitelist$Name <- str_replace_all(whitelist$Name, "CLocus", ">MACA") #rename names
for (l in 1:nrow(whitelist)){
  out <- append(out, whitelist$Name[l])
  out <- append(out, tolower(whitelist$Sequence[l]))
}#for

# check
head(out)
tail(out)

# save to fasta
write(out, "SpecificityCheckTemplates_plusWhitelist.fa")

# write whitelist to separate fasta (for optimize_primers step)
out <- c()
for (l in 1:nrow(whitelist)){
  out <- append(out, whitelist$Name[l])
  out <- append(out, tolower(whitelist$Sequence[l]))
}#for

head(out)
length(out)/2

write(out, "Marten_Panel1_whitelist.fa")
