library(stringr)

# adding template info to primer pairs
pairs <- read.csv('/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/SNPpanel/Panel2_200initialpairs/multiplex_wormhole/4_OptimizedSets/Run0_150Loci_allLoci_dimers.csv',
                  colClasses=c("PrimerPairID"="character"))
templates <- read.csv('/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/SNPpanel/Panel2_200initialpairs/multiplex_wormhole/0_Inputs/CoastalMartenTemplates_MAF10_CENSOR_Martesmartes_trimmed.csv')
names(templates)
names(pairs)

pairs$SEQUENCE_ID <- unlist(lapply(pairs$PrimerPairID, function(x) strsplit(x, "\\.")[[1]][1]))

merged <- merge(pairs, templates, by="SEQUENCE_ID")
head(merged)

# add primer info
primers <- read.csv('/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/SNPpanel/Panel2_200initialpairs/multiplex_wormhole/4_OptimizedSets/Run0_150Loci_allLoci_primers.csv')
design <- read.csv('/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/SNPpanel/Panel2_200initialpairs/multiplex_wormhole/2_FilteredPrimers/SpecificityCheckTemplates_passed.csv')
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
  merged$FW[merged$PrimerPairID==pair] <- primersub$PrimerSeq[which(primersub$Direction.x=="FW")]
  merged$FWtemp[merged$PrimerPairID==pair] <- primersub$AnnealingTempC[which(primersub$Direction.x=="FW")]
  merged$FWpropbound[merged$PrimerPairID==pair] <- primersub$PropBound[which(primersub$Direction.x=="FW")]
  merged$REV[merged$PrimerPairID==pair] <- primersub$PrimerSeq[which(primersub$Direction=="REV")]
  merged$REVtemp[merged$PrimerPairID==pair] <- primersub$AnnealingTempC[which(primersub$Direction.x=="REV")]
  merged$REVpropbound[merged$PrimerPairID==pair] <- primersub$PropBound[which(primersub$Direction.x=="REV")]
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

View(merged)

write.csv(merged, '/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/SNPpanel/Panel2_200initialpairs/multiplex_wormhole/4_OptimizedSets/FullPrimerInfo_Run0_150allloci.csv')
