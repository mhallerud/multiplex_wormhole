#!/usr/bin/env Rscript

# input user-defined parameters
args=commandArgs(trailingOnly=TRUE)

if(length(args)<2) stop("Must provide input PrimerDimerReport.txt file & Delta G threshold")

IN = args[1] 
DELTAG_THRESHOLD = args[2]

#IN='/Users/maggiehallerud/Marten_Primer_Design/TEST/PrimerDimerReport_25jul2023.txt'
#DELTAG_THRESHOLD=-6

# load dependency
library(stringr)


# define prefix to use for outputs (match inputs to avoid confusion)
OUT=strsplit(basename(IN), split="[.]")[[1]][1]
OUTDIR=dirname(IN)

# read in file
infile=readLines(IN)

# make table with primer1, primer2, deltaG
dimers_txt=unique(infile[grep('kcal/mol', infile)])
primer1=unlist(lapply(1:length(dimers_txt), function(X) str_replace_all(str_split(dimers_txt[X], " ")[[1]][8], " ", "")))
primer2=unlist(lapply(1:length(dimers_txt), function(X) str_replace_all(str_replace(str_split(dimers_txt[X], " ")[[1]][10], ":", ""), " ","")))
deltaG=as.numeric(unlist(lapply(1:length(dimers_txt), function(X) str_split(dimers_txt[X]," ")[[1]][11])))
dimers <- data.frame(Primer1=primer1, Primer2=primer2, DeltaG=deltaG)
#write.table(dimers, paste0(paste(OUTDIR, OUT, sep="/"), ".tsv"), sep="\t", row.names=FALSE)
#nrow(dimers)

# only keep dimers below on deltaG threshold
filtered_dimers <- dimers[dimers$DeltaG<DELTAG_THRESHOLD,]
#write.table(filtered_dimers, paste0(paste(OUTDIR, OUT, sep="/"), "_Filtered.tsv"), sep="\t", row.names=FALSE)

# This function will ensure equivalent labels for e.g. Primer1xPrimer2 and Primer2xPrimer1 comparisons
# This doesn't matter for PrimerSuite, but could matter when primer interactions are duplicated in a report.
uniquePairLabels <- function(df, field1, field2){
  pairs <- sort(c(df[field1], df[field2]))
  return(paste(pairs[1],"x",pairs[2]))
}#uniquePairs

# export # dimers per primer
filtered_dimers$DeltaG <- 1 # convert delta G to binary (1=dimer present)
#remove self-dimers since we've already screened for these
filtered_dimers <- filtered_dimers[!filtered_dimers$Primer1==filtered_dimers$Primer2,] 
# check that there's only one interaction per primer pair...
#filtered_dimers$PrimersCompared <- apply(filtered_dimers, 1, function(X) uniquePairLabels(X, "Primer1","Primer2"))
#dimers_agg <- aggregate(DeltaG~PrimersCompared, data=filtered_dimers, FUN="sum")
#unique(dimers_agg$DeltaG)
# set up matrix to store wide format results
primers=sort(unique(c(filtered_dimers$Primer1, filtered_dimers$Primer2)))
primer_interactions=matrix(NA, nrow=length(primers), ncol=length(primers))
colnames(primer_interactions) <- primers
rownames(primer_interactions) <- primers
for (primer in primers){
  for (primer2 in primers){
    dimer_sub <- filtered_dimers[which(filtered_dimers$Primer1%in%c(primer,primer2) & filtered_dimers$Primer2%in%c(primer,primer2)),]
    primer_interactions[primer, primer2] <- nrow(dimer_sub)
    primer_interactions[primer2,primer] <- nrow(dimer_sub)
  }#primer2
}#primer
# export # interactions per primer
write.table(primer_interactions, paste0(paste(OUTDIR, OUT, sep="/"), "_PrimerInteractions_wide.csv"), sep=",", row.names=TRUE)

# aggregate interactions by primer pair
# make DF with pair IDs
pair_dimers <- filtered_dimers
pair_dimers$Pair1 <- unlist(lapply(pair_dimers$Primer1, function(X) str_split(str_split(X, "_FW")[[1]][1], "_REV")[[1]][1] ))
pair_dimers$Pair2 <- unlist(lapply(pair_dimers$Primer2, function(X) str_split(str_split(X, "_FW")[[1]][1], "_REV")[[1]][1] ))
# set up n_pairs x n_pairs matrix to store wide format results
pairs=sort(unique(c(pair_dimers$Pair1, pair_dimers$Pair2)))
primer_pair_interactions <- matrix(NA, nrow=length(pairs), ncol=length(pairs))
colnames(primer_pair_interactions) <- pairs
rownames(primer_pair_interactions) <- pairs
# loop through each combination of primer pairs and tally the number of unique primer-primer interactions (possible 0-4)
for (pair in pairs){
  for (pair2 in pairs){
    pair_sub <- pair_dimers[which(pair_dimers$Pair1%in%c(pair,pair2) & pair_dimers$Pair2%in%c(pair,pair2)),]
    primer_pair_interactions[pair,pair2] <- nrow(pair_sub) 
    primer_pair_interactions[pair2,pair] <- nrow(pair_sub)
  }#for pair2
}#for pair

# export wide format matrix
write.table(primer_pair_interactions, paste0(paste(OUTDIR, OUT, sep="/"), "_PrimerPairInteractions_wide.csv"), sep=",", row.names=TRUE)

# convert interactions by primer to long format (1 record per primer w/ # interactions)
primer_interaction_sums <- as.data.frame(rowSums(primer_interactions))
write.table(primer_interaction_sums, paste0(paste(OUTDIR, OUT, sep="/"), "_PrimerInteractions_sum.csv"), sep=",", row.names=TRUE)

# convert interactions by primer pair to long format (1 record per primer pair w. # interactions)
pair_interaction_sums <- as.data.frame(rowSums(primer_pair_interactions))
write.table(pair_interaction_sums, paste0(paste(OUTDIR, OUT, sep="/"), "_PrimerPairInteractions_sum.csv"), sep=",", row.names=TRUE)

# convert to binary primer pair interactions
binary_pair_interactions <- primer_pair_interactions
binary_pair_interactions[binary_pair_interactions>0] <- 1
write.table(binary_pair_interactions, paste0(paste(OUTDIR, OUT, sep="/"), "_PrimerPairInteractions_wide_binary.csv"), sep=",", row.names=TRUE)

# convert binary to long format
binary_pairs_sum <- as.data.frame(rowSums(binary_pair_interactions))
write.table(binary_pairs_sum, paste0(paste(OUTDIR, OUT, sep="/"), "_PrimerPairInteractions_sum_binary.csv"), sep=",", row.names=TRUE)
