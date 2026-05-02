##--------RUN PRIMER-TREE FOR MULTIPLEX WORMHOLE OUTPUTS-------------##
library(primerTree)
library(plyr)
library(dplyr)

# function:
runPrimerTree <- function(primers, organisms, 
                          fwd_adapter="tcgtcggcagcgtcagatgtgtataagagacag",
                          rev_adapter="gtctcgtgggctcggagatgtgtataagagacag",
                          ...){
  #---------READ INPUTS ----------#
  # load multiplex_wormhole output
  primers <- read.csv(primers)
  
  # remove Illumina Nextera adapters from primer sequences
  primers$Sequence <- gsub(fwd_adapter, "", primers$Sequence)
  primers$Sequence <- gsub(rev_adapter, "", primers$Sequence)
  
  # extract forward sequences
  forwards <- primers$Sequence[endsWith(primers$PrimerID, ".FWD")]
  reverses <- primers$Sequence[endsWith(primers$PrimerID, ".REV")]
  names <- gsub(".FWD", "", primers$PrimerID[endsWith(primers$PrimerID,".FWD")])
  
  #---------RUN PRIMER-BLAST SPECIFICITY CHECK FOR ALL PRIMERS-----------#
  ## NOTE: This will consider mispriming from all combinations of FWD/REV primer pairs!
  ## only one "organism" argument allowed, so run multiple times per "organism" group
  # set up function with global settings:
  run_primer_search <- function(fwd, rev, organism, ...){
    primerTree::primer_search(forward=fwd, reverse=rev, organism=organism, ...)
  }#run_primer_search
  all_hits <- data.frame()
  # loop through forwards
  for (fwd in 1:length(forwards)){
    # loop through all reverses
    for (rev in 1:length(reverses)){
      print(paste("Running PRIMER-BLAST for",names[fwd],names[rev]))
      # loop through organisms
      for (o in organisms){
        search <- run_primer_search(forwards[fwd], reverses[rev], organism=o, ...)
        for (l in 1:length(search)){
          hits <- primerTree::parse_primer_hits(search[[l]])
          if (!is.null(hits)){
            hits$FWD <- paste0(names[fwd], ".FWD")
            hits$REV <- paste0(names[rev], ".REV")
            all_hits <- rbind(all_hits, hits)
          }#if
        }#for l
      }#for o
    }#for rev
  }#for fwd
  
  #--------PULL TAXONOMY INFORMATION FOR ALL HITS--------------#
  print("Pulling taxonomy information from NCBI.....")
  taxa <- primerTree::get_taxonomy(all_hits$accession)
  all_hits <- merge(all_hits, taxa, by="accession")
  
  #--------PULL SEQUENCE INFORMATION FOR ALL HITS--------------#
  print("Retrieving sequences from NCBI.....")
  all_hits$Sequence <- NA
  for (row in 1:nrow(all_hits)){
    seq <- get_sequence(accession = all_hits$accession[row],
                        start = all_hits$product_start[row],
                        stop = all_hits$product_stop[row])
    all_hits$Sequence[row] <- paste(unlist(as.character(seq)), collapse="")
  }#for
  return(all_hits)
}#runPrimerTree


## The above function will PRIMER-BLAST all FWD/REV combinations of primers 
## against the set organisms. 
## Details on PRIMER-BLAST options found at bottom of script
## example usage:
# in this example, the target organism and relatives are checked (carnivora), as well as
# diet items (rodents and rabbits, aves [i.e. birds], plethodontidae, Ericaceae) and possible
# lab contaminants (human, bacteria)
primerblast <- runPrimerTree(primers="Final_Primers/MACA100_N20_Run01_primers.csv", #path to CSV of primers (PrimerID, Sequence)
                             organisms=c("Carnivora","Rodents and rabbits","Aves","Plethodontidae",
                                         "Ericacea","Homo sapiens","Bacteria"), #vector/list
                             fwd_adapter="tcgtcggcagcgtcagatgtgtataagagacag", #Nextera FWD adapter
                             rev_adapter="gtctcgtgggctcggagatgtgtataagagacag", #Nextera REV adapter
                             primer_specificity_database="core_nt", # use the core nucleotide database from NCBI
                             exclude_env="False", # don't exclude environmental/uncultured samples
                             MAX_TARGET_SIZE=1000 #max amplicon size
                             #add additional PRIMER-BLAST settings here if desired
)

#----------------PLOT TREE---------------#
library(DECIPHER)
library(Biostrings)
## NOTE: Consider adding target species amplicons into the primerblast CSV before proceeding

# configure plotting environment- these will be saved to a PDF since we're making lots
pdf("PRIMERBLAST_Trees.pdf")
par(mar=c(14,3,3,3))
# run new tree for each primer pair combination
names = unique(gsub(".REV","", gsub(".FWD", "", c(primerblast$FWD, primerblast$REV))))
for (name in names){
  # (for simplicity, only running actual primer pairs)
  sub <- primerblast[which(startsWith(primerblast$FWD,name) & startsWith(primerblast$REV,name)),]
  if(nrow(sub)>1){
    # convert sequences to DNAstringset
    dnaseqs <- sub$Sequence
    names(dnaseqs) <- paste(sub$species, sub$accession)
    dnas <- Biostrings::DNAStringSet(dnaseqs, use.names=TRUE)
    # align sequences
    align <- DECIPHER::AlignSeqs(dnas)
    # run maximum likelihood tree with ancestral state reconstruction
    mltree <- DECIPHER::TreeLine(align, method="ML", reconstruct=TRUE, type="dendrogram")
    # plot trees with number of state transitions
    plottree <- MapCharacters(mltree, labelEdges=TRUE)
    plot(plottree, edgePar=list(p.col=NA, p.border=NA, t.col="#55CC99", t.cex=0.7))
    attr(plottree[[1]], "change") #status changes on first branch left of (virtual) root
  }#if
}#for

dev.off()


### FULL PRIMER-BLAST OPTIONS WITH DEFAULTS SHOWN:
### PRIMER-BLAST WEBSITE: https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi
### To determine which option corresponds to which webpage setting, right-click
### on a website setting and select "Inspect element".
  #SEQFILE=NA,
  #PRIMER5_START=NA,
  #PRIMER5_END=NA,
  #PRIMER3_START=NA,
  #PRIMER3_END=NA,
  #PRIMER_LEFT_INPUT=NA,
  #PRIMER_RIGHT_INPUT=NA,
  #PRIMER_PRODUCT_MIN=70,
  #PRIMER_PRODUCT_MAX=1000,
  #PRIMER_NUM_RETURN=10,
  #PRIMER_MIN_TM=57.0,
  #PRIMER_OPT_TM=60.0,
  #PRIMER_MAX_TM=63.0,
  #PRIMER_MAX_DIFF_TM=3,
  #PREFER_3END=,
  #PRIMER_ON_SPLICE_SITE=0,
  #SPLICE_SITE_OVERLAP_5END=7,
  #SPLICE_SITE_OVERLAP_3END=4,
  #SPLICE_SITE_OVERLAP_3END_MAX=8,
  #SPAN_INTRON=,
  #MIN_INTRON_SIZE=1000,
  #MAX_INTRON_SIZE=1000000,
  #SEARCH_SPECIFIC_PRIMER="on",
  #SEARCHMODE=0,
  #PRIMER_SPECIFICITY_DATABASE="nt",
  #CUSTOMSEQFILE=NA,
  #EXCLUDE_XM=NA,
  #EXCLUDE_ENV=NA,
  #ORGANISM="Homo sapiens", 
  #ORGANISM2="rodents and rabbits 314147",
  #slctOrg=,
  #ENTREZ_QUERY=NA,
  #TOTAL_PRIMER_SPECIFICITY_MISMATCH=1,
  #PRIMER_3END_SPECIFICITY_MISMATCH=1,
  #MISMATCH_REGION_LENGTH=5,
  #TOTAL_MISMATCH_IGNORE=6,
  #MAX_TARGET_SIZE=1000,
  #ALLOW_TRANSCRIPT_VARIANTS=NA,
  #NEWWIN=NA,
  #SHOW_SVIEWER=on,
  #HITSIZE=50000,
  #EVALUE=30000,
  #WORD_SIZE=7,
  #MAX_CANDIDATE_PRIMER=500,
  #NUM_TARGETS=20,
  #NUM_TARGETS_WITH_PRIMERS=1000,
  #MAX_TARGET_PER_TEMPLATE=100
  #PRODUCT_MIN_TM=,
  #PRODUCT_OPT_TM=,
  #PRODUCT_MAX_TM=,
  #PRIMER_MIN_SIZE=15,
  #PRIMER_OPT_SIZE=20,
  #PRIMER_MAX_SIZE=25,
  #PRIMER_MIN_GC=20.0,
  #PRIMER_MAX_GC=80.0,
  #GC_CLAMP=0,
  #POLYX=5,
  #PRIMER_MAX_END_STABILITY=9,
  #PRIMER_MAX_END_GC=5,
  #TH_OLOGO_ALIGNMENT=,
  #TH_TEMPLATE_ALIGNMENT=,
  #PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00,
  #PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00,
  #PRIMER_MAX_SELF_ANY_TH=45.0,
  #PRIMER_MAX_SELF_END_TH=35.0,
  #PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0,
  #PRIMER_PAIR_MAX_COMPL_END_TH=35.0,
  #PRIMER_MAX_HAIRPIN_TH=24.0,
  #PRIMER_MAX_TEMPLATE_MISPRIMING=12.00,
  #PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00,
  #SELF_ANY=8.00,
  #SELF_END=3.00,
  #PRIMER_PAIR_MAX_COMPL_ANY=8.00,
  #PRIMER_PAIR_MAX_COMPL_END=3.00,
  #EXCLUDED_REGIONS=NA,
  #OVERLAP=NA,
  #OVERLAP_5END=7,
  #OVERLAP_3END=4,
  #MONO_CATIONS=50.0,
  #DIVA_CATIONS=1.5,
  #CON_DNTPS=0.6,
  #SALT_FORMULAR=1,
  #TM_METHOD=1,
  #CON_ANEAL_OLIGO=50.0,
  #NO_SNP=,
  #PRIMER_MISPRIMING_LIBRARY=AUTO,
  #LOW_COMPLEXITY_FILTER=on,
  #PICK_HYB_PROBE=,
  #PRIMER_INTERNAL_OLIGO_MIN_SIZE=18,
  #PRIMER_INTERNAL_OLIGO_OPT_SIZE=20,
  #PRIMER_INTERNAL_OLIGO_MAX_SIZE=27,
  #PRIMER_INTERNAL_OLIGO_MIN_TM=57.0,
  #PRIMER_INTERNAL_OLIGO_OPT_TM=60.0,
  #PRIMER_INTERNAL_OLIGO_MAX_TM=63.0,
  #PRIMER_INTERNAL_OLIGO_MIN_GC=20.0,
  #PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT=50,
  #PRIMER_INTERNAL_OLIGO_MAX_GC=80.0,
  #NA=NA,
  #NEWWIN=,
  #SHOW_SVIEWER=on