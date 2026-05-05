##--------RUN PRIMER-TREE FOR MULTIPLEX WORMHOLE OUTPUTS-------------##
# function:
runPrimerTree <- function(primers, organisms, 
                          fwd_adapter="tcgtcggcagcgtcagatgtgtataagagacag",
                          rev_adapter="gtctcgtgggctcggagatgtgtataagagacag",
                          all_combos=FALSE,
                          ...){
  # load packages
  library(primerTree)
  library(plyr)
  library(dplyr)
  
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
      if(all_combos | fwd==rev){
        print(paste("Running PRIMER-BLAST for",names[fwd],names[rev]))
        # loop through organisms
        for (o in organisms){
          skip_to_next <- FALSE
          tryCatch({
            search <- run_primer_search(forwards[fwd], reverses[rev], organism=o, ...)
            for (l in 1:length(search)){
              hits <- primerTree::parse_primer_hits(search[[l]])
              if (!is.null(hits)){
                hits$FWD <- paste0(names[fwd], ".FWD")
                hits$REV <- paste0(names[rev], ".REV")
                all_hits <- rbind(all_hits, hits)
              }#if
            }#for l
          },#try, 
          error = function(cond) {
              message(paste("PRIMER-BLAST failed for this combo:", names[fwd], names[rev]))
              message("Here's the original error message:")
              message(conditionMessage(cond))
              skip_to_next <- TRUE
          })#tryCatch
          # bypass errors
          if(skip_to_next) { next }
        }#for o
      }#if
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
    seq <- orimerTree::get_sequence(accession = all_hits$accession[row],
                        start = all_hits$product_start[row],
                        stop = all_hits$product_stop[row])
    all_hits$Sequence[row] <- paste(unlist(as.character(seq)), collapse="")
  }#for
  return(all_hits)
}#runPrimerTree





#---------------EXTRACT AMPLICON SEQUENCES FOR DESIGNED PRIMERS-------------#
extractPrimerInfo <- function(templates, filtprimers, finalprimers,
                                fwd_adapter="tcgtcggcagcgtcagatgtgtataagagacag",
                                rev_adapter="gtctcgtgggctcggagatgtgtataagagacag"){
  # First, load all datasets
  templates <- read.csv(templates)
  filtprimers <- read.csv(filtprimers)
  finalprimers <- read.csv(finalprimers)
  # ensure all names are character format
  templates$SEQUENCE_ID <- as.character(templates$SEQUENCE_ID)
  filtprimers$LocusID <- as.character(filtprimers$LocusID)
  filtprimers$PrimerID <- as.character(filtprimers$PrimerID)

  # extract primer pair ID from primer pair names formatted SeqID.#.FWD
  filtprimers$PairID <- unlist(lapply(filtprimers$PrimerID, function(X){
      x <- strsplit(X, "\\.")[[1]]
      return (paste(x[1],x[2], sep="."))
    }))#lapply

  # subset primer info for designed primers
  pairnames = unique(gsub(".REV","", gsub(".FWD", "", finalprimers$PrimerID)))
  primersub <- filtprimers[filtprimers$PairID %in% pairnames,]
  
  # copy over template info for each primer pair
  primersub$TemplateSeq <- templates$SEQUENCE_TEMPLATE[match(primersub$LocusID, 
                                                             templates$SEQUENCE_ID)]
  
  # grad amplicon info for each primer pair
  primersub$AmpliconSeq <- NA
  primersub$AmpliconTarget <- NA
  for (id in unique(primersub$LocusID)){
    # extract template sequence based on FWD-start + amplicon length
    row <- which(primersub$LocusID==id & endsWith(primersub$PrimerID,".FWD"))
    START <- primersub$StartBP[row]
    LEN <- primersub$AmpliconSize[row] - nchar(fwd_adapter) - nchar(rev_adapter)
    TEMPLATE <- primersub$TemplateSeq[row]
    primersub$AmpliconSeq[which(primersub$LocusID==id)] <- substr(TEMPLATE, START, START+LEN-1)
    # extract target BP- readjusted within amplicon
    TARGET <- templates$SEQUENCE_TARGET[match(id, templates$SEQUENCE_ID)]
    TARGET <- strsplit(TARGET, ",")[[1]]
    TARGET_START <- as.numeric(TARGET[1]); TARGET_LEN <- as.numeric(TARGET[2])
    NEW_TARGET <- paste(TARGET_START-START+1, TARGET_LEN, sep=",")
    primersub$AmpliconTarget[which(primersub$LocusID==id)] <- NEW_TARGET
  }#for pair

  return(primersub)
}#extractPrimerInfo


  


#-------------------PLOT PRIMER-BLAST RESULTS ALIGNED TO TARGETS----------------#
plotPrimerBlast <- function(primerblast, primerinfo, species="TARGET"){
  # load dependencies
  library(DECIPHER)
  library(Biostrings)
  # configure plotting environment to allow space for labels
  par(mar=c(1,1,1,16))
  # run new tree for each primer pair combination
  names <- unique(gsub(".FWD", "", primerblast$FWD))
  for (id in names){
    # subset results to those associated with this primer pair  
    sub <- primerblast[which(startsWith(primerblast$FWD,id) | startsWith(primerblast$REV,id)),]
    # extract amplicon sequence for pair
    amp <- primerinfo[primerinfo$PairID==id, c("LocusID","AmpliconSeq")][1,]#keep FWD record
    names(amp) <- c("accession","Sequence")
    amp$species <- species
    # merge primerblast and amplicon info for this primer pair
    sub <- plyr::rbind.fill(sub, amp)
    if(nrow(sub)>1){
      # convert sequences to DNAstringset
      dnaseqs <- tbl$Sequence
      names(dnaseqs) <- paste(tbl$species, tbl$accession, tbl)
      dnas <- Biostrings::DNAStringSet(dnaseqs, use.names=TRUE)
      # align (unique) sequences
      align <- DECIPHER::AlignSeqs(unique(dnas))
      # run maximum likelihood tree with ancestral state reconstruction
      mltree <- DECIPHER::TreeLine(align, method="ML", reconstruct=TRUE, type="dendrogram")
      # plot trees with number of state transitions
      plottree <- MapCharacters(mltree, labelEdges=TRUE)
      #plottree <- dendrapply(plottree, function(x){
      #  attr(x, "edgetext") <- paste(attr(x,"edgetext"),"\n")
      #})#dendrapply
      plot(plottree, edgePar=list(p.col=NA, p.border=NA, t.col="#55CC99", t.cex=0.8, t.font=2), 
           main=id, yaxt="n", horiz=TRUE)
      #attr(plottree[[1]], "change") #status changes on first branch left of (virtual) root
    }#if
  }#for
  # reset plot margins
  par(mar=c(5.1,4.1,4.1,2.1))
}#plotPrimerBlast  

