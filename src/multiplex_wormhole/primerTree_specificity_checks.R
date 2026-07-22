##--------RUN PRIMER-TREE FOR MULTIPLEX WORMHOLE OUTPUTS-------------##
# function:
runPrimerTree <- function(primers, organisms, outcsv,
                          fwd_adapter="tcgtcggcagcgtcagatgtgtataagagacag",
                          rev_adapter="gtctcgtgggctcggagatgtgtataagagacag",
                          all_combos=FALSE,
                          THREADS=1, 
                          MAX_TARGET_SIZE=600, #max amplicon size is 600 bp
                          EXCLUDE_ENV="$0", #checks box for "exclude environmental samples"
                          PRIMER_SPECIFICITY_DATABASE="PRIMERDB/genome_selected_species", 
                          ...){
  # check inputs
  if (!file.exists(primers)) stop(paste0(primers," file could not be found!"))
  if (!all_combos%in%c(TRUE,FALSE)) stop("all combos must be TRUE or FALSE")
  if (!is.numeric(MAX_TARGET_SIZE)) stop("MAX_TARGET_SIZE is not numeric!")
  if (!is.numeric(THREADS)) stop("THREADS is not numeric!")
  if (!EXCLUDE_ENV%in%c("$0","")) stop("EXCLUDE_ENV must be '$0' or ''")
  # load packages
  library(primerTree)
  library(plyr)
  library(dplyr)
  #---------READ INPUTS ----------#
  # load multiplex_wormhole output
  primers <- read.csv(primers)
  # make everything lowercase
  primers$Sequence <- tolower(primers$Sequence)
  fwd_adapter <- tolower(fwd_adapter)
  rev_adapter <- tolower(rev_adapter)
  if (sum(grepl(fwd_adapter, primers$Sequence))<1) stop("fwd_adapter not found in primer sequences!")
  if (sum(grepl(rev_adapter, primers$Sequence))<1) stop("rev_adapter not found in primer sequences!")
  # remove Illumina Nextera adapters from primer sequences
  primers$Sequence <- gsub(fwd_adapter, "", primers$Sequence)
  primers$Sequence <- gsub(rev_adapter, "", primers$Sequence)
  # extract FWD, REV seqs & pair IDs
  orig_pairs <- data.frame(ID = gsub(".FWD", "", primers$PrimerID[endsWith(primers$PrimerID,".FWD")]),
                          FWDseq = primers$Sequence[endsWith(primers$PrimerID, ".FWD")],
                          REVseq = primers$Sequence[endsWith(primers$PrimerID, ".REV")])
  #---------RUN PRIMER-BLAST SPECIFICITY CHECK FOR ALL PRIMERS-----------#
  # set up df of primer pairs to check
  if(all_combos){
    pairs <- expand.grid(FWDid=orig_pairs$ID, REVid=orig_pairs$ID)
    pairs$FWDseq <- orig_pairs$FWDseq[match(pairs$FWDid, orig_pairs$ID)]
    pairs$REVseq <- orig_pairs$REVseq[match(pairs$REVid, orig_pairs$ID)]
  }else{
    pairs <- orig_pairs
    pairs$FWDid <- pairs$REVid <- pairs$ID
  }#ifelse
  #------SUB-FUNCTION TO RUN PRIMERTREE FOR EACH PRIMER PAIR ---------
  runPair <- function(pair, organisms, MAX_TARGET_SIZE, EXCLUDE_ENV, 
                      PRIMER_SPECIFICITY_DATABASE, ...){
    print(paste("Running PRIMER-BLAST for",pair$FWDid,pair$REVid))
    all_hits <- data.frame()
    # loop through organisms
    for (o in organisms){
        skip_to_next <- FALSE
        tryCatch({
            search <- primerTree::primer_search(pair$FWDseq, pair$REVseq,
                                                organism=o, 
                                                max_target_size=MAX_TARGET_SIZE,
                                                exclude_env=EXCLUDE_ENV,
                                                primer_specificity_database=PRIMER_SPECIFICITY_DATABASE,
                                                ...)
            Sys.sleep(0.1)#avoiding overloading NCBI
            # parse hits
            hits_list <- vector("list", length(search))
            if(length(search)>0){
              for (l in 1:length(search)){
                hits_list[[l]] <- primerTree::parse_primer_hits(search[[l]])
              }#for l
              new_hits <- do.call(rbind, Filter(Negate(is.null), hits_list))
              if (!is.null(new_hits)){
                new_hits$FWD <- paste0(pair$FWDid, ".FWD")
                new_hits$REV <- paste0(pair$REVid, ".REV")
                new_hits$FWDseq <- pair$FWDseq
                new_hits$REVseq <- pair$REVseq
                all_hits <- rbind(all_hits, new_hits)
              }#if
            }#if
          },#try, 
          error = function(cond) {
            message(paste("PRIMER-BLAST failed for this combo:", 
                          pair$FWDid, pair$REVid, "organism:", o))
            message("Here's the original error message:")
            message(conditionMessage(cond))
            skip_to_next <- TRUE
            NULL
          })#tryCatch
          # bypass errors in for loop
          if(skip_to_next) { next }
        }#for o in organisms
    return(all_hits)
    }#runPair
  #----------LOOP THROUGH PRIMER PAIRS & RUN PRIMERBLAST---------#
  # use multi-threading if available
  if (THREADS>1){
    # set up multi-threading
    library(parallel)
    library(doParallel)
    # check cores- adjust based on availability / needs
    nc = parallel::detectCores() #check available cores
    if (THREADS>nc) THREADS <- nc-1
    if (THREADS>nrow(pairs)) THREADS <- nrow(pairs)
    if (THREADS>10){
      print("WARNING: Running more than 10 threads is not recommended because it's likely to cause issues with the NCBI API!")
    }#if
    # initialize thread cluster
    cluster <- parallel::makeCluster(THREADS, outfile="")
    # register the cluster
    doParallel::registerDoParallel(cluster)
    # set up tmp dir
    tmp_dir <- file.path(tempdir(), "mw_primerblast")
    dir.create(tmp_dir, showWarnings=FALSE)
    # multi-thread primerTree across pairs
    all_hits <- foreach::foreach(i=1:nrow(pairs), .combine=rbind) %dopar% {
      result <- runPair(pairs[i,], organisms, MAX_TARGET_SIZE, EXCLUDE_ENV, PRIMER_SPECIFICITY_DATABASE, ...)
      # save progress...
      tmp_file <- paste0("blast__pair", i, "_tmp.csv")
      if(nrow(result)>0) write.csv(result, file.path(tmp_dir, tmp_file), row.names=FALSE)
      result
    }#foreach
    # close cluster
    parallel::stopCluster(cluster)
    # save checkpoint file, remove temps
    if (nrow(all_hits)>0) write.csv(all_hits, outcsv, row.names=FALSE)
    unlink(tmp_dir, recursive=TRUE)
  # if no multi-threading, loop through each primer pair instead,
  # saving progress as it runs
  }else{
      all_hits <- data.frame()
      for (i in 1:nrow(pairs)){
        new_hits <- runPair(pairs[i,], organisms, MAX_TARGET_SIZE, EXCLUDE_ENV, 
                            PRIMER_SPECIFICITY_DATABASE, ...)
        all_hits <- rbind(all_hits, new_hits)
        if(nrow(all_hits)>0) write.csv(all_hits, outcsv, row.names=FALSE)
      }#for
    }#ifelse
  # alternative cluster-friendly multi-threading
  #   library(future)
  #   library(future.apply)
  #   library(future.batchtools)
  #   # Detect environment and set plan accordingly
  #   if (Sys.getenv("SLURM_JOB_ID") != ""){
  #       #Sys.getenv("PBS_JOBID") != "" |
  #       #Sys.getenv("LSB_JOBID") != "" ){
  #   # Running on SLURM cluster
  #   plan(batchtools_slurm,
  #        resources = list(ncpus=THREADS, walltime="24:00:00", memory="8gb"))
  #   }else if(THREADS > 1){
  #     # Running locally with multiple cores
  #     plan(multisession, workers=THREADS)
  #   }else{
  #     # Single-threaded fallback
  #     plan(sequential)
  #   }#ifelse
  # # run pairs with parallel processing
  # all_hits <- do.call(rbind, future_lapply(1:nrow(pairs), function(i){
  #   runPair(pairs[i,], organisms, MAX_TARGET_SIZE, EXCLUDE_ENV,
  #             PRIMER_SPECIFICITY_DATABASE, ...)
  #   }))#do.call
  # # save progress
  # write.csv(all_hits, outcsv, row.names=FALSE)
  if (nrow(all_hits)>0){
    unique_accs <- unique(all_hits$accession)
    #--------PULL TAXONOMY INFORMATION FOR ALL HITS--------------#
    print("Pulling taxonomy information from NCBI.....")
    taxa <- primerTree::get_taxonomy(unique_accs)
    Sys.sleep(0.1)#avoiding overloading NCBI
    all_hits <- merge(all_hits, taxa, by="accession")
    write.csv(all_hits, outcsv, row.names=FALSE)
    #--------PULL SEQUENCE INFORMATION FOR ALL HITS--------------#
    print("Retrieving sequences from NCBI.....")
    # run with rentrez::entrez_fetch
    seqs <- primerTree::get_sequences(accession = all_hits$accession,
                                      start = all_hits$product_start,
                                      stop = all_hits$product_stop)
    Sys.sleep(0.1)#avoiding overloading NCBI
    seqs <- lapply(as.character(seqs), function(x) paste(unlist(x), collapse=""))
    all_hits$Sequence <- as.character(seqs[match(all_hits$accession, names(seqs))])
  }#if
  #-----SAVE FINAL OUTPUT-------#
  write.csv(all_hits, outcsv, row.names=FALSE)
  #return(all_hits)
}#runPrimerTree





#---------------EXTRACT AMPLICON SEQUENCES FOR DESIGNED PRIMERS-------------#
extractPrimerInfo <- function(templates, filtprimers, finalprimers,
                              fwd_adapter="tcgtcggcagcgtcagatgtgtataagagacag",
                              rev_adapter="gtctcgtgggctcggagatgtgtataagagacag"){
  # check inputs
  if (!file.exists(templates)) stop(paste0(templates," file could not be found!"))
  if (!file.exists(filtprimers)) stop(paste0(filtprimers," file could not be found!"))
  if (!file.exists(finalprimers)) stop(paste0(finalprimers," file could not be found!"))
  # First, load all datasets
  templates <- read.csv(templates)
  filtprimers <- read.csv(filtprimers)
  finalprimers <- read.csv(finalprimers)
  if (!sum(grepl(fwd_adapter, tolower(finalprimers$SEQUENCE)))>0) stop("fwd_adapter not found in finalprimer sequences!")
  if (!sum(grepl(rev_adapter, tolower(finalprimers$SEQUENCE)))>0) stop("rev_adapter not found in finalprimer sequences!")
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
plotAmpliconTrees <- function(primerblast, primerinfo=NA, species="TARGET", dG=0, 
                              dG_end=NA, MAX_AMPLICON_SIZE=500, THREADS=1, ...){
  #-------READ INPUTS----------#
  # check for proper inputs...
  if(!is.data.frame(primerblast)) stop("primerblast must be a dataframe!")
  if(!is.na(primerinfo)) if(!is.data.frame(primerinfo)) stop("primerinfo must be a dataframe!")
  if(!is.na(dG)) if(!is.numeric(dG)) stop("dG must be numeric!")
  if(!is.na(dG_end)) if(!is.numeric(dG_end)) stop("dG must be numeric!")
  if(!is.character(species)) stop("species must be characterformat!")
  if(!is.numeric(MAX_AMPLICON_SIZE)) stop("MAX_AMPLICON_SIZE must be numeric!")
  if(!is.numeric(THREADS)) stop("THREADS must be numeric!")
  # load dependencies
  library(DECIPHER)
  library(Biostrings)
  print(paste("# Off-target sequences:", nrow(primerblast)))
  # filter primerblast outputs based on thermodynamics
  if(!is.na(dG)){
    primerblast <- primerblast[which(primerblast$dG_FWD<dG & primerblast$dG_REV<dG),]
  }#ifdG
  if(!is.na(dG_end)){
    primerblast <- primerblast[which(primerblast$dG_FWD_END<dG_end & primerblast$dG_REV_END<dG_end),]
  }#ifdG_end
  print(paste("# Off-target sequences after delta G filtering:", nrow(primerblast)))
  # filter based on amplicon size
  if(!is.na(MAX_AMPLICON_SIZE)){
    primerblast <- primerblast[which(primerblast$product_length<=MAX_AMPLICON_SIZE),]
  }#if
  print(paste("# Off-target sequences after amplicon size filtering:", nrow(primerblast)))
  
  #-------SUB-FUNCTION TO PLOTTREE FOR SINGLE PRIMER PAIR-------#
  plotTree <- function(id, ...){
    # configure plotting environment to allow space for labels
    # subset results to those associated with this primer pair  
    sub <- primerblast[which(startsWith(primerblast$FWD,id) | startsWith(primerblast$REV,id)),]
    # extract amplicon sequence for pair
    if(!is.na(primerinfo) && is.data.frame(primerinfo)){
      amp <- primerinfo[primerinfo$PairID==id, c("LocusID","AmpliconSeq")][1,]#keep FWD record
      names(amp) <- c("accession","Sequence")
      amp$species <- species
      # merge primerblast and amplicon info for this primer pair
      sub <- plyr::rbind.fill(sub, amp)
    }#if
    if(nrow(sub)>1){
      # convert sequences to DNAstringset
      dnaseqs <- sub$Sequence
      names(dnaseqs) <- paste(sub$species, sub$accession)
      dnas <- Biostrings::DNAStringSet(dnaseqs, use.names=TRUE)
      if(length(unique(dnas))>1){
        # align (unique) sequences
        align <- DECIPHER::AlignSeqs(unique(dnas))
        # run maximum likelihood tree with ancestral state reconstruction
        mltree <- DECIPHER::TreeLine(align, method="ML", reconstruct=TRUE, type="dendrogram")
        # plot trees with number of state transitions
        plottree <- DECIPHER::MapCharacters(mltree, labelEdges=TRUE)
        #plottree <- dendrapply(plottree, function(x){
        #  attr(x, "edgetext") <- paste(attr(x,"edgetext"),"\n")
        #})#dendrapply
        plot(plottree, edgePar=list(p.col=NA, p.border=NA, t.col="#55CC99", t.cex=0.8, t.font=2), 
             main=id, yaxt="n", horiz=TRUE, ...)
        #attr(plottree[[1]], "change") #status changes on first branch left of (virtual) root
        return(plottree)
      }else{
        return(NULL)
      }#ifelse
    }else{
      print("No sequences passed filtering to plot!")
      return(NULL)
    }#ifelse
  }#plotTree
  #--------LOOP THROUGH PRIMER PAIRS TO PLOT TREES----------#
  names <- unique(gsub(".FWD", "", primerblast$FWD))
  if (length(names)>0){
    if(THREADS>1){
      library(parallel)
      library(doParallel)
      # check cores- adjust based on availability / needs
      nc = parallel::detectCores() #check available cores
      if (THREADS>nc) THREADS <- nc-1
      if (THREADS>length(names)) THREADS <- length(names)
      # set up cluster
      cluster <- parallel::makeCluster(THREADS, outfile="")
      doParallel::registerDoParallel(cluster)
      # multi-thread plotting for each primer pair
      plots <- foreach::foreach(i=1:length(names)) %dopar% {
        plotTree(names[i])
      }#foreach
      # now plot...
      for (i in 1:length(plots)){
        if(!is.null(plots[[i]])){
          plot(plots[[i]], edgePar=list(p.col=NA, p.border=NA, t.col="#55CC99", t.cex=0.8, t.font=2), 
             main=names[i], yaxt="n", horiz=TRUE, ...)
        }else{
          print(paste("No plot for ", names[i]))
        }#ifelse
      }#for
      # close cluster
      parallel::stopCluster(cluster)
    }else{
        for (id in names){
          plotTree(id, ...)
        }#for
      }#ifelse
  }#else
}#plotAmpliconTrees





plotMismatches <- function(primerblast, group, title=element_blank()){
  if(!is.data.frame(primerblast)) stop("primerblast must be a data frame!")
  if(!group%in%names(primerblast)) stop("'group' is not a field in primerblast!")
  library(ggplot2)
  # calc overall mismatches
  primerblast$mismatch_forward <- as.numeric(primerblast$mismatch_forward)
  primerblast$mismatch_reverse <- as.numeric(primerblast$mismatch_reverse)
  primerblast$mismatches <- sapply(1:nrow(primerblast), 
                                   function(x) (primerblast[x,"mismatch_forward"]+
                                                  primerblast[x,"mismatch_reverse"]) / 2)
  # find "worst" off-target for each group
  group <- primerblast[,group]
  names <- sapply(unique(primerblast$FWD), function(x) strsplit(x,"\\.")[[1]][1])
  mm_out <- expand.grid(Primer=names, Group=unique(group))
  mm_out$Primer <- as.character(mm_out$Primer)
  # mm_out$MinMean <- sapply(1:nrow(mm_out), function(x){
  #   sub <- primerblast$mismatches[which(startsWith(primerblast$FWD, mm_out$Primer[x])
  #                                       & group==mm_out$Group[x])]
  #   ifelse(length(sub)>0, return(min(sub)), return(NA))
  # })#sapply
  mm_out$MinFWD <- sapply(1:nrow(mm_out), function(x){
    sub <- primerblast$mismatch_forward[which(startsWith(primerblast$FWD, mm_out$Primer[x])
                                        & group==mm_out$Group[x])]
    ifelse(length(sub)>0, return(min(sub)), return(NA))
  })#sapply
  mm_out$MinREV <- sapply(1:nrow(mm_out), function(x){
    sub <- primerblast$mismatch_reverse[which(startsWith(primerblast$FWD, mm_out$Primer[x])
                                        & group==mm_out$Group[x])]
    ifelse(length(sub)>0, return(min(sub)), return(NA))
  })#sapply
  # reshape FWD/REV mismatches into rows
  mm_out_fwd <- mm_out[,c("Primer","Group","MinFWD")]
  mm_out_rev <- mm_out[,c("Primer","Group","MinREV")]
  mm_out_fwd$Primer <- paste0(mm_out_fwd$Primer, "-FWD")
  mm_out_rev$Primer <- paste0(mm_out_rev$Primer, "-REV")
  names(mm_out_fwd)[3] <- names(mm_out_rev)[3] <- "Mismatches"
  mm_out <- rbind(mm_out_fwd, mm_out_rev)
  mm_out$Primer <- factor(mm_out$Primer, levels=sort(unique(mm_out$Primer)))
  # plot mean mismatches per primer pair x group
  plot <- ggplot2::ggplot(mm_out)+
    geom_tile(aes(x=Primer, y=Group, fill=Mismatches),col='black',lwd=0.2)+
    theme(axis.text.x=element_text(angle=90, hjust=1))+
    scale_fill_gradient2(low="darkred", mid="orange", high="yellow", na.value="white",
                         midpoint=2.5)+#median(mm_out$Mismatches,na.rm=T))+
    xlab(element_blank())+ylab(element_blank())+
    ggtitle(title)
  print(plot)
}#plotMismatches





# function to convert FASTA to CSV
FASTA2CSV <- function(infa, outcsv){
  lines <- readLines(infa)
  names <- gsub("^>", "", lines[startsWith(lines, ">")])
  seqs <- lines[!startsWith(lines, ">")]
  out <- data.frame(PrimerID=names, Sequence=seqs)
  write.csv(out, outcsv)
}#FAST2CSV


