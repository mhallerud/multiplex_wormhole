# multiplex_wormhole
Optimizing PCR primer design for multiplex amplicon sequencing by minimizing predicted dimer loads. The target audience is GT-seq panel development, however the process is transferable to any application where multiple distinct amplicons are being targeted for multiplex PCR.

**Important: This pipeline is still in development. Use at your own risk!**

# Dependencies
This pipeline relies on primer3 for primer design and MFEprimer for dimer calculations. To setup your files to run:
- Primer3 can be downloaded [here](https://github.com/primer3-org/primer3/releases). The path to primer3_core will need to be updated on line 18 of scripts/primer3.sh will need to be updated to reflect your local version. Primer3 version
- MFEprimer can be downloaded [here](https://www.mfeprimer.com/mfeprimer-3.1/#2-command-line-version). The path to MFEprimer-*-awd will need to be updated on lines 80-81 of multiplex_primer_design.sh to reflect your local version. MFEprimer Version 3.2.7 on Darwin was used during development. 

The pipeline was built and tested on MacOS with Python 3.9.13 and bash 3.2.57(1)-release.

# Running the pipeline
1. Create a CSV file with a row for each target locus and columns for locus ID, DNA sequence, and target position (following the primer3 <start bp>,<length> format). For SNPs, the create_in_templates.R file takes the input VCF and matching FASTA (i.e., FASTA sequence IDs match VCF CHROM IDs) and outputs the templates CSV. See example inputs in the [examples folder](https://github.com/mhallerud/multiplex_wormhole/examples). 

2. If you have a pre-existing primer set that you are looking to expand, create a FASTA file containing the adapter-ligated primer sequences.
  
3. Run the multiplex_primer_design.py script to design primers for your loci, filter primers based on thermodynamic settings, and optimize a primer set for multiplex PCR.

   To specify primer design settings:
   - The primer3 settings used for primer design are found in primer3_settings/primer3_Base_NoSecondaryFilters.txt and primer3_settings/primer3_Broad_NoSecondaryFilters.txt. These settings can be edited directly, or you can create your own Broad and Base setting files and change the filepaths in the primer3_batch_design.py script.
   - The base primer3 settings add Nextera Illumina adapters to all primers by default. If your lab workflow does not involve using adapter-ligated primers or if you are using a different set of adapters, adjust the 'SEQUENCE_OVERHANG_LEFT' and 'SEQUENCE_OVERHANG_RIGHT' settings in the primer3 settings files.

   To specify primer filtering criteria:
   - Lines 88-93 in the multiplex_primer_design.py file can be edited to specify different filtering parameters. 
   
   To specify dimer prediction criteria:
   - Lines 150-151 in the multiplex_primer_design.py file can be edited to specify different delta G and score thresholds for reporting dimers.

# Optimization Workflow & Information on Default Settings
1. Primers are designed for each template sequence using primer3, including calculating secondary structures using Illumina i5 and i7 Nextera adapters as overhangs.
  Defaults:
   - Annealing Temp: 52 C
   - Amplicon size: 70-120 bp
   - Primer size: 18-26 bp (optimal: 20)
   - Primer Tm: 57-63 Celsius (optimal: 60)
   - Primer GC content: 30-70% (optimal: 50)
   - GC Clamp: 1
   - Max End GC: 4
   - Max Poly X: 4
   - dNTP concentration: 0.25 mM
   - template concentration: 50 nM
   - divalent cation concentration: 3.8 mM
    - monovalent cation concentration: 50 mM
  If primers can't be found for the above settings, the following constraints are relaxed:
    - Primer GC content: 20-80%
    - GC Clamp: 0
    - Max End GC: 5
    - Max Poly X: 5
2. Primer pairs are filtered to avoid secondary structures. By default, primer pairs are discarded if the following criteria are met:
    - secondary structure Tm > 45
        AND
    - delta G < threshold
        - hairpins: -2 kcal/mol
        - homodimers or heterodimers at primer ends: -5 kcal/mol
        - homodimers or heterodimers not at primer ends: -10 kcal/mol
3. Specificity of primers are checked against all templates. Any primer pairs that aren't specific to one locus are discarded.
   Note: I also recommend checking amplicon sequences against a genome or larger dataset to ensure that you aren't amplifying multiple genomic regions (unless this is the goal). This can be done using short-read sequence alignment software such as bwa mem or by BLAST-ing sequences against the genome(s) of interest.
5. Primer dimer formation, that is chemical binding between primers, is the largest issue with designing large multiplex PCR primer sets- especially when working with degraded samples with lower quality DNA. Primer dimers are predicted using MFEprimer, with one MFEprimer output including all primer dimers expected to form and a second output including only primer dimers forming on the 3' end (the least stable end) of primers. The following defaults are used for reporting predicted dimer formation:
    - delta G threshold for structures at the 3' end: -3 kcal/mol
    - delta G threshold for any other structures: -6 kcal/mol
    - score limit: 3 (score calculated with +1 for each bp match and -1 for each bp mismatch)
    - max mismatches in dimer: 40 bp
    - dNTP concentration: 0.25 mM
    - oligo concentration: 50 nM
    - divalent cation concentration: 3.8 mM
    - monovalent cation concentration: 50 mM
6. Primer dimer reports are translated into the following tables:
    - A N_PRIMERS x N_PRIMERS table which includes the total number of primer dimers estimated for all pairwise combinations of primer pairs.
    - A N_PRIMERS x 2 table which includes the total number of primer dimers contributed by each primer pair.
    - A N_PRIMERS X N_PRIMERS table which specifies pairwise primer pair interactions in binary (0 = no secondary structures between primer 1 and primer 2, 1 = at least 1 secondary structure between primer 1 and primer 2).
    - A N_PRIMERS X 2 table which inclues the total number of primer pairs with which each primer pair interacts.
   These tables don't include interactions between primers within the same pair (which were already handled in step 2).
7. A set of primers for "N" loci is selected that minimizes the number of dimers between selected primer pairs. This algorithm was largely inspired by SADDLE and Marxan Conservation Planning software. Specifically, the following steps are followed:
     - The best primer pairs (i.e., the primer pairs with the lowest dimer load) are identified for each locus.
           NOTE: Dimer load can be calculated either as the total number of primer dimers that loci contributes, or as the number of binary primer pair interactions (1 = primer pair has at least one primer dimer, 0 = primer pair has no interactions) the loci contributes. By default, the algorithm considers the binary primer pair interactions. To optimize based on total primer dimers instead, use the *_PrimerPair_interactions_sum.csv and *_PrimerPair_interactions_wide.csv as the DIMER_SUMS and DIMER_TABLE inputs. 
     - An initial set of primer pairs is selected by choosing the best (i.e., lowest dimer load) "N" primer pairs from the locus-specific bests.
     - The algorithm then iteratively improves this initial set of primer pairs by switching the worst (i.e., highest dimer load) primer pair with a random alternative (either for the same locus, or from a locus that is not currently included in the set).
         - If the change resulted in an improvement (i.e., fewer TOTAL dimer load for the set), then this new set is used as the starting point for the next iteration.
         - If the change did not improve the set, then other alternatives for this worst pair are tested until a replacement that improves the set is found. If no alternatives improve the set, then this primer pair is blacklisted (i.e., kept in the set but a message will be returned suggesting to possibly remove this pair as an option) and the algorithm moves on to the next worst primer pair.
      - The algorithm stops when one of the following conditions are met:
         - The total primer dimer load is 0, and therefore there is no longer a need to optimize.
         - No further improvements can be made to the set (i.e., all remaining primer pairs have been blacklisted- no alternatives to these will result in improvements).
         - A user-defined number of iterations is reached.
Because the randomness in this process leads to slightly different outputs in each run, we recommend running 10-100 times (depending on the difficulty of your problem) and using the run with the lowest overall dimer load.

# Notes on the Optimization Process
This pipeline was created for designing multiplex PCR primers for SNP genotyping from reduced-representation sequencing data. The pipeline functions best when there are many potential targets relative to the number of desired loci (i.e., there are many alternatives that can be tested). Specifically, we have been using ~200 potential loci to design primer sets for GT-seq with a target of ~50 loci. The process will likely perform poorly on problems with few targets relative to the number of desired loci (e.g., 50 input loci with a desired assay of 40 loci).

# Using other dimer calculation tools
The pipeline is built to use MFEprimer dimer to calculate dimer formation, however the optimization process will accept any input tables as long as the 2 input tables specify 1) pairwise dimer loads between primer pairs and 2) the total dimer load per primer pair, with primer pair IDs matching between the input templates and both tables. See example inputs in the [examples folder](https://github.com/mhallerud/multiplex_wormhole/examples).

One alternative to MFEprimer dimer is Primer Suite's Primer Dimer tool. To use this tool to calculate dimers:
1. Copy the *_specificityCheck_passed.csv file into [primer-dimer.com](https://primer-dimer.com), select 'Multiplex Analysis' and 'Dimer Structure Report', and click Submit. Depending on how many loci you provided, this step may take awhile (~20 minutes for primers for 50 loci, ~xx hours for 1200 loci).

2. Run scripts/translate_primerSuite_report.R on the PrimerDimerReport.txt file downloaded from primer-dimer.com to convert this to a CSV. This can either be run directly in R or via the command line following the syntax:
   Rscript --vanilla scripts/translate_primerSuite_report.R <PATH_TO_PrimerDimerReport.txt> <DELTA_G_THRESHOLD>
The delta G threshold specified will filter out any primer dimers with delta G above this value.

Another alternative is [ThermoFisher's Multiplex Primer Design](https://www.thermofisher.com/us/en/home/brands/thermo-scientific/molecular-biology/molecular-biology-learning-center/molecular-biology-resource-library/thermo-scientific-web-tools/multiple-primer-analyzer.html), however multiplex_wormhole currently has no support for translating these outputs into tables.


## Citations
Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M and Rozen SG. 2012. Primer3--new capabilities and interfaces. *Nucleic Acids Research* 40(15):e115. [http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3424584/](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3424584/)

Wang K, Li H, Xu Y, Shao Q, Yi J, Wang R, Cai W, Hang X, Zhang C, Cai H, and Qu W. 2019. MFEprimer-3.0: quality control for PCR primers.
*Nucleic Acids Research* 47(W1): W610–W613. [https://doi.org/10.1093/nar/gkz351](https://doi.org/10.1093/nar/gkz351)

Xie NG, Wang MX, Song P, Mao S, Wang Y, Yang Y, Luo J, Ren S, and Zhang DY. 2022. Designing highly multiplex PCR primer sets with Simulated Annealing Design using Dimer Likelihood Estimation (SADDLE). *Nature Communications* 13(1): 1881. doi: [10.1038/s41467-022-29500-4](https://doi.org/10.1038/s41467-022-29500-4). 
