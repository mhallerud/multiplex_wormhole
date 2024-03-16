# multiplex_wormhole
Optimizing PCR primer design for multiplex amplicon sequencing.

**Important: This pipeline is still in development. Use at your own risk!**

# Dependencies
This pipeline relies on primer3 for primer design and MFEprimer for dimer calculations. To setup your files to run:
- Primer3 can be downloaded [here](https://github.com/primer3-org/primer3/releases). The path to primer3_core will need to be updated on line 18 of scripts/primer3.sh will need to be updated to reflect your local version. Primer3 version
- MFEprimer can be downloaded [here](https://www.mfeprimer.com/mfeprimer-3.1/#2-command-line-version). The path to MFEprimer-*-awd will need to be updated on lines 80-81 of multiplex_primer_design.sh to reflect your local version. MFEprimer Version 3.2.7 on Darwin was used during development. 

The pipeline was built and tested on MacOS with Python 3.9.13 in Spyder and bash 3.2.57(1)-release. The following (normally pre-installed) Python dependencies are required:
- pandas: install by running "pip install pandas" in the command line
General Python modules required (these normally come pre-installed): os, sys, csv, random, math, signal, gc, itertools


# Running pipeline
1. Create a CSV file with a row for each target and columns for locus ID, template sequence, and target position (following primer3 <start bp>,<length> format). For SNPs, the create_in_templates.R file takes a paired VCF and FASTA (FASTA loci match VCF CHROM field) and outputs the templates CSV. See example inputs in the [examples folder](https://github.com/mhallerud/multiplex_wormhole/examples).

2. Run multiplex_primer_design.py to design and filter primers for your targets:
   ./multiplex_primer_design.sh <TEMPLATE_CSV> <OUTPUT_DIRECTORY> <RESULTS_PREFIX>

   To specify primer design settings:
   - The primer3 settings used for primer design are found in primer3_settings/primer3_Base_NoSecondaryFilters.txt and primer3_settings/primer3_Broad_NoSecondaryFilters.txt. These settings can be edited directly, or you can create your own Broad and Base setting files and change the filepaths in the primer3_batch_design.sh script.
   - If you want to adjust PCR specifications, make sure to also change the parameters used for estimating dimer formation (lines 151-152 in multiplex_primer_design.py)

   To specify filtering criteria:
   - Lines 91-96 in the multiplex_primer_design.py file can be edited to specify different filtering parameters. 
   - Lines 151-152 in the multiplex_primer_design.py file can be edited to specify different delta G and score thresholds for considering dimers.


# Steps in pipeline
For more details on arguments and defaults for each function, go to the documentation for that function (multiplex_wormhole/docs).

0. Set up a folder structure for storing inputs and outputs
1. primer3_batch_design: Primers are designed for each template sequence using primer3, including predicting secondary structures (hairpins, homodimers, and heterodimers) within the primer pair.
2. filter_primers: Primer pairs are filtered to avoid likely secondary structures based on Gibbs free energy (deltaG) and annealing temperatures. 
3. check_primer_specificity: Specificity of primers are checked against all templates. Any primer pairs that aren't specific to one locus are discarded to avoid off-target amplification.
4. MFEprimer dimer: Primer dimers are predicted using MFEprimer, with one MFEprimer output including all primer dimers expected to form and a second output including only primer dimers forming on the 3' end of primers.
5. tabulate_MFEprimer_dimers: Primer dimer reports are translated into tables counting pairwise primer pair interactions and total interactions per primer pair. 
6. optimize_primers: A set of primers for "N" loci is selected that minimizes the number of negative interactions between primer pairs. An initial primer set is selected using a pseudo-greedy algorithm where the primer pairs with the cumulative lowest dimer load (across all loci provided) are selected, then adaptive simulated annealing is used to explore the optimization space around this initial primer set by randomly swapping out primer pairs and keeping improvements while allowing for 'mistakes' that may improve the primer set in the long run, and finally the best primer set found during adaptive simulated annealing is entered into a simple iterative improvement algorithm where the worst loci are swapped for better alternatives.
7. Specifically, the following steps are followed:
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

# Notes on the Optimization Process
Due to the random nature of the optimization algorithm, the returned primer set will only be an estimate of the ideal solution and will vary with each run. I recommend running the optimization 10-20 times and keeping the best set of these runs, unless a set with 0 dimers is found or the problem is simple (few desired # loci vs. many options).

This pipeline was created for designing multiplex PCR primers for SNP genotyping from reduced-representation sequencing data. The pipeline functions best when there are many potential targets relative to the number of desired loci (i.e., there are many alternatives that can be tested). The process will likely perform poorly on problems with few targets relative to the number of desired loci (e.g., 50 input loci with a desired SNP panel of 40 loci).

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
