# multiplex_wormhole
*In silico* optimization for multiplex PCR assays by minimizing predicted primer dimer loads. The target audience is GT-seq panel developers, however the process is transferable to any application where many amplicons are targeted in multiplex PCR.

multiplex_wormhole was created for the purpose of optimizing large multiplex PCR panels that are applied to noninvasive wildlife genetic samples for the purposes of SNP-based individual identification and relatedness analysis.  Default settings are therefore tailored towards amplifying degraded DNA from low quality samples. 

**Important: This pipeline is still in development. Use at your own risk!**
**Please report any problems or potential enhancements in the GitHub Issues page.**


# Dependencies
- Primer3 is used for primer design and can be downloaded [here](https://github.com/primer3-org/primer3/releases). Keep track of where this downloads- the path to primer3_core will need to be added at line 55 in multiplex_primer_design.py.
- MFEprimer is used for dimer calculations and can be downloaded [here](https://www.mfeprimer.com/mfeprimer-3.1/#2-command-line-version). The path to MFEprimer-*-awd will need to be updated on line 54 of multiplex_primer_design.py. MFEprimer Version 3.2.7 on Darwin was used during development. 

The following (normally pre-installed) Python dependencies are required:
- pandas: install by running "pip install pandas" in the command line. Developed with pandas version 1.4.4
General Python modules required (these come pre-installed with most Python installations): os, sys, csv, random, math, signal, gc, itertools, shutil, mathplotlib

multiplex_wormhole was built and tested on MacOS with Python 3.9.13 in Spyder and bash 3.2.57(1)-release. 


# Quick-Start Guide
1. Create a CSV file with a row for each target and columns for locus ID, template sequence, and target position following primer3 format: <start bp>,<length>. The [create_in_templates](scripts/create_in_templates.R) R script can be used to create this file using a VCF and FASTA as inputs. See example inputs in the [examples folder](https://github.com/mhallerud/multiplex_wormhole/examples). Each DNA sequence in the FASTA file will be treated as a unique target for PCR amplification (i.e., sequences should each include one "target" region). The R script handles two input types:

* *de novo*: Matches formatting output of *de novo* Stacks pipeline where the FASTA was created using `populations --fasta-loci`. Assumes the VCF 'CHROM' field matches the FASTA sequence headers with a prefix such as "CLocus_".

* *reference-aligned: Assumes the FASTA sequence header is in the format `CHROM`:`startBP`-`endBP` where `CHROM` matches the reference sequence denoted as 'CHROM' field in the VCF, `startBP` is the position in the reference sequence where the FASTA sequence starts, and `endBP` is the position in the reference sequence where the FASTA sequence ends. This format can be extracted from reference-aligned SNP data in VCF format (e.g., output from Stacks, bcftools mpileup, etc.) using the following commands with REF.FASTA as the reference genome:

```
# NOTE: Flanking regions of 100-bp are used here as this is roughly the amplicon length targeted and provides good primer binding space on either side of the target SNP, but this can be adjusted if longer or shorter amplicons are desired. 
# thin to 1 SNP per 100-bp to avoid overlapping template sequences
# this step is not necessary if linkage disequilibrium pruning has already been performed on the input VCF
vcftools --vcf <IN.VCF> --thin 100 --out Thinned100bp --recode

# convert SNPs into bed format 
bcftools query -f '%CHROM\t%POS\n' Thinned100bp.recode.vcf | awk '{print $1"\t"$2-1"\t"$2}' | awk '{print $0"\t"$1":"$3}' > Thinned100bp.bed

# make tab-delimited file with reference sequence IDs in the first field and sequence lengths in the second field
seqkit fx2tab --length --name <REF.FASTA> chr_lengths.txt 

# grab 100-bp flanking regions around thinned SNPs
bedtools slop -b 100 -i Thinned100bp.bed -g chr_lengths.txt > Thinned100bp_Flanking100bp.bed

# make fasta of these
bedtools getfasta -fi <REF.FASTA> -bed Thinned100bp_Flanking100bp.bed -fo <OUT.FA>
```

2. Run multiplex_primer_design.py by changing the filepaths on lines 43, 54-61. Default parameters are provided but can be adjusted within this script.

***This pipeline has very specific file requirements and naming conventions. Problems can be avoided by running the full pipeline from start to finish for each new dataset.***


# Steps in the pipeline
The [multiplex_primer_design](multiplex_primer_design.py) script provides a workflow for running the pipeline. For more details on arguments and defaults for each function, click on the function links to see documentation.

0. Set up a folder structure for storing inputs and outputs
1. [Batch Primer Design](docs/1_BatchPrimerDesign.md): Primers are designed for each template sequence using primer3, including predicting secondary structures (hairpins, homodimers, and heterodimers) within the primer pair.
2. [Filter Primers](docs/2_FilterPrimers.md): Primer pairs are filtered to avoid likely secondary structures based on Gibbs free energy (deltaG) and annealing temperatures. 
3. [Check Primer Specificity](docs/3_CheckPrimerSpecificity.md): Specificity of primers are checked against all templates. Any primer pairs that aren't specific to one locus are discarded to avoid off-target amplification.

   *Keeplist primers should be added to the previous step's output before proceeding to the next step. This is automatically handled in the multiplex_primer_design script.*

4. [Predict Cross-Primer Dimers](docs/4_DimerPrediction.md): Primer dimers are predicted using MFEprimer, with one MFEprimer output including all primer dimers expected to form and a second output including only primer dimers forming on the 3' end of primers.
5. [Tabulate Dimers](docs/5_TabulateDimers.md): Primer dimer reports are translated into tables counting pairwise primer pair interactions and total interactions per primer pair.

   *Simulated annealing parameter space may be explored before proceeding to optimization. See [plot_SA_parameters](docs/6A_ExploreOptimParameters.md) for details.*

6. [Optimize Multiplex Primer Set](docs/6_OptimizeMultiplexPrimerSet.md): A set of primers for "N" loci is selected that minimizes the number of secondary interactions (i.e., dimer load) between primer pairs. An initial primer set is selected using a pseudo-greedy algorithm where the primer pairs with the cumulative lowest dimer load (across all loci provided) are selected, then adaptive simulated annealing is used to explore the optimization space around this initial primer set by randomly swapping out primer pairs and keeping improvements while allowing for 'mistakes' that may improve the primer set in the long run, and finally the best primer set found during adaptive simulated annealing is entered into a simple iterative improvement algorithm where the worst loci are swapped for better alternatives. Because these algorithms have an element of randomness, try running this optimization step multiple times and choosing the best outcome.

   *The multiple_run_optimization script automates multiple runs of the optimization process. Code for running this can be found within the multiplex_primer_design script.*

# Notes on the optimization process
Optimization works best when there are >4 times as many candidates relative to the number of target loci. Problems with few candidates relative to targets (e.g., 50 inputs loci for a 40-locus target) are likely to underperform because relatively few combinations are possible. The number of candidate loci needed to build dimer-free primer sets increases substantially with the number of target loci because finding primers that do not interact becomes more difficult. For example, we have had success designing 50-plexes from 200 input sequences, but designing a 150-plex with minimal primer interactions required >2000 candidate sequences.

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
