# multiplex_wormhole
*In silico* optimization for multiplex PCR assays by minimizing predicted primer dimer loads. The target audience is GT-seq panel developers, however the process is transferable to any application where many amplicons are targeted in multiplex PCR.

multiplex_wormhole was created for the purpose of optimizing large multiplex PCR panels that are applied to noninvasive wildlife genetic samples for the purposes of SNP-based individual identification and relatedness analysis.  Default settings are therefore tailored towards amplifying degraded DNA from low quality samples. 

**Important: This pipeline is still in development. Use at your own risk!**
**Please report any problems or potential enhancements in the GitHub Issues page.**


# Dependencies
- Primer3 is used for primer design and can be downloaded [here](https://github.com/primer3-org/primer3/releases). Keep track of where this downloads- the path to primer3_core will need to be updated at line 55 in multiplex_primer_design.py and line 89 in multiplex_wormhole.py.
- MFEprimer is used for dimer calculations and can be downloaded [here](https://www.mfeprimer.com/mfeprimer-3.1/#2-command-line-version). The path to MFEprimer-*-awd will need to be updated on line 54 of multiplex_primer_design.py and line 88 in multiplex_wormhole.py. MFEprimer Version 3.2.7 on Darwin was used during development. 

The following (normally pre-installed) Python dependencies are required:
- pandas: install by running "pip install pandas" in the command line. Developed with pandas version 1.4.4
General Python modules required (these come pre-installed with most Python installations): os, sys, csv, random, math, signal, gc, itertools, shutil, mathplotlib, datetime

multiplex_wormhole was built and tested on MacOS with Python 3.9.13 in the Spyder IDE.


# Inputs and Outputs
## Inputs
The only required input to multiplex wormhole is a table with DNA sequences including target and flanking regions (see step 1 under quick-start guide below for details). Optionally, users can also add a "keeplist" FASTA file of primers which must be included in the final multiplex, for example primers from a pre-existing panel that is being added onto. *Note that multiplex wormhole designs primers based on standardized primer3 PCR conditions including a melting temperature between 58-62 degrees. If keeplist primers were designed using different PCR conditions and primer design settings, these should be adjusted in the files found in the primer3_settings directory.*

I use the following steps to prepare input data:
1. *SNP discovery* using data from reduced representation sequencing (e.g., RADseq) or whole genome sequencing (WGS) projects
2. *Initial SNP filtering* following standard genomics guidelines to remove erroneous SNPs (e.g., removing SNPs with low quality scores, low depths or extremely high depths, and high missingness; O'Leary et al. 2018)
3. *Target SNP filtering* where SNPs with high information content based on the objective(s) of the panel are selected. For example, SNPs with high minor allele frequencies are selected for individual identification, while SNPs with high Fst or delta values are selected for population assignment. Multiplex wormhole treats all candidates equally, it is therefore the responsibility of the user to select informative candidate DNA sequences based on their objectives. 
4. *Extract flanking regions around SNPs* using --fasta-loci for RADseq data and bedtools getfasta for WGS data
5. *Remove candidates in known repetitive regions* by checking alignments to known repetitive elements using (CENSOR](https://www.girinst.org/censor/index.php)

## Outputs
Multiplex wormhole is an **in silico** design tool and additional testing of primers is needed. For example, here is my process for further checking and optimizing multiplex wormhole output:
1. *Specificity check against prey species* using [PRIMER-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)- this is especially important when genotyping scat samples
2. *Initial lab testing in equimolar multiplex* to assess amplification success in target samples, test specificity against prey and closely related co-occurring species, and estimate genotyping error rates based on PCR replicates and/or comparison of high-quality and low-quality samples from the same known individuals

# Quick-Start Guide
## 1. Set up input data
Create a CSV file with a row for each candidate and columns named SEQUENCE_ID containing sequence names, SEQUENCE_TEMPLATE containing the template DNA sequence in 5'-->3' direction, and SEQUENCE_TARGET containing the base pairs targeted for PCR amplification following primer3 format: *startBP*,*length*. See example inputs in the [examples folder](https://github.com/mhallerud/multiplex_wormhole/examples). 

Each DNA sequence in the FASTA file will be treated as a unique target for PCR amplification (i.e., sequences should each include one "target" region). 

The [create_in_templates](scripts/create_in_templates.R) R script can be used to create this file using a VCF and FASTA as inputs. The R script handles two input types:

* *de novo*: Matches formatting output of *de novo* Stacks pipeline with FASTA created by `populations --fasta-loci`. Assumes the VCF 'CHROM' field matches the FASTA sequence headers with a set prefix such as "CLocus_".

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

## 2. Run the multiplex_wormhole pipeline
### Running the full pipeline
Edit filepaths to primer3 and MFEprimer dependencies on lines 88-89 of multiplex_wormhole.py (this only needs to be done once), then the script can be run in python:
```
# add multiplex wormhole folder to your PATH file
import sys
sys.path.append("~/multiplex_wormhole")
# load the module
from multiplex_wormhole import main as multiplex_wormhole
# run multiplex wormhole with defaults
multiplex_wormhole(TEMPLATES, N_LOCI, OUTDIR, PREFIX=None, KEEPLIST_FA=None, N_RUNS=10, ITERATIONS=5000, SIMPLE=2000)
```
or from the linux command line:
```
# run multiplex wormhole with all options
python3 multiplex_wormhole.py TEMPLATES N_LOCI OUTDIR PREFIX KEEPLIST_FA N_RUNS ITERATIONS SIMPLE
# run multiplex wormhole with defaults
python3 multiplex_wormhole.py TEMPLATES N_LOCI OUTDIR 'None' 'None' 10 5000 2000
```

### Running piecewise steps of the pipeline
The multiplex_primer_design.py script is set up for running multiplex_wormhole one step at a time within a Python environment. Update the paths to your multiplex_wormhole scripts on line 44, dependencies on lines 54-55, then inputs can be specified on lines 58-62. Default parameters are provided for all steps but can be adjusted within the script. 
### Parameters
* **TEMPLATES**: Input CSV containing candidate DNA template IDs, sequences, and targets.
* **N_LOCI**: Number of target sequences in final multiplex primer set.
* **OUTDIR**: Directory where outputs and intermediates will be stored.
* **PREFIX**: Prefix for output optimization files.
* **KEEPLIST_FA**: FASTA containing primers that are required to be included in the final multiplex. These should be in the 5'-->3' direction and should include adapters, if primers are ordered with adapter sequences. Primer names must be in the format `LOCUS_ID`.`Pair`.`FWD/REV` to avoid designing multiple primer pairs for the same LOCUS_ID.
* **N_RUNS**: Number of times to run the full optimization process.
* **ITERATIONS**: Number of iterations to run adaptive simulated annealing algorithm (see [Optimize Multiplex Primer Set](docs/6_OptimizeMultiplexPrimerSet.md) for details).
* **SIMPLE**: Number of iterations to run simple iterative improvement algorithm (see [Optimize Multiplex Primer Set](docs/6_OptimizeMultiplexPrimerSet.md) for details).
Note that optimization will stop before the full number of iterations is run if a primer set with 0 dimers is found.
Multiplex_wormhole uses many additional parameters including filtering thresholds for including primers or dimers, simulated annealing algorithm parameters, and output names. Defaults for these are set in the multiplex_wormhole.py script but can be adjusted as desired, or you can run the pipeline one step at a time to further explore these parameters.


# Workflows for Various Applications
A) Designing a new panel --> Start at Step 1

B) Test dimer load of a pre-existing panel --> Steps 4-5

C) Combining pre-existing panels or primer sets --> Start at Step 3

D) Expand a current panel with newly designed primers --> Start at Step 1

# Steps in the Pipeline
See [multiplex_primer_design](multiplex_primer_design.py) to run the full pipeline.

0. Set up a folder structure for storing inputs and outputs
1. [Batch Primer Design](docs/1_BatchPrimerDesign.md)
Primers are designed for each template sequence using primer3, including predicting secondary structures (hairpins, homodimers, and heterodimers) within the primer pair.
2. [Filter Primers](docs/2_FilterPrimers.md) 
Primer pairs are filtered to avoid likely secondary structures based on Gibbs free energy (deltaG) and annealing temperatures. 
3. [Check Primer Specificity](docs/3_CheckPrimerSpecificity.md)
Specificity of primers are checked against all templates. Any primer pairs that aren't specific to one locus are discarded to avoid off-target amplification.

   *Keeplist primers should be added to the previous step's output before proceeding to the next step. This is automatically handled in the multiplex_primer_design script.*

4. [Predict Cross-Primer Dimers](docs/4_DimerPrediction.md): Primer dimers are predicted using MFEprimer, with one MFEprimer output including all primer dimers expected to form and a second output including only primer dimers forming on the 3' end of primers.
5. [Tabulate Dimers](docs/5_TabulateDimers.md): Primer dimer reports are translated into tables counting pairwise primer pair interactions and total interactions per primer pair.

   *Simulated annealing parameter space may be explored before proceeding to optimization. See [plot_SA_parameters](docs/6A_ExploreOptimParameters.md) for details.*

6. [Optimize Multiplex Primer Set](docs/6_OptimizeMultiplexPrimerSet.md): A set of primers for "N" loci is selected that minimizes the number of secondary interactions (i.e., dimer load) between primer pairs. An initial primer set is selected using a pseudo-greedy algorithm where the primer pairs with the cumulative lowest dimer load (across all loci provided) are selected, then adaptive simulated annealing is used to explore the optimization space around this initial primer set by randomly swapping out primer pairs and keeping improvements while allowing for 'mistakes' that may improve the primer set in the long run, and finally the best primer set found during adaptive simulated annealing is entered into a simple iterative improvement algorithm where the worst loci are swapped for better alternatives. Because these algorithms have an element of randomness, try running this optimization step multiple times and choosing the best outcome.

   *The multiple_run_optimization script automates multiple runs of the optimization process. The **multiple_run_optimization** script automates this process and is included in the multiplex_primer_design script.*

# A Note on the Optimization Process
Optimization works best when there are >4 times as many candidates relative to the number of target loci. Problems with few candidates relative to targets (e.g., 50 inputs loci for a 40-locus target) are likely to underperform because relatively few combinations are possible. The number of candidate loci needed to build dimer-free primer sets increases substantially with the number of target loci because finding primers that do not interact becomes more difficult. For example, we have had success designing 50-plexes from 200 input sequences, but designing a 150-plex with minimal primer interactions required >2000 candidate sequences.

# Alternative dimer calculation tools
The pipeline is built to use MFEprimer dimer to calculate dimer formation, however the optimization process will accept any input tables as long as the 2 input tables specify 1) pairwise dimer loads between primer pairs and 2) the total dimer load per primer pair, with primer pair IDs matching between the input templates and both tables. See example inputs in the [examples folder](https://github.com/mhallerud/multiplex_wormhole/examples).

One alternative to MFEprimer dimer is Primer Suite's Primer Dimer tool, though we have found dimers to be overestimated using this tool. To use this tool to calculate dimers:
1. Copy the *_specificityCheck_passed.csv file into [primer-dimer.com](https://primer-dimer.com), select 'Multiplex Analysis' and 'Dimer Structure Report', and click Submit. Depending on how many loci you provided, this step may take awhile (~20 minutes for primers for 50 loci, multiple hours for 1200 loci).

2. Run scripts/translate_primerSuite_report.R on the PrimerDimerReport.txt file downloaded from primer-dimer.com to convert this to a CSV. The delta G threshold specified will filter out any primer dimers with delta G above this value.

Another alternative is [ThermoFisher's Multiplex Primer Design](https://www.thermofisher.com/us/en/home/brands/thermo-scientific/molecular-biology/molecular-biology-learning-center/molecular-biology-resource-library/thermo-scientific-web-tools/multiple-primer-analyzer.html), however multiplex_wormhole currently has no support for translating these outputs into tables.


## Citations
Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M and Rozen SG. 2012. Primer3--new capabilities and interfaces. *Nucleic Acids Research* 40(15):e115. [http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3424584/](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3424584/)

Wang K, Li H, Xu Y, Shao Q, Yi J, Wang R, Cai W, Hang X, Zhang C, Cai H, and Qu W. 2019. MFEprimer-3.0: quality control for PCR primers.
*Nucleic Acids Research* 47(W1): W610–W613. [https://doi.org/10.1093/nar/gkz351](https://doi.org/10.1093/nar/gkz351)

Xie NG, Wang MX, Song P, Mao S, Wang Y, Yang Y, Luo J, Ren S, and Zhang DY. 2022. Designing highly multiplex PCR primer sets with Simulated Annealing Design using Dimer Likelihood Estimation (SADDLE). *Nature Communications* 13(1): 1881. doi: [10.1038/s41467-022-29500-4](https://doi.org/10.1038/s41467-022-29500-4). 
