# multiplex_wormhole
Optimizing PCR primer design for multiplex amplicon sequencing by minimizing predicted dimer loads. The target audience is GT-seq panel development, however the process is transferable to any application where multiple distinct amplicons are being targeted for multiplex PCR.

multiplex_wormhole was created for optimizing large multiplex PCR assays for the purpose of SNP-based genotyping wildlife from noninvasive genetic samples. Defaults are tailored towards amplifying degraded DNA from low quality samples. 

**Important: This pipeline is still in development. Use at your own risk!**

# Dependencies
This pipeline relies on primer3 for primer design and MFEprimer for dimer calculations. To setup your files to run:
- Primer3 can be downloaded [here](https://github.com/primer3-org/primer3/releases). Keep track of where this downloads- the path to primer3_core will need to be added at line 55 in multiplex_primer_design.py.
- MFEprimer can be downloaded [here](https://www.mfeprimer.com/mfeprimer-3.1/#2-command-line-version). The path to MFEprimer-*-awd will need to be updated on line 54 of multiplex_primer_design.py. MFEprimer Version 3.2.7 on Darwin was used during development. 

The following (normally pre-installed) Python dependencies are required:
- pandas: install by running "pip install pandas" in the command line. Developed with pandas version 1.4.4
General Python modules required (these come pre-installed with many Python installations): os, sys, csv, random, math, signal, gc, itertools, shutil, mathplotlib

multiplex_wormhole was built and tested on MacOS with Python 3.9.13 in Spyder and bash 3.2.57(1)-release. 

# Running the pipeline
1. Create a CSV file with a row for each target and columns for locus ID, template sequence, and target position (following primer3 <start bp>,<length> format). For SNPs, the create_in_templates.R file takes a paired VCF and FASTA (FASTA loci match VCF CHROM field) and outputs the templates CSV. See example inputs in the [examples folder](https://github.com/mhallerud/multiplex_wormhole/examples).

2. Run multiplex_primer_design.py by changing the filepaths on lines 43, 54-60. Default parameters are provided throughout the script but may be changed as desired.

***This pipeline has very specific file requirements and naming conventions. For simplicity, it is recommended to run the pipeline from the beginning for each new run to avoid problems.***


# Steps in the pipeline
The [multiplex_primer_design](multiplex_primer_design.py) script provides a workflow for running the pipeline. For more details on arguments and defaults for each function, click on the function links to see documentation.

0. Set up a folder structure for storing inputs and outputs
1. [Batch Primer Design](docs/1_BatchPrimerDesign.md): Primers are designed for each template sequence using primer3, including predicting secondary structures (hairpins, homodimers, and heterodimers) within the primer pair.
2. [Filter Primers](docs/2_FilterPrimers.md): Primer pairs are filtered to avoid likely secondary structures based on Gibbs free energy (deltaG) and annealing temperatures. 
3. [Check Primer Specificity](docs/3_CheckPrimerSpecificity.md): Specificity of primers are checked against all templates. Any primer pairs that aren't specific to one locus are discarded to avoid off-target amplification.

   *Whitelist primers should be added to the previous step's output before proceeding to the next step. This is automatically handled in the multiplex_primer_design script.*

4. [Predict Cross-Primer Dimers](docs/4_PrimerPredictions.md): Primer dimers are predicted using MFEprimer, with one MFEprimer output including all primer dimers expected to form and a second output including only primer dimers forming on the 3' end of primers.
5. [Tabulate Dimers](docs/5_TabulateDimers.md): Primer dimer reports are translated into tables counting pairwise primer pair interactions and total interactions per primer pair.

   *Simulated annealing parameter space may be explored before proceeding to optimization. See [plot_SA_parameters](docs/6A_ExploreOptimParameters.md) for details.*

6. [Optimize Multiplex Primer Set](docs/6_OptimizeMultiplexPrimerSet.md): A set of primers for "N" loci is selected that minimizes the number of secondary interactions (i.e., dimer load) between primer pairs. An initial primer set is selected using a pseudo-greedy algorithm where the primer pairs with the cumulative lowest dimer load (across all loci provided) are selected, then adaptive simulated annealing is used to explore the optimization space around this initial primer set by randomly swapping out primer pairs and keeping improvements while allowing for 'mistakes' that may improve the primer set in the long run, and finally the best primer set found during adaptive simulated annealing is entered into a simple iterative improvement algorithm where the worst loci are swapped for better alternatives. Because these algorithms have an element of randomness, try running this optimization step multiple times and choosing the best outcome.

   *The multiple_run_optimization script automates multiple runs of the optimization process. Code for running this can be found within the multiplex_primer_design script.*

# Notes on the optimization process
Optimization functions best when there are 3-5 times as many potential templates relative to the number of desired loci. The process will likely perform poorly on problems with few targets relative to the number of desired loci (e.g., 50 input loci with a desired SNP panel of 40 loci). The number of potential templates needed to build dimer-free sets increases geometrically with the number of loci desired in the final panel. We have had success with designing 50-plexes from 200 input templates and 150-plexes from 2000 input templates.

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
