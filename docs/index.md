---
layout: default
nav_order: 2
has_children: true
---
# multiplex wormhole
*In silico* multiplex PCR primer design for noninvasive wildlife genetics.

![Multiplex wormhole logo with DNA entering black hole](assets/images/logo.png)

# Contents
1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Input File Format](#input-file-format)
4. [Multiplex Wormhole Functions](#multiplex-wormhole-functions)
5. [Recommended Workflows](#recommended-workflows)
6. [Handling Outputs](#handling-outputs)
7. [Optimization Details](#optimization-details)
8. [Contact](#problems-questions)
9. [Success Stories](#success-stories)
10. [Citations](#citations)


**multiplex wormhole is largely an automation of pre-existing protocols for noninvasive SNP panel development (Eriksson et al. 2020). All default parameters for primer design, PCR parameters, and dimer prediction are based on this protocol and assume that the two-step PCR with Illumina Nextera adapters is being used. We recommend that the indexing step use Illumina unique dual indexes to minimize error rates due to index hopping, which is likely to be especially problematic with low-quality samples.**


## Installation
1.  **Set up a Python virtual environment**
If you already use Python, we recommend setting up a virtual environment to avoid affecting your pre-existing scripts. [Anaconda](https://www.anaconda.com/products/navigator) is the recommended virtual env manager for Python. We recommend the Spyder Python interpreter for running code line-by-line, other options are PyCharm and VScode. 

To set up a unix virtual environment with conda or pixi, run `pixi init` or `conda init` within the folder from which you plan to run multiplex_wormhole. Then, to enter your python environment, run `pixi shell` or `conda activate`. 

2. **Install dependencies**
multiplex_wormhole depends on the Python libraries primer3-py for primer design, pandas and numpy for data frame handling, and matplotlib for plotting. The required versions can be installed on the command line/terminal using:
```
pip install pandas==1.4.4
pip install numpy==1.24.4
pip install matplotlib==3.5.2
pip install primer3-py==2.0.0
```

3. **Download multiplex wormhole**
multiplex wormhole scripts can be downloaded from Github. On the command line, run:
```
git clone https://github.com/mhallerud/multiplex_wormhole
```

4. **Configure MFEprimer dependency**
MFEprimer (Wang et al. 2019) is used for dimer prediction in multiplex wormhole. To install and configure the MFEprimer binaries, run from the command line:
```
python3 multiplex_wormhole/src/scripts/setup_mfeprimer.py
```
or, from within a Python interpreter:
```
# load multiplex wormhole functions
import sys
sys.path.append("/multiplex_wormhole/src")
from setup_mfeprimer import main as setup_mfeprimer
setup_mfeprimer()
```

If this fails, MFEprimer3.2.7 can be downloaded manually [here](https://github.com/quwubin/MFEprimer-3.0/releases), then unzip, ensure permissions allow execution (chmod +x), and save to your ~/multiplex_wormhole/src directory. 

Now you are ready to run multiplex wormhole!

## Quick Start
Command line usage:
```
# panel design
# -d: deltaG optimization (details in Optimization section below
# -v: verbose
python3.9 multiplex_wormhole.py -t TEMPLATES -n NLOCI -o OUTDIR [-p PREFIX] [-k KEEPLIST] [-r RUNS] [-i ITER] [-s SIMPLE] [-d] [-v]

# panel assessment
panel_assessment.py [PRIMERFASTA | PRIMERCSV]
```

Python usage:
```
# setup
import sys
sys.path.append("YOURPATHTO/multiplex_wormhole")

# panel design
from multiplex_wormhole import main as multiplex_wormhole
multiplex_wormhole(TEMPLATES,
                   N_LOCI,
                   OUTDIR,
                   KEEPLIST,
                   N_RUNS,
                   ITERATIONS,
                   SIMPLE,
                   deltaG,
                   VERBOSE)

# panel assessment
from panel_assessment import main as assess_panel
assess_panel(PRIMERFASTA)
# or
assess_panel(PRIMERCSV)
```

**Using the multiplex_wormhole function will use defaults for all other parameters. Check out the workflow below and specific functions for higher flexibility in settings.**


## Input file format
### Templates file
The primary input to multiplex wormhole is a CSV format containing information on template DNA sequences and target regions. multiplex wormhole was developed for SNP-based genotyping with data originating from new or previously published genomic sequencing data or pre-existing SNP panels, but the method can be easily adjusted to accommodate multiplex PCR design for other targeted sequencing applications such as microsatellites, methylation regions, or functional regions. 

Each DNA sequence in the FASTA file will be treated as a unique target for PCR amplification (i.e., sequences should be non-overlapping and unique). The file should have 3 columns named:
* SEQUENCE_ID : template names, without punctuation. All punctuation will be automatically removed. Names should also be unique. 
* SEQUENCE_TEMPLATE : DNA template sequence in the 5'-->3' direction.
* SEQUENCE_TARGET : identifies the base pairs targeted for PCR ampflication following primer3 format: <startBP,length>. For example, a SNP at the 100th base pair in the sequence would be denoted as 100,1 in this field. If there are 2 SNPs, for example at the 50th and 90th base pairs in the sequence, 50,40 would be the target. See the [example input CSV](https://github.com/mhallerud/multiplex_wormhole/blob/main/examples/Input_Templates.csv). 


| SEQUENCE_ID   | SEQUENCE_TEMPLATE    | SEQUENCE_TARGET    |
| ------------- | -------------------- | ------------------ |
| CLocus_704    | TCAGAGAC...          | 53,1               |
| ...           | ...                  | ...                |


The [create_in_templates](https://github.com/mhallerud/multiplex_wormhole/blob/main/src/scripts/create_in_templates.R) R script can be used to create this file using VCF and FASTA as inputs. The R script handles two input types:

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

### Keeplist file (Optional)
Users can also add a "keeplist" FASTA file of primers which must be included in the final multiplex, for example primers from a pre-existing panel or representing functionally important regions (e.g., sexing loci). This option can be invoked for augmenting existing panels. Importantly, primer design settings in the batch primer design step (details below) should be adjusted to match the PCR settings used to design the keeplist primers. 

### A Note on Candidate Loci
Multiplex wormhole treats all candidates equally, it is therefore the responsibility of the user to select informative candidate DNA sequences based on their objectives. 

For example, I use the following steps to prepare input data for applications focused on individual identification:
1. *SNP discovery*: Identifying SNPs from reduced representation sequencing (e.g., RADseq) or whole genome sequencing (WGS) projects, if available. Alternatively, SNP data can be pulled from pre-existing panels, including Fluidigm panels.
2. *Initial SNP filtering*: If using genomic data, follow standard genomics guidelines to remove erroneous SNPs (e.g., removing SNPs with low quality scores, low depths or extremely high depths, and high missingness; O'Leary et al. 2018). 
3. *Candidate SNP filtering*: SNPs with high information content based on the objective(s) of the panel are selected. For example, SNPs with high minor allele frequencies are selected for individual identification, while SNPs with high Fst or delta values are selected for population assignment. 
4. *Extract flanking regions around SNPs*: Using --fasta-loci for de novo RADseq data and bedtools getfasta for reference-aligned RADseq or WGS data. For reference-aligned data, masking repetitive regions (e.g., using RepeatMasker) is recommended prior to this step.
5. *Remove candidates in known repetitive regions*: Check alignments against known repetitive elements using (CENSOR](https://www.girinst.org/censor/index.php).


## Multiplex Wormhole Functions
For maximum flexibility, the full multiplex wormhole workflow is available at: [multiplex_primer_design.py](https://github.com/mhallerud/multiplex_wormhole/blob/main/src/multiplex_primer_design.py). Click on the functions below for detailed information on inputs, outputs, and settings (including defaults). 

1. [Batch Primer Design](1_BatchPrimerDesign.md) with `primer3_batch_design.py`
   
2. [Dimer Prediction](2_DimerPrediction.md) with `MFEprimer dimer`
  
3. [Tabulate Dimers](3_TabulateDimers.md) with `tabulate_dimers.py`

   (Optional): [Explore Optimization Parameters](4A_ExploreOptimParameters.md) with `plot_ASA_temps.py`
  
4. [Optimize Multiplex Primers](4_OptimizeMultiplexPrimerSet.md) with `optimize_multiplex.py`

5. [Multiple Run Optimization](5_MultipleRunOptimization.md) with `multiple_run_optimization.py`

6. [Panel Assessment](6_AssessPanel.md) with `panel_assessment.py`

Helper Functions:
* [Convert CSV to FASTA](CSVtoFASTA.md) with `CSVtoFASTA.py`
* [Add Keeplist to FASTA](AddKeeplistToFASTA.md) with `add_keeplist_to_fasta.py`

![diagram showing mw workflow](assets/images/diagram.png)
Diagram for multiplex wormhole workflow, with black arrows showing the panel design workflow and blue arrows showing the panel assessment workflow.


## Recommended Workflows
### Designing a Novel Panel
This workflow will optimize a novel panel from data including the full template DNA sequences (i.e., target + flanking regions).

TEMPLATES --> primer3_batch_design.py --> MFEprimer dimer --> tabulate_dimers.py --> multiple_run_optimization.py

### Augmenting an Existing Panel
This workflow will optimally expand an existing multiplex primer set (KEEPLIST) using new DNA template data (TEMPLATES).

TEMPLATES + KEEPLIST --> primer3_batch_design.py --> MFEprimer dimer --> tabulate_dimers.py --> multiple_run_optimization.py

### Assessing an Existing Panel
This workflow assesses dimer load of an existing multiplex primer set (PRMERS).

PRIMERS --> panel_assessment.py

### Re-designing an Existing Panel
An existing panel (or combined panels) can be redesigned, but the approach will depend on whether only primer sequences available or the full templates are available. If template data are available, ensure that there is one template per target and use the "novel panel" or "augmented panel" (if unequal data availability) workflows. If only primer sequences are available, the only approach for optimization is to reduce panel size. To do so, use the following workflow with the original PRIMERS as the input fasta at the optimization step and set N_LOCI smaller than the current panel:

PRIMERS --> MFEprimer dimer --> tabulate_dimers.py --> multiple_run_optimization.py


## Handling Outputs
Multiplex wormhole is an ***in silico*** design tool, and although it facilitates the primer design process, additional checks and lab optimization remain critical to panel success. Importantly, multiplex wormhole does not do specificity checking within the design process. The following steps should be taken prior to ordering primers.

### Specificity Checks
1. *Specificity check against the target species*: [PRIMER-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) can be used to check each set of primers against the genome of the target species (or a closely related species) to ensure that targets only amplify one target in the genome. This is particularly important in noninvasive samples because sequence variation due to paralogs (i.e., regions present multiple times across the genome) is not differentiable from sequence dropout, which is the most common form of genotyping error. Paralogs in highly variable regions (e.g., mRNA transcript variants) can be particularly problematic because dropout may be linked to variation in the primer binding region which causes biased amplification failure. 
2. *Specificity check against diet species / potential field contaminators*: Ensuring species specificity is critical to accuracy if the panel is to be applied to noninvasive samples. Specificity checking uses the same process as described above, but against genomes and GenBank for species present in the diet. 
3. *Secondary dimer check*: Dimer prediction is an imperfect science, and primers should be checked for particularly bad dimers e.g.  by rerunning with MFEprimer v3 and/or other software such as [MFEprimer v4](https://m4.igenetech.com/dimer), [primer-dimer.com](www.primer-dimer.com), and/or [ThermoFisher's Multiple Primer Analyzer](https://www.thermofisher.com/us/en/home/brands/thermo-scientific/molecular-biology/molecular-biology-learning-center/molecular-biology-resource-library/thermo-scientific-web-tools/multiple-primer-analyzer.html).

### Lab Testing
We recommend the protocol in Eriksson et al. (2020) for lab testing:
1. Initial lab test in equimolar multiplex (0.2 uM per primer) on a range of sample qualities. Primer pairs that fail to amplify at this step and/or primers forming strong secondary structures that consume sequencing depth should be removed. Primer concentrations can also be adjusted to improve evenness of sequencing depth across targets.
2. The second round of optimization should test specificity and genotyping error rates. These can be assessed using matched samples from the same individual (e.g., scat, hair, and tissue from the same animal) and testing on off-target species (e.g., diet items and closely related co-occurring species). Note that technical replicates are better than nothing, but are insufficient for assessing specificity and error issues specific to low-quality samples.


## Optimization Details
The multiplex wormhole `optimize_multiplex.py` function uses a combination of simple iterative improvement and simulated annealing algorithms to minimize dimer load in the multiplex primer set, where dimer load can be measured by a) tallying pairwise dimers or b) maximizing mean Gibbs free energy (deltaG; a meeasure of dimer strength) of pairwise dimers. The optimization process heavily relies on having abundant candidates relative to the number of target loci, and more templates will be needed to design larger panels while keeping dimer load low. See [Optimization Algorithm](https://mhallerud.github.io/multiplexwormhole/optimization-process) for details. 


## Problems? Ideas?
Problems, bugs, or thoughts for enhancement can be reported on the [GitHub Issues page](https://github.com/mhallerud/multiplex_wormhole/issues). You can also contact Maggie Hallerud (hallerum@oregonstate.edu).


## Success Stories
Multiplex wormhole has been used to develop noninvasive genotyping panels for:

![Humboldt marten](assets/images/marten.png)
![Gray wolf](assets/images/graywolf.png)
![Pacific fisher](assets/images/fisher.JPG)
![Sea_otter](assets/images/seaotter.png)
![Cougar](assets/images/cougar.png)
![Bobcat](assets/images/bobcat.png)
![American_black_bear](assets/images/blackbear.png)

Contact us if you want to add your species to the list!


## Citations
Eriksson, CE, Ruprecht J, Levi T. 2020. More affordable and effective noninvasive SNP genotyping using high-throughput amplicon sequencing. Molecular Ecology Resources 20(4): [doi: 10.1111/1755-0998.13208](https://doi.org/10.1111/1755-0998.13208).

O'Leary et al. 

Untergasser, A, Cutcutache, I, Koressaar, T, Ye, J, Faircloth, BC, Remm, M, Rozen SG. 2012. Primer3--new capabilities and interfaces. Nucleic Acids Research: e115. [doi: 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596). 

Wang et al. 2019
