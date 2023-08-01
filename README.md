# multiplex_wormhole
Optimizing PCR primer design for multiplex amplicon sequencing.

# setup
This pipeline relies on primer3 for primer design. Primer3 can be downloaded (here)[https://github.com/primer3-org/primer3/releases]. The path to primer3_core on line will need to be updated on line 18 of scripts/primer3.sh will need to be updated to reflect your local version.

# Steps to run the pipeline
1. Create a CSV file with a row for each target and columns for locus ID, template sequence, and target position (following primer3 <start bp>,<length> format). For SNPs, the create_in_templates.R file takes the input VCF and matching FASTA and outputs the templates CSV.

2. Run primer_design_STEP1.sh to design and filter primers for your targets:
   ./primer_design_STEP1.sh <TEMPLATE_CSV> <OUTPUT_DIRECTORY> <RESULTS_PREFIX>

   To specify input settings: The primer3 settings used for primer design are found in primer3_settings/primer3_Base_NoSecondaryFilters.txt and primer3_settings/primer3_Broad_NoSecondaryFilters.txt. These settings can be edited directly, or you can create your own Broad and Base setting files and change the filepaths in the primer3_batch_design.sh script.

   To specify filtering criteria: Line 30 in the primer_design_STEP1.sh file can be edited to specify different filtering parameters. The format is:
   ./scripts/filter_primers_Tm_dG.sh <max_Tm> <min_hairpin_dG> <min_ends_dG> <min_self_dG> <OUTPUT_DIRECTORY> <RESULTS_PREFIX> > <LOGFILE>
   
3. Copy the *_specificityCheck_passed.csv file into (primer-dimer.com)[https://primer-dimer.com], select 'Multiplex Analysis' and 'Dimer Structure Report', and click Submit. Depending on how many loci you provided, this step may take awhile (~20 minutes for primers for 50 loci, ~xx hours for 1200 loci).

4. Run translate_primerSuite_report.R on the PrimerDimerReport.txt file downloaded from primer-dimer.com to convert this to a CSV. This can either be run directly in R or via the command line following the syntax:
   Rscript --vanilla scripts/translate_primerSuite_report.R <PATH_TO_PrimerDimerReport.txt> <DELTA_G_THRESHOLD>
   
5. Run optimize_primers.py in Python: Edit user-defined parameters at the top of the file, then you can run the whole script.


# Process
1. Primers are designed for each template sequence using primer3, including calculating secondary structures using Illumina i5 and i7 adapters as overhangs.
  Defaults:
   - Annealing Temp: 52 C
   - Amplicon size: 70-120 bp
   - Primer size: 18-26 bp (optimal: 20)
   - Primer Tm: 57-63 Celsius (optimal: 60)
   - Primer GC content: 30-70% (optimal: 50)
   - GC Clamp: 1
   - Max End GC: 4
   - Max Poly X: 4
  If primers can't be found for the above settings, constraints are relaxed:
    - Primer GC content: 20-80%
    - GC Clamp: 0
    - Max End GC: 5
    - Max Poly X: 5
2. Primer pairs are filtered to avoid secondary structures. By default, primer pairs are discarded if the following criteria are met:
    - secondary structure Tm > 45
        AND
    - delta G < threshold
        - hairpins: -3 kcal/mol
        - homodimers or heterodimers at primer ends: -5 kcal/mol
        - homodimers or heterodimers not at primer ends: -10 kcal/mol
3. Specificity of primers are checked against all templates. Any primer pairs that aren't specific to one locus are discarded.
4. A primer dimer table is calculated which includes the number of primer dimers estimated for all combinations of primer pairs. Currently this is being done using PrimerSuite's Primer Dimer function, but in the future the goal will be to add downstream capability for ThermoFisher's Multiplex Primer Design and SADDLE DimerScores.
  Primer dimers are only considered if their delta G is below a given threshold (default: -6 kcal/mol).
  In the future, functionality may be added to filter primer dimers based on Tm and location (end vs. middle).
5. A set of primers for "N" loci is selected that minimizes the number of negative interactions between primer pairs. This algorithm was largely inspired by SADDLE and Marxan Conservation Planning software. Specifically, the following steps are followed:
     - The best primer pairs (i.e., the primer pairs with the lowest dimer load) are identified for each locus.
           NOTE: Dimer load can be calculated either as the total number of primer dimers that loci contributes, or as the number of binary primer pair interactions (1 = primer pair has at least one primer dimer, 0 = primer pair has no interactions) the loci contributes. To optimize based on total primer dimers, use *_PrimerPair_interactions_sum.csv and *_PrimerPair_interactions_wide.csv (output from translate_primerSuite_report.R) as the DIMER_SUMS and DIMER_TABLE inputs. To optimize based on binary primer pair interactions, use *_PrimerPair_interactions_sum_binary.csv and *_PrimerPair_interactions_wide_binary.csv instead.
     - An initial set of primer pairs is selected by choosing the best (i.e., lowest dimer load) "N" primer pairs from the locus-specific bests.
     - The algorithm then iteratively improves this initial set of primer pairs by switching the worst (i.e., highest dimer load) primer pair with a random alternative (either for the same locus, or from a locus that is not currently included in the set).
         - If the change resulted in an improvement (i.e., fewer TOTAL dimer load for the set), then this new set is used as the starting point for the next iteration.
         - If the change did not improve the set, then other alternatives for this worst pair are tested until a replacement that improves the set is found. If no alternatives improve the set, then this primer pair is blacklisted (i.e., kept in the set but a message will be returned suggesting to possibly remove this pair as an option) and the algorithm moves on to the next worst primer pair.
      - The algorithm stops when one of the following conditions are met:
         - The total primer dimer load is 0, and therefore there is no longer a need to optimize.
         - No further improvements can be made to the set (i.e., all remaining primer pairs have been blacklisted- no alternatives to these will result in improvements).
         - A user-defined number of iterations is reached.

# Notes on the Optimization Process
Due to the random nature of the optimization algorithm, the returned primer set will only be an estimate of the ideal solution and will vary with each run. I recommend running the optimization 10-20 times and keepign the best set of these runs, unless a set with 0 dimers is found or the problem is simple (few desired # loci vs. many options).

This pipeline was created for designing multiplex PCR primers for SNP genotyping from reduced-representation sequencing data. The pipeline functions best when there are many potential targets relative to the number of desired loci (i.e., there are many alternatives that can be tested). The process will likely perform poorly on problems with few targets relative to the number of desired loci (e.g., 50 input loci with a desired SNP panel of 40 loci).
