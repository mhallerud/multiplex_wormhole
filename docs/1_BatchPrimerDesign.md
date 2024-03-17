# Batch Primer Design with `primer3_batch_design`

## Purpose
Primers are designed for each template sequence using primer3, including calculating secondary structures using Illumina i5 and i7 adapters as overhangs.

## Defaults:
Narrow primer3 settings (see --- file):
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

If primers can't be found for the above settings, constraints are relaxed (see --- file):
    - Primer GC content: 20-80%
    - GC Clamp: 0
    - Max End GC: 5
    - Max Poly X: 5

## Usage

