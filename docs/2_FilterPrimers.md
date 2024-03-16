# Filtering Primer Pairs <filter_primers>

## Purpose
Primer pairs are filtered to avoid within-pair secondary structures. 

## Usage

## Defaults 
Thresholds for discarding primer pairs are set directly in the function. By default in <multiplex_primer_design>, primer pairs are discarded if the following criteria are met:
- secondary structure Tm > 45 C
  AND
- delta G < threshold where the threshold varies depending on the structure type:
  - hairpins: -2 kcal/mol
  - homodimers or heterodimers at primer ends: -5 kcal/mol
  - homodimers or heterodimers not at primer ends: -10 kcal/mol
