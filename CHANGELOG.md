# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.1.2 - 2026-07-21
### Fixed
- Added cap to inner loop of GLS in the optimize_multiplex step. This is required to avoid getting stuck on one primer swap.

## 1.1.1 - 2026-07-21
### Fixed
- Added cap to sampling iterations for SetTemps in the optimize_multiplex step. This is required to avoid stalling in instances where few options are available or changes are almost always good.
- Typo in plotAmpliconTrees that caused function to crash

## [1.1.0] - 2026-07-16
### Added
- plotMismatches function in primerTree_specificity_checks.R
- multi-threading in runPrimerTree and plotAmpliconTrees
- multi-processing in multiplex-wormhole, which is carried into mw-primer-design, MFEprimer dimer, and mw-mult-optimization steps
### Fixed
- Fixed erroneous thermodynamic calculations in mw-specificity
- Fixed issues with KEEPLIST in mw-optimize-multiplex
- Improved functionality of inner loop in simple iterative improvement
- primerTree R script writes to an output file to save progress intermittently
- Improved compatibility with Windows
- Other minor bugs that occurred in niche cases
### Removed
- multiplex_primer_design.py script which was a deprecated duplicate of multiplexWormhole.py

## [1.0.5] - 2026-07-03
## Fixed
- Lots of minor bug fixes to improve documentation and stability across various use cases.

## [1.0.1] - 2026-06-25
- Initial public release.
