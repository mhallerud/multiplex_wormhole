from scripts.primer3_batch_design import main as primer3BatchDesign
from scripts.tabulate_MFEprimer_dimers import main as tabulateDimers
from scripts.optimize_primers import main as optimizeMultiplex
from scripts.multiple_run_optimization import main as multipleOptimizations
from scripts.CSVtoFasta import main as CSVtoFASTA

import importlib
plotASAtemps = importlib.import_module("plot_ASA_temps")
