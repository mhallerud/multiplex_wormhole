import sys, os
sys.path.append(os.path.dirname(__file__))

from multiplex_wormhole import main as multiplex_wormhole
from panel_assessment import main as assessPanel

from scripts.setup_mfeprimer import main as setup_mfeprimer
from scripts.primer3_batch_design import main as primer3BatchDesign
from scripts.tabulate_MFEprimer_dimers import main as tabulateDimers
from scripts.optimize_primers import main as optimizeMultiplex
from scripts.multiple_run_optimization import main as multipleOptimizations
from scripts.CSVtoFasta import main as CSVtoFASTA

import importlib
plotASAtemps = importlib.import_module("plot_ASA_temps")

import glob
MFEprimer_PATH = glob.glob(os.path.dirname(__file__)+"/*mfeprimer*")[0]
