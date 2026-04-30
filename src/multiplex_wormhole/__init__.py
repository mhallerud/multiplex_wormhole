# LOAD ALL SUBMODULES
import sys, os, importlib, glob
sys.path.append(os.path.dirname(__file__))

import multiplexWormhole
import panel_assessment
import batch_primer3_design
import tabulate_dimers
import optimize_multiplex 
import multiple_run_optimization 
from helpers import setup_mfeprimer
from helpers import CSVtoFasta
from helpers import add_keeplist_to_fasta
plot_ASA_temps = importlib.import_module("plot_ASA_temps")



# SET NAMES FOR MAIN FUNCTIONS
assessPanel = panel_assessment.main
primer3BatchDesign = batch_primer3_design.main
tabulateDimers = tabulate_dimers.main
plotASAtemps = plot_ASA_temps.main
optimizeMultiplex = optimize_multiplex.main
multipleOptimizations = multiple_run_optimization.main
CSVtoFASTA = CSVtoFasta.main
addKeeplist2FASTA = add_keeplist_to_fasta.main
multiplexWormhole = multiplexWormhole.main


# SET DEPENDENCY PATH FOR MFEPRIMER
MFEprimer_PATH = glob.glob(os.path.dirname(__file__)+"/*mfeprimer*")[0]



__all__ = ["multiplexWormhole", "panel_assessment", "setup_mfeprimer", "batch_primer3_design", "tabulate_dimers", "optimize_multiplex", "multiple_run_optimization", "CSVtoFasta", "add_keeplist_to_fasta", "plot_ASA_temps"]
#, "assessPanel", "primer3BatchDesign", "tabulateDimers", "plotASAtemps", "optimizeMultiplex", "multipleOptimizations", "CSVtoFASTA", "addKeeplist2FASTA", "MFEprimer_PATH"]
