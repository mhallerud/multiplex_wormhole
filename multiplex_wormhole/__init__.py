# LOAD ALL SUBMODULES
import sys, os, importlib, glob, types
sys.path.append(os.path.dirname(__file__))

from ./multiplex_wormhole import main as multiplex_wormhole
from . import panel_assessment
from . import setup_mfeprimer
from . import batch_primer3_design
from . import tabulate_dimers
from . import optimize_multiplex 
from . import multiple_run_optimization 
from . CSVtoFasta
from . import add_keeplist_to_fasta
plot_ASA_temps = importlib.import_module("plot_ASA_temps")



# SET NAMES FOR MAIN FUNCTIONS
assessPanel = panel_assessment.main
primer3BatchDesign = batch_primer3_design.main
tabulateDimers = tabulate_dimers.main
plotASAtemps = plot_ASA_temps.main
optimizeMultiplex = optimize_multiplex.main
multipleOptimizations = multiple_run_optimization.main
CSVtoFASTA = CSVtoFASTA.main
addKeeplist2FASTA = add_keeplist_to_fasta.main



# SET DEPENDENCY PATH FOR MFEPRIMER
MFEprimer_PATH = glob.glob(os.path.dirname(__file__)+"/*mfeprimer*")[0]



__all__ = ["multiplex_wormhole", "panel_assessment", "setup_mfeprimer", "batch_primer3_design", "tabulate_dimers", "optimize_multiplex", "multiple_run_optimization", "CSVtoFasta", "add_keeplist_to_fasta", "plot_ASA_temps"]
#, "assessPanel", "primer3BatchDesign", "tabulateDimers", "plotASAtemps", "optimizeMultiplex", "multipleOptimizations", "CSVtoFASTA", "addKeeplist2FASTA", "MFEprimer_PATH"]