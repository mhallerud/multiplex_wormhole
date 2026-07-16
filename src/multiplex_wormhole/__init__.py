# LOAD ALL SUBMODULES
import sys, os, importlib, glob
sys.path.append(os.path.dirname(__file__))

from . import multiplexWormhole
from . import panel_assessment
from . import batch_primer3_design
from . import tabulate_dimers
from . import optimize_multiplex 
from . import multiple_run_optimization 
from . import offtarget_thermodynamics
from .helpers import _setup_mfeprimer
from .helpers import CSVtoFasta
from .helpers import add_keeplist_to_fasta
from . import plot_ASA_temps as plot_ASA_temps
#plot_ASA_temps = importlib.import_module("./plot_ASA_temps")



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
offtargetThermodynamics = offtarget_thermodynamics.main
setup_mfeprimer = _setup_mfeprimer.main


# SET DEPENDENCY PATH FOR MFEPRIMER
_hits = glob.glob(os.path.join(os.path.dirname(__file__), "*mfeprimer*"))
MFEprimer_PATH = _hits[0] if _hits else None



__all__ = ["multiplexWormhole", "panel_assessment", "setup_mfeprimer", "batch_primer3_design", 
           "tabulate_dimers", "optimize_multiplex", "multiple_run_optimization", "CSVtoFasta", 
           "add_keeplist_to_fasta", "plot_ASA_temps", "offtarget_thermodynamics"]
