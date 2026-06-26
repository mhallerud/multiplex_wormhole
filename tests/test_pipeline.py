# tests/test_pipeline.py
import os
import sys
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src/multiplex_wormhole"))
from multiplexWormhole import main as multiplexWormhole
from helpers import setup_mfeprimer
from panel_assessment import main as assessPanel
EXAMPLES = os.path.join(os.path.dirname(__file__), "../examples/Input_Templates.csv")
FASTA = os.path.join(os.path.dirname(__file__), "../examples/example_keeplist.fasta")

#import shutil
#HAS_MFEPRIMER = shutil.which("mfeprimer") is not None
HAS_MFEPRIMER = setup_mfeprimer.main() is not None



@pytest.mark.skipif(not HAS_MFEPRIMER, reason="MFEprimer binary not installed")
def test_multiplex_wormhole(tmp_path, templates=EXAMPLES):
    # test multiplex wormhole - primer design pipeline
    multiplexWormhole(TEMPLATES=templates, N_LOCI=30, OUTDIR=tmp_path, PREFIX='pytest', 
                      KEEPLIST_FA=None, N_RUNS=1, SIMPLE=30)
    primer_out = os.path.join(tmp_path, "1_PrimerDesign/FilteredPrimers.fa")
    dimer_out = os.path.join(tmp_path, "2_PredictedDimers/MFEprimerDimers.txt")
    table_out = os.path.join(tmp_path, "2_PredictedDimers/PrimerPairInteractions_wide.csv")
    opt_out = os.path.join(tmp_path, "3_OptimizedMultiplexes/pytest_Run01_primers.fasta")
    assert os.path.exists(primer_out), "Batch primer design failed!"
    assert os.path.exists(dimer_out), "MFEprimer dimer failed!"
    assert os.path.exists(table_out), "Tabulate dimers failed!"
    assert os.path.exists(opt_out), "Optimization failed!"

    # test multiple run summary
    #results = pd.read_csv(os.path.join(tmp_path, "3_OptimizedMultiplexed/pytest_RunSummary.csv"))
    #assert len(results) > 0, "No primer pairs in output"
    
    # test multiplex wormhole - primer design with keeplist
    multiplexWormhole(TEMPLATES=templates, N_LOCI=40, OUTDIR=tmp_path, PREFIX='pytest_keeplist',
                      KEEPLIST_FA=opt_out, N_RUNS=1, SIMPLE=30)
    keeplist_out = os.path.join(tmp_path, "3_OptimizedMultiplexes/pytest_keeplist_Run01_primers.fasta")
    assert os.path.exists(keeplist_out), "Optimization with keeplist failed!"
    
    # clean up
    os.remove(os.path.join(os.path.dirname(__file__), "../multiplex_wormhole_pytest.log"))
    os.remove(os.path.join(os.path.dirname(__file__), "../multiplex_wormhole_pytest_keeplist.log"))
    


def test_panel_assessment(fasta=FASTA):
    counts = assessPanel(PRIMERS=fasta)
    assert counts, "Panel assessment failed!"
    # clean up
    os.remove(os.path.join(os.path.dirname(__file__), "../examples/example_keeplist_MFEdimers.txt"))
    os.remove(os.path.join(os.path.dirname(__file__), "../examples/example_keeplist_MFEdimers_ends.txt"))
    os.remove(os.path.join(os.path.dirname(__file__), "../examples/example_keeplist_PrimerPairDeltaG.log"))
    os.remove(os.path.join(os.path.dirname(__file__), "../examples/example_keeplist_PrimerPairDeltaG_mean.csv"))
    os.remove(os.path.join(os.path.dirname(__file__), "../examples/example_keeplist_PrimerPairDeltaG_wide.csv"))
    os.remove(os.path.join(os.path.dirname(__file__), "../examples/example_keeplist_PrimerPairDimers.log"))
    os.remove(os.path.join(os.path.dirname(__file__), "../examples/example_keeplist_PrimerPairDimers_binary_sum.csv"))
    os.remove(os.path.join(os.path.dirname(__file__), "../examples/example_keeplist_PrimerPairDimers_binary_wide.csv"))
    os.remove(os.path.join(os.path.dirname(__file__), "../examples/example_keeplist_PrimerPairDimers_sum.csv"))
    os.remove(os.path.join(os.path.dirname(__file__), "../examples/example_keeplist_PrimerPairDimers_wide.csv"))
    os.remove(os.path.join(os.path.dirname(__file__), ".log"))


