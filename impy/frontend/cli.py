"""
A commandline interface that mimics CRMC to ease transition from CRMC.

CRMC (Cosmic Ray Monte Carlo package) https://web.ikp.kit.edu/rulrich/crmc.html
"""

import argparse
from collections import namedtuple

ModelData = namedtuple("ModelData", "name module lib run_class evt_class")

models = {
    # numbers must match definition in CRMC
    0: ModelData("EPOS-LHC", "epos", "eposlhc", "EPOSMCEvent", "EPOSMCRun"),
    2: ModelData("QGSJet01", "qgsjet", "qgs01", "QGSJet01Run", "QGSJet01Event"),
    4: ModelData("Pythia6", "pythia6", "pythia6", "PYTHAMCRun", "PYTHAEvent"),
    6: ModelData("SIBYLL-2.3", "sibyll", "sib23", "SibyllMCRun", "SibyllMCEvent"),
    14: ModelData("SIBYLL-2.3c", "sibyll", "sib23", "SibyllMCRun", "SibyllMCEvent"),
    11: ModelData("QGSJetII.03", "gqsjet", "QGSJetIIMCRun", "QGSJetIIEvent"),
    7: ModelData("QGSJetII.04", "qgsjet", "QGSJetIIMCRun", "QGSJetIIEvent"),
    13: ModelData("DPMJetII.55", "dpmjetII", "DpmjetIIMCRun", "DpmjetIIMCEvent"),
}

parser = argparse.ArgumentParser()
parser.add_argument("-v,--version", action="version", version="0.1")
parser.add
parser.add_argument("-m,--model", type=int, default=0,
    help="model [0=EPOS_LHC, 1=EPOS_1.99, 2=QGSJET01,"
         "4=Pythia_6.115, 6=Sibyll_2.3, 7=QGSJETII-04,"
         "12=DPMJet 3.0-6]")
