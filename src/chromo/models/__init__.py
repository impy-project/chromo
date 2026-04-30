from chromo.models.dpmjetIII import DpmjetIII19, DpmjetIII307
from chromo.models.epos import EposLHC, EposLHCR, EposLHCRHadrRescattering
from chromo.models.fluka import Fluka
from chromo.models.fluka_decay import (
    STABLE_DEFAULT,
    DecayChainHandler,
    DecayChannel,
    DecayLine,
    FlukaDecay,
    FlukaIsotope,
)
from chromo.models.phojet import Phojet19, Phojet112
from chromo.models.pythia6 import Pythia6
from chromo.models.pythia8 import Pythia8, Pythia8Angantyr, Pythia8Cascade
from chromo.models.qgsjet import QGSJet01d, QGSJetII03, QGSJetII04, QGSJetIII
from chromo.models.sibyll import (
    Sibyll21,
    Sibyll23c,
    Sibyll23d,
    Sibyll23dStarMixed,
    Sibyll23e,
    Sibyll23eStarMixed,
)
from chromo.models.sophia import Sophia20
from chromo.models.urqmd import UrQMD34
from chromo.util import deprecated_aliases

__all__ = (
    "STABLE_DEFAULT",
    "DecayChainHandler",
    "DecayChannel",
    "DecayLine",
    "DpmjetIII19",
    "DpmjetIII307",
    "EposLHC",
    "EposLHCR",
    "EposLHCRHadrRescattering",
    "Fluka",
    "FlukaDecay",
    "FlukaIsotope",
    "Phojet19",
    "Phojet112",
    "Pythia6",
    "Pythia8",
    "Pythia8Angantyr",
    "Pythia8Cascade",
    "QGSJet01d",
    "QGSJetII03",
    "QGSJetII04",
    "QGSJetIII",
    "Sibyll21",
    "Sibyll23c",
    "Sibyll23d",
    "Sibyll23dStarMixed",
    "Sibyll23e",
    "Sibyll23eStarMixed",
    "Sophia20",
    "UrQMD34",
)

__getattr__ = deprecated_aliases(
    __name__, {"DpmjetIII193": "DpmjetIII19", "Phojet193": "Phojet19"}
)
