from chromo.models.dpmjetIII import DpmjetIII19_DEV
from chromo.models.sibyll import (
    Sibyll23c00,
    Sibyll23c01,
    Sibyll23c02,
    Sibyll23c03,
    Sibyll23dStarBar,
    Sibyll23dStarNoEnh,
    Sibyll23dStarRho,
    Sibyll23dStarStrange,
    Sibyll23e_DEV,
    SibyllRS_DEV,
    SibyllNext,
    SibyllNextFullyResolvedHMD,
    SibyllNextNoExcitationNucleonScaling,
    Sibyll23eStarBar,
    Sibyll23eStarNoEnh,
    Sibyll23eStarRho,
    Sibyll23eStarStrange,
)
from chromo.util import deprecated_aliases

__all__ = (
    "DpmjetIII19_DEV",
    "Sibyll23c00",
    "Sibyll23c01",
    "Sibyll23c02",
    "Sibyll23c03",
    "Sibyll23dStarBar",
    "Sibyll23dStarNoEnh",
    "Sibyll23dStarRho",
    "Sibyll23dStarStrange",
    "Sibyll23e_DEV",
    "SibyllRS_DEV",
    "SibyllNext",
    "SibyllNextFullyResolvedHMD",
    "SibyllNextNoExcitationNucleonScaling",
    "Sibyll23eStarBar",
    "Sibyll23eStarNoEnh",
    "Sibyll23eStarRho",
    "Sibyll23eStarStrange",
)

__getattr__ = deprecated_aliases(__name__, {"DpmjetIII193_DEV": "DpmjetIII19_DEV"})
