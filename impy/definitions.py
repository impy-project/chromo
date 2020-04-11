'''This module initializes lists of namedtuple that link the definitions
of models and their various versions to the existing wrapper classes.
'''

from collections import namedtuple
from impy.models import (sibyll, dpmjetIII, epos, phojet, urqmd, pythia6,
                         pythia8, qgsjet)

# Objects of this type contain all default initialization directions
# for each interaction model and create dictionaries that link the
# these settings to either the 'tag' or the 'crmc_id'
InteractionModelDef = namedtuple('InteractionModelDef', [
    'tag', 'name', 'version', 'crmc_id', 'library_name', 'RunClass',
    'EventClass', 'output_frame'
])

# The list of available interaction models. Extend with new models or
# versions when available
interaction_model_nt_init = [
    [
        'SIBYLL23C', 'SIBYLL', '2.3c03', -1, 'sib23c03', sibyll.SIBYLLRun,
        sibyll.SibyllEvent, 'center-of-mass'
    ],
    [
        'SIBYLL23C00', 'SIBYLL', '2.3c00', -1, 'sib23c00', sibyll.SIBYLLRun,
        sibyll.SibyllEvent, 'center-of-mass'
    ],
    [
        'SIBYLL23C01', 'SIBYLL', '2.3c01', -1, 'sib23c01', sibyll.SIBYLLRun,
        sibyll.SibyllEvent, 'center-of-mass'
    ],
    [
        'SIBYLL23C02', 'SIBYLL', '2.3c02', -1, 'sib23c02', sibyll.SIBYLLRun,
        sibyll.SibyllEvent, 'center-of-mass'
    ],
    [
        'SIBYLL23C03', 'SIBYLL', '2.3c03', -1, 'sib23c03', sibyll.SIBYLLRun,
        sibyll.SibyllEvent, 'center-of-mass'
    ],
    [
        'SIBYLL23C04', 'SIBYLL', '2.3c04', -1, 'sib23c04', sibyll.SIBYLLRun,
        sibyll.SibyllEvent, 'center-of-mass'
    ],
    [
        'SIBYLL23', 'SIBYLL', '2.3', -1, 'sib23', sibyll.SIBYLLRun,
        sibyll.SibyllEvent, 'center-of-mass'
    ],
    [
        'SIBYLL21', 'SIBYLL', '2.1', -1, 'sib21', sibyll.SIBYLLRun,
        sibyll.SibyllEvent, 'center-of-mass'
    ],
    [
        'DPMJETIII191', 'DPMJET-III', '19.1', -1, 'dpmjetIII191',
        dpmjetIII.DpmjetIIIRun, dpmjetIII.DpmjetIIIEvent, 'center-of-mass'
    ],
    [
        'DPMJETIII306', 'DPMJET-III', '3.0-6', -1, 'dpmjet306',
        dpmjetIII.DpmjetIIIRun, dpmjetIII.DpmjetIIIEvent, 'center-of-mass'
    ],
    [
        'EPOSLHC', 'EPOS', 'LHC', -1, 'eposlhc', epos.EPOSRun, epos.EPOSEvent,
        'center-of-mass'
    ],
    [
        'PHOJET112', 'PHOJET', '1.12-35', -1, 'dpmjet306',
        phojet.PHOJETRun, phojet.PhojetEvent, 'center-of-mass'
    ],
    [
        'PHOJET191', 'PHOJET', '19.1', -1, 'dpmjetIII191', phojet.PHOJETRun,
        phojet.PhojetEvent, 'center-of-mass'
    ],
    [
        'URQMD34', 'UrQMD', '3.4', -1, 'urqmd34', urqmd.UrQMDRun,
        urqmd.UrQMDEvent, 'center-of-mass'
    ],
    [
        'PYTHIA6', 'Pythia', '6.428', -1, 'pythia6', pythia6.PYTHIA6Run,
        pythia6.PYTHIA6Event, 'center-of-mass'
    ],
    [
        'PYTHIA8', 'Pythia', '8.240', -1, 'pythia8', pythia8.PYTHIA8Run,
        pythia8.PYTHIA8Event, 'center-of-mass'
    ],
    [
        'QGSJET01C', 'QGSJet', '01c', -1, 'qgs01', qgsjet.QGSJet01Run,
        qgsjet.QGSJETEvent, 'laboratory'
    ],
    [
        'QGSJETII03', 'QGSJet', 'II-03', -1, 'qgsII03', qgsjet.QGSJetIIRun,
        qgsjet.QGSJETEvent, 'laboratory'
    ],
    [
        'QGSJETII04', 'QGSJet', 'II-04', -1, 'qgsII04', qgsjet.QGSJetIIRun,
        qgsjet.QGSJETEvent, 'laboratory'
    ]
]

# Different kinds of lookup-maps/dictionaries to refer to models by name
# or ID
interaction_model_by_tag = {}
interaction_model_by_crmc_id = {}

for arg_tup in interaction_model_nt_init:
    nt = InteractionModelDef(*arg_tup)
    interaction_model_by_tag[nt.tag] = nt
    interaction_model_by_crmc_id[nt.crmc_id] = nt


def make_generator_instance(int_model_def):
    """Returns instance of a <Model>Run.
    """
    return int_model_def.RunClass(int_model_def)
