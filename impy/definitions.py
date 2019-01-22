'''This module initializes lists of namedtuple that link the definitions
of models and their various versions to the existing wrapper classes.
'''

from collections import namedtuple
from impy.models import (sibyll, dpmjetIII, epos, phojet)

# Objects of this type contain all default initialization directions
# for each interaction model and create dictionaries that link the
# these settings to either the 'tag' or the 'crmc_id'
InteractionModelDef = namedtuple('InteractionModelDef', [
    'tag', 'name', 'version', 'crmc_id', 'library_name', 'RunClass',
    'EventClass', 'output_frame'
])

# The list of available interaction models. Extend with new models or
# versions when available
interaction_model_nt_init = [[
    'SIBYLL23C', 'SIBYLL', '2.3c', -1, 'sib23c', sibyll.SIBYLLRun,
    sibyll.SibyllEvent, 'center-of-mass'
], [
    'SIBYLL23', 'SIBYLL', '2.3', -1, 'sib23', sibyll.SIBYLLRun,
    sibyll.SibyllEvent, 'center-of-mass'
], [
    'SIBYLL21', 'SIBYLL', '2.1', -1, 'sib21', sibyll.SIBYLLRun,
    sibyll.SibyllEvent, 'center-of-mass'
], [
    'DPMJETIII171', 'DPMJET-III', '17.1', -1, 'dpmjetIII171',
    dpmjetIII.DpmjetIIIRun, dpmjetIII.DpmjetIIIEvent, 'center-of-mass'
], [
    'DPMJETIII306', 'DPMJET-III', '3.0-6', -1, 'dpmjet306',
    dpmjetIII.DpmjetIIIRun, dpmjetIII.DpmjetIIIEvent, 'center-of-mass'
], [
    'EPOSLHC', 'EPOS', 'LHC', -1, 'eposlhc', epos.EPOSRun, epos.EPOSEvent,
    'center-of-mass'
], [
    'PHOJET112', 'PHOJET', '1.12-35', -1, 'dpmjet306', dpmjetIII.DpmjetIIIRun,
    dpmjetIII.DpmjetIIIEvent, 'center-of-mass'
], [
    'PHOJET171', 'PHOJET', '17.1', -1, 'dpmjetIII171', phojet.PHOJETRun,
    phojet.PhojetEvent, 'center-of-mass'
]]

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
