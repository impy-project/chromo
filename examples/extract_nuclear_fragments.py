"""
Select nuclear remnants in p + Pb collisions with DPMJET-III.

See doc/nuclear_fragments.md for the status codes and generator support.
"""

from collections import Counter

from chromo.kinematics import FixedTarget
from chromo.models import DpmjetIII307
from chromo.util import pdg2AZ

kin = FixedTarget(100, "p", "Pb208")
run = DpmjetIII307(kin, seed=1)

for event in run(1):
    frags = event.final_state_with_nucl_frag()
    nuclei = frags[frags.status == 4]
    nucleons = frags[frags.status == 5]

    print(f"final state:         {len(event.final_state())}")
    print(f"with remnants:       {len(frags)}")
    print(f"residual nuclei:     {[pdg2AZ(pid) for pid in nuclei.pid]}")
    print(f"spectator nucleons:  {dict(Counter(int(p) for p in nucleons.pid))}")
