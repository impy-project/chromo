"""
Inspect nuclear remnants and placeholder records in a DPMJET h+A run.

Nuclear fragments and beam remnants never pass ``event.final_state()``:
that filter keeps only status 1 entries, and everything nuclear in the
event stack carries other status codes. This script prints what is
actually in the raw event stack for p + O at 100 TeV beam momentum and
how to select the relevant entries with numpy masks. See
doc/nuclear_fragments.md for the details.

Run with::

    python examples/extract_nuclear_fragments.py
"""

import numpy as np

from chromo.constants import GeV
from chromo.kinematics import FixedTarget, Momentum
from chromo.models import DpmjetIII307

kin = FixedTarget(Momentum(1e5 * GeV), "p", "O16")
run = DpmjetIII307(kin, seed=1)

for event in run(1):
    codes, counts = np.unique(event.status, return_counts=True)
    print("status codes in the raw event stack:")
    for c, n in zip(codes, counts):
        print(f"  status {c:5d}: {n} entries")
    print()

    # chromo stores the incoming beam particles at indices 0 and 1 with
    # status 4. For nuclear beams/targets this includes the beam nucleus.
    print("beam records (indices 0, 1):")
    for i in (0, 1):
        print(f"  index {i}: pdgid={event.pid[i]}  status={event.status[i]}")
    print()

    # DPMJET leaves the records of hadronization chains/strings in the
    # stack with PDG ID 99999 once their content has been converted into
    # final-state hadrons (these are NOT nuclear fragments, see
    # doc/nuclear_fragments.md). Those with status 2 stand for chains
    # which were directly emitted as one hadron, the ones with status
    # >= 100 are the fragmented strings:
    chain = event.pid == 99999
    direct = chain & (event.status == 2)
    fragmented = chain & (event.status >= 100)
    print(f"chain placeholders (pid 99999): {np.sum(chain)} total,")
    print(f"  {np.sum(direct)} collapsed to single hadrons (status 2),")
    print(
        f"  {np.sum(fragmented)} fragmented strings "
        f"(statuses {sorted({int(s) for s in event.status[fragmented]})})"
    )
    print()

    # The final state contains no nuclei at all:
    final = event.final_state()
    n_nuclei = np.sum(np.abs(final.pid) > 1000000000)
    print(f"final_state(): {len(final)} particles, nuclei among them: {n_nuclei}")
