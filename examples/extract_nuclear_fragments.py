"""
Extract nuclear fragments and remnants from a DPMJET h+A run.

Nuclear fragments and beam remnants never pass ``event.final_state()``:
that filter keeps only status 1 entries, and nuclei in the event stack
carry other status codes. This script prints what is actually in the raw
event stack for p + O at 100 TeV beam momentum and shows how to select
the interesting entries with numpy masks. See doc/nuclear_fragments.md
for the details.

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

    # DPMJET stores nuclear fragments as placeholder entries with PDG ID
    # 99999. Entries with status 2 are hadronization chains that are
    # already absorbed into the final state; the nuclear prefragments
    # carry status codes >= 100 (DPMJET uses 103, 104, 106, and higher
    # 5xx/6xx/7xx codes for these):
    prefrag = (event.pid == 99999) & (event.status >= 100)
    print(f"nuclear prefragments (pid 99999, status >= 100): {np.sum(prefrag)}")
    for i in np.where(prefrag)[0]:
        print(
            f"  status={event.status[i]}  m={event.m[i]:.1f} GeV"
            f"  E={event.en[i]:.1f} GeV"
        )
    print()

    # The final state contains no nuclei at all:
    final = event.final_state()
    n_nuclei = np.sum(np.abs(final.pid) > 1000000000)
    print(f"final_state(): {len(final)} particles, nuclei among them: {n_nuclei}")
