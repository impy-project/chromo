"""
Extract nuclear beam records and nucleon-level remnants from a DPMJET h+A run.

Nuclear fragments and beam remnants never pass ``event.final_state()``: that
filter keeps only status 1 entries. ``event.final_state_with_nucl_frag()``
additionally returns the nucleus records (chromo status 4, proper 10LZZZAAAI
PDG codes) and the nucleon-level remnants of the intranuclear cascade
(chromo status 5). This script shows the selection on p + Pb events, where
the cascade bookkeeping is abundant. See doc/nuclear_fragments.md for the
generator-by-generator details and the caveats.

Run with::

    python examples/extract_nuclear_fragments.py
"""

from collections import Counter

from chromo.kinematics import FixedTarget
from chromo.models import DpmjetIII307

kin = FixedTarget(100, "p", "Pb208")
run = DpmjetIII307(kin, seed=1)

for event in run(1):
    codes = Counter(event.status)
    print("status codes in the raw event stack:")
    for c, n in sorted(codes.items()):
        print(f"  status {c:5d}: {n} entries")
    print()

    # chromo normalizes every generator: the incoming beam particles are
    # records 0 and 1 with status 4 (for a nuclear target this is the beam
    # nucleus with a proper PDG code), nucleon-level remnants get status 5.
    with_frags = event.final_state_with_nucl_frag()
    nuclei = with_frags[with_frags.status == 4]
    remnants = with_frags[with_frags.status == 5]
    print(f"final_state(): {len(event.final_state())} hadrons")
    print(
        f"final_state_with_nucl_frag(): {len(with_frags)} records = "
        f"{len(with_frags) - len(nuclei) - len(remnants)} hadrons + "
        f"{len(nuclei)} nucleus records + {len(remnants)} remnant nucleons"
    )
    print("nucleus records (status 4):")
    for pid in nuclei.pid:
        pid = int(pid)
        if abs(pid) > 1000000000:
            print(
                f"  pdgid={pid}  A={abs(pid) // 10 % 1000}  Z={abs(pid) // 10000 % 10000}"
            )
        else:
            print(f"  pdgid={pid}  (beam nucleon)")
    print(f"remnant nucleons (status 5): {Counter(int(p) for p in remnants.pid)}")
