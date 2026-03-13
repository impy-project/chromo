"""Quick verification script for Angantyr with precomputed initialization files.

Usage:
    python scripts/test_angantyr.py
"""

import sys
from pathlib import Path

# Ensure chromo is importable from the source tree
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from chromo.kinematics import CenterOfMass
from chromo.models import Pythia8Angantyr
from chromo.util import name2pdg


def main():
    print("Testing Pythia8Angantyr with p + O16 at 100 GeV CMS ...")
    kin = CenterOfMass(100, "p", "O16")
    gen = Pythia8Angantyr(kin, seed=42)

    n_events = 5
    n_ok = 0
    for event in gen(n_events):
        hi = gen._pythia.info.hiInfo
        b = hi.b if hi is not None else float("nan")
        npp = hi.nPartProj if hi is not None else 0
        npt = hi.nPartTarg if hi is not None else 0
        fs = event.final_state_charged()
        print(
            f"  event: b={b:.2f} fm, nPartProj={npp}, nPartTarg={npt}, "
            f"nFinalCharged={len(fs)}"
        )
        assert hi is not None, "hiInfo is None — Angantyr not active"
        n_ok += 1

    print(f"\nAll {n_ok}/{n_events} events OK.")

    cs = gen.cross_section()
    print(f"\nCross sections for p+O16 at 100 GeV CMS:")
    print(f"  total    = {cs.total:.1f} mb")
    print(f"  inelastic= {cs.inelastic:.1f} mb")
    print(f"  elastic  = {cs.elastic:.1f} mb")

    print("\nVerification PASSED.")


if __name__ == "__main__":
    main()
