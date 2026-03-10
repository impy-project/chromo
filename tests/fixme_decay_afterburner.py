import os
import sys
from collections import Counter

from chromo import chromo_config
from chromo.constants import GeV
from chromo.kinematics import EventKinematics
from chromo.models import DpmjetIII191

# This class will go through the event and decay all particles that should be
# unstable but did not decay in some other generator


# FIXME this should not be defined in a test
# if it works, it should be part of the library and tested here
class Pythia8DecayAfterburner:
    def __init__(self, stable_list):
        self.stable_list = stable_list
        self._init_pythia()

    def _init_pythia(self):
        # FIXME Anatoly test cannot contain special paths
        pythia_dir = os.path.join(
            os.path.expanduser("~"), "devel", "apps", "pythia8240"
        )
        sys.path.append(os.path.join(pythia_dir, "lib"))
        import pythia8

        self.pythia = pythia8.Pythia()
        self.pythia.readString("ProcessLevel:all = off")
        self.pythia.readString("ParticleDecays:tau0Max = 1e100")
        # Set muons unstable
        self.pythia.particleData.mayDecay(13, True)
        self.pythia.init()

    def __call__(self, event):
        nappend = 0
        pid, status, px, py, pz, en, m = [], [], [], [], [], [], []

        for ip in range(event.npart):
            if event.status[ip] != 1 or abs(event.pid[ip]) in self.stable_list:
                continue

            # Set particle to not final and simulate decay
            event.status[ip] = 2
            self.pythia.event.reset()
            # put decaying particle
            self.pythia.event.append(
                int(event.pid[ip]),
                91,
                0,
                0,
                event.px[ip],
                event.py[ip],
                event.pz[ip],
                event.en[ip],
                event.m[ip],
            )
            self.pythia.particleData.mayDecay(int(event.pid[ip]), True)
            # Decay it
            self.pythia.forceHadronLevel()
            for p in self.pythia.event:
                if not p.isFinal():
                    continue
                pid.append(p.id())
                status.append(1)
                px.append(p.px())
                py.append(p.py())
                pz.append(p.pz())
                en.append(p.e())
                m.append(p.m())
                nappend += 1

        append = slice(event.npart, event.npart + nappend)

        event.en[append] = en
        event.pid[append] = pid
        event.status[append] = status
        event.px[append] = px
        event.py[append] = py
        event.pz[append] = pz
        event.m[append] = m

        event.npart = nappend
        event.selection = slice(None, event.npart)
        event._apply_slicing()


def test_decay_afterburner():
    event_kinematics = EventKinematics(
        ecm=200 * GeV,
        particle1=2212,
        particle2=2212,
        # particle2=(14,7)
    )

    # Watch out this setting!
    chromo_config["pre_slice"] = False

    # The rest is pretty standard
    generator = DpmjetIII191(event_kinematics)

    # Here provide the list of particles which you want to retain as stable
    pythia_afterburner = Pythia8DecayAfterburner(
        stable_list=[2212, 11, 12, 14, 15, 16, 22]
    )

    before = Counter()
    after = Counter()
    for event in generator(200):
        before.update(event.pid)
        # This has to be the first call after an event is generated. The event object
        # will be modified and finalized by this call
        pythia_afterburner(event)
        # Here filter only the particles that are remaining as stable
        event.filter_final_state()
        # Enjoy the result
        after.update(event.pid)
    assert len(after) < len(before)
    # TODO check list of specific particles that are expected to
    # be decayed by Pythia but not Dpmjet
