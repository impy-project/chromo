from impy.kinematics import CenterOfMass
from impy.constants import TeV
import impy.models as im


evt_kin = CenterOfMass(13 * TeV, "gamma", "proton")
generator = im.Sophia20(evt_kin)

for event in generator(10):
    print(event.pid, event.en)
