import impy
from impy.kinematics import CompositeTarget
from impy.kinematics import EventKinematics
from impy.constants import TeV, GeV


def recognize_particle_input_type(arg):
    if isinstance(arg, int):
        return "pdg_id"
    elif isinstance(arg, str):
        return "string"
    elif isinstance(arg, tuple):
        if len(arg) != 2:
            raise ValueError(
                "tuple (A, Z) should have len == 2, but given = {0}".format(arg)
            )
        if not isinstance(arg[0], int):
            raise ValueError(
                "1st entry should be 'int' type, but it is {0} = {1}".format(
                    type(arg[0]), arg[0]
                )
            )
        if not isinstance(arg[1], int):
            raise ValueError(
                "1st entry should be 'int' type, but it is {0} = {1}".format(
                    type(arg[1]), arg[1]
                )
            )
        return "tuple"
    elif isinstance(arg, CompositeTarget):
        return "composite_target"
    else:
        raise ValueError("Unmaintained parameter type {0} = {1}".format(type(arg), arg))


air = CompositeTarget("Air")

air.add_component("Nitrogen", 14, 7, 2 * 0.78084)
air.add_component("Oxygen", 16, 8, 2 * 0.20946)
# air.add_component("Argon", 40, 18, 0.00934)
air.add_component("Oxygen(Vapor)", 16, 8, 0.004)
air.add_component("Hydrogen(Vapor)", 1, 1, 2 * 0.004)

print(recognize_particle_input_type(2121))
print(recognize_particle_input_type("proton"))
print(recognize_particle_input_type((5, 5)))
print(recognize_particle_input_type(air))


impy.impy_config["user_frame"] = "laboratory"

ekin = impy.kinematics.CompositeTargetKinematics(
    elab=13 * TeV, p1pdg=2212, composite_target=air
)

# # ekin = EventKinematics(elab = 13 * TeV, p1pdg = 2212,
# #                        nuc2_prop = air.get_random_AZ())

# generator = impy.models.DpmjetIII191(ekin)
# # generator = impy.models.Sibyll23d(ekin)

# for event in generator(1000):
#     print("Nevent = ", generator.nevents)
