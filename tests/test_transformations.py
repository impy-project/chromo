import numpy as np
from impy.kinematics import (
    matrix_rotate_back_hep,
    matrix_lorentz_boost_hep,
    from_cms,
    from_laboratory,
)


# random unit vector
def get_random_3vect():
    v = []
    for _ in range(3):
        v.append(np.random.random())
    v = np.array(v)
    v = v / np.sqrt(np.sum(v[:] ** 2))
    return v


# random 4-momentum
def p_vec(energy, mass):
    n3 = get_random_3vect()
    p3 = np.sqrt((energy - mass) * (energy + mass)) * n3
    return np.array([*p3, energy])


# Create npart random 4-momenta
def create_ps(npart, mass):
    pcol = []
    for _ in range(npart):
        en = np.random.uniform(1e2, 1e3)
        pcol.append(p_vec(en, mass))
    return np.array(pcol).transpose()


# Get mass
def get_mass(p):
    pmod = np.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
    emod = p[3]
    return np.sqrt((emod - pmod) * (emod + pmod))


def test_rotate_back():
    """Test matrix_rotate_back_hep()"""
    for _ in range(1000):
        # Get arbitrary 3-vector with module vec_mod
        vec_mod = 235
        vec = vec_mod * get_random_3vect()
        vec4 = np.array([*vec[:], 0], dtype=np.float64)
        # Transform 3-vector with module vec_mod along Z
        # using obtrained matrix
        vec4z = np.array([0, 0, vec_mod, 0], dtype=np.float64)
        vec4z_back = np.dot(matrix_rotate_back_hep(vec4), vec4z)

        for i in range(3):
            rel_error = abs(vec4z_back[i] - vec4[i]) / vec4[i]
            assert (
                rel_error < 1e-15
            ), "Too large error = {0} in back transformation for i = {1}".format(
                rel_error, i
            )


def test_lorentz_boost():
    for _ in range(1000):
        energy = 1e3  # GeV
        mass = 0.938
        # Random 4-momentum
        p4 = p_vec(energy, mass)
        # Transform to the frame with
        gamma = 1000
        n_direction = get_random_3vect()

        # Lorentz boost
        p4_new = np.dot(matrix_lorentz_boost_hep(gamma, n_direction), p4)
        # Transform back
        p4_trans = np.dot(matrix_lorentz_boost_hep(gamma, -n_direction), p4_new)

        for i in range(4):
            rel_error = abs(p4_trans[i] - p4[i]) / p4[i]
            tolerance_level = 1e-15 * gamma**2 * (energy / mass) ** 2
            assert (
                rel_error < tolerance_level
            ), "Error = {0} > tolerance level {2} in back transformation for i = {1}".format(
                rel_error, i, tolerance_level
            )


def test_from_cms():
    for _ in range(1000):
        mass1 = 0.938
        mass2 = 0.938
        energy1 = 1e3
        energy2 = 1e2
        p1 = p_vec(energy1, mass1)
        p2 = p_vec(energy2, mass2)

        s = p1 + p2
        ecm = np.sqrt(s[3] ** 2 - np.sum(s[:3] ** 2))
        s2 = ecm**2
        p_cms = np.sqrt((s2 - (mass1 + mass2) ** 2) * (s2 - (mass1 - mass2) ** 2)) / (
            2 * ecm
        )
        e_cms = np.sqrt(p_cms**2 + mass1**2)
        # Along z direction
        p1_cms = np.array([0, 0, p_cms, e_cms], dtype=np.float64)
        p2_cms = np.array([0, 0, -p_cms, e_cms], dtype=np.float64)
        # Get set of vectors[4, :]
        p_array = np.array([p1_cms, p2_cms]).transpose()
        p_original = from_cms(p1, p2, p_array).transpose()
        p1_trans = p_original[0]
        p2_trans = p_original[1]

        gamma = (energy1 + energy2) / ecm
        tolerance_level = 1e-15 * ((energy1 / mass1) * gamma) ** 2

        for i in range(4):
            rel_error1 = abs(p1_trans[i] - p1[i]) / p1[i]
            rel_error2 = abs(p2_trans[i] - p2[i]) / p2[i]
            assert (
                rel_error1 < tolerance_level
            ), "Error1 = {0} > tolerance level {2} in back transformation for i = {1}".format(
                rel_error1, i, tolerance_level
            )
            assert (
                rel_error2 < tolerance_level
            ), "Error2 = {0} > tolerance level {2} in back transformation for i = {1}".format(
                rel_error2, i, tolerance_level
            )


def test_from_laboratory():
    for _ in range(1000):
        mass1 = 0.938
        mass2 = 0.938
        energy1 = 1e3
        energy2 = 1e2
        p1 = p_vec(energy1, mass1)
        p2 = p_vec(energy2, mass2)

        s = p1 + p2
        ecm = np.sqrt(s[3] ** 2 - np.sum(s[:3] ** 2))
        elab = (ecm**2 - mass1**2 - mass2**2) / (2 * mass2)
        plab = np.sqrt((elab - mass1) * (elab + mass1))
        # Along z direction
        p1_rest = np.array([0, 0, plab, elab], dtype=np.float64)
        p2_rest = np.array([0, 0, 0, mass2], dtype=np.float64)
        # Get set of vectors[4, :]
        p_array = np.array([p1_rest, p2_rest]).transpose()
        p_original = from_laboratory(p1, p2, p_array).transpose()
        p1_trans = p_original[0]
        p2_trans = p_original[1]

        gamma = energy2 / mass2
        tolerance_level = 1e-15 * ((energy1 / mass1) * gamma) ** 2

        for i in range(4):
            rel_error1 = abs(p1_trans[i] - p1[i]) / p1[i]
            rel_error2 = abs(p2_trans[i] - p2[i]) / p2[i]
            assert (
                rel_error1 < tolerance_level
            ), "Error1 = {0} > tolerance level {2} in back transformation for i = {1}".format(
                rel_error1, i, tolerance_level
            )
            assert (
                rel_error2 < tolerance_level
            ), "Error2 = {0} > tolerance level {2} in back transformation for i = {1}".format(
                rel_error2, i, tolerance_level
            )
