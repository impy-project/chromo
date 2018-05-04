import pyhepmc

def test_basic():
    #
    # In this example we will place the following event into HepMC "by hand"
    #
    #     name status pdg_id  parent Px       Py    Pz       Energy      Mass
    #  1  !p+!    3   2212    0,0    0.000    0.000 7000.000 7000.000    0.938
    #  3  !p+!    3   2212    0,0    0.000    0.000-7000.000 7000.000    0.938
    #=========================================================================
    #  2  !d!     3      1    1,1    0.750   -1.569   32.191   32.238    0.000
    #  4  !u~!    3     -2    2,2   -3.047  -19.000  -54.629   57.920    0.000
    #  5  !W-!    3    -24    1,2    1.517   -20.68  -20.605   85.925   80.799
    #  6  !gamma! 1     22    1,2   -3.813    0.113   -1.833    4.233    0.000
    #  7  !d!     1      1    5,5   -2.445   28.816    6.082   29.552    0.010
    #  8  !u~!    1     -2    5,5    3.962  -49.498  -26.687   56.373    0.006

    # now we build the graph, which will looks like
    #                       p7                         #
    # p1                   /                           #
    #   \v1__p3      p5---v4                           #
    #         \_v3_/       \                           #
    #         /    \        p8                         #
    #    v2__p4     \                                  #
    #   /            p6                                #
    # p2                                               #
    #                                                  #
    evt = pyhepmc.GenEvent(pyhepmc.Units.GEV, pyhepmc.Units.MM)

    #                           px      py        pz       e      pdgid status
    p1 = pyhepmc.GenParticle( ( 0.0,    0.0,   7000.0,  7000.0  ), 2212,  3 )
    p2 = pyhepmc.GenParticle( ( 0.750, -1.569,   32.191,  32.238),    1,  3 )
    p3 = pyhepmc.GenParticle( ( 0.0,    0.0,  -7000.0,  7000.0  ), 2212,  3 )
    p4 = pyhepmc.GenParticle( (-3.047,-19.0,    -54.629,  57.920),   -2,  3 )

    v1 = pyhepmc.GenVertex();
    v1.add_particle_in (p1)
    v1.add_particle_out(p2)
    evt.add_vertex(v1)

    # Set vertex status if needed
    v1.status = 4

    v2 = pyhepmc.GenVertex()
    v2.add_particle_in (p3)
    v2.add_particle_out(p4)
    evt.add_vertex(v2)

    v3 = pyhepmc.GenVertex()
    v3.add_particle_in(p2)
    v3.add_particle_in(p4)
    evt.add_vertex(v3)

    p5 = pyhepmc.GenParticle( (-3.813,  0.113, -1.833, 4.233),  22, 1 )
    p6 = pyhepmc.GenParticle( ( 1.517,-20.68, -20.605,85.925), -24, 3 )

    v3.add_particle_out(p5)
    v3.add_particle_out(p6)

    v4 = pyhepmc.GenVertex()
    v4.add_particle_in (p6)
    evt.add_vertex(v4)

    p7 = pyhepmc.GenParticle( (-2.445, 28.816,  6.082,29.552),  1, 1 )
    p8 = pyhepmc.GenParticle( ( 3.962,-49.498,-26.687,56.373), -2, 1 )

    v4.add_particle_out(p7)
    v4.add_particle_out(p8)

    oss = pyhepmc.ostringstream()
    f = pyhepmc.WriterAscii(oss)
    f.write_event(evt)
    f.close()

    assert oss.str() == """HepMC::Version 3.0.0
HepMC::IO_GenEvent-START_EVENT_LISTING
E 0 4 8
U GEV MM
P 1 0 2212 0.0000000000000000e+00 0.0000000000000000e+00 7.0000000000000000e+03 7.0000000000000000e+03 0.0000000000000000e+00 3
V -1 4 [1]
P 2 -1 1 7.5000000000000000e-01 -1.5690000000000000e+00 3.2191000000000003e+01 3.2238000000000000e+01 6.2465990744549081e-02 3
P 3 0 2212 0.0000000000000000e+00 0.0000000000000000e+00 -7.0000000000000000e+03 7.0000000000000000e+03 0.0000000000000000e+00 3
P 4 3 -2 -3.0470000000000002e+00 -1.9000000000000000e+01 -5.4628999999999998e+01 5.7920000000000002e+01 3.3845236001575724e-01 3
V -3 0 [2,4]
P 5 -3 22 -3.8130000000000002e+00 1.1300000000000000e-01 -1.8330000000000000e+00 4.2329999999999997e+00 8.1621075709617186e-02 1
P 6 -3 -24 1.5169999999999999e+00 -2.0680000000000000e+01 -2.0605000000000000e+01 8.5924999999999997e+01 8.0799603408680156e+01 3
P 7 6 1 -2.4449999999999998e+00 2.8815999999999999e+01 6.0819999999999999e+00 2.9552000000000000e+01 -9.9503768772913739e-02 1
P 8 6 -2 3.9620000000000002e+00 -4.9497999999999998e+01 -2.6687000000000001e+01 5.6372999999999998e+01 -1.7403447934355551e-01 1
HepMC::IO_GenEvent-END_EVENT_LISTING

"""


def test_sequence_access():
    evt = pyhepmc.GenEvent()
    evt.particles = (pyhepmc.GenParticle(),)
    evt.particles[0].momentum = (1, 2, 3, 4)
    evt.particles[0].pid = 5
    evt.vertices = (pyhepmc.GenVertex(),)
    evt.vertices[0].position = (1, 2, 3, 4)
    assert evt.particles == [pyhepmc.GenParticle((1, 2, 3, 4), 5),]
    assert evt.vertices == [pyhepmc.GenVertex((1, 2, 3, 4)),]


def test_read_write():
    evt1 = pyhepmc.GenEvent()
    p1 = pyhepmc.GenParticle( (1, 2, 3, 4), 1, 1)
    p2 = pyhepmc.GenParticle( (5, 6, 7, 8), 2, 2)
    p3 = pyhepmc.GenParticle( (9, 10, 11, 12), 3, 3)
    p4 = pyhepmc.GenParticle( (13, 14, 15, 16), 4, 4)
    v1 = pyhepmc.GenVertex()
    v1.add_particle_in(p1)
    v1.add_particle_in(p2)
    v1.add_particle_out(p3)
    v1.add_particle_out(p4)
    evt1.add_vertex(v1)

    f = pyhepmc.WriterAscii("testfile.dat")
    f.write_event(evt1)
    f.close()

    f = pyhepmc.ReaderAscii("testfile.dat")
    evt2 = pyhepmc.GenEvent()

    assert evt1 != evt2
    f.read_event(evt2)
    f.close()

    assert evt1.particles == evt2.particles
    assert evt1.vertices == evt2.vertices
    assert evt1 == evt2


if __name__ == '__main__':
    test_basic()
    test_sequence_access()
    test_read_write()