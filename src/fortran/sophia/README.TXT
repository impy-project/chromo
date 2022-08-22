

                           SOPHIA 2.0
                           ==========



A.M"ucke, Ralph Engel, J.P.Rachen, R.J.Protheroe, and Todor Stanev

(astro-ph/9903478, to appear in Comp.Phys.Commun.)


Files:
-----


SOPHIA20.f           (main program)

eventgen.f           (hadronic event generator)

sampling.f           (sampling of photon energies)

inpoutput.f          (histogramming, input/output)

jetset74dp.f         (Lund fragmentation program, changed to double prec.)

rndm.f               (random number generator)

testPL.inp            (sample input for straight power law calculation)

testBP.inp            (sample input for broken power law calculation)

testBB.inp            (sample input for black body radiation calculation)

README.TXT           (this file)





Particle record:
---------------


After event generation the final state particles are found in 

       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb

NP is the number of particles produced in the last event.
The particle momenta of particle I in x, y, z directions are 
P(I,1), P(I,2), P(I,3) and P(I,4) is the energy (in GeV). 
The particle mass is stored in P(I,5).

The particle identity is given by LLIST with the following numbering
scheme:

 -----------------------
 code particle     mass
 -----------------------
   1  gam          0.0000
   2  e+           0.0005
   3  e-           0.0005
   4  mu+          0.1057
   5  mu-          0.1057
   6  pi0          0.1350
   7  pi+          0.1396
   8  pi-          0.1396
   9  k+           0.4937
  10  k-           0.4937
  11  k0l          0.4977
  12  k0s          0.4977
  13  p            0.9383
  14  n            0.9396
  15  nue          0.0000
  16  nueb         0.0000
  17  num          0.0000
  18  numb         0.0000
  21  k0           0.4977
  22  k0b          0.4977
  23  eta          0.5488
  24  etap         0.9576
  25  rho+         0.7714
  26  rho-         0.7714
  27  rho0         0.7717
  28  k*+          0.8921
  29  k*-          0.8921
  30  k*0          0.8965
  31  k*0b         0.8965
  32  omeg         0.7826
  33  phi          1.1020
  34  SIG+         1.1894
  35  SIG0         1.1925
  36  SIG-         1.1973
  37  XI0          1.3149
  38  XI-          1.3213
  39  LAM          1.1156
  40  DELT++       1.2300
  41  DELT+        1.2310
  42  DELT0        1.2320
  43  DELT-        1.2330
  44  SIG*+        1.3828
  45  SIG*0        1.3837
  46  SIG*-        1.3872
  47  XI*0         1.5318
  48  XI*-         1.5350
  49  OME*-        1.6724

Antibaryons have negative ID numbers correspondingly.
Decayed particles are marked by adding 10000 to their ID code.

The decay of a particle with code ID can be turned off via

      IDB(ID) = -ABS(IDB(ID))

and turned on via

      IDB(ID) = ABS(IDB(ID))

This has to be done only once at the beginning of event generation or
might be changed on an event-by-event basis.



The particles  produced in the last EVENTGEN call can be listed 
conveniently by 
      CALL print_event(1)




Related publications:
--------------------

- M"ucke A., et al. 1999, astro-ph/9903478, to appear in Comp.Phys.Commun.

- M"ucke A., et al. 1999, astro-ph/9808279, to appear in PASA.

- M"ucke A., et al. 1999, to appear in: Proc. of the
  19th Texas Symposium on Relativistic Astrophysics, Paris,
  France, Dec. 1998. Eds.: J.~Paul, T.~Montmerle \& E.~Aubourg 
  (CEA Saclay)

- M"ucke A., et al. 1999, astro-ph/9905153, to appear in: Proc. of
  19th Texas Symposium on Relativistic Astrophysics, Paris, France,
  Dec. 1998. Eds.: J.~Paul, T.~Montmerle \& E.~Aubourg 
  (CEA Saclay)

- M"ucke A., et al 1999, to appear in: Proc. of 26th Int.Cosmic
  Ray Conf. (Salt Lake City, Utah) 
