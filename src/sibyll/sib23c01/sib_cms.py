#!/usr/bin/env python

import sib23ccrs as sib
import numpy as np
import matplotlib.pyplot as plt

# initialize random number generator
sib.rnd_ini()

# initialize sibyll
sib.sibyll_ini()

# call separate decay ini
# declares stable/unstable
sib.dec_ini()

# number of events
Nevt = 10000

# 7TeV center-of mass
ecm = 7000.

# proton beam
ibeam = 13

# proton target (one nucleon)
itarget = 1

# oxygen target (16 nucleons)
#itarget = 16

nRapBins = 201
rap_hist = np.zeros(nRapBins-1)

for evt in xrange(Nevt):

    # create event
    sib.sibyll(ibeam,itarget,ecm)

    # call decay routine
    sib.decsib()

    # read particle stack

    # number of particles
    numberOfParticles = sib.s_plist.np
    
    # list of particles in event
    idList = sib.s_plist.llist[ : numberOfParticles ]

    # longitudinal momentum
    pz = sib.s_plist.p[ : numberOfParticles , 2 ]
    pts = sib.s_plist.p[ : numberOfParticles , 0 ] ** 2 + sib.s_plist.p[ : numberOfParticles , 1 ] ** 2
    mts = pts + sib.s_plist.p[ : numberOfParticles, 4 ] ** 2
    # do something...
    # calc rapidity
    id_pos = sib.s_plist.p[ :numberOfParticles , 2] > 0
    id_neg = sib.s_plist.p[ :numberOfParticles , 2] < 0
    rap = np.zeros( numberOfParticles ,dtype=pz.dtype)
    try:
        rap_pos = \
            np.log( ( np.sqrt( pz[ id_pos ]**2 + mts[ id_pos ] ) + pz[ id_pos ] ) / np.sqrt(mts)[ id_pos ] )
    except FloatingPointError:
        print 'calculation of rapidity failed!! (pz,mts):',pz[ id_pos ],mts[ id_pos ]

    try:
        rap_neg = \
            np.log(  np.sqrt(mts)[ id_neg ] \
                         / (np.sqrt(pz[ id_neg ]**2 \
                                        + mts[ id_neg ]) - pz[ id_neg ]))
    except FloatingPointError:
        print 'calculation of rapidity failed!! (pz,mts):',pz[ id_pos ],mts[ id_pos ]

    rap[ id_pos ] = rap_pos
    rap[ id_neg ] = rap_neg

    rap_hist += np.histogram( rap , bins = np.linspace(-10,10,nRapBins) )[0]
    # e.g call sib_list, which prints event
#    sib.sib_list(6)



plt.plot(np.linspace(-10,10,nRapBins)[:-1], rap_hist , ls='steps-post')
plt.xlabel('Pseudorapidity')
plt.ylabel('dN/deta')
plt.show()
