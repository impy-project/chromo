C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

      SUBROUTINE FASTJETPPgenkt(P,NPART,R,PALG,F77JETS,NJETS
     .,ARRAY1,ARRAY2,IPRI)
      DOUBLE PRECISION P(4,*), R, PALG, F77JETS(4,*)
      DOUBLE PRECISION dP, dR, dPALG
      INTEGER          NPART, NJETS,Nd,ARRAY1(*),ARRAY2(*),IPRI
      F77JETS(1,1)=0d0
      NJETS=0
      dP=P(1,1)
      Nd=NPART
      dR=R
      dPALG=PALG
      array1(1)=0
      array2(1)=0
      idum=ipri
c      write(*,*)"FastJet called with :",F77JETS(1,1),NJETS,R,PALG
c      stop"But FastJet not installed !"
      END
