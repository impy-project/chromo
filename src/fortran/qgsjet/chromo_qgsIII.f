C This file are aux functions for stand-alone use
C of CORSIKA interaction models

      SUBROUTINE CQGSINI( DATDIR, LUN, DEBUGNUM )
C-----------------------------------------------------------------------
C  C(ompact) Q(UARK) G(LUON) S(TRING JET MODEL) INI(TIALZATION)
C
C  INITIALIZES QGSJET-II MODEL quitely and fully.
C  THIS SUBROUTINE IS CALLED FROM START.
C-----------------------------------------------------------------------
         IMPLICIT NONE
         CHARACTER*256 DATDIR
         COMMON /AREA40/  JDIFR
         INTEGER          JDIFR
         INTEGER        VERSION
         INTEGER    LUN, moniou
         INTEGER    DEBUGNUM, debug
         common /qgarr43/ moniou
         common /qgdebug/ debug
         SAVE
C-----------------------------------------------------------------------

C  COMMON MODEL PARAMETERS SETTING
         CALL QGSET
         moniou = LUN
         debug = DEBUGNUM
C  Call global initialization subroutine
         CALL QGAINI( DATDIR )
C  Call cross-section initialization subroutine
C      CALL QGSSIGINI

      END

      SUBROUTINE CHEPEVT
C-----------------------------------------------------------------------
C  Convert to HEPEVT common block
C
C
C  FOLLOWING NOTATIONS FOR PARTICLES TYPES ARE USED WITHIN QGSJET:
C             0 - PI0,
C             1 - PI+,
C            -1 - PI-,
C             2 - P,
C            -2 - P-BAR,
C             3 - N,
C            -3 - N-BAR,
C             4 - K+,
C            -4 - K-,
C             5 - K0S,
C            -5 - K0L
C             6 - LAMBDA
C            -6 - LAMBDA-BAR
C             7 - DELTA++
C            -7 - DELTA++-BAR
C             8 - DELTA-
C            -8 - DELTA--BAR
C             9 - nothing
C            -9 - nothing
C            10 - ETA
C            16 - RHO+
C           -16 - RHO-
C            17 - RHO0
C            18 - K*+
C           -18 - K*-
C            19 - K*0
C           -19 - K*0-BAR
C            20 - sigma+
C           -20 - sigma+bar
C            21 - sigma-
C           -21 - sigma-bar
C-----------------------------------------------------------------------
         IMPLICIT NONE

         INTEGER NPTMAX, ICH, NSP
         DOUBLE PRECISION ESP
         PARAMETER(NPTMAX=95000)
         COMMON /QGARR12/ NSP
         COMMON /QGARR14/ ESP(4,NPTMAX),ICH(NPTMAX)


         INTEGER NEVHEP,NMXHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
         DOUBLE PRECISION PHEP,VHEP
         PARAMETER (NMXHEP=NPTMAX)
         COMMON /HEPEVT/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &      JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),
     &      VHEP(4,NMXHEP)
         INTEGER ICHG
         COMMON /QGCHG/  ICHG(NMXHEP)
C     Particle tables start with the ID -10(rho0) going through 0 (pi0).
         CHARACTER*12 NAME(-21:21)
         DATA NAME /
     &   'Sigmabar+   ', 'Sigmabar-   ', 'K*bar0      ', 'K*-         ',
     &   '            ', 'rho-        ', '            ', '            ',
     &   '            ', '            ', '            ', '            ',
     &   '            ', 'Deltabar--  ', 'Deltabar0   ', 'Lambdabar0  ',
     &   'K_L0        ', 'K-          ', 'nbar0       ', 'pbar-       ',
     &   'pi-         ', 'pi0         ', 'pi+         ', 'p+          ',
     &   'n0          ', 'K+          ', 'K_S0        ', 'Lambda0     ',
     &   'Delta++     ', 'Delta-      ', '            ', 'eta         ',
     &   '            ', '            ', '            ', '            ',
     &   '            ', 'rho+        ', 'rho0        ', 'K*+         ',
     &   'K*0         ', 'Sigma+      ', 'Sigma-      ' /

         INTEGER IPDGID(-21:21)
         DATA IPDGID /
     &   -3112, -3222, -313, -323, 0, -213, 0, 0, 0, 0,
     &   0, 0, 0, -2224, -2114, -3122, 130, -321, -2112, -2212,
     &   -211, 111, 211, 2212, 2112, 321, 310, 3122, 2224, 1114,
     &   0, 221, 0, 0, 0, 0, 0, 213, 113, 323, 313, 3222, 3112 /

         DOUBLE PRECISION QMASS(-21:21)
         DATA QMASS /
     &   1.19d0, 1.19d0, 0.496d0, 0.496d0, 0.0d0, 0.76d0,
     &   0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.19d0,
     &   1.19d0, 1.116d0, 0.496d0, 0.496d0, 0.93827999d0,
     &   0.93827999d0, 0.14d0,0.14d0, 0.14d0, 0.93827999d0,
     &   0.93827999d0, 0.496d0, 0.496d0,1.116d0, 1.23d0, 1.19d0, 0.0d0,
     &   0.548d0, 0.0d0, 0.0d0, 0.0d0,0.0d0, 0.0d0, 0.76d0, 0.76d0,
     &   0.496d0, 0.496d0, 1.19d0, 1.19d0 /

         INTEGER ICHRG(-21:21)
         DATA ICHRG /
     &   1, -1, 0, -1, 0, -1, 0, 0, 0, 0,
     &   0, 0, 0, -2, 0, 0, 0, -1, 0, -1,
     &   -1, 0, 1, 1, 0, 1, 0, 0, 2, -1,
     &   0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
     &   0, 1, -1 /

         INTEGER I

         NHEP = nsp

         DO I=1,nsp
C         WRITE(6,*) I, ich(I), esp(:,I)
            NHEP = NSP
            ISTHEP(I) = 1
            IDHEP(I) = IPDGID(ich(I))
            PHEP(1,I) = esp(3,I)
            PHEP(2,I) = esp(4,I)
            PHEP(3,I) = esp(2,I)
            PHEP(4,I) = esp(1,I)
            PHEP(5,I) = QMASS(ich(I))
            ICHG(I) = ICHRG(ich(I))
         END DO


      END

*-- Author :    T. Pierog IKP KIT KARLSRUHE       25/10/2012
C=======================================================================

      subroutine  LzmaOpenFile(name)

C-----------------------------------------------------------------------
C DUMMY FUNCTION TO BE COMPATIBLE WITH CRMC
C-----------------------------------------------------------------------

         IMPLICIT NONE
         character*256 name,name2
         name2=name
      end

*-- Author :    T. Pierog IKP KIT KARLSRUHE       25/10/2012
C=======================================================================

      subroutine  LzmaCloseFile()

C-----------------------------------------------------------------------
C DUMMY FUNCTION TO BE COMPATIBLE WITH CRMC
C-----------------------------------------------------------------------

         IMPLICIT NONE
      end

*-- Author :    T. Pierog IKP KIT KARLSRUHE       25/10/2012
C=======================================================================

      subroutine LzmaFillArray(dum,idum)

C-----------------------------------------------------------------------
C DUMMY FUNCTION TO BE COMPATIBLE WITH CRMC
C-----------------------------------------------------------------------

         IMPLICIT NONE
         double precision dum,dum2
         integer idum,idum2
         dum2=dum
         idum2=idum
      end

*-- Author :    T. Pierog IKP KIT KARLSRUHE       25/10/2012
C=======================================================================

      integer function size(array)

C-----------------------------------------------------------------------
C DUMMY FUNCTION TO BE COMPATIBLE WITH CRMC
C-----------------------------------------------------------------------

         IMPLICIT NONE
         double precision array(*)
         size=int(array(1))
      end

*-- Author :    A. Fedynitch ASIoP       14/03/2024
C=======================================================================

      subroutine CQGSHH_HA_CS(gtot,gin,gel,gdp,gdt,bel)

C-----------------------------------------------------------------------
C Calculate QGSJET cross sections using native functions
C call qgini before calling this function
C-----------------------------------------------------------------------

         implicit double precision (a-h,o-z)
         common /qgarr1/  ia(2),icz,icp
         common /qgarr2/  scm,wp0,wm0
         common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
         common /qgarr15/ fp(3),rq(3),cd(2,3)
         common /qgarr49/ trnuc(56),twsnuc(56),twbnuc(56)
         dimension gz0(5),gz1(3)
Cf2py    intent(out) :: gtot,gin,gel,gdp,gdt,bel

         if (ia(1).ne.1) then
            write(*,*) 'Error: chromo_qgsect called with ia(1) ne 1'
         endif
         ! icz = hadron class (N, pi, K)
         !e0n interaction energy, scm center of mass
         rp=(rq(icz)+rq(2)+alfp*log(scm))*4.d0*.0389d0
         if(ia(2).eq.1)then                        !hadron-proton interaction
            call qgfz(0.d0,gz0,0,0)
            gtot=gz0(1)                              !total cross-section
            gin=(gz0(2)+gz0(3)+gz0(4))*.5d0          !inelastic cross section
            bel=gz0(5)                               !elastic scattering slope
            gel=gtot-gin                             !elastic cross-section
            gdp=gz0(3)*.5d0                          !proj. diffr. cross section
            gdt=gz0(4)*.5d0                          !targ. diffr. cross section
            if(debug.ge.1)write (moniou,219)gtot,gin,gel,gdp,gdt,bel

         else    !hadron-nucleus interaction
            ! nuclear densities set in qgini
            !**anfe: different from qgini!!
            ! bm=rnuc(2)+dlog(29.d0)*wsnuc(2)          !impact parameter for integr.
            ! bm=rnuc(2)+5.d0*wsnuc(2)
            anorm=qganrm(rnuc(2),wsnuc(2),wbnuc(2))*rp    !density normalization
            call qggau(gz1)                          !integration over b< bm

            !integration over b> bm
            gin=10.d0*(gz1(1)+gz1(2)+gz1(3) )        !inelastic cross section
            gdp=10.d0*gz1(2)                               !proj. diffr. cross section
            gdt=10.d0*gz1(3)                               !targ. diffr. cross section
            ! Need to multiply by 10.d0
            if(debug.ge.1)write (moniou,220) gin,gz1(3),gz1(2)
         endif
  219    format(2x,'qgaini: hadron-proton cross sections:'
     *     /4x,'gtot=',e10.3,2x,'gin=',e10.3,2x,'gel=',e10.3/4x
     *     ,'gdifrp=',e10.3,2x,'gdifrt=',e10.3,2x,'b_el=',e10.3,2x)
  220    format(2x,'qgaini: hadron-nucleus cross sections:'
     *     /4x,'gin=',e10.3,2x,'gdifr_targ=',e10.3,2x
     *     ,'gdifr_proj=',e10.3)
      end subroutine
