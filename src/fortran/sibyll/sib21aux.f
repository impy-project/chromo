      SUBROUTINE PDG_INI
C----------------------------------------------------------------
C     PDG conversion blocks \FR'13
C----------------------------------------------------------------
      SAVE
      COMMON /S_DEBUG/ Ncall, Ndebug, lun
      PARAMETER ( ID_PDG_MAX = 99 )
      COMMON /S_PDG2PID/ ID_PDG_LIST(ID_PDG_MAX),ID_LIST(577)
      DATA ID_PDG_LIST /22,-11,11,-13,13,111,211,-211,321,-321, !10
     &     130,310,2212,2112,12,-12,14,-14,-2212,-2112,         !20
     &     311,-311,221,331,213,-213,113,10321,-10321,10311,    !30
     &     -10311,223,333,3222,3212,3112,3322,3312,3122,2224,   !40
     &     2214,2114,1114,3224,3214,3114,3324,3314,3334,0,      !50
     &     8*0,411,-411,10*0,                                   !70
     &     421,-421,441,431,-431,433,-433,413,-413,423,         !80
     &     -423,0,443,4222,4212,4112,4232,4132,4122,0,          !90
     &     0,0,0,4224,4214,4114,4324,4314,4332/

      IF(Ndebug.gt.2)
     & WRITE(lun,*) 'INITIALIZING PDG TABLES..'
      CALL pho_cpcini(ID_pdg_max,ID_pdg_list,ID_list)
      
      END

      INTEGER FUNCTION ISIB_PDG2PID(Npdg)
C----------------------------------------------------------------
C     conversion of PDG standard particle code to SIBYLL internal
C
C     input:     Npdg        PDG particle number
C     output:    sib_pdg2pid internal particle id
C
C     based on similar phojet function \FR'13
C----------------------------------------------------------------
      SAVE
      COMMON /S_PDG2PID/ IPID_PDG_LIST(99),ID_LIST(577)
      COMMON /S_DEBUG/ Ncall, Ndebug, lun
      DIMENSION LBAR(99)
      DATA LBAR /1,3,2,5,4,6,8,7,10,9,11,12,-13,-14,16,15,18,17,13,14,
     +  22,21,23,24,26,25,27,29,28,31,30,32,33,-34,-35,-36,-37,-38,-39,
     +  -40,-41,-42,-43,-44,-45,-46,-47,-48,-49,9*0,60,59,10*0,72,71,
     +  73,75,74,77,76,79,78,81,80,0,83,-84,-85,-86,-87,-88,-89,4*0,-94,
     +  -95,-96,-97,-98,-99 /

      Nin = abs(Npdg)
      if((Nin.gt.99999).or.(Nin.eq.0)) then
C  invalid particle number
        if(ndebug.gt.5) write(lun,'(1x,A,I10)')
     &    'isib_pdg2pid: invalid PDG ID number ',Npdg
        isib_pdg2pid = 0
        return
      else If(Nin.le.577) then
C  simple case
        Nout = Nin
      else
C  use hash algorithm
        Nout = mod(Nin,577)
      endif
 100  continue

C  particle not in table
      if(ID_list(Nout).Eq.0) then
        if(ndebug.ge.0) write(lun,'(1x,A,I10)')
     &    'isib_pdg2pid: particle not in table ',Npdg
        isib_pdg2pid = 0
        return
      endif

      if(IPID_pdg_list(ID_list(Nout)).eq.Nin) then
C  particle ID found
        isib_pdg2pid = ID_list(Nout)
        if (NPDG.lt.0) isib_pdg2pid = lbar( isib_pdg2pid )
        return
      else
C  increment and try again
        Nout = Nout + 5
        If(Nout.gt.577) Nout = Mod(Nout,577)
        goto 100
      endif

      END
      
      SUBROUTINE SIB_PARTPR(LUN)
C----------------------------------------------------------------
C     prints the particles known to SIBYLL with their internal
C     and PDG labels \FR'13
C----------------------------------------------------------------
      COMMON /S_CNAM/ NAMP (0:49)
      CHARACTER NAMP*6

      WRITE(LUN,50)
 50   FORMAT(2X,'SIBYLL PARTICLE TABLE:',/,2x,48('-'))
      WRITE(LUN,100)
 100  FORMAT(2X,'Particle',4X,'SIB PID',6x,'SIB2PDG',6x,'SIB2PDG^-1', 
     &     /, 2X,48('-'))

      DO J=1,99
         IA = ISIB_PID2PDG( j )         
         IF(IA.ne.0)
     &        WRITE (LUN,120) NAMP(J), J, IA, ISIB_PDG2PID( IA )
      ENDDO
 120  FORMAT(4X,A6,4X,I4,8X,I6,8X,I4)

      END


      INTEGER FUNCTION ISIB_PID2PDG(Npid)
C----------------------------------------------------------------
C     conversion of SIBYLL internal particle code to PDG standard
C
C     input:     Npid        internal particle number
C     output:    sib_pid2pdg  PDG particle number
C
C     based on similar phojet function \FR'13
C----------------------------------------------------------------
      SAVE
      COMMON /S_PDG2PID/ ID_PDG_LIST(99),ID_LIST(577)
      INTEGER NPIDA,NPID

      Npida = iabs(Npid)
      isib_pid2pdg = ID_PDG_LIST(Npida)
      IF(NPID.lt.0)isib_pid2pdg = isign(isib_pid2pdg,Npid)
      RETURN
      END



      SUBROUTINE pho_cpcini(Nrows,Number,List)
C***********************************************************************
C
C     initialization of particle hash table
C
C     input:   Number     vector with Nrows entries according to PDG
C                         convention
C
C     output:  List       vector with hash table
C
C     (this code is based on the function initpns written by
C      Gerry Lynch, LBL, January 1990)
C
C***********************************************************************
      IMPLICIT NONE
      SAVE
      integer Number(*),List(*),Nrows
      Integer Nin,Nout,Ip,I

      do I = 1,577
        List(I) = 0
      enddo
C    Loop over all of the elements in the Number vector

        Do 500 Ip = 1,Nrows
            Nin = Number(Ip)

C    Calculate a list number for this particle id number
            If(Nin.Gt.99999.or.Nin.Le.0) Then
                 Nout = -1
            Else If(Nin.Le.577) Then
                 Nout = Nin
            Else
                 Nout = Mod(Nin,577)
            End If

 200        continue

            If(Nout.Lt.0) Then
C    Count the bad entries
c                Write(6,'(1x,a,i10)')
c     &            'pho_cpcini: invalid particle ID',Nin
                Go to 500
            End If
            If(List(Nout).eq.0) Then
                List(Nout) = Ip
            Else
                If(Nin.eq.Number(List(Nout))) Then
c                  Write(6,'(1x,a,i10)')
c     &              'pho_cpcini: double particle ID',Nin
                End If
                Nout = Nout + 5
                If(Nout.Gt.577) Nout = Mod(Nout, 577)

                Go to 200
            End If
 500      Continue

      END

      SUBROUTINE SIBHEP1
C-----------------------------------------------------------------------
C  Convert to HEPEVT common block
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL P
      INTEGER NP,LLIST,LLIST1,NP_max,NEVSIB
      DATA NEVSIB /0/
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      COMMON /S_PLIST1/ LLIST1(8000)
      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(49), ISTR(49), IBAR(49)

      INTEGER NEVHEP,NMXHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      DOUBLE PRECISION PHEP,VHEP
      PARAMETER (NMXHEP=NP_max)
      COMMON /HEPEVT/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &                JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),
     &                VHEP(4,NMXHEP)

      INTEGER ICHG
      COMMON /SCHG/  ICHG(NMXHEP)

      INTEGER I, ISIB_PID2PDG
      EXTERNAL ISIB_PID2PDG

      SAVE NEVSIB

      NHEP = NP
      NEVHEP = NEVSIB

      DO I=1,NHEP
         IF (ABS(LLIST(I)).LT.10000) THEN
            ISTHEP(I) = 1
         ELSE
            ISTHEP(I) = 2
         END IF
         ICHG(I) = ICHP(ABS(LLIST(I)))
         IDHEP(I) = ISIB_PID2PDG(MOD(LLIST(I),10000))
         JMOHEP(1,I) = LLIST1(I)
         JMOHEP(2,I) = LLIST1(I)
         PHEP(1,I) = DBLE(P(I,1))
         PHEP(2,I) = DBLE(P(I,2))
         PHEP(3,I) = DBLE(P(I,3))
         PHEP(4,I) = DBLE(P(I,4))
         PHEP(5,I) = DBLE(P(I,5))
      END DO

      NEVSIB = NEVSIB + 1
      END
