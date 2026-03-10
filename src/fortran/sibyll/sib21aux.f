C=======================================================================

      SUBROUTINE SIB_PARTPR(LUN)

C----------------------------------------------------------------
C     prints the particles known to SIBYLL with their internal
C     and PDG labels \FR'13
C----------------------------------------------------------------
         COMMON /S_MASS1/ AM(49), AM2(49)

         INTEGER ICHP,ISTR,IBAR
         COMMON /S_CHP/ ICHP(49), ISTR(49), IBAR(49)

         INTEGER ICHM
         COMMON /S_CHM/ ICHM(49)

         CHARACTER*6 NAMP
         COMMON /S_CNAM/ NAMP (0:49)
         SAVE

         WRITE(LUN,50)
   50    FORMAT(/,2X,16X,'SIBYLL PARTICLE TABLE:',/,2x,80('-'))
         WRITE(LUN,100)
  100    FORMAT(2X,'Particle',4X,'SIB PID',6x,'SIB2PDG',6x,'SIB2PDG^-1',
     &      4x,'MASS',4x,'STRG',4x,'CHRM',4x,'BRYN'/, 2X,80('-'))

         DO J=1,49
            IA = ISIB_PID2PDG( j )
            IF(IA.ne.0)THEN
               ISIBPDG2PIDIA=ISIB_PDG2PID( IA )
            ELSE
               WRITE(LUN,'(1X,A,I2)') 'PDG conversion not found! pid=',j
            ENDIF
            WRITE (LUN,120)  NAMP(J), J, IA, ISIBPDG2PIDIA, AM(J),
     &        ISTR(J), ICHM(J), IBAR(J)
         ENDDO
  120    FORMAT(4X,A6,4X,I4,7X,I7,8X,I4,5X,F9.3,3(6X,I2))

      END

C=======================================================================

      INTEGER FUNCTION ISIB_PID2PDG(Npid)

C----------------------------------------------------------------
C     conversion of SIBYLL internal particle code to PDG standard
C
C     input:     Npid        internal particle number
C     output:    sib_pid2pdg  PDG particle number
C
C     based on similar phojet function \FR'13
C----------------------------------------------------------------
         COMMON /S_PDG2PID/ ID_PDG_LIST(99),ID_LIST(577)
         INTEGER NPIDA,NPID
         SAVE

         Npida = iabs(Npid)
         ISIB_PID2PDG = ID_PDG_LIST(Npida)
         IF(NPID.lt.0)ISIB_PID2PDG = isign(ISIB_PID2PDG,Npid)
         RETURN
      END

C=======================================================================

      INTEGER FUNCTION ISIB_PDG2PID(Npdg)

C-----------------------------------------------------------------------
C     conversion of PDG standard particle code to SIBYLL internal
C
C     input:     Npdg        PDG particle number
C     output:    sib_pdg2pid internal particle id
C
C     based on similar phojet function \FR'13
C----------------------------------------------------------------
         COMMON /S_PDG2PID/ IPID_PDG_LIST(99),ID_LIST(577)

         INTEGER NCALL, NDEBUG, LUN
         COMMON /S_DEBUG/ NCALL, NDEBUG, LUN

         COMMON /S_CSYDEC/ CBR(102), KDEC(612), LBARP(49), IDB(49)
         SAVE

         Nin = abs(Npdg)
         if((Nin.gt.999999).or.(Nin.eq.0)) then
C  invalid particle number
            if(ndebug.gt.5) write(6,'(1x,A,I10)')
     &      ' ISIB_PDG2PID: invalid PDG ID number ',Npdg
            ISIB_PDG2PID = 0
            return
         else If(Nin.le.577) then
C  simple case
            Nout = Nin
         else
C  use hash algorithm
            Nout = mod(Nin,577)
         endif

  100    continue

C  particle not in table
         if(ID_list(Nout).Eq.0) then
            if(ndebug.gt.0) write(6,'(1x,A,I10)')
     &      ' ISIB_PDG2PID: particle not in table ',Npdg
            ISIB_PDG2PID = 0
            return
         endif
         ID_out = ID_list(Nout)
! #ifdef SIBYLLSP
         ! IF(abs(ID_out).gt.49)then
! #else
         IF(abs(ID_out).gt.99)then
! #endif
            ISIB_PDG2PID = 0
            return
         else

            if(IPID_PDG_LIST(ID_list(Nout)).eq.Nin) then
C     particle ID found
               ISIB_PDG2PID = ID_list(Nout)
               if (NPDG.lt.0) ISIB_PDG2PID = lbarp( ISIB_PDG2PID )
               return
            else
C     increment and try again
               Nout = Nout + 5
               If(Nout.gt.577) Nout = Mod(Nout,577)
               goto 100
            endif
         endif
      END

C=======================================================================

      SUBROUTINE PDG_INI

C-----------------------------------------------------------------------
C     PDG conversion blocks \FR'13
C----------------------------------------------------------------
         INTEGER NCALL, NDEBUG, LUN
         COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
         PARAMETER ( ID_PDG_MAX = 260 )
         COMMON /S_PDG2PID/ ID_PDG_LIST(99),ID_LIST(577)
         SAVE
         DATA ID_PDG_LIST /22,-11,11,-13,13,111,211,-211,321,-321, !10
     &     130,310,2212,2112,12,-12,14,-14,-2212,-2112,         !20
     &     311,-311,221,331,213,-213,113,323,-323,313,          !30
     &     -313,223,333,3222,3212,3112,3322,3312,3122,2224,     !40
     &     2214,2114,1114,3224,3214,3114,3324,3314,3334,0,      !50
     &     202212,202112,212212,212112,4*0,411,-411,            !60
     &     900111,900211,-900211,7*0,                           !70
     &     421,-421,441,431,-431,433,-433,413,-413,423,         !80
     &     -423,0,443,4222,4212,4112,4232,4132,4122,-15,        !9
     &     15,-16,16,4224,4214,4114,4324,4314,4332/

         IF(Ndebug.gt.2)
     &   WRITE(LUN,*) ' INITIALIZING PDG TABLES..'
         CALL SIB_CPCINI(ID_pdg_max,ID_pdg_list,ID_list)

      END

C=======================================================================

      SUBROUTINE SIB_CPCINI(Nrows,Number,List)

C-----------------------------------------------------------------------
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

         INTEGER NCALL, NDEBUG, LUN
         COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
         integer Number(*),List(*),Nrows
         Integer Nin,Nout,Ip,I
         SAVE

         do I = 1,577
            List(I) = 0
         enddo

C    Loop over all of the elements in the Number vector

         Do 500 Ip = 1,Nrows
            Nin = Number(Ip)

C    Calculate a list number for this particle id number
            If(Nin.Gt.999999.or.Nin.Le.0) Then
               Nout = -1
            Else If(Nin.Le.577) Then
               Nout = Nin
            Else
               Nout = Mod(Nin,577)
            End If

  200       continue

            If(Nout.Lt.0) Then
C    Count the bad entries
               IF(Ndebug.gt.3) Write(LUN,'(1x,a,i10)')
     &            ' SIB_CPCINI: invalid particle ID',Nin
               Go to 500
            End If
            If(List(Nout).eq.0) Then
               List(Nout) = Ip
            Else
               If(Nin.eq.Number(List(Nout))) Then
                  IF(Ndebug.gt.3)Write(LUN,'(1x,a,i10)')
     &              ' SIB_CPCINI: double particle  ID',Nin
               End If
               Nout = Nout + 5
               If(Nout.Gt.577) Nout = Mod(Nout, 577)

               Go to 200
            End If
  500    Continue

      END
