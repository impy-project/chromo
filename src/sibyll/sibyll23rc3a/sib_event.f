      SUBROUTINE SIB_START_EV (SQS, L, IA, IAFLG, NW, JDIF)
C-----------------------------------------------------------------------
C...Beginning of a SIBYLL interaction 
C.
C.  add l.m. Glauber SD cross section for pAir  13/FR
C.
C.  INPUT : SQS = c.m.s. energy (GeV)
C.          L = 1:proton, 2:charged pion
C.          IA = mass of target nucleon
C.          IAFLG = target is air
C. 
C.  OUTPUT: NW    = number of wounded nucleons
C.          JDIF(JW)  = diffraction code    !!!! changed to field !!!!
C.                  (0 : non-diffractive interaction)
C.                  (1 : forward diffraction)
C.                  (2 : backward diffraction)
C.                  (3 : double diffraction)
C.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE

c     external type declarations
      INTEGER NW_max,JDIF,IA,L,IAFLG,NW
      DOUBLE PRECISION SQS
      PARAMETER (NW_max = 20)
      DIMENSION JDIF(NW_max)
      
      DOUBLE PRECISION B, BMAX
      INTEGER NTRY, NA
      COMMON /S_CNCM0/ B, BMAX, NTRY, NA
      INCLUDE 'sib_difmass_cmmn.inc'
      DOUBLE PRECISION XI_MAX, ALAM
      COMMON /GLAUB_SCR/ XI_MAX, ALAM(61)
      include 'sib_cflafr_cmmn.inc'

c     local type declarations
      DOUBLE PRECISION SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO,
     &     SIGPROD,SIGBDIF,S_RNDM,S,PF,PB,PD,P0,P1,P2,R
      DIMENSION SIGDIF(3)
      INTEGER K

C...sample number of wounded nucleons
c     read hadron-nucleon cross section from table
      CALL SIB_SIGMA_HP(L,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO) 
      IF (IA .GT. 1)  THEN
         IF(IPAR(12).NE.0)THEN
            IF(IPAR(12).eq.3)THEN
c     distinguish between nuclear cross sections..
               IF(IAFLG.eq.0)THEN
c     if target is nucleus calc. hadron-nucleus cross section (slow)
                  CALL SIB_SIGMA_HNUC(L,IA,SQS,SIGprod,SIGbdif)
               ELSE
c     if target is air read hadron-air cross section from table
                  CALL SIB_SIGMA_HAIR(L,SQS,SIGprod,SIGbdif)
               ENDIF
            ELSE
c     always use air cross section...
               CALL SIB_SIGMA_HAIR(L,SQS,SIGprod,SIGbdif)
            ENDIF
C     2channel low-mass (coherent) diffraction?
            IF(S_RNDM(L).LT.SIGbdif/SIGprod)THEN
               NW = 1
               JDIF(1) = 1
               RETURN
            ENDIF
         ENDIF
c     sample number of wounded nucleons
         CALL INT_H_NUC (IA, SIGT, SLOPE, RHO) 
      ELSE
         NA = 1
      ENDIF      
      NW = NA

C...new treatment of diffraction 
      IF(IA.GT.1) THEN
c     hadron-nucleus case
         IF(NW.eq.1)THEN
            IF(IPAR(12).NE.0)THEN
c     high mass (incoherent) diffraction?
               S = SQS ** 2
               PF =(1.D0-dLOG(S*XI_MAX/XM2MIN(L))/
     &              dLOG(S*PAR(13)/XM2MIN(L)))*SIGDIF(1)/SIGINEL
               PB = SIGDIF(2)/SIGINEL
               PD = SIGDIF(3)/SIGINEL
            ELSE
               PF = SIGDIF(1)/SIGINEL
               PB = SIGDIF(2)/SIGINEL
               PD = SIGDIF(3)/SIGINEL
            ENDIF
         ELSE
c     Nw>1:
            IF(IPAR(12).EQ.1)THEN
c     all interactions with Nw>1 are non-diff.
               DO K=1, NW
                  JDIF(K) = 0
               ENDDO
               RETURN
            ELSE
c     some Nw>1 are attached by diff. 
               PF = PAR(124)*SIGDIF(1)/SIGINEL
               PB = PAR(124)*SIGDIF(2)/SIGINEL
               PD = PAR(124)*SIGDIF(3)/SIGINEL
            ENDIF
         ENDIF
      ELSE
c     hadron-nucleon case
         PF = SIGDIF(1)/SIGINEL
         PB = SIGDIF(2)/SIGINEL
         PD = SIGDIF(3)/SIGINEL
      ENDIF
      P0 = 1.D0-PF-PB-PD
      P1 = P0 + PF
      P2 = P1 + PB
      DO K=1, NW
         R = S_RNDM(0)
         IF (R .LT. P0)  THEN
            JDIF(K) = 0
         ELSE IF (R .LT. P1)  THEN
            JDIF(K) = 1
         ELSE IF (R .LT. P2)  THEN
            JDIF(K) = 2
         ELSE 
            JDIF(K) = 3
         ENDIF
      ENDDO
      
      END
