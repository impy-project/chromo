      SUBROUTINE INT_H_NUC (IA, SIGT, SLOPE, RHO) 
C...Compute with a montecarlo method the "multiple interaction structure"
C.  of an hadron-nucleus collision.
C.  
C.
C.  INPUT : IA               = mass of target nucleus
C.          SIGT (mbarn)     = total hp cross section
C.          SLOPE (GeV**-2)  = slope of hp elastic scattering
C.          RHO              = real/imaginary part of forward elastic
C.                             scattering amplitude
C.
C.  OUTPUT : in COMMON block /CNCMS0/
C.           B = impact parameter (fm)
C.           BMAX = maximum impact parameter for generation
C.           NTRY = number of "trials" before one interaction
C.           NA = number of wounded nucleons in A
C. Author : P.Lipari  (may 1993)
C---------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (IAMAX=56)
      COMMON /S_CNCM0/ B, BMAX, NTRY, NA
      DIMENSION XA(IAMAX), YA(IAMAX)
      INCLUDE 'sib_utl_cmmn.inc'      
      CC = SIGT/(FOUR*PI*SLOPE*CMBARN)         
      DEN = TWO*SLOPE*CMBARN*0.1D0
      BMAX = ONE*10             ! fm
      NTRY = 0
      CALL NUC_CONF (IA, XA, YA)
1000  B = BMAX*dSQRT(S_RNDM(0))
      PHI = TWOPI*S_RNDM(0)
      BX = B*dCOS(PHI)
      BY = B*dSIN(PHI)
      NTRY = NTRY+1
      NA = 0
      DO JA=1,IA
         S = (XA(JA)-BX)**2 + (YA(JA)-BY)**2
         F = dEXP(-S/DEN)
         PEL = CC*CC*(ONE+RHO*RHO)*F*F
         PINEL  = TWO*CC*F-PEL
         R = S_RNDM(0)
         IF (R .LT. PINEL)  THEN
            NA = NA + 1
         ENDIF
      ENDDO
      IF (NA .EQ. 0)  GOTO 1000
      RETURN
      END

