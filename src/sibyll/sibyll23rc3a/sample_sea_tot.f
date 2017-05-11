
      SUBROUTINE SAMPLE_SEA_TOT
     &     (KRMNT,NINT,NSEA,XGAM,XJET,STR_MASS,XSJ,XX)
      IMPLICIT NONE
      SAVE
Cf2py double precision, intent(in) :: xgam,xjet,str_mass
Cf2py integer, intent(in) :: Nsea,nint,krmnt
Cf2py double precision, intent(out) :: xsj,xx

c     include COMMON blocks
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      
c     input/output type definitions
      DOUBLE PRECISION XGAM,XJET,STR_MASS,XSEA,XX,XSJ
      DIMENSION XX(2*NW_max+2)

      INTEGER NSEA,NINT,KRMNT

c     local type definitions
      DOUBLE PRECISION AC,GAMMA,S_RNDM,XA,XREM,R,Z,Z1,Z2,XMINA
      DATA AC /-0.2761856692D0/ ! log(2) - gamma(Eulero)
      INTEGER j,jj,ilast

      GAMMA = xgam
      XMINA = TWO*STR_mass/SQS
      IF(IPAR(73).eq.1.and.NINT.gt.1) GAMMA = PAR(119)
      IF(ndebug.gt.3) THEN
         WRITE(LUN,*)'SAMPLE_SEA_TOT: called with ',
     &        '(KRMNT,NINT,NSEA,XGAM,XJET,STR_MASS):', 
     &        KRMNT,NINT,NSEA,XGAM,XJET,STR_MASS
         
         WRITE(LUN,*)'SAMPLE_SEA_TOT: XMIN,XMIN*N,XREM:',
     &        XMINA,NSEA*XMINA,ONE-XJET
      ENDIF
c     sample total fraction for sea partons..
      Z1 = LOG(DBLE(NSEA))
 50   Z2 = LOG(HALF*SQS*(ONE-XJET)/STR_MASS-TWO)
      R = S_RNDM(0)
      Z=(Z1+AC)*(ONE+R*(((Z2+AC)/(Z1+AC))**NSEA-ONE))
     &     **(ONE/DBLE(NSEA))-AC
      XSEA = XMINA*EXP(Z)
      IF(ndebug.gt.3) WRITE(LUN,*) ' total SEA fraction:' , xsea
      IF ( (ONE-XSEA)**GAMMA .LT. S_RNDM(0)) GOTO 50
c     maximal fraction remaining for valence..
 60   XREM = XSEA - DBLE(Nsea)*XMINA
      IF(ndebug.gt.3) 
     &     WRITE(LUN,*) ' Xsea,xval,xjet:',
     &     xsea,ONE-XSEA-XJET,xjet
      
C...  Split the energy of sea partons among the different partons
      DO j=1,Nsea-1
         jj = 2+j
         IF(KRMNT.eq.0) jj = 4+j
c     fraction for first parton
         XA = XREM*S_RNDM(J)
c     for interactions other than first decrease energy fraction
c     (beam side hadron can participate in multiple binary collisions)
c     IF(NINT.gt.1.and.j.gt.2*KRMNT) XA=SIGN(ABS(XA)**PAR(116),XA)
         XX(jj) = XMINA + XA
c     new remainder
         XREM = XREM - XA
         IF(ndebug.gt.3) write(lun,*)'x1,j,rem,xa',xX(jj),jj,xrem,xa
      enddo
c     last parton..
      ilast = 2+Nsea
      IF(KRMNT.eq.0) ilast = 4+Nsea
      XX(ILAST) = XMINA + XREM

c     break symmetry between nucleon interactions
c     first interaction takes most energy
      IF(NINT.gt.1.and.IPAR(71).eq.1)THEN
         JJ = 3
         IF(KRMNT.eq.0) JJ = 5
         if(ndebug.gt.4) write(lun,*) 'x1+x2,p*xeq:',
     &        XX(JJ)+XX(JJ+1),PAR(117)*XSEA/NINT
         IF(XX(JJ)+XX(JJ+1).lt.PAR(117)*XSEA/NINT) GOTO 60
      ENDIF

      XSJ = XSJ + XSEA
      IF(ndebug.gt.3)THEN  
         write(lun,*)'x1,N,rem',xx(ilast),ilast,xrem
         write(lun,*) 'xseajet',xsj
      endif

      END
