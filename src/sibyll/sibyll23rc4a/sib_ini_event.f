      SUBROUTINE INI_EVENT(ECM,KBEAM,IATARG,IMOD)
C-----------------------------------------------------------------------
C     initializes the stacks and event info common
c     if Imod : 0 - initiate subevent in recursive call
c                  ( keeps the final hadron stack intact )
C             : 1 - initiate entire new event
C-----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
c     external type declarations
      DOUBLE PRECISION ECM
      INTEGER KBEAM,IATARG,IMOD

c     COMMONs
      include 'sib_nw_prm.inc'
      include 'sib_int_prm.inc'
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_chist_cmmn.inc'
      INCLUDE 'sib_indx_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_rand_cmmn.inc'

c     local types
      DOUBLE PRECISION PZ,E1,PAWT,S_RNDM,R,FOX
      INTEGER KK,JJ,KBA,IREFout,JN
      DATA FOX /0.257/
      
      IF(NDEBUG.gt.0.and.IMOD.eq.1) 
     &     WRITE(LUN,'(A48,F8.2,I3,I3,I3)')
     &     'INI_EVENT: called with (ECM,KBEAM,IATARG,NCALL):',
     &     ECM,KBEAM,IATARG,NCALL

c     set final particle stack to zero
      IF(IMOD.eq.1)then
         NP = 0
         NWD = 0
         NJET = 0
         NSOF = 0
c     keep rand gen state at beginning of current event
c     in common sib_rand
!          CALL PHO_RNDSO(U2,C2,CD2,CM2,II2,JJ2)
      endif

      call ini_prtn_stck(0,0)

c     clear index cache
      do kk=1,3
         IBMRDX(kk) = 0
      ENDDO
      do jj=1,NW_max
         do kk=1,3
            ICSTDX(jj,kk) = 0
            ICSTDX(jj+1,kk) = 0
            ITGRDX(jj,kk) = 0
            IINTDX(jj) = 0
         ENDDO
      ENDDO

      SQS   = Ecm
      S     = SQS*SQS
      
      KB = KBEAM
      KBA = IABS(KBEAM)
c     add beam particles to parton stack, lvl -2
      PZ = PAWT(SQS,AM(KBA),AM(13))
      E1 = SQRT(PZ**2+AM2(KBA))
      call add_prtn(ZERO,ZERO,PZ,E1,AM(KBA),KB,-2,0,IREFout)
      IF(IMOD.eq.1)THEN
         IAT = IATARG
         IF(IATARG.EQ.1)THEN
            KT(1) = 13
         ELSE
            IF(IATARG.eq.0)THEN
C...  Generate an 'air' interaction by choosing Nitrogen or Oxygen
               R = S_RNDM(0)
               IATARG = 14
               IF (R .LT. FOX)  IATARG = 16
            ENDIF
            DO JN=1,IATARG
c     for nuclear target: proton (13) or neutron (14)
               KT(JN) = 13 + INT((TWO-EPS8)*S_RNDM(JN))
            ENDDO
         ENDIF
      ELSE
         KT(1) = IATARG
      ENDIF

C...energy-dependent transverse momentum cutoff
c...EJA correction 2007.03.27
      IF(IPAR(27).eq.1)THEN
         PTmin = PAR(10)+PAR(11)*EXP(PAR(12)*SQRT(LOG(SQS)))
      else
         PTmin = PAR(10)+PAR(11)*EXP(PAR(12)*SQRT(LOG(S)))
      endif
      XMIN = FOUR*PTmin**2/S
      ZMIN = LOG(XMIN)
      IF(ndebug.gt.0)then
         write(lun,*) 'INI_EVNT: ncall:', ncall
         write(lun,'(A33,F8.2,F8.2,F8.5,E8.3,F8.5)')
     &        'INI_EVNT: (SQS,S,PTmin,Xmin,Zmin)',
     &        SQS,S,PTmin,Xmin,Zmin
         write(lun,*) 'INI_EVNT: KB,IAT,KT',KB,IAT,(KT(jj),jj=1,IATARG)
      endif
      CALL PTSETUP_4FLV(ECM)

      END
