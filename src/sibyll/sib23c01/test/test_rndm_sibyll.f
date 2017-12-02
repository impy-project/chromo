      PROGRAM MAIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER LDIFF
      COMMON /S_CLDIF/ LDIFF
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      CHARACTER *(20) FILENA
      DATA FILENA /'test_rndm_sibyll.rnd'/

c      ndebug = 6
c      ndebug = 1
      lun = 7
c      ldiff = 1

      imode = 9
      
      call sibyll_ini
      call rnd_ini
      
      call test_parameters(0)

 100  FORMAT(/,35('#'),2X,I8,2x,'/',I8,2x,'/',F8.3,2X,'/',2X,F8.5,
     &2X,'/',F8.5,2X,'/',F8.5,2X,'/',F8.5,2X,35('#'),/)
      
      nevt = 0
      nprt = 10000
      rtime = 0.D0
      rtime_max = 0.D0
      rtime_tot_max = 0.D0
      rtime_tot = 0.D0

      call CPU_TIME(start)

      DO WHILE (nevt>-1)
c$$$         IF(0.5.lt.S_RNDM(nevt))then
c$$$            ibeam = 6 + int(8.99*s_rndm(nevt))
c$$$            IS = -1 + 2*INT((2.D0-EPS8)*S_RNDM(nevt))
c$$$            if(ibeam.gt.12) ibeam = is*ibeam
c$$$         else
c$$$            ibeam = 34 + int(5.99*s_rndm(nevt))
c$$$            IS = -1 + 2*INT((2.D0-EPS8)*S_RNDM(nevt))
c$$$            ibeam = isign(ibeam,is)
c$$$         endif
c$$$c         ibeam = -13
         xaenergy = 1. + 5 * s_rndm(nevt)
         xenergy = 10.**xaenergy
c         xenergy = 139.37054329156945
c$$$         itarget = int(1.99D0*s_rndm(nevt))
c     itarget = 0
         call sel_config(imode,ibeam,itarget)
         call CPU_TIME(start)
         if(nevt.lt.10)
     &        write(6,*)'(beam,target,energy):',ibeam,itarget,xenergy
         ncall0 = ncall
         CALL sibyll(ibeam,itarget,xenergy)
         call CPU_TIME(finish)
         nevt = nevt + 1
         rtime = rtime + (finish-start)
         rtime_tot = rtime_tot + (finish-start)
         rtime_max = max((finish-start),rtime_max)
         if((finish-start).gt.rtime_tot_max)then
            rtime_tot_max = (finish-start)
            write(6,*)'(beam,target,energy,rtime,ncall,nrej):',
     &           ibeam,itarget,xenergy,rtime_tot_max,ncall,ncall-ncall0
            write(lun,*)'(beam,target,energy,rtime,nevt,ncall,nrej):',
     &           ibeam,itarget,xenergy,rtime_tot_max,
     &           nevt,ncall,ncall-ncall0
            call prnt_prtn_stck
            call sib_list(lun)
            call pho_RNDST(2,FILENA)
         endif
         start = finish
         if(mod(nevt,nprt).eq.0) then
            write(6,100) , nevt , ncall, nevt/real(ncall),
     &           rtime/dble(nprt),rtime_max,
     &           rtime_tot/real(nevt),rtime_tot_max
            rtime = 0.D0
            rtime_max = 0.D0
         endif
c     call sib_list(6)
         call test_parameters(1)
      ENDDO
      write(6,*)'ncall,efficiency:' , ncall,real(nevt)/real(ncall)
     &     , 1/sqrt(real(NEVT))
      END


      subroutine test_parameters(imod)
      implicit none
      save
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER LDIFF
      COMMON /S_CLDIF/ LDIFF
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10
      
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      integer i,ipar0,imod
      double precision par0
      dimension ipar0(nipar_max), par0(npar_max)
      select case(imod)
      case(0)
         do i=1,NPAR_max
            par0(i)=par(i)
         enddo
         do i=1,NIPAR_max
            ipar0(i)=ipar(i)
         enddo
         return
      case(1)
      do i=1,NPAR_max
         if(PAR(I)-PAR0(i).gt.eps3)then
            write(lun,*) 'inconsistent parameter! (par,ncall)',i,ncall
            stop
         endif
      ENDDO
      do i=1,NIPAR_max
         if(IPAR(I)-IPAR0(i).ne.0)then
            write(lun,*) 'inconsistent parameter! (ipar,ncall)',i,ncall
            stop
         endif
      ENDDO
      return
      case default
         return
      end select
      end

      SUBroutine sel_config(IMOD,IBM,ITG)
      IMPLICIT NONE
      INTEGER IMOD,IBM,ITG
      INTEGER IS
      DOUBLE PRECISION TWO,S_RNDM,EPS8
      DATA TWO,EPS8 /2.D0,1.D-8/
      INTEGER iCHMBM
      DIMENSION iCHMBM(10)
      DATA iCHMBM /59,60,71,72,74,75,87,88,89,99/
c      WRITE(6,*) ' SEL_CONFIG: input (IMOD,IBM,ITG)',IMOD,IBM,ITG
      
      SELECT CASE(IMOD)
      CASE(1)
         IF(0.5.lt.S_RNDM(IMOD))then
            ibm = 6 + int(8.99*s_rndm(IMOD))
            IS = -1 + 2*INT((TWO-EPS8)*S_RNDM(IMOD))
            if(ibm.gt.12) ibm = is*ibm
         else
            ibm = 34 + int(5.99*s_rndm(IMOD))
            IS = -1 + 2*INT((TWO-EPS8)*S_RNDM(IMOD))
            ibm = isign(ibm,is)
         endif
         itg = int((TWO-EPS8)*s_rndm(IMOD))
         
      CASE(2)
c     proton - proton interactions                           
         ibm = 13
         itg = 1

      CASE(3)
c     proton - air interactions                           
         ibm = 13
         itg = 0

      CASE(4)
c     pion - proton interactions                           
         ibm = 7
         itg = 1

      CASE(5)
c     pion - air interactions                  
         ibm = 7
         itg = 0
         
      CASE(6)
c     hyperon - proton interactions         
         ibm = 34 + int(5.99*s_rndm(IMOD))
         IS = -1 + 2*INT((TWO-EPS8)*S_RNDM(IMOD))
         ibm = isign(ibm,is)
         itg = 1

      CASE(7)
c     hyperon - air interactions         
         ibm = 34 + int(5.99*s_rndm(IMOD))
         IS = -1 + 2*INT((TWO-EPS8)*S_RNDM(IMOD))
         ibm = isign(ibm,is)
         itg = 0

      CASE(8)
c     charmed - proton interactions
         ibm = ichmbm( 1+int(9.99*s_rndm(IMOD) ) )
         IS = 1
         IF(ibm.gt.83)
     &        IS = -1 + 2*INT((TWO-EPS8)*S_RNDM(IMOD))
         ibm = isign(ibm,is)
         itg = 1

      CASE(9)
c     charmed - air interactions
         ibm = ichmbm( 1+int(9.99*s_rndm(IMOD) ) )
         IS = 1
         IF(ibm.gt.83)
     &        IS = -1 + 2*INT((TWO-EPS8)*S_RNDM(IMOD))
         ibm = isign(ibm,is)
         itg = 0
         
      CASE default
         print *,'wrong setup!'
c         stop
         
      END SELECT
      END
      
