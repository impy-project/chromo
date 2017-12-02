      program main
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      CHARACTER *(18) FILENA
      DATA FILENA /'sib-dminus-bug.rnd'/

c      ndebug = 8
      lun = 7

      call rnd_ini
      CALL SIBYLL_INI
      call dec_ini

c     set random number generator to specific state
c     determined by previous runs of test_sibyll2..
c      call pho_RNDST(1,FILENA)

c      lun = 6
      
      ii = 0
 100  FORMAT(/,35('#'),2X,I8,2x,'/',I8,2x,'/',F8.3,2X,35('#'),/)


      rtime = 0.D0
      rtime_max = 0.D0
      rtime_tot_max = 0.D0
      rtime_tot = 0.D0

      call CPU_TIME(start)
      rtime_last = 1.d0
      
      ecm = 30.
      NEVT = 100000
      ibeam = 13
      itarget = 1
      
      do i = 1, NEVT
         ii = i
         call CPU_TIME(start)
         if(mod(ii,NEVT/10).eq.0) then
             write(6,100) , ii , ncall, ii/real(ncall)
             write(6,*) 'rtime intv.,rtime/rtime_last,rtime max,',
     &            'rtime tot,rtime max tot',
     &            rtime/dble(NEVT/10),rtime/rtime_last,rtime_max,
     &            rtime_tot/real(i),rtime_tot_max
           rtime_last = rtime
           rtime = 0.D0
           rtime_max = 0.D0           
         endif
         call sibyll(ibeam,itarget,ecm)
         call decsib
         call CPU_TIME(finish)
         rtime = rtime + (finish-start)
         rtime_tot = rtime_tot + (finish-start)
         rtime_max = max((finish-start),rtime_max)
         if((finish-start).gt.rtime_tot_max)then
            rtime_tot_max = (finish-start)
         endif
         start = finish
c         call sib_list(6)
      enddo

      write(6,*) 'Ncall, MC efficiency:' , Ncall , NEVT/real(NCALL)
     &     , 1/sqrt(real(NEVT))
      write(6,*) ,'total. runtime:', rtime_tot
      end

