      program main
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CLDIF/ LDIFF
 
      ndebug = 0
      lun = 6
      ldiff = 0

      call rnd_ini
      call sibyll_ini

c     set random number generator to specific state
c     determined by previous runs of test_sibyll2..
c      call pho_RNDST(1,FILENA)
      ii = 0
 100  FORMAT(/,35('#'),I4,2X,35('#'),/)

      ecm = 7000.
      NEVT = 15
      ibeam = 13
      itarget = 1
      do i = 1, NEVT
         ii = i
         write(6,100) , ii
         print *, 'rndm:' , s_rndm(0)
         call sibyll(ibeam,itarget,ecm)
         call prnt_prtn_stck
c         call decsib
         call sib_list(6)
c         do ii=1,NP
c            IF(niorig(ii).gt.nwd) stop
c         enddo
      enddo
      write(6,*) 'Ncall, MC efficiency:' , Ncall , NEVT/real(NCALL)
     &     , 1/sqrt(real(NEVT))
      end


