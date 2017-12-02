      PROGRAM MAIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER ( NEVTMAX=800000 )
      DOUBLE PRECISION XRND
      DIMENSION XRND(NEVTMAX)
      COMMON /RNDMGAS/ ISET
      
c      ndebug = 7
c      lun = 6
      
c     call rnd_ini
c     21         172          56         150
      ISET = 0
      CALL PHO_RNDIN(21,172,56,150)
C       NA1,NA2,NA3,NB1  - values for initializing the generator
C                          NA? must be in 1..178 and not all 1;
C                          12,34,56  are the standard values
C                          NB1 must be in 1..168;
C                          78  is the standard value      
      call sibyll_ini
      
 100  FORMAT(/,35('#'),F8.5,2X,35('#'),/)
 101  FORMAT(/,35('#'),I6,F8.5,2X,35('#'),/)
 102  FORMAT(/,35('#'),I6,F8.5,F8.5,2X,35('#'),/)            

c     set beam and target mode:
c     random: 1,
c     proton- proton: 2
c     proton-air: 3
c     pion-proton: 4      
c     pion-air: 5
c     hyperon-proton: 6
c     hyperon-air: 7      
c     charm-proton: 8
c     charm-air   : 9
      IMODE = 1
      
      EPS = 1.e-3
      NEVT = 50000
      write(6,*)'creating ', nevt , ' events in config:', imode
      
c     first run through random number sequence
      write(6,100) S_RNDM(0)
      DO II=1,NEVT

      call sel_config(imode,ibeam,itarget)
         
c     select energy
         xaenergy = 1. + 5 * s_rndm(0)
         xenergy = 10.**xaenergy         
         
         call sibyll(ibeam,itarget,xenergy)
c         call sib_list(6)
c     store next random number
         XRND(II) = S_RNDM(II)
         if(mod(ii,NEVT/10).eq.0) then
            write(6,101) ii, XRND(II)
            write(6,*) 'current config: ibeam,itarget,xenergy',
     &           ibeam,itarget,xenergy
         endif
      ENDDO
c     print last event
c     call prnt_prtn_stck
      call sib_list(6)
      write(6,*) 'next random number: ' , s_rndm(0)
      
c     reinitialize random number generator and run again
c     call rnd_ini
      ISET = 0
      CALL PHO_RNDIN(21,172,56,150)
c      CALL PHO_RNDIN(112,36,56,78)      
      write(6,100) S_RNDM(0)
      DO II=1,NEVT
         
         call sel_config(imode,ibeam,itarget)
         
         xaenergy = 1. + 5 * s_rndm(0)
         xenergy = 10.**xaenergy         

         call sibyll(ibeam,itarget,xenergy)
c     get next random number
         XNEXT = S_RNDM(II)         
         if(mod(ii,NEVT/10).eq.0) then
            write(6,102) ii, xrnd(ii), xnext
            write(6,*) 'current config: ibeam,itarget,xenergy',
     &           ibeam,itarget,xenergy
         endif
c     check if random numbers have diverged
         IF(ABS(XRND(II)-XNEXT).gt.EPS)THEN
            write(6,*) 'sequences have diverged! position' , II,
     &           XRND(II), xnext
            call sib_list(6)
            stop
         ENDIF
      ENDDO
c     print last event
c      call prnt_prtn_stck
      call sib_list(6)
      write(6,*) 'next random number: ' , s_rndm(0)
      write(6,*) '**** success! ****'
      END      
     
      SUBroutine sel_config(IMOD,IBM,ITG)
      IMPLICIT NONE
      INTEGER IMOD,IBM,ITG
      INTEGER IS
      DOUBLE PRECISION TWO,S_RNDM,EPS8
      DATA TWO,EPS8 /2.D0,1.D-8/
      INTEGER iCHMBM
      DIMENSION iCHMBM(10)
      DATA iCHMBM /59,60,71,72,74,75,87,88,89,99/

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
         stop
         
      END SELECT
      END
