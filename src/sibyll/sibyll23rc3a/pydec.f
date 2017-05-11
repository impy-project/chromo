      SUBROUTINE PYDEC
C--------------------------------------------------
c     modified version of DECSIB using
c     PYTHIA decay routine PYDECY \FR'14
C--------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_csydec_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_plist1_cmmn.inc'
      INCLUDE 'sib_parto_cmmn.inc'

c$$$      COMMON /S_CSYDEC/ CBR(223+22), KDEC(1338+6*22), LBARP(99), IDB(99)
c$$$      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
c$$$      COMMON /S_PLIST1/ LLIST1(8000)
c$$$      COMMON /S_MASS1/ AM(99), AM2(99)
c$$$      COMMON /S_PARTO/ NFORIG(8000), NPORIG(8000), IPFLAG, NINT
c$$$      COMMON /S_DEBUG/ Ncall, Ndebug, LUN
      DIMENSION  LL(10), PD(10,5)  
      COMMON/PYJETS/N,NPAD,Klist(4000,5),dPpY(4000,5),dV(4000,5)

      NN = 1
      DO J=1,NP
         LLIST1(J) = 0
      ENDDO
      DO WHILE (NN .LE. NP)
         L= LLIST(NN)
         LA = IABS(L)
         if(LA.lt.100) then
           IF (IDB(LA) .GT. 0)  THEN
              DO K=1,5
                 dppy(1,K) = P(NN,K)
              ENDDO
              if(llist1(NN).eq.0)then
                 klist(1,1) = 1
              else
                 klist(1,1) = 11
              endif
              N = 1
              klist(1,2) = isib_pid2pdg( l )
              CALL PYDECY (1)
              LLIST(NN) = LLIST(NN)+ISIGN(10000,LLIST(NN))
c              print *, 'decay products:',N
c             call pylist(1)
              jj = 2
              DO WHILE ( jj.le.N )
 100             continue
c                 print *,'jj',jj
c                 print *,'decay product:',jj,klist(jj,2)
c                 print *,'py id:',klist(jj,2)
c                 print *,'sib id:',isib_pdg2pid( klist(jj,2) )
                 lldecprod = isib_pdg2pid( klist(jj,2) )
                 IF(lldecprod .eq. 0)THEN
c     intermediate state unknow to sibyll, must decay
c     first store rest of decay products at an offset
c                    print *,'unknown intermediate!',jj,klist(jj,2)
c                    call pylist(1)
                    if(ndebug.ge.5)then
                       write(lun,*) 'unknown intermediate!',
     +                      jj,klist(jj,2)
                       call pylist(1)
                    endif
c     decay intermediate state
                    CALL PYDECY (JJ)
c     overwrite unknown intermediate with one of its decay products
c     N is last such final state particle
                    do nn=1,3
                       klist(jj,nn) = klist(N,nn)
                    enddo
                    do kk=1,5
                       dppy(jj,kk)=dppy(N,kk)
                    enddo
                    N=N-1
                    if(ndebug.ge.5) call pylist(1)
c                    call pylist(1)
c                    print *,'N:',N
                    goto 100
                 ENDIF
                 ll(jj-1) = lldecprod
                 do kk=1,5
                    PD(jj-1,kk) = dppy(jj,kk)
                 enddo
c                 print *,'wrote to cache: jj,ll,n',jj-1,ll(jj-1),n
                 jj = jj + 1
              ENDDO
              ND = N-1
              DO J=1,ND
                NP = NP+1
                if(NP.gt.8000) then
                  write(LUN,'(1x,a,2i8)') 
     &              'PYDEC: no space left in S_PLIST (NP,ND):',NP,ND
                  NP = NP-1
                  return
                endif
                DO K=1,5
                  P(NP,K) = PD(J,K)
                ENDDO
c                print *,'writing to stack:',J,LL(J)
c                print*,'parent:',NN,LLIST(NN)
c                print *,'next on stack:',NN+1,LLIST(NN+1)
c                print *,'total:np,nd',np,nd
                LLIST(NP)=LL(J)
                LLIST1(NP)=NN
                NPORIG(NP)= NPORIG(NN)
                niorig(NP)= NIORIG(NN)
                NFORIG(NP) = L
c                if(ll(J).eq.0) stop
              ENDDO
           ENDIF
         endif
         NN = NN+1
      ENDDO

      END


      SUBROUTINE PYDEC_INI
C****************************************
C     routine to initialize particle 
C     decay tables in PYTHIA       \FR'14
C****************************************
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun

      CALL PDG_INI

      write(LUN,*) '-----------------------------------------'
      write(LUN,*) 'PYDECINI: setting particles stable '
      write(LUN,*) '-----------------------------------------'

c     define stable particles
      do i=1,14
         call py_set_stable( i )
      enddo

c     K0s stable
      call py_set_stable( 12 )

c     Lambda/anti lambda stable
      call py_set_stable( 39 )

c     Sigmas stable
      do i=34,36
         call py_set_stable( i )
      enddo

C     Eta stable
c      call py_set_stable( 23 )
      
      end


      subroutine py_set_stable( ipid )
C****************************************
C     routine to define particle as stable 
C     within PYTHIA.
C     uses sibyll particle labels
C     requires pdg_ini             \FR'14
C****************************************

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      COMMON /S_DEBUG/ Ncall, Ndebug, LUN
      COMMON /S_CNAM/ NAMP (0:99)
      CHARACTER NAMP*6
      INTEGER PYCOMP
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)

      ipdg = isib_pid2pdg( ipid )
      ipyth = pycomp( ipdg )
      mdcy( ipyth , 1 ) = 0
      IF(ndebug.gt.0)
     +     write(LUN,*) 'setting ',namp(i),ipdg,ipyth, 
     +     ' stable in PYTHIA'

      end
