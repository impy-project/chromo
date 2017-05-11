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

      DIMENSION  LL(10), PD(10,5)  
      COMMON/PYJETS/N,NPAD,Klist(4000,5),dPpY(4000,5),dV(4000,5)
      integer pycomp
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
c     check if particle known to pythia
              LPY = PYCOMP(klist(1,2))
c     if unknown, skip and decay with sibyll
              IF(LPY.eq.0) THEN
                 IF(ndebug.ge.5)
     &                write(lun,*) 'particle unknow to pythia! ' ,
     &                lpy, klist(1,2) , 'skipping..'
                 NN = NN + 1
                 cycle
              endif
              IF(ndebug.gt.0)THEN
                 write(lun,*) 'PYDEC: pythia stack prior to decay:'
                 call pylist(1)
              ENDIF
              CALL PYDECY (1)
              LLIST(NN) = LLIST(NN)+ISIGN(10000,LLIST(NN))
              IF(ndebug.gt.0)THEN
                 write(lun,*) 'initial pythia decay done. pylist:'
                 call pylist(1)
              ENDIF
              IF(ndebug.ge.5)
     &             write(lun,*)
     &             'checking final state for unknown..'
              jj = 2
              DO WHILE ( jj.le.N )
 100             continue
                 lldecprod = isib_pdg2pid( klist(jj,2) )
c     intermediate state unknow to sibyll, must decay
c     first store rest of decay products at an offset
                 IF(lldecprod .eq. 0)THEN
                    if(ndebug.ge.5)then
                       write(lun,*) ' unknown intermediate! (j,id)',
     +                      jj,klist(jj,2)
                       call pylist(1)
                       write(lun,*) ' trying decay..'
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
                    if(ndebug.ge.5) THEN
                       write(lun,*) 'stack after decay of intermediate:'
                       call pylist(1)
                    ENDIF
                    goto 100
                 ENDIF
                 ll(jj-1) = lldecprod
                 do kk=1,5
                    PD(jj-1,kk) = dppy(jj,kk)
                 enddo
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
                LLIST(NP)=LL(J)
                LLIST1(NP)=NN
                NPORIG(NP)= NPORIG(NN)
                niorig(NP)= NIORIG(NN)
                NFORIG(NP) = L
              ENDDO
           ENDIF
         endif
         NN = NN+1
      ENDDO

      END


      SUBROUTINE PYDEC_INI
C****************************************
C     routine to initialize particle 
C     decay tables in PYTHIA       \FR'15
C****************************************
      IMPLICIT NONE
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_csydec_cmmn.inc'
      INTEGER i
      CALL PDG_INI

      write(LUN,*) '-----------------------------------------'
      write(LUN,*) 'PYDECINI: using definitions from DECINI! '
      write(LUN,*) '-----------------------------------------'

c     define stable particles
      do i=1,99
         if(idb(i).lt.0)then
            call py_set_stable( i )
         endif
      enddo
      
      end


      subroutine py_set_stable( ipid )
C****************************************
C     routine to define particle as stable 
C     within PYTHIA.
C     uses sibyll particle labels
C     requires pdg_ini             \FR'14
C****************************************
      IMPLICIT NONE
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cnam_cmmn.inc'
c     external
      INTEGER ipid
c     internal
      INTEGER ipdg,ipyth,isib_pid2pdg      
      CHARACTER CODE*18
c     pythia link
      INTEGER PYCOMP,MDCY,MDME,KFDP
      DOUBLE PRECISION BRAT
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)

      ipdg = isib_pid2pdg( ipid )
      ipyth = pycomp( ipdg )
      if(ndebug.gt.0)then
         code = '                  '
         code(1:6) = namp(iabs(ipid))
         if (ipid .lt. 0) code(7:9) = 'bar'
         write(lun,*) 'setting ',code,ipdg,ipyth, 
     &        ' stable in PYTHIA'
      endif
      mdcy( ipyth , 1 ) = 0
      end
