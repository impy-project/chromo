C*********************************************************************
 
      subroutine upyth(id1,iso31,id2,iso32,sqrts)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER id1,iso31,id2,iso32,iso3,kfa
      integer pdgid
      double precision sqrts,mm_to_fmc

      real*8 ranf
      external ranf

c      character*4 reffram 

      include 'newpart.f'
      include 'options.f'

C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYDATR/MRPY(6),RRPY(100)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/,
     &     /PYINT1/,/PYDATR/
      integer pyk,pychge,pycomp
C...Local arrays.
c      DIMENSION PSUM(5),PINI(6),PFIN(6)

C...List of string ends besides those quarks or diquarks from nucleons      
      DIMENSION NSTREND(20)
      DIMENSION NENTRY(4000),SUPPFACT(4000)
      character chaup_pro*20,chaup_tar*20
 
                 
ccccccccccccccccccccccccc section 1: set up PYTHIA cccccccccccccccccccccc

c conversion from mm to fm/c
      mm_to_fmc=1d3/(1d-15/2.99792458d8)

c... redo pythia if scattering was soft
 437  continue

C...Reset process type, kinematics cuts, and the flags used.
      DO 190 ISUB=1,500
         MSUB(ISUB)=0
 190  CONTINUE
      
      msel=0
      msub(11)=1 !default processes
      msub(12)=1
      msub(13)=1
      msub(28)=1
      msub(53)=1
      msub(68)=1
      msub(95)=1
c      msub(81)=1 !extra
c      msub(82)=1
      msub(86)=0                ! turn on J/psi
      msub(87)=0
      msub(88)=0
      msub(89)=0
      msub(106)=0

c... use PYTHIA standard

      MSTP(2)=1  ! 1st order alpha_s

CC use default values
c      MSTP(33)=1 ! K-factor
c      PARP(31)=1.5
c      MSTP(81)=1 ! allow multiple interactions
c      MSTP(82)=1 ! gaussian matter distribution in nucleon
c      MSTP(111)=1 ! do fragmentation (LUND)
c      MSTP(131)=0 ! generate only one event at a time
c      MSTP(133)=0

c... no decays
      mstj(21)=0

c... suppress output
      mstp(122)=0

c... set the random number seed to urqmd random number
      mrpy(1)=int(900000000*ranf(0))

ccccccccccccccccccc section 2: convert urqmd input to PDGID ccccccccccccccc

c      write(6,*) 'id and iso ',id1,iso31,id2,iso32

      kfpro=pdgid(id1,iso31)
      kftar=pdgid(id2,iso32)

      call hepnam(kfpro,chaup_pro)
      call hepnam(kftar,chaup_tar)
ccccccccccccccccccc section 3: call PYTHIA cccccccccccccccccccccccccccccccc

c      write(6,*) 'calling pythia with ',chaup_pro,chaup_tar,sqrts

      CALL PYINIT('CMS',chaup_pro,chaup_tar,sqrts)
 
C...Generate nevent events of each required type.

      CALL PYEVNT

c ... hard scattering? Q>1.5 GeV
      if (vint(51).le.1.5) goto 437

c... define string ends for leading particle cross sections
C...Get the potential string end from the hard process      
      istrend=0
      if(ABS(K(3,2)).le.3)then
         if(K(5,2).ne.K(3,2))then
            istrend = istrend + 1
            nstrend(istrend)=3
         else if(K(7,2).eq.K(3,2))then
            istrend = istrend + 1
            nstrend(istrend)=7
         endif
      endif
      if(ABS(K(4,2)).le.3)then
         if(K(6,2).ne.K(4,2))then
            istrend = istrend + 1
            nstrend(istrend)=4
         else if(K(8,2).eq.K(4,2))then
            istrend = istrend + 1
            nstrend(istrend)=8
         endif
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...Some other cases:      
C...1). gluon fusion or quark fusion? ~1% 
C       or quark-antiquark annilation ~0.1%
c      if(K(5,2).ne.K(7,2).or.K(6,2).ne.K(8,2))then
c         write(*,*)i
c         call pylist(1)
c      endif
C...2). gluon shower to a quark. ~2%
c      nsig = 0
c      if(istrend.eq.0)then
c         do kk=3,8
c            if(K(kk,2).ne.21)nsig = 1
c         enddo
c      endif
c      if(K(5,2).ne.K(7,2).or.K(6,2).ne.K(8,2))nsig=0
c      if(nsig.eq.1)then
c         write(*,*)i
c         call pylist(1)
c      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C...Find the entries of the output hadrons (entries with K(I,1)=1) after
C   the call of pyedit(2).
      nhadron=0
      do j=1,N
         if(K(j,1).eq.1)then
            nhadron=nhadron+1
            NENTRY(j)=nhadron
         else
            NENTRY(j)=0
         endif
      enddo

C...Initialize of the array of leading factor.
      do j=1,nhadron
         SUPPFACT(j)=0
      enddo

C...Get the strings
      istrsig=0
      nenstrbegin=0
      nenstrend=0
      do j=msti(4)+1,N
         if(K(j,1).eq.12.and.istrsig.eq.0)then
            !begin of a string
            istrsig=1
            nenstrbegin=j

         else if(K(j,1).eq.11.and.istrsig.eq.1)then
            !end of a string
            istrsig=0
            nenstrend=j

            if(K(j,3).le.msti(4))then
            !consider only the original string
            !firtly, get the entry of the string
            ndt=K(j,4)
 1099       if(K(ndt,2).ne.92.and.K(ndt,2).ne.91)then
               ndt=K(ndt,4)
               goto 1099
            endif

            !find the leading parton                     
            ileadsigl=0
            if(K(nenstrbegin,3).le.2.and.K(nenstrbegin,3).gt.0)then
               ileadsigl=1
            else
               do jj=1,istrend
               if(K(nenstrbegin,3).eq.nstrend(jj)
     $            .and.K(nenstrbegin,2).eq.K(nstrend(jj),2))then
                  ileadsigl=1
               endif
               enddo
            endif
            ileadsigr=0
            if(K(nenstrend,3).le.2.and.K(nenstrend,3).gt.0)then
               ileadsigr=1
            else
               do jj=1,istrend
               if(K(nenstrend,3).eq.nstrend(jj)
     $            .and.K(nenstrend,2).eq.K(nstrend(jj),2))then
                  ileadsigr=1
               endif
               enddo
            endif
            
            if(K(ndt,2).eq.91)then
               !cluster
               !has 1 or 2 daughters
               !ASSUME:both quarks in a di-quark are from the incident hadrons

               if(abs(K(nenstrbegin,2)).lt.1000
     $            .and.abs(K(nenstrend,2)).lt.1000)then
                  !meson string (T)
                  if(K(ndt,4).eq.K(ndt,5))then
                     !single meson
                     if(ileadsigl.eq.1)then
                        SUPPFACT(NENTRY(K(ndt,4)))=
     $                    SUPPFACT(NENTRY(K(ndt,4)))+0.5d0
                     endif
                     if(ileadsigr.eq.1)then
                        SUPPFACT(NENTRY(K(ndt,4)))=
     $                    SUPPFACT(NENTRY(K(ndt,4)))+0.5d0
                     endif
                  else
                     !double meson
                     if(ileadsigl.eq.1)SUPPFACT(NENTRY(K(ndt,4)))=0.5d0
                     if(ileadsigr.eq.1)SUPPFACT(NENTRY(K(ndt,5)))=0.5d0
                  endif
               else
                  !baryon string (T)
                  if(abs(K(nenstrend,2)).gt.1000)then
                     nenstrtmp=nenstrbegin
                     nenstrbegin=nenstrend
                     nenstrend=nenstrtmp
                     ileadsigtmp=ileadsigl
                     ileadsigl=ileadsigr
                     ileadsigr=ileadsigtmp
                  endif

                  if(K(ndt,4).eq.K(ndt,5))then
                     !single baryon
                     if(ileadsigl.eq.1)then
                        SUPPFACT(NENTRY(K(ndt,4)))=
     $                    SUPPFACT(NENTRY(K(ndt,4)))+0.66667d0
                     endif
                     if(ileadsigr.eq.1)then
                        SUPPFACT(NENTRY(K(ndt,4)))=
     $                    SUPPFACT(NENTRY(K(ndt,4)))+0.33333d0
                     endif
                  else
                     !baryon+meson
                     if(ileadsigl.eq.1)then
                        if(abs(K(K(ndt,4),2)).gt.1000)then
                           SUPPFACT(NENTRY(K(ndt,4)))=0.66667d0
                        else
                           SUPPFACT(NENTRY(K(ndt,4)))=0.5d0
                           SUPPFACT(NENTRY(K(ndt,5)))=0.33333d0
                        endif
                     endif
                     if(ileadsigr.eq.1)then
                        if(abs(K(K(ndt,5),2)).gt.1000)then
                           SUPPFACT(NENTRY(K(ndt,5)))=
     $                       SUPPFACT(NENTRY(K(ndt,5)))+0.33333d0
                        else
                           SUPPFACT(NENTRY(K(ndt,5)))=0.5d0
                        endif
                     endif
                  endif
               endif
            else if(K(ndt,2).eq.92)then
               !string
               !has >= 2 daughters
               !ASSUME:both quarks in a di-quark are from the incident hadrons

               if(ileadsigl.eq.1)then
                  if(abs(K(nenstrbegin,2)).lt.10)then
                     !leading quark or anti-quark
                     if(abs(K(K(ndt,4),2)).lt.1000)then
                        SUPPFACT(NENTRY(K(ndt,4)))=0.5d0
                     else
                        SUPPFACT(NENTRY(K(ndt,4)))=0.33333d0
                     endif
                  else if(abs(K(nenstrbegin,2)).gt.1000)then
                     !leading di-quark
                     if(abs(K(K(ndt,4),2)).lt.1000)then
                        SUPPFACT(NENTRY(K(ndt,4)))=0.5d0
                        SUPPFACT(NENTRY(K(ndt,4)+1))=0.33333d0
                     else
                        SUPPFACT(NENTRY(K(ndt,4)))=0.66667d0
                     endif
                  endif
               endif
               
               if(ileadsigr.eq.1)then
                  if(abs(K(nenstrend,2)).lt.10)then
                     !leading quark or anti-quark
                     if(abs(K(K(ndt,5),2)).lt.1000)then
                        SUPPFACT(NENTRY(K(ndt,5)))=
     $                     SUPPFACT(NENTRY(K(ndt,5)))+0.5d0
                     else
                        SUPPFACT(NENTRY(K(ndt,5)))=
     $                     SUPPFACT(NENTRY(K(ndt,5)))+0.33333d0
                     endif
                  else if(abs(K(nenstrend,2)).gt.1000)then
                     !leading di-quark
                     if(abs(K(K(ndt,5),2)).lt.1000)then
                        SUPPFACT(NENTRY(K(ndt,5)))=
     $                     SUPPFACT(NENTRY(K(ndt,5)))+0.5d0
                        SUPPFACT(NENTRY(K(ndt,5)-1))=
     $                     SUPPFACT(NENTRY(K(ndt,5)-1))+0.33333d0
                     else
                        SUPPFACT(NENTRY(K(ndt,5)))=
     $                     SUPPFACT(NENTRY(K(ndt,5)))+0.66667d0
                     endif
                  endif
               endif

            else 
               write(*,*)'ERROR in PYTHIA: ',i,j,ndt,K(ndt,2)
               call pylist(2)
            endif
            endif
      
         else if(K(j,1).eq.1.and.K(j,3).le.2)then
            !direct hadrons
            if(abs(K(j,2)).lt.1000)then
               !meson
               SUPPFACT(NENTRY(j))=0.5d0
            else
               !baryon
               !Need to check whether there are exceptions(only one l.q.)
               SUPPFACT(NENTRY(j))=0.66667d0
            endif

         else if(K(j,1).eq.11.and.K(j,2).eq.92.and.istrsig.eq.0)then
            !the string entry. the end of the treatment.
            goto 1101
         endif
      enddo

 1101 continue

c... check out the suppression factor before the removement
c      call pylist(2)
c      do j=1,N
c         if(K(j,1).eq.1.and.SUPPFACT(NENTRY(j)).ne.0)then
c         write(*,*)j,SUPPFACT(NENTRY(j))
c         endif
c      enddo
      
c... remove unwanted particles from event
      call pyedit(2)

 
C...List statistics for each process type.
c     CALL PYSTAT(1)


ccccccccccccccccccc section 4: convert PYTHIA arrays to UrQMD newpart array
c-
      nexit=n
      nstring1=nexit
      nstring2=0

c      write(6,*) 'exit pythia with ',nexit

      do 10 i=1,nexit

c         write(6,*) ' particle ',k(i,2)
         kfa=k(i,2)
         call pdg2id(ityp,iso3,kfa)
         itypnew(i)=ityp
         i3new(i)=iso3

c set leadfac
c.. suppress some part of the cross section to mimic coherence
         SUPPFACT(i)=SUPPFACT(i)*CTParam(59)
         leadfac(i)=SUPPFACT(i)

         do 20 j=1,4
            pnew(j,i)=p(i,j)
            xnew(j,i)=v(i,j)*mm_to_fmc !useless
 20      continue
         pnew(5,i)=p(i,5)
c set formation time
         if(pnew(5,i).gt.0d0) then
            xnew(4,i)=0.8d0*pnew(4,i)/pnew(5,i)
         endif

 10   continue


 
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pdg2id (ityp,iso3,pdgid)
c
c     Author   : Steffen A. Bass
c     Date     : 06/08/98
c     Revision : 1.0
c
c based on ityp2pdg.f from Henning Weber
c 
c     converts PDG-Id to  UrQMD-Id 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none


      integer pdgid,i
      integer ityp
      integer iso3

      integer tab_size
      parameter (TAB_SIZE = 172)

      integer idtab(3,TAB_SIZE)

      character chau*20

      data idtab/
c Neutron
     .       1, -1,  2112,  
c Proton
     .       1,  1,  2212,
c N*
     .       2, -1, 12112, 
     .       2,  1, 12212,
     .       3, -1,  1214, 
     .       3,  1,  2124, 
     .       4, -1, 22112, 
     .       4,  1, 22212,
     .       5, -1, 32112,
     .       5,  1, 32212,
     .       6, -1,  2116,
     .       6,  1,  2216,
     .       7, -1, 12116,
     .       7,  1, 12216,
     .       8, -1, 21214,
     .       8,  1, 22124,
     .       9, -1, 42112, 
     .       9,  1, 42212, 
     .      10, -1, 31214, 
     .      10,  1, 32124, 
     .      14, -1,  1218, 
     .      14,  1,  2128, 
c Delta
     .      17, -3,  1114,
     .      17, -1,  2114,
     .      17,  1,  2214,
     .      17,  3,  2224,
     .      18, -3, 31114,
     .      18, -1, 32114,
     .      18,  1, 32214,
     .      18,  3, 32224,
     .      19, -3,  1112,
     .      19, -1,  1212,
     .      19,  1,  2122,
     .      19,  3,  2222,
     .      20, -3, 11114,
     .      20, -1, 12114,
     .      20,  1, 12214,
     .      20,  3, 12224,
     .      21, -3, 11112,
     .      21, -1, 11212,
     .      21,  1, 12122,
     .      21,  3, 12222,
     .      22, -3,  1116,
     .      22, -1,  1216,
     .      22,  1,  2126,
     .      22,  3,  2226,
     .      23, -3, 21112,
     .      23, -1, 21212,
     .      23,  1, 22122,
     .      23,  3, 22222,
     .      24, -3, 21114,
     .      24, -1, 22114,
     .      24,  1, 22214,
     .      24,  3, 22224,
     .      25, -3, 11116,
     .      25, -1, 11216,
     .      25,  1, 12126,
     .      25,  3, 12226,
     .      26, -3,  1118,
     .      26, -1,  2118,
     .      26,  1,  2218,
     .      26,  3,  2228,
c Lambda
     .      27,  0,  3122,
     .      28,  0, 13122,   
     .      29,  0,  3124,   
     .      30,  0, 23122,   
     .      31,  0, 33122,
     .      32,  0, 13124,
     .      33,  0, 43122,   
     .      34,  0, 53122,   
     .      35,  0,  3126,   
     .      36,  0, 13126,   
     .      37,  0, 23124,   
     .      38,  0,  3128,   
     .      39,  0, 23126,   
c Sigma
     .      40, -2,  3112,
     .      40,  0,  3212,
     .      40,  2,  3222,
     .      41, -2,  3114,
     .      41,  0,  3214,
     .      41,  2,  3224,
     .      42, -2, 13112,
     .      42,  0, 13212,
     .      42,  2, 13222,
     .      43, -2, 13114,
     .      43,  0, 13214,
     .      43,  2, 13224,
     .      44, -2, 23112,
     .      44,  0, 23212,
     .      44,  2, 23222,
     .      45, -2,  3116,
     .      45,  0,  3216,
     .      45,  2,  3226,
     .      46, -2, 13116,
     .      46,  0, 13216,
     .      46,  2, 13226,
     .      47, -2, 23114,
     .      47,  0, 23214,
     .      47,  2, 23224,
     .      48, -2,  3118,
     .      48,  0,  3218,
     .      48,  2,  3228,
c Xi
     .      49, -1,  3312,
     .      49,  1,  3322,
     .      50, -1,  3314,
     .      50,  1,  3324,
     .      52, -1, 13314,
     .      52,  1, 13324,
c Omega
     .      55,  0,  3334,
c gamma
     .     100,  0,    22, 
c pion
     .     101, -2,  -211,
     .     101,  0,   111,
     .     101,  2,   211,
c eta
     .     102,  0,   221,
c omega
     .     103,  0,   223,
c rho
     .     104, -2,  -213,
     .     104,  0,   113, 
     .     104,  2,   213,
c f0(980)
     .     105,  0, 10221,
c kaon
     .     106, -1,   311,
     .     106,  1,   321,
c eta'
     .     107,  0,   331,
c k*(892)
     .     108, -1,   313,
     .     108,  1,   323,
c phi
     .     109,  0,   333,
c k0*(1430)
     .     110, -1, 10313,
     .     110,  1, 10323,
c a0(980)
     .     111, -2,-10211,
     .     111,  0, 10111,
     .     111,  2, 10211,
c f0(1370)
     .     112,  0, 20221,
c k1(1270)
     .     113, -1, 10313,
     .     113,  1, 10323,
c a1(1260)
     .     114, -2,-20213,
     .     114,  0, 20113,
     .     114,  2, 20213,
c f1(1285)
     .     115,  0, 20223,
c f1'(1510)
     .     116,  0, 40223,
c k2*(1430)
     .     117, -1,   315,
     .     117,  1,   325,
c a2(1329)
     .     118, -2,  -215,
     .     118,  0,   115,
     .     118,  2,   215,
c f2(1270)
     .     119,  0,   225,
c f2'(1525)
     .     120,  0,   335,
c k1(1400)
     .     121, -1, 20313,
     .     121,  1, 20323,
c b1
     .     122, -2,-10213,
     .     122,  0, 10113,
     .     122,  2, 10213,
c h1
     .     123,  0, 10223,
c K* (1410)
     .     125, -1, 30313,
     .     125,  1, 30323,
c rho (1450)
     .     126, -2,-40213,
     .     126,  0, 40113,
     .     126,  2, 40213,
c omega (1420)
     .     127,  0, 50223,
c phi(1680)
     .     128,  0, 10333,
c k*(1680)
     .     129, -1, 40313,
     .     129,  1, 40323,
c rho(1700)
     .     130, -2,-30213,
     .     130,  0, 30113,
     .     130,  2, 30213,
c omega(1600)
     .     131,  0, 60223,
c phi(1850)     
     .     132,  0,   337,
c D
     .     133,  1,   411,    
     .     133, -1,   421,  
c D*
     .     134,  1, 10411,    
     .     134, -1, 10421,  
c J/Psi
     .     135,  0, 443,
c Psi'
     .     136,  0, 100443,
c Chi_C
     .     137,  0, 10441/ 


c search for the ITYP in IDTAB


      do 99 i=1,TAB_SIZE
         if(idtab(3,i).eq.pdgid) then
            ityp=idtab(1,i)
            iso3=idtab(2,i)
            return
         elseif((pdgid.lt.0).and.(idtab(3,i).eq.abs(pdgid))) then
            ityp=(-1)*idtab(1,i)
            iso3=(-1)*idtab(2,i)
            return
         endif
 99   continue



c handle all unknown ityps from PYTHIA
      ityp=1000+abs(pdgid)
      if(pdgid.lt.0) ityp=(-1)*ityp
      iso3=0


c      call hepnam(pdgid,chau)  
c      write(6,*) '(Info) pdg2ityp: ',chau,' PDG-Id: ',pdgid,' 
c     &  -> ityp',ityp
c      stop 137

      return
      end




 
