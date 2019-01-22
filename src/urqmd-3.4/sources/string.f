c $Id: string.f,v 1.23 2007/01/30 14:50:27 bleicher Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stringdec(ityp,iz2,smass,part,ident2,npart)
c
cinput smass   : Stringmass
cinput ityp    : Particle ID
cinput iz2     : Isospin$_3\cdot 2$
c
coutput part   : 4-momenta, 4-position, masses (array)
coutput ident2 : ityp, iz2 (array)
coutput npart  : number of outgoing particles
c
c     This subroutine performs string fragmentation.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      PARAMETER(MXPTCL=200)
      COMMON/PARTCL/ PPTCL(9,MXPTCL),nptcl,IDENT(MXPTCL),IDCAY(MXPTCL)

      include 'comstr.f'

      dimension part(9,mxptcl)
      dimension ident2(2,mxptcl)


c we call the translation routine
        call ityp2id(ityp,iz2,ifa,ifb)

        goto 1
cspl... string fragmentation called with quark id's and energy as arguments
         entry qstring(ifanew,ifbnew,smass,part,ident2,npart)
         ifa=ifanew
         ifb=ifbnew
 1      continue


        smem=smass
c here we call the fragmentation routine. the produced hadrons and their
c properties are returned via the pptcl- and ident-array in the common-block
       call string(ifa,ifb,smass)
c here the array pptcl has been filled with nptcl entries (1->nptcl)
c now we translate to uqmd and shift the pptcl- and ident-info to the
c corresponding part- and ident2-arrays of the newpart-common-block:
       npart=nptcl       
       do 2 i=1,nptcl
        call id2ityp(ident(i),pptcl(5,i),itypout,iz2out)
        ident2(1,i)=itypout
        ident2(2,i)=iz2out
        smem=smem-pptcl(4,i)
        do 3 j=1,9
         part(j,i)=pptcl(j,i)
 3       continue
 2      continue    
c check for energy conservation:
       if(abs(smem).gt.1.0D-5)then
        write(*,*)'! stringdec: energy difference=',smem 
        write(*,*)'ifa,ifb,smass=',ifa,ifb,smass
       endif
       return

       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine strini
c
c output     : via common blocks
c
c     {\tt strini} calculates mixing angles for the meson-multipletts
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      include 'options.f'
        include 'comres.f'
      include 'comstr.f'
        real*8 m3
        real*8 massit
        integer jit
c mixing angles of meson multiplets according to flavor SU(3) quark model:
cspl-0795 these parameters assign the pure u/ubar,d/dbar,s/sbar states
c  (e.g. 110,220,330) to the physical particles according to the su(3)
c  quark model. The flavor mixing angles are chosen according to quadratic
c  Gell-Mann-Okubo mass formula (Review of Particle Properties,
c  Phys. Rev D50 (1994) 1319). For the scalar mesons this formula is not 
c  applicable. We assume an ideal mixing angle (tan(theta)=1/sqrt(2)).
c
c pseudoscalar: theta=-10 deg
c vector      : theta= 39 deg
c pseudovector: theta= 51 deg
c tensor      : theta= 28 deg
c 
      real*8 mixang(njspin)
c ideal mixing angles assumed for the last three multiplets
      data mixang/-10d0,39d0,35.3d0,51d0,28d0,35.3d0,35.3d0,35.3d0/
c 

      pi = 4d0*datan2(1d0,1d0)

c calculate 'singlet shift probabilities', e.g. a 11 (u-ubar) state
c can be changed to a 22 (d-dbar) or 33 (s-sbar) state with a certain
c probability. THEN they can be identified with physical hadrons!
      do 3 i=1,njspin
         mixang(i)=mixang(i)/36d1*2d0*3.1416d0
         PMIX1S(1,i)=(dcos(mixang(i))/sqrt(6d0)
     &        -dsin(mixang(i))/sqrt(3d0))**2
         PMIX1S(2,i)=PMIX1S(1,i)
         PMIX1S(3,i)=(-(dcos(mixang(i))*2d0/dsqrt(6d0))
     &        -dsin(mixang(i))/dsqrt(3d0))**2
         PMIX2S(1,i)=0.5d0
         PMIX2S(2,i)=0.5d0
         PMIX2S(3,i)=1d0
         
ce calculate probabilities of the meson multipletts
ce according to parm=(spin degeneracy)/(average mass) *ctp(50 ff.)
           parm(i)=0d0
           m3=0d0
           do 102 j=0,3
             itp=mlt2it(4*(i-1)+j+1)
             m3=m3+massit(itp)
             jpc=jit(itp)/2
102        continue
           parm(i)=parm(i)+(2*jpc+1)/m3*4*CTParam(49+i)
         
c the mixing-angles are the same for 'string' and 'cluster':
         do 2 k=1,3            
            PMIX1C(k,i)=PMIX1S(k,i)
            PMIX2C(k,i)=PMIX2S(k,i)
c            write(6,*)'#',pmix1c(k,i)
 2       continue
 3    continue
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE GAUSPT(PT0,SIGQT)
c
cinput   sigqt  : Width of Gaussian
c
coutput  pt0    : transverse momentum 
c
C     generate pt with Gaussian
c     distribution $\propto pt \exp(-pt^2/sigqt^2)$
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real*8 pt0,sigqt,rnd,ranf

      RND=ranf(0)
      PT0=SIGQT*SQRT(-DLOG(1.d0-RND))
      RETURN
      END


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN)
c
cinput ID     : quarkcode
c
coutput ifl1  : single quarks id 
coutput ifl2  : single quarks id 
coutput ifl3  : single quarks id 
coutput jspin : spin id 
c
C          THIS SUBROUTINE UNPACKS THE IDENT CODE ID=+/-IJKL
c
C          -MESONS:
C          I=0, J<=K, +/- IS SIGN FOR J,
C          ID=110 FOR PI0, ID=220 FOR ETA, ETC.
c
C          -BARYONS:
C          I<=J<=K IN GENERAL,
C          J<I<K FOR SECOND STATE ANTISYMMETRIC IN (I,J), EG. L = 2130
c
C          -DIQUARKS:
C          ID=+/-IJ00, I<J FOR DIQUARK COMPOSED OF I,J.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      IDABS=IABS(ID)
c extract the single (anti-)quark ids from the hadron id:
      I=IDABS/1000
      J=MOD(IDABS/100,10)
      K=MOD(IDABS/10,10)
      JSPIN=MOD(IDABS,10)
c diquarks:
      IF(ID.NE.0.AND.MOD(ID,100).EQ.0) GO TO 300
c quarks oder so:
      IF(J.EQ.0) GO TO 200
c mesonen:
      IF(I.EQ.0) GO TO 100

C..BARYONS:
C..ONLY X,Y BARYONS ARE QQX, QQY, Q=U,D,S.
      IFL1=ISIGN(I,ID)
      IFL2=ISIGN(J,ID)
      IFL3=ISIGN(K,ID)
      RETURN
C          MESONS
100   CONTINUE
      IFL1=0
      IFL2=ISIGN(J,ID)
      IFL3=ISIGN(K,-ID)
      RETURN
200   CONTINUE
      IFL1=0  
      IFL2=0  
      IFL3=0  
      JSPIN=0 
      return
300   IFL1=ISIGN(I,ID)
      IFL2=ISIGN(J,ID)
      IFL3=0
      JSPIN=0
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE STRING(IFL1,IFL2,AMSTR)
c
cinput   amstr     : stringmass
cinput   ifl1      : leading (di)quark (along +Z)
cinput   ifl2      : leading (di)quark
c
c output : produced particles via common block ({\tt pptcl}) 
c
C     Hadron production via string fragmentation. masses acc. to Breit
c     Wigner distr., incl. production of mesonic and baryonic resonances
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      
      COMMON/COMTRY/ NTRIES
      PARAMETER(MXPTCL=200)
      COMMON/PARTCL/ PPTCL(9,MXPTCL),nptcl,IDENT(MXPTCL),IDCAY(MXPTCL)

      include 'comstr.f'

      DIMENSION PX1L(2),PY1L(2),PX1(2),PY1(2),PMTS(2),W(2),IFL(2)
      LOGICAL DIQBR,SPINT,BACK
      DIMENSION VC(3)
      DIMENSION PPTL(5,MXPTCL),PPTR(5,MXPTCL),
     *IDENTL(MXPTCL),IDENTR(MXPTCL)
      COMMON/KAPPA/ XAP
      COMMON/COLRET/ LRET
      LOGICAL LRET
      COMMON/CONSTI/ CONSTI
      LOGICAL CONSTI

      LOGICAL leading

      include 'options.f'

      IFLL=0
      PXL=0
      PYL=0
c..strange quark and charm quark probabilities:
      prbs=ctparam(6)
      prbc=ctparam(7)
cspl the diquark-suppresion parameter is reduced for small
c    string masses (finite size effect) see A. Jahns diploma thesis
      if(amstr.le.2.8d0)then
        pbars=0.d0
      elseif (amstr.le.5d0)then
        pbars=((amstr-2.8d0)/2.2d0)**3*ctparam(8)
      else
        pbars=ctparam(8)
      endif

c..some parameters (see input.f)
      dmas=ctparam(9)
      dmass=ctparam(10)
      pardbs=ctparam(25)
      sigqts=ctparam(42)
C
C  PRBS STRANGENESS SUPPRESSION PARAMETER
C  PRBC CHARM SUPPRESSION PARAMETER
      PU=1./(2.+PRBS+PRBC)
C
      LRET=.FALSE.
C
      NREP = 0
      DIQBR=.TRUE.
      NFIX=0
      BACK=.TRUE.

c.. here starts everything (if the string break up didn't work
c..                         we start again at this point)
 100  I=NFIX
      ilead=0
      NPR=0
      NPL=0
      NPTCL=NFIX
c.. nrep=number of tries to break string
      NREP=NREP+1
      IF(NREP.LT.NTRIES) GO TO 102
      LRET=.TRUE.
      write(*,*)'!! STRING: no fragmentation.'
      return
102   CONTINUE
      IFL(1)=IFL1
      IFL(2)=IFL2
      DO 110 J=1,2
 110  W(J)=AMSTR
      DO 120 J=1,2
      PX1L(J)=0.
      PY1L(J)=0.
      PX1(J)=0.
 120  PY1(J)=0.
C  WILL THERE BE ONLY ONE BREAK OR NOT ?
      SPINT=.TRUE.
      KSPIN=1
      IF(MOD(IFL(1),100).EQ.0.AND.MOD(IFL(2),100).EQ.0) GOTO 131
      IDR=IDPARS(IFL(1),IFL(2),SPINT,KSPIN)
      WEND=(AMASS(IDR)+DMAS)**2
      GO TO 151
131   IFCN=1
      IF(ranf(0).GT.0.5) IFCN=2
      IFLC1=IFCN
      IF(IFL(1).LT.0) IFLC1=-IFCN
      IDR1=IDPARS(IFL(1),IFLC1,SPINT,KSPIN)
      IDR2=IDPARS(IFL(2),-IFLC1,SPINT,KSPIN)
      WEND=(AMASS(IDR1)+AMASS(IDR2)+DMAS)**2
151   SPINT=.FALSE.
      KSPIN=0
c.. only one break goto 225 is the end of the fragmentation
      IF(W(1)*W(2).LE.WEND) GO TO 225

c..the main iteration loop for the fragmentation:
 130  I=I+1
      IF(I.GT.MXPTCL) GO TO 9999
C  CHOOSE SIDE OF BREAK
      JSIDE=INT(1.+2.*ranf(0))
      IF(JSIDE.EQ.1) NPR=NPR+1
      IF(JSIDE.EQ.2) NPL=NPL+1
      IF(NPR.GT.MXPTCL.OR.NPL.GT.MXPTCL) GO TO 9999
C  IS IFL(JSIDE) A QUARK OR A DIQUARK ? (di-quark-->150)
      IF(MOD(IFL(JSIDE),100).EQ.0) GO TO 150
C  IFL(JSIDE) IS A QUARK
C  NOW WE SELECT Q,QBAR PAIR OR QQ,QQBAR PAIR
      DIQBR=.FALSE.
      DRND=ranf(0)
c.. do a qq-qqbar pair with certain prob. 
      IF(DRND.LT.PBARS) GO TO 140
C  Q,QBAR PAIR
      IFLN=ISIGN(IFLAV(PU,PRBS),-IFL(JSIDE))
      GO TO 200
C  QQ,QQBAR PAIR
140   IQ1=IFLAV(PU,PRBS)
      IQ2=IFLAV(PU,PRBS)
c.. suppr. double strange di-quarks(ss) with certain prob. (acc. to ctp 29)
      if((IQ1.eq.3.and.IQ2.eq.3)
     &     .and.ranf(0).gt.CTParam(29))goto 140
c..  additional qs diquark suppression needed for anti-lambda suppression 
c... single strange di-quark suppresion: 
       if((IQ1.eq.3.and.IQ2.ne.3).or.(IQ1.ne.3.and.IQ2.eq.3)
     &  .and.ranf(0).gt.CTParam(49))goto 140
c
      IF(IQ1.LE.IQ2) GO TO 145
      ISWAP=IQ1
      IQ1=IQ2
      IQ2=ISWAP
145   IFQQ=1000*IQ1+100*IQ2
      IFLN=ISIGN(IFQQ,IFL(JSIDE))
      GO TO 200
c..the di-quark part:
C  IFL(JSIDE) IS A DIQUARK
C  CAN DIQUARK BREAK OR NOT
 150  DRND=ranf(0)
      IF(DRND.LE. PARDBS) GO TO 190
C DIQUARK BREAK (prob. in PARDBS)
      CALL FLAVOR(IFL(JSIDE),IFLD1,IFLD2,IFLD3,JSPIN)
      IFLL=IFLD1
      IFL(JSIDE)=IFLD2
      DRND=ranf(0)
      IF(DRND.GE.PARQLS) GO TO 160
      IFLL=IFLD2
      IFL(JSIDE)=IFLD1
 160  DIQBR=.TRUE.
C  LEADING QUARK TRANSVERSE MOMENTUM
      CALL GAUSPT(PTL0,SIGQTS)
      PHI=2.*PI*ranf(0)
      PXL=PTL0*COS(PHI)
      PYL=PTL0*SIN(PHI)
      PX1L(JSIDE)=PX1(JSIDE)
      PY1L(JSIDE)=PY1(JSIDE)
      PX1(JSIDE)=-PXL
      PY1(JSIDE)=-PYL
C  Q,QBAR PAIR
      IFLN=ISIGN(IFLAV(PU,PRBS),-IFL(JSIDE))
      GO TO 200
C  DIQUARK DOES NOT BREAK
C  Q,QBAR PAIR
 190  IFLN=ISIGN(IFLAV(PU,PRBS),IFL(JSIDE))
      DIQBR=.FALSE.
C  IDENT,MASS AND TRANSVERSE MOMENTUM OF PARTICLE
 200  IDENT(I)=IDPARS(IFL(JSIDE),IFLN,SPINT,KSPIN)
      PPTCL(5,I)=AMASS(IDENT(I))
      SIGQTSN=SIGQTS
      IF(MOD(IFLN,100).EQ.0) SIGQTSN=sigqts*ctparam(38) 
      if(iabs(ifln).eq.3.or.iabs(ifl(jside)).eq.3)
     &                  sigqtsn=sigqts*ctparam(39) 
      CALL GAUSPT(PT2,SIGQTSN)
c no pt for leading hadron:
      leading=.false.
      if((JSIDE.EQ.1.and.NPR.eq.1).and.
     & (abs(IDENT(I)).eq.1120 .or.abs(IDENT(I)).eq.1220) )then
c     & (abs(IDENT(I)).ge.1110))then
        leading=.true.
cblu        pt2=pt2/2d0
        ilead=ilead+1
      endif
      if((JSIDE.EQ.2.and.NPL.eq.1).and.
     & (abs(IDENT(I)).eq.1120 .or.abs(IDENT(I)).eq.1220) )then
c     & (abs(IDENT(I)).ge.1110))then
        leading=.true.
cblu        pt2=pt2/2d0
        ilead=ilead+1
      endif
c..transverse momentum choosen for the newly produced hadron
      PHI=2.*PI*ranf(0)
      PX2=PT2*COS(PHI)
      PY2=PT2*SIN(PHI)
      PPTCL(1,I)=PX1(JSIDE)+PX2
      PPTCL(2,I)=PY1(JSIDE)+PY2
C  GENERATE Z-momentum
      PMTS(3-JSIDE)=AMASS(IABS(IFL(3-JSIDE)))**2
      PTS=PPTCL(1,I)**2+PPTCL(2,I)**2
      PMTS(JSIDE)=PPTCL(5,I)**2+PTS
      IF(PMTS(JSIDE)+PMTS(3-JSIDE).GE.PARRS*W(1)*W(2)) GO TO 100
      ZMIN=PMTS(JSIDE)/(W(1)*W(2))
      ZMAX=1.-PMTS(3-JSIDE)/(W(1)*W(2))
      IF(ZMIN.GE.ZMAX) GO TO 100
C..WARNING: VERY IMPORTANT THE ORDER OF IFL AND IFLN IN ZFRAGS
c.. fraction of momentum acc. to the fragmentation fct.
      Z=ZFRAGS(IFL(JSIDE),IFLN,PTS,ZMIN,ZMAX,leading)

      PPTCL(3,I)=0.5*(Z*W(JSIDE)-PMTS(JSIDE)/
     *(Z*W(JSIDE)))*(-1.)**(JSIDE+1)
      PPTCL(4,I)=0.5*(Z*W(JSIDE)+PMTS(JSIDE)/(Z*W(JSIDE)))
      IDCAY(I)=0
      IF(.NOT.(JSIDE.EQ.1)) GO TO 282
      IDENTR(NPR)=IDENT(I)
      PPTR(1,NPR)=PPTCL(1,I)
      PPTR(2,NPR)=PPTCL(2,I)
      PPTR(3,NPR)=PPTCL(3,I)
      PPTR(4,NPR)=PPTCL(4,I)
      PPTR(5,NPR)=PPTCL(5,I)
 282  IF(.NOT.(JSIDE.EQ.2)) GO TO 283
      IDENTL(NPL)=IDENT(I)
      PPTL(1,NPL)=PPTCL(1,I)
      PPTL(2,NPL)=PPTCL(2,I)
      PPTL(3,NPL)=PPTCL(3,I)
      PPTL(4,NPL)=PPTCL(4,I)
      PPTL(5,NPL)=PPTCL(5,I)
 283  IF(DIQBR) GO TO 210
      IFL(JSIDE)=-IFLN
      PX1(JSIDE)=-PX2
      PY1(JSIDE)=-PY2
      GO TO 220
C  NEW DIQUARK CREATION
210   ID1=IABS(IFLL)
      ID2=IABS(IFLN)
      IF(ID1.LE.ID2) GO TO 215
      ISWAP=ID1
      ID1=ID2
      ID1=ISWAP
215   IFL(JSIDE)=ISIGN(1000*ID1+100*ID2,IFLL)
      PX1L(JSIDE)=PX1L(JSIDE)+PXL-PX2
      PY1L(JSIDE)=PY1L(JSIDE)+PYL-PY2
      PX1(JSIDE)=PX1L(JSIDE)
      PY1(JSIDE)=PY1L(JSIDE)
 220  W(1)=W(1)-PPTCL(4,I)-PPTCL(3,I)
      W(2)=W(2)-PPTCL(4,I)+PPTCL(3,I)
      SPINT=.TRUE.
      KSPIN=2
      IF(MOD(IFL(1),100).EQ.0.AND.MOD(IFL(2),100).EQ.0) GO TO 240
      IDB=IDPARS(IFL(1),IFL(2),SPINT,KSPIN)
      AMB=AMASS(IDB)+dmass
      GO TO 211
240   IFCN=1
      IF(ranf(0).GT.0.5) IFCN=2
      IFLC1=-IFCN
      IF(IFL(1).GT.0) IFLC1=IFCN
      IFLC2=-IFLC1
      IKH1=IDPARS(IFL(1),IFLC1,SPINT,KSPIN)
      IKH2=IDPARS(IFL(2),IFLC2,SPINT,KSPIN)
      AMB=AMASS(IKH1)+AMASS(IKH2)+DMASs
211   P1X=PX1(1)+PX1(2)
      P1Y=PY1(1)+PY1(2)
      PT12=P1X**2+P1Y**2
      W12=W(1)*W(2)
      AMS2=W12-PT12
      IF(AMS2.LT.AMB**2) GO TO 100
      SPINT=.TRUE.
      KSPIN=1
      IF(MOD(IFL(1),100).EQ.0.AND.MOD(IFL(2),100).EQ.0) GO TO 231
      IDR=IDPARS(IFL(1),IFL(2),SPINT,KSPIN)
      WEND=(AMASS(IDR)+DMAS)**2
      GO TO 232
231   IKHR1=IDPARS(IFL(1),IFLC1,SPINT,KSPIN)
      IKHR2=IDPARS(IFL(2),IFLC2,SPINT,KSPIN)
      WEND=(AMASS(IKHR1)+AMASS(IKHR2)+DMAS)**2
232   SPINT=.FALSE.
      KSPIN=0
      IF(W(1)*W(2).GE.WEND) GO TO 130
      GO TO 230
225   P1X=PX1(1)+PX1(2)
      P1Y=PY1(1)+PY1(2)
      PT12=P1X**2+P1Y**2
      W12=W(1)*W(2)
      AMS2=W12-PT12
C  LAST BREAK OF STRING
 230  NPTCL=I
      AMC=SQRT(AMS2)
      EC=(W(1)+W(2))/2.0
      VC(1)=P1X/EC
      VC(2)=P1Y/EC
      VC(3)=(W(1)-W(2))/(2.0*EC)
      NIN1=NPTCL+1
c.. the last break of the string will be done in clustr
      CALL CLUSTR(IFL(1),IFL(2),AMC,ilead)
      IF(LRET) GO TO 100
      NFIN1=NPTCL
      CALL LORTR(VC,NIN1,NFIN1,BACK)
      NPR=NPR+1
      NPL=NPL+1
      IF(NPR.GT.MXPTCL.OR.NPL.GT.MXPTCL) GO TO 9999
c..the hadron from the left and the right side of the string
c..are copied to the final pptcl array
      IDENTL(NPL)=IDENT(NFIN1)
      PPTL(1,NPL)=PPTCL(1,NFIN1)
      PPTL(2,NPL)=PPTCL(2,NFIN1)
      PPTL(3,NPL)=PPTCL(3,NFIN1)
      PPTL(4,NPL)=PPTCL(4,NFIN1)
      PPTL(5,NPL)=PPTCL(5,NFIN1)
      IDENTR(NPR)=IDENT(NIN1)
      PPTR(1,NPR)=PPTCL(1,NIN1)
      PPTR(2,NPR)=PPTCL(2,NIN1)
      PPTR(3,NPR)=PPTCL(3,NIN1)
      PPTR(4,NPR)=PPTCL(4,NIN1)
      PPTR(5,NPR)=PPTCL(5,NIN1)
      JJ=NFIX
      DO 284 J=1,NPR
      JJ=JJ+1
      IDENT(JJ)=IDENTR(J)
      PPTCL(1,JJ)=PPTR(1,J)
      PPTCL(2,JJ)=PPTR(2,J)
      PPTCL(3,JJ)=PPTR(3,J)
      PPTCL(4,JJ)=PPTR(4,J)
      PPTCL(5,JJ)=PPTR(5,J)
284   CONTINUE
      JJ=NFIX+NPR
      DO 285 J=1,NPL
      JJ=JJ+1
      K=NPL-J+1
      IDENT(JJ)=IDENTL(K)
      PPTCL(1,JJ)=PPTL(1,K)
      PPTCL(2,JJ)=PPTL(2,K)
      PPTCL(3,JJ)=PPTL(3,K)
      PPTCL(4,JJ)=PPTL(4,K)
      PPTCL(5,JJ)=PPTL(5,K)
285   CONTINUE
      N1=NFIX+1
      N2=NFIX+NPR+NPL-1
c.. we choose the LUND scheme
      consti=.false.

      IF(CONSTI) THEN
C------------------------------------------------------C
C----- CONSTITUENT      TIME           ----------------C
C------------------------------------------------------C
      DO 1286 J=N1,N2
      P3S=0.
      ES=0.
      DO 1287 L=N1,J
      P3S=P3S+PPTCL(3,L)
1287  ES=ES+PPTCL(4,L)
c.. TI is the formation time of the particle
c.. ZI is the z coordinate
      TI=(AMSTR-2.*P3S)/(2.*XAP)
      ZI=(AMSTR-2.*ES)/(2.*XAP)
      IF(J.NE.N2) GO TO 1288
      TII=TI
      ZII=ZI
1288  PPTCL(6,J)=0.
      PPTCL(7,J)=0.
      PPTCL(8,J)=ZI
      PPTCL(9,J)=TI
1286  CONTINUE
C
      PPTCL(6,N2+1)=0.
      PPTCL(7,N2+1)=0.
      PPTCL(8,N2+1)=ZII
      PPTCL(9,N2+1)=TII
C
      GO TO 1253
      ENDIF
C------------------------------------------------------C
C-----  INSIDE-OUTSIDE  TIME (LUND)    ----------------C
C------------------------------------------------------C
      DO 286 J=N1,NPTCL
      P3S=0.
      ES=0.
      NJ=J-1
      IF(NJ.EQ.0) GO TO 289
      DO 287 L=N1,NJ
      P3S=P3S+PPTCL(3,L)
 287  ES=ES+PPTCL(4,L)
c.. TI is the formation time of the particle
c.. ZI is the z coordinate
 289  TI=(AMSTR-2.*P3S+PPTCL(4,J)-PPTCL(3,J))/(2.*XAP)
      ZI=(AMSTR-2.*ES-PPTCL(4,J)+PPTCL(3,J))/(2.*XAP)
      PPTCL(6,J)=0.
      PPTCL(7,J)=0.
      PPTCL(8,J)=ZI
      PPTCL(9,J)=TI
 286  CONTINUE
1253  RETURN

c.. warning if to many hadrons are produced in string
c.. increase the particle arrays to avoid this
9999  WRITE(6,9998) I
9998  FORMAT(//10X,40H...stop IN STRING..NPTCL TOO HIGH NPTCL=,I5)
      stop 137
      END


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE CLUSTR(IFL1,IFL2,AMCTR,ilead)
c
cinput   amctr     : stringmass, 
cinput   ifl1      : leading quark (or diquark) along $+Z$ axis 
cinput   ifl2      : 2nd leading quark (or diquark)
cinput   ilead     : $2-ilead=$ number of leading (di-)quarks
c
c output : produced particles via common block ({\tt pptcl})  
c
C  HADRONS PRODUCTION BY MEANS CLUSTER BREAKING
C  WITH QUARK AND ANTIQUARK OR QUARK AND DIQUARK OR DIQUARK AND
C  ANTIDIQUARK IFL1 AND IFL2 ON ENDS.
c  Only the final 2 particles are created in {\tt clustr}! 
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      COMMON/COMTRY/ NTRIES
      
      PARAMETER(MXPTCL=200)
      COMMON/PARTCL/ PPTCL(9,MXPTCL),nptcl,IDENT(MXPTCL),IDCAY(MXPTCL)

      include 'comstr.f'

      COMMON/COLRET/ LRET
      LOGICAL LRET
      DIMENSION IFL(2),U(3)
      LOGICAL SPINT
      real*8 valint(1)
      common /values/ valint

      include 'options.f'

c.. strange and charm suppression in clustr (see string)
      prbs=ctparam(6)
      prbc=ctparam(7)
c..the diquark-suppresion parameter is reduced for small
c    string masses (finite size effect) see A. Jahns diploma thesis
      if(amctr.le.2.8d0)then
        pbarc=0.d0
      elseif (amctr.le.5d0)then
        pbarc=((amctr-2.8d0)/2.2d0)**3*ctparam(8)
      else
        pbarc=ctparam(8)
      endif

      dmas=ctparam(9)
      dmass=ctparam(10)

C
C  PRBS STRANGENESS SUPPRESSION PARAMETER
C  PRBC CHARM SUPPRESSION PARAMETER
      PU=1./(2.+PRBS+PRBC)
C
      NFIX=NPTCL
      NREP=0
      LRET=.FALSE.
 100  I=NFIX
      IF(NREP.LT.NTRIES) GO TO 101
      LRET=.TRUE.
      RETURN
101   CONTINUE
      KSPIN=0
      IFL(1)=IFL1
      IFL(2)=IFL2
      SPINT=.FALSE.
      I=I+2
      IF(I.GT.MXPTCL) GO TO 9999
C  CHOOSE SIDE OF BREAK
      JSIDE=1
C  IF ANY IFL IS A DIQUARK
      IF(MOD(IFL(1),100).EQ.0.OR.MOD(IFL(2),100).EQ.0) GO TO 150
C  IFL(1) AND IFL(2) ARE QUARKS
C  SELECT Q,QBARPAIR OR QQ,QQBAR PAIR
      DRND=ranf(0)
      IF(DRND.LT.PBARC.and.valint(1).eq.0.d0) GO TO 140
C  Q,QBAR PAIR
      IFLN=ISIGN(IFLAV(PU,PRBS),-IFL(JSIDE))
      GO TO 200
C  QQ,QQBAR PAIR
140   IQ1=IFLAV(PU,PRBS)
      IQ2=IFLAV(PU,PRBS)
      IF(IQ1.LE.IQ2) GO TO 145
      ISWAP=IQ1
      IQ1=IQ2
      IQ2=ISWAP
145   IFQQ=1000*IQ1+100*IQ2
      IFLN=ISIGN(IFQQ,IFL(JSIDE))
      GO TO 200
C  IFL(1) OR IFL(2) IS DIQUARK
C  Q,QBAR PAIR
150   IPSIGN=IFL(JSIDE)
      IF(MOD(IFL(JSIDE),100).EQ.0) GO TO 130
      IPSIGN=-IFL(JSIDE)
130   IFLN=ISIGN(IFLAV(PU,PRBS),IPSIGN)
C  IDENTS AND MASSES OF PARTICLES
 200  continue 
c..quark-flip included (to describe some phi data) 
      if(CTParam(5).gt.ranf(0).and.mod(ifln,100).ne.0.and.
     &   mod(ifl(jside),100).ne.0.and.mod(ifl(3-jside),100).ne.0)then 
c quark-flip
        IDENT(I-1)=IDPARC(IFL(JSIDE),IFL(3-JSIDE),SPINT,KSPIN) 
        IDENT(I)=IDPARC(-IFLN,IFLN,SPINT,KSPIN)
      else 
        IDENT(I-1)=IDPARC(IFL(JSIDE),IFLN,SPINT,KSPIN)
        IDENT(I)=IDPARC(IFL(3-JSIDE),-IFLN,SPINT,KSPIN)
      end if
c..for special bbar-b annihilation reactions for conservation 
c  of total quantum numbers 
      if(valint(1).ne.0.d0)then
        ifq1=int(valint(1)/10.)
        ifq2=-mod(int(valint(1)),10)
        if(isign(1,ifln).eq.isign(1,ifq1))then
          IDENT(I-1)=IDPARC(IFL(JSIDE),ifq1,SPINT,KSPIN)
          IDENT(I)=IDPARC(IFL(3-JSIDE),ifq2,SPINT,KSPIN)
        else
          IDENT(I-1)=IDPARC(IFL(JSIDE),ifq2,SPINT,KSPIN)
          IDENT(I)=IDPARC(IFL(3-JSIDE),ifq1,SPINT,KSPIN)
        endif
      endif

      PPTCL(5,I-1)=AMASS(IDENT(I-1))
      PPTCL(5,I)=AMASS(IDENT(I))
C  IF TOO LOW MASS,START ALL OVER (i.e. goto 100)
      DEMAS=0.15
      IF(IFLN.LT.3) DEMAS=0.
c      IF(AMCTR.GT.PPTCL(5,I-1)+PPTCL(5,I)+DEMAS)  GO TO 102
      IF(AMCTR.GT.PPTCL(5,I-1)+PPTCL(5,I)+DEMAS) then
       if(mod(ifl1,100).eq.0.or.mod(ifl2,100).eq.0)goto 102
c..maximum kinetic energy cutoff for meson-clustr:
c..we want a lot of energy in the particle mass in this last break
       IF(AMCTR-PPTCL(5,I-1)-PPTCL(5,I).lt.ctparam(43)) GO TO 102
      endif
      NREP=NREP+1
c.. 100 starts all over
      GO TO 100

102   CONTINUE

c.. isotropic px py pz distribution
      PA=DBLPCM(AMCTR,PPTCL(5,I-1),PPTCL(5,I))
      U(3)=1.-2.*ranf(0)
      PHI=2.*PI*ranf(0)
      ST=SQRT(1.-U(3)**2)
      U(1)=ST*COS(PHI)
      U(2)=ST*SIN(PHI)
      PPTCL(1,I-1)=PA*U(1)
      PPTCL(1,I)=-(PA*U(1))
      PPTCL(2,I-1)=PA*U(2)
      PPTCL(2,I)=-(PA*U(2))
      PPTCL(3,I-1)=PA*U(3)
      PPTCL(3,I)=-(PA*U(3))
      PA2=PA**2
      PPTCL(4,I-1)=SQRT(PA2+PPTCL(5,I-1)**2)
      PPTCL(4,I)=SQRT(PA2+PPTCL(5,I)**2)
      IDCAY(I-1)=0
      IDCAY(I)=0
      NPTCL=I
      
c..forward/backward distribution in clustr for baryons 
c..(no pt in the last string break!)
c..pt for the baryon comes from parton kick in the excitation
c      if(abs(ident(i)).ge.1000.or.abs(ident(i-1)).ge.1000)then
c      PPTCL(1,I-1)=0.d0
c      PPTCL(1,I)=0.d0
c      PPTCL(2,I-1)=0.d0
c      PPTCL(2,I)=0.d0
c      PPTCL(3,I-1)=PA
c      PPTCL(3,I)=-PA
c      PA2=PA**2     
c      PPTCL(4,I-1)=SQRT(PA2+PPTCL(5,I-1)**2)
c      PPTCL(4,I)=SQRT(PA2+PPTCL(5,I)**2)
c      IDCAY(I-1)=0
c      IDCAY(I)=0    
c      NPTCL=I       
c      endif

c..if baryon number=+-1, force the (anti-)baryon in positive
c..z-direction (just pick the right hemisphere)
      if(ctoption(29).gt.0)then
      if( (iabs(ident(i)/1000).ne.0.and.iabs(ident(i-1)/1000).eq.0
     &    .and.pptcl(3,i).lt.0.d0).or.
     &  (iabs(ident(i-1)/1000).ne.0.and.iabs(ident(i)/1000).eq.0
     &    .and.pptcl(3,i-1).lt.0.d0)) then
            pptcl(3,i)  =-pptcl(3,i)
            pptcl(3,i-1)=-pptcl(3,i-1)
      endif
      endif
        

      RETURN
c.. particle array to small warning:
9999  WRITE(6,9998) I
9998  FORMAT(//10X,40H...stop IN CLUSTR..NPTCL TOO HIGH NPTCL=,I5)
      stop 137
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer FUNCTION IFLAV(PU,PRBS)
c
cinput PU      : 1-{\tt PU}= up (down, resp.) probability
cinput PRBS    : Strange quark suppression
c
c output : {\tt iflav}:  flavor of created quark
c
c     Returns quark flavor acc. to suppression prob's:
c     1=up, 2=down, 3=strange, 4=charm
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

C
      RNDOM=ranf(0)
C
      IF(RNDOM.GT.PU) GO TO 1
c..create up quark
      IFLAV=1
      RETURN
  1   IF(RNDOM.GT.2.0*PU) GO TO 2
c..create down quark
      IFLAV=2
      RETURN
  2   IF(RNDOM.GT.PU*(2.0+PRBS)) GO TO 3
c..create strange quark
      IFLAV=3
      RETURN
c..create charm quark
  3   IFLAV=4
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 FUNCTION ZFRAGS(IFL,IFLN,PT2,ZMIN,ZMAX,leading)
c
cinput IFL       : ID of existing quark 
cinput IFLN      : ID of newly created quark
cinput PT2       : $p_t$ of newly created hadron 
cinput ZMIN      : lowest allowed longitudinal momentum fraction
cinput ZMAX      : highest allowed longitudinal momentum fraction
cinput leading   : flag for leading particle 
c
coutput ZFRAGS   : longitudinal momentum fraction of created hadron
c
c     According to the fragmentation function(s), longitudinal momentum 
c     is assigned to the hadron.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      LOGICAL leading

      include 'options.f'

      COMMON/INPRNT/ ITDKY,ITLIS
      PARAMETER(ALFT=0.5,ARHO=0.5,APHI=0.,APSI=-2.)
      PARAMETER(AN=-0.5,ALA=-0.75,ALAC=-1.75)
      PARAMETER(AKSI=-1.0,AUSC=-2.0,AUCC=-2.0)

c.. cto 21 chooses the fragmentation fct.
      if ((.not.leading).or.(ifln.eq.3)) then
         affm=ctparam(47)
         bffm=ctparam(48)
 5108  ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
       YF=((1d0-ZFRAGS)**bffm*(bffm+1)*affm+1-affm)/3d0
       if(yf.gt.1d0.or.yf.lt.0d0)then 
        write(6,*)'ZFRAGS: wrong norm:',yf,zmin,zmax,zfrags
       end if
       IF(ranf(0).LE.YF) RETURN
       GO TO 5108
      endif

c.. GAUSSIAN fragmentation fct.
      if(CTOption(21).eq.0)then  
         affm=CTParam(36)
         bffm=CTParam(37)
c.. suppress low momentum particles
      deltaz=zmax-zmin
      zmin=zmin+deltaz*0.25d0
 108  ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
c.. gaussian distr.
      yf=exp(-(((zfrags-bffm)**2)/(2.*affm**2)))
c.. field-feynmann fragmentation fct.
c      YF=((1d0-ZFRAGS)**bffm*(bffm+1)*affm+1-affm)/3d0     
      if(yf.gt.1d0.or.yf.lt.0d0)then
       write(6,*)'ZFRAGS: wrong norm:',yf,zmin,zmax,zfrags
      end if
c      return
      IF(ranf(0).LE.YF) RETURN                    
      GO TO 108                                   

      else if(CTOption(21).eq.1)then
c..lund-fragmentation fct.
 1008 ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1-zfrags)*exp(-(0.7*pt2/zfrags))/3d0
      if(yf.gt.1d0.or.yf.lt.0d0)then
       write(6,*)'ZFRAGS: wrong norm:',yf,zmin,zmax,zfrags
      end if
      IF(ranf(0).LE.YF) RETURN                    
      GO TO 1008 

      else if(CTOption(21).eq.2)then 
c.. kaidalov's fragmentation fct. 
      ID1=IABS(IFL)
      ID2=IABS(IFLN)
      IF(MOD(ID2,100).EQ.0) GO TO 15
      GO TO(1,2,3,4),ID2
C  UU-TRAJECTORY
 1    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.0-ZFRAGS)**(ALFT-ARHO)
      IF(ranf(0).LE.YF) RETURN
      GO TO  1
C DD-TRAJECTORY
 2    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-ARHO)
      IF(ranf(0).LE.YF) RETURN
      GO TO  2
C SS-TRAJECTORY
 3    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-APHI)
      IF(ranf(0).LE.YF) RETURN
      GO TO  3
C CC-TRAJECTORY
 4     ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-APSI)
      IF(ranf(0).LE.YF) RETURN
      GO TO  4
C
15    Continue
      CALL FLAVOR(ID2,IFL2,IFL3,IFL1,ISPIN)
      IF((IFL2.EQ.1.AND.IFL3.EQ.1)) GO TO 16
      IF((IFL2.EQ.1.AND.IFL3.EQ.2)) GO TO 17
      IF((IFL2.EQ.1.AND.IFL3.EQ.3)) GO TO 18
      IF((IFL2.EQ.1.AND.IFL3.EQ.4)) GO TO 19
      IF((IFL2.EQ.2.AND.IFL3.EQ.2)) GO TO 20
      IF((IFL2.EQ.2.AND.IFL3.EQ.3)) GO TO 21
      IF((IFL2.EQ.2.AND.IFL3.EQ.4)) GO TO 22
      IF((IFL2.EQ.3.AND.IFL3.EQ.3)) GO TO 23
      IF((IFL2.EQ.3.AND.IFL3.EQ.4)) GO TO 24
      IF((IFL2.EQ.4.AND.IFL3.EQ.4)) GO TO 25
C UUUU-TRAJECTORY
16    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*AN-ARHO))
      IF(ranf(0).LE.YF) RETURN
      GO TO 16
C UDUD-TRAJECTORY
17    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*AN-ARHO))
      IF(ranf(0).LE.YF) RETURN
      GO TO 17
C USUS-TRAJECTORY
18    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*ALA-ARHO))
      IF(ranf(0).LE.YF) RETURN
      GO TO 18
C UCUC-TRAJECTORY
19    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*ALAC-ARHO))
      IF(ranf(0).LE.YF) RETURN
      GO TO 19
C DDDD-TRAJECTORY
20    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*AN-ARHO))
      IF(ranf(0).LE.YF) RETURN
      GO TO 16
C DSDS-TRAJECTORY
21    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*ALA-ARHO))
      IF(ranf(0).LE.YF) RETURN
      GO TO 21
C DCDC-TRAJECTORY
22    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*ALAC-ARHO))
      IF(ranf(0).LE.YF) RETURN
      GO TO 22
C SSSS-TRAJECTORY
23    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*AKSI-ARHO))
      IF(ranf(0).LE.YF) RETURN
      GO TO 23
C SCSC-TRAJECTORY
24    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*AUSC-ARHO))
      IF(ranf(0).LE.YF) RETURN
      GO TO 24
C CCCC-BARYON
25    ZFRAGS=ZMIN+ranf(0)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*AUCC-ARHO))
      IF(ranf(0).LE.YF) RETURN
      GO TO 25
      else
       write(6,*)'string.f: cto(21)=',ctoption(21),' not valid'
       stop 137
      end if

      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE LORTR(V,NIN,NFIN,BACK)
c
cinput V         : boost velocity (3-vector)
cinput NIN       : lower boundary in the pptcl array
cinput NFIN      : upper boundary in the pptcl array
cinput back      : inversion flag for transformation 
c
c output : via common block ({\tt pptcl})
c
c     Performs a Lorentz-boost of a part of the pptcl array
c     (from {\tt nin} to {\tt nfin})
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      PARAMETER(MXPTCL=200)
      COMMON/PARTCL/ PPTCL(9,MXPTCL),nptcl,IDENT(MXPTCL),IDCAY(MXPTCL)
      DIMENSION V(3),VD(3)
      LOGICAL BACK
      L=1
      IF(BACK) L=-1
      DO 3 I=1,3
3     VD(I)=V(I)
      VVD=VD(1)*VD(1)+VD(2)*VD(2)+VD(3)*VD(3)
      GAD=1.D0/DSQRT(DABS(1.D0-VVD))
      GA=GAD
      DO 100 J=NIN,NFIN
      VP=V(1)*PPTCL(1,J)+V(2)*PPTCL(2,J)+V(3)*PPTCL(3,J)
      GAVP=GA*(GA*VP/(1.+GA)-FLOAT(L)*PPTCL(4,J))
      PPTCL(1,J)=PPTCL(1,J)+GAVP*V(1)
      PPTCL(2,J)=PPTCL(2,J)+GAVP*V(2)
      PPTCL(3,J)=PPTCL(3,J)+GAVP*V(3)
      PMAS=PPTCL(5,J)
      PPTCL(4,J)=SQRT(PPTCL(1,J)**2+PPTCL(2,J)**2+PPTCL(3,J)**2+
     +PMAS**2)
100   CONTINUE
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      BLOCK DATA D2
c
c     Initial values for several common blocks
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
C
C    INPUT QUARK DISTRIBUTION PARAMETERS
C

      include 'comstr.f'

      COMMON/KAPPA/ XAP
      COMMON/CONSTI/ CONSTI
      LOGICAL CONSTI
       DATA CONSTI/.FALSE./
C
C
C   DATA FOR COSPAR 
c this parameter controls the min. energy for ending string fragmentation
c default      DATA DMAS/0.35/
c now: ctparam(9) and ctparam(10)
c      DATA DMAS/1./
c      data dmass/0.25/
C
C
C  INPUT STRING TENSION (FERMI/GEV)
      DATA XAP/1.0/
C
C  INPUT FRAGMENTATION PARAMETERS:
C
      DATA PARQLS/0.5/,PARRS/0.9/
C  PRBS STRANGENESS SUPPRESSION PARAMETER
C  PRBC CHARM SUPPRESSION PARAMETER
C  PBARC AND PBARS BARYON PAIRS SUPPRESION PARAMETERS
C      DATA PRBS/0.3/,PRBC/0.0005/,PBARC/0.09/,PBARS/0.09/
c these are now ctparam(6),ctparam(7) and ctparam(8)
c      DATA PRBS/0.33/,PRBC/0.0/,PBARC/0.09/,PBARS/0.09/
C
C  THESE PARAMETERS ARE IMPORTANT FOR RESONANCE PRODUCTION
      DATA PJSPNC/.50/
      DATA PJSPNS/.50/
c      DATA PJSPNC/1./
c      DATA PJSPNS/1./
C      DATA PMIX1C/.25,.25,.5,.5,.5,1./,PMIX2C/.5,.5,1.,0.,0.,1./
c      DATA PMIX1S/.25,.25,.5,.5,.5,1./,PMIX2S/.5,.5,1.,0.,0.,1./
c QUARK MIXING PARAMETERS
cspl-0795 these parameters assign the pure u/ubar,d/dbar,s/sbar states
c  (e.g. 110,220,330) to the physical particles according to the su(3)
c  quark model. The flavor mixing angles are chosen according to quadratic
c  Gell-Mann-Okubo mass formula (Review of Particle Properties,
c  Phys. Rev D50 (1994) 1319). For the scalar mesons this formula is not 
c  applicable. We assume an ideal mixing angle (tan(theta)=1/sqrt(2)).
c
c pseudoscalar: theta=-10 deg
c vector      : theta= 39 deg
c pseudovector: theta= 51 deg
c tensor      : theta= 28 deg
c
c PMIX1C/S(IF1,JSPIN) gives the quark content of the heavy isosinglet:
c  e.g. the eta' (330) is 25% u/ubar 25% d/dbar and 50% s/sbar
c PMIX2C/S(IF1,JSPIN) chooses for nonstrange qqbars between the isosinglet 
c  the and the isotriplet 
c  quark content of the heavy singlett
c      DATA PMIX1C/.25,.25,.5,   ! eta'
c     &     0.,0.,1.,    ! phi
c     &     0.,0.,1.,    ! f_0
c     &     .04,.04,.92, ! f_1'
c     &     .01,.01,.98/ ! f_2'
c      DATA PMIX2C/.5 ,.5 ,1.,.5,.5,1.,.5,.5,1.,.5 ,.5 ,1. ,.5 ,.5 , 1./
c      DATA PMIX1S/.25,.25,.5,0.,0.,1.,0.,0.,1.,.04,.04,.92,.01,.01,.98/,
c     &     PMIX2S/.5 ,.5 ,1.,.5,.5,1.,.5,.5,1.,.5 ,.5 ,1. ,.5 ,.5 , 1./
C
C  SIGQTS IS IMPORTANT PARAMETER FOR PRODUCED HADRON TRANSVERSE MOMENTA
C
c       DATA SIGQTS/0.65/ -> ctparam(42)
       END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      BLOCK DATA D4
c
cinput ()      : NONE
c
coutput ()     : via common blocks
c
c     Initial values for several common blocks
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      
      COMMON/COMTRY/ NTRIES
C
C
C          DATA COMTRY
      DATA NTRIES/1000/
C

      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer FUNCTION IDPARS(IFL01,IFL02,SPINT,IR)
c
cinput IFL01  : ID of (anti-)quarks/di-quarks 
cinput IFL02  : ID of (anti-)quarks/di-quarks 
cinput SPINT        : flag for spin assignment
cinput IR           : Determines particle spin
c
c output : Quark code of the hadron
c
C   CONSTRUCT MESON FROM QUARK AND ANTIQUARK WITH FLAVORS IFL01,IFL02
C   OR CONSTRUCT BARYON FROM DIQUARK AND QUARK OR ANTIDIQUARK AND
C   ANTIQUARK WITH FLAVORS IFL01,IFL02.
c   THE MESON MULTIPLETT IS CHOSEN ACC. TO SUPPRESSION PARAM'S:
c parm gives the probability for different meson multiplets according
c to spin degeneracy and average mass ratios
c spin-parity 0- : 1- : 0+ : 1+ : 2+ = parm(1):parm(2)...:parm(5)
c If SPINT=.t., IR will be used to assign particle spin
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)


      LOGICAL SPINT
      include 'options.f'

      include 'comstr.f'

C
      IFL1=IFL01
      IFL2=IFL02
C  CONSTRUCT MESON WITH ACCOUNT FLAVOR MIXING
      IF(MOD(IFL1,100).EQ.0) GO TO 420
      IF(MOD(IFL2,100).EQ.0) GO TO 425
c
      call getbran(parm,1,njspin,psdum,1,njspin,jspin)
      jspin=jspin-1

      IF(SPINT.AND.IR.EQ.2) JSPIN=0
      IF(SPINT.AND.IR.EQ.1) JSPIN=1
      ID1=IFL1
      ID2=IFL2
      IF(ID1+ID2.NE.0) GO TO 400
C
      ID1=IABS(ID1)
      IF(ID1.GE.4) GO TO 401
      RND=ranf(0)
      ID1=INT(PMIX1S(ID1,JSPIN+1)+RND)+INT(PMIX2S(ID1,JSPIN+1)+RND)+1
 401  ID2=-ID1
C
 400  IF(IABS(ID1).LE.IABS(ID2)) GO TO 410
      ISAVE=ID1
      ID1=ID2
      ID2=ISAVE
 410  IDHAD=ISIGN(100*IABS(ID1)+10*IABS(ID2)+JSPIN,ID1)
      GO TO 470
C  CONSTRUCT BARYON IDENT
 420  ID3=ISIGN(MOD(IFL1/100,10),IFL1)
      ID2=IFL1/1000
      ID1=IFL2
      GO TO 430
 425  ID3=ISIGN(MOD(IFL2/100,10),IFL2)
      ID2=IFL2/1000
      ID1=IFL1
 430  IF(IABS(ID1).LE.IABS(ID2)) GO TO 431
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
 431  IF(IABS(ID2).LE.IABS(ID3)) GO TO 432
      ISWAP=ID2
      ID2=ID3
      ID3=ISWAP
 432  IF(IABS(ID1).LE.IABS(ID2)) GO TO 440
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
 440  JSPIN=1
      IF(ID1.EQ.ID2.AND.ID2.EQ.ID3) GO TO 450
      JSPIN=INT(ranf(0)+PJSPNS)
      IF(SPINT.AND.IR.EQ.2) JSPIN=0
      IF(SPINT.AND.IR.EQ.1) JSPIN=1

 450  IF(JSPIN.EQ.1.OR.ID1.EQ.ID2.OR.ID2.EQ.ID3) GO TO 460
      DRND=ranf(0)
      IF(DRND.GT.PJSPNS) GO TO 460
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
 460  IDHAD=1000*IABS(ID1)+100*IABS(ID2)+10*IABS(ID3)+JSPIN
      IDHAD=ISIGN(IDHAD,IFL1)
 470  IDPARS=IDHAD
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer FUNCTION IDPARC(IFL01,IFL02,SPINT,IR)
c
cinput IFL01  : ID of (anti-)quarks/di-quarks 
cinput IFL02  : ID of (anti-)quarks/di-quarks 
cinput SPINT        : flag for spin assignment
cinput IR           : Determines particle spin
c
c output : Quark code of the hadron
c
C   CONSTRUCT MESON FROM QUARK AND ANTIQUARK WITH FLAVORS IFL01,IFL02
C   OR CONSTRUCT BARYON FROM DIQUARK AND QUARK OR ANTIDIQUARK AND
C   ANTIQUARK WITH FLAVORS IFL01,IFL02.
c   THE MESON MULTIPLETT IS CHOSEN ACC. TO SUPPRESSION PARAM'S:
c parm gives the probability for different meson multiplets according
c to spin degeneracy and average mass ratios
c spin-parity 0- : 1- : 0+ : 1+ : 2+ = parm(1):parm(2)...:parm(5)
c If SPINT=.t., IR will be used to assign particle spin
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      LOGICAL SPINT
      include 'options.f'

      include 'comstr.f'

C
      IFL1=IFL01
      IFL2=IFL02
C  CONSTRUCT MESON WITH ACCOUNT FLAVOR MIXING
      IF(MOD(IFL1,100).EQ.0) GO TO 420
      IF(MOD(IFL2,100).EQ.0) GO TO 425

c.. choose multiplett by its probability
      call getbran(parm,1,njspin,dummy,1,njspin,jspin)
      jspin=jspin-1

      IF(SPINT.AND.IR.EQ.2) JSPIN=0
      IF(SPINT.AND.IR.EQ.1) JSPIN=1
      ID1=IFL1
      ID2=IFL2
      IF(ID1+ID2.NE.0) GO TO 400
C
      ID1=IABS(ID1)
      IF(ID1.GE.4) GO TO 401
c.. singlet mixing acc. to mixing angles
      RND=ranf(0)
      ID1=INT(PMIX1C(ID1,JSPIN+1)+RND)+INT(PMIX2C(ID1,JSPIN+1)+RND)+1
 401  ID2=-ID1
C
 400  IF(IABS(ID1).LE.IABS(ID2)) GO TO 410
      ISAVE=ID1
      ID1=ID2
      ID2=ISAVE
 410  IDHAD=ISIGN(100*IABS(ID1)+10*IABS(ID2)+JSPIN,ID1)
      GO TO 470
C  CONSTRUCT BARYON IDENT
 420  ID3=ISIGN(MOD(IFL1/100,10),IFL1)
      ID2=IFL1/1000
      ID1=IFL2
      GO TO 430
 425  ID3=ISIGN(MOD(IFL2/100,10),IFL2)
      ID2=IFL2/1000
      ID1=IFL1
 430  IF(IABS(ID1).LE.IABS(ID2)) GO TO 431
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
 431  IF(IABS(ID2).LE.IABS(ID3)) GO TO 432
      ISWAP=ID2
      ID2=ID3
      ID3=ISWAP
 432  IF(IABS(ID1).LE.IABS(ID2)) GO TO 440
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
 440  JSPIN=1
      IF(ID1.EQ.ID2.AND.ID2.EQ.ID3) GO TO 450
      JSPIN=INT(ranf(0)+PJSPNC)
      IF(SPINT.AND.IR.EQ.2) JSPIN=0
      IF(SPINT.AND.IR.EQ.1) JSPIN=1
 450  IF(JSPIN.EQ.1.OR.ID1.EQ.ID2.OR.ID2.EQ.ID3) GO TO 460
      DRND=ranf(0)
      IF(DRND.GT.PJSPNC) GO TO 460
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
 460  IDHAD=1000*IABS(ID1)+100*IABS(ID2)+10*IABS(ID3)+JSPIN
      IDHAD=ISIGN(IDHAD,IFL1)
 470  IDPARC=IDHAD
       RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 FUNCTION AMASS(ID)
c
c
cinput ID      : Quark code
c
c
C          THIS FUNCTION RETURNS THE MASS OF THE PARTICLE WITH
C          IDENT CODE ID. (QUARK-BASED IDENT CODE)
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      include 'comres.f'
      include 'comstr.f' 

      real*8 mmin,mmax,m0,mminit, massit, widit
      integer isoit
      dimension amq(4)
      include 'options.f'

c.. quark masses (u,d,s,c)
      DATA AMq/.15,.15,.45,1.6/

c fraction of baryonresonances
      bresfrac=ctparam(11)

      CALL FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN)
      idabs=iabs(id)

c get quark masses
      if (idabs.le.4) then
        amass=amq(idabs)
        return
      endif
c get diquark masses 
      IF(ID.NE.0.AND.MOD(ID,100).EQ.0) then
        AMASS=AMq(IABS(IFL1))+AMq(IABS(IFL2))
        return
      endif

c get hadron masses

c get baryon masses
c (anti-)nucleon ? 
      if ((idabs.eq.1120).or.(idabs.eq.1220)) then
        if (ranf(0).lt.bresfrac)then
c.. N*
          amass=getmass(mresmax,1)
          return
        else
c.. N
          amass=0.938
          return
        endif
      endif
c (anti-)delta ?
      if ((idabs.eq.1111).or.(idabs.eq.1121)
     &   .or.(idabs.eq.1221).or.(idabs.eq.2221)) then
        if (ranf(0).lt.bresfrac)then
c.. D*
          amass=getmass(mresmax,2)
          return
        else
c.. Delta(1232)
          m0=massit(mindel)
          w0=widit(mindel)        
c get meson mass accord. to breit wigner distr. 
          mmin=mminit(mindel)
          mmax=m0+3d0*w0  
          call getmas(m0,w0,mindel,isoit(mindel),mmin,mmax,-1d0,amass)
          return 
        endif
      endif
c (other baryons)
      if(idabs-1000.ge.0) then 
        call id2ityp(id,0.d0,itypin,iz2)
c..check range and avoid double counting for explicitely treated resonances
        if ((abs(itypin).ge.minlam.and.abs(itypin).le.maxcas).and.
     &   (idabs.ne.1231.and.idabs.ne.1131.and.idabs.ne.2231.and.
     &    idabs.ne.1331.and.idabs.ne.2331).and.
     &   (ctoption(31).eq.1)) then
c.. generate also non-groundstate ityps for strange baryons
         call probitypres(itypin,ityp)
        else
           ityp=itypin
        endif
         amass=massit(iabs(ityp))
       return
      endif

c mesons
      if(idabs.le.330+njspin-1) then
        call id2ityp(id,0.d0,ityp,iz2)
        m0=massit(iabs(ityp))
        w0=widit(iabs(ityp))        
c get meson mass accord. to breit wigner distr. 
        mmin=max(mminit(iabs(ityp)),m0-3d0*w0)
        mmax=m0+3d0*w0  
        call getmas(m0,w0,ityp,iz2,mmin,mmax,-1d0,amass)
        return
      endif
      if(idabs.ge.440.and.idabs.le.447)then
        amass=massit(135)
        return
      else if(idabs.ge.340.and.idabs.le.347)then
        amass=massit(138)
        return
      endif
      write(6,*)'! amass: no mass for part.id:',id
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 FUNCTION DBLPCM(A,B,C)
c
cinput A      : Mass of particle A
cinput B      : Mass of particle B
cinput C      : Mass of particle C
c
c
c In the rest frame of a particle of mass A, decaying into particles
c of masses B and C, {\tt dblpcm} returns the momenta of the outgoing 
c particles.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      DA=A
      DB=B
      DC=C
      DVAL=(DA**2-DB**2-DC**2)**2-(2.D0*DB*DC)**2
      DBLPCM=0.
      IF(DVAL.GT.0.D0)DBLPCM=DSQRT(DVAL)/(2.D0*DA)
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ityp2id(ityp,iz2,ifa,ifb)
c
cinput ityp    : UrQMD particle ID 
cinput iz2     : UrQMD $2\cdot I_3$ 
c
coutput ifa : quarkcode (diquark)
coutput ifb : quarkcode (quark)
c
c  returns quark id from uqmd-ityp and isospin z-component (times 2) 
c
c  quark ids are:  1 up,  2 down,  3 strange,  4 charm
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'comres.f'
      integer itypabs,ityp,iz2,ifa,ifb,sf
      integer t3,if(3)
      itypabs=iabs(ityp)
      t3=iz2*isign(1,ityp)
            
c nucleons ?
      if (itypabs.lt.minmes) then
       sf=strres(itypabs)   
c uuu
       if(t3.eq.3) then
        if(1)=1
        if(2)=1
        if(3)=1
        call mquarks(if,ifa,ifb)
        ifa=ifa*isign(1,ityp)
        ifb=ifb*isign(1,ityp)
        return
       endif
c uus
       if(t3.eq.2) then
        if(1)=1
        if(2)=1
        if(3)=3
        call mquarks(if,ifa,ifb)
        ifa=ifa*isign(1,ityp)
        ifb=ifb*isign(1,ityp)
        return
       endif
c uud
       if((t3.eq.1).and.(sf.eq.0)) then
        if(1)=1
        if(2)=1
        if(3)=2
        call mquarks(if,ifa,ifb)
        ifa=ifa*isign(1,ityp)
        ifb=ifb*isign(1,ityp)
        return
       endif
c uss
       if((t3.eq.1).and.(sf.eq.2)) then
        if(1)=1
        if(2)=3
        if(3)=3
        call mquarks(if,ifa,ifb)
        ifa=ifa*isign(1,ityp)
        ifb=ifb*isign(1,ityp)
        return
       endif
c uds
       if((t3.eq.0).and.(sf.eq.1)) then
        if(1)=1
        if(2)=2
        if(3)=3
        call mquarks(if,ifa,ifb)
        ifa=ifa*isign(1,ityp)
        ifb=ifb*isign(1,ityp)
        return
       endif
c sss
       if((t3.eq.0).and.(sf.eq.3)) then
        if(1)=3
        if(2)=3
        if(3)=3
        call mquarks(if,ifa,ifb)
        ifa=ifa*isign(1,ityp)
        ifb=ifb*isign(1,ityp)
        return
       endif
c udd
       if((t3.eq.-1).and.(sf.eq.0)) then
        if(1)=1
        if(2)=2
        if(3)=2
        call mquarks(if,ifa,ifb)
        ifa=ifa*isign(1,ityp)
        ifb=ifb*isign(1,ityp)
        return
       endif
c dss
       if((t3.eq.-1).and.(sf.eq.2)) then
        if(1)=2
        if(2)=3
        if(3)=3
        call mquarks(if,ifa,ifb)
        ifa=ifa*isign(1,ityp)
        ifb=ifb*isign(1,ityp)
        return
       endif
c dds
       if (t3.eq.-2) then
        if(1)=2
        if(2)=2
        if(3)=3
        call mquarks(if,ifa,ifb)
        ifa=ifa*isign(1,ityp)
        ifb=ifb*isign(1,ityp)
        return
       endif
c ddd
       if (t3.eq.-3) then
        if(1)=2
        if(2)=2
        if(3)=2
        call mquarks(if,ifa,ifb)
        ifa=ifa*isign(1,ityp)
        ifb=ifb*isign(1,ityp)
        return
       endif

      endif
c-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --  
c bosons
      if (itypabs.ge.minmes) then
       sf=strmes(itypabs)   
c d ubar
       if(t3.eq.-2) then
        ifa=2
        ifb=-1
        return
       endif
c d sbar
       if((t3.eq.-1).and.(sf*isign(1,ityp).eq.-1)) then
        ifa=2
        ifb=-3
        return
       endif
c s ubar
       if((t3.eq.-1).and.(sf*isign(1,ityp).eq.+1)) then
        ifa=3
        ifb=-2
        return
       endif
c q qbar
       if(t3.eq.0) then
c u ubar
c all neutral triplet states:
        if ((itypabs.eq.minmes+1)   .or.      !pi0                           
     &       (itypabs.eq.minmes+4)  .or.                                 
     &       (itypabs.eq.minmes+11)  .or.                                 
     &       (itypabs.eq.minmes+14)  .or.                                 
     &       (itypabs.eq.minmes+18)  .or.                                 
     &       (itypabs.eq.minmes+22)  .or.                                 
     &       (itypabs.eq.minmes+26)  .or.                                 
     &       (itypabs.eq.minmes+30)  .or.                                 
c the gamma gets also a quark content
     &       (itypabs.eq.minmes)) then
         ifa=1
         ifb=-1
         return
        endif
c d dbar
c light singlet states (generally less strange quark content)
        if ((itypabs.eq.minmes+2) .or.        !eta
     &      (itypabs.eq.minmes+3)  .or.
     &      (itypabs.eq.minmes+5)  .or.    
     &      (itypabs.eq.minmes+15)  .or.  
     &      (itypabs.eq.minmes+19)  .or. 
     &      (itypabs.eq.minmes+23)  .or. 
     &      (itypabs.eq.minmes+27)  .or.
     &      (itypabs.eq.minmes+31)) then 
         ifa=2
         ifb=-2
         return
        endif
c s sbar
        if ((itypabs.eq.minmes+7) .or.        !eta' (958)
     &      (itypabs.eq.minmes+9)  .or.        !phi (1020)
     &      (itypabs.eq.minmes+12) .or.        !f_0  (980)
     &      (itypabs.eq.minmes+16) .or.        !f_1 (1510)
     &      (itypabs.eq.minmes+20) .or.        !f_2'(1525)
     &      (itypabs.eq.minmes+24) .or.        !f_2'(1525)
     &      (itypabs.eq.minmes+28) .or.        !f_2'(1525)
     &      (itypabs.eq.minmes+32))then
         ifa=3
         ifb=-3
         return
        endif
       endif
c u sbar
       if((t3.eq.1).and.(sf*isign(1,ityp).eq.-1)) then
        ifa=1
        ifb=-3
        return
       endif
c s dbar
       if((t3.eq.1).and.(sf*isign(1,ityp).eq.+1)) then
        ifa=3
        ifb=-1
        return
       endif
c u dbar
       if(t3.eq.2) then
        ifa=1
        ifb=-2
        return
       endif
c charm:
       if(itypabs.ge.133.and.itypabs.lt.135) then
c D mesons
         if(iz2.eq.1.and.(ityp.eq.-133.or.ityp.eq.-134))then
           ifa=-4
           ifb=1
           return
         else if(iz2.eq.-1.and.(ityp.eq.133.or.ityp.eq.134))then
           ifa=4
           ifb=-1
           return
         else if(iz2.eq.1.and.(ityp.eq.133.or.ityp.eq.134))then
           ifa=4
           ifb=-2
           return
         else if(iz2.eq.-1.and.(ityp.eq.-133.or.ityp.eq.-134))then
           ifa=-4
           ifb=2
           return
         endif
       elseif(itypabs.gt.135.and.itypabs.le.137) then
c charmonium
          ifa=4
          ifb=-4
          return
       else if(itypabs.ge.138.and.itypabs.le.139) then
          if(ityp.eq.138.or.ityp.eq.139) then
            ifa=4
            ifb=-3
            return
          else
            ifa=-4
            ifb=3
            return
          end if
       endif


      endif

c  in any other case we will do a mesonic string
c.. and a warning !
      write(*,*)'! ityp2id: ityp',ityp,',iz2',iz2,
     &          ' can not be converted into id. Please check.'
      ifa=1
      ifb=-1      
      RETURN
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mquarks(if,ifa,ifb)
c
cinput if     : single quarks (array)
c
coutput ifa : diquark
coutput ifb : quark
c    
c this routine adds randomly the single quarks of the
c baryon to become a diquark and a quark. 
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

      include 'options.f'

      real*8 ranf
      integer if(0:2),ifa,ifb,ir

      integer i,prod,inons
 

      ir=int(3.*ranf(0))

      ifa=1000*if(ir)+100*if(mod(ir+1,3))
      ifb=if(mod(ir+2,3))
c.. check if heavy quark clusters are switched off (cto 37=0)
      if (CTOption(37).eq.0)  return

c.. not switched off -> clusters strange quarks to
c.. diquark molecule, i.e. keep ss-diquark together
c. are there 2 strange quarks
      prod=if(0)*if(1)*if(2) 
      if (prod.eq.9.or.prod.eq.18.or.prod.eq.27) then
c. find the non-strange quark inons
       inons=0
       do i=0,2
        if(if(i).lt.3) inons=i
       enddo
       ifb=if(inons)
       ifa=1000*if(mod(inons+1,3))+100*if(mod(inons+2,3))
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine id2ityp(id,mass,it,iz)
c
cinput id    : quarkcode 
cinput mass  : mass of resonance 
c
coutput: it : UrQMD particle ID
coutput: iz : $2\cdot I_3$ of particle
c
c     returns UrQMD ID ({\tt it}) and isospin z-component (times 2) ({\tt IZ}) 
c     from quark id
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'comres.f'
      include 'options.f'
      integer IT,IZ,id,idabs,i,j,k,jspin,idloc,whichres
      real*8 mass,pardelt,mminit,dm,dmold,massit
      integer iit,irun

c.. extract quark id's and spin
      IDABS=IABS(ID)
      I=IDABS/1000
      J=MOD(IDABS/100,10)
      K=MOD(IDABS/10,10)
      JSPIN=MOD(IDABS,10)

c single diquarks are not popular in urqmd:
      IF(ID.NE.0.AND.MOD(ID,100).EQ.0) GO TO 222
c single quarks are not popular in urqmd:
      IF(J.EQ.0) GO TO 222

c mesons:
       if(I.EQ.0) then
c calculate isospin from quark content:
       IZ=0
       if(j.lt.3) IZ=3-2*j
       if(k.lt.3) IZ=IZ+(2*k-3)
       IZ=isign(1,id)*IZ
 
       idloc=idabs-jspin
       if((idloc.eq.110).or.(idloc.eq.120))then
c..triplett (e.g. pion)
          IT= mlt2it(jspin*4+1)
          return
       else if(idloc.eq.220)then
c..1st singlett (e.g. eta)
          IT= mlt2it(jspin*4+3)
          return
       else if(idloc.eq.130.or.idloc.eq.230)then
c..strange doublett (e.g. (anti-)kaon)
          IT=isign(mlt2it(jspin*4+2),id)
          return
       else if(idloc.eq.330)then
c..2nd singlett (e.g. eta')
          IT=mlt2it(jspin*4+4)
          return
c...charmed mesons
      else if(abs(j).eq.4.and.abs(k).eq.4)then
      it=135
      iz=0
      return
      else if(abs(j).eq.4.or.abs(k).eq.4)then        
         if((j.eq.4.and.k.eq.1).or.(k.eq.4.and.j.eq.1))then
          if(id.gt.0)then
c..D0
           it=133       
           iz=-1
           return
          else
c..D0Bar
           it=-133      
           iz=1
           return
          endif
         else if((j.eq.4.and.k.eq.2).or.(k.eq.4.and.j.eq.2))then
          if(id.gt.0)then
c..D+     
          it=133      
          iz=1
          return
         else
c..D-
          it=-133      
          iz=-1
          return
         endif
        endif
      it=133
      iz=0
      return
c

       else
          goto 222
       endif

      else
c.. baryons:

c calculate isospin from quark content:
       IZ=0
       if(i.lt.3) IZ=3-2*i
       if(j.lt.3) IZ=IZ+3-2*j
       if(k.lt.3) IZ=IZ+3-2*k
       IZ=isign(1,id)*IZ
c spin 1/2 baryons (all nucleon-resonances are treated here!)
      if(jspin.eq.0)then
c (anti-)nucleons and resonances
       if(idabs.eq.1120.or.idabs.eq.1220)then
c mass below parnuc -> nucleon (ground state), mass above parnuc ->N*
         if(mass.lt.mminit(minnuc+1))then
           IT=isign(minnuc,id)
           return
         else
           IT=isign(whichres(dble(mass),1),id)
           return
         endif
       else if(idabs.eq.2130)then
c lambda
           iit=minlam
           if (mass.gt.1d0)then
           dmold=1d30
           do irun = minlam,maxlam
              dm=abs(massit(irun)-mass)
              if (dm.le.dmold)then
                 dmold=dm
                 iit=irun
              end if
           end do
           end if
           IT=isign(iit,id)
           return
       else if(idabs.eq.1230.or.idabs.eq.1130.or.idabs.eq.2230)then
c sigma
           iit=minsig
           if (mass.gt.1d0)then
           dmold=1d30
           do irun = minsig,maxsig
              dm=abs(massit(irun)-mass)
              if (dm.le.dmold)then
                 dmold=dm
                 iit=irun
              end if
           end do
           end if
           IT=isign(iit,id)
           return
       else if(idabs.eq.1330.or.idabs.eq.2330)then
c cascade
           iit=mincas
           if (mass.gt.1d0)then
           dmold=1d30
           do irun = mincas,maxcas
              dm=abs(massit(irun)-mass)
              if (dm.le.dmold)then
                 dmold=dm
                 iit=irun
              end if
           end do
           endif
           IT=isign(iit,id)
           return
       else
           goto 222
       endif
c spin 3/2 baryons (all delta-resonances are treated here!)
      else if(jspin.eq.1)then
       if(idabs.eq.1111.or.idabs.eq.1121.or.
     &    idabs.eq.1221.or.idabs.eq.2221)then
c if mass below pardelt -> Delta1232, otherwise: Delta-resonance
         pardelt=1.45
         if(mass.lt.pardelt)then
c delta 1232
           IT=isign(whichres(dble(mass),0),id)
           return
         else
c delta resonance
           IT=isign(whichres(dble(mass),2),id)          
           return
         endif
       else if(idabs.eq.1231.or.idabs.eq.1131.or.idabs.eq.2231)then
c sigma*
           IT=isign(minsig+1,id)
           return
       else if(idabs.eq.1331.or.idabs.eq.2331)then
c cascade*
           IT=isign(mincas+1,id)
           return
       else if(idabs.eq.3331)then
c omega
           IT=isign(minome,id)
           return
       else
        goto 222
       endif
c higher spin baryons include here:

      else
       goto 222
      endif
      endif
      
 222  continue
      write(6,*)'! ID=',id,' can not be converted into ityp'
      write(6,*)'I=',i,'J=',j,'K=',k,'spin=',jspin
      RETURN
      END

chp new subroutine for id2ityp

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine id2itypnew(id,mass,it,iz)
c
cinput id    : quarkcode 
cinput mass  : mass of resonance 
c
coutput: it : UrQMD particle ID
coutput: iz : $2\cdot I_3$ of particle
c
c     returns UrQMD ID ({\tt it}) and isospin z-component (times 2) ({\tt IZ}) 
c     from quark id
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'comres.f'
      include 'options.f'
      integer IT,IZ,id,idabs,i,j,k,jspin,idloc,itclass
      real*8 mass,dm,dmold,massit
      integer irun, idmin,idmax
      
      real*8 ranf
      external ranf
            
c.. extract quark id's and spin
      IDABS=IABS(ID)
      I=IDABS/1000
      J=MOD(IDABS/100,10)
      K=MOD(IDABS/10,10)
      JSPIN=MOD(IDABS,10)

c single diquarks are not popular in urqmd:
      IF(ID.NE.0.AND.MOD(ID,100).EQ.0) GO TO 222
c single quarks are not popular in urqmd:
      IF(J.EQ.0) GO TO 222

      idloc=idabs-jspin

c mesons:
      if(I.EQ.0) then
c calculate isospin from quark content:
       IZ=0
       if(j.lt.3) IZ=3-2*j
       if(k.lt.3) IZ=IZ+(2*k-3)
       IZ=isign(1,id)*IZ
 
       if((idloc.eq.110).or.(idloc.eq.120))then
c..triplett (e.g. pion)
          itclass=1
c          IT= mlt2it(jspin*4+1)
       else if(idloc.eq.220)then
c..1st singlett (e.g. eta)
          itclass=3
c          IT= mlt2it(jspin*4+3)
       else if(idloc.eq.130.or.idloc.eq.230)then
c..strange doublett (e.g. (anti-)kaon)
          itclass=2
c          IT=isign(mlt2it(jspin*4+2),id)
       else if(idloc.eq.330)then
c..2nd singlett (e.g. eta')
          itclass=4
c          IT=mlt2it(jspin*4+4)
       else
          goto 222
       endif
       
c..pick up the mass-closest resonance
       dm=abs(mass-massit(mlt2it(itclass)))
       dmold=dm
       do jspin=0,7
          irun = mlt2it(jspin*4+itclass)
          if(itclass.eq.2) irun=isign(irun,id)
          dm=abs(mass-massit(irun))
          if(dm.le.dmold)then
             IT=irun
             dmold=dm
          endif
       enddo
       
       return

      else
c.. baryons:

c calculate isospin from quark content:
       IZ=0
       if(i.lt.3) IZ=3-2*i
       if(j.lt.3) IZ=IZ+3-2*j
       if(k.lt.3) IZ=IZ+3-2*k
       IZ=isign(1,id)*IZ
      
       if(idloc.eq.1120.or.idloc.eq.1220)then
c N* or D*             !FIXME
          if(ranf(0).lt.0.5)then
c N*, ASSUME there are equal numbers of N* and D* resonances
             idmin=minnuc
             idmax=maxnuc
          else
c D+,D0
             idmin=mindel
             idmax=maxdel
          endif
       else if(idloc.eq.1110.or.idloc.eq.2220)then
c D++,D-
          idmin=mindel
          idmax=maxdel
       else if(idloc.eq.1230)then
c L or S0              !FIXME
          if(ranf(0).lt.0.3)then
c Lambda !lambda is suppressed respect to sigma0?
             idmin=minlam
             idmax=maxlam
          else
c Sigma0
             idmin=minsig
             idmax=maxsig
          endif
       else if(idloc.eq.1130.or.idloc.eq.2230)then
c Sigam+,Sigma-
          idmin=minsig
          idmax=maxsig
       else if(idloc.eq.1330.or.idloc.eq.2330)then
c Cascade
          idmin=mincas
          idmax=maxcas
       else if(idloc.eq.3330)then
c Omega
          idmin=minome
          idmax=minome
       else
          goto 222
       endif
 
c..pick up the mass-closest resonance
       dm=abs(mass-massit(idmin))
       dmold=dm
       do irun=idmin,idmax
          dm=abs(mass-massit(irun))
          if(dm.le.dmold)then
             IT=isign(irun,id)
             dmold=dm
          endif
       enddo
       
       return
      endif
               
 222  continue
      write(6,*)'! ID=',id,' can not be converted into ityp'
      write(6,*)'I=',i,'J=',j,'K=',k,'spin=',jspin
      RETURN
      END

chp end of changes



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine probitypres(itypin,itypout)
c
cinput itypin : ityp of groundstate
c
coutput: itypout : ityp of groundstate plus resonances 
c
c     returns a new ityp which also includes resonce ityps, this is
c     necessary to include hyperon resonces as long as no proper getmass
c     for hyperons existst.
c     the probabilities are chooses according to a exp(delta m/b) distribution.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'comres.f'

      integer itypin,itypout,ss,i,itypinn
      real*8 prob(1:200),probres,norm,deltam,massit

      ss=itypin/abs(itypin)
      itypinn=abs(itypin)
     
c.. the lambda 
      if (itypinn.ge.minlam.and.itypinn.le.maxlam)then
c get probabilities and norm
         norm=0d0
         do i=minlam+1,maxlam
            deltam=abs(massit(i)-massit(minlam))
            prob(i)=probres(deltam,massit(minlam),i)
            norm=norm+prob(i)
         end do
         do i=minlam+1,maxlam
            prob(i)=prob(i)/norm
         end do
         call findityp(prob,itypout,minlam+1,maxlam)
         itypout=itypout*ss
         return
      end if

c.. the sigma
      if (itypinn.ge.minsig.and.itypinn.le.maxsig)then
c get probabilities and norm
         norm=0d0
         do i=minsig+1,maxsig
            deltam=abs(massit(i)-massit(minsig))
            prob(i)=probres(deltam,massit(minsig),i)
c... do not generate masses for minsig+1 they are treated ecplicitely
            if (i.eq.minsig+1) prob(i)=0d0
            norm=norm+prob(i)
         end do
         do i=minsig+1,maxsig
            prob(i)=prob(i)/norm
         end do
         call findityp(prob,itypout,minsig+1,maxsig)
         itypout=itypout*ss
         return
      end if

c.. the cascades
      if (itypinn.ge.mincas.and.itypinn.le.maxcas)then
c get probabilities and norm
            norm=0d0
         do i=mincas+1,maxcas
            deltam=abs(massit(i)-massit(mincas))
            prob(i)=probres(deltam,massit(mincas),i)
c... do not generate masses for mincas+1 they are treated ecplicitely
            if (i.eq.mincas+1) prob(i)=0d0
            norm=norm+prob(i)
         end do
         do i=mincas+1,maxcas
            prob(i)=prob(i)/norm
         end do
         call findityp(prob,itypout,mincas+1,maxcas)
         itypout=itypout*ss
         return
      end if  
      write(*,*)'itypin, itypout',itypin,itypout
      write(*,*)'Error in probitypres!'
      stop 137 
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function probres(dm,minmass,it)
c
cinput dm : mass difference to groundstate
cinput minmass : mass of groundstate
cinput it : ityp
c
coutput: probres: unnormalized probability for this state
c
c     returns probabilty for a higher mass state according to
c     exponential distribution with degeneracy
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real*8 dm,minmass,j,g,T
      integer it, getspin

c.. assume temperature T=170 MeV (from statistical model, e.g. becattini,heinz)      
      T=0.170d0

c.. get spin*2
      j=1d0*getspin(it,1)

      if (j.lt.0) then
c.. take spin into account via J= m^2 and deg. g=2j+1 (regge theory)     
       j=2d0*(minmass+dm)**2
      end if

      g=(j+1d0)

      probres=g*exp(-(dm/T))
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine findityp (p,it,mini,maxi)
c
cinput p : array with normalized probabilities
cinput mini : lowest index 
cinput maxi : largest index
c
coutput: it: ityp according to probability
c
c     returns a new ityp according to probabilties defined in p
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real*8 p(1:200),y,ranf
      integer it,mini,maxi,ix
      
      it=mini
      ix=0
      y=0d0

 1    continue
      ix=int(mini+int((maxi-mini+1)*ranf(0)))
      y=ranf(0)
      if (y.lt.p(ix)) then 
        it=ix
        return
      end if
      goto 1
      end

C---------------------------------------------------------------------------
C                        THE END








