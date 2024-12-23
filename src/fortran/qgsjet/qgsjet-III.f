C======================================================================C 
C                                                                      C
C    QQQ        GGG      SSSS    JJJJJJJ   EEEEEEE   TTTTTTT     III   C
C   Q   Q      G   G    S    S         J   E            T        III   C
C  Q     Q    G         S              J   E            T        III   C
C  Q     Q    G   GGG    SSSS          J   EEEEE        T    ==  III   C
C  Q   Q Q    G     G        S         J   E            T        III   C
C   Q   Q      G   G    S    S    J   J    E            T        III   C
C    QQQ        GGG      SSSS      JJJ     EEEEEEE      T        III   C
C                                                                      C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C                  QUARK - GLUON - STRING - JET - III MODEL            C
C                                                                      C
C                HIGH ENERGY HADRON INTERACTION PROGRAM                C
C                                                                      C
C                                  BY                                  C
C                                                                      C
C                           S. OSTAPCHENKO                             C
C                                                                      C
C                        Hamburg University                            C
C                 II Institute for Theoretical Physics                 C
C                 e-mail: sergey.ostapchenko@desy.de                   C
C=======================================================================

c Publications to cite: PRD 109 (2024) 034002; PRD 109 (2024) 094019
c Last modified: 04.10.2024
c Regarding secondary particle types, see the procedure qgreg
c NB: because of a final model tuning, certain parameter values differ
c from what was specified in PRD 109 (2024) 094019

c=============================================================================
      subroutine qgset
c-----------------------------------------------------------------------------
c common model parameters setting
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      character*7 ty
      character*2 tyq
      parameter(iapmax=208,nfock=3)
      dimension rratio(nfock),rq0(nfock)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr3/  rmin,emax,eev
      common /qgarr6/  pi,bm,amws
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),b
      common /qgarr8/  pt2w,be(4),dc(5),deta,drho,almpt,ptdif
      common /qgarr10/ am0,amn,amk,amsig,amlam,ameta
      common /qgarr11/ b10
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr19/ wpiex(2,3)
      common /qgarr20/ spmax
      common /qgarr21/ dmmin(3),dmres(3),wdres(3)
      common /qgarr26/ factk,fqscal
      common /qgarr28/ arr(5),alpq
      common /qgarr32/ pt2str,pt2rem(3)
      common /qgarr38/ htfac
      common /qgarr40/ apipr,alphaf,bpi,breg
      common /qgarr41/ ty(6)
      common /qgarr42/ tyq(16)
      common /qgarr43/ moniou
      common /qgarr51/ epsxmn
      common /opt/     jopt
      common /qgdebug/ debug
      common /qgsIIInex1/xan(iapmax,3),xbn(iapmax,3)  !used to link with nexus
     *,bqgs,bmaxqgs,bmaxnex,bminnex

      moniou=6             !output channel for debugging
      debug=0              !debugging level
                           !(0 - no debugging, 1 - very geheral,
                           !2 - more detailed, 3 - all function calls,
                           !4 - all returned values, 5 - technical)
      if(debug.ge.1)write (moniou,210)
      
      bqgs=0.d0            !used to link with nexus
      bmaxqgs=0.d0         !used to link with nexus
      bmaxnex=-1.d0        !used to link with nexus
      bminnex=0.d0         !used to link with nexus

      jopt=1               !parameter option

      if(jopt.eq.1)then    !tunable parameters !htwn1
c soft Pomeron parameters
       dels=0.21d0         !overcriticality
       alfp=0.21d0         !trajectory slope
c multi-Pomeron vertex parameters
       r3p=0.024d0         !triple-Pomeron coupling (/4/pi)
       sgap=exp(1.5d0)     !minimal rap-gap between 3P-vertices
c coupling to DGLAP
       qt0=2.d0            !q**2 cutoff
       rr=0.11d0           !normalization of parton density in the soft pomeron
       beth(1)=2.d0        !gluon distr. hardness (pion)
       beth(2)=4.d0        !gluon distr. hardness (proton)
       beth(3)=2.d0        !gluon distr. hardness (kaon)
       dgqq=0.28d0         !sea quark/gluon relative weight
       htfac=3.d0          !normalization for ht-effects  
       factk=1.5d0         !k-factor value
c hadron-Pomeron vertices
       frg=2.02d0          !Pomeron coupling per unit transverse area
       rq0(1)=0.85d0       !vertex slope for largest diffr. eigenst. (pion)
       rq0(2)=1.86d0       !vertex slope for largest diffr. eigenst. (proton)
       rq0(3)=0.59d0       !vertex slope for largest diffr. eigenst. (kaon)
       rratio(1)=0.2d0     !ratio smallest/largest diffr. eigenst. (pion)
       rratio(2)=0.11d0    !ratio smallest/largest diffr. eigenst. (proton)
       rratio(3)=0.28d0    !ratio smallest/largest diffr. eigenst. (kaon)

c parameters for soft/hard fragmentation:

       qtf=0.15d0          !q**2 cutoff for timelike cascades
       almpt=1.55d0        !string fragmentation parameter
       pt2w=0.45d0         !switching to 2-particle decay (kinematic parameter)
       alpq=0.7d0          !constituent (anti-)quark exponent
c relative probabilities for pion/f exchange processes
       wpiex(1,1)=0.37d0   !pion-pi
       wpiex(2,1)=0.15d0   !pion-f
       wpiex(1,2)=0.37d0   !proton-pi
       wpiex(2,2)=0.08d0   !proton-f
       wpiex(1,3)=0.11d0   !kaon-pi
       wpiex(2,3)=0.d0     !kaon-f
c dc(i) - relative probabilities for qq~(qqq~q~)-pair creation from vacuum
       dc(1)=0.095d0       !udu~d~
       dc(2)=0.125d0       !ss~
       dc(3)=0.22d0        !ss~ (intrinsic)
       dc(4)=0.038d0       !us/ds
c parameters for pt-distributions of string ends & hadron remnants
       pt2str=0.15d0       !string ends
       pt2rem(1)=0.2d0     !pion
       pt2rem(2)=0.72d0    !proton
       pt2rem(3)=0.2d0     !kaon
c be(i) - parameters for pt-distributions
       be(1)=0.21d0        !uu~(dd~)
       be(2)=0.33d0        !qqq~q~
       be(3)=0.28d0        !ss~
       bpi=0.15d0          !slope for pion-nucleon form factor
       breg=4.d0           !slope for f-nucleon form factor   
       ptdif=0.43d0        !diffractive momentum transfer

c parameters for nuclear spectator part fragmentation:
       rmin=3.35d0    !coupling radius squared (fm^2)
       emax=0.11d0    !relative critical energy ( / <E_ex>, <E_ex>~12.5 MeV )
       eev=.025d0     !relative evaporation energy ( / <E_ex>, <E_ex>~12.5 MeV )
      else
       stop'wrong option!!!'
      endif
      g3p=(4.d0*3.1416d0*r3p)/dels !corresponds to critical renormalized pomeron

c Pomeron-hadron couplings for diffractive eigenstates
      do icz=1,3
       do ifock=1,nfock
        cc(ifock,icz)=1.d0/dble(nfock)   !equal weigts for diffr. eigenstates
        rq(ifock,icz)=rq0(icz)*rratio(icz)**(dble(ifock-1)
     *  /dble(nfock-1))
        fp(ifock,icz)=frg*rq(ifock,icz)
       enddo
      enddo

!other parameters and constants:

      spmax=1.d11          !max energy squared for tabulations
      delh=0.19d0          !effective exponent for MC sampling (technical)
      if(delh.eq.dels)stop'delh=dels!!!'
      epsxmn=0.01d0        !pt-resolution scale (technical)
      alm=0.04d0           !lambda_qcd squared
      fqscal=4.d0          !factor for fact. scale (Mf^2=p_t^2/fqscal)
      deta=0.11111d0       !ratio of etas production to all pions (1/9)
      drho=0.3333d0        !proportion of rho (~1/3)
      dc(5)=0.0d0          !to switch off charmed particles set to 0.000
c auxiliary constants
      b10=.43876194d0      !initial value of the pseudorandom sequence
      pi=3.1416d0          !pi-value
      amws=0.523d0         !diffusive radius for saxon-wood density
c regge intercepts for the uu~, qqq~q~, us~, uc~ trajectories
      arr(1)=0.5d0         !qq~-trajectory
      arr(2)=-0.5d0        !qqq~q~-trajectory
      arr(3)=0.d0          !us~-trajectory
      alphaf=0.7d0         !intercept for f-trajectory
      apipr=0.9d0          !reggeon slope
c lowest resonance masses for low-mass excitations
      dmmin(1)=0.76d0      !rho
      dmmin(2)=1.23d0      !delta
      dmmin(3)=0.89d0      !K*
c mass and width for resonance contribution to low mass diffraction
      dmres(1)=1.23d0      !pion
      dmres(2)=1.44d0      !proton
      dmres(3)=1.27d0      !kaon
      wdres(1)=0.3d0       !pion
      wdres(2)=0.3d0       !proton
      wdres(3)=0.1d0       !kaon
c proton, kaon, pion, lambda, sigma, eta masses
      amn=0.93827999d0
      amk=0.496d0
      am0=0.14d0
      amlam=1.116d0
      amsig=1.19d0
      ameta=0.548d0
c initial particle classes
      ty(1)='pion   '
      ty(2)='nucleon'
      ty(3)='kaon   '
c parton types
      tyq(1)='DD'
      tyq(2)='UU'
      tyq(3)='C '
      tyq(4)='S '
      tyq(5)='UD'
      tyq(6)='D '
      tyq(7)='U '
      tyq(8)='g '
      tyq(9)='u '
      tyq(10)='d '
      tyq(11)='ud'
      tyq(12)='s '
      tyq(13)='c '
      tyq(14)='uu'
      tyq(15)='dd'
      if(debug.ge.2)write (moniou,202)

210   format(2x,'qgset - common model parameters setting')
202   format(2x,'qgset - end')
      return
      end
      
c=============================================================================
      subroutine qgaini( DATDIR )
c-----------------------------------------------------------------------------
c common initialization procedure
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      CHARACTER DATDIR*(132)
      real qggamfun
      integer debug
      character *7 ty
      logical lcalc
      parameter(iapmax=208,nfock=3)
      dimension evs(40,100,3,2),ixemax(40,3,2),gz0(5),gz1(3),fann(17)
     *,qfan0(11,13),qnorm(3),v1p(5),qxpt(7)
     *,sig3ht(8),sig3htm(5),vht(5),vhtqp(4),vhtqt(4),vht0(5)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr10/ am(6)
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr20/ spmax
      common /qgarr24/ qpomr(11,11,288,nfock**2,77)
     *,dhteik(11,11,216,nfock**2,60)
      common /qgarr25/ ahv(3)
      common /qgarr26/ factk,fqscal
      common /qgarr27/ qlegi(31,11,nfock,3,6),qfanu(31,11,11,3*nfock,5)
     *,qfanc(31,11,11,39,13*nfock)
      common /qgarr28/ arr(5),alpq
      common /qgarr29/ cstot(40,40,240)
      common /qgarr30/ cs0(40,40,240)
      common /qgarr31/ csborn(40,240)
      common /qgarr33/ fsud(10,2)
      common /qgarr34/ qrt(10,101,2)
      common /qgarr38/ htfac
      common /qgarr39/ qpomi(31,11,18),qpomis(31,11,11,11,5)
      common /qgarr41/ ty(6)
      common /qgarr43/ moniou
      common /qgarr47/ gsect(10,3,6)
      common /qgarr48/ qgsasect(10,6,6)
      common /qgarr51/ epsxmn
      common /qgarr52/ evk(40,40,100,3,2)
      common /qgarr56/ pdfr(31,11,11,9,4*nfock)
      common /qgarr57/ qloopr(31,11,11,2,2)
      common /qgarr73/ htfacm
      common /qgarr74/ feikht(20,10,11,66,24*nfock),ffhtm(20,10,12,6,8)
     *,flhtm(20,12,6,2,2)
c auxiliary common blocks to calculate hadron-nucleus cross-sections
      common /arr1/    trnuc(56),twsnuc(56),twbnuc(56)
      common /arr3/    x1(7),a1(7)
      common /arr6/    x3(8),a3(8)
      common /opt/     jopt
      common /qgdebug/ debug
      character*500 fnIIIdat,fnIIIncs                         !to link with nexus
      common /version/   version                            !to link with nexus
      common/qgsIIIfname/ fnIIIdat, fnIIIncs, ifIIIdat, ifIIIncs !to link with nexus
      common/qgsIIInfname/nfnIIIdat, nfnIIIncs                 !to link with nexus
      common/producetab/ producetables                      !to link with CRMC
      logical producetables

      if(debug.ge.1)write (moniou,210)
      version = 301

c-------------------------------------------------
      write(*,100)
 100  format(' ',
     *           '====================================================',
     *     /,' ','|                                                  |',
     *     /,' ','|         QUARK GLUON STRING JET - III MODEL       |',
     *     /,' ','|                                                  |',
     *     /,' ','|         HADRONIC INTERACTION MONTE CARLO         |',
     *     /,' ','|                        BY                        |',
     *     /,' ','|                 S. OSTAPCHENKO                   |',
     *     /,' ','|                                                  |',
     *     /,' ','|             e-mail: sergei@tf.phys.ntnu.no       |',
     *     /,' ','|                                                  |',
     *     /,' ','|                   Version III-04                 |',
     *     /,' ','|                                                  |',
     *     /,' ','| Publications to be cited when using this program:|',
     *     /,' ','| S. Ostapchenko, PRD 109 (2024) 034002; 094019    |',
     *     /,' ','|                                                  |',
     *     /,' ','| Any modification has to be approved by the author|',
     *     /,' ','====================================================',
     *     /)

c large-x behavior for valence quarks (GRS-pion and GRV-proton)
      ahv(1)=.383d0+.624d0*dlog(dlog(qt0/.204d0**2)
     */dlog(.26d0/.204d0**2))
      ahv(3)=ahv(1)
      sq=dlog(dlog(qt0/.232d0**2)/dlog(.23d0/.232d0**2))
      ahv(2)=2.997d0+.753d0*sq-.076d0*sq*sq
c valence quark momentum share 
      do icz=1,2      !only pions and nucleons
       qnorm(icz)=0.d0
       do i=1,7
       do m=1,2
        tp=1.d0-(.5d0+x1(i)*(m-1.5d0))**(2.d0/3.d0)
        xp=1.d0-tp**(1.d0/(1.d0+ahv(icz)))       
        qnorm(icz)=qnorm(icz)+a1(i)*(qggrv(xp,qt0,icz,1)
     *  +qggrv(xp,qt0,icz,2))/dsqrt(1.d0-tp)
       enddo
       enddo
       qnorm(icz)=qnorm(icz)/(1.d0+ahv(icz))/3.d0
      enddo
      qnorm(3)=qnorm(1)
c fix gluon PDFs by momentum conservation
      do icz=1,3
       do icdp=1,nfock
        bbbi(icdp,icz)=0.d0
       enddo
      enddo
      bbbpom=0.d0
      iii=0
1     iii=iii+1
      gnorm=0.d0
      seanrm=0.d0
      do i=1,8
      do m=1,2
       xxg=(.5d0+x3(i)*(m-1.5d0))**(1.d0/(1.d0-dels))
       gnorm=gnorm+a3(i)*qgppdi(xxg,0)
       seanrm=seanrm+a3(i)*qgppdi(xxg,1)
      enddo
      enddo
      gnorm=gnorm/(1.d0-dels)*g3p*rr*2.d0*pi
      seanrm=seanrm/(1.d0-dels)*g3p*rr*2.d0*pi
      dbbb=(1.d0-gnorm-seanrm)/(4.d0*pi*rr*g3p)/qggamfun(real(2.d0
     *-dels))*qggamfun(real(3.d0-dels+beth(1)))/qggamfun(real(1.d0
     *+beth(1)))/(1.d0+bbbpom*(2.d0-dels)/(3.d0-dels+beth(1))
     *+.5d0*bbbpom*(bbbpom-1.d0)*(2.d0-dels)
     **(3.d0-dels)/(3.d0-dels+beth(1))/(4.d0-dels+beth(1)))
      bbbpom=bbbpom+dbbb
      if(abs(1.d0-gnorm-seanrm).gt.1.d-4)goto 1

      do icz=1,3
       do ic=1,nfock
        iii=0
2       iii=iii+1
        gnorm=0.d0
        seanrm=0.d0
        do i=1,8
        do m=1,2
         xxg=(.5d0+x3(i)*(m-1.5d0))**(1.d0/(1.d0-dels))
         gnorm=gnorm+a3(i)*qgppdc(xxg,0,ic,icz)
         seanrm=seanrm+a3(i)*qgppdc(xxg,1,ic,icz)
        enddo
        enddo
        gnorm=gnorm/(1.d0-dels)*fp(ic,icz)*rr*2.d0*pi
        seanrm=seanrm/(1.d0-dels)*fp(ic,icz)*rr*2.d0*pi
        bbbi(ic,icz)=bbbi(ic,icz)+(1.d0-qnorm(icz)-gnorm-seanrm)
     *   /4.d0/pi/rr/fp(ic,icz)*qggamfun(real(3.d0-dels+beth(icz)))
     *   /qggamfun(real(2.d0-dels))/qggamfun(real(1.d0+beth(icz)))
     *   /(1.d0+bbbi(ic,icz)*(2.d0-dels)/(3.d0-dels+beth(icz))
     *   +.5d0*bbbi(ic,icz)*(bbbi(ic,icz)-1.d0)*(2.d0-dels)
     *   *(3.d0-dels)/(3.d0-dels+beth(icz))/(4.d0-dels+beth(icz)))
        bbbi(ic,icz)=max(-20.d0,bbbi(ic,icz))
        bbbi(ic,icz)=min(10.d0,bbbi(ic,icz))
        if(abs(qnorm(icz)+gnorm+seanrm-1.d0).gt.1.d-4)goto 2
       enddo
      enddo

c-----------------------------------------------------------------------------
c     reading cross sections from the file
c      print *,ifIIIdat,DATDIR(1:INDEX(DATDIR,' ')-1)
c     *       ,fnIIIdat(1:nfnIIIdat)
      if(ifIIIdat.ne.1)then
       inquire(file=DATDIR(1:INDEX(DATDIR,' ')-1)//'qgsjetIII.dat'
     *        ,exist=lcalc)
      else
       inquire(file=fnIIIdat(1:nfnIIIdat),exist=lcalc) !used to link with nexus
      endif
      lzmaUse=0
      if(lcalc)then
        if(debug.ge.2)write (moniou,205)
         if(ifIIIdat.ne.1)then
            open(1,file=DATDIR(1:INDEX(DATDIR,' ')-1)//'qgsjetIII.dat'
     *           ,status='old')
         else                   !used to link with nexus
            if (LEN(fnIIIdat).gt.6.and.
     *           fnIIIdat(nfnIIIdat-4:nfnIIIdat) .eq. ".lzma") then
               lzmaUse=1
               call LzmaOpenFile(fnIIIdat(1:nfnIIIdat))
            else
               open(ifIIIdat,file=fnIIIdat(1:nfnIIIdat),status='old')
            endif
         endif

         if (lzmaUse.ne.0) then

          if(debug.ge.0)write (moniou,214) 'qgsjetIII.dat.lzma'

          call LzmaFillArray(csborn,size(csborn))
          call LzmaFillArray(cs0,size(cs0))
          call LzmaFillArray(cstot,size(cstot))
          call LzmaFillArray(evk,size(evk))
          call LzmaFillArray(qpomi,size(qpomi))
          call LzmaFillArray(qpomis,size(qpomis))
          call LzmaFillArray(qloopr,size(qloopr))
          call LzmaFillArray(qlegi,size(qlegi))
          call LzmaFillArray(qfanu,size(qfanu))
          call LzmaFillArray(qfanc,size(qfanc))
          call LzmaFillArray(pdfr,size(pdfr))
          call LzmaFillArray(qpomr,size(qpomr))
          call LzmaFillArray(dhteik,size(dhteik))
          call LzmaFillArray(feikht,size(feikht))
          call LzmaFillArray(ffhtm,size(ffhtm))
          call LzmaFillArray(flhtm,size(flhtm))
          call LzmaFillArray(gsect,size(gsect))
          call LzmaFillArray(fsud,size(fsud))
          call LzmaFillArray(qrt,size(qrt))
          call LzmaCloseFile()
        else
          if(debug.ge.0)write (moniou,214) 'qgsjetIII.dat'
          read (1,*)csborn,cs0,cstot,evk,qpomi,qpomis,qloopr,qlegi,qfanu
     * ,qfanc,pdfr,qpomr,dhteik,feikht,ffhtm,flhtm,gsect,fsud,qrt
          close(1)
        endif
      
       if(debug.ge.0)write (moniou,202)
       
c fix gluon PDF by momentum conservation
       do icz=1,3
        do icdp=1,nfock
         iii=0
3        iii=iii+1
         gnorm=0.d0
         seanrm=0.d0
         do i=1,8
         do m=1,2
          xxg=(.5d0+x3(i)*(m-1.5d0))**(1.d0/(1.d0-dels))
          gnorm=gnorm+a3(i)*qgpdff(xxg,icdp,icz,1)*xxg**dels
          seanrm=seanrm+a3(i)*qgpdff(xxg,icdp,icz,2)*xxg**dels
         enddo
         enddo
         gnorm=gnorm/(1.d0-dels)/2.d0
         seanrm=seanrm/(1.d0-dels)/2.d0
         bbbi(icdp,icz)=bbbi(icdp,icz)+(1.d0-gnorm-seanrm)
     *   /4.d0/pi/rr/fp(icdp,icz)*qggamfun(real(3.d0-dels+beth(icz)))
     *   /qggamfun(real(2.d0-dels))/qggamfun(real(1.d0+beth(icz)))
     *   /(1.d0+bbbi(icdp,icz)*(2.d0-dels)/(3.d0-dels+beth(icz))
     *   +.5d0*bbbi(icdp,icz)*(bbbi(icdp,icz)-1.d0)*(2.d0-dels)
     *   *(3.d0-dels)/(3.d0-dels+beth(icz))/(4.d0-dels+beth(icz)))
         bbbi(icdp,icz)=max(-20.d0,bbbi(icdp,icz))
         bbbi(icdp,icz)=min(10.d0,bbbi(icdp,icz))
         if(abs(1.d0-gnorm-seanrm).gt.1.d-4)goto 3
        enddo
       enddo
       goto 15

      elseif(.not.producetables)then
        write(moniou,*) "Missing qgsjetIII.dat file !"
        write(moniou,*) "Please correct the defined path ",
     &"or force production ..."
        stop
      endif

c--------------------------------------------------
c qcd evolution and qcd ladder cross sections
      if(debug.ge.0)write (moniou,201)
      do i=1,40
      do m=1,3
      do k=1,2
       ixemax(i,m,k)=99
      do j=1,40
      do l=1,100
       evk(i,j,l,m,k)=0.d0
      enddo
      enddo
      enddo
      enddo
      enddo

      n=1
4     n=n+1
      do m=1,3
      do k=1,2
       if(m.ne.3.or.k.ne.1)then
        do i=1,39
         if(ixemax(i,m,k).gt.0)then
          qi=spmax**((i-1)/39.d0)
          qq=qi*(spmax/qi)**(1.d0/39.d0)
          do l=1,99
           if(l.le.37)then
            xx=.1d0/(.1d0*spmax)**((37.d0-l)/36.d0)
           elseif(l.le.69)then
            xx=.1d0+.8d0*(l-37.d0)/32.d0
           else
            xx=1.d0-.1d0*(10.d0*epsxmn)**((l-69.d0)/31.d0)
           endif

           ev=qgev(qi,qq,qq,xx,m,k)/qgfap(xx,m,k)
           if(m.eq.1.and.k.eq.1.or.m.ne.1.and.k.ne.1)then
            evs(i,l,m,k)=dlog(1.d0+ev*4.5d0*qgsudx(qi,m)/qgsudx(qq,m)
     *      /dlog(dlog(qq/alm)/dlog(qi/alm)))
           else
            evs(i,l,m,k)=dlog(1.d0+ev/.3d0*(dlog(epsxmn)+.75d0)
     *      /(qgsudx(qq,1)/qgsudx(qi,1)-qgsudx(qq,2)/qgsudx(qi,2)))
           endif
          enddo
         endif
        enddo
       endif
      enddo
      enddo

      jec=0
      do m=1,3
      do k=1,2
       if(m.ne.3.or.k.ne.1)then
        do i=1,39
         if(ixemax(i,m,k).gt.0)then
          qi=spmax**((i-1)/39.d0)
          qq=qi*(spmax/qi)**(1.d0/39.d0)
          imx=ixemax(i,m,k)
          do l=imx,1,-1
           if(l.le.37)then
            xx=.1d0/(.1d0*spmax)**((37.d0-l)/36.d0)
           elseif(l.le.69)then
            xx=.1d0+.8d0*(l-37.d0)/32.d0
           else
            xx=1.d0-.1d0*(10.d0*epsxmn)**((l-69.d0)/31.d0)
           endif

           if(abs(evs(i,l,m,k)-evk(i,2,l,m,k)).gt.1.d-3)then
            evk(i,2,l,m,k)=evs(i,l,m,k)
            jec=1
           elseif(ixemax(i,m,k).eq.l)then
            ixemax(i,m,k)=l-1
           endif
          enddo
         endif
        enddo
       endif
      enddo
      enddo

      do i=1,39
       qi=spmax**((i-1)/39.d0)
       qj=qi*(spmax/qi)**(1.d0/39.d0)
       qq=qi*(spmax/qi)**(2.d0/39.d0)
       do l=99,1,-1
        if(l.le.37)then
         xx=.1d0/(.1d0*spmax)**((37.d0-l)/36.d0)
        elseif(l.le.69)then
         xx=.1d0+.8d0*(l-37.d0)/32.d0
        else
         xx=1.d0-.1d0*(10.d0*epsxmn)**((l-69.d0)/31.d0)
        endif
        do m=1,3
        do k=1,2
         if(m.ne.3.or.k.ne.1)then
          ev=(qgev(qi,qj,qq,xx,m,k)
     *    +qgevi(qi,qj,xx,m,k)*qgsudx(qq,k)/qgsudx(qj,k)
     *    +qgevi(qj,qq,xx,m,k)*qgsudx(qj,m)/qgsudx(qi,m))/qgfap(xx,m,k)
          if(m.eq.1.and.k.eq.1.or.m.ne.1.and.k.ne.1)then
           evk(i,3,l,m,k)=dlog(ev*4.5d0*qgsudx(qi,m)/qgsudx(qq,m)
     *     /dlog(dlog(qq/alm)/dlog(qi/alm)))
          else
           evk(i,3,l,m,k)=dlog(ev/.3d0*(dlog(epsxmn)+.75d0)
     *     /(qgsudx(qq,1)/qgsudx(qi,1)-qgsudx(qq,2)/qgsudx(qi,2)))
          endif
         endif
        enddo
        enddo
       enddo
      enddo
      if(jec.ne.0)goto 4

      do i=1,39
       qi=spmax**((i-1)/39.d0)
      do j=4,40
       qj=qi*(spmax/qi)**((j-2)/39.d0)
       qq=qi*(spmax/qi)**((j-1)/39.d0)
       do l=99,1,-1
        if(l.le.37)then
         xx=.1d0/(.1d0*spmax)**((37.d0-l)/36.d0)
        elseif(l.le.69)then
         xx=.1d0+.8d0*(l-37.d0)/32.d0
        else
         xx=1.d0-.1d0*(10.d0*epsxmn)**((l-69.d0)/31.d0)
        endif
        do m=1,3
        do k=1,2
         if(m.ne.3.or.k.ne.1)then
          ev=(qgev(qi,qj,qq,xx,m,k)
     *    +qgevi(qi,qj,xx,m,k)*qgsudx(qq,k)/qgsudx(qj,k)
     *    +qgevi(qj,qq,xx,m,k)*qgsudx(qj,m)/qgsudx(qi,m))/qgfap(xx,m,k)
          if(m.eq.1.and.k.eq.1.or.m.ne.1.and.k.ne.1)then
           evk(i,j,l,m,k)=dlog(ev*4.5d0*qgsudx(qi,m)/qgsudx(qq,m)
     *     /dlog(dlog(qq/alm)/dlog(qi/alm)))
          else
           evk(i,j,l,m,k)=dlog(ev/.3d0*(dlog(epsxmn)+.75d0)
     *     /(qgsudx(qq,1)/qgsudx(qi,1)-qgsudx(qq,2)/qgsudx(qi,2)))
          endif
         endif
        enddo
        enddo
       enddo
      enddo
      enddo

c--------------------------------------------------
c qcd ladder cross sections
      do i=1,40
       qi=(spmax/4.d0/fqscal)**((i-1)/39.d0)  !q^2 cutoff for born process
       s2min=qi*4.d0*fqscal                !energy threshold for 2->2 subprocess
      do m=1,4                                !parton types (1-g, 2,3,4-q)
       if(m.gt.2)then
        lmin=2
       else
        lmin=1
       endif
      do l=lmin,2                             !parton types (1-g, 2-q)
      do k=1,40
       sk=s2min*(spmax/s2min)**((k-1)/39.d0)  !c.m. energy squared
       k1=k+40*(m-1)+80*(l-1)
       csborn(i,k1)=dlog(qgborn(qi,qi,sk,m-1,l-1)) !born cross-section (2->2)
       if(.not.(csborn(i,k1).ge.0.d0.or.csborn(i,k1).lt.0.d0))stop
      enddo
      enddo
      enddo
      enddo

      do i=1,40
       qi=(spmax/4.d0/fqscal)**((i-1)/39.d0)
      do j=1,40
       qj=qi*(spmax/4.d0/fqscal/qi)**((j-1)/39.d0)
       s2min=qj*4.d0*fqscal
       smin=s2min/(1.d0-epsxmn)
      do m=1,4
       if(m.gt.2)then
        lmin=2
       else
        lmin=1
       endif
      do l=lmin,2                             !parton types (1-g, 2-q)
      do k=1,40
       sk=s2min*(spmax/s2min)**((k-1)/39.d0)
       k1=k+40*(m-1)+80*(l-1)
       tmin=qj*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qj*fqscal/sk)))
       sjtot=qgjett(qi,qj,sk,m-1,l-1)
       sjord1=qgjeto(qi,qj,sk,m-1,l-1)
       sjord2=qgjeto(qj,qi,sk,l-1,m-1)
       born=qgborn(qi,qj,sk,m-1,l-1)
       if(k.eq.1.or.j.eq.40.or.i.eq.40.or.sk.le.smin)then
        cstot(i,j,k1)=dlog(born)
        cs0(i,j,k1)=cstot(i,j,k1)
       else
        cstot(i,j,k1)=dlog(born+(sjtot+sjord1+sjord2)
     *  /(1.d0/tmin-2.d0/sk))
        cs0(i,j,k1)=dlog(born+sjord1/(1.d0/tmin-2.d0/sk))
       endif
       if(.not.(cstot(i,j,k1).ge.0.d0.or.cstot(i,j,k1).lt.0.d0))stop
       if(.not.(cs0(i,j,k1).ge.0.d0.or.cs0(i,j,k1).lt.0.d0))stop
      enddo
      enddo
      enddo
      enddo
      enddo
      
c--------------------------------------------------
c enhanced diagrams
      jrep=1
      if(debug.ge.1)write (moniou,211)
      do iy=1,31
      do iz=1,11
      do iqq=1,18
       qpomi(iy,iz,iqq)=0.d0
      enddo
      enddo
      enddo
      do iy=2,31
       sy=sgap**2*(spmax/sgap**4)**((iy-1)/30.d0)
       rp=alfp*log(sy)*4.d0*.0389d0
      do iz=2,11
       if(iz.gt.6)then
        z=.2d0*(iz-6)
        b=sqrt(-log(z)*rp)
       elseif(iz.gt.1)then
        b=sqrt(-rp*(log(0.2d0)+2.5d0*(iz-7)))
        z=exp(-b*b/rp)
       endif
       qpomi(iy,iz,1)=dlog(qgpint(sy,b*b)/qgpini(sy,b*b,0.d0,0.d0,0)
     * +1.d0)
      enddo
      enddo

c-------------------------------------------------
c loop contribution
      do iy=2,31
       sy=sgap**2*(spmax/sgap**4)**((iy-1)/30.d0)
       rp=alfp*log(sy)*4.d0*.0389d0
       do iz=2,11
       do iqq=2,7
        qpomi(iy,iz,iqq)=qpomi(iy-1,iz,iqq)
        if(iy.lt.31)qpomi(iy+1,iz,iqq)=qpomi(iy,iz,iqq)
       enddo
       enddo
       n=0
5      n=n+1
       nrep=0
       do iz=2,11
        if(iz.gt.6)then
         z=.2d0*(iz-6)
         b=sqrt(-log(z)*rp)
        elseif(iz.gt.1)then
         b=sqrt(-rp*(log(0.2d0)+2.5d0*(iz-7)))
         z=exp(-b*b/rp)
        endif
        call qgloop(sy,b*b,fann,1)
        do iqq=1,6
         if(fann(iqq).gt.0.d0)then
          qfan0(iz,iqq)=dlog(fann(iqq)/qgpini(sy,b*b,0.d0,0.d0,0))
         elseif(iy.gt.2)then
          qfan0(iz,iqq)=min(2.d0*qpomi(iy-1,iz,iqq+1)
     *    -qpomi(iy-2,iz,iqq+1),qpomi(iy-1,iz,iqq+1))
         else
          stop'loop<0: iy=2'
         endif
         if(abs(qfan0(iz,iqq)-qpomi(iy,iz,iqq+1)).gt.1.d-3)nrep=1
        enddo
       enddo
       do iz=2,11
        do iqq=2,7
         qpomi(iy,iz,iqq)=qfan0(iz,iqq-1)
        enddo
       enddo
       if(nrep.eq.1.and.n.lt.100)goto 5
      enddo
      
c-------------------------------------------------
c cut loops
      do iy=2,31
       sy=sgap**2*(spmax/sgap**4)**((iy-1)/30.d0)
       rp=alfp*log(sy)*4.d0*.0389d0
       do iz=2,11
       do iqq=8,18
        qpomi(iy,iz,iqq)=qpomi(iy-1,iz,iqq)
        if(iy.lt.31)qpomi(iy+1,iz,iqq)=qpomi(iy,iz,iqq)
       enddo
       enddo
       n=0
6      n=n+1
       nrep=0
       do iz=2,11
        if(iz.gt.6)then
         z=.2d0*(iz-6)
         b=sqrt(-log(z)*rp)
        else
         b=sqrt(-rp*(log(0.2d0)+2.5d0*(iz-7)))
         z=exp(-b*b/rp)
        endif
        call qgloop(sy,b*b,fann,2)
        do iqq=8,18
         if(fann(iqq-1).gt.0.d0)then
          qfan0(iz,iqq-7)=dlog(fann(iqq-1)/qgpini(sy,b*b,0.d0,0.d0,0))
         elseif(iy.gt.2)then
          qfan0(iz,iqq-7)=min(2.d0*qpomi(iy-1,iz,iqq)
     *    -qpomi(iy-2,iz,iqq),qpomi(iy-1,iz,iqq))
         else
          stop'loop<0: iy=2'
         endif
         if(abs(qfan0(iz,iqq-7)-qpomi(iy,iz,iqq)).gt.1.d-3)nrep=1
        enddo
       enddo
       do iz=2,11
       do iqq=8,18
        qpomi(iy,iz,iqq)=qfan0(iz,iqq-7)
       enddo
       enddo
       if(nrep.eq.1.and.n.lt.50)goto 6
      enddo
      
c-------------------------------------------------
c loop pdfs
      do iv=1,11  
       vvx=dble(iv-1)/10.d0
       if(iv.eq.1)then
        do iy=1,31
        do iz=1,11
         do iqq=1,2
         do jj=1,2
          qloopr(iy,iz,iv,iqq,jj)=0.d0
         enddo
         enddo
        enddo
        enddo
       else
        do iy=1,31
        do iz=1,11
         do iqq=1,2
         do jj=1,2
          qloopr(iy,iz,iv,iqq,jj)=qloopr(iy,iz,iv-1,iqq,jj)
         enddo
         enddo
        enddo
        enddo
       endif
            
       do iy=2,31
        sy=sgap*(spmax/sgap**2)**((iy-1)/30.d0)
        rp=alfp*dlog(sy)*4.d0*.0389d0
        do iz=2,11
         do iqq=1,2
         do jj=1,2
          qloopr(iy,iz,iv,iqq,jj)=qloopr(iy-1,iz,iv,iqq,jj)
         enddo
         enddo
        enddo
        n=0
7       n=n+1
        nrep=0
        do iz=2,11
         if(iz.gt.6)then
          z=.2d0*(iz-6)
          b=sqrt(-log(z)*rp)
         else
          b=sqrt(-rp*(log(0.2d0)+2.5d0*(iz-7)))
          z=exp(-b*b/rp)
         endif
         do iqq=1,2
         do jj=1,2
          ff=qgloopr(sy,b*b,vvx,iqq,jj)
          if(ff.gt.0.d0)then
           qfan0(iz,iqq+2*(jj-1))=dlog(ff)
          elseif(iy.gt.2)then
           qfan0(iz,iqq+2*(jj-1))=min(qloopr(iy-1,iz,iv,iqq,jj)
     *     ,2.d0*qloopr(iy-1,iz,iv,iqq,jj)-qloopr(iy-2,iz,iv,iqq,jj))
          else
           qfan0(iz,iqq+2*(jj-1))=qloopr(iy-1,iz,iv,iqq,jj)
          endif
          if(abs(qfan0(iz,iqq+2*(jj-1))-qloopr(iy,iz,iv,iqq,jj)) 
     *    .gt.1.d-3)nrep=1
         enddo
         enddo
        enddo

        do iz=2,11
         do iqq=1,2
         do jj=1,2
          qloopr(iy,iz,iv,iqq,jj)=qfan0(iz,iqq+2*(jj-1))
         enddo
         enddo
        enddo
        if(iv.gt.1.and.nrep.eq.1.and.n.lt.50)goto 7
       enddo
      enddo

c-------------------------------------------------
c cut loops with proj/targ screening corrections
      do iv=1,11
       vvx=dble(iv-1)/10.d0
      do iv1=1,11
       vvxt=dble(iv1-1)/10.d0
       do iy=1,31
       do iz=1,11
       do iqq=1,5
         qpomis(iy,iz,iv,iv1,iqq)=0.d0
       enddo
       enddo
       enddo

       do iy=2,31
        sy=sgap**2*(spmax/sgap**4)**((iy-1)/30.d0)
        rp=alfp*log(sy)*4.d0*.0389d0
        do iz=2,11
        do iqq=1,5
         qpomis(iy,iz,iv,iv1,iqq)=qpomis(iy-1,iz,iv,iv1,iqq)
         if(iy.lt.31)qpomis(iy+1,iz,iv,iv1,iqq)
     *   =qpomis(iy,iz,iv,iv1,iqq)
        enddo
        enddo

        n=0
8       n=n+1
        nrep=0
        do iz=2,11
         if(iz.gt.6)then
          z=.2d0*(iz-6)
          b=sqrt(-log(z)*rp)
         else
          b=sqrt(-rp*(log(0.2d0)+2.5d0*(iz-7)))
          z=exp(-b*b/rp)
         endif
         call qgloos(sy,b*b,vvx,vvxt,fann)  
         vi0=qgpini(sy,b*b,0.d0,0.d0,0)
         do iqq=1,5
          if(fann(iqq).gt.0.d0)then
           qfan0(iz,iqq)=dlog(fann(iqq)/vi0)
          elseif(iy.gt.2)then
           qfan0(iz,iqq)=min(2.d0*qpomis(iy-1,iz,iv,iv1,iqq)
     *     -qpomis(iy-2,iz,iv,iv1,iqq),qpomis(iy-1,iz,iv,iv1,iqq))
          else
           qfan0(iz,iqq)=qpomis(iy-1,iz,iv,iv1,iqq)
          endif
          if(abs(qfan0(iz,iqq)-qpomis(iy,iz,iv,iv1,iqq)).gt.1.d-3)
     *    nrep=1
         enddo
        enddo
        do iz=2,11
         do iqq=1,5
          qpomis(iy,iz,iv,iv1,iqq)=qfan0(iz,iqq)
         enddo
        enddo
        if(nrep.eq.1.and.n.lt.50)goto 8
       enddo
      enddo
      enddo

c-------------------------------------------------
c integrated Pomeron leg eikonals
9     if(debug.ge.1)write (moniou,212)
      do icz=1,3
      do ifock=1,nfock
      do iy=1,31
      do iz=1,11  
      do iqq=1,6
       qlegi(iy,iz,ifock,icz,iqq)=0.d0
      enddo
      enddo
      enddo
      enddo
      enddo
      do icz=1,3
      do ifock=1,nfock
       do iy=2,31
        sy=sgap**2*(spmax/sgap**3)**((iy-1)/30.d0)
        rp=(rq(ifock,icz)+alfp*log(sy))*4.d0*.0389d0
        do iz=2,11  
         do iqq=1,2
          if(iz.gt.6)then
           z=.2d0*(iz-6)
           b=sqrt(-log(z)*rp)
          elseif(iz.gt.1)then
           b=sqrt(-rp*(log(0.2d0)+2.5d0*(iz-7)))
           z=exp(-b*b/rp)
          endif
          qlegi(iy,iz,ifock,icz,iqq)=log(qglsh(sy,b*b,ifock,icz,iqq-1)
     *    /qglegi(sy,b*b,ifock,icz,0)+1.d0)
         enddo
        enddo
       enddo
      enddo
      enddo

c-------------------------------------------------
c loop-legs 
      do icz=1,3
      do ifock=1,nfock
       do iy=2,31
        sy=sgap**2*(spmax/sgap**3)**((iy-1)/30.d0)
        rp=(rq(ifock,icz)+alfp*log(sy))*4.d0*.0389d0
        do iz=2,11  
         if(iz.gt.6)then
          z=.2d0*(iz-6)
          b=sqrt(-log(z)*rp)
         elseif(iz.ne.1)then
          b=sqrt(-rp*(log(0.2d0)+2.5d0*(iz-7)))
          z=exp(-b*b/rp)
         endif
         call qglool(sy,b*b,ifock,icz,fann)
         do iqq=3,6       
          if(fann(iqq-2).gt.0.d0)then
           qlegi(iy,iz,ifock,icz,iqq)=log(fann(iqq-2)
     *     /qglegi(sy,b*b,ifock,icz,0))
          elseif(iy.gt.2)then
           qlegi(iy,iz,ifock,icz,iqq)=2.d0*qlegi(iy-1,iz,ifock,icz,iqq)
     *     -qlegi(iy-2,iz,ifock,icz,iqq)
          else
           stop'lool: iy=2!'
          endif
         enddo
        enddo
       enddo
      enddo
      enddo
      
c-------------------------------------------------
c uncut fan-contributions
      if(debug.ge.1)write (moniou,213)
      do icz=1,3
      do iv=1,11
      do ifock=1,nfock
      do iy=1,31
      do iz=1,11 
      do iqq=1,5
       qfanu(iy,iz,iv,ifock+nfock*(icz-1),iqq)=0.d0
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
              
      do icz=1,3
      do iv=1,11
       vvx=dble(iv-1)/10.d0
      do ifock=1,nfock
       do iy=2,31
        sy=sgap**2*(spmax/sgap**3)**((iy-1)/30.d0)
        rp=(rq(ifock,icz)+alfp*dlog(sy))*4.d0*.0389d0
        do iz=2,11 
        do iqq=1,5
         qfanu(iy,iz,iv,ifock+nfock*(icz-1),iqq)
     *   =qfanu(iy-1,iz,iv,ifock+nfock*(icz-1),iqq)
         if(iy.lt.31)qfanu(iy+1,iz,iv,ifock+nfock*(icz-1),iqq)
     *   =qfanu(iy,iz,iv,ifock+nfock*(icz-1),iqq)
        enddo
        enddo

        n=1
10      n=n+1
        nrep=0
        do iz=2,11  
         if(iz.gt.6)then
          z=.2d0*dble(iz-6)
          b=dsqrt(-dlog(z)*rp)
         else
          b=dsqrt(-rp*(dlog(0.2d0)+2.5d0*dble(iz-7)))
          z=dexp(-b*b/rp)
         endif
         call qgfan(sy,b*b,vvx,ifock,icz,fann)
         fanmin=qglegi(sy,b*b,ifock,icz,0)*qgfact(sy)
         do iqq=1,4
          fann(iqq)=max(fanmin,fann(iqq))
         enddo
         fann(1)=min(fann(1),fann(2))
         fann(4)=min(fann(4),fann(2))
         fann(3)=min(fann(4),fann(3))
 
         do iqq=1,5
          if(fann(iqq).gt.0.d0)then
           qfan0(iz,iqq)=dlog(fann(iqq)/qglegi(sy,b*b,ifock,icz,0))
          elseif(iy.gt.2)then
           qfan0(iz,iqq)=min(qfanu(iy-1,iz,iv,ifock+nfock*(icz-1),iqq)
     *     ,2.d0*qfanu(iy-1,iz,iv,ifock+nfock*(icz-1),iqq)
     *     -qfanu(iy-2,iz,iv,ifock+nfock*(icz-1),iqq))
          else
           stop'qfanu: iy=2'
          endif
          if(abs(qfan0(iz,iqq)-qfanu(iy,iz,iv,ifock+nfock*(icz-1),iqq))
     *    .gt.1.d-3)nrep=1
         enddo
        enddo
       
        do iz=2,11 
        do iqq=1,5
         qfanu(iy,iz,iv,ifock+nfock*(icz-1),iqq)=qfan0(iz,iqq)
        enddo
        enddo
        if(nrep.eq.1)goto 10
       enddo
      enddo
      enddo
      enddo

c-------------------------------------------------
c b-dependent pdfs (hp-case)
      if(jrep.ne.0)then
       do icz=1,3
       do ifock=1,nfock
        do iv=1,1
         vvx=dble(iv-1)/10.d0
        do iv1=1,1
         vvxp=dble(iv1-1)/5.d0
        do iy=1,31
         xr=(spmax/sgap**2)**(-(iy-1)/30.d0)/sgap
         rp=(rq(ifock,icz)-alfp*dlog(xr))*4.d0*.0389d0
         do iz=1,11
          if(iy.eq.1.or.iz.eq.1)then
           do jj=1,1
           do iqq=1,2
            pdfr(iy,iz,iv,icz+(icz-1)*(3-icz)*(iv1+1)
     *      ,ifock+nfock*(iqq-1)+2*nfock*(jj-1))=0.d0
           enddo
           enddo
          else
           if(iz.gt.6)then
            z=.2d0*dble(iz-6)
            b=dsqrt(-dlog(z)*rp)
           else
            b=dsqrt(-rp*(dlog(0.2d0)+2.5d0*dble(iz-7)))
            z=dexp(-b*b/rp)
           endif
           do jj=1,1 
            do iqq=1,2
             ffmin=max(1.d-5,qgpdfm(xr,ifock,icz,iqq))
             ff=max(ffmin,qgpdfb(xr,b*b,vvx,vvxp,ifock,icz,iqq,jj))
             if(ff.gt.0.d0)then
              pdfr(iy,iz,iv,icz+(icz-1)*(3-icz)*(iv1+1)
     *        ,ifock+nfock*(iqq-1)+2*nfock*(jj-1))=dlog(ff)
             else
              stop'pdfr<0'
             endif
            enddo
           enddo
          endif  
         enddo
        enddo
        enddo
        enddo
       enddo
       enddo
       
c fix gluon PDF by momentum conservation
       do icz=1,3
        do icdp=1,nfock
         iii=0
11       iii=iii+1
         gnorm=0.d0
         seanrm=0.d0
         do i=1,7
         do m=1,2
          xxg=(.5d0+x1(i)*(m-1.5d0))**(1.d0/(1.d0-dels))
          gnorm=gnorm+a1(i)*qgpdff(xxg,icdp,icz,1)*xxg**dels
          seanrm=seanrm+a1(i)*qgpdff(xxg,icdp,icz,2)*xxg**dels
         enddo
         enddo
         gnorm=gnorm/(1.d0-dels)/2.d0
         seanrm=seanrm/(1.d0-dels)/2.d0
         bbbi(icdp,icz)=bbbi(icdp,icz)+(1.d0-gnorm-seanrm)
     *   /4.d0/pi/rr/fp(icdp,icz)*qggamfun(real(3.d0-dels+beth(icz)))
     *   /qggamfun(real(2.d0-dels))/qggamfun(real(1.d0+beth(icz)))
     *   /(1.d0+bbbi(icdp,icz)*(2.d0-dels)/(3.d0-dels+beth(icz))
     *   +.5d0*bbbi(icdp,icz)*(bbbi(icdp,icz)-1.d0)*(2.d0-dels)
     *   *(3.d0-dels)/(3.d0-dels+beth(icz))/(4.d0-dels+beth(icz)))
         bbbi(icdp,icz)=max(-20.d0,bbbi(icdp,icz))
         bbbi(icdp,icz)=min(10.d0,bbbi(icdp,icz))
         if(abs(1.d0-gnorm-seanrm).gt.1.d-4)goto 11
        enddo
       enddo
       jrep=0
       goto 9
      endif

c-------------------------------------------------
c cut fan contributions
      if(debug.ge.1)write (moniou,215)
      do icz=1,3
      do ifock=1,nfock
      do iv=1,11
      do iv1=1,1+5*(icz-1)*(3-icz)
      do iv2=1,1+5*(icz-1)*(3-icz)
      do iy=1,31
      do iz=1,11 
      do iqq=1,13
       qfanc(iy,iz,iv,icz+(icz-1)*(3-icz)*(iv1+1+6*(iv2-1))
     * ,ifock+nfock*(iqq-1))=0.d0
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
     
      do icz=1,3                                  !hadron class
      do ifock=1,nfock
c vvx,vvxp,vvxpl - screening corrections from targ. and nuclear proj. fans
      do iv=1,11  
       vvx=dble(iv-1)/10.d0
      do iv1=1,1+5*(icz-1)*(3-icz)
       vvxp=dble(iv1-1)/5.d0
      do iv2=1,1+5*(icz-1)*(3-icz)
       vvxpl=vvx*dble(iv2-1)/5.d0
       do iy=2,31
        sy=sgap**2*(spmax/sgap**3)**((iy-1)/30.d0)
        rp=(rq(ifock,icz)+alfp*dlog(sy))*4.d0*.0389d0
        do iz=2,11 
        do iqq=1,13
         qfanc(iy,iz,iv,icz+(icz-1)*(3-icz)*(iv1+1+6*(iv2-1))
     *   ,ifock+nfock*(iqq-1))=qfanc(iy-1,iz,iv,icz+(icz-1)*(3-icz)
     *   *(iv1+1+6*(iv2-1)),ifock+nfock*(iqq-1))
         if(iy.lt.31)qfanc(iy+1,iz,iv,icz+(icz-1)*(3-icz)
     *   *(iv1+1+6*(iv2-1)),ifock+nfock*(iqq-1))=qfanc(iy,iz,iv
     *   ,icz+(icz-1)*(3-icz)*(iv1+1+6*(iv2-1)),ifock+nfock*(iqq-1))
        enddo
        enddo

        n=1
12      n=n+1                          !number of t-channel iterations
        nrep=0
        do iz=2,11  
         if(iz.gt.6)then
          z=.2d0*dble(iz-6)
          b=dsqrt(-dlog(z)*rp)
         else
          b=dsqrt(-rp*(dlog(0.2d0)+2.5d0*dble(iz-7)))
          z=dexp(-b*b/rp)
         endif
         call qgfanc(sy,b*b,vvx,vvxp,vvxpl,ifock,icz,fann)
         fanmin=qglegi(sy,b*b,ifock,icz,0)*qgfact(sy)
         do iqq=1,12
          fann(iqq)=max(fanmin,fann(iqq))
         enddo
         fann(2)=min(fann(1),fann(2))
         fann(3)=min(fann(3),fann(2))
         fann(6)=min(fann(5),fann(6))
         fann(7)=min(fann(7),fann(6))
         fann(8)=min(fann(7),fann(8))
         fann(9)=min(fann(9),fann(10))
         fann(11)=min(fann(11),fann(12))

         do iqq=1,13
          if(fann(iqq).gt.0.d0)then
           qfan0(iz,iqq)=dlog(fann(iqq)/qglegi(sy,b*b,ifock,icz,0))
          elseif(iy.gt.2)then
           qfan0(iz,iqq)=min(2.d0*qfanc(iy-1,iz,iv,icz+(icz-1)*(3-icz)
     *     *(iv1+1+6*(iv2-1)),ifock+nfock*(iqq-1))-qfanc(iy-2,iz,iv
     *     ,icz+(icz-1)*(3-icz)*(iv1+1+6*(iv2-1)),ifock+nfock*(iqq-1))
     *     ,qfanc(iy-1,iz,iv,icz+(icz-1)*(3-icz)*(iv1+1+6*(iv2-1))
     *     ,ifock+nfock*(iqq-1)))
          else
           stop'qfanc: iy=2'
          endif
          if(abs(qfan0(iz,iqq)-qfanc(iy,iz,iv,icz+(icz-1)*(3-icz)
     *    *(iv1+1+6*(iv2-1)),ifock+nfock*(iqq-1))).gt.1.d-3)nrep=1
         enddo
        enddo 

        do iz=2,11 
        do iqq=1,13 
         qfanc(iy,iz,iv,icz+(icz-1)*(3-icz)*(iv1+1+6*(iv2-1))
     *   ,ifock+nfock*(iqq-1))=qfan0(iz,iqq)
        enddo
        enddo
        if(nrep.eq.1.and.n.lt.50)goto 12
       enddo
      enddo
      enddo
      enddo
      enddo
      enddo

c-------------------------------------------------
c b-dependent pdfs
      do icz=1,3
      do ifock=1,nfock
       do iv=1,11
        vvx=dble(iv-1)/10.d0
       do iv1=1,1+5*(icz-1)*(3-icz)
        vvxp=dble(iv1-1)/5.d0
       do iy=1,31
        xr=(spmax/sgap**2)**(-(iy-1)/30.d0)/sgap
        rp=(rq(ifock,icz)-alfp*dlog(xr))*4.d0*.0389d0
        do iz=1,11
         if(iy.eq.1.or.iz.eq.1)then
          do jj=1,2
          do iqq=1,2
           pdfr(iy,iz,iv,icz+(icz-1)*(3-icz)*(iv1+1)
     *     ,ifock+nfock*(iqq-1)+2*nfock*(jj-1))=0.d0
          enddo
          enddo
         else
          if(iz.gt.6)then
           z=.2d0*dble(iz-6)
           b=dsqrt(-dlog(z)*rp)
          else
           b=dsqrt(-rp*(dlog(0.2d0)+2.5d0*dble(iz-7)))
           z=dexp(-b*b/rp)
          endif
          do jj=1,2 
           do iqq=1,2
            ffmin=max(1.d-5,qgpdfm(xr,ifock,icz,iqq))
            ff=max(ffmin,qgpdfb(xr,b*b,vvx,vvxp,ifock,icz,iqq,jj))
            if(ff.gt.0.d0)then
             pdfr(iy,iz,iv,icz+(icz-1)*(3-icz)*(iv1+1)
     *       ,ifock+nfock*(iqq-1)+2*nfock*(jj-1))=dlog(ff)
            else
             stop'pdfr<0'
            endif
           enddo
          enddo
         endif  
        enddo
       enddo
       enddo
       enddo
      enddo
      enddo

c-------------------------------------------------
c integrands for HT corrections
      if(debug.ge.1)write (moniou,204)
      s2min=4.d0*fqscal*qt0
      do icz=1,3
      do iy=1,20
       sy=s2min*(spmax/s2min)**(dble(iy)/20.d0)
      do ix=1,10
       xp=(sy/s2min)**((ix-11)/11.d0)
      do iq1=1,2
      do iq2=1,2
       tmin=qt0*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qt0*fqscal/sy/xp)))
       sigg=qgborn(qt0,qt0,xp*sy,iq1-1,iq2-1)*(1.d0/tmin-2.d0/sy/xp)
     * +qgjett(qt0,qt0,xp*sy,iq1-1,iq2-1)
     * +qgjeto(qt0,qt0,xp*sy,iq1-1,iq2-1)
     * +qgjeto(qt0,qt0,xp*sy,iq2-1,iq1-1)

      do ifock=1,nfock
       rp=(rq(ifock,icz)-alfp*dlog(xp))*4.d0*.0389d0
      do iz=1,11 
       if(iz.gt.6)then
        z=.2d0*dble(iz-6)
        bbp=max(0.d0,-dlog(z)*rp)
       else
        bbp=-rp*(dlog(0.2d0)+2.5d0*dble(iz-7))
        z=dexp(-bbp/rp)
       endif
      do iv=1,6
       vvx=dble(iv-1)/5.d0
      do ig=1,11
       if(ig.le.6)then
        genh=1.d0+.1d0*(ig-1)
       else
        genh=1.5d0*exp(dble(ig-6)/2.d0)
       endif
       htfacm=htfac*genh

       do jj=1,2
        if(iz.eq.1)then
         feikht(iy,ix,iz,iv+6*(ig-1),iq1+2*(iq2-1)+4*(ifock-1)
     *   +4*nfock*(icz-1)+12*nfock*(jj-1))=0.d0
        else
         ff=qgfht(sy,xp,bbp,vvx,iq1,iq2,ifock,icz,jj)
         ff1=sigg*qgpdfbi(xp,bbp,vvx,0.d0,ifock,icz,iq1,jj)
         if(ff1.gt.0.d0)then
          feikht(iy,ix,iz,iv+6*(ig-1),iq1+2*(iq2-1)+4*(ifock-1)
     *   +4*nfock*(icz-1)+12*nfock*(jj-1))=dlog(max(1.d-10,ff/ff1))
         else
          feikht(iy,ix,iz,iv+6*(ig-1),iq1+2*(iq2-1)+4*(ifock-1)
     *   +4*nfock*(icz-1)+12*nfock*(jj-1))=0.d0
         endif
        endif
       enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      do iy=1,20
       sy=s2min*(spmax/s2min/sgap)**(dble(iy)/20.d0)
      do ix=1,10
       xp=(sy/s2min)**((ix-11)/11.d0)
      do iv=1,6
       vvxi=dble(iv-1)/5.d0
      do iqp=1,2
      do iqt=1,2
       tmin=qt0*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qt0*fqscal/sy/xp)))
       sigg=qgborn(qt0,qt0,xp*sy,iqp-1,iqt-1)*(1.d0/tmin-2.d0/sy/xp)
     * +qgjett(qt0,qt0,xp*sy,iqp-1,iqt-1)
     * +qgjeto(qt0,qt0,xp*sy,iqp-1,iqt-1)
     * +qgjeto(qt0,qt0,xp*sy,iqt-1,iqp-1)
      do jj=1,2
       do ifact=1,12
        if(ifact.le.7)then
         facht=dble(ifact-2)
        else
         facht=5.d0*exp(dble(ifact-7)/2.d0)
        endif
        if(ifact.eq.1)then
         ff0=sigg*qgloopi(1.d0/xp,vvxi,iqp,jj)                  !no 1/2
         ffhtm(iy,ix,ifact,iv,iqp+2*(iqt-1)+4*(jj-1))=dlog(ff0)
        else
         htfm=qgfhtm(xp*sy,xp,facht,vvxi,iqp,iqt,jj)
         if(ifact.eq.2)ff1=htfm
         ffhtm(iy,ix,ifact,iv,iqp+2*(iqt-1)+4*(jj-1))
     *   =dlog(max(1.d-10,htfm/ff1))
        endif
       enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      do iy=1,20
       sy=s2min*(spmax/s2min/sgap)**(dble(iy)/20.d0)
      do iv=1,6
       vvxi=dble(iv-1)/5.d0
      do iqt=1,2
      do jj=1,2
       do ifact=1,12
        if(ifact.le.7)then
         facht=dble(ifact-2)
        else
         facht=5.d0*exp(dble(ifact-7)/2.d0)
        endif
        if(ifact.eq.1)then
         ff0=qglm0(sy,vvxi,iqt,jj)
         flhtm(iy,ifact,iv,iqt,jj)=dlog(ff0)
        else
         htlm=qglhtm(sy,facht,vvxi,iqt,jj)
         if(ifact.eq.2)ff1=htlm
         flhtm(iy,ifact,iv,iqt,jj)=dlog(max(1.d-10,htlm/ff1))
        endif
       enddo
      enddo
      enddo
      enddo
      enddo
     
c-------------------------------------------------
c total interaction eikonals
      do icz=1,3
      do ifockp=1,nfock
      do ifockt=1,nfock
       iic=ifockp+nfock*(ifockt-1)
       do iy=1,11
        e0n=10.d0**iy
        sy=2.d0*e0n*am(2)+am(2)**2+am(icz)**2
        rp=(rq(ifockp,icz)+rq(ifockt,2)+alfp*log(sy))*4.d0*.0389d0
        do iz=1,11  
         if(iz.eq.1)then
          do iv=1,6
          do iv1=1,1+5*(icz-1)*(3-icz)
          do iv2=1,6
           ivv=iv+6*(iv2-1)+18*(icz-1)*(4-icz+2*(3-icz)*iv1)
           do iqq=1,7
           do ig=1,11
            qpomr(iy,iz,ivv,iic,iqq+7*(ig-1))=0.d0
           enddo
           enddo
           if(icz.eq.2)then
            do iqd=1,5
             dhteik(iy,iz,iv+6*(iv1-1)+36*(iv2-1),iic,iqd)=1.d-20
             do ig=1,11
              dhteik(iy,iz,iv+6*(iv1-1)+36*(iv2-1),iic,iqd+5*ig)=0.d0
             enddo
            enddo
           endif
          enddo
          enddo
          enddo
         else
          if(iz.gt.6)then
           z=.2d0*(iz-6)
           b=sqrt(-log(z)*rp)
          else
           b=sqrt(-rp*(log(0.2d0)+2.5d0*(iz-7)))
           z=exp(-b*b/rp)
          endif
          vsoft=qgpomi(sy,b*b,0.d0,0.d0,0.d0,1.d0,1.d0
     *    ,ifockp,ifockt,icz,0)
 
          vgg=qgfsh(sy,b*b,ifockp,ifockt,icz,0)
          vqg=qgfsh(sy,b*b,ifockp,ifockt,icz,1)
          vgq=qgfsh(sy,b*b,ifockp,ifockt,icz,2)
          vqq=qghard(sy,b*b,ifockp,ifockt,icz)
          qxp=vsoft+vgg
          qxpq=vqg+vgq+vqq
    
          do iv=1,6
           vvx=(iv-1)/5.d0
          do iv1=1,1+5*(icz-1)*(3-icz)
           vvxp=(iv1-1)/5.d0
          do iv2=1,6
           vvxt=(iv2-1)/5.d0
           ivv=iv+6*(iv2-1)+18*(icz-1)*(4-icz+2*(3-icz)*iv1)

           v3p=qg3pom(sy,b,vvx,vvxp,vvxt,ifockp,ifockt,icz,1)
           v3pq=qg3pom(sy,b,vvx,vvxp,vvxt,ifockp,ifockt,icz,2)
           call qgpcut(sy,b,vvx,vvxp,vvxt,ifockp,ifockt,icz,v1p)
            
           do ig=1,11           !nuclear (target) enhancement of HT corrections
            if(ig.le.6)then
             genh=1.d0+.1d0*(ig-1)
            else
             genh=1.5d0*exp(dble(ig-6)/2.d0)
            endif
            call qg3pht(sig3ht,sy,b,vvx,vvxp,vvxt,1.d0
     *      ,ifockp,ifockt,icz,2)
            do i=1,2                    !single hard (1 - uncut, 2 - cut)
             vht(i)=sig3ht(2*i)-sig3ht(2*i-1)
             vhtqt(i)=sig3ht(2*i+4)-sig3ht(2*i+3)
            enddo
           
            call qghtm(sig3htm,sy,b,vvx,vvxp,vvxt,1.d0
     *      ,ifockp,ifockt,icz,2)
            vht(3)=sig3htm(1)           !mult-hard (gg-uncut)
            vht(4)=sig3htm(3)           !mult-hard (gg-cut)
            vht(5)=sig3htm(5)           !mult-hard (soft-cut)
            vhtqt(3)=sig3htm(2)         !mult-hard (qv-uncut)
            vhtqt(4)=sig3htm(4)         !mult-hard (qv-cut)

            call qg3pht(sig3ht,sy,b,vvx,vvxt,vvxp,genh
     *      ,ifockt,ifockp,2,icz)
            do i=1,2
             if(sig3ht(2*i-1).ne.0.d0)vht(i)=vht(i)
     *       *sig3ht(2*i)/sig3ht(2*i-1)+sig3ht(2*i)-sig3ht(2*i-1)
             vhtqp(i)=sig3ht(2*i+4)-sig3ht(2*i+3)
            enddo
            call qghtm(sig3htm,sy,b,vvx,vvxt,vvxp,genh
     *      ,ifockt,ifockp,2,icz)
            vht(3)=vht(3)+sig3htm(1)
            vht(4)=vht(4)+sig3htm(3)
            vht(5)=vht(5)+sig3htm(5)
            vhtqp(3)=sig3htm(2)
            vhtqp(4)=sig3htm(4)
            
            if(icz.eq.2)then
             do iqd=1,5
              if(iqd.eq.1)then           !gg-uncut
               vhtn=sig3ht(1)-sig3ht(2)-sig3htm(1)
              elseif(iqd.eq.2)then       !q_v-uncut
               vhtn=sig3ht(5)-sig3ht(6)-sig3htm(2)
              elseif(iqd.eq.3)then       !gg-cut
               vhtn=sig3ht(3)-sig3ht(4)-sig3htm(3)
              elseif(iqd.eq.4)then       !soft-cut
               vhtn=sig3htm(5)
              elseif(iqd.eq.5)then       !q_v-cut
               vhtn=sig3ht(7)-sig3ht(8)-sig3htm(4)
              endif
              if(ig.eq.1)then
               vht0(iqd)=vhtn
               dhteik(iy,iz,iv+6*(iv2-1)+36*(iv1-1)
     *         ,ifockt+nfock*(ifockp-1),iqd)=dlog(max(1.d-20,vhtn))
              else
               if(vht0(iqd).gt.0.d0)then
                dhteik(iy,iz,iv+6*(iv2-1)+36*(iv1-1)
     *         ,ifockt+nfock*(ifockp-1),iqd+5*ig)=dlog(vhtn/vht0(iqd))
               else
                dhteik(iy,iz,iv+6*(iv2-1)+36*(iv1-1)
     *         ,ifockt+nfock*(ifockp-1),iqd+5*ig)=0.d0
               endif
              endif
             enddo
            endif
           
            do iqq=1,7
             if(iqq.eq.1)then            !(g+soft)-uncut
              qxpt=qxp+v3p+vht(1)+vht(3)
             elseif(iqq.eq.2)then        !q_v-uncut
              qxpt=qxpq+v3pq+vsoft+vhtqp(1)+vhtqt(1)+vhtqp(3)+vhtqt(3)
             elseif(iqq.eq.3)then        !(g+soft)-cut
              qxpt=qxp+v1p(1)+vht(2)+vht(4)+vht(5)
             elseif(iqq.eq.4)then        !soft-cut
              qxpt=vsoft+v1p(2)+vht(5)
             elseif(iqq.eq.5)then        !q_v-cut
              qxpt=qxpq+v1p(3)+vsoft+vhtqp(2)+vhtqt(2)
     *        +vhtqp(4)+vhtqt(4)
             elseif(iqq.eq.6)then        !q_v_p-cut
              qxpt(iqq)=vqg+v1p(4)+vsoft+vhtqp(2)+vhtqp(4)
             elseif(iqq.eq.7)then        !q_v_t-cut
              qxpt(iqq)=vgq+v1p(5)+vsoft+vhtqt(2)+vhtqt(4)
             else
              stop'pomr: iqq?!'
             endif
             
             if(qxpt(iqq).gt.0.d0)then
              qpomr(iy,iz,ivv,iic,iqq+7*(ig-1))=log(qxpt(iqq)/vsoft)
             elseif(iy.gt.2)then
              qpomr(iy,iz,ivv,iic,iqq+7*(ig-1))
     *        =min(qpomr(iy-1,iz,ivv,iic,iqq+7*(ig-1))
     *        ,2.d0*qpomr(iy-1,iz,ivv,iic,iqq+7*(ig-1))
     *        -qpomr(iy-2,iz,ivv,iic,iqq+7*(ig-1)))
             else
              qpomr(iy,iz,ivv,iic,iqq+7*(ig-1))
     *        =qpomr(iy-1,iz,ivv,iic,iqq+7*(ig-1))
             endif
            enddo
           enddo
          enddo
          enddo
          enddo
         endif
        enddo
       enddo
      enddo
      enddo
      enddo

c-------------------------------------------------
c interaction cross sections
      if(debug.ge.1)write (moniou,203)
      ia(1)=1
      do iy=1,10
       e0n=10.d0**iy                     !interaction energy
       do iiz=1,3  
        icz=iiz                          !hadron class
        scm=2.d0*e0n*am(2)+am(2)**2+am(icz)**2
        do iia=1,6
         if(iia.le.4)then
          ia(2)=4**(iia-1)               !target mass number
         elseif(iia.eq.5)then
          ia(2)=14
         else
          ia(2)=40
         endif
         if(debug.ge.1)write (moniou,206)e0n,ty(icz),ia(2)

         if(ia(2).lt.10)then             !gaussian for light nucleus density
          rnuc(2)=trnuc(ia(2))           !nuclear radius (from data tables)
          wsnuc(2)=amws                  !not used
          wbnuc(2)=0.d0                  !not used
         elseif(ia(2).le.56)then   !3-parameter Fermi distr. for heavier nuclei
          rnuc(2)=trnuc(ia(2))           !nuclear radius
          wsnuc(2)=twsnuc(ia(2))         !diffuseness
          wbnuc(2)=twbnuc(ia(2))         !'wine-bottle' parameter
         else                            !2-parameter Fermi for heaviest
          rnuc(2)=1.19d0*dble(ia(2))**(1.d0/3.d0)
     *    -1.38d0*dble(ia(2))**(-1.d0/3.d0) !nuclear radius
          wsnuc(2)=amws                  !diffuseness
          wbnuc(2)=0.d0                  !not used
         endif

         if(ia(2).eq.1)then
          call qgfz(0.d0,gz0,0,0)        !hadron-proton cross sections (mb)
          gtot=gz0(1)                    !total cross-section
          gin=(gz0(2)+gz0(3)+gz0(4))/2.d0!inelastic cross section
          gel=gtot-gin                   !elastic cross section
          gdp=gz0(3)/2.d0                !proj. low mass diffr. (+double LMD)
          gdt=gz0(4)/2.d0                !target low mass diffraction
          bel=gz0(5)                     !elastic scattering slope
          if(iy.le.10)gsect(iy,icz,iia)=log(gin)
          if(debug.ge.1)write (moniou,225)gtot,gin,gel,gdp,gdt,bel
         else
          call qggau(gz1)                !hadron-nucleus cross sections (fm^2)
          gin=gz1(1)*10.d0               !inelastic cross section (mb)
          if(iy.le.10)gsect(iy,icz,iia)=log(gin)
          if(debug.ge.1)write (moniou,224)gin,gz1(3)*10.d0,gz1(2)*10.d0
         endif
         if(.not.(gsect(iy,icz,iia).le.0.d0
     *   .or.gsect(iy,icz,iia).gt.0.d0))stop'qpomr-nan'
        enddo
       enddo
      enddo
         
c-----------------------------------------------------------------------------
c timelike Sudakov formfactor
      if(debug.ge.1)write (moniou,221)
      do m=1,2                     !parton type (1-g, 2-q)
       fsud(1,m)=0.d0
      do k=2,10
       qmax=qtf*4.d0**(1.d0+k)     !effective virtuality (qt**2/z**2/(1-z)**2)
       fsud(k,m)=qgsudt(qmax,m)
      enddo
      enddo
      if(debug.ge.1)write (moniou,222)
      do m=1,2                     !parton type (1-g, 2-q)
      do k=1,10
       qlmax=1.38629d0*(k-1)
       qrt(k,1,m)=0.d0
       qrt(k,101,m)=qlmax
      do i=1,99                    !bins in Sudakov formfactor
       if(k.eq.1)then
        qrt(k,i+1,m)=0.d0
       else
        qrt(k,i+1,m)=qgroot(qlmax,.01d0*i,m)
       endif
      enddo
      enddo
      enddo
     

c-----------------------------------------------------------------------------
c writing cross sections to the file
      if(debug.ge.1)write (moniou,220)
      if(ifIIIdat.ne.1)then
       open(1,file=DATDIR(1:INDEX(DATDIR,' ')-1)//'qgsjetIII.dat'
     * ,status='unknown')
      else                                              !used to link with nexus
       open(ifIIIdat,file=fnIIIdat(1:nfnIIIdat),status='unknown')
      endif
      write (1,*)csborn,cs0,cstot,evk,qpomi,qpomis,qloopr,qlegi,qfanu
     * ,qfanc,pdfr,qpomr,dhteik,feikht,ffhtm,flhtm,gsect,fsud,qrt
      close(1)

c-----------------------------------------------------------------------------
c nuclear cross sections
15    if(ifIIIncs.ne.2)then
       inquire(file=DATDIR(1:INDEX(DATDIR,' ')-1)//'sectnu-III'
     * ,exist=lcalc)
      else                                                  !ctp
       inquire(file=fnIIIncs(1:nfnIIIncs),exist=lcalc)
      endif
      if(lcalc)then
       if(debug.ge.0)write (moniou,207)
       if(ifIIIncs.ne.2)then
        open(2,file=DATDIR(1:INDEX(DATDIR,' ')-1)//'sectnu-III'
     *  ,status='old')
       else                                                  !ctp
        open(ifIIIncs,file=fnIIIncs(1:nfnIIIncs),status='old')
       endif
       read (2,*)qgsasect
       close(2)

      elseif(.not.producetables)then
        write(moniou,*) "Missing sectnu-III file !"
        write(moniou,*) "Please correct the defined path ",
     &"or force production ..."
        stop

      else
       niter=5000                   !number of iterations
       do ie=1,10
        e0n=10.d0**ie               !interaction energy (per nucleon)
       do iia1=1,6
        iap=2**iia1                 !proj. mass number
       do iia2=1,6
        if(iia2.le.4)then
         iat=4**(iia2-1)            !targ. mass number
        elseif(iia2.eq.5)then
         iat=14
        else
         iat=40
        endif
        if(debug.ge.1)write (moniou,208)e0n,iap,iat

        call qgini(e0n,2,iap,iat)
        call qgcrossc(niter,gtot,gprod,gabs,gdd,gqel,gcoh)
        qgsasect(ie,iia1,iia2)=log(gprod)
        if(debug.ge.1)write (moniou,209)gtot,gprod,gabs,gdd,gqel,gcoh
       enddo
       enddo
       enddo
       if(ifIIIncs.ne.2)then
        open(2,file=DATDIR(1:INDEX(DATDIR,' ')-1)//'sectnu-III'
     *  ,status='unknown')
       else                                                  !ctp
        open(ifIIIncs,file=fnIIIncs(1:nfnIIIncs),status='unknown')
       endif
       write (2,*)qgsasect
       close(2)
      endif
 
      if(debug.ge.3)write (moniou,218)
201   format(2x,'qgaini: hard cross sections calculation')
202   format(2x,'qgaini: hard cross sections read out')
203   format(2x,'qgaini: hadron-proton/nucleus cross sections')
204   format(2x,'qgaini: integrands for HT corrections')
205   format(2x,'qgaini: pretabulation of the interaction eikonals')
206   format(2x,'qgaini: initial particle energy:',e10.3,2x
     *,'its type:',a7,2x,'target mass number:',i2)
207   format(2x,'qgaini: nuclear cross sections readout from the file'
     *,' sectnu-III')
208   format(2x,'qgaini: initial nucleus energy:',e10.3,2x
     *,'projectile mass:',i2,2x,'target mass:',i2)
209   format(2x,'gtot',d10.3,'  gprod',d10.3,'  gabs',d10.3
     */2x,'gdd',d10.3,'  gqel',d10.3,' gcoh',d10.3)
210   format(2x,'qgaini - main initialization procedure')
211   format(2x,'qgaini: Pomeron loop contributions')
212   format(2x,'qgaini: integrated Pomeron leg eikonals')
213   format(2x,'qgaini: integrated fan contributions')
214   format(2x,'qgaini: cross sections readout from the file: ', A,2x)
c     *,' qgsjetIII.dat')
215   format(2x,'qgaini: integrated cut fan contributions')
c216   format(2x,'qgaini: integrated cut Pomeron eikonals')
218   format(2x,'qgaini - end')
220   format(2x,'qgaini: cross sections are written to the file'
     *,' qgsjetIII.dat')
221   format(2x,'qgaini: timelike Sudakov formfactor')
222   format(2x,'qgaini: effective virtuality for inversion')
224   format(2x,'qgaini: hadron-nucleus cross sections:'
     */4x,'gin=',e10.3,2x,'gdifr_targ=',e10.3,2x
     *,'gdifr_proj=',e10.3)
225   format(2x,'qgaini: hadron-proton cross sections:'
     */4x,'gtot=',e10.3,2x,'gin=',e10.3,2x,'gel=',e10.3/4x
     *,'gdifrp=',e10.3,2x,'gdifrt=',e10.3,2x,'b_el=',e10.3,2x)
c226   format(2x,'qgaini: ht-corrections to cut fans')
      return
      end

c=============================================================================
      subroutine qgini(e0n,icp0,iap,iat)
c-----------------------------------------------------------------------------
c additional initialization procedure
c e0n  - interaction energy (per hadron/nucleon),
c icp0 - hadron type (+-1 - pi+-, +-2 - p(p~), +-3 - n(n~),
c                     +-4 - K+-, +-5 - K_l/s),
c iap  - projectile mass number (1 - for a hadron),
c iat  - target mass number
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208,nfock=3)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr4/  ey0(3)
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),bcoll
      common /qgarr10/ am(6)
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /arr1/    trnuc(56),twsnuc(56),twbnuc(56)
      common /qgdebug/ debug
      common /qgsIIInex1/xan(iapmax,3),xbn(iapmax,3)  !used to link with nexus
     *,bqgs,bmaxqgs,bmaxnex,bminnex

      if(debug.ge.1)write (moniou,201)icp0,iap,iat,e0n
      icp=icp0
      ia(1)=iap
      ia(2)=iat

      icz=iabs(icp)/2+1  !particle class (1 - pion, 2 - nucleon, 3 - kaon)

      scm=2.d0*e0n*am(2)+am(2)**2+am(icz)**2   !c.m. energy squared
      wplab=e0n+dsqrt(e0n**2-am(icz)**2)       !proj. LC momentum
      ey0(1)=dsqrt(scm)/(am(2)+wplab)          !Lorentz boost to lab. frame
      ey0(2)=1.d0
      ey0(3)=1.d0

c-------------------------------------------------
c nuclear radii and weights for nuclear configurations - procedure qggea
      do i=1,2
       if(ia(i).lt.10.and.ia(i).ne.1)then      !gaussian density
        rnuc(i)=trnuc(ia(i))                   !nuclear radius
        rnuc(i)=rnuc(i)*dsqrt(2.d0*ia(i)/(ia(i)-1.d0)) !rnuc -> rnuc*a/(a-1)
       elseif(ia(i).ne.1)then
        if(ia(i).le.56)then                    !3-parameter Fermi
         rnuc(i)=trnuc(ia(i))                  !nuclear radius
         wsnuc(i)=twsnuc(ia(i))                !diffuseness
         wbnuc(i)=twbnuc(ia(i))                !'wine-bottle' parameter
        else                                   !2-parameter Fermi
         rnuc(i)=1.19*float(ia(i))**(1./3.)-1.38*float(ia(i))**(-1./3.)
         wsnuc(i)=amws                         !diffuseness
         wbnuc(i)=0.d0                         !not used
        endif
        cr1(i)=1.d0+3.d0/rnuc(i)*wsnuc(i)+6.d0/(rnuc(i)/wsnuc(i))**2
     *  +6.d0/(rnuc(i)/wsnuc(i))**3
        cr2(i)=3.d0/rnuc(i)*wsnuc(i)
        cr3(i)=3.d0/rnuc(i)*wsnuc(i)+6.d0/(rnuc(i)/wsnuc(i))**2
       endif
      enddo
      if(ia(1).ne.1.and.ia(2).ne.1)then        !primary nucleus
       bm=rnuc(1)+rnuc(2)+3.d0*max(wsnuc(1),wsnuc(2))
     & +dsqrt(alfp*log(scm))                   !b-cutoff
      elseif(ia(2).ne.1)then                   !hadron-nucleus
       bm=rnuc(2)+3.d0*wsnuc(2)+dsqrt(rq(1,icz)+alfp*log(scm)) !b-cutoff
      elseif(ia(1).ne.1)then                   !nucleus-proton
       bm=rnuc(1)+3.d0*wsnuc(1)+dsqrt(rq(1,2)+alfp*log(scm))   !b-cutoff
      else                                     !hadron-proton
       bm=dsqrt(rq(1,icz)+rq(1,2)+alfp*log(scm))               !b-cutoff
      endif

      bmaxqgs=bm                                      !used to link with nexus

      if(debug.ge.3)write (moniou,202)
201   format(2x,'qgini - miniinitialization: particle type icp0=',
     *i2,2x,'projectile mass number iap=',i2/4x,
     *'target mass number iat=',i2,' interaction energy e0n=',e10.3)
202   format(2x,'qgini - end')
      return
      end

c=============================================================================
      double precision function qgpint(sy,bb)
c-----------------------------------------------------------------------------
c qgpint - interm. Pomeron eikonal (semihard)
c sy - pomeron mass squared,
c bb - impact parameter squared
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)sy,bb

      qgpint=0.d0      
      s2min=4.d0*fqscal*qt0
      xmin=s2min/sy
      if(xmin.ge.1.d0)return
      
      xmin=xmin**(delh-dels)
c numerical integration over z1
      do i=1,7
      do m=1,2
       z1=(.5d0*(1.d0+xmin-(2*m-3)*x1(i)*(1.d0-xmin)))
     * **(1.d0/(delh-dels))
       ww=z1*sy
       sjqq=qgjit(qt0,qt0,ww,2,2)  !inclusive qq cross-section
       sjqg=qgjit(qt0,qt0,ww,1,2)  !inclusive qg cross-section
       sjgg=qgjit(qt0,qt0,ww,1,1)  !inclusive gg cross-section
       sjqqq=qgjit(qt0,qt0,ww,3,2) !inclusive qq cross-section (same flavor)
       sjqqa=qgjit(qt0,qt0,ww,4,2) !inclusive qbar-q cross-section (same flavor)
        
       st2=0.d0
       do j=1,7
       do k=1,2
        xx=.5d0*(1.d0+x1(j)*(2*k-3))
        xp=z1**xx
        xm=z1/xp
        glu1=qgppdi(xp,0)
        sea1=qgppdi(xp,1)
        glu2=qgppdi(xm,0)
        sea2=qgppdi(xm,1)
        st2=st2+a1(j)*(glu1*glu2*sjgg+(glu1*sea2+glu2*sea1)*sjqg
     *  +sea1*sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0))
       enddo
       enddo
       rh=-alfp*dlog(z1)
       qgpint=qgpint-a1(i)*dlog(z1)/z1**delh*st2
     * *exp(-bb/(4.d0*.0389d0*rh))/rh      
      enddo
      enddo
      qgpint=qgpint*(1.d0-xmin)/(delh-dels)*factk*rr**2*g3p**2/2.d0*pi
      
      if(debug.ge.3)write (moniou,202)qgpint
201   format(2x,'qgpint - interm. Pomeron eikonal:'
     */4x,'sy=',e10.3,2x,'bb=',e10.3)
202   format(2x,'qgpint=',e10.3)
      return
      end
    
c------------------------------------------------------------------------
      double precision function qgpini(sy,bb,vvxp,vvxt,iqq)
c-----------------------------------------------------------------------
c qgpini - intermediate Pomeron eikonal
c sy   - pomeron mass squared,
c bb   - impact parameter squared,
c vvxp - projectile screening factor,
c vvxt - target screening factor
c iqq=0  - single soft Pomeron
c iqq=1  - single Pomeron (s+g)
c iqq=2  - general loop contribution
c iqq=3  - single Pomeron end on one side
c iqq=4  - single Pomeron ends on both sides
c iqq=5  - general soft loop 
c iqq=6  - soft loop with single Pomeron end on one side
c iqq=7  - soft loop with single Pomeron ends on both sides
c iqq=8  - single cut Pomeron (s+g)
c iqq=9  - single cut Pomeron with single end
c iqq=10 - single cut Pomeron with 2 single ends
c iqq=11 - single cut soft Pomeron
c iqq=12 - single cut soft Pomeron with single end
c iqq=13 - single cut soft Pomeron with 2 single ends
c iqq=14 - any cuts except the complete rap-gap
c iqq=15 - no rap-gap at one side
c iqq=16 - no rap-gap at one side and single Pomeron on the other
c iqq=17 - single cut Pomeron end at one side
c iqq=18 - single cut Pomeron end at one side and single Pomeron on the other
c  with proj/targ screening corrections:
c iqq=19 - single cut Pomeron
c iqq=20 - diffractive cut, Puu
c iqq=21 - diffractive cut, Puu-Puc
c iqq=22 - diffractive cut, Pcc
c iqq=23 - diffractive cut, Pcc+Pcu
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wz(3),wi(3),wj(3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr20/ spmax
      common /qgarr39/ qpomi(31,11,18),qpomis(31,11,11,11,5)
      common /qgarr43/ moniou
      common /qgdebug/ debug
      
      qgpini=0.d0
      rp=alfp*dlog(sy)*4.d0*.0389d0
      z=exp(-bb/rp)
      if(sy.le.sgap**2.or.iqq.eq.0)goto 1

      yl=log(sy/sgap**2)/log(spmax/sgap**4)*30.d0+1.d0
      k=max(1,int(yl-1.d0))
      k=min(k,29)
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)

      z2=.2d0*exp(-12.5d0)
      if(z.gt.z2.or.vvxp+vvxt.ne.0.d0)then
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=(-bb/rp-dlog(0.2d0))/2.5d0+7.d0
       endif
       jz=min(9,int(zz))
       jz=max(2,jz)
       if(jz.eq.6)jz=5
       wz(2)=zz-jz
       wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
       wz(1)=1.d0-wz(2)+wz(3)
       wz(2)=wz(2)-2.d0*wz(3)
      else
       z3=z2*exp(2.5d0)
       jz=1
       wz(1)=(z-z2)*(z-z3)/z2/z3
       wz(2)=z*(z-z3)/z2/(z2-z3)
       wz(3)=z*(z-z2)/z3/(z3-z2)
      endif

      if(iqq.le.18)then
       do l1=1,3
        l2=jz+l1-1
       do k1=1,3
        k2=k+k1-1
        qgpini=qgpini+qpomi(k2,l2,iqq)*wk(k1)*wz(l1)
       enddo
       enddo
       
      else
       vi=vvxp*10.d0+1.d0
       i=max(1,int(vi))
       i=min(i,9)
       wi(2)=vi-i
       wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
       wi(1)=1.d0-wi(2)+wi(3)
       wi(2)=wi(2)-2.d0*wi(3)
       
       vj=vvxt*10.d0+1.d0
       j=max(1,int(vj))
       j=min(j,9)
       wj(2)=vj-j
       wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
       wj(1)=1.d0-wj(2)+wj(3)
       wj(2)=wj(2)-2.d0*wj(3)
       jmax=3

       do j1=1,jmax
        j2=j+j1-1
       do i1=1,3
        i2=i+i1-1
       do l1=1,3
        l2=jz+l1-1
       do k1=1,3
        k2=k+k1-1
        qgpini=qgpini+qpomis(k2,l2,i2,j2,iqq-18)
     *  *wk(k1)*wz(l1)*wi(i1)*wj(j1)
       enddo
       enddo
       enddo
       enddo
      endif
      if(z.lt.z2/2.d0.and.vvxp+vvxt.ne.0.d0)qgpini=min(0.d0,qgpini)

1     qgpini=exp(qgpini)*sy**dels*g3p**2*z/rp*4.d0*.0389d0
      return 
      end

c------------------------------------------------------------------------
      double precision function qglegi(sy,bb,icdp,icz,iqq)
c-----------------------------------------------------------------------
c qglegi - integrated Pomeron leg eikonal
c sy   - pomeron mass squared,
c bb   - impact parameter squared,
c icdp - diffractive state for the hadron,
c icz  - hadron class
c iqq=0  - single soft leg
c iqq=1  - single (s+g)-leg Pomeron
c iqq=2  - single q-leg Pomeron
c iqq=3  - all loops
c iqq=4  - all loops with single Pomeron end
c iqq=5  - all soft loops
c iqq=6  - all soft loops with single Pomeron end
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wz(3)
      parameter(nfock=3)
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr20/ spmax
      common /qgarr27/ qlegi(31,11,nfock,3,6),qfanu(31,11,11,3*nfock,5)
     *,qfanc(31,11,11,39,13*nfock)
      common /qgarr43/ moniou
      common /qgdebug/ debug
      
      qglegi=0.d0
      rp=(rq(icdp,icz)+alfp*log(sy))*4.d0*.0389d0
      z=exp(-bb/rp)
      if(sy.le.sgap**2.or.iqq.eq.0)goto 1

      yl=log(sy/sgap**2)/log(spmax/sgap**3)*30.d0+1.d0
      k=max(1,int(yl))
      k=min(k,29)     
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)
      
      z2=.2d0*exp(-12.5d0)
      if(z.gt.z2)then
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=(-bb/rp-dlog(0.2d0))/2.5d0+7.d0
       endif
       jz=min(9,int(zz))
       jz=max(2,jz)
       if(jz.eq.6)jz=5
       wz(2)=zz-jz
       wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
       wz(1)=1.d0-wz(2)+wz(3)
       wz(2)=wz(2)-2.d0*wz(3)
      else
       z3=z2*exp(2.5d0)
       jz=1
       wz(1)=(z-z2)*(z-z3)/z2/z3
       wz(2)=z*(z-z3)/z2/(z2-z3)
       wz(3)=z*(z-z2)/z3/(z3-z2)
      endif
      
      do l1=1,3
       l2=jz+l1-1
      do k1=1,3
       k2=k+k1-1
       qglegi=qglegi+qlegi(k2,l2,icdp,icz,iqq)*wk(k1)*wz(l1)
      enddo
      enddo
1     qglegi=exp(qglegi)
      if(iqq.eq.2)qglegi=max(0.d0,qglegi-1.d0)
      qglegi=qglegi*z*sy**dels*fp(icdp,icz)*g3p/rp*4.d0*.0389d0
      return 
      end
      
c=============================================================================
      double precision function qglsh(sy,bb,icdp,icz,iqq)
c-----------------------------------------------------------------------------
c qglsh - unintegrated Pomeron leg eikonal
c sy  - pomeron mass squared,
c xp  - light cone momentum share,
c bb  - impact parameter squared,
c icz - hadron class
c iqq=0 - gluon/sea quark contribution,
c iqq=1 - valence quark contribution
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)sy,bb,icz

      qglsh=0.d0
      s2min=4.d0*fqscal*qt0
      if(sy.lt.1.001d0*s2min)return

      xmin=(s2min/sy)**(delh-dels)
      alh=.5d0+dels
c numerical integration over zh
      do i1=1,7
      do m1=1,2
       zh=(.5d0*(1.d0+xmin-(2*m1-3)*x1(i1)*(1.d0-xmin)))
     * **(1.d0/(delh-dels))
       ww=zh*sy         !c.m. energy squared for hard interaction
       sjqq=qgjit(qt0,qt0,ww,2,2)
       sjqg=qgjit(qt0,qt0,ww,1,2)
       sjgg=qgjit(qt0,qt0,ww,1,1)
       sjqqq=qgjit(qt0,qt0,ww,3,2) !inclusive qq cross-section (same flavor)
       sjqqa=qgjit(qt0,qt0,ww,4,2) !inclusive qbar-q cross-section (same flavor)
       if(debug.ge.3)write (moniou,203)ww,sjgg
     * ,sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0
        
       if(iqq.eq.0)then
        stg=0.d0
        do i2=1,7
        do m2=1,2
         xx=.5d0*(1.d0+x1(i2)*(2*m2-3))
         xph=zh**xx
         xmh=zh/xph
         glu1=qgppdc(xph,0,icdp,icz)
         sea1=qgppdc(xph,1,icdp,icz)
         glu2=qgppdi(xmh,0)
         sea2=qgppdi(xmh,1)
         rh=rq(icdp,icz)-alfp*dlog(zh)
         stsum=(glu1*glu2*sjgg+(glu1*sea2+glu2*sea1)*sjqg
     *   +sea1*sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0))
     *   *exp(-bb/(4.d0*.0389d0*rh))/rh
         stg=stg+a1(i2)*stsum
        enddo
        enddo
        qglsh=qglsh-a1(i1)*dlog(zh)/zh**delh*stg

       elseif(iqq.eq.1)then
        stq=0.d0
        xam=zh**alh
        do i2=1,7
        do m2=1,2
         xph=(.5d0*(1.d0+xam+x1(i2)*(2*m2-3)*(1.d0-xam)))**(1.d0/alh)
         xmh=zh/xph
         glu2=qgppdi(xmh,0)
         sea2=qgppdi(xmh,1)
         rh=rq(icdp,icz)-alfp*dlog(xmh)
         stq=stq+a1(i2)*qgvpdf(xph,icz)/dsqrt(xph)
     *   *(glu2*sjqg+sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0))
     *   *exp(-bb/(4.d0*.0389d0*rh))/rh
        enddo
        enddo
        qglsh=qglsh+a1(i1)/zh**delh*stq*(1.d0-xam)
       endif
      enddo
      enddo
      if(iqq.eq.0)then
       qglsh=qglsh*rr**2*(1.d0-xmin)/(delh-dels)*fp(icdp,icz)*g3p*factk
     * /2.d0*pi
      elseif(iqq.eq.1)then
       qglsh=qglsh*rr*g3p*(1.d0-xmin)/(delh-dels)*factk/alh/8.d0
      endif
     
      if(debug.ge.3)write (moniou,202)qglsh
201   format(2x,'qglsh - unintegrated Pomeron leg eikonal:'
     */4x,'s=',e10.3,2x,'b^2=',e10.3,2x,'icz=',i1)
202   format(2x,'qglsh=',e10.3)
203   format(2x,'qglsh:',2x,'s_hard(gg)=',e10.3
     *,2x,'sigma_hard(qq)=',e10.3)
      return
      end    
  
c------------------------------------------------------------------------
      subroutine qgloop(sy,bb,fann,jj)
c-----------------------------------------------------------------------
c qgloop - intermediate Pomeron eikonal with loops
c sy   - pomeron mass squared,
c bb   - impact parameter squared,
c jj=1 - uncut loops (iqq=1,...6)
c jj=2 - cut loops (iqq=7,...17)
c iqq=1  - general loop contribution
c iqq=2  - single Pomeron end on one side
c iqq=3  - single Pomeron ends on both sides
c iqq=4  - general soft loop 
c iqq=5  - soft loop with single Pomeron end on one side
c iqq=6  - soft loop with single Pomeron ends on both sides
c iqq=7  - single cut Pomeron
c iqq=8  - single cut Pomeron with single end
c iqq=9  - single cut Pomeron with 2 single ends
c iqq=10 - single cut soft Pomeron
c iqq=11 - single cut soft Pomeron with single end
c iqq=12 - single cut soft Pomeron with 2 single ends
c iqq=13 - any cuts except the complete rap-gap
c iqq=14 - no rap-gap at one side
c iqq=15 - no rap-gap at one side and single Pomeron on the other
c iqq=16 - single cut Pomeron at one side
c iqq=17 - single cut Pomeron at one side and single Pomeron on the other
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension fann(17)
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)
      
      do iqq=1,17
       fann(iqq)=0.d0
      enddo
      if(sy.le.sgap**2)goto 1

      do ix1=1,7
      do mx1=1,2
       xpomr=(sy/sgap**2)**(-.5d0-x1(ix1)*(mx1-1.5d0))/sgap
       rp=-alfp*log(xpomr)*4.d0*.0389d0
       rp1=alfp*log(xpomr*sy)*4.d0*.0389d0
       rp2=rp*rp1/(rp+rp1)
      do ix2=1,7
      do mx2=1,2
       z=.5d0+x1(ix2)*(mx2-1.5d0)
       bb0=-rp2*log(z)
      do ix3=1,7
      do mx3=1,2
       phi=pi*(.5d0+x1(ix3)*(mx3-1.5d0))
       bb1=(dsqrt(bb)*rp1/(rp+rp1)-dsqrt(bb0)*cos(phi))**2
     * +bb0*sin(phi)**2
       bb2=(dsqrt(bb)*rp/(rp+rp1)+dsqrt(bb0)*cos(phi))**2
     * +bb0*sin(phi)**2
       
       vis=qgpini(xpomr*sy,bb1,0.d0,0.d0,0)           !single soft pomeron
       vist=min(vis,qgpini(xpomr*sy,bb1,0.d0,0.d0,5)) !soft loop sequence
       vi=max(vis,qgpini(xpomr*sy,bb1,0.d0,0.d0,1))   !single pomeron (s+g)
       vit=min(vi,qgpini(xpomr*sy,bb1,0.d0,0.d0,2))   !loop sequence         

       vp0=qgpini(1.d0/xpomr,bb2,0.d0,0.d0,4)        !loop sequence (2 1P-ends)
       vp1=min(vp0,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,3))  !loop sequence (1P-end)
       vp=min(vp1,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,2))   !loop sequence 

       vps0=qgpini(1.d0/xpomr,bb2,0.d0,0.d0,7)          !soft loops (2 1P-ends)
       vps1=min(vps0,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,6))  !soft loops (1P-end)
       vps=min(vps1,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,5))   !soft loop sequence

       if(jj.eq.1)then
        do iqq=1,6
         if(iqq.eq.1)then     !general loop sequence 
          dpx=vi*(min(0.d0,1.d0-exp(-vp)-vp)+vp-vp1)
     *    +min(0.d0,1.d0-exp(-vit)-vit)*(1.d0-exp(-vp))
         elseif(iqq.eq.2)then !loop sequence with single Pomeron end on one side
          dpx=vi*(min(0.d0,1.d0-exp(-vp)-vp)+vp-vp1)
         elseif(iqq.eq.3)then !loop sequence with single P ends on both sides
          dpx=vi*(vp1-vp0)
         elseif(iqq.eq.4)then !soft loop sequence
          dpx=vis*(min(0.d0,1.d0-exp(-vps)-vps)+vps-vps1)
     *    +min(0.d0,1.d0-exp(-vist)-vist)*(1.d0-exp(-vps))
         elseif(iqq.eq.5)then !soft loop sequence with 1P end on one side
          dpx=vis*(min(0.d0,1.d0-exp(-vps)-vps)+vps-vps1)
         elseif(iqq.eq.6)then !soft loop sequence with 1P ends on both sides
          dpx=vis*(vps1-vps0)
         else 
          stop'qgloop: iqq?!'
         endif
         fann(iqq)=fann(iqq)+a1(ix1)*a1(ix2)*a1(ix3)*dpx/z*rp2
        enddo
        
       else
        vpc0=min(vp0,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,10))  !cut P (two 1P-ends)
        vpc1=min(vpc0,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,9))  !cut P (1P-end)
        vpc=min(vpc1,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,8))   !single cut Pomeron
        vpcn=min(vp,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,14))   !no gap through
        victn=min(vit,qgpini(xpomr*sy,bb1,0.d0,0.d0,14))   !no gap through
        victg=min(victn,qgpini(xpomr*sy,bb1,0.d0,0.d0,15)) !no gap at the end
        vict1=min(victg,qgpini(xpomr*sy,bb1,0.d0,0.d0,17)) !1 cut P at the end
        vict=min(vict1,qgpini(xpomr*sy,bb1,0.d0,0.d0,8))   !single cut Pomeron
        
        vpc0s=min(vpc0,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,13)) !cut soft P (2 ends)
        vpc1s=min(vpc0s,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,12))!cut soft P (1P-end)
        vpcs=min(vpc1s,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,11)) !cut soft Pomeron
        victs=min(vict,qgpini(xpomr*sy,bb1,0.d0,0.d0,11))   !cut soft Pomeron
        do iqq=7,17
         if(iqq.eq.7)then       !single cut Pomeron
          dpx=vi*(vpc*exp(-2.d0*vpcn)-vpc1)
     *    +vict*(exp(-2.d0*victn)-1.d0)*vpc*exp(-2.d0*vpcn)
         elseif(iqq.eq.8)then   !single cut Pomeron with 1P end on one side
          dpx=vi*(vpc*exp(-2.d0*vpcn)-vpc1)
         elseif(iqq.eq.9)then   !single cut Pomeron with 1P end on both sides
          dpx=vi*(vpc1-vpc0)
         elseif(iqq.eq.10)then  !single cut soft Pomeron
          dpx=vis*(vpcs*exp(-2.d0*vpcn)-vpc1s)
     *    +victs*(exp(-2.d0*victn)-1.d0)*vpcs*exp(-2.d0*vpcn)
         elseif(iqq.eq.11)then !single cut soft Pomeron (1P-end on one side)
          dpx=vis*(vpcs*exp(-2.d0*vpcn)-vpc1s)
         elseif(iqq.eq.12)then !single cut soft Pomeron (1P-ends on both sides)
          dpx=vis*(vpc1s-vpc0s)
         elseif(iqq.eq.13)then !any cuts except the complete rap-gap
          dpx=vi*(min(0.d0,1.d0-exp(-vp)-vp)+vp-vp1)
     *    +.5d0*min(0.d0,1.d0-exp(-vit)-vit)*(1.d0-exp(-2.d0*vpcn))
     *    +.5d0*min(0.d0,1.d0-exp(-2.d0*victn)-2.d0*victn)
     *    *max(0.d0,1.d0-exp(-vp)-.5d0*(1.d0-exp(-2.d0*vpcn)))
         elseif(iqq.eq.14)then !no rap-gap at one side
          dpx=vi*(min(0.d0,1.d0-exp(-vp)-vp)+vp-vp1)
     *    +(.5d0*max(0.d0,1.d0-exp(-2.d0*victn)*(1.d0+2.d0*victn))
     *    +victg*(exp(-2.d0*victn)-1.d0))*(1.d0-exp(-vp))
         elseif(iqq.eq.15)then !no rap-gap at one side and 1P-end on the other
          dpx=vi*(vp1-vp0)
     *    +(.5d0*max(0.d0,1.d0-exp(-2.d0*victn)*(1.d0+2.d0*victn))
     *    +victg*(exp(-2.d0*victn)-1.d0))*vp1
         elseif(iqq.eq.16)then !single cut Pomeron at one side
          dpx=vi*(min(0.d0,1.d0-exp(-vp)-vp)+vp-vp1)
     *    +vict1*(exp(-2.d0*victn)-1.d0)*(1.d0-exp(-vp))
         elseif(iqq.eq.17)then !single cut P at one side and 1P-end on the other
          dpx=vi*(vp1-vp0)+vict1*(exp(-2.d0*victn)-1.d0)*vp1
         else 
          stop'qgloop: iqq?!'
         endif
         fann(iqq)=fann(iqq)+a1(ix1)*a1(ix2)*a1(ix3)*dpx/z*rp2
        enddo
       endif
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
1     continue
      dpin=qgpini(sy,bb,0.d0,0.d0,1)
      dpins=qgpini(sy,bb,0.d0,0.d0,0)
      do iqq=1,17
       fann(iqq)=fann(iqq)*log(sy/sgap**2)
     * /8.d0*pi*r3p/.0389d0/g3p**3
       if(iqq.gt.3.and.iqq.lt.7.or.iqq.gt.9.and.iqq.lt.13)then 
        fann(iqq)=fann(iqq)+dpins
       else
        fann(iqq)=fann(iqq)+dpin
       endif
      enddo
      return
      end  

c=============================================================================
      double precision function qgloopr(sy,bb,vvx,iqq,jj)
c-----------------------------------------------------------------------------
c qgloopr - loop sequence coupled to hard evolution (auxilliary procedure)
c sy  - energy scale,
c bb  - impact parameter squared,
c vvx - absorptive factor,
c iqq=1 - g,
c iqq=2 - q,
c jj=1  - uncut,
c jj=2  - cut
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      qgloopr=0.d0
      if(sy.gt.sgap)then
       do ix1=1,7
       do mx1=1,2
        xpomr=(sy/sgap)**(-.5d0-x1(ix1)*(mx1-1.5d0))/sgap
        rp=-alfp*log(xpomr)*2.d0*.0389d0
        rp1=alfp*log(xpomr*sy)*4.d0*.0389d0
        rp2=rp*rp1/(rp+rp1)
        
        do ix2=1,7
        do mx2=1,2
         z=.5d0+x1(ix2)*(mx2-1.5d0)
         bb0=-rp2*log(z)
        do ix3=1,7
        do mx3=1,2
         phi=pi*(.5d0+x1(ix3)*(mx3-1.5d0))
         bb1=(dsqrt(bb)*rp1/(rp+rp1)-dsqrt(bb0)*cos(phi))**2
     *   +bb0*sin(phi)**2
         bb2=(dsqrt(bb)*rp/(rp+rp1)+dsqrt(bb0)*cos(phi))**2
     *   +bb0*sin(phi)**2
       
         vi0=rr*g3p*(xpomr*sy)**dels*qgppdi(1.d0/xpomr/sy,iqq-1)/rp1
     *   *exp(-bb1/rp1)*4.d0*.0389d0
         vis=qgloopri(xpomr*sy,bb1,vvx,iqq,jj)
         if(jj.eq.1)then
          vip1=qgpini(1.d0/xpomr,bb2,0.d0,0.d0,6)
          vip=min(vip1,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,5))
          dpx=vi0*(min(0.d0,1.d0-exp(-vip)-vip)+vip-vip1)
          if(vvx.ne.0.d0)dpx=dpx-vis*(1.d0-exp(-vip))*vvx
         else
          vicn=qgpini(1.d0/xpomr,bb2,0.d0,0.d0,14)
          vipc1=qgpini(1.d0/xpomr,bb2,0.d0,0.d0,12)
          vipc=min(vipc1,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,11))
          dpx=vi0*(vipc*exp(-2.d0*vicn)-vipc1)
          if(vvx.ne.0.d0)dpx=dpx-vis*(1.d0-(1.d0-vvx)**2)
     *    *vipc*exp(-2.d0*vicn)
         endif
         qgloopr=qgloopr+a1(ix1)*a1(ix2)*a1(ix3)*dpx/z*rp2
        enddo
        enddo
        enddo
        enddo
       enddo
       enddo
       qgloopr=qgloopr*dlog(sy/sgap)*pi*r3p/g3p**3/8.d0/.0389d0   !GeV^2
      endif
      rp=alfp*log(sy)*4.d0*.0389d0
      qgloopr=qgloopr/(rr*g3p*sy**dels*qgppdi(1.d0/sy,iqq-1)
     **exp(-bb/rp)/rp*4.d0*.0389d0)+1.d0                            !GeV^2
      return
      end
      
c=============================================================================
      double precision function qgloopri(sy,bb,vvx,iqq,jj)
c-----------------------------------------------------------------------------
c qgloopr - loop sequence coupled to hard evolution (auxilliary procedure)
c sy  - energy scale,
c bb  - impact parameter squared,
c vvx - absorptive factor,
c iqq=1 - g,
c iqq=2 - q, 
c jj=1  - uncut,
c jj=2  - cut
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wz(3),wj(3)
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr20/ spmax
      common /qgarr43/ moniou
      common /qgarr57/ qloopr(31,11,11,2,2)
      common /qgdebug/ debug

      qgloopri=0.d0
      rp=alfp*dlog(sy)*4.d0*.0389d0
      z=exp(-bb/rp)
      if(sy.le.sgap)goto 1
      
      yl=log(sy/sgap)/log(spmax/sgap**2)*30.d0+1.d0
      k=max(1,int(yl-1.d0))
      k=min(k,29)
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)

      z2=.2d0*exp(-12.5d0)
      if(z.gt.z2.or.vvx.ne.0.d0)then
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=(-bb/rp-dlog(0.2d0))/2.5d0+7.d0
       endif
       jz=min(9,int(zz))
       jz=max(2,jz)
       if(jz.eq.6)jz=5
       wz(2)=zz-jz
       wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
       wz(1)=1.d0-wz(2)+wz(3)
       wz(2)=wz(2)-2.d0*wz(3)
      else
       z3=z2*exp(2.5d0)
       jz=1
       wz(1)=(z-z2)*(z-z3)/z2/z3
       wz(2)=z*(z-z3)/z2/(z2-z3)
       wz(3)=z*(z-z2)/z3/(z3-z2)
      endif

      vl=max(1.d0,vvx*10.d0+1.d0)
      if(vvx.eq.0.d0)then
       ivmax=1
       j=1
       wj(1)=1.d0
      else
       j=min(int(vl),9)    
       wj(2)=vl-dble(j)
       wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
       wj(1)=1.d0-wj(2)+wj(3)
       wj(2)=wj(2)-2.d0*wj(3)
       ivmax=3
      endif

      do j1=1,ivmax
       j2=j+j1-1
      do l1=1,3
       l2=jz+l1-1
      do k1=1,3
       k2=k+k1-1
       qgloopri=qgloopri+qloopr(k2,l2,j2,iqq,jj)*wk(k1)*wz(l1)*wj(j1)
      enddo
      enddo
      enddo
     
1     if(z.lt.z2/2.d0.and.vvx.ne.0.d0)qgloopri=min(0.d0,qgloopri)
      qgloopri=exp(qgloopri)*4.d0*rr*sy**dels*qgppdi(1.d0/sy,iqq-1)
     **g3p*z/rp*.0389d0
      return
      end

c------------------------------------------------------------------------
      subroutine qgloos(sy,bb,vvxp,vvxt,fann)
c-----------------------------------------------------------------------
c qgloos - intermediate Pomeron eikonal with screening corrections
c sy   - pomeron mass squared,
c bb   - impact parameter squared,
c vvxp = 1 - exp[-sum_{i} chi_proj(i)] 
c vvxt = 1 - exp[-sum_j chi_targ(j)] 
c iqq=1  - single cut Pomeron
c iqq=2  - diffractive cut, Puu
c iqq=3  - diffractive cut, Puu-Puc
c iqq=4  - diffractive cut, Pcc
c iqq=5  - diffractive cut, Pcc+Pcu
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension fann(17)
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)
      
      do iqq=1,5
       fann(iqq)=0.d0
      enddo
      if(sy.le.sgap**2)goto 1

      do ix1=1,7
      do mx1=1,2
       xpomr=(sy/sgap**2)**(-.5d0-x1(ix1)*(mx1-1.5d0))/sgap
       rp=-alfp*log(xpomr)*4.d0*.0389d0
       rp1=alfp*log(xpomr*sy)*4.d0*.0389d0
       rp2=rp*rp1/(rp+rp1)
      do ix2=1,7
      do mx2=1,2
       z=.5d0+x1(ix2)*(mx2-1.5d0)
       bb0=-rp2*log(z)
      do ix3=1,7
      do mx3=1,2
       phi=pi*(.5d0+x1(ix3)*(mx3-1.5d0))
       bb1=(dsqrt(bb)*rp1/(rp+rp1)-dsqrt(bb0)*cos(phi))**2
     * +bb0*sin(phi)**2
       bb2=(dsqrt(bb)*rp/(rp+rp1)+dsqrt(bb0)*cos(phi))**2
     * +bb0*sin(phi)**2
       
       vit=qgpini(xpomr*sy,bb1,0.d0,0.d0,2)
       vicn=min(vit,qgpini(xpomr*sy,bb1,0.d0,0.d0,14))
       vic1=min(vicn,qgpini(xpomr*sy,bb1,0.d0,0.d0,8))
       
       vu=qgpini(1.d0/xpomr,bb2,0.d0,0.d0,2)
       vcn=min(vu,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,14))
       vc1=min(vcn,qgpini(1.d0/xpomr,bb2,0.d0,0.d0,8))
       
       vpc1=max(0.d0,qgpini(1.d0/xpomr,bb2,vvxp,vvxt,19)
     * +vc1*(exp(-2.d0*vcn)-1.d0))
       vduu=qgpini(1.d0/xpomr,bb2,vvxp,vvxt,20)
       vduc=max(0.d0,vduu-qgpini(1.d0/xpomr,bb2,vvxp,vvxt,21))
       vduu=max(0.d0,vduu+min(0.d0,1.d0-exp(-vu)-vu))
       vdcc=qgpini(1.d0/xpomr,bb2,vvxp,vvxt,22)
       vdcu=max(0.d0,qgpini(1.d0/xpomr,bb2,vvxp,vvxt,23)-vdcc)
       vdcc=max(0.d0,vdcc+.5d0*(1.d0-exp(-vu))**2-vu        !vuncut being added
     * +.5d0*(exp(2.d0*(vu-vcn))-1.d0)*exp(-2.d0*vu))

       do iqq=1,5
        if(iqq.eq.1)then       !single cut Pomeron
         dpx=vpc1*vic1*exp(-2.d0*vicn)
     *   *((1.d0-vvxp)**2*(1.d0-vvxt)**2-1.d0)
        elseif(iqq.eq.2)then   !Puu
         dpx=(1.d0-exp(-vit))
     *   *(vduu*((1.d0-vvxp)*(1.d0-vvxt)*(1.d0-vvxp*vvxt)-1.d0)
     *   -vdcu*(1.d0-vvxp)**2*(1.d0-vvxt)*vvxt)
        elseif(iqq.eq.3)then   !Puu-Puc
         dpx=(1.d0-exp(-vit))
     *   *((vduu-vduc)*((1.d0-vvxp)*(1.d0-vvxt)*(1.d0-vvxp*vvxt)-1.d0)
     *   -(vdcc+vdcu)*(1.d0-vvxp)**2*(1.d0-vvxt)*vvxt)
        elseif(iqq.eq.4)then   !Pcc
         dpx=.5d0*((1.d0-exp(-vit))**2
     *   +(exp(2.d0*(vit-vicn))-1.d0)*exp(-2.d0*vit))
     *   *(vdcc*((1.d0-vvxp)**2*(1.d0-vvxt)**2-1.d0)
     *   -vduc*(1.d0-vvxp)*(1.d0-vvxt)**2*vvxp)
        elseif(iqq.eq.5)then   !Pcc+Pcu
         dpx=.5d0*((1.d0-exp(-vit))**2
     *   +(exp(2.d0*(vit-vicn))-1.d0)*exp(-2.d0*vit))
     *   *((vdcc+vdcu)*((1.d0-vvxp)**2*(1.d0-vvxt)**2-1.d0)
     *   +(vduu-vduc)*(1.d0-vvxp)*(1.d0-vvxt)**2*vvxp)
        endif
        fann(iqq)=fann(iqq)+a1(ix1)*a1(ix2)*a1(ix3)*dpx/z*rp2
       enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
1     vit=qgpini(sy,bb,0.d0,0.d0,2)
      vic1=qgpini(sy,bb,0.d0,0.d0,8)
      do iqq=1,5
       fann(iqq)=fann(iqq)*log(sy/sgap**2)/8.d0*pi*r3p/.0389d0/g3p**3
       if(iqq.eq.1)then
        fann(iqq)=fann(iqq)+vic1
       else
        fann(iqq)=fann(iqq)+vit
       endif
      enddo
      return
      end  

c------------------------------------------------------------------------
      subroutine qglool(sy,bb,icdp,icz,fann)
c-----------------------------------------------------------------------
c qglool - integrated Pomeron leg eikonal with loops
c sy   - pomeron mass squared,
c bb   - impact parameter squared,
c icz  - hadron class
c iqq=1  - all
c iqq=2  - single Pomeron end
c iqq=3  - all soft
c iqq=4  - soft with single Pomeron end
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension fann(17)
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      do iqq=1,4
       fann(iqq)=0.d0
      enddo
      if(sy.le.sgap**2)goto 1
      
      do ix1=1,7
      do mx1=1,2
       xpomr=(sy/sgap**2)**(-.5d0-x1(ix1)*(mx1-1.5d0))/sgap
       rp=(rq(icdp,icz)-alfp*log(xpomr))*4.d0*.0389d0
       rp1=alfp*log(xpomr*sy)*4.d0*.0389d0
       rp2=rp*rp1/(rp+rp1)
      do ix2=1,7
      do mx2=1,2
       z=.5d0+x1(ix2)*(mx2-1.5d0)
       bb0=-rp2*log(z)
      do ix3=1,7
      do mx3=1,2
       phi=pi*(.5d0+x1(ix3)*(mx3-1.5d0))
       bb1=(dsqrt(bb)*rp1/(rp+rp1)-dsqrt(bb0)*cos(phi))**2
     * +bb0*sin(phi)**2
       bb2=(dsqrt(bb)*rp/(rp+rp1)+dsqrt(bb0)*cos(phi))**2
     * +bb0*sin(phi)**2
      
       vpls=qglegi(1.d0/xpomr,bb2,icdp,icz,0)
       vpl=max(vpls,qglegi(1.d0/xpomr,bb2,icdp,icz,1))
       vi0=qgpini(xpomr*sy,bb1,0.d0,0.d0,4)
       vi1=min(vi0,qgpini(xpomr*sy,bb1,0.d0,0.d0,3))
       vi=min(vi1,qgpini(xpomr*sy,bb1,0.d0,0.d0,2))
       vi0s=qgpini(xpomr*sy,bb1,0.d0,0.d0,7)
       vi1s=min(vi0s,qgpini(xpomr*sy,bb1,0.d0,0.d0,6))
       vis=min(vi1s,qgpini(xpomr*sy,bb1,0.d0,0.d0,5))
       
       do iqq=1,4
        if(iqq.eq.1)then
         dpx=vpl*(min(0.d0,1.d0-exp(-vi)-vi)+vi-vi1)
        elseif(iqq.eq.2)then
         dpx=vpl*(vi1-vi0)
        elseif(iqq.eq.3)then
         dpx=vpls*(min(0.d0,1.d0-exp(-vis)-vis)+vis-vi1s)
        elseif(iqq.eq.4)then
         dpx=vpls*(vi1s-vi0s)
        else
         stop'qglool: iqq?!'
        endif
        fann(iqq)=fann(iqq)+a1(ix1)*a1(ix2)*a1(ix3)*dpx/z*rp2
       enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
1     dlool=qglegi(sy,bb,icdp,icz,1)
      dlools=qglegi(sy,bb,icdp,icz,0)
      do iqq=1,4
       fann(iqq)=fann(iqq)*log(sy/sgap**2)/8.d0*pi*r3p/.0389d0/g3p**3
       if(iqq.gt.2)then
        fann(iqq)=fann(iqq)+dlools
       else
        fann(iqq)=fann(iqq)+dlool
       endif
      enddo
      return
      end 

c------------------------------------------------------------------------
      subroutine qgfan(sy,bb,vvx,icdp,icz,fann)
c-----------------------------------------------------------------------
c qgfan - integrated fan-contributions
c sy    - c.m. energy squared,
c bb    - impact parameter squared,
c icdp - diffractive state for the projectile,
c icz  - hadron class
c vvx  = 1 - exp[-sum_j chi_targ(j) - sum_{i<I} chi_proj(i)]
c iqq=1 - general fan with loops
c iqq=2 - general fan with single Pomeron end
c iqq=3 - soft fan with loops
c iqq=4 - soft fan with single Pomeron end
c iqq=5 - q-fan with loops
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension fann(17)
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)
      common /arr6/    x3(8),a3(8)
      
      do iqq=1,5
       fann(iqq)=0.d0
      enddo
      if(sy.le.sgap**2)goto 1
      
      do ix1=1,7
      do mx1=1,2
       xpomr1=(sy/sgap**2)**(-.5d0-x1(ix1)*(mx1-1.5d0))/sgap
       rp=(rq(icdp,icz)-alfp*log(xpomr1))*4.d0*.0389d0
       rp1=alfp*log(xpomr1*sy)*4.d0*.0389d0
       rp2=rp*rp1/(rp+rp1)
       do ix2=1,7
       do mx2=1,2
        z=.5d0+x1(ix2)*(mx2-1.5d0)
        bb0=-rp2*log(z)
       do ix3=1,7
       do mx3=1,2
        phi=pi*(.5d0+x1(ix3)*(mx3-1.5d0))
        bb1=(dsqrt(bb)*rp1/(rp+rp1)-dsqrt(bb0)*cos(phi))**2
     *  +bb0*sin(phi)**2
        bb2=(dsqrt(bb)*rp/(rp+rp1)+dsqrt(bb0)*cos(phi))**2
     *  +bb0*sin(phi)**2
       
        vpls=qglegi(1.d0/xpomr1,bb2,icdp,icz,0)          !single soft Pomeron
        vpl=qglegi(1.d0/xpomr1,bb2,icdp,icz,1)           !single Pomeron (s+g)
        vplq=qglegi(1.d0/xpomr1,bb2,icdp,icz,2)          !single q-Pomeron 

        vpfq=min(vplq,qgfani(1.d0/xpomr1,bb2,vvx,0.d0,0.d0,icdp,icz,5))!q-fan
        vpfs=min(vpls,qgfani(1.d0/xpomr1,bb2,vvx,0.d0,0.d0,icdp,icz,3))!soft fan
        vpf=min(vpl,qgfani(1.d0/xpomr1,bb2,vvx,0.d0,0.d0,icdp,icz,1))  !fan
      
        vi0=qgpini(xpomr1*sy,bb1,0.d0,0.d0,4)          !Pomeron loop (2 1P-ends)
        vi1=min(vi0,qgpini(xpomr1*sy,bb1,0.d0,0.d0,3)) !Pomeron loop (1P-end)
        vi=min(vi1,qgpini(xpomr1*sy,bb1,0.d0,0.d0,2))  !general Pomeron loop
        vi0s=qgpini(xpomr1*sy,bb1,0.d0,0.d0,7)         !soft P loop (2 1P-ends)
        vi1s=min(vi0s,qgpini(xpomr1*sy,bb1,0.d0,0.d0,6)) !soft P loop (1P-end)
        vis=min(vi1s,qgpini(xpomr1*sy,bb1,0.d0,0.d0,5))  !soft Pomeron loop
     
        do iqq=1,5
         if(iqq.eq.1)then         !general fan with loops
          dpx=(1.d0-exp(-vi))*(min(0.d0,1.d0-exp(-vpf)-vpf)
     *    *(1.d0-vvx)-vpf*vvx+vpl)-vpl*vi1
         elseif(iqq.eq.2)then     !general fan with single Pomeron end
          dpx=vi1*(min(0.d0,1.d0-exp(-vpf)-vpf)*(1.d0-vvx)-vpf*vvx)
     *    +vpl*(vi1-vi0)
         elseif(iqq.eq.3)then     !soft fan with loops
          dpx=(1.d0-exp(-vis))*(min(0.d0,1.d0-exp(-vpfs)-vpfs)
     *    *(1.d0-vvx)-vpfs*vvx+vpls)-vpls*vi1s
         elseif(iqq.eq.4)then     !soft fan with single Pomeron end
          dpx=vi1s*(min(0.d0,1.d0-exp(-vpfs)-vpfs)*(1.d0-vvx)-vpfs*vvx)
     *    +vpls*(vi1s-vi0s)
         elseif(iqq.eq.5)then     !q-fan with loops
          dpx=(1.d0-exp(-vis))*(vpfq*(exp(-vpf)*(1.d0-vvx)-1.d0)+vplq)
     *    -vplq*vi1s
         endif
         fann(iqq)=fann(iqq)+a1(ix1)*a1(ix2)*a1(ix3)*dpx/z*rp2
        enddo
       enddo
       enddo
       enddo
       enddo
      enddo
      enddo

1     flegs=qglegi(sy,bb,icdp,icz,0)
      fleg=qglegi(sy,bb,icdp,icz,1)
      flegq=qglegi(sy,bb,icdp,icz,2)+flegs
      do iqq=1,5
       fann(iqq)=fann(iqq)*dlog(sy/sgap**2)/8.d0*pi*r3p/.0389d0/g3p**3
       if(iqq.le.2)then
        fann(iqq)=fann(iqq)+fleg
       elseif(iqq.le.4)then
        fann(iqq)=fann(iqq)+flegs
       else
        fann(iqq)=fann(iqq)+flegq
       endif
      enddo
      return
      end  
    
c------------------------------------------------------------------------
      subroutine qgfanc(sy,bb,vvx,vvxp,vvxpl,icdp,icz,fann)
c-----------------------------------------------------------------------
c qgfan - cut fan-contributions
c sy    - c.m. energy squared,
c bb    - impact parameter squared,
c icdp  - diffractive state for the projectile,
c icz   - hadron class,
c vvx  = 1 - exp[-sum_j chi_targ(j) - sum_{i<I} chi_proj(i)]
c vvxpl= 1 - exp[-sum_{i<I} chi_proj(i)] 
c vvxp = 1 - exp[-sum_{i>I} chi^(6)_proj(i)] (iqq=1,2,3)
c vvxp = 1 - exp[-sum_{i>I} chi^(9)_proj(i)] (iqq=4)
c vvxp = 1 - exp[-sum_{i>I} chi_proj(i)]  (iqq=5-13)
c iqq=1  - cut handle fan
c iqq=2  - no rap-gap at the end
c iqq=3  - single cut Pomeron end
c iqq=4  - total fan-like contribution
c iqq=5  - leg-like cut (ND)
c iqq=6  - leg-like cut with cut handle
c iqq=7  - leg-like cut without a rap-gap at the end
c iqq=8  - leg-like cut with single cut Pomeron end
c iqq=9  - single cut Pomeron
c iqq=10 - single cut Pomeron with single Pomeron end
c iqq=11 - single cut soft Pomeron
c iqq=12 - single cut soft Pomeron with single Pomeron end
c iqq=13 - single cut q-Pomeron
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension fann(17)
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)
      
      do iqq=1,13
       fann(iqq)=0.d0
      enddo
      if(sy.le.sgap**2)goto 1
      
      if(vvx.gt..999d0)then
       vvxs=0.d0
      else
       vvxs=(1.d0-vvx)**2/(1.d0-vvxpl)
      endif
      
      do ix1=1,7
      do mx1=1,2
       xpomr1=(sy/sgap**2)**(-.5d0-x1(ix1)*(mx1-1.5d0))/sgap
       rp=(rq(icdp,icz)-alfp*log(xpomr1))*4.d0*.0389d0
       rp1=alfp*log(xpomr1*sy)*4.d0*.0389d0
       rp2=rp*rp1/(rp+rp1)
       do ix2=1,7
       do mx2=1,2
        z=.5d0+x1(ix2)*(mx2-1.5d0)
        bb0=-rp2*log(z)
       do ix3=1,7
       do mx3=1,2
        phi=pi*(.5d0+x1(ix3)*(mx3-1.5d0))
        bb1=(dsqrt(bb)*rp1/(rp+rp1)-dsqrt(bb0)*cos(phi))**2
     *  +bb0*sin(phi)**2
        bb2=(dsqrt(bb)*rp/(rp+rp1)+dsqrt(bb0)*cos(phi))**2
     *  +bb0*sin(phi)**2
       
        vi1=qgpini(xpomr1*sy,bb1,0.d0,0.d0,3)
        vi=min(vi1,qgpini(xpomr1*sy,bb1,0.d0,0.d0,2))
        vicn=min(vi,qgpini(xpomr1*sy,bb1,0.d0,0.d0,14))
        vicgap1=qgpini(xpomr1*sy,bb1,0.d0,0.d0,16)
        vicgap=min(vicgap1,qgpini(xpomr1*sy,bb1,0.d0,0.d0,15))
        vic1p0=min(vicgap1,qgpini(xpomr1*sy,bb1,0.d0,0.d0,18))
        vic1p=min(vic1p0,qgpini(xpomr1*sy,bb1,0.d0,0.d0,17))
        vic0=qgpini(xpomr1*sy,bb1,0.d0,0.d0,10)
        vic1=min(vic0,qgpini(xpomr1*sy,bb1,0.d0,0.d0,9))
        vic=min(vic1,qgpini(xpomr1*sy,bb1,0.d0,0.d0,8))

        vic0s=qgpini(xpomr1*sy,bb1,0.d0,0.d0,13)
        vic1s=min(vic0s,qgpini(xpomr1*sy,bb1,0.d0,0.d0,12))
        vics=min(vic1s,qgpini(xpomr1*sy,bb1,0.d0,0.d0,11))

        vpls=qglegi(1.d0/xpomr1,bb2,icdp,icz,0)
        vpl=qglegi(1.d0/xpomr1,bb2,icdp,icz,1)
        vplq=qglegi(1.d0/xpomr1,bb2,icdp,icz,2)
       
        vpf=min(vpl,qgfani(1.d0/xpomr1,bb2,vvx,0.d0,0.d0,icdp,icz,1))
        vpfq=qgfani(1.d0/xpomr1,bb2,vvx,0.d0,0.d0,icdp,icz,5)
        vpfc0=min(vpf,qgfani(1.d0/xpomr1,bb2,vvx,vvxp,vvxpl 
     *  ,icdp,icz,6))
        vpfct=qgfani(1.d0/xpomr1,bb2,vvx,vvxp,vvxpl,icdp,icz,9)
        vpfl=min(vpf,qgfani(1.d0/xpomr1,bb2,vvx,vvxp,vvxpl
     *  ,icdp,icz,10))
        vpfl0=min(vpfl,qgfani(1.d0/xpomr1,bb2,vvx,vvxp,vvxpl
     *  ,icdp,icz,11))
        vpfcq=qgfani(1.d0/xpomr1,bb2,vvx,vvxp,vvxpl,icdp,icz,18)
        vpfc=qgfani(1.d0/xpomr1,bb2,vvx,vvxp,vvxpl,icdp,icz,14)
        vpfcs=qgfani(1.d0/xpomr1,bb2,vvx,vvxp,vvxpl,icdp,icz,16)
     
        do iqq=1,13
         if(iqq.eq.1)then          !cut handle
          dpx=(1.d0-exp(-vi))*(vvxs*(min(0.d0,1.d0-exp(-vpfc0)-vpfc0)
     *    *(1.d0-vvxp)+vvxp*(1.d0-exp(-vpf)))
     *    +vpfc0*(vvxs*(1.d0-vvxp)-1.d0)+vpl)-vpl*vi1
         elseif(iqq.eq.2)then      !no rap-gap at the end
          dpx=(.5d0*max(0.d0,1.d0-exp(-2.d0*vicn)*(1.d0+2.d0*vicn))
     *    +vicgap*exp(-2.d0*vicn))*(vvxs*(min(0.d0
     *    ,1.d0-exp(-vpfc0)-vpfc0)*(1.d0-vvxp)+vvxp*(1.d0-exp(-vpf)))
     *    +vpfc0*(vvxs*(1.d0-vvxp)-1.d0)+vpl)-vpl*vicgap1
         elseif(iqq.eq.3)then      !single cut Pomeron end
          dpx=vic1p*exp(-2.d0*vicn)*(vvxs*(min(0.d0
     *    ,1.d0-exp(-vpfc0)-vpfc0)*(1.d0-vvxp)+vvxp*(1.d0-exp(-vpf)))
     *    +vpfc0*(vvxs*(1.d0-vvxp)-1.d0)+vpl)-vpl*vic1p0
         elseif(iqq.eq.4)then      !total fan-like contribution
          dpx=(1.d0-exp(-vi))*((1.d0-vvxpl)*(min(0.d0
     *    ,1.d0-exp(-vpfct)-vpfct)*(1.d0-vvxp)+vvxp*(1.d0-exp(-vpf)))
     *    -vpfct*(1.d0-(1.d0-vvxp)*(1.d0-vvxpl))+vpl)-vpl*vi1
         elseif(iqq.eq.5)then      !leg-like cut
          dpx=(1.d0-exp(-vi))*(vpfl*((1.d0-vvx)*(1.d0-vvxpl)
     *    *(1.d0-vvxp)**2*exp(-2.d0*vpf)-1.d0)+vpl)-vpl*vi1
         elseif(iqq.eq.6)then      !leg-like cut with cut handle
          dpx=(1.d0-exp(-vi))
     *    *(vpfl0*((1.d0-vvx)**2*(1.d0-vvxp)**2*exp(-2.d0*vpf)-1.d0)
     *    -(vpfl-vpfl0)*vvxs*(1.d0-vvxp)*exp(-vpf)
     *    *(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpf))+vpl)-vpl*vi1
         elseif(iqq.eq.7)then      !leg-like cut without a rap-gap at the end
          dpx=(.5d0*max(0.d0,1.d0-exp(-2.d0*vicn)*(1.d0+2.d0*vicn))
     *    +vicgap*exp(-2.d0*vicn))
     *    *(vpfl0*((1.d0-vvx)**2*(1.d0-vvxp)**2*exp(-2.d0*vpf)-1.d0)
     *    -(vpfl-vpfl0)*vvxs*(1.d0-vvxp)*exp(-vpf)
     *    *(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpf))+vpl)-vpl*vicgap1
         elseif(iqq.eq.8)then      !leg-like cut with single cut Pomeron end
          dpx=vic1p*exp(-2.d0*vicn)
     *    *(vpfl0*((1.d0-vvx)**2*(1.d0-vvxp)**2*exp(-2.d0*vpf)-1.d0)
     *    -(vpfl-vpfl0)*vvxs*(1.d0-vvxp)*exp(-vpf)
     *    *(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpf))+vpl)-vpl*vic1p0
         elseif(iqq.eq.9)then      !single cut Pomeron
          dpx=vic*exp(-2.d0*vicn)*(vpfc*((1.d0-vvx)**2
     *    *(1.d0-vvxp)**2*exp(-2.d0*vpf-2.d0*vpfq)-1.d0)+vpl)-vpl*vic1
         elseif(iqq.eq.10)then     !single cut Pomeron with single Pomeron end
          dpx=vic1*vpfc*((1.d0-vvx)**2*(1.d0-vvxp)**2
     *    *exp(-2.d0*vpf-2.d0*vpfq)-1.d0)+vpl*(vic1-vic0)
         elseif(iqq.eq.11)then     !single cut soft Pomeron
          dpx=vics*exp(-2.d0*vicn)*(vpfcs*((1.d0-vvx)**2*(1.d0-vvxp)**2
     *    *exp(-2.d0*vpf-2.d0*vpfq)-1.d0)+vpls)-vpls*vic1s
         elseif(iqq.eq.12)then     !single cut soft Pomeron with single 1P-end
          dpx=vic1s*vpfcs*((1.d0-vvx)**2*(1.d0-vvxp)**2
     *    *exp(-2.d0*vpf-2.d0*vpfq)-1.d0)+vpls*(vic1s-vic0s)
         elseif(iqq.eq.13)then     !single cut q-Pomeron
          dpx=vics*exp(-2.d0*vicn)*(vpfcq*((1.d0-vvx)**2
     *    *(1.d0-vvxp)**2*exp(-2.d0*vpf)-1.d0)+vplq)-vplq*vic1s
         endif
         fann(iqq)=fann(iqq)+a1(ix1)*a1(ix2)*a1(ix3)*dpx/z*rp2
        enddo
       enddo
       enddo
       enddo
       enddo
      enddo
      enddo

1     dfan=qglegi(sy,bb,icdp,icz,1)
      dfans=qglegi(sy,bb,icdp,icz,0)
      dfanq=qglegi(sy,bb,icdp,icz,2)+dfans
      do iqq=1,13
       fann(iqq)=fann(iqq)*dlog(sy/sgap**2)/8.d0*pi*r3p/.0389d0/g3p**3
       if(iqq.le.10)then
        fann(iqq)=fann(iqq)+dfan
       elseif(iqq.le.12)then
        fann(iqq)=fann(iqq)+dfans
       elseif(iqq.eq.13)then
        fann(iqq)=fann(iqq)+dfanq
       else
        stop'qgfanc: iqq?!'
       endif
      enddo
      return
      end 
       
c------------------------------------------------------------------------
      double precision function qgfani(sy,bb,vvx,vvxp,vvxpl,icdp,icz
     *,iqq)
c-----------------------------------------------------------------------
c qgfani - integrated fan-contributions
c sy   - c.m. energy squared,
c bb   - impact parameter squared,
c icdp - diffractive state for the projectile,
c icz  - hadron class,
c vvx  = 1 - exp[-sum_j chi_targ(j) - sum_{i<I} chi_proj(i)]
c vvxp = vvxpl = 0                           (iqq=1-5)
c vvxpl= 1 - exp[-sum_{i<I} chi_proj(i)] 
c vvxp = 1 - exp[-sum_{i>I} chi^(6)_proj(i)] (iqq=6,7,8)
c vvxp = 1 - exp[-sum_{i>I} chi^(9)_proj(i)] (iqq=9)
c vvxp = 1 - exp[-sum_{i>I} chi_proj(i)]     (iqq=10-18)
c iqq=1  - general fan with loops
c iqq=2  - general fan with single Pomeron end
c iqq=3  - soft fan with loops
c iqq=4  - soft fan with single Pomeron end
c iqq=5  - q-fan with loops
c iqq=6  - cut handle fan
c iqq=7  - cut handle fan (no rap-gap at the end)
c iqq=8  - cut handle fan (single cut Pomeron end)
c iqq=9  - total fan-like contribution
c iqq=10 - leg-like cut
c iqq=11 - leg-like cut with cut handle
c iqq=12 - leg-like cut without a rap-gap at the end
c iqq=13 - leg-like cut with single Pomeron end
c iqq=14 - single cut Pomeron (with screening corrections)
c iqq=15 - single cut Pomeron with single Pomeron end
c iqq=16 - single cut soft Pomeron (with screening corrections)
c iqq=17 - single cut soft Pomeron with single Pomeron end
c iqq=18 - single cut q-Pomeron (with screening corrections)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wz(3),wj(3),wi(3),wn(3)
      parameter(nfock=3)
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr20/ spmax
      common /qgarr27/ qlegi(31,11,nfock,3,6),qfanu(31,11,11,3*nfock,5)
     *,qfanc(31,11,11,39,13*nfock)
      common /qgarr43/ moniou
      common /qgdebug/ debug
            
      qgfani=0.d0
      if(sy.lt..999d0*sgap)stop'qgfani:sy!!!'
      if(sy.le.sgap**2)goto 1

      yl=dlog(sy/sgap**2)/dlog(spmax/sgap**3)*30.d0+1.d0
      k=max(1,int(yl-1.d0))
      k=min(k,29)
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)
      
      rp=(rq(icdp,icz)+alfp*dlog(sy))*4.d0*.0389d0
      z=dexp(-bb/rp)
      z2=.2d0*exp(-12.5d0)
      z3=z2*exp(2.5d0)
      if(z.gt.z2.or.vvx+vvxp+vvxpl.ne.0.d0)then
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=(-bb/rp-dlog(0.2d0))/2.5d0+7.d0
       endif
       jz=min(9,int(zz))
       jz=max(2,jz)
       if(jz.eq.6)jz=5
       wz(2)=zz-jz
       wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
       wz(1)=1.d0-wz(2)+wz(3)
       wz(2)=wz(2)-2.d0*wz(3)
      else
       jz=1
       wz(1)=(z-z2)*(z-z3)/z2/z3
       wz(2)=z*(z-z3)/z2/(z2-z3)
       wz(3)=z*(z-z2)/z3/(z3-z2)
      endif

      vl=max(1.d0,vvx*10.d0+1.d0)
      if(vvx.eq.0.d0)then
       ivmax=1
       j=1
       wj(1)=1.d0
      else
       j=min(int(vl),9)    
       wj(2)=vl-dble(j)
       wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
       wj(1)=1.d0-wj(2)+wj(3)
       wj(2)=wj(2)-2.d0*wj(3)
       ivmax=3
      endif

      if(iqq.le.5)then
       ii=icdp+nfock*(icz-1)
       do j1=1,ivmax
        j2=j+j1-1
       do l1=1,3
        l2=jz+l1-1
       do k1=1,3
        k2=k+k1-1
        qgfani=qgfani+qfanu(k2,l2,j2,ii,iqq)*wk(k1)*wz(l1)*wj(j1)
       enddo
       enddo
       enddo
       
      elseif(icz.ne.2.or.vvxp+vvxpl.eq.0.d0)then !hadron (no proj. nucl. corr.)
       ll=icz+(icz-1)*(3-icz)*2
       ii=icdp+nfock*(iqq-6)
       do j1=1,ivmax
        j2=j+j1-1
       do l1=1,3
        l2=jz+l1-1
       do k1=1,3
        k2=k+k1-1
        qgfani=qgfani+qfanc(k2,l2,j2,ll,ii)*wk(k1)*wz(l1)*wj(j1)
       enddo
       enddo
       enddo

      else
       if(vvxp.eq.0.d0)then         !hadron (no nuclear graphs)
        iv1max=1
        i=1
        wi(1)=1.d0
       else                         !nuclear effects 
        iv1max=2
        vl1=max(1.d0,vvxp*5.d0+1.d0)
        i=min(int(vl1),5)    
        wi(2)=vl1-i
        wi(1)=1.d0-wi(2)
       endif
        
       if(vvx.lt..01d0.or.vvxpl.eq.0.d0)then    !weak (no) screening
        iv2max=1
         n=1
        wn(1)=1.d0
       else                                    !nuclear effects 
        iv2max=2
        vl2=max(1.d0,vvxpl/vvx*5.d0+1.d0)
        n=min(int(vl2),5)    
        wn(2)=vl2-n
        wn(1)=1.d0-wn(2)
       endif

       ii=icdp+nfock*(iqq-6)
       do n1=1,iv2max
        n2=n+n1-2
       do i1=1,iv1max
        i2=i+i1+2
       do j1=1,ivmax
        j2=j+j1-1
       do l1=1,3
        l2=jz+l1-1
       do k1=1,3
        k2=k+k1-1
        qgfani=qgfani+qfanc(k2,l2,j2,i2+6*n2,ii)
     *  *wk(k1)*wz(l1)*wj(j1)*wi(i1)*wn(n1)
       enddo
       enddo
       enddo
       enddo
       enddo
      endif
      if(z.lt.z2/2.d0.and.vvx+vvxp+vvxpl.ne.0.d0)
     *qgfani=min(0.d0,qgfani)
      
1     qgfani=dexp(qgfani)
      if(iqq.eq.5.or.iqq.eq.18)qgfani=max(0.d0,qgfani-1.d0)
      qgfani=qgfani*qglegi(sy,bb,icdp,icz,0)
      return 
      end

c=============================================================================
      double precision function qg3pom(sy,b,vvx,vvxp,vvxt
     *,icdp,icdt,icz,iqq)
c-----------------------------------------------------------------------
c qg3pom - integrated 3p-contributions to the interaction eikonal
c sy   - pomeron mass squared,
c b    - impact parameter,
c vvx  = 1 - exp[-sum_{j<J} chi_targ(j) - sum_{i<I} chi_proj(i)]
c vvxp = 1 - exp[-sum_{i>I} chi_proj(i)] 
c vvxt = 1 - exp[-sum_{j>J} chi_targ(j)]
c icdp - diffractive state for the projectile,
c icdt - diffractive state for the target,
c icz  - projectile class,
c iqq=1 - q_v not involved,
c iqq=2 - q_v involved
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)
      
      qg3pom=0.d0      
      if(sy.le.sgap**2)return

      do ix1=1,7
      do mx1=1,2
       xpomr1=(sy/sgap**2)**(-(.5+x1(ix1)*(mx1-1.5)))/sgap
       rp1=(rq(icdp,icz)-alfp*log(xpomr1))*4.d0*.0389d0
       rp2=(rq(icdt,2)+alfp*log(xpomr1*sy))*4.d0*.0389d0
       rp=rp1*rp2/(rp1+rp2)
      do ib1=1,7
      do mb1=1,2
       z=.5d0+x1(ib1)*(mb1-1.5d0)
       bb0=-rp*dlog(z)
      do ib2=1,7
      do mb2=1,2
       phi=pi*(.5d0+x1(ib2)*(mb2-1.5d0))
       bb1=(b*rp1/(rp1+rp2)+dsqrt(bb0)*cos(phi))**2+bb0*sin(phi)**2
       bb2=(b*rp2/(rp1+rp2)-dsqrt(bb0)*cos(phi))**2+bb0*sin(phi)**2
               
       v1p0=qglegi(1.d0/xpomr1,bb1,icdp,icz,1)           !single Pomeron (s+g)
       v1t0=qglegi(xpomr1*sy,bb2,icdt,2,1)               !single Pomeron (s+g)
       v1p1=min(v1p0,qglegi(1.d0/xpomr1,bb1,icdp,icz,4)) !single leg with 1P-end
       v1t1=min(v1t0,qglegi(xpomr1*sy,bb2,icdt,2,4))     !single leg with 1P-end
       v1p=min(v1p1,qglegi(1.d0/xpomr1,bb1,icdp,icz,3))  !single leg
       v1t=min(v1t1,qglegi(xpomr1*sy,bb2,icdt,2,3))      !single leg

       vpf0=min(v1p,qgfani(1.d0/xpomr1,bb1
     * ,1.d0-(1.d0-vvx)*(1.d0-vvxt),0.d0,0.d0,icdp,icz,1))
       vtf0=min(v1t,qgfani(xpomr1*sy,bb2
     * ,1.d0-(1.d0-vvx)*(1.d0-vvxp),0.d0,0.d0,icdt,2,1))
       
       n=1
1      n=n+1
       vpf=min(v1p,qgfani(1.d0/xpomr1,bb1
     * ,1.d0-(1.d0-vvx)*(1.d0-vvxt)*exp(-vtf0),0.d0,0.d0,icdp,icz,1)) !fan
       vtf=min(v1t,qgfani(xpomr1*sy,bb2
     * ,1.d0-(1.d0-vvx)*(1.d0-vvxp)*exp(-vpf0),0.d0,0.d0,icdt,2,1))   !fan
       if(abs(1.d0-vpf/vpf0)+abs(1.d0-vtf/vtf0).gt.1.d-2.and.n.lt.100)
     * then
        vpf0=vpf
        vtf0=vtf
        goto 1
       endif

       if(iqq.eq.1)then
        dpx=(1.d0-vvx)*(min(0.d0,1.d0-exp(-vpf)-vpf)
     *  *min(0.d0,1.d0-exp(-vtf)-vtf)
     *  +vpf*min(0.d0,1.d0-exp(-vtf)-vtf)
     *  +vtf*min(0.d0,1.d0-exp(-vpf)-vpf))-vvx*vpf*vtf
     *  -.5d0*(vtf-v1t)*(min(0.d0,1.d0-exp(-vpf)-vpf)
     *  *(1.d0-vvx)*(1.d0-vvxt)*exp(-vtf)
     *  -vpf*(1.d0-(1.d0-vvx)*(1.d0-vvxt)*exp(-vtf)))
     *  -.5d0*(vpf-v1p)*(min(0.d0,1.d0-exp(-vtf)-vtf)
     *  *(1.d0-vvx)*(1.d0-vvxp)*exp(-vpf)
     *  -vtf*(1.d0-(1.d0-vvx)*(1.d0-vvxp)*exp(-vpf)))
     *  +.5d0*(v1t-v1t1)*v1p0+.5d0*(v1p-v1p1)*v1t0
       else
        vpq=qglegi(1.d0/xpomr1,bb1,icdp,icz,2)  !single q-Pomeron
        vtq=qglegi(xpomr1*sy,bb2,icdt,2,2)      !single q-Pomeron
        vpf1=max(vpf,qgfani(1.d0/xpomr1,bb1
     *  ,1.d0-(1.d0-vvx)*(1.d0-vvxt)*exp(-vtf0),0.d0,0.d0,icdp,icz,2))
        vpf1=min(v1p,vpf1)
        vtf1=max(vtf,qgfani(xpomr1*sy,bb2
     *  ,1.d0-(1.d0-vvx)*(1.d0-vvxp)*exp(-vpf0),0.d0,0.d0,icdt,2,2))
        vtf1=min(v1t,vtf1)

        dpx=vpq*((1.d0-exp(-vtf))*exp(-vpf)*(1.d0-vvx)*(1.d0-vvxp)
     *  -vtf1)
     *  +vtq*((1.d0-exp(-vpf))*exp(-vtf)*(1.d0-vvx)*(1.d0-vvxt)-vpf1)
       endif
       dpx=min(1.d0,dpx)
       qg3pom=qg3pom+a1(ib1)*a1(ib2)*a1(ix1)/z*rp*dpx
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      qg3pom=qg3pom/8.d0*log(sy/sgap**2)*(r3p*pi/.0389d0)/g3p**3
      return 
      end

c------------------------------------------------------------------------
      subroutine qgpcut(sy,b,vvx,vvxp,vvxt,icdp,icdt,icz,v1p)
c-----------------------------------------------------------------------
c qgpcut - all screening corrections for single cut Pomeron contribution
c sy   - pomeron mass squared,
c bb   - impact parameter squared,
c vvx  = 1 - exp[-sum_{j<J} chi_targ(j) - sum_{i<I} chi_proj(i)]
c vvxp = 1 - exp[-sum_{i>I} chi_proj(i)] 
c vvxt = 1 - exp[-sum_{j>J} chi_targ(j)] 
c icdp - diffractive state for the projectile,
c icdt - diffractive state for the target,
c icz  - projectile class,
c iqq=1 - general Pomeron (s+g),
c iqq=2 - soft Pomeron
c iqq=3 - q_v-Pomeron (qg+gq)
c iqq=4 - qg-Pomeron
c iqq=5 - gq-Pomeron
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      dimension v1p(5)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      do i=1,5
       v1p(i)=0.d0
      enddo
      if(sy.le.sgap**2)return

      do ix1=1,7
      do mx1=1,2
       xpomr1=(sy/sgap**2)**(-(.5+x1(ix1)*(mx1-1.5)))/sgap
       rp1=(rq(icdp,icz)-alfp*log(xpomr1))*4.d0*.0389d0
       rp2=(rq(icdt,2)+alfp*log(xpomr1*sy))*4.d0*.0389d0
       rp=rp1*rp2/(rp1+rp2)
      do ib1=1,7
      do mb1=1,2
       z=.5d0+x1(ib1)*(mb1-1.5d0)
       bb0=-rp*dlog(z)
      do ib2=1,7
      do mb2=1,2
       phi=pi*(.5d0+x1(ib2)*(mb2-1.5d0))
       bb1=(b*rp1/(rp1+rp2)+dsqrt(bb0)*cos(phi))**2+bb0*sin(phi)**2
       bb2=(b*rp2/(rp1+rp2)-dsqrt(bb0)*cos(phi))**2+bb0*sin(phi)**2
               
       vpls=qglegi(1.d0/xpomr1,bb1,icdp,icz,0)             !single soft Pomeron
       vtls=qglegi(xpomr1*sy,bb2,icdt,2,0)                 !single soft Pomeron
       vpl=max(vpls,qglegi(1.d0/xpomr1,bb1,icdp,icz,1))    !single Pomeron (s+g)
       vtl=max(vtls,qglegi(xpomr1*sy,bb2,icdt,2,1))        !single Pomeron (s+g)
       vplq=qglegi(1.d0/xpomr1,bb1,icdp,icz,2)             !single q-Pomeron
       vtlq=qglegi(xpomr1*sy,bb2,icdt,2,2)                 !single q-Pomeron

       vpf0=qgfani(1.d0/xpomr1,bb1,1.d0-(1.d0-vvx)*(1.d0-vvxt)
     * ,0.d0,0.d0,icdp,icz,1)
       vtf0=qgfani(xpomr1*sy,bb2,1.d0-(1.d0-vvx)*(1.d0-vvxp)
     * ,0.d0,0.d0,icdt,2,1)
       n=1
1      n=n+1
       vpfq=qgfani(1.d0/xpomr1,bb1,1.d0-(1.d0-vvx)*(1.d0-vvxt)
     * *exp(-vtf0),0.d0,0.d0,icdp,icz,5)               !q-fan contribution
       vtfq=qgfani(xpomr1*sy,bb2,1.d0-(1.d0-vvx)*(1.d0-vvxp)
     * *exp(-vpf0),0.d0,0.d0,icdt,2,5)                 !q-fan contribution
       vpf=qgfani(1.d0/xpomr1,bb1,1.d0-(1.d0-vvx)*(1.d0-vvxt)
     * *exp(-vtf0-vtfq),0.d0,0.d0,icdp,icz,1)          !general fan contribution
       vtf=qgfani(xpomr1*sy,bb2,1.d0-(1.d0-vvx)*(1.d0-vvxp)
     * *exp(-vpf0-vpfq),0.d0,0.d0,icdt,2,1)            !general fan contribution
       if(abs(1.d0-vpf/vpf0)+abs(1.d0-vtf/vtf0).gt.1.d-2.and.n.lt.100)
     * then
        vpf0=vpf
        vtf0=vtf
        goto 1
       endif

       vplc1=qgfani(1.d0/xpomr1,bb1,1.d0-(1.d0-vvx)*(1.d0-vvxt)
     * *exp(-vtf-vtfq),vvxp,0.d0,icdp,icz,15)      !single cut Pomeron (1P-end)
       vplc=min(vplc1,qgfani(1.d0/xpomr1,bb1,1.d0-(1.d0-vvx)
     * *(1.d0-vvxt)*exp(-vtf-vtfq),vvxp,0.d0,icdp,icz,14)) !single cut P
       vtlc1=qgfani(xpomr1*sy,bb2,1.d0-(1.d0-vvx)*(1.d0-vvxp)
     * *exp(-vpf-vpfq),vvxt,0.d0,icdt,2,15)        !single cut Pomeron (1P-end)
       vtlc=min(vtlc1,qgfani(xpomr1*sy,bb2,1.d0-(1.d0-vvx)*(1.d0-vvxp)
     * *exp(-vpf-vpfq),vvxt,0.d0,icdt,2,14))       !single cut Pomeron in a fan
       vplc1s=min(vplc1,qgfani(1.d0/xpomr1,bb1,1.d0-(1.d0-vvx)
     * *(1.d0-vvxt)*exp(-vtf-vtfq),vvxp,0.d0,icdp,icz,17)) !cut soft P (1P-end)
       vplcs=min(vplc1s,qgfani(1.d0/xpomr1,bb1,1.d0-(1.d0-vvx)
     * *(1.d0-vvxt)*exp(-vtf-vtfq),vvxp,0.d0,icdp,icz,16)) !single cut soft P
       vtlc1s=min(vtlc1,qgfani(xpomr1*sy,bb2,1.d0-(1.d0-vvx)
     * *(1.d0-vvxp)*exp(-vpf-vpfq),vvxt,0.d0,icdt,2,17))   !cut soft P (1P-end)
       vtlcs=min(vtlc1s,qgfani(xpomr1*sy,bb2,1.d0-(1.d0-vvx)
     * *(1.d0-vvxp)*exp(-vpf-vpfq),vvxt,0.d0,icdt,2,16))   !single cut soft P

       do iqq=1,5
        if(iqq.eq.1)then                 !general Pomeron (s+g)
         dpx=(vplc*vtl+vtlc*vpl)*exp(-2.d0*vpf-2.d0*vtf
     *   -2.d0*vpfq-2.d0*vtfq)*((1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt))**2
     *   -vplc1*vtl-vtlc1*vpl
        elseif(iqq.eq.2)then             !soft Pomeron
         dpx=(vplcs*vtls+vtlcs*vpls)*exp(-2.d0*vpf-2.d0*vtf
     *   -2.d0*vpfq-2.d0*vtfq)*((1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt))**2
     *   -vplc1s*vtls-vtlc1s*vpls
        elseif(iqq.eq.3)then             !q_v Pomeron
         dpx=2.d0*((vplcs*vtlq+vtlcs*vplq)*exp(-2.d0*vpf-2.d0*vtf)
     *   *((1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt))**2
     *   -vplc1s*vtlq-vtlc1s*vplq)
        elseif(iqq.eq.4)then             !qg Pomeron
         dpx=2.d0*vplq*(vtlcs*exp(-2.d0*vpf-2.d0*vtf)
     *   *((1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt))**2-vtlc1s)
        elseif(iqq.eq.5)then             !gq Pomeron
         dpx=2.d0*vtlq*(vplcs*exp(-2.d0*vpf-2.d0*vtf)
     *   *((1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt))**2-vplc1s)
        else
         stop'qgpcut: iqq?!'
        endif
        v1p(iqq)=v1p(iqq)+a1(ib1)*a1(ib2)*a1(ix1)/z*rp*dpx
       enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      do iqq=1,5
       v1p(iqq)=v1p(iqq)/16.d0*log(sy/sgap**2)*(r3p*pi/.0389d0)/g3p**3
      enddo
      return 
      end

c------------------------------------------------------------------------
      double precision function qgpomi(sy,bb,vvx,vvxp0,vvxt0
     *,genhp0,genht0,icdp0,icdt0,icz,iqq0)
c-----------------------------------------------------------------------
c qgpomi - integrated  eikonal contributions, including screening corrections
c sy   - pomeron mass squared,
c bb   - impact parameter squared,
c vvx  = 1 - exp[-sum_{j<J} chi_targ(j) - sum_{i<I} chi_proj(i)]
c vvxp = 1 - exp[-sum_{i>I} chi_proj(i)] 
c vvxt = 1 - exp[-sum_{j>J} chi_targ(j)] 
c icdp - diffractive state for the projectile,
c icdt - diffractive state for the target,
c icz  - projectile class
c iqq=0 - soft Pomeron,
c iqq=1 - 3p-contributions (q_v not involved),
c iqq=2 - 3p-contributions (q_v involved),
c iqq=3 - screening corrections to single cut Pomeron (s+g),
c iqq=4 - screening corrections to single cut soft Pomeron,
c iqq=5 - screening corrections to single cut q_v-Pomeron (qg+gq)
c iqq=6 - screening corrections to single cut qg-Pomeron
c iqq=7 - screening corrections to single cut gq-Pomeron
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wz(3),wi(3),wj(3),wm(3),wn(3)
      parameter(nfock=3)
      common /qgarr10/ am(6)
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr24/ qpomr(11,11,288,nfock**2,77)
     *,dhteik(11,11,216,nfock**2,60)
      common /qgarr43/ moniou
      common /qgdebug/ debug
      
      if(genhp0.ne.1.d0.and.icz.ne.2)stop'not pion-nucleus'
      if(genhp0.ne.1.d0.and.genht0.eq.1.d0)then
       genhp=1.d0
       genht=max(1.d0,genhp0)
       vvxp=vvxt0
       vvxt=vvxp0
       icdp=icdt0
       icdt=icdp0
       if(iqq0.gt.5)then
        iqq=13-iqq0
       else
        iqq=iqq0
       endif
      else
       genhp=max(1.d0,genhp0)
       genht=max(1.d0,genht0)
       vvxp=vvxp0
       vvxt=vvxt0
       icdp=icdp0
       icdt=icdt0
       iqq=iqq0
      endif
      
      qgpomi=0.d0
      rp=(rq(icdp,icz)+rq(icdt,2)+alfp*log(sy))*4.d0*.0389d0
      z=exp(-bb/rp)
      vsoft=fp(icdp,icz)*fp(icdt,2)*sy**dels*z/rp*4.d0*.0389d0
      if(iqq.eq.0)goto 1
      
      yl=dlog10((sy-am(2)**2-am(icz)**2)/2.d0/am(2))
      k=max(1,int(yl))
      k=min(k,9)     
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)

      z2=.2d0*exp(-12.5d0)
      z3=z2*exp(2.5d0)
      if(z.gt.z2.or.vvx+vvxp+vvxt.ne.0.d0.or.genhp*genht.ne.1.d0)then
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=(-bb/rp-dlog(0.2d0))/2.5d0+7.d0
       endif
       jz=min(9,int(zz))
       jz=max(2,jz)
       if(jz.eq.6)jz=5
       wz(2)=zz-jz
       wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
       wz(1)=1.d0-wz(2)+wz(3)
       wz(2)=wz(2)-2.d0*wz(3)
      else
       jz=1
       wz(1)=(z-z2)*(z-z3)/z2/z3
       wz(2)=z*(z-z3)/z2/(z2-z3)
       wz(3)=z*(z-z2)/z3/(z3-z2)
      endif

      ml=icdp+nfock*(icdt-1)
      if(vvx+vvxp+vvxt.eq.0.d0.and.genhp*genht.eq.1.d0)then
       ivv=1+18*(icz-1)*(10-3*icz)
       do l1=1,3
        l2=jz+l1-1
       do k1=1,3
        k2=k+k1-1
        qgpomi=qgpomi+qpomr(k2,l2,ivv,ml,iqq)*wk(k1)*wz(l1)
       enddo
       enddo
      else
       vl=max(1.d0,vvx*5.d0+1.d0)
       j=min(int(vl),4)    
       wj(2)=vl-j
       wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
       wj(1)=1.d0-wj(2)+wj(3)
       wj(2)=wj(2)-2.d0*wj(3)

       if(vvxp.eq.0.d0)then
        i=1
        i1max=1
        wi(1)=1.d0
       else
        if(icz.ne.2)stop'qgpomi: nuclear screening not for pions/kaons'
        vl1=max(1.d0,vvxp*5.d0+1.d0)
        i=min(int(vl1),4)    
        wi(2)=vl1-i
        wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
        wi(1)=1.d0-wi(2)+wi(3)
        wi(2)=wi(2)-2.d0*wi(3)
        i1max=3
       endif

       vl2=max(1.d0,vvxt*5.d0+1.d0)
       m=min(int(vl2),4)    
       wm(2)=vl2-m
       wm(3)=wm(2)*(wm(2)-1.d0)*.5d0
       wm(1)=1.d0-wm(2)+wm(3)
       wm(2)=wm(2)-2.d0*wm(3)

       if(genht.eq.1.d0)then
        do m1=1,3
         m2=m+m1-2
        do i1=1,i1max
         i2=i+i1-1
        do j1=1,3
         j2=j+j1-1
         ivv=j2+6*m2+18*(icz-1)*(4-icz+2*(3-icz)*i2)
        do l1=1,3
         l2=jz+l1-1
        do k1=1,3
         k2=k+k1-1
         qgpomi=qgpomi+qpomr(k2,l2,ivv,ml,iqq)
     *   *wk(k1)*wz(l1)*wj(j1)*wi(i1)*wm(m1)
        enddo
        enddo
        enddo
        enddo
        enddo
       else
        if(genht.le.1.5d0)then
         yg=(genht-1.d0)*10.d0+1.d0
        else
         yg=dlog(genht/1.5d0)*2.d0+6.d0
        endif
        n=min(9,int(yg))
        n=max(1,n)
        if(n.eq.5)n=4
        wn(2)=yg-n
        wn(3)=wn(2)*(wn(2)-1.d0)*.5d0
        wn(1)=1.d0-wn(2)+wn(3)
        wn(2)=wn(2)-2.d0*wn(3)
        do n1=1,3
         n2=n+n1-2
        do m1=1,3
         m2=m+m1-2
        do i1=1,i1max
         i2=i+i1-1
        do j1=1,3
         j2=j+j1-1
         ivv=j2+6*m2+18*(icz-1)*(4-icz+2*(3-icz)*i2)
        do l1=1,3
         l2=jz+l1-1
        do k1=1,3
         k2=k+k1-1
         qgpomi=qgpomi+qpomr(k2,l2,ivv,ml,iqq+7*n2)
     *   *wk(k1)*wz(l1)*wj(j1)*wi(i1)*wm(m1)*wn(n1)
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
       endif
      endif
      if(z.lt.z2/2.d0.and.(vvx+vvxp+vvxt.ne.0.d0
     *.or.genhp*genht.ne.1.d0))qgpomi=min(0.d0,qgpomi)

1     qgpomi=exp(qgpomi)
      if(iqq.eq.2.or.iqq.ge.5)qgpomi=max(0.d0,qgpomi-1.d0)
      qgpomi=qgpomi*vsoft
      
      if(genhp.ne.1.d0.and.iqq.ne.0.and.iqq.ne.6)then
       dht0=0.d0
       dht=0.d0
       if(genhp.le.1.5d0)then
        yg=(genhp-1.d0)*10.d0+1.d0
       else
        yg=dlog(genhp/1.5d0)*2.d0+6.d0
       endif
       n=min(9,int(yg))
       n=max(1,n)
       if(n.eq.5)n=4
       wn(2)=yg-n
       wn(3)=wn(2)*(wn(2)-1.d0)*.5d0
       wn(1)=1.d0-wn(2)+wn(3)
       wn(2)=wn(2)-2.d0*wn(3)
       iqd=min(5,iqq)
       do j1=1,3
        j2=j+j1-1
       do i1=1,i1max
        i2=i+i1-2
       do m1=1,3
        m2=m+m1-2
        ivv=j2+6*i2+36*m2
       do l1=1,3
        l2=jz+l1-1
       do k1=1,3
        k2=k+k1-1
        dht0=dht0+dhteik(k2,l2,ivv,ml,iqd)
     *  *wk(k1)*wz(l1)*wj(j1)*wi(i1)*wm(m1)
        do n1=1,3
         n2=n+n1-1
         dht=dht+dhteik(k2,l2,ivv,ml,iqd+5*n2)
     *   *wk(k1)*wz(l1)*wj(j1)*wi(i1)*wm(m1)*wn(n1)
        enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       dht=max(0.d0,exp(dht)-1.d0)*exp(dht0)
       if(iqq.ne.4)then
        qgpomi=max(0.d0,qgpomi-dht)
       else
        qgpomi=qgpomi+dht
       endif
      endif
      return 
      end 
  
c------------------------------------------------------------------------
      double precision function qgppdi(xp,iqq)
c-----------------------------------------------------------------------
c qgppdi - parton distributions in the Pomeron (*xp^dels)
c xp    - parton LC momentum share,
c iqq=0 - gluon
c iqq=1 - sea quark
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr43/ moniou
      common /qgdebug/ debug   
      
      if(iqq.lt.0.or.iqq.gt.1)stop'qgppdi:iqq out of range!!!'   
      if(debug.ge.3)write (moniou,201)xp,iqq
      if(xp.ge..9999999d0)then
       qgppdi=0.d0
      else
       if(iqq.eq.0)then                             !gluon
        qgppdi=(1.d0-xp)**beth(1)*(1.d0+xp)**bbbpom*(1.d0-dgqq)
       elseif(iqq.eq.1)then                         !quark
        qgppdi=qgftlf(xp,beth(1),bbbpom)*dgqq
       endif
      endif
      if(debug.ge.4)write (moniou,202)qgppdi
     
201   format(2x,'qgppdi - parton distr. in the Pomeron:'
     */4x,'xp=',e10.3,2x,'iqq=',i1)
202   format(2x,'qgppdi=',e10.3)
      return 
      end      
       
c------------------------------------------------------------------------
      double precision function qgppdc(xp,iqq,icdp,icz)
c-----------------------------------------------------------------------
c qgppdc - parton distr. in the Pomeron for different Fock states (*xp^dels)
c xp    - parton LC momentum share,
c iqq=0 - gluon
c iqq=1 - sea quark
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr43/ moniou
      common /qgdebug/ debug   
      
      if(iqq.lt.0.or.iqq.gt.1)stop'qgppdc:iqq out of range!!!'   
      if(debug.ge.3)write (moniou,201)xp,iqq
      if(xp.ge..9999999d0)then
       qgppdc=0.d0
      else
       if(iqq.eq.0)then                             !gluon
        qgppdc=(1.d0-xp)**beth(icz)*(1.d0+xp)**bbbi(icdp,icz)
     *  *(1.d0-dgqq)
       elseif(iqq.eq.1)then                         !quark
        qgppdc=qgftlf(xp,beth(icz),bbbi(icdp,icz))*dgqq
       endif
      endif
      if(debug.ge.4)write (moniou,202)qgppdc
     
201   format(2x,'qgppdi - parton distr. in the Pomeron:'
     */4x,'xp=',e10.3,2x,'iqq=',i1)
202   format(2x,'qgppdi=',e10.3)
      return 
      end      

c=============================================================================
      double precision function qgpdfb0(x,bb,icdp,icz,iqq)
c-----------------------------------------------------------------------------
c qgpdfb0 - generalized parton momentum distributions (without absorption) 
c [xf(x,b,qt0), GeV^2]
c x    - Feinman x,
c bb   - impact parameter squared,
c icdp - diffractive eigenstate,
c icz  - hadron class
c iqq=1 - g,
c iqq=2 - q_sea,
c iqq=3 - q_v
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(iqq.le.2)then
       rp=(rq(icdp,icz)-alfp*dlog(x))*4.d0*.0389d0
       qgpdfb0=rr*fp(icdp,icz)*qgppdc(x,iqq-1,icdp,icz)/x**dels
     * /rp*exp(-bb/rp)*4.d0*.0389d0
      elseif(iqq.eq.3)then
       rh=rq(icdp,icz)*4.d0*.0389d0
       qgpdfb0=qgvpdf(x,icz)/rh*exp(-bb/rh)*.0389d0/pi
      else
       stop'qgpdfb0: iqq?!'
      endif
      return
      end
      
c=============================================================================
      double precision function qgvpdf(x,icz)
c-----------------------------------------------------------------------------
c qgvpdf - valence quark structure function
c x   - Feinman x,
c icz - hadron class
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr25/ ahv(3)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      qgvpdf=(qggrv(x,qt0,icz,1)+qggrv(x,qt0,icz,2))*(1.d0-x)**ahv(icz)
      return
      end
      
c=============================================================================
      double precision function qgpdfb(x,bb,vvx,vvxp,icdp,icz,iqq,jj)
c-----------------------------------------------------------------------------
c qgpdfb - reaction-dependent generalized parton momentum distributions 
c [xf(x,b,qt0), GeV^2]
c x    - Feinman x,
c bb   - impact parameter squared,
c vvx  - screening factor,
c icdp - diffractive eigenstate,
c icz  - hadron class
c iqq=1 - g,
c iqq=2 - q,
c jj=1  - uncut,
c jj=2  - 1-P cut
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      qgpdfb=0.d0
      if(x*sgap.lt.1.d0)then      
       do ix=1,7
       do mx=1,2
        xpomr=(x*sgap)**(.5d0+x1(ix)*(mx-1.5d0))/sgap
        rp=(rq(icdp,icz)-alfp*log(xpomr))*4.d0*.0389d0
        rp1=alfp*dlog(xpomr/x)*4.d0*.0389d0
        rp2=rp1*rp/(rp1+rp)
        do ix2=1,7
        do mx2=1,2
         zz=.5d0+x1(ix2)*(mx2-1.5d0)
         bb0=-rp2*dlog(zz)
        do ix3=1,7
        do mx3=1,2
         phi=pi*(.5d0+x1(ix3)*(mx3-1.5d0))
         bb1=(dsqrt(bb)*rp1/(rp+rp1)-dsqrt(bb0)*cos(phi))**2
     *   +bb0*sin(phi)**2
         bb2=(dsqrt(bb)*rp/(rp+rp1)+dsqrt(bb0)*cos(phi))**2
     *   +bb0*sin(phi)**2
       
         vpf1=qgfani(1.d0/xpomr,bb2,vvx,0.d0,0.d0,icdp,icz,4) !soft fan (1P-end)
         vpf=min(vpf1,qgfani(1.d0/xpomr,bb2,vvx,0.d0,0.d0,icdp,icz,3)) !soft fan
         vi=qgppdi(x/xpomr,iqq-1)*(xpomr/x)**dels*exp(-bb1/rp1)
     *   *rr*g3p/rp1*4.d0*.0389d0
        
         if(jj.eq.1)then
          dpx=min(0.d0,1.d0-exp(-vpf)-vpf)*(1.d0-vvx)-vpf*vvx
     *    +(vpf-vpf1)
         else
          vpfc1=qgfani(1.d0/xpomr,bb2,vvx,vvxp,0.d0,icdp,icz,17) !cut P (1P-end)
          vpfc=min(vpfc1,qgfani(1.d0/xpomr,bb2,vvx,vvxp,0.d0
     *    ,icdp,icz,16))                           !cut fan (single cut Pomeron)
          dpx=vpfc*exp(-2.d0*vpf)*(1.d0-vvx)**2*(1.d0-vvxp)**2-vpfc1
         endif
         qgpdfb=qgpdfb+a1(ix)*a1(ix2)*a1(ix3)*dpx*vi/zz*rp2
        enddo
        enddo
        enddo
        enddo
       enddo
       enddo
       qgpdfb=-qgpdfb*dlog(x*sgap)*(r3p*pi/.0389d0)/g3p**3/8.d0 !GeV^2
      endif
      qgpdfb=(1.d0+qgpdfb/qgpdfb0(x,bb,icdp,icz,iqq))
      return
      end
      
c=============================================================================
      double precision function qgpdfm(x,icdp,icz,iqq)
c-----------------------------------------------------------------------------
c qgpdfm - minimum for GPDs
c x    - Feinman x,
c icdp - diffractive eigenstate,
c icz  - hadron class
c iqq=1 - g,
c iqq=2 - q_sea
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      qgpdfm=0.d0
      if(x*sgap.lt.1.d0)then      
       do ix=1,7
       do mx=1,2
        xpomr=(x*sgap)**(.5d0+x1(ix)*(mx-1.5d0))/sgap
         qgpdfm=qgpdfm+a1(ix)*qgppdi(x/xpomr,iqq-1)*qgfact(1.d0/xpomr)
       enddo
       enddo
       qgpdfm=qgpdfm*dlog(x*sgap)*dels/2.d0/qgppdc(x,iqq-1,icdp,icz)
      endif
      qgpdfm=1.d0+qgpdfm
      return
      end
      
c=============================================================================
      double precision function qgpdfbi(x,bb,vvx,vvxp,icdp,icz,iqq,jj)
c-----------------------------------------------------------------------------
c qgpdfb - reaction-dependent GPDs
c [xf(x,b,qt0), GeV^2]
c x   - Feinman x,
c bb  - impact parameter squared,
c vvx  - screening factor,
c icdp - diffractive eigenstate,
c icz - hadron class
c iqq=1 - g,
c iqq=2 - q_sea,
c iqq=3 - q_v,
c jj=1  - uncut,
c jj=2  - 1-P cut
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wz(3),wi(3),wj(3)
      parameter(nfock=3)
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr20/ spmax
      common /qgarr43/ moniou
      common /qgarr56/ pdfr(31,11,11,9,4*nfock)
      common /qgdebug/ debug

      qgpdfbi=0.d0
      if(x.ge.1.d0/sgap.or.iqq.eq.3)goto 1
      
      yl=-dlog(x*sgap)/dlog(spmax/sgap**2)*30.d0+1.d0
      k=max(1,int(yl-1.d0))
      k=min(k,29)
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)

      rp=(rq(icdp,icz)-alfp*dlog(x))*4.d0*.0389d0
      z=dexp(-bb/rp)
      z2=.2d0*exp(-12.5d0)
      z3=z2*exp(2.5d0)
      if(z.gt.z2)then
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=(-bb/rp-dlog(0.2d0))/2.5d0+7.d0
       endif
       jz=min(9,int(zz))
       jz=max(2,jz)
       if(jz.eq.6)jz=5
       wz(2)=zz-jz
       wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
       wz(1)=1.d0-wz(2)+wz(3)
       wz(2)=wz(2)-2.d0*wz(3)
      else
       jz=1
       wz(1)=(z-z2)*(z-z3)/z2/z3
       wz(2)=z*(z-z3)/z2/(z2-z3)
       wz(3)=z*(z-z2)/z3/(z3-z2)
      endif

      if(vvx+vvxp.eq.0.d0)then
       ii=icdp+nfock*(iqq-1)+2*nfock*(jj-1)
       iic=icz+2*(icz-1)*(3-icz)
       do l1=1,3
        l2=jz+l1-1
       do k1=1,3
        k2=k+k1-1
        qgpdfbi=qgpdfbi+pdfr(k2,l2,1,iic,ii)*wk(k1)*wz(l1)
       enddo
       enddo
      else
       if(vvx.eq.0.d0)then
        ivmax=1
        j=1
        wj(1)=1.d0
       else
        vl=max(1.d0,vvx*10.d0+1.d0)
        j=min(int(vl),9)    
        wj(2)=vl-dble(j)
        wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
        wj(1)=1.d0-wj(2)+wj(3)
        wj(2)=wj(2)-2.d0*wj(3)
        ivmax=3
       endif

       if(icz.ne.2.or.vvxp.eq.0.d0)then
        iv1max=1
        i=1
        wi(1)=1.d0
       else
        iv1max=2
        vl1=max(1.d0,vvxp*5.d0+1.d0)
        i=min(int(vl1),5)    
        wi(2)=vl1-i
        wi(1)=1.d0-wi(2)
       endif
       ii=icdp+nfock*(iqq-1)+2*nfock*(jj-1)
       do i1=1,iv1max
        i2=i+i1
        iic=icz+i2*(icz-1)*(3-icz)
       do j1=1,ivmax
        j2=j+j1-1
       do l1=1,3
        l2=jz+l1-1
       do k1=1,3
        k2=k+k1-1
        qgpdfbi=qgpdfbi+pdfr(k2,l2,j2,iic,ii)
     *  *wk(k1)*wz(l1)*wj(j1)*wi(i1)
       enddo
       enddo
       enddo
       enddo
      endif
      if(z.lt.z2/2.d0.and.vvx+vvxp.ne.0.d0)qgpdfbi=min(0.d0,qgpdfbi)

1     qgpdfbi=dexp(qgpdfbi)*qgpdfb0(x,bb,icdp,icz,iqq)
      return
      end

c=============================================================================
      double precision function qgpdfi(x,icz,iqq)
c-----------------------------------------------------------------------------
c qgpdfi - integrated pdfs at q0-scale [x f(x,qt0)]
c x   - Feinman x,
c icz - hadron class
c iqq=1 - g,
c iqq=2 - q
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208,nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      qgpdfi=0.d0
      do icdp=1,nfock     
       rp=(rq(icdp,icz)-alfp*log(x))*4.d0*.0389d0
       dpb=0.d0
       do ib=1,7
       do mb=1,2
        z=.5d0+x1(ib)*(mb-1.5d0)
        bb=-rp*log(z)
        pdf=qgpdfbi(x,bb,0.d0,0.d0,icdp,icz,iqq,1)
        if(iqq.eq.2)pdf=pdf+qgpdfbi(x,bb,0.d0,0.d0,icdp,icz,3,1)
        dpb=dpb+a1(ib)*pdf/z*rp
       enddo
       enddo
       qgpdfi=qgpdfi+cc(icdp,icz)*dpb*pi/.0389d0/2.d0
      enddo
      return
      end

c=============================================================================
      double precision function qgpdff(x,icdp,icz,iqq)
c-----------------------------------------------------------------------------
c qgpdff - integrated pdfs at q0-scale for given Fock state
c x   - Feinman x,
c icz - hadron class
c iqq=1 - g,
c iqq=2 - q
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208,nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      qgpdff=0.d0
      rp=(rq(icdp,icz)-alfp*log(x))*4.d0*.0389d0
      dpb=0.d0
      do ib=1,7
      do mb=1,2
       z=.5d0+x1(ib)*(mb-1.5d0)
       bb=-rp*log(z)
       pdf=qgpdfbi(x,bb,0.d0,0.d0,icdp,icz,iqq,1)
       if(iqq.eq.2)pdf=pdf+qgpdfbi(x,bb,0.d0,0.d0,icdp,icz,3,1)
       dpb=dpb+a1(ib)*pdf/z*rp
      enddo
      enddo
      qgpdff=dpb*pi/.0389d0/2.d0
      return
      end

c=============================================================================
      double precision function qgfsh(sy,bb,icdp,icdt,icz,iqq)
c-----------------------------------------------------------------------------
c qgfsh - semihard interaction eikonal
c sy   - pomeron mass squared,
c bb   - impact parameter squared,
c icdp - diffractive state for the projectile,
c icdt - diffractive state for the target,
c icz  - hadron class
c iqq  - type of the hard interaction (0-gg, 1-q_vg, 2-gq_v)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)sy,bb,iqq,icz

      qgfsh=0.d0      
      s2min=4.d0*fqscal*qt0
      xmin=s2min/sy
      if(xmin.ge.1.d0)return
      xmin=xmin**(delh-dels)
      if(iqq.eq.1)then
       icv=icz
       icq=2
       icdq=icdt
      elseif(iqq.eq.2)then
       icv=2
       icq=icz
       icdq=icdp
      endif
      if(debug.ge.3)write (moniou,205)xmin,iqq

c numerical integration over z1
      do i=1,7
      do m=1,2
       z1=(.5d0*(1.d0+xmin-(2*m-3)*x1(i)*(1.d0-xmin)))
     * **(1.d0/(delh-dels))
       ww=z1*sy
       sjqq=qgjit(qt0,qt0,ww,2,2)
       sjqg=qgjit(qt0,qt0,ww,1,2)
       sjgg=qgjit(qt0,qt0,ww,1,1)
       sjqqq=qgjit(qt0,qt0,ww,3,2) !inclusive qq cross-section (same flavor)
       sjqqa=qgjit(qt0,qt0,ww,4,2) !inclusive qbar-q cross-section (same flavor)
       if(debug.ge.3)write (moniou,203)ww,sjgg
     * ,sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0
        
       if(iqq.eq.0)then
        st2=0.d0
        do j=1,7
        do k=1,2
         xx=.5d0*(1.d0+x1(j)*(2*k-3))
         xp=z1**xx
         xm=z1/xp
         glu1=qgppdc(xp,0,icdp,icz)
         sea1=qgppdc(xp,1,icdp,icz)
         glu2=qgppdc(xm,0,icdt,2)
         sea2=qgppdc(xm,1,icdt,2)
         st2=st2+a1(j)*(glu1*glu2*sjgg+(glu1*sea2+glu2*sea1)*sjqg
     *   +sea1*sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0))
        enddo
        enddo
        rh=rq(icdp,icz)+rq(icdt,2)-alfp*dlog(z1)
        qgfsh=qgfsh-a1(i)*dlog(z1)/z1**delh*st2
     *  *exp(-bb/(4.d0*.0389d0*rh))/rh
      
       else
        st2=0.d0
        alh=.5d0+dels
        xam=z1**alh
        do j=1,7
        do k=1,2
         xp=(.5d0*(1.d0+xam+x1(j)*(2*k-3)*(1.d0-xam)))**(1.d0/alh)
         xm=z1/xp
         glu=qgppdc(xm,0,icdq,icq)
         sea=qgppdc(xm,1,icdq,icq)
         rh=rq(icdp,icz)+rq(icdt,2)-alfp*dlog(xm)
         fst=(glu*sjqg+sea*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0))
     *   *qgvpdf(xp,icv)/dsqrt(xp)*exp(-bb/(4.d0*.0389d0*rh))/rh
         st2=st2+a1(j)*fst
        enddo
        enddo
        st2=st2*(1.d0-xam)/alh
        qgfsh=qgfsh+a1(i)/z1**delh*st2
       endif
      enddo
      enddo
      if(iqq.eq.0)then
       qgfsh=qgfsh*rr**2*fp(icdp,icz)*fp(icdt,2)*factk
     * *(1.d0-xmin)/(delh-dels)/2.d0*pi
      else
       qgfsh=qgfsh*rr*fp(icdq,icq)*factk*(1.d0-xmin)/(delh-dels)/8.d0
      endif
      
      if(debug.ge.3)write (moniou,202)qgfsh
201   format(2x,'qgfsh - semihard interaction eikonal:'
     */4x,'sy=',e10.3,2x,'bb=',e10.3,2x,'iqq=',i1,2x,'icz=',i1)
202   format(2x,'qgfsh=',e10.3)
203   format(2x,'qgfsh:',2x,'sigma_hard(gg)=',e10.3
     *,2x,'sigma_hard(qq)=',e10.3)
205   format(2x,'qgfsh:',2x,'xmin=',e10.3,2x,'iqq=',i3)
      return
      end     
                                                
c------------------------------------------------------------------------
      double precision function qgftlf(zz,bett,bbb)
c-----------------------------------------------------------------------
c qgftlf - auxilliary function for semihard eikonals calculation
c zz - ratio of the quark and pomeron light cone x (zz=x_g/x_P)
c integration over quark to gluon light cone momentum ratio (z=x/x_g):
c qgftlf=int(dz) z^dels * (1-zz/z)^bett * (1+zz/z)^bbb * P_qg(z)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)zz
201   format(2x,'qgftlf:',2x,'zz=',e10.3)

      qgftlf=0.d0
      zmin=zz**(1.d0+dels)
      do i=1,7
      do m=1,2
       z=(.5d0*(1.d0+zmin+(2*m-3)*x1(i)*(1.d0-zmin)))
     * **(1.d0/(1.d0+dels))
       qgftlf=qgftlf+a1(i)*max(1.d-9,(1.d0-zz/z))**bett
     * *(1.d0+zz/z)**bbb*(z**2+(1.d0-z)**2)
       enddo
      enddo
      qgftlf=qgftlf*1.5d0*(1.d0-zmin)/(1.d0+dels)   !1.5=naflav/2 at Q0

      if(debug.ge.3)write (moniou,202)qgftlf
202   format(2x,'qgftlf=',e10.3)
      return
      end

c=============================================================================
      subroutine qgfz(b,gz,iddp1,iddp2)
c----------------------------------------------------------------------------
c hadron-hadron and hadron-nucleus cross sections calculation
c----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208,nfock=3)
      dimension gz(5),wt1(3),wt2(3)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)b,iddp1,iddp2
      do l=1,5
       gz(l)=0.d0
      enddo
      rp=(rq(1,icz)+rq(1,2)+alfp*log(scm))*4.d0*.0389d0

      do i1=1,7
      do m=1,2
       z=.5d0+x1(i1)*(m-1.5d0)
       bb1=rp*z
       bb2=rp*(1.d0-dlog(z))       
       do l=1,3
        wt1(l)=0.d0
        wt2(l)=0.d0
       enddo

       if(ia(2).eq.1)then
        do idd1=1,nfock
        do idd2=1,nfock
         vv1=exp(-qgpomi(scm,bb1,0.d0,0.d0,0.d0,1.d0,1.d0,idd1,idd2,icz
     *   ,1)-qgpomi(scm,bb1,0.d0,0.d0,0.d0,1.d0,1.d0,idd1,idd2,icz,2))
         vv2=exp(-qgpomi(scm,bb2,0.d0,0.d0,0.d0,1.d0,1.d0,idd1,idd2,icz
     *   ,1)-qgpomi(scm,bb2,0.d0,0.d0,0.d0,1.d0,1.d0,idd1,idd2,icz,2))
         do l=1,2
          wt1(l)=wt1(l)+cc(idd1,icz)*cc(idd2,2)*vv1**l
          wt2(l)=wt2(l)+cc(idd1,icz)*cc(idd2,2)*vv2**l
         enddo
         do idd3=1,nfock
          wt1(3)=wt1(3)+cc(idd1,icz)*cc(idd2,2)*cc(idd3,icz)*vv1
     *    *exp(-qgpomi(scm,bb1,0.d0,0.d0,0.d0,1.d0,1.d0,idd3,idd2,icz
     *    ,1)-qgpomi(scm,bb1,0.d0,0.d0,0.d0,1.d0,1.d0,idd3,idd2,icz,2))
          wt2(3)=wt2(3)+cc(idd1,icz)*cc(idd2,2)*cc(idd3,icz)*vv2
     *    *exp(-qgpomi(scm,bb2,0.d0,0.d0,0.d0,1.d0,1.d0,idd3,idd2,icz
     *    ,1)-qgpomi(scm,bb2,0.d0,0.d0,0.d0,1.d0,1.d0,idd3,idd2,icz,2))
         enddo
        enddo
        enddo
        do l=1,2
         gz(l)=gz(l)+a1(i1)*((1.d0-wt1(l))+(1.d0-wt2(l))/z)
        enddo
        gz(3)=gz(3)+a1(i1)*((wt1(2)-wt1(3))+(wt2(2)-wt2(3))/z)
        gz(4)=gz(4)+a1(i1)*((wt1(3)-wt1(1)**2)+(wt2(3)-wt2(1)**2)/z)
        gz(5)=gz(5)+a1(i1)*((1.d0-wt1(1))*bb1+(1.d0-wt2(1))/z*bb2)

       else
        do idd1=1,nfock
        do idd2=1,nfock
         vv1=exp(-qgpomi(scm,bb1,0.d0,0.d0,0.d0,1.d0,1.d0,iddp1,idd1
     *   ,icz,1)-qgpomi(scm,bb1,0.d0,0.d0,0.d0,1.d0,1.d0,iddp1,idd1,icz
     *   ,2)-qgpomi(scm,bb1,0.d0,0.d0,0.d0,1.d0,1.d0,iddp2,idd2,icz,1)
     *   -qgpomi(scm,bb1,0.d0,0.d0,0.d0,1.d0,1.d0,iddp2,idd2,icz,2))
         vv2=exp(-qgpomi(scm,bb2,0.d0,0.d0,0.d0,1.d0,1.d0,iddp1,idd1
     *   ,icz,1)-qgpomi(scm,bb2,0.d0,0.d0,0.d0,1.d0,1.d0,iddp1,idd1,icz
     *   ,2)-qgpomi(scm,bb2,0.d0,0.d0,0.d0,1.d0,1.d0,iddp2,idd2,icz,1)
     *   -qgpomi(scm,bb2,0.d0,0.d0,0.d0,1.d0,1.d0,iddp2,idd2,icz,2))

         if(idd1.eq.idd2)then
          wt1(1)=wt1(1)+cc(idd1,2)*vv1
          wt2(1)=wt2(1)+cc(idd1,2)*vv2
         endif
         wt1(2)=wt1(2)+cc(idd1,2)*cc(idd2,2)*vv1
         wt2(2)=wt2(2)+cc(idd1,2)*cc(idd2,2)*vv2
        enddo
        enddo
        cg1=qgrot(b,dsqrt(bb1))
        cg2=qgrot(b,dsqrt(bb2))
        do l=1,2
         gz(l)=gz(l)+a1(i1)*(cg1*(1.d0-wt1(l))+cg2*(1.d0-wt2(l))/z)
        enddo
       endif
      enddo
      enddo
      if(ia(2).eq.1.and.iddp1.eq.0.and.iddp2.eq.0)then     !hadron-proton
       do l=1,5
        gz(l)=gz(l)*pi*rp*10.d0          !in mb
       enddo
       gz(5)=gz(5)/gz(1)/2.d0/.0389d0
      else
       do l=1,2
        gz(l)=gz(l)*rp/2.d0              !in fm^2
       enddo
      endif

      if(debug.ge.2)write (moniou,203)gz
      if(debug.ge.3)write (moniou,202)
201   format(2x,'qgfz - hadronic cross-sections calculation'
     */4x,'b=',e10.3,2x,'iddp=',2i3)
202   format(2x,'qgfz - end')
203   format(2x,'qgfz: gz=',5e10.3)
      return
      end

c=============================================================================
      double precision function qghard(sy,bb,icdp,icdt,icz)
c-----------------------------------------------------------------------------
c qghard - hard quark-quark interaction cross-section
c s - energy squared for the interaction (hadron-hadron),
c icz - type of the primaty hadron (nucleon)
c----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr25/ ahv(3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)sy,icz

      qghard=0.d0
      s2min=4.d0*fqscal*qt0
      xmin=s2min/sy
      if(xmin.ge.1.d0)return
      xmin=xmin**(delh+.5d0)

      do i=1,7
      do m=1,2
       z1=(.5d0*(1.d0+xmin-(2*m-3)*x1(i)*(1.d0-xmin)))
     * **(1.d0/(delh+.5d0))
       sj=qgjit(qt0,qt0,z1*sy,2,2)
       sjqq=qgjit(qt0,qt0,z1*sy,3,2)
       sjqa=qgjit(qt0,qt0,z1*sy,4,2)
       if(debug.ge.3)write (moniou,203)z1*sy,sj

       st2=0.d0
       do j=1,7
       do k=1,2
        xx=.5d0*(1.d0+x1(j)*(2*k-3))
        xp=z1**xx
        xm=z1/xp
        ut=qggrv(xm,qt0,2,1)
        vt=qggrv(xm,qt0,2,2)
        up=qggrv(xp,qt0,icz,1)
        vp=qggrv(xp,qt0,icz,2)
        if(icz.eq.1)then
         sigg=sj*(up*dt+dp*ut)+sjqq*up*ut+sjqa*dp*dt
        elseif(icz.eq.2)then
         sigg=sj*(up*dt+dp*ut)+sjqq*(up*ut+dp*dt)
        elseif(icz.eq.3)then
         sigg=sj*dp*(ut+dt)+sjqq*up*ut+sj*up*dt
        else
         stop'qghard: wrong icz'
        endif
        st2=st2+a1(j)*sigg*(1.d0-xp)**ahv(icz)*(1.d0-xm)**ahv(2)
       enddo
       enddo
       qghard=qghard-a1(i)*st2/z1**(delh+.5d0)*dlog(z1)
      enddo
      enddo
      qghard=qghard*(1.d0-xmin)/(.5d0+delh)*.25d0*factk
      rh=rq(icdp,icz)+rq(icdt,2)
      qghard=qghard/(8.d0*pi*rh)*exp(-bb/(4.d0*.0389d0*rh))
      
      if(debug.ge.2)write (moniou,202)qghard
201   format(2x,'qghard - hard quark-quark interaction eikonal:'
     */2x,'s=',e10.3,2x,'icz=',i1)
202   format(2x,'qghard=',e10.3)
203   format(2x,'qghard:',2x,'s_hard=',e10.3)
      return
      end

c=============================================================================
      subroutine qgbdef(bba,bbb,xxa,yya,xxb,yyb,xxp,yyp,jb)
c-----------------------------------------------------------------------
c qgbdef - defines coordinates (xxp,yyp) of a multi-pomeron vertex
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      
      xx=xxa-xxb
      yy=yya-yyb
      bb=xx**2+yy**2
      if(bb.lt.1.d-5)then
       xxp=xxb+dsqrt(bba)
       yyp=yyb
      elseif(abs(yy).lt.1.d-8)then
       xxp=(bba-bbb+xxb**2-xxa**2)/2.d0/(xxb-xxa)
       yyp=yyb+(2*jb-3)*dsqrt(max(0.d0,bbb-(xxb-xxp)**2))
      else
       bbd=bb+bbb-bba
       discr=max(0.d0,4.d0*bb*bbb-bbd**2)
       xxp=(xx*bbd+(2*jb-3)*abs(yy)*dsqrt(discr))/2.d0/bb
       yyp=(bbd-2.d0*xx*xxp)/2.d0/yy
       xxp=xxp+xxb
       yyp=yyp+yyb
      endif      
      return
      end

c=============================================================================
      subroutine qgv(x,y,xb,vin,vdd,vabs)
c xxv - eikonal dependent factor for hadron-nucleus interaction
c (used for total and diffractive hadron-nucleus cross-sections calculation)
c----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208,nfock=3)
      dimension xb(iapmax,3),vabs(nfock)  !,fhard(3)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)x,y

      vin=0.d0
      vdd=0.d0
      do iddp1=1,nfock
       dv=0.d0
       do m=1,ia(2)
        bb=(x-xb(m,1))**2+(y-xb(m,2))**2
        dv=dv+qgpomi(scm,bb,0.d0,0.d0,0.d0,1.d0,1.d0,iddp1,iddt(m),icz
     *  ,1)+qgpomi(scm,bb,0.d0,0.d0,0.d0,1.d0,1.d0,iddp1,iddt(m),icz,2)
       enddo
       dv=exp(-dv)
       vabs(iddp1)=1.d0-dv**2       !1-exp(-2 * chi_i)
       vdd=vdd+cc(iddp1,icz)*dv**2  !sum_i cc(i) exp(-2 * chi_i)
       vin=vin+cc(iddp1,icz)*dv     !sum_i cc(i) exp(-chi_i)
      enddo
      vin=1.d0-vin**2               !1-sum_ij cc(i) cc(j) exp(-chi_i-chi_j)
      vdd=vdd+vin-1.d0 !sum_i cc(i)*exp(-2*chi_i)-sum_ij cc(i)*cc(j)*exp(-chi_i-chi_j)

      if(debug.ge.3)write (moniou,202)vin,vdd,vabs
201   format(2x,'qgv - eikonal factor: nucleon coordinates x='
     *  ,e10.3,2x,'y=',e10.3)
202   format(2x,'vin=',e10.3,2x,'vdd=',e10.3,2x,'vabs=',2e10.3)
      return
      end

c=============================================================================
      subroutine qgfdf(xxp,yyp,xpomr,vpac,vtac,vpht,vtht,vpacq,vtacq
     *,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)
c-----------------------------------------------------------------------
c qgfdf - configuration of fan contributions (cut and uncut fans)
c xxp, yyp -  coordinates of the multi-Pomeron vertex,
c xpomr    - LC momentum share of the multi-Pomeron vertex,
c ip       - proj. index,
c it       - targ. index
c vvx   = 1 - exp[-sum_{j<J} chi_targ(j) - sum_{i<I} chi_proj(i)]
c vvxp  = 1 - exp[-sum_{i>I} chi_proj(i)]
c vvxt  = 1 - exp[-sum_{j>J} chi_targ(j)]
c vvxpl = 1 - exp[-sum_{i<I} chi_proj(i)]
c vvxtl = 1 - exp[-sum_{j<J} chi_targ(j)]
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208,nfock=3)
      dimension vpac(iapmax),vtac(iapmax),vpacs(iapmax),vtacs(iapmax)
     *,vpacq(iapmax),vtacq(iapmax),vphtq(iapmax,2),vthtq(iapmax,2)
     *,vpht(iapmax,2),vtht(iapmax,2),fanht(2,2),fanhtm(2,2)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),bcoll
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr35/ iqvp(iapmax),iqvt(iapmax)
      common /qgarr43/ moniou
      common /qgarr46/ iconab(iapmax,iapmax),icona(iapmax)
     *,iconb(iapmax)
      common /qgdebug/ debug

      if(debug.ge.3)write (moniou,201)xxp,yyp,xpomr,ip,it
      vvx=0.d0
      vvxp=0.d0
      vvxt=0.d0
      vvxpl=0.d0
      vvxtl=0.d0
      genhp=1.d0
      genht=1.d0
      shmin=4.d0*fqscal*qt0
      if(scm.lt.1.001d0*sgap**2)return
      
      bbp0=(xa(ip,1)+bcoll-xxp)**2+(xa(ip,2)-yyp)**2
      bbt0=(xb(it,1)-xxp)**2+(xb(it,2)-yyp)**2
      sumup0=0.d0
      sumut0=0.d0
      pdfgp=0.d0
      pdfgt=0.d0
      gtotp=0.d0
      gtott=0.d0
                     
      nn=0
1     nn=nn+1 
      sumup=0.d0                      !proj. fans without targ. screening
      sumut=0.d0                      !targ. fans without proj. screening
      do ipp=1,ia(1)
       if(iconab(ipp,it).eq.0)then
        vpacs(ipp)=0.d0
       else
        bbp=(xa(ipp,1)+bcoll-xxp)**2+(xa(ipp,2)-yyp)**2
        vpacs(ipp)=qgfani(1.d0/xpomr,bbp,1.d0-exp(-sumup-sumut0)
     *  ,0.d0,0.d0,iddp(ipp),icz,3)
        if(nn.eq.1)then
         pdfgi=qgpdfbi(xpomr,bbp,1.d0-exp(-sumup),0.d0
     *   ,iddp(ipp),icz,1,1)
         if(ipp.eq.ip)pdfgp=pdfgi
         gtotp=gtotp+pdfgi
        endif
        sumup=sumup+vpacs(ipp)
       endif
      enddo

      do itt=1,ia(2)
       if(iconab(ip,itt).eq.0)then
        vtacs(itt)=0.d0
       else
        bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
        vtacs(itt)=qgfani(xpomr*scm,bbt,1.d0-exp(-sumut-sumup0)
     *  ,0.d0,0.d0,iddt(itt),2,3)
        if(nn.eq.1)then
         pdfgi=qgpdfbi(1.d0/xpomr/scm,bbt,1.d0-exp(-sumut),0.d0
     *   ,iddt(itt),2,1,1)
         if(itt.eq.it)pdfgt=pdfgi
         gtott=gtott+pdfgi
        endif
        sumut=sumut+vtacs(itt)
       endif
      enddo
      if(abs(sumup-sumup0)+abs(sumut-sumut0).gt..01d0.and.nn.lt.50)then
       sumup0=sumup
       sumut0=sumut
       goto 1
      endif
      
      if(ia(1).eq.1.or.pdfgp.eq.0.d0)then
       genhp=1.d0
      else
       genhp=max(1.d0,gtotp/pdfgp)
      endif
      if(ia(2).eq.1.or.pdfgt.eq.0.d0)then
       genht=1.d0
      else
       genht=max(1.d0,gtott/pdfgt)
      endif
      
      do ipp=1,ia(1)
       if(iconab(ipp,it).eq.0)then     !no connection (ipp far away from it)
        vpac(ipp)=0.d0
        vpacq(ipp)=0.d0
        do jj=1,2
         vpht(ipp,jj)=0.d0
         vphtq(ipp,jj)=0.d0
        enddo
       elseif(xpomr*shmin.ge.1.d0)then
        do jj=1,2
         vpht(ipp,jj)=0.d0
         vphtq(ipp,jj)=0.d0
        enddo
       else
        vvxs=sumut
        vvxps=0.d0
        do ipi=1,ia(1)
         if(ipi.lt.ipp)then
          vvxs=vvxs+vpacs(ipi)
         elseif(ipi.gt.ipp)then
          vvxps=vvxps+vpacs(ipi)
         endif
        enddo
        vvxs=1.d0-exp(-vvxs)
        vvxps=1.d0-exp(-vvxps)
        bbp=(xa(ipp,1)+bcoll-xxp)**2+(xa(ipp,2)-yyp)**2
        call qgfanht(fanht,scm,1.d0/xpomr,bbt0,bbp,vvxs,vvxps,genht
     *  ,iddt(it),iddp(ipp),2,icz)
        call qgfanhtm(fanhtm,scm,1.d0/xpomr/scm,bbt0,bbp,vvxs,vvxps
     *  ,genht,iddt(it),iddp(ipp),2,icz)
        do jj=1,2
         vpht(ipp,jj)=min(0.d0,fanht(1,jj)+fanhtm(1,jj))
         vphtq(ipp,jj)=min(0.d0,fanht(2,jj)+fanhtm(2,jj))
        enddo
       endif
      enddo

      do itt=1,ia(2)
       if(iconab(ip,itt).eq.0)then     !no connection
        vtac(itt)=0.d0
        vtacq(itt)=0.d0
        do jj=1,2
         vtht(itt,jj)=0.d0
         vthtq(itt,jj)=0.d0
        enddo
       elseif(xpomr*scm.le.shmin)then
        do jj=1,2
         vtht(itt,jj)=0.d0
         vthtq(itt,jj)=0.d0
        enddo
       else
        vvxs=sumup
        vvxts=0.d0
        do iti=1,ia(2)
         if(iti.lt.itt)then
          vvxs=vvxs+vtacs(iti)
         elseif(iti.gt.itt)then
          vvxts=vvxts+vtacs(iti)
         endif
        enddo
        vvxs=1.d0-exp(-vvxs)
        vvxts=1.d0-exp(-vvxts)
        bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
        call qgfanht(fanht,scm,xpomr*scm,bbp0,bbt,vvxs,vvxts,genhp
     *  ,iddp(ip),iddt(itt),icz,2)
        call qgfanhtm(fanhtm,scm,xpomr,bbp0,bbt,vvxs,vvxts,genhp
     *  ,iddp(ip),iddt(itt),icz,2)
        do jj=1,2
         vtht(itt,jj)=min(0.d0,fanht(1,jj)+fanhtm(1,jj))
         vthtq(itt,jj)=min(0.d0,fanht(2,jj)+fanhtm(2,jj))
        enddo 
       endif
      enddo
      
      sumup0=sumup
      sumut0=sumut
      nn=0
2     nn=nn+1 
      sumup=0.d0                      !proj. fans without targ. screening
      sumut=0.d0                      !targ. fans without proj. screening
      do ipp=1,ia(1)
       if(iconab(ipp,it).ne.0)then
        bbp=(xa(ipp,1)+bcoll-xxp)**2+(xa(ipp,2)-yyp)**2
        vps=qgfani(1.d0/xpomr,bbp,1.d0-exp(-sumup-sumut0)
     *  ,0.d0,0.d0,iddp(ipp),icz,3)
        vpac(ipp)=max(vps,vpht(ipp,1)+qgfani(1.d0/xpomr,bbp
     *  ,1.d0-exp(-sumup-sumut0),0.d0,0.d0,iddp(ipp),icz,1))
        sumup=sumup+vpac(ipp)
        if(iqvp(ipp).eq.0)then
         vpacq(ipp)=max(0.d0,vphtq(ipp,1)+qgfani(1.d0/xpomr,bbp
     *   ,1.d0-exp(-sumup-sumut0),0.d0,0.d0,iddp(ipp),icz,5))
        else
         vpacq(ipp)=0.d0
        endif
       endif
      enddo

      do itt=1,ia(2)
       if(iconab(ip,itt).ne.0)then
        bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
        vts=qgfani(xpomr*scm,bbt,1.d0-exp(-sumut-sumup0)
     *  ,0.d0,0.d0,iddt(itt),2,3)
        vtac(itt)=max(vts,vtht(itt,1)+qgfani(xpomr*scm,bbt
     *  ,1.d0-exp(-sumut-sumup0),0.d0,0.d0,iddt(itt),2,1))
        sumut=sumut+vtac(itt)
        if(iqvt(itt).eq.0)then
         vtacq(itt)=max(0.d0,vthtq(itt,1)+qgfani(xpomr*scm,bbt
     *  ,1.d0-exp(-sumut-sumup0),0.d0,0.d0,iddt(itt),2,5))
        else
         vtacq(itt)=0.d0
        endif
       endif
      enddo
      if(abs(sumup-sumup0)+abs(sumut-sumut0).gt..01d0.and.nn.lt.50)then
       sumup0=sumup
       sumut0=sumut
       goto 2
      endif
      
      if(ia(1).ne.1)then
       do ipp=1,ia(1)
        if(ipp.lt.ip)then
         vvxpl=vvxpl+vpac(ipp)
        elseif(ipp.gt.ip)then
         vvxp=vvxp+vpac(ipp)
        endif
       enddo
      endif
      if(ia(2).ne.1)then
       do itt=1,ia(2)
        if(itt.lt.it)then
         vvxtl=vvxtl+vtac(itt)
        elseif(itt.gt.it)then
         vvxt=vvxt+vtac(itt)
        endif
       enddo
      endif
      
      vvx=1.d0-exp(-vvxpl-vvxtl)
      vvxp=1.d0-exp(-vvxp)
      vvxpl=1.d0-exp(-vvxpl)
      vvxt=1.d0-exp(-vvxt)
      vvxtl=1.d0-exp(-vvxtl)
      if(debug.ge.4)write (moniou,202)
      
201   format(2x,'qgfdf - configuration of fan contributions:'
     */2x,'xxp=',e10.3,2x,'yyp=',e10.3,2x,'xpomr=',e10.3
     *,2x,'ip=',i3,2x,'it=',i3)
202   format(2x,'qgfdf - end')
      return 
      end

c=============================================================================
      subroutine qgconf
c-----------------------------------------------------------------------------
c interaction (cut Pomeron) configuration:
c b - impact parameter,
c xa(1-iap,3), xb(1-iat,3) - proj. and targ. nucleon coordinates,
c iddp(1-iap), iddt(1-iat) - proj. and targ. nucleon diffractive eigenstates,
c icona(1-iap) - connection for proj. nucleons (0 if too far from the target),
c iconab(1-iap,1-iat) - connection for proj.-targ. nucleons (0 if too far from
c each other),
c nwp, nwt - numbers of wounded proj. and targ. nucleons (inelastic or diff.),
c iwp(1-iap), iwt(1-iat) - indexes for wounded proj. and targ. nucleons
c (0 - intact, 1 - inel., 2,3 - diffr., -1 - recoiled from diffraction),
c ncola(1-iap), ncolb(1-iat) - index for inel.-wounded proj. and targ. nucleons,
c nbpom  - total number of Pomeron blocks,
c ias(k) (ibs(k)) - index of the proj. (targ.) nucleon for k-th Pomeron block,
c bbpom(k) - squared impact parameter (between proj. and targ.) for k-th block,
c nqs(k) - number of single Pomerons in k-th block (without cut 3P-vertexes),
c npompr(k) - number of proj. leg Pomerons in k-th block,
c npomtg(k) - number of targ. leg Pomerons in k-th block,
c npomin(k) - number of interm. Pomerons (between 2 3P-vertexes) in k-th block,
c xpopin(n,k) - LC momentum of the upper 3P-vertex for n-th interm. Pomeron
c in k-th block,
c xpomin(n,k) - LC momentum of the lower 3P-vertex for n-th interm. Pomeron
c in k-th block,
c nnpr(i,k) - proj. participant index for i-th single Pomeron in k-th block,
c nntg(i,k) - targ. participant index for i-th single Pomeron in k-th block,
c ilpr(i,k) - proj. index for i-th proj. leg Pomeron in k-th block,
c iltg(i,k) - proj. index for i-th targ. leg Pomeron in k-th block,
c lnpr(i,k) - proj. participant index for i-th proj. leg Pomeron in k-th block,
c lntg(i,k) - targ. participant index for i-th targ. leg Pomeron in k-th block,
c lqa(ip) - number of cut Pomerons connected to ip-th proj. nucleon (hadron),
c lqb(it) - number of cut Pomerons connected to it-th targ. nucleon (hadron),
c nbpi(n,ip) - block index for n-th Pomeron connected to ip-th proj. nucleon,
c nbti(n,it) - block index for n-th Pomeron connected to it-th targ. nucleon,
c idnpi(n,ip) - type of n-th Pomeron (0 - single, 1 - leg) connected to ip-th
c proj. nucleon,
c idnti(n,it) - type of n-th Pomeron (0 - single, 1 - leg) connected to it-th
c targ. nucleon,
c nppi(n,ip) - index in the block of n-th Pomeron connected to ip-th proj.
c nucleon (for single Pomerons),
c npti(n,it) - index in the block of n-th Pomeron connected to it-th targ.
c nucleon (for single Pomerons),
c nlpi(n,ip) - index in the block of n-th Pomeron connected to ip-th proj.
c nucleon (for leg Pomerons),
c nlti(n,it) - index in the block of n-th Pomeron connected to it-th targ.
c nucleon (for leg Pomerons),
c iprcn(ip) - index of the recoiled targ. nucleon for ip-th proj. nucleon
c (undergoing diffraction),
c itgcn(it) - index of the recoiled proj. nucleon for it-th targ. nucleon
c (undergoing diffraction),
c bpompr(n,ip) - squared impact parameter for n-th leg Pomeron connected
c to ip-th proj. nucleon,
c bpomtg(n,it) - squared impact parameter for n-th leg Pomeron connected
c to it-th targ. nucleon,
c xpompr(n,ip) - LC momentum of the 3P-vertex for n-th leg Pomeron connected
c to ip-th proj. nucleon,
c xpomtg(n,it) - LC momentum of the 3P-vertex for n-th leg Pomeron connected
c to it-th targ. nucleon
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208,npbmax=1000,npmax=200,legmax=900,nfock=3
     *,nptmax=95000)
      dimension xas(iapmax,3),vabsi(nfock,iapmax),wdifi(iapmax)
     *,vpac(iapmax),vtac(iapmax),vpht(iapmax,2),vtht(iapmax,2)
     *,vpacq(iapmax),vtacq(iapmax)
     *,bpomim(npmax),xpompi(legmax),xpomti(legmax),xpomip(npmax)
     *,xpomim(npmax),bpompi(legmax,2),bpomti(legmax,2),ipompi(legmax)
     *,ipomti(legmax),idpompi(legmax),idpomti(legmax),ncola(iapmax)
     *,ncolb(iapmax),wdp(nfock,iapmax),wdt(nfock,iapmax)
     *,wabs(nfock,nfock),xrapmin(100),xrapmax(100)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr4/  ey0(3)
      common /qgarr6/  pi,bm,amws
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),bcoll
      common /qgarr9/  iwp(iapmax),iwt(iapmax),lqa(iapmax),lqb(iapmax)
     *,iprcn(iapmax),itgcn(iapmax),ias(npbmax),ibs(npbmax),nqs(npbmax)
     *,npompr(npbmax),npomtg(npbmax),npomin(npbmax),nnpr(npmax,npbmax)
     *,nntg(npmax,npbmax),ilpr(legmax,npbmax),iltg(legmax,npbmax)
     *,lnpr(legmax,npbmax),lntg(legmax,npbmax),idpom(npmax,npbmax)
     *,idpomp(legmax,npbmax),idpomt(legmax,npbmax),nbpi(legmax,iapmax)
     *,nbti(legmax,iapmax),idnpi(legmax,iapmax),idnti(legmax,iapmax)
     *,nppi(legmax,iapmax),npti(legmax,iapmax),nlpi(legmax,iapmax)
     *,nlti(legmax,iapmax)
      common /qgarr10/ am(6)
      common /qgarr11/ b10
      common /qgarr12/ nsp
      common /qgarr13/ nsf,iaf(iapmax)
      common /qgarr14/ esp(4,nptmax),ich(nptmax)
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr22/ xppr0(iapmax),xmtg0(iapmax)
      common /qgarr23/ bbpom(npbmax),bpompr(legmax,iapmax,2)
     *,bpomtg(legmax,iapmax,2),xpompr(legmax,iapmax)
     *,xpomtg(legmax,iapmax),xpopin(npmax,npbmax),xpomin(npmax,npbmax)
     *,bpomin(npmax,npbmax)
      common /qgarr35/ iqvp(iapmax),iqvt(iapmax)
      common /qgarr43/ moniou
      common /qgarr46/ iconab(iapmax,iapmax),icona(iapmax)
     *,iconb(iapmax)
      common /qgarr55/ nwt,nwp           !N of wounded targ.(proj.) nucleons
      common /qgdebug/ debug
      common /qgsIIInex1/xan(iapmax,3),xbn(iapmax,3) !used to link with nexus
     *,bqgs,bmaxqgs,bmaxnex,bminnex
      common /jdiff/   jdiff             !diffr. type (external use)
ctp from epos
      integer ng1evt,ng2evt,ikoevt
      real    rglevt,sglevt,eglevt,fglevt,typevt
      common/c2evt/ng1evt,ng2evt,rglevt,sglevt,eglevt,fglevt,ikoevt
     *,typevt            !in epos.inc
      common /ebal/ ebal0(4),ebal(4)
      external qgran

      if(debug.ge.1)write (moniou,201)
      nsp=0
      nsf=0
      nsp0=nsp

c initialization
1     do i=1,ia(1)
       aks=qgran(b10)
       do ic=1,nfock
        aks=aks-cc(ic,icz)
        if(aks.lt.0.d0)goto 2
       enddo
       ic=nfock
2      iddp(i)=ic                          !diffractive eigenstates for proj.
      enddo
      do i=1,ia(2)
       aks=qgran(b10)
       do ic=1,nfock
        aks=aks-cc(ic,2)
        if(aks.lt.0.d0)goto 3
       enddo
       ic=nfock
3      iddt(i)=ic                          !diffractive eigenstates for targ.
      enddo

c-------------------------------------------------
c squared impact parameter is sampled uniformly (b**2<bm**2)
      bcoll=bm*dsqrt(qgran(b10))
      if(debug.ge.1)write (moniou,202)bcoll

c      if(bmaxnex.ge.0.d0)then              !used to link with nexus
c       b1=bminnex
c       b2=min(bm,bmaxnex)
c       if(b1.ge.b2)stop'bmin > bmax in qgsjet'
c       bcoll=dsqrt(b1*b1+(b2*b2-b1*b1)*qgran(b10))
c       bqgs=bcoll
c      endif
       bqgs=bcoll

c-------------------------------------------------
c nuclear configurations
      if(debug.ge.1)write (moniou,203)
      if(ia(1).gt.1)then          !projectile nucleon coordinates
       call qggea(ia(1),xa,1)     !xa(n,i), i=1,2,3 - x,y,z for n-th nucleon
      else
       do i=1,3
        xa(1,i)=0.d0              !projectile hadron
       enddo
      endif
      if(ia(2).gt.1)then          !target nucleon coordinates
       call qggea(ia(2),xb,2)     !xb(n,i), i=1,2,3 - x,y,z for n-th nucleon
      else
       do i=1,3
        xb(1,i)=0.d0              !target proton
       enddo
      endif

c-------------------------------------------------
c check connections
      if(debug.ge.1)write (moniou,204)
      do it=1,ia(2)
       iconb(it)=0
      enddo
      nconts=0
      do ip=1,ia(1)
       icdp=iddp(ip)
       icona(ip)=0
       do it=1,ia(2)
        icdt=iddt(it)
        bbp=(xa(ip,1)+bcoll-xb(it,1))**2+(xa(ip,2)-xb(it,2))**2
        vv1p=qgpomi(scm,bbp,0.d0,0.d0,0.d0,1.d0,1.d0,icdp,icdt,icz,1)
     *  +qgpomi(scm,bbp,0.d0,0.d0,0.d0,1.d0,1.d0,icdp,icdt,icz,2)
        if(vv1p.gt.1.d-3)then
         if(debug.ge.2)write (moniou,205)ip,it
         iconab(ip,it)=1
         icona(ip)=1
         iconb(it)=1
         nconts=nconts+1
         if(debug.ge.2)write (moniou,206)ip
         if(debug.ge.2)write (moniou,207)it
        else
         iconab(ip,it)=0
        endif
       enddo
      enddo
      if(nconts.eq.0)goto 1

      nrej=0
4     nrej=nrej+1
      if(debug.ge.2)write (moniou,208)nrej
      if(nrej.gt.100.or.nrej.gt.10.and.scm.gt.1.d3)then   !reject rare configs.
       if(debug.ge.1)write (moniou,209)
       goto 1
      endif
      nsp=nsp0
      nbpom=0
      nhard=0
      nwp=0
      nwt=0
      do i=1,ia(1)
       lqa(i)=0
       iwp(i)=0
       ncola(i)=0
       iqvp(i)=0
       do ic=1,nfock
        wdp(ic,i)=0.d0
       enddo
       xppr0(i)=1.d0
      enddo
      do i=1,ia(2)
       lqb(i)=0
       iwt(i)=0
       ncolb(i)=0
       iqvt(i)=0
       do ic=1,nfock
        wdt(ic,i)=0.d0
       enddo
       xmtg0(i)=1.d0
      enddo
      nqs(1)=0
      npomin(1)=0
      npompr(1)=0
      npomtg(1)=0

c-------------------------------------------------
c Pomeron configuration
      if(debug.ge.1)write (moniou,210)
      do 10 ip=1,ia(1)                   !loop over all projectile nucleons
       if(debug.ge.2)write (moniou,211)ip
       if(icona(ip).eq.0)goto 10
       x=xa(ip,1)+bcoll            !proj. x is shifted by the impact parameter b
       y=xa(ip,2)
       icdp=iddp(ip)                     !diffr. eigenstate for ip

       do 9 it=1,ia(2)                   !loop over all target nucleons
        if(debug.ge.2)write (moniou,212)it
        if(iconab(ip,it).eq.0)goto 9
        icdt=iddt(it)                    !diffr. eigenstate for it
        bbp=(x-xb(it,1))**2+(y-xb(it,2))**2 !distance squared between ip, it

c calculate nuclear screening factors for "middle point" -> eikonals
        xpomr=1.d0/dsqrt(scm)
        fff=(rq(icdt,2)+alfp*log(xpomr*scm))
     *  /(rq(icdp,icz)+rq(icdt,2)+alfp*log(scm))
        xxp=fff*x+(1.d0-fff)*xb(it,1)
        yyp=fff*y+(1.d0-fff)*xb(it,2)
        call qgfdf(xxp,yyp,xpomr,vpac,vtac,vpht,vtht,vpacq,vtacq
     *  ,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)
        vv=qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht,icdp,icdt,icz,1) !total eik.
        if(iqvp(ip)*iqvt(it).ne.0)then
         vvq=0.d0
        else
         vvq=qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht,icdp,icdt,icz,2) !q_v-eik.
        endif

        if(qgran(b10).gt.1.d0-exp(-2.d0*vv-2.d0*vvq))goto 9   !if no int.
        iwt(it)=1
        iwp(ip)=1
        ncola(ip)=ncola(ip)+1            !N of binary collisions for ip
        ncolb(it)=ncolb(it)+1            !N of binary collisions for it
        nbpom=nbpom+1                    !new Pomeron block
        if(nbpom.gt.npbmax)goto 4
        ias(nbpom)=ip         !proj. index for current elementary interaction
        ibs(nbpom)=it         !targ. index for current elementary interaction
        bbpom(nbpom)=bbp      !distance squared between ip, it
        if(debug.ge.2)write (moniou,214)nbpom,ip,it,n

        nqs(nbpom)=0
        npomin(nbpom)=0
        npompr(nbpom)=0
        npomtg(nbpom)=0
        vv1p=min(vv,qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht
     *  ,icdp,icdt,icz,3))               !1 cut P
        vv1ps=min(vv1p,qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht
     *  ,icdp,icdt,icz,4))               !soft P
        if(iqvp(ip)*iqvt(it).ne.0)then   !both q_v quarks already involved
         vvqc=0.d0
        else
         vvqc=min(vvq,qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht
     *   ,icdp,icdt,icz,5))              !1 cut q-P
        endif
        if(debug.ge.2)write (moniou,213)vv,vv1p
        if(qgran(b10).lt.2.d0*vvqc/(1.d0-exp(-2.d0*vv-2.d0*vvq)))then
         iqvc=1
         n=npgen(2.d0*vv,0,npmax)  !number of elem. inter. for (ip-it) collision
        else
         iqvc=0
         vvm=2.d0*vv+2.d0*max(0.d0,vvq-vvqc)
         if(vvm.lt..1d0)then
          n=npgen(vvm,1,npmax)     !number of elem. inter. for (ip-it) collision
         else
5         n=npgen(2.d0*vv,0,npmax)     !number of elem. inter. for (ip-it) coll.
          n=n+min(1,int(qgran(b10)+2.d0*max(0.d0,vvq-vvqc)))
          if(n.eq.0)goto 5
         endif
        endif
        
        do i=1-iqvc,n
         if(i.eq.0.or.qgran(b10).le.vv1p/vv.or.vv.lt.1.d-2
     *   .or.scm.lt.1.001d0*sgap**2)then !single cut P
          if(debug.ge.2)write (moniou,215)i
          np=nqs(nbpom)+1
          if(np.gt.npmax)goto 4
          nqs(nbpom)=np                  !update Pomeron number in the block
          l0=lqa(ip)+1
          if(l0.gt.legmax)goto 4
          lqa(ip)=l0                     !update number of connections for proj.
          nnpr(np,nbpom)=l0              !index for connected proj. participant
          nbpi(l0,ip)=nbpom
          idnpi(l0,ip)=0
          nppi(l0,ip)=np
          l0=lqb(it)+1
          if(l0.gt.legmax)goto 4
          lqb(it)=l0
          nntg(np,nbpom)=l0              !index for connected targ. participant
          nbti(l0,it)=nbpom
          idnti(l0,it)=0
          npti(l0,it)=np
          
          if(i.eq.0)then
           nhard=nhard+1
           if(nhard.gt.npmax)goto 4
           if(iqvt(it).ne.0)then
            idpom(np,nbpom)=2
            iqvp(ip)=1
           elseif(iqvp(ip).ne.0)then
            idpom(np,nbpom)=3
            iqvt(it)=1
           else
            vvqg=min(vvqc,qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht
     *      ,icdp,icdt,icz,6))!cut qg-P
            vvgq=min(vvqc,qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht
     *      ,icdp,icdt,icz,7))!cut gq-P
            aks=max(vvqc,vvqg+vvgq)*qgran(b10)
            if(aks.lt.vvqg)then
             idpom(np,nbpom)=2
             iqvp(ip)=1
            elseif(aks.lt.vvqg+vvgq)then
             idpom(np,nbpom)=3
             iqvt(it)=1
            else
             idpom(np,nbpom)=4
             iqvp(ip)=1
             iqvt(it)=1
            endif
           endif
          else
           if(qgran(b10).le.vv1ps/vv1p.or.vv.lt.1.d-2
     *     .or.scm.lt.1.001d0*sgap**2)then
            idpom(np,nbpom)=0
           else
            idpom(np,nbpom)=1
            nhard=nhard+1
            if(nhard.gt.npmax)goto 4
           endif
          endif

         else                           !multi-Pomeron vertex
          if(debug.ge.2)write (moniou,219)
          call qg3pdf(xpompi,xpomti,bpompi,bpomti,xpomip,xpomim,bpomim
     *    ,npompi,npomti,npin,ipompi,ipomti,idpompi,idpomti,wdp,wdt
     *    ,ip,it,iret)
          if(iret.ne.0)goto 4

          if(npompi.ne.0)then
           if(debug.ge.2)write (moniou,221)i,npompi
           lpr=0
6          lpr=lpr+1
           xmi=1.d0/xpompi(lpr)/scm
           ipp=ipompi(lpr)
           if(xmi.lt.xmtg0(it))then
            xmtg0(it)=xmtg0(it)-xmi
            np=npompr(nbpom)+1
            if(np.gt.legmax)goto 4
            npompr(nbpom)=np
            iwp(ipp)=1
            ilpr(np,nbpom)=ipp
            l0=lqa(ipp)+1
            if(l0.gt.legmax)goto 4
            lqa(ipp)=l0
            if(idpompi(lpr).eq.0)then
             idpomp(np,nbpom)=idpompi(lpr)
            else
             idpomp(np,nbpom)=idpompi(lpr)+4
             nhard=nhard+1
             if(nhard.gt.npmax)goto 4
            endif
            lnpr(np,nbpom)=l0
            nbpi(l0,ipp)=nbpom
            idnpi(l0,ipp)=1
            nlpi(l0,ipp)=np
            xpompr(l0,ipp)=xmi
            do l=1,2
             bpompr(l0,ipp,l)=bpompi(lpr,l)
            enddo
           else
            npompi=npompi-1
            if(npompi.eq.0.and.lqa(ipp).eq.0)goto 4
            if(lpr.le.npompi)then
             do l=lpr,npompi
              xpompi(l)=xpompi(l+1)
              ipompi(l)=ipompi(l+1)
              idpompi(l)=idpompi(l+1)
              do m=1,2
               bpompi(l,m)=bpompi(l+1,m)
              enddo
             enddo
            endif
            lpr=lpr-1
           endif
           if(lpr.lt.npompi)goto 6
          endif
          
          if(npomti.ne.0)then
           if(debug.ge.2)write (moniou,222)i,npomti
           ltg=0
7          ltg=ltg+1
           xpi=xpomti(ltg)
           itt=ipomti(ltg)
           if(xpi.lt.xppr0(ip))then
            xppr0(ip)=xppr0(ip)-xpi
            np=npomtg(nbpom)+1
            if(np.gt.legmax)goto 4
            npomtg(nbpom)=np
            iwt(itt)=1
            iltg(np,nbpom)=itt
            l0=lqb(itt)+1
            if(l0.gt.legmax)goto 4
            lqb(itt)=l0
            lntg(np,nbpom)=l0
            if(idpomti(ltg).eq.0)then
             idpomt(np,nbpom)=idpomti(ltg)
            else
             idpomt(np,nbpom)=idpomti(ltg)+6
             nhard=nhard+1
             if(nhard.gt.npmax)goto 4
            endif
            nbti(l0,itt)=nbpom
            idnti(l0,itt)=1
            nlti(l0,itt)=np
            xpomtg(l0,itt)=xpomti(ltg)
            do l=1,2
             bpomtg(l0,itt,l)=bpomti(ltg,l)
            enddo
           else
            npomti=npomti-1
            if(npomti.eq.0.and.lqb(itt).eq.0)goto 4
            if(ltg.le.npomti)then
             do l=ltg,npomti
              xpomti(l)=xpomti(l+1)
              ipomti(l)=ipomti(l+1)
              idpomti(l)=idpomti(l+1)
              do m=1,2
               bpomti(l,m)=bpomti(l+1,m)
              enddo
             enddo
            endif
            ltg=ltg-1
           endif
           if(ltg.lt.npomti)goto 7
          endif

          if(npin.ne.0)then
           if(debug.ge.2)write (moniou,220)i,npin
           lint=0
8          lint=lint+1
           xpi=xpomip(lint)
           xmi=xpomim(lint)
           if(xpi.lt.xppr0(ip).and.xmi.lt.xmtg0(it))then
            xppr0(ip)=xppr0(ip)-xpi
            xmtg0(it)=xmtg0(it)-xmi
            npomin(nbpom)=npomin(nbpom)+1
            if(npomin(nbpom).gt.npmax)goto 4
            l1=npomin(nbpom)
            xpopin(l1,nbpom)=xpi
            xpomin(l1,nbpom)=xmi
            bpomin(l1,nbpom)=bpomim(lint)
           else
            npin=npin-1
            if(lint.le.npin)then
             do l=lint,npin
              xpomip(l)=xpomip(l+1)
              xpomim(l)=xpomim(l+1)
              bpomim(l)=bpomim(l+1)
             enddo
            endif
            lint=lint-1
           endif
           if(lint.lt.npin)goto 8
          endif
         endif
        enddo                   !end of Pomeron loop
9      continue                 !end of it-loop
10    continue                  !end of ip-loop

c-------------------------------------------------
c low mass diffraction (hadron-hadron case)
      if(ia(1).eq.1.and.ia(2).eq.1.and.iwp(1).eq.0.and.iwt(1).eq.0)then
       wel=0.d0
       wnorm=0.d0
       do icdp=1,nfock
       do icdt=1,nfock
        vv=qgpomi(scm,bcoll**2,0.d0,0.d0,0.d0,1.d0,1.d0
     *  ,icdp,icdt,icz,1)  !total eik. (s+g)
        vvq=qgpomi(scm,bcoll**2,0.d0,0.d0,0.d0,1.d0,1.d0
     *  ,icdp,icdt,icz,2) !q_v-eikonal
        wabs(icdp,icdt)=exp(-vv-vvq)
        wnorm=wnorm+cc(icdp,icz)*cc(icdt,2)*wabs(icdp,icdt)**2
        wel=wel+cc(icdp,icz)*cc(icdt,2)*wabs(icdp,icdt)
       enddo
       enddo
       wnd=wel**2/wnorm
       if(qgran(b10).le.wnd)then
        if(debug.ge.1)write (moniou,231)
        goto 1
       endif
       
       wdtot=max(0.d0,wnorm-wel**2)
       wdifp=0.d0
       wdift=0.d0
       do icdp=1,nfock
       do icdt=1,nfock
       do ic=1,nfock
        wdifp=wdifp+cc(icdp,icz)*cc(icdt,2)*cc(ic,2)*wabs(icdp,icdt)
     *  *wabs(icdp,ic)
        wdift=wdift+cc(icdp,icz)*cc(icdt,2)*cc(ic,icz)*wabs(icdp,icdt)
     *  *wabs(ic,icdt)
       enddo
       enddo
       enddo
       wdifp=max(0.d0,wdifp-wel**2)
       wdift=max(0.d0,wdift-wel**2)
       aks=wdtot*qgran(b10)
       if(aks.lt.wdifp)then
        nwp=nwp+1
        iwp(1)=2
        iprcn(1)=1
        iwt(1)=-1
       elseif(aks.lt.wdifp+wdift)then
        nwt=nwt+1
        iwt(1)=2
        itgcn(1)=1
        iwp(1)=-1
       else
        nwp=nwp+1
        nwt=nwt+1
        iwp(1)=2
        iwt(1)=2
        iprcn(1)=1
        itgcn(1)=1
       endif
       goto 15
      endif

c-------------------------------------------------
c   diffraction (hadron-nucleus & nucleus-nucleus)
      do ip=1,ia(1)             !loop over all projectile nucleons
       x=xa(ip,1)+bcoll         !proj. x is shifted by b
       y=xa(ip,2)
       if(iwp(ip).eq.1)then
        nwp=nwp+1               !one more wounded proj. nucleon
        dptot=0.d0
        do i=1,nfock
         dptot=dptot+wdp(i,ip)
        enddo
        if(lqa(ip).eq.0.and.dptot.ne.0.d0)then
         icdps=iddp(ip)
         xpomr=1.d0/dsqrt(scm)
         do it=1,ia(2)
          if(iconab(ip,it).ne.0)then
           bbp=(x-xb(it,1))**2+(y-xb(it,2))**2
           icdt=iddt(it)
           do icdp=1,nfock
            iddp(ip)=icdp
            fff=(rq(icdt,2)+alfp*log(xpomr*scm))
     *      /(rq(icdp,icz)+rq(icdt,2)+alfp*log(scm))
            xxp=fff*x+(1.d0-fff)*xb(it,1)
            yyp=fff*y+(1.d0-fff)*xb(it,2)
            call qgfdf(xxp,yyp,xpomr,vpac,vtac,vpht,vtht,vpacq,vtacq
     *      ,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)
            vv=qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht,icdp,icdt,icz,1
     *      )+qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht,icdp,icdt,icz,2) !tot-eik
            wdp(icdp,ip)=wdp(icdp,ip)*exp(-vv)
           enddo
          endif
         enddo
         iddp(ip)=icdps
         wnorm=0.d0
         wel=0.d0
         do i=1,nfock
          wnorm=wnorm+cc(i,icz)*wdp(i,ip)**2
          wel=wel+cc(i,icz)*wdp(i,ip)
         enddo
         if(wnorm.gt.0.d0.and.qgran(b10).gt.wel**2/wnorm)iwp(ip)=3   !LMD excit.
        endif

       elseif(icona(ip).ne.0)then
        if(debug.ge.2)write (moniou,223)ip
        icdps=iddp(ip)
        xpomr=1.d0/dsqrt(scm)
        do it=1,ia(2)
         bbp=(x-xb(it,1))**2+(y-xb(it,2))**2
         icdt=iddt(it)
         do icdp=1,nfock
          if(iconab(ip,it).eq.0)then
           vabsi(icdp,it)=0.d0
          else
           iddp(ip)=icdp
           fff=(rq(icdt,2)+alfp*log(xpomr*scm))
     *     /(rq(icdp,icz)+rq(icdt,2)+alfp*log(scm))
           xxp=fff*x+(1.d0-fff)*xb(it,1)
           yyp=fff*y+(1.d0-fff)*xb(it,2)
           call qgfdf(xxp,yyp,xpomr,vpac,vtac,vpht,vtht,vpacq,vtacq
     *     ,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)
           vv=qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht,icdp,icdt,icz,1)
     *     +qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht,icdp,icdt,icz,2) !tot-eik.
           vabsi(icdp,it)=vv
          endif
         enddo
        enddo
        iddp(ip)=icdps
        wnorm=0.d0
        wel=0.d0
        do i=1,nfock
         xel=0.d0
         do it=1,ia(2)
          xel=xel+vabsi(i,it)
         enddo
         wnorm=wnorm+cc(i,icz)*exp(-2.d0*xel)
         wel=wel+cc(i,icz)*exp(-xel)
        enddo
        if(qgran(b10).gt.wel**2/wnorm)then !proj. diffr.
         wdift=0.d0
         do it=1,ia(2)
          if(iwt(it).ne.-1)then
           wnorm=0.d0
           wel=0.d0
           do i=1,nfock
            wnorm=wnorm+cc(i,icz)*exp(-2.d0*vabsi(i,it))
            wel=wel+cc(i,icz)*exp(-vabsi(i,it))
           enddo
           wdifi(it)=max(0.d0,1.d0-wel**2/wnorm)
           wdift=wdift+wdifi(it)
          else
           wdifi(it)=0.d0
          endif
         enddo
         if(wdift.gt.0.d0)then
          nwp=nwp+1
          iwp(ip)=2
          aks=qgran(b10)*wdift
          do it=1,ia(2)
           aks=aks-wdifi(it)
           if(aks.lt.0.d0)goto 11
          enddo
          stop'p-diffr: it>ia(2)'
11        iprcn(ip)=it
          if(iwt(it).eq.0)iwt(it)=-1
          if(debug.ge.2)write (moniou,224)ip,it
         endif
        endif
       endif
      enddo                     !end of ip-loop

      do 14 it=1,ia(2)          !check target diffraction
       if(iwt(it).eq.1)then
        nwt=nwt+1                         !one more wounded targ. nucleon
        dttot=0.d0
        do i=1,nfock
         dttot=dttot+wdt(i,it)
        enddo
        if(lqb(it).eq.0.and.dttot.ne.0.d0)then
         icdts=iddt(it)
         xpomr=1.d0/dsqrt(scm)
         do ip=1,ia(1)
          if(iconab(ip,it).ne.0)then
           bbp=(xa(ip,1)+bcoll-xb(it,1))**2+(xa(ip,2)-xb(it,2))**2
           icdp=iddp(ip)
           do icdt=1,nfock
            iddt(it)=icdt
            fff=(rq(icdt,2)+alfp*log(xpomr*scm))
     *      /(rq(icdp,icz)+rq(icdt,2)+alfp*log(scm))
            xxp=fff*(xa(ip,1)+bcoll)+(1.d0-fff)*xb(it,1)
            yyp=fff*xa(ip,2)+(1.d0-fff)*xb(it,2)
            call qgfdf(xxp,yyp,xpomr,vpac,vtac,vpht,vtht,vpacq,vtacq
     *      ,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)
            vv=qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht,icdp,icdt,icz,1
     *      )+qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht,icdp,icdt,icz,2) !tot-eik
            wdt(icdt,it)=wdt(icdt,it)*exp(-vv)
           enddo
          endif
         enddo
         iddt(it)=icdts
         wnorm=0.d0
         wel=0.d0
         do i=1,nfock
          wnorm=wnorm+cc(i,2)*wdt(i,it)**2
          wel=wel+cc(i,2)*wdt(i,it)
         enddo
         if(qgran(b10).gt.wel**2/wnorm)iwt(it)=3
        endif

       elseif(iconb(it).ne.0)then
        if(debug.ge.2)write (moniou,225)it
        icdts=iddt(it)
        xpomr=1.d0/dsqrt(scm)
        do ip=1,ia(1)
         bbp=(xa(ip,1)+bcoll-xb(it,1))**2+(xa(ip,2)-xb(it,2))**2
         icdp=iddp(ip)
         do icdt=1,nfock
          if(iconab(ip,it).eq.0)then
           vabsi(icdt,ip)=0.d0
          else
           iddt(it)=icdt
           fff=(rq(icdt,2)+alfp*log(xpomr*scm))
     *     /(rq(icdp,icz)+rq(icdt,2)+alfp*log(scm))
           xxp=fff*(xa(ip,1)+bcoll)+(1.d0-fff)*xb(it,1)
           yyp=fff*xa(ip,2)+(1.d0-fff)*xb(it,2)
           call qgfdf(xxp,yyp,xpomr,vpac,vtac,vpht,vtht,vpacq,vtacq
     *     ,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)
           vv=qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht,icdp,icdt,icz,1)
     *     +qgpomi(scm,bbp,vvx,vvxp,vvxt,genhp,genht,icdp,icdt,icz,2) !tot-eik.
           vabsi(icdt,ip)=vv
          endif
         enddo
        enddo
        iddt(it)=icdts
        wnorm=0.d0
        wel=0.d0
        do i=1,nfock
         xel=0.d0
         do ip=1,ia(1)
          xel=xel+vabsi(i,ip)
         enddo
         wnorm=wnorm+cc(i,2)*exp(-2.d0*xel)
         wel=wel+cc(i,2)*exp(-xel)
        enddo
        if(qgran(b10).gt.wel**2/wnorm)then !targ. diffr.
         wdift=0.d0
         do ip=1,ia(1)
          if(iwp(ip).eq.-1)then
           wdifi(ip)=0.d0
          else
           if(iwp(ip).eq.2)then
            itt=iprcn(ip)
            if(itt.eq.it)goto 13
            if(iwt(itt).eq.2)then
             wdifi(ip)=0.d0
             goto 12
            endif
           endif
           wnorm=0.d0
           wel=0.d0
           do i=1,nfock
            wnorm=wnorm+cc(i,2)*exp(-2.d0*vabsi(i,ip))
            wel=wel+cc(i,2)*exp(-vabsi(i,ip))
           enddo
           wdifi(ip)=max(0.d0,1.d0-wel**2/wnorm)
          endif
12        wdift=wdift+wdifi(ip)
         enddo
         if(wdift.eq.0.d0)goto 14
         aks=qgran(b10)*wdift
         do ip=1,ia(1)
          aks=aks-wdifi(ip)
          if(aks.lt.0.d0)goto 13
         enddo
         stop't-diffr: ip>ia(1)'
13       nwt=nwt+1
         iwt(it)=2
         itgcn(it)=ip
         if(debug.ge.2)write (moniou,226)it,ip
         if(iwp(ip).eq.0)then
          iwp(ip)=-1
         elseif(iwp(ip).eq.2)then
          itt=iprcn(ip)
          iprcn(ip)=it
          if(itt.ne.it.and.iwt(itt).eq.-1)iwt(itt)=0
         endif
        endif
       endif
14    continue
15    if(nwp+nwt.eq.0)then                !no interaction
       if(debug.ge.1)write (moniou,231)
       goto 1  
      endif

c check diffractive cross sections 
      jdiff=0                             !non-diffractive
      nqst=0
      nint=0
      if(nbpom.ne.0)then
       do i=1,nbpom
        nqst=nqst+nqs(i)
        nint=nint+npomin(i)
       enddo
      endif
      if(nqst.eq.0)then   !not elastic nor ND
       lqat=0
       do ip=1,ia(1)
        lqat=lqat+lqa(ip)
       enddo
       lqbt=0
       do it=1,ia(2)
        lqbt=lqbt+lqb(it)
       enddo
       iwpt=0
       do ip=1,ia(1)
        if(iwp(ip).eq.1)then
         iwpt=1
         goto 16
        elseif(iwp(ip).ge.2)then
         iwpt=2
        endif
       enddo
16     iwtt=0
       do it=1,ia(2)
        if(iwt(it).eq.1)then
         iwtt=1
         goto 17
        elseif(iwt(it).ge.2)then
         iwtt=2
        endif
       enddo
17     if(lqat.eq.0.and.lqbt.eq.0)then
        if(nbpom.eq.0.or.nint.eq.0)then
         if(iwpt.eq.2.and.iwtt.ne.2)then
          jdiff=6                         !SD(LM)-proj
         elseif(iwpt.ne.2.and.iwtt.eq.2)then
          jdiff=7                         !SD(LM)-targ
         elseif(iwpt.eq.2.and.iwtt.eq.2)then
          jdiff=8                         !DD(LM)
         else
          goto 1
         endif
        else
         if(iwpt.ne.2.and.iwtt.ne.2)then
          jdiff=9                         !CD(DPE)
         else
          jdiff=10                        !CD+LMD
         endif
        endif
       elseif(lqat.gt.0.and.lqbt.eq.0.and.iwtt.ne.2)then
        jdiff=1                          !SD(HM)-proj
       elseif(lqat.eq.0.and.lqbt.gt.0.and.iwpt.ne.2)then
        jdiff=2                          !SD(HM)-targ
       elseif(lqat.gt.0.and.lqbt.eq.0.and.iwtt.eq.2)then
        jdiff=3                          !DD(LHM)-proj
       elseif(lqat.eq.0.and.lqbt.gt.0.and.iwpt.eq.2)then
        jdiff=4                          !DD(LHM)-targ

       elseif(lqat.gt.0.and.lqbt.gt.0)then
        if(nbpom.eq.0)stop'problem with nbpom!!!'
        xrapmax(1)=1.d0
        xrapmin(1)=1.d0/scm
        do ibpom=1,nbpom
         if(npompr(ibpom).gt.0)then    !check all prog. legs
          do i=1,npompr(ibpom)
           ip=ilpr(i,ibpom)
           lpom=lnpr(i,ibpom)
           xrapmax(1)=min(xrapmax(1),1.d0/xpompr(lpom,ip)/scm)
          enddo
         endif
         if(npomtg(ibpom).gt.0)then    !check all targ. legs
          do i=1,npomtg(ibpom)
           it=iltg(i,ibpom)
           lpom=lntg(i,ibpom)
           xrapmin(1)=max(xrapmin(1),xpomtg(lpom,it))
          enddo
         endif
        enddo
        if(xrapmin(1).gt..999d0*xrapmax(1))goto 19
        nraps=1
18      if(nraps.gt.90)stop'nraps>90'
        do ibpom=1,nbpom
         if(npomin(ibpom).gt.0)then
          do i=1,npomin(ibpom)
           if(nraps.eq.1)then
            if(1.d0/scm/xpomin(i,ibpom).lt..999d0*xrapmax(1)
     *      .and.xpopin(i,ibpom).gt.1.001d0*xrapmin(1))then !rap-gaps changed
             if(1.d0/scm/xpomin(i,ibpom).lt.1.001d0*xrapmin(1)
     *       .and.xpopin(i,ibpom).gt..999d0*xrapmax(1))then !no rap-gap (filled)
               goto 19
             elseif(xpopin(i,ibpom).gt..999d0*xrapmax(1))then
              xrapmax(1)=1.d0/scm/xpomin(i,ibpom)
             elseif(1.d0/scm/xpomin(i,ibpom).lt.1.001d0*xrapmin(1))then
              xrapmin(1)=xpopin(i,ibpom)
             else
              xrapmin(2)=xrapmin(1)
              xrapmin(1)=xpopin(i,ibpom)
              xrapmax(2)=1.d0/scm/xpomin(i,ibpom)
              nraps=2
              goto 18
             endif
            endif
           else
            if(1.d0/scm/xpomin(i,ibpom).lt..999d0*xrapmax(1)
     *      .and.xpopin(i,ibpom).gt.1.001d0*xrapmin(nraps))then !rapgaps changed
             if(1.d0/scm/xpomin(i,ibpom).lt.1.001d0*xrapmin(nraps)
     *       .and.xpopin(i,ibpom).gt..999d0*xrapmax(1))then !no rap-gap (filled)
              goto 19
             else
              do irap=1,nraps
               if(xpopin(i,ibpom).gt..999d0*xrapmax(irap).and.1.d0/scm
     *         /xpomin(i,ibpom).lt.1.001d0*xrapmin(irap))then !gap filled
                if(irap.lt.nraps)then
                 do j=irap,nraps-1
                  xrapmax(j)=xrapmax(j+1)
                  xrapmin(j)=xrapmin(j+1)
                 enddo
                endif
                nraps=nraps-1
                goto 18
               elseif(xpopin(i,ibpom).gt..999d0*xrapmax(irap))then
                xrapmax(irap)=min(1.d0/scm/xpomin(i,ibpom)
     *          ,xrapmax(irap))
               elseif(1.d0/scm/xpomin(i,ibpom)
     *         .lt.1.001d0*xrapmin(irap))then
                xrapmin(irap)=max(xpopin(i,ibpom),xrapmin(irap))
               elseif(1.d0/scm/xpomin(i,ibpom).gt.xrapmin(irap)
     *         .and.xpopin(i,ibpom).lt.xrapmax(irap))then
                xrapmin(irap)=max(xpopin(i,ibpom),xrapmin(irap))
                if(irap.lt.nraps)then
                 do j=1,nraps-irap
                  xrapmax(nraps-j+2)=xrapmax(nraps-j+1)
                  xrapmin(nraps-j+2)=xrapmin(nraps-j+1)
                 enddo
                endif
                xrapmin(irap+1)=xrapmin(irap)
                xrapmax(irap+1)=1.d0/scm/xpomin(i,ibpom)
                nraps=nraps+1
                goto 18
               endif
              enddo                       !end of irap-loop
             endif
            endif
           endif
          enddo                           !end of npin-loop
         endif
        enddo                             !end of ibpom-loop
        jdiff=5                          !DD(HM)
       endif
      endif                              !end of diffr. check

ctp define collision type
19    typevt=0                      !no interaction
      if(jdiff.eq.0.or.jdiff.gt.10)then                     !ND (no rap-gaps)
       typevt=1
      elseif(jdiff.eq.8.or.jdiff.eq.10.or.jdiff.eq.3.or.jdiff.eq.4
     *.or.jdiff.eq.5)then                                   !DD + (CD+LMD)
       typevt=2
      elseif(jdiff.eq.1.or.jdiff.eq.6)then                  !SD pro
       typevt=4
      elseif(jdiff.eq.2.or.jdiff.eq.7)then                  !SD tar
       typevt=-4
      elseif(jdiff.eq.9)then                                !CD
       typevt=3
      else
       stop'problem with typevt!'
      endif

c form projectile spectator part
      if(debug.ge.1)write (moniou,227)
      nspec=0
      do ip=1,ia(1)
       if(iwp(ip).eq.0)then
        if(debug.ge.2)write (moniou,228)ip
        nspec=nspec+1
        do l=1,3
         xas(nspec,l)=xa(ip,l)
        enddo
       endif
      enddo
      nspect=0
      do it=1,ia(2)
       if(iwt(it).eq.0)nspect=nspect+1
      enddo

c inelastic interaction: energy sharing and particle production
      if(ia(1).eq.nspec.or.ia(2).eq.nspect)stop'ia(1)=nspec!!!'
      if(debug.ge.1)write (moniou,229)

      wp0=dsqrt(scm)*(1.d0+am(2)/(am(2)+wplab)
     **(dble(ia(2)-nspect)/dble(ia(1)-nspec)-1.d0))
      wm0=dsqrt(scm)+am(icz)**2/dsqrt(scm)*(1.d0+am(2)/wplab)
     **(dble(ia(1)-nspec)/dble(ia(2)-nspect)-1.d0)

      ebal0(1)=.5d0*(wp0*(ia(1)-nspec)+wm0*(ia(2)-nspect))
      ebal0(2)=.5d0*(wp0*(ia(1)-nspec)-wm0*(ia(2)-nspect))
      ebal0(3)=0.d0
      ebal0(4)=0.d0
       
      call qgsha(nbpom,ncola,ncolb,iret)
      if(iret.ne.0)goto 4
      if(nsp.le.nsp0+2)then
       if(debug.ge.1)write (moniou,230)nsp,lqa(1),lqb(1),nwp,nwt
       goto 4
      endif
      if(debug.ge.1)write (moniou,232)nsp

      if(icp.eq.1.or.icp.eq.2.or.icp.eq.4)then
       ibal=-2
      elseif(icp.eq.-1.or.icp.eq.-2.or.icp.eq.-4)then
       ibal=0
      else
       ibal=-1
      endif
      do i=1,nsp
       if(ich(i).eq.1.or.ich(i).eq.2.or.ich(i).eq.4.or.ich(i).eq.16
     * .or.ich(i).eq.-8.or.ich(i).eq.20.or.ich(i).eq.-21
     * .or.ich(i).eq.18)ibal=ibal+1
       if(ich(i).eq.-1.or.ich(i).eq.-2.or.ich(i).eq.-4.or.ich(i).eq.-16
     * .or.ich(i).eq.8.or.ich(i).eq.-20.or.ich(i).eq.21
     * .or.ich(i).eq.-18)ibal=ibal-1
       if(ich(i).eq.7)ibal=ibal+2
       if(ich(i).eq.-7)ibal=ibal-2
      enddo

      if(ibal.ne.0.and.ia(1)+ia(2).eq.2)then
       write(*,*)'charge balance:',ibal,icp,ia
       do i=1,nsp
        write(*,*)i,ich(i)
       enddo
      endif
      
      if(abs(ebal(1))/wp0.gt.1.d-6.or.abs(ebal(3)).gt.1.d-8)then
       write(*,*)'ebal-conf',ebal,nsp,nspec,nspect
      endif

c fragmentation of the projectile spectator part
      if(debug.ge.1)write (moniou,233)
      call qgfrgm(nspec,xas)
      if(debug.ge.1)write (moniou,234)nsf
      if(debug.ge.1)write (moniou,235)

201   format(2x,'qgconf - configuration of the interaction')
202   format(2x,'qgconf: impact parameter b=',e10.3,' fm')
203   format(2x,'qgconf: nuclear configurations')
204   format(2x,'qgconf: check connections')
205   format(2x,'qgconf: ',i3,'-th proj. nucleon may interact with '
     *,i3,'-th target nucleon')
206   format(2x,'qgconf: ',i3,'-th projectile nucleon may interact')
207   format(2x,'qgconf: ',i3,'-th target nucleon may interact')
208   format(2x,'qgconf: ',i3,'-th rejection,'
     *,' redo Pomeron configuration')
209   format(2x,'qgconf: too many rejections,'
     *,' redo nuclear configuartions')
210   format(2x,'qgconf: Pomeron configuration')
211   format(2x,'qgconf: check ',i3,'-th projectile nucleon')
212   format(2x,'qgconf: interaction with ',i3,'-th target nucleon?')
213   format(2x,'qgconf: eikonals - total: ',e10.3,2x,'single: ',e10.3)
214   format(2x,'qgconf: ',i4,'-th Pomeron block connected to ',i3
     *,'-th proj. nucleon and'/4x,i3,'-th targ. nucleon;'
     *,' number of element. processes in the block: ',i3)
215   format(2x,'qgconf: ',i3
     *,'-th process in the block is single cut Pomeron')
219   format(2x,'qgconf: configuration of multi-Pomeron vertexes')
220   format(2x,'qgconf: ',i3,'-th process in the block contains '
     *,i3,' interm. Pomerons')
221   format(2x,'qgconf: ',i3,'-th process in the block contains '
     *,i3,' proj. legs')
222   format(2x,'qgconf: ',i3,'-th process in the block contains '
     *,i3,' targ. legs')
223   format(2x,'qgconf: check diffraction for ',i3,'-th proj. nucleon')
224   format(2x,'qgconf: diffr. of ',i3,'-th proj. nucleon,'
     *,' recoil of ',i3,'-th targ. nucleon')
225   format(2x,'qgconf: check diffraction for ',i3,'-th targ. nucleon')
226   format(2x,'qgconf: diffr. of ',i3,'-th targ. nucleon,'
     *,' recoil of ',i3,'-th proj. nucleon')
227   format(2x,'qgconf: projectile spectator part')
228   format(2x,'qgconf: ',i3,'-th proj. nucleon stays idle')
229   format(2x,'qgconf: inelastic interaction: energy sharing'
     *,' and particle production')
230   format(2x,'qgconf: no particle produced - rejection',5i4)
231   format(2x,'qgconf: no interaction - rejection')
232   format(2x,'qgconf: ',i5,' particles have been produced')
233   format(2x,'qgconf: fragmentation of the proj. spectator part')
234   format(2x,'qgconf: ',i3,' proj. fragments have been produced')
235   format(2x,'qgconf - end')
      return
      end
      
c=============================================================================
      subroutine qg3pdf(xpompi,xpomti,bpompi,bpomti,xpomip,xpomim
     *,bpomim,nppr,nptg,npin,ipompi,ipomti,idpompi,idpomti,wdp,wdt
     *,ip,it,iret)
c-----------------------------------------------------------------------
c qg3pdf - configuration for multi-Pomeron/diffractive contributions
c ip,it - indexes of proj. and targ. nucleons for current collision
c to determine:
c nppr - number of proj. leg Pomerons in the process,
c nptg - number of targ. leg Pomerons in the process,
c npin - number of interm. Pomerons (between 2 3P-vertexes) in the process,
c xpomip(i) - LC momentum of the upper 3P-vertex for i-th interm. Pomeron
c in the process,
c xpomim(i) - LC momentum of the lower 3P-vertex for i-th interm. Pomeron
c in the process,
c ipompi(i) - proj. index for i-th proj. leg Pomeron in the process,
c ipomti(i) - proj. index for i-th targ. leg Pomeron in the process,
c bpompi(i) - squared impact param. for i-th proj. leg Pomeron in the process,
c bpomti(i) - squared impact param. for i-th targ. leg Pomeron in the process,
c xpompi(i) - LC momentum of the 3P-vertex for i-th proj. leg Pomeron,
c xpomti(i) - LC momentum of the 3P-vertex for i-th targ. leg Pomeron
c iret=1 - reject configuration
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208,npbmax=1000,npmax=200,legmax=900,levmax=20
     *,nfock=3)
      dimension vpac(iapmax),vtac(iapmax),vpht(iapmax,2),vtht(iapmax,2)
     *,vpacq(iapmax),vtacq(iapmax),vpac0(iapmax),vtac0(iapmax)
     *,vpact(iapmax),vtact(iapmax),xpomip(npmax),xpomim(npmax)
     *,bpomim(npmax),xpompi(legmax),xpomti(legmax)
     *,bpompi(legmax,2),bpomti(legmax,2),ipompi(legmax),ipomti(legmax)
     *,idpompi(legmax),idpomti(legmax),ippr0(legmax),iptg0(legmax)
     *,nppm(levmax),ippm(legmax,levmax),ii(levmax),xpomm(levmax)
     *,wgpm(levmax),xxm(levmax),yym(levmax),itypr0(legmax)
     *,itytg0(legmax),itypm(legmax,levmax),vv(18),wdp(nfock,iapmax)
     *,wdt(nfock,iapmax)
     *,xpmmin(npmax),xpmmax(npmax),bbm(npmax),ityppm(npmax),iorm(npmax)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),bcoll
      common /qgarr9/  iwp(iapmax),iwt(iapmax),lqa(iapmax),lqb(iapmax)
     *,iprcn(iapmax),itgcn(iapmax),ias(npbmax),ibs(npbmax),nqs(npbmax)
     *,npompr(npbmax),npomtg(npbmax),npomin(npbmax),nnpr(npmax,npbmax)
     *,nntg(npmax,npbmax),ilpr(legmax,npbmax),iltg(legmax,npbmax)
     *,lnpr(legmax,npbmax),lntg(legmax,npbmax),idpom(npmax,npbmax)
     *,idpomp(legmax,npbmax),idpomt(legmax,npbmax),nbpi(legmax,iapmax)
     *,nbti(legmax,iapmax),idnpi(legmax,iapmax),idnti(legmax,iapmax)
     *,nppi(legmax,iapmax),npti(legmax,iapmax),nlpi(legmax,iapmax)
     *,nlti(legmax,iapmax)
      common /qgarr11/ b10
      common /qgarr12/ nsp
      common /qgarr13/ nsf,iaf(iapmax)
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr23/ bbpom(npbmax),bpompr(legmax,iapmax,2)
     *,bpomtg(legmax,iapmax,2),xpompr(legmax,iapmax)
     *,xpomtg(legmax,iapmax),xpopin(npmax,npbmax),xpomin(npmax,npbmax)
     *,bpomin(npmax,npbmax)
      common /qgarr26/ factk,fqscal
      common /qgarr35/ iqvp(iapmax),iqvt(iapmax)
      common /qgarr43/ moniou
      common /qgarr46/ iconab(iapmax,iapmax),icona(iapmax)
     *,iconb(iapmax)
      common /qgdebug/ debug
      external qgran

      if(debug.ge.2)write (moniou,201)ip,it
      
      if(scm.le.sgap**2)stop'qg3pdf: scm<sgap**2!' 
      iret=0
      vpacng=0.d0
      vtacng=0.d0
      vpacpe=0.d0
      vtacpe=0.d0
      vimp=0.d0
      viuc=0.d0
      viuu=0.d0
      vip=0.d0
      vicc=0.d0
      vicu=0.d0
      s2min=4.d0*fqscal*qt0
      
c normalization of rejection function      
      xpomr=1.d0/dsqrt(scm)
      bpt=dsqrt((xa(ip,1)+bcoll-xb(it,1))**2+(xa(ip,2)-xb(it,2))**2)
      rp1=(rq(iddp(ip),icz)-alfp*dlog(xpomr))*4.d0*.0389d0  
      rp2=(rq(iddt(it),2)+alfp*dlog(xpomr*scm))*4.d0*.0389d0 
      rp0=rp1*rp2/(rp1+rp2)
      bbpr=(bpt*rp1/(rp1+rp2))**2 
      bbtg=(bpt*rp2/(rp1+rp2))**2 
      call qgbdef(bbpr,bbtg,xa(ip,1)+bcoll,xa(ip,2),xb(it,1),xb(it,2)
     *,xxp,yyp,1)

      call qgfdf(xxp,yyp,xpomr,vpac,vtac,vpht,vtht,vpacq,vtacq
     *,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)     
      vvxts=1.d0-(1.d0-vvx)*(1.d0-vvxt)*exp(-vtac(it))
      vplcps=qgfani(1.d0/xpomr,bbpr,vvxts,vvxp,vvxpl,iddp(ip),icz,16)
      vplcp=max(vplcps,vpht(ip,2)
     *+qgfani(1.d0/xpomr,bbpr,vvxts,vvxp,vvxpl,iddp(ip),icz,14))
      vplcpe=max(vplcp,vpht(ip,1)
     *+qgfani(1.d0/xpomr,bbpr,vvxts,vvxp,vvxpl,iddp(ip),icz,13))
      vplcng=max(vplcpe,vpht(ip,1)
     *+qgfani(1.d0/xpomr,bbpr,vvxts,vvxp,vvxpl,iddp(ip),icz,12))
      vplc0=max(vplcng,vpht(ip,1)
     *+qgfani(1.d0/xpomr,bbpr,vvxts,vvxp,vvxpl,iddp(ip),icz,11))
      vplc=max(vplc0,vpht(ip,1)
     *+qgfani(1.d0/xpomr,bbpr,vvxts,vvxp,vvxpl,iddp(ip),icz,10))

      vvxps=1.d0-(1.d0-vvx)*(1.d0-vvxp)*exp(-vpac(ip))
      vtlcps=qgfani(xpomr*scm,bbtg,vvxps,vvxt,vvxtl,iddt(it),2,16)
      vtlcp=max(vtlcps,vtht(it,2)
     *+qgfani(xpomr*scm,bbtg,vvxps,vvxt,vvxtl,iddt(it),2,14))
      vtlcpe=max(vtlcp,vtht(it,1)
     *+qgfani(xpomr*scm,bbtg,vvxps,vvxt,vvxtl,iddt(it),2,13))
      vtlcng=max(vtlcpe,vtht(it,1)
     *+qgfani(xpomr*scm,bbtg,vvxps,vvxt,vvxtl,iddt(it),2,12))
      vtlc0=max(vtlcng,vtht(it,1)
     *+qgfani(xpomr*scm,bbtg,vvxps,vvxt,vvxtl,iddt(it),2,11))
      vtlc=max(vtlc0,vtht(it,1)
     *+qgfani(xpomr*scm,bbtg,vvxps,vvxt,vvxtl,iddt(it),2,10))

      sumcp0=0.d0
      sumup=0.d0
      do i=1,ia(1)
       sumup=sumup+vpac(i)
      enddo
      vvxs=(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
      do i=1,ia(1)-ip+1
       ipp=ia(1)-i+1
       bbp=(xa(ipp,1)+bcoll-xxp)**2+(xa(ipp,2)-yyp)**2
       sumup=sumup-vpac(ipp)
       vps=qgfani(1.d0/xpomr,bbp,1.d0-vvxs*exp(-sumup)
     * ,1.d0-exp(-sumcp0),1.d0-exp(-sumup),iddp(ipp),icz,3)
       vpac0(ipp)=max(vps,vpht(ipp,1)+qgfani(1.d0/xpomr,bbp
     * ,1.d0-vvxs*exp(-sumup),1.d0-exp(-sumcp0),1.d0-exp(-sumup)
     * ,iddp(ipp),icz,6))
       vpac0(ipp)=min(vpac(ipp),vpac0(ipp))
       if(ipp.gt.ip)sumcp0=sumcp0+vpac0(ipp)
      enddo
      sumct0=0.d0
      sumut=0.d0
      do i=1,ia(2)
       sumut=sumut+vtac(i)
      enddo
      vvxs=(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
      do i=1,ia(2)-it+1
       itt=ia(2)-i+1
       bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
       sumut=sumut-vtac(itt)
       vts=qgfani(xpomr*scm,bbt,1.d0-vvxs*exp(-sumut)
     * ,1.d0-exp(-sumct0),1.d0-exp(-sumut),iddt(itt),2,3)
       vtac0(itt)=max(vts,vtht(itt,1)+qgfani(xpomr*scm,bbt
     * ,1.d0-vvxs*exp(-sumut),1.d0-exp(-sumct0),1.d0-exp(-sumut)
     * ,iddt(itt),2,6))
       vtac0(itt)=min(vtac(itt),vtac0(itt))
       if(itt.gt.it)sumct0=sumct0+vtac0(itt)
      enddo
      vvxp0=1.d0-exp(-sumcp0)
      vvxt0=1.d0-exp(-sumct0)
      
      sumvp=0.d0
      sumvpl=0.d0
      do ipp=1,ia(1)
       if(iqvp(ipp).eq.0)then
        sumvp=sumvp+vpacq(ipp)
        if(ipp.lt.ip)sumvpl=sumvpl+vpacq(ipp)
       else
        vpacq(ipp)=0.d0
       endif
      enddo
      sumvt=0.d0
      sumvtl=0.d0
      do itt=1,ia(2)
       if(iqvt(itt).eq.0)then
        sumvt=sumvt+vtacq(itt)
        if(itt.lt.it)sumvtl=sumvtl+vtacq(itt)
       else
        vtacq(itt)=0.d0
       endif
      enddo
      do i=1,18
       vv(i)=0.d0
      enddo       
                                  
c weights for vertex contr.: vv(1): >=2 proj. legs & >=2 targ. legs (N1 in PRD)
      vv(1)=(max(0.d0,1.d0-exp(-2.d0*vpac(ip))*(1.d0+2.d0*vpac(ip)))
     *+2.d0*vpac(ip)*exp(-2.d0*vpac(ip))*(1.d0-(1.d0-vvxp)**2))
     **(max(0.d0,1.d0-exp(-2.d0*vtac(it))*(1.d0+2.d0*vtac(it)))
     *+2.d0*vtac(it)*exp(-2.d0*vtac(it))*(1.d0-(1.d0-vvxt)**2))
     **(1.d0-vvx)**2*exp(-2.d0*sumvpl-2.d0*sumvtl)
     *-2.d0*(max(0.d0,1.d0-exp(-2.d0*vtac(it))*(1.d0+2.d0*vtac(it)))
     *+2.d0*vtac(it)*exp(-2.d0*vtac(it))*(1.d0-(1.d0-vvxt)**2))
     **(max(0.d0,exp(vpac(ip)-vpac0(ip))-1.d0-(vpac(ip)-vpac0(ip)))
     **(1.d0-vvxp0)+(vpac(ip)-vpac0(ip))*(vvxp-vvxp0))*exp(-vpac(ip))
     **(1.d0-vvx)*(1.d0-vvxtl)*exp(-sumvp-2.d0*sumvtl)
     *-2.d0*(max(0.d0,1.d0-exp(-2.d0*vpac(ip))*(1.d0+2.d0*vpac(ip)))
     *+2.d0*vpac(ip)*exp(-2.d0*vpac(ip))*(1.d0-(1.d0-vvxp)**2))
     **(max(0.d0,exp(vtac(it)-vtac0(it))-1.d0-(vtac(it)-vtac0(it)))
     **(1.d0-vvxt0)+(vtac(it)-vtac0(it))*(vvxt-vvxt0))*exp(-vtac(it))
     **(1.d0-vvx)*(1.d0-vvxpl)*exp(-2.d0*sumvpl-sumvt)
c vv(2): 0 proj. legs and 0 targ. legs (N8 in PRD)
      vv(2)=(((1.d0-exp(-vpac(ip)))**2*(1.d0-vvxpl)
     *+2.d0*(1.d0-exp(-vpac(ip)))*vvxpl)*exp(-2.d0*sumvp)*(1.d0-vvxpl)
     *+2.d0*vpacq(ip)*(1.d0-exp(-vpac(ip)-sumvpl)*(1.d0-vvxpl)
     **(1.d0-vvxp))*exp(-sumvpl))
     **(((1.d0-exp(-vtac(it)))**2*(1.d0-vvxtl)
     *+2.d0*(1.d0-exp(-vtac(it)))*vvxtl)*exp(-2.d0*sumvt)*(1.d0-vvxtl)
     *+2.d0*vtacq(it)*(1.d0-exp(-vtac(it)-sumvtl)*(1.d0-vvxtl)
     **(1.d0-vvxt))*exp(-sumvtl))
c vv(3): 0 proj. legs and >=2 targ. legs (N2 in PRD)
      vv(3)=(((1.d0-exp(-vpac(ip)))**2*(1.d0-vvxpl)
     *+2.d0*(1.d0-exp(-vpac(ip)))*vvxpl)*exp(-2.d0*sumvp)*(1.d0-vvxpl)
     *+2.d0*vpacq(ip)*(1.d0-exp(-vpac(ip)-sumvpl)*(1.d0-vvxpl)
     **(1.d0-vvxp))*exp(-sumvpl))
     **((max(0.d0,1.d0-exp(-2.d0*vtac(it))*(1.d0+2.d0*vtac(it)))
     *+2.d0*vtac(it)*exp(-2.d0*vtac(it))*(1.d0-(1.d0-vvxt)**2))
     **(1.d0-vvxtl)**2*exp(-2.d0*sumvtl)
     *-2.d0*(max(0.d0,exp(vtac(it)-vtac0(it))-1.d0-(vtac(it)
     *-vtac0(it)))*(1.d0-vvxt0)+(vtac(it)-vtac0(it))*(vvxt-vvxt0))
     **exp(-vtac(it))*(1.d0-vvxtl)*exp(-sumvt))
c vv(4): >=2 proj. legs and 0 targ. legs (N3 in PRD)
      vv(4)=((max(0.d0,1.d0-exp(-2.d0*vpac(ip))*(1.d0+2.d0*vpac(ip)))
     *+2.d0*vpac(ip)*exp(-2.d0*vpac(ip))*(1.d0-(1.d0-vvxp)**2))
     **(1.d0-vvxpl)**2*exp(-2.d0*sumvpl)
     *-2.d0*(max(0.d0,exp(vpac(ip)-vpac0(ip))-1.d0-(vpac(ip)
     *-vpac0(ip)))*(1.d0-vvxp0)+(vpac(ip)-vpac0(ip))*(vvxp-vvxp0))
     **exp(-vpac(ip))*exp(-sumvp)*(1.d0-vvxpl))
     **(((1.d0-exp(-vtac(it)))**2*(1.d0-vvxtl)
     *+2.d0*(1.d0-exp(-vtac(it)))*vvxtl)*exp(-2.d0*sumvt)*(1.d0-vvxtl)
     *+2.d0*vtacq(it)*(1.d0-exp(-vtac(it)-sumvtl)*(1.d0-vvxtl)
     **(1.d0-vvxt))*exp(-sumvtl))
c vv(5): 0 proj. legs and >=2 targ. (handle) legs (N7 in PRD)
      vv(5)=4.d0*((1.d0-exp(-vpac(ip)))*exp(-vpacq(ip))+vpacq(ip))
     **exp(-sumvpl)*(max(0.d0,exp(vtac(it)-vtac0(it))-1.d0
     *-(vtac(it)-vtac0(it)))*(1.d0-vvxt0)+(vtac(it)-vtac0(it))
     **(vvxt-vvxt0))*exp(-vtac(it)-sumvt)*(1.d0-vvx)
      if(xpomr*scm.lt.1.1d0*sgap**2)vv(5)=0.d0
c vv(6): >=2 proj. (handle) legs and 0 targ. legs (N6 in PRD)
      vv(6)=4.d0*((1.d0-exp(-vtac(it)))*exp(-vtacq(it))+vtacq(it))
     **exp(-sumvtl)*(max(0.d0,exp(vpac(ip)-vpac0(ip))-1.d0
     *-(vpac(ip)-vpac0(ip)))*(1.d0-vvxp0)+(vpac(ip)-vpac0(ip))
     **(vvxp-vvxp0))*exp(-vpac(ip)-sumvp)*(1.d0-vvx)
      if(xpomr*sgap**2.gt..9d0)vv(6)=0.d0
c vv(7): >=2 proj. legs and 1 targ. leg (N4+N9 in PRD)
      vv(7)=(max(0.d0,1.d0-exp(-2.d0*vpac(ip))*(1.d0+2.d0*vpac(ip)))
     *+2.d0*vpac(ip)*exp(-2.d0*vpac(ip))*(1.d0-(1.d0-vvxp)**2))
     **((vtac0(it)+vtlc0)*exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
     **exp(-sumvt)-(vtac(it)+vtlc-vtac0(it)-vtlc0)*(1.d0-exp(-vtac(it))
     **(1.d0-vvxt)*(1.d0-vvxtl)*exp(-sumvt)))*(1.d0-vvxt)*(1.d0-vvxpl)
     **(1.d0-vvx)*exp(-vtac(it)-sumvt-2.d0*sumvpl)
     *-2.d0*(max(0.d0,exp(vpac(ip)-vpac0(ip))-1.d0-(vpac(ip)
     *-vpac0(ip)))*(1.d0-vvxp0)+(vpac(ip)-vpac0(ip))*(vvxp-vvxp0))
     **(vtac(it)+vtlc)*(1.d0-vvxt)**2*(1.d0-vvxtl)*(1.d0-vvx)
     **exp(-vpac(ip)-2.d0*vtac(it)-sumvp-2.d0*sumvt)
c vv(8): 1 proj. leg and >=2 targ. legs (N4+N9 in PRD)
      vv(8)=(max(0.d0,1.d0-exp(-2.d0*vtac(it))*(1.d0+2.d0*vtac(it)))
     *+2.d0*vtac(it)*exp(-2.d0*vtac(it))*(1.d0-(1.d0-vvxt)**2))
     **((vpac0(ip)+vplc0)*exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
     **exp(-sumvp)-(vpac(ip)+vplc-vpac0(ip)-vplc0)*(1.d0-exp(-vpac(ip))
     **(1.d0-vvxp)*(1.d0-vvxpl)*exp(-sumvp)))*(1.d0-vvxp)*(1.d0-vvxtl)
     **(1.d0-vvx)*exp(-vpac(ip)-sumvp-2.d0*sumvtl)   
     *-2.d0*(vpac(ip)+vplc)*exp(-2.d0*vpac(ip)-vtac(it))
     **(max(0.d0,exp(vtac(it)-vtac0(it))-1.d0-(vtac(it)-vtac0(it)))
     **(1.d0-vvxt0)+(vtac(it)-vtac0(it))*(vvxt-vvxt0))*(1.d0-vvxpl)
     **(1.d0-vvxp)**2*(1.d0-vvx)*exp(-2.d0*sumvp-sumvt)
c vv(9): 0 proj. legs and 1 targ. leg (N5+N10 in PRD)
      vv(9)=(((1.d0-exp(-vpac(ip)))**2*(1.d0-vvxpl)
     *+2.d0*(1.d0-exp(-vpac(ip)))*vvxpl)*exp(-2.d0*sumvp)*(1.d0-vvxpl)
     *+2.d0*vpacq(ip)*(1.d0-exp(-vpac(ip)-sumvpl)*(1.d0-vvxpl)
     **(1.d0-vvxp))*exp(-sumvpl))
     **((vtac0(it)+vtlc0)*exp(-vtac(it)-sumvt)*(1.d0-vvxt)*(1.d0-vvxtl)
     *-(vtac(it)+vtlc-vtac0(it)-vtlc0)
     **(1.d0-exp(-vtac(it)-sumvt)*(1.d0-vvxt)*(1.d0-vvxtl)))
     **(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it)-sumvt)
c vv(10): 1 proj. leg and 0 targ. legs (N5+N10 in PRD)
      vv(10)=((vpac0(ip)+vplc0)*exp(-vpac(ip)-sumvp)*(1.d0-vvxp)
     **(1.d0-vvxpl)-(vpac(ip)+vplc-vpac0(ip)-vplc0)*(1.d0
     *-exp(-vpac(ip)-sumvp)*(1.d0-vvxp)*(1.d0-vvxpl)))*exp(-vpac(ip))
     **(((1.d0-exp(-vtac(it)))**2*(1.d0-vvxtl)
     *+2.d0*(1.d0-exp(-vtac(it)))*vvxtl)*exp(-2.d0*sumvt)*(1.d0-vvxtl)
     *+2.d0*vtacq(it)*(1.d0-exp(-vtac(it)-sumvtl)*(1.d0-vvxtl)
     **(1.d0-vvxt))*exp(-sumvtl))*(1.d0-vvxpl)*(1.d0-vvxp)*exp(-sumvp)
c vv(11): 1 cut proj. leg and 1 targ. leg (N11 in PRD minus 1P-cut)
      vv(11)=((vtlc0-vtlcpe)*exp(-vtac(it)-sumvt)*(1.d0-vvxt)
     **(1.d0-vvxtl)-(vtlc-vtlc0)*(1.d0-exp(-vtac(it)-sumvt)
     **(1.d0-vvxt)*(1.d0-vvxtl)))
     **vplcp*exp(-2.d0*vpac(ip)-vtac(it)-2.d0*sumvp-sumvt)
     **(1.d0-vvxt)*(1.d0-vvxp)**2*(1.d0-vvxpl)*(1.d0-vvx)*2.d0
      if(xpomr*scm.lt.1.1d0*sgap**2)vv(11)=0.d0
c vv(12): 1 proj. leg and 1 cut targ. leg (N11 in PRD minus 1P-cut)
      vv(12)=((vplc0-vplcpe)*exp(-vpac(ip)-sumvp)*(1.d0-vvxp)
     **(1.d0-vvxpl)-(vplc-vplc0)*(1.d0-exp(-vpac(ip)-sumvp)
     **(1.d0-vvxp)*(1.d0-vvxpl)))
     **vtlcp*exp(-2.d0*vtac(it)-vpac(ip)-sumvp-2.d0*sumvt)
     **(1.d0-vvxp)*(1.d0-vvxt)**2*(1.d0-vvxtl)*(1.d0-vvx)*2.d0
      if(xpomr*sgap**2.gt..9d0)vv(12)=0.d0
c vv(13): >=1 proj. quark-leg and >=2 targ. legs
      vv(13)=((max(0.d0,1.d0-exp(-2.d0*vtac(it))*(1.d0+2.d0*vtac(it)))
     *+2.d0*vtac(it)*exp(-2.d0*vtac(it))*(1.d0-(1.d0-vvxt)**2))
     **(1.d0-vvxtl)
     *-2.d0*(max(0.d0,exp(vtac(it)-vtac0(it))-1.d0-(vtac(it)
     *-vtac0(it)))*(1.d0-vvxt0)+(vtac(it)-vtac0(it))*(vvxt-vvxt0))
     **exp(-vtac(it)))*vpacq(ip)*exp(-2.d0*vpac(ip)-2.d0*sumvpl)
     **(1.d0-vvxp)**2*(1.d0-vvxpl)*(1.d0-vvx)*2.d0
c vv(14): >=1 targ. quark-leg and >=2 proj. legs
      vv(14)=((max(0.d0,1.d0-exp(-2.d0*vpac(ip))*(1.d0+2.d0*vpac(ip)))
     *+2.d0*vpac(ip)*exp(-2.d0*vpac(ip))*(1.d0-(1.d0-vvxp)**2))
     **(1.d0-vvxpl)
     *-2.d0*(max(0.d0,exp(vpac(ip)-vpac0(ip))-1.d0-(vpac(ip)
     *-vpac0(ip)))*(1.d0-vvxp0)+(vpac(ip)-vpac0(ip))*(vvxp-vvxp0))
     **exp(-vpac(ip)))*vtacq(it)*exp(-2.d0*vtac(it)-2.d0*sumvtl)
     **(1.d0-vvxt)**2*(1.d0-vvxtl)*(1.d0-vvx)*2.d0
c vv(15): >=1 proj. quark-leg and 0 targ. legs
      vv(15)=2.d0*vpacq(ip)*exp(-2.d0*vpac(ip)-2.d0*sumvpl)
     **(1.d0-vvxp)**2*(1.d0-vvxpl)*(1.d0-vvx)
     **((1.d0-exp(-vtac(it)))**2*(1.d0-vvxtl)
     *+2.d0*(1.d0-exp(-vtac(it)))*vvxtl)
c vv(16): >=1 proj. quark-leg and 0 targ. legs
      vv(16)=2.d0*vtacq(it)*exp(-2.d0*vtac(it)-2.d0*sumvtl)
     **(1.d0-vvxt)**2*(1.d0-vvxtl)*(1.d0-vvx)
     **((1.d0-exp(-vpac(ip)))**2*(1.d0-vvxpl)
     *+2.d0*(1.d0-exp(-vpac(ip)))*vvxpl)
c vv(17): >=1 proj. quark-leg and 1 targ. leg
      vv(17)=2.d0*vpacq(ip)*exp(-2.d0*vpac(ip)-2.d0*sumvpl)
     **(1.d0-vvxp)**2*(1.d0-vvxpl)*(1.d0-vvx)
     **((vtlc0-vtlcpe)*exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
     *-(vtlc-vtlc0)*(1.d0-exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)))
     **exp(-vtac(it))*(1.d0-vvxt)
      if(xpomr*scm.lt.1.1d0*sgap**2)vv(17)=0.d0
c vv(18): >=1  targ. quark-leg and 1 targ. leg
      vv(18)=2.d0*vtacq(it)*exp(-2.d0*vtac(it)-2.d0*sumvtl)
     **(1.d0-vvxt)**2*(1.d0-vvxtl)*(1.d0-vvx)
     **((vplc0-vplcpe)*exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
     *-(vplc-vplc0)*(1.d0-exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)))
     **exp(-vpac(ip))*(1.d0-vvxp)
      if(xpomr*sgap**2.gt..9d0)vv(18)=0.d0

      gb0=0.d0
      do i=1,18
       gb0=gb0+max(0.d0,vv(i))
      enddo       
      if(gb0.eq.0.d0)then
       if(debug.ge.3)write (moniou,202)
       iret=1
       goto 31
      endif
      if(debug.ge.3)write (moniou,203)gb0
      
1     xpomr=(scm/sgap**2)**(-qgran(b10))/sgap   !LC momentum for 3P-vertex
      rp1=(rq(iddp(ip),icz)-alfp*dlog(xpomr))*4.d0*.0389d0  
      rp2=(rq(iddt(it),2)+alfp*dlog(xpomr*scm))*4.d0*.0389d0 
      rp=rp1*rp2/(rp1+rp2)
      z=qgran(b10)
      phi=pi*qgran(b10)
      b0=dsqrt(-rp*dlog(z))
      bbpr=(bpt*rp1/(rp1+rp2)+b0*cos(phi))**2+(b0*sin(phi))**2 
      bbtg=(bpt*rp2/(rp1+rp2)-b0*cos(phi))**2+(b0*sin(phi))**2 
      call qgbdef(bbpr,bbtg,xa(ip,1)+bcoll,xa(ip,2),xb(it,1),xb(it,2)
     *,xxp,yyp,int(1.5d0+qgran(b10)))   !determine coordinates for the vertex
           
      call qgfdf(xxp,yyp,xpomr,vpac,vtac,vpht,vtht,vpacq,vtacq
     *,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)     
      vvxts=1.d0-(1.d0-vvx)*(1.d0-vvxt)*exp(-vtac(it))
      vplcps=qgfani(1.d0/xpomr,bbpr,vvxts,vvxp,vvxpl,iddp(ip),icz,16)
      vplcp=max(vplcps,vpht(ip,2)
     *+qgfani(1.d0/xpomr,bbpr,vvxts,vvxp,vvxpl,iddp(ip),icz,14))
      vplcpe=max(vplcp,vpht(ip,1)
     *+qgfani(1.d0/xpomr,bbpr,vvxts,vvxp,vvxpl,iddp(ip),icz,13))
      vplcng=max(vplcpe,vpht(ip,1)
     *+qgfani(1.d0/xpomr,bbpr,vvxts,vvxp,vvxpl,iddp(ip),icz,12))
      vplc0=max(vplcng,vpht(ip,1)
     *+qgfani(1.d0/xpomr,bbpr,vvxts,vvxp,vvxpl,iddp(ip),icz,11))
      vplc=max(vplc0,vpht(ip,1)
     *+qgfani(1.d0/xpomr,bbpr,vvxts,vvxp,vvxpl,iddp(ip),icz,10))

      vvxps=1.d0-(1.d0-vvx)*(1.d0-vvxp)*exp(-vpac(ip))
      vtlcps=qgfani(xpomr*scm,bbtg,vvxps,vvxt,vvxtl,iddt(it),2,16)
      vtlcp=max(vtlcps,vtht(it,2)
     *+qgfani(xpomr*scm,bbtg,vvxps,vvxt,vvxtl,iddt(it),2,14))
      vtlcpe=max(vtlcp,vtht(it,1)
     *+qgfani(xpomr*scm,bbtg,vvxps,vvxt,vvxtl,iddt(it),2,13))
      vtlcng=max(vtlcpe,vtht(it,1)
     *+qgfani(xpomr*scm,bbtg,vvxps,vvxt,vvxtl,iddt(it),2,12))
      vtlc0=max(vtlcng,vtht(it,1)
     *+qgfani(xpomr*scm,bbtg,vvxps,vvxt,vvxtl,iddt(it),2,11))
      vtlc=max(vtlc0,vtht(it,1)
     *+qgfani(xpomr*scm,bbtg,vvxps,vvxt,vvxtl,iddt(it),2,10))

      sumcp0=0.d0
      sumup=0.d0
      do i=1,ia(1)
       sumup=sumup+vpac(i)
      enddo
      vvxs=(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
      do i=1,ia(1)-ip+1
       ipp=ia(1)-i+1
       bbp=(xa(ipp,1)+bcoll-xxp)**2+(xa(ipp,2)-yyp)**2
       sumup=sumup-vpac(ipp)
       vps=qgfani(1.d0/xpomr,bbp,1.d0-vvxs*exp(-sumup)
     * ,1.d0-exp(-sumcp0),1.d0-exp(-sumup),iddp(ipp),icz,3)
       vpac0(ipp)=max(vps,vpht(ipp,1)+qgfani(1.d0/xpomr,bbp
     * ,1.d0-vvxs*exp(-sumup),1.d0-exp(-sumcp0),1.d0-exp(-sumup)
     * ,iddp(ipp),icz,6))
       vpac0(ipp)=min(vpac(ipp),vpac0(ipp))
       if(ipp.gt.ip)sumcp0=sumcp0+vpac0(ipp)
      enddo
      sumct0=0.d0
      sumut=0.d0
      do i=1,ia(2)
       sumut=sumut+vtac(i)
      enddo
      vvxs=(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
      do i=1,ia(2)-it+1
       itt=ia(2)-i+1
       bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
       sumut=sumut-vtac(itt)
       vts=qgfani(xpomr*scm,bbt,1.d0-vvxs*exp(-sumut)
     * ,1.d0-exp(-sumct0),1.d0-exp(-sumut),iddt(itt),2,3)
       vtac0(itt)=max(vts,vtht(itt,1)+qgfani(xpomr*scm,bbt
     * ,1.d0-vvxs*exp(-sumut),1.d0-exp(-sumct0),1.d0-exp(-sumut)
     * ,iddt(itt),2,6))
       vtac0(itt)=min(vtac(itt),vtac0(itt))
       if(itt.gt.it)sumct0=sumct0+vtac0(itt)
      enddo
      vvxp0=1.d0-exp(-sumcp0)
      vvxt0=1.d0-exp(-sumct0)
      
      sumvp=0.d0
      sumvpl=0.d0
      do ipp=1,ia(1)
       if(iqvp(ipp).eq.0)then
        sumvp=sumvp+vpacq(ipp)
        if(ipp.lt.ip)sumvpl=sumvpl+vpacq(ipp)
       else
        vpacq(ipp)=0.d0
       endif
      enddo
      sumvt=0.d0
      sumvtl=0.d0
      do itt=1,ia(2)
       if(iqvt(itt).eq.0)then
        sumvt=sumvt+vtacq(itt)
        if(itt.lt.it)sumvtl=sumvtl+vtacq(itt)
       else
        vtacq(itt)=0.d0
       endif
      enddo
      do i=1,18
       vv(i)=0.d0
      enddo       
                                  
c weights for vertex contr.: vv(1): >=2 proj. legs & >=2 targ. legs (N1 in PRD)
      vv(1)=(max(0.d0,1.d0-exp(-2.d0*vpac(ip))*(1.d0+2.d0*vpac(ip)))
     *+2.d0*vpac(ip)*exp(-2.d0*vpac(ip))*(1.d0-(1.d0-vvxp)**2))
     **(max(0.d0,1.d0-exp(-2.d0*vtac(it))*(1.d0+2.d0*vtac(it)))
     *+2.d0*vtac(it)*exp(-2.d0*vtac(it))*(1.d0-(1.d0-vvxt)**2))
     **(1.d0-vvx)**2*exp(-2.d0*sumvpl-2.d0*sumvtl)
     *-2.d0*(max(0.d0,1.d0-exp(-2.d0*vtac(it))*(1.d0+2.d0*vtac(it)))
     *+2.d0*vtac(it)*exp(-2.d0*vtac(it))*(1.d0-(1.d0-vvxt)**2))
     **(max(0.d0,exp(vpac(ip)-vpac0(ip))-1.d0-(vpac(ip)-vpac0(ip)))
     **(1.d0-vvxp0)+(vpac(ip)-vpac0(ip))*(vvxp-vvxp0))*exp(-vpac(ip))
     **(1.d0-vvx)*(1.d0-vvxtl)*exp(-sumvp-2.d0*sumvtl)
     *-2.d0*(max(0.d0,1.d0-exp(-2.d0*vpac(ip))*(1.d0+2.d0*vpac(ip)))
     *+2.d0*vpac(ip)*exp(-2.d0*vpac(ip))*(1.d0-(1.d0-vvxp)**2))
     **(max(0.d0,exp(vtac(it)-vtac0(it))-1.d0-(vtac(it)-vtac0(it)))
     **(1.d0-vvxt0)+(vtac(it)-vtac0(it))*(vvxt-vvxt0))*exp(-vtac(it))
     **(1.d0-vvx)*(1.d0-vvxpl)*exp(-2.d0*sumvpl-sumvt)
c vv(2): 0 proj. legs and 0 targ. legs (N8 in PRD)
      vv(2)=(((1.d0-exp(-vpac(ip)))**2*(1.d0-vvxpl)
     *+2.d0*(1.d0-exp(-vpac(ip)))*vvxpl)*exp(-2.d0*sumvp)*(1.d0-vvxpl)
     *+2.d0*vpacq(ip)*(1.d0-exp(-vpac(ip)-sumvpl)*(1.d0-vvxpl)
     **(1.d0-vvxp))*exp(-sumvpl))
     **(((1.d0-exp(-vtac(it)))**2*(1.d0-vvxtl)
     *+2.d0*(1.d0-exp(-vtac(it)))*vvxtl)*exp(-2.d0*sumvt)*(1.d0-vvxtl)
     *+2.d0*vtacq(it)*(1.d0-exp(-vtac(it)-sumvtl)*(1.d0-vvxtl)
     **(1.d0-vvxt))*exp(-sumvtl))
c vv(3): 0 proj. legs and >=2 targ. legs (N2 in PRD)
      vv(3)=(((1.d0-exp(-vpac(ip)))**2*(1.d0-vvxpl)
     *+2.d0*(1.d0-exp(-vpac(ip)))*vvxpl)*exp(-2.d0*sumvp)*(1.d0-vvxpl)
     *+2.d0*vpacq(ip)*(1.d0-exp(-vpac(ip)-sumvpl)*(1.d0-vvxpl)
     **(1.d0-vvxp))*exp(-sumvpl))
     **((max(0.d0,1.d0-exp(-2.d0*vtac(it))*(1.d0+2.d0*vtac(it)))
     *+2.d0*vtac(it)*exp(-2.d0*vtac(it))*(1.d0-(1.d0-vvxt)**2))
     **(1.d0-vvxtl)**2*exp(-2.d0*sumvtl)
     *-2.d0*(max(0.d0,exp(vtac(it)-vtac0(it))-1.d0-(vtac(it)
     *-vtac0(it)))*(1.d0-vvxt0)+(vtac(it)-vtac0(it))*(vvxt-vvxt0))
     **exp(-vtac(it))*(1.d0-vvxtl)*exp(-sumvt))
c vv(4): >=2 proj. legs and 0 targ. legs (N3 in PRD)
      vv(4)=((max(0.d0,1.d0-exp(-2.d0*vpac(ip))*(1.d0+2.d0*vpac(ip)))
     *+2.d0*vpac(ip)*exp(-2.d0*vpac(ip))*(1.d0-(1.d0-vvxp)**2))
     **(1.d0-vvxpl)**2*exp(-2.d0*sumvpl)
     *-2.d0*(max(0.d0,exp(vpac(ip)-vpac0(ip))-1.d0-(vpac(ip)
     *-vpac0(ip)))*(1.d0-vvxp0)+(vpac(ip)-vpac0(ip))*(vvxp-vvxp0))
     **exp(-vpac(ip))*exp(-sumvp)*(1.d0-vvxpl))
     **(((1.d0-exp(-vtac(it)))**2*(1.d0-vvxtl)
     *+2.d0*(1.d0-exp(-vtac(it)))*vvxtl)*exp(-2.d0*sumvt)*(1.d0-vvxtl)
     *+2.d0*vtacq(it)*(1.d0-exp(-vtac(it)-sumvtl)*(1.d0-vvxtl)
     **(1.d0-vvxt))*exp(-sumvtl))
c vv(5): 0 proj. legs and >=2 targ. (handle) legs (N7 in PRD)
      vv(5)=4.d0*((1.d0-exp(-vpac(ip)))*exp(-vpacq(ip))+vpacq(ip))
     **exp(-sumvpl)*(max(0.d0,exp(vtac(it)-vtac0(it))-1.d0
     *-(vtac(it)-vtac0(it)))*(1.d0-vvxt0)+(vtac(it)-vtac0(it))
     **(vvxt-vvxt0))*exp(-vtac(it)-sumvt)*(1.d0-vvx)
      if(xpomr*scm.lt.1.1d0*sgap**2)vv(5)=0.d0
c vv(6): >=2 proj. (handle) legs and 0 targ. legs (N6 in PRD)
      vv(6)=4.d0*((1.d0-exp(-vtac(it)))*exp(-vtacq(it))+vtacq(it))
     **exp(-sumvtl)*(max(0.d0,exp(vpac(ip)-vpac0(ip))-1.d0
     *-(vpac(ip)-vpac0(ip)))*(1.d0-vvxp0)+(vpac(ip)-vpac0(ip))
     **(vvxp-vvxp0))*exp(-vpac(ip)-sumvp)*(1.d0-vvx)
      if(xpomr*sgap**2.gt..9d0)vv(6)=0.d0
c vv(7): >=2 proj. legs and 1 targ. leg (N4+N9 in PRD)
      vv(7)=(max(0.d0,1.d0-exp(-2.d0*vpac(ip))*(1.d0+2.d0*vpac(ip)))
     *+2.d0*vpac(ip)*exp(-2.d0*vpac(ip))*(1.d0-(1.d0-vvxp)**2))
     **((vtac0(it)+vtlc0)*exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
     **exp(-sumvt)-(vtac(it)+vtlc-vtac0(it)-vtlc0)*(1.d0-exp(-vtac(it))
     **(1.d0-vvxt)*(1.d0-vvxtl)*exp(-sumvt)))*(1.d0-vvxt)*(1.d0-vvxpl)
     **(1.d0-vvx)*exp(-vtac(it)-sumvt-2.d0*sumvpl)
     *-2.d0*(max(0.d0,exp(vpac(ip)-vpac0(ip))-1.d0-(vpac(ip)
     *-vpac0(ip)))*(1.d0-vvxp0)+(vpac(ip)-vpac0(ip))*(vvxp-vvxp0))
     **(vtac(it)+vtlc)*(1.d0-vvxt)**2*(1.d0-vvxtl)*(1.d0-vvx)
     **exp(-vpac(ip)-2.d0*vtac(it)-sumvp-2.d0*sumvt)
c vv(8): 1 proj. leg and >=2 targ. legs (N4+N9 in PRD)
      vv(8)=(max(0.d0,1.d0-exp(-2.d0*vtac(it))*(1.d0+2.d0*vtac(it)))
     *+2.d0*vtac(it)*exp(-2.d0*vtac(it))*(1.d0-(1.d0-vvxt)**2))
     **((vpac0(ip)+vplc0)*exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
     **exp(-sumvp)-(vpac(ip)+vplc-vpac0(ip)-vplc0)*(1.d0-exp(-vpac(ip))
     **(1.d0-vvxp)*(1.d0-vvxpl)*exp(-sumvp)))*(1.d0-vvxp)*(1.d0-vvxtl)
     **(1.d0-vvx)*exp(-vpac(ip)-sumvp-2.d0*sumvtl)   
     *-2.d0*(vpac(ip)+vplc)*exp(-2.d0*vpac(ip)-vtac(it))
     **(max(0.d0,exp(vtac(it)-vtac0(it))-1.d0-(vtac(it)-vtac0(it)))
     **(1.d0-vvxt0)+(vtac(it)-vtac0(it))*(vvxt-vvxt0))*(1.d0-vvxpl)
     **(1.d0-vvxp)**2*(1.d0-vvx)*exp(-2.d0*sumvp-sumvt)
c vv(9): 0 proj. legs and 1 targ. leg (N5+N10 in PRD)
      vv(9)=(((1.d0-exp(-vpac(ip)))**2*(1.d0-vvxpl)
     *+2.d0*(1.d0-exp(-vpac(ip)))*vvxpl)*exp(-2.d0*sumvp)*(1.d0-vvxpl)
     *+2.d0*vpacq(ip)*(1.d0-exp(-vpac(ip)-sumvpl)*(1.d0-vvxpl)
     **(1.d0-vvxp))*exp(-sumvpl))
     **((vtac0(it)+vtlc0)*exp(-vtac(it)-sumvt)*(1.d0-vvxt)*(1.d0-vvxtl)
     *-(vtac(it)+vtlc-vtac0(it)-vtlc0)
     **(1.d0-exp(-vtac(it)-sumvt)*(1.d0-vvxt)*(1.d0-vvxtl)))
     **(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it)-sumvt)
c vv(10): 1 proj. leg and 0 targ. legs (N5+N10 in PRD)
      vv(10)=((vpac0(ip)+vplc0)*exp(-vpac(ip)-sumvp)*(1.d0-vvxp)
     **(1.d0-vvxpl)-(vpac(ip)+vplc-vpac0(ip)-vplc0)*(1.d0
     *-exp(-vpac(ip)-sumvp)*(1.d0-vvxp)*(1.d0-vvxpl)))
     **exp(-vpac(ip)-sumvp)*(((1.d0-exp(-vtac(it)))**2*(1.d0-vvxtl)
     *+2.d0*(1.d0-exp(-vtac(it)))*vvxtl)*exp(-2.d0*sumvt)*(1.d0-vvxtl)
     *+2.d0*vtacq(it)*(1.d0-exp(-vtac(it)-sumvtl)*(1.d0-vvxtl)
     **(1.d0-vvxt))*exp(-sumvtl))*(1.d0-vvxpl)*(1.d0-vvxp)
c vv(11): 1 cut proj. leg and 1 targ. leg (N11 in PRD minus 1P-cut)
      vv(11)=((vtlc0-vtlcpe)*exp(-vtac(it)-sumvt)*(1.d0-vvxt)
     **(1.d0-vvxtl)-(vtlc-vtlc0)*(1.d0-exp(-vtac(it)-sumvt)
     **(1.d0-vvxt)*(1.d0-vvxtl)))
     **vplcp*exp(-2.d0*vpac(ip)-vtac(it)-2.d0*sumvp-sumvt)
     **(1.d0-vvxt)*(1.d0-vvxp)**2*(1.d0-vvxpl)*(1.d0-vvx)*2.d0
      if(xpomr*scm.lt.1.1d0*sgap**2)vv(11)=0.d0
c vv(12): 1 proj. leg and 1 cut targ. leg (N11 in PRD minus 1P-cut)
      vv(12)=((vplc0-vplcpe)*exp(-vpac(ip)-sumvp)*(1.d0-vvxp)
     **(1.d0-vvxpl)-(vplc-vplc0)*(1.d0-exp(-vpac(ip)-sumvp)
     **(1.d0-vvxp)*(1.d0-vvxpl)))
     **vtlcp*exp(-2.d0*vtac(it)-vpac(ip)-sumvp-2.d0*sumvt)
     **(1.d0-vvxp)*(1.d0-vvxt)**2*(1.d0-vvxtl)*(1.d0-vvx)*2.d0
      if(xpomr*sgap**2.gt..9d0)vv(12)=0.d0
c vv(13): >=1 proj. quark-leg and >=2 targ. legs
      vv(13)=((max(0.d0,1.d0-exp(-2.d0*vtac(it))*(1.d0+2.d0*vtac(it)))
     *+2.d0*vtac(it)*exp(-2.d0*vtac(it))*(1.d0-(1.d0-vvxt)**2))
     **(1.d0-vvxtl)
     *-2.d0*(max(0.d0,exp(vtac(it)-vtac0(it))-1.d0-(vtac(it)
     *-vtac0(it)))*(1.d0-vvxt0)+(vtac(it)-vtac0(it))*(vvxt-vvxt0))
     **exp(-vtac(it)))*vpacq(ip)*exp(-2.d0*vpac(ip)-2.d0*sumvpl)
     **(1.d0-vvxp)**2*(1.d0-vvxpl)*(1.d0-vvx)*2.d0
c vv(14): >=1 targ. quark-leg and >=2 proj. legs
      vv(14)=((max(0.d0,1.d0-exp(-2.d0*vpac(ip))*(1.d0+2.d0*vpac(ip)))
     *+2.d0*vpac(ip)*exp(-2.d0*vpac(ip))*(1.d0-(1.d0-vvxp)**2))
     **(1.d0-vvxpl)
     *-2.d0*(max(0.d0,exp(vpac(ip)-vpac0(ip))-1.d0-(vpac(ip)
     *-vpac0(ip)))*(1.d0-vvxp0)+(vpac(ip)-vpac0(ip))*(vvxp-vvxp0))
     **exp(-vpac(ip)))*vtacq(it)*exp(-2.d0*vtac(it)-2.d0*sumvtl)
     **(1.d0-vvxt)**2*(1.d0-vvxtl)*(1.d0-vvx)*2.d0
c vv(15): >=1 proj. quark-leg and 0 targ. legs
      vv(15)=2.d0*vpacq(ip)*exp(-2.d0*vpac(ip)-2.d0*sumvpl)
     **(1.d0-vvxp)**2*(1.d0-vvxpl)*(1.d0-vvx)
     **((1.d0-exp(-vtac(it)))**2*(1.d0-vvxtl)
     *+2.d0*(1.d0-exp(-vtac(it)))*vvxtl)
c vv(16): >=1 proj. quark-leg and 0 targ. legs
      vv(16)=2.d0*vtacq(it)*exp(-2.d0*vtac(it)-2.d0*sumvtl)
     **(1.d0-vvxt)**2*(1.d0-vvxtl)*(1.d0-vvx)
     **((1.d0-exp(-vpac(ip)))**2*(1.d0-vvxpl)
     *+2.d0*(1.d0-exp(-vpac(ip)))*vvxpl)
c vv(17): >=1 proj. quark-leg and 1 targ. leg
      vv(17)=2.d0*vpacq(ip)*exp(-2.d0*vpac(ip)-2.d0*sumvpl)
     **(1.d0-vvxp)**2*(1.d0-vvxpl)*(1.d0-vvx)
     **((vtlc0-vtlcpe)*exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
     *-(vtlc-vtlc0)*(1.d0-exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)))
     **exp(-vtac(it))*(1.d0-vvxt)
      if(xpomr*scm.lt.1.1d0*sgap**2)vv(17)=0.d0
c vv(18): >=1  targ. quark-leg and 1 targ. leg
      vv(18)=2.d0*vtacq(it)*exp(-2.d0*vtac(it)-2.d0*sumvtl)
     **(1.d0-vvxt)**2*(1.d0-vvxtl)*(1.d0-vvx)
     **((vplc0-vplcpe)*exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
     *-(vplc-vplc0)*(1.d0-exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)))
     **exp(-vpac(ip))*(1.d0-vvxp)
      if(xpomr*sgap**2.gt..9d0)vv(18)=0.d0

      gb=0.d0
      do i=1,18
       vv(i)=max(0.d0,vv(i))
       gb=gb+vv(i)
      enddo       
      gb=gb/gb0/z*rp/rp0  /max(8.1d0,1.3d0*dlog10(scm))
      if(debug.ge.5)write (moniou,204)xpomr,bbpr,bbtg,gb
      if(qgran(b10).gt.gb)goto 1
      if(debug.ge.3)write (moniou,205)xpomr,bbpr,bbtg,xxp,yyp

      sumcpt=0.d0
      sumup=0.d0
      do i=1,ia(1)
       sumup=sumup+vpac(i)
      enddo
      vvxs=(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
      do i=1,ia(1)-ip+1
       ipp=ia(1)-i+1
       bbp=(xa(ipp,1)+bcoll-xxp)**2+(xa(ipp,2)-yyp)**2
       sumup=sumup-vpac(ipp)
       vpact(ipp)=max(vpac(ipp),vpht(ipp,1)
     * +qgfani(1.d0/xpomr,bbp,1.d0-vvxs*exp(-sumup)
     * ,1.d0-exp(-sumcpt),1.d0-exp(-sumup),iddp(ipp),icz,9))
       if(ipp.gt.ip)sumcpt=sumcpt+vpact(ipp)
      enddo
      sumctt=0.d0
      sumut=0.d0
      do i=1,ia(2)
       sumut=sumut+vtac(i)
      enddo
      vvxs=(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
      do i=1,ia(2)-it+1
       itt=ia(2)-i+1
       bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
       sumut=sumut-vtac(itt)
       vtact(itt)=max(vtac(itt),vtht(itt,1)
     * +qgfani(xpomr*scm,bbt,1.d0-vvxs*exp(-sumut)
     * ,1.d0-exp(-sumctt),1.d0-exp(-sumut),iddt(itt),2,9))
       if(itt.gt.it)sumctt=sumctt+vtact(itt)
      enddo
      vvxpt=1.d0-exp(-sumcpt)
      vvxtt=1.d0-exp(-sumctt)

      vvt=0.d0
      do i=1,18
       vvt=vvt+max(0.d0,vv(i))
      enddo
      if(.not.(vvt.gt.0.d0))stop'vvt<0'  !?????????????

      aks=qgran(b10)*vvt
      do jt=1,18
       aks=aks-vv(jt)
       if(aks.lt.0.d0)goto 2
      enddo
      stop'jt>18!'
      
2     if(xpomr*scm.gt.sgap**2)then
       wzgp=max(0.d0,min(0.d0,1.d0-exp(-vtac0(it))-vtac0(it)) 
     * -min(0.d0,1.d0-exp(-vtact(it))-vtact(it)))
     * -max(0.d0,vtact(it)-vtac0(it))*(1.d0-(1.d0-vvxt)*(1.d0-vvxtl)
     * *exp(-vtac(it)))
     * +vvxtt*(1.d0-exp(-vtact(it)))-vvxt0*(1.d0-exp(-vtac0(it)))
       wzgp=2.d0*max(0.d0,wzgp)
       wzgp1=2.d0*vpac(ip)*exp(-2.d0*vpac(ip))*(1.d0-vvxp)**2
     * *(1.d0-vvxpl)**2*wzgp                                !1 cut fan
       wzgp=(max(0.d0,1.d0-exp(-2.d0*vpac(ip))*(1.d0+2.d0*vpac(ip)))
     * +2.d0*vpac(ip)*exp(-2.d0*vpac(ip))*(1.d0-(1.d0-vvxp)**2))
     * *(1.d0-vvxpl)**2*wzgp                                !>1 cut fans
      else
       wzgp=0.d0
       wzgp1=0.d0
      endif
      if(xpomr*sgap**2.lt.1.d0)then
       wzgt=max(0.d0,min(0.d0,1.d0-exp(-vpac0(ip))-vpac0(ip)) 
     * -min(0.d0,1.d0-exp(-vpact(ip))-vpact(ip)))
     * -max(0.d0,vpact(ip)-vpac0(ip))*(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)
     * *exp(-vpac(ip)))
     * +vvxpt*(1.d0-exp(-vpact(ip)))-vvxp0*(1.d0-exp(-vpac0(ip)))
       wzgt=2.d0*max(0.d0,wzgt)
       wzgt1=2.d0*vtac(it)*exp(-2.d0*vtac(it))*(1.d0-vvxt)**2
     * *(1.d0-vvxtl)**2*wzgt                                !1 cut fan
       wzgt=(max(0.d0,1.d0-exp(-2.d0*vtac(it))*(1.d0+2.d0*vtac(it)))
     * +2.d0*vtac(it)*exp(-2.d0*vtac(it))*(1.d0-(1.d0-vvxt)**2))
     * *(1.d0-vvxtl)**2*wzgt                                !1 cut fan
      else
       wzgt=0.d0
       wzgt1=0.d0
      endif

      nppr0=0
      nptg0=0
      npprh0=0
      nptgh0=0
      wgpr0=0.d0
      wgtg0=0.d0
      if(jt.eq.1.or.jt.eq.4.or.jt.eq.7.or.jt.eq.14)then         !>=2 proj. fans
       ntry=0
3      ntry=ntry+1
       npprh0=0
       if(ip.eq.ia(1).or.ntry.gt.100)then
        nppr0=npgen(2.d0*vpac(ip),2,npmax)                !number of proj. fans
        do i=1,nppr0
         if(qgran(b10).le.vpac0(ip)/vpac(ip).or.xpomr*sgap**2.gt..9d0)
     *   then
          itypr0(i)=0                 !cut handle fan
         else
          npprh0=npprh0+1
          itypr0(i)=1                 !uncut handle fan
         endif
         ippr0(i)=ip
        enddo
        if(npprh0.ne.0)wh=(vpac(ip)/vpac0(ip)-1.d0)/npprh0
       else
        nppr0=npgen(2.d0*vpac(ip),1,npmax)
        do i=1,nppr0
         if(qgran(b10).le.vpac0(ip)/vpac(ip).or.xpomr*sgap**2.gt..9d0)
     *   then
          itypr0(i)=0
         else
          npprh0=npprh0+1
          itypr0(i)=1
         endif
         ippr0(i)=ip
        enddo
        if(npprh0.ne.0)wh=(vpac(ip)/vpac0(ip)-1.d0)/npprh0
        do ipp=ip+1,ia(1)
         ninc=npgen(2.d0*vpac(ipp),0,npmax)
         if(ninc.ne.0)then
          nppr0=nppr0+ninc
          nh0=npprh0
          if(nppr0.gt.legmax)then
           iret=1
           goto 31
          endif
          do i=nppr0-ninc+1,nppr0
           if(qgran(b10).le.vpac0(ipp)/vpac(ipp)
     *     .or.xpomr*sgap**2.gt..9d0)then
            itypr0(i)=0
           else
            npprh0=npprh0+1
            itypr0(i)=1
           endif
           ippr0(i)=ipp
          enddo
          if(npprh0.gt.nh0)wh=(vpac(ipp)/vpac0(ipp)-1.d0)/(npprh0-nh0)
         endif
        enddo
        if(nppr0.eq.1)goto 3
       endif
       if(nppr0.le.npprh0+1)then
        if((1.d0-vvxp)*(1.d0-vvxpl).le.0.d0)stop'qg3pdf: vvxp=0!'
        if(jt.ne.7)then
         wh0=1.d0-exp(vpac(ip)+(1.d0-nppr0)*dlog(2.d0))
     *   /(1.d0-vvxp)/(1.d0-vvxpl)
        else
         wh0=1.d0-exp(vpac(ip)+(1.d0-nppr0)*dlog(2.d0))
     *   /(1.d0-vvxp)/(1.d0-vvxpl)
     *   *(vtac(it)+vtlc)*exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
     *   /((vtac0(it)+vtlc0)*exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
     *   -(vtac(it)+vtlc-vtac0(it)-vtlc0)
     *   *(1.d0-exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)))
        endif
        if(wh0.lt.0.d0.and.(nppr0.eq.npprh0.or.nppr0.eq.npprh0+1
     *  .and.qgran(b10).gt.1.d0+wh*wh0))goto 3
       endif
      endif
        
      if(jt.eq.1.or.jt.eq.3.or.jt.eq.8.or.jt.eq.13)then
       ntry=0
4      ntry=ntry+1
       nptgh0=0
       if(it.eq.ia(2).or.ntry.gt.100)then
        nptg0=npgen(2.d0*vtac(it),2,npmax)
        do i=1,nptg0
         if(qgran(b10).le.vtac0(it)/vtac(it)
     *   .or.xpomr*scm.lt.1.1d0*sgap**2)then
          itytg0(i)=0
         else
          nptgh0=nptgh0+1
          itytg0(i)=1
         endif
         iptg0(i)=it
        enddo
        if(nptgh0.ne.0)wh=(vtac(it)/vtac0(it)-1.d0)/nptgh0
       else
        nptg0=npgen(2.d0*vtac(it),1,npmax)
        do i=1,nptg0
         if(qgran(b10).le.vtac0(it)/vtac(it)
     *   .or.xpomr*scm.lt.1.1d0*sgap**2)then
          itytg0(i)=0
         else
          nptgh0=nptgh0+1
          itytg0(i)=1
         endif
         iptg0(i)=it
        enddo
        if(nptgh0.ne.0)wh=(vtac(it)/vtac0(it)-1.d0)/nptgh0
        do itt=it+1,ia(2)
         ninc=npgen(2.d0*vtac(itt),0,npmax)
         if(ninc.ne.0)then
          nptg0=nptg0+ninc
          nh0=nptgh0
          if(nptg0.gt.legmax)then
           iret=1
           goto 31
          endif
          do i=nptg0-ninc+1,nptg0
           if(qgran(b10).le.vtac0(itt)/vtac(itt)
     *     .or.xpomr*scm.lt.1.1d0*sgap**2)then
            itytg0(i)=0
           else
            nptgh0=nptgh0+1
            itytg0(i)=1
           endif
           iptg0(i)=itt
          enddo
          if(nptgh0.gt.nh0)wh=(vtac(itt)/vtac0(itt)-1.d0)/(nptgh0-nh0)
         endif
        enddo
        if(nptg0.eq.1)goto 4
       endif
       if(nptg0.le.nptgh0+1)then
        if((1.d0-vvxt)*(1.d0-vvxtl).le.0.d0)stop'qg3pdf: vvxt=0!'
        if(jt.ne.8)then
         wh0=1.d0-exp(vtac(it)+(1.d0-nptg0)*dlog(2.d0))
     *   /(1.d0-vvxt)/(1.d0-vvxtl)
        else
         wh0=1.d0-exp(vtac(it)+(1.d0-nptg0)*dlog(2.d0))
     *   /(1.d0-vvxt)/(1.d0-vvxtl)
     *   *(vpac(ip)+vplc)*exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
     *   /((vpac0(ip)+vplc0)*exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
     *   -(vpac(ip)+vplc-vpac0(ip)-vplc0)
     *   *(1.d0-exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)))
        endif
        if(wh0.lt.0.d0.and.(nptg0.eq.nptgh0.or.nptg0.eq.nptgh0+1
     *  .and.qgran(b10).gt.1.d0+wh*wh0))goto 4
       endif
      endif
      
      if(jt.eq.6)then
       ntry=0
5      ntry=ntry+1
       if(ip.eq.ia(1).or.ntry.gt.100)then
        nppr0=npgen(vpac(ip)-vpac0(ip),2,npmax)
        do i=1,nppr0
         itypr0(i)=1
         ippr0(i)=ip
        enddo
       else
        nppr0=npgen(vpac(ip)-vpac0(ip),1,npmax)
        do i=1,nppr0
         itypr0(i)=1
         ippr0(i)=ip
        enddo
        do ipp=ip+1,ia(1)
         ninc=npgen(vpac(ipp)-vpac0(ipp),0,npmax)
         if(ninc.ne.0)then
          nppr0=nppr0+ninc
          if(nppr0.gt.legmax)then
           iret=1
           goto 31
          endif
          do i=nppr0-ninc+1,nppr0
           itypr0(i)=1
           ippr0(i)=ipp
          enddo
         endif
        enddo
        if(nppr0.eq.1)goto 5
       endif
      endif
        
      if(jt.eq.5)then
       ntry=0
6      ntry=ntry+1
       if(it.eq.ia(2).or.ntry.gt.100)then
        nptg0=npgen(vtac(it)-vtac0(it),2,npmax)
        do i=1,nptg0
         itytg0(i)=1
         iptg0(i)=it
        enddo
       else
        nptg0=npgen(vtac(it)-vtac0(it),1,npmax)
        do i=1,nptg0
         itytg0(i)=1
         iptg0(i)=it
        enddo
        do itt=it+1,ia(2)
         ninc=npgen(vtac(itt)-vtac0(itt),0,npmax)
         if(ninc.ne.0)then
          nptg0=nptg0+ninc
          if(nptg0.gt.legmax)then
           iret=1
           goto 31
          endif
          do i=nptg0-ninc+1,nptg0
           itytg0(i)=1
           iptg0(i)=itt
          enddo
         endif
        enddo
        if(nptg0.eq.1)goto 6
       endif
      endif
      
      gbt=1.d0
      if((jt.eq.1.and.nptgh0.lt.nptg0.or.jt.eq.4.or.jt.eq.14)
     *.and.npprh0.eq.nppr0)then
       gbt=1.d0-exp(vpac(ip)+(1.d0-nppr0)*dlog(2.d0))
     * /(1.d0-vvxp)/(1.d0-vvxpl)
      elseif((jt.eq.1.and.npprh0.lt.nppr0.or.jt.eq.3.or.jt.eq.13)
     *.and.nptgh0.eq.nptg0)then
       gbt=1.d0-exp(vtac(it)+(1.d0-nptg0)*dlog(2.d0))
     * /(1.d0-vvxt)/(1.d0-vvxtl)
      elseif(jt.eq.1.and.nptgh0.eq.nptg0.and.npprh0.eq.nppr0)then
       gbt=1.d0-exp(vpac(ip)+(1.d0-nppr0)*dlog(2.d0))
     * /(1.d0-vvxp)/(1.d0-vvxpl)
     * -exp(vtac(it)+(1.d0-nptg0)*dlog(2.d0))/(1.d0-vvxt)/(1.d0-vvxtl)
      elseif(jt.eq.7.and.npprh0.eq.nppr0)then
       gbt=1.d0-exp(vpac(ip)+(1.d0-nppr0)*dlog(2.d0))
     * /(1.d0-vvxp)/(1.d0-vvxpl)
     * *(vtac(it)+vtlc)*exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
     * /((vtac0(it)+vtlc0)*exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
     * -(vtac(it)+vtlc-vtac0(it)-vtlc0)
     * *(1.d0-exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)))
      elseif(jt.eq.8.and.nptgh0.eq.nptg0)then
       gbt=1.d0-exp(vtac(it)+(1.d0-nptg0)*dlog(2.d0))
     * /(1.d0-vvxt)/(1.d0-vvxtl)
     * *(vpac(ip)+vplc)*exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
     * /((vpac0(ip)+vplc0)*exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
     * -(vpac(ip)+vplc-vpac0(ip)-vplc0)
     * *(1.d0-exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)))
      endif
      if(qgran(b10).gt.gbt)goto 2
       
      if(jt.eq.7.or.jt.eq.9.or.jt.eq.11.or.jt.eq.12.or.jt.eq.17)then
       nptg0=1
       iptg0(1)=it
      endif
      if(jt.eq.8.or.jt.eq.10.or.jt.eq.11.or.jt.eq.12.or.jt.eq.18)then
       nppr0=1
       ippr0(1)=ip
      endif

c less important part of 'zigzag' cuts (sub-per cent effect)
      if(jt.eq.1.and.qgran(b10).lt.wzgp/vv(1))then
       nppr0=0
       npprh0=0
      endif
      if(jt.eq.8.and.qgran(b10).lt.wzgp1/vv(8))then
       nppr0=0
       npprh0=0
      endif
      if(jt.eq.1.and.qgran(b10).lt.wzgt/vv(1))then
       nptg0=0
       nptgh0=0
      endif
      if(jt.eq.7.and.qgran(b10).lt.wzgt1/vv(7))then
       nptg0=0
       nptgh0=0
      endif

      if(jt.eq.8.and.nptgh0.lt.nptg0.or.jt.eq.10)then !'fan' from cut vertex
       vpacng=min(vpac0(ip),vpht(ip,1)
     * +qgfani(1.d0/xpomr,bbpr,vvxts,vvxp0,vvxpl,iddp(ip),icz,7))
       vpacng=max(vplcng,vpacng)
       factor=exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
       wng=(vpacng+vplcng)*factor
       wgap=max(0.d0,(vpac0(ip)+vplc0)*factor
     * -(vpac(ip)+vplc-vpac0(ip)-vplc0)*(1.d0-factor)-wng)
       if(qgran(b10).ge.wgap/(wgap+wng).or.xpomr*sgap**2.gt..9d0)then
        if(qgran(b10).lt.vpacng/(vpacng+vplcng))then
         itypr0(1)=2            !cut 'fan' (no gap at the end)
        else
         itypr0(1)=4            !cut 'leg' (no gap at the end)
        endif
       else
        wfg=max(0.d0,(vpac0(ip)-vpacng)*factor
     *  -(vpac(ip)-vpac0(ip))*(1.d0-factor))
        wlg=max(0.d0,(vplc0-vplcng)*factor-(vplc-vplc0)*(1.d0-factor))
        if(qgran(b10).lt.wfg/(wfg+wlg))then
         itypr0(1)=3            !cut 'fan' (gap at the end)
         wgpr0=1.d0-wfg/factor/(vpac0(ip)-vpacng)
        else
         itypr0(1)=5            !cut 'leg' (gap at the end)
         wgpr0=1.d0-wlg/factor/(vplc0-vplcng)
        endif
       endif
       
      elseif(jt.eq.8.and.nptgh0.eq.nptg0)then !'fan' from cut/uncut vertex
       vpacng=min(vpac0(ip),vpht(ip,1)
     * +qgfani(1.d0/xpomr,bbpr,vvxts,vvxp0,vvxpl,iddp(ip),icz,7))
       vpacng=max(vplcng,vpacng)
       factor=exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
       wng=(vpacng+vplcng)*factor*(1.d0-exp(vtac(it)
     * +(1.d0-nptg0)*dlog(2.d0))/(1.d0-vvxt)/(1.d0-vvxtl))
       wgap=max(0.d0,(vpac0(ip)+vplc0)*factor
     * -(vpac(ip)+vplc-vpac0(ip)-vplc0)*(1.d0-factor)
     * -exp(vtac(it)+(1.d0-nptg0)*dlog(2.d0))/(1.d0-vvxt)/(1.d0-vvxtl)
     * *(vpac(ip)+vplc)*factor-wng)
       if(qgran(b10).ge.wgap/(wgap+wng).or.xpomr*sgap**2.gt..9d0)then
        if(qgran(b10).lt.vpacng/(vpacng+vplcng))then
         itypr0(1)=2            !cut 'fan' (no gap at the end)
        else
         itypr0(1)=4            !cut 'leg' (no gap at the end)
        endif
       else
        wfg=max(0.d0,(vpac0(ip)-vpacng)*factor
     *  -(vpac(ip)-vpac0(ip))*(1.d0-factor)
     *  -exp(vtac(it)+(1.d0-nptg0)*dlog(2.d0))
     *  /(1.d0-vvxt)/(1.d0-vvxtl)*(vpac(ip)-vpacng)*factor)
        wlg=max(0.d0,(vplc0-vplcng)*factor-(vplc-vplc0)*(1.d0-factor)
     *  -exp(vtac(it)+(1.d0-nptg0)*dlog(2.d0))
     *  /(1.d0-vvxt)/(1.d0-vvxtl)*(vplc-vplcng)*factor)
        if(qgran(b10).lt.wfg/(wfg+wlg))then
         itypr0(1)=3            !cut 'fan' (gap at the end)
         wgpr0=(vpac(ip)-vpac0(ip))*(1.d0-factor)
     *   /(wfg+(vpac(ip)-vpac0(ip))*(1.d0-factor))
        else
         itypr0(1)=5            !cut 'leg' (gap at the end)
         wgpr0=(vplc-vplc0)*(1.d0-factor)
     *   /(wlg+(vplc-vplc0)*(1.d0-factor))
        endif
       endif
    
      elseif(jt.eq.11)then
       itypr0(1)=6
      elseif(jt.eq.12.or.jt.eq.18)then
       factor=exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
       wng=max(0.d0,vplcng-vplcpe)*factor
     * /((vplc0-vplcpe)*factor-(vplc-vplc0)*(1.d0-factor))
       if(qgran(b10).le.wng)then
        itypr0(1)=7            !cut 'leg' (>1 cut Poms at the end)
       else
        itypr0(1)=5            !cut 'leg' (gap at the end)
        wgpr0=(1.d0-factor)/factor*(vplc-vplc0)/(vplc0-vplcng)
       endif
      endif
    
      if(jt.eq.7.and.npprh0.lt.nppr0.or.jt.eq.9)then !'fan' from cut vertex
       vtacng=min(vtac0(it),vtht(it,1)
     * +qgfani(xpomr*scm,bbtg,vvxps,vvxt0,vvxtl,iddt(it),2,7))
       vtacng=max(vtlcng,vtacng)
       factor=exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
       wng=(vtacng+vtlcng)*factor
       wgap=max(0.d0,(vtac0(it)+vtlc0)*factor
     * -(vtac(it)+vtlc-vtac0(it)-vtlc0)*(1.d0-factor)-wng)
       if(qgran(b10).ge.wgap/(wgap+wng)
     * .or.xpomr*scm.lt.1.1d0*sgap**2)then
        if(qgran(b10).lt.vtacng/(vtacng+vtlcng))then
         itytg0(1)=2            !cut 'fan' (no gap at the end)
        else
         itytg0(1)=4            !cut 'leg' (no gap at the end)
        endif
       else
        wfg=max(0.d0,(vtac0(it)-vtacng)*factor
     *  -(vtac(it)-vtac0(it))*(1.d0-factor))
        wlg=max(0.d0,(vtlc0-vtlcng)*factor-(vtlc-vtlc0)*(1.d0-factor))
        if(qgran(b10).lt.wfg/(wfg+wlg))then
         itytg0(1)=3            !cut 'fan' (gap at the end)
         wgtg0=(1.d0-factor)/factor
     *   *(vtac(it)-vtac0(it))/(vtac0(it)-vtacng)
        else
         itytg0(1)=5            !cut 'leg' (gap at the end)
         wgtg0=(1.d0-factor)/factor*(vtlc-vtlc0)/(vtlc0-vtlcng)
        endif
       endif
       
      elseif(jt.eq.7.and.npprh0.eq.nppr0)then !'fan' from cut/uncut vertex
       vtacng=min(vtac0(it),vtht(it,1)
     * +qgfani(xpomr*scm,bbtg,vvxps,vvxt0,vvxtl,iddt(it),2,7))
       vtacng=max(vtlcng,vtacng)
       factor=exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
       wng=(vtacng+vtlcng)*factor*(1.d0-exp(vpac(ip)
     * +(1.d0-nppr0)*dlog(2.d0))/(1.d0-vvxp)/(1.d0-vvxpl))
       wgap=max(0.d0,(vtac0(it)+vtlc0)*factor
     * -(vtac(it)+vtlc-vtac0(it)-vtlc0)*(1.d0-factor)
     * -exp(vpac(ip)+(1.d0-nppr0)*dlog(2.d0))/(1.d0-vvxp)/(1.d0-vvxpl)
     * *(vtac(it)+vtlc)*factor-wng)
       if(qgran(b10).ge.wgap/(wgap+wng)
     * .or.xpomr*scm.lt.1.1d0*sgap**2)then
        if(qgran(b10).lt.vtacng/(vtacng+vtlcng))then
         itytg0(1)=2            !cut 'fan' (no gap at the end)
        else
         itytg0(1)=4            !cut 'leg' (no gap at the end)
        endif
       else
        wfg=max(0.d0,(vtac0(it)-vtacng)*factor
     *  -(vtac(it)-vtac0(it))*(1.d0-factor)
     *  -exp(vpac(ip)+(1.d0-nppr0)*dlog(2.d0))
     *  /(1.d0-vvxp)/(1.d0-vvxpl)*(vtac(it)-vtacng)*factor)
        wlg=max(0.d0,(vtlc0-vtlcng)*factor-(vtlc-vtlc0)*(1.d0-factor)
     *  -exp(vpac(ip)+(1.d0-nppr0)*dlog(2.d0))
     *  /(1.d0-vvxp)/(1.d0-vvxpl)*(vtlc-vtlcng)*factor)
        if(qgran(b10).lt.wfg/(wfg+wlg))then
         itytg0(1)=3            !cut 'fan' (gap at the end)
         wgtg0=(1.d0-factor)*(vtac(it)-vtac0(it))
     *   /(wfg+(1.d0-factor)*(vtac(it)-vtac0(it)))
        else
         itytg0(1)=5            !cut 'leg' (gap at the end)
         wgtg0=(1.d0-factor)*(vtlc-vtlc0)
     *   /(wlg+(1.d0-factor)*(vtlc-vtlc0))
        endif
       endif
    
      elseif(jt.eq.12)then
       itytg0(1)=6
      elseif(jt.eq.11.or.jt.eq.17)then
       factor=exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
       wng=max(0.d0,vtlcng-vtlcpe)*factor
     * /((vtlc0-vtlcpe)*factor-(vtlc-vtlc0)*(1.d0-factor))
       if(qgran(b10).le.wng)then
        itytg0(1)=7            !cut 'leg' (>1 cut Poms at the end)
       else
        itytg0(1)=5            !cut 'leg' (gap at the end)
        wgtg0=(1.d0-factor)/factor*(vtlc-vtlc0)/(vtlc0-vtlcng)
       endif
      endif
      if(debug.ge.3)write (moniou,206)nppr0,nptg0
      
      nppr=0
      nptg=0      
      npin=0
      nzpi=0
      
      if(nppr0.eq.1.and.itypr0(1).eq.6)then     !single cut Pomeron
       nppr=1
       if(qgran(b10).le.vplcps/vplcp.or.xpomr*s2min.gt..9d0)then
        idpompi(nppr)=0
       else
        idpompi(nppr)=1
       endif
       xpompi(nppr)=xpomr
       ipompi(nppr)=ip
       bpompi(nppr,1)=bbpr
       bpompi(nppr,2)=bbtg
       if(debug.ge.4)write (moniou,209)nppr,ip,bbpr,xpompi(nppr)
       nppr0=0
      endif
      if(nptg0.eq.1.and.itytg0(1).eq.6)then     !single cut Pomeron
       nptg=1
       if(qgran(b10).le.vtlcps/vtlcp.or.xpomr*scm.lt.1.1d0*s2min)then
        idpomti(nptg)=0
       else
        idpomti(nptg)=1
       endif
       xpomti(nptg)=xpomr
       ipomti(nptg)=it
       bpomti(nptg,1)=bbtg
       bpomti(nptg,2)=bbpr
       if(debug.ge.4)write (moniou,217)nptg,it,bbtg,xpomti(nptg)
       nptg0=0
      endif
      
      vvxps=vvxp
      vvxpls=vvxpl
      vvxp0s=vvxp0
      if(nppr0.ne.0)then
       i=0
7      i=i+1
       ityp=itypr0(i)
       if(ityp.eq.0.or.ityp.eq.2.or.ityp.eq.4)then
        ipp=ippr0(i)
        bbp=(xa(ipp,1)+bcoll-xxp)**2+(xa(ipp,2)-yyp)**2
        vvxp=0.d0
        vvxpl=0.d0
        vvxp0=0.d0
        if(ia(1).gt.1)then
         do l=1,ia(1)
          if(l.lt.ipp)then
           vvxpl=vvxpl+vpac(l)
          elseif(l.gt.ipp)then
           vvxp=vvxp+vpac(l)
           vvxp0=vvxp0+vpac0(l)
          endif
         enddo
        endif
        vvxp=1.d0-exp(-vvxp)
        vvxpl=1.d0-exp(-vvxpl)
        vvxp0=1.d0-exp(-vvxp0)
        vvxts=1.d0-(1.d0-vvxt)*(1.d0-vvxtl)*(1.d0-vvxpl)*exp(-vtac(it))
        vplcps=qgfani(1.d0/xpomr,bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,16)
        vplcp=max(vplcps,vpht(ipp,2)
     *  +qgfani(1.d0/xpomr,bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,14))
        if(ityp.ne.4)then
         vpacpe=max(vplcp,vpht(ipp,1)
     *   +qgfani(1.d0/xpomr,bbp,vvxts,vvxp0,vvxpl,iddp(ipp),icz,8))
         vpacng=max(vpacpe,vpht(ipp,1)
     *   +qgfani(1.d0/xpomr,bbp,vvxts,vvxp0,vvxpl,iddp(ipp),icz,7))
        else
         vplcpe=max(vplcp,vpht(ipp,1)
     *   +qgfani(1.d0/xpomr,bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,13))
         vplcng=max(vplcpe,vpht(ipp,1)
     *   +qgfani(1.d0/xpomr,bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,12))
        endif 
          
        if(ityp.eq.0)then
         aks=qgran(b10)*vpac0(ipp)
         if(aks.le.vplcp.or.xpomr*sgap**2.gt..9d0)then  
          itypr0(i)=6        !single cut Pomeron
         elseif(aks.lt.vpacpe)then  
          itypr0(i)=-1       !'fan' (cut Pomeron end)
         elseif(aks.lt.vpacng)then  
          itypr0(i)=2        !'fan' (>1 cut Poms at the end)
         endif
        elseif(ityp.eq.2)then
         aks=qgran(b10)*vpacng
         if(aks.le.vplcp.or.xpomr*sgap**2.gt..9d0)then  
          itypr0(i)=6        !single cut Pomeron
         elseif(aks.lt.vpacpe)then  
          itypr0(i)=-1       !'fan' (cut Pomeron end)
         endif
        elseif(ityp.eq.4)then
         aks=qgran(b10)*vplcng
         if(aks.le.vplcp.or.xpomr*sgap**2.gt..9d0)then  
          itypr0(i)=6        !single cut Pomeron
         elseif(aks.gt.vplcpe.or.xpomr*sgap**3.gt..9d0)then  
          itypr0(i)=7        !'leg' (>1 cut Poms at the end)
         endif
        endif
        
        if(itypr0(i).eq.6)then        !single cut Pomeron
         nppr=nppr+1
         if(qgran(b10).le.vplcps/vplcp.or.xpomr*s2min.gt..9d0)then
          idpompi(nppr)=0
         else
          idpompi(nppr)=1
         endif
         xpompi(nppr)=xpomr
         ipompi(nppr)=ipp
         bpompi(nppr,1)=bbp
         bpompi(nppr,2)=(xb(it,1)-xxp)**2+(xb(it,2)-yyp)**2
         if(debug.ge.4)write (moniou,209)nppr,ipp,bbp,xpompi(nppr)
         nppr0=nppr0-1
         if(nppr0.ge.i)then
          do l=i,nppr0
           ippr0(l)=ippr0(l+1)
           itypr0(l)=itypr0(l+1)
          enddo
         endif
         i=i-1
        endif
       endif
       if(i.lt.nppr0)goto 7
      endif
         
      vvxp=vvxps
      vvxpl=vvxpls
      vvxp0=vvxp0s
      vvxts=vvxt
      vvxtls=vvxtl
      vvxt0s=vvxt0
      if(nptg0.ne.0)then
       i=0
8      i=i+1
       ityt=itytg0(i)
       if(ityt.eq.0.or.ityt.eq.2.or.ityt.eq.4)then
        itt=iptg0(i)
        bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
        vvxt=0.d0
        vvxtl=0.d0
        vvxt0=0.d0
        if(ia(2).gt.1)then
         do l=1,ia(2)
          if(l.lt.itt)then
           vvxtl=vvxtl+vtac(l)
          elseif(l.gt.itt)then
           vvxt=vvxt+vtac(l)
           vvxt0=vvxt0+vtac0(l)
          endif
         enddo
        endif
        vvxt=1.d0-exp(-vvxt)
        vvxtl=1.d0-exp(-vvxtl)
        vvxt0=1.d0-exp(-vvxt0)
        vvxps=1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*(1.d0-vvxtl)*exp(-vpac(ip))
        vtlcps=qgfani(xpomr*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,16)
        vtlcp=max(vtlcps,vtht(itt,2)
     *  +qgfani(xpomr*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,14))
        if(ityt.ne.4)then
         vtacpe=max(vtlcp,vtht(itt,1)
     *   +qgfani(xpomr*scm,bbt,vvxps,vvxt0,vvxtl,iddt(itt),2,8))
         vtacng=max(vtacpe,vtht(itt,1)
     *   +qgfani(xpomr*scm,bbt,vvxps,vvxt0,vvxtl,iddt(itt),2,7))
        else
         vtlcpe=max(vtlcp,vtht(itt,1)
     *   +qgfani(xpomr*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,13))
         vtlcng=max(vtlcpe,vtht(itt,1)
     *   +qgfani(xpomr*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,12))
        endif 
          
        if(ityt.eq.0)then
         aks=qgran(b10)*vtac0(itt)
         if(aks.le.vtlcp.or.xpomr*scm.lt.1.1d0*sgap**2)then  
          itytg0(i)=6        !single cut Pomeron
         elseif(aks.lt.vtacpe)then  
          itytg0(i)=-1       !'fan' (cut Pomeron end)
         elseif(aks.lt.vtacng)then  
          itytg0(i)=2        !'fan' (>1 cut Poms at the end)
         endif
        elseif(ityt.eq.2)then
         aks=qgran(b10)*vtacng
         if(aks.le.vtlcp.or.xpomr*scm.lt.1.1d0*sgap**2)then  
          itytg0(i)=6        !single cut Pomeron
         elseif(aks.lt.vtacpe)then  
          itytg0(i)=-1       !'fan' (cut Pomeron end)
         endif
        elseif(ityt.eq.4)then
         aks=qgran(b10)*vtlcng
         if(aks.le.vtlcp.or.xpomr*scm.lt.1.1d0*sgap**2)then  
          itytg0(i)=6
         elseif(aks.gt.vtlcpe.or.xpomr*scm.lt.1.1d0*sgap**3)then
          itytg0(i)=7        !'leg' (>1 cut Poms at the end)
         endif
        endif
        
        if(itytg0(i).eq.6)then        !single cut Pomeron
         nptg=nptg+1
         if(qgran(b10).le.vtlcps/vtlcp.or.xpomr*scm.lt.1.1d0*s2min)then
          idpomti(nptg)=0
         else
          idpomti(nptg)=1
         endif
         xpomti(nptg)=xpomr
         ipomti(nptg)=itt
         bpomti(nptg,1)=bbt
         bpomti(nptg,2)=(xa(ip,1)+bcoll-xxp)**2+(xa(ip,2)-yyp)**2
         if(debug.ge.4)write (moniou,217)nptg,itt,bbt,xpomti(nptg)
         nptg0=nptg0-1
         if(nptg0.ge.i)then
          do l=i,nptg0
           iptg0(l)=iptg0(l+1)
           itytg0(l)=itytg0(l+1)
          enddo
         endif
         i=i-1
        endif
       endif
       if(i.lt.nptg0)goto 8
      endif            
      vvxt=vvxts
      vvxtl=vvxtls
      vvxt0=vvxt0s
      
      if((jt.eq.13.or.jt.eq.15.or.jt.eq.17.or.(jt.eq.1.or.jt.eq.4
     *.or.jt.eq.7).and.qgran(b10).lt.2.d0*vpacq(ip))
     *.and.iqvp(ip).eq.0)then
       nppr=nppr+1
       idpompi(nppr)=2
       xpompi(nppr)=xpomr
       ipompi(nppr)=ip
       bpompi(nppr,1)=bbpr
       bpompi(nppr,2)=bbtg
       iqvp(ip)=1
       if(debug.ge.4)write (moniou,209)nppr,ip,bbpr,xpompi(nppr)
       
       if(ip.lt.ia(1))then
        do ipp=ip+1,ia(1)
         if(qgran(b10).lt.2.d0*vpacq(ipp).and.iqvp(ipp).eq.0)then
          nppr=nppr+1
          idpompi(nppr)=2
          xpompi(nppr)=xpomr
          ipompi(nppr)=ipp
          bbp=(xa(ipp,1)+bcoll-xxp)**2+(xa(ipp,2)-yyp)**2
          bpompi(nppr,1)=bbp
          bpompi(nppr,2)=(xb(it,1)-xxp)**2+(xb(it,2)-yyp)**2
          iqvp(ipp)=1
          if(debug.ge.4)write (moniou,209)nppr,ipp,bbp,xpompi(nppr)
         endif
        enddo
       endif
      endif
      
      if((jt.eq.14.or.jt.eq.16.or.jt.eq.18
     *.or.(jt.eq.1.or.jt.eq.3.or.jt.eq.8)
     *.and.qgran(b10).lt.2.d0*vtacq(it)).and.iqvt(it).eq.0)then
       nptg=nptg+1
       idpomti(nptg)=2
       xpomti(nptg)=xpomr
       ipomti(nptg)=it
       bpomti(nptg,1)=bbtg
       bpomti(nptg,2)=bbpr
       iqvt(it)=1
       if(debug.ge.4)write (moniou,217)nptg,it,bbtg,xpomti(nptg)

       if(it.lt.ia(2))then
        do itt=it+1,ia(2)
         if(qgran(b10).lt.2.d0*vtacq(itt).and.iqvt(itt).eq.0)then
          nptg=nptg+1
          idpomti(nptg)=2
          xpomti(nptg)=xpomr
          ipomti(nptg)=itt
          bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
          bpomti(nptg,1)=bbt
          bpomti(nptg,2)=(xa(ip,1)+bcoll-xxp)**2+(xa(ip,2)-yyp)**2
          iqvt(itt)=1
          if(debug.ge.4)write (moniou,217)nptg,itt,bbt,xpomti(nptg)
         endif
        enddo
       endif
      endif

c zigzag contributions
      if(jt.eq.1.and.xpomr*sgap**2.lt..9d0)then
       vvxps=1.d0-((1.d0-vvxp)*(1.d0-vvxpl))**2*exp(-2.d0*vpac(ip))
       vzzp=2.d0*max(0.d0,vpact(ip)-vpac0(ip))
       vzzt=(1.d0-exp(-vtac0(it)))*(1.d0-vvxt0)*vvxps
     * +min(0.d0,1.d0-exp(-vtact(it))-vtact(it))
     * -min(0.d0,1.d0-exp(-vtac0(it))-vtac0(it))
     * +vvxt0*(1.d0-exp(-vtac0(it)))-vvxtt*(1.d0-exp(-vtact(it)))
     
       wzzp=vzzp*max(0.d0,vzzt)
       nzzp=npgen(wzzp/vv(1),0,npmax)
      else
       nzzp=0
      endif
      
      if(jt.eq.1.and.xpomr*scm.gt.1.1d0*sgap**2)then
       vvxts=1.d0-((1.d0-vvxt)*(1.d0-vvxtl))**2*exp(-2.d0*vtac(it))
       vzzt=2.d0*max(0.d0,vtact(it)-vtac0(it))
       vzzp=(1.d0-exp(-vpac0(ip)))*(1.d0-vvxp0)*vvxts
     * +min(0.d0,1.d0-exp(-vpact(ip))-vpact(ip))
     * -min(0.d0,1.d0-exp(-vpac0(ip))-vpac0(ip))
     * +vvxp0*(1.d0-exp(-vpac0(ip)))-vvxpt*(1.d0-exp(-vpact(ip)))
     
       wzzt=vzzt*max(0.d0,vzzp)
       nzzt=npgen(wzzt/vv(1),0,npmax)
      else
       nzzt=0
      endif
      
      if(nzzp.ne.0)then
       bpm=(xa(ip,1)+bcoll-xxp)**2+(xa(ip,2)-yyp)**2
       xpomr0=dsqrt(xpomr)
       rp1=(rq(iddp(ip),icz)-alfp*dlog(xpomr0))*4.d0*.0389d0  
       rp2=alfp*dlog(xpomr0/xpomr)*4.d0*.0389d0 
       rp0=rp1*rp2/(rp1+rp2)
       bbp=bpm*(rp1/(rp1+rp2))**2 
       bbi=bpm*(rp2/(rp1+rp2))**2 
       call qgbdef(bbp,bbi,xa(ip,1)+bcoll,xa(ip,2),xxp,yyp,xxp0,yyp0,1)
       call qgfdf(xxp0,yyp0,xpomr0,vpac,vtac,vpht,vtht,vpacq,vtacq
     * ,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)

       sumcp0=0.d0
       sumcpt=0.d0
       sumup=0.d0
       do i=1,ia(1)
        sumup=sumup+vpac(i)
       enddo
       vvxs=(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
       do i=1,ia(1)-ip+1
        ipp=ia(1)-i+1
        bbpi=(xa(ipp,1)+bcoll-xxp0)**2+(xa(ipp,2)-yyp0)**2
        sumup=sumup-vpac(ipp)
        vps=qgfani(1.d0/xpomr0,bbpi,1.d0-vvxs*exp(-sumup)
     *  ,1.d0-exp(-sumcp0),1.d0-exp(-sumup),iddp(ipp),icz,3)
        vpac0(ipp)=max(vps,vpht(ipp,1)
     *  +qgfani(1.d0/xpomr0,bbpi,1.d0-vvxs*exp(-sumup)
     *  ,1.d0-exp(-sumcp0),1.d0-exp(-sumup),iddp(ipp),icz,6))
        vpac0(ipp)=min(vpac(ipp),vpac0(ipp))
        vpact(ipp)=max(vpac(ipp),vpht(ipp,1)
     *  +qgfani(1.d0/xpomr0,bbpi,1.d0-vvxs*exp(-sumup)
     *  ,1.d0-exp(-sumcpt),1.d0-exp(-sumup),iddp(ipp),icz,9))
        if(ipp.gt.ip)then
         sumcp0=sumcp0+vpac0(ipp)
         sumcpt=sumcpt+vpact(ipp)
        endif
       enddo
       vvxpt=1.d0-exp(-sumcpt)
       vvxp0=1.d0-exp(-sumcp0)

       sumut=0.d0
       sumct0=0.d0
       do i=1,ia(2)
        sumut=sumut+vtac(i)
       enddo
       vvxs=(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
       do i=1,ia(2)-it+1
        itt=ia(2)-i+1
        bbti=(xb(itt,1)-xxp0)**2+(xb(itt,2)-yyp0)**2
        sumut=sumut-vtac(itt)
        vts=qgfani(xpomr0*scm,bbti,1.d0-vvxs*exp(-sumut)
     *  ,1.d0-exp(-sumct0),1.d0-exp(-sumut),iddt(itt),2,3)
        vtac0(itt)=max(vts,vtht(itt,1)
     *  +qgfani(xpomr0*scm,bbti,1.d0-vvxs*exp(-sumut)
     *  ,1.d0-exp(-sumct0),1.d0-exp(-sumut),iddt(itt),2,6))
        vtac0(itt)=min(vtac(itt),vtac0(itt))
        if(itt.gt.it)sumct0=sumct0+vtac0(itt)
       enddo
       vvxt0=1.d0-exp(-sumct0)
       viu=qgpini(xpomr0/xpomr,bbi,0.d0,0.d0,2)

       gb0=(1.d0-exp(-viu))*((1.d0-exp(-vpac0(ip)))*(1.d0-vvxp0)
     * *(1.d0-exp(-2.d0*vtac(it))*((1.d0-vvxt)*(1.d0-vvxtl))**2)
     * +min(0.d0,1.d0-exp(-vpact(ip))-vpact(ip))
     * -min(0.d0,1.d0-exp(-vpac0(ip))-vpac0(ip))
     * +vvxp0*(1.d0-exp(-vpac0(ip)))-vvxpt*(1.d0-exp(-vpact(ip))))
       if(gb0.le.0.d0)then
        nzzp=0
        goto 40
       endif
       
       do in=1,nzzp
        nrej=0
9       xpomri=(xpomr*sgap**2)**qgran(b10)/sgap
        rp1=(rq(iddp(ip),icz)-alfp*dlog(xpomri))*4.d0*.0389d0  
        rp2=alfp*dlog(xpomri/xpomr)*4.d0*.0389d0 
        rp=rp1*rp2/(rp1+rp2)
        z=qgran(b10)
        phi=pi*qgran(b10)
        b0=dsqrt(-rp*dlog(z))
        bbp=(dsqrt(bpm)*rp1/(rp1+rp2)+b0*cos(phi))**2+(b0*sin(phi))**2 
        bbi=(dsqrt(bpm)*rp2/(rp1+rp2)-b0*cos(phi))**2+(b0*sin(phi))**2 
        call qgbdef(bbp,bbi,xa(ip,1)+bcoll,xa(ip,2),xxp,yyp
     *  ,xxi,yyi,int(1.5d0+qgran(b10)))   !coordinates for the vertex
        call qgfdf(xxi,yyi,xpomri,vpac,vtac,vpht,vtht,vpacq,vtacq
     *  ,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)
     
        sumcp0=0.d0
        sumcpt=0.d0
        sumup=0.d0
        do i=1,ia(1)
         sumup=sumup+vpac(i)
        enddo
        vvxs=(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
        do i=1,ia(1)-ip+1
         ipp=ia(1)-i+1
         bbpi=(xa(ipp,1)+bcoll-xxi)**2+(xa(ipp,2)-yyi)**2
         sumup=sumup-vpac(ipp)
         vps=qgfani(1.d0/xpomri,bbpi,1.d0-vvxs*exp(-sumup)
     *   ,1.d0-exp(-sumcp0),1.d0-exp(-sumup),iddp(ipp),icz,3)
         vpac0(ipp)=max(vps,vpht(ipp,1)
     *   +qgfani(1.d0/xpomri,bbpi,1.d0-vvxs*exp(-sumup)
     *   ,1.d0-exp(-sumcp0),1.d0-exp(-sumup),iddp(ipp),icz,6))
         vpac0(ipp)=min(vpac(ipp),vpac0(ipp))
         vpact(ipp)=max(vpac(ipp),vpht(ipp,1)
     *   +qgfani(1.d0/xpomri,bbpi,1.d0-vvxs*exp(-sumup)
     *   ,1.d0-exp(-sumcpt),1.d0-exp(-sumup),iddp(ipp),icz,9))
         if(ipp.gt.ip)then
          sumcp0=sumcp0+vpac0(ipp)
          sumcpt=sumcpt+vpact(ipp)
         endif
        enddo
        vvxpt=1.d0-exp(-sumcpt)
        vvxp0=1.d0-exp(-sumcp0)
       
        sumut=0.d0
        sumct0=0.d0
        do i=1,ia(2)
         sumut=sumut+vtac(i)
        enddo
        vvxs=(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
        do i=1,ia(2)-it+1
         itt=ia(2)-i+1
         bbti=(xb(itt,1)-xxi)**2+(xb(itt,2)-yyi)**2
         sumut=sumut-vtac(itt)
         vts=qgfani(xpomri*scm,bbti,1.d0-vvxs*exp(-sumut)
     *   ,1.d0-exp(-sumct0),1.d0-exp(-sumut),iddt(itt),2,3)
         vtac0(itt)=max(vts,vtht(itt,1)
     *   +qgfani(xpomri*scm,bbti,1.d0-vvxs*exp(-sumut)
     *   ,1.d0-exp(-sumct0),1.d0-exp(-sumut),iddt(itt),2,6))
         vtac0(itt)=min(vtac(itt),vtac0(itt))
         if(itt.gt.it)sumct0=sumct0+vtac0(itt)
        enddo
        vvxt0=1.d0-exp(-sumct0)
        viu=qgpini(xpomri/xpomr,bbi,0.d0,0.d0,2)
        vim=min(viu,qgpini(xpomri/xpomr,bbi,0.d0,0.d0,14))

        gb=(1.d0-exp(-viu))*((1.d0-exp(-vpac0(ip)))*(1.d0-vvxp0)
     *  *(1.d0-exp(-2.d0*vtac(it))*((1.d0-vvxt)*(1.d0-vvxtl))**2)
     *  +min(0.d0,1.d0-exp(-vpact(ip))-vpact(ip))
     *  -min(0.d0,1.d0-exp(-vpac0(ip))-vpac0(ip))
     *  +vvxp0*(1.d0-exp(-vpac0(ip)))-vvxpt*(1.d0-exp(-vpact(ip))))

        gb=gb/gb0/z*rp/rp0/30.d0
        nrej=nrej+1
        if(qgran(b10).gt.gb.and.nrej.lt.200)goto 9
                        
        if(qgran(b10).le.(1.d0-exp(-2.d0*vim))/2.d0/(1.d0-exp(-viu)))
     *  then
         ninc=npgen(2.d0*vim,1,npmax)
         v1p=min(vim,qgpini(xpomri/xpomr,bbi,0.d0,0.d0,8))
         nzpi=0
         do i=1,ninc
          if(qgran(b10).le.v1p/vim.or.xpomri/xpomr.lt.1.1d0*sgap**2)
     *    then
           npin=npin+1
           if(npin.gt.npmax)then
            iret=1
            goto 31
           endif
           xpomim(npin)=1.d0/xpomr/scm
           xpomip(npin)=xpomri
           bpomim(npin)=bbi
           if(debug.ge.4)write (moniou,211)npin,xpomip(npin)
     *     ,xpomim(npin),bpomim(npin)
          else
           nzpi=nzpi+1
           if(nzpi.gt.npmax)then
            iret=1
            goto 31
           endif
           xpmmin(nzpi)=xpomr
           xpmmax(nzpi)=xpomri
           bbm(nzpi)=bbi
           ityppm(nzpi)=1
           iorm(nzpi)=1
          endif
         enddo
        endif
       enddo   !in=1,nzzp
      endif    !nzzp.ne.0

40    if(nzzt.ne.0)then
       btm=(xb(it,1)-xxp)**2+(xb(it,2)-yyp)**2
       xpomr0=dsqrt(xpomr/scm)
       rp1=(rq(iddt(it),2)+alfp*dlog(xpomr0*scm))*4.d0*.0389d0  
       rp2=alfp*dlog(xpomr/xpomr0)*4.d0*.0389d0 
       rp0=rp1*rp2/(rp1+rp2)
       bbt=btm*(rp1/(rp1+rp2))**2 
       bbi=btm*(rp2/(rp1+rp2))**2 
       call qgbdef(bbt,bbi,xb(it,1),xb(it,2),xxp,yyp,xxp0,yyp0,1)     
       call qgfdf(xxp0,yyp0,xpomr0,vpac,vtac,vpht,vtht,vpacq,vtacq
     * ,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)

       sumct0=0.d0
       sumctt=0.d0
       sumut=0.d0
       do i=1,ia(2)
        sumut=sumut+vtac(i)
       enddo
       vvxs=(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
       do i=1,ia(2)-it+1
        itt=ia(2)-i+1
        bbti=(xb(itt,1)-xxp0)**2+(xb(itt,2)-yyp0)**2
        sumut=sumut-vtac(itt)
        vts=qgfani(xpomr0*scm,bbti,1.d0-vvxs*exp(-sumut)
     *  ,1.d0-exp(-sumct0),1.d0-exp(-sumut),iddt(itt),2,3)
        vtac0(itt)=max(vts,vtht(itt,1)
     *  +qgfani(xpomr0*scm,bbti,1.d0-vvxs*exp(-sumut)
     *  ,1.d0-exp(-sumct0),1.d0-exp(-sumut),iddt(itt),2,6))
        vtac0(itt)=min(vtac(itt),vtac0(itt))
        vtact(itt)=max(vtac(itt),vtht(itt,1)
     *  +qgfani(xpomr0*scm,bbti,1.d0-vvxs*exp(-sumut)
     *  ,1.d0-exp(-sumctt),1.d0-exp(-sumut),iddt(itt),2,9))
        if(itt.gt.it)then
         sumct0=sumct0+vtac0(itt)
         sumctt=sumctt+vtact(itt)
        endif
       enddo
       vvxtt=1.d0-exp(-sumctt)
       vvxt0=1.d0-exp(-sumct0)

       sumcp0=0.d0
       sumup=0.d0
       do i=1,ia(1)
        sumup=sumup+vpac(i)
       enddo
       vvxs=(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
       do i=1,ia(1)-ip+1
        ipp=ia(1)-i+1
        bbpi=(xa(ipp,1)+bcoll-xxp0)**2+(xa(ipp,2)-yyp0)**2
        sumup=sumup-vpac(ipp)
        vps=qgfani(1.d0/xpomr0,bbpi,1.d0-vvxs*exp(-sumup)
     *  ,1.d0-exp(-sumcp0),1.d0-exp(-sumup),iddp(ipp),icz,3)
        vpac0(ipp)=max(vps,vpht(ipp,1)
     *  +qgfani(1.d0/xpomr0,bbpi,1.d0-vvxs*exp(-sumup)
     *  ,1.d0-exp(-sumcp0),1.d0-exp(-sumup),iddp(ipp),icz,6))
        vpac0(ipp)=min(vpac(ipp),vpac0(ipp))
        if(ipp.gt.ip)sumcp0=sumcp0+vpac0(ipp)
       enddo
       vvxp0=1.d0-exp(-sumcp0)
       viu=qgpini(xpomr/xpomr0,bbi,0.d0,0.d0,2)

       gb0=(1.d0-exp(-viu))*((1.d0-exp(-vtac0(it)))*(1.d0-vvxt0)
     * *(1.d0-exp(-2.d0*vpac(ip))*((1.d0-vvxp)*(1.d0-vvxpl))**2)
     * +min(0.d0,1.d0-exp(-vtact(it))-vtact(it))
     * -min(0.d0,1.d0-exp(-vtac0(it))-vtac0(it))
     * +vvxt0*(1.d0-exp(-vtac0(it)))-vvxtt*(1.d0-exp(-vtact(it))))
       if(gb0.le.0.d0)then
        nzzt=0
        goto 41
       endif
       
       do in=1,nzzt
        nrej=0
10      xpomri=xpomr/sgap/(xpomr*scm/sgap**2)**qgran(b10)
        rp1=(rq(iddt(it),2)+alfp*dlog(xpomri*scm))*4.d0*.0389d0  
        rp2=alfp*dlog(xpomr/xpomri)*4.d0*.0389d0 
        rp=rp1*rp2/(rp1+rp2)
        z=qgran(b10)
        phi=pi*qgran(b10)
        b0=dsqrt(-rp*dlog(z))
        bbt=(dsqrt(btm)*rp1/(rp1+rp2)+b0*cos(phi))**2+(b0*sin(phi))**2 
        bbi=(dsqrt(btm)*rp2/(rp1+rp2)-b0*cos(phi))**2+(b0*sin(phi))**2 
        call qgbdef(bbt,bbi,xb(it,1),xb(it,2),xxp,yyp,xxi,yyi     
     *  ,int(1.5d0+qgran(b10)))   !coordinates for the vertex
        call qgfdf(xxi,yyi,xpomri,vpac,vtac,vpht,vtht,vpacq,vtacq
     *  ,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)
     
        sumct0=0.d0
        sumctt=0.d0
        sumut=0.d0
        do i=1,ia(2)
         sumut=sumut+vtac(i)
        enddo
        vvxs=(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
        do i=1,ia(2)-it+1
         itt=ia(2)-i+1
         bbti=(xb(itt,1)-xxi)**2+(xb(itt,2)-yyi)**2
         sumut=sumut-vtac(itt)
         vts=qgfani(xpomri*scm,bbti,1.d0-vvxs*exp(-sumut)
     *   ,1.d0-exp(-sumct0),1.d0-exp(-sumut),iddt(itt),2,3)
         vtac0(itt)=max(vts,vtht(itt,1)
     *   +qgfani(xpomri*scm,bbti,1.d0-vvxs*exp(-sumut)
     *   ,1.d0-exp(-sumct0),1.d0-exp(-sumut),iddt(itt),2,6))
         vtac0(itt)=min(vtac(itt),vtac0(itt))
         vtact(itt)=max(vtac(itt),vtht(itt,1)
     *   +qgfani(xpomri*scm,bbti,1.d0-vvxs*exp(-sumut)
     *   ,1.d0-exp(-sumctt),1.d0-exp(-sumut),iddt(itt),2,9))
         if(itt.gt.it)then
          sumct0=sumct0+vtac0(itt)
          sumctt=sumctt+vtact(itt)
         endif
        enddo
        vvxtt=1.d0-exp(-sumctt)
        vvxt0=1.d0-exp(-sumct0)       
       
        sumcp0=0.d0
        sumup=0.d0
        do i=1,ia(1)
         sumup=sumup+vpac(i)
        enddo
        vvxs=(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
        do i=1,ia(1)-ip+1
         ipp=ia(1)-i+1
         bbpi=(xa(ipp,1)+bcoll-xxi)**2+(xa(ipp,2)-yyi)**2
         sumup=sumup-vpac(ipp)
         vps=qgfani(1.d0/xpomri,bbpi,1.d0-vvxs*exp(-sumup)
     *   ,1.d0-exp(-sumcp0),1.d0-exp(-sumup),iddp(ipp),icz,3)
         vpac0(ipp)=max(vps,vpht(ipp,1)
     *   +qgfani(1.d0/xpomri,bbpi,1.d0-vvxs*exp(-sumup)
     *   ,1.d0-exp(-sumcp0),1.d0-exp(-sumup),iddp(ipp),icz,6))
         vpac0(ipp)=min(vpac(ipp),vpac0(ipp))
         if(ipp.gt.ip)sumcp0=sumcp0+vpac0(ipp)
        enddo
        vvxp0=1.d0-exp(-sumcp0)      
        viu=qgpini(xpomr/xpomri,bbi,0.d0,0.d0,2)
        vim=min(viu,qgpini(xpomr/xpomri,bbi,0.d0,0.d0,14))

        gb=(1.d0-exp(-viu))*((1.d0-exp(-vtac0(it)))*(1.d0-vvxt0)
     *  *(1.d0-exp(-2.d0*vpac(ip))*((1.d0-vvxp)*(1.d0-vvxpl))**2)
     *  +min(0.d0,1.d0-exp(-vtact(it))-vtact(it))
     *  -min(0.d0,1.d0-exp(-vtac0(it))-vtac0(it))
     *  +vvxt0*(1.d0-exp(-vtac0(it)))-vvxtt*(1.d0-exp(-vtact(it))))

        gb=gb/gb0/z*rp/rp0/30.d0
        nrej=nrej+1
        if(qgran(b10).gt.gb.and.nrej.lt.200)goto 10
                        
        if(qgran(b10).le.(1.d0-exp(-2.d0*vim))/2.d0/(1.d0-exp(-viu)))
     *  then
         ninc=npgen(2.d0*vim,1,npmax)
         v1p=min(vim,qgpini(xpomr/xpomri,bbi,0.d0,0.d0,8))
         do i=1,ninc
          if(qgran(b10).le.v1p/vim.or.xpomr/xpomri.lt.1.1d0*sgap**2)
     *    then
           npin=npin+1
           if(npin.gt.npmax)then
            iret=1
            goto 31
           endif
           xpomim(npin)=1.d0/xpomri/scm
           xpomip(npin)=xpomr
           bpomim(npin)=bbi
           if(debug.ge.4)write (moniou,211)npin,xpomip(npin)
     *     ,xpomim(npin),bpomim(npin)
          else
           nzpi=nzpi+1
           if(nzpi.gt.npmax)then
            iret=1
            goto 31
           endif
           xpmmin(nzpi)=xpomri
           xpmmax(nzpi)=xpomr
           bbm(nzpi)=bbi
           ityppm(nzpi)=1
           iorm(nzpi)=-1
          endif
         enddo
        endif
       enddo   !in=1,nzzt
      endif    !nzzt.ne.0

41    if(nzpi.ne.0)then
11     xpimin=xpmmin(nzpi)
       xpimax=xpmmax(nzpi)
       bbi=bbm(nzpi)
       itypom=ityppm(nzpi)
       ipor=iorm(nzpi)
       if(xpimax/xpimin.lt.sgap**2)stop'qg3pdf:syz<sgap^2!?'
       
       bbzu=bbi/4.d0
       bbzd=bbzu
       rp0=alfp*dlog(xpimax/xpimin)*.0389d0
       syp=dsqrt(xpimax/xpimin)
       viu=qgpini(syp,bbzu,0.d0,0.d0,2)
       vimu=min(viu,qgpini(syp,bbzu,0.d0,0.d0,14))
       vingu=min(vimu,qgpini(syp,bbzu,0.d0,0.d0,15))
       vipeu=min(vingu,qgpini(syp,bbzu,0.d0,0.d0,17))
       v1pu=min(vipeu,qgpini(syp,bbzu,0.d0,0.d0,8))

       vid=viu
       vimd=vimu
       vingd=vingu
       viped=vipeu
       v1pd=v1pu
       
       if(itypom.eq.2)then    !all except 1P-end
        vv1=0.d0
       else
        vv1=4.d0*v1pd*exp(-2.d0*vimd)
     *  *max(0.d0,1.d0-exp(-viu)-vipeu*exp(-2.d0*vimu))
       endif
       if(itypom.eq.3)then    !no rapgap at the end
        vv2=0.d0
       else
        vv2=(2.d0*min(0.d0,1.d0-exp(-vid)-vid)
     *  -min(0.d0,1.d0-exp(-2.d0*vimd)-2.d0*vimd)+2.d0*(vid-vimd))
     *  *(max(0.d0,1.d0-exp(-2.d0*vimu)*(1.d0+2.d0*vimu))
     *  +2.d0*vingu*exp(-2.d0*vimu))
       endif
       vv3=2.d0*max(0.d0,1.d0-exp(-2.d0*vimd)*(1.d0+2.d0*vimd))
     * *(1.d0-exp(-viu))
       gbz0=(vv1+vv2+vv3)

       nrej=0
12     xpi=xpimin*sgap*(xpimax/xpimin/sgap**2)**qgran(b10)
       rpd=alfp*dlog(xpi/xpimin)*4.d0*.0389d0  
       rpu=alfp*dlog(xpimax/xpi)*4.d0*.0389d0  
       rp=rpd*rpu/(rpd+rpu)
       z=qgran(b10)
       phi=pi*qgran(b10)
       b0=dsqrt(-rp*dlog(z))
       bbzu=(dsqrt(bbi)*rpu/(rpd+rpu)+b0*cos(phi))**2
     * +(b0*sin(phi))**2
       bbzd=(dsqrt(bbi)*rpd/(rpd+rpu)-b0*cos(phi))**2
     * +(b0*sin(phi))**2

       syu=xpimax/xpi
       viu=qgpini(syu,bbzu,0.d0,0.d0,2)
       vimu=min(viu,qgpini(syu,bbzu,0.d0,0.d0,14))
       vingu=min(vimu,qgpini(syu,bbzu,0.d0,0.d0,15))
       vipeu=min(vingu,qgpini(syu,bbzu,0.d0,0.d0,17))
       v1pu=min(vipeu,qgpini(syu,bbzu,0.d0,0.d0,8))

       syd=xpi/xpimin
       vid=qgpini(syd,bbzd,0.d0,0.d0,2)
       vimd=min(vid,qgpini(syd,bbzd,0.d0,0.d0,14))
       vingd=min(vimd,qgpini(syd,bbzd,0.d0,0.d0,15))
       viped=min(vingd,qgpini(syd,bbzd,0.d0,0.d0,17))
       v1pd=min(viped,qgpini(syd,bbzd,0.d0,0.d0,8))
       
       if(itypom.eq.2)then
        vv1=0.d0
       elseif(ipor.eq.1)then
        vv1=4.d0*v1pd*exp(-2.d0*vimd)
     *  *max(0.d0,1.d0-exp(-viu)-vipeu*exp(-2.d0*vimu))
       else
        vv1=4.d0*v1pu*exp(-2.d0*vimu)
     *  *max(0.d0,1.d0-exp(-vid)-viped*exp(-2.d0*vimd))
       endif
       if(itypom.eq.3)then
        vv2=0.d0
       elseif(ipor.eq.1)then
        vv2=(2.d0*min(0.d0,1.d0-exp(-vid)-vid)
     *  -min(0.d0,1.d0-exp(-2.d0*vimd)-2.d0*vimd)+2.d0*(vid-vimd))
     *  *(max(0.d0,1.d0-exp(-2.d0*vimu)*(1.d0+2.d0*vimu))
     *  +2.d0*vingu*exp(-2.d0*vimu))
       else
        vv2=(2.d0*min(0.d0,1.d0-exp(-viu)-viu)
     *  -min(0.d0,1.d0-exp(-2.d0*vimu)-2.d0*vimu)+2.d0*(viu-vimu))
     *  *(max(0.d0,1.d0-exp(-2.d0*vimd)*(1.d0+2.d0*vimd))
     *  +2.d0*vingd*exp(-2.d0*vimd))
       endif
       if(ipor.eq.1)then
        vv3=2.d0*max(0.d0,1.d0-exp(-2.d0*vimd)*(1.d0+2.d0*vimd))
     *  *(1.d0-exp(-viu))
       else
        vv3=2.d0*max(0.d0,1.d0-exp(-2.d0*vimu)*(1.d0+2.d0*vimu))
     *  *(1.d0-exp(-vid))
       endif
       vvt=vv1+vv2+vv3
       gbz=vvt/gbz0/z*rp/rp0   /15.
       nrej=nrej+1
       if(qgran(b10).gt.gbz.and.nrej.lt.200)goto 12

       nzpi=nzpi-1
       aks=vvt*qgran(b10)
       if(aks.lt.vv1)then
        jtd=1
       elseif(aks.lt.vv1+vv2)then
        jtd=2
       else
        jtd=3
       endif
       
       if(jtd.eq.1)then
        npin=npin+1
        if(npin.gt.npmax)then
         iret=1
         goto 31
        endif
        if(ipor.eq.1)then
         xpomim(npin)=1.d0/xpimin/scm
         xpomip(npin)=xpi
         bpomim(npin)=bbzd
        else
         xpomim(npin)=1.d0/xpi/scm
         xpomip(npin)=xpimax
         bpomim(npin)=bbzu
        endif
        if(debug.ge.4)write (moniou,211)npin,xpomip(npin)
     *  ,xpomim(npin),bpomim(npin)

       elseif(jtd.eq.3)then
        if(ipor.eq.1)then
         ninc=npgen(2.d0*vimd,2,npmax)
         w1p=v1pd/vimd
         syp=xpi/xpimin
        else
         ninc=npgen(2.d0*vimu,2,npmax)
         w1p=v1pu/vimu
         syp=xpimax/xpi
        endif
        do i=1,ninc
         if(qgran(b10).le.w1p.or.syp.lt.1.1d0*sgap**2)then
          npin=npin+1
          if(npin.gt.npmax)then
           iret=1
           goto 31
          endif
          if(ipor.eq.1)then
           xpomim(npin)=1.d0/xpimin/scm
           xpomip(npin)=xpi
           bpomim(npin)=bbzd
          else
           xpomim(npin)=1.d0/xpi/scm
           xpomip(npin)=xpimax
           bpomim(npin)=bbzu
          endif
          if(debug.ge.4)write (moniou,211)npin,xpomip(npin)
     *    ,xpomim(npin),bpomim(npin)
         else
          nzpi=nzpi+1
          if(nzpi.gt.npmax)then
           iret=1
           goto 31
          endif
          ityppm(nzpi)=1
          if(ipor.eq.1)then
           xpmmin(nzpi)=xpimin
           xpmmax(nzpi)=xpi
           bbm(nzpi)=bbzd
           iorm(nzpi)=1
          else
           xpmmin(nzpi)=xpi
           xpmmax(nzpi)=xpimax
           bbm(nzpi)=bbzu
           iorm(nzpi)=-1
          endif
         endif
        enddo
       endif

       if(ipor.eq.1)then
        v1p=v1pu
        vim=vimu
        syp=xpimax/xpi
        if(jtd.eq.1)then
         vvv=2.d0*max(0.d0,1.d0-exp(-viu)-vipeu*exp(-2.d0*vimu))
        elseif(jtd.eq.2)then
         vvv=max(0.d0,1.d0-exp(-2.d0*vimu)*(1.d0+2.d0*vimu))
     *   +2.d0*vingu*exp(-2.d0*vimu)
        elseif(jtd.eq.3)then
         vvv=2.d0*(1.d0-exp(-viu))
        else
         stop'qg3pdf:jtd?!'
        endif
       else
        v1p=v1pd
        vim=vimd
        syp=xpi/xpimin
        if(jtd.eq.1)then
         vvv=2.d0*max(0.d0,1.d0-exp(-vid)-viped*exp(-2.d0*vimd))
        elseif(jtd.eq.2)then
         vvv=max(0.d0,1.d0-exp(-2.d0*vimd)*(1.d0+2.d0*vimd))
     *   +2.d0*vingd*exp(-2.d0*vimd)
        elseif(jtd.eq.3)then
         vvv=2.d0*(1.d0-exp(-vid))
        else
         stop'qg3pdf:jtd?!'
        endif
       endif
       
       vnp=max(0.d0,1.d0-exp(-2.d0*vim)*(1.d0+2.d0*vim))
       if(qgran(b10).lt.vnp/vvv)then
        ninc=npgen(2.d0*vim,2,npmax)
        do i=1,ninc
         if(qgran(b10).le.v1p/vim.or.syp.lt.1.1d0*sgap**2)then
          npin=npin+1
          if(npin.gt.npmax)then
           iret=1
           goto 31
          endif
          if(ipor.eq.1)then
           xpomim(npin)=1.d0/xpi/scm
           xpomip(npin)=xpimax
           bpomim(npin)=bbzu
          else
           xpomim(npin)=1.d0/xpimin/scm
           xpomip(npin)=xpi
           bpomim(npin)=bbzd
          endif
          if(debug.ge.4)write (moniou,211)npin,xpomip(npin)
     *    ,xpomim(npin),bpomim(npin)
         else
          nzpi=nzpi+1
          if(nzpi.gt.npmax)then
           iret=1
           goto 31
          endif
          ityppm(nzpi)=1
          iorm(nzpi)=ipor
          if(ipor.eq.1)then
           xpmmin(nzpi)=xpi
           xpmmax(nzpi)=xpimax
           bbm(nzpi)=bbzu
          else
           xpmmin(nzpi)=xpimin
           xpmmax(nzpi)=xpi
           bbm(nzpi)=bbzd
          endif
         endif
        enddo
        
       else
        if(jtd.eq.1.and.syp.gt.1.1d0*sgap**2.and.(ipor.eq.1
     *  .and.qgran(b10).lt.2.d0*(vimu-vipeu)*exp(-2.d0*vimu)/(vvv-vnp)
     *  .or.ipor.eq.-1.and.qgran(b10).lt.2.d0*(vimd-viped)
     *  *exp(-2.d0*vimd)/(vvv-vnp)))then
         nzpi=nzpi+1
         if(nzpi.gt.npmax)then
          iret=1
          goto 31
         endif
         ityppm(nzpi)=2
         iorm(nzpi)=ipor
         if(ipor.eq.1)then
          xpmmin(nzpi)=xpi
          xpmmax(nzpi)=xpimax
          bbm(nzpi)=bbzu
         else
          xpmmin(nzpi)=xpimin
          xpmmax(nzpi)=xpi
          bbm(nzpi)=bbzd
         endif

        elseif(jtd.eq.2)then
         if(syp.lt.1.1d0*sgap**2
     *   .or.qgran(b10).lt.2.d0*v1p*exp(-2.d0*vim)/(vvv-vnp))then
          npin=npin+1
          if(npin.gt.npmax)then
           iret=1
           goto 31
          endif
          if(ipor.eq.1)then
           xpomim(npin)=1.d0/xpi/scm
           xpomip(npin)=xpimax
           bpomim(npin)=bbzu
          else
           xpomim(npin)=1.d0/xpimin/scm
           xpomip(npin)=xpi
           bpomim(npin)=bbzd
          endif
          if(debug.ge.4)write (moniou,211)npin,xpomip(npin)
     *    ,xpomim(npin),bpomim(npin)
         else
          nzpi=nzpi+1
          if(nzpi.gt.npmax)then
           iret=1
           goto 31
          endif
          ityppm(nzpi)=3
          iorm(nzpi)=ipor
          if(ipor.eq.1)then
           xpmmin(nzpi)=xpi
           xpmmax(nzpi)=xpimax
           bbm(nzpi)=bbzu
          else
           xpmmin(nzpi)=xpimin
           xpmmax(nzpi)=xpi
           bbm(nzpi)=bbzd
          endif
         endif

        elseif(jtd.eq.3.and.qgran(b10)
     *  .lt.2.d0*vim*exp(-2.d0*vim)/(vvv-vnp))then
         if(syp.lt.1.1d0*sgap**2.or.qgran(b10).lt.v1p/vim)then
          npin=npin+1
          if(npin.gt.npmax)then
           iret=1
           goto 31
          endif
          if(ipor.eq.1)then
           xpomim(npin)=1.d0/xpi/scm
           xpomip(npin)=xpimax
           bpomim(npin)=bbzu
          else
           xpomim(npin)=1.d0/xpimin/scm
           xpomip(npin)=xpi
           bpomim(npin)=bbzd
          endif
          if(debug.ge.4)write (moniou,211)npin,xpomip(npin)
     *    ,xpomim(npin),bpomim(npin)
         else
          nzpi=nzpi+1
          if(nzpi.gt.npmax)then
           iret=1
           goto 31
          endif
          ityppm(nzpi)=1
          iorm(nzpi)=ipor
          if(ipor.eq.1)then
           xpmmin(nzpi)=xpi
           xpmmax(nzpi)=xpimax
           bbm(nzpi)=bbzu
          else
           xpmmin(nzpi)=xpimin
           xpmmax(nzpi)=xpi
           bbm(nzpi)=bbzd
          endif
         endif
        endif
       endif
       if(nzpi.gt.0)goto 11
      endif   !end of 'zigzags'

      call qgfdf(xxp,yyp,xpomr,vpac,vtac,vpht,vtht,vpacq,vtacq
     *,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)
      if((jt.eq.2.or.jt.eq.3.or.jt.eq.9.or.jt.eq.16).and.qgran(b10).gt
     *.2.d0*vvxpl/((1.d0-exp(-vpac(ip)))*(1.d0-vvxpl)+2.d0*vvxpl))then
       icdps=iddp(ip)
       do icdp=1,nfock
        iddp(ip)=icdp
        call qgfdf(xxp,yyp,xpomr,vpac,vtac,vpht,vtht,vpacq,vtacq
     *  ,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)
        wdp(icdp,ip)=(1.d0-exp(-vpac(ip)))*(1.d0-vvxpl)
       enddo
       iddp(ip)=icdps
      endif
      call qgfdf(xxp,yyp,xpomr,vpac,vtac,vpht,vtht,vpacq,vtacq
     *,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)
      if((jt.eq.2.or.jt.eq.4.or.jt.eq.10.or.jt.eq.15).and.qgran(b10).gt
     *.2.d0*vvxtl/((1.d0-exp(-vtac(it)))*(1.d0-vvxtl)+2.d0*vvxtl))then
       icdts=iddt(it)
       do icdt=1,nfock
        iddt(it)=icdt
        call qgfdf(xxp,yyp,xpomr,vpac,vtac,vpht,vtht,vpacq,vtacq
     *  ,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,it)
        wdt(icdt,it)=(1.d0-exp(-vtac(it)))*(1.d0-vvxtl)
       enddo
       iddt(it)=icdts
      endif

c projectile 'fans'      
      if(nppr0.eq.0)goto 20
      m=0
      nppm(1)=nppr0
      xpomm(1)=xpomr
      wgpm(1)=wgpr0
      xxm(1)=xxp
      yym(1)=yyp
      do i=1,nppr0
       ippm(i,1)=ippr0(i)
       itypm(i,1)=itypr0(i)
      enddo
            
13    m=m+1                                 !next level multi-Pomeron vertex
      if(m.gt.levmax)then
       iret=1
       goto 31
      endif
      ii(m)=0
14    ii(m)=ii(m)+1                         !next cut fan in the vertex
      if(ii(m).gt.nppm(m))then              !all fans at the level considered
       m=m-1                                !one level down
       if(m.eq.0)goto 20                    !all proj. fans considered 
       goto 14
      endif 
      l=ii(m)
      ipp=ippm(l,m)                         !proj. index for the leg      
      itypom=itypm(l,m)                     !type of the cut 
      bpm=(xa(ipp,1)+bcoll-xxm(m))**2+(xa(ipp,2)-yym(m))**2  !b^2 for the leg
      if(debug.ge.4)write (moniou,208)ii(m),m,ipp,bpm
      if(xpomm(m)*sgap**2.gt.1.d0)stop'xpomm(m)*sgap**2>1!'
      if(itypom.eq.4.and.xpomm(m)*sgap**3.gt.1.d0)
     *stop'4:xpomm(m)*sgap**3>1!'
      
      if(debug.ge.4)write (moniou,210)m
      xpomr0=dsqrt(xpomm(m))
      if(itypom.eq.-1)xpomr0=max(xpomr0,.9d0/sgap)
      if(itypom.eq.4)xpomr0=dsqrt(xpomm(m)/sgap)
      rp1=(rq(iddp(ipp),icz)-alfp*dlog(xpomr0))*4.d0*.0389d0  
      rp2=alfp*dlog(xpomr0/xpomm(m))*4.d0*.0389d0 
      rp0=rp1*rp2/(rp1+rp2)
      bbp=bpm*(rp1/(rp1+rp2))**2 
      bbi=bpm*(rp2/(rp1+rp2))**2 
      call qgbdef(bbp,bbi,xa(ipp,1)+bcoll,xa(ipp,2),xxm(m),yym(m)
     *,xxp0,yyp0,1)      
       
      call qgfdf(xxp0,yyp0,xpomr0,vpac,vtac,vpht,vtht,vpacq,vtacq
     *,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ipp,it)
      vvxts=1.d0-(1.d0-vvx)*(1.d0-vvxt)*exp(-vtac(it))
      viu=qgpini(xpomr0/xpomm(m),bbi,0.d0,0.d0,2)
      vim=2.d0*min(viu,qgpini(xpomr0/xpomm(m),bbi,0.d0,0.d0,14))
      if(itypom.eq.-1.or.itypom.eq.4)then         !single cut Pomeron at the end
       vvxi=1.d0-(1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt)
     * *exp(-vpac(ipp)-vtac(it))
       vip=qgpini(xpomr0/xpomm(m),bbi,vvxi,0.d0,19)
     * +qgpini(xpomr0/xpomm(m),bbi,0.d0,0.d0,8)*(exp(-vim)-1.d0)
       if(vip.lt.0.d0)vip=qgpini(xpomr0/xpomm(m),bbi,0.d0,0.d0,8)
     * *exp(-vim)
      elseif(itypom.eq.2.or.itypom.eq.7)then       !>1 cut Poms at the end
       vimp=max(0.d0,1.d0-exp(-vim)*(1.d0+vim))
      else                                         !rap-gap
       vvxpin=1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ipp))
       vvxtin=1.d0-(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
       viuu=qgpini(xpomr0/xpomm(m),bbi,vvxpin,vvxtin,20)
       viuc=max(0.d0,viuu-qgpini(xpomr0/xpomm(m),bbi,vvxpin,vvxtin,21))
       viuu=max(0.d0,viuu+min(0.d0,1.d0-exp(-viu)-viu))
       vicc=qgpini(xpomr0/xpomm(m),bbi,vvxpin,vvxtin,22)
       vicu=max(0.d0,qgpini(xpomr0/xpomm(m),bbi,vvxpin,vvxtin,23)-vicc)
       vicc=max(0.d0,vicc+.5d0*(1.d0-exp(-viu))**2-viu
     * +.5d0*(exp(2.d0*viu-vim)-1.d0)*exp(-2.d0*viu))
      endif

      vplcps=qgfani(1.d0/xpomr0,bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,16)
      vplcp=max(vplcps,vpht(ipp,2)
     *+qgfani(1.d0/xpomr0,bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,14))
      if(itypom.le.3)then
       sumup=0.d0
       vvxp0=0.d0
       do i=1,ia(1)
        sumup=sumup+vpac(i)
       enddo
       vvxs=(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
       do i=1,ia(1)-ipp+1
        ipi=ia(1)-i+1
        bbl=(xa(ipi,1)+bcoll-xxp0)**2+(xa(ipi,2)-yyp0)**2
        sumup=sumup-vpac(ipi)
        vps=qgfani(1.d0/xpomr0,bbl,1.d0-vvxs*exp(-sumup)
     *  ,1.d0-exp(-vvxp0),1.d0-exp(-sumup),iddp(ipi),icz,3)
        vpac0(ipi)=max(vps,vpht(ipi,1)
     *  +qgfani(1.d0/xpomr0,bbl,1.d0-vvxs*exp(-sumup)
     *  ,1.d0-exp(-vvxp0),1.d0-exp(-sumup),iddp(ipi),icz,6))
        vpac0(ipi)=min(vpac(ipi),vpac0(ipi))
        if(ipi.gt.ipp)vvxp0=vvxp0+vpac0(ipi)
       enddo
       vvxp0=1.d0-exp(-vvxp0)
       vpacpe=max(vplcp,vpht(ipp,1)
     * +qgfani(1.d0/xpomr0,bbp,vvxts,vvxp0,vvxpl,iddp(ipp),icz,8))
       vpacng=max(vpacpe,vpht(ipp,1)
     * +qgfani(1.d0/xpomr0,bbp,vvxts,vvxp0,vvxpl,iddp(ipp),icz,7))
      else
       vplcpe=max(vplcp,vpht(ipp,1)
     * +qgfani(1.d0/xpomr0,bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,13))
       vplcng=max(vplcpe,vpht(ipp,1)
     * +qgfani(1.d0/xpomr0,bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,12))
       vplc0=max(vplcng,vpht(ipp,1)
     * +qgfani(1.d0/xpomr0,bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,11))
       vplc=max(vplc0,vpht(ipp,1)
     * +qgfani(1.d0/xpomr0,bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,10))
      endif

      if(itypom.eq.-1)then          !'fan' (single cut Pomeron at the end)
       gb0=vip*((max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))
     * +((1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl)
     * +2.d0*((vpac0(ipp)-vpacpe)*exp(-vpac(ipp))*(1.d0-vvxp)
     * *(1.d0-vvxpl)-(vpac(ipp)-vpac0(ipp))*(1.d0-exp(-vpac(ipp))
     * *(1.d0-vvxp)*(1.d0-vvxpl)))*exp(-vpac(ipp))*(1.d0-vvxp))
     * *(1.d0-vvx)*(1.d0-vvxt)**2*(1.d0-vvxtl)*exp(-2.d0*vtac(it))
       gb0=gb0*50.d0
      elseif(itypom.eq.0)then      !'fan' (cut loop at the end - rapgap)
       gb0=((max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))
     * +((1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl)
     * +2.d0*vpacng*exp(-2.d0*vpac(ipp))*(1.d0-vvxp)**2*(1.d0-vvxpl))
     * *(vicc*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
     * -vicu*(1.d0-(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))))
     * *(1.d0-vvx)*(1.d0-vvxt)*exp(-vtac(it))
     * -2.d0*vicu*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp)-vtac(it))
     * *(1.d0-vvx)*(1.d0-vvxt)
      elseif(itypom.eq.1)then      !'fan' (uncut end - rapgap)
       gb0=((max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))
     * +((1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl)
     * +2.d0*vpacng*exp(-2.d0*vpac(ipp))*(1.d0-vvxp)**2*(1.d0-vvxpl))
     * *(viuc*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
     * +viuu*(1.d0-(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))))
     * *(1.d0-vvx)*(1.d0-vvxt)*exp(-vtac(it))
     * +2.d0*viuu*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp)-vtac(it))
     * *(1.d0-vvx)*(1.d0-vvxt)
      elseif(itypom.eq.2)then        !'fan' (>1 cut Poms at the end)
       gb0=vimp*((max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))
     * +((1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl)
     * +2.d0*(vpac0(ipp)*exp(-vpac(ipp))*(1.d0-vvxp)
     * *(1.d0-vvxpl)-(vpac(ipp)-vpac0(ipp))*(1.d0-exp(-vpac(ipp))
     * *(1.d0-vvxp)*(1.d0-vvxpl)))*exp(-vpac(ipp))*(1.d0-vvxp))
     * *(1.d0-vvx)*(1.d0-vvxt)**2*(1.d0-vvxtl)*exp(-2.d0*vtac(it))
      elseif(itypom.eq.3)then      !'fan' (cut/uncut end - rapgap)
       gb0=((max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))
     * +((1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl)
     * +2.d0*vpacng*exp(-2.d0*vpac(ipp))*(1.d0-vvxp)**2*(1.d0-vvxpl))
     * *(vicc-wgpm(m)*viuc)*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
     * *(1.d0-vvx)*(1.d0-vvxt)*exp(-vtac(it))
      elseif(itypom.eq.4)then          !'leg' (single cut Pomeron at the end)
       gb0=vip*((vplc0-vplcpe)*exp(-vpac(ipp))*(1.d0-vvxp)
     * *(1.d0-vvxpl))*exp(-vpac(ipp)-2.d0*vtac(it))*(1.d0-vvxp)
     * *(1.d0-vvx)*(1.d0-vvxt)**2*(1.d0-vvxtl)
      elseif(itypom.eq.5)then      !'leg' (cut/uncut end - rapgap)
       gb0=vplcng*exp(-2.d0*vpac(ipp)-vtac(it))
     * *(1.d0-vvxp)**2*(1.d0-vvxpl)*(1.d0-vvx)*(1.d0-vvxt)
     * *(vicc-wgpm(m)*viuc)*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
      elseif(itypom.eq.7)then      !'leg' (>1 cut Poms at the end)
       gb0=vimp*(vplc0*exp(-vpac(ipp))*(1.d0-vvxp)*(1.d0-vvxpl)
     * -(vplc-vplc0)*(1.d0-exp(-vpac(ipp))*(1.d0-vvxp)*(1.d0-vvxpl)))
     * *exp(-vpac(ipp)-2.d0*vtac(it))*(1.d0-vvxp)
     * *(1.d0-vvx)*(1.d0-vvxt)**2*(1.d0-vvxtl)
      endif
      if(gb0.le.0.d0)then
       nppr=nppr+1
       if(nppr.gt.legmax)then
        iret=1
        goto 31
       endif
       if(qgran(b10).le.vplcps/vplcp)then
        idpompi(nppr)=0
       else
        idpompi(nppr)=1
       endif
       xpompi(nppr)=xpomm(m)
       ipompi(nppr)=ipp
       bpompi(nppr,1)=bpm
       bpompi(nppr,2)=(xb(it,1)-xxm(m))**2+(xb(it,2)-yym(m))**2
       goto 14
      endif
      nrej=0

15    xpomm(m+1)=(xpomm(m)*sgap**2)**qgran(b10)/sgap
      if(itypom.eq.4)xpomm(m+1)=(xpomm(m)*sgap**3)**qgran(b10)/sgap**2
      rp1=(rq(iddp(ipp),icz)-alfp*dlog(xpomm(m+1)))*4.d0*.0389d0  
      rp2=alfp*dlog(xpomm(m+1)/xpomm(m))*4.d0*.0389d0 
      rp=rp1*rp2/(rp1+rp2)
      z=qgran(b10)
      phi=pi*qgran(b10)
      b0=dsqrt(-rp*dlog(z))
      bbp=(dsqrt(bpm)*rp1/(rp1+rp2)+b0*cos(phi))**2+(b0*sin(phi))**2
      bbi=(dsqrt(bpm)*rp2/(rp1+rp2)-b0*cos(phi))**2+(b0*sin(phi))**2
      call qgbdef(bbp,bbi,xa(ipp,1)+bcoll,xa(ipp,2),xxm(m),yym(m)
     *,xxm(m+1),yym(m+1),int(1.5d0+qgran(b10)))   !coordinates for the vertex

      call qgfdf(xxm(m+1),yym(m+1),xpomm(m+1),vpac,vtac,vpht,vtht
     *,vpacq,vtacq,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ipp,it)
      vvxts=1.d0-(1.d0-vvx)*(1.d0-vvxt)*exp(-vtac(it))
      viu=qgpini(xpomm(m+1)/xpomm(m),bbi,0.d0,0.d0,2)
      vim=2.d0*min(viu,qgpini(xpomm(m+1)/xpomm(m),bbi,0.d0,0.d0,14))
      if(itypom.eq.-1.or.itypom.eq.4)then         !single cut Pomeron at the end
       vvxi=1.d0-(1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt)
     * *exp(-vpac(ipp)-vtac(it))
       vip=qgpini(xpomm(m+1)/xpomm(m),bbi,vvxi,0.d0,19)
     * +qgpini(xpomm(m+1)/xpomm(m),bbi,0.d0,0.d0,8)*(exp(-vim)-1.d0)
       if(vip.lt.0.d0)vip=qgpini(xpomm(m+1)/xpomm(m),bbi,0.d0,0.d0,8)
     * *exp(-vim)
      elseif(itypom.eq.2.or.itypom.eq.7)then       !>1 cut Poms at the end
       vimp=max(0.d0,1.d0-exp(-vim)*(1.d0+vim))
      else                                         !rap-gap
       vvxpin=1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ipp))
       vvxtin=1.d0-(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
       viuu=qgpini(xpomm(m+1)/xpomm(m),bbi,vvxpin,vvxtin,20)
       viuc=max(0.d0,viuu-qgpini(xpomm(m+1)/xpomm(m),bbi,vvxpin,vvxtin
     * ,21))
       viuu=max(0.d0,viuu+min(0.d0,1.d0-exp(-viu)-viu))
       vicc=qgpini(xpomm(m+1)/xpomm(m),bbi,vvxpin,vvxtin,22)
       vicu=max(0.d0,qgpini(xpomm(m+1)/xpomm(m),bbi,vvxpin,vvxtin,23)
     * -vicc)
       vicc=max(0.d0,vicc+.5d0*(1.d0-exp(-viu))**2-viu
     * +.5d0*(exp(2.d0*viu-vim)-1.d0)*exp(-2.d0*viu))
      endif

      vplcps=qgfani(1.d0/xpomm(m+1),bbp,vvxts,vvxp,vvxpl,iddp(ipp)
     *,icz,16)
      vplcp=max(vplcps,vpht(ipp,2)
     *+qgfani(1.d0/xpomm(m+1),bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,14))
      if(itypom.le.3)then
       sumup=0.d0
       vvxp0=0.d0
       do i=1,ia(1)
        sumup=sumup+vpac(i)
       enddo
       vvxs=(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
       do i=1,ia(1)-ipp+1
        ipi=ia(1)-i+1
        bbl=(xa(ipi,1)+bcoll-xxm(m+1))**2+(xa(ipi,2)-yym(m+1))**2
        sumup=sumup-vpac(ipi)
        vps=qgfani(1.d0/xpomm(m+1),bbl,1.d0-vvxs*exp(-sumup)
     *  ,1.d0-exp(-vvxp0),1.d0-exp(-sumup),iddp(ipi),icz,3)
        vpac0(ipi)=max(vps,vpht(ipi,1)
     *  +qgfani(1.d0/xpomm(m+1),bbl,1.d0-vvxs*exp(-sumup)
     *  ,1.d0-exp(-vvxp0),1.d0-exp(-sumup),iddp(ipi),icz,6))
        vpac0(ipi)=min(vpac(ipi),vpac0(ipi))
        if(ipi.gt.ipp)vvxp0=vvxp0+vpac0(ipi)
       enddo
       vvxp0=1.d0-exp(-vvxp0)
       vpacpe=max(vplcp,vpht(ipp,1)
     * +qgfani(1.d0/xpomm(m+1),bbp,vvxts,vvxp0,vvxpl,iddp(ipp),icz,8))
       vpacng=max(vpacpe,vpht(ipp,1)
     * +qgfani(1.d0/xpomm(m+1),bbp,vvxts,vvxp0,vvxpl,iddp(ipp),icz,7))
      else
       vplcpe=max(vplcp,vpht(ipp,1)
     * +qgfani(1.d0/xpomm(m+1),bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,13))
       vplcng=max(vplcpe,vpht(ipp,1)
     * +qgfani(1.d0/xpomm(m+1),bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,12))
       vplc0=max(vplcng,vpht(ipp,1)
     * +qgfani(1.d0/xpomm(m+1),bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,11))
       vplc=max(vplc0,vpht(ipp,1)
     * +qgfani(1.d0/xpomm(m+1),bbp,vvxts,vvxp,vvxpl,iddp(ipp),icz,10))
      endif

      if(itypom.eq.-1)then          !'fan' (single cut Pomeron at the end)
       gb=vip*((max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))
     * +((1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl)
     * +2.d0*((vpac0(ipp)-vpacpe)*exp(-vpac(ipp))*(1.d0-vvxp)
     * *(1.d0-vvxpl)-(vpac(ipp)-vpac0(ipp))*(1.d0-exp(-vpac(ipp))
     * *(1.d0-vvxp)*(1.d0-vvxpl)))*exp(-vpac(ipp))*(1.d0-vvxp))
     * *(1.d0-vvx)*(1.d0-vvxt)**2*(1.d0-vvxtl)*exp(-2.d0*vtac(it))
      elseif(itypom.eq.0)then      !'fan' (cut loop at the end - rapgap)
       gb=((max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))
     * +((1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl)
     * +2.d0*vpacng*exp(-2.d0*vpac(ipp))*(1.d0-vvxp)**2*(1.d0-vvxpl))
     * *(vicc*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
     * -vicu*(1.d0-(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))))
     * *(1.d0-vvx)*(1.d0-vvxt)*exp(-vtac(it))
     * -2.d0*vicu*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp)-vtac(it))
     * *(1.d0-vvx)*(1.d0-vvxt)
      elseif(itypom.eq.1)then      !'fan' (uncut end - rapgap)
       gb=((max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))
     * +((1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl)
     * +2.d0*vpacng*exp(-2.d0*vpac(ipp))*(1.d0-vvxp)**2*(1.d0-vvxpl))
     * *(viuc*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
     * +viuu*(1.d0-(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))))
     * *(1.d0-vvx)*(1.d0-vvxt)*exp(-vtac(it))
     * +2.d0*viuu*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp)-vtac(it))
     * *(1.d0-vvx)*(1.d0-vvxt)
      elseif(itypom.eq.2)then        !'fan' (>1 cut Poms at the end)
       gb=vimp*((max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))
     * +((1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl)
     * +2.d0*(vpac0(ipp)*exp(-vpac(ipp))*(1.d0-vvxp)
     * *(1.d0-vvxpl)-(vpac(ipp)-vpac0(ipp))*(1.d0-exp(-vpac(ipp))
     * *(1.d0-vvxp)*(1.d0-vvxpl)))*exp(-vpac(ipp))*(1.d0-vvxp))
     * *(1.d0-vvx)*(1.d0-vvxt)**2*(1.d0-vvxtl)*exp(-2.d0*vtac(it))
      elseif(itypom.eq.3)then      !'fan' (cut/uncut end - rapgap)
       gb=((max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))
     * +((1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl)
     * +2.d0*vpacng*exp(-2.d0*vpac(ipp))*(1.d0-vvxp)**2*(1.d0-vvxpl))
     * *((vicc-wgpm(m)*viuc)*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
     * -(vicu+wgpm(m)*viuu)*(1.d0-(1.d0-vvxt)*(1.d0-vvxtl)
     * *exp(-vtac(it))))*(1.d0-vvx)*(1.d0-vvxt)*exp(-vtac(it))
     * -2.d0*(vicu+wgpm(m)*viuu)*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))
     * -1.d0-(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp)-vtac(it))
     * *(1.d0-vvx)*(1.d0-vvxt)
      elseif(itypom.eq.4)then          !'leg' (single cut Pomeron at the end)
       gb=vip*((vplc0-vplcpe)*exp(-vpac(ipp))*(1.d0-vvxp)
     * *(1.d0-vvxpl)-(vplc-vplc0)*(1.d0-exp(-vpac(ipp))*(1.d0-vvxp)
     * *(1.d0-vvxpl)))*exp(-vpac(ipp)-2.d0*vtac(it))*(1.d0-vvxp)
     * *(1.d0-vvx)*(1.d0-vvxt)**2*(1.d0-vvxtl)
      elseif(itypom.eq.5)then      !'leg' (cut/uncut end - rapgap)
       gb=vplcng*exp(-2.d0*vpac(ipp)-vtac(it))
     * *(1.d0-vvxp)**2*(1.d0-vvxpl)*(1.d0-vvx)*(1.d0-vvxt)
     * *((vicc-wgpm(m)*viuc)*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
     * -(vicu+wgpm(m)*viuu)*(1.d0-(1.d0-vvxt)*(1.d0-vvxtl)
     * *exp(-vtac(it))))
      elseif(itypom.eq.7)then      !'leg' (>1 cut Poms at the end)
       gb=vimp*(vplc0*exp(-vpac(ipp))*(1.d0-vvxp)*(1.d0-vvxpl)
     * -(vplc-vplc0)*(1.d0-exp(-vpac(ipp))*(1.d0-vvxp)*(1.d0-vvxpl)))
     * *exp(-vpac(ipp)-2.d0*vtac(it))*(1.d0-vvxp)
     * *(1.d0-vvx)*(1.d0-vvxt)**2*(1.d0-vvxtl)
      endif
      gb=gb/gb0/z*rp/rp0 /15.d0
      nrej=nrej+1
      if(qgran(b10).gt.gb.and.nrej.lt.1000)goto 15

      if(itypom.eq.-1.or.itypom.eq.4)then  !single cut Pomeron in the handle
       npin=npin+1
       if(npin.gt.npmax)then
        iret=1
        goto 31
       endif
       xpomim(npin)=1.d0/xpomm(m)/scm
       xpomip(npin)=xpomm(m+1)
       bpomim(npin)=bbi
       if(debug.ge.4)write (moniou,211)npin,xpomip(npin),xpomim(npin)
     * ,bpomim(npin)
      elseif(itypom.eq.2.or.itypom.eq.7)then   !>1 cut Pomerons in the handle
       ninc=npgen(vim,2,npmax)
       npin=npin+ninc
       if(npin.gt.npmax)then
        iret=1
        goto 31
       endif
       do i=npin-ninc+1,npin
        xpomim(i)=1.d0/xpomm(m)/scm
        xpomip(i)=xpomm(m+1)
        bpomim(i)=bbi
        if(debug.ge.4)write (moniou,211)i,xpomip(i),xpomim(i)
     *  ,bpomim(i)
       enddo
      endif

      if(itypom.eq.-1)then      !single cut Pomeron in the 'handle'
       vv1=(max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))     
       vv2=(1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl
       vv3=2.d0*((vpac0(ipp)-vpacpe)*exp(-vpac(ipp))*(1.d0-vvxp)
     * *(1.d0-vvxpl)-(vpac(ipp)-vpac0(ipp))*(1.d0-exp(-vpac(ipp))
     * *(1.d0-vvxp)*(1.d0-vvxpl)))*exp(-vpac(ipp))*(1.d0-vvxp)
       if(xpomm(m+1)*sgap**2.gt..9d0.or.vv3.lt.0.d0)vv3=0.d0
       aks=(vv1+vv2+vv3)*qgran(b10)
       if(aks.lt.vv1)then
        jt=1                     !>1 cut fans
       elseif(aks.lt.vv1+vv2)then
        jt=2                     !diffr. cut
       else
        jt=3                     !1 cut fan (without 1P-end)
       endif
      elseif(itypom.eq.0)then    !cut 'loop' in the 'handle' (rap-gap)
       vv1=(max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp)-vtac(it))
     * *(1.d0-vvxt)*(1.d0-vvxtl)*(vicc+vicu)
     * /(vicc*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
     * -vicu*(1.d0-(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))))
       vv1=max(0.d0,vv1)
       vv2=(1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl
       vv3=2.d0*vpacng*exp(-2.d0*vpac(ipp))*(1.d0-vvxp)**2*(1.d0-vvxpl)
       aks=(vv1+vv2+vv3)*qgran(b10)
       if(aks.lt.vv1)then
        jt=1                     !>1 cut fans
       elseif(aks.lt.vv1+vv2)then
        jt=2                     !diffr. cut
       else
        jt=3                     !1 cut fan (no rapgap at the end)
       endif
      elseif(itypom.eq.1)then    !uncut 'handle' (rap-gap)
       vv1=(max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))
       vv2=(1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl
       vv3=2.d0*vpacng*exp(-2.d0*vpac(ipp))*(1.d0-vvxp)**2*(1.d0-vvxpl)
       vv4=2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0-(vpac(ipp)
     * -vpac0(ipp)))*(1.d0-vvxp0)+(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))
     * *exp(-vpac(ipp))*viuu/(viuu*(1.d0-(1.d0-vvxt)*(1.d0-vvxtl)
     * *exp(-vtac(it)))+viuc*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it)))
       if(xpomm(m+1)*sgap**2.gt..9d0.or.vv4.lt.0.d0)vv4=0.d0
       aks=(vv1+vv2+vv3+vv4)*qgran(b10)
       if(aks.lt.vv1)then
        jt=1                     !>1 cut fans
       elseif(aks.lt.vv1+vv2)then
        jt=2                     !diffr. cut
       elseif(aks.lt.vv1+vv2+vv3)then
        jt=3                     !1 cut fan (no rapgap at the end)
       else
        jt=4                     !>1 cut 'handle' fans
       endif
      elseif(itypom.eq.2)then    !>1 cut Pomerons in the 'handle'
       vv1=(max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))*exp(-vpac(ipp))
       vv2=(1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl     
       vv3=2.d0*(vpac0(ipp)*exp(-vpac(ipp))*(1.d0-vvxp)
     * *(1.d0-vvxpl)-(vpac(ipp)-vpac0(ipp))*(1.d0-exp(-vpac(ipp))
     * *(1.d0-vvxp)*(1.d0-vvxpl)))*exp(-vpac(ipp))*(1.d0-vvxp)
       aks=(vv1+vv2+vv3)*qgran(b10)
       if(aks.lt.vv1)then
        jt=1                     !>1 cut fans
       elseif(aks.lt.vv1+vv2)then
        jt=2                     !diffr. cut
       else
        jt=3                     !1 cut fan
       endif
      elseif(itypom.eq.3)then    !rap-gap in the 'handle'
       vv1=(max(0.d0,1.d0-exp(-2.d0*vpac(ipp))
     * *(1.d0+2.d0*vpac(ipp)))+2.d0*vpac(ipp)*exp(-2.d0*vpac(ipp))
     * *(1.d0-(1.d0-vvxp)**2))*(1.d0-vvxpl)
     * -2.d0*(max(0.d0,exp(vpac(ipp)-vpac0(ipp))-1.d0
     * -(vpac(ipp)-vpac0(ipp)))*(1.d0-vvxp0)
     * +(vpac(ipp)-vpac0(ipp))*(vvxp-vvxp0))
     * *exp(-vpac(ipp)-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
     * *(vicc+vicu+wgpm(m)*(viuu-viuc))
     * /((vicc-wgpm(m)*viuc)*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
     * -(vicu+wgpm(m)*viuu)*(1.d0-(1.d0-vvxt)*(1.d0-vvxtl)
     * *exp(-vtac(it))))
       vv1=max(0.d0,vv1)
       vv2=(1.d0-exp(-vpac(ipp)))**2*(1.d0-vvxpl)
     * +2.d0*(1.d0-exp(-vpac(ipp)))*vvxpl
       vv3=2.d0*vpacng*exp(-2.d0*vpac(ipp))*(1.d0-vvxp)**2
     * *(1.d0-vvxpl)
       aks=(vv1+vv2+vv3)*qgran(b10)
       if(aks.lt.vv1)then
        jt=1                     !>1 cut fans
       elseif(aks.lt.vv1+vv2)then
        jt=2                     !diffr. cut
       else
        jt=3                     !1 cut fan (no rapgap at the end)
       endif
      else
       jt=5                      !cut leg
      endif

      nppm(m+1)=0
      wgpm(m+1)=0.d0
      if(jt.eq.1)then                        !>1 cut fans
       ntry=0
16     ntry=ntry+1
       nphm=0
       if(ipp.eq.ia(1).or.ntry.gt.100)then
        nppm(m+1)=npgen(2.d0*vpac(ipp),2,npmax)
        do i=1,nppm(m+1)
         if(qgran(b10).le.vpac0(ipp)/vpac(ipp)
     *   .or.xpomm(m+1)*sgap**2.gt..9d0)then
          itypm(i,m+1)=0
         else
          itypm(i,m+1)=1
          nphm=nphm+1
         endif
         ippm(i,m+1)=ipp
        enddo
        if(nphm.ne.0)wh=(vpac(ipp)/vpac0(ipp)-1.d0)/nphm
       else
        nppm(m+1)=npgen(2.d0*vpac(ipp),1,npmax)
        do i=1,nppm(m+1)
         if(qgran(b10).le.vpac0(ipp)/vpac(ipp)
     *   .or.xpomm(m+1)*sgap**2.gt..9d0)then
          itypm(i,m+1)=0
         else
          itypm(i,m+1)=1
          nphm=nphm+1
         endif
         ippm(i,m+1)=ipp
        enddo
        if(nphm.ne.0)wh=(vpac(ipp)/vpac0(ipp)-1.d0)/nphm
        do ipi=ipp+1,ia(1)
         ninc=npgen(2.d0*vpac(ipi),0,npmax)
         if(ninc.ne.0)then
          nppm(m+1)=nppm(m+1)+ninc
          nh0=nphm
          if(nppm(m+1).gt.legmax)then
           iret=1
           goto 31
          endif
          do i=nppm(m+1)-ninc+1,nppm(m+1)
           if(qgran(b10).le.vpac0(ipi)/vpac(ipi)
     *     .or.xpomm(m+1)*sgap**2.gt..9d0)then
            itypm(i,m+1)=0
           else
            itypm(i,m+1)=1
            nphm=nphm+1
           endif
           ippm(i,m+1)=ipi
          enddo
          if(nphm.gt.nh0)wh=(vpac(ipi)/vpac0(ipi)-1.d0)/(nphm-nh0)
         endif
        enddo
        if(nppm(m+1).eq.1)goto 16
       endif
       
       if(nphm+1.ge.nppm(m+1))then
        if(itypom.eq.-1.or.itypom.eq.1.or.itypom.eq.2)then
         gbt=1.d0-exp(vpac(ipp)+(1.d0-nphm)*dlog(2.d0))
     *   /(1.d0-vvxp)/(1.d0-vvxpl)
        elseif(itypom.eq.0)then
         gbt=1.d0-(vicc+vicu)*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
     *   /(vicc*(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
     *   -vicu*(1.d0-(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))))
     *   *exp(vpac(ipp)+(1.d0-nphm)*dlog(2.d0))
     *   /(1.d0-vvxp)/(1.d0-vvxpl)
        elseif(itypom.eq.3)then
         gbt=1.d0-(vicc+vicu+wgpm(m)*(viuu-viuc))
     *   *(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
     *   /((vicc-wgpm(m)*viuc)*(1.d0-vvxt)*(1.d0-vvxtl)
     *   *exp(-vtac(it))-(vicu+wgpm(m)*viuu)
     *   *(1.d0-(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))))
     *   *exp(vpac(ipp)+(1.d0-nphm)*dlog(2.d0))
     *   /(1.d0-vvxp)/(1.d0-vvxpl)
        else
         stop'unknown itypom'
        endif
        if(nphm.eq.nppm(m+1).and.qgran(b10).gt.gbt
     *  .or.nphm+1.eq.nppm(m+1).and.qgran(b10).gt.1.d0+wh*gbt)then
         ntry=0
          goto 16
        endif
       endif
        
      elseif(jt.eq.4)then                    !>1 cut 'handle' fans
       ntry=0
17     ntry=ntry+1
       if(ipp.eq.ia(1).or.ntry.gt.100)then
        nppm(m+1)=npgen(vpac(ipp)-vpac0(ipp),2,npmax)
        do i=1,nppm(m+1)
         itypm(i,m+1)=1
         ippm(i,m+1)=ipp
        enddo
       else
        nppm(m+1)=npgen(vpac(ipp)-vpac0(ipp),1,npmax)
        do i=1,nppm(m+1)
         itypm(i,m+1)=1
         ippm(i,m+1)=ipp
        enddo
        do ipi=ipp+1,ia(1)
         ninc=npgen(vpac(ipi)-vpac0(ipi),0,npmax)
         if(ninc.ne.0)then
          nppm(m+1)=nppm(m+1)+ninc
          if(nppm(m+1).gt.legmax)then
           iret=1
           goto 31
          endif
          do i=nppm(m+1)-ninc+1,nppm(m+1)
           itypm(i,m+1)=1
           ippm(i,m+1)=ipi
          enddo
         endif
        enddo
        if(nppm(m+1).eq.1)goto 17
       endif
       
      elseif(jt.eq.3)then                    !1 cut fan
       nppm(m+1)=1
       ippm(1,m+1)=ipp
       if(itypom.eq.-1)then             !single cut Pomeron in the 'handle'
        factor=exp(-vpac(ipp))*(1.d0-vvxp)*(1.d0-vvxpl)
        wng=(vpacng-vpacpe)*factor/((vpac0(ipp)-vpacpe)*factor
     *  -(vpac(ipp)-vpac0(ipp))*(1.d0-factor))
        if(qgran(b10).le.wng.or.wng.lt.0.d0)then
         itypm(1,m+1)=2          !>1 cut Pomerons in the 'handle'
        else
         itypm(1,m+1)=3          !rap-gap in the 'handle'
         wgpm(m+1)=(1.d0-factor)/factor*(vpac(ipp)-vpac0(ipp))
     *   /(vpac0(ipp)-vpacng)
        endif
       elseif(itypom.eq.2)then          !>1 cut Pomerons in the 'handle'
        factor=exp(-vpac(ipp))*(1.d0-vvxp)*(1.d0-vvxpl)
        wng=vpacng*factor/(vpac0(ipp)*factor
     *  -(vpac(ipp)-vpac0(ipp))*(1.d0-factor))
        if(qgran(b10).le.wng.or.wng.lt.0.d0
     *  .or.xpomm(m+1)*sgap**2.gt..9d0)then
         if(qgran(b10).le.vpacpe/vpacng
     *   .or.xpomm(m+1)*sgap**2.gt..9d0)then
          itypm(1,m+1)=-1        !single cut Pomeron at the end
         else
          itypm(1,m+1)=2         !>1 cut Pomerons in the 'handle'
         endif
        else
         itypm(1,m+1)=3          !rap-gap in the 'handle'
         wgpm(m+1)=(1.d0-factor)/factor*(vpac(ipp)-vpac0(ipp))
     *   /(vpac0(ipp)-vpacng)
        endif
       else                             !rap-gap in the 'handle' (itypom=0,1,3)
        if(qgran(b10).le.vpacpe/vpacng
     *  .or.xpomm(m+1)*sgap**2.gt..9d0)then
         itypm(1,m+1)=-1         !single cut Pomeron at the end
        else
         itypm(1,m+1)=2          !>1 cut Pomerons in the 'handle'
        endif
       endif

       if(itypm(1,m+1).eq.-1)then     !single cut Pomeron in the 'handle'
        if(qgran(b10).le.vplcp/vpacpe
     *  .or.xpomm(m+1)*sgap**2.gt..9d0)itypm(1,m+1)=6 !single cut Pomeron
       endif

      elseif(jt.eq.5)then                    !cut 'leg'
       nppm(m+1)=1
       ippm(1,m+1)=ipp
       if(itypom.eq.4)then              !single cut Pomeron at the end
        if(xpomm(m+1)*sgap**2.ge.1.d0)stop'=4:xpomm(m+1)*sgap**2>1'
        factor=exp(-vpac(ipp))*(1.d0-vvxp)*(1.d0-vvxpl)
        wng=(vplcng-vplcpe)*factor/((vplc0-vplcpe)*factor
     *  -(vplc-vplc0)*(1.d0-factor))
        if(qgran(b10).le.wng.or.wng.lt.0.d0)then
         itypm(1,m+1)=7          !>1 cut Pomerons in the 'handle'
        else
         itypm(1,m+1)=5          !rap-gap at the end
         wgpm(m+1)=(1.d0-factor)/factor*(vplc-vplc0)/(vplc0-vplcng)
        endif
       elseif(itypom.eq.5)then          !rap-gap at the end (cut or uncut loop)
        if(qgran(b10).le.vplcpe/vplcng
     *  .or.xpomm(m+1)*sgap**2.gt..9d0)then
         itypm(1,m+1)=4          !single cut Pomeron at the end
        else
         itypm(1,m+1)=7          !>1 cut Pomerons in the 'handle'
        endif
       elseif(itypom.eq.7)then          !>1 cut Pomerons in the 'handle'
        factor=exp(-vpac(ipp))*(1.d0-vvxp)*(1.d0-vvxpl)
        wng=vplcng*factor/(vplc0*factor-(vplc-vplc0)*(1.d0-factor))
        if(qgran(b10).le.wng.or.wng.lt.0.d0
     *  .or.xpomm(m+1)*sgap**2.gt..9d0)then
         if(qgran(b10).le.vplcpe/vplcng
     *   .or.xpomm(m+1)*sgap**2.gt..9d0)then
          itypm(1,m+1)=4         !single cut Pomeron at the end
         else
          itypm(1,m+1)=7         !>1 cut Pomerons in the 'handle'
         endif
        else
         itypm(1,m+1)=5          !rap-gap at the end
         wgpm(m+1)=(1.d0-factor)/factor*(vplc-vplc0)/(vplc0-vplcng)
        endif
       endif

       if(itypm(1,m+1).eq.4)then        !single cut Pomeron at the end
        if(qgran(b10).le.vplcp/vplcpe
     *  .or.xpomm(m+1)*sgap**3.gt..9d0)itypm(1,m+1)=6 !single cut Pomeron
       endif
      endif

      if(nppm(m+1).eq.1.and.itypm(1,m+1).eq.6)then  !record single cut Pomeron 
       nppr=nppr+1
       if(nppr.gt.legmax)then
        iret=1
        goto 31
       endif
       if(qgran(b10).le.vplcps/vplcp.or.xpomm(m+1)*s2min.gt..9d0)then
        idpompi(nppr)=0
       else
        idpompi(nppr)=1
       endif
       xpompi(nppr)=xpomm(m+1)
       ipompi(nppr)=ipp
       bpompi(nppr,1)=bbp
       bpompi(nppr,2)=(xb(it,1)-xxm(m+1))**2+(xb(it,2)-yym(m+1))**2
       nppm(m+1)=0
       if(debug.ge.4)write (moniou,209)nppr,ipp,bbp,xpompi(nppr)

      elseif(nppm(m+1).gt.1)then
       vvxpls=vvxpl
       vvxs=(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(it))
       i=0
18     i=i+1
       ityp=itypm(i,m+1)
       if(ityp.eq.0)then    !cut 'handle' fan
        ipi=ippm(i,m+1)
        bbi=(xa(ipi,1)+bcoll-xxm(m+1))**2+(xa(ipi,2)-yym(m+1))**2
        vvxp=0.d0
        vvxpl=0.d0
        vvxp0=0.d0
        if(ia(1).gt.1)then
         do l=1,ia(1)
          if(l.lt.ipi)then
           vvxpl=vvxpl+vpac(l)
          elseif(l.gt.ipi)then
           vvxp=vvxp+vpac(l)
           vvxp0=vvxp0+vpac0(l)
          endif
         enddo
        endif
        vvxp=1.d0-exp(-vvxp)
        vvxpl=1.d0-exp(-vvxpl)
        vvxp0=1.d0-exp(-vvxp0)
        vvxts=1.d0-vvxs*(1.d0-vvxpl)
        vplcps=qgfani(1.d0/xpomm(m+1),bbi,vvxts,vvxp,vvxpl,iddp(ipi)
     *  ,icz,16)
        vplcp=max(vplcps,vpht(ipi,2)
     *  +qgfani(1.d0/xpomm(m+1),bbi,vvxts,vvxp,vvxpl,iddp(ipi),icz,14))
        vpacpe=max(vplcp,vpht(ipi,1)
     *  +qgfani(1.d0/xpomm(m+1),bbi,vvxts,vvxp0,vvxpl,iddp(ipi),icz,8))
        vpacng=max(vpacpe,vpht(ipi,1)
     *  +qgfani(1.d0/xpomm(m+1),bbi,vvxts,vvxp0,vvxpl,iddp(ipi),icz,7))

        aks=qgran(b10)*vpac0(ipi)
        if(aks.le.vplcp.or.xpomm(m+1)*sgap**2.gt..9d0)then
         itypm(i,m+1)=6          !single cut Pomeron
        elseif(aks.lt.vpacpe)then  
         itypm(i,m+1)=-1         !single cut Pomeron in the 'handle'
        elseif(aks.lt.vpacng)then  
         itypm(i,m+1)=2          !>1 cut Pomerons in the 'handle'
        endif

        if(itypm(i,m+1).eq.6)then      !record single cut Pomeron
         nppr=nppr+1
         if(nppr.gt.legmax)then
          iret=1
          goto 31
         endif
         if(qgran(b10).le.vplcps/vplcp.or.xpomm(m+1)*s2min.gt..9d0)then
          idpompi(nppr)=0
         else
          idpompi(nppr)=1
         endif
         xpompi(nppr)=xpomm(m+1)
         ipompi(nppr)=ipi
         bpompi(nppr,1)=bbi
         bpompi(nppr,2)=(xb(it,1)-xxm(m+1))**2+(xb(it,2)-yym(m+1))**2
         if(debug.ge.4)write (moniou,209)nppr,ipi,bbi,xpompi(nppr)
         nppm(m+1)=nppm(m+1)-1
         if(nppm(m+1).ge.i)then
          do l=i,nppm(m+1)
           ippm(l,m+1)=ippm(l+1,m+1)
           itypm(l,m+1)=itypm(l+1,m+1)
          enddo
         endif
         i=i-1
        endif
       endif
       if(i.lt.nppm(m+1))goto 18
       vvxpl=vvxpls
      endif
       
      if(jt.eq.2.and.qgran(b10).gt.2.d0*vvxpl
     */((1.d0-exp(-vpac(ipp)))*(1.d0-vvxpl)+2.d0*vvxpl))then
       if(debug.ge.4)write (moniou,212)
       icdps=iddp(ipp)
       do icdp=1,nfock
        iddp(ipp)=icdp
        call qgfdf(xxm(m+1),yym(m+1),xpomm(m+1),vpac,vtac,vpht,vtht
     *  ,vpacq,vtacq,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht
     *  ,ipp,it)
        wdp(icdp,ipp)=(1.d0-exp(-vpac(ipp)))*(1.d0-vvxpl)
       enddo
       iddp(ipp)=icdps
      endif
       
      if(nppm(m+1).ne.0)then
       goto 13
      else
       goto 14
      endif

20    continue
      if(debug.ge.3)write (moniou,214)nppr      
      if(nptg0.eq.0)goto 31
      
c target 'fans'      
      m=0
      nppm(1)=nptg0
      xpomm(1)=xpomr
      wgpm(1)=wgtg0
      xxm(1)=xxp
      yym(1)=yyp
      do i=1,nptg0
       ippm(i,1)=iptg0(i)
       itypm(i,1)=itytg0(i)
      enddo
                        
21    m=m+1                                   !next level multi-Pomeron vertex
      if(m.gt.levmax)then
       iret=1
       goto 31
      endif
      ii(m)=0
22    ii(m)=ii(m)+1                           !next cut fan in the vertex
      if(ii(m).gt.nppm(m))then                !all fans at the level considered
       m=m-1                                  !one level down
       if(m.eq.0)goto 31                      !all targ. fans considered
       goto 22
      endif 
      l=ii(m)
      itt=ippm(l,m)                           !targ. index for the leg      
      itypom=itypm(l,m)                       !type of the cut      
      btm=(xb(itt,1)-xxm(m))**2+(xb(itt,2)-yym(m))**2  !b^2 for the leg  
      if(debug.ge.4)write (moniou,216)ii(m),m,itt,btm
      if(xpomm(m)*scm.lt.sgap**2)stop'xpomm(m)*scm<sgap**2!'
      if(itypom.eq.4.and.xpomm(m)*scm.lt.sgap**3)
     *stop'4t:xpomm(m)*scm<sgap**3!'
       
      if(debug.ge.4)write (moniou,210)m
      xpomr0=dsqrt(xpomm(m)/scm) 
      if(itypom.eq.-1)xpomr0=min(xpomr0,1.1d0*sgap/scm)
      if(itypom.eq.4)xpomr0=dsqrt(xpomm(m)*sgap/scm)
      rp1=(rq(iddt(itt),2)+alfp*dlog(xpomr0*scm))*4.d0*.0389d0  
      rp2=alfp*dlog(xpomm(m)/xpomr0)*4.d0*.0389d0 
      rp0=rp1*rp2/(rp1+rp2)
      bbt=btm*(rp1/(rp1+rp2))**2 
      bbi=btm*(rp2/(rp1+rp2))**2  
      call qgbdef(bbt,bbi,xb(itt,1),xb(itt,2),xxm(m),yym(m)
     *,xxp0,yyp0,1)     
       
      call qgfdf(xxp0,yyp0,xpomr0,vpac,vtac,vpht,vtht,vpacq,vtacq
     *,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,itt)
      vvxps=1.d0-(1.d0-vvx)*(1.d0-vvxp)*exp(-vpac(ip))
      viu=qgpini(xpomm(m)/xpomr0,bbi,0.d0,0.d0,2)
      vim=2.d0*min(viu,qgpini(xpomm(m)/xpomr0,bbi,0.d0,0.d0,14))
      if(itypom.eq.-1.or.itypom.eq.4)then      !single cut Pomeron at the end
       vvxi=1.d0-(1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt)
     * *exp(-vpac(ip)-vtac(itt))
       vip=qgpini(xpomm(m)/xpomr0,bbi,vvxi,0.d0,19)
     * +qgpini(xpomm(m)/xpomr0,bbi,0.d0,0.d0,8)*(exp(-vim)-1.d0)
       if(vip.lt.0.d0)vip=qgpini(xpomm(m)/xpomr0,bbi,0.d0,0.d0,8)
     * *exp(-vim)
      elseif(itypom.eq.2.or.itypom.eq.7)then   !>1 cut Pomerons at the end
       vimp=max(0.d0,1.d0-exp(-vim)*(1.d0+vim))
      else                                     !rap-gap at the end
       vvxpin=1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
       vvxtin=1.d0-(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(itt))
       viuu=qgpini(xpomm(m)/xpomr0,bbi,vvxtin,vvxpin,20)
       viuc=max(0.d0,viuu-qgpini(xpomm(m)/xpomr0,bbi,vvxtin,vvxpin,21))
       viuu=max(0.d0,viuu+min(0.d0,1.d0-exp(-viu)-viu))
       vicc=qgpini(xpomm(m)/xpomr0,bbi,vvxtin,vvxpin,22)
       vicu=max(0.d0,qgpini(xpomm(m)/xpomr0,bbi,vvxtin,vvxpin,23)-vicc)
       vicc=max(0.d0,vicc+.5d0*(1.d0-exp(-viu))**2-viu
     * +.5d0*(exp(2.d0*viu-vim)-1.d0)*exp(-2.d0*viu))
      endif

      vtlcps=qgfani(xpomr0*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,16)
      vtlcp=max(vtlcps,vtht(itt,2)
     *+qgfani(xpomr0*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,14))
      if(itypom.le.3)then                         !cut 'fan'
       sumut=0.d0
       vvxt0=0.d0
       do i=1,ia(2)
        sumut=sumut+vtac(i)
       enddo
       vvxs=(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
       do i=1,ia(2)-itt+1
        iti=ia(2)-i+1
        bbl=(xb(iti,1)-xxp0)**2+(xb(iti,2)-yyp0)**2
        sumut=sumut-vtac(iti)
        vts=qgfani(xpomr0*scm,bbl,1.d0-vvxs*exp(-sumut)
     *  ,1.d0-exp(-vvxt0),1.d0-exp(-sumut),iddt(iti),2,3)
        vtac0(iti)=max(vts,vtht(iti,1)
     *  +qgfani(xpomr0*scm,bbl,1.d0-vvxs*exp(-sumut)
     *  ,1.d0-exp(-vvxt0),1.d0-exp(-sumut),iddt(iti),2,6))
        vtac0(iti)=min(vtac(iti),vtac0(iti))
        if(iti.gt.itt)vvxt0=vvxt0+vtac0(iti)
       enddo
       vvxt0=1.d0-exp(-vvxt0)
       vtacpe=max(vtlcp,vtht(itt,1)
     * +qgfani(xpomr0*scm,bbt,vvxps,vvxt0,vvxtl,iddt(itt),2,8))
       vtacng=max(vtacpe,vtht(itt,1)
     * +qgfani(xpomr0*scm,bbt,vvxps,vvxt0,vvxtl,iddt(itt),2,7))    
      else                                        !cut 'leg'
       vtlcpe=max(vtlcp,vtht(itt,1)
     * +qgfani(xpomr0*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,13))
       vtlcng=max(vtlcpe,vtht(itt,1)
     * +qgfani(xpomr0*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,12))
       vtlc0=max(vtlcng,vtht(itt,1)
     * +qgfani(xpomr0*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,11))
       vtlc=max(vtlc0,vtht(itt,1)
     * +qgfani(xpomr0*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,10))
      endif

      if(itypom.eq.-1)then         !'fan' (single cut Pomeron at the end)
       gb0=vip*((max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))
     * +((1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl)
     * +2.d0*((vtac0(itt)-vtacpe)*exp(-vtac(itt))*(1.d0-vvxt)
     * *(1.d0-vvxtl)-(vtac(itt)-vtac0(itt))*(1.d0-exp(-vtac(itt))
     * *(1.d0-vvxt)*(1.d0-vvxtl)))*exp(-vtac(itt))*(1.d0-vvxt))
     * *(1.d0-vvx)*(1.d0-vvxp)**2*(1.d0-vvxpl)*exp(-2.d0*vpac(ip))
       gb0=gb0*50.d0
      elseif(itypom.eq.0)then      !'fan' (cut loop at the end - rapgap)
       gb0=((max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))
     * +((1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl)
     * +2.d0*vtacng*exp(-2.d0*vtac(itt))*(1.d0-vvxt)**2*(1.d0-vvxtl))
     * *(vicc*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
     * -vicu*(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))))
     * *(1.d0-vvx)*(1.d0-vvxp)*exp(-vpac(ip))
     * -2.d0*vicu*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt)-vpac(ip))
     * *(1.d0-vvx)*(1.d0-vvxp)
      elseif(itypom.eq.1)then      !'fan' (uncut end - rapgap)
       gb0=((max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))
     * +((1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl)
     * +2.d0*vtacng*exp(-2.d0*vtac(itt))*(1.d0-vvxt)**2*(1.d0-vvxtl))
     * *(viuc*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
     * +viuu*(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))))
     * *(1.d0-vvx)*(1.d0-vvxp)*exp(-vpac(ip))
     * +2.d0*viuu*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt)-vpac(ip))
     * *(1.d0-vvx)*(1.d0-vvxp)
      elseif(itypom.eq.2)then      !'fan' (>1 cut Poms at the end)
       gb0=vimp*((max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))
     * +((1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl)
     * +2.d0*(vtac0(itt)*exp(-vtac(itt))*(1.d0-vvxt)
     * *(1.d0-vvxtl)-(vtac(itt)-vtac0(itt))*(1.d0-exp(-vtac(itt))
     * *(1.d0-vvxt)*(1.d0-vvxtl)))*exp(-vtac(itt))*(1.d0-vvxt))
     * *(1.d0-vvx)*(1.d0-vvxp)**2*(1.d0-vvxpl)*exp(-2.d0*vpac(ip))
      elseif(itypom.eq.3)then      !'fan' (cut/uncut end - rapgap)
       gb0=((max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))
     * +((1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl)
     * +2.d0*vtacng*exp(-2.d0*vtac(itt))*(1.d0-vvxt)**2*(1.d0-vvxtl))
     * *(vicc-wgpm(m)*viuc)*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
     * *(1.d0-vvx)*(1.d0-vvxp)*exp(-vpac(ip))
      elseif(itypom.eq.4)then      !'leg' (single cut Pomeron at the end)
       gb0=vip*((vtlc0-vtlcpe)*exp(-vtac(itt))*(1.d0-vvxt)
     * *(1.d0-vvxtl))*exp(-vtac(itt)-2.d0*vpac(ip))*(1.d0-vvxt)
     * *(1.d0-vvx)*(1.d0-vvxp)**2*(1.d0-vvxpl)
      elseif(itypom.eq.5)then      !'leg' (cut/uncut end - rapgap)
       gb0=vtlcng*exp(-2.d0*vtac(itt)-vpac(ip))
     * *(1.d0-vvxt)**2*(1.d0-vvxtl)*(1.d0-vvx)*(1.d0-vvxp)
     * *(vicc-wgpm(m)*viuc)*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
      elseif(itypom.eq.7)then      !'leg' (>1 cut Poms at the end)
       gb0=vimp*(vtlc0*exp(-vtac(itt))*(1.d0-vvxt)*(1.d0-vvxtl)
     * -(vtlc-vtlc0)*(1.d0-exp(-vtac(itt))*(1.d0-vvxt)*(1.d0-vvxtl)))
     * *exp(-vtac(itt)-2.d0*vpac(ip))*(1.d0-vvxt)
     * *(1.d0-vvx)*(1.d0-vvxp)**2*(1.d0-vvxpl)
      endif
      if(gb0.le.0.d0)then
       nptg=nptg+1
       if(nptg.gt.legmax)then
        iret=1
        goto 31
       endif
       if(qgran(b10).le.vtlcps/vtlcp)then
        idpomti(nptg)=0
       else
        idpomti(nptg)=1
       endif
       xpomti(nptg)=xpomm(m)
       ipomti(nptg)=itt
       bpomti(nptg,1)=btm
       bpomti(nptg,2)=(xa(ip,1)+bcoll-xxm(m))**2+(xa(ip,2)-yym(m))**2
       goto 22
      endif
      nrej=0

23    xpomm(m+1)=xpomm(m)/sgap/(xpomm(m)*scm/sgap**2)**qgran(b10)
      if(itypom.eq.4)xpomm(m+1)=xpomm(m)/sgap
     */(xpomm(m)*scm/sgap**3)**qgran(b10)     
      rp1=(rq(iddt(itt),2)+alfp*dlog(xpomm(m+1)*scm))*4.d0*.0389d0
      rp2=alfp*dlog(xpomm(m)/xpomm(m+1))*4.d0*.0389d0 
      rp=rp1*rp2/(rp1+rp2)
      z=qgran(b10)
      phi=pi*qgran(b10)
      b0=dsqrt(-rp*dlog(z))
      bbt=(dsqrt(btm)*rp1/(rp1+rp2)+b0*cos(phi))**2+(b0*sin(phi))**2 
      bbi=(dsqrt(btm)*rp2/(rp1+rp2)-b0*cos(phi))**2+(b0*sin(phi))**2 
      call qgbdef(bbt,bbi,xb(itt,1),xb(itt,2),xxm(m),yym(m)
     *,xxm(m+1),yym(m+1),int(1.5d0+qgran(b10)))   !coordinates for the vertex

      call qgfdf(xxm(m+1),yym(m+1),xpomm(m+1),vpac,vtac,vpht,vtht
     *,vpacq,vtacq,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,itt)
      vvxps=1.d0-(1.d0-vvx)*(1.d0-vvxp)*exp(-vpac(ip))
      viu=qgpini(xpomm(m)/xpomm(m+1),bbi,0.d0,0.d0,2)
      vim=2.d0*min(viu,qgpini(xpomm(m)/xpomm(m+1),bbi,0.d0,0.d0,14))
      if(itypom.eq.-1.or.itypom.eq.4)then      !single cut Pomeron at the end
       vvxi=1.d0-(1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt)
     * *exp(-vpac(ip)-vtac(itt))
       vip=qgpini(xpomm(m)/xpomm(m+1),bbi,vvxi,0.d0,19)
     * +qgpini(xpomm(m)/xpomm(m+1),bbi,0.d0,0.d0,8)*(exp(-vim)-1.d0)
       if(vip.lt.0.d0)vip=qgpini(xpomm(m)/xpomm(m+1),bbi,0.d0,0.d0,8)
     * *exp(-vim)
      elseif(itypom.eq.2.or.itypom.eq.7)then   !>1 cut Pomerons at the end
       vimp=max(0.d0,1.d0-exp(-vim)*(1.d0+vim))
      else                                     !rap-gap at the end
       vvxpin=1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
       vvxtin=1.d0-(1.d0-vvxt)*(1.d0-vvxtl)*exp(-vtac(itt))
       viuu=qgpini(xpomm(m)/xpomm(m+1),bbi,vvxtin,vvxpin,20)
       viuc=max(0.d0,viuu-qgpini(xpomm(m)/xpomm(m+1),bbi,vvxtin,vvxpin
     * ,21))
       viuu=max(0.d0,viuu+min(0.d0,1.d0-exp(-viu)-viu))
       vicc=qgpini(xpomm(m)/xpomm(m+1),bbi,vvxtin,vvxpin,22)
       vicu=max(0.d0,qgpini(xpomm(m)/xpomm(m+1),bbi,vvxtin,vvxpin,23)
     * -vicc)
       vicc=max(0.d0,vicc+.5d0*(1.d0-exp(-viu))**2-viu
     * +.5d0*(exp(2.d0*viu-vim)-1.d0)*exp(-2.d0*viu))
      endif

      vtlcps=qgfani(xpomm(m+1)*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,16)
      vtlcp=max(vtlcps,vtht(itt,2)
     *+qgfani(xpomm(m+1)*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,14))
      if(itypom.le.3)then                         !cut 'fan'
       sumut=0.d0
       vvxt0=0.d0
       do i=1,ia(2)
        sumut=sumut+vtac(i)
       enddo
       vvxs=(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
       do i=1,ia(2)-itt+1
        iti=ia(2)-i+1
        bbl=(xb(iti,1)-xxm(m+1))**2+(xb(iti,2)-yym(m+1))**2
        sumut=sumut-vtac(iti)
        vts=qgfani(xpomm(m+1)*scm,bbl,1.d0-vvxs*exp(-sumut)
     *  ,1.d0-exp(-vvxt0),1.d0-exp(-sumut),iddt(iti),2,3)
        vtac0(iti)=max(vts,vtht(iti,1)
     *  +qgfani(xpomm(m+1)*scm,bbl,1.d0-vvxs*exp(-sumut)
     *  ,1.d0-exp(-vvxt0),1.d0-exp(-sumut),iddt(iti),2,6))
        vtac0(iti)=min(vtac(iti),vtac0(iti))
        if(iti.gt.itt)vvxt0=vvxt0+vtac0(iti)
       enddo
       vvxt0=1.d0-exp(-vvxt0)
       vtacpe=max(vtlcp,vtht(itt,1)
     * +qgfani(xpomm(m+1)*scm,bbt,vvxps,vvxt0,vvxtl,iddt(itt),2,8))
       vtacng=max(vtacpe,vtht(itt,1)
     * +qgfani(xpomm(m+1)*scm,bbt,vvxps,vvxt0,vvxtl,iddt(itt),2,7))
      else                                        !cut 'leg'
       vtlcpe=max(vtlcp,vtht(itt,1)
     * +qgfani(xpomm(m+1)*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,13))
       vtlcng=max(vtlcpe,vtht(itt,1)
     * +qgfani(xpomm(m+1)*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,12))
       vtlc0=max(vtlcng,vtht(itt,1)
     * +qgfani(xpomm(m+1)*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,11))
       vtlc=max(vtlc0,vtht(itt,1)
     * +qgfani(xpomm(m+1)*scm,bbt,vvxps,vvxt,vvxtl,iddt(itt),2,10))
      endif

      if(itypom.eq.-1)then         !'fan' (single cut Pomeron at the end)
       gb=vip*((max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))
     * +((1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl)
     * +2.d0*((vtac0(itt)-vtacpe)*exp(-vtac(itt))*(1.d0-vvxt)
     * *(1.d0-vvxtl)-(vtac(itt)-vtac0(itt))*(1.d0-exp(-vtac(itt))
     * *(1.d0-vvxt)*(1.d0-vvxtl)))*exp(-vtac(itt))*(1.d0-vvxt))
     * *(1.d0-vvx)*(1.d0-vvxp)**2*(1.d0-vvxpl)*exp(-2.d0*vpac(ip))
      elseif(itypom.eq.0)then      !'fan' (cut loop at the end - rapgap)
       gb=((max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))
     * +((1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl)
     * +2.d0*vtacng*exp(-2.d0*vtac(itt))*(1.d0-vvxt)**2*(1.d0-vvxtl))
     * *(vicc*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
     * -vicu*(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))))
     * *(1.d0-vvx)*(1.d0-vvxp)*exp(-vpac(ip))
     * -2.d0*vicu*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt)-vpac(ip))
     * *(1.d0-vvx)*(1.d0-vvxp)
      elseif(itypom.eq.1)then      !'fan' (uncut end - rapgap)
       gb=((max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))
     * +((1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl)
     * +2.d0*vtacng*exp(-2.d0*vtac(itt))*(1.d0-vvxt)**2*(1.d0-vvxtl))
     * *(viuc*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
     * +viuu*(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))))
     * *(1.d0-vvx)*(1.d0-vvxp)*exp(-vpac(ip))
     * +2.d0*viuu*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt)-vpac(ip))
     * *(1.d0-vvx)*(1.d0-vvxp)
      elseif(itypom.eq.2)then      !'fan' (>1 cut Poms at the end)
       gb=vimp*((max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))
     * +((1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl)
     * +2.d0*(vtac0(itt)*exp(-vtac(itt))*(1.d0-vvxt)
     * *(1.d0-vvxtl)-(vtac(itt)-vtac0(itt))*(1.d0-exp(-vtac(itt))
     * *(1.d0-vvxt)*(1.d0-vvxtl)))*exp(-vtac(itt))*(1.d0-vvxt))
     * *(1.d0-vvx)*(1.d0-vvxp)**2*(1.d0-vvxpl)*exp(-2.d0*vpac(ip))
      elseif(itypom.eq.3)then      !'fan' (cut/uncut end - rapgap)
       gb=((max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))
     * +((1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl)
     * +2.d0*vtacng*exp(-2.d0*vtac(itt))*(1.d0-vvxt)**2*(1.d0-vvxtl))
     * *((vicc-wgpm(m)*viuc)*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
     * -(vicu+wgpm(m)*viuu)*(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)
     * *exp(-vpac(ip))))*(1.d0-vvx)*(1.d0-vvxp)*exp(-vpac(ip))
     * -2.d0*(vicu+wgpm(m)*viuu)*(max(0.d0,exp(vtac(itt)-vtac0(itt))
     * -1.d0-(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt)-vpac(ip))
     * *(1.d0-vvx)*(1.d0-vvxp)
      elseif(itypom.eq.4)then      !'leg' (single cut Pomeron at the end)
       gb=vip*((vtlc0-vtlcpe)*exp(-vtac(itt))*(1.d0-vvxt)
     * *(1.d0-vvxtl)-(vtlc-vtlc0)*(1.d0-exp(-vtac(itt))*(1.d0-vvxt)
     * *(1.d0-vvxtl)))*exp(-vtac(itt)-2.d0*vpac(ip))*(1.d0-vvxt)
     * *(1.d0-vvx)*(1.d0-vvxp)**2*(1.d0-vvxpl)
      elseif(itypom.eq.5)then      !'leg' (cut/uncut end - rapgap)
       gb=vtlcng*exp(-2.d0*vtac(itt)-vpac(ip))
     * *(1.d0-vvxt)**2*(1.d0-vvxtl)*(1.d0-vvx)*(1.d0-vvxp)
     * *((vicc-wgpm(m)*viuc)*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
     * -(vicu+wgpm(m)*viuu)*(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)
     * *exp(-vpac(ip))))    
      elseif(itypom.eq.7)then      !'leg' (>1 cut Poms at the end)
       gb=vimp*(vtlc0*exp(-vtac(itt))*(1.d0-vvxt)*(1.d0-vvxtl)
     * -(vtlc-vtlc0)*(1.d0-exp(-vtac(itt))*(1.d0-vvxt)*(1.d0-vvxtl)))
     * *exp(-vtac(itt)-2.d0*vpac(ip))*(1.d0-vvxt)
     * *(1.d0-vvx)*(1.d0-vvxp)**2*(1.d0-vvxpl)
      endif
      nrej=nrej+1
      gb=gb/gb0/z*rp/rp0 /15.d0
      if(qgran(b10).gt.gb.and.nrej.lt.1000)goto 23

      if(itypom.eq.-1.or.itypom.eq.4)then    !'single cut Pomeron in the handle
       npin=npin+1
       if(npin.gt.npmax)then
        iret=1
        goto 31
       endif
       xpomim(npin)=1.d0/xpomm(m+1)/scm
       xpomip(npin)=xpomm(m)
       bpomim(npin)=bbi
       if(debug.ge.4)write (moniou,211)npin,xpomip(npin),xpomim(npin)
     * ,bpomim(npin)
      elseif(itypom.eq.2.or.itypom.eq.7)then !>1 cut Pomerons in the handle
       ninc=npgen(vim,2,npmax)
       npin=npin+ninc
       if(npin.gt.npmax)then
        iret=1
        goto 31
       endif
       do i=npin-ninc+1,npin
        xpomim(i)=1.d0/xpomm(m+1)/scm
        xpomip(i)=xpomm(m)
        bpomim(i)=bbi
        if(debug.ge.4)write (moniou,211)i,xpomip(i),xpomim(i),bpomim(i)
       enddo
      endif
        
      if(itypom.eq.-1)then      !single cut Pomeron in the 'handle'
       vv1=(max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))    
       vv1=max(0.d0,vv1)
       vv2=(1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl    
       vv3=2.d0*((vtac0(itt)-vtacpe)*exp(-vtac(itt))*(1.d0-vvxt)
     * *(1.d0-vvxtl)-(vtac(itt)-vtac0(itt))*(1.d0-exp(-vtac(itt))
     * *(1.d0-vvxt)*(1.d0-vvxtl)))*exp(-vtac(itt))*(1.d0-vvxt)
       if(xpomm(m+1)*scm.lt.1.1d0*sgap**2.or.vv3.lt.0.d0)vv3=0.d0
       aks=(vv1+vv2+vv3)*qgran(b10)
       if(aks.lt.vv1)then
        jt=1                     !>1 cut fans
       elseif(aks.lt.vv1+vv2)then
        jt=2                     !diffr. cut
       else
        jt=3                     !1 cut fan (without 1P-end)
       endif
      elseif(itypom.eq.0)then      !cut 'loop' in the 'handle'
       vv1=(max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt)-vpac(ip))
     * *(1.d0-vvxp)*(1.d0-vvxpl)*(vicc+vicu)
     * /(vicc*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
     * -vicu*(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))))
       vv1=max(0.d0,vv1)
       vv2=(1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl
       vv3=2.d0*vtacng*exp(-2.d0*vtac(itt))*(1.d0-vvxt)**2*(1.d0-vvxtl)
       aks=(vv1+vv2+vv3)*qgran(b10)
       if(aks.lt.vv1)then
        jt=1                     !>1 cut fans
       elseif(aks.lt.vv1+vv2)then
        jt=2                     !diffr. cut
       else
        jt=3                     !1 cut fan (no rapgap at the end)
       endif
      elseif(itypom.eq.1)then    !uncut 'handle' (rap-gap)
       vv1=(max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))
       vv2=(1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl
       vv3=2.d0*vtacng*exp(-2.d0*vtac(itt))*(1.d0-vvxt)**2*(1.d0-vvxtl)
       vv4=2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0-(vtac(itt)
     * -vtac0(itt)))*(1.d0-vvxt0)+(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))
     * *exp(-vtac(itt))*viuu/(viuu*(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)
     * *exp(-vpac(ip)))+viuc*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip)))
       if(xpomm(m+1)*scm.lt.1.1d0*sgap**2.or.vv4.lt.0.d0)vv4=0.d0
       aks=(vv1+vv2+vv3+vv4)*qgran(b10)
       if(aks.lt.vv1)then
        jt=1                     !>1 cut fans
       elseif(aks.lt.vv1+vv2)then
        jt=2                     !diffr. cut
       elseif(aks.lt.vv1+vv2+vv3)then
        jt=3                     !1 cut fan (no rapgap at the end)
       else
        jt=4                     !>1 cut 'handle' fans
       endif
      elseif(itypom.eq.2)then    !>1 cut Pomerons in the 'handle'
       vv1=(max(0.d0,1.d0-exp(-2.d0*vtac(itt))
     * *(1.d0+2.d0*vtac(itt)))+2.d0*vtac(itt)*exp(-2.d0*vtac(itt))
     * *(1.d0-(1.d0-vvxt)**2))*(1.d0-vvxtl)
     * -2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)
     * +(vtac(itt)-vtac0(itt))*(vvxt-vvxt0))*exp(-vtac(itt))
       vv2=(1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl
       vv3=2.d0*(vtac0(itt)*exp(-vtac(itt))*(1.d0-vvxt)
     * *(1.d0-vvxtl)-(vtac(itt)-vtac0(itt))*(1.d0-exp(-vtac(itt))
     * *(1.d0-vvxt)*(1.d0-vvxtl)))*exp(-vtac(itt))*(1.d0-vvxt)
       aks=(vv1+vv2+vv3)*qgran(b10)
       if(aks.lt.vv1)then
        jt=1                     !>1 cut fans
       elseif(aks.lt.vv1+vv2)then
        jt=2                     !diffr. cut
       else
        jt=3                     !1 cut fan
       endif
      elseif(itypom.eq.3)then    !rap-gap in the 'handle'
       vv1=(max(0.d0,1.d0-exp(-2.d0*vtac(itt))*(1.d0+2.d0*vtac(itt)))
     * +2.d0*vtac(itt)*exp(-2.d0*vtac(itt))*(1.d0-(1.d0-vvxt)**2))
     * *(1.d0-vvxtl)-2.d0*(max(0.d0,exp(vtac(itt)-vtac0(itt))-1.d0
     * -(vtac(itt)-vtac0(itt)))*(1.d0-vvxt0)+(vtac(itt)-vtac0(itt))
     * *(vvxt-vvxt0))*exp(-vtac(itt)-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
     * *(vicc+vicu+wgpm(m)*(viuu-viuc))
     * /((vicc-wgpm(m)*viuc)*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
     * -(vicu+wgpm(m)*viuu)*(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)
     * *exp(-vpac(ip))))
       vv1=max(0.d0,vv1)
       vv2=(1.d0-exp(-vtac(itt)))**2*(1.d0-vvxtl)
     * +2.d0*(1.d0-exp(-vtac(itt)))*vvxtl
       vv3=2.d0*vtacng*exp(-2.d0*vtac(itt))*(1.d0-vvxt)**2
     * *(1.d0-vvxtl)
       aks=(vv1+vv2+vv3)*qgran(b10)
       if(aks.lt.vv1)then
        jt=1                     !>1 cut fans
       elseif(aks.lt.vv1+vv2)then
        jt=2                     !diffr. cut
       else
        jt=3                     !1 cut fan (no rapgap at the end)
       endif
      else
       jt=5                      !cut leg
      endif

      nppm(m+1)=0
      wgpm(m+1)=0.d0
      if(jt.eq.1)then                        !>1 cut fans
       ntry=0
24     ntry=ntry+1
       nphm=0
       if(itt.eq.ia(2).or.ntry.gt.100)then
        nppm(m+1)=npgen(2.d0*vtac(itt),2,npmax)
        do i=1,nppm(m+1)
         if(qgran(b10).le.vtac0(itt)/vtac(itt)
     *   .or.xpomm(m+1)*scm.lt.1.1d0*sgap**2)then
          itypm(i,m+1)=0
         else
          nphm=nphm+1
          itypm(i,m+1)=1
         endif
         ippm(i,m+1)=itt
        enddo
        if(nphm.ne.0)wh=(vtac(itt)/vtac0(itt)-1.d0)/nphm
       else
        nppm(m+1)=npgen(2.d0*vtac(itt),1,npmax)
        do i=1,nppm(m+1)
         if(qgran(b10).le.vtac0(itt)/vtac(itt)
     *   .or.xpomm(m+1)*scm.lt.1.1d0*sgap**2)then
          itypm(i,m+1)=0
         else
          nphm=nphm+1
          itypm(i,m+1)=1
         endif
         ippm(i,m+1)=itt
        enddo
        if(nphm.ne.0)wh=(vtac(itt)/vtac0(itt)-1.d0)/nphm
        do iti=itt+1,ia(2)
         ninc=npgen(2.d0*vtac(iti),0,npmax)
         if(ninc.ne.0)then
          nppm(m+1)=nppm(m+1)+ninc
          nh0=nphm
          if(nppm(m+1).gt.legmax)then
           iret=1
           goto 31
          endif
          do i=nppm(m+1)-ninc+1,nppm(m+1)
           if(qgran(b10).le.vtac0(iti)/vtac(iti)
     *     .or.xpomm(m+1)*scm.lt.1.1d0*sgap**2)then
            itypm(i,m+1)=0
           else
            nphm=nphm+1
            itypm(i,m+1)=1
           endif
           ippm(i,m+1)=iti
          enddo
          if(nphm.gt.nh0)wh=(vtac(iti)/vtac0(iti)-1.d0)/(nphm-nh0)
         endif
        enddo
        if(nppm(m+1).eq.1)goto 24
       endif
       
       if(nphm+1.ge.nppm(m+1))then
        if(itypom.eq.-1.or.itypom.eq.1.or.itypom.eq.2)then
         gbt=1.d0-exp(vtac(itt)+(1.d0-nphm)*dlog(2.d0))
     *   /(1.d0-vvxt)/(1.d0-vvxtl)
        elseif(itypom.eq.0)then
         gbt=1.d0-(vicc+vicu)*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
     *   /(vicc*(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
     *   -vicu*(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))))
     *   *exp(vtac(itt)+(1.d0-nphm)*dlog(2.d0))
     *   /(1.d0-vvxt)/(1.d0-vvxtl)
        elseif(itypom.eq.3)then
         gbt=1.d0-(vicc+vicu+wgpm(m)*(viuu-viuc))
     *   *(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
     *   /((vicc-wgpm(m)*viuc)*(1.d0-vvxp)*(1.d0-vvxpl)
     *   *exp(-vpac(ip))-(vicu+wgpm(m)*viuu)
     *   *(1.d0-(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))))
     *   *exp(vtac(itt)+(1.d0-nphm)*dlog(2.d0))
     *   /(1.d0-vvxt)/(1.d0-vvxtl)
        else
         stop'unknown itypom'
        endif
        if(nphm.eq.nppm(m+1).and.qgran(b10).gt.gbt
     *  .or.nphm+1.eq.nppm(m+1).and.qgran(b10).gt.1.d0+wh*gbt)then
         ntry=0
          goto 24
        endif
       endif
        
      elseif(jt.eq.4)then                    !>1 cut 'handle' fans
       ntry=0
25     ntry=ntry+1
       if(itt.eq.ia(2).or.ntry.gt.100)then
        nppm(m+1)=npgen(vtac(itt)-vtac0(itt),2,npmax)
        do i=1,nppm(m+1)
         itypm(i,m+1)=1
         ippm(i,m+1)=itt
        enddo
       else
        nppm(m+1)=npgen(vtac(itt)-vtac0(itt),1,npmax)
        do i=1,nppm(m+1)
         itypm(i,m+1)=1
         ippm(i,m+1)=itt
        enddo
        do iti=itt+1,ia(2)
         ninc=npgen(vtac(iti)-vtac0(iti),0,npmax)
         if(ninc.ne.0)then
          nppm(m+1)=nppm(m+1)+ninc
          if(nppm(m+1).gt.legmax)then
           iret=1
           goto 31
          endif
          do i=nppm(m+1)-ninc+1,nppm(m+1)
           itypm(i,m+1)=1
           ippm(i,m+1)=iti
          enddo
         endif
        enddo
        if(nppm(m+1).eq.1)goto 25
       endif
       
      elseif(jt.eq.3)then                    !1 cut fan
       nppm(m+1)=1
       ippm(1,m+1)=itt
       if(itypom.eq.-1)then             !single cut Pomeron in the 'handle'
        factor=exp(-vtac(itt))*(1.d0-vvxt)*(1.d0-vvxtl)
        wng=(vtacng-vtacpe)*factor/((vtac0(itt)-vtacpe)*factor
     *  -(vtac(itt)-vtac0(itt))*(1.d0-factor))
        if(qgran(b10).le.wng.or.wng.lt.0.d0)then
         itypm(1,m+1)=2          !>1 cut Pomerons in the 'handle'
        else
         itypm(1,m+1)=3          !rap-gap in the 'handle'
         wgpm(m+1)=(1.d0-factor)/factor*(vtac(itt)-vtac0(itt))
     *   /(vtac0(itt)-vtacng)
        endif
       elseif(itypom.eq.2)then          !>1 cut Pomerons in the 'handle'
        factor=exp(-vtac(itt))*(1.d0-vvxt)*(1.d0-vvxtl)
        wng=vtacng*factor/(vtac0(itt)*factor
     *  -(vtac(itt)-vtac0(itt))*(1.d0-factor))
        if(qgran(b10).le.wng.or.wng.lt.0.d0
     *  .or.xpomm(m+1)*scm.lt.1.1d0*sgap**2)then
         if(qgran(b10).le.vtacpe/vtacng
     *   .or.xpomm(m+1)*scm.lt.1.1d0*sgap**2)then
          itypm(1,m+1)=-1        !single cut Pomeron at the end
         else
          itypm(1,m+1)=2         !>1 cut Pomerons in the 'handle'
         endif
        else
         itypm(1,m+1)=3          !rap-gap in the 'handle'
         wgpm(m+1)=(1.d0-factor)/factor*(vtac(itt)-vtac0(itt))
     *   /(vtac0(itt)-vtacng)
        endif
       else                             !rap-gap in the 'handle'
        if(qgran(b10).le.vtacpe/vtacng
     *  .or.xpomm(m+1)*scm.lt.1.1d0*sgap**2)then
         itypm(1,m+1)=-1         !single cut Pomeron at the end
        else
         itypm(1,m+1)=2          !>1 cut Pomerons in the 'handle'
        endif
       endif
        
       if(itypm(1,m+1).eq.-1)then     !single cut Pomeron at the end
        if(qgran(b10).le.vtlcp/vtacpe
     *  .or.xpomm(m+1)*scm.lt.1.1d0*sgap**2)itypm(1,m+1)=6 !single cut Pomeron
       endif
        
      elseif(jt.eq.5)then                    !cut 'leg'
       nppm(m+1)=1
       ippm(1,m+1)=itt
       if(itypom.eq.4)then              !single cut Pomeron at the end
        if(xpomm(m+1)*scm.le.sgap**2)stop'=4:xpomm(m+1)*scm<sgap**2'
        factor=exp(-vtac(itt))*(1.d0-vvxt)*(1.d0-vvxtl)
        wng=(vtlcng-vtlcpe)*factor/((vtlc0-vtlcpe)*factor
     *  -(vtlc-vtlc0)*(1.d0-factor))
        if(qgran(b10).le.wng.or.wng.lt.0.d0)then
         itypm(1,m+1)=7          !>1 cut Pomerons in the 'handle'
        else
         itypm(1,m+1)=5          !rap-gap at the end
         wgpm(m+1)=(1.d0-factor)/factor*(vtlc-vtlc0)/(vtlc0-vtlcng)
        endif
       elseif(itypom.eq.5)then          !rap-gap at the end (cut or uncut loop)
        if(qgran(b10).le.vtlcpe/vtlcng
     *  .or.xpomm(m+1)*scm.lt.1.1d0*sgap**2)then
         itypm(1,m+1)=4          !single cut Pomeron at the end
        else
         itypm(1,m+1)=7          !>1 cut Pomerons in the 'handle'
        endif
       elseif(itypom.eq.7)then          !>1 cut Pomerons in the 'handle'
        factor=exp(-vtac(itt))*(1.d0-vvxt)*(1.d0-vvxtl)
        wng=vtlcng*factor/(vtlc0*factor-(vtlc-vtlc0)*(1.d0-factor))
        if(qgran(b10).le.wng.or.wng.lt.0.d0
     *  .or.xpomm(m+1)*scm.lt.1.1d0*sgap**2)then
         if(qgran(b10).le.vtlcpe/vtlcng
     *   .or.xpomm(m+1)*scm.lt.1.1d0*sgap**2)then
          itypm(1,m+1)=4         !single cut Pomeron at the end
         else
          itypm(1,m+1)=7         !>1 cut Pomerons in the 'handle'
         endif
        else
         itypm(1,m+1)=5          !rap-gap at the end
         wgpm(m+1)=(1.d0-factor)/factor*(vtlc-vtlc0)/(vtlc0-vtlcng)
        endif
       endif
        
       if(itypm(1,m+1).eq.4)then        !single cut Pomeron at the end
        if(qgran(b10).le.vtlcp/vtlcpe
     *  .or.xpomm(m+1)*scm.lt.1.1d0*sgap**3)itypm(1,m+1)=6 !single cut Pomeron
       endif
      endif

      if(nppm(m+1).eq.1.and.itypm(1,m+1).eq.6)then  !record single cut Pomeron 
       nptg=nptg+1
       if(nptg.gt.legmax)then
        iret=1
        goto 31
       endif
       if(qgran(b10).le.vtlcps/vtlcp
     * .or.xpomm(m+1)*scm.lt.1.1d0*s2min)then
        idpomti(nptg)=0
       else
        idpomti(nptg)=1
       endif
       xpomti(nptg)=xpomm(m+1)
       ipomti(nptg)=itt
       bpomti(nptg,1)=bbt
       bpomti(nptg,2)=(xa(ip,1)+bcoll-xxm(m+1))**2
     * +(xa(ip,2)-yym(m+1))**2
       nppm(m+1)=0
       if(debug.ge.4)write (moniou,217)nptg,itt,bbt,xpomti(nptg)

      elseif(nppm(m+1).gt.1)then  
       vvxtls=vvxtl
       vvxs=(1.d0-vvxp)*(1.d0-vvxpl)*exp(-vpac(ip))
       i=0
26     i=i+1
       ityp=itypm(i,m+1)
       if(ityp.eq.0)then
        iti=ippm(i,m+1)
        bbi=(xb(iti,1)-xxm(m+1))**2+(xb(iti,2)-yym(m+1))**2
        vvxt=0.d0
        vvxtl=0.d0
        vvxt0=0.d0
        if(ia(2).gt.1)then
         do l=1,ia(2)
          if(l.lt.iti)then
           vvxtl=vvxtl+vtac(l)
          elseif(l.gt.iti)then
           vvxt=vvxt+vtac(l)
           vvxt0=vvxt0+vtac0(l)
          endif
         enddo
        endif
        vvxt=1.d0-exp(-vvxt)
        vvxtl=1.d0-exp(-vvxtl)
        vvxt0=1.d0-exp(-vvxt0)
        vvxps=1.d0-vvxs*(1.d0-vvxtl)
        vtlcps=qgfani(xpomm(m+1)*scm,bbi,vvxps,vvxt,vvxtl,iddt(iti)
     *  ,2,16)
        vtlcp=max(vtlcps,vtht(iti,2)
     *  +qgfani(xpomm(m+1)*scm,bbi,vvxps,vvxt,vvxtl,iddt(iti),2,14)) 
        vtacpe=max(vtlcp,vtht(iti,1)
     *  +qgfani(xpomm(m+1)*scm,bbi,vvxps,vvxt0,vvxtl,iddt(iti),2,8))
        vtacng=max(vtacpe,vtht(iti,1)
     *  +qgfani(xpomm(m+1)*scm,bbi,vvxps,vvxt0,vvxtl,iddt(iti),2,7))
         
        aks=qgran(b10)*vtac0(iti)
        if(aks.le.vtlcp.or.xpomm(m+1)*scm.lt.1.1d0*sgap**2)then  
         itypm(i,m+1)=6          !single cut Pomeron
        elseif(aks.lt.vtacpe)then  
         itypm(i,m+1)=-1         !single cut Pomeron in the 'handle'
        elseif(aks.lt.vtacng)then  
         itypm(i,m+1)=2          !>1 cut Pomerons in the 'handle'
        endif
        
        if(itypm(i,m+1).eq.6)then      !record single cut Pomeron
         nptg=nptg+1
         if(nptg.gt.legmax)then
          iret=1
          goto 31
         endif
         if(qgran(b10).le.vtlcps/vtlcp
     *   .or.xpomm(m+1)*scm.lt.1.1d0*s2min)then
          idpomti(nptg)=0
         else
          idpomti(nptg)=1
         endif
         xpomti(nptg)=xpomm(m+1)
         ipomti(nptg)=iti
         bpomti(nptg,1)=bbi
         bpomti(nptg,2)=(xa(ip,1)+bcoll-xxm(m+1))**2
     *   +(xa(ip,2)-yym(m+1))**2
         if(debug.ge.4)write (moniou,217)nptg,iti,bbi,xpomti(nptg)
         nppm(m+1)=nppm(m+1)-1
         if(nppm(m+1).ge.i)then
          do l=i,nppm(m+1)
           ippm(l,m+1)=ippm(l+1,m+1)
           itypm(l,m+1)=itypm(l+1,m+1)
          enddo
         endif
         i=i-1
        endif
       endif
       if(i.lt.nppm(m+1))goto 26
       vvxtl=vvxtls
      endif
       
      if(jt.eq.2.and.qgran(b10).gt.2.d0*vvxtl
     */((1.d0-exp(-vtac(itt)))*(1.d0-vvxtl)+2.d0*vvxtl))then
       if(debug.ge.4)write (moniou,212)
       icdts=iddt(itt)
       do icdt=1,nfock
        iddt(itt)=icdt
        call qgfdf(xxm(m+1),yym(m+1),xpomm(m+1),vpac,vtac,vpht,vtht
     *  ,vpacq,vtacq,vvx,vvxp,vvxt,vvxpl,vvxtl,genhp,genht,ip,itt)
        wdt(icdt,itt)=(1.d0-exp(-vtac(itt)))*(1.d0-vvxtl)
       enddo
       iddt(itt)=icdts
      endif
       
      if(nppm(m+1).ne.0)then
       goto 21
      else
       goto 22
      endif     
31    continue
      if(debug.ge.2)write (moniou,219)nppr,nptg,npin,iret

201   format(2x,'qg3pdf - configuration for multi-Pomeron'
     *,'/diffractive contributions'                        
     */4x,i2,'-th proj. nucleon',2x,i2,'-th targ. nucleon')
202   format(2x,'qg3pdf: problem with initial normalization'
     *,' -> rejection')
203   format(2x,'qg3pdf: normalization of rejection function - ',e10.3)
204   format(2x,'qg3pdf: xpomr=',e10.3,2x,'bbpr=',e10.3,2x,'bbtg=',e10.3
     *,2x,'gb=',e10.3)
205   format(2x,'qg3pdf: xpomr=',e10.3,2x,'bbpr=',e10.3,2x,'bbtg=',e10.3
     *,2x,'xxp=',e10.3,2x,'yyp=',e10.3)
206   format(2x,'qg3pdf: main vertex, nppr0=',i3,2x,'nptg0=',i3)
208   format(2x,'qg3pdf: check',i3,'-th cut fan at ',i2,'-th level,'
     *,' proj. index - ',i3,2x,'b^2=',e10.3)
209   format(2x,'qg3pdf: ',i3,'-th proj. leg, proj. index - ',i3
     *,2x,'b^2=',e10.3,2x,'xpomr=',e10.3,2x,'vvx=',e10.3)
210   format(2x,'qg3pdf: new vertex at ',i3,'-th level')
211   format(2x,'qg3pdf: ',i3,'-th interm. Pomeron'
     */4x,'xpomip=',e10.3,2x,'xpomim=',e10.3,2x,'bpomim=',e10.3)
212   format(2x,'qg3pdf: diffractive cut')
214   format(2x,'qg3pdf: total number of proj. legs - ',i3)
216   format(2x,'qg3pdf: check',i3,'-th cut fan at ',i2,'-th level,'
     *,' targ. index - ',i3,2x,'b^2=',e10.3)
217   format(2x,'qg3pdf: ',i3,'-th targ. leg, targ. index - ',i3
     *,2x,'b^2=',e10.3,2x,'xpomr=',e10.3,2x,'vvx=',e10.3)
219   format(2x,'qg3pdf - end',2x,'number of proj. legs:',i3
     *,2x,'number of targ. legs:',i3
     */4x,'number of interm. Pomerons:',i3,'return flag:',i2)
      return
      end

c=============================================================================
      subroutine qgsha(nbpom,ncola,ncolb,iret)
c-----------------------------------------------------------------------------
c qgsha - inelastic interaction (energy sharing and particle production)
c nbpom - number of Pomeron blocks (nucleon(hadron)-nucleon collisions),
c ncola - number of inel.-wounded proj. nucleons,
c ncolb - number of inel.-wounded targ. nucleons
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208,npbmax=1000,npmax=200,legmax=900,njmax=50000
     *,nfock=3)
      dimension wppr(iapmax),wmtg(iapmax),izp(iapmax),izt(iapmax)
     *,ncola(iapmax),ncolb(iapmax),ncola0(iapmax),ncolb0(iapmax)
     *,lva(iapmax),lvb(iapmax),ipnh(legmax,iapmax),itnh(legmax,iapmax)
     *,wph(npmax),wmh(npmax),xxh(npmax),yyh(npmax),bbh(npmax),wch(3)
     *,wchard(3,npmax),vvxph(npmax),vvxth(npmax)
     *,wppri(legmax,iapmax,2),wmtgi(legmax,iapmax,2)
     *,ptprs(legmax,iapmax,2,2),pttgs(legmax,iapmax,2,2),ptps(2,2)
     *,ptts(2,2),ptprr(iapmax,2),pttgr(iapmax,2),ptpr(2),pttr(2)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr9/  iwp(iapmax),iwt(iapmax),lqa(iapmax),lqb(iapmax)
     *,iprcn(iapmax),itgcn(iapmax),ias(npbmax),ibs(npbmax),nqs(npbmax)
     *,npompr(npbmax),npomtg(npbmax),npomin(npbmax),nnpr(npmax,npbmax)
     *,nntg(npmax,npbmax),ilpr(legmax,npbmax),iltg(legmax,npbmax)
     *,lnpr(legmax,npbmax),lntg(legmax,npbmax),idpom(npmax,npbmax)
     *,idpomp(legmax,npbmax),idpomt(legmax,npbmax),nbpi(legmax,iapmax)
     *,nbti(legmax,iapmax),idnpi(legmax,iapmax),idnti(legmax,iapmax)
     *,nppi(legmax,iapmax),npti(legmax,iapmax),nlpi(legmax,iapmax)
     *,nlti(legmax,iapmax)
      common /qgarr11/ b10
      common /qgarr12/ nsp
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr22/ xppr0(iapmax),xmtg0(iapmax)
      common /qgarr23/ bbpom(npbmax),bpompr(legmax,iapmax,2)
     *,bpomtg(legmax,iapmax,2),xpompr(legmax,iapmax)
     *,xpomtg(legmax,iapmax),xpopin(npmax,npbmax),xpomin(npmax,npbmax)
     *,bpomin(npmax,npbmax)
      common /qgarr26/ factk,fqscal
      common /qgarr28/ arr(5),alpq
      common /qgarr37/ eqj(4,njmax),iqj(njmax),ncj(2,njmax),nj
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /ebal/    ebal0(4),ebal(4)
      external qgran

      if(debug.ge.1)write (moniou,201)nbpom
      
      nsp0=nsp
      iret=0
      nret=0
      s2min=4.d0*fqscal*qt0                      !threshold energy

1     nsp=nsp0   
      nj=0
      nhard=0
      do i=1,4
       ebal(i)=ebal0(i)
      enddo

      if(iret.ne.0)then                          !energy-sharing rejected
       nret=nret+1
       if(nret.gt.20.and.scm.gt.1.d3.or.nret.gt.100)return !overrej.->redo conf.
      endif

c-------------------------------------------------
c initial nucleon (hadron) types
      if(ia(1).ne.1)then
       do i=1,ia(1)
        izp(i)=int(2.5d0+qgran(b10))             !i-th projectile nucleon type
       enddo
      else
       izp(1)=icp                                !projectile hadron type
      endif
      if(ia(2).ne.1)then
       do i=1,ia(2)
        izt(i)=int(2.5d0+qgran(b10))             !i-th target nucleon type
       enddo
      else
       izt(1)=2                                  !target proton
      endif

      do i=1,ia(1)
       lva(i)=0
       ncola0(i)=ncola(i)
       wppr(i)=xppr0(i)*wp0
      enddo
      do i=1,ia(2)
       lvb(i)=0
       ncolb0(i)=ncolb(i)
       wmtg(i)=xmtg0(i)*wm0
      enddo

c-------------------------------------------------
c energy-momentum sharing for hard processes
      if(nbpom.ne.0)then
       do npb=1,nbpom                            !loop over collisions
        ip=ias(npb)                              !proj. index
        it=ibs(npb)                              !targ. index
        icdp=iddp(ip)                            !proj. diffr. eigenstate
        icdt=iddt(it)                            !targ. diffr. eigenstate
        bb=bbpom(npb)                            !b^2 between proj. and targ.
        if(debug.ge.1)write (moniou,211)npb,ip,it,nqs(npb)
     *  ,npomin(npb),npompr(npb),npomtg(npb)

        if(nqs(npb).ne.0)then
         do np=1,nqs(npb)                        !loop over single Pomerons
          idp=idpom(np,npb)                      !Pomeron type
          if(idp.ne.0)then                       !hard rescattering
           nhard=nhard+1
           if(nhard.gt.npmax)stop'nhard>npmax'
           lnp=nnpr(np,npb)                      !proj. constituent index
           lnt=nntg(np,npb)                      !targ. constituent index
           ipnh(lnp,ip)=nhard                    !hard scatt. index
           itnh(lnt,it)=nhard                    !hard scatt. index
           wpi=wppr(ip)                          !LC+ available
           wmi=wmtg(it)                          !LC- available
           call qglchard(wpi,wmi,xph,xmh,bb,0.d0,wch,xxp,yyp,bbi
     *     ,vvxps,vvxts,idp,icdp,icdt,izp(ip),izt(it),ip,it,icz,iret)
           if(iret.ne.0.or.xph*wp0.ge.wppr(ip)
     *     .or.xmh*wm0.ge.wmtg(it).or.xph*xmh*scm.le.s2min)then
            iret=1
            goto 1
           endif
           wph(nhard)=xph*wp0                    !LC+ for hard rescattering
           wmh(nhard)=xmh*wm0                    !LC+ for hard rescattering
           wppr(ip)=wppr(ip)-wph(nhard)          !remaining LC+
           wmtg(it)=wmtg(it)-wmh(nhard)          !remaining LC-
           if(idp.eq.2.or.idp.eq.4)lva(ip)=lnp
           if(idp.eq.3.or.idp.eq.4)lvb(it)=lnt
           do i=1,3
            wchard(i,nhard)=wch(i)
           enddo
           xxh(nhard)=xxp
           yyh(nhard)=yyp
           bbh(nhard)=bbi
           vvxph(nhard)=vvxps
           vvxth(nhard)=vvxts
          endif
         enddo                                   !np-loop
        endif

        if(npompr(npb).ne.0)then
         do np=1,npompr(npb)                     !loop over proj. leg Pomerons
          idp=idpomp(np,npb)                     !Pomeron type
          if(idp.ne.0)then                       !hard rescattering
           if(idp.ne.5.and.idp.ne.6)stop'qgsha: idpomp?!'
           nhard=nhard+1 
           if(nhard.gt.npmax)stop'nhard>npmax'
           ipp=ilpr(np,npb)                      !proj. index
           icdpp=iddp(ipp)                       !proj. diffr. eigenstate
           lnp=lnpr(np,npb)                      !index for proj. constituent
           ipnh(lnp,ipp)=nhard                   !hard scatt. index
           wpi=wppr(ipp)                         !LC+ available
           wmi=xpompr(lnp,ipp)*wm0               !LC- available
           bbp=bpompr(lnp,ipp,1)                 !b^2 to the proj.
           bbt=bpompr(lnp,ipp,2)                 !b^2 to the target
           call qglchard(wpi,wmi,xph,xmh,bbp,bbt,wch,xxp,yyp,bbi,vvxps
     *     ,vvxts,idp,icdpp,icdt,izp(ipp),izt(it),ipp,it,icz,iret)
           if(iret.ne.0.or.xph*wp0.ge.wppr(ipp)
     *     .or.xmh*wm0.ge.wmtg(it).or.xph*xmh*scm.le.s2min)then
            iret=1
            goto 1
           endif
           wph(nhard)=xph*wp0                    !LC+ for hard rescattering
           wmh(nhard)=xmh*wm0                    !LC+ for hard rescattering
           wppr(ipp)=wppr(ipp)-wph(nhard)        !remaining LC+
           wmtg(it)=wmtg(it)-wmh(nhard)          !remaining LC-
           if(idp.eq.6)lva(ipp)=lnp
           do i=1,3
            wchard(i,nhard)=wch(i)
           enddo
           xxh(nhard)=xxp
           yyh(nhard)=yyp
           bbh(nhard)=bbi
           vvxph(nhard)=vvxps
           vvxth(nhard)=vvxts
          endif
         enddo                                   !np-loop
        endif

        if(npomtg(npb).ne.0)then
         do np=1,npomtg(npb)                     !loop over targ. leg Pomerons
          idp=idpomt(np,npb)                     !Pomeron type
          if(idp.ne.0)then                       !hard rescattering
           if(idp.ne.7.and.idp.ne.8)stop'qgsha: idpomt?!'
           nhard=nhard+1 
           if(nhard.gt.npmax)stop'nhard>npmax'
           itt=iltg(np,npb)                      !targ. index
           icdtt=iddt(itt)                       !targ. diffr. eigenstate
           lnt=lntg(np,npb)                      !index for targ. constituent
           itnh(lnt,itt)=nhard                   !hard scatt. index
           wpi=xpomtg(lnt,itt)*wp0               !LC+ available
           wmi=wmtg(itt)                         !LC- available
           bbt=bpomtg(lnt,itt,1)                 !b^2 to the target
           bbp=bpomtg(lnt,itt,2)                 !b^2 to the proj.
           call qglchard(wpi,wmi,xph,xmh,bbp,bbt,wch,xxp,yyp,bbi,vvxps
     *     ,vvxts,idp,icdp,icdtt,izp(ip),izt(itt),ip,itt,icz,iret)
           if(iret.ne.0.or.xph*wp0.ge.wppr(ip)
     *     .or.xmh*wm0.ge.wmtg(itt).or.xph*xmh*scm.le.s2min)then
            iret=1
            goto 1
           endif
           wph(nhard)=xph*wp0                    !LC+ for hard rescattering
           wmh(nhard)=xmh*wm0                    !LC+ for hard rescattering
           wmtg(itt)=wmtg(itt)-wmh(nhard)        !remaining LC-
           wppr(ip)=wppr(ip)-wph(nhard)          !remaining LC+
           if(idp.eq.8)lvb(itt)=lnt
           do i=1,3
            wchard(i,nhard)=wch(i)
           enddo
           xxh(nhard)=xxp
           yyh(nhard)=yyp
           bbh(nhard)=bbi
           vvxph(nhard)=vvxps
           vvxth(nhard)=vvxts
          endif
         enddo                                   !np-loop
        endif
       enddo                                     !npb-loop
      endif

c-------------------------------------------------
c energy-momentum sharing for constituent partons (string ends)
      do ip=1,ia(1)
       if(lqa(ip).gt.1)then
        if(lva(ip).eq.0)then
         nlmax=lqa(ip)-1
        else
         nlmax=lqa(ip)
        endif
        do np=1,nlmax
         if(np.ne.lva(ip))then
          arest=2.d0*(1.d0-alpq)*(nlmax-np+1)
          if(np.lt.lva(ip))arest=arest-2.d0*(1.d0-alpq)
          do i=1,2
           arest=arest-(1.d0-alpq)
           wppri(np,ip,i)=wppr(ip)*qglcgen(arest,lva(ip),icz)
           wppr(ip)=wppr(ip)-wppri(np,ip,i)
          enddo
         endif
        enddo
       endif
      enddo
      
      do it=1,ia(2)
       if(lqb(it).gt.1)then
        if(lvb(it).eq.0)then
         nlmax=lqb(it)-1
        else
         nlmax=lqb(it)
        endif
        do np=1,nlmax
         if(np.ne.lvb(it))then
          arest=2.d0*(1.d0-alpq)*(nlmax-np+1)
          if(np.lt.lvb(it))arest=arest-2.d0*(1.d0-alpq)
          do i=1,2
           arest=arest-(1.d0-alpq)
           wmtgi(np,it,i)=wmtg(it)*qglcgen(arest,lvb(it),2)
           wmtg(it)=wmtg(it)-wmtgi(np,it,i)
          enddo
         endif
        enddo
       endif
      enddo

c-------------------------------------------------
c transverse momenta for constituent partons (string ends)
      do ip=1,ia(1)
       do i=1,2
        ptpr(i)=0.d0
       enddo
       if(lqa(ip).ne.0)then
        if(lva(ip).eq.0)then
         nstr=2*lqa(ip)-1
        else
         nstr=2*lqa(ip)-2
        endif
        do np=1,lqa(ip)
        do istr=1,2
         if(np.ne.lva(ip).and.(np.lt.lqa(ip).or.lva(ip).ne.0
     *   .or.istr.eq.1))then
          nstr=nstr-1
          do i=1,2
           call qgptstr(pts,ptpr(i),nstr,icz)
           ptprs(np,ip,istr,i)=pts
          enddo
         endif
        enddo
        enddo
       endif
       do i=1,2
        ptprr(ip,i)=-ptpr(i)
       enddo
      enddo

      do it=1,ia(2)
       do i=1,2
        pttr(i)=0.d0
       enddo
       if(lqb(it).ne.0)then
        if(lvb(it).eq.0)then
         nstr=2*lqb(it)-1
        else
         nstr=2*lqb(it)-2
        endif
        do np=1,lqb(it)
        do istr=1,2
         if(np.ne.lvb(it).and.(np.lt.lqb(it).or.lvb(it).ne.0
     *   .or.istr.eq.1))then
          nstr=nstr-1
          do i=1,2
           call qgptstr(pts,pttr(i),nstr,2)
           pttgs(np,it,istr,i)=pts
          enddo
         endif
        enddo
        enddo
       endif
       do i=1,2
        pttgr(it,i)=-pttr(i)
       enddo
      enddo

c-------------------------------------------------
c treatment of low mass diffraction
      if(debug.ge.1)write (moniou,207)
      do ip=1,ia(1)                              !loop over proj. nucleons
       if(iwp(ip).eq.2)then                      !diffraction dissociation
        it=iprcn(ip)
        if(debug.ge.2)write (moniou,208)ip,it
        if(iwt(it).eq.2)then
         call qgdifr(wppr(ip),wmtg(it),ptprr(ip,1),ptprr(ip,2)
     *   ,pttgr(it,1),pttgr(it,2),izp(ip),izt(it),2,2,iret)
        elseif(iwt(it).eq.-1)then
         call qgdifr(wppr(ip),wmtg(it),ptprr(ip,1),ptprr(ip,2)
     *   ,pttgr(it,1),pttgr(it,2),izp(ip),izt(it),2,0,iret)
        elseif(iwt(it).gt.0)then
         call qgdifr(wppr(ip),wmtg(it),ptprr(ip,1),ptprr(ip,2)
     *   ,pttgr(it,1),pttgr(it,2),izp(ip),izt(it),2,-1,iret)
        else
         stop'wrong connection for diffraction'
        endif
        if(iret.ne.0)goto 1
       endif
      enddo

      do it=1,ia(2)                              !loop over targ. nucleons
       if(iwt(it).eq.2)then                      !diffraction dissociation
        ip=itgcn(it)
        if(debug.ge.2)write (moniou,209)it,ip
        if(iwp(ip).eq.-1)then
         call qgdifr(wppr(ip),wmtg(it),ptprr(ip,1),ptprr(ip,2)
     *   ,pttgr(it,1),pttgr(it,2),izp(ip),izt(it),0,2,iret)
        elseif(iwp(ip).gt.0.and.iwp(ip).ne.2)then
         call qgdifr(wppr(ip),wmtg(it),ptprr(ip,1),ptprr(ip,2)
     *   ,pttgr(it,1),pttgr(it,2),izp(ip),izt(it),-1,2,iret)
        endif
        if(iret.ne.0)goto 1
       endif
      enddo

c-------------------------------------------------
c treatment of particle production
      if(nbpom.ne.0)then
       if(debug.ge.1)write (moniou,210)
       do npb=1,nbpom                            !loop over collisions
        ip=ias(npb)                              !proj. index
        it=ibs(npb)                              !targ. index
        icdp=iddp(ip)                            !proj. diffr. eigenstate
        icdt=iddt(it)                            !targ. diffr. eigenstate
        ncola0(ip)=ncola0(ip)-1
        ncolb0(it)=ncolb0(it)-1
        if(debug.ge.1)write (moniou,211)npb,ip,it,nqs(npb),npomin(npb)
     *  ,npompr(npb),npomtg(npb)

        if(npomin(npb).ne.0)then
         do n=1,npomin(npb)                      !loop over interm. Pomerons
          wpi=xpopin(n,npb)*wp0                  !LC+ for the Pomeron
          wmi=xpomin(n,npb)*wm0                  !LC- for the Pomeron
          if(debug.ge.2)write (moniou,212)n,wpi,wmi
          do i=1,2
           do j=1,2
            ptps(j,i)=0.d0
            ptts(j,i)=0.d0
           enddo
           ptpr(i)=0.d0
           pttr(i)=0.d0
          enddo
          call qgstr(wpi,0.d0,wmi,0.d0,wppr(ip),wmtg(it)
     *    ,ptps,ptts,ptpr,pttr,0,0)
         enddo
        endif

        if(npompr(npb).ne.0)then
         do l=1,npompr(npb)                      !loop over proj. leg Pomerons
          ipp=ilpr(l,npb)                        !proj. index
          icdpp=iddp(ipp)                        !proj. diffr. eigenstate
          lnp=lnpr(l,npb)                        !index for proj. constituent
          if(lva(ipp).eq.0.and.lnp.lt.lqa(ipp)
     *    .or.lva(ipp).ne.0.and.lnp.ne.lva(ipp))then
           wp1=wppri(lnp,ipp,1)                  !LC+ for 1st string end
           wp2=wppri(lnp,ipp,2)                  !LC+ for 2nd string end
           wmi=wm0*xpompr(lnp,ipp)               !LC- for the Pomeron
           idp=idpomp(l,npb)
           do i=1,2
            do j=1,2
             ptps(j,i)=ptprs(lnp,ipp,j,i)
             ptts(j,i)=0.d0
            enddo
            ptpr(i)=ptprr(ipp,i)
            pttr(i)=0.d0
           enddo
           if(idp.eq.0)then
            call qgstr(wp1,wp2,wmi,0.d0,wppr(ipp),wmtg(it)
     *      ,ptps,ptts,ptpr,pttr,1,0)
           else        !QCD evolution and hadronization for semi-hard Pomeron
            nhard=ipnh(lnp,ipp)
            do i=1,3
             wch(i)=wchard(i,nhard)
            enddo
            call qghot(wph(nhard),wmh(nhard),wp1,wp2,wmi,0.d0,wch
     *      ,ptps,ptts,xxh(nhard),yyh(nhard),bbh(nhard),vvxph(nhard)
     *      ,vvxth(nhard),idp,izp(ipp),izt(it),icdpp,icdt,0,0,0,0
     *      ,ipp,it,iret)
            if(iret.ne.0)goto 1
           endif
           do i=1,2
            ptprr(ipp,i)=ptpr(i)
           enddo
          endif
         enddo
        endif

        if(npomtg(npb).ne.0)then
         do l=1,npomtg(npb)                      !loop over targ. leg Pomerons
          itt=iltg(l,npb)                        !targ. index
          icdtt=iddt(itt)                        !targ. diffr. eigenstate
          lnt=lntg(l,npb)                        !index for targ. constituent
          if(lvb(itt).eq.0.and.lnt.lt.lqb(itt)
     *    .or.lvb(itt).ne.0.and.lnt.ne.lvb(itt))then
           wpi=wp0*xpomtg(lnt,itt)               !LC+ for the Pomeron
           wm1=wmtgi(lnt,itt,1)                  !LC- for 1st string end
           wm2=wmtgi(lnt,itt,2)                  !LC- for 2nd string end
           idp=idpomt(l,npb)
           do i=1,2
            do j=1,2
             ptps(j,i)=0.d0
             ptts(j,i)=pttgs(lnt,itt,j,i)
            enddo
            ptpr(i)=0.d0
            pttr(i)=pttgr(itt,i)
           enddo
           if(idp.eq.0)then
            call qgstr(wpi,0.d0,wm1,wm2,wppr(ip),wmtg(itt)
     *      ,ptps,ptts,ptpr,pttr,0,1)
           else        !QCD evolution and hadronization for semi-hard Pomeron
            nhard=itnh(lnt,itt)
            do i=1,3
             wch(i)=wchard(i,nhard)
            enddo
            call qghot(wph(nhard),wmh(nhard),wpi,0.d0,wm1,wm2,wch
     *      ,ptps,ptts,xxh(nhard),yyh(nhard),bbh(nhard),vvxph(nhard)
     *      ,vvxth(nhard),idp,izp(ip),izt(itt),icdp,icdtt,0,0,0,0
     *      ,ip,itt,iret)
            if(iret.ne.0)goto 1
           endif
           do i=1,2
            pttgr(itt,i)=pttr(i)
           enddo
          endif
         enddo
        endif
        
        if(nqs(npb).ne.0)then
         do n=1,nqs(npb)                         !loop over single Pomerons
          lnp=nnpr(n,npb)                        !index for proj. constituent
          lnt=nntg(n,npb)                        !index for targ. constituent
          if((lva(ip).eq.0.and.lnp.lt.lqa(ip).or.lva(ip).ne.0
     *    .and.lnp.ne.lva(ip)).and.(lvb(it).eq.0.and.lnt.lt.lqb(it)
     *    .or.lvb(it).ne.0.and.lnt.ne.lvb(it)))then
           wp1=wppri(lnp,ip,1)                   !LC+ for 1st string end
           wp2=wppri(lnp,ip,2)                   !LC+ for 2nd string end
           wm1=wmtgi(lnt,it,1)                   !LC- for 1st string end
           wm2=wmtgi(lnt,it,2)                   !LC- for 2nd string end
           idp=idpom(n,npb)
           do i=1,2
            do j=1,2
             ptps(j,i)=ptprs(lnp,ip,j,i)
             ptts(j,i)=pttgs(lnt,it,j,i)
            enddo
            ptpr(i)=ptprr(ip,i)
            pttr(i)=pttgr(it,i)
           enddo
           if(idp.eq.0)then
            call qgstr(wp1,wp2,wm1,wm2,wppr(ip),wmtg(it)
     *      ,ptps,ptts,ptpr,pttr,1,1)
           else        !QCD evolution and hadronization for semi-hard Pomeron
            nhard=ipnh(lnp,ip)
            do i=1,3
             wch(i)=wchard(i,nhard)
            enddo
            call qghot(wph(nhard),wmh(nhard),wp1,wp2,wm1,wm2,wch
     *      ,ptps,ptts,xxh(nhard),yyh(nhard),bbh(nhard),vvxph(nhard)
     *      ,vvxth(nhard),idp,izp(ip),izt(it),icdp,icdt,0,0,0,0
     *      ,ip,it,iret)
            if(iret.ne.0)goto 1
           endif
           do i=1,2
            ptprr(ip,i)=ptpr(i)
            pttgr(it,i)=pttr(i)
           enddo
          endif
         enddo
        endif

        if(lqa(ip).eq.0.and.ncola0(ip).eq.0
     *  .and.lqb(it).eq.0.and.ncolb0(it).eq.0)then
         call qgdifr(wppr(ip),wmtg(it),ptprr(ip,1),ptprr(ip,2)
     *   ,pttgr(it,1),pttgr(it,2),izp(ip),izt(it)
     *   ,iwp(ip)-1,iwt(it)-1,iret)
        elseif(lqa(ip).eq.0.and.ncola0(ip).eq.0
     *  .and.(lqb(it).ne.0.or.ncolb0(it).ne.0))then
         call qgdifr(wppr(ip),wmtg(it),ptprr(ip,1),ptprr(ip,2)
     *   ,pttgr(it,1),pttgr(it,2),izp(ip),izt(it),iwp(ip)-1,-1,iret)
        elseif((lqa(ip).ne.0.or.ncola0(ip).ne.0)
     *  .and.lqb(it).eq.0.and.ncolb0(it).eq.0)then
         call qgdifr(wppr(ip),wmtg(it),ptprr(ip,1),ptprr(ip,2)
     *   ,pttgr(it,1),pttgr(it,2),izp(ip),izt(it),-1,iwt(it)-1,iret)
        endif
        if(iret.ne.0)goto 1
       enddo                                     !end of collision loop
      endif
      
      do ip=1,ia(1)
       if(lqa(ip).ne.0)then
        if(lva(ip).ne.0)then
         lnp=lva(ip)
        else
         lnp=lqa(ip)
        endif
        wp1=wppr(ip)
        wp2=0.d0
        
        if(idnpi(lnp,ip).eq.0)then
         npb=nbpi(lnp,ip)                        !coll. index
         np=nppi(lnp,ip)                         !Pomeron index
         idp=idpom(np,npb)
         it=ibs(npb)                             !targ. index
         lnt=nntg(np,npb)                        !index for targ. constituent
          if(lvb(it).eq.0.and.lnt.eq.lqb(it)
     *   .or.lvb(it).ne.0.and.lnt.eq.lvb(it))then
          iqvt=1
          wm1=wmtg(it)
          wm2=0.d0
         else
          iqvt=0
          wm1=wmtgi(lnt,it,1)                    !LC- for 1st string end
          wm2=wmtgi(lnt,it,2)                    !LC- for 2nd string end
         endif
        else
         iqvt=0
         npb=nbpi(lnp,ip)                        !coll. index
         np=nlpi(lnp,ip)                         !Pomeron index
         idp=idpomp(np,npb)
         it=ibs(npb)                             !targ. index
         wm1=wm0*xpompr(lnp,ip)                  !LC- for the Pomeron
         wm2=0.d0
        endif
        icdp=iddp(ip)                           !proj. diffr. eigenstate
        icdt=iddt(it)                           !targ. diffr. eigenstate

        if(lva(ip).eq.0)then
         do i=1,2
          ptps(1,i)=ptprs(lnp,ip,1,i)
          ptps(2,i)=ptprr(ip,i)
         enddo
        else
         do i=1,2
          ptps(1,i)=0.d0
          ptps(2,i)=ptprr(ip,i)
         enddo
        endif
        if(idnpi(lnp,ip).eq.0)then
         if(iqvt.eq.1.and.lvb(it).eq.0)then
          do i=1,2
           ptts(1,i)=pttgs(lnt,it,1,i)
           ptts(2,i)=pttgr(it,i)
          enddo
         elseif(iqvt.eq.1.and.lvb(it).ne.0)then
          do i=1,2
           ptts(1,i)=0.d0
           ptts(2,i)=pttgr(it,i)
          enddo
         else
          do i=1,2
           do j=1,2
            ptts(j,i)=pttgs(lnt,it,j,i)
           enddo
          enddo
         endif
        else
         do i=1,2
          do j=1,2
           ptts(j,i)=0.d0
          enddo
         enddo
        endif

        if(idp.eq.0)then
         call qgvstr(wp1,wp2,wm1,wm2,ptps,ptts,izp(ip),izt(it),1,iqvt
     *   ,lqa(ip)+lva(ip),lqb(it)+lvb(it),idnpi(lnp,ip),icz,iret)
        else        !QCD evolution and hadronization for semi-hard Pomeron
         nhard=ipnh(lnp,ip)
         do i=1,3
          wch(i)=wchard(i,nhard)
         enddo
         call qghot(wph(nhard),wmh(nhard),wp1,wp2,wm1,wm2,wch,ptps,ptts
     *   ,xxh(nhard),yyh(nhard),bbh(nhard),vvxph(nhard),vvxth(nhard)
     *   ,idp,izp(ip),izt(it),icdp,icdt,1,iqvt,lqa(ip)+lva(ip)
     *   ,lqb(it)+lvb(it),ip,it,iret)
        endif
        if(iret.ne.0)goto 1
       endif
      enddo

      do 2 it=1,ia(2)
       if(lqb(it).ne.0)then
        if(lvb(it).ne.0)then
         lnt=lvb(it)
        else
         lnt=lqb(it)
        endif
        wm1=wmtg(it)
        wm2=0.d0
        iqvp=0
        
        if(idnti(lnt,it).eq.0)then
         npb=nbti(lnt,it)                        !coll. index
         np=npti(lnt,it)                         !Pomeron index
         idp=idpom(np,npb)
         ip=ias(npb)                             !targ. index
         lnp=nnpr(np,npb)                        !index for targ. constituent
          if(lva(ip).ne.0.and.lnp.eq.lva(ip)
     *   .or.lva(ip).eq.0.and.lnp.eq.lqa(ip))goto 2
         wp1=wppri(lnp,ip,1)                     !LC+ for 1st string end
         wp2=wppri(lnp,ip,2)                     !LC+ for 2nd string end
        else
         npb=nbti(lnt,it)                        !coll. index
         np=nlti(lnt,it)                         !Pomeron index
         idp=idpomt(np,npb)
         ip=ias(npb)                             !targ. index
         wp1=wp0*xpomtg(lnt,it)                  !LC+ for the Pomeron
         wp2=0.d0
        endif
        icdp=iddp(ip)                           !proj. diffr. eigenstate
        icdt=iddt(it)                           !targ. diffr. eigenstate

        if(lvb(it).eq.0)then
         do i=1,2
          ptts(1,i)=pttgs(lnt,it,1,i)
          ptts(2,i)=pttgr(it,i)
         enddo
        else
         do i=1,2
          ptts(1,i)=0.d0
          ptts(2,i)=pttgr(it,i)
         enddo
        endif
        if(idnti(lnt,it).eq.0)then
         do i=1,2
          do j=1,2
           ptps(j,i)=ptprs(lnp,ip,j,i)
          enddo
         enddo
        else
         do i=1,2
          do j=1,2
           ptps(j,i)=0.d0
          enddo
         enddo
        endif

        if(idp.eq.0)then
         call qgvstr(wp1,wp2,wm1,wm2,ptps,ptts,izp(ip),izt(it),iqvp,1
     *   ,lqa(ip)+lva(ip),lqb(it)+lvb(it),idnti(lnt,it),icz,iret)
        else        !QCD evolution and hadronization for semi-hard Pomeron
         nhard=itnh(lnt,it)
         do i=1,3
          wch(i)=wchard(i,nhard)
         enddo
         call qghot(wph(nhard),wmh(nhard),wp1,wp2,wm1,wm2,wch,ptps,ptts
     *   ,xxh(nhard),yyh(nhard),bbh(nhard),vvxph(nhard),vvxth(nhard)
     *   ,idp,izp(ip),izt(it),icdp,icdt,iqvp,1,lqa(ip)+lva(ip)
     *   ,lqb(it)+lvb(it),ip,it,iret)
        endif
        if(iret.ne.0)goto 1
       endif
2     continue

      if(nj.ne.0)then                            !parton color connections
       if(debug.ge.1)write (moniou,216)nj
       call qgjarr(jfl)
       if(jfl.eq.0)then
        iret=1
        goto 1
       endif
       if(debug.ge.1)write (moniou,217)
       call qgxjet(iret)                         !jet hadronization
       if(iret.ne.0)goto 1
      endif
      if(debug.ge.1)write (moniou,218)

201   format(2x,'qgsha - inelastic interaction, N of Pomeron blocks:'
     *,i4)
207   format(2x,'qgsha: treatment of low mass diffraction')
208   format(2x,'qgsha: diffraction of ',i3,'-th proj. nucleon,'
     *,' recoil of ',i3,'-th targ. nucleon')
209   format(2x,'qgsha: diffraction of ',i3,'-th targ. nucleon,'
     *,' recoil of ',i3,'-th proj. nucleon')
210   format(2x,'qgsha: particle production for all cut Pomerons')
211   format(2x,'qgsha: ',i4,'-th collision,  proj. index - ',i3
     *,2x,'targ. index - ',i3/4x,'N of single Pomerons - ',i3
     *,2x,' N of interm. Pomerons - ',i3
     */4x,'N of proj. legs - ',i3,2x,'N of targ. legs - ',i3)
212   format(2x,'qgsha: particle production for '
     *,i3,'-th interm. Pomeron'
     */4x,'light cone momenta for the Pomeron:',2e10.3)
216   format(2x,'qgsha: arrangement of color connections for '
     *,i5,' final partons')
217   format(2x,'qgsha: jet hadronization')
218   format(2x,'qgsha - end')
      return
      end

c=============================================================================
      subroutine qglchard(wpi,wmi,xph,xmh,bbp0,bbt0,wch,xxh,yyh
     *,bbi,vvxps,vvxts,idp,icdp,icdt,izp,izt,ip,it,icz,iret)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208,nfock=3)
      dimension vpac(iapmax),vtac(iapmax),vpht(iapmax,2),vtht(iapmax,2)
     *,vpacq(iapmax),vtacq(iapmax),wch(3),dvhtp(2,2),dvhtt(2,2)
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),bcoll
      common /qgarr11/ b10
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr25/ ahv(3)
      common /qgarr26/ factk,fqscal
      common /qgarr38/ htfac
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      iret=0
      s2min=4.d0*fqscal*qt0                      !threshold for hard scattering
      xph=0.d0
      xmh=0.d0
      xxh=0.d0
      yyh=0.d0
      bbi=0.d0
      vvxps=0.d0
      vvxts=0.d0
      if(wpi*wmi.le.s2min)then
       iret=1
       return
      endif
      
      xpmax=wpi/wp0
      xmmax=wmi/wm0
      zmin=s2min/wp0/wm0
      zmax=xpmax*xmmax
      xpmin=xpmax*zmin/zmax
      xmmin=xmmax*zmin/zmax
      if(idp.eq.4)then
       sjqq=qgjit(qt0,qt0,zmax*wp0*wm0,2,2)      !inclusive qq cross-section
       sjqqq=qgjit(qt0,qt0,zmax*wp0*wm0,3,2)
       sjqqa=qgjit(qt0,qt0,zmax*wp0*wm0,4,2)
       qvp1=qggrv(xpmax,qt0,icz,1)*(1.d0-xpmin)**ahv(icz)
       qvp2=qggrv(xpmax,qt0,icz,2)*(1.d0-xpmin)**ahv(icz)
       qvt1=qggrv(xmmax,qt0,2,1)*(1.d0-xmmin)**ahv(2)
       qvp2=qggrv(xmmax,qt0,2,2)*(1.d0-xmmin)**ahv(2)
       if(icz.eq.1)then
        if(izp.eq.1.and.izt.eq.2.or.izp.eq.-1.and.izt.eq.3)then
         gb0=qvp1*qvt1*sjqqq+qvp2*qvt2*sjqqa+(qvp1*qvt2+qvp2*qvt1)*sjqq
        else
         gb0=qvp1*qvt2*sjqqq+qvp2*qvt1*sjqqa+(qvp1*qvt1+qvp2*qvt2)*sjqq
        endif
       elseif(icz.eq.2)then
        if(izp.eq.izt)then
         gb0=(qvp1*qvt1+qvp2*qvt2)*sjqqq+(qvp1*qvt2+qvp2*qvt1)*sjqq
        elseif(izp.gt.0)then
         gb0=(qvp1*qvt2+qvp2*qvt1)*sjqqq+(qvp1*qvt1+qvp2*qvt2)*sjqq
        elseif(izp+izt.eq.0)then
         gb0=(qvp1*qvt1+qvp2*qvt2)*sjqqa+(qvp1*qvt2+qvp2*qvt1)*sjqq
        else
         gb0=(qvp1*qvt2+qvp2*qvt1)*sjqqa+(qvp1*qvt1+qvp2*qvt2)*sjqq
        endif
       elseif(icz.eq.3)then
        if(izp.eq.4.and.izt.eq.2.or.izp.eq.5.and.izt.eq.3)then
         gb0=qvp1*qvt1*sjqqq+(qvp1*qvt2+qvp2*qvt1+qvp2*qvt2)*sjqq
        elseif(izp.gt.0)then
         gb0=qvp1*qvt2*sjqqq+(qvp1*qvt1+qvp2*qvt1+qvp2*qvt2)*sjqq
        elseif(izp.eq.-4.and.izt.eq.2.or.izp.eq.-5.and.izt.eq.3)then
         gb0=qvp1*qvt1*sjqqa+(qvp1*qvt2+qvp2*qvt1+qvp2*qvt2)*sjqq
        else
         gb0=qvp1*qvt2*sjqqa+(qvp1*qvt1+qvp2*qvt1+qvp2*qvt2)*sjqq
        endif
       endif
     
1      zh=zmin*(1.d0-qgran(b10)*(1.d0-(zmax/zmin)**(delh+.4d0)))
     * **(1.d0/(delh+.4d0))
       if(qgran(b10).gt.dlog(zmax/zh)/dlog(zmax/zmin))goto 1
       xph=xpmax*(zh/zmax)**qgran(b10)           !LC+ momentum fraction
       xmh=zh/xph                                !LC- momentum fraction
       sjqq=qgjit(qt0,qt0,zh*wp0*wm0,2,2)        !inclusive qq cross-section
       sjqqq=qgjit(qt0,qt0,zh*wp0*wm0,3,2)
       sjqqa=qgjit(qt0,qt0,zh*wp0*wm0,4,2)
       qvp1=qggrv(xph,qt0,icz,1)*(1.d0-xph)**ahv(icz)
       qvp2=qggrv(xph,qt0,icz,2)*(1.d0-xph)**ahv(icz)
       qvt1=qggrv(xmh,qt0,2,1)*(1.d0-xmh)**ahv(2)
       qvt2=qggrv(xmh,qt0,2,2)*(1.d0-xmh)**ahv(2)
       if(icz.eq.1)then
        if(izp.eq.1.and.izt.eq.2.or.izp.eq.-1.and.izt.eq.3)then
         wch(1)=qvp1*qvt1*sjqqq
         wch(2)=qvp2*qvt1*sjqq
         wch(3)=qvp1*qvt2*sjqq
         tot=wch(1)+wch(2)+wch(3)+qvp2*qvt2*sjqqa
        else
         wch(1)=qvp1*qvt1*sjqq
         wch(2)=qvp2*qvt1*sjqqa
         wch(3)=qvp1*qvt2*sjqqq
         tot=wch(1)+wch(2)+wch(3)+qvp2*qvt2*sjqq
        endif
       elseif(icz.eq.2)then
        if(izp.eq.izt)then
         wch(1)=qvp1*qvt1*sjqqq
         wch(2)=qvp2*qvt1*sjqq
         wch(3)=qvp1*qvt2*sjqq
         tot=wch(1)+wch(2)+wch(3)+qvp2*qvt2*sjqqq
        elseif(izp.gt.0)then
         wch(1)=qvp1*qvt1*sjqq
         wch(2)=qvp2*qvt1*sjqqq
         wch(3)=qvp1*qvt2*sjqqq
         tot=wch(1)+wch(2)+wch(3)+qvp2*qvt2*sjqq
        elseif(izp+izt.eq.0)then
         wch(1)=qvp1*qvt1*sjqqa
         wch(2)=qvp2*qvt1*sjqq
         wch(3)=qvp1*qvt2*sjqq
         tot=wch(1)+wch(2)+wch(3)+qvp2*qvt2*sjqqa
        else
         wch(1)=qvp1*qvt1*sjqq
         wch(2)=qvp2*qvt1*sjqqa
         wch(3)=qvp1*qvt2*sjqqa
         tot=wch(1)+wch(2)+wch(3)+qvp2*qvt2*sjqq
        endif
       elseif(icz.eq.3)then
        if(izp.eq.4.and.izt.eq.2.or.izp.eq.5.and.izt.eq.3)then
         wch(1)=qvp1*qvt1*sjqqq
         wch(2)=qvp2*qvt1*sjqq
         wch(3)=qvp1*qvt2*sjqq
         tot=wch(1)+wch(2)+wch(3)+qvp2*qvt2*sjqq
        elseif(izp.gt.0)then
         wch(1)=qvp1*qvt1*sjqq
         wch(2)=qvp2*qvt1*sjqq
         wch(3)=qvp1*qvt2*sjqqq
         tot=wch(1)+wch(2)+wch(3)+qvp2*qvt2*sjqq
        elseif(izp.eq.-4.and.izt.eq.2.or.izp.eq.-5.and.izt.eq.3)then
         wch(1)=qvp1*qvt1*sjqqa
         wch(2)=qvp2*qvt1*sjqq
         wch(3)=qvp1*qvt2*sjqq
         tot=wch(1)+wch(2)+wch(3)+qvp2*qvt2*sjqq
        else
         wch(1)=qvp1*qvt1*sjqq
         wch(2)=qvp2*qvt1*sjqq
         wch(3)=qvp1*qvt2*sjqqa
         tot=wch(1)+wch(2)+wch(3)+qvp2*qvt2*sjqq
        endif
       endif
       gb=tot/(zh/zmax)**(delh+.4d0)/gb0/1.4d0
       if(qgran(b10).gt.gb)goto 1
       do i=1,3
        wch(i)=wch(i)/tot
       enddo
       
      else
       xxp0=xa(ip,1)+bcoll
       yyp0=xa(ip,2)
       xxt0=xb(it,1)
       yyt0=xb(it,2)
       if(idp.gt.4)then
        call qgbdef(bbp0,bbt0,xxp0,yyp0,xxt0,yyt0,xxp,yyp
     *  ,int(1.5d0+qgran(b10)))
        if(idp.eq.5.or.idp.eq.6)then
         bb0=bbp0
         xxt0=xxp
         yyt0=yyp
        elseif(idp.eq.7.or.idp.eq.8)then
         bb0=bbt0
         xxp0=xxp
         yyp0=yyp
        endif
       else
        bb0=bbp0
       endif

       if(idp.eq.1)then
        rpmax=(rq(icdp,icz)+rq(icdt,2)-alfp*dlog(zmin))*4.d0*.0389d0
        rpmin=(rq(icdp,icz)+rq(icdt,2)-alfp*dlog(zmax))*4.d0*.0389d0
       elseif(idp.eq.2)then
        rpmax=(rq(icdp,icz)+rq(icdt,2)-alfp*dlog(xmmin))
     *  *4.d0*.0389d0
        rpmin=(rq(icdp,icz)+rq(icdt,2)-alfp*dlog(xmmax))*4.d0*.0389d0
       elseif(idp.eq.3)then
        rpmax=(rq(icdp,icz)+rq(icdt,2)-alfp*dlog(xpmin))
     *  *4.d0*.0389d0
        rpmin=(rq(icdp,icz)+rq(icdt,2)-alfp*dlog(xpmax))*4.d0*.0389d0
       elseif(idp.eq.5)then
        rpmax=(rq(icdp,icz)-alfp*dlog(xpmin))*4.d0*.0389d0
        rpmin=(rq(icdp,icz)-alfp*dlog(xpmax))*4.d0*.0389d0
       elseif(idp.eq.6)then
        rpmax=(rq(icdp,icz)-alfp*dlog(zmin/zmax))*4.d0*.0389d0
        rpmin=rq(icdp,icz)*4.d0*.0389d0
       elseif(idp.eq.7)then
        rpmax=(rq(icdt,2)-alfp*dlog(xmmin))*4.d0*.0389d0
        rpmin=(rq(icdt,2)-alfp*dlog(xmmax))*4.d0*.0389d0
       elseif(idp.eq.8)then
        rpmax=(rq(icdt,2)-alfp*dlog(zmin/zmax))*4.d0*.0389d0
        rpmin=rq(icdt,2)*4.d0*.0389d0
       else
        stop'qglchard: wrong idp!'
       endif
       gb0=exp(-bb0/rpmax)/rpmin
       
       sjqq=qgjit(qt0,qt0,zmax*wp0*wm0,2,2)      !inclusive qq cross-section
       sjqg=qgjit(qt0,qt0,zmax*wp0*wm0,1,2)      !inclusive qg cross-section
       sjgg=qgjit(qt0,qt0,zmax*wp0*wm0,1,1)      !inclusive gg cross-section
       sjqqq=qgjit(qt0,qt0,zmax*wp0*wm0,3,2)
       sjqqa=qgjit(qt0,qt0,zmax*wp0*wm0,4,2)
       if(idp.eq.1.or.idp.eq.3.or.idp.eq.5)then
        glu1=qgppdc(xpmin,0,icdp,icz)*fp(icdp,icz)*rr*4.d0*.0389d0
        sea1=qgppdc(xpmin,1,icdp,icz)*fp(icdp,icz)*rr*4.d0*.0389d0
       elseif(idp.eq.2.or.idp.eq.6)then
        qv1=(qggrv(xpmax,qt0,icz,1)+qggrv(xpmax,qt0,icz,2))
     *  *(1.d0-xpmin)**ahv(icz)/xpmax**(delh+.4d0)
       elseif(idp.eq.7.or.idp.eq.8)then
        glu1=qgppdi(xpmin/xpmax,0)*g3p*rr*4.d0*.0389d0
        sea1=qgppdi(xpmin/xpmax,1)*g3p*rr*4.d0*.0389d0
       endif
       if(idp.eq.1.or.idp.eq.2.or.idp.eq.7)then
        glu2=qgppdc(xmmin,0,icdt,2)*fp(icdt,2)*rr*4.d0*.0389d0
        sea2=qgppdc(xmmin,1,icdt,2)*fp(icdt,2)*rr*4.d0*.0389d0
       elseif(idp.eq.3.or.idp.eq.8)then
        qv2=(qggrv(xmmax,qt0,2,1)+qggrv(xmmax,qt0,2,2))
     *  *(1.d0-xmmin)**ahv(2)/xmmax**(delh+.4d0)
       elseif(idp.eq.5.or.idp.eq.6)then
        glu2=qgppdi(xmmin/xmmax,0)*g3p*rr*4.d0*.0389d0
        sea2=qgppdi(xmmin/xmmax,1)*g3p*rr*4.d0*.0389d0
       endif
       if(idp.eq.1.or.idp.eq.5.or.idp.eq.7)then
        gb0=gb0*(glu1*glu2*sjgg+sea1*sea2*(sjqq/1.5d0+sjqqq/6.d0
     *  +sjqqa/6.d0)+(sea1*glu2+glu1*sea2)*sjqg)
       elseif(idp.eq.2.or.idp.eq.6)then
        gb0=gb0*qv1*(glu2*sjqg+sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0))
       elseif(idp.eq.3.or.idp.eq.8)then
        gb0=gb0*qv2*(glu1*sjqg+sea1*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0))
       endif 
       nret=0
       nren=0
       gbmax=0.d0
     
2      zh=zmin*(1.d0-qgran(b10)*(1.d0-(zmax/zmin)**(delh-dels)))
     * **(1.d0/(delh-dels))
       if(idp.eq.1.or.idp.eq.5.or.idp.eq.7)then
        if(qgran(b10).gt.dlog(zmax/zh)/dlog(zmax/zmin))goto 2
        xph=xpmax*(zh/zmax)**qgran(b10)          !LC+ momentum fraction
        xmh=zh/xph                               !LC- momentum fraction
        gb=(zh/zmax)**(dels-delh)
       else
        gbz=(1.d0-(zh/zmax)**(dels+.4d0))
     *  /(1.d0-(zmin/zmax)**(dels+.4d0))
        if(qgran(b10).gt.gbz)goto 2
        xph=(1.d0-qgran(b10)*(1.d0-(zh/zmax)**(dels+.4d0)))
     *  **(1.d0/(dels+.4d0))
        if(idp.eq.2.or.idp.eq.6)then
         xph=xpmax*xph                           !LC+ momentum fraction
         xmh=zh/xph                              !LC- momentum fraction
         gb=(zh/zmax)**(dels-delh)/xph**(dels+.4d0)
        elseif(idp.eq.3.or.idp.eq.8)then
         xmh=xmmax*xph                           !LC- momentum fraction
         xph=zh/xmh                              !LC+ momentum fraction   
         gb=(zh/zmax)**(dels-delh)/xmh**(dels+.4d0)
        endif
       endif
       
       if(idp.eq.1.or.idp.eq.3.or.idp.eq.5)then
        rp1=(rq(icdp,icz)-alfp*dlog(xph))*4.d0*.0389d0
       elseif(idp.eq.2.or.idp.eq.6)then
        rp1=rq(icdp,icz)*4.d0*.0389d0
       elseif(idp.eq.7.or.idp.eq.8)then
        rp1=alfp*dlog(xpmax/xph)*4.d0*.0389d0
       else
        stop'qglchard: wrong idp!'
       endif
       if(idp.eq.1.or.idp.eq.2.or.idp.eq.7)then
        rp2=(rq(icdt,2)-alfp*dlog(xmh))*4.d0*.0389d0
       elseif(idp.eq.3.or.idp.eq.8)then
        rp2=rq(icdt,2)*4.d0*.0389d0
       elseif(idp.eq.5.or.idp.eq.6)then
        rp2=alfp*dlog(xmmax/xmh)*4.d0*.0389d0
       else
        stop'qglchard: wrong idp!'
       endif
       
       rp=rp1*rp2/(rp1+rp2)
       zb=qgran(b10)
       gb=gb*rp/zb/gb0
       phi=2.d0*pi*qgran(b10)
       bh=dsqrt(-rp*dlog(zb))
       bbp=(dsqrt(bb0)*rp1/(rp1+rp2)+bh*cos(phi))**2+(bh*sin(phi))**2 
       bbt=(dsqrt(bb0)*rp2/(rp1+rp2)-bh*cos(phi))**2+(bh*sin(phi))**2 
       call qgbdef(bbp,bbt,xxp0,yyp0,xxt0,yyt0,xxh,yyh
     * ,int(1.5d0+qgran(b10)))
       do iqp=1,2
       do iqt=1,2
        dvhtp(iqp,iqt)=1.d0
        dvhtt(iqp,iqt)=1.d0
       enddo
       enddo

       if(xph*sgap.lt..9d0)then
        call qgfdf(xxh,yyh,xph,vpac,vtac,vpht,vtht,vpacq,vtacq
     *  ,vvxh,vvxp,vvxt,vvxpl,vvxtl,genhpp,genhtp,ip,it)
        vvxps=1.d0-exp(-vtac(it))*(1.d0-vvxt)*(1.d0-vvxtl)
     *  *(1.d0-vvxpl)
       else
        genhpp=1.d0
        genhtp=1.d0
        vvxps=0.d0
       endif
       if(xmh*sgap.lt..9d0)then
        call qgfdf(xxh,yyh,1.d0/xmh/scm,vpac,vtac,vpht,vtht,vpacq,vtacq
     *  ,vvxh,vvxp,vvxt,vvxpl,vvxtl,genhpm,genhtm,ip,it)
        vvxts=1.d0-exp(-vpac(ip))*(1.d0-vvxp)*(1.d0-vvxpl)
     *  *(1.d0-vvxtl)
       else
        genhpm=1.d0
        genhtm=1.d0
        vvxts=0.d0
       endif
       if(idp.eq.1.or.idp.eq.3.or.idp.eq.5)then
        glu1=qgpdfbi(xph,bbp,vvxps,vvxp,icdp,icz,1,2)
        sea1=qgpdfbi(xph,bbp,vvxps,vvxp,icdp,icz,2,2)
        vvxps=1.d0-(1.d0-vvxps)*(1.d0-vvxp)
       elseif(idp.eq.2.or.idp.eq.6)then
        qv1=qgvpdf(xph,icz)*exp(-bbp/rp1)/rp1
        vvxps=1.d0-(1.d0-vvxps)*(1.d0-vvxp)
       elseif(idp.eq.7.or.idp.eq.8)then
        vvxps=1.d0-(1.d0-vvxps)*exp(-vpac(ip))*(1.d0-vvxp)
        glu1=qgloopri(xpmax/xph,bbp,vvxps,1,2)
        sea1=qgloopri(xpmax/xph,bbp,vvxps,2,2)
        bbi=bbp
       endif
       if(idp.eq.1.or.idp.eq.2.or.idp.eq.7)then
        glu2=qgpdfbi(xmh,bbt,vvxts,vvxt,icdt,2,1,2)
        sea2=qgpdfbi(xmh,bbt,vvxts,vvxt,icdt,2,2,2)
        vvxts=1.d0-(1.d0-vvxts)*(1.d0-vvxt)
       elseif(idp.eq.3.or.idp.eq.8)then
        qv2=qgvpdf(xmh,2)*exp(-bbt/rp2)/rp2
        vvxts=1.d0-(1.d0-vvxts)*(1.d0-vvxt)
       elseif(idp.eq.5.or.idp.eq.6)then
        vvxts=1.d0-(1.d0-vvxts)*exp(-vtac(it))*(1.d0-vvxt)
        glu2=qgloopri(xmmax/xmh,bbt,vvxts,1,2)
        sea2=qgloopri(xmmax/xmh,bbt,vvxts,2,2)
        bbi=bbt
       endif

       sjqq=qgjit(qt0,qt0,zh*wp0*wm0,2,2)        !inclusive qq cross-section
       sjqg=qgjit(qt0,qt0,zh*wp0*wm0,1,2)        !inclusive qg cross-section
       sjgg=qgjit(qt0,qt0,zh*wp0*wm0,1,1)        !inclusive gg cross-section
       sjqqq=qgjit(qt0,qt0,zh*wp0*wm0,3,2)
       sjqqa=qgjit(qt0,qt0,zh*wp0*wm0,4,2)
       
       if(idp.eq.1.or.idp.eq.3)then
        do iq1=1,2
        do iq2=1,2
         dvhtp(iq1,iq2)=qgfhti(xmh*wp0*wm0,xph,bbp,vvxps,genhpm
     *   ,iq1,iq2,icdp,icz,2)
        enddo
        enddo
       elseif(idp.eq.7.or.idp.eq.8)then
        bbg=(xa(ip,1)+bcoll-xxh)**2+(xa(ip,2)-yyh)**2
        xg=qt0/xmh/wp0/wm0
        pdfg=qgpdfbi(xg,bbg,0.d0,0.d0,icdp,icz,1,1)
        facht=htfac*pi**3*2.d0*pdfg*genhpm
        do iq1=1,2
        do iq2=1,2
         dvhtp(iq1,iq2)=qgfhtmi(zh*wp0*wm0,xph/xpmax,facht,vvxps
     *   ,iq1,iq2,2)
        enddo
        enddo
       endif

       if(idp.eq.1.or.idp.eq.2)then
        do iq1=1,2
        do iq2=1,2
         dvhtt(iq1,iq2)=qgfhti(xph*wp0*wm0,xmh,bbt,vvxts,genhtp
     *   ,iq1,iq2,icdt,2,2)
        enddo
        enddo
       elseif(idp.eq.5.or.idp.eq.6)then
        bbg=(xb(it,1)-xxh)**2+(xb(it,2)-yyh)**2
        xg=qt0/xph/wp0/wm0
        pdfg=qgpdfbi(xg,bbg,0.d0,0.d0,icdt,2,1,1)
        facht=htfac*pi**3*2.d0*pdfg*genhtp
        do iq1=1,2
        do iq2=1,2
         dvhtt(iq1,iq2)=qgfhtmi(zh*wp0*wm0,xmh/xmmax,facht,vvxts
     *   ,iq1,iq2,2)
        enddo
        enddo
       endif
       
       if(idp.eq.1)then
        gb=gb*(glu1*glu2*sjgg*dvhtp(1,1)*dvhtt(1,1)
     *  +sea1*sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0)*dvhtp(2,2)
     *  *dvhtt(2,2)+sea1*glu2*sjqg*dvhtp(2,1)*dvhtt(1,2)
     *  +glu1*sea2*sjqg*dvhtp(1,2)*dvhtt(2,1))
       elseif(idp.eq.2.or.idp.eq.6)then
        gb=gb*qv1*(glu2*sjqg*dvhtt(1,2)
     *  +sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0)*dvhtt(2,2))
       elseif(idp.eq.3.or.idp.eq.8)then
        gb=gb*qv2*(glu1*sjqg*dvhtp(1,2)
     *  +sea1*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0)*dvhtp(2,2))
       elseif(idp.eq.5)then
        gb=gb*(glu1*glu2*sjgg*dvhtt(1,1)
     *  +sea1*sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0)*dvhtt(2,2)
     *  +sea1*glu2*sjqg*dvhtt(1,2)+glu1*sea2*sjqg*dvhtt(2,1))
       elseif(idp.eq.7)then
        gb=gb*(glu1*glu2*sjgg*dvhtp(1,1)
     *  +sea1*sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0)*dvhtp(2,2)
     *  +sea1*glu2*sjqg*dvhtp(2,1)+glu1*sea2*sjqg*dvhtp(1,2))
       endif
       gb=gb/2.d0
       gbmax=max(gb,gbmax)
       if(qgran(b10).gt.gb)then
        nret=nret+1
        if(nret.gt.30)then
         if(gbmax.le.0.d0.or.nren.gt.10)then
          iret=1
          return
         endif
         nren=nren+1
         gb0=gb0*gbmax/nren*1.5
         gbmax=0.d0
         nret=0
        endif
        goto 2
       endif

       if(idp.eq.2.or.idp.eq.6)then
        qvp1=qggrv(xph,qt0,icz,1)
        qvp2=qggrv(xph,qt0,icz,2)
        wch(1)=qvp1*glu2*sjqg*dvhtt(1,2)
        wch(2)=qvp2*glu2*sjqg*dvhtt(1,2)
        wch(3)=qvp1*sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0)*dvhtt(2,2)
        tot=wch(1)+wch(2)+wch(3)
     *  +qvp2*sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0)*dvhtt(2,2)
       elseif(idp.eq.3.or.idp.eq.8)then
        qvt1=qggrv(xmh,qt0,2,1)
        qvt2=qggrv(xmh,qt0,2,2)
        wch(1)=glu1*sjqg*qvt1*dvhtp(1,2)
        wch(2)=sea1*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0)*qvt1*dvhtp(2,2)
        wch(3)=glu1*sjqg*qvt2*dvhtp(1,2)
        tot=wch(1)+wch(2)+wch(3)
     *  +sea1*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0)*qvt2*dvhtp(2,2)
       elseif(idp.eq.1)then
        wch(1)=glu1*glu2*sjgg*dvhtp(1,1)*dvhtt(1,1)
        wch(2)=sea1*glu2*sjqg*dvhtp(2,1)*dvhtt(1,2)
        wch(3)=glu1*sea2*sjqg*dvhtp(1,2)*dvhtt(2,1)
        tot=wch(1)+wch(2)+wch(3)+sea1*sea2
     *  *(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0)*dvhtp(2,2)*dvhtt(2,2)
       elseif(idp.eq.5)then
        wch(1)=glu1*glu2*sjgg*dvhtt(1,1)
        wch(2)=sea1*glu2*sjqg*dvhtt(1,2)
        wch(3)=glu1*sea2*sjqg*dvhtt(2,1)
        tot=wch(1)+wch(2)+wch(3)
     *  +sea1*sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0)*dvhtt(2,2)
       elseif(idp.eq.7)then
        wch(1)=glu1*glu2*sjgg*dvhtp(1,1)
        wch(2)=sea1*glu2*sjqg*dvhtp(2,1)
        wch(3)=glu1*sea2*sjqg*dvhtp(1,2)
        tot=wch(1)+wch(2)+wch(3)
     *  +sea1*sea2*(sjqq/1.5d0+sjqqq/6.d0+sjqqa/6.d0)*dvhtp(2,2)
       endif
       do i=1,3
        wch(i)=wch(i)/tot
       enddo
      endif
      return
      end
    
c=============================================================================
      double precision function qglcgen(arest,lva,icz)
c------------------------------------------------------------------------
c qglcgen - light cone momentum for a string end
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr28/ arr(5),alpq
      common /qgarr43/ moniou
      common /qgdebug/ debug
      
      if(lva.eq.0)then
       ahl=1.d0-arr(1)-arr(icz)
      else
       ahl=-arr(icz)
      endif
      ahl=ahl+arest
      qglcgen=qgslc(alpq,ahl)
      return
      end

c=============================================================================
      double precision function qgslc(alpq,ahl)
c------------------------------------------------------------------------
c qgslc - technical procedure for light cone momentum sharing
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran
      
      if(ahl.lt.0.d0)then
1      z=1.d0-qgran(b10)**(1.d0/(1.d0+ahl))
       z=min(z,.9999999d0)
       x=z**(1.d0/(1.d0-alpq))
       gb=((1.d0-x)/(1.d0-z))**ahl
       if(qgran(b10).gt.gb)goto 1
      else
2      x=qgran(b10)**(1.d0/(1.d0-alpq))
       gb=(1.d0-x)**ahl
       if(qgran(b10).gt.gb)goto 2
      endif
      qgslc=x
      return
      end

c=============================================================================
      subroutine qgptstr(pts,ptr,nstr,icz)
c------------------------------------------------------------------------
c qgptstr - generation of string end transverse momenta
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr11/ b10
      common /qgarr32/ pt2str,pt2rem(3)
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran
      
      betk=dsqrt(pt2str*(pt2rem(icz)+nstr*pt2str)
     */(pt2rem(icz)+(nstr+1)*pt2str))
     
      pts=betk*(qgran(b10)+qgran(b10)+qgran(b10)-1.5d0)
     *-ptr*pt2str/(pt2rem(icz)+(nstr+1)*pt2str)
      ptr=ptr+pts
      return
      end

c=============================================================================
      subroutine qgstr(wp1,wp2,wm1,wm2,wpr,wmr
     *,ptps,ptts,ptpr,pttr,iqsp,iqst)
c-----------------------------------------------------------------------------
c qgstr - string fragmentation process
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wps(2),wms(2),icp(2),icm(2),ep1(4),ep2(4)
     *,ptps(2,2),ptts(2,2),ptp(2,2),ptt(2,2),ptpr(2),pttr(2)
      common /qgarr8/  pt2w,be(4),dc(5),deta,drho,almpt,ptdif
      common /qgarr10/ am(6)
      common /qgarr11/ b10
      common /qgarr28/ arr(5),alpq
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      if(debug.ge.2)write (moniou,201)wp1,wp2,wm1,wm2,wpr,wmr

      nrejmax=100
      smin=4.d0*am(1)**2
      if(iqsp.ne.0.and.iqst.ne.0.and.wp1*wm1-(ptps(1,1)+ptts(1,1))**2
     *-(ptps(1,2)+ptts(1,2))**2.le.smin.and.wp2*wm2
     *-(ptps(2,1)+ptts(2,1))**2-(ptps(2,2)+ptts(2,2))**2.le.smin
     *.or.iqsp.eq.0.and.iqst.ne.0.and.wp1*wm1-ptts(1,1)**2-ptts(1,2)**2
     *.le.smin.and.wp1*wm2-ptts(2,1)**2-ptts(2,2)**2.le.smin
     *.or.iqsp.ne.0.and.iqst.eq.0.and.wp1*wm1-ptps(1,1)**2-ptps(1,2)**2
     *.le.smin.and.wp2*wm1-ptps(2,1)**2-ptps(2,2)**2.le.smin
     *.or.iqsp.eq.0.and.iqst.eq.0.and.wp1*wm1.le.smin)goto 2

      nrej=0
1     nrej=nrej+1
      if(nrej.gt.nrejmax)goto 2 

      icp(1)=int(1.5+qgran(b10))             !sampling string end types
      icp(2)=-icp(1)
      icm(2)=int(1.5+qgran(b10))
      icm(1)=-icm(2)
      if(iqsp.ne.0.and.iqst.ne.0)then
       if(wp1*wm1-(ptps(1,1)+ptts(1,1))**2-(ptps(1,2)+ptts(1,2))**2
     * .gt.smin.and.wp2*wm2-(ptps(2,1)+ptts(2,1))**2
     * -(ptps(2,2)+ptts(2,2))**2.gt.smin)then
        nstr=2
        wps(1)=wp1
        wps(2)=wp2
        wms(1)=wm1
        wms(2)=wm2
        do j=1,2
        do i=1,2
         ptp(j,i)=ptps(j,i)
         ptt(j,i)=ptts(j,i)
        enddo
        enddo
       elseif(wp1*wm1-(ptps(1,1)+ptts(1,1))**2-(ptps(1,2)+ptts(1,2))**2
     * .gt.smin.and.icp(1).eq.icm(2))then
        nstr=1
        wps(1)=wp1
        wms(1)=wm1
        wpr=wpr+wp2
        wmr=wmr+wm2
        do i=1,2
         ptp(1,i)=ptps(1,i)
         ptt(1,i)=ptts(1,i)
         ptpr(i)=ptpr(i)+ptps(2,i)
         pttr(i)=pttr(i)+ptts(2,i)
        enddo
       elseif(wp2*wm2-(ptps(2,1)+ptts(2,1))**2-(ptps(2,2)+ptts(2,2))**2
     * .gt.smin.and.icp(1).eq.icm(2))then
        nstr=1
        wps(1)=wp2
        wms(1)=wm2
        wpr=wpr+wp1
        wmr=wmr+wm1
        do i=1,2
         ptp(1,i)=ptps(2,i)
         ptt(1,i)=ptts(2,i)
         ptpr(i)=ptpr(i)+ptps(1,i)
         pttr(i)=pttr(i)+ptts(1,i)
        enddo
        icp(1)=icp(2)
        icm(1)=icm(2)
       else
        goto 1
       endif

      elseif(iqsp.ne.0.and.iqst.eq.0)then
       x=qgslc(alpq,-alpq)
       if(wp1*wm1*x-ptps(1,1)**2-ptps(1,2)**2.gt.smin
     * .and.wp2*wm1*(1.d0-x)-ptps(2,1)**2-ptps(2,2)**2.gt.smin)then
        nstr=2
        wps(1)=wp1
        wps(2)=wp2
        wms(1)=wm1*x
        wms(2)=wm1*(1.d0-x)
        do j=1,2
        do i=1,2
         ptp(j,i)=ptps(j,i)
         ptt(j,i)=0.d0
        enddo
        enddo
       elseif(wp1*wm1-ptps(1,1)**2-ptps(1,2)**2.gt.smin
     * .and.icp(1).eq.icm(2))then
        nstr=1
        wps(1)=wp1
        wms(1)=wm1
        wpr=wpr+wp2
        do i=1,2
         ptp(1,i)=ptps(1,i)
         ptt(1,i)=0.d0
         ptpr(i)=ptpr(i)+ptps(2,i)
        enddo
       elseif(wp2*wm1-ptps(2,1)**2-ptps(2,2)**2.gt.smin
     * .and.icp(1).eq.icm(2))then
        nstr=1
        wps(1)=wp2
        wms(1)=wm1
        wpr=wpr+wp1
        do i=1,2
         ptp(1,i)=ptps(2,i)
         ptt(1,i)=0.d0
         ptpr(i)=ptpr(i)+ptps(1,i)
        enddo
        icp(1)=icp(2)
        icm(1)=icm(2)
       else
        goto 1
       endif

      elseif(iqsp.eq.0.and.iqst.ne.0)then
       x=qgslc(alpq,-alpq)
       if(wp1*wm1*x-ptts(1,1)**2-ptts(1,2)**2.gt.smin
     * .and.wp1*wm2*(1.d0-x)-ptts(2,1)**2-ptts(2,2)**2.gt.smin)then
        nstr=2
        wps(1)=wp1*x
        wps(2)=wp1*(1.d0-x)
        wms(1)=wm1
        wms(2)=wm2
        do j=1,2
        do i=1,2
         ptp(j,i)=0.d0
         ptt(j,i)=ptts(j,i)
        enddo
        enddo
       elseif(wp1*wm1-ptts(1,1)**2-ptts(1,2)**2.gt.smin
     * .and.icp(1).eq.icm(2))then
        nstr=1
        wps(1)=wp1
        wms(1)=wm1
        wmr=wmr+wm2
        do i=1,2
         ptp(1,i)=0.d0
         ptt(1,i)=ptts(1,i)
         pttr(i)=pttr(i)+ptts(2,i)
        enddo
       elseif(wp1*wm2-ptts(2,1)**2-ptts(2,2)**2.gt.smin
     * .and.icp(1).eq.icm(2))then
        nstr=1
        wps(1)=wp1
        wms(1)=wm2
        wmr=wmr+wm1
        do i=1,2
         ptp(1,i)=0.d0
         ptt(1,i)=ptts(2,i)
         pttr(i)=pttr(i)+ptts(1,i)
        enddo
        icp(1)=icp(2)
        icm(1)=icm(2)
       else
        goto 1
       endif
       
      else
       xp=qgslc(alpq,-alpq)
       xm=qgslc(alpq,-alpq)
       if(wp1*wm1*xp*xm.gt.smin
     * .and.wp1*wm1*(1.d0-xp)*(1.d0-xm).gt.smin)then
        nstr=2
        wps(1)=wp1*xp
        wps(2)=wp1*(1.d0-xp)
        wms(1)=wm1*xm
        wms(2)=wm1*(1.d0-xm)
        do j=1,2
        do i=1,2
         ptp(j,i)=0.d0
         ptt(j,i)=0.d0
        enddo
        enddo
       elseif(wp1*wm1.gt.smin.and.icp(1).eq.icm(2))then
        nstr=1
        wps(1)=wp1
        wms(1)=wm1
        do i=1,2
         ptp(1,i)=0.d0
         ptt(1,i)=0.d0
        enddo
        if(qgran(b10).gt..5d0)then
         icp(1)=-icp(1)
         icm(1)=-icm(1)
        endif
       else
        goto 1
       endif
      endif

      do i=1,nstr
       ww=wps(i)*wms(i)
       pt2p=ptp(i,1)**2+ptp(i,2)**2
       pt2t=ptt(i,1)**2+ptt(i,2)**2
       zz=qgtwd(ww,pt2p,pt2t)
       wpp=zz*wps(i)
       wpm=pt2p/wpp
       ep1(1)=.5d0*(wpp+wpm)
       ep1(2)=.5d0*(wpp-wpm)
       ep1(3)=ptp(i,1)
       ep1(4)=ptp(i,2)
       ep2(1)=.5d0*(wps(i)+wms(i))-ep1(1)
       ep2(2)=.5d0*(wps(i)-wms(i))-ep1(2)
       ep2(3)=ptt(i,1)
       ep2(4)=ptt(i,2)
       call qggene(ep1,ep2,icp(i),icm(i),0)      !string fragmentation
      enddo
      if(debug.ge.3)write (moniou,202)wpr,wmr
      return
      
2     wpr=wpr+wp1+wp2
      wmr=wmr+wm1+wm2
      do i=1,2
       ptpr(i)=ptpr(i)+ptps(1,i)+ptps(2,i)
       pttr(i)=pttr(i)+ptts(1,i)+ptts(2,i)
      enddo
      if(debug.ge.3)write (moniou,202)wpr,wmr
      
201   format(2x,'qgstr: wp1=',e10.3,2x,'wp2=',e10.3,2x,'wm1=',e10.3
     *,2x,'wm2=',e10.3,2x,'wpr=',e10.3,2x,'wmr=',e10.3)
202   format(2x,'qgstr - returned light cone momenta:'
     *,2x,'wpr=',e10.3,2x,'wmr=',e10.3)
      return
      end
      
c=============================================================================
      subroutine qgvstr(wp1,wp2,wm1,wm2,ptps,ptts
     *,izp,izt,iqvp,iqvt,lqa,lqb,idpi,icz,iret)
c-----------------------------------------------------------------------------
c qgvstr - string fragmentation process (involving valence diquarks)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wps(2),wms(2),icps(2),icms(2),ep1(4),ep2(4)
     *,ptps(2,2),ptts(2,2),ptpf(2,2),pttf(2,2)
      common /qgarr8/  pt2w,be(4),dc(5),deta,drho,almpt,ptdif
      common /qgarr10/ am(6)
      common /qgarr11/ b10
      common /qgarr19/ wpiex(2,3)
      common /qgarr28/ arr(5),alpq
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /jdiff/   jdiff             !diffr. type (external use)
      external qgran

      if(debug.ge.2)write (moniou,201)wp1,wp2,wm1,wm2

      nrejmax=100
      jdiff0=jdiff
      if(iqvp.eq.0.and.iqvt.eq.0)stop'no valence?!'

      nrej=0
1     nrej=nrej+1
      jdiff=jdiff0
      jexp=0
      jext=0
      if(nrej.gt.nrejmax)goto 3
      am1=am(1)
      am2=am(1)
      if(iqvp.eq.0)then
       am12=am(1)
      else
       am12=am(icz)
      endif
      if(iqvt.eq.0)then
       am21=am(1)
      else
       am21=am(2)
      endif
      if((wp1+wp2)*(wm1+wm2)-(ptps(1,1)+ptps(2,1)+ptts(1,1)
     *+ptts(2,1))**2-(ptps(1,2)+ptps(2,2)+ptts(1,2)+ptts(2,2))**2
     *.le.(am1+am12+am21+am2)**2)goto 3
      
      if(lqa.eq.1.and.iqvp.eq.1
     *.and.qgran(b10).lt.wpiex(1,icz)+wpiex(2,icz))then
       if(qgran(b10).lt.wpiex(1,icz)/(wpiex(1,icz)+wpiex(2,icz)))then
        jexp=1
        call qgpiex(xnp,ptp,am2hp,izp,icpr,icps(1),icps(2),icz)
       else
        jexp=2
        call qgrex(xnp,ptp,am2hp,izp,icpr,icps(1),icps(2),icz,jexp)
       endif
       xp=qgslc(arr(1),-arr(1))
       am12=am(1)
       if(icps(1).lt.0.d0)then
        xp=1.d0-xp
        icd=icps(1)
        icps(1)=icps(2)
        icps(2)=icd
        amd=am1
        am1=am12
        am12=amd
       endif
       call qgcs(cpr,spr)
       pt2pr=(ptps(1,1)+ptps(2,1)+ptp*cpr)**2
     * +(ptps(1,2)+ptps(2,2)+ptp*spr)**2
       wprp=xnp*wp1
       wprm=(am2hp+pt2pr)/wprp
       jdiff=10+jexp                              !pion exch. for forward prod.

       if(lqb.eq.1.and.iqvt.eq.1
     * .and.qgran(b10).lt.wpiex(1,2)+wpiex(2,2))then
        if(qgran(b10).lt.wpiex(1,2)/(wpiex(1,2)+wpiex(2,2)))then     
         jext=1
         call qgpiex(xnt,ptt,am2ht,izt,ictg,icms(2),icms(1),2)
        else
         jext=2
         call qgrex(xnt,ptt,am2ht,izt,ictg,icms(2),icms(1),2,jext)
        endif
        xm=qgslc(arr(1),-arr(1))
         am21=am(1)
        call qgcs(ctg,stg)
        pt2tg=(ptts(1,1)+ptts(2,1)+ptt*ctg)**2
     *  +(ptts(1,2)+ptts(2,2)+ptt*stg)**2
        wtgm=xnt*wm1
        wtgp=(am2ht+pt2tg)/wtgm

        wps(1)=xp*(wp1-wprp-wtgp)
        wps(2)=wp1-wprp-wtgp-wps(1)
        wms(1)=xm*(wm1-wtgm-wprm)
        wms(2)=wm1-wtgm-wprm-wms(1)
        ptpf(1,1)=-ptp*cpr*xp
        ptpf(1,2)=-ptp*spr*xp
        ptpf(2,1)=-ptp*cpr*(1.d0-xp)
        ptpf(2,2)=-ptp*spr*(1.d0-xp)
        pttf(1,1)=-ptt*ctg*xm
        pttf(1,2)=-ptt*stg*xm
        pttf(2,1)=-ptt*ctg*(1.d0-xm)
        pttf(2,2)=-ptt*stg*(1.d0-xm)

        if(wps(1).gt.0.d0.and.wps(2).gt.0.d0.and.wps(1)*wms(1)
     *  -(ptpf(1,1)+pttf(1,1))**2-(ptpf(1,2)+pttf(1,2))**2
     *  .gt.(am1+am21)**2.and.wps(2)*wms(2)-(ptpf(2,1)+pttf(2,1))**2
     *  -(ptpf(2,2)+pttf(2,2))**2.gt.(am12+am2)**2)then
         nstr=2
        else
         goto 1
        endif
        ep1(1)=.5d0*(wprp+wprm)
        ep1(2)=.5d0*(wprp-wprm)
        ep1(3)=ptps(1,1)+ptps(2,1)+ptp*cpr
        ep1(4)=ptps(1,2)+ptps(2,2)+ptp*spr
        call qgreg(ep1,icpr)
        ep2(1)=.5d0*(wtgp+wtgm)
        ep2(2)=.5d0*(wtgp-wtgm)
        ep2(3)=ptts(1,1)+ptts(2,1)+ptt*ctg
        ep2(4)=ptts(1,2)+ptts(2,2)+ptt*stg
        call qgreg(ep2,ictg)
        
       else
        if(iqvt.eq.1)then
         call qgvlc(wm1-wprm,wms(2),wms(1),izt,icms(2),icms(1))
        else
         icms(2)=int(1.5+qgran(b10))             !sampling string end types
         icms(1)=-icms(2)
         if(idpi.eq.0)then
          wms(1)=wm1*(wm1+wm2-wprm)/(wm1+wm2)
          wms(2)=wm2*(wm1+wm2-wprm)/(wm1+wm2)
         else
          wms(1)=(wm1-wprm)*qgslc(alpq,-alpq)
          wms(2)=wm1-wprm-wms(1)
         endif
        endif
        wps(1)=xp*(1.d0-xnp)*wp1
        wps(2)=(1.d0-xp)*(1.d0-xnp)*wp1
        do j=1,2
        do i=1,2
         pttf(j,i)=ptts(3-j,i)
        enddo
        enddo
        ptpf(1,1)=-ptp*cpr*xp
        ptpf(1,2)=-ptp*spr*xp
        ptpf(2,1)=-ptp*cpr*(1.d0-xp)
        ptpf(2,2)=-ptp*spr*(1.d0-xp)
        
        if(wps(1)*wms(1)-(ptpf(1,1)+pttf(1,1))**2-(ptpf(1,2)
     *  +pttf(1,2))**2.gt.(am1+am21)**2.and.wps(2)*wms(2)-(ptpf(2,1)
     *  +pttf(2,1))**2-(ptpf(2,2)+pttf(2,2))**2.gt.(am12+am2)**2)then
         nstr=2
        else
         goto 1
        endif
        ep1(1)=.5d0*(wprp+wprm)
        ep1(2)=.5d0*(wprp-wprm)
        ep1(3)=ptps(1,1)+ptps(2,1)+ptp*cpr
        ep1(4)=ptps(1,2)+ptps(2,2)+ptp*spr
        call qgreg(ep1,icpr)
       endif
       goto 2
      endif

      do j=1,2
      do i=1,2
       ptpf(j,i)=ptps(j,i)
      enddo
      enddo
      if(lqb.eq.1.and.iqvt.eq.1
     *.and.qgran(b10).lt.wpiex(1,2)+wpiex(2,2))then
       if(qgran(b10).lt.wpiex(1,2)/(wpiex(1,2)+wpiex(2,2)))then     
        jext=1
        call qgpiex(xnt,ptt,am2ht,izt,ictg,icms(2),icms(1),2)
       else
        jext=2
        call qgrex(xnt,ptt,am2ht,izt,ictg,icms(2),icms(1),2,jext)
       endif
       xm=qgslc(arr(1),-arr(1))
       am21=am(1)
       call qgcs(ctg,stg)
       pt2tg=(ptts(1,1)+ptts(2,1)+ptt*ctg)**2
     * +(ptts(1,2)+ptts(2,2)+ptt*stg)**2
       wtgm=xnt*wm1
       wtgp=(am2ht+pt2tg)/wtgm

       if(iqvp.eq.1)then
        call qgvlc(wp1-wtgp,wps(1),wps(2),izp,icps(1),icps(2))
        if(icps(1).lt.0)then
         wpd=wps(1)
         wps(1)=wps(2)
         wps(2)=wpd
         icd=icps(1)
         icps(1)=icps(2)
         icps(2)=icd
         amd=am1
         am1=am12
         am12=amd
         do j=1,2
         do i=1,2
          ptpf(j,i)=ptps(3-j,i)
         enddo
         enddo
        endif
       else
        icps(1)=int(1.5+qgran(b10))              !sampling string end types
        icps(2)=-icps(1)
        if(idpi.eq.0)then
         wps(1)=wp1*(wp1+wp2-wtgp)/(wp1+wp2)
         wps(2)=wp2*(wp1+wp2-wtgp)/(wp1+wp2)
        else
         wps(1)=(wp1-wtgp)*qgslc(alpq,-alpq)
         wps(2)=wp1-wtgp-wps(1)
        endif
       endif
       wms(1)=xm*(1.d0-xnt)*wm1
       wms(2)=(1.d0-xm)*(1.d0-xnt)*wm1
       pttf(1,1)=-ptt*ctg*xm
       pttf(1,2)=-ptt*stg*xm
       pttf(2,1)=-ptt*ctg*(1.d0-xm)
       pttf(2,2)=-ptt*stg*(1.d0-xm)
        
       if(wps(1)*wms(1)-(ptpf(1,1)+pttf(1,1))**2-(ptpf(1,2)
     * +pttf(1,2))**2.gt.(am21+am1)**2.and.wps(2)*wms(2)-(ptpf(2,1)
     * +pttf(2,1))**2-(ptpf(2,2)+pttf(2,2))**2.gt.(am2+am12)**2)then
        nstr=2
       else
        goto 1
       endif
       ep2(1)=.5d0*(wtgp+wtgm)
       ep2(2)=.5d0*(wtgp-wtgm)
       ep2(3)=ptts(1,1)+ptts(2,1)+ptt*ctg
       ep2(4)=ptts(1,2)+ptts(2,2)+ptt*stg
       call qgreg(ep2,ictg)
       goto 2
      endif

      if(iqvp.eq.1)then
       call qgvlc(wp1,wps(1),wps(2),izp,icps(1),icps(2))
       if(icps(1).lt.0)then
        wpd=wps(1)
        wps(1)=wps(2)
        wps(2)=wpd
        icd=icps(1)
        icps(1)=icps(2)
        icps(2)=icd
        amd=am1
        am1=am12
        am12=amd
        do j=1,2
        do i=1,2
         ptpf(j,i)=ptps(3-j,i)
        enddo
        enddo
       endif
      else
       icps(1)=int(1.5+qgran(b10))               !sampling string end types
       icps(2)=-icps(1)
       if(idpi.eq.0)then
        wps(1)=wp1
        wps(2)=wp2
       else
        wps(1)=wp1*qgslc(alpq,-alpq)
        wps(2)=wp1-wps(1)
       endif
      endif
      if(iqvt.eq.1)then
       call qgvlc(wm1,wms(2),wms(1),izt,icms(2),icms(1))
      else
       icms(2)=int(1.5+qgran(b10))               !sampling string end types
       icms(1)=-icms(2)
       if(idpi.eq.0)then
        wms(1)=wm1
        wms(2)=wm2
       else
        wms(1)=wm1*qgslc(alpq,-alpq)
        wms(2)=wm1-wms(1)
       endif
      endif
      do j=1,2
      do i=1,2
       pttf(j,i)=ptts(3-j,i)
      enddo
      enddo
      
      do i=1,2
       if(iabs(icps(i)).gt.5.and.iabs(icms(i)).gt.5)goto 1
      enddo
      
      if(wps(1)*wms(1)-(ptpf(1,1)+pttf(1,1))**2-(ptpf(1,2)
     *+pttf(1,2))**2.gt.(am1+am21)**2.and.wps(2)*wms(2)-(ptpf(2,1)
     *+pttf(2,1))**2-(ptpf(2,2)+pttf(2,2))**2.gt.(am12+am2)**2)then
       nstr=2
      elseif(idpi.ne.0.and.iqvp.eq.1.and.(icps(1)+icms(1).eq.0
     *.and.wp1*wm1-(ptpf(1,1)+ptpf(2,1)+pttf(1,1)+pttf(2,1))**2
     *-(ptpf(1,2)+ptpf(2,2)+pttf(1,2)+pttf(2,2))**2.gt.(am12+am2)**2
     *.or.icps(2)+icms(2).eq.0.and.wp1*wm1-(ptpf(1,1)+ptpf(2,1)
     *+pttf(1,1)+pttf(2,1))**2-(ptpf(1,2)+ptpf(2,2)+pttf(1,2)
     *+pttf(2,2))**2.gt.(am1+am21)**2))then
       nstr=1
       wps(1)=wp1
       wms(1)=wm1
       if(icps(1)+icms(1).eq.0)then
        icps(1)=icps(2)
        icms(1)=icms(2)
       endif
       do i=1,2
        ptpf(1,i)=ptpf(1,i)+ptpf(2,i)
        pttf(1,i)=pttf(1,i)+pttf(2,i)
       enddo
      elseif(idpi.ne.0.and.iqvt.eq.1.and.icps(2)+icms(2).eq.0.and.wp1
     **wm1-(ptpf(1,1)+ptpf(2,1)+pttf(1,1)+pttf(2,1))**2-(ptpf(1,2)
     *+ptpf(2,2)+pttf(1,2)+pttf(2,2))**2.gt.(am1+am21)**2)then
       nstr=1
       wps(1)=wp1
       wms(1)=wm1
       do i=1,2
        ptpf(1,i)=ptpf(1,i)+ptpf(2,i)
        pttf(1,i)=pttf(1,i)+pttf(2,i)
       enddo
      else
       goto 1
      endif

2     do i=1,nstr
       ww=wps(i)*wms(i)
       wwi=ww-(ptpf(i,1)+pttf(i,1))**2-(ptpf(i,2)+pttf(i,2))**2
       pt2p=ptpf(i,1)**2+ptpf(i,2)**2
       pt2t=pttf(i,1)**2+pttf(i,2)**2
       zz=qgtwd(ww,pt2p,pt2t)
       wpp=zz*wps(i)
       wpm=pt2p/wpp
       ep1(1)=.5d0*(wpp+wpm)
       ep1(2)=.5d0*(wpp-wpm)
       ep1(3)=ptpf(i,1)
       ep1(4)=ptpf(i,2)
       ep2(1)=.5d0*(wps(i)+wms(i))-ep1(1)
       ep2(2)=.5d0*(wps(i)-wms(i))-ep1(2)
       ep2(3)=pttf(i,1)
       ep2(4)=pttf(i,2)
       if((iqvt.eq.1.and.jext.eq.0.or.iqvp.eq.1.and.jexp.eq.0)
     * .and.(i.eq.1.and.nstr.eq.2.and.wwi.gt.1.5d0*(am1+am21+am(1))**2
     * .or.i.eq.2.and.wwi.gt.1.5d0*(am12+am2+am(1))**2))then
        j3p=1
       else
        j3p=0
       endif
       call qggene(ep1,ep2,icps(i),icms(i),j3p)    !string fragmentation
      enddo
      iret=0
      if(debug.ge.3)write (moniou,202)iret
      return
      
3     iret=1
      if(debug.ge.3)write (moniou,202)iret
      
201   format(2x,'qgvstr: wp1=',e10.3,2x,'wp2=',e10.3,2x,'wm1=',e10.3
     *,2x,'wm2=',e10.3)
202   format(2x,'qgvstr: error code=',i4)
      return
      end
       
c===========================================================================
      subroutine qgpiex(xh,pt,am2hf,izp,icr,ic1,ic2,icz)
c---------------------------------------------------------------------------
c qgpiex - pion exchange process
c xh  - leading hadron LC momentum fraction;
c pt  - leading hadron pt;
c icr - leading hadron type;
c ic1, ic2 - (anti)quark flavours for string ends
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr10/ am(6)
      common /qgarr11/ b10
      common /qgarr21/ dmmin(3),dmres(3),wdres(3)
      common /qgarr40/ apipr,alphaf,bpi,breg
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran
      
      am2=am(1)**2
      am2h=am(icz)**2
      if(icz.eq.2)then
       am2hf=am2h
      else
       am2hf=dmmin(icz)**2
      endif
      is=iabs(izp)/izp
      if(icz.eq.1)then
       if(qgran(b10).lt..5d0)then
        icr=17          !rho0
        ic1=(3-izp)/2
        ic2=ic1-3
       else
        icr=16*is       !rho+- (+-16)
        ic1=int(1.5d0+qgran(b10))
        ic2=-ic1
       endif
      elseif(icz.eq.2)then
       aks=9.d0*qgran(b10)
       if(aks.lt.3.d0)then
        icr=izp
        ic1=izp-is
        if(qgran(b10).lt..3333d0)ic1=3*is-ic1
        ic2=-ic1
       elseif(aks.lt.7.d0)then
        icr=5*is-izp
        ic1=izp-is
        ic2=ic1-3*is
       else
        am2hf=dmmin(2)**2
        icr=izp+5*is
        ic1=4*is-izp
        ic2=ic1-3*is
       endif
       if(is.eq.-1)then
        ics=ic1
        ic1=ic2
        ic2=ics
       endif
      elseif(icz.eq.3)then
       if(qgran(b10).lt..3333d0)then
        icr=14*is+izp   !K* (+-18,19)
        ic1=is*(izp-3*is)
        ic2=-ic1
       else
        icr=23*is-izp   !K*
        if(izp.eq.4.or.izp.eq.-5)then
         ic1=1
        else
         ic1=2
        endif
        ic2=ic1-3
       endif
      else
       stop'qgpiex: icz???'
      endif

      alf=2.d0+dels+2.d0*apipr*am2
1     xh=1.d0-qgran(b10)**(1.d0/alf)
      t=-dlog(qgran(b10))/bpi/2.d0
      pt2=xh*t-(1.d0-xh)*(am2hf-xh*am2h)
      if(pt2.lt.0.d0)goto 1
      gb=(1.d0-xh)**(2.d0*apipr*t)*t*am2/(t+am2)**2*4.d0
      if(gb.gt.1.d0)stop'qgpiex: gb>1?!'
      if(qgran(b10).gt.gb)goto 1
      pt=dsqrt(pt2)
      return
      end
      
c===========================================================================
      subroutine qgrex(xh,pt,am2hf,izp,icr,ic1,ic2,icz,jex)
c---------------------------------------------------------------------------
c qgrex - reggeon exchange process
c xh  - leading hadron LC momentum fraction;
c pt  - leading hadron pt;
c icr - leading hadron type;
c ic1, ic2 - (anti)quark flavours for string ends;
c icz - incident hadron class;
c jex - reggeon type (2 - f)
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr10/ am(6)
      common /qgarr11/ b10
      common /qgarr21/ dmmin(3),dmres(3),wdres(3)
      common /qgarr28/ arr(5),alpq
      common /qgarr40/ apipr,alphaf,bpi,breg
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran
      
      am2h=am(icz)**2
      if(jex.eq.2)then                  !f-reggeon
       am2hf=am2h
       icr=izp                          !unchanged
       ic1=int(1.5d0+qgran(b10))
       ic2=-ic1
      else
       stop'qgrex: jex???'
      endif
      alf=2.d0-2.d0*alphaf+dels

1     xh=1.d0-qgran(b10)**(1.d0/alf)
      t=-dlog(qgran(b10))/breg/2.d0
      pt2=xh*t-(1.d0-xh)*(am2hf-xh*am2h)
      if(pt2.lt.0.d0)goto 1
      gb=(1.d0-xh)**(2.d0*apipr*t)
      if(qgran(b10).gt.gb)goto 1
      pt=dsqrt(pt2)
      return
      end

c===========================================================================
      subroutine qgvlc(wp,wp1,wp2,izp,ic1,ic2)
c---------------------------------------------------------------------------
c qgvlc - valence quark (diquark) flavours & LC momenta for string ends
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr11/ b10
      common /qgarr28/ arr(5),alpq
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      if(debug.ge.2)write (moniou,201)wp,izp

      iza=iabs(izp)
      if(iza.ne.0)is=izp/iza
      if(iza.le.1)then
       wp1=wp*qgslc(arr(1),-arr(1))
       if(iza.eq.0)then
        ic1=int(1.5d0+qgran(b10))*(1-2*int(.5d0+qgran(b10)))
        ic2=-ic1
       else
        ic=int(1.5d0+qgran(b10))
        ic1=ic*(3-2*ic)*is
        ic2=(3-ic)*(2*ic-3)*is
       endif
      elseif(iza.le.3)then
       zq=qgslc(arr(1),-arr(2))
       wp1=wp*zq
       if(qgran(b10).lt.2.d0/(3.d0-zq))then
        ic1=(iza-1)*is
        ic2=3*is
       else
        ic1=(4-iza)*is
        ic2=(4+iza)*is
       endif
      elseif(iza.le.5)then
       wp1=wp*qgslc(arr(1),-arr(3))
       ic1=(iza-3)*is
       ic2=-4*is
      else
       stop'qgvlc: izp?'
      endif
      wp2=wp-wp1
 
      if(debug.ge.3)write (moniou,202)ic1,ic2,wp1,wp2
201   format(2x,'qgvlc: remnant LC momentum wp=',e10.3
     *,2x,' remnant type izp=',i2)
202   format(2x,'qgvlc-end: parton flavors ic1=',i2,' ic2=',i2
     *,' LC momenta: ',2e10.3)
      return
      end
      
c=============================================================================
      subroutine qgvdef(ic,ic1,ic2,izp)
c-----------------------------------------------------------------------------
c qgvdef - valence quark flavour for hard scattering
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      if(debug.ge.2)write (moniou,201)ic,izp

      iza=iabs(izp)
      if(iza.ne.0)is=izp/iza
      if(iza.eq.0)then
       ic1=ic*(1-2*int(.5d0+qgran(b10)))
       ic2=-ic1
      elseif(iza.eq.1)then
       ic1=ic*(3-2*ic)*is
       ic2=(3-ic)*(2*ic-3)*is
      elseif(iza.eq.2)then
       ic1=ic*is
       ic2=3*ic*is
      elseif(iza.eq.3)then
       ic1=(3-ic)*is
       ic2=(4*ic-1)*is
      elseif(iza.eq.4)then
       if(ic.eq.1)then
        ic1=is
        ic2=-4*is
       else
        ic1=-4*is
        ic2=is
       endif
      elseif(iza.eq.5)then
       if(ic.eq.1)then
        ic1=2*is
        ic2=-4*is
       else
        ic1=-4*is
        ic2=2*is
       endif
      else
       stop'qgvdef: izp?'
      endif

      if(debug.ge.3)write (moniou,202)ic1,ic2
201   format(2x,'qgvdef: sampled parton ic=',i2,' hadron type izp=',i1)
202   format(2x,'qgvdef-end: parton flavors ic1=',i2,'ic2=',i2)
      return
      end

c=============================================================================
      subroutine qghot(wph,wmh,wp1,wp2,wm1,wm2,wch,ptps,ptts,xxp,yyp
     *,bbi,vvxps,vvxts,idp,izp,izt,icdp,icdt,iqvp,iqvt,lqa,lqb
     *,ip,it,iret)
c---------------------------------------------------------------------------
c qghot - semihard process
c wph,wmh - LC momenta for ladder leg partons,
c wp1,wp2,wm1,wm2 - LC momenta for constituent partons,
c izp, izt  - types of proj. and targ. remnants,
c icdp,icdt - proj. and targ.  diffractive eigenstates,
c iqq - type of semihard process: 1,5,7 - gg, 2 - q_vg, 3 - gq_v, 4 - q_vq_v
c iqvp=1 - if proj. diquark involved (0 otherwise),
c iqvt=1 - if targ. diquark involved (0 otherwise)
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      character*2 tyq
      parameter(iapmax=208,njmax=50000,nfock=3)
      dimension ept(4),ep1(4),ep2(4),ep3(4),ebal(4),ebal2(4),ey(3)
     *,wsub(3),qmin(2),wp(2),iqc(2),iqp(2),nqc(2),ncc(2,2),qv1(30,50)
     *,zv1(30,50),qm1(30,50),iqv1(30,50),ldau1(30,49),lpar1(30,50)
     *,qv2(30,50),zv2(30,50),qm2(30,50),iqv2(30,50),ldau2(30,49),iqb(2)
     *,lpar2(30,50),ptps(2,2),ptts(2,2),ptpf(2,2),pttf(2,2),wch(3)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),bcoll
      common /qgarr8/  pt2w,be(4),dc(5),deta,drho,almpt,ptdif
      common /qgarr10/ am(6)
      common /qgarr11/ b10
      common /qgarr12/ nsp
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr19/ wpiex(2,3)
      common /qgarr26/ factk,fqscal
      common /qgarr28/ arr(5),alpq
      common /qgarr37/ eqj(4,njmax),iqj(njmax),ncj(2,njmax),nj
      common /qgarr38/ htfac
      common /qgarr42/ tyq(16)
      common /qgarr43/ moniou
      common /qgarr46/ iconab(iapmax,iapmax),icona(iapmax)
     *,iconb(iapmax)
      common /qgarr51/ epsxmn
      common /qgdebug/ debug
      common /jdiff/   jdiff             !diffr. type (external use)
      common /ebal/    ebal0(4),ebal1(4)
      external qgran

      if(debug.ge.1)write (moniou,201)idp,wph,wmh,izp,izt,icdp,icdt
     *,iqvp,iqvt,nj

      nj0=nj                                     !save number of final partons
      nsp0=nsp                                   !save number of final particles
      jdiff0=jdiff
      nrejmax=500
      iret=0
      nrej=0
      do i=1,4
       ebal2(i)=ebal1(i)
      enddo

1     nrej=nrej+1
      if(nrej.gt.nrejmax)goto 10
      ebal(1)=.5d0*(wph+wp1+wp2+wmh+wm1+wm2)     !balans of 4-momentum
      ebal(2)=.5d0*(wph+wp1+wp2-wmh-wm1-wm2)
      ebal(3)=ptps(1,1)+ptps(2,1)+ptts(1,1)+ptts(2,1)
      ebal(4)=ptps(1,2)+ptps(2,2)+ptts(1,2)+ptts(2,2)
      nj=nj0
      nsp=nsp0
      jdiff=jdiff0
      wpi=wph                                    !LC+ for the hard interaction
      wmi=wmh                                    !LC- for the hard interaction
      if(debug.ge.2)write (moniou,202)wpi*wmi
      s2min=4.d0*fqscal*qt0                      !threshold energy
      if(wph*wmh.lt.s2min)stop'qghot: sy<s2min!!!'
      do i=1,4
       ebal1(i)=ebal2(i)
      enddo

      aks=qgran(b10)
      do jc=1,3
       aks=aks-wch(jc)
       if(aks.lt.0.d0)goto 2
      enddo
      jc=4

c-------------------------------------------------
c color connections to constituent partons
2     if(idp.eq.2.or.idp.eq.4.or.idp.eq.6)then   !q_v proj. involved
       ic=jc
       if(ic.gt.2)ic=ic-2   
       call qgvdef(ic,ic1,ic2,izp)               !leading state flavor
       iqc(1)=ic1                                !upper leg parton
       nj=nj+1
       if(nj.gt.njmax)stop'increase njmax!!!'
       iqj(nj)=ic2                               !leading jet parton
       ncc(1,1)=nj                               !color conn. with leading jet
       ncc(2,1)=0
       wmq=(ptps(2,1)**2+ptps(2,2)**2)/wp1
       wmi=wmi-wmq
       if(wmi.lt.0.d0.or.wpi*wmi.le.s2min)goto 10
       eqj(1,nj)=.5d0*(wp1+wmq)
       eqj(2,nj)=.5d0*(wp1-wmq)
       eqj(3,nj)=ptps(2,1)
       eqj(4,nj)=ptps(2,2)
       
      else                                       !g or q_sea as upper leg parton
       do i=1,2
       do j=1,2
        ptpf(j,i)=ptps(j,i)
       enddo
       enddo
       if(idp.eq.1.or.idp.eq.3.or.idp.eq.5)then
        if(iqvp.ne.0)then
         if(lqa.ne.1.or.qgran(b10).gt.wpiex(1,icz)+wpiex(2,icz))then
          call qgvlc(wp1,wps1,wps2,izp,icp1,icp2)
          if(icp1.lt.0)then
           icd=icp1
           wpd=wps1
           icp1=icp2
           icp2=icd
           wps1=wps2
           wps2=wpd
           do i=1,2
           do j=1,2
            ptpf(j,i)=ptps(3-j,i)
           enddo
           enddo
          endif
         else
          if(qgran(b10).lt.wpiex(1,icz)
     *    /(wpiex(1,icz)+wpiex(2,icz)))then
           jex=1
           call qgpiex(xnp,ptp,am2hp,izp,icpr,icp1,icp2,icz)
          else
           jex=2
           call qgrex(xnp,ptp,am2hp,izp,icpr,icp1,icp2,icz,jex)
          endif
          xp=qgslc(arr(1),-arr(1))
          if(icp1.lt.0.d0)then
           xp=1.d0-xp
           icd=icp1
           icp1=icp2
           icp2=icd
          endif
          call qgcs(cpr,spr)
          pt2pr=(ptps(1,1)+ptps(2,1)+ptp*cpr)**2
     *    +(ptps(1,2)+ptps(2,2)+ptp*spr)**2
          wprp=xnp*wp1
          wprm=(am2hp+pt2pr)/wprp
          wmi=wmi-wprm
          if(wmi.lt.0.d0.or.wpi*wmi.le.s2min)goto 1
          ep1(1)=.5d0*(wprp+wprm)
          ep1(2)=.5d0*(wprp-wprm)
          ep1(3)=ptps(1,1)+ptps(2,1)+ptp*cpr
          ep1(4)=ptps(1,2)+ptps(2,2)+ptp*spr
          do i=1,4
           ebal(i)=ebal(i)-ep1(i)
          enddo
          call qgreg(ep1,icpr)
          wps1=xp*(wp1-wprp)
          wps2=wp1-wprp-wps1
          ptpf(1,1)=-ptp*cpr*xp
          ptpf(1,2)=-ptp*spr*xp
          ptpf(2,1)=-ptp*cpr*(1.d0-xp)
          ptpf(2,2)=-ptp*spr*(1.d0-xp)
          jdiff=10+jex                           !pion exch. for forward prod.
         endif
        else
         icp1=int(1.5d0+qgran(b10))
         icp2=-icp1
         wps1=wp1
         wps2=wp2
        endif
       elseif(idp.eq.7.or.idp.eq.8)then          !g or q_sea from 3P-vertex
        icp1=int(1.5d0+qgran(b10))
        icp2=-icp1
        wps1=wp1*qgslc(alpq,-alpq)
        wps2=wp1-wps1
       else
        stop'qghot: idp?'
       endif
       
       if(jc.eq.1.or.jc.eq.3)then                !g as upper leg parton
        iqc(1)=0                                 !upper leg parton
        nj=nj+1
        if(nj.gt.njmax)stop'increase njmax!!!'
        iqj(nj)=icp2                             !leading jet parton
        ncc(1,1)=nj                              !color conn. to leading jet
        wmq2=(ptpf(2,1)**2+ptpf(2,2)**2)/wps2
        wmq1=(ptpf(1,1)**2+ptpf(1,2)**2)/wps1
        if(wmi-wmq1-wmq2.le.0.d0.or.wpi*(wmi-wmq1-wmq2).le.s2min)then
         if(wps1.gt.wps2)then
          ptpf(1,1)=ptpf(1,1)+ptpf(2,1)
          ptpf(1,2)=ptpf(1,2)+ptpf(2,2)
          ptpf(2,1)=0.d0
          ptpf(2,2)=0.d0
          wmq1=(ptpf(1,1)**2+ptpf(1,2)**2)/wps1
          wmq2=0.d0
         else
          ptpf(2,1)=ptpf(2,1)+ptpf(1,1)
          ptpf(2,2)=ptpf(2,2)+ptpf(1,2)
          ptpf(1,1)=0.d0
          ptpf(1,2)=0.d0
          wmq2=(ptpf(2,1)**2+ptpf(2,2)**2)/wps2
          wmq1=0.d0
         endif
        endif
        wmi=wmi-wmq1-wmq2
        if(wmi.lt.0.d0.or.wpi*wmi.le.s2min)goto 1
        eqj(1,nj)=.5d0*(wps2+wmq2)
        eqj(2,nj)=.5d0*(wps2-wmq2)
        eqj(3,nj)=ptpf(2,1)
        eqj(4,nj)=ptpf(2,2)
        nj=nj+1
        if(nj.gt.njmax)stop'increase njmax!!!'
        iqj(nj)=icp1                             !leading jet parton
        ncc(2,1)=nj
        eqj(1,nj)=.5d0*(wps1+wmq1)
        eqj(2,nj)=.5d0*(wps1-wmq1)
        eqj(3,nj)=ptpf(1,1)
        eqj(4,nj)=ptpf(1,2)
       else
        if(qgran(b10).lt..3333d0)then
         iqc(1)=3
         icst=-4
         am2=am(3)
        else
         iqc(1)=int(1.5d0+qgran(b10))
         icst=-iqc(1)
         am2=am(1)
        endif
        if(qgran(b10).lt..5d0)then
         iqc(1)=-iqc(1)
         icst=-icst
         icd=icp1
         wpd=wps1
         icp1=icp2
         icp2=icd
         wps1=wps2
         wps2=wpd
         do i=1,2
          ptd=ptpf(1,i)
          ptpf(1,i)=ptpf(2,i)
          ptpf(2,i)=ptd
         enddo
        endif
        
        nj=nj+1
        if(nj.gt.njmax)stop'increase njmax!!!'
        iqj(nj)=icp2                             !leading jet parton
        ncc(1,1)=nj                              !color conn. with leading jet
        ncc(2,1)=0
        wmq=(ptpf(2,1)**2+ptpf(2,2)**2)/wps2
        wmi=wmi-wmq
        eqj(1,nj)=.5d0*(wps2+wmq)
        eqj(2,nj)=.5d0*(wps2-wmq)
        eqj(3,nj)=ptpf(2,1)
        eqj(4,nj)=ptpf(2,2)
        icsp=icp1
        if(iabs(icp1).le.2)then
         icq=1
        elseif(iabs(icp1).eq.4)then
         icq=3
        else
         icq=2
        endif
        am1=am(icq)
        wmq=(ptpf(1,1)**2+ptpf(1,2)**2)/wps1
        wmi=wmi-wmq
        ep1(1)=.5d0*(wps1+wmq)
        ep1(2)=.5d0*(wps1-wmq)
        ep1(3)=ptpf(1,1)
        ep1(4)=ptpf(1,2)

        if(idp.eq.7.or.idp.eq.8)then
         bett=beth(1)
         bbbt=bbbpom
         xph=wph/wp1
        else
         bett=beth(icz)
         bbbt=bbbi(icdp,icz)
         xph=wph/wp0
        endif
3       zg=1.d0-qgran(b10)*(1.d0-xph)            !gluon splitting into qq~
        if(qgran(b10).gt.zg**dels*((1.d0-xph/zg)/(1.d0-xph))**bett
     *  *(.5d0+xph/zg/2.d0)**bbbt*(zg**2+(1.d0-zg)**2))goto 3
        wmq=zg/(1.d0-zg)/wph
        wmi=wmi-wmq
        if(wmi.lt.0.d0.or.wmi*wpi.le.s2min)goto 1
        ep2(1)=.5d0*wmq
        ep2(2)=-ep2(1)
        ep2(3)=0.d0
        ep2(4)=0.d0
        do i=1,4
         ept(i)=ep1(i)+ep2(i)
         ebal(i)=ebal(i)-ept(i)
        enddo
        if(qgnrm(ept).lt.(am1+am2)**2)goto 1
        call qggene(ep1,ep2,icsp,icst,0)         !string fragmentation
       endif
      endif

      if(idp.eq.3.or.idp.eq.4.or.idp.eq.8)then   !q_v targ. involved
       ic=(jc+1)/2  
       call qgvdef(ic,ic1,ic2,izt)               !leading state flavor
       iqc(2)=ic1                                !lower leg parton
       nj=nj+1
       if(nj.gt.njmax)stop'increase njmax!!!'
       iqj(nj)=ic2                               !leading jet parton
       ncc(1,2)=nj                               !color conn. with leading jet
       ncc(2,2)=0
       wpq=(ptts(2,1)**2+ptts(2,2)**2)/wm1
       wpi=wpi-wpq
       eqj(1,nj)=.5d0*(wpq+wm1)
       eqj(2,nj)=.5d0*(wpq-wm1)
       eqj(3,nj)=ptts(2,1)
       eqj(4,nj)=ptts(2,2)
       if(wpi.lt.0.d0.or.wpi*wmi.le.s2min)goto 10

      else                                       !g or q_sea as lower leg parton
       do i=1,2
       do j=1,2
        pttf(j,i)=ptts(j,i)
       enddo
       enddo
       if(idp.eq.1.or.idp.eq.2.or.idp.eq.7)then
        if(iqvt.ne.0)then
         if(lqb.ne.1.or.qgran(b10).gt.wpiex(1,2)+wpiex(2,2))then
          call qgvlc(wm1,wms1,wms2,izt,ict1,ict2)
         else
          if(qgran(b10).lt.wpiex(1,2)/(wpiex(1,2)+wpiex(2,2)))then     
           jex=1
           call qgpiex(xnt,ptt,am2ht,izt,ictg,ict1,ict2,2)
          else
           jex=2
           call qgrex(xnt,ptt,am2ht,izt,ictg,ict1,ict2,2,jex)
          endif
          xm=qgslc(arr(1),-arr(1))
          call qgcs(ctg,stg)
          pt2tg=(ptts(1,1)+ptts(2,1)+ptt*ctg)**2
     *    +(ptts(1,2)+ptts(2,2)+ptt*stg)**2
          wtgm=xnt*wm1
          wtgp=(am2ht+pt2tg)/wtgm
          wpi=wpi-wtgp
          if(wpi.lt.0.d0.or.wpi*wmi.le.s2min)goto 1
          ep2(1)=.5d0*(wtgp+wtgm)
          ep2(2)=.5d0*(wtgp-wtgm)
          ep2(3)=ptts(1,1)+ptts(2,1)+ptt*ctg
          ep2(4)=ptts(1,2)+ptts(2,2)+ptt*stg
          do i=1,4
           ebal(i)=ebal(i)-ep2(i)
          enddo
          call qgreg(ep2,ictg)
          wms1=xm*(wm1-wtgm)
          wms2=wm1-wtgm-wms1
          pttf(1,1)=-ptt*ctg*xm
          pttf(1,2)=-ptt*stg*xm
          pttf(2,1)=-ptt*ctg*(1.d0-xm)
          pttf(2,2)=-ptt*stg*(1.d0-xm)
         endif
        else
         ict1=int(1.5d0+qgran(b10))
         ict2=-ict1
         wms1=wm1
         wms2=wm2
        endif
       elseif(idp.eq.5.or.idp.eq.6)then          !g or q_sea from 3P-vertex
        ict1=int(1.5d0+qgran(b10))
        ict2=-ict1
        wms1=wm1*qgslc(alpq,-alpq)
        wms2=wm1-wms1
       endif

       if(jc.le.2)then                          !g as lower leg parton
        iqc(2)=0                                !lower leg parton
        nj=nj+1
        if(nj.gt.njmax)stop'increase njmax!!!'
        iqj(nj)=ict2                            !leading jet parton
        ncc(1,2)=nj                             !color connect. with leading jet
        wpq2=(pttf(2,1)**2+pttf(2,2)**2)/wms2
        wpq1=(pttf(1,1)**2+pttf(1,2)**2)/wms1
         if(wpi-wpq1-wpq2.le.0.d0.or.wmi*(wpi-wpq1-wpq2).le.s2min)then
         if(wms1.gt.wms2)then
          pttf(1,1)=pttf(1,1)+pttf(2,1)
          pttf(1,2)=pttf(1,2)+pttf(2,2)
          pttf(2,1)=0.d0
          pttf(2,2)=0.d0
          wpq1=(pttf(1,1)**2+pttf(1,2)**2)/wms1
          wpq2=0.d0
         else
          pttf(2,1)=pttf(2,1)+pttf(1,1)
          pttf(2,2)=pttf(2,2)+pttf(1,2)
          pttf(1,1)=0.d0
          pttf(1,2)=0.d0
          wpq2=(pttf(2,1)**2+pttf(2,2)**2)/wms2
          wpq1=0.d0
         endif
        endif
        wpi=wpi-wpq1-wpq2
        if(wpi.lt.0.d0.or.wpi*wmi.le.s2min)goto 1
        eqj(1,nj)=.5d0*(wpq2+wms2)
        eqj(2,nj)=.5d0*(wpq2-wms2)
        eqj(3,nj)=pttf(2,1)
        eqj(4,nj)=pttf(2,2)
        nj=nj+1
        if(nj.gt.njmax)stop'increase njmax!!!'
        iqj(nj)=ict1                            !leading jet parton
        ncc(2,2)=nj
        eqj(1,nj)=.5d0*(wpq1+wms1)
        eqj(2,nj)=.5d0*(wpq1-wms1)
        eqj(3,nj)=pttf(1,1)
        eqj(4,nj)=pttf(1,2)
       else
        if(qgran(b10).lt..3333d0)then
         iqc(2)=3
         icsp=-4
         am2=am(3)
        else
         iqc(2)=int(1.5d0+qgran(b10))
         icsp=-iqc(2)
         am2=am(1)
        endif
        if(qgran(b10).lt..5d0)then
         iqc(2)=-iqc(2)
         icsp=-icsp
         icd=ict1
         wmd=wms1
         ict1=ict2
         ict2=icd
         wms1=wms2
         wms2=wmd
         do i=1,2
          ptd=pttf(1,i)
          pttf(1,i)=pttf(2,i)
          pttf(2,i)=ptd
         enddo
        endif
        
        nj=nj+1
        if(nj.gt.njmax)stop'increase njmax!!!'
        iqj(nj)=ict2                            !leading jet parton
        ncc(1,2)=nj                             !color connect. with leading jet
        ncc(2,2)=0
        wpq=(pttf(2,1)**2+pttf(2,2)**2)/wms2
        wpi=wpi-wpq
        eqj(1,nj)=.5d0*(wpq+wms2)
        eqj(2,nj)=.5d0*(wpq-wms2)
        eqj(3,nj)=pttf(2,1)
        eqj(4,nj)=pttf(2,2)
        icst=ict1
        if(iabs(ict1).le.2)then
         icq=1
        elseif(iabs(ict1).eq.4)then
         icq=3
        else
         icq=2
        endif
        am1=am(icq)
        wpq=(pttf(1,1)**2+pttf(1,2)**2)/wms1
        wpi=wpi-wpq
        ep1(1)=.5d0*(wpq+wms1)
        ep1(2)=.5d0*(wpq-wms1)
        ep1(3)=pttf(1,1)
        ep1(4)=pttf(1,2)

        if(idp.eq.5.or.idp.eq.6)then
         bett=beth(1)
         bbbt=bbbpom
         xmh=wmh/wm1
        else
         bett=beth(2)
         bbbt=bbbi(icdt,2)
         xmh=wmh/wm0
        endif
4       zg=1.d0-qgran(b10)*(1.d0-xmh)           !gluon splitting into qq~
        if(qgran(b10).gt.zg**dels*((1.d0-xmh/zg)/(1.d0-xmh))**bett
     *  *(.5d0+xmh/zg/2.d0)**bbbt*(zg**2+(1.d0-zg)**2))goto 4
        wpq=zg/(1.d0-zg)/wmh
        wpi=wpi-wpq
        if(wpi.lt.0.d0.or.wmi*wpi.le.s2min)goto 1
        ep2(1)=.5d0*wpq
        ep2(2)=ep2(1)
        ep2(3)=0.d0
        ep2(4)=0.d0
        do i=1,4
         ept(i)=ep1(i)+ep2(i)
         ebal(i)=ebal(i)-ept(i)
        enddo
        if(qgnrm(ept).lt.(am1+am2)**2)goto 1
        call qggene(ep2,ep1,icsp,icst,0)        !string fragmentation
       endif
      endif
      
      if(idp.ne.4.and.iqc(1).ne.0.and.iqc(2).ne.0)then
       sjq=qgjit(qt0,qt0,wph*wmh,2,2)
       sjqq=qgjit(qt0,qt0,wph*wmh,3,2)
       sjqa=qgjit(qt0,qt0,wph*wmh,4,2)
       if(iqc(1).eq.iqc(2))then
        gbf=sjqq/(sjq/1.5d0+sjqq/6.d0+sjqa/6.d0)
       elseif(iqc(1)+iqc(2).eq.0)then
        gbf=sjqa/(sjq/1.5d0+sjqq/6.d0+sjqa/6.d0)
       else
        gbf=sjq/(sjq/1.5d0+sjqq/6.d0+sjqa/6.d0)
       endif
       gbf=gbf/1.7d0
       if(qgran(b10).gt.gbf)goto 1
      endif

      if(wpi.lt.wph.or.wmi.lt.wmh)then
       ic1=min(1,iabs(iqc(1)))+1
       ic2=min(1,iabs(iqc(2)))+1
       if(ic1+ic2.eq.4.and.iqc(1).eq.iqc(2))then
        gb=qgjit(qt0,qt0,wpi*wmi,3,2)/qgjit(qt0,qt0,wph*wmh,3,2)
       elseif(ic1+ic2.eq.4.and.iqc(1)+iqc(2).eq.0)then
        gb=qgjit(qt0,qt0,wpi*wmi,4,2)/qgjit(qt0,qt0,wph*wmh,4,2)
       else
        gb=qgjit(qt0,qt0,wpi*wmi,ic1,ic2)
     *  /qgjit(qt0,qt0,wph*wmh,ic1,ic2)
       endif
       if(qgran(b10).gt.gb)goto 1
      endif

c-------------------------------------------------
c hard parton cascade
      ept(1)=.5d0*(wpi+wmi)                     !ladder 4-momentum
      ept(2)=.5d0*(wpi-wmi)
      ept(3)=0.d0
      ept(4)=0.d0
      qmin(1)=qt0                               !q^2 cutoff for the upper leg
      qmin(2)=qt0                               !q^2 cutoff for the lower leg
      qminn=max(qmin(1),qmin(2))                !overall q^2 cutoff
      si=qgnrm(ept)
      jini=1
      jj=int(1.5d0+qgran(b10))      !1st upper (jj=1) or lower (jj=2) leg parton

5     if(debug.ge.3)write (moniou,203)si,iqc,ept
      pt2=ept(3)**2+ept(4)**2
      pt=dsqrt(pt2)
      ww=si+pt2
      iqp(1)=min(1,iabs(iqc(1)))+1
      iqp(2)=min(1,iabs(iqc(2)))+1
      wp(1)=ept(1)+ept(2)                       !LC+ for the ladder
      wp(2)=ept(1)-ept(2)                       !LC- for the ladder
      s2min=4.d0*fqscal*qminn    !minimal energy squared for 2-parton production
      if(jini.eq.1)then                         !general ladder
       if(iqp(1)+iqp(2).eq.4.and.iqc(1).eq.iqc(2))then
        sj=qgjit(qmin(jj),qmin(3-jj),si,3,2)    !total ladder
        sj1=qgjit1(qmin(3-jj),qmin(jj),si,3,2)  !one-way ordered
        sjb=qgbit(qmin(1),qmin(2),si,3,2)       !born contribution
       elseif(iqp(1)+iqp(2).eq.4.and.iqc(1)+iqc(2).eq.0)then
        sj=qgjit(qmin(jj),qmin(3-jj),si,4,2)    !total ladder
        sj1=qgjit1(qmin(3-jj),qmin(jj),si,4,2)  !one-way ordered
        sjb=qgbit(qmin(1),qmin(2),si,4,2)       !born contribution
       else
        sj=qgjit(qmin(jj),qmin(3-jj),si,iqp(jj),iqp(3-jj))   !total ladder
        sj1=qgjit1(qmin(3-jj),qmin(jj),si,iqp(3-jj),iqp(jj)) !one-way ordered
        sjb=qgbit(qmin(1),qmin(2),si,iqp(1),iqp(2))          !born contribution
       endif
       aks=qgran(b10)
       if(aks.lt.sjb/sj)then
        goto 7                                  !born process sampled
       elseif(aks.lt.sj1/sj)then                !change to ordered ladder
        jj=3-jj
        sj=sj1
        jini=0
       endif
      else                                      !one-way ordered ladder
       if(iqp(1)+iqp(2).eq.4.and.iqc(1).eq.iqc(2))then
        sj=qgjit1(qmin(jj),qmin(3-jj),si,3,2)   !one-way ordered
        sjb=qgbit(qmin(1),qmin(2),si,3,2)       !born contribution
       elseif(iqp(1)+iqp(2).eq.4.and.iqc(1)+iqc(2).eq.0)then
        sj=qgjit1(qmin(jj),qmin(3-jj),si,4,2)   !one-way ordered
        sjb=qgbit(qmin(1),qmin(2),si,4,2)       !born contribution
       else
        sj=qgjit1(qmin(jj),qmin(3-jj),si,iqp(jj),iqp(3-jj)) !one-way ordered
        sjb=qgbit(qmin(1),qmin(2),si,iqp(1),iqp(2))         !born contribution
       endif
       if(qgran(b10).lt.sjb/sj)goto 7           !born process sampled
      endif
      wwmin=(s2min+qmin(jj)+pt2-2.d0*pt*dsqrt(qmin(jj)*epsxmn))
     */(1.d0-epsxmn)             !minimal energy squared for 3-parton production
      if(debug.ge.3)write (moniou,204)s2min,wwmin,sj,sjb
      if(ww.lt.1.1d0*wwmin)goto 7               !energy too low -> born process

      xxx=pt*dsqrt(qmin(jj))/ww
      xmin=(s2min+qmin(jj)+pt2)/ww
      xmin=xmin-2.d0*xxx*(xxx+dsqrt(xxx**2+1.d0-xmin))
      xmax=1.d0-epsxmn
      if(debug.ge.3)write (moniou,205)xmin,xmax

      qqmax=(pt*dsqrt(epsxmn)+dsqrt(max(0.d0,pt2*epsxmn
     *+(1.d0+4.d0*fqscal)*(xmax*ww-pt2))))/(1.d0+4.d0*fqscal)
      qqmin=qmin(jj)        !minimal parton virtuality in the current rung
      if(debug.ge.3)write (moniou,206)qqmin,qqmax

      qm0=qqmin
      xm0=xmax
      s2max=xm0*ww

      if(jini.eq.1)then
       sj0=qgjit(qm0,qmin(3-jj),s2max,1,iqp(3-jj))*qgfap(xm0,iqp(jj),1)
       if(iqp(3-jj).eq.2.and.iqp(jj).eq.1)then
        sj0=sj0+qgfap(xm0,iqp(jj),2)*(qgjit(qm0,qmin(3-jj),s2max,2,2)
     *  /1.5d0+qgjit(qm0,qmin(3-jj),s2max,3,2)/6.d0
     *  +qgjit(qm0,qmin(3-jj),s2max,4,2)/6.d0)
       elseif(iqp(3-jj).eq.2.and.iqc(1).eq.iqc(2))then
        sj0=sj0+qgfap(xm0,iqp(jj),2)*qgjit(qm0,qmin(3-jj),s2max,3,2)
       elseif(iqp(3-jj).eq.2.and.iqc(1)+iqc(2).eq.0)then
        sj0=sj0+qgfap(xm0,iqp(jj),2)*qgjit(qm0,qmin(3-jj),s2max,4,2)
       else
        sj0=sj0+qgfap(xm0,iqp(jj),2)
     *  *qgjit(qm0,qmin(3-jj),s2max,2,iqp(3-jj))
       endif
      else
       sj0=qgjit1(qm0,qmin(3-jj),s2max,1,iqp(3-jj))
     * *qgfap(xm0,iqp(jj),1)
       if(iqp(3-jj).eq.2.and.iqp(jj).eq.1)then
        sj0=sj0+qgfap(xm0,iqp(jj),2)*(qgjit1(qm0,qmin(3-jj),s2max,2,2)
     *  /1.5d0+qgjit1(qm0,qmin(3-jj),s2max,3,2)/6.d0
     *  +qgjit1(qm0,qmin(3-jj),s2max,4,2)/6.d0)
       elseif(iqp(3-jj).eq.2.and.iqc(1).eq.iqc(2))then
        sj0=sj0+qgfap(xm0,iqp(jj),2)*qgjit1(qm0,qmin(3-jj),s2max,3,2)
       elseif(iqp(3-jj).eq.2.and.iqc(1)+iqc(2).eq.0)then
        sj0=sj0+qgfap(xm0,iqp(jj),2)*qgjit1(qm0,qmin(3-jj),s2max,4,2)
       else
        sj0=sj0+qgfap(xm0,iqp(jj),2)
     *  *qgjit1(qm0,qmin(3-jj),s2max,2,iqp(3-jj))
       endif
      endif
      gb0=sj0*qm0*qgalf(qm0/alm)*qgsudx(qm0,iqp(jj)) *4.5d0  !normal. of accept.
      if(xm0.le..5d0)then
       gb0=gb0*(2.d0*xm0)**(1.d0-delh)
      else
       gb0=gb0*(1.d0-xm0)*2.d0
      endif
      if(debug.ge.3)write (moniou,208)xm0,xmin,xmax,gb0

      xmin2=max(.5d0,xmin)
      xmin1=xmin**delh
      xmax1=min(xmax,.5d0)**delh
      if(xmin.ge..5d0)then                             !choose proposal function
       djl=1.d0
      elseif(xmax.lt..5d0)then
       djl=0.d0
      else
       djl=1.d0/(1.d0+((2.d0*xmin)**delh-1.d0)/delh
     * /dlog(2.d0*(1.d0-xmax)))
      endif

c-------------------------------------------------
c propose x, q^2
6     if(qgran(b10).gt.djl)then
       x=(xmin1+qgran(b10)*(xmax1-xmin1))**(1.d0/delh) !parton LC share
      else
       x=1.d0-(1.d0-xmin2)*((1.d0-xmax)/(1.d0-xmin2))**qgran(b10)
      endif
      qq=qqmin/(1.d0+qgran(b10)*(qqmin/qqmax-1.d0))    !parton virtuality
      qt2=qq*(1.d0-x)                                  !parton p_t^2
      if(debug.ge.4)write (moniou,209)qq,qqmin,qqmax,x,qt2

      qmin2=max(qq,qminn)
      qt=dsqrt(qt2)
      call qgcs(c,s)
      ep3(3)=qt*c                               !final parton p_x, p_y
      ep3(4)=qt*s
      pt2new=(ept(3)-ep3(3))**2+(ept(4)-ep3(4))**2!p_t^2 for the remained ladder
      s2min2=max(s2min,4.d0*fqscal*qmin2)  !new ladder kinematic limit
      s2=x*ww-qt2*x/(1.d0-x)-pt2new        !mass squared for the remained ladder
      if(s2.lt.s2min2)goto 6           !ladder mass below threshold -> rejection

      if(jini.eq.1)then                    !weights for g- and q-legs
       sj1=qgjit(qq,qmin(3-jj),s2,1,iqp(3-jj))*qgfap(x,iqp(jj),1)
       sq=qgjit(qq,qmin(3-jj),s2,2,iqp(3-jj))
       sqq=0.d0
       sqa=0.d0
       if(iqp(3-jj).eq.2.and.iqp(jj).eq.1)then
        sqq=qgjit(qq,qmin(3-jj),s2,3,2)
        sqa=qgjit(qq,qmin(3-jj),s2,4,2)
        sj2=qgfap(x,iqp(jj),2)*(sq/1.5d0+sqq/6.d0+sqa/6.d0)
       elseif(iqp(3-jj).eq.2.and.iqc(1).eq.iqc(2))then
        sqq=qgjit(qq,qmin(3-jj),s2,3,2)
        sj2=qgfap(x,iqp(jj),2)*sqq
       elseif(iqp(3-jj).eq.2.and.iqc(1)+iqc(2).eq.0)then
        sqa=qgjit(qq,qmin(3-jj),s2,4,2)
        sj2=qgfap(x,iqp(jj),2)*sqa
       else
        sj2=qgfap(x,iqp(jj),2)*sq
       endif
      else
       sj1=qgjit1(qq,qmin(3-jj),s2,1,iqp(3-jj))*qgfap(x,iqp(jj),1)
       sq=qgjit1(qq,qmin(3-jj),s2,2,iqp(3-jj))
       sqq=0.d0
       sqa=0.d0
       if(iqp(3-jj).eq.2.and.iqp(jj).eq.1)then
        sqq=qgjit1(qq,qmin(3-jj),s2,3,2)
        sqa=qgjit1(qq,qmin(3-jj),s2,4,2)
        sj2=qgfap(x,iqp(jj),2)*(sq/1.5d0+sqq/6.d0+sqa/6.d0)
       elseif(iqp(3-jj).eq.2.and.iqc(1).eq.iqc(2))then
        sqq=qgjit1(qq,qmin(3-jj),s2,3,2)
        sj2=qgfap(x,iqp(jj),2)*sqq
       elseif(iqp(3-jj).eq.2.and.iqc(1)+iqc(2).eq.0)then
        sqa=qgjit1(qq,qmin(3-jj),s2,4,2)
        sj2=qgfap(x,iqp(jj),2)*sqa
       else
        sj2=qgfap(x,iqp(jj),2)*sq
       endif
      endif
      gb7=(sj1+sj2)*qgalf(qq/alm)*qq*qgsudx(qq,iqp(jj))/gb0  /3.1d0
                               !acceptance probability for x and q**2 simulation
      if(x.le..5d0)then
       gb7=gb7*(2.d0*x)**(1.d0-delh)
      else
       gb7=gb7*(1.d0-x)*2.d0
      endif
      if(debug.ge.4)write (moniou,210)gb7,s2,sj1,sj2,jj,jini
      if(qgran(b10).gt.gb7)goto 6               !rejection

c-------------------------------------------------
c define color flow for the emitted jet; perform final state emission
      nqc(2)=0
      if(qgran(b10).lt.sj1/(sj1+sj2))then       !new gluon-leg ladder
       if(iqc(jj).eq.0)then                     !g -> gg
        jt=1
        jq=int(1.5d0+qgran(b10))
        nqc(1)=ncc(jq,jj)                       !color connection for the jet
       else                                     !q -> qg
        jt=2
        if(iqc(jj).gt.0)then                    !orientation of color flow
         jq=1
        else
         jq=2
        endif
        nqc(1)=0
        ncc(jq,jj)=ncc(1,jj)                    !color connection for the jet
       endif
       iq1=iqc(jj)                              !jet flavor (type)
       iqc(jj)=0                                !new ladder leg flavor (type)

      else                                      !new quark-leg ladder
       if(iqc(jj).ne.0)then                     !q -> gq
        iq1=0
        jt=3
        if(iqc(jj).gt.0)then                    !orientation of color flow
         jq=1
        else
         jq=2
        endif
        nqc(1)=ncc(1,jj)                        !color connection for the jet
       else                                     !g -> qq~
        jq=int(1.5d0+qgran(b10))                !orientation of color flow
66      iq1=int(3.d0*qgran(b10)+1.d0)*(3-2*jq)  !jet flavor (type)
        if(iqp(3-jj).ne.1)then
         if(jini.eq.1)then
          sjq=qgjit(qq,qmin(3-jj),s2,2,2)
          sjqq=qgjit(qq,qmin(3-jj),s2,3,2)
          sjqa=qgjit(qq,qmin(3-jj),s2,4,2)
         else
          sjq=qgjit1(qq,qmin(3-jj),s2,2,2)
          sjqq=qgjit1(qq,qmin(3-jj),s2,3,2)
          sjqa=qgjit1(qq,qmin(3-jj),s2,4,2)
         endif
         if(iq1.eq.iqc(2))then
          gbf=sjqa/(sjq/1.5d0+sjqq/6.d0+sjqa/6.d0)
         elseif(iq1+iqc(2).eq.0)then
          gbf=sjqq/(sjq/1.5d0+sjqq/6.d0+sjqa/6.d0)
         else
          gbf=sjq/(sjq/1.5d0+sjqq/6.d0+sjqa/6.d0)
         endif
         gbf=gbf/1.7d0
         if(qgran(b10).gt.gbf)goto 66
        endif
        iqc(jj)=-iq1                            !new ladder leg flavor (type)
        jt=4
        nqc(1)=ncc(jq,jj)                       !color connections for the jet
        ncc(1,jj)=ncc(3-jq,jj)
       endif
      endif
      if(debug.ge.3)write (moniou,211)jt

      call qgcjet(qt2,iq1,qv1,zv1,qm1,iqv1,ldau1,lpar1,jq) !final state emission
      si=x*ww-(qt2+qm1(1,1))*x/(1.d0-x)-pt2new  !mass squared for the new ladder
      if(si.gt.s2min2)then
       iq=min(1,iabs(iqc(jj)))+1
       if(jini.eq.1)then
        if(iq+iqp(3-jj).eq.4.and.iqc(1).eq.iqc(2))then
         gb=qgjit(qq,qmin(3-jj),si,3,2)/qgjit(qq,qmin(3-jj),s2,3,2)
        elseif(iq+iqp(3-jj).eq.4.and.iqc(1)+iqc(2).eq.0)then
         gb=qgjit(qq,qmin(3-jj),si,4,2)/qgjit(qq,qmin(3-jj),s2,4,2)
        else
         gb=qgjit(qq,qmin(3-jj),si,iq,iqp(3-jj))
     *   /qgjit(qq,qmin(3-jj),s2,iq,iqp(3-jj))
        endif
       else
        if(iq+iqp(3-jj).eq.4.and.iqc(1).eq.iqc(2))then
         gb=qgjit1(qq,qmin(3-jj),si,3,2)/qgjit1(qq,qmin(3-jj),s2,3,2)
        elseif(iq+iqp(3-jj).eq.4.and.iqc(1)+iqc(2).eq.0)then
         gb=qgjit1(qq,qmin(3-jj),si,4,2)/qgjit1(qq,qmin(3-jj),s2,4,2)
        else
         gb=qgjit1(qq,qmin(3-jj),si,iq,iqp(3-jj))
     *   /qgjit1(qq,qmin(3-jj),s2,iq,iqp(3-jj))
        endif
       endif
       if(qgran(b10).gt.gb)goto 1               !weight for jet mass correction
      else                                      !below threshold -> rejection
       goto 1
      endif

      wp3=wp(jj)*(1.d0-x)
      wm3=(qt2+qm1(1,1))/wp3
      ep3(1)=.5d0*(wp3+wm3)                     !jet 4-momentum
      ep3(2)=.5d0*(wp3-wm3)*(3-2*jj)
      call qgrec(ep3,nqc,qv1,zv1,qm1,iqv1,ldau1,lpar1,jq)
                              !reconstruction of 4-momenta for all final partons
c-------------------------------------------------
c define color connections for the new ladder
      if(jt.eq.1)then
       if(ncc(1,jj).eq.0.and.ncc(2,jj).eq.0)ncc(3-jq,jj)=nqc(1)
       ncc(jq,jj)=nqc(2)
      elseif(jt.eq.2)then
       ncc(3-jq,jj)=nqc(1)
      elseif(jt.eq.3)then
       ncc(1,jj)=nqc(2)
      elseif(jt.eq.4.and.ncc(1,jj).eq.0.and.ncc(2,jj).eq.0)then
       ncc(1,jj)=nqc(1)
      endif

      if(iabs(iq1).eq.3)then
       iqqq=8+iq1/3*4
      else
       iqqq=8+iq1
      endif
      if(debug.ge.3)write (moniou,212)tyq(iqqq),qt2,ep3
      do i=1,4
       ept(i)=ept(i)-ep3(i)                     !new ladder 4-momentum
      enddo
      qmin(jj)=qq                               !new virtuality cutoffs
      qminn=qmin2
      goto 5                                    !consider next parton emission

c------------------------------------------------
c born process - last parton pair production in the ladder
7     if(debug.ge.2)write (moniou,214)si,qminn,iqc
      tmin=qminn*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qminn*fqscal/si)))
      qtmin=tmin*(1.d0-tmin/si)
      fb=0.d0
      if(iqc(1).ne.0.or.iqc(2).ne.0)then
       do n=1,3
        fb=fb+qgfbor(si,tmin,iqc(1),iqc(2),n)
       enddo
       gb0=tmin**2*qgalf(qtmin/fqscal/alm)**2*fb*1.15d0 !norm. of acceptance
      else
       do n=1,3
        fb=fb+qgfbor(si,.5d0*si,iqc(1),iqc(2),n)
       enddo
       gb0=.25d0*si**2*qgalf(qtmin/fqscal/alm)**2*fb
      endif
      gb0=gb0*qgsudx(qtmin/fqscal,iqp(1))*qgsudx(qtmin/fqscal,iqp(2))
      if(debug.ge.3)write (moniou,215)gb0

8     q2=tmin/(1.d0-qgran(b10)*(1.d0-2.d0*tmin/si)) !proposed q^2
      z=q2/si                                       !parton LC momentum share
      qt2=q2*(1.d0-z)                               !parton p_t^2
      if(qgran(b10).lt..5d0)then
       jm=2
       tq=si-q2
      else
       jm=1
       tq=q2
      endif
      wtot=0.d0
      do n=1,3
       wsub(n)=qgfbor(si,tq,iqc(1),iqc(2),n)
       wtot=wtot+wsub(n)
      enddo
      gb=q2**2*qgalf(qt2/fqscal/alm)**2*wtot*qgsudx(qt2/fqscal,iqp(1))
     **qgsudx(qt2/fqscal,iqp(2))/gb0                !acceptance probabilty
      if(debug.ge.4)write (moniou,216)gb,q2,z,qt2
      if(qgran(b10).gt.gb)goto 8                    !rejection

      nsub=1                                        !subprocess
      do i=1,2
       iqb(i)=iqc(i)
      enddo
      aks=wtot*qgran(b10)
      if(aks.gt.wsub(1))then
       if(aks.gt.wsub(1)+wsub(2))then
        nsub=3
        do i=1,2
         iqb(i)=0
        enddo
       else
        nsub=2
        if(iqc(1).ne.0)then
         iqb(1)=int(3.d0*qgran(b10)+1.d0)*iabs(iqc(1))/iqc(1)
        else
         iqb(1)=int(3.d0*qgran(b10)+1.d0)*(1-2*int(.5d0+qgran(b10)))
        endif
        iqb(2)=-iqb(1)
       endif
      endif
      
c HT-corrections (rejection)
      if(idp.eq.2.or.idp.eq.4.or.idp.eq.5.or.idp.eq.6)then
       dhtp=1.d0
      else
       xg=qt0/scm*wm0/wp(2)
       if(ia(1).eq.1)then
        bbg=(xa(ip,1)+bcoll-xxp)**2+(xa(ip,2)-yyp)**2
        bbp=bbg
        pdfg=qgpdfbi(xg,bbg,0.d0,0.d0,icdp,icz,1,1)
       else
        pdfg=0.d0
        sumps=0.d0
        do ipp=1,ia(1)
         if(iconab(ipp,it).ne.0)then 
          bbg=(xa(ipp,1)+bcoll-xxp)**2+(xa(ipp,2)-yyp)**2
          if(ipp.eq.ip)bbp=bbg
          pdfg=pdfg+qgpdfbi(xg,bbg,1.d0-exp(-sumps),0.d0
     *    ,iddp(ipp),icz,1,1)
          if(xg*sgap.lt.1.d0)sumps=sumps+qgfani(1.d0/xg,bbg
     *    ,1.d0-exp(-sumps),0.d0,0.d0,iddp(ipp),icz,3)
         endif
        enddo
       endif
       facht=1.d0+pdfg/q2*qgalf(qt2/fqscal/alm)*2.d0*pi**3*htfac
     * *(3.d0*(1-min(1,iabs(iqb(1))))+4.d0/3.d0*min(1,iabs(iqb(1))))
       if(jc.eq.1.or.jc.eq.3)then
        iqpp=1
       else
        iqpp=2
       endif
       if(idp.eq.1.or.idp.eq.3)then
        xpp=wpi/wp0
        pdfp0=qgpdfbi(xpp,bbp,vvxps,0.d0,icdp,icz,iqpp,2)
        pdfp=qgpdfbi(xpp*facht,bbp,vvxps,0.d0,icdp,icz,iqpp,2)
       else
        xpp=wpi/wp1
        pdfp0=qgloopri(1.d0/xpp,bbi,vvxps,iqpp,2)
        pdfp=qgloopri(1.d0/xpp/facht,bbi,vvxps,iqpp,2)
       endif
       dhtp=pdfp/pdfp0/facht
      endif

      if(idp.eq.3.or.idp.eq.4.or.idp.eq.7.or.idp.eq.8)then
       dhtt=1.d0
      else
       xg=qt0/scm*wp0/wp(1)
       if(ia(2).eq.1)then
        bbg=(xb(it,1)-xxp)**2+(xb(it,2)-yyp)**2
        bbt=bbg
        pdfg=qgpdfbi(xg,bbg,0.d0,0.d0,icdt,2,1,1)
       else
        pdfg=0.d0
        sumts=0.d0
        do itt=1,ia(2)
         if(iconab(ip,itt).ne.0)then
          bbg=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
          if(itt.eq.it)bbt=bbg
          pdfg=pdfg+qgpdfbi(xg,bbg,1.d0-exp(-sumts),0.d0
     *    ,iddt(itt),2,1,1)
          if(xg*sgap.lt.1.d0)sumts=sumts+qgfani(1.d0/xg,bbg
     *    ,1.d0-exp(-sumts),0.d0,0.d0,iddt(itt),2,3)
         endif
        enddo
       endif
       facht=1.d0+pdfg/q2*qgalf(qt2/fqscal/alm)*2.d0*pi**3*htfac
     * *(3.d0*(1-min(1,iabs(iqb(2))))+4.d0/3.d0*min(1,iabs(iqb(2))))
       if(jc.le.2)then
        iqtt=1
       else
        iqtt=2
       endif
       if(idp.eq.1.or.idp.eq.2)then
        xmt=wmi/wm0
        pdft0=qgpdfbi(xmt,bbt,vvxts,0.d0,icdt,2,iqtt,2)
        pdft=qgpdfbi(xmt*facht,bbt,vvxts,0.d0,icdt,2,iqtt,2)
       else
        xmt=wmi/wm1
        pdft0=qgloopri(1.d0/xmt,bbi,vvxts,iqtt,2)
        pdft=qgloopri(1.d0/xmt/facht,bbi,vvxts,iqtt,2)
       endif
       dhtt=pdft/pdft0/facht
      endif
      gbht=dhtp*dhtt
      if(qgran(b10).gt.gbht)goto 1

c-------------------------------------------------
c define color connections for the 1st emitted jet
      nqc(2)=0
      if(nsub.eq.1)then
       if(iqc(1).eq.0.and.iqc(2).eq.0)then          !gg-process
        jq=int(1.5d0+qgran(b10))                    !orientation of color flow
        nqc(1)=ncc(jq,jm)
        if(qgran(b10).lt..5d0)then
         jt=1                                       !gg -> gg
         njc1=ncc(3-jq,jm)                        !color connections for 1st jet
         njc2=ncc(jq,3-jm)
         if(iqj(njc1).ne.0)then
          ncj(1,njc1)=njc2
         else
          ncj(jq,njc1)=njc2
         endif
         if(iqj(njc2).ne.0)then
          ncj(1,njc2)=njc1
         else
          ncj(3-jq,njc2)=njc1
         endif
        else                                !gg -> gg (inverse color connection)
         jt=2
         nqc(2)=ncc(3-jq,3-jm)
        endif

       elseif(iqc(1)*iqc(2).eq.0)then               !qg -> qg
        if(iqc(1)+iqc(2).gt.0)then                  !orientation of color flow
         jq=1
        else
         jq=2
        endif
        if(qgran(b10).lt..5d0)then
         if(iqc(jm).eq.0)then
          jt=3
          nqc(1)=ncc(jq,jm)
          njc1=ncc(3-jq,jm)
          njc2=ncc(1,3-jm)
          if(iqj(njc1).ne.0)then
           ncj(1,njc1)=njc2
          else
           ncj(jq,njc1)=njc2
          endif
          if(iqj(njc2).ne.0)then
           ncj(1,njc2)=njc1
          else
           ncj(3-jq,njc2)=njc1
          endif
         else
          jt=4
          nqc(1)=0
          njc1=ncc(1,jm)
          njc2=ncc(3-jq,3-jm)
          if(njc2.ne.0)then
           if(iqj(njc1).ne.0)then
            ncj(1,njc1)=njc2
           else
            ncj(3-jq,njc1)=njc2
           endif
           if(iqj(njc2).ne.0)then
            ncj(1,njc2)=njc1
           else
            ncj(jq,njc2)=njc1
           endif
          endif
         endif
        else
         if(iqc(jm).eq.0)then
          jt=5
          nqc(2)=ncc(3-jq,jm)
          nqc(1)=ncc(1,3-jm)
         else
          jt=6
          nqc(1)=ncc(jq,3-jm)
         endif
        endif

       elseif(iqc(1)*iqc(2).gt.0)then               !qq (q~q~) -> qq (q~q~)
        jt=7
        if(iqc(1).gt.0)then
         jq=1
        else
         jq=2
        endif
        nqc(1)=ncc(1,3-jm)
       elseif(iqc(1)*iqc(2).lt.0)then               !qq~ -> qq~
        jt=8
        if(iqc(jm).gt.0)then
         jq=1
        else
         jq=2
        endif
        nqc(1)=0
        njc1=ncc(1,jm)
        njc2=ncc(1,3-jm)
        if(iqj(njc1).ne.0)then
         ncj(1,njc1)=njc2
        else
         ncj(3-jq,njc1)=njc2
        endif
        if(iqj(njc2).ne.0)then
         ncj(1,njc2)=njc1
        else
         ncj(jq,njc2)=njc1
        endif
       endif

      elseif(nsub.eq.2.and.iqc(1)+iqc(2).eq.0
     *.and.iqc(1)*iqc(2).ne.0)then                  !qq~ -> q'q'~
       jt=9
       if(iqc(jm).gt.0)then
        jq=1
       else
        jq=2
       endif
       nqc(1)=ncc(1,jm)
      elseif(nsub.eq.2.and.iqc(1).eq.0.and.iqc(2).eq.0)then !gg -> q'q'~
       jt=10
       if(iqb(jm).gt.0)then
        jq=1
       else
        jq=2
       endif
       nqc(1)=ncc(jq,jm)
       njc1=ncc(3-jq,jm)
       njc2=ncc(jq,3-jm)
       if(iqj(njc1).ne.0)then
        ncj(1,njc1)=njc2
       else
        ncj(jq,njc1)=njc2
       endif
       if(iqj(njc2).ne.0)then
        ncj(1,njc2)=njc1
       else
        ncj(3-jq,njc2)=njc1
       endif
      elseif(nsub.eq.3.and.iqc(1)+iqc(2).eq.0
     *.and.iqc(1)*iqc(2).ne.0)then                 !qq~ -> gg
       jt=11
       if(iqc(jm).gt.0)then
        jq=1
       else
        jq=2
       endif
       nqc(1)=ncc(1,jm)
      else
       stop'qghot:wrong subprocess'
      endif

      if(jt.lt.8.or.jt.eq.11)then
       jq2=jq
      else
       jq2=3-jq
      endif
      if(debug.ge.3)write (moniou,211)jt
      call qgcjet(qt2,iqb(jm),qv1,zv1,qm1,iqv1,ldau1,lpar1,jq) !final state.conf
      call qgcjet(qt2,iqb(3-jm),qv2,zv2,qm2,iqv2,ldau2,lpar2,jq2)
      amt1=qt2+qm1(1,1)
      amt2=qt2+qm2(1,1)
      if(dsqrt(si).gt.dsqrt(amt1)+dsqrt(amt2))then
       z=qgtwd(si,amt1,amt2)
      else
       if(debug.ge.4)write (moniou,217)dsqrt(si),dsqrt(amt1),dsqrt(amt2)
       goto 1                                      !below threshold -> rejection
      endif
      call qgdeft(si,ept,ey)
      wp3=z*dsqrt(si)
      wm3=(qt2+qm1(1,1))/wp3
      ep3(1)=.5d0*(wp3+wm3)                         !1st jet 4-momentum
      ep3(2)=.5d0*(wp3-wm3)
      qt=dsqrt(qt2)
      call qgcs(c,s)
      ep3(3)=qt*c
      ep3(4)=qt*s
      call qgtran(ep3,ey,1)
      call qgrec(ep3,nqc,qv1,zv1,qm1,iqv1,ldau1,lpar1,jq) !jet reconstruction
      if(iabs(iqb(jm)).eq.3)then
       iqqq=8+iqb(jm)/3*4
      else
       iqqq=8+iqb(jm)
      endif
      if(debug.ge.3)write (moniou,212)tyq(iqqq),qt2,ep3
      wp3=(1.d0-z)*dsqrt(si)
      wm3=(qt2+qm2(1,1))/wp3
      ep3(1)=.5d0*(wp3+wm3)                         !2nd jet 4-momentum
      ep3(2)=.5d0*(wp3-wm3)
      ep3(3)=-qt*c
      ep3(4)=-qt*s
      call qgtran(ep3,ey,1)

c-------------------------------------------------
c define color connections for the 2nd emitted jet
      if(jt.eq.1)then
       nqc(1)=nqc(2)
       nqc(2)=ncc(3-jq,3-jm)
      elseif(jt.eq.2)then
       nqc(2)=ncc(3-jq,jm)
       nqc(1)=ncc(jq,3-jm)
      elseif(jt.eq.3)then
       nqc(1)=nqc(2)
      elseif(jt.eq.4)then
       nqc(2)=nqc(1)
       nqc(1)=ncc(jq,3-jm)
      elseif(jt.eq.5)then
       nqc(1)=ncc(jq,jm)
      elseif(jt.eq.6)then
       nqc(2)=ncc(3-jq,3-jm)
       nqc(1)=ncc(1,jm)
      elseif(jt.eq.7)then
       nqc(1)=ncc(1,jm)
      elseif(jt.eq.9)then
       nqc(1)=ncc(1,3-jm)
      elseif(jt.eq.10)then
       nqc(1)=ncc(3-jq,3-jm)
      elseif(jt.eq.11)then
       nqc(1)=nqc(2)
       nqc(2)=ncc(1,3-jm)
      endif
      call qgrec(ep3,nqc,qv2,zv2,qm2,iqv2,ldau2,lpar2,jq2) !jet reconstruction
      if(iabs(iqb(3-jm)).eq.3)then
       iqqq=8+iqb(3-jm)/3*4
      else
       iqqq=8+iqb(3-jm)
      endif
      if(debug.ge.3)write (moniou,212)tyq(iqqq),qt2,ep3

      if(nj.eq.nj0)stop'qghot: nj=nj0?'
      do i=nj0+1,nj
       do j=1,4
        ebal(j)=ebal(j)-eqj(j,i)
       enddo
      enddo
      if(ebal(1).gt.1.d-5)write(*,*)'ebal-hot',ebal,nj0,nj,idp

      if(debug.ge.2)write (moniou,218)nj
      if(debug.ge.5)write (moniou,219)ebal
      if(debug.ge.1)write (moniou,220)
      return
10    iret=1      

201   format(2x,'qghot - semihard interaction:'/
     *4x,'type of the interaction - ',i2/
     *4x,'initial light cone momenta - ',2e10.3/
     *4x,'remnant types - ',2i3,2x,'diffr. eigenstates - ',2i2/
     *4x,'proj. valence flag - ',i2,2x,'targ. valence flag - ',i2/
     *4x,'initial number of final partons - ',i4)
202   format(2x,'qghot: mass squared for parton ladder - ',e10.3)
203   format(2x,'qghot: ',' mass squared for the laddder:',e10.3/
     *4x,'ladder end flavors:',2i3/4x,'ladder 4-momentum: ',4e10.3)
204   format(2x,'qghot: kinematic bounds s2min=',e10.3,
     *2x,'wwmin=',e10.3/4x,'jet cross section sj=',e10.3,
     *2x,'born cross section sjb=',e10.3)
205   format(2x,'qghot: xmin=',e10.3,2x,'xmax=',e10.3)
206   format(2x,'qghot: qqmin=',e10.3,2x,'qqmax=',e10.3)
208   format(2x,'qghot: xm0=',e10.3,2x,'xmin=',e10.3,2x,
     *'xmax=',e10.3,2x,'gb0=',e10.3)
209   format(2x,'qghot: qq=',e10.3,2x,'qqmin=',e10.3,2x,
     *'qqmax=',e10.3,2x,'x=',e10.3,2x,'qt2=',e10.3)
210   format(2x,'qghot: gb7=',e10.3,2x,'s2=',e10.3,2x,'sj1=',e10.3
     *,2x,'sj2=',e10.3,2x,'jj=',i2,2x,'jini=',i2)
211   format(2x,'qghot: colour connection jt=:',i1)
212   format(2x,'qghot: new jet flavor:',a2,
     *' pt squared for the jet:',e10.3/4x,'jet 4-momentum:',4e10.3)
214   format(2x,'qghot - highest virtuality subprocess in the ladder:'/
     *4x,'mass squared for the process:',e10.3/4x,'q^2-cutoff:',e10.3
     *,2x,'iqc=',2i3)
215   format(2x,'qghot - normalization of acceptance:',' gb0=',e10.3)
216   format(2x,'qghot - acceptance probabilty:'/
     *4x,'gb=',e10.3,2x,'q2=',e10.3,2x,'z=',e10.3,2x,'qt2=',e10.3)
217   format(2x,'qghot: ecm=',e10.3,2x,'mt1=',e10.3,2x,'mt2=',e10.3)
218   format(2x,'qghot: total number of jets - ',i4)
219   format(2x,'qghot: 4-momentum balans - ',4e10.3)
220   format(2x,'qghot - end')
      return
      end

c------------------------------------------------------------------------
      function npgen(vv,npmin,npmax)
c-----------------------------------------------------------------------
c npgen -  Poisson distribution
c vv    - average number
c npmin - minimal number
c npmax - maximal number
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      if(npmin.eq.0)then
       aks=qgran(b10)
       vvn=exp(-vv)
       do n=1,npmax
        aks=aks-vvn
        if(aks.lt.0.d0)goto 1
        vvn=vvn*vv/dble(n)
       enddo
      elseif(npmin.eq.1)then
       aks=qgran(b10)*(1.d0-exp(-vv))
       vvn=exp(-vv)
       do n=1,npmax
        vvn=vvn*vv/dble(n)
        aks=aks-vvn
        if(aks.lt.0.d0)goto 2
       enddo
      elseif(npmin.eq.2)then
       aks=qgran(b10)*(1.d0-exp(-vv)*(1.d0+vv))
       vvn=vv*exp(-vv)
       do n=2,npmax
        vvn=vvn*vv/dble(n)
        aks=aks-vvn
        if(aks.lt.0.d0)goto 2
       enddo
      else
       stop'npgen'
      endif
1     n=n-1
2     npgen=n
      return
      end

c=============================================================================
      double precision function qgbit(qi,qj,s,m,l)
c------------------------------------------------------------------------
c qgbit - born cross-section interpolation
c qi,qj - effective momentum cutoffs for the scattering,
c s - total c.m. energy squared for the scattering,
c m - parton type at current end of the ladder (1 - g, 2 - q)
c l - parton type at opposite end of the ladder (1 - g, 2 - q)
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wi(3),wk(3)
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr20/ spmax
      common /qgarr26/ factk,fqscal
      common /qgarr31/ csj(40,240)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)qi,qj,s,m,l
      qgbit=0.d0
      qq=max(qi,qj)
      s2min=qq*4.d0*fqscal
      if(s.le..99d0*s2min)then
       if(debug.ge.3)write (moniou,202)qgbit
       return
      endif

      tmin=qq*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qq*fqscal/s)))
      qli=dlog(qq)/dlog(spmax/4.d0/fqscal)*39.d0+1.d0
      sl=dlog(s/s2min)/dlog(spmax/s2min)*39.d0+1.d0
      i=min(38,int(qli))
      k=min(38,int(sl))
      mm=max(m,l)
      if(mm.le.2)then
       ml=40*(m-1)+80*(l-1)
      else
       ml=40*(mm+1)
      endif

      wk(2)=sl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)
      wi(2)=qli-i
      wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
      wi(1)=1.d0-wi(2)+wi(3)
      wi(2)=wi(2)-2.d0*wi(3)
      do k1=1,3
       k2=k+k1-1+ml
      do i1=1,3
       qgbit=qgbit+csj(i+i1-1,k2)*wi(i1)*wk(k1)
      enddo
      enddo
      qgbit=exp(qgbit)*(1.d0/tmin-2.d0/s)
      if(qi.lt.qq)qgbit=qgbit*qgsudx(qq,m)/qgsudx(qi,m)
      if(qj.lt.qq)qgbit=qgbit*qgsudx(qq,l)/qgsudx(qj,l)

      if(debug.ge.3)write (moniou,202)qgbit
201   format(2x,'qgbit: qi=',e10.3,2x,'qj=',e10.3
     *,2x,'s= ',e10.3,2x,'m= ',i1,2x,'l= ',i1)
202   format(2x,'qgbit=',e10.3)
      return
      end

c=============================================================================
      double precision function qgfbor(s,t,iq1,iq2,n)
c---------------------------------------------------------------------------
c qgfbor - integrand for the born cross-section (matrix element squared)
c s - total c.m. energy squared for the scattering,
c t - invariant variable for the scattering abs[(p1-p3)**2],
c iq1 - parton type at current end of the ladder (0 - g, 1,2 - q)
c iq2 - parton type at opposite end of the ladder (0 - g, 1,2 - q)
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)s,t,iq1,iq2

      u=s-t
      if(n.eq.1)then
       if(iq1.eq.0.and.iq2.eq.0)then        !gluon-gluon
        qgfbor=(3.d0-t*u/s**2+s*u/t**2+s*t/u**2)*4.5d0
       elseif(iq1*iq2.eq.0)then             !gluon-quark
        qgfbor=(s**2+u**2)/t**2+(s/u+u/s)/2.25d0
       elseif(iq1.eq.iq2)then               !quark-quark (same flavor)
        qgfbor=((s**2+u**2)/t**2+(s**2+t**2)/u**2)/2.25d0
     *  -s**2/t/u/3.375d0
       elseif(iq1+iq2.eq.0)then             !quark-antiquark (same flavor)
        qgfbor=((s**2+u**2)/t**2+(u**2+t**2)/s**2)/2.25d0
     *  +u**2/t/s/3.375d0
       else                                 !quark-antiquark (different flavors)
        qgfbor=(s**2+u**2)/t**2/2.25d0
       endif
      elseif(n.eq.2)then
       if(iq1.eq.0.and.iq2.eq.0)then        !gluon-gluon->quark-antiquark
        qgfbor=.5d0*(t/u+u/t)-1.125d0*(t*t+u*u)/s**2
       elseif(iq1+iq2.eq.0)then             !quark-antiquark->quark-antiquark
        qgfbor=(t*t+u*u)/s**2/1.125d0       !(different flavor)
       else
        qgfbor=0.d0
       endif
      elseif(n.eq.3)then
       if(iq1.ne.0.and.iq1+iq2.eq.0)then    !quark-antiquark->gluon-gluon
        qgfbor=32.d0/27.d0*(t/u+u/t)-(t*t+u*u)/s**2/.375d0
       else
        qgfbor=0.d0
       endif
      endif

      if(debug.ge.2)write (moniou,202)qgfbor
201   format(2x,'qgfbor - hard scattering matrix element squared:'/
     *4x,'s=',e10.3,2x,'|t|=',e10.3,2x,'iq1=',i1,2x,'iq2=',i1)
202   format(2x,'qgfbor=',e10.3)
      return
      end

c=============================================================================
      double precision function qgborn(qi,qj,s,iq1,iq2)
c-----------------------------------------------------------------------------
c qgborn - hard 2->2 parton scattering born cross-section
c s is the c.m. energy square for the scattering process,
c iq1 - parton type at current end of the ladder (0 - g, 1,2 etc. - q)
c iq2 - parton type at opposite end of the ladder (0 - g, 1,2 etc. - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)qi,qj,s,iq1,iq2

      qgborn=0.d0
      qq=max(qi,qj)
      tmin=qq*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qq*fqscal/s)))
      mm=iq1
      ll=iq2
      if(mm*ll.ne.0)then
       if(max(mm,ll).eq.1)then
        mm=1
        ll=2
       elseif(max(mm,ll).eq.2)then
        mm=1
        ll=1
       elseif(max(mm,ll).eq.3)then
        mm=1
        ll=-1
       else
        stop'wrong parton types in qgborn'
       endif
      endif
      do i=1,7
      do m=1,2
       t=2.d0*tmin/(1.d0+2.d0*tmin/s-x1(i)*(2*m-3)*(1.d0-2.d0*tmin/s))
       qt=t*(1.d0-t/s)
       fb=0.d0
       do n=1,3
        fb=fb+qgfbor(s,t,mm,ll,n)+qgfbor(s,s-t,ll,mm,n)
       enddo
       fb=fb*qgsudx(qt/fqscal,min(2,iabs(iq1)+1))
     * *qgsudx(qt/fqscal,min(2,iabs(iq2)+1))
       qgborn=qgborn+a1(i)*fb*qgalf(qt/fqscal/alm)**2*t**2
      enddo
      enddo
      qgborn=qgborn*2.d0*pi**3/s**2
      qgborn=qgborn/qgsudx(qi,min(2,iabs(iq1)+1))
     */qgsudx(qj,min(2,iabs(iq2)+1))
      if(mm.eq.ll)qgborn=qgborn*.5d0

      if(debug.ge.3)write (moniou,202)qgborn
201   format(2x,'qgborn: qi=',e10.3,2x,'qj=',e10.3,2x,
     *'s= ',e10.3,2x,'iq1= ',i1,2x,'iq2= ',i1)
202   format(2x,'qgborn=',e10.3)
      return
      end

c=============================================================================
      subroutine qgcjet(qq,iq1,qv,zv,qm,iqv,ldau,lpar,jq)
c-----------------------------------------------------------------------------
c final state emission process (all branchings as well as parton masses
c are determined)
c qq - maximal effective momentum transfer for the first branching
c iq1 - initial jet flavour (0 - for gluon)
c qv(i,j) - effective momentum for the branching of the parton in i-th row
c on j-th level (0 - in case of no branching)  - to be determined
c zv(i,j) - z-value for the branching of the parton in i-th row
c on j-th level - to be determined
c qm(i,j) - mass squared for the parton in i-th row
c on j-th level - to be determined
c iqv(i,j) - flavour for the parton in i-th row on j-th level
c - to be determined
c ldau(i,j) - first daughter row for the branching of the parton in i-th row
c on j-th level - to be determined
c lpar(i,j) - the parent row for the parton in i-th row
c on j-th level - to be determined
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension qmax(30,50),iqm(2),lnv(50),
     *qv(30,50),zv(30,50),qm(30,50),iqv(30,50),
     *ldau(30,49),lpar(30,50)
      parameter(nfock=3)
      common /qgarr11/ b10
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      if(debug.ge.2)write (moniou,201)qq,iq1,jq

      do i=2,50
       lnv(i)=0
      enddo
      lnv(1)=1
      qmax(1,1)=qq
      iqv(1,1)=iq1
      nlev=1
      nrow=1

2     qlmax=dlog(qmax(nrow,nlev)/qtf/16.d0)
      iq=min(1,iabs(iqv(nrow,nlev)))+1

      if(qgran(b10).gt.qgsudi(qlmax,iq))then   !parton branching
       q=qgqint(qlmax,qgran(b10),iq)
       z=qgzsim(q,iq)
       ll=lnv(nlev+1)+1
       ldau(nrow,nlev)=ll
       lpar(ll,nlev+1)=nrow
       lpar(ll+1,nlev+1)=nrow
       lnv(nlev+1)=ll+1

       if(iq.ne.1)then
        if((3-2*jq)*iqv(nrow,nlev).gt.0)then
         iqm(1)=0
         iqm(2)=iqv(nrow,nlev)
        else
         iqm(2)=0
         iqm(1)=iqv(nrow,nlev)
         z=1.d0-z
        endif
       else
        wg=qgfap(z,1,1)
        wg=wg/(wg+qgfap(z,1,2))
        if(qgran(b10).lt.wg)then
         iqm(1)=0
         iqm(2)=0
        else
         iqm(1)=int(3.d0*qgran(b10)+1.d0)*(3-2*jq)
         iqm(2)=-iqm(1)
        endif
        if(qgran(b10).lt..5d0)z=1.d0-z
       endif
       qv(nrow,nlev)=q
       zv(nrow,nlev)=z
       nrow=ll
       nlev=nlev+1
       qmax(nrow,nlev)=q*z**2
       qmax(nrow+1,nlev)=q*(1.d0-z)**2
       iqv(nrow,nlev)=iqm(1)
       iqv(nrow+1,nlev)=iqm(2)
       if(debug.ge.3)write (moniou,203)nlev,nrow,q,z
       goto 2
      else
       qv(nrow,nlev)=0.d0
       zv(nrow,nlev)=0.d0
       qm(nrow,nlev)=0.d0
       if(debug.ge.3)write (moniou,204)nlev,nrow
      endif

3     continue
      if(nlev.eq.1)then
       if(debug.ge.3)write (moniou,202)
       return
      endif

      lprow=lpar(nrow,nlev)
      if(ldau(lprow,nlev-1).eq.nrow)then
       nrow=nrow+1
       goto 2
      else
       z=zv(lprow,nlev-1)
       qm(lprow,nlev-1)=z*(1.d0-z)*qv(lprow,nlev-1)
     * +qm(nrow-1,nlev)/z+qm(nrow,nlev)/(1.d0-z)
       nrow=lprow
       nlev=nlev-1
       if(debug.ge.3)write (moniou,205)nlev,nrow,qm(lprow,nlev)
       goto 3
      endif

201   format(2x,'qgcjet: qq=',e10.3,2x,'iq1= ',i1,2x,'jq=',i1)
202   format(2x,'qgcjet - end')
203   format(2x,'qgcjet: new branching at level nlev=',i2,' nrow=',i2
     */4x,' effective momentum q=',e10.3,2x,' z=',e10.3)
204   format(2x,'qgcjet: new final jet at level nlev=',i2,' nrow=',i2)
205   format(2x,'qgcjet: jet mass at level nlev=',i2,' nrow=',i2
     *,' - qm=',e10.3)
      end

c===========================================================================
      subroutine qgcs(c,s)
c---------------------------------------------------------------------------
c c,s - cos and sin generation for uniformly distributed angle 0<fi<2*pi
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      if(debug.ge.3)write (moniou,201)
1     s1=2.d0*qgran(b10)-1.d0
      s2=2.d0*qgran(b10)-1.d0
      s3=s1*s1+s2*s2
      if(s3.gt.1.d0)goto 1
      s3=dsqrt(s3)
      c=s1/s3
      s=s2/s3

      if(debug.ge.3)write (moniou,202)c,s
201   format(2x,'qgcs - cos(fi) and sin(fi) are generated',
     *' (0<fi<2*pi)')
202   format(2x,'qgcs: c=',e10.3,2x,'s=',e10.3)
      return
      end

c===========================================================================
      subroutine qgdeft(s,ep,ey)
c---------------------------------------------------------------------------
c determination of the parameters for the lorentz transform to the rest frame
c system for 4-vector ep
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ey(3),ep(4)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)ep,s

      do i=1,3
       if(ep(i+1).eq.0.d0)then
        ey(i)=1.d0
       else
        wp=ep(1)+ep(i+1)
        wm=ep(1)-ep(i+1)
        if(wm/wp.lt.1.d-8)then
         ww=s
         do l=1,3
          if(l.ne.i)ww=ww+ep(l+1)**2
         enddo
         wm=ww/wp
        endif
        ey(i)=dsqrt(wm/wp)
        ep(1)=wp*ey(i)
        ep(i+1)=0.d0
       endif
      enddo

      if(debug.ge.3)write (moniou,202)ey
201   format(2x,'qgdeft - lorentz boost parameters:'
     */4x,'4-vector ep=',4e10.3/4x,'4-vector squared s=',e10.3)
202   format(2x,'qgdeft: lorentz boost parameters ey(i)=',2x,3e10.3)
      return
      end

c=============================================================================
      subroutine qgdefr(ep,s0x,c0x,s0,c0)
c-----------------------------------------------------------------------------
c determination of the parameters the spatial rotation to the lab. system
c for 4-vector ep
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ep(4)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)ep

c transverse momentum square for the current parton (ep)
      pt2=ep(3)**2+ep(4)**2
      if(pt2.ne.0.d0)then
       pt=dsqrt(pt2)
c system rotation to get pt=0 - euler angles are determined (c0x = cos theta,
c s0x = sin theta, c0 = cos phi, s0 = sin phi)
       c0x=ep(3)/pt
       s0x=ep(4)/pt
c total momentum for the gluon
       pl=dsqrt(pt2+ep(2)**2)
       s0=pt/pl
       c0=ep(2)/pl
      else
       c0x=1.d0
       s0x=0.d0
       pl=abs(ep(2))
       s0=0.d0
       c0=ep(2)/pl
      endif
      ep(2)=pl
      ep(3)=0.d0
      ep(4)=0.d0

      if(debug.ge.3)write (moniou,202)s0x,c0x,s0,c0,ep
201   format(2x,'qgdefr - spatial rotation parameters'/4x,
     *'4-vector ep=',2x,4(e10.3,1x))
202   format(2x,'qgdefr: spatial rotation parameters'/
     *4x,'s0x=',e10.3,2x,'c0x=',e10.3,2x,'s0=',e10.3,2x,'c0=',e10.3/
     *4x,'rotated 4-vector ep=',4(e10.3,1x))
      return
      end

c=============================================================================
      double precision function qgfap(x,j,l)
c------------------------------------------------------------------------
c qgfap - altarelli-parisi function (multiplied by x)
c x - light cone momentum share value,
c j - type of the parent parton (1-g,2-q)
c l - type of the daughter parton (1-g,2-q)
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)x,j,l

      if(j.eq.1)then
       if(l.eq.1)then
        qgfap=((1.d0-x)/x+x/(1.d0-x)+x*(1.d0-x))*6.d0
       else
        qgfap=(x**2+(1.d0-x)**2)*3.d0
       endif
      else
       if(l.eq.1)then
        qgfap=(1.d0+(1.d0-x)**2)/x/.75d0
       else
        qgfap=(x**2+1.d0)/(1.d0-x)/.75d0
       endif
      endif

      if(debug.ge.3)write (moniou,202)qgfap
201   format(2x,'qgfap - altarelli-parisi function:'
     *,2x,'x=',e10.3,2x,'j=',i1,2x,'l=',i1)
202   format(2x,'qgfap=',e10.3)
      return
      end

c=============================================================================
      subroutine qggea(ia,xa,jj)
c-----------------------------------------------------------------------------
c qggea - nuclear configuration simulation (nucleons positions)
c ia - number of nucleons to be considered
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208)
      dimension xa(iapmax,3)
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      if(debug.ge.2)write (moniou,201)jj,ia

      if(ia.ge.10)then
       do i=1,ia
1       zuk=qgran(b10)*cr1(jj)-1.d0
        if(zuk.le.0.d0)then
         tt=rnuc(jj)/wsnuc(jj)*(qgran(b10)**.3333d0-1.d0)
        elseif(zuk.le.cr2(jj))then
         tt=-dlog(qgran(b10))
        elseif(zuk.le.cr3(jj))then
         tt=-dlog(qgran(b10))-dlog(qgran(b10))
        else
         tt=-dlog(qgran(b10))-dlog(qgran(b10))-dlog(qgran(b10))
        endif
        rim=tt*wsnuc(jj)+rnuc(jj)
        
        if(qgran(b10).gt.(1.d0+wbnuc(jj)*rim**2/rnuc(jj)**2)
     *  /(1.d0+5.d0*abs(wbnuc(jj)))/(1.d0+exp(-abs(tt))))goto 1
       z=rim*(2.d0*qgran(b10)-1.d0)
        rim=dsqrt(rim*rim-z*z)
        xa(i,3)=z
        call qgcs(c,s)
        xa(i,1)=rim*c
        xa(i,2)=rim*s
       enddo
      else !Gaussian distr. (conserving center of mass)
       do l=1,3
        summ=0.d0
        do i=1,ia-1
         j=ia-i
         aks=rnuc(jj)*(qgran(b10)+qgran(b10)+qgran(b10)-1.5d0)
         k=j+1
         xa(k,l)=summ-aks*sqrt(float(j)/k)
         summ=summ+aks/sqrt(float(j*k))
        enddo
        xa(1,l)=summ
       enddo
      endif

      if(debug.ge.3)then
       write (moniou,203)
       do i=1,ia
        write (moniou,204)i,(xa(i,l),l=1,3)
       enddo
       write (moniou,202)
      endif
201   format(2x,'qggea - configuration of the nucleus ',i1,';',2x,
     *'coordinates for ',i2,' nucleons')
202   format(2x,'qggea - end')
203   format(2x,'qggea:  positions of the nucleons')
204   format(2x,'qggea: ',i2,' - ',3(e10.3,1x))
      return
      end

c=============================================================================
      double precision function qgapi(x,j,l)
c-----------------------------------------------------------------------------
c qgapi - integrated altarelli-parisi function
c x - light cone momentum share value,
c j - type of initial parton (1 - g, 2 - q)
c l - type of final parton (1 - g, 2 - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)x,j,l

      if(j.eq.1)then
       if(l.eq.1)then
        qgapi=6.d0*(dlog(x/(1.d0-x))-x**3/3.d0+x**2/2.d0-2.d0*x)
       else
        qgapi=3.d0*(x+x**3/1.5d0-x*x)
       endif
      else
       if(l.eq.1)then
        qgapi=(dlog(x)-x+.25d0*x*x)/.375d0
       else
        z=1.d0-x
        qgapi=-(dlog(z)-z+.25d0*z*z)/.375d0
       endif
      endif

      if(debug.ge.2)write (moniou,202)qgapi
201   format(2x,'qgapi: x=',e10.3,2x,'j= ',i1,2x,'l= ',i1)
202   format(2x,'qgapi=',e10.3)
      return
      end

c=============================================================================
      subroutine qgjarr(jfl)
c-----------------------------------------------------------------------------
c final jets rearrangement according to their colour connections
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(njmax=50000)
      dimension mark(njmax),ept(4)
      common /qgarr10/ am(6)
      common /qgarr36/ epjet(4,njmax),ipjet(njmax),njtot
      common /qgarr37/ eqj(4,njmax),iqj(njmax),ncj(2,njmax),nj
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)nj
      if(debug.ge.2.and.nj.ne.0)then
       do i=1,nj
        write (moniou,203)i,iqj(i),(eqj(l,i),l=1,4)
        if(iqj(i).eq.0)then
         write (moniou,204)ncj(1,i),ncj(2,i)
        else
         ncdum=0
         write (moniou,204)ncj(1,i),ncdum
        endif
       enddo
      endif

      njpar=0
      jfl=0
      do i=1,nj
       mark(i)=1
      enddo
      njtot=0

1     continue
      do ij=1,nj
       if(mark(ij).ne.0.and.iqj(ij).ne.0)goto 2
      enddo
      stop'problem in qgjarr?'
2     continue

      jfirst=1
      if(iabs(iqj(ij)).le.2)then
       am1=am(1)
      elseif(iabs(iqj(ij)).eq.4)then
       am1=am(3)
      else
       am1=am(2)
      endif
      do i=1,4
       ept(i)=0.d0
      enddo

3     mark(ij)=0
      njtot=njtot+1
      ipjet(njtot)=iqj(ij)
      do i=1,4
       ept(i)=ept(i)+eqj(i,ij)
       epjet(i,njtot)=eqj(i,ij)
      enddo

      if(iqj(ij).ne.0)then
       if(jfirst.ne.1)then
        if(iabs(iqj(ij)).le.2)then
         am2=am(1)
        elseif(iabs(iqj(ij)).eq.4)then
         am2=am(3)
        else
         am2=am(2)
        endif
        amj=(am1+am2)**2
        if(amj.gt.qgnrm(ept))then
         if(debug.ge.3)write (moniou,202)jfl
         return
        endif

        if(njtot.lt.nj)then
         goto 1
        else
         jfl=1
         nj=0
         if(debug.ge.3)write (moniou,202)jfl
         return
        endif
       else
        jfirst=0
        njpar=ij
        ij=ncj(1,ij)
        goto 3
       endif
      else
       if(ncj(1,ij).eq.njpar)then
        njdau=ncj(2,ij)
       else
        njdau=ncj(1,ij)
       endif
       njpar=ij
       ij=njdau
       goto 3
      endif

201   format(2x,'qgjarr: total number of jets nj=',i4)
202   format(2x,'qgjarr - end,jfl=',i2)
203   format(2x,'qgjarr: ij=',i3,2x,'iqj=',i2,2x,'eqj=',4e10.3)
204   format(2x,'qgjarr: ncj=',2i3)
      end

c=============================================================================
      double precision function qgjit(q1,q2,s,m,l)
c-----------------------------------------------------------------------------
c qgjit - inclusive hard cross-section interpolation - for any ordering
c in the ladder
c q1 - effective momentum cutoff for current end of the ladder,
c q2 - effective momentum cutoff for opposide end of the ladder,
c s - total c.m. energy squared for the ladder,
c m - parton type at current end of the ladder (1 - g, 2,3,4 - q)
c l - parton type at opposite end of the ladder (1 - g, 2,3,4 - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wi(3),wj(3),wk(3)
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr20/ spmax
      common /qgarr26/ factk,fqscal
      common /qgarr29/ csj(40,40,240)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)s,q1,q2,m,l

      qgjit=0.d0
      qq=max(q1,q2)
      s2min=qq*4.d0*fqscal
      if(s.lt..9999d0*s2min)then
       write (*,202)qgjit
       if(debug.ge.3)write (moniou,202)qgjit
       return
      endif

      mm=max(m,l)
      if(q1.le.q2)then
       qi=q1
       qj=q2
       ml=40*(m-1)+80*(l-1)
      else
       qi=q2
       qj=q1
       ml=40*(l-1)+80*(m-1)
      endif
      if(mm.gt.2)ml=40*(mm+1)

      tmin=qq*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qq*fqscal/s)))
      qli=dlog(qi)/dlog(spmax/4.d0/fqscal)*39.d0+1.d0
      if(qi.lt..99d0*spmax/4.d0/fqscal)then
       qlj=dlog(qj/qi)/dlog(spmax/4.d0/fqscal/qi)*39.d0+1.d0
      else
       qlj=1.d0
      endif
      sl=dlog(s/s2min)/dlog(spmax/s2min)*39.d0+1.d0
      i=min(38,int(qli))
      j=min(38,int(qlj))
      k=min(38,int(sl))

      wk(2)=sl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)
      wi(2)=qli-i
      wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
      wi(1)=1.d0-wi(2)+wi(3)
      wi(2)=wi(2)-2.d0*wi(3)
      wj(2)=qlj-j
      wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
      wj(1)=1.d0-wj(2)+wj(3)
      wj(2)=wj(2)-2.d0*wj(3)
      do k1=1,3
       k2=k+k1-1+ml
      do i1=1,3
      do j1=1,3
       qgjit=qgjit+csj(i+i1-1,j+j1-1,k2)*wi(i1)*wj(j1)*wk(k1)
      enddo
      enddo
      enddo
      qgjit=exp(qgjit)*(1.d0/tmin-2.d0/s)

      if(debug.ge.3)write (moniou,202)qgjit
201   format(2x,'qgjit - unordered ladder cross section interpol.:'/4x,
     *'s=',e10.3,2x,'q1=',e10.3,2x,'q2=',e10.3,2x,2x,'m=',i1,2x,'l=',i1)
202   format(2x,'qgjit=',e10.3)
      return
      end

c=============================================================================
      double precision function qgjit1(q1,q2,s,m,l)
c-----------------------------------------------------------------------------
c qgjit1 - inclusive hard cross-section interpolation - for strict ordering
c in the ladder
c q1 - effective momentum cutoff for current end of the ladder,
c q2 - effective momentum cutoff for opposide end of the ladder,
c s - total c.m. energy squared for the ladder,
c m - parton type at current end of the ladder (1 - g, 2,3,4 - q)
c l - parton type at opposite end of the ladder (1 - g, 2,3,4 - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wi(3),wj(3),wk(3)
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr20/ spmax
      common /qgarr26/ factk,fqscal
      common /qgarr30/ csj(40,40,240)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)s,q1,q2,m,l

      qgjit1=0.d0
      qq=max(q1,q2)
      s2min=qq*4.d0*fqscal
      if(s.le.s2min)then
       if(debug.ge.3)write (moniou,202)qgjit1
       return
      endif

      tmin=qq*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qq*fqscal/s)))
      qli=dlog(q1)/dlog(spmax/4.d0/fqscal)*39.d0+1.d0
      if(q1.lt..99d0*spmax/4.d0/fqscal)then
       qlj=dlog(qq/q1)/dlog(spmax/4.d0/fqscal/q1)*39.d0+1.d0
      else
       qlj=1.d0
      endif
      sl=dlog(s/s2min)/dlog(spmax/s2min)*39.d0+1.d0
      i=min(38,int(qli))
      j=min(38,int(qlj))
      k=min(38,int(sl))
      mm=max(m,l)
      if(mm.le.2)then
       ml=40*(m-1)+80*(l-1)
      else
       ml=40*(mm+1)
      endif

      wk(2)=sl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)
      wi(2)=qli-i
      wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
      wi(1)=1.d0-wi(2)+wi(3)
      wi(2)=wi(2)-2.d0*wi(3)
      wj(2)=qlj-j
      wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
      wj(1)=1.d0-wj(2)+wj(3)
      wj(2)=wj(2)-2.d0*wj(3)

      do k1=1,3
       k2=k+k1-1+ml
      do i1=1,3
      do j1=1,3
       qgjit1=qgjit1+csj(i+i1-1,j+j1-1,k2)*wi(i1)*wj(j1)*wk(k1)
      enddo
      enddo
      enddo
      qgjit1=exp(qgjit1)*(1.d0/tmin-2.d0/s)
      if(q2.lt.q1)qgjit1=qgjit1*qgsudx(q1,min(2,l))/qgsudx(q2,min(2,l))

      if(debug.ge.3)write (moniou,202)qgjit1
201   format(2x,'qgjit1 - ordered ladder cross section interpol.:'/4x,
     *'s=',e10.3,2x,'q1=',e10.3,2x,'q2=',e10.3,2x,2x,'m=',i1,2x,'l=',i1)
202   format(2x,'qgjit1=',e10.3)
      return
      end

c=============================================================================
      double precision function qglam(s,a,b)
c-----------------------------------------------------------------------------
c kinematical function for two particle decay - maximal pt-value
c a - first particle mass squared,
c b - second particle mass squared,
c s - two particle invariant mass
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)s,a,b

      s1=dsqrt(s)
      a1=dsqrt(a)
      b1=dsqrt(b)
      qglam=max(0.d0,s1-a1-b1)*(s1-a1+b1)*(s1+a1-b1)*(s1+a1+b1)/4.d0/s

      if(debug.ge.3)write (moniou,202)qglam
201   format(2x,'qglam - kinematical function s=',e10.3,2x,'a='
     *,e10.3,2x,'b=',e10.3)
202   format(2x,'qglam=',e10.3)
      return
      end

c=============================================================================
      double precision function qgnrm(ep)
c-----------------------------------------------------------------------------
c 4-vector squared calculation
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ep(4)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)ep
      qgnrm=(ep(1)-ep(2))*(ep(1)+ep(2))-ep(3)**2-ep(4)**2

      if(debug.ge.3)write (moniou,202)qgnrm
201   format(2x,'qgnrm - 4-vector squared for ','ep=',4(e10.3,1x))
202   format(2x,'qgnrm=',e10.3)
      return
      end

c===========================================================================
      subroutine qgrec(ep,nqc,qv,zv,qm,iqv,ldau,lpar,jq)
c---------------------------------------------------------------------------
c jet reconstructuring procedure - 4-momenta for all final jets are determ.
c ep(i) - jet 4-momentum
c---------------------------------------------------------------------------
c qv(i,j) - effective momentum for the branching of the parton in i-th row
c on j-th level (0 - in case of no branching)
c zv(i,j) - z-value for the branching of the parton in i-th row
c on j-th level
c qm(i,j) - mass squared for the parton in i-th row
c on j-th level
c iqv(i,j) - flavours for the parton in i-th row on j-th level
c ldau(i,j) - first daughter row for the branching of the parton in i-th row
c on j-th level
c lpar(i,j) - the parent row for the parton in i-th row on j-th level
c----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(njmax=50000)
      dimension ep(4),ep3(4),epv(4,30,50),nqc(2),ncc(2,30,50),
     *qv(30,50),zv(30,50),qm(30,50),iqv(30,50),
     *ldau(30,49),lpar(30,50)
c eqj(i,nj) - 4-momentum for the final jet nj
c iqj(nj) - flavour for the final jet nj
c ncj(m,nj) - colour connections for the final jet nj
      common /qgarr37/ eqj(4,njmax),iqj(njmax),ncj(2,njmax),nj
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)jq,ep,iqv(1,1),nqc

      do i=1,4
       epv(i,1,1)=ep(i)
      enddo
      ncc(1,1,1)=nqc(1)
      if(iqv(1,1).eq.0)ncc(2,1,1)=nqc(2)
      nlev=1
      nrow=1

2     continue
      if(qv(nrow,nlev).eq.0.d0)then            !final parton
       nj=nj+1
       do i=1,4
        eqj(i,nj)=epv(i,nrow,nlev)
       enddo
       iqj(nj)=iqv(nrow,nlev)
       if(iabs(iqj(nj)).eq.3)iqj(nj)=iqj(nj)*4/3

       if(iqj(nj).ne.0)then                    !(anti)quark
        njc=ncc(1,nrow,nlev)
        if(njc.ne.0)then
         ncj(1,nj)=njc
         iqc=iqj(njc)
         if(iqc.ne.0)then
          ncj(1,njc)=nj
         else
          if(iqj(nj).gt.0)then
           ncj(2,njc)=nj
          else
           ncj(1,njc)=nj
          endif
         endif
        else
         ncc(1,nrow,nlev)=nj
        endif
        
       else                                    !gluon
        do m=1,2
         if(jq.eq.1)then
          m1=m
         else
          m1=3-m
         endif
         njc=ncc(m1,nrow,nlev)
         if(njc.ne.0)then
          ncj(m,nj)=njc
          iqc=iqj(njc)
          if(iqc.ne.0)then
           ncj(1,njc)=nj
          else
           ncj(3-m,njc)=nj
          endif
         else
          ncc(m1,nrow,nlev)=nj
         endif
        enddo
       endif
       if(debug.ge.3)write (moniou,204)
     * nj,nlev,nrow,iqj(nj),(eqj(i,nj),i=1,4)

      else
       do i=1,4
        ep3(i)=epv(i,nrow,nlev)
       enddo
       call qgdefr(ep3,s0x,c0x,s0,c0)
       z=zv(nrow,nlev)
       qt2=(z*(1.d0-z))**2*qv(nrow,nlev)
       ldrow=ldau(nrow,nlev)

       wp0=ep3(1)+ep3(2)
       wpi=z*wp0
       wmi=(qt2+qm(ldrow,nlev+1))/wpi
       ep3(1)=.5d0*(wpi+wmi)
       ep3(2)=.5d0*(wpi-wmi)
       qt=dsqrt(qt2)
       call qgcs(c,s)
       ep3(3)=qt*c
       ep3(4)=qt*s
       call qgrota(ep3,s0x,c0x,s0,c0)
       do i=1,4
        epv(i,ldrow,nlev+1)=ep3(i)
       enddo
       if(debug.ge.3)write (moniou,206)nlev+1,ldrow,ep3

       wpi=(1.d0-z)*wp0
       wmi=(qt2+qm(ldrow+1,nlev+1))/wpi
       ep3(1)=.5d0*(wpi+wmi)
       ep3(2)=.5d0*(wpi-wmi)
       ep3(3)=-qt*c
       ep3(4)=-qt*s
       call qgrota(ep3,s0x,c0x,s0,c0)
       do i=1,4
        epv(i,ldrow+1,nlev+1)=ep3(i)
       enddo
       if(debug.ge.3)write (moniou,206)nlev+1,ldrow+1,ep3

       if(iqv(nrow,nlev).eq.0)then
        if(iqv(ldrow,nlev+1).ne.0)then         !g->qq~
         ncc(1,ldrow,nlev+1)=ncc(1,nrow,nlev)
         ncc(1,ldrow+1,nlev+1)=ncc(2,nrow,nlev)
        else                                   !g->gg
         ncc(1,ldrow,nlev+1)=ncc(1,nrow,nlev)
         ncc(2,ldrow,nlev+1)=0
         ncc(1,ldrow+1,nlev+1)=0
         ncc(2,ldrow+1,nlev+1)=ncc(2,nrow,nlev)
        endif
       else
        if(iqv(ldrow,nlev+1).eq.0)then
         ncc(1,ldrow,nlev+1)=ncc(1,nrow,nlev)
         ncc(2,ldrow,nlev+1)=0
         ncc(1,ldrow+1,nlev+1)=0
        else
         ncc(1,ldrow,nlev+1)=0
         ncc(1,ldrow+1,nlev+1)=0
         ncc(2,ldrow+1,nlev+1)=ncc(1,nrow,nlev)
        endif
       endif

       nrow=ldrow
       nlev=nlev+1
       goto 2
      endif

8     continue
      if(nlev.eq.1)then
       if(nqc(1).eq.0)nqc(1)=ncc(1,1,1)
       if(iqv(1,1).eq.0.and.nqc(2).eq.0)nqc(2)=ncc(2,1,1)
       if(debug.ge.3)write (moniou,202)
       return
      endif

      lprow=lpar(nrow,nlev)
      if(ldau(lprow,nlev-1).eq.nrow)then
       if(iqv(nrow,nlev).eq.0)then
        if(ncc(1,lprow,nlev-1).eq.0)ncc(1,lprow,nlev-1)=ncc(1,nrow,nlev)
        ncc(1,nrow+1,nlev)=ncc(2,nrow,nlev)
       else
        if(iqv(lprow,nlev-1).eq.0)then
         if(ncc(1,lprow,nlev-1).eq.0)
     *   ncc(1,lprow,nlev-1)=ncc(1,nrow,nlev)
        else
         ncc(1,nrow+1,nlev)=ncc(1,nrow,nlev)
        endif
       endif
       nrow=nrow+1
       goto 2
      else
       if(iqv(nrow,nlev).eq.0)then
        if(iqv(lprow,nlev-1).eq.0)then
         if(ncc(2,lprow,nlev-1).eq.0)
     *   ncc(2,lprow,nlev-1)=ncc(2,nrow,nlev)
        else
         if(ncc(1,lprow,nlev-1).eq.0)
     *   ncc(1,lprow,nlev-1)=ncc(2,nrow,nlev)
        endif
       else
        if(iqv(lprow,nlev-1).eq.0.and.ncc(2,lprow,nlev-1).eq.0)
     *  ncc(2,lprow,nlev-1)=ncc(1,nrow,nlev)
       endif
       nrow=lprow
       nlev=nlev-1
       goto 8
      endif

201   format(2x,'qgrec - jet reconstructuring: jq=',i1
     */4x,'jet 4-momentum ep=',4(e10.3,1x)
     */4x,'jet flavor: ',i2,2x,'colour connections: ',2i3)
202   format(2x,'qgrec - end')
204   format(2x,'qgrec: ',i3,'-th final jet at level nlev=',i2,' nrow='
     *,i2/4x,'jet flavor: ',i3,2x,'jet 4-momentum:',4(e10.3,1x))
206   format(2x,'qgrec: jet at level nlev='
     *,i2,' nrow=',i2/4x,'jet 4-momentum:',4(e10.3,1x))
      end

c=============================================================================
      double precision function qgroot(qlmax,g,j)
c-----------------------------------------------------------------------------
c qgroot - effective momentum tabulation for given set of random number
c values and maximal effective momentum qmax values - according to the
c probability of branching: (1 - timelike sudakov formfactor)
c qlmax - ln qmax/16/qtf,
c g - dzeta number (some function of ksi)
c j - type of the parton (1-g,2-q)
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)qlmax,g,j

      ql0=0.d0
      ql1=qlmax
      f0=-g
      f1=1.d0-g
      sud0=-dlog(qgsudi(qlmax,j))

1     ql2=ql1-(ql1-ql0)*f1/(f1-f0)
      if(ql2.lt.0.d0)then
       ql2=0.d0
       f2=-g
      elseif(ql2.gt.qlmax)then
       ql2=qlmax
       f2=1.d0-g
      else
       f2=-dlog(qgsudi(ql2,j))/sud0-g
      endif
      if(abs(f2).gt.1.d-3)then
       if(f1*f2.lt.0.d0)then
        ql0=ql1
        f0=f1
       endif
       ql1=ql2
       f1=f2
       goto 1
      else
       qgroot=ql2
      endif

      if(debug.ge.3)write (moniou,202)qgroot
201   format(2x,'qgqint - branching momentum tabulation:'
     */4x,'qlmax=',e10.3,2x,'g=',e10.3,2x,'j=',i1)
202   format(2x,'qgroot=',e10.3)
      return
      end

c=============================================================================
      subroutine qgrota(ep,s0x,c0x,s0,c0)
c-----------------------------------------------------------------------------
c spatial rotation to the lab. system for 4-vector ep
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ep(4),ep1(3)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)ep,s0x,c0x,s0,c0

      ep1(3)=ep(4)
      ep1(2)=ep(2)*s0+ep(3)*c0
      ep1(1)=ep(2)*c0-ep(3)*s0
      ep(2)=ep1(1)
      ep(4)=ep1(2)*s0x+ep1(3)*c0x
      ep(3)=ep1(2)*c0x-ep1(3)*s0x

      if(debug.ge.3)write (moniou,202)ep
201   format(2x,'qgrota - spatial rotation:'/4x,'4-vector ep=',4(e10.3
     *,1x)/4x,'s0x=',e10.3,'c0x=',e10.3,2x,'s0=',e10.3,'c0=',e10.3)
202   format(2x,'qgrota: rotated 4-vector ep=',2x,4e10.3)
      return
      end

c=============================================================================
      double precision function qgqint(qlmax,g,j)
c-----------------------------------------------------------------------------
c qgqint - effective momentum interpolation for given random number g
c and maximal effective momentum qmax
c qlmax - ln qmax/16/qtf,
c g - random number (0<g<1)
c j - type of the parton (1-g,2-q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wi(3),wk(3)
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr34/ qrt(10,101,2)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)qlmax,g,j

      qli=qlmax/1.38629d0
      sud0=1.d0/qgsudi(qlmax,j)
      sl=100.d0*dlog(1.d0-g*(1.d0-sud0))/dlog(sud0)
      i=int(qli)
      k=int(sl)
      if(k.gt.98)k=98
      wk(2)=sl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)
      qgqint=0.d0
      if(i.gt.7)i=7
      wi(2)=qli-i
      wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
      wi(1)=1.d0-wi(2)+wi(3)
      wi(2)=wi(2)-2.d0*wi(3)
      do k1=1,3
      do i1=1,3
       qgqint=qgqint+qrt(i+i1,k+k1,j)*wi(i1)*wk(k1)
      enddo
      enddo
      if(qgqint.le.0.d0)qgqint=0.d0
      qgqint=16.d0*qtf*exp(qgqint)

      if(debug.ge.3)write (moniou,202)qgqint
201   format(2x,'qgqint - branching momentum interpolation:'
     */4x,'qlmax=',e10.3,2x,'g=',e10.3,2x,'j=',i1)
202   format(2x,'qgqint=',e10.3)
      return
      end

c=============================================================================
      double precision function qgalf(qq)
c-----------------------------------------------------------------------------
c qgalf - alpha_s(qq)/2/pi
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr43/ moniou
      common /qgdebug/ debug

      qgalf=2.d0/9.d0/dlog(qq)
      return
      end

c=============================================================================
      subroutine qgtran(ep,ey,jj)
c-----------------------------------------------------------------------------
c lorentz transform according to parameters ey ( determining lorentz shift
c along the z,x,y-axis respectively (ey(1),ey(2),ey(3)))
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ey(3),ep(4)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)ep,ey

      if(jj.eq.1)then
c lorentz boost to lab. frame according to 1/ey(i) parameters
       do i=1,3
        if(ey(4-i).ne.1.d0)then
         wp=(ep(1)+ep(5-i))/ey(4-i)
         wm=(ep(1)-ep(5-i))*ey(4-i)
         ep(1)=.5d0*(wp+wm)
         ep(5-i)=.5d0*(wp-wm)
        endif
       enddo
      else
c lorentz boost to c.m. frame according to ey(i) parameters
       do i=1,3
        if(ey(i).ne.1.d0)then
         wp=(ep(1)+ep(i+1))*ey(i)
         wm=(ep(1)-ep(i+1))/ey(i)
         ep(1)=.5d0*(wp+wm)
         ep(i+1)=.5d0*(wp-wm)
        endif
       enddo
      endif

      if(debug.ge.3)write (moniou,202)ep
201   format(2x,'qgtran - lorentz boost for 4-vector'/4x,'ep='
     *,2x,4(e10.3,1x)/4x,'boost parameters ey=',3e10.3)
202   format(2x,'qgtran: transformed 4-vector ep=',2x,4(e10.3,1x))
      return
      end

c=============================================================================
      double precision function qgsudi(qlmax,j)
c-----------------------------------------------------------------------------
c qgsudi - timelike sudakov formfactor interpolation
c qlmax - ln qmax/16/qtf,
c j - type of the parton (1-g,2-q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3)
      common /qgarr33/ fsud(10,2)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)j,qlmax

      ql=qlmax/1.38629d0
      if(ql.le.0.d0)then
       qgsudi=1.d0
      else
       k=int(ql)
       if(k.gt.7)k=7
       wk(2)=ql-k
       wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
       wk(1)=1.d0-wk(2)+wk(3)
       wk(2)=wk(2)-2.d0*wk(3)

       qgsudi=0.d0
       do k1=1,3
        qgsudi=qgsudi+fsud(k+k1,j)*wk(k1)
       enddo
       if(qgsudi.le.0.d0)qgsudi=0.d0
       qgsudi=exp(-qgsudi)
      endif

      if(debug.ge.3)write (moniou,202)qgsudi
201   format(2x,'qgsudi - spacelike form factor interpolation:'
     */4x,'parton type j=',i1,2x,'momentum logarithm qlmax=',e10.3)
202   format(2x,'qgsudi=',e10.3)
      return
      end

c=============================================================================
      double precision function qgsudx(q,j)
c-----------------------------------------------------------------------------
c qgsudx - spacelike sudakov formfactor
c q - maximal value of the effective momentum,
c j - type of parton (1 - g, 2 - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr43/ moniou
      common /qgarr51/ epsxmn
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)j,q

      if(q.gt.1.d0)then
       qgsudx=dlog(dlog(q/alm)/dlog(1.d0/alm))*(.75d0+dlog(epsxmn))
       if(j.eq.1)then
        qgsudx=exp(qgsudx/.75d0)
       else
        qgsudx=exp(qgsudx*16.d0/27.d0)
       endif
      else
       qgsudx=1.d0
      endif

      if(debug.ge.3)write (moniou,202)qgsudx
201   format(2x,'qgsudx - spacelike form factor: parton type j='
     *,i1,2x,'momentum q=',e10.3)
202   format(2x,'qgsudx=',e10.3)
      return
      end

c=============================================================================
      double precision function qgsudt(qmax,j)
c-----------------------------------------------------------------------------
c qgsudt - timelike sudakov formfactor
c qmax - maximal value of the effective momentum,
c j - type of parton (1 - g, 2 - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr43/ moniou
      common /arr3/    x1(7),a1(7)
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)j,qmax

      qgsudt=0.d0
      qlmax=dlog(dlog(qmax/16.d0/alm))
      qfl=dlog(dlog(qtf/alm))
      do i=1,7
      do m=1,2
       qtl=.5d0*(qlmax+qfl+(2*m-3)*x1(i)*(qlmax-qfl))
       qt=alm*exp(exp(qtl))
       if(qt.ge.qmax/16.d0)qt=qmax/16.0001d0
       zmin=.5d0-dsqrt((.25d0-dsqrt(qt/qmax)))
       zmax=1.d0-zmin

       if(j.eq.1)then
        ap=(qgapi(zmax,1,1)-qgapi(zmin,1,1)+
     *  qgapi(zmax,1,2)-qgapi(zmin,1,2))*.5d0
       else
        ap=qgapi(zmax,2,1)-qgapi(zmin,2,1)
       endif
       qgsudt=qgsudt+a1(i)*ap
      enddo
      enddo
      qgsudt=qgsudt*(qlmax-qfl)/9.d0

      if(debug.ge.3)write (moniou,202)qgsudt
201   format(2x,'qgsudt - timelike form factor: parton type j='
     *,i1,2x,'momentum qmax=',e10.3)
202   format(2x,'qgsudt=',e10.3)
      return
      end

c=============================================================================
      double precision function qgtwd(s,a,b)
c-----------------------------------------------------------------------------
c kinematical function for two particle decay - light cone momentum share
c for the particle of mass squared a,
c b - partner's mass squared,
c s - two particle invariant mass
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)s,a,b

      x=.5d0*(1.d0+(a-b)/s)
      dx=x-dsqrt(a/s)
      if(dx.gt.0.d0)then
       x=x+dsqrt(dx)*dsqrt(x+dsqrt(a/s))
      else
       x=dsqrt(a/s)
      endif
      qgtwd=x

      if(debug.ge.3)write (moniou,202)qgtwd
201   format(2x,'qgtwd: s=',e10.3,2x,'a=',e10.3,2x,'b=',e10.3)
202   format(2x,'qgtwd=',e10.3)
      return
      end

c=============================================================================
      double precision function qgzsim(qq,j)
c-----------------------------------------------------------------------------
c qgzsim - light cone momentum share simulation (for the timelike
c branching)
c qq - effective momentum value,
c j - type of the parent parton (1-g,2-q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr11/ b10
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      if(debug.ge.2)write (moniou,201)qq,j

      zmin=.5d0-dsqrt(.25d0-dsqrt(qtf/qq))
      qlf=dlog(qtf/alm)
1     continue
      if(j.eq.1)then
       qgzsim=.5d0*(2.d0*zmin)**qgran(b10)
       gb=qgzsim*(qgfap(qgzsim,1,1)+qgfap(qgzsim,1,2))/7.5d0
      else
       qgzsim=zmin*((1.d0-zmin)/zmin)**qgran(b10)
       gb=qgzsim*qgfap(qgzsim,2,1)*.375d0
      endif
      qt=qq*(qgzsim*(1.d0-qgzsim))**2
      gb=gb/dlog(qt/alm)*qlf
      if(debug.ge.3)write (moniou,203)qt,gb
      if(qgran(b10).gt.gb)goto 1

      if(debug.ge.3)write (moniou,202)qgzsim
201   format(2x,'qgzsim - z-share simulation: qq=',e10.3,2x,'j=',i1)
202   format(2x,'qgzsim=',e10.3)
203   format(2x,'qgzsim: qt=',e10.3,2x,'gb=',e10.3)
      return
      end

c=============================================================================
      subroutine qgdifr(wppr,wmtg,ptxpr,ptypr,ptxtg,ptytg
     *,izp,izt,jexpr,jextg,iret)
c-----------------------------------------------------------------------------
c qgdifr - treatment of diffractive excitations
c wppr - LC momentum for projectile remnant;
c wptg - LC momentum for target remnant;
c izp  - projectile remnant type;
c izt  - target remnant type;
c jexpr/jextg = -1 - more collisions to follow;
c             =  0 - no excitation;
c             =  2 - low mass excitation
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ep1(4),ep2(4),ept(4),ey(3)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wplab,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr8/  pt2w,be(4),dc(5),deta,drho,almpt,ptdif
      common /qgarr10/ am(6)
      common /qgarr11/ b10
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr21/ dmmin(3),dmres(3),wdres(3)
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /npn/npp,nnn
      external qgran

      if(debug.ge.2)write (moniou,201)izp,izt,wppr,wmtg,jexpr,jextg

      iret=0
      jexip=0
      jexit=0
      if(jexpr.eq.2)jexip=1                     !low mass excitation for proj.
      if(jextg.eq.2)jexit=1                     !low mass excitation for targ.
      if(jexpr*(jexpr+1)*(jexpr-2).ne.0.or.jextg*(jextg+1)*(jextg-2)
     *.ne.0)stop'qgdifr: wrong states?!'
      if(jexpr.eq.-1.and.jextg.eq.-1)stop'qgdifr: no diffraction?!'
      
      ddmin1=0.d0
      ddmax1=0.d0
      sd0=wppr*wmtg-(ptxpr+ptxtg)**2-(ptypr+ptytg)**2 !energy squared available
      if(sd0.le.0.d0)goto 5

      jj=int(1.5d0+qgran(b10))                   !start from proj./targ. side
      if(jextg.eq.-1)then                        !more collisions to follow
       dmass2=0.d0
       ddmin2=0.d0
      elseif(jexit.eq.0)then                     !no excitation
       if(iabs(izt).eq.7.or.iabs(izt).eq.8)then  !delta++/-
        dmass2=dmmin(2)
       else
        dmass2=am(2)
       endif
       ddmin2=dmass2
      else                                       !low mass diffraction
       ddmin2=dmmin(2)+am(1)
      endif
      if(jexpr.eq.-1)then                        !more collisions to follow
       dmass1=0.d0
       ddmin1=0.d0
      elseif(jexip.eq.0)then                     !no excitation
       if(iabs(izp).eq.7.or.iabs(izp).eq.8)then  !delta++/-
        dmass1=dmmin(2)
       else
        dmass1=am(icz)
       endif
       ddmin1=dmass1
      else                                       !low mass diffraction
       ddmin1=dmmin(icz)+am(1)
      endif
      if(jj.eq.2)goto 2

1     if(jexpr.eq.2)then                         !mass for excited proj. remnant
       if(izp.eq.0)stop'LMD:izp=0?!'
       if(jj.eq.1)then
        ddmax1=dsqrt(sd0)-ddmin2
       else
        ddmax1=dsqrt(sd0)-dmass2
       endif
       if(ddmax1.le.ddmin1)goto 5

       ddmax=min(ddmax1,dmres(icz)*dsqrt(1.d0+wdres(icz)/dmres(icz)))
       if(ddmax1.gt.ddmax)then
        wres=1.d0/(1.d0+ddmax**2/dmres(icz)/wdres(icz)
     *  *(1.d0-ddmax/ddmax1)/(pi/4.d0-atan((ddmin1**2
     *  -dmres(icz)**2)/dmres(icz)/wdres(icz))))
       else
        wres=1.d0
       endif
       if(qgran(b10).gt.wres)then !PPR contribution
        dmass1=ddmax/(1.d0-qgran(b10)*(1.d0-ddmax/ddmax1))
       else                                      !resonance contribution
        dmass1=dmres(icz)*dsqrt(1.d0+wdres(icz)/dmres(icz)*tan(
     *  atan((ddmax**2-dmres(icz)**2)/dmres(icz)/wdres(icz))-qgran(b10)
     *  *(atan((ddmax**2-dmres(icz)**2)/dmres(icz)/wdres(icz))
     *  -atan((ddmin1**2-dmres(icz)**2)/dmres(icz)/wdres(icz)))))
        jexip=0
        izp=izp+10*izp/iabs(izp)
       endif
      endif
      if(jj.eq.2)goto 3

2     if(jextg.eq.2)then                         !mass for excited targ. remnant
       if(jj.eq.1)then
        ddmax2=dsqrt(sd0)-dmass1
       else
        ddmax2=dsqrt(sd0)-ddmin1
       endif
       if(ddmax2.le.ddmin2)goto 5

       ddmax=min(ddmax2,dmres(2)*dsqrt(1.d0+wdres(2)/dmres(2)))
       if(ddmax2.gt.ddmax)then
        wres=1.d0/(1.d0+ddmax**2/dmres(2)/wdres(2)*(1.d0-ddmax/ddmax2)
     *  /(pi/4.d0-atan((ddmin2**2-dmres(2)**2)/dmres(2)/wdres(2))))
       else
        wres=1.d0
       endif
       if(qgran(b10).gt.wres)then !PPR contribution
        dmass2=ddmax/(1.d0-qgran(b10)*(1.d0-ddmax/ddmax2))
       else                                      !resonance contribution
        dmass2=dmres(2)*dsqrt(1.d0+wdres(2)/dmres(2)
     *  *tan(atan((ddmax**2-dmres(2)**2)/dmres(2)/wdres(2))
     *  -qgran(b10)*(atan((ddmax**2-dmres(2)**2)/dmres(2)/wdres(2))
     *  -atan((ddmin2**2-dmres(2)**2)/dmres(2)/wdres(2)))))
        izt=izt+10*izt/iabs(izt)
        jexit=0
       endif
      endif
      if(jj.eq.2)goto 1

3     if(sd0.le.(dmass1+dmass2)**2)goto 5
      dmass1=dmass1**2
      dmass2=dmass2**2
      ptmax=max(0.d0,qglam(sd0,dmass1,dmass2))
      nret=0
      nren=0
      
4     nret=nret+1
      if(nret.gt.10)then
       if(nren.gt.10)goto 5
       ptmax=ptmax/2.d0
       nret=0
       nren=nren+1
      endif
      pti=-ptdif**2*dlog(1.d0-qgran(b10)*(1.d0-exp(-ptmax/ptdif**2)))
      pt=dsqrt(pti)
      call qgcs(c,s)
      
      if(jexpr.eq.-1)then
       ep2(3)=ptxtg+pt*c
       ep2(4)=ptytg+pt*s
       wpd2=(dmass2+ep2(3)**2+ep2(4)**2)/wmtg
       if(wpd2.gt.wppr)goto 4
       ep2(1)=.5d0*(wpd2+wmtg)
       ep2(2)=.5d0*(wpd2-wmtg)
       wmtg=0.d0
       wppr=wppr-wpd2
       ptxpr=ptxpr-pt*c
       ptypr=ptypr-pt*s
      elseif(jextg.eq.-1)then
       ep1(3)=ptxpr+pt*c
       ep1(4)=ptypr+pt*s
       wmd1=(dmass1+ep1(3)**2+ep1(4)**2)/wppr
       if(wmd1.gt.wmtg)goto 4
       ep1(1)=.5d0*(wppr+wmd1)
       ep1(2)=.5d0*(wppr-wmd1)
       wppr=0.d0
       wmtg=wmtg-wmd1
       ptxtg=ptxtg-pt*c
       ptytg=ptytg-pt*s

      else
       ept(1)=.5d0*(wppr+wmtg)
       ept(2)=.5d0*(wppr-wmtg)
       ept(3)=ptxpr+ptxtg
       ept(4)=ptypr+ptytg
       call qgdeft(sd0,ept,ey)               !boost to c.m. for the two remnants
       amt1=dmass1+pti
       amt2=dmass2+pti
       wpd1=dsqrt(sd0)*qgtwd(sd0,amt1,amt2)
       wpd2=dsqrt(sd0)-wpd1
       wmd1=amt1/wpd1
       wmd2=amt2/wpd2

       ep1(1)=.5d0*(wpd1+wmd1)
       ep1(2)=.5d0*(wpd1-wmd1)
       ep1(3)=pt*c
       ep1(4)=pt*s
       ep2(1)=.5d0*(wpd2+wmd2)
       ep2(2)=.5d0*(wpd2-wmd2)
       ep2(3)=-pt*c
       ep2(4)=-pt*s
       call qgtran(ep1,ey,1)
       call qgtran(ep2,ey,1)
       wppr=0.d0
       wmtg=0.d0
      endif

      if(jexpr.ne.-1)then
       if(jexip.eq.0)then
        call qgreg(ep1,izp)
       else
        if(izp.ne.0)is=iabs(izp)/izp     
        if(icz.eq.1)then
         if(iabs(izp).ge.4)then
          stop'qgdifr?'
          ic2=-4*is
          ic1=izp-3*is
         elseif(izp.ne.0)then
          ic1=izp*(1-3*int(.5d0+qgran(b10)))
          ic2=-izp-ic1
         else
          stop'qgdifr?'
          ic1=int(1.5d0+qgran(b10))*(2*int(.5d0+qgran(b10))-1)
          ic2=-ic1
         endif
        elseif(icz.eq.2)then
         if(iabs(izp).lt.7)then
          if(qgran(b10).gt..33333d0)then
           ic1=3*is
           ic2=izp-is
          else
           ic1=izp+4*is
           ic2=4*is-izp
          endif
         else
          ic1=izp-is
          ic2=izp-6*is
         endif
        elseif(icz.eq.3)then
         ic1=-4*is
         ic2=izp-3*is
        endif   
        call qgdeft(dmass1,ep1,ey)               !boost to the remnant rest frame
        ept(1)=.5d0*dsqrt(dmass1)
        ept(2)=ept(1)
        ept(3)=0.d0
        ept(4)=0.d0
        ep1(1)=ept(1)
        ep1(2)=-ept(1)
        ep1(3)=0.d0
        ep1(4)=0.d0
        call qgtran(ep1,ey,1)
        call qgtran(ept,ey,1)
        call qggene(ept,ep1,ic1,ic2,1)
       endif
      endif
      
      if(jextg.ne.-1)then
       if(jexit.eq.0)then
        call qgreg(ep2,izt)
       else
        is=iabs(izt)/izt
        if(iabs(izt).lt.7)then
         if(qgran(b10).gt..33333d0)then
          ic1=3*is
          ic2=izt-is
         else
          ic1=izt+4*is
          ic2=4*is-izt
         endif
        else
         stop'qgdifr?'
         ic1=izt-is
         ic2=izt-6*is
        endif
        call qgdeft(dmass2,ep2,ey)               !boost to the remnant rest frame
        ep1(1)=.5d0*dsqrt(dmass2)
        ep1(2)=ep1(1)
        ep1(3)=0.d0
        ep1(4)=0.d0
        ep2(1)=ep1(1)
        ep2(2)=-ep1(1)
        ep2(3)=0.d0
        ep2(4)=0.d0
        call qgtran(ep1,ey,1)
        call qgtran(ep2,ey,1)
        call qggene(ep1,ep2,ic2,ic1,1)
       endif
      endif
      return
5     iret=1
   
      if(debug.ge.3)write (moniou,202)
201   format(2x,'qgdifr - leading clusters hadronization:'
     */4x,'cluster types izp=',i2,2x,
     *'izt=',i2/4x,'available light cone momenta: wppr=',e10.3,
     *' wmtg=',e10.3,' jexpr',i3,' jextg',i3)
202   format(2x,'qgdifr - end')
      return
      end

c=============================================================================
      subroutine qgfau(b,gz)
c-----------------------------------------------------------------------------
c integrands for hadron-hadron and hadron-nucleus cross-sections calculation
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208,nfock=3)
      dimension gz(3),gz0(5)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)b

      do l=1,3
       gz(l)=0.d0
      enddo

      ab=float(ia(2))
      do iddp1=1,nfock
      do iddp2=1,nfock
       call qgfz(b,gz0,iddp1,iddp2)
       if(iddp1.eq.iddp2)gz(1)=gz(1)+(1.d0-gz0(1)*anorm)**ab
     * *cc(iddp1,icz)
       do l=2,3
        gz(l)=gz(l)+(1.d0-gz0(l-1)*anorm)**ab
     *  *cc(iddp1,icz)*cc(iddp2,icz)
       enddo
      enddo
      enddo
      gprod=1.d0-gz(3)
      gz(3)=gz(2)-gz(3)
      gz(2)=gz(1)-gz(2)
      gz(1)=gprod

      if(debug.ge.2)write (moniou,203)gz
      if(debug.ge.3)write (moniou,202)
201   format(2x,'qgfau - integrands for hadron-hadron and hadron'
     *,'-nucleus cross-sections calculation'/4x,'b=',e10.3)
202   format(2x,'qgfau - end')
203   format(2x,'qgfau: gz=',3e10.3)
      return
      end

c=============================================================================
      subroutine qgfrag(sa,na,rc)
c-----------------------------------------------------------------------------
c connected nucleon clasters extraction - used for the nuclear spectator part
c multifragmentation
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=208)
      dimension sa(iapmax,3)
      common /qgarr13/ nsf,iaf(iapmax)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.2)write (moniou,201)na
      if(debug.ge.3)then
       write (moniou,203)
       do i=1,na
        write (moniou,204)(sa(i,l),l=1,3)
       enddo
      endif

      ni=1
      ng=1
      j=0
1     j=j+1
      j1=ni+1

      do 4 i=j1,na
       ri=0.d0
       do m=1,3
        ri=ri+(sa(j,m)-sa(i,m))**2
       enddo
       if(ri.gt.rc)goto 4

       ni=ni+1
       ng=ng+1
       if(i.eq.ni)goto 4
       do m=1,3
        s0=sa(ni,m)
        sa(ni,m)=sa(i,m)
        sa(i,m)=s0
       enddo
4     continue

      if(j.lt.ni.and.na-ni.gt.0)goto 1
      nsf=nsf+1
      iaf(nsf)=ng
      if(debug.ge.3)write (moniou,206)nsf,iaf(nsf)

      ng=1
      j=ni
      ni=ni+1
      if(na.eq.ni)then
       nsf=nsf+1
       iaf(nsf)=1
       if(debug.ge.3)write (moniou,206)nsf,iaf(nsf)
      elseif(na.gt.ni)then
       goto 1
      endif

      if(debug.ge.3)write (moniou,202)
201   format(2x,'qgfrag-multifragmentation: nucleus mass number: na='
     *,i2)
202   format(2x,'qgfrag - end')
203   format(2x,'nucleons coordinates:')
204   format(2x,3e10.3)
206   format(2x,'qgfrag: fragment n',i2,2x,'fragment mass - ',i2)
      return
      end

c=============================================================================
      subroutine qgfrgm(ns,xa)
c-----------------------------------------------------------------------------
c qgfrgm - fragmentation of the spectator part of the nucleus
c xa     - array for spectator nucleon positions;
c ns     - total number of spectators;
c nsf    - number of secondary fragments;
c iaf(i) - mass of i-th fragment
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter(iapmax=208)
      dimension xa(iapmax,3)
      integer debug
      common /qgarr1/  ia(2),icz,icp
      common /qgarr3/  rmin,emax,eev
      common /qgarr11/ b10
      common /qgarr13/ nsf,iaf(iapmax)
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      if(debug.ge.2)write (moniou,201)ns

      nsf=0
      if(ns.eq.0)then                  !no fragments
       return
      elseif(ns.eq.1)then              !single spectator nucleon recorded
       nsf=nsf+1
       iaf(nsf)=1
       if(debug.ge.3)write (moniou,205)
       return
      endif

      eex=0.d0                         !excitation energy for spectator part
           !sum of excitations due to wounded nucleons (including diffractive)
      do i=1,ia(1)-ns
c partial excitation according to f(e) ~ 1/sqrt(e) * exp(-e/(2*<e>))
       eex=eex+(qgran(b10)+qgran(b10)+qgran(b10)+
     * qgran(b10)+qgran(b10)-2.5d0)**2*2.4d0
      enddo
      if(debug.ge.3)write (moniou,203)eex

      if(eex/ns.gt.emax)then    !if eex>emax -> multifragmentation
       call qgfrag(xa,ns,rmin)  !multifragmentation (percolation algorithm)
      else                      !otherwise eveporation
       nf=npgen(eex/eev,0,ns-1) !number of eveporated nucleons (mean=eex/eev)
       nsf=nsf+1
       iaf(nsf)=ns-nf           !recording of the fragment produced
       if(debug.ge.3)write (moniou,206)iaf(nsf)

       nal=nf/4                 !number of evapotared alphas (taken as nf/4)
       if(nal.ne.0)then
        do i=1,nal              !recording the evaporated alphas
         nsf=nsf+1
         iaf(nsf)=4
        enddo
       endif
       nf=nf-4*nal

       if(nf.ne.0)then
        do i=1,nf               !recording the evaporated nucleons
         nsf=nsf+1
         iaf(nsf)=1
        enddo
       endif
       if(debug.ge.3)write (moniou,204)nf,nal
      endif

      if(debug.ge.3)write (moniou,202)
201   format(2x,'qgfrgm: number of spectators: ns=',i2)
202   format(2x,'qgfrgm - end')
203   format(2x,'qgfrgm: excitation energy: eex=',e10.3)
204   format(2x,'qgfrgm - evaporation: number of nucleons nf='
     *,i2,'number of alphas nal=',i2)
205   format(2x,'qgfrgm - single spectator')
206   format(2x,'qgfrgm - evaporation: mass number of the fragment:',i2)
      return
      end

c=============================================================================
      subroutine qggau(gz)
c-----------------------------------------------------------------------------
c impact parameter integration for impact parameters <bm -
c for hadron-hadron and hadron-nucleus cross-sections calculation
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension gz(3),gz0(3)
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)

      do i=1,3
       gz(i)=0.d0
      enddo
      anorm=qganrm(rnuc(2),wsnuc(2),wbnuc(2))  !density normalization
      bm=rnuc(2)+2.d0*wsnuc(2)
      do i=1,7
      do m=1,2
       b=bm*dsqrt(.5d0+x1(i)*(m-1.5d0))
       call qgfau(b,gz0)
       do l=1,3
        gz(l)=gz(l)+gz0(l)*a1(i)
       enddo
      enddo
      enddo
      do l=1,3
       gz(l)=gz(l)*bm**2*pi*.5d0
      enddo

      do i=1,7
      do m=1,2
       b=bm-wsnuc(2)*dlog(.5d0+x1(i)*(m-1.5d0))
       call qgfau(b,gz0)
       do l=1,3
        gz(l)=gz(l)+gz0(l)*a1(i)*exp((b-bm)/wsnuc(2))*b*pi*wsnuc(2)
       enddo
      enddo
      enddo

      if(debug.ge.3)write (moniou,202)
201   format(2x,'qggau - nuclear cross-sections calculation')
202   format(2x,'qggau - end')
      return
      end

c=============================================================================
      double precision function qganrm(rnuc,wsnuc,wbnuc)
c-----------------------------------------------------------------------------
c impact parameter integration for impact parameters <bm -
c for hadron-hadron and hadron-nucleus cross-sections calculation
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)

      qganrm=0.d0
      do i=1,7
      do m=1,2
       r=rnuc*(.5d0+x1(i)*(m-1.5d0))**(1.d0/3.d0)
       quq=(r-rnuc)/wsnuc
       if(quq.lt.1.d80)qganrm=qganrm+a1(i)/(1.d0+exp(quq))
     * *(1.d0+wbnuc*(r/rnuc)**2)
      enddo
      enddo
      qganrm=qganrm*rnuc**3*pi/1.5d0

      dnrm=0.d0
      do i=1,7
      do m=1,2
       t=.5d0+x1(i)*(m-1.5d0)
       r=rnuc-wsnuc*log(t)
       dnrm=dnrm+a1(i)/(1.d0+t)*r*r
     * *(1.d0+wbnuc*(r/rnuc)**2)
      enddo
      enddo
      qganrm=1.d0/(qganrm+dnrm*2.d0*pi*wsnuc)

      if(debug.ge.3)write (moniou,202)qganrm
201   format(2x,'qganrm - nuclear density normalization')
202   format(2x,'qganrm=',e10.3)
      return
      end
      
c=============================================================================
      double precision function qgptgen(ptmax,bet)
c-----------------------------------------------------------------------------
c pt-generation for string fragmentation
c ptmax - kinematic limit,
c bet   - parameter (related to average pt)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran
      
      if(ptmax.lt.bet)then
1      pti=ptmax*dsqrt(qgran(b10))
       if(qgran(b10).gt.exp(-pti/bet))goto 1
      else
2      pti=-bet*dlog(qgran(b10)*qgran(b10))
       if(pti.gt.ptmax)goto 2
      endif
      qgptgen=pti
      return
      end
      
c=============================================================================
      subroutine qggene(epp0,epm0,icp0,icm0,jregg)
c-----------------------------------------------------------------------------
c qggene - string fragmentation into secondary hadrons
c wp0, wm0 - initial light cone momenta ( e+p, e-p ) of string ends;
c ic1, ic2 - their types (1 - u, -1 - U, 2 - d, -2 - D, 3 - ud, -3 - UD, 4 - s,
c -4 - S, 6 - uu, -6 - UU, 7 - dd, -7 - DD, 8 - us, -8 - US, 9 - ds, -9 - DS)
c jregg=1 for Reggeon string (at least 3 particles to produce)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      character *2 tyq
      dimension wp(2),ey(3),ept(4),ep(4),ep1(4),ep2(4),amj(2)
     *,epp0(4),epm0(4),ic(2),jend(2)
      common /qgarr8/  pt2w,bep,ben,bek,bec,dc(5),deta,drho,almpt,ptdif
      common /qgarr10/ am0,amn,amk,amsig,amlam,ameta
      common /qgarr11/ b10
      common /qgarr21/ dmmin(3),dmres(3),wdres(3)
      common /qgarr28/ arr(5),alpq
      common /qgarr42/ tyq(16)
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      if(debug.ge.2)write (moniou,201)tyq(8+icp0),tyq(8+icm0)
     *,epp0,epm0
     
      do i=1,4
       ep1(i)=epp0(i)
       ep2(i)=epm0(i)
      enddo
      
      ic(1)=icp0                        !parton types at string ends
      ic(2)=icm0
      j3part=jregg
      jend(1)=1
      jend(2)=1

      if(iabs(ic(1)).eq.5.or.iabs(ic(2)).eq.5.or.iabs(ic(1)).ge.8
     *.or.iabs(ic(2)).ge.8)stop'qggene: problem with parton types'

1     do i=1,4
       ept(i)=ep1(i)+ep2(i)
      enddo
      ww=qgnrm(ept)                     !mass squared for the string
      sww=dsqrt(ww)

      if(ic(1)*ic(2).gt.0.and.(iabs(ic(1)).le.2.or.iabs(ic(1)).eq.4)
     *.and.(iabs(ic(2)).le.2.or.iabs(ic(2)).eq.4).or.ic(1)*ic(2).lt.0
     *.and.((iabs(ic(1)).le.2.or.iabs(ic(1)).eq.4).and.(iabs(ic(2)).eq.3
     *.or.iabs(ic(2)).ge.6).or.(iabs(ic(1)).eq.3.or.iabs(ic(1)).ge.6)
     *.and.(iabs(ic(2)).le.2.or.iabs(ic(2)).eq.4)))
     *stop'qggene: problem with parton types'
      
      do i=1,2
       iab=iabs(ic(i))
       if(iab.le.2)then
        amj(i)=am0
       elseif(iab.eq.4)then
        amj(i)=amk
       elseif(iab.eq.6.or.iab.eq.7)then
        amj(i)=dmmin(2)
       elseif(iab.ge.8)then
        amj(i)=amsig
       else
        amj(i)=amn
       endif
      enddo

      j=int(1.5d0+qgran(b10))           !choose string end to start
      iab=iabs(ic(j))
      iaj=iabs(ic(3-j))
      is=ic(j)/iab
      
      if(debug.ge.3)then
       iqt=8+ic(j)
       iqt2=8+ic(3-j)
       if(iqt.gt.0.and.iqt2.gt.0)
     * write (moniou,203)j,tyq(iqt),tyq(iqt2),ww
      endif

      call qgdeft(ww,ept,ey)            !boost to c.m. for the string
      call qgtran(ep1,ey,-1)
      call qgdefr(ep1,s0x,c0x,s0,c0)
      
c      if(j3part.ne.0.or.sww.gt.amj(1)+amj(2)+am0+wwm)then !>2 particles
      if(j3part.ne.0.or.sww.gt.dsqrt(amj(1)**2+pt2w)
     *+dsqrt(amj(2)**2+pt2w)+dsqrt(am0**2+pt2w))then        !>2 particles
       blf=0.d0
       bet=0.d0
       alf=0.d0
       if(iab.le.2)then                 !light (anti-)quark at the string end
        wqs=0.d0                        !us/ds,US,DS-diquark
        wud=0.d0                        !ud/UD-diquark
        wss=0.d0                        !s/S-quark
        if(sww.gt.dsqrt(amj(3-j)**2+pt2w)+2.d0*dsqrt(amlam**2+pt2w))
     *  wqs=dc(4)
        if(sww.gt.dsqrt(amj(3-j)**2+pt2w)+2.d0*dsqrt(amn**2+pt2w))
     *  wud=dc(1)
        if(sww.gt.dsqrt(amj(3-j)**2+pt2w)+2.d0*dsqrt(amk**2+pt2w))
     *  wss=dc(2)
        aks=qgran(b10)*(1.d0+wqs+wud+wss)
        
        if(aks.lt.wqs)then                                  !lambda/sigma
         alf=almpt-arr(2)+arr(1)-arr(3)
         blf=1.d0-arr(2)-arr(3)
         bet=ben
         s2min=(amj(3-j)+amlam)**2
         if(qgran(b10).lt..5d0.and.sww.gt.amj(3-j)+amlam+amsig)then
          ami=amsig**2
          ic0=ic(j)+19*is               !sigma+-
          icc=-7*is-ic(j)               !US(DS/us/ds)
         else
          ami=amlam**2
          ic0=6*is                      !(anti-)lambda
          icc=-10*is+ic(j)              !US(DS/us/ds)
         endif
        elseif(aks.lt.wqs+wud)then                          !nucleon
         alf=almpt-arr(2)
         blf=1.d0-arr(1)-arr(2)
         bet=ben
         ami=amn**2
         s2min=(amj(3-j)+amn)**2
         ic0=ic(j)+is                   !p,n(P,N)
         icc=-3*is                      !UD(ud)
        elseif(aks.lt.wqs+wud+wss)then                      !kaon
         alf=almpt-arr(3)
         blf=1.d0-arr(1)-arr(3)
         bet=bek
         ami=amk**2
         s2min=(amj(3-j)+amk)**2
         ic0=ic(j)+3*is                 !K+(K-,K0,K0~)
         icc=4*is                       !s(S)
        else
         weta=0.d0
         wrho=0.d0
         if(sww.gt.dsqrt(amj(3-j)**2+pt2w)+dsqrt(ameta**2+pt2w)
     *   +dsqrt(am0**2+pt2w))weta=deta
         if(sww.gt.dsqrt(amj(3-j)**2+pt2w)+dsqrt(dmmin(1)**2+pt2w)
     *   +dsqrt(am0**2+pt2w))wrho=drho
         aks=qgran(b10)*(1.d0+weta+wrho)
         alf=almpt-arr(1)
         blf=1.d0-2.d0*arr(1)
         bet=bep
         s2min=(amj(3-j)+am0)**2
         if(aks.lt.weta)then                                !eta
          ami=ameta**2
          ic0=10                        !eta
          icc=ic(j)                     !unchanged
         elseif(aks.lt.weta+wrho)then                       !rho
          ami=dmmin(1)**2
          if(qgran(b10).lt..3333d0)then
           ic0=17                       !rho0
           icc=ic(j)                    !unchanged
          else
           ic0=16*(3*is-2*ic(j))        !rho+-
           icc=3*is-ic(j)               !d,u(D,U)
          endif
         else                                               !pion
          ami=am0**2
          if(qgran(b10).lt..3333d0)then
           ic0=0                        !pi0
           icc=ic(j)                    !unchanged
          else
           ic0=3*is-2*ic(j)             !pi+-
           icc=3*is-ic(j)               !d,u(D,U)
          endif
         endif
        endif

       elseif(iab.ge.8)then             !us(ds/US/DS)
        alf=almpt-arr(1)
        blf=1.d0-arr(2)-arr(3)
        bet=bep
        s2min=(amj(3-j)+am0)**2
        if(qgran(b10).lt..5d0)then
         ami=amsig**2
         ic0=ic(j)+12*is                !sigma+-
         icc=7*is-ic(j)                 !U(D/u/d)
        else
         ami=amlam**2
         ic0=6*is                       !(anti-)lambda
         icc=ic(j)-10*is                !D(U/d/u)
        endif
       
       elseif(iab.eq.4)then             !strange (anti-)quark as the string end
        wud=0.d0                        !ud/UD-diquark
        if(sww.gt.dsqrt(amj(3-j)**2+pt2w)+dsqrt(amn**2+pt2w)
     *  +dsqrt(amlam**2+pt2w))wud=dc(1)
        aks=qgran(b10)*(1.d0+wud)
        if(aks.lt.wud)then                                  !lambda
         alf=almpt-arr(2)
         blf=1.d0-arr(2)-arr(3)
         bet=ben
         ami=amlam**2
         s2min=(amj(3-j)+amn)**2
         ic0=6*is                       !(anti-)lambda
         icc=-3*is                      !UD(ud)
        else                                                !kaon
         alf=almpt-arr(1)
         blf=1.d0-arr(1)-arr(3)
         bet=bep
         ami=amk**2
         s2min=(amj(3-j)+am0)**2
         icc=is*int(1.5d0+qgran(b10))   !u,d(U,D)
         ic0=-3*is-icc                  !K-(K+,K0,K0~)
        endif

       else                             !ud(uu,dd,UD,UU,DD) diquark string end
        wss=0.d0
        if(iab.eq.3.and.sww.gt.dsqrt(amj(3-j)**2+pt2w)
     *  +dsqrt(amk**2+pt2w)+dsqrt(amlam**2+pt2w).or.iab.ge.6
     *  .and.sww.gt.dsqrt(amj(3-j)**2+pt2w)+dsqrt(amk**2+pt2w)
     *  +dsqrt(amsig**2+pt2w))wss=dc(2)*(1-jend(j))+dc(3)*jend(j) !'intr. str.'
        aks=qgran(b10)*(1.d0+wss)
        if(aks.lt.wss)then                                  !lambda
         alf=almpt-arr(3)
         blf=1.d0-arr(2)-arr(3)
         if(iab.ne.3)blf=blf+1.d0       !uu-diquark!
         bet=bek
         s2min=(amj(3-j)+amk)**2
         if(iab.eq.3)then
          ami=amlam**2
          ic0=6*is                      !(anti-)lambda
         else
          ami=amsig**2
          ic0=ic(j)+14*is               !sigma+-
         endif
         icc=-4*is                      !S(s)
        else                                                !nucleon
         alf=almpt-arr(1)
         blf=1.d0-arr(1)-arr(2)
         if(iab.ne.3)blf=blf+1.d0       !uu-diquark!
         bet=ben
         ami=amn**2
         s2min=(amj(3-j)+am0)**2
         if(iab.eq.3)then
          ic0=is*int(2.5d0+qgran(b10))  !p,n(P,N)
          icc=is-ic0                    !U,D(u,d)
         elseif(qgran(b10).lt..5d0.or.sww.le.amj(3-j)+am0+dmmin(2))then
          ic0=ic(j)-4*is                !p,n(P,N)
          icc=ic0-4*is                  !D,U(d,u)
         else
          ami=dmmin(2)**2
          ic0=ic(j)+is                                      !delta++/-
          icc=6*is-ic0
         endif
        endif
       endif

       ptmax=dsqrt(max(0.d0,qglam(ww,s2min,ami)))
       pti=qgptgen(ptmax,bet)
       amt=ami+pti**2
       restm=s2min+pti**2
       zmax=qgtwd(ww,amt,restm)
       zmin=1.d0-qgtwd(ww,restm,amt)
       z1=(1.d0-zmax)**alf
       z2=(1.d0-zmin)**alf
2      z=1.-(z1+(z2-z1)*qgran(b10))**(1./alf)
       gba=exp(2.d0*amt*(1.d0-z)/z*dlog(1.d0-z)
     * -2.d0*amt*(1.d0-zmax)/zmax*dlog(1.d0-zmax))
       if(gba.gt.1.d0)write(*,*)'gba',gba,z/zmax
       if(qgran(b10).gt.(z/zmax)**blf*gba)goto 2
       
c       if(qgran(b10).gt.(z/zmax)**blf)goto 2

       wp(j)=z*sww
       wp(3-j)=amt/wp(j)
       ep(1)=.5d0*(wp(1)+wp(2))
       ep(2)=.5d0*(wp(1)-wp(2))
       call qgcs(c,s)
       ep(3)=pti*c
       ep(4)=pti*s
       if(s0x.ne.0.d0.or.s0.ne.0.d0)call qgrota(ep,s0x,c0x,s0,c0)
       call qgtran(ep,ey,1)
       call qgreg(ep,ic0)
      
       wp(j)=sww-wp(j)
       wpm=pti*pti/wp(j)
       wp(3-j)=sww-wp(3-j)-wpm
       if(j.eq.1)then
        ep1(1)=.5d0*(wp(1)+wpm)
        ep1(2)=.5d0*(wp(1)-wpm)
        ep1(3)=-pti*c
        ep1(4)=-pti*s
        ep2(1)=.5d0*wp(2)
        ep2(2)=-ep2(1)
        ep2(3)=0.d0
        ep2(4)=0.d0
       else
        ep2(1)=.5d0*(wpm+wp(2))
        ep2(2)=.5d0*(wpm-wp(2))
        ep2(3)=-pti*c
        ep2(4)=-pti*s
        ep1(1)=.5d0*wp(1)
        ep1(2)=ep1(1)
        ep1(3)=0.d0
        ep1(4)=0.d0
       endif
       if(s0x.ne.0.d0.or.s0.ne.0.d0)then
        call qgrota(ep1,s0x,c0x,s0,c0)
        call qgrota(ep2,s0x,c0x,s0,c0)
       endif
       call qgtran(ep1,ey,1)
       call qgtran(ep2,ey,1)
       ic(j)=icc
       j3part=0
       jend(j)=0
       goto 1
       
      else
       bet=bep
       if(iab.ge.8.or.iaj.ge.8)then     !us,ds,US,DS
        if(iab.lt.8)then
         j=3-j
         iab=iabs(ic(j))
         iaj=iabs(ic(3-j))
         is=ic(j)/iabs(ic(j))
        endif
        if(iaj.ge.8)then
         if(ic(j)+ic(3-j).eq.0)then
          if(qgran(b10).lt..5d0.or.sww.le.2.d0*amsig)then
           ic(j)=6*is
          else
           ic(j)=ic(j)+12*is
          endif
          ic(3-j)=-ic(j)
         else
          if(qgran(b10).lt..5d0)then
           ic(j)=6*is
           ic(3-j)=ic(3-j)-12*is
          else
           ic(j)=ic(j)+12*is
           ic(3-j)=-6*is
          endif
         endif
        elseif(iaj.ge.6)then
         if(iaj.eq.6)then
          ic(3-j)=-2*is
          if(iabs(ic(j)).eq.8)then
           ic(j)=6*is
          else
           ic(j)=21*is
          endif
         else
          ic(3-j)=-3*is
          if(iabs(ic(j)).eq.8)then
           ic(j)=20*is
          else
           ic(j)=6*is
          endif
         endif
        elseif(iaj.eq.4)then
         if(qgran(b10).lt..5d0.or.sww.le.amsig+amk)then
          ic(3-j)=ic(j)-13*is
          ic(j)=6*is
         else
          ic(3-j)=4*is-ic(j)
          ic(j)=ic(j)+12*is
         endif
        elseif(iaj.eq.3)then
         if(qgran(b10).lt..5d0.or.sww.le.amsig+amn)then
          ic(3-j)=ic(j)-11*is
          ic(j)=6*is
         else
          ic(3-j)=-ic(j)+6*is
          ic(j)=ic(j)+12*is
         endif
        else
         if(sww.le.amsig+am0.or.(ic(j)-ic(3-j)-7*is.eq.0
     *   .and.qgran(b10).gt..3333d0.or.ic(j)-ic(3-j)-7*is.ne.0
     *   .and.qgran(b10).lt..3333d0))then
          if(ic(j)-ic(3-j)-7*is.eq.0)then
           ic(3-j)=3*is-2*ic(3-j)
          else
           ic(3-j)=0
          endif
          ic(j)=6*is
         else
          ic(3-j)=ic(j)-7*is-ic(3-j)
          ic(j)=ic(j)+12*is
         endif
        endif
       
       elseif(iab.gt.5.or.iaj.gt.5)then !uu, dd, UU, DD
        if(iab.lt.6)then
         j=3-j
         iab=iabs(ic(j))
         iaj=iabs(ic(3-j))
         is=ic(j)/iabs(ic(j))
        endif
        if(iaj.lt.6)then
         if(iaj.eq.3)then
          if(qgran(b10).gt..5d0.or.sww.le.amn+dmmin(2))then
           ic(j)=ic(j)-4*is
           ic(3-j)=ic(j)-5*is
          else 
           ic(j)=ic(j)+is
           ic(3-j)=5*is-ic(j)
          endif
         elseif(iaj.eq.4)then
          if(qgran(b10).gt..5d0.or.sww.le.amk+dmmin(2))then
           ic(j)=ic(j)-4*is
           ic(3-j)=ic(j)-7*is
          else 
           ic(j)=ic(j)+is
           ic(3-j)=3*is-ic(j)
          endif
         else
          if(ic(j).eq.5*is+ic(3-j))then
           if(qgran(b10).gt..3333d0.or.sww.le.am0+dmmin(2))then
            ic(j)=ic(j)-4*is
            ic(3-j)=3*is-2*ic(3-j)
           else 
            ic(j)=ic(j)+is
            ic(3-j)=0
           endif
          else
           if(qgran(b10).lt..3333d0.or.sww.le.am0+dmmin(2))then
            ic(j)=ic(j)-4*is
            ic(3-j)=0
           else 
            ic(j)=ic(j)+is
            ic(3-j)=3*is-2*ic(3-j)
           endif
          endif
         endif
        elseif(ic(1)+ic(2).eq.0)then
         ic(j)=ic(j)-4*is
         ic(3-j)=-ic(j)
        elseif(qgran(b10).lt..5d0)then
         ic(j)=ic(j)-4*is
         ic(3-j)=ic(3-j)-is
        else
         ic(j)=ic(j)+is
         ic(3-j)=ic(3-j)+4*is
        endif

       elseif(iab.le.2.and.iaj.le.2)then
        wud=0.d0                        !ud/UD-diquark
        wss=0.d0                        !s/S-quark
        if(sww.gt.2.d0*amn)wud=dc(1)
        if(sww.gt.2.d0*amk)wss=dc(2)
        aks=qgran(b10)*(1.d0+wud+wss)
        if(aks.lt.wud)then
         bet=ben
         ic(j)=ic(j)+is
         ic(3-j)=ic(3-j)-is
        elseif(aks.lt.wud+wss)then
         bet=bek
         ic(j)=ic(j)+3*is
         ic(3-j)=ic(3-j)-3*is
        else
         weta=0.d0
         wrho=0.d0
         ic0=-ic(1)-ic(2)
         if(sww.gt.ameta+am0)weta=deta/(3-2*iabs(ic0))
         if(sww.gt.dmmin(1)+am0)wrho=drho
         aks=qgran(b10)*(1.d0+weta+wrho)
         if(aks.lt.weta)then
          ic(j)=10              !eta
          ic(3-j)=ic0
         elseif(aks.lt.weta+wrho)then
          if(ic0.ne.0)then               !different flavors
           if(qgran(b10).lt..5d0)then
            ic(j)=16*ic0                 !rho+-
            ic(3-j)=0
           else
            ic(j)=17                     !rho0
            ic(3-j)=ic0
           endif
          else                           !same flavor
           if(qgran(b10).lt..2d0)then
            ic(j)=17                     !rho0
            ic(3-j)=0
           else
            ic(j)=3*is-2*ic(j)
            ic(3-j)=-ic(j)
            ic(j)=ic(j)*16
           endif
          endif
         else
          if(ic0.ne.0)then                !different flavors
           ic(j)=ic0
           ic(3-j)=0
          else                            !same flavor
           if(qgran(b10).lt..2d0)then
            ic(j)=0
            ic(3-j)=0
           else
            ic(j)=3*is-2*ic(j)
            ic(3-j)=-ic(j)
           endif
          endif
         endif
        endif

       elseif(iab.eq.4.or.iaj.eq.4)then  !s(S)
        if(iab.ne.4)then
         j=3-j
         iab=iabs(ic(j))
         iaj=iabs(ic(3-j))
         is=ic(j)/iabs(ic(j))
        endif
        if(iaj.gt.2)then
         ic(j)=-is*int(4.5d0+qgran(b10))
         if(iaj.eq.4)then
          ic(3-j)=-ic(j)
         else
          ic(3-j)=-ic(j)-2*is
         endif
        elseif(sww.gt.amlam+amn.and.qgran(b10)
     *  .lt.dc(1)/(1.d0+dc(1)))then
         bet=ben
         ic(j)=6*is
         ic(3-j)=ic(3-j)-is
        else
         if(qgran(b10).lt..3333d0)then
          ic(j)=ic(3-j)-3*is
          ic(3-j)=0
         else
          ic(j)=-6*is-ic(3-j)
          ic(3-j)=-2*ic(3-j)-3*is
         endif
        endif
       
       else                   !ud(UD)
        if(iab.ne.3)then
         j=3-j
         iab=iabs(ic(j))
         iaj=iabs(ic(3-j))
         is=ic(j)/iabs(ic(j))
        endif
        wss=dc(2)*(1-jend(j))+dc(3)*jend(j)
        if(iab.eq.iaj)then
         if(sww.gt.2.d0*amlam.and.qgran(b10).lt.wss/(1.d0+wss))then
          bet=bek
          ic(j)=6*is
          ic(3-j)=-ic(j)
         else
          ic(j)=is*int(2.5d0+qgran(b10))
          ic(3-j)=-ic(j)
         endif
        elseif(sww.gt.amlam+amk.and.qgran(b10).lt.wss/(1.d0+wss))then
         bet=bek
         ic(j)=6*is
         ic(3-j)=ic(3-j)+3*is
        else
         if(qgran(b10).lt..3333d0)then
          ic(j)=ic(3-j)+is
          ic(3-j)=0
         else
          ic(j)=4*is-ic(3-j)
          ic(3-j)=3*is-2*ic(3-j)
         endif
        endif
       endif

       do i=1,2
        iab=iabs(ic(i))
        if(iab.le.1)then                !pion
         amj(i)=am0**2
        elseif(iab.le.3)then            !nucleon
         amj(i)=amn**2
        elseif(iab.le.5)then            !kaon
         amj(i)=amk**2
        elseif(iab.eq.6)then            !lambda
         amj(i)=amlam**2
        elseif(iab.le.8)then            !delta
         amj(i)=dmmin(2)**2
        elseif(iab.eq.10)then           !eta
         amj(i)=ameta**2
        elseif(iab.eq.17.or.iab.eq.16)then !rho
         amj(i)=dmmin(1)**2
        elseif(iab.eq.20.or.iab.eq.21)then !sigma
         amj(i)=amsig**2
        else
          stop'qggene: problem with id'
        endif
       enddo

       ptmax=dsqrt(qglam(ww,amj(1),amj(2)))
       if(ptmax.lt.0.d0)stop'qggene: ptmax=0?!'
       pti=qgptgen(ptmax,bet)
       amt1=amj(1)+pti**2
       amt2=amj(2)+pti**2
       z=qgtwd(ww,amt1,amt2)
       wp(1)=z*sww
       wp(2)=amt1/wp(1)
       ep1(1)=.5d0*(wp(1)+wp(2))
       ep1(2)=.5d0*(wp(1)-wp(2))
       call qgcs(c,s)
       ep1(3)=pti*c
       ep1(4)=pti*s
       ep2(1)=sww-ep1(1)
       do i=2,4
        ep2(i)=-ep1(i)
       enddo
      
       if(s0x.ne.0.d0.or.s0.ne.0.d0)then
        call qgrota(ep1,s0x,c0x,s0,c0)
        call qgrota(ep2,s0x,c0x,s0,c0)
       endif
       call qgtran(ep1,ey,1)
       call qgtran(ep2,ey,1)
       call qgreg(ep1,ic(1))
       call qgreg(ep2,ic(2))
      endif
      if(debug.ge.3)write (moniou,202)
      return

201   format(2x,'qggene: parton flavors at the ends of the string:'
     *,2x,a2,2x,a2/4x,'light cone momenta of the string: ',e10.3
     *,2x,e10.3/4x,'ey0=',3e10.3/4x,'s0x=',e10.3,2x,'c0x=',e10.3
     *,2x,'s0=',e10.3,2x,'c0=',e10.3)
202   format(2x,'qggene - end')
203   format(2x,'qggene: current parton flavor at the end '
     *,i1,' of the string: ',a2/4x,' string mass: ',e10.3)
      end

c=============================================================================
      subroutine qgxjet(iret)
c-----------------------------------------------------------------------------
c procedure for jet hadronization
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(njmax=50000)
      dimension ep(4),ept(4),ept1(4),ey(3),ep2(4)
     *,epj(4,2,2*njmax),ipj(2,2*njmax)
      common /qgarr8/  pt2w,be(4),dc(5),deta,drho,almpt,ptdif
      common /qgarr10/ am(6)
      common /qgarr11/ b10
      common /qgarr21/ dmmin(3),dmres(3),wdres(3)
      common /qgarr36/ epjet(4,njmax),ipjet(njmax),njtot
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      if(debug.ge.2)write (moniou,201)njtot
201   format(2x,'qgxjet: total number of jets njtot=',i4)

      nj0=1
      njet0=0
      nrej=0
      iret=0
      nrejmax=1000

1     njet=njet0
      do i=1,4
       ept(i)=epjet(i,nj0)
       epj(i,1,njet+1)=ept(i)
      enddo
      iq1=ipjet(nj0)
      ipj(1,njet+1)=iq1
      
      if(iq1.eq.0)stop'qgxjet: iq1=0???'
      if(iabs(iq1).le.2)then
       am1=am(1)
       if(iq1.gt.0)then
        jq=1
       else
        jq=2
       endif
      elseif(iabs(iq1).eq.4)then
       am1=am(3)
       if(iq1.gt.0)then
        jq=1
       else
        jq=2
       endif
      else
       am1=am(2)
       if(iq1.gt.0)then
        jq=2
       else
        jq=1
       endif
      endif

      ij=nj0
2     ij=ij+1
      njet=njet+1
      iq2=ipjet(ij)

      if(iq2.eq.0)then
6      aks=qgran(b10)
       if(qgran(b10).gt.(4.d0*aks*(1.d0-aks)))goto 6
       do i=1,4
        epi=epjet(i,ij)*aks
        epj(i,2,njet)=epi
        ept(i)=ept(i)+epi
       enddo
       if(qgran(b10).lt.dc(2))then
        ipj(2,njet)=4*(2*jq-3)
        amj=am(3)
       else
        ipj(2,njet)=int(1.5d0+qgran(b10))*(2*jq-3)
        amj=am(1)
       endif

       if(qgnrm(ept).gt.(am1+amj)**2)then
        if(debug.ge.3)write (moniou,211)njet,ipj(1,njet),ipj(2,njet)
     *  ,qgnrm(ept),ept
        ipj(1,njet+1)=-ipj(2,njet)
        do i=1,4
         ept(i)=epjet(i,ij)-epj(i,2,njet)
         epj(i,1,njet+1)=ept(i)
        enddo
        am1=amj
        goto 2
       elseif(nrej.lt.nrejmax)then
        nrej=nrej+1
        goto 1
       else
3       do i=1,4
         ept(i)=epjet(i,ij)+epjet(i,ij-1)+epjet(i,ij+1)
         ep(i)=epjet(i,ij-1)
         ept1(i)=ept(i)
        enddo
        ww=qgnrm(ept1)
        if(ww.le.0.d0)then
         if(ij.gt.nj0+1)then
          ij=ij-1
          goto 3
         elseif(ipjet(ij+1).eq.0)then
          ij=ij+1
          goto 3
         else
          iret=1
          return
         endif
        endif
        ipjet(ij)=ipjet(ij+1)
        sww=sqrt(ww)
        call qgdeft(ww,ept1,ey)
        call qgtran(ep,ey,-1)
        call qgdefr(ep,s0x,c0x,s0,c0)
        ep(1)=.5d0*sww
        ep(2)=.5d0*sww
        ep(3)=0.d0
        ep(4)=0.d0
        call qgrota(ep,s0x,c0x,s0,c0)
        call qgtran(ep,ey,1)
        do i=1,4
         epjet(i,ij-1)=ep(i)
         epjet(i,ij)=ept(i)-ep(i)
        enddo

        if(njtot.gt.ij+1)then
         do j=ij+1,njtot-1
          ipjet(j)=ipjet(j+1)
         do i=1,4
          epjet(i,j)=epjet(i,j+1)
         enddo
         enddo
        endif
        nrej=0
        njtot=njtot-1
        goto 1
       endif

      else
       ipj(2,njet)=iq2
       do i=1,4
        epi=epjet(i,ij)
        epj(i,2,njet)=epi
        ept(i)=ept(i)+epi
       enddo
       if(iabs(iq2).le.2)then
        am2=am(1)
       elseif(iabs(iq2).eq.4)then
        am2=am(3)
       else
        am2=am(2)
       endif

       if(qgnrm(ept).gt.(am1+am2)**2)then
        if(debug.ge.3)write (moniou,211)njet,ipj(1,njet),ipj(2,njet)
     *  ,qgnrm(ept),ept

        nj0=ij+1
        njet0=njet
        nrej=0
        if(ij.lt.njtot)then
         goto 1
        else
         goto 5
        endif
       elseif(nrej.lt.nrejmax.and.ij.gt.nj0+1)then
        nrej=nrej+1
        goto 1
       else
4       if(ij.le.nj0+1)then
         iret=1
         return
        endif
        do i=1,4
         ept(i)=epjet(i,ij)+epjet(i,ij-1)+epjet(i,ij-2)
         ep(i)=epjet(i,ij-2)
         ept1(i)=ept(i)
        enddo
        ww=qgnrm(ept1)
        if(ww.le.0.d0)then
         ij=ij-1
         goto 4
        endif
        ipjet(ij-1)=ipjet(ij)
        sww=sqrt(ww)
        call qgdeft(ww,ept1,ey)
        call qgtran(ep,ey,-1)
        call qgdefr(ep,s0x,c0x,s0,c0)
        ep(1)=.5d0*sww
        ep(2)=.5d0*sww
        ep(3)=0.d0
        ep(4)=0.d0
        call qgrota(ep,s0x,c0x,s0,c0)
        call qgtran(ep,ey,1)
        do i=1,4
         epjet(i,ij-2)=ep(i)
         epjet(i,ij-1)=ept(i)-ep(i)
        enddo

        if(ij.lt.njtot)then
         do j=ij,njtot-1
          ipjet(j)=ipjet(j+1)
         do i=1,4
          epjet(i,j)=epjet(i,j+1)
         enddo
         enddo
        endif
        nrej=0
        njtot=njtot-1
        goto 1
       endif
      endif

5     do ij=1,njet
       do i=1,4
        ep(i)=epj(i,1,ij)
        ep2(i)=epj(i,2,ij)
        ept(i)=ep(i)+ep2(i)
       enddo
       ww=qgnrm(ept)  !invariant mass squared for the jet
       if(debug.ge.3)write (moniou,208)
     * ij,njet,ww,ipj(1,ij),ipj(2,ij)
       if(iabs(ipj(1,ij)).gt.5.and.iabs(ipj(2,ij)).gt.5.and.ipj(1,ij)
     * +ipj(2,ij).ne.0.and.ww.le.1.1d0*(am(2)+dmmin(2))**2)then
        iret=1
        return
       endif
       call qggene(ep,ep2,ipj(1,ij),ipj(2,ij),0)
      enddo
      
      if(debug.ge.3)write (moniou,202)      
202   format(2x,'qgxjet - end')
208   format(2x,'qgxjet: ij=',i2,2x,'njet=',i3,2x,'ww=',e10.3
     *,2x,'ic=',2i3)
211   format(2x,'qgxjet: njet=',i3,2x,'ic=',2i2,2x,'mass=',e10.3
     *,2x,'ep=',4e10.3)
      return
      end

c=============================================================================
      double precision function qgrot(b,s)
c-----------------------------------------------------------------------------
c qgrot - convolution of nuclear profile functions (axial angle integration)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)b,s

      qgrot=0.d0
      do i=1,7
      do m=1,2
       phi=pi*(.5d0+x1(i)*(m-1.5d0))
       bb=b**2+s**2+2.d0*b*s*cos(phi)
       qgrot=qgrot+a1(i)*qgt(bb)
      enddo
      enddo
      qgrot=qgrot*pi/2.d0

      if(debug.ge.2)write (moniou,202)qgrot
201   format(2x,'qgrot - axial angle integration of the ',
     *'nuclear profile function'/4x,
     *'impact parameter b=',e10.3,2x,'nucleon coordinate s=',e10.3)
202   format(2x,'qgrot=',e10.3)
      return
      end

c===========================================================================
      double precision function qgt(b)
c---------------------------------------------------------------------------
c qgt - nuclear profile function value at impact parameter squared b
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)b

      qgt=0.
      zm=rnuc(2)**2-b
      if(zm.gt.4.*b)then
       zm=dsqrt(zm)
      else
       zm=2.*dsqrt(b)
      endif

      do i=1,7
      do m=1,2
       z1=zm*(.5d0+x1(i)*(m-1.5d0))
       r=dsqrt(b+z1**2)
       quq=(r-rnuc(2))/wsnuc(2)
       if (quq.lt.85.)qgt=qgt+a1(i)/(1.+exp(quq))
     * *(1.d0+wbnuc(2)*(r/rnuc(2))**2)
      enddo
      enddo
      qgt=qgt*zm

      dt=0.
      do i=1,7
      do m=1,2
       z1=zm-wsnuc(2)*log(.5d0+x1(i)*(m-1.5d0))
       r=dsqrt(b+z1**2)
       quq=(r-rnuc(2)-z1+zm)/wsnuc(2)
       if (quq.lt.85.)dt=dt+a1(i)/(exp((zm-z1)/wsnuc(2))+exp(quq))
     * *(1.d0+wbnuc(2)*(r/rnuc(2))**2)
      enddo
      enddo
      qgt=qgt+dt*wsnuc(2)

      if(debug.ge.3)write (moniou,202)qgt
201   format(2x,'qgt - nuclear profile function value at impact'
     *,' parameter squared b=',e10.3)
202   format(2x,'qgt=',e10.3)
      return
      end

c=============================================================================
      block data qgdata
c-----------------------------------------------------------------------------
c constants for numerical integration (gaussian weights)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /arr1/ trnuc(56),twsnuc(56),twbnuc(56)
      common /arr3/ x1(7),a1(7)
      common /arr4/ x4(2),a4(2)
      common /arr5/ x5(2),a5(2)
      common /arr6/ x3(8),a3(8)
      common /arr8/ x2(4),a2
      common /arr9/ x9(3),a9(3)
      data x1/.9862838d0,.9284349d0,.8272013d0,.6872929d0,.5152486d0,
     *.3191124d0,.1080549d0/
      data a1/.03511946d0,.08015809d0,.1215186d0,.1572032d0,
     *.1855384d0,.2051985d0,.2152639d0/
      data x2/.00960736d0,.0842652d0,.222215d0,.402455d0/
      data a2/.392699d0/
      data x3/.0950125098d0,.2816035507d0,.4580167776d0,.6178762444d0
     +,.7554044083d0,.8656312023d0,.9445750230d0,.9894009349d0/
      data a3/.1894506104d0,.1826034150d0,.1691565193d0,.1495959888d0
     +,.1246289712d0,.0951585116d0,.0622535239d0, .0271524594d0/
      data x4/ 0.339981,0.861136/
      data a4/ 0.652145,0.347855/
      data x5/.585786d0,3.41421d0/
      data a5/.853553d0,.146447d0/
      data x9/.93247d0,.661209d0,.238619d0/
      data a9/.171324d0,.360762d0,.467914d0/
      data trnuc/0.69d0,1.71d0,1.53d0,1.37d0,1.37d0,2.09d0,1.95d0
     *,1.95d0,2.06d0,1.76d0,1.67d0,1.74d0,1.66d0,2.57d0,2.334d0
     *,2.608d0,2.201d0,2.331d0,2.58d0,2.791d0,2.791d0,2.782d0,2.74d0
     *,3.192d0,3.22d0,3.05d0,3.07d0,3.34d0,3.338d0,3.252d0
     *,3.369d0,3.244d0,3.244d0,3.313d0,3.476d0,3.54d0,3.554d0
     *,3.554d0,3.743d0,3.73d0,3.744d0,3.759d0,3.774d0,3.788d0
     *,3.802d0,3.815d0,3.829d0,3.843d0,3.855d0,3.941d0
     *,3.94d0,3.984d0,4.d0,4.074d0,3.89d0,4.111d0/
      data twsnuc/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
     *,0.55d0,0.55d0,0.56d0,0.56d0,0.5052d0,0.498d0,0.513d0
     *,0.55d0,0.55d0,0.567d0,0.698d0,0.698d0,0.549d0,0.55d0
     *,0.604d0,0.58d0,0.523d0,0.519d0,0.58d0,0.547d0,0.553d0
     *,0.582d0,0.55d0,0.55d0,0.7d0,0.599d0,0.507d0,0.588d0
     *,0.588d0,0.585d0,0.62d0,0.55d0,0.55d0,0.55d0,0.55d0
     *,0.55d0,0.55d0,0.55d0,0.588d0,0.588d0
     *,0.566d0,0.505d0,0.542d0,0.557d0,0.536d0,0.567d0,0.558d0/
      data twbnuc/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
     *,0.d0,0.d0,0.d0,0.d0,-0.18d0,0.139d0,-0.051d0,0.d0,0.d0
     *,0.d0,-0.168d0,0.d0,0.d0,0.d0,-0.249d0,-0.236d0,0.d0,0.d0
     *,0.233d0,-0.203d0,-0.078d0,-0.173d0,0.d0,0.d0,0.d0,-0.1d0
     *,0.d0,-0.13d0,-0.13d0,-0.201d0,-0.19d0,0.d0,0.d0,0.d0,0.d0
     *,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
     *,0.d0,0.d0/
      end

c-----------------------------------------------------------------------
      real function qggamfun(x)
c-----------------------------------------------------------------------
c     gamma fctn
c-----------------------------------------------------------------------
      dimension c(13)
      data c
     1/ 0.00053 96989 58808, 0.00261 93072 82746, 0.02044 96308 23590,
     2  0.07309 48364 14370, 0.27964 36915 78538, 0.55338 76923 85769,
     3  0.99999 99999 99998,-0.00083 27247 08684, 0.00469 86580 79622,
     4  0.02252 38347 47260,-0.17044 79328 74746,-0.05681 03350 86194,
     5  1.13060 33572 86556/
      qggamfun=0
      z=x
      if(x .gt. 0.0) goto1
      if(x .eq. aint(x)) goto5
      z=1.0-z
    1 f=1.0/z
      if(z .le. 1.0) goto4
      f=1.0
    2 continue
      if(z .lt. 2.0) goto3
      z=z-1.0
      f=f*z
      goto2
    3 z=z-1.0
    4 qggamfun=
     1 f*((((((c(1)*z+c(2))*z+c(3))*z+c(4))*z+c(5))*z+c(6))*z+c(7))/
     2   ((((((c(8)*z+c(9))*z+c(10))*z+c(11))*z+c(12))*z+c(13))*z+1.0)
      if(x .gt. 0.0) return
      qggamfun=3.141592653589793/(sin(3.141592653589793*x)*qggamfun)
      return
    5 write(*,10)x
   10 format(1x,'argument of gamma fctn = ',e20.5)
      stop
      end

c-------------------------------------------------------------------------------
      subroutine qgcrossc(niter,gtot,gprod,gabs,gdd,gqel,gcoh)
c-------------------------------------------------------------------------------
c nucleus-nucleus (nucleus-hydrogen) interaction cross sections
c gtot  - total cross section
c gprod - production cross section (projectile diffraction included)
c gabs  - cut pomerons cross section
c gdd   - projectile diffraction cross section
c gqel  - quasielastic (projectile nucleon knock-out) cross section
c gcoh  - coherent (elastic with respect to the projectile) cross section
c (target diffraction is not treated explicitely and contributes to
c gdd, gqel, gcoh).
c-------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter(iapmax=208,nfock=3)
      dimension wabs(28),wdd(28),wqel(28),wcoh(28)
     *,b0(28),ai(28),xa(iapmax,3),xb(iapmax,3)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr11/ b10
      common /qgarr16/ cc(nfock,3),iddp(iapmax),iddt(iapmax)
      common /arr3/    x1(7),a1(7)
      external qgran

      e1=exp(-1.d0)
      bm=rnuc(1)+rnuc(2)+2.d0*max(wsnuc(1),wsnuc(2))
      do i=1,7
       b0(15-i)=bm*sqrt((1.d0+x1(i))/2.d0)
       b0(i)=bm*sqrt((1.d0-x1(i))/2.d0)
       ai(i)=a1(i)*bm**2*5.d0*pi    !change to mb
       ai(15-i)=ai(i)
      enddo
      do i=1,7
       tp=(1.d0+x1(i))/2.d0
       tm=(1.d0-x1(i))/2.d0
       b0(14+i)=bm-log(tp)*max(wsnuc(1),wsnuc(2))
       b0(29-i)=bm-log(tm)*max(wsnuc(1),wsnuc(2))
       ai(14+i)=a1(i)*b0(14+i)/tp*10.d0*max(wsnuc(1),wsnuc(2))*pi
       ai(29-i)=a1(i)*b0(29-i)/tm*10.d0*max(wsnuc(1),wsnuc(2))*pi
      enddo
      do i=1,28
       wabs(i)=0.
       wdd(i)=0.
       wqel(i)=0.
       wcoh(i)=0.
      enddo

      do nc=1,niter
       do i=1,ia(2)
        aks=qgran(b10)
        do ic=1,nfock
         aks=aks-cc(ic,2)
         if(aks.lt.0.d0)goto 1
        enddo
        ic=nfock
1       iddt(i)=ic                          !diffractive eigenstates for targ.
       enddo
       if(ia(1).eq.1)then
        xa(1,1)=0.d0
        xa(1,2)=0.d0
        xa(1,3)=0.d0
       else
        call qggea(ia(1),xa,1)
       endif
       if(ia(2).eq.1)then
        xb(1,1)=0.d0
        xb(1,2)=0.d0
        xb(1,3)=0.d0
       else
        call qggea(ia(2),xb,2)
       endif

       do i=1,28
        call qggcr(b0(i),gabs,gdd,gqel,gcoh,xa,xb,ia(1))
        wabs(i)=wabs(i)+gabs
        wdd(i)=wdd(i)+gdd
        wqel(i)=wqel(i)+gqel
        wcoh(i)=wcoh(i)+gcoh
       enddo
      enddo

      gabs=0.
      gdd=0.
      gqel=0.
      gcoh=0.
      do i=1,28
       gabs=gabs+ai(i)*wabs(i)/niter
       gdd=gdd+ai(i)*wdd(i)/niter
       gqel=gqel+ai(i)*wqel(i)/niter
       gcoh=gcoh+ai(i)*wcoh(i)/niter
      enddo
      gprod=gabs+gdd
      gtot=gprod+gqel+gcoh
      return
      end

c-------------------------------------------------------------------------------
      subroutine qggcr(bcoll,gabs,gdd,gqel,gcoh,xa,xb,ia)
c-------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter(iapmax=208,nfock=3)
      dimension xa(iapmax,3),xb(iapmax,3),vabs(nfock)

      gabs=1.
      gdd=1.
      gqel=1.
      gcoh=1.
      do n=1,ia
       call qgv(xa(n,1)+bcoll,xa(n,2),xb,vin,vdd,vabs)
       gabs=gabs*(vdd-vin+1.d0)          !prod_n^A [sum_i c_i exp(-2chi_i(n))]
       gdd=gdd*(1.-vin)                  !prod_n^A [sum_i c_i exp(-chi_i(n))]^2
       gqel=gqel*(2.d0*dsqrt(1.d0-vin)-1.d0)
                                       !prod_n^A [sum_i c_i exp(-chi_i(n)) - 1]
       gcoh=gcoh*dsqrt(1.d0-vin)
      enddo
      gcoh=1.-2.*gcoh+gqel
      gqel=gdd-gqel
      gdd=gabs-gdd
      gabs=1.-gabs
      return
      end

c-------------------------------------------------------------------------------
      double precision function qgsect(e0n,icz,iap0,iat0)    !so18032013
c-------------------------------------------------------------------------------
c qgsect - hadron-nucleus (hadron-nucleus) particle production cross section
c e0n - lab. energy per projectile nucleon (hadron),
c icz - hadron class,
c iap - projectile mass number (1=<iap<=iapmax),
c iat - target mass number     (1=<iat<=iapmax)
c-------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wa(3),wb(3)
      common /qgarr47/ gsect(10,3,6)
      common /qgarr48/ qgsasect(10,6,6)
      common /qgarr43/ moniou
      common /qgdebug/ debug

      if(debug.ge.3)write (moniou,201)e0n,icz,iap0,iat0
      qgsect=0.d0

      iap=iap0                                              !so18032013-beg
      iat=iat0
      if(iat.eq.1.and.iap.ne.1)then
       iap=iat0
       iat=iap0
      endif                                                 !so18032013-end

      ye=max(1.d0,dlog10(e0n))
      je=min(8,int(ye))
      wk(2)=ye-je
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)

      if(iat.eq.14)then
       jb=5
       jmax=1
       wb(1)=1.d0
      elseif(iat.eq.40)then
       jb=6
       jmax=1
       wb(1)=1.d0
      else
       jmax=3
       yb=max(1.d0,dlog(dble(iat))/1.38629d0+1.d0)
       jb=min(int(yb),2)
       wb(2)=yb-jb
       wb(3)=wb(2)*(wb(2)-1.d0)*.5d0
       wb(1)=1.d0-wb(2)+wb(3)
       wb(2)=wb(2)-2.d0*wb(3)
      endif

      if(iap.eq.1)then
       do i=1,3
       do l=1,jmax
        qgsect=qgsect+gsect(je+i-1,icz,jb+l-1)*wk(i)*wb(l)
       enddo
       enddo
      else
       ya=max(1.d0,dlog(dble(iap)/2.d0)/.69315d0+1.d0)
       ja=min(int(ya),4)
       wa(2)=ya-ja
       wa(3)=wa(2)*(wa(2)-1.d0)*.5d0
       wa(1)=1.d0-wa(2)+wa(3)
       wa(2)=wa(2)-2.d0*wa(3)
       do i=1,3
       do m=1,3
       do l=1,jmax
        qgsect=qgsect+qgsasect(je+i-1,ja+m-1,jb+l-1)*wk(i)*wa(m)*wb(l)
       enddo
       enddo
       enddo
      endif
      qgsect=exp(qgsect)
      if(debug.ge.4)write (moniou,202)

201   format(2x,'qgsect - nucleus-nucleus production cross section'
     */4x,'lab. energy per nucleon - ',e10.3,2x,'hadron class - ',i2
     */4x,'proj. mass N - ',i3,2x,'targ. mass N - ',i3)
202   format(2x,'qgsect=',e10.3)
      return
      end

c=============================================================================
      subroutine qgreg(ep0,ic)
c-----------------------------------------------------------------------
c qgreg - registration of produced hadron
c ep0 - 4-momentum,
c ic  - hadron type
c output:
c nsh - number of secondary hadrons;
c esp(i,j) - 4-vectors (E,p_z,p_x,p_y) for secondary hadrons (j=1,nsh);
c ich(j) - hadron types (0 - pi0, 1 - pi+, -1 - pi-, 2 - p, -2 - pbar, 3 - n,
c -3 - nbar, 4 - K+, -4 - K-, 5 - Ks, -5 - Kl, 6 - lambda, -6 - lambdabar,
c 7 - delta++, -7 - delta++bar, 8 - delta-, -8 - delta-bar, 10 - eta,
c 16 - rho+, -16 - rho-, 17 - rho0, 18 - K*+, -18 - K*-, 19 - K*0, -19 - K*0bar,
c 20 - sigma+, -20 - sigma+bar, 21 - sigma-, -21 - sigma-bar)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nptmax=95000)
      dimension ep(4),ep0(4),ep1(4),ep2(4),ep3(4)
      common /qgarr4/  ey0(3)
      common /qgarr10/ am0,amn,amk,amsig,amlam,ameta
      common /qgarr11/ b10
      common /qgarr12/ nsh
      common /qgarr14/ esp(4,nptmax),ich(nptmax)
      common /qgarr21/ dmmin(3),dmres(3),wdres(3)
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /ebal/    ebal0(4),ebal(4)
      external qgran

      if(debug.ge.3)write (moniou,201)ic,ep0,nsh
      nsh=nsh+1

      nstprev = nsh

      if(.not.ep0(1).gt.0.d0)then
       write (moniou,201)ic,ep0,nsh
       stop'NaN!!!'
      endif
      if(ic.eq.-10.or.ic.eq.9)stop'ic=-10?!!!'
      if(nsh.gt.nptmax)stop'increase nptmax!!!'
      iab=iabs(ic)
      do i=1,4
       ep(i)=ep0(i)
      enddo

      if(iab.eq.11)then                  !pi* -> rho + pi
       am2=qgnrm(ep)
       call qgdec2(ep,ep1,ep2,am2,dmmin(1)**2,am0**2)
       if(qgran(b10).lt..5d0)then  !rho0 + pi+/-
        ich(nsh)=17   !-10
        ich(nsh+1)=ic/iab
        do i=1,4
          esp(i,nsh)=ep1(i)
          ep(i)=ep2(i)
        enddo
        nsh=nsh+1
       else      !rho+/- + pi0 -> pi+/- + 2 pi0
        ich(nsh)=16*ic/iab    !rho+/-
        ich(nsh+1)=0
        do i=1,4
          esp(i,nsh)=ep1(i)
          ep(i)=ep2(i)
        enddo
        nsh=nsh+1
       endif

      elseif(iab.eq.12.or.iab.eq.13)then       !N*
       am2=qgnrm(ep)
       if(6.d0*qgran(b10).lt.1.d0.and.dsqrt(am2).gt.dmmin(2)+am0)then !delta+pi
        call qgdec2(ep,ep1,ep2,am2,dmmin(2)**2,am0**2)
        call qgdec2(ep1,ep3,ep,dmmin(2)**2,amn**2,am0**2)
        ich(nsh)=2*ic-25*ic/iab
        ich(nsh+1)=ic-10*ic/iab
        ich(nsh+2)=-ich(nsh)
        do i=1,4
         esp(i,nsh)=ep2(i)
         esp(i,nsh+1)=ep3(i)
        enddo
        nsh=nsh+2
       else                                    !N + pi
        if(dsqrt(am2).le.amn+am0)stop'N* - mass too low!!!'
        call qgdec2(ep,ep1,ep2,am2,amn**2,am0**2)
        do i=1,4
         esp(i,nsh)=ep1(i)
         ep(i)=ep2(i)
        enddo
        if(qgran(b10).lt..4d0)then
         ich(nsh)=ic-10*ic/iab
         ich(nsh+1)=0
        else
         ich(nsh)=15*ic/iab-ic
         ich(nsh+1)=25*ic/iab-2*ic
        endif
        nsh=nsh+1
       endif

      elseif(iab.eq.14.or.iab.eq.15)then       !K1
       am2=qgnrm(ep)
       if(dsqrt(am2).gt.dmmin(1)+amk)then      !rho + K
        call qgdec2(ep,ep1,ep2,am2,dmmin(1)**2,amk**2)
        if(3.d0*qgran(b10).lt.1.d0)then  !rho0
         ich(nsh)=ic-10*ic/iab
         ich(nsh+1)=17   !-10
         do i=1,4
           esp(i,nsh)=ep2(i)
           ep(i)=ep1(i)
         enddo
         nsh=nsh+1
        else 
         ich(nsh)=19*ic/iab-ic
         ich(nsh+1)=16*ic/iab    !rho+/-
         if(iab.eq.15)ich(nsh+1)=-ich(nsh+1)
         do i=1,4
          esp(i,nsh)=ep2(i)
          ep(i)=ep1(i)
         enddo
         nsh=nsh+1
        endif
       else                                    !K* + pi
        call qgdec2(ep,ep1,ep2,am2,dmmin(3)**2,am0**2)
        if(3.d0*qgran(b10).lt.1.d0)then
         ich(nsh)=0
         ich(nsh+1)=ic+4*ic/iab    !K*
        else
         ich(nsh)=29*ic/iab-2*ic
         ich(nsh+1)=33*ic/iab-ic    !K*
        endif
        do i=1,4
         esp(i,nsh)=ep2(i)
         ep(i)=ep1(i)
        enddo
        nsh=nsh+1
       endif

      elseif(iab.eq.5)then                     !K0,K0~
       ich(nsh)=10*int(.5d0+qgran(b10))-5

      else
       ich(nsh)=ic
      endif

      do i=1,4
       esp(i,nsh)=ep(i)
      enddo

      do n=nstprev,nsh
        do i=1,4
          ep(i)=esp(i,n)
          ebal(i)=ebal(i)-ep(i)
        enddo
        call qgtran(ep,ey0,1)
        do i=1,4
          esp(i,n)=ep(i)
        enddo
      enddo

      if(debug.ge.4)write (moniou,202)

201   format(2x,'qgreg: ic=',i2,2x,'c.m. 4-momentum:',2x,4(e10.3,1x)/
     * 4x,'number of particles in the storage: ',i5)
202   format(2x,'qgreg - end')
      return
      end

c-----------------------------------------------------------------------------
      subroutine qgdec2(ep,ep1,ep2,ww,a,b)
c two particle decay
      implicit double precision (a-h,o-z)
      integer debug
      dimension ep(4),ep1(4),ep2(4),ey(3)
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /qgdebug/ debug
      external qgran

      if(debug.ge.2)write (moniou,201)ep,ww,a,b
201   format(2x,'qgdec2: 4-momentum:',2x,4(e10.3,1x)
     */4x,'ww=',e10.3,2x,'a=',e10.3,2x,'b=',e10.3)

      pl=qglam(ww,a,b)
      ep1(1)=dsqrt(pl+a)
      ep2(1)=dsqrt(pl+b)
      pl=dsqrt(pl)
      cosz=2.d0*qgran(b10)-1.d0
      pt=pl*dsqrt(1.d0-cosz**2)
      ep1(2)=pl*cosz
      call qgcs(c,s)
      ep1(3)=pt*c
      ep1(4)=pt*s
      do i=2,4
       ep2(i)=-ep1(i)
      enddo
      call qgdeft(ww,ep,ey)
      call qgtran(ep1,ey,1)
      call qgtran(ep2,ey,1)
      if(debug.ge.3)write (moniou,203)
203   format(2x,'qgdec2 - end')
      return
      end

c------------------------------------------------------------------------
      double precision function qggrv(x,qqs,icq,iq)
c------------------------------------------------------------------------
c qggrv - GRV PDFs
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr25/ ahv(3)

      qggrv=0.
      if(x.gt..99999d0.and.(qqs.ne.qt0.or.iq.ne.1.and.iq.ne.2))return

      if(icq.eq.2)then
       sq=dlog(dlog(qqs/.232d0**2)/dlog(.23d0/.232d0**2))
       if(iq.eq.0)then                                 !gluon
        alg=.524d0
        betg=1.088d0
        aag=1.742d0-.93d0*sq
        bbg=-.399d0*sq**2
        ag=7.486d0-2.185d0*sq
        bg=16.69d0-22.74d0*sq+5.779d0*sq*sq
        cg=-25.59d0+29.71d0*sq-7.296d0*sq*sq
        dg=2.792d0+2.215d0*sq+.422d0*sq*sq-.104d0*sq*sq*sq
        eg=.807d0+2.005d0*sq
        eeg=3.841d0+.361d0*sq
        qggrv=(1.d0-x)**dg*(x**aag*(ag+bg*x+cg*x**2)*log(1.d0/x)**bbg
     *  +sq**alg*exp(-eg+sqrt(eeg*sq**betg*log(1.d0/x))))
       elseif(iq.eq.1.or.iq.eq.2)then                  !u_v or d_v
        aau=.59d0-.024d0*sq
        bbu=.131d0+.063d0*sq
        auu=2.284d0+.802d0*sq+.055d0*sq*sq
        au=-.449d0-.138d0*sq-.076d0*sq*sq
        bu=.213d0+2.669d0*sq-.728d0*sq*sq
        cu=8.854d0-9.135d0*sq+1.979d0*sq*sq
        du=2.997d0+.753d0*sq-.076d0*sq*sq
        uv=auu*x**aau*(1.d0+au*x**bbu+bu*x+cu*x**1.5d0)
        if(qqs.ne.qt0)uv=uv*(1.d0-x)**du

        aad=.376d0
        bbd=.486d0+.062d0*sq
        add=.371d0+.083d0*sq+.039d0*sq*sq
        ad=-.509d0+3.31d0*sq-1.248d0*sq*sq
        bd=12.41d0-10.52d0*sq+2.267d0*sq*sq
        ccd=6.373d0-6.208d0*sq+1.418d0*sq*sq
        dd=3.691d0+.799d0*sq-.071d0*sq*sq
        dv=add*x**aad*(1.d0+ad*x**bbd+bd*x+ccd*x**1.5d0)
        if(qqs.ne.qt0)then
         dv=dv*(1.d0-x)**dd
        elseif(x.gt..99999d0)then
         dv=0.d0
        else
         dv=dv*(1.d0-x)**(dd-ahv(2))
        endif
        if(iq.eq.1)then                              !u_v
         qggrv=uv
        elseif(iq.eq.2)then                          !d_v
         qggrv=dv
        endif

       elseif(iq.eq.-3)then                           !s_sea
        als=.914
        bets=.577
        aas=1.798-.596*sq
        as=-5.548+3.669*sqrt(sq)-.616*sq
        bs=18.92-16.73*sqrt(sq)+5.168*sq
        ds=6.379-.35*sq+.142*sq*sq
        es=3.981+1.638*sq
        ees=6.402
        qggrv=(1.-x)**ds*sq**als/log(1./x)**aas*(1.+as*sqrt(x)
     *  +bs*x)*exp(-es+sqrt(ees*sq**bets*log(1./x)))
       elseif(iabs(iq).lt.3)then                      !u_sea or d_sea
        aadel=.409-.005*sq
        bbdel=.799+.071*sq
        addel=.082+.014*sq+.008*sq*sq
        adel=-38.07+36.13*sq-.656*sq*sq
        bdel=90.31-74.15*sq+7.645*sq*sq
        ccdel=0.
        ddel=7.486+1.217*sq-.159*sq*sq
        delv=addel*x**aadel*(1.-x)**ddel
     *  *(1.+adel*x**bbdel+bdel*x+ccdel*x**1.5)

        alud=1.451
        betud=.271
        aaud=.41-.232*sq
        bbud=.534-.457*sq
        aud=.89-.14*sq
        bud=-.981
        cud=.32+.683*sq
        dud=4.752+1.164*sq+.286*sq*sq
        eud=4.119+1.713*sq
        eeud=.682+2.978*sq
        udsea=(1.-x)**dud*(x**aaud*(aud+bud*x+cud*x**2)
     *  *log(1./x)**bbud+sq**alud*exp(-eud+sqrt(eeud*sq**betud
     *  *log(1./x))))

        if(iq.eq.-1)then                           !u_sea
         qggrv=(udsea-delv)/2.
        elseif(iq.eq.-2)then                       !d_sea
         qggrv=(udsea+delv)/2.
        endif
       else
        qggrv=0.
       endif

      elseif(icq.eq.1.or.icq.eq.3)then
       sq=dlog(dlog(qqs/.204d0**2)/dlog(.26d0/.204d0**2))
       if(iq.eq.1.or.iq.eq.2)then
        aapi=.517-.02*sq
        api=-.037-.578*sq
        bpi=.241+.251*sq
        dpi=.383+.624*sq
        anorm=1.212+.498*sq+.009*sq**2
        qggrv=.5*anorm*x**aapi*(1.+api*sqrt(x)+bpi*x)
        if(qqs.ne.qt0)qggrv=qggrv*(1.d0-x)**dpi
       elseif(iq.eq.0)then
          alfpi=.504
          betpi=.226
          aapi=2.251-1.339*sqrt(sq)
          api=2.668-1.265*sq+.156*sq**2
          bbpi=0.
          bpi=-1.839+.386*sq
          cpi=-1.014+.92*sq-.101*sq**2
          dpi=-.077+1.466*sq
          epi=1.245+1.833*sq
          eppi=.51+3.844*sq
          qggrv=(1.-x)**dpi*(x**aapi*(api+bpi*sqrt(x)+cpi*x)*
     *    log(1./x)**bbpi+sq**alfpi*
     *    exp(-epi+sqrt(eppi*sq**betpi*log(1./x))))
        elseif(iq.eq.-3)then
          alfpi=.823
          betpi=.65
          aapi=1.036-.709*sq
          api=-1.245+.713*sq
          bpi=5.58-1.281*sq
          dpi=2.746-.191*sq
          epi=5.101+1.294*sq
          eppi=4.854-.437*sq
          qggrv=sq**alfpi/log(1./x)**aapi*(1.-x)**dpi*
     *    (1.+api*sqrt(x)+bpi*x)*
     *    exp(-epi+sqrt(eppi*sq**betpi*log(1./x)))
        elseif(iabs(iq).lt.3)then
          alfpi=1.147
          betpi=1.241
          aapi=.309-.134*sqrt(sq)
          api=.219-.054*sq
          bbpi=.893-.264*sqrt(sq)
          bpi=-.593+.24*sq
          cpi=1.1-.452*sq
          dpi=3.526+.491*sq
          epi=4.521+1.583*sq
          eppi=3.102
          qggrv=(1.-x)**dpi*(x**aapi*(api+bpi*sqrt(x)+cpi*x)*
     *    log(1./x)**bbpi+sq**alfpi*
     *    exp(-epi+sqrt(eppi*sq**betpi*log(1./x))))
        else
          qggrv=0.
        endif
      else
       qggrv=0.
      endif
      return
      end

c------------------------------------------------------------------------
      double precision function qgev(q1,qj,qq,xx,j,l)
c------------------------------------------------------------------------
c qgev - PDF evolution
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr51/ epsxmn
      common /arr3/    x1(7),a1(7)

      qgev=0.d0
      zmax=1.d0-epsxmn
      zmin=xx/zmax
      if(zmin.ge.zmax)return

      if(qj.eq.qq)then
       do i1=1,7
       do m1=1,2
        qi=q1*(qq/q1)**(.5d0+x1(i1)*(m1-1.5d0))

        fz1=0.d0
        fz2=0.d0
        fz3=0.d0
        zmin1=max(.2d0,zmin)
        zmax1=min(.2d0,zmax)
        zmax1=min(5.d0*xx,zmax1)
        zmax2=min(zmin1,zmax)
        zmin2=max(zmax1,zmin)

        if(zmax1.gt.zmin)then
         do i=1,7
         do m=1,2
          z=xx+(zmin-xx)*((zmax1-xx)/(zmin-xx))**(.5d0+(m-1.5d0)*x1(i))
          do k=1,2
           if(j.ne.3.or.k.ne.1)then
            fz1=fz1+a1(i)*qgevi(q1,qi,xx/z,j,k)*qgfap(z,k,l)*(1.d0-xx/z)
           endif
          enddo
         enddo
         enddo
         fz1=fz1*dlog((zmax1-xx)/(zmin-xx))
        endif
        if(zmin1.lt.zmax)then
         do i=1,7
         do m=1,2
          z=1.d0-(1.d0-zmax)*((1.d0-zmin1)/(1.d0-zmax))
     *    **(.5d0+x1(i)*(m-1.5d0))
          do k=1,2
           if(j.ne.3.or.k.ne.1)then
            fz2=fz2+a1(i)*qgevi(q1,qi,xx/z,j,k)*qgfap(z,k,l)
     *      *(1.d0/z-1.d0)
           endif
          enddo
         enddo
         enddo
         fz2=fz2*dlog((1.d0-zmin1)/(1.d0-zmax))
        endif
        if(zmax2.gt.zmin2)then
         do i=1,7
         do m=1,2
          z=zmin2*(zmax2/zmin2)**(.5d0+x1(i)*(m-1.5d0))
          do k=1,2
           if(j.ne.3.or.k.ne.1)then
            fz3=fz3+a1(i)*qgevi(q1,qi,xx/z,j,k)*qgfap(z,k,l)
           endif
          enddo
         enddo
         enddo
         fz3=fz3*dlog(zmax2/zmin2)
        endif
        qgev=qgev+a1(i1)*(fz1+fz2+fz3)/qgsudx(qi,l)*qgalf(qi/alm)
       enddo
       enddo
       qgev=qgev*dlog(qq/q1)/4.d0*qgsudx(qq,l)

      else
       fz1=0.d0
       fz2=0.d0
       fz3=0.d0
       zmin1=max(.2d0,zmin)
       zmax1=min(.2d0,zmax)
       zmax1=min(5.d0*xx,zmax1)
       zmax2=min(zmin1,zmax)
       zmin2=max(zmax1,zmin)

       if(zmax1.gt.zmin)then
        do i=1,7
        do m=1,2
         z=xx+(zmin-xx)*((zmax1-xx)/(zmin-xx))**(.5d0+(m-1.5d0)*x1(i))
         do k=1,2
          if(j.ne.3)then
           fz1=fz1+a1(i)*qgevi(q1,qj,xx/z,j,k)*qgevi(qj,qq,z,k,l)
     *     *(1.d0-xx/z)
          elseif(k.ne.1)then
           fz1=fz1+a1(i)*qgevi(q1,qj,xx/z,3,2)*qgevi(qj,qq,z,3,2)
     *     *(1.d0-xx/z)
          endif
         enddo
        enddo
        enddo
        fz1=fz1*dlog((zmax1-xx)/(zmin-xx))
       endif
       if(zmin1.lt.zmax)then
        do i=1,7
        do m=1,2
         z=1.d0-(1.d0-zmax)*((1.d0-zmin1)/(1.d0-zmax))
     *   **(.5d0+x1(i)*(m-1.5d0))
         do k=1,2
          if(j.ne.3)then
           fz2=fz2+a1(i)*qgevi(q1,qj,xx/z,j,k)*qgevi(qj,qq,z,k,l)
     *     *(1.d0/z-1.d0)
          elseif(k.ne.1)then
           fz2=fz2+a1(i)*qgevi(q1,qj,xx/z,3,2)*qgevi(qj,qq,z,3,2)
     *     *(1.d0/z-1.d0)
          endif
         enddo
        enddo
        enddo
        fz2=fz2*dlog((1.d0-zmin1)/(1.d0-zmax))
       endif
       if(zmax2.gt.zmin2)then
        do i=1,7
        do m=1,2
         z=zmin2*(zmax2/zmin2)**(.5d0+x1(i)*(m-1.5d0))
         do k=1,2
          if(j.ne.3)then
           fz2=fz2+a1(i)*qgevi(q1,qj,xx/z,j,k)*qgevi(qj,qq,z,k,l)
          elseif(k.ne.1)then
           fz2=fz2+a1(i)*qgevi(q1,qj,xx/z,3,2)*qgevi(qj,qq,z,3,2)
          endif
         enddo
        enddo
        enddo
        fz3=fz3*dlog(zmax2/zmin2)
       endif
       qgev=(fz1+fz2+fz3)/2.d0
      endif
      return
      end

c------------------------------------------------------------------------
      double precision function qgevi(q1,qq,xx,m,l)
c------------------------------------------------------------------------
c qgevi - PDF evolution - interpolation
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension wi(3),wj(3),wk(3)
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr20/ spmax
      common /qgarr51/ epsxmn
      common /qgarr52/ evk(40,40,100,3,2)

      qgevi=0.d0
      if(q1.ge..9999d0*spmax)goto 1

      if(xx.le..1d0)then
       yx=37.d0-dlog(.1d0/xx)/dlog(.1d0*spmax)*36.d0
       k=max(1,int(yx))
       k=min(k,35)
      elseif(xx.le..9d0)then
       yx=(xx-.1d0)*40.d0+37.d0
       k=max(37,int(yx))
       k=min(k,67)
      else
       yx=dlog(10.d0*(1.d0-xx))/log(10.d0*epsxmn)*31.d0+69.d0
       k=max(69,int(yx))
       k=min(k,98)
      endif
      wk(2)=yx-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)

      qli=log(q1)/dlog(spmax)*39.d0+1.d0
      qlj=log(qq/q1)/dlog(spmax/q1)*39.d0+1.d0
      i=max(1,int(1.0001d0*qli))
      i=min(i,38)
      wi(2)=qli-i
      wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
      wi(1)=1.d0-wi(2)+wi(3)
      wi(2)=wi(2)-2.d0*wi(3)

      j=max(1,int(1.0001d0*qlj))
      j=min(j,38)
      wj(2)=qlj-j
      wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
      wj(1)=1.d0-wj(2)+wj(3)
      wj(2)=wj(2)-2.d0*wj(3)

      do i1=1,3
      do j1=1,3
      do k1=1,3
       k2=k+k1-1
       qgevi=qgevi+evk(i+i1-1,j+j1-1,k2,m,l)*wi(i1)*wj(j1)*wk(k1)
      enddo
      enddo
      enddo
1     qgevi=exp(qgevi)*qgfap(xx,m,l)
      if(m.eq.1.and.l.eq.1.or.m.ne.1.and.l.ne.1)then
       qgevi=qgevi/4.5d0/qgsudx(q1,m)*qgsudx(qq,m)
     * *dlog(dlog(qq/alm)/dlog(q1/alm))
      else
       qgevi=qgevi*.3d0/(dlog(epsxmn)+.75d0)
     * *(qgsudx(qq,1)/qgsudx(q1,1)-qgsudx(qq,2)/qgsudx(q1,2))
      endif
      return
      end

c=============================================================================
      double precision function qgjeto(qi,qj,s,iq1,iq2)
c-----------------------------------------------------------------------------
c qgjeto - hard 2->2 parton scattering born cross-section
c s is the c.m. energy square for the scattering process,
c iq1 - parton type at current end of the ladder (0 - g, 1,2 etc. - q)
c iq2 - parton type at opposite end of the ladder (0 - g, 1,2 etc. - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr51/ epsxmn
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)qi,qj,s,iq1,iq2

      qgjeto=0.d0
      qq=max(qi,qj)

      zmin=qq*fqscal*4.d0/s
      zmax=1.d0-epsxmn
      if(zmin.ge.zmax)return

      dpx1=0.d0
      zmin1=min(.2d0,1.d0-zmin)
      do i1=1,7
      do m1=1,2
       z=1.d0-epsxmn*(zmin1/epsxmn)**(.5d0+x1(i1)*(m1-1.5d0))

       si=z*s
       fb=qgjeti(qi,qj,si,z,1.d0,iq1,iq2,1)
       dpx1=dpx1+a1(i1)*fb*(1.d0-z)
      enddo
      enddo
      dpx1=dpx1*dlog(zmin1/epsxmn)

      dpx2=0.d0
      if(zmin.lt..8d0)then
       zmin1=zmin**(-delh)
       zmax1=.8d0**(-delh)
       do i1=1,7
       do m1=1,2
        z=(.5d0*(zmax1+zmin1+(zmax1-zmin1)*x1(i1)*(2*m1-3)))
     *  **(-1.d0/delh)

        si=z*s
        fb=qgjeti(qi,qj,si,z,1.d0,iq1,iq2,1)
        dpx2=dpx2+a1(i1)*fb*z**(1.d0+delh)
       enddo
       enddo
       dpx2=dpx2*(zmin1-zmax1)/delh
      endif
      qgjeto=(dpx1+dpx2)/qgsudx(qj,min(2,iabs(iq2)+1))*pi**3

      if(debug.ge.3)write (moniou,202)qgjeto
201   format(2x,'qgjeto: qi=',e10.3,2x,'qj=',e10.3,2x,
     *'s= ',e10.3,2x,'iq1= ',i1,2x,'iq2= ',i1)
202   format(2x,'qgjeto=',e10.3)
      return
      end

c=============================================================================
      double precision function qgjett(qi,qj,s,iq1,iq2)
c-----------------------------------------------------------------------------
c qgjett - hard 2->2 parton scattering born cross-section
c s is the c.m. energy square for the scattering process,
c iq1 - parton type at current end of the ladder (0 - g, 1,2 etc. - q)
c iq2 - parton type at opposite end of the ladder (0 - g, 1,2 etc. - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr51/ epsxmn
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(debug.ge.2)write (moniou,201)qi,qj,s,iq1,iq2

      qgjett=0.d0
      qq=max(qi,qj)

      zmin=qq*fqscal*4.d0/s
      zmax=(1.d0-epsxmn)**2
      if(zmin.ge.zmax)return
      zmin1=zmin**(-delh)
      zmax1=zmax**(-delh)
      do i1=1,7
      do m1=1,2
       z=(.5d0*(zmax1+zmin1+(zmax1-zmin1)*x1(i1)*(2*m1-3)))
     * **(-1.d0/delh)

       si=z*s
       fb1=0.d0
       zmin2=min(.2d0,1.d0-dsqrt(z))
       do i2=1,7
       do m2=1,2
        z1=1.d0-epsxmn*(zmin2/epsxmn)**(.5d0+x1(i2)*(m2-1.5d0))
        z2=z/z1
        fb1=fb1+a1(i2)*(qgjeti(qi,qj,si,z1,z2,iq1,iq2,2)
     *  +qgjeti(qi,qj,si,z2,z1,iq1,iq2,2))*(1.d0/z1-1.d0)
       enddo
       enddo
       fb1=fb1*dlog(zmin2/epsxmn)

       fb2=0.d0
       if(z.lt..64d0)then
        do i2=1,7
        do m2=1,2
         z1=.8d0*(dsqrt(z)/.8d0)**(.5d0+x1(i2)*(m2-1.5d0))
         z2=z/z1
         fb2=fb2+a1(i2)*(qgjeti(qi,qj,si,z1,z2,iq1,iq2,2)
     *   +qgjeti(qi,qj,si,z2,z1,iq1,iq2,2))
        enddo
        enddo
        fb2=fb2*dlog(.64d0/z)/2.d0
       endif

       qgjett=qgjett+a1(i1)*(fb1+fb2)*z**(1.d0+delh)
      enddo
      enddo
      qgjett=qgjett*(zmin1-zmax1)/delh*pi**3/2.d0

      if(debug.ge.3)write (moniou,202)qgjett
201   format(2x,'qgjett: qi=',e10.3,2x,'qj=',e10.3,2x,
     *'s= ',e10.3,2x,'iq1= ',i1,2x,'iq2= ',i1)
202   format(2x,'qgjett=',e10.3)
      return
      end

c=============================================================================
      double precision function qgjeti(qi,qj,si,z1,z2,iq1,iq2,jj)
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      dimension sigsi(2)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr51/ epsxmn
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      qgjeti=0.d0
      qq=max(qi,qj)
      tmin=qq*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qq*fqscal/si)))
      if(tmin.ge.si/2.d0)return
      do i=1,7
      do m=1,2
       t=2.d0*tmin/(1.d0+2.d0*tmin/si
     *   -x1(i)*(2*m-3)*(1.d0-2.d0*tmin/si))
       qt=t*(1.d0-t/si)

       fb=0.d0
       if(jj.eq.1)then
        do iql=1,2
         iq=2*iql-2
         sigs=0.d0
         do n=1,3
          sigs=sigs+qgfbor(si,t,iq,min(1,iabs(iq2)),n)
     *    +qgfbor(si,si-t,iq,min(1,iabs(iq2)),n)
         enddo
         if(iq.eq.iq2)sigs=sigs/2.d0
         evs=qgevi(qi,qt/fqscal,z1,min(2,iabs(iq1)+1),iql)
         fb=fb+evs*sigs
        enddo
        if(iq2.ne.0)then
         do iql=1,2
          iq=3-2*iql
          sigsi(iql)=0.d0
          do n=1,3
           sigsi(iql)=sigsi(iql)+qgfbor(si,t,iq,1,n)
     *     +qgfbor(si,si-t,iq,1,n)
          enddo
         enddo
         fb=fb+evs*(sigsi(1)/12.d0+sigsi(2)/6.d0-sigs/3.d0)
         if(iq1.ne.0)then
          evns=qgevi(qi,qt/fqscal,z1,3,2)
          if(max(iq1,iq2).lt.2)then
           fb=fb+evns*(sigs/3.d0-sigsi(1)/12.d0-sigsi(2)/6.d0)
          elseif(max(iq1,iq2).eq.2)then
           fb=fb+evns*(sigsi(1)/2.4d0-sigsi(2)/6.d0-sigs/1.5d0)
          elseif(max(iq1,iq2).eq.3)then
           fb=fb+evns*(sigsi(2)/1.2d0-sigsi(1)/12.d0-sigs/1.5d0)
          else
           stop'problem with parton types in qgjeti'
          endif
         endif
        endif
        fb=fb*qgsudx(qt/fqscal,min(2,iabs(iq2)+1))

       else
        do iql=1,2
         iq=2*iql-2
         evsl=qgevi(qi,qt/fqscal,z1,min(2,iabs(iq1)+1),iql)
        do iqr=1,2
         evsr=qgevi(qj,qt/fqscal,z2,min(2,iabs(iq2)+1),iqr)
         sigs=0.d0
         do n=1,3
          sigs=sigs+qgfbor(si,t,iq,iqr-1,n)+qgfbor(si,si-t,iq,iqr-1,n)
         enddo
         if(iq.eq.iqr-1)sigs=sigs/2.d0
         fb=fb+evsl*evsr*sigs
        enddo
        enddo
        do iql=1,2
         iq=3-2*iql
         sigsi(iql)=0.d0
         do n=1,3
          sigsi(iql)=sigsi(iql)+qgfbor(si,t,iq,1,n)
     *    +qgfbor(si,si-t,iq,1,n)
         enddo
        enddo
        fb=fb+evsl*evsr*(sigsi(1)/12.d0+sigsi(2)/6.d0-sigs/3.d0)
        if(iq1.ne.0.and.iq2.ne.0)then
         evnsl=qgevi(qi,qt/fqscal,z1,3,2)
         evnsr=qgevi(qj,qt/fqscal,z2,3,2)
         if(max(iq1,iq2).lt.2)then
          fb=fb+evnsl*evnsr*(sigs/3.d0-sigsi(1)/12.d0-sigsi(2)/6.d0)
         elseif(max(iq1,iq2).eq.2)then
          fb=fb+evnsl*evnsr*(sigsi(1)/2.4d0-sigsi(2)/6.d0-sigs/1.5d0)
         elseif(max(iq1,iq2).eq.3)then
          fb=fb+evns*(sigsi(2)/1.2d0-sigsi(1)/12.d0-sigs/1.5d0)
         else
          stop'problem with parton types in qgjeti'
         endif
        endif
       endif
       qgjeti=qgjeti+a1(i)*fb*qgalf(qt/fqscal/alm)**2*t**2
      enddo
      enddo
      qgjeti=qgjeti*(1.d0/tmin-2.d0/si)/si**2
      return
      end

c------------------------------------------------------------------------
      double precision function qgpdf(xx,qq,icz,jj)
c-----------------------------------------------------------------------
c qgpdf - parton distribution function for proton
c qq  - virtuality scale,
c xx  - light cone x,
c icz - hadron type,
c jj  - parton type (0 - gluon, 1 - u_v, 2 - d_v, -1 - q_sea)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr25/ ahv(3)
      common /qgarr51/ epsxmn
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      if(jj.eq.0)then
       qgpdf=qgpdfi(xx,icz,1)
      elseif(jj.eq.1.or.jj.eq.2)then
       qgpdf=qggrv(xx,qt0,icz,jj)*(1.d0-xx)**ahv(icz)
      else
       qgpdf=qgpdfi(xx,icz,2)-(qggrv(xx,qt0,icz,1)
     * +qggrv(xx,qt0,icz,2))*(1.d0-xx)**ahv(icz)
      endif
      qgpdf=qgpdf*qgsudx(qq,iabs(jj)+1)/qgsudx(qt0,iabs(jj)+1)

      xmin=xx/(1.d0-epsxmn)
      if(xmin.lt.1.d0.and.qq.gt.qt0)then
       dpd1=0.d0
       dpd2=0.d0
       xm=max(xmin,.3d0)
       do i=1,7         !numerical integration over zx
       do m=1,2
        zx=1.d0-(1.d0-xm)*(.5d0+(m-1.5d0)*x1(i))**.25d0
        z=xx/zx

        gl=qgpdfi(zx,icz,1)
        uv=qggrv(zx,qt0,icz,1)*(1.d0-zx)**ahv(icz)
        dv=qggrv(zx,qt0,icz,2)*(1.d0-zx)**ahv(icz)
        sea=qgpdfi(zx,icz,2)-uv-dv
        if(jj.eq.0)then
         fz=qgevi(qt0,qq,z,1,1)*gl+qgevi(qt0,qq,z,2,1)*(uv+dv+sea)
        elseif(jj.eq.1)then
         fz=qgevi(qt0,qq,z,3,2)*uv
        elseif(jj.eq.2)then
         fz=qgevi(qt0,qq,z,3,2)*dv
        else
         akns=qgevi(qt0,qq,z,3,2)              !nonsinglet contribution
         aks=(qgevi(qt0,qq,z,2,2)-akns)        !singlet contribution
         fz=(qgevi(qt0,qq,z,1,2)*gl+aks*(uv+dv+sea)+akns*sea)
        endif
        dpd1=dpd1+a1(i)*fz/zx**2/(1.d0-zx)**3
       enddo
       enddo
       dpd1=dpd1*(1.d0-xm)**4/8.d0*xx

       if(xm.gt.xmin)then
        do i=1,7         !numerical integration
        do m=1,2
         zx=xx+(xm-xx)*((xmin-xx)/(xm-xx))**(.5d0-(m-1.5d0)*x1(i))
         z=xx/zx

         gl=qgpdfi(zx,icz,1)
         uv=qggrv(zx,qt0,icz,1)*(1.d0-zx)**ahv(icz)
         dv=qggrv(zx,qt0,icz,2)*(1.d0-zx)**ahv(icz)
         sea=qgpdfi(zx,icz,2)-uv-dv
         if(jj.eq.0)then
          fz=qgevi(qt0,qq,z,1,1)*gl+qgevi(qt0,qq,z,2,1)*(uv+dv+sea)
         elseif(jj.eq.1)then
          fz=qgevi(qt0,qq,z,3,2)*uv
         elseif(jj.eq.2)then
          fz=qgevi(qt0,qq,z,3,2)*dv
         else
          akns=qgevi(qt0,qq,z,3,2)              !nonsinglet contribution
          aks=(qgevi(qt0,qq,z,2,2)-akns)        !singlet contribution
          fz=(qgevi(qt0,qq,z,1,2)*gl+aks*(uv+dv+sea)+akns*sea)
         endif
         dpd2=dpd2+a1(i)*fz*(1.d0-xx/zx)/zx
        enddo
        enddo
        dpd2=dpd2*dlog((xm-xx)/(xmin-xx))*.5d0*xx
       endif
       qgpdf=qgpdf+dpd2+dpd1
      endif
      return
      end

c-------------------------------------------------------------------------------
      double precision function qgfact(sy)
c-------------------------------------------------------------------------------
c qgsect - hadron-nucleus (hadron-nucleus) particle production cross section
c e0n - lab. energy per projectile nucleon (hadron),
c icz - hadron class,
c iap - projectile mass number (1=<iap<=iapmax),
c iat - target mass number     (1=<iat<=iapmax)
c-------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      
      k=0
      f=1.d0
      qgfact=1.d0
      
1     k=k+1
      if(dlog(sy).gt.dlog(sgap)*(k+1))then
       f=-f*dels/k
       qgfact=qgfact+f*(dlog(sy)-dlog(sgap)*(k+1))**k
       goto 1
      endif
      return
      end

c=============================================================================
      double precision function qgloopi(sy,vvx,iq,jj)
c-----------------------------------------------------------------------
c qgloopi - integrated loop contribution
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr4/    x4(2),a4(2)
      
      qgloopi=0.d0
      if(sy.le.1.d0)return
      rh=4.d0*.0389d0*alfp*dlog(sy)
      do ib=1,2
      do mb=1,2
       zz=.5d0+x4(ib)*(mb-1.5d0)
       bb=-rh*dlog(zz)
       qgloopi=qgloopi+a4(ib)/zz*qgloopri(sy,bb,vvx,iq,jj)
      enddo
      enddo
      qgloopi=qgloopi*pi*rh/.0389d0/2.d0
      return 
      end
 
c------------------------------------------------------------------------
      subroutine qg3pht(sig3ht,sy,b,vvx,vvxp,vvxt,genh
     *,icdp,icdt,iczp,iczt)
c-----------------------------------------------------------------------
c qg3pht - eikonals for single hard scattering contribution (ht corrections)
c sy   - pomeron mass squared,
c b    - impact parameter squared,
c vvx  = 1 - exp[-sum_{j<J} chi_targ(j) - sum_{i<I} chi_proj(i)]
c vvxp = 1 - exp[-sum_{i>I} chi_proj(i)] 
c vvxt = 1 - exp[-sum_{j>J} chi_targ(j)] 
c icdp - diffractive state for the projectile,
c icdt - diffractive state for the target,
c iczp - projectile class,
c iczt - target class,
c sig3ht: 1 - g-uncut(no HT), 2 - g-uncut(HT), 3 - g-cut(no HT), 4 - g-cut(HT),
c 5 - qv-uncut(no HT), 6 - qv-uncut(HT), 7 - qv-cut(no HT), 8 - qv-cut(HT)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      dimension sigh(2,2),pdfp(2,2),pdft(3,2),sig3ht(8)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)
      common /arr9/    x9(3),a9(3)

      do i=1,8
       sig3ht(i)=0.d0
      enddo
      s2min=4.d0*fqscal*qt0
      xmin=s2min/sy
      if(xmin.ge.1.d0)return
      xmin=xmin**(delh-dels)
      
      if(vvx+vvxp+vvxt.gt.0.d0.or.genh.ne.1.d0)then
       ixxm=3
      else
       ixxm=7
      endif
      do ix1=1,ixxm
      do mx1=1,2
       if(ixxm.eq.7)then
        zpm=(.5d0*(1.d0+xmin-(2*mx1-3)*x1(ix1)*(1.d0-xmin)))
     *  **(1.d0/(delh-dels))
       else
        zpm=(.5d0*(1.d0+xmin-(2*mx1-3)*x9(ix1)*(1.d0-xmin)))
     *  **(1.d0/(delh-dels))
       endif
       ww=zpm*sy
       do iq1=1,2
       do iq2=1,2
        sigh(iq1,iq2)=qgjit(qt0,qt0,ww,iq1,iq2)
       enddo
       enddo
       do ix2=1,ixxm
       do mx2=1,2
        if(ixxm.eq.7)then
         xx=.5d0*(1.d0+x1(ix2)*(2*mx2-3))
        else
         xx=.5d0*(1.d0+x9(ix2)*(2*mx2-3))
        endif
        xp=zpm**xx
        xm=zpm/xp
        rp1=(rq(icdp,iczp)-alfp*log(xp))*4.d0*.0389d0
        rp2=(rq(icdt,iczt)-alfp*log(xm))*4.d0*.0389d0
        rp=rp1*rp2/(rp1+rp2)
        
        do ib1=1,ixxm
        do mb1=1,2
         if(ixxm.eq.7)then
          z=.5d0+x1(ib1)*(mb1-1.5d0)
         else
          z=.5d0+x9(ib1)*(mb1-1.5d0)
         endif
         bb0=-rp*dlog(z)
        do ib2=1,ixxm
        do mb2=1,2
         if(ixxm.eq.7)then
          phi=pi*(.5d0+x1(ib2)*(mb2-1.5d0))
         else
          phi=pi*(.5d0+x9(ib2)*(mb2-1.5d0))
         endif
         bbp=(b*rp1/(rp1+rp2)+dsqrt(bb0)*cos(phi))**2+bb0*sin(phi)**2
         bbt=(b*rp2/(rp1+rp2)-dsqrt(bb0)*cos(phi))**2+bb0*sin(phi)**2
               
         if(xp*sgap.lt.1.d0)then
          vpf0=qgfani(1.d0/xp,bbp,1.d0-(1.d0-vvx)*(1.d0-vvxt)
     *    ,0.d0,0.d0,icdp,iczp,3)
         else
          vpf0=0.d0
         endif
         if(xp*sy.gt.sgap)then
          vtf0=qgfani(xp*sy,bbt,1.d0-(1.d0-vvx)*(1.d0-vvxp)
     *    ,0.d0,0.d0,icdt,iczt,3)
         else
          vtf0=0.d0
         endif
         n=1
1        n=n+1
         if(xp*sgap.lt.1.d0)then
          vpf=qgfani(1.d0/xp,bbp,1.d0-(1.d0-vvx)*(1.d0-vvxt)*exp(-vtf0)
     *    ,0.d0,0.d0,icdp,iczp,3)
         else
          vpf=0.d0
         endif
         if(xp*sy.gt.sgap)then
          vtf=qgfani(xp*sy,bbt,1.d0-(1.d0-vvx)*(1.d0-vvxp)*exp(-vpf0)
     *    ,0.d0,0.d0,icdt,iczt,3)
         else
          vtf=0.d0
         endif
         if(vpf0.gt.0.d0.and.vtf0.gt.0.d0.and.abs(1.d0-vpf/vpf0)
     *   +abs(1.d0-vtf/vtf0).gt.1.d-2.and.n.lt.100)then
          vpf0=vpf
          vtf0=vtf
          goto 1
         endif
         vvxtg=1.d0-exp(-vtf)*(1.d0-vvx)*(1.d0-vvxt)
         do jj=1,2
         do iq1=1,2
          pdfp(iq1,jj)=qgpdfbi(xp,bbp,vvxtg,vvxp,icdp,iczp,iq1,jj)
         enddo
         enddo
         
         if(xm*sy.gt.sgap)then
          vpf0=qgfani(xm*sy,bbp,1.d0-(1.d0-vvx)*(1.d0-vvxt)
     *    ,0.d0,0.d0,icdp,iczp,3)
         else
          vpf0=0.d0
         endif
         if(xm*sgap.lt.1.d0)then
          vtf0=qgfani(1.d0/xm,bbt,1.d0-(1.d0-vvx)*(1.d0-vvxp)
     *    ,0.d0,0.d0,icdt,iczt,3)
         else
          vtf0=0.d0
         endif
         n=1
2        n=n+1
         if(xm*sy.gt.sgap)then
          vpf=qgfani(xm*sy,bbp,1.d0-(1.d0-vvx)*(1.d0-vvxt)*exp(-vtf0)
     *    ,0.d0,0.d0,icdp,iczp,3)
         else
          vpf=0.d0
         endif
         if(xm*sgap.lt.1.d0)then
          vtf=qgfani(1.d0/xm,bbt,1.d0-(1.d0-vvx)*(1.d0-vvxp)
     *    *exp(-vpf0),0.d0,0.d0,icdt,iczt,3)
         else
          vtf=0.d0
         endif
         if(vpf0.gt.0.d0.and.vtf0.gt.0.d0.and.abs(1.d0-vpf/vpf0)
     *   +abs(1.d0-vtf/vtf0).gt.1.d-2.and.n.lt.100)then
          vpf0=vpf
          vtf0=vtf
          goto 2
         endif
         vvxpr=1.d0-exp(-vpf)*(1.d0-vvx)*(1.d0-vvxp)
         do jj=1,2
         do iq2=1,3
          pdft(iq2,jj)=qgpdfbi(xm,bbt,vvxpr,vvxt,icdt,iczt,iq2,jj)
         enddo
         enddo

         do jj=1,2
          if(jj.eq.1)then
           vvxs=vvxtg
          else
           vvxs=1.d0-(1.d0-vvxtg)*(1.d0-vvxp)
          endif
         do jht=1,2
         do jv=1,2
          sigt=0.d0
          do iq1=1,2
          do iq2=jv,2
           sigg=sigh(iq1,iq2)*pdfp(iq1,jj)*pdft(iq2+jv-1,jj)
           if(jht.eq.2)sigg=sigg*qgfhti(xm*sy,xp,bbp,vvxs,genh
     *     ,iq1,iq2,icdp,iczp,jj)
           sigt=sigt+sigg
          enddo
          enddo
          if(ixxm.eq.7)then
           sig3ht(jht+2*(jj-1)+4*(jv-1))=sig3ht(jht+2*(jj-1)+4*(jv-1))
     *     -a1(ib1)*a1(ib2)*a1(ix1)*a1(ix2)/z*rp*sigt
     *     *dlog(zpm)/zpm**(delh-dels)
          else
           sig3ht(jht+2*(jj-1)+4*(jv-1))=sig3ht(jht+2*(jj-1)+4*(jv-1))
     *     -a9(ib1)*a9(ib2)*a9(ix1)*a9(ix2)/z*rp*sigt
     *     *dlog(zpm)/zpm**(delh-dels)
          endif
         enddo
         enddo
         enddo
        enddo
        enddo
        enddo
        enddo
       enddo
       enddo
      enddo
      enddo
      do i=1,8
       sig3ht(i)=sig3ht(i)/32.d0*(1.d0-xmin)/(delh-dels)*factk*pi !1/2 included
     * /.0389d0
      enddo
      return 
      end

c=============================================================================
      double precision function qgfht(sm,xp,bbp,vvx,iq1,iq2,icdp,icz
     *,jj)
c-----------------------------------------------------------------------------
c qgfht - HT correction to hard scattering cross-section
c sm = xm*s,
c iq1 - projectile parton type for the scattering
c iq2 - target parton type for the scattering
c jj=1 - uncut,
c jj=2 - cut
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr51/ epsxmn
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)
      common /arr9/    x9(3),a9(3)

      qgfht=0.d0     
      zmin=4.d0*qt0*fqscal/sm/xp
      if(zmin.ge.1.d0)return
      
      qgfht=qgjetht(sm,xp*sm,xp,bbp,vvx,1.d0,1.d0,iq1,iq2,icdp,icz,jj)         
      zmax=1.d0-epsxmn
      if(zmin.ge.zmax)return

      dpx1=0.d0
      zmin1=min(.2d0,1.d0-zmin)
      do i1=1,3 !7
      do m1=1,2
       z=1.d0-epsxmn*(zmin1/epsxmn)**(.5d0+x9(i1)*(m1-1.5d0))       
       si=xp*z*sm
       fb=qgjetht(sm,si,xp,bbp,vvx,z,1.d0,iq1,iq2,icdp,icz,jj)        
     * +qgjetht(sm,si,xp,bbp,vvx,1.d0,z,iq1,iq2,icdp,icz,jj)
       dpx1=dpx1+a9(i1)*fb*(1.d0-z)
      enddo
      enddo
      qgfht=qgfht+dpx1*dlog(zmin1/epsxmn)/2.d0
      
      if(zmin.lt..8d0)then
       dpx2=0.d0
       zmin1=zmin**(-delh)
       zmax1=.8d0**(-delh)
       do i1=1,3 !7
       do m1=1,2
        z=(.5d0*(zmax1+zmin1+(zmax1-zmin1)*x9(i1)*(2*m1-3)))
     *  **(-1.d0/delh)
        si=xp*z*sm
        fb=qgjetht(sm,si,xp,bbp,vvx,z,1.d0,iq1,iq2,icdp,icz,jj)
     *  +qgjetht(sm,si,xp,bbp,vvx,1.d0,z,iq1,iq2,icdp,icz,jj)
        dpx2=dpx2+a9(i1)*fb*z**(1.d0+delh)
       enddo
       enddo
       qgfht=qgfht+dpx2*(zmin1-zmax1)/delh/2.d0
      endif
      
      zmax=(1.d0-epsxmn)**2
      if(zmin.ge.zmax)return
      
      dpx3=0.d0
      zmin1=zmin**(-delh)
      zmax1=zmax**(-delh)
      do i1=1,3 !7
      do m1=1,2
       z=(.5d0*(zmax1+zmin1+(zmax1-zmin1)*x9(i1)*(2*m1-3)))
     * **(-1.d0/delh)
       si=xp*z*sm
       
       fb1=0.d0
       zmin2=min(.2d0,1.d0-dsqrt(z))
       do i2=1,3 !7
       do m2=1,2
        z1=1.d0-epsxmn*(zmin2/epsxmn)**(.5d0+x9(i2)*(m2-1.5d0))
        z2=z/z1
        fb1=fb1+a9(i2)*(qgjetht(sm,si,xp,bbp,vvx,z1,z2,iq1,iq2,icdp
     *  ,icz,jj)+qgjetht(sm,si,xp,bbp,vvx,z2,z1,iq1,iq2,icdp,icz,jj))
     *  *(1.d0/z1-1.d0)
       enddo
       enddo
       fb1=fb1*dlog(zmin2/epsxmn)
       
       fb2=0.d0
       if(z.lt..64d0)then
        do i2=1,3 !7
        do m2=1,2
         z1=.8d0*(dsqrt(z)/.8d0)**(.5d0+x9(i2)*(m2-1.5d0))
          z2=z/z1
         fb2=fb2+a9(i2)*(qgjetht(sm,si,xp,bbp,vvx,z1,z2,iq1,iq2,icdp
     *   ,icz,jj)+qgjetht(sm,si,xp,bbp,vvx,z2,z1,iq1,iq2,icdp,icz,jj))
        enddo
        enddo
        fb2=fb2*dlog(.8d0/dsqrt(z))
       endif
       dpx3=dpx3+a9(i1)*(fb1+fb2)*z**(1.d0+delh)
      enddo
      enddo
      qgfht=qgfht+dpx3*(zmin1-zmax1)/delh/4.d0
      return 
      end

c=============================================================================
      double precision function qgjetht(sm,si,xp,bbp,vvx,zp,zm
     *,iq1,iq2,icdp,icz,jj)
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr73/ htfacm
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)
     
      qgjetht=0.d0
      tmin=qt0*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qt0*fqscal/si)))
      if(tmin.ge.si/2.d0)return
      
      xg=qt0/zm/sm
      pdfg=qgpdfbi(xg,bbp,0.d0,0.d0,icdp,icz,1,1)

      do i=1,7
      do m=1,2
       t=2.d0*tmin/(1.d0+2.d0*tmin/si
     * -x1(i)*(2*m-3)*(1.d0-2.d0*tmin/si))
       qt=t*(1.d0-t/si)
       if(qt.lt.qt0)stop'qgjetht: qt<qt0'
       facht=htfacm*pi**3*2.d0/t*qgalf(qt/fqscal/alm)*pdfg
       
       fb=0.d0
       if(zp*zm.gt..99999d0)then
        xpm=xp*(1.d0+facht*(3.d0*(2-iq1)+4.d0/3.d0*(iq1-1)))
        sim=si
        if(xpm.lt.1.d0)then
         pdfp=qgpdfbi(xpm,bbp,vvx,0.d0,icdp,icz,iq1,jj)*xp/xpm
         iq=2*iq1-2
         do n=1,3
          fb=fb+qgfbor(sim,t,iq,iq2-1,n)
          if(iq.ne.iq2-1)fb=fb+qgfbor(sim,sim-t,iq,iq2-1,n)
         enddo
         fb=fb*qgsudx(qt/fqscal,iq1)*qgsudx(qt/fqscal,iq2)
     &   /qgsudx(qt0,iq1)/qgsudx(qt0,iq2)*pdfp/sim**2
        endif

       elseif(zm.eq.1.d0)then
        do iql=1,2
         iq=2*iql-2
         dfb=0.d0
         xpm=xp*(1.d0+facht*(3.d0*(2-iql)+4.d0/3.d0*(iql-1)))
         sim=si
         if(xpm.lt.1.d0)then
          pdfp=qgpdfbi(xpm,bbp,vvx,0.d0,icdp,icz,iq1,jj)*xp/xpm
          do n=1,3
           dfb=dfb+qgfbor(sim,t,iq,iq2-1,n)
           if(iq.ne.iq2-1)dfb=dfb+qgfbor(sim,sim-t,iq,iq2-1,n)
          enddo
          fb=fb+dfb*qgevi(qt0,qt/fqscal,zp,iq1,iql)*pdfp/sim**2
         endif
        enddo
        fb=fb*qgsudx(qt/fqscal,iq2)/qgsudx(qt0,iq2)

       elseif(zp.eq.1.d0)then
        xpm=xp*(1.d0+facht*(3.d0*(2-iq1)+4.d0/3.d0*(iq1-1)))
        sim=si
        if(xpm.lt.1.d0)then
         pdfp=qgpdfbi(xpm,bbp,vvx,0.d0,icdp,icz,iq1,jj)*xp/xpm
         do iqr=1,2
          iq=2*iqr-2
          dfb=0.d0
          do n=1,3
           dfb=dfb+qgfbor(sim,t,iq1-1,iq,n)
           if(iq.ne.iq1-1)dfb=dfb+qgfbor(sim,sim-t,iq1-1,iq,n)
          enddo
          fb=fb+dfb*qgevi(qt0,qt/fqscal,zm,iq2,iqr)
         enddo
         fb=fb*qgsudx(qt/fqscal,iq1)/qgsudx(qt0,iq1)*pdfp/sim**2
        endif

       else
        do iql=1,2
         iq=2*iql-2
         xpm=xp*(1.d0+facht*(3.d0*(2-iql)+4.d0/3.d0*(iql-1)))
         sim=si
         if(xpm.lt.1.d0)then
          pdfp=qgpdfbi(xpm,bbp,vvx,0.d0,icdp,icz,iq1,jj)*xp/xpm
          do iqr=1,2
           dfb=0.d0
           do n=1,3
            dfb=dfb+qgfbor(sim,t,iq,iqr-1,n)
            if(iq.ne.iqr-1)dfb=dfb+qgfbor(sim,sim-t,iq,iqr-1,n)
           enddo
           fb=fb+dfb*qgevi(qt0,qt/fqscal,zp,iq1,iql)
     *     *qgevi(qt0,qt/fqscal,zm,iq2,iqr)*pdfp/sim**2
          enddo
         endif
        enddo
       endif

       qgjetht=qgjetht+a1(i)*fb*qgalf(qt/fqscal/alm)**2*t**2
      enddo
      enddo
      qgjetht=qgjetht*(1.d0/tmin-2.d0/si)*pi**3*2.d0
      return
      end

c-----------------------------------------------------------------------
      double precision function qgfhti(sy,xp,bb,vvx,genh
     *,iq1,iq2,icdp,icz,jj)
c-----------------------------------------------------------------------
c qgfhti - integrand for HT correction to the eikonal
c sy = xm*s,
c bb   - impact parameter squared,
c xp   - LC+ momentum fraction,
c vvx  - screening factor,
c iq1 - projectile parton type for the scattering,
c iq2 - target parton type for the scattering
c icdp - diffractive state for the hadron,
c icz  - hadron class
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wz(3),wi(3),wj(3),wm(3)
      parameter(nfock=3)
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr20/ spmax
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr74/ feikht(20,10,11,66,24*nfock),ffhtm(20,10,12,6,8)
     *,flhtm(20,12,6,2,2)
      common /qgdebug/ debug
      
      qgfhti=0.d0
      s2min=4.d0*fqscal*qt0
      if(xp*sy.le.s2min)goto 1

      yl=log(sy/s2min)/log(spmax/s2min)*20.d0
      k=max(1,int(yl))
      k=min(k,18)     
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)

      yx=dlog(xp)/log(sy/s2min)*11.d0+11.d0
      i=max(1,int(yx))
      i=min(i,8)     
      wi(2)=yx-i
      wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
      wi(1)=1.d0-wi(2)+wi(3)
      wi(2)=wi(2)-2.d0*wi(3)
      
      yv=vvx*5.d0+1.d0
      j=max(1,int(yv))
      j=min(j,4)     
      wj(2)=yv-j
      wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
      wj(1)=1.d0-wj(2)+wj(3)
      wj(2)=wj(2)-2.d0*wj(3)
       
      rp=(rq(icdp,icz)-alfp*log(xp))*4.d0*.0389d0
      z=exp(-bb/rp)
      z2=.2d0*exp(-12.5d0)
      if(z.gt.z2)then
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=(-bb/rp-dlog(0.2d0))/2.5d0+7.d0
       endif
       jz=min(9,int(zz))
       jz=max(2,jz)
       if(jz.eq.6)jz=5
       wz(2)=zz-jz
       wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
       wz(1)=1.d0-wz(2)+wz(3)
       wz(2)=wz(2)-2.d0*wz(3)
      else
       z3=z2*exp(2.5d0)
       jz=1
       wz(1)=(z-z2)*(z-z3)/z2/z3
       wz(2)=z*(z-z3)/z2/(z2-z3)
       wz(3)=z*(z-z2)/z3/(z3-z2)
      endif
      
      ii=iq1+2*(iq2-1)+4*(icdp-1)+4*nfock*(icz-1)+12*nfock*(jj-1)
      if(genh.eq.1.d0)then
       do j1=1,3
        j2=j+j1-1
       do l1=1,3
        l2=jz+l1-1
       do i1=1,3
        i2=i+i1-1
       do k1=1,3
        k2=k+k1-1
        qgfhti=qgfhti+feikht(k2,i2,l2,j2,ii)
     *  *wk(k1)*wz(l1)*wj(j1)*wi(i1)
       enddo
       enddo
       enddo
       enddo
      else
       if(genh.le.1.5d0)then
        yg=(genh-1.d0)*10.d0+1.d0
       else
        yg=dlog(genh/1.5d0)*2.d0+6.d0
       endif
       m=min(9,int(yg))
       m=max(1,m)
       if(m.eq.5)m=4
       wm(2)=yg-m
       wm(3)=wm(2)*(wm(2)-1.d0)*.5d0
       wm(1)=1.d0-wm(2)+wm(3)
       wm(2)=wm(2)-2.d0*wm(3)
       do m1=1,3
        m2=m+m1-2
       do j1=1,3
        j2=j+j1-1
       do l1=1,3
        l2=jz+l1-1
       do i1=1,3
        i2=i+i1-1
       do k1=1,3
        k2=k+k1-1
        qgfhti=qgfhti+feikht(k2,i2,l2,j2+6*m2,ii)
     *  *wk(k1)*wz(l1)*wj(j1)*wi(i1)*wm(m1)
       enddo
       enddo
       enddo
       enddo
       enddo
      endif
1     qgfhti=exp(qgfhti)
      return 
      end

c=============================================================================
      subroutine qghtm(sig3htm,sy,b,vvx,vvxp,vvxt,genh
     *,icdp,icdt,iczp,iczt)
c-----------------------------------------------------------------------
c qghtm - ht-fan ccorrections (>=2 hard legs)
c sy   - c.m. energy squared,
c b    - impact parameter,
c icdp - projectile diffractive state,
c icdt - target diffractive state,
c iczp - projectile class,
c iczt - target class
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      double precision dpeik(2,2,2),feik(2,2,2),pdft(3,2),sig3htm(5)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(nfock,3),rq(nfock,3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr38/ htfac
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)
      common /arr9/    x9(3),a9(3)
      common /arr4/    x4(2),a4(2)
      
      do i=1,5
       sig3htm(i)=0.d0
      enddo
      shmin=4.d0*fqscal*qt0  
      if(sy.le.shmin*sgap)return

      if(vvx+vvxp+vvxt.gt.0.d0.or.genh.ne.1.d0)then
       ixxm=3
      else
       ixxm=7
      endif
      do ix1=1,ixxm
      do mx1=1,2
       if(ixxm.eq.7)then
        xpomr=(sy/shmin/sgap)**(-.5d0-x1(ix1)*(mx1-1.5d0))/sgap
       else
        xpomr=(sy/shmin/sgap)**(-.5d0-x9(ix1)*(mx1-1.5d0))/sgap
       endif
       rp1=(rq(icdp,iczp)-alfp*log(xpomr))*4.d0*.0389d0
       rp2=(rq(icdt,iczt)+alfp*log(xpomr*sy))*4.d0*.0389d0
       rp=rp1*rp2/(rp1+rp2)

       do ib1=1,ixxm
       do mb1=1,2
        if(ixxm.eq.7)then
         z=.5d0+x1(ib1)*(mb1-1.5d0)
        else
         z=.5d0+x9(ib1)*(mb1-1.5d0)
        endif
        bb0=-rp*dlog(z)
       do ib2=1,ixxm
       do mb2=1,2
        if(ixxm.eq.7)then
         phi=pi*(.5d0+x1(ib2)*(mb2-1.5d0))
        else
         phi=pi*(.5d0+x9(ib2)*(mb2-1.5d0))
        endif
        bb1=(b*rp1/(rp1+rp2)+dsqrt(bb0)*cos(phi))**2+bb0*sin(phi)**2
        bb2=(b*rp2/(rp1+rp2)-dsqrt(bb0)*cos(phi))**2+bb0*sin(phi)**2

        vpf0=qgfani(1.d0/xpomr,bb1,1.d0-(1.d0-vvx)*(1.d0-vvxt)
     *  ,0.d0,0.d0,icdp,iczp,3)
        vtf0=qgfani(xpomr*sy,bb2,1.d0-(1.d0-vvx)*(1.d0-vvxp)
     *  ,0.d0,0.d0,icdt,iczt,3)
       
        n=1
1       n=n+1
        vpf=qgfani(1.d0/xpomr,bb1,1.d0-exp(-vtf0)*(1.d0-vvx)
     *  *(1.d0-vvxt),0.d0,0.d0,icdp,iczp,3)
        vtf=qgfani(xpomr*sy,bb2,1.d0-exp(-vpf0)*(1.d0-vvx)*(1.d0-vvxp)
     *  ,0.d0,0.d0,icdt,iczt,3)
        if(vpf0.gt.0.d0.and.vtf0.gt.0.d0.and.abs(1.d0-vpf/vpf0)
     *  +abs(1.d0-vtf/vtf0).gt.1.d-2.and.n.lt.100)then
         vpf0=vpf
         vtf0=vtf
         goto 1
        endif
        fanp=vpf
        fant=vtf
        fanpc=qgfani(1.d0/xpomr,bb1,1.d0-exp(-vtf0)*(1.d0-vvx)
     *  *(1.d0-vvxt),vvxp,0.d0,icdp,iczp,16)
        fantc=qgfani(xpomr*sy,bb2,1.d0-exp(-vpf0)*(1.d0-vvx)
     *  *(1.d0-vvxp),vvxt,0.d0,icdt,iczt,16)
        vvxi=1.d0-exp(-fanp-fant)*(1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt)
        dpx=(1.d0-exp(-fanp))*exp(-fant)*(1.d0-vvx)*(1.d0-vvxt)
        dpxc=fanpc*exp(-2.d0*fanp-2.d0*fant)
     *  *(1.d0-vvx)**2*(1.d0-vvxt)**2*(1.d0-vvxp)**2

        do jj=1,2
        do jv=1,2
        do ieik=1,2
         dpeik(ieik,jj,jv)=0.d0
        enddo
        enddo
        enddo
        do ix2=1,7
        do mx2=1,2
         xm=(xpomr*sy/shmin)**(-(.5+x1(ix2)*(mx2-1.5)))
         xg=qt0/xm/sy
         pdfg=qgpdfbi(xg,bb1,0.d0,0.d0,icdp,iczp,1,1)
         facht=htfac*pi**3*2.d0*pdfg*genh
       
         vpf0=qgfani(xm*sy,bb1,1.d0-(1.d0-vvx)*(1.d0-vvxt),0.d0,0.d0
     *   ,icdp,iczp,3)
         if(xm*sgap.lt.1.d0)then
          vtf0=qgfani(1.d0/xm,bb2,1.d0-(1.d0-vvx)*(1.d0-vvxp),0.d0,0.d0
     *    ,icdt,iczt,3)
         else
          vtf0=0.d0
         endif
        
         n=1
2        n=n+1
         vpf=qgfani(xm*sy,bb1,1.d0-exp(-vtf0)*(1.d0-vvx)*(1.d0-vvxt)
     *   ,0.d0,0.d0,icdp,iczp,3)
         vvxpr=1.d0-exp(-vpf0)*(1.d0-vvx)*(1.d0-vvxp)
         if(xm*sgap.lt.1.d0)then
          vtf=qgfani(1.d0/xm,bb2,vvxpr,0.d0,0.d0,icdt,iczt,3)
         else
          vtf=0.d0
         endif
         if(vpf0.gt.0.d0.and.vtf0.gt.0.d0.and.abs(1.d0-vpf/vpf0)
     *   +abs(1.d0-vtf/vtf0).gt.1.d-2.and.n.lt.100)then
          vpf0=vpf
          vtf0=vtf
          goto 2
         endif

         do jj=1,2
         do iqt=1,3
          pdft(iqt,jj)=qgpdfbi(xm,bb2,vvxpr,vvxt,icdt,iczt,iqt,jj)
         enddo
         enddo
         do iqt=1,2
         do jj=1,2
          feik(1,iqt,jj)=qglhtmi(xm*xpomr*sy,-1.d0,vvxi,iqt,jj)
          feik(2,iqt,jj)=feik(1,iqt,jj)
     *    *qglhtmi(xm*xpomr*sy,facht,vvxi,iqt,jj)
         enddo
         enddo
         do ieik=1,2
         do jj=1,2
         do jv=1,2
          do iqt=jv,2
           dpeik(ieik,jj,jv)=dpeik(ieik,jj,jv)+a1(ix2)
     *     *feik(ieik,iqt,jj)*pdft(iqt+jv-1,jj)
          enddo
         enddo
         enddo
         enddo
        enddo
        enddo
        do ieik=1,2
        do jj=1,2
        do jv=1,2
         dpeik(ieik,jj,jv)=dpeik(ieik,jj,jv)
     *   *dlog(xpomr*sy/shmin)/4.d0*factk          !1/2 included
        enddo
        enddo
        enddo

        do iqq=1,5
         if(iqq.eq.1)then               !g-uncut
          dpxi=dpx*(min(0.d0,1.d0-exp(-dpeik(2,1,1))-dpeik(2,1,1))
     *    -min(0.d0,1.d0-exp(-dpeik(1,1,1))-dpeik(1,1,1)))
         elseif(iqq.eq.2)then           !q-uncut
          dpxi=dpx*(dpeik(1,1,2)*(1.d0-exp(-dpeik(1,1,1)))
     *    -dpeik(2,1,2)*(1.d0-exp(-dpeik(2,1,1))))
         elseif(iqq.eq.3)then           !g-cut
          dpxi=dpxc*(dpeik(1,2,1)*((1.d0-exp(-2.d0*dpeik(1,1,1)))
     *    +2.d0*dpeik(1,1,2)*exp(-2.d0*dpeik(1,1,1)))
     *    -dpeik(2,2,1)*((1.d0-exp(-2.d0*dpeik(2,1,1)))
     *    +2.d0*dpeik(2,1,2)*exp(-2.d0*dpeik(2,1,1))))
         elseif(iqq.eq.4)then           !q-cut
          dpxi=dpxc*(dpeik(1,2,2)*(1.d0-exp(-2.d0*dpeik(1,1,1)))
     *    -dpeik(2,2,2)*(1.d0-exp(-2.d0*dpeik(2,1,1))))
         elseif(iqq.eq.5)then           !soft-cut
          dpxi=dpxc*fantc*((1.d0-exp(-2.d0*dpeik(1,1,1)))
     *    +2.d0*dpeik(1,1,2)*exp(-2.d0*dpeik(1,1,1))
     *    -(1.d0-exp(-2.d0*dpeik(2,1,1)))
     *    -2.d0*dpeik(2,1,2)*exp(-2.d0*dpeik(2,1,1)))
         endif

         if(ixxm.eq.7)then
          sig3htm(iqq)=sig3htm(iqq)+a1(ix1)*a1(ib1)*a1(ib2)*dpxi/z*rp
         else
          sig3htm(iqq)=sig3htm(iqq)+a9(ix1)*a9(ib1)*a9(ib2)*dpxi/z*rp
         endif
        enddo
       enddo
       enddo
       enddo
       enddo
      enddo
      enddo
      do iqq=1,5
       sig3htm(iqq)=sig3htm(iqq)*dlog(sy/shmin/sgap)
     * *pi*r3p/g3p**3/.0389d0/8.d0
      enddo
      return 
      end

c=============================================================================
      double precision function qglhtm(sm,facht,vvxi,iqt,jj)
c-----------------------------------------------------------------------------
c qglhtm - integrand for HT-fan contribution (>=2 hard legs)
c sm = xm*xpomr*s,
c vvxi - screeing correction for intermediate pomeron,
c iq2 - target parton type for the scattering
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      qglhtm=0.d0     
      xpmin=4.d0*qt0*fqscal/sm
      if(xpmin.ge.1.d0)return
      
      do ix1=1,7
      do mx1=1,2
       xp=xpmin**(.5d0+x1(ix1)*(mx1-1.5d0))
      do iqp=1,2
       dpx=qgfhtmi(xp*sm,xp,-1.d0,vvxi,iqp,iqt,jj)
       dpx=dpx*qgfhtmi(xp*sm,xp,facht,vvxi,iqp,iqt,jj)
       qglhtm=qglhtm+a1(ix1)*dpx
      enddo
      enddo
      enddo
      qglhtm=-qglhtm*dlog(xpmin)/2.d0
      return 
      end
      
c-----------------------------------------------------------------------
      double precision function qglhtmi(sm,facht,vvxi,iqt,jj)
c-----------------------------------------------------------------------
c qglhtmi - integrand for HT-fan contribution (interpolation)
c sm = xm*xpomr*s,
c vvxi  - screeing correction for intermediate pomeron,
c facht - strength of HT correction,
c iqt   - target parton type for the scattering
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wi(3),wj(3)
      parameter(nfock=3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr20/ spmax
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr74/ feikht(20,10,11,66,24*nfock),ffhtm(20,10,12,6,8)
     *,flhtm(20,12,6,2,2)
      common /qgdebug/ debug
      
      qglhtmi=0.d0
      s2min=4.d0*fqscal*qt0
      if(sm.le.s2min)return

      yl=log(sm/s2min)/log(spmax/s2min/sgap)*20.d0
      k=max(1,int(yl))
      k=min(k,18)     
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)
      
      yv=vvxi*5.d0+1.d0
      i=max(1,int(yv))
      i=min(i,4)     
      wi(2)=yv-i
      wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
      wi(1)=1.d0-wi(2)+wi(3)
      wi(2)=wi(2)-2.d0*wi(3)

      if(facht.le.0.d0)then
       if(facht.lt.0.d0)then
        jf=1
       else
        jf=2
       endif
       do i1=1,3
        i2=i+i1-1
       do k1=1,3
        k2=k+k1-1
        qglhtmi=qglhtmi+flhtm(k2,jf,i2,iqt,jj)*wk(k1)*wi(i1)
       enddo
       enddo
      else
       if(facht.le.5.d0)then
        yf=facht+2.d0
       else
        yf=dlog(facht/5.d0)*2.d0+7.d0
       endif
       j=max(2,int(yf))
       j=min(j,10)
       if(j.eq.6)j=5     
       wj(2)=yf-j
       wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
       wj(1)=1.d0-wj(2)+wj(3)
       wj(2)=wj(2)-2.d0*wj(3)

       do i1=1,3
        i2=i+i1-1
       do j1=1,3
        j2=j+j1-1
       do k1=1,3
        k2=k+k1-1
        qglhtmi=qglhtmi+flhtm(k2,j2,i2,iqt,jj)*wk(k1)*wj(j1)*wi(i1)
       enddo
       enddo
       enddo
      endif
      qglhtmi=exp(qglhtmi)
      return 
      end

c=============================================================================
      double precision function qglm0(sm,vvxi,iqt,jj)
c-----------------------------------------------------------------------------
c qglm0 - convolution of sig-jet with intermediate pomeron
c sm = xm*xpomr*s,
c vvxi - screeing correction for intermediate pomeron,
c iq2 - target parton type for the scattering
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)

      qglm0=0.d0     
      xpmin=4.d0*qt0*fqscal/sm
      if(xpmin.ge.1.d0)return
      
      do ix1=1,7
      do mx1=1,2
       xp=xpmin**(.5d0+x1(ix1)*(mx1-1.5d0))
      do iqp=1,2
       sigj=qgjit(qt0,qt0,xp*sm,iqp,iqt)*qgloopi(1.d0/xp,vvxi,iqp,jj)
       qglm0=qglm0+a1(ix1)*sigj
      enddo
      enddo
      enddo
      qglm0=-qglm0*dlog(xpmin)/2.d0
      return 
      end

c=============================================================================
      double precision function qgfhtm(sy,xp,facht,vvxi,iqp,iqt,jj)
c-----------------------------------------------------------------------------
c qgfhtm - integrand for HT-fan contribution (>=2 hard legs)
c sy = xp*xm*xpomr*s,
c vvxi - screeing correction for intermediate pomeron,
c iq1 - projectile parton type for the scattering,
c iq2 - target parton type for the scattering
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr51/ epsxmn
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)
      common /arr9/    x9(3),a9(3)

      qgfhtm=0.d0     
      s2min=4.d0*qt0*fqscal
      if(sy.le.s2min)return
      
      dxpt=qgjhtm(sy,xp,1.d0,1.d0,facht,vvxi,iqp,iqt,jj)
      zmin=s2min/sy
      zmax=1.d0-epsxmn
      if(zmin.lt.zmax)then
       dpx1=0.d0
       zmin1=min(.2d0,1.d0-zmin)
       do i1=1,3 !7
       do m1=1,2
        z=1.d0-epsxmn*(zmin1/epsxmn)**(.5d0+x9(i1)*(m1-1.5d0))
        si=z*sy
        fb=qgjhtm(si,xp,z,1.d0,facht,vvxi,iqp,iqt,jj)
     *  +qgjhtm(si,xp,1.d0,z,facht,vvxi,iqp,iqt,jj)
        dpx1=dpx1+a9(i1)*fb*(1.d0-z)
       enddo
       enddo
       dxpt=dxpt+dpx1*dlog(zmin1/epsxmn)/2.d0
      
       if(zmin.lt..8d0)then
        dpx2=0.d0
        zmin1=zmin**(-delh)
        zmax1=.8d0**(-delh)
        do i1=1,3 !7
        do m1=1,2
         z=(.5d0*(zmax1+zmin1+(zmax1-zmin1)*x9(i1)*(2*m1-3)))
     *   **(-1.d0/delh)
         si=z*sy
         fb=qgjhtm(si,xp,z,1.d0,facht,vvxi,iqp,iqt,jj)
     *   +qgjhtm(si,xp,1.d0,z,facht,vvxi,iqp,iqt,jj)
         dpx2=dpx2+a9(i1)*fb*z**(1.d0+delh)
        enddo
        enddo
        dxpt=dxpt+dpx2*(zmin1-zmax1)/delh/2.d0
       endif
      endif
      
      zmax=(1.d0-epsxmn)**2
      if(zmin.lt.zmax)then
       dpx3=0.d0
       zmin1=zmin**(-delh)
       zmax1=zmax**(-delh)
       do i1=1,3 !7
       do m1=1,2
        z=(.5d0*(zmax1+zmin1+(zmax1-zmin1)*x9(i1)*(2*m1-3)))
     *  **(-1.d0/delh)
        si=z*sy
       
        fb1=0.d0
        zmin2=min(.2d0,1.d0-dsqrt(z))
        do i2=1,3 !7
        do m2=1,2
         z1=1.d0-epsxmn*(zmin2/epsxmn)**(.5d0+x9(i2)*(m2-1.5d0))
         z2=z/z1
         fb1=fb1+a9(i2)*(qgjhtm(si,xp,z1,z2,facht,vvxi,iqp,iqt,jj)
     *   +qgjhtm(si,xp,z2,z1,facht,vvxi,iqp,iqt,jj))*(1.d0/z1-1.d0)
        enddo
        enddo
        fb1=fb1*dlog(zmin2/epsxmn)
       
        fb2=0.d0
        if(z.lt..64d0)then
         do i2=1,3 !7
         do m2=1,2
          z1=.8d0*(dsqrt(z)/.8d0)**(.5d0+x9(i2)*(m2-1.5d0))
          z2=z/z1
          fb2=fb2+a9(i2)*(qgjhtm(si,xp,z1,z2,facht,vvxi,iqp,iqt,jj)
     *    +qgjhtm(si,xp,z2,z1,facht,vvxi,iqp,iqt,jj))
         enddo
         enddo
         fb2=fb2*dlog(.8d0/dsqrt(z))
        endif
        dpx3=dpx3+a9(i1)*(fb1+fb2)*z**(1.d0+delh)
       enddo
       enddo
       dxpt=dxpt+dpx3*(zmin1-zmax1)/delh/4.d0
      endif
      qgfhtm=dxpt
      return 
      end

c=============================================================================
      double precision function qgjhtm(si,xp,zp,zm,facht,vvxi,iq1,iq2
     *,jj)
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      common /qgarr6/  pi,bm,amws
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)
     
      qgjhtm=0.d0
      tmin=qt0*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qt0*fqscal/si)))
      if(tmin.ge.si/2.d0)return
      
      do i=1,7
      do m=1,2
       t=2.d0*tmin/(1.d0+2.d0*tmin/si
     * -x1(i)*(2*m-3)*(1.d0-2.d0*tmin/si))
       qt=t*(1.d0-t/si)
      
       fb=0.d0
       if(zp*zm.gt..99999d0)then
        iq=2*iq1-2
        xpm=xp*(1.d0+facht*(3.d0*(2-iq1)+4.d0/3.d0*(iq1-1))/t
     *  *qgalf(qt/fqscal/alm))
        if(xpm.lt.1.d0)then
         do n=1,3
          fb=fb+qgfbor(si,t,iq,iq2-1,n)
          if(iq.ne.iq2-1)fb=fb+qgfbor(si,si-t,iq,iq2-1,n)
         enddo
         fb=fb*qgsudx(qt/fqscal,iq1)*qgsudx(qt/fqscal,iq2)
     &   /qgsudx(qt0,iq1)/qgsudx(qt0,iq2)/si**2
     *   *qgloopi(1.d0/xpm,vvxi,iq1,jj)*xp/xpm
        endif

       elseif(zm.eq.1.d0)then
        do iql=1,2
         iq=2*iql-2
         xpm=xp*(1.d0+facht*(3.d0*(2-iql)+4.d0/3.d0*(iql-1))/t
     *   *qgalf(qt/fqscal/alm))
         if(xpm.lt.1.d0)then
          dfb=0.d0
          do n=1,3
           dfb=dfb+qgfbor(si,t,iq,iq2-1,n)
           if(iq.ne.iq2-1)dfb=dfb+qgfbor(si,si-t,iq,iq2-1,n)
          enddo
          fb=fb+dfb*qgevi(qt0,qt/fqscal,zp,iq1,iql)/si**2
     *    *qgloopi(1.d0/xpm,vvxi,iq1,jj)*xp/xpm
         endif
        enddo
        fb=fb*qgsudx(qt/fqscal,iq2)/qgsudx(qt0,iq2)

       elseif(zp.eq.1.d0)then
        xpm=xp*(1.d0+facht*(3.d0*(2-iq1)+4.d0/3.d0*(iq1-1))/t
     *  *qgalf(qt/fqscal/alm))
        if(xpm.lt.1.d0)then
         do iqr=1,2
          iq=2*iqr-2
          dfb=0.d0
          do n=1,3
           dfb=dfb+qgfbor(si,t,iq1-1,iq,n)
           if(iq.ne.iq1-1)dfb=dfb+qgfbor(si,si-t,iq1-1,iq,n)
          enddo
          fb=fb+dfb*qgevi(qt0,qt/fqscal,zm,iq2,iqr)
         enddo
         fb=fb*qgsudx(qt/fqscal,iq1)/qgsudx(qt0,iq1)/si**2
     *   *qgloopi(1.d0/xpm,vvxi,iq1,jj)*xp/xpm
        endif

       else
        do iql=1,2
         iq=2*iql-2
         xpm=xp*(1.d0+facht*(3.d0*(2-iql)+4.d0/3.d0*(iql-1))/t
     *   *qgalf(qt/fqscal/alm))
         if(xpm.lt.1.d0)then
          do iqr=1,2
           dfb=0.d0
           do n=1,3
            dfb=dfb+qgfbor(si,t,iq,iqr-1,n)
            if(iq.ne.iqr-1)dfb=dfb+qgfbor(si,si-t,iq,iqr-1,n)
           enddo
           fb=fb+dfb*qgevi(qt0,qt/fqscal,zp,iq1,iql)
     *     *qgevi(qt0,qt/fqscal,zm,iq2,iqr)/si**2
     *     *qgloopi(1.d0/xpm,vvxi,iq1,jj)*xp/xpm
          enddo
         endif
        enddo
       endif
       qgjhtm=qgjhtm+a1(i)*fb*qgalf(qt/fqscal/alm)**2*t**2
      enddo
      enddo
      qgjhtm=qgjhtm*(1.d0/tmin-2.d0/si)*pi**3*2.d0
      return
      end
      
c-----------------------------------------------------------------------
      double precision function qgfhtmi(sy,xp,facht,vvxi,iqp,iqt,jj)
c-----------------------------------------------------------------------
c qgfhtmi - integrand for HT-fan contribution (interpolation)
c sy = xp*xm*xpomr*s,
c vvxi  - screeing correction for intermediate pomeron,
c facht - strength of HT correction,
c iqt   - target parton type for the scattering
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wi(3),wj(3),wm(3)
      parameter(nfock=3)
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr20/ spmax
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr74/ feikht(20,10,11,66,24*nfock),ffhtm(20,10,12,6,8)
     *,flhtm(20,12,6,2,2)
      common /qgdebug/ debug
      
      qgfhtmi=0.d0
      s2min=4.d0*fqscal*qt0
      if(sy.le.s2min)return

      yl=log(sy/xp/s2min)/log(spmax/s2min/sgap)*20.d0
      k=max(1,int(yl))
      k=min(k,18)     
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)
      
      yx=dlog(xp)/log(sy/xp/s2min)*11.d0+11.d0
      m=max(1,int(yx))
      m=min(m,8)     
      wm(2)=yx-m
      wm(3)=wm(2)*(wm(2)-1.d0)*.5d0
      wm(1)=1.d0-wm(2)+wm(3)
      wm(2)=wm(2)-2.d0*wm(3)
      
      yv=vvxi*5.d0+1.d0
      i=max(1,int(yv))
      i=min(i,4)     
      wi(2)=yv-i
      wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
      wi(1)=1.d0-wi(2)+wi(3)
      wi(2)=wi(2)-2.d0*wi(3)

      if(facht.le.0.d0)then
       if(facht.lt.0.d0)then
        jf=1
       else
        jf=2
       endif
       do i1=1,3
        i2=i+i1-1
       do m1=1,3
        m2=m+m1-1
       do k1=1,3
        k2=k+k1-1
        qgfhtmi=qgfhtmi+ffhtm(k2,m2,jf,i2,iqp+2*(iqt-1)+4*(jj-1))
     *  *wk(k1)*wi(i1)*wm(m1)
       enddo
       enddo
       enddo
      else
       if(facht.le.5.d0)then
        yf=facht+2.d0
       else
        yf=dlog(facht/5.d0)*2.d0+7.d0
       endif
       j=max(2,int(yf))
       j=min(j,10)
       if(j.eq.6)j=5     
       wj(2)=yf-j
       wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
       wj(1)=1.d0-wj(2)+wj(3)
       wj(2)=wj(2)-2.d0*wj(3)

       do i1=1,3
        i2=i+i1-1
       do j1=1,3
        j2=j+j1-1
       do m1=1,3
        m2=m+m1-1
       do k1=1,3
        k2=k+k1-1
        qgfhtmi=qgfhtmi+ffhtm(k2,m2,j2,i2,iqp+2*(iqt-1)+4*(jj-1))
     *  *wk(k1)*wj(j1)*wi(i1)*wm(m1)
       enddo
       enddo
       enddo
       enddo
      endif
      qgfhtmi=exp(qgfhtmi)
      return 
      end

c=============================================================================
      subroutine qgfanht(fanht,sy,sm,bbp,bbt,vvx,vvxt,genh
     *,icdp,icdt,iczp,iczt)
c-----------------------------------------------------------------------
c qgfanht - ht-fan correction to energy sharing
c sy   - c.m. energy squared,
c sm = xpomr*s,
c bbp  - impact parameter squared to the projectile,
c bbt  - impact parameter squared to the target,
c genh - nuclear enhancement of the soft gluon density,
c icdp - projectile diffractive state,
c icdt - target diffractive state,
c iczp - projectile class,
c iczt - target class
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      double precision fanht(2,2),pdft(3,2),feik(2,2)
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr38/ htfac
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr3/    x1(7),a1(7)
      common /arr9/    x9(3),a9(3)
      common /arr4/    x4(2),a4(2)
      
      do i=1,2
      do jv=1,2
       fanht(jv,i)=0.d0
      enddo
      enddo
      shmin=4.d0*fqscal*qt0  
      if(sm.le.shmin)return

c      if(vvxt.gt.0.d0.or.genh.ne.1.d0)then
       ixxm=2
c      else
c       ixxm=7
c      endif
      vtf=qgfani(sm,bbt,vvx,0.d0,0.d0,icdt,iczt,3)
      vvxi=1.d0-exp(-vtf)*(1.d0-vvx)*(1.d0-vvxt)
      do ix1=1,ixxm
      do mx1=1,2
       if(ixxm.eq.7)then
        xm=(sm/shmin)**(-.5d0-x1(ix1)*(mx1-1.5d0))
       else
        xm=(sm/shmin)**(-.5d0-x4(ix1)*(mx1-1.5d0))
       endif
       xg=qt0/xm/sy
       pdfg=qgpdfbi(xg,bbp,0.d0,0.d0,icdp,iczp,1,1)
       facht=htfac*pi**3*2.d0*pdfg*genh
       si=sm*xm
       
       do jj=1,2
        do iqt=1,3
         pdft(iqt,jj)=qgpdfbi(xm,bbt,vvx,vvxt,icdt,iczt,iqt,jj)
        enddo
       enddo
       do jj=1,2
        do iqt=1,2
         feik(iqt,jj)=qglhtmi(si,-1.d0,vvxi,iqt,jj)
     *   *max(0.d0,1.d0-qglhtmi(si,facht,vvxi,iqt,jj))
        enddo
       enddo
       do jj=1,2
        do jv=1,2
         do iqt=jv,2
          if(ixxm.eq.7)then
           fanht(jv,jj)=fanht(jv,jj)+a1(ix1)*feik(iqt,jj)
     *     *pdft(iqt+jv-1,jj)
          else
           fanht(jv,jj)=fanht(jv,jj)+a4(ix1)*feik(iqt,jj)
     *     *pdft(iqt+jv-1,jj)
          endif
         enddo
        enddo
       enddo
      enddo
      enddo
      do jv=1,2
      do jj=1,2
       fanht(jv,jj)=-fanht(jv,jj)*dlog(sm/shmin)*factk/4.d0  !1/2 included
      enddo
      enddo
      return 
      end

c=============================================================================
      subroutine qgfanhtm(fanhtm,sy,xpomr,bbp,bbt,vvx,vvxt,genhp
     *,icdp,icdt,iczp,iczt)
c-----------------------------------------------------------------------
c qgfanht - ht-fan correction to energy sharing
c sy = xpomr*s,
c bbp  - impact parameter squared to the projectile,
c bbt  - impact parameter squared to the target,
c genh - nuclear enhancement of the soft gluon density,
c icdp - projectile diffractive state,
c icdt - target diffractive state,
c iczp - projectile class,
c iczt - target class,
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nfock=3)
      double precision fanhtm(2,2),fanhtt(2,2,2),feik(2,2,2),vint(2)
     *,dpx(2,2),pdft(3,2)
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,rr,r3p,g3p,delh,sgap
      common /qgarr18/ alm,qt0,qtf,bbbpom,dgqq,beth(3),bbbi(nfock,3)
      common /qgarr26/ factk,fqscal
      common /qgarr38/ htfac
      common /qgarr43/ moniou
      common /qgdebug/ debug
      common /arr4/    x4(2),a4(2)
      
      do jv=1,2
      do jj=1,2
       fanhtm(jv,jj)=0.d0
      enddo
      enddo
      return !subleading (insignificant) HT corrections neglected to speed up
      shmin=4.d0*fqscal*qt0  
      if(xpomr*sy.le.shmin*sgap)return

      do ix1=1,2
      do mx1=1,2
       xpomrm=xpomr/sgap*(shmin*sgap/sy/xpomr)
     * **(.5d0+x4(ix1)*(mx1-1.5d0))

       vtf=qgfani(xpomrm*sy,bbt,vvx,0.d0,0.d0,icdt,iczt,3)
       vvxi=1.d0-exp(-vtf)*(1.d0-vvx)*(1.d0-vvxt)
       rh=4.d0*.0389d0*alfp*dlog(xpomr/xpomrm)
       do jj=1,2
        vint(jj)=0.d0
       enddo
       do ib=1,2
       do mb=1,2
        zz=.5d0+x4(ib)*(mb-1.5d0)
        bbi=-rh*dlog(zz)
        vint(1)=vint(1)+a4(ib)/zz*qgpini(xpomr/xpomrm,bbi,0.d0,0.d0,5)
        vint(2)=vint(2)+a4(ib)/zz*qgpini(xpomr/xpomrm,bbi,0.d0,0.d0,11)
       enddo
       enddo
       do jj=1,2
        vint(jj)=vint(jj)*pi*rh/.0389d0/2.d0
       enddo

       do jv=1,2
       do jj=1,2
       do ieik=1,2
        fanhtt(ieik,jv,jj)=0.d0
       enddo
       enddo
       enddo
       do ix2=1,2
       do mx2=1,2
        xm=(shmin/xpomrm/sy)**(.5d0+x4(ix2)*(mx2-1.5d0))
        xg=qt0/xm/sy
        pdfg=qgpdfbi(xg,bbp,0.d0,0.d0,icdp,iczp,1,1)*genhp
        facht=htfac*pi**3*2.d0*pdfg
        si=sy*xm*xpomrm
       
        do jj=1,2
         do iqt=1,3
          pdft(iqt,jj)=qgpdfbi(xm,bbt,vvx,vvxt,icdt,iczt,iqt,jj)
         enddo
        enddo
        do jj=1,2
         do iqt=1,2
          feik(1,iqt,jj)=qglhtmi(si,-1.d0,vvxi,iqt,jj)
          feik(2,iqt,jj)=feik(1,iqt,jj)*qglhtmi(si,facht,vvxi,iqt,jj)
         enddo
        enddo
        do jj=1,2
         do jv=1,2
          do ieik=1,2
           do iqt=jv,2
            fanhtt(ieik,jv,jj)=fanhtt(ieik,jv,jj)
     *      +a4(ix2)*feik(ieik,iqt,jj)*pdft(iqt+jv-1,jj)
           enddo
          enddo
         enddo
        enddo
       enddo
       enddo
       do jj=1,2
       do jv=1,2
       do ieik=1,2
        fanhtt(ieik,jv,jj)=-fanhtt(ieik,jv,jj)*dlog(shmin/xpomrm/sy)
     *  *factk/4.d0                  !1/2 included
       enddo
       enddo
       enddo

       dpx(1,1)=(min(0.d0,1.d0-exp(-fanhtt(2,1,1))-fanhtt(2,1,1))
     * -min(0.d0,1.d0-exp(-fanhtt(1,1,1))-fanhtt(1,1,1)))*exp(-vtf)
     * *(1.d0-vvx)
       dpx(2,1)=(fanhtt(2,2,1)*(exp(-fanhtt(2,1,1))-1.d0)
     * -fanhtt(1,2,1)*(exp(-fanhtt(1,1,1))-1.d0))*exp(-vtf)*(1.d0-vvx)
       do jv=1,2
        dpx(jv,2)=(fanhtt(1,jv,2)*(1.d0-exp(-2.d0*fanhtt(1,1,1)))
     *  -fanhtt(2,jv,2)*(1.d0-exp(-2.d0*fanhtt(2,1,1))))*(1.d0-vvxi)**2
       enddo
       do jv=1,2
       do jj=1,2
        fanhtm(jv,jj)=fanhtm(jv,jj)+a4(ix1)*dpx(jv,jj)*vint(jj)
       enddo
       enddo
      enddo
      enddo
      do jv=1,2
      do jj=1,2
       fanhtm(jv,jj)=-fanhtm(jv,jj)*dlog(shmin*sgap/sy/xpomr)
     * /2.d0*r3p/g3p**3
      enddo
      enddo
      return 
      end

