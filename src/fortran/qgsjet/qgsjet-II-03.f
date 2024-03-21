C======================================================================C
C                                                                      C
C     QQQ        GGG      SSSS    JJJJJJJ   EEEEEEE   TTTTTTT     I I  C
C    Q   Q      G   G    S    S         J   E            T        I I  C
C   Q     Q    G         S              J   E            T        I I  C
C   Q     Q    G   GGG    SSSS          J   EEEEE        T    ==  I I  C
C   Q   Q Q    G     G        S         J   E            T        I I  C
C    Q   Q      G   G    S    S    J   J    E            T        I I  C
C     QQQQQ      GGG      SSSS      JJJ     EEEEEEE      T        I I  C
C                                                                      C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C                  QUARK - GLUON - STRING - JET - II MODEL             C
C                                                                      C
C                HIGH ENERGY HADRON INTERACTION PROGRAM                C
C                                                                      C
C                                  BY                                  C
C                                                                      C
C                           S. OSTAPCHENKO                             C
C                                                                      C
C          Institut fuer Kernphysik, Forschungszentrum Karlsruhe       C
C               Moscow State University,  Moscow, Russia               C
C                   e-mail: serguei@ik.fzk.de                          C
C                           serg@eas10.eas.npi.msu.ru                  C
C----------------------------------------------------------------------C
C               Publication to be cited when using this program:       C
C         S. Ostapchenko, Phys.Lett. B636 (2006) 40                    C
C         S. Ostapchenko, Phys.Rev. D 74 (2006) 014026                 C
C----------------------------------------------------------------------C
C                        LIST OF MODIFICATIONS                         C
C                                                                      C
C (Any modification of this program has to be approved by the author)  C
C                                                                      C
C 24.01.2005 - beta-version completed (qgsjet-II-01)                   C
C 21.02.2005, 25.02.2005, 03.03.2005 - bugs removed                    C 
c (after debugging by D.Heck)                                          C
C 12.04.2005 - final version                                           C
C 12.12.2005 - technical update -  version II-03:                      C
C    improved treatment of Pomeron cuts (all "net" cuts included);     C
C    improved treatment of nuclear config. (more consistent diffr.);   C
C    "baryon junction" mechanism included (motivated by RHIC data);    C
C    better parameter calibration, e.g. including RHIC data            C
C 16.12.2005 - bug removed                                             C
C 21.02.2006 - some commons enlarged to avoid frequent rejects  D.H.   C
C 26.04.2006 - reduce unnecessary looping in qgsha              D.H.   C
C                                                                      C
C                 last modification:  21.02.2006                       C
C                        Version qgsjet-II-03                          C
C=======================================================================
c=============================================================================
      subroutine qgset
c-----------------------------------------------------------------------------
c model parameters&costants
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      character*7 ty
      character*2 tyq
      parameter(iapmax=207)
      common /qgarr3/  rmin,emax,eev
      common /qgarr6/  pi,bm,amws
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),b
      common /qgarr8/  wwm,be(4),dc(5),deta,almpt,ptdif,ptndi
      common /qgarr10/ am0,amn,amk,amc,amlamc,amlam,ameta,ammu
      common /qgarr11/ b10
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr19/ ahl(3)
      common /qgarr20/ spmax
      common /qgarr21/ dmmin(3),wex(3)
      common /qgarr26/ factk,fqscal
      common /qgarr28/ arr(5)
      common /qgarr32/ epsxmn
      common /qgarr41/ ty(6)
      common /qgarr42/ tyq(16)
      common /qgarr43/ moniou
      common /qgarr44/ jopt
      common /debug/   debug

      moniou=6             !output channel for debugging
      debug=0              !debugging level 
                           !(0 - no debugging, 1 - very geheral,
                           !2 - more detailed, 3 - all function calls,
			   !4 - all returned values, 5 - technical)
      if(debug.ge.1)write (moniou,201)
      
      spmax=1.d10          !maximal energy for tabulations
      jopt=1               !parameter option
      if(jopt.eq.1)then
c soft Pomeron parameters
       dels=.175d0         !overcriticality
       alfp=.105d0         !trajectory slope
       rq(1)=1.1d0         !vertex slope (pi)
       cd(1,1)=1.65d0      !relative strenth of the 1st diffr. eigenst. (pi)
       cd(2,1)=.35d0       !relative strenth of the 2nd diffr. eigenst. (pi)
       rq(2)=2.d0          !vertex slope (p)
       cd(1,2)=1.6d0       !relative strenth of the 1st diffr. eigenst. (p)
       cd(2,2)=.4d0        !relative strenth of the 2nd diffr. eigenst. (p)
       rq(3)=.65d0         !vertex slope (kaon)
       cd(1,3)=1.7d0       !relative strenth of the 1st diffr. eigenst. (kaon)
       cd(2,3)=.3d0        !relative strenth of the 2nd diffr. eigenst. (kaon)
       sigs=1.3d0          !soft parton cross section
       alpd=0.d0           !Pomeron end exponent
       r3p=.0135d0         !triple-Pomeron vertex
       g3p=.25d0           !slope of multi-Pomeron vertex         
       almpt=1.6d0         !string fragmentation parameter
c hard Pomeron parameters
       qt0=2.5d0           !q**2 cutoff
       factk=1.5d0         !k-factor value
       betp=2.d0           !gluon structure function parameter for soft P
       fqscal=4.d0         !factor for Q^2-scale (Q^2=p_t^2/fqscal)
       dgqq=.19d0          !sea quark contribution
       alm=.04d0           !lambda_QCD squared
       qtf=.1d0            !p_t**2 resolution scale
       delh=0.25d0         !effective sampling exponent for the hard Pomeron
       epsxmn=.01d0        !x-cutoff for soft gluon emission (only in pdf-s)
c weigts of diffractive eigenstates      
       cc(1,1)=.5d0        !pion          
       cc(2,1)=.5d0
       cc(1,2)=.5d0        !proton
       cc(2,2)=.5d0
       cc(1,3)=.5d0        !kaon
       cc(2,3)=.5d0
c leading state exponents
       ahl(1)=-0.05d0      !pion
       ahl(2)=1.3d0        !proton
       ahl(3)=-0.3d0       !kaon
c remnant excitation probabilities
       wex(1)=1.d0         !pion
       wex(2)=.7d0         !proton
       wex(3)=1.d0         !kaon                                                                    
                                   
      else
       stop'wrong option!!!'
      endif
           
c Regge intercepts for secondary trajectories
      arr(1)=0.5d0         !qq~-trajectory
      arr(2)=-0.5d0        !qqq~q~-trajectory
      arr(3)=0.d0          !us~-trajectory
      arr(4)=-2.d0         !uc~-trajectory
c string fragmentation - relative probabilities for qq~(qqq~q~)-pairs 
c creation from vacuum:
c (1-3) - for quark (u,d,u~,d~) fragmentation;
c (4-5) - for diquark (ud, u~d~) fragmentation
      dc(1)=.1d0           !udu~d~
      dc(2)=.07d0          !ss~
      dc(3)=.000d0         !charmed particles switched off
      dc(4)=.2d0           !ss~ (for di-quark fragmentation)
      dc(5)=.0d0           !charmed particles switched off
      deta=.11111d0        !ratio of etas to all pions (1/9)
      wwm=1.2d0            !mass threshold for 3-particle string decay  
c be(i) - parameters for pt-distributions 
      be(1)=.3d0           !uu~(dd~)
      be(2)=.3d0           !ss~
      be(3)=.3d0           !qqq~q~
      be(4)=.40d0          !cc~
      ptdif=.12d0          !diffractive momentum transfer
      ptndi=.22d0          !non-diffractive momentum transfer
c parameters for nuclear spectator part fragmentation:
      rmin=3.35d0          !coupling radius squared (fm^2)
      emax=.11d0           !relative critical energy ( /<E_ex>, <E_ex>~12.5 MeV )
      eev=.25d0            !elative evaporation energy ( /<E_ex>, <E_ex>~12.5 MeV )
c auxiliary constants
      b10=.43876194d0      !initial value for random generator (not used)
      pi=3.1416d0 
      amws=.523d0          !diffusion radius for saxon-wood density (A>56)
c particle masses
      amn=.939d0
      amk=.496d0
      am0=.14d0
      amc=1.868d0
      amlam=1.116d0
      amlamc=2.27d0
      ameta=.548d0
      ammu=.1057d0
c minimal low-mass diffractive mass
      dmmin(1)=.76d0       !pion
      dmmin(2)=1.23d0      !proton
      dmmin(3)=.89d0       !kaon
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
      
201   format(2x,'qgset - setting model parameters')
202   format(2x,'qgset - end')
      return
      end

c=============================================================================
      subroutine qgaini( DATDIR )
c-----------------------------------------------------------------------------
c qgaini - common initialization procedure
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      CHARACTER DATDIR*(256)
      real gamfun
      integer debug
      character *7 ty
      logical lcalc
      parameter(iapmax=207)
      dimension mij(40,40,4),nij(40,40,4,2),cs1(40,40,160,2)
     *,evs(40,100,3,2),ixemax(40,3,2),qfan0(51,11,6),fann(6)
     *,gz0(5),gz1(3)
cdh  *,gz0(5),gz1(3),wa(3)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr10/ am(7),ammu
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr19/ ahl(3)
      common /qgarr20/ spmax
      common /qgarr24/ qpomu(11,11,2,2,3),qpom1(11,11,6,4,3)
     *,qpomr(11,11,6,216,12)
      common /qgarr25/ ahv(3)
      common /qgarr26/ factk,fqscal
      common /qgarr27/ qlegi(51,11,11,3),qfanu(51,11,11,6)
     *,qfan(51,11,11,6,120)
      common /qgarr28/ arr(5)
      common /qgarr29/ cstot(40,40,160)
      common /qgarr30/ cs0(40,40,160,2)
      common /qgarr31/ csborn(40,160)
      common /qgarr32/ epsxmn
      common /qgarr33/ fsud(10,2)
      common /qgarr34/ qrt(10,101,2)
      common /qgarr35/ qlegc(11,10,11,11,24)
      common /qgarr38/ qpomc(11,100,11,11,60)
      common /qgarr39/ qpomi(51,11,11),qpomi0(11,11)
      common /qgarr41/ ty(6)
      common /qgarr43/ moniou
      common /qgarr44/ jopt
      common /qgarr45/ qppdi(10,11,2)
      common /qgarr47/ gsect(10,5,6)
      common /qgarr48/ asect(10,6,6)
      common /qgarr49/ trnuc(56),twsnuc(56),twbnuc(56)
      common /qgarr50/ x1(7),a1(7)
      common /qgarr54/ evk(40,40,100,3,2)
      common /debug/   debug

      if(debug.ge.1)write (moniou,201)jopt

c-------------------------------------------------
      write(*,100)
 100  format(' ',
     *           '====================================================',
     *     /,' ','|                                                  |',
     *     /,' ','|         QUARK GLUON STRING JET -II MODEL         |',
     *     /,' ','|                                                  |',
     *     /,' ','|         HADRONIC INTERACTION MONTE CARLO         |',
     *     /,' ','|                        BY                        |',
     *     /,' ','|                 S. OSTAPCHENKO                   |',
     *     /,' ','|                                                  |',
     *     /,' ','|            e-mail: serguei@ik.fzk.de             |',
     *     /,' ','|                                                  |',
     *     /,' ','|                   Version II-03                  |',
     *     /,' ','|                                                  |',
     *     /,' ','| Publication to be cited when using this program: |',
     *     /,' ','| S.Ostapchenko, Phys.Lett. B636 (2006) 40         |',
     *     /,' ','| S. Ostapchenko, Phys.Rev. D 74 (2006) 014026     |',
     *     /,' ','|                                                  |',
     *     /,' ','| last modification:  26.04.2006                   |',
     *     /,' ','|                                                  |',
     *     /,' ','| Any modification has to be approved by the author|',
     *     /,' ','====================================================',
     *     /)

c-----------------------------------------------------------------------------
c normalization of parton density in the soft Pomeron
      rr=gamfun(real(2.d0+betp-dels))/gamfun(real(1.d0-dels))
     */gamfun(real(1.d0+betp))/4.d0/pi
c valence quark x->1 behavior (GRV)
      ahv(1)=.383d0+.624d0*dlog(dlog(qt0/.204d0**2)
     */dlog(.26d0/.204d0**2))
      ahv(3)=ahv(1)
      sq=dlog(dlog(qt0/.232d0**2)/dlog(.23d0/.232d0**2))
      ahv(2)=2.997d0+.753d0*sq-.076d0*sq*sq

c valence quark momentum share (pion, kaon)      
      qnorm1=0.d0
      do i=1,7
      do m=1,2
       tp=1.d0-(.5d0+x1(i)*(m-1.5d0))**(2.d0/3.d0)
       xp=1.d0-tp**(1.d0/(1.d0+ahv(1)))       
       qnorm1=qnorm1+a1(i)*(qggrv(xp,qt0,1,1)+qggrv(xp,qt0,1,2))
     * /dsqrt(1.d0-tp)
      enddo
      enddo
      qnorm1=qnorm1/(1.d0+ahv(1))/3.d0
c valence quark momentum share (proton)      
      qnorm2=0.d0
      do i=1,7
      do m=1,2
       tp=1.d0-(.5d0+x1(i)*(m-1.5d0))**(2.d0/3.d0)
       xp=1.d0-tp**(1.d0/(1.d0+ahv(2)))       
       qnorm2=qnorm2+a1(i)*(qggrv(xp,qt0,2,1)+qggrv(xp,qt0,2,2))
     * /dsqrt(1.d0-tp)
      enddo
      enddo
      qnorm2=qnorm2/(1.d0+ahv(2))/3.d0     

c Pomeron-hadron vertex
      fp(2)=(1.d0-qnorm2)*gamfun(real(3.d0-alpd+ahl(2)))
     */gamfun(real(1.d0+ahl(2)))/gamfun(real(2.d0-alpd))
      do icz=1,3,2 
       fp(icz)=gamfun(real(3.d0-alpd+ahl(icz)))/gamfun(real(2.d0-alpd))
     * /gamfun(real(1.d0+ahl(icz)))*(1.d0-qnorm1)
      enddo
c check normalization
      do icz=1,3
       gnorm=0.d0
       seanrm=0.d0
       do i=1,7
       do m=1,2
        xxg=(.5d0+x1(i)*(m-1.5d0))**(1.d0/(1.d0-dels))
        gnorm=gnorm+a1(i)*qgftld(xxg,icz)
        seanrm=seanrm+a1(i)*qgftle(xxg,icz)
       enddo
       enddo
       gnorm=gnorm/(1.d0-dels)*fp(icz)*rr*2.d0*pi
       seanrm=seanrm/(1.d0-dels)*fp(icz)*rr*2.d0*pi
       if(icz.eq.2)then
        xxnorm=qnorm2+gnorm+seanrm
       else
        xxnorm=qnorm1+gnorm+seanrm
       endif      
       if(debug.ge.1)write (moniou,202)ty(icz),xxnorm,gnorm,seanrm
     * ,fp(icz),rr
      enddo
      
      inquire(file=DATDIR(1:INDEX(DATDIR,' ')-1)//'qgsdat-II-03',
     *        exist=lcalc)
      if(lcalc)then
       if(debug.ge.1)write (moniou,203)
       open(1,file=DATDIR(1:INDEX(DATDIR,' ')-1)//'qgsdat-II-03',
     *   form='unformatted',status='old')
       read (1) csborn,cs0,cstot,evk,qppdi,qlegi,qfanu,qfan
     * ,qpomu,qpom1,qpomr,qlegc,qpomc,qpomi,qpomi0,gsect,fsud,qrt
       close(1)	      
       goto 5
      endif
      
c--------------------------------------------------
c qcd evolution and qcd ladder cross sections
c--------------------------------------------------
      if(debug.ge.1)write (moniou,204)   

c qcd evolution
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
1     n=n+1
      if(debug.ge.1)write (moniou,232)n

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
      if(jec.ne.0)goto 1
      
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
      
c qcd ladder cross sections
      do i=1,40
       qi=(spmax/4.d0/fqscal)**((i-1)/39.d0)  !q^2 cutoff for born process
       s2min=qi*4.d0*fqscal          !energy threshold for 2->2 subprocess
      do m=1,2                                !parton types (1-g, 2-q)
      do l=1,2                                !parton types (1-g, 2-q) 
       l1=2*l-1
      do k=1,40
       sk=s2min*(spmax/s2min)**((k-1)/39.d0)  !c.m. energy squared
       k1=k+40*(m-1)+80*(l-1)
       csborn(i,k1)=dlog(qgborn(qi,qi,sk,m-1,l1-1,2)) !born cross-section (2->2)	  
      enddo
      enddo
      enddo
      enddo

c initialization of ladder cross-sections
      do i=1,40
       qi=(spmax/4.d0/fqscal)**((i-1)/39.d0)
      do j=1,40
       qj=qi*(spmax/4.d0/fqscal/qi)**((j-1)/39.d0)
       s2min=qj*4.d0*fqscal
      do m=1,2
      do l=1,2  
       l1=2*l-1
       ml=m+2*l-2
      do k=1,40
       sk=s2min*(spmax/s2min)**((k-1)/39.d0)
       k1=k+40*(m-1)+80*(l-1)
       cstot(i,j,k1)=dlog(qgborn(qi,qj,sk,m-1,l1-1,1))
       cs0(i,j,k1,1)=cstot(i,j,k1)
       cs0(i,j,k1,2)=dlog(qgborn(qi,qj,sk,m-1,l1-1,3))
       mij(i,j,ml)=2
       nij(i,j,ml,1)=2	  
       nij(i,j,ml,2)=2	  
      enddo
      enddo
      enddo
      enddo
      enddo

      n=2                             !number of ladder rungs considered
2     continue
      if(debug.ge.1)write (moniou,205)n,mij(1,1,1)
     *,nij(1,1,1,1),nij(1,1,1,2)

      do i=1,39
       qi=(spmax/4.d0/fqscal)**((i-1)/39.d0)       !q^2 for upper parton
      do j=1,39
       qj=qi*(spmax/4.d0/fqscal/qi)**((j-1)/39.d0) !q^2 for downer parton
       s2min=qj*4.d0*fqscal                !energy threshold for 2->2 subprocess
       smin=s2min/(1.d0-epsxmn)            !energy threshold for 2->3 subprocess
      do m=1,2                                     !parton types (1-g, 2-q)
      do l=1,2                                     !parton types (1-g, 2-q)
       l1=2*l-1
       ml=m+2*l-2
       do jj=1,2
        kmin=nij(i,j,ml,jj)                  !lowest energy bin for another rung
        if(kmin.le.40)then
         do k=kmin,40
          sk=s2min*(spmax/s2min)**((k-1)/39.d0)
          if(sk.le.smin)then
           nij(i,j,ml,jj)=nij(i,j,ml,jj)+1         !skip the bin
          else
           k1=k+40*(m-1)+80*(l-1)
           tmin=qj*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qj*fqscal/sk)))
                              !cross-section for one-way ordered ladder
           cs1(i,j,k1,jj)=dlog(qgjet1(qi,qj,sk,s2min,m,l,jj)
     *     /(1.d0/tmin-2.d0/sk)+qgborn(qi,qj,sk,m-1,l1-1,2*jj-1))
          endif
         enddo
        endif	  
       enddo
      enddo
      enddo
      enddo
      enddo

c check convergence
      do i=1,39
      do j=1,39
      do m=1,2
      do l=1,2   
       ml=m+2*l-2
       do jj=1,2
        kmin=nij(i,j,ml,jj)
        if(kmin.le.40)then
         do k=40,kmin,-1
          k1=k+40*(m-1)+80*(l-1)
          if(abs(cs1(i,j,k1,jj)-cs0(i,j,k1,jj)).gt.1.d-2)then
	   cs0(i,j,k1,jj)=cs1(i,j,k1,jj)
          elseif(k.eq.nij(i,j,ml,jj))then
           nij(i,j,ml,jj)=nij(i,j,ml,jj)+1
	  endif
         enddo
        endif
       enddo
      enddo
      enddo
      enddo
      enddo

      do i=1,39
       qi=(spmax/4.d0/fqscal)**((i-1)/39.d0)
      do j=1,39
       qj=qi*(spmax/4.d0/fqscal/qi)**((j-1)/39.d0)
       s2min=qj*4.d0*fqscal         !min energy squared for 2->2 subprocess
       smin=s2min/(1.d0-epsxmn)     !min energy squared for 2->3 subprocess
      do m=1,2
      do l=1,2
       ml=m+2*l-2
       kmin=mij(i,j,ml)             !min energy bin for more ladder rungs
       if(kmin.le.40)then
        do k=kmin,40
         sk=s2min*(spmax/s2min)**((k-1)/39.d0)
         if(sk.le.smin)then
          mij(i,j,ml)=mij(i,j,ml)+1
         else
          k1=k+40*(m-1)+80*(l-1)
          tmin=qj*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qj*fqscal/sk)))
                                    !cross-section for general ladder
          cs1(i,j,k1,1)=dlog((qgjet(qi,qj,sk,s2min,m,l)
     *    +qgjit1(qj,qi,sk,l,m,1))/(1.d0/tmin-2.d0/sk))
         endif
        enddo
       endif
      enddo
      enddo
      enddo
      enddo

c check convergence
      do i=1,39
      do j=1,39
      do m=1,2
      do l=1,2
       ml=m+2*l-2
       kmin=mij(i,j,ml)             !min energy bin for more ladder rungs
       if(kmin.le.40)then
        do k=40,kmin,-1
         k1=k+40*(m-1)+80*(l-1)
         if(abs(cs1(i,j,k1,1)-cstot(i,j,k1)).gt.1.d-2)then
          cstot(i,j,k1)=cs1(i,j,k1,1)
         elseif(k.eq.mij(i,j,ml))then     
          mij(i,j,ml)=mij(i,j,ml)+1
	 endif
        enddo
       endif
      enddo
      enddo
      enddo
      enddo

      n=n+1                         !one more rung
      do i=1,39
      do j=1,39
      do l=1,4
       if(mij(i,j,l).le.40.or.nij(i,j,l,1).le.40.or.nij(i,j,l,2).le.40)
     * goto 2
      enddo
      enddo
      enddo

c parton distributions in the Pomeron
      if(debug.ge.1)write (moniou,209)
      do iqq=1,2                               !parton type (1 - g, 2 - q_s)
      do ix=1,10                               !x-binning
       if(ix.le.5)then
        xp=.2d0*(.2d0*spmax)**((ix-6)/5.d0)
       else
        xp=.2d0*(ix-5)
       endif
      do iv=1,11
       vvx=(iv-1)/10.d0     !relative strenth of screening corrections (0<vvx<1),
       
       if(ix.eq.10.or.iv.eq.1)then
        qppdi(ix,iv,iqq)=0.d0
       else
        if(iqq.eq.1)then
         qppdi(ix,iv,1)=dlog(qgppd(xp,vvx,0)/(1.d0-xp)**betp
     *   /(1.d0-dgqq))       
        elseif(iqq.eq.2)then
         qppdi(ix,iv,2)=dlog(qgppd(xp,vvx,1)/qgftlf(xp)/dgqq)
        endif
       endif
      enddo
      enddo
      enddo

c eikonal for an itermediate gg-Pomeron
      if(debug.ge.1)write (moniou,210)
      s2min=4.d0*fqscal*qt0
      do iy=1,51
       sy=spmax**((iy-1)/50.d0)                   !Pomeron mass squared
       rp=alfp*dlog(sy)*4.d0*.0389d0
       if(iy.eq.1)then
        do iz=1,11
        do iv=1,11
         qpomi(iy,iz,iv)=0.d0
        enddo 
        enddo 
       else 
        do iz=1,11  
         if(iz.gt.6)then
          z=.2d0*(iz-6)
          b=sqrt(-log(z)*rp)                      !impact parameter between ends
         elseif(iz.gt.1)then
          b=sqrt(-rp*(log(0.2d0)+.8d0*(iz-7)))
          z=exp(-b*b/rp)
         endif      
        do iv=1,11
         vvx=(iv-1)/10.d0                         !relative screening strenth
       
         alf=4.d0*pi*r3p/g3p*sigs*vvx
         if(iz.eq.1.or.sy.le.s2min)then
          qpomi(iy,iz,iv)=-alf*dlog(sy)
         else
          qpomi(iy,iz,iv)=log(qgpint(sy,b*b,vvx)
     *    /sy**dels/z*rp/sigs/4.d0/.0389d0+sy**(-alf))
         endif
        enddo
        enddo
       endif
      enddo
      
c integrated eikonal for an itermediate gg-Pomeron
      if(debug.ge.1)write (moniou,211)
      s2min=4.d0*fqscal*qt0
      do iy=1,11
       sy=2.d0*s2min*(spmax/s2min/2.d0)**((iy-1)/10.d0)    !Pomeron mass squared
      do iv=1,11
       vvx=(iv-1)/10.d0       
       qpomi0(iy,iv)=log(qgpint0(sy,vvx)/sy**delh)
      enddo
      enddo

c integrated Pomeron leg eikonals
      if(debug.ge.1)write (moniou,212)
      do icz=1,3                                  !hadron class
      do iy=1,51
       sy=spmax**((iy-1)/50.d0)                   !Pomeron mass squared
       rp=(rq(icz)+alfp*log(sy))*4.d0*.0389d0
      do iz=1,11  
       if(iz.gt.6)then
        z=.2d0*(iz-6)
        b=sqrt(-log(z)*rp)                        !impact parameter between ends
       elseif(iz.gt.1)then
        b=sqrt(-rp*(log(0.2d0)+.8d0*(iz-7)))
        z=exp(-b*b/rp)
       else
        b=1.d5
       endif
      do iv=1,11
       vvx=(iv-1)/10.d0                           !relative screening strenth
       
       qxl=qgleg(sy,b*b,vvx,1,icz)
       if(iz.eq.1)then
        qlegi(iy,iz,iv,icz)=dlog(qxl)
       else
        if(qxl.gt.0.d0)then
         qlegi(iy,iz,iv,icz)=log(qxl/z)
        else
         qlegi(iy,iz,iv,icz)=qlegi(iy,iz,iv-1,icz)
        endif
       endif
      enddo
      enddo
      enddo
      enddo
      
c integrated fan contributions
      if(debug.ge.1)write (moniou,213)
      do icz=1,3                                  !hadron class
      do icdp=1,2                                 !diffractive eigenstate    
      do iv=1,11
       vvx=(iv-1)/10.d0                           !relative screening strenth

       do iy=1,51                                 !initialization
       do iz=1,11 
        qfanu(iy,iz,iv,icdp+2*(icz-1))=qlegi(iy,iz,iv,icz)
     *  +log(cd(icdp,icz))
       enddo
       enddo      
       do iy=2,51
        sy=spmax**(dble(iy-1)/50.d0)              !Pomeron mass squared
        rp=(rq(icz)+alfp*dlog(sy))*4.d0*.0389d0	
        do iz=2,11 
         qfanu(iy,iz,iv,icdp+2*(icz-1))=qfanu(iy-1,iz,iv,icdp+2*(icz-1))
        enddo
	
        n=1
3       n=n+1                                    !number of t-channel iterations
        nrep=0
        do iz=2,11  
         if(iz.gt.6)then
          z=.2d0*(iz-6)
          b=sqrt(-log(z)*rp)                      !impact parameter between ends
         else
          b=sqrt(-rp*(log(0.2d0)+.8d0*(iz-7)))
          z=exp(-b*b/rp)
         endif
         qxf=qgfan(sy,b*b,vvx,icdp,icz)
	 if(qxf.gt.0.d0)then
          qfan0(iy,iz,1)=dlog(qxf/z)
	 else
          qfan0(iy,iz,1)=qfan0(iy,iz-1,1)
	 endif
	 if(abs(qfan0(iy,iz,1)-qfanu(iy,iz,iv,icdp+2*(icz-1)))
     *   .gt.1.d-3)nrep=1
        enddo       
        do iz=2,11 
         qfanu(iy,iz,iv,icdp+2*(icz-1))=qfan0(iy,iz,1)
        enddo
        if(nrep.eq.1)goto 3
       enddo
      enddo
      enddo
      enddo
      
c-------------------------------------------------
c integrated cut fan contributions
      if(debug.ge.1)write (moniou,215)
      do icz=1,3                                  !hadron class
      do icdp=1,2                                 !diffractive eigenstate
c vvx,vvxp,vvxt - screening corrections from targ. and nuclear proj. fans
      do iv=1,11  
       vvx=dble(iv-1)/10.d0
      do iv1=1,6
       vvxp=dble(iv1-1)/5.d0
      do iv2=1,6
       vvxt=vvx*dble(iv2-1)/5.d0

       do iy=1,51                                !initialization
       do iz=1,11 
        if(iv2.eq.1)then
         do iqq=2,3
	  qfan(iy,iz,iv,iv1,icdp+2*(icz-1)+6*(iqq-2))
     *    =qlegi(iy,iz,iv,icz)+dlog(cd(icdp,icz))
	 enddo
        endif
        qfan(iy,iz,iv,iv1,12+iv2+6*(icdp+2*(icz-1)-1))
     *  =qlegi(iy,iz,iv,icz)+dlog(cd(icdp,icz))
	qfan(iy,iz,iv,iv1,12+iv2+6*(icdp+2*(icz-1)+5))
     *  =qfanu(iy,iz,iv,icdp+2*(icz-1))
	qfan(iy,iz,iv,iv1,12+iv2+6*(icdp+2*(icz-1)+11))
     *  =qfanu(iy,iz,iv,icdp+2*(icz-1))
       enddo
       enddo
       
       do iy=2,51
        sy=spmax**((iy-1)/50.d0)                  !Pomeron mass squared
        rp=(rq(icz)+alfp*log(sy))*4.d0*.0389d0	
        do iz=2,11 
         do iqq=4,6
          qfan(iy,iz,iv,iv1,12+iv2+6*(icdp-1+2*(icz-1)+6*(iqq-4)))
     *    =qfan(iy-1,iz,iv,iv1,12+iv2+6*(icdp-1+2*(icz-1)+6*(iqq-4)))
         enddo
	 if(iv2.eq.1)then
          do iqq=2,3
           qfan(iy,iz,iv,iv1,icdp+2*(icz-1)+6*(iqq-2))
     *     =qfan(iy-1,iz,iv,iv1,icdp+2*(icz-1)+6*(iqq-2))
          enddo
 	 endif
        enddo
	
        n=1
4       n=n+1                                    !number of t-channel iterations
        nrep=0
        do iz=2,11  
         if(iz.gt.6)then
          z=.2d0*(iz-6)
          b=sqrt(-log(z)*rp)                      !impact parameter between ends
         else
          b=sqrt(-rp*(log(0.2d0)+.8d0*(iz-7)))
          z=exp(-b*b/rp)
         endif
         call qgfan1(fann,sy,b*b,vvx,vvxp,vvxt,icdp,icz)
         do iqq=2,6
          if(fann(iqq).gt.0.d0)then
           qfan0(iy,iz,iqq)=dlog(fann(iqq)/z)
	  else
           qfan0(iy,iz,iqq)=qfan0(iy,iz-1,iqq)
	  endif
	 enddo
	 do iqq=4,6
	  if(abs(qfan0(iy,iz,iqq)-qfan(iy,iz,iv,iv1,12+iv2
     *    +6*(icdp-1+2*(icz-1)+6*(iqq-4)))).gt.1.d-3)nrep=1
         enddo
         if(iv2.eq.1)then
	  do iqq=2,3
           if(abs(qfan0(iy,iz,iqq)-qfan(iy,iz,iv,iv1
     *     ,icdp+2*(icz-1)+6*(iqq-2))).gt.1.d-3)nrep=1
          enddo
	 endif
        enddo      
        do iz=2,11 
         do iqq=4,6 
          qfan(iy,iz,iv,iv1,12+iv2+6*(icdp-1+2*(icz-1)+6*(iqq-4)))
     *    =qfan0(iy,iz,iqq)
         enddo
         if(iv2.eq.1)then
          do iqq=2,3 
           qfan(iy,iz,iv,iv1,icdp+2*(icz-1)+6*(iqq-2))=qfan0(iy,iz,iqq)
          enddo
         endif
        enddo
        if(nrep.eq.1)goto 4
		
        do iqq=2,6
        do iz=2,11  
         if(iz.gt.6)then
          z=.2d0*(iz-6)
          b=sqrt(-log(z)*rp)
         else
          b=sqrt(-rp*(log(0.2d0)+.8d0*(iz-7)))
          z=exp(-b*b/rp)
         endif
        enddo
        enddo
       enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      
c integrated cut Pomeron eikonals
      if(debug.ge.1)write (moniou,216)
      do icz=1,3                                    !hadron class
      do iy=1,11
       e0n=10.d0**iy                                !interaction energy
       sy=2.d0*e0n*am(2)+am(2)**2+am(icz)**2        
       rp=(rq(icz)+rq(2)+alfp*log(sy))*4.d0*.0389d0
      do iz=1,11  
       if(iz.gt.6)then
        z=.2d0*(iz-6)
        b=sqrt(-log(z)*rp)                          !impact parameter
       elseif(iz.gt.1)then
        b=sqrt(-rp*(log(0.2d0)+.8d0*(iz-7)))
        z=exp(-b*b/rp)
       endif

       vsoft=fp(icz)*fp(2)*sigs*sy**dels/rp*4.d0*.0389d0
     * *gamfun(real(1.d0-alpd+dels))**2*gamfun(real(1.d0+ahl(icz)))
     * *gamfun(real(1.d0+ahl(2)))/gamfun(real(2.d0-alpd+dels+ahl(icz)))
     * /gamfun(real(2.d0-alpd+dels+ahl(2)))
                                       !soft Pomeron eikonal (without screening)
       if(iz.ne.1)then
        vgg=qgfsh(sy,b*b,icz,0)        !gg-Pomeron eikonal (without screening)
        vqg=qgfsh(sy,b*b,icz,1)        !qg-Pomeron eikonal (without screening)
        vgq=qgfsh(sy,b*b,icz,2)        !gq-Pomeron eikonal (without screening)
        vqq=qghard(sy,b*b,icz)         !qq-Pomeron eikonal (without screening)
       endif
       
       do icdp=1,2                              !proj.diffractive eigenstate       
       do icdt=1,2                              !targ.diffractive eigenstate       
        if(iz.eq.1)then
         qpomu(iy,iz,icdp,icdt,icz)=log(vsoft*cd(icdp,icz)*cd(icdt,2))
         do iv=1,6
          qpom1(iy,iz,iv,icdp+2*(icdt-1),icz)
     *    =log(vsoft*cd(icdp,icz)*cd(icdt,2))
         do iv1=1,6
         do iv2=1,6
         do iv3=1,6
          qpomr(iy,iz,iv,iv1+6*(iv2-1)+36*(iv3-1)
     *    ,icdp+2*(icdt-1)+4*(icz-1))=-100.d0
	 enddo
	 enddo
	 enddo
	 enddo
	else
         do iv=1,6  
          vvxp=(iv-1)/5.d0                       !relative screening strenth	  	  

          vtot=qgftot(sy,b,vvxp,icdp,icdt,icz)
          qxp1=vtot+vqq*cd(icdp,icz)*cd(icdt,2)  !1-Pomeron eikonal (with scr.)
          qpom1(iy,iz,iv,icdp+2*(icdt-1),icz)=log(qxp1/z)
	 
         do iv1=1,6                 !screening factors for h-A (AA) interactions
          vvxt=(iv1-1)/5.d0
         do iv2=1,6
          vvxpa=(iv2-1)/5.d0
         do iv3=1,6
          vvxta=(iv3-1)/5.d0

          v3p=qg3pa(sy,b,vvxp,vvxt,vvxpa,vvxta,icdp,icdt,icz) 
	                                 !contribution of multi-Pomeron vertexes
          if(v3p.gt.0.d0)then
           qpomr(iy,iz,iv,iv1+6*(iv2-1)+36*(iv3-1)
     *     ,icdp+2*(icdt-1)+4*(icz-1))=log(v3p/z)
          else
           qpomr(iy,iz,iv,iv1+6*(iv2-1)+36*(iv3-1),icdp+2*(icdt-1)
     *     +4*(icz-1))=qpomr(iy,iz-1,iv,iv1+6*(iv2-1)+36*(iv3-1)
     *     ,icdp+2*(icdt-1)+4*(icz-1))
	  endif
     
          if(iv.eq.1.and.iv1.eq.1.and.iv2.eq.1.and.iv3.eq.1)then
	   qxp2=qxp1+v3p                           !total cut eikonal
           qpomu(iy,iz,icdp,icdt,icz)=log(qxp2/z)
	   
           if(debug.ge.2)write(moniou,217)ty(icz)
     *     ,e0n,icdp,icdt,qxp2,qxp1,v3p
	  endif
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

c interaction cross sections
      ia(1)=1
      do iy=1,11        
       e0n=10.d0**iy                                      !interaction energy
       scm=2.d0*e0n*am(2)+am(2)**2+am(icz)**2

       do iiz=1,3  
        icz=iiz                                           !hadron class
        rp=(rq(icz)+rq(2)+alfp*log(scm))*4.d0*.0389d0
        do iia=1,6
         if(iia.le.4)then
	  ia(2)=4**(iia-1)                                !target mass number
	 elseif(iia.eq.5)then
	  ia(2)=14
	 else
	  ia(2)=40
	 endif	 
         if(debug.ge.1)write (moniou,218)e0n,ty(icz),ia(2)
c nuclear densities
         if(ia(2).lt.10)then                       !light nuclei - gaussian
          rnuc(2)=.9d0*float(ia(2))**.3333         !nuclear radius
          wsnuc(2)=amws                            !not used
	  wbnuc(2)=0.d0                            !not used
         elseif(ia(2).le.56)then                   !3-parameter Fermi
          rnuc(2)=trnuc(ia(2))                     !nuclear radius
          wsnuc(2)=twsnuc(ia(2))                   !diffuseness
	  wbnuc(2)=twbnuc(ia(2))                   !"wine-bottle" parameter
	 else                                      !2-parameter Fermi
          rnuc(2)=1.19*float(ia(2))**(1./3.)-1.38*float(ia(2))**(-1./3.) !radius
          wsnuc(2)=amws                            !diffuseness
	  wbnuc(2)=0.d0
         endif

         if(ia(2).eq.1)then                        !hadron-proton interaction
          call qgfz(0.d0,gz0,0,0)
          gtot=gz0(1)                              !total cross-section
          gin=(gz0(2)+gz0(3)+gz0(4))*.5d0          !inelastic cross section
          bel=gz0(5)                               !elastic scattering slope
          gel=gtot-gin                             !elastic cross-section
          gdp=gz0(3)*.5d0                          !proj. diffr. cross section
          gdt=gz0(4)*.5d0                          !targ. diffr. cross section
          if(iy.le.10)gsect(iy,icz,iia)=log(gin)
          if(debug.ge.1)write (moniou,219)gtot,gin,gel,gdp,gdt,bel
	  
         else                                      !hadron-nucleus interaction
          bm=rnuc(2)+dlog(29.d0)*wsnuc(2)          !impact parameter for integr.
	  anorm=qganrm(rnuc(2),wsnuc(2),wbnuc(2))*rp    !density normalization
          call qggau(gz1)                          !integration over b< bm
          call qggau1(gz1)                         !integration over b> bm
          gin=gz1(1)+gz1(2)+gz1(3)                 !inelastic cross section
          if(iy.le.10)gsect(iy,icz,iia)=log(gin*10.d0)
          if(debug.ge.1)write (moniou,220)
     *    gin*10.d0,gz1(3)*10.d0,gz1(2)*10.d0
         endif
        enddo
       enddo
      enddo

c-----------------------------------------------------------------------------
c timelike Sudakov formfactor
      if(debug.ge.1)write (moniou,221)
      do m=1,2                       !parton type (1-g, 2-q)
       fsud(1,m)=0.d0
      do k=2,10
       qmax=qtf*4.d0**(1.d0+k)       !effective virtuality (qt**2/z**2/(1-z)**2)
       fsud(k,m)=qgsudt(qmax,m)
      enddo
      enddo
c-----------------------------------------------------------------------------
c effective virtuality (used for inversion in timelike branching)
      if(debug.ge.1)write (moniou,222)
      do m=1,2                       !parton type (1-g, 2-q)
      do k=1,10
       qlmax=1.38629d0*(k-1)
       qrt(k,1,m)=0.d0
       qrt(k,101,m)=qlmax
      do i=1,99                      !bins in Sudakov formfactor
       if(k.eq.1)then
        qrt(k,i+1,m)=0.d0
       else
        qrt(k,i+1,m)=qgroot(qlmax,.01d0*i,m)
       endif
      enddo
      enddo
      enddo
      
c-------------------------------------------------
c cut Pomeron leg eikonals (semi-hard)
      if(debug.ge.1)write (moniou,223)
      s2min=4.d0*fqscal*qt0
      do icz=1,3                                           !hadron class   
      do iy=1,11
       sy=2.d0*s2min*(spmax/s2min/2.d0)**((iy-1)/10.d0)    !Pomeron mass squared
      do icdp=1,2                                          !diffr. eigenstate      
      do iv=1,11
       vvx=(iv-1)/10.d0                                    !relative scr. strenth
       
       do ix=1,10
        if(ix.le.5)then
         xp=.2d0*(sy/spmax/2.d0)**((6-ix)/5.d0)            !Pomeron end momentum
        else
         xp=.2d0*(ix-5)
        endif
        rp=(rq(icz)+alfp*log(sy/xp))*4.d0*.0389d0
       do iz=2,11  
        if(iz.gt.6)then
         z=.2d0*(iz-6)
         b=sqrt(-log(z)*rp)                                !impact parameter
        elseif(iz.gt.1)then
         b=sqrt(-rp*(log(0.2d0)+.8d0*(iz-7)))
         z=exp(-b*b/rp)
        endif
       
        qxg=qglsh(sy,xp,b*b,vvx,icdp,icz,0,1)              !gg-leg
        if(qxg.gt.0.d0)then
         qlegc(iy,ix,iz,iv,icdp+2*(icz-1)+12)=log(qxg/sy**delh/z) 
	else
         qlegc(iy,ix,iz,iv,icdp+2*(icz-1)+12)
     *   =qlegc(iy,ix,iz,iv-1,icdp+2*(icz-1)+12) 
	endif  
 	qxq=qglsh(sy,xp,b*b,vvx,icdp,icz,1,1)              !qg-leg    
        if(qxq.gt.0.d0)then
         qlegc(iy,ix,iz,iv,icdp+2*(icz-1)+18)=log(qxq/sy**delh/z)       
	else
         qlegc(iy,ix,iz,iv,icdp+2*(icz-1)+18)
     *   =qlegc(iy,ix,iz,iv-1,icdp+2*(icz-1)+18)
	endif  
       enddo
       enddo
      enddo
      enddo
      enddo
      enddo

c-------------------------------------------------
c cut Pomeron leg eikonals (soft and total)
      if(debug.ge.1)write (moniou,224)
      do icz=1,3                                          !hadron class 
      do iy=1,11
       sy=spmax**((iy-1)/10.d0)                           !Pomeron mass squared
      do icdp=1,2                                         !diffr. eigenstate       
      do iv=1,11
       vvx=(iv-1)/10.d0                                   !relative scr. strenth      
      do ix=1,10
       if(ix.le.5)then
        xp=.2d0*(sy/spmax/2.d0)**((6-ix)/5.d0)            !Pomeron end momentum
       else
        xp=.2d0*(ix-5)
       endif
       rp=(rq(icz)+alfp*log(sy/xp))*4.d0*.0389d0
      do iz=1,11  
       if(iz.gt.6)then
        z=.2d0*(iz-6)
        b=sqrt(-log(z)*rp)                                !impact parameter
       elseif(iz.gt.1)then
        b=sqrt(-rp*(log(0.2d0)+.8d0*(iz-7)))
        z=exp(-b*b/rp)
       endif
       
       if(iy.eq.1.or.iz.eq.1)then
        qlegc(iy,ix,iz,iv,icdp+2*(icz-1)) 
     *  =log(fp(icz)*sigs*g3p/rp*4.d0*.0389d0*cd(icdp,icz))
        qlegc(iy,ix,iz,iv,icdp+2*(icz-1)+6)
     *  =qlegc(iy,ix,iz,iv,icdp+2*(icz-1)) 
       else
        vg=qglegc(sy,xp,b*b,vvx,icdp,icz,1)
        vq=qglegc(sy,xp,b*b,vvx,icdp,icz,2)

        qxs=qgls(sy,xp,b*b,vvx,icdp,icz,1)
        if(qxs.gt.0.d0)then
         qlegc(iy,ix,iz,iv,icdp+2*(icz-1)+6)=log(qxs/sy**dels/z)  !soft eikonal 
         qlegc(iy,ix,iz,iv,icdp+2*(icz-1))                        !total eikonal
     *   =log((qxs+(vg+vq)*xp**alpd)/sy**dels/z)
	else
         qlegc(iy,ix,iz,iv,icdp+2*(icz-1)+6)
     *   =qlegc(iy,ix,iz,iv-1,icdp+2*(icz-1)+6) 
         qlegc(iy,ix,iz,iv,icdp+2*(icz-1)) 
     *   =log((vg+vq)*xp**alpd/sy**dels/z)
	endif
       endif
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      
c-------------------------------------------------
c cut Pomeron eikonals (semi-hard)
      if(debug.ge.1)write (moniou,225)
      do icz=1,3                                          !proj. class
      do iy=1,11
       sy=2.d0*s2min*(spmax/s2min/2.d0)**((iy-1)/10.d0)   !Pomeron mass squared
      do ix1=1,10
       if(ix1.le.5)then
        xp=.2d0*(sy/spmax/2.d0)**((6-ix1)/5.d0)           !Pomeron LC+ momentum
       else
        xp=.2d0*(ix1-5)
       endif
      do ix2=1,10
       if(ix2.le.5)then
        xm=.2d0*(sy/spmax/2.d0)**((6-ix2)/5.d0)           !Pomeron LC- momentum
       else
        xm=.2d0*(ix2-5)
       endif
       rp=(rq(icz)+rq(2)+alfp*log(sy/xp/xm))*4.d0*.0389d0
      do iz=2,11  
       if(iz.gt.6)then
        z=.2d0*(iz-6)
        b=sqrt(-log(z)*rp)                                !impact parameter
       elseif(iz.gt.1)then
        b=sqrt(-rp*(log(0.2d0)+.8d0*(iz-7)))
        z=exp(-b*b/rp)
       endif
      do iv=1,11
       vvx=(iv-1)/10.d0                                   !relative scr. strenth
       
       do icdp=1,2       
       do icdt=1,2 
        vgg=qgpsh(sy,xp,xm,b,vvx,icdp,icdt,icz,0)          !gg-Pomeron
        vqg=qgpsh(sy,xp,xm,b,vvx,icdp,icdt,icz,1)          !qg-Pomeron
        vgq=qgpsh(sy,xp,xm,b,vvx,icdp,icdt,icz,2)          !gq-Pomeron
	if(vgg.gt.0.d0)then         
         qpomc(iy,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+24)
     *   =log(vgg/sy**delh/z)
        else
         qpomc(iy,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+24)
     *   =qpomc(iy-1,ix1+10*(ix2-1),iv,iz,icdp+2*(icdt-1)+4*(icz-1)+24)
	endif
        if(vqg.gt.0.d0)then         
         qpomc(iy,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+36)
     *   =log(vqg/sy**delh/z)
        else
         qpomc(iy,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+36)
     *   =qpomc(iy-1,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+36)
	endif
	if(vgq.gt.0.d0)then         
         qpomc(iy,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+48)
     *   =log(vgq/sy**delh/z)	
        else
         qpomc(iy,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+48)
     *   =qpomc(iy-1,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+48)
	endif
       enddo
       enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo      

c-------------------------------------------------
c cut Pomeron eikonals (soft and total)
      if(debug.ge.1)write (moniou,226)
      do icz=1,3                                          !proj. class
      do iy=1,11
       sy=spmax**((iy-1)/10.d0)                           !Pomeron mass squared
      do ix1=1,10
       if(ix1.le.5)then
        xp=.2d0*(sy/spmax/2.d0)**((6-ix1)/5.d0)           !Pomeron LC+ momentum
       else
        xp=.2d0*(ix1-5)
       endif
      do ix2=1,10
       if(ix2.le.5)then
        xm=.2d0*(sy/spmax/2.d0)**((6-ix2)/5.d0)           !Pomeron LC- momentum
       else
        xm=.2d0*(ix2-5)
       endif
       rp=(rq(icz)+rq(2)+alfp*log(sy/xp/xm))*4.d0*.0389d0
      do iz=1,11  
       if(iz.gt.6)then
        z=.2d0*(iz-6)
        b=sqrt(-log(z)*rp)                                !impact parameter
       elseif(iz.gt.1)then
        b=sqrt(-rp*(log(0.2d0)+.8d0*(iz-7)))
        z=exp(-b*b/rp)
       endif
      do iv=1,11
       vvx=(iv-1)/10.d0                                   !relative scr. strenth
       
       vsoft=sy**dels*fp(icz)*fp(2)*sigs/(xp*xm)**alpd/rp*4.d0*.0389d0 !soft Pom.
       if(iz.ne.1)vqq=qgpomc(sy,xp,xm,b*b,0.d0,1,1,icz,4) !qq-Pomeron
     * /cd(1,icz)/cd(1,2)
      
       do icdp=1,2       
       do icdt=1,2
        if(iy.eq.1.or.iz.eq.1)then
         qpomc(iy,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+12)
     *   =log(vsoft*cd(icdp,icz)*cd(icdt,2)/sy**dels*(xp*xm)**alpd)
         qpomc(iy,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1))
     *   =qpomc(iy,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+12)
	else
         vgg=qgpomc(sy,xp,xm,b*b,vvx,icdp,icdt,icz,1)     !gg-Pomeron
	 vqg=qgpomc(sy,xp,xm,b*b,vvx,icdp,icdt,icz,2)     !qg-Pomeron
	 vgq=qgpomc(sy,xp,xm,b*b,vvx,icdp,icdt,icz,3)     !gq-Pomeron
	 qps=qgpsoft(sy,xp,xm,b,vvx,icdp,icdt,icz)
	 if(qps.gt.0.d0)then         
          qpomc(iy,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+12)
     *    =log(qps/sy**dels*(xp*xm)**alpd/z)
         else
          qpomc(iy,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+12)
     *    =qpomc(iy-1,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1)+12)
	 endif
	 qpt=vqq*cd(icdp,icz)*cd(icdt,2)+max(0.d0,qps)+vgg+vqg+vgq !total eikonal
	 qpomc(iy,ix1+10*(ix2-1),iz,iv,icdp+2*(icdt-1)+4*(icz-1))
     *   =log(qpt/sy**dels*(xp*xm)**alpd/z)
        endif
       enddo
       enddo
      enddo
      enddo
      enddo
      enddo
      enddo      
      enddo      

      if(debug.ge.1)write (moniou,227)      
      open(1,file='qgsdat-II-03',form='unformatted',status='unknown')
      write (1) csborn,cs0,cstot,evk,qppdi,qlegi,qfanu,qfan
     *,qpomu,qpom1,qpomr,qlegc,qpomc,qpomi,qpomi0,gsect,fsud,qrt
      close(1)
      
5     continue
c Nuclear cross sections
      inquire(file=DATDIR(1:INDEX(DATDIR,' ')-1)//'sectnu-II-03',
     *        exist=lcalc)
      if(lcalc)then
       if(debug.ge.1)write (moniou,228)       
       open(2,file=DATDIR(1:INDEX(DATDIR,' ')-1)//'sectnu-II-03',
     *       status='old')
       read (2,*)asect
       close(2)

      else
       niter=5000                               !number of iterations
       do ie=1,10
        e0n=10.d0**ie                           !interaction energy (per nucleon)
       do iia1=1,6
        iap=2**iia1                             !proj. mass number
       do iia2=1,6
        if(iia2.le.4)then
         iat=4**(iia2-1)                        !targ. mass number
	elseif(iia2.eq.5)then
	 iat=14
	else
	 iat=40
	endif
        if(debug.ge.1)write (moniou,229)e0n,iap,iat
	
        call qgini(e0n,2,iap,iat)
        call qgcrossc(niter,gtot,gprod,gabs,gdd,gqel,gcoh)
        if(debug.ge.1)write (moniou,230)gtot,gprod,gabs,gdd,gqel,gcoh
        asect(ie,iia1,iia2)=log(gprod)
       enddo
       enddo
       enddo
       open(2,file='sectnu-II-03',status='unknown')
       write (2,*)asect
       close(2)
      endif
      if(debug.ge.1)write (moniou,231)

201   format(2x,'qgaini - main initialization procedure'
     *,' option=',i1)
202   format(2x,'qgaini: normalization of parton densities'
     */4x,a7,2x,'norm=',e10.3,2x,'<x_g>=',e10.3,2x,'<x_q_s>=',e10.3
     */4x,'fp(icz)=',e10.3,2x,'rr=',e10.3)
203   format(2x,'qgaini: hard cross section ratios readout from the'
     *,' file qgsdat-II-03')
204   format(2x,'qgaini: hard cross sections calculation')
205   format(2x,'qgaini: number of rungs considered:',i2
     */4x,'starting energy bin for ordered and general ladders:',3i4)
206   format(2x,'qgaini: new and old cross sections'
     *,' for ordered ladder:'/4x,e10.3,3x,e10.3)
207   format(2x,'qgaini: new and old cross sections'
     *,' for general ladder:'/4x,e10.3,3x,e10.3)
208   format(2x,'qgaini: hard cross sections finished')
209   format(2x,'qgaini: parton distributions in the Pomeron')
210   format(2x,'qgaini: eikonal for an itermediate gg-Pomeron')
211   format(2x,'qgaini: integrated eikonal for an itermediate'
     *,' gg-Pomeron')
212   format(2x,'qgaini: integrated Pomeron leg eikonals')
213   format(2x,'qgaini: integrated fan contributions')
214   format(2x,'qgaini: integrated cut Pomeron leg eikonals')
215   format(2x,'qgaini: integrated cut fan contributions')
216   format(2x,'qgaini: integrated cut Pomeron eikonals')
217   format(2x,'qgspsaini: ',a7,'-proton interaction'
     */4x,' e0=',e10.3,' icdp=',i1,' icdt=',i1,2x,'total eikonal: '
     *,e10.3/4x,'1-Pomeron eikonal: ',e10.3,2x,'multi-Pomeron: '
     *,e10.3)
218   format(2x,'qgaini: initial particle energy:',e10.3,2x
     *,'its type:',a7,2x,'target mass number:',i2)
219   format(2x,'qgaini: hadron-proton cross sections:'
     */4x,'gtot=',e10.3,2x,'gin=',e10.3,2x,'gel=',e10.3/4x
     *,'gdifrp=',e10.3,2x,'gdifrt=',e10.3,2x,'b_el=',e10.3,2x)
220   format(2x,'qgaini: hadron-nucleus cross sections:'
     */4x,'gin=',e10.3,2x,'gdifr_targ=',e10.3,2x
     *,'gdifr_proj=',e10.3)
221   format(2x,'qgaini: timelike Sudakov formfactor')
222   format(2x,'qgaini: effective virtuality for inversion')
223   format(2x,'qgaini: cut Pomeron leg eikonals (semi-hard)')
224   format(2x,'qgaini: cut Pomeron leg eikonals (soft and total)')
225   format(2x,'qgaini: cut Pomeron eikonals (semi-hard)')
226   format(2x,'qgaini: cut Pomeron eikonals (soft and total)')
227   format(2x,'qgaini: cross sections are written to the file'
     *,' qgsdat-II-03')
228   format(2x,'qgaini: nuclear cross sections readout from the file'
     *,' sectnu-II-03')
229   format(2x,'qgaini: initial nucleus energy:',e10.3,2x
     *,'projectile mass:',i2,2x,'target mass:',i2)
230   format(2x,'gtot',d10.3,'  gprod',d10.3,'  gabs',d10.3
     */2x,'gdd',d10.3,'  gqel',d10.3,' gcoh',d10.3)
231   format(2x,'qgaini - end')
232   format(2x,'qgaini: qcd evolution - ',i2,'-th order contribution')
      return
      end

c=============================================================================
      subroutine qgini(e0n,icp0,iap,iat)
c-----------------------------------------------------------------------------
c qgini - additional initialization procedure (before each run)
c e0n  - interaction energy (per hadron/nucleon),
c icp0 - hadron type (+-1 - pi+-, +-2 - p(p~), +-3 - n(n~),
c                     +-4 - K+-, +-5 - K_l/s),
c iap  - projectile mass number (1 - for a hadron),
c iat  - target mass number 
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207)
cdh   dimension wk(3),wa(3)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr4/  ey0(3)
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),b
      common /qgarr10/ am(7),ammu
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr43/ moniou
      common /qgarr49/ trnuc(56),twsnuc(56),twbnuc(56)
      common /debug/   debug

      if(debug.ge.1)write (moniou,201)icp0,iap,iat,e0n
      icp=icp0
      ia(1)=iap
      ia(2)=iat
      
c icz - primary particle class (1- pion, 2 - nucleon, 3 - kaon)
      if(iabs(icp).eq.12)then
       icz=6
      elseif(iabs(icp).lt.6)then
       icz=iabs(icp)/2+1
      else
       icz=(iabs(icp)+1)/2
      endif

      scm=2.d0*e0n*am(2)+am(2)**2+am(icz)**2                !c.m. energy squared
      ey0(1)=dsqrt(scm)/(e0n+am(2)+dsqrt(e0n**2-am(icz)**2))!Lorentz boost to lab.
      ey0(2)=1.d0
      ey0(3)=1.d0
      wp0=(e0n+dsqrt(e0n**2-am(icz)**2))*ey0(1)             !initial LC+ mometum
      wm0=am(2)/ey0(1)                                      !initial LC- mometum

c nuclear radii and weights for sampling nuclear configurations
      do i=1,2
       if(ia(i).lt.10.and.ia(i).ne.1)then     !gaussian density for light nuclei
        rnuc(i)=.9d0*float(ia(i))**.3333      !nuclear radius
        if(ia(i).eq.2)rnuc(i)=3.16d0
        rnuc(i)=rnuc(i)*dsqrt(2.d0*ia(i)/(ia(i)-1.d0))
                                  !rnuc -> rnuc*a/(a-1) - to use Van-Hove method
       elseif(ia(i).ne.1)then
        if(ia(i).le.56)then                   !3-parameter Fermi
         rnuc(i)=trnuc(ia(i))                 !nuclear radius
         wsnuc(i)=twsnuc(ia(i))               !diffuseness
	 wbnuc(i)=twbnuc(ia(i))               !"wine-bottle" parameter
	else                                  !2-parameter Fermi
         rnuc(i)=1.19*float(ia(i))**(1./3.)-1.38*float(ia(i))**(-1./3.) !radius
         wsnuc(i)=amws                        !diffuseness
	 wbnuc(i)=0.d0
	endif
        cr1(i)=1.d0+3.d0/rnuc(i)*wsnuc(i)+6.d0/(rnuc(i)/wsnuc(i))**2
     *  +6.d0/(rnuc(i)/wsnuc(i))**3
        cr2(i)=3.d0/rnuc(i)*wsnuc(i)
        cr3(i)=3.d0/rnuc(i)*wsnuc(i)+6.d0/(rnuc(i)/wsnuc(i))**2
       endif
      enddo

      if(ia(1).ne.1)then                              !primary nucleus
       bm=rnuc(1)+rnuc(2)+5.d0*max(wsnuc(1),wsnuc(2)) !impact parameter cutoff
      elseif(ia(2).ne.1)then                          !hadron-nucleus
       bm=rnuc(2)+5.d0*wsnuc(2)
      else                                            !hadron-proton
       bm=2.d0*dsqrt((rq(icz)+rq(2)+alfp*log(scm))*4.d0*.0398d0)
      endif
      if(debug.ge.2)write (moniou,202)
      
201   format(2x,'qgini - miniinitialization: particle type icp0=',
     *i2,2x,'projectile mass number iap=',i2/4x,
     *'target mass number iat=',i2,' interaction energy e0n=',e10.3)
202   format(2x,'qgini - end')
      return
      end

c=============================================================================
      double precision function qgftot(sy,b,vvx,icdp,icdt,icz)
c-----------------------------------------------------------------------------
c qgftot - 1-Pomeron eikonal (with screening)
c sy   - Pomeron mass squared,
c b    - impact parameter,
c vvx  - relative screening strength,
c icdp - projectile diffractive eigenstate,
c icdt - target diffractive eigenstate,
c icz  - hadron class
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr19/ ahl(3)
      common /qgarr25/ ahv(3)
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /qgarr51/ x4(2),a4(2)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)sy,b,icdp,icdt,icz
      qgftot=0.d0      
      do i1=1,2
      do m1=1,2
       tp=1.d0-(.5d0+x4(i1)*(m1-1.5d0))**(1./(1.+dels-alpd))
       xp=1.d0-tp**(1.d0/(1.d0+ahl(icz)))                          !LC+
      do i2=1,2
      do m2=1,2
       tm=1.d0-(.5d0+x4(i2)*(m2-1.5d0))**(1./(1.+dels-alpd))
       xm=1.d0-tm**(1.d0/(1.d0+ahl(2)))                            !LC-
       
       ws=qgpsoft(xp*xm*sy,xp,xm,b,vvx,icdp,icdt,icz)              !soft Pomeron
       wgg=qgpsh(xp*xm*sy,xp,xm,b,vvx,icdp,icdt,icz,0)/(xp*xm)**alpd !gg-Pomeron
       wqg=qgpsh(xp*xm*sy,xp,xm,b,vvx,icdp,icdt,icz,1)               !qg-Pomeron
     * /xm**alpd/dsqrt(xp)*(1.d0-xp)**(ahv(icz)-ahl(icz))            !gq-Pomeron
       wgq=qgpsh(xp*xm*sy,xp,xm,b,vvx,icdp,icdt,icz,2)
     * /xp**alpd/dsqrt(xm)*(1.d0-xm)**(ahv(2)-ahl(2))

       qgftot=qgftot+a4(i1)*a4(i2)*(ws+wgg+wqg+wgq)
     * /((1.d0-tp)*(1.d0-tm))**(dels-alpd)
      enddo
      enddo
      enddo
      enddo
      qgftot=qgftot/4.d0/(1.+ahl(icz))/(1.+ahl(2))/(1.d0+dels-alpd)**2
      if(debug.ge.4)write (moniou,202)qgftot
          
201   format(2x,'qgftot - 1-Pomeron eikonal:'/4x,'sy=',e10.3
     *,2x,'b=',e10.3,2x,'icdp=',i1,2x,'icdt=',i1,2x,'icz=',i1)
202   format(2x,'qgftot=',e10.3)
      return
      end
      
c=============================================================================
      double precision function qgfsh(sy,bb,icz,iqq)
c-----------------------------------------------------------------------------
c qgfsh - semihard interaction eikonal (without screening)
c sy  - Pomeron mass squared,
c bb  - squared impact parameter,
c icz - hadron class,
c iqq - type of the hard interaction (0-gg, 1-q_vg, 2-gq_v)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr25/ ahv(3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)sy,bb,iqq,icz
      qgfsh=0.d0      
      s2min=4.d0*fqscal*qt0                !energy threshold for 2->2 process
      xmin=s2min/sy
      if(xmin.ge.1.d0)return
      xmin=xmin**(delh-dels)
      if(iqq.eq.1)then
       icv=icz
       icq=2
      elseif(iqq.eq.2)then
       icv=2
       icq=icz
      endif
      if(debug.ge.5)write (moniou,202)xmin,iqq

c numerical integration over z1
      do i=1,7
      do m=1,2
       z1=(.5d0*(1.d0+xmin-(2*m-3)*x1(i)*(1.d0-xmin)))
     * **(1.d0/(delh-dels))
       ww=z1*sy                            !energy squared for hard process
       sjqq=qgjit(qt0,qt0,ww,2,2)          !qq-ladder cross section
       sjqg=qgjit(qt0,qt0,ww,1,2)          !qg-ladder cross section
       sjgg=qgjit(qt0,qt0,ww,1,1)          !gg-ladder cross section
       if(debug.ge.5)write (moniou,203)ww,sjqq+sjqg+sjgg
        
       if(iqq.eq.0)then                    !gg-ladder
        st2=0.d0
        do j=1,7
        do k=1,2
         xx=.5d0*(1.d0+x1(j)*(2*k-3))
         xp=z1**xx
         xm=z1/xp
         glu1=qgftld(xp,icz)              !gluon distribution in the Pomeron
         sea1=qgftle(xp,icz)              !sea quark distribution in the Pomeron
         glu2=qgftld(xm,2)
         sea2=qgftle(xm,2)
         st2=st2+a1(j)*(glu1*glu2*sjgg+(glu1*sea2+glu2*sea1)*sjqg
     *   +sea1*sea2*sjqq)
        enddo
        enddo
        rh=rq(icz)+rq(2)-alfp*dlog(z1)
        qgfsh=qgfsh-a1(i)*dlog(z1)/z1**delh*st2
     *  *exp(-bb/(4.d0*.0389d0*rh))/rh
      
       else                                !gq_v-ladder
        st2=0.d0
        alh=.5d0+dels
        xam=z1**alh

        do j=1,7
        do k=1,2
         xp=(.5d0*(1.d0+xam+x1(j)*(2*k-3)*(1.d0-xam)))**(1.d0/alh)
         xm=z1/xp
         glu=qgftld(xm,icq)
         sea=qgftle(xm,icq)
         rh=rq(icz)+rq(2)-alfp*dlog(xm)
         fst=(glu*sjqg+sea*sjqq)*(1.d0-xp)**ahv(icv)
     *   *(qggrv(xp,qt0,icv,1)+qggrv(xp,qt0,icv,2))/dsqrt(xp)
     *   *exp(-bb/(4.d0*.0389d0*rh))/rh
         st2=st2+a1(j)*fst
        enddo
        enddo
        st2=st2*(1.d0-xam)/alh
        qgfsh=qgfsh+a1(i)/z1**delh*st2
       endif
      enddo
      enddo

      if(iqq.eq.0)then
       qgfsh=qgfsh*rr**2*(1.d0-xmin)/(delh-dels)*fp(icz)*fp(2)*factk
     * /2.d0*pi
      else
       qgfsh=qgfsh*rr*fp(icq)*(1.d0-xmin)/(delh-dels)*factk/8.d0
      endif
      if(debug.ge.4)write (moniou,204)qgfsh
      
201   format(2x,'qgfsh - semihard interaction eikonal:'
     */4x,'sy=',e10.3,2x,'bb=',e10.3,2x,'iqq=',i1,2x,'icz=',i1)
202   format(2x,'qgfsh:',2x,'xmin=',e10.3,2x,'iqq=',i3)
203   format(2x,'qgfsh:',2x,'s_hard=',e10.3,2x,'sigma_hard=',e10.3)
204   format(2x,'qgfsh=',e10.3)
      return
      end
     
c=============================================================================
      double precision function qgleg(sy,bb,vvx,icdp,icz)
c-----------------------------------------------------------------------------
c qgleg - integrated Pomeron leg eikonal
c sy   - Pomeron mass squared,
c bb   - squared impact parameter,
c vvx  - relative strenth of screening corrections (0<vvx<1),
c icdp - diffractive eigenstate,
c icz  - hadron class,
c jj=0 - screening from current proj. hadron (nucleon) included in vvx,
c jj=1 - treated explicitely
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr19/ ahl(3)
      common /qgarr25/ ahv(3)
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)sy,bb,vvx,icdp,icz         !so161205
      qgleg=0.d0
      if(sy.lt.1.001d0)then
       tmin=1.d0
      else
       tmin=(1.d0-(1.d0-1.d0/sy)**(1.+ahl(icz)))**(1.+dels-alpd)
      endif
      if(debug.ge.5)write (moniou,202)tmin
      do i1=1,7
      do m1=1,2
       tp=1.d0-(.5d0*(1.d0+tmin)+x1(i1)*(m1-1.5d0)*(1.d0-tmin))
     * **(1./(1.+dels-alpd))
       if(tp.gt.1.d-9)then
        xp=1.d0-tp**(1.d0/(1.d0+ahl(icz)))                !Pomeron end momentum
       else
        xp=1.d0
       endif       
       ws=qgls(xp*sy,xp,bb,vvx,icdp,icz,0)/xp**alpd       !soft Pomeron
       if(bb.gt.1.d9)then    !asymptotics for b->infty
        wg=0.d0
	wq=0.d0
       else
        wg=qglsh(xp*sy,xp,bb,vvx,icdp,icz,0,0)/xp**alpd   !g-Pomeron
        wq=qglsh(xp*sy,xp,bb,vvx,icdp,icz,1,0)/dsqrt(xp)  !q-Pomeron
     *  *(1.d0-xp)**(ahv(icz)-ahl(icz))
       endif
       qgleg=qgleg+a1(i1)*(ws+wg+wq)/(1.d0-tp)**(dels-alpd)
      enddo
      enddo
      qgleg=qgleg/2.d0/(1.+ahl(icz))/(1.d0+dels-alpd)
      if(debug.ge.4)write (moniou,203)qgleg
     
201   format(2x,'qgleg - Pomeron leg eikonal:'/4x,'s=',e10.3,2x,'b^2='
     *,e10.3,2x,'vvx=',e10.3,2x,'icdp=',i1,2x,'icz=',i1)        !so161205
202   format(2x,'qgleg:',2x,'tmin=',e10.3)
203   format(2x,'qgleg=',e10.3)
      return
      end
      
c------------------------------------------------------------------------
      double precision function qglegi(sy,bb,vvx,icdp,icz)
c-----------------------------------------------------------------------
c qglegi - integration of Pomeron leg eikonal
c sy   - Pomeron mass squared,
c bb   - squared impact parameter,
c vvx  - relative strenth of screening corrections (0<vvx<1),
c icdp - diffractive eigenstate,
c icz  - hadron class,
c jj=0 - screening from current proj. hadron (nucleon) included in vvx,
c jj=1 - treated explicitely
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
cdh   dimension wk(3),wj(3),wz(3),wi(3)
      dimension wk(3),wj(3),wz(3)
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr19/ ahl(3)
      common /qgarr20/ spmax
      common /qgarr27/ qlegi(51,11,11,3),qfanu(51,11,11,6)
     *,qfan(51,11,11,6,120)
      common /qgarr43/ moniou
      common /debug/   debug
     
      if(debug.ge.3)write (moniou,201)sy,bb,vvx,icdp,icz
      qglegi=0.d0
      if(sy.le.1.d0)return

      rp=(rq(icz)+alfp*log(max(1.d0,sy)))*4.d0*.0389d0
      z=exp(-bb/rp)
      if(z.lt..2d0*exp(-4.d0))then
       izmax=2
       jz=1
       wz(2)=5.d0*z*exp(4.d0)
       wz(1)=1.d0-wz(2)
      else
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=(-bb/rp-log(0.2d0))/.8d0+7.d0
       endif
       jz=min(10,int(zz))
       jz=max(2,jz)
       wz(2)=zz-jz
       wz(1)=1.d0-wz(2)
       izmax=2
      endif
      
      yl=log(sy)/log(spmax)*50.d0+1.d0
      k=max(1,int(yl))
      k=min(k,50)     
      wk(2)=yl-k
      wk(1)=1.d0-wk(2)
      iymax=2

      vl=max(1.d0,vvx*10.d0+1.d0)
      if(vvx.eq.0.d0)then
       j=1
       ivmax=1
       wj(1)=1.d0
      elseif(vl.lt.2.d0)then
       j=1    
       wj(2)=vl-j
       wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
       wj(1)=1.d0-wj(2)+wj(3)
       wj(2)=wj(2)-2.d0*wj(3)
       ivmax=3
      else
       j=min(int(vl),10)    
       wj(2)=vl-j
       wj(1)=1.d0-wj(2)
       ivmax=2
      endif
      
      do j1=1,ivmax
       j2=j+j1-1
      do l1=1,izmax
       l2=jz+l1-1
      do k1=1,iymax
       k2=k+k1-1
       qglegi=qglegi+qlegi(k2,l2,j2,icz)*wk(k1)*wz(l1)*wj(j1)
      enddo
      enddo
      enddo
      qglegi=exp(qglegi)*z*cd(icdp,icz)
     **(1.d0-(1.d0-(1.d0-1.d0/sy)**(1.+ahl(icz)))**(1.+dels-alpd))
      
      if(debug.ge.4)write (moniou,202)qglegi
      
201   format(2x,'qglegi - interpolation of Pomeron leg eikonal:'
     */4x,'s=',e10.3,2x,'b^2=',e10.3,2x,'vvx=',e10.3
     *,2x,'icdp=',i1,2x,'icz=',i1,2x)
202   format(2x,'qglegi=',e10.3)
      return 
      end
     
c------------------------------------------------------------------------
      double precision function qgfan(sy,bb,vvx,icdp,icz)
c-----------------------------------------------------------------------
c qgfan - integrated uncut fan-contributions
c sy   - Pomeron mass squared,
c bb   - squared impact parameter,
c vvx  - relative strenth of screening corrections (0<vvx<1),
c icdp - diffractive eigenstate,
c icz  - hadron class
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr19/ ahl(3)
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /qgarr51/ x4(2),a4(2)
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)sy,bb,vvx,icdp,icz
      if(sy.le.1.d0)then                         !no kinematic space
       qgfan=0.d0
       if(debug.ge.4)write (moniou,202)qgfan
       return
      endif
      qgfan=qglegi(sy,bb,vvx,icdp,icz)
      
      dps=0.d0
      do ix1=1,7
      do mx1=1,2
       xpomr1=sy**(-(.5d0+x1(ix1)*(mx1-1.5d0)))  !-log(rap.) for multi-P vertex
	 	
       rp1=alfp*log(xpomr1*sy)*4.d0*.0389d0	
       do ix2=1,2  !7
       do mx2=1,2
        bb1=-rp1*log(.5d0+x4(ix2)*(mx2-1.5d0))   !impact parameter for unterm. P 
       do ix3=1,2  !7
       do mx3=1,2
        phi=2.d0*pi*(.5d0+x4(ix3)*(mx3-1.5d0))   !polar angle
	bb2=(dsqrt(bb)-dsqrt(bb1)*cos(phi))**2+bb1*sin(phi)**2  
                                            !impact parameter for multi-P vertex
        v1p=qgfani(1.d0/xpomr1,bb2,vvx,0.d0,0.d0,icdp,icz,1)   !fan contribution
        v1i=qgpini(xpomr1*sy,bb1,vvx)                          !intermediate P
        dps=dps+a1(ix1)*a4(ix2)*a4(ix3)*v1i*(1.d0-exp(-v1p)-v1p)
       enddo
       enddo
       enddo
       enddo
      enddo
      enddo
      qgfan=qgfan+dps*log(sy)/2.d0*pi*r3p/g3p*sigs*(1.d0-vvx)
      qgfan=qgfan/(1.d0-(1.d0-(1.d0-1.d0/sy)**(1.+ahl(icz)))  !divided by factor
     ***(1.+dels-alpd))                                       !for interpolation
      if(debug.ge.4)write (moniou,202)qgfan
     
201   format(2x,'qgfan - integrated uncut fan-contributions:'
     */4x,'s=',e10.3,2x,'b^2=',e10.3,2x,'vvx=',e10.3
     *,2x,'icdp=',i1,2x,'icz=',i1)
202   format(2x,'qgfan=',e10.3)
      return
      end  

c------------------------------------------------------------------------
      subroutine qgfan1(fann,sy,bb,vvx,vvxp,vvxt,icdp,icz)
c-----------------------------------------------------------------------
c qgfan1 - integrated cut fan-contributions
c sy    - Pomeron mass squared,
c bb    - squared impact parameter,
c vvx,vvxp,vvxt - screening corrections from targ. and nuclear proj. fans
c icdp  - diffractive eigenstate,
c icz   - hadron class,
c iqq=2 - 1-Pomeron contribution (screened); "holes" included,
c iqq=3 - 1-Pomeron contribution without "holes"
c iqq=4 - 1-Pomeron contribution without uncut end
c iqq=5 - total "fan" cut without uncut end,
c iqq=6 - total cut without uncut end
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension fann(6)
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr19/ ahl(3)
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /qgarr51/ x4(2),a4(2)
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)sy,bb,vvx,vvxp,vvxt,icdp,icz,iqq
      do iqq=2,6
       fann(iqq)=0.d0
      enddo      
      if(sy.le.1.d0)goto 1                                  !no kinematic space
      
      do ix1=1,7
      do mx1=1,2
       xpomr1=sy**(-(.5d0+x1(ix1)*(mx1-1.5d0)))  !-log(rap.) for multi-P vertex	 	
       rp1=alfp*log(xpomr1*sy)*4.d0*.0389d0	
       do ix2=1,2
       do mx2=1,2
        bb1=-rp1*log(.5d0+x4(ix2)*(mx2-1.5d0))   !impact parameter for unterm. P
       do ix3=1,2
       do mx3=1,2
        phi=2.d0*pi*(.5d0+x4(ix3)*(mx3-1.5d0))   !polar angle
	bb2=(dsqrt(bb)-dsqrt(bb1)*cos(phi))**2+bb1*sin(phi)**2
                                            !impact parameter for multi-P vertex       
        v1p=qgfani(1.d0/xpomr1,bb2,vvx,0.d0,0.d0,icdp,icz,1) !different cut fans
	v1p0=min(v1p,qgfani(1.d0/xpomr1,bb2,vvx,vvxp,vvxt,icdp,icz,5))
	v1pt=qgfani(1.d0/xpomr1,bb2,vvx,vvxp,vvxt,icdp,icz,6)
	v1pl=qgfani(1.d0/xpomr1,bb2,vvx,vvxp,0.d0,icdp,icz,2)
	v1p1=qgfani(1.d0/xpomr1,bb2,vvx,vvxp,0.d0,icdp,icz,3)
        v1pl0=min(v1pl,qgfani(1.d0/xpomr1,bb2,vvx,vvxp,vvxt,icdp,icz,4))
	
	do iqq=2,6
         if(iqq.eq.5)then
          vvxi=1.d0-(1.d0-vvx)*(1.d0-vvxt)
	  v1i=qgpini(xpomr1*sy,bb1,vvxi)
	  dpx=v1i*(1.d0-vvxi)*(1.d0-exp(-v1p)-v1p0
     *    -(1.d0-vvxp)*exp(-v1p)*(exp(v1p-v1p0)-1.d0))
         elseif(iqq.eq.6)then
          v1i=qgpini(xpomr1*sy,bb1,vvxt)
	  dpx=v1i*((1.d0-vvxt)*(1.d0-exp(-v1p)-v1pt)
     *    +(exp(-v1p)-exp(-v1pt))*(1.d0-vvxp))
         elseif(iqq.eq.2)then
          vvxi=1.d0-(1.d0-vvx)*(1.d0-vvxp)
	  v1i=qgpini(xpomr1*sy,bb1,vvxi)
	  dpx=v1i*(1.d0-vvxi)*v1pl*(exp(-2.d0*v1p)-1.d0)
         elseif(iqq.eq.3)then
          vvxi=1.d0-((1.d0-vvx)*(1.d0-vvxp))**2
	  v1i=qgpini(xpomr1*sy,bb1,vvxi)
	  dpx=v1i*(1.d0-vvxi)*v1p1*(exp(-2.d0*v1p)-1.d0)
         elseif(iqq.eq.4)then
          vvxi=1.d0-dsqrt((1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt)**3)
	  v1i=qgpini(xpomr1*sy,bb1,vvxi)
          dpx=v1i*((1.d0-vvxi)*v1pl0*(exp(-v1p)-1.d0)+v1pl*exp(-v1p)
     *    *((1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt)*exp(-v1p)-1.d0+vvxi))
         endif     
         fann(iqq)=fann(iqq)+a1(ix1)*a4(ix2)*a4(ix3)*dpx
	enddo
       enddo
       enddo
       enddo
       enddo
      enddo
      enddo
1     continue
      do iqq=2,6                                           !add leg contribution
       fann(iqq)=fann(iqq)*dlog(sy)*pi*r3p/g3p*sigs/2.d0
       if(iqq.eq.5)then
        vvxi=1.d0-(1.d0-vvx)*(1.d0-vvxt)
       elseif(iqq.eq.6)then
        vvxi=vvxt
       elseif(iqq.eq.2)then
        vvxi=1.d0-(1.d0-vvx)*(1.d0-vvxp)
       elseif(iqq.eq.3)then
        vvxi=1.d0-((1.d0-vvx)*(1.d0-vvxp))**2
       elseif(iqq.eq.4)then
        vvxi=1.d0-dsqrt((1.d0-vvx)*(1.d0-vvxp)*(1.d0-vvxt)**3)
       endif
       fann(iqq)=(fann(iqq)+qglegi(sy,bb,vvxi,icdp,icz))
     * /(1.d0-(1.d0-(1.d0-1.d0/sy)**(1.+ahl(icz)))**(1.+dels-alpd))
      enddo
      if(debug.ge.4)write (moniou,202)
     
201   format(2x,'qgfan1 - integrated cut fan-contributions:'
     */4x,'s=',e10.3,2x,'b^2=',e10.3,2x,'vvx=',e10.3,2x
     */4x,'vvxp=',e10.3,'vvxt=',e10.3,2x,'icdp=',i1
     *,2x,'icz=',i1,2x,'iqq=',i1)
202   format(2x,'qgfan1-end')
      return
      end  

c------------------------------------------------------------------------
      double precision function qgfani(sy,bb,vvx,vvxp,vvxt,icdp,icz,iqq)
c-----------------------------------------------------------------------
c qgfani - integrated fan-contributions (interpolation)
c sy    - Pomeron mass squared,
c bb    - squared impact parameter,
c vvx,vvxp,vvxt - screening corrections from targ. and nuclear proj. fans
c icdp  - diffractive eigenstate,
c icz   - hadron class,
c iqq=1 - total cut,
c iqq=2 - 1-Pomeron contribution (screened); "holes" included,
c iqq=3 - 1-Pomeron contribution without "holes"
c iqq=4 - 1-Pomeron contribution without uncut end
c iqq=5 - total "fan" cut without uncut end,
c iqq=6 - total cut without uncut end
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wz(3),wj(3),wi(3),wn(3)
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr19/ ahl(3)
      common /qgarr20/ spmax
      common /qgarr27/ qlegi(51,11,11,3),qfanu(51,11,11,6)
     *,qfan(51,11,11,6,120)
      common /qgarr43/ moniou
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)sy,bb,vvx,vvxp,vvxt,icdp,icz,iqq  !so161205
      qgfani=0.d0
      if(sy.le.1.d0)then                          !no kinematic space
       if(debug.ge.4)write (moniou,202)qgfani
       return
      endif

      rp=(rq(icz)+alfp*dlog(max(1.d0,sy)))*4.d0*.0389d0
      z=dexp(-bb/rp)
      if(z.lt..2d0*dexp(-4.d0))then
       izmax=2
       jz=1
       wz(2)=5.d0*z*dexp(4.d0)
       wz(1)=1.d0-wz(2)
      else
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=(-bb/rp-dlog(0.2d0))/.8d0+7.d0
       endif
       jz=min(10,int(zz))
       jz=max(2,jz)
       wz(2)=zz-dble(jz)
       wz(1)=1.d0-wz(2)
       izmax=2
      endif

      yl=dlog(sy)/dlog(spmax)*50.d0+1.d0
      k=max(1,int(1.00001d0*yl-1.d0))
      k=min(k,49)
      if(k.eq.1)then
       wk(2)=yl-dble(k)
       wk(1)=1.d0-wk(2)
       iymax=2
      else
       wk(2)=yl-dble(k)
       wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
       wk(1)=1.d0-wk(2)+wk(3)
       wk(2)=wk(2)-2.d0*wk(3)
       iymax=3
      endif

      vl=max(1.d0,vvx*10.d0+1.d0)
      if(vvx.eq.0.d0)then
       j=1
       ivmax=1
       wj(1)=1.d0
      elseif(vl.lt.2.d0)then
       j=1    
       wj(2)=vl-dble(j)
       wj(1)=1.d0-wj(2)
       ivmax=2
      else
       j=min(int(vl),10)    
       wj(2)=vl-dble(j)
       wj(1)=1.d0-wj(2)
       ivmax=2
      endif
      
      if(iqq.eq.1)then                               !uncut fan
       do j1=1,ivmax
        j2=j+j1-1
       do l1=1,izmax
        l2=jz+l1-1
       do k1=1,iymax
        k2=k+k1-1
        qgfani=qgfani+qfanu(k2,l2,j2,icdp+2*(icz-1))
     *  *wk(k1)*wz(l1)*wj(j1)
       enddo
       enddo
       enddo
      else
       if(vvx.eq.1.d0.and.iqq.ne.6)then              !full screening
        ivmax=1
        iv1max=1
        iv2max=1
	i=1
	n=1
	j=11
        wi(1)=1.d0
        wj(1)=1.d0
        wn(1)=1.d0
       elseif(vvxp.eq.0.d0.and.vvxt.eq.0.d0)then     !hadron (no nuclear graphs)
        iv1max=1
        iv2max=1
	i=1
	n=1
        wi(1)=1.d0
        wn(1)=1.d0
       else         !nuclear effects (Pomerons may couple to different nucleons)
        iv1max=2
        vl1=max(1.d0,vvxp*5.d0+1.d0)
        i=min(int(vl1),5)    
        wi(2)=vl1-i
        wi(1)=1.d0-wi(2)
	
        if(vvx.lt..01d0)then
         iv2max=1
	 n=1
         wn(1)=1.d0
	else
	 iv2max=2
         vl2=max(1.d0,vvxt/vvx*5.d0+1.d0)
         n=min(int(vl2),5)    
         wn(2)=vl2-n
         wn(1)=1.d0-wn(2)
	endif
       endif

       if(iqq.le.3)then
        do i1=1,iv1max
        do j1=1,ivmax
         j2=j+j1-1
        do l1=1,izmax
         l2=jz+l1-1
        do k1=1,iymax
         k2=k+k1-1
         qgfani=qgfani+qfan(k2,l2,j2,i+i1-1,icdp+2*(icz-1)+6*(iqq-2))
     *   *wk(k1)*wz(l1)*wj(j1)*wi(i1)
        enddo
        enddo
        enddo
        enddo
       else
        do n1=1,iv2max
        do i1=1,iv1max
        do j1=1,ivmax
         j2=j+j1-1
        do l1=1,izmax
         l2=jz+l1-1
        do k1=1,iymax
         k2=k+k1-1
         qgfani=qgfani+qfan(k2,l2,j2,i+i1-1,n+n1+11+6*(icdp-1
     *   +2*(icz-1)+6*(iqq-4)))*wk(k1)*wz(l1)*wj(j1)*wi(i1)*wn(n1)
        enddo
        enddo
        enddo
        enddo
        enddo
       endif
      endif
      qgfani=dexp(qgfani)*(1.d0-(1.d0-(1.d0-1.d0/sy)**(1.d0+ahl(icz)))
     ***(1.d0+dels-alpd))*z
      if(debug.ge.4)write (moniou,202)qgfani

201   format(2x,'qgfani - integrated fan-contributions (interpolation):'
     */4x,'s=',e10.3,2x,'b^2=',e10.3,2x,'vvx=',e10.3,2x      !so161205
     */4x,'vvxp=',e10.3,'vvxt=',e10.3                        !so161205
     *,2x,'icdp=',i1,2x,'icz=',i1,2x,'iqq=',i1)
202   format(2x,'qgfani=',e10.3)
      return 
      end
                                                
c------------------------------------------------------------------------
      double precision function qgppd(xph,vvx,iqq)
c-----------------------------------------------------------------------
c qgppd - parton distributions in the Pomeron (with screening)
c xph   - parton LC momentum share,
c vvx   - relative strenth of screening corrections,
c iqq=0 - gluon,
c iqq=1 - sea quark
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)xph,vvx,iqq
      qgppd=0.d0
      alf=4.d0*pi*r3p/g3p*sigs*vvx
      if(alf.gt..02d0)xpmin=xph**alf
      do i=1,7
      do m=1,2
       if(alf.gt.0.02d0)then
        xpomr=(.5d0*(1.d0+xpmin-x1(i)*(2*m-3)*(1.d0-xpmin)))**(1.d0/alf)
       else
        xpomr=xph**(.5d0+x1(i)*(m-1.5d0))
       endif
       if(1.d0-xph/xpomr.gt.1.d-5)then             !otherwise =0
        if(alf.gt.0.02d0)then                      !weak screening
	 if(iqq.eq.0)then                          !gluon
          qgppd=qgppd+a1(i)*(1.d0-xph/xpomr)**betp*(1.d0-xpmin)
	 else                                      !quark
          qgppd=qgppd+a1(i)*qgftlf(xph/xpomr)*(1.d0-xpmin)
	 endif
        else                                       !strong screening
	 if(iqq.eq.0)then                          !gluon
          qgppd=qgppd-a1(i)*(1.d0-xph/xpomr)**betp*xpomr**alf
     *    *alf*log(xph)
 	 else                                      !quark
          qgppd=qgppd-a1(i)*qgftlf(xph/xpomr)*xpomr**alf*alf*log(xph)
         endif
        endif
       endif
      enddo
      enddo
      if(iqq.eq.0)then                             !gluon
       qgppd=max(((1.d0-xph)**betp-qgppd/2.d0)*(1.d0-dgqq)
     * ,(1.d0-xph)**betp*(1.d0-dgqq)*xph**alf)
      else                                         !quark
       qgppd=max((qgftlf(xph)-qgppd/2.d0)*dgqq
     * ,qgftlf(xph)*dgqq*xph**alf)
      endif
      if(debug.ge.4)write (moniou,202)qgppd
     
201   format(2x,'qgppd - parton distributions in the Pomeron:'
     */4x,'xph=',e10.3,2x,'vvx=',e10.3,2x,'iqq=',i1)
202   format(2x,'qgppd=',e10.3)
      return
      end
     
c------------------------------------------------------------------------
      double precision function qgppdi(xp,vvx,iqq)
c-----------------------------------------------------------------------
c qgppdi - parton distributions in the Pomeron (interpolation)
c xp    - parton LC momentum share,
c vvx   - relative strenth of screening corrections,
c iqq=0 - gluon
c iqq=1 - sea quark
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wi(3),wj(3)
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr20/ spmax
      common /qgarr43/ moniou
      common /qgarr45/ qppdi(10,11,2)
      common /debug/   debug
       
      if(debug.ge.3)write (moniou,201)xp,vvx,iqq
      if(xp.lt..2d0)then
       xl=5.d0*log(5.d0*xp)/log(.2d0*spmax)+6.d0
      else
       xl=5.d0*xp+5.d0
      endif
      i=min(8,int(xl))
      i=max(1,i)
      if(i.eq.5)i=4
      wi(2)=xl-i
      wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
      wi(1)=1.d0-wi(2)+wi(3)
      wi(2)=wi(2)-2.d0*wi(3)
      
      vl=max(1.d0,vvx*10.d0+1.d0)
      j=min(int(vl),9)    
      wj(2)=vl-j
      wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
      wj(1)=1.d0-wj(2)+wj(3)
      wj(2)=wj(2)-2.d0*wj(3)

      qgppdi=0.d0
      do j1=1,3
       j2=j+j1-1
      do i1=1,3
       i2=i+i1-1
       qgppdi=qgppdi+qppdi(i2,j2,iqq+1)*wi(i1)*wj(j1)
      enddo
      enddo
      qgppdi=exp(qgppdi)
      if(iqq.eq.0)then                             !gluon
       qgppdi=qgppdi*(1.d0-xp)**betp*(1.d0-dgqq)
      elseif(iqq.eq.1)then                         !quark
       qgppdi=qgppdi*qgftlf(xp)*dgqq
      endif
      if(debug.ge.4)write (moniou,202)qgppdi
     
201   format(2x,'qgppdi - parton distr. in the Pomeron (interpol.):'
     */4x,'xp=',e10.3,2x,'vvx=',e10.3,2x,'iqq=',i1)
202   format(2x,'qgppdi=',e10.3)
      return 
      end      
      
c=============================================================================
      double precision function qgftld(z,icz)
c-----------------------------------------------------------------------------
c qgftld - auxilliary function for calculations of semihard eikonals
c (proportional to gluon pdf: g(z)*z^(1+dels))
c z   - LC x,
c icz - hadron class
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr19/ ahl(3)
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)z,icz
      qgftld=0.d0
      xpmin=z**(1.d0-alpd+dels)
      do i1=1,7
      do m1=1,2
       tp=1.d0-(1.d0-xpmin)*(.5d0+x1(i1)*(m1-1.5d0))
     * **(1.d0/(1.d0+ahl(icz)))
       xp=tp**(1.d0/((1.d0-alpd+dels)))
       qgftld=qgftld+a1(i1)*((1.d0-xp)/(1.d0-tp))**ahl(icz)
     * *(1.d0-z/xp)**betp
      enddo
      enddo
      qgftld=qgftld*.5d0*(1.d0-xpmin)**(ahl(icz)+1.d0)
     */(ahl(icz)+1.d0)/(1.d0-alpd+dels)*(1.d0-dgqq)
      if(debug.ge.4)write (moniou,202)qgftld
      
201   format(2x,'qgftld:',2x,'z=',e10.3,2x,'icz=',i1)
202   format(2x,'qgftld=',e10.3)
      return
      end

c------------------------------------------------------------------------
      double precision function qgftle(z,icz)
c-----------------------------------------------------------------------
c qgftle - auxilliary function calculations of semihard eikonals
c (proportional to sea quark pdf: q_s(z)*z^(1+dels))
c z   - LC x,
c icz - hadron class
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr19/ ahl(3)
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)z,icz
      qgftle=0.d0
      xpmin=z**(1.d0-alpd+dels)
      do i1=1,7
      do m1=1,2
       tp=1.d0-(1.d0-xpmin)*(.5d0+x1(i1)*(m1-1.5d0))
     * **(1.d0/(1.d0+ahl(icz)))
       xp=tp**(1.d0/((1.d0-alpd+dels)))
       qgftle=qgftle+a1(i1)*((1.d0-xp)/(1.d0-tp))**ahl(icz)
     * *qgftlf(z/xp)
      enddo
      enddo
      qgftle=qgftle*.5d0*(1.d0-xpmin)**(ahl(icz)+1.d0)
     */(ahl(icz)+1.d0)/(1.d0-alpd+dels)*dgqq
      if(debug.ge.4)write (moniou,202)qgftle
      
201   format(2x,'qgftle:',2x,'z=',e10.3,2x,'icz=',i1)
202   format(2x,'qgftle=',e10.3)
      return
      end

c------------------------------------------------------------------------
      double precision function qgftlf(zz)
c-----------------------------------------------------------------------
c qgftlf - auxilliary function for calculations of semihard eikonals
c (qgftlf=int(dz) z^dels * (1-zz/z)^betp * P_qG(z))
c zz - ratio of the quark and Pomeron light cone x (zz=x_G/x_P)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)zz
      qgftlf=0.d0
      zmin=zz**(1.d0+dels)
      do i=1,7
      do m=1,2
       z=(.5d0*(1.d0+zmin+(2*m-3)*x1(i)*(1.d0-zmin)))
     * **(1.d0/(1.d0+dels))
       qgftlf=qgftlf+a1(i)*max(1.d-9,(1.d0-zz/z))**betp
     * *(z**2+(1.d0-z)**2)
      enddo
      enddo
      qgftlf=qgftlf*1.5d0*(1.d0-zmin)/(1.d0+dels)   !1.5=n_flav/2 
      if(debug.ge.4)write (moniou,202)qgftlf

201   format(2x,'qgftlf:',2x,'zz=',e10.3)
202   format(2x,'qgftlf=',e10.3)
      return
      end

c=============================================================================
      subroutine qgfz(b,gz,iddp1,iddp2)
c----------------------------------------------------------------------------
c qgfz - hadron-hadron and hadron-nucleus cross sections
c b            - impact parameter (in case of hadron-nucleus),
c iddp1, iddp2 - diffractive eigenstates (in case of hadron-nucleus),
c gz(i), i=1,5 - contributions of different cuts (elastic, inel., and diffr.)
c----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207)
      dimension gz(5),wt1(3),wt2(3)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug
Cf2py intent(in) b,iddp1,iddp2
Cf2py intent(out) gz

      if(debug.ge.3)write (moniou,201)b,iddp1,iddp2
      do l=1,5
       gz(l)=0.d0
      enddo
      rp=(rq(icz)+rq(2)+alfp*log(scm))*4.d0*.0389d0
      if(ia(2).eq.1.and.iddp1.eq.0.and.iddp2.eq.0)then
       g0=pi*rp*10.d0                     !normalization factor (in mb)
       bm=2.d0*dsqrt(rp)                  !impact parameter for exp. fall-down
      endif

      do i1=1,7
      do m=1,2
       z=.5d0+x1(i1)*(m-1.5d0)
       bb1=rp*z
       bb2=rp*(1.d0-dlog(z))
       
       do l=1,3
        wt1(l)=0.d0
        wt2(l)=0.d0
       enddo
       if(ia(2).eq.1)then                 !hadron-proton
        do idd1=1,2
        do idd2=1,2
         vv1=exp(-qgpomi(scm,bb1,0.d0,0.d0,0.d0,0.d0,idd1,idd2,icz,1))
	 vv2=exp(-qgpomi(scm,bb2,0.d0,0.d0,0.d0,0.d0,idd1,idd2,icz,1))
                                          !eikonals for given diffr. eigenstates
         do l=1,2
          wt1(l)=wt1(l)+cc(idd1,icz)*cc(idd2,2)*vv1**l
          wt2(l)=wt2(l)+cc(idd1,icz)*cc(idd2,2)*vv2**l
         enddo
         do idd3=1,2
          wt1(3)=wt1(3)+cc(idd1,icz)*cc(idd2,2)*cc(idd3,icz)*vv1
     *    *exp(-qgpomi(scm,bb1,0.d0,0.d0,0.d0,0.d0,idd3,idd2,icz,1))
          wt2(3)=wt2(3)+cc(idd1,icz)*cc(idd2,2)*cc(idd3,icz)*vv2
     *    *exp(-qgpomi(scm,bb2,0.d0,0.d0,0.d0,0.d0,idd3,idd2,icz,1))
         enddo
        enddo
        enddo
        do l=1,2                        !gz(i) - contributions of different cuts
         gz(l)=gz(l)+a1(i1)*((1.d0-wt1(l))+(1.d0-wt2(l))/z)
        enddo
        gz(3)=gz(3)+a1(i1)*((wt1(2)-wt1(3))+(wt2(2)-wt2(3))/z)
        gz(4)=gz(4)+a1(i1)*((wt1(3)-wt1(1)**2)+(wt2(3)-wt2(1)**2)/z)
        gz(5)=gz(5)+a1(i1)*((1.d0-wt1(1))*z+(1.d0-wt2(1))/z*(1.-log(z)))

       else                              !hadron-nucleus
        do idd1=1,2
        do idd2=1,2
         vv1=exp(-qgpomi(scm,bb1,0.d0,0.d0,0.d0,0.d0,iddp1,idd1,icz,1)
     *   -qgpomi(scm,bb1,0.d0,0.d0,0.d0,0.d0,iddp2,idd2,icz,1))
         vv2=exp(-qgpomi(scm,bb2,0.d0,0.d0,0.d0,0.d0,iddp1,idd1,icz,1)
     *   -qgpomi(scm,bb2,0.d0,0.d0,0.d0,0.d0,iddp2,idd2,icz,1))

         if(idd1.eq.idd2)then
          wt1(1)=wt1(1)+cc(idd1,2)*vv1
          wt2(1)=wt2(1)+cc(idd1,2)*vv2
         endif
         wt1(2)=wt1(2)+cc(idd1,2)*cc(idd2,2)*vv1
         wt2(2)=wt2(2)+cc(idd1,2)*cc(idd2,2)*vv2
        enddo
        enddo
        cg1=qgrot(b,dsqrt(bb1))          !convolution with nuclear profile
        cg2=qgrot(b,dsqrt(bb2))
        do l=1,2
         gz(l)=gz(l)+a1(i1)*(cg1*(1.d0-wt1(l))+cg2*(1.d0-wt2(l))/z)
        enddo
       endif
      enddo
      enddo
      if(ia(2).eq.1.and.iddp1.eq.0.and.iddp2.eq.0)then     !hadron-proton
       do l=1,5
        gz(l)=gz(l)*g0
       enddo
       gz(5)=gz(5)/gz(1)*(rq(icz)+rq(2)+alfp*log(scm))*2.d0
      endif
      if(debug.ge.3)write (moniou,202)gz
      if(debug.ge.4)write (moniou,203)
      
201   format(2x,'qgfz - hadronic cross-sections calculation'
     */4x,'b=',e10.3,2x,'iddp=',2i3)
202   format(2x,'qgfz: gz=',5e10.3)      
203   format(2x,'qgfz - end')
      return
      end

c=============================================================================
      double precision function qghard(sy,bb,icz)
c-----------------------------------------------------------------------------
c qghard - hard quark-quark eikonal
c s   - energy squared for the interaction (hadron-hadron),
c bb  - squared impact parameter,
c icz - type of the primaty hadron (nucleon)
c----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr25/ ahv(3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)sy,icz
      qghard=0.d0
      s2min=4.d0*fqscal*qt0                   !energy threshold for hard process
      xmin=s2min/sy
      if(xmin.ge.1.d0)return
      xmin=xmin**(delh+.5d0)

      do i=1,7
      do m=1,2
       z1=(.5d0*(1.d0+xmin-(2*m-3)*x1(i)*(1.d0-xmin)))
     * **(1.d0/(delh+.5d0))

       st2=0.d0
       do j=1,7
       do k=1,2
        xx=.5d0*(1.d0+x1(j)*(2*k-3))
	xp=z1**xx
	xm=z1/xp
        st2=st2+a1(j)*(1.d0-xp)**ahv(icz)*(1.d0-xm)**ahv(2)
     *  *(qggrv(xp,qt0,icz,1)+qggrv(xp,qt0,icz,2))
     *  *(qggrv(xm,qt0,2,1)+qggrv(xm,qt0,2,2))/dsqrt(z1)
       enddo
       enddo
       sj=qgjit(qt0,qt0,z1*sy,2,2)
       st2=-st2*dlog(z1)*sj
       if(debug.ge.5)write (moniou,202)z1*sy,sj

       qghard=qghard+a1(i)/z1**delh*st2
      enddo
      enddo
      qghard=qghard*(1.d0-xmin)/(.5d0+delh)*.25d0*factk
      qghard=qghard/(8.d0*pi*(rq(icz)+rq(2)))        !gaussian distrobition in b
     **exp(-bb/(4.d0*.0389d0*(rq(icz)+rq(2))))      
      if(debug.ge.4)write (moniou,203)qghard
      
201   format(2x,'qghard - hard quark-quark interaction eikonal:'
     */2x,'s=',e10.3,2x,'icz=',i1)
202   format(2x,'qghard:',2x,'s_hard=',e10.3,2x,'sigma_hard=',e10.3)
203   format(2x,'qghard=',e10.3)
      return
      end

c=============================================================================
      double precision function qg3pom(sy,b,icdp,icdt,icz)
c-----------------------------------------------------------------------
c qg3pom - integrated uncut multi-Pomeron contributions 
c to the interaction eikonal
c sy   - Pomeron mass squared,
c b    - impact parameter,
c icdp - projectile diffractive eigenstate,
c icdp - target diffractive eigenstate,
c icz  - projectile class
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /qgarr53/ x9(3),a9(3)
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)sy,b,icdp,icdt,icz
      qg3pom=0.d0      
      if(sy.le.1.d0)return                        !below threshold

      rp=(rq(icz)+rq(2)+alfp*log(sy))*4.d0*.0389d0
      do ib1=1,7
      do mb1=1,2
       z=.5+x1(ib1)*(mb1-1.5)
       b1=sqrt(-rp/4.*log(z))
       dbb=b1**2+.25*b**2
      do ib2=1,7
      do mb2=1,2
       phi=pi*(.5+x1(ib2)*(mb2-1.5))
       bb1=dbb+b*b1*cos(phi)
       bb2=2.*dbb-bb1  
      do ix1=1,7
      do mx1=1,2
       xpomr1=sy**(-(.5+x1(ix1)*(mx1-1.5)))
       	
       v1p1=qglegi(1.d0/xpomr1,bb1,0.d0,icdp,icz)          !unscreened proj. leg
       v1t1=qglegi(xpomr1*sy,bb2,0.d0,icdt,2)              !unscreened targ. leg
       v1pf=min(v1p1,qgfani(1.d0/xpomr1,bb1,0.d0,0.d0,0.d0,icdp,icz,1))
                                                           !proj. fan
       v1tf=min(v1t1,qgfani(xpomr1*sy,bb2,0.d0,0.d0,0.d0,icdt,2,1))
                                                           !targ. fan
       v1pt=min(v1pf,qgfani(1.d0/xpomr1,bb1,v1tf,0.d0,0.d0,icdp,icz,1))
                                                           !proj. gener. fan
       v1tt=min(v1tf,qgfani(xpomr1*sy,bb2,v1pf,0.d0,0.d0,icdt,2,1))
                                                           !targ. gener. fan
       v1pr=v1pt-v1pf                                      !proj. reversed fan
       v1tr=v1tt-v1tf                                      !targ. reversed fan
       
       dpx=min(0.d0,1.d0-exp(-v1pt)-v1pt)*min(0.d0,1.d0-exp(-v1tt)-v1tt)
     * +v1pt*min(0.d0,1.d0-exp(-v1tt)-v1tt)
     * +v1tt*min(0.d0,1.d0-exp(-v1pt)-v1pt)
     * -.5d0*v1pr*min(0.d0,(1.d0-exp(-v1tt))*exp(-v1pf)-v1tt)
     * -.5d0*(v1pt-v1p1)*min(0.d0,1.d0-exp(-v1tf)-v1tf)
     * -.5d0*v1tr*min(0.d0,(1.d0-exp(-v1pt))*exp(-v1tf)-v1pt)
     * -.5d0*(v1tt-v1t1)*min(0.d0,1.d0-exp(-v1pf)-v1pf)
       
       qg3pom=qg3pom+a1(ib1)*a1(ib2)*a1(ix1)/z*dpx
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      qg3pom=qg3pom*rp/32.d0*log(sy)*(r3p*pi/.0389d0)/g3p**3        
      if(debug.ge.3)write (moniou,202)qg3pom
      
201   format(2x,'qg3pom - multi-Pomeron contributions to the eikonal:'
     */2x,'sy=',e10.3,2x,'b=',e10.3,2x,'icdp=',i1,2x,'icdt=',i1
     *,2x,'icz=',i1)
202   format(2x,'qg3pom=',e10.3)
      return 
      end
      
c=============================================================================
      double precision function qg3pa(sy,b,vvxp,vvxt,vvxpa,vvxta
     *,icdp,icdt,icz)
c-----------------------------------------------------------------------
c qg3pa - integrated cut multi-Pomeron contributions to the interaction eikonal
c sy          - Pomeron mass squared,
c b           - impact parameter,
c vvxp, vvxpa - proj. nuclear screening (cut and uncut legs to other nucleons),
c vvxt, vvxta - targ. nuclear screening,
c icdp        - projectile diffractive eigenstate,
c icdp        - target diffractive eigenstate,
c icz         - projectile class
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /qgarr51/ x4(2),a4(2)
      common /qgarr53/ x9(3),a9(3)
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)sy,b,vvxp,vvxt,vvxpa,vvxta
     *,icdp,icdt,icz
      qg3pa=0.d0      
      if(sy.le.1.d0)return

      rp=(rq(icz)+rq(2)+alfp*log(sy))*4.d0*.0389d0
      do ib1=1,2
      do mb1=1,2
       z=.5+x4(ib1)*(mb1-1.5)
       b1=sqrt(-rp/4.*log(z))
       dbb=b1**2+.25*b**2
      do ib2=1,2
      do mb2=1,2
       phi=pi*(.5+x4(ib2)*(mb2-1.5))
       bb1=dbb+b*b1*cos(phi)
       bb2=2.*dbb-bb1 
      do ix1=1,7
      do mx1=1,2
       xpomr1=sy**(-(.5+x1(ix1)*(mx1-1.5)))
       
       vpau0=qgfani(1.d0/xpomr1,bb1,1.d0-(1.d0-vvxp)*(1.d0-vvxt)
     * *(1.d0-vvxta),0.d0,0.d0,icdp,icz,1)
       vtau0=qgfani(xpomr1*sy,bb2,1.d0-(1.d0-vvxt)*(1.d0-vvxp)
     * *(1.d0-vvxpa),0.d0,0.d0,icdt,2,1)
     
       nn=0
1      nn=nn+1
       vpac=qgfani(1.d0/xpomr1,bb1,1.d0-exp(-vtau0)*(1.d0-vvxp)
     * *(1.d0-vvxt)*(1.d0-vvxta),0.d0,0.d0,icdp,icz,1)
       vtac=qgfani(xpomr1*sy,bb2,1.d0-exp(-vpau0)*(1.d0-vvxt)
     * *(1.d0-vvxp)*(1.d0-vvxpa),0.d0,0.d0,icdt,2,1)
       if((abs(vpac-vpau0).gt.1.d-2.or.abs(vtac-vtau0).gt.1.d-2)
     * .and.nn.lt.100)then
        vpau0=vpac
	vtau0=vtac
	goto 1
       endif
     
       v1p=min(vpac,qgfani(1.d0/xpomr1,bb1
     * ,1.d0-exp(-vtac)*(1.d0-vvxt)*(1.d0-vvxta)*(1.d0-vvxp)
     * ,1.d0-(1.d0-vvxp)*(1.d0-vvxpa)**2,0.d0,icdp,icz,2))
       v1t=min(vtac,qgfani(xpomr1*sy,bb2
     * ,1.d0-exp(-vpac)*(1.d0-vvxp)*(1.d0-vvxpa)*(1.d0-vvxt)
     * ,1.d0-(1.d0-vvxt)*(1.d0-vvxta)**2,0.d0,icdt,2,2))

       dpx=max(0.d0,1.d0-exp(-vpac)-vpac*exp(-vpac)*(1.d0-vvxpa))
     * *max(0.d0,1.d0-exp(-vtac)-vtac*exp(-vtac)*(1.d0-vvxta))-.25d0
     * *((1.d0-exp(-vpac))**2*(1.d0-vvxp)+2.d0*(1.d0-exp(-vpac))*vvxp)
     * *((1.d0-exp(-vtac))**2*(1.d0-vvxt)+2.d0*(1.d0-exp(-vtac))*vvxt)
     * +.5d0*(vtac+v1t)*exp(-vtac)*(1.d0-vvxta)
     * *max(0.d0,1.d0-exp(-vpac)-vpac*exp(-vpac)*(1.d0-vvxpa))
     * +.5d0*(vpac+v1p)*exp(-vpac)*(1.d0-vvxpa)
     * *max(0.d0,1.d0-exp(-vtac)-vtac*exp(-vtac)*(1.d0-vvxta))
     * +.5d0*exp(-vpac-vtac)*(1.d0-vvxpa)*(1.d0-vvxta)
     * *(v1p*(1.d0-exp(-vpac)*(1.d0-vvxp)*(1.d0-vvxpa))
     * *(vtac+v1t*exp(-vtac)*(1.d0-vvxt)*(1.d0-vvxta))
     * +v1t*(1.d0-exp(-vtac)*(1.d0-vvxt)*(1.d0-vvxta))
     * *(vpac+v1p*exp(-vpac)*(1.d0-vvxp)*(1.d0-vvxpa)))

       qg3pa=qg3pa+a4(ib1)*a4(ib2)*a1(ix1)/z*max(0.d0,min(1.d0,dpx))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      qg3pa=qg3pa*rp/32.d0*log(sy)*(r3p*pi/.0389d0)/g3p**3  

      if(debug.ge.3)write (moniou,202)qg3pa
      
201   format(2x,'qg3pa - cut multi-Pomeron eikonal contributions:'
     */2x,'sy=',e10.3,2x,'b=',e10.3,2x,'vvxp=',e10.3,2x,'vvxt=',e10.3
     */2x,'vvxpa=',e10.3,2x,'vvxta=',e10.3,2x,'icdp=',i1,2x,'icdt=',i1
     *,2x,'icz=',i1)
202   format(2x,'qg3pa=',e10.3)
      return 
      end
      
c=============================================================================
      subroutine qgbdef(bba,bbb,xxa,yya,xxb,yyb,xxp,yyp,jb)
c-----------------------------------------------------------------------
c qgbdef - defines coordinates (xxp,yyp) of a multi-Pomeron vertex
c bba - squared distance between the vertex and the projectile,
c bbb - squared distance between the vertex and the target,
c xxa, yya - projectile coordinates,
c xxb, yyb - target coordinates,
c jb=1,2 - mirror choice (left- or right-positioned)
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr43/ moniou
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)bba,bbb,xxa,yya,xxb,yyb,jb

      xx=xxa-xxb
      yy=yya-yyb
      bb=xx**2+yy**2
      if(bb.lt.1.d-3)then
       xxp=xxb+dsqrt(bba)
       yyp=yyb
      elseif(abs(yy).lt.1.d-6)then
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
      if(debug.ge.4)write (moniou,202)xxp,yyp
      if(debug.ge.4)write (moniou,203)
      
201   format(2x,'qgbdef - coordinates of a multi-Pomeron vertex:'
     */2x,'bba=',e10.3,2x,'bbb=',e10.3,2x,'xxa=',e10.3,2x,'yya=',e10.3
     */2x,'xxb=',e10.3,2x,'yyb=',e10.3,2x,'jb=',i1)
202   format(2x,'qgbdef: xxp=',e10.3,2x,'yyp=',e10.3)      
203   format(2x,'qgbdef - end')
      return
      end

c=============================================================================
      subroutine qgfdf(xxp,yyp,xpomr,vpac,vtac,sumup,sumut,sumcp,sumct
     *,ip,it,imp,imt)
c-----------------------------------------------------------------------
c qgfdf - configuration of fan contributions (cut and uncut fans)
c xxp, yyp -  coordinates of the multi-Pomeron vertex,
c xpomr    - LC momentum share of the multi-Pomeron vertex,
c ip       - proj. index,
c it       - targ. index
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207)
cdh   dimension vpau0(iapmax),vtau0(iapmax),vpau(iapmax),vtau(iapmax)
      dimension vpau0(iapmax),vtau0(iapmax)
     *,vpac(iapmax),vtac(iapmax)
cdh  *,vpac(iapmax),vtac(iapmax),vpac0(iapmax),vtac0(iapmax)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),b
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr43/ moniou
      common /qgarr46/ iconab(iapmax,iapmax),icona(iapmax)
     *,iconb(iapmax)
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)xxp,yyp,xpomr,ip,it,imp,imt

      sumup0=0.d0                      !proj. fans without targ. screening
      do ipp=1,ia(1)
       if(ipp.ne.imp)then
        if(iconab(ipp,it).eq.0)then    !no connection 
	                               !(nucleon too far from the vertex)
         vpau0(ipp)=0.d0
	else
         bbp=(xa(ipp,1)+b-xxp)**2+(xa(ipp,2)-yyp)**2
         vpau0(ipp)=qgfani(1.d0/xpomr,bbp,1.d0-exp(-sumup0)
     *   ,0.d0,0.d0,iddp(ipp),icz,1)
         sumup0=sumup0+vpau0(ipp)
	endif
       endif
      enddo
      if(iconab(imp,it).eq.0)then
       vpau0(imp)=0.d0
      else
       bbp=(xa(imp,1)+b-xxp)**2+(xa(imp,2)-yyp)**2
       vpau0(imp)=qgfani(1.d0/xpomr,bbp,1.d0-exp(-sumup0)
     * ,0.d0,0.d0,iddp(imp),icz,1)
       sumup0=sumup0+vpau0(imp)
      endif
             	
      sumut0=0.d0                      !targ. fans without proj. screening
      do itt=1,ia(2)
       if(itt.ne.imt)then
        if(iconab(ip,itt).eq.0)then    !no connection
         vtau0(itt)=0.d0
	else
         bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
         vtau0(itt)=qgfani(xpomr*scm,bbt,1.d0-exp(-sumut0)
     *   ,0.d0,0.d0,iddt(itt),2,1)
         sumut0=sumut0+vtau0(itt)
	endif
       endif
      enddo
      if(iconab(ip,imt).eq.0)then
       vtau0(imt)=0.d0
      else
       bbt=(xb(imt,1)-xxp)**2+(xb(imt,2)-yyp)**2
       vtau0(imt)=qgfani(xpomr*scm,bbt,1.d0-exp(-sumut0)
     * ,0.d0,0.d0,iddt(imt),2,1)
       sumut0=sumut0+vtau0(imt)
      endif
             	
      nn=0
1     nn=nn+1 
      sumup=0.d0                       !proj. fans with targ. screening
      do ipp=1,ia(1)
       if(ipp.ne.imp)then
        if(iconab(ipp,it).eq.0)then    !no connection
         vpac(ipp)=0.d0
	else
         bbp=(xa(ipp,1)+b-xxp)**2+(xa(ipp,2)-yyp)**2
         vpac(ipp)=qgfani(1.d0/xpomr,bbp,1.d0-exp(-sumup-sumut0)
     *   ,0.d0,0.d0,iddp(ipp),icz,1)
         sumup=sumup+vpac(ipp)
	endif
       endif
      enddo
      if(iconab(imp,it).eq.0)then
       vpac(imp)=0.d0
      else
       bbp=(xa(imp,1)+b-xxp)**2+(xa(imp,2)-yyp)**2
       vpac(imp)=qgfani(1.d0/xpomr,bbp,1.d0-exp(-sumup-sumut0)
     * ,0.d0,0.d0,iddp(imp),icz,1)
       sumup=sumup+vpac(imp)
      endif
             	
      sumut=0.d0                      !targ. uncut fans with proj. screening
      do itt=1,ia(2)
       if(itt.ne.imt)then
        if(iconab(ip,itt).eq.0)then
         vtac(itt)=0.d0
	else
         bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
         vtac(itt)=qgfani(xpomr*scm,bbt,1.d0-exp(-sumut-sumup0)
     *   ,0.d0,0.d0,iddt(itt),2,1)
         sumut=sumut+vtac(itt)
	endif
       endif
      enddo
      if(iconab(ip,imt).eq.0)then
       vtac(imt)=0.d0
      else
       bbt=(xb(imt,1)-xxp)**2+(xb(imt,2)-yyp)**2
       vtac(imt)=qgfani(xpomr*scm,bbt,1.d0-exp(-sumut-sumup0)
     * ,0.d0,0.d0,iddt(imt),2,1)
       sumut=sumut+vtac(imt)
      endif
      
      if((abs(sumup-sumup0).gt..01d0.or.abs(sumut-sumut0).gt..01d0)
     *.and.nn.lt.100)then
       do i=1,ia(1)
        vpau0(i)=vpac(i)
       enddo
       do i=1,ia(2)
        vtau0(i)=vtac(i)
       enddo
       sumup0=sumup
       sumut0=sumut
       goto 1
      endif     
      
      sumcp=0.d0
      if(ip.gt.1)then
       do ipp=1,ip-1
	sumcp=sumcp+vpac(ipp)
       enddo
      endif            	
      sumct=0.d0
      if(it.gt.1)then
       do itt=1,it-1
        sumct=sumct+vtac(itt)
       enddo
      endif
      if(debug.ge.4)write (moniou,202)
      
201   format(2x,'qgfdf - configuration of fan contributions:'
     */2x,'xxp=',e10.3,2x,'yyp=',e10.3,2x,'xpomr=',e10.3
     *,2x,'ip=',i3,2x,'it=',i3,2x,'imp=',i3,2x,'imt=',i3)
202   format(2x,'qgfdf - end')
      return 
      end

c------------------------------------------------------------------------
      double precision function qgpomi(sy,bb,vvxp,vvxt,vvxpa,vvxta
     *,icdp,icdt,icz,iqq)
c-----------------------------------------------------------------------
c qgpomi - integrated  eikonal contributions
c sy          - Pomeron mass squared,
c bb          - squared impact parameter,
c vvxp, vvxpa - proj. nuclear screening (cut and uncut legs to other nucleons),
c vvxt, vvxta - targ. nuclear screening,
c icdp        - diffractive state for the projectile,
c icdt        - diffractive state for the target,
c icz         - projectile class
c iqq=1 - total uncut,
c iqq=2 - 1-cut,
c iqq=3 - multi-cut
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wz(3),wi(3),wj(3),wm(3),wn(3)
      common /qgarr10/ am(7),ammu
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr24/ qpomu(11,11,2,2,3),qpom1(11,11,6,4,3)
     *,qpomr(11,11,6,216,12)
      common /qgarr43/ moniou
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)sy,bb,vvxp,vvxt,vvxpa,vvxta
     *,icdp,icdt,icz,iqq

      qgpomi=0.d0
      if(sy.le.1.d0.and.iqq.eq.3)return                     !below threshold
      rp=(rq(icz)+rq(2)+alfp*log(max(1.d0,sy)))*4.d0*.0389d0
      z=exp(-bb/rp)
      if(z.lt..2d0*exp(-4.d0).and.iqq.le.2)then
       izmax=2
       jz=1
       wz(2)=5.d0*z*exp(4.d0)
       wz(1)=1.d0-wz(2)
      else
       izmax=3
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=max(2.d0,(-bb/rp-log(0.2d0))/.8d0+7.d0)
       endif
       jz=min(9,int(zz))
       jz=max(2,jz)
       if(jz.eq.6)jz=5
       wz(2)=zz-jz
       wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
       wz(1)=1.d0-wz(2)+wz(3)
       wz(2)=wz(2)-2.d0*wz(3)
      endif

      yl=dlog10((sy-am(2)**2-am(icz)**2)/2.d0/am(2))
      k=max(1,int(yl))
      k=min(k,9)     
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)

      if(iqq.eq.1)then                                !uncut eikonal
       do l1=1,izmax
        l2=jz+l1-1
       do k1=1,3
        k2=k+k1-1
        qgpomi=qgpomi+qpomu(k2,l2,icdp,icdt,icz)*wk(k1)*wz(l1)
       enddo
       enddo
      else                                            !cut eikonals

       if(iqq.eq.2)then
        vl=max(1.d0,vvxp*5.d0+1.d0)
        j=min(int(vl),4)    
        wj(2)=vl-j
        wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
        wj(1)=1.d0-wj(2)+wj(3)
        wj(2)=wj(2)-2.d0*wj(3)
        do j1=1,3
         j2=j+j1-1
        do l1=1,izmax
         l2=jz+l1-1
        do k1=1,3
         k2=k+k1-1
         qgpomi=qgpomi+qpom1(k2,l2,j2,icdp+2*(icdt-1),icz)
     *   *wk(k1)*wz(l1)*wj(j1)
        enddo
        enddo
        enddo       
       else
        ml=icdp+2*(icdt-1)+4*(icz-1)
        vl=max(1.d0,vvxp*5.d0+1.d0)
        j=min(int(vl),3)    
        wj(2)=vl-j
        wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
        wj(1)=1.d0-wj(2)+wj(3)
        wj(2)=wj(2)-2.d0*wj(3)

        vl1=max(1.d0,vvxt*5.d0+1.d0)
        i=min(int(vl1),3)    
        wi(2)=vl1-i
        wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
        wi(1)=1.d0-wi(2)+wi(3)
        wi(2)=wi(2)-2.d0*wi(3)

        vl2=max(1.d0,vvxpa*5.d0+1.d0)
        m=min(int(vl2),5)    
        wm(2)=vl2-m
        wm(1)=1.d0-wm(2)

        vl3=max(1.d0,vvxta*5.d0+1.d0)
        n=min(int(vl3),5)    
        wn(2)=vl3-n
        wn(1)=1.d0-wn(2)

        do n1=1,2
         n2=n+n1-2
        do m1=1,2
         m2=m+m1-2
        do i1=1,3
         i2=i+i1-1
        do j1=1,3
         j2=j+j1-1
        do l1=1,izmax
         l2=jz+l1-1
        do k1=1,3
         k2=k+k1-1
         qgpomi=qgpomi+qpomr(k2,l2,j2,i2+6*m2+36*n2,ml)
     *   *wk(k1)*wz(l1)*wj(j1)*wi(i1)*wm(m1)*wn(n1)
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
       endif
      endif
      qgpomi=exp(qgpomi)*z
      if(iqq.eq.3)qgpomi=qgpomi*(1.d0-vvxp)*(1.d0-vvxt)      
      if(debug.ge.4)write (moniou,202)qgpomi
      
201   format(2x,'qgpomi - integrated  eikonal contributions:'
     */2x,'sy=',e10.3,2x,'bb=',e10.3,2x,'vvxp=',e10.3,2x,'vvxt=',e10.3
     */2x,'vvxpa=',e10.3,2x,'vvxta=',e10.3,2x,'icdp=',i1,2x,'icdt=',i1
     *,2x,'icz=',i1,2x,'iqq=',i1)
202   format(2x,'qgpomi=',e10.3)
      return 
      end  
                                                     
c=============================================================================
      double precision function qgpint(sy,bb,vvx)
c-----------------------------------------------------------------------------
c qgpint - intermediate gg-Pomeron eikonal
c sy  - Pomeron mass squared,
c bb - squared impact parameter,
c vvx - relative strenth of screening corrections (0<vvx<1)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)sy,bb,vvx

      qgpint=0.d0
      s2min=4.d0*fqscal*qt0
      if(sy.lt.1.001d0*s2min)return

      alf=4.d0*pi*r3p/g3p*sigs*vvx
      xmin=(s2min/sy)**(delh+alf-dels)

      do i1=1,7
      do m1=1,2
       zh=(.5d0*(1.d0+xmin-(2*m1-3)*x1(i1)*(1.d0-xmin)))
     * **(1.d0/(delh+alf-dels))
       ww=zh*sy                        !c.m. energy squared for hard interaction
       sjqq=qgjit(qt0,qt0,ww,2,2)      !qq-ladder cross section
       sjqg=qgjit(qt0,qt0,ww,1,2)      !qg-ladder cross section
       sjgg=qgjit(qt0,qt0,ww,1,1)      !gg-ladder cross section
        
       stg=0.d0
       do i2=1,7
       do m2=1,2
        xx=.5d0*(1.d0+x1(i2)*(2*m2-3))
        xph=zh**xx
        xmh=zh/xph
        glu1=qgppdi(xph,vvx,0)  !gluon distribution (screened) for upper Pomeron
	sea1=qgppdi(xph,vvx,1)  !quark distribution (screened) for upper Pomeron
        glu2=qgppdi(xmh,vvx,0)  !gluon distribution (screened) for lower Pomeron
	sea2=qgppdi(xmh,vvx,1)  !quark distribution (screened) for lower Pomeron

        stg=stg+a1(i2)*(glu1*glu2*sjgg+(glu1*sea2+glu2*sea1)*sjqg
     *  +sea1*sea2*sjqq)
       enddo
       enddo
       rh=-alfp*dlog(zh)
       qgpint=qgpint-a1(i1)*dlog(zh)/zh**(delh+alf)*stg
     * *exp(-bb/(4.d0*.0389d0*rh))/rh
      enddo
      enddo
      qgpint=qgpint*rr**2*(1.d0-xmin)/(delh+alf-dels)*factk
     */2.d0*pi     
      if(debug.ge.4)write (moniou,202)qgpint
      
201   format(2x,'qgpint - intermediate gg-Pomeron eikonal:'
     */4x,'sy=',e10.3,2x,'bb=',e10.3,2x,'vvx=',e10.3)
202   format(2x,'qgpint=',e10.3)
      return
      end

c=============================================================================
      double precision function qgpint0(sy,vvx)
c-----------------------------------------------------------------------------
c qgpint0 - integrated (over b) gg-Pomeron eikonal
c sy  - Pomeron mass squared,
c vvx - relative strenth of screening corrections (0<vvx<1)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)sy,vvx

      qgpint0=0.d0
      s2min=4.d0*fqscal*qt0
      if(sy.lt.1.001d0*s2min)return

      alf=4.d0*pi*r3p/g3p*sigs*vvx
      xmin=(s2min/sy)**(delh+alf-dels)
c numerical integration over zh
      do i1=1,7
      do m1=1,2
       zh=(.5d0*(1.d0+xmin-(2*m1-3)*x1(i1)*(1.d0-xmin)))
     * **(1.d0/(delh+alf-dels))
       ww=zh*sy                        !c.m. energy squared for hard interaction
       sjqq=qgjit(qt0,qt0,ww,2,2)
       sjqg=qgjit(qt0,qt0,ww,1,2)
       sjgg=qgjit(qt0,qt0,ww,1,1)
        
       stg=0.d0
       do i2=1,7
       do m2=1,2
        xx=.5d0*(1.d0+x1(i2)*(2*m2-3))
        xph=zh**xx
        xmh=zh/xph
        glu1=qgppdi(xph,vvx,0)
	sea1=qgppdi(xph,vvx,1)
        glu2=qgppdi(xmh,vvx,0)
	sea2=qgppdi(xmh,vvx,1)

        stg=stg+a1(i2)*(glu1*glu2*sjgg+(glu1*sea2+glu2*sea1)*sjqg
     *  +sea1*sea2*sjqq)
       enddo
       enddo
       qgpint0=qgpint0-a1(i1)*dlog(zh)/zh**(delh+alf)*stg
      enddo
      enddo
      qgpint0=qgpint0*rr**2*(1.d0-xmin)/(delh+alf-dels)*factk
     **2.d0*.0389d0*pi     
      if(debug.ge.4)write (moniou,202)qgpint0
      
201   format(2x,'qgpint0 - integrated gg-Pomeron eikonal:'
     */4x,'sy=',e10.3,2x,'vvx=',e10.3)
202   format(2x,'qgpint0=',e10.3)
      return
      end

c------------------------------------------------------------------------
      double precision function qgpin0(sy,vvx)
c-----------------------------------------------------------------------
c qgpin0 - intermediate gg-Pomeron eikonal
c sy   - Pomeron mass squared,
c vvx - relative strenth of screening corrections (0<vvx<1),
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wj(3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr20/ spmax
      common /qgarr26/ factk,fqscal
      common /qgarr39/ qpomi(51,11,11),qpomi0(11,11)
      common /qgarr43/ moniou
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)sy,vvx

      qgpin0=0.d0
      s2min=4.d0*fqscal*qt0
      if(sy.lt.1.001d0*s2min)then
       if(debug.ge.4)write (moniou,202)qgpin0
       return
      endif
      
      yl=log(sy/s2min/2.d0)/log(spmax/s2min/2.d0)*10.d0+1.d0
      k=max(1,int(yl))
      k=min(k,9)     
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)

      vl=max(1.d0,vvx*10.d0+1.d0)
      j=min(int(vl),9)    
      wj(2)=vl-j
      wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
      wj(1)=1.d0-wj(2)+wj(3)
      wj(2)=wj(2)-2.d0*wj(3)

      do j1=1,3
       j2=j+j1-1
      do k1=1,3
       k2=k+k1-1
       qgpin0=qgpin0+qpomi0(k2,j2)*wk(k1)*wj(j1)
      enddo
      enddo
      qgpin0=exp(qgpin0)*sy**delh      
      if(debug.ge.4)write (moniou,202)qgpin0
      
201   format(2x,'qgpin0 - interpolation of gg-Pomeron eikonal:'
     */4x,'sy=',e10.3,2x,'vvx=',e10.3)
202   format(2x,'qgpin0=',e10.3)
      return 
      end

c------------------------------------------------------------------------
      double precision function qgpini(sy,bb,vvx)
c-----------------------------------------------------------------------
c qgpini - intermediate gg-Pomeron eikonal
c sy   - Pomeron mass squared,
c bb - squared impact parameter,
c vvx - relative strenth of screening corrections (0<vvx<1)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wj(3),wz(3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr20/ spmax
      common /qgarr39/ qpomi(51,11,11),qpomi0(11,11)
      common /qgarr43/ moniou
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)sy,bb,vvx

      qgpini=0.d0
      yl=log(sy)/log(spmax)*50.d0+1.d0
      k=max(1,int(yl))
      k=min(k,49)     
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)

      rp=alfp*dlog(max(1.d0,sy))*4.d0*.0389d0
      if(rp.le.1.d-10)then
       z=1.d0
      else
       z=exp(-bb/rp)
      endif
      if(z.lt..2d0*exp(-4.d0))then
       izmax=2
       jz=1
       wz(2)=5.d0*z*exp(4.d0)
       wz(1)=1.d0-wz(2)
      else
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=(-bb/rp-log(0.2d0))/.8d0+7.d0
       endif
       jz=min(9,int(zz))
       jz=max(2,jz)
       if(jz.eq.6)jz=5
       wz(2)=zz-jz
       wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
       wz(1)=1.d0-wz(2)+wz(3)
       wz(2)=wz(2)-2.d0*wz(3)
       izmax=3
      endif

      vl=max(1.d0,vvx*10.d0+1.d0)
      j=min(int(vl),9)    
      wj(2)=vl-j
      wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
      wj(1)=1.d0-wj(2)+wj(3)
      wj(2)=wj(2)-2.d0*wj(3)

      do j1=1,3
       j2=j+j1-1
      do l1=1,izmax
       l2=jz+l1-1
      do k1=1,3
       k2=k+k1-1
       qgpini=qgpini+qpomi(k2,l2,j2)*wk(k1)*wj(j1)*wz(l1)
      enddo
      enddo
      enddo
      qgpini=exp(qgpini)*sy**dels    
      if(debug.ge.4)write (moniou,202)qgpini
      
201   format(2x,'qgpini - interpolation of gg-Pomeron eikonal:'
     */4x,'sy=',e10.3,2x,'bb=',e10.3,2x,'vvx=',e10.3)
202   format(2x,'qgpini=',e10.3)
      return 
      end
   
c=============================================================================
      double precision function qglsh(sy,xp,bb,vvx,icdp,icz,iqq,jj)
c-----------------------------------------------------------------------------
c qglsh - unintegrated Pomeron leg eikonal (semihard)
c sy   - Pomeron mass squared,
c xp   - Pomeron LC momentum,
c bb   - squared impact parameter,
c vvx  - relative strenth of screening corrections (0<vvx<1),
c icdp -diffractive eigenstate for current nucleon (hadron),
c icz  - hadron class,
c iqq=0 - gluon contribution,
c iqq=1 - sea quark contribution,
c jj=0 - screening corrections from current nucleon (hadron) included in vvx,
c jj=1 - screening corrections from current nucleon (hadron) treated explicitely
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr19/ ahl(3)
      common /qgarr25/ ahv(3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)sy,xp,bb,vvx,icdp,icz,iqq,jj
      qglsh=0.d0
      s2min=4.d0*fqscal*qt0
      if(sy.lt.1.001d0*s2min)then
       if(debug.ge.4)write (moniou,202)qglsh
       return
      endif
       
      if(jj.eq.0)then
       vvx0=vvx
      else
       vvx0=vvx*(2.d0-vvx)
      endif
      alf=4.d0*pi*r3p/g3p*sigs*vvx0
      xmin=(s2min/sy)**(delh+alf-dels)

      do i1=1,7
      do m1=1,2
       zh=(.5d0*(1.d0+xmin-(2*m1-3)*x1(i1)*(1.d0-xmin)))
     * **(1.d0/(delh+alf-dels))
       ww=zh*sy                        !c.m. energy squared for hard interaction
       sjqq=qgjit(qt0,qt0,ww,2,2)      !qq-ladder cross section
       sjqg=qgjit(qt0,qt0,ww,1,2)      !qg-ladder cross section
       sjgg=qgjit(qt0,qt0,ww,1,1)      !gg-ladder cross section
        
       if(iqq.eq.0)then                !gluon contribution
        stg=0.d0
        do i2=1,7
        do m2=1,2
         xx=.5d0*(1.d0+x1(i2)*(2*m2-3))
         xph=zh**xx
         xmh=zh/xph
	 if(jj.eq.0)then
          glu1=qgppdi(xph,vvx,0)!gluon distribution (screened) for upper Pomeron
	  sea1=qgppdi(xph,vvx,1)!quark distribution (screened) for upper Pomeron
          glu2=qgppdi(xmh,vvx,0)!gluon distribution (screened) for lower Pomeron
	  sea2=qgppdi(xmh,vvx,1)!quark distribution (screened) for lower Pomeron
	 else
	  v1p1=qgfani(1.d0/xp/dsqrt(xph),bb,vvx,0.d0,0.d0,icdp,icz,1)
	  v1p2=qgfani(sy/xp*dsqrt(xmh),bb,vvx,0.d0,0.d0,icdp,icz,1)
	  vvx1=1.d0-(1.d0-vvx)**2*exp(-2.d0*v1p1) !new screening factors include
	  vvx2=1.d0-(1.d0-vvx)**2*exp(-2.d0*v1p2) ! fans from current hadron
          glu1=qgppdi(xph,vvx1,0)
	  sea1=qgppdi(xph,vvx1,1)
          glu2=qgppdi(xmh,vvx2,0)
	  sea2=qgppdi(xmh,vvx2,1)
	 endif
         rh=rq(icz)-alfp*dlog(zh*xp)
	 
         stsum=(glu1*glu2*sjgg+(glu1*sea2+glu2*sea1)*sjqg
     *   +sea1*sea2*sjqq)*exp(-bb/(4.d0*.0389d0*rh))/rh
         if(stsum.lt.0.d0)then
	  stsum=0.d0
	 endif
     
         stg=stg+a1(i2)*stsum
        enddo
        enddo
        qglsh=qglsh-a1(i1)*dlog(zh)/zh**(delh+alf)*stg

       elseif(iqq.eq.1)then                !quark contribution
        xmh=zh
	if(jj.eq.0)then
         glu2=qgppdi(xmh,vvx,0)
	 sea2=qgppdi(xmh,vvx,1)
	else
	 v1p2=qgfani(sy/xp*dsqrt(xmh),bb,vvx,0.d0,0.d0,icdp,icz,1)
	 vvx2=1.d0-(1.d0-vvx)**2*exp(-2.d0*v1p2)
         glu2=qgppdi(xmh,vvx2,0)
	 sea2=qgppdi(xmh,vvx2,1)
	endif
        rh=rq(icz)-alfp*dlog(zh)
	 
	stq=(glu2*sjqg+sea2*sjqq)*exp(-bb/(4.d0*.0389d0*rh))/rh
        if(stq.lt.0.d0)then
	 stq=0.d0
	endif
        qglsh=qglsh+a1(i1)/zh**(delh+alf)*stq
     *  *(qggrv(xp,qt0,icz,1)+qggrv(xp,qt0,icz,2))/dsqrt(xp)
       endif
      enddo
      enddo
      if(iqq.eq.0)then
       qglsh=qglsh*rr**2*(1.d0-xmin)/(delh+alf-dels)*fp(icz)*factk*g3p
     * /2.d0*pi
      elseif(iqq.eq.1)then
       qglsh=qglsh*rr*(1.d0-xmin)/(delh+alf-dels)*factk*g3p/4.d0
      endif
      if(jj.eq.1)qglsh=qglsh*cd(icdp,icz)     
      if(debug.ge.4)write (moniou,202)qglsh
      
201   format(2x,'qglsh - unintegrated Pomeron leg eikonal:'
     */4x,'sy=',e10.3,2x,'xp=',e10.3,2x,'b^2=',e10.3,2x,'vvx=',e10.3
     *,2x,'icdp=',i1,2x,'icz=',i1,2x,'iqq=',i1,2x,'jj=',i1)
202   format(2x,'qglsh=',e10.3)
      return
      end
      
c------------------------------------------------------------------------
      double precision function qglegc(sy,xp,bb,vvx,icdp,icz,iqq)
c-----------------------------------------------------------------------
c qglegc - interpolation of Pomeron leg eikonal
c sy   - Pomeron mass squared,
c xp   - Pomeron LC momentum,
c bb   - squared impact parameter,
c vvx - relative strenth of screening corrections (0<vvx<1),
c icdp - diffractive eigenstate for the hadron,
c icz  - hadron class
c iqq=-1 - total contribution,
c iqq=0  - soft contribution,
c iqq=1  - gluon contribution,
c iqq=2  - sea quark contribution
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wj(3),wi(3),wz(3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr19/ ahl(3)
      common /qgarr20/ spmax
      common /qgarr25/ ahv(3)
      common /qgarr26/ factk,fqscal
      common /qgarr35/ qlegc(11,10,11,11,24)
      common /qgarr43/ moniou
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)sy,xp,bb,vvx,icdp,icz,iqq

      qglegc=0.d0
      s2min=4.d0*fqscal*qt0
      if(iqq.gt.0.and.sy.lt.1.001d0*s2min.or.iqq.eq.2.and.xp.gt..99d0)
     *then
       if(debug.ge.4)write (moniou,202)qglegc
       return
      endif
      
      rp=(rq(icz)+alfp*log(max(1.d0,sy/xp)))*4.d0*.0389d0
      z=exp(-bb/rp)
      if(sy.le.1.d0)then
       qglegc=sy**dels*fp(icz)*sigs*g3p/rp*4.d0*.0389d0
     * /xp**alpd*z*cd(icdp,icz)       
       if(debug.ge.4)write (moniou,202)qglegc
       return
      endif

      if(iqq.le.0.and.z.lt..2d0*exp(-4.d0))then
       izmax=2
       jz=1
       wz(2)=5.d0*z*exp(4.d0)
       wz(1)=1.d0-wz(2)
      else
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=(-bb/rp-log(0.2d0))/.8d0+7.d0
       endif
       if(zz.lt.2.d0)then
        jz=2
        wz(2)=zz-jz
        wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
        wz(1)=1.d0-wz(2)+wz(3)
        wz(2)=wz(2)-2.d0*wz(3)
        izmax=3
       else
        jz=min(10,int(zz))
        jz=max(2,jz)
        wz(2)=zz-jz
        wz(1)=1.d0-wz(2)
        izmax=2
       endif
      endif

      if(iqq.le.0)then
       yl=max(0.d0,dlog(sy)/dlog(spmax))*10.d0+1.d0
      else
       yl=dlog(sy/s2min/2.d0)/dlog(spmax/s2min/2.d0)*10.d0+1.d0
      endif
      k=max(1,int(yl))
      k=min(k,10)     
      wk(2)=yl-k
      wk(1)=1.d0-wk(2)
      iymax=2

      vl=max(1.d0,vvx*10.d0+1.d0)
      if(vl.lt.2.d0)then
       j=1    
       wj(2)=vl-j
       wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
       wj(1)=1.d0-wj(2)+wj(3)
       wj(2)=wj(2)-2.d0*wj(3)
       ivmax=3
      else
       j=min(int(vl),10)    
       wj(2)=vl-j
       wj(1)=1.d0-wj(2)
       ivmax=2
      endif

      if(xp.lt..2d0.and.sy.lt.1.5d0*spmax)then
       xl=6.d0-5.d0*log(5.d0*xp)/log(sy/spmax/2.d0)
      else
       xl=5.d0*xp+5.d0
      endif
      if(xl.lt.1.d0.or.sy.ge.1.5d0*spmax)then
       i=min(8,int(xl))
       i=max(1,i)
       if(i.eq.5)i=4
       if(sy.ge.1.5d0*spmax)i=max(i,6)
       wi(2)=xl-i
       wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
       wi(1)=1.d0-wi(2)+wi(3)
       wi(2)=wi(2)-2.d0*wi(3)
       ixmax=3
      else
       i=min(9,int(xl))
       i=max(1,i)
       wi(2)=xl-i
       wi(1)=1.d0-wi(2)
       ixmax=2
      endif
      
      do j1=1,ivmax
       j2=j+j1-1
      do l1=1,izmax
       l2=jz+l1-1
      do i1=1,ixmax
       i2=i+i1-1
      do k1=1,iymax
       k2=k+k1-1
       qglegc=qglegc+qlegc(k2,i2,l2,j2,icdp+2*(icz-1)+6*(iqq+1))
     * *wk(k1)*wi(i1)*wz(l1)*wj(j1)
      enddo
      enddo
      enddo
      enddo
      if(iqq.le.0)then
       qglegc=exp(qglegc)*sy**dels*z/xp**alpd
      else
       qglegc=exp(qglegc)*sy**delh*z
       if(iqq.eq.1)then
        qglegc=qglegc/xp**alpd
       else
        qglegc=qglegc/dsqrt(xp)*(1.d0-xp)**(ahv(icz)-ahl(icz))
       endif
      endif
      if(debug.ge.4)write (moniou,202)qglegc
      
201   format(2x,'qglegc - interpolation of Pomeron leg eikonal:'
     */4x,'sy=',e10.3,2x,'xp=',e10.3,2x,'b^2=',e10.3,2x,'vvx=',e10.3
     *,2x,'icdp=',i1,2x,'icz=',i1,2x,'iqq=',i1)
202   format(2x,'qglegc=',e10.3)
      return 
      end
      
c=============================================================================
      double precision function qgpsh(sy,xpp,xpm,b,vvx0
     *,icdp,icdt,icz,iqq)
c-----------------------------------------------------------------------------
c qgpsh - unintegrated semihard Pomeron eikonal
c sy         - Pomeron mass squared,
c xpp, xpm   - Pomeron LC momenta,
c b          - impact parameter,
c vvx0       - relative strenth of nuclear screening corrections (0<vvx0<1),
c icdp, icdt - proj. and targ. diffractive eigenstates,
c icz        - hadron class,
c iqq        - type of the hard interaction (0-gg, 1-q_vg, 2-gq_v)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr25/ ahv(3)
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /qgarr51/ x4(2),a4(2)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)sy,xpp,xpm,b,vvx0,icdp,icdt
     *,icz,iqq
      qgpsh=0.d0      
      s2min=4.d0*fqscal*qt0               !energy threshold for hard interaction
      if(s2min/sy.ge.1.d0)then
       if(debug.ge.4)write (moniou,202)qgpsh
       return
      endif
      
      if(iqq.ne.2)then
       icv=icz
       icq=2
       xp=xpp
       xm=xpm
       icdv=icdp
       icdq=icdt
      else
       icv=2
       icq=icz
       xp=xpm
       xm=xpp
       icdq=icdp
       icdv=icdt
      endif
      if(iqq.eq.0)then
       rp=(rq(icz)+rq(2)+alfp*log(max(1.d0,sy/xp/xm/s2min)))
     * *4.d0*.0389d0
      else
       rp=(rq(icz)+rq(2)+alfp*log(max(1.d0,sy/xm/s2min)))*4.d0*.0389d0
      endif
      
      do ib1=1,2
      do mb1=1,2
       z=.5+x4(ib1)*(mb1-1.5)
       b1=sqrt(-rp/4.*log(z))
       dbb=b1**2+.25*b**2
      do ib2=1,2
      do mb2=1,2
       phi=pi*(.5+x4(ib2)*(mb2-1.5))
       bb1=dbb+b*b1*cos(phi)
       bb2=2.*dbb-bb1  

       v1pnu0=qgfani(1.d0/xpp*dsqrt(sy),bb1,vvx0,0.d0,0.d0,icdp,icz,1)
       v1tnu0=qgfani(1.d0/xpm*dsqrt(sy),bb2,vvx0,0.d0,0.d0,icdt,2,1)
       nn=0
1      nn=nn+1
       vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx0)
       vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx0)
       v1tnu=qgfani(1.d0/xpm*dsqrt(sy),bb2,vvxt,0.d0,0.d0,icdt,2,1)
       v1pnu=qgfani(1.d0/xpp*dsqrt(sy),bb1,vvxp,0.d0,0.d0,icdp,icz,1)
       if((abs(v1pnu0-v1pnu).gt.1.d-2.or.abs(v1tnu0-v1tnu).gt.1.d-2)
     * .and.nn.lt.100)then
        v1pnu0=v1pnu
        v1tnu0=v1tnu
	goto 1
       endif
       vvx=1.d0-exp(-2.d0*(v1pnu+v1tnu))*(1.d0-vvx0)**2
       alf=4.d0*pi*r3p*sigs/g3p*vvx
       
       xmin=(s2min/sy)**(delh+alf-dels)
       dqg=0.d0

       do i=1,7 
       do m=1,2
        z1=(.5d0*(1.d0+xmin-(2*m-3)*x1(i)*(1.d0-xmin)))
     *  **(1.d0/(delh+alf-dels))
        ww=z1*sy
        sjqq=qgjit(qt0,qt0,ww,2,2)
        sjqg=qgjit(qt0,qt0,ww,1,2)
        sjgg=qgjit(qt0,qt0,ww,1,1)
        
        if(iqq.eq.0)then                                !gg-Pomeron
         st2=0.d0
         do j=1,2
         do k=1,2
          xx=.5d0*(1.d0+x4(j)*(2*k-3))
          xph=z1**xx
          xmh=z1/xph
	  
          v1pnu0=qgfani(1.d0/xpp/dsqrt(xph),bb1,vvx0,0.d0,0.d0
     *    ,icdp,icz,1)
          v1tnu0=qgfani(sy/xpm*dsqrt(xph),bb2,vvx0,0.d0,0.d0,icdt,2,1)
          nn=0
2         nn=nn+1
          vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx0)
          vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx0)
          v1tnu=qgfani(sy/xpm*dsqrt(xph),bb2,vvxt,0.d0,0.d0,icdt,2,1)
          v1pnu=qgfani(1.d0/xpp/dsqrt(xph),bb1,vvxp,0.d0,0.d0
     *    ,icdp,icz,1)
          if((abs(v1pnu0-v1pnu).gt.1.d-2.or.abs(v1tnu0-v1tnu).gt.1.d-2)
     *    .and.nn.lt.100)then
           v1pnu0=v1pnu
           v1tnu0=v1tnu
	   goto 2
          endif
          vvx1=1.d0-exp(-2.d0*(v1pnu+v1tnu))*(1.d0-vvx0)**2
	  
          v1pnu0=qgfani(sy/xpp*dsqrt(xmh),bb1,vvx0,0.d0,0.d0,icdp,icz,1)
          v1tnu0=qgfani(1.d0/xpm/dsqrt(xmh),bb2,vvx0,0.d0,0.d0,icdt,2,1)
          nn=0
3         nn=nn+1
          vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx0)
          vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx0)
          v1tnu=qgfani(1.d0/xpm/dsqrt(xmh),bb2,vvxt,0.d0,0.d0,icdt,2,1)
          v1pnu=qgfani(sy/xpp*dsqrt(xmh),bb1,vvxp,0.d0,0.d0,icdp,icz,1)
          if((abs(v1pnu0-v1pnu).gt.1.d-2.or.abs(v1tnu0-v1tnu).gt.1.d-2)
     *    .and.nn.lt.100)then
           v1pnu0=v1pnu
           v1tnu0=v1tnu
	   goto 3
          endif
          vvx2=1.d0-exp(-2.d0*(v1pnu+v1tnu))*(1.d0-vvx0)**2
	  
          glu1=qgppdi(xph,vvx1,0)
          sea1=qgppdi(xph,vvx1,1)
          glu2=qgppdi(xmh,vvx2,0)
          sea2=qgppdi(xmh,vvx2,1)
          st2=st2+a4(j)*(glu1*glu2*sjgg+(glu1*sea2+glu2*sea1)*sjqg
     *    +sea1*sea2*sjqq)
         enddo
         enddo
         rh1=rq(icz)-alfp*dlog(xpp*xph)
         rh2=rq(2)-alfp*dlog(xpm*xmh)
         dqg=dqg-a1(i)*dlog(z1)/z1**(delh+alf)*st2
     *   *exp(-(bb1/rh1+bb2/rh2)/4.d0/.0389d0)/rh1/rh2
      
        else                                !qg-Pomeron
         xmh=z1
	 if(iqq.eq.1)then
	  bbv=bb1
	  bbq=bb2
	 else
	  bbv=bb2
	  bbq=bb1
	 endif
         v1pnu0=qgfani(sy/xp*dsqrt(xmh),bbv,vvx0,0.d0,0.d0,icdv,icv,1)
         v1tnu0=qgfani(1.d0/xm/dsqrt(xmh),bbq,vvx0,0.d0,0.d0,icdq,icq,1)
         nn=0
4        nn=nn+1
         vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx0)
         vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx0)
         v1tnu=qgfani(1.d0/xm/dsqrt(xmh),bbq,vvxt,0.d0,0.d0,icdq,icq,1)
         v1pnu=qgfani(sy/xp*dsqrt(xmh),bbv,vvxp,0.d0,0.d0,icdv,icv,1)
         if((abs(v1pnu0-v1pnu).gt.1.d-2.or.abs(v1tnu0-v1tnu).gt.1.d-2)
     *   .and.nn.lt.100)then
          v1pnu0=v1pnu
          v1tnu0=v1tnu
	  goto 4
         endif
         vvx2=1.d0-exp(-2.d0*(v1pnu+v1tnu))*(1.d0-vvx0)**2
	  
         glu=qgppdi(xmh,vvx2,0)
         sea=qgppdi(xmh,vvx2,1)
         rh1=rq(icv)
         rh2=rq(icq)-alfp*dlog(xm*xmh)
         fst=(glu*sjqg+sea*sjqq)
     *   *(qggrv(xp,qt0,icv,1)+qggrv(xp,qt0,icv,2))/dsqrt(xp)
     *   *exp(-(bbv/rh1+bbq/rh2)/4.d0/.0389d0)/rh1/rh2
         dqg=dqg+a1(i)/z1**(delh+alf)*fst
        endif
       enddo
       enddo
       dqg=dqg*(1.d0-xmin)/(delh+alf-dels)
       qgpsh=qgpsh+a4(ib1)*a4(ib2)*dqg/z
      enddo
      enddo
      enddo
      enddo

      if(iqq.eq.0)then
       qgpsh=qgpsh*rr**2*fp(icz)*fp(2)*factk/128.d0*pi*rp/.0389d0
     * *cd(icdp,icz)*cd(icdt,2)
      else
       qgpsh=qgpsh*rr*fp(icq)*factk/256.d0*rp/.0389d0
     * *cd(icdp,icz)*cd(icdt,2)
      endif     
      if(debug.ge.4)write (moniou,202)qgpsh
      
201   format(2x,'qgpsh - unintegrated semihard Pomeron eikonal:'
     */4x,'sy=',e10.3,2x,'xpp=',e10.3,2x,'xpm=',e10.3,2x,'b=',e10.3
     */4x,'vvx0=',e10.3,2x,'icdp=',i1,2x,'icdt=',i1,2x,'icz=',i1
     *,2x,'iqq=',i1)
202   format(2x,'qgpsh=',e10.3)
      return
      end	
      
c=============================================================================
      double precision function qgpsoft(sy,xpp,xpm,b,vvx0,icdp,icdt,icz)
c-----------------------------------------------------------------------------
c qgpsoft - soft Pomeron eikonal
c sy         - Pomeron mass squared,
c xpp, xpm   - Pomeron LC momenta,
c b          - impact parameter,
c vvx0       - relative strenth of nuclear screening corrections (0<vvx0<1),
c icdp, icdt - proj. and targ. diffractive eigenstates,
c icz        - hadron class
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /qgarr51/ x4(2),a4(2)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)sy,xpp,xpm,b,vvx0,icdp,icdt,icz
      qgpsoft=0.d0      
      rp=(rq(icz)+rq(2)+alfp*log(max(1.d0,sy/xpp/xpm)))*4.d0*.0389d0
      rh1=rq(icz)+alfp*log(max(1.d0,dsqrt(sy)/xpp))
      rh2=rq(2)+alfp*log(max(1.d0,dsqrt(sy)/xpm))
      
      do ib1=1,7 
      do mb1=1,2
       z=.5+x1(ib1)*(mb1-1.5)
       b1=sqrt(-rp/4.*log(z))
       dbb=b1**2+.25*b**2
      do ib2=1,7
      do mb2=1,2
       phi=pi*(.5+x1(ib2)*(mb2-1.5))
       bb1=dbb+b*b1*cos(phi)
       bb2=2.*dbb-bb1  

       v1pnu0=qgfani(1.d0/xpp*dsqrt(sy),bb1,vvx0,0.d0,0.d0,icdp,icz,1)
       v1tnu0=qgfani(1.d0/xpm*dsqrt(sy),bb2,vvx0,0.d0,0.d0,icdt,2,1)
       nn=0
1      nn=nn+1
       vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx0)
       vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx0)
       v1tnu=qgfani(1.d0/xpm*dsqrt(sy),bb2,vvxt,0.d0,0.d0,icdt,2,1)
       v1pnu=qgfani(1.d0/xpp*dsqrt(sy),bb1,vvxp,0.d0,0.d0,icdp,icz,1)
       if((abs(v1pnu0-v1pnu).gt.1.d-2.or.abs(v1tnu0-v1tnu).gt.1.d-2)
     * .and.nn.lt.100)then
        v1pnu0=v1pnu
        v1tnu0=v1tnu
	goto 1
       endif
       vvx=1.d0-exp(-2.d0*(v1pnu+v1tnu))*(1.d0-vvx0)**2
       alf=4.d0*pi*r3p*sigs/g3p*vvx

       qgpsoft=qgpsoft+a1(ib1)*a1(ib2)*sy**dels/max(sy,1.d0)**alf
     * /z*exp(-(bb1/rh1+bb2/rh2)/4.d0/.0389d0)/rh1/rh2 
      enddo
      enddo   
      enddo
      enddo   
      qgpsoft=qgpsoft*fp(icz)*fp(2)*sigs*rp/64.d0/.0389d0
     **cd(icdp,icz)*cd(icdt,2)/(xpp*xpm)**alpd
      if(debug.ge.4)write (moniou,202)qgpsoft
      
201   format(2x,'qgpsoft - soft Pomeron eikonal:'
     */4x,'sy=',e10.3,2x,'xpp=',e10.3,2x,'xpm=',e10.3,2x,'b=',e10.3
     */4x,'vvx0=',e10.3,2x,'icdp=',i1,2x,'icdt=',i1,2x,'icz=',i1)
202   format(2x,'qgpsoft=',e10.3)
      return
      end	
     
c=============================================================================
      double precision function qgls(sy,xp,bb,vvx,icdp,icz,jj)
c-----------------------------------------------------------------------------
c qgls - soft Pomeron leg eikonal
c sy   - Pomeron mass squared,
c xp   - Pomeron light cone momentum,
c bb   - squared impact parameter to the 3p-vertex,
c vvx  - relative strenth of screening corrections,
c icdp - diffractive eigenstate for the connected hadron,
c icz  - hadron class,
c jj=0 - screening corrections from current nucleon (hadron) included in vvx,
c jj=1 - screening corrections from current nucleon (hadron) treated explicitely
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /qgarr51/ x4(2),a4(2)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)sy,xp,bb,vvx,icdp,icz,jj
      rp=rq(icz)+alfp*log(max(1.d0,sy/xp))
      if(jj.eq.0.or.bb.gt.1.d9)then                !analytic result
       alf=4.d0*pi*r3p/g3p*sigs*vvx
       qgls=sy**dels/max(1.d0,sy)**alf*fp(icz)*sigs*g3p/rp
       if(bb.lt.1.d9)qgls=qgls*exp(-bb/(4.d0*.0389d0*rp))
       if(jj.eq.1)qgls=qgls*cd(icdp,icz)
       if(debug.ge.4)write (moniou,202)qgls
       return       
      endif
      
      alf=4.d0*pi*r3p/g3p*sigs*vvx*(2.d0-vvx)
      qgls=0.d0      
      do ix1=1,7
      do mx1=1,2
       xpomr1=xp/sy**(.5d0+x1(ix1)*(mx1-1.5d0))
       rp1=alfp*log(xpomr1/xp*sy)*4.d0*.0389d0	
       do ix2=1,2
       do mx2=1,2
        bb1=-rp1*log(.5d0+x4(ix2)*(mx2-1.5d0))
       do ix3=1,2
       do mx3=1,2
        phi=2.d0*pi*(.5d0+x4(ix3)*(mx3-1.5d0))
	bb2=(dsqrt(bb)-dsqrt(bb1)*cos(phi))**2+bb1*sin(phi)**2
       
        rp2=rq(icz)-alfp*log(xpomr1)
        v1s=max(1.d0,xp/xpomr1)**(-alf)
     *  /rp2*exp(-bb2/(4.d0*.0389d0*rp2))
        v1p=qgfani(1.d0/xpomr1,bb2,vvx,0.d0,0.d0,icdp,icz,1)
        v1pi=qgfani(dsqrt(sy/xp/xpomr1),dsqrt(bb*bb2),vvx,0.d0,0.d0
     *  ,icdp,icz,1)
        alf1=4.d0*pi*r3p/g3p*sigs*(1.d0-(1.d0-vvx)**2*exp(-2.d0*v1pi))
        v1i=max(1.d0,xpomr1/xp*sy)**(-alf1)

        qgls=qgls+a1(ix1)*a4(ix2)*a4(ix3)*v1i*v1s*(exp(-2.d0*v1p)-1.d0)
       enddo
       enddo
       enddo
       enddo
      enddo
      enddo
      qgls=(qgls/2.d0*pi*r3p/g3p*sigs*(1.d0-vvx)**2*dlog(sy)
     *+max(1.d0,sy)**(-alf)/rp*exp(-bb/(4.d0*.0389d0*rp)))
     **sy**dels*fp(icz)*sigs*g3p*cd(icdp,icz)
      if(debug.ge.4)write (moniou,202)qgls
      
201   format(2x,'qgls - soft Pomeron leg eikonal:'
     */4x,'sy=',e10.3,2x,'xp=',e10.3,2x,'b^2=',e10.3,2x,'vvx=',e10.3
     *,2x,'icdp=',i1,2x,'icz=',i1,2x,'jj=',i1)
202   format(2x,'qgls=',e10.3)
      return
      end	
      
c=============================================================================
      double precision function qgpomc(sy,xp,xm,bb,vvx
     *,icdp,icdt,icz,iqq)
c-----------------------------------------------------------------------
c qgpomc - unintegrated cut Pomeron eikonal
c sy         - Pomeron mass squared,
c xp,xm      - Pomeron light cone momenta,
c bb         - squared impact parameter,
c vvx        - relative strenth of nuclear screening corrections,
c icdp, icdt - proj. and targ. diffractive eigenstates,
c icz        - hadron class
c iqq=-1 - total,
c iqq=0  - soft contribution,
c iqq=1  - gg contribution,
c iqq=2  - qg contribution
c iqq=3  - gq contribution
c iqq=4  - qq contribution
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3),wi(3),wj(3),wz(3),wm(3)
      common /qgarr6/  pi,bm,amws
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr19/ ahl(3)
      common /qgarr20/ spmax
      common /qgarr25/ ahv(3)
      common /qgarr26/ factk,fqscal
      common /qgarr38/ qpomc(11,100,11,11,60)
      common /qgarr43/ moniou
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)sy,xp,xm,bb,vvx
     *,icdp,icdt,icz,iqq

      qgpomc=0.d0
      s2min=4.d0*fqscal*qt0
      if(iqq.gt.0.and.sy.lt.1.001d0*s2min.or.iqq.eq.4
     *.and.(xp.gt..99d0.or.xm.gt..99d0).or.iqq.eq.2.and.xp.gt..99d0
     *.or.iqq.eq.3.and.xm.gt..99d0)then
       if(debug.ge.4)write (moniou,202)qgpomc
       return
      endif
      
      if(iqq.eq.4)then                          !qq contribution
       sj=qgjit(qt0,qt0,sy,2,2)
       qgpomc=sj*factk*(qggrv(xp,qt0,icz,1)+qggrv(xp,qt0,icz,2))
     * *(qggrv(xm,qt0,2,1)+qggrv(xm,qt0,2,2))/xp/xm
     * *(1.d0-xp)**(ahv(icz)-ahl(icz))*(1.d0-xm)**(ahv(2)-ahl(2))
     * *exp(-bb/(4.d0*.0389d0*(rq(icz)+rq(2))))
     * /(8.d0*pi*(rq(icz)+rq(2)))*cd(icdp,icz)*cd(icdt,2)
       return
      endif
      
      rp=(rq(icz)+rq(2)+alfp*log(sy/xp/xm))*4.d0*.0389d0
      z=exp(-bb/rp)
      if(sy.le.1.d0)then
       qgpomc=sy**dels*fp(icz)*fp(2)*sigs/(xp*xm)**alpd*z/rp
     * *4.d0*.0389d0*cd(icdp,icz)*cd(icdt,2)
       return
      endif

      if(iqq.le.0.and.z.lt..2d0*exp(-4.d0))then
       izmax=2
       jz=1
       wz(2)=5.d0*z*exp(4.d0)
       wz(1)=1.d0-wz(2)
      else
       if(z.gt..2d0)then
        zz=5.d0*z+6.d0
       else
        zz=(-bb/rp-log(0.2d0))/.8d0+7.d0
       endif
       if(zz.lt.2.d0)then
        jz=2
        wz(2)=zz-jz
        wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
        wz(1)=1.d0-wz(2)+wz(3)
        wz(2)=wz(2)-2.d0*wz(3)
        izmax=3
       else
        jz=min(10,int(zz))
        jz=max(2,jz)
        wz(2)=zz-jz
        wz(1)=1.d0-wz(2)
        izmax=2
       endif
      endif

      if(iqq.le.0)then                      !total or soft
       yl=max(0.d0,dlog(sy)/dlog(spmax))*10.d0+1.d0
      else                                  !semihard
       yl=dlog(sy/s2min/2.d0)/dlog(spmax/s2min/2.d0)*10.d0+1.d0
      endif
      k=max(1,int(yl))
      k=min(k,10)     
      wk(2)=yl-k
      wk(1)=1.d0-wk(2)
      iymax=2
      
      if(xp.lt..2d0.and.sy.lt.1.5d0*spmax)then
       xl1=6.d0-5.d0*log(5.d0*xp)/log(sy/spmax/2.d0)
      else
       xl1=5.d0*xp+5.d0
      endif
      if(xl1.lt.1.d0.or.sy.ge.1.5d0*spmax)then
       i=min(8,int(xl1))
       i=max(1,i)
       if(i.eq.5)i=4
       if(sy.ge.1.5d0*spmax)i=max(i,6)
       wi(2)=xl1-i
       wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
       wi(1)=1.d0-wi(2)+wi(3)
       wi(2)=wi(2)-2.d0*wi(3)
       ix1max=3
      else
       i=min(9,int(xl1))
       i=max(1,i)
       wi(2)=xl1-i
       wi(1)=1.d0-wi(2)
       ix1max=2
      endif
      
      if(xm.lt..2d0.and.sy.lt.1.5d0*spmax)then
       xl2=6.d0-5.d0*log(5.d0*xm)/log(sy/spmax/2.d0)
      else
       xl2=5.d0*xm+5.d0
      endif
      if(xl2.lt.1.d0.or.sy.ge.1.5d0*spmax)then
       j=min(8,int(xl2))
       j=max(1,j)
       if(j.eq.5)j=4
       if(sy.ge.1.5d0*spmax)j=max(j,6)
       wj(2)=xl2-j
       wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
       wj(1)=1.d0-wj(2)+wj(3)
       wj(2)=wj(2)-2.d0*wj(3)
       ix2max=3
      else
       j=min(9,int(xl2))
       j=max(1,j)
       wj(2)=xl2-j
       wj(1)=1.d0-wj(2)
       ix2max=2
      endif

      ml=icdp+2*(icdt-1)+4*(icz-1)+12*(iqq+1)
      if(vvx.eq.0.d0)then                     !hadron-proton collision
       do l1=1,izmax
        l2=jz+l1-1
       do j1=1,ix2max
        j2=j+j1-1
       do i1=1,ix1max
        i2=i+i1-1
       do k1=1,iymax
        k2=k+k1-1
        qgpomc=qgpomc+qpomc(k2,i2+10*(j2-1),l2,1,ml)
     *  *wk(k1)*wi(i1)*wj(j1)*wz(l1)
       enddo
       enddo
       enddo
       enddo
      else                                    !hA (AA) collision
       vl=max(1.d0,vvx*10.d0+1.d0)
       if(vl.lt.2.d0)then
        m=1    
        wm(2)=vl-m
        wm(3)=wm(2)*(wm(2)-1.d0)*.5d0
        wm(1)=1.d0-wm(2)+wm(3)
        wm(2)=wm(2)-2.d0*wm(3)
        ivmax=3
       else
        m=min(int(vl),10)    
        wm(2)=vl-m
        wm(1)=1.d0-wm(2)
        ivmax=2
       endif

       do m1=1,ivmax
        m2=m+m1-1
       do l1=1,izmax
        l2=jz+l1-1
       do j1=1,ix2max
        j2=j+j1-1
       do i1=1,ix1max
        i2=i+i1-1
       do k1=1,iymax
        k2=k+k1-1
        qgpomc=qgpomc+qpomc(k2,i2+10*(j2-1),l2,m2,ml)
     *  *wk(k1)*wi(i1)*wj(j1)*wz(l1)*wm(m1)
       enddo
       enddo
       enddo
       enddo
       enddo
      endif
      qgpomc=exp(qgpomc)*z
      if(iqq.le.0)then
       qgpomc=qgpomc/(xp*xm)**alpd*sy**dels
      elseif(iqq.eq.1)then
       qgpomc=qgpomc/(xp*xm)**alpd*sy**delh
      elseif(iqq.eq.2)then
       qgpomc=qgpomc/xm**alpd/dsqrt(xp)*(1.d0-xp)**(ahv(icz)-ahl(icz))
     * *sy**delh
      elseif(iqq.eq.3)then
       qgpomc=qgpomc/xp**alpd/dsqrt(xm)*(1.d0-xm)**(ahv(2)-ahl(2))
     * *sy**delh
      endif     
      if(debug.ge.4)write (moniou,202)qgpomc
      
201   format(2x,'qgpomc - unintegrated cut Pomeron eikonal:'
     */4x,'sy=',e10.3,2x,'xp=',e10.3,2x,'xm=',e10.3,2x,'b^2=',e10.3
     */4x,'vvx=',e10.3,2x,'icdp=',i1,2x,'icdt=',i1,2x,'icz=',i1
     *,2x,'iqq=',i1)
202   format(2x,'qgpomc=',e10.3)
      return 
      end

c=============================================================================
      subroutine qgv(x,y,xb,vin,vdd,vabs)
c qgv - eikonal factors for nucleus-nucleus interaction
c (used for cross-section calculation)
c x,y - projectile nucleon coordinates
c----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207)
cdh   dimension xb(iapmax,3),fhard(3),vabs(2)
      dimension xb(iapmax,3),         vabs(2)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)x,y

      vin=0.d0
      vdd=0.d0
      do iddp1=1,2
       dv=0.d0
       do m=1,ia(2)
        bb=(x-xb(m,1))**2+(y-xb(m,2))**2
        dv=dv+qgpomi(scm,bb,0.d0,0.d0,0.d0,0.d0,iddp1,iddt(m),icz,1)
       enddo
       dv=exp(-dv)
       vabs(iddp1)=1.d0-dv**2       !1-exp(-2 * chi_i)
       vdd=vdd+cc(iddp1,icz)*dv**2  !sum_i cc(i) exp(-2 * chi_i)
       vin=vin+cc(iddp1,icz)*dv     !sum_i cc(i) exp(-chi_i)
      enddo
      vin=1.d0-vin**2               !1-sum_ij cc(i) cc(j) exp(-chi_i-chi_j)
      vdd=vdd+vin-1.d0 
          !sum_i cc(i) exp(-2*chi_i) - sum_ij cc(i) cc(j) exp(-chi_i-chi_j)

      if(debug.ge.4)write (moniou,202)vin,vdd,vabs
      if(debug.ge.4)write (moniou,203)
      
201   format(2x,'qgv - eikonal factor: nucleon coordinates x=',
     *e10.3,2x,'y=',e10.3)
202   format(2x,'vin=',e10.3,2x,'vdd=',e10.3,2x,'vabs=',2e10.3)
203   format(2x,'qgv - end')
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
c (0 - intact, 1 - inel., 2 - diffr., -1 - recoiled from diffraction),
c ncola(1-iap), ncolb(1-iat) - index for inel.-wounded proj. and targ. nucleons,
c nbpom  - total number of Pomeron blocks,
c ias(k) (ibs(k)) - index of the proj. (targ.) nucleon for k-th Pomeron block,
c bbpom(k) - squared impact parameter (between proj. and targ.) for k-th block,
c vvxpom(k) - relative strenth of A-screening corrections for k-th block,
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
c vvxpr(n,ip) - relative strenth of A-screening corrections for n-th leg 
c Pomeron connected to ip-th proj. nucleon,
c vvxtg(n,it) - relative strenth of A-screening corrections for n-th leg 
c Pomeron connected to it-th targ. nucleon,
c xpompr(n,ip) - LC momentum of the 3P-vertex for n-th leg Pomeron connected 
c to ip-th proj. nucleon,
c xpomtg(n,it) - LC momentum of the 3P-vertex for n-th leg Pomeron connected 
c to it-th targ. nucleon
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207,npbmax=1000,npnmax=1000,npmax=5000
     *,legmax=900)
      dimension xas(iapmax,3),vabs(2),vabsi(2,iapmax),wdifi(iapmax)
     *,vpac(iapmax),vtac(iapmax),xpomip(npmax),xpomim(npmax)
     *,xpompi(legmax),xpomti(legmax),vvxpi(legmax),vvxti(legmax)
     *,bpompi(legmax),bpomti(legmax),ipompi(legmax),ipomti(legmax)
     *,vvcont(iapmax),ncola(iapmax),ncolb(iapmax)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),b
      common /qgarr9/  iwp(iapmax),iwt(iapmax),lqa(iapmax),lqb(iapmax)
     *,iprcn(iapmax),itgcn(iapmax),ias(npbmax),ibs(npbmax),nqs(npbmax)
     *,npompr(npbmax),npomtg(npbmax),npomin(npbmax),nnpr(npmax,npbmax)
     *,nntg(npmax,npbmax),ilpr(legmax,npbmax),iltg(legmax,npbmax)
     *,lnpr(legmax,npbmax),lntg(legmax,npbmax)
     *,nbpi(npnmax,iapmax),nbti(npnmax,iapmax),idnpi(npnmax,iapmax)
     *,idnti(npnmax,iapmax),nppi(npnmax,iapmax),npti(npnmax,iapmax)
     *,nlpi(npnmax,iapmax),nlti(npnmax,iapmax)
      common /qgarr11/ b10
c nsp - number of secondary particles
      common /qgarr12/ nsp
      common /qgarr13/ nsf,iaf(iapmax)   
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr22/ wppr(iapmax),wmtg(iapmax)
      common /qgarr23/ bbpom(npbmax),vvxpom(npbmax)
     *,bpompr(npnmax,iapmax),bpomtg(npnmax,iapmax),vvxpr(npnmax,iapmax)
     *,vvxtg(npnmax,iapmax),xpompr(npnmax,iapmax),xpomtg(npnmax,iapmax)
     *,xpopin(npmax,npbmax),xpomin(npmax,npbmax)
      common /qgarr43/ moniou
      common /qgarr46/ iconab(iapmax,iapmax),icona(iapmax)
     *,iconb(iapmax)
      common /qgarr55/ nwt,nwp  !number of wounded targ.(proj.) nucleons
      common /debug/   debug
      external psran

      if(debug.ge.1)write (moniou,201)     
      nsp=0
      nsf=0
      nsp0=nsp  
      
c initialization
1     continue
      do i=1,ia(1)
       iddp(i)=1+int(psran(b10)+cc(2,icz)) !diffractive eigenstates for proj.
      enddo
      do i=1,ia(2)
       iddt(i)=1+int(psran(b10)+cc(2,2))   !diffractive eigenstates for targ.
      enddo

c-------------------------------------------------
c squared impact parameter is sampled uniformly (b**2<bm**2)
      b=bm*dsqrt(psran(b10))
      if(debug.ge.1)write (moniou,202)b

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
       vvcont(it)=0.d0
      enddo
      
      do ip=1,ia(1)
       icdp=iddp(ip)             
       vvconp=0.d0
       do it=1,ia(2)
        icdt=iddt(it)                        
	bbp=(xa(ip,1)+b-xb(it,1))**2+(xa(ip,2)-xb(it,2))**2   
        vv1p=qgpomi(scm,bbp,0.d0,0.d0,0.d0,0.d0,icdp,icdt,icz,1)
	if(vv1p.gt.1.d-2)then
         if(debug.ge.2)write (moniou,205)ip,it      
	 iconab(ip,it)=1
	else
	 iconab(ip,it)=0
	endif
        vvcont(it)=vvcont(it)+vv1p
        vvconp=vvconp+vv1p
       enddo
       if(vvconp.gt.1.d-2)then
        if(debug.ge.2)write (moniou,206)ip      
	icona(ip)=1
       else
	icona(ip)=0
       endif
      enddo
      
      do it=1,ia(2)
       if(vvcont(it).gt.1.d-2)then
        if(debug.ge.2)write (moniou,207)it      
	iconb(it)=1
       else
	iconb(it)=0
       endif
      enddo
      nrej=0

2     nrej=nrej+1
      if(debug.ge.2)write (moniou,208)nrej     
cdh 2006.02.21
cdh   if(nrej.gt.100)then
      if(nrej.gt.5)then
       if(debug.ge.1)write (moniou,209)
       goto 1
      endif
      nsp=nsp0
      nbpom=0
      nwp=0
      nwt=0
      do i=1,ia(1)
       lqa(i)=0
       iwp(i)=0
       ncola(i)=0
       wppr(i)=wp0                         !LC+ for projectile nucleons
      enddo
      do i=1,ia(2)
       lqb(i)=0
       iwt(i)=0
       ncolb(i)=0
       wmtg(i)=wm0                         !LC- for target nucleons
      enddo

c-------------------------------------------------
c Pomeron configuration
      if(debug.ge.1)write (moniou,210)
      do 5 ip=1,ia(1)             !loop over all projectile nucleons
       if(debug.ge.2)write (moniou,211)ip
       if(icona(ip).eq.0)goto 5
       x=xa(ip,1)+b               !proj. x is shifted by the impact parameter b
       y=xa(ip,2)
       icdp=iddp(ip)              !diffr. eigenstate for ip

       do 3 it=1,ia(2)            !loop over all target nucleons
        if(debug.ge.2)write (moniou,212)it
        if(iconab(ip,it).eq.0)goto 3
        icdt=iddt(it)                         !diffr. eigenstate for it
	bbp=(x-xb(it,1))**2+(y-xb(it,2))**2   !distance squared between ip, it

c calculate nuclear screening factors for "middle point" -> eikonals
        xpomr=1.d0/dsqrt(scm)
	xxp=.5d0*(x+xb(it,1))
	yyp=.5d0*(y+xb(it,2))
        call qgfdf(xxp,yyp,xpomr,vpac,vtac,sumup,sumut,sumcp,sumct
     *  ,ip,it,ip,it)
        vvx0=1.d0-exp(min(0.d0,vpac(ip)+vtac(it)-sumup-sumut))
        vv1p=qgpomi(scm,bbp,vvx0,0.d0,0.d0,0.d0,icdp,icdt,icz,2) !1-Pomeron eikonal
        call qgfdf(xxp,yyp,xpomr,vpac,vtac,sumup,sumut,sumcp,sumct
     *  ,ip,it,ia(1),ia(2))
        vvxp=1.d0-exp(-sumcp)
        vvxt=1.d0-exp(-sumct)
        vvxpa=1.d0-exp(min(0.d0,vpac(ip)+sumcp-sumup))
        vvxta=1.d0-exp(min(0.d0,vtac(it)+sumct-sumut))
	vv3p=qgpomi(scm,bbp,vvxp,vvxt,vvxpa,vvxta,icdp,icdt,icz,3) !multi-Pomeron
        vv=vv1p+vv3p                                        !total eikonal
        if(debug.ge.2)write (moniou,213)vv,vv1p,vv3p
	
        if(psran(b10).gt.1.d0-exp(-2.d0*vv))goto 3  !1.-exp(-2*vv) - probability
	                                            !for inelastic interaction
        iwt(it)=1
	iwp(ip)=1
	ncola(ip)=ncola(ip)+1
	ncolb(it)=ncolb(it)+1
	
        n=npgen(2.d0*vv,1,50) !number of cut Pomerons for (ip-it) interaction	
        nbpom=nbpom+1         !new Pomeron block
        if(nbpom.gt.npbmax)then
	 goto 2
	endif
        ias(nbpom)=ip         !proj. index for current elementary interaction
        ibs(nbpom)=it         !targ. index for current elementary interaction
        bbpom(nbpom)=bbp      !distance squared between ip, it
        vvxpom(nbpom)=vvx0    !relative strenth of A-screening corrections
        if(debug.ge.2)write (moniou,214)nbpom,ip,it,n
	
        nqs(nbpom)=0
	npomin(nbpom)=0
	npompr(nbpom)=0
	npomtg(nbpom)=0
	do i=1,n
	 aks=psran(b10)*vv
	 if(aks.gt.vv3p)then               !single Pomeron
          if(debug.ge.2)write (moniou,215)i
          np=nqs(nbpom)+1
          if(np.gt.legmax)then
	   goto 2
	  endif
          nqs(nbpom)=np                    !update Pomeron number in the block
	  l0=lqa(ip)+1 
          if(l0.gt.npnmax)then
	   goto 2
	  endif
          lqa(ip)=l0                       !update number of connections for proj.
	  nnpr(np,nbpom)=l0                !index for connected proj. participant
	  nbpi(l0,ip)=nbpom
	  idnpi(l0,ip)=0
	  nppi(l0,ip)=np
	  l0=lqb(it)+1 
          if(l0.gt.npnmax)then
	   goto 2
	  endif
          lqb(it)=l0
	  nntg(np,nbpom)=l0                !index for connected targ. participant
	  nbti(l0,it)=nbpom
	  idnti(l0,it)=0
	  npti(l0,it)=np
	  
	 else                              !multi-Pomeron vertex
c configuration of multi-Pomeron vertexes	  
          if(debug.ge.2)write (moniou,219)
          call qg3pdf(vvxpi,vvxti,xpompi,xpomti,bpompi,bpomti,xpomip
     *    ,xpomim,npompi,npomti,npin,ipompi,ipomti,ip,it,iret)
          if(iret.ne.0)goto 2

	  if(npin.ne.0)then
           if(debug.ge.2)write (moniou,220)i,npin
	   npomin(nbpom)=npomin(nbpom)+npin
	   if(npomin(nbpom).gt.npmax)then
	    goto 2
	   endif
	   do l=1,npin
	    l1=npomin(nbpom)+l-npin
	    xpopin(l1,nbpom)=xpomip(l)
	    xpomin(l1,nbpom)=xpomim(l)
	   enddo
	  endif
	  if(npompi.ne.0)then
           if(debug.ge.2)write (moniou,221)i,npompi
	   do m=1,npompi
	    np=npompr(nbpom)+1
            if(np.gt.legmax)then
	     goto 2
	    endif
	    npompr(nbpom)=np
	    ipp=ipompi(m)
            iwp(ipp)=1
	    ilpr(np,nbpom)=ipp
	    l0=lqa(ipp)+1 
	    if(l0.gt.npnmax)then
	     goto 2
            endif
            lqa(ipp)=l0
	    lnpr(np,nbpom)=l0
	    nbpi(l0,ipp)=nbpom
	    idnpi(l0,ipp)=1
	    nlpi(l0,ipp)=np
	    vvxpr(l0,ipp)=vvxpi(m)
	    xpompr(l0,ipp)=1.d0/xpompi(m)/scm
	    bpompr(l0,ipp)=bpompi(m)
	   enddo
	  endif
	  if(npomti.ne.0)then  
           if(debug.ge.2)write (moniou,222)i,npomti
	   do m=1,npomti
	    np=npomtg(nbpom)+1
            if(np.gt.legmax)then
	     goto 2
	    endif
	    npomtg(nbpom)=np
	    itt=ipomti(m)
            iwt(itt)=1
	    iltg(np,nbpom)=itt
	    l0=lqb(itt)+1 
	    if(l0.gt.npnmax)then
	     goto 2
            endif
            lqb(itt)=l0
	    lntg(np,nbpom)=l0
	    nbti(l0,itt)=nbpom
	    idnti(l0,itt)=1
	    nlti(l0,itt)=np
	    vvxtg(l0,itt)=vvxti(m)
	    xpomtg(l0,itt)=xpomti(m)
	    bpomtg(l0,itt)=bpomti(m)
	   enddo
	  endif
	 endif	
	enddo             !end of Pomeron loop
3      continue           !end of it-loop
      
       if(iwp(ip).ne.0)then
        nwp=nwp+1                         !one more wounded proj. nucleon
       else                               !check projectile diffraction
        if(debug.ge.2)write (moniou,223)ip
        vabs(1)=0.d0
        vabs(2)=0.d0
	icdps=iddp(ip)
        do it=1,ia(2)
 	 bbp=(x-xb(it,1))**2+(y-xb(it,2))**2
         icdt=iddt(it)
	 do icdp=1,2
          if(iconab(ip,it).eq.0)then
	   vabsi(icdp,it)=0.d0
	  else
	   iddp(ip)=icdp	 
           xpomr=1.d0/dsqrt(scm)
	   xxp=.5d0*(x+xb(it,1))
	   yyp=.5d0*(y+xb(it,2))
           call qgfdf(xxp,yyp,xpomr,vpac,vtac,sumup,sumut,sumcp,sumct
     *     ,ip,it,ip,it)
           vvx=1.d0-exp(min(0.d0,vpac(ip)+vtac(it)-sumup-sumut))
           vv1p=qgpomi(scm,bbp,vvx,0.d0,0.d0,0.d0,icdp,icdt,icz,2) !1-Pomeron eikonal
           call qgfdf(xxp,yyp,xpomr,vpac,vtac,sumup,sumut,sumcp,sumct
     *     ,ip,it,ia(1),ia(2))
           vvxp=1.d0-exp(-sumcp)
           vvxt=1.d0-exp(-sumct)
           vvxpa=1.d0-exp(min(0.d0,vpac(ip)+sumcp-sumup))
           vvxta=1.d0-exp(min(0.d0,vtac(it)+sumct-sumut))
	   vv3p=qgpomi(scm,bbp,vvxp,vvxt,vvxpa,vvxta,icdp,icdt,icz,3) !multi-Pomeron
           vv=vv1p+vv3p                                     !total eikonal

	   vabsi(icdp,it)=vv
           vabs(icdp)=vabs(icdp)+vv
	  endif
	 enddo
        enddo
	iddp(ip)=icdps	
        wdifr=cc(1,icz)*cc(2,icz)*(exp(-vabs(1))-exp(-vabs(2)))**2
     *  /(cc(1,icz)*exp(-2.d0*vabs(1))+cc(2,icz)*exp(-2.d0*vabs(2)))
     
        if(psran(b10).lt.wdifr)then       !projectile diffraction
	 wdift=0.d0
         do it=1,ia(2)
	  if(iwt(it).ne.-1)then
	   wdifi(it)=cc(1,icz)*cc(2,icz)*(exp(-vabsi(1,it))
     *     -exp(-vabsi(2,it)))**2/(cc(1,icz)*exp(-2.d0*vabsi(1,it))
     *     +cc(2,icz)*exp(-2.d0*vabsi(2,it)))
	   wdift=wdift+wdifi(it)
          else
	   wdifi(it)=0.d0
	  endif
	 enddo
	 if(wdift.ne.0.d0)then
          nwp=nwp+1
	  iwp(ip)=2
	  aks=psran(b10)*wdift
          do it=1,ia(2)
	   aks=aks-wdifi(it)
	   if(aks.lt.0.d0)goto 4
	  enddo
4	  continue
          iprcn(ip)=it
	  if(iwt(it).eq.0)iwt(it)=-1
          if(debug.ge.2)write (moniou,224)ip,it
	 endif
        endif
       endif
5     continue                            !end of ip-loop

      do 8 it=1,ia(2)                     !check target diffraction
       if(iwt(it).gt.0)then 
        nwt=nwt+1                         !one more wounded targ. nucleon
       elseif(iconb(it).ne.0)then
        if(debug.ge.2)write (moniou,225)it
        vabs(1)=0.d0
	vabs(2)=0.d0
	icdts=iddt(it)
        do ip=1,ia(1)
	 bbp=(xa(ip,1)+b-xb(it,1))**2+(xa(ip,2)-xb(it,2))**2
         icdp=iddp(ip)
	 do icdt=1,2
          if(iconab(ip,it).eq.0)then
	   vabsi(icdt,ip)=0.d0
	  else
	   iddt(it)=icdt	 
           xpomr=1.d0/dsqrt(scm)
	   xxp=.5d0*(xa(ip,1)+b+xb(it,1))
	   yyp=.5d0*(xa(ip,2)+xb(it,2))
           call qgfdf(xxp,yyp,xpomr,vpac,vtac,sumup,sumut,sumcp,sumct
     *     ,ip,it,ip,it)
           vvx=1.d0-exp(min(0.d0,vpac(ip)+vtac(it)-sumup-sumut))
           vv1p=qgpomi(scm,bbp,vvx,0.d0,0.d0,0.d0,icdp,icdt,icz,2) !1-Pomeron eikonal
           call qgfdf(xxp,yyp,xpomr,vpac,vtac,sumup,sumut,sumcp,sumct
     *     ,ip,it,ia(1),ia(2))
           vvxp=1.d0-exp(-sumcp)
           vvxt=1.d0-exp(-sumct)
           vvxpa=1.d0-exp(min(0.d0,vpac(ip)+sumcp-sumup))
           vvxta=1.d0-exp(min(0.d0,vtac(it)+sumct-sumut))
	   vv3p=qgpomi(scm,bbp,vvxp,vvxt,vvxpa,vvxta,icdp,icdt,icz,3) !multi-Pomeron
           vv=vv1p+vv3p                                        !total eikonal

	   vabsi(icdt,ip)=vv
           vabs(icdt)=vabs(icdt)+vv
	  endif
	 enddo
	enddo
	iddt(it)=icdts	
	wdifr=cc(1,2)*cc(2,2)*(exp(-vabs(1))-exp(-vabs(2)))**2
     *  /(cc(1,2)*exp(-2.d0*vabs(1))+cc(2,2)*exp(-2.d0*vabs(2)))
     
        if(psran(b10).lt.wdifr)then       !target diffraction
	 wdift=0.d0
         do ip=1,ia(1)
	  if(iwp(ip).eq.-1)then
	   wdifi(ip)=0.d0
	  else
	   if(iwp(ip).eq.2)then
	    itt=iprcn(ip)
	    if(itt.eq.it)goto 7
	    if(iwt(itt).eq.2)then
	     wdifi(ip)=0.d0
	     goto 6
	    endif
	   endif	    
	   wdifi(ip)=cc(1,2)*cc(2,2)*(exp(-vabsi(1,ip))
     *     -exp(-vabsi(2,ip)))**2/(cc(1,2)*exp(-2.d0*vabsi(1,ip))
     *     +cc(2,2)*exp(-2.d0*vabsi(2,ip)))
          endif
6	  wdift=wdift+wdifi(ip)
	 enddo
	 if(wdift.eq.0.d0)goto 8
         nwt=nwt+1
	 iwt(it)=2
	 aks=psran(b10)*wdift
         do ip=1,ia(1)
	  aks=aks-wdifi(ip)
	  if(aks.lt.0.d0)goto 7
	 enddo
7	 continue
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
8     continue 

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

c inelastic interaction: energy sharing and particle production 
      if(nwp.ne.0.or.nwt.ne.0)then
       if(debug.ge.1)write (moniou,229)
       call qgsha(nbpom,ncola,ncolb)
       if(nsp.le.nsp0+2)then
        if(debug.ge.1)write (moniou,230)
        goto 1
       endif
      else                                 !no interaction
       if(debug.ge.1)write (moniou,231)
       goto 1
      endif
      if(debug.ge.1)write (moniou,232)nsp

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
213   format(2x,'qgconf: eikonals - total: ',e10.3,2x,'single: ',e10.3
     *,2x,'multi-P: ',e10.3)
214   format(2x,'qgconf: ',i4,'-th Pomeron block connected to ',i3
     *,'-th proj. nucleon and'/4x,i3,'-th targ. nucleon;'
     *,' number of element. processes in the block: ',i3)
215   format(2x,'qgconf: ',i3
     *,'-th process in the block is single cut Pomeron')
216   format(2x,'qgconf: ',i3,'-th process in the block contains'
     *,' multi-Pomeron vertex')
217   format(2x,'qgconf: ',i3,'-th process in the block has a'
     *,' rap. gap with the proj.')
218   format(2x,'qgconf: ',i3,'-th process in the block has a'
     *,' rap. gap with the targ.')
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
230   format(2x,'qgconf: no particle produced - rejection')
231   format(2x,'qgconf: no interaction - rejection')
232   format(2x,'qgconf: ',i5,' particles have been produced')
233   format(2x,'qgconf: fragmentation of the proj. spectator part')
234   format(2x,'qgconf: ',i3,' proj. fragments have been produced')
235   format(2x,'qgconf - end')
      return
      end

c=============================================================================
      subroutine qg3pdf(vvxpi,vvxti,xpompi,xpomti,bpompi,bpomti
     *,xpomip,xpomim,nppr,nptg,npin,ipompi,ipomti,ip,it,iret)
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
c vvxpi(i) - relative strenth of scr. corrections for i-th proj. leg Pomeron,
c vvxti(i) - relative strenth of scr. corrections for i-th targ. leg Pomeron,
c xpompi(i) - LC momentum of the 3P-vertex for i-th proj. leg Pomeron,
c xpomti(i) - LC momentum of the 3P-vertex for i-th targ. leg Pomeron
c iret=1 - reject configuration
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207,npmax=5000,levmax=20,legmax=900)
      dimension vpac(iapmax),vtac(iapmax)
     *,vpac0(iapmax),vtac0(iapmax),vpact0(iapmax),vtact0(iapmax)
     *,xpomip(npmax),xpomim(npmax),xpompi(legmax),xpomti(legmax)
     *,vvxpi(legmax),vvxti(legmax),bpompi(legmax),bpomti(legmax)
     *,ipompi(legmax),ipomti(legmax),ippr0(legmax),iptg0(legmax)
     *,nppm(levmax),ippm(npmax,levmax),ii(levmax),xpomm(levmax)
     *,itypr0(legmax),itytg0(legmax),itypm(legmax,levmax)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),b
      common /qgarr11/ b10
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr19/ ahl(3)
      common /qgarr43/ moniou
      common /qgarr46/ iconab(iapmax,iapmax),icona(iapmax)
     *,iconb(iapmax)
      common /debug/   debug
      external psran

      if(debug.ge.2)write (moniou,201)ip,it               !so161205
      
      iret=0
c normalization of rejection function      
      alf=4.d0*pi*r3p*sigs/g3p
      xpomr=1.d0/dsqrt(scm)	   
      rp=(rq(icz)+rq(2)+alfp*log(scm))*4.d0*.0389d0  
      xxp=.5d0*(xa(ip,1)+b+xb(it,1))    
      yyp=.5d0*(xa(ip,2)+xb(it,2))
      bpt=dsqrt((xa(ip,1)+b-xb(it,1))**2+(xa(ip,2)-xb(it,2))**2)
      call qgfdf(xxp,yyp,xpomr,vpac,vtac,sumup,sumut,sumcp,sumct
     *,ip,it,ia(1),ia(2))
      bbp=(xa(ip,1)+b-xxp)**2+(xa(ip,2)-yyp)**2    
      bbt=(xxp-xb(it,1))**2+(yyp-xb(it,2))**2    
      v1p=min(vpac(ip),qgfani(1.d0/xpomr,bbp,1.d0-exp(-sumut-sumcp)
     *,1.d0-exp(min(0.d0,sumcp+2.d0*vpac(ip)-2.d0*sumup))
     *,0.d0,iddp(ip),icz,2))
      v1t=min(vtac(it),qgfani(xpomr*scm,bbt,1.d0-exp(-sumup-sumct)
     *,1.d0-exp(min(0.d0,sumct+2.d0*vtac(it)-2.d0*sumut))
     *,0.d0,iddt(it),2,2))

      gb0
     *=max(0.d0,(1.d0-exp(-vpac(ip)))*exp(-sumcp)-vpac(ip)*exp(-sumup))
     **max(0.d0,(1.d0-exp(-vtac(it)))*exp(-sumct)-vtac(it)*exp(-sumut))
     *-.25d0*((1.d0-exp(-vpac(ip)))**2*exp(-2.d0*sumcp)
     *+2.d0*(1.d0-exp(-vpac(ip)))*(1.d0-exp(-sumcp))*exp(-sumcp))
     **((1.d0-exp(-vtac(it)))**2*exp(-2.d0*sumct)
     *+2.d0*(1.d0-exp(-vtac(it)))*(1.d0-exp(-sumct))*exp(-sumct))
     *+.5d0*(vtac(it)+v1t)*exp(-sumut)
     **max(0.d0,(1.d0-exp(-vpac(ip)))*exp(-sumcp)-vpac(ip)*exp(-sumup))
     *+.5d0*(vpac(ip)+v1p)*exp(-sumup)
     **max(0.d0,(1.d0-exp(-vtac(it)))*exp(-sumct)-vtac(it)*exp(-sumut))
     *+.5d0*exp(-sumup-sumut)
     **(v1p*(1.d0-exp(-sumup))*(vtac(it)+v1t*exp(-sumut))
     *+v1t*(1.d0-exp(-sumut))*(vpac(ip)+v1p*exp(-sumup)))
     
      gbcor10=(1.d0-(1.d0-(1.d0-xpomr)**(1.d0+ahl(icz)))
     ***(1.d0+dels-alpd))
      gbcor20=(1.d0-(1.d0-(1.d0-1.d0/xpomr/scm)
     ***(1.d0+ahl(2)))**(1.d0+dels-alpd))	  
      gb0=min(1.d0,gb0)*max(1.d0,dlog10(scm)*2.d0)
      if(gb0.le.0.d0.or.gbcor10.le.0.d0.or.gbcor20.le.0.d0)then	  
       if(debug.ge.3)write (moniou,202)
       iret=1
       goto 21
      endif
      gb0=gb0/gbcor10/gbcor20/(gbcor10+gbcor20)
      if(debug.ge.3)write (moniou,203)gb0
      
1     continue
      xpomr=scm**(-psran(b10))               !proposed LC momentum for 3P-vertex    
      gbcor1=(1.d0-(1.d0-(1.d0-xpomr)**(1.d0+ahl(icz)))
     ***(1.d0+dels-alpd))
      gbcor2=(1.d0-(1.d0-(1.d0-1.d0/xpomr/scm)
     ***(1.d0+ahl(2)))**(1.d0+dels-alpd))	  
      if(psran(b10).gt.gbcor1*gbcor2*(gbcor1+gbcor2)
     */gbcor10/gbcor20/(gbcor10+gbcor20))goto 1
      
      z=psran(b10)
      b1=sqrt(-rp/4.*log(z))
      dbb=b1**2+.25*bpt**2
      phi=pi*psran(b10)
      bbpr=dbb+bpt*b1*cos(phi)            !squared impact parameter to the proj.
      bbtg=2.*dbb-bbpr                    !squared impact parameter to the targ.
      call qgbdef(bbpr,bbtg,xa(ip,1)+b,xa(ip,2),xb(it,1),xb(it,2)
     *,xxp,yyp,int(1.5d0+psran(b10)))     !determine coordinates for the vertex      
      call qgfdf(xxp,yyp,xpomr,vpac,vtac,sumup,sumut,sumcp,sumct
     *,ip,it,ia(1),ia(2))
                          	          !calculate fan contributions	
      v1p=min(vpac(ip),qgfani(1.d0/xpomr,bbpr,1.d0-exp(-sumut-sumcp)
     *,1.d0-exp(min(0.d0,sumcp+2.d0*vpac(ip)-2.d0*sumup))
     *,0.d0,iddp(ip),icz,2))
      v1t=min(vtac(it),qgfani(xpomr*scm,bbtg,1.d0-exp(-sumup-sumct)
     *,1.d0-exp(min(0.d0,sumct+2.d0*vtac(it)-2.d0*sumut))
     *,0.d0,iddt(it),2,2))
     
      gb
     *=max(0.d0,(1.d0-exp(-vpac(ip)))*exp(-sumcp)-vpac(ip)*exp(-sumup))
     **max(0.d0,(1.d0-exp(-vtac(it)))*exp(-sumct)-vtac(it)*exp(-sumut))
     *-.25d0*((1.d0-exp(-vpac(ip)))**2*exp(-2.d0*sumcp)
     *+2.d0*(1.d0-exp(-vpac(ip)))*(1.d0-exp(-sumcp))*exp(-sumcp))
     **((1.d0-exp(-vtac(it)))**2*exp(-2.d0*sumct)
     *+2.d0*(1.d0-exp(-vtac(it)))*(1.d0-exp(-sumct))*exp(-sumct))
     *+.5d0*(vtac(it)+v1t)*exp(-sumut)
     **max(0.d0,(1.d0-exp(-vpac(ip)))*exp(-sumcp)-vpac(ip)*exp(-sumup))
     *+.5d0*(vpac(ip)+v1p)*exp(-sumup)
     **max(0.d0,(1.d0-exp(-vtac(it)))*exp(-sumct)-vtac(it)*exp(-sumut))
     *+.5d0*exp(-sumup-sumut)
     **(v1p*(1.d0-exp(-sumup))*(vtac(it)+v1t*exp(-sumut))
     *+v1t*(1.d0-exp(-sumut))*(vpac(ip)+v1p*exp(-sumup)))
      gbs=gb
      gb=min(gb,1.d0)/gbcor1/gbcor2/(gbcor1+gbcor2)/gb0/z      
      if(debug.ge.5)write (moniou,204)xpomr,bbpr,bbtg,gb
     
      if(psran(b10).gt.gb)goto 1
      if(debug.ge.3)write (moniou,205)xpomr,bbpr,bbtg,xxp,yyp

      sumcp0=0.d0
      sumcpt=0.d0
      sumcptf=0.d0
      do i=1,ia(1)
       ipp=ia(1)-i+1
       bbp=(xa(ipp,1)+b-xxp)**2+(xa(ipp,2)-yyp)**2
       if(ipp.ge.ip)vpac0(ipp)=min(vpac(ipp),qgfani(1.d0/xpomr,bbp
     * ,1.d0-exp(-sumut-sumcp),1.d0-exp(-sumcp0),1.d0-exp(-sumut)
     * ,iddp(ipp),icz,5))
       vpact0(ipp)=qgfani(1.d0/xpomr,bbp,1.d0-exp(-sumut-sumcp)
     * ,1.d0-exp(-sumcptf),1.d0-exp(-sumcp),iddp(ipp),icz,6)
       sumcptf=sumcptf+vpact0(ipp)
       if(ipp.gt.ip)then
        sumcp0=sumcp0+vpac0(ipp)
        sumcpt=sumcpt+vpact0(ipp)
       endif
      enddo
      sumct0=0.d0
      sumctt=0.d0
      sumcttf=0.d0
      do i=1,ia(2)
       itt=ia(2)-i+1
       bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
       if(itt.ge.it)vtac0(itt)=min(vtac(itt),qgfani(xpomr*scm,bbt
     * ,1.d0-exp(-sumup-sumct),1.d0-exp(-sumct0),1.d0-exp(-sumup)
     * ,iddt(itt),2,5))
       vtact0(itt)=qgfani(xpomr*scm,bbt,1.d0-exp(-sumup
     * -sumct),1.d0-exp(-sumcttf),1.d0-exp(-sumct),iddt(itt),2,6)
       sumcttf=sumcttf+vtact0(itt)
       if(itt.gt.it)then
        sumct0=sumct0+vtac0(itt)
        sumctt=sumctt+vtact0(itt)
       endif
      enddo
      v1p0=min(v1p,qgfani(1.d0/xpomr,bbpr,1.d0-exp(-sumut-sumcp)
     *,1.d0-exp(min(0.d0,sumcp+2.d0*vpac(ip)-2.d0*sumup))
     *,1.d0-exp(-sumut),iddp(ip),icz,4))
      v1t0=min(v1t,qgfani(xpomr*scm,bbtg,1.d0-exp(-sumup-sumct)
     *,1.d0-exp(min(0.d0,sumct+2.d0*vtac(it)-2.d0*sumut))
     *,1.d0-exp(-sumup),iddt(it),2,4))
           	       	
c weights for vertex contributions: 
c vv1: 1 proj. leg and >1 targ. legs
       vv1=max(0.d0,(1.d0-exp(-2.d0*vpac(ip)))*exp(-2.d0*sumcp)
     *-2.d0*vpac(ip)*exp(-2.d0*sumup))
     **max(0.d0,(1.d0-exp(-2.d0*vtac(it)))*exp(-2.d0*sumct)
     *-2.d0*vtac(it)*exp(-2.d0*sumut))
     *-2.d0*max(0.d0,(exp(-vpac0(ip))-exp(-vpac(ip)))
     **exp(-sumcp0-sumcp)-(vpac(ip)-vpac0(ip))*exp(-sumup))
     **max(0.d0,(1.d0-exp(-2.d0*vtac(it)))*exp(-2.d0*sumct)
     *-2.d0*vtac(it)*exp(-2.d0*sumut))
     *-2.d0*max(0.d0,(exp(-vtac0(it))-exp(-vtac(it)))
     **exp(-sumct0-sumct)-(vtac(it)-vtac0(it))*exp(-sumut))
     **max(0.d0,(1.d0-exp(-2.d0*vpac(ip)))*exp(-2.d0*sumcp)
     *-2.d0*vpac(ip)*exp(-2.d0*sumup))
c vv2: >1 proj. legs and 1 targ. leg
      vv2=max(0.d0,(1.d0-exp(-2.d0*vpac(ip)))*exp(-2.d0*sumcp)
     *-2.d0*vpac(ip)*exp(-2.d0*sumup))
     **((vtac0(it)+v1t0)*exp(-2.d0*sumut)
     *-(vtac(it)-vtac0(it)+v1t-v1t0)*exp(-sumut)*(1.d0-exp(-sumut)))
     *-2.d0*max(0.d0,(exp(-vpac0(ip))-exp(-vpac(ip)))
     **exp(-sumcp0-sumcp)-(vpac(ip)-vpac0(ip))*exp(-sumup))
     **(vtac(it)+v1t)*exp(-2.d0*sumut)     
c vv3: >1 proj. legs and >1 targ. legs
      vv3=((vpac0(ip)+v1p0)*exp(-2.d0*sumup)
     *-(vpac(ip)-vpac0(ip)+v1p-v1p0)*exp(-sumup)*(1.d0-exp(-sumup)))
     **max(0.d0,(1.d0-exp(-2.d0*vtac(it)))*exp(-2.d0*sumct)
     *-2.d0*vtac(it)*exp(-2.d0*sumut))
     *-2.d0*max(0.d0,(exp(-vtac0(it))-exp(-vtac(it)))
     **exp(-sumct0-sumct)-(vtac(it)-vtac0(it))*exp(-sumut))
     **(vpac(ip)+v1p)*exp(-2.d0*sumup)            
c vv4: 0 proj. legs and >1 targ. legs (targ. diffr.)
      vv4=max(0.d0,(1.d0-exp(-2.d0*vtac(it)))*exp(-2.d0*sumct)
     *-2.d0*vtac(it)*exp(-2.d0*sumut))
     **((1.d0-exp(-vpac(ip)))**2*exp(-2.d0*sumcp)
     *+2.d0*(1.d0-exp(-vpac(ip)))*(1.d0-exp(-sumcp))*exp(-sumcp))     
c vv5: 0 proj. legs and 1 targ. leg (targ. diffr.)
      vv5=((vtac0(it)+v1t0)*exp(-2.d0*sumut)
     *-(vtac(it)-vtac0(it)+v1t-v1t0)*exp(-sumut)*(1.d0-exp(-sumut)))
     **((1.d0-exp(-vpac(ip)))**2*exp(-2.d0*sumcp)
     *+2.d0*(1.d0-exp(-vpac(ip)))*(1.d0-exp(-sumcp))*exp(-sumcp))
c vv6: 0 proj. legs and >1 targ. legs with "holes" (targ. diffr.)
      vv6=2.d0*(1.d0-exp(-2.d0*vpac(ip)))*exp(-2.d0*sumcp)
     **max(0.d0,(exp(-vtac0(it))-exp(-vtac(it)))*exp(-sumct0-sumct)
     *-(vtac(it)-vtac0(it))*exp(-sumut))
c vv7: 0 targ. legs and >1 proj. legs (proj. diffr.)       
      vv7=max(0.d0,(1.d0-exp(-2.d0*vpac(ip)))*exp(-2.d0*sumcp)
     *-2.d0*vpac(ip)*exp(-2.d0*sumup))
     **((1.d0-exp(-vtac(it)))**2*exp(-2.d0*sumct)
     *+2.d0*(1.d0-exp(-vtac(it)))*(1.d0-exp(-sumct))*exp(-sumct))     
c vv8: 0 targ. legs and 1 proj. leg (proj. diffr.)       
      vv8=((vpac0(ip)+v1p0)*exp(-2.d0*sumup)
     *-(vpac(ip)-vpac0(ip)+v1p-v1p0)*exp(-sumup)*(1.d0-exp(-sumup)))
     **((1.d0-exp(-vtac(it)))**2*exp(-2.d0*sumct)
     *+2.d0*(1.d0-exp(-vtac(it)))*(1.d0-exp(-sumct))*exp(-sumct))     
c vv9: 0 targ. legs and >1 proj. legs with "holes" (proj. diffr.)       
      vv9=2.d0*(1.d0-exp(-2.d0*vtac(it)))*exp(-2.d0*sumct)
     **max(0.d0,(exp(-vpac0(ip))-exp(-vpac(ip)))*exp(-sumcp0-sumcp)
     *-(vpac(ip)-vpac0(ip))*exp(-sumup))
     
      vvt=vv1+vv2+vv3+vv4+vv5+vv6+vv7+vv8+vv9
      if(vvt.gt.0.d0)then
       aks=psran(b10)*vvt
       if(aks.lt.vv1)then
        jt=1
       elseif(aks.lt.vv1+vv2)then
        jt=2
       elseif(aks.lt.vv1+vv2+vv3)then
        jt=3
       elseif(aks.lt.vv1+vv2+vv3+vv4)then
        jt=4
       elseif(aks.lt.vv1+vv2+vv3+vv4+vv5)then
        jt=5
       elseif(aks.lt.vv1+vv2+vv3+vv4+vv5+vv6)then
        jt=6
       elseif(aks.lt.vv1+vv2+vv3+vv4+vv5+vv6+vv7)then
        jt=7
       elseif(aks.lt.vv1+vv2+vv3+vv4+vv5+vv6+vv7+vv8)then
        jt=8
       else
        jt=9
       endif
      elseif(vv1.gt.0.d0)then
       jt=1
      elseif(vv2.gt.0.d0)then
       jt=2
      elseif(vv3.gt.0.d0)then
       jt=3
      elseif(vv4.gt.0.d0)then
       jt=4
      elseif(vv5.gt.0.d0)then
       jt=5
      elseif(vv6.gt.0.d0)then
       jt=6
      elseif(vv7.gt.0.d0)then
       jt=7
      elseif(vv8.gt.0.d0)then
       jt=8
      elseif(vv9.gt.0.d0)then
       jt=9
      else
       iret=1
       goto 21
      endif
       
2     continue
      if(jt.eq.1.or.jt.eq.2.or.jt.eq.7)then
       ntry=0
3      ntry=ntry+1
       npprh0=0
       if(ip.eq.ia(1).or.ntry.gt.100)then
        nppr0=npgen(2.d0*vpac(ip),2,20)
        do i=1,nppr0
	 if(psran(b10).lt.vpac0(ip)/vpac(ip))then
	  itypr0(i)=0
	 else
          npprh0=npprh0+1
	  itypr0(i)=1
	 endif
	 ippr0(i)=ip
	enddo	 
       else
        nppr0=npgen(2.d0*vpac(ip),1,20)
        do i=1,nppr0
	 if(psran(b10).lt.vpac0(ip)/vpac(ip))then
	  itypr0(i)=0
	 else
          npprh0=npprh0+1
	  itypr0(i)=1
	 endif
	 ippr0(i)=ip
	enddo
	do ipp=ip+1,ia(1)
         ninc=npgen(2.d0*vpac(ipp),0,20)
	 if(ninc.ne.0)then
          nppr0=nppr0+ninc
          if(nppr0.gt.legmax)then
	   iret=1
	   goto 21
	  endif
	  do i=nppr0-ninc+1,nppr0
	   if(psran(b10).lt.vpac0(ipp)/vpac(ipp))then
	    itypr0(i)=0
	   else
            npprh0=npprh0+1
	    itypr0(i)=1
	   endif
	   ippr0(i)=ipp
	  enddo
	 endif
	enddo
	if(nppr0.eq.1)goto 3
       endif
      endif
	
      if(jt.eq.1.or.jt.eq.3.or.jt.eq.4)then
       ntry=0
4      ntry=ntry+1
       nptgh0=0
       if(it.eq.ia(2).or.ntry.gt.100)then
        nptg0=npgen(2.d0*vtac(it),2,20)
        do i=1,nptg0
	 if(psran(b10).lt.vtac0(it)/vtac(it))then
	  itytg0(i)=0
	 else
          nptgh0=nptgh0+1
	  itytg0(i)=1
	 endif
	 iptg0(i)=it
	enddo	 
       else
        nptg0=npgen(2.d0*vtac(it),1,20)
        do i=1,nptg0
	 if(psran(b10).lt.vtac0(it)/vtac(it))then
	  itytg0(i)=0
	 else
          nptgh0=nptgh0+1
	  itytg0(i)=1
	 endif
	 iptg0(i)=it
	enddo
	do itt=it+1,ia(2)
         ninc=npgen(2.d0*vtac(itt),0,20)
	 if(ninc.ne.0)then
          nptg0=nptg0+ninc
          if(nptg0.gt.legmax)then
	   iret=1
	   goto 21
	  endif
	  do i=nptg0-ninc+1,nptg0
	   if(psran(b10).lt.vtac0(itt)/vtac(itt))then
	    itytg0(i)=0
	   else
            nptgh0=nptgh0+1
	    itytg0(i)=1
	   endif
	   iptg0(i)=itt
	  enddo
	 endif
	enddo
	if(nptg0.eq.1)goto 4
       endif
      endif
      
      if(jt.eq.9)then
       ntry=0
5      ntry=ntry+1
       if(ip.eq.ia(1).or.ntry.gt.100)then
        nppr0=npgen(vpac(ip)-vpac0(ip),2,20)
        do i=1,nppr0
	 itypr0(i)=1
	 ippr0(i)=ip
	enddo	 
       else
        nppr0=npgen(vpac(ip)-vpac0(ip),1,20)
        do i=1,nppr0
	 itypr0(i)=1
	 ippr0(i)=ip
	enddo
	do ipp=ip+1,ia(1)
         ninc=npgen(vpac(ipp)-vpac0(ipp),0,20)
	 if(ninc.ne.0)then
          nppr0=nppr0+ninc
          if(nppr0.gt.legmax)then
	   iret=1
	   goto 21
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
	
      if(jt.eq.6)then
       ntry=0
6      ntry=ntry+1
       if(it.eq.ia(2).or.ntry.gt.100)then
        nptg0=npgen(vtac(it)-vtac0(it),2,20)
        do i=1,nptg0
	 itytg0(i)=1
	 iptg0(i)=it
	enddo	 
       else
        nptg0=npgen(vtac(it)-vtac0(it),1,20)
        do i=1,nptg0
	 itytg0(i)=1
	 iptg0(i)=it
	enddo
	do itt=it+1,ia(2)
         ninc=npgen(vtac(itt)-vtac0(itt),0,20)
	 if(ninc.ne.0)then
          nptg0=nptg0+ninc
          if(nptg0.gt.legmax)then
	   iret=1
	   goto 21
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
      if(jt.eq.1.and.nptgh0.lt.nptg0.and.npprh0.eq.nppr0)then
       gbt=1.d0-exp(sumup+(1.d0-nppr0)*dlog(2.d0))
      elseif(jt.eq.1.and.nptgh0.eq.nptg0.and.npprh0.lt.nppr0)then
       gbt=1.d0-exp(sumut+(1.d0-nptg0)*dlog(2.d0))
      elseif(jt.eq.1.and.nptgh0.eq.nptg0.and.npprh0.eq.nppr0)then
       gbt=1.d0-exp(sumut+(1.d0-nptg0)*dlog(2.d0))
     * -exp(sumup+(1.d0-nppr0)*dlog(2.d0))
      elseif(jt.eq.2.and.npprh0.eq.nppr0)then
       gbt=((vtac0(it)+v1t0)*exp(-2.d0*sumut)
     * -(vtac(it)-vtac0(it)+v1t-v1t0)*exp(-sumut)*(1.d0-exp(-sumut))
     * -exp(sumup+(1.d0-nppr0)*dlog(2.d0))
     * *(vtac(it)+v1t)*exp(-2.d0*sumut))
     * /((vtac0(it)+v1t0)*exp(-2.d0*sumut)
     * -(vtac(it)-vtac0(it)+v1t-v1t0)*exp(-sumut)*(1.d0-exp(-sumut)))
      elseif(jt.eq.3.and.nptgh0.eq.nptg0)then
       gbt=((vpac0(ip)+v1p0)*exp(-2.d0*sumup)
     * -(vpac(ip)-vpac0(ip)+v1p-v1p0)*exp(-sumup)*(1.d0-exp(-sumup))
     * -exp(sumut+(1.d0-nptg0)*dlog(2.d0))
     * *(vpac(ip)+v1p)*exp(-2.d0*sumup))
     * /((vpac0(ip)+v1p0)*exp(-2.d0*sumup)
     * -(vpac(ip)-vpac0(ip)+v1p-v1p0)*exp(-sumup)*(1.d0-exp(-sumup)))
      endif
      if(psran(b10).gt.gbt)goto 2

      if(jt.eq.2.or.jt.eq.5)then
       nptg0=1
       itytg0(1)=0
       iptg0(1)=it
      elseif(jt.eq.3.or.jt.eq.8)then
       nppr0=1
       itypr0(1)=0
       ippr0(1)=ip
      endif
      
      if(jt.eq.4.or.jt.eq.5.or.jt.eq.6)then
       nppr0=0
      elseif(jt.eq.7.or.jt.eq.8.or.jt.eq.9)then
       nptg0=0
      endif
      if(debug.ge.3)write (moniou,206)nppr0,nptg0
      
      nppr=0
      nptg=0      
      npin=0
      if((jt.eq.5.or.jt.eq.2.and.npprh0.lt.nppr0).and.psran(b10)
     *.lt.(vtac(it)-vtac0(it)+v1t-v1t0)*exp(-sumut)*(1.d0-exp(-sumut))
     */((vtac0(it)+v1t0)*exp(-2.d0*sumut)-(vtac(it)-vtac0(it)+v1t-v1t0)
     **exp(-sumut)*(1.d0-exp(-sumut)))
     *.or.jt.eq.2.and.npprh0.eq.nppr0.and.psran(b10)
     *.lt.(vtac(it)-vtac0(it)+v1t-v1t0)*(exp(-sumut)*(1.d0-exp(-sumut))
     **2.d0**nppr0*exp(-sumup)+2.d0*exp(-2.d0*sumut))
     */(((vtac0(it)+v1t0)*exp(-2.d0*sumut)-(vtac(it)-vtac0(it)+v1t-v1t0)
     **exp(-sumut)*(1.d0-exp(-sumut)))*2.d0**nppr0*exp(-sumup)
     *-2.d0*(vtac(it)+v1t)*exp(-2.d0*sumut)))then
       vvx=max(0.d0,1.d0+exp(-2.d0*(sumup+sumut))
     * -exp(-2.d0*sumup-sumut)-exp(-2.d0*sumut-sumup))
7      xpomr1=xpomr/(xpomr*scm)**psran(b10)
       v1i=((xpomr/xpomr1)**(dels-alf*vvx)
     * +qgpin0(xpomr/xpomr1,vvx)/sigs)*.7d0
       if(psran(b10).gt.v1i)goto 7
       npin=npin+1
       if(npin.gt.npmax)then
	iret=1
	goto 21
       endif
       xpomim(npin)=1.d0/xpomr1/scm
       xpomip(npin)=xpomr		
      endif
      
      if((jt.eq.8.or.jt.eq.3.and.nptgh0.lt.nptg0).and.psran(b10)
     *.lt.(vpac(ip)-vpac0(ip)+v1p-v1p0)*exp(-sumup)*(1.d0-exp(-sumup))
     */((vpac0(ip)+v1p0)*exp(-2.d0*sumup)-(vpac(ip)-vpac0(ip)+v1p-v1p0)
     **exp(-sumup)*(1.d0-exp(-sumup)))
     *.or.jt.eq.3.and.nptgh0.eq.nptg0.and.psran(b10)
     *.lt.(vpac(ip)-vpac0(ip)+v1p-v1p0)*(exp(-sumup)*(1.d0-exp(-sumup))
     **2.d0**nptg0*exp(-sumut)+2.d0*exp(-2.d0*sumup))
     */(((vpac0(ip)+v1p0)*exp(-2.d0*sumup)-(vpac(ip)-vpac0(ip)+v1p-v1p0)
     **exp(-sumup)*(1.d0-exp(-sumup)))*2.d0**nptg0*exp(-sumut)
     *-2.d0*(vpac(ip)+v1p)*exp(-2.d0*sumup)))then
       vvx=max(0.d0,1.d0+exp(-2.d0*(sumup+sumut))
     * -exp(-2.d0*sumup-sumut)-exp(-2.d0*sumut-sumup))
8      xpomr1=xpomr**psran(b10)
       v1i=((xpomr1/xpomr)**(dels-alf*vvx)
     * +qgpin0(xpomr1/xpomr,vvx)/sigs)*.7d0
       if(psran(b10).gt.v1i)goto 8
       npin=npin+1
       if(npin.gt.npmax)then
	iret=1
	goto 21
       endif
       xpomim(npin)=1.d0/xpomr/scm
       xpomip(npin)=xpomr1
      endif

      gbipr=exp(-sumct)*((1.d0-exp(-2.d0*sumup))*(1.d0-exp(-vtac(it)))
     *+(exp(-vtac0(it))-exp(-vtac(it)))*exp(-sumct0)
     **(exp(-2.d0*sumup)-max(0.d0,sumctt-sumct0))+max(0.d0
     *,exp(-vtac(it))-exp(-vtact0(it)))*exp(-max(sumctt,sumct0))
     *-max(0.d0,vtact0(it)-vtac0(it))*exp(-sumct0-vtac0(it)))
     
      gbitg=exp(-sumcp)*((1.d0-exp(-2.d0*sumut))*(1.d0-exp(-vpac(ip)))
     *+(exp(-vpac0(ip))-exp(-vpac(ip)))*exp(-sumcp0)
     **(exp(-2.d0*sumut)-max(0.d0,sumcpt-sumcp0))+max(0.d0
     *,exp(-vpac(ip))-exp(-vpact0(ip)))*exp(-max(sumcpt,sumcp0))
     *-max(0.d0,vpact0(ip)-vpac0(ip))*exp(-sumcp0-vpac0(ip)))
     
      gbi=(gbipr*max(0.d0,vpact0(ip)-vpac0(ip))
     *+gbitg*max(0.d0,vtact0(it)-vtac0(it)))/4.d0
      if(gbi.lt.0.d0)then
       gbi=0.d0
      endif
cdh 2006.02.21
cdh   nprin=npgen(gbi/gbs,0,npmax)
cdh   ntgin=npgen(gbi/gbs,0,npmax)
      nprin=npgen(gbi/gbs,0,50)
      ntgin=npgen(gbi/gbs,0,50)
      
      if(npin+nprin+ntgin.gt.npmax)then
	iret=1
	goto 21
      endif

      if(nprin.ne.0.or.ntgin.ne.0)then
       vvxt0=1.d0-exp(-sumcttf)
       gb0pr=exp(-sumcp)*((1.d0-exp(-2.d0*sumut))*(1.d0-exp(-vpac(ip)))
     * +exp(-2.d0*sumut)*(exp(-vpac0(ip))-exp(-vpac(ip)))*exp(-sumcp0)
     * +max(0.d0,exp(-vpac(ip))-exp(-vpact0(ip)))*exp(-sumcpt))
     * -max(0.d0,vpact0(ip)-vpac0(ip))*(1.d0-vvxt0)
       vvxp0=1.d0-exp(-sumcptf)
       gb0tg=exp(-sumct)*((1.d0-exp(-2.d0*sumup))*(1.d0-exp(-vtac(it)))
     * +exp(-2.d0*sumup)*(exp(-vtac0(it))-exp(-vtac(it)))*exp(-sumct0)
     * +max(0.d0,exp(-vtac(it))-exp(-vtact0(it)))*exp(-sumctt))
     * -max(0.d0,vtact0(it)-vtac0(it))*(1.d0-vvxp0)
       gb0=gb0pr*gbipr+gb0tg*gbitg
       if(gb0.le.0.d0)then
	nprin=0
	ntgin=0
       endif
      endif

      if(nprin.ne.0)then
       gbcor10=1.d0-(1.d0-(1.d0-xpomr)**(1.d0+ahl(icz)))
     * **(1.d0+dels-alpd)
       do in=1,nprin
	n=0
22      xpomr1=xpomr**psran(b10)
        gbcor1=1.d0-(1.d0-(1.d0-xpomr1)**(1.d0+ahl(icz)))
     *  **(1.d0+dels-alpd)
        if(psran(b10).gt.gbcor1/gbcor10)goto 22
	 
	call qgfdf(xxp,yyp,xpomr1,vpac,vtac,sumup,sumut,sumcp,sumct
     *  ,ip,it,ia(1),ia(2))
        sumcp0=0.d0
        sumcpt=0.d0
        sumcptf=0.d0
        do i=1,ia(1)
         ipp=ia(1)-i+1
         bbp=(xa(ipp,1)+b-xxp)**2+(xa(ipp,2)-yyp)**2
         if(ipp.ge.ip)vpac0(ipp)=min(vpac(ipp),qgfani(1.d0/xpomr1,bbp
     *   ,1.d0-exp(-sumut-sumcp),1.d0-exp(-sumcp0),1.d0-exp(-sumut)
     *   ,iddp(ipp),icz,5))
         vpact0(ipp)=qgfani(1.d0/xpomr1,bbp,1.d0-exp(-sumut-sumcp)
     *   ,1.d0-exp(-sumcptf),1.d0-exp(-sumcp),iddp(ipp),icz,6)
         sumcptf=sumcptf+vpact0(ipp)
         if(ipp.gt.ip)then
          sumcp0=sumcp0+vpac0(ipp)
          sumcpt=sumcpt+vpact0(ipp)
         endif
        enddo
        sumct0=0.d0
        sumctt=0.d0
        sumcttf=0.d0
        do i=1,ia(2)
         itt=ia(2)-i+1
         bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
         if(itt.ge.it)vtac0(itt)=min(vtac(itt),qgfani(xpomr1*scm,bbt
     *   ,1.d0-exp(-sumup-sumct),1.d0-exp(-sumct0),1.d0-exp(-sumup)
     *   ,iddt(itt),2,5))
         vtact0(itt)=qgfani(xpomr1*scm,bbt,1.d0-exp(-sumup-sumct)
     *   ,1.d0-exp(-sumcttf),1.d0-exp(-sumct),iddt(itt),2,6)
         sumcttf=sumcttf+vtact0(itt)
         if(itt.gt.it)then
          sumct0=sumct0+vtac0(itt)
          sumctt=sumctt+vtact0(itt)
         endif
        enddo
	
        vvxt=1.d0-dsqrt((1.d0-vvxt0)*exp(-sumcttf))
        vvxp=1.d0-dsqrt((1.d0-vvxp0)*exp(-sumcptf))
        v1ip=(xpomr1/xpomr)**(dels-alf*vvxt)
     *  +qgpin0(xpomr1/xpomr,vvxt)/sigs
	v1it=(xpomr1/xpomr)**(dels-alf*vvxp)
     *  +qgpin0(xpomr1/xpomr,vvxp)/sigs
    
	gbpr=v1ip*gbipr
     *  *(exp(-sumcp)*((1.d0-exp(-2.d0*sumut))*(1.d0-exp(-vpac(ip)))
     *  +exp(-2.d0*sumut)*(exp(-vpac0(ip))-exp(-vpac(ip)))*exp(-sumcp0)
     *  +max(0.d0,exp(-vpac(ip))-exp(-vpact0(ip)))*exp(-sumcpt))
     *  -max(0.d0,vpact0(ip)-vpac0(ip))*(1.d0-vvxt))
     *  +v1it*gb0tg
     *  *exp(-sumcp)*((1.d0-exp(-2.d0*sumut))*(1.d0-exp(-vpac(ip)))
     *  +(exp(-vpac0(ip))-exp(-vpac(ip)))*exp(-sumcp0)
     *  *(exp(-2.d0*sumut)-max(0.d0,sumcpt-sumcp0))+max(0.d0
     *  ,exp(-vpac(ip))-exp(-vpact0(ip)))*exp(-max(sumcpt,sumcp0))
     *  -max(0.d0,vpact0(ip)-vpac0(ip))*exp(-sumcp0-vpac0(ip)))
        gbpr=gbpr/gb0/gbcor1*gbcor10    *.05
        if(psran(b10).gt.gbpr)then
	 n=n+1
	 if(n.gt.1000)then
	  goto 23
	 endif
	 goto 22
        endif
        npin=npin+1
        xpomim(npin)=1.d0/xpomr/scm
        xpomip(npin)=xpomr1
       enddo
      endif
23    continue

      if(ntgin.ne.0)then
       gbcor20=1.d0-(1.d0-(1.d0-1.d0/xpomr/scm)
     * **(1.d0+ahl(2)))**(1.d0+dels-alpd)
       do in=1,ntgin
	n=0
24      xpomr1=xpomr/(xpomr*scm)**psran(b10)
        gbcor2=1.d0-(1.d0-(1.d0-1.d0/xpomr1/scm)
     *  **(1.d0+ahl(2)))**(1.d0+dels-alpd)
        if(psran(b10).gt.gbcor2/gbcor20)goto 24
	 
        call qgfdf(xxp,yyp,xpomr1,vpac,vtac,sumup,sumut,sumcp,sumct
     *  ,ip,it,ia(1),ia(2))
        sumcp0=0.d0
        sumcpt=0.d0
        sumcptf=0.d0
        do i=1,ia(1)
         ipp=ia(1)-i+1
         bbp=(xa(ipp,1)+b-xxp)**2+(xa(ipp,2)-yyp)**2
         if(ipp.ge.ip)vpac0(ipp)=min(vpac(ipp),qgfani(1.d0/xpomr1,bbp
     *   ,1.d0-exp(-sumut-sumcp),1.d0-exp(-sumcp0),1.d0-exp(-sumut)
     *   ,iddp(ipp),icz,5))
         vpact0(ipp)=qgfani(1.d0/xpomr1,bbp,1.d0-exp(-sumut-sumcp)
     *   ,1.d0-exp(-sumcptf),1.d0-exp(-sumcp),iddp(ipp),icz,6)
         sumcptf=sumcptf+vpact0(ipp)
         if(ipp.gt.ip)then
          sumcp0=sumcp0+vpac0(ipp)
          sumcpt=sumcpt+vpact0(ipp)
         endif
        enddo
        sumct0=0.d0
        sumctt=0.d0
        sumcttf=0.d0
        do i=1,ia(2)
         itt=ia(2)-i+1
         bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
         if(itt.ge.it)vtac0(itt)=min(vtac(itt),qgfani(xpomr1*scm,bbt
     *   ,1.d0-exp(-sumup-sumct),1.d0-exp(-sumct0),1.d0-exp(-sumup)
     *   ,iddt(itt),2,5))
         vtact0(itt)=qgfani(xpomr1*scm,bbt,1.d0-exp(-sumup-sumct)
     *   ,1.d0-exp(-sumcttf),1.d0-exp(-sumct),iddt(itt),2,6)
         sumcttf=sumcttf+vtact0(itt)
         if(itt.gt.it)then
          sumct0=sumct0+vtac0(itt)
          sumctt=sumctt+vtact0(itt)
         endif
        enddo

        vvxt=1.d0-dsqrt((1.d0-vvxt0)*exp(-sumcttf))
        vvxp=1.d0-dsqrt((1.d0-vvxp0)*exp(-sumcptf))
        v1ip=(xpomr/xpomr1)**(dels-alf*vvxt)
     *  +qgpin0(xpomr/xpomr1,vvxt)/sigs
	v1it=(xpomr/xpomr1)**(dels-alf*vvxp)
     *  +qgpin0(xpomr/xpomr1,vvxp)/sigs
     
	gbtg=v1it*gbitg
     *  *(exp(-sumct)*((1.d0-exp(-2.d0*sumup))*(1.d0-exp(-vtac(it)))
     *  +exp(-2.d0*sumup)*(exp(-vtac0(it))-exp(-vtac(it)))*exp(-sumct0)
     *  +max(0.d0,exp(-vtac(it))-exp(-vtact0(it)))*exp(-sumctt))
     *  -max(0.d0,vtact0(it)-vtac0(it))*(1.d0-vvxp))
     *  +v1ip*gb0pr
     *  *exp(-sumct)*((1.d0-exp(-2.d0*sumup))*(1.d0-exp(-vtac(it)))
     *  +(exp(-vtac0(it))-exp(-vtac(it)))*exp(-sumct0)
     *  *(exp(-2.d0*sumup)-max(0.d0,sumctt-sumct0))+max(0.d0
     *  ,exp(-vtac(it))-exp(-vtact0(it)))*exp(-max(sumctt,sumct0))
     *  -max(0.d0,vtact0(it)-vtac0(it))*exp(-sumct0-vtac0(it)))
        gbtg=gbtg/gb0/gbcor2*gbcor20   *.05
	
        if(psran(b10).gt.gbtg)then
	 n=n+1
	 if(n.gt.1000)then
	  goto 25
	 endif
	 goto 24
	endif
        npin=npin+1
        xpomim(npin)=1.d0/xpomr1/scm
        xpomip(npin)=xpomr
       enddo
      endif    
25    continue     

      call qgfdf(xxp,yyp,xpomr,vpac,vtac,sumup,sumut,sumcp,sumct
     *,ip,it,ia(1),ia(2))
      v1p=min(vpac(ip),qgfani(1.d0/xpomr,bbpr,1.d0-exp(-sumut-sumcp)
     *,1.d0-exp(min(0.d0,sumcp+2.d0*vpac(ip)-2.d0*sumup))
     *,0.d0,iddp(ip),icz,2))
      v1p0=min(v1p,qgfani(1.d0/xpomr,bbpr,1.d0-exp(-sumut-sumcp)
     *,1.d0-exp(min(0.d0,sumcp+2.d0*vpac(ip)-2.d0*sumup))
     *,1.d0-exp(-sumut),iddp(ip),icz,4))
      sumcp0=0.d0
      do i=1,ia(1)-ip+1
       ipp=ia(1)-i+1
       bbp=(xa(ipp,1)+b-xxp)**2+(xa(ipp,2)-yyp)**2
       vpac0(ipp)=min(vpac(ipp),qgfani(1.d0/xpomr,bbp,1.d0-exp(-sumut
     * -sumcp),1.d0-exp(-sumcp0),1.d0-exp(-sumut),iddp(ipp),icz,5))
       if(ipp.gt.ip)sumcp0=sumcp0+vpac0(ipp)
      enddo

      if(nppr0.eq.1.and.((jt.eq.8.or.jt.eq.3.and.nptgh0.lt.nptg0)
     *.and.psran(b10).lt.(v1p0*exp(-2.d0*sumup)
     *-(v1p-v1p0)*exp(-sumup)*(1.d0-exp(-sumup)))
     */max(.001d0,(vpac0(ip)+v1p0)*exp(-2.d0*sumup)
     *-(vpac(ip)-vpac0(ip)+v1p-v1p0)*exp(-sumup)*(1.d0-exp(-sumup)))
     *.or.jt.eq.3.and.nptgh0.eq.nptg0.and.psran(b10)
     *.lt.((v1p0*exp(-2.d0*sumup)-(v1p-v1p0)*exp(-sumup)
     **(1.d0-exp(-sumup)))*2.d0**nptg0*exp(-2.d0*sumut)
     *-2.d0*v1p*exp(-2.d0*sumup-sumut))
     */max(.001d0,((vpac0(ip)+v1p0)*exp(-2.d0*sumup)
     *-(vpac(ip)-vpac0(ip)+v1p-v1p0)*exp(-sumup)*(1.d0-exp(-sumup)))
     **2.d0**nptg0*exp(-2.d0*sumut)
     *-2.d0*(vpac(ip)+v1p)*exp(-2.d0*sumup-sumut))))then     !just one proj. leg
       nppr=1                               !number of proj. legs in the process
       xpompi(1)=xpomr                      !LC momentum for the 3P-vertex
       ipompi(1)=ip                         !proj. index for the leg
       vvxpi(1)=1.d0-exp(-sumup-sumut)      !screening factor
       bpompi(1)=bbpr                       !b^2 for the leg Pomeron
       if(debug.ge.3)write (moniou,207)
       goto 15
       
      elseif(nppr0.ne.0)then
       m=0
       nppm(1)=nppr0
       xpomm(1)=xpomr
       do i=1,nppr0
	ippm(i,1)=ippr0(i)
        itypm(i,1)=itypr0(i)
       enddo
            
9      m=m+1                                 !next level multi-Pomeron vertex
       ii(m)=0
10     ii(m)=ii(m)+1                         !next cut fan in the vertex
       if(ii(m).gt.nppm(m))then              !all fans at the level considered
        m=m-1                                !one level down
	if(m.eq.0)goto 15                    !all proj. fans considered 
	goto 10
       endif 
       l=ii(m)
       ipp=ippm(l,m)                         !proj. index for the leg      
       itypom=itypm(l,m)      
       bbp=(xa(ipp,1)+b-xxp)**2+(xa(ipp,2)-yyp)**2      !b^2 for the leg    
       if(debug.ge.4)write (moniou,208)ii(m),m,ipp,bbp
       call qgfdf(xxp,yyp,xpomm(m),vpac,vtac,sumup,sumut,sumcp,sumct
     * ,ipp,it,ia(1),ia(2))
       v1p=min(vpac(ipp),qgfani(1.d0/xpomm(m),bbp,1.d0-exp(-sumut-sumcp)
     * ,1.d0-exp(min(0.d0,sumcp+2.d0*vpac(ipp)-2.d0*sumup))
     * ,0.d0,iddp(ipp),icz,2))
       v1p0=min(v1p,qgfani(1.d0/xpomm(m),bbp,1.d0-exp(-sumut-sumcp)
     * ,1.d0-exp(min(0.d0,sumcp+2.d0*vpac(ipp)-2.d0*sumup))
     * ,1.d0-exp(-sumut),iddp(ipp),icz,4))
       sumcp0=0.d0
       do i=1,ia(1)-ipp+1
        ip1=ia(1)-i+1
        bbp1=(xa(ip1,1)+b-xxp)**2+(xa(ip1,2)-yyp)**2
        vpac0(ip1)=min(vpac(ip1),qgfani(1.d0/xpomm(m),bbp1
     *  ,1.d0-exp(-sumut-sumcp),1.d0-exp(-sumcp0)
     *  ,1.d0-exp(-sumut),iddp(ip1),icz,5))
        if(ip1.gt.ipp)sumcp0=sumcp0+vpac0(ip1)
       enddo
       vvxt=1.d0-exp(-sumut)

       if(itypom.eq.0.and.psran(b10).lt.v1p0/vpac0(ipp)
     * .or.m.eq.levmax)then                                        !single leg
        nppr=nppr+1
        if(nppr.gt.legmax)then
	 iret=1
	 goto 21
        endif
        xpompi(nppr)=xpomm(m)
        vvxpi(nppr)=1.d0-exp(-sumup-sumut)
        ipompi(nppr)=ipp
	bpompi(nppr)=bbp
        if(debug.ge.4)write (moniou,209)nppr,ipp,bbp,xpompi(nppr)
     *  ,vvxpi(nppr)
        goto 10        

       else                                            !new multi-Pomeron vertex
        if(debug.ge.4)write (moniou,210)m
        if(itypom.eq.0)then
         v3pt0=(max(0.d0,(1.d0-exp(-vpac(ipp)))*exp(-sumcp)
     *   -vpac0(ipp)*exp(-2.d0*sumup)
     *   -(exp(-vpac0(ipp))-exp(-vpac(ipp)))*exp(-sumcp0-sumcp))
     *   +(v1p-v1p0)*exp(-sumup)*(1.d0-exp(-sumup)))*(1.d0-vvxt)**2
          vvx=1.d0-exp(-2.d0*(sumup+sumut))
	else
         v3pt0=(1.d0-exp(-vpac(ipp)))*exp(-sumcp)*vvxt*(1.d0-vvxt)     
     *   +max(0.d0,(exp(-vpac0(ipp))-exp(-vpac(ipp)))*exp(-sumcp0-sumcp)
     *   -(vpac(ipp)-vpac0(ipp))*exp(-sumup))*(1.d0-vvxt)**2
         vvx=1.d0-exp(-2.d0*sumut-sumup)
	endif
	if(v3pt0.le.0.d0)then
	 iret=1
	 goto 21
	endif
        gbcor10=1.d0-(1.d0-(1.d0-xpomm(m))**(1.d0+ahl(icz)))
     *  **(1.d0+dels-alpd)

11      xpomm(m+1)=xpomm(m)**psran(b10)
        gbcor1=1.d0-(1.d0-(1.d0-xpomm(m+1))**(1.d0+ahl(icz)))
     *  **(1.d0+dels-alpd)
        if(psran(b10).gt.gbcor1/gbcor10)goto 11

        call qgfdf(xxp,yyp,xpomm(m+1),vpac,vtac,sumup,sumut,sumcp,sumct
     *  ,ipp,it,ia(1),ia(2))
        v1p=min(vpac(ipp),qgfani(1.d0/xpomm(m+1),bbp
     *  ,1.d0-exp(-sumut-sumcp),1.d0-exp(min(0.d0,sumcp+2.d0*vpac(ipp)
     *  -2.d0*sumup)),0.d0,iddp(ipp),icz,2))
        v1p0=min(v1p,qgfani(1.d0/xpomm(m+1),bbp,1.d0-exp(-sumut-sumcp)
     *  ,1.d0-exp(min(0.d0,sumcp+2.d0*vpac(ipp)-2.d0*sumup))
     *  ,1.d0-exp(-sumut),iddp(ipp),icz,4))
        sumcp0=0.d0
        do i=1,ia(1)-ipp+1
         ip1=ia(1)-i+1
         bbp1=(xa(ip1,1)+b-xxp)**2+(xa(ip1,2)-yyp)**2
         vpac0(ip1)=min(vpac(ip1),qgfani(1.d0/xpomm(m+1),bbp1
     *   ,1.d0-exp(-sumut-sumcp),1.d0-exp(-sumcp0)
     *   ,1.d0-exp(-sumut),iddp(ip1),icz,5))
         if(ip1.gt.ipp)sumcp0=sumcp0+vpac0(ip1)
        enddo
        vvxt=1.d0-exp(-sumut)

	v1i=(xpomm(m+1)/xpomm(m))**(dels-alf*vvx)
     *  +qgpin0(xpomm(m+1)/xpomm(m),vvx)/sigs
	if(itypom.eq.0)then
         v3pt=(max(0.d0,(1.d0-exp(-vpac(ipp)))*exp(-sumcp)
     *   -vpac0(ipp)*exp(-2.d0*sumup)
     *   -(exp(-vpac0(ipp))-exp(-vpac(ipp)))*exp(-sumcp0-sumcp))
     *   +(v1p-v1p0)*exp(-sumup)*(1.d0-exp(-sumup)))*(1.d0-vvxt)**2
	else
         v3pt=(1.d0-exp(-vpac(ipp)))*exp(-sumcp)*vvxt*(1.d0-vvxt)     
     *   +max(0.d0,(exp(-vpac0(ipp))-exp(-vpac(ipp)))*exp(-sumcp0-sumcp)
     *   -(vpac(ipp)-vpac0(ipp))*exp(-sumup))*(1.d0-vvxt)**2
	endif
        gb3pd=v3pt*v1i/v3pt0/gbcor1*gbcor10 *.2d0     
        if(psran(b10).gt.gb3pd)goto 11	
        
        if(itypom.eq.0)then
         npin=npin+1
         if(npin.gt.npmax)then
	  iret=1
	  goto 21
	 endif
         xpomim(npin)=1.d0/xpomm(m)/scm
         xpomip(npin)=xpomm(m+1)
         if(debug.ge.4)write (moniou,211)npin,xpomip(npin),xpomim(npin)
        endif
	
	aks=psran(b10)*v3pt	
	v3pd=.5d0*(1.d0-exp(-vpac(ipp)))**2*exp(-2.d0*sumcp)      !diffr. cut probability
     *  +(1.d0-exp(-vpac(ipp)))*exp(-sumcp)*(1.d0-exp(-sumcp))
        if(itypom.eq.0)then
         v3pd=v3pd*(1.d0-vvxt)**2
	else
         v3pd=v3pd*vvxt*(1.d0-vvxt)
	endif
	if(aks.lt.v3pd)then
         if(debug.ge.4)write (moniou,212)
	 goto 10
	endif
	
	if(itypom.eq.1)then
	 v3p1=v1p0*exp(-2.d0*sumup)*vvxt*(1.d0-vvxt)
	 if(aks.lt.v3pd+v3p1)then
          nppr=nppr+1
          if(nppr.gt.legmax)then
	   iret=1
	   goto 21
          endif
          xpompi(nppr)=xpomm(m+1)
          vvxpi(nppr)=1.d0-exp(-sumup-sumut)
          ipompi(nppr)=ipp
	  bpompi(nppr)=bbp
          goto 10 
	 endif
	       
	 v3ph=max(0.d0,(exp(-vpac0(ipp))-exp(-vpac(ipp)))
     *   *exp(-sumcp0-sumcp)-(vpac(ipp)-vpac0(ipp))*exp(-sumup))
     *   *(1.d0-vvxt)**2
         if(aks.lt.v3pd+v3p1+v3ph)then
	  ntry=0
12	  ntry=ntry+1
	  if(ipp.eq.ia(1).or.ntry.gt.100)then
           nppm(m+1)=npgen(vpac(ipp)-vpac0(ipp),2,20)
           do i=1,nppm(m+1)
	    itypm(i,m+1)=1
	    ippm(i,m+1)=ipp
	   enddo
	  else
           nppm(m+1)=npgen(vpac(ipp)-vpac0(ipp),1,20)
           do i=1,nppm(m+1)
	    itypm(i,m+1)=1
	    ippm(i,m+1)=ipp
	   enddo
	   do ip1=ipp+1,ia(1)
            ninc=npgen(vpac(ip1)-vpac0(ip1),0,20)
	    if(ninc.ne.0)then
             nppm(m+1)=nppm(m+1)+ninc
             if(nppm(m+1).gt.legmax)then
	      iret=1
	      goto 21
	     endif
	     do i=nppm(m+1)-ninc+1,nppm(m+1)
	      itypm(i,m+1)=1
	      ippm(i,m+1)=ip1
	     enddo
	    endif
	   enddo
	   if(nppm(m+1).eq.1)goto 12
	  endif
	  goto 9
         endif
        endif
	
        ntry=0
14      ntry=ntry+1
        npprh=0
	if(ipp.eq.ia(1).or.ntry.gt.100)then
         nppm(m+1)=npgen(2.d0*vpac(ipp),2,20)
         do i=1,nppm(m+1)
	  if(psran(b10).lt.vpac0(ipp)/vpac(ipp))then
	   itypm(i,m+1)=0
	  else
           npprh=npprh+1
	   itypm(i,m+1)=1
	  endif
	  ippm(i,m+1)=ipp
	 enddo
	else
         nppm(m+1)=npgen(2.d0*vpac(ipp),1,20)
         do i=1,nppm(m+1)
	  if(psran(b10).lt.vpac0(ipp)/vpac(ipp))then
	   itypm(i,m+1)=0
	  else
           npprh=npprh+1
	   itypm(i,m+1)=1
	  endif
	  ippm(i,m+1)=ipp
	 enddo
	 do ip1=ipp+1,ia(1)
          ninc=npgen(2.d0*vpac(ip1),0,20)
	  if(ninc.ne.0)then
           nppm(m+1)=nppm(m+1)+ninc
           if(nppm(m+1).gt.legmax)then
	    iret=1
	    goto 21
	   endif
	   do i=nppm(m+1)-ninc+1,nppm(m+1)
	    if(psran(b10).lt.vpac0(ip1)/vpac(ip1))then
	     itypm(i,m+1)=0
	    else
             npprh=npprh+1
	     itypm(i,m+1)=1
	    endif
	    ippm(i,m+1)=ip1
	   enddo
	  endif
	 enddo
	 if(nppm(m+1).eq.1)goto 14
	endif
	if(itypom.eq.0.and.npprh.eq.nppm(m+1).and.psran(b10)
     *  .gt.1.d0-exp(sumup+(1.d0-npprh)*dlog(2.d0)))goto 14
        if(debug.ge.4)write (moniou,213)nppm(m+1)
	goto 9
       endif
      endif

15    continue
      if(debug.ge.3)write (moniou,214)nppr
      call qgfdf(xxp,yyp,xpomr,vpac,vtac,sumup,sumut,sumcp,sumct
     *,ip,it,ia(1),ia(2))
      v1t=min(vtac(it),qgfani(xpomr*scm,bbtg,1.d0-exp(-sumup-sumct)
     *,1.d0-exp(min(0.d0,sumct+2.d0*vtac(it)-2.d0*sumut))
     *,0.d0,iddt(it),2,2))
      v1t0=min(v1t,qgfani(xpomr*scm,bbtg,1.d0-exp(-sumup-sumct)
     *,1.d0-exp(min(0.d0,sumct+2.d0*vtac(it)-2.d0*sumut))
     *,1.d0-exp(-sumup),iddt(it),2,4))
      sumct0=0.d0
      do i=1,ia(2)-it+1
       itt=ia(2)-i+1
       bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2
       vtac0(itt)=min(vtac(itt),qgfani(xpomr*scm,bbt,1.d0-exp(-sumup
     * -sumct),1.d0-exp(-sumct0),1.d0-exp(-sumup),iddt(itt),2,5))
       if(itt.gt.it)sumct0=sumct0+vtac0(itt)
      enddo
      
      if(nptg0.eq.1.and.((jt.eq.5.or.jt.eq.2.and.npprh0.lt.nppr0)
     *.and.psran(b10).lt.(v1t0*exp(-2.d0*sumut)
     *-(v1t-v1t0)*exp(-sumut)*(1.d0-exp(-sumut)))
     */max(.001d0,(vtac0(it)+v1t0)*exp(-2.d0*sumut)
     *-(vtac(it)-vtac0(it)+v1t-v1t0)*exp(-sumut)*(1.d0-exp(-sumut)))
     *.or.jt.eq.2.and.npprh0.eq.nppr0.and.psran(b10)
     *.lt.((v1t0*exp(-2.d0*sumut)-(v1t-v1t0)*exp(-sumut)
     **(1.d0-exp(-sumut)))*2.d0**nppr0*exp(-2.d0*sumup)
     *-2.d0*v1t*exp(-2.d0*sumut-sumup))
     */max(.001d0,((vtac0(it)+v1t0)*exp(-2.d0*sumut)
     *-(vtac(it)-vtac0(it)+v1t-v1t0)*exp(-sumut)*(1.d0-exp(-sumut)))
     **2.d0**nppr0*exp(-2.d0*sumup)
     *-2.d0*(vtac(it)+v1t)*exp(-2.d0*sumut-sumup))))then     !just one targ. leg
       nptg=1
       xpomti(1)=xpomr                            !LC momentum for the 3P-vertex
       ipomti(1)=it                               !targ. index for the leg
       vvxti(1)=1.d0-exp(-sumut-sumup)            !screening factor
       bpomti(1)=bbtg                             !b^2 for the leg Pomeron      
       if(debug.ge.3)write (moniou,215)
       goto 21
       
      elseif(nptg0.ne.0)then
       m=0
       nppm(1)=nptg0
       xpomm(1)=xpomr
       do i=1,nptg0
	ippm(i,1)=iptg0(i)
        itypm(i,1)=itytg0(i)
       enddo
                        
16     m=m+1                                   !next level multi-Pomeron vertex
       ii(m)=0
17     ii(m)=ii(m)+1                           !next cut fan in the vertex
       if(ii(m).gt.nppm(m))then                !all fans at the level considered
        m=m-1                                  !one level down
	if(m.eq.0)goto 21                      !all targ. fans considered
	goto 17
       endif 
       l=ii(m)
       itt=ippm(l,m)                           !targ. index for the leg      
       itypom=itypm(l,m)      
       bbt=(xb(itt,1)-xxp)**2+(xb(itt,2)-yyp)**2  !b^2 for the leg  
       if(debug.ge.4)write (moniou,216)ii(m),m,itt,bbt
       
       call qgfdf(xxp,yyp,xpomm(m),vpac,vtac,sumup,sumut,sumcp,sumct
     * ,ip,itt,ia(1),ia(2))
       v1t=min(vtac(itt),qgfani(xpomm(m)*scm,bbt,1.d0-exp(-sumup-sumct)
     * ,1.d0-exp(min(0.d0,sumct+2.d0*vtac(itt)-2.d0*sumut))
     * ,0.d0,iddt(itt),2,2))
       v1t0=min(v1t,qgfani(xpomm(m)*scm,bbt,1.d0-exp(-sumup-sumct)
     * ,1.d0-exp(min(0.d0,sumct+2.d0*vtac(itt)-2.d0*sumut))
     * ,1.d0-exp(-sumup),iddt(itt),2,4))
       sumct0=0.d0
       do i=1,ia(2)-itt+1
        it1=ia(2)-i+1
        bbt1=(xb(it1,1)-xxp)**2+(xb(it1,2)-yyp)**2
        vtac0(it1)=min(vtac(it1),qgfani(xpomm(m)*scm,bbt1
     *  ,1.d0-exp(-sumup-sumct),1.d0-exp(-sumct0),1.d0-exp(-sumup)
     *  ,iddt(it1),2,5))
        if(it1.gt.itt)sumct0=sumct0+vtac0(it1)
       enddo
       vvxp=1.d0-exp(-sumup)

       if(itypom.eq.0.and.psran(b10).lt.v1t0/vtac0(itt)            !single leg
     * .or.m.eq.levmax)then
        nptg=nptg+1
        if(nptg.gt.legmax)then
	 iret=1
	 goto 21
	endif
        xpomti(nptg)=xpomm(m)
        vvxti(nptg)=1.d0-exp(-sumut-sumup)
        ipomti(nptg)=itt
	bpomti(nptg)=bbt
        if(debug.ge.4)write (moniou,217)nptg,itt,bbt,xpomti(nptg)
     *  ,vvxti(nptg)
        goto 17        

       else                                            !new multi-Pomeron vertex	
        if(debug.ge.4)write (moniou,210)m
        if(itypom.eq.0)then
         v3pt0=(max(0.d0,(1.d0-exp(-vtac(itt)))*exp(-sumct)
     *   -vtac0(itt)*exp(-2.d0*sumut)
     *   -(exp(-vtac0(itt))-exp(-vtac(itt)))*exp(-sumct0-sumct))
     *   +(v1t-v1t0)*exp(-sumut)*(1.d0-exp(-sumut)))*(1.d0-vvxp)**2
          vvx=1.d0-exp(-2.d0*(sumup+sumut))
	else
         v3pt0=(1.d0-exp(-vtac(itt)))*exp(-sumct)*vvxp*(1.d0-vvxp)     
     *   +max(0.d0,(exp(-vtac0(itt))-exp(-vtac(itt)))*exp(-sumct0-sumct)
     *   -(vtac(itt)-vtac0(itt))*exp(-sumut))*(1.d0-vvxp)**2
         vvx=1.d0-exp(-2.d0*sumup-sumut)
	endif
        gbcor20=1.d0-(1.d0-(1.d0-1.d0/xpomm(m)/scm)
     *  **(1.d0+ahl(2)))**(1.d0+dels-alpd)

18      xpomm(m+1)=xpomm(m)/(xpomm(m)*scm)**psran(b10)
        gbcor2=1.d0-(1.d0-(1.d0-1.d0/xpomm(m+1)/scm)
     *  **(1.d0+ahl(2)))**(1.d0+dels-alpd)
        if(psran(b10).gt.gbcor2/gbcor20)goto 18

        call qgfdf(xxp,yyp,xpomm(m+1),vpac,vtac,sumup,sumut,sumcp,sumct
     *  ,ip,itt,ia(1),ia(2))
        v1t=min(vtac(itt),qgfani(xpomm(m+1)*scm,bbt
     *  ,1.d0-exp(-sumup-sumct),1.d0-exp(min(0.d0,sumct+2.d0*vtac(itt)
     *  -2.d0*sumut)),0.d0,iddt(itt),2,2))
        v1t0=min(v1t,qgfani(xpomm(m+1)*scm,bbt,1.d0-exp(-sumup-sumct)
     *  ,1.d0-exp(min(0.d0,sumct+2.d0*vtac(itt)-2.d0*sumut))
     *  ,1.d0-exp(-sumup),iddt(itt),2,4))
        sumct0=0.d0
        do i=1,ia(2)-itt+1
         it1=ia(2)-i+1
         bbt1=(xb(it1,1)-xxp)**2+(xb(it1,2)-yyp)**2
         vtac0(it1)=min(vtac(it1),qgfani(xpomm(m+1)*scm,bbt1
     *   ,1.d0-exp(-sumup-sumct),1.d0-exp(-sumct0),1.d0-exp(-sumup)
     *   ,iddt(it1),2,5))
         if(it1.gt.itt)sumct0=sumct0+vtac0(it1)
        enddo
        vvxp=1.d0-exp(-sumup)       
	v1i=(xpomm(m)/xpomm(m+1))**(dels-alf*vvx)
     *  +qgpin0(xpomm(m)/xpomm(m+1),vvx)/sigs
     
	if(itypom.eq.0)then
         v3pt=(max(0.d0,(1.d0-exp(-vtac(itt)))*exp(-sumct)
     *   -vtac0(itt)*exp(-2.d0*sumut)
     *   -(exp(-vtac0(itt))-exp(-vtac(itt)))*exp(-sumct0-sumct))
     *   +(v1t-v1t0)*exp(-sumut)*(1.d0-exp(-sumut)))*(1.d0-vvxp)**2
	else
         v3pt=(1.d0-exp(-vtac(itt)))*exp(-sumct)*vvxp*(1.d0-vvxp)
     *   +max(0.d0,(exp(-vtac0(itt))-exp(-vtac(itt)))*exp(-sumct0-sumct)
     *   -(vtac(itt)-vtac0(itt))*exp(-sumut))*(1.d0-vvxp)**2
	endif
        gb3pd=v3pt*v1i/v3pt0/gbcor2*gbcor20   *.2d0
        if(psran(b10).gt.gb3pd)goto 18
	
        if(itypom.eq.0)then
         npin=npin+1
         if(npin.gt.npmax)then
	  iret=1
	  goto 21
	 endif
         xpomim(npin)=1.d0/xpomm(m+1)/scm
         xpomip(npin)=xpomm(m)
         if(debug.ge.4)write (moniou,211)npin,xpomip(npin),xpomim(npin)
        endif
		
	aks=psran(b10)*v3pt
        v3pd=.5d0*(1.d0-exp(-vtac(itt)))**2*exp(-2.d0*sumct)     !diffr. cut probability
     *  +(1.d0-exp(-vtac(itt)))*exp(-sumct)*(1.d0-exp(-sumct))
        if(itypom.eq.0)then
         v3pd=v3pd*(1.d0-vvxp)**2
	else
         v3pd=v3pd*vvxp*(1.d0-vvxp)
	endif
	if(aks.lt.v3pd)then
         if(debug.ge.4)write (moniou,212)
	 goto 17
        endif

	if(itypom.eq.1)then
	 v3p1=v1t0*exp(-2.d0*sumut)*vvxp*(1.d0-vvxp)
	 if(aks.lt.v3pd+v3p1)then
          nptg=nptg+1
          if(nptg.gt.legmax)then
	   iret=1
	   goto 21
	  endif
          xpomti(nptg)=xpomm(m+1)
          vvxti(nptg)=1.d0-exp(-sumut-sumup)
          ipomti(nptg)=itt
	  bpomti(nptg)=bbt
          goto 17
	 endif
	       
	 v3ph=max(0.d0,(exp(-vtac0(itt))-exp(-vtac(itt)))
     *   *exp(-sumct0-sumct)-(vtac(itt)-vtac0(itt))*exp(-sumut))
     *   *(1.d0-vvxp)**2
	 if(aks.lt.v3pd+v3p1+v3ph)then
	  ntry=0
19	  ntry=ntry+1
	  if(itt.eq.ia(2).or.ntry.gt.100)then
           nppm(m+1)=npgen(vtac(itt)-vtac0(itt),2,20)
           do i=1,nppm(m+1)
	    itypm(i,m+1)=1
	   ippm(i,m+1)=itt
	   enddo
	  else
           nppm(m+1)=npgen(vtac(itt)-vtac0(itt),1,20)
           do i=1,nppm(m+1)
	    itypm(i,m+1)=1
	   ippm(i,m+1)=itt
	   enddo
	   do it1=itt+1,ia(2)
            ninc=npgen(vtac(it1)-vtac0(it1),0,20)
	    if(ninc.ne.0)then
             nppm(m+1)=nppm(m+1)+ninc
             if(nppm(m+1).gt.legmax)then
	      iret=1
	      goto 21
	     endif
	     do i=nppm(m+1)-ninc+1,nppm(m+1)
	      itypm(i,m+1)=1
	      ippm(i,m+1)=it1
	     enddo
	    endif
	   enddo
	   if(nppm(m+1).eq.1)goto 19
	  endif
	  goto 16
         endif
        endif
	
	ntry=0
20	ntry=ntry+1
        nptgh=0
	if(itt.eq.ia(2).or.ntry.gt.100)then
         nppm(m+1)=npgen(2.d0*vtac(itt),2,20)
	 do i=1,nppm(m+1)
	  if(psran(b10).lt.vtac0(itt)/vtac(itt))then
	   itypm(i,m+1)=0
	  else
           nptgh=nptgh+1
	   itypm(i,m+1)=1
	  endif
	  ippm(i,m+1)=itt
	 enddo
	else
         nppm(m+1)=npgen(2.d0*vtac(itt),1,20)
         do i=1,nppm(m+1)
	  if(psran(b10).lt.vtac0(itt)/vtac(itt))then
	   itypm(i,m+1)=0
	  else
           nptgh=nptgh+1
	   itypm(i,m+1)=1
	  endif
	  ippm(i,m+1)=itt
	 enddo
	 do it1=itt+1,ia(2)
          ninc=npgen(2.d0*vtac(it1),0,20)
	  if(ninc.ne.0)then
           nppm(m+1)=nppm(m+1)+ninc
           if(nppm(m+1).gt.legmax)then
	    iret=1
	    goto 21
	   endif
	   do i=nppm(m+1)-ninc+1,nppm(m+1)
	    if(psran(b10).lt.vtac0(it1)/vtac(it1))then
	     itypm(i,m+1)=0
	    else
             nptgh=nptgh+1
	     itypm(i,m+1)=1
	    endif
	    ippm(i,m+1)=it1
	   enddo
	  endif
	 enddo
	 if(nppm(m+1).eq.1)goto 20
	endif
	if(itypom.eq.0.and.nptgh.eq.nppm(m+1).and.psran(b10)
     *  .gt.1.d0-exp(sumut+(1.d0-nptgh)*dlog(2.d0)))goto 20
        if(debug.ge.4)write (moniou,218)nppm(m+1)
	goto 16
       endif
      endif
21    continue
      if(debug.ge.2)write (moniou,219)nppr,nptg,npin,iret

201   format(2x,'qg3pdf - configuration for multi-Pomeron'
     *,'/diffractive contributions'                              !so161205
     */4x,i2,'-th proj. nucleon',2x,i2,'-th targ. nucleon')      !so161205
202   format(2x,'qg3pdf: problem with initial normalization'
     *,' -> rejection')
203   format(2x,'qg3pdf: normalization of rejection function - ',e10.3)
204   format(2x,'qg3pdf: xpomr=',e10.3,2x,'bbpr=',e10.3,2x,'bbtg=',e10.3
     *,2x,'gb=',e10.3)
205   format(2x,'qg3pdf: xpomr=',e10.3,2x,'bbpr=',e10.3,2x,'bbtg=',e10.3
     *,2x,'xxp=',e10.3,2x,'yyp=',e10.3)
206   format(2x,'qg3pdf: main vertex, nppr0=',i3,2x,'nptg0=',i3)
207   format(2x,'qg3pdf: single proj. leg Pomeron')
208   format(2x,'qg3pdf: check',i3,'-th cut fan at ',i2,'-th level,'
     *,' proj. index - ',i3,2x,'b^2=',e10.3)
209   format(2x,'qg3pdf: ',i3,'-th proj. leg, proj. index - ',i3
     *,2x,'b^2=',e10.3,2x,'xpomr=',e10.3,2x,'vvx=',e10.3)
210   format(2x,'qg3pdf: new vertex at ',i3,'-th level')
211   format(2x,'qg3pdf: ',i3,'-th interm. Pomeron, xpomip=',e10.3
     *,2x,'xpomim=',e10.3)
212   format(2x,'qg3pdf: diffractive cut')
213   format(2x,'qg3pdf: new multi-Pomeron vertex with',i3,'proj. legs')
214   format(2x,'qg3pdf: total number of proj. legs - ',i3)
215   format(2x,'qg3pdf: single targ. leg Pomeron')
216   format(2x,'qg3pdf: check',i3,'-th cut fan at ',i2,'-th level,'
     *,' targ. index - ',i3,2x,'b^2=',e10.3)
217   format(2x,'qg3pdf: ',i3,'-th targ. leg, targ. index - ',i3
     *,2x,'b^2=',e10.3,2x,'xpomr=',e10.3,2x,'vvx=',e10.3)
218   format(2x,'qg3pdf: new multi-Pomeron vertex with',i3,'targ. legs')
219   format(2x,'qg3pdf - end',2x,'number of proj. legs:',i3
     *,2x,'number of targ. legs:',i3
     */4x,'number of interm. Pomerons:',i3,'return flag:',i2)
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
      common /debug/   debug
      external psran

      if(debug.ge.3)write (moniou,201)vv,npmin,npmax
      if(npmin.eq.0)then
       aks=psran(b10)
       vvn=exp(-vv)
       do n=1,npmax
 	aks=aks-vvn
        if(aks.lt.0.d0)goto 1
 	vvn=vvn*vv/dble(n)
       enddo
      elseif(npmin.eq.1)then
       aks=psran(b10)*(1.d0-exp(-vv))
       vvn=exp(-vv)
       do n=1,npmax
 	vvn=vvn*vv/dble(n)
 	aks=aks-vvn
        if(aks.lt.0.d0)goto 2
       enddo
      elseif(npmin.eq.2)then
       aks=psran(b10)*(1.d0-exp(-vv)*(1.d0+vv))
       vvn=vv*exp(-vv)
       do n=2,npmax
 	vvn=vvn*vv/dble(n)
 	aks=aks-vvn
        if(aks.lt.0.d0)goto 2
       enddo
      else
       stop'npgen - wrong input'
      endif
1     n=n-1
2     npgen=n
      if(debug.ge.4)write (moniou,202)n
      
201   format(2x,'npgen -  Poisson distribution',2x,'average=',e10.3
     *,2x,'min=',i1,2x,'max=',i3)
202   format(2x,'npgen=',e10.3)
      return
      end
                       
c=============================================================================
      subroutine qgsha(nbpom,ncola,ncolb)
c-----------------------------------------------------------------------------
c qgsha - inelastic interaction (energy sharing and particle production)
c nbpom - number of Pomeron blocks (nucleon(hadron)-nucleon collisions),
c ncola - number of inel.-wounded proj. nucleons,
c ncolb - number of inel.-wounded targ. nucleons
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207,npbmax=1000,npnmax=1000,npmax=5000
     *,legmax=900,njmax=50000)
      dimension wppr0(iapmax),wmtg0(iapmax),wppr1(iapmax),wmtg1(iapmax)
     *,wppr2(iapmax),wmtg2(iapmax),izp(iapmax),izt(iapmax)
     *,ila(iapmax),ilb(iapmax),lva(iapmax),lvb(iapmax)
     *,lqa0(iapmax),lqb0(iapmax),ncola(iapmax),ncolb(iapmax)
     *,ncola0(iapmax),ncolb0(iapmax)
     *,xpomp0(npnmax,iapmax),xpomt0(npnmax,iapmax)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr9/  iwp(iapmax),iwt(iapmax),lqa(iapmax),lqb(iapmax)
     *,iprcn(iapmax),itgcn(iapmax),ias(npbmax),ibs(npbmax),nqs(npbmax)
     *,npompr(npbmax),npomtg(npbmax),npomin(npbmax),nnpr(npmax,npbmax)
     *,nntg(npmax,npbmax),ilpr(legmax,npbmax),iltg(legmax,npbmax)
     *,lnpr(legmax,npbmax),lntg(legmax,npbmax)
     *,nbpi(npnmax,iapmax),nbti(npnmax,iapmax),idnpi(npnmax,iapmax)
     *,idnti(npnmax,iapmax),nppi(npnmax,iapmax),npti(npnmax,iapmax)
     *,nlpi(npnmax,iapmax),nlti(npnmax,iapmax)
      common /qgarr11/ b10
      common /qgarr12/ nsp
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr22/ wppr(iapmax),wmtg(iapmax)
      common /qgarr23/ bbpom(npbmax),vvxpom(npbmax)
     *,bpompr(npnmax,iapmax),bpomtg(npnmax,iapmax),vvxpr(npnmax,iapmax)
     *,vvxtg(npnmax,iapmax),xpompr(npnmax,iapmax),xpomtg(npnmax,iapmax)
     *,xpopin(npmax,npbmax),xpomin(npmax,npbmax)
      common /qgarr37/ eqj(4,njmax),iqj(njmax),ncj(2,njmax),nj
      common /qgarr40/ xppr(npnmax,iapmax),xmtg(npnmax,iapmax)
      common /qgarr43/ moniou
      common /debug/   debug
      external psran

      if(debug.ge.1)write (moniou,201)nbpom             !so161205
      nsp0=nsp
      
      do j=1,ia(1)
       if(lqa(j).ne.0)then
        do i=1,lqa(j)
	 if(idnpi(i,j).ne.0)xpomp0(i,j)=xpompr(i,j)
        enddo
       endif
      enddo
      do j=1,ia(2)
       if(lqb(j).ne.0)then
        do i=1,lqb(j)
	 if(idnti(i,j).ne.0)xpomt0(i,j)=xpomtg(i,j)
        enddo
       endif
      enddo
      
1     nsp=nsp0
      nj=0

      do j=1,ia(1)
       if(lqa(j).ne.0)then
        do i=1,lqa(j)
	 if(idnpi(i,j).ne.0)xpompr(i,j)=xpomp0(i,j)
        enddo
       endif
      enddo
      do j=1,ia(2)
       if(lqb(j).ne.0)then
        do i=1,lqb(j)
	 if(idnti(i,j).ne.0)xpomtg(i,j)=xpomt0(i,j)
        enddo
       endif
      enddo
      
c-------------------------------------------------
c initial nucleon (hadron) types
      if(ia(1).ne.1)then
       do i=1,ia(1)                
        izp(i)=int(2.5d0+psran(b10))   !i-th projectile nucleon type
       enddo
      else
       izp(1)=icp                      !projectile hadron type
      endif
      if(ia(2).ne.1)then
       do i=1,ia(2)
        izt(i)=int(2.5d0+psran(b10))   !i-th target nucleon type
       enddo
      else
       izt(1)=2                        !target proton
      endif
      
      do i=1,ia(1)
       lqa0(i)=lqa(i)
       lva(i)=0
       ncola0(i)=ncola(i)
      enddo
      do i=1,ia(2)
       lqb0(i)=lqb(i)
       lvb(i)=0
       ncolb0(i)=ncolb(i)
      enddo
      
c-------------------------------------------------
c energy-momentum sharing between Pomerons
      if(nbpom.ne.0)then
       if(debug.ge.1)write (moniou,202)
       call qgprox(0)        !initial x-configuration
       gbl0=qgweix(nbpom)    !log-weight for the initial x-configuration
       nrej=0
       nchange=0
       gbnorm=.1d0
cdh    gbhmax=-1000.d0           !so 26.04.06
       gbhmax=-10.d10

2      continue
       call qgprox(1)        !proposed x-configuration
       gbl=qgweix(nbpom)     !log-weight for the proposed x-configuration       
       gbh=gbl-gbl0-gbnorm   !log of acceptance probability
       gbhmax=max(gbhmax,gbh)

       if(debug.ge.5)write (moniou,203)gbh,nrej,nchange
       if(gbh.lt.-50.d0.or.psran(b10).gt.exp(gbh))then
        nrej=nrej+1
cdh     if(nrej.gt.100)then              !so 26.04.06
        if(nrej.gt.50)then               !too many rejections
	 nrej=0
	 nchange=nchange+1
	 gbnorm=gbnorm+gbhmax+.5d0        !new normalization of acceptance
cdh      gbhmax=-1000.d0                  !so 26.04.06
         gbhmax=-10.d10
         if(debug.ge.4)write (moniou,204)nchange
	endif
        goto 2                            !rejection
       endif
      endif
       
c-------------------------------------------------
c leading remnant LC momenta
      if(debug.ge.1)write (moniou,205)
      do i=1,ia(1)                        !loop over proj. nucleons
       wppr0(i)=wppr(i)
       wppr1(i)=0.d0
       wppr2(i)=0.d0
       if(lqa(i).ne.0)then
        do l=1,lqa(i)                     !loop over constituent partons
	 wppr0(i)=wppr0(i)-wp0*xppr(l,i)  !subtract Pomeron LC momentum
	 if(wppr0(i).le.0.d0)stop'w^+<0!!!'
	enddo
       endif
      enddo
      do i=1,ia(2)                        !loop over targ. nucleons
       wmtg0(i)=wmtg(i)
       wmtg1(i)=0.d0
       wmtg2(i)=0.d0
       if(lqb(i).ne.0)then
        do l=1,lqb(i)                     !loop over constituent partons
	 wmtg0(i)=wmtg0(i)-wm0*xmtg(l,i)  !subtract Pomeron LC momentum
	 if(wmtg0(i).le.0.d0)stop'w^-<0!!!'
	enddo
       endif
      enddo

c-------------------------------------------------
c momentum conservation (correction for 3p-vertexes)
      if(debug.ge.1)write (moniou,206)
      if(nbpom.ne.0)then
       do nb=1,nbpom                            !loop over collisions
        ip=ias(nb)                              !proj. index
	it=ibs(nb)                              !targ. index
	if(nqs(nb).ne.0)then
	 do np=1,nqs(nb)             !loop over single Pomerons in the collision
	  lnp=nnpr(np,nb)                       !proj. constituent parton index
	  lnt=nntg(np,nb)                       !targ. constituent parton index
          wppr1(ip)=wppr1(ip)+xppr(lnp,ip)*wp0  !count Pomeron LC momentum	
          wmtg1(it)=wmtg1(it)+xmtg(lnt,it)*wm0  !count Pomeron LC momentum
	 enddo
	endif		
	if(npomin(nb).ne.0)then
	 do np=1,npomin(nb)         !loop over interm. Pomerons in the collision
	  xpp=xpopin(np,nb)
          xpm=xpomin(np,nb)
	  if(xpp*xpm*scm.gt.1.d0)then
           wppr2(ip)=wppr2(ip)+xpp*wp0          !count Pomeron LC momentum	
           wmtg2(it)=wmtg2(it)+xpm*wm0          !count Pomeron LC momentum
	  else
	   xpopin(np,nb)=0.d0
	   xpomin(np,nb)=0.d0
	  endif
	 enddo
	endif 	 
	if(npompr(nb).ne.0)then
	 do np=1,npompr(nb)       !loop over proj. leg Pomerons in the collision
	  ipp=ilpr(np,nb)                       !proj. index
	  lnp=lnpr(np,nb)                       !proj. constituent parton index
	  xpp=xppr(lnp,ipp)
	  xpm=xpompr(lnp,ipp)
          if(xpp*xpm*scm.gt.1.d0)then
           wppr1(ipp)=wppr1(ipp)+xpp*wp0        !count Pomeron LC momentum	
           wmtg2(it)=wmtg2(it)+xpm*wm0          !count Pomeron LC momentum
	  else
	   xppr(lnp,ipp)=0.d0
	   xpompr(lnp,ipp)=0.d0
	  endif
	 enddo
	endif
	if(npomtg(nb).ne.0)then
	 do np=1,npomtg(nb)       !loop over targ. leg Pomerons in the collision
	  itt=iltg(np,nb)                       !targ. index
	  lnt=lntg(np,nb)                       !targ. constituent parton index
	  xpp=xpomtg(lnt,itt)
	  xpm=xmtg(lnt,itt)
	  if(xpp*xpm*scm.gt.1.d0)then
           wppr2(ip)=wppr2(ip)+xpp*wp0	        !count Pomeron LC momentum
           wmtg1(itt)=wmtg1(itt)+xpm*wm0        !count Pomeron LC momentum
	  else
	   xmtg(lnt,itt)=0.d0
	   xpomtg(lnt,itt)=0.d0
	  endif
	 enddo
	endif
       enddo
      endif

      do ip=1,ia(1)
       if(wppr1(ip)+wppr2(ip).ne.0.d0)then
        if(lqa(ip).ne.0)then
	 do i=1,lqa(ip)
	  xppr(i,ip)=xppr(i,ip)*(wppr(ip)-wppr0(ip)) !renorm. for const. partons
     *    /(wppr1(ip)+wppr2(ip))
         enddo

         do nb=1,nbpom
	  if(ias(nb).eq.ip.and.(npomtg(nb).ne.0.or.npomin(nb).ne.0))then
	   if(npomin(nb).ne.0)then
	    do np=1,npomin(nb)
	     xpopin(np,nb)=xpopin(np,nb)*(wppr(ip)-wppr0(ip)) 
     *       /(wppr1(ip)+wppr2(ip))
	    enddo
	   endif
	   if(npomtg(nb).ne.0)then
	    do np=1,npomtg(nb)
	     itt=iltg(np,nb)
	     lnt=lntg(np,nb)
	     xpomtg(lnt,itt)=xpomtg(lnt,itt)*(wppr(ip)-wppr0(ip))
     *       /(wppr1(ip)+wppr2(ip))
	    enddo
	   endif
	  endif
	 enddo
	 
        elseif(wppr2(ip).gt.wppr(ip))then
         wppr0(ip)=wppr(ip)/2.d0
         do nb=1,nbpom
	  if(ias(nb).eq.ip.and.(npomtg(nb).ne.0.or.npomin(nb).ne.0))then
	   if(npomin(nb).ne.0)then
	    do np=1,npomin(nb)
	     xpopin(np,nb)=xpopin(np,nb)*wppr(ip)/wppr2(ip)/2.d0
	    enddo
	   endif
	   if(npomtg(nb).ne.0)then
	    do np=1,npomtg(nb)
	     itt=iltg(np,nb)
	     lnt=lntg(np,nb)
	     xpomtg(lnt,itt)=xpomtg(lnt,itt)*wppr(ip)/wppr2(ip)/2.d0
	    enddo
	   endif
	  endif
	 enddo
	else
         wppr0(ip)=wppr(ip)-wppr2(ip)
	endif
       else
        wppr0(ip)=wp0
       endif
      enddo
	 
      do it=1,ia(2)
       if(wmtg1(it)+wmtg2(it).ne.0.d0)then
        if(lqb(it).ne.0)then
	 do i=1,lqb(it)
	  xmtg(i,it)=xmtg(i,it)*(wmtg(it)-wmtg0(it))
     *    /(wmtg1(it)+wmtg2(it))
         enddo

         do nb=1,nbpom
	  if(ibs(nb).eq.it.and.(npompr(nb).ne.0.or.npomin(nb).ne.0))then
	   if(npomin(nb).ne.0)then
	    do np=1,npomin(nb)
	     xpomin(np,nb)=xpomin(np,nb)*(wmtg(it)-wmtg0(it))
     *       /(wmtg1(it)+wmtg2(it))
	    enddo
	   endif
	   if(npompr(nb).ne.0)then
	    do np=1,npompr(nb)
	     ipp=ilpr(np,nb)
	     lnp=lnpr(np,nb)
	     xpompr(lnp,ipp)=xpompr(lnp,ipp)*(wmtg(it)-wmtg0(it))
     *       /(wmtg1(it)+wmtg2(it))
	    enddo
	   endif
	  endif
	 enddo
	 
        elseif(wmtg2(it).gt.wmtg(it))then
         wmtg0(it)=wmtg(it)/2.d0
         do nb=1,nbpom
	  if(ibs(nb).eq.it.and.(npompr(nb).ne.0.or.npomin(nb).ne.0))then
	   if(npomin(nb).ne.0)then
	    do np=1,npomin(nb)
	     xpomin(np,nb)=xpomin(np,nb)*wmtg(it)/wmtg2(it)/2.d0
	    enddo
	   endif
	   if(npompr(nb).ne.0)then
	    do np=1,npompr(nb)
	     ipp=ilpr(np,nb)
	     lnp=lnpr(np,nb)
	     xpompr(lnp,ipp)=xpompr(lnp,ipp)*wmtg(it)/wmtg2(it)/2.d0
	    enddo
	   endif
	  endif
	 enddo
	else
         wmtg0(it)=wmtg(it)-wmtg2(it)
	endif
       else
        wmtg0(it)=wm0
       endif
      enddo	 
	 
c-------------------------------------------------
c treatment of low mass diffraction
      if(debug.ge.1)write (moniou,207)
      do ip=1,ia(1)                        !loop over proj. nucleons
       if(iwp(ip).eq.2)then                !diffraction dissociation
        it=iprcn(ip)
        if(debug.ge.2)write (moniou,208)ip,it
	if(iwt(it).eq.2)then
         call qgdifr(wppr0(ip),wmtg0(it),izp(ip),izt(it),-2,-2,iret)
	elseif(iwt(it).eq.-1)then
         call qgdifr(wppr0(ip),wmtg0(it),izp(ip),izt(it),-2,0,iret)
	elseif(iwt(it).eq.1)then
         call qgdifr(wppr0(ip),wmtg0(it),izp(ip),izt(it),-2,-1,iret)
	else
	 stop'wrong connection for diffraction'
	endif
	if(iret.eq.1)goto 1
       endif
      enddo
      
      do it=1,ia(2)                        !loop over targ. nucleons
       if(iwt(it).eq.2)then                !diffraction dissociation
        ip=itgcn(it)
        if(debug.ge.2)write (moniou,209)it,ip
	if(iwp(ip).eq.-1)then
         call qgdifr(wppr0(ip),wmtg0(it),izp(ip),izt(it),0,-2,iret)
	elseif(iwp(ip).eq.1)then
         call qgdifr(wppr0(ip),wmtg0(it),izp(ip),izt(it),-1,-2,iret)
	endif
	if(iret.eq.1)goto 1
       endif
      enddo

c-------------------------------------------------
c particle production for all cut Pomerons
      if(nbpom.ne.0)then
       if(debug.ge.1)write (moniou,210)
       do npb=1,nbpom                            !loop over collisions
        ip=ias(npb)                              !proj. index
        it=ibs(npb)                              !targ. index
        icdp=iddp(ip)                            !proj. diffr. eigenstate
        icdt=iddt(it)                            !targ. diffr. eigenstate
        bbp=bbpom(npb)                           !b^2 between proj. and targ.
	vvx=vvxpom(npb)                          !nuclear screening factor
        if(debug.ge.1)write (moniou,211)npb,ip,it,bbp,vvx,nqs(npb)
     *  ,npomin(npb),npompr(npb),npomtg(npb)
	
	if(npomin(npb).ne.0)then
	 do n=1,npomin(npb)                      !loop over interm. Pomerons
	  wpi=xpopin(n,npb)*wp0                  !LC+ for the Pomeron
	  wmi=xpomin(n,npb)*wm0                  !LC- for the Pomeron
          if(debug.ge.2)write (moniou,212)n,wpi,wmi
          if(wpi*wmi.ne.0.d0)then
           ic11=0
           ic12=0
           ic21=0
           ic22=0
           call qgstr(wpi,wmi,wppr0(ip),wmtg0(it)
     *     ,ic11,ic12,ic22,ic21,0,0)             !string hadronization
          endif	
	 enddo
	endif

	if(nqs(npb).ne.0)then
         do n=1,nqs(npb)                         !loop over single Pomerons
	  lnp=nnpr(n,npb)                        !index for proj. constituent
	  lnt=nntg(n,npb)                        !index for targ. constituent
	  lqa0(ip)=lqa0(ip)-1
	  lqb0(it)=lqb0(it)-1
	  xpi=xppr(lnp,ip)
	  xmi=xmtg(lnt,it)
          wpi=wp0*xpi                            !LC+ for the Pomeron
          wmi=wm0*xmi                            !LC- for the Pomeron
          sy=wpi*wmi
	  wsoft=qgpomc(sy,xpi,xmi,bbp,vvx,icdp,icdt,icz,0)!weight of soft inter.
	  wgg=qgpomc(sy,xpi,xmi,bbp,vvx,icdp,icdt,icz,1) !weight of gg-hard int.
	  wqg=qgpomc(sy,xpi,xmi,bbp,vvx,icdp,icdt,icz,2) !weight of qg-hard int.
	  wgq=qgpomc(sy,xpi,xmi,bbp,vvx,icdp,icdt,icz,3) !weight of gq-hard int.
	  wqq=qgpomc(sy,xpi,xmi,bbp,vvx,icdp,icdt,icz,4) !weight of qq-hard int.
	  aks=psran(b10)*(wsoft+wgg+wqg+wgq+wqq)
	  if(debug.ge.2)write (moniou,213)n,wpi,wmi

	  if(aks.lt.wsoft)then            !string hadronization for soft Pomeron 
           if(lqa0(ip).eq.0.and.lva(ip).eq.0)then
            call qgixxd(izp(ip),ic11,ic12,icz)
           else
            ic11=0
            ic12=0
           endif
           if(lqb0(it).eq.0.and.lvb(it).eq.0)then
            call qgixxd(izt(it),ic21,ic22,2)
           else
            ic21=0
            ic22=0
           endif
           call qgstr(wpi,wmi,wppr0(ip),wmtg0(it),ic11,ic12,ic22,ic21
     *     ,1,1)	  
	  else            !QCD evolution and hadronization for semi-hard Pomeron
	   if(lva(ip).eq.0.and.lvb(it).eq.0.and.aks.lt.wsoft+wqq)then
            iqq=3
	    lva(ip)=1
	    lvb(it)=1
	   elseif(lva(ip).eq.0.and.aks.lt.wsoft+wqq+wqg)then
            iqq=1
	    lva(ip)=1
	   elseif(lvb(it).eq.0.and.aks.lt.wsoft+wqq+wqg+wgq)then
            iqq=2
	    lvb(it)=1
	   else
            iqq=0
	   endif	   

           call qghot(wpi,wmi,dsqrt(bbp),vvx,nva,nvb,izp(ip),izt(it)
     *     ,icdp,icdt,icz,iqq,0)            !QCD evolution + jet hadronization
           if(iqq.eq.1.or.iqq.eq.3)ila(ip)=nva
           if(iqq.eq.2.or.iqq.eq.3)ilb(it)=nvb
	  endif
	 enddo
	endif
	 
        if(npompr(npb).ne.0)then
	 do l=1,npompr(npb)                !loop over proj. leg Pomerons
	  ipp=ilpr(l,npb)                  !proj. index
	  lnp=lnpr(l,npb)                  !index for proj. constituent
	  bbpr=bpompr(lnp,ipp)             !b^2 for the Pomeron
	  vvxp=vvxpr(lnp,ipp)              !screening factor
	  lqa0(ipp)=lqa0(ipp)-1
	  xpi=xppr(lnp,ipp)
	  xmi=xpompr(lnp,ipp)
          wpi=wp0*xpi                      !LC+ for the Pomeron
	  wmi=wm0*xmi                      !LC- for the Pomeron
	  sy=wpi*wmi
	  if(sy.ne.0.d0)then
	   wsoft=qglegc(sy,xpi,bbpr,vvxp,iddp(ipp),icz,0)!weight of soft inter.
	   wgg=qglegc(sy,xpi,bbpr,vvxp,iddp(ipp),icz,1)  !weight of gg-hard int.
	   wqg=qglegc(sy,xpi,bbpr,vvxp,iddp(ipp),icz,2)  !weight of qg-hard int.
	  else
	   wsoft=1.d0
	   wgg=0.d0
	   wqg=0.d0
	  endif
          aks=psran(b10)*(wsoft+wgg+wqg)
	  if(debug.ge.2)write (moniou,214)l,wpi,wmi

	  if(aks.le.wsoft)then            !string hadronization for soft Pomeron 
           if(lqa0(ipp).eq.0.and.lva(ipp).eq.0.and.sy.ne.0.d0)then
            call qgixxd(izp(ipp),ic11,ic12,icz)
           else
            ic11=0
            ic12=0
           endif
           ic21=0
           ic22=0
           call qgstr(wpi,wmi,wppr0(ipp),wmtg0(it),ic11,ic12,ic22,ic21
     *     ,1,0)	  

	  else            !QCD evolution and hadronization for semi-hard Pomeron
	   if(lva(ipp).eq.0.and.aks.lt.wsoft+wqg)then
            iqq=1
	    lva(ipp)=1
	   else
            iqq=0
	   endif	   
           call qghot(wpi,wmi,dsqrt(bbpr),vvxp,nva,nvb,izp(ipp),izt(it)
     *     ,iddp(ipp),icdt,icz,iqq,1)         !QCD evolution + jet hadronization
           if(iqq.eq.1)ila(ipp)=nva
	  endif	    
          call qglead(wppr0(ipp),wmtg0(it),lqa(ipp),lqb(it)
     *    ,lqa0(ipp)+ncola0(ipp),lqb0(it)+ncolb0(it),lva(ipp),lvb(it)
     *    ,izp(ipp),izt(it),ila(ipp),ilb(it),iret)            !remnant treatment
          if(iret.ne.0)goto 1
	 enddo
	endif	
	  
	if(npomtg(npb).ne.0)then
	 do l=1,npomtg(npb)                !loop over targ. leg Pomerons
	  itt=iltg(l,npb)                  !targ. index
	  lnt=lntg(l,npb)                  !index for targ. constituent
	  bbtg=bpomtg(lnt,itt)             !b^2 for the Pomeron
	  vvxt=vvxtg(lnt,itt)              !screening factor
	  lqb0(itt)=lqb0(itt)-1
	  xmi=xmtg(lnt,itt)
          wmi=wm0*xmi                      !LC- for the Pomeron
	  wpi=xpomtg(lnt,itt)*wp0          !LC+ for the Pomeron
	  sy=wpi*wmi
	  if(sy.ne.0.d0)then
	   wsoft=qglegc(sy,xmi,bbtg,vvxt,iddt(itt),2,0)!weight of soft inter.
	   wgg=qglegc(sy,xmi,bbtg,vvxt,iddt(itt),2,1)  !weight of gg-hard int.
	   wqg=qglegc(sy,xmi,bbtg,vvxt,iddt(itt),2,2)  !weight of qg-hard int.
	  else
	   wsoft=1.d0
	   wgg=0.d0
	   wqg=0.d0
	  endif
          aks=psran(b10)*(wsoft+wgg+wqg)
	  if(debug.ge.2)write (moniou,215)l,wpi,wmi

	  if(aks.le.wsoft)then            !string hadronization for soft Pomeron 
           ic11=0
           ic12=0
           if(lqb0(itt).eq.0.and.lvb(itt).eq.0.and.sy.ne.0.d0)then
            call qgixxd(izt(itt),ic21,ic22,2)
           else
            ic21=0
            ic22=0
           endif
           call qgstr(wpi,wmi,wppr0(ip),wmtg0(itt),ic11,ic12,ic22,ic21
     *     ,0,1)	  
	   	  
	  else            !QCD evolution and hadronization for semi-hard Pomeron
	   if(lvb(itt).eq.0.and.aks.lt.wsoft+wqg)then
            iqq=2
	    lvb(itt)=1
	   else
            iqq=0
	   endif	   
           call qghot(wpi,wmi,dsqrt(bbtg),vvxt,nva,nvb,izp(ip),izt(itt)
     *     ,icdp,iddt(itt),icz,iqq,2)         !QCD evolution + jet hadronization
           if(iqq.eq.2)ilb(itt)=nvb
	  endif
	  call qglead(wppr0(ip),wmtg0(itt),lqa(ip),lqb(itt)
     *    ,lqa0(ip)+ncola0(ip),lqb0(itt)+ncolb0(itt),lva(ip),lvb(itt)
     *    ,izp(ip),izt(itt),ila(ip),ilb(itt),iret)    !remnant treatment
          if(iret.ne.0)goto 1
	 enddo
	endif	  
        ncola0(ip)=ncola0(ip)-1
        ncolb0(it)=ncolb0(it)-1
        call qglead(wppr0(ip),wmtg0(it),lqa(ip),lqb(it)
     *  ,lqa0(ip)+ncola0(ip),lqb0(it)+ncolb0(it),lva(ip),lvb(it)
     *  ,izp(ip),izt(it),ila(ip),ilb(it),iret)        !remnant treatment
        if(iret.ne.0)goto 1
       enddo     	                              !end of collision loop
      endif
	 
      if(nj.ne.0)then                   !arrangement of parton color connections
       if(debug.ge.1)write (moniou,216)nj
       call qgjarr(jfl)
       if(jfl.eq.0)goto 1
       if(debug.ge.1)write (moniou,217)
       call qgxjet                      !jet hadronization
      endif      
      if(debug.ge.1)write (moniou,218)

201   format(2x,'qgsha - inelastic interaction, N of Pomeron blocks:'    !so161205
     *,i4)                                                               !so161205
202   format(2x,'qgsha: energy-momentum sharing between Pomerons')
203   format(2x,'qgsha: log of acceptance probability - ',e10.3
     */4x,'N of rejections - ',i4,2x,'N of renorm. - ',i3)
204   format(2x,'qgsha:  new normalization of acceptance,'
     *,' N of renorm. - ',i3)
205   format(2x,'qgsha: leading remnant LC momenta')
206   format(2x,'qgsha: momentum conservation '
     *,'(correction for 3p-vertexes)')
207   format(2x,'qgsha: treatment of low mass diffraction')
208   format(2x,'qgsha: diffraction of ',i3,'-th proj. nucleon,'
     *,' recoil of ',i3,'-th targ. nucleon')
209   format(2x,'qgsha: diffraction of ',i3,'-th targ. nucleon,'
     *,' recoil of ',i3,'-th proj. nucleon')
210   format(2x,'qgsha: particle production for all cut Pomerons')
211   format(2x,'qgsha: ',i4,'-th collision,  proj. index - ',i3,2x
     *,'targ. index - ',i3
     */4x,'b^2=',e10.3,2x,'vvx=',e10.3,2x,'N of single Pomerons - ',i3
     *,2x,' N of interm. Pomerons - ',i3
     */4x,'N of proj. legs - ',i3,2x,'N of targ. legs - ',i3)
212   format(2x,'qgsha: particle production for '
     *,i3,'-th interm. Pomeron'
     */4x,'light cone momenta for the Pomeron:',2e10.3)
213   format(2x,'qgsha: particle production for '
     *,i3,'-th single Pomeron'
     */4x,'light cone momenta for the Pomeron:',2e10.3)
214   format(2x,'qgsha: particle production for '
     *,i3,'-th proj. leg Pomeron'
     */4x,'light cone momenta for the Pomeron:',2e10.3)
215   format(2x,'qgsha: particle production for '
     *,i3,'-th targ. leg Pomeron'
     */4x,'light cone momenta for the Pomeron:',2e10.3)
216   format(2x,'qgsha: arrangement of color connections for '
     *,i5,' final partons')
217   format(2x,'qgsha: jet hadronization')
218   format(2x,'qgsha - end')
      return
      end

c=============================================================================
      subroutine qglead(wppr0,wmtg0,lqa,lqb,lqa0,lqb0,lva,lvb
     *,izp,izt,ila,ilb,iret)
c-------------------------------------------------------------------------
c qglead - treatment of leading hadron states
c-------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(njmax=50000)
      common /qgarr37/ eqj(4,njmax),iqj(njmax),ncj(2,njmax),nj
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)
      
      iret=0
      if(lqa0.eq.0.and.lqb0.eq.0)then
       if(lva.eq.0.and.lvb.eq.0)then
	call qgdifr(wppr0,wmtg0,izp,izt,lqa,lqb,iret)
       elseif(lva.eq.0)then
        call qgdifr(wppr0,wmtg0,izp,izt,lqa,-1,iret)
       elseif(lvb.eq.0)then
        call qgdifr(wppr0,wmtg0,izp,izt,-1,lqb,iret)
       endif
       if(lva.eq.1)then
        eqj(1,ila)=.5d0*wppr0
        eqj(2,ila)=eqj(1,ila)
        eqj(3,ila)=0.d0
        eqj(4,ila)=0.d0
       endif
       if(lvb.eq.1)then
        eqj(1,ilb)=.5d0*wmtg0
        eqj(2,ilb)=-eqj(1,ilb)
        eqj(3,ilb)=0.d0
        eqj(4,ilb)=0.d0
       endif
      elseif(lqa0.eq.0)then
       if(lva.eq.0)then
        call qgdifr(wppr0,wmtg0,izp,izt,lqa,-1,iret)
       else
        eqj(1,ila)=.5d0*wppr0
        eqj(2,ila)=eqj(1,ila)
        eqj(3,ila)=0.d0
        eqj(4,ila)=0.d0
       endif
      elseif(lqb0.eq.0)then
       if(lvb.eq.0)then
        call qgdifr(wppr0,wmtg0,izp,izt,-1,lqb,iret)
       else
        eqj(1,ilb)=.5d0*wmtg0
        eqj(2,ilb)=-eqj(1,ilb)
        eqj(3,ilb)=0.d0
        eqj(4,ilb)=0.d0
       endif
      endif
      if(debug.ge.4)write (moniou,202)
      
201   format(2x,'qglead - treatment of leading hadron states')
202   format(2x,'qglead - end: iret=',i2)
      return
      end
	
c=============================================================================
      subroutine qgprox(imode)
c-------------------------------------------------------------------------
c qgprox - propose Pomeron end LC momenta 
c imod = 0 - to define normalization
c imod = 1 - propose values according to x^(delf-alpd) * (1 - sum_i x_i)^ahl
c-------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207,npbmax=1000,npnmax=1000,npmax=5000
     *,legmax=900)
cdh   dimension alpi(npnmax,iapmax),alti(npnmax,iapmax)
cdh  *,alpt(iapmax),altt(iapmax)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr9/  iwp(iapmax),iwt(iapmax),lqa(iapmax),lqb(iapmax)
     *,iprcn(iapmax),itgcn(iapmax),ias(npbmax),ibs(npbmax),nqs(npbmax)
     *,npompr(npbmax),npomtg(npbmax),npomin(npbmax),nnpr(npmax,npbmax)
     *,nntg(npmax,npbmax),ilpr(legmax,npbmax),iltg(legmax,npbmax)
     *,lnpr(legmax,npbmax),lntg(legmax,npbmax)
     *,nbpi(npnmax,iapmax),nbti(npnmax,iapmax),idnpi(npnmax,iapmax)
     *,idnti(npnmax,iapmax),nppi(npnmax,iapmax),npti(npnmax,iapmax)
     *,nlpi(npnmax,iapmax),nlti(npnmax,iapmax)
      common /qgarr11/ b10
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr19/ ahl(3)
      common /qgarr22/ wppr(iapmax),wmtg(iapmax)
      common /qgarr23/ bbpom(npbmax),vvxpom(npbmax)
     *,bpompr(npnmax,iapmax),bpomtg(npnmax,iapmax),vvxpr(npnmax,iapmax)
     *,vvxtg(npnmax,iapmax),xpompr(npnmax,iapmax),xpomtg(npnmax,iapmax)
     *,xpopin(npmax,npbmax),xpomin(npmax,npbmax)
      common /qgarr40/ xppr(npnmax,iapmax),xmtg(npnmax,iapmax)
      common /qgarr43/ moniou
      common /debug/   debug
      external psran
       
      if(debug.ge.3)write (moniou,201)imode
      
      delf=dels
      if(imode.eq.0)then                    !0-configuration (for normalization)
       do ip=1,ia(1)                        !loop over proj. nucleons
        if(lqa(ip).ne.0)then
         do n=1,lqa(ip)                     !loop over proj. constituents
	  if(idnpi(n,ip).eq.0)then
	   xppr(n,ip)=1.d0/wp0              !LC+ for single Pomeron
	  else
	   xppr(n,ip)=1.d0/xpompr(n,ip)/scm !LC+ for leg Pomeron
	  endif
 	 enddo
        endif
       enddo		    
       do it=1,ia(2)                        !loop over targ. nucleons
        if(lqb(it).ne.0)then
         do n=1,lqb(it)                     !loop over targ. constituents
	  if(idnti(n,it).eq.0)then
	   xmtg(n,it)=1.d0/wm0              !LC- for single Pomeron
	  else
	   xmtg(n,it)=1.d0/xpomtg(n,it)/scm !LC- for leg Pomeron
	  endif
	 enddo
        endif
       enddo	

      else                                  !proposed configuration 
       do ip=1,ia(1)                        !loop over proj. nucleons
        if(lqa(ip).ne.0)then
	 xpt=wppr(ip)/wp0
         do n=1,lqa(ip)                     !loop over proj. constituents
	  nrej=0
	  alfl=ahl(icz)+(lqa(ip)-n)*(1.d0-alpd+delf)
	  if(icz.eq.2)alfl=alfl-float(lqa(ip)-1)/lqa(ip)  !baryon "junction"
	  gb0=(1.d0-.11d0**(1.d0/(1.d0-alpd+delf)))**alfl
     *    *exp(alfl*(1.d0-alpd+delf)*.11d0)*2.d0
1         continue
c proposal functions are chosen depending on the parameters
c to assure an efficient procedure
          if(delf-alpd.ge.0.d0.and.alfl.ge.0.d0
     *    .or.delf-alpd.lt.0.d0.and.alfl.le.0.d0)then
           up=1.d0-psran(b10)**(1.d0/(1.d0-alpd+delf))
	   if(1.d0-up.lt.1.d-20)goto 1
           tp=1.d0-up**(1.d0/(1.d0+alfl))
	   gb=(tp/(1.d0-up))**(-alpd+delf)
          elseif(delf-alpd.lt.0.d0.and.alfl.gt.0.d0)then
           up=-log(1.d0-psran(b10)*(1.d0-exp(-alfl*(1.d0-alpd+delf))))
     *     /alfl/(1.d0-alpd+delf)
           tp=up**(1.d0/(1.d0-alpd+delf))
	   gb=(1.d0-tp)**alfl*exp(alfl*(1.d0-alpd+delf)*up)/gb0
	  else
           tp=1.d0-psran(b10)**(1.d0/(1.d0+alfl))
	   gb=tp**(delf-alpd)
          endif
	  if(psran(b10).gt.gb)then
	   nrej=nrej+1
	   goto 1
          endif
	  xppr(n,ip)=tp*xpt                 !proposed LC+ for the constituent 
          xpt=xpt-xppr(n,ip)                !LC+ of the remnant
 	 enddo
        endif
       enddo	
	    
       do it=1,ia(2)                        !loop over targ. nucleons
        if(lqb(it).ne.0)then
	 xmt=wmtg(it)/wm0
         do n=1,lqb(it)                     !loop over targ. constituents
	  nrej=0
	  alfl=ahl(2)+(lqb(it)-n)*(1.d0-alpd+delf)
     *    -float(lqb(it)-1)/lqb(it)                       !baryon "junction"
	  gb0=(1.d0-.11d0**(1.d0/(1.d0-alpd+delf)))**alfl
     *    *exp(alfl*(1.d0-alpd+delf)*.11d0)*2.d0
2         continue
          if(delf-alpd.ge.0.d0.and.alfl.ge.0.d0
     *    .or.delf-alpd.lt.0.d0.and.alfl.le.0.d0)then
           up=1.d0-psran(b10)**(1.d0/(1.d0-alpd+delf))
	   if(1.d0-up.lt.1.d-20)goto 2
           tp=1.d0-up**(1.d0/(1.d0+alfl))
	   gb=(tp/(1.d0-up))**(-alpd+delf)
          elseif(delf-alpd.lt.0.d0.and.alfl.gt.0.d0)then
           up=-log(1.d0-psran(b10)*(1.d0-exp(-alfl*(1.d0-alpd+delf))))
     *     /alfl/(1.d0-alpd+delf)
           tp=up**(1.d0/(1.d0-alpd+delf))
	   gb=(1.d0-tp)**alfl*exp(alfl*(1.d0-alpd+delf)*up)/gb0
	  else
           tp=1.d0-psran(b10)**(1.d0/(1.d0+alfl))
	   gb=tp**(delf-alpd)
          endif
          if(psran(b10).gt.gb)then
	   nrej=nrej+1
	   goto 2
          endif
          xmtg(n,it)=tp*xmt                 !proposed LC- for the constituent
          xmt=xmt-xmtg(n,it)                !LC- of the remnant
 	 enddo
        endif
       enddo   	    
      endif
      if(debug.ge.4)write (moniou,202)

201   format(2x,'qgprox - propose Pomeron end LC momenta, imode=',i2)
202   format(2x,'qgprox - end')
      return        
      end
      
c=============================================================================
      double precision function qgweix(nbpom)
c-------------------------------------------------------------------------
c qgweix - log-weight of x-configuration 
c imod = 0 - to define normalization
c imod = 1 - propose values according to x^(delf-alpd) * (1 - sum_i x_i)^ahl
c-------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207,npbmax=1000,npnmax=1000,npmax=5000
     *,legmax=900)
cdh   dimension alpi(npnmax,iapmax),alti(npnmax,iapmax)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr9/  iwp(iapmax),iwt(iapmax),lqa(iapmax),lqb(iapmax)
     *,iprcn(iapmax),itgcn(iapmax),ias(npbmax),ibs(npbmax),nqs(npbmax)
     *,npompr(npbmax),npomtg(npbmax),npomin(npbmax),nnpr(npmax,npbmax)
     *,nntg(npmax,npbmax),ilpr(legmax,npbmax),iltg(legmax,npbmax)
     *,lnpr(legmax,npbmax),lntg(legmax,npbmax)
     *,nbpi(npnmax,iapmax),nbti(npnmax,iapmax),idnpi(npnmax,iapmax)
     *,idnti(npnmax,iapmax),nppi(npnmax,iapmax),npti(npnmax,iapmax)
     *,nlpi(npnmax,iapmax),nlti(npnmax,iapmax)
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr23/ bbpom(npbmax),vvxpom(npbmax)
     *,bpompr(npnmax,iapmax),bpomtg(npnmax,iapmax),vvxpr(npnmax,iapmax)
     *,vvxtg(npnmax,iapmax),xpompr(npnmax,iapmax),xpomtg(npnmax,iapmax)
     *,xpopin(npmax,npbmax),xpomin(npmax,npbmax)
      common /qgarr40/ xppr(npnmax,iapmax),xmtg(npnmax,iapmax)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)nbpom

      delf=dels
      qgweix=0.d0
      do npb=1,nbpom                              !loop over collisions
       ip=ias(npb)                                !proj. index
       it=ibs(npb)                                !targ. index
       icdp=iddp(ip)                              !proj. diffr. eigenstate
       icdt=iddt(it)                              !targ. diffr. eigenstate
       bbp=bbpom(npb)                             !b^2 between proj. and targ. 
       vvx=vvxpom(npb)                            !nuclear screening factor
       if(nqs(npb).ne.0)then
        do n=1,nqs(npb)                           !loop over single Pomerons
	 lnp=nnpr(n,npb)                          !proj. constituent index
	 lnt=nntg(n,npb)                          !targ. constituent index
	 xpp=xppr(lnp,ip)                         !LC+ for the Pomeron
	 xpm=xmtg(lnt,it)                         !LC- for the Pomeron
         qgweix=qgweix+dlog(qgpomc(scm*xpp*xpm,xpp,xpm,bbp,vvx
     *   ,icdp,icdt,icz,-1)/(xpp*xpm)**(delf-alpd)) !add single Pomeron contrib.
        enddo
       endif
       if(npompr(npb).ne.0)then
	do l=1,npompr(npb)                        !loop over proj. leg Pomerons
	 ipp=ilpr(l,npb)                          !proj. index
	 lnp=lnpr(l,npb)                          !proj. constituent index
	 xpp=xppr(lnp,ipp)                        !LC+ for the Pomeron
	 xpomr=1.d0/xpompr(lnp,ipp)/scm           !LC+ for the 3P vertex
	 vvxp=vvxpr(lnp,ipp)                      !screening factor
	 bbpr=bpompr(lnp,ipp)	                  !b^2 for the Pomeron
         qgweix=qgweix+dlog(qglegc(xpp/xpomr,xpp,bbpr,vvxp
     *   ,iddp(ipp),icz,-1)/xpp**(-alpd+delf))    !add leg Pomeron contrib.
	enddo
       endif
       if(npomtg(npb).ne.0)then
	do l=1,npomtg(npb)                        !loop over targ. leg Pomerons
	 itt=iltg(l,npb)                          !targ. index
	 lnt=lntg(l,npb)                          !targ. constituent index
	 xpm=xmtg(lnt,itt)                        !LC- for the Pomeron
	 xpomr=xpomtg(lnt,itt)                    !LC+ for the 3P vertex
	 vvxt=vvxtg(lnt,itt)                      !screening factor
	 bbtg=bpomtg(lnt,itt)	                  !b^2 for the Pomeron
         qgweix=qgweix+dlog(qglegc(xpomr*scm*xpm,xpm,bbtg,vvxt
     *   ,iddt(itt),2,-1)/xpm**(-alpd+delf))      !add leg Pomeron contrib.
	enddo
       endif
      enddo
      if(debug.ge.4)write (moniou,202)qgweix

201   format(2x,'qgweix - log-weight of x-configuration,'
     *,' N of collisions - ',i4)
202   format(2x,'qgweix=',e10.3)
      return        
      end

c=============================================================================
      subroutine qghot(wpp,wpm,b,vvx,nva,nvb,izp,izt,icdp,icdt,icz,iqq
     *,jpt)
c---------------------------------------------------------------------------
c qghot - semi-hard process
c wpp,wpm   - LC momenta for the constituent partons,
c b         - impact parameter for the semi-hard Pomeron,
c izp, izt  - types of proj. and targ. remnants,
c icdp,icdt - proj. and targ.  diffractive eigenstates,
c iqq - type of the semi-hard process: 0 - gg, 1 - q_vg, 2 - gq_v, 3 - q_vq_v
c jpt=0 - single Pomeron,
c jpt=1 - proj. leg Pomeron,
c jpt=2 - targ. leg Pomeron
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      character*2 tyq
      parameter(njmax=50000)
      dimension ept(4),ep3(4),ey(3),ebal(4),
     *qmin(2),wp(2),iqc(2),iqp(2),nqc(2),ncc(2,2),
     *qv1(30,50),zv1(30,50),qm1(30,50),iqv1(30,50),
     *ldau1(30,49),lpar1(30,50),
     *qv2(30,50),zv2(30,50),qm2(30,50),iqv2(30,50),
     *ldau2(30,49),lpar2(30,50)
      common /qgarr2/  scm,wp0,wm0
      common /qgarr6/  pi,bm,amws
      common /qgarr8/  wwm,be(4),dc(5),deta,almpt,ptdif,ptndi
      common /qgarr10/ am(7),ammu
      common /qgarr11/ b10
      common /qgarr12/ nsp
      common /qgarr15/ fp(3),rq(3),cd(2,3)
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr26/ factk,fqscal
      common /qgarr37/ eqj(4,njmax),iqj(njmax),ncj(2,njmax),nj
      common /qgarr42/ tyq(16)
      common /qgarr43/ moniou
      common /debug/   debug
      external psran

      if(debug.ge.1)write (moniou,201)iqq,wpp,wpm,izp,izt,icdp,icdt
     *,icz,jpt,nj
     
      nj0=nj                       !store number of final partons
      nsp0=nsp                     !store number of final particles

1     sy=wpp*wpm  !energy squared for semi-hard inter. (including preevolution)
      nj=nj0
      nsp=nsp0

      if(iqq.eq.3)then             !q_vq_v-ladder
       wpi=wpp                     !LC+ for the hard interaction
       wmi=wpm                     !LC- for the hard interaction
      else

c-------------------------------------------------
c normalization of acceptance
       s2min=4.d0*fqscal*qt0       !threshold energy
       if(sy.lt.s2min)stop'qghot: sy<s2min!!!'
       xmin=s2min/sy
       iq=(iqq+1)/2+1              !auxilliary type of parton (1 - g, 2 - q(q~))
       sj=qgjit(qt0,qt0,sy,1,iq)   !inclusive parton-parton cross-sections
       if(iqq.eq.0)then
        gb0=-dlog(xmin)*(1.d0-dsqrt(xmin))**(2.d0*betp)*sj
       else
        gb0=(1.d0-xmin)**betp*sj
       endif
       if(jpt.eq.0)then            !single Pomeron
        if(iqq.eq.0)then
         gb0=gb0/(rq(icz)+rq(2)+alfp*dlog(scm/sy))
     *   *exp(-b*b/(4.d0*.0389d0*(rq(icz)+rq(2)+alfp*dlog(scm/s2min))))
        elseif(iqq.eq.1)then
         gb0=gb0/(rq(icz)+rq(2)+alfp*dlog(wm0/wpm))*exp(-b*b
     *   /(4.d0*.0389d0*(rq(icz)+rq(2)+alfp*dlog(wpp*wm0/s2min))))
        elseif(iqq.eq.2)then
         gb0=gb0/(rq(icz)+rq(2)+alfp*dlog(wp0/wpp))*exp(-b*b
     *   /(4.d0*.0389d0*(rq(icz)+rq(2)+alfp*dlog(wpm*wp0/s2min))))
        endif
       elseif(jpt.eq.1)then        !proj. leg Pomeron
        if(iqq.eq.0)then
         gb0=gb0/(rq(icz)+alfp*dlog(wp0/wpp))
     *   *exp(-b*b/(4.d0*.0389d0*(rq(icz)+alfp*dlog(wp0*wpm/s2min))))
        elseif(iqq.eq.1)then
         gb0=gb0/rq(icz)
     *   *exp(-b*b/(4.d0*.0389d0*(rq(icz)+alfp*dlog(sy/s2min))))
        endif
       elseif(jpt.eq.2)then        !targ. leg Pomeron
        if(iqq.eq.0)then
         gb0=gb0/(rq(2)+alfp*dlog(wm0/wpm))
     *   *exp(-b*b/(4.d0*.0389d0*(rq(2)+alfp*dlog(wm0*wpp/s2min))))
        elseif(iqq.eq.2)then
         gb0=gb0/rq(2)
     *   *exp(-b*b/(4.d0*.0389d0*(rq(2)+alfp*dlog(sy/s2min))))
        endif
       endif
	    
c-------------------------------------------------
c sharing of LC momenta between soft preevolution and hard ladder
2      zpm=(1.d0-psran(b10)*(1.d0-xmin**(delh-dels)))
     * **(1.d0/(delh-dels))
       sjqq=qgjit(qt0,qt0,zpm*sy,2,2)  !inclusive qq cross-section
       sjqg=qgjit(qt0,qt0,zpm*sy,1,2)  !inclusive qg cross-section
       sjgg=qgjit(qt0,qt0,zpm*sy,1,1)  !inclusive gg cross-section

       if(iqq.eq.0)then              !gg-ladder
        xp=zpm**psran(b10)           !LC+ momentum share
        xm=zpm/xp                    !LC- momentum share
        wpi=wpp*xp                   !LC+ for the hard interaction
        wmi=wpm*xm                   !LC- for the hard interaction
	if(jpt.eq.0)then             !single Pomeron
	 v1pnu0=qgfani(wp0/wpp/xp,b*b/4.d0,vvx,0.d0,0.d0,icdp,icz,1)
         v1tnu0=qgfani(wpp/wp0*xp*scm,b*b/4.d0,vvx,0.d0,0.d0,icdt,2,1)
         nn=0
21       nn=nn+1
         vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx)
         vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx)
         v1tnu=qgfani(wpp/wp0*xp*scm,b*b/4.d0,vvxt,0.d0,0.d0,icdt,2,1)
         v1pnu=qgfani(wp0/wpp/xp,b*b/4.d0,vvxp,0.d0,0.d0,icdp,icz,1)
         if((abs(v1pnu0-v1pnu).gt.1.d-1.or.abs(v1tnu0-v1tnu).gt.1.d-1)
     *   .and.nn.lt.100)then
          v1pnu0=v1pnu
          v1tnu0=v1tnu
	  goto 21
         endif
         vvx1=1.d0-exp(-2.d0*(v1pnu+v1tnu))*(1.d0-vvx)**2
       
         v1pnu0=qgfani(wpm/wm0*xm*scm,b*b/4.d0,vvx,0.d0,0.d0,icdp,icz,1)
         v1tnu0=qgfani(wm0/wpm/xm,b*b/4.d0,vvx,0.d0,0.d0,icdt,2,1)
         nn=0
22       nn=nn+1
         vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx)
         vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx)
         v1tnu=qgfani(wm0/wpm/xm,b*b/4.d0,vvxt,0.d0,0.d0,icdt,2,1)
         v1pnu=qgfani(wpm/wm0*xm*scm,b*b/4.d0,vvxp,0.d0,0.d0,icdp,icz,1)
         if((abs(v1pnu0-v1pnu).gt.1.d-1.or.abs(v1tnu0-v1tnu).gt.1.d-1)
     *   .and.nn.lt.100)then
          v1pnu0=v1pnu
          v1tnu0=v1tnu
	  goto 22
         endif
         vvx2=1.d0-exp(-2.d0*(v1pnu+v1tnu))*(1.d0-vvx)**2

	 glu1=qgppdi(xp,vvx1,0)      !gluon PDF in the upper Pomeron
	 sea1=qgppdi(xp,vvx1,1)      !quark PDF in the upper Pomeron
	 glu2=qgppdi(xm,vvx2,0)      !gluon PDF in the lower Pomeron
	 sea2=qgppdi(xm,vvx2,1)      !quark PDF in the lower Pomeron
	 wwgg=glu1*glu2*sjgg
	 wwqg=sea1*glu2*sjqg
	 wwgq=glu1*sea2*sjqg
	 wwqq=sea1*sea2*sjqq
	else                         !leg Pomeron
	 glu1=qgppdi(xp,vvx,0)
	 sea1=qgppdi(xp,vvx,1)
	 glu2=qgppdi(xm,vvx,0)
	 sea2=qgppdi(xm,vvx,1)
	 wwgg=glu1*glu2*sjgg
	 wwqg=sea1*glu2*sjqg
	 wwgq=glu1*sea2*sjqg
	 wwqq=sea1*sea2*sjqq
	endif
        gbyj=-dlog(zpm)*(wwgg+wwqg+wwgq+wwqq)
	if(jpt.eq.0)then
	 rh=rq(icz)+rq(2)-alfp*dlog(wpp/wp0*wpm/wm0*zpm)
	elseif(jpt.eq.1)then
	 rh=rq(icz)-alfp*dlog(wpp/wp0*zpm)
	elseif(jpt.eq.2)then
	 rh=rq(2)-alfp*dlog(wpm/wm0*zpm)
	endif
	gbyj=gbyj/rh*exp(-b*b/(4.d0*.0389d0*rh))

       else                          !q_vg-(gq_v-)ladder
        if(iqq.eq.1)then             !q_vg-ladder
         wpi=wpp
         wmi=wpm*zpm
	 xm=zpm
	 if(jpt.eq.0)then            !single Pomeron
	  v1pnu0=qgfani(wpm/wm0*xm*scm,b*b/4.d0,vvx,0.d0,0.d0
     *    ,icdp,icz,1)
          v1tnu0=qgfani(wm0/wpm/xm,b*b/4.d0,vvx,0.d0,0.d0,icdt,2,1)
          nn=0
23        nn=nn+1
          vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx)
          vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx)
	  v1tnu=qgfani(wm0/wpm/xm,b*b/4.d0,vvxt,0.d0,0.d0,icdt,2,1)
          v1pnu=qgfani(wpm/wm0*xm*scm,b*b/4.d0,vvxp,0.d0,0.d0
     *    ,icdp,icz,1)
          if((abs(v1pnu0-v1pnu).gt.1.d-1.or.abs(v1tnu0-v1tnu).gt.1.d-1)
     *    .and.nn.lt.100)then
           v1pnu0=v1pnu
           v1tnu0=v1tnu
	   goto 23
          endif
          vvx2=1.d0-exp(-2.d0*(v1pnu+v1tnu))*(1.d0-vvx)**2

	  wwqg=qgppdi(zpm,vvx2,0)*sjqg
	  wwqq=qgppdi(zpm,vvx2,1)*sjqq
	 else                        !leg Pomeron
	  wwqg=qgppdi(zpm,vvx,0)*sjqg
	  wwqq=qgppdi(zpm,vvx,1)*sjqq
	 endif
        elseif(iqq.eq.2)then         !gq_v-ladder
         wpi=wpp*zpm
         wmi=wpm
	 xp=zpm
	 if(jpt.eq.0)then            !single Pomeron
	  v1pnu0=qgfani(wp0/wpp/xp,b*b/4.d0,vvx,0.d0,0.d0,icdp,icz,1)
          v1tnu0=qgfani(wpp/wp0*xp*scm,b*b/4.d0,vvx,0.d0,0.d0,icdt,2,1)
          nn=0
24        nn=nn+1
          vvxt=1.d0-exp(-v1pnu0)*(1.d0-vvx)
          vvxp=1.d0-exp(-v1tnu0)*(1.d0-vvx)
	  v1tnu=qgfani(wpp/wp0*xp*scm,b*b/4.d0,vvxt,0.d0,0.d0,icdt,2,1)
          v1pnu=qgfani(wp0/wpp/xp,b*b/4.d0,vvxp,0.d0,0.d0,icdp,icz,1)
	  if((abs(v1pnu0-v1pnu).gt.1.d-1.or.abs(v1tnu0-v1tnu).gt.1.d-1)
     *    .and.nn.lt.100)then
           v1pnu0=v1pnu
           v1tnu0=v1tnu
	   goto 24
          endif
	  vvx1=1.d0-exp(-2.d0*(v1pnu+v1tnu))*(1.d0-vvx)**2
       
	  wwqg=qgppdi(zpm,vvx1,0)*sjqg
	  wwqq=qgppdi(zpm,vvx1,1)*sjqq
	 else                        !leg Pomeron
	  wwqg=qgppdi(zpm,vvx,0)*sjqg
	  wwqq=qgppdi(zpm,vvx,1)*sjqq
	 endif
        endif
        gbyj=wwqg+wwqq
	if(jpt.eq.0)then
         if(iqq.eq.1)then
	  rh=rq(icz)+rq(2)-alfp*dlog(wpm/wm0*zpm)
	 else
	  rh=rq(icz)+rq(2)-alfp*dlog(wpp/wp0*zpm)
	 endif
	elseif(jpt.eq.1)then
	 rh=rq(icz)-alfp*dlog(zpm)
	elseif(jpt.eq.2)then
	 rh=rq(2)-alfp*dlog(zpm)
	endif
	gbyj=gbyj/rh*exp(-b*b/(4.d0*.0389d0*rh))
       endif

       gbyj=gbyj/gb0/zpm**delh
       if(psran(b10).gt.gbyj)goto 2
      endif
      if(debug.ge.2)write (moniou,202)wpi*wmi

11    wpi1=wpi
      wmi1=wmi
      wpq=0.d0
      wmq=0.d0
      nj=nj0                     !initialization for the number of final partons
      rrr=psran(b10)
      jqq=0                                  !gg-ladder
      if(iqq.eq.1.or.iqq.eq.2)then
       if(rrr.lt.wwqq/(wwqg+wwqq))jqq=1      !q_vq_s-laddder
      elseif(iqq.eq.0)then
       if(rrr.lt.wwqg/(wwgg+wwqg+wwgq+wwqq))then
        jqq=1                                !q_sg-ladder
       elseif(rrr.lt.(wwqg+wwgq)/(wwgg+wwqg+wwgq+wwqq))then
        jqq=2                                !gq_s-ladder
       elseif(rrr.lt.(wwqg+wwgq+wwqq)/(wwgg+wwqg+wwgq+wwqq))then
        jqq=3                                !q_sq_s-ladder
       endif
      endif

c-------------------------------------------------
c parton types for the ladder legs and for the leading jets
c iqc(1) - flavor for the upper quark (0 in case of gluon),
c iqc(2) - the same for the lower one
      if(iqq.ne.0.and.iqq.ne.2)then          !q_v from the proj.
       call qgvdef(izp,ic1,ic2,icz)          !leading state flavor
       iqc(1)=ic1                            !upper leg parton
       nj=nj+1
       if(nj.gt.njmax)stop'increase njmax!!!'
       nva=nj
       iqj(nj)=ic2                           !leading jet parton
       ncc(1,1)=nj                           !color connection with leading jet
       ncc(2,1)=0
      else                                   !g(q_s) from the proj.
       nj=nj+1
       if(nj.gt.njmax)stop'increase njmax!!!'
       if(psran(b10).lt.dc(2))then
        iqj(nj)=-4
       else
        iqj(nj)=-int(2.d0*psran(b10)+1.d0)
       endif
       iqj(nj+1)=-iqj(nj)
       wp1=wpp-wpi
       wp2=wp1*psran(b10)
       wp1=wp1-wp2
       eqj(1,nj)=.5d0*wp1
       eqj(2,nj)=eqj(1,nj)
       eqj(3,nj)=0.d0
       eqj(4,nj)=0.d0
       eqj(1,nj+1)=.5d0*wp2
       eqj(2,nj+1)=eqj(1,nj+1)
       eqj(3,nj+1)=0.d0
       eqj(4,nj+1)=0.d0
       if(jqq.eq.0.or.iqq.eq.0.and.jqq.eq.2)then
        iqc(1)=0
        ncc(1,1)=nj
        ncc(2,1)=nj+1
        nj=nj+1
        if(nj.gt.njmax)stop'increase njmax!!!'
       else
        if(psran(b10).lt.dc(2))then
         iqc(1)=3*(2.d0*int(.5d0+psran(b10))-1.d0)
	else
         iqc(1)=int(2.d0*psran(b10)+1.d0)
     *   *(2.d0*int(.5d0+psran(b10))-1.d0)
        endif
12      zg=xp+psran(b10)*(1.d0-xp)           !gluon splitting into qq~
        if(psran(b10).gt.zg**dels*((1.d0-xp/zg)/ (1.d0-xp))**betp)
     *  goto 12
        xg=xp/zg
        wpq0=wpp*(xg-xp)
        wmq=1.d0/wpq0
        wmi1=wmi1-wmq
        if(wmi1*wpi1.le.s2min)goto 11
        nj=nj+2
        if(nj.gt.njmax)stop'increase njmax!!!'
	iqj(nj)=-iqc(1)
	if(iabs(iqc(1)).eq.3)iqj(nj)=iqj(nj)*4/3
        eqj(1,nj)=.5d0*wmq
        eqj(2,nj)=-.5d0*wmq
        eqj(3,nj)=0.d0
        eqj(4,nj)=0.d0
        if(iqc(1).gt.0)then
         ncj(1,nj)=nj-1
         ncj(1,nj-1)=nj
         ncj(2,nj)=0
         ncj(2,nj-1)=0
         ncc(1,1)=nj-2
         ncc(2,1)=0
        else
         ncj(1,nj)=nj-2
         ncj(1,nj-2)=nj
         ncj(2,nj)=0
         ncj(2,nj-2)=0
         ncc(1,1)=nj-1
         ncc(2,1)=0
        endif
       endif
      endif

      if((iqq-2)*(iqq-3)*(iqq-4).eq.0)then     !q_v from the targ.
       call qgvdef(izt,ic1,ic2,2)              !leading state flavor
       iqc(2)=ic1                              !lower leg parton
       nj=nj+1
       if(nj.gt.njmax)stop'increase njmax!!!'
       nvb=nj
       iqj(nj)=ic2
       ncc(1,2)=nj
       ncc(2,2)=0
      else
       nj=nj+1
       if(nj.gt.njmax)stop'increase njmax!!!'
       if(psran(b10).lt.dc(2))then
        iqj(nj)=-4
       else
        iqj(nj)=-int(2.d0*psran(b10)+1.d0)
       endif
       iqj(nj+1)=-iqj(nj)
       wm1=wpm-wmi
       wm2=wm1*psran(b10)
       wm1=wm1-wm2
       eqj(1,nj)=.5d0*wm1
       eqj(2,nj)=-eqj(1,nj)
       eqj(3,nj)=0.d0
       eqj(4,nj)=0.d0
       eqj(1,nj+1)=.5d0*wm2
       eqj(2,nj+1)=-eqj(1,nj+1)
       eqj(3,nj+1)=0.d0
       eqj(4,nj+1)=0.d0
       if(jqq.eq.0.or.iqq.eq.0.and.jqq.eq.1)then
        iqc(2)=0
        ncc(1,2)=nj
        ncc(2,2)=nj+1
        nj=nj+1
        if(nj.gt.njmax)stop'increase njmax!!!'
       else
        if(psran(b10).lt.dc(2))then
         iqc(2)=3*(2.d0*int(.5d0+psran(b10))-1.d0)
	else
         iqc(2)=int(2.d0*psran(b10)+1.d0)
     *   *(2.d0*int(.5d0+psran(b10))-1.d0)
        endif
14      zg=xm+psran(b10)*(1.d0-xm)           !gluon splitting into qq~
        if(psran(b10).gt.zg**dels*((1.d0-xm/zg)/ (1.d0-xm))**betp)
     *  goto 14
        xg=xm/zg
        wmq0=wpm*(xg-xm)
        wpq=1.d0/wmq0
        wpi1=wpi1-wpq
        if(wmi1*wpi1.le.s2min)goto 11
        nj=nj+2
        if(nj.gt.njmax)stop'increase njmax!!!'
	iqj(nj)=-iqc(2)
	if(iabs(iqc(2)).eq.3)iqj(nj)=iqj(nj)*4/3
        eqj(1,nj)=.5d0*wpq
        eqj(2,nj)=.5d0*wpq
        eqj(3,nj)=0.d0
        eqj(4,nj)=0.d0
        if(iqc(2).gt.0)then
         ncj(1,nj)=nj-1
         ncj(1,nj-1)=nj
         ncj(2,nj)=0
         ncj(2,nj-1)=0
         ncc(1,2)=nj-2
         ncc(2,2)=0
        else
         ncj(1,nj)=nj-2
         ncj(1,nj-2)=nj
         ncj(2,nj)=0
         ncj(2,nj-2)=0
         ncc(1,2)=nj-1
         ncc(2,2)=0
        endif
       endif
      endif

      if(jqq.ne.0)then
       if(iqq.ne.0.or.iqq.eq.0.and.jqq.eq.3)then
        sjqq1=qgjit(qt0,qt0,wpi1*wmi1,2,2)
        gbs=sjqq1/sjqq
       else
        sjqg1=qgjit(qt0,qt0,wpi1*wmi1,1,2)
        gbs=sjqg1/sjqg
       endif
       if(psran(b10).gt.gbs)goto 11
      endif
      wpi=wpi1
      wmi=wmi1

      ept(1)=.5d0*(wpi+wmi)      !ladder 4-momentum
      ept(2)=.5d0*(wpi-wmi)
      ept(3)=0.d0
      ept(4)=0.d0
      qmin(1)=qt0                !q^2 cutoff for the upper leg
      qmin(2)=qt0                !q^2 cutoff for the downer leg
      qminn=max(qmin(1),qmin(2)) !overall q^2 cutoff
      si=qgnrm(ept)
      jini=1
      jj=int(1.5d0+psran(b10)) !1st parton at upper (jj=1) or downer (jj=2) leg

3     continue

      aaa=qgnrm(ept)             !ladder mass squared
      if(debug.ge.3)write (moniou,203)si,iqc,ept,aaa

      pt2=ept(3)**2+ept(4)**2
      pt=dsqrt(pt2)
      ww=si+pt2

      iqp(1)=min(1,iabs(iqc(1)))+1
      iqp(2)=min(1,iabs(iqc(2)))+1
      wp(1)=ept(1)+ept(2)                 !LC+ for the ladder
      wp(2)=ept(1)-ept(2)                 !LC- for the ladder
      s2min=4.d0*fqscal*qminn   !minimal energy squared for 2-parton production
      wwmin=(s2min+qminn+pt2-2.d0*pt*dsqrt(qtf))/(1.d0-qtf/qminn)
                                !minimal energy squared for 3-parton production

      if(ww.lt.1.2d0*wwmin)goto 6         !energy too low -> born process

      if(jini.eq.1)then                   !general ladder
       sj=qgjit(qmin(1),qmin(2),si,iqp(1),iqp(2))   !total ladder contribution
       sj1=qgjit1(qmin(3-jj),qmin(jj),si,iqp(3-jj),iqp(jj),2) !one-way ordered
       sjb=0.d0
       if(psran(b10).lt.sj1/sj)then       !change to one-way ordered ladder
        jj=3-jj
	sj=sj1
        sjb=qgbit(qmin(1),qmin(2),si,iqp(1),iqp(2))          !born contribution
	jini=0
       endif
      else                                !one-way ordered ladder
       sj=qgjit1(qmin(jj),qmin(3-jj),si,iqp(jj),iqp(3-jj),2) !one-way ordered
       sjb=qgbit(qmin(1),qmin(2),si,iqp(1),iqp(2))         !born contribution
      endif
      if(debug.ge.3)write (moniou,204)s2min,wwmin,sj,sjb

      if(psran(b10).lt.sjb/sj)goto 6      !born process sampled

c xmin is the minimal parton LC share (for minimal virtuality qmin(jj))
      xxx=1.d0+(2.d0*pt*dsqrt(qtf)-pt2)/ww
      xmin=(xxx+dsqrt(xxx**2-4.d0*qtf*(1.d0+4.d0*fqscal)/ww))/2.d0
      xmax=qtf*(1.d0+4.d0*fqscal)/ww/xmin

      xxx=1.d0+(2.d0*pt*dsqrt(qtf)-pt2-s2min)/ww
      xmin=min(xmin,(xxx+dsqrt(xxx**2-4.d0*qtf/ww))/2.d0)

      xmax=1.d0-xmax
      xmin=1.d0-xmin
      if(debug.ge.3)write (moniou,205)xmin,xmax

      qqmax=qtf/(1.d0-xmax)        !maximal parton virtuality in the current run
      qqmin=qtf/(1.d0-xmin)        !minimal parton virtuality in the current run
      if(debug.ge.3)write (moniou,206)qqmin,qqmax

      if(qqmin.lt.qmin(jj))then
       qqmin=qmin(jj)
       xxx=pt*dsqrt(qqmin)/ww
       xmin=max(xmin,1.d0-(xxx+dsqrt(xxx**2+1.d0
     * -(s2min+qqmin+pt2)/ww))**2)
      endif
      if(debug.ge.3)write (moniou,207)xmin,qqmin

      qm0=qminn
      xm0=1.d0-qtf/qm0
      s2max=xm0*ww

      if(jini.eq.1)then
       sj0=qgjit(qm0,qmin(3-jj),s2max,1,iqp(3-jj))*qgfap(xm0,iqp(jj),1)
     * +qgjit(qm0,qmin(3-jj),s2max,2,iqp(3-jj))*qgfap(xm0,iqp(jj),2)
      else
       sj0=qgjit1(qm0,qmin(3-jj),s2max,1,iqp(3-jj),2)
     * *qgfap(xm0,iqp(jj),1)
     * +qgjit1(qm0,qmin(3-jj),s2max,2,iqp(3-jj),2)*qgfap(xm0,iqp(jj),2)
      endif

      gb0=sj0*qm0*qgalf(qm0/alm)*qgsuds(qm0,iqp(jj)) *4.5d0  !normal. of accept.  
      if(xm0.le..5d0)then
       gb0=gb0*xm0**(1.d0-delh)
      else
       gb0=gb0*(1.d0-xm0)*2.d0**delh
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
4     continue
      if(psran(b10).gt.djl)then
       x=(xmin1+psran(b10)*(xmax1-xmin1))**(1.d0/delh) !parton LC share
      else
       x=1.d0-(1.d0-xmin2)*((1.d0-xmax)/(1.d0-xmin2))**psran(b10)
      endif
      qq=qqmin/(1.d0+psran(b10)*(qqmin/qqmax-1.d0))    !parton virtuality
      qt2=qq*(1.d0-x)                                  !parton p_t^2
      if(debug.ge.4)write (moniou,209)qq,qqmin,qqmax,x,qt2

      if(qt2.lt.qtf)goto 4                 !too collinear parton -> rejection
      if(qq.gt.qminn)then                  !update virtuality cutoff
       qmin2=qq
      else
       qmin2=qminn
      endif
      qt=dsqrt(qt2)
      call qgcs(c,s)
      ep3(3)=qt*c                          !final parton p_x, p_y
      ep3(4)=qt*s
      pt2new=(ept(3)-ep3(3))**2+(ept(4)-ep3(4))**2!p_t^2 for the remained ladder
      s2min2=max(s2min,4.d0*fqscal*qmin2)  !new ladder kinematic limit 
      s2=x*ww-qt2*x/(1.d0-x)-pt2new        !mass squared for the remained ladder
      if(s2.lt.s2min2)goto 4           !ladder mass below threshold -> rejection

      if(jini.eq.1)then                    !weights for g- and q-legs
       sj1=qgjit(qq,qmin(3-jj),s2,1,iqp(3-jj))*qgfap(x,iqp(jj),1)
       sj2=qgjit(qq,qmin(3-jj),s2,2,iqp(3-jj))*qgfap(x,iqp(jj),2)
      else
       sj1=qgjit1(qq,qmin(3-jj),s2,1,iqp(3-jj),2)*qgfap(x,iqp(jj),1)
       sj2=qgjit1(qq,qmin(3-jj),s2,2,iqp(3-jj),2)*qgfap(x,iqp(jj),2)
      endif
      gb7=(sj1+sj2)*qgalf(qq/alm)*qq*qgsuds(qq,iqp(jj))/gb0
                               !acceptance probability for x and q**2 simulation
      if(x.le..5d0)then
       gb7=gb7*x**(1.d0-delh)
      else
       gb7=gb7*(1.d0-x)*2.d0**delh
      endif
      if(debug.ge.4)write (moniou,210)gb7,s2,sj1,sj2,jj,jini
      if(psran(b10).gt.gb7)goto 4          !rejection

c-------------------------------------------------
c define color flow for the emitted jet; perform final state emission
      nqc(2)=0
      if(psran(b10).lt.sj1/(sj1+sj2))then         !new gluon-leg ladder
       if(iqc(jj).eq.0)then                       !g -> gg
        jt=1
        jq=int(1.5d0+psran(b10))
        nqc(1)=ncc(jq,jj)                         !color connection for the jet
        nqc(2)=0
       else                                       !q -> qg
        jt=2
        if(iqc(jj).gt.0)then                      !orientation of color flow
         jq=1
        else
         jq=2
        endif
        nqc(1)=0
        ncc(jq,jj)=ncc(1,jj)                      !color connection for the jet
       endif
       iq1=iqc(jj)                                !jet flavor (type)
       iqc(jj)=0                                  !new ladder leg flavor (type)

      else                                        !new quark-leg ladder
       if(iqc(jj).ne.0)then                       !q -> gq
        iq1=0
        jt=3
        if(iqc(jj).gt.0)then                      !orientation of color flow
         jq=1
        else
         jq=2
        endif
        nqc(1)=ncc(1,jj)                          !color connection for the jet
        nqc(2)=0

       else                                       !g -> qq~
        jq=int(1.5d0+psran(b10))                  !orientation of color flow
        iq1=int(3.d0*psran(b10)+1.d0)*(3-2*jq)    !jet flavor (type)
        iqc(jj)=-iq1                              !new ladder leg flavor (type)
        jt=4
        nqc(1)=ncc(jq,jj)                         !color connections for the jet
        ncc(1,jj)=ncc(3-jq,jj)
       endif
      endif
      if(debug.ge.3)write (moniou,211)jt

      call qgcjet(qt2,iq1,qv1,zv1,qm1,iqv1,ldau1,lpar1,jq) !final state emission
      si=x*ww-(qt2+qm1(1,1))*x/(1.d0-x)-pt2new  !mass squared for the new ladder
      if(si.gt.s2min2)then
       iq=min(1,iabs(iqc(jj)))+1
       if(jini.eq.1)then
        gb=qgjit(qq,qmin(3-jj),si,iq,iqp(3-jj))
     *  /qgjit(qq,qmin(3-jj),s2,iq,iqp(3-jj))
       else
        gb=qgjit1(qq,qmin(3-jj),si,iq,iqp(3-jj),2)
     *  /qgjit1(qq,qmin(3-jj),s2,iq,iqp(3-jj),2)
       endif
       if(psran(b10).gt.gb)goto 1        !jet mass correction for the acceptance
      else                                        !below threshold -> rejection
       goto 1
      endif

      wp3=wp(jj)*(1.d0-x)
      wm3=(qt2+qm1(1,1))/wp3
      ep3(1)=.5d0*(wp3+wm3)                       !jet 4-momentum
      ep3(2)=.5d0*(wp3-wm3)*(3-2*jj)
      call qgrec(ep3,nqc,qv1,zv1,qm1,iqv1,ldau1,lpar1,jq)
                               !reconstruction of 4-momenta of all final partons
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
       ept(i)=ept(i)-ep3(i)                       !new ladder 4-momentum
      enddo
      qmin(jj)=qq                                 !new virtuality cutoffs
      qminn=qmin2
      goto 3                                      !consider next parton emission

c------------------------------------------------
c born process - last parton pair production in the ladder
6     continue
      if(debug.ge.2)write (moniou,214)si,qminn,iqc
      tmin=qminn*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qminn*fqscal/si)))
      qtmin=tmin*(1.d0-tmin/si)
      if(iqc(1).ne.0.or.iqc(2).ne.0)then
       gb0=tmin**2*qgalf(qtmin/fqscal/alm)**2
     * *qgfbor(si,tmin,iqc(1),iqc(2),1)
      else
       gb0=.25d0*si**2*qgalf(qtmin/fqscal/alm)**2
     * *qgfbor(si,.5d0*si,iqc(1),iqc(2),1)
      endif
      gb0=gb0*qgsuds(qtmin/fqscal,iqp(1))*qgsuds(qtmin/fqscal,iqp(2))
                                                    !normalization of acceptance	
      if(debug.ge.3)write (moniou,215)gb0
     
7     q2=tmin/(1.d0-psran(b10)*(1.d0-2.d0*tmin/si))   !proposed q^2
      z=q2/si                                         !parton LC momentum share
      qt2=q2*(1.d0-z)                                 !parton p_t^2
      if(psran(b10).lt..5d0)then
       jm=2
       tq=si-q2
      else
       jm=1
       tq=q2
      endif
      gb=q2**2*qgalf(qt2/fqscal/alm)**2*qgfbor(si,tq,iqc(1),iqc(2),1)
     **qgsuds(qt2/fqscal,iqp(1))*qgsuds(qt2/fqscal,iqp(2))/gb0
                                                      !acceptance probabilty
      if(debug.ge.4)write (moniou,216)gb,q2,z,qt2     
      if(psran(b10).gt.gb)goto 7                      !rejection

c-------------------------------------------------
c define color connections for the 1st emitted jet
      nqc(2)=0
      if(iqc(1).eq.0.and.iqc(2).eq.0)then             !gg-process
       jq=int(1.5d0+psran(b10))                       !orientation of color flow
       nqc(1)=ncc(jq,jm)

       if(psran(b10).lt..5d0)then
        jt=1                                          !gg -> gg
        nqc(2)=0
        njc1=ncc(3-jq,jm)                         !color connections for 1st jet
        njc2=ncc(jq,3-jm)
        if(ncc(1,1).eq.0.and.ncc(2,1).eq.0)then
         if(jm.eq.1)nqc(1)=njc2
        else
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
        endif
       else                                 !gg -> gg (inverse color connection)
        jt=2
        nqc(2)=ncc(3-jq,3-jm)
       endif

      elseif(iqc(1)*iqc(2).eq.0)then                  !qg -> qg
       if(iqc(1)+iqc(2).gt.0)then                     !orientation of color flow
        jq=1
       else
        jq=2
       endif
       if(psran(b10).lt..5d0)then
        if(iqc(jm).eq.0)then
         jt=3
         nqc(1)=ncc(jq,jm)
         nqc(2)=0
         njc1=ncc(3-jq,jm)
         njc2=ncc(1,3-jm)
         if(ncc(1,jm).eq.0.and.ncc(2,jm).eq.0)then
          nqc(1)=njc2
         else
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
      
      elseif(iqc(1)*iqc(2).gt.0)then                  !qq (q~q~) -> qq (q~q~)
       jt=7
       if(iqc(1).gt.0)then
        jq=1
       else
        jq=2
       endif
       nqc(1)=ncc(1,3-jm)
      else                                            !qq~ -> qq~
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
      if(jt.ne.8)then
       jq2=jq
      else
       jq2=3-jq
      endif
      if(debug.ge.3)write (moniou,211)jt
      call qgcjet(qt2,iqc(jm),qv1,zv1,qm1,iqv1,ldau1,lpar1,jq)!final state emis.
      call qgcjet(qt2,iqc(3-jm),qv2,zv2,qm2,iqv2,ldau2,lpar2,jq2)
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
      ep3(1)=.5d0*(wp3+wm3)                        !1st jet 4-momentum
      ep3(2)=.5d0*(wp3-wm3)
      qt=dsqrt(qt2)
      call qgcs(c,s)
      ep3(3)=qt*c
      ep3(4)=qt*s

      call qgtran(ep3,ey,1)
      call qgrec(ep3,nqc,qv1,zv1,qm1,iqv1,ldau1,lpar1,jq)
                               !reconstruction of 4-momenta of all final partons
      if(iabs(iqc(jm)).eq.3)then
       iqqq=8+iqc(jm)/3*4
      else
       iqqq=8+iqc(jm)
      endif
      if(debug.ge.3)write (moniou,212)tyq(iqqq),qt2,ep3

      wp3=(1.d0-z)*dsqrt(si)
      wm3=(qt2+qm2(1,1))/wp3
      ep3(1)=.5d0*(wp3+wm3)                        !2nd jet 4-momentum
      ep3(2)=.5d0*(wp3-wm3)
      ep3(3)=-qt*c
      ep3(4)=-qt*s
      call qgtran(ep3,ey,1)

c-------------------------------------------------
c define color connections for the 2nd emitted jet
      if(jt.eq.1)then
       nqc(1)=nqc(2)
       if(ncc(1,3-jm).eq.0.and.ncc(2,3-jm).eq.0)then
        nqc(2)=ncc(3-jq,jm)
       else
        nqc(2)=ncc(3-jq,3-jm)
       endif
      elseif(jt.eq.2)then
       if(ncc(1,1).eq.0.and.ncc(2,1).eq.0)then
        if(jm.eq.1)then
         nqc(2)=nqc(1)
         nqc(1)=ncc(jq,3-jm)
        else
         nqc(1)=nqc(2)
         nqc(2)=ncc(3-jq,jm)
        endif
       else
        nqc(2)=ncc(3-jq,jm)
        nqc(1)=ncc(jq,3-jm)
       endif
      elseif(jt.eq.3)then
       nqc(1)=nqc(2)
      elseif(jt.eq.4)then
       nqc(2)=nqc(1)
       if(ncc(1,1).eq.0.and.ncc(2,1).eq.0)then
        nqc(1)=ncc(1,jm)
       else
        nqc(1)=ncc(jq,3-jm)
       endif
      elseif(jt.eq.5)then
       if(ncc(1,jm).eq.0.and.ncc(2,jm).eq.0)then
        nqc(1)=nqc(2)
       else
        nqc(1)=ncc(jq,jm)
       endif
      elseif(jt.eq.6)then
       if(ncc(1,3-jm).eq.0.and.ncc(2,3-jm).eq.0)then
        nqc(2)=nqc(1)
       else
        nqc(2)=ncc(3-jq,3-jm)
       endif
       nqc(1)=ncc(1,jm)
      elseif(jt.eq.7)then
       nqc(1)=ncc(1,jm)
      endif
      call qgrec(ep3,nqc,qv2,zv2,qm2,iqv2,ldau2,lpar2,jq2)
                               !reconstruction of 4-momenta of all final partons
      if(iabs(iqc(3-jm)).eq.3)then
       iqqq=8+iqc(3-jm)/3*4
      else
       iqqq=8+iqc(3-jm)
      endif
      if(debug.ge.3)write (moniou,212)tyq(iqqq),qt2,ep3

      ebal(1)=.5d0*(wpp+wpm)                          !balans of 4-momentum
      ebal(2)=.5d0*(wpp-wpm)
      ebal(3)=0.d0
      ebal(4)=0.d0
      do i=nj0+1,nj
       if(iqq.eq.0.or.iqq.eq.1.and.i.ne.nva.or.iqq.eq.2
     * .and.i.ne.nvb.or.iqq.eq.3.and.i.ne.nva.and.i.ne.nvb)then
        do j=1,4
         ebal(j)=ebal(j)-eqj(j,i)
        enddo
       endif
      enddo
      if(ebal(1).gt.1.d-5)write (*,*)'ebal',ebal,iqq,jqq,nva,nvb     
      if(debug.ge.2)write (moniou,218)nj
      if(debug.ge.5)write (moniou,219)ebal
      if(debug.ge.1)write (moniou,220)

201   format(2x,'qghot - semihard interaction:'/
     *4x,'type of the interaction - ',i2/
     *4x,'initial light cone momenta - ',2e10.3/
     *4x,'remnant types - ',2i3,2x,'diffr. eigenstates - ',2i2/
     *4x,'proj. class - ',i2,2x,'Pomeron type - ',i2/
     *4x,'initial number of final partons - ',i4)
202   format(2x,'qghot: mass squared for parton ladder - ',e10.3)
203   format(2x,'qghot: ',' mass squared for the laddder:',e10.3/
     *4x,'ladder end flavors:',2i3/4x,'ladder 5-momentum: ',5e10.3)
204   format(2x,'qghot: kinematic bounds s2min=',e10.3,
     *2x,'wwmin=',e10.3/4x,'jet cross section sj=',e10.3,
     *2x,'born cross section sjb=',e10.3)
205   format(2x,'qghot: xmin=',e10.3,2x,'xmax=',e10.3)
206   format(2x,'qghot: qqmin=',e10.3,2x,'qqmax=',e10.3)
207   format(2x,'qghot: xmin=',e10.3,2x,'qqmin=',e10.3)
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

c=============================================================================
      double precision function qgbit(qi,qj,s,m,l)
c------------------------------------------------------------------------
c qgbit - born cross-section interpolation
c qi,qj - virtuality cutoffs for the scattering,
c s     - total c.m. energy squared for the scattering,
c m     - parton type at current end of the ladder (1 - g, 2 - q)
c l     - parton type at opposite end of the ladder (1 - g, 2 - q)
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wi(3),wk(3)
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr20/ spmax
      common /qgarr26/ factk,fqscal
      common /qgarr31/ csj(40,160)
      common /qgarr43/ moniou
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)qi,qj,s,m,l
      qgbit=0.d0
      qq=max(qi,qj)
      s2min=qq*4.d0*fqscal
      if(s.le..99d0*s2min)then
       if(debug.ge.4)write (moniou,202)qgbit
       return
      endif

      tmin=qq*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qq*fqscal/s)))
      ml=40*(m-1)+80*(l-1)
      qli=dlog(qq)/dlog(spmax/4.d0/fqscal)*39.d0+1.d0
      sl=dlog(s/s2min)/dlog(spmax/s2min)*39.d0+1.d0
      i=min(38,int(qli))
      k=min(38,int(sl))

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
      if(qi.lt.qq)qgbit=qgbit*qgsuds(qq,m)/qgsuds(qi,m)
      if(qj.lt.qq)qgbit=qgbit*qgsuds(qq,l)/qgsuds(qj,l)
      if(debug.ge.4)write (moniou,202)qgbit
      
201   format(2x,'qgbit - born cross-section interpolation:'
     */4x,'qi=',e10.3,2x,'qj=',e10.3,2x,'s= ',e10.3,2x,'m= ',i1
     *,2x,'l= ',i1)
202   format(2x,'qgbit=',e10.3)
      return
      end

c=============================================================================
      double precision function qgfbor(s,t,iq1,iq2,n)
c---------------------------------------------------------------------------
c qgfbor - integrand for the born cross-section (matrix element squared)
c s   - total c.m. energy squared for the scattering,
c t   - invariant variable for the scattering abs[(p1-p3)**2],
c iq1 - parton type at current end of the ladder (0 - g, 1,2 - q)
c iq2 - parton type at opposite end of the ladder (0 - g, 1,2 - q)
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)s,t,iq1,iq2
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
      if(debug.ge.4)write (moniou,202)qgfbor
      
201   format(2x,'qgfbor - hard scattering matrix element squared:'/
     *4x,'s=',e10.3,2x,'|t|=',e10.3,2x,'iq1=',i1,2x,'iq2=',i1)
202   format(2x,'qgfbor=',e10.3)
      return
      end

c=============================================================================
      double precision function qgborn(qi,qj,s,iq1,iq2,jj)
c-----------------------------------------------------------------------------
c qgborn - hard 2->2 parton scattering born cross-section
c s   - c.m. energy squared for the scattering process,
c iq1 - parton type at current end of the ladder (0 - g, 1,2 etc. - q)
c iq2 - parton type at opposite end of the ladder (0 - g, 1,2 etc. - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr26/ factk,fqscal
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)qi,qj,s,iq1,iq2
      
      qgborn=0.d0
      qq=max(qi,qj)
      tmin=qq*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qq*fqscal/s)))
      do i=1,7
      do m=1,2
       t=2.d0*tmin/(1.d0+2.d0*tmin/s-x1(i)*(2*m-3)*(1.d0-2.d0*tmin/s))
       qt=t*(1.d0-t/s)
       
       fb=0.d0
       do n=1,3
        fb=fb+qgfbor(s,t,iq1,iq2,n)+qgfbor(s,s-t,iq1,iq2,n)
       enddo
       if(jj.eq.1)then
        fb=fb*qgsudx(qt/fqscal,iabs(iq1)+1)
     *  *qgsudx(qt/fqscal,iabs(iq2)+1)
       elseif(jj.eq.2)then
        fb=fb*qgsuds(qt/fqscal,iabs(iq1)+1)
     *  *qgsuds(qt/fqscal,iabs(iq2)+1)
       elseif(jj.eq.3)then
        fb=fb*qgsudx(qt/fqscal,iabs(iq1)+1)
     *  *qgsuds(qt/fqscal,iabs(iq2)+1)
       endif
       qgborn=qgborn+a1(i)*fb*qgalf(qt/fqscal/alm)**2*t**2
      enddo
      enddo
      qgborn=qgborn*2.d0*pi**3/s**2
      if(jj.eq.1)then
       qgborn=qgborn/qgsudx(qi,iabs(iq1)+1)/qgsudx(qj,iabs(iq2)+1)
      elseif(jj.eq.2)then
       qgborn=qgborn/qgsuds(qi,iabs(iq1)+1)/qgsuds(qj,iabs(iq2)+1)
      elseif(jj.eq.3)then
       qgborn=qgborn/qgsudx(qi,iabs(iq1)+1)/qgsuds(qj,iabs(iq2)+1)
      endif
      if(iq1.eq.iq2)qgborn=qgborn*.5d0
      if(debug.ge.4)write (moniou,202)qgborn
      
201   format(2x,'qgborn - 2->2 parton scattering cross-section:'
     */4x,'qi=',e10.3,2x,'qj=',e10.3,2x,'s= ',e10.3,2x,'iq1= ',i1
     *,2x,'iq2= ',i1)
202   format(2x,'qgborn=',e10.3)
      return
      end

c=============================================================================
      subroutine qgcjet(qq,iq1,qv,zv,qm,iqv,ldau,lpar,jq)
c-----------------------------------------------------------------------------
c qgcjet - final state emission process
c (all branchings as well as parton masses are determined)
c qq  - maximal virtuality for the first branching,
c iq1 - initial parton type (flavour),
c jq  - orientation of the color flow (1,2)
c to determine:
c qv(i,j)  - parton branching virt. in i-th row on j-th level (0 - no branching)
c zv(i,j)  - z-value for parton branching in i-th row on j-th level
c qm(i,j)  - virtual parton mass squared in i-th row on j-th level
c iqv(i,j) - parton type (flavour) in i-th row on j-th level
c ldau(i,j) - first daughter row for parton branching in i-th row on j-th level
c lpar(i,j) - parent row for the parton in i-th row on j-th level
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension qmax(30,50),iqm(2),lnv(50),
     *qv(30,50),zv(30,50),qm(30,50),iqv(30,50),
     *ldau(30,49),lpar(30,50)
      common /qgarr11/ b10
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr43/ moniou
      common /debug/   debug
      external psran

      if(debug.ge.3)write (moniou,201)qq,iq1,jq
      do i=2,20
       lnv(i)=0
      enddo
      lnv(1)=1
      qmax(1,1)=qq
      iqv(1,1)=iq1
      nlev=1
      nrow=1

2     qlmax=dlog(qmax(nrow,nlev)/qtf/16.d0)
      iq=min(1,iabs(iqv(nrow,nlev)))+1

      if(psran(b10).gt.qgsudi(qlmax,iq))then
       q=qgqint(qlmax,psran(b10),iq)
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
        if(psran(b10).lt.wg)then
         iqm(1)=0
         iqm(2)=0
        else
         iqm(1)=int(3.d0*psran(b10)+1.d0)*(3-2*jq)
         iqm(2)=-iqm(1)
        endif
        if(psran(b10).lt..5d0)z=1.d0-z
       endif
       qv(nrow,nlev)=q
       zv(nrow,nlev)=z
       nrow=ll
       nlev=nlev+1
       qmax(nrow,nlev)=q*z**2
       qmax(nrow+1,nlev)=q*(1.d0-z)**2
       iqv(nrow,nlev)=iqm(1)
       iqv(nrow+1,nlev)=iqm(2)
       if(debug.ge.4)write (moniou,202)nlev,nrow,q,z   
       goto 2
      else
       qv(nrow,nlev)=0.d0
       zv(nrow,nlev)=0.d0
       qm(nrow,nlev)=0.d0
       if(debug.ge.4)write (moniou,203)nlev,nrow
      endif

3     continue
      if(nlev.eq.1)then
       if(debug.ge.4)write (moniou,205)
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
       if(debug.ge.4)write (moniou,204)nlev,nrow,qm(lprow,nlev)
       goto 3
      endif
      
201   format(2x,'qgcjet - final state emission process:'
     */4x,'qq=',e10.3,2x,'iq1= ',i1,2x,'jq=',i1)
202   format(2x,'qgcjet: new branching at level nlev=',i2,' nrow=',i2
     */4x,' effective momentum q=',e10.3,2x,' z=',e10.3)
203   format(2x,'qgcjet: new final jet at level nlev=',i2,' nrow=',i2)
204   format(2x,'qgcjet: jet mass at level nlev=',i2,' nrow=',i2
     *,' - qm=',e10.3)
205   format(2x,'qgcjet - end')
      end

c===========================================================================
      subroutine qgcs(c,s)
c---------------------------------------------------------------------------
c qgcs - cos and sin generation for uniformly distributed angle 0<fi<2*pi
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /debug/   debug
      external psran

      if(debug.ge.3)write (moniou,201)    
1     s1=2.d0*psran(b10)-1.d0
      s2=2.d0*psran(b10)-1.d0
      s3=s1*s1+s2*s2
      if(s3.gt.1.d0)goto 1
      s3=dsqrt(s3)
      c=s1/s3
      s=s2/s3      
      if(debug.ge.4)write (moniou,202)c,s
      
201   format(2x,'qgcs - cos(fi) and sin(fi) are generated',
     *' (0<fi<2*pi)')
202   format(2x,'qgcs: c=',e10.3,2x,'s=',e10.3)
      return
      end

c===========================================================================
      subroutine qgdeft(s,ep,ey)
c---------------------------------------------------------------------------
c qgdeft - Lorentz boost to the rest frame for 4-vector ep
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ey(3),ep(4)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)ep,s
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
      if(debug.ge.4)write (moniou,202)ey
      
201   format(2x,'qgdeft - lorentz boost parameters:'
     */4x,'4-vector ep=',4e10.3/4x,'4-vector squared s=',e10.3)
202   format(2x,'qgdeft: lorentz boost parameters ey(i)=',2x,3e10.3)
      return
      end

c=============================================================================
      subroutine qgdefr(ep,s0x,c0x,s0,c0)
c-----------------------------------------------------------------------------
c qgdefr - spacial rotation to z-axis for 4-vector ep
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ep(4)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)ep
      pt2=ep(3)**2+ep(4)**2               !squared p_t
      if(pt2.ne.0.d0)then
       pt=dsqrt(pt2)
       c0x=ep(3)/pt                       !cos phi
       s0x=ep(4)/pt                       !sin phi
       pl=dsqrt(pt2+ep(2)**2)             !total momentum
       s0=pt/pl                           !sin theta
       c0=ep(2)/pl                        !cos theta
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
      if(debug.ge.4)write (moniou,202)s0x,c0x,s0,c0,ep
      
201   format(2x,'qgdefr - spacial rotation parameters'/4x,
     *'4-vector ep=',2x,4(e10.3,1x))
202   format(2x,'qgdefr: spacial rotation parameters'/
     *4x,'s0x=',e10.3,2x,'c0x=',e10.3,2x,'s0=',e10.3,2x,'c0=',e10.3/
     *4x,'rotated 4-vector ep=',4(e10.3,1x))
      return
      end
	
c=============================================================================
      double precision function qgfap(x,j,l)
c------------------------------------------------------------------------
c qgfap - Altarelli-Parisi kernel (multiplied by x)
c x - light cone momentum share,
c j - type of the parent parton (1-g,2-q),
c l - type of the daughter parton (1-g,2-q)
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)x,j,l
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
      if(debug.ge.4)write (moniou,202)qgfap
      
201   format(2x,'qgfap - Altarelli-Parisi kernel:'
     *,2x,'x=',e10.3,2x,'j=',i1,2x,'l=',i1)
202   format(2x,'qgfap=',e10.3)
      return
      end

c=============================================================================
      subroutine qggea(ia,xa,jj)
c-----------------------------------------------------------------------------
c qggea - nuclear configuration (nucleon positions)
c ia - nuclear mass number,
c jj=1 - projectile
c jj=2 - target
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207)
      dimension xa(iapmax,3)
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /debug/   debug
      external psran

      if(debug.ge.3)write (moniou,201)jj,ia
      if(ia.ge.10)then                 !light nucleus -> Van-Hove method
       do i=1,ia
1       zuk=psran(b10)*cr1(jj)-1.d0
        if(zuk)2,2,3
2       tt=rnuc(jj)/wsnuc(jj)*(psran(b10)**.3333d0-1.d0)
        goto 6
3       if(zuk.gt.cr2(jj))goto 4
        tt=-dlog(psran(b10))
        goto 6
4       if(zuk.gt.cr3(jj))goto 5
        tt=-dlog(psran(b10))-dlog(psran(b10))
        goto 6
5       tt=-dlog(psran(b10))-dlog(psran(b10))-dlog(psran(b10))
6       rim=tt*wsnuc(jj)+rnuc(jj)
        if(psran(b10).gt.(1.d0+wbnuc(jj)*rim**2/rnuc(jj)**2)
     *  /(1.d0+exp(-abs(tt))))goto 1
        z=rim*(2.d0*psran(b10)-1.d0)
        rim=dsqrt(rim*rim-z*z)
        xa(i,3)=z
        call qgcs(c,s)
        xa(i,1)=rim*c
        xa(i,2)=rim*s
       enddo
       
      else                               !heavy nucleus -> 3-parameter Fermi
       do l=1,3
        summ=0.d0
        do i=1,ia-1
         j=ia-i
         aks=rnuc(jj)*(psran(b10)+psran(b10)+psran(b10)-1.5d0)
         k=j+1
         xa(k,l)=summ-aks*sqrt(float(j)/k)
         summ=summ+aks/sqrt(float(j*k))
        enddo
        xa(1,l)=summ
       enddo
      endif
      if(debug.ge.5)then
       write (moniou,202)
       do i=1,ia
        write (moniou,203)i,(xa(i,l),l=1,3)
       enddo
      endif
      if(debug.ge.4)write (moniou,204)
      
201   format(2x,'qggea - configuration of the nucleus ',i1,';',2x,
     *'coordinates for ',i2,' nucleons')
202   format(2x,'qggea:  positions of the nucleons')
203   format(2x,'qggea: ',i2,' - ',3(e10.3,1x))
204   format(2x,'qggea - end')
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
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)x,j,l
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
      if(debug.ge.4)write (moniou,202)qgapi
      
201   format(2x,'qgapi - integrated AP kernel: x=',e10.3
     *,2x,'j= ',i1,2x,'l= ',i1)
202   format(2x,'qgapi=',e10.3)
      return
      end

c=============================================================================
      subroutine qgjarr(jfl)
c-----------------------------------------------------------------------------
c qgjarr - final jets rearrangement according to their colour connections
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(njmax=50000)
      dimension mark(njmax),ept(4)
      common /qgarr10/ am(7),ammu
      common /qgarr36/ epjet(4,njmax),ipjet(njmax),njtot
      common /qgarr37/ eqj(4,njmax),iqj(njmax),ncj(2,njmax),nj
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)nj
      if(debug.ge.3.and.nj.ne.0)then
       do i=1,nj
        write (moniou,202)i,iqj(i),(eqj(l,i),l=1,4)
        if(iqj(i).eq.0)then
         write (moniou,203)ncj(1,i),ncj(2,i)
        else
         ncdum=0
         write (moniou,203)ncj(1,i),ncdum
        endif
       enddo
      endif

      jfl=0
      do i=1,nj
       mark(i)=1
      enddo
      njtot=0

2     continue
      do ij=1,nj
       if(mark(ij).ne.0.and.iqj(ij).ne.0)goto 4
      enddo
4     continue

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

6     mark(ij)=0
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
         if(debug.ge.3)write (moniou,204)jfl
         return
        endif

        if(njtot.lt.nj)then
         goto 2
        else
         jfl=1
         nj=0
         if(debug.ge.3)write (moniou,204)jfl
         return
        endif
       else
        jfirst=0
        njpar=ij
        ij=ncj(1,ij)
        goto 6
       endif
      else
       if(ncj(1,ij).eq.njpar)then
        njdau=ncj(2,ij)
       else
        njdau=ncj(1,ij)
       endif
       njpar=ij
       ij=njdau
       goto 6
      endif
      
201   format(2x,'qgjarr - color arrangement: total number of jets nj='
     *,i4)
202   format(2x,'qgjarr: ij=',i3,2x,'iqj=',i2,2x,'eqj=',4e10.3)
203   format(2x,'qgjarr: ncj=',2i3)
204   format(2x,'qgjarr - end, return flag: ',i2)
      end

c=============================================================================
      double precision function qgjet(q1,q2,s,s2min,j,l)
c-----------------------------------------------------------------------------
c qgjet - inclusive ladder cross-section (one more rung added)
c q1    - virtuality cutoff for current end of the ladder,
c q2    - virtuality cutoff for opposite end of the ladder,
c s     - total c.m. energy squared for the ladder,
c s2min - minimal c.m. energy squared for born process (above q1 and q2),
c j     - parton type at current end of the ladder (1 - g, 2 - q),
c l     - parton type at opposite end of the ladder (1 - g, 2 - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr26/ factk,fqscal
      common /qgarr32/ epsxmn
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)s,q1,q2,s2min,j,l
      qgjet=0.d0
      qmax=s/4.d0/fqscal*(1.d0-epsxmn)
      qmin=q1
      if(debug.ge.5)write (moniou,202)qmin,qmax

      if(qmax.gt.qmin)then
       do i=1,7
       do m=1,2
        qi=2.d0*qmin/(1.d0+qmin/qmax+(2*m-3)*x1(i)*(1.d0-qmin/qmax))
        zmax=(1.d0-epsxmn)**delh
        zmin=(max(4.d0*fqscal*qi,s2min)/s)**delh
        fsj=0.d0	  
        if(debug.ge.5)write (moniou,203)qi,zmin,zmax

        if(zmax.gt.zmin)then
         do i1=1,7
         do m1=1,2
          z=(.5d0*(zmax+zmin+(2*m1-3)*x1(i1)*(zmax-zmin)))**(1.d0/delh)
          s2=z*s
          sj=0.d0
          do k=1,2
           sj=sj+qgjit(qi,q2,s2,k,l)*qgfap(z,j,k)*z
          enddo
          fsj=fsj+a1(i1)*sj/z**delh
         enddo
         enddo
         fsj=fsj*(zmax-zmin)
        endif
        qgjet=qgjet+a1(i)*fsj*qi*qgsudx(qi,j)*qgalf(qi/alm)
       enddo
       enddo
       qgjet=qgjet*(1.d0/qmin-1.d0/qmax)/qgsudx(q1,j)/delh/4.d0
      endif     
      if(debug.ge.4)write (moniou,204)qgjet
      
201   format(2x,'qgjet - general ladder cross section:'
     */4x,'s=',e10.3,2x,'q1=',e10.3,2x,'q2=',e10.3,2x,'s2min=',
     *e10.3,2x,'j=',i1,2x,'l=',i1)
202   format(2x,'qgjet:',2x,'qmin=',e10.3,2x,'qmax=',e10.3)
203   format(2x,'qgjet:',2x,'qi=',e10.3,2x,'zmin=',e10.3
     *,2x,'zmax=',e10.3)
204   format(2x,'qgjet=',e10.3)
      return
      end

c=============================================================================
      double precision function qgjet1(q1,q2,s,s2min,j,l,jj)
c-----------------------------------------------------------------------------
c qgjet1 - inclusive one-way ordered ladder cross-section (one more rung added)
c q1    - virtuality cutoff for current end of the ladder,
c q2    - virtuality cutoff for opposite end of the ladder,
c s     - total c.m. energy squared for the ladder,
c s2min - minimal c.m. energy squared for born process (above q1 and q2),
c j     - parton type at current end of the ladder (1 - g, 2 - q),
c l     - parton type at opposite end of the ladder (1 - g, 2 - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr26/ factk,fqscal
      common /qgarr32/ epsxmn
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)s,q1,q2,s2min,j,l
      qgjet1=0.d0
      qmax=s/4.d0/fqscal*(1.d0-epsxmn)
      qmin=q1
      if(debug.ge.5)write (moniou,202)qmin,qmax

      if(qmax.gt.qmin)then
       do i=1,7
       do m=1,2
        qi=2.d0*qmin/(1.d0+qmin/qmax+(2*m-3)*x1(i)*(1.d0-qmin/qmax))
        zmax=(1.d0-epsxmn)**delh
        zmin=(max(4.d0*fqscal*qi,s2min)/s)**delh
        fsj=0.d0
        if(debug.ge.5)write (moniou,203)qi,zmin,zmax

        if(zmax.gt.zmin)then
         do i1=1,7
         do m1=1,2
          z=(.5d0*(zmax+zmin+(2*m1-3)*x1(i1)*(zmax-zmin)))**(1.d0/delh)
          s2=z*s

          sj=0.d0
          do k=1,2
           sj=sj+qgjit1(qi,q2,s2,k,l,jj)*qgfap(z,j,k)*z
          enddo
          fsj=fsj+a1(i1)*sj/z**delh
         enddo
         enddo
         fsj=fsj*(zmax-zmin)
        endif
        qgjet1=qgjet1+a1(i)*fsj*qi*qgsudx(qi,j)*qgalf(qi/alm)
       enddo
       enddo
       qgjet1=qgjet1*(1.d0/qmin-1.d0/qmax)/qgsudx(q1,j)/delh/4.d0
      endif
      if(debug.ge.4)write (moniou,204)qgjet1
      
201   format(2x,'qgjet1 - one-way ordered ladder cross section:'
     */4x,'s=',e10.3,2x,'q1=',e10.3,2x,'q2=',e10.3,2x,'s2min=',
     *e10.3,2x,'j=',i1,2x,'l=',i1)
202   format(2x,'qgjet1:',2x,'qmin=',e10.3,2x,'qmax=',e10.3)
203   format(2x,'qgjet1:',2x,'qi=',e10.3,2x,'zmin=',e10.3
     *,2x,'zmax=',e10.3)
204   format(2x,'qgjet1=',e10.3)
      return
      end

c=============================================================================
      double precision function qgjit(q1,q2,s,m,l)
c-----------------------------------------------------------------------------
c qgjit - interpolation of general ladder cross-section
c q1    - virtuality cutoff for current end of the ladder,
c q2    - virtuality cutoff for opposite end of the ladder,
c s     - total c.m. energy squared for the ladder,
c m - parton type at current end of the ladder (1 - g, 2 - q),
c l - parton type at opposite end of the ladder (1 - g, 2 - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wi(3),wj(3),wk(3)
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr20/ spmax
      common /qgarr26/ factk,fqscal
      common /qgarr29/ csj(40,40,160)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)s,q1,q2,m,l
      qgjit=0.d0
      qq=max(q1,q2)
      s2min=qq*4.d0*fqscal
      if(s.le..99d0*s2min)then
       if(debug.ge.4)write (moniou,202)qgjit
       return
      endif

      if(q1.le.q2)then
       qi=q1
       qj=q2
       ml=40*(m-1)+80*(l-1)
      else
       qi=q2
       qj=q1
       ml=40*(l-1)+80*(m-1)
      endif

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
      if(debug.ge.4)write (moniou,202)qgjit
      
201   format(2x,'qgjit - interpolation of general ladder cross-section:'
     */4x,'s=',e10.3,2x,'q1=',e10.3,2x,'q2=',e10.3,2x,2x,'m=',i1
     *,2x,'l=',i1)
202   format(2x,'qgjit=',e10.3)
      return
      end

c=============================================================================
      double precision function qgjit1(q1,q2,s,m,l,jj)
c-----------------------------------------------------------------------------
c qgjit1 - interpolation of one-way ordered ladder cross-section
c q1    - virtuality cutoff for current end of the ladder,
c q2    - virtuality cutoff for opposite end of the ladder,
c s     - total c.m. energy squared for the ladder,
c m - parton type at current end of the ladder (1 - g, 2 - q),
c l - parton type at opposite end of the ladder (1 - g, 2 - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wi(3),wj(3),wk(3)
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr20/ spmax
      common /qgarr26/ factk,fqscal
      common /qgarr30/ csj(40,40,160,2)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)s,q1,q2,m,l
      qgjit1=0.d0
      qq=max(q1,q2)
      s2min=qq*4.d0*fqscal
      if(s.le.s2min)then
       if(debug.ge.4)write (moniou,202)qgjit1
       return
      endif

      tmin=qq*fqscal/(.5d0+dsqrt(max(0.d0,.25d0-qq*fqscal/s)))
      ml=40*(m-1)+80*(l-1)
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
       qgjit1=qgjit1+csj(i+i1-1,j+j1-1,k2,jj)*wi(i1)*wj(j1)*wk(k1)
      enddo
      enddo
      enddo
      qgjit1=exp(qgjit1)*(1.d0/tmin-2.d0/s)
      if(q2.lt.q1)then
       if(jj.eq.1)then
        qgjit1=qgjit1*qgsudx(q1,l)/qgsudx(q2,l)
       else
        qgjit1=qgjit1*qgsuds(q1,l)/qgsuds(q2,l)
       endif
      endif
      if(debug.ge.4)write (moniou,202)qgjit1
      
201   format(2x,'qgjit1 - one-way ordered ladder cross section:'/4x,
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
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)s,a,b     
      qglam=.25d0/s*(s+a-b)**2-a 
           
      if(debug.ge.4)write (moniou,202)qglam      
201   format(2x,'qglam - kinematical function, s=',e10.3,2x,'a='
     *,e10.3,2x,'b=',e10.3)
202   format(2x,'qglam=',e10.3)
      return
      end

c=============================================================================
      double precision function qgnrm(ep)
c-----------------------------------------------------------------------------
c invariant mass for the 4-vector ep
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ep(4)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)ep
      qgnrm=(ep(1)-ep(2))*(ep(1)+ep(2))-ep(3)**2-ep(4)**2
      
      if(debug.ge.4)write (moniou,202)qgnrm      
201   format(2x,'qgnrm - invariant mass for ep=',4(e10.3,1x))
202   format(2x,'qgnrm=',e10.3)
      return
      end

c===========================================================================
      subroutine qgrec(ep,nqc,qv,zv,qm,iqv,ldau,lpar,jq)
c---------------------------------------------------------------------------
c qgrec - jet reconstructing (4-momenta for all final partons)
c ep(i)    - parent parton 4-momentum
c nqc      - color connections of the parent parton
c qv(i,j)  - parton branching virt. in i-th row on j-th level (0 - no branching)
c zv(i,j)  - z-value for parton branching in i-th row on j-th level
c qm(i,j)  - virtual parton mass squared in i-th row on j-th level
c iqv(i,j) - parton type (flavour) in i-th row on j-th level
c ldau(i,j) - first daughter row for parton branching in i-th row on j-th level
c lpar(i,j) - parent row for the parton in i-th row on j-th level
c jq        - orientation of the color flow (1,2)
c to determine:
c eqj(i,nj) - 4-momentum for the final jet nj
c iqj(nj)   - type (flavour) for the final jet nj
c ncj(m,nj) - colour connections for the final jet nj
c----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(njmax=50000)
      dimension ep(4),ep3(4),epv(4,30,50),nqc(2),ncc(2,30,50),
     *qv(30,50),zv(30,50),qm(30,50),iqv(30,50),
     *ldau(30,49),lpar(30,50)
      common /qgarr37/ eqj(4,njmax),iqj(njmax),ncj(2,njmax),nj
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)jq,ep,iqv(1,1),nqc
      do i=1,4
       epv(i,1,1)=ep(i)
      enddo
      ncc(1,1,1)=nqc(1)
      if(iqv(1,1).eq.0)ncc(2,1,1)=nqc(2)
      nlev=1
      nrow=1

2     continue
      if(qv(nrow,nlev).eq.0.d0)then
       nj=nj+1
       if(nj.gt.njmax)stop'increase njmax!!!'
       do i=1,4
        eqj(i,nj)=epv(i,nrow,nlev)
       enddo
       iqj(nj)=iqv(nrow,nlev)
       if(iabs(iqj(nj)).eq.3)iqj(nj)=iqj(nj)*4/3

       if(iqj(nj).ne.0)then
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
       else

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
       if(debug.ge.4)write (moniou,202)
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
       if(debug.ge.4)write (moniou,203)nlev+1,ldrow,ep3

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
       if(debug.ge.4)write (moniou,203)nlev+1,ldrow+1,ep3

       if(iqv(nrow,nlev).eq.0)then
        if(iqv(ldrow,nlev+1).ne.0)then
         ncc(1,ldrow,nlev+1)=ncc(1,nrow,nlev)
         ncc(1,ldrow+1,nlev+1)=ncc(2,nrow,nlev)
        else
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
       if(debug.ge.4)write (moniou,204)
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
      
201   format(2x,'qgrec - jet reconstructing: jq=',i1
     */4x,'jet 4-momentum ep=',4(e10.3,1x)
     */4x,'jet flavor: ',i2,2x,'colour connections: ',2i3)
202   format(2x,'qgrec: ',i3,'-th final jet at level nlev=',i2,' nrow='
     *,i2/4x,'jet flavor: ',i3,2x,'jet 4-momentum:',4(e10.3,1x))
203   format(2x,'qgrec: jet at level nlev='
     *,i2,' nrow=',i2/4x,'jet 4-momentum:',4(e10.3,1x))
204   format(2x,'qgrec - end')
      end

c=============================================================================
      double precision function qgroot(qlmax,g,j)
c-----------------------------------------------------------------------------
c qgroot - effective momentum tabulation for s-channel branching
c qlmax - ln qmax/16/qtf (qmax - maximal effective momentum),
c g - dzeta number (function of random number),
c j - type of the parton (1-g,2-q)
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)qlmax,g,j
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
       ql0=ql1
       ql1=ql2
       f0=f1
       f1=f2
       goto 1
      else
       qgroot=ql2
      endif
      if(debug.ge.4)write (moniou,202)qgroot
      
201   format(2x,'qgqint - branching momentum tabulation:'
     */4x,'qlmax=',e10.3,2x,'g=',e10.3,2x,'j=',i1)
202   format(2x,'qgroot=',e10.3)
      return
      end

c=============================================================================
      subroutine qgrota(ep,s0x,c0x,s0,c0)
c-----------------------------------------------------------------------------
c qgrota - spacial rotation to the lab. system for 4-vector ep
c s0x,c0x,s0,c0 - rotation sines and cosines
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ep(4),ep1(3)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)ep,s0x,c0x,s0,c0
      ep1(3)=ep(4)
      ep1(2)=ep(2)*s0+ep(3)*c0
      ep1(1)=ep(2)*c0-ep(3)*s0
      ep(2)=ep1(1)
      ep(4)=ep1(2)*s0x+ep1(3)*c0x
      ep(3)=ep1(2)*c0x-ep1(3)*s0x
      if(debug.ge.4)write (moniou,202)ep
      
201   format(2x,'qgrota - spacial rotation:'/4x,'4-vector ep=',4(e10.3
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
c g - random number (0<g<1),
c j - type of the parton (1-g,2-q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wi(3),wk(3)
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr34/ qrt(10,101,2)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)qlmax,g,j
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
      if(debug.ge.4)write (moniou,202)qgqint
      
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
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)
      qgalf=2.d0/9.d0/dlog(qq)
      
      if(debug.ge.4)write (moniou,202)qgalf
201   format(2x,'qgalf - alpha_s/2/pi')
202   format(2x,'qgalf=',e10.3)
      return
      end

c=============================================================================
      subroutine qgtran(ep,ey,jj)
c-----------------------------------------------------------------------------
c qgtran - Lorentz boost according to parameters ey
c (along z,x,y-axis respectively - ey(1),ey(2),ey(3))
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ey(3),ep(4)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)ep,ey
      if(jj.eq.1)then
c lorentz transform to lab. system according to 1/ey(i) parameters
       do i=1,3
        if(ey(4-i).ne.1.d0)then
         wp=(ep(1)+ep(5-i))/ey(4-i)
         wm=(ep(1)-ep(5-i))*ey(4-i)
         ep(1)=.5d0*(wp+wm)
         ep(5-i)=.5d0*(wp-wm)
        endif
       enddo
      else
c lorentz transform to lab. system according to ey(i) parameters
       do i=1,3
        if(ey(i).ne.1.d0)then
         wp=(ep(1)+ep(i+1))*ey(i)
         wm=(ep(1)-ep(i+1))/ey(i)
         ep(1)=.5d0*(wp+wm)
         ep(i+1)=.5d0*(wp-wm)
        endif
       enddo
      endif
      if(debug.ge.4)write (moniou,202)ep
      
201   format(2x,'qgtran - lorentz boost for 4-vector'/4x,'ep='
     *,2x,4(e10.3,1x)/4x,'boost parameters ey=',3e10.3)
202   format(2x,'qgtran: transformed 4-vector ep=',2x,4(e10.3,1x))
      return
      end

c=============================================================================
      double precision function qgsudi(qlmax,j)
c-----------------------------------------------------------------------------
c qgsudi - timelike Sudakov formfactor interpolation
c qlmax - ln qmax/16/qtf,
c j - type of the parton (1-g,2-q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wk(3)
      common /qgarr33/ fsud(10,2)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)j,qlmax
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
      if(debug.ge.4)write (moniou,202)qgsudi
      
201   format(2x,'qgsudi - timelike form factor interpolation:'
     */4x,'parton type j=',i1,2x,'momentum logarithm qlmax=',e10.3)
202   format(2x,'qgsudi=',e10.3)
      return
      end

c=============================================================================
      double precision function qgsuds(q,j)
c-----------------------------------------------------------------------------
c qgsuds - spacelike Sudakov formfactor
c q - virtuality,
c j - type of parton (1 - g, 2 - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)j,q
      if(q.gt.1.d0)then
       qgsuds=(dlog(q)-dlog(dlog(q/alm)/dlog(1.d0/alm))
     * *(.75d0+dlog(qtf/alm)))/4.5d0
       if(j.eq.1)then
        qgsuds=qgsuds*6.d0
       else
        qgsuds=qgsuds/.375d0
       endif
       qgsuds=exp(-qgsuds)
      else
       qgsuds=1.d0
      endif     
      if(debug.ge.4)write (moniou,202)qgsuds
      
201   format(2x,'qgsuds - spacelike form factor: parton type j='
     *,i1,2x,'momentum q=',e10.3)
202   format(2x,'qgsuds=',e10.3)
      return
      end

c=============================================================================
      double precision function qgsudx(q,j)
c-----------------------------------------------------------------------------
c qgsudx - spacelike Sudakov formfactor
c q - virtuality,
c j - type of parton (1 - g, 2 - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr32/ epsxmn
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)j,q
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
      if(debug.ge.4)write (moniou,202)qgsudx
      
201   format(2x,'qgsudx - spacelike form factor: parton type j='
     *,i1,2x,'momentum q=',e10.3)
202   format(2x,'qgsudx=',e10.3)
      return
      end

c=============================================================================
      double precision function qgsudt(qmax,j)
c-----------------------------------------------------------------------------
c qgsudt - timelike Sudakov formfactor
c qmax - maximal value of the effective momentum,
c j - type of parton (1 - g, 2 - q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)j,qmax
      qgsudt=0.d0
      qlmax=dlog(dlog(qmax/16.d0/alm))
      qfl=dlog(dlog(qtf/alm))
c numerical integration over transverse momentum square;
c gaussian integration is used
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
      if(debug.ge.4)write (moniou,202)qgsudt
      
201   format(2x,'qgsudt - timelike form factor: parton type j='
     *,i1,2x,'momentum qmax=',e10.3)
202   format(2x,'qgsudt=',e10.3)
      return
      end

c=============================================================================
      double precision function qgtwd(s,a,b)
c-----------------------------------------------------------------------------
c qgtwd - kinematic function for two-particle decay 
c (light cone momentum share of the particle of mass squared a)
c b - partner's mass squared,
c s - two particle invariant mass
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)s,a,b
      x=.5d0*(1.d0+(a-b)/s)
      dx=x-dsqrt(a/s)
      if(dx.gt.0.d0)then
       x=x+dsqrt(dx)*dsqrt(x+dsqrt(a/s))
      else
       x=dsqrt(a/s)
      endif
      qgtwd=x
      if(debug.ge.4)write (moniou,202)qgtwd
      
201   format(2x,'qgtwd: s=',e10.3,2x,'a=',e10.3,2x,'b=',e10.3)
202   format(2x,'qgtwd=',e10.3)
      return
      end

c=============================================================================
      subroutine qgvdef(ich,ic1,ic2,icz)
c-----------------------------------------------------------------------------
c qgvdef - determination of valence quark flavour
c (for valence quark hard scattering)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /debug/   debug
      external psran

      if(debug.ge.3)write (moniou,201)ich,icz
      is=iabs(ich)/ich
      if(icz.eq.1)then
       ic1=ich*(1-3*int(.5d0+psran(b10)))
       ic2=-ic1-ich
      elseif(icz.eq.2)then
       if(psran(b10).gt..33333d0.or.ich.lt.0)then
        ic1=ich-is
        ic2=3*is
       else
        ic1=4*is-ich
        ic2=ich+4*is
       endif
      elseif(icz.eq.3)then
       ic1=ich-3*is
       ic2=-4*is
      elseif(icz.eq.4)then
       ic1=ich-9*is
       ic2=5*is
      endif
      if(debug.ge.4)write (moniou,202)ic1,ic2
      
201   format(2x,'qgvdef: hadron type ich=',i2,' auxilliary type icz='
     *,i1)
202   format(2x,'qgvdef-end: parton flavors ic1=',i2,
     *'ic2=',i2)
      return
      end

c=============================================================================
      double precision function qgzsim(qq,j)
c-----------------------------------------------------------------------------
c qgzsim - light cone momentum share (for timelike branching)
c qq - effective momentum value,
c j - type of the parent parton (1-g,2-q)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr11/ b10
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr43/ moniou
      common /debug/   debug
      external psran

      if(debug.ge.3)write (moniou,201)qq,j
      zmin=.5d0-dsqrt(.25d0-dsqrt(qtf/qq))
      qlf=dlog(qtf/alm)
1     continue
      if(j.eq.1)then
       qgzsim=.5d0*(2.d0*zmin)**psran(b10)
       gb=qgzsim*(qgfap(qgzsim,1,1)+qgfap(qgzsim,1,2))/7.5d0
      else
       qgzsim=zmin*((1.d0-zmin)/zmin)**psran(b10)
       gb=qgzsim*qgfap(qgzsim,2,1)*.375d0
      endif
      qt=qq*(qgzsim*(1.d0-qgzsim))**2
      gb=gb/dlog(qt/alm)*qlf
      if(psran(b10).gt.gb)goto 1
      if(debug.ge.4)write (moniou,202)qgzsim
      
201   format(2x,'qgzsim - z-share simulation: qq=',e10.3,2x,'j=',i1)
202   format(2x,'qgzsim=',e10.3)
      return
      end

c===========================================================================
      subroutine qgixxd(ich,ic1,ic2,icz)
c---------------------------------------------------------------------------
c qgixxd - parton flavours for valence quark soft interaction (charge exchange)
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr8/  wwm,be(4),dc(5),deta,almpt,ptdif,ptndi
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /debug/   debug
      external psran
        
      if(debug.ge.3)write (moniou,201)ich,icz
      is=iabs(ich)/ich
      if(icz.eq.1)then                      !pion
       ic1=ich*(1-3*int(.5d0+psran(b10)))
       if(psran(b10).lt.dc(2))then
        ic2=-4*ic1/iabs(ic1)
	if(iabs(ic1).eq.1)then
	 ich1=-5*is
	else
	 ich1=4*is
	endif
       else
        ich1=ich*int(.5d0+psran(b10))
        ic2=-ic1*iabs(ich1)-(ich+ic1)*iabs(ich-ich1)
       endif
      elseif(icz.eq.2)then
       ic1=int(1.3333d0+psran(b10))   !valence quark type (for proton)
       if(ic1.eq.1)then
        ich1=int(psran(b10)+.5d0)+2   !leading nucleon type
        ic2=1-ich1
       elseif(psran(b10).lt..5d0)then
        ich1=2
	ic2=-2
       else
        ich1=7                        !uuu
	ic2=-1
       endif
       if(iabs(ich).eq.3)then         !neutron
        ic1=3-ic1
        ic2=-3-ic2
	if(ich1.eq.7)then
	 ich1=8                       !ddd
	else
         ich1=5-ich1
	endif
       endif
       if(ich.lt.0)then
        ic1=-ic1
        ic2=-ic2
        ich1=-ich1
       endif
      elseif(icz.eq.3)then
       ic1=ich-3*is
       ic2=-is*int(1.5d0+psran(b10))
       ich1=3*is-ic2
      elseif(icz.eq.4)then
       ic1=ich-9*is
       ic2=is*int(1.5d0+psran(b10))
       ich1=9*is-ic2
      elseif(icz.eq.5)then
       ic1=is*int(1.5d0+psran(b10))
       ic2=-ic1
       ich1=ich
      endif
      ich=ich1
      if(debug.ge.4)write (moniou,202)ic1,ic2,ich
      
201   format(2x,'qgixxd: hadron type ich=',i2,' auxilliary type icz='
     *,i1)
202   format(2x,'qgixxd-end: parton flavors ic1=',i2,' ic2='
     *,i2,'new hadron type ich=',i2)
      return
      end

c=============================================================================
      subroutine qgdifr(wppr,wmtg,izp,izt,jexpr,jextg,iret)
c-----------------------------------------------------------------------------
c qgdifr - treatment of diffraction dissociation / leading hadron states
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ey(3),ep(4)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr2/  scm,wp0,wm0
      common /qgarr8/  wwm,be(4),dc(5),deta,almpt,ptdif,ptndi
      common /qgarr10/ am(7),ammu
      common /qgarr11/ b10
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr21/ dmmin(3),wex(3)
      common /qgarr43/ moniou
      common /debug/   debug
      external psran

      if(debug.ge.3)write (moniou,201)izp,izt,wppr,wmtg
      iret=0
      jexip=0
      jexit=0
      if(jexpr.eq.-2)jexip=1                             !so161205
      if(jexpr.gt.0)then                                 !so161205
       if(psran(b10).gt.(1.d0-wex(icz))**jexpr)jexip=1   !so161205
      endif                                              !so161205
      if(jextg.eq.-2.or.jextg.gt.0.and.psran(b10)
     *.gt.(1.d0-wex(2))**jextg)jexit=1
      
      if(wppr.eq.wp0.and.jexpr.gt.0.and.jexip.eq.0)jexip=1
      if(wmtg.eq.wm0.and.jextg.gt.0.and.jexit.eq.0)jexit=1

      sd0=wppr*wmtg
      if(jextg.eq.-1)then
       dmass2=0.d0
       ddmin2=0.d0
      elseif(jexit.eq.0)then
       if(iabs(izt).ge.7)then
        dmass2=dmmin(2)
       else
        dmass2=am(2)
       endif
       ddmin2=dmass2
      else
       ddmin2=dmmin(2)
      endif
      if(jexpr.eq.-1)then
       dmass1=0.d0
       ddmin1=0.d0
      elseif(jexip.eq.0)then
       if(iabs(izp).ge.7)then
        dmass1=dmmin(2)
       else
        dmass1=am(icz)
       endif
       ddmin1=dmass1
      else
       ddmin1=dmmin(icz)
       ddmax1=dsqrt(sd0)-ddmin2
      endif
      ddmax2=dsqrt(sd0)-ddmin1

      if(psran(b10).lt..5d0)then
       if(jexip.eq.1)then
        if(ddmax1.gt.ddmin1)then
         dmass1=ddmin1*(1.d0-psran(b10)*(1.d0-(ddmax1/ddmin1)
     *   **(alpd-1.d0)))**(1.d0/(alpd-1.d0))
        else
         if(iabs(izp).ge.7)then
          dmass1=dmmin(2)
         else
          dmass1=am(icz)
 	 endif
 	 jexip=0
        endif
       endif

       if(jexit.eq.1)then
        ddmax2=dsqrt(sd0)-dmass1
        if(ddmax2.gt.ddmin2)then
         dmass2=ddmin2*(1.d0-psran(b10)*(1.d0-(ddmax2/ddmin2)
     *   **(alpd-1.d0)))**(1.d0/(alpd-1.d0))
        else
         if(iabs(izt).ge.7)then
          dmass2=dmmin(2)
         else
          dmass2=am(2)
 	endif
	 jexit=0
        endif
       endif
   
      else
       if(jexit.eq.1)then
        if(ddmax2.gt.ddmin2)then
         dmass2=ddmin2*(1.d0-psran(b10)*(1.d0-(ddmax2/ddmin2)
     *   **(alpd-1.d0)))**(1.d0/(alpd-1.d0))
        else
         if(iabs(izt).ge.7)then
          dmass2=dmmin(2)
         else
          dmass2=am(2)
 	endif
	 jexit=0
        endif
       endif
   
       if(jexip.eq.1)then
        ddmax1=dsqrt(sd0)-dmass2
        if(ddmax1.gt.ddmin1)then
         dmass1=ddmin1*(1.d0-psran(b10)*(1.d0-(ddmax1/ddmin1)
     *   **(alpd-1.d0)))**(1.d0/(alpd-1.d0))
        else
         if(iabs(izp).ge.7)then
          dmass1=dmmin(2)
         else
          dmass1=am(icz)
 	 endif
 	 jexip=0
        endif
       endif
      endif
      
      wpp=wppr
      wpm=wmtg
      if(sd0.lt.(dmass1+dmass2)**2)then
       sd0=1.1d0*(dmass1+dmass2)**2
       if(wpp.lt.dsqrt(sd0/scm)*am(2))wpp=dsqrt(sd0/scm)*am(2)
       wpm=sd0/wpp
      endif
      dmass1=dmass1**2
      dmass2=dmass2**2
      
      if(jexpr.ne.-1.and.jextg.ne.-1)then
       ptmax=max(0.d0,qglam(sd0,dmass1,dmass2))
       if(jexpr.eq.-2.or.jextg.eq.-2)then
        ptmean=ptdif
       else
        ptmean=ptndi
       endif
       if(ptmax.lt.ptmean**2)then
1       pti=ptmax*psran(b10)
        if(psran(b10).gt.exp(-dsqrt(pti)/ptmean))goto 1
       else
2       pti=(ptmean*dlog(psran(b10)*psran(b10)))**2
        if(pti.gt.ptmax)goto 2
       endif
      else
       pti=0.d0
      endif
      amt1=dmass1+pti
      amt2=dmass2+pti
      wpd1=wpp*qgtwd(sd0,amt1,amt2)
      if(wpd1.gt.0.d0)then
       wmd1=amt1/wpd1
      else
       wmd1=0.d0
      endif
      wmd2=wpm-wmd1
      if(wmd2.gt.0.d0)then
       wpd2=amt2/wmd2
      else
       wpd2=0.d0
      endif
      pt=dsqrt(pti)
      call qgcs(c,s)
      
      if(jexpr.eq.-1)then
       wppr=wpd1
       if(wmd1.ne.0.d0)stop'wmd1.ne.0!!!'
      else
       ep(1)=.5d0*(wpd1+wmd1)
       ep(2)=.5d0*(wpd1-wmd1)
       ep(3)=pt*c
       ep(4)=pt*s
       wppr=0.d0
       if(jexip.eq.0)then
        call qgreg(ep,izp)
       else
        if(izp.ne.0)is=iabs(izp)/izp     
        if(icz.eq.1)then
	 if(iabs(izp).ge.4)then
	  ic2=-4*is
	  ic1=izp-3*is
	  if(psran(b10).lt..5d0)then
	   ics=ic1
	   ic1=ic2
	   ic2=ics
	  endif
         elseif(izp.ne.0)then
          ic1=izp*(1-3*int(.5d0+psran(b10)))
          ic2=-izp-ic1
         else
          ic1=int(1.5d0+psran(b10))*(2*int(.5d0+psran(b10))-1)
          ic2=-ic1
         endif
        elseif(icz.eq.2)then
         if(iabs(izp).ge.7)then
	  ic1=izp-is
	  ic2=izp-6*is
	 else
	  if(psran(b10).gt..33333d0)then
           ic1=3*is
           ic2=izp-is
          else
           ic1=izp+4*is
           ic2=4*is-izp
	  endif
	 endif
        elseif(icz.eq.3)then
         ic1=-4*is
         ic2=izp-3*is
        endif      
        call qgdeft(dmass1,ep,ey)
        call qggene(dsqrt(dmass1),dsqrt(dmass1),ey
     *  ,0.d0,1.d0,0.d0,1.d0,ic1,ic2)
       endif
      endif

      if(jextg.eq.-1)then
       wmtg=wmd2
       if(wpd2.ne.0.d0)stop'wpd2.ne.0!!!'
      else
       ep(1)=.5d0*(wpd2+wmd2)
       ep(2)=.5d0*(wpd2-wmd2)
       ep(3)=-pt*c
       ep(4)=-pt*s
       wmtg=0.d0
       if(jexit.eq.0)then
        call qgreg(ep,izt)
       else
        is=iabs(izt)/izt
        if(iabs(izt).ge.7)then
	 ic1=izt-is
	 ic2=izt-6*is
	else
	 if(psran(b10).gt..33333d0)then
          ic1=3*is
          ic2=izt-is
         else
          ic1=izt+4*is
          ic2=4*is-izt
	 endif
	endif
        call qgdeft(dmass2,ep,ey)
        call qggene(dsqrt(dmass2),dsqrt(dmass2),ey
     *  ,0.d0,1.d0,0.d0,1.d0,ic2,ic1)
       endif
      endif    
      if(debug.ge.4)write (moniou,202)
      
201   format(2x,'qgdifr - leading clusters hadronization:'
     */4x,'cluster types izp=',i2,2x,
     *'izt=',i2/4x,'available light cone momenta: wppr=',e10.3,
     *' wmtg=',e10.3)
202   format(2x,'qgdifr - end')
      return
      end

c=============================================================================
      subroutine qgfau(b,gz)
c-----------------------------------------------------------------------------
c qgfau - integrands for hadron-hadron and hadron-nucleus cross-sections calc.
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207)
      dimension gz(3),gz0(5)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)b
      do l=1,3
       gz(l)=0.d0
      enddo

      ab=float(ia(2))
      do iddp1=1,2
      do iddp2=1,2
       call qgfz(b,gz0,iddp1,iddp2)
       if(iddp1.eq.iddp2)gz(1)=gz(1)+(1.d0-gz0(1)*anorm)**ab
     * *cc(iddp1,icz)
       do l=2,3
        gz(l)=gz(l)+(1.d0-gz0(l-1)*anorm)**ab
     *  *cc(iddp1,icz)*cc(iddp2,icz)
       enddo
      enddo
      enddo
                                                                       
      gz(3)=gz(2)-gz(3)
      gz(2)=gz(1)-gz(2)
      gz(1)=1.d0-gz(1)
      if(debug.ge.4)write (moniou,202)gz
      if(debug.ge.4)write (moniou,203)
      
201   format(2x,'qgfau - integrands for hadron-hadron and hadron'
     *,'-nucleus cross-sections'/4x,'b=',e10.3)
202   format(2x,'qgfau: gz=',3e10.3)
203   format(2x,'qgfau - end')
      return
      end

c=============================================================================
      subroutine qgfrag(sa,na,rc)
c-----------------------------------------------------------------------------
c qgfrag - multifragmentation (search for connected nucleon clasters)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207)
      dimension sa(iapmax,3)
      common /qgarr13/ nsf,iaf(iapmax)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)na
      if(debug.ge.4)then
       write (moniou,202)
       do i=1,na
        write (moniou,203)(sa(i,l),l=1,3)
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
      if(debug.ge.4)write (moniou,204)nsf,iaf(nsf)

      ng=1
      j=ni
      ni=ni+1
      if(na-ni)6,5,1
5     nsf=nsf+1
      iaf(nsf)=1
      if(debug.ge.4)write (moniou,204)nsf,iaf(nsf)
6     continue
      if(debug.ge.4)write (moniou,205)
      
201   format(2x,'qgfrag-multifragmentation: nucleus mass number: na='
     *,i2)
202   format(2x,'nucleons coordinates:')
203   format(2x,3e10.3)
204   format(2x,'qgfrag: fragment n',i2,2x,'fragment mass - ',i2)
205   format(2x,'qgfrag - end')
      return
      end

c=============================================================================
      subroutine qgfrgm(ns,xa)
c-----------------------------------------------------------------------------
c qgfrgm - fragmentation of nuclear spectator part
c xa - spectator nucleon coordinates,
c ns - total number of spectators
c to determine:
c nsf    - number of secondary fragments;
c iaf(i) - mass of the i-th fragment
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter(iapmax=207)
      dimension xa(iapmax,3)
      integer debug
      common /qgarr1/  ia(2),icz,icp
      common /qgarr3/  rmin,emax,eev
      common /qgarr11/ b10
      common /qgarr13/ nsf,iaf(iapmax)
      common /qgarr43/ moniou
      common /debug/   debug
      external psran

      if(debug.ge.3)write (moniou,201)ns

      nsf=0
      if(ns-1)6,1,2
1     nsf=nsf+1                               !single spectator nucleon recorded
      iaf(nsf)=1
      if(debug.ge.4)write (moniou,202)
      goto 6

c excitation energy is calculated summing contributions of all wounded nucleons
c partial excitation is simulated as f(e) ~ 1/sqrt(e)* exp(-e/(2*<e>))
2     eex=0.d0         
      do i=1,ia(1)-ns
       eex=eex+(psran(b10)+psran(b10)+psran(b10)+
     * psran(b10)+psran(b10)-2.5d0)**2*2.4d0
      enddo     
      if(debug.ge.4)write (moniou,203)eex

c multifragmentation if excitation energy per spectator is larger than emax
      if(eex/ns.gt.emax)then
       call qgfrag(xa,ns,rmin)       !multifragmentation - percolation algorithm
      else
c otherwise average number of eveporated nucleons equals eex/eev
c (eev - mean excitation energy carried out by one nucleon)
       nf=npgen(eex/eev,0,ns-1)
       nsf=nsf+1
       iaf(nsf)=ns-nf                !recording the produced fragment
       if(debug.ge.4)write (moniou,204)iaf(nsf)

c some part of excitation energy is carried out by alphas; we determine the
c number of alphas simply as nf/4
       nal=nf/4
       if(nal.ne.0)then
        do i=1,nal                   !recording of the evaporated alphas
         nsf=nsf+1
         iaf(nsf)=4
        enddo
       endif
       nf=nf-4*nal
       if(nf.ne.0)then
        do i=1,nf
         nsf=nsf+1
         iaf(nsf)=1                   !recording of the evaporated nucleons
        enddo
       endif
       if(debug.ge.4)write (moniou,205)nf,nal
      endif
6     continue
      if(debug.ge.4)write (moniou,206)
      
201   format(2x,'qgfrgm: number of spectators: ns=',i2)
202   format(2x,'qgfrgm - single spectator')
203   format(2x,'qgfrgm: excitation energy: eex=',e10.3)
204   format(2x,'qgfrgm - evaporation: mass number of the fragment:',i2)
205   format(2x,'qgfrgm - evaporation: number of nucleons nf='
     *,i2,'number of alphas nal=',i2)
206   format(2x,'qgfrgm - end')
      return
      end

c=============================================================================
      subroutine qggau(gz)
c-----------------------------------------------------------------------------
c qggau - impact parameter integration for impact parameters < bm 
c (for hadron-hadron and hadron-nucleus cross-sections)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension gz(3),gz0(3)
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)
      do i=1,3
       gz(i)=0.d0
      enddo
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
      if(debug.ge.4)write (moniou,202)
      
201   format(2x,'qggau - nuclear cross-sections calculation')
202   format(2x,'qggau - end')
      return
      end

c=============================================================================
      subroutine qggau1(gz)
c-----------------------------------------------------------------------------
c qggau1 - impact parameter integration for impact parameters > bm
c (for hadron-hadron and hadron-nucleus cross-sections)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension gz(3),gz0(3)
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)
      do i=1,7
      do m=1,2
       b=bm-wsnuc(2)*dlog(.5d0+x1(i)*(m-1.5d0))
       call qgfau(b,gz0)
       do l=1,3
        gz(l)=gz(l)+gz0(l)*a1(i)*exp((b-bm)/wsnuc(2))*b*pi*wsnuc(2)
       enddo
      enddo
      enddo     
      if(debug.ge.4)write (moniou,202)
      
201   format(2x,'qggau1 - nuclear cross-sections calculation')
202   format(2x,'qggau1 - end')
      return
      end

c=============================================================================
      double precision function qganrm(rnuc,wsnuc,wbnuc)
c-----------------------------------------------------------------------------
c qganrm - normalization for 3-parameter Fermi nuclear density
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr6/  pi,bm,amws
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)
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
      if(debug.ge.4)write (moniou,202)qganrm
      
201   format(2x,'qganrm - nuclear density normalization')
202   format(2x,'qganrm=',e10.3)
      return
      end

c=============================================================================
      subroutine qggene(wp0,wm0,ey0,s0x,c0x,s0,c0,ic1,ic2)
c-----------------------------------------------------------------------------
c qggene - string fragmentation into secondary hadrons
c wp0, wm0 - initial LC momenta ( e+p, e-p ) of the (di-)quarks - string ends,
c ic1, ic2 - their types:
c 1 - u, -1 - U, 2 - d, -2 - D, 3 - ud, -3 - UD, 4 - s, -4 - S,
c 5 - c, -5 - C, 6 - uu, 7 - dd, -6 - UU, -7 - DD
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      character *2 tyq
cdh   dimension wp(2),ic(2),ept(4),ep(4),ep1(4),ey(3),ey0(3)
      dimension wp(2),ic(2),ept(4),ep(4),       ey(3),ey0(3)
      common /qgarr8/  wwm,bep,ben,bek,bec,dc(5),deta,almpt,ptdif,ptndi
      common /qgarr10/ am0,amn,amk,amc,amlamc,amlam,ameta,ammu
      common /qgarr11/ b10
      common /qgarr19/ ahl(3)
      common /qgarr28/ arr(5)
      common /qgarr42/ tyq(16)
      common /qgarr43/ moniou
      common /debug/   debug
      external psran

      if(debug.ge.3)write (moniou,201)tyq(8+ic1),tyq(8+ic2)
     *,wp0,wm0,ey0,s0x,c0x,s0,c0     
      ww=wp0*wm0
      ept(1)=.5d0*(wp0+wm0)
      ept(2)=.5d0*(wp0-wm0)
      ept(3)=0.d0
      ept(4)=0.d0
      ic(1)=ic1
      ic(2)=ic2

1     sww=dsqrt(ww)
      call qgdeft(ww,ept,ey)
      j=int(2.d0*psran(b10))+1

      if(debug.ge.5)then
       iqt=8+ic(j)
       write (moniou,202)j,tyq(iqt),ww
      endif

      iab=iabs(ic(j))
      is=ic(j)/iab
      if(iab.gt.5)iab=3
      iaj=iabs(ic(3-j))
      if(iaj.gt.5)iaj=3
      if(iaj.eq.3)then
       restm=amn
      elseif(iaj.eq.4)then
       restm=amk
      elseif(iaj.eq.5)then
       restm=amc
      else
       restm=am0
      endif

      if(iab.le.2.and.sww.gt.restm+2.d0*am0+wwm.or.
     *iab.eq.3.and.sww.gt.restm+am0+amn+wwm.or.
     *iab.eq.4.and.sww.gt.restm+am0+amk+wwm.or.
     *iab.eq.5.and.sww.gt.restm+am0+amc+wwm)then
       if(iab.le.2)then
        if(sww.gt.restm+2.d0*amn.and.psran(b10).lt.dc(1))then
c nucleon generation
         restm=(restm+amn)**2
         bet=ben
         ami=amn**2
         alf=almpt-arr(2)
         blf=1.D0-ARR(1)-ARR(2)
         ic0=ic(j)+is
         ic(j)=-3*is
        elseif(sww.gt.restm+2.d0*amk.and.psran(b10).lt.dc(2))then
c kaon generation
         restm=(restm+amk)**2
         bet=bek
         ami=amk**2
         alf=almpt-arr(3)
         blf=1.D0-ARR(1)-ARR(3)
         ic0=ic(j)+3*is
         ic(j)=4*is
        elseif(sww.gt.restm+ameta+am0.and.psran(b10).lt.deta)then
c eta generation
         restm=(restm+am0)**2
         bet=bek
         ami=ameta**2
         alf=almpt-arr(1)
         blf=1.D0-2.D0*ARR(1)
         ic0=10
        else
c pion generation
         restm=(restm+am0)**2
         bet=bep
         ami=am0**2
         alf=almpt-arr(1)
         blf=1.D0-2.D0*ARR(1)
         if(psran(b10).lt..3333d0)then
          ic0=0
         else
          ic0=3*is-2*ic(j)
          ic(j)=3*is-ic(j)
         endif
        endif

       elseif(iab.eq.3)then
        if(sww.gt.restm+amk+amlam.and.psran(b10).lt.dc(4)
     *  .and.iabs(ic(j)).eq.3)then
c lambda generation
         restm=(restm+amk)**2
         bet=bek
         ami=amlam**2
         alf=almpt-arr(3)
         blf=1.D0-ARR(1)-ARR(2)+arr(1)-arr(3)
         ic0=6*is
         ic(j)=-4*is
        else
c nucleon generation
         restm=(restm+am0)**2
         bet=ben
         ami=amn**2
         alf=almpt-arr(1)
         blf=1.D0-ARR(1)-ARR(2)
         if(iabs(ic(j)).eq.3)then
          ic0=is*int(2.5d0+psran(b10))
          ic(j)=is-ic0
         else
          ic0=ic(j)-4*is
          ic(j)=ic0-4*is
         endif
        endif

       elseif(iab.eq.4)then
        if(sww.gt.restm+amn+amlam.and.psran(b10).lt.dc(1))then
c lambda generation
         restm=(restm+amn)**2
         bet=ben
         ami=amlam**2
         alf=almpt-arr(2)
         blf=1.D0-ARR(1)-ARR(2)+arr(1)-arr(3)
         ic0=6*is
         ic(j)=-3*is
        else
c kaon generation
         restm=(restm+am0)**2
         bet=bep
         ami=amk**2
         alf=almpt-arr(1)
         blf=1.D0-ARR(1)-ARR(3)
         ic(j)=is*int(1.5d0+psran(b10))
         ic0=-3*is-ic(j)
        endif
       endif

       ptmax=qglam(ww,restm,ami)
       if(ptmax.lt.0.)ptmax=0.
       if(ptmax.lt.bet**2)then
2       pti=ptmax*psran(b10)
        if(psran(b10).gt.exp(-dsqrt(pti)/bet))goto 2
       else
3       pti=(bet*dlog(psran(b10)*psran(b10)))**2
        if(pti.gt.ptmax)goto 3
       endif

       amt=ami+pti
       restm1=restm+pti
       zmin=dsqrt(amt/ww)
       zmax=qgtwd(ww,amt,restm1)
       z1=(1.-zmax)**alf
       z2=(1.-zmin)**alf
4      z=1.-(z1+(z2-z1)*psran(b10))**(1./alf)
       if(psran(b10).gt.(z/zmax)**blf)goto 4
       wp(j)=z*sww
       wp(3-j)=amt/wp(j)
       ep(1)=.5d0*(wp(1)+wp(2))
       ep(2)=.5d0*(wp(1)-wp(2))
       pti=dsqrt(pti)
       call qgcs(c,s)
       ep(3)=pti*c
       ep(4)=pti*s
       ept(1)=sww-ep(1)
       do i=2,4
        ept(i)=-ep(i)
       enddo
       ww=qgnrm(ept)
       if(ww.lt.restm)goto 4

       call qgtran(ep,ey,1)
       call qgtran(ept,ey,1)
       if(s0x.ne.0.d0.or.s0.ne.0.d0)then
        call qgrota(ep,s0x,c0x,s0,c0)
       endif
       if(ey0(1)*ey0(2)*ey0(3).ne.1.d0)then
        call qgtran(ep,ey0,1)
       endif
       call qgreg(ep,ic0)
       goto 1
      endif

      ami2=restm**2
      bet=bep
      if(iab.le.2.and.iaj.le.2)then
       if(sww.gt.2.d0*amk.and.psran(b10).lt.dc(2))then
        bet=bek
        ami=amk**2
        ami2=ami
        ic(j)=ic(j)+3*is
        ic(3-j)=ic(3-j)-3*is
       else
        ami=am0**2
        ic0=-ic(1)-ic(2)
        if(ic0.ne.0)then
         ic(j)=ic0*int(.5d0+psran(b10))
         ic(3-j)=ic0-ic(j)
        else
         if(psran(b10).lt..2d0)then
          ic(j)=0
          ic(3-j)=0
         else
          ic(j)=3*is-2*ic(j)
          ic(3-j)=-ic(j)
         endif
        endif
       endif
	
      elseif(iab.eq.3.or.iaj.eq.3)then
       if(iab.eq.3)then
        ami=amn**2
        if(iabs(ic(j)).eq.3)then
         if(iaj.eq.3)then
          if(iabs(ic(3-j)).eq.3)then
           if(sww.gt.2.d0*amlam.and.psran(b10).lt.dc(4))then
            bet=bek
            ami=amlam**2
            ami2=ami
            ic(j)=6*is
            ic(3-j)=-6*is
	   else
            ic(j)=is*int(2.5d0+psran(b10))
            ic(3-j)=-ic(j)
	   endif
          else
           ic(3-j)=ic(3-j)+4*is
           ic(j)=5*is+ic(3-j)
          endif
         elseif(iaj.lt.3)then
          if(sww.gt.amlam+amk.and.psran(b10).lt.dc(4))then
           bet=bek
           ami=amlam**2
           ami2=amk**2
           ic(j)=6*is
           ic(3-j)=ic(3-j)+3*is
	  else
           if(psran(b10).lt..3333d0)then
            ic(j)=ic(3-j)+is
            ic(3-j)=0
           else
            ic(j)=is*(4-iaj)
            ic(3-j)=is*(3-2*iaj)
           endif
          endif
         elseif(iaj.eq.4)then
          ic(j)=is*int(2.5d0+psran(b10))
          ic(3-j)=-ic(j)-2*is
         elseif(iaj.eq.5)then
          ic(j)=is*int(2.5d0+psran(b10))
          ic(3-j)=-ic(j)+10*is
         endif
        else
         ic(j)=ic(j)-4*is
         ic0=ic(j)-4*is
         if(iaj.eq.3)then
          ic(3-j)=ic0-is
         elseif(iaj.lt.3)then
          ic(3-j)=-ic(3-j)-ic0
         elseif(iaj.eq.4)then
          ic(3-j)=ic0-3*is
         elseif(iaj.eq.5)then
          ic(3-j)=ic0+9*is
         endif
        endif
       else
        if(iabs(ic(3-j)).eq.3)then
         if(iab.lt.3)then
          if(sww.gt.amlam+amk.and.psran(b10).lt.dc(4))then
           bet=bek
           ami2=amlam**2
           ami=amk**2
           ic(j)=ic(j)+3*is
           ic(3-j)=6*is
	  else
           ami=am0**2
           if(psran(b10).lt..3333d0)then
            ic(3-j)=ic(j)+is
            ic(j)=0
           else
            ic(3-j)=is*(4-iab)
            ic(j)=is*(3-2*iab)
           endif
          endif
	 elseif(iab.eq.4)then
          ami=amk**2
          ic(3-j)=is*int(2.5d0+psran(b10))
          ic(j)=-ic(3-j)-2*is
         elseif(iab.eq.5)then
          ami=amc**2
          ic(3-j)=is*int(2.5d0+psran(b10))
          ic(j)=-ic(3-j)+10*is
         endif
        else
         ic(3-j)=ic(3-j)-4*is
         ic0=ic(3-j)-4*is
         if(iab.lt.3)then
          ami=am0**2
          ic(j)=-ic0-ic(j)
         elseif(iab.eq.4)then
          ami=amk**2
          ic(j)=ic0-3*is
         elseif(iab.eq.5)then
          ami=amc**2
          ic(j)=ic0+9*is
         endif
        endif
       endif
      elseif(iab.eq.4.or.iaj.eq.4)then
       if(iab.eq.4)then
        ami=amk**2
        if(iaj.eq.4)then
         ic(j)=-is*int(4.5d0+psran(b10))
         ic(3-j)=-ic(j)
        elseif(iaj.eq.5)then
         ic(j)=-is*int(4.5d0+psran(b10))
         ic(3-j)=-ic(j)-12*is
        else
         ic0=ic(3-j)+int(.6667d0+psran(b10))*(-3*is-2*ic(3-j))
         ic(j)=ic0-3*is
         ic(3-j)=ic0-ic(3-j)
        endif
       else
        if(iab.le.2)then
         ami=am0**2
         ic0=ic(j)+int(.6667d0+psran(b10))*(3*is-2*ic(j))
         ic(j)=ic0-ic(j)
         ic(3-j)=ic0+3*is
        elseif(iab.eq.5)then
         ami=amc**2
         ic(3-j)=is*int(4.5d0+psran(b10))
         ic(j)=-ic(3-j)+12*is
        endif
       endif
      elseif(iab.eq.5.or.iaj.eq.5)then
       if(iab.eq.5)then
        ami=amc**2
        if(iaj.eq.5)then
         ic(j)=is*int(7.5d0+psran(b10))
         ic(3-j)=-ic(j)
        else
         ic0=ic(3-j)+int(.6667d0+psran(b10))*(-3*is-2*ic(3-j))
         ic(j)=ic0+9*is
         ic(3-j)=ic0-ic(3-j)
        endif
       else
        ami=am0**2
        ic0=ic(j)+int(.6667d0+psran(b10))*(3*is-2*ic(j))
        ic(j)=ic0-ic(j)
        ic(3-j)=ic0-9*is
       endif
      endif

      ptmax=qglam(ww,ami2,ami)
      if(ptmax.lt.0.)ptmax=0.
      if(ptmax.lt.bet**2)then
5      pti=ptmax*psran(b10)
       if(psran(b10).gt.exp(-dsqrt(pti)/bet))goto 5
      else
6      pti=(bet*dlog(psran(b10)*psran(b10)))**2
       if(pti.gt.ptmax)goto 6
      endif
      amt1=ami+pti
      amt2=ami2+pti
      z=qgtwd(ww,amt1,amt2)
      wp(j)=z*sww
      wp(3-j)=amt1/wp(j)
      ep(1)=.5d0*(wp(1)+wp(2))
      ep(2)=.5d0*(wp(1)-wp(2))
      pti=dsqrt(pti)
      call qgcs(c,s)
      ep(3)=pti*c
      ep(4)=pti*s
      ept(1)=sww-ep(1)
      do i=2,4
       ept(i)=-ep(i)
      enddo
      call qgtran(ep,ey,1)
      call qgtran(ept,ey,1)
      if(s0x.ne.0.d0.or.s0.ne.0.d0)then
        call qgrota(ep,s0x,c0x,s0,c0)
        call qgrota(ept,s0x,c0x,s0,c0)
      endif
      if(ey0(1)*ey0(2)*ey0(3).ne.1.d0)then
        call qgtran(ep,ey0,1)
        call qgtran(ept,ey0,1)
      endif	
      call qgreg(ep,ic(j))
      call qgreg(ept,ic(3-j))
      
      if(debug.ge.4)write (moniou,203)
      return
      
201   format(2x,'qggene: parton flavors at the ends of the string:'
     *,2x,a2,2x,a2/4x,'light cone momenta of the string: ',e10.3
     *,2x,e10.3/4x,'ey0=',3e10.3/4x,'s0x=',e10.3,2x,'c0x=',e10.3
     *,2x,'s0=',e10.3,2x,'c0=',e10.3)
202   format(2x,'qggene: current parton flavor at the end '
     *,i1,' of the string: ',a2/4x,' string mass: ',e10.3)
203   format(2x,'qggene - end')
      end

c=============================================================================
      subroutine qgxjet
c-----------------------------------------------------------------------------
c qgxjet - jet hadronization
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(njmax=50000)
      dimension ep(4),ept(4),ept1(4),ey(3)
     *,epj(4,2,2*njmax),ipj(2,2*njmax)
      common /qgarr8/  wwm,bep,ben,bek,bec,dc(5),deta,almpt
     *,ptdif,ptndi
      common /qgarr10/ am(7),ammu
      common /qgarr11/ b10
      common /qgarr36/ epjet(4,njmax),ipjet(njmax),njtot
      common /qgarr43/ moniou
      common /debug/   debug
      external psran

      if(debug.ge.3)write (moniou,201)njtot
      nj0=1
      njet0=0
      nrej=0

1     njet=njet0
      do i=1,4
       ept(i)=epjet(i,nj0)
       epj(i,1,njet+1)=ept(i)
      enddo
      iq1=ipjet(nj0)
      ipj(1,njet+1)=iq1
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
       aks=psran(b10)
       do i=1,4
        epi=epjet(i,ij)*aks
        epj(i,2,njet)=epi
        ept(i)=ept(i)+epi
       enddo
       if(psran(b10).lt.dc(2))then
        ipj(2,njet)=4*(2*jq-3)
	am2=am(3)
       else
        ipj(2,njet)=int(1.5d0+psran(b10))*(2*jq-3)
	am2=am(1)
       endif

       if(qgnrm(ept).gt.(am1+am2)**2)then
        if(debug.ge.5)write (moniou,202)njet,ipj(1,njet),ipj(2,njet)
     *  ,qgnrm(ept),ept
        ipj(1,njet+1)=-ipj(2,njet)
        do i=1,4
         ept(i)=epjet(i,ij)-epj(i,2,njet)
         epj(i,1,njet+1)=ept(i)
        enddo
        am1=am2
        goto 2
       elseif(nrej.lt.100000)then
        nrej=nrej+1
        goto 1
       else
3       continue
        do i=1,4
         ept(i)=epjet(i,ij)+epjet(i,ij-1)+epjet(i,ij+1)
         ep(i)=epjet(i,ij-1)
         ept1(i)=ept(i)
        enddo
        ww=qgnrm(ept1)
        if(ww.le.0.)then
         if(ij.gt.nj0+1)then
          ij=ij-1
          goto 3
         else
          ij=ij+1
          goto 3
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
        if(debug.ge.5)write (moniou,202)njet,ipj(1,njet),ipj(2,njet)
     *  ,qgnrm(ept),ept

        nj0=ij+1
        njet0=njet
        nrej=0
        if(ij.lt.njtot)then
         goto 1
        else
         goto 5
        endif
       elseif(nrej.lt.100000)then
        nrej=nrej+1
        goto 1
       else
4       continue
        do i=1,4
         ept(i)=epjet(i,ij)+epjet(i,ij-1)+epjet(i,ij-2)
         ep(i)=epjet(i,ij-2)
         ept1(i)=ept(i)
        enddo
        ww=qgnrm(ept1)
        if(ww.lt.0.d0)then
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

5     continue
      do ij=1,njet
       do i=1,4
        ep(i)=epj(i,1,ij)
        ept(i)=ep(i)+epj(i,2,ij)
       enddo
       ww=qgnrm(ept)                     !invariant mass squared for the jet
       if(debug.ge.4)write (moniou,203)ij,njet,ww,ipj(1,ij),ipj(2,ij)
       sww=dsqrt(ww)
       call qgdeft(ww,ept,ey)
       call qgtran(ep,ey,-1)
       call qgdefr(ep,s0x,c0x,s0,c0)
       call qggene(sww,sww,ey,s0x,c0x,s0,c0,ipj(1,ij),ipj(2,ij))
      enddo      
      if(debug.ge.4)write (moniou,204) 
           
201   format(2x,'qgxjet - total number of jets njtot=',i4)
202   format(2x,'qgxjet: njet=',i3,2x,'ic=',2i2,2x,'mass=',e10.3
     *,2x,'ep=',4e10.3)
203   format(2x,'qgxjet: ij=',i2,2x,'njet=',i3,2x,'ww=',e10.3
     *,2x,'ic=',2i3)
204   format(2x,'qgxjet - end')
      return
      end

c=============================================================================
      subroutine qgreg(ep0,ic)
c-----------------------------------------------------------------------
c qgreg - registration of produced hadron
c ep0 - 4-momentum,
c ic  - hadron type
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(nptmax=95000)
      dimension ep(4),ep0(4),ep1(4),ep2(4)
      common /qgarr4/  ey0(3)
      common /qgarr10/ am0,amn,amk,amc,amlamc,amlam,ameta,ammu
      common /qgarr11/ b10
      common /qgarr12/ nsh
      common /qgarr14/ esp(4,nptmax),ich(nptmax)
      common /qgarr21/ dmmin(3),wex(3)
      common /qgarr43/ moniou
      common /debug/   debug
      external psran
      
      if(debug.ge.3)write (moniou,201)ic,ep0,nsh
      nsh=nsh+1
      if(nsh.gt.nptmax)stop'increase nptmax!!!'
      iab=iabs(ic)
      do i=1,4
       ep(i)=ep0(i)
      enddo
      call qgtran(ep,ey0,1)

      if(iab.eq.7.or.iab.eq.8)then
       call qgdec2(ep,ep1,ep2,dmmin(2)**2,amn**2,am0**2)
       ich(nsh)=ic-5*ic/iab
       do i=1,4
        esp(i,nsh)=ep1(i)
        ep(i)=ep2(i)
       enddo
       nsh=nsh+1
       ich(nsh)=15*ic/iab-2*ic
      elseif(iab.eq.5)then
        ich(nsh)=10*int(.5d0+psran(b10))-5
      elseif(ic.eq.16)then
        ich(nsh)=12
      else
        ich(nsh)=ic
      endif
      do i=1,4
       esp(i,nsh)=ep(i)
      enddo
      if(debug.ge.4)write (moniou,202)
      
201   format(2x,'qgreg: ic=',i2,2x,'c.m. 4-momentum:',2x,4(e10.3,1x)/
     * 4x,'number of particles in the storage: ',i5)
202   format(2x,'qgreg - end')
      return
      end

c=============================================================================
      double precision function qgrot(b,s)
c-----------------------------------------------------------------------------
c qgrot - convolution of nuclear profile functions (axial angle integration)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr52/ x2(4),a2
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)b,s
      qgrot=0.d0
      do i=1,4
       sb1=b**2+s**2-2.*b*s*(2.*x2(i)-1.)
       sb2=b**2+s**2-2.*b*s*(1.-2.*x2(i))
       qgrot=qgrot+(qgt(sb1)+qgt(sb2))
      enddo
      qgrot=qgrot*a2
      if(debug.ge.4)write (moniou,202)qgrot
      
201   format(2x,'qgrot - axial angle integration of the ',
     *'nuclear profile function'/4x,
     *'impact parameter b=',e10.3,2x,'nucleon coordinate s=',e10.3)
202   format(2x,'qgrot=',e10.3)
      return
      end

c=============================================================================
      subroutine qgstr(wpi0,wmi0,wp0,wm0,ic10,ic120,ic210,ic20,jp,jt)
c-----------------------------------------------------------------------------
c qgstr - fragmentation process for the Pomeron 
c (quarks and antiquarks types at the ends of the two strings,
c energy-momentum sharing, string fragmentation)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ey(3)
      common /qgarr6/  pi,bm,amws
      common /qgarr8/  wwm,be(4),dc(5),deta,almpt,ptdif,ptndi
      common /qgarr10/ am(7),ammu
      common /qgarr11/ b10
      common /qgarr17/ dels,alfp,sigs,rr,r3p,g3p,alpd,delh
      common /qgarr43/ moniou
      common /debug/   debug
      external psran

      if(debug.ge.3)write (moniou,201)wpi0,wmi0,wp0,wm0
      do i=1,3
       ey(i)=1.d0
      enddo
      wpi=wpi0
      wmi=wmi0
c quark-antiquark types (1 - u, 2 - d, -1 - u~, -2 - d~); s-quarks are
c taken into consideration at the fragmentation step
      if(ic10.eq.0)then
       if(psran(b10).lt.dc(2))then
        ic1=4
	ic12=-4
       else
        ic1=int(1.5+psran(b10))
        ic12=-ic1
       endif
      elseif(ic10.gt.0)then
       ic1=ic10
       ic12=ic120
      else
       ic1=ic120
       ic12=ic10
      endif

      if(ic20.eq.0)then
       if(psran(b10).lt.dc(2))then
        ic2=4
	ic21=-4
       else
        ic2=int(1.5+psran(b10))
        ic21=-ic2
       endif
      elseif(ic20.gt.0)then
       ic2=ic20
       ic21=ic210
      else
       ic2=ic210
       ic21=ic20
      endif

c longitudinal momenta for the strings
      if(jp.eq.0)then
       wp1=wpi*cos(pi*psran(b10))**2
      else
1      xp=.5d0*psran(b10)**(2.d0/(1.d0-alpd))
       if(psran(b10).gt.(2.d0*(1.d0-xp))**(-.5d0*(1.d0+alpd)))goto 1
       wp1=wpi*xp
       if(psran(b10).lt..5d0)wp1=wpi-wp1
      endif
      if(jt.eq.0)then
       wm1=wmi*cos(pi*psran(b10))**2
      else
2      xm=.5d0*psran(b10)**(2.d0/(1.d0-alpd))
       if(psran(b10).gt.(2.d0*(1.d0-xm))**(-.5d0*(1.d0+alpd)))goto 2
       wm1=wmi*xm
       if(psran(b10).lt..5d0)wm1=wmi-wm1
      endif
      wpi=wpi-wp1
      wmi=wmi-wm1
c string masses
      sm1=wp1*wm1
      sm2=wpi*wmi

c mass thresholds
      if(iabs(ic1).le.2)then
       am1=am(1)
      elseif(iabs(ic1).eq.3)then
       am1=am(2)
      elseif(iabs(ic1).eq.4)then
       am1=am(3)
      endif
      if(iabs(ic2).le.2)then
       am2=am(1)
      elseif(iabs(ic2).eq.3)then
       am2=am(2)
      elseif(iabs(ic2).eq.4)then
       am2=am(3)
      endif
      if(iabs(ic12).le.2)then
       am12=am(1)
      elseif(iabs(ic12).eq.3)then
       am12=am(2)
      elseif(iabs(ic12).eq.4)then
       am12=am(3)
      endif
      if(iabs(ic21).le.2)then
       am21=am(1)
      elseif(iabs(ic21).eq.3)then
       am21=am(2)
      elseif(iabs(ic21).eq.4)then
       am21=am(3)
      endif

c too short strings are neglected (energy is given to partner string 
c or to the hadron (nucleon) to which the Pomeron is connected)
      if(sm1.gt.am1+am21.and.sm2.gt.am2+am12)then
c strings fragmentation is simulated - gener
       call qggene(wp1,wm1,ey,0.d0,1.d0,0.d0,1.d0,ic1,ic21)
       call qggene(wpi,wmi,ey,0.d0,1.d0,0.d0,1.d0,ic12,ic2)
      elseif((wpi+wp1)*(wmi+wm1).gt.am1+am21)then
       call qggene(wp1+wpi,wm1+wmi,ey,0.d0,1.d0,0.d0,1.d0,ic1,ic21)
      elseif((wpi+wp1)*(wmi+wm1).gt.am2+am12)then
       call qggene(wp1+wpi,wm1+wmi,ey,0.d0,1.d0,0.d0,1.d0,ic12,ic2)
      else
       wp0=wp0+wp1+wpi
       wm0=wm0+wm1+wmi
      endif
      if(debug.ge.4)write (moniou,202)wp0,wm0
      
201   format(2x,'qgstr: wpi0=',e10.3,2x,'wmi0=',e10.3
     *,2x,'wp0=',e10.3,2x,'wm0=',e10.3)
202   format(2x,'qgstr - returned light cone momenta:'
     *,2x,'wp0=',e10.3,2x,'wm0=',e10.3)
      return
      end

c===========================================================================
      double precision function qgt(b)
c---------------------------------------------------------------------------
c qgt - nuclear profile function at squared impact parameter b
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)b     
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
      qgt=qgt*zm*0.5d0

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
      qgt=qgt+dt*wsnuc(2)/2.d0
      if(debug.ge.4)write (moniou,202)qgt
      
201   format(2x,'qgt - nuclear profile function at impact'
     *,' parameter b^2=',e10.3)
202   format(2x,'qgt=',e10.3)
      return
      end
      
c-------------------------------------------------------------------------------
**anfe: renamed to qgcrossc to comply with qgsII-04
      subroutine qgcrossc(niter,gtot,gprod,gabs,gdd,gqel,gcoh)
c-------------------------------------------------------------------------------
c qgcrossc - nucleus-nucleus (nucleus-hydrogen) interaction cross sections
c gtot  - total cross section
c gprod - production cross section (projectile diffraction included)
c gabs  - cut Pomerons cross section
c gdd   - projectile diffraction cross section
c gqel  - quasielastic (projectile nucleon knock-out) cross section
c gcoh  - coherent (elastic with respect to the projectile) cross section
c (target diffraction is not treated explicitely and contributes to
c gdd, gqel, gcoh).
c-------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207)
cdh   dimension wabs(28),wdd(28),wqel(28),wcoh(28),wtot(28)
      dimension wabs(28),wdd(28),wqel(28),wcoh(28)
     *,wprod(28),b0(28),ai(28),xa(iapmax,3),xb(iapmax,3)
      common /qgarr1/  ia(2),icz,icp
      common /qgarr5/  rnuc(2),wsnuc(2),wbnuc(2),anorm
     *,cr1(2),cr2(2),cr3(2)
      common /qgarr6/  pi,bm,amws
      common /qgarr11/ b10
      common /qgarr16/ cc(2,3),iddp(iapmax),iddt(iapmax)
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug
      external psran
cf2py intent(in) :: niter
cf2py intent(out) :: gtot,gprod,gabs,gdd,gqel,gcoh

      if(debug.ge.3)write (moniou,201)niter
      do i=1,7
       b0(15-i)=bm*sqrt((1.d0+x1(i))/2.d0)
       b0(i)=bm*sqrt((1.d0-x1(i))/2.d0)
       ai(i)=a1(i)*bm**2*5.d0*pi
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
        iddt(i)=1+int(psran(b10)+cc(2,2))
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
       wabs(i)=wabs(i)/niter
       wdd(i)=wdd(i)/niter
       wqel(i)=wqel(i)/niter
       wcoh(i)=wcoh(i)/niter
       wprod(i)=wabs(i)+wdd(i)
       gabs=gabs+ai(i)*wabs(i)
       gdd=gdd+ai(i)*wdd(i)
       gqel=gqel+ai(i)*wqel(i)
       gcoh=gcoh+ai(i)*wcoh(i)
      enddo
      gprod=gabs+gdd
      gtot=gprod+gqel+gcoh
      if(debug.ge.4)write (moniou,202)
      
201   format(2x,'crossc - nucleus-nucleus interaction cross sections,'
     *,' N of iter.:',i5)
202   format(2x,'crossc: gtot=',e10.3,2x,'gprod=',e10.3,2x,'gabs=',e10.3
     */4x,'gdd=',e10.3,2x,'gqel=',e10.3,2x,'gcoh=',e10.3)
      return
      end

c-------------------------------------------------------------------------------
      subroutine qggcr(b,gabs,gdd,gqel,gcoh,xa,xb,ia)
c-------------------------------------------------------------------------------
c qggcr - integrands (b-profiles) for nucleus-nucleus cross sections
c b - impact parameter
c-------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      parameter(iapmax=207)
      dimension xa(iapmax,3),xb(iapmax,3),vabs(2)
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)b
      gabs=1.
      gdd=1.
      gqel=1.
      gcoh=1.
      do n=1,ia
       call qgv(xa(n,1)+b,xa(n,2),xb,vin,vdd,vabs)
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
      if(debug.ge.4)write (moniou,202)
      
201   format(2x,'qggcr - integrands for nucleus-nucleus cross sections,'
     *,' b=',e10.3)
202   format(2x,'qggcr: gabs=',e10.3,2x,'gdd=',e10.3,2x,'gqel=',e10.3
     *,2x,'gcoh=',e10.3)
      return
      end

c-------------------------------------------------------------------------------
      double precision function qgsect(e0n,icz,iap,iat)
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
      common /qgarr43/ moniou
      common /qgarr47/ gsect(10,5,6)
      common /qgarr48/ asect(10,6,6)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)e0n,icz,iap,iat
      qgsect=0.d0
      ye=dlog10(e0n)
      if(ye.lt.1.d0)ye=1.d0
      je=int(ye)
      if(je.gt.8)je=8

      wk(2)=ye-je
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)

      yb=iat
      yb=dlog(yb)/1.38629d0+1.d0
      jb=min(int(yb),2)
      wb(2)=yb-jb
      wb(3)=wb(2)*(wb(2)-1.d0)*.5d0
      wb(1)=1.d0-wb(2)+wb(3)
      wb(2)=wb(2)-2.d0*wb(3)

      if(iap.eq.1)then
       if(iat.eq.14)then
        do i=1,3
         qgsect=qgsect+gsect(je+i-1,icz,5)*wk(i)
        enddo
       elseif(iat.eq.40)then
        do i=1,3
         qgsect=qgsect+gsect(je+i-1,icz,6)*wk(i)
        enddo
       else
        do i=1,3
        do l=1,3
         qgsect=qgsect+gsect(je+i-1,icz,jb+l-1)*wk(i)*wb(l)
        enddo
        enddo
       endif
      else
       ya=iap
       ya=dlog(ya/2.d0)/.69315d0+1.d0
       ja=min(int(ya),4)
       wa(2)=ya-ja
       wa(3)=wa(2)*(wa(2)-1.d0)*.5d0
       wa(1)=1.d0-wa(2)+wa(3)
       wa(2)=wa(2)-2.d0*wa(3)
       if(iat.eq.14)then
        do i=1,3
        do m=1,3
         qgsect=qgsect+asect(je+i-1,ja+m-1,5)*wk(i)*wa(m)
        enddo
        enddo
       elseif(iat.eq.40)then
        do i=1,3
        do m=1,3
         qgsect=qgsect+asect(je+i-1,ja+m-1,6)*wk(i)*wa(m)
        enddo
        enddo
       else
        do i=1,3
        do m=1,3
        do l=1,3
         qgsect=qgsect+asect(je+i-1,ja+m-1,jb+l-1)*wk(i)*wa(m)*wb(l)
        enddo
        enddo
        enddo
       endif
      endif
      qgsect=exp(qgsect)
      if(debug.ge.4)write (moniou,202)
      
201   format(2x,'qgsect - nucleus-nucleus production cross section'
     */4x,'lab. energy per nucleon - ',e10.3,2x,'hadron class - ',i2
     */4x,'proj. mass N - ',i3,2x,'targ. mass N - ',i3)
202   format(2x,'qgsect=',e10.3)
      return
      end
            
c-----------------------------------------------------------------------------
      subroutine qgdec2(ep,ep1,ep2,ww,a,b)
c-----------------------------------------------------------------------------
c qgdec2 - two particle decay
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension ep(4),ep1(4),ep2(4),ey(3)
      common /qgarr11/ b10
      common /qgarr43/ moniou
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)ep,ww,a,b
      pl=qglam(ww,a,b)
      ep1(1)=dsqrt(pl+a)
      ep2(1)=dsqrt(pl+b)
      pl=dsqrt(pl)
      cosz=2.d0*psran(b10)-1.d0
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
      if(debug.ge.4)write (moniou,202)
      
201   format(2x,'qgdec2: 4-momentum:',2x,4(e10.3,1x)
     */4x,'ww=',e10.3,2x,'a=',e10.3,2x,'b=',e10.3)
202   format(2x,'qgdec2 - end')
      return
      end
      
c=============================================================================
      block data qgdata
c-----------------------------------------------------------------------------
c constants for numerical integration (gaussian weights), 
c parameters for nuclear profiles (3-parameter Fermi)
c-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /qgarr49/ trnuc(56),twsnuc(56),twbnuc(56)
      common /qgarr50/ x1(7),a1(7)
      common /qgarr51/ x4(2),a4(2)
      common /qgarr52/ x2(4),a2
      common /qgarr53/ x9(3),a9(3)
      data x1/.9862838d0,.9284349d0,.8272013d0,.6872929d0,.5152486d0,
     *.3191124d0,.1080549d0/
      data a1/.03511946d0,.08015809d0,.1215186d0,.1572032d0,
     *.1855384d0,.2051985d0,.2152639d0/
      data x2/.00960736d0,.0842652d0,.222215d0,.402455d0/
      data a2/.392699d0/
      data x4/ 0.339981,0.861136/
      data a4/ 0.652145,0.347855/
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

c------------------------------------------------------------------------
      double precision function qggrv(x,qqs,icq,iq)
c------------------------------------------------------------------------
c qggrv - GRV PDFs
c x   - parton LC momentum share,
c qqs - parton virtuality,
c iq  - parton type (flavor),
c icq - hadron type
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr25/ ahv(3)
      common /qgarr43/ moniou
      common /debug/   debug
      
      if(debug.ge.3)write (moniou,201)x,qqs,iq,icq
      if(x.gt..99999d0.and.(qqs.ne.qt0.or.iq.ne.1.and.iq.ne.2))then
       qggrv=0.d0
       return
      endif
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
      if(debug.ge.4)write (moniou,202)qggrv
      
201   format(2x,'qggrv - GRV PDFS: parton x =',e10.3,2x,'q^2=',e10.3
     *,2x,'type=',i2,2x,'hadron type=',i2)
202   format(2x,'qggrv=',e10.3)
      return
      end
      
c-----------------------------------------------------------------------
      real function gamfun(x)
c-----------------------------------------------------------------------
c     gamma function
c-----------------------------------------------------------------------
      dimension c(13)
      data c
     1/ 0.00053 96989 58808, 0.00261 93072 82746, 0.02044 96308 23590,
     2  0.07309 48364 14370, 0.27964 36915 78538, 0.55338 76923 85769,
     3  0.99999 99999 99998,-0.00083 27247 08684, 0.00469 86580 79622,
     4  0.02252 38347 47260,-0.17044 79328 74746,-0.05681 03350 86194,
     5  1.13060 33572 86556/
      gamfun=0
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
    4 gamfun=
     1 f*((((((c(1)*z+c(2))*z+c(3))*z+c(4))*z+c(5))*z+c(6))*z+c(7))/
     2   ((((((c(8)*z+c(9))*z+c(10))*z+c(11))*z+c(12))*z+c(13))*z+1.0)
      if(x .gt. 0.0) return
      gamfun=3.141592653589793/(sin(3.141592653589793*x)*gamfun)
      return
    5 write(*,10)x
   10 format(1x,'argument of gamma fctn = ',e20.5)
      stop
      end
      
c------------------------------------------------------------------------
      double precision function qgev(q1,qj,qq,xx,j,l)
c------------------------------------------------------------------------
c qgev - qcd evolution factor
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr32/ epsxmn
      common /qgarr43/ moniou
      common /qgarr50/ x1(7),a1(7)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)xx,q1,qj,qq,m,l

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
      if(debug.ge.4)write (moniou,202)qgev
      
201   format(2x,'qgev - qcd evolution factor:'
     */4x,'xx=',e10.3,2x,'q1=',e10.3,2x,'qj=',e10.3,2x,'qq=',e10.3
     *,2x,'m=',i1,2x,'l=',i1)
202   format(2x,'qgev=',e10.3)
      return
      end
      
c------------------------------------------------------------------------
      double precision function qgevi(q1,qq,xx,m,l)
c------------------------------------------------------------------------
c qgevi - qcd evolution factor - interpolation 
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer debug
      dimension wi(3),wj(3),wk(3)
      common /qgarr18/ alm,qt0,qtf,betp,dgqq
      common /qgarr20/ spmax
      common /qgarr32/ epsxmn
      common /qgarr43/ moniou
      common /qgarr54/ evk(40,40,100,3,2)
      common /debug/   debug

      if(debug.ge.3)write (moniou,201)xx,q1,qq,m,l
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
      if(debug.ge.4)write (moniou,202)qgevi
      
201   format(2x,'qgevi - interpolation of qcd evolution factor:'
     */4x,'xx=',e10.3,2x,'q1=',e10.3,2x,'qq=',e10.3,2x,2x,'m=',i1
     *,2x,'l=',i1)
202   format(2x,'qgevi=',e10.3)
      return
      end
