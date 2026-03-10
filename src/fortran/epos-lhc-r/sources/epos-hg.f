C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
      subroutine hgcaaa
c-----------------------------------------------------------------------
c hadronic resonance gas in grand canonical treatment
c returns T, chemical potentials and hadronic yield
c (hadron chemical potentials as combinations of quark chemical potentials)
c
c input:
c   iostat: 1: Boltzmann approximation, 0: quantum statistics  /metr3/
c   tecm:                    droplet energy      /confg/
c   volu:                    droplet volume      /confg/
c   keu ked kes kec keb ket: net flavor number   /drop5/
c
c output:
c   tem    : temperature [GeV]                            /cgchg/
c   chem(1:nflav): quark chem. pot. [GeV]                 /cflav/
c   chemgc(1:nspecs): hadron chem. pot. [GeV]             /cgchg/
c   ptlngc(1:nspecs): hadron number                       /cgchg/
c   rmsngc(1:nspecs): standard deviation of hadron number /cgchg/
c
c exact treatment (iostat=0):
c for massive hadrons     : first in Boltzmann approximation with analytical
c                           expressions for particle and energy densities,
c                           then by using quantum statistics in integral form,
c                           extracting mu and T using numerical integration
c                           and an iterative procedure to solve for mu, T
c for massless hadrons    : using analytic expressions for massles particles
c                           and employing the  same algorithm as for massive
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cbol/rmsbol(mspecs),ptlbol(mspecs),chebol(mspecs),tembol
      common/cflavs/nflavs,kef(nflav),chem(nflav)
      common/ciakt/gen,iafs,ians,genm
      common/cnrit/nrit
      gen=10.0**(-epsgc)
      genm=gen/10.

      isho=ish
      if(ishsub/100.eq.51)ish=mod(ishsub,100)

      iug=(1+iospec)/2*2-1


c     initialization
c     --------------
      kef(1)=keu
      kef(2)=ked
      kef(3)=kes
      kef(4)=kec
      kef(5)=keb
      kef(6)=ket

      nflavs=3
      if(iug.eq.1)nflavs=1
      if(iug.eq.3)nflavs=2
      if(iug.eq.5)nflavs=2
      if(iug.eq.7)nflavs=3
      if(iug.eq.9)nflavs=3
      if(iug.eq.11)nflavs=3
      tem=0.0
      do i=1,nflavs
      chem(i)=0.0
      enddo
      call hgchac(0)
      do i=1,nspecs
      ptlngc(i)=0.0
      rmsngc(i)=0.0
      enddo
      nrit=0

      if(ish.ge.5)then
        write(ifch,*)('-',l=1,10)
     *,' entry sr hgcaaa ',('-',l=1,30)
        write(ifch,'(1x,a,2x,3i3)')
     *'>>> grand canonical hadron gas for droplet with u d s content:'
     *,keu,ked,kes
        write(ifch,'(1x,a,2x,f7.3,2x,a,2x,f7.3)')
     *'mass [GeV]:',tecm,'volume [fm^3]:',volu
      endif

      if(iug.eq.1.and.keu.ne.0.and.ish.ge.5)then
      write(ifch,*)'inversion impossible !!!'
      write(ifch,*)'keu=0 required for this option'
      write(ifch,*)'T = mu(i) = 0 returned'
      if(ish.ge.5)write(ifch,*)('-',i=1,30)
     *,' exit sr hgcaaa ',('-',i=1,10)
      return
      endif
      if(iug.eq.3.and.(keu+ked).ne.0.and.ish.ge.5)then
      write(ifch,*)'inversion impossible !!!'
      write(ifch,*)'keu+ked=0 required for this option'
      write(ifch,*)'T = mu(i) = 0 returned'
      if(ish.ge.5)write(ifch,*)('-',i=1,30)
     *,' exit sr hgcaaa ',('-',i=1,10)
      return
      endif
      kf=keu+ked+kes+kec+keb+ket
      kf=abs(kf)
      if(kf.ne.0)then
      if(mod(kf,3).ne.0.and.ish.ge.5)then
      write(ifch,*)'inversion impossible !!!'
      write(ifch,*)'sum must be multiple of three'
      write(ifch,*)'T = mu(i) = 0 returned'
      if(ish.ge.5)write(ifch,*)('-',i=1,30)
     *,' exit sr hgcaaa ',('-',i=1,10)
      return
      endif
      endif


c     initial T (m=0, baryon free)
c     -------------------------------
      gfac=0.0

       if(iostat.eq.0.and.iospec.eq.iug)then
      do i=1,nspecs
      igsp=int(gspecs(i))
      if(mod(igsp,2).eq.0)then
      gfac=gfac+7.*gspecs(i)/8.
      else
      gfac=gfac+gspecs(i)
      endif
      enddo
      if(iabs(ispecs(nspecs)).lt.10)gfac=gfac+16.
      tem=(tecm/volu*hquer**3*30./pi**2/gfac)**.25
       else
      do i=1,nspecs
      gfac=gfac+gspecs(i)
      enddo
      if(iabs(ispecs(nspecs)).lt.10)gfac=gfac+16.
      tem=(tecm/volu*hquer**3*pi**2/gfac/3.)**.25
      tem=2.*tem
       endif

      if(ish.ge.5)write(ifch,1)'initial T :',tem
1     format(1x,a,3x,f9.6)

      if(ish.ge.5)write(ifch,*)'iospec: ',iospec

       if(ish.ge.5.and.iospec.ne.iug)then
      write(ifch,*)'inversion in Boltzmann approx. :'
       elseif(ish.ge.5.and.iospec.eq.iug)then
      write(ifch,*)'inversion for massless hadrons :'
       endif

       if(ish.ge.5)then
      if(nflavs.eq.1)write(ifch,'(3x,a,8x,a)')
     *'T:','chem u:'
      if(nflavs.eq.2)write(ifch,'(3x,a,8x,a,5x,a)')
     *'T:','chem u:','chem d:'
      if(nflavs.eq.3)write(ifch,'(3x,a,8x,a,5x,a,5x,a)')
     *'T:','chem u:','chem d:','chem s:'
       endif

      k=1
10    continue
      if(ish.ge.9.and.mod(k,10).eq.0)
     *write(ifch,*)'hgc iteration:',k
      if(ish.ge.9)call hgccch(1)

c     search for temperature (chem=const)
c     -----------------------------------
      idt=0
      temo=tem

       if(iospec.eq.iug)then

c     massless particles
c     ------------------
      if(iostat.eq.0)then
      if(ish.ge.9)
     *write(ifch,*)'iteration (massless):',k
      call hgctm0
      elseif(iostat.eq.1)then
      if(ish.ge.9)
     *write(ifch,*)'iteration (Boltzmann, massless):',k
      call hgctbo(ibna)
      if(ibna.eq.1)then
      tem=temo
      goto20
      endif
      endif

       else

c     Boltzmann approxiamtion (massive particles)
c     -------------------------------------------
      if(ish.ge.9)
     *write(ifch,*)'iteration (Boltzmann, massive):',k
      call hgctbo(ibna)
      if(ibna.eq.1)then
      tem=temo
      goto20
      endif

       endif

      if(tem.le.1.e-6.and.ish.ge.5)then
      write(ifch,*)'inversion imposssible'
      write(ifch,*)'T:',tem
      if(ioinco.ge.1)call hnbmin(iospec,keu,ked,kes,kec)
      if(ish.ge.5)write(ifch,*)('-',i=1,30)
     *,' exit sr hgcaaa ',('-',i=1,10)
      ish=isho
      return
      endif

      dt=abs(temo-tem)
      if(dt.le.gen*temo.or.dt.le.genm)idt=1

c     search for chemical potentials (tem=const)
c     ------------------------------------------
      idch=0
      ibna=0

        do iafs=1,nflavs
      chemo=chem(iafs)

       if(iospec.eq.iug)then

c     massless particles
c     ------------------
      if(iostat.eq.0)then
      call hgccm0
      elseif(iostat.eq.1)then
      call hgccbo(ibna)
      endif

       else

c     Boltzmann approxiamtion (massive particles)
c     -------------------------------------------
      call hgccbo(ibna)

       endif

      dch=abs(chemo-chem(iafs))
      if(ish.ge.9)write(ifch,*)'dch:',dch
      if(dch.le.abs(gen*chemo).or.dch.le.genm)idch=idch+1
      if(ibna.eq.1)then
      chem(iafs)=chemo
      call hgchac(0)
      goto20
      endif

        enddo


c     new hadron chem. potentials
c     ---------------------------
      call hgchac(0)


      if(ish.ge.5.and.nflavs.eq.1)
     *write(ifch,'(1x,f8.6,2x,f9.6)')
     *tem,chem(1)
      if(ish.ge.5.and.nflavs.eq.2)
     *write(ifch,'(1x,f8.6,2x,f9.6,2x,f9.6)')
     *tem,chem(1),chem(2)
      if(ish.ge.5.and.nflavs.eq.3)
     *write(ifch,'(1x,f8.6,2x,f9.6,2x,f9.6,2x,f9.6)')
     *tem,chem(1),chem(2),chem(3)
      if(idch.eq.nflavs.and.idt.eq.1)goto20


      k=k+1

       if(k.gt.300)then
       if(ish.ge.5)
     *write(ifch,*)'failure in approximate solution'
      goto20
       endif

      goto10

20    continue
      if(ish.ge.9)call hgccch(0)
      if(ish.ge.5)write(ifch,'(1x,a,1x,f9.6)')'  T  :',tem
      do i=1,nflavs
      if(i.eq.1.and.ish.ge.5)
     *write(ifch,'(1x,a,1x,f9.6)')'chem u:',chem(1)
      if(i.eq.2.and.ish.ge.5)
     *write(ifch,'(1x,a,1x,f9.6)')'chem d:',chem(2)
      if(i.eq.3.and.ish.ge.5)
     *write(ifch,'(1x,a,1x,f9.6)')'chem s:',chem(3)
      enddo


c     checking results
c     ----------------
      if(ish.ge.5)call hgcchb

c     particle yield
c     --------------
      call hgcpyi(1)

c     checking flavor conservation
c     ----------------------------
      if(ish.ge.5)call hgccfc

      if(iug.eq.iospec.and.iostat.eq.0)then
      if(ish.ge.5)write(ifch,*)
     *'approximation and exact treatment equal'
      if(ish.ge.5)write(ifch,*)('-',i=1,30)
     *,' exit sr hgcaaa ',('-',i=1,10)
      ish=isho
      return
      endif

c     continue or return approximate values
c     -------------------------------------
      do i=1,nspecs
      rmsbol(i)=rmsngc(i)
      ptlbol(i)=ptlngc(i)
      chebol(i)=chemgc(i)
      enddo
      tembol=tem
      if(iostat.eq.1)then
      if(ish.ge.5)write(ifch,*)('-',i=1,30)
     *,' exit sr hgcaaa ',('-',i=1,10)
      ish=isho
      return
      endif


c     quantum statistics
c     ------------------
      if(ish.ge.5)write(ifch,*)'quantum statistics:'
      if(ish.ge.5.and.nflavs.eq.1)write(ifch,'(3x,a,8x,a)')
     *'T:','chem u:'
      if(ish.ge.5.and.nflavs.eq.2)write(ifch,'(3x,a,8x,a,6x,a)')
     *'T:','chem u:','chem d:'
      if(ish.ge.5.and.nflavs.eq.3)write(ifch,'(3x,a,8x,a,6x,a,6x,a)')
     *'T:','chem u:','chem d:','chem s:'
      k=1

30    continue
      if(ish.ge.9.and.mod(k,10).eq.0)
     *write(ifch,*)'hgc iteration:',k

c     new temperature
c     ---------------
      idt=0
      temo=tem
      call hgctex
      if(ish.ge.5.and.nflavs.eq.1)
     *write(ifch,'(1x,f10.8,2x,f10.7)')
     *tem,chem(1)
      if(ish.ge.5.and.nflavs.eq.2)
     *write(ifch,'(1x,f10.8,2x,f10.7,2x,f10.7)')
     *tem,chem(1),chem(2)
      if(ish.ge.5.and.nflavs.eq.3)
     *write(ifch,'(1x,f10.8,2x,f10.7,2x,f10.7,2x,f10.7)')
     *tem,chem(1),chem(2),chem(3)

      if(tem.le.1.e-6.and.ish.ge.5)then
      write(ifch,*)'inversion imposssible'
      write(ifch,*)'T:',tem
      call hnbmin(iospec,keu,ked,kes,kec)
      if(ish.ge.5)write(ifch,*)('-',i=1,30)
     *,' exit sr hgcaaa ',('-',i=1,10)
      ish=isho
      return
      endif

      dt=abs(temo-tem)
      if(dt.le.gen*temo.or.dt.le.genm)idt=1
      if(ish.ge.9)write(ifch,*)'dtem:',dt

c     new quark chem. potentials
c     --------------------------
      idch=0
      do iafs=1,nflavs
      chemo=chem(iafs)
      call hgccex
      dch=abs(chemo-chem(iafs))
      if(ish.ge.9)write(ifch,*)'dche:',dch
      if(dch.le.abs(gen*chemo).or.dch.le.genm)idch=idch+1
      enddo

c     new hadron chem. potentials
c     ---------------------------
      call hgchac(0)

       if(idch.eq.nflavs.and.idt.eq.1)then

      if(ish.ge.5)write(ifch,*)'results:'
      if(ish.ge.5)write(ifch,51)'  T  :',tem
      if(nflavs.ge.1.and.ish.ge.5)write(ifch,51)'chem u:',chem(1)
      if(nflavs.ge.2.and.ish.ge.5)write(ifch,51)'chem d:',chem(2)
      if(nflavs.ge.3.and.ish.ge.5)write(ifch,51)'chem s:',chem(3)
51    format(1x,a,3x,f9.6)

c     checking results
c     ----------------
      if(ish.ge.5)call hgcchh(i)

c     particle yield
c     --------------
      call hgcpyi(0)

c     checking flavor conservation
c     ----------------------------
      call hgccfc

      if(ish.ge.5)write(ifch,*)('-',i=1,30)
     *,' exit sr hgcaaa ',('-',i=1,10)
      ish=isho
      return
       endif

       if(k.gt.300)then
       if(ish.ge.5)
     *write(ifch,*)'failure in exact solution'
      if(ish.ge.5)write(ifch,*)'results:'
      if(ish.ge.5)write(ifch,51)'  T  :',tem
      if(nflavs.ge.1.and.ish.ge.5)write(ifch,51)'chem u:',chem(1)
      if(nflavs.ge.2.and.ish.ge.5)write(ifch,51)'chem d:',chem(2)
      if(nflavs.ge.3.and.ish.ge.5)write(ifch,51)'chem s:',chem(3)

c     particle yield
c     --------------
      call hgcpyi(0)

      if(ish.ge.5)write(ifch,*)('-',i=1,30)
     *,' exit sr hgcaaa ',('-',i=1,10)
      ish=isho
      return

       endif

      k=k+1
      goto30

      end


c---------------------------------------------------------------------
      function hgcbi0(x)
c---------------------------------------------------------------------
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     *1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     *0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,
     *-0.1647633d-1,0.392377d-2/
      if (abs(x).lt.3.75) then
        y=dble((x/3.75)**2)
        hgcbi0=sngl(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        ax=abs(x)
        y=dble(3.75/ax)
        hgcbi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *(q7+y*(q8+y*q9))))))))
      endif
      return
      end


c------------------------------------------------------------------------
      function hgcbi1(x)
c------------------------------------------------------------------------
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     *0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     *-0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,-0.2895312d-1,
     *0.1787654d-1,-0.420059d-2/
      if (abs(x).lt.3.75) then
        y=dble((x/3.75)**2)
        hgcbi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        ax=abs(x)
        y=dble(3.75/ax)
        hgcbi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *(q7+y*(q8+y*q9))))))))
        if(x.lt.0.)hgcbi1=-hgcbi1
      endif
      return
      END


c---------------------------------------------------------------------
      function hgcbk0(x)
c------------------------------------------------------------------------
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     *0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     *-0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
      if (x.le.2.0) then
        y=dble(x*x/4.0)
        hgcbk0=(-log(x/2.0)*hgcbi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*
     *(p6+y*p7))))))
      else
        y=dble(2.0/x)
        hgcbk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *q7))))))
      endif
      return
      END


c---------------------------------------------------------------
      function hgcbk1(x)
c--------------------------------------------------------------------
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     *-0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     *0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
      if (x.le.2.0) then
        y=dble(x*x/4.0)
        hgcbk1=(log(x/2.0)*hgcbi1(x))+(1.0/x)*(p1+y*(p2+y*(p3+y*(p4+y*
     *(p5+y*(p6+y*p7))))))
      else
        y=dble(2.0/x)
        hgcbk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *q7))))))
      endif
      return
      END


c-------------------------------------------------------------------
      function hgcbk(n,x)
c------------------------------------------------------------------
      tox=2.0/x
      bkm=hgcbk0(x)
      bk=hgcbk1(x)
      do 11 j=1,n-1
        bkp=bkm+j*tox*bk
        bkm=bk
        bk=bkp
11    continue
      hgcbk=bk
      return
      END


c----------------------------------------------------------------
      subroutine hgccbo(iba)
c----------------------------------------------------------------
c returns new chem(iafs) for boltzmann statistics
c  input:
c    tem
c    kef/volu
c  output:
c    chem(iafs)
c-----------------------------------------------------------------------
      common/cnsta/pi,pii,hquer,prom,piom,ainfin
      common/drop6/tecm,volu
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      parameter(nflav=6)
      common/cflavs/nflavs,kef(nflav),chem(nflav)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)
      common/ciakt/gen,iafs,ians,genm
      external hgcbk
      k=1
      iba=0
      c1=-0.5
      c2=0.5
      goto11

c     new chemical potential
c     ----------------------
10    chem(iafs)=c1+0.5*(c2-c1)
11    continue
      fd=0.0
      call hgchac(0)

        do i=1,nspecs

        if(ifok(iafs,i).ne.0)then
       if((chemgc(i)/tem).gt.70.)then
      hpd=1.e30
       else
      hpd=exp(chemgc(i)/tem)
       endif
       if(aspecs(i).ne.0.)then
      fk2=hgcbk(2,aspecs(i)/tem)
      hpd=hpd*gspecs(i)*aspecs(i)**2*tem*fk2
     */2./pi**2/hquer**3
       else
      hpd=hpd*gspecs(i)*tem**3/pi**2/hquer**3
       endif
      hfd=ifok(iafs,i)*hpd
      fd=fd+hfd
        endif

        enddo

      dfd=abs(fd-(kef(iafs)/volu))
      if(dfd.le.abs(gen*(kef(iafs)/volu)).or.dfd.le.genm)return
c     if(abs(fd).ge.100.)then
c     iba=1
c     return
c     endif


       if(fd.gt.(kef(iafs)/volu))then
      c2=chem(iafs)
      else
      c1=chem(iafs)
       endif

      k=k+1
      if(k.gt.300)return

      goto10

      end


c----------------------------------------------------------------------
      subroutine hgccch(iii)
c----------------------------------------------------------------------
c checks convergence of iterative algorithm
c plots iteration values for T and mu_i
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cflavs/nflavs,kef(nflav),chem(nflav)
      parameter (nbin=500)
      common/cdatc/data(nbin),datb(nbin),datc(nbin),datd(nbin)
     *,date(nbin),datf(nbin),datg(nbin),dath(nbin),dati(nbin)
      common/cnrit/nrit
      character cen*4,cvol*4,cu*3,cd*3,cs*3

           if(iii.gt.0)then

      nrit=nrit+1
      data(nrit)=nrit
      datb(nrit)=tem
      datc(nrit)=chem(1)
      datd(nrit)=chem(2)
      date(nrit)=chem(3)

           elseif(iii.eq.0)then

      nrit=nrit+1
      data(nrit)=nrit
      datb(nrit)=tem
      datc(nrit)=chem(1)
      datd(nrit)=chem(2)
      date(nrit)=chem(3)
      do i=1,nrit
      datf(i)=datb(nrit)
      datg(i)=datc(nrit)
      dath(i)=datd(nrit)
      dati(i)=date(nrit)
      enddo

      x1=data(1)
      x2=data(nrit)
      write(cen,'(f4.1)')tecm
      write(cvol,'(f4.1)')volu
      write(cu,'(i3)')keu
      write(cd,'(i3)')ked
      write(cs,'(i3)')kes


      write(ifhi,'(a)')       'newpage zone 1 4 1 openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')       'text 0 0 "xaxis Iteration"'
      write(ifhi,'(a)')       'text 0 0 "yaxis T (GeV)"'
      write(ifhi,'(a)')       'text 0.15 0.9 "E= '//cen//'"'
      write(ifhi,'(a)')       'text 0.4 0.9 "V= '//cvol//'"'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(3a)')'yrange',' auto',' auto'
      write(ifhi,'(a)')       'array 2'
      do j=1,nrit
      write(ifhi,'(2e12.4)')data(j),datb(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0-'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(3a)')'yrange',' auto',' auto'
      write(ifhi,'(a)')       'array 2'
      do j=1,nrit
      write(ifhi,'(2e12.4)')data(j),datf(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')       'text 0 0 "xaxis Iteration"'
      write(ifhi,'(a)')       'text 0 0 "yaxis [m]^1! (GeV)"'
      write(ifhi,'(a)')       'text 0.15 0.9 "Q^1!= '//cu//'"'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(3a)')'yrange',' auto',' auto'
      write(ifhi,'(a)')       'array 2'
      do j=1,nrit
      write(ifhi,'(2e12.4)')data(j),datc(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0-'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(3a)')'yrange',' auto',' auto'
      write(ifhi,'(a)')       'array 2'
      do j=1,nrit
      write(ifhi,'(2e12.4)')data(j),datg(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')       'text 0 0 "xaxis Iteration"'
      write(ifhi,'(a)')       'text 0 0 "yaxis [m]^2! (GeV)"'
      write(ifhi,'(a)')       'text 0.15 0.9 "Q^2!= '//cd//'"'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(3a)')'yrange',' auto',' auto'
      write(ifhi,'(a)')       'array 2'
      do j=1,nrit
      write(ifhi,'(2e12.4)')data(j),datd(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0-'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(3a)')'yrange',' auto',' auto'
      write(ifhi,'(a)')       'array 2'
      do j=1,nrit
      write(ifhi,'(2e12.4)')data(j),dath(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')       'text 0 0 "xaxis Iteration"'
      write(ifhi,'(a)')       'text 0 0 "yaxis [m]^3! (GeV)"'
      write(ifhi,'(a)')       'text 0.15 0.9 "Q^3!= '//cs//'"'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(3a)')'yrange',' auto',' auto'
      write(ifhi,'(a)')       'array 2'
      do j=1,nrit
      write(ifhi,'(2e12.4)')data(j),date(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0-'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(3a)')'yrange',' auto',' auto'
      write(ifhi,'(a)')       'array 2'
      do j=1,nrit
      write(ifhi,'(2e12.4)')data(j),dati(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0'

           endif

       return

       end

c-----------------------------------------------------------------------
      subroutine hgccex
c-----------------------------------------------------------------------
c returns new chem(iafs) for massive quantum statistics
c  input:
c    tem
c    kef/volu
c  output:
c    chem(iafs)
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cflavs/nflavs,kef(nflav),chem(nflav)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)
      common/ciakt/gen,iafs,ians,genm
      external hgcfhn

      k=1

      c1=-0.5
      c2=0.5
      goto11

c     new chemical potential
c     ----------------------
10    chem(iafs)=c1+0.5*(c2-c1)
11    continue

      fd=0.0
        do ians=1,nspecs
       if(ifok(iafs,ians).ne.0)then

      call hgchac(0)
      call hgclim(a,b)
      if(b.eq.0.0)then
      hpd=0.0
      else
      call uttraq(hgcfhn,a,b,hpd)
      endif
      hpd=hpd*gspecs(ians)/2./pi**2/hquer**3
      fd=fd+hpd*ifok(iafs,ians)

       endif
        enddo

      dfd=abs(fd-(kef(iafs)/volu))
      if(dfd.le.abs(gen*(kef(iafs)/volu)).or.dfd.le.genm)return

       if(fd.gt.(kef(iafs)/volu))then
      c2=chem(iafs)
      else
      c1=chem(iafs)
       endif

      k=k+1
      if(k.gt.300)then
      if(ish.ge.5)
     *write(ifch,*)'failure at cex at iafs:',iafs
      return
      endif

      goto10

      end


c------------------------------------------------------------------
      subroutine hgccfc
c------------------------------------------------------------------
c checks flavor conservation in particle yield
c------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cflavs/nflavs,kef(nflav),chem(nflav)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)

      if(ish.ge.5)write(ifch,*)'checking flavor conservation'
      do i=1,nflavs
      ckef=0.0
      do ii=1,nspecs
      ckef=ckef+ifok(i,ii)*ptlngc(ii)
      enddo
      dkef=abs(ckef-kef(i))
      if(dkef.le.1.e-2)then
      if(i.eq.1.and.ish.ge.5)write(ifch,*)'u conserved'
      if(i.eq.2.and.ish.ge.5)write(ifch,*)'d conserved'
      if(i.eq.3.and.ish.ge.5)write(ifch,*)'s conserved'
      else
      if(i.eq.1.and.ish.ge.5)write(ifch,*)'u not conserved'
      if(i.eq.2.and.ish.ge.5)write(ifch,*)'d not conserved'
      if(i.eq.3.and.ish.ge.5)write(ifch,*)'s not conserved'
      if(ish.ge.5)write(ifch,*)'df:',dkef
      endif
      enddo

      return
      end

c----------------------------------------------------------------
      subroutine hgcchb
c----------------------------------------------------------------
c checks results by numerical integration
c----------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cflavs/nflavs,kef(nflav),chem(nflav)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)
      common/ciakt/gen,iafs,ians,genm
      external hgcfbe
      external hgcfbn
      if(ish.ge.5)write(ifch,*)
     *'check by numer. calc. of expect. values:'
      iced=0
      ceden=0.0
       do ians=1,nspecs
      call hgclim(a,b)
      if(b.eq.0.0)then
      cedh=0.0
      else
      call uttraq(hgcfbe,a,b,cedh)
      endif
      if(ish.ge.9)write(ifch,*)'cedh:',cedh
      ced=cedh*gspecs(ians)/2./pi**2/hquer**3
      ceden=ceden+ced
       enddo

      if(iabs(ispecs(nspecs)).lt.10)
     *ceden=ceden+(8.*pi**2*tem**4/15.+bag4rt**4)/hquer**3

      if(ish.ge.5)write(ifch,*)'energy density   :',ceden
      ded=abs((tecm/volu)-ceden)
      if((tecm/volu)*gen.ge.ded.or.ded.le.gen)iced=1
      icfd=0

       do i=1,nflavs
      cfd=0.0
      do ians=1,nspecs
      call hgclim(a,b)
      if(b.eq.0.0)then
      hpd=0.0
      else
      call uttraq(hgcfbn,a,b,hpd)
      endif
      hfd=ifok(i,ians)*hpd*gspecs(ians)/2./pi**2/hquer**3
      if(ish.ge.9)write(ifch,*)'hfd:',hfd
      cfd=cfd+hfd
      enddo
      if(i.eq.1.and.ish.ge.5)write(ifch,5)'flavor density u :',cfd
      if(i.eq.2.and.ish.ge.5)write(ifch,5)'flavor density d :',cfd
      if(i.eq.3.and.ish.ge.5)write(ifch,5)'flavor density s :',cfd
5     format(1x,a,1x,f12.6)
      dfd=abs(cfd-(kef(i)/volu))
      if(abs(gen*(kef(i)/volu)).ge.dfd.or.dfd.le.gen)
     *icfd=icfd+1
       enddo

       if(iced.eq.1.and.icfd.eq.nflavs)then
      if(ish.ge.5)write(ifch,*)'results agree'
      else
      if(ish.ge.5)write(ifch,*)'results disagree'
       endif

      return
      end

c----------------------------------------------------------------
      subroutine hgcchh(icorr)
c----------------------------------------------------------------
c checks results by numerical integration
c----------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cflavs/nflavs,kef(nflav),chem(nflav)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)
      common/ciakt/gen,iafs,ians,genm
      external hgcfhe
      external hgcfhn
      icorr=0
      if(ish.ge.5)write(ifch,*)
     *'check by numer. calc. of expect. values:'

      iced=0
      ceden=0.0
       do ians=1,nspecs
      call hgclim(a,b)
      if(b.eq.0.0)then
      cedh=0.0
      else
      call uttraq(hgcfhe,a,b,cedh)
      endif
      if(ish.ge.9)write(ifch,*)'cedh:',cedh
      ced=cedh*gspecs(ians)/2./pi**2/hquer**3
      ceden=ceden+ced
       enddo

      if(iabs(ispecs(nspecs)).lt.10)
     *ceden=ceden+(8.*pi**2*tem**4/15.+bag4rt**4)/hquer**3

      if(ish.ge.5)write(ifch,*)'energy density   :',ceden
      ded=abs((tecm/volu)-ceden)
      if((tecm/volu)*gen.ge.ded.or.ded.le.gen)iced=1

      icfd=0

       do i=1,nflavs
      cfd=0.0
      do ians=1,nspecs
      call hgclim(a,b)
      if(b.eq.0.0)then
      hpd=0.0
      else
      call uttraq(hgcfhn,a,b,hpd)
      endif
      hfd=ifok(i,ians)*hpd*gspecs(ians)/2./pi**2/hquer**3
      if(ish.ge.9)write(ifch,*)'hfd:',hfd
      cfd=cfd+hfd
      enddo
      if(i.eq.1.and.ish.ge.5)write(ifch,5)'flavor density u :',cfd
      if(i.eq.2.and.ish.ge.5)write(ifch,5)'flavor density d :',cfd
      if(i.eq.3.and.ish.ge.5)write(ifch,5)'flavor density s :',cfd
5     format(1x,a,1x,f9.6)
      dfd=abs(cfd-(kef(i)/volu))
      if(abs(gen*(kef(i)/volu)).ge.dfd.or.dfd.le.gen)
     *icfd=icfd+1
       enddo

       if(iced.eq.1.and.icfd.eq.nflavs)then
      if(ish.ge.5)write(ifch,*)'results agree'
      icorr=1
      else
      if(ish.ge.5)write(ifch,*)'results disagree'
       endif

      return
      end


c--------------------------------------------------------------------
      subroutine hgccm0
c--------------------------------------------------------------------
c returns new quark chemical potentials for massless quantum statistics
c input:
c  tem
c  kef/volu
c output:
c  chem
c---------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cflavs/nflavs,kef(nflav),chem(nflav)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)
      common/ciakt/gen,iafs,ians,genm
      external hgcfhn
      k=1
      z3=1.2020569

      c1=-0.5
      c2=0.5
      goto11

c     new chemical potential
c     ----------------------
10    chem(iafs)=c1+0.5*(c2-c1)
11    continue

      fd=0.0
      call hgchac(0)

               do i=1,nspecs
              if(ifok(iafs,i).ne.0)then

           igsp=int(gspecs(i))
          if(mod(igsp,2).eq.0)then

       if(ispecs(i).gt.0)then
      hpd=gspecs(i)*(chemgc(i)*tem**2+chemgc(i)**3/pi**2)/6./hquer**3
       else
      hpd=0.0
       endif

c            else
c      if(ispecs(i).gt.0)then
c     hpd=gspecs(i)*(chemgc(i)*tem**2/3.-chemgc(i)**3/pi**2/6.)/hquer**3
c      else
c     hpd=0.0
c      endif
c        endif

c     n=1
c0    xx=n*abs(chemgc(i))/tem
c     if(xx.le.60.)then
c     hpd=hpd+(-1.)**(n+1)/n**3/exp(xx)
c     n=n+1
c     goto20
c     endif
c     hpd=hpd*gspecs(i)*tem**3/pi**2/hquer**3
c     if(chemgc(i).eq.abs(chemgc(i)))then
c     hpd=gspecs(i)*(chemgc(i)*tem**2+chemgc(i)**3/pi**2)/6./hquer**3
c    *-hpd
c     endif

c      else
c     hpd=3.*gspecs(i)*tem**3*z3/4./pi**2/hquer**3
c      endif

         else

      hpd=gspecs(i)*tem**3*z3/pi**2/hquer**3

         endif

      hfd=hpd*ifok(iafs,i)
      fd=fd+hfd

       endif
        enddo

      dfd=abs(fd-(kef(iafs)/volu))
      if(dfd.le.abs(gen*(kef(iafs)/volu)).or.dfd.le.genm)return

       if(fd.gt.(kef(iafs)/volu))then
      c2=chem(iafs)
      else
      c1=chem(iafs)
       endif

      k=k+1
      if(k.gt.300)then
      if(ish.ge.5)
     *write(ifch,*)'failure at cm0 at iafs:',iafs
      return
      endif
      goto10
      end

c-----------------------------------------------------------------------
      function hgcfbe(x)
c-----------------------------------------------------------------------
c integrand of energy density
c------------------------------------------------------------------------
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/ciakt/gen,iafs,ians,genm
      eex=81.
      hgcfbe=0.0
      sq=sqrt(x**2+aspecs(ians)**2)
      if(tem.ne.0.0)eex=(sq-chemgc(ians))/tem
      if(eex.gt.60.)return
      if(eex.lt.-60)then
      hgcfbe=1.e25
      return
      endif

      hgcfbe=sq*x**2*exp(-eex)

      return
      end

c-----------------------------------------------------------------
      function hgcfbf(x)
c-----------------------------------------------------------------
c integrand of mean square variance of  energy
c----------------------------------------------------------------
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/ciakt/gen,iafs,ians,genm
      eex=61
      hgcfbf=0.0

      sq=sqrt(x**2+aspecs(ians)**2)
      if(tem.ne.0.0)eex=(sq-chemgc(ians))/tem
      if(eex.gt.60.)return
      if(eex.lt.-60)then
      hgcfbf=1.e25
      return
      endif

      hgcfbf=(aspecs(ians)**2+x**2)*x**2*exp(-eex)

      return
      end

c-----------------------------------------------------------------
      function hgcfbn(x)
c-----------------------------------------------------------------
c integrand of hadron density
c-----------------------------------------------------------------
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/ciakt/gen,iafs,ians,genm
      eex=81.
      hgcfbn=0.0

      sq=sqrt(x**2+aspecs(ians)**2)
      if(tem.ne.0.0)eex=(sq-chemgc(ians))/tem
      if(eex.gt.80.)return
      if(eex.lt.-60)then
      hgcfbn=1.e25
      return
      endif

      hgcfbn=x**2*exp(-eex)

      return
      end

c-----------------------------------------------------------------------
      function hgcfhe(x)
c-----------------------------------------------------------------------
c integrand of energy density
c------------------------------------------------------------------------
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/ciakt/gen,iafs,ians,genm
      eex=81.
      hgcfhe=0.0
      igsp=int(gspecs(ians))

      sq=sqrt(x**2+aspecs(ians)**2)
      if(tem.ne.0.0)eex=(sq-chemgc(ians))/tem
      if(eex.gt.80.)return

       if(mod(igsp,2).ne.0)then
      d=-1.0
      if(eex.lt.1.e-10)return
       else
      d=1.0
       endif

      hgcfhe=sq*x**2/(exp(eex)+d)

      return
      end

c-----------------------------------------------------------------
      function hgcfhf(x)
c-----------------------------------------------------------------
c integrand of mean square variance of  energy
c----------------------------------------------------------------
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/ciakt/gen,iafs,ians,genm
      eex=61
      hgcfhf=0.0
      igsp=int(gspecs(ians))

      sq=sqrt(x**2+aspecs(ians)**2)
      if(tem.ne.0.0)eex=(sq-chemgc(ians))/tem
      if(eex.gt.60.)return
      if(eex.lt.(-60.))return

       if(mod(igsp,2).ne.0)then
      d=-1.0
      if(eex.lt.1.0e-10.and.eex.gt.(-1.0e-10))return
       else
      d=1.0
       endif

      hgcfhf=(aspecs(ians)**2+x**2)*x**2/(exp(eex)+2.0*d+exp(-eex))

      return
      end

c-----------------------------------------------------------------
      function hgcfhn(x)
c-----------------------------------------------------------------
c integrand of hadron density
c-----------------------------------------------------------------
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/ciakt/gen,iafs,ians,genm
      eex=81.
      hgcfhn=0.0
      igsp=int(gspecs(ians))

      sq=sqrt(x**2+aspecs(ians)**2)
      if(tem.ne.0.0)eex=(sq-chemgc(ians))/tem
      if(eex.gt.80.)return

       if(mod(igsp,2).ne.0)then
      d=-1.0
      if(eex.lt.1.e-10)return
       else
      d=1.0
       endif

      hgcfhn=x**2/(exp(eex)+d)

      return
      end

c-----------------------------------------------------------------
      function hgcfhw(x)
c-----------------------------------------------------------------
c integrand of mean square variance of hadron yield
c----------------------------------------------------------------
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/ciakt/gen,iafs,ians,genm
      eex=61
      hgcfhw=0.0
      igsp=int(gspecs(ians))

      sq=sqrt(x**2+aspecs(ians)**2)
      if(tem.ne.0.0)eex=(sq-chemgc(ians))/tem
      if(eex.gt.60.)return
      if(eex.lt.(-60.))return

       if(mod(igsp,2).ne.0)then
      d=-1.0
      if(eex.lt.1.0e-10.and.eex.gt.(-1.0e-10))return
       else
      d=1.0
       endif

      hgcfhw=x**2/(exp(eex)+2.0*d+exp(-eex))

      return
      end


c-----------------------------------------------------------------
      subroutine hgchac(iboco)
c------------------------------------------------------------------
c returns hadronic chemical potentials as combinations of quark
c chemical potentials
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cflavs/nflavs,kef(nflav),chem(nflav)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)

       do i=1,nspecs
      chemgc(i)=0.0
      do ii=1,nflavs
      chemgc(i)=chemgc(i)+ifok(ii,i)*chem(ii)
      if(ish.ge.9)write(ifch,*)'mu_i:',chem(ii),' k_i:',ifok(ii,i)
      enddo
      if(ish.ge.9)write(ifch,*)'mu_nu:',chemgc(i)
      igsp=int(gspecs(i))
      if(mod(igsp,2).ne.0.and.chemgc(i).gt.aspecs(i).and.iboco.eq.0)
     *chemgc(i)=aspecs(i)
       enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine hgclim(a,b)
c----------------------------------------------------------------------
c returns integration limits for numerical evaluation of particle
c and energy densities using quantum statistics
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/ciakt/gen,iafs,ians,genm

      igsp=int(gspecs(ians))

       if(mod(igsp,2).ne.0)then
      a=0.001
       else
      a=0.0
       endif

      b=0.0
      bb=(chemgc(ians)+tem*80.)**2-aspecs(ians)**2
      if(ish.ge.9)write(ifch,*)'bb:',bb
      if(bb.ge.0.0)b=sqrt(bb)
      if(bb.lt.0.0)then
      if(ish.ge.9)write(ifch,*)'failure at hgclim, bb=',bb
      if(ish.ge.9)write(ifch,'(1x,a,i5,a,2x,f12.6,1x,a,2x,f9.6)')
     *'mu(',ispecs(ians),'):',chemgc(ians),' T:',tem
      endif
      if(ish.ge.9)write(ifch,*)'ians:',ians,' a:',a,' b:',b
      return
      end

c------------------------------------------------------------------------
      subroutine hgcnbi(iret)
c-----------------------------------------------------------------------
c uses hgcaaa results to generate initial hadron set, nlattc, iozero
c input:
c    ptlngc(1:nspecs): particle number expectation values  /cgchg/
c output:
c     nump:           number of hadrons   /chnbin/
c     ihadro(1:nump): hadron ids          /chnbin/
c     nlattc:         lattice size        /clatt/
c     iozero:         zero weight         /metr1/
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/chnbin/nump,ihadro(maxp)
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cgctot/rmstot,ptltot
      common/camgc/amgc,samgc,amtot
      common/cflavs/nflavs,kef(nflav),chem(nflav)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)
      common/clatt/nlattc,npmax
      common/cgcnb/nptlgc(mspecs)
      common/ctaue/taue
      common/cgck/k(nflav),kp(nflav),kps(nflav)
     *,idp(maxp),ida(mspecs),idb(mspecs)
      integer hgcndn

      iret=0
      isho=ish
      if(ishsub/100.eq.50)ish=mod(ishsub,100)

      if(ish.ge.7)write(ifch,*)('-',l=1,10)
     *,' entry sr hgcnbi ',('-',l=1,30)


      nh=nint(ptltot)
      iug=(1+iospec)/2*2-1
      if(iug.lt.9)call utstop('hgcnbi: iospec < 9&')

c     determine nlattc
c     ----------------
        if(ionlat.eq.1)then
      s1=ptltot+2.*rmstot
      s2=1.3*ptltot
      s=max(s1,s2,6.)
      nlattc=nint(s)
       elseif(ionlat.eq.2)then
      s1=ptltot+3.*rmstot
      s2=1.5*ptltot
      s=max(s1,s2,6.)
      nlattc=nint(s)
       elseif(ionlat.eq.3)then
      s1=ptltot+4.*rmstot
      s2=2.*ptltot
      s=max(s1,s2,6.)
      nlattc=nint(s)
        elseif(ionlat.eq.0)then
      nlattc=8*(tecm/10)*(1/(tecm/volu))**0.2*(nspecs/3.)**0.3
      if(aspecs(1).lt.0.010)nlattc=nlattc*3
      nlattc=max(nlattc,20)
        endif

      if(ish.ge.7)write(ifch,*)'nlattc:',nlattc

c     determine iozero
c     ----------------
      if(iozero.eq.-1)then
      iozero=nspecs
      elseif(iozero.eq.-2)then
      iozero=nspecs*int(sqrt(volu/tecm))
      endif

c     modify iozero for testing
c     -------------------------
      if(iozevt.gt.0)then
      iozero=(nrevt/iozevt+1)*iozinc   !nrevt=event number - 1 !!
      write(ifch,*)'nrevt+1:',nrevt+1,'   iozero:',iozero
      endif

c     initial hadron set
c     ------------------
      ammin=2.*aspecs(1)
      if(tecm.lt.ammin)then
      write(ifch,*)'impossible to generate hadron configuration'
      call utstop('hgcnbi: tecm less than two pi0 masses&')
      endif

      kk=1
100   continue

       if(kk.gt.20)then
       iret=1
      if(ish.ge.7)then
      write(ifch,*)'failed to generate hadron set for'
     *,' event:',nrevt+1
      write(ifch,*)'u d s :',keu,ked,kes,' E:',tecm
      write(ifch,*)('-',i=1,30)
     *,' exit sr hgcnbi ',('-',i=1,10)
      endif
      ish=isho
      return
        endif

      amtot=0.0
      do i=1,nspecs
      nptlgc(i)=0
      enddo
      do ii=1,nflavs
      k(ii)=kef(ii)
      enddo

      if(ish.ge.7)write(ifch,*)
     *'sample hadron multiplicities and total mass:'

      kbar=keu+ked+kes
      kpar=iabs(keu)+iabs(ked)+iabs(kes)
      nbar=kbar/3.
      if(ish.ge.7)write(ifch,*)'baryon number:',nbar,' parton number:'
     *,kpar

      nn=2
      if(ioinco.ne.2)then
      nn=hgcndn(0)
      else
      nn=nh
      endif
      nb=iabs(nbar)
      if(ish.ge.7)write(ifch,*)'<n>:',nh,' n_sample:',nn,' n_bar:',nb
      if(nn.gt.nb.and.nb.ne.0.and.nb.ge.nh)nn=nb
      if(nn.lt.nb.and.nb.ne.0)nn=nb
      km=kpar-iabs(kbar)
      nt=km/2+nb
      if(nt.gt.nn)nn=nt
      nn=max(nn,2)

      if(ioinco.eq.2)then
      nit=15*taue
      else
      itpn=100
      nit=nn*itpn
      endif
      nbb=0
      n=0

c     start with nb protons
      nptlgc(19)=nptlgc(19)+nb
      n=nb
      amtot=amtot+nb*aspecs(19)
      do ii=1,nflavs
      k(ii)=k(ii)-ifok(ii,19)*nb
      enddo
      nbb=nbb+nb


       do it=1,nit

      xsp=nspecs
      x0=0.5
      xib=x0+xsp*rangen()
      ib=nint(xib)
      if(ib.gt.nspecs)ib=nspecs
      if(ib.lt.1)ib=1
      kb=ifok(1,ib)+ifok(2,ib)+ifok(3,ib)
      if(rangen().lt.0.5.and.nptlgc(ib).ge.1)then
      ni=-1
      else
      ni=1
      endif
      as=1.0
      if(nptlgc(ib).eq.0)as=0.5
      if(nptlgc(ib).eq.1.and.ni.eq.(-1))as=2.0
      if(ish.ge.9)write(ifch,*)
     *'id:',ispecs(ib),' <i>:',ptlngc(ib),' ni:',ni

         if(ni.ne.0)then

       if(ptlngc(ib).gt.5.0)then

      pnla=hgcpnl(ib,0)
      pnlb=hgcpnl(ib,ni)
      pnlog=-pnla+pnlb
      if(ish.ge.9)write(ifch,*)'pnlog:',pnlog
      if(pnlog.lt.60)then
      pn=exp(pnlog)
      else
      pn=1.1
      endif

       else

      if(ni.eq.1)then
      pn=ptlngc(ib)/(nptlgc(ib)+1)
      elseif(ni.eq.(-1).and.ptlngc(ib).gt.1.e-20)then
      pn=nptlgc(ib)/ptlngc(ib)
      elseif(nptlgc(ib).gt.0)then
      pn=1.1
      else
      pn=0.0
      endif

       endif

       pm=1.0
       if(ioinfl.ge.0)then
      pmla=hgcpml(ib,0,ib,0)
      pmlb=hgcpml(ib,ni,ib,0)
      pmlog=-pmla+pmlb
      if(ish.ge.9)write(ifch,*)'pmlog:',pmlog
      if(pmlog.lt.60)then
      pm=exp(pmlog)
      else
      pm=1.1
      endif
       endif

      p=pn*pm*as
      r=rangen()
      if(r.le.p)then
      nptlgc(ib)=nptlgc(ib)+ni
      n=n+ni
      amtot=amtot+ni*aspecs(ib)
      do ii=1,nflavs
      k(ii)=k(ii)-ifok(ii,ib)*ni
      enddo
      if(kb.ne.0)nbb=nbb+ni
      if(ish.ge.7.and.ni.gt.0)write(ifch,*)'add:'
      if(ish.ge.7.and.ni.lt.0)write(ifch,*)'remove:'
      if(ish.ge.7)write(ifch,*)'id:',ispecs(ib),' <n_i>:',ptlngc(ib)
     *,' n_i:',nptlgc(ib)
      if(ish.ge.7)write(ifch,*)'<n>:',nn,' it:',it
      if(ish.ge.7)write(ifch,*)'<M>:',amgc,' M:',amtot
      if(ish.ge.7)write(ifch,*)'p:',p,' r:',r
      if(ish.ge.7)write(ifch,*)'flav defect: u:',k(1),' d:'
     *,k(2),' s:',k(3)
      if(n.ge.nn.and.ioinco.ne.2)goto102
      endif

       endif

       enddo


102   continue

       ndd=0
c      if(nbb.lt.nb)then
c      nba=nb-nbb
c     if(nbar.gt.0)then
c     if(ish.ge.7)write(ifch,*)'add protons: nba:',nba
c     nptlgc(19)=nptlgc(19)+nba
c     n=n+nba
c     amtot=amtot+aspecs(19)*nba
c     elseif(nbar.lt.0)then
c     if(ish.ge.7)write(ifch,*)'add aprotons: nba:',nba
c     nptlgc(20)=nptlgc(20)+nba
c     n=n+nba
c     amtot=amtot+aspecs(20)*nba
c     endif
c      endif
       if(n.lt.nn.and.ioinco.ne.2)then
      ndd=nn-n
      nd=mod(ndd,4)
      xn=n
      xnn=nn
      xl=(xnn-xn)/4.
      l=aint(xl)
      if(ish.ge.7)write(ifch,*)'add pions/etas: ndd:',ndd
     *,' l:',l,' nd:',nd
        if(l.ge.1)then
       do j=1,l
      nptlgc(1)=nptlgc(1)+1
      nptlgc(2)=nptlgc(2)+1
      nptlgc(3)=nptlgc(3)+1
      nptlgc(8)=nptlgc(8)+1
      amtot=amtot+aspecs(1)+aspecs(2)+aspecs(3)+aspecs(8)
       enddo
        endif
      if(nd.eq.1)then
      nptlgc(1)=nptlgc(1)+1
      amtot=amtot+aspecs(1)
      elseif(nd.eq.2)then
      nptlgc(2)=nptlgc(2)+1
      nptlgc(3)=nptlgc(3)+1
      amtot=amtot+aspecs(2)+aspecs(3)
      elseif(nd.eq.3)then
      nptlgc(2)=nptlgc(2)+1
      nptlgc(3)=nptlgc(3)+1
      nptlgc(1)=nptlgc(1)+1
      amtot=amtot+aspecs(2)+aspecs(3)+aspecs(1)
      endif
       endif

       if(n.eq.0.and.ioinco.eq.2)then
      nptlgc(2)=nptlgc(2)+1
      nptlgc(3)=nptlgc(3)+1
      amtot=amtot+aspecs(2)+aspecs(3)
       elseif(n.eq.1.and.ioinco.eq.2)then
      nptlgc(1)=nptlgc(1)+1
      amtot=amtot+aspecs(1)
       endif

      if(amtot.ge.tecm.and.ioinfl.ge.0)then
      if(ish.ge.7)write(ifch,*)
     *'total mass exceeded , redo configuration'
      kk=kk+1
      goto100
      endif


      iii=0
      if(ish.ge.7)then
        write(ifch,*)'u d s :',keu,ked,kes,' E:',tecm
        write(ifch,*)
     *'hadron set without flavor conservation:'
      endif
      do i=1,nspecs
      n=nptlgc(i)
      if(n.ge.1)then
      do j=1,n
      iii=iii+1
      if(iii.gt.maxp)stop'iii>maxp in hgcnbi'
      idp(iii)=ispecs(i)
      enddo
      endif
      enddo
      if(ish.ge.7)then
        write(ifch,'(1x,10i6)')(idp(i),i=1,iii)
        write(ifch,*)'flav defect: u:',k(1),' d:'
     *,k(2),' s:',k(3)
        write(ifch,*)'M:',amtot,' <M>:',amgc
      endif
      if(ioinfl.le.0)goto1000

      ll=1
      llmax=nn*25
      ior=1

120        if(k(1).ne.0.or.k(2).ne.0.or.k(3).ne.0)then

        if(kk.gt.6)ior=0

      if(ish.ge.7)write(ifch,*)
     *'remaining flavor defect before operation:',ll
      if(ish.ge.7)write(ifch,*)'flav defect: u:',k(1),' d:'
     *,k(2),' s:',k(3)

      nida=0
      do i=1,nspecs
      if(nptlgc(i).gt.0)then
      nida=nida+1
      ida(nida)=i
      endif
      enddo

      if(nida.eq.0)then
      if(ish.ge.7)write(ifch,*)'no proposals in a , redo'
      kk=kk+1
      goto100
      endif


      xna=0.5+nida*rangen()
      na=nint(xna)
      if(na.gt.nida)na=nida
      if(na.lt.1)na=1
      ia=ida(na)
      if(ish.ge.7)write(ifch,*)'nida:',nida,' ia:',ia

      nidb=0
      do ii=1,nflavs
      kp(ii)=k(ii)+ifok(ii,ia)
      kps(ii)=isign(1,kp(ii))
      enddo
      if(ish.ge.7)write(ifch,*)
     *'   assemble: u:',kp(1),' d:',kp(2),' s:',kp(3)
      do i=1,nspecs
      iacc=0
      naccsp=0
      naccmi=1
      do ii=1,nflavs
      naccsp=naccsp+iabs(ifok(ii,i))
      if(kp(ii).ne.0)then
      if(kps(ii)*ifok(ii,i).le.kps(ii)*kp(ii)
     *.and.kps(ii)*ifok(ii,i).gt.0)iacc=iacc+iabs(ifok(ii,i))
      endif
      enddo
      if(kp(1).eq.0.and.kp(2).eq.0.and.kp(3).eq.0)naccmi=0
      if(iacc.eq.naccsp.and.naccsp.ge.naccmi)then
      nidb=nidb+1
      idb(nidb)=i
      endif
      enddo

      if(nidb.eq.0)then
      if(ish.ge.7)write(ifch,*)'no proposals in b , redo'
      kk=kk+1
      goto100
      endif

      xnb=0.5+nidb*rangen()
      nb=nint(xnb)
      if(nb.gt.nidb)nb=nidb
      if(nb.lt.1)nb=1
      ib=idb(nb)
      if(ish.ge.7)write(ifch,*)'nidb:',nidb,' ib:',ib
      if(ish.ge.7)write(ifch,*)
     *'proposal:',ispecs(ia),' --> ',ispecs(ib)

      asym=1.0

c      if(asym.gt.0.0)then

       if(ptlngc(ia).gt.5.0)then
      pnali=hgcpnl(ia,0)
      pnalf=hgcpnl(ia,-1)
      pnalog=-pnali+pnalf
      if(ish.ge.7)write(ifch,*)'pnalog:',pnalog
      if(pnalog.lt.60)then
      pna=exp(pnalog)
      else
      pna=1.1
      endif
       else
      if(ptlngc(ia).gt.1.e-20)then
      pna=nptlgc(ia)/ptlngc(ia)
      elseif(nptlgc(ia).gt.0)then
      pna=1.1
      else
      pna=0.0
      endif
       endif

       if(ptlngc(ib).gt.5.0)then
      pnbli=hgcpnl(ib,0)
      pnblf=hgcpnl(ib,1)
      pnblog=-pnbli+pnblf
      if(ish.ge.7)write(ifch,*)'pnblog:',pnblog
      if(pnblog.lt.60)then
      pnb=exp(pnblog)
      else
      pnb=1.1
      endif
       else
      pnb=ptlngc(ib)/(nptlgc(ib)+1)
       endif


      pmli=hgcpml(ia,0,ib,0)
      pmlf=hgcpml(ia,-1,ib,1)
      pmlog=-pmli+pmlf
      if(ish.ge.7)write(ifch,*)'pmlog:',pmlog
      if(pmlog.lt.60)then
      pm=exp(pmlog)
      else
      pm=1.1
      endif

      p=pna*pnb*pm*asym
      if(ior.eq.0)then
      r=0.0
      else
      r=rangen()
      endif

c      else

c     r=1.0
c     p=0.0

c      endif

       if(r.lt.p)then
      if(ish.ge.7)write(ifch,*)'p:',p,' r:',r,' asymmetry:',asym
      if(ish.ge.7)write(ifch,*)'remove ',ispecs(ia),'  add ',ispecs(ib)
     *,'  proposal accepted'
      nptlgc(ia)=nptlgc(ia)-1
      nptlgc(ib)=nptlgc(ib)+1
      amtot=amtot-aspecs(ia)+aspecs(ib)
      do ii=1,nflavs
      k(ii)=k(ii)+ifok(ii,ia)-ifok(ii,ib)
      enddo
       endif

       
        if(k(1).ne.0.or.k(2).ne.0.or.k(3).ne.0)then
       ll=ll+1
      if(ll.le.llmax)then
      goto120
      else
      if(ish.ge.7)write(ifch,*)'failed to remove defect, redo'
      kk=kk+1
      goto100
       endif
        endif

         endif

1000  continue

      nump=0
      kcu=0
      kcd=0
      kcs=0
      do i=1,nspecs
      n=nptlgc(i)
      if(n.ge.1)then
      do j=1,n
      nump=nump+1
      ihadro(nump)=ispecs(i)
      kcu=kcu+ifok(1,i)
      kcd=kcd+ifok(2,i)
      kcs=kcs+ifok(3,i)
      enddo
      endif
      enddo

          if(ioinfl.gt.0)then
        if(kcu.ne.keu.or.kcd.ne.ked.or.kcs.ne.kes)then
      if(ish.ge.7)write(ifch,*)
     *'failed to remove flavor defect, redo configuration'
      kk=kk+1
      goto100
        endif
          endif

          nutot=0
          chitot=0.
      if(ioinct.ge.1)then
        chitot=0.0
        nutot=nspecs
        do i=1,nspecs
        chi=0.0
        if(rmsngc(i).gt.1.e-10)chi=(ptlngc(i)-nptlgc(i))/rmsngc(i)
        chitot=chitot+chi**2
        enddo
        call xhgccc(chitot)

        u=0
        d=0
        s=0
        do i=1,nspecs
        u=u+ifok(1,i)*nptlgc(i)
        d=d+ifok(2,i)*nptlgc(i)
        s=s+ifok(3,i)*nptlgc(i)
        enddo
        call xhgcfl(u,d,s,0)
        call xhgcam(amtot,0)
      endif

      if(ish.ge.7)then
        write(ifch,*)
     *'initial hadron set for droplet decay:'
        write(ifch,'(1x,10i6)')(ihadro(i),i=1,nump)
      endif
       if(nump.ge.nlattc)then
         nlattc=nump+1
         if(ish.ge.7)then
           write(ifch,*)'initial set > nlattc !'
           write(ifch,*)'new nlattc:',nlattc
         endif
       endif
       if(ish.ge.7)then
         write(ifch,*)'keu:',kef(1),' kcu:',kcu,' ku:',k(1)
         write(ifch,*)'ked:',kef(2),' kcd:',kcd,' kd:',k(2)
         write(ifch,*)'kes:',kef(3),' kcs:',kcs,' ks:',k(3)
         write(ifch,*)' nh:',nh,' nump:',nump
         write(ifch,*)' nu:',nutot,'  chi^2:',chitot
         write(ifch,*)'iozero:',iozero,'  iomom:',iomom
         write(ifch,*)
     *'total mass:',amtot,' droplet mass:',tecm
         write(ifch,*)'trials needed:',kk
     *,' operations needed:',ll
         write(ifch,*)'iterations:',it,' pions added:',ndd
         write(ifch,*)('-',i=1,30)
     *,' exit sr hgcnbi ',('-',i=1,10)
       endif
      ish=isho
      return

      end

c--------------------------------------------------------------------
      integer function hgcndn(i)
c--------------------------------------------------------------------
c returns random multiplicity from gaussian distribution for species i
c---------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cgctot/rmstot,ptltot
      common/clatt/nlattc,npmax
      a=iowidn
      kk=0

       if(i.eq.0)then

1     continue
      kk=kk+1
      p=0.0
      nmin=2
      nh=nint(ptltot)
      nmax=nlattc
      xn=1.5+(nmax-nmin)*rangen()
      n=nint(xn)
      x=(n-ptltot)**2/2.0
      y=-70.
      if(rmstot.gt.1.e-15)y=-x/rmstot**2*a**2
      if(y.lt.70.)p=exp(y)
      if(rmstot.gt.1.e-15.and.iowidn.lt.0)p=p/sqrt(2.*pi)/rmstot
      if(p.ge.rangen())then
      hgcndn=n
      if(ish.ge.9)write(ifch,*)'hgcndn: k:',kk,' n:',hgcndn
      return
      else
      if(kk.le.25)goto1
      hgcndn=max(2,nh)
      if(ish.ge.9)write(ifch,*)'hgcndn: k:',kk,' n:',hgcndn
      return
      endif

       else

2     continue
      kk=kk+1
      p=0.0
      nmin=0
      nh=nint(ptlngc(i))
      nmax=2*nh
      nmax=max(2,nmax)
      xn=-0.5+(nmax-nmin)*rangen()
      n=nint(xn)
      x=(n-ptlngc(i))**2/2.0
      if(x.lt.1.e-30)then
      p=1.
      else
      y=-70.
      if(rmsngc(i).gt.1.e-15)y=-x/rmsngc(i)**2
      if(y.lt.70.)p=exp(y)
      if(rmsngc(i).gt.1.e-15.and.iowidn.lt.0)
     *p=p/sqrt(2.*pi)/rmsngc(i)
      endif
      if(p.ge.rangen())then
      hgcndn=n
      if(ish.ge.9)write(ifch,*)'hgcndn: k:',kk,' n:',hgcndn
      return
      else
      if(kk.le.25)goto2
      hgcndn=nh
      if(ish.ge.9)write(ifch,*)'hgcndn: k:',kk,' n:',hgcndn
      return
      endif

       endif

      end

c--------------------------------------------------------------------
      function hgcpml(i1,n1,i2,n2)
c--------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/camgc/amgc,samgc,amtot
      common/cgcnb/nptlgc(mspecs)
      if(ish.ge.9)write(ifch,*)'i1:',i1,' i2:',i2
      if(ish.ge.9)write(ifch,*)'n1:',n1,' n2:',n2
      hgcpml=-1.e30
      ampr=n1*aspecs(i1)+n2*aspecs(i2)
      if((amtot+ampr).lt.tecm.and.(amtot+ampr).ge.0
     *.and.nptlgc(i1).ge.(-n1).and.nptlgc(i2).ge.(-n2))then
      hgcpml=0.0
      pl=(amtot-amgc+ampr)**2/2.0
      if(pl.lt.1.e-30)then
      hgcpml=0.0
      return
      endif
      if(samgc.gt.1.e-15)hgcpml=-pl/samgc**2
      endif
      if(ish.ge.9)write(ifch,*)'hgcpml:',hgcpml
      return
      end

c--------------------------------------------------------------------
      function hgcpnl(i,n)
c--------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cgcnb/nptlgc(mspecs)
      if(ish.ge.9)write(ifch,*)'i:',i,' n:',n
      hgcpnl=-1.e30
      if(nptlgc(i).ge.(-n))then
      pl=(nptlgc(i)-ptlngc(i)+n)**2/2.0
      if(pl.lt.1.e-30)then
      hgcpnl=0.0
      return
      endif
      if(rmsngc(i).gt.1.e-15)hgcpnl=-pl/rmsngc(i)**2
      endif
      if(ish.ge.9)write(ifch,*)'hgcpnl:',hgcpnl
      return
      end


c--------------------------------------------------------------------
      subroutine hgcpen
c--------------------------------------------------------------------
c returns array for twodimensional plot of energy- and flavor-
c density
c--------------------------------------------------------------------
c xpar1,xpar2 temperature range
c xpar3       # of bins for temperature
c xpar4,xpar5 chem.pot. range
c xpar6       # of bins for chem.pot.
c xpar7       max. density
c xpar8       strange chem.pot.
c--------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cflavs/nflavs,kef(nflav),chem(nflav)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)
      common/ciakt/gen,iafs,ians,genm
      parameter (nbin=100)
      real edensi(nbin,nbin),qdensi(nbin,nbin)
      external hgcfhe
      external hgcfhn
      external hgcfbe
      external hgcfbn

      iug=(1+iospec)/2*2-1

c     initialization
c     --------------

      if(iug.eq.1)nflavs=1
      if(iug.eq.3)nflavs=2
      if(iug.eq.5)nflavs=2
      if(iug.eq.7)nflavs=3
      if(iug.eq.9)nflavs=3
      if(iug.eq.11)nflavs=3
      tem=0.0
      do i=1,nflavs
      chem(i)=0.0
      enddo
      call hgchac(0)
      do i=1,nspecs
      ptlngc(i)=0.0
      rmsngc(i)=0.0
      enddo

      nbt=nint(xpar3)
      nbc=nint(xpar6)
      nbc=min(nbc,100)
      nbt=min(nbt,100)
      dt=(xpar2-xpar1)/nbt
      dc=(xpar5-xpar4)/nbc
      ymax=xpar7
      cs=xpar8


      t0=xpar1+dt/2.
      c0=xpar4+dc/2
      do i=1,nbc
      chem(1)=c0+(i-1)*dc
      chem(2)=chem(1)
      chem(3)=cs
      chem(4)=0.0
      chem(5)=0.0
      chem(6)=0.0
      call hgchac(0)
      do ii=1,nbt
      tem=t0+(ii-1)*dt
      if(ish.ge.5)write(ifch,*)' mu:',chem(1),' T:',tem

       qd=0.0
       ed=0.0

      do ians=1,nspecs

      call hgclim(a,b)

      if(b.eq.0.0)then
      hden=0.0
      elseif(iostat.eq.0)then
      call uttraq(hgcfhn,a,b,hden)
      elseif(iostat.eq.1)then
      call uttraq(hgcfbn,a,b,hden)
      endif
      hd=hden*gspecs(ians)/2./pi**2/hquer**3

      if(ish.ge.7)write(ifch,*)'i:',ians,' n_u:',ifok(1,ians),' hd:',hd

      qd=qd+ifok(1,ians)*hd+ifok(2,ians)*hd
      if(qd.gt.ymax)qd=ymax
c     if(qd.gt.ymax)qd=0.0
      if(qd.lt.-ymax)qd=-ymax
c     if(qd.lt.-ymax)qd=0.0


      if(b.eq.0.0)then
      edi=0.0
      elseif(iostat.eq.0)then
      call uttraq(hgcfhe,a,b,edi)
      elseif(iostat.eq.1)then
      call uttraq(hgcfbe,a,b,edi)
      endif
      edi=edi*gspecs(ians)/2./pi**2/hquer**3

      if(ish.ge.7)write(ifch,*)'i:',ians,' mu:',chemgc(ians)
     *                        ,' edi:',edi

      ed=ed+edi
      if(ed.gt.ymax)ed=ymax
c     if(ed.gt.ymax)ed=0.0
      enddo

      if(ish.ge.5)write(ifch,*)' ed:',ed,' qd:',qd
      edensi(i,ii)=ed
      qdensi(i,ii)=qd

      enddo
      enddo

      write(ifhi,'(a)')      'openhisto'
      write(ifhi,'(a,2e11.3)')'xrange',xpar1,xpar2
      write(ifhi,'(a,2e11.3)')'yrange',xpar4,xpar5
      write(ifhi,'(a)')      'set ityp2d 5'
      write(ifhi,'(a,i4)')   'array2d',nbt
      do j=1,nbc
      do jj=1,nbt
      write(ifhi,'(e11.3)') edensi(j,jj)
      enddo
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot2d'

      write(ifhi,'(a)')      'openhisto'
      write(ifhi,'(a,2e11.3)')'xrange',xpar1,xpar2
      write(ifhi,'(a,2e11.3)')'yrange',xpar4,xpar5
      write(ifhi,'(a)')      'set ityp2d 5'
      write(ifhi,'(a,i4)')   'array2d',nbt
      do j=1,nbc
      do jj=1,nbt
      write(ifhi,'(e11.3)') qdensi(j,jj)
      enddo
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot2d'

        return
        end

c--------------------------------------------------------------------
      subroutine hgcpfl
c--------------------------------------------------------------------
c returns array for twodimensional plot of energy- and flavor-
c density fluctuations
c--------------------------------------------------------------------
c xpar1,xpar2 temperature range
c xpar3       # of bins for temperature
c xpar4,xpar5 chem.pot. range
c xpar6       # of bins for chem.pot.
c xpar7       max. density
c xpar8       strange chem.pot.
c--------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cflavs/nflavs,kef(nflav),chem(nflav)
      common/ciakt/gen,iafs,ians,genm
      parameter (nbin=100)
      real efl(nbin,nbin),qfl(nbin,nbin),v(nbin),wn(nbin),we(nbin)
      external hgcfhf
      external hgcfhe
      external hgcfhn
      external hgcfhw
      external hgcfbf
      external hgcfbe
      external hgcfbn

      iug=(1+iospec)/2*2-1

c     initialization
c     --------------

      if(iug.eq.1)nflavs=1
      if(iug.eq.3)nflavs=2
      if(iug.eq.5)nflavs=2
      if(iug.eq.7)nflavs=3
      if(iug.eq.9)nflavs=3
      if(iug.eq.11)nflavs=3
      tem=0.0
      do i=1,nflavs
      chem(i)=0.0
      enddo
      call hgchac(0)
      do i=1,nspecs
      ptlngc(i)=0.0
      rmsngc(i)=0.0
      enddo

      nbt=nint(xpar3)
      nbv=nint(xpar6)
      nbv=min(nbv,100)
      nbt=min(nbt,100)
      dt=(xpar2-xpar1)/nbt
      dv=(xpar5-xpar4)/nbv
      ymax=1.e20
      chem(1)=xpar7
      chem(2)=xpar7
      chem(3)=xpar8
      call hgchac(0)


      t0=xpar1+dt/2.
      v0=xpar4
      do i=1,nbv
      volu=v0+(i-1)*dv
      do ii=1,nbt
      tem=t0+(ii-1)*dt
      if(ish.ge.5)write(ifch,*)'volu:',volu,' tem:',tem

       ev=0.0
       ee=0.0
       qv=0.0
       qe=0.0

      do ians=1,nspecs

      call hgclim(a,b)

      if(b.eq.0.0)then
      hn=0.0
      hv=0.0
      elseif(iostat.eq.0)then
      call uttraq(hgcfhn,a,b,hn)
      call uttraq(hgcfhw,a,b,hv)
      elseif(iostat.eq.1)then
      call uttraq(hgcfbn,a,b,hn)
      hv=hn
      endif
      hn=hn*volu*gspecs(ians)/2./pi**2/hquer**3
      hv=hv*volu*gspecs(ians)/2./pi**2/hquer**3
      if(ish.ge.5)write(ifch,*)'hn:',hn,' hv:',hv

      hn=max(hn,1.e-15)
      qv=qv+hv
      qe=qe+hn


      if(qv.gt.ymax)qv=ymax
      if(qe.gt.ymax)qe=ymax


      if(b.eq.0.0)then
      eei=0.0
      evi=0.0
      elseif(iostat.eq.0)then
      call uttraq(hgcfhe,a,b,eei)
      call uttraq(hgcfhf,a,b,evi)
      elseif(iostat.eq.1)then
      call uttraq(hgcfbe,a,b,eei)
      call uttraq(hgcfbf,a,b,evi)
      endif
      eei=eei*volu*gspecs(ians)/2./pi**2/hquer**3
      evi=evi*volu*gspecs(ians)/2./pi**2/hquer**3
      if(ish.ge.5)write(ifch,*)'eei:',eei,' evi:',evi


      eei=max(eei,1.e-15)
      ev=ev+evi
      ee=ee+eei
      if(ev.gt.ymax)ev=ymax
      if(ee.gt.ymax)ee=ymax
      enddo
      if(ish.ge.5)write(ifch,*)'qv:',qv,' ev:',ev

      qfl(i,ii)=0.
      efl(i,ii)=0.
      if(ev.gt.0.0.and.ee.gt.1.e-15)efl(i,ii)=sqrt(ev)/ee
      if(qv.gt.0.0.and.ee.gt.1.e-15)qfl(i,ii)=sqrt(qv)/qe
      if(tem.eq.0.195)then
      we(i)=efl(i,ii)
      wn(i)=qfl(i,ii)
      v(i)=volu
      endif

      enddo
      enddo

      write(ifhi,'(a)')      'openhisto'
      write(ifhi,'(a,2e11.3)')'xrange',xpar1,xpar2
      write(ifhi,'(a,2e11.3)')'yrange',xpar4,xpar5
      write(ifhi,'(a)')      'set ityp2d 5'
      write(ifhi,'(a,i4)')   'array2d',nbt
      do j=1,nbv
      do jj=1,nbt
      write(ifhi,'(e11.3)') efl(j,jj)
      enddo
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot2d'

      write(ifhi,'(a)')      'openhisto'
      write(ifhi,'(a,2e11.3)')'xrange',xpar1,xpar2
      write(ifhi,'(a,2e11.3)')'yrange',xpar4,xpar5
      write(ifhi,'(a)')      'set ityp2d 5'
      write(ifhi,'(a,i4)')   'array2d',nbt
      do j=1,nbv
      do jj=1,nbt
      write(ifhi,'(e11.3)') qfl(j,jj)
      enddo
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot2d'

      write(ifhi,'(a)')      'newpage zone 1 2 1'
      write(ifhi,'(a)')      'openhisto'
      write(ifhi,'(a,2e11.3)')'xrange',xpar4,xpar5
      write(ifhi,'(a)')      'htyp lfu xmod lin ymod lin'
      write(ifhi,'(a,i4)')   'array 2'
      do j=1,nbv
      write(ifhi,'(2e13.5)')v(j),we(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0'

      write(ifhi,'(a)')      'openhisto'
      write(ifhi,'(a,2e11.3)')'xrange',xpar4,xpar5
      write(ifhi,'(a)')      'htyp lfu xmod lin ymod lin'
      write(ifhi,'(a,i4)')   'array 2'
      do j=1,nbv
      write(ifhi,'(2e13.5)')v(j),wn(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0'


        return
        end


c------------------------------------------------------------------
      subroutine hgcpyi(ist)
c------------------------------------------------------------------
c returns particle yield
c input:
c   tem   : temperature
c   chemgc: chemical potentials
c output:
c   ptlngc: expectation value of particle number for each species
c   rmsngc: standard deviation of ptlngc
c   ptltot: total particle number
c   rmstot: standard deviation of ptltot
c works for hadrons and partons
c  ist=1 boltzmann statistics
c  ist=0 quantum statistics
c--------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cgctot/rmstot,ptltot
      common/camgc/amgc,samgc,amtot
      common/ciakt/gen,iafs,ians,genm
      external hgcfhw
      external hgcfhn

       if(iabs(ispecs(nspecs)).lt.10)then

c     parton yield
c     ------------
      if(ish.ge.5)write(ifch,*)'parton yield:'
      gln=16.*1.20206*tem**3/pi**2*volu/hquer**3
      sdg=sqrt(gln)   !!???
      if(ish.ge.5)write(ifch,'(1x,a,f10.4,2x,a,f9.4,a)')
     *'<N(    0)> :',gln,' sigma :',sdg,' (qm-statistics!)'
      ptltot=gln
      rmstot=0.0
      vartot=gln

       else

      if(ish.ge.5)write(ifch,*)'hadronic yield:'
      ptltot=0.0
      rmstot=0.0
      vartot=0.0

       endif

      amgc=0.0
      samgc=0.0

       do ians=1,nspecs

c     hadronic yield
c     --------------
       if(ist.eq.0)then

      call hgclim(a,b)
      if(b.eq.0.0)then
      hden=0.0
      else
      call uttraq(hgcfhn,a,b,hden)
      endif
      ptlngc(ians)=hden*volu*gspecs(ians)/2./pi**2/hquer**3

       else

       if((chemgc(ians)/tem).gt.70.)then
      hpd=1.e30
       else
      hpd=exp(chemgc(ians)/tem)
       endif
       if(aspecs(ians).ne.0.)then
      fk2=hgcbk(2,aspecs(ians)/tem)
      hpd=hpd*gspecs(ians)*aspecs(ians)**2*tem*fk2
     */2./pi**2/hquer**3
       else
      hpd=hpd*gspecs(ians)*tem**3/pi**2/hquer**3
       endif
      ptlngc(ians)=hpd*volu

       endif

      ptltot=ptltot+ptlngc(ians)
      amgc=amgc+ptlngc(ians)*aspecs(ians)
      if(amgc.ge.tecm)amgc=tecm*0.9

c     standard deviation
c     ------------------
      rmsngc(ians)=0.0

       if(ist.eq.0)then

      call uttraq(hgcfhw,a,b,var)
      var=var*gspecs(ians)*volu/2./pi**2/hquer**3
      vartot=vartot+var
      if(var.ge.0.0)rmsngc(ians)=sqrt(var)
      samgc=samgc+var*aspecs(ians)

       else

      if(ptlngc(ians).ge.0.0)rmsngc(ians)=sqrt(ptlngc(ians))
      vartot=vartot+ptlngc(ians)
      samgc=samgc+ptlngc(ians)*aspecs(ians)

       endif


      if(ish.ge.7)write(ifch,'(2x,a,i5,a,2x,f8.4,5x,a,3x,f8.4)')
     *'m(',ispecs(ians),')  :',aspecs(ians),'mu :',chemgc(ians)
      if(ish.ge.5)write(ifch,'(1x,a,i5,a,2x,f8.4,2x,a,2x,f10.4)')
     *'<N(',ispecs(ians),')> :',ptlngc(ians),'sigma :',rmsngc(ians)

       enddo

      if(vartot.ge.0.0)rmstot=sqrt(vartot)
      if(samgc.ge.0.0)samgc=sqrt(samgc)
      if(amgc.ge.tecm)samgc=sqrt(amgc)
      if(ish.ge.5)write(ifch,'(1x,a,2x,f8.4,2x,a,2x,f10.4)')
     *'<N(  all)> :',ptltot,'sigma :',rmstot
      if(ish.ge.5)write(ifch,'(1x,a,2x,f8.4,2x,a,2x,f10.4)')
     *'<M_tot>    :',amgc,'sigma :',samgc

      return
      end

c------------------------------------------------------------------------
      subroutine hgctbo(iba)
c------------------------------------------------------------------------
c returns new tem using boltzmann statistics in analytic form
c  input:
c    chemgc
c    tecm/volu
c  output:
c    tem
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/ciakt/gen,iafs,ians,genm
      external hgcbk
      external hgcbk1
      iba=0
      k=1
      t1=0.0
      t2=1.0

      goto15

10    tem=t1+.5*(t2-t1)
      if(tem.le.1.e-7)return
15    eden=0.0

        do i=1,nspecs

       if(aspecs(i).ne.0)then
      if(tem.ne.0.)arr=aspecs(i)/tem
      cba=(aspecs(i)/tem+12.*tem/aspecs(i)-3.*chemgc(i)/aspecs(i))
     **hgcbk(2,arr)+(3.-chemgc(i)/tem)*hgcbk1(arr)
       else
      cba=4.*tem-chemgc(i)
       endif

      if(cba.lt.0.0)then
      iba=1
      return
      endif

      x=71.
      if(tem.ne.0.)x=chemgc(i)/tem

       if(x.le.70.)then
      y=exp(x)
       else
      y=1.e30
       endif

       if(aspecs(i).ne.0.)then
      edi=y*(3./arr*hgcbk(2,arr)+hgcbk1(arr))
     **gspecs(i)*aspecs(i)**3*tem/2./pi**2/hquer**3
       else
      edi=y*3.*gspecs(i)*tem**4/pi**2/hquer**3
       endif

      eden=eden+edi

        enddo

      if(iabs(ispecs(nspecs)).lt.10)
     *eden=eden+(8.*pi**2*tem**4/15.+bag4rt**4)/hquer**3

      de=abs(eden-(tecm/volu))
      if(de.le.gen*(tecm/volu).or.de.le.genm)return
c     if(eden.ge.100.)return

       if(eden.gt.(tecm/volu))then
      t2=tem
      else
      t1=tem
       endif

       if(k.gt.300)return

      k=k+1
      goto10
      end

c----------------------------------------------------------------------
      subroutine hgctex
c----------------------------------------------------------------------
c returns new tem using massive quantum statistics in integral form
c  input:
c    chemgc
c    tecm/volu
c  output:
c    tem
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/ciakt/gen,iafs,ians,genm
      external hgcfhe
      k=1
      t1=0.0
      t2=tem+0.1
      goto15

c     new temperature
c     ---------------
10    tem=t1+.5*(t2-t1)
15    continue
      if(tem.le.1.e-6)return
      eden=0.0

       do ians=1,nspecs
      call hgclim(a,b)
      if(b.eq.0.0)then
      edi=0.0
      else
      call uttraq(hgcfhe,a,b,edi)
      endif
      edi=edi*gspecs(ians)/2./pi**2/hquer**3
      eden=eden+edi
       enddo

      if(iabs(ispecs(nspecs)).lt.10)
     *eden=eden+(8.*pi**2*tem**4/15.+bag4rt**4)/hquer**3

      de=abs(eden-(tecm/volu))
      if(de.le.gen*(tecm/volu).or.de.le.genm)return

       if(eden.gt.(tecm/volu))then
      t2=tem
      else
      t1=tem
       endif

       if(k.gt.300)then
       if(ish.ge.5)
     *write(ifch,*)'failure in tex'
      return
       endif

      k=k+1
      goto10
      end

c-----------------------------------------------------------------
      subroutine hgctm0
c-----------------------------------------------------------------
c returns new tem using massless quantum statistics in analytic form
c  input:
c    chemgc
c    tecm/volu
c  output:
c    tem
c----------------------------------------------------------------------

      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/ciakt/gen,iafs,ians,genm

      k=1

      t1=0.0
      t2=1.0
10    tem=t1+.5*(t2-t1)
      if(tem.le.1.e-6)return
      eden=0.0

        do i=1,nspecs

      igsp=int(gspecs(i))
      if(mod(igsp,2).eq.0)then
      edhm0=7./240.*pi**2*tem**4+chemgc(i)**2*tem**2/8.
     *+chemgc(i)**4/pi**2/16.
      else
      edhm0=pi**2*tem**4/30.+chemgc(i)**2*tem**2/4.
     *-chemgc(i)**4/pi**2/16.
      endif
      edi=edhm0*gspecs(i)/hquer**3


      eden=eden+edi
        enddo

      if(iabs(ispecs(nspecs)).lt.10)
     *eden=eden+(8.*pi**2*tem**4/15.+bag4rt**4)/hquer**3

      de=abs(eden-(tecm/volu))
      if(de.le.gen*(tecm/volu).or.de.le.genm)return

       if(eden.gt.(tecm/volu))then
      t2=tem
      else
      t1=tem
       endif

       if(k.gt.300)then
       if(ish.ge.5)
     *write(ifch,*)'failure in tm0'
      return
       endif

      k=k+1
      goto10
      end








c-----------------------------------------------------------------------
      subroutine xhgcam(amt,iii)
c-----------------------------------------------------------------------
c creates unnormalized histogram for total mass of grand
c canonically generated sample
c xpar1: nr. of bins
c xpar2: m_1 (lower boundary)
c xpar3: m_2 (upper boundary)
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter(nbmx=200)
      common/camdat/data(nbmx),datb(nbmx)
      parameter(mxclu=10000)
      real am(mxclu)
      character cen*6,cvol*6
      save am

      if(iii.eq.0)then

      am(nrclu)=amt

      return

      elseif(iii.lt.0)then

      nbin=nint(xpar3)
      x1=xpar1
      x2=xpar2
      dam=(x2-x1)/nbin
      write(cen,'(f6.1)')tecm
      write(cvol,'(f6.1)')volu

      do i=1,nbin
      data(i)=x1+(i-1)*dam
      datb(i)=0.0
      enddo

      do i=1,nrclu
      xnb=(am(i)-x1)/dam+1.
      nb=nint(xnb)
      if(nb.le.nbin.and.nb.ge.1)datb(nb)=datb(nb)+1
      enddo

      write(ifhi,'(a)')       'newpage zone 1 2 1'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp his'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis total mass"'
      write(ifhi,'(a)')    'text 0 0 "yaxis N"'
      write(ifhi,'(a,a)')     'text 0.4 0.91 "V='//cvol//'"'
      write(ifhi,'(a,a)')     'text 0.15 0.91 "E='//cen//'"'
      write(ifhi,'(a)')       'array 2'

         do j=1,nbin
      write(ifhi,'(2e13.5)')data(j),datb(j)
         enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'


      return

           endif

       end

c-----------------------------------------------------------------------
      subroutine xhgccc(chi)
c-----------------------------------------------------------------------
c creates unnormalized histogram for chi-squared test of initial
c configuration (grand-canonical results are used)
c for chi>0: chi-squared for each droplet configuration is written
c            to /cchi/
c for chi<0: creates histogram
c            xpar1 specifies lower limit
c            xpar2 specifies upper limit
c            xpar3 specifies bin width
c  newpage, zone and plot commands not included !!!
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter(nbin=200)
      common/chidat/data(nbin),datb(nbin)
      parameter(mxclu=10000)
      common/cchi/chi2(mxclu)
      character cnu*2,cinco*1,cen*6,cvol*6
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)

         if(chi.ge.0.0)then

      nrclu=nrclu+1
      chi2(nrclu)=chi

      return

         elseif(chi.lt.0.0)then

      x1=nint(xpar1)
      x2=nint(xpar2)
      da=xpar3
      write(cnu,'(i2)')nspecs
      write(cinco,'(i1)')ioinco
      write(cen,'(f6.1)')tecm
      write(cvol,'(f6.1)')volu

      if(x2.eq.0)x2=50.0
      da=max(0.1,da)
      a0=x1

      do i=1,nbin
      data(i)=a0+(i-1)*da
      datb(i)=0.0
      enddo

      do i=1,nrclu
      nb=(chi2(i)+da/2.-a0)/da
      if(nb.le.nbin.and.nb.ge.1)datb(nb)=datb(nb)+1
      enddo

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp his'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis [V]^2"'
      write(ifhi,'(a)')    'text 0 0 "yaxis f([V]^2,n?eff!)"'
      if(iappl.eq.4)write(ifhi,'(a,a)')'text 0.4 0.91 "V='//cvol//'"'
      if(iappl.eq.4)write(ifhi,'(a,a)')'text 0.15 0.91 "E='//cen//'"'
      write(ifhi,'(a)')       'array 2'

         do j=1,nbin
      dat=datb(j)/nevent/da
      write(ifhi,'(2e13.5)')data(j),dat
         enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto'

      return

           endif

       end

c-----------------------------------------------------------------------
      subroutine xhgcen
c-----------------------------------------------------------------------
c  creates energy spectrum plot for decayed QM-droplet
c  using grand canonical results
c input:
c  xpar1 specifies particle species by paige id, 0 for all
c  xpar2 and xpar3 specify xrange of plot
c  xpar4 specifies line type : dashed (0), dotted (1), full (2) dado (3)
c  xpar5 specifies statistics to be used ,(0) same as iostat
c                                         (1) boltzmann
c output:
c  histo-file
c  newpage, zone and plot commands not included !!!
c-----------------------------------------------------------------------
      include 'epos.inc'
      common/citer/iter,itermx
      parameter (nbin=200)
      real datx(nbin),daty(nbin)
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cbol/rmsbol(mspecs),ptlbol(mspecs),chebol(mspecs),tembol
      character ctem*5,cit*5,cen*6,cvo*6,chem*5

      idpa=nint(xpar1)
      x1=xpar2
      x2=xpar3
      ltyp=nint(xpar4)
      ist=nint(xpar5)
      if(ist.eq.0.and.iostat.eq.1)ist=1

      id=0
      jx=100
      do i=1,nspecs
      if(ispecs(i).eq.idpa)id=i
      enddo

      dx=(x2-x1)/2./jx
      x0=x1+dx

         do j=1,jx
         datx(j)=x0+(j-1)*dx*2.
         daty(j)=0.0

       if(id.eq.0)then

      do 10 i=1,nspecs
      dnde=0.0
        if(datx(j).ge.aspecs(i))then
      x=100.
      if(tem.ne.0.0.and.ist.eq.0)x=(datx(j)-chemgc(i))/tem
      if(tem.ne.0.0.and.ist.eq.1)x=(datx(j)-chebol(i))/tembol
      igsp=gspecs(i)
       if(x.ge.60)goto10
       if(mod(igsp,2).eq.0.and.ist.eq.0)then
      dnde=1./(exp(x)+1.)
       elseif(x.le.1.e-7.and.ist.eq.0)then
      dnde=1.e7
       elseif(ist.eq.0)then
      dnde=1./(exp(x)-1.)
       elseif(ist.eq.1)then
      dnde=exp(-x)
       endif
        endif
      daty(j)=daty(j)+dnde*gspecs(i)*volu/hquer**3/8./pi**3
10    continue

       else

      dnde=0.0
        if(datx(j).ge.aspecs(id))then
      x=100.
      if(tem.ne.0.0.and.ist.eq.0)x=(datx(j)-chemgc(id))/tem
      if(tem.ne.0.0.and.ist.eq.1)x=(datx(j)-chebol(id))/tembol
      igsp=gspecs(id)
       if(x.ge.60)goto11
       if(mod(igsp,2).eq.0.and.ist.eq.0)then
      dnde=1./(exp(x)+1.)
       elseif(x.le.1.e-7.and.ist.eq.0)then
      dnde=1.e7
       elseif(ist.eq.0)then
      dnde=1./(exp(x)-1.)
       elseif(ist.eq.1)then
      dnde=exp(-x)
       endif
        endif
11    daty(j)=dnde*gspecs(id)*volu/hquer**3/8./pi**3

       endif

         enddo

      ctem='     '
      chem='     '
      if(tem.gt.0.)write(ctem,'(f5.3)')tem
      write(cen,'(f6.1)')tecm
      write(cvo,'(f6.1)')volu
      if(id.gt.0)write(chem,'(f5.3)')chemgc(id)
      write(cit,'(i5)')itermx
      write(ifhi,'(a)')       'openhisto'
      if(ltyp.eq.0)then
      write(ifhi,'(a)')       'htyp lda'
      elseif(ltyp.eq.1)then
      write(ifhi,'(a)')       'htyp ldo'
      elseif(ltyp.eq.2)then
      write(ifhi,'(a)')       'htyp lfu'
      elseif(ltyp.eq.3)then
      write(ifhi,'(a)')       'htyp ldd'
      endif
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis E?[n]! (GeV)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dN?[n]!/d^3!p"'
      write(ifhi,'(a,a)')     'text 0.3 0.10 "T='//ctem//'"'
      write(ifhi,'(a,a)')     'text 0.3 0.20 "[m]?[n]!='//chem//'"'
      write(ifhi,'(a,a)')     'text 0.3 0.20 "i?max!='//cit//'"'
      if(iocite.ne.1)then
      write(ifhi,'(a,a)')     'text 0.4 0.91 "V='//cvo//'"'
      write(ifhi,'(a,a)')     'text 0.15 0.91 "E='//cen//'"'
      endif
      write(ifhi,'(a)')       'array 2'

         do j=1,jx
      write(ifhi,'(2e12.4)')datx(j),daty(j)
         enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto'

      return
      end

c-----------------------------------------------------------------------
      subroutine xhgcfl(u,d,s,iii)
c-----------------------------------------------------------------------
c creates unnormalized histogram for net flavor content of grand
c canonically generated sample
c xpar1: specifies width of plot, netflavor centered
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter(nb=200)
      common/cfldat/data(nb),datb(nb),datc(nb),datu(nb)
     *,datd(nb),dats(nb)
      parameter(mxclu=10000)
      integer ku(mxclu),kd(mxclu),ks(mxclu)
      character cfl*3,cen*6,cvol*6
      save ku,kd,ks

      if(iii.eq.0)then

      ku(nrclu)=u
      kd(nrclu)=d
      ks(nrclu)=s

      return

      elseif(iii.lt.0)then

      kwid=nint(xpar1)
      nbin=2*kwid+1
      x1u=keu-kwid
      x2u=keu+kwid
      x1d=ked-kwid
      x2d=ked+kwid
      x1s=kes-kwid
      x2s=kes+kwid
      write(cen,'(f6.1)')tecm
      write(cvol,'(f6.1)')volu

      do i=1,nbin
      data(i)=x1u+(i-1)
      datb(i)=x1d+(i-1)
      datc(i)=x1s+(i-1)
      datu(i)=0.0
      datd(i)=0.0
      dats(i)=0.0
      enddo

      do i=1,nrclu
      nbu=(ku(i)-x1u+1)
      nbd=(kd(i)-x1d+1)
      nbs=(ks(i)-x1s+1)
      if(nbu.le.nbin.and.nbu.ge.1)datu(nbu)=datu(nbu)+1
      if(nbd.le.nbin.and.nbd.ge.1)datd(nbd)=datd(nbd)+1
      if(nbs.le.nbin.and.nbs.ge.1)dats(nbs)=dats(nbs)+1
      enddo

      write(ifhi,'(a)')       'newpage zone 1 3 1'

      write(cfl,'(i3)')keu
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp his'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1u,x2u
      write(ifhi,'(a)')    'text 0 0 "xaxis net u content"'
      write(ifhi,'(a)')    'text 0 0 "yaxis N"'
      write(ifhi,'(a,a)')     'text 0.4 0.91 "V='//cvol//'"'
      write(ifhi,'(a,a)')     'text 0.15 0.91 "E='//cen//'"'
      write(ifhi,'(a,a)')     'text 0.65 0.91 "N?u!='//cfl//'"'
      write(ifhi,'(a)')       'array 2'

         do j=1,nbin
      write(ifhi,'(2e13.5)')data(j),datu(j)
         enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      write(cfl,'(i3)')ked
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp his'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1d,x2d
      write(ifhi,'(a)')    'text 0 0 "xaxis net d content"'
      write(ifhi,'(a)')    'text 0 0 "yaxis N"'
      write(ifhi,'(a,a)')     'text 0.4 0.91 "V='//cvol//'"'
      write(ifhi,'(a,a)')     'text 0.15 0.91 "E='//cen//'"'
      write(ifhi,'(a,a)')     'text 0.65 0.91 "N?d!='//cfl//'"'
      write(ifhi,'(a)')       'array 2'

         do j=1,nbin
      write(ifhi,'(2e13.5)')datb(j),datd(j)
         enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      write(cfl,'(i3)')kes
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp his'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1s,x2s
      write(ifhi,'(a)')    'text 0 0 "xaxis net s content"'
      write(ifhi,'(a)')    'text 0 0 "yaxis N"'
      write(ifhi,'(a,a)')     'text 0.4 0.91 "V='//cvol//'"'
      write(ifhi,'(a,a)')     'text 0.15 0.91 "E='//cen//'"'
      write(ifhi,'(a,a)')     'text 0.65 0.91 "N?s!='//cfl//'"'
      write(ifhi,'(a)')       'array 2'

         do j=1,nbin
      write(ifhi,'(2e13.5)')datc(j),dats(j)
         enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      return

           endif

       end

c-----------------------------------------------------------------------
      subroutine xhgcmt
c-----------------------------------------------------------------------
c creates transverse mass spectrum for QM-droplet decay
c according to grand canonical results
c input:
c  xpar1 specifies particle species by paige id, 0 for all
c  xpar2 and xpar3 specify xrange of plot
c  xpar4 specifies line type : dashed (0), dotted (1), full (2)
c output:
c  histo-file
c  newpage, zone and plot commands not included !!!
c-----------------------------------------------------------------------
      include 'epos.inc'
      common/citer/iter,itermx
      parameter (nbin=200)
      real datx(nbin),daty(nbin)
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      character cen*6,cvo*6,cit*5,ctem*5

      idpa=nint(xpar1)
      x1=xpar2
      x2=xpar3
      ltyp=nint(xpar4)

      id=0
      jx=100
      do i=1,nspecs
      if(ispecs(i).eq.idpa)id=i
      enddo

      dx=(x2-x1)/2./jx
      x0=x1+dx

         do j=1,jx
         datx(j)=x0+(j-1)*dx*2.
         daty(j)=0.0

       if(id.eq.0)then

      do 10 i=1,nspecs
      dndmt=0.0
      if(datx(j).ge.aspecs(i))then
      x=100.
      xx=100.
      if(tem.ne.0.)x=datx(j)/tem
      if(tem.ne.0.)xx=chemgc(i)/tem
      if(abs(xx).le.60)dndmt=gspecs(i)*volu/hquer**3*exp(xx)*datx(j)
     */4./pi**3*hgcbk1(x)
      endif
      daty(j)=daty(j)+dndmt
10    continue

       else

      dndmt=0.0
      if(datx(j).ge.aspecs(id))then
      x=100.
      xx=100.
      if(tem.ne.0.)x=datx(j)/tem
      if(tem.ne.0.)xx=chemgc(id)/tem
      if(abs(xx).le.60)dndmt=gspecs(id)*volu/hquer**3*exp(xx)*datx(j)
     */4./pi**3*hgcbk1(x)
      endif
      daty(j)=dndmt

       endif

         enddo

      write(cit,'(i5)')itermx
      write(cen,'(f6.1)')tecm
      write(cvo,'(f6.1)')volu
      write(ctem,'(f5.3)')tem
      write(ifhi,'(a)')       'openhisto'
      if(ltyp.eq.0)then
      write(ifhi,'(a)')       'htyp lda'
      elseif(ltyp.eq.1)then
      write(ifhi,'(a)')       'htyp ldo'
      elseif(ltyp.eq.2)then
      write(ifhi,'(a)')       'htyp lfu'
      endif
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis m?t! (GeV)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dN?[n]!/d^2!m?t! "'
      write(ifhi,'(a,a)')     'text 0.3 0.10 "T='//ctem//'"'
      write(ifhi,'(a,a)')     'text 0.3 0.20 "i?max!='//cit//'"'
      write(ifhi,'(a,a)')     'text 0.4 0.91 "V='//cvo//'"'
      write(ifhi,'(a,a)')     'text 0.15 0.91 "E='//cen//'"'
      write(ifhi,'(a)')       'array 2'

         do j=1,jx
      write(ifhi,'(2e12.4)')datx(j),daty(j)
         enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto'

      return
      end

c-----------------------------------------------------------------------
      subroutine xhgcmu
c-----------------------------------------------------------------------
c creates multiplicity plot for decayed QM-droplet
c according to grand canonical results
c input:
c  xpar1 specifies species by paige id, 0 for total multiplicity
c  xpar2 specifies xrange to be set automatically (0) or by hand (1)
c  xpar3 and xpar4 xrange if xpar2 ne 0
c  xpar5 xrange = average+-sigma*xpar5
c  xpar6 specifies line type : dashed (0), dotted (1), full (2)
c  xpar7 specifies statistics : same as iostat (0)
c                               boltzmann (1)
c output:
c  histo-file
c  newpage, zone and plot commands not included !!!
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (nbin=200)
      real datx(nbin),daty(nbin)
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cbol/rmsbol(mspecs),ptlbol(mspecs),chebol(mspecs),tembol
      common/cgctot/rmstot,ptltot
      character cyield*8,cen*6,cvo*6,cinco*1


      idpa=nint(xpar1)
      ixra=nint(xpar2)
      iwid=nint(xpar5)
      ltyp=nint(xpar6)
      ist=nint(xpar7)
      if(ist.eq.0.and.iostat.eq.1)ist=1


      pn=0.0
      id=0
      jx=100
      ymin=1./nevent/10.
      if(nevent.le.10)ymin=ymin/10.
      do i=1,nspecs
      if(ispecs(i).eq.idpa)id=i
      enddo

       if(ixra.eq.1)then
      x1=anint(xpar3)
      x2=anint(xpar4)
       else
      if(id.eq.0)then
      x1=anint(ptltot-iwid*rmstot)
      x2=anint(ptltot+iwid*rmstot)
      else
      x1=anint(ptlngc(id)-iwid*rmsngc(id))
      x2=anint(ptlngc(id)+iwid*rmsngc(id))
      endif
      x2=max(x2,3.0)
       endif

      x1=max(x1,0.0)
      dx=(x2-x1)/2./jx
      x0=x1+dx

      do j=1,jx
      datx(j)=x0+(j-1)*dx*2.
      if(id.eq.0)then

c     total multiplicity
c     ------------------
      x=100.
      if(rmstot.ge.1.e-10)x=(datx(j)-ptltot)**2/rmstot**2/2.

       if(x.ge.60)then
      pn=0.0
       else
      pn=exp(-x)/rmstot/sqrt(2.*pi)
       endif

      daty(j)=pn

         else

c     one species (specified by id)
c     ------------------------------
      x=100.
      if(rmsngc(id).ge.1.e-10.and.ist.eq.0)
     *x=(datx(j)-ptlngc(id))**2/rmsngc(id)**2/2.
      if(rmsbol(id).ge.1.e-10.and.ist.eq.1)
     *x=(datx(j)-ptlbol(id))**2/rmsbol(id)**2/2.

       if(x.ge.60)then
      pn=0.0
       else
      pn=0.0
      if(ist.eq.0)pn=exp(-x)/rmsngc(id)/sqrt(2*pi)
      if(ist.eq.1)pn=exp(-x)/rmsbol(id)/sqrt(2*pi)
       endif

      daty(j)=pn

         endif
         enddo

      if(id.eq.0)then
      write(cyield,'(f8.3)')ptltot
      else
      write(cyield,'(f8.3)')ptlngc(id)
      endif
      write(cinco,'(i1)')ioinco
      write(cen,'(f6.1)')tecm
      write(cvo,'(f6.1)')volu
      write(ifhi,'(a)')       'openhisto'
      if(ltyp.eq.0)then
      write(ifhi,'(a)')       'htyp lda'
      elseif(ltyp.eq.1)then
      write(ifhi,'(a)')       'htyp ldo'
      elseif(ltyp.eq.2)then
      write(ifhi,'(a)')       'htyp lfu'
      elseif(ltyp.eq.3)then
      write(ifhi,'(a)')       'htyp ldd'
      endif
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a,e11.3,a)')'yrange',ymin,'  auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis N?[n]!"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(N?[n]!)"'
      write(ifhi,'(a,a)')'text 0.3 0.10 "" "L#N?[n]!"G#='//cyield//'""'
      write(ifhi,'(a,a)')     'text 0.3 0.2 "conf?in!='//cinco//'"'
      if(iocite.ne.1)then
      write(ifhi,'(a,a)')     'text 0.4 0.91 "V='//cvo//'"'
      write(ifhi,'(a,a)')     'text 0.15 0.91 "E='//cen//'"'
      endif
      write(ifhi,'(a)')       'array 2'

         do j=1,jx
      write(ifhi,'(2e12.4)')datx(j),daty(j)
         enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto'


      return
      end


c-----------------------------------------------------------------------
      subroutine xhgcmx
c-----------------------------------------------------------------------
c creates multiplicity plot for decayed QM-droplet
c according to grand canonical results POISSON DISTRIB.!!!!
c input:
c  xpar1 specifies species by paige id, 0 for total multiplicity
c  xpar2 specifies xrange to be set automatically (0) or by hand (1)
c  xpar3 and xpar4 xrange if xpar2 ne 0
c  xpar5 xrange = average+-sigma*xpar5
c  xpar6 specifies line type : dashed (0), dotted (1), full (2) dado (3)
c  xpar7 specifies statistics : same as iostat (0)
c                               boltzmann (1)
c output:
c  histo-file
c  newpage, zone and plot commands not included !!!
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (nbin=200)
      real datx(nbin),daty(nbin)
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cbol/rmsbol(mspecs),ptlbol(mspecs),chebol(mspecs),tembol
      common/cgctot/rmstot,ptltot
      character cyield*8,cen*6,cvo*6,cinco*1


      idpa=nint(xpar1)
      ixra=nint(xpar2)
      iwid=nint(xpar5)
      ltyp=nint(xpar6)
      ist=nint(xpar7)
      if(ist.eq.0.and.iostat.eq.1)ist=1
      pn=0.


      id=0
      ymin=1./nevent/10.
      if(nevent.le.10)ymin=ymin/10.
      do i=1,nspecs
      if(ispecs(i).eq.idpa)id=i
      enddo

       if(ixra.eq.1)then
      n1=nint(xpar3)
      n2=nint(xpar4)
       else
      if(id.eq.0)then
      n1=nint(ptltot-iwid*rmstot)
      n2=nint(ptltot+iwid*rmstot)
      else
      n1=nint(ptlngc(id)-iwid*rmsngc(id))
      n2=nint(ptlngc(id)+iwid*rmsngc(id))
      endif
      n2=max(n2,3)
       endif

      n1=max(n1,0)
      jx=n2+1

      do j=1,jx
      datx(j)=j-1
      jf=1
      if(j.gt.1)then
      do i=1,j-1
      jf=jf*i
      enddo
      endif
      if(id.eq.0)then

c     total multiplicity
c     ------------------

      daty(j)=1./jf*ptltot**(j-1)*exp(-ptltot)

         else

c     one species (specified by id)
c     ------------------------------

      if(ist.eq.0)pn=1./jf*ptlngc(id)**(j-1)*exp(-ptlngc(id))
      if(ist.eq.1)pn=1./jf*ptlbol(id)**(j-1)*exp(-ptlbol(id))

      daty(j)=pn

         endif
         enddo

      if(id.eq.0)then
      write(cyield,'(f8.3)')ptltot
      else
      write(cyield,'(f8.3)')ptlngc(id)
      endif
      write(cinco,'(i1)')ioinco
      write(cen,'(f6.1)')tecm
      write(cvo,'(f6.1)')volu
      write(ifhi,'(a)')       'openhisto'
      if(ltyp.eq.0)then
      write(ifhi,'(a)')       'htyp lda'
      elseif(ltyp.eq.1)then
      write(ifhi,'(a)')       'htyp ldo'
      elseif(ltyp.eq.2)then
      write(ifhi,'(a)')       'htyp lfu'
      elseif(ltyp.eq.3)then
      write(ifhi,'(a)')       'htyp ldd'
      endif
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2i3)')'xrange',n1,n2
      write(ifhi,'(a,e11.3,a)')'yrange',ymin,'  auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis N?[n]!"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(N?[n]!)"'
      write(ifhi,'(a,a)')'text 0.3 0.10 "" "L#N?[n]!"G#='//cyield//'""'
      write(ifhi,'(a,a)')     'text 0.3 0.2 "conf?in!='//cinco//'"'
      if(iocite.ne.1)then
      write(ifhi,'(a,a)')     'text 0.4 0.91 "V='//cvo//'"'
      write(ifhi,'(a,a)')     'text 0.15 0.91 "E='//cen//'"'
      endif
      write(ifhi,'(a)')       'array 2'

         do j=1,jx
      write(ifhi,'(2e12.4)')datx(j),daty(j)
         enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto'


      return
      end

c-----------------------------------------------------------------------
      subroutine xhgcpt
c-----------------------------------------------------------------------
c creates transverse momentum spectrum for decayed QM-droplet
c according to grand canonical results
c input:
c  xpar1 specifies particle species by paige id, 0 for all
c  xpar2 rapidity window
c  xpar3 and xpar4 specify xrange of plot
c  xpar5 specifies line type : dashed (0), dotted (1), full (2)
c output:
c  histo-file
c  newpage, zone and plot commands not included !!!
c-----------------------------------------------------------------------
      include 'epos.inc'
      common/citer/iter,itermx
      parameter (nbin=200)
      real datx(nbin),daty(nbin)
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      character crap*5,cen*6,cvo*6,cit*5

      idpa=nint(xpar1)
      y=xpar2
      x1=xpar3
      x2=xpar4
      ltyp=xpar5

      write(crap,'(f5.1)')y
      id=0
      jx=100
      do i=1,nspecs
      if(ispecs(i).eq.idpa)id=i
      enddo

      dx=(x2-x1)/2./jx
      x0=x1+dx

         do j=1,jx
         datx(j)=x0+(j-1)*dx*2.
         daty(j)=0.0

       if(id.eq.0)then

      do 10 i=1,nspecs
      x=100.
      if(tem.ne.0.)
     *x=(sqrt(aspecs(i)**2+datx(j)**2)*cosh(y)-chemgc(i))/tem
       if(x.ge.60)then
      dndpt=0.0
       else
      dndpt=exp(-x)
       endif
      dndpt=dndpt*gspecs(i)*volu/hquer**3*cosh(y)
     **sqrt(aspecs(i)**2+datx(j)**2)/8./pi**3
      daty(j)=daty(j)+dndpt
10    continue

       else

      x=100.
      if(tem.ne.0.)
     *x=(sqrt(aspecs(id)**2+datx(j)**2)*cosh(y)-chemgc(id))/tem
       if(x.ge.60)then
      dndpt=0.0
       else
      dndpt=exp(-x)
       endif
      dndpt=dndpt*gspecs(id)*volu/hquer**3*cosh(y)
     **sqrt(aspecs(id)**2+datx(j)**2)/8./pi**3
      daty(j)=dndpt

       endif

         enddo

      write(cit,'(i5)')itermx
      write(cen,'(f6.1)')tecm
      write(cvo,'(f6.1)')volu
      write(ifhi,'(a)')       'openhisto'
      if(ltyp.eq.0)then
      write(ifhi,'(a)')       'htyp lda'
      elseif(ltyp.eq.1)then
      write(ifhi,'(a)')       'htyp ldo'
      elseif(ltyp.eq.2)then
      write(ifhi,'(a)')       'htyp lfu'
      endif
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis p?t! (GeV/c)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dN?[n]!/dyd^2!p?t!"'
      write(ifhi,'(a)')    'text 0.10 0.10 "y = '//crap//'"'
      write(ifhi,'(a)')    'text 0.10 0.30 "i?max! = '//cit//'"'
      write(ifhi,'(a,a)')     'text 0.4 0.91 "V='//cvo//'"'
      write(ifhi,'(a,a)')     'text 0.15 0.91 "E='//cen//'"'
      write(ifhi,'(a)')       'array 2'

         do j=1,jx
      write(ifhi,'(2e12.4)')datx(j),daty(j)
         enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto'

      return
      end

c-----------------------------------------------------------------------
      subroutine xhgcra
c-----------------------------------------------------------------------
c creates rapidity distribution for decayed QM-droplet
c according to grand canonical results
c input:
c  xpar1 specifies particle species by paige id, 0 for all
c  xpar2 and xpar3 specify xrange of plot
c  xpar4 specifies line type : dashed (0), dotted (1), full (2)
c output:
c  histo-file
c  newpage, zone and plot commands not included !!!
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (nbin=200)
      real datx(nbin),daty(nbin)
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cgctot/rmstot,ptltot
      character cen*6,cvo*6,cng*8

      idpa=nint(xpar1)
      x1=nint(xpar2)
      x2=nint(xpar3)
      ltyp=nint(xpar4)

      id=0
      jx=100
      ymin=1./nevent/10.
      if(nevent.le.10)ymin=ymin/10.
      do i=1,nspecs
      if(ispecs(i).eq.idpa)id=i
      enddo

      dx=(x2-x1)/2./jx
      x0=x1+dx

         do j=1,jx

         datx(j)=x0+(j-1)*dx*2.
         daty(j)=0.0
         y=datx(j)
         if(ish.ge.9)write(ifch,*)'cosh y:',cosh(y)

       if(id.eq.0)then

      do 10 i=1,nspecs
      dndy=0.0
      sum=aspecs(i)**2*tem+2.*aspecs(i)*tem**2/cosh(y)
     *+2.*tem**3/cosh(y)**2
      x=100.
      if(tem.ne.0.0)
     *x=(aspecs(i)*cosh(y)-chemgc(i))/tem

       if(x.ge.60.)then
      pro=0.0
       else
      pro=exp(-x)
      endif

      pro=pro*gspecs(i)*volu/hquer**3/4./pi**2

      if(pro.ge.(1.e-30).and.sum.ge.(1.e-30))then
      che=alog(pro)+alog(sum)
      else
      che=-61.0
      endif
      if(che.le.60.0.and.che.ge.(-60.0))dndy=pro*sum
c     if(che.le.60.0.and.che.ge.(-60.0))dndy=exp(che)

      daty(j)=daty(j)+dndy

10    continue

       else

      dndy=0.0
      sum=aspecs(id)**2*tem+2.*aspecs(id)*tem**2/cosh(y)
     *+2.*tem**3/cosh(y)**2
      x=100.
      if(tem.ne.0.0)
     *x=(aspecs(id)*cosh(y)-chemgc(id))/tem

       if(x.ge.60.)then
      pro=0.0
       else
      pro=exp(-x)
      endif

      pro=pro*gspecs(id)*volu/hquer**3/4./pi**2

      if(pro.ge.(1.e-30).and.sum.ge.(1.e-30))then
      che=alog(pro)+alog(sum)
      else
      che=-61.0
      endif
      if(che.le.60..and.che.ge.-60.)dndy=pro*sum

      daty(j)=dndy

       endif

         enddo

      write(cen,'(f6.1)')tecm
      write(cvo,'(f6.1)')volu
      if(id.eq.0)then
      write(cng,'(f8.3)')ptltot
      else
      write(cng,'(f8.3)')ptlngc(id)
      endif
      write(ifhi,'(a)')       'openhisto'
      if(ltyp.eq.0)then
      write(ifhi,'(a)')       'htyp lda'
      elseif(ltyp.eq.1)then
      write(ifhi,'(a)')       'htyp ldo'
      elseif(ltyp.eq.2)then
      write(ifhi,'(a)')       'htyp lfu'
      endif

      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a,e11.3,a)')'yrange',ymin,'  auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis y"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dN?[n]!/dy"'
      write(ifhi,'(a,a)')     'text 0.4 0.91 "V='//cvo//'"'
      write(ifhi,'(a,a)')     'text 0.15 0.91 "E='//cen//'"'
      write(ifhi,'(a,a)')     'text 0.3 0.10 "N?[n]!='//cng//'"'
      write(ifhi,'(a)')       'array 2'

         do j=1,jx
      write(ifhi,'(2e12.4)')datx(j),daty(j)
         enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto'

      return
      end

















