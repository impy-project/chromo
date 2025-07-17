C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c----------------------------------------------------------------------
      subroutine hnbaaa(jcasi,ip,iret)  !former hnbaaa156 from epos-yyy
c----------------------------------------------------------------------
c  Decay of object ip, where the type of the object is:
c        
c  jcase = 3 ... cluster from parametrized flow approach (iorsdf=3)  
c          6 ... cluster from rope approach for initial conditions 
c          7 ... volume elements at FO hypersurface in local rest frame
c          8 ... cluster decay without mass rescaling for testing 
c          9 ... remnant droplets
c----------------------------------------------------------------------
c  two methods :
c     NRPS via LIPS (if yiespecs <= iocova)             
c     NRPS via Hagedorn method (if yiespecs > iocova) 
c----------------------------------------------------------------------
c  Structure: 
c   1) Rescaling mass depending on jcase
c   2) N-body decay (loop over hnbmet) 
c         & flow given to particles
c----------------------------------------------------------------------
c  Mass is reduced to allow longitudinal (to smear the matrix effect)
c  and radial flow but in case of low mass clusters, flow is forced to 
c  get the energy given by fradflo*mass (not enough particles)
c  particles are still there after the boost. So we should not try to
c  restore the original mass with radial flow, but define radial flow 
c  first and then complete the energy by rescaling pz.
c  Multiplicity of droplet is very high with narrow distribution around
c  eta=0 if nothing done. Long flow reduce multiplicity and elongate 
c  distribution in eta.
c----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incico'
      include 'epos.incho'
      real bglaub
      common/cbglaub/bglaub
      common/cxyzt/xptl(mxptl),yptl(mxptl),zptl(mxptl)
     * ,tptl(mxptl),optl(mxptl),uptl(mxptl),sptl(mxptl)
     *,rptl(mxptl,3)
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      common/citer/iter,itermx
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common /cttaus/  tpro,zpro,ttar,ztar,ttaus,detap,detat
      integer jc(nflav,2)
      parameter (mspecs=400)
      parameter (mxpair=mspecs**2*4)
      double precision p(5),c(5),pa(5),yrmax
      real piii(5)
      parameter(maxit=500000)
      common/count/nacc,nrej,naccit(maxit),nptot,npit(maxit)
      dimension be(4),pe(5)
      double precision pall(4),am,ptm
      common/cspecs2/pspecs(0:mspecs+1),yie2volu,yiespecs
      common/yradx/yrad(maxp),phirad(maxp)
      common/cylng/ylng(maxp)
c      real Y3(maxp)
      common/xxxspecsy/ndrop(-4:4,-4:4,-4:4)
      common/cdelzet/delzet,delsce
      common/cvocell/vocell,xlongcell
      real xcell,scell,zcell,delxce
      integer jjj,m1cell,m3cell,nptlb,nptla
      common/jintpoc2/xcell,scell,zcell,delxce
     &,jjj,m1cell,m3cell,nptlb,nptla
      real aminclust
      common/jintpoc4/aminclust
      
      call utpri('hnbaaa',ish,ishini,4)
      jcase=jcasi
      ntry=0
 10   continue
      ntry=ntry+1
      iret=0
      if(ish.ge.4)then
      write(ifch,140)sngl(ttaus)
  140 format(/' ----------------------------------'/
     *        '      NB decay  at tau =',f6.2/
     *        ' ----------------------------------')
      write(ifch,*)'droplet:'
      call clist('&',ip,ip,60,60)
      endif

    !~~~~~~~~~~~~~~ basic settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      ianisotropy = 0
      if(jcase.eq.3)ianisotropy = 1

      epsi = efrout !energy density used to determine volu

    !~~~~~~~~~~~~~~ useful definitions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do j=1,5
      c(j)=pptl(j,ip)
      enddo

      call idquac(ip,nqi,nsi,nai,jc)
      keu=jc(1,1)-jc(1,2)
      ked=jc(2,1)-jc(2,2)
      kes=jc(3,1)-jc(3,2)
      kec=jc(4,1)-jc(4,2)
      keb=jc(5,1)-jc(5,2)
      ket=jc(6,1)-jc(6,2)
c      write(ifch,*)'hnb :  uds=',keu,ked,kes,'   E=',pptl(5,ip)

      amin0=utamnu(keu,ked,kes,kec,keb,ket,4)   !utamnu(...,4) and not utamnu(...,5)
                                               !could be too light after flow
c      print *,'hnb :  uds=',keu,ked,kes,' E=',pptl(5,ip), ' amin=',amin
      tecmor=pptl(5,ip)
      tecm=pptl(5,ip)
      dens=radptl(ip)
      ityip=ityptl(ip)
      if(tecm.le.amin0)then
        iret=1
        if(ish.ge.4)write(ifch,*)'Decay skipped (M too low) !'
     &                           ,tecm,amin0
        goto 1000
      endif
      tecmx=tecm
      tecmxx=tecm

      yco=0.
      corrco=1.
      ncx=1
      yrmax=0d0
      fradflo=1.
      fecc=0.
      taufo=abs(tauzer2)
      Z=0.
      aa=1.
      bb=0.
      cc=0.
      dd=1.

    !###################################################################
    !###################################################################
    !                   Rescaling for different cases
    !###################################################################
    !###################################################################        

        ! redefine energy to account for collective flow
        ! \int_0^yrmax f(y) d2y = E_new (effective mass)
        ! \int_0^yrmax cosh(y) f(y) d2y = E_old

    !~~~~~~~~~~~~~~~~~~~~~~~~~~ 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if(jcase.eq.9)then !for droplets from remnants

        !split droplet in ncx subdroplets 
c       amin=amin0!+max(0.,aminclust)
        amin=max(1.2*amin0,amdrmin)
        xmxms=max(amin,100.)
        if(tecm.gt.2.*xmxms)then
          ncx=nint(tecm/xmxms)
          tecm=tecmx/float(ncx)
          do while(tecm.le.amin)
            ncx=ncx-1
            if(ncx.eq.0)call utstop("Problem in hnbaa with min mass!&")
            tecm=tecmx/float(ncx)
          enddo
        else
          ncx=1
        endif

        yco=10.*abs(ylongmx)  !10 is to have the same minimum yco than for core
        tecmx=tecm

        if(yco.gt.0.)then
          corr=yco/sinh(yco)
          if(tecm*corr.le.amin)then
            do while(tecm*corr.lt.amin.and.yco.gt.1e-2)
              yco=yco*0.9
              corr=1.
              if(yco.gt.1e-2)corr=yco/sinh(yco)
c             print *,'corr',yco,corr,tecm*corr,amin
            enddo
          endif
          tecm=tecm*corr
        endif
        
    !~~~~~~~~~~~~~~~~~~~~~~~~~~ 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      elseif(jcase.eq.3)then ! Effective flow option

c        Z=max(0.,ng1evt-2.)/400. !ng1evt=Npart  
c        taufo=taufo*(1+Z*3.)
        taufo=taufo*ratiomass**0.66 !0.33    !ratiodens**0.33 vary from ~0.75 to  ~5
c        taufo=taufo*(1.+ratiodens/50.)    !ratiodens**0.33 vary from ~0.75 to  ~5
c        print *,taufo,Z,ratiomass**0.33,ratiomass,ratiomass*ratiodens
        amin=max(1.2*amin0,max(amdrmin,abs(aminclust)))

c        yrmax=0d0
c        yrmax=dble(yrmaxi*dens/avgdens)
c       yrmax=dble(yrmaxi*min(dens/avgdens

        yrmax=dble(yrmaxi)
        if(yradpx.gt.0.)then
          yrmax=yrmax*dble(min(max(dens/avgdens,0.)
     .                    ,max(1.,yradpx/(ng1evt*yrmax))))
c        .                       ,1.+yradpx*exp(avgdens)))
        else
          if(rangen().lt.(dens/avgdens)**max(-5.,yradpx
     .                                    *max(1,ng1evt-1)))
     .        yrmax=yrmax*dble(dens)/dble(avgdens)
        endif

c        print *,ng1evt,'yrmax',yrmaxi,yrmax,dens/avgdens,yradpx
c     .,(dens/avgdens)**(max(-5.,yradpx*max(1,ng1evt-1))),yfgeo
c        fradflo=fradflii
ccc        dens=min(tecmor,yrmaxi)/amuseg
c        density=tecmor
c        volr=amuseg!*float(max(1,koievt))**yradpx!*max(bglaub**3,0.)!vocell*ydsrd*dezptl(ipo)**yradpx)
cc        voll=max(0.,ydslg*dezptl(ipo)/m1cell**2*m3cell)**yradpi
c        voll=amuseg !max(1.,vocell*ydslg*dezptl(ipo)**yradpi)!(0.5*float(npjevt+ntgevt))**yradpi
cc        voll=float(max(1,koievt))**yradmi
ccc        voll=1.+yradpi*log(max(1.,vocell*ydslg*dezptl(ipo)))
ccc        if(voll.gt.bglaub**3)print *,'----->',tecmor
c        dens=density/volr
c        if(dens.gt.1.)then
c          fgeo=1.+ydsrd*float(abs(npjevt-ntgevt))/float(npjevt+ntgevt)
c          yrmax=dble(yradpp/max(bglaub,1.))
c          yrmax=max(yrmax,dble(fgeo*yradmx*log(dens)))
cc     .                        ,dble(yrmaxi))
c          yrmax=yrmax/dble(0.5*float(npjevt+ntgevt))**yradpx
cc        yrmax=yrmax*(2.5d0/dble(xmaxico))**yradpx
c        endif
        if(yrmax.gt.1d-5)then
          fradflo=sngl(1d0/
     &         ((sinh(yrmax)*yrmax-cosh(yrmax)+1d0)/(yrmax**2/2d0)))
        else
          fradflo=1.
          yrmax=0d0
        endif
c
cc local definition of yco 
cc using local definition increase multiplicity at large eta : less good for LHCb mult distributions but for pPb it is necessary because otherwise the pseudorapidity shape is deformed (small and large cluster have the same reduction of mass)
c        if(ylongmx.lt.0.)then
cc          yco=yfgeo            !* 1.75
ccc          dens=min(tecmor,yfgeo)
cc use a mix of local and global definition : if ydslg=0, use tecmor, if ydslg=1 us tecmor but with reduction due to participant number. In between a mixed is used where the pseudorapidity dep of the global mass is flattened
c          density=max(0.,density/voll)
cc          dens=max(1.,tecmor)
cc          yco=max(0.,abs(ylongmx)+yradmi*log(dens))
ccc          yco=max(0.,abs(ylongmx)+yradmi*log(min(tecmor,yfgeo)))
cc minimum should increase to reproduce pPb
c          ycoi=max(abs(ylongmx)*nint(dezptl(ip))/10000
cc     .       *(1.+0.15*log(float(max(1,ikoevt))/float(max(1,koievt))))
c     .       ,yradmi*log(density))*xlongcell
cc          print *,'2 ->',bimevt,koievt,ikoevt,yfgeo,density,ycoi
cc     .   ,yradmi*log(density),abs(ylongmx)*nint(dezptl(ip))/10000
cc     .       *(1.+0.15*log(float(max(1,ikoevt))/float(max(1,koievt))))
cc          yco=abs(ylongmx)*xlongcell*yfgeo*log(dens)
cc          yco=yco/(0.5*float(npjevt+ntgevt))**yradpi
c        else
c          ycoi=ylongmx*xlongcell*nint(dezptl(ip))/10000
c        endif
        yco=sign(ycori,dezptl(ip))*nint(dezptl(ip))/10000 !*tecmor/avgmass    !if it is too large it kills the tail of the multiplicity distribution in pp
        if(yco.lt.0.)yco=5.*abs(ylongmx)  !5 is to have the same minimum yco than for core but here is a droplet from remnant, do not reduce mass too much

c
        if(ish.ge.3)write(ifch,*)'initial yrmax,yco,vol='
cc      print *,'yrmax,yco,vol=',dezptl(ipo),yrmaxi,abs(ylongmx)*xlongcell
     .       ,yrmax,yco,dens,yfgeo,ycori,nint(dezptl(ip))/10000
c
c        aumin=amin*(1d0+yrmax)  !for rad and long flow
c        if(tecm.lt.aumin)then
c         print *,'ici',aumin,tecm,amin,yrmax,0.5d0*dble(tecm/amin-1d0)
c         yrmax=0.5d0*dble(tecm/amin-1.) !to have always some flow
c         aumin=amin*(1.+sngl(yrmax)) !for rad and long flow  
c         if(yrmax.gt.1d-5)fradflo=sngl(1d0/
c     &         ((sinh(yrmax)*yrmax-cosh(yrmax)+1d0)/(yrmax**2/2d0)))
c        endif
c        fradflo=min(1.,max(fradflo,1.01*amin/tecm)) !if flow too large, do something anyway (saturation of flow)
c
c        if(yrmax.gt.1d-5)then
c          fradflo=sngl(1d0/
c     &    ((sinh(yrmax)*yrmax-cosh(yrmax)+1d0)/(yrmax**2/2d0)))
c        else
c          fradflo=1.
c          yrmax=0d0
c        endif
c        if(tecm*fradflo.lt.amin)then
c          fradflo=1.
c          yrmax=0d0
c        endif
        corrco=1.
        if(yco.gt.0.)corrco=sinh(yco)/yco
c        corrco=corrcoi
        
        tecmxx=tecm
        corr=fradflo/corrco
        if(tecm*corr.le.amin.and.tecm.gt.0.)then
          do while(tecm*corr.lt.amin.and.yco.gt.1e-2)
            yco=yco*0.9
            corrco=1.
            if(yco.gt.0.)corrco=sinh(yco)/yco
            corr=fradflo/corrco
c            write(ifch,*)'corco',yco,corrco,corr,tecm*corr,amin
          enddo
        endif
c
        if(tecm*corr.ge.amin) then
          tecm=tecm*fradflo
c          if(tecm.lt.amin)stop'aaahnb: small mass. should not happen.'
        else
          do while(tecm*corr.lt.amin.and.yrmax.gt.1d-3)
            yrmax=yrmax*0.9d0
            fradflo=sngl(1d0/
     &    ((sinh(yrmax)*yrmax-cosh(yrmax)+1d0)/(yrmax**2/2d0)))
            corr=fradflo/corrco
c            write(ifch,*)'yrmax',yrmax,corr,tecm*corr,amin
          enddo
          if(yrmax.le.1d-3)then
            if(ish.ge.4)
     &      write(ifch,*)'no more flow',yrmax,corr,tecm*corr,amin,tecmxx
            if(ityip.eq.60)then
              iret=1
              goto 1000
            else        ! if remnant droplet, redo without flow
              jcase=9
              goto 10
            endif
          else
            tecm=tecm*fradflo
          endif
        endif
        
        tecmx=tecm
c        if(yco.gt.0.) then
c          if(tecm.ge.amin)then
            tecm=tecm/corrco
c            do while(tecm.lt.amin)
c              yco=yco*0.5
c              corrco=1.
c              if(yco.gt.0.)corrco=sinh(yco)/yco
c              tecm=tecmx/corrco
c            enddo
c          else
c            tecm=tecmx
c            yco=0.
c          endif
c        else
c          yco=0.
c        endif

        if(ish.ge.2)
     &  write(ifch,'(a,5f7.2,i7,i12,i7)')  '==> cluster energy: '
     &      ,pptl(5,ip),tecmx,tecm,yrmax,yco,ip,idptl(ip),ityptl(ip)

    !###################################################################
    !###################################################################
      else
        call utstop('jcase&')
      endif                     !end jcase 
    !###################################################################
    !###################################################################

      nptlc=nptl

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do nc=1,ncx  ! ncx iteration in case of jcase=9
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if(nc.gt.1)then
        ! emit first droplet with flavors of original droplet 
        ! and then complete with droplets with null net quark content
        keu=0
        ked=0
        kes=0
        kec=0
        keb=0
        ket=0
      endif

      if(ish.ge.4)write(ifch,*)'===== nc / ncx: ',nc,'/',ncx

    !###################################################################
    !###################################################################
    !###################################################################
    !###################################################################
    !           N-body decay & flow to particles
    !-------------------------------------------------------------------
    !  The following code uses
    !  tecm .......... rescaled CMS energy of the decaying object
    !  tecmx ......... tecm before rescaling (long)
    !  tecmxx ........ tecm before rescaling 
    !  epsi .......... energy density used to define volume volu 
    !  yco ........... collective long boost                        
    !  yrmax ......... radial boost
    !  ianisotropy ... anisotropic boost (if 1)
    !###################################################################
    !###################################################################
    !###################################################################
    !###################################################################

    !~~~~~~~~~define volume~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      volu=tecm/epsi
       !print*,tecm,volu
      yiespecs=yie2volu*volu

    !~~~~~~~~~decay~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      call hnbinit(iretx)
    
      if(iretx.ne.0)then
        !write(ifch,*)'***** unsucessfull hnbinit *****'
        iret=2
        goto 1000
      endif
      if(ioinct.ge.1)goto 1

      do iter=1,nint( itermx * fitermet )
        if(iter.le.maxit)naccit(iter)=0
        call hnbmet
      enddo
c check
      eee=0
      do i=1,np
      eee=eee+pcm(4,i)
      enddo
      if(abs(eee-tecm)/tecm.gt.1e-4)write(ifmt,'(a,2e10.2,i3,3x,3i4)')
     .'hnbaaa (1) energy pb',eee,tecm,jcase,keu,ked,kes

1     continue

      if(ioceau.eq.1.and.iappl.eq.1)call xhnbte(ip)
      scal=0.

    !~~~~~~~~~~long coll flow -> particles~~~~~~~~~~~~~~~~

      tecmxxx=tecm
      do i=1,np
        ylng(i)=0
      enddo
c      if(yco.gt.0.) then
c        errlim=0.0001
c        tecm=tecmx
c        niter=0
c 611    energ=0.
c        niter=niter+1
c        do i=1,np
c          pt2=0.!pcm(1,i)**2+pcm(2,i)**2+amass(i)**2
c          ylng(i)=(2*rangen()-1)*yco
c          be(3)=sinh(ylng(i))
c          be(4)=cosh(ylng(i))
c          energ=energ+be(4)*pcm(4,i)-be(3)*pcm(3,i)
c        enddo
c        if(abs(energ-tecm).gt.0.1.and.niter.lt.1000)then
c          goto 611
c        elseif(niter.ge.1000)then
c          if(ish.ge.1)write(ifch,*)'Long Flow failed:'
c     &                             ,energ,tecm
c          yco=0
c          tecm=tecmxxx
c          goto 400
c        endif
c                !print*,'===== energy after flow boosts',energ,'   soll: ',tecm
c        do j=1,4
c          pe(j)=0.
c        enddo
c        do i=1,np
c          be(1)= 0
c          be(2)= 0
c          be(3)= sinh(ylng(i))
c          be(4)= cosh(ylng(i))
c          call utlob3(-1,be(1),be(2),be(3),be(4),1e0
c     *         , pcm(1,i), pcm(2,i), pcm(3,i), pcm(4,i))  !KW2012 first argu -1 !!!
c          do j=1,4
c          pe(j)=pe(j)+pcm(j,i)
c          enddo
c        enddo
c        pe(5)=sqrt(pe(4)**2-pe(3)**2-pe(2)**2-pe(1)**2)
c       !write(6,'(a,5e11.3)')'flow boosts',pe
c        do j=1,4
c          pa(j)=0.d0
c        enddo
c        do i=1,np
c          v=pcm(3,i)/pcm(4,i) 
c          Y2=0.5 * log ((1+v)/(1-v))
c          call utlob3(1,pe(1),pe(2),pe(3),pe(4),pe(5)
c     *         , pcm(1,i), pcm(2,i), pcm(3,i), pcm(4,i))
c          do j=1,4
c            pa(j)=pa(j)+pcm(j,i)
c          enddo
c          v=pcm(3,i)/pcm(4,i) 
c          Y3(i)=0.5 * log ((1+v)/(1-v))
c          ylng(i)=ylng(i)+Y3(i)-Y2  !should be equal Y3-Y1
c        enddo
c        pa(5)=sqrt(pa(4)**2-pa(3)**2-pa(2)**2-pa(1)**2)
c                !write(6,'(a,5e11.3)')' cms boost ',pa
cc rescale all momenta to conserve initial energy
c        esoll=tecm
c        scal=1.
c        do ipass=1,200
c          sum=0.
c          do  j=1,np
c            k1=3
cc            pt2=(pcm(1,j)**2+pcm(2,j)**2)
cc            scal=1.+(scal-1.)/max(1.,pt2)       !large pt are peaked around eta=0
cc            if(pcm(3,j)**2.gt.amass(j)**2/pt2)k1=3
c            do k=k1,3  
c              pcm(k,j)=scal*pcm(k,j)
c            enddo
c            pcm(4,j)=sqrt(pcm(1,j)**2+pcm(2,j)**2+pcm(3,j)**2
c     *           +amass(j)**2)
c            sum=sum+pcm(4,j)
c          enddo
c          scal=esoll/sum
c          !write(6,*)'ipass,scal,e,esoll:'
c          !    $         ,ipass,scal,sum,esoll
c          if(abs(scal-1.).le.errlim) goto 301
c        enddo
c 301    continue
c        do j=1,4
c          pa(j)=0.d0
c        enddo
c        do i=1,np
c          do j=1,4
c            pa(j)=pa(j)+pcm(j,i)
c          enddo
c          v=pcm(3,i)/pcm(4,i) 
c          Y4=0.5 * log ((1+v)/(1-v))
c          ylng(i)=ylng(i)+Y4-Y3(i)  !should be equal Y4-Y1
c        enddo
c        pa(5)=sqrt(pa(4)**2-pa(3)**2-pa(2)**2-pa(1)**2)
c        !write(6,'(a,5e11.3)')' rescaling ',pa
c      endif
      if(yco.gt.0.) then
        errlim=0.0001
        tecm=tecmx
        niter=0
 611    energ=0.
        niter=niter+1
        do i=1,np
          ylng(i)=(2*rangen()-1)*yco
          be(3)=sinh(ylng(i))
          be(4)=cosh(ylng(i))
          energ=energ+be(4)*pcm(4,i)-be(3)*pcm(3,i)
        enddo
        if(abs(energ-tecm).gt.0.1.and.niter.lt.1000)then
          goto 611
        elseif(niter.ge.1000)then
          if(ish.ge.1)write(ifch,*)'Long Flow failed:'
     &                             ,energ,tecm
          yco=0
          tecm=tecmxxx
          goto 400
        endif
                !print*,'===== energy after flow boosts',energ,'   soll: ',tecm
        do j=1,4
          pe(j)=0.
        enddo
        do i=1,np
          be(1)= 0
          be(2)= 0
          be(3)= sinh(ylng(i))
          be(4)= cosh(ylng(i))
          call utlob3(1,be(1),be(2),be(3),be(4),1e0
     *         , pcm(1,i), pcm(2,i), pcm(3,i), pcm(4,i))
          do j=1,4
          pe(j)=pe(j)+pcm(j,i)
          enddo
        enddo
        pe(5)=sqrt(pe(4)**2-pe(3)**2-pe(2)**2-pe(1)**2)
       !write(6,'(a,5e11.3)')'flow boosts',pe
        do j=1,4
          pa(j)=0.
        enddo
        do i=1,np
          call utlob3(1,pe(1),pe(2),pe(3),pe(4),pe(5)
     *         , pcm(1,i), pcm(2,i), pcm(3,i), pcm(4,i))
          do j=1,4
            pa(j)=pa(j)+pcm(j,i)
          enddo
        enddo
        pa(5)=sqrt(pa(4)**2-pa(3)**2-pa(2)**2-pa(1)**2)
                !write(6,'(a,5e11.3)')' cms boost ',pa
        esoll=tecm
        scal=1.
        do ipass=1,200
          sum=0.
          do  j=1,np
            do k=3,3
              pcm(k,j)=scal*pcm(k,j)
            enddo
            pcm(4,j)=sqrt(pcm(1,j)**2+pcm(2,j)**2+pcm(3,j)**2
     *           +amass(j)**2)
            sum=sum+pcm(4,j)
          enddo
          scal=esoll/sum
          !write(6,*)'ipass,scal,e,esoll:'
          !    $         ,ipass,scal,sum,esoll
          if(abs(scal-1.).le.errlim) goto301
        enddo
 301    continue
        do j=1,4
          pa(j)=0.
        enddo
        do i=1,np
          do j=1,4
            pa(j)=pa(j)+pcm(j,i)
          enddo
        enddo
        pa(5)=sqrt(pa(4)**2-pa(3)**2-pa(2)**2-pa(1)**2)
        !write(6,'(a,5e11.3)')' rescaling ',pa
      endif

 400  continue

    !~~~~~~~~~~radial flow -> particles~~~~~~~~~~~~~~~~~~

      ipo=ip
      if(iorptl(ipo).ne.0)ipo=iorptl(ipo)
      if(yrmax.gt.0.d0) then
        fecc=0.
        if(ianisotropy.eq.1)then  !ckw
          xx=uptl(ipo)   ! <x**2>
          yy=optl(ipo)   ! <y**2>
          xy=desptl(ipo) ! <x*y>
          dta=0.5*abs(xx-yy)
          ev1=(xx+yy)/2+sqrt(dta**2+xy**2)
          ev2=(xx+yy)/2-sqrt(dta**2+xy**2)
          yy=ev1
          xx=ev2
          if(abs(xx+yy).gt.0.)then
            ecc=(yy-xx)/(yy+xx)
            fecc=min(facecc,ecc) !be careful : fecc change <pt> since it is the elliptical deformation of the sub cluster(give strength of v2)
c          print *,fecc,xx,yy,idptl(ip)
            phiclu=phievt
            aa=cos(phiclu)
            bb=sin(phiclu)
            cc=-sin(phiclu)
            dd=cos(phiclu)
          endif
        endif
        if(fecc.lt.1e-5)then
          fecc=0
          aa=1
          bb=0
          cc=0
          dd=1
        endif
        errlim=0.0001
        tecm=tecmxx
        niter=0
 610    energ=0.
        niter=niter+1
        do i=1,np
          do j=1,4
            pa(j)=pcm(j,i)
          enddo
          yrad(i)=sqrt(rangen())
          phirad(i)=2.*pi*rangen()
          pt2=(pcm(1,i)**2+pcm(2,i)**2) !+amass(i)**2)
          bex=dsinh(dble(yrad(i)*yrmax))*cos(phirad(i))
     *       *(1+fecc/(1.+pt2))
          bey=dsinh(dble(yrad(i)*yrmax))*sin(phirad(i))
     *       *(1-fecc/(1.+pt2))
          ! 
          ! **** attention **** the 4velocity be is the velocity of the 
          !                     lab seen from the flowing fluid element
          be(1)=aa*bex+cc*bey
          be(2)=bb*bex+dd*bey
          be(3)=0d0
          be(4)=sqrt(1+be(1)**2+be(2)**2)
c          bp=0d0
c          do k=1,3
c            bp=bp-pcm(k,i)*be(k)
c          enddo
c         en=be(4)*pcm(4,i)+bp
          bet=-(pa(1)*be(1)+pa(2)*be(2))/(be(4)+1.)
c         bet=bet+pcm(4,i) !would be the correct Lorentz trafo
          aemass=max(sqrt(pt2),min(ydsrd,amass(i)))!min(ydsrd,sqrt(pt2+amass(i)**2)) !sqrt(pt2+amass(i)**2) !min(1.,amass(i)) !min(1.,sqrt(pt2+amass(i)**2)) !max(sqrt(pt2),amass(i)) !max(0.,amass(i)-0.2)**0.0002
          if(abs(ident(i)).eq.17)aemass=0      !to be checked for deuterium 
          bet=bet+aemass !sqrt(pt2+aemass**2) !simplified version without p(3) and limited mass
          pa(1)=pa(1)-bet*be(1)
          pa(2)=pa(2)-bet*be(2)
          en=sqrt(pa(1)**2+pa(2)**2+pa(3)**2+amass(i)**2)
cc        print *,i,pa(3),en,pa(4),en/pa(4)
          energ=energ+en
        enddo
        if(energ-tecm.gt.0..and.niter.lt.1000)then
          goto 610
        elseif(niter.ge.1000)then
          if(ish.ge.1)write(ifch,*)'Radial Flow failed:'
c          print *,'Radial Flow failed:'
     &                             ,yrmax,energ,tecm
          iret=3
          goto 1000
        endif
        do j=1,4
          pa(j)=0.d0
        enddo
        do i=1,np
          pt2=(pcm(1,i)**2+pcm(2,i)**2)!+amass(i)**2)
          bex=dsinh(dble(yrad(i)*yrmax))*cos(phirad(i))
     *       *(1+fecc/(1.+pt2))
          bey=dsinh(dble(yrad(i)*yrmax))*sin(phirad(i))
     *       *(1-fecc/(1.+pt2))
          ! 
          ! **** attention **** the 4velocity be is the velocity of the 
          !                     lab seen from the flowing fluid element
          be(1)=aa*bex+cc*bey
          be(2)=bb*bex+dd*bey
          be(3)=0d0
          be(4)=sqrt(1+be(1)**2+be(2)**2)
c mimic boost transformation but protect against too high values of p(3) (p(3)~p(4)) 
            bet=-(pcm(1,i)*be(1)+pcm(2,i)*be(2))/(be(4)+1.)

            aemass=max(sqrt(pt2),min(ydsrd,amass(i)))!min(ydsrd,sqrt(pt2+amass(i)**2)) !min(1.,amass(i)) !min(1.,sqrt(pt2+amass(i)**2))  !max(sqrt(pt2),amass(i))!sqrt(pt2+amass(i)**2) !max(sqrt(pt2),amass(i))!max(0.,amass(i)-0.2)**0.0002
            if(abs(ident(i)).eq.17)aemass=0. !to be checked for deuterium 
            bet=bet+aemass !sqrt(pt2+aemass**2) !simplified version without p(3) and limited mass
            pcm(1,i)=pcm(1,i)-bet*be(1)
            pcm(2,i)=pcm(2,i)-bet*be(2)
c          call utlob3(1,be(1),be(2),be(3),be(4),1e0
c     *         , pcm(1,i), pcm(2,i), pcm(3,i), pcm(4,i))
          do j=1,2
            pa(j)=pa(j)+pcm(j,i)
          enddo
        enddo
c rescale momenta to have pt conservation
        p1shift=-sngl(pa(1)/dble(np))
        p2shift=-sngl(pa(2)/dble(np))
        do j=1,np
          pcm(1,j)=pcm(1,j)+p1shift
          pcm(2,j)=pcm(2,j)+p2shift
          pcm(4,j)=sqrt(pcm(1,j)**2+pcm(2,j)**2
     &         +pcm(3,j)**2+amass(j)**2)
        enddo
c rescale all momenta to conserve initial energy
        esoll=tecm
        scal=1.
        do ipass=1,500
          sum=0.
          do  j=1,np
            k1=3
            pt2=max(float(min(0,5-np)),pcm(1,j)**2+pcm(2,j)**2)
            scal=1.+(scal-1.)!     !large pt are peaked around eta=0
c            if(pcm(3,j)**2.gt.amass(j)**2/pt2)k1=3
            do k=k1,3  
              pcm(k,j)=scal*pcm(k,j)
            enddo
            pcm(4,j)=sqrt(pcm(1,j)**2+pcm(2,j)**2+pcm(3,j)**2
     *           +amass(j)**2)
            sum=sum+pcm(4,j)
          enddo
          scal=esoll/sum
c          write(6,*)'ipass,scal,e,esoll:'
c     $         ,ipass,scal,sum,esoll,p1shift,p2shift
          if(abs(scal-1.).le.errlim) goto 300
        enddo
 300    continue
      else !yrmax.le.0.d0 (no rad flow)
        do n=1,np
          yrad(n)=0.
          phirad(n)=0.
        enddo
      endif

    !~~~~~~~~end radial flow stuff ~~~~~~~~~~~~~~~~~~~~~~~~

      rini=0.
      if(jcase.eq.9)then
        ityip=ityip+1
        rini=3./5.*((3.*radptl(ip)/4./pi)**0.3333)**2 !<r**2>=3/5*R**2 for sphere of radius R
      else
        ityip=60
        xx=uptl(ipo)            ! <x**2>
        yy=optl(ipo)            ! <y**2>
        rini=sqrt(5./3.*(xx+yy)) !<r**2>=3/5*R**2 for sphere of radius R
      endif
      do n=1,np
        nptl=nptl+1
        if(nptl.gt.mxptl)call utstop('hnbaaa: mxptl too small&')
        idptl(nptl)=ident(n)
        do j=1,4
          p(j)=pcm(j,n)
        enddo
        p(5)=amass(n)
        call utlob2(-1,c(1),c(2),c(3),c(4),c(5),p(1),p(2),p(3),p(4),10)
        call setpptl(nptl
     .  ,sngl(p(1)),sngl(p(2)),sngl(p(3)),sngl(p(4)),sngl(p(5)))
        call setityptl(nptl,ityip)
        call getiorptl(ip,ipo)
        call setiorptl(nptl,ip)
        call setjorptl(nptl,ipo)
        call setistptl(nptl,0)
        call getzpaptl(ipo,zpa1,zpa2)
        call setzpaptl(nptl,zpa1,zpa2)
        radptl(nptl)=0.
        dezptl(nptl)=0.
        itsptl(nptl)=0
c        if(jcase.eq.6)then !rope
c          call setrinptl(nptl,ylng(n))
c        endif 
        !protection against very high momentum particle (it can happen very very 
        !boosted cluster (which do no really make sense anyway))
        if(jcase.eq.3)then   !ckw
          if(p(4).ge.0.5*engy)then
            if(ish.ge.4)call clist('&',nptlc+1,nptl,60,60)
            nptl=nptlc
            iret=4
            if(ish.ge.4)write(ifch,*)'Decay skipped (p4 high) !',ntry
            if(ntry.lt.10)goto 10
            goto 1000
          endif
        endif
        if(jcase.eq.3.or.jcase.eq.6)then   !ckw
          call getityptl(nptl,ityn)
          if(ityn.eq.60)then
            r=1.15*rini*yrad(n) !yrad=y/ymax
            ptm=max(p(5), 0.0001d0)
c            ptm=sqrt(p(5)**2+p(1)**2+p(2)**2)
            if(tauzer2.gt.0.)then
              rndm=rangen()
              tau=-taufo*log(rndm)*p(4)/ptm
            else
              tau=taufo*p(4)/ptm
            endif
c            tau=tauzer+rangen()*taufo
cc            tau=tauzer+min(abs(ylng(n))*4.,taufo*10.)   !change pseudorapidity shape
cc            tau=tauzer+taufo*(2.25/sqrt(yrad(n)**2+0.04)-1.)   !do not give good result with hacas (probably too large)
c            !print *,tau,taufo
c            !z=xorptl(3,ipo)
c            !t=xorptl(4,ipo)
c            !zeta=0.5*log((t+z)/(t-z))-0.5*delzet+2*0.5*delzet*rangen()
c            zeta=sign(ainfin,sngl(p(3)))
c            if((p(4)-abs(p(3))).gt.1d-12)then
c              zeta=0.5*log((p(4)+p(3))/(p(4)-p(3)))
c            endif
c            z=tau*sinh(zeta)
c            t=tau*cosh(zeta)
            bex=cos(phirad(n))
            bey=sin(phirad(n))
            be(1)=-(aa*bex+cc*bey)
            be(2)=-(bb*bex+dd*bey)
            call getxorptl(ipo,xo1,xo2,xo3,xo4)
c            call setxorptl(nptl
c     .           ,xo1+r*be(1),xo2+r*be(2),z,t)
c            call updateXor(ip,nptl)
c            call getxorptl(nptl,xo1,xo2,xo3,xo4)
            xo1=xo1+p(1)/p(4)*tau
            xo2=xo2+p(2)/p(4)*tau
            xo3=xo3+p(3)/p(4)*tau
            xo4=xo4+p(4)/p(4)*tau
            call setxorptl(nptl
     .           ,xo1+r*be(1),xo2+r*be(2),xo3,xo4)
            tivptl(1,nptl)=xorptl(4,nptl)
            call updateXor(1,nptl)     !just for tivptl(2)
          else                  !distribute position to avoid problems with UrQMD
            yr=sqrt(rangen())
            phir=2.*pi*rangen()
            r=1.15*rini*yr
            tau=2.25/sqrt(yr**2+0.04)-0.75
            zeta=sign(ainfin,sngl(p(3)))
            if((p(4)-abs(p(3))).gt.1d-12)then
              zeta=0.5*log((p(4)+p(3))/(p(4)-p(3)))
            endif
            z=tau*sinh(zeta)
            t=tau*cosh(zeta)
            call getxorptl(ipo,xo1,xo2,xo3,xo4)
            call setxorptl(nptl,xo1+r*cos(phir),xo2+r*sin(phir),z,t)
          endif
c          call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
c          r=rangen()
c          call settivptl(nptl,t,t+taugm*(-alog(r)))
        elseif(jcase.eq.9)then  !same as in ProReF
          yr=sqrt(rangen())
          phi=2.*pi*rangen()
          radius=1.15*rini*yr
          bex=cos(phi)
          bey=sin(phi)
          call updateXor(nptl,nptl)
          call getxorptl(nptl,xo1,xo2,xo3,xo4)
          call setxorptl(nptl,xo1-radius*bex,xo2-radius*bey,xo3,xo4)
        endif
      enddo

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      enddo         !nc loop
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if(ish.ge.4)then
        write(ifch,*)'decay products:'
        call clist('&',nptlc+1,nptl,40,60)
      endif
      
c check
      do j=1,4
        pall(j)=0d0
      enddo
      do i=nptlc+1,nptl
        call getpptl(i,piii(1),piii(2),piii(3),piii(4),piii(5))
        do j=1,4
          pall(j)=pall(j)+dble(piii(j))
        enddo
      enddo
      am=(pall(4)-pall(3))*(pall(4)+pall(3))-pall(2)**2-pall(1)**2
      IF (am.ge.0.d0) THEN
        am=sqrt(am)
      ELSE
        am=0.d0    ! problems!!!
      END IF
      if(ish.ge.1.and.abs(sngl(am)-tecmor)/tecmor.gt.1e-1)
     .write(ifmt,'(a,3e11.3,i3,3x,3i4)')
     .'hnbaaa (2) energy pb',tecm,tecmor,am,jcase,keu,ked,kes

 1000 continue
      call utprix('hnbaaa',ish,ishini,4)
      return
      end


c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################
c#########                                                                  #########
c#########                     hnb  routines                                #########
c#########                                                                  #########
c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################

c-----------------------------------------------------------------------
      subroutine hnbcreate
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter(mxpairst=145000)
      common/chnbcreate/ihnbcreate
      parameter (literm=500)
      parameter (mspecs=400)
      parameter (mxpair=mspecs**2*4)  
      parameter (mxlk=375)
      parameter (mxidx=3**6)
      if(ihnbcreate.eq.1)return
      ihnbcreate=1
      call memo(1,'create hnb objects (1) ;')
      call wgtpairstcreate(mxpairst)
      call memo(2,';')
      call memo(1,'create hnb objects (2) ;')
      call idpairstcreate(2,mxpairst)
      call memo(2,';')
      call memo(1,'create hnb objects (3) ;')
      call lkfokcreate(mxlk,7,7,7,7)
      call memo(2,';')
      call hnbspd(iospec)
      call hnbspd2
      call hnbpajini
      call memo(1,'create hnb objects (4) ;')
      call lspecscreate(literm,nspecs)  
      call memo(2,';')
      end

c-----------------------------------------------------------------------
      subroutine hnbdestroy
c-----------------------------------------------------------------------
      include 'epos.inc'
      common/chnbcreate/ihnbcreate
      if(ihnbcreate.eq.0)return
      ihnbcreate=0
      call memo(1,'destroy hnb objects ;')
      call lspecsdestroy()
      call wgtpairstdestroy()
      call idpairstdestroy()
      call lkfokdestroy()
      call memo(2,';')
      end

c------------------------------------------------------------------------------
      subroutine hnbcor(mode)
c------------------------------------------------------------------------------
c determines(mode=1) and plots (mode=2) two particle  correlations
c for the configurations /confg/
c------------------------------------------------------------------------------
      include 'epos.inc'
      integer bns
      parameter (maxp=6000,bns=100)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      dimension zwei(bns),zz(bns)!,phi(bns),yy(bns)
      common/cor/wert(bns),cwert(bns)
      character*6 cen,cvol

           if(mode.eq.1)then

      nctcor=nctcor+1

      if(nctcor.eq.1)then
      do nn=1,bns
      wert(nn)=0
      cwert(nn)=0
      enddo
      endif

      ll=0

      do ii=1,np-1
      do jj=ii+1,np

      ll=ll+1
      prod=0

      do kk=1,3
      prod=prod+pcm(kk,ii)*pcm(kk,jj)
      enddo

      cs=prod/pcm(5,ii)/pcm(5,jj)

      if(abs(cs).gt.1.)then
      cs=aint(cs)
      ang=acos(cs)
      else
      ang=acos(cs)
      endif

      if(cs.eq.1.)then
      nk=bns
      nw=1
      elseif(ang.eq.pi)then
      nk=1
      nw=bns
      else
      nw=1+aint(ang/pi*bns)
      nk=1+aint((cs+1.)/2.*bns)
      endif
      nw=min(nw,bns)
      nk=min(nk,bns)

      wert(nw)=wert(nw)+1
      cwert(nk)=cwert(nk)+1

      enddo
      enddo

           elseif(mode.eq.2)then

      do mm=1,bns
c      phi(mm)=.5*pi/bns+(mm-1)*pi/bns
      zwei(mm)=.5*2./bns+(mm-1)*2./bns-1.
c      yy(mm)=wert(mm)/nctcor
      zz(mm)=cwert(mm)/nctcor
      enddo

      write(cen,'(f6.1)')tecm
      write(cvol,'(f6.1)')volu

      write(ifhi,'(a)')    'newpage zone 1 1 1 openhisto'
      write(ifhi,'(a)')    'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')    'xrange -1 1'
      write(ifhi,'(a)')    'text 0 0 "xaxis cosine"'
      write(ifhi,'(a)')    'text 0 0 "yaxis counts"'
      write(ifhi,'(a)')    'text 0.4 0.91 "V='//cvol//'"'
      write(ifhi,'(a)')    'text 0.15 0.91 "E='//cen//'"'
      write(ifhi,'(a)')    'array 2'
         do mm=1,bns
      write(ifhi,'(2e13.5)')zwei(mm),zz(mm)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

           endif

      return
      end

c----------------------------------------------------------------------
      subroutine hnbfac(faclog)
c----------------------------------------------------------------------
c  returns log of factor for phase space weight
c  faclog= log{ prod[ m_i*(2*s_i+1)*volu/4/pi**3/hquer**3/(n_l+1-i) ] }
c      ~~~~~~~~~~~~~~
c  corresponds to eq. 67 of micro paper :
c         Cvol * Cdeg * Cident * Cmicro
c    the factors partly compensate each other !!
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
c      integer ii(maxp)
      common /clatt/nlattc,npmax

      faclog=0

c sum_i log m_i*g_i*volu/4/pi**3/hquer**3/(n_l+1-i) -> flog
      flog=0
      do i=1,np
      call hnbfaf(i,gg,am)
      flog=flog+alog(gg*am*volu/4/pi**3/hquer**3/(nlattc+1-i))
      enddo
      faclog=faclog+flog

      return
      end

c----------------------------------------------------------------------
      subroutine hnbfaf(i,gg,am)
c----------------------------------------------------------------------
c  returns degeneracy gg and mass am  for factor f5
c----------------------------------------------------------------------
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/metr1/iospec,iocova,iopair,iozero,ioflac,iomom
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      common/drop6/tecm,volu

      hquer=0.197327
      cc=0.216416
      dd=13.773935

ckw   call hnbspi(ident(i),spideg)
      do ii=1,nspecs
        if(ident(i).eq.ispecs(ii))then
          gg=gspecs(ii)
          goto 777 
        endif
      enddo
      !print*,i,ident(i)
      stop'####### ERROR 18042018 #######'
 777  continue
      !print*,'hnbfaf',i,ident(i),gg

      am=0.5                ! 1 / 2     (no dimension)

      return
      end

c----------------------------------------------------------------------
      subroutine hnbiiw(x,f,df)
c----------------------------------------------------------------------
c returns fctn value and first derivative at x of the
c i-th integrated weight fctn minus random number
c for the asympotic phase space integral.
c input:
c   x:   x-value
c   iii: i-value (via common/ciiw/iii,rrr)
c   rrr: random number   ( " )
c output:
c   f:   fctn value
c   df:  first derivative
c----------------------------------------------------------------------
      common/ciiw/iii,rrr
      i=iii
      f=x**(2*i-2)*(i-(i-1)*x**2)-rrr
      df=2*i*(i-1)*(x**(2*i-3)-x**(2*i-1))
      return
      end

c----------------------------------------------------------------------
      subroutine hnbinit(iret)
c----------------------------------------------------------------------
c  generates initial configuration
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cspecs2/pspecs(0:mspecs+1),yie2volu,yiespecs
      common/crnoz/rnoz(maxp-1)
      common/crnoy/rnoy(3,maxp-1)
      common/citer/iter,itermx
      common/cfact/faclog
      common/chnbin/nump,ihadro(maxp)
      common /clatt/nlattc,npmax
      parameter(maxit=500000)
      common/count/nacc,nrej,naccit(maxit),nptot,npit(maxit)
      common/ctaue/taue
      if(ish.ge.7)write(ifch,*)('-',i=1,10)
     *,' entry sr hnbinit ',('-',i=1,30)

      iret=0
 
      iter=0

      ktnbod=0

      nlattc=8*(tecm/10)*(1/(tecm/volu))**0.2*(nspecs/3.)**0.3
      if(aspecs(1).lt.0.010)nlattc=nlattc*3
      nlattc=max(nlattc,20)
      if(iternc.lt.0)iternc=1.500*nlattc

      itermx=iterma
      if(itermx.le.0)then
        e=tecm/volu
        b=1.1*(e+0.33)**0.66
        a=13.*(e+0.13)**(-0.65)
        tm=34.*(e+0.65)**(-0.61)
        t=a+b*volu
        taue=max(t,tm)
        itermx=(-itermx)*taue
      else
        taue=0
      endif
      if(ish.ge.5)write(ifch,*)'itermx:',itermx

      if(iternc.gt.itermx/2)iternc=itermx/2

      ntry=0
 77   continue
      ntry=ntry+1

      keuz=keu
      kedz=ked
      kesz=kes
      kecz=kec
      iret1=0
      iret2=0
      if(ioinco.eq.0)then
            call hnbmin(iospec,keuz,kedz,kesz,kecz)
            if(iograc.eq.1)call hgcaaa
      elseif(ioinco.ge.1)then
            nk=keu+ked+ked+kec
            if(tecm.lt.1.5.and.nk.eq.0)then
                   call hnbmin(iospec,keuz,kedz,kesz,kecz)
            elseif(tecm.lt.2.0.and.nk.ne.0)then
                   call hnbmin(iospec,keuz,kedz,kesz,kecz)
            else
                   iret1=1
                   if(ntry.le.3)then
                     call hgcaaa
                     call hgcnbi(iret2)
                   else
                     iret2=1
                   endif
                   if(iret2.eq.1)then
                     call hnbmin(iospec,keuz,kedz,kesz,kecz)
                     iret=0
                     if(ish.ge.5)then
                       write(ifch,*)'hadron set from hnbmin:'
                       write(ifch,'(10i6)')(ihadro(k),k=1,nump)
                     endif
                   endif
           endif
      endif

      np=nump+nadd
      if(np.gt.maxp)then
        print*,'ERROR for tecm =',tecm 
        stop'ERROR np too large'
      endif
 
      nlattc=max(nlattc,1+int(np*1.2))
c      print *,np,nlattc
      if(nlattc-1.gt.maxp)then
        print*,'ERROR for tecm =',tecm 
        stop'ERROR maxp too small'
      endif
 
      do i= 1, nlattc-1
      rnoz(i)=rangen()
      enddo

      do i= 1, nlattc-1
      do k=1,3
      rnoy(k,i)=rangen()
      enddo
      enddo
      !call flpsore(rno,ntm2) not needed any more

      if(nadd.gt.0)then
      do i=nump+1,np
      ihadro(i)=110
      enddo
      endif

      do i=1,np
        ident(i)=ihadro(i)
        amass(i)=-1
        do j=1,nspecs
          if(ident(i).eq.ispecs(j))then
            amass(i)=aspecs(j)
            goto1
          endif
        enddo
    1   continue
        if(amass(i).lt.0.)then
          if(ntry.lt.10)goto 77           
          write(ifch,*)'***** hnbinit'
          write(ifch,*)'***** invalid particle species '
          write(ifch,*)'***** iocova =',iocova,'   ntry=',ntry
          write(ifch,*)'***** droplet mass:',tecm,'   ioinco =',ioinco 
          write(ifch,*)'***** flavour:'
          write(ifch,*)'*****',keu,ked,kes,kec,keb,ket
          write(ifch,*)'***** np =',np,'   nump =',nump,'   nadd =',nadd
          write(ifch,*)'***** list of hadrons (min hadron set):'
          write(ifch,*)(ihadro(ii),ii=1,nump)
          write(ifch,*)'***** exit with iret = 1 (invalid event) A'
          iret=1
          goto1000
        endif
      enddo
      if(yiespecs.le.iocova)call hnbody    !covariant
      if(yiespecs.gt.iocova)call hnbodz    !noncovariant
      call hnbfac(faclog)
      wtlog=wtxlog+faclog

      iret=0
      if(wtlog.le.-0.99999E+35)then
          if(ntry.lt.10)goto 77 
cc        if(ish.ge.1) then
          write(ifch,*)'***** hnbinit'
          write(ifch,*)'***** covariant : ',yiespecs.le.iocova
          write(ifch,*)'***** wtlog for initl config < -1E+35'
          write(ifch,*)'***** wtlog,faclog:',wtlog,faclog
          write(ifch,*)'***** iret1,2:',iret1,iret2
          write(ifch,*)'***** droplet mass:',tecm
          write(ifch,*)'***** flavour:'
          write(ifch,*)'*****',keu,ked,kes,kec,keb,ket
          write(ifch,*)'***** list of hadrons (min hadron set):'
          write(ifch,*)(ihadro(i),i=1,nump)
          write(ifch,*)'***** exit with iret = 1 (invalid event) B'
cc        endif
        iret=1
        goto1000
      endif

      if(ish.ge.7)then
        write(ifch,*)'initial configuration:'
        call xhnbwri
      endif

      itermx=iterma
      if(itermx.le.0)then
      e=tecm/volu
      b=1.1*(e+0.33)**0.66
      a=13.*(e+0.13)**(-0.65)
      tm=34.*(e+0.65)**(-0.61)
      t=a+b*volu
      taue=max(t,tm)
      itermx=(-itermx)*taue
      else
      taue=0
      endif
      if(ish.ge.5)write(ifch,*)'itermx:',itermx

      if(iternc.gt.itermx/2)iternc=itermx/2

      nacc=0
      nrej=0

 1000 continue

      if(ish.ge.7)write(ifch,*)('-',i=1,30)
     *,' exit sr hnbinit ',('-',i=1,10)

      return
      end

cc----------------------------------------------------------------------
cc already commented subroutine hnbint(tecmx,nevtxx,nsho) removed in 3241
cc----------------------------------------------------------------------

cc----------------------------------------------------------------------
      subroutine hnbmet
c----------------------------------------------------------------------
c  change (or not) configuration via metropolis
c  configuration=np,tecm,amass(),ident(),pcm(),volu,wtlog
c    (common /confg/)
c  nlattc (in /clatt/) must be set before calling this routine
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      common/crnoz/rnoz(maxp-1)
      common/crnoy/rnoy(3,maxp-1)
      real rnozo(maxp-1)
      real rnoyo(3,maxp-1)
      common/cfact/faclog
      dimension amasso(maxp),idento(maxp),pcmo(5,maxp)
      integer jc1(nflav,2),jc2(nflav,2)
     .,jc12old(nflav,2),jc12new(nflav,2),jc12miss(nflav,2)
     .,jc34old(nflav,2),jc34miss(nflav,2)
      common/citer/iter,itermx
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cspecs2/pspecs(0:mspecs+1),yie2volu,yiespecs
      parameter (literm=500)
      common/cmet/kspecs(mspecs),liter
     *,iterl(literm),iterc(literm)
      parameter (mxpair=mspecs**2*4)
      common /clatt/nlattc,npmax
      parameter (nhise=1) !(nhise=100)
      common/chise/hise(mspecs,nhise)
      integer id3old,id4old!,id3new,id4new
      parameter(maxit=500000)
      common/count/nacc,nrej,naccit(maxit),nptot,npit(maxit)
      if(ish.ge.7)then
        write(ifch,*)('-',i=1,10)
     *,' entry sr hnbmet ',('-',i=1,30)
        write(ifch,'(1x,a,i4)')'iteration:',iter
      endif
      ielast=1 !inelastic
      !if(rangen().lt.felamet)ielast=2 
      if( iter .gt. itermx * fitermet * felamet )ielast=2 
      if(mod(iter,iterpr).eq.0)write(ifmt,'(a,i9,4i5,2i8,f11.1)')
     .        'iteration:',iter,np,maxp,nlattc
     .        ,ielast,nacc,nrej,yiespecs
      if(maxp.gt.np)then
      do n=np+1,maxp
      ident(n)=0
      enddo
      endif

c     for iter=1
c     ----------
           if(iter.eq.1)then
      liter=1
      do i=1,nspecs
      kspecs(i)=0
      nptot=0
      enddo
      do li=1,literm
      iterc(li)=0
      enddo
      do j=1,mspecs
      do i=1,nhise
      hise(j,i)=0
      enddo
      enddo
      call xhnbzmu(-1)
           endif

c     remember old configuration
c     --------------------------
      wtlo=wtlog
      wtlox=wtxlog
      faclo=faclog
      npo=np
      if(np-1.gt.0)then
       do i=1,np-1
        rnozo(i)=rnoz(i)
        do k=1,3
         rnoyo(k,i)=rnoy(k,i)
        enddo
       enddo
      endif
      if(np.gt.0)then
      do i=1,np
      amasso(i)=amass(i)
      idento(i)=ident(i)
      do j=1,5
      pcmo(j,i)=pcm(j,i)
      enddo
      enddo
      endif

c      if(mod(iter,5).ne.0)goto 777

c     --------------------------------------------------------------
           if(iopair.eq.2.or.iopair.eq.4)then
                 !determine 2 pairs, construct 2 new pairs, update ident
                 !(double pair method)
c     --------------------------------------------------------------

      if(ielast.eq.1)then
        !~~~~~12~~~~~ 
        call hnbpad(1,n1,n2,n3dummy,n4dummy,jc12old) 
        id1old=ident(n1)
        id2old=ident(n2)
        call specs_getp(id1old,p1old)
        call specs_getp(id2old,p2old)
        call specs_getrid(id1)
        call specs_getp(id1,p1)
        call idtr7(id1,jc1)
        call specs_getrid(id2)
        call specs_getp(id2,p2)
        call idtr7(id2,jc2)
        call jcadd(1,jc1,jc2,jc12new) !add
        ident(n1)=id1
        ident(n2)=id2
        call jcadd(-1,jc12old,jc12new,jc12miss) !subtract 
        !~~~~~34~~~~~
        ncnt=0
    2   ncnt=ncnt+1
c        if(mod(ncnt,1000000).eq.0)call utstop('hnbmet&')
        if(mod(ncnt,1000000).eq.0)write(ifmt,*)
     .                           'Hnbmet loop :',ncnt
        call hnbpad(2,n1,n2,n3,n4,jc34old)
        !call jcprint(1,jc34old)
        id3old=ident(n3)
        id4old=ident(n4)
        call jcadd(1,jc12miss,jc34old,jc34miss)
        !write(ifch,*)'<---',id3old,id4old,jc34miss
       call hnbpaj(1,jc34miss,id3,id4,p34)
        !write(ifch,*)'--->',id3,id4,p34
        if(p34.lt.1e-6)goto 2
        xab=p1*p2*p34
        call hnbpaj(2,jc34old,id3old,id4old,p34old)
        xba=p1old*p2old*p34old
        if(iopair.eq.4)then
          call specs_getp(id3old,p3old)
          call specs_getp(id4old,p4old)
          call specs_getrid(id3)
          call specs_getp(id3,p3)
          call specs_getrid(id4)
          call specs_getp(id4,p4)
          xab=p1*p2*p3*p4
          xba=p1old*p2old*p3old*p4old
        endif
        ident(n3)=id3
        ident(n4)=id4
      elseif(ielast.eq.2)then
        xab=1
        xba=1
      else
        stop'ERROR 03062018c'
      endif 
      call hnbrmz(ielast)

c     -----------------------------------------------------------
           elseif(iopair.eq.3)then !no flavor consideration, 
                                   !for testing, using iospec=18
c     -----------------------------------------------------------

      if(ielast.eq.1)then
        !flavor unchanged
        n1=min( nlattc , int(1+rangen()*nlattc) )
        io=ident(n1)
        call specs_getp(io,xba)
        call specs_getrid(id)
        ident(n1)=id
        call specs_getp(id,xab)
      elseif(ielast.eq.2)then
        xab=1
        xba=1
      else
        stop'ERROR 03062018b'
      endif 
      call hnbrmz(ielast)

c     -----------------------------------------------------
           else  !wrong option
      call utstop('hnbmet: invalid choice for iopair&')
           endif
c     -----------------------------------------------------


c     determine masses/momenta/weight of trial configuration
c     ------------------------------------------------------

c 777  continue

      if(np.ge.2)then
        do i=1,np
          amass(i)=-1
          do j=1,nspecs
            if(ident(i).eq.ispecs(j))then
            amass(i)=aspecs(j)
            goto1
            endif
          enddo
    1     continue
          if(amass(i).lt.0.)
     *    call utstop('hnbmet: invalid particle species&')
        enddo
        !ckw2018 keepr=0
        !ckw2018 call hnbolo(1000) !instead of "call hnbody" for testing
        !ckw2018 keepr=1
        !~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(yiespecs.le.iocova)call hnbody
        if(yiespecs.gt.iocova)call hnbodz
        !~~~~~~~~~~~~~~~~~~~~~~~~~~
      else
        wtxlog=-1e35
      endif
      call hnbfac(faclog)
      wtlog=wtxlog+faclog
      if(ish.ge.7)then
        write(ifch,*)'trial configuration:'
        call xhnbwri
      endif

c     accept or not trial configuration (metropolis)
c     ----------------------------------------------

      if(ish.ge.7)write(ifch,'(1x,a,4i5,a,4i5,a)')
     *'metropolis decision for '
     *,id1old,id2old,id3old,id4old,'   -->  '
     *,id1,id2,id3,id4,' :'
      iacc=0
           if(wtlog-wtlo.lt.30.)then
      q=exp(wtlog-wtlo)*xba/xab
      r=rangen()
      if(r.le.q)iacc=1
      if(ish.ge.7)write(ifch,*)'new weight / old weight:',q,'    '
     *,'random number:',r
           else
      iacc=1
      if(ish.ge.7)write(ifch,*)'log new weight / old weight:'
     *,wtlog-wtlo
           endif
           if(iacc.eq.1)then
      !print*,np
      if(ish.ge.7)write(ifch,*)'new configuration accepted'
      nacc=nacc+1
      if(iter.le.maxit)naccit(iter)=1
           else
      !print*,'               ',np
      if(ish.ge.7)write(ifch,*)'old configuration kept'
      nrej=nrej+1
      wtlog=wtlo
      wtxlog=wtlox
      faclog=faclo
      np=npo
      if(np-1.gt.0)then
       do i=1,np-1
        rnoz(i)=rnozo(i)
        do k=1,3
         rnoy(k,i)=rnoyo(k,i)
        enddo
       enddo
      endif
      if(np.gt.0)then
      do i=1,np
      amass(i)=amasso(i)
      ident(i)=idento(i)
      do j=1,5
      pcm(j,i)=pcmo(j,i)
      enddo
      enddo
      endif
           endif
           if(ioobsv.eq.0)then
      if(iter.le.maxit)npit(iter)=np
      if(iter.gt.iternc)nptot=nptot+np
      else
      npob=0
      do i=1,np
      if(ioobsv.eq.ident(i))npob=npob+1
      enddo
      npit(iter)=npob
      if(iter.gt.iternc)nptot=nptot+npob
           endif
      if(ish.ge.7)then
        write(ifch,*)'actual configuration:'
        call xhnbwri
        if(ish.eq.27)stop'change this?????????????' !call hnbcor(1)
      endif

c     printout/return
c     ---------------
      if(iosngl.ne.nrevt+1.and.iocite.ne.1)goto1000
      npmax=max(npmax,np)
           if(liter.le.literm)then
      iterc(liter)=iterc(liter)+1
      do i=1,np
      do j=1,nspecs
      if(ident(i).eq.ispecs(j))then
      call lspecsincrement(liter,j)
      goto8
      endif
      enddo
    8 continue
      enddo
      if(mod(iter,iterpl).eq.0)then
      iterl(liter)=iter
      liter=liter+1
c     if(liter.le.literm)then
c     iterc(liter)=iterc(liter-1)
c     do j=1,nspecs
c     lspecs(liter,j)=lspecs(liter-1,j)
c     enddo
c     endif
      endif
           endif
      if(iter.le.iternc)return

           do i=1,np
      call xhnbzen(i)  !fill energy histogram
      do j=1,nspecs
      if(ident(i).eq.ispecs(j))then
      kspecs(j)=kspecs(j)+1
      goto7
      endif
      enddo
    7 continue
           enddo
      call xhnbzmu(1)  !fill multiplicity histogram

           if(iter.eq.itermx.and.npmax.ge.nlattc.and.ish.ge.1)then
      call utmsg('hnbmet&')
      write(ifch,*)'*****  nlattc too small'
      write(ifch,*)'nlattc:',nlattc,'   npmax:',npmax
      call utmsgf
           endif

1000  continue
      if(ish.ge.7)then
        write(ifch,*)'accepted proposals:',nacc
     *,'  rejected proposals:',nrej
        write(ifch,*)('-',i=1,30)
     *,' exit sr hnbmet ',('-',i=1,10)
      endif
      return
      end

c----------------------------------------------------------------------
      subroutine hnbmin(iospec,keux,kedx,kesx,kecx)
c----------------------------------------------------------------------
c  returns min hadron set with given u,d,s,c content
c  input:
c     keux: net u quark number
c     kedx: net d quark number
c     kesx: net s quark number
c     kecx: net c quark number
c  output (written to /chnbin/):
c     nump: number of hadrons
c     ihadro(n): hadron id for n'th hadron
c----------------------------------------------------------------------
      !no include aaa.h to avoid conflicts
      integer      iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      parameter(maxp=6000)
      common/chnbin/nump,ihadro(maxp)
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      logical wri
      character f1*11
      wri=.false.
      if(ish.ge.7)wri=.true.
      if(wri)then
        write(ifch,*)('-',i=1,10)
     *       ,' entry sr hnbmin ',('-',i=1,30)
        write(ifch,*)'iospec,keux,kedx,kesx,kecx'
     *               ,iospec,keux,kedx,kesx,kecx
      endif

      if(iospec.eq.18)then
        nump=1
        ihadro(nump)=ispecs(1)
        return
      endif

      nump=0
      f1='(4i3,i7,i6)'
      ke=iabs(keux+kedx+kesx+kecx)

      if(keux+kedx+kesx+kecx.ge.0)then
      keu=keux
      ked=kedx
      kes=kesx
      kec=kecx
      isi=1
      else
      keu=-keux
      ked=-kedx
      kes=-kesx
      kec=-kecx
      isi=-1
      endif
      if(wri)write(ifch,'(4i3)')keux,kedx,kesx,kecx
      if(wri)write(ifch,'(4i3)')keu,ked,kes,kec

c get rid of anti-c and c (140, 240, -140, -240)
      if(kec.ne.0)then
   10 continue
      if(kec.lt.0)then
      kec=kec+1
      if(keu.gt.ked)then
      keu=keu-1
      nump=nump+1
      ihadro(nump)=140
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      else
      ked=ked-1
      nump=nump+1
      ihadro(nump)=240
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      endif
      goto10
      endif
   11 continue
      if(kec.gt.0)then
      kec=kec-1
      if(keu.lt.ked)then
      keu=keu+1
      nump=nump+1
      ihadro(nump)=-140
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      else
      ked=ked+1
      nump=nump+1
      ihadro(nump)=-240
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      endif
      goto11
      endif
      endif

c get rid of anti-s (130,230)
    5 continue
      if(kes.lt.0)then
      kes=kes+1
      if(keu.ge.ked)then
      keu=keu-1
      nump=nump+1
      ihadro(nump)=130
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      else
      ked=ked-1
      nump=nump+1
      ihadro(nump)=230
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      endif
      goto5
      endif

c get rid of anti-d (120, -230)
   6  continue
      if(ked.lt.0)then
      ked=ked+1
      if(keu.ge.kes)then
      keu=keu-1
      nump=nump+1
      ihadro(nump)=120
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      else
      kes=kes-1
      nump=nump+1
      ihadro(nump)=-230
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      endif
      goto6
      endif

c get rid of anti-u (-120, -130)
    7 continue
      if(keu.lt.0)then
      keu=keu+1
      if(ked.ge.kes)then
      ked=ked-1
      nump=nump+1
      ihadro(nump)=-120
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      else
      kes=kes-1
      nump=nump+1
      ihadro(nump)=-130
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      endif
      goto7
      endif

      if(keu+ked+kes+kec.ne.ke)call utstop('hnbmin: sum_kei /= ke&')

      keq=keu+ked

c get rid of s (3331, x330, xx30)
      i=4
    2 i=i-1
    3 continue
      if((4-i)*kes.gt.(i-1)*keq)then
      kes=kes-i
      keq=keq-3+i
      nump=nump+1
      if(i.eq.3)ihadro(nump)=3331
      if(i.eq.2)ihadro(nump)=0330
      if(i.eq.1)ihadro(nump)=0030
           if(i.lt.3)then
      do j=1,3-i
      l=1+2*rangen()
      if(keu.gt.ked)l=1
      if(keu.lt.ked)l=2
      if(l.eq.1)keu=keu-1
      if(l.eq.2)ked=ked-1
      ihadro(nump)=ihadro(nump)+l*10**(4-j)
      enddo
           endif
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      if(kes.lt.0)call utstop('hnbmin: negative kes&')
      if(keq.lt.0)call utstop('hnbmin: negative keq&')
      goto3
      endif
      if(i.gt.1)goto2

      if(keu+ked.ne.keq)call utstop('hnbmin: keu+ked /= keq&')

c get rid of d (2221, 1220, 1120)
      i=4
   12 i=i-1
   13 continue
      if((4-i)*ked.gt.(i-1)*keu)then
      ked=ked-i
      keu=keu-3+i
      if(i.eq.3)then
      nump=nump+2
      ihadro(nump)=1220
      ihadro(nump-1)=-120
      else
      nump=nump+1
      if(i.eq.2)ihadro(nump)=1220
      if(i.eq.1)ihadro(nump)=1120
      endif
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      if(ked.lt.0)call utstop('hnbmin: negative ked&')
      if(keu.lt.0)call utstop('hnbmin: negative keu&')
      goto13
      endif
      if(i.gt.1)goto12

      if(ked.ne.0)call utstop('hnbmin: ked .ne. 0&')

c get rid of u (1111)
    9 continue
      if(keu.gt.0)then
      keu=keu-3
      nump=nump+2
      ihadro(nump)=1120
      ihadro(nump-1)=120
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      if(keu.lt.0)call utstop('hnbmin: negative keu&')
      goto9
      endif

      if(keu.ne.0)call utstop('hnbmin: keu .ne. 0&')

      if(isi.eq.-1)then
      do i=1,nump
      ihadro(i)=isi*ihadro(i)
      enddo
      endif

      do lo=1,2
      if(nump.lt.2)then
      nump=nump+1
      ihadro(nump)=110
      if(wri)write(ifch,f1)keu,ked,kes,kec,nump,ihadro(nump)
      endif
      enddo

      if(wri)write(ifch,*)('-',i=1,30)
     *,' exit sr hnbmin ',('-',i=1,10)
      return
      end

c-------------------------------------------------------------
      subroutine hnbody
c-------------------------------------------------------------
c KW 2018/05  Small but crucial modification: 
c Whereas sofar Lorentz-invariant phase space (LIPS) was 
c considered (see below), we add a factor 
c           \prod{ 2*E_i }
c which amounts to considering non-relativistic phase space
c (NRPS). The event generation is the same, but the weight 
c changes (by this factor). 
c-------------------------------------------------------------
c   Formerly subr genbod from genlib (cernlib).
c   Modified by K. Werner, march 94.
c   Subroutine to generate n-body event
c   according to Lorentz-invariant phase space (LIPS).
c   Adapted from fowl (cern w505) sept. 1974 by f. james.
c   Events are generated in their own center-of-mass.
c-------------------------------------------------------------
c   input to and output from subr thru common block config.
c   input:
c             np=number of outgoing particles
c             tecm=total energy in center-of-mass
c             amass(i)=mass of ith outgoing particle
c   output:
c             pcm(1,i)=x-momentum if ith particle
c             pcm(2,i)=y-momentum if ith particle
c             pcm(3,i)=z-momentum if ith particle
c             pcm(4,i)=energy of ith particle
c             pcm(5,i)=momentum of ith particle
c             wtxlog=log of weight of event
c--------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      dimension emm(maxp)
      common/crnoy/rnoy(3,maxp-1)
c     !pcm1 is linear equiv. of pcm to avoid double indices
      dimension em(maxp),pd(maxp),ems(maxp),sm(maxp)
     *,pcm1(5*maxp)
      dimension ri(maxp),xi(maxp),zi(maxp)
      common/cffq/ffqlog(maxp)
      common/ciiw/iii,rrr
      equivalence (nt,np),(amass(1),em(1)),(pcm1(1),pcm(1,1))
      logical wri
      data twopi/6.2831853073/
      external hnbiiw
ctp060829      nas=5 !must be at least 3
      wri=.false.
      if(ish.ge.7)wri=.true.
      if(wri)then
        write(ifch,*)('-',i=1,10)
     *,' entry sr hnbody ',('-',i=1,30)
        write(ifch,1200)np,tecm
        write(ifch,*)'particle masses:'
        write(ifch,'(1x,10f9.3)')(amass(n),n=1,np)
      endif

c..... initialization

      ktnbod=ktnbod + 1
      if(ktnbod.le.1)then
        !... ffq(n) = pi * (twopi)**(n-2) / (n-2)!
        ffqlog(1)=-1e35
        ffqlog(2)=alog(pi)
        do n=3,maxp
        ffqlog(n)=ffqlog(n-1)+log(twopi/(n-2))
        enddo
      endif

      if(nt.lt.2) goto 1001
      if(nt.gt.maxp) goto 1002
      ntm1=nt-1
      ntm2=nt-2
      ntnm4=3*nt - 4
      emm(1)=em(1)
      tm=0.0
      do 2 i=1,nt
      ems(i)=em(i)**2
      tm=tm+em(i)
    2 sm(i)=tm
      tecmtm=tecm-tm
      if(tecmtm.le.0.0) goto 1000
      emm(nt)=tecm
      wtmlog=alog(tecmtm)*ntm2 + ffqlog(nt) - alog(tecm)

c...special cases

      if(ntm2.lt.0)then
       goto 9
      elseif(ntm2.eq.0)then
       goto 5
      endif

c...calculate z_i distributed as i*z*(i-1)
      do i= 1, ntm2
        ri(i)= rnoy(1,i+1)
        zi(i)=ri(i)**(1./i)
      enddo

c...calculate x_i
      xi(ntm2)=zi(ntm2)
      do i=ntm2-1,1,-1
      xi(i)=zi(i)*xi(i+1)
      enddo

c...calculate emm().......M_i

      do 6 j=2,ntm1
    6 emm(j)=xi(j-1)*tecmtm+sm(j) !rnoy(1,j=1) not used

c...calculate wtlog

    5 continue
      wtxlog=wtmlog
      ir=ntm2
      do 7 i=1,ntm1
        pd(i)=hnbpdk(emm(i+1),emm(i),em(i+1))
        if(pd(i).gt.0.)then
          pdlog=alog(pd(i))
        else
          pdlog=-1e35
        endif
        wtxlog=wtxlog+pdlog
    7 continue

c...complete specification of event (raubold-lynch method)

      pcm(1,1)=0.0
      pcm(2,1)=pd(1)
      pcm(3,1)=0.0
      do i=2,nt
        pcm(1,i)=0.0
        pcm(2,i)=-pd(i-1)
        pcm(3,i)=0.0
        bang=twopi*rnoy(2,i-1)
        cb=cos(bang)
        sb=sin(bang)
        c=2.0*rnoy(3,i-1)-1.0
        s=sqrt(1.0-c*c)
        if(i.ne.nt)then
          esys=sqrt(pd(i)**2+emm(i)**2)
          beta=pd(i)/esys
          gama=esys/emm(i)
          do j=1,i
            ndx=5*j - 5
            aa= pcm1(ndx+1)**2 + pcm1(ndx+2)**2 + pcm1(ndx+3)**2
            pcm1(ndx+5)=sqrt(aa)
            pcm1(ndx+4)=sqrt(aa+ems(j))
            call hnbrt2(c,s,cb,sb,pcm,j)
            psave=gama*(pcm(2,j)+beta*pcm(4,j))
            pcm(2,j)=psave
          enddo
        else !(i.eq.nt)
          do j=1,i
            aa=pcm(1,j)**2 + pcm(2,j)**2 + pcm(3,j)**2
            pcm(5,j)=sqrt(aa)
            pcm(4,j)=sqrt(aa+ems(j))
            call hnbrt2(c,s,cb,sb,pcm,j)
          enddo
        endif
      enddo

c...add factor NRPS / LIPS  (sofar we considered LIPS)

      do i=1,nt
        wtxlog=wtxlog+alog(2*pcm(4,i))  ! \prod{ 2*E_i }
      enddo

c...returns

  9   continue
      goto1111

 1000 continue
      if(wri)
     *write(ifch,*)'available energy zero or negative -> wtxlog=-1e35'
      wtxlog=-1e35
      goto1111

 1001 continue
      if(wri)
     *write(ifch,*)'less than 2 outgoing particles -> wtxlog=-1e35'
      wtxlog=-1e35
      goto1111

 1002 continue
      write(ifch,*)'too many outgoing particles'
      write(ifch,1150) ktnbod
 1150 format(47h0 above error detected in hnbody at call number,i7)
      write(ifch,1200) np,tecm
 1200 format(' np:',i6/' tecm:',f10.5)
      write(ifch,*)'particle masses:'
      write(ifch,'(1x,10f6.3)')(amass(jk),jk=1,np)
      stop

1111  continue
      if(wri)write(ifch,*)('-',i=1,30)
     *,' exit sr hnbody ',('-',i=1,10)
      return
      end

c---------------------------------------------------------------------------------------------------------
      SUBROUTINE FLPSORE(A,N)
C---------------------------------------------------------------------------------------------------------
C CERN PROGLIB# M103    FLPSOR          .VERSION KERNFOR  3.15  820113
C ORIG. 29/04/78
C
C   SORT THE ONE-DIMENSIONAL FLOATING POINT ARRAY A(1),...,A(N) BY
C   INCREASING VALUES
C
C-    PROGRAM  M103  TAKEN FROM CERN PROGRAM LIBRARY,  29-APR-78
C----------------------------------------------------------------------------------------------------------
      DIMENSION A(N)
      COMMON /SLATE/ LT(20),RT(20)
      INTEGER R,RT
C
      LEVEL=1
      LT(1)=1
      RT(1)=N
   10 L=LT(LEVEL)
      R=RT(LEVEL)
      LEVEL=LEVEL-1
   20 IF(R.GT.L) GO TO 200
      IF(LEVEL.lt.0)then
       goto 50
      ELSEIF(LEVEL.eq.0)then
       goto 50
      ELSE
       goto 10
      ENDIF 
C
C   SUBDIVIDE THE INTERVAL L,R
C     L : LOWER LIMIT OF THE INTERVAL (INPUT)
C     R : UPPER LIMIT OF THE INTERVAL (INPUT)
C     J : UPPER LIMIT OF LOWER SUB-INTERVAL (OUTPUT)
C     I : LOWER LIMIT OF UPPER SUB-INTERVAL (OUTPUT)
C
  200 I=L
      J=R
      M=(L+R)/2
      X=A(M)
  220 IF(A(I).GE.X) GO TO 230
      I=I+1
      GO TO 220
  230 IF(A(J).LE.X) GO TO 231
      J=J-1
      GO TO 230
C
  231 IF(I.GT.J) GO TO 232
      W=A(I)
      A(I)=A(J)
      A(J)=W
      I=I+1
      J=J-1
      IF(I.LE.J) GO TO 220
C
  232 LEVEL=LEVEL+1
      IF((R-I).GE.(J-L)) GO TO 30
      LT(LEVEL)=L
      RT(LEVEL)=J
      L=I
      GO TO 20
   30 LT(LEVEL)=I
      RT(LEVEL)=R
      R=J
      GO TO 20
   50 continue

      do i=1,n-1
        if(a(i).gt.a(i+1))stop'FLPSORE: ERROR.                    '
      enddo

      RETURN
      END





c-------------------------------------------------------------
      subroutine hnbodz
c-------------------------------------------------------------
c   subr to generate n-body event
c   according to non-invariant phase space.
c   the phase space integral is the sum over the weights exp(wtxlog)
c   divided by the number of events.
c   ref.: hagedorn, nuov. cim. suppl ix, x (1958) 646.
c   events are generated in their own center-of-mass.
c
c   input to and output from subr is thru common block config.
c   input:
c             np=number of outgoing particles
c             tecm=total energy in center-of-mass
c             amass(i)=mass of ith outgoing particle
c   output:
c             pcm(1,i)=x-momentum of ith particle
c             pcm(2,i)=y-momentum of ith particle
c             pcm(3,i)=z-momentum of ith particle
c             pcm(4,i)=energy of ith particle
c             pcm(5,i)=momentum of ith particle
c             wtxlog=log of weight of event
c--------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      common /clatt/nlattc,npmax
      common/cffq/ffqlog(maxp)
      dimension ti(maxp),xi(maxp),si(maxp),zi(maxp)
      common/crnoz/rnoz(maxp-1)
      double precision ps(5)

      call utpri('hnbodz',ish,ishini,6)
      if(ish.ge.6)write(ifch,1200)np,tecm
      if(ish.ge.6)write(ifch,*)'particle masses:'
      if(ish.ge.6)write(ifch,'(1x,10f6.3)')(amass(n),n=1,np)

c initialization ktnbod=1
      ktnbod=ktnbod + 1
      if(ktnbod.gt.1) goto 1
c     !ffqlog(n) = log{ (4*pi)**n  / (n-1)! }
      ffqlog(1)=alog(4*pi)
      do n=2,maxp
      ffqlog(n)=ffqlog(n-1)+alog(4*pi/(n-1))
      enddo
    1 continue
c set wtxlog -infinity for np<2
      if(np.lt.2) goto 1001
c special treatment for np=2
      if(np.eq.2)then
      if(tecm.lt.amass(1)+amass(2)+0.00001)goto1000
      p0=utpcm(tecm,amass(1),amass(2))
      wtxlog=alog( 4*pi*p0
     */(1/sqrt(amass(1)**2+p0**2)+1/sqrt(amass(2)**2+p0**2)) )
      if(ish.ge.7)
     *write(ifch,*)'wtxlog:',wtxlog,'   (np=2 treatment)'
      bang=2*pi*rangen()
      cb=cos(bang)
      sb=sin(bang)
      c=2.0*rangen()-1.0
      s=sqrt(1.0-c*c)
      do 9 i=1,2
      is=2*i-3
      pcm(5,i)=p0
      pcm(1,i)=is*pcm(5,i)*s*cb
      pcm(2,i)=is*pcm(5,i)*s*sb
      pcm(3,i)=is*pcm(5,i)*c
      pcm(4,i)=sqrt(amass(i)**2+p0**2)
    9 continue
      goto1111
      endif
c stop if np too large
      if(np.gt.maxp) goto 1002
c initialization all ktnbod
      tm=0.0
      do 2 i=1,np
      tm=tm+amass(i)
    2 continue
      tt=tecm-tm
      if(tt.le.0.0) goto 1000
c prefactor
      wtxlog=alog(tt)*(np-1) + ffqlog(np)
      if(ish.ge.7)
     *write(ifch,*)'wtxlog:',wtxlog,'   (prefactor)'
c fill rnoz with np-1 random numbers
      if(keepr.eq.0)then
        do i= 1, np-1
          rnoz(i)=rangen()
        enddo
c update rnoz
      else
      !done in hnbrmz, look for "update rnoz"
      !  do lo=1,iomom
      !    j=1+rangen()*nlattc
      !    rnoz(j)=rangen()
      !  enddo
      endif
c calculate z_i distributed as i*z*(i-1)
      do i= 1, np-1
      zi(i)=rnoz(i)**(1./i)
      enddo
c calculate x_i
      xi(np)=1
      do i=np-1,1,-1
      xi(i)=zi(i)*xi(i+1)
      enddo
c calculate t_i, e_i, p_i
      if(ish.ge.9)write(ifch,*)'calculate t_i, e_i, p_i ...'
      si(1)=xi(1)*tt
      do i=2,np-1
      si(i)=xi(i)*tt
      enddo
      ti(1)=si(1)
      if(ti(1).le.0.)ti(1)=1e-10
      ti(np)=tt-si(np-1)
      if(ti(np).le.0.)ti(np)=1e-10
      do i=np-1,2,-1
      ti(i)=si(i)-si(i-1)
      if(ti(i).le.0.)ti(i)=1e-10
      enddo
      do i=1,np
      pcm(1,i)=0
      pcm(2,i)=0
      pcm(3,i)=0
      pcm(4,i)=ti(i)+amass(i)
      p52=ti(i)*(ti(i)+2*amass(i))
      if(p52.gt.0)then
      pcm(5,i)=sqrt(p52)
      else
      pcm(5,i)=ti(i)*sqrt(1+2*amass(i)/ti(i))
      endif
      enddo
c calculate wtxlog
      !~~~~~~~~~~~~~~~~~~~~~
      call hnbraw(7,20000,w)  !7,200   !ckw2018
      !~~~~~~~~~~~~~~~~~~~~~
      if(w.gt.0.)then
      wtxlog=wtxlog+alog(w)
      else
      wtxlog=wtxlog-1e+30
      endif
      do 7 i=1,np
      wtxlog=wtxlog+alog(pcm(5,i))+alog(ti(i)+amass(i))
    7 continue
      if(ish.ge.7)
     *write(ifch,*)'wtxlog:',wtxlog
c check
      eee=0
      do i=1,np
      eee=eee+pcm(4,i)
      enddo
      if(abs(eee-tecm)/tecm.gt.1e-4)print*,'hnbodz (1) energy pb'
     .,eee,tecm
c print
      if(ish.ge.7)then
      write(ifch,*)'momenta:'
      do j=1,4
      ps(j)=0
      enddo
      do i=1,np
      do j=1,4
      ps(j)=ps(j)+pcm(j,i)
      enddo
      write(ifch,'(1x,i3,5x,5f12.5)')i,(pcm(j,i),j=1,5)
      enddo
      ps(5)=dsqrt(ps(1)**2+ps(2)**2+ps(3)**2)
      write(ifch,'(1x,a4,8x,5f12.5)')'sum:',(sngl(ps(j)),j=1,5)
      endif
      if(w.le.0.)goto1111
c complete specification of event (random rotations and then deformations)
      call hnbrot
      if(ish.ge.7)write(ifch,*)'momenta after rotations:'
      call hnbrop(96,0)
      call hnbrod
      if(ish.ge.7)write(ifch,*)'momenta after deformations:'
      call hnbrop(96,1)
c check
      eee=0
      do i=1,np
      eee=eee+pcm(4,i)
      enddo
      if(abs(eee-tecm)/tecm.gt.1e-4)print*,'hnbodz (2) energy pb'
     .,eee,tecm

      goto1111

c error returns
 1000 continue
      if(ish.ge.6)
     *write(ifch,*)'available energy zero or negative -> wtxlog=-1e35'
      wtxlog=-1e35
      goto1111

 1001 continue
      if(ish.ge.6)
     *write(ifch,*)'less than 2 outgoing particles -> wtxlog=-1e35'
      wtxlog=-1e35
      goto1111

 1002 continue
      write(ifch,*)'too many outgoing particles'
      write(ifch,1150) ktnbod
 1150 format(47h0 above error detected in hnbody at call number,i7)
      write(ifch,1200) np,tecm
 1200 format(' np:',i6/' tecm:',f10.5)
      write(ifch,*)'particle masses:'
      write(ifch,'(1x,10f6.3)')(amass(jk),jk=1,np)
      stop

1111  continue
      call utprix('hnbodz',ish,ishini,6)
      return
      end

c-----------------------------------------------------------------------
      subroutine hnbolo(loops)
c-----------------------------------------------------------------------
c  loop over hnbody
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      parameter (mspecs=400)
      common/cspecs2/pspecs(0:mspecs+1),yie2volu,yiespecs
      a=0
      k=0
      do j=1,loops
c-c   if(mod(j,iterpr).eq.0)write(ifmt,*)'     iteration:',iter,j
      if(yiespecs.le.iocova)call hnbody
      if(yiespecs.gt.iocova)call hnbodz
      if(ish.ge.8)write(ifch,*)'j:',j,'   wtxlog:',wtxlog
           if(wtxlog.gt.-1e30)then
      k=k+1
      if(k.eq.1)c=wtxlog
           if(a.gt.0.)then
      if(alog(a).lt.wtxlog-c-20)then
      a=0
      c=wtxlog
      endif
           endif
      a=a+exp(wtxlog-c)
           endif
      if(ish.ge.8)write(ifch,*)'k:',k,'   c:',c
      enddo
      a=a/loops
      wtxlog=alog(a)+c
      return
      end

c-----------------------------------------------------------------------
      function hnbpdk(a,b,c)
c-----------------------------------------------------------------------
c  formerly pdk from cernlib
c  returns momentum p for twobody decay  a --> b + c
c           a, b, c are the three masses
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c  this p is related to twobody phase space as R2 = pi * p /a
c-----------------------------------------------------------------------
      double precision aa,bb,cc,a2,b2,c2
      aa=a
      bb=b
      cc=c
      a2=aa*aa
      b2=bb*bb
      c2=cc*cc
      if(a2 + (b2-c2)**2/a2-2.0*(b2+c2).le.0.)then
      hnbpdk = 0
      else
      hnbpdk = 0.5*dsqrt(a2 + (b2-c2)**2/a2 - 2.0*(b2+c2))
      endif
      return
      end

c----------------------------------------------------------------------
      subroutine hnbpad(k,n1,n2,n3,n4,jc)
c----------------------------------------------------------------------
c  k=1: determ pair indices n1,n2
c  k=2: determ pair indices n3,n4 (.ne. n1,n2)
c  k=1 and k=2 jc: flavour of pair
c----------------------------------------------------------------------
      include 'epos.inc'
      integer jc(nflav,2),jc1(nflav,2),jc2(nflav,2)
      common /clatt/nlattc,npmax
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog

      k1=0
      k2=0
      if(k.eq.2)then
      k1=n1
      k2=n2
      endif

c     determine n1,n2 and mm
c     ----------------------
    1 continue
      n1=1+rangen()*nlattc
      n1=min(n1,nlattc)
    2 continue
      n2=1+rangen()*nlattc
      n2=min(n2,nlattc)
      if(n2.eq.n1)goto2
      if(n2.lt.n1)then
      n1r=n1
      n1=n2
      n2=n1r
      endif
      if(k.eq.2)then
      if(n1.eq.k1.or.n1.eq.k2.or.n2.eq.k1.or.n2.eq.k2)goto1
      endif
      if(ish.ge.7)then
        write(ifch,*)'to be replaced 1 :',n1,ident(n1)
        write(ifch,*)'to be replaced 2 :',n2,ident(n2)
      endif

      call idtr7( ident(n1) ,jc1)
      call idtr7( ident(n2) ,jc2)
      call jcadd(1,jc1,jc2,jc) !add

      if(k.eq.2)then
      n3=n1
      n4=n2
      endif

      return
      end

c----------------------------------------------------------------------
      subroutine jcadd(isi,jc1,jc2,jc)
c----------------------------------------------------------------------
      include 'epos.inc'
      integer jc(nflav,2),jc1(nflav,2),jc2(nflav,2)
      do i=1,nflav
        do j=1,2
          jc(i,j)=jc1(i,j)+isi*jc2(i,j)
        enddo
      enddo
      do i=1,nflav
        j12=jc(i,1)-jc(i,2)
        if(j12.ge.0)then
          jc(i,1)=j12
          jc(i,2)=0
        else
          jc(i,1)=0
          jc(i,2)=-j12
        endif
      enddo
      end

c----------------------------------------------------------------------
      subroutine jcprint(ii,jc)
c----------------------------------------------------------------------
      include 'epos.inc'
      dimension jc(nflav,2)
      write(ifch,'(5x,i5,5x,3i5)')ii,
     .jc(1,1)-jc(1,2),jc(2,1)-jc(2,2),jc(3,1)-jc(3,2)
      end

c----------------------------------------------------------------------
      subroutine hnbpaj(iii,jc,id1,id2,wei)
c----------------------------------------------------------------------
c iii=1: For given flavour jc, returns hadron pair id1,id2 and weight
c iii=2: For given flavour jc, id1,id2, returns  weight
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(mspecs=400,mxids=381)
      parameter(mxpair=mspecs**2*4)
      parameter (mxidx=3**6)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cspec2/jspecs(2,nflav,mspecs+1)
      parameter (mxlk=375)
      common/cspec5/idxpair(-6:6,-6:6,-6:6),ipairst(2,mxidx)
      common/cspec6/wpairst(mxidx)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)
      dimension jc(nflav,2)!,jc2(nflav,2)

      wei=1
      if(iopair.eq.4)return

      wei=0

      iqu=jc(1,1)-jc(1,2)
      iqd=jc(2,1)-jc(2,2)
      iqs=jc(3,1)-jc(3,2)
      
      if(abs(iqu).gt.6)return
      if(abs(iqd).gt.6)return
      if(abs(iqs).gt.6)return

      idx=idxpair(iqu,iqd,iqs)
      !write(ifch,*)'idx',iqu,iqd,iqs,idx
      if(idx.eq.0)return

      ipair=ipairst(1,idx)
      !write(ifch,*)'ipair',ipair
      if(ipair.eq.0)return
 
      wtot=wpairst(idx) 

      if(iii.eq.1)then !~~~~~~~~~~~~~ random choice ~~~~~~~~~~~~~~~~

      r=rangen()
      wtt=0
      do i=1,ipair
        call idpairstget(1,i+ipairst(2,idx),id1)
        call idpairstget(2,i+ipairst(2,idx),id2)
        call wgtpairstget(i+ipairst(2,idx),wt)
        wei=wt/wtot
        wtt=wtt+wei
        !write(ifmt,'(a,i6,i6,2x,3i3,3x,2i7,2f9.4)')'hnbpaj',i
        !.  ,idx,iqu,iqd,iqs,id1,id2,wtt,r
        if(wtt.ge.r)goto 4
      enddo
      i=min(i,ipair)
  4   continue

      !write(ifmt,'(a,3i6,f9.4)')'hnbpaj',i,id1,id2,wei

      else !~~~~~~~~~~~~~~~ given ids ~~~~~~~~~~~~~~~~~~~~
 
      do i=1,ipair
        call idpairstget(1,i+ipairst(2,idx),id1xx)
        call idpairstget(2,i+ipairst(2,idx),id2xx)
        call wgtpairstget(i+ipairst(2,idx),wt)
        wei=wt/wtot
        if(id1.eq.id1xx.and.id2.eq.id2xx)goto 5
      enddo
  5   continue
 
      endif

      end

c----------------------------------------------------------------------
      subroutine hnbpajini
c----------------------------------------------------------------------
c  initialize array to speed up hnbpaj calculation
c  store sum of weights iwpair of possible pairs in an array
c  for any combinations of quarks
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(mspecs=400,mxids=381)
      parameter(mxpair=mspecs**2*4)
      parameter(mxpairst=145000)
      parameter (mxidx=3**6)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cspec2/jspecs(2,nflav,mspecs+1)
      common/cspecs2/pspecs(0:mspecs+1),yie2volu,yiespecs
      parameter (mxlk=375)
      common/cspec5/idxpair(-6:6,-6:6,-6:6),ipairst(2,mxidx)
      common/cspec6/wpairst(mxidx)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)
      dimension ids(mxids),wgts(mxids)
      if(iopair.eq.4)return

c      write(ifmt,*)' Initialize droplet decay ...'

c     construct possible pairs id1,id2
c     --------------------------------
      idx=0
      jpair=0
      do i=1,mxpairst
        call idpairstset(1,i,0)  !id1 of pair
        call idpairstset(2,i,0)  !id2 of pair
        call wgtpairstset(i,0.) !weight pair
      enddo

      do iqu=-6,6
      do iqd=-6,6
      do iqs=-6,6

      !iqu,iqd,iqs reprents needed flavor to compensate removed one

      itot=iqu+iqd+iqs
      iatot=abs(iqu)+abs(iqd)+abs(iqs)
      if(mod(itot,3).ne.0.or.iatot.gt.12)then
         idxpair(iqu,iqd,iqs)=0
         !write(ifmt,'(20x,3i4,a)')iqu,iqd,iqs,' skip'
         goto4
      endif

c initialize

      idx=idx+1
      idxpair(iqu,iqd,iqs)=idx

      ipairst(2,idx)=jpair
      ipair=0
      wsum=0
      wpairst(idx)=0!weight sum
      do i=1,mxids
        ids(i)=0
        wgts(i)=0
      enddo
      !write(ifmt,'(a,i4,2x,3i3)')'hnbpajini',idx,iqu,iqd,iqs

c --------

      do i1=1,nspecs+1 !first pair element 

      id1=ispecs(i1)
      prob1 = pspecs(i1)-pspecs(i1-1)

      iquu=iqu-ifok(1,i1)
      iqdd=iqd-ifok(2,i1)
      iqss=iqs-ifok(3,i1)
      iqcc=0

      nids=0

      if(abs(iquu).gt.3)goto3
      if(abs(iqdd).gt.3)goto3
      if(abs(iqss).gt.3)goto3

      call lkfokget(1,4+iquu,4+iqdd,4+iqss,4+iqcc, lkfok1 )
      !write(ifmt,'(i7,3x,3i4,f8.3,i7)')id1,iquu,iqdd,iqss,prob1,lkfok1
      if(lkfok1.gt.0)then
        if(lkfok1.gt.mxlk)stop'HNBPAJINI: mxlk too small'
        do ii=1,lkfok1
          nids=nids+1
          call lkfokget(1+ii,4+iquu,4+iqdd,4+iqss,4+iqcc, i2 )
          ids(nids)=ispecs(i2)
          prob2 = pspecs(i2)-pspecs(i2-1)
          wgts(nids) = prob1 * prob2
          !print*,ii,nids,i2,ids(nids),prob2
        enddo
      endif
      if(nids.eq.0)goto3

      if(nids.gt.mxpair)call utstop('hnbpajini: mxpair too small&')

      call idtr8(id1,ju1,jd1,js1)
      do k=1,nids
        ipair=ipair+1
        jpair=jpair+1
        wsum=wsum+wgts(k)
        if(ipair+ipairst(2,idx).ne.jpair)stop'#### ERROR 10072018 ####'
        call idpairstset(1,ipair+ipairst(2,idx),id1)
        call idpairstset(2,ipair+ipairst(2,idx),ids(k))
        if(ipair+ipairst(2,idx).gt.mxpairst)
     .  stop'#### ERROR 19022019 ####'
        call wgtpairstset( ipair+ipairst(2,idx),wgts(k) )
        id2=ids(k)
        call idtr8(id2,ju2,jd2,js2)
        !write(ifmt,'(i4,2x,3i4,3x,3i4,2x,3i4,2x,2i9,2f9.4)')ipair,
        !.  iqu,iqd,iqs, ju1,jd1,js1, ju2,jd2,js2, id1,id2, wgts(k)*1000
        !.  ,wsum*1000
        if(iqu.ne.ju1+ju2)stop'ERROR 27062018a'
        if(iqd.ne.jd1+jd2)stop'ERROR 27062018b'
        if(iqs.ne.js1+js2)stop'ERROR 27062018b'
      enddo

    3 continue

      enddo !i1

      ipairst(1,idx)=ipair
      wpairst(idx) = wsum  
      !print*,'hnbpajini ipairst',idx,ipairst(1,idx),ipairst(2,idx)+1
      !.,ipair+ipairst(2,idx)

      wtt=0
      do i=1,ipair
        call wgtpairstget(i+ipairst(2,idx),wt)
        wtt=wtt+wt
      enddo
      if(abs(wtt-wsum).gt.1e-4*wsum)stop'####### ERROR 29062018 #######'
      !print*,wtt,wsum

      !write(ifmt,'(a,i4,2x,3i3,i8,f9.3)')'hnbpajini',idx,iqu,iqd,iqs
      !.,ipair,wsum*1000

c     no pair found
c     -------------
      if(ipair.eq.0)then
        if(wsum.gt.1e-5)call utstop('hnbpajini: iwpair.ne.0&')
      endif

  4   continue

      enddo
      enddo
      enddo
      
      return
      end

c--------------------------------------------------------------------
      subroutine hnbraw(npx,npy,w)
c--------------------------------------------------------------------
c returns random walk fctn w=w(0,p_1,p_2,...,p_n) for noncovariant
c phase space integral (see hagedorn, suppl nuov cim ix(x) (1958)646)
c input: dimension np and momenta p_i=pcm(5,i) via /confg/
c    1   < np <= npx : hagedorn method
c    npx < np <= npy : integral method
c    npy < np        : asymptotic method
c--------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      integer ii(maxp),isi(maxp)
      double precision ppcm(maxp),ww,ppsum,ppmax
      common /arh6/xh6(6),wh6(6)
      external hnbrax
      common/cepsr/nepsr
      if(ish.ge.9)write(ifch,*)('-',i=1,10)
     *,' entry sr hnbraw ',('-',i=1,30)

      if(np.lt.3)call utstop('hnbraw: np must be at least 3&')

      kper=5
      pi=3.1415927
      pmax=0
      do i=1,np
      pmax=pmax+pcm(5,i)
      enddo
      wio=0
      win=0
      whd=0

      npxx=0
      pp2=0
      do i=1,np
        pp2=pp2+pcm(5,i)**2
        if(pcm(5,i).gt.0.1)npxx=npxx+1
      enddo
      pp26=pp2/6
      rpp26=sqrt(pp26)
      !the gaussian approx of the integrand of the RW function is
      !      f_G = 1/(2*pi**2) * x**2 * exp( -(P*x)**2 )
      ! with with P = rpp26
      xmx=1/rpp26 !location of the maximum of f_G 
      xmx3=xmx*3 !estimate of the upper bound  of trapezoidal integration

c     sum p_i - 2*p_max not positive
c     ------------------------------
      px=0
      ps=0
      do i=1,np
      px=max(px,pcm(5,i))
      ps=ps+pcm(5,i)
      enddo
      if(ps-2*px.le.0.)then
      w=0
      if(ish.ge.7)write(ifch,'(1x,a,e12.5,4x)')
     *'sum p_i - 2*p_max not positive -->  w:',w
      goto1000
      endif

      if(np.le.npx)goto 9
      if(np.le.npy)goto 22

c     asymptotic method based on f_G (see above)
c     ------------------------------
ckw  asymptotic formula not precise enough even for np ~ 200 - 500
ckw  BUT the approximate integrand f_G can be employed for numerical
ckw  Gauss Hermite integration (very efficient), see next paragraph

      was=0
      do i=1,np
      was=was+pcm(5,i)**2
      enddo
      was=(was*2*pi/3)**(-1.5)
      if(ish.ge.7)write(ifch,'(1x,a,e12.5,4x)')
     *'asymptotic method: was:',was

      if(np.gt.npy)then
      w=was
      goto1000
      endif

c     integral method
c     ---------------
  22  continue
      if(ish.ge.9)write(ifch,*)'integral method...'
      !~~~~trapeze method
      eps=epsr
      b=xmx3*2 !to be on the safe side
      !b=1.5*np/pmax  !pi*np*kper/pmax
      !call uttrap2(hnbrax,0.,b,eps,win,iok,np,it)
      !w=win
      !~~~~Gauss-Hermite method
      !Gauss Hermite based on f ~ f_G ~ exp( -(P*x)**2 ),  P = rpp26
      !beyond np=30 very precise (>1/1000) 
      !beyond np=10 still very precise (>1/1000) in almost all cases
      ww=0
      do k=1,6
        ww=ww+hnbray(rpp26,xh6(k)/rpp26)/rpp26*wh6(k)
      enddo
      iok=1 !Gauss-Hermite always ok by construction
      !~~~~~Compare both methods 
      !if(abs(w-ww)/ww.gt.1e-3)
      !.write(*,'(a,2i5,f8.3,3x,3f8.3)')'trap, gauss hermite:',np,npxx
      !.,rpp26, w*1000, ww*1000,abs(w-ww)/ww
      w=ww

      if(iok.eq.1)goto1000

      !if(ish.ge.1)then
        write(ifch,*)'WARNING hnbraw - iok not 1 (should be 1)'
      !endif
      goto1000

c     hagedorn method (double)
c     ------------------------
    9 continue
      ppmax=0
      do i=1,np
      ppcm(i)=pcm(5,i)
      ppmax=ppmax+ppcm(i)
      enddo
      ww=0
      do i=1,np
      ii(i)=0
      isi(i)=1
      enddo
      ppsum=ppmax
      i=0
      iprosi=1
      ww=iprosi*(ppsum/ppmax)**(np-3)
      if(ish.ge.8)
     *write(ifch,'(4x,i5,12x,f7.2,i5,f11.2)')np,sngl(ppsum)
     *,iprosi,sngl(ww)
    5 continue
      i=i+1
      if(i.gt.np)goto6
      if(ii(i).eq.1)goto5
      iprosi=-iprosi
      isi(i)=-isi(i)
      ppsum=ppsum+2*isi(i)*ppcm(i)
           if(ppsum.gt.0.or.ppsum.eq.0..and.isi(i).gt.0)then
      ww=ww+iprosi*(ppsum/ppmax)**(np-3)
      if(ish.ge.8)
     *write(ifch,'(4x,2i5,2f7.2,i5,f11.2)')
     *np,i,sngl(2*isi(i)*ppcm(i)),sngl(ppsum),iprosi,sngl(ww)
           else
      if(ish.ge.8)
     *write(ifch,'(4x,2i5,2f7.2,i5,4x,a)')
     *np,i,sngl(2*isi(i)*ppcm(i)),sngl(ppsum),iprosi,'not counted'
           endif
      ii(i)=1
      if(i.gt.1)then
      do j=1,i-1
      ii(j)=0
      enddo
      endif
      i=0
      goto5
    6 continue
      do i=1,np
      ww=ww*pmax/ppcm(i)/2./i
      enddo
      ww=-ww/pmax**3/pi/2.*np*(np-1)*(np-2)
      whd=ww
      if(ish.ge.7)write(ifch,'(1x,a,e12.5,4x,a)')
     *'hagedorn method:   whd:',whd,'double precision'

      w=whd

1000  continue
      if(ish.ge.9)write(ifch,*)('-',i=1,30)
     *,' exit sr hnbraw ',('-',i=1,10)
      return
      end

c--------------------------------------------------------------------
      function hnbrax(x)
c--------------------------------------------------------------------
c returns integrand for random walk fctn w=w(0,p_1,p_2,...,p_n):
c 1./(2*pi**2) * x**2 * prod[sin(p_i*x)/(p_i*x)]
c input: dimension np and momenta p_i=pcm(5,i) via /confg/
c--------------------------------------------------------------------
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      common/cnsta/pi,pii,hquer,prom,piom,ainfin
      hnbrax= pii * x**2
      do i=1,np
      px=pcm(5,i)*x
      if(px.ne.0.)hnbrax=hnbrax*sin(px)/px
      enddo
      return
      end

c--------------------------------------------------------------------
      function hnbray(p,x)
c--------------------------------------------------------------------
      hnbray = hnbrax(x) * exp( (p*x)**2 )
      return
      end

c-----------------------------------------------------------------
      subroutine uttrap2(func,a,b,eps,s,iok,np,it)
c-----------------------------------------------------------------
c trapezoidal method for integration.
c input: fctn func and limits a,b
c output: value s of the integral
c-----------------------------------------------------------------
      include 'epos.inc'

      INTEGER JMAX
      REAL a,b,func,s
      EXTERNAL func
      PARAMETER (JMAX=20)
CU    USES uttras
      INTEGER j
      REAL olds

c      if(np.lt.20)then
c        fac=100
c      elseif(np.lt.50)then
c        fac=10
c      else
        fac=1
c      endif 

      itmax=8
      it=0
      b=b*3./5.
    3 continue
      it=it+1
      b=b*5./3.

      olds=-1.e30
      do j=1,JMAX
        call uttras2(func,a,b,j,eps,s1,s2,iok)
        s=s1+s2
        ds=abs(s-olds)
        !if(it.ge.2)write(ifch,*)it,np,j,b,s,ds/abs(olds),iok
        if (ds.lt.fac*eps*abs(olds)) goto 11
        olds=s
      enddo
  11  continue

      if(iok.eq.0.and.it.lt.itmax)goto3

      if(iok.eq.1)then
        write(ifch,*)'WARNING uttrap2 - accuracy np int ok it:'
     *  ,ds/abs(olds),np,s,iok,it
        !do i=1, 2**(JMAX-1)
        !  x=a+i*(b-a)/2**(JMAX-1)
        !  write(ifch,*)x,func(x)
        !enddo
      endif
      END

c-----------------------------------------------------------------
      subroutine uttras2(func,a,b,n,eps,s1,s2,iok)
c-----------------------------------------------------------------
c performs one iteration of the trapezoidal method for integration
c-----------------------------------------------------------------
      EXTERNAL func
      if(n.eq.1)then
        s1=0.5*(b-a)*func(a)
        s2=0.5*(b-a)*func(b)
      elseif(n.eq.2)then
        del=(b-a)
        x=a+0.5*(b-a)
        w=func(x)
        s1=0.5*s1 + del*w 
        s2=0.5*s2 + del*w 
      else
        it=2**(n-2)
        del=(b-a)/it
        c=(a+b)/2 
        x1=a+0.5*del
        x2=c+0.5*del
        sum1=0.
        sum2=0.
        do j=1,it/2
          sum1=sum1+func(x1)
          sum2=sum2+func(x2)
          x1=x1+del
          x2=x2+del
        enddo
        s1=0.5*( s1 + del*sum1 )
        s2=0.5*( s2 + del*sum2 )
      endif
      s=s1+s2
      iok=0
      if(abs(s-s1).lt.eps*abs(s))iok=1
      end

c----------------------------------------------------------------------
      subroutine hnbrmz(ielast)
c----------------------------------------------------------------------
c  removes intermediate zeros from ident
c  updates np
c  update rnoz and rnoy
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
c      integer identx(maxp)
      common /clatt/nlattc,npmax
      common/crnoz/rnoz(maxp-1)
      common/crnoy/rnoy(3,maxp-1)
      if(ish.ge.9)write(ifch,*)('-',i=1,10)
     *,' entry sr hnbrmz ',('-',i=1,30)
      if(np.eq.0)goto1000
      if(ielast.gt.2)stop'ERROR 03062018'

      if(ielast.eq.2)then !elastic
        m=min( np , int(1+rangen()*np) )
        rnoz(m)=rangen()
        do k=1,3
         rnoy(k,m)=rangen()
        enddo        
        return
      endif

c      do i=1,np
c      identx(i)=ident(i)
c      enddo
c      npx=np

      i=0
      np=nlattc+1
      npp=0

      !----------------
      !find next zero i
      !----------------
    1 i=i+1
      if(i.gt.nlattc)then
        np=nlattc
        goto1000
      endif
      if(ident(i).ne.0)then
        npp=i
        goto1
      endif

c    2 np=np-1
c      if(np.eq.0)goto1000
c      if(ident(np).eq.0)goto2
c      if(i.eq.np+1)goto1000
c      ident(i)=ident(np)
c      ident(np)=0
c      !update rnoz 
c      if(np.lt.nlattc)then
c        rnoz(i)=rnoz(np)
c        rnoz(np)=rangen()
c      else  
c        rnoz(i)=rangen()
c      endif

      !--------------------------------------------------
      ! find next nonzero m -- move m to i -- make m zero 
      !--------------------------------------------------
      m=i
    3 m=m+1
      if(m.gt.nlattc)then
        np=npp
        goto1000
      endif
      if(ident(m).eq.0)goto3
      npp=i
      ident(i)=ident(m)
      ident(m)=0
      !update rnoz 
      if(m.lt.nlattc)then
        rnoz(i)=rnoz(m)
        rnoz(m)=rangen()
        do k=1,3
         rnoy(k,i)=rnoy(k,m)
         rnoy(k,m)=rangen()
        enddo
      else  
        rnoz(i)=rangen()
        do k=1,3
         rnoy(k,i)=rangen()
        enddo
      endif

      goto1

1000  continue
      if(ish.ge.9)write(ifch,*)('-',i=1,30)
     *,' exit sr hnbrmz ',('-',i=1,10)
      end

c----------------------------------------------------------------------
      subroutine hnbrod
c----------------------------------------------------------------------
c deformes polygon of a sequence of arbitrarily rotated momentum
c vectors such that the polygon gets closed
c    input: pcm(1-3,i) representing polygon
c    output: pcm(1-3,i) representing closed polygon
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      real x(3),y(3),z(3),w(3)
      if(ish.ge.8)write(ifch,*)'sr hnbrod: polygon deformation:'

      err=0.01

      kmax=1000
      fac=0.30
      x2max=(err*tecm)**2

      if(ish.ge.8)write(ifch,'(a,i4,a,f12.6)')
     *' kmax:',kmax,'   x2max:',x2max

      x(1)=0
      x(2)=0
      x(3)=0
      do i=1,np
      x(1)=x(1)+pcm(1,i)
      x(2)=x(2)+pcm(2,i)
      x(3)=x(3)+pcm(3,i)
      enddo ! i

      k=0
   1  continue

      x2=x(1)**2+x(2)**2+x(3)**2
      if(ish.ge.8)write(ifch,'(a,i3,a,3f9.3,a,f12.6)')
     *' it',k,':   x:',x,'      x2:',x2
      if(x2.le.x2max)goto1000
      if(k.gt.kmax)goto1001

      k=k+1
      ir=1+rangen()*np
      ir=min(ir,np)

      z(1)=-x(1)
      z(2)=-x(2)
      z(3)=-x(3)
      x(1)=x(1)-pcm(1,ir)
      x(2)=x(2)-pcm(2,ir)
      x(3)=x(3)-pcm(3,ir)
      y(1)=pcm(1,ir)
      y(2)=pcm(2,ir)
      y(3)=pcm(3,ir)
      if(ish.ge.9)write(ifch,'(a,i3,a,3f9.3,a,3f9.3,a,i4)')
     *' it',k,':   x:',x,'   y:',y,'  ir:',ir
      xxx=x(1)**2+x(2)**2+x(3)**2
      yyy=y(1)**2+y(2)**2+y(3)**2
      zzz=z(1)**2+z(2)**2+z(3)**2
         if(xxx.gt.0..and.yyy.gt.0..and.zzz.gt.0.)then
c      xx=sqrt(xxx)
      yy=sqrt(yyy)
      zz=sqrt(zzz)
      a=min(fac,fac*yy/zz)
      w(1)=y(1)+a*z(1)
      w(2)=y(2)+a*z(2)
      w(3)=y(3)+a*z(3)
      www=w(1)**2+w(2)**2+w(3)**2
         if(www.gt.0.)then
      ww=sqrt(www)
      y(1)=yy/ww*w(1)
      y(2)=yy/ww*w(2)
      y(3)=yy/ww*w(3)
      pcm(1,ir)=y(1)
      pcm(2,ir)=y(2)
      pcm(3,ir)=y(3)
         endif
         endif
      x(1)=x(1)+y(1)
      x(2)=x(2)+y(2)
      x(3)=x(3)+y(3)
      if(ish.ge.9)write(ifch,'(a,i3,a,3f9.3,a,3f9.3,a,i4)')
     *' it',k,':   x:',x,'   y:',y,'  ir:',ir

      goto1

 1001 continue
      call utmsg('hnbrod&')
      write(ifch,*)'*****  total 3-momentum nonzero'
      write(ifch,'(3f12.5,5x,2f12.5)')(x(j),j=1,3),x2,x2max
      call utmsgf

 1000 continue
      return

      end

c----------------------------------------------------------------------
      subroutine hnbrop(ishx,ichk)
c----------------------------------------------------------------------
c  prints momenta of configuration (essentially to check rotation procedure)
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      double precision ps(5)
      err=0.01
      do j=1,4
      ps(j)=0
      enddo
      do i=1,np
      do j=1,4
      ps(j)=ps(j)+pcm(j,i)
      enddo
      if(ish.ge.ishx)write(ifch,'(1x,i3,5x,5f12.5)')i,(pcm(j,i),j=1,3)
     *,sqrt(pcm(1,i)**2+pcm(2,i)**2+pcm(3,i)**2),pcm(5,i)
      enddo
      ps(5)=dsqrt(ps(1)**2+ps(2)**2+ps(3)**2)
      if(ish.ge.ishx)write(ifch,'(1x,a4,8x,5f12.5)')
     *'sum:',(sngl(ps(j)),j=1,5)
           if(ichk.eq.1)then
           if(dabs(ps(1)).gt.err*tecm.or.dabs(ps(2)).gt.err*tecm
     *.or.dabs(ps(3)).gt.err*tecm)then
      call utmsg('hnbrop&')
      write(ifch,*)'*****  total 3-momentum nonzero'
      write(ifch,'(9x,5f12.5)')(sngl(ps(j)),j=1,5)
      call utmsgf
           endif
           endif
      return
      end

c----------------------------------------------------------------------
      subroutine hnbrot
c----------------------------------------------------------------------
c rotates momenta of /confg/ randomly
c   input: pcm(5,i)
c   output: pcm(1-3,i)
c----------------------------------------------------------------------
      common/cnsta/pi,pii,hquer,prom,piom,ainfin
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      real u(3)

      do i=1,np
      u(3)=2.*rangen()-1.
      phi=2.*pi*rangen()
      u(1)=sqrt(1.-u(3)**2)*cos(phi)
      u(2)=sqrt(1.-u(3)**2)*sin(phi)
      pcm(1,i)=pcm(5,i)*u(1)
      pcm(2,i)=pcm(5,i)*u(2)
      pcm(3,i)=pcm(5,i)*u(3)
      enddo

      return
      end

cc-------------------------------------------------------------------
c      subroutine hnbrt2old(c,s,c2,s2,pr,i)
cc-------------------------------------------------------------------
cc  formerly subr rotes2 from cernlib
cc  this subr now does two rotations (xy and xz)
cc-------------------------------------------------------------------
c      parameter(maxp=6000)
c      dimension pr(5*maxp)
c      k1 = 5*i - 4
c      k2 = k1 + 1
c      sa = pr(k1)
c      sb = pr(k2)
c      a      = sa*c - sb*s
c      pr(k2) = sa*s + sb*c
c      k2 = k2 + 1
c      b = pr(k2)
c      pr(k1) = a*c2 - b*s2
c      pr(k2) = a*s2 + b*c2
c      return
c      end
c
c-------------------------------------------------------------------
      subroutine hnbrt2(c,s,c2,s2,pr,i)
c-------------------------------------------------------------------
c  formerly subr rotes2 from cernlib
c  this subr now does two rotations (xy and xz)
c-------------------------------------------------------------------
      parameter(maxp=6000)
      dimension pr(5,maxp)
      k1 = 5*i - 4
      k2 = k1 + 1
      sa = pr(1,i)
      sb = pr(2,i)
      a      = sa*c - sb*s
      pr(2,i) = sa*s + sb*c
      k2 = k2 + 1
      b = pr(3,i)
      pr(1,i) = a*c2 - b*s2
      pr(3,i) = a*s2 + b*c2
      return
      end

cc-----------------------------------------------------------------------
c      subroutine hnbsor(a,n)
cc-----------------------------------------------------------------------
cc cern proglib# m103    flpsor          .version kernfor  3.15  820113
cc orig. 29/04/78
cc-----------------------------------------------------------------------
cc   sort the one-dimensional floating point array a(1),...,a(n) by
cc   increasing values
cc-----------------------------------------------------------------------
c      dimension a(*)
c      common /slate/ lt(20),rt(20)
c      integer r,rt
cc
c      level=1
c      lt(1)=1
c      rt(1)=n
c   10 l=lt(level)
c      r=rt(level)
c      level=level-1
c   20 if(r.gt.l) go to 200
c      if(level) 50,50,10
cc
cc   subdivide the interval l,r
cc     l : lower limit of the interval (input)
cc     r : upper limit of the interval (input)
cc     j : upper limit of lower sub-interval (output)
cc     i : lower limit of upper sub-interval (output)
cc
c  200 i=l
c      j=r
c      m=(l+r)/2
c      x=a(m)
c  220 if(a(i).ge.x) go to 230
c      i=i+1
c      go to 220
c  230 if(a(j).le.x) go to 231
c      j=j-1
c      go to 230
cc
c  231 if(i.gt.j) go to 232
c      w=a(i)
c      a(i)=a(j)
c      a(j)=w
c      i=i+1
c      j=j-1
c      if(i.le.j) go to 220
cc
c  232 level=level+1
c      if(level.gt.20)stop'level too large'
c      if((r-i).ge.(j-l)) go to 30
c      lt(level)=l
c      rt(level)=j
c      l=i
c      go to 20
c   30 lt(level)=i
c      rt(level)=r
c      r=j
c      go to 20
c   50 return
c      end
c
c-----------------------------------------------------------------------
      subroutine hnbspd(iopt)
c-----------------------------------------------------------------------
c  defines particle species and masses and degeneracies.
c  input:
c    iopt=odd number: massless
c    iopt=even number: same as iopt-1, but massive
c    iopt= 1: pi0 (massless)
c    iopt= 2: pi0
c    iopt= 3: pi-,pi0,pi+ (massless)
c    iopt= 4: pi-,pi0,pi+
c    iopt= 5: pi-,pi0,pi+,prt,aprt,ntr,antr (massless)
c    iopt= 6: pi-,pi0,pi+,prt,aprt,ntr,antr
c    iopt= 7: 25 hadrons (massless)
c    iopt= 8: 25 hadrons
c    iopt= 9: 54 hadrons (massless)
c    iopt=10: 54 hadrons
c    iopt=11:  3 quarks  (massless)
c    iopt=12:  3 quarks
c    iopt=13:  54 hadrons + J/psi   (massless)
c    iopt=14:  54 hadrons + J/psi
c    iopt=15:  54 hadrons + J/psi + H  (massless)
c    iopt=16:  54 hadrons + J/psi + H
c    iopt=18:  54 hadrons (mass,g equal case 10, flavorwise pi0)
c    iopt=19:  54 hadrons (mass=0.5,g=1, flavorwise equal case 9)
c    iopt=21:  54 hadrons from ptl table (massless)
c    iopt=22:  54 hadrons from ptl table
c    iopt=23: 371 hadrons from ptl table (massless)
c    iopt=24: 371 hadrons from ptl table 
c    iopt=25: 371 hadrons from ptl table (flavorless)
c  output:
c    nspecs: nr of species
c    ispecs: id's
c    aspecs: masses
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      double precision das
      character dname*30
      parameter(nspe01=1,nspe03=3,nspe05=7,nspe07=25,nspe09=54)
      parameter(nspe11=6,nspe13=55,nspe15=56)
      real jspe01(nspe01),jspe03(nspe03),jspe05(nspe05),jspe07(nspe07)
     *,jspe09(nspe09),jspe11(nspe11),jspe13(nspe13),jspe15(nspe15)
      data jspe01/   110 /
      data jspe03/   110,  120, -120 /
      data jspe05/   110,  120, -120, 1120,-1120, 1220,-1220 /
      data jspe07/
     *   110,  120, -120,  130, -130,  230, -230,  220,  330
     *, 1120,-1120, 1220,-1220, 1130,-1130, 2130,-2130
     *, 1230,-1230, 2230,-2230, 1330,-1330, 2330,-2330 /
      data jspe09/
     *   110,  120, -120,  130, -130,  230, -230,  220,  330
     *,  111,  121, -121,  131, -131,  231, -231,  221,  331
     *, 1120,-1120, 1220,-1220, 1130,-1130, 2130,-2130
     *, 1230,-1230, 2230,-2230, 1330,-1330, 2330,-2330
     *, 1111,-1111, 1121,-1121, 1221,-1221, 2221,-2221, 1131,-1131
     *, 1231,-1231, 2231,-2231, 1331,-1331, 2331,-2331, 3331,-3331 /
      data jspe11/
     *     1,   -1,    2,   -2,    3,   -3   /
      data jspe13/
     *   110,  120, -120,  130, -130,  230, -230,  220,  330
     *,  111,  121, -121,  131, -131,  231, -231,  221,  331
     *, 1120,-1120, 1220,-1220, 1130,-1130, 2130,-2130
     *, 1230,-1230, 2230,-2230, 1330,-1330, 2330,-2330
     *, 1111,-1111, 1121,-1121, 1221,-1221, 2221,-2221, 1131,-1131
     *, 1231,-1231, 2231,-2231, 1331,-1331, 2331,-2331, 3331,-3331
     *, 441 /
      data jspe15/
     *   110,  120, -120,  130, -130,  230, -230,  220,  330
     *,  111,  121, -121,  131, -131,  231, -231,  221,  331
     *, 1120,-1120, 1220,-1220, 1130,-1130, 2130,-2130
     *, 1230,-1230, 2230,-2230, 1330,-1330, 2330,-2330
     *, 1111,-1111, 1121,-1121, 1221,-1221, 2221,-2221, 1131,-1131
     *, 1231,-1231, 2231,-2231, 1331,-1331, 2331,-2331, 3331,-3331
     *, 441 , 30 /
      data ncnthnb/0/
      save ncnthnb
      ncnthnb=ncnthnb+1
      gg=0.

      if(iopt.le.19)then

        ioptx=(1+iopt)/2*2-1
        if(ioptx.eq.1)nspecs=nspe01
        if(ioptx.eq.3)nspecs=nspe03
        if(ioptx.eq.5)nspecs=nspe05
        if(ioptx.eq.7)nspecs=nspe07
        if(ioptx.eq.9)nspecs=nspe09
        if(ioptx.eq.11)nspecs=nspe11
        if(ioptx.eq.13)nspecs=nspe13
        if(ioptx.eq.15)nspecs=nspe15
        if(iopt.eq.18)nspecs=nspe09
        if(iopt.eq.19)nspecs=nspe09
        do i=1,nspecs
          if(ioptx.eq.1)ispecs(i)=jspe01(i)
          if(ioptx.eq.3)ispecs(i)=jspe03(i)
          if(ioptx.eq.5)ispecs(i)=jspe05(i)
          if(ioptx.eq.7)ispecs(i)=jspe07(i)
          if(ioptx.eq.9)ispecs(i)=jspe09(i)
          if(ioptx.eq.11)ispecs(i)=jspe11(i)
          if(ioptx.eq.13)ispecs(i)=jspe13(i)
          if(ioptx.eq.15)ispecs(i)=jspe15(i)
          if(iopt.eq.17)ispecs(i)=jspe09(i)
          if(iopt.eq.18)ispecs(i)=jspe09(i)
          if(iopt.eq.19)then
            aspecs(i)=0.5
          elseif(ioptx.eq.iopt)then !odd iopt
            aspecs(i)=0
          else
            id=ispecs(i)
            call idmass(id,am)
            aspecs(i)=am
          endif
          call hnbspi(ispecs(i),gg)
          if(iopt.eq.19)gg=1
          gspecs(i)=gg
          if(iopt.eq.18)then
            id=ispecs(i)
            if(id.gt.0)ispecs(i)=990000+id
            if(id.lt.0)ispecs(i)=-990000+id
          endif
          !print*,i,ispecs(i),aspecs(i),gspecs(i)
        enddo
        ispecs(nspecs+1)=0
        aspecs(nspecs+1)=0
        gspecs(nspecs+1)=0

      elseif(iopt.le.20)then

        !not used
        stop'####### ERROR in hnbspd (2) #######'

      elseif(iopt.le.25)then

        if(iopt.eq.21.or.iopt.eq.22)iflagmx=1
        if(iopt.eq.23.or.iopt.eq.24)iflagmx=2
        if(iopt.eq.25)iflagmx=2
        m=0
        last=0
        do
          call getNextParticle(dname,das,iQ,idxx,igs,js,iB,iS
     .         ,last,iflag)
          !--------------------------------------------------
          !      iflag=1 : 54 "basic" hadrons
          !      iflag=2 : charm=bottom=top=0 
          !      iflag=3 : all
          !-------------------------------------------------- 
          ! old method : call igetparticle(dname
          !.,das,iQ,idPDG,igs,js,iB,iS,last,iflag) 
          !--------------------------------------- 
          if(abs(last).eq.1)then
           if(iflag.gt.0.and.iflag.le.iflagmx)then
            m=m+1
            aa=das
            gg=igs
ctp2406 introduce some corrections for fine tuning (tuned to LHC) ???????????????
            fac=1
            call idflav(idxx,ifl1,ifl2,ifl3,jspin,index)
            ifls=0
            iflb=0
            if(abs(ifl1).eq.3)ifls=ifls+1
            if(abs(ifl2).eq.3)ifls=ifls+1
            if(abs(ifl3).eq.3)ifls=ifls+1
            if(ifls.eq.1)then
              if(abs(idxx).lt.1000)then
                fac=facts       !kaons
c                if(abs(idxx)/10.eq.23.or.abs(idxx).eq.20)fac=fac*0.85 !help in some places but worth for underlying events with strange particles and not really justified (the difference between produced and observed may come from hadronic rescattering "hidding" part of the particles ?
              else
                fac=factb    !lambda,sigma      !in underlying event with strange particles, it seems to be NOT necessary ? (but then MB yield is too low ?)
              endif
c            elseif(ifls.eq.0.and.abs(idxx).gt.1000)then
c              fac=1./factb         !reduce Omega
            endif
            gg=gg*fac
ctp2406 end of correction ??????????????
            ! old method ! idxx=idtrafo('pdg','nxs',idPDG)
            if(m.gt.mspecs)then
              print*,'####### m > mspecs :',m,mspecs 
              stop'####### ERROR in hnbspd #######'
            endif
            if(iopt.eq.25)then
              id=idxx
              if(id.gt.0)idxx=990000+id
              if(id.lt.0)idxx=-990000+id
            endif
            ispecs(m)=idxx
            aspecs(m)=aa
            gspecs(m)=gg 
c            print*,'particle',m,idxx,aa,gg,iflag 
           endif
          endif
          if(last.lt.0)goto1
        enddo
   1    continue
        nspecs=m
        ispecs(nspecs+1)=0
        aspecs(nspecs+1)=0
        gspecs(nspecs+1)=0

      else
        call utstop('ERROR hnbspd: invalid option&')
      endif 

      if(ncnthnb.eq.1)call hnbgra

      return
      end

c-----------------------------------------------------------------------
      subroutine hnbspd2
c-----------------------------------------------------------------------
      include 'epos.inc'
c-in
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
c-out
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)
      common/cspec2/jspecs(2,nflav,mspecs+1)
      parameter (mxlk=375)
c-local
      integer jc(nflav,2)

      do nf=1,nflav
        ifoa(nf)=0
      enddo
      do iic=-3, 3              
       do iis=-3, 3
        do iid=-3, 3
         do iiu=-3, 3
          do ii=1,mxlk
           call lkfokset(ii,4+iiu,4+iid,4+iis,4+iic, 0 )  
          enddo
         enddo
        enddo
       enddo
      enddo

      do i=1,nspecs+1
        id=ispecs(i)
        call idtr7(id,jc)
        do nf=1,nflav
        ifok(nf,i)=jc(nf,1)-jc(nf,2)
        ifoa(nf)=ifoa(nf)+iabs(ifok(nf,i))
        jspecs(1,nf,i)=jc(nf,1)
        jspecs(2,nf,i)=jc(nf,2)
        enddo
        iiu=ifok(1,i)
        iid=ifok(2,i)
        iis=ifok(3,i)
        iic=ifok(4,i)  !-charm
        if(abs(iiu).gt.3)stop'HNBSPD: u-dimension of lkfok too small'
        if(abs(iid).gt.3)stop'HNBSPD: d-dimension of lkfok too small'
        if(abs(iis).gt.3)stop'HNBSPD: s-dimension of lkfok too small'
        if(abs(iic).gt.3)stop'HNBSPD: c-dimension of lkfok too small'  
        if(ifok(5,i).ne.0)stop'HNBSPD: lkfok needs index for b'
        if(ifok(6,i).ne.0)stop'HNBSPD: lkfok needs index for t'
        call lkfokincrement(1,4+iiu,4+iid,4+iis,4+iic)            
        call lkfokget(1,4+iiu,4+iid,4+iis,4+iic, ii )            
        if(ii.gt.mxlk)stop'HNBSPD: ii-dimension of lkfok too small'
        call lkfokset(1+ii,4+iiu,4+iid,4+iis,4+iic, i )           
      enddo

      end

c-------------------------------------------------------------
      subroutine specs_getrid(id)
c-------------------------------------------------------------
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cspecs2/pspecs(0:mspecs+1),yie2volu,yiespecs
      r=rangen()
      do i=1,nspecs+1
        !print*,i,r,pspecs(i),ispecs(i)
        if(pspecs(i).ge.r)then
          id=ispecs(i)
          !print*,'------->',id
          return 
        endif
      enddo
      stop'ERROR 24062018b'
      end

c-------------------------------------------------------------
      subroutine specs_getp(id,p)
c-------------------------------------------------------------
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cspecs2/pspecs(0:mspecs+1),yie2volu,yiespecs
      if(id.eq.0)then
        p=pspecs(nspecs+1)-pspecs(nspecs)
        return
      endif
      do i=1,nspecs
        if(ispecs(i).eq.id)then
          p=pspecs(i)-pspecs(i-1)
          return
        endif
      enddo
      stop'ERROR 24062018'
      end

c-------------------------------------------------------------
      subroutine hnbgraee(T,ee)
c-------------------------------------------------------------
c input: 
c    T = temperature in GeV
c output: 
c    ee = Energy density of grand canonical HG for muB = 0 in GeV/fm3
c -------------------------------------------------------
c using (with a = g / (2*pi*hbar)^3 ):
c    dn = a*V * exp(-E/T) d3p
c the energy density is:
c    ee = \sum  V^-1 * \int E dn
c We find the formula below by using (from Textbooks) with z=am/T :
c K1(z) = z * \int_1^\infty (x^2-1)^0.5 * exp(-zx) dx
c K2(z) = 1/3 * z^2 * \int_1^\infty (x^2-1)^1.5 * exp(-zx) dx 
c --------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      hbarc=0.197327
       ee=0
        do i=1,nspecs
          g=gspecs(i)
          am=aspecs(i)
          a=g/(hbarc*2*pi)**3
          ee=ee + a * 4 * pi * am * am * T
     .       * ( 3 * T * sbessk(2,am/T) + am * sbessk1(am/T) )
        enddo
        end

c-------------------------------------------------------------
      subroutine hnbgra
c-------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incho'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cspecs2/pspecs(0:mspecs+1),yie2volu,yiespecs
      real aff(0:mspecs+1)
      parameter (nbins=500)
      real far(6,2,nbins/5)
      real ysumi(6,2)
      integer iplots(6,2),jplots(6,2)
      data iplots/ 990120, 990130, 991120, 992130, 992330, 993331
     .           ,-990120,-990130,-991120,-992130,-992330,-993331 /
      data jplots/ 120, 130, 1120, 2130, 2330, 3331
     .           ,-120,-130,-1120,-2130,-2330,-3331 / 

c determine T corresponding to a HG with given energy density efrout
      hbarc=0.197327
      mxee=25
      a=0.010
      b=0.400
      do i=1,mxee !bisection method for solving eps(T)=efrout
        call hnbgraee(a,fa) !returns eps=fa for given T=a
        fa=fa-efrout
        call hnbgraee(b,fb)
        fb=fb-efrout
        c=(a+b)/2
        call hnbgraee(c,fc)
        fc=fc-efrout
        !print*,a,b,c,'     ',fa,fb,fc,'    ',efrout
        if(fc.eq.0)goto 77
        if(fa*fc.lt.0)then
          b=c
        else
          a=c
        endif    
        if(b-a.lt.1e-5)goto 77
      enddo
  77  continue
      T=(a+b)/2
      call hnbgraee(T,ee)

c determine hadron yields / V acc to HG with temp T (muB=0)

      call gaulag(xlag,wlag,nlag,0.)
      do n=1,nlag
      wlag(n)=wlag(n)*exp(xlag(n))
      enddo
      eesum=0
      aff(0)=0
      do i=1,nspecs
        g=gspecs(i)
        am=aspecs(i)
        a=g/(hbarc*2*pi)**3
        fsum=0
        esum=0
        do n=1,nlag
          x=xlag(n)  ! p/T
          e=sqrt(am**2+x**2*T**2)
          w=exp(-sqrt(am**2/T**2+x**2)) 
          fsum=fsum+wlag(n)*x**2*w
          esum=esum+wlag(n)*x**2*w*e
        enddo
        fsum=fsum * a * 4 * pi * T**3  
        esum=esum * a * 4 * pi * T**3 
        aff(i)=aff(i-1)+fsum
        eesum=eesum+esum
        !write(*,'(i5,3f10.5)')i,esum,fsum,aff(i)
      enddo
      yie2volu=aff(nspecs) !yield per volume
      write(ifmt,'(a,4f9.5,f9.3)')
     .'hnbgra: T eps(3x) yie/V:',T,ee,eesum,efrout,yie2volu
       
      wzero=aff(nspecs) !yield / V  of zero (per def) 
      aff(nspecs+1)=aff(nspecs)+wzero !nspecs+1 = zero ptl
      pspecs(0)=0
      do i=1,nspecs+1
        pspecs(i)=aff(i)/aff(nspecs+1)
      enddo
      pspecs(nspecs+1)=1 !to avoid numerical error

      !do i=1,nspecs+1
      !  id=0 
      !  if(i.le.nspecs)id=ispecs(i)
      !  call specs_getp(id,p)
      !  write(ifmt,'(2i9,2f10.2)')i,id,(pspecs(i)-pspecs(i-1))*1000
      !.   ,p*1000
      !enddo
      
      if(iappl.ne.4)return

c-prepare plot and check for given volu and facmicro
      del=3.0/nbins
      ee=0
      yy=0
      do i=1,nspecs
        g=gspecs(i)
        am=aspecs(i)
        a=g/(hbarc*2*pi)**3
        esum=0
        ysum=0
        do j=1,nbins
          p=(j-0.5)*del
          e=sqrt(p*p+am*am)
          f=a*volu*4*pi*p**2*exp(-e/T)  * facmicro ! facmicro to compare different energies
          if(mod(j,5).eq.0)then
            do n=1,6
            do m=1,2
              if(ispecs(i).eq.iplots(n,m)
     .       .or.ispecs(i).eq.jplots(n,m))far(n,m,j/5)=f
            enddo 
            enddo
          endif
          esum=esum+f*e
          ysum=ysum+f
        enddo
        esum=esum*del
        ysum=ysum*del
        do n=1,6
        do m=1,2
          if(ispecs(i).eq.iplots(n,m)
     .   .or.ispecs(i).eq.jplots(n,m))ysumi(n,m)=ysum
        enddo 
        enddo
        ee=ee+esum
        yy=yy+ysum
      enddo
      write(ifmt,'(a,i5,f8.3,3x,2f8.3)')
     .'hnbgra: nspecs eps yie/V(2x): ',nspecs, ee/(volu*facmicro)
     .,aff(nspecs),yy/(volu*facmicro)
      
c-plot    
      do n=1,6
        write(ifhi,'(a,i1,$)')'openhisto name testp',n
        write(ifhi,'(a,$)')' xrange 0 3 yrange 1e-4 30  htyp lkv'
        write(ifhi,'(a)')' xmod lin ymod log'
        write(ifhi,'(a)')'text 0 0 ""xaxis p   "" '
        write(ifhi,'(a)')'text 0 0 ""yaxis dn/dp  ""' 
        write(ifhi,'(a)')'histoweight 1'
        write(ifhi,'(a)')'array  2'
        do j=1,nbins/5
          p=(j*5-0.5)*del
          write(ifhi,'(2f12.5)') p,far(n,1,j)+far(n,2,j)
        enddo
        xpo=-0.28+n*0.20 
        !print*,esum,ee,ii,xpo
        write(ifhi,'(a)')'endarray'
        write(ifhi,'(a,f4.2,a,f8.3,a)')
     .  'text ',xpo,' 1.02 "" ',ysumi(n,1)+ysumi(n,2),'  "" '
        write(ifhi,'(a,$)')'closehisto plot 0'
        if(n.ne.6)write(ifhi,'(a)')'-'
        if(n.eq.6)write(ifhi,'(a)')' '
      enddo
      end

c-------------------------------------------------------------
      subroutine idtr7(id,jc)
c-------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2),ifl(3)
      do i=1,nflav
        jc(i,1)=0
        jc(i,2)=0
      enddo
      if(id.eq.0)return
      if(iabs(id).gt.16.and.iabs(id).lt.20)then
        if(id.eq.17)then
          jc(1,1)=3
          jc(2,1)=3
        elseif(id.eq.-17)then
          jc(1,2)=3
          jc(2,2)=3
        elseif(id.eq.18)then
          jc(1,1)=4
          jc(2,1)=5
        elseif(id.eq.-18)then
          jc(1,2)=4
          jc(2,2)=5
        elseif(id.eq.19)then
          jc(1,1)=6
          jc(2,2)=6
        elseif(id.eq.-19)then
          jc(1,2)=6
          jc(2,2)=6
        endif
        return
      endif
      if(id.eq.30)then
         jc(1,1)=2
         jc(2,1)=2
         jc(3,1)=2
         return
      elseif(id.eq.-30)then
         jc(1,2)=2
         jc(2,2)=2
         jc(3,2)=2
         return
      endif
      if(abs(id).gt.990000)then
        idxx=110
      else
        idxx=id
      endif
      call idflav(idxx,ifl(1),ifl(2),ifl(3),jspin,idu)
      do j=1,3
         if(ifl(j).gt.0)jc(ifl(j),1)=jc(ifl(j),1)+1
         if(ifl(j).lt.0)jc(-ifl(j),2)=jc(-ifl(j),2)+1
      enddo
      end

c-------------------------------------------------------------
      subroutine idtr8(id,ju,jd,js)
c-------------------------------------------------------------
      integer ifl(3)
      ju=0
      jd=0
      js=0 
      if(id.eq.0)return
      if(iabs(id).gt.16.and.iabs(id).lt.20)then
        if(abs(id).eq.17)then
          ju=sign(3,id)
          jd=sign(3,id)
        elseif(abs(id).eq.18)then
          ju=sign(4,id)
          jd=sign(5,id)
        elseif(abs(id).eq.19)then
          ju=sign(6,id)
          jd=sign(6,id)
        endif
        return
      endif
      if(abs(id).eq.30)then
         ju=sign(2,id)
         jd=sign(2,id)
         js=sign(2,id)
         return
      endif
      if(abs(id).gt.990000)then
        idxx=110
      else
        idxx=id
      endif
      call idflav(idxx,ifl(1),ifl(2),ifl(3),jspin,idu)
      do j=1,3
         if(abs(ifl(j)).eq.1)ju=ju+sign(1,ifl(j))
         if(abs(ifl(j)).eq.2)jd=jd+sign(1,ifl(j))
         if(abs(ifl(j)).eq.3)js=js+sign(1,ifl(j))
      enddo
      end

c----------------------------------------------------------------------
      subroutine hnbspi(id,spideg)
c----------------------------------------------------------------------
c  returns spin degeneracy spideg for particle id-code id
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter (nspec=62)
      dimension ispec(nspec),spid(nspec)
      data ispec/
     *     1,   -1,    2,   -2,    3,   -3
     *,  110,  120, -120,  220,  130, -130,  230, -230,  330
     *,  111,  121, -121,  221,  131, -131,  231, -231,  331
     *, 1120, 1220, 1130, 2130, 1230, 2230, 1330, 2330
     *, 1111, 1121, 1221, 2221, 1131, 1231, 2231, 1331, 2331, 3331
     *,-1120,-1220,-1130,-2130,-1230,-2230,-1330,-2330
     *,-1111,-1121,-1221,-2221,-1131,-1231,-2231,-1331,-2331,-3331
     *,441,30/
      data spid/
     *  6*6.
     *, 9*1.
     *, 9*3.
     *, 8*2.
     *,10*4.
     *, 8*2.
     *,10*4.
     *,1*3
     *,1*3/
      do i=1,nspec
      if(id.eq.ispec(i))then
      spideg=spid(i)
c      fac=1
      !factb ... not used
      !factq ... not used
c      call idflav(id,ifl1,ifl2,ifl3,jspin,index)
c      ifls=0
c      iflb=0
c      if(abs(ifl1).eq.3)ifls=ifls+1
c      if(abs(ifl2).eq.3)ifls=ifls+1
c      if(abs(ifl3).eq.3)ifls=ifls+1
c      if((i.gt.10.and.i.lt.15).or.(i.gt.19.and.i.lt.24))then
c        fac=facts               !reduce kaons (increase lambda)
c      elseif(i.eq.60)then
c        fac=factb               !reduce Omega
c      endif
      !-----------------
ckw2018      spideg=spideg*fac 
ckw2018 fac must be 1 to correpond to microcanonical hadronization
      !-----------------
      goto 1
      endif
      enddo
      call utstop('hnbspi: id not found&')
    1 continue
      return
      end

c----------------------------------------------------------------------
      subroutine hnbtst(iof12)
c----------------------------------------------------------------------
c  calculates logs of prefactors and phase space integral
c  for ultrarelativistic limit (massless particles) and (2*s_i+1)=1
c  f12log and w15log=w35log+f12log not calculated calculated for iof12=0
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      common/ctst/psulog,wtulog
      integer ii(maxp)
      common /clatt/nlattc,npmax

      pi=3.1415927
      hquer=0.197327
      ish0=ish
      if(ishsub/100.eq.23)ish=mod(ishsub,100)
      do i=1,np
      ii(i)=1
      enddo

      if(ish.ge.7)write(ifch,*)('-',i=1,10)
     *,' entry sr hnbtst ',('-',i=1,30)
      if(ish.ge.7)write(ifch,*)'configuration:'
      if(ish.ge.7)write(ifch,*)(ident(i),i=1,np)
      if(ish.ge.7)write(ifch,*)'n_l:',nlattc,'   n_0:',nlattc-np

c log of prod m_i*volu/4/pi**3/hquer**3 -> f5log
      f5log=0
      do i=1,np
      call hnbfaf(i,gg,am)
      f5log=f5log+alog(gg*am*volu/4/pi**3/hquer**3)
      enddo
      if(ish.ge.7)write(ifch,*)'log(f5):',f5log

c log f4log=0
      f4log=0
      if(ish.ge.7)write(ifch,*)'log(f4):',f4log

c log of 1/prod n_alpha! -> f3log
      dbllog=0
      n1=1
      nx=1
    1 continue
      i=0
      x=0
      do n2=n1,np
      if(ident(n2).eq.ident(n1))then
      ii(n2)=0
      i=i+1
      x=x+alog(i*1.)
      endif
      if(ii(n2).ne.0.and.n2.gt.n1.and.nx.eq.n1
     *.and.ident(n2).ne.ident(n1))nx=n2
      enddo
      dbllog=dbllog+x
      if(nx.gt.n1)then
      n1=nx
      goto1
      endif
      f3log=-dbllog
      if(ish.ge.7)write(ifch,*)'log(f3):'
     *,f3log

c log of f3 * f4 * f5
      f35log=f5log+f4log+f3log
      if(ish.ge.7)write(ifch,*)'log(f3*f4*f5):',f35log

c log of phase space integral --> psilog
      psilog=alog(2.*np*np*(np-1)/tecm**4/pi)
      do i=1,np
      psilog=psilog+alog(tecm**2*pi/2./i/i)
      enddo
      if(ish.ge.7)write(ifch,*)'method 1 log(psi):',psilog
      psilog=-alog(2.*np-1)
      psilog=psilog+(np-1)*alog(pi/2.)
      do i=1,2*np-2
      psilog=psilog+alog((2.*np+i-2)/i)
      enddo
      do i=1,3*np-4
      psilog=psilog+alog(tecm/i)
      enddo
      if(ish.ge.7)write(ifch,*)'method 2 log(psi):',psilog

c log of phase space integral * f3 * f4 * f5
      w35log=f35log+psilog
      if(ish.ge.7)write(ifch,*)'log(f35*psi):',w35log

           if(iof12.ne.0)then

c log of macro/micro factor (f1*f2) --> f12log
      deglog=0
      do i=1,np
      deglog=deglog+alog(1.*i)
      enddo
      deglog=deglog+f3log
      do i=1,np
      deglog=deglog+alog(nlattc+1.-i)-alog(1.*i)
      enddo
      f12log=-deglog

      w15log=w35log+f12log
      if(ish.ge.7)then
        write(ifch,*)'log(f1*f2):',f12log
        write(ifch,*)'log(f15*psi):',w15log
        write(ifch,'(1x,4(a,3x))')
     *'log(fac):','log(psi):',' log(wt):','log(wta):'
        write(ifch,'(1x,4(f9.3,3x))')
     *f12log+f35log,psilog,w15log,w15log-f12log
      endif

           endif

      psulog=psilog
      wtulog=w35log

      if(ish.ge.7)write(ifch,*)('-',i=1,30)
     *,' exit sr hnbtst ',('-',i=1,10)
      ish=ish0
      return
      end



c#######################################################################
c#######################################################################
c#######################################################################
c###########################  xhnbroutines ############################# 
c#######################################################################
c#######################################################################
c#######################################################################


c-------------------------------------------------------------
      subroutine xhnbspf(ku,kd,ks,kc,kb,kt,j,n,spelog)
c-------------------------------------------------------------
c  returns spelog = log of factor for consid. different species
c  spelog is double precision
c  option ioflac determines the method:
c     ioflac=1: ignore flavour conservation
c     ioflac=2: flavour conservation implemented straightforward
c                 (only for nspecs=3,7)
c     ioflac=3: flavour conservation via generating fctn
c  further input:
c     ku,...,kt (integer) : flavour
c     j (integer) : excluded species
c     n (integer) : multiplicity
c-------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)
      integer m(7),l(7),ifot(nflav)
      common/csph/ifox(nflav),ifoy(nflav),jx,nx,ifom(nflav,mspecs)
      parameter(mxfacu=200)
      double precision faci(0:mxfacu)
      double precision utgam2,spelog,spe
c      parameter(numax=100,kqmax=100)
c      parameter(mxhh=200)
      if(ish.ge.9)write(ifch,*)('-',i=1,10)
     *,' entry sr xhnbspf ',('-',i=1,30)
      if(ish.ge.9)write(ifch,'(1x,a,9x,a,4x,a)')
     *' ku kd ks kc kb kt','j','n'
      if(ish.ge.9)write(ifch,'(1x,6i3,5x,2i5)')
     *ku,kd,ks,kc,kb,kt,j,n
      k=nspecs
      jx=j
      nx=n
      ifot(1)=ku
      ifot(2)=kd
      ifot(3)=ks
      ifot(4)=kc
      ifot(5)=kb
      ifot(6)=kt

           if(ioflac.eq.1)then

      if(ish.ge.9)write(ifch,'(1x,a,i1)')'ioflac=',ioflac
      g=0
      do i=1,nspecs
      if(i.ne.j)g=g+gspecs(i)
      enddo
      spelog=n*dlog(1.d0*g)

           elseif(ioflac.eq.2)then

      if(ish.ge.9)write(ifch,'(1x,a,i2)')'ioflac:',ioflac
           if(k.eq.3)then
      if(ish.ge.9)write(ifch,'(1x,a,i2)')'nspecs:',nspecs
      spe=0d0
           if(j.lt.1.or. j.gt.k)then
      do 1 n1=0,n
      do 2 n2=0,n-n1
      n3=n-n1-n2
      do 5 nf=1,nflav
      if(ifoa(nf).eq.0.and.ifot(nf).eq.0)goto5
      if(n1*ifok(nf,1)+n2*ifok(nf,2)+n3*ifok(nf,3).ne.ifot(nf))goto2
    5 continue
      spe=spe+utgam2(1.d0+n)
     &/utgam2(1.d0+n1)/utgam2(1.d0+n2)/utgam2(1.d0+n3)
     &*gspecs(1)**n1*gspecs(2)**n2*gspecs(3)**n3
    2 continue
    1 continue
           else
      do 3 i1=0,n
      i2=n-i1
      m(1)=0
      m(2)=i1
      m(3)=i2
      do i=1,3
      ii=1+mod(j-2+i,3)
      l(ii)=m(i)
      enddo
      n1=l(1)
      n2=l(2)
      n3=l(3)
      do 6 nf=1,nflav
      if(ifoa(nf).eq.0.and.ifot(nf).eq.0)goto6
      if(n1*ifok(nf,1)+n2*ifok(nf,2)+n3*ifok(nf,3).ne.ifot(nf))goto3
    6 continue
      spe=spe+utgam2(1.d0+n)
     &/utgam2(1.d0+n1)/utgam2(1.d0+n2)/utgam2(1.d0+n3)
     &*gspecs(1)**n1*gspecs(2)**n2*gspecs(3)**n3
    3 continue
           endif
      if(ish.ge.9)write(ifch,*)'spe:',spe
      spelog=-1000
      if(spe.gt.0.d0)spelog=dlog(spe)
      if(ish.ge.9)write(ifch,*)'spelog:',spelog
           elseif(k.eq.7)then
      if(ish.ge.9)write(ifch,'(1x,a,i2)')'nspecs:',nspecs
      if(n.gt.mxfacu)call utstop('xhnbspf: mxfacu too small&')
      do lf=0,n
      faci(lf)=1.d0/utgam2(1d0+lf)
      enddo
      spe=0
           if(j.lt.1.or. j.gt.k)then
      do n1=0,n
      do n2=0,n-n1
      do n3=0,n-n1-n2
      do n4=0,n-n1-n2-n3
      do n5=0,n-n1-n2-n3-n4
      do 12 n6=0,n-n1-n2-n3-n4-n5
      n7=n-n1-n2-n3-n4-n5-n6
      do 15 nf=1,nflav
      if(ifoa(nf).eq.0.and.ifot(nf).eq.0)goto15
      if(n1*ifok(nf,1)+n2*ifok(nf,2)+n3*ifok(nf,3)+n4*ifok(nf,4)
     *+n5*ifok(nf,5)+n6*ifok(nf,6)+n7*ifok(nf,7).ne.ifot(nf))goto12
   15 continue
      spe=spe+1d0/faci(n)*faci(n1)*faci(n2)*faci(n3)*faci(n4)
     &*faci(n5)*faci(n6)*faci(n7)
     &*gspecs(1)**n1*gspecs(2)**n2*gspecs(3)**n3*gspecs(4)**n4
     &*gspecs(5)**n5*gspecs(6)**n6*gspecs(7)**n7
   12 continue
      enddo
      enddo
      enddo
      enddo
      enddo
           else
      do i1=0,n
      do i2=0,n-i1
      do i3=0,n-i1-i2
      do i4=0,n-i1-i2-i3
      do 13 i5=0,n-i1-i2-i3-i4
      i6=n-i1-i2-i3-i4-i5
      m(1)=0
      m(2)=i1
      m(3)=i2
      m(4)=i3
      m(5)=i4
      m(6)=i5
      m(7)=i6
      do i=1,7
      ii=1+mod(j-2+i,7)
      l(ii)=m(i)
      enddo
      n1=l(1)
      n2=l(2)
      n3=l(3)
      n4=l(4)
      n5=l(5)
      n6=l(6)
      n7=l(7)
      do 16 nf=1,nflav
      if(ifoa(nf).eq.0.and.ifot(nf).eq.0)goto16
      if(n1*ifok(nf,1)+n2*ifok(nf,2)+n3*ifok(nf,3)+n4*ifok(nf,4)
     *+n5*ifok(nf,5)+n6*ifok(nf,6)+n7*ifok(nf,7).ne.ifot(nf))goto13
   16 continue
      spe=spe+1d0/faci(n)*faci(n1)*faci(n2)*faci(n3)*faci(n4)
     &*faci(n5)*faci(n6)*faci(n7)
     &*gspecs(1)**n1*gspecs(2)**n2*gspecs(3)**n3*gspecs(4)**n4
     &*gspecs(5)**n5*gspecs(6)**n6*gspecs(7)**n7
   13 continue
      enddo
      enddo
      enddo
      enddo
           endif
      if(ish.ge.9)write(ifch,*)'spe:',spe
      spelog=-1000
      if(spe.gt.0.d0)spelog=dlog(spe)
      if(ish.ge.9)write(ifch,*)'spelog:',spelog
           else
      call utstop('xhnbspf: ioflac=2 only for nspecs=3,7&')
           endif

           elseif(ioflac.eq.3)then

      call utstop('xhnbspf: ioflac must be 1 or 2&')

           endif

      if(ish.ge.9)write(ifch,*)('-',i=1,30)
     *,' exit sr xhnbspf ',('-',i=1,10)
      return
      end

c-------------------------------------------------------------
      subroutine xhnbspg(ku,kd,ks,kc,kb,kt,j,n,spelog)
c-------------------------------------------------------------
      include 'epos.inc'
      double precision spelog,spalog
      if(ioflac.ne.0)return
      ioflac=2
      call xhnbspf(ku,kd,ks,kc,kb,kt,j,n,spalog)
      ioflac=3
      call xhnbspf(ku,kd,ks,kc,kb,kt,j,n,spelog)
      ioflac=0
      write(ifch,*)'ioflac=2/3:',spalog,spelog
      return
      end

c----------------------------------------------------------------------
      subroutine xhnbwri
c----------------------------------------------------------------------
c  writes (to ifch) an configuration
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      common/cfact/faclog
      write(ifch,'(1x,a,i5)')'np:',np
      write(ifch,'(1x,3(a,3x))')
     *'log(fac):','log(psi):',' log(wt):'
      if(wtlog.gt.-1e30.and.wtxlog.gt.-1e30)then
      write(ifch,'(1x,3(f9.3,3x))')faclog,wtxlog,wtlog
      else
      write(ifch,*)faclog,wtxlog,wtlog
      endif
      if(np.le.1)return
      call hnbtst(1)
      write(ifch,*)'particle id codes:'
      write(ifch,'(1x,10i6)')(ident(n),n=1,np)
      write(ifch,*)'particle masses:'
      write(ifch,'(1x,10f6.3)')(amass(n),n=1,np)
      end

c----------------------------------------------------------------------
      subroutine xhnbzen(iii)
c----------------------------------------------------------------------
c analysis of events. energy spectra.
c for iii>0: filling histogram considering ptl iii
c----------------------------------------------------------------------
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      parameter (nhise=1) !(nhise=100)
      common/chise/hise(mspecs,nhise)
      de=2./nhise/2.

      j=0

           if(iii.gt.0)then

      i=iii
      do l=1,nspecs
      if(ident(i).eq.ispecs(l))then
      j=l
      goto1
      endif
      enddo
    1 continue
      am=aspecs(j)
      e=pcm(4,i)
      ke=1+int((e-am)/(2*de))
      if(ke.ge.1.and.ke.le.nhise)hise(j,ke)=hise(j,ke)+1
      return

           else

      stop'STOP in xhnbzen: iii=0'

           endif

      end

c----------------------------------------------------------------------
      subroutine xhnbzmu(iii)
c----------------------------------------------------------------------
c analysis of events. multiplicity spectra.
c for iii<0: settting histograms to zero (should be first call)
c for iii>0: filling histogram considering ptl iii
c----------------------------------------------------------------------
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      parameter (nhismu=1) !(nhismu=500)
      common/chismu/hismu(mspecs,0:nhismu),hismus(nhismu)

           if(iii.lt.0)then

      do i=1,nhismu
      hismus(i)=0
      enddo
      do j=1,nspecs
      do i=0,nhismu
      hismu(j,i)=0
      enddo
      enddo
      goto1000

           elseif(iii.gt.0)then

      if(np.ge.1.and.np.le.nhismu)hismus(np)=hismus(np)+1
      do j=1,nspecs
      mu=0
      do i=1,np
      if(ident(i).eq.ispecs(j))mu=mu+1
      enddo
      if(mu.ge.0.and.mu.le.nhismu)hismu(j,mu)=hismu(j,mu)+1
      enddo
      goto1000

           else

      stop'STOP in sr xhnbzmu: iii must not be 0'

           endif

1000  continue
      return
      end


c-----------------------------------------------------------------------
      subroutine xhnben
c-----------------------------------------------------------------------
c produces histogram of energy spectrum (after metropolis run)
c complete histogram: openhisto ... closehisto
c iocite=1 required
c-----------------------------------------------------------------------
c xpar1: particle species (venus id-code)
c xpar2: 1: actual spectrum 2: fit
c xpar3: 1: de/d3p 2: ede/d3e
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      parameter (nhise=1) !(nhise=100)
      common/chise/hise(mspecs,nhise)
      parameter (literm=500)
      common/cmet/kspecs(mspecs),liter
     *,iterl(literm),iterc(literm)
      real datx(nhise),daty(nhise),dats(nhise)
      common/citer/iter,itermx
      character ch*1,chid*5,cyield*9,ctem*5
      de=2./nhise/2.

      if(iocite.ne.1)stop'STOP: xhnben: iocite=1 required'

      idcode=nint(xpar1)
      mode=nint(xpar2)
      kind=nint(xpar3)

           do j=1,nspecs
           if(idcode.eq.ispecs(j))then

      id=idcode
      am=aspecs(j)
      yield=1.*kspecs(j)/(itermx-iternc)
      if(kind.eq.1)ch=' '
      if(kind.eq.2)ch='e'
      ll=kind-1
      e0=am+de
      nebins=0
        do i=1,nhise
      e=e0+(i-1)*2*de
      p1=sqrt((e-de)**2-am**2)
      p2=sqrt((e+de)**2-am**2)
      d3p=4*pi*(p2**3-p1**3)/3
      datx(i)=e
      y=(1-ll+ll*e)*hise(j,i)/(itermx-iternc)/d3p
      if(y.gt.0.)then
      nebins=nebins+1
      daty(i)=alog(y)
      d=y/sqrt(hise(j,i))
      dats(i)=1e10
      if(y-d.gt.0.)dats(i)=alog(y+d)-alog(y-d)
      else
      daty(i)=-100
      dats(i)=1e10
      endif
c-c   if(e.lt.0.2)dats(i)=1e10
        enddo
      a=0.
      b=0.
        if(nebins.ge.3)then
      call utfit(datx,daty,nhise,dats,1,a,b,siga,sigb,chi2,q)
      tem=-1./b
      if(tem.lt.0.050.or.tem.gt.10.)then
      tem=0.
      a=0.
      b=0.
      endif
        endif
      do i=1,nhise
      daty(i)=exp(daty(i))
      enddo
      write(chid,'(i5)')id
      write(cyield,'(f9.4)')yield
      ctem='     '
      if(tem.gt.0.)write(ctem,'(f5.3)')tem
      write(ifhi,'(a)')    'openhisto xrange 0 3'
      write(ifhi,'(a)')    'htyp lin xmod lin ymod log'
      write(ifhi,'(a,a)')  'text 0 0 "title id='//chid
     *                           ,'   N='//cyield//'   T='//ctem//'"'
      write(ifhi,'(a)')    'text 0 0 "xaxis energy (GeV)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis '//ch//' dn/d3p (GeV-3)"'
      write(ifhi,'(a)')    'array 2'
      do i=1,nhise
      if(mode.eq.1)write(ifhi,'(2e12.4)')datx(i),daty(i)
      if(mode.eq.2)write(ifhi,'(2e12.4)')datx(i),exp(a+b*datx(i))
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto'

           endif
           enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine xhnbit
c-----------------------------------------------------------------------
c produces histogram of multiplicity versus iterations (after metropolis run)
c complete histogram: openhisto ... closehisto
c iocite=1 required
c-----------------------------------------------------------------------
c xpar1: particle species (0=all, else venus id-code)
c xpar2: 1:actual multiplicity 2:average multiplicity 3:grand canonical
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      parameter (literm=500)
      common/cmet/kspecs(mspecs),liter
     *,iterl(literm),iterc(literm)
      real datlx(literm),datly(literm)
      common/citer/iter,itermx
      character chid*5,ctecm*5,cvolu*6
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      common/cgctot/rmstot,ptltot

      if(iocite.ne.1)stop'STOP: xhnbit: iocite=1 required'

      idcode=nint(xpar1)
      mode=nint(xpar2)

           if(idcode.eq.0)then

      yield=0
      do j=1,nspecs
      yield=yield+1.*kspecs(j)/(itermx-iternc)
      enddo
      datlx(1)=(iterl(1)+1)/2.
      do li=2,liter-1
      datlx(li)=(iterl(li)+iterl(li-1)+1)/2.
      enddo
      x1=0
      x2=iterl(liter-1)
      do li=1,liter-1
      y=0
      do j=1,nspecs
      call lspecsget(li,j,ival)
      y=y+ival
      enddo
      if(mode.eq.1)datly(li)=y/iterc(li)
      if(mode.eq.2)datly(li)=yield
      if(mode.eq.3)datly(li)=ptltot
      enddo
      write(ctecm,'(f5.1)')tecm
      write(cvolu,'(f6.1)')volu
      write(ifhi,'(a,2e11.3)')'openhisto xrange',x1,x2
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,a)')     'text 0 0 "title E = '//ctecm//'   V = '
     *                                 ,cvolu//'"'
      write(ifhi,'(a)')       'text 0 0 "xaxis iterations"'
      write(ifhi,'(a)')       'text 0 0 "yaxis multiplicity"'
      write(ifhi,'(a)')       'array 2'
      do i=1,liter-1
      write(ifhi,'(2e12.4)')   datlx(i),datly(i)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'

           else

           do j=1,nspecs
           if(idcode.eq.ispecs(j))then

      yield=1.*kspecs(j)/(itermx-iternc)
      write(chid,'(i5)')idcode
      do li=1,liter-1
      datlx(li)=iterl(li)
      enddo
      x1=0
      x2=datlx(liter-1)
      do li=1,liter-1
      if(mode.eq.1)then
        call lspecsget(li,j,ival)
        datly(li)=ival*1./iterc(li)
      endif 
      if(mode.eq.2)datly(li)=yield
      if(mode.eq.3)datly(li)=ptlngc(j)
      enddo
      write(ifhi,'(a,2e11.3)')'openhisto xrange',x1,x2
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')       'text 0 0 "title id='//chid//'"'
      write(ifhi,'(a)')       'text 0 0 "xaxis iterations "'
      write(ifhi,'(a)')       'text 0 0 "yaxis multiplicity"'
      write(ifhi,'(a)')       'array 2'
      do i=1,liter-1
      write(ifhi,'(2e12.4)')   datlx(i),datly(i)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'

           endif
           enddo

           endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xhnbmu
c-----------------------------------------------------------------------
c produces histogram of multiplicity distribution (after metropolis run)
c complete histogram: openhisto ... closehisto
c iocite=1 required
c-----------------------------------------------------------------------
c xpar1: particle species (0=all, else venus id-code)
c xpar2: xrange automatic (0) or given via xpar3,4 (else)
c xpar3,4: xrange
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      parameter (nhismu=1) !(nhismu=500)
      common/chismu/hismu(mspecs,0:nhismu),hismus(nhismu)
      parameter (literm=500)
      common/cmet/kspecs(mspecs),liter
     *,iterl(literm),iterc(literm)
      real datx(nhismu),daty(nhismu)
      common/citer/iter,itermx
      common /clatt/nlattc,npmax
      character chid*5,cyield*9,ctecm*5,cvolu*6

      if(iocite.ne.1)stop'STOP: xhnbmu: iocite=1 required'

      idcode=nint(xpar1)
      ixr=nint(xpar2)
      xx1=xpar3
      xx2=xpar4

      write(ctecm,'(f5.1)')tecm
      write(cvolu,'(f6.1)')volu

           if(idcode.eq.0)then

      yield=0
      do j=1,nspecs
      yield=yield+1.*kspecs(j)/(itermx-iternc)
      enddo
      write(cyield,'(f9.4)')yield
      i1=0
      i2=nlattc
      mus=0
      do i=1,nhismu
      if(i1.eq.0.and.nint(hismus(i)).gt.0)i1=i
      if(nint(hismus(i)).gt.0)i2=i
      mus=mus+hismus(i)
      enddo
      ij=0.5*(i1+i2)*0.20
      if(itermx.le.1000)ij=0.5*(i1+i2)*0.40
      if(itermx.le.100)ij=0.5*(i1+i2)*0.80
      i1=i1-ij
      i1=max(i1,2)
      i2=i2+ij
      ii=10
      if(i1.le.50)ii=5
      if(i1.le.20)ii=2
      i1=i1/ii*ii
      i2=i2/ii*ii+ii
           do i=i1,i2
      l=1+i-i1
      datx(l)=i
      daty(l)=hismus(i)/mus
           enddo
      jx=1+i2-i1
      if(ixr.eq.0)then
      x1=i1
      x2=i2
      else
      x1=xx1
      x2=xx2
      endif
      write(ifhi,'(a,2e11.3)')'openhisto xrange',x1,x2
      write(ifhi,'(a)')       'htyp lin xmod lin ymod log'
      write(ifhi,'(a,a)')     'text 0 0 "title E = '//ctecm//'   V = '
     *                              ,cvolu//'"'
      write(ifhi,'(a)')       'text 0 0 "xaxis multiplicity n  "'
      write(ifhi,'(a)')       'text 0 0 "yaxis dN/dn"'
      write(ifhi,'(a)')       'text 0.30 0.25 "N?MC!='//cyield//'"'
      write(ifhi,'(a)')       'array 2'
      do i=1,jx
      write(ifhi,'(2e12.4)')   datx(i),daty(i)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'

           else

           do j=1,nspecs
           if(idcode.eq.ispecs(j))then

      yield=1.*kspecs(j)/(itermx-iternc)
      write(cyield,'(f9.4)')yield
      write(chid,'(i5)')idcode
      i1=0
      i2=nlattc
      mus=0
      do i=0,nhismu
      if(i1.eq.0.and.nint(hismu(j,i)).gt.0)i1=i
      if(nint(hismu(j,i)).gt.0)i2=i
      mus=mus+hismu(j,i)
      enddo
      ij=0.5*(i1+i2)*0.30
      if(itermx.le.1000)ij=0.5*(i1+i2)*0.60
      if(itermx.le.100)ij=0.5*(i1+i2)*1.20
      i1=i1-ij
      i1=max(i1,0)
      i2=i2+ij
      ii=10
      if(i1.le.50)ii=5
      if(i1.le.20)ii=2
      i1=i1/ii*ii
      i2=i2/ii*ii+ii
           do i=i1,i2
      l=1+i-i1
      datx(l)=i
      daty(l)=hismu(j,i)/mus
           enddo
      jx=1+i2-i1
      if(ixr.eq.0)then
      x1=i1
      x2=i2
      else
      x1=xx1
      x2=xx2
      endif
      write(ifhi,'(a,2e11.3)')'openhisto xrange',x1,x2
      write(ifhi,'(a)')       'htyp lin xmod lin ymod log'
      write(ifhi,'(a)')       'text 0 0 "title id='//chid//'"'
      write(ifhi,'(a)')       'text 0 0 "xaxis multiplicity n  "'
      write(ifhi,'(a)')       'text 0 0 "yaxis dN/dn"'
      write(ifhi,'(a)')       'text 0.30 0.25 "N?MC!='//cyield//'"'
      write(ifhi,'(a)')       'array 2'
      do i=1,jx
      write(ifhi,'(2e12.4)')   datx(i),daty(i)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'

           endif
           enddo

           endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xhnbmz
c-----------------------------------------------------------------------
c produces histogram of multiplicity distribution from droplet decay
c or average multiplicity versus iterations
c for massless hadrons
c complete histogram: openhisto ... closehisto
c-----------------------------------------------------------------------
c xpar1: particle species (0=all, else venus id-code)
c xpar2: lower limit multiplicity
c xpar3: upper limit multiplicity
c xpar4: lower limit total multiplicity   (also necc for xpar1.ne.0)
c xpar5: upper limit  "      "            (also necc for xpar1.ne.0)
c xpar6: sets htyp: 1->lfu, 2->ldo, 3->lda, 4->ldd
c xpar7: 0: multiplicity distribution
c        >0: av multiplicity vs iterations (itermx=xpar7)
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      common/ctst/psulog,wtulog
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      parameter (nhismu=1) !(nhismu=500)
      common/cflac/ifok(nflav,mspecs+1),ifoa(nflav)
      real datx(nhismu),datyu(nhismu)
      character cyieur*9
      real pzlog(nhismu)
      double precision spelog,cc,bb,dsu
      common/cyield/yield
      character*3 htyp

      idcode=nint(xpar1)
      x1=xpar2
      x2=xpar3
      i1=nint(xpar2)
      i2=nint(xpar3)
      ii1=nint(xpar4)
      ii2=nint(xpar5)
      ih=nint(xpar6)
      htyp='lin'
      if(ih.eq.1)htyp='lfu'
      if(ih.eq.2)htyp='ldo'
      if(ih.eq.3)htyp='lda'
      if(ih.eq.4)htyp='ldd'
      itmax=nint(xpar7)
      pzlog(1)=0.

      wtrlog=-1e30
           do i=ii1,ii2
      if(i.ge.2)then
      np=i
      do k=1,np
      ident(k)=110
      enddo
      call hnbtst(0)
      wtzlog=wtulog
      if(ioflac.eq.0)call xhnbspg(keu,ked,kes,kec,keb,ket,0,np,spelog)
      if(ioflac.ne.0)call xhnbspf(keu,ked,kes,kec,keb,ket,0,np,spelog)
      wtulog=wtulog+spelog
      else
      wtzlog=-1000
      wtulog=-1000
      endif
      pzlog(1+i-ii1)=wtzlog
      datyu(1+i-ii1)=wtulog
      wtrlog=max(wtrlog,wtulog)
           enddo
      yield=0
      su=0
           do i=ii1,ii2
      l=1+i-ii1
      pzlog(l)=pzlog(l)-wtrlog
      datyu(l)=datyu(l)-wtrlog
      if(datyu(l).gt.-50.)then
      datyu(l)=exp(datyu(l))
      else
      datyu(l)=exp(-50.)
      endif
      yield=yield+i*datyu(l)
      su=su+datyu(l)
           enddo
      yield=yield/su
           do i=ii1,ii2
      l=1+i-ii1
      datx(l)=i
      datyu(l)=datyu(l)/su
           enddo
      jx=1+ii2-ii1
      write(cyieur,'(f9.4)')yield
c     ---
        if(idcode.eq.0.and.itmax.eq.0)then
      write(ifhi,'(a,2e11.3)')'openhisto xrange',x1,x2
      write(ifhi,'(a)')       'htyp '//htyp//' xmod lin ymod log'
      write(ifhi,'(a)')       'text 0.30 0.15 "N?ana!='//cyieur//'"'
      write(ifhi,'(a)')       'array 2'
      do i=1,jx
      write(ifhi,'(2e12.4)')   datx(i),datyu(i)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'
        elseif(idcode.eq.0)then
      write(ifhi,'(a,2e11.3)')'openhisto xrange',0.,itmax*1.
      write(ifhi,'(a)')       'htyp '//htyp//' xmod lin ymod lin'
      write(ifhi,'(a)')       'array 2'
      itm=20
      do i=1,itm
      write(ifhi,'(2e12.4)')   (i-1.)*itmax/(itm-1.),yield
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'
        endif
c     ---
      if(idcode.eq.0)return

           do j=1,nspecs
           if(idcode.eq.ispecs(j))then

      wtrlog=-1e30
           do i=i1,i2
      l=1+i-i1
      datx(l)=i
           enddo
      yield=0
      suj=0
      dsu=su
           do i=i1,i2
      l=1+i-i1
      bb=0
      nfi=0
      do ntot=max(i+1,ii1),min(i2*nspecs,ii2)
      nfi=nfi+1
      cc=1d0
      do kc=1,i
      cc=cc*(1.+ntot-kc)/kc*gspecs(j)
      enddo
      ku=keu-i*ifok(1,j)
      kd=ked-i*ifok(2,j)
      ks=kes-i*ifok(3,j)
      kc=kec-i*ifok(4,j)
      kb=keb-i*ifok(5,j)
      kt=ket-i*ifok(6,j)
      if(ioflac.eq.0)call xhnbspg(ku,kd,ks,kc,kb,kt,j,ntot-i,spelog)
      if(ioflac.ne.0)call xhnbspf(ku,kd,ks,kc,kb,kt,j,ntot-i,spelog)
      cc=cc*dexp(spelog)
      bb=bb+cc*dexp(1.d0*pzlog(1+ntot-ii1))/dsu
      enddo
      datyu(l)=bb
      yield=yield+i*datyu(l)
      suj=suj+datyu(l)
           enddo
      yield=yield/suj
      jx=1+i2-i1
      write(cyieur,'(f9.4)')yield
c     ---
        if(itmax.eq.0)then
      write(ifhi,'(a,2e11.3)')'openhisto xrange',x1,x2
      write(ifhi,'(a)')       'htyp '//htyp//' xmod lin ymod log'
      write(ifhi,'(a)')       'text 0.30 0.15 "N?ana!='//cyieur//'"'
      write(ifhi,'(a)')       'array 2'
      do i=1,jx
      write(ifhi,'(2e12.4)')   datx(i),datyu(i)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'
        else
      write(ifhi,'(a,2e11.3)')'openhisto xrange',0.,itmax*1.
      write(ifhi,'(a)')       'htyp '//htyp//' xmod lin ymod lin'
      write(ifhi,'(a)')       'array 2'
      itm=20
      do i=1,itm
      write(ifhi,'(2e12.4)')   (i-1.)*itmax/(itm-1.),yield
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'
        endif
c     ---
      return

           endif
           enddo

      end

c-----------------------------------------------------------------------
      subroutine xhnbte(iii)
c-----------------------------------------------------------------------
c fills histograms (iii>=0) or writes histogram to histo-file (iii<0)
c regarding exponential autocorrelation time and acceptance rate
c
c input:
c   requires complete run with application hadron (iappl=1)
c   or application metropolis (iappl=4)
c   ioceau=1 necessary
c
c  output:
c   for iii=0 (only valid for iappl=4):
c     data(nrevt): nrevt  (event number)               /cdat/
c     datb(nrevt): taui   (calculated corr time)       /cdat/
c     datc(nrevt): accrat (acceptance rate)            /cdat/
c     datd(nrevt): taue   (parametrized corr time)     /cdat/
c   for iii>0 (only valid for iappl=1):
c     nrclu=nrclu+1                                    /cnrclu/
c     data(nrclu): nrclu  (droplet number)             /cdat/
c     datb(nrclu): taui-taue (calc - param corr time)  /cdat/
c     datc(nrclu): accrat (acceptance rate)            /cdat/
c     datd(nrclu): avnp (average particle number)      /cdat/
c   for iii<0:
c     writes complete histogram (openhisto ... closehisto) to histofile
c       for iappl=4:                for iappl=1:
c         xpar1=1: (data,datb,datd) xpar1=1: (data,datb)
c         xpar1=2: (data,datc)      xpar1=2: (data,datd)
c                                   xpar1=3: (data,datc)
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxit=500000)
      common/count/nacc,nrej,naccit(maxit),nptot,npit(maxit)
      common/citer/iter,itermx
      common /clatt/nlattc,npmax
      common/cgctot/rmstot,ptltot
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      parameter (nbin=500)
      common/cdat/ data(nbin),datb(nbin),datc(nbin),datd(nbin)
      real dev(maxit)
      character cobs*5,cnc*5,cdz*5,czer*5
     *,cmom*5,cnp*7,cen*7,cvol*7,clatt*5,cit*5
      common/ctaue/taue

      if(ioceau.ne.1)stop'STOP: ioceau=1 required'
      if(iii.eq.0.and.iappl.ne.4)stop'STOP: iappl=4 required'
      if(iii.gt.0.and.iappl.ne.1)stop'STOP: iappl=1 required'

      if(iii.lt.0)jjj=nint(xpar1)

      id=0
      ish0=ish
c     ish=98

c          ----------------
           if(iii.ge.0)then
c          ----------------

      if(iii.gt.0)nrclu=nrclu+1
      if(nrclu.gt.500)return

c     mean
c     ----
      xnptot=nptot
      avnp=xnptot/(itermx-iternc)
      if(ish.ge.9)write(ifch,*)'event:',nrevt,'   droplet:',nrclu
     *,'   avnp:',avnp

c     calculate corfct_0
c     ------------------
      corzer=0.0
      do i=iternc+1,itermx
      dev(i)=npit(i)-avnp
      corzer=corzer+dev(i)**2
      enddo
      corzer=corzer/(itermx-iternc)
      if(ish.ge.9)write(ifch,*)'c_0:',corzer

c     calculate corfct_1
c     ------------------
      corone=0.0
      do i=iternc+1,itermx-1
      corone=corone+dev(i)*dev(i+1)
      enddo
      corone=corone/(itermx-iternc-1)

c     calculate initial autocorrelation time
c     -----------------------------------------
      if(corone.gt.1.e-30.and.corzer.gt.1.e-30)then
      r=alog(corone)-alog(corzer)
      if(ish.ge.9)write(ifch,*)'log rho_1:',r
      taui=(-1.)/r
      else
      taui=0.
      endif
      if(ish.ge.9)write(ifch,*)'tau_init:',taui

c     calculate parametrized autocorrelation time (if necessary)
c     ----------------------------------------------------------
      if(taue.eq.0.0)then
      e=tecm/volu
      b=1.1*(e+0.33)**0.66
      a=13.*(e+0.13)**(-0.65)
      tm=34.*(e+0.65)**(-0.61)
      t=a+b*volu
      taue=max(t,tm)
      endif

c     calculate acceptance rate
c     -------------------------
      xa=nacc
      ya=itermx
      accrat=xa/ya

c     write to data/b/c/d
c     -------------------
       if(iii.eq.0)then
      if(iozevt.gt.0)then
      data(nrevt)=iozero
      else
      data(nrevt)=nrevt
      endif
      datb(nrevt)=taui
      datc(nrevt)=accrat
      datd(nrevt)=taue
       else
      data(nrclu)=nrclu
      datb(nrclu)=taui-taue
      datc(nrclu)=accrat
      datd(nrclu)=avnp
       endif

c          -----------------------------------
           elseif(iii.lt.0.and.iappl.eq.4)then
c          -----------------------------------

      write(cmom,'(i3)')iomom
      write(cen,'(f7.3)')tecm
       if(ioobsv.eq.0)then
      write(cnp,'(f7.3)')ptltot
       else
       do i=1,nspecs
       if(ioobsv.eq.ispecs(i))id=i
       enddo
      write(cnp,'(f7.3)')ptlngc(id)
       endif
      write(cvol,'(f7.3)')volu
      write(clatt,'(i3)')nlattc
      write(cit,'(i5)')itermx
      if(ioobsv.eq.0)then
      write(cobs,'(a)')'all'
      else
      write(cobs,'(i5)')ioobsv
      endif
      write(cnc,'(i5)')iternc
      if(iozevt.eq.0)write(czer,'(i5)')iozero
      if(iozevt.gt.0)write(cdz,'(i5)')iozinc

      x1=1
      x2=nevent

      if(jjj.eq.1)then

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      if(iozevt.gt.0)then
      write(ifhi,'(a)')       'text 0 0 "xaxis iozero"'
      else
      write(ifhi,'(a)')       'text 0 0 "xaxis event"'
      endif
      write(ifhi,'(a)')       'text 0 0 "yaxis [t]?exp!"'
      write(ifhi,'(a)')       'text 0.05 0.95 "E='//cen//'"'
      write(ifhi,'(a)')       'text 0.2  0.95 "V='//cvol//'"'
      write(ifhi,'(a)')       'text 0.35 0.95 "N?g!='//cnp//'"'
      write(ifhi,'(a)')       'text 0.55 0.95 "observable  '//cobs//'"'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')       'array 2'
      do j=1,nevent
      write(ifhi,'(2e12.4)')data(j),datb(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0-'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')       'array 2'
      do j=1,nevent
      write(ifhi,'(2e12.4)')data(j),datd(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'

      elseif(jjj.eq.2)then

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      if(iozevt.gt.0)then
      write(ifhi,'(a)')       'text 0 0 "xaxis iozero"'
      else
      write(ifhi,'(a)')       'text 0 0 "xaxis event"'
      endif
      write(ifhi,'(a)')       'text 0 0 "yaxis acceptence rate"'
      write(ifhi,'(a)')       'text 0.05 0.95 "iomom= '//cmom//'"'
      write(ifhi,'(a)')       'text 0.2  0.95 "nlattc= '//clatt//'"'
      if(iozevt.eq.0)
     *write(ifhi,'(a)')       'text 0.35 0.95 "iozero= '//czer//'"'
      write(ifhi,'(a)')       'text 0.55 0.95 "itermx= '//cit//'"'
      write(ifhi,'(a)')       'text 0.75 0.95 "iternc= '//cnc//'"'
      if(iozevt.gt.0)
     *write(ifhi,'(a)')       'text 0.35  0.95 "dzero= '//cdz//'"'
      if(iorejz.eq.1)
     *write(ifhi,'(a)')    'text 0.25 0.05 "zeros rejected !"'
      if(ioinco.ge.1)then
      write(ifhi,'(a)')    'text 0.05 0.05 "hot start"'
      else
      write(ifhi,'(a)')    'text 0.05 0.05 "cold start"'
      endif
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')       'array 2'
      do j=1,nevent
      write(ifhi,'(2e12.4)')data(j),datc(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'

      endif

c          -----------------------------------
           elseif(iii.lt.0.and.iappl.eq.1)then
c          -----------------------------------

      if(ioobsv.eq.0)then
      write(cobs,'(a)')'all'
      else
      write(cobs,'(i5)')ioobsv
      endif

      x1=1
      x2=nrclu

      if(jjj.eq.1)then

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')       'text 0 0 "xaxis droplet"'
      write(ifhi,'(a)')       'text 0 0 "yaxis [D][t]?exp!"'
      write(ifhi,'(a)')       'text 0.05 0.91 "[D][t]?exp!=[t]?measured!
     *-[t]?parametrized"'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a,a,a)')'yrange',' auto',' auto'
      write(ifhi,'(a)')       'array 2'
      do j=1,nrclu
      write(ifhi,'(2e12.4)')data(j),datb(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'

      elseif(jjj.eq.2)then

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')       'text 0 0 "xaxis droplet"'
      write(ifhi,'(a)')       'text 0 0 "yaxis N?obs!"'
      write(ifhi,'(a)')       'text 0.05 0.95 "observable  '//cobs//'"'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')       'array 2'
      do j=1,nrclu
      write(ifhi,'(2e12.4)')data(j),datd(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'

      elseif(jjj.eq.3)then

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')       'text 0 0 "xaxis droplet"'
      write(ifhi,'(a)')       'text 0 0 "yaxis accep. rate"'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')       'array 2'
      do j=1,nrclu
      write(ifhi,'(2e12.4)')data(j),datc(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'

      endif

c          -----
           endif
c          -----

      ish=ish0
      return
      end

c-------------------------------------------------------------------------
      subroutine xhnbti(iii)
c-------------------------------------------------------------------------
c fills histograms (iii=0) or writes histogram to histo-file (iii<0)
c regarding integrated autocorrelation time and corresponding multiplicity
c and variance
c
c input:
c   requires complete run with application metropolis (iappl=4)
c   iociau=1 necessary
c   iompar (parameter for windowing algorithm by  a.d.sokal) must
c   be set to 3 < c_M < 11
c
c  output:
c   for iii=0 (only valid for iappl=4):
c     data(nrevt): nrevt (event number)              /cdat/
c     datb(nrevt): tau   (calculated int corr time)  /cdat/
c     datc(nrevt): stau  (variance tau)              /cdat/
c     datd(nrevt): avnp  (multiplicity)              /cdat/
c     date(nrevt): sobs  (variance multiplicity)     /cdat/
c     datf(nrevt):       (gc multiplicity)           /cdat/
c   for iii=0 and iosngl>0:
c     writes complete set of histograms (newpage zone 1 3 1
c     openhisto ... closehisto plot0 ... openhisto ... closehisto plot 0)
c     concerning acceptance rate, rejection rate, correlation function
c     for specific event, specified by value of iosngl (=nrevt+1)
c   for iii<0:
c     writes complete histogram (openhisto ... closehisto) to histofile
c       xpar1=1: (data,datb,datc)
c       xpar1=2: (data,datd,date,datf)
c------------------------------------------------------------------------
      include 'epos.inc'
      parameter(maxit=500000)
      common/count/nacc,nrej,naccit(maxit),nptot,npit(maxit)
      common/citer/iter,itermx
      parameter(maxp=6000)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      common /clatt/nlattc,npmax
      common/cgctot/rmstot,ptltot
      parameter (mspecs=400)
      common/cspecs/nspecs
     .,ispecs(mspecs+1),aspecs(mspecs+1),gspecs(mspecs+1)
      common/cgchg/rmsngc(mspecs),ptlngc(mspecs),chemgc(mspecs),tem
      parameter (nbin=500)
      common/cdat2/data(nbin),datb(nbin),datc(nbin),datd(nbin)
     *,date(nbin),datf(nbin),datg(nbin),dath(nbin)
      common/cdat3/datx(nbin),daty(nbin),datz(nbin),datr(nbin)
     *,dats(nbin)
      real corfct(maxit),dev(maxit)
      character cobs*5,cdz*5,ccuev*5,cmpar*3,ctau*7
      character cmom*5,cnp*7,cen*7,cvol*7,clatt*5,cit*5,cavnp*7
      character cnacc*10,cnrej*10,caver*10,cioz*5,ciom*3,cnlat*5

      if(iociau.ne.1)stop'STOP: iociau=1 required'
      if(iii.eq.0.and.iappl.ne.4)stop'STOP: iappl=4 required'
      if(iii.gt.0)stop'STOP: iii>0 not supported'

      jjj=nint(xpar1)
      id=0

c          ----------------
           if(iii.eq.0)then
c          ----------------

c     mean
c     ----
      xnptot=nptot
      avnp=xnptot/(itermx-iternc)
      if(ish.ge.9)write(ifch,*)'event:',nrevt,'   avnp:',avnp

c     normalization of corfct_i
c     -------------------------
      corzer=0.0
      do i=iternc+1,itermx
      dev(i)=npit(i)-avnp
      if(ish.ge.9)write(ifch,*)'i:',i,'  dev_i:',dev(i)
      corzer=corzer+dev(i)**2
      enddo
      corzer=corzer/(itermx-iternc)
      if(ish.ge.9)write(ifch,*)'c_0:',corzer

c     calculate corfct_i
c     ------------------
      nt=itermx-iternc-1
      do it=1,nt
      corfct(it)=0.0
      do i=iternc+1,itermx-it
      corfct(it)=corfct(it)+dev(i)*dev(i+it)
      enddo
      corfct(it)=corfct(it)/(itermx-iternc-it)
      if(it.le.10.and.ish.ge.9)
     *write(ifch,*)'t:',it,'  c_t:',corfct(it)
      enddo

c     calculate initial autocorrelation time
c     -----------------------------------------
      if(corfct(1).gt.1.e-30.and.corzer.gt.1.e-30)then
      r=alog(corfct(1))-alog(corzer)
      if(ish.ge.9)write(ifch,*)'log rho_1:',r
      taui=(-1.)/r
      else
      taui=0.
      endif
      if(ish.ge.9)write(ifch,*)'tau_init:',taui

c     calculate integrated autocorrelation time
c     -----------------------------------------
      k=1
      mpar=iompar
      tau=taui
      taux=taui
      taum=0.0
      mcut=0
      if(ish.ge.9)write(ifch,*)'initial tau:',tau,'   c_M:',mpar

        if(corzer.gt.1.e-30)then

5     mcut=mpar*abs(taux)
      tauo=tau
      tau=.5
      do it=1,mcut
      tau=tau+corfct(it)/corzer
      enddo
      taum=taum+tau
      taux=taum/k
      if(ish.ge.9)write(ifch,*)'iteration:',k,'   M:',mcut,'  tau:',tau
      if(mcut.lt.(mpar*tau).or.mcut.gt.(10.*tau))then
      dt=abs(tau-tauo)
      if(k.lt.20.and.dt.gt.0.2)then
      k=k+1
      goto5
      endif
      endif
      mcut=mpar*abs(taux)
      if(ish.ge.9)write(ifch,*)'tau_mean:',taux,'   M:',mcut
      tau=0.5
      do it=1,mcut
      tau=tau+corfct(it)/corzer
      enddo

       endif

      vtau=(2.*mcut+1.)*2./(itermx-iternc)*tau**2
      stau=0.0
      if(vtau.ge.0.0)stau=sqrt(vtau)
      if(ish.ge.9)
     *write(ifch,*)'tau_int:',tau,'   var:',vtau,'   sig:',stau

c     calculate variance of observable
c     --------------------------------
      vobs=2.*tau*corzer/(itermx-iternc)
      sobs=0.0
      if(vobs.ge.0.0)sobs=sqrt(vobs)

c     write to data-f
c     ---------------
       if(ioobsv.eq.0)then
      datf(nrevt)=ptltot
       else
      do j=1,np
      if(ioobsv.eq.ispecs(j))id=j
      enddo
      datf(nrevt)=ptlngc(id)
       endif
      datb(nrevt)=tau
      datc(nrevt)=stau
      date(nrevt)=sobs
      datd(nrevt)=avnp
      if(iozevt.gt.0)then
      data(nrevt)=iozero
      else
      data(nrevt)=nrevt
      endif

c          -------------------------
           if(iosngl.eq.nrevt+1)then
c          -------------------------

      nb=itermx/iterpl
      if(nb.gt.nbin)nb=nbin

      datx(1)=iterpl/2
      daty(1)=naccit(1)
      datz(1)=1-naccit(1)
      if(iterpl.ge.2)then
      do j=1,iterpl-1
      daty(1)=daty(1)+naccit(1+j)
      datz(1)=datz(1)+1-naccit(1+j)
      enddo
      endif
      datr(1)=daty(1)/iterpl
      dats(1)=datz(1)/iterpl
      do i=2,nb
      datx(i)=datx(i-1)+iterpl
      daty(i)=daty(i-1)
      datz(i)=datz(i-1)
      do j=1,iterpl
      daty(i)=daty(i)+naccit((i-1)*iterpl+j)
      datz(i)=datz(i)+1-naccit((i-1)*iterpl+j)
      enddo
      datr(i)=daty(i)/i/iterpl
      dats(i)=datz(i)/i/iterpl
      enddo
      b=nacc
      c=itermx
      avrate=b/c
      write(cnacc,'(i6)')nacc
      write(cnrej,'(i6)')nrej
      write(caver,'(f5.3)')avrate
      write(cioz,'(i5)')iozero
      write(ciom,'(i3)')iomom
      write(cnlat,'(i5)')nlattc
      x1=datx(1)
      x2=datx(nb)

      write(ifhi,'(a)')       'newpage zone 1 3 1 openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')       'text 0 0 "xaxis iterations"'
      write(ifhi,'(a)')       'text 0 0 "yaxis acceptence rate"'
      write(ifhi,'(a)')       'text 0.6 0.5 "accepted '//cnacc//'"'
      write(ifhi,'(a)')       'text 0.6 0.4 "rejected  '//cnrej//'"'
      write(ifhi,'(a)')       'text 0.6 0.3 "aver. rate  '//caver//'"'
      write(ifhi,'(a)')       'text 0.4 0.5 "nlattc='//cnlat//'"'
      write(ifhi,'(a)')       'text 0.4 0.4 "iozero='//cioz//'"'
      write(ifhi,'(a)')       'text 0.4 0.3 "iomom='//ciom//'"'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')       'array 2'
      do j=1,nb
      write(ifhi,'(2e12.4)')datx(j),datr(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0-'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')       'array 2'
      do j=1,nb
      write(ifhi,'(2e12.4)')datx(j),dats(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0'

      m=min(mcut,500)
      do i=1,m
      datg(i)=i
      dath(i)=1000.
      if(corzer.gt.1.e-30)dath(i)=corfct(i)/corzer
      enddo
      write(ccuev,'(i5)')nrevt+1
      write(cmpar,'(i3)')mpar
      write(ctau,'(i7)')tau
      x1=1.
      x2=m

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')       'text 0 0 "xaxis t"'
      write(ifhi,'(a)')       'text 0 0 "yaxis correl. func."'
      write(ifhi,'(a)')       'text 0.8 0.95 "event '//ccuev//'"'
      write(ifhi,'(a)')'text 0.05 0.95  "window parameter= '//cmpar//'"'
      write(ifhi,'(a)')       'text 0.35 0.95  "tau= '//ctau//'"'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a,a,a)')'yrange',' auto',' auto'
      write(ifhi,'(a)')       'array 2'
      do j=1,m
      write(ifhi,'(2e12.4)')datg(j),dath(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0'

c          -----
           endif
c          -----

c          --------------------
           elseif(iii.lt.0)then
c          --------------------

      write(cmom,'(i3)')iomom
       if(ioobsv.eq.0)then
      write(cnp,'(f7.3)')ptltot
       else
      do j=1,np
      if(ioobsv.eq.ispecs(j))id=j
      enddo
      write(cnp,'(f7.3)')ptlngc(id)
       endif
      write(cen,'(f7.3)')tecm
      write(cvol,'(f7.3)')volu
      write(clatt,'(i3)')nlattc
      write(cit,'(i5)')itermx
      write(cavnp,'(f7.3)')avnp
      if(iozevt.gt.0)
     *write(cdz,'(i5)')iozinc
      write(cmpar,'(i3)')mpar
      if(ioobsv.eq.0)then
      write(cobs,'(a)')'all'
      else
      write(cobs,'(i5)')ioobsv
      endif

      x1=data(1)
      x2=data(nevent)

      if(jjj.eq.1)then

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp pnt xmod lin ymod lin'
      if(iozevt.gt.0)then
      write(ifhi,'(a)')       'text 0 0 "xaxis iozero"'
      else
      write(ifhi,'(a)')       'text 0 0 "xaxis event"'
      endif
      write(ifhi,'(a)')       'text 0 0 "yaxis [t]?int!"'
      write(ifhi,'(a)')'text 0.05 0.95  "window parameter '//cmpar//'"'
      if(iozevt.gt.0)
     *write(ifhi,'(a)')       'text 0.8  0.95 "dzero= '//cdz//'"'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')       'array 3'
      do j=1,nevent
      write(ifhi,'(3e12.4)')data(j),datb(j),datc(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'

      elseif(jjj.eq.2)then

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp pnt xmod lin ymod lin'
      if(iozevt.gt.0)then
      write(ifhi,'(a)')       'text 0 0 "xaxis iozero"'
      else
      write(ifhi,'(a)')       'text 0 0 "xaxis event"'
      endif
      write(ifhi,'(a)')       'text 0 0 "yaxis multiplicity"'
      write(ifhi,'(a)')       'text 0.05 0.95 "E='//cen//'"'
      write(ifhi,'(a)')       'text 0.2 0.95 "V='//cvol//'"'
      write(ifhi,'(a)')       'text 0.35 0.95 "N?g!='//cnp//'"'
      write(ifhi,'(a)')       'text 0.55 0.95 "observable  '//cobs//'"'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a,a,a)')'yrange',' auto',' auto'
      write(ifhi,'(a)')       'array 3'
      do j=1,nevent
      write(ifhi,'(3e12.4)')data(j),datd(j),date(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto   plot 0-'


      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lda xmod lin ymod lin'
      write(ifhi,'(a)')       'array 2'
      do j=1,nevent
      write(ifhi,'(2e12.4)')data(j),datf(j)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto'

      endif

c          -----
           endif
c          -----

      return
      end
      
