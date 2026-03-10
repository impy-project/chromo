C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
      subroutine decayall(n,i99)
c-----------------------------------------------------------------------
c  decay of objects n to nptl, including their children
c-----------------------------------------------------------------------
      include 'epos.inc'
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cnptlbur/nptlbur
      common/cij99/ij99
      character*24 txt
      logical go
      !-------------------------------
      !call ranfini(222222222,2,1) !to have same random nr for 
      !call ranfini(111111111,1,2) !  mx and mt run
      !print*,'TEST core decay ',rangen()
      !do i=1,nptl
      !if(ityptl(i)/10.eq.6)istptl(i)=7
      !enddo
      !-------------------------------
      call utpri('decaya',ish,ishini,5)
      ij99=i99
      if(i99.eq.99)then
        ndecay_save=ndecay
        idecay_save=idecay
        ndecay=0
        idecay=ihdecay
        istdec=0
      elseif(i99.eq.999)then
        ndecay_save=ndecay
        idecay_save=idecay
        ndecay=111111110
        idecay=1
        istdec=0
      elseif(i99.eq.3)then
        ndecay_save=ndecay
        idecay_save=idecay
        if(nptl.gt.mxptl-(nptl-nbdky))goto 999 !skip copy of particles used for analysis only
        ndecay=0
        idecay=ihdecay
        do ip=n,nptlbur
          if(istptl(ip).eq.3.or.istptl(ip).eq.4)then
            nptl=nptl+1
            call checkcccptl(nptl)
            if(nptl.le.mxptl)then
              call utrepl(nptl,ip)
              istptl(nptl)=8
              iorptl(nptl)=ip
            else
              call utstop('decayall: mxptl too small&')
            endif
          endif
        enddo
        istdec=8
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! nptlfff=149500
      ! nptlxxx=nptl
      ! write(ifmt,*)'#####  CHECK ##### nptl artificially set to',nptlfff
      ! do ij=nptl+1,nptlfff
      ! istptl(ij)=999
      ! enddo
      ! nptl=nptlfff
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      else
        istdec=0
        ndecay_save=ndecay
        idecay_save=idecay      
      endif
      nptlbd=nptl
      if(ndecay.eq.1)goto 999
      if(idecay.eq.0)goto 999
      if(iappl.eq.4.or.iappl.eq.7.or.iappl.eq.9)then
        nptli=1
      else
        nptli=maproj+matarg+1
      endif
      ttaus=1
      np1=n
 1    np2=nptl
      ip=np1-1
      nnnpt=0
      do while (ip.lt.np2)
        ip=ip+1
         call getistptl(ip,ist)
         call getidptl(ip,id)
         if(ist.eq.istdec)then  !consider last generation particles
           if(id.ne.0)then !from spherio: holes in the list (id=0)
             nptlb=nptl
             ntry=0
 2           ntry=ntry+1
             if(i99.eq.3)call setistptl(ip,0) 
             call hdecas(ip,iret)
             if(i99.eq.3)then
               call getistptl(ip,ist)
               if(ist.eq.0)ival=8
               if(ist.eq.1)ival=6
               call setistptl(ip,ival)
             endif  
             if(iret.eq.1)then
               call getidptl(ip,idip)
               call getistptl(ip,istip)
               call getityptl(ip,ityip)
               call getpptl(ip,p1,p2,p3,p4,p5) 
               if(ish.ge.2)then
                  idor=0
                  call getiorptl(ip,iorip)
                  if(iorip.gt.0)then
                    call getidptl(iorip,idor)
                  endif
                  write(ifch,*)('*',k=1,60)
                  write(ifch,*)'***** hdecas_iret=1'
                  write(ifch,*)'***** i=',ip,'    id=',idip
     .            ,'    ist=',istip,'    ity=',ityip
     .            ,'    idor=',idor
                  write(ifch,*)'***** mass=',p5
               endif
               jerr(9)=jerr(9)+1
               if(ntry.eq.2)then
                 write(ifmt,*)'####### ERROR 21072016 ; id ist ity = '
     .           ,idip,istip,ityip
                 stop
               endif 
               call idmass(idip,am)
               p4=sqrt(p1*p2+p2*p2+p3*p3+am*am)  !put particle on-shell to avoid problems later
               call setpptl(ip,p1,p2,p3,p4,am)
               if(ish.ge.2)then
                  write(ifch,*)'***** changed to ',am
                  write(ifch,*)('*',k=1,60)
               endif
               nptl=nptlb
               goto 2
             endif
           endif
         endif
c remove useless particles if not enough space
         if(nclean.gt.0.and.nptl.gt.mxptl/4.and.nnnpt.eq.0)then
           nnnip=0
           do iii=maproj+matarg+1,np2
             go=.true.
             if(nclean.eq.1.and.istptl(iii).le.istmax)go=.false.
             if(nclean.eq.2.and.istptl(iii).lt.10)go=.false.
             if(go.and.mod(istptl(iii),10).ne.0)then
               istptl(iii)=99
               nnnpt=nnnpt+1
               if(iii.le.ip)nnnip=nnnip-1
             endif
           enddo
           if(nnnpt.gt.0)then
             nptl0=nptl
             call utclea(1,nptl0)
             np2=np2-nnnpt
             ip=ip-nnnip
             nptli=nptl
             nbdky=nbdky-nnnpt
             if(ish.ge.1)
     .       write(ifch,*)"clean decayall ->",nptl,nnnpt,ip,np2
           else
             if(ish.ge.1)
     .       write(ifch,*)"no part. for clean decayall !"
             nnnpt=-1
           endif
         endif
      enddo
      if(i99.eq.3)then
        do ip=np2+1,nptl
          call setistptl(ip,8)
        enddo
      endif
      np1=np2+1
      if(np1.le.nptl)goto1
      if(ish.ge.1.and.i99.eq.3)
     . write(ifch,'(a,i6)')'final number of slots used:',nptl
      if(nptl.eq.nptlbd)goto 999
      if(ish.ge.2)then
      txt=              ' &     '
      if(i99.eq.99)txt= ' (99)& '
      if(i99.eq.999)txt=' (999)& '
      if(i99.eq.3)txt=' (3)& '
      call alist('list after decayall'//txt,1,nptl)
      endif
      call clop(3)
  999 continue
      if(i99.eq.99.or.i99.eq.999.or.i99.eq.3)then
        ndecay=ndecay_save
        idecay=idecay_save
      endif
      if(ish.ge.4)write(ifch,'(a,i8,i5)')
     .'end of decayall - n i99 = ',n,i99
      call utprix('decaya',ish,ishini,5)
      end

c-----------------------------------------------------------------------
      subroutine hdecas(i,iret)
c-----------------------------------------------------------------------
c  decay of object i  (main decay routine)
c-----------------------------------------------------------------------

      include 'epos.inc'
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat!,zor,tor
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      double precision ttaux,ttauz
      !integer jcdu(nflav,2)

      !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
      !  write(ifch,'(a,3i8,i4,5f10.2,4x,4f10.2)')
      !. ' WWWWWW     hacas mother ',iorptl(i),i,idptl(i),ityptl(i)
      !. ,sqrt(pptl(1,i)**2+pptl(2,i)**2)
      !. , xorptl(4,i),xorptl(3,i),xorptl(1,i),xorptl(2,i)
      !. ,tivptl(2,i)
      !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

      iret=0
      nptlb=nptl

c no last generation -> no decay

      if(istptl(i).ne.0)return

      if(nptl.gt.mxptl-10)then
        if(model.eq.12)then
          write(ifmt,*)
     .    'Number of particles too large after decay with DPMJET !'
          nptl=mxptl-10
        else
          call alist('end&',1,nptl)
          call utstop('hdecas: mxptl too small&')
        endif
      endif
c entry

      call utpri('hdecas',ish,ishini,5)
      ttauz=ttaus


c skip nuclei
      if(idptl(i).gt.1000000000)return

c small droplet decay

      if(iabs(idptl(i)).gt.100000000)then
        print *,'hdecas',idptl(i)
        stop'hdecas: no longer supported (2).       '
      endif

c  ordinary decay


c      if(model.eq.1)then
c        call idmass(111,amrho0)
c        call idmass(221,amomeg)
c        ioi=iorptl(i)
c        if(ioi.gt.0.and.(idptl(i).eq.111.or.idptl(i).eq.221))then
c          if(.not.(iabs(idptl(ioi)).lt.10000
c     *         .and.jorptl(i).eq.0))then
c            
c            if(iLHC.eq.1.and.((ityptl(i).ge.20.and.ityptl(i).le.39)
c     *           .or.ityptl(i).eq.42.or.ityptl(i).eq.52
c     *           .or.ityptl(i).eq.53
c     *           .or.(ityptl(i).eq.43.and..not.(abs(idproj).lt.200
c     *                                 .and.pptl(5,ioi).lt.2.75))))then
cc mix rho and omegas only from string production (with remnant string mass high enough not to be low mass diffraction (fine tuning from NA22) and if not decay product
c              if(idptl(i).eq.111)idptl(i)=221
c              if(idptl(i).eq.221.and.rangen().gt.0.5)idptl(i)=111
c            elseif(iLHC.eq.0.and..not.(ityptl(i).eq.60))then
c              if(idptl(i).eq.111)idptl(i)=221
c              if(idptl(i).eq.221.and.rangen().gt.0.5)idptl(i)=111
c            endif
c            
c          endif
c        endif
c      endif

      if(ctaumin.gt.0.)then
        call idtau(idptl(i),1.,1.,ctau)       !ctau in fm
        if(ctau*1.e-13.gt.ctaumin)goto 1000   !ctaumin in cm
      endif

      ida=iabs(idptl(i))

      if(.not.(iappl.eq.7.and.i.eq.1))then
      if(mod(ndecay        ,10).eq.1
     *.and.ida.ne.0.and.ida.lt.10000)goto1000
      if(mod(ndecay/10     ,10).eq.1.and.ida.eq.  20)goto1000
      if(mod(ndecay/100    ,10).eq.1.and.ida.eq.2130)goto1000
      if(mod(ndecay/1000   ,10).eq.1.and.ida.eq.1130)goto1000
      if(mod(ndecay/1000   ,10).eq.1.and.ida.eq.2230)goto1000
      if(mod(ndecay/10000  ,10).eq.1.and.ida.eq.2330)goto1000
      if(mod(ndecay/10000  ,10).eq.1.and.ida.eq.1330)goto1000
      if(mod(ndecay/100000 ,10).eq.1.and.ida.eq.3331)goto1000
      if(mod(ndecay/1000000,10).eq.1.and.ida.eq. 110)goto1000
      if(mod(ndecay/1000000,10).eq.1.and.ida.eq. 220)goto1000

      if(nrnody.gt.0)then
      do nod=1,nrnody
      if(idptl(i).eq.nody(nod))goto 1000
      enddo
      endif


      endif

      call hdecay(i,iret)
      if(iret.eq.1)goto1000
      if(nptl.le.nptlb)then
        iret=-1
        goto 1000
      endif

c ---successful decay---

      istptl(i)=1
      ifrptl(1,i)=nptlb+1
      ifrptl(2,i)=nptl

      t=tivptl(2,i)
      x=xorptl(1,i)+(t-xorptl(4,i))*pptl(1,i)/pptl(4,i)
      y=xorptl(2,i)+(t-xorptl(4,i))*pptl(2,i)/pptl(4,i)
      z=xorptl(3,i)+(t-xorptl(4,i))*pptl(3,i)/pptl(4,i)
      call jtaux(t,z,ttaux)
      ttaus=ttaux
      if( ttaus.gt.0d0 ) then
        call jtauin
        call jtaus(z,ttest,sz)
        if (abs(t-ttest).gt.1e-5*ttest.and.ish.ge.1) then
          call utmsg('hdecas&')
          write(ifch,*)'*****  t /= ttest'
          write(ifch,*)t,ttest,i,z,t,xorptl(3,i),xorptl(4,i)
     $         ,pptl(3,i),pptl(4,i)
          call utmsgf
        endif
      endif

c loop over decay products

      do 20 n=nptlb+1,nptl
      iorptl(n)=i
      jorptl(n)=0
      istptl(n)=0
      ifrptl(1,n)=0
      ifrptl(2,n)=0
      rad=0
      phi=0
      ti=t
      zi=z
      xorptl(1,n)=x + rad*cos(phi)
      xorptl(2,n)=y + rad*sin(phi)
      xorptl(3,n)=zi
      xorptl(4,n)=ti
c      io=n
c1     io=iorptl(io)
c      if(ish.ge.4)write(ifch,*)'io = ',io,'  origin: ',iorptl(io)
c      if(io.eq.iorptl(io))call utmsg("Strange iorptl in hdecas&")
cc security to avoid infinite loop
c      if(iorptl(io).gt.0.and.io.ne.iorptl(io))goto 1  
c      if(ish.ge.4)write(ifch,*)'origin: ',io,idptl(io)
c      zor=xorptl(3,io)
c      tor=xorptl(4,io)
c      call idquac(io,nq,ndummy1,ndummy2,jcdu)
      tivptl(1,n)=ti
c      r=rangen()
c      tauran=-taurea*alog(r)
c      call jtaix(n,tauran,zor,tor,zis,tis)
c      tivptl(1,n)=amax1(ti,tis)
      call idtau(idptl(n),pptl(4,n),pptl(5,n),taugm)
      ityptl(n)=ityptl(i)
      r=rangen()
      tivptl(2,n)=ti+taugm*(-alog(r))
      radptl(n)=0.
      dezptl(n)=0.
      itsptl(n)=itsptl(i)
      rinptl(n)=rinptl(i)
      zpaptl(1,n)=zpaptl(1,i)
      zpaptl(2,n)=zpaptl(2,i)
      !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
      !  write(ifch,'(a,3i8,i4,5f10.2,4x,4(f10.2,1x))')
      !. ' WWWWWW     hdecas       ',iorptl(n),n,idptl(n),ityptl(n)
      !. ,sqrt(pptl(1,n)**2+pptl(2,n)**2)
      !. , xorptl(4,n),xorptl(3,n),xorptl(1,n),xorptl(2,n)
      !. ,tivptl(1,n),tivptl(2,n)
      !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
20    continue

      if(iabs(idptl(nptlb+1)).le.6) then
        call gakli2(0,0)
        write (*,*) 'nptlb+1,nptl:',nptlb+1,nptl
        istptl(nptlb+1)=1
        do n=nptlb+2,nptl
          istptl(n)=20
        enddo
        call gakfra(0,iret)
        if(iret.eq.1)goto1000
        call gakli2(0,0)
      endif

1000  continue
      ttaus=ttauz
      call jtauin
      call utprix('hdecas',ish,ishini,5)
      return
      end

c-----------------------------------------------------------------------
      subroutine hdecay(ip,iret)
c-----------------------------------------------------------------------
c  decays particle ip from /cptl/
c  for ip being no resonance: call StaHad
c  for ip being resonance: standard resonance decay  procedure
c-----------------------------------------------------------------------
      include 'epos.inc'
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/wco/wmass2,wgam2
      parameter (mxlook=10000,mxdky=4000)
      common/dkytab/look(mxlook),cbr(mxdky),mode(5,mxdky)
      dimension pgen(5,10),rnd(10),u(3),beta(3)
     1     ,reduce(10)
      dimension prest(4,10),kno(10)
      data reduce/1.,1.,2.,5.,15.,60.,250.,1500.,1.2E4,1.2E5/
      data twome/1.022006e-3/

c          fctn definitions
      dot(i1,i2)=prest(4,i1)*prest(4,i2)-prest(1,i1)*prest(1,i2)
     *-prest(2,i1)*prest(2,i2)-prest(3,i1)*prest(3,i2)
c          charged w propagator.
      wprop(z)=(z-wmass2**2)**2+(wmass2*wgam2)**2

      call utpri('hdecay',ish,ishini,5)

      ipp=ip
      iret=0
      nptlb=nptl

      if(ish.ge.4.or.idptl(ip).eq.0)write(ifch,*)'ip,id,mass: ',ip
     &                                ,idptl(ip),pptl(5,ip),ibreit
      if(idptl(ip).eq.0)call utstop("Pb in hdecay with id=0!&")


      if(model.eq.4.and.iappl.eq.7)then
        if(abs(idptl(ipp)).gt.13.and.abs(idptl(ipp)).ne.1120
     &.and.abs(idptl(ipp)).ne.15)call decaymod(ipp,iret)
        if(iret.gt.0)goto 1000
        naddptl=0
        goto 900 
      endif

c     no k_long decay
c     ---------------
c     if(idptl(ip).eq.-20)goto1000

c     select decay mode
c     -----------------
      ntry=0
      wimax=-1.
      amss0=0.
2     ntry=ntry+1
           if(ntry.gt.100)then
      if(ish.ge.1)then
      call utmsg('hdecay&')
      write(ifch,*)'*****  decay not possible. iret = 1.'
      call utmsgf
      endif
      iret=1
      goto1000
           endif
      idlv1=idptl(ip)
      amss=pptl(5,ip)

c Decay of deuteron

      if(abs(idlv1).eq.17)then
        amss=1.01*amss
        naddptl=2
        call idmass(1120,amnew)
        pptl(5,nptl+1)=amnew
        idptl(nptl+1)=1120
        sum=amnew
        call idmass(1220,amnew)
        pptl(5,nptl+2)=amnew
        idptl(nptl+2)=1220
        sum=sum+amnew
        goto 111
      endif

c Decay of triton

      if(abs(idlv1).eq.18)then
        amss=1.01*amss
        naddptl=3
        call idmass(1120,amnew)
        pptl(5,nptl+1)=amnew
        idptl(nptl+1)=1120
        sum=amnew
        call idmass(1220,amnew)
        pptl(5,nptl+2)=amnew
        idptl(nptl+2)=1220
        sum=sum+amnew
        call idmass(1220,amnew)
        pptl(5,nptl+3)=amnew
        idptl(nptl+3)=1220
        sum=sum+amnew
         goto 111
      endif

c Decay of alpha

      if(abs(idlv1).eq.19)then
        amss=1.01*amss
        naddptl=4
        call idmass(1120,amnew)
        pptl(5,nptl+1)=amnew
        idptl(nptl+1)=1120
        sum=amnew
        call idmass(1220,amnew)
        pptl(5,nptl+2)=amnew
        idptl(nptl+2)=1220
        sum=sum+amnew
        call idmass(1120,amnew)
        pptl(5,nptl+3)=amnew
        idptl(nptl+3)=1120
        sum=sum+amnew
        call idmass(1220,amnew)
        pptl(5,nptl+4)=amnew
        idptl(nptl+4)=1220
        sum=sum+amnew
        goto 111
      endif

      if(ibreit.eq.1)then
c  Introduce mass broadening due to life time (Breit-Wigner distribution)
      ctau=tivptl(2,ip)-tivptl(1,ip)     !ctau*gamma (fm)
c      ctau=ctau*pptl(5,ip)/pptl(4,ip)            !ctau (fm)
      if(ctau.gt.0..and.ctau.lt.1.e16)then
        call idtau(idlv1,1.,1.,ctau0)       !nominal ctau in fm
        wi=0.197/ctau0                      !width (GeV)
        call idmass(idlv1,amss0)
        if(abs(idlv1).eq.20.or.idlv1.eq.111.or.idlv1.eq.221)
     &  amss=amss0    !because of particle mixing
c only if mean mass is used
        if(abs(amss-amss0).lt.1.e-5.or.wimax.gt.0.)then
          if(wimax.lt.0.)then
            atwimax=atan(20.)       !use large limit (10 sigma) refined later
          else
            atwimax=atan(2.*wimax/wi)
          endif
          amss=amss0+0.5*wi*tan((2.*rangen()-1.)
     &                          *atwimax) !maximum broadening
c      write(*,*)'broadening:',ip,idlv1,abs(pptl(5,ip)-amss)/wi,amss,ntry
        endif
      endif
      endif

c  select one of the decay channel
      ipoint=look(iabs(idlv1))-1
      if(idlv1.eq.-20)ipoint=look(320)-1
      if(ipoint.lt.0) goto1000
      try=rangen()
      !-------------------------------
      idforce=0
      if(ifoele.eq.1)then !charmed, bottomed -> only electrons
        jd=abs(idlv1)
        if(   jd.eq.240.or.jd.eq.140.or.jd.eq.340.or.jd.eq.2140
     .    .or.jd.eq.150.or.jd.eq.250.or.jd.eq.350.or.jd.eq.2150 ) then
          idforce=12
        endif
      endif
      if(ifoele.eq.2)then !charmed, bottomed -> only muons
        jd=abs(idlv1)
        if(   jd.eq.240.or.jd.eq.140.or.jd.eq.340.or.jd.eq.2140
     .    .or.jd.eq.150.or.jd.eq.250.or.jd.eq.350.or.jd.eq.2150 ) then
          idforce=14
        endif
      endif
      if(ifoele.eq.3)then !JPsi -> only muons
        jd=abs(idlv1)
        if(   jd.eq.441 ) then
          idforce=14
        endif
      endif
      !-----------------------------------
      if(idforce.gt.0)then
        sumbranch=0
        cbra=0.
        do kfo=1,2 ! 1 -> sumbranch, 2 -> ipoint  
          ipt=ipoint+1
          sumbr=0
          igo=1
          do while((ipt.eq.ipoint+1.or.cbra.lt.0.9999).and.igo.eq.1)
            iok=0
            do i=1,5
              if(abs(mode(i,ipt)).eq.idforce)iok=1
            enddo
            if(iok.eq.1)then
              if(ipt.eq.ipoint+1)then
                delbr=cbr(ipt)
              else
                delbr=cbr(ipt)-cbr(ipt-1)
              endif
              sumbr=sumbr+delbr
              !print*,'CHECK dky',kfo,delbr,sumbr,try*sumbranch,igo
              if(kfo.eq.2.and.try*sumbranch.le.sumbr)then
                igo=0 !last cycle
                ipoint=ipt
              endif 
            endif
            cbra=cbr(ipt)
            ipt=ipt+1
          enddo!while
          sumbranch=sumbr
        enddo!kfo
        !print*,'CHECK dky',ipoint,idlv1,(mode(i,ipoint),i=1,5)
      else
  100   ipoint=ipoint+1
        if(ish.ge.4)write(ifch,*)'ipoint,cbr,try',ipoint,cbr(ipoint),try
        if(try.gt.cbr(ipoint)) goto100
      endif 
      !-----------------------------------
      if(idforce.gt.0)then 
        iele=0
        do i=1,5
          if(abs(mode(i,ipoint)).eq.idforce)iele=1
        enddo
        if(iele.eq.0)then
          write(ifmt,*)'ERROR 24042016 ',idlv1,try,ipoint
          stop
        endif
      endif
      !----------------------------------
      naddptl=0
      sum=0.
c      nstart=nptl+1                  !?????????????????unused
      new=0
      do 110 i=1,5         !store id and mass of products
        if(mode(i,ipoint).eq.0) goto 110
        if(nptl+naddptl+1.gt.mxptl) goto 9999
        if(iabs( mode(1,ipoint)) .le. 6.and.i.eq.2)then   !decay into quark ???
          call vedi(mode(1,ipoint),mode(2,ipoint),k3,idlv1)
          idptl(new)=idlv1
          call idmass(idlv1,amnew)
          pptl(5,new)=amnew
          sum=pptl(5,new)
        else                                 !decay into particles
          naddptl=naddptl+1
          new=nptl+naddptl
          idptl(new)=mode(i,ipoint)
          idlv1=idptl(new)
          call idmass(idlv1,pptl(5,new))
          sum=sum+pptl(5,new)
        endif
        wimax=max(wimax,1.01*(amss0-sum))   !adjust max width to minimum mass
 110  continue
 111  continue
      if(ish.ge.4)write(ifch,*)'naddptl,sum,amss',naddptl,sum,amss
      if(naddptl.ne.1.and.sum.ge.amss)goto 2
 112  naddptl1=naddptl-1
      pgen(4,1)=0.
      do 120 j=1,3
      pgen(j,1)=pptl(j,ip)
      pgen(4,1)=pgen(4,1)+pgen(j,1)**2
120   continue
      pgen(5,1)=amss !needed because of deuteron, triton and alpha decay and OK
      pgen(4,1)=sqrt(pgen(4,1)+pgen(5,1)**2)

      pgen(5,naddptl)=pptl(5,nptl+naddptl)
      if(naddptl.eq.1) goto 700            !one body decay
      if(naddptl.eq.2) goto 400            !two body decay

      if(ish.ge.4)write(ifch,*)'>= 3 body decay'

c     use kroll-wada distribution for pi0 and eta dalitz decays.
c     ----------------------------------------------
      if(.not.((idptl(ip).eq.110.or.idptl(ip).eq.220).and.
     1iabs(idptl(nptl+2)).eq.12)) goto 130
      ntry=0             !decay of pi0 or eta into electron
125   ntry=ntry+1
           if(ntry.gt.10)then
      if(ish.ge. 1)then
      call utmsg('hdecay&')
      write(ifch,*)'*****  ntry > 10. iret = 1.'
      write(ifch,*)'***** amee,ree,wtee',amee,ree,wtee
      call utmsgf
      endif
      iret=1
      goto1000
           endif
      amee=twome*(pptl(5,ip)/twome)**rangen()
      ree=(twome/amee)**2
      wtee=(1.-(amee/pptl(5,ip))**2)**3*sqrt(1.-ree)*(1.+.5*ree)
      if(wtee.lt.rangen()) goto 125
      pgen(5,2)=amee
      goto 400
130   continue

c     calculate maximum phase-space weight
c     ------------------------------------
      wtmax=1./reduce(naddptl)
      sum1=pgen(5,1)
      sum2=sum-pptl(5,nptl+1)
      do 200 i=1,naddptl1
      wtmax=wtmax*utpcm(sum1,sum2,pptl(5,nptl+i))
      sum1=sum1-pptl(5,nptl+i)
      sum2=sum2-pptl(5,nptl+i+1)
200   continue

c     generate uniform naddptl-body phase space
c     --------------------------------------
      ntry=0
300   ntry=ntry+1
           if(ntry.gt.10000)then
      if(ish.ge. 0)then
      call utmsg('hdecay&')
      write(ifch,*)'*****  infinite loop (2). iret = 1.'
      write(ifch,*)'***** ip,idptl(ip),pptl(5,ip):'
     *,ip,idptl(ip),pptl(5,ip)
      write(ifch,*)'***** wt,wtmax:',wt,wtmax
      write(ifch,*)'***** i,pgen(5,i),pptl(5,nptl+i),idptl(nptl+i):'
      do i=1,naddptl
      write(ifch,*)i,pgen(5,i),pptl(5,nptl+i),idptl(nptl+i)
      enddo
      call utmsgf
      endif
      iret=1
      goto 1000
           endif
      rnd(1)=1.
      jsave=1
      do 310 i=2,naddptl1
      rnew=rangen()
      i1=i-1
      do 320 jj1=1,i1
      j=i-jj1
      jsave=j+1
      if(rnew.le.rnd(j)) goto310
      rnd(jsave)=rnd(j)
320   continue
310   rnd(jsave)=rnew
      rnd(naddptl)=0.
      wt=1.
      sum1=sum
      do 330 i=2,naddptl
      sum1=sum1-pptl(5,nptl+i-1)
      pgen(5,i)=sum1+rnd(i)*(pgen(5,1)-sum)
      a=pgen(5,i-1)
      b=pgen(5,i)
      c=pptl(5,nptl+i-1)
      wt=wt*utpcm(a,b,c)
330   continue
      if(wt.lt.rangen()*wtmax) goto 300

c     carry out two-body decays in pgen frames
c     ----------------------------------------
400   continue
      if(ish.ge.4)write(ifch,*)'2 body decay'
      do 410 i=1,naddptl1
      !print*,'DKY i id id1-4',i,idptl(ip),mode(1,ipoint),mode(2,ipoint)
      !.,mode(3,ipoint),mode(4,ipoint)
      qcm=utpcmi(pgen(5,i),pgen(5,i+1),pptl(5,nptl+i),iret)
      if(iret.ne.0)then
        write(ifmt,*)'ERROR DKY i id id1-4',i,idptl(ip),mode(1,ipoint)
     .  ,mode(2,ipoint),mode(3,ipoint),mode(4,ipoint)
        stop  
      endif
      u(3)=2.*rangen()-1.
      phi=2.*pi*rangen()
      u(1)=sqrt(1.-u(3)**2)*cos(phi)
      u(2)=sqrt(1.-u(3)**2)*sin(phi)
      do 420 j=1,3
      pptl(j,nptl+i)=qcm*u(j)
      pgen(j,i+1)=-pptl(j,nptl+i)
420   continue
      pptl(4,nptl+i)=sqrt(qcm**2+pptl(5,nptl+i)**2)
      pgen(4,i+1)=sqrt(qcm**2+pgen(5,i+1)**2)
410   continue
      do 430 j=1,4
      pptl(j,nptl+naddptl)=pgen(j,naddptl)
430   continue

c     boost pgen frames to lab frame
c          also save momenta in rest frame (last frame)
c     -------------------------------------------------
      do 500 ii=1,naddptl1
      i=naddptl-ii
      do 510 j=1,3
      beta(j)=pgen(j,i)/pgen(4,i)
510   continue
      gamma=pgen(4,i)/pgen(5,i)
      do 520 k=i,naddptl
      k1=nptl+k
      b=pptl(1,k1)
      if(ish.ge.10)write(ifch,*)'hdecay 500',k1,b
      bp=beta(1)*pptl(1,k1)+beta(2)*pptl(2,k1)+beta(3)*pptl(3,k1)
      do 530 j=1,3
      prest(j,k)=pptl(j,k1)
      a=pptl(j,k1)
      b=gamma*beta(j)*(pptl(4,k1)+bp*gamma/(gamma+1.))
      pptl(j,k1)=a+b
530   continue
      prest(4,k)=pptl(4,k1)
      pptl(4,k1)=gamma*(pptl(4,k1)+bp)
      if(pptl(4,k1).lt.1.d-5)then
        pptl(4,k1)=sqrt(pptl(1,k1)*pptl(1,k1)+pptl(2,k1)*pptl(2,k1)
     &                 +pptl(3,k1)*pptl(3,k1))
      endif
520   continue
500   continue

c     matrix elements
c     ---------------
        if(iabs(idptl(ip)).eq.14)then                  !muon decay
          goto 650
        elseif(naddptl.eq.3.and..not.((idptl(ip).eq.110
     1         .or.idptl(ip).eq.220).and.iabs(idptl(nptl+2)).eq.12))then
          if(idptl(ip).eq.221.or.idptl(ip).eq.331)then  !omeg and phi decay
            goto 610
          elseif(iabs(idptl(ip)).eq.130.or.       !Kl and K decay
     1       idptl(ip).eq.-20)then
            if(iabs(idptl(nptl+2)).lt.20)then   !semi-leptonic
              goto 630
            else                                !hadronic
              goto 640
            endif
          elseif(iabs(idptl(nptl+1)).lt.20.and. !other semi-leptonic decay
     1       idptl(nptl+1).ne.10)then
            goto 620
          elseif(iabs(idptl(nptl+2)).le.6)then
            goto 605            !decay into quark
          else
            goto 800
          endif
        else
         goto 800
        endif

 605    wt=pptl(5,ip)*pptl(5,nptl+1)*dot(2,3)
        IF(wt.LT.rangen()*pptl(5,ip)**4/16.) goto 300
        ams=sqrt(dot(2,2)+dot(3,3)+2.*dot(2,3))
        kno(1)=idptl(nptl+2)
        kno(2)=idptl(nptl+3)
        if(ammin(kno(1),kno(2)).gt.ams)then
          call vedi(kno(1),kno(2),iddum,idlv2)
          idptl(nptl+2)=idlv2
          call idmass(idlv2,amnew2)
          pptl(5,nptl+2)=amnew2
          naddptl=2
          goto 112
        endif
c......multiplicity
        PS=0
        if(dot(2,2).gt.0.)PS=sqrt(dot(2,2))
        psq=0
        if(dot(3,3).gt.0.)psq=sqrt(dot(3,3))
c        PSP=PS                  !!???????????????unused
        np=0                    !!!!?????
        nq=2
        CNDE=4.5*LOG(MAX((ams-PS-PSQ)/0.7,1.1))
c        IF(MMAT.EQ.12) CNDE=CNDE+PARJ(63)
 769    NTRY=NTRY+1
        IF(NTRY.GT.1000) THEN
          write(*,*)'hdecay caught in infinite loop'
          write(ifch,*)'hdecay caught in infinite loop'
          iret=1
          goto 1000
        ENDIF
        GAUSS=SQRT(-2.*CNDE*LOG(MAX(1E-10,rangen())))*
     &       SIN(2.*pi*rangen())
        ND=0.5+0.5*NP+0.25*NQ+CNDE+GAUSS
        IF(ND.LT.NP+NQ/2.OR.ND.LT.2.OR.ND.GT.10) GOTO 769


c......choose hadrons


        kno(3)=kno(1)
        kno(4)=kno(2)

        CONTINUE
        IF(ND.EQ.NP+NQ/2) GOTO 773
        DO I=nptl+2,nptl+2+nd-nq/2-1
          JT=2+1+INT((NQ-1) * rangen() )
          CALL vedi(kno(JT),0,KFL2,idlv3)
          idptl(i)=idlv3
c          IF(K(I,2).EQ.0) GOTO 769
          kno(JT)=-KFL2
        enddo
 773    CONTINUE
        CALL vedi(kno(3),kno(4),KFLDMP,idlv4)
        idptl(nptl+2+nd-nq/2)=idlv4
        sum=0.
        do i=nptl+2,nptl+2+nd-nq/2
          call idmass(idptl(i),am)
          pptl(5,i)=am
          sum=sum+am
        enddo
        if(sum.gt.ams) goto 769
c......goto phase space dis....
        ipxx=nptl+2+nd-nq/2+1
        do j=1,4
          pptl(j,ipxx)=pptl(j,ipp)-pptl(j,nptl+1)
        enddo
        pptl(5,ipxx)=ams
        idptl(ipxx)=sign(80,idptl(ipp))
        nptl=nptl+1
        naddptl=nd
        goto 112


c     omeg and phi decay
c          use vectors in rest frame
c     ------------------------------
610   wt=(pptl(5,nptl+1)*pptl(5,nptl+2)*pptl(5,nptl+3))**2
     1-(pptl(5,nptl+1)*dot(2,3))**2
     2-(pptl(5,nptl+2)*dot(1,3))**2
     3-(pptl(5,nptl+3)*dot(1,2))**2
     4+2.*dot(1,2)*dot(2,3)*dot(1,3)
      if(wt.lt.rangen()*pptl(5,ip)**6/108.) goto300
      goto 800

c     semileptonic and quark decays
c          use vectors in rest frame, where ip has (m,0,0,0)
c          include w propagator
c     ------------------------------------------------------
620   wt=(pptl(5,ip)*prest(4,2))*dot(1,3)
      s12=pptl(5,nptl+1)**2+pptl(5,nptl+2)**2+2.*dot(1,2)
      s12max=pptl(5,ip)**2
      wt=wt*wprop(s12max)/wprop(s12)
      if(wt.lt.rangen()*pptl(5,ip)**4/16.) goto 300
      goto 800

c     semileptonic kaon decays
c          use vectors in rest frame, where ip has (m,0,0,0)
c          include form factor FML
c     ------------------------------------------------------
630   if(iabs(idptl(ip)).eq.130)then
        if(iabs(idptl(nptl+2)).eq.12)then
          ncha=1          !K   -> Pi0 + e + Nu
        else
          ncha=2          !K   -> Pi0 + Mu + Nu
        endif
      else
        if(iabs(idptl(nptl+2)).eq.12)then
          ncha=3          !K0  -> Pi + e + Nu
        else
          ncha=4          !K0  -> Pi + Mu + Nu
        endif
      endif

      wt=FML(ncha,pptl(5,ip),pptl(5,nptl+1),pptl(5,nptl+2)
     &       ,prest(4,1),prest(4,2),prest(4,3))
      if(wt.lt.rangen()) goto 300
      goto 800

c     hadronic kaon decays
c          use vectors in rest frame, where ip has (m,0,0,0)
c          include form factor FM
c     ------------------------------------------------------
640   if(iabs(idptl(ip)).eq.130)then
        if(iabs(idptl(nptl+3)).eq.120)then
          ncha=1          !K   -> 3 Pi
        else
          ncha=2          !K   ->  Pi + 2 Pi0
        endif
      else
        if(iabs(idptl(nptl+1)).eq.110)then
          ncha=3          !K0  -> 3 Pi0
        else
          ncha=4          !K0  -> 2 Pi + Pi0
        endif
      endif
      S0=(pptl(5,ip)**2+pptl(5,nptl+1)**2+pptl(5,nptl+2)**2
     &   +pptl(5,nptl+3)**2)/3.d0
      S1=pptl(5,ip)**2+pptl(5,nptl+1)**2-2.*prest(4,1)*pptl(5,ip)
      S2=pptl(5,ip)**2+pptl(5,nptl+2)**2-2.*prest(4,2)*pptl(5,ip)
      S3=pptl(5,ip)**2+pptl(5,nptl+3)**2-2.*prest(4,3)*pptl(5,ip)
      wt=FM(ncha,S0,S1,S2,S3)
      if(wt.lt.rangen()) goto 300
      goto 800

c     muon decays
c          use vectors in rest frame, where ip has (m,0,0,0)
c          include form factor FMU
c     ------------------------------------------------------
650   xxx=2.*prest(4,1)/pptl(5,ip)            !reduced energy of electron
      if(xxx.gt.1.) goto 300
      wt=FMU(xxx)
      rrr=rangen()
      if(wt.lt.rrr) goto 300
      goto 800

c     one-particle decays
c     -------------------
700   continue
      do 710 j=1,5
      pptl(j,nptl+1)=pptl(j,ip)
710   continue

c     swap particles and antiparticles if idptl(ip)<0
c     -----------------------------------------------
 800    continue
        if(iabs(idptl(ip)).eq.80)then
          nptl=nptl-1
          naddptl=naddptl+1
        endif
        !update values for printout
        do i=1,naddptl
          istptl(nptl+i)=0
          iorptl(nptl+i)=ipp
          jorptl(nptl+i)=0
          ifrptl(1,nptl+i)=0
          ifrptl(2,nptl+i)=0
          ityptl(nptl+i)=ityptl(ipp)
        enddo
        if(idptl(ipp).ge.0.or.iabs(idptl(ipp)).eq.20) goto 900
        do 810 i=1,naddptl
          idabs=iabs(idptl(nptl+i))
          !ifl1=idabs/1000  !cKW2108 This is old!!!
          !ifl2=mod(idabs/100,10)
          !ifl3=mod(idabs/10,10)
          call idflav(idabs,ifl1,ifl2,ifl3,jdummy,idummy)  !cKW2108 New version
          if(ifl1.eq.0.and.ifl2.ne.0.and.ifl2.eq.-ifl3) goto 810
          if(idabs.eq.9.or.idabs.eq.10.or.idabs.eq.20) goto 810
          if(idabs.eq.29.or.idabs.eq.30.or.idabs.eq.40) goto 810
          idptl(nptl+i)=-idptl(nptl+i)
 810    continue

 900    continue
c save mass (and energy) as used for decay (mass shift due to lifetime)
        pptl(4,ip)=pgen(4,1)
        pptl(5,ip)=pgen(5,1)
        nptl=nptl+naddptl
        if(nptl.gt.mxptl)call utstop('hdecay: nptl>mxptl&')
c        nqk=0           !???????????????????unused
        if(iabs(idptl(nptl)).lt.10.or.mod(idptl(nptl),100).eq.0)then
c          call utstop('hdecay: decay ptcl is parton&')
        endif

c     print
c     -----

      if(ish.ge.3)then
      write(ifch,140)sngl(ttaus)
  140 format(/' ----------------------------'/
     *'    decay  at tau =',f6.2/
     *' ----------------------------')
      write(ifch,*)'decaying object:'
      call alist('&',ip,ip)
      write(ifch,*)'decay products:'
      call alist('&',nptlb+1,nptl)
      endif
      if(ish.ge.5)then
      write(ifch,*)'momentum sum:'
      do kk=1,5
      pptl(kk,nptl+1)=0
      do ii=nptlb+1,nptl
      pptl(kk,nptl+1)=pptl(kk,nptl+1)+pptl(kk,ii)
      enddo
      enddo
      call alist('&',nptl+1,nptl+1)
      endif

c     exit
c     ----

 1000 continue
      ip=ipp
      if(iret.ne.0.and.ish.ge.1)then
        write(ifmt,'(a)')'hdecay: redo event'
        write(ifch,'(a)')'hdecay: redo event'
c        stop
      endif
      call utprix('hdecay',ish,ishini,5)
      return

 9999   call utstop('hdecay: mxptl too small&')
        end

c---------------------------------------------------------------------
      subroutine vedi(k1,k2,k3,id)
c---------------------------------------------------------------------
      include 'epos.inc'
      if(k2.eq.0)then
        if(rangen().lt.pdiqua.and.iabs(k1).lt.6)then
          ifl1=int(rangen()/pud)+1
          ifl2=int(rangen()/pud)+1
          k3=-min(ifl1,ifl2)*1000-max(ifl1,ifl2)*100
        else
          k3=int(rangen()/pud)+1
        endif
        if(k1.gt.0.and.k1.le.6)k3=-k3
        if(k1.lt.-1000)k3=-k3
      else
        k3=k2
      endif
      id=idsp(k1,k3)
      if(iabs(id).le.999) then
        ids=max(mod(iabs(id)/100,10),mod(iabs(id)/10,10))
        if(ids.le.2)then
          idr=sign(iabs(id)+int(rangen()+0.5),id)
        elseif(ids.eq.3)then
          idr=sign(iabs(id)+int(rangen()+0.6),id)
        else
          idr=sign(iabs(id)+int(rangen()+0.75),id)
        endif
      else
        idr=sign(iabs(id)+int(0.5+rangen()),id)
      endif
      id=idr
      if(ish.ge.5)write(ifch,*) 'Flavor:',k1,k2,k3,id
      end

c-----------------------------------------------------------------------
      subroutine hdecin(lprint)
c-----------------------------------------------------------------------
c     sets up /dkytab/
c-----------------------------------------------------------------------
      include 'epos.inc'
      common/wco/wmass2,wgam2
      dimension imode(6)
      character*8 idlabl,lmode(6),lres
      character*8 iblank
      logical lprint
      parameter (mxlook=10000,mxdky=4000)
      common/dkytab/look(mxlook),cbr(mxdky),mode(5,mxdky)
      common/nodcay/nodcay,noeta,nopi0,nonunu,noevol,nohadr
      logical nodcay,noeta,nopi0,nonunu,noevol,nohadr
      parameter (ndectbmax=4000)
      real dectab(7,ndectbmax)

      write(ifmt,'(2a)')'load ',fnnx(1:nfnnx)//'/idky5.dt'
      open(95,file= fnnx(1:nfnnx)//'/idky5.dt',STATUS='UNKNOWN')

      j=1
      idlast=0
      do
        read(95,*, end=999)dectab(1,j),dectab(2,j),dectab(3,j)
     .  ,dectab(4,j),dectab(5,j),dectab(6,j),dectab(7,j)
        id=nint(dectab(1,j))
        if(id.eq.idlast)then
          dectab(2,j)=dectab(2,j)+dectab(2,j-1)
        elseif(j.gt.1)then
          if(abs(dectab(2,j-1)-1.0).gt.1e-4)then
            write(ifmt,'(a,$)')'WARNING DKY'
            write(ifmt,'(i5,f8.5,$)')nint(dectab(1,j-1)),dectab(2,j-1)
            write(ifmt,'(1x,a1,i4,f8.5,a1,$)')
     .      '(',nint(dectab(1,j-2)),dectab(2,j-2),')'
            dectab(2,j-1)=1.000000
            write(ifmt,'(a,$)')' changed to'
            write(ifmt,'(i5,f8.5)')nint(dectab(1,j-1)),dectab(2,j-1)
          endif
        endif 
        idlast=id
        j=j+1
        if(j.gt.ndectbmax)stop'ndectbmax too small !'
      enddo      
999   continue
      ndectb=j-1
      !print*,ndectb !Check ndectb
      close(95)

c     determine wmass2,wgam2
c     ----------------------
      alfa=1./137.036
      gf=1.16570e-5
      sin2w=.215
      sinw=sqrt(sin2w)
c      cosw=sqrt(1.-sin2w)           !?????????????????unused
      amw=sqrt(pi*alfa/(.9304*sqrt(2.)*gf))/sinw
      wmass2=amw
      call idmass(5,amlep5)
      call idmass(6,amlep6)
      ngam=12
      if(amlep5+amlep6.gt.amw) ngam=9
      wgam2=gf*amw**3/(6.*pi*sqrt(2.))*ngam

      data iblank/' '/
      ird=0
      do 1 i=1,mxlook
1     look(i)=0
      do 2 i=1,mxdky
      do 3 j=1,5
3     mode(j,i)=0
2     cbr(i)=0.
      nodcay=.false.
      noeta=.false.
      nopi0=.false.
      nonunu=.false.
      noevol=.false.
      nohadr=.false.
      if(lprint) write(ifch,10)
10    format('1',30('*')/' *',28x,'*'/
     1' *',5x,'isajet decay table',5x,'*'/
     2' *',28x,'*'/' ',30('*')//
     36x,'part',18x,'decay mode',19x,'cum br',15x,'ident',17x,
     4'decay ident')
      loop=0
      iold=0
      if(nodcay) return

200   loop=loop+1
      if(loop.gt.mxdky) goto9999
220   do 210 i=1,5
      imode(i)=0
      lmode(i)=iblank
210   continue
      ird=ird+1
      if(ird.gt.ndectb)return
c      if(ird.gt.1171)return   ! ??????????????????????????
      ires=nint(dectab(1,ird))
      br=dectab(2,ird)
      do 215 i=1,5
215   imode(i)=nint(dectab(2+i,ird))
      if(nopi0.and.ires.eq.110) goto220
      if(noeta.and.ires.eq.220) goto220
      if(ires.eq.iold) goto230
      if(ires.lt.0.or.ires.gt.mxlook)
     *call utstop('hdecin: ires out of range&')
      look(ires)=loop
230   iold=ires
      cbr(loop)=br
      do 240 i=1,5
      mode(i,loop)=imode(i)
      if(imode(i).ne.0) lmode(i)=idlabl(imode(i))
240   continue
      lres=idlabl(ires)
      if(lprint) write(ifch,20) lres,(lmode(k),k=1,5),
     1br,ires,(imode(k),k=1,5)
20    format(6x,a5,6x,5(a5,2x),3x,f8.5,15x,i5,4x,5(i5,2x))
      goto200

9999  write(ifch,*)'loop=', loop
      call utstop('hdecin: loop > mxdky&')

      end

C -----------------------------------------------
      FUNCTION FM(NQ,S0,S1,S2,S3)
C -----------------------------------------------
C Normalized TRANSITION MATRIX FOR THE DALIZT PLOT DISTRI.
C OF K -> 3 PIONS. PARAMETRIZATION OF WEINBERG
C AS DESCRIBE IN PARTICLE DATA BOOK.
C G IS THE LINEAR COEFFICIENT (SLOPE g)
C H IS THE QUADRATIC COEFFICIENT h
C D IS THE QUADRATIC COEFFICIENT k
C Amax is the maximum of this amplitude (taken from Corsika by D. Heck)
C NQ is the decay channel :
C   1 - K -> 3 Pi
C   2 - K -> Pi + 2 Pi0
C   3 - K0 -> 3 Pi0
C   4 - K0 -> 2 Pi + Pi0
C -----------------------------------------------
      DIMENSION G(4),H(4),D(4),Amax(4)
      PARAMETER (PIM=139.57E-3)
      DATA G/-0.2154,0.594,0.,0.67/
      DATA H/0.01,0.035,0.,0.079/
      DATA D/-0.01,0.,0.,0.0098/
      DATA Amax/1.27,1.84,1.,2.22/

      FM=1.+G(NQ)*(S3-S0)/(PIM*PIM)+H(NQ)*((S3-S0)/(PIM*PIM))**2
     *+D(NQ)*((S2-S1)/(PIM*PIM))**2
      FM=FM/Amax(NQ)

      RETURN
      END
C -----------------------------------------------
      FUNCTION FML(N,AM,RM1,RM2,E1S,E2S,E3S)
C -----------------------------------------------
C Normalized DALITZ PLOT DENSITY (RHO)
C OF K -> 1 PION + 2 LEPTONS
C AS DESCRIBE IN PARTICLE DATA BOOK.
C CLP IS THE LAMBDA + FORM FACTOR COEFFICIENT
C CLN IS THE LAMBDA 0 FORM FACTOR COEFFICIENT
C EEP IS E'pion
C GP IS THE F+(t) FORM FACTOR (t=AM*AM+SM1-2.D0*AM*E1S)
C H IS EPS(t)=F-(t)/F+(t) WHERE F- IS CALCULATED FROM F0
C Amax is the maximum of this density (taken from Corsika by D. Heck)
C N is the decay channel :
C   1 - K -> Pi0 + e + Nu
C   2 - K -> Pi0 + Mu + Nu
C   3 - K0 -> Pi + e + Nu
C   4 - K0 -> Pi + Mu + Nu
C -----------------------------------------------
      DIMENSION CLP(4),CLN(4),Amax(4)
      DATA CLP/0.0276,0.031,0.0288,0.034/
      DATA CLN/0.0,0.006,0.,0.025/
      DATA Amax/1.28e-2,1.194e-2,1.31e-2,1.241e-2/

      SM1=RM1*RM1
      SM2=RM2*RM2
      EEP=0.5D0*(AM*AM+SM1-SM2)/AM-E1S
      GP=1.+CLP(N)*(AM*AM+SM1-2.*AM*E1S)/SM1
      H=(AM*AM-SM1)/SM1*(CLN(N)-CLP(N))/GP
      FML=GP*GP*(AM*(2.*E2S*E3S-AM*EEP)+
     *SM2*(0.25*EEP-E3S)+H*SM2*(E3S-0.5*EEP)+
     *0.25*H*H*SM2*EEP)
      FML=FML/Amax(N)
      RETURN
      END
C -----------------------------------------------
      FUNCTION FMU(X)
C -----------------------------------------------
C PROBABILITY DISTRI. FOR ELECTRON ENERGY FROM MUON DECAY :
C MU -> 2NU + E. DESCRIBE IN PARTICLE DATA BOOK.
C (SIMPLIFY DIFFERENTIAL DECAY RATE INTEGRATED)
C X REDUCED ENERGY OF PARTICLE
C -----------------------------------------------

      FMU=2.*(3.-2.*X)*X*X

      RETURN
      END
