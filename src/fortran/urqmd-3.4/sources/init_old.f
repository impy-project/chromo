c $Id: init.f,v 1.22 2007/01/30 14:50:25 bleicher Exp $
C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine uinit(io)
c
c Revision : 1.0
c
cinput io : flag for call to {\tt input(io)} 
c
c This subroutine calls initialization procedures for uqmd
c i.e. random generater, etc. 
c Routines called before the first (physical) event of uqmd should
c enter here. an exception is the subroutine {\tt init}. {\tt init} should NOT
c be included in uinit. 
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit none
      include 'coms.f'
      include 'options.f'
      include 'boxinc.f'
      include 'inputs.f'
      integer io

      call set0
      call params

c.. display logo
      call urqmdlogo

c read input file
      call input(io)

      call strini

      call output(19)

      firstseed=.true.
      fixedseed=ranseed.gt.0
      if(fixedseed)write(6,*)'fixed random number:',ranseed
      call sseed(ranseed)
      call loginit
      if(CTOption(33).eq.0.or.CTOption(9).eq.0) call loadwtab(io)
c in case of CASCADE mode, the potentials need not be calculated
      if(EoS.ne.0) then
         if(logSky) call potdww
         if(logYuk) call potYuk
         if(logPau) call potPau
         if(logCb)  call potCb
      endif
c
c  calculate the normalization of resonances distribution...
c
      call norm_init

      if(io.ne.0) return

c     do not initialize projectile and target if old event is read in
      if(CTOption(40).ne.0) return
c SEH: in urqlnk in corsika
      if(boxflag.eq.0) then
c
c initialize nuclei (projectile and target) and store them
c
c initialize normal projectile
         if(prspflg.eq.0) then
c         if(eos.eq.0) then
            call cascinit(Zp,Ap,1)
c         else
c            write(6,*)'illegal EOS in init.'
c            stop 137
c         endif
         endif
c initialize normal target
         if(At.ne.0) then
            if(trspflg.eq.0) then
c            if(eos.eq.0) then
               call cascinit(Zt,At,2)
c            else
c               write(6,*)'illegal EOS in init.'
c               stop 137
c            endif
            endif
         endif
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine init
c
c     Revision : 1.0
c     This subroutine calls initialization procedures for different
c     equations of state and calculation modi.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'coms.f'
      include 'options.f'
      include 'comres.f'
      include 'inputs.f'
      include 'freezeout.f'
      include 'newpart.f'
      include 'boxinc.f'
      include 'colltab.f'

      integer j,k,icount,npold,getspin,fchg,indsp(2),isrt,ipbm
      real*8 nucrad,dstt,dstp,pcm,eb,embeam,emtarget
      real*8 massit,ranf,pcms,dectim
      real*8 gaeq,beeq,galab,belab,ppeq,pteq
      real*8 pboost
      real*8 ratio
      integer AAp, AAt
        integer nnuc
        parameter (nnuc=11)
        save isrt,ipbm
        logical bcorr
        common /ini/ bcorr

        integer i,flagy
        real*8 alf,regula
c momenta
       real*8 mx,my,mz  
c important: never change!
      npart = 0
      npold = 0
      nbar=0
      nmes=0
      apt=0
      uid_cnt=0
c reset counters
c     all collisions/decays
      ctag  = 0
c     all decays
      dectag = 0
c     number of prod. hard resonances
      NHardRes=0
c     number of prod. soft resonances
      NSoftRes=0
c     number of prod. resonances via decay
      NDecRes=0
c     number of blocked collisions
      NBlColl=0
c     number of elastic collisions
      NElColl=0
c
      eb=0D0
c icount is the number of EXTRAordinary pro/tar combinations (i.e. pion ...)
      icount = 0
c reset particle vectors
      do 20 j=1,nmax
        spin(j)  = 0
        ncoll(j) = 0
        lstcoll(j)=0
        r0(j) = 0.0
        rx(j)    = 0.0
        ry(j)    = 0.0
        rz(j)    = 0.0
        p0(j)    = 0.0
        px(j)    = 0.0
        py(j)    = 0.0
        pz(j)    = 0.0
        frr0(j) = 0.0
        frrx(j)    = 0.0
        frry(j)    = 0.0
        frrz(j)    = 0.0
        frp0(j)    = 0.0
        frpx(j)    = 0.0
        frpy(j)    = 0.0
        frpz(j)    = 0.0
        fmass(j) = 0.0
        charge(j)= 0
        iso3(j)  = 0
        ityp(j)  = 0
        dectime(j)= 0.0
        origin(j)=0
        tform(j)=0.0
        xtotfac(j)=1.0
        uid(j)=0
         ffermpx(j) = 0.0
         ffermpy(j) = 0.0
         ffermpz(j) = 0.0
ctd
         do 21 k=1,2
            p0td(k,j)=0.d0
            pxtd(k,j)=0.d0
            pytd(k,j)=0.d0
            pztd(k,j)=0.d0
            fmasstd(k,j)=0.d0
            ityptd(k,j)=0
            iso3td(k,j)=0
 21      continue
 20   continue


      if(CTOption(40).ne.0) then
         call getoldevent
         return
      endif


      if (boxflag.eq.1) then 
         mbpx=0
         mbpy=0
         mbpz=0
         mx=0
         my=0
         mz=0

         nbar=0
         nmes=0
         flagy=edensflag
         ctoption(30)=0         ! no frozen fermi for box
         do 100 cbox=1,mbox
            if (flagy.ge.1) then
               bptpmax(cbox)=edens/mbox
            endif       
            call bptinit(cbox)
 100     Continue

c prevents a collective motion in the box
         do 143 i=1,npart
            mbpx=mbpx+px(i)
            mbpy=mbpy+py(i)
            mbpz=mbpz+pz(i)
 143     continue
         do 142 i=1,npart
            px(i)=px(i)-mbpx/npart
            py(i)=py(i)-mbpy/npart
            pz(i)=pz(i)-mbpz/npart
            call setonshell(i)
 142     continue
         
         if (flagy.ge.1) then
            alf=Regula(edens)
            do 42  i=1,npart   
               px(i) = alf*px(i)
               py(i) = alf*py(i)
               pz(i) = alf*pz(i)
               call setonshell(i)
 42         continue
         endif
         Write(*,*) 'walls selected'
         mbflag=2
         if (solid.gt.0) Write(*,*)'solid walls selected'
         
         Return
      EndIf

      if(At.ne.0) then
         if (CTParam(21).eq.0.0) then
            dstp = nucrad(Ap)+CTParam(41)
            dstt = nucrad(At)+CTParam(41)
         else
            ratio = sqrt((1 + 4.0*CTParam(21)/3.0) / 
     $           (1 - 2.0*CTParam(21)/3.0) )
            dstp = nucrad(Ap)*ratio**(2.0/3.0)+CTParam(41)
            dstt = nucrad(At)*ratio**(2.0/3.0)+CTParam(41)
         endif
c add radius offset
         dstp=dstp+CTParam(30)
         dstt=dstt+CTParam(30)      
c         
c         dst=0.5d0*(dstt+dstp)
      else
c            dst=0.d0
            dstp=0d0
            dstt=0d0
      endif

ce For anti nuclei:
      AAp = abs(Ap)
      AAt = abs(At)

c
c fix masses of projectile and target for calculation of pbeam,ecm,pcm
      if(prspflg.eq.0) then
         embeam=AAp*EMNUC
      elseif(prspflg.eq.1) then
         icount=icount+1
c!!!  sofar only pro/tar with fixed masses allowed
         embeam=massit(spityp(1))
      endif
      if(trspflg.eq.0) then
         emtarget=AAt*EMNUC
      elseif(trspflg.eq.1) then
         icount=icount+1
         emtarget=massit(spityp(2))
      endif
c         
c p(equal_speed) with given elab  cccccccccccccccccccccccccccccccccccccc

      if(srtflag.eq.0) then
         
         ebeam=AAp*ebeam
         eb=ebeam+embeam
         pbeam=sqrt(ebeam*(ebeam+2.0d0*embeam))
         
c     p(equal_speed) with given sqrt(s) ccccccccccccccccccccccccccccccccccccc
         
      elseif(srtflag.eq.1) then
         
c     in this mode, everything has to calculated on a per particle basis
         embeam=embeam/dble(AAp)
         emtarget=emtarget/dble(AAt)
         
         if(emtarget+embeam.gt.srtmin)then
            srtmin=emtarget+embeam+1d-2
            write(6,*)' *** error:initial energy below treshold'
            write(6,*)'     c.m. energy will be increased to:',
     &           srtmin
            srtmax=max(srtmax,srtmin)
         end if
         if(efuncflag.eq.0) then
            ecm=srtmin
         else                   ! if(efuncflag.eq.1) then
c     excitation function
            if(mod((event-1)*nsrt,nevents).eq.0
     &           .and.firstev.gt.0)then
               if(efuncflag.eq.1)then 
                  ecm=ecm+(srtmax-srtmin)/dble(nsrt-1)
               else if(efuncflag.eq.2)then
                  isrt=isrt+1
                  ecm=srtmin*exp( 
     &                 (dlog(srtmax/srtmin)) 
     &                 *isrt/(nsrt-1))
               end if
            elseif(firstev.eq.0) then
               isrt=0
               firstev=1
               ecm=srtmin
            endif
         endif
c     
c     this is all on a per particle basis
         pcm=pcms(ecm,embeam,emtarget)
         ebeam=sqrt(embeam**2 + (pcm*ecm/emtarget)**2) - embeam
         pbeam=pcm*ecm/emtarget
c     now revert to full quantities
         ebeam=AAp*ebeam
         pbeam=AAp*pbeam
         embeam=embeam*AAp
         emtarget=emtarget*AAt
         eb=sqrt(embeam**2+pbeam**2)
         

         
c     p(equal_speed) with given plab ccccccccccccccccccccccccccccccccccc
      elseif(srtflag.eq.2) then
         if(efuncflag.gt.0) then
c     excitation function
c     calculate the next pbeam if number of events at current pbeam exceeds nevents/npb 
            if (mod((event-1)*npb,nevents).eq.0
     &           .and.firstev.gt.0) then
c     excitation function linear in pbeam
               if (efuncflag.eq.1) then 
                  pbeam=pbeam+(pbmax-pbmin)/dble(npb-1)
c     else use a logaritmic excitation function
               else if(efuncflag.eq.2) then
                  ipbm=ipbm+1
                  pbeam=pbmin*exp( 
     &                 (dlog(pbmax/pbmin)) 
     &                 *ipbm/(npb-1))
               end if
            else if (firstev.eq.0) then
               ipbm=0
               firstev=1
               pbeam=pbmin
            endif
         endif
c     input was pbeam per particle
         pbeam=AAp*pbeam
         eb=sqrt(embeam**2+pbeam**2)
         ebeam=eb-embeam
      endif   
      
c     now do the calculation of equal_speed quantities
      galab=eb/embeam           ! gamma_lab
      belab=pbeam/eb            ! beta_lab
      gaeq=sqrt(0.5*(1+galab))
      beeq=belab*galab/(1+galab)
      ppeq=gaeq*beeq*embeam     ! p_projectile(eq)
      pteq=-(gaeq*beeq*emtarget)
      
c     reduce to per particle quantities
      ppeq=ppeq/dble(AAp)
      if(AAt.ne.0) then
         pteq=pteq/dble(AAt)
         emtarget=emtarget/dble(AAt)
      endif
      embeam=embeam/dble(AAp)
      pbeam=pbeam/dble(AAp)
      ebeam=ebeam/dble(AAp)
c     the following is the NN sqrt(s)
      ecm=sqrt(embeam**2+2*eb/dble(AAp)*emtarget+emtarget**2)  
ccccccccccccccccccccccccccccccccccccccccccccc
c compute transformation betas for output

         pcm=max(1d-10,pbeam*emtarget/ecm)

         if(CTOption(27).eq.0) then
            betann=0.d0
            betatar=pcm/sqrt(pcm**2+emtarget**2)
            betapro=-(1.*pcm/sqrt(pcm**2+embeam**2))
         elseif(CTOption(27).eq.1) then
            betann=-(1*pcm/sqrt(pcm**2+emtarget**2))
            betatar=0.d0
            betapro=-(1*pbeam/sqrt(pbeam**2+embeam**2))
         elseif(CTOption(27).eq.2) then
            betann=pcm/sqrt(pcm**2+emtarget**2)
            betatar=pbeam/sqrt(pbeam**2+embeam**2)
            betapro=0.d0
         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c determine impact parameter
         if(CTOption(5).eq.0) then
            bimp=bdist
         elseif(CTOption(5).eq.1) then
            if(bdist.gt.(nucrad(Ap)+nucrad(At)+2*CTParam(30))) 
     &           bdist=nucrad(Ap)+nucrad(At)+2*CTParam(30)
            bimp=sqrt(bmin**2 + ranf(0) * (bdist**2 - bmin**2))
         elseif(CTOption(5).eq.2) then
            if(bdist.gt.(nucrad(Ap)+nucrad(At)+2*CTParam(30))) 
     &           bdist=nucrad(Ap)+nucrad(At)+2*CTParam(30)
            bimp=bmin+ranf(0)*(bdist-bmin)
         else
            write(6,*)'illegal CTOption(5) :',CTOption(5)
            stop 137
         endif

         if(At.eq.0) bimp=0.d0

c initialize normal projectile
         if(prspflg.eq.0) then
            if(mod(event,nnuc).eq.0)then
               call cascinit(Zp,Ap,1)
            endif
            call getnucleus(1,npart)
            npart=npart+AAp
c change reference frame
            if (CTOption(27).eq.1) then   
               pboost = -pbeam
            elseif (CTOption(27).eq.2) then
               pboost = 0.d0
            else
               pboost = -ppeq
            endif  
            call boostnuc(npold+1,npold+AAp,
     &                    pboost,0.5*bimp,-dstp)
c save fermi motion
            if (CTOption(30).eq.1) then
               call savefermi(npold+1,npold+AAp,-pboost)
            endif
            npold=npart
            nbar=nbar+AAp
            if (CTParam(20).ne.0) then
               call getnucleus(1,npart)
               npart=npart+AAp
               call boostnuc(npold+1,npold+AAp,
     &              pboost,0.5*bimp,-dstp+CTParam(20))
               if (CTOption(30).eq.1) then
                  call savefermi(npold+1,npold+AAp,-pboost)
               endif
               npold=npart
               nbar=nbar+AAp
            endif
         endif


c initialize normal target
         if(At.ne.0) then
            if(trspflg.eq.0) then
               if(mod(event,nnuc).eq.0)then
                  call cascinit(Zt,At,2)
               endif
               call getnucleus(2,npart)
               npart=npart+AAt
c change ref. frame
               if(CTOption(27).eq.1) then  
                  pboost = 0.d0
               elseif(CTOption(27).eq.2) then
                  pboost = pbeam
               else
                  pboost = -pteq
               endif
               call boostnuc(npold+1,npold+AAt,
     &                       pboost,-(0.5*bimp),dstt)
c save fermi motion
               if (CTOption(30).eq.1) then
                  call savefermi(npold+1,npold+AAt,-pboost)
               endif
               npold=npart
               nbar = nbar + AAt
            endif
         endif

c set unique ID-tag counter (is not initialized with getnucleus calls)
         uid_cnt=npart

         if(icount.eq.0) then
c set counter for collupd 
            apt=Ap
            return
c initialize special PRO/TAR combinations
         elseif(icount.eq.1) then
            if(prspflg.eq.1) then
               indsp(1)=1
c the "regular" target sits first in the arrays
               apt=At
            else
               indsp(1)=2
               apt=Ap
            endif
         elseif(icount.eq.2) then
            if(abs(spityp(1)).le.abs(spityp(2))) then
               indsp(1)=1
               indsp(2)=2
               apt=Ap
            else
               indsp(1)=2
               indsp(2)=1
               apt=At
            endif
         endif
         do 40 j=1,icount
            npart=npart+1
            if(abs(spityp(indsp(j))).lt.minmes) then
               nbar=nbar+1
            else
               nmes=nmes+1
            endif
            iso3(npart) = spiso3(indsp(j))
            ityp(npart) = spityp(indsp(j))
            uid_cnt=uid_cnt+1
            uid(npart)=uid_cnt
            spin(npart) = getspin(ityp(npart),-1)
            charge(npart)=fchg(iso3(npart),ityp(npart))
            fmass(npart) = massit(ityp(npart))
            rx(npart) = 0.d0
            ry(npart) = 0.d0
            rz(npart) = 0.d0
            px(npart) = 0.d0 
            py(npart) = 0.d0
c       pz ist stored in pbeam,p?eq! 
            pz(npart) = 0.d0
            p0(npart)=sqrt(px(npart)**2+py(npart)**2+pz(npart)**2
     &           +fmass(npart)**2)
            if(indsp(j).eq.1) then
               if(CTOption(27).eq.1) then   
                  call boostnuc(npart,npart,-pbeam,0.5*bimp,-dstp)
               elseif(CTOption(27).eq.2) then
                  call boostnuc(npart,npart,0.d0,0.5*bimp,-dstp)
               else
                  call boostnuc(npart,npart,-ppeq,0.5*bimp,-dstp)
               endif
            elseif(indsp(j).eq.2) then
               if(CTOption(27).eq.1) then
                  call boostnuc(npart,npart,0.d0,-(0.5*bimp),dstt)
               elseif(CTOption(27).eq.2) then
                  call boostnuc(npart,npart,pbeam,-(0.5*bimp),dstt)
               else
                  call boostnuc(npart,npart,-pteq,-(0.5*bimp),dstt)
               endif
            endif
            dectime(npart) = dectim(npart,1)

 40      continue
      return
      end

      subroutine savefermi(i1,i2,p)
c
c     Revision : 1.0
c     Store fermi momentum in fferm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'coms.f'
      
      integer i,i1,i2
      real*8 p

      if (i1.eq.0) return

      if (ncoll(i1).gt.0) return        

      do i=i1,i2
         ffermpx(i)=px(i)
         ffermpy(i)=py(i)
         ffermpz(i)=pz(i)-p
         px(i)=0.0
         py(i)=0.0
         pz(i)=p
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine addfermi(ind,p)
c
c     Revision : 1.0
c     Restore fermi momentum from fferm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'coms.f'
      
      integer ind
      real*8 p

      if (ind.eq.0) return
      
      p = pz(ind)
      px(ind) = px(ind)+ffermpx(ind)
      py(ind) = py(ind)+ffermpy(ind)
      pz(ind) = pz(ind)+ffermpz(ind)
      ffermpx(ind) = 0.0
      ffermpy(ind) = 0.0
      ffermpz(ind) = 0.0

      return
      end

