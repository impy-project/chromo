c$Id: scatter.f,v 1.26 2007/01/30 14:50:27 bleicher Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine scatter(id1,id2,ssigtot,sqrts,ccolfac)
c
c  Revision : 1.0
c
cinput id1   : index of particle 1
cinput id2   : index of particle 2
cinput ssigtot: total cross section 
cinput sqrts : $sqrt{s}$ of collision
cinput ccolfac : scale factor for color fluctuations
c
c  This subroutine performs the scattering/annihilation
c  of two particles  or the decay
c  of one particle in the incoming channel.
c
c  Structure of this routine:
c  \begin{enumerate}
c       \item transform to NN system for proper kinematics
c       \item get collision class and number of exit-channels
c       \item loop over exit-channels and get partial cross sections
c       \item select exit-channel
c       \item save information in case of pauli-blocking
c       \item call kinematics routines
c       \item call output routine for collision file
c  \end{enumerate}
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

      include 'coms.f'
      include 'comres.f'
      include 'newpart.f'
      include 'options.f'
      include 'boxinc.f'


c local variables
      real*8 sqrts,ssigtot,sigpart,sigma(0:maxpsig),sigsum,sigfac
      real*8 e1,e2,lambda,colldens,ccolfac
      integer ind1,ind2,ityp1,ityp2,isigline,nCh,ii
      integer i,j,itot1,itot2,iso31,iso32
      integer id1,id2
chp counter for new loop
      integer ict
c     functions and subroutines
      integer collclass,isoit

c
c     make local copies of particle indices
      ind1=id1
      ind2=id2
c     increment collision-counter
      ctag=ctag+1

c     save total cross section in sigma(0)
      sigma(0)=ssigtot


c     initialize some arrays (definitions are listed in newpart.f)
      do 13 j=1,mprt
         do 12 i=1,5
            pnew(i,j)=0.0
 12      continue
         inew(j)=0
         leadfac(j)=1.d0
 13   continue

c if ind2 is less than 0, a collision with a wall takes place 
      if (ind2.lt.0) then

         call file15out(ind1,ind2,sqrts,ssigtot,sigpart)

c solid < 1 then periodic boundary conditions
            if (solid.lt.1) then
                iline = 80
                if (ind2.eq.-1) then
                        rx(ind1)=rx(ind1)-lbox
                elseif (ind2.eq.-4) then
                        rx(ind1)=rx(ind1)+lbox
                elseif (ind2.eq.-2) then
                        ry(ind1)=ry(ind1)-lbox
                elseif (ind2.eq.-5) then
                        ry(ind1)=ry(ind1)+lbox
                elseif (ind2.eq.-3) then
                        rz(ind1)=rz(ind1)-lbox
                else
                        rz(ind1)=rz(ind1)+lbox
                endif
            else
c solid walls
                iline = 81
                if (ind2.eq.-1) then
                        px(ind1)=-px(ind1)
                elseif (ind2.eq.-4) then
                        px(ind1)=-px(ind1)
                elseif (ind2.eq.-2) then
                        py(ind1)=-py(ind1)
                elseif (ind2.eq.-5) then
                        py(ind1)=-py(ind1)
                elseif (ind2.eq.-3) then
                        pz(ind1)=-pz(ind1)
                else
                        pz(ind1)=-pz(ind1)
                endif           
            Endif       

                nexit=1         
                inew(1)=ind1

        call f15outch(0.d0)

        return
        endif
cc      

c add fermi motion
      if (CTOption(30).eq.1) then
         call addfermi(ind1,peq1)
         call addfermi(ind2,peq2)
      endif

c     clear array of partial cross sections
      do 14 i=1,maxpsig
         sigma(i)=0.d0
 14   continue

c     save particle indices into pslot for later use (Danielewicz-Pratt delays)
      pslot(1)=ind1
      pslot(2)=ind2

c some abbreviations
      ityp1=ityp(ind1)
      itot1=isoit(ityp1)
      iso31=iso3(ind1)
      if(ind2.eq.0) then
         ityp2=0
         itot2=0
         iso32=0
      else
         ityp2=ityp(ind2)
         itot2=isoit(ityp2)
         iso32=iso3(ind2)
      endif




c     transform in to NN system for proper kinematics

      e1 = p0(ind1)

      if(ind2.eq.0) then 
cccc DECAY cccccccccccccccccccccccc

c     increment decay counter
         dectag=dectag+1

c     prepare output to decay file
         call file16entry(ind1)

         sqrts=fmass(ind1)
         betax=px(ind1)/e1
         betay=py(ind1)/e1
         betaz=pz(ind1)/e1
         p0nn=p0(ind1)
         pxnn=px(ind1)
         pynn=py(ind1)
         pznn=pz(ind1)
         pnn=sqrt(pxnn*pxnn+pynn*pynn+pznn*pznn)

c     set tag for decay:
         iline=20
         
      else   
cccc SCATTERING/ANNIHILATION  cccc
         e2 = p0(ind2)
c     transformation betas into two particle CMS system
         betax=(px(ind1)+px(ind2))/(e1+e2)
         betay=(py(ind1)+py(ind2))/(e1+e2)
         betaz=(pz(ind1)+pz(ind2))/(e1+e2)

c     calculate momenta in two particle CMS system
         pxnn=px(ind1)
         pynn=py(ind1)
         pznn=pz(ind1)
         p0nn=e1
c     call to Lorentz transformation
         call rotbos(0d0,0d0,-betax,-betay,-betaz,
     &        pxnn,pynn,pznn,p0nn)
         pnn=sqrt(pxnn*pxnn+pynn*pynn+pznn*pznn)

c     reduced cross section for leading hadrons of string fragmentation
c     sigfac is scaling factor for cross section
         sigfac=1.d0
         if(tform(ind1).le.acttime
     &           .and.tform(ind2).gt.acttime) then
            sigfac=xtotfac(ind2)
         else if(tform(ind2).le.acttime
     &           .and.tform(ind1).gt.acttime) then
            sigfac=xtotfac(ind1)
         elseif(tform(ind1).gt.acttime
     &        .and.tform(ind2).gt.acttime) then 
            sigfac=xtotfac(ind1)*xtotfac(ind2)
         endif

c modify sigfac due to color fluctuations
         sigfac = sigfac*ccolfac

c now get line-number for sigmaLN array in blockres
         isigline=collclass(ityp1,iso31,ityp2,iso32)

c  number of exit-channels:
         nCh=SigmaLn(1,1,isigline)

         do 10 ii=3,nCh+2 ! loop over exit-channels
c           get  process-id (iline) from sigmaLN array in blockres
            iline=SigmaLn(ii,1,isigline)
            if(iline.gt.0) then ! normal cross sections:

               call crossx(iline,sqrts,ityp1,iso31,fmass(ind1),
     &                     ityp2,iso32,fmass(ind2),sigma(ii-2))

            else   !  detailed balance:
               call crossz(iline,sqrts,ityp1,iso31,fmass(ind1),
     &                     ityp2,iso32,fmass(ind2),sigma(ii-2))
            endif  ! end of detailed balance
c     ensure reduction of part. cross sections within formation time
            sigma(ii-2)=sigfac*sigma(ii-2)


 10      continue               ! end of exit channel loop
   
chp label if elastic scattering is disabled        
         ict=0 
 101     continue
chp
c     unitarize partial cross sections (rescale sum to match the total cross section)
         call normit (sigma,isigline)

c     select partial cross section ii
         call getbran(sigma,0,maxpsig,sigsum,1,nCh,ii)

         if(sigsum.lt.1d-10)then 
c           write(6,*)'***(W) scatter: ',
c     ,     'no entry found for fin. state -> forced elastic scat. ',
c     ,     'particles:line=',isigline,'ecm=',sqrts
c           write(6,*)'     collision of    : ',
c     &                        ityp1,iso31,fmass(ind1),
c     &                        ityp2,iso32,fmass(ind2)
c           write(6,*)'                sigma=',sigma,ii
c     force elastic scattering
           sigpart=sigma(0)
           iline=13
        else
           sigpart=sigma(ii)
           iline=SigmaLN(ii+2,1,isigline)
        end if
chp no elastic scattering of initial pp pn and nn
       ict=ict+1 
       if(CTOption(7).eq.1.and.(iline.eq.17.or.iline.eq.19))then
        if(ict.le.10)then
          goto 101
        else 
c         write(6,*) 'in scatter.f: too many tries',
c     &    'elastic scattering happens despite cto 7 1'  
        end if
       end if   
chp      
 
      endif     
c.. take care sigpart in case of decay (no sigpart defined)
      if (iline.eq.20) sigpart=0d0

cccc end of scatter/decay if
      call file15out(ind1,ind2,sqrts,ssigtot,sigpart)

c     save old particle information in case of pauli blocking
c     or rejection due to energy non-conservation
      call saveinfo(ind1,1)
      call saveinfo(ind2,2)
      
c... prepare exit channel
      call scatprep(ind1,ind2,sqrts,sigpart)

      lambda=1.0d0
      call prescale(lambda)

      call scatfinal(colldens)

c     output to collision file
      call f15outch(colldens)
      call osc99_coll

c     write output to decay file
      if(ind2.eq.0) call file16write
c     in case of ctoption(13).ne.0 more output is written with the next call
      call f16outch

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine scatprep(in1,in2,sqrts,sigpart)
c
c     Revision : 1.0
c
cinput ind1   : index of ingoing particle 1 
cinput ind2   : index of ingoing particle 2
cinput sqrts  :  $sqrt{s}$ of collision
cinput sigpart: partial cross section for process
c
c {\tt scatprep, prescale} and {\tt scatfinal} handle the collision/decay
c kinematics and the book-keeping of the global particle vectors.
c In {\tt scatprep} the exit channel is generated via a call to
c {\tt make22}, the information is stored in common blocks of the 
c {\tt newpart.f} file. Furthermore the new particle order including
c new slots for the exit channel are generated in {\tt scatprep}
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'coms.f'
      include 'options.f'
      include 'comres.f'
      include 'newpart.f'

      integer i,ind1,ind2,k,in1,in2
      integer itmp(mprt),ipmp(mprt)
      real*8 sqrts,sigpart
      real*8 phi1,phi2
      real*8 pzi1,pzi2,pxi1,pxi2,pyi1,pyi2
      real*8 rstringx(2),rstringy(2),rstringz(2),tstring(2),tformold(2)
      real*8 rpott(2),rpotx(2),rpoty(2),rpotz(2)
      real*8 th,rdum, ctheta1
      integer bar,nb1,nb2,nm1,nm2,meslist2(mprt)
      integer barlist1(mprt),barlist2(mprt),meslist1(mprt)
      logical dnscall

c     functions:

      real*8 pcms
      common /scatcomr/rstringx,rstringy,rstringz,tstring,
     &                rpott,rpotx,rpoty,rpotz,
     &                pzi1,pzi2,pxi1,pxi2,pyi1,pyi2,
     &                ctheta1,phi1,th,phi2,tformold
      common /scatcomi/itmp,ipmp,ind1,ind2,nb1,nb2,
     &                 bar,nm1,nm2
      common /scatcoml/dnscall

c     call to paulibl for baryon-density
      dnscall=CTOption(39).ne.0


c     reset some pointer arrays
      do 10 i=1,mprt
         itmp(i)=0
         ipmp(i)=0
         barlist1(i)=0
         barlist2(i)=0
         meslist1(i)=0
 10   continue

      ind1=in1
      ind2=in2

c     reset xtotfacs
      if(tform(ind1).le.acttime) xtotfac(ind1)=1.d0
      if(ind2.gt.0) then
        if(tform(ind2).le.acttime) xtotfac(ind2)=1.d0
      endif

c     1. determine mass, ityp's, iso3's of outgoing channel
c
      if(ind2.ne.0) then ! scattering/annihilation
         call make22(iline,sqrts,
     &               ityp(ind1),iso3(ind1),fmass(ind1),xtotfac(ind1),
     &               ityp(ind2),iso3(ind2),fmass(ind2),xtotfac(ind2))


      elseif(ind2.eq.0) then ! decay
         call make22(iline,sqrts,
     &               ityp(ind1),iso3(ind1),fmass(ind1),xtotfac(ind1),
     &               0,0,0d0,1d0)
      endif

c set here the true nexit (number of particles in the exit channel): 
      nexit=nstring1+nstring2

        
c
c     2a). annihilation/soft resonance production: store properties
c
      if(nstring2.eq.0) then
c     store average  time and position of particles (only one part/string coming out...)
         tstring(1)=(r0(ind1)+r0(ind2))/2. 
         rstringx(1)=(rx(ind1)+rx(ind2))/2.
         rstringy(1)=(ry(ind1)+ry(ind2))/2.
         rstringz(1)=(rz(ind1)+rz(ind2))/2.
c     do the same for the MD trajectories
         rpott(1)=(r0_t(ind1)+r0_t(ind2))/2.
         rpotx(1)=(rx_t(ind1)+rx_t(ind2))/2.
         rpoty(1)=(ry_t(ind1)+ry_t(ind2))/2.
         rpotz(1)=(rz_t(ind1)+rz_t(ind2))/2.

         pnnout=pnn


c     store the incoming momenta for further use before the slots get erased/replaced
c     by outgoing particles
         pxi1=px(ind1)
         pxi2=px(ind2)
         pyi1=py(ind1)
         pyi2=py(ind2)
         pzi1=pz(ind1)
         pzi2=pz(ind2)

c     store old formation times
         tformold(1)=max(tform(ind1),tform(ind2))
         tformold(2)=max(tform(ind1),tform(ind2))
c two leading hadrons form s-channel -> the larger reduction factor
c of the incoming hadrons is maintained:
         xtotfacold(1)=max(xtotfac(ind1),xtotfac(ind2))
         xtotfacold(2)=max(xtotfac(ind1),xtotfac(ind2))
      else

c     2b) 1. scattering/decay/strings: store both locations for later use
c        store time and position of both incoming particles
         tstring(1)=r0(ind1)
         rstringx(1)=rx(ind1)
         rstringy(1)=ry(ind1)
         rstringz(1)=rz(ind1)
c        likewise for the MD trajectory arrays
         rpott(1)=r0_t(ind1)
         rpotx(1)=rx_t(ind1)
         rpoty(1)=ry_t(ind1)
         rpotz(1)=rz_t(ind1)

c        store old formation time and reduction factor
         tformold(1)=tform(ind1)
         xtotfacold(1)=xtotfac(ind1)

         if(ind2.eq.0) then !decay
c           likewise for particles 2,3,4...
            tstring(2)=r0(ind1)
            rstringx(2)=rx(ind1)
            rstringy(2)=ry(ind1)
            rstringz(2)=rz(ind1)
            rpott(2)=r0_t(ind1)
            rpotx(2)=rx_t(ind1)
            rpoty(2)=ry_t(ind1)
            rpotz(2)=rz_t(ind1)

            tformold(2)=tform(ind1)
            xtotfacold(2)=xtotfac(ind1)

         else
            tstring(2)=r0(ind2)
            rstringx(2)=rx(ind2)
            rstringy(2)=ry(ind2)
            rstringz(2)=rz(ind2)

            rpott(2)=r0_t(ind2)
            rpotx(2)=rx_t(ind2)
            rpoty(2)=ry_t(ind2)
            rpotz(2)=rz_t(ind2)

            tformold(2)=tform(ind2)
            xtotfacold(2)=xtotfac(ind2)

         endif
c
c     2b) 2. get relative momentum in 2 particle cms (only for ang. distrib.!)
c
         if(iline.gt.0) then
            pnnout=pcms(sqrts,mstring(1),mstring(2))
         endif

      endif    ! nstring2.eq.0
c
c     2b) 3.  get angular distribution
c

      if (ind2.ne.0) then
       call angdisnew(sqrts,fmass(ind1),fmass(ind2),iline,ctheta1,phi1)
      else 
       call angdisnew(sqrts,fmass(ind1),0.0d0,iline,ctheta1,phi1)
      end if  


c   for CTOption(2)<>0 phi1=0 and 2-part. scattering is in the scattering plane
      if(ind2.ne.0.and.CTOption(2).ne.0) then
         phi1=0.d0
      endif

c     get angles between two particle CMS and computational frame
      call getang(pxnn,pynn,pznn,th,phi2,rdum)

c     option for pure forward streaming
      if(CTOption(14).ne.0)then 
         ctheta1=1d0
         th=0d0
      end if
c
c     3. generate slots for new particles
c
c     count baryons AND anti-baryons in old slots (NOT net-baryons!!)
      bar=0
      if(abs(ityp(ind1)).le.maxbar)then
         bar=bar+1
      endif
c
cccc in case of baryon-boson scattering, inew(1) contains the baryon,
cccc otherwise it does not matter
c
      if(ind2.ne.0) then
         if(abs(ityp(ind2)).le.maxbar) then
            bar=bar+1
         endif
         if(abs(ityp(ind2)).le.maxbar.and.
     &        abs(ityp(ind1)).ge.minmes) then
            inew(1)=ind2
            inew(2)=ind1
         else
            inew(1)=ind1
            inew(2)=ind2
         endif
      else
cccc decay: must create a second slot
         inew(1)=ind1
c     increment meson counter
         nmes=nmes+1
c     set new slot number
         inew(2)=nbar+nmes
c     create the slot
         call addpart(inew(2))
c     this is needed for delpart (to esure correct counters):
         ityp(inew(2))=999
c     three or four body decays
         if(nexit.gt.2) then
            do 91 i=3,nexit
               nmes=nmes+1
               inew(i)=nbar+nmes
               call addpart(inew(i))
               ityp(inew(i))=999
 91         continue
         endif
      endif

cccc soft resonance production (second slot must be deleted)
      if(nexit.eq.1) then
c        the smaller index automatically belongs to a baryon if there was one
         inew(1)=min0(ind1,ind2)
c        get index to be deleted 
         inew(2)=max0(ind1,ind2)
c        delete second slot
         call delpart(inew(2))
      endif

c     generate sorting-order for new particles: 
c      1. leading baryon of particle/string 1
c      2. leading baryon of particle/string 2
c      3. all other baryons of particle/string 1
c      4. all other baryons of particle/string 2
c      5. all mesons of particle/string 1
c      6. all mesons of particle/string 2
c     the itypnew-indices are stored into the itmp-array
c     the location-index (do the new particles "belong" to incident particle 1 or 2) 
c      is stored in the ipmp-array

      nb1=0  ! number of baryons in string 1
      nb2=0  ! number of baryons in string 2
      nm1=0  ! number of mesons in string 1
      nm2=0  ! number of mesons in string 2
c     now get the number of baryons and mesons and their ityps from string/particle 1
      call instring(1,nstring1,nb1,nm1,barlist1,meslist1)
      if(nstring2.ne.0) 
c     likewise from string/particle 2
     &     call instring(nstring1+1,nexit,nb2,nm2,barlist2,meslist2)

c now really generate the above described sorting-order
      if(  (bar.eq.0)
     & .or.(bar.eq.1.and.nb1.gt.0)
     & .or.(bar.eq.2.and.nb1.gt.0.and.nb2.eq.0)) then
c     loop over baryons in string/particle 1
         do 381 i=1,nb1
            itmp(i)=barlist1(i)
            ipmp(i)=1
 381     continue
c     loop over baryons in string/particle 2
         do 382 i=1,nb2
            itmp(nb1+i)=barlist2(i)
            ipmp(nb1+i)=2
 382     continue
      elseif((bar.eq.1.and.nb2.gt.0.and.nb1.eq.0)
     &   .or.(bar.eq.2.and.nb2.gt.0.and.nb1.eq.0)) then
c     loop over baryons in string/particle 2
         do 391 i=1,nb2
            itmp(i)=barlist2(i)
            ipmp(i)=2
 391     continue
      elseif(bar.eq.2.and.nb1.gt.0.and.nb2.gt.0) then
         itmp(1)=barlist1(1)
         ipmp(1)=1
         itmp(2)=barlist2(1)
         ipmp(2)=2
c     loop over baryons in string/particle 1
         do 371 i=2,nb1
            itmp(1+i)=barlist1(i)
            ipmp(1+i)=1
 371     continue
c     loop over baryons in string/particle 2
         do 372 i=2,nb2
            itmp(nb1+i)=barlist2(i)
            ipmp(nb1+i)=2
 372     continue
      endif
c     loop over mesons in string/particle 1
      do 383 i=1,nm1
         itmp(nb1+nb2+i)=meslist1(i)
         ipmp(nb1+nb2+i)=1
 383  continue
c     loop over mesons in string/particle 2
      do 384 i=1,nm2
         itmp(nb1+nb2+nm1+i)=meslist2(i)
         ipmp(nb1+nb2+nm1+i)=2
 384  continue


c in case of annihilation without baryon production
c or if a baryon-antibaryon pair is created via a meson-meson collision,
c delete old slots
c
c note: this sequence is needed to keep all baryons at the top
c       of the particle table. However, it destroys distinctive
c       projectile/target areas in the particle arrays.
c

      if((bar.eq.2.and.(nb1+nb2).eq.0)
     &   .or.
     &   (bar.eq.0.and.(nb1+nb2).gt.0)
     &  ) then

c     store average  time and position of particles 
         tstring(1)=(r0(ind1)+r0(ind2))/2.
         rstringx(1)=(rx(ind1)+rx(ind2))/2.
         rstringy(1)=(ry(ind1)+ry(ind2))/2.
         rstringz(1)=(rz(ind1)+rz(ind2))/2.

         tstring(2)=tstring(1) 
         rstringx(2)=rstringx(1)
         rstringy(2)=rstringy(1)
         rstringz(2)=rstringz(1)

c     do likewise for MD trajectories
         rpott(1)=(r0_t(ind1)+r0_t(ind2))/2.
         rpotx(1)=(rx_t(ind1)+rx_t(ind2))/2.
         rpoty(1)=(ry_t(ind1)+ry_t(ind2))/2.
         rpotz(1)=(rz_t(ind1)+rz_t(ind2))/2.
         rpott(2)=rpott(1)
         rpotx(2)=rpotx(1)
         rpoty(2)=rpoty(1)
         rpotz(2)=rpotz(1)

         call delpart(inew(1))
         call delpart(inew(2))
         inew(1)=0
         inew(2)=0
c in case of a meson-string with baryon-antibaryon creation
c the meson-slot must be deleted
      elseif(bar.eq.1.and.(nb1+nb2).gt.1) then
         call delpart(inew(2))
         inew(2)=0
      endif

c now create new slots:
      do 307 i=1,nexit 
         if(inew(i).lt.1) then
c     the new particle ID is stored in itypnew
            if(iabs(itypnew(itmp(i))).le.maxbar) then
c     the particle is a baryon
               do 385 k=1,i
c     make sure that mesons in the inew array are shifted upwards
                  if(inew(k).gt.nbar)inew(k)=inew(k)+1
 385           continue
c     increment baryon counter
               nbar=nbar+1
c     this is the new particle slot
               inew(i)=nbar
            elseif(iabs(itypnew(itmp(i))).ge.minmes) then
c     the particle is a meson
               nmes=nmes+1
               inew(i)=nbar+nmes
            endif
c     create the slot
            call addpart(inew(i))
         endif
 307  continue

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prescale(lambda)
c
c     Revision : 1.0
c
cinput lambda : scaling factor for momenta of outgoing particles
c
c {\tt scatprep, prescale} and {\tt scatfinal} handle the collision/decay
c kinematics and the book-keeping of the global particle vectors.
c In {\tt prescale} the actual kinematics for the exit-channel, including
c a call to {\tt angdis} for the angular distribution and the transformation
c from the two particle restframe to the computational frame is being 
c performed. Momenta and locations of exit channel particles are written
c to the global particle vectors. \\
c The scaling factor {\tt lambda} is needed to ensure global energy conservation
c in the case of momentum dependent interactions (MDI). In that case {\tt prescale}
c can be called several succesive times with different values of {\tt lambda}
c to minimize $E_{tot,in} - E_{tot,out}$.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'coms.f'
      include 'options.f'
      include 'newpart.f'

      integer i,j,itmp(mprt),ipmp(mprt),ind1,ind2,nb1,nb2
      real*8 lambda,th,phi2,ctheta1,phi1,theta3,phi3,pabs
      real*8 pzi1,pzi2,pxi1,pxi2,pyi1,pyi2
      real*8 rstringx(2),rstringy(2),rstringz(2),tstring(2)
      real*8 rpott(2),rpotx(2),rpoty(2),rpotz(2)
      real*8 tauf(mprt),tformold(2)
      integer getspin,bar,nm1,nm2

      common /scatcomr/rstringx,rstringy,rstringz,tstring,
     &                 rpott,rpotx,rpoty,rpotz,
     &                 pzi1,pzi2,pxi1,pxi2,pyi1,pyi2,
     &                 ctheta1,phi1,th,phi2,tformold
      common /scatcomi/itmp,ipmp,ind1,ind2,nb1,nb2,
     &                 bar,nm1,nm2


      if(nexit.eq.1) then ! soft resonance production

c new momenta are stored in pnew
         pnew(1,1)=lambda*(pxi1+pxi2)
         pnew(2,1)=lambda*(pyi1+pyi2)
         pnew(3,1)=lambda*(pzi1+pzi2)
c this is the energy
         pnew(4,1)=dsqrt(pnew(5,1)**2+
     &             pnew(1,1)**2+pnew(2,1)**2+pnew(3,1)**2)
c formation time
         tauf(1)=0.d0         
      else                      ! scattering/decay
         do 205 j=1,nexit
c compute formation time (as a eigentime)
            tauf(j)=xnew(4,j)*pnew(5,j)/pnew(4,j)

c     rescale momenta of particles
            call getang(pnew(1,j),pnew(2,j),pnew(3,j),theta3,phi3,pabs)
            pabs=lambda*pabs
            call putang(pnew(1,j),pnew(2,j),pnew(3,j),theta3,phi3,pabs)


c check for forward time-delay
c in case of delay the momenta are already in the comp. frame
            if(.not.(CTOption(34).eq.2.and.iline.eq.20.and.
     &         ityptd(1,pslot(1)).ne.0)) then

c     rotate in the from z-axis to th&phi given by angdis 
               call rotbos(dacos(ctheta1),phi1,0d0,0d0,0d0,
     &                     pnew(1,j),pnew(2,j),pnew(3,j),pnew(4,j))

c     rotate from the NN to the comp. sys. and
c     transform particles to computational system

c     decays should not be rotated
               if (iline.eq.20) then 
                  th = 0.d0
                  phi2 = 0.d0
               end if

               call rotbos(th,phi2,betax,betay,betaz,
     &                     pnew(1,j),pnew(2,j),pnew(3,j),pnew(4,j))

            endif
c end of delay-if


 205     continue
      endif

c     write coordinates for energy-check and pauli-blocking
c     correct particle/string-locations from relative to absolute

      do 215 i=1,nexit
c     write locations to global arrays
c     the ipmp values determine wether the new particle belongs to incoming slot 1 or 2
         r0(inew(i))=tstring(ipmp(i))
         rx(inew(i))=rstringx(ipmp(i))
         ry(inew(i))=rstringy(ipmp(i))
         rz(inew(i))=rstringz(ipmp(i))
cpot
c     likewise for MD trajectories
         r0_t(inew(i))=rpott(ipmp(i))
         rx_t(inew(i))=rpotx(ipmp(i))
         ry_t(inew(i))=rpoty(ipmp(i))
         rz_t(inew(i))=rpotz(ipmp(i))

c     write momenta to global arrays
         p0(inew(i))=pnew(4,itmp(i))
         px(inew(i))=pnew(1,itmp(i))
         py(inew(i))=pnew(2,itmp(i))
         pz(inew(i))=pnew(3,itmp(i))


c store formation time and leading hadron 
         if(pnew(5,itmp(i)).gt.0d0)then
            tform(inew(i))=tstring(1)+tauf(itmp(i))*
     &                  pnew(4,itmp(i))/pnew(5,itmp(i))
         else
           tform(inew(i))=tstring(ipmp(i))
         endif

c     cross section reduction factor and string ID
         xtotfac(inew(i))=leadfac(itmp(i))

c additional reduction of the leading hadrons' cross section
c in case the incoming hadron x-sec was already reduced
c otherwise an elastic scattering would erase formation time...
         if(xtotfacold(ipmp(i)).lt.1d0) then
            xtotfac(inew(i))=xtotfacold(ipmp(i))*xtotfac(inew(i))
         endif
         if(tform(inew(i)).lt.tformold(ipmp(i))
     &      .and.xtotfac(inew(i)).gt.0) then
            tform(inew(i))=tformold(ipmp(i))
         endif


c     write mass, ID, I3 and spin to global arrays
         fmass(inew(i))=pnew(5,itmp(i))
         ityp(inew(i))=itypnew(itmp(i))
         iso3(inew(i))=i3new(itmp(i))
         spin(inew(i))=getspin(itypnew(itmp(i)),-1)

 215  continue

c set lstcoll:
c     lstcoll relates the outgoing particles of this scattering/decay interaction
c     to each other via the number of the collision - it is used to prevent them
c     from directly interacting again with each other because this would be unphysical.
      do i=1,nexit
         lstcoll(inew(i))=ctag
      end do
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine scatFinal(colldens)
c
c     Revision : 1.0
c
coutput colldens : baryon density at point of interaction
c
c {\tt scatprep, prescale} and {\tt scatfinal} handle the collision/decay
c kinematics and the book-keeping of the global particle vectors.
c In {\tt scatfinal} the interaction is checked for Pauli-blocking and all
c particle {\em quantum numbers} which so far have not been set are written
c to the global arrays. Some collision counters are incremented here, too.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'coms.f'
      include 'comres.f'
      include 'newpart.f'
      include 'freezeout.f'
      include 'options.f'

      integer i,ind1,ind2,itmp(mprt),ipmp(mprt),nb1,nb2,j
      integer ikill1,ikill2,bar,nm1,nm2
      integer iloffset
      real*8 colldens
      logical dnscall,dnsdum

c     functions:
      integer fchg
      real*8 dectim
      logical paulibl

      common /scatcomi/itmp,ipmp,ind1,ind2,nb1,nb2,
     &                 bar,nm1,nm2
      common /scatcoml/dnscall

c     check pauli-blocking (only for nucleons and delta(1232)
c     loop over all baryons in exit-channel
      do 306 i=1,nb1+nb2
         if(itypnew(itmp(i)).eq.nucleon) then
            dnscall=.false.
            ikill1=0
            ikill2=0
c     call to Pauli-Blocker
            if(paulibl(inew(i),colldens,-1000)) then
c     pauli-blocked collision, restore information and return

c     first delete unnecessary slots
               if(nexit.gt.2) then
                  do 317 j=3,nexit
                     call delpart(inew(j))
 317              continue
               endif

c in case of blocked decay, delete second entry
               if(ind2.eq.0) then
                  call delpart(inew(2))
c only delete/create slots if something more than a mere
c swapping has happened
               elseif(.not.(((inew(1).eq.ind1).or.(inew(1).eq.ind2))
     &                .and.((inew(2).eq.ind1).or.(inew(2).eq.ind2)))) 
     &                 then
                  ikill1=max0(inew(1),inew(2))
                  ikill2=min0(inew(1),inew(2))
                  call delpart(ikill1)          
                  if(ikill2.gt.0) call delpart(ikill2)
cgg make sure slot with the lower slotnumber is created first
                  if(ind2.gt.0) then
                    call addpart(min0(ind1,ind2))
                    call addpart(max0(ind1,ind2))
                  else 
                    call addpart(ind1)
                  endif
               endif
c     restore old contents of slots ind1 and ind2
               call saveinfo(ind1,-1)
               call saveinfo(ind2,-2)
c     in case of blocked decay, sample new decay time
               if(ind2.eq.0) then
                  dectime(ind1)=dectim(ind1,iline)+acttime
               endif
c     restore frozen fermi state if necessary
               if (CTOption(30).eq.1) then
                  call savefermi(ind1,ind1,peq1)
                  call savefermi(ind2,ind2,peq2)
               endif
c     set lstcoll
               if(ind2.ne.0) then
                  lstcoll(ind1)=ctag
                  lstcoll(ind2)=ctag
               else
                  lstcoll(ind1)=ctag
               endif
               
c     set baryon and meson counters 
               if(ikill1.ne.0) then
                  if(abs(ityp(ind1)).lt.minmes) then
                     nbar=nbar+1
                  else
                     nmes=nmes+1
                  endif
               endif
               if(ind2.gt.0.and.ikill2.ne.0) then
                  if(abs(ityp(ind2)).lt.minmes) then
                     nbar=nbar+1
                  else
                     nmes=nmes+1
                  endif
               endif

c     increment blocking counter
               NBlColl=NBlColl+1
c     set number of particles in exit-channel to zero
               nexit=0
               return
            endif
         endif
 306  continue

c     the collision was not blocked, proceed...

c     4. rewrite arrays, clear entries etc.
c     fill the rest of the slot

c     in case the Pauli-Blocker was not called, do it now in order to 
c     get the baryon density at the interaction point, colldens
      if(dnscall) dnsdum=paulibl(inew(1),colldens,-1000)

c     write rest of the quantum numbers to global vectors
      do 312 i=1,nexit
c     number of collisions
         ncoll(inew(i))=ncoll(inew(i))+1
c     charge
         charge(inew(i))=fchg(iso3(inew(i)),ityp(inew(i)))
c     time of decay
         dectime(inew(i))=dectim(inew(i),iline)+tform(inew(i))
c     last process the particle was involved in
         if ( (iline.ne.13).and.(iline.ne.17).and.(iline.ne.19).and.
     $        (iline.ne.22).and.(iline.ne.26).and.(iline.ne.38)) then
            if (ind2.gt.0) then
               iloffset=0
               if ( (iline.eq.27).or.(iline.eq.28) ) then
                  if ( (CTOption(41).gt.1) ) then
                     if ( (abs(itypt(1)).gt.100).and.
     $                    (abs(itypt(2)).gt.100) ) then
                        iloffset=20
                     endif
                     if ( (itypt(1).eq.ityp(inew(i))).and.
     $                    (iso3t(1).eq.iso3(inew(i))) ) then
                        origin(inew(i))=origint(1)+100
                        uid(inew(i))=uidt(1)
                     elseif ( (itypt(2).eq.ityp(inew(i))).and.
     $                       (iso3t(2).eq.iso3(inew(i))) ) then
                        origin(inew(i))=origint(2)+100
                        uid(inew(i))=uidt(2)
                     else
                        origin(inew(i))=iline+iloffset
     $                       +1000*(iabs(itypt(1))+1000*iabs(itypt(2)))
                     endif
                  else
                     origin(inew(i))=iline+iloffset
     $                    +1000*(iabs(itypt(1))+1000*iabs(itypt(2)))
                  endif
               else
                  origin(inew(i))=iline
     $                 +1000*(iabs(itypt(1))+1000*iabs(itypt(2)))
               endif
            else
               origin(inew(i))=iline
     $              +1000*(iabs(itypt(1)))
            endif
c     unique ID tag
            uid_cnt=uid_cnt+1
            uid(inew(i))=uid_cnt
         else
            origin(inew(i))=origint(i)+100
            uid(inew(i))=uidt(i)
            if (ityp(inew(1)).ne.itypt(1)) then
               if (nexit.ne.2) then
                  write(6,*) "fatal error in scatter: nexit.ne.2!"
                  stop 137
               endif
               origin(inew(i))=origint(3-i)+100
               uid(inew(i))=uidt(3-i)
            endif
         endif
c     freeze-out coordinates
c     freeze-out coordinates
chp modify freeze-out coordinates
chp particles should not freeze out before they are formed
chp add tform and propagate accordingly
chp not used when visualization output is wanted, then freeze-out 
chp coordinates are production points
        if(tform(inew(i)).gt.r0(inew(i)).and.CTOption(55).eq.0)then
         frr0(inew(i))=tform(inew(i))
         frrx(inew(i))=rx(inew(i))
     & +px(inew(i))/p0(inew(i))*(tform(inew(i))-r0(inew(i)))
         frry(inew(i))=ry(inew(i))
     & +py(inew(i))/p0(inew(i))*(tform(inew(i))-r0(inew(i)))
         frrz(inew(i))=rz(inew(i))
     & +pz(inew(i))/p0(inew(i))*(tform(inew(i))-r0(inew(i)))
chp otherwise just proceed normally
        else
         frr0(inew(i))=r0(inew(i))
         frrx(inew(i))=rx(inew(i))
         frry(inew(i))=ry(inew(i))
         frrz(inew(i))=rz(inew(i))
        end if
         frp0(inew(i))=p0(inew(i))
         frpx(inew(i))=px(inew(i))
         frpy(inew(i))=py(inew(i))
         frpz(inew(i))=pz(inew(i))

c forward time-delay
      if(CTOption(34).eq.2.and.(iline.eq.36.or.iline.eq.37)) then
         do 307 j=1,2
            p0td(j,inew(i))=pold(4,j)
            pxtd(j,inew(i))=pold(1,j)
            pytd(j,inew(i))=pold(2,j)
            pztd(j,inew(i))=pold(3,j)
            fmasstd(j,inew(i))=pold(5,j)
            ityptd(j,inew(i))=itypold(j)
            iso3td(j,inew(i))=iso3old(j)
 307     continue
      else
         do 308 j=1,2
            p0td(j,inew(i))=0.d0
            pxtd(j,inew(i))=0.d0
            pytd(j,inew(i))=0.d0
            pztd(j,inew(i))=0.d0
            fmasstd(j,inew(i))=0.d0
            ityptd(j,inew(i))=0
            iso3td(j,inew(i))=0
 308     continue
      endif
               
c     in a N particle exit channel the produced resonance 
c     either comes from BB->B* ? ? or a String -> B* ? ? 
         if(ityp(inew(i)).ge.minres.and.ityp(inew(i)).le.maxres) then
            if(iline.ne.20.or.iline.ne.10) then
c     B B -> ? ? 
               NHardRes=NHardRes+1
            elseif(iline.eq.20) then
               NDecRes=NDecRes+1
            elseif(iline.eq.10) then
               NSoftRes=NSoftRes+1   
            endif
         endif
 312  continue
c other counters for B-strings/M-strings could be added here...
      if(iline.eq.13.or.iline.eq.17
     &     .or.iline.eq.19.or.iline.eq.26) then
         NElColl=NElColl+1
      endif

      return 
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function collclass(ityp1,iso31,ityp2,iso32)
c
c     Revision : 1.0
c
c     This function links the ingoing collision channel to the 
c     line number of the sigmaLN array in {\tt blockres.f} in which the
c     types of cross section parametrizations to perform are defined.
c
cinput ityp1  : ityp of particle 1
cinput iso31  : $I_3$ of particle 1
cinput ityp2  : ityp of particle 2
cinput iso32  : $I_3$ of particle 2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'comres.f'
      include 'options.f'
      integer ityp1,ityp2,iso31,iso32,i1,i2,iz1,iz2
      real*8 d1,d2

c     copy to new variables MUST be done because of swpizm
c     (otherwise scatter will crash due to swapped particle properties)
      i1=iabs(ityp1)
      i2=iabs(ityp2)
      iz1=iso31
      iz2=iso32
      d1=0.d0
      d2=0.d0


c... blank out PYHTIA pdg-id's
c... define collclasses for specific reactions if needed
      if(i1.gt.1000 .or. i2.gt.1000) then
         collclass=0
         return
      endif


cJS    No rescattering of photons

      if((i1.eq.minmes).or.(i2.eq.minmes))then
         collclass=0
         return
      endif

              
      if(i1.lt.i2)call swpizm(i1,iz1,d1,i2,iz2,d2)

c baryon-antibaryon
      if(i1.le.maxbar.and.i2.le.maxbar.and.ityp1*ityp2.lt.0) then
         collclass=11
         if(CTOption(19).ne.0)collclass=0
         return
      end if 

c nucleon nucleon
      if(i1.eq.nucleon.and.i2.eq.nucleon) then
         if(iz1.eq.iz2) then
c proton proton or neutron-neutron
            collclass=2
            return
         else
c proton neutron
            collclass=1
            return
         endif
      elseif(i1.eq.mindel.and.i2.eq.nucleon) then
c Delta(1232) nucleon
         collclass=3
         return
      elseif(i1.gt.minnuc.and.i1.le.maxnuc.and.
     &       i2.eq.nucleon) then
c N+ nucleon
         collclass=4
         return
      elseif(i1.ge.mindel.and.
     &       i1.le.maxdel.and.
     &       i2.eq.nucleon) then
c D* nucleon
         collclass=5
         return
      elseif(i1.eq.mindel.and.i2.eq.mindel) then
c Delta(1232)-Delta(1232) 
         collclass=6
         return
      elseif(i1.eq.mindel.and.
     &       i2.gt.minnuc.and.
     &       i2.le.maxnuc) then
c Delta(1232)-N*
         collclass=7
         return
      elseif(i2.eq.mindel.and.
     &       i1.gt.mindel.and.
     &       i1.le.maxdel) then
c Delta(1232)-D*
         collclass=8
         return
      elseif(i1.ge.minmes.and.i1.lt.133.and.i2.le.maxbar) then
c Boson-Baryon
         collclass=9
         return
      elseif(i1.ge.minmes.and.i1.lt.133
     &        .and.i2.ge.minmes.and.i2.lt.133) then
c Boson-Boson
         collclass=10
         return
      elseif(i1.gt.minnuc.and.
     &       i1.le.maxdel.and.
     &       i2.gt.minnuc.and.
     &       i2.le.maxdel) then
c D*-D* or D*-N* or N*-N*
         collclass=12
         return
      elseif(i1.le.maxbar.and.i2.le.maxbar) then
c all remaining BB-collisions
         collclass=13
         return
      elseif(i1.ge.133.and.i2.ge.minmes) then
c M_c - Meson
         collclass=14
      elseif(i1.ge.133.and.i2.le.maxbar) then
c M_c - Baryon
         collclass=15
         return
c no charmonium-charmonium collisions
c (should not happen, see conditions for M_c - Meson scattering)
      elseif(i1.ge.135.and.i1.le.137.and.i2.ge.135.and.i2.le.137)then
         collclass=0
      else
c class not defined (sets sigtot to zero)
      write(*,*)'scatter: collclass of ',i1,i2,' not yet defined!'
         collclass=0
      endif
      return
      end


C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine getang(x,y,z,th,ph,r)
c
c gives spherical coordinates of cartesian 3-vector $(x,y,z)$
c
c input : 3-vector x,y,z 
c output: angles {\tt th}($\vartheta$), {\tt ph}($\varphi$) and radius {\tt r}
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      real*8 x,y,z,th,ph,cut,r
      parameter (cut=1d-9)
      if(abs(x).lt.cut.and.abs(y).lt.cut) then
         ph=0d0
      else
         ph=datan2(y,x)
      endif
      r=sqrt(x*x+y*y+z*z)
      th=dacos(z/max(r,cut))
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine putang(x,y,z,th,ph,r)
c
c  creates 3-vector $(x,y,z)$ out of  spherical coordinates 
c  input: angles {\tt th}($\vartheta$), {\tt ph}($\varphi$) and radius {\tt r}
c  output: 3-vector x,y,z 
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      real*8 x,y,z,th,ph,r
      x=r*sin(th)*cos(ph)
      y=r*sin(th)*sin(ph)
      z=r*cos(th)
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      SUBROUTINE rotbos(THE,PHI,BEX,BEY,BEZ,p1,p2,p3,p4)
c
c INPUT: the,phi,bex,bey,bez,p  
c OUTPUT: p
c 1)rotate 4-vector p according to the and phi 2/ boost 4-vector p
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      real*8 P(4),BEX,BEY,BEZ,GA,BEP,GABEP,rot(3,3),the,phi
      real*8 p1,p2,p3,p4,bb2
      integer j


      IF(THE**2+PHI**2.GT.1E-20) THEN
C...ROTATE 
        ROT(1,1)=COS(THE)*COS(PHI)
        ROT(1,2)=-SIN(PHI)
        ROT(1,3)=SIN(THE)*COS(PHI)
        ROT(2,1)=COS(THE)*SIN(PHI)
        ROT(2,2)=COS(PHI)
        ROT(2,3)=SIN(THE)*SIN(PHI)
        ROT(3,1)=-SIN(THE)
        ROT(3,2)=0.
        ROT(3,3)=COS(THE)
        DO 108 J=1,3
 108       P(J)=ROT(J,1)*P1+ROT(J,2)*P2+ROT(J,3)*P3
        p(4)=p4
      else
        p(1)=p1
        p(2)=p2
        p(3)=p3
        p(4)=p4       
      ENDIF

      bb2=BEX**2+BEY**2+BEZ**2
      IF(bb2.GT.1E-20) THEN
C...LORENTZ BOOST (TYPICALLY FROM REST TO MOMENTUM/ENERGY=BETA)
        GA=1D0/DSQRT(1D0-bb2)
        BEP=BEX*P(1)+BEY*P(2)+BEZ*P(3)
        GABEP=GA*(GA*BEP/(1D0+GA)+P(4))
        P(1)=P(1)+GABEP*BEX
        P(2)=P(2)+GABEP*BEY
        P(3)=P(3)+GABEP*BEZ
        P(4)=GA*(P(4)+BEP)
      ENDIF

      p1=p(1)
      p2=p(2)
      p3=p(3)
      p4=p(4)

      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine instring(low,high,nb,nm,barlist,meslist)
c
c     version: 1.0
c
c     this subroutine returns arrays with the indices of the baryons/mesons
c     found in the itypnew-array of {\tt newpart.f} 
c
cinput low:  lower search boundary
cinput high: upper search boundary
c
coutput nb : number of baryons in range
coutput nm : number of mesons in range
coutput barlist: list of indices for baryons
coutput meslist: list of indices for mesons
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'newpart.f'
      include 'comres.f'

      integer i,low,high,nb,nm,barlist(mprt),meslist(mprt)

      nb=0
      nm=0

      do 10 i=low,high
         if(abs(itypnew(i)).le.maxbar) then
          nb=nb+1
          barlist(nb)=i
       else
          nm=nm+1
          meslist(nm)=i
       endif
 10   continue
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine leadhad(n1l,n2l,nbl)
c
cinput n1l : First index of this string in the newpart-arrays
cinput n2l : Last index of this string in the newpart-arrays
cinput nbl : Baryon number of string system, nbl=3 is e+e-
c
c output : via common-block {\tt newpart.f}
c
c This subroutine determines the position of the leading hadrons
c of a string fragmentation in the newpart-arrays and assigns the
c reduction factor for the total cross sections (for all non-leading
c hadrons this factor is zero, i.e. they have no cross-section during
c their formation time)
c
C####C##1#########2#########3#########4#########5#########6#########7##

      implicit none
      include 'newpart.f'
      include 'comres.f'
      include 'coms.f'
      integer n1l,n2l,nbl,ll

c default: no leading hadron at all
       do 1 ll=n1l,n2l
         leadfac(ll)=0.d0
 1    continue

c no leading quarks in e+e-
      if(nbl.eq.3) return

c full x-section, if there is only one particle
      if(n1l.eq.n2l)then
        leadfac(n1l)=1.d0
        return
      endif

c the last hadron contains always one leading quark
      if(iabs(itypnew(n2l)).le.maxbar)then
        leadfac(n2l)=0.33d0
      else
        leadfac(n2l)=0.5d0
      endif
c baryonic string: look for the first (=leading) baryon
      if(nbl.gt.0)then
        do 10 ll=n1l,n2l
         if(iabs(itypnew(ll)).le.maxbar)then
          leadfac(ll)=0.66d0
          goto 20
         endif        
 10    continue
 20    continue
      else
c mesonic string: the first hadron contains one leading quark
        if(iabs(itypnew(n1l)).le.maxbar)then
          leadfac(n1l)=0.33d0
        else
          leadfac(n1l)=0.5d0
        endif
      endif

      return
      end

