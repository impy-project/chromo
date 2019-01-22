c $Id: urqmd.f,v 1.25 2007/01/30 14:50:32 bleicher Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     subroutine UrQMD
c
c      Authors : The UrQMD collaboration 
c                S.A. Bass, M. Belkacem, M. Bleicher, M. Brandstetter,
c                L. Bravina, C. Ernst, L. Gerland, M. Hofmann, 
c                S. Hofmann, J. Konopka, G. Mao, L. Neise, S. Soff,
c                C. Spieles, H. Weber, L.A. Winckelmann, H. Stoecker
c                and W. Greiner
c
c     Revision: 1.2
c
cc    contact address:
cc
cc                     urqmd@th.physik.uni-frankfurt.de
cc
c               
c This is the main module of {\tt urqmd}. It servers as a connection between
c the initialization, the propagation (including the real part of the 
c optical potential) and the collision term (imaginary part of the optical
c potential).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'coms.f'
      include 'comres.f'
      include 'options.f'
      include 'colltab.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'boxinc.f'
c
      integer i,j,k,steps,ii,ocharge,ncharge, mc, mp, noc, it1,it2
      real*8 sqrts,otime,xdummy,st
      logical isstable
      integer stidx,CTOsave
      real*8 Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      common /energies/ Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      integer cti1sav,cti2sav
chp hydro variables
      real*8 thydro_start,thydro,nucrad
      logical lhydro
  
c
c     numerical/technical initialisation
c
      call uinit(0)
c
c  Main program
c

      mc=0
      mp=0
      noc=0
c
c loop over all events
c
      do 10 event=1,nevents
c     start event here
c    

c     time is the system time at the BEGINNING of every timestep
      time = 0.0
chp hydro flag, hydro should be called only once
      lhydro=.true.

c     initialize random number generator
c     call auto-seed generator only for first event and if no seed was fixed
      if(.not.firstseed.and.(.not.fixedseed)) then
         ranseed=-(1*abs(ranseed))
         call sseed(ranseed)
      else
         firstseed=.false.
      endif

      write(6,*)'event# ',event,ranseed

c
c     initialisation of physics quantities 
c
      call init
cbb if we are reading old events, check the success of the read-in:
      if (CTOption(40).ne.0.and.(.not.success)) then
        exit
      endif

chp hydro switch      
      if (CToption(45).eq.1)then
chp hydro start time (nuclei have passed each other)
chp ebeam is only the kinetic energy 
chp CTParam(65) is useful for the variation of the start time 
chp default value is one
       thydro_start=CTParam(65)*2.d0*nucrad(Ap)*sqrt(2.d0*emnuc/ebeam)
       write(*,300) 'hydro starts after',thydro_start
chp lower limit for hydro start time
       if(thydro_start.lt.CTParam(63)) then
        thydro_start=CTParam(63)
        write(6,300) '... extended to',CTParam(63)
       end if
      end if
 300  format(a18,x,f5.2,' fm/c')

c old time if an old fort.14 is used 
      if(CTOption(40).ne.0)time=acttime

c output preparation

c write headers to file
      call output(13)
      call output(14)
      call output(15)
      call output(16)
      if(event.eq.1) then
        call output(17)
        call osc_header
        call osc99_header
      endif
      call osc99_event(-1)

c for CTOption(4)=1 : output of initialization configuration
      if(CTOption(4).eq.1)call file14out(0)

c     participant/spectator model:
      if(CTOption(28).ne.0) call rmspec(0.5d0*bimp,-(0.5d0*bimp))

c     compute time of output
      otime = outsteps*dtimestep
      
c reset time step counter
      steps = 0

c  loop over all timesteps

      do 20  steps=1,nsteps

c store coordinates in arrays with *_t
c this is needed for MD type propagation
         if (eos.ne.0) then
            do 23 j=1,npart
               r0_t(j) = r0(j)
               rx_t(j) = rx(j)
               ry_t(j) = ry(j)
               rz_t(j) = rz(j)
 23         continue
         end if

c we are at the beginning of the timestep, set current time (acttime) 
         acttime = time
c  option for MD without collision term
         if(CTOption(16).ne.0) goto 103

c  Load collision table with next collisions in current timestep
         call colload
c     check for collisions in time-step, nct = # of collisions in table
         if (nct.gt.0) then

c     entry-point for collision loop in case of full colload after every coll. 
 101        continue
            k = 0
c     normal entry-point for collision loop 
 100        continue
c     get next collision
            call getnext(k)

c     exit collision loop if no collisions are left
            if (k.eq.0) goto 102
chp call hydro if start time is reached
            if(CTOption(45).eq.1)then
           
             if(cttime(k).gt.thydro_start.and.lhydro)then
              st=thydro_start-acttime
              call cascstep(acttime,st)
chp all particle arrays will be modified by hydro
              call hydro(thydro_start,thydro)
              acttime=thydro_start
              lhydro=.false.
              if(thydro.gt.1.d-8.or.CTOption(48).eq.1)then
chp full update of collision table
              call colload
                      
              go to 101
              end if
             end if    
            end if 

c  propagate all particles to next collision time
c  store actual time in acttime, propagation time st=cttime(k)-acttime
            st=cttime(k)-acttime
            call cascstep(acttime,st)
c  new actual time (for upcoming collision)
            acttime = cttime(k)

c  perform collision 

            if(cti2(k).gt.0.and.
     .           abs(sqrts(cti1(k),cti2(k))-ctsqrts(k)).gt.1d-3)then
               write(6,*)' ***(E) wrong collision update (col) ***'
               write(6,*)cti1(k),cti2(k),
     .              ctsqrts(k),sqrts(cti1(k),cti2(k))
            else if(cti2(k).eq.0.and.
     .              abs(fmass(cti1(k))-ctsqrts(k)).gt.1d-3) then
               write(6,*)' *** main(W) wrong collision update (decay)'
               write(6,*)ctag,cti1(k),ityp(cti1(k)),dectime(cti1(k)),
     .              fmass(cti1(k)),ctsqrts(k)
            endif

            ocharge=charge(cti1(k))
            if(cti2(k).gt.0) ocharge=ocharge+charge(cti2(k))

c     store quantities in local variables for charge conservation check
            it1= ityp(cti1(k))
            if(cti2(k).gt.0)it2= ityp(cti2(k))

c increment "dirty" collision counter
            if(cti2(k).gt.0)then !scatter
               mc=mc+1
            endif
c     perform scattering/decay
            cti1sav = cti1(k)               
            cti2sav = cti2(k)    
            call scatter(cti1(k),cti2(k),ctsigtot(k),ctsqrts(k),
     &                   ctcolfluc(k))

c
c  update collision table 
c
c     normal update mode
            if(CTOption(17).eq.0) then
               if(nexit.eq.0) then
c     new collision partners for pauli-blocked states (nexit=0)
                  if (cti1(k).ne.cti1sav.or.cti2(k).ne.cti2sav) then
                   cti1(k) = cti1sav 
                   cti2(k) = cti2sav 
                  endif
                  call collupd(cti1(k),1)
                  if(cti2(k).gt.0) call collupd(cti2(k),1)
               else
                  ncharge=0
c     new collision partners for scattered/produced particles (nexit><0)
                  do 30 i=1,nexit
c     ncharge is used for charge conservation check
                     ncharge=ncharge+charge(inew(i))
                     call collupd(inew(i),1)
 30               continue

c     charge conservation check
                  if(ocharge.ne.ncharge) then
                     write(6,*)'ch-conservation error coll/dec ',ctag
                     write(6,*)'   it1:',it1,'   it2:',it2
                     write(6,*)'   ch:',ocharge,ncharge
                     write(6,*)'cti1(k),cti2(k),ctsigtot(k),ctsqrts(k)'
                     write(6,*)cti1(k),cti2(k),ctsigtot(k),ctsqrts(k)
                  endif
               endif

c     update collisions for partners of annihilated particles
               do 55 ii=1,nsav
                  call collupd(ctsav(ii),1)
 55            continue
               nsav=0

            else ! (CTOption(17).ne.0)
c     full collision load
               call colload
            endif

            if (CTOption(17).eq.0) goto 100
            goto 101

c     this is the point to jump to after all collisions in the timestep
c     have been taken care of
 102        continue

         endif ! (nct.gt.0)

c  After all collisions in the timestep are done, propagate to end of 
c  the timestep.

c     point to jump to in case of MD without collision term
 103     continue

c     increment timestep
         time = time+dtimestep

c  After all collisions in the timestep are done, propagate to end of 
c  the timestep.
         call cascstep(acttime,time-acttime)

c     in case of potential interaction, do MD propagation step
         if (eos.ne.0) then

c set initial conditions for MD propagation-step
            do 24 j=1,npart
               r0(j) = r0_t(j)
               rx(j) = rx_t(j)
               ry(j) = ry_t(j)
               rz(j) = rz_t(j)
 24         continue

c now molecular dynamics trajectories
            call proprk(time,dtimestep)

         end if ! (eos.ne.0)

c     perform output if desired
         if(mod(steps,outsteps).eq.0.and.steps.lt.nsteps)then 
            if(CTOption(28).eq.2)call spectrans(otime)
            call file14out(steps)
           if(CTOption(55).eq.1)then
            call osc_vis(steps)
           endif  
         endif ! output handling

 20   continue ! time step loop


ce
        acttime=time
c optional decay of all unstable particles before final output
c DANGER: pauli-blocked decays are not performed !!!
         if(CTOption(18).eq.0) then
c no do-loop is used because npart changes in loop-structure
            i=0
            nct=0
            actcol=0
c disable Pauli-Blocker for final decays
            CTOsave=CTOption(10)
            CTOption(10)=1
c decay loop structure starts here
 40         continue
            i=i+1

c is particle unstable
            if(dectime(i).lt.1.d30) then
 41            continue
               isstable = .false.
               do 44 stidx=1,nstable
                  if (ityp(i).eq.stabvec(stidx)) then
c                     write (6,*) 'no decay of particle ',ityp(i)
                     isstable = .true.
                  endif
 44            enddo
               if (.not.isstable) then
c     perform decay
                  call scatter(i,0,0.d0,fmass(i),xdummy)
c     backtracing if decay-product is unstable itself
                  if(dectime(i).lt.1.d30) goto 41
               endif
            endif
c     check next particle
            if(i.lt.npart) goto 40
         endif ! final decay
         CTOption(10)=CTOsave
c final output

           if(CTOption(28).eq.2)call spectrans(otime)

         call file13out(nsteps)
         if(CTOption(50).eq.0)then
          call file14out(nsteps)
         end if
         call file16out
         if(CTOption(50).eq.0.and.CTOption(55).eq.0)then
          call osc_event
         end if
         if(CTOption(50).eq.0.and.CTOption(55).eq.1)then
          call osc_vis(nsteps)
         end if
         call osc99_event(1)
         call osc99_eoe
      
         mp=mp+npart
         if(ctag.eq.0)then
           write(*,*)'(W) No collision in event ',event
           noc=noc+1
         endif

c     end of event loop
 10   continue

      write(6,301)'no. of collisions = ',mc/dble(nevents),' (per event)'
      write(6,301)'final particles   = ',mp/dble(nevents),' (per event)'
      write(6,302)'empty events      : ', noc,noc*1d2/dble(nevents)
 301  format(a19,f8.1,a12)
 302  format(a19,i8,' = ',f5.1,'%')
      end

 

