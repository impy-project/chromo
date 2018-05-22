c $Id: dectim.f,v 1.9 2007/01/30 14:50:24 bleicher Exp $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function dectim(ind,iproc)
c
c     Revision : 1.0
C
cinput ind : ID of particle
cinput iproc: process ID for resonance creation
couput dectim: time of decay
c
c     This function computes a random choice for the time at which
c     a resonance will decay and transformes it to the computational
c     frame.
c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      implicit none
      include 'coms.f'   
      include 'options.f'
      include 'colltab.f'

      integer ind,iproc
      real*8 gg,wid,tau,ranf,fwidth,widit,fbwnorm,mr,massit
      real*8 tmp,factor


      if(ityp(ind).eq.134.and.iso3(ind).eq.-1)then
       dectim=10.5d0*p0(ind)/fmass(ind)*hqc
       return
      endif
      if(ityp(ind).eq.-134.and.iso3(ind).eq.1)then
       dectim=10.5d0*p0(ind)/fmass(ind)*hqc
       return
      endif

c... keep unknown PYTHIA particles stable
      if(abs(ityp(ind)).gt.1000) then
         DECTIM=1.d34
         return
      endif

c..D-Mesons
      if(abs(ityp(ind)).ge.133)then
         DECTIM=1.d34
         return
      endif

c
c first determine width of resonace
c
      if(CTOption(1).eq.0.and.CTOption(34).ne.1) then
c     mass dependent width
         wid=fwidth(ityp(ind),iso3(ind),fmass(ind))*CTParam(1)
      else
c     fixed width
         wid=widit(ityp(ind))*CTParam(1)
      end if

C  ...  REST FRAME DECAY TIME
      if(WID .GT. 1.d-10) then
         if(CTOption(34).lt.2) then
c     "normal" life time tau=1/gamma
            TAU=-(dLOG(1.d0-RANF(0))/WID)
         else
c     use Danielewicz delay
            if(iproc.ne.36.and.iproc.ne.37) then
c     delay for scattering wave
               TAU=-(dLOG(1.d0-RANF(0))*
     &              fbwnorm(fmass(ind),ityp(ind),iso3(ind))*pi/2.d0)    
            else
c     delay for forward wave
               if(CTOption(34).eq.2) then
                  factor=1.d0/CTParam(58)
               elseif(CTOption(34).eq.3) then
                  factor=(ctsigtot(actcol)-CTParam(58))/CTParam(58)
               elseif(CTOption(34).eq.4) then
                  tmp=dsqrt(2.d0/(wid*3.14d0*
     &                 fbwnorm(fmass(ind),ityp(ind),iso3(ind))))
                  factor=1.d0/(CTParam(58)*tmp)
               else
                  factor=1.d0
               endif
               mr=massit(ityp(ind))
               TAU=-(dLOG(1.d0-RANF(0))*factor*
     &              2.d0*pi*fbwnorm(fmass(ind),ityp(ind),iso3(ind))*
     &              (fmass(ind)-mr)**2/wid**2)
            endif
         endif
      ELSE          
c     stable particle
         DECTIM=1.d34
         RETURN
      END IF
C  ...  APPLY TIME DILATION
C  ...  GAMMA FOR THE RESONANCE RESTFRAME <-> COMP. FRAME TRAFO  
      gg=p0(ind)/fmass(ind)
      DECTIM=TAU*GG*hqc

      RETURN
      END
