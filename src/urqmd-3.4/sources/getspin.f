c $Id: getspin.f,v 1.5 2007/01/30 14:50:24 bleicher Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      integer function getspin(iityp,itag)
c
c     Revision : 1.0
c
cinput ityp   : ID of particle
cinput itag   : flag for return value
c
c output: $2*J_{tot}$ of particle
c
c     This subroutine converts global ityp to maximum spin and optionally
c     chooses a random projection ({\tt itag=-1}).
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      implicit none
      include 'comres.f'
c
      integer ityp,iityp,jtot,itag
      real*8 ranf
c
      ityp=abs(iityp)

      if(ityp.ge.nucleon.and.ityp.le.maxbar) then
         jtot=Jres(ityp)
      elseif(ityp.ge.offmeson.and.ityp.le.maxmeson) then
         jtot=Jmes(ityp)
      elseif(ityp.gt.1000.and.ityp.lt.1e8) then
c... quick and dirty fix for PDGID+1000 codes: 
c... needs to be modified if spin is needed
         getspin=0
         return                
      else
         write(6,*)'undefined total isospin in getspin:'
         write(6,*)'ityp: ',iityp
         stop 137
      endif
c
      if(itag.eq.1) then
         getspin=jtot
      elseif(itag.eq.-1) then
         getspin=jtot-2*int(ranf(0)*(jtot+1))
      else
         write(6,*)'itag-error in getspin.f'
         stop 137
      endif

      return
      end




