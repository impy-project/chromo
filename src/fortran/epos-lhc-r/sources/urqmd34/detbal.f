c $Id: detbal.f,v 1.12 1999/01/18 09:57:00 ernst Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine detbal(sqrts,ityp1,ityp2,iso31,iso32,
     &                  em1,em2,itnew1,itnew2,dbfact)
c
c     Revision : 1.0 
c
cinput sqrts   : sqrt(s)
cinput ityp1   : ityp of incoming particle 1
cinput ityp2   : ityp of incoming particle 2
cinput iso31   : 2*I3 of incoming particle 1
cinput iso32   : 2*I3 of incoming particle 2
cinput em1     : mass of incoming particle 1
cinput em2     : mass of incoming particle 2
cinput itnew1  : ityp of outgoing particle 1
cinput itnew2  : ityp of outgoing particle 2
c
coutput dbfact     : correction factor for cross section
c
c     This subroutine calculates a correction factor for the 
c     partial crosssection based on the principle of detailed balance.
C
c     For {\tt CTOption(3)=0} a modified detailed balance (default) is used
c     which takes finite resonance widths into account. For
c     {\tt CTOption(3)=1} the old standard detailed balance relation is used.
c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      implicit none

      include 'coms.f'
      include 'options.f'
      include 'comres.f'
      include 'newpart.f'
c

      real*8 sqrts, dbfact
      integer iso31, iso32
      real*8  em1, em2, clebweight
      integer ityp1, ityp2, itnew1, itnew2
c     local vars for integration of Breit-Wigner:
      real*8 mepsilon
      real*8 oq, q, minwid, factor
      integer inres, outres,idum1,idum2,idum3,idum4
c     called functions
      real*8 pcms, dbweight, massit
      real*8 pmean, widit, dgcgkfct
      integer isoit
      real*8 detbalin
      external detbalin 

c
c     0.1 MeV shift for integrator-maxvalue
      parameter(mepsilon=0.0001)
c     minimal width for "unstable" particle
      parameter( minwid=1.d-4 )

c
      idum1=0
      idum2=0
      idum3=0
      idum4=0
c
c     fix itypes, iso3 and  phase-space for outgoing particles
c
c     a) set up call to isocgk and getmass: determine outgoing isospins
c        and masses


c
c clebweight: actually areduction of given isospin_summed cross_section 
c             to actual incoming channel - here it is used to probe
c             wether the process in question is isospin allowed or not
c
      clebweight=dbweight(isoit(ityp1),iso31,isoit(ityp2),iso32,
     &     isoit(itnew1),isoit(itnew2))
      if(clebweight.lt.0.00001) then
         dbfact=0.d0
         return
      endif


c
c     b) determine momenta
c
      pnnout=pcms(sqrts,massit(itnew1),massit(itnew2))
c
c     d) now calculate correction factor
c

c
c     call to dgcgkfct which calculates degeneracy factors and clebsches
c     
      factor=dgcgkfct(ityp1,ityp2,iso31,iso32,itnew1,itnew2)


      if(CTOption(3).eq.0) then

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c modified detailed balance 
c
c reference: Danielewicz and Bertsch: Nuclear Physics A533(1991) 712.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c
c     count resonances in the incoming channel:
c
         inres=0
         if(widit(ityp1).gt.minwid) inres=inres+1
         if(widit(ityp2).gt.minwid) inres=inres+1

c
c     count resonances in the outgoing channel:
         outres=0
         if(widit(itnew1).gt.minwid) outres=outres+1
         if(widit(itnew2).gt.minwid) outres=outres+1

         if(inres.eq.0) then
c     in this case detbal without resonances, original prescription
            dbfact=factor*pnnout**2/(pnn**2)
c
cccccccccccccccccccccccccccccc
         elseif(inres.eq.1) then
c modified det-bal for one resonance
c

c     now generate the correction factor
         dbfact=factor*pnnout**2/
     &          pmean(sqrts,ityp1,iso31,ityp2,iso32,
     &                idum1,idum2,idum3,idum4,2)

cccccccccccccccccccccccccc
      else
c     modified det-bal for two resonances
c     reference: S.A. Bass, private calculation
c
ccccccccccccccccccccccccccccccccc
         oq=0D0
         if(outres.gt.0) then

c here we have B* B* to B* N
            oq=pmean(sqrts,itnew1,-99,itnew2,-99,
     &               idum1,idum2,idum3,idum4,2)

         endif
ccccccccccccccccccccccccccccccccc

         q=pmean(sqrts,ityp1,iso31,ityp2,iso32,
     &           idum1,idum2,idum3,idum4,2)
c
c     now generate the correction factor
            if(outres.eq.0) then
               dbfact=factor*pnnout**2/(max(1.d-12,q))
            else
               dbfact=factor*oq/(max(1.d-12,q))
            endif
         endif

         return
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c original detailed balance
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif(CTOption(3).eq.1) then
         dbfact=factor*pnnout**2/(pnn**2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c error processing
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
         write(6,*)'undefined detailed balance mode in DETBAL'
         dbfact=1.
      endif
c
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      real*8 function ppiso(pid,ityp1,iso31,ityp2,iso32,itnew1,itnew2)
c
c     Revision : 1.0 
c
cinput pid     : ID of process
cinput ityp1   : ityp of incoming particle 1
cinput ityp2   : ityp of incoming particle 2
cinput iso31   : 2*I3 of incoming particle 1
cinput iso32   : 2*I3 of incoming particle 2
cinput itnew1  : ityp of outgoing particle 1
cinput itnew2  : ityp of outgoing particle 2
c
coutput nniso     : isospin-factor for the reaction $p p \to B B$
c
c     This subroutine calculates the isospin-factor for resonance
c     excitation in inelastic proton proton collisions.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      implicit none

      real*8 dbweight,factor,dgcgkfct
      integer isoit,pid,itag,iz1,iz2,i1,i2,im,jm,ityp1,ityp2
      integer iso31,iso32,itnew1,itnew2,itmp1,itmp2

      include 'comres.f'
      include 'newpart.f'


      i1=ityp1
      i2=ityp2
      iz1=iso31
      iz2=iso32
      im=itnew1
      jm=itnew2

      if(pid.gt.0) then
         ppiso=dbweight(i1,iz1,i2,iz2,isoit(im),isoit(jm))/
     /        dbweight(1,1,1,1,isoit(im),isoit(jm))
      else

         factor=dgcgkfct(i1,i2,iz1,iz2,nucleon,nucleon)
         if(factor.le.1.d-8) then
            ppiso=0.d0
            return
         endif

         nexit=2
         itot(1)=isoit(nucleon)
         itot(2)=isoit(nucleon)
         call isocgk4(isoit(i1),iz1,isoit(i2),iz2,itot,i3new,itag)
         itmp1=i3new(1)
         itmp2=i3new(2)
         ppiso=dbweight(nucleon,itmp1,nucleon,itmp2,
     &           isoit(im),isoit(jm))/
     /           dbweight(1,1,1,1,isoit(im),isoit(jm))
      endif

      return

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      real*8 function dgcgkfct(ityp1,ityp2,iso31,iso32,itnew1,itnew2)
c
c     Revision : 1.0 
c
cinput ityp1   : ityp of incoming particle 1
cinput ityp2   : ityp of incoming particle 2
cinput iso31   : 2*I3 of incoming particle 1
cinput iso32   : 2*I3 of incoming particle 2
cinput itnew1  : ityp of outgoing particle 1
cinput itnew2  : ityp of outgoing particle 2
c
coutput dgcgkfct     : product of degeneracy and cgk factor for detailed bal.
c
c     This subroutine calculates the product of the spin and isospin
c     degeneracy factors and the a isospin correction factor
c     (isospin dependence of cross section) for the detailed balance
c     cross section.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      implicit none

      include 'comres.f'
c

      integer iso31,iso32
      real*8  clebweight
      integer ityp1, ityp2,stot(4),itnew1,itnew2
c     components of degeneracy factor
      integer gin1,gin2,gout1,gout2
      real*8 dgfact
c     called functions
      real*8 dbweight
      integer jit,isoit
c
c     a) set up call to isocgk and getmass: determine outgoing isospins
c        and masses

c
c clebweight: reduction of given isospin_summed cross_section to actual
c             incoming channel
c
c 1) reduction:
      clebweight=dbweight(isoit(ityp1),iso31,isoit(ityp2),iso32,
     &     isoit(itnew1),isoit(itnew2))
      if(clebweight.lt.0.00001) then
         dgcgkfct=0.d0
         return
      endif


c     c) calculate degeneracy factors 
c        reference: S. Bass, GSI-Report 93-13 p. 25 and references therein
c
c     get spins: in-channel stot(1 and 2), out-channel stot(3 and 4)
      stot(1)=jit(ityp1)
      stot(2)=jit(ityp2)
      stot(3)=jit(itnew1)
      stot(4)=jit(itnew2)
c

      gout1=(stot(3)+1)
      gout2=(stot(4)+1)
      gin1=(stot(1)+1)
      gin2=(stot(2)+1)
c
c
c     the degeneracy factor is 
      dgfact=dble(gout1*gout2)/dble(gin1*gin2)
c
c     d) now calculate correction factor
c
c
      dgcgkfct=dgfact*clebweight
c 
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function pmean(sqrts,itp1,iz1,itp2,iz2,
     &                      itp3,iz3,itp4,iz4,ipwr)
c
c     Revision : 1.0
c
cinput sqrts  : $\sqrt{s}$
cinput itp1   : ityp of particle 1
cinput iz1    : $2 \cdot I_3$ of particle 1
cinput itp2   : ityp of particle 2
cinput iz2    : $2 \cdot I_3$ of particle 2
cinput itp3   : ityp of particle 3
cinput iz3    : $2 \cdot I_3$ of particle 3
cinput itp4   : ityp of particle 4
cinput iz4    : $2 \cdot I_3$ of particle 4
cinput ipwr   : power of $p_{mean}$ to integrate
c
c     This function returns the value of the following integral:
c     \begin{displaymath}
c     \int\limits_{m_1= {\tt mmin}}^{\tt mmax} 
c      p_{CMS}^{\tt ipwr}(\sqrt{s},m_1,m_2) A_1(m_1) A_2(m_2) \; dm_1 dm_2 
c     \end{displaymath}
c     with $A_r(m)$ being the spectral function of the resonance:
c     \begin{displaymath}
c      A(m) = \displaystyle\frac{\Gamma(m)/2}{(m-m_0)^2+\Gamma(m)^2/4}
c     \end{displaymath}
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

c      include "comres.f"

      real*8 sqrts,minwid,q1,q2
      real*8 mmin1,mfix,mmin2,maxs1,maxs2,mepsilon,diverg1
      real*8 smass
      integer itp1,ipwr,iz1,itp2,iz2,itp3,iz3,itp4,iz4
      integer inres,izz1,ind1,ind2

cfunctions
      real*8 detbalin,widit,massit,detbalin2,mminit,pcms
c      integer
      external detbalin,detbalin2

c     minimal width for "unstable" particle
      parameter( minwid=1.d-3 )
c     0.1 MeV shift for integrator-maxvalue
      parameter(mepsilon=0.0001)


      if(sqrts.le.mminit(itp1)+mminit(itp2)) then
         pmean=0.d0
         return
      endif

c     count broad particles and store number in inres
c     NOTE: only particles 1 and 2 may be broad!!!!
c
      inres=0
      if(widit(itp1).gt.minwid) inres=inres+1
      if(widit(itp2).gt.minwid) inres=inres+1

      if(inres.eq.0) then
c     in this case the Breit-Wigner distributions are Delta functions,
c     no integrations necessary
         smass=massit(itp2)
         if(itp3.ne.0) smass=smass+massit(itp3)
         if(itp4.ne.0) smass=smass+massit(itp4)

         pmean=pcms(sqrts,massit(itp1),smass)**ipwr

         return
c
cccccccccccccccccccccccccccccc
      elseif(inres.eq.1) then
c modified det-bal for one resonance
c
c first determine which particle is the resonance and store ityp in ind1
         if(widit(itp1).gt.minwid) then
            ind2=itp2
            ind1=itp1
            izz1=iz1
         else
            ind2=itp1
            ind1=itp2
            izz1=iz2
         endif

c     now set integration boundaries
         mmin1=mminit(ind1)
         mfix=mminit(ind2)

         if(itp3.ne.0) mfix=mfix+massit(itp3)
         if(itp4.ne.0) mfix=mfix+massit(itp4)

         maxs1=sqrts-mfix-mepsilon
c     the integration might be divided by the pole of the Breit-Wigner
c     then two integrations with diverg1 as upper or lower boundary
c     respectively are necessary
         diverg1=massit(ind1)
c
c     now perform integration in function detbalin
c     integrate f(m1)=pcms(sqrts,m1,m2)**ipwr*fbwnorm(m1,ityp1)
         q1=0.d0
         q2=0.d0
         if(mmin1.le.diverg1) then
            if(maxs1.gt.diverg1) then
               call qsimp(detbalin,mmin1,diverg1,
     &              ind1,izz1,mfix,sqrts,ipwr,q1,-1)
               call qsimp(detbalin,diverg1,maxs1,
     &              ind1,izz1,mfix,sqrts,ipwr,q2,1)
            else
               call qsimp(detbalin,mmin1,maxs1,
     &              ind1,izz1,mfix,sqrts,ipwr,q1,-1)
            endif
         else
               call qsimp(detbalin,mmin1,maxs1,
     &              ind1,izz1,mfix,sqrts,ipwr,q2,1)
         endif

         pmean=(q1+q2)

         return
c
cccccccccccccccccccccccccc
      else
c 2 resonances to integrate over
c
c
         if(itp3.ne.0) then
            write(6,*) 'ERROR in pmean: only one broad particle allowed'
            write(6,*) '                in case of 3 or 4 body decays!!'
            stop 137
         endif


c     outer integration:
c     set integration boundaries
         mmin1=mminit(itp1)
         mmin2=mminit(itp2)
         maxs2=sqrts-mmin1
         diverg1=massit(itp2)
         q1=0.d0
         q2=0.d0
         if(mmin2.le.diverg1) then
            if(maxs2.gt.diverg1) then             
               call qsimp2(detbalin2,mmin2,diverg1,
     &              itp2,iz2,mmin1,itp1,iz1,ipwr,sqrts,q1,-1)
               call qsimp2(detbalin2,diverg1,maxs2,
     &              itp2,iz2,mmin1,itp1,iz1,ipwr,sqrts,q1,1)
            else
               call qsimp2(detbalin2,mmin2,maxs2,
     &              itp2,iz2,mmin1,itp1,iz1,ipwr,sqrts,q1,-1)
            endif
         else
            call qsimp2(detbalin2,mmin2,maxs2,
     &           itp2,iz2,mmin1,itp1,iz1,ipwr,sqrts,q1,1)
         endif

         pmean=(q1+q2)

         return
      endif


      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function detbalin(m1,ityp1,iz1,m2,sqrts,ipwr)
c
c     Revision : 1.0
c
cinput m1     : mass of resonance (integration variable)
cinput ityp1  : ityp of Delta/N* resonance for fbwnorm()
cinput iz1    : $2\cdot I_3$ of resonance
cinput m2     : second mass for call to pcms
cinput sqrts  : sqrt(s)
cinput ipwr   : power for $p_mean$
c
c     This function is an integrand for the modified detailed balance: 
c     \begin{displaymath}
c     detbalin(m1)=\, p_{CMS}^{\tt ipwr}(\sqrt{s},m_1,m_2) A_1(m_1)
c     \end{displaymath}
c     with $A_r(M)$ being the spectral function of the resonance:
c     \begin{displaymath}
c      A(m) = \displaystyle\frac{\Gamma(m)/2}{(m-m_0)^2+\Gamma(m)^2/4}
c     \end{displaymath}
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     arguments
      real*8 m1,m2,sqrts
      integer ityp1,iz1,ipwr
c     called functions
      real*8 fbwnorm,pcms

      detbalin=pcms(sqrts,m1,m2)**ipwr*fbwnorm(m1,ityp1,iz1)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function detbalin2(m2,ityp2,iz2,min1,
     &                          ityp1,iz1,ipwr,sqrts)
c
c     Revision : 1.0
c
c     This function represents the integrand
c     \begin{displaymath}
c     detbalin2\,=\,A_2(m2)\; 
c     \int \limits_{m_N + m_{\pi}}^{\sqrt{s} - m_2}
c     p_{rel}(\sqrt{s},m_1,m_2) A_1(m_1)  \, d m_1
c     \end{displaymath}
c     for the modified detailed balance with two resonances in the
c     incoming channel.
c     $A_r(M)$ is the spectral function of the resonance (see {\tt detbalin}).
c
c
cinput m2     : mass of resonance2 (outer integration in detbal)
cinput ityp1  : ityp of resonance1
cinput ityp2  : ityp of resonance2
cinput min1   : lower boundary for integration via {\tt qsimp}
ccinput max1   : upper boundary for integration via {\tt qsimp}
cinput iz1    : $2\cdot I_3$ of resonance 1
cinput iz2    : $2\cdot I_3$ of resonance 2
cinput ipwr   : power for $p_mean$
cinput sqrts  : $\sqrt{s}$
c
c     output: 
c             detbalin2  : value of function
c
c     function and subroutine calls:
c                      detbalin (referenced as external)
c                      fbwnorm
c                      qsimp
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ityp1,ityp2,iz1,iz2,ipwr
      real*8 m2,min1,max1,sqrts,q1,q2,diverg1,fbwnorm,mepsilon
      real*8 massit
      real*8 detbalin
      external detbalin
c     0.1 MeV shift for integrator-maxvalue
      parameter(mepsilon=0.0001)

      max1=sqrts-m2-mepsilon

      diverg1=massit(ityp1)
      q1=0.d0
      q2=0.d0
      if(min1.le.diverg1) then 
         if(max1.gt.diverg1) then
            call qsimp(detbalin,min1,diverg1,
     &           ityp1,iz1,m2,sqrts,ipwr,q1,-1)
            call qsimp(detbalin,diverg1,max1,
     &           ityp1,iz1,m2,sqrts,ipwr,q2,1)
         else
            call qsimp(detbalin,min1,max1,
     &           ityp1,iz1,m2,sqrts,ipwr,q1,-1)
         endif
      else
            call qsimp(detbalin,min1,max1,
     &           ityp1,iz1,m2,sqrts,ipwr,q2,1)
      endif
            
      detbalin2=fbwnorm(m2,ityp2,iz2)*(q1+q2)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function fbwnorm(m,ires,iz1)
c
c     Revision : 1.0 
c
cinput m   : mass of resonance
cinput ires: ityp of resonance
cinput iz1 : $2\cdot I_3$ of resonance
c
c     {\tt fbwnorm} returns a Breit-Wigner distribution (non-relativistic)
c     which is normalized to 1 in the limit of mass-independent
c     decay widths. However this function uses mass-dependent decay
c     widths when available. The function only uses widths down to 
c     a lower boundary of 1 MeV, smaller widths are automatically set
c     to 1 MeV. For {\tt iz=-99} fixed widths are used instead of
c     a call to {\tt fwidth}. You should use {\tt fbrwig} for standard purpose
c       since in case of mass dependent widths fbwnorm() is not very well 
c       defined for widths smaller than 1 MeV.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit none

      real*8 m,gam2,mres,fwidth,massit,widit,gam,minwid
      integer ires,ires1,iz1
      include 'comres.f'
      include 'coms.f'
      include 'comwid.f'

c     minimal width for "unstable" particle
      parameter( minwid=1.d-3 )

      ires1 = ires
      mres = massit(ires1)
      if(iz1.eq.-99.or.wtabflg.eq.0)then
        gam = widit(ires1)
      else
        gam = fwidth(ires1,iz1,m)
      end if
c     cutoff for small widths
      gam=max(gam,minwid)
      gam2=gam**2
      fbwnorm = 0.5*gam/(pi*((m-mres)**2+gam2/4.0))!*norm
      return 
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  here start the numerical receipies routines for numerical integration
c  no more physics beyond this point in the file!!!
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c(c) numerical receipies, adapted for f(x,idum,dum,dum)
      SUBROUTINE qsimp(func,a,b,idum1,idum2,dum2,dum3,idum3,s,flag)
c
c     Simpson integration via Numerical Receipies.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'options.f'

      INTEGER JMAX,j,idum1,idum2,idum3,flag
      REAL*8 a,b,func,s,EPS
      REAL*8 os,ost,st,dum2,dum3
      EXTERNAL func
      PARAMETER (JMAX=20)

      if(b-a.le.1.d-4) then
         s=0.d0
         return
      endif

      EPS = 6.d-2
      if (CTOption(35).eq.1) EPS=6.d-3

      ost=-1.d30
      os= -1.d30
      do 11 j=1,JMAX
         if(flag.eq.-1) then
            call midsqu1(func,a,b,idum1,idum2,dum2,dum3,idum3,st,j)
         elseif(flag.eq.1) then
            call midsql1(func,a,b,idum1,idum2,dum2,dum3,idum3,st,j)
         endif
        s=(9.*st-ost)/8.
        if (abs(s-os).le.EPS*abs(os)) return
        os=s
        ost=st
11    continue
      write(6,*)  'too many steps in qsimp, increase JMAX!'

      return      
      END

C  (C) Copr. 1986-92 Numerical Recipes Software.



      SUBROUTINE midsqu1(funk,aa,bb,idum1,idum2,dum2,dum3,idum3,s,n)
c     modified midpoint rule; allows singuarity at upper limit
      implicit none
      integer idum1,idum2,idum3
      real*8 dum2,dum3
      INTEGER n
      REAL*8 aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x,func,a,b,xx

      func(x)=2.*x*funk(bb-x**2,idum1,idum2,dum2,dum3,idum3)

      b=sqrt(bb-aa)
      a=0.d0
      if (n.eq.1) then
      xx=0.5d0*(a+b)

      s=(b-a)*func(0.5d0*(a+b))

      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END

      SUBROUTINE midsql1(funk,aa,bb,idum1,idum2,dum2,dum3,idum3,s,n)
c     modified midpoint rule; allows singularity at lower limit
      implicit none
      integer idum1,idum2,idum3
      real*8 dum2,dum3
      INTEGER n
      REAL*8 aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x,func,a,b
      func(x)=2.*x*funk(aa+x**2,idum1,idum2,dum2,dum3,idum3)
      b=sqrt(bb-aa)
      a=0.
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c(c) numerical receipies, adapted for f(x,idum,dum,dum)
      SUBROUTINE qsimp2(func,a,b,idum1,idum2,dum1,idum3,idum4,
     &                  idum5,dum2,s,flag)
c
c     Simpson integration via Numerical Receipies.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      include 'options.f'

      INTEGER JMAX,j,idum1,idum2,idum3,idum4,idum5,flag
      REAL*8 a,b,s,EPS
      REAL*8 os,ost,st,dum1,dum2
      REAL*8 func
      PARAMETER (JMAX=10)
      external func

      if(b-a.le.1.d-4) then
         s=0.d0
         return
      endif

      EPS = 6.d-2
      if (CTOption(35).eq.1) EPS=6.d-3

      ost=-1.d30
      os= -1.d30
      do 11 j=1,JMAX
         if(flag.eq.-1) then
            call midsqu2(func,a,b,idum1,idum2,dum1,idum3,idum4,
     &                  idum5,dum2,st,j)
         elseif(flag.eq.1) then
            call midsql2(func,a,b,idum1,idum2,dum1,idum3,idum4,
     &                  idum5,dum2,st,j)
         endif
        s=(9.*st-ost)/8.

        if (abs(s-os).le.EPS*abs(os)) return
        os=s
        ost=st
11    continue

      write(6,*)  'too many steps in qsimp2, increase JMAX!'

      return      
      END


      SUBROUTINE midsqu2(funk,aa,bb,idum1,idum2,dum1,idum3,idum4,
     &                  idum5,dum2,s,n)
c     modified midpoint rule; allows singuarity at upper limit
      implicit none
      integer idum1,idum2,idum3,idum4,idum5
      real*8 dum1,dum2
      INTEGER n
      REAL*8 aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x,func,a,b
      func(x)=2.*x*funk(bb-x**2,
     &                  idum1,idum2,dum1,idum3,idum4,idum5,dum2)
      b=sqrt(bb-aa)
      a=0.
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END

      SUBROUTINE midsql2(funk,aa,bb,idum1,idum2,dum1,idum3,idum4,
     &                  idum5,dum2,s,n)
c     modified midpoint rule; allows singularity at lower limit
      implicit none
      integer idum1,idum2,idum3,idum4,idum5
      real*8 dum1,dum2
      INTEGER n
      REAL*8 aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x,func,a,b
      func(x)=2.*x*funk(aa+x**2,
     &                  idum1,idum2,dum1,idum3,idum4,idum5,dum2)
      b=sqrt(bb-aa)
      a=0.
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END







