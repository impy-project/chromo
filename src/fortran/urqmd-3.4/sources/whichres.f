c $Id: whichres.f,v 1.7 2007/01/30 14:50:32 bleicher Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function whichres(m,class)
c
c     Revision : 1.0  
cinput   m          : Mass of the resonance
cinput   class      : class of resonance: see subr. getrange
coutput  whichres   : itype of the resonance
c
c     DETERMINES ID OF A RESONANCE WITH MASS M ACCORDING TO THE GLOBAL
c     TYPE OF THE RESONANCE. THE TYPE OF THE RESONANCE IS SELECTED
c     RANDOMLY BETWEEN ALL TYPES PERMITTED.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real*8 m
      integer class,i,im,ip

      call getrange(class,im,ip)

        call whichi(i,im,ip,m)

        whichres=i

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine whichi(i,im,ip,m)
c
c     Revision : 1.0  
cinput   m         : Mass of the resonance
cinput   im,ip   : lower and upper limit of itypes
coutput  i               : itype of the resonance
c
c     DETERMINES ID OF A RESONANCE WITH MASS M ACCORDING TO  
c     the range defined by {\tt im} and {\tt ip}. 
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'comres.f'
      real*8 m,fmax,f(-itmax:itmax),fbrwig,bwnorm
      integer i,im,ip,jit,isoit
      include 'options.f'


      do 101 i=im,ip
c        f(i) = breitwig(m,i)
         f(i) = fbrwig(i,isoit(i),m,1)/bwnorm(i)*dble(jit(i)+1)

 101  continue

      call getbran(f,-itmax,itmax,fmax,im,ip,i)

      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getrange(class,i1,i2)
c
c     Revision : 1.0 
c
cinput  class  : class  of resonance:
c                   0  = Delta(1232)
c                   1  = N* 
c                   2  = Delta* (EXcludind D(1232))
c                   3  = all nonstrange resonances allowed
c                   4  = nucleon
c                        11  = N and N*
c                        12  = Delta(1232) and Delta*
c                        13  = Lambda, Lambda*
c                        14  = Sigma, Sigma*
c                        15  = Cascade, Cascade*
c                        16  = Omega(s)
c
coutput i1,i2  : range of resonance IDs (from,to) 
c
C     INTERFACE BETWEEN class AND paticle ID
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer class,i1,i2
      include 'comres.f'

      if (class.eq.0) then
c  inspect only Delta(1232)
         i1 = mindel
         i2 = mindel
      else if (class.eq.1) then
c  all N*-Resonances
         i1 = minnuc+1
         i2 = maxnuc
      else if (class.eq.2) then
c  all Delta-Resonances (EXcluding D(1232))
         i1 = mindel+1
         i2 = maxdel
      else if (class.eq.3) then
c  all non strange Resonances 
         i1 = minres
         i2 = maxres
      else if (class.eq.4) then
         i1 = nucleon
         i2 = nucleon
      else if (class.eq.11) then
         i1 = minnuc
         i2 = maxnuc
      else if (class.eq.12) then
         i1 = mindel
         i2 = maxdel
      else if (class.eq.13) then
         i1 = minlam
         i2 = maxlam
      else if (class.eq.14) then
         i1 = minsig
         i2 = maxsig
      else if (class.eq.15) then
         i1 = mincas
         i2 = maxcas
      else if (class.eq.16) then
         i1 = minome
         i2 = maxome
      else
c  something went wrong
         write(6,*) 'getrange: class=',class,' not valid...'
         stop 137
      endif
     
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getinw(i,m,class,mmax)
c
c Revision : 1.0 
c
cinput    i   : Resonance ID
cinput   mmax : upper limit for mass to create
coutput   i   : Resonance ID of particle in same class
coutput   m   : Mass of particle in same class
coutput class : Class of resonance i (see {\tt getrange})        
c
c Does in a way the inverse of subr. getrange: look in which class
c resonance i is and determine mass m and itype i of same class with
c maximal mass of mmax.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      include 'comres.f'
      integer i,j,im,ip,isav,class,ii,iim,iip
      real*8 m,mmax

      real*8 massit,getmass
      integer whichres,isign

      isav=i
      do 108 j=0,2
        call getrange(j,im,ip) 
        if(i.ge.im.and.i.le.ip)then
          m=getmass(mmax,j)
          i=whichres(m,j)
          class=j
          return
        end if
 108  continue
c if failed keep old value of i and standard value for m
claw next line is for debug purpose
        write(6,*)'getinw: itype not in resonance range:itype=',isav
      i=isav
      m=massit(i)
      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry getirg(ii,iim,iip)
c
cinput ii       : particle ID
coutput iim,iip : minimal and maximal itype of particle's class
c
c Determine which range (acc. to {\tt getrange}) particle 
c {\tt ii} belongs to.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do 1008 j=11,16 
        call getrange(j,iim,iip) 
        if(iabs(ii).ge.iim.and.iabs(ii).le.iip)then
           iim=min(isign(iim,ii),isign(iip,ii))
           iip=max(isign(iim,ii),isign(iip,ii))
           return
        endif
 1008   continue
c only one particle (ii) within range
      iim=ii
      iip=ii
      return

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function massdist(m,class)
c
c Revision : 1.0 
c
cinput  m        : mass to look at
cinput  class    : resonance class (see {\tt getrange}
coutput massdist : value of mass distribution at mass m
c       
C MASS DISTRIBUTION AT MASS M FOR GIVEN class OF RESONANCES.
C DISTRIBUTION IS CALCULATED BY SUPERPOSING BREIT-WIGNER-
C DISTRIBUTIONS FOR ALL RESONANCES OF A KIND.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real*8 m
      integer class,i1,i2,i,isoit,jit
      include 'comres.f'
      real*8 fbrwig,bwnorm
      include 'options.f'

      call getrange(class,i1,i2)
c  Delta(1232) is included via i1=i2=mindel
      massdist = 0.0
      do 101 i=i1,i2
c         massdist=massdist+breitwig(m,i)
         massdist=massdist+fbrwig(i,isoit(i),m,1)
     &        /bwnorm(i)*dble(jit(i)+1)
 101  continue
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function pcms(ecm,m1,m2)
c calculates the CM-momentum in a 2-body decay/coll. depending on ecm
      implicit none
      real*8 ecm,m1,m2,s

      if (ecm.le.m1+m2) then
         pcms = 0.0
         return
      endif

      s = ecm*ecm
      pcms = sqrt( (s-(m1+m2)**2)*(s-(m1-m2)**2)/(4.0*s))
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function bcms(ecm,m1,m2)
c calculates the CM-velocity in a 2-body decay/coll. depending on ecm

      implicit none
      real*8 ecm,m1,m2,s

      if (ecm.le.m1+m2) then
         bcms = 0.0
         return
      endif

      s = ecm*ecm
      bcms = sqrt( (1d0-(m1+m2)**2/s) * (1d0-(m1-m2)**2/s) )
      return
      end



