c $Id: gnuranf.f,v 1.7 1999/11/24 19:47:50 ssoff Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ranf(idummy)
c
c     wrapper for the local random number generator
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'coms.f'

      integer idummy, iseed, oldseed
      real*8 ranf, ran2
      common /seed/iseed,oldseed

      ranf = ran2(iseed)
      if (check) then
         if (ranf.ge.1.0D0.or.ranf.le.0.0D0) then
            call error ('ranf','ranf gives endpoint values',ranf,3)
         endif
      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sseed(ranseed)
c
c     reset the random number generato
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      real*8 dummy, ran2
      integer iseed, ranseed, oldseed, time, timeseed
      common /seed/iseed,oldseed

      if (ranseed.le.0) then
         timeseed = abs(time())
         if (timeseed.eq.oldseed) return
         ranseed = timeseed
      endif
      oldseed = ranseed
      iseed = -ranseed
      dummy = ran2(iseed)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ran2(idum)
c
c Long period (>2E18) random number generator of L'Ecuyer with
c Bays-Durham shuffle and added safeguards. Returns a uniform random
c deviate between 0.0 and 1.0 (exclusive of the endpoint values).
c Call with idum a negative integer to initialize; thereafter, do
c not alter idum between successive deviates in a sequence. RNMX
c should approximate the largest floating value that is less than 1.
c
C  (C) Copr. 1986-92 Numerical Recipes Software.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real*8 ran2,AM,EPS,RNMX

      parameter (IM1=2147483563,IM2=2147483399,AM=1.0D0/IM1,IMM1=IM1-1,
     $     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     $     NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.0D0-16,RNMX=1.0D0-EPS)
      integer idum2, j, k, iv(NTAB), iy
      save iv, iy, idum2
      data idum2/123456789/, iv/NTAB*0/, iy/0/
      
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
 11      continue
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if (iy.lt.1) iy=iy+IMM1
      ran2=min(AM*iy,RNMX)

      return
      end

C  (C) Copr. 1986-92 Numerical Recipes Software 5&40.
