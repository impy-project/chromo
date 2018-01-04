c $Id: comwid.f,v 1.15 2007/01/30 14:50:24 bleicher Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Revision : 1.0
c
c     common block for the tabulated branching ratios
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c number of points for the interpolation
      integer widnsp
c the branching ratios are interpolated with high precision in the range
c from mintab to maxtab1 and from maxtab1 to maxtab2 with a smaller precision

c version number of the table
      integer tabver
      real*8 mintab, maxtab1, maxtab2

c set the default values 
      parameter (widnsp=120)
      parameter (mintab=0.1d0)
      parameter (maxtab1=5.0d0)
      parameter (maxtab2=50.0d0)
c increase this parameter, if you make changes, which require a new table 
      parameter (tabver=9)

c tabulated x-values (i.e. sqrt(s) of the collision)
      real*8 tabx (1:widnsp)
c tabulated y-values (i.e. the branching ratios) and the second
c derivatives of the function. 
c full baryon ratio 
      real*8 fbtaby (1:widnsp,minbar:maxbar,1:2)
c partial baryon ratios
      real*8 pbtaby (1:widnsp,1:2,minbar:maxbar,0:maxbra)
c full meson ratio
      real*8 fmtaby (1:widnsp,minmes:maxmes,1:2)
c partial meson ratios
      real*8 pmtaby (1:widnsp,1:2,minmes:maxmes,0:maxbrm)

c Breit-Wigner norms
c norm of Breit-Wigner with mass dependent widths baryons/mesons
        real*8 bwbarnorm(minbar:maxbar),bwmesnorm(minmes:maxmes)

c tabulated fppfit()
c tabulated x-values (i.e. sqrt(s) of the collision)
      real*8 tabxnd (1:widnsp)
c 2-resonance channels
c                      x     deriv ND N*..D*
      real*8 frrtaby(1:widnsp,1:2,1:2,2:maxdel)

c this flag indicates the progress of tabulating the function
      integer wtabflg
c name of file containing the tables
      character*77 tabname

      common /decaywidth/ tabx,fbtaby,pbtaby,fmtaby,pmtaby,wtabflg
        common /brwignorm/ bwbarnorm,bwmesnorm
        common /xsections/ tabxnd,frrtaby
      common /tabnames/ tabname

