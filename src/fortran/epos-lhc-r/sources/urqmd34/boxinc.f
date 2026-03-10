c $Id: boxinc.f,v 1.3 1999/01/18 09:56:53 ernst Exp $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Unit     : all modules using the box
c     Version  : 1.0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       max count of different species
        integer bptmax
        parameter(bptmax=20)

c      counter
        integer cbox
c      Flags
        integer boxflag
        integer edensflag, mbflag
        integer para,solid,mtest
c      number of different species
        integer mbox
c      particle counters
        integer bptpart(bptmax),bptityp(bptmax),bptiso3(bptmax)
        real*8 bptpmax(bptmax)
        real*8 edens
c      edge length of a cube in fm
        real*8  lbox
c      half edge lengt of a cube 
        real*8  lboxhalbe
c      double edge length of a cube
        real*8  lboxd
c      momenta
        real*8 mbp0, mbpx, mbpy, mbpz
        
        common /boxic/ cbox,boxflag, mbox, bptityp, bptiso3, bptpart
        common /boxic/ edensflag, para, solid, mbflag, mtest
        common /boxrc/ lbox, lboxhalbe, lboxd, bptpmax, edens
        common /boxrc/ mbp0, mbpx, mbpy, mbpz
        
