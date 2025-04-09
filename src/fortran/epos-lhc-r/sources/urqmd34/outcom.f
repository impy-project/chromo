c $Id: outcom.f,v 1.3 2000/01/12 16:02:39 bass Exp $
c temporary storage for in-channel  
      real*8 tsqrts,tstot,tsigpart
      real*8 tr0(3),trx(3),try(3),trz(3),ttform(3),txtotfac(3),
     @     tp0(3),tpx(3),tpy(3),tpz(3),tm(3)
      integer tind(3),tityp(3),tiso3(3),tcoll(3),tstrange(3)
      integer tlcoll(3),tcharge(3),torigin(3),tuid(3)

      common/outcom/tsqrts,tstot,tsigpart,
     &     tr0,trx,try,trz,ttform,txtotfac,
     @     tp0,tpx,tpy,tpz,tm,
     &     tind,tityp,tiso3,tcoll,tstrange,
     &     tlcoll,tcharge,torigin,tuid
