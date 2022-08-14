c $Id: inputs.f,v 1.4 2000/01/12 16:02:36 bass Exp $
c     include file for data from the input subroutine
      integer nevents,spityp(2),prspflg
      integer trspflg,spiso3(2),outsteps,bflag,srtflag,efuncflag
      integer nsrt,npb,firstev
      real*8  srtmin,srtmax,pbeam,betann,betatar,betapro
      real*8  pbmin,pbmax

      common /inputs/nevents,spityp,prspflg,trspflg,
     &               spiso3,outsteps,bflag,srtflag,efuncflag,nsrt,
     &               firstev, npb
      common /input2/ srtmin,srtmax,pbeam,betann,betatar,betapro,
     &                pbmin, pbmax

c maximum number of test particle per nucleon
      integer MaxNTest
      parameter (MaxNTest=20)
c
c maximum mass for projectile and target
      integer AAmax
      parameter(AAmax=300)
c storage arrays for projectile and target nuclei
      integer PT_iso3(AAmax*MaxNTest,2),PT_ityp(AAmax*MaxNTest,2)
      integer PT_spin(AAmax*MaxNTest,2),PT_charge(AAmax*MaxNTest,2)
      integer PT_AA(2),PT_uid(AAmax*MaxNTest,2)
      real*8 PT_r0(AAmax*MaxNTest,2),PT_rx(AAmax*MaxNTest,2)
      real*8 PT_ry(AAmax*MaxNTest,2),PT_rz(AAmax*MaxNTest,2)
      real*8 PT_p0(AAmax*MaxNTest,2),PT_px(AAmax*MaxNTest,2)
      real*8 PT_py(AAmax*MaxNTest,2),PT_pz(AAmax*MaxNTest,2)
      real*8 PT_dectime(AAmax*MaxNTest,2),PT_fmass(AAmax*MaxNTest,2)
      real*8 PT_rho(AAmax*MaxNTest,2),PT_pmax(AAmax*MaxNTest,2)

c common blocks (data transfer between cascinit and getinit)
      common/ProTarInts/PT_iso3,PT_ityp,PT_spin,PT_charge,PT_AA,PT_uid
      common/ProTarReals/PT_r0,PT_rx,PT_ry,PT_rz,PT_fmass,PT_dectime,
     &                   PT_p0,PT_px,PT_py,PT_pz,PT_rho,PT_pmax
c
