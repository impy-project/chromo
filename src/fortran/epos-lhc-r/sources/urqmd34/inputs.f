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

c
c maximum mass for projectile and target
      integer AAmax
      parameter(AAmax=300)
c storage arrays for projectile and target nuclei
      integer PT_iso3(AAmax,2),PT_ityp(AAmax,2),PT_spin(AAmax,2)
      integer PT_charge(AAmax,2),PT_AA(2),PT_uid(AAmax,2)
      real*8 PT_r0(AAmax,2),PT_rx(AAmax,2),PT_ry(AAmax,2),PT_rz(AAmax,2)
      real*8 PT_p0(AAmax,2),PT_px(AAmax,2),PT_py(AAmax,2),PT_pz(AAmax,2)
      real*8 PT_dectime(AAmax,2),PT_fmass(AAmax,2),PT_rho(AAmax,2)
      real*8 PT_pmax(AAmax,2)
c common blocks (data transfer between cascinit and getinit)
      common/ProTarInts/PT_iso3,PT_ityp,PT_spin,PT_charge,PT_AA,PT_uid
      common/ProTarReals/PT_r0,PT_rx,PT_ry,PT_rz,PT_fmass,PT_dectime,
     &                   PT_p0,PT_px,PT_py,PT_pz,PT_rho,PT_pmax
c
