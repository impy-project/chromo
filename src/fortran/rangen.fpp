c=======================================================================
      subroutine npyrng( rval )
c-----------------------------------------------------------------------
c  interface to C code
c-----------------------------------------------------------------------
      implicit none
      external npynxt
      double precision rval
      integer*8 bitgen
      common /npy/bitgen
      call npynxt(rval, bitgen)
      end

c=======================================================================
      subroutine ranfgt(seed)
c-----------------------------------------------------------------------
c     used by EPOS
c     Stash the state of the random number generator
c     interface to C code
c-----------------------------------------------------------------------           
      implicit none
      external npynxt_get_state
      integer*8 bitgen
      integer*8 state_arr(4)
      double precision seed
      common /npy/ bitgen
      common /npystash/ state_arr
c     seed is dummy for compiler
      seed = 0       
      call npynxt_get_state(state_arr, bitgen)   
      end


c=======================================================================
      subroutine ranfst(seed)
c-----------------------------------------------------------------------
c     used by EPOS
c     Restore the state of the random number generator from stash
c     interface to C code
c-----------------------------------------------------------------------           
      implicit none
      external npynxt_set_state
      integer*8 bitgen
      integer*8 state_arr(4)
      double precision seed
      common /npy/ bitgen
      common /npystash/ state_arr
c     seed is dummy for compiler     
      seed = 0     
      call npynxt_set_state(state_arr, bitgen)  
      end         

c=======================================================================
      function rangen()
c-----------------------------------------------------------------------
c  used by EPOS
c-----------------------------------------------------------------------
      double precision rval
 20   call npyrng(rval)
      if(rval.le.0.or.rval.ge.1) goto 20
      rangen=sngl(rval)
      end

c=======================================================================
      function drangen( dummy )
c-----------------------------------------------------------------------
c  used by EPOS
c-----------------------------------------------------------------------
      double precision drangen, dummy
      call npyrng(drangen)
      end

c=======================================================================
      double precision function simrnd()
c-----------------------------------------------------------------------
c  alternative interface to npyrng
c-----------------------------------------------------------------------
      call npyrng(simrnd)
      end

c=======================================================================
      function gasdev( dummy )
c-----------------------------------------------------------------------
c  used by SIBYLL-2.3x
c-----------------------------------------------------------------------
      implicit none
      external npygas
      double precision gasdev
      integer dummy
      integer*8 bitgen
      common /npy/bitgen
      call npygas(gasdev, bitgen)
      end

c=======================================================================
      function spgasdev( dummy )
c-----------------------------------------------------------------------
c  used by SIBYLL-2.1
c-----------------------------------------------------------------------
      implicit none
      external npygas
      real spgasdev
      integer dummy
      integer*8 bitgen
      common /npy/bitgen
      double precision rval
      call npygas(rval, bitgen)
      spgasdev = sngl(rval)
      end

c=======================================================================
      subroutine rm48( rvec,lenv )
c-----------------------------------------------------------------------
c  r(ando)m (number generator for evaporation module/dpmjet)
c
c  this subroutine is called from routines of evaporation module.
c  arguments:
c   rvec   = doubl.prec. vector field to be filled with random numbers
c   lenv   = length of vector (# of randnumbers to be generated)
c-----------------------------------------------------------------------
      implicit none
      double precision rvec(*), rval
      integer lenv,ivec
      do ivec = 1, lenv
        call npyrng(rval)
        rvec(ivec) = rval
      enddo
      end

c=======================================================================
      subroutine rmmard( rvec,lenv,iseq )
c-----------------------------------------------------------------------
c  used by epos
c-----------------------------------------------------------------------
      implicit none
      double precision rvec(*)
      integer lenv,iseq
      call rm48(rvec,lenv)
      end

c=======================================================================
      function dranf(dummy)
c-----------------------------------------------------------------------
c  used by epos
c-----------------------------------------------------------------------
      implicit none
      double precision dranf,dummy
      call npyrng(dranf)
      end

c=======================================================================
c  sibyll random generator
c-----------------------------------------------------------------------
      function s_rndm(dummy)
c-----------------------------------------------------------------------
      implicit none
      integer dummy
#ifdef SIBYLL_21
      real s_rndm
      double precision rval
555   call npyrng(rval)
      if ((rval.le.0e0).or.(rval.ge.1e0)) goto 555
      s_rndm = sngl(rval)
#else
      double precision s_rndm
      call npyrng(s_rndm)
#endif
      end

c=======================================================================
      double precision function pyr()
c-----------------------------------------------------------------------
c  pythia random generator
c-----------------------------------------------------------------------
      implicit none
      call npyrng(pyr)
      end

c=======================================================================
      double precision function rndm()
c-----------------------------------------------------------------------
c  random generator for dpmjet
c-----------------------------------------------------------------------
      implicit none
      call npyrng(rndm)
      end

c=======================================================================
      double precision function psran()
c-----------------------------------------------------------------------
c  random generator for qgsjet
c-----------------------------------------------------------------------
      implicit none
      call npyrng(psran)
      end

c=======================================================================
      double precision function ranf()
c-----------------------------------------------------------------------
c  random generator for urqmd
c-----------------------------------------------------------------------
      implicit none
      call npyrng(ranf)
      end

c=======================================================================
      double precision function rlu()
c-----------------------------------------------------------------------
c  random generator for jetset
c-----------------------------------------------------------------------
      implicit none
      call npyrng(rlu)
      end

c=======================================================================
      double precision function dt_rndm(vdummy)
c-----------------------------------------------------------------------
c  this functon is called from dpm_jet306 routines.
c-----------------------------------------------------------------------
      implicit none
      double precision vdummy
      call npyrng(dt_rndm)
      end
