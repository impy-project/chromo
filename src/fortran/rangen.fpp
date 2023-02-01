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
      double precision function simrnd()
c-----------------------------------------------------------------------
c  alternative interface to npyrng
c-----------------------------------------------------------------------
      call npyrng(simrnd)
      end

c=======================================================================
      double precision function gasdev( dummy )
c-----------------------------------------------------------------------
c  used by SIBYLL-2.3x
c-----------------------------------------------------------------------
      implicit none
      external npygas
      integer dummy
      integer*8 bitgen
      common /npy/bitgen
      call npygas(gasdev, bitgen)
      end

c=======================================================================
      real function spgasdev( dummy )
c-----------------------------------------------------------------------
c  used by SIBYLL-2.1
c-----------------------------------------------------------------------
      implicit none
      external npygas
      integer dummy
      integer*8 bitgen
      common /npy/bitgen
      double precision rval
      call npygas(rval, bitgen)
      spgasdev = real(rval)
      end

c=======================================================================
      subroutine rmmard( rvec,lenv,iseq )
c-----------------------------------------------------------------------
c  used by epos
c-----------------------------------------------------------------------
      implicit none
      double precision rvec(*), rval
      integer iseq,lenv,ivec
      do ivec = 1, lenv
        call npyrng(rval)
        rvec(ivec) = rval
      enddo
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
      double precision function dranf(dummy)
c-----------------------------------------------------------------------
c  used by epos
c-----------------------------------------------------------------------
      implicit none
      double precision dummy
      call npyrng(dranf)
      end

c=======================================================================
c  sibyll random generator
c-----------------------------------------------------------------------
#ifdef SIBYLL_21
      real function s_rndm(dummy)
#else
      double precision function s_rndm(dummy)
#endif 
c-----------------------------------------------------------------------
      implicit none
      integer dummy
      double precision simrnd
#ifdef SIBYLL_21
555   s_rndm = real(simrnd())
      if ((s_rndm.le.0e0).or.(s_rndm.ge.1e0)) goto 555     
#else
      s_rndm = simrnd()
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
