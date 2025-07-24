c=======================================================================
c     Below: Dummy implementations of PRNG initialization functions
c     used by EPOS. We use the Numpy PRNG so these functions are
c     replaced with dummies that do nothing.
c=======================================================================

c=======================================================================
      subroutine rmmaqd( iseed, iseq, chopt )
c-----------------------------------------------------------------------
      implicit none
      integer          iseed(3), iseq
      character        chopt*(*), cchopt*12
      end


c=======================================================================
      subroutine ranfini(seed,iseq,iqq)
c-----------------------------------------------------------------------
      implicit none
      double precision seed
      integer iseq, iqq
      end

c=======================================================================
      subroutine ranfcv(seed)
c-----------------------------------------------------------------------
      implicit none
      double precision seed
      end

c=======================================================================
      subroutine ranflim(seed)
c-----------------------------------------------------------------------
      implicit none
      double precision seed
      end
