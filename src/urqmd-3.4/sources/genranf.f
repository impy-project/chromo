
      function ranf(dummy)

      implicit none

      real*8 ranf
      integer dummy
      real rand

      ranf = dble(rand())

      return
      end

      subroutine sseed(ranseed)

      implicit none

      integer ranseed
      integer oldseed
      integer time

      save oldseed

      if (ranseed.gt.0) then
         call srand(abs(ranseed))
         return
      endif
      ranseed = time()
      if (abs(oldseed).eq.ranseed) then
         ranseed = -abs(ranseed)
         return
      else
         oldseed = ranseed
      endif
      call srand(ranseed)

      return
      end
