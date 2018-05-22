c $Id: alpharanf.f,v 1.2 1999/01/18 09:56:50 ernst Exp $
C*****************   R A N F   *******************************
      real*8 FUNCTION RANF(ix)
      integer ix
      real*8 rand
      external rand
       ranf=dble(rand())
c       write(7,*) ranf
      RETURN
      end
c
c
      subroutine SSEED(ranseed)
      integer ranseed,oldseed,time
      external time,srand
      save

      itt=0


      if(ranseed.gt.0) then
         call srand(abs(ranseed))
         return
      endif
      ranseed=time()
      if(abs(oldseed).eq.ranseed) then
         ranseed=-1*abs(ranseed)
         return
      else
         oldseed=ranseed
      endif
      call srand(ranseed)
      RETURN
      END
