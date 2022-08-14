cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine error (function_name, message, value, level)
c
c     Revision : 1.0
c
c     output of errors and warnings
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      character function_name*(*)
      character message*(*)
      real*8 value
      integer level

      include 'inputs.f'

      integer errdev

      errdev=6

      if ((level.lt.1).or.(level.gt.3)) then
         write (errdev,*) '*** Error in Subroutine error:'
         write (errdev,*) '*** Message: Wrong Errorlevel'
         write (errdev,*) '*** Related Value: ', level
      endif

      if (level.eq.1) then 
         write (errdev,*) '*** Warning in Subroutine ',function_name,':'
      elseif (level.eq.2) then
         write (errdev,*) '*** Error in Subroutine ',function_name,':'
      else
         write (errdev,*) '*** Fatal Error in Subroutine ',
     $        function_name,':'
      endif
      write (errdev,*) '*** Message: ',message
      write (errdev,*) '*** Related Value: ',value 

      if (level.ge.3) then
         write (errdev,*)
         write (errdev,*) '*** Program stopped.'
         stop 137
      endif

      return
      end





