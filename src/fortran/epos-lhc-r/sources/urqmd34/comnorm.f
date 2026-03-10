c $Id: comnorm.f,v 1.5 1999/01/18 09:56:55 ernst Exp $
      integer n
      parameter (n = 800)
      real*8 x_norm(0:3,1:n),y_norm(0:3,1:n)
      real*8 dx
      common /normsplin/ x_norm,y_norm,dx

