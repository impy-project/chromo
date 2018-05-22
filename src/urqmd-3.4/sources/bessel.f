C  (C) Copr. 1986-92 Numerical Recipes Software.

       real*8 function bessi0(x)
       real*8 x
c Returns the modified bessel function I0(x) for any real x
       real*8 ax
       double precision p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,
     &  q8,q9,y
       save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
       data p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     &   1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
       data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     &  0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,
     &  0.2635537d-1,-0.1647633d-1,0.392377d-2/
       if(abs(x).lt.3.75d0)then
         y=(x/3.75d0)**2
         bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
       else
        ax=abs(x)
        y=3.75d0/ax
        bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4
     &     +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))        
       endif
       return 
       end


       real*8 FUNCTION bessk0(x)
   
       double precision x
C USES bessi0
C Returns the modified Bessel function K0(x) for positive real x.
       REAL*8 bessi0
       DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,
     * q2,q3,q4,q5,q6,q7,y

       SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7

       DATA p1,p2,p3,p4,p5,p6,p7/-.57721566d0,.42278420d0,.23069756d0,
     * .3488590d-1,.262698d-2,.10750d-3,.74d-5/

       DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-.7832358d-1,.2189568d-1,
     * -.1062446d-1,.587872d-2,-.251540d-2,.53208d-3/

       if (x.le.2.0d0) then
        y=x*x/4.0d0
        bessk0=(-log(x/2.0d0)*bessi0(x))+(p1+y*(p2+y*(p3+
     * y*(p4+y*(p5+y*(p6+y*p7))))))
       else
        y=(2.0d0/x)
       bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+
     * y*(q4+y*(q5+y*(q6+y*q7))))))

       endif
       return
       END

        

       real*8 FUNCTION bessi1(x)
   
       double precision x

       REAL*8 ax
       DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,
     * q8,q9,y

       SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
       DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     * 0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
       DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     * -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,
     * -0.2895312d-1,0.1787654d-1,-0.420059d-2/
       if (abs(x).lt.3.75d0) then
        y=(x/3.75d0)**2
        bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
       else
        ax=abs(x)
        y=3.75d0/ax
        bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+
     * y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
       if(x.lt.0.0d0)bessi1=-bessi1
       endif
       return
       END
        
       real*8 function bessk1(x)
       real*8 x
c     uses bessi1
c     returns the modified Bessel function Ki(x) for positive real x
       real*8 bessi1
       double precision p1,p2,p3,p4,p5,p6,p7,q1,
     &  q2,q3,q4,q5,q6,q7,y
       save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
       data p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     &  -0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
       data q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498618d0,
     &  -0.3655620d-1,
     &  0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
       if (x.le.2.0d0) then
          y=x*x/4.0d0
          bessk1=(log(x/2.0d0)*bessi1(x))+(1.0d0/x)*(p1+y*(p2+
     &      y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
       else
          y=2.0d0/x
          bessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+
     &      y*(q6+y*q7))))))
       endif
       return
       end


       real*8 FUNCTION bessk(n,x)
       INTEGER n
       REAL*8 x
C      USES bessk0,bessk1
c      Returns the modified Bessel function Kn(x) for positive x and n >= 2.
       INTEGER j
       REAL*8 bk,bkm,bkp,tox,bessk0,bessk1

       if (n.lt.2.0d0) then
        write(6,*) 'n < 2 in bessel.f'
       end if 
       tox=2.0d0/x

       bkm=bessk0(x)
       bk=bessk1(x)

       do j=1,n-1
        bkp=bkm+j*tox*bk
        bkm=bk
        bk=bkp
       end do
       bessk=bk
       return
       END

C  (C) Copr. 1986-92 Numerical Recipes Software.
