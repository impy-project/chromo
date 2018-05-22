c $Id: iso.f,v 1.7 1999/01/18 09:57:06 ernst Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ISOCGK4(J1,M1,J2,M2,Jnew,Mnew,ITAG)
c
c     Revision : 1.0
c
c     This subroutine determines according to probabilities given by
c     Clebsch Gordan cefficients the total and 3-component of the
c     isospin of up to 4 outgoing particles
C
c     input:  
c     (isocgk: two part. in and two part. out) 
c              J1    : 2*I  of ingoing particle 1
c              M1    : 2*I3 of ingoing particle 1
c              J2    : 2*I  of ingoing particle 2
c              M2    : 2*I3 of ingoing particle 2
c              Jnew  : 2*I  of outgoing particles (array)
c              
c
c     (isonew: one part. in and two part. out) 
c              J     : 2*I  of ingoing particle
c              M     : 2*I3 of ingoing particle
c              Jnew  : 2*I  of outgoing particles (array)
c              ITAG=-50 << necessary for correct functioning of routine
c
c     input/output:
c              Mnew  : 2*I3 of outgoing particles (array)
c                      Mnew(i)=-9 to determine the I3 component of the
c                          i-th particle randomly
c              ITAG  : = -1 then no possible isospin combination has been found
c
c
c     function calls:
c                     ranf()
c                     clebsch
c
c important global variables:
c              nexit : number of outgoing particles
c
c important local variables:
C     JMINOL/JMINNW  THE MINIMAL POSSIBLE TOTAL ISOSPIN IN IN-/OUT-STATE
C     M       TOTAL I3 OF IN-/OUT-STATE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'newpart.f'
      integer m1,j1,m2,j2,itag,Jnew,Mnew
      integer Jtot,M,Jminol,Jmaxol,Jminnw,Jmaxnw,i,j,k,l,il,jp,jpr
      integer jmin,jmax,Nj,ifind
      integer m1pr,m1p,m1pos,m1out
      integer m2pr,m2p,m2pos,m2out
      integer m3pr,m3p,m3pos,m3out
      integer m4pr,m4p,m4pos,m4out,Mmin,Mmax
      integer m12,m34,j12,j34
      integer Jin,Min
      real*8 pjin,prbout,prbsum,zrand,c12,c34,c_tot
      DIMENSION pjin(20),prbout(20,20,20,20)
      DIMENSION m1out(20),m2out(20),m3out(20),m4out(20)
      DIMENSION JNEW(mprt),Mnew(mprt),Mmin(mprt),Mmax(mprt)
      real*8 ranf,clebsch

      ITAG=0
      M=M1+M2


c 1.) first treat  some special cases
        if(nexit.gt.4)then
          write(6,*)'ISOCGK: only <=4 outgoing Isospins can be coupled'
          itag=-1
          stop 137
          return
        endif

        if(nexit.eq.3)Jnew(4)=0

      if(nexit.eq.2)then       !nexit.eq.2
        Jnew(3)=0
        Jnew(4)=0
c check for zero in out-channel:
       if(Jnew(1).eq.0) then
          Mnew(1)=0
          Mnew(2)=m
          return
       elseif(Jnew(2).eq.0) then
          Mnew(2)=0
          Mnew(1)=m
          return
       endif
      endif                   !nexit.eq.2

c 2.) determine possible min and max isospins for in/out state
c
c determine number of possible in-states
       Jminol=  MAX0(IABS(J1-J2),IABS(M))
       Jmaxol= J1+J2

c determine number of possible out-states
c JMINNW=  MAX0(IABS(J1NEW-J2NEW),IABS(M))
        Jminnw=1000
        do 1 i=-1,1,2
         do 2 j=-1,1,2
          do 3 k=-1,1,2
           do 4 l=-1,1,2
            jp=IABS(i*Jnew(1)+j*Jnew(2)+k*Jnew(3)+l*Jnew(4))
            if(jp.lt.Jminnw)Jminnw=jp
4          continue
3         continue
2        continue
1       continue
        Jminnw=MAX0(Jminnw,IABS(M))

c JMAXNW= J1NEW+J2NEW
        Jmaxnw=0
         do 5 i=1,nexit
          Jmaxnw=Jmaxnw+Jnew(i)
5       continue

c  check which possible states match (are common for in AND out state)
       Jmin  =  MAX0(Jminol,Jminnw)
       Jmax  =  MIN0(Jmaxol,Jmaxnw)
c error check for unphysical input
       if(Jmin.gt.Jmax) then
          itag=-1
          write(6,*)'isocgk: jmin > jmax : unphysical input!'
          write(6,*) J1,M1,J2,M2,jnew(1),jnew(2),jmin,jmax
          return
       endif

c 3.) calculate number of possible isospins
       nj = (Jmax-Jmin)/2 +1
       if(J1.eq.0.or.J2.eq.0)then
          if(Jmin.ne.Jmax) then
             itag=-1
             write(6,*) 'J1(2)=0,Jmin.ne.Jmax IN ISOCGK - check calling'
             return
          endif
          if(J1.eq.0.and.J2.eq.0) then
             write(6,*) "J1,J2=0 IN ISOCGK - can't couple this"
             itag=-1
               return
          endif
c here only one total isospin is possible (probability is unity)
          pjin(1)=1.
          goto 310
       END IF      !J1.EQ.0.OR.J2.EQ.0
c if no overlap between in and out state, return with itag=-1
       if(nj.le.0) then
          itag=-1
          write(6,*)'Isocgk: nj.le.0 - no combination possible'
          return
       endif

       ifind=0
c 4) loops over all possible combinations of J1,J2,M1,M2,Jtot
c    to get the probabilities of the in-channel couplings
             DO 6 jpr=Jmin,Jmax,2
                      ifind=ifind+1
                      pjin(ifind)=clebsch(J1,J2,m1,m2,jpr)
6          CONTINUE

C     
c error message, if not all possible Jtot's have been found
c       if(ifind.ne.nj) then
c          write(6,*)'ERROR IN ISOCGK IFIND.NE.NJ'
c          stop 137
c       endif
c sum CGKs over all possible Jtots (-> probabilities)
       prbsum=0.
       do 7 il=1,nj
          prbsum=prbsum+pjin(il)
7      continue

c check for nonsense
       IF(prbsum.le.0.) THEN
          write(6,*)'ERROR IN ISOCGK 30:PRBSUM.LE.0.'
          stop 137
       END IF
c normalize PJIN(.) to 1 
c now PJIN contains CGK-based probabilities for the different possible Jtots
       Do 8 il=1,nj
          pjin(il)=pjin(il)/prbsum
8      Continue
310   continue
c 5) now throw dice to determine one of the possible Jtots
       zrand=ranf(0)
       Do 9 il=1,nj
          if(zrand.lt.pjin(il))then
c this is now the "real Jtot"
             Jtot= Jmin +2*(il-1)
               goto 11
          else
             zrand=zrand-pjin(il)
          endif
9      Continue
11       Continue


       goto 111

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c this is the entry for one in and >= two out particles
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ENTRY ISONEW4(JIN,MIN,JNEW,MNEW,ITAG)

       IF(ITAG.EQ.-50) THEN

          Jtot=Jin
          M=Min
          itag=0

       END IF                   !itag=-50

c here now both cases (one/two-in particles) together
c now the in-channel is determined -> get out-channel
 111   continue

c special cases
       if(nexit.gt.4)then
          write(6,*)'ISONEW: only <=4 outgoing Isospins can be coupled'
          itag=-1
          stop 137
          return
       endif
       
       if(nexit.eq.3)then
          Jnew(4)=0
          Mnew(4)=0
       endif

       if(nexit.eq.2)then
          Jnew(3)=0
          Jnew(4)=0
          Mnew(3)=0
          Mnew(4)=0
c check for zero in out-channel:
          if(Jnew(1).eq.0) then
             mnew(1)=0
             mnew(2)=m
             return
          elseif(Jnew(2).eq.0) then
             mnew(2)=0
             mnew(1)=m
             return
          endif
       endif                    !nexit=2


c reset counters
        m1pos=2*Jnew(1)+1
        m2pos=2*Jnew(2)+1
        m3pos=2*Jnew(3)+1
        m4pos=2*Jnew(4)+1
       Do 161 m1p=1,20
          m1out(m1p)=0
          m2out(m1p)=0
          m3out(m1p)=0
          m4out(m1p)=0
          Do 162 m2p=1,m2pos
           Do 163 m3p=1,m3pos
            Do 164 m4p=1,m4pos
           prbout(m1p,m2p,m3p,m4p)=0d0
164         Continue
163      Continue
162     Continue
161    Continue

c get min/maximal M
        Do 100 i=1,4
            Mmin(i)=-Jnew(i)
            Mmax(i)=Jnew(i)
100     Continue

c calculate the possible |J_i M_i> combinations and their probability
        do 112 j12=Iabs(Jnew(1)-Jnew(2)),(Jnew(1)+Jnew(2)),2
        do 134 j34=Iabs(Jnew(3)-Jnew(4)),(Jnew(3)+Jnew(4)),2
        m1pos=0
        Do 41 m1pr=Mmin(1),Mmax(1),2
         m1pos=m1pos+1
         m1out(m1pos)=m1pr
           m2pos=0
         Do 42 m2pr=Mmin(2),Mmax(2),2
          m2pos=m2pos+1
          m2out(m2pos)=m2pr
            m3pos=0
          Do 43 m3pr=Mmin(3),Mmax(3),2
           m3pos=m3pos+1
           m3out(m3pos)=m3pr
             m4pos=0
             Do 44 m4pr=Mmin(4),Mmax(4),2
              m4pos=m4pos+1
              m4out(m4pos)=m4pr
c            m4pr=m-m1pr-m2pr-m3pr
              If(m1pr+m2pr+m3pr+m4pr.ne.m)goto 44
              m12=m1pr+m2pr
              m34=m3pr+m4pr

                c12=clebsch(Jnew(1),Jnew(2),m1pr,m2pr,J12)
                c34=clebsch(Jnew(3),Jnew(4),m3pr,m4pr,J34)
                c_tot= clebsch(J12,J34,m12,m34,Jtot)

                prbout(m1pos,m2pos,m3pos,m4pos)=
     +         prbout(m1pos,m2pos,m3pos,m4pos)+c12*c34*c_tot

                
 44          Continue
 43       Continue
 42    Continue
 41   Continue
 134  continue
 112  continue

c error check
       if(m1pos.eq.0.or.m2pos.eq.0.or.m3pos.eq.0.or.m4pos.eq.0)then
          write(6,*)'IN ISOCGK/ISONEW: MPOS=0 ERROR'
            write(6,*)"Can't couple Jin, Min=",Jtot,M
            write(6,*)'To J1,J2,J3,J4=',Jnew(1),Jnew(2),Jnew(3),Jnew(4)
          itag=-1
          return
       endif

c sum up all CGKs
       prbsum=0.
       Do 51 m1p=1,m1pos
          Do 52 m2p=1,m2pos
           Do 53 m3p=1,m3pos 
            Do 54 m4p=1,m4pos
           prbsum=prbsum+prbout(m1p,m2p,m3p,m4p)
54          Continue
53       continue
52      continue
51     continue

c error check
       IF(prbsum.le.0.) then
          write(6,*)'ERROR IN ISOCGK/ISONEW:PRBSUM.LE.0.'
            write(6,*)"Can't couple Jin, Min=",Jtot,M
            write(6,*)'To J1,J2,J3,J4=',Jnew(1),Jnew(2),Jnew(3),Jnew(4)
          stop 137
       endif

c normalize to 1 (now we have real probabilities for different Mout combis)
       Do 61 m1p=1,m1pos
          Do 62 m2p=1,m2pos
           Do 63 m3p=1,m3pos 
            Do 64 m4p=1,m4pos
           prbout(m1p,m2p,m3p,m4p)=prbout(m1p,m2p,m3p,m4p)/prbsum
c      write(*,*)'!p= ',m1out(m1p),m2out(m2p),m3out(m3p),
c     & m4out(m4p),prbout(m1p,m2p,m3p,m4p)
64          Continue
63       Continue
62      Continue
61     Continue

c now determine according to the PRBOUT values the outgoing M combination
       zrand=ranf(0)
       Do 71 m1p=1,m1pos
          Do 72 m2p=1,m2pos
           Do 73 m3p=1,m3pos 
            Do 74 m4p=1,m4pos
           if(zrand.lt.prbout(m1p,m2p,m3p,m4p)) then
             Mnew(1)= M1out(m1p)
             Mnew(2)= M2out(m2p)
             Mnew(3)= M3out(m3p)
               Mnew(4)= M4out(m4p)
               goto 70
           else
             zrand=zrand-prbout(m1p,m2p,m3p,m4p)
           endif
74          Continue
73       Continue
72      Continue
71     Continue
70       Continue
       RETURN
       END


C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function fcgk(i1,i2,iz1,iz2,i)
c returns the normalized clebsch gorden factor also for combinations 
c involving strange mesons and antibaryons
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      include 'comres.f'
      real*8 c
      integer i1,i2,iz1,iz2,i,iz,isoit,i12,i12a,ir,strit,icnt
      logical nombbb

      fcgk=0D0
      icnt=0
      c=0d0
      iz=iz1+iz2
      if(isoit(i).lt.iabs(iz))goto 1008
      if(isoit(i1)*isoit(i2).eq.0)then
         c=1d0
         goto 1008
      end if

      call cgknrm(isoit(i),iz,isoit(i1),isoit(i2),iz1,iz2,ir,c) 
      if(i1.eq.i2.and.iz1.ne.iz2)c=2d0*c
c... particle exchange 

      if(ir.ne.0)then
         icnt=icnt+1
         if(icnt.le.1)then 
           write(6,*)'fcgk: no iso-spin decomposition found for:',
     @       i,iz,' to ',i1,iz1,'+',i2,iz2
           write(6,*)'      please check this channel'
        end if
        return  
      end if

      if(strit(i).eq.0)then
c... this is now for particle+antiparticle (except nonstrange mesons)
         i12=i1*i2
c... the charge conjugated states have the same weight
csab: new condition, the old one violated antibaryon-meson symmetry
         if(i12.lt.0.and.min(iabs(i1),iabs(i2)).gt.maxbar)then
c... for example anti-K* + K 
            if(i1.ne.-i2)c=c*5d-1
         end if
      end if
1008   fcgk=c  
      return
      end
C####C##1#########2#########3#########4#########5#########6#########7##
       subroutine cgknrm(JIN,MIN,J1NEW,J2NEW,M1IN,M2IN,ierr,cf)
C gives the normalized cg-factor i.e. poosibility into a given
C iso-spin decomposition of JIN,MIN into J1NEW,J2NEW,M1IN,M2IN
C ierr equals 0 if there is any alowed J1,J2,M1,M2 (not necessaryly
C equal to J1NEW,J2NEW,M1IN,M2). 
C ierr is not equal 0 if all channels are iso-spin forbidden  
C for specific couplings possibly involving strange particles or
C anti-particles function fcgk should be used (see beyond)  
Coutput cf :  normalized cg-factor
C####C##1#########2#########3#########4#########5#########6#########7##
       implicit integer (i - n)
       implicit real*8 (a - h , o - z)
       DIMENSION PRBOUT(20),M1OUT(20)
       real*8 clebsch

c the in particle of course defines Jtot and Mtot 
       cf =0d0
       ierr = 0
 1     JTOT=JIN
       M=MIN
       ITAG=0
c check for zero in out-channel:
          if(j1new.eq.0) then
             m1new=0
             m2new=m
             return
          elseif(j2new.eq.0) then
             m2new=0
             m1new=m
             return
          endif

c here now both cases (one/two in particles) together
c reset counters
       M1POS=0
c loop over all J1,J2,Jtot,M1,M2 combinations
       do 39 m1pr=-j1new,j1new,2
          m2pr=m-m1pr
c if J1new and J2new and Jtot and Mtot create a match then store M1new
c inM1OUT array and the CGK in PRBOUT array; the counter for possible
c Mnew combinations is M1POS 
          M1POS=M1POS+1
          M1OUT(M1POS)=M1PR
          PRBOUT(M1POS)=clebsch(j1new,j2new,m1pr,m2pr,jtot)
          if( ! (m1pr.eq.m2in.and.m2pr.eq.m1in).or.
     @       (m2pr.eq.m2in.and.m1pr.eq.m1in))cf=cf+PRBOUT(M1POS)
 39    continue
c error check
       IF(M1POS.EQ.0) then
          write(6,*)'IN ISOCGK: M1POS=0 ERROR'
          write(6,*)'jtot,j1new,j2new,m',jtot,j1new,j2new,m
          itag=-1
          return
       endif
c sum over all CGKs
       PRBSUM=0.
       DO 50 M1P=1,M1POS
          PRBSUM=PRBSUM+PRBOUT(M1P)
 50    continue

c error check
       IF(PRBSUM.LT.1d-3) then
          ierr = 1
          cf = 0d0 
          return
       endif
c normalize to 1 (now we have real probabilities for different Mout combis)
       cf=cf/PRBSUM

       RETURN
       END

          
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function dbweight(j1,m1,j2,m2,j1new,j2new)
c
c     Revision : 1.0
c
c     This function delivers a weight, based on a Clebsch Gordan
c     coefficient, for detailed balance cross sections
C
c     input:  
c              J1    : 2*I  of ingoing particle 1
c              J2    : 2*I  of ingoing particle 2
c              M1    : 2*I3 of ingoing particle 1
c              M2    : 2*I3 of ingoing particle 2
c              J1new : 2*I  of outgoing particle 1
c              J2new : 2*I  of outgoing particle 2
c
c     output:
c              weight : weight factor
c
c     function calls:
c                     clebsch
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer j1,m1,j2,m2,j1new,j2new,jminol,jmaxol,jminnw,jmaxnw
      integer nj,jpr,m,jmax,jmin,ind
      real*8 clebsch,weight(10)

      dbweight=0.d0
      m=m1+m2
c determine number of possible states
c fist the in state
      JMINOL=  MAX0(IABS(J1-J2),IABS(M))
      JMAXOL= J1+J2
c now the out state
      JMINNW=  MAX0(IABS(J1NEW-J2NEW),IABS(M))
      JMAXNW= J1NEW+J2NEW
c  which possible states match (are common for in AND out state)
      JMIN  =  MAX0(JMINOL,JMINNW)
      JMAX  =  MIN0(JMAXOL,JMAXNW)
      NJ = (JMAX-JMIN)/2 +1
      if(nj.lt.1) return
      ind=0
      do 18 jpr=jmin,jmax,2
         ind=ind+1
         weight(ind)=clebsch(j1,j2,m1,m2,jpr)
         dbweight=dbweight+weight(ind)
 18   continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function clebsch(j1,j2,m1,m2,j3)
c
c     Revision : 1.0 
c
c     This function delivers a Clebsch Gordan Coefficient, which has been
c     calculated by w3j.
C
c     input:  
c              J1    : 2*I  of ingoing particle 1
c              J2    : 2*I  of ingoing particle 2
c              M1    : 2*I3 of ingoing particle 1
c              M2    : 2*I3 of ingoing particle 2
c              J3    : 2*I  of projection requested 
c                           (i.e. resonance to be formed)
c
c     output:
c              clebsch : CGK**2
c
c     function calls:
c                     w3j
c                     !!! first call function after intialization with loginit
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      real*8 LogFak( 0 : 100 ),dj1,dj2,dj3,dm1,dm2,dm3,w3j
      real*8 cfct,cgk
      integer j1,j2,j3,m1,m2,ipot,jmm,jmm1
      common /FACTORIALS/LogFak

      integer jmax
      parameter(jmax=7) 
      real*8 cgktab(0:jmax,0:jmax,-jmax:jmax,-jmax:jmax,0:jmax)
      common /cgks/cgktab

c     each cgk for j's in the range up to jmax is calculated only once
c     and then stored in the cgktab table for further use
      jmm1=max(j1,j2)
      jmm=max(j3,jmm1)
      if(jmm.gt.jmax.or.(cgktab(j1,j2,m1,m2,j3).lt.-8.d0)) then

         dj1=dble(j1)/2.d0
         dj2=dble(j2)/2.d0
         dj3=dble(j3)/2.d0
         dm1=dble(m1)/2.d0
         dm2=dble(m2)/2.d0
         dm3=-(dm1+dm2)
         ipot=(j1+m1+j2-m2)/2
         cfct=sqrt(2*dj3+1.d0)/(-(1.d0**ipot))
         cgk=cfct*w3j(dj1,dj2,dj3,dm1,dm2,dm3)
         clebsch=cgk**2
         if(jmm.le.jmax) then
            cgktab(j1,j2,m1,m2,j3)=clebsch
         endif
      else
         clebsch=cgktab(j1,j2,m1,m2,j3)
      endif
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function W3j( J1, J2, J3, M1, M2, M3 )
c
c
c  This program calculates the 3-j wigner symbols according to the
c  representation of A. Lindner. 
c
c  Reference:
c  A. Lindner, Drehimpulse in der Quantenmechanik, Teubner 1984, P.39
c
c====================================================================== 
      implicit none
c  Program input

      real*8       J1, J2, J3, M1, M2, M3

c  Program returns

      real*8       W3j

c  Global variables

      real*8       LogFak( 0 : 100 )
      common /FACTORIALS/   LogFak

c  Program variables

      real*8       R1, R2, R3,  R4, R5, R6, R7, R8, R9
      real*8       N( 1 : 3, 1 : 3 )
      real*8       Sum1, Sum2
      real*8       Sigma
      real*8       LF_R1, LF_R2, LF_R3, LF_R4, LF_R5, LF_R6
      real*8       LF_R7, LF_R8, LF_R9
      real*8       LF_Sigma
      real*8       Hlp1, Hlp2, Pre
      real*8       Summe, S( 0 : 100 )
      integer*4    Signum
      integer*4    in
      real*8       minimal
      integer*4    imin, jmin
      integer*4    i, j
      real*8       dn
      real*8       dummy

c  Start of calculation
                                                                                
c  Evaluation due to equivalence with Regge symbol

c     call LogInit
C
      Sigma = J1 + J2 + J3
                                                                                
      N( 1, 1 ) = -J1 + J2 + J3
      N( 1, 2 ) =  J1 - J2 + J3
      N( 1, 3 ) =  J1 + J2 - J3
      N( 2, 1 ) =  J1 - M1
      N( 2, 2 ) =  J2 - M2
      N( 2, 3 ) =  J3 - M3
      N( 3, 1 ) =  J1 + M1
      N( 3, 2 ) =  J2 + M2
      N( 3, 3 ) =  J3 + M3

      do 20 i = 1, 3
         do 10 j = 1, 3
            if ( nint( N( i, j ) ) .lt. 0 ) goto 99999
 10      continue
         Sum1 = N( i, 1 ) + N( i, 2 ) + N( i, 3 )
         Sum2 = N( 1, i ) + N( 2, i ) + N( 3, i )
         if ( nint( Sum1 ) .ne. nint( Sigma ) )   goto 99999
         if ( nint( Sum2 ) .ne. nint( Sigma ) )   goto 99999
 20   continue

c      do 101 i=1, 3
c         write(6,'(3f14.5)') (N(i,j), j=1, 3 )
c 101  continue

      imin = 1
      jmin = 1
      Signum = 1
      minimal = N( 1, 1 ) 

c  Looking for the smallest N( i, j )

      do 40 i = 1, 3
         do 30 j = 1, 3
            if ( N(i,j) .lt. minimal )   then
               minimal = N( i, j )
               imin = i
               jmin = j 
            endif
 30      continue
 40   continue

      Signum = 1

      if ( imin .gt. 1 )   then 
         do 50 j = 1, 3 
            dummy = N( 1, j )
            N( 1, j ) = N( imin, j )
            N( imin, j ) = dummy
 50      continue
         Signum = (-1)**nint( Sigma )
      endif 

      if ( jmin .gt. 1 )   then
         do 60 i = 1, 3
            dummy = N( i, 1 )
            N( i, 1 ) = N( i, jmin )
            N( i, jmin ) = dummy
 60      continue
         Signum = (-1)**nint( Sigma ) * Signum
      endif 


      R1 = N( 1, 1 ) 
      R2 = N( 1, 2 )
      R3 = N( 1, 3 )
      R4 = N( 2, 1 )
      R5 = N( 2, 2 )
      R6 = N( 2, 3 )
      R7 = N( 3, 1 )
      R8 = N( 3, 2 )
      R9 = N( 3, 3 )

      LF_R1 = LogFak( nint( R1 ) ) 
      LF_R2 = LogFak( nint( R2 ) )
      LF_R3 = LogFak( nint( R3 ) )
      LF_R4 = LogFak( nint( R4 ) ) 
      LF_R5 = LogFak( nint( R5 ) )
      LF_R6 = LogFak( nint( R6 ) )
      LF_R7 = LogFak( nint( R7 ) )
      LF_R8 = LogFak( nint( R8 ) )
      LF_R9 = LogFak( nint( R9 ) ) 
      LF_Sigma = LogFak( nint( Sigma+1.d0 ) )

      Hlp1 = ( LF_R2 + LF_R3 + LF_R4 + LF_R7 - LF_Sigma - 
     &         LF_R1 - LF_R5 - LF_R9 - LF_R6 - LF_R8 ) / 2.d0

      Pre = dexp( Hlp1 ) * (-1)**( nint( R6 + R8 ) )

      Hlp2 = LF_R6 - LogFak( nint(R6-R1) )
     &       + LF_R8 - LogFak( nint(R8-R1) )
      S( 0 ) = dexp( Hlp2 )
      Summe = S( 0 )

      do 70 in = 1, nint( R1 )
         dn = dble( in ) 
         S( in ) = (-1)*S( in-1 ) * ( R1+1.d0-dn ) * ( R5+1.d0-dn )
     &        * ( R9+1.d0-dn ) / dn / ( R6-R1+dn ) / ( R8-R1+dn )
         Summe = Summe + S( in )
 70   continue

      W3j = Pre * Summe * Signum
      return

99999 W3j = 0.d0
      return

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine LogInit
c=====================================================================
c
c  This function computes the logarithm of the factorials and
c  stores it in the array LogFak.
c
c=====================================================================

      implicit none

c  Program output
      integer jmax
      parameter(jmax=7)
      real*8 cgktab(0:jmax,0:jmax,-jmax:jmax,-jmax:jmax,0:jmax)
      common /cgks/cgktab

      real*8      LogFak( 0 : 100 )

      common / FACTORIALS /   LogFak

c  Program variables

      integer*4   i,j1,j2,j3,m1,m2

c  Program start
      do 1 j1=0,jmax
         do 1 j2=0,jmax
            do 1 m1=-jmax,jmax
               do 1 m2=-jmax,jmax
                  do 1 j3=0,jmax
                     cgktab(j1,j2,m1,m2,j3)=-9.d0
 1                continue


      LogFak( 0 ) = 0.d0 
      do 10 i = 1, 100
         LogFak( i ) = LogFak( i-1 ) + dlog( dble( i ) )  
 10   continue

      end
































