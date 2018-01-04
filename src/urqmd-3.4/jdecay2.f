ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine nbodydec(rm)
c
c     Revision : 1.0
c
c input:  rm: Resonance mass
c
c     {\tt nbodydec} performs the decay of a resonance with mass rm in 
c     its local rest  frame into nexit particles  with  4-momenta  and 
c     masses  stored  in the array  pnew (see  comment to SUB jdecay).
c     The accessible many-body  phase-space is homogenously populated,
c     i.e. each configuration has equal probability. The theory behind
c     this  approach  can be  found in  M.M. Block and J.D. Jackson, Z.
c     Phys. C 3, 255 (1980). The original  routine is contained in  CPC
c     (Code ACGJ). It has been modified for uQMD purposes.
c     More documentation and better readability are to follow.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        implicit none
        integer j,i,imin1
        include 'newpart.f'
        real*8 rm
        real*8 p4loc(0:3), p1loc(3), ploc(3)
        real*8 M(mprt), MEFFloc(mprt), MASS
        real*8 M1,F1,M11                   


        real*8 z3, v2, ptot, costheta, z4, v1, reci1, z9, delx2, xi,
     +   z10, pp1dot,ener1,p2, p1s, z5, sintheta, phii,psin,
     +   p1sq, u, u1, ximin,ximax,delu,energy,pi,wmax,w,esys,z2,s,
     +   z8,reci,b, delxi,a,delz,ranf

        integer ntry

        LOGICAL MASSLESS                                  

        p4loc(0) = rm
        p4loc(1) = 0.0
        p4loc(2) = 0.0
        p4loc(3) = 0.0
        wmax = 1.0            
                                                     
        PI=4.*ATAN(1.)                                                          
C                                                                               
        M1=0. !Initialize M1=sum of all masses.                                 
        MASSLESS=.FALSE.                                                        
C       Read masses.                                                            
        DO I=1,nexit
          m(i) = pnew(5,i)
          M1=M1+M(I)
        ENDDO                                                                   
        IF (M1.EQ.0.) THEN  ! If massless.                                      
          MASSLESS=.TRUE.                                                    
          WMAX=1.                                                            
        ENDIF                                                                   
        MEFFloc(1)=M(1)
C       Initialize ESYS.                                                        
        ESYS=SQRT(p4loc(0)*p4loc(0)+p4loc(1)*p4loc(1)+
     +            p4loc(2)*p4loc(2)+p4loc(3)*p4loc(3))              
C                                                                               

C       Main Calculation                                                        
C                                                                               
        ntry=0
60        W=1.    !Initial weight, for each new event.                          
          ntry=ntry+1
          MASS=M1 !Initial M=total mass of ALL particles.                       
          ENERGY=p4loc(0)

          DO I=nexit-1,2,-1   !Loop over all N-2 effective masses needed.   
            U1=M(I+1)/ENERGY                                                    
            U=U1**2                                                             
            MASS=MASS-M(I+1)                                                    
C           MASS=SUM of all rest masses of the REMAINING particles.             
            XIMIN=(MASS/ENERGY)**2          !This is Xi,minimum.                
            DELU=1.-U1                                                          
            XIMAX=DELU**2           !This is Xi,maximum,where                   
                        !XIMIN <= XI <= XIMAX.                                  
            DELXI=XIMAX-XIMIN       ! DELXI=delta (XI)=XIMAX-XIMIN.             
            B=(1.+U-XIMIN)
            A=dble(I)*B            ! A=commonly used factor.                   
            IMIN1=I-1               ! IMIN1=commonly used factor.               
            DELZ=(A-dble(IMIN1)*DELXI)*DELXI**IMIN1        ! DELZ=Zmax-Zmin    
C                                                                               
C           Here, we introduce the FAST generators.                             
C                                                                               
            RECI=1./dble(I)                                                    
            S=1./(dble(I)-dble(IMIN1)*DELXI/B)
                        ! for the distribution i*(i-1)*y**(i-2)*(1-y).          
100         Z2=ranf(0)                                                      
            IF (Z2.LT.S) THEN  ! The distribution i*(i-1)*y**(i-2)*(1-y).       
                Z8=ranf(0)                                                  
                Z9=ranf(0)                                                  
                RECI1=1./dble(IMIN1)                                           
                DELX2=DELXI*Z8**RECI*Z9**RECI1                                  
              ELSE                                                              
                Z10=ranf(0)                                                 
                DELX2=DELXI*Z10**RECI   ! The probability distribution          
                                ! i*y**(i-1)                                    
            ENDIF                                                               
            IF (DELX2.EQ.0.) GOTO 100 ! Guards against division by 0.           
110         XI=XIMIN+DELX2 ! XI=XImin+deltaXI.                                  
                        ! We reweight (multiply) W by                           
                                ! DELZ*F1*[(XI/(XI-XIMIN))**(I-2)]/(1+U-XI) .   
            W=W*DELZ*F1(U,XI)*((XI/DELX2)**(I-2))/(1.+U-XI)                     
            MEFFloc(I)=ENERGY*SQRT(XI)
                                ! update E for next effective mass.             
            ENERGY=MEFFloc(I)                                                      
          ENDDO                                                                 
          V1=(M(1)/ENERGY)**2     ! Set up final weight, with particles 1 and 2.
          V2=(M(2)/ENERGY)**2                                                   
          W=W*F1(V1,V2)
                !We find WMAX, the max weight, here.                            
c          IF (W.GT.WMAX) WMAX=W   ! Update WMAX.                                
                        ! This routine selects W=1 (unweighted events).         
          Z3=ranf(0)                                                        
          IF (W.LT.WMAX*Z3.and.ntry.le.1000) THEN                                                
             GOTO 60                                                            
          ENDIF                                                                 
        ! We have accepted event, so see if we Lorentz transform it.    

        M11=p4loc(0)                                                               
        p1loc(1)=p4loc(1)                                                             
        p1loc(2)=p4loc(2)                                                             
        p1loc(3)=p4loc(3)                                                             

        !Iterate over all blob masses, MEFF(I), where MEFF(1)=M(1),MEFF(N)=E*.  
        DO 2500 I=nexit,2,-1                                                             
            ENERGY=.5*(M11+(M(I)**2-MEFFloc(I-1)**2)/M11)                          
            PTOT=SQRT(ENERGY**2-M(I)**2)                                        
                !Find RANDOM cos(theta*)=COSTHETA, random PHI*=PHI              
                ! SINTHETA=SIN(THETHA*)                                         
            Z4=ranf(0)                                                      
            COSTHETA=2.*Z4-1. ! -PI <= THETA* <= PI                             
            SINTHETA=SQRT(1.-COSTHETA**2)                                       
            Z5=ranf(0)                                                      
            PHII=2.*PI*Z5   ! 0 <= PHI* <= 2*PI, random PHII                    
            PSIN=PTOT*SINTHETA !Commonly used combination.                      
                ! Calculate momentum compon. of particle I, ploc(k), k=1 to 3.     
            ploc(1)=PSIN*COS(PHII)
            ploc(2)=PSIN*SIN(PHII)
            ploc(3)=PTOT*COSTHETA ! z-component.                                   
            P1SQ=p1loc(1)**2+p1loc(2)**2+p1loc(3)**2                                     
            P1S=SQRT(P1SQ)                                                      
            ENER1=SQRT(P1SQ+M11**2)                                             
                ! Calculate Plab(i) =                                           
                        !P*(i) + betagamma(i)*                                  
                        ! [Energy + betagamma(j).ploc(j)/(gamma+1)],               
                                !where . means DOT product, i,j=x,y,z.          
            PP1DOT=ploc(1)*p1loc(1)+ploc(2)*p1loc(2)+ploc(3)*p1loc(3)
            A=(ENERGY+PP1DOT/M11/(1.+ENER1/M11))/M11                            
                ! Plab=P1 for particle I;store in matrix OUT(K,I,3), update     
                                ! new M11 and new p1loc()=ploc()-p1loc().                
            P2=0.                                                               
            DO J=1,3                                                            
                 ploc(J)=ploc(J)+A*p1loc(J)
                 P2=P2+ploc(J)*ploc(J)
                 p1loc(J)=p1loc(J)-ploc(J)
            ENDDO                                                               
            ENERGY=SQRT(P2+M(I)**2)
            M11=MEFFloc(I-1)
c            WRITE (5,2600) M(I),ENERGY,ploc(1),ploc(2),ploc(3)                           
            pnew(5,i) = m(i)
            pnew(4,i) = sqrt(m(i)**2+ploc(1)**2+ploc(2)**2+ploc(3)**2)
            pnew(1,i) = ploc(1) 
            pnew(2,i) = ploc(2) 
            pnew(3,i) = ploc(3) 
2500    ENDDO                                                                   
2600    FORMAT(0P,F7.4,2X,G13.7,5X,G13.7,3X,G13.7,3X,G13.7)                     
        P2=0. ! Do LAST particle here.                                          
        DO J=1,3                                                                
            ploc(J)=p1loc(J)
            P2=P2+ploc(J)*ploc(J)                                                     
        ENDDO                                                                   
        ENERGY=SQRT(P2+M(1)**2)
c               WRITE (5,2600) M(1),ENERGY,ploc(1),ploc(2),ploc(3)                        
c               WRITE (5,*)                                                      
            pnew(5,1) = m(1)
            pnew(4,1) = sqrt(m(1)**2+ploc(1)**2+ploc(2)**2+ploc(3)**2)
            pnew(1,1) = ploc(1) 
            pnew(2,1) = ploc(2) 
            pnew(3,1) = ploc(3) 

        RETURN                                                                  
        END                                                                     
!-------------------------------------------------------------------------      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        FUNCTION F1(V1,V2)                                                 
c        ! Function F1(V1,V2)=SQR(1+(V1-V2)**2-2*(V1+V2))=2*(P*)/(E*).
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none
        REAL*8 F1, F2, V1, V2
        F2=1.+(V1-V2)**2-2.*(V1+V2)                                             
        IF (F2.LE.0.) THEN                                                      
             F1=0.                                                              
           ELSE                                                                 
             F1=SQRT(F2)
        ENDIF                                                                   
        END                                                                     

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function M_inv_2(v01,vx1,vy1,vz1,
     +                 v02,vx2,vy2,vz2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 M_inv_2,v01,vx1,vy1,vz1,
     +               v02,vx2,vy2,vz2 

      M_inv_2 = sqrt((v01+v02)**2
     +              -(vx1+vx2)**2
     +              -(vy1+vy2)**2
     +              -(vz1+vz2)**2)
      return 
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function M_inv_3(v01,vx1,vy1,vz1,
     +                 v02,vx2,vy2,vz2,
     +                 v03,vx3,vy3,vz3)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 M_inv_3,v01,vx1,vy1,vz1,
     +               v02,vx2,vy2,vz2, 
     +               v03,vx3,vy3,vz3 

      M_inv_3 = sqrt((v01+v02+v03)**2
     +              -(vx1+vx2+vx3)**2
     +              -(vy1+vy2+vy3)**2
     +              -(vz1+vz2+vz3)**2)
      return 
      end





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine jdecay(rm)
C        input px,py,pz : CM-momenta of total system                   
C              rm:        Mass of resonance (sqrt(s))
c for pnew and pgen : 
c      first index: 1=px, 2=py, 3=pz, 4=E, 5=m0
c      second index: particle number
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
       include 'newpart.f'
       real*8 pgen(5,mprt),rnd(mprt),u(3),beta(3),wt,tmp
       real*8 wtmax,rm,sum,pi,sum1,sum2,pcms,ranf,gamma,bp,phi,qcm,r1234

       parameter(pi=3.141592654)
       integer n,nadd1,i,j,ii,k
c
       pgen(1,1)=0.d0
       pgen(2,1)=0.d0
       pgen(3,1)=0.d0
       pgen(5,1)=rm                                                    
       pgen(4,1)=rm
c
       nadd1=nexit-1
c                                                    
       pgen(5,nexit)=pnew(5,nexit)                                       
                                                                       
c Two body decay                                                       
c ---------------                                                      
       if(nexit.eq.2) goto 400                                          
       sum=0.
c sum: sum of masses in the outgoing channel 
       do 20 n=1,nexit                                                     
          sum=sum+pnew(5,n)                                             
 20    continue                                                        
                                                                       
c     calculate maximum phase-space weight wtmax 
c     ------------------------------------ 
       wtmax=0.5                                                  
       sum1=pgen(5,1)                                             
       sum2=sum-pnew(5,1)                                         
       do 200 i=1,nadd1                                           
          wtmax=wtmax*pcms(sum1,sum2,pnew(5,i))                          
          sum1=sum1-pnew(5,i)                                             
          sum2=sum2-pnew(5,i+1)                                           
 200   continue                                                        
                                                                       
c     generate uniform nexit-body phase space                           
c     --------------------------------------                           
300   continue 
c first generate nexit random numbers with decreasing value
c as excess energy distribution weights
      rnd(1)=ranf(1)
      do 110 i=2,nexit
         rnd(i)=ranf(1)
         do 120 j=i,2,-1
            if(rnd(j).gt.rnd(j-1)) then
               tmp=rnd(j-1)
               rnd(j-1)=rnd(j)
               rnd(j)=tmp
            endif
 120      continue
 110   continue
c last weight has to be zero
      rnd(nexit)=0.                                                     
c now ?
      wt=1.                                                            
      sum1=sum                                                         
      do 330 i=2,nexit                                                  
         sum1=sum1-pnew(5,i-1)                                           
         pgen(5,i)=sum1+rnd(i)*(pgen(5,1)-sum)
         if(pgen(5,1)-sum.lt.0.0) write(6,*)'glrrrrrp'
         wt=wt*pcms(pgen(5,i-1),pgen(5,i),pnew(5,i-1))                  
 330  continue                                                         
      r1234=ranf(1)                                                   
      if(wt.lt.r1234*wtmax) goto 300                                    
      
c     carry out two-body decays in pgen frames                         
c     ----------------------------------------                         
 400  continue                                                         
      do 410 i=1,nadd1                                                 
         qcm=pcms(pgen(5,i),pgen(5,i+1),pnew(5,i))                     
c        u(3) is cos(theta)
         u(3)=2.*ranf(1)-1.                                            
         phi=2.*pi*ranf(1)                                             
         u(1)=sqrt(1.-u(3)**2)*cos(phi)                                 
         u(2)=sqrt(1.-u(3)**2)*sin(phi)                                 
         do 420 j=1,3                                                   
            pnew(j,i)=qcm*u(j)                                           
            pgen(j,i+1)=-pnew(j,i)                                       
 420     continue                                                       
         pnew(4,i)=sqrt(qcm**2+pnew(5,i)**2)                            
         pgen(4,i+1)=sqrt(qcm**2+pgen(5,i+1)**2)                        
 410  continue                                                         
      do 430 j=1,4                                                     
         pnew(j,nexit)=pgen(j,nexit)                                      
 430  continue                                                         
      
c     boost pgen frames to lab frame                                   
c     -------------------------------------------------                
      do 500 ii=1,nadd1                                                
         i=nexit-ii                                                      
         do 510 j=1,3                                                   
            beta(j)=pgen(j,i)/pgen(4,i)                                  
 510     continue                                                       
         gamma=pgen(4,i)/pgen(5,i)                                      
         do 520 k=i,nexit                                                
            bp=beta(1)*pnew(1,k)+beta(2)*pnew(2,k)+beta(3)*pnew(3,k)     
            do 530 j=1,3                                                 
               pnew(j,k)=pnew(j,k)+gamma*beta(j)*(pnew(4,k)               
     &              +bp*gamma/(gamma+1.))                                      
 530        continue                                                     
            pnew(4,k)=gamma*(pnew(4,k)+bp)                               
 520     continue                                                       
 500  continue                                                         

      return                                                           
      end                                                              
                                                                       
