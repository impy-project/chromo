c $Id: paulibl.f,v 1.4 1999/01/18 09:57:11 ernst Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       logical function paulibl(i,rhob,pseudorap)
c
c     Revision : 1.0
c
c     This function determines wether the final state of particle i 
c     is pauli-blocked ({\tt .true.}) or not ({\tt .false.}).
c     The baryon-density at the location of particle i is returned in {\tt rhob}
c
c                     
cinput   i :     Index of particle to be added
c
coutput rhob :   baryon density at location of i
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       implicit none
       integer i, j
       real*8 afit, bfit, rho, f, rhob, test, r2, p2, pgw, rhob0
       parameter (afit=1.49641d0, bfit=0.208736d0)

       real*8 qij(4),pij(4),pr2,q2,qp

c.. variables used for density calculation

       real*8 sigma,rhobcf,rhobrs,betaxj,betayj,betazj,norm
       parameter (sigma=1.5d0)
       real*8 vx,vy,vz,vz2,gammaz,g,b
       real*8 betax,betay,betaz,j0cf,j1cf,j2cf,j3cf
       real*8 beta2,gboost,j0rs,j1rs,j2rs,j3rs  
       real*8 pseudorap,pseudorapj,momtotj,ptransj,diff

       include 'coms.f'
       include 'comres.f'
       include 'options.f'

       rho = 0.0d0
       rhob = 0.0d0
       rhob0 = 0d0
       f = 0.0d0
       pgw = 1.0d0/hqc/hqc/gw

c..   if Lorentz-contraction of projectile and target is enabled
       if(CTOption(23).eq.0)goto 3

       do 1 j=1,nbar
          r2 = (rx(i)-rx(j))**2+(ry(i)-ry(j))**2+(rz(i)-rz(j))**2
          if ((ityp(i).eq.ityp(j)).and.(iso3(i).eq.iso3(j))) then
             p2 = (px(i)+ffermpx(i)-px(j)-ffermpx(j))**2
     &           +(py(i)+ffermpy(i)-py(j)-ffermpy(j))**2
     &           +(pz(i)+ffermpz(i)-pz(j)-ffermpz(j))**2
             p2 = 0.25d0*p2
             rho = rho + dexp(-(2.0d0*gw*r2))
             f = f + dexp(-(gw*r2)-pgw*p2)
          end if
          rhob = rhob + dexp(-(2.0d0*gw*r2))
 1     continue
       paulibl = .true.
       test = afit + bfit*rho
       if (test.gt.f) paulibl = .false.
       if (CTOption(10).eq.1) paulibl=.false.

       rhob = rhob*(2.0d0*gw/pi)**1.5/rho0 
                

       return

 3     continue

c  first the pauliblocking
       do 108 j=1,nbar
          qij(4)=r0(i)-r0(j)
          qij(1)=rx(i)-rx(j)
          qij(2)=ry(i)-ry(j)
          qij(3)=rz(i)-rz(j)
          pij(4)=p0(i)+p0(j)
          pij(1)=px(i)+ffermpx(i)+px(j)+ffermpx(j)
          pij(2)=py(i)+ffermpy(i)+py(j)+ffermpy(j)
          pij(3)=pz(i)+ffermpz(i)+pz(j)+ffermpz(j)  
          q2=qij(4)**2-qij(1)**2-qij(2)**2-qij(3)**2
          p2=pij(4)**2-pij(1)**2-pij(2)**2-pij(3)**2
          qp=qij(4)*pij(4)-qij(1)*pij(1)-qij(2)*pij(2)
     .       -qij(3)*pij(3)
          r2=qp**2/p2 - q2
          if(r2.lt.0) then
             write(6,*)'***(E) negative transverse distance !!',r2
             write(6,*)r0(j),rx(j),ry(j),rz(j),p0(j),px(j),py(j),pz(j)
             write(6,*)r0(i),rx(i),ry(i),rz(i),p0(i),px(i),py(i),pz(i)
             r2=1000.d0
          endif
          pr2 = (px(i)-px(j))**2+(py(i)-py(j))**2+(pz(i)-pz(j))**2
          if ((ityp(i).eq.ityp(j)).and.(iso3(i).eq.iso3(j))) then
             rho=rho+dexp(-(2.0d0*gw*r2))
             f=f+dexp(-(gw*r2)-.25d0*pgw*pr2)
          end if
c baryon density in rest frame of particle
c         if(j.ne.lstcoll(i))rhob=rhob+dexp(-(2.0d0*gw*r2))
 108  continue
      paulibl=.true.
      test=afit+bfit*rho
      if (test.gt.f) paulibl=.false.
      if (CTOption(10).eq.1) paulibl=.false.

c      rhob=max(0d0,min(.1d3,rhob*(2.0d0*gw/pi)**1.5/rho0)) 

      if(CTOption(46).eq.1) goto 4
c calculate the baryon density via the four-current
c.. put the necessary variables to zero

        rhobcf=0.0d0
        rhobrs=0.0d0
        betaxj=0.0d0
        betayj=0.0d0
        betazj=0.0d0
        norm=0.0d0
        betax=0.0d0
        betay=0.0d0
        betaz=0.0d0
        
c loop over all baryons   
      do 109 j=1,nbar     

               
c calculate the velocities       
        vx=(px(j)+ffermpx(j))/p0(j)
        vy=(py(j)+ffermpy(j))/p0(j)
        vz=(pz(j)+ffermpz(j))/p0(j)
    
        vz2=vz**2
      
c calculate the gamma for the Lorentz-contraction of the Gaussian
       if((1.0d0-vz2).gt.1.0d-16)then
        gammaz=1.0d0/(sqrt(1.0d0-vz2))
       else
        gammaz=1.0d0
       end if  
c the density in the computational frame
c addition of brayons and subtraction of antibaryons

c gaussian distribution including sign of the charge     
       b=ityp(j)/abs(ityp(j))
       g=b*gammaz*exp(-(((rx(j)-rx(i))**2+(ry(j)-ry(i))**2
     &    +(rz(j)-rz(i))**2*gammaz**2)/(2.0d0*sigma**2)))    
        
      rhobcf=rhobcf+(1.0d0/(2.0d0*pi))**1.5*(1.0d0/sigma)**3*g
              
c mean values of velocities    
        betaxj=betaxj+vx*g
        betayj=betayj+vy*g
        betazj=betazj+vz*g
c corresponding normalization
        norm=norm+abs(g)
        
 109  continue   
c velocity of the cell around i         
       if(norm.gt.1.d-16)then
        betax=betaxj/norm
        betay=betayj/norm
        betaz=betazj/norm
       else
        rhob=0.0d0
        return
       end if
          
c the components of jmu in the computational frame 
        j0cf=rhobcf
        j1cf=rhobcf*betax
        j2cf=rhobcf*betay
        j3cf=rhobcf*betaz
     
        beta2=betax**2+betay**2+betaz**2
      if((1.0d0-beta2).gt.1.0d-16)then 
       gboost=1.0d0/(sqrt(1.0d0-beta2))
      else
       gboost=1.0d0
      end if        

c transform jmu to the frame where j-vector is zero
        j0rs=gboost*j0cf-betax*gboost*j1cf
     &    -betay*gboost*j2cf-betaz*gboost*j3cf
        j1rs=(-betax)*gboost*j0cf
     &    +(1.0d0+(gboost-1.0d0)*betax**2/beta2)*j1cf
     &    +(gboost-1.0d0)*betax*betay/beta2*j2cf
     &    +(gboost-1.0d0)*betax*betaz/beta2*j3cf
        j2rs=(-betay)*gboost*j0cf
     &    +(gboost-1.0d0)*betay*betax/beta2*j1cf
     &    +(1.0d0+(gboost-1.0d0)*betay**2/beta2)*j2cf
     &    +(gboost-1.0d0)*betay*betaz/beta2*j3cf
        j3rs=(-betaz)*gboost*j0cf
     &    +(gboost-1.0d0)*betaz*betax/beta2*j1cf
     &    +(gboost-1.0d0)*betaz*betay/beta2*j2cf
     &    +(1.0d0+(gboost-1.0d0)*betaz**2/beta2)*j3cf

c the zero component is the baryon density at position i                   
      rhob=j0rs
      
      return

c calculate quark density instead of baryon density
 4    continue
c calculate the quark density via the four-current
c.. put the necessary variables to zero

        rhobcf=0.0d0
        rhobrs=0.0d0
        betaxj=0.0d0
        betayj=0.0d0
        betazj=0.0d0
        norm=0.0d0
        betax=0.0d0
        betay=0.0d0
        betaz=0.0d0
        
c loop over all particles   
        do 110 j=1,nbar+nmes        
           momtotj=dsqrt((pz(j)+ffermpz(j))**2+(py(j)+ffermpy(j))**2
     $          +(px(j)+ffermpx(j))**2)
           ptransj=dsqrt((py(j)+ffermpy(j))**2
     $          +(px(j)+ffermpx(j))**2)
           pseudorapj=0.5d0*dlog((momtotj+ptransj)/(momtotj-ptransj))     
           diff=dabs(pseudorapj-pseudorap)
c           write(*,*)j, pseudorap, diff
           if((diff.le.CTParam(70)).or.(pseudorap.lt.-999.d0))then
 
c calculate the velocities       
        vx=(px(j)+ffermpx(j))/p0(j)
        vy=(py(j)+ffermpy(j))/p0(j)
        vz=(pz(j)+ffermpz(j))/p0(j)
    
        vz2=vz**2
      
c calculate the gamma for the Lorentz-contraction of the Gaussian
       if((1.0d0-vz2).gt.1.0d-16)then
        gammaz=1.0d0/(sqrt(1.0d0-vz2))
       else
        gammaz=1.0d0
       end if  
c the density in the computational frame
c addition of valence quarks and antiquarks        
       if(abs(ityp(j)).lt.100) then 
          b=3.0d0
       else
          b=2.0d0             
       end if

      g=b*gammaz*exp(-(((rx(j)-rx(i))**2+(ry(j)-ry(i))**2
     &    +(rz(j)-rz(i))**2*gammaz**2)/(2*sigma**2)))   

      rhobcf=rhobcf+(1.0d0/(2.0d0*pi))**1.5*(1.0d0/sigma)**3*g
    
c mean values of velocities  
        
        betaxj=betaxj+vx*g
        betayj=betayj+vy*g
        betazj=betazj+vz*g
       
c corresponding normalization
        norm=norm+g
       endif 
 110  continue   
c velocity of the cell around i         
        betax=betaxj/norm
        betay=betayj/norm
        betaz=betazj/norm
c the components of jmu in the computational frame 
        j0cf=rhobcf
        j1cf=rhobcf*betax
        j2cf=rhobcf*betay
        j3cf=rhobcf*betaz
     
        beta2=betax**2+betay**2+betaz**2
      if((1.0d0-beta2).gt.1.0d-16)then 
       gboost=1.0d0/(sqrt(1.0d0-beta2))
      else
       gboost=1.0d0
      end if        

c transform jmu to the frame where j-vector is zero
        j0rs=gboost*j0cf-betax*gboost*j1cf
     &    -betay*gboost*j2cf-betaz*gboost*j3cf
        j1rs=(-betax)*gboost*j0cf
     &    +(1.0d0+(gboost-1.0d0)*betax**2/beta2)*j1cf
     &    +(gboost-1.0d0)*betax*betay/beta2*j2cf
     &    +(gboost-1.0d0)*betax*betaz/beta2*j3cf
        j2rs=(-betay)*gboost*j0cf
     &    +(gboost-1.0d0)*betay*betax/beta2*j1cf
     &    +(1.0d0+(gboost-1.0d0)*betay**2/beta2)*j2cf
     &    +(gboost-1.0d0)*betay*betaz/beta2*j3cf
        j3rs=(-betaz)*gboost*j0cf
     &    +(gboost-1.0d0)*betaz*betax/beta2*j1cf
     &    +(gboost-1.0d0)*betaz*betay/beta2*j2cf
     &    +(1.0d0+(gboost-1.0d0)*betaz**2/beta2)*j3cf

c the zero component is the quark density at position i                    
      rhob=j0rs




      return
      
      end




