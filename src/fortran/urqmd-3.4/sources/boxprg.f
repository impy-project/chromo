c $Id: boxprg.f,v 1.8 1999/01/18 09:56:53 ernst Exp $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bptinit(ibox)
c
c     Unit     : Initis all the particles setted by the bpt command
c     Version  : 1.0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      
      include 'coms.f'
      include 'comres.f'
      include 'boxinc.f'
      include 'options.f'

c Var
c       counter, spin 
      integer i,ibox,fchg,getspin
c       randomnumbergenerator, particlemass, deacy times
      real*8 ranf,massit,dectim
c       momentum, angle distribution
      real*8 P,cost,sint,phi    

      ecm=0.d0
      ebeam=0.d0
        
c  main program

c loop over all particles
      do 42 i=npart+1,npart+bptpart(ibox)
c          configuration space 
         r0(i)=0.d0     
         rx(i)=lboxhalbe*(1-2*ranf(0)) 
         ry(i)=lboxhalbe*(1-2*ranf(0)) 
         rz(i)=lboxhalbe*(1-2*ranf(0))                  
 
c          isospin and ityp
         iso3(i)=bptiso3(ibox)
         ityp(i)=bptityp(ibox)
         uid_cnt = uid_cnt + 1
         uid(i)=uid_cnt

c set baryon and   meson numbers                
         if(abs(ityp(i)).le.maxbar) then
            nbar=nbar+1
         else
            nmes=nmes+1
         endif  
         
c               charge
         charge(i)=fchg(iso3(i),ityp(i))
c               massarray
         fmass(i)=massit(ityp(i))
c               Spin
         spin(i)=getspin(ityp(i),-1)
c               decaytime
         dectime(i)=dectim(i,1)
        
 42   continue
c       End of loop

      if (edensflag.le.0) then
c       homogenious momentum distribution, randomly distributed
c       max momentum is a parameter
         do 45 i=npart+1,npart+bptpart(ibox)
            P=bptpmax(ibox)*ranf(0)**(1./3.)
            cost = 1.-2.*ranf(0)
            sint = sqrt(1.-cost**2)
            phi = 2.*Pi*ranf(0)
            px(i) = P*sint*cos(phi)
            py(i) = P*sint*sin(phi)
            pz(i) = P*cost
            call setonshell(i)
 45      continue
      elseif (edensflag.ge.1) then
           
c       energiedensity

c loop over all particles
         do 60 i=npart+1,npart+bptpart(ibox)
            P=bptpmax(ibox)/bptpart(ibox)*ranf(0)**(1./3.)
              
            cost = 1.-2.*ranf(0)
            sint = sqrt(1.-cost**2)
            phi = 2.*Pi*ranf(0)
 
c different initialisations
            
            if (para.eq.0) then
                           
c Boxmode                                               
               if (i.eq.1) write(*,*) 'Boxmode'
               px(i) = P*sint*cos(phi)
               py(i) = P*sint*sin(phi)
               pz(i) = P*cost
                 
            elseif(para.eq.1) then

c stream over stream 
               if (i.eq.1) write(*,*) 'streammode'      
               px(i) = 0.d0
               py(i) = 0.d0
               pz(i) = bptpmax(ibox)/bptpart(ibox)*(-1.d0)**i
            elseif(para.eq.2) then

c slab on slab
               if (i.eq.1) write(*,*) 'slabmode'
               px(i)=0.d0
               py(i)=0.d0
               if (rz(i).gt.0) then
                  pz(i)=(-1.0d0)*bptpmax(ibox)/bptpart(ibox)
               else
                  pz(i)=bptpmax(ibox)/bptpart(ibox)
               endif
            endif
 60      continue
      endif
c       sum of particles
      npart=npart+bptpart(ibox)
      Write(*,*) 'Particles = ',npart
cc
      return

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function swapi(x,dx)
c
c Version: 1.0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real*8 swapi, x, dx  

      swapi = x
 1    if (swapi.lt.-dx) then
         swapi = swapi + 2.0d0*dx
         goto 1
      end if  
 2    if (swapi.gt.dx) then
         swapi = swapi - 2.0d0*dx
         goto 2
      end if  
      end     


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Function Energie(alpha,max)
c
c     Unit     : calculate the energy
c     Version  : 1.0 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      
      include 'coms.f'
      integer i
      
      real*8 alpha,max
      real*8 Energie, E

      E=0
      Do 42 i=1,npart
         E=E+sqrt((alpha**2)*(px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))+
     $        fmass(i)**2)
 42   continue

      Energie=E-max     
  
      Return
      End



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Function Regula(me)
c
c     Unit     : Searches for the zero of the function
c     Version  : 1.0 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      implicit none
        
      include 'coms.f'
        
      real*8 Regula,xn,xu,x0,me,Energie
      integer i
      real*8 E1,E2,E3

c init
      xu=0.d0
      x0=3.d0     
c main
      Write(*,*) 'Regula is running!'
      i=0
 10   Continue
      i=i+1
      E1=Energie(x0,me)        
      E2=Energie(xu,me) 
        
      xn=x0-(E1*(x0-xu))/(E1-E2)
      E3=Energie(xn,me)
      IF ((E2*E3).LE.0) then
         x0=xn
      else
         xu=xn
      EndIF

      IF ((ABS(x0-xu).GE.1.D-12).and.(i.le.1000).and.(
     &     ((E3.ge.1.D-12).or.(-E3.ge.1.D-12)))) goto 10 
      Regula=xn 
      End


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function wallcoll(ipart,wall)
c
c     Unit     : Collisions with an imaginary wall
c     Version  : 1.0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      
      include 'coms.f'
      include 'boxinc.f'
      include 'options.f'
      
c var
      real*8  betax,betay,betaz
      real*8  ty,tz
      real*8    tn
      integer wall,ipart
      integer wally,wallz
        
        
c Mainprogram
c init the variables
      wall=0
      tn=0

c velocity 
      betax=px(ipart)/p0(ipart)
      betay=py(ipart)/p0(ipart)
      betaz=pz(ipart)/p0(ipart)

c check which wall is reached next and sort by impact time
c wall presents the wall and tn the time

      if (betax.lt.0) then
         wall=-4
         tn=(-lboxhalbe-rx(ipart))/(-max(-betax,1.d-13))
      else
         wall=-1
         tn=((lboxhalbe-rx(ipart))/max(betax,1.d-13))
      endif
        
      if (betay.lt.0) then
         wally=-5
         ty=(-lboxhalbe-ry(ipart))/(-max(-betay,1.d-13))
      else
         wally=-2
         ty=((lboxhalbe-ry(ipart))/max(betay,1.d-13))
      endif
        
      if(ty.lt.tn) then 
         tn=ty
         wall=wally
      endif
      
      if (betaz.lt.0) then
         wallz=-6
         tz=(-lboxhalbe-rz(ipart))/(-max(-betaz,1.d-13))
      else
         wallz=-3
         tz=((lboxhalbe-rz(ipart))/max(betaz,1.d-13))
      endif
        
      if(tz.lt.tn) then
         tn=tz
         wall=wallz
      endif


c sets the time for the earliest wall collision
      wallcoll=tn
      
      return
      End


