

      subroutine hydro(thydro_start,thydro)

      implicit none


      include 'coms.f'
      include 'freezeout.f'
      INCLUDE 'defs.f'
      include 'options.f'

      real*8 mx(ngr,ngr,ngr),my(ngr,ngr,ngr),mz(ngr,ngr,ngr)
      real*8 elab(ngr,ngr,ngr),r(ngr,ngr,ngr)
      real*8 dx,dt_hydro,t,vol,cut
      real*8 thydro,thydro_start,starttime
      real*8 e0,n0,Bag,hc3,pi2
      integer n
      real*4 ifr(ngr*ngr*ngr),jfr(ngr*ngr*ngr),kfr(ngr*ngr*ngr)
      integer mmax,m
      real*4 tf(ngr*ngr*ngr),tfr(ngr*ngr*ngr),muqfr(ngr*ngr*ngr)
      real*4 musfr(ngr*ngr*ngr)
      real*8 vyfr(ngr*ngr*ngr),vxfr(ngr*ngr*ngr)
      real*8 vzfr(ngr*ngr*ngr),gfr(ngr*ngr*ngr)
      real*8 dstfr(ngr*ngr*ngr),dsxfr(ngr*ngr*ngr)
      real*8 dsyfr(ngr*ngr*ngr),dszfr(ngr*ngr*ngr)       

c hydro common blocks
      common /grstate/ e0,n0,Bag,hc3,pi2
      common /gitter/ dx,dt_hydro,t,vol,cut
      common /onefgrid/elab,mx,my,mz,r
c common block necessary for gradual freeze-out
      common /start/ starttime 
      common /line/ n

      integer j,k,testflag,ithydro
c.. variables for final decays
      integer stidx, CTOsave,i
      logical isstable
      real*8 xdummy
c
c     Define constants. (for 1f hydro)
c     
      pi2 = dacos(-1d0)**2
      hc3 = 197.327053d0**3d0 
      e0 = 146.51751415742d0 
      n0 = 0.15891d0 
      Bag = 235d0**4 
      n=0

c.... IMPORTANT HYDRO CELL SIZE !!!!!!

      dx=CTParam(61)

      if(Ap.ne.At)then
            write(*,*) 'cannot use hydro for asymmetric systems'
            stop 137
      end if

      call f15hyin(thydro_start)      
c.. this subroutine puts the urqmd particles on the grid      
      call urqmdtohydro
      
      testflag=CTOption(48)
      starttime=thydro_start


c.. hydro call      
      write(*,*) 'starting hydro'

c      write (*,*)
c     $   "#############################################################"
c      write (*,*)
c     $   "##     The SHASTA algorithm by D. Rischke is used          ##"
c      write (*,*)
c     $   "##     for the 3d-hydro evolution                          ##"
c      write (*,*)
c     $   "#############################################################"
       
      
      call onefluid(testflag,ithydro,thydro_start,ifr,jfr,kfr,tf,dstfr,
     &   dsxfr,dsyfr,dszfr,tfr,muqfr,musfr,vxfr,vyfr,vzfr,mmax)
      write (*,'(x,a20,x,f6.3," fm")') 'hydro finished after',t-dt_hydro
    
      if(ithydro.eq.1.and.testflag.eq.0)then
       thydro=0.0d0
       write(*,*) 'hydro not called'
       return
      endif 

      thydro=t-dt_hydro

c.. perform the freezeout      
c      call cooperfrye(thydro_start)
      call coopertime(ifr,jfr,kfr,tf,dstfr,dsxfr,dsyfr,dszfr,
     &           tfr,muqfr,musfr,vxfr,vyfr,vzfr,mmax)
      call arraysort(thydro_start)

c.. write time information of hydro evolution to f15
      call f15outhy(thydro_start,thydro)
c.. copy the particles back into urqmd  
      call inithydro
      call f15hyout(thydro_start,thydro)

c.. write ouput directly after hydro freeze-out      
      if(CTOption(50).eq.1)then
c.. perform final decays
        if(CTOption(18).eq.0) then
c no do-loop is used because npart changes in loop-structure
            i=0
c disable Pauli-Blocker for final decays
            CTOsave=CTOption(10)
            CTOption(10)=1
c decay loop structure starts here
 40         continue
            i=i+1
c is particle unstable
            if(dectime(i).lt.1.d30) then
 41            continue
               isstable = .false.
               do 44 stidx=1,nstable
                  if (ityp(i).eq.stabvec(stidx)) then
c     write (6,*) 'no decay of particle ',ityp(i)
                     isstable = .true.
                  endif
 44            enddo
               if (.not.isstable) then
c     perform decay
                  call scatter(i,0,0.d0,fmass(i),xdummy)
c     backtracing if decay-product is unstable itself
                  if(dectime(i).lt.1.d30) goto 41
               endif
            endif
c     check next particle
            if(i.lt.npart) goto 40
         endif ! final decay
         CTOption(10)=CTOsave
c     final output
         call osc_event
         call file14out(2)
      end if
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine urqmdtohydro
      implicit none


      include 'coms.f'
      INCLUDE 'defs.f'
      include 'freezeout.f'    
      include 'options.f'  

      real*8 mx(ngr,ngr,ngr),my(ngr,ngr,ngr),mz(ngr,ngr,ngr)
      real*8 elab(ngr,ngr,ngr),r(ngr,ngr,ngr)
      real*8 dx,dt_hydro,t,vol,cut,ecor
      real*8 h_etot,h_btot,enorm
      integer mmax


      common /gitter/ dx,dt_hydro,t,vol,cut
      common /onefgrid/elab,mx,my,mz,r
      common /grstate/ e0,n0,Bag,hc3,pi2
      
      common /h_etot/ h_etot,h_btot

c  This program initializes the hydrodynamic fields using urqmd output.
c  This program operates in MeV like the hydrocode

c.. particle x^mu and p^mu from urqmd
      real*8 tpos,xpos,ypos,zpos,energy,phx,phy,phz
c.. baryon number
      real*8 barnum
c.. ityp, cell positions and counting integers      
      integer ptype,i,j,k,line,maxline,index
c.. grid positions      
      real*8 xposgr,yposgr,zposgr
c.. counting variables to reduce computing power
      integer ir,jr,kr,d,irr,jrr,krr,ien,tail,numspeco
c.. parameters for the gaussians      
      real*8 sigma, norm,rq
c.. total energies and momenta      
      real*8 etotu,etothy,pxtotu,pxtothy,pytotu,pytothy
      real*8 pztotu,pztothy,bartot,bartotu
c.. gamma in z-direction      
      real*8 gamma_gauss,momtot,ptrans,pseudorap

c.. grid parameter
      real*8 pi2,n0,e0,Bag,hc3,enerspec

c.. rapidty cut for spectators
      real*8 rapypar,rapyzpar,bdum,quarkdens,paulibl
      integer oldcto
c.. spectator array
      real*8 specr0(nmax), specrx(nmax), specry(nmax), specrz(nmax),
     +     specp0(nmax), specpx(nmax), specpy(nmax), specpz(nmax),
     +     specfmass(nmax), spectform(nmax),specxtotfac(nmax),
     $     specdectime(nmax)
      integer specncoll(nmax),speccharge(nmax),specityp(nmax),
     &     speclstcoll(nmax),speciso3(nmax),specorigin(nmax)
      common /speccoor/specr0,specrx,specry,specrz,specp0, specpx, 
     &     specpy,specpz,specfmass,spectform,specxtotfac,
     &     specncoll,speccharge,specityp,speclstcoll,speciso3,specorigin
     $     ,specdectime

      real*8 specffermpx(nmax), specffermpy(nmax), specffermpz(nmax)
      common /specffermi/ specffermpx, specffermpy, specffermpz

      real*8 specfrr0(nmax), specfrrx(nmax), specfrry(nmax),
     &       specfrrz(nmax),specfrp0(nmax), specfrpx(nmax), 
     &       specfrpy(nmax), specfrpz(nmax)

      common /specfrcoor/ specfrr0, specfrrx, specfrry, specfrrz,
     &                    specfrp0, specfrpx, specfrpy, specfrpz        

      integer numspec,Zspec,Aspec,strangebe,strit
      common /numspec/ numspec,Zspec,Aspec,strangebe

chp.. make Gaussian width optional
      sigma=CTParam(68)
      

c.. initialization of variables 
      do i=1,ngr
       do j=1,ngr
        do k=1,ngr
         elab(i,j,k)=0
         mx(i,j,k)=0
         my(i,j,k)=0
         mz(i,j,k)=0
         r(i,j,k)=0
        end do
       end do
      end do
      
      etotu=0
      etothy=0
      pxtotu=0
      pxtothy=0
      pytotu=0
      pytothy=0
      pztotu=0
      pztothy=0
      bartotu=0
      bartot=0
      j=0      
      Zspec=0  
      numspec=0
      strangebe=0
      Aspec=0
      enerspec=0.0d0
      numspeco=0
      if(CTOption(49).eq.0)then
c.. reset spectator arrays
       do j=1,nmax
          
        specr0(j)=0.0d0
        specrx(j)=0.0d0
        specry(j)=0.0d0
        specrz(j)=0.0d0
        specp0(j)=0.0d0
        specpx(j)=0.0d0
        specpy(j)=0.0d0
        specpz(j)=0.0d0
        specfmass(j)=0.0d0
        specityp(j)=0
        speciso3(j)=0
        speccharge(j)=0
        speclstcoll(j)=0
        specncoll(j)=0
        specorigin(j)=0
        specdectime(j)=0.0d0
        spectform(j)=0.0d0
        specxtotfac(j)=1.0d0
c.. get the freezeout coordinates right
        specfrr0(j) = 0.0d0
        specfrrx(j) = 0.0d0
        specfrry(j) = 0.0d0
        specfrrz(j) = 0.0d0
        specfrp0(j) = 0.0d0
        specfrpx(j) = 0.0d0
        specfrpy(j) = 0.0d0
        specfrpz(j) = 0.0d0
        specffermpx(j)=0.0d0
        specffermpy(j)=0.0d0
        specffermpz(j)=0.0d0
        
       end do   

       j=0
c.. sort out the spectators
       do i=1,npart
          if(ncoll(i).gt.0)then
             rapyzpar=0.5d0*dabs(dlog((p0(i)+pz(i))
     $            /(p0(i)-pz(i))))
             
             rapypar=0.5d0*dabs(dlog((p0(i)+dsqrt(pz(i)*pz(i)+
     $            px(i)*px(i)+py(i)*py(i)))
     $            /(p0(i)-dsqrt(pz(i)*pz(i)+
     $            px(i)*px(i)+py(i)*py(i)))))
          else
             rapypar=1.0d10
          endif 

c    Calculate the local density
            quarkdens=10.0d0
           
          if(CTParam(69).gt.0.001)then
             momtot=dsqrt((pz(i)+ffermpz(i))**2+(py(i)+ffermpy(i))**2
     $            +(px(i)+ffermpx(i))**2)
             ptrans=dsqrt((py(i)+ffermpy(i))**2
     $            +(px(i)+ffermpx(i))**2)
             pseudorap=0.5d0*dlog((momtot+ptrans)/(momtot-ptrans))     
             oldcto=CTOption(46)
             CTOption(46)=1
             bdum=paulibl(i,quarkdens,pseudorap)
             quarkdens=quarkdens/0.16d0/3.0d0
             CToption(46)=oldcto
             
          endif 

           if ((ncoll(i).lt.1).or.(rapyzpar.gt.CTParam(66))
     $        .or.(quarkdens.lt.ctparam(69)))then
             j=j+1
             if(ncoll(i).lt.1)numspec=numspec+1
             enerspec=enerspec+p0(i)
             specr0(j)=r0(i)
             specrx(j)=rx(i)
             specry(j)=ry(i)
             specrz(j)=rz(i)
             specp0(j)=p0(i)
             specpx(j)=px(i)
             specpy(j)=py(i)
             specpz(j)=pz(i)
             specfmass(j)=fmass(i)
             specityp(j)=ityp(i)
             speciso3(j)=iso3(i)
             speccharge(j)=charge(i)
             speclstcoll(j)=lstcoll(i)
             specncoll(j)=ncoll(i)
             specorigin(j)=origin(i)
             specdectime(j)=dectime(i)
             spectform(j)=tform(i)
             specxtotfac(j)=xtotfac(i)
c..   get the freezeout coordinates right
             specfrr0(j) = frr0(i)
             specfrrx(j) = frrx(i)
             specfrry(j) = frry(i)
             specfrrz(j) = frrz(i)
             specfrp0(j) = frp0(i)
             specfrpx(j) = frpx(i)
             specfrpy(j) = frpy(i)
             specfrpz(j) = frpz(i)
             specffermpx(j)=ffermpx(i)
             specffermpy(j)=ffermpy(i)
             specffermpz(j)=ffermpz(i)
             Zspec=Zspec+charge(i)
             if(abs(ityp(i)).lt.100)then
                Aspec=Aspec+ityp(i)/abs(ityp(i))
             endif
             strangebe=strangebe+strit(ityp(i))
c             write(25,*)rx(i),ry(i),rz(i),rapyzpar,rapypar,ityp(i),0
          end if   
       end do   
       write(*,*)'number of spectators: ',numspec
       numspec=j
        
       
c.. total energy and momenta      
       do ien=1,npart

          rapyzpar=0.5d0*dabs(dlog((p0(ien)+pz(ien))
     $         /(p0(ien)-pz(ien))))            

          rapypar=0.5d0*dabs(dlog((p0(ien)+dsqrt(pz(ien)*pz(ien)+
     $         px(ien)*px(ien)+py(ien)*py(ien)))
     $         /(p0(ien)-dsqrt(pz(ien)*pz(ien)+
     $         px(ien)*px(ien)+py(ien)*py(ien)))))
           quarkdens=10.0d0

          if(CTParam(69).gt.0.001)then
             momtot=dsqrt((pz(ien)+ffermpz(ien))**2
     $            +(py(ien)+ffermpy(ien))**2
     $            +(px(ien)+ffermpx(ien))**2)
             ptrans=dsqrt((py(ien)+ffermpy(ien))**2
     $            +(px(ien)+ffermpx(ien))**2)
             pseudorap=0.5d0*dlog((momtot+ptrans)/(momtot-ptrans)) 
c             write(*,*)ien, pseudorap
             oldcto=CTOption(46)
             CTOption(46)=1
             bdum=paulibl(ien,quarkdens,pseudorap)
             quarkdens=quarkdens/0.16d0/3.0d0
c     write(*,*)quarkdens, ien
             CToption(46)=oldcto
          endif
          if((quarkdens.le.ctparam(69)).and.(ncoll(ien).gt.0))then
             ecor=ecor+p0(ien)*1000.0d0
          endif

          if((ncoll(ien).gt.0).and.(rapyzpar.le.CTParam(66))
     $         .and.(quarkdens.gt.ctparam(69)))then
             
             etotu=etotu+p0(ien)*1000
             pxtotu=pxtotu+dabs(px(ien))*1000
             pytotu=pytotu+dabs(py(ien))*1000
             pztotu=pztotu+dabs(pz(ien))*1000
             if(abs(ityp(ien)).lt.100)then
                bartotu=bartotu+abs(ityp(ien))/ityp(ien)
             end if
          end if
       end do

       write(*,300) 'core fraction: ', (1.0-(ecor/(ecor+etotu)))*100.0
 300   format(a15,f5.1,'%')

       do index=1,npart
          rapyzpar=0.5d0*dabs(dlog((p0(index)+pz(index))
     $         /(p0(index)-pz(index))))
          
          rapypar=0.5d0*dabs(dlog((p0(index)+dsqrt(pz(index)*pz(index)+
     $         px(index)*px(index)+py(index)*py(index)))
     $         /(p0(index)-dsqrt(pz(index)*pz(index)+
     $         px(index)*px(index)+py(index)*py(index)))))          

           quarkdens=10.0d0
          if(CTParam(69).gt.0.001)then
             momtot=dsqrt((pz(index)+ffermpz(index))**2
     $            +(py(index)+ffermpy(index))**2
     $            +(px(index)+ffermpx(index))**2)
             ptrans=dsqrt((py(index)+ffermpy(index))**2
     $            +(px(index)+ffermpx(index))**2)
             pseudorap=0.5d0*dlog((momtot+ptrans)/(momtot-ptrans)) 
             oldcto=CTOption(46)
             CTOption(46)=1
             bdum=paulibl(index,quarkdens,pseudorap)
             quarkdens=quarkdens/0.16d0/3.0d0
             CToption(46)=oldcto
          endif

          if((ncoll(index).gt.0).and.(rapyzpar.le.CTParam(66))
     $            .and.(quarkdens.gt.ctparam(69)))then

             xpos=rx(index)
             ypos=ry(index)
             zpos=rz(index)


c..   if the particle is on the grid
             if(xpos.lt.dx*ngr/2.0d0.and.xpos.gt.-dx*ngr/2.0d0.and.
     &            ypos.lt.dx*ngr/2.0d0.and.ypos.gt.-dx*ngr/2.0d0.and.
     &            zpos.lt.dx*ngr/2.0d0.and.zpos.gt.-dx*ngr/2.0d0)then
                


c.. cell index defintion
                i=anint(xpos/dx+0.5d0+0.5d0*ngr)
                j=anint(ypos/dx+0.5d0+0.5d0*ngr)
                k=anint(zpos/dx+0.5d0+0.5d0*ngr)
c..   to save computing time the tails of the Gaussians are cut off
                tail=IDINT(4.0d0/dx+0.5d0)   
                if(i-tail.ge.0.and.i+tail.lt.ngr.and.
     &               j-tail.ge.0.and.j+tail.lt.ngr.and.
     &               k-tail.ge.0.and.k+tail.lt.ngr)then
                   irr=i-tail+1
                   jrr=j-tail+1
                   krr=k-tail+1
                   d=2*tail-1
c..   if the particle is on the boundary nothing is cut off
                else
                   irr=1
                   jrr=1
                   krr=1
                   d=ngr-1
                end if      
        
c..   gamma in z-direction
                gamma_gauss=sqrt(p0(index)**2/(p0(index)**2-
     &               (pz(index)+ffermpz(index))**2))
               if (p0(index)-abs(pz(index)+ffermpz(index)).le.1.d-8)then
                   
                   write(*,*) 'urqmdtohydro:',gamma_gauss, p0(index), 
     &                  pz(index)+ffermpz(index)
                endif
                
c..   calculate norm numerically
                enorm = 0d0
        do 702 ir=irr,irr+d
           do 602 jr=jrr,jrr+d
              do 502 kr=krr,krr+d
c..   grid position where something is put
                 xposgr=(ir-ngr/2.0d0)*dx
                 yposgr=(jr-ngr/2.0d0)*dx
                 zposgr=(kr-ngr/2.0d0)*dx
c..   exponent (x_c-x)^2
                 rq= (-xpos+xposgr)**2+(-ypos+yposgr)**2
     &                +(-zpos+zposgr)**2*gamma_gauss**2              
c..   integrating gauss
          
                 
                 enorm=enorm
     &                + exp(-rq/(2.0d0*sigma**2))*dx*dx*dx
                 
                 
 502          continue    
 602       continue  
 702    continue         
        

        
c..   general norm for the Gaussian distributions
        norm = 1d0/enorm
c..   normalized densities
        energy= p0(index)*norm*1000/e0
        phx=(px(index)+ffermpx(index))*norm*1000/e0
        phy=(py(index)+ffermpy(index))*norm*1000/e0
        phz=(pz(index)+ffermpz(index))*norm*1000/e0
c..   baryon number         
        barnum=0.0d0
           if(abs(ityp(index)).lt.100)then
              barnum=ityp(index)/abs(ityp(index))*norm/n0
           end if        
           
c..   put the distributions on the grid        
           do 700 ir=irr,irr+d
              do 600 jr=jrr,jrr+d
                 do 500 kr=krr,krr+d
c..   grid position where something is put
                    xposgr=(ir-ngr/2.0d0)*dx
                    yposgr=(jr-ngr/2.0d0)*dx
                    zposgr=(kr-ngr/2.0d0)*dx
c..   exponent (x_c-x)^2
                   rq= (-xpos+xposgr)**2+(-ypos+yposgr)**2
     &                  +(-zpos+zposgr)**2*gamma_gauss**2              
c..   addition of densities in the cell
                   if(energy*exp(-rq/(2.0d0*sigma**2)).ge.1.0d-12)then
                      elab(ir,jr,kr)=elab(ir,jr,kr)
     &                     + energy*exp(-rq/(2.0d0*sigma**2))            
                      mx(ir,jr,kr)=mx(ir,jr,kr)
     &                     + phx*exp(-rq/(2.0d0*sigma**2))            
                      my(ir,jr,kr)=my(ir,jr,kr)
     &                     + phy*exp(-rq/(2.0d0*sigma**2))            
                      mz(ir,jr,kr)=mz(ir,jr,kr)
     &                     + phz*exp(-rq/(2.0d0*sigma**2))            
                      r(ir,jr,kr)=r(ir,jr,kr)
     &                     + barnum*exp(-rq/(2.0d0*sigma**2))
                      
                      etothy=etothy+energy*exp(-rq/(2.0d0*sigma**2))            
                      pxtothy=pxtothy+abs(phx)*exp(-rq/(2.0d0*sigma**2))            
                      pytothy=pytothy+abs(phy)*exp(-rq/(2.0d0*sigma**2))            
                      pztothy=pztothy+abs(phz)*exp(-rq/(2.0d0*sigma**2))            
                      bartot=bartot+barnum*exp(-rq/(2.0d0*sigma**2))           
                      
                   end if
 500            continue    
 600         continue  
 700      continue      
       end if           
      end if
      end do   
      
      else 

c.. put also spectators on hydro grid
c.. total energy and momenta      
       do ien=1,npart
        etotu=etotu+p0(ien)*1000
        pxtotu=pxtotu+abs(px(ien))*1000
        pytotu=pytotu+abs(py(ien))*1000
        pztotu=pztotu+abs(pz(ien))*1000
         if(abs(ityp(ien)).lt.100)then
          bartotu=bartotu+abs(ityp(ien))/ityp(ien)
         end if
       end do 
     
       do index=1,npart
        xpos=rx(index)
        ypos=ry(index)
        zpos=rz(index)
c.. if the particle is on the grid
        if(xpos.lt.dx*ngr/2.0d0.and.xpos.gt.-dx*ngr/2.0d0.and.
     &     ypos.lt.dx*ngr/2.0d0.and.ypos.gt.-dx*ngr/2.0d0.and.
     &     zpos.lt.dx*ngr/2.0d0.and.zpos.gt.-dx*ngr/2.0d0)then



c.. cell index defintion
         i=anint(xpos/dx+0.5d0+0.5d0*ngr)
         j=anint(ypos/dx+0.5d0+0.5d0*ngr)
         k=anint(zpos/dx+0.5d0+0.5d0*ngr)
c.. to save computing time the tails of the Gaussians are cut off
         tail=ngr/10   
         if(i-tail.ge.0.and.i+tail.lt.ngr.and.
     &      j-tail.ge.0.and.j+tail.lt.ngr.and.
     &      k-tail.ge.0.and.k+tail.lt.ngr)then
          irr=i-tail
          jrr=j-tail
          krr=k-tail
          d=2*tail
c.. if the particle is on the boundary nothing is cut off
         else
          irr=0
          jrr=0
          krr=0
          d=ngr       
         end if      

c.. gamma in z-direction
          gamma_gauss=sqrt(p0(index)**2/(p0(index)**2-
     &              (pz(index)+ffermpz(index))**2))
         if (p0(index)-abs(pz(index)+ffermpz(index)).le.1.d-8)then

          write(*,*) 'urqmdtohydro:',gamma_gauss, p0(index), 
     &    pz(index)+ffermpz(index)
         endif
c.. general norm for the Gaussian distributions
          norm = sqrt(1/(2.0d0*pi*sigma**2))**3*gamma_gauss
c.. normalized densities
          energy= p0(index)*norm*1000/e0
          phx=(px(index)+ffermpx(index))*norm*1000/e0
          phy=(py(index)+ffermpy(index))*norm*1000/e0
          phz=(pz(index)+ffermpz(index))*norm*1000/e0
c.. baryon number         
          barnum=0.0d0
          if(abs(ityp(index)).lt.100)then
           barnum=ityp(index)/abs(ityp(index))*norm/n0
          end if         
        
c.. put the distributions on the grid        
          do 701 ir=irr,irr+d
           do 601 jr=jrr,jrr+d
            do 501 kr=krr,krr+d
c.. grid position where something is put
             xposgr=(ir-ngr/2.0d0)*dx
             yposgr=(jr-ngr/2.0d0)*dx
             zposgr=(kr-ngr/2.0d0)*dx
c.. exponent (x_c-x)^2
             rq= (-xpos+xposgr)**2+(-ypos+yposgr)**2
     &        +(-zpos+zposgr)**2*gamma_gauss**2              
c.. addition of densities in the cell
             if(energy*exp(-rq/(2.0d0*sigma**2)).ge.1.0d-12)then
              elab(ir,jr,kr)=elab(ir,jr,kr)
     &            + energy*exp(-rq/(2.0d0*sigma**2))            
              mx(ir,jr,kr)=mx(ir,jr,kr)
     &            + phx*exp(-rq/(2.0d0*sigma**2))            
              my(ir,jr,kr)=my(ir,jr,kr)
     &            + phy*exp(-rq/(2.0d0*sigma**2))            
              mz(ir,jr,kr)=mz(ir,jr,kr)
     &            + phz*exp(-rq/(2.0d0*sigma**2))            
              r(ir,jr,kr)=r(ir,jr,kr)
     &            + barnum*exp(-rq/(2.0d0*sigma**2))

              etothy=etothy+energy*exp(-rq/(2.0d0*sigma**2))            
              pxtothy=pxtothy+abs(phx)*exp(-rq/(2.0d0*sigma**2))            
              pytothy=pytothy+abs(phy)*exp(-rq/(2.0d0*sigma**2))            
              pztothy=pztothy+abs(phz)*exp(-rq/(2.0d0*sigma**2))            
              bartot=bartot+barnum*exp(-rq/(2.0d0*sigma**2))           
          
             end if
 501       continue    
 601      continue  
 701     continue      
        end if           
       end do   
      end if   
     

c.. calculation of conserved quantities and check output versus input

      etothy=etothy*e0*dx**3
      pxtothy=pxtothy*e0*dx**3       
      pytothy=pytothy*e0*dx**3       
      pztothy=pztothy*e0*dx**3       
      bartot=bartot*n0*dx**3

      h_etot=etotu/1.d3
      h_btot=bartotu
      
      write(*,301) 'particlization','urqmd','fluid','diff'
      write(*,302) 'energy:           '
     +              ,h_etot,etothy/1.d3,(etothy/1.d3-h_etot)/h_etot*100
      write(*,302) 'baryons:          '
     +              ,h_btot,bartot, (bartot-h_btot)/h_btot*100
      write(*,303) 'strangeness:', strangebe
      write(*,303) 'charge:     ', -zspec+zt+zp

      if (abs((etothy/1.d3-h_etot)/h_etot).gt.0.1d0)then
         write(*,*) 'energy difference is greater than 10 %'
         write(*,*) 'evolution is stopped'
c         stop 137
      end if   
 301  format(a14,x,a5,x,a7,x,a6,x,'(conservation check)')
 302  format(a12,x,f7.1,x,f7.1,x,f5.1,'%')
 303  format(a12,9x,i7)

      h_btot=anint(bartot)
c      write(*,*) 'energy comparison (etotu vs. etothy)',
c     &      etotu/1.d3,etothy/1.d3
c      write(*,*) 'x-momentum comparison (pxtotu vs. pxtothy)',
c     &      pxtotu/1.d3,pxtothy/1.d3
c      write(*,*) 'y-momentum comparison (pytotu vs. pytothy)',
c     &      pytotu/1.d3,pytothy/1.d3
c      write(*,*) 'z-momentum comparison (pztotu vs. pztothy)',
c     &      pztotu/1.d3,pztothy/1.d3
c      write(*,*) 'bayon number comparison (bartotu vs. bartothy):',
c     &      bartotu,bartot      

c      write(*,*) 'Energy conserved by: ',
c     &  ((etothy-etotu)/etotu)*100,'%'
c      write(*,*) 'X-momentum conserved by: ',
c     &  (pxtothy-pxtotu)/(pxtotu)*100,'%'
c      write(*,*) 'Y-momentum conserved by: ',
c     &  (pytothy-pytotu)/(pytotu)*100,'%'
c      write(*,*) 'Z-momentum conserved by: ',
c     &  (pztothy-pztotu)/(pztotu)*100,'%'
c      write(*,*) 'Baryon number conserved by ',
c     &   (bartot-bartotu)/(bartotu)*100,'%'
            
        
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coopertime(ifr,jfr,kfr,tf,dstfr,dsxfr,dsyfr,dszfr,
     &           tfr,muqfr,musfr,vxfr,vyfr,vzfr,mmax)
      implicit none

      
      include 'comres.f'
      include 'coms.f'
      INCLUDE 'defs.f'

      real*4 ifr(ngr*ngr*ngr),jfr(ngr*ngr*ngr),kfr(ngr*ngr*ngr)
      integer mmax,m
      real*4 tf(ngr*ngr*ngr),tfr(ngr*ngr*ngr),muqfr(ngr*ngr*ngr)
      real*4 musfr(ngr*ngr*ngr)
      real*8 vyfr(ngr*ngr*ngr),vxfr(ngr*ngr*ngr)
      real*8 vzfr(ngr*ngr*ngr),gfr(ngr*ngr*ngr)
      real*8 dstfr(ngr*ngr*ngr),dsxfr(ngr*ngr*ngr)
      real*8 dsyfr(ngr*ngr*ngr),dszfr(ngr*ngr*ngr)       
c.. hydro parameters
      real*8 dx,dt_hydro,t,vol,cut
      real*8 e0,n0,Bag,hc3,pi2
c.. arrays to transfer hydro particles to urqmd
      real*8 h_r0(nmax), h_rx(nmax), h_ry(nmax), h_rz(nmax),
     +     h_p0(nmax), h_px(nmax), h_py(nmax), h_pz(nmax),
     +     h_fmass(nmax), h_tform(nmax),h_xtotfac(nmax)
      integer h_ncoll(nmax),h_charge(nmax),h_ityp(nmax),
     &     h_lstcoll(nmax),h_iso3(nmax),h_origin(nmax)
      real*8 h_etot,h_btot,starttime
      integer n
      real*8 Peaks 
c.. hydro common blocks
      common /gitter/ dx,dt_hydro,t,vol,cut
c      common /onefgrid/ elab,mx,my,mz,r,mmax
      common /grstate/ e0,n0,Bag,hc3,pi2
c.. thermodynamic quantities
      common /start/ starttime 
      common /line/ n  
c.. transfer array
      common /hcoor/h_r0,h_rx,h_ry,h_rz,h_p0, h_px, h_py, h_pz,
     +     h_fmass,h_tform,h_xtotfac,
     &     h_ncoll,h_charge,h_ityp,h_lstcoll,h_iso3,h_origin
c.. energy conservation
      common /h_etot/ h_etot,h_btot

c.. this program operates in GeV      
      
c.. volume in GeV^(-3) 
      real*8 volume,beta2
c.. masses of particles and id      
      real*8 massit,mass
      integer iso,ii,it(0:20),barit,isospin,ch    
c.. degeneracy factors
      integer isoit,fchg,iso3it,strit  
c.. variables for momentum simulation      
      real*8 Temp,mubar,mustr,gam,ranf,mu
      real*8 vxb,vyb,vzb,lp0,lpx,lpy,lpz
      real*8 dst,dsx,dsy,dsz
      integer q,iitmax,netbaryons,tot_strange
      integer netbaryons2,mode
c..   charge conservation
      real*8 diffold,des_ch,diffnew,tot_charge,i,j,k
      
      real*8 econ,econa,econ2,econa2,econs,econ3,econa3,econ4
      

      
c..   Z/A of the spectators        
      integer numspec,Zspec,Aspec,strangebe
      common /numspec/ numspec,Zspec,Aspec,strangebe
      
      econ=0.0d0
      econ2=0.0d0
      econs=0.0d0
      econ3=0.0d0
      econ4=0.0d0
 
      netbaryons=0
      netbaryons2=0
      tot_charge=0

      tot_strange=0
      


      
      
c..   mode=1 produce all strange particles
c..   mode=2 produce the equivalent number of antistrange particles
c..   mode=3 produce all nonstrange antibaryons
c..   mode=4 produce all nonstrange baryons
c..   mode=5 produce negatively charged non-strange mesons
c..   mode=6 produce positively charged non-strange mesons 
c..   mode=7 Energy conservation by non-strange neutral mesons

      write(*,'(a1,x,a4,x,a6,x,a5,x,a6,x,a6)')
     %      'm','bary','charge','stran','energy','ch_des'
      do mode=1,7
         
         
c..   randomly generate particles until baryon number is reached        
 10      continue
         
c..   choose randomly a cell
         m=int(ranf(0)*(mmax-1)+1)
         
         
         if(tfr(m).gt.3.d0) then
c     .. transformation to GeV and usual definitions of the chemical potentials
            Temp=tfr(m)*1.0d-3
            mubar=3.0d0*muqfr(m)*1.0d-3         
            mustr=(muqfr(m)-musfr(m))*1.0d-3
            
            vxb=vxfr(m)
            vyb=vyfr(m)
            vzb=vzfr(m)             
            beta2=vxb**2+vyb**2+vzb**2 
            gam=1.d0/(sqrt(1.d0-beta2)) 
            i=ifr(m)
            j=jfr(m)
            k=kfr(m)
            dst=dstfr(m)
            dsx=dsxfr(m)
            dsy=dsyfr(m)
            dsz=dszfr(m)  

c           avoid acausal transport

	    if((tf(m)+starttime).lt.dabs(k)) goto 10
c..  volume of element via scalar product            
            volume=((1/hqc)**3)*gam*(dst+dsx*vxb+dsy*vyb+dsz*vzb)    
             

       if(volume.lt.0.d0)then
c             write(*,*) 'negative volume element, omitted so far'
             go to 10   
       end if
  
c..   subroutine mkityp decides if and which particle is produced
            call mkityp(Temp,mubar,mustr,gam,it,iitmax,volume)

c..   itmax equals zero, then go to next cell         
            if(iitmax.eq.0) go to 10
c..   loop over all produced particles 
            do q=1,iitmax             

               
               
c..   count baryons
               


c..   produce only the wanted particles in this mode             
               if (mode.eq.2.and.strit(it(q)).ge.0) goto 12
               if (mode.eq.4.and.(it(q).lt.0.or.it(q).ge.100.or.
     $              strit(it(q)).ne.0))goto 12
               
               des_ch=netbaryons*(Zt+Zp-Zspec+0d0)/(At+Ap-Aspec+0d0)
               iso=isoit(it(q)) 
               isospin=iso3it(iso)
               ch=fchg(isospin,it(q)) 
         if (mode.eq.6.and.(abs(it(q)).lt.100.or.strit(it(q))
     $              .ne.0.or.ch.le.0))goto 12

         if (mode.eq.7.and.(abs(it(q)).lt.100.or.strit(it(q))
     $              .ne.0.or.ch.ne.0))goto 12               

               

               

               call cooper(it(q),Temp,mubar,mustr,gam,vxb,vyb,vzb,
     &              dst,dsx,dsy,dsz,lp0,lpx,lpy,lpz) 
                     
               if(lp0.lt.1.d-8) then
                  if(mode.eq.1.and.abs(it(q)).lt.100)then
                     netbaryons2=netbaryons2-it(q)/abs(it(q))
                  end if   
                  go to 10
               end if

               if(mode.eq.1)then
                  econ2=econ2+lp0
               end if
               if(mode.eq.3)then
                  econ3=econ3+lp0
               end if
               if(mode.eq.5)then
                  econ4=econ4+lp0  
               end if  
      
               
c..   produce only the wanted particles in this mode             
               if (mode.eq.1.and.strit(it(q)).le.0) goto 12 
               if (mode.eq.3.and.(it(q).gt.0.or.it(q).le.-100.or.
     $              strit(it(q)).ne.0))goto 12

               if (mode.eq.5.and.(abs(it(q)).lt.100.or.
     &            strit(it(q)).ne.0.or.ch.ge.0))goto 12
               
c..   store all necessary particle information in array         
               n=n+1
c     
c     Calculate position vectors (single precision).
c    
   
               
               
               h_r0(n)=tf(m)+starttime
               h_rx(n)=i
               h_ry(n)=j
               h_rz(n)=k
               

               h_p0(n)=lp0
               h_px(n)=lpx
               h_py(n)=lpy
               h_pz(n)=lpz
             
               h_fmass(n)=massit(it(q))
               h_ityp(n)=it(q)
               h_iso3(n)=isospin
               h_charge(n)=ch   

               h_lstcoll(n)=0
               h_ncoll(n)=0
               h_origin(n)=96
               h_tform(n)=0.0d0
               h_xtotfac(n)=1.0d0  
c..   count the total charge and strange quarks   
               tot_charge=tot_charge+h_charge(n)
               tot_strange=tot_strange+strit(h_ityp(n))          

c..   count baryons
               if(abs(h_ityp(n)).lt.100)then
                  netbaryons=netbaryons+h_ityp(n)/abs(h_ityp(n))
               end if   
               econa=econ
               econ=econ+h_p0(n)

               
 12            continue
            end do


         end if           
c..   produce particles until baryon number is reached
c..   initial baryon number on the hydro grid 
        if (mode.eq.1.and.econ2.lt.h_etot) goto 10 
        if (mode.eq.2.and.tot_strange.gt.strangebe) goto 10
        if (mode.eq.3.and.econ3.le.h_etot) goto 10
        if (mode.eq.4.and.netbaryons.le.h_btot) goto 10
        if (mode.eq.5.and.econ4.le.h_etot) goto 10
        if (mode.eq.6.and.tot_charge.lt.des_ch) goto 10    

        if(mode.eq.7.and.econ.lt.h_etot*0.999) goto 10 
        if(mode.eq.7.and.(dabs(econa-h_etot)
     &   .gt.dabs(econ-h_etot))) goto 10
      
      netbaryons=0
      econs=0.0d0
      do ii=1,n-1
         econs=econs+h_p0(ii)
         if(abs(h_ityp(ii)).lt.100)then
            netbaryons=netbaryons+h_ityp(ii)/abs(h_ityp(ii))
         end if         
      enddo
      
      if(mode.eq.4)then
         netbaryons=netbaryons+1 
      endif

      write(*,'(i1,x,i4,x,f6.1,x,i4,x,f7.1,x,f6.1)')
     &     mode,netbaryons,tot_charge,tot_strange,econs,des_ch
      
      end do                    ! end of mode loop
      
      
       
      npart=n-1
      
      return 
      
      end
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine arraysort(thydro_start)
      implicit none

      include 'coms.f'
      include 'comres.f'

      integer i,j,n,b,m
      real*8 thydro_start



      real*8 br0(nmax),brx(nmax),bry(nmax),brz(nmax)
      real*8 bp0(nmax),bpx(nmax),bpy(nmax),bpz(nmax)
      real*8 bfmass(nmax),btform(nmax),bxtotfac(nmax),bdectime(nmax)
      integer bityp(nmax),biso3(nmax),bcharge(nmax),blstcoll(nmax)
      integer borigin(nmax),bncoll(nmax)

      real*8 mr0(nmax),mrx(nmax),mry(nmax),mrz(nmax)
      real*8 mp0(nmax),mpx(nmax),mpy(nmax),mpz(nmax)
      real*8 mfmass(nmax),mtform(nmax),mxtotfac(nmax),mdectime(nmax)
      integer mityp(nmax),miso3(nmax),mcharge(nmax),mlstcoll(nmax)
      integer morigin(nmax),mncoll(nmax)


c.. arrays to transfer hydro particles to urqmd
      real*8 h_r0(nmax), h_rx(nmax), h_ry(nmax), h_rz(nmax),
     +     h_p0(nmax), h_px(nmax), h_py(nmax), h_pz(nmax),
     +     h_fmass(nmax), h_tform(nmax),h_xtotfac(nmax)
      integer h_ncoll(nmax),h_charge(nmax),h_ityp(nmax),
     &     h_lstcoll(nmax),h_iso3(nmax),h_origin(nmax)

c.. transfer array
      common /hcoor/h_r0,h_rx,h_ry,h_rz,h_p0, h_px, h_py, h_pz,
     +     h_fmass,h_tform,h_xtotfac,
     &     h_ncoll,h_charge,h_ityp,h_lstcoll,h_iso3,h_origin
 
      common /line/ n 
      
      do i=1,n-1
        h_tform(i)=h_r0(i) 
        h_rx(i)=h_rx(i)-h_px(i)/h_p0(i)*(h_tform(i)-thydro_start)
        h_ry(i)=h_ry(i)-h_py(i)/h_p0(i)*(h_tform(i)-thydro_start)
        h_rz(i)=h_rz(i)-h_pz(i)/h_p0(i)*(h_tform(i)-thydro_start)
        h_r0(i)=thydro_start
        if(h_tform(i).gt.0.0d0)then
           h_xtotfac(i)=0.0d0
        endif 
      end do   

      b=0
      m=0

      do i=1,n-1 
       if(abs(h_ityp(i)).le.maxbar)then
        b=b+1
        br0(b)=h_r0(i)
        brx(b)=h_rx(i)
        bry(b)=h_ry(i)
        brz(b)=h_rz(i)
        bp0(b)=h_p0(i)
        bpx(b)=h_px(i)
        bpy(b)=h_py(i)
        bpz(b)=h_pz(i)
       
        bityp(b)=h_ityp(i)
        biso3(b)=h_iso3(i)
        bfmass(b)=h_fmass(i)        
        bcharge(b)=h_charge(i)
        blstcoll(b)=h_lstcoll(i)
        bncoll(b)=h_ncoll(i)
        borigin(b)=h_origin(i)
        bxtotfac(b)=h_xtotfac(i)  
        btform(b)=h_tform(i)
       else
        m=m+1
        mr0(m)=h_r0(i)
        mrx(m)=h_rx(i)
        mry(m)=h_ry(i)
        mrz(m)=h_rz(i)
        mp0(m)=h_p0(i)
        mpx(m)=h_px(i)
        mpy(m)=h_py(i)
        mpz(m)=h_pz(i)
       
        mityp(m)=h_ityp(i)
        miso3(m)=h_iso3(i)
        mfmass(m)=h_fmass(i)        
        mcharge(m)=h_charge(i)
        mlstcoll(m)=h_lstcoll(i)
        mncoll(m)=h_ncoll(i)
        morigin(m)=h_origin(i)
        mxtotfac(m)=h_xtotfac(i)  
        mtform(m)=h_tform(i)
       endif 
      end do
      do i=1,b
       h_r0(i)=br0(i) 
       h_rx(i)=brx(i)
       h_ry(i)=bry(i)
        h_rz(i)=brz(i)
        h_p0(i)=bp0(i)
        h_px(i)=bpx(i)
        h_py(i)=bpy(i)
        h_pz(i)=bpz(i)
        h_ityp(i)=bityp(i)
        h_iso3(i)=biso3(i)
        h_fmass(i)=bfmass(i)
        h_charge(i)=bcharge(i)
        h_lstcoll(i)=blstcoll(i)
        h_ncoll(i)=bncoll(i)
        h_origin(i)=borigin(i)
        h_xtotfac(i)=bxtotfac(i)
        h_tform(i)=btform(i)
      end do
      do i=b+1,b+m
        h_r0(i)=mr0(i-b) 
        h_rx(i)=mrx(i-b)
        h_ry(i)=mry(i-b)
        h_rz(i)=mrz(i-b)
        h_p0(i)=mp0(i-b)
        h_px(i)=mpx(i-b)
        h_py(i)=mpy(i-b)
        h_pz(i)=mpz(i-b)
        h_ityp(i)=mityp(i-b)
        h_iso3(i)=miso3(i-b)
        h_fmass(i)=mfmass(i-b)
        h_charge(i)=mcharge(i-b)
        h_lstcoll(i)=mlstcoll(i-b)
        h_ncoll(i)=mncoll(i-b)
        h_origin(i)=morigin(i-b)
        h_xtotfac(i)=mxtotfac(i-b)
        h_tform(i)=mtform(i-b)
      end do

     
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      subroutine mkityp(Temp,mubar,mustr,gam,it,iitmax,volume)
c.. this subroutine decides if and which particle is produced

      implicit none
      
      include 'comres.f'
      include 'itypdata.f'
      include 'coms.f'
     
      real*8 Temp,mubar,mustr,gam
      integer it(0:20),ii,id,degfac,iitmax
      integer strit,barit,jit,isoit
      real*8 massit,bessk,volume,mass
      integer b,idum,factorial,kmax,q,i
      real*8 mu,bosesum,pnd(265),nop,partdens
      real*8 prob,probid,ranf,dx,poisson,y

c.. loop over all urqmd particles to calculate the particle 
c.. number densities per cell
      bosesum=0.0d0         
      nop=0.0d0
      do 100 ii=1,265
c.. if it is a valid urqmd ityp
         if(itypes(ii).ne.0)then
            id=itypes(ii)
c..   call the necessary quantum numbers of this particle type           
            mass=massit(id)
c..   degeneracy factor           
        degfac=(jit(id)+1)*(isoit(id)+1)
c..   complete chemical potential for all species           
        mu=-strit(id)*mustr+barit(id)*mubar
        
c..   Bose distribution for pions         
      if((id.eq.101).or.(iabs(id).eq.106).or.(iabs(id).eq.108))then
	   bosesum=0.0d0
           do b=1,10
              bosesum=bosesum+(1.0d0/b)*bessk(2,(b*(mass/Temp)))  
     $        *dexp(b*mu/temp)       
           end do        
           pnd(ii)=volume*degfac
     &       *mass**2*Temp/(2.0d0*pi**2)*bosesum         
           
c..   J�ttner distribution for all other particle species         
        else 
           pnd(ii)=volume*degfac
     &          *mass**2*Temp/(2.0d0*pi**2)
     &          *exp(mu/Temp)
     &          *bessk(2,(mass/Temp))
        end if
         nop=nop+pnd(ii)
      end if
 100  continue          
      
c      write(*,*) 'in mkityp.. all densities are calculated'
c.. check if Poisson approximation is good enough (error <=1%)
      if(nop.lt.0.01)then
c.. generate random number to decide if a particle is produced
       prob=ranf(0)         
       if(prob.lt.nop)then
c.. generate random number to decide the ityp         
        probid=ranf(0)*nop
        partdens=0.0d0
c.. loop over all valid urqmd ityps
        do 200 ii=1,265
           if(itypes(ii).ne.0)then
              partdens=partdens+pnd(ii)
              if(probid.lt.partdens)then
                 it(1)=itypes(ii)
                 iitmax=1 
                 go to 201
              endif
           endif
 200    continue            
 201    continue 
      else
         it(1)=0
         iitmax=0
      end if
      else
c.. full Poisson distribution is used

         kmax=5
 937     iitmax=int(ranf(0)*kmax)
         y=ranf(0)
         poisson=nop**iitmax/factorial(iitmax)*exp(-nop)
         if (y.gt.poisson) goto 937
       if (iitmax.gt.0) then
        do i=1,iitmax
c.. generate random number to decide the ityp         
           probid=ranf(0)*nop
           partdens=0.0d0
c.. loop over all valid urqmd ityps
         do 700 ii=1,265
          if(itypes(ii).ne.0)then
c           write(*,*) 'in poisson', ii, pnd(ii)
             partdens=partdens+pnd(ii)
           if(probid.lt.partdens)then
              it(i)=itypes(ii)
              go to 701
           endif
          endif
 700     continue            
 701     continue 
        end do   
       else
          it(1)=0
          iitmax=0
       end if
      end if
      return
      end 
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function factorial(i)
      
      implicit none
      integer i,k
      
      factorial=1
      if (i.eq.0) return
      do k=1,i
         factorial=factorial*k
      end do
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                           
      integer function barit(ia)
      implicit none
c.. this function gives the baryon number of each hadron species      
      integer ia
                                                                                                                                                                    
      if(abs(ia).gt.100)then
       barit=0
      else if(ia.lt.0)then
       barit=-1
      else
       barit=1
      end if              
       
      return
      end                                                                                                                                                                     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      integer function iso3it(is)
      implicit none 
c.. this function defines the iso3 component via Monte Carlo      
     
     
      integer is
      real*8 probiso,ranf
      
c.. isospin division            
      if (is.eq.0)then
       iso3it=0
      else if (is.eq.1) then
       probiso=ranf(0)*2.0d0
       if(probiso.lt.1.0d0) then
        iso3it=1
       else 
        iso3it=-1
       endif
      else if (is.eq.2) then
       probiso=ranf(0)*3.0d0           
       if(probiso.lt.1.0d0) then
        iso3it=2 
       else if(probiso.lt.2.0d0) then
        iso3it=0
       else
        iso3it=-2
       endif
      else if (is.eq.3)then  
       probiso=ranf(0)*4.0d0
       if(probiso.lt.1.0d0) then
        iso3it=-3    
       else if(probiso.lt.2.0d0) then
        iso3it=-1
       else if(probiso.lt.3.0d0) then
        iso3it=1
       else 
        iso3it=3
       end if
      else
        write(6,*) 'unknown isospin', is
      end if 

      return 
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine cooper(id,T,mub,mus,g,vx,vy,vz,
     &      dst,dsx,dsy,dsz,p0,px,py,pz) 
c.. this subroutine generates the momenta
    
       implicit none

       include 'options.f' 
        
       real*8 T,mub,mus,g,vx,vy,vz,p0,px,py,pz,mu
       integer id,a,i
       real*8 massit,m,fmax,f,fp,ptmax,pzmax
       integer strit,barit
       real*8 ranf,den,energy
       real*8 dst,dsx,dsy,dsz           
       real*8  time,  acttime, bdist, ebeam, bimp,bmin,ecm
       common /rsys/ time,acttime,bdist,bimp,bmin,ebeam,ecm
       real*8 pxmax,pxmin,pymax,pymin,pzmin,lhs,rsph
       
   
       real*8 intpxmin,intpymin,intpzmin,fpxmax,fpymax,fpzmax
       real*8 dist,dp
       integer ix,iy,iz
       
       

       m=massit(id)
       mu=-strit(id)*mus+barit(id)*mub

c.. maximum pt and pz range for particle generation       
       ptmax=ecm
       pzmax=2.d0*ecm
c.. bose or fermi distribution       
       a=1
       if (abs(id).gt.100)a=-1

c.. maximum of the distribution is calculated via analytical formula 
c.. for the zeroth of the derivative
c        den=vx**2+vy**2+vz**2-1.0d0

c        px=abs(vx)/vx*sqrt(-vx**2*m**2/den)
c        py=abs(vy)/vy*sqrt(-vy**2*m**2/den)
c        pz=abs(vz)/vz*sqrt(-vz**2*m**2/den)
c        p0=sqrt(m**2+px**2+py**2+pz**2)
c        fmax=fp(p0,px,py,pz,vx,vy,vz,g,a,T,mu)

cc       write(*,*) 'momenta at max:',px,py,pz,p0
c        i=0

c        if(CTOption(53).eq.0) then

cc       define "radius" for cut at exp(-30)*fmax
cc        rsph=sqrt((30.d0*T+m)**2-m**2)
c        lhs=exp(30.d0)*(exp(m/T)+a*exp(mu/T))-a*exp(mu/T)
c        if(lhs.lt.0.d0) then
c          write(*,*)'ERROR:negative argument for log!'
c          write(*,*) 'm, mu,T,a,lhs', m, mu,T,a,lhs
c          stop 137
c        endif
c        if((T*log(lhs))**2-m**2.lt.0) then
c          write(*,*)'ERROR:sqrt of negative nuber!'
c          stop 137
c        endif

cc..     radius
c        rsph=sqrt((T*log(lhs))**2-m**2)

cc        write(*,*)'T,mu,radius',T,mu,rsph

c        call lorsphere(rsph,g,vx,vy,vz,m,pxmax,pxmin,
c     & pymax,pymin,pzmax,pzmin)

cccccccccccccccccc       check     ccccccccccccccccccccccccccc
c        if(px.gt.pxmax.or.px.lt.pxmin) then
c            write(*,*)'smt wrong with max/min values of px'
c            write(*,*)'px at max:',px,'pxmax,pxmin:',pxmax,pxmin
c            stop 137
c        endif
c        if(py.gt.pymax.or.py.lt.pymin) then
c            write(*,*)'smt wrong with max/min values of py'
c            write(*,*)'py at max:',py,'pymax,pymin:',pymax,pymin
c            stop 137
c        endif
c        if(pz.gt.pzmax.or.pz.lt.pzmin) then
c            write(*,*)'smt wrong with max/min values of pz'
c            write(*,*)'pz at max:',pz,'pzmax,pzmin:',pzmax,pzmin
c            stop 137
c        endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      

        call Peaks(fmax,fpxmax,fpymax,fpzmax,m,a,T,mu,  
     &             vx,vy,vz,dst,dsx,dsy,dsz)

c.. minimum/maximum momentum --- check from time to time
        pxmin=fpxmax-5.d0
        pxmax=fpxmax+5.d0
        pymin=fpymax-5.d0
        pymax=fpymax+5.d0
        pzmin=fpzmax-10.d0
        pzmax=fpzmax+10.d0

        i=0
   
   
 1     continue
       i=i+1

       if(CTOption(53).eq.0) then 
          px=pxmin+ranf(0)*(pxmax-pxmin)
          py=pymin+ranf(0)*(pymax-pymin)
          pz=pzmin+ranf(0)*(pzmax-pzmin)
       else
          px=-ptmax+2.*ranf(0)*ptmax
          py=-ptmax+2.*ranf(0)*ptmax
          pz=-pzmax+2.*ranf(0)*pzmax
       endif
   
       p0=sqrt(m**2+px**2+py**2+pz**2)
c       write(*,*) 'generated momenta:',px,py,pz,p0,i 
       f=ranf(0)*fmax
       
       if(fp(p0,px,py,pz,vx,vy,vz,g,dst,dsx,dsy,dsz,a,T,mu)
     $   .gt.fmax)then
        write(*,*)'fmax too small'
        write(*,*) 'fmax',fmax
       write(*,*) 'fp',fp(p0,px,py,pz,vx,vy,vz,g,dst,dsx,dsy,dsz,a,T,mu)
       endif 
       
       if(f.gt.fp(p0,px,py,pz,vx,vy,vz,g,dst,dsx,dsy,dsz,a,T,mu)
     $  .and.i.le.10000000) then
c       if(f.gt.fp(p0,px,py,pz,vx,vy,vz,g,a,T,mu)) then
        go to 1
       end if
       if(i.gt.10000000)then
        p0=0.0d0
       end if 
c       write(*,*) 'number of trials in cooper:', i

       return 
       end
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       real*8 function fp(p0,px,py,pz,vx,vy,vz,g,
     $                      dst,dsx,dsy,dsz,a,T,mu)
 
       implicit none
       real*8 p0,px,py,pz,vx,vy,vz,g,T,mu
       integer a
       real*8 pmum,dsmupmu
       real*8 dst,dsx,dsy,dsz
c.. the hypersurface vector is covariant disg_mu       
       dsmupmu=dst*p0+dsx*px+dsy*py+dsz*pz
       pmum=g*(p0-px*vx-py*vy-pz*vz)
       fp=1.d0/(exp((pmum-mu)/T)+a)
       fp=dsmupmu*fp/p0      
 
       return
       end
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine inithydro
c
c     Revision : 1.0
c     This subroutine calls initialization procedures for different
c     equations of state and calculation modi.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'coms.f'
      include 'options.f'
      include 'comres.f'
      include 'inputs.f'
      include 'freezeout.f'
      include 'newpart.f'
      include 'boxinc.f'
      include 'colltab.f'

      integer j,k,icount,npold,getspin,fchg,indsp(2),isrt,ipbm
      real*8 nucrad,dstt,dstp,pcm,eb,embeam,emtarget
      real*8 massit,ranf,pcms,dectim
      real*8 gaeq,beeq,galab,belab,ppeq,pteq
      real*8 pboost
      real*8 ratio
      integer AAp, AAt
        integer nnuc
        parameter (nnuc=11)
        save isrt,ipbm
        logical bcorr
        common /ini/ bcorr

        integer i,flagy
        real*8 alf,regula
c momenta
       real*8 mx,my,mz  
c important: never change!
      npold = 0
      nbar=0
      nmes=0
      apt=0
      uid_cnt=0
c reset counters
c
      eb=0D0
c icount is the number of EXTRAordinary pro/tar combinations (i.e. pion ...)
      icount = 0
c reset particle vectors
      do 20 j=1,nmax
        spin(j)  = 0
        ncoll(j) = 0
        lstcoll(j)=0
        r0(j) = 0.0
        rx(j)    = 0.0
        ry(j)    = 0.0
        rz(j)    = 0.0
        p0(j)    = 0.0
        px(j)    = 0.0
        py(j)    = 0.0
        pz(j)    = 0.0
        frr0(j) = 0.0
        frrx(j)    = 0.0
        frry(j)    = 0.0
        frrz(j)    = 0.0
        frp0(j)    = 0.0
        frpx(j)    = 0.0
        frpy(j)    = 0.0
        frpz(j)    = 0.0
        fmass(j) = 0.0
        charge(j)= 0
        iso3(j)  = 0
        ityp(j)  = 0
        dectime(j)= 0.0
        origin(j)=0
        tform(j)=0.0
        xtotfac(j)=1.0
        uid(j)=0
         ffermpx(j) = 0.0
         ffermpy(j) = 0.0
         ffermpz(j) = 0.0
ctd
         do 21 k=1,2
            p0td(k,j)=0.d0
            pxtd(k,j)=0.d0
            pytd(k,j)=0.d0
            pztd(k,j)=0.d0
            fmasstd(k,j)=0.d0
            ityptd(k,j)=0
            iso3td(k,j)=0
 21      continue
 20   continue

c.. initialize urqmd with hydro output
        call gethydrooutput

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine gethydrooutput

      implicit none

      include 'comres.f'
      include 'coms.f'
      include 'options.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'freezeout.f'
      include 'boxinc.f'

      real*8 h_r0(nmax), h_rx(nmax), h_ry(nmax), h_rz(nmax),
     +     h_p0(nmax), h_px(nmax), h_py(nmax), h_pz(nmax),
     +     h_fmass(nmax), h_tform(nmax),h_xtotfac(nmax)
      integer h_ncoll(nmax),h_charge(nmax),h_ityp(nmax),
     &     h_lstcoll(nmax),h_iso3(nmax),h_origin(nmax)
      common /hcoor/h_r0,h_rx,h_ry,h_rz,h_p0, h_px, h_py, h_pz,
     +     h_fmass,h_tform,h_xtotfac,
     &     h_ncoll,h_charge,h_ityp,h_lstcoll,h_iso3,h_origin

      
      common /h_etot/ h_etot,h_btot

      integer i,j,iinelcoll,ttime,itotcoll,bartotu
      real*8 dectim,etotu,pxtotu,pytotu,pztotu,h_etot,h_btot

c.. spectator array
      real*8 specr0(nmax), specrx(nmax), specry(nmax), specrz(nmax),
     +     specp0(nmax), specpx(nmax), specpy(nmax), specpz(nmax),
     +     specfmass(nmax), spectform(nmax),specxtotfac(nmax),
     $     specdectime(nmax)
      integer specncoll(nmax),speccharge(nmax),specityp(nmax),
     &     speclstcoll(nmax),speciso3(nmax),specorigin(nmax)
      common /speccoor/specr0,specrx,specry,specrz,specp0, specpx, 
     &     specpy,specpz,specfmass,spectform,specxtotfac,
     &     specncoll,speccharge,specityp,speclstcoll,speciso3,specorigin
     $     ,specdectime
      real*8 specffermpx(nmax), specffermpy(nmax), specffermpz(nmax)
      common /specffermi/ specffermpx, specffermpy, specffermpz

      real*8 specfrr0(nmax), specfrrx(nmax), specfrry(nmax),
     &       specfrrz(nmax),specfrp0(nmax), specfrpx(nmax), 
     &       specfrpy(nmax), specfrpz(nmax)

      common /specfrcoor/ specfrr0, specfrrx, specfrry, specfrrz,
     &                    specfrp0, specfrpx, specfrpy, specfrpz        

      integer numspec,Zspec,Aspec,strangebe
      common /numspec/ numspec,Zspec,Aspec,strangebe

      real*8 specthydro

c now read particle-output
      nbar=0
      etotu=0.0d0
      pxtotu=0.0d0
      pytotu=0.0d0
      pztotu=0.0d0
      bartotu=0


        
c.. if spectators are on the hydro grid    
      if(CTOption(49).eq.1)then 
       do 39 i=1,npart
        r0(i)=h_r0(i)
        rx(i)=h_rx(i)
        ry(i)=h_ry(i)
        rz(i)=h_rz(i)
        p0(i)=h_p0(i)
        px(i)=h_px(i)
        py(i)=h_py(i)
        pz(i)=h_pz(i)
        fmass(i)=h_fmass(i)
        ityp(i)=h_ityp(i)
        iso3(i)=h_iso3(i)
        charge(i)=h_charge(i)
        lstcoll(i)=h_lstcoll(i)
        ncoll(i)=h_ncoll(i)
        origin(i)=h_origin(i)
        dectime(i)=h_tform(i)+dectim(i,0)
        tform(i)=h_tform(i)
        xtotfac(i)=h_xtotfac(i)
c.. get the freezeout coordinates right
        frr0(i) = h_tform(i)
        frrx(i) = h_rx(i)+(h_tform(i)-r0(i))*h_px(i)/h_p0(i)
        frry(i) = h_ry(i)+(h_tform(i)-r0(i))*h_py(i)/h_p0(i)
        frrz(i) = h_rz(i)+(h_tform(i)-r0(i))*h_pz(i)/h_p0(i)
        frp0(i) = h_p0(i)
        frpx(i) = h_px(i)
        frpy(i) = h_py(i)
        frpz(i) = h_pz(i)

        etotu=etotu+p0(i)
        pxtotu=pxtotu+abs(px(i))
        pytotu=pytotu+abs(py(i))
        pztotu=pztotu+abs(pz(i))
        if(abs(ityp(i)).lt.100)then
         bartotu=bartotu+ityp(i)/abs(ityp(i))      
        end if 
 
        if(abs(ityp(i)).le.maxbar)nbar=nbar+1
 39    continue
c.. default is seperate spectator propagation
      else 

c.. first the spectators
         do 41 j=1,numspec
            r0(j)=specr0(j)
            rx(j)=specrx(j)
            ry(j)=specry(j)
            rz(j)=specrz(j)
            p0(j)=specp0(j)
            px(j)=specpx(j)
            py(j)=specpy(j)
            pz(j)=specpz(j)
            fmass(j)=specfmass(j)
            ityp(j)=specityp(j)
            iso3(j)=speciso3(j)
            charge(j)=speccharge(j)
            lstcoll(j)=speclstcoll(j)
            ncoll(j)=specncoll(j)
            origin(j)=specorigin(j)
            dectime(j)=specdectime(j)
            tform(j)=spectform(j)
            xtotfac(j)=specxtotfac(j)
c..   get the freezeout coordinates right
            frr0(j) = specr0(j)
            frrx(j) = specrx(j)
            frry(j) = specry(j)
            frrz(j) = specrz(j)
            frp0(j) = specp0(j)
            frpx(j) = specpx(j)
            frpy(j) = specpy(j)
            frpz(j) = specpz(j)
c..   get the old fermi momenta
            ffermpx(j) = specffermpx(j)
            ffermpy(j) = specffermpy(j)
            ffermpz(j) = specffermpz(j)
            
            if(abs(ityp(j)).le.maxbar)nbar=nbar+1
            
 41      continue
         
c..   then all particles from cooperfrye
         do 40 i=numspec+1,numspec+npart
            r0(i)=h_r0(i-numspec)
            rx(i)=h_rx(i-numspec)
            ry(i)=h_ry(i-numspec)
            rz(i)=h_rz(i-numspec)
            p0(i)=h_p0(i-numspec)
            px(i)=h_px(i-numspec)
            py(i)=h_py(i-numspec)
            pz(i)=h_pz(i-numspec)
            fmass(i)=h_fmass(i-numspec)
            ityp(i)=h_ityp(i-numspec)
            iso3(i)=h_iso3(i-numspec)
            charge(i)=h_charge(i-numspec)
            lstcoll(i)=h_lstcoll(i-numspec)
            ncoll(i)=h_ncoll(i-numspec)
            origin(i)=h_origin(i-numspec)
            dectime(i)=h_tform(i-numspec)+dectim(i,0)
            tform(i)=h_tform(i-numspec)
            xtotfac(i)=h_xtotfac(i-numspec)
c..   get the freezeout coordinates right
            frr0(i) = h_tform(i-numspec)
            frrx(i) = h_rx(i-numspec)+(h_tform(i-numspec)-r0(i))
     &           *h_px(i-numspec)/h_p0(i-numspec)
            frry(i) = h_ry(i-numspec)+(h_tform(i-numspec)-r0(i))
     &           *h_py(i-numspec)/h_p0(i-numspec)
            frrz(i) = h_rz(i-numspec)+(h_tform(i-numspec)-r0(i))
     &           *h_pz(i-numspec)/h_p0(i-numspec)
            frp0(i) = h_p0(i-numspec)
            frpx(i) = h_px(i-numspec)
            frpy(i) = h_py(i-numspec)
            frpz(i) = h_pz(i-numspec)
 40      continue
c..   conservation check only for "hydro" part      
         do 42 i=numspec+1,numspec+npart
            etotu=etotu+p0(i)
            pxtotu=pxtotu+abs(px(i))
            pytotu=pytotu+abs(py(i))
            pztotu=pztotu+abs(pz(i))
            if(abs(ityp(i)).lt.100)then
               bartotu=bartotu+ityp(i)/abs(ityp(i))      
            end if 
            
            if(abs(ityp(i)).le.maxbar)nbar=nbar+1
 42      continue
      end if
c..   add spectators to npart
      npart=npart+numspec
      nmes=npart-nbar
      acttime=r0(1) 
      

      call arraysort2(1.0d0)
      
c     read options-file
c     call getparams
      write(*,301) 'fluidization  ','fluid','urqmd','diff'
      write(*,302) 'energy:           '
     +              ,h_etot,etotu,(etotu-h_etot)/h_etot*100
      write(*,304) 'baryons:          '
     +              ,h_btot,bartotu,(bartotu-h_btot)/h_btot*100
 301  format(a14,x,a5,x,a7,x,a6,x,'(conservation check)')
 302  format(a12,x,f7.1,x,f7.1,x,f5.1,'%')
 304  format(a12,x,f7.1,x,i7,x,f5.1,'%')

      if (abs((etotu-h_etot)/h_etot).gt.0.1d0)then
           write(*,*) 'energy difference is greater than 10 %'
           write(*,*) 'evolution is stopped'
      end if   

      return
c stop in case of EoF
 666  stop 137
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine arraysort2(thydro_start)
      implicit none

      include 'comres.f'
      include 'coms.f'
      include 'options.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'freezeout.f'
      include 'boxinc.f'
      integer i,j,b,m
      real*8 thydro_start

      
      
      real*8 br0(nmax),brx(nmax),bry(nmax),brz(nmax)
      real*8 bp0(nmax),bpx(nmax),bpy(nmax),bpz(nmax)
      real*8 bfmass(nmax),btform(nmax),bxtotfac(nmax),bdectime(nmax)
      integer bityp(nmax),biso3(nmax),bcharge(nmax),blstcoll(nmax)
      integer borigin(nmax),bncoll(nmax)

      real*8 mr0(nmax),mrx(nmax),mry(nmax),mrz(nmax)
      real*8 mp0(nmax),mpx(nmax),mpy(nmax),mpz(nmax)
      real*8 mfmass(nmax),mtform(nmax),mxtotfac(nmax),mdectime(nmax)
      integer mityp(nmax),miso3(nmax),mcharge(nmax),mlstcoll(nmax)
      integer morigin(nmax),mncoll(nmax)
       

      b=0
      m=0

      do i=1,npart 
         if(abs(ityp(i)).le.maxbar)then
            b=b+1
            br0(b)=r0(i)
            brx(b)=rx(i)
            bry(b)=ry(i)
            brz(b)=rz(i)
            bp0(b)=p0(i)
            bpx(b)=px(i)
            bpy(b)=py(i)
            bpz(b)=pz(i)
            
            bityp(b)=ityp(i)
            biso3(b)=iso3(i)
            bfmass(b)=fmass(i)        
            bcharge(b)=charge(i)
            blstcoll(b)=lstcoll(i)
            bncoll(b)=ncoll(i)
            borigin(b)=origin(i)
            bxtotfac(b)=xtotfac(i)  
            btform(b)=tform(i)
            bdectime(b)=dectime(i)
         else
            m=m+1
            mr0(m)=r0(i)
            mrx(m)=rx(i)
            mry(m)=ry(i)
            mrz(m)=rz(i)
            mp0(m)=p0(i)
            mpx(m)=px(i)
            mpy(m)=py(i)
            mpz(m)=pz(i)
            
            mityp(m)=ityp(i)
            miso3(m)=iso3(i)
            mfmass(m)=fmass(i)        
            mcharge(m)=charge(i)
            mlstcoll(m)=lstcoll(i)
            mncoll(m)=ncoll(i)
            morigin(m)=origin(i)
            mxtotfac(m)=xtotfac(i)  
            mtform(m)=tform(i)
            mdectime(m)=dectime(i)
         endif 
      end do
      do i=1,b
         r0(i)=br0(i) 
         rx(i)=brx(i)
         ry(i)=bry(i)
         rz(i)=brz(i)
         p0(i)=bp0(i)
         px(i)=bpx(i)
         py(i)=bpy(i)
         pz(i)=bpz(i)
         ityp(i)=bityp(i)
         iso3(i)=biso3(i)
         fmass(i)=bfmass(i)
         charge(i)=bcharge(i)
         lstcoll(i)=blstcoll(i)
         ncoll(i)=bncoll(i)
         origin(i)=borigin(i)
         xtotfac(i)=bxtotfac(i)
         tform(i)=btform(i)
         dectime(i)=bdectime(i)
      end do
      do i=b+1,b+m
         r0(i)=mr0(i-b) 
         rx(i)=mrx(i-b)
         ry(i)=mry(i-b)
         rz(i)=mrz(i-b)
         p0(i)=mp0(i-b)
         px(i)=mpx(i-b)
         py(i)=mpy(i-b)
         pz(i)=mpz(i-b)
         ityp(i)=mityp(i-b)
         iso3(i)=miso3(i-b)
         fmass(i)=mfmass(i-b)
         charge(i)=mcharge(i-b)
         lstcoll(i)=mlstcoll(i-b)
         ncoll(i)=mncoll(i-b)
         origin(i)=morigin(i-b)
         xtotfac(i)=mxtotfac(i-b)
         tform(i)=mtform(i-b)
         dectime(i)=mdectime(i-b)
      end do
  
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lorsphere(r,g,vx,vy,vz,m,pxmax,pxmin,
     &     pymax,pymin,pzmax,pzmin)
      implicit none
      real*8 r,g,vx,vy,vz,m,pxmax,pxmin,pymax,pymin,pzmax,pzmin
      real*8 p1x,p2x,p1y,p2y,p1z,p2z
      real*8 bqx,bqy,bqz
      real*8 pcx,pcy,pcz
      real*8 k,v2,qxbar,qybar,qzbar
      
      v2=vx**2+vy**2+vz**2
      
c     boost point (qx,qy,qz)=(0,0,r)
      
      p1x=bqx(0.d0,0.d0,r,g,vx,vy,vz,m)
      p1y=bqy(0.d0,0.d0,r,g,vx,vy,vz,m)
      p1z=bqz(0.d0,0.d0,r,g,vx,vy,vz,m)
      
c     boost point (qx,qy,qz)=(0,0,-r)
      p2x=bqx(0.d0,0.d0,-r,g,vx,vy,vz,m)
      p2y=bqy(0.d0,0.d0,-r,g,vx,vy,vz,m)
      p2z=bqz(0.d0,0.d0,-r,g,vx,vy,vz,m)
      
c     write(*,*)'p1z,p2z',p1z,p2z
      
      pzmin=p2z
      pzmax=p1z
      
c     find if bqz has an intrinsic local minimum
c     or maximum along the (0,0,qz) axis
      k=g*vz*v2/(v2+vz**2*(g-1.d0))
      if(abs(k).gt.1) then
         if(vz.lt.0) then
            qzbar=m/sqrt(k**2-1)
            if(qzbar.lt.r) then
               p1z=bqz(0.d0,0.d0,qzbar,g,vx,vy,vz,m)
               pzmax=p1z
            endif
         elseif(vz.gt.0) then
            qzbar=-m/sqrt(k**2-1)
            if(qzbar.gt.-r) then
               p2z=bqz(0.d0,0.d0,qzbar,g,vx,vy,vz,m)
               pzmin=p2z
            endif
         endif
      endif
      
c     boost point (qx,qy,qz)=(0,r,0)
      p1y=bqy(0.d0,r,0.d0,g,vx,vy,vz,m)


c     boost point (qx,qy,qz)=(0,-r,0)
      p2y=bqy(0.d0,-r,0.d0,g,vx,vy,vz,m)
      
      
      pymin=p2y                 !min(p2y,pymin)
      pymax=p1y                 !max(p1y,pymax)
      
c     find if bqy has an intrinsic local minimum
c     or maximum along the (0,qy,0) axis
      k=g*vy*v2/(v2+vy**2*(g-1.d0))
      if(abs(k).gt.1) then
         if(vy.lt.0) then
            qybar=m/sqrt(k**2-1)
            if(qybar.lt.r) then
               p1y=bqy(0.d0,qybar,0.d0,g,vx,vy,vz,m)
               pymax=max(pymax,p1y)
            endif
         elseif(vy.gt.0) then
            qybar=-m/sqrt(k**2-1)
            if(qybar.gt.-r) then
               p2y=bqy(0.d0,qybar,0.d0,g,vx,vy,vz,m)
               pymin=min(pymin,p2y)
            endif
         endif
      endif
      
c     boost point (qx,qy,qz)=(r,0,0)
      
      p1x=bqx(r,0.d0,0.d0,g,vx,vy,vz,m)
c     !      p1y=bqy(r,0.d0,0.d0,g,vx,vy,vz,m)
c     !      p1z=bqz(r,0.d0,0.d0,g,vx,vy,vz,m)
      
c     boost point (qx,qy,qz)=(-r,0,0)
      p2x=bqx(-r,0.d0,0.d0,g,vx,vy,vz,m)
c     !      p2y=bqy(-r,0.d0,0.d0,g,vx,vy,vz,m)
c     !      p2z=bqz(-r,0.d0,0.d0,g,vx,vy,vz,m)
      
      
      pxmin=p2x                 !min(p2x,pxmin)
c     !      pymin=min(min(p1y,p2y),pymin)
c     !      pzmin=min(min(p1z,p2z),pzmin)
      pxmax=p1x                 !max(p1x,pxmax)
c     !      pymax=max(max(p1y,p2y),pymax)
c     !      pzmax=max(max(p1z,p2z),pzmax)
      
c     find if bqx has an intrinsic local minimum
c     or maximum along the (qx,0,0) axis
      k=g*vx*v2/(v2+vx**2*(g-1.d0))
      if(abs(k).gt.1) then
         if(vx.lt.0) then
            qxbar=m/sqrt(k**2-1)
            if(qxbar.lt.r) then
               p1x=bqx(qxbar,0.d0,0.d0,g,vx,vy,vz,m)
               pxmax=max(pxmax,p1x)
            endif
         elseif(vx.gt.0) then
            qxbar=-m/sqrt(k**2-1)
            if(qxbar.gt.-r) then
               p2x=bqx(qxbar,0.d0,0.d0,g,vx,vy,vz,m)
               pxmin=min(pxmin,p2x)
            endif
         endif
      endif
      
      
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     returns boosted qx
      real*8 function bqx(qx,qy,qz,g,vx,vy,vz,m)
      implicit none
      real*8 qx,qy,qz,g,vx,vy,vz,m
      real*8 q0,v2,vsq
      
      q0=sqrt(m**2+qx**2+qy**2+qz**2)
      v2=vx**2+vy**2+vz**2
      vsq=vx*qx+vy*qy+vz*qz

      bqx=qx+(g-1.d0)/v2*vsq*vx+g*vx*q0
      
      
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       returns boosted qy
      real*8 function bqy(qx,qy,qz,g,vx,vy,vz,m)
      implicit none
      real*8 qx,qy,qz,g,vx,vy,vz,m
      real*8 q0,v2,vsq

      q0=sqrt(m**2+qx**2+qy**2+qz**2)
      v2=vx**2+vy**2+vz**2
      vsq=vx*qx+vy*qy+vz*qz

      bqy=qy+(g-1.d0)/v2*vsq*vy+g*vy*q0
      
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       returns boosted qz
      real*8 function bqz(qx,qy,qz,g,vx,vy,vz,m)
      implicit none
      real*8 qx,qy,qz,g,vx,vy,vz,m
      real*8 q0,v2,vsq
      
      q0=sqrt(m**2+qx**2+qy**2+qz**2)
      v2=vx**2+vy**2+vz**2
      vsq=vx*qx+vy*qy+vz*qz
      
      bqz=qz+(g-1.d0)/v2*vsq*vz+g*vz*q0
      
      return
      end
      
