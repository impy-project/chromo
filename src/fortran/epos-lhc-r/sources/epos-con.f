c-----------------------------------------------------------------------
      subroutine conaa(iret)
c-----------------------------------------------------------------------
c  determines interaction configuration
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      include 'epos.incsem'
      include 'epos.incpar'

      common/geom/rmproj,rmtarg,bmax,bkmx
      common/nucl3/phi,bimp
      common/cfacmss/facmss
c      double precision yyrmax
      common/scrangle/ phik3(kollmx),thetak3(kollmx)

      call utpri('conaa ',ish,ishini,4)

      iret=0

c     initialisations
c     ---------------

      vel=tanh(ypjtl-yhaha)+tanh(yhaha)
      rpro=2.*(r2had(iclpro)+slopom*log(engy*engy))*4*0.0389  !x2 to have ~1 for pp 7 TeV
      rtrg=2.*(r2had(icltar)+slopom*log(engy*engy))*4*0.0389

c     determine phi, bimp, coll, iproj, itarg, x/y/zproj, x/y/ztarg
c     ---------------------------------------------------------------

           if(iokoll.eq.1)then

      koll=matarg
      do k=1,koll
        do n=1,4
          coord(n,k)=0.
        enddo
      enddo
      bimp=0
      phi=0
      xproj(1)=0
      yproj(1)=0
      zproj(1)=0
      lproj(1)=koll
      lproj3(1)=0
      coltrp(1)=0d0
      ang=2*pi*rangen()
      diproj(1,1)=disize*cos(ang)*rpro
      diproj(1,2)=disize*sin(ang)*rpro
      ang=2*pi*rangen()
      ditarg(1,1)=disize*cos(ang)*rtrg
      ditarg(1,2)=disize*sin(ang)*rtrg
      do k=1,koll
        bij=bkmx*sqrt(rangen())
        bk(k)=bij
        iproj(k)=1
        itarg(k)=k
        coltrt(k)=0d0
        phi=2.*pi*rangen()
        xtarg(k)=bij*cos(phi)
        ytarg(k)=bij*sin(phi)
        ztarg(k)=0
        ltarg(k)=1
        kproj(1,k)=k
        ktarg(k,1)=k
        ltarg3(k)=0
        ktarg3(k,1)=0
        kproj3(k,1)=0
        if(iscreen.ne.0.and.bij.le.bkmxndif)then
c         if(zbrmax.le.0..or. bij.lt.zbcutx+zbrmax*rangen())then
c          if(zbrmax.le.0..or. rangen().lt.zbcutx/bij*zbrmax)then
          if(zbrads.gt.0.)then
            lproj3(1)=lproj3(1)+1
            ltarg3(k)=1
            kproj3(1,lproj3(1))=k
            ktarg3(k,1)=k
          endif
        endif
      enddo

           elseif(maproj.eq.1.and.matarg.eq.1)then

      b1=bminim
      b2=amin1(bkmx,bmaxim)
      if(b1.gt.b2)call utstop('conaa: bmin > bmax&')
      bimp=sqrt(b1*b1+(b2*b2-b1*b1)*rangen())
      koll=1
      do n=1,4
        coord(n,1)=0.
      enddo
      bk(1)=bimp
      iproj(1)=1
      itarg(1)=1
      phi=2.*pi*rangen()
      xproj(1)=bimp*cos(phi)
      yproj(1)=bimp*sin(phi)
      zproj(1)=0
      xtarg(1)=0
      ytarg(1)=0
      ztarg(1)=0
      lproj(1)=1
      ltarg(1)=1
      lproj3(1)=1
      ltarg3(1)=1
      kproj3(1,1)=1
      ktarg3(1,1)=1
      kproj(1,1)=1
      ktarg(1,1)=1
      coltrp(1)=dble(rangen()-zbrmax(iclpro))
      coltrt(1)=dble(rangen()-zbrmax(icltar))
c      coltrp(1)=dble(max(0.0001,1.+(rangen()-0.5)*zbrmax))
c      coltrt(1)=dble(max(0.0001,1.+(rangen()-0.5)*zbrmax))
      ang=2*pi*rangen()
      diproj(1,1)=disize*cos(ang)*rpro
      diproj(1,2)=disize*sin(ang)*rpro
      ang=2*pi*rangen()
      ditarg(1,1)=disize*cos(ang)*rtrg
      ditarg(1,2)=disize*sin(ang)*rtrg

           else

      call conxyz('p',mamx,xproj,yproj,zproj,ypjtl-yhaha)
      call conxyz('t',mamx,xtarg,ytarg,ztarg,yhaha)

      bx=0
      by=0
      if(maproj.gt.0)then
      if(bimevt.lt.0)then
        b1=bminim
        b2=amin1(rmproj+rmtarg,bmaxim)
        if(b1.gt.b2)call utstop('conaa: bmin > bmax&')
        bimp=sqrt(b1**2+(b2**2-b1**2)*rangen())
        if(nbarray.gt.0)bimp=barray(mod(nrevt,nbarray)+1)
        if(jpsi.gt.0)then
          bimp=b1+(b2-b1)*(float(mod(nrevt,12))+rangen())/12.
          bimevt=bimp
        endif
        phi=phimin+rangen()*(phimax-phimin)
      else
        phi=phievt
        bimp=bimevt
      endif
      bx=cos(phi)*bimp
      by=sin(phi)*bimp
      endif
      if(jpsi.lt.0)then !modify b
        bx=xtarg(1)
        by=ytarg(1)
      endif
      if(maproj.eq.0)goto1000
      koll=0
      do i=1,maproj
        lproj(i)=0
        lproj3(i)=0
        coltrp(i)=dble(rangen()-zbrmax(iclpro))!/maproj
c      coltrp(i)=dble(max(0.0001,1.+(rangen()-0.5)*zbrmax))
        !dipole:
        ang=2*pi*rangen()
        diproj(i,1)=disize*cos(ang)*rpro
        diproj(i,2)=disize*sin(ang)*rpro
      enddo
      do j=1,matarg
        ltarg(j)=0
        ltarg3(j)=0
        coltrt(j)=dble(rangen()-zbrmax(icltar))!/matarg
c      coltrt(j)=dble(max(0.0001,1.+(rangen()-0.5)*zbrmax))
        !dipole:
        ang=2*pi*rangen()
        ditarg(j,1)=disize*cos(ang)*rtrg
        ditarg(j,2)=disize*sin(ang)*rtrg
      enddo
c      rrr=rangen()
      do 12 i=1,maproj
      do 11 j=1,matarg
      if(jpsi.lt.0.and.ztarg(j).le.ztarg(1))goto 11
      bij=sqrt((xproj(i)+bx-xtarg(j))**2+(yproj(i)+by-ytarg(j))**2)
      if(ish.ge.7)write(ifch,*)'i_p:',i,' i_t:',j,' b_ij:',bij
      if(bij.gt.bkmx)goto 11

      koll=koll+1
      if(koll.gt.kollmx)call utstop('conaa: kollmx too small&')
      bk(koll)=bij
      bkx(koll)=xproj(i)+bx-xtarg(j)
      bky(koll)=yproj(i)+by-ytarg(j)
      iproj(koll)=i
      itarg(koll)=j
      lproj(i)=lproj(i)+1
      ltarg(j)=ltarg(j)+1
      kproj(i,lproj(i))=koll
      ktarg(j,ltarg(j))=koll
      phik3(koll)=0.
      thetak3(koll)=0.
      if(iscreen.ne.0.and.bij.le.bkmxndif)then
       if(zbrads.gt.0.)then
c        if(zbrmax.gt.0..and.bij.gt.zbcutx+zbrmax*rangen())goto 11
c          print *,'zbrmax',zbrmax,zbcutx,bij, zbcutx*zbrmax
c        if(zbrmax.gt.0..and. rrr*bij.ge.zbcutx*zbrmax)goto 11
        lproj3(i)=lproj3(i)+1
        ltarg3(j)=ltarg3(j)+1
        kproj3(i,lproj3(i))=koll
        ktarg3(j,ltarg3(j))=koll
c define angle for anti-shadowing
        if(abs(bky(koll)).gt.1.e-6)then
          if(abs(bkx(koll)).gt.1.e-6)then
            phik3(koll)=atan(bky(koll)/bkx(koll))
          else
            phik3(koll)=sign(0.5*pi,bky(koll))
          endif
        elseif(bkx(koll).lt.0.)then
          phik3(koll)=pi
        endif
        if(bk(koll).gt.0.)then
          thetak3(koll)=atan(bglaubx/bk(koll))
        else
          thetak3(koll)=0.5*pi
        endif
       endif
      endif

11    continue
12    continue

      do k=1,koll
        do n=1,4
          coord(n,k)=0.
        enddo
      enddo


          endif

      if(ish.ge.3)write(ifch,*)'koll=',koll
      if(koll.eq.0)goto 1001


c     determine coord
c     ---------------
      do kl=1,koll
      i=iproj(kl)
      j=itarg(kl)
      dist=ztarg(j)-zproj(i)
      coord(1,kl)=(xproj(i)+xtarg(j))*0.5
      coord(2,kl)=(yproj(i)+ytarg(j))*0.5
      coord(3,kl)=(zproj(i)+ztarg(j))*0.5
      coord(4,kl)=dist/vel
      enddo

      if(iscreen.ne.0.and.zbrads.gt.0.)call CalcScrPair(bimp)

    !~~~~~redefine energy in case of imposed radial flow~~~~~~~~~~~~~~~~
c      yrmaxi=max(0.,yradmx+yradmi*log(1.+engy/sqrt(float(koll))))
c      if(yrmaxi.gt.1e-5)then
c        yyrmax=dble(yrmaxi)
c        fradflii=sngl(1d0/
c     &  ((sinh(yyrmax)*yyrmax-cosh(yyrmax)+1d0)/(yyrmax**2/2d0)))
c      else
c        fradflii=1.
c      endif
c      if(ish.ge.3)write(ifch,*)'yrmaxi=',yrmaxi


c     exit
c     ----
1000  continue
      if(ish.ge.5)then
      write(ifch,*)'ia,x/y/zproj/coltrp:'
      do mm=1,maproj
      write(ifch,*)mm,xproj(mm),yproj(mm),zproj(mm),coltrp(mm)
      enddo
      write(ifch,*)'ia,x/y/ztarg/coltrt:'
      do mm=1,matarg
      write(ifch,*)mm,xtarg(mm),ytarg(mm),ztarg(mm),coltrt(mm)
      enddo
      write(ifch,*)'iret',iret
      endif
      call utprix('conaa ',ish,ishini,4)
      return

1001  continue !iret=1 causes redo of whole collision
      iret=1
      if(ish.ge.3)then
      write(ifch,*)
      write(ifch,*)'***** subroutine conaa:'
      write(ifch,*)'***** no nucleon pair found --> no interaction'
      write(ifch,*)
      endif
      goto 1000

      end

c-----------------------------------------------------------------------
      subroutine conigr(n,k)
c-----------------------------------------------------------------------
c  determines position of parton pair at position (n,k)
c-----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incems"      

      !KW210403
      !These are "effective" Pomeron positons, which may include initial random
      !walks between valence quarks and Pomeron end.
      !Parton positions are  very important for initial energy distribution asymmetry 
      !and therefore for everything related to flow, which is used to fix parameters.
      !Should be completely independent from the Metroplis stuff, and most importantly:
      !parameters should not be changed to improve the "non-hydro" part!!! 
      !It takes weeks to run full hydro to optimize these parameters.   

      !random:
      xi=0.
      yi=0.
      xj=0.
      yj=0.
c      call ranpospair(xi,yi,ri,  xj,yj,rj)
      i=iproj(k)
      j=itarg(k)
      !dipole:
      isi=1
      if(rangen().le.0.5)isi=-1
      dix=isi*diproj(i,1)*0.5
      diy=isi*diproj(i,2)*0.5
      isi=1 
      if(rangen().le.0.5)isi=-1
      djx=isi*ditarg(j,1)*0.5
      djy=isi*ditarg(j,2)*0.5
      !add random and dipole
      coordpr(1,n,k)=(xi+xj+dix+djx)*0.5+coord(1,k)
      coordpr(2,n,k)=(yi+yj+diy+djy)*0.5+coord(2,k)
      end

c-----------------------------------------------------------------------
      function frdmzcut(b)
c-----------------------------------------------------------------------
c  determines interaction configuration
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incpar'

c      frdmzcut=zbcutx!*exp(-(b*rangen()))
      frdmzcut=max(0.,b+zbcutx)
c      bdum=b
      return
      end

c-----------------------------------------------------------------------
      subroutine CalcScrPair(b)
c-----------------------------------------------------------------------
c  determines interaction configuration
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      include 'epos.incpar'

      common/scrangle/ phik3(kollmx),thetak3(kollmx)
      logical lascr(kollmx),cont


c     order pairs for splitting as a function of b
c     ---------------
      do i=1,maproj
        if(lproj3(i).gt.1)then
          do l=2,lproj3(i)
            m=l
 100        continue
              kp=kproj3(i,m-1)
              kl=kproj3(i,m)
              if(bk(kl).lt.bk(kp))then
                kproj3(i,m-1)=kl
                kproj3(i,m)=kp
                m=m-1
            if(m.gt.1)goto 100
              endif
          enddo
        endif
      enddo
      do j=1,matarg
        if(ltarg3(j).gt.1)then
          do l=2,ltarg3(j)
            m=l
 200        continue
              kp=ktarg3(j,m-1)
              kl=ktarg3(j,m)
              if(bk(kl).lt.bk(kp))then
                ktarg3(j,m-1)=kl
                ktarg3(j,m)=kp
                m=m-1
            if(m.gt.1)goto 200
              endif
          enddo
        endif
      enddo
      return   !put shadowing in GfunPark

      
      if(koll.gt.1)then
c Define anti-shadowing as a consequence of geometrical screening of
c nucleons at large b by the one at small b

c Projectile
c Check phi not to be in the range of a nucleon with smaller b
        do i=1,maproj
          if(lproj3(i).gt.1)then
            kl=kproj3(i,1)
            lascr(kl)=.true.
            do 300 l=lproj3(i),2,-1
              kl=kproj3(i,l)
              lascr(kl)=.true.
              if(bk(kl).ge.frdmzcut(b))then
                do m=1,l-1
                  km=kproj3(i,m)
                  if(kl.ne.km.and.bk(km).le.frdmzcut(b)
     &          .and.phik3(kl).ge.phik3(km)-thetak3(km)
     &          .and.phik3(kl).le.phik3(km)+thetak3(km))then
                    lascr(kl)=.false.
                    goto 300
                  endif
                enddo
              endif
 300        continue
          endif
        enddo
c suppress screened pair from the list
        do i=1,maproj
          if(lproj3(i).gt.1)then
            do l=1,lproj3(i)
              kl=kproj3(i,l)
            enddo
          endif
        enddo
        do i=1,maproj
          if(lproj3(i).gt.1)then
            n=2
            cont=lproj3(i).gt.1
            do while(cont)
              l=lproj3(i)
              kl=kproj3(i,l)
c suppress end of the list in order to have the last pair active
              do while(.not.lascr(kl).and.l.gt.1)
                kproj3(i,l)=0
                lproj3(i)=lproj3(i)-1
                l=lproj3(i)
                kl=kproj3(i,l)
              enddo
              cont=lproj3(i).gt.n
              if(cont)then
c compress list
                kn=kproj3(i,n)
                if(.not.lascr(kn))then
                  m=lproj3(i)   !last pair always active
                  km=kproj3(i,m)
                  kproj3(i,n)=km
                  kproj3(i,m)=0
                  lproj3(i)=lproj3(i)-1
                endif
              endif
              n=min(lproj3(i)-1,n)+1
              cont=lproj3(i).ne.n
            enddo
          endif
        enddo

c Target
c Check phi not to be in the range of a nucleon with smaller b
        do j=1,matarg
          if(ltarg3(j).gt.1)then
            kl=ktarg3(j,1)
            lascr(kl)=.true.
            do 400 l=ltarg3(j),2,-1
              kl=ktarg3(j,l)
              lascr(kl)=.true.
              if(bk(kl).le.frdmzcut(b))then
                do m=1,l-1
                  km=ktarg3(j,m)
                  if(km.ne.kl.and.bk(km).le.frdmzcut(b)
     &          .and.phik3(kl).ge.phik3(km)-thetak3(km)
     &          .and.phik3(kl).le.phik3(km)+thetak3(km))then
                    lascr(kl)=.false.
                    goto 400
                  endif
                enddo
              endif
 400        continue
          endif
        enddo
c suppress screened pair from the list
        do j=1,matarg
          if(ltarg3(j).gt.1)then
            n=2
            cont=ltarg3(j).gt.1
            do while(cont)
              l=ltarg3(j)
              kl=ktarg3(j,l)
c suppress end of the list in order to have the last pair active
              do while(.not.lascr(kl).and.l.gt.1)
                ktarg3(j,l)=0
                ltarg3(j)=ltarg3(j)-1
                l=ltarg3(j)
                kl=ktarg3(j,l)
              enddo
              cont=ltarg3(j).gt.n
              if(cont)then
c compress list
                kn=ktarg3(j,n)
                if(.not.lascr(kn))then
                  m=ltarg3(j)   !last pair always active
                  km=ktarg3(j,m)
                  ktarg3(j,n)=km
                  ktarg3(j,m)=0
                  ltarg3(j)=ltarg3(j)-1
                endif
              endif
              n=min(ltarg3(j)-1,n)+1
              cont=ltarg3(j).ne.n
            enddo
          endif
        enddo
      endif


      end

c-----------------------------------------------------------------------
      subroutine xGeometry(iii)
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      include 'epos.incsem'
      common/xgeom/nnn,naa(kollmx),nbb(kollmx)
      character*5 fmt1,fmt2

      if(iii.eq.0)then

      do k=1,kollmx
        naa(k)=0
      enddo
      nnn=0


      elseif(iii.eq.1)then

      ngl=0
      do k=1,koll
        r=bk(k)
        if(r.le.sqrt(sigine/10./pi))ngl=ngl+1
      enddo
      if(ngl.ne.0)then
        nnn=nnn+1
        naa(ngl)=naa(ngl)+1
      endif

      elseif(iii.eq.2)then

      if(xpar1.eq.0..and.xpar2.eq.0.)then
       print*,'---------- xGeometry -----------------------'
       return
      endif
      x1=1-0.01*xpar2
      x2=1-0.01*xpar1
      kmx=0
      nbb(1)=naa(1)
      do k=2,kollmx
        if(naa(k).ne.0.)kmx=k
        nbb(k)=nbb(k-1)+naa(k)
      enddo
      k1=0
      k2=0
      do k=1,kmx
        x=nbb(k)/float(nnn)
        if(x.lt.x1)k1=k
        if(x.lt.x2)k2=k
      enddo
      k1=k1+1
      k2=k2+1
      x=0
      av=0
      su=0
      p1=0.
      p2=0.
      do k=1,kmx
        xb=x
        x=nbb(k)/float(nnn)
        y=naa(k)/float(nnn)
        dx=x-xb
        p=0.
        if(k.eq.k1)then
          p=(x-x1)/dx
          p1=p
        elseif(k.eq.k2)then
          p=(x2-xb)/dx
          p2=p
        elseif(k.gt.k1.and.k.lt.k2)then
          p=1
        endif
      av=av+y*p*k
      su=su+y*p
      enddo
      av=av/su
      n1=nint(100*p1)
      n2=nint(100*p2)
      if(n1.eq.0)then
        k1=k1+1
        n1=100
      endif
      if(k1.le.9)fmt1='i1,4x'
      if(k1.gt.9.and.k1.le.99)fmt1='i2,3x'
      if(k1.gt.99.and.k1.le.999)fmt1='i3,2x'
      if(k1.gt.999.and.k1.le.9999)fmt1='i4,1x'
      if(k2.le.9)fmt2='i1,4x'
      if(k2.gt.9.and.k2.le.99)fmt2='i2,3x'
      if(k2.gt.99.and.k2.le.999)fmt2='i3,2x'
      if(k2.gt.999.and.k2.le.9999)fmt2='i4,1x'
      write(6,'(i4,a,i5,a,'//fmt1//',i6,a,i5,a,'//fmt2//',5x,a,f8.2)')
     &       nint(xpar2),':MIN',n1,'%',k1
     &      ,nint(xpar1),':MAX',n2,'%',k2   ,'av:',av
      endif

      end

c-----------------------------------------------------------------------
      function conbmx()
c-----------------------------------------------------------------------
      double precision om1intbc,p,eps
      include 'epos.inc'
      include 'epos.incsem'

      conbmx=0.
      b1=0.
      b2=10.
      eps=5.0d-3
      p=1.d0-dexp(-om1intbc(b2))
      if(p.gt.2.d0*eps)return

      ntry=0

10    ntry=ntry+1
      b=b1+.5*(b2-b1)

      p=(1.d0-dexp(-om1intbc(b)))
      
      if(p.gt.eps)then
       if(p.gt.2.d0*eps)then
        b1=b
       else
        conbmx=b
       return
       endif
      else
       if(p.lt.eps/5.d0)then
        b2=b
       else
        conbmx=b
        return
       endif
      endif

      if(ntry.le.1000)goto 10
      write(ifmt,*)'Too much try in conbmx ... bmax=',b
      conbmx=b
      return

      end

c-----------------------------------------------------------------------
      function conbmxndif()
c-----------------------------------------------------------------------
      double precision om1intbc,p,eps
      include 'epos.inc'
      include 'epos.incsem'


      iomegasave=iomega
      iomega=2
      conbmxndif=0.
      b1=0.
      b2=10.
      conbmxndif=b2
      eps=1d-10
      p=1.d0-dexp(-om1intbc(b2))
      if(p.gt.2.d0*eps)goto 100

      ntry=0

10    ntry=ntry+1
      b=b1+.5*(b2-b1)

      p=(1.d0-dexp(-om1intbc(b)))

      if(p.gt.eps)then
       if(p.gt.2.d0*eps)then
        b1=b
       else
        conbmxndif=b
       goto 100
       endif
      else
       if(p.lt.eps/5.d0)then
        b2=b
       else
        conbmxndif=b
        goto 100
       endif
      endif

      if(ntry.le.1000)goto 10
      write(ifmt,*)'Too much try in conbmxndif ... bkmxndif=',b
      conbmxndif=b
 100  iomega=iomegasave
      return

      end

c-----------------------------------------------------------------------
      function conbmxdif()
c-----------------------------------------------------------------------
c find b to have (1-exp(-om))pmax=pdiff
c-----------------------------------------------------------------------
      double precision om1intbc,pmax,drootom,pdiff
      include 'epos.inc'
      include 'epos.incsem'

      conbmxdif=0.
      b1=0.
      bmax=7.
      iomegasave=iomega
      iomega=2

      eps=1.e-5
      pmax=1.d0-dexp(-om1intbc(b1))
      pdiff=facdif
      if(pmax.gt.eps)then
        conbmxdif=drootom(pdiff,pmax,bmax,eps)
      endif
      iomega=iomegasave

      return

      end

c-----------------------------------------------------------------------
      subroutine conre
c-----------------------------------------------------------------------
c     initializes remnants taking into account the fact that neutron and
c     proton should be closest neighbor
c-----------------------------------------------------------------------
      include "epos.incems"
      include 'epos.inc'
      logical lneighbor,lset(mamx)
      dimension distlk(mamx,mamx)

      call utpri('conre ',ish,ishini,6)

      lneighbor=.true.

c     proj
c     ----
      la=laproj
      ma=iabs(maproj)
      las=0
      mas=0
      do l=1,ma
        lset(l)=.false.
      enddo
      if(la.gt.0)then
        if(lneighbor)then
          do l=1,ma
            do k=1,ma
              distlk(l,k)=-1.
            enddo
          enddo
          mstart=ma-2*la
          if(mstart.gt.0)then
            idstart=1220
          elseif(mstart.lt.0)then
            idstart=1120
            mstart=abs(mstart)
          else
            mstart=1
            idstart=1220
            if(rangen().le.(la-las)*1./(ma-mas))id=1120
          endif
        else
          idstart=1220
          mstart=ma
        endif
      else
        mstart=1
        idstart=idproj
        las=1
        mas=ma
      endif
      do l=1,mstart
      if(la.lt.0)then
        if(iabs(idproj).lt.20)then
          id=idproj
        else
          ia=iabs(idproj/10)
          is=idproj/iabs(idproj)
          if(ia.ne.111.and.ia.ne.222.and.ia.ne.333)id=idproj/10*10
          if(ia.eq.111.or. ia.eq.222.or. ia.eq.333)id=idproj/10*10+is
          if(ia.eq.213)id=1230*is
          if(ia.eq.214)id=1240*is
        endif
      else
c assign to the first (randomly distributed) nucleons some id
        id=sign(idstart,idproj)
        if(lneighbor)then
          lset(l)=.true.
        else
          if(rangen().le.(la-las)*1./(ma-mas))id=sign(1120,idproj)
        endif
        if(abs(id).eq.1120)las=las+1
        mas=mas+1
      endif
      ic1=idtrai(1,id,1)
      ic2=idtrai(2,id,1)
      icproj(1,l)=ic1
      icproj(2,l)=ic2
      enddo

      if(lneighbor)then
c       look for closest neighbor to assign different id
        l=0
        do while(mas.lt.ma)
          l=l+1
          if(lset(l))then
            k=l+1
            if(k.gt.ma)k=1
            distmin=ainfin
            kmin=0
            do while(k.ne.l)
              if(.not.lset(k))then
                dist=distlk(l,k)
                if(dist.lt.0.)dist=(xproj(k)-xproj(l))**2
     &             +(yproj(k)-yproj(l))**2+(zproj(k)-zproj(l))**2
                if(dist.le.distmin)then
                  distmin=dist
                  kmin=k
                endif
              endif
              k=k+1
              if(k.gt.ma)k=1
            enddo
            if(kmin.eq.0)call utstop("kmin should not be undefined !&")
            lset(kmin)=.true.
            mas=mas+1
            ic1=0
            ic2=0
            if(icproj(1,l).eq.210000)then
              ic1=120000
            elseif(icproj(1,l).ne.0)then
              ic1=210000
              las=las+1
            endif
            if(icproj(2,l).eq.210000)then
              ic2=120000
            elseif(icproj(2,l).ne.0)then
              ic2=210000
              las=las+1
            endif
            icproj(1,kmin)=ic1
            icproj(2,kmin)=ic2
            l=kmin-1
          endif
          if(l.eq.ma)l=0
        enddo
      endif
c     targ
c     ----
      la=latarg
      ma=iabs(matarg)
      las=0
      mas=0
      do l=1,ma
        lset(l)=.false.
      enddo
      if(la.gt.0)then
        if(lneighbor)then
          do l=1,ma
            do k=1,ma
              distlk(l,k)=-1.
            enddo
          enddo
          mstart=ma-2*la
          if(mstart.gt.0)then
            idstart=1220
          elseif(mstart.lt.0)then
            idstart=1120
            mstart=abs(mstart)
          else
            mstart=1
            idstart=1220
            if(rangen().le.(la-las)*1./(ma-mas))id=1120
          endif
        else
          idstart=1220
          mstart=ma
        endif
      else
        mstart=1
        idstart=idtarg
        las=1
        mas=ma
      endif
      do l=1,mstart
      if(la.lt.0)then
        if(iabs(idtarg).lt.20)then
          id=idtarg
        else
          ia=iabs(idtarg/10)
          is=idtarg/iabs(idtarg)
          if(ia.ne.111.and.ia.ne.222.and.ia.ne.333)id=idtarg/10*10
          if(ia.eq.111.or. ia.eq.222.or. ia.eq.333)id=idtarg/10*10+is
          if(ia.eq.213)id=1230*is
        endif
      else
        id=idstart
        if(lneighbor)then
          lset(l)=.true.
        else
          if(rangen().le.(la-las)*1./(ma-mas))id=1120
        endif
        if(id.eq.1120)las=las+1
        mas=mas+1
      endif
      ic1=idtrai(1,id,1)
      ic2=idtrai(2,id,1)
      ictarg(1,l)=ic1
      ictarg(2,l)=ic2
      enddo

      if(lneighbor)then
c       look for closest neighbor to assign different id
        l=0
        ic2=0
        do while(mas.lt.ma)
          l=l+1
          if(lset(l))then
            k=l+1
            if(k.gt.ma)k=1
            distmin=ainfin
            kmin=0
            do while(k.ne.l)
              if(.not.lset(k))then
                dist=distlk(l,k)
                if(dist.lt.0.)dist=(xtarg(k)-xtarg(l))**2
     &             +(ytarg(k)-ytarg(l))**2+(ztarg(k)-ztarg(l))**2
                if(dist.le.distmin)then
                  distmin=dist
                  kmin=k
                endif
              endif
              k=k+1
              if(k.gt.ma)k=1
            enddo
            if(kmin.eq.0)call utstop("kmin should not be undefined !&")
            lset(kmin)=.true.
            mas=mas+1
            if(ictarg(1,l).eq.210000)then
              ic1=120000
            else
              ic1=210000
              las=las+1
            endif
            ictarg(1,kmin)=ic1
            ictarg(2,kmin)=ic2
            l=kmin-1
          endif
          if(l.eq.ma)l=0
        enddo
      endif

      call utprix('conre ',ish,ishini,6)
      return
      end

c-----------------------------------------------------------------------
      subroutine conrl
c-----------------------------------------------------------------------
c  initializes target remnant in case of appl lepton
c-----------------------------------------------------------------------
      include "epos.incems"
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
      common/hadr2/iomodl,idproj,idtarg,wexcit

c     targ
c     ----
      la=latarg
      ma=iabs(matarg)
      las=0
      mas=0
      do l=1,ma
      if(la.lt.0)then
      id=idtarg
      else
      id=1220
      if(rangen().le.(la-las)*1./(ma-mas))id=1120
      if(id.eq.1120)las=las+1
      mas=mas+1
      endif
      ic1=idtrai(1,id,1)
      ic2=idtrai(2,id,1)
      ictarg(1,l)=ic1
      ictarg(2,l)=ic2
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine conwr
c-----------------------------------------------------------------------
c     writes /cptl/
c-----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incems"
      double precision XA(64,3),XB(64,3),BQGS,BMAXQGS,BMAXNEX,BMINNEX
      COMMON /Q_QGSNEX1/ XA,XB,BQGS,BMAXQGS,BMAXNEX,BMINNEX
      common/nucl3/phi,bimp
      parameter(iapmax=208)
      double precision bqgs2,bmaxqgs2,bmaxnex2,bminnex2,xan,xbn
      common /qgsIInex1/xan(iapmax,3),xbn(iapmax,3)
     *,bqgs2,bmaxqgs2,bmaxnex2,bminnex2
      double precision bqgs3,bmaxqgs3,bmaxnex3,bminnex3,xan3,xbn3
      common /qgsIIInex1/xan3(iapmax,3),xbn3(iapmax,3)
     *,bqgs3,bmaxqgs3,bmaxnex3,bminnex3
      common/photrans/phoele(4),ebeam
      integer ic(2)

      call utpri('conwr ',ish,ishini,6)

      bx=cos(phi)*bimp
      by=sin(phi)*bimp

c     write /cptl/
c     ------------
      nptl=0

           if(iokoll.eq.1)then   ! precisely matarg collisions

      nptl=nptl+1
      do 3 i=1,4
3     xorptl(i,nptl)=0
      tivptl(1,nptl)=-ainfin
      tivptl(2,nptl)=0
      istptl(nptl)=1
      iorptl(nptl)=-1
      jorptl(nptl)=0
      do 1 k=1,koll
      nptl=nptl+1
      do 4 i=1,4
4     xorptl(i,nptl)=0
      tivptl(1,nptl)=-ainfin
      tivptl(2,nptl)=0
      istptl(nptl)=1
      iorptl(nptl)=-1
      jorptl(nptl)=0
1     continue

           elseif(iappl.ne.7)then

c             print *,'proj'

      do 6 i=1,maproj
      nptl=nptl+1
      istptl(nptl)=0
      iorptl(nptl)=0
      jorptl(nptl)=0
      if(model.eq.2)then       !QGSJet
      xproj(i)=XA(i,1)
      yproj(i)=XA(i,2)
      zproj(i)=XA(i,3)
      istptl(nptl)=1
      iorptl(nptl)=-1
      elseif(model.eq.7.or.model.eq.11)then       !QGSJetII
      xproj(i)=xan(i,1)
      yproj(i)=xan(i,2)
      zproj(i)=xan(i,3)
      istptl(nptl)=1
      iorptl(nptl)=-1
      elseif(model.eq.13)then       !QGSJetIII
      xproj(i)=xan3(i,1)
      yproj(i)=xan3(i,2)
      zproj(i)=xan3(i,3)
      istptl(nptl)=1
      iorptl(nptl)=-1
      elseif(model.ge.3)then       !Gheisha, ...
      istptl(nptl)=1
      iorptl(nptl)=-1
      endif
      xorptl(1,nptl)=xproj(i)+bx/2
      xorptl(2,nptl)=yproj(i)+by/2
      xorptl(3,nptl)=zproj(i)
      xorptl(4,nptl)=0
      tivptl(1,nptl)=-ainfin
c for visualisation uncomment
c-c   tivptl(1,nptl)=-100
      tivptl(2,nptl)= ainfin
c             print *,i,xorptl(1,nptl),xorptl(2,nptl),xorptl(3,nptl)
6     continue
c             print *,'targ'
      do 7 i=1,matarg
      nptl=nptl+1
      istptl(nptl)=0
      iorptl(nptl)=0
      jorptl(nptl)=0
      if(model.eq.2)then       !QGSJet
      xtarg(i)=XB(i,1)
      ytarg(i)=XB(i,2)
      ztarg(i)=XB(i,3)
      istptl(nptl)=1
      iorptl(nptl)=-1
      elseif(model.eq.7.or.model.eq.11)then       !QGSJetII
      xtarg(i)=xbn(i,1)
      ytarg(i)=xbn(i,2)
      ztarg(i)=xbn(i,3)
      istptl(nptl)=1
      iorptl(nptl)=-1
      elseif(model.eq.13)then       !QGSJetIII
      xtarg(i)=xbn3(i,1)
      ytarg(i)=xbn3(i,2)
      ztarg(i)=xbn3(i,3)
      istptl(nptl)=1
      iorptl(nptl)=-1
      elseif(model.ge.3)then       !Gheisha, ...
      istptl(nptl)=1
      iorptl(nptl)=-1
      endif
      xorptl(1,nptl)=xtarg(i)-bx/2
      xorptl(2,nptl)=ytarg(i)-by/2
      xorptl(3,nptl)=ztarg(i)
      xorptl(4,nptl)=0
      tivptl(1,nptl)=-ainfin
c for visualisation uncomment
c-c   tivptl(1,nptl)=-100
      tivptl(2,nptl)= ainfin
c             print *,i,xorptl(1,nptl),xorptl(2,nptl),xorptl(3,nptl)
7     continue
      if(abs(idprojin).eq.12.and..not.chargex)then   !electron for fake DIS
c electron projectile
        nptl=nptl+1
        istptl(nptl)=41
        iorptl(nptl)=-1
        jorptl(nptl)=-1
        iorptl(1)=nptl         !projectile coming from lepton
        xorptl(1,nptl)=bx/2
        xorptl(2,nptl)=by/2
        xorptl(3,nptl)=0.
        xorptl(4,nptl)=0.
        tivptl(1,nptl)=-ainfin
        tivptl(2,nptl)=0.
c target nucleons (in lab frame)
        do i=1,matarg
          nptl=nptl+1
          istptl(nptl)=41
          iorptl(nptl)=-1
          jorptl(nptl)=-1
          xorptl(1,nptl)=xtarg(i)-bx/2
          xorptl(2,nptl)=ytarg(i)-by/2
          xorptl(3,nptl)=ztarg(i)
          xorptl(4,nptl)=0
          tivptl(1,nptl)=-ainfin
c         for visualisation uncomment
c         -c   tivptl(1,nptl)=-100
          tivptl(2,nptl)= ainfin
        enddo
c electron remnant
        nptl=nptl+1
        istptl(nptl)=0
        iorptl(nptl)=maproj+matarg+1
        jorptl(nptl)=-1
        xorptl(1,nptl)=bx/2
        xorptl(2,nptl)=by/2
        xorptl(3,nptl)=0.
        xorptl(4,nptl)=0.
        tivptl(1,nptl)=0.
        tivptl(2,nptl)= ainfin
      endif

          endif

      nptl=0
      if(iappl.le.2)then
      do i=1,maproj
      nptl=nptl+1
      ic(1)=icproj(1,i)
      ic(2)=icproj(2,i)
      id=idtra(ic,0,0,0)
      if(laproj.eq.-1)id=idproj
c      id=idtra(ic,0,0,3)      !tp071107 imix=3 ??????????
      call idmass(id,ams)
      idptl(nptl)=id
      pptl(1,nptl)=0.
      pptl(2,nptl)=0.
      pptl(3,nptl)=pnullx
      pptl(4,nptl)=sqrt(pnullx**2+ams**2)
      pptl(5,nptl)=ams
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      ityptl(nptl)=0
      enddo
      endif
      if(iappl.ne.7)then
      do i=1,matarg
      nptl=nptl+1
      ic(1)=ictarg(1,i)
      ic(2)=ictarg(2,i)
      id=idtra(ic,0,0,0)
      if(latarg.eq.-1)id=idtarg
c      id=idtra(ic,0,0,3)      !tp071107 imix=3 ??????????
      call idmass(id,ams)
      idptl(nptl)=id
      pptl(1,nptl)=0.
      pptl(2,nptl)=0.
      pptl(3,nptl)=-pnullx
      pptl(4,nptl)=sqrt(pnullx**2+ams**2)
      pptl(5,nptl)=ams
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      ityptl(nptl)=0
      enddo
      if(abs(idprojin).eq.12.and..not.chargex)then   !electron for fake DIS
c electron projectile
        nptl=nptl+1
        id=idprojin
        call idmass(id,ams)
        idptl(nptl)=id
        pptl(1,nptl)=0.
        pptl(2,nptl)=0.
        pptl(3,nptl)=sqrt(max(0.,(elepti+ams)*(elepti-ams)))
        pptl(4,nptl)=elepti
        pptl(5,nptl)=ams
        ifrptl(1,nptl)=1
        ifrptl(2,nptl)=1
        ityptl(nptl)=40
c target nucleons (in lab frame)
        do i=1,matarg
          nptl=nptl+1
          idptl(nptl)=idptl(maproj+i)
          pptl(1,nptl)=0.
          pptl(2,nptl)=0.
          pptl(3,nptl)=-pnll
          pptl(4,nptl)=ebeam
          pptl(5,nptl)=pptl(5,maproj+i)
          ifrptl(1,nptl)=maproj+i
          ifrptl(2,nptl)=maproj+i
          ityptl(nptl)=50
        enddo
c electron remnant
        nptl=nptl+1
        idptl(nptl)=id
        pptl(1,nptl)=phoele(1)
        pptl(2,nptl)=phoele(2)
        pptl(3,nptl)=phoele(3)
        pptl(4,nptl)=phoele(4)
        pptl(5,nptl)=ams
        ifrptl(1,nptl)=0
        ifrptl(2,nptl)=0
        ityptl(nptl)=40
      endif

      else

      nptl=nptl+1
      id=idproj
      call idmass(id,ams)
      idptl(nptl)=id
      pptl(1,nptl)=0.
      pptl(2,nptl)=0.
      pptl(3,nptl)=pnullx
      pptl(4,nptl)=sqrt(pnullx**2+ams**2)
      pptl(5,nptl)=ams
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      ityptl(nptl)=0
      iorptl(nptl)=-1
      jorptl(nptl)=0
      istptl(nptl)=0
      do 5 i=1,4
 5      xorptl(i,nptl)=0
      tivptl(1,nptl)=0
      tivptl(2,nptl)=0
      endif

c     exit
c     ----

      call utprix('conwr ',ish,ishini,6)
      return
      end

c------------------------------------------------------------------------
      subroutine conxyz(ch,n,x,y,z,ynuc)
c-----------------------------------------------------------------------
      include 'epos.inc'

      real x(n),y(n),z(n)
      character ch*1

      massnr=0
      iii=0
      if(ch.eq.'p')then
      massnr=maproj
      iii=1
      elseif(ch.eq.'t')then
      massnr=matarg
      iii=2
      else
      call utstop('conxyz: nucleus neither proj nor targ&')
      endif

      if(massnr.eq.0)return
      if(massnr.gt.n)call utstop('conxyz: massnr.gt.n&')
      if(massnr.eq.1)then
      x(1)=0
      y(1)=0
      z(1)=0
      return
      endif

      rad=radnuc(massnr)

      if(massnr.ge.10)then !---wood-saxon density---

        rad=rad/difnuc(massnr)
        cr1=1.+3./rad+6./rad**2+6./rad**3
        cr2=3./rad
        cr3=3./rad+6./rad**2
        do i=1,massnr
   1      zuk=rangen()*cr1-1.
          if(zuk.le.0.)then
            tt=rad*(rangen()**.3333-1.)
          elseif(zuk.le.cr2 )then
            tt=-log(rangen())
          elseif(zuk.lt.cr3 )then
            tt=-log(rangen())-log(rangen())
          else
            tt=-log(rangen())-log(rangen())-log(rangen())
          endif
          if(rangen().gt.1./(1.+exp(-abs(tt))))goto 1
          rim=tt+rad
          zz=rim*(2.*rangen()-1.)
          rim=sqrt(rim*rim-zz*zz)
          z(i)=zz*difnuc(massnr)
          call pscs(c,s)
          x(i)=rim*c*difnuc(massnr)
          y(i)=rim*s*difnuc(massnr)
        enddo

      elseif(massnr.ge.3)then  ! ---gaussian density---

        rad=rad*sqrt(2.*massnr/(massnr-1.))   !van hove simulation
        do l=1,3
          summ=0.
          do i=1,massnr-1
            j=massnr-i
            aks=rad *(rangen()+rangen()+rangen()-1.5)
            k=j+1
            if(l.eq.1)x(k)=summ-aks*sqrt(float(j)/k)
            if(l.eq.2)y(k)=summ-aks*sqrt(float(j)/k)
            if(l.eq.3)z(k)=summ-aks*sqrt(float(j)/k)
            summ=summ+aks/sqrt(float(j*k))
          enddo
          if(l.eq.1)x(1)=summ
          if(l.eq.2)y(1)=summ
          if(l.eq.3)z(1)=summ
        enddo

      elseif(massnr.eq.2)then  ! ---deuteron---

        !.........r=t*difnuc(massnr), t~exp(-2*t)*(1-exp(-a*t))
        a=radnuc(massnr)
  2     t=-0.5*alog(rangen())  !~exp(-2*t)
        if(rangen().gt.(1-exp(-a*t))**2)goto2
        r=t*difnuc(massnr)
        zz=r*(2.*rangen()-1.)
        call pscs(c,s)
        rxy=sqrt(r*r-zz*zz)
        z(1)=0.5*zz
        x(1)=0.5*rxy*c
        y(1)=0.5*rxy*s
        z(2)=-z(1)
        x(2)=-x(1)
        y(2)=-y(1)

      else

        stop'conxyz: wrong massnr.     '

      endif

c...plot preparation

      rmax=(radnuc(massnr)+3)
      drnucl(iii)=rmax/mxnucl
      nrnucl(iii)=nrnucl(iii)+1
      do i=1,massnr
        r=sqrt(x(i)**2+y(i)**2+z(i)**2)
        k=1+int(r/drnucl(iii))
        if(k.le.mxnucl)rnucl(k,iii)=rnucl(k,iii)+1
      enddo

c...lorentz trafo

      do i=1,massnr
      z(i)=z(i)/cosh(ynuc)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine conini
c-----------------------------------------------------------------------
      include 'epos.inc'

      imax=max(maproj,matarg)
      if(idtargin.eq.0)imax=max(imax,40)
      do massnr=1,mamxx
        dif=0.54
        rad=0.
        if(massnr.gt.imax.or.massnr.eq.1)then
          dif=0
          rad=0
        elseif(massnr.eq.197)then
          dif=0.562
          rad=6.5
        elseif(massnr.ge.10)then
          rad=1.12*massnr**0.33333-0.86*massnr**(-0.33333)
        elseif(massnr.ge.3)then
          rad=.9*float(massnr)**.3333
        elseif(massnr.eq.2)then
          dif=4.316
          rad=4.68
        endif
        difnuc(massnr)=dif
        radnuc(massnr)=rad
      enddo

      end

c-----------------------------------------------------------------------
      subroutine xConNuclDens(iii)
c-----------------------------------------------------------------------
c plots distribution of nucleons in nuclei
c  iii = 1 (proj) or 2 (targ)
c-----------------------------------------------------------------------
      include 'epos.inc'
      if(model.ne.1)return
      massnr=1
      if(iii.eq.1)then
        massnr=maproj
      elseif(iii.eq.2)then
        massnr=matarg
      endif
      if(massnr.eq.1)return
      a=1./4.316
      b=4.68
      write(ifhi,'(a)') '!-----------------------------------------'
      write(ifhi,'(a)') '!          nuclear density          '
      write(ifhi,'(a)') '!-----------------------------------------'
      write(ifhi,'(a)')       'openhisto'
      if(massnr.ge.10)write(ifhi,'(a)')'htyp lin xmod lin ymod lin'
      if(massnr.lt.10)write(ifhi,'(a)')'htyp lin xmod lin ymod log'
      write(ifhi,'(a)')       'text 0 0 "title nuclear density"'
      write(ifhi,'(a)')       'text 0.99 -0.15 " r (fm) "'
      write(ifhi,'(a)')       'text 0 0 "yaxis rho(r)"'
      write(ifhi,'(a,2e11.3)')'xrange',0.,mxnucl*drnucl(iii)
      write(ifhi,'(3a)')'yrange',' 0 ',' auto'
      write(ifhi,'(a)')       'array 2'
      do j=1,mxnucl
      r=(j-0.5)*drnucl(iii)
      d=0.5*drnucl(iii)
      write(ifhi,'(2e12.4)')  r,rnucl(j,iii)/nrnucl(iii)/
     *                     (4./3.*pi*((r+d)**3-(r-d)**3))
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lbo '
      write(ifhi,'(a)')       'array 2'
      do j=1,mxnucl
      r=(j-0.5)*drnucl(iii)
      rr=2*r
      rho=1
      if(massnr.eq.2)then
        rho=1.00*((1-exp(-b*a*rr))*exp(-a*rr)/rr)**2
      elseif(massnr.eq.197)then
        rho=0.16/(1+exp((r-6.5)/0.562))
      elseif(massnr.ge.10)then
        rho=0.16/(1+exp((r-radnuc(massnr))/difnuc(massnr)))
      endif
      write(ifhi,'(2e12.4)')  r,rho
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0'
      end

c-----------------------------------------------------------------------
      subroutine xConThick(iii)
c-----------------------------------------------------------------------
      ! plots sigma_pp *T_A (b)  (=average number of collisions)
      ! T_A = thickness function
      !  iii = 1 (proj) or 2 (targ)
      !----------------------------------------------------------------
      include 'epos.inc'
      parameter(iconimax=20,iconkmax=100)
      real thick(2,0:iconimax)
      imax=iconimax
      kmax=iconkmax
      if(model.ne.1)return
      ramx=mxnucl*drnucl(iii)
      do i=0,imax
        x=i/float(imax)*ramx
        sum=0
        rho0=conrho(iii,0.)
        h=ramx/kmax
        do k=1,kmax
          z=k*h
          r=sqrt(x*x+z*z)
          rho2=conrho(iii,r)
          z=(k-0.5)*h
          r=sqrt(x*x+z*z)
          rho1=conrho(iii,r)
          sum=sum+h/6.*(rho0+4*rho1+rho2)
          rho0=rho2
        enddo
        sum=sum*2 ! integral fro -infty to +infty
        thick(1,i)=x
        thick(2,i)=sum
      enddo
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lbo '
      write(ifhi,'(a)')       'txt "xaxis b (fm)" '
      write(ifhi,'(a)')       'txt "yaxis [s]?pp! T?A! (b) " '
      write(ifhi,'(a)')       'array 2'
      do i=0,imax
        write(ifhi,'(2e12.4)') thick(1,i),sigine/10.*thick(2,i)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0'

      end

c-----------------------------------------------------------------------
      function conrho(iii,r)
c-----------------------------------------------------------------------
      ! nuclear density
      !  iii = 1 (proj) or 2 (targ)
      !----------------------------------------------------------------
      include 'epos.inc'
      conrho=1.
      if(model.ne.1)return
      massnr=1
      if(iii.eq.1)then
        massnr=maproj
      elseif(iii.eq.2)then
        massnr=matarg
      endif
      if(massnr.eq.1)return
      a=1./4.316
      b=4.68
      rr=2*r
      rho=1
      if(massnr.eq.2.and.rr.gt.0.)then
        rho=1.00*((1-exp(-b*a*rr))*exp(-a*rr)/rr)**2
      elseif(massnr.eq.197)then
        rho=0.16/(1+exp((r-6.5)/0.562))
      elseif(massnr.ge.10)then
        rho=0.16/(1+exp((r-radnuc(massnr))/difnuc(massnr)))
      endif
      conrho=rho
      end


c----------------------------------------------------------------------
      subroutine geoglauber
c----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incems"
      common/geom/rmproj,rmtarg,bmax,bkmx
      common/nucl3/phi,bimp
      common/ckopjtg/kopj(mamx),kotg(mamx)
      common/cbglaub/bglaub
      logical phi_zero

      phi_zero=.false.


! calculate collision number according to Glauber ------------------

c      rs=r2had(iclpro)+r2had(icltar)+slopom*log(sy)
c      bglaub2=4.*.0389*rs
      xs=sigine   !xsectionpar()   !?????????????????????????????????
      bglaub=max(0.,xs/10./pi)        !10= fm^2 -> mb
c      bglhard=bglaub*min(1.,xfrachard)
      bglaub=sqrt(bglaub)
      nglevt=0
      nglndif=0
      nglhard=0
      bglndif=max(0.,sigcut/10./pi) !10= fm^2 -> mb
c      bglhard=bglndif*xfrachard
      bglndif=sqrt(bglndif)
c      bglhard=sqrt(bglhard)
      do i=1,maproj
        kopj(i)=0
      enddo
      do i=1,matarg
        kotg(i)=0
      enddo
      if(abs(phimin).lt.1e-5.and.abs(phimax).lt.1e-5)phi_zero=.true.
      if(phi_zero)then
        x2av=0
        y2av=0
        xyav=0
        xav=0
        yav=0
        bx=cos(phi)*bimp
        by=sin(phi)*bimp
      else     !not used
        bx=0.
        by=0.
      endif  
      do ko=1,koll
        ip=iproj(ko)
        it=itarg(ko)
        r=bk(ko)
        if(phi_zero)then
          r2=(xproj(ip)+bx-xtarg(it))**2+(yproj(ip)+by-ytarg(it))**2
          if((maproj.gt.1.or.matarg.gt.1).and.
     .       abs(r**2-r2).gt.1e-3)then
            print*,'geoglauber: r**2,r2:',r**2,r2
            stop
          endif
        endif
        if(r.le.bglaub)then
          nglevt=nglevt+1
          kopj(ip)=kopj(ip)+1
          kotg(it)=kotg(it)+1
c pair dependent hard interactions
c          bglhard=max(0.,xfrachard*sngl(sigmak(1,ko))/10./pi) !10= fm^2 -> mb
c          bglhard=sqrt(bglhard)
c          if(r.le.bglhard)nglhard=nglhard+1
        endif
        if(r.le.bglndif)nglndif=nglndif+1
      enddo
      ng1evt=0
      ng11evt=0
      ng12evt=0
      ng2evt=0
      ngspecp=laproj+latarg
      ngspecn=maproj+matarg-ngspecp
      do i=1,maproj
        if(kopj(i).ge.1)then
          ng1evt=ng1evt+1
          ng11evt=ng11evt+1
          if(abs(idptl(i)).eq.1120)ngspecp=ngspecp-1
          if(abs(idptl(i)).eq.1220)ngspecn=ngspecn-1
        endif 
        if(kopj(i).ge.2)ng2evt=ng2evt+1
      enddo   
      do i=1,matarg
        if(kotg(i).ge.1)then
          ng1evt=ng1evt+1
          ng12evt=ng12evt+1
          if(idptl(maproj+i).eq.1120)ngspecp=ngspecp-1
          if(idptl(maproj+i).eq.1220)ngspecn=ngspecn-1
        endif
        if(kotg(i).ge.2)ng2evt=ng2evt+1
      enddo   
      eglevt=0
      fglevt=0
      rglevt=0
      sglevt=0
      if(phi_zero)then
        do i=1,maproj
          if(kopj(i).ge.1)then
            x=xproj(i)+bx/2
            y=yproj(i)+by/2
            xav=xav+x
            yav=yav+y
            x2av=x2av+x*x
            y2av=y2av+y*y
            xyav=xyav+x*y
            !print*,x,y
          endif
        enddo
        do i=1,matarg
          if(kotg(i).ge.1)then
            x=xtarg(i)-bx/2
            y=ytarg(i)-by/2
            xav=xav+x
            yav=yav+y
            x2av=x2av+x*x
            y2av=y2av+y*y
            xyav=xyav+x*y
            !print*,x,y
          endif
        enddo
        if(ng1evt.gt.0)then
        x2av=x2av/float(ng1evt)
        y2av=y2av/float(ng1evt)
        xyav=xyav/float(ng1evt)
        xav= xav /float(ng1evt)
        yav= yav /float(ng1evt)
        eglevt=(y2av-x2av)/(y2av+x2av)
        fglevt=sqrt( (y2av-x2av)**2+4*(xyav-xav*yav)**2 )  /(y2av+x2av)
        rglevt=ng2evt/float(ng1evt)
        sglevt=pi*sqrt(x2av)*sqrt(y2av)
        !~~~~~~~~~~~~~~~~~~~~
        xx=x2av ! <x**2>
        yy=y2av ! <y**2>
        xy=xyav ! <x*y>
        dta=0.5*abs(xx-yy)
        ev1=(xx+yy)/2+sqrt(dta**2+xy**2)
        ev2=(xx+yy)/2-sqrt(dta**2+xy**2)
        yy=ev1
        xx=ev2
        fglevt=(yy-xx)/(yy+xx)
        !~~~~~~~~~~~~~~~~~~~~
        endif
        !if(ng1evt.gt.0)print*,bimp,eglevt,bx,by,phi
      endif  
  
      end

c-----------------------------------------------------------------------
      function percentile()
c-----------------------------------------------------------------------
      include "epos.inc"  
      common/ccc20/icc20 
      common/nucl3/phi,bimp
      real wi(2)
      if(izmode.ne.1)stop'\n\n  ERROR 24022013c \n\n'
      if(icc20.ne.1)stop'\n\n  ERROR 24022013b \n\n'
      do i=1,20
        ii=i
        if(bimp.ge.zclass(1,i).and.bimp.le.zclass(2,i))goto 1
      enddo
  1   cmin=(ii-1)*0.05
      cmax= ii*0.05
      bmin=zclass(1,ii)
      bmax=zclass(2,ii)
      frac=(bimp-bmin)/(bmax-bmin)
      wi(1)=1-frac
      wi(2)=frac
      percentile   = wi(1)*cmin +wi(2)*cmax
      end
      





