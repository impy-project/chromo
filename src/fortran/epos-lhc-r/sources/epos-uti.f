C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
      subroutine psmirror(ep3)
c-----------------------------------------------------------------------
      common/cmirror/mirror
      dimension ep3(4)
      if(mirror.eq.1)ep3(2)=-ep3(2)
      end
c-----------------------------------------------------------------------
      subroutine psmirror1(ep)
c-----------------------------------------------------------------------
      common/cmirror/mirror
      if(mirror.eq.1)ep=-ep
      end
c-----------------------------------------------------------------------
      subroutine setmirror(ivalue)
c-----------------------------------------------------------------------
      common/cmirror/mirror
      mirror=ivalue
      end

c-----------------------------------------------------------------------
      subroutine swap(x,y)
c-----------------------------------------------------------------------
      xold=x
      x=y
      y=xold
      end
c-----------------------------------------------------------------------
      subroutine iswap(ix,iy)
c-----------------------------------------------------------------------
      ixold=ix
      ix=iy
      iy=ixold
      end

c-----------------------------------------------------------------------
      logical function leq(x,y)
c-----------------------------------------------------------------------
      if(nint(x*10).eq.nint(y*10))then
        leq=.true.
      else
        leq=.false.
      endif
      end 
c-----------------------------------------------------------------------
      subroutine getSystemABE(ia,ib,e)
c-----------------------------------------------------------------------
      include "epos.inc"
      ia=maproj
      ib=matarg 
      e=engy
      end
c-----------------------------------------------------------------------
      subroutine getSystemType(isys)
c-----------------------------------------------------------------------
      include "epos.inc"
      isys=0
      if(iappl.ne.1)return ! not hadron
      if(iabs(maproj)*iabs(matarg).eq.1)then !pp
         isys=1
      elseif(iabs(maproj).le.3.or.iabs(matarg).le.3)then!aA (a=p,d,He)
         isys=2
      elseif(max(iabs(maproj),iabs(matarg)).lt.170)then!mid AA
         isys=3
      else !big AA
         isys=4
      endif
      end

c-----------------------------------------------------------------------
      integer function igetNpom()
c-----------------------------------------------------------------------
      include "epos.inc"
      n=0
      do i=1,nptl
        call getistptl(i,istx)
        if(istx.eq.30.or.istx.eq.31)n=n+1
      enddo
      igetNpom=n
      end

c-----------------------------------------------------------------------
      subroutine getNhpom(i) !hard Poms, known after call ProPoTy in emsaa
c-----------------------------------------------------------------------
      include "epos.incems"
      i=npr(3,1) !3 means "hard", 1 means "collision 1"
      end

c----------------------------------------------------------------------
      function polar(x,y)
c----------------------------------------------------------------------
      pi=3.1415927
      if(abs(x).gt.1.e-6)then
        phi=atan(y/x)
        if(x.lt.0.)phi=pi+phi
        if(phi.lt.0)phi=2*pi+phi
      else
        phi=0.5*pi
        if(y.lt.0)phi=phi+pi
      endif
      polar=phi
      end

c-----------------------------------------------------------------------
      subroutine setDefine(line,i,j,k)
c-----------------------------------------------------------------------
      include "epos.inc"
      character line*1000
      ne=0
      call utword(line,i,j,ne)
      if(line(i:j).eq.'bim3')stop'****** bim3->bim03 ****** '
      if(line(i:j).eq.'bim5')stop'****** bim5->bim05 ****** '
      if(line(i:j).eq.'bim6')stop'****** bim6->bim06 ****** '
      if(ndefine+1.gt.mxdefine)stop'too many `define` statements.      '
      l1=j-i+1
      if(l1.gt.99)stop'`define` argument 1 too long.            '
      w1define(ndefine+1)(1:l1)=line(i:j)
      w1define(ndefine+1)(l1+1:l1+1)=' '
      call utword(line,i,j,ne)
      l2=j-i+1
      if(l2.gt.99)stop'`define` argument 2 too long.            '
      w2define(ndefine+1)(1:l2)=line(i:j)
      w2define(ndefine+1)(l2+1:l2+1)=' '
      ndefine=ndefine+1
      l1define(ndefine)=l1
      l2define(ndefine)=l2
      if(k.eq.2)then
        call utword(line,i,j,ne)
        call utword(line,i,j,ne)
      endif
      end
      
c-----------------------------------------------------------------------
      integer function igetYieldCharged()
c-----------------------------------------------------------------------
      include "epos.inc"
      n=0
      do i=1,nptl
        call getistptl(i,istx)
        if(istx.eq.0)then
          call getidptl(i,id)
          call idchrg( 999 ,id,ch)
          call idflav(id,i1,i2,i3,jdu,idu)
          call getpptl(i,p1,p2,p3,p4,p5)
          if(abs(ch).gt.0.1.and.i2.ne.0.and.i3.ne.0)then
            pt2=p1**2+p2**2
            pp=sqrt(pt2+p3**2)
            if(pt2.gt.0)then  
              pt=sqrt(pt2)
              eta=sign(1.,p3)*alog((pp+abs(p3))/pt)
            else
              eta=1e10                  !
            endif
            if(abs(eta).le.1.0)then
             n=n+1
            endif
          endif
        endif
      enddo
      igetYieldCharged=n
      end

c-----------------------------------------------------------------------
      integer function igetYieldCharm12()
c-----------------------------------------------------------------------
      include "epos.inc"
      n=0
      do i=1,nptl
        call getistptl(i,istx)
        call getidptl(i,id)
        ida=abs(id)
        !if(istx.eq.21.and.ida.eq.4)print*,'CHARMQUA',id
        call getpptl(i,p1,p2,p3,p4,p5)
        if(ida.eq.140.or.ida.eq.240.or.ida.eq.241)then !used by ALICE
          amt=p5**2+p1**2+p2**2
          if(amt.gt.0..and.p4+abs(p3).gt.0.)then  
            amt=sqrt(amt)
            rap=sign(1.,p3)*log((p4+abs(p3))/amt) 
            pt=sqrt(p1**2+p2**2)
            !print*,'CHARMMES',id,rap,pt
          else
            rap=1e10                  !
          endif
          if(abs(rap).le.0.5.and.pt.gt.1.and.pt.lt.2)n=n+1
        endif
      enddo
      igetYieldCharm12=n
      end

c-----------------------------------------------------------------------
      function s1234(x)
c-----------------------------------------------------------------------
      double precision x(4) 
      s1234= (x(1)-x(2)) * (x(1)+x(2)) - x(3)**2 - x(4)**2
      end

c-----------------------------------------------------------------------
      subroutine getEposPath(path,length)
c-----------------------------------------------------------------------
      include "epos.inc"
      character*500 path 
      path=fnnx
      length=nfnnx
      end

c-----------------------------------------------------------------------
      subroutine getMonitorFileIndex(ifmtx)
c-----------------------------------------------------------------------
      include "epos.inc"
      ifmtx=ifmt
      end

c-----------------------------------------------------------------------
      subroutine getCheckFileIndex(ifchx)
c-----------------------------------------------------------------------
      include "epos.inc"
      ifchx=ifch
      end

c-----------------------------------------------------------------------
      integer function igetIsh()
c-----------------------------------------------------------------------
      include "epos.inc"
      igetIsh=ish
      end
      subroutine setIsh(ishx)
      include "epos.inc"
      ish=ishx
      end

c-----------------------------------------------------------------------
      subroutine setIshy(ishyx)
c-----------------------------------------------------------------------
      common/cishy/ishy
      ishy=ishyx
      end
      subroutine getIshy(ishyx)
      common/cishy/ishy
      ishyx=ishy
      end

c-----------------------------------------------------------------------
      subroutine storeValues(a,b,c,d,e)
c-----------------------------------------------------------------------
      common/cParams/ax,bx,cx,dx,ex
      ax=a
      bx=b
      cx=c
      dx=d
      ex=e
      end
      subroutine restoreValues(a,b,c,d,e)
      common/cParams/ax,bx,cx,dx,ex
      a=ax
      b=bx
      c=cx
      d=dx
      e=ex
      end

c-----------------------------------------------------------------------
      subroutine setVparam(i,val)
c-----------------------------------------------------------------------
      include "epos.inc"
      vparam(i)=val
      end
      subroutine getVparam(i,val)
      include "epos.inc"
      val=vparam(i)
      end
      subroutine incVparam(i,val)
      include "epos.inc"
      vparam(i)=vparam(i)+val
      end

c-----------------------------------------------------------------------
      subroutine setLadderindex(ladderindexx)
c-----------------------------------------------------------------------
      common/cladderindex/ladderindex
      ladderindex=ladderindexx
      end
      subroutine getLadderindex(ladderindexx)
      common/cladderindex/ladderindex
      ladderindexx=ladderindex
      end

c-----------------------------------------------------------------------
      subroutine setKomg(ivalue)
c-----------------------------------------------------------------------
      integer komg
      common/ckomg/komg
      komg=ivalue
      end
      subroutine getKomg(ivalue)
      integer komg
      common/ckomg/komg
      ivalue=komg
      end

c-----------------------------------------------------------------------
      subroutine updateXor(n,nn)
c-----------------------------------------------------------------------
      include "epos.inc"
      real xori(4)
      call getiorptl(n,ior1)  !father
      ist1=0
      if(ior1.gt.0)call getistptl(ior1,ist1) !status
      call getityptl(nn,ity)
      if(ity.eq.4.or.ity.eq.5)ist1=41    !for spectators
      if (ist1.eq.29)then        !string
        taux=taustr
      elseif(ist1.eq.41)then     !remnant
        taux=taurem
      elseif(ist1.eq.10.or.ist1.eq.11.or.ist1.eq.12)then     !core
        taux=tauzer2
      else
        goto 100
      endif
      xori(1)=xorptl(1,ior1)
      xori(2)=xorptl(2,ior1)
      xori(3)=xorptl(3,ior1)
      xori(4)=xorptl(4,ior1)
      ptm=1.
      if(pptl(5,nn).gt.0.)ptm=pptl(5,nn)
c      ptm=sqrt(pptl(5,nn)**2+pptl(1,nn)**2+pptl(2,nn)**2)
      if(taux.gt.0.)then
        r=rangen()
        tauran=-abs(taux)*log(r)*pptl(4,nn)/ptm
      else
        tauran=abs(taux)*pptl(4,nn)/ptm
      endif
c      write(ifch,*)'update',nn,ptm,pptl(4,nn),xori(4),tauran
      xorptl(1,nn)=xori(1)+pptl(1,nn)/pptl(4,nn)*tauran
      xorptl(2,nn)=xori(2)+pptl(2,nn)/pptl(4,nn)*tauran
      xorptl(3,nn)=xori(3)+pptl(3,nn)/pptl(4,nn)*tauran
      xorptl(4,nn)=xori(4)+pptl(4,nn)/pptl(4,nn)*tauran
      tivptl(1,nn)=xorptl(4,nn)
  100 continue
      call idtau(idptl(nn),pptl(4,nn),pptl(5,nn),taugm)
      tivptl(2,nn)=tivptl(1,nn)+taugm*(-log(rangen()))
      end

!---------------------------------------------------------------------------------
!> @brief
!> get ity and jor of 25-parton associated to given hard parton or hadron
!>  in case of soft origin:               return ity = jor = -992
!>  in case of projectile remnant origin: return ity = jor = -994
!>  in case of target remnant origin:     return ity = jor = -995
!> @author Klaus WERNER 
c---------------------------------------------------------------------------------
      subroutine getItyJor25(j,ityxx,jorxx) !ity,jor of 25-parton
c---------------------------------------------------------------------------------
      integer, intent(in) :: j       !< particle index (of parton or hadron)
      integer, intent(out) :: ityxx  !< ity value of 25-parton
      integer, intent(out) :: jorxx  !< jor value of 25-parton
      include "epos.inc"
      data ncount/0/
      save ncount
      if(j.le.0.or.model.ne.1)return
      idepos=ideposf( 3 ,j)
      ityxx=-999
      jorxx=-999
      iok=0
      jj=0
      if((abs(idepos).ge.1.and.abs(idepos).le.6         !quarks
     .   .or.abs(idepos).eq.9                                !gluons
     .   .or.(abs(idepos).lt.9999.and.mod(abs(idepos),100).eq.0) !diquarks
     .  ).and.j.gt.1)then    
        call getityptl(j,ityjj)
        if(ityjj/10.eq.3)then !hard
          jj=j
          iok=1
        elseif(ityjj/10.eq.2)then !soft
          if(istptl(j).eq.21)then
           ityxx=-992
           jorxx=-992
          endif
        elseif(ityjj/10.eq.4)then !proj remnant
          ityxx=-994
          jorxx=-994
        elseif(ityjj/10.eq.5)then !targ remnant
          ityxx=-995
          jorxx=-995
        endif
      elseif(ihacas.eq.0)then !otherwise makes no sense
        call index2521(j,i25first,i21last)
        if(i25first.ne.0)then                 !hadrons with ity/10 = 3
          jj=i25first
          iok=1
        endif
      endif
      if(iok.eq.1.and.jj.le.nptl)then
       ityxx=-1
       jorxx=-1
       if(istptl(jj).ne.21.and.istptl(jj).ne.25)then
         ncount=ncount+1 
         if(ncount.eq.1.and.ish.ge.1)then
           write(ifmt,*)'WARNING in getItyJor25() : parton found'
           write(ifch,*)'WARNING in getItyJor25() : parton found'
           do jik=1,jj
           write(ifch,*)'PTL',iorptl(jik),jorptl(jik),jik
     .      ,ifrptl(1,jik),ifrptl(2,jik),idptl(jik),istptl(jik)
     .      ,ityptl(jik)
           enddo 
           !stop'ERROR 210806'
         endif
       else
         return
       endif
       jchk=jj
       !print*,'CHK    ',jchk,idptl(jchk),istptl(jchk),istptl(jchk-1)
       do while(istptl(jchk).ne.25
     .  .and.istptl(jchk-1).ge.21.and.istptl(jchk-1).le.28)
         jchk=jchk-1
       !print*,'CHK    ',jchk,idptl(jchk),istptl(jchk),istptl(jchk-1)
       enddo
       if(istptl(jchk).eq.25)then
         jorxx=jorptl(jchk)
         ityxx=ityptl(jchk)
         !print*,'CHK===>',jchk,idptl(jchk),istptl(jchk),jorxx
       endif
       !print*,'CHK-----------------------------------------',jorxx
      endif
      end

c-----------------------------------------------------------------------
      function utgenBW(a,g)
c-----------------------------------------------------------------------
      ! generates mass x according to 
      ! Relativistic Breit–Wigner distribution
      ! f(x2) ~ 1 / ( (x2-a**2)**2 + a**2*g**2 )
      !  x2 = mass squared x**2
      !  a = average mass
      !  g = gamma = width
      ! variable transform: 
      !    x2=a*g*tan(y)+a**2
      ! |dx2/dy|=m*g/cos(y)**2
      ! gives prob(y)=flat
      !-------------------------------------------
      ! x2 between 2*m_pi and \infty
      ami=2*0.14 
      ymin=-atan((a**2-ami**2)/g/a)
      ymax=3.14159/2.
      y=ymin+rangen()*(ymax-ymin)
      x2=a*g*tan(y)+a**2
      x=sqrt(x2)
      utgenBW=x
      end

c-----------------------------------------------------------------------
      subroutine utphiBW(i)
c-----------------------------------------------------------------------
      include "epos.inc"
      common/creswi/reswi 
      call getidptl(i,id)
      call getityptl(i,ity)
      idepos=id
      if(ity.eq.61)idepos=idtrafo('pdg','nxs',id)
      if(idepos.ne.331)stop'####### ERROR 20112015 #######'
      !old:
      call getpptl(i,p1,p2,p3,p4,p5)
      amt=sqrt(p1**2+p2**2+p5**2)
      rap=sign(1.,p3)*log((p4+abs(p3))/amt)
      !new:
      call idtau(idepos,1.,1.,dummy)
      g=reswi
      call idmass(idepos,a)
c      print*,g,a,engy
      ntry=0
 10   ntry=ntry+1
      if(ntry.gt.999)call utstop('problem in utphiBW !&')
      am=utgenBW(a,g)
c      print*,'utphiBW',idepos,p5,'--->',a,g,'--->',am
      if(am.ge.engy.or.abs(a-am).gt.10.*g)goto 10  !limit mass to available energy and may 10.*width
      p5=am
      amt=sqrt(p1**2+p2**2+p5**2)
      p3=amt*sinh(rap)
      p4=amt*cosh(rap)
      call setpptl(i,p1,p2,p3,p4,p5)
      end

c-----------------------------------------------------------------------
      subroutine testBW
c-----------------------------------------------------------------------
      parameter (jmax=20,n=10000)
      real w(jmax)
      do j=1,jmax
        w(j)=0
      enddo
      amin=1.
      amax=1.04
      dam=(amax-amin)/jmax
      call idmass(331,amt)
      do i=1,n
       rap=rangen()*5 - 2.5
       p3=amt*sinh(rap)
       p4=amt*cosh(rap)
       p1=rangen()
       p2=rangen()
       p5=1
       call setidptl(i,331)
       call setpptl(i,p1,p2,p3,p4,p5)
      enddo
      do i=1,n
       call utphiBW(i)
       call getpptl(i,p1,p2,p3,p4,p5)
       j= (p5-amin)/dam+1
       if(j.ge.1.and.j.le.jmax)w(j)=w(j)+1
      enddo
      do j=1,jmax
        am=amin+(j-0.5)*dam
        a=1.0195
        g=0.0043 
        f= 1 / ( (am**2-a**2)**2 + a**2*g**2 )
        print*, am , w(j)/n/dam , f*0.0027
      enddo
      stop'Normal stop in testBW'
      end 

c-----------------------------------------------------------------------
      logical function nequal(a,b)
c-----------------------------------------------------------------------
      nequal=.false.
      if(abs(a-b).gt.1e-4*max(a,b))nequal=.true.
      end

c-----------------------------------------------------------------------
      subroutine memo(i,text)
c-----------------------------------------------------------------------
#ifndef CHROMO
      include "epos.inc"
      character text*(*)
      n=index(text,';')
      if(n.gt.1)write(ifmt,'(a,$)')text(1:n)
      call checkmemory(memory)
      if(i.eq.0)write(ifmt,'(a,i8)')'  memory',memory
      if(i.eq.1)write(ifmt,'(a,i8,$)')'  memory',memory
      if(i.eq.2)write(ifmt,'(a,i8)')'  ==>',memory
#endif
      end

c-----------------------------------------------------------------------
      function AngleEllipsoid(xx,yy,xy)
c-----------------------------------------------------------------------
      common/cecc/ecct,rrrt
      pi=3.1415927
      ranphi=0
      xpxp=xx
      ypyp=yy
      dta=0.5*(yy-xx)
      eba=0.5*(xx+yy)
      ww=-xy
      if(dta.ne.0.)then
       !------------------------------
       !make rotation of axes (angle +phi)
       ! xp=  x*cosphi + y*sinphi
       ! yp= -x*sinphi + y*cosphi
       !chose phi such that xpyp=0
       !-------------------------------
       ranphi=0.5*atan(ww/dta)
       coco=cos(ranphi)**2
       sisi=sin(ranphi)**2
       cosi=cos(ranphi)*sin(ranphi)
       xpxp=xx*coco+yy*sisi+2*xy*cosi
       ypyp=xx*sisi+yy*coco-2*xy*cosi
       if(ypyp.lt.xpxp)then
        zzzz=xpxp
        xpxp=ypyp
        ypyp=zzzz
        ranphi=ranphi+pi/2
       endif
      endif
      AngleEllipsoid=ranphi
      ecct=(ypyp-xpxp)/(ypyp+xpxp)
      rrrt=sqrt((ypyp+xpxp)/2.)
      end
      
c-----------------------------------------------------------------------
      subroutine utcopy
c-----------------------------------------------------------------------
      include "epos.inc"
      common/cutcopy/nptl1xx,nptl2xx
      if(irescl.eq.3)then
      if(ispherio.eq.1.or.ikolmn.gt.10)
     .write(ifmt,'(a)')'copy for rescale'
      nptlx=nptl
      do iloo=1,nptlx
        ip=iloo   
c        if(iloo.gt.mxptl)then
c          call restorecccptl(iloo,mxptl+2)
c          ip=mxptl+2
c        endif
        if(ip.le.nptlpt.or.istptl(ip).eq.0)then
          ist=-1
          if(istptl(ip).eq.0)ist=-2
          nptl=nptl+1
          call checkcccptl(nptl)
c          if(nptl.le.mxptl)then
            call utrepl(nptl,ip)
            istptl(nptl)=ist
c          else
c            call utrepl(mxptl+1,ip)
c            istptl(mxptl+1)=ist
c            call dumpcccptl(mxptl+1,nptl)
c          endif
        endif
      enddo  
      nptl1xx=nptlx+1
      nptl2xx=nptlx+nptlpt
      endif
      end
      
c-----------------------------------------------------------------------
      subroutine rescaleRap2(scal,asym , nptl2,ist_0,ii , p3,p4)
c-----------------------------------------------------------------------
      include "epos.inc"
      real p(5)
      sum3=0.
      sum4=0.
      do  j=nptl2+1,nptl 
        call getistptl(j,istx)
        call getityptl(j,ityx)
        if(mod(istx,10).eq.ist_0.and.ityx.ne.4.and.ityx.ne.5)then
          call getpptl(j,p(1),p(2),p(3),p(4),p(5))
          if(p(4).gt.1d-4)then
            !p(3)=scal*p(3)
            !p(4)=sqrt(p(1)**2+p(2)**2+p(3)**2+p(5)**2)
            amt=sqrt(p(1)**2+p(2)**2+p(5)**2)
            rap=sign(1.,p(3))*alog((p(4)+abs(p(3)))/amt) 
            if(rap.gt.0.)then
              rap=rap*scal*(1+asym)
            else
              rap=rap*scal*(1-asym)
            endif 
            p(3)=amt*sinh(rap)
            p(4)=amt*cosh(rap)
            if(ii.eq.1)call setpptl(j,p(1),p(2),p(3),p(4),p(5))
            sum3=sum3+p(3)
            sum4=sum4+p(4)
          endif
        endif
      enddo
      p3=sum3
      p4=sum4
      end

c-----------------------------------------------------------------------
      subroutine rescaleRap1(scal , nptl2,ist_0 , asym,p4,iret) 
c-----------------------------------------------------------------------
      include "epos.inc"
      iret=0
      ii=0
      fab=1
      scalmax=0
      do while(scalmax.le.0.85.and.fab.gt.0)
        scalmax=scalmax+0.10
        a=-1 
        asym=a*scalmax
        call rescaleRap2(scal,asym , nptl2,ist_0,ii , p3,p4)
        fa=p3
        b=1
        asym=b*scalmax
        call rescaleRap2(scal,asym , nptl2,ist_0,ii , p3,p4)
        fb=p3 
        fab=fa*fb
      enddo
      if(fb*fa.ge.0)then        !should be < 0
        iF(ish.ge.1)then
        write(ifmt,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
        write(ifmt,*)'ERROR 08042021a in rescaleRap1'
        write(ifmt,*)' fa fb = ',fa,fb,'  (should have diff sign)' 
        write(ifmt,*)'scalmax probably too small. CHECK:'
        do kk=-9,9 
         asym=kk/5.*scalmax
         call rescaleRap2(scal,asym , nptl2,ist_0,ii , p3,p4)
         write(ifmt,*)'CHECK rescaleRap1',scal,asym,p3
        enddo
        write(ifmt,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
        endif
        iret=1
        return
c       stop
      endif
      fc=1e20
      ipass=0
      do while(ipass.lt.100.and.abs(fc).gt.0.00001*p4)
        ipass=ipass+1
        c=(a+b)/2 
        asym=c*scalmax
        call rescaleRap2(scal,asym , nptl2,ist_0,ii , p3,p4)
        fc=p3 
        if(ish.ge.4)write(ifch,'(a,i5,2f12.2,4x,2f9.4)')
     .  'rescaleRap1: iter p3 p4 =',ipass,p3,p4,scal,asym
        if(fa*fc.le.0.)then
          b=c
          fb=fc
        else
          a=c
          fa=fc
        endif
      enddo
      if(abs(fc).gt.0.0001*p4.and.ish.ge.1)then
          write(ifmt,*) 
     .    'WARNING rescaleRap1: fc = ',fc,' -> too big!'
      endif
      end

c-----------------------------------------------------------------------
      subroutine rescaleRap(nptl2,ist_0,esoll,iret) 
c-----------------------------------------------------------------------
      include "epos.inc"
      double precision esoll,fab
      iret=0
      fab=1
      fa=1.
      fb=1.
      scalmax=0
      do while(scalmax.le.0.70.and.fab.gt.0)
        scalmax=scalmax+0.10
        a=-1 
        scal=1+a*scalmax
        call rescaleRap1(scal , nptl2,ist_0 , asym,p4,iret)
        if(iret.ne.0)return
        fa=p4-esoll
        b=1
        scal=1+b*scalmax
        call rescaleRap1(scal , nptl2,ist_0 , asym,p4,iret) 
        if(iret.ne.0)return
        fb=p4-esoll
        fab=fa*fb
        call checkFaFb(a,b,fa,fb,scalmax,nptl2,ist_0  
     .    , c,fc,scal,asym,iretab,iret)
        if(iret.ne.0)return
        if(iretab.eq.1)goto 888
      enddo
      do while(b.lt.3.and.fab.gt.0)
        b=b+0.1
        scal=1+b*scalmax
        call rescaleRap1(scal , nptl2,ist_0 , asym,p4,iret) 
        if(iret.ne.0)return
        fb=p4-esoll
        fab=fa*fb
        call checkFaFb(a,b,fa,fb,scalmax,nptl2,ist_0  
     .    , c,fc,scal,asym,iretab,iret)
        if(iret.ne.0)return
        if(iretab.eq.1)goto 888
      enddo
      kk=0 
      do while(kk.lt.20.and.fab.gt.0)
        kk=kk+1
        scal=1+a*scalmax
        scalnew=scal/2.
        a=(scalnew-1)/scalmax
        scal=1+a*scalmax
        if(abs(scal-scalnew).gt.1e-5*scal
     .  .and.abs(scal-scalnew).gt.1e-6)then
          write(ifmt,*)'++++++++++++++++++++++++'
          write(ifmt,*)'ERROR 11062022 in rescaleRap'
          write(ifmt,*)'scal,scalnew,scalmax',scal,scalnew,scalmax
          write(ifmt,*)'++++++++++++++++++++++++'
          stop
        endif
        call rescaleRap1(scal , nptl2,ist_0 , asym,p4,iret) 
        if(iret.ne.0)return
        fa=p4-esoll
        fab=fa*fb
        call checkFaFb(a,b,fa,fb,scalmax,nptl2,ist_0  
     .    , c,fc,scal,asym,iretab,iret)
        if(iret.ne.0)return
        if(iretab.eq.1)goto 888
      enddo
      if(fa*fb.ge.0)then !should be < 0
        if(ish.ge.1)write(ifmt,*)
     .  'WARNING rescaleRap No scale Redo event' 
        iret=1
        return
        !write(ifmt,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
        !write(ifmt,*)'ERROR 08042021b in rescaleRap'
        !write(ifmt,*)' fa fb = ',fa,fb,'  (should have diff sign)' 
        !write(ifmt,*)'scalmax probably too small. CHECK:'
        !do kk=-9,9 
        ! scal=1+kk/5.*scalmax
        ! if(scal.gt.0.)then
        !   call rescaleRap1(scal, nptl2,ist_0 , asym,p4,iret)
        !   write(ifmt,*)'CHECK rescaleRap',scal,asym,p4-esoll
        ! endif
        !enddo
        !write(ifmt,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
        !stop
      endif
      fc=1e20
      ipass=0
      do while(ipass.lt.20.and.abs(fc).gt.0.0001*esoll)
        ipass=ipass+1
        c=(a+b)/2 
        scal=1+c*scalmax
        call rescaleRap1(scal , nptl2,ist_0 , asym,p4,iret) 
        if(iret.ne.0)return
        fc=p4-esoll
        if(ish.ge.4)write(ifch,'(a,i3,2f12.2,4x,f9.4,f9.5)')
     .  'rescaleRap: it p4 Esoll =',ipass,p4,esoll,scal,asym
        if(fa*fc.le.0.)then
          b=c
          fb=fc
        else
          a=c
          fa=fc
        endif
      enddo
  888 if(abs(fc).gt.0.001*esoll.and.ish.ge.1)then
          write(ifmt,*) 
     .    'WARNING rescaleRap: fc = ',fc,' -> too big!'
      endif
      call rescaleRap2(scal,asym , nptl2,ist_0, 1  , p3,p4)
      if(ish.ge.1)then
      if(abs(scal-1).gt.1e-3.and.abs((p4-esoll)/esoll).gt.1e-3)
     . write(ifmt,'(2a,i3,f8.2,f10.4,f10.4,f10.5)')'WARNING rescaleRap:'
     .     ,' it p3 dE/E scales =',ipass,p3,(p4-esoll)/esoll,scal-1,asym
      endif
      end
c---------------------------
      subroutine checkFaFb(a,b,fa,fb,scalmax,nptl2,ist_0
     .    ,c,fc,scal,asym,iretab,iret)
      iretab=0
      if(fa.eq.0..or.fb.eq.0.)then
        if(fa.eq.0.)then
          iretab=1
          c=a
        elseif(fb.eq.0.)then
          iretab=1
          c=b
        endif
        fc=0 
        scal=1+c*scalmax
        call rescaleRap1(scal , nptl2,ist_0 , asym,p4,iret) 
      endif
      end

c-----------------------------------------------------------------------
      subroutine utrescxx(iret,icopy)
c-----------------------------------------------------------------------
c New procedure based on rapidity rescaling, should always work 
c-----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incems"
      common/cutcopy/nptl1xx,nptl2xx
      double precision p1,p2,p3,esoll,utdble,xx1,xx2,xx3,xx4
      real p(5)
      dimension p1(5),p2(5) ,p3(5) 
      call utpri('utrexx',ish,ishini,4)
      
      iret=0
      if(irescl.eq.0)return
      nptlpt=max(2,iabs(maproj)+iabs(matarg))

      imin=nptlpt+1
      if(iappl.eq.1)then
        imin=1
      endif

      nptl1=1
      nptl2=nptlpt
      ist_0=0
      ist_1=1
      
      if(icopy.eq.1)then
        call utcopy
        imin=nptl1xx
        nptl1=nptl1xx
        nptl2=nptl2xx
        ist_0=-2
        ist_1=-1
      endif

      if(nptl.le.nptl2)then
        write(ifmt,'(80a)')('-',n=1,80)
        do n=1,nptl
          write(ifmt,'(7i8)')iorptl(n),jorptl(n),n,(ifrptl(i,n),i=1,2)
     .    ,idptl(n),istptl(n)   
        enddo
        write(ifmt,'(a)')'  ERROR 08042021c'
        write(ifmt,'(80a)')('-',n=1,80)
        call utstop("utrescxx (1) &")
      endif

c     compute p1
c     ----------
      
      esoll=0.d0
      p1(1)=0.d0
      p1(2)=0.d0
      p1(3)=0.d0
      p1(4)=0.d0  
      do i=nptl1,nptl2
        call getistptl(i,istx)
        call getityptl(i,ityx)
        if(istx.eq.ist_1.and.ityx.eq.0)then !participant
          call getpptl(i,p(1),p(2),p(3),p(4),p(5))
          do j=1,4
            p1(j)=p1(j)+utdble(p(j))
          enddo
        endif  
      enddo
      p1(5)=dsqrt((p1(4)+p1(3))*(p1(4)-p1(3))-p1(2)*p1(2)-p1(1)*p1(1))
      esoll=p1(5)
      if(ish.ge.4) write (ifch,'(a,5g13.6)') 'boost-vector p1',p1

c     trafo + compute p2
c     ------------------
 
      p2(1)=0.d0
      p2(2)=0.d0
      p2(3)=0.d0
      p2(4)=0.d0
      do i=nptl2+1,nptl
        call getistptl(i,istx)
        call getityptl(i,ityx)
        if(mod(istx,10).eq.ist_0.and.ityx.ne.4.and.ityx.ne.5)then
          call getpptl(i,p(1),p(2),p(3),p(4),p(5))
          !check that particles are on-shell and force it if not
          am=(p(4)+p(3))*(p(4)-p(3))-p(1)*p(1)-p(2)*p(2)-p(5)*p(5)
          if(abs(am)/max(1e-4,p(4)).gt.5.e-4.or.p(4).le.1e-4)then      !particle off-shell, redefine energy
            if(p(5).le.0.)then
              call getidptl(i,idx)
              call idmass(idx,am)
              p(5)=am
            endif
            p(4)=sqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3)+p(5)*p(5))
          endif
          if(p(4).gt.1e-4)then
            xx1=utdble(p(1))
            xx2=utdble(p(2))
            xx3=utdble(p(3))
            xx4=sqrt(xx1**2+xx2**2+xx3**2+utdble(p(5))**2)
            call utlob2(1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $           ,xx1,xx2,xx3,xx4,71)
            p(4)=sngl(sqrt(xx1**2+xx2**2+xx3**2+utdble(p(5))**2))
            if(abs(p(4)-sngl(xx4)).gt.1e-4*p(4))
     $                           call utstop('Pb in utrescxx !&')
            p(3)=sngl(xx3)
            p(2)=sngl(xx2)
            p(1)=sngl(xx1)
            call setpptl(i,p(1),p(2),p(3),p(4),p(5))
            do j=1,4
              p2(j)=p2(j)+utdble(p(j))
            enddo
          endif
        endif
      enddo
      p2(5)=dsqrt((p2(4)+p2(3))*(p2(4)-p2(3))-p2(2)*p2(2)-p2(1)*p2(1))
      if(ish.ge.4) write (ifch,'(a,5g13.6)') 'boost-vector p2',p2
      
c     trafo + compute p3
c     ------------------
 
      p3(1)=0.d0
      p3(2)=0.d0
      p3(3)=0.d0
      p3(4)=0.d0
      do i=nptl2+1,nptl
        call getistptl(i,istx)
        call getityptl(i,ityx)
        if(mod(istx,10).eq.ist_0.and.ityx.ne.4.and.ityx.ne.5)then
          call getpptl(i,p(1),p(2),p(3),p(4),p(5))
          if(p(4).gt.1e-4)then
            call utlob4(1,p2(1),p2(2),p2(3),p2(4),p2(5)
     $    ,p(1),p(2),p(3),p(4))
            call setpptl(i,p(1),p(2),p(3),p(4),p(5))
            do j=1,4
              p3(j)=p3(j)+utdble(p(j))
            enddo
          endif
       endif
      enddo
      p3(5)=dsqrt((p3(4)+p3(3))*(p3(4)-p3(3))-p3(2)*p3(2)-p3(1)*p3(1))
      if(ish.ge.4) write (ifch,'(a,5g13.6)') 'boost-vector p3',p3

c     rescale momenta in rest frame
c     -----------------------------

      call rescaleRap(nptl2,ist_0,esoll,iret)
      if(iret.eq.1)return
 
c     trafo
c     -----
      do i=nptl2+1,nptl
        call getistptl(i,istx)
        call getpptl(i,p(1),p(2),p(3),p(4),p(5))
        call getityptl(i,ityx)
        if(mod(istx,10).eq.ist_0.and.p(4).gt.1e-4
     $     .and.ityx.ne.4.and.ityx.ne.5)then
          call utlob4(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $         ,p(1),p(2),p(3),p(4))
          call setpptl(i,p(1),p(2),p(3),p(4),p(5))
          if(icopy.eq.0)call updateXor(i,i)  
        endif
      enddo
      
      call utprix('utrexx',ish,ishini,4)

      end

c-----------------------------------------------------------------------
      subroutine utrfrag(ir,nmin,iret)
c-----------------------------------------------------------------------
c     rescaling to conserve energy and momentum for particles from nmin to nptl
c     ir=1 for projectile
c     ir=2 for target
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      double precision p1(5),esoll,ppp,psoll,pt1soll,pt2soll
      integer nmin,iret
      logical lim
      call utpri('utrfrag',ish,ishini,4)

      errlim=0.005

      iret=0
      if(nptl.le.nptlpt) goto 9999

      esoll=0.d0
      psoll=0.d0
      p1(1)=0.d0
      p1(2)=0.d0
      p1(3)=0.d0
      p1(4)=0.d0
      imin=nmin
      if(ir.eq.1)then
        i0=1
        i1=maproj
      else
        i0=maproj+1
        i1=maproj+matarg
      endif
c store projectile and target in p0 and sum pz an E in p1(3) and p1(4)
      do i=i0,i1
        if(istptl(i).eq.1.and.iorptl(i) .eq. 0 ) then
c calculate total energy of primaries
          do j=1,4
            p1(j)=p1(j)+dble(pptl(j,i))
          enddo
        endif
      enddo
      psoll=max(dble(pnullx),abs(p1(3)))
c calculate boost vector to have sum of pt=0 and sum of pz=0
      ppp=(p1(4)+p1(3))*(p1(4)-p1(3))-p1(2)*p1(2)-p1(1)*p1(1)
      if(ppp.gt.0.d0)then
        p1(5)=sqrt(ppp)
      elseif(abs(ppp/dble((i1-i0)*engy)).gt.1e-6)then
        iret=1
        write(ifch,*)'p1',p1(1),p1(2),p1(3),p1(4),ppp
        if(ish.ge.2)write (ifch,*) 'Problem in utrfrag : redo event'
        write (ifmt,*) 'Problem in utrfrag : redo event'
c       call utstop('utresc&')
        goto 9999
      else
        goto 9999
      endif
      esoll=p1(5)
      if(ish.ge.4) write (ifch,'(a,5g14.6)') 'boost-vector: ',p1

c     trafo
c     -----
      npart=0
      pt1soll=0.d0
      pt2soll=0.d0
      do i=imin,nptl
        call utlob4(1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $         ,pptl(1,i),pptl(2,i),pptl(3,i),pptl(4,i))
        npart=npart+1
        pt1soll=pt1soll+dble(pptl(1,i))
        pt2soll=pt2soll+dble(pptl(2,i))
      enddo

      if(ish.ge.6)then
        call clist('list after boost&',nmin,nptl,0,0)
      endif

      if(ish.ge.5)write(ifch,'(a)')'--------rescale momenta----------'

      if(ish.ge.6)write(ifch,*)'ptsoll:',pt1soll,pt2soll,npart
      pt1soll=pt1soll/dble(npart)
      pt2soll=pt2soll/dble(npart)
      do i=imin,nptl
        pptl(1,i)=pptl(1,i)-sngl(pt1soll)
        pptl(2,i)=pptl(2,i)-sngl(pt2soll)
        pptl(4,i)=sqrt(pptl(1,i)**2+pptl(2,i)**2
     &       +pptl(3,i)**2+pptl(5,i)**2)
      enddo

      if(ish.ge.6)write(ifch,*)'esoll,psoll,pmax:',esoll,psoll

c     rescale momenta in rest frame
c     -----------------------------
      scal=1.
      diff0=0.
      npart=nptl-imin+1
      do ipass=1,300
        sum=0.
        sum3=0.
        difft=diff0
        do  i=0,npart-1
          j=imin+i
c modified particles
          if(ipass.gt.1)then
            pptl3new=0.
            diff=difft/npart/pptl(5,j)
            if(ir.eq.1)then   !limit maximum momentum not to be too much above initial pmax
              if(diff.lt.0.)
     *        diff=max(diff, pptl(3,j)-exp(-4.*abs(pptl(3,j))))
            else
              if(diff.gt.0.)
     *        diff=min(diff, pptl(3,j)+exp(-4.*abs(pptl(3,j))))
            endif
            pptl(1,j)=scal*pptl(1,j)
            pptl(2,j)=scal*pptl(2,j)
            pptl(3,j)=scal*(pptl(3,j)-diff)
            pptl(4,j)=sqrt(pptl(1,j)**2+pptl(2,j)**2
     *           +pptl(3,j)**2+pptl(5,j)**2)
c           write(ifch,*)'par',j,pptl(4,j),pptl(3,j),diff,difft,scal
          endif
c sum over all particles
          sum=sum+pptl(4,j)
          sum3=sum3+pptl(3,j)
        enddo

        diff=sum3
        scal=sngl(esoll)/(sum-diff)
        if(ish.ge.6)write(ifch,*)
     $       'ipass,scal,diff/psoll,e,pz:'
     $   ,ipass,scal,diff/psoll,sum,sum3
        if(abs(scal-1.).le.errlim.and.abs(diff/psoll).lt.errlim)
     $  goto 300
        if(ipass.lt.200)then
          diff0=diff
        else
          diff0=0.
        endif
      enddo
      lim=.not.(abs(scal-1.).le.errlim.and.abs(diff/psoll).lt.errlim)
 
      if(ish.ge.1)then
        call utmsg('utrfrag&')
        write(ifch,*)'*****  scal=',scal,diff/psoll,lim
        call utmsgf
      endif

      if(lim)then
        if(ish.ge.1)then
          write(ifmt,*)'Warning in utrfrag !'
          write(ifch,'(a,i10,d25.15)')'redo momenta ...'
        endif
        iret=1
        goto 9999
      endif

c     trafo
c     -----
 300  continue

      if(ish.ge.6)then
        call clist('list after rescaling&',nmin,nptl,0,0)
      endif

      do i=imin,nptl
        call utlob4(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $           ,pptl(1,i),pptl(2,i),pptl(3,i),pptl(4,i))
      enddo

      if(ish.ge.4)call clist('list after utrfrag &',nmin,nptl,0,0)

 9999 continue
      call utprix('utrfrag',ish,ishini,4)

      end


c-----------------------------------------------------------------------
      subroutine utresc(iret)
c-----------------------------------------------------------------------
c  if irescl=1 rescaling is done, otherwise the purpose of going through
c  this routine is to not change the seed in case of no rescaling
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      double precision p1,esoll,ppp,seedp,psoll,pt1soll,pt2soll
     .,xx1,xx2,xx3,xx4,utdble
      dimension p1(5),p2(4),p0(5,mamx+mamx),pini(mxptl)
      logical force,nolead(mxptl),lspec(mxptl),lim
      data scalmean/0./scalevt/0./
      save scalmean,scalevt
      call utpri('utresc',ish,ishini,4)

      errlim=max(0.00005,0.05/engy)

      iret=0
      call ranfgt(seedp)        !not to change the seed ...
      if(nptl.le.nptlpt) goto 9999

      if(ish.ge.8)then
        call alistf('list before boost&')
      endif
      esoll=0.d0
      psoll=0.d0
      p1(1)=0.d0
      p1(2)=0.d0
      p1(3)=0.d0
      p1(4)=0.d0
      p2(1)=0.d0
      p2(2)=0.d0
      p2(3)=0.d0
      p2(4)=0.d0
      ipmax=4
      imin=nptlpt+1
      if(iappl.eq.1)then
        imin=1
        ipmax=0
      endif
c store projectile and target in p0 and sum pz an E in p1(3) and p1(4)
      do i=1,nptlpt
        nolead(i)=.true.
        do j=1,5
          p0(j,i)=dble(pptl(j,i))
        enddo
c calculate total energy of primaries
c       if(mod(istptl(i),10).eq.1)then
          do j=ipmax+1,4
            p1(j)=p1(j)+dble(pptl(j,i))
          enddo
c       endif
c calculate total energy of secondaries
        if(mod(istptl(i),10).eq.0)then
          do j=ipmax+1,4
            p2(j)=p2(j)+dble(pptl(j,i))
          enddo
        endif
      enddo
c fix secondaries counted in the system
      do i=nptlpt+1,nptl
       if(mod(istptl(i),10).eq.0)then
c check maximum energy
         if(pptl(4,i).gt.engy*0.5+3.
     &     .and.ityptl(i).ne.47.and.ityptl(i).ne.57.and.ityptl(i).lt.90
     &     .and.ityptl(i).ne.4.and.ityptl(i).ne.5)then !exclude active spectators
           if(ish.ge.1)write(ifch,*)'Energy of particle too high !'
     &                                     ,i,ityptl(i),pptl(4,i)
c           call utstop('utresc&')
           if(ityptl(i).eq.48.or.ityptl(i).eq.58)then  !diffractive resonance
             pptl(4,i)=0.5*engy
             amt=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(5,i)**2)
             si=sign(1.,pptl(3,i))
             pptl(3,i)=(pptl(4,i)+amt)*(pptl(4,i)-amt)
             if(pptl(3,i).gt.0.)then
               pptl(3,i)=si*sqrt(pptl(3,i))
             else
               iret=1
             endif
           else
             iret=1
           endif
           if(iret.eq.1)then
             if(ish.ge.1)write(ifch,*)'Warning in utresc: redo event...'
c             call utstop('Energy of particle too high !&')
             goto 9999
           endif
         endif
c fix pt (p1(1) and p2(1)) from secondaries
         do j=1,ipmax
           p1(j)=p1(j)+dble(pptl(j,i))
         enddo
c calculate total energy of secondaries
         do j=ipmax+1,4
           p2(j)=p2(j)+dble(pptl(j,i))
         enddo
         lspec(i)=.false.
c         if(((ityptl(i).eq.45.or.ityptl(i).eq.48).and.maproj.ge.100)
c     &  .or.((ityptl(i).eq.55.or.ityptl(i).eq.58).and.matarg.ge.100))
c     &   lspec(i)=.true.
c         if((abs(pptl(3,i)/pnullx).le.0.9
c        & .and.abs(pptl(3,i)).gt.pptl(5,i)).or.lspec(i))then
         if((abs(pptl(3,i)/pnullx).le.1.1
     &        .and.abs(pptl(3,i)).gt.pptl(5,i))
     &       .or.ityptl(i).eq.47.or.ityptl(i).eq.47
     &       .or.ityptl(i).eq.94.or.ityptl(i).eq.95)then
           nolead(i)=.true.
           if(ish.ge.8)write(ifch,*)'nolead',i
         else
           nolead(i)=.false.
           if(ish.ge.8)write(ifch,*)'lead',i
         endif
       endif
      enddo
      psoll=max(dble(pnullx),abs(p1(3)))
c check if energy is conserved before boost
      if(iappl.eq.1)then
        diff=abs(p2(3)-p1(3))
        scal=p2(4)/p1(4)
        if(abs(scal-1.).le.errlim.and.abs(diff/psoll).lt.errlim
     &     .and.(iLHC.eq.0.or.
     &    (abs(p2(1)).lt.errlim.and.abs(p2(2)).lt.errlim)))then
          if(ish.ge.4) 
     & write (ifch,'(2x,a,2g14.6)') 'Energy OK: ',scal,abs(diff/psoll)
          goto 9999
        else
         diff=0.
         scal=1.
        endif
      endif
c calculate boost vector to have sum of pt=0 and sum of pz=0
      ppp=(p1(4)+p1(3))*(p1(4)-p1(3))-p1(2)*p1(2)-p1(1)*p1(1)
      if(ppp.gt.0.d0)then
        p1(5)=sqrt(ppp)
      else
        iret=1
        write(ifch,*)'p1',p1(1),p1(2),p1(3),p1(4),ppp
        if(ish.ge.2)write (ifch,*) 'Problem in utresc (1): redo event'
        write (ifmt,*) 'Problem in utresc (1): redo event'
c       call utstop('utresc&')
        goto 9999
      endif
      esoll=p1(5)
      if(ish.ge.4) write (ifch,'(a,5g14.6)') 'boost-vector: ',p1

c     trafo
c     -----
      pmax=0.
      npart=0
      pt1soll=0.d0
      pt2soll=0.d0
      do i=imin,nptl
        if(mod(istptl(i),10).le.1)then
          call utlob4(1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $         ,pptl(1,i),pptl(2,i),pptl(3,i),pptl(4,i))
          if(mod(istptl(i),10).eq.0.and.i.gt.nptlpt.and.nolead(i).and.
     & .not.((ityptl(i).eq.4.or.ityptl(i).eq.5).and.istptl(i).eq.0))then
            npart=npart+1
            pt1soll=pt1soll+dble(pptl(1,i))
            pt2soll=pt2soll+dble(pptl(2,i))
          endif
        endif
        if(i.gt.nptlpt.and.nolead(i))pmax=max(pmax,abs(pptl(3,i)))
      enddo

      if(ish.ge.6)then
        call alistf('list after boost&')
      endif

      if(ish.ge.5)write(ifch,'(a)')'--------rescale momenta----------'

        if(ish.ge.6)write(ifch,*)'ptsoll:',pt1soll,pt2soll,npart,errlim
        pt1soll=pt1soll/dble(npart)
        pt2soll=pt2soll/dble(npart)
        do i=nptlpt+1,nptl
          if(mod(istptl(i),10).eq.0.and.nolead(i).and.
     & .not.((ityptl(i).eq.4.or.ityptl(i).eq.5).and.istptl(i).eq.0))then
            pptl(1,i)=pptl(1,i)-sngl(pt1soll)
c     .               +2.*max(0.,pptl(1,i)-0.5)
            pptl(2,i)=pptl(2,i)-sngl(pt2soll)
c     .               +2.*max(0.,pptl(2,i)-0.5)
            pptl(4,i)=sqrt(pptl(1,i)**2+pptl(2,i)**2
     &                    +pptl(3,i)**2+pptl(5,i)**2)
          endif
        enddo

      if(ish.ge.6)write(ifch,*)'esoll,psoll,pmax:',esoll,psoll,pmax

c     rescale momenta in rest frame
c     -----------------------------
      scal=1.
      diff0=0.
      ndif0=1
      ferr=0.1
      force=.false.
      npart=nptl-imin+1
      do ipass=1,300
        sum=0.
        sum3=0.
        difft=diff0
        ndif=0
        nfirst=int(rangen()*float(npart))    !start list randomly
        do  i=0,npart-1
          j=imin+i+nfirst
          if(j.gt.nptl)j=imin+i+nfirst-npart
c          write(ifch,*)j,istptl(j),nolead(j),lspec(j)
          if(mod(istptl(j),10).eq.0)then
c modified particles
            if(nolead(j))then
c            if(j.gt.nptlpt)then
c            if(abs(pptl(3,j))/pnullx.lt.0.9)then  !not spectator or diffraction
c              if(scal.eq.1..and.abs(diff0).lt.1.e-6)then
              if(ipass.eq.1)then
                ndif=ndif+1
                pini(j)=pptl(3,j)
              else
                pptl3new=0.
                if( force .or.(
     &            ityptl(j)/10.eq.4.or.ityptl(j)/10.eq.5
     &               ))then     !try just remnant first
                  ndif=ndif+1
                  diff=sign(min(0.3*abs(pini(j)),
     &                   2.*rangen()*abs(difft)/max(1,ndif0)),difft)
                  pptl3new=scal*(pptl(3,j)-diff)
c                  write(ifch,*)'par',j,pptl3new,pptl(3,j),diff,difft
c     &                 ,ndif,pmax,scal
                  if(abs(pptl3new).lt.pmax)then
c particle should not be too fast or too modified
                    if(abs(pptl3new-pini(j)).lt.ferr*abs(pini(j))
     &         .or.(lspec(j).and.abs(pptl3new).lt.abs(0.8*pini(j))))then
c                  write(ifch,*)'used'
                      difft=difft-diff
                      pptl(3,j)=scal*(pptl(3,j)-diff)
                      pptl(4,j)=sqrt(pptl(1,j)**2+pptl(2,j)**2
     *                     +pptl(3,j)**2+pptl(5,j)**2)
                    endif
                  endif
                endif
              endif
            endif
c sum over all particles
            sum=sum+pptl(4,j)
            sum3=sum3+pptl(3,j)
          endif
        enddo

        diff=sum3
        scal=sngl(esoll)/(sum-diff)
        if(ish.ge.6)write(ifch,*)
     $       'ipass,scal,diff/psoll,e,pz,ndif,f:'
     $   ,ipass,scal,diff/psoll,sum,sum3,ndif,force
        if(abs(scal-1.).le.errlim.and.abs(diff/psoll).lt.errlim)
     $  goto 300
        if(ndif.gt.0.and.(force.or.ipass.lt.150))then
          ndif0=ndif
          diff0=diff
        elseif(force)then
          if(ish.ge.2)
     $    write(ifmt,*)'Warning in utresc: no more allowed particle'
          goto 302
        else
          force=.true.
          ferr=0.5
          diff0=0.
          ndif0=ndif
        endif
      enddo
 302  if(iLHC.eq.1)then
        lim=.not.(abs(scal-1.).le.errlim.and.abs(diff/psoll).lt.errlim)
      else
        lim=abs(scal)+abs(diff/psoll).gt.2.5
      endif
 
      if(ish.ge.1)then
        call utmsg('utrescl&')
        write(ifch,*)'*****  scal=',scal,diff/psoll,lim
        call utmsgf
      endif

      if(lim)then
        if(ish.ge.1)then
          write(ifmt,*)'Warning in utresc !'
          write(ifch,'(a,i10,d25.15)')'Warning: redo EPOS event ...'
     &                                ,nint(seedj),seedc
        endif
        iret=1
        goto 9999
      endif

c     trafo
c     -----
 300  continue

      do i=1,nptl
        if(i.le.nptlpt)then
          do j=1,5
            pptl(j,i)=p0(j,i)
          enddo
        else
          if(mod(istptl(i),10).le.1)then
c          write(ifch,*)'av',i,idptl(i),istptl(i)
c     $  ,pptl(1,i),pptl(2,i),pptl(3,i),pptl(4,i)
c     $  ,(pptl(4,i)+pptl(3,i))*(pptl(4,i)-pptl(3,i))
c     $  -pptl(1,i)*pptl(1,i)-pptl(2,i)*pptl(2,i)
            xx1=utdble(pptl(1,i))
            xx2=utdble(pptl(2,i))
            xx3=utdble(pptl(3,i))
            xx4=sqrt(xx1**2+xx2**2+xx3**2+utdble(pptl(5,i))**2)
            call utlob2(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $                                 ,xx1,xx2,xx3,xx4,70)
            pptl(4,i)=sngl(sqrt(xx1**2+xx2**2+xx3**2
     $                         +utdble(pptl(5,i))**2))
            if(abs(pptl(4,i)-sngl(xx4)).gt.1e-4*pptl(4,i))
     $                      call utstop('Pb in utresc !&')
            pptl(3,i)=sngl(xx3)
            pptl(2,i)=sngl(xx2)
            pptl(1,i)=sngl(xx1)
c          write(ifch,*)'ap',i,idptl(i),istptl(i)
c     $  ,pptl(1,i),pptl(2,i),pptl(3,i),pptl(4,i)
c     $  ,(pptl(4,i)+pptl(3,i))*(pptl(4,i)-pptl(3,i))
c     $  -pptl(1,i)*pptl(1,i)-pptl(2,i)*pptl(2,i)
          endif
        endif
      enddo

      if(ish.ge.4)call alist('list after rescaling&',1,nptl)

 9999 continue
      if(ish.ge.2)then
        scalevt=scalevt+1.
        scalmean=scalmean+scal
        write(ifch,*)' average rescaling factor: ',scalmean
     &                                            /scalevt
      endif
      call ranfst(seedp)        ! ... after this subroutine
      call utprix('utresc',ish,ishini,4)

      end

c-----------------------------------------------------------------------
      subroutine utghost(iret)
c-----------------------------------------------------------------------
c  if irescl=1 make particle on-shell if not
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      double precision seedp
      call utpri('ughost',ish,ishini,4)

      iret=0
      call ranfgt(seedp)        !not to change the seed ...
      if(nptl.le.nptlpt) goto 9999

      if(ish.ge.5)write(ifch,'(a)')'---------mark ghosts---------'

c     mark ghosts
c     -----------
      do  j=nptlpt+1,nptl
        if(istptl(j).le.1.and.pptl(4,j).gt.0.d0)then
          if(istptl(j).eq.0)then !fix max of final particles only
            amass=pptl(5,j)
            call idmass(idptl(j),amass)
            call idwidth(idptl(j),wi)
            del=max(1.e-3,2.*wi)
            if(abs(idptl(j)).gt.100.and.idptl(j).ne.221.and.
     &       idptl(j).ne.702.and.abs(pptl(5,j)-amass).gt.del)then
              if(ish.ge.5)write(ifch,*)'wrong particle mass',j,idptl(j)
     &                                           ,pptl(5,j),amass,del
              amass=pptl(5,j)
              call idres(idptl(j),amass,idr,iadj,0)
              if(idr.ne.0)then
                if(iadj.eq.1)pptl(5,j)=amass
                idptl(j)=idr
              else
                call idmass(idptl(j),amass)
                pptl(5,j)=amass
              endif
              call idtau(idptl(j),pptl(4,j),pptl(5,j),taugm)
              tivptl(2,j)=tivptl(1,j)+taugm*(-alog(rangen()))
            elseif(idptl(j).eq.221.or.idptl(j).eq.702)then    !particles coming from mix with rho0 or eta'   (omega and f0 width smaller than rho0 or eta' width, so fix the mass tothe average one)
              pptl(5,j)=amass
            endif
          endif
          if(abs((pptl(4,j)+pptl(3,j))*(pptl(4,j)-pptl(3,j))
     $         -pptl(2,j)**2-pptl(1,j)**2-pptl(5,j)**2).gt.0.3
     $       .and.abs(1.-abs(pptl(3,j))/pptl(4,j)).gt.0.01)then
        !print*,'ghost',ityptl(j),idptl(j)
           if(ish.ge.1)write(ifmt,*)'ghost:',j,idptl(j),ityptl(j)
           if(ish.ge.5)then
              write(ifch,'(a,$)')'ghost:'
              call alistc("&",j,j)
            endif
            ityptl(j)=100+ityptl(j)/10
          else
c ensure that all particles are really on-shell
            pptl(4,j)=sqrt(pptl(1,j)**2+pptl(2,j)**2
     *                    +pptl(3,j)**2+pptl(5,j)**2)
          endif
        elseif(mod(istptl(j),10).eq.0)then
c if not droplet with fusion
          if(istptl(j).ne.10.or.iorsdf.ne.3)then
            if(ish.ge.1)then
              write(ifmt,*)'Lost particle (E=0)'
              write(ifch,*)'Lost particle (E=0) :'
              call alistc("utghost&",j,j)
            endif
            istptl(j)=istptl(j)+14
          endif
        endif
      enddo

      if(ish.ge.5)write(ifch,'(a)')'---------treat ghosts---------'

c     treat ghosts
c     ------------
      ifirst=1
      scal=1.
      pfif=0.
      efif=0.
      ntry=0
 132  nfif=0
      psum=0
      esum=0.
      ntry=ntry+1
      do  j=nptlpt+1,nptl
        if(mod(istptl(j),10).eq.0)then
          if(ityptl(j).gt.100)then
            nfif=nfif+1
            if(ifirst.eq.1)then
              pfif=pfif+pptl(3,j)
              if(pptl(4,j).gt.0.)efif=efif+pptl(4,j)
            endif
            if(irescl.ge.1) then
              if(ifirst.gt.1)then
                if(pptl(4,j).gt.0.)then
                  Einv=1./pptl(4,j)
                  amt=1.-(pptl(5,j)*Einv)**2+(pptl(1,j)*Einv)**2
     $                +(pptl(2,j)*Einv)**2
                else
                  amt=-1.
                endif
                if(amt.gt.0.)then
                  pptl(3,j)=sign(pptl(4,j),pptl(3,j))*sqrt(amt)
                else
                  y=(rangen()+rangen()+rangen()+rangen()-2.)/2.*yhaha
                  y=sign(abs(y),pptl(3,j))
                  pptl(3,j)
     $                 =sqrt(pptl(5,j)**2+pptl(1,j)**2
     $                 +pptl(2,j)**2)*sinh(y)
                  pptl(4,j)
     $                 =sqrt(pptl(5,j)**2+pptl(1,j)**2
     $                 +pptl(2,j)**2)*cosh(y)
                  efif=efif+pptl(4,j)
                endif
                ifirst=0
              else
c                do k=1,3
                do k=3,3
                  pptl(k,j)=pptl(k,j)*scal
                enddo
                pptl(4,j)=sqrt(pptl(1,j)**2+pptl(2,j)**2+pptl(3,j)**2
     *                 +pptl(5,j)**2)
              endif
            endif
            psum=psum+pptl(3,j)
            esum=esum+pptl(4,j)
            if(ish.ge.5)
     $           write (ifch,*) 'nrevt,psum,esum,pfif,efif,nfif,scal'
     $           ,nrevt,psum,esum,pfif,efif,nfif,scal
          endif
        endif
      enddo
      if ( ish.ge.5 )  write (ifch,*) 'tot',nfif,efif,pfif,esum,psum


      if(nfif.gt.5.or.(esum.gt.0.05*engy.and.nfif.ne.1))then
        if(ifirst.eq.0)then
          do  j=nptlpt+1,nptl
            if ( ityptl(j).ge.101 .and. ityptl(j).le.105 )then
              if((psum-pfif)*(1.-scal).ge.0)
     &             pptl(3,j)=pptl(3,j)-(psum-pfif)/nfif
            endif
          enddo
        else
          ifirst=2
          goto 132
        endif
        scal=efif/esum
        if ( ish.ge.5 )  write (ifch,*) 'scal',scal
        if ( abs(scal-1.) .gt. 0.05 ) then
          if(ntry.le.1000)then
            goto 132
          else
            iret=1
            if(ish.ge.2)write (ifch,*) 'Problem in utghost : redo event'
            if(ish.ge.1)write (ifmt,*) 'Problem in utghost : redo event'
            goto 9999
         endif
        endif
      endif

c finally put ghost on-shell
      do  j=nptlpt+1,nptl
        if ( ityptl(j).ge.101 .and. ityptl(j).le.105 )then
          pptl(4,j)=sqrt(pptl(1,j)**2+pptl(2,j)**2+pptl(3,j)**2
     *                 +pptl(5,j)**2)
        endif
      enddo
      

      if(ish.ge.5)write(ifch,'(a)')'---------Check Ghost list---------'

c Check Ghost list

      if(ish.ge.5)then
        do  j=nptlpt+1,nptl
          if(mod(istptl(j),10).eq.0)then
            if(ityptl(j).le.105.and.ityptl(j).ge.101)then
              write(ifch,'(a,$)')'ghost:'
              call alistc("&",j,j)
            endif
          endif
        enddo
      endif


 9999 continue
      call ranfst(seedp)        ! ... after this subroutine
      call utprix('ughost',ish,ishini,4)

      end

c-----------------------------------------------------------------------
      subroutine utrsph(iret)
c-----------------------------------------------------------------------
c  if irescl=1 and ispherio=1 rescaling is done for particle used by
c  spherio as initial condition.
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      double precision p1,esoll,utdble
      dimension p1(5),p0(5,mamx+mamx)
      call utpri('utrsph',ish,ishini,4)

      errlim=0.0001

      iret=0
      if(nptl.le.nptlpt) goto 9999

      esoll=0.d0
      p1(1)=0.d0
      p1(2)=0.d0
      p1(3)=0.d0
      p1(4)=0.d0
      do i=nptlpt+1,nptl
        if((istptl(i).le.11
     $   .and.(iorptl(i).ge.1.and.istptl(iorptl(i)).eq.41))
     $   .or.istptl(i).eq.20.or.istptl(i).eq.21)then
         do j=1,2
           p1(j)=p1(j)+utdble(pptl(j,i))
         enddo
       endif
      enddo
      do i=1,nptlpt
         do j=1,5
           p0(j,i)=pptl(j,i)
         enddo
         do j=3,4
           p1(j)=p1(j)+utdble(pptl(j,i))
         enddo
      enddo
      p1(5)=dsqrt((p1(4)+p1(3))*(p1(4)-p1(3))-p1(2)**2.d0-p1(1)**2.d0)
      esoll=p1(5)
      if(ish.ge.4) write (ifch,'(a,5g13.6)') 'boost-vector',p1

c     trafo
c     -----
      do i=1,nptl
        if((istptl(i).le.11
     $   .and.(iorptl(i).ge.1.and.istptl(iorptl(i)).eq.41))
     $   .or.istptl(i).eq.20.or.istptl(i).eq.21
     $   .or.(istptl(i).eq.0.and.i.le.nptlpt))then
          call utlob4(1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $         ,pptl(1,i),pptl(2,i),pptl(3,i),pptl(4,i))
        endif
      enddo


      if(ish.ge.5)write(ifch,'(a)')'------------------'

c     rescale momenta in rest frame
c     -----------------------------

      scal=1.
      diff=0.
      do ipass=1,1000
        sum=0.
        sum3=0.
        ndif=0
        do  j=1,nptl
        if((istptl(j).le.11
     $   .and.(iorptl(j).ge.1.and.istptl(iorptl(j)).eq.41))
     $   .or.istptl(j).eq.20.or.istptl(j).eq.21
     $   .or.(istptl(j).eq.0.and.j.le.nptlpt))then
            if(j.gt.nptlpt)then
              ndif=ndif+1
              pptl(3,j)=scal*(pptl(3,j)-diff)
              pptl(4,j)=sqrt(pptl(1,j)**2+pptl(2,j)**2+pptl(3,j)**2
     *           +pptl(5,j)**2)
            endif
            sum=sum+pptl(4,j)
            sum3=sum3+pptl(3,j)
          endif
        enddo

        diff=sum3/real(ndif)
        scal=real(esoll)/sum
        if(ish.ge.6)write(ifch,*)'ipass,scal,diff,e,esoll,pz,ndif:'
     $       ,ipass,scal,diff,sum,esoll,sum3,ndif
        if(abs(scal-1.).le.errlim.and.abs(diff).lt.10.*errlim) goto300
      enddo
      if(ish.ge.1)then
        call utmsg('hresph&')
        write(ifch,*)'*****  scal=',scal,diff
        call utmsgf
      endif


c     trafo
c     -----
 300  continue
c      do i=nptlpt+1,nptl
      do i=1,nptl
        if((istptl(i).le.11
     $   .and.(iorptl(i).ge.1.and.istptl(iorptl(i)).eq.41))
     $   .or.istptl(i).eq.20.or.istptl(i).eq.21
     $   .or.(istptl(i).eq.0.and.i.le.nptlpt))then
          call utlob4(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $         ,pptl(1,i),pptl(2,i),pptl(3,i),pptl(4,i))
        endif
        if(i.le.nptlpt)then
          do j=1,5
            pptl(j,i)=p0(j,i)
          enddo
        endif
      enddo

 9999 call utprix('utrsph',ish,ishini,4)

      end

cc-----------------------------------------------------------------------
c      double precision function dddlog(xxx)
cc-----------------------------------------------------------------------
c      double precision xxx
c      dddlog=-1d50
c      if(xxx.gt.0d0)dddlog=dlog(xxx)
c      end
c
ccc-----------------------------------------------------------------------
c      subroutine randfl(jc,iqa0,iflav,ic,isame)
cc-----------------------------------------------------------------------
cc     returns random flavour ic(2) (iqa0=1:quark,2:antiquark,11:diquark)
cc-----------------------------------------------------------------------
c      include 'epos.inc'
c      real probab(nflav),probsu(nflav+1)
c      integer jc(nflav,2),jc0(nflav,2),ic(2)
c      if(ish.ge.6)then
c      write(ifch,*)('-',i=1,10)
c     *,' entry sr randfl ',('-',i=1,30)
c      write(ifch,*)'iqa0:',iqa0
c      write(ifch,*)'jc:'
c      write(ifch,*)jc
c      endif
c      iflav=0
c      ic(1)=0
c      ic(2)=0
c      do 10 n=1,nflav
c      do 10 i=1,2
c10    jc0(n,i)=0
c      iqa1=iqa0*10
c9999  iqa1=iqa1/10
c      if(iqa1.eq.0)goto9998
c      iqa=mod(iqa1,10)
c      su=0
c      do 20 i=1,nflav
c      probab(i)=jc(i,iqa)-jc0(i,iqa)
c      if(isame.eq.1)probab(i)=probab(i)*(jc(i,3-iqa)-jc0(i,3-iqa))
c20    su=su+probab(i)
c      if(su.lt..5)then
c      iflav=0
c      ic(1)=0
c      ic(2)=0
c      goto9998
c      endif
c      probsu(1)=0.
c      do 30 i=1,nflav
c      probsu(i+1)=probsu(i)+probab(i)/su
c      if(probsu(i+1)-probsu(i).lt.1e-5)probsu(i+1)=probsu(i)
c30    continue
c      r=rangen()*probsu(nflav+1)
c      do 50 i=1,nflav
c      if(probsu(i).le.r.and.r.lt.probsu(i+1))iflav=i
c50    continue
c      jc0(iflav,iqa)=jc0(iflav,iqa)+1
c      if(isame.eq.1)jc0(iflav,3-iqa)=jc0(iflav,3-iqa)+1
c      call idenco(jc0,ic,ireten)
c      if(ireten.eq.1)call utstop('randfl: idenco ret code = 1&')
c      if(ish.ge.6)then
c      write(ifch,*)'probab:'
c      write(ifch,*)probab
c      write(ifch,*)'probsu:'
c      write(ifch,*)probsu
c      write(ifch,*)'ran#:',r,'   flav:',iflav
c      endif
c      goto9999
c9998  continue
c      if(ish.ge.6)write(ifch,*)('-',i=1,30)
c     *,' exit sr randfl ',('-',i=1,10)
c      return
c      end
c
c
cc-----------------------------------------------------------------------
c      subroutine ranhvy(x,eps)
cc-----------------------------------------------------------------------
cc     generates x for heavy particle fragmentation according to
cc     the peterson form
cc          d(x)=1/(x*(1-1/x-eps/(1-x))**2)
cc              =d0(x)*d1(x)*d2(x)
cc          d0(x)=(1-x)**2/((1-x)**2+eps)**2
cc          d1(x)=x
cc          d2(x)=(((1-x)**2+eps)/((1-x)**2+eps*x))**2
cc     using x=1-y**pow
cc     generates flat in x if eps>1.
cc-----------------------------------------------------------------------
c      data aln4/1.3863/
c      if(eps.lt.1.) then
c        pow=alog((3.+eps)/eps)/aln4
c        ymx=(eps*(3.*pow-1.)/(pow+1.))**(.5/pow)
c        zmx=1-ymx**pow
c        d0mx=(1-zmx)**2/((1.-zmx)**2+eps)**2*pow*ymx**(pow-1.)
c        d2mx=2./(2.-sqrt(eps))
c      else
c        pow=1.
c        zmx=0.
c        d0mx=(1.-zmx)**2/((1.-zmx)**2+eps)**2
c        d2mx=1.+eps
c      endif
cc
cc          generate z according to (1-z)**2/((1-z)**2+eps*z)**2
c1     continue
c      y=rangen()
c      z=1.-y**pow
cc
c      d0z=(1.-z)**2/((1.-z)**2+eps)**2*pow*y**(pow-1.)
c      if(d0z.lt.rangen()*d0mx) goto1
cc
cc          check remaining factors
c      d1=z
c      d2=(((1.-z)**2+eps)/((1.-z)**2+eps*z))**2
c      if(d1*d2.lt.rangen()*d2mx) goto1
cc
cc          good x
c      x=z
c      return
c      end
c
c=============================================================================
      subroutine qIIgfrag(sa,na,rc,nsf,iaf)
c-----------------------------------------------------------------------------
c connected nucleon clusters extraction - used for the nuclear spectator part
c multifragmentation
c     adapted code from QGSJETII from Sergey Ostapchenko
c     based on : N.N. Kalmykov and S.S. Ostapchenko Yad. Fiz. 56, 105-119 (March 1993)
c-----------------------------------------------------------------------------
      implicit none
      include 'epos.inc'
      real sa(mamxx,4),rc,ri,s0
      integer na,nsf,iaf(mamxx),i,j,l,ni,ng,j1,m

      if(ish.ge.2)write (ifch,201)na
      if(ish.ge.6)then
       write (ifch,203)
       do i=1,na
        write (ifch,204)(sa(i,l),l=1,4)
       enddo
      endif

      ni=1
      ng=nint(sa(ni,4))
      j=0
1     j=j+1
      j1=ni+1

      do 4 i=j1,na
       ri=0.d0
       do m=1,3
        ri=ri+(sa(j,m)-sa(i,m))**2
       enddo
       if(ri.gt.rc)goto 4

       ni=ni+1
       if(i.ne.ni)then
         do m=1,3
           s0=sa(ni,m)
           sa(ni,m)=sa(i,m)
           sa(i,m)=s0
         enddo
         s0=sa(ni,4)
         sa(ni,4)=sa(i,4)
         sa(i,4)=s0
       endif
       ng=ng+nint(sa(ni,4))
4     continue

      if(j.lt.ni.and.na-ni.gt.0)goto 1
      nsf=nsf+1
      iaf(nsf)=ng
      if(ish.ge.3)write (ifch,206)nsf,iaf(nsf)

      j=ni
      ni=ni+1
      ng=nint(sa(ni,4))
      if(na.eq.ni)then
       nsf=nsf+1
       iaf(nsf)=ng
       if(ish.ge.3)write (ifch,206)nsf,iaf(nsf)
      elseif(na.gt.ni)then
       goto 1
      endif

      if(ish.ge.3)write (ifch,202)
201   format(2x,'qIIgfrag-multifragmentation: nucleus mass number: na='
     *,i2)
202   format(2x,'qIIgfrag - end')
203   format(2x,'nucleons coordinates:')
204   format(2x,4e10.3)
206   format(2x,'qIIgfrag: fragment n',i2,2x,'fragment mass : ',i5)
      return
      end

*-- Author :    D. HECK IK FZK KARLSRUHE       27/04/1994
C=======================================================================

      SUBROUTINE EPOVAPOR(MAPRO,INEW,ITHEO,XA,JFIN,ITYP,PFRX,PFRY,PFRZ)

C-----------------------------------------------------------------------
C  (E)VAPOR(ATION OF NUCLEONS AND ALPHA PARTICLES FROM FRAGMENT)
C
C  TREATES THE REMAINING UNFRAGMENTED NUCLEUS
C  EVAPORATION FOLLOWING CAMPI APPROXIMATION.
C  SEE: X. CAMPI AND J. HUEFNER, PHYS.REV. C24 (1981) 2199
C  AND  J.J. GAIMARD, THESE UNIVERSITE PARIS 7, (1990)
C  THIS SUBROUTINE IS CALLED FROM SDPM, DPMJST, NSTORE, AND VSTORE.
C  ARGUMENTS INPUT:
C   MAPRO        = NUMBER OF NUCLEONS OF PROJECTILE
C   INEW         = PARTICLE TYPE OF SPECTATOR FRAGMENT
C  ARGUMENTS OUTPUT:
C   JFIN         = NUMBER OF FRAGMENTS
C   ITYP(1:JFIN) = NATURE (PARTICLE CODE) OF FRAGMENTS
C   PFRX(1:JFIN) = TRANSVERSE MOMENTUM OF FRAGMENTS IN X-DIRECTION
C   PFRY(1:JFIN) = TRANSVERSE MOMENTUM OF FRAGMENTS IN Y-DIRECTION
C   PFRZ(1:JFIN) = LONGITUDINAL MOMENTUM OF FRAGMENTS IN Z-DIRECTION
C
C  FROM CORSIKA AND ADAPTED BY T. PIEROG TO INCLUDE LONG MOMENTUM AND
C  MORE REALISTIC FRAGMENTS IN PARTICULAR WITH MULTIFRAGMENTATION
C-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      common/eporansto2/irndmseq
      integer irndmseq
      common /cfragm/rminfrg,emaxfrg,facgrey,p3grey
      real rminfrg,emaxfrg,facgrey,p3grey
      real  xproj,yproj,zproj,xtarg,ytarg,ztarg
      real  xa(mamxx,4)

      DOUBLE PRECISION PFR(mamxx),PFRX(mamxx),PFRY(mamxx),PFRZ(mamxx)
     *                ,RD(2*mamxx),SPFRY,SPFRZ,drangen,zmin
      DOUBLE PRECISION AFIN,AGLH,APRF,BGLH,EEX,PHIFR,RANNORM,SPFRX,FAC
      INTEGER          ITYP(mamxx),IARM,INEW,ITYPRM,INRM,IZRM,JC,lseq
     *     ,JFIN,K,L,MAPRO,MF,NFIN,NINTA,NNUC,NPRF,NNSTEP
     *     ,ITHEO,I!,IS,LS
      SAVE
      EXTERNAL         RANNORM
C-----------------------------------------------------------------------

      IF(ish.ge.7)WRITE(ifch,*)'EPOVAPOR : MAPRO,INEW=',MAPRO,INEW,ITHEO

      IF(MAPRO.LT.0)THEN
        NFIN = JFIN
        JFIN = 0
        GOTO 25
      ENDIF

      JC     = 0
      JFIN   = 0
      lseq   = irndmseq
      ITYPRM = INEW
      NPRF   = INEW/100
      NINTA  = MAPRO - ITHEO  !wounded nucleons
c      print *,mapro,nprf,ninta,itheo

c Code adapted from QGSJETII
      eex=0.d0                         !excitation energy for spectator part
           !sum of excitations due to wounded nucleons (including diffractive)
      do i=1,NINTA
c partial excitation according to f(e) ~ 1/sqrt(e) * exp(-e/(2*<e>))
       eex=eex+(drangen(eex)+drangen(eex)+drangen(eex)+
     * drangen(eex)+drangen(eex)-2.5d0)**2*2.4d0
      enddo
      if(ish.ge.3)write (ifch,203)eex,nprf
203   format(2x,'EPOVAPOR: excitation energy: eex=',e10.3,' nprf=',i4)

      if(eex/NPRF.gt.emaxfrg)then !if eex>emaxqII -> multifragmentation
        
        call qIIgfrag(xa,NPRF,rminfrg,JFIN,ITYP) !multifragmentation (percolation algorithm)
        
      else

      IF ( NINTA .EQ. 0 ) THEN
C  NO NUCLEON HAS INTERACTED
        JFIN    = 1
        PFRX(1)  = 0.D0
        PFRY(1)  = 0.D0
        PFRZ(1)  = 0.D0
        ITYP(1) = ITYPRM
        IF(ish.ge.7)WRITE(ifch,*) 'EPOVAPOR : JFIN,NINTA=',JFIN,NINTA
        GOTO 50
      ENDIF

C  EXCITATION ENERGY EEX OF PREFRAGMENT
C  SEE: J.J. GAIMARD, THESE UNIVERSITE PARIS 7, (1990), CHPT. 4.2
      EEX = 0.D0
      CALL RMMARD( RD,2*NINTA,lseq )
      DO  L = 1, NINTA
        IF ( RD(NINTA+L) .LT. RD(L) ) RD(L) = 1.D0 - RD(L)
        EEX = EEX + RD(L)
      ENDDO
C  DEPTH OF WOODS-SAXON POTENTIAL TO FERMI SURFACE IS 0.040 GEV
      IF(ish.ge.7)WRITE(ifch,*)'EPOVAPOR : EEX=',SNGL(EEX*0.04D0),
     &                                            ' GEV'
C  EVAPORATION: EACH EVAPORATION STEP NEEDS ABOUT 0.020 GEV, THEREFORE
C  NNSTEP IS EEX * 0.04/0.02 = EEX * 2.
      NNSTEP = INT( EEX*2.5D0 )      !fine tuned to NA61 data

      
      IF ( NNSTEP .LE. 0 ) THEN
C  EXCITATION ENERGY TOO SMALL, NO MULTIFRAGMENTATION
        JFIN    = 1
        PFRX(1)  = 0.D0
        PFRY(1)  = 0.D0
        PFRZ(1)  = 0.D0
        ITYP(1) = ITYPRM
        IF(ish.ge.7)WRITE(ifch,*) 'EPOVAPOR : JFIN,EEX=',JFIN,SNGL(EEX)
        GOTO 50

      ENDIF

C  AFIN IS ATOMIC NUMBER OF FINAL NUCLEUS
      APRF = DBLE(NPRF)
      AFIN = APRF - 1.6D0 * DBLE(NNSTEP)
      NFIN = MAX( 0, INT( AFIN+0.5D0 ) )
C  CORRESPONDS TO DEFINITION; FRAGMENTATION-EVAPORATION
C  CONVOLUTION EMU07 /MODEL ABRASION EVAPORATION (JNC FZK APRIL 94)
C  NNUC IS NUMBER OF EVAPORATING NUCLEONS
      NNUC = NPRF - NFIN
      IF(ish.ge.7)WRITE(ifch,*) 'EPOVAPOR : NFIN,NNUC=',NFIN,NNUC

      IF     ( NNUC .LE. 0 ) THEN
C  NO EVAPORATION
        JFIN    = 1
        PFRX(1)  = 0.D0
        PFRY(1)  = 0.D0
        PFRZ(1)  = 0.D0
        ITYP(1) = ITYPRM
        GOTO 50

      ELSE

        IARM = ITYPRM/100
        DO WHILE (NNUC .GT. 0)

          IF ( IARM .LE. 0 .OR. NNUC .LE. 0 ) GOTO 100
          IZRM = MOD(ITYPRM,100)
          INRM = IARM - IZRM
          JC   = JC + 1
          CALL RMMARD( RD,2,lseq )
          IF ( NNUC.GE.4 .AND. IZRM.GE.2 .AND. INRM.GE.2
     .         .AND. RD(1) .LT. 1.2D0 ) THEN
            ITYP(JC) = 402          !alpha
            NNUC     = NNUC - 4
          ELSEIF ( NNUC.GE.3 .AND.IZRM .GE.1 .AND. INRM.GE.3 ) THEN
            ITYP(JC) = 301          !TRITIUM
            NNUC     = NNUC - 3
          ELSEIF ( NNUC.GE.2 .AND. IZRM.GE.1 .AND. INRM.GE.2 ) THEN
            ITYP(JC) = 201          !DEUTERIUM
            NNUC     = NNUC - 2
          ELSEIF ( (RD(2) .LT. 0.5 .AND. IZRM .GE. 1)
     .        .OR. INRM.EQ.0  ) THEN
            ITYP(JC) = 101          !PROTRON
            NNUC     = NNUC - 1
          ELSE
            ITYP(JC) = 100          !NEUTRON
            NNUC     = NNUC - 1
          ENDIF

          ITYPRM   = ITYPRM - ITYP(JC)
          IARM = ITYPRM/100

        ENDDO

        
c      ELSEIF ( NNUC .GE. 4 ) THEN
cC  EVAPORATION WITH FORMATION OF ALPHA PARTICLES POSSIBLE
cC  IARM, IZRM, INRM ARE NUMBER OF NUCLEONS, PROTONS, NEUTRONS OF
cC  REMAINDER
c        DO  LS = 1, NNSTEP
c          IARM = ITYPRM/100
c          IF ( IARM .LE. 0 .OR. NNUC .LE. 0 ) GOTO 100
c          IZRM = MOD(ITYPRM,100)
c          INRM = IARM - IZRM
c          JC   = JC + 1
c          CALL RMMARD( RD,2,lseq )
c          IF ( RD(1) .LT. 0.2D0  .AND.  IZRM .GE. 2
cc          IF ( IZRM .GE. 2
c     *                           .AND.  INRM .GE. 2 ) THEN
c            ITYP(JC) = 402          !alpha
c            NNUC     = NNUC - 4
c            ITYPRM   = ITYPRM - 402
c          ELSE
c            IF ( IZRM .EQ. 1 .AND. INRM .GT. IZRM ) THEN
c              ITYP(JC) = 100
c              ITYPRM   = ITYPRM - ITYP(JC)              
c            ELSEIF ( INRM .EQ. 1 .AND. IZRM .GT. INRM ) THEN
c              ITYP(JC) = 101
c              ITYPRM   = ITYPRM - ITYP(JC)         
c            ELSEIF(RD(2)*(IZRM+INRM).LT.IZRM.AND.IZRM.GE.INRM)THEN
c              ITYP(JC) = 101
c              ITYPRM   = ITYPRM - ITYP(JC)
c            ELSE
c              ITYP(JC) = 100
c              ITYPRM   = ITYPRM - ITYP(JC)
c            ENDIF
c            NNUC = NNUC - 1
c          ENDIF
c        ENDDO
c      ENDIF
c
c      IF ( NNUC .LT. 4 ) THEN
cC  EVAPORATION WITHOUT FORMATION OF ALPHA PARTICLES
c        CALL RMMARD( RD,NNUC,lseq )
c        DO  IS = 1, NNUC
c          IARM = ITYPRM/100
c          IF ( IARM .LE. 0 ) GOTO 100
c          IZRM = MOD(ITYPRM,100)
c          INRM = IARM - IZRM
c          JC   = JC + 1
c          IF ( IZRM .EQ. 1 .AND. INRM .GT. IZRM ) THEN
c            ITYP(JC) = 100
c            ITYPRM   = ITYPRM - ITYP(JC)              
c          ELSEIF ( INRM .EQ. 1 .AND. IZRM .GT. IZRM ) THEN
c            ITYP(JC) = 101
c            ITYPRM   = ITYPRM - ITYP(JC)             
c          ELSEIF ( RD(IS)*IARM .LT. IZRM .AND. IZRM .GE. INRM ) THEN
c            ITYP(JC) = 101
c            ITYPRM   = ITYPRM - ITYP(JC)
c          ELSE
c            ITYP(JC) = 100
c            ITYPRM   = ITYPRM - ITYP(JC)
c          ENDIF
c        ENDDO
      ENDIF

 100   CONTINUE

      IF(ITYPRM.GT.0)THEN
        JC = JC +1
        ITYP(JC) = ITYPRM
      ENDIF

      JFIN=JC

      ENDIF

 50   CONTINUE

      JC=JFIN
      
      DO MF=1,JFIN
      IARM = ITYP(MF)/100
      IF ( IARM .EQ. 0 )call utstop("IARM shouldn't be 0 here !&")
      IZRM = MOD(ITYP(MF),100)
      INRM = IARM - IZRM
      JC = JC + 1
C CLEAN FRAGMENT TO HAVE ALWAYS LESS N THAN 5+1.6*P  (NUCLEI TOO UNSTABLE ABOVE THIS)
c     CALL RMMARD( RD,2,lseq )
      ZMIN = MAX(0D0,MIN(DBLE(IZRM),5D0))
      DO WHILE ( INRM .GT. INT( ZMIN + 1.6D0 * DBLE(IZRM) ) )
        ITYP(JC) = 1220
        ITYP(MF)   = ITYP(MF) - 100              
        IARM = ITYP(MF)/100
        IZRM = MOD(ITYP(MF),100)
        INRM = IARM - IZRM
        JC = JC + 1
      ENDDO
      IF ( IZRM .GE. 35 ) THEN
C CLEAN FRAGMENT NOT TO HAVE MORE P THAN N WHEN Z > 35
        DO WHILE ( IZRM .GT. INRM )
          ITYP(JC) = 1120
          ITYP(MF)   = ITYP(MF) - 101              
          IARM = ITYP(MF)/100
          IZRM = MOD(ITYP(MF),100)
          INRM = IARM - IZRM
          JC = JC + 1
        ENDDO
      ELSE
C CLEAN FRAGMENT NOT TO HAVE MORE P THAN  SQRT(35*N) WHEN Z < 35
        DO WHILE ( IZRM .GT. NINT( ( SQRT( DBLE(35*INRM) ) ) )) 
          ITYP(JC) = 1120
          ITYP(MF)   = ITYP(MF) - 101              
          IARM = ITYP(MF)/100
          IZRM = MOD(ITYP(MF),100)
          INRM = IARM - IZRM
          JC = JC + 1
        ENDDO
      ENDIF
      IF ( INRM .EQ. 4 .AND. IZRM.EQ.5) THEN     !EXCLUDED
        ITYP(JC) = 1120
        ITYP(MF)   = ITYP(MF) - 101              
        IARM = ITYP(MF)/100
        IZRM = MOD(ITYP(MF),100)
        INRM = IARM - IZRM
        JC = JC + 1
      ENDIF
      IF ( (INRM .EQ. 3 .OR. INRM .EQ.5) .AND. IZRM.EQ.2) THEN     !EXCLUDED
        ITYP(JC) = 1220
        ITYP(MF)   = ITYP(MF) - 100              
        IARM = ITYP(MF)/100
        IZRM = MOD(ITYP(MF),100)
        INRM = IARM - IZRM
        JC = JC + 1
      ENDIF
      IF ( INRM .EQ. 4 .AND. IZRM.EQ.4 ) THEN     !EXCLUDED
        DO WHILE ( IZRM .GE. 2 .AND. INRM .GE. 2 )
          ITYP(JC) = -402
          ITYP(MF)   = ITYP(MF) - 402              
          IARM = ITYP(MF)/100
          IZRM = MOD(ITYP(MF),100)
          INRM = IARM - IZRM
          JC = JC + 1
        ENDDO
        IF(ITYP(MF).NE.0)call utstop("ITYP should be 0 here !&")
        JC = JC - 1
      ELSEIF     ( ITYP(MF) .GT. 200 ) THEN
        ITYP(JC) = -ITYP(MF)
      ELSEIF ( ITYP(MF) .EQ. 101 ) THEN
        ITYP(JC) = 1120
      ELSEIF ( ITYP(MF) .EQ. 100 ) THEN
        ITYP(JC) = 1220
      ELSE
        JC = JC - 1
        IF ( ITYP(MF) .NE. 0 )THEN
          WRITE(ifch,*)'EPOVAPOR : ILLEGAL PARTICLE ITYP(MF) =',ITYP(MF)
          call utstop("ITYP didn't reach 0 !&")
        ENDIF
      ENDIF
      

      ENDDO

      NFIN = JC-JFIN
      
 25   CONTINUE
      NNUC = 0
      IF ( NFIN .GT. 1 ) THEN
       FAC = MIN(DBLE(PNULLX),10D0)
       IF(ish.ge.7)WRITE(ifch,*) 
     *   'EPOVAPOR :  NO        ITYP     PFR       PFL'
       IF ( infragm .EQ. 2 ) THEN
C  EVAPORATION WITH PT AFTER PARAMETRIZED JACEE DATA
         DO  MF = 1, NFIN
          ITYP(MF) = ITYP(MF+JFIN)
          IF(ITYP(MF).LT.0)THEN
            IARM=-ITYP(MF)/100
            PFR(MF) = 0d0
            PFRZ(MF)= 0d0
          ELSE
            NNUC = NNUC + 1 
            IARM=1
            PFR(MF) = RANNORM(0.0888D0,0.0444D0) !RANNORM(0.088D0,0.044D0)
            PFRZ(MF)= RANNORM(0.D0,0.450D0*FAC) !Fermi motion about 300 MeV
          ENDIF
          IF(ish.ge.7)WRITE(ifch,*) MF,ITYP(MF),SNGL(PFR(MF))
     &                                         ,SNGL(PFRZ(MF))
        ENDDO
       ELSEIF ( infragm .EQ. 3 ) THEN
C  EVAPORATION WITH PT AFTER GOLDHABER''S MODEL (PHYS.LETT.53B(1974)306)
        DO  MF = 1, NFIN
          ITYP(MF) = ITYP(MF+JFIN)
          K    = MAX( 1, -ITYP(MF)/100 )
          BGLH = K * (MAPRO - K) / DBLE(MAPRO-1)
C  THE VALUE 0.103 [GEV] IS SIGMA(0)=P(FERMI)/SQRT(5.)
*         AGLH = 0.103D0 * SQRT( BGLH )
C  THE VALUE 0.090 [GEV] IS EXPERIMENTALLY DETERMINED SIGMA(0)
          AGLH = 0.090D0 * SQRT( BGLH )
          PFR(MF) = RANNORM(0.D0,AGLH)
          PFRZ(MF)= RANNORM(0.000D0,0.500D0)    !from pAg at 100 GeV
          IF(ish.ge.7)WRITE(ifch,*) MF,ITYP(MF),SNGL(PFR(MF))
     &                                         ,SNGL(PFRZ(MF))
        ENDDO
      ELSE
C  EVAPORATION WITHOUT TRANSVERSE MOMENTUM
        DO  MF = 1, NFIN
          ITYP(MF) = ITYP(MF+JFIN)
          PFR(MF) = 0.D0
          PFRZ(MF)= 0.D0
          IF(ish.ge.7)WRITE(ifch,*) MF,ITYP(MF),SNGL(PFR(MF))
     &                                         ,SNGL(PFRZ(MF))
        ENDDO
       ENDIF
       JFIN=NFIN
C  CALCULATE RESIDUAL TRANSVERSE MOMENTUM
       SPFRX = 0.D0
       SPFRY = 0.D0
       SPFRZ = 0.D0
       CALL RMMARD( RD,JFIN,lseq )
       DO  MF = 1, JFIN
        PHIFR = PI * RD(MF)
        PFRX(MF) = PFR(MF) * COS( PHIFR )
        PFRY(MF) = PFR(MF) * SIN( PHIFR )
        SPFRY = SPFRY + PFRY(MF)
        SPFRX = SPFRX + PFRX(MF)
        SPFRZ = SPFRZ + PFRZ(MF)
       ENDDO
C  CORRECT ALL MOMENTA FOR MOMENTUM CONSERVATION
       SPFRX = SPFRX / NNUC !JFIN
       SPFRY = SPFRY / NNUC !JFIN
       SPFRZ = SPFRZ / NNUC !JFIN
       DO  MF = 1, JFIN
         IF(ITYP(MF).GT.0D0)THEN
           PFRX(MF) = PFRX(MF) - SPFRX
           PFRY(MF) = PFRY(MF) - SPFRY
           PFRZ(MF) = PFRZ(MF) - SPFRZ
         ENDIF
       ENDDO
      ELSE
C ONLY ONE FRAGMENT
        JFIN=1
        MF = 1
        ITYP(MF) = ITYP(MF+JFIN)
        PFR(MF) = 0.D0
        PFRZ(MF)= 0.D0
        IF(ish.ge.7)WRITE(ifch,*) MF,ITYP(MF),SNGL(PFR(MF))
     &                                         ,SNGL(PFRZ(MF))
      ENDIF

      IF(ish.ge.7)WRITE(ifch,*) 'EPOVAPOR : NINTA,JFIN=',NINTA,JFIN

      RETURN
      END

*-- Author :    The CORSIKA development group   21/04/1994
C=======================================================================

      DOUBLE PRECISION FUNCTION RANNORM( A,B )

C-----------------------------------------------------------------------
C  RAN(DOM NUMBER) NOR(MALLY DISTRIBUTED)
C
C  GENERATES NORMAL DISTRIBUTED RANDOM NUMBER
C  DELIVERS 2 UNCORRELATED RANDOM NUMBERS,
C  THEREFORE RANDOM CALLS ARE ONLY NECESSARY EVERY SECOND TIME.
c  but to be used with CONEX/CORSIKA we should always use 2 new number
c  to be able to reproduce the shower
C  REFERENCE : NUMERICAL RECIPES, W.H. PRESS ET AL.,
C              CAMBRIDGE UNIVERSITY PRESS, 1992  ISBN 0 521 43064 X
C  ARGUMENTS:
C   A      = MEAN VALUE
C   B      = STANDARD DEVIATION
C
C  FROM CORSIKA
C-----------------------------------------------------------------------

      IMPLICIT NONE
      double precision facrdm,u1rdm,u2rdm,drangen
c      logical knordm
c      data knordm/.true./

      DOUBLE PRECISION A,B,RR
      SAVE facrdm,u1rdm,u2rdm!,knordm
C-----------------------------------------------------------------------

c      IF ( KNORdm ) THEN
  1     CONTINUE
        U1rdm = 2.D0*drangen(a) - 1.D0
        U2rdm = 2.D0*drangen(b) - 1.D0
        RR = U1rdm**2 + U2rdm**2
        IF ( RR .GE. 1.D0  .OR.  RR .EQ. 0.D0 ) GOTO 1
        FACrdm = SQRT( (-2.D0) * LOG(RR) / RR )

        RANNORM = FACrdm * U1rdm * B + A
c        KNORdm   = .FALSE.
c      ELSE
c        RANNORM = FACrdm * U2rdm * B + A
c        KNORdm   = .TRUE.
c      ENDIF

      RETURN
      END

c-----------------------------------------------------------------------
      function ransig()
c-----------------------------------------------------------------------
c     returns randomly +1 or -1
c-----------------------------------------------------------------------
      ransig=1
      if(rangen().gt.0.5)ransig=-1
      return
      end

cc-----------------------------------------------------------------------
c      function ranxq(n,x,q,xmin)
cc-----------------------------------------------------------------------
cc     returns random number according to x(i) q(i) with x>=xmin
cc-----------------------------------------------------------------------
c      include 'epos.inc'
c      real x(n),q(n)
c      imin=1
c      if(xmin.eq.0.)goto3
c      i1=1
c      i2=n
c1     i=i1+(i2-i1)/2
c      if(x(i).lt.xmin)then
c      i1=i
c      elseif(x(i).gt.xmin)then
c      i2=i
c      else
c      imin=i
c      goto3
c      endif
c      if(i2-i1.gt.1)goto1
c      imin=i2
c3     continue
c      if(q(imin).gt.q(n)*.9999)then
c      ranxq=xmin
c      goto4
c      endif
c      qran=q(imin)+rangen()*(q(n)-q(imin))
c      ranxq=utinvt(n,x,q,qran)
c4     continue
c
c      if(ranxq.lt.xmin)then
c      call utmsg('ranxq&')
c      write(ifch,*)'*****  ranxq=',ranxq,' <       xmin=',xmin
c      write(ifch,*)'q(imin) q q(n):',q(imin),qran,q(n)
c      write(ifch,*)'x(imin) x x(n):',x(imin),ranxq,x(n)
c      call utmsgf
c      ranxq=xmin
c      endif
c
c      return
c      end
c
cc  ***** end r-routines
cc  ***** beg s-routines
c
cc-----------------------------------------------------------------------
c      function sbet(z,w)
cc-----------------------------------------------------------------------
c      sbet=utgam1(z)*utgam1(w)/utgam1(z+w)
c      return
c      end
c
cc-----------------------------------------------------------------------
c      function smass(a,y,z)
cc-----------------------------------------------------------------------
cc     returns droplet mass (in gev) (per droplet, not (!) per nucleon)
cc     according to berger/jaffe mass formula, prc35(1987)213 eq.2.31,
cc     see also c. dover, BNL-46322, intersections-meeting, tucson, 91.
cc     a: massnr, y: hypercharge, z: charge,
cc-----------------------------------------------------------------------
c      common/cmass/thet,epsi,as,ac,dy,dz,ym,cz,zm,sigma,rzero
c      ymin=ym*a
c      zmin=cz/(dz/a+zm/a**(1./3.))
c      smass=epsi*a+as*a**(2./3.)+(ac/a**(1./3.)+dz/a/2.)*(z-zmin)**2
c     *+dy/a/2.*(y-ymin)**2
c      return
c      end
c
cc-----------------------------------------------------------------------
c      subroutine smassi(theta)
cc-----------------------------------------------------------------------
cc     initialization for smass.
cc     calculates parameters for berger/jaffe mass formula
cc     (prc35(1987)213 eq.2.31, see also c. dover, BNL-46322).
cc     theta: parameter that determines all parameters in mass formula.
cc-----------------------------------------------------------------------
c      common/cmass/thet,epsi,as,ac,dy,dz,ym,cz,zm,sigma,rzero
c      thet=theta
c
c      astr=.150
c      pi=3.14159
c      alp=1./137.
c
c      co=cos(theta)
c      si=sin(theta)
c      bet=(1+co**3)/2.
c      rzero=si/astr/(  2./3./pi*(1+co**3)  )**(1./3.)
cctp060829      cs=astr/si
c      cz=-astr/si*(  (  .5*(1+co**3)  )**(1./3.)-1  )
c      sigma=6./8./pi*(astr/si)**3*(co**2/6.-si**2*(1-si)/3.-
c     *1./3./pi*(pi/2.-theta-sin(2*theta)+si**3*alog((1+co)/si)))
c
c      epsi=astr*((.5*(1+co**3))**(1./3.)+2)/si
c      as=4*pi*sigma*rzero**2
c      ac=3./5.*alp/rzero
c      dz=astr/si*bet**(1./3.)*co**2*
c     *(co**4*(1+bet**(2./3.))+(1+bet)**2)/
c     *(  (2*co**2+bet**(1./3.))*(co**4*(1+bet**(2./3.))+(1+bet)**2)-
c     *(co**4+bet**(1./3.)*(1+bet))*((2*bet**(2./3.)-1)*co**2+1+bet)  )
c      dy=astr/6.*(1+co**3)**3/si*
c     *(  1+(1+co)/(4*(1+co**3))**(2./3.)  )/
c     *(co**6+co+co*(.5*(1+co**3))**(4./3.))
c      zm=6*alp/(5*rzero)
c      ym=(1-co**3)/(1+co**3)
c
c      return
c      end
c
cc-----------------------------------------------------------------------
c      subroutine smassp
cc-----------------------------------------------------------------------
cc     prints smass.
cc-----------------------------------------------------------------------
c      include 'epos.inc'
c      common/cmass/thet,epsi,as,ac,dy,dz,ym,cz,zm,sigma,rzero
c      real eng(14),ymi(14),zmi(14)
c      pi=3.14159
c      write(ifch,*)'parameters of mass formula:'
c      write(ifch,*)'theta=',thet,'   epsi=',epsi
c      write(ifch,*)'as=',as,'   ac=',ac
c      write(ifch,*)'dy=',dy,'   dz=',dz
c      write(ifch,*)'ym=',ym
c      write(ifch,*)'cz dz zm=',cz,dz,zm
c      write(ifch,*)'sigma**1/3=',sigma**(1./3.),'   rzero=',rzero
c      write(ifch,*)'mass:'
c      write(ifch,5000)(j,j=1,14)
c5000  format(5x,'a:',14i5)
c      do 4 j=1,14
c      a=j
c      ymi(j)=ym*a
c4     zmi(j)=cz/(dz/a+zm/a**(1./3.))
c      write(ifch,5002)(ymi(j),j=1,14)
c5002  format(1x,'ymin: ',14f5.2)
c      write(ifch,5003)(zmi(j),j=1,14)
c5003  format(1x,'zmin: ',14f5.2)
c      do 2 i=1,15
c      ns=11-i
c      do 3 j=1,14
c      a=j
c      y=a-ns
c      z=0.
c3     eng(j)=smass(a,y,z)/a
c      write(ifch,5001)ns,(eng(j),j=1,14)
c5001  format(1x,'s=',i2,2x,14f5.2)
c2     continue
c      write(ifch,*)'mass-mass(free):'
c      write(ifch,5000)(j,j=1,14)
c      do 5 i=1,15
c      ns=11-i
c      do 6 j=1,14
c      a=j
c      y=a-ns
c      z=0.
c      call smassu(a,y,z,ku,kd,ks,kc)
c6     eng(j)=(smass(a,y,z)-utamnu(ku,kd,ks,kc,0,0,3))/a
c      write(ifch,5001)ns,(eng(j),j=1,14)
c5     continue
c
c      stop
c      end
c
cc-----------------------------------------------------------------------
c      subroutine smasst(kux,kdx,ksx,kcx,a,y,z)
cc-----------------------------------------------------------------------
cc     input: kux,kdx,ksx,kcx = net quark numbers (for u,d,s,c quarks).
cc     output: massnr a, hypercharge y and charge z.
cc-----------------------------------------------------------------------
c      sg=1
c      if(kux+kdx+ksx+kcx.lt.0.)sg=-1
c      ku=sg*kux
c      kd=sg*kdx
c      ks=sg*ksx
c      kc=sg*kcx
c      k=ku+kd+ks+kc
c      if(mod(k,3).ne.0)stop'noninteger baryon number'
c      a=k/3
c      y=a-ks
c      nz=2*ku-kd-ks+2*kc
c      if(mod(nz,3).ne.0)stop'noninteger charge'
c      z=nz/3
c      return
c      end
c
cc-----------------------------------------------------------------------
c      subroutine smassu(ax,yx,zx,ku,kd,ks,kc)
cc-----------------------------------------------------------------------
cc     input: massnr ax, hypercharge yx and charge zx.
cc     output: ku,kd,ks,kc = net quark numbers (for u,d,s,c quarks).
cc-----------------------------------------------------------------------
c      sg=1
c      if(ax.lt.0.)sg=-1
c      a=sg*ax
c      y=sg*yx
c      z=sg*zx
c      ku=nint(a+z)
c      kd=nint(a-z+y)
c      ks=nint(a-y)
c      kc=0
c      return
c      end
c
cc-----------------------------------------------------------------------
c      function spoc(a,b,c,d,x)
cc-----------------------------------------------------------------------
cc     power fctn with cutoff
cc-----------------------------------------------------------------------
c      spoc=0
c      if(a.eq.0..and.b.eq.0.)return
c      spoc =a+b*x**c
c      spoc0=a+b*d**c
c      spoc=amin1(spoc,spoc0)
c      spoc=amax1(0.,spoc)
c      return
c      end
c
c-----------------------------------------------------------------------
      function utacos(x)
c-----------------------------------------------------------------------
c     returns acos(x) for -1 <= x <= 1 , acos(+-1) else
c-----------------------------------------------------------------------
      include 'epos.inc'
      argum=x
      if(x.lt.-1.)then
      if(ish.ge.1)then
      call utmsg('utacos&')
      write(ifch,*)'*****  argum = ',argum,' set -1'
      call utmsgf
      endif
      argum=-1.
      elseif(x.gt.1.)then
      if(ish.ge.1)then
      call utmsg('utacos&')
      write(ifch,*)'*****  argum = ',argum,' set 1'
      call utmsgf
      endif
      argum=1.
      endif
      utacos=acos(argum)
      return
      end

c----------------------------------------------------------------------
      function utamnu(keux,kedx,kesx,kecx,kebx,ketx,modus)
c----------------------------------------------------------------------
c     returns min mass of droplet with given u,d,s,c content
c     keux: net u quark number
c     kedx: net d quark number
c     kesx: net s quark number
c     kecx: net c quark number
c     kebx: net b quark number
c     ketx: net t quark number
c     modus: 4=two lowest multiplets; 5=lowest multiplet
c----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/csjcga/amnull,asuha(7)
      common/drop4/asuhax(7),asuhay(7)

      if(modus.lt.4.or.modus.gt.5)stop'UTAMNU: not supported'
c 1    format(' flavours:',6i5 )
c 100  format(' flavours+mass:',6i5,f8.2 )
c      write(ifch,1)keux,kedx,kesx,kecx,kebx,ketx

      amnull=0.

      do i=1,7
      if(modus.eq.4)asuha(i)=asuhax(i)    !two lowest multiplets
      if(modus.eq.5)asuha(i)=asuhay(i)    !lowest multiplet
      enddo

      ke=iabs(keux+kedx+kesx+kecx+kebx+ketx)

      if(keux+kedx+kesx+kecx+kebx+ketx.ge.0)then
      keu=keux
      ked=kedx
      kes=kesx
      kec=kecx
      keb=kebx
      ket=ketx
      else
      keu=-keux
      ked=-kedx
      kes=-kesx
      kec=-kecx
      keb=-kebx
      ket=-ketx
      endif

c      write(ifch,*)keu,ked,kes,kec,keb,ket

c   removing top mesons  to remove t quarks or antiquarks
      if(ket.ne.0)then
12    continue
      ii=sign(1,ket)
      ket=ket-ii
      if(ii*keu.le.ii*ked)then
      keu=keu+ii
      else
      ked=ked+ii
      endif
      amnull=amnull+200.    ! ???????
      if(ket.ne.0)goto12
      endif

c   removing bottom mesons  to remove b quarks or antiquarks
      if(keb.ne.0)then
11    continue
      ii=sign(1,keb)
      keb=keb-ii
      if(ii*keu.le.ii*ked)then
      keu=keu+ii
      else
      ked=ked+ii
      endif
      amnull=amnull+6. !5.28   ! (more than B-meson)
      if(keb.ne.0)goto11
      endif

c   removing charm mesons  to remove c quarks or antiquarks
      if(kec.ne.0)then
10    continue
      ii=sign(1,kec)
      kec=kec-ii
      if(keu*ii.le.ked*ii)then
      keu=keu+ii
      else
      ked=ked+ii
      endif
      amnull=amnull+2.2 !1.87  ! (more than D-meson)
      if(kec.ne.0)goto10
      endif

c      write(ifch,100)keu,ked,kes,kec,keb,ket,amnull

c   removing mesons to remove s antiquarks
5     continue
      if(kes.lt.0)then
      amnull=amnull+asuha(6)
      if(keu.ge.ked)then
      keu=keu-1
      else
      ked=ked-1
      endif
      kes=kes+1
      goto5
      endif

c   removing mesons to remove d antiquarks
6     continue
      if(ked.lt.0)then
      if(keu.ge.kes)then
      amnull=amnull+asuha(5)
      keu=keu-1
      else
      amnull=amnull+asuha(6)
      kes=kes-1
      endif
      ked=ked+1
      goto6
      endif

c   removing mesons to remove u antiquarks
7     continue
      if(keu.lt.0)then
      if(ked.ge.kes)then
      amnull=amnull+asuha(5)
      ked=ked-1
      else
      amnull=amnull+asuha(6)
      kes=kes-1
      endif
      keu=keu+1
      goto7
      endif

c      write(ifch,100)keu,ked,kes,kec,keb,ket,amnull
c      print*,keu,ked,kes,kec,keb,ket,amnull

      if(keu+ked+kes+kec+keb+ket.ne.ke)
     *call utstop('utamnu: sum_kei /= ke&')
      keq=keu+ked
      keqx=keq
      amnux=0

c   removing strange baryons
      i=4
2     i=i-1
3     continue
      if((4-i)*kes.gt.(i-1)*keq)then
      amnux=amnux+asuha(1+i)
      kes=kes-i
      keq=keq-3+i
      if(kes.lt.0)call utstop('utamnu: negative kes&')
      if(keq.lt.0)call utstop('utamnu: negative keq&')
      goto3
      endif
      if(i.gt.1)goto2
      if(keqx.gt.keq)then
      do 8 k=1,keqx-keq
      if(keu.ge.ked)then
      keu=keu-1
      else
      ked=ked-1
      endif
8     continue
      endif

      if(keu+ked.ne.keq)call utstop('utamnu: keu+ked /= keq&')
c      write(ifch,100)keu,ked,kes,kec,keb,ket,amnull+amnux
c      print*,keu,ked,kes,kec,keb,ket,amnull+amnux

c   removing nonstrange baryons
9     continue
      if(keu.gt.2*ked)then
      amnux=amnux+asuha(7)
      keu=keu-3
      if(keu.lt.0)call utstop('utamnu: negative keu&')
      goto9
      endif
      if(ked.gt.2*keu)then
      amnux=amnux+asuha(7)
      ked=ked-3
      if(ked.lt.0)call utstop('utamnu: negative ked&')
      goto9
      endif
      keq=keu+ked

c      write(ifch,100)keu,ked,kes,kec,keb,ket,amnull+amnux
c      print*,keu,ked,kes,kec,keb,ket,amnull+amnux

      if(mod(keq,3).ne.0)call utstop('utamnu: mod(keq,3) /= 0&')
      amnux=amnux+asuha(1)*keq/3

c      write(ifch,100)keu,ked,kes,kec,keb,ket,amnull+amnux
c      print*,keu,ked,kes,kec,keb,ket,amnull+amnux

      amnull=amnull+amnux

      if(amnull.eq.0)amnull=asuha(5)

      utamnu=amnull
      return
      end

c-----------------------------------------------------------------------
      function utamnx(jcp,jcm)
c-----------------------------------------------------------------------
c returns minimum mass for the decay of jcp---jcm (by calling utamnu).
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jcp(nflav,2),jcm(nflav,2)

      do i=1,nflav
      do j=1,2
      if(jcp(i,j).ne.0)goto1
      enddo
      enddo
      keu=jcm(1,1)-jcm(1,2)
      ked=jcm(2,1)-jcm(2,2)
      kes=jcm(3,1)-jcm(3,2)
      kec=jcm(4,1)-jcm(4,2)
      keb=jcm(5,1)-jcm(5,2)
      ket=jcm(6,1)-jcm(6,2)
      utamnx=utamnu(keu,ked,kes,kec,keb,ket,5)
      return
1     continue

      do i=1,nflav
      do j=1,2
      if(jcm(i,j).ne.0)goto2
      enddo
      enddo
      keu=jcp(1,1)-jcp(1,2)
      ked=jcp(2,1)-jcp(2,2)
      kes=jcp(3,1)-jcp(3,2)
      kec=jcp(4,1)-jcp(4,2)
      keb=jcp(5,1)-jcp(5,2)
      ket=jcp(6,1)-jcp(6,2)
      utamnx=utamnu(keu,ked,kes,kec,keb,ket,5)
      return
2     continue

      keu=jcp(1,1)-jcp(1,2)
      ked=jcp(2,1)-jcp(2,2)
      kes=jcp(3,1)-jcp(3,2)
      kec=jcp(4,1)-jcp(4,2)
      keb=jcp(5,1)-jcp(5,2)
      ket=jcp(6,1)-jcp(6,2)
      ke=keu+ked+kes+kec+keb+ket
      if(mod(ke+1,3).eq.0)then
        keu=keu+1
        amms1=utamnu(keu,ked,kes,kec,keb,ket,5)
        keu=keu-1
        ked=ked+1
        amms2=utamnu(keu,ked,kes,kec,keb,ket,5)
      elseif(mod(ke-1,3).eq.0)then
        keu=keu-1
        amms1=utamnu(keu,ked,kes,kec,keb,ket,5)
        keu=keu+1
        ked=ked-1
        amms2=utamnu(keu,ked,kes,kec,keb,ket,5)
      else
        amms1=0
        amms2=0
        amms3=0
        amms4=0
        call utstop('utamnx: no singlet possible (1)&')
      endif
      keu=jcm(1,1)-jcm(1,2)
      ked=jcm(2,1)-jcm(2,2)
      kes=jcm(3,1)-jcm(3,2)
      kec=jcm(4,1)-jcm(4,2)
      keb=jcm(5,1)-jcm(5,2)
      ket=jcm(6,1)-jcm(6,2)
      ke=keu+ked+kes+kec+keb+ket
      if(mod(ke+1,3).eq.0)then
        keu=keu+1
        amms3=utamnu(keu,ked,kes,kec,keb,ket,5)
        keu=keu-1
        ked=ked+1
        amms4=utamnu(keu,ked,kes,kec,keb,ket,5)
      elseif(mod(ke-1,3).eq.0)then
        keu=keu-1
        amms3=utamnu(keu,ked,kes,kec,keb,ket,5)
        keu=keu+1
        ked=ked-1
        amms4=utamnu(keu,ked,kes,kec,keb,ket,5)
      else
        call utstop('utamnx: no singlet possible (2)&')
      endif
      utamnx=min(amms1+amms3,amms2+amms4)
c       print *,amms1,amms3,amms2,amms4,jcp,jcm
      return
      end



cc-----------------------------------------------------------------------
c      function utamny(jcp,jcm)
cc-----------------------------------------------------------------------
cc returns minimum mass of jcp+jcm (by calling utamnu).
cc-----------------------------------------------------------------------
c      parameter (nflav=6)
c      integer jcp(nflav,2),jcm(nflav,2),jc(nflav,2)
c      do 7 nf=1,nflav
c      jc(nf,1)=jcp(nf,1)+jcm(nf,1)
c7     jc(nf,2)=jcp(nf,2)+jcm(nf,2)
c      keu=jc(1,1)-jc(1,2)
c      ked=jc(2,1)-jc(2,2)
c      kes=jc(3,1)-jc(3,2)
c      kec=jc(4,1)-jc(4,2)
c      keb=jc(5,1)-jc(5,2)
c      ket=jc(6,1)-jc(6,2)
c      utamny=utamnu(keu,ked,kes,kec,keb,ket,5)
c      return
c      end
c
c-----------------------------------------------------------------------
      function utamnz(jc,modus)
c-----------------------------------------------------------------------
c returns minimum mass of jc (by calling utamnu).
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2)
      keu=jc(1,1)-jc(1,2)
      ked=jc(2,1)-jc(2,2)
      kes=jc(3,1)-jc(3,2)
      kec=jc(4,1)-jc(4,2)
      keb=jc(5,1)-jc(5,2)
      ket=jc(6,1)-jc(6,2)
      utamnz=utamnu(keu,ked,kes,kec,keb,ket,modus)
      return
      end

c-----------------------------------------------------------------------
      subroutine utar(i1,i2,i3,x0,x1,x2,x3,xx)
c-----------------------------------------------------------------------
c     returns the array xx with xx(1)=x0 <= xx(i) <= xx(i3)=x3
c-----------------------------------------------------------------------
      real xx(i3)
      do 1 i=1,i1-1
  1   xx(i)=x0+(i-1.)/(i1-1.)*(x1-x0)
      do 2 i=i1,i2-1
  2   xx(i)=x1+(i-i1*1.)/(i2-i1*1.)*(x2-x1)
      do 3 i=i2,i3
  3   xx(i)=x2+(i-i2*1.)/(i3-i2*1.)*(x3-x2)
      return
      end

cc---------------------------------------------------------------------
c      subroutine utaxis(i,j,a1,a2,a3)
cc-----------------------------------------------------------------------
cc     calculates the axis defined by the ptls i,j in the i,j cm system
cc---------------------------------------------------------------------
c      include 'epos.inc'
c      double precision pi1,pi2,pi3,pi4,pj1,pj2,pj3,pj4,p1,p2,p3,p4,p5
c     *,err,a
c      a1=0
c      a2=0
c      a3=1
c      pi1=dble(pptl(1,i))
c      pi2=dble(pptl(2,i))
c      pi3=dble(pptl(3,i))
c      pi4=dble(pptl(4,i))
c      pj1=dble(pptl(1,j))
c      pj2=dble(pptl(2,j))
c      pj3=dble(pptl(3,j))
c      pj4=dble(pptl(4,j))
c      p1=pi1+pj1
c      p2=pi2+pj2
c      p3=pi3+pj3
c      p4=pi4+pj4
c      p5=dsqrt(p4**2-p3**2-p2**2-p1**2)
c      call utlob2(1,p1,p2,p3,p4,p5,pi1,pi2,pi3,pi4,50)
c      call utlob2(1,p1,p2,p3,p4,p5,pj1,pj2,pj3,pj4,51)
c           err=(pi1+pj1)**2+(pi2+pj2)**2+(pi3+pj3)**2
c           if(err.gt.1d-3)then
c      call utmsg('utaxis&')
c      write(ifch,*)'*****  err=',err
c      write(ifch,*)'pi:',pi1,pi2,pi3,pi4
c      write(ifch,*)'pj:',pj1,pj2,pj3,pj4
c      call utmsgf
c           endif
c      a=dsqrt( (pj1-pi1)**2 + (pj2-pi2)**2 + (pj3-pi3)**2 )
c      if(a.eq.0.d0)return
c      a1=sngl((pi1-pj1)/a)
c      a2=sngl((pi2-pj2)/a)
c      a3=sngl((pi3-pj3)/a)
c      return
c      end
c
cc---------------------------------------------------------------------
c      subroutine uthalf(i,j,zz1,zz2,iret)
cc-----------------------------------------------------------------------
cc     give equal energy (E_i+E_j)/2 to particle i+j in their cm system
cc---------------------------------------------------------------------
c      include 'epos.inc'
c      double precision pi1,pi2,pi3,pi4,pi5,pj1,pj2,pj3,pj4,pj5
c     *,p1,p2,p3,p4,p5,err,pt,pti,sinp,cosp,pmax,phi,drangen!,rrr
c      iret=0
c      pi1=dble(pptl(1,i))
c      pi2=dble(pptl(2,i))
c      pi3=dble(pptl(3,i))
c      pi5=dble(pptl(5,i))
c      pi4=sqrt(pi1**2+pi2**2+pi3**2+pi5**2)
c      pj1=dble(pptl(1,j))
c      pj2=dble(pptl(2,j))
c      pj3=dble(pptl(3,j))
c      pj5=dble(pptl(5,j))
c      pj4=sqrt(pj1**2+pj2**2+pj3**2+pj5**2)
c      if(ish.ge.6)then
c        write(ifch,*)'uthalf for ',i,' and ',j
c        write(ifch,*)'in ',idptl(i),pi1,pi2,pi3,pi4,pi5
c        write(ifch,*)'<> ',idptl(j),pj1,pj2,pj3,pj4,pj5
c      endif
c      p1=pi1+pj1
c      p2=pi2+pj2
c      p3=pi3+pj3
c      p4=pi4+pj4
c      p5=(p4-p3)*(p4+p3)-p2**2-p1**2
c      if(p5.lt.0d0.or.(pi3.lt.0.99d0*pj3.and.mod(ityptl(j)/10,10).eq.4)
c     &       .or.(pi3.gt.0.99d0*pj3.and.mod(ityptl(j)/10,10).eq.5))then
cc      if(p5.lt.0d0)then
c        if(ish.ge.7)write(ifch,*)'Inversion not possible (1)',p5,pi3,pj3
c        iret=1
c        return
c      else
c        p5=sqrt(p5)
c      endif
c      call utlob2(1,p1,p2,p3,p4,p5,pi1,pi2,pi3,pi4,50)
c      call utlob2(1,p1,p2,p3,p4,p5,pj1,pj2,pj3,pj4,51)
c           err=(pi1+pj1)**2+(pi2+pj2)**2+(pi3+pj3)**2
c           if(err.gt.1d-3)then
c      call utmsg('uthalf&')
c      write(ifch,*)'*****  err=',err
c      write(ifch,*)'pi:',pi1,pi2,pi3,pi4
c      write(ifch,*)'pj:',pj1,pj2,pj3,pj4
c      call utmsgf
c           endif
c      if(ish.ge.8)then
c        write(ifch,*)'pi 1:',pi1,pi2,pi3,pi4,pi5
c     &       ,sqrt(pi1**2+pi2**2+pi3**2+pi5**2)
c        write(ifch,*)'pj 1:',pj1,pj2,pj3,pj4,pj5
c     &       ,sqrt(pj1**2+pj2**2+pj3**2+pj5**2)
c      endif
c
c      phi=drangen(p5)*2d0*dble(pi)
c      pti=sqrt(pi1*pi1+pi2*pi2)
c
cc      sinp=abs(asin(pi2/pti))
c
c      cosp=cos(phi)
c      sinp=sin(phi)
cc      cosp=pi1/pti
cc      sinp=pi2/pti
c      pini=abs(pi3)
c      pmax=pini
c
c      ntry=0
c 10   ntry=ntry+1
c      pi1=-pj1
c      pi2=-pj2
c      pi3=-pj3
c
c
cc      rrr=dble(ranptcut(4.))*max(0.1d0,dble(zz))+0.01d0
cc      rrr=dble(ranptcut(zz*max(1.,float(ntry/10)**3)))        !flip if pt too large
cc      rrr=dble(min(1.,ranptcut(zz/float(max(1,ntry/10)))))        !return if pt too large
c
c      rrr=rangen()
c      phi=dble(pi*zz1*(1.-rrr*zz2/float(max(1,ntry/10))))
c      call utroa2(phi,cosp,sinp,0d0,pi1,pi2,pi3)    !rotation around an axis perpendicular to p3
c      pt=pi1*pi1+pi2*pi2
c
ccc      pi3=pini*(1d0-(pmax/pini+1d0)*min(1d0,rrr))
ccc     &   *sign(1.d0,pi3)
cc      pi3=pini*(1d0-(pmax/pini+1d0)*(1d0-min(1d0,rrr)))
cc     &   *sign(1.d0,pi3)
cccc      pi3=(1d0-exp(-drangen(pti)))*pi3
cccc      pi3=min(1d0,(drangen()+exp(-min(100d0,dble(zz)))))*pi3
cccc      pi3=(dble(reminv)+exp(-min(100d0,dble(zz))))*pi3
cc      pt=((pi4+pi3)*(pi4-pi3)-pi5*pi5)
c      if(ish.ge.8)write(ifch,*)'ut',ntry,zz,rrr,-pi3/pj3,pt,pti*pti
c     &,idptl(i),idptl(j)
c      if((pt.lt.0d0.or.pt.gt.max(2.*pti*pti,1.d0))
cc      if(pt.lt.0d0
c     &                               .and.ntry.lt.1000)then
c        goto 10
c      elseif(ntry.ge.1000)then
c        if(ish.ge.7)write(ifch,*)'Inversion not possible (2)',pt,pti*pti
c        iret=1
c        return         !pion distribution fall at xf=1
cc        pi3=pj3       !if flip all particle with large pt
cc        pt=pti*pti    !then very very hard pi+ spectra (flat at xf=1 like in pp400 !)
c                       !but to hard for K, rho, eta, etc ..
c      endif
cc      print *,'ut',ntry,phi/(pi-sinp),zz,rrr,-pi3/pj3,pt/pti/pti
cc      pt=sqrt(pt)
cc      pi1=cosp*pt
cc      pi2=sinp*pt
c
c      pj1=-pi1
c      pj2=-pi2
c      pj3=-pi3
cc      pi4=sqrt(pi3*pi3+pt*pt+pi5*pi5)
cc      pj4=sqrt(pj3*pj3+pt*pt+pj5*pj5)
c      pi4=sqrt(pi3*pi3+pt+pi5*pi5)
c      pj4=sqrt(pj3*pj3+pt+pj5*pj5)
c      if(ish.ge.8)then
c        write(ifch,*)'pi 2:',pi1,pi2,pi3,pi4,pi5
c     &       ,sqrt(pi1**2+pi2**2+pi3**2+pi5**2)
c        write(ifch,*)'pj 2:',pj1,pj2,pj3,pj4,pj5
c     &       ,sqrt(pj1**2+pj2**2+pj3**2+pj5**2)
c      endif
c      call utlob2(-1,p1,p2,p3,p4,p5,pi1,pi2,pi3,pi4,-50)
c      call utlob2(-1,p1,p2,p3,p4,p5,pj1,pj2,pj3,pj4,-51)
c      if(pi3/dble(pptl(3,i)).gt.1.00001d0
c     &.or.(pi3.gt.1.001d0*pj3.and.mod(ityptl(j)/10,10).eq.4)
c     &.or.(pi3.lt.1.001d0*pj3.and.mod(ityptl(j)/10,10).eq.5))then
c        if(ish.ge.7)write(ifch,*)'Inversion not possible (3)',pi3,pj3
c        iret=1
c        return
c      endif
c      id=idptl(i)
c      pptl(1,i)=sngl(pi1)
c      pptl(2,i)=sngl(pi2)
c      pptl(3,i)=sngl(pi3)
c      pptl(4,i)=sngl(pi4)
c      pptl(5,i)=sngl(pi5)
cc      pptl(5,i)=sngl(pj5)
cc      idptl(i)=idptl(j)
c      pptl(1,j)=sngl(pj1)
c      pptl(2,j)=sngl(pj2)
c      pptl(3,j)=sngl(pj3)
c      pptl(4,j)=sngl(pj4)
c      pptl(5,j)=sngl(pj5)
cc      pptl(5,j)=sngl(pi5)
cc      idptl(j)=id
c      if(ish.ge.6)then
c        write(ifch,*)'out ',idptl(i),pi1,pi2,pi3,pi4
c        write(ifch,*)' <> ',idptl(j),pj1,pj2,pj3,pj4
c      endif
c
c      return
c      end
c
cc-----------------------------------------------------------------------
c      subroutine utchm(arp,arm,ii)
cc-----------------------------------------------------------------------
cc     checks whether arp**2=0 and arm**2=0.
cc-----------------------------------------------------------------------
c      include 'epos.inc'
c      double precision arp(4),arm(4),difp,difm
c      difp=arp(4)**2-arp(1)**2-arp(2)**2-arp(3)**2
c      difm=arm(4)**2-arm(1)**2-arm(2)**2-arm(3)**2
c      if(dabs(difp).gt.1e-3*arp(4)**2
c     *.or.dabs(difm).gt.1e-3*arm(4)**2)then
c      call utmsg('utchm&')
c      write(ifch,*)'*****  mass non zero  -  ',ii
c      write(ifch,*)'jet-mass**2`s:    ',difp,difm
c      write(ifch,*)'energy**2`s:      ',arp(4)**2,arm(4)**2
c      write(ifch,*)(sngl(arp(i)),i=1,4)
c      write(ifch,*)(sngl(arm(i)),i=1,4)
c      call utmsgf
c      endif
c      return
c      end
c
c-----------------------------------------------------------------------
      subroutine utclea(nptlii,nptl0)
c-----------------------------------------------------------------------
c     starting from nptlii
c     overwrites istptl=99 particles in /cptl/, reduces so nptl
c     and update minfra and maxfra
c-----------------------------------------------------------------------
      include 'epos.inc'
      integer newptl(mxptl)!,oldptl(mxptl),ii(mxptl)

      ish0=ish
      if(ishsub/100.eq.18)ish=mod(ishsub,100)

      call utpri('utclea',ish,ishini,2)

      nptli=max(maproj+matarg+1,nptlii)
      minfra0=minfra
      maxfra0=maxfra
      minfra1=maxfra
      maxfra1=minfra
      if(ish.ge.2)write(ifch,*)'entering subr utclea:',nptl
     &                                                ,minfra,maxfra
      if(ish.ge.7)then
      write(ifch,*)('-',l=1,68)
      write(ifch,*)'sr utclea. initial.'
      write(ifch,*)('-',l=1,68)
      do 34 n=nptli,nptl
      write(ifch,116)iorptl(n),jorptl(n),n,ifrptl(1,n),ifrptl(2,n)
     *,idptl(n),sqrt(pptl(1,n)**2+pptl(2,n)**2),pptl(3,n),pptl(5,n)
     *,istptl(n),ityptl(n)
34    continue
116   format(1x,i6,i6,4x,i6,4x,i6,i6,i12,3x,3(e8.2,1x),i3,i3)
      endif

c      ish=ish0
c      ish0=ish
c      if(ishsub/100.eq.18)ish=mod(ishsub,100)

      i=nptli-1
1     i=i+1
      if(i.gt.nptl)goto 1000
      if(istptl(i).eq.99)goto 2
      newptl(i)=i
c      oldptl(i)=i
      goto 1

2     i=i-1
      j=i
3     i=i+1
4     j=j+1
      if(j.gt.nptl)goto 5
      newptl(j)=0
      if(istptl(j).eq.99)goto 4
      newptl(j)=i
c      oldptl(i)=j
c      write(ifch,*)'move',j,' to ',i
c       write(ifch,*)idptl(i),ityptl(i),idptl(j),ityptl(j),minfra,maxfra
      call utrepl(i,j)
      if(j.ge.minfra0.and.j.le.maxfra0)then
        minfra1=min(minfra1,i)
        maxfra1=max(maxfra1,i)
      endif
      goto 3

5     nptl=i-1
      if(nptl.eq.0)then
        nptl0=0
        goto 1000
      endif

20    n0=newptl(nptl0)
      if(n0.gt.0)then
      nptl0=n0
      else
      nptl0=nptl0-1
      if(nptl0.gt.0)goto 20
      endif


c      do 11 k=1,nptl
c      io=iorptl(k)
c      if(io.le.0)ii(k)=io
c      if(io.gt.0)ii(k)=newptl(io)
c11    continue
c      do 12 k=1,nptl
c12    iorptl(k)=ii(k)
c
c      do 13 k=1,nptl
c      jo=jorptl(k)
c      if(jo.le.0)ii(k)=jo
c      if(jo.gt.0)ii(k)=newptl(jo)
c13    continue
c      do 14 k=1,nptl
c14    jorptl(k)=ii(k)
c
c      do 15 k=1,nptl
c      if1=ifrptl(1,k)
c      if(if1.le.0)ii(k)=if1
c      if(if1.gt.0)ii(k)=newptl(if1)
c15    continue
c      do 16 k=1,nptl
c16    ifrptl(1,k)=ii(k)
c
c      do 17 k=1,nptl
c      if2=ifrptl(2,k)
c      if(if2.le.0)ii(k)=if2
c      if(if2.gt.0)ii(k)=newptl(if2)
c17    continue
c      do 18 k=1,nptl
c18    ifrptl(2,k)=ii(k)
c
c      do 19 k=1,nptl
c      if(ifrptl(1,k).eq.0.and.ifrptl(2,k).gt.0)ifrptl(1,k)=ifrptl(2,k)
c      if(ifrptl(2,k).eq.0.and.ifrptl(1,k).gt.0)ifrptl(2,k)=ifrptl(1,k)
c19    continue

1000  continue

      if(minfra1.lt.minfra0)minfra=minfra1
      if(maxfra1.ge.minfra1)maxfra=maxfra1

      if(ish.ge.2)then
      write(ifch,*)'exiting subr utclea:'
      do 35 n=1,nptl
      write(ifch,116)iorptl(n),jorptl(n),n,ifrptl(1,n),ifrptl(2,n)
     *,idptl(n),sqrt(pptl(1,n)**2+pptl(2,n)**2),pptl(3,n),pptl(5,n)
     *,istptl(n),ityptl(n)
35    continue
      endif

      if(ish.ge.2)write(ifch,*)'exiting subr utclea:',nptl
     &                                                ,minfra,maxfra

      call utprix('utclea',ish,ishini,2)
      ish=ish0
      return
      end

c---------------------------------------------------------------------
      subroutine utfit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
c---------------------------------------------------------------------
c linear fit to data
c input:
c    ndata: nr of data points
c    x(),y(),sig(): data
c    mwt: unweighted (0) or weighted (else) data points
c output:
c    a,b: parameters of linear fit a+b*x
c---------------------------------------------------------------------
      INTEGER mwt,ndata
      REAL a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
CU    USES utgmq
      INTEGER i
      REAL sigdat,ss,st2,sx,sxoss,sy,t,wt,utgmq
      sx=0.
      sy=0.
      st2=0.
      b=0.
      if(mwt.ne.0) then
        ss=0.
        do 11 i=1,ndata
          wt=1./(sig(i)**2)
          ss=ss+wt
          sx=sx+x(i)*wt
          sy=sy+y(i)*wt
11      continue
      else
        do 12 i=1,ndata
          sx=sx+x(i)
          sy=sy+y(i)
12      continue
        ss=float(ndata)
      endif
      sxoss=sx/ss
      if(mwt.ne.0) then
        do 13 i=1,ndata
          t=(x(i)-sxoss)/sig(i)
          st2=st2+t*t
          b=b+t*y(i)/sig(i)
13      continue
      else
        do 14 i=1,ndata
          t=x(i)-sxoss
          st2=st2+t*t
          b=b+t*y(i)
14      continue
      endif
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1.+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1./st2)
      chi2=0.
      if(mwt.eq.0) then
        do 15 i=1,ndata
          chi2=chi2+(y(i)-a-b*x(i))**2
15      continue
        q=1.
        sigdat=sqrt(chi2/(ndata-2))
        siga=siga*sigdat
        sigb=sigb*sigdat
      else
        do 16 i=1,ndata
          chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
16      continue
        q=utgmq(0.5*(ndata-2),0.5*chi2)
      endif
      return
      END

c-----------------------------------------------------------------------
      function utgam1(x)
c-----------------------------------------------------------------------
c  gamma fctn tabulated
c  single precision
c-----------------------------------------------------------------------
      double precision utgamtab,utgam,al,dl
      common/gamtab/utgamtab(10000)

      if(x.gt.0.01.and.x.lt.99.99)then
        al=100.d0*dble(x)
        k1=int(al)
        k2=k1+1
        dl =al-dble(k1)
        utgam1=real(utgamtab(k2)*dl+utgamtab(k1)*(1.d0-dl))
      elseif(x.eq.0.)then
        utgam1=0.
      else
        utgam1=real(utgam(dble(x)))
      endif

      end

c-----------------------------------------------------------------------
      double precision function utgam2(x)
c-----------------------------------------------------------------------
c  gamma fctn tabulated
c  double precision
c-----------------------------------------------------------------------
      double precision utgamtab,x,al,dl,utgam
      common/gamtab/utgamtab(10000)

      if(x.gt.0.01d0.and.x.le.99.99d0)then
        al=100.d0*x
        k1=int(al)
        k2=k1+1
        dl =al-dble(k1)
        utgam2=utgamtab(k2)*dl+utgamtab(k1)*(1.d0-dl)
      elseif(x.eq.0.d0)then
        utgam2=0.d0
      else
        utgam2=utgam(x)
      endif

      end

c-----------------------------------------------------------------------
      double precision function utgam(x)
c-----------------------------------------------------------------------
c  gamma fctn
c  double precision
c-----------------------------------------------------------------------
      include 'epos.inc'
      double precision c(13),x,z,f
      data c
     1/ 0.00053 96989 58808, 0.00261 93072 82746, 0.02044 96308 23590,
     2  0.07309 48364 14370, 0.27964 36915 78538, 0.55338 76923 85769,
     3  0.99999 99999 99998,-0.00083 27247 08684, 0.00469 86580 79622,
     4  0.02252 38347 47260,-0.17044 79328 74746,-0.05681 03350 86194,
     5  1.13060 33572 86556/
      utgam=0d0
      z=x
      if(x .gt. 170.d0) goto6
      if(x .gt. 0.0d0) goto1
      if(x .eq. int(x)) goto5
      z=1.0d0-z
    1 f=1.0d0/z
      if(z .le. 1.0d0) goto4
      f=1.0d0
    2 continue
      if(z .lt. 2.0d0) goto3
      z=z-1.0d0
      f=f*z
      goto2
    3 z=z-1.0d0
    4 utgam=
     1 f*((((((c(1)*z+c(2))*z+c(3))*z+c(4))*z+c(5))*z+c(6))*z+c(7))/
     2 ((((((c(8)*z+c(9))*z+c(10))*z+c(11))*z+c(12))*z+c(13))*z+1.0d0)
      if(x .gt. 0.0d0) return
      utgam=3.141592653589793d0/(sin(3.141592653589793d0*x)*utgam)
      return
    5 write(ifch,10)sngl(x)
   10 format(1x,'argument of gamma fctn = ',e20.5)
      call utstop('utgam : negative integer argument&')
    6 write(ifch,11)sngl(x)
   11 format(1x,'argument of gamma fctn = ',e20.5)
      call utstop('utgam : argument too large&')
      end

c---------------------------------------------------------------------
      subroutine utgcf(gammcf,a,x,gln)
c---------------------------------------------------------------------
      INTEGER ITMAX
      REAL a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES utgmln
      INTEGER i
      REAL an,b,c,d,del,h,utgmln
      gln=utgmln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      call utstop("a too large, ITMAX too small in utgcf&")
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END

c---------------------------------------------------------------------
      function utgmln(xx)
c---------------------------------------------------------------------
      REAL utgmln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      utgmln=tmp+log(stp*ser/x)
      return
      END

c---------------------------------------------------------------------
      function utgmq(a,x)
c---------------------------------------------------------------------
      REAL a,utgmq,x
CU    USES utgcf,utgser
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.) call utstop("bad arguments in utgmq&")
      if(x.lt.a+1.)then
        call utgser(gamser,a,x,gln)
        utgmq=1.-gamser
      else
        call utgcf(gammcf,a,x,gln)
        utgmq=gammcf
      endif
      return
      END

c---------------------------------------------------------------------
      subroutine utgser(gamser,a,x,gln)
c---------------------------------------------------------------------
      INTEGER ITMAX
      REAL a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES utgmln
      INTEGER n
      REAL ap,del,sum,utgmln
      gln=utgmln(a)
      if(x.le.0.)then
        if(x.lt.0.)call utstop("x < 0 in utgser&")
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      call utstop("a too large, ITMAX too small in utgser&")
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END

c-------------------------------------------------------------------------
      subroutine uticpl(ic,ifla,iqaq,iret)
c-------------------------------------------------------------------------
c  adds a quark (iqaq=1) or antiquark (iqaq=2) of flavour ifla
c  to 2-id ic
c-------------------------------------------------------------------------
      include 'epos.inc'
      integer jc(nflav,2),ic(2)
      iret=0
      if(ifla.eq.0)return
      call iddeco(ic,jc)
      if(ish.ge.8)write(ifch,'(2i8,12i3)')ic,jc
      jqaq=3-iqaq
      if(jc(ifla,jqaq).gt.0)then
      jc(ifla,jqaq)=jc(ifla,jqaq)-1
      else
      jc(ifla,iqaq)=jc(ifla,iqaq)+1
      endif
      call idcomj(jc)
      call idenco(jc,ic,ireten)
      if(ish.ge.8)write(ifch,'(2i8,12i3)')ic,jc
      if(ireten.eq.1)iret=1
      if(ic(1).eq.0.and.ic(2).eq.0.and.ireten.eq.0)then
      ic(1)=100000
      ic(2)=100000
      endif
      return
      end

cc-----------------------------------------------------------------------
c      subroutine utindx(n,xar,x,i)
cc-----------------------------------------------------------------------
cc  input:  dimension n
cc          array xar(n) with xar(i) > xar(i-1)
cc          some number x between xar(1) and xar(n)
cc  output: the index i such that x is between xar(i)  and xar(i+1)
cc-----------------------------------------------------------------------
c      include 'epos.inc'
c      real xar(n)
c           if(x.lt.xar(1))then
c      if(ish.ge.5)then
c      call utmsg('utindx&')
c      write(ifch,*)'*****  x=',x,' < xar(1)=',xar(1)
c      call utmsgf
c      endif
c      i=1
c      return
c           elseif(x.gt.xar(n))then
c      if(ish.ge.5)then
c      call utmsg('utindx&')
c      write(ifch,*)'*****  x=',x,' > xar(n)=',xar(n)
c      call utmsgf
c      endif
c      i=n
c      return
c           endif
c      lu=1
c      lo=n
c1     lz=(lo+lu)/2
c      if((xar(lu).le.x).and.(x.le.xar(lz)))then
c      lo=lz
c      elseif((xar(lz).lt.x).and.(x.le.xar(lo)))then
c      lu=lz
c      else
c      call utstop('utindx: no interval found&')
c      endif
c      if((lo-lu).ge.2) goto1
c      if(lo.le.lu)call utstop('utinvt: lo.le.lu&')
c      i=lu
c      return
c      end
c
c-----------------------------------------------------------------------
      function utinvt(n,x,q,y)
c-----------------------------------------------------------------------
c     returns x with y=q(x)
c-----------------------------------------------------------------------
      include 'epos.inc'
      real x(n),q(n)
      if(q(n).eq.0.)call utstop('utinvt: q(n)=0&')
           if(y.lt.0.)then
      if(ish.ge.1)then
      call utmsg('utinvt&')
      write(ifch,*)'*****  y=',y,' < 0'
      call utmsgf
      endif
      y=0.
           elseif(y.gt.q(n))then
      if(ish.ge.1)then
      call utmsg('utinvt&')
      write(ifch,*)'*****  y=',y,' > ',q(n)
      call utmsgf
      endif
      y=q(n)
           endif
      lu=1
      lo=n
1     lz=(lo+lu)/2
      if((q(lu).le.y).and.(y.le.q(lz)))then
      lo=lz
      elseif((q(lz).lt.y).and.(y.le.q(lo)))then
      lu=lz
      else
      write(ifch,*)'q(1),y,q(n):',q(1),y,q(n)
      write(ifch,*)'lu,lz,lo:',lu,lz,lo
      write(ifch,*)'q(lu),q(lz),q(lo):',q(lu),q(lz),q(lo)
      call utstop('utinvt: no interval found&')
      endif
      if((lo-lu).ge.2) goto1
      if(lo.le.lu)call utstop('utinvt: lo.le.lu&')
      utinvt=x(lu)+(y-q(lu))*(x(lo)-x(lu))/(q(lo)-q(lu))
      return
      end

c-----------------------------------------------------------------------
      subroutine utlob2(isig,p1,p2,p3,p4,p5,x1,x2,x3,x4,idi)
c-----------------------------------------------------------------------
c  performs a lorentz boost, double prec.
c  isig=+1 is to boost the four vector x1,x2,x3,x4 such as to obtain it
c  in the frame specified by the 5-vector p1...p5 (5-vector=4-vector+mass).
c  isig=-1: the other way round, that means,
c  if the 4-vector x1...x4 is given in some frame characterized by
c  p1...p5 with respect to to some lab-frame, utlob2 returns the 4-vector
c  x1...x4  in the lab frame.
c  idi is a call identifyer (integer) to identify the call in case of problem
c-----------------------------------------------------------------------
      include 'epos.inc'
      double precision beta(4),z(4),p1,p2,p3,p4,p5,pp,bp,x1,x2,x3,x4
     *,xx0,x10,x20,x30,x40,x4x,x0123
           if(ish.ge.2)then
      if(ish.ge.9)then
      write(ifch,101)x1,x2,x3,x4,(x4-x3)*(x4+x3)-x2*x2-x1*x1
      write(ifch,301)p1,p2,p3,p4,p5,(p4-p3)*(p4+p3)-p2*p2-p1*p1
101   format(' utlob2: x =  ',5e13.5)
301   format('         p =  ',6e13.5)
      endif
      pp=(p4-p3)*(p4+p3)-p2*p2-p1*p1
      if(dabs(pp-p5*p5).gt.1d-3*p4*p4.and.dabs(pp-p5*p5).gt.1d-3)then
      call utmsg('utlob2&')
      write(ifch,*)'*****  p**2 .ne. p5**2'
      write(ifch,*)'call identifyer:',idi
      write(ifch,*)'p**2,p5**2: ',pp,p5*p5
      write(ifch,*)'p: ',p1,p2,p3,p4,p5
      call utmsgf
      endif
      x10=x1
      x20=x2
      x30=x3
      x40=x4
           endif
      xx0=(x4-x3)*(x4+x3)-x2*x2-x1*x1
      if(p5.le.0.)then
      call utmsg('utlob2&')
      write(ifch,*)'*****  p5 negative.'
      write(ifch,*)'call identifyer:',idi
      write(ifch,*)'p(5): ',p1,p2,p3,p4,p5
      write(ifmt,*)'call identifyer:',idi
      write(ifmt,*)'p(5): ',p1,p2,p3,p4,p5
      call utmsgf
c      force_crash=z(int(p5)) !to get the Backtrace after crash
      call utstop('utlob2: p5 negative.&')
      endif
      z(1)=x1
      z(2)=x2
      z(3)=x3
      z(4)=x4
      beta(1)=-p1/p5
      beta(2)=-p2/p5
      beta(3)=-p3/p5
      beta(4)= p4/p5
      bp=0.
      do 220 k=1,3
220   bp=bp+z(k)*isig*beta(k)
      do 230 k=1,3
230   z(k)=z(k)+isig*beta(k)*z(4)
     *+isig*beta(k)*bp/(beta(4)+1.)
      z(4)=beta(4)*z(4)+bp
      x1=z(1)
      x2=z(2)
      x3=z(3)
      x4=z(4)
      if(ish.ge.9)write(ifch,103)
c      write(ifmt,103)
     *x1,x2,x3,x4,(x4-x3)*(x4+x3)-x2*x2-x1*x1,xx0
 103  format(' utlob2: x =  ',5e13.5, ' M = ',e13.5)
      x4x=x4
      x0123=xx0+x1*x1+x2*x2+x3*x3
      if(x0123.gt.0.)then
      x4=sign( dsqrt(x0123) , x4x )
      else
      x4=0
      endif
      if(ish.ge.9)then
      write(ifch,101)x1,x2,x3,x4,(x4-x3)*(x4+x3)-x2*x2-x1*x1
      endif
           if(ish.ge.2)then
      if(ish.ge.9)write(ifch,*)'check x**2_ini -- x**2_fin'
      if(dabs(x4-x4x).gt.1d-2*dabs(x4).and.dabs(x4-x4x).gt.1d-2)then
      call utmsg('utlob2&')
      write(ifch,*)'*****  x**2_ini .ne. x**2_fin.'
      write(ifch,*)'call identifyer:',idi
      write(ifch,*)'x1 x2 x3 x4 x**2 (initial/final/corrected):'
102   format(5e13.5)
      write(ifch,102)x10,x20,x30,x40,(x40-x30)*(x40+x30)-x20*x20-x10*x10
      write(ifch,102)x1,x2,x3,x4x,(x4x-x3)*(x4x+x3)-x2*x2-x1*x1
      write(ifch,102)x1,x2,x3,x4,(x4-x3)*(x4+x3)-x2*x2-x1*x1
      call utmsgf
      endif
           endif
      if(ish.ge.9)write(ifch,*)'return from utlob2'
      return
      end

c-----------------------------------------------------------------------
      subroutine utlob3(isig,p1,p2,p3,p4,p5,x1,x2,x3,x4)
c-----------------------------------------------------------------------
c  performs a lorentz boost, double prec.
c  but arguments are single precision
c-----------------------------------------------------------------------
      double precision xx1,xx2,xx3,xx4,utdble
      xx1=utdble(x1)
      xx2=utdble(x2)
      xx3=utdble(x3)
      xx4=utdble(x4)
      call utlob2(isig
     *,utdble(p1),utdble(p2),utdble(p3),utdble(p4),utdble(p5)
     *,xx1,xx2,xx3,xx4,52)
      x1=sngl(xx1)
      x2=sngl(xx2)
      x3=sngl(xx3)
      x4=sngl(xx4)
      return
      end

c-----------------------------------------------------------------------
      subroutine utlob5(yboost,x1,x2,x3,x4,x5)
c-----------------------------------------------------------------------
      amt=sqrt(x5**2+x1**2+x2**2)
      y=sign(1.,x3)*alog((x4+abs(x3))/amt)
      y=y-yboost
      x4=amt*cosh(y)
      x3=amt*sinh(y)
      return
      end

c-----------------------------------------------------------------------
      subroutine utlob4(isig,pp1,pp2,pp3,pp4,pp5,x1,x2,x3,x4)
c-----------------------------------------------------------------------
c  performs a lorentz boost, double prec.
c  but arguments are partly single precision
c-----------------------------------------------------------------------
      double precision xx1,xx2,xx3,xx4,pp1,pp2,pp3,pp4,pp5,utdble
      xx1=utdble(x1)
      xx2=utdble(x2)
      xx3=utdble(x3)
      xx4=utdble(x4)
      call utlob2(isig,pp1,pp2,pp3,pp4,pp5,xx1,xx2,xx3,xx4,53)
      x1=sngl(xx1)
      x2=sngl(xx2)
      x3=sngl(xx3)
      x4=sngl(xx4)
      return
      end


c-----------------------------------------------------------------------
      subroutine utlobo(isig,p1,p2,p3,p4,p5,x1,x2,x3,x4)
c-----------------------------------------------------------------------
c     performs a lorentz boost
c-----------------------------------------------------------------------
      include 'epos.inc'
      real beta(4),z(4)
      if(p5.le.0.)then
      call utmsg('utlobo&')
      write(ifch,*)'*****  mass <= 0.'
      write(ifch,*)'p(5): ',p1,p2,p3,p4,p5
      call utmsgf
      call utstop('utlobo: mass <= 0.&')
      endif
      z(1)=x1
      z(2)=x2
      z(3)=x3
      z(4)=x4
      beta(1)=-p1/p5
      beta(2)=-p2/p5
      beta(3)=-p3/p5
      beta(4)= p4/p5
      bp=0.
      do 220 k=1,3
220   bp=bp+z(k)*isig*beta(k)
      do 230 k=1,3
230   z(k)=z(k)+isig*beta(k)*z(4)
     *+isig*beta(k)*bp/(beta(4)+1.)
      z(4)=beta(4)*z(4)+bp
      x1=z(1)
      x2=z(2)
      x3=z(3)
      x4=z(4)
      return
      end

c-----------------------------------------------------------------------
      subroutine utloc(ar,n,a,l)
c-----------------------------------------------------------------------
      real ar(n)
      do 1 i=1,n
      l=i-1
      if(a.lt.ar(i))return
1     continue
      l=n
      return
      end

cc-----------------------------------------------------------------------
c      subroutine utlow(cone)
cc-----------------------------------------------------------------------
c      character*1 cone
c      if(cone.eq.'A')cone='a'
c      if(cone.eq.'B')cone='b'
c      if(cone.eq.'C')cone='c'
c      if(cone.eq.'D')cone='d'
c      if(cone.eq.'E')cone='e'
c      if(cone.eq.'F')cone='f'
c      if(cone.eq.'G')cone='g'
c      if(cone.eq.'H')cone='h'
c      if(cone.eq.'I')cone='i'
c      if(cone.eq.'J')cone='j'
c      if(cone.eq.'K')cone='k'
c      if(cone.eq.'L')cone='l'
c      if(cone.eq.'M')cone='m'
c      if(cone.eq.'N')cone='n'
c      if(cone.eq.'O')cone='o'
c      if(cone.eq.'P')cone='p'
c      if(cone.eq.'Q')cone='q'
c      if(cone.eq.'R')cone='r'
c      if(cone.eq.'S')cone='s'
c      if(cone.eq.'T')cone='t'
c      if(cone.eq.'U')cone='u'
c      if(cone.eq.'V')cone='v'
c      if(cone.eq.'W')cone='w'
c      if(cone.eq.'X')cone='x'
c      if(cone.eq.'Y')cone='y'
c      if(cone.eq.'Z')cone='z'
c      return
c      end
c
cc-----------------------------------------------------------------------
c      subroutine utlow3(cthree)
cc-----------------------------------------------------------------------
c      character cthree*3
c      do 1 i=1,3
c1     call utlow(cthree(i:i))
c      return
c      end
c
cc-----------------------------------------------------------------------
c      subroutine utlow6(csix)
cc-----------------------------------------------------------------------
c      character csix*6
c      do 1 i=1,6
c1     call utlow(csix(i:i))
c      return
c      end
c
cc-----------------------------------------------------------------------
c      function utmom(k,n,x,q)
cc-----------------------------------------------------------------------
cc     calculates kth moment for f(x) with q(i)=int[0,x(i)]f(z)dz
cc-----------------------------------------------------------------------
c      real x(n),q(n)
c      if(n.lt.2)call utstop('utmom : dimension too small&')
c      utmom=0
c      do 1 i=2,n
c1     utmom=utmom+((x(i)+x(i-1))/2)**k*(q(i)-q(i-1))
c      utmom=utmom/q(n)
c      return
c      end
c
c-----------------------------------------------------------------------
      function utpcm(a,b,c)
c-----------------------------------------------------------------------
c     calculates cm momentum for a-->b+c
c-----------------------------------------------------------------------
      val=(a*a-b*b-c*c)*(a*a-b*b-c*c)-(2.*b*c)*(2.*b*c)
      if(val.lt.0..and.val.gt.-1e-4)then
      utpcm=0
      return
      endif
      utpcm=sqrt(val)/(2.*a)
      return
      end

c-----------------------------------------------------------------------
      function utpcmi(a,b,c,iret)
c-----------------------------------------------------------------------
c     calculates cm momentum for a-->b+c
c-----------------------------------------------------------------------
      iret=0
      utpcmi=0
      val=(a*a-b*b-c*c)*(a*a-b*b-c*c)-(2.*b*c)*(2.*b*c)
      if(val.lt.0..and.val.gt.-1e-4)then
        return
      endif
      if(val.lt.0)then
        iret=1
        return
      endif
      utpcmi=sqrt(val)/(2.*a)
      return
      end

c-----------------------------------------------------------------------
      double precision function utpcmd(a,b,c,iret)
c-----------------------------------------------------------------------
c     calculates cm momentum for a-->b+c
c-----------------------------------------------------------------------
      double precision a,b,c,val
      iret=0
      val=(a*a-b*b-c*c)*(a*a-b*b-c*c)-(2.*b*c)*(2.*b*c)
      utpcmd=0d0
      if(val.lt.0d0.and.val.gt.-1d-4)then
        return
      elseif(val.lt.0d0)then
        iret=1
        return
      endif
      utpcmd=sqrt(val)/(2.d0*a)
      return
      end

c-----------------------------------------------------------------------
      subroutine utpri(text,ishi,ishini,ishx)
c-----------------------------------------------------------------------
      include 'epos.inc'
      character*(*) text
c      double precision seedx                               !!!
      ishini=ishi
      if(ishevt.ne.0.and.nrevt+1.ne.ishevt)return
      if(jprint.ge.ishx)then
        write(ifmt,'(1x,43a)')
     *  ('-',i=1,10),' entry ',text,' ',('-',i=1,30)
        if(ifmt.ne.6)close(ifmt)
        if(ifmt.ne.6)open(ifmt,file=fnmt(1:nfnmt),access='append')
      endif
      if(nrpri.gt.0)then
        do nr=1,nrpri
        if(subpri(nr)(1:6).eq.text)then
        ishi=ishpri(nr)
        endif
        enddo
      endif
      if(ishi.ge.ishx)then
        write(ifch,'(1x,43a)')
     *  ('-',i=1,10),' entry ',text,' ',('-',i=1,30)
c       call ranfgt(seedx)                                   !!!
c       if(ishi.ge.ishx)write(ifch,*)'seed:',seedx            !!!
      endif
      return
      end

c-----------------------------------------------------------------------
      subroutine utprix(text,ishi,ishini,ishx)
c-----------------------------------------------------------------------
      include 'epos.inc'
      character*(*) text
      if(ishevt.ne.0.and.nrevt+1.ne.ishevt)return
      if(jprint.ge.ishx)then
        write(ifmt,'(1x,44a)')
     *  ('-',i=1,30),' exit ',text,' ',('-',i=1,11)
        if(ifmt.ne.6)close(ifmt)
        if(ifmt.ne.6)open(ifmt,file=fnmt(1:nfnmt),access='append')
      endif
      if(ish.ge.ishx)write(ifch,'(1x,44a)')
     *('-',i=1,30),' exit ',text,' ',('-',i=1,11)
      ishi=ishini
      return
      end

c-----------------------------------------------------------------------
      subroutine utprj(text,ishi,ishini,ishx)
c-----------------------------------------------------------------------
      include 'epos.inc'
      character*(*) text
c      double precision seedx                               !!!
      idx=index(text,' ')-1
      ishini=ishi
      if(ishevt.ne.0.and.nrevt+1.ne.ishevt)return
      if(nrpri.gt.0)then
      do nr=1,nrpri
      if(subpri(nr)(1:idx).eq.text(1:idx))then
      ishi=ishpri(nr)
      endif
      enddo
      endif
      if(ish.ge.ishx)then
        write(ifch,'(1x,43a)')
     *  ('-',i=1,10),' entry ',text(1:idx),' ',('-',i=1,30)
c       call ranfgt(seedx)                                   !!!
c       if(ish.ge.ishx)write(ifch,*)'seed:',seedx            !!!
      endif
      return
      end

c-----------------------------------------------------------------------
      subroutine utprjx(text,ishi,ishini,ishx)
c-----------------------------------------------------------------------
      include 'epos.inc'
      character*(*) text
      idx=index(text,' ')-1
      if(ishevt.ne.0.and.nrevt+1.ne.ishevt)return
      if(ish.ge.ishx)write(ifch,'(1x,44a)')
     *('-',i=1,30),' exit ',text(1:idx),' ',('-',i=1,11)
      ishi=ishini
      return
      end

c-----------------------------------------------------------------------
      function utquad(m,x,f,k)
c-----------------------------------------------------------------------
c     performs an integration according to simpson
c-----------------------------------------------------------------------
      real x(m),f(m)
      utquad=0
      do 1 i=1,k-1
  1   utquad=utquad+(f(i)+f(i+1))/2*(x(i+1)-x(i))
      return
      end

c-----------------------------------------------------------------------
      subroutine utquaf(fu,n,x,q,x0,x1,x2,x3)
c-----------------------------------------------------------------------
c     returns q(i) = integral [x(1)->x(i)] fu(x) dx
c-----------------------------------------------------------------------
      include 'epos.inc'
      real x(n),q(n)
      parameter (m=10)
      real xa(m),fa(m)
      external fu
      if(x1.lt.x0.or.x2.lt.x1.or.x3.lt.x2)then
      if(ish.ge.1)then
      call utmsg('utquaf&')
      write(ifch,*)'   xi=',x0,x1,x2,x3
      call utmsgf
      endif
      endif
      call utar(n/3,n*2/3,n,x0,x1,x2,x3,x)
      q(1)=0
      do 2 i=2,n
      do 3 k=1,m
      z=x(i-1)+(k-1.)/(m-1.)*(x(i)-x(i-1))
      xa(k)=z
3     fa(k)=fu(z)
      q(i)=q(i-1)+utquad(m,xa,fa,m)
2     continue
      return
      end

c-----------------------------------------------------------------------
      subroutine utrepl(i,j)
c-----------------------------------------------------------------------
c     i is replaced by j in /cptl/
c-----------------------------------------------------------------------
      include 'epos.inc'
      do 1 k=1,5
1     pptl(k,i)  =pptl(k,j)
      iorptl(i)  = 0 !iorptl(j)
      idptl(i)   =idptl(j)
      istptl(i)  =istptl(j)
      do 2 k=1,2
2     tivptl(k,i)=tivptl(k,j)
      do 3 k=1,2
3     ifrptl(k,i)= 0 !ifrptl(k,j)
      jorptl(i)  = 0 !jorptl(j)
      do 4 k=1,4
4     xorptl(k,i)=xorptl(k,j)
      do 5 k=1,4
5     ibptl(k,i) =ibptl(k,j)
      ityptl(i)  =ityptl(j)
      iaaptl(i)  =iaaptl(j)
      radptl(i)  =radptl(j)
      desptl(i)  =desptl(j)
      dezptl(i)  =dezptl(j)
      qsqptl(i)  =qsqptl(j)
      zpaptl(1,i)=zpaptl(1,j)
      zpaptl(2,i)=zpaptl(2,j)
      itsptl(i)  =itsptl(j)
      rinptl(i)  =rinptl(j)
      return
      end

c-----------------------------------------------------------------------
      subroutine utrepla(i,j)
c-----------------------------------------------------------------------
c     i is replaced by j in /cptl/
c-----------------------------------------------------------------------
      include 'epos.inc'
      do 1 k=1,5
1     pptl(k,i)  =pptl(k,j)
      iorptl(i)  = iorptl(j)
      idptl(i)   =idptl(j)
      istptl(i)  =istptl(j)
      do 2 k=1,2
2     tivptl(k,i)=tivptl(k,j)
      do 3 k=1,2
3     ifrptl(k,i)= ifrptl(k,j)
      jorptl(i)  = jorptl(j)
      do 4 k=1,4
4     xorptl(k,i)=xorptl(k,j)
      do 5 k=1,4
5     ibptl(k,i) =ibptl(k,j)
      ityptl(i)  =ityptl(j)
      iaaptl(i)  =iaaptl(j)
      radptl(i)  =radptl(j)
      desptl(i)  =desptl(j)
      dezptl(i)  =dezptl(j)
      qsqptl(i)  =qsqptl(j)
      zpaptl(1,i)=zpaptl(1,j)
      zpaptl(2,i)=zpaptl(2,j)
      itsptl(i)  =itsptl(j)
      rinptl(i)  =rinptl(j)
      return
      end

cc-----------------------------------------------------------------------
c      subroutine utresm(icp1,icp2,icm1,icm2,amp,idpr,iadj,ireten)
cc-----------------------------------------------------------------------
c      parameter (nflav=6)
c      integer icm(2),icp(2),jcm(nflav,2),jcp(nflav,2)
c      icm(1)=icm1
c      icm(2)=icm2
c      icp(1)=icp1
c      icp(2)=icp2
c      CALL IDDECO(ICM,JCM)
c      CALL IDDECO(ICP,JCP)
c      do 37 nf=1,nflav
c      do 37 k=1,2
c37    jcP(nf,k)=jcp(nf,k)+jcm(nf,k)
c      CALL IDENCO(JCP,ICP,IRETEN)
c      IDP=IDTRA(ICP,0,0,3)
c      call idres(idp,amp,idpr,iadj,1)
c      return
c      end
c
cc-----------------------------------------------------------------------
c      subroutine utroa1(phi,a1,a2,a3,x1,x2,x3)
cc-----------------------------------------------------------------------
cc  rotates x by angle phi around axis a.
cc  normalization of a is irrelevant.
cc-----------------------------------------------------------------------
c      double precision aaa,aa(3),xxx,xx(3),e1(3),e2(3),e3(3),xp,xt,dphi
c      dphi=phi
c      xx(1)=x1
c      xx(2)=x2
c      xx(3)=x3
c      aa(1)=a1
c      aa(2)=a2
c      aa(3)=a3
c      aaa=0
c      xxx=0
c      do i=1,3
c      aaa=aaa+aa(i)**2
c      xxx=xxx+xx(i)**2
c      enddo
c      if(xxx.eq.0d0)return
c      if(aaa.eq.0d0)call utstop('utroa1: zero rotation axis&')
c      aaa=dsqrt(aaa)
c      xxx=dsqrt(xxx)
cc e3 = a / !a!
c      do i=1,3
c      e3(i)=aa(i)/aaa
c      enddo
cc x_parallel
c      xp=0
c      do i=1,3
c      xp=xp+xx(i)*e3(i)
c      enddo
cc x_transverse
c      if(xxx**2-xp**2.le.0.)return
c      xt=dsqrt(xxx**2-xp**2)
cc e1 = vector x_transverse / absolute value x_transverse
c      do i=1,3
c      e1(i)=(xx(i)-e3(i)*xp)/xt
c      enddo
cc e2 orthogonal e3,e1
c      call utvec2(e3,e1,e2)
cc rotate x
c      do i=1,3
c      xx(i)=xp*e3(i)+xt*dcos(dphi)*e1(i)+xt*dsin(dphi)*e2(i)
c      enddo
cc back to single precision
c      x1=xx(1)
c      x2=xx(2)
c      x3=xx(3)
c      return
c      end
c
c-----------------------------------------------------------------------
      subroutine utroa1(phi,a1,a2,a3,x1,x2,x3)
c-----------------------------------------------------------------------
c  rotates x by angle phi around axis a (argument single precision)
c  normalization of a is irrelevant.
c-----------------------------------------------------------------------
      double precision aa(3),xx(3),dphi
      dphi=phi
      xx(1)=x1
      xx(2)=x2
      xx(3)=x3
      aa(1)=a1
      aa(2)=a2
      aa(3)=a3
      call utroa2(dphi,aa(1),aa(2),aa(3),xx(1),xx(2),xx(3))
c back to single precision
      x1=sngl(xx(1))
      x2=sngl(xx(2))
      x3=sngl(xx(3))
      return
      end

c-----------------------------------------------------------------------
      subroutine utroa2(phi,a1,a2,a3,x1,x2,x3)
c-----------------------------------------------------------------------
c  rotates x by angle phi around axis a.
c  normalization of a is irrelevant.
c  double precision phi,a1,a2,a3,x1,x2,x3
c-----------------------------------------------------------------------
      double precision phi,a1,a2,a3,x1,x2,x3
      double precision aaa,aa(3),xxx,xx(3),e1(3),e2(3),e3(3),xp,xt,dphi
      dphi=phi
      xx(1)=x1
      xx(2)=x2
      xx(3)=x3
      aa(1)=a1
      aa(2)=a2
      aa(3)=a3
      aaa=0d0
      xxx=0d0
      do i=1,3
      aaa=aaa+aa(i)**2
      xxx=xxx+xx(i)**2
      enddo
      if(xxx.eq.0d0)return
      if(aaa.eq.0d0)call utstop('utroa1: zero rotation axis&')
      aaa=1.0/dsqrt(aaa)
c e3 = a / !a!
      do i=1,3
      e3(i)=aa(i)*aaa
      enddo
c x_parallel
      xp=0
      do i=1,3
      xp=xp+xx(i)*e3(i)
      enddo
c x_transverse
      if(xxx-xp**2.le.0.)return
      xt=dsqrt(xxx-xp**2)
c e1 = vector x_transverse / absolute value x_transverse
      do i=1,3
      e1(i)=(xx(i)-e3(i)*xp)/xt
      enddo
c e2 orthogonal e3,e1
      call utvec2(e3,e1,e2)
c rotate x
      do i=1,3
      xx(i)=xp*e3(i)+xt*cos(dphi)*e1(i)+xt*sin(dphi)*e2(i)
      enddo
      xxx=0d0
      do i=1,3
      xxx=xxx+xx(i)**2
      enddo
c back to single precision
      x1=xx(1)
      x2=xx(2)
      x3=xx(3)
      return
      end

cc--------------------------------------------------------------------
c      function utroot(funcd,x1,x2,xacc)
cc--------------------------------------------------------------------
cc combination of newton-raphson and bisection method for root finding
cc input:
cc   funcd: subr returning fctn value and first derivative
cc   x1,x2: x-interval
cc   xacc:  accuracy
cc output:
cc   utroot: root
cc--------------------------------------------------------------------
c      include 'epos.inc'
c      INTEGER MAXIT
c      REAL utroot,x1,x2,xacc
c      EXTERNAL funcd
c      PARAMETER (MAXIT=100)
c      INTEGER j
c      REAL df,dx,dxold,f,fh,fl,temp,xh,xl
c      call funcd(x1,fl,df)
c      call funcd(x2,fh,df)
c      if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.))
c     *call utstop('utroot: root must be bracketed&')
c      if(fl.eq.0.)then
c        utroot=x1
c        return
c      else if(fh.eq.0.)then
c        utroot=x2
c        return
c      else if(fl.lt.0.)then
c        xl=x1
c        xh=x2
c      else
c        xh=x1
c        xl=x2
c      endif
c      utroot=.5*(x1+x2)
c      dxold=abs(x2-x1)
c      dx=dxold
c      call funcd(utroot,f,df)
c      do 11 j=1,MAXIT
c        if(((utroot-xh)*df-f)*((utroot-xl)*df-f).ge.0..or. abs(2.*
c     *f).gt.abs(dxold*df) ) then
c          dxold=dx
c          dx=0.5*(xh-xl)
c          utroot=xl+dx
c          if(xl.eq.utroot)return
c        else
c          dxold=dx
c          dx=f/df
c          temp=utroot
c          utroot=utroot-dx
c          if(temp.eq.utroot)return
c        endif
c        if(abs(dx).lt.xacc) return
c        call funcd(utroot,f,df)
c        if(f.lt.0.) then
c          xl=utroot
c        else
c          xh=utroot
c        endif
c11    continue
c      call utmsg('utroot&')
c      write(ifch,*)'*****  exceeding maximum iterations'
c      write(ifch,*)'dx:',dx
c      call utmsgf
c      return
c      END
c
c-----------------------------------------------------------------------
      subroutine utrot2(isig,ax,ay,az,x,y,z)
c-----------------------------------------------------------------------
c     performs a rotation, double prec.
c-----------------------------------------------------------------------
      include 'epos.inc'
      double precision ax,ay,az,x,y,z,rx,ry,rz
     *,alp,bet,cosa,sina,cosb,sinb,xs,ys,zs
         if(ax**2.eq.0.and.ay**2.eq.0.and.az**2.eq.0.)then
      write(ifch,*)'ax**2,ay**2,az**2:',ax**2,ay**2,az**2
      write(ifch,*)'ax,ay,az:',ax,ay,az
      call utstop('utrot2: zero vector.&')
         endif
         if(az.ge.0.)then
      rx=ax
      ry=ay
      rz=az
         else
      rx=-ax
      ry=-ay
      rz=-az
         endif
      if(rz**2+ry**2.ne.0.)then
      alp=dabs(dacos(rz/dsqrt(rz**2+ry**2)))*sign(1.,sngl(ry))
      bet=
     *dabs(dacos(dsqrt(rz**2+ry**2)/dsqrt(rz**2+ry**2+rx**2)))*
     *sign(1.,sngl(rx))
      else
      alp=3.1415927d0/2d0
      bet=3.1415927d0/2d0
      endif
      cosa=dcos(alp)
      sina=dsin(alp)
      cosb=dcos(bet)
      sinb=dsin(bet)
           if(isig.ge.0)then
      xs=x*cosb-y*sina*sinb-z*cosa*sinb
      ys=       y*cosa     -z*sina
      zs=x*sinb+y*sina*cosb+z*cosa*cosb
           else     !if(isig.lt.0)then
      xs= x*cosb            +z*sinb
      ys=-x*sinb*sina+y*cosa+z*cosb*sina
      zs=-x*sinb*cosa-y*sina+z*cosb*cosa
           endif
      x=xs
      y=ys
      z=zs
      return
      end


c-----------------------------------------------------------------------
      subroutine utrot4(isig,ax,ay,az,x,y,z)
c-----------------------------------------------------------------------
c     performs a rotation, double prec.
c     arguments partly single
c-----------------------------------------------------------------------
      double precision ax,ay,az,xx,yy,zz,utdble
      xx=utdble(x)
      yy=utdble(y)
      zz=utdble(z)
      call utrot2(isig,ax,ay,az,xx,yy,zz)
      x=sngl(xx)
      y=sngl(yy)
      z=sngl(zz)
      return
      end

c-----------------------------------------------------------------------
      subroutine utrota(isig,ax,ay,az,x,y,z)
c-----------------------------------------------------------------------
c     performs a rotation
c-----------------------------------------------------------------------
         if(az.ge.0.)then
      rx=ax
      ry=ay
      rz=az
         else
      rx=-ax
      ry=-ay
      rz=-az
         endif
      if(rz.eq.0..and.ry.eq.0.)then
        alp=0.
        stop
      else
        alp=abs(utacos(rz/sqrt(rz**2+ry**2)))*sign(1.,ry)
      endif
      bet=
     *abs(utacos(sqrt(rz**2+ry**2)/sqrt(rz**2+ry**2+rx**2)))*sign(1.,rx)
      cosa=cos(alp)
      sina=sin(alp)
      cosb=cos(bet)
      sinb=sin(bet)
           if(isig.ge.0)then
      xs=x*cosb-y*sina*sinb-z*cosa*sinb
      ys=       y*cosa     -z*sina
      zs=x*sinb+y*sina*cosb+z*cosa*cosb
           else        !if(isig.lt.0)then
      xs= x*cosb            +z*sinb
      ys=-x*sinb*sina+y*cosa+z*cosb*sina
      zs=-x*sinb*cosa-y*sina+z*cosb*cosa
           endif
      x=xs
      y=ys
      z=zs
      return
      end

c-----------------------------------------------------------------------
      subroutine utrotb(alp,bet,x,y,z)
c-----------------------------------------------------------------------
c     performs a rotation
c-----------------------------------------------------------------------
      cosa=cos(alp)
      sina=sin(alp)
      cosb=cos(bet)
      sinb=sin(bet)
      xs=x*cosb-y*sina*sinb-z*cosa*sinb
      ys=       y*cosa     -z*sina
      zs=x*sinb+y*sina*cosb+z*cosa*cosb
      x=xs
      y=ys
      z=zs
      return
      end


c-----------------------------------------------------------------------
      subroutine utstop(text)
c-----------------------------------------------------------------------
c  returns error message and stops execution.
c  text is an optonal text to appear in the error message.
c  text is a character string of length 40;
c     for shorter text, it has to be terminated by &;
c        example: call utstop('error in subr xyz&')
c-----------------------------------------------------------------------
      include 'epos.inc'
c      parameter(itext=40)
      character  text*(*)
      imax=index(text,'&')
      do 1 j=1,2
      if(j.eq.1)then
        ifi=ifch
      else        !if(j.eq.2)
        ifi=ifmt
      endif
      if(imax.gt.1)then
      write(ifi,101)('*',k=1,72),text(1:imax-1)
     *,nrevt+1,nint(seedj),seedc,('*',k=1,72)
      else
      write(ifi,101)('*',k=1,72),' '
     *,nrevt+1,nint(seedj),seedc,('*',k=1,72)
      endif
101   format(
     *1x,72a1
     */1x,'***** STOP in ',a
     */1x,'***** current event number: ',i12
     */1x,'***** initial seed for current run:',i10
     */1x,'***** initial seed for current event:',d25.15
     */1x,72a1)
1     continue
c      c=0.
c      b=a/c
      stop
      entry utmsg(text)
      imax=index(text,'&')
      if(imax.eq.0)imax=6
      imsg=imsg+1
      write(ifch,'(1x,74a1)')('*',j=1,72)
      write(ifch,100)text(1:imax),nrevt+1,nint(seedj),seedc
100   format(1x,'***** msg from ',a,'.   es:',i7,2x,i9,2x,d23.17)
      return
      entry utmsgf
      if(ish.eq.1)return
      write(ifch,'(1x,74a1)')('*',j=1,72)
      end

c-----------------------------------------------------------------
      subroutine uttrap(func,a,b,s)
c-----------------------------------------------------------------
c trapezoidal method for integration.
c input: fctn func and limits a,b
c output: value s of the integral
c-----------------------------------------------------------------
      include 'epos.inc'

      INTEGER JMAX
      REAL a,b,func,s
      EXTERNAL func
      PARAMETER (JMAX=10)
CU    USES uttras
      INTEGER j
      REAL olds

      ierr=0
      olds=-1.e30
      do 11 j=1,JMAX
        if(ish.ge.9)write(ifch,*)'sr uttrap:   j:',j
        call uttras(func,a,b,s,j)
        ds=abs(s-olds)
        if (ds.lt.epsr*abs(olds)) return
        olds=s
11    continue

c-c   nepsr=nepsr+1
      if(ish.ge.9)then
      call utmsg('uttrap&')
      write(ifch,*)
     *'*****  requested accuracy could not be achieved'
      write(ifch,*)'achieved accuracy: ',ds/abs(olds)
      write(ifch,*)'requested accuracy:',epsr
      call utmsgf
      endif

      END

c-----------------------------------------------------------------
      subroutine uttraq(func,a,b,s)
c-----------------------------------------------------------------
c trapezoidal method for integration.
c input: function func and limits a,b
c output: value s of the integral
c-----------------------------------------------------------------

      REAL a,b,func,s
      EXTERNAL func
      PARAMETER (eps=5.e-6)
CU    USES uttras
      INTEGER j
      REAL olds
      olds=-1.e30
      j=1
10      call uttras(func,a,b,s,j)
        ds=abs(s-olds)
        if (ds.le.eps*abs(olds)) return
        olds=s
        if(j.ge.15)then
          print *,"precision not reached",ds/olds
          return
        endif
        j=j+1
      goto 10
      END

c-----------------------------------------------------------------
      subroutine uttras(func,a,b,s,n)
c-----------------------------------------------------------------
c performs one iteration of the trapezoidal method for integration
c-----------------------------------------------------------------
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END

cc-----------------------------------------------------------------------
c      subroutine utvec1(a,b,c)
cc-----------------------------------------------------------------------
cc  returns vector product c = a x b .
cc-----------------------------------------------------------------------
c      real a(3),b(3),c(3)
c      c(1)=a(2)*b(3)-a(3)*b(2)
c      c(2)=a(3)*b(1)-a(1)*b(3)
c      c(3)=a(1)*b(2)-a(2)*b(1)
c      return
c      end
c
c-----------------------------------------------------------------------
      subroutine utvec2(a,b,c)
c-----------------------------------------------------------------------
c  returns vector product c = a x b .
c  a,b,c double precision.
c-----------------------------------------------------------------------
      double precision a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end

c-------------------------------------------------------------------
      subroutine utword(line,i,j,iqu)
c-------------------------------------------------------------------
c  finds the first word of the character string line(j+1:1000).
c  the word is line(i:j) (with new i and j).
c  if j<0 or if no word found --> new line read.
c  a text between quotes "..." is considered one word;
c  stronger: a text between double quotes ""..."" is consid one word
c  stronger: a text between "{ and }" is considered one word
c-------------------------------------------------------------------
c  input:
c    line: character string (*1000)
c    i: integer between 1 and 1000
c    iqu: for iqu=1 a "?" is written to output before reading a line,
c         otherwise (iqu/=1) nothing is typed
c  output:
c    i,j: left and right end of word (word=line(i:j))
c-------------------------------------------------------------------
      include 'epos.inc'
      parameter(mempty=1)
      character*1 empty(mempty),mk
      character line*(*)
      character*200 linex
      character*2 mrk
      data empty/' '/
      save

      imax=len(line)

      if(j.ge.0)then
        i=j
        goto 1
      endif

    5 continue
      if(macrocurrent.ge.1.and.macrocurrentline.ge.1)then
        if(macrocurrentline.le.macrotextlines(macrocurrent))goto 4
      endif

      !----------------------------------
      !     read line
      !----------------------------------
      if(iqu.eq.1.and.iprmpt.gt.0)write(ifmt,'(a)')'?'
      if(nopen.eq.0)then
        ifopx=ifop
      elseif(nopen.gt.0)then
        ifopx=20+nopen
      else !if(nopen.lt.0)
        ifopx=ifcp
      endif
      read(ifopx,'(a)',end=9999)line
c      print *,ifopx,line
      imx=imax
      do while(imx.gt.2.and.line(imx:imx).eq.' ')
        imx=imx-1
      enddo
      imx=imx+1

      !----------------------------------------------------------------
      !   #do 
      !----------------------------------------------------------------
      !do loop to be used in optns as
      !#do ( val-1-1 ... val-1-N , val-2 ... val-2-N , ... )
      ! ... several lines containg place holders {{i}}
      !#done
      !where each {{i}} will be replaced in the jth cycle by val-i-j  
      !----------------------------------------------------------------
      if(line(1:4).eq.'#do ')then
        i=5
        do while(line(i:i).eq.' ')
          i=i+1
        enddo
        if(line(i:i).ne.'(')stop'ERROR 26042021b'
        i=i+1
        k1=1
        k2=1
        igo=1 
        do while(igo.eq.1)
          do while(line(i:i).eq.' ')
            i=i+1
            if(i.gt.imax)then
              print*,'ERROR utword i > ',imax
              print*,line
              stop
            endif
          enddo
          iquote=0
          if(line(i:i).eq.'"')then !quotes may be used "..." to accomodate spaces
            iquote=1
            i=i+1
          endif 
          j=i
          do while(line(j:j).ne.' '.or.iquote.eq.1)
            if(line(j:j).eq.'"')iquote=0 !expresson ends with "
            j=j+1
          enddo
          !print*,'CHECK #do parameters: ',k1,k2,i,j,line(i:j)
          if(line(i:j).ne.') ')then
            if(line(i:j).ne.', ')then
              dooptnsargtext(k1,k2)=line(i:j)
            else
              k2=k2+1
              k1=0 
              if(k2.gt.k2domax)stop'ERROR 26042021c'
            endif
            i=j+1
            k1=k1+1
            if(k1.gt.k1domax)stop'ERROR 26042021d'
          else
            igo=0
          endif
        enddo
        k1dooptns=k1
        k2dooptns=k2
        idooptns=1
        if(nopen+2.gt.9)stop'ERROR 26042021a'
        open(unit=20+nopen+2,file=fncp(1:nfncp)//'Do2',status='unknown')
        open(unit=20+nopen+1,file=fncp(1:nfncp)//'Do1',status='unknown')
        nrdox=0 
        goto 5 
      endif
      !--------------#done----------------------------------------
      !End of do domain, close tmp file, make it current read file
      !-----------------------------------------------------------  
      if(line(1:6).eq.'#done ')then
        rewind(20+nopen+2)
        nopen=nopen+2 
        idooptns=2 
        k1=1
        nrdo=0
        goto 5
      endif
      !--------------do/done text to tmp---------------
      !Write text between #do and #done to tmp file 
      !------------------------------------------------
      if(idooptns.eq.1)then
        nrdox=nrdox+1  
        write(20+nopen+2,'(a)')line(1:imx)
        !print*,'CHECK #do write',20+nopen+2,'  ',line(1:80)
        goto 5
      endif
      !--------------do/done text expansion---------------------
      !Read from tmp file, expand and write to another tmp file 
      !--------------------------------------------------------
      if(idooptns.eq.2)then
        nrdo=nrdo+1
        i=1
        j=1
        do while(i.le.imx) 
          if(j.gt.len(linex))then
            if(line(i:i).ne.' ')stop'ERROR 27042021d'
            jx=j
            i=i+1
            j=j+1
          elseif(line(i:i+1).ne.'{{'
     .      .and.line(i+3:i+4).ne.'}}')then
            linex(j:j)=line(i:i)
            jx=j
            i=i+1
            j=j+1
          else
            read(line(i+2:i+2),*)k2
            idx=index(dooptnsargtext(k1,k2),' ')-1
            idx2=index(dooptnsargtext(k1,k2),'"')-1
            if(idx2.gt.0)idx=idx2
            linex(j:j+idx-1)=dooptnsargtext(k1,k2)(1:idx)
            jx=j+idx-1
            i=i+5
            j=j+idx
          endif
        enddo
        !print*,'CHECK #do write ',20+nopen-1,'  ',linex(1:jx)
        write(20+nopen-1,'(a)')linex(1:jx)
        if(nrdo.eq.nrdox)then
          k1=k1+1
          nrdo=0
          rewind(20+nopen)
          if(k1.eq.k1dooptns)then
            idooptns=0
            close(20+nopen)
            nopen=nopen-1
            rewind(20+nopen)
          endif 
        endif
        goto 5 
      endif

      !----------------------------------------------------------------
      !   #macro
      !----------------------------------------------------------------
      !macros, to be used  in optns as
      !#macro <name> <number of arguments>
      ! ... text with placeholders <<i>>
      !#endmacro
      ! ...
      !#expandmacro <name> <argument1> <argement2> ...
      !---------------------------------------------------------------- 
      !macronr .............. Number of macros
      !macroname(m) ......... Name of macro m
      !macrotext(m,l) ....... Text line l of macro m
      !macroargs(m) ......... Number of arguments of macro m
      !macroargtext(m,k) ... Value (text) of argument k of macro m
      !macrocurrent ......... Current macro number (to be used for read)
      !macrocurrentline ..... Next line nr of Current macro (for read)  
      !----------------------------------------------------------------
      if(line(1:11).eq.'#macroreset')then
        macronr=0
        if(ibeginwrite.eq.1)stop'#macroreset in write environment'
        goto 5
      endif
      !-------------#macro-------------------------------
      !Read and store the text between  #macro and #endmacro
      !---------------------------------------------------
      if(line(1:7).eq.'#macro ')then
        if(ibeginwrite.eq.1)stop'#macro in write environment'
        macronr=macronr+1
        if(macronr.gt.macromax)stop'ERROR 26032021'
        i=7
        do while(line(i:i).eq.' ')
          i=i+1
        enddo
        j=i  
        do while(line(j:j).ne.' ')
          j=j+1
        enddo
        !print*,'macro: ',line(i:j)
        macroname(macronr)(1:j-i+1)=line(i:j)
        macroname(macronr)(j-i+2:j-i+2)=' '
        i=j+1
        do while(line(i:i).eq.' ')
          i=i+1
        enddo
        j=i  
        do while(line(j:j).ne.' ')
          j=j+1
        enddo
        nrline=0
        read(line(i:j),*)macroargs(macronr) !number of arguments (show up as $1,$2...)
        if(macroargs(macronr).gt.macroargmax)stop'ERROR 26032021c'
        do while( line(1:9).ne.'#endmacro')
          read(ifopx,'(a)',end=9999)line
          do i=161,imax
            if(line(i:i).ne.' ')stop'#macro line too long'
          enddo  
          if(line(1:9).ne.'#endmacro')then
            nrline=nrline+1
            if(nrline.gt.macrolines)stop'ERROR 26032021b'
            macrotext(macronr,nrline)(1:160)=line(1:160)
            !write(ifmt,'(a,i3,$)')'CHECK macrotext:',nrline
            !write(ifmt,'(3x,a50,2i5)')line(1:50),i,j
          endif
        enddo
        macrotextlines(macronr)=nrline
        goto 5
      endif
      !----------#expandmacro---------------------
      !get macro name and arguments for expansion
      !------------------------------------------- 
      if(line(1:12).eq.'#expandmacro')then
        if(ibeginwrite.eq.1)stop'#expandmacro in write environment'
        i=13
        do while(line(i:i).eq.' ')
          i=i+1
        enddo
        j=i  
        do while(line(j:j).ne.' ')
          j=j+1
        enddo
        !print*,'expandmacro name: ',line(i:j)
        macro=1
        do while(macroname(macro)(1:j-i+1).ne.line(i:j))
          macro=macro+1
          if(macro.gt.macromax)stop'Macro unknown' !no match
        enddo 
        do k=1,macroargs(macro)
          i=j+1
          do while(line(i:i).eq.' ')
            i=i+1
          enddo
          j=i  
          do while(line(j:j).ne.' ')
            j=j+1
          enddo
          if(j-i+1.gt.macroargtextlength)then
            print*,'---------------------------------------------------' 
            print*,'ERROR macroargtext length > ', macroargtextlength
            print*,'string: ','-->'//line(i:j)//'<--'
            print*,'---------------------------------------------------' 
          endif
          macroargtext(macro,k)(1:j-i+1)=line(i:j)
          macroargtext(macro,k)(j-i+2:j-i+2)=' '
          !write(ifmt,'(a,$)')'CHECK expandmacro arg: '
          !write(ifmt,'(i5,3x,a)')k,line(i:j)
        enddo
        if(imx.gt.j)then
          print*,'-----------------------------------------------------' 
          print*,'ERROR too many #expandmacro arguments' 
          print*,'unmatched arguments: ',line(j+1:imx)
          print*,'-----------------------------------------------------' 
          stop
        endif 
        macrocurrent=macro
        macrocurrentline=1
        goto 5
      endif
      !-------------end expandmacro---------------
   4  continue
      !-------------expand macro text---------------
      !do the macro expansion line by line based on
      ! current macro and its next (not yet read) line
      !---------------------------------------------
      if(macrocurrent.ge.1.and.macrocurrentline.ge.1)then
      if(macrocurrentline.le.macrotextlines(macrocurrent))then
        nr=macrocurrentline
        i=1
        j=1
        do while(i.lt.160-4) 
          if(macrotext(macro,nr)(i:i+1).ne.'<<'
     .      .or.(macrotext(macro,nr)(i+3:i+4).ne.'>>'
     .           .and.macrotext(macro,nr)(i+4:i+5).ne.'>>'))then
            line(j:j)=macrotext(macro,nr)(i:i)
            i=i+1
            j=j+1
          elseif(macrotext(macro,nr)(i:i+1).eq.'<<'
     .      .and.macrotext(macro,nr)(i+3:i+4).eq.'>>')then !1-9
            read(macrotext(macro,nr)(i+2:i+2),*)k
            idx=index(macroargtext(macro,k),' ')-1
            line(j:j+idx-1)=macroargtext(macro,k)(1:idx)
            i=i+5
            j=j+idx
          elseif(macrotext(macro,nr)(i:i+1).eq.'<<'
     .      .and.macrotext(macro,nr)(i+4:i+5).eq.'>>')then !10-999
            read(macrotext(macro,nr)(i+2:i+3),*)k
            idx=index(macroargtext(macro,k),' ')-1
            line(j:j+idx-1)=macroargtext(macro,k)(1:idx)
            i=i+6
            j=j+idx
          else
            stop'ERROR 03072021'
          endif
        enddo 
        !write(ifmt,'(a,i3,$)')'macrotext:',nr
        !write(ifmt,'(3x,a80,2i5)')line(1:80),i,j
        macrocurrentline=macrocurrentline+1
        if(macrocurrentline.gt.macrotextlines(macrocurrent))then
          macrocurrent=0
          macrocurrentline=0
        endif
      endif
      endif
      !-------------end expand macro text---------------

      !----------------------------------------------------------------
      ! Write environment
      !   beginwrite ... endwrite
      ! to write sequences of text lines into histo file.
      ! The text lines may contain escape objects of the form 
      !        [[command]] 
      ! which allows to interpret "command" as normal EPOS
      ! command line input (frequently used commands are
      ! "histoweight" and "writearray 3") 
      !----------------------------------------------------------------
      if(line(1:10).eq.'beginwrite')then
        if(ibeginwrite.eq.1)stop'beginwrite in write environment'
        ibeginwrite=1
        do i=11,imax
          if(line(i:i).ne.' ')stop'character found after beginwrite'
        enddo  
        iwritecount=0
        goto 5
      endif
      if(line(1:8).eq.'endwrite')then
        if(ibeginwrite.eq.0)stop'endwrite outside write environment'
        ibeginwrite=0
        do i=9,imax
          if(line(i:i).ne.' ')stop'character found after endwrite'
        enddo  
        goto 5
      endif
      if(ibeginwrite.eq.1)then
        ii4=imax
        do while(ii4.gt.1.and.line(ii4:ii4).eq.' ')
          ii4=ii4-1
        enddo
        if(ii4.eq.1.and.line(1:1).eq.' ')then !empty line
          line(1:9)='write ~ ~' 
          goto 3
        endif
        ii1=1
        do while(ii1.lt.ii4.and.line(ii1:ii1).eq.' ')
          ii1=ii1+1
        enddo
        if(iwritecount.gt.0)ii1=min(ii1,ii1x)
        ii3=ii4-1
        do while(ii3.gt.ii1.and.line(ii3:ii3+1).ne.']]')
          ii3=ii3-1
        enddo
        ii2=ii1
        do while(ii2.lt.ii4-1.and.line(ii2:ii2+1).ne.'[[')
          ii2=ii2+1
        enddo
        if(ii4.lt.ii1)then!empty line
          line(1:9)='write ~ ~' 
          goto 3
        elseif(ii3.le.ii1)then!line without [[...]]  
          kk=7-ii1+1 !to insert "write ~"
          ii3=-2
        else!line with [[...]]  
          if(ii2+2.gt.ii3-1)stop'ERROR 20042021a cannot be'
          do i=ii2+2,ii3-1
            if(line(i:i+1).eq.'[['.or.line(i:i+1).eq.']]')
     .      stop'Only one [[...]] expression per line'
          enddo
          kk=14-ii1+1 !to insert twice "write ~"
          if(line(ii2:ii2+1).ne.'[[')stop'ERROR 20042021b cannot be'
          if(line(ii3:ii3+1).ne.']]')stop'ERROR 20042021c cannot be'
          line(ii2:ii2+1)='~ '
          line(ii3:ii3+1)='  '
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !  ii1        ii2-1    ii2+2       ii3-1  ii3+2        ii4
          !   |           |       |           |      |            |
          !   *****Text****     ~ ***Command***      *****Text*****
          ! |                                       |
          !write ~                                write ~           
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        endif
        kk4=ii4+kk
        if(kk4+1.gt.imax)stop'line in write environment too long'
        line(kk4+1:kk4+1)='~'
        i=kk4
        do while(i-kk.ge.ii3+2.and.i-kk.ge.ii1)
          line(i:i)=line(i-kk:i-kk)
          i=i-1
        enddo 
        if(i.lt.kk4)i=i+1
        line(i-7:i-1)='write ~'
        if(i-8.ge.1)line(i-8:i-8)=' ' 
        i=i-9
        if(i.ge.0)then
          kk=kk-7
          do while(i-kk.ge.ii1)
            line(i:i)=line(i-kk:i-kk)
            i=i-1
          enddo 
          i=i+1
          line(i-7:i-1)='write ~' 
        endif
        do i=1,kk4
          if(line(i:i+7).eq.'write ~~')then
            line(i:i+7)='        '
          endif
          if(line(i:i+8).eq.'write ~ ~')then
            line(i:i+8)='         '
          endif
        enddo 
        iwritecount=iwritecount+1
        if(iwritecount.eq.1)ii1x=ii1
      endif

      !----------
      !print line
      !----------
    3 continue
      if(iecho.eq.1.or.(nopen.ge.0.and.kcpopen.eq.1))then
       kmax=2
       do k=3,imax
       if(line(k:k).ne.' ')kmax=k
       enddo
      else
        kmax=0
      endif
      if(nopen.ge.0.and.kcpopen.eq.1)
     &  write(ifcp,'(a)')line(1:kmax)
      if(iecho.eq.1.and.line(1:3).ne.'not'.and.line(1:3).ne.'#if')
     &  write(ifmt,'(a)')line(1:kmax)
      i=0

      !----------------------------------------------------
      !find next non-empty i, request for next line in case 
      !----------------------------------------------------
    1 i=i+1
      if(i.gt.imax)goto 5
      if(line(i:i).eq.'!')goto 5
      do ne=1,mempty
      if(line(i:i).eq.empty(ne))goto 1
      enddo

      !----------------------------------------------
      !----------------------------------------------

      nbla=1
      mrk='  '
      mk=' '
      if(line(i:i).eq.'~')mk='~'
      if(line(i:i+1).eq.'"{')mrk='}"'
      if(line(i:i+1).eq.'""')mrk='""'
      if(mrk.ne.'  ')goto 10
      if(line(i:i).eq.'"')mk='"'
      if(mk.ne.' ')goto 8
      j=i-1
    6 j=j+1
      if(j.gt.imax)goto 7
      if(line(j:j).eq.'!')goto 7
      do ne=1,mempty
      if(line(j:j).eq.empty(ne))goto 7
      enddo
      goto 6

    8 continue
      if(i.ge.imax-1)stop'utword: make line shorter!!!         '
      i=i+1
      j=i
      if(line(j:j).eq.mk)stop'utword: empty string!!!           '
    9 j=j+1
      if(j.gt.imax)then                 !reach the end of the line
        j=j-nbla+2
        goto 7
      endif
      if(line(j:j).eq.' ')then
        nbla=nbla+1
      else
        nbla=2
      endif
      if(line(j:j).eq.mk)then
      line(i-1:i-1)=' '
      line(j:j)=' '
      goto 7
      endif
      goto 9

   10 continue
      if(i.ge.imax-3)stop'utword: make line shorter!!!!          '
      i=i+2
      j=i
      if(line(j:j+1).eq.mrk)stop'utword: empty string!!!!        '
   11 j=j+1
      if(j.gt.imax-1)then                 !reach the end of the line
        j=j-nbla+2
        goto 7
      endif
      if(line(j:j+1).eq.mrk)then
      line(i-2:i-1)='  '
      line(j:j+1)='  '
      goto 7
      endif
      if(line(j:j).eq.' ')then
        nbla=nbla+1
      else
        nbla=2
      endif
      goto 11

    7 j=j-1

      if(line(i:i+1).eq.'::'.and.line(j-1:j).eq.'::')then
        line(i:i+1)='  '
        line(j-1:j)='  '
        i=i+2
        j=j-2
        goto 7777
      endif

      if(iqu.ne.-1)call expand(line,i,j)
      do k=i,j
        if(k.gt.1)then
          if(line(k-1:k).eq.'::')then
            line(k-1:k)='  '
            j=k-2
            goto 12
          endif
        endif
      enddo
 12   continue
      if(iqu.ne.-1)call expand(line,i,j)

7777  continue
      return

9999  close(ifopx)
      nopen=nopen-1
      if(nopen.eq.0.and.iprmpt.eq.-1)iprmpt=1
      goto 5
      end

      subroutine expand(line,i,j)
      include "epos.inc"
      character line*(*)
      if(ndefine.le.0)return
      imax=len(line) 
      imx=imax
      do while(imx.gt.2.and.line(imx:imx).eq.' ')
        imx=imx-1
      enddo
      imx=imx+1
      j0=0
      jj=0
      do ndf=1,ndefine
        l1=l1define(ndf)
        l2=l2define(ndf)
        do i0=i,j+1-l1
          if(line(i0:i0-1+l1).eq.w1define(ndf)(1:l1))then
           if(line(i0-1:i0-1).ne.'\'.and.line(i0-1:i0-1).ne.'/')then
            jj=1 
            if(l2.eq.l1)then
              line(i0:i0-1+l1)=w2define(ndf)(1:l2)
            elseif(l2.lt.l1)then
              line(i0:i0+l2-1)=w2define(ndf)(1:l2)
              do k=i0+l2,i0-1+l1
                line(k:k)=' '
              enddo
            elseif(l2.gt.l1)then
              ier=0
              do k=i0+l1,i0+l2
                if(line(k:k).ne.' ')ier=1
              enddo
              if(ier.eq.1)then
                !write(ifmt,'(9a)')('-----------------',m=1,4)
                !write(ifmt,*)'uti.f:expand'
                !write(ifmt,*)'  i j i0 = ',i,j,i0
                !write(ifmt,*)'  no space for `define` replacement.'
                !write(ifmt,*)'  ',line(i0:i0-1+l2+1),l1
                !write(ifmt,*)'  ',w2define(ndf)(1:l2),l2
                k1=l2-l1
                do k=imx,i0+l1,-1
                  line(k+k1:k+k1)=line(k:k)
                  line(k:k)=' '
                enddo
                !write(ifmt,*)'  ',line(i0:i0-1+l2+1+k1)
                !write(ifmt,'(9a)')('-----------------',m=1,4)
                !stop'exit'
              endif
              line(i0:i0+l2-1)=w2define(ndf)(1:l2)
              j=i0+l2-1
            endif
           elseif(line(i0-1:i0-1).eq.'\')then
            line(i0-1:i0-1)='/'
           elseif(line(i0-1:i0-1).eq.'/')then
            line(i0-1:i0-1)=' '
           endif
          endif
        enddo
      enddo
      if(jj.eq.1)then
        do k=i,j
          if(line(k:k).ne.' ')j0=j
        enddo
        j=j0
      endif
      end

c--------------------------------------------------------------------
      subroutine dollartext(text,line,i,j)
c--------------------------------------------------------------------
      character line*(*)
      character text*(*)
      character*30 hdtext(8)
      parameter(leng=30)
      common/chdtext/hdtext
      character cfmt*5
      cfmt='(a  )'
      write(cfmt(3:4),'(i2)')leng
      n=index(text,';')-1
      n1=n-1
      do l=i,j-n1
      if(line(l:l+n1).eq.text(1:n))then
        read(text(n:n),'(i1)')iii
        do jj=l+n,l+leng-1
          if(line(jj:jj).ne.' ')stop'###### ERROR 26052014 #####'
        enddo
        write(line(l:l+leng-1),cfmt)hdtext(iii)
        return
      endif
      enddo
      end

c--------------------------------------------------------------------
      subroutine dollar(text,yield,line,i,j)
c--------------------------------------------------------------------
      character line*(*)
      character text*(*)
      character*9 cfmt
      n=index(text,';')-1
      n1=n-1
      do l=i,j-n1
      if(line(l:l+n1).eq.text(1:n))then
        nd=0
        do while(line(l+n1+nd+1:l+n1+nd+1).eq.' ')
          nd=nd+1
        enddo
        cfmt='(f_._,_x)'
        ayield=abs(yield)
        if(ayield.eq.0.)then
          write(line(l:l+n1),'(a)')'0'
        else
          fk=log10(ayield)
          k=1+fk
          if(float(k).eq.1+fk)k=k+1
          m2=4-k
          m2=max(m2,0)
          m2=min(m2,n-1)
          m1=1
          if(k.gt.0)m1=k
          nn=n
          do ik=1,nd
            if(nn-(m1+1+m2).lt.0)nn=nn+1
          enddo
          mx=nn-(m1+1+m2)
          write(cfmt(3:3),'(i1)')m1+1+m2
          write(cfmt(5:5),'(i1)')m2
          write(cfmt(7:7),'(i1)')mx
          if(mx.eq.0)write(cfmt(6:8),'(a)')'   '
             !write(*,'(40x,5(a,10x))')'nn',' k','m2','m1','mx'
             !print*,yield,log10(ayield),nn,k,m2,m1,mx,cfmt,fk
             !print*,' '
          if(mx.ge.0)then
            write(line(l:l+nn-1),cfmt)yield
          else
            do ii=l,l+nn-1
            write(line(ii:ii),'(a)')'?'
            enddo
          endif
        endif
      endif
      enddo
      end

c--------------------------------------------------------------------
      subroutine utworn(line,j,ne)
c--------------------------------------------------------------------
c  returns number ne of nonempty characters of line(j+1:1000)
c--------------------------------------------------------------------
      character line*(*)
      imax=len(line)
      ne=0
      do l=j+1,imax
      if(line(l:l).ne.' ')ne=ne+1
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine getairmol(iz,ia)
c-----------------------------------------------------------------------
      include 'epos.inc'
      i=0
      r=rangen()
      do while(r.gt.0.)  ! choose air-molecule
        i=i+1
        r=r-airwnxs(i)
      enddo
      iz = nint(airznxs(i))
      ia = nint(airanxs(i))
      end

c----------------------------------------------------------------------

      subroutine factoriel

c----------------------------------------------------------------------
c tabulation of fctrl(n)=n!, facto(n)=1/n! and utgamtab(x) for x=0 to 50
c----------------------------------------------------------------------
      include 'epos.incems'
      double precision utgamtab,utgam,x
      common/gamtab/utgamtab(10000)

      nfctrl=100
      fctrl(0)=1.D0
      facto(0)=1.D0
      do i=1,min(npommx,nfctrl)
        fctrl(i)=fctrl(i-1)*dble(i)
        facto(i)=1.d0/fctrl(i)
      enddo

      do k=1,10000
        x=dble(k)/100.d0
        utgamtab(k)=utgam(x)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine fremnu(amin,ca,cb,ca0,cb0,ic1,ic2,ic3,ic4)
c-----------------------------------------------------------------------
      common/hadr2/iomodl,idproj,idtarg,wexcit
      real pnll,ptq
      common/hadr1/pnll,ptq,exmass,cutmss,wproj,wtarg
      parameter (nflav=6)
      integer ic(2),jc(nflav,2)
      ic(1)=ca
      ic(2)=cb
      call iddeco(ic,jc)
      keu=jc(1,1)-jc(1,2)
      ked=jc(2,1)-jc(2,2)
      kes=jc(3,1)-jc(3,2)
      kec=jc(4,1)-jc(4,2)
      keb=jc(5,1)-jc(5,2)
      ket=jc(6,1)-jc(6,2)
      amin=utamnu(keu,ked,kes,kec,keb,ket,5)  !???4=2mults, 5=1mult
        if(ca-ca0.eq.0..and.cb-cb0.eq.0..and.rangen().gt.wexcit)then
      ic3=0
      ic4=0
      ic1=ca
      ic2=cb
       else
      amin=amin+exmass
      n=0
      do i=1,4
      do j=1,2
      n=n+jc(i,j)
      enddo
      enddo
      k=1+rangen()*n
      do i=1,4
      do j=1,2
      k=k-jc(i,j)
      if(k.le.0)goto 1
      enddo
      enddo
1     if(j.eq.1)then
      ic3=10**(6-i)
      ic4=0
      else
      ic3=0
      ic4=10**(6-i)
      endif
      ic1=int(ca)-ic3
      ic2=int(cb)-ic4
        endif
      return
      end


c-----------------------------------------------------------------------
      function fremnux(jc)
c-----------------------------------------------------------------------
      real bidiba,amhdibar
      integer iodiba
      common /ems12/bidiba,amhdibar,iodiba  ! defaut iodiba=0. if iodiba=1, study H-Dibaryon
      real pnll,ptq
      common/hadr1/pnll,ptq,exmass,cutmss,wproj,wtarg
      parameter (nflav=6)
      integer jc(nflav,2)!,ic(2)
c      ic(1)=210000
c      ic(2)=0
c      call iddeco(ic,jc)
      if(iodiba.eq.1)then
        if((jc(1,1).eq.2.and.jc(2,1).eq.2.and.jc(3,1).eq.2
     .          .and.jc(1,2)+jc(2,2)+jc(3,2).eq.0)
     . .or.(jc(1,2).eq.2.and.jc(2,2).eq.2.and.jc(3,2).eq.2
     .       .and.jc(1,1)+jc(2,1)+jc(3,1).eq.0))then
          fremnux=amhdibar-bidiba
          return
        endif
      endif
      keu=jc(1,1)-jc(1,2)
      ked=jc(2,1)-jc(2,2)
      kes=jc(3,1)-jc(3,2)
      kec=jc(4,1)-jc(4,2)
      keb=jc(5,1)-jc(5,2)
      ket=jc(6,1)-jc(6,2)
      fremnux=utamnu(keu,ked,kes,kec,keb,ket,4)+exmass  !???4=2mults, 5=1mult
      return
      end

c-----------------------------------------------------------------------
      function fremnux2(jc)
c-----------------------------------------------------------------------
      real pnll,ptq
      common/hadr1/pnll,ptq,exmass,cutmss,wproj,wtarg
      parameter (nflav=6)
      integer jc(nflav,2)!,ic(2)
c      ic(1)=210000
c      ic(2)=0
c      call iddeco(ic,jc)
      keu=jc(1,1)-jc(1,2)
      ked=jc(2,1)-jc(2,2)
      kes=jc(3,1)-jc(3,2)
      kec=jc(4,1)-jc(4,2)
      keb=jc(5,1)-jc(5,2)
      ket=jc(6,1)-jc(6,2)
      fremnux2=utamnu(keu,ked,kes,kec,keb,ket,4) !+exmass  !???4=2mults, 5=1mult
      return
      end

c-----------------------------------------------------------------------
      function fremnux3(jci)
c-----------------------------------------------------------------------
c minimum mass from ic counting all quarks
c-----------------------------------------------------------------------
      include 'epos.inc'
      integer jc(nflav,2),jci(nflav,2)!,ic(2)
c      ic(1)=210000
c      ic(2)=0
c      print *,'start',ic
      fremnux3=0.
      do j=1,2
        do i=1,nflav
          jc(i,j)=jci(i,j)
        enddo
      enddo
c      call iddeco(ic,jc)
      call idquacjc(jc,nqua,naqu)
        do ii=1,2
      if(ii.eq.1)then
        nqu=nqua
      else
        nqu=naqu
      endif
      if(nqu.ge.3)then
        do while(jc(3,ii).ne.0.and.nqu.ge.3)  !count baryons with s quark
          jc(3,ii)=jc(3,ii)-1
          if(jc(3,ii).gt.0)then
            jc(3,ii)=jc(3,ii)-1
            if(jc(3,ii).gt.0)then
              jc(3,ii)=jc(3,ii)-1
              fremnux3=fremnux3+asuhax(4)
            elseif(jc(2,ii).gt.0)then
              jc(2,ii)=jc(2,ii)-1
              fremnux3=fremnux3+asuhax(4)
            elseif(jc(1,ii).gt.0)then
              jc(1,ii)=jc(1,ii)-1
              fremnux3=fremnux3+asuhax(4)
            endif
          elseif(jc(2,ii).gt.0)then
            jc(2,ii)=jc(2,ii)-1
            if(jc(1,ii).gt.0)then
              jc(1,ii)=jc(1,ii)-1
              fremnux3=fremnux3+asuhax(3)
            elseif(jc(2,ii).gt.0)then
              jc(2,ii)=jc(2,ii)-1
              fremnux3=fremnux3+asuhax(3)
            endif
          elseif(jc(1,ii).gt.0)then
            jc(1,ii)=jc(1,ii)-2
            fremnux3=fremnux3+asuhay(3)
          endif
          nqu=nqu-3
        enddo
        do while(jc(2,ii).ne.0.and.nqu.ge.3)  !count baryons with d quark
          jc(2,ii)=jc(2,ii)-1
          if(jc(1,ii).gt.0)then
            jc(1,ii)=jc(1,ii)-1
            if(jc(2,ii).gt.0)then
              jc(2,ii)=jc(2,ii)-1
              fremnux3=fremnux3+asuhay(2)
            elseif(jc(1,ii).gt.0)then
              jc(1,ii)=jc(1,ii)-1
              fremnux3=fremnux3+asuhay(2)
            endif
          elseif(jc(2,ii).gt.0)then
            jc(2,ii)=jc(2,ii)-2
            fremnux3=fremnux3+asuhay(3)
          endif
          nqu=nqu-3
        enddo
        do while(jc(1,ii).ne.0.and.nqu.ge.3)  !count baryons with s quark
          jc(1,ii)=jc(1,ii)-3
          fremnux3=fremnux3+asuhay(3)
          nqu=nqu-3
        enddo
        if(ii.eq.1)then
          nqua=nqu
        else
          naqu=nqu
        endif
      endif
c      print *,ii,nqua,naqu,jc,fremnux3
      enddo
      if(nqua+naqu.ne.0)then
      do while(jc(3,1).ne.0)    !count mesons with s quark
        jc(3,1)=jc(3,1)-1
        if(jc(3,2).gt.0)then
          jc(3,2)=jc(3,2)-1
          fremnux3=fremnux3+asuhax(6)
        elseif(jc(2,2).gt.0)then
          jc(2,2)=jc(2,2)-1
          fremnux3=fremnux3+asuhay(6)
        elseif(jc(1,2).gt.0)then
          jc(1,2)=jc(1,2)-1
          fremnux3=fremnux3+asuhay(6)
        endif
      enddo
      do while(jc(2,1).ne.0)    !count mesons with d quark
        jc(2,1)=jc(2,1)-1
        if(jc(2,2).gt.0)then
          jc(2,2)=jc(2,2)-1
          fremnux3=fremnux3+asuhay(5)
        elseif(jc(1,2).gt.0)then
          jc(1,2)=jc(1,2)-1
          fremnux3=fremnux3+asuhay(5)
        endif
      enddo
      do while(jc(1,1).ne.0)    !count mesons with s quark
        jc(1,1)=jc(1,1)-1
        if(jc(1,2).gt.0)then
          jc(1,2)=jc(1,2)-1
          fremnux3=fremnux3+asuhay(5)
        endif
      enddo
      endif
c      fremnux3=fremnux3+0.5
c      print *,'stop',nqua,naqu,fremnux3

      return
      end

c-----------------------------------------------------------------------
      subroutine fremnx(ammax,amin,sm,ic3,ic4,iret)
c-----------------------------------------------------------------------
      common/psar9/ alpr
      include 'epos.inc'
      iret=0
      if(ic3.eq.0.and.ic4.eq.0)then
        if(ammax.lt.amin**2)then
          iret=1
          return
        endif
        sm=amin**2
      else
c       ammax1=min(ammax,(engy/4.)**2)
        ammax1=ammax
        if(ammax1.lt.amin**2)then
          iret=1
          return
        endif
        if(alpr.eq.-1.)then
          sm=amin**2*(ammax1/amin**2)**rangen()
        else
          sm=amin**2*(1.+((ammax1/amin**2)**(1.+alpr)-1.)
     *    *rangen())**(1./(1.+alpr))
        endif
      endif
      return
      end

      SUBROUTINE gaulag(x,w,n,alf)
      INTEGER n,MAXIT
      REAL alf,w(n),x(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.D-14,MAXIT=10)
CU    USES gammln
      INTEGER i,its,j
      REAL ai,gammln
      DOUBLE PRECISION p1,p2,p3,pp,z,z1
      z=0.
      do 13 i=1,n
        if(i.eq.1)then
          z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
        else if(i.eq.2)then
          z=z+(15.+6.25*alf)/(1.+.9*alf+2.5*n)
        else
          ai=i-2
          z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.+3.5*ai))*
     *(z-x(i-2))/(1.+.3*alf)
        endif
        do 12 its=1,MAXIT
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
11        continue
          pp=(n*p1-(n+alf)*p2)/z
          z1=z
          z=z1-p1/pp
          if(abs(z-z1).le.EPS)goto 1
12      continue
        call utstop("too many iterations in gaulag")
1       x(i)=z
        w(i)=-exp(gammln(alf+n)-gammln(float(n)))/(pp*n*p2)
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 4+1$!].

c-----------------------------------------------------------------------
      FUNCTION gammln(xx)
c-----------------------------------------------------------------------
c log(Gamma(x)) single precision
c-----------------------------------------------------------------------
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END


c-----------------------------------------------------------------------
      double precision FUNCTION dgammln(x)
c-----------------------------------------------------------------------
c log(Gamma(x)) double precision
c-----------------------------------------------------------------------
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      dgammln=tmp+log(stp*ser/x)
      return
      END

C     (C) Copr. 1986-92 Numerical Recipes Software 4+1$!].

c-----------------------------------------------------------------------
      double precision function digamma(ex)
c-----------------------------------------------------------------------
c   Digamma function (derivative of logarithm of Gamma function)
c   from "An Atlas of functions" (Spanier and Oldham ch. 44, p. 423)
c   T. Pierog, 22/10/2010 (from a previous work 12/02/1999)
c-----------------------------------------------------------------------
      implicit none
      double precision ex,x,g,f

c Initialization

      g=0d0
      x=ex
      digamma=0d0

c Digamma function calculation

      do while (x.lt.5d0.and.abs(x).gt.1d-15)
        g=g+1d0/x
        x=x+1d0
      enddo

      if(abs(x).gt.1d-15)then
        f=1d0+(0.46d0/(x*x)-1d0)/(10d0*x*x)
        f=-((f/(6d0*x)+1d0)/(2d0*x))+log(x)-g

c get out function value at point x

        digamma=f

      else

c infinite value for digamma
        call utstop('infinity in digamma function!&')

      endif

      return

      end


c-----------------------------------------------------------------------
      double precision FUNCTION betai(a,b,x)
c-----------------------------------------------------------------------
c   Incomplete beta function I_x(a, b)=B_x(a,b)/B(a,b)
c   with B_x(a,b)=int from 0 to x of (dt t**(a-1)*(1-t)**(b-1))
c   from "NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING"
c   (ISBN 0-521-43064-X)
c   T. Pierog, 20/03/2017
c-----------------------------------------------------------------------
      double precision a,b,x
C USES betacf,dgammln
      double precision bt,betacf,dgammln
      if(x.lt.0.d0.or.x.gt.1.d0)then
        print *,a,b,x
        call utstop('bad argument x in betai&')
      endif
      if(x.eq.0.d0.or.x.eq.1.d0)then
        bt=0.d0
      else          !Factors in front of the continued fraction.
        bt=exp(dgammln(a+b)-dgammln(a)-dgammln(b)
     *       +a*log(x)+b*log(1.-x))
      endif
      if(x.lt.(a+1.d0)/(a+b+2.d0))then   !Use continued fraction directly.
        betai=bt*betacf(a,b,x)/a
        return
      else
        betai=1.d0-bt*betacf(b,a,1.d0-x)/b     !Use continued fraction after making the symmetry transformation.
        return
      endif
      END

c-----------------------------------------------------------------------
      double precision function betacf(a,b,x)
c-----------------------------------------------------------------------
c   continued fraction evaluation routine
c   from "NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING"
c   (ISBN 0-521-43064-X)
c   T. Pierog, 20/03/2017
c-----------------------------------------------------------------------
      INTEGER MAXIT
      double precision a,b,x,EPS,FPMIN
      PARAMETER (MAXIT=100,EPS=3.d-7,FPMIN=1.d-30)
c Used by betai: Evaluates continued fraction for incomplete beta function by modified Lentz's method (Â§5.2).
      INTEGER m,m2
      double precision aa,c,d,del,h,qab,qam,qap
      qab=a+b       !These q's will be used in factors that occur in the
      qap=a+1.d0      !coefficients (6.4.6).
      qam=a-1.d0
      c=1.d0          !First step of Lentz's method.
      d=1.d0-qab*x/qap
      if(abs(d).lt.FPMIN)d=FPMIN
      d=1.d0/d
      h=d
      do m=1,MAXIT
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.d0+aa*d     !One step (the even one) of the recurrence.
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.d0+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        h=h*d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.d0+aa*d     !Next step of the recurrence (the odd one).
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.d0+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=d*c
        h=h*del
        if(abs(del-1.d0).lt.EPS)goto 1 !Are we done?
      enddo
      call utstop("a or b too big, or MAXIT too small in betacf&")
 1    betacf=h
      return
      END


c-----------------------------------------------------------------------
      double precision function utdble(x)
c-----------------------------------------------------------------------
c     define x as double precision with 0 in additionnal decimals
c-----------------------------------------------------------------------
      implicit none
      real x,xa
      integer ie,im,ip
      parameter(ip=6) !7 digit precision

      xa=abs(x)
      if(xa.gt.10.**(-ip))then
        ie=ip-int(log10(xa))
        im=nint(xa*10.**ie)
        utdble=sign(dble(im)/10d0**ie,dble(x))
      else
        utdble=0d0
      endif

c      print *,x,dble(x),im,utdble
      
      return
      end
      

c ------- modified Bessel function K_n(x) for n.ge.2
      
      FUNCTION sbessk(n,x)
      INTEGER n
      REAL sbessk,x
CU    USES sbessk0,sbessk1
      INTEGER j
      REAL bk,bkm,bkp,tox,sbessk0,sbessk1
      if (n.lt.2) stop 'bad argument n in sbessk'
      tox=2.0/x
      bkm=sbessk0(x)
      bk=sbessk1(x)
      do 11 j=1,n-1
        bkp=bkm+j*tox*bk
        bkm=bk
        bk=bkp
11    continue
      sbessk=bk
      return
      END

c-----------------------------------------------------------------------
      subroutine setInput(line,i,j,nopen,iprmpt,iopcnt)
c-----------------------------------------------------------------------
      character line*1000, cfmt*4
      common/ciexhd/iexhd
      call utword(line,i,j,0)
      if(nopen.ge.0)then
       nopen=nopen+1
       if(nopen.gt.9)stop'too many nested input commands'
c       print*,line(i:j)
       open(unit=20+nopen,file=line(i:j),status='old')
       if(iprmpt.eq.1)iprmpt=-1
       if(line(j-5:j).eq.'.optns')then
          n=0
        if(line(j-7:j-7).eq.'-')then
          n=1
        elseif(line(j-8:j-8).eq.'-')then
          n=2
        elseif(line(j-9:j-9).eq.'-')then
          n=3
        elseif(line(j-10:j-10).eq.'-')then
          n=4
        elseif(line(j-11:j-11).eq.'-')then
          n=5
        elseif(line(j-12:j-12).eq.'-')then
          n=6
        endif
        if(n.gt.0.and.line(j-8-n:j-6-n).eq.'/z-')n=0
        iopcnt=0
        if(n.gt.0)then
          cfmt='(i )'
          write(cfmt(3:3),'(i1)')n
          jj=j-5-n
          read(line(jj:j-6),cfmt)iopcnt
        endif
        do k=i,j-4
          if(line(k:k+4).eq.'ppihd')iexhd=1
        enddo
       endif
      endif
      end

c-----------------------------------------------------------------------
      subroutine setXinput(line,i,j,nopen,iprmpt,fnnx,nfnnx)
c-----------------------------------------------------------------------
      character line*1000, fnnx*500
      common/ciexhd/iexhd
      call utword(line,i,j,0)
      if(nopen.ge.0)then
       nopen=nopen+1
       if(nopen.gt.9)stop'too many nested input commands'
       open(unit=20+nopen,file=fnnx(1:nfnnx)//line(i:j)
     .   ,status='old')
       if(iprmpt.eq.1)iprmpt=-1
      endif
      end
   
c-----------------------------------------------------------------------
      subroutine setIf1(line,i,j) !#if1 definition
c-----------------------------------------------------------------------
      include "epos.inc"
      character*1000 line
      character cext1*10
      common/ccext1/cext1
      character*2 celse(20)
      integer nelse,inotelse
      common/ccelse/nelse,inotelse,celse
      lcext1=index(cext1,' ')-1
      inot=0
      call utword(line,i,j,0)
      if(line(i:j).eq.'#not')then
        inot=1
        call utword(line,i,j,0)
      endif
      nif1=1
      if(line(i:i).eq.'#')then
        read(line(i+1:j),*)nif1
        call utword(line,i,j,0)
      endif
      nelse=nif1
      iskip=2 !skip
      celse(1)=line(i:j)
      if(line(i:j).eq.cext1(1:lcext1))iskip=1 !noskip
      do nuw=2,nif1
        call utword(line,i,j,0)
        celse(nuw)=line(i:j)
        if(line(i:j).eq.cext1(1:lcext1))iskip=1 !noskip
      enddo
      if(inot.eq.1)iskip=3-iskip !invert
      if(iskip.eq.2)then !skip
        iechosave=iecho
        iecho=0
        do while(line(i:j).ne.'#else1'
     .   .and.line(i:j).ne.'#fi'.and.line(i:j).ne.'#fi1')
          call utword(line,i,j,-1)
        enddo
        iecho=iechosave
        if(line(i:j).eq.'#else1')j=max(0,j-7)
      else
        kmax=2
        do k=3,1000
        if(line(k:k).ne.' ')kmax=k
        enddo
        if(iecho.eq.1)write(ifmt,'(a)')line(1:kmax)
      endif
      end

c-----------------------------------------------------------------------
      subroutine setElse1(line,i,j) !#else1 definition
c-----------------------------------------------------------------------
      include "epos.inc"
      character*1000 line
      character cext1*10
      common/ccext1/cext1
      character*2 celse(20)
      integer nelse,inotelse
      common/ccelse/nelse,inotelse,celse
      lcext1=index(cext1,' ')-1
      inot=inotelse
      nif1=nelse
      iskip=2 !skip
      if(celse(1).eq.cext1(1:lcext1))iskip=1 !noskip
      do nuw=2,nif1
        if(celse(nuw).eq.cext1(1:lcext1))iskip=1 !noskip
      enddo
      if(inot.eq.0)iskip=3-iskip !invert
      if(iskip.eq.2)then !skip
        iechosave=iecho
        iecho=0
        do while(line(i:j).ne.'#fi'.and.line(i:j).ne.'#fi1')
          call utword(line,i,j,-1)
        enddo
        iecho=iechosave
      else
        kmax=2
        do k=3,1000
        if(line(k:k).ne.' ')kmax=k
        enddo
        if(iecho.eq.1)write(ifmt,'(a)')line(1:kmax)
      endif
      end

c-----------------------------------------------------------------------
      subroutine setCentralityClass(line,i,j,nopen)
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/hydr4/zclass(5,100),izmode,maxsurf,iofrout,jzmode(7)
      common/ccc20/icc20
      character line*1000
      call utword(line,i,j,0)
      if(line(i:j).eq.'cc20')then
        icc20=1
        return
      endif
      read(line(i:j),*)val
      ival=nint(val)
      jval=ival
      call utword(line,i,j,0)
      if(line(i:j).eq.'-')then
        call utword(line,i,j,0)
        read(line(i:j),*)val
        jval=nint(val)
        call utword(line,i,j,0)
      endif
      read(line(i:j),*)val1
      call utword(line,i,j,0)
      read(line(i:j),*)val2
      call utword(line,i,j,0)
      read(line(i:j),*)val3
      if(izmode.eq.2)then
        call utword(line,i,j,0)
        read(line(i:j),*)val4
        call utword(line,i,j,0)
        read(line(i:j),*)val5
      endif
      do i=ival,jval
      if(nopen.ne.-1)then   !only first read
        if(zclass(3,i).gt.0.001)then
          write(ifmt,'(a,i3)') 'ERROR 24022013 zclass3 in use, i=',i
          stop'24022013'
        endif
      endif
      zclass(1,i)=val1
      zclass(2,i)=val2
      zclass(3,i)=val3
      if(izmode.eq.2)then
        zclass(4,i)=val4
        zclass(5,i)=val5
      endif
      enddo
      end

c-----------------------------------------------------------------------
      subroutine setCentrality(line,i,j,iopcnt,ifmt
     .               ,icentrality,jcentrality,ffac,imax,ival,ishxxx)
c-----------------------------------------------------------------------
      common/hydr4/zclass(5,100),izmode,maxsurf,iofrout,jzmode(7)
      character line*1000
      
      if(line(i:j).eq.'within')then
        if(iopcnt.eq.0)then
          write(ifmt,'(//a/a/a)')
     .     'ERROR in aread:'
     .    ,'centrality within  and  iopcnt = 0  incompatible'
     .    ,'centrality within  requires batch mode!'
          stop'in aread \n\n '
        endif
        call utword(line,i,j,0)
        if(line(i:j).ne.'{')stop'\n\n ERROR 19112011\n\n'
        k=0
        i0=i
        j0=j
        do while(line(i:j).ne.'}')
          call utword(line,i,j,0)
          k=k+1
        enddo
        ival=k-1
        i=i0
        j=j0
        ij=1+mod((iopcnt-1),ival)
        do k=1,ival
          call utword(line,i,j,0)
          read(line(i:j),*)val1
          if(k.eq.ij)then     
            icentrality=nint(val1)
          endif
        enddo
        call utword(line,i,j,0)
        if(line(i:j).ne.'}')stop'\n\n ERROR 19112011b\n\n'
        jcentrality=icentrality
        ffac=zclass(3,icentrality)
        imax=1
        if(ishxxx.eq.1)
     .    write(ifmt,'(a,i3)')'centrality set to',icentrality
      else
        read(line(i:j),*)val
        ival=nint(val)
        if(ival.lt.0)then
          jcentrality=ival
          if(iopcnt.eq.0)then
            write(ifmt,'(//a/a/a)')
     .       'ERROR in aread:'
     .      ,'centrality < 0  and  iopcnt = 0  incompatible'
     .      ,'centrality < 0  requires batch mode!'
            stop'in aread \n\n '
          endif
          imax=0
          do i=1,100
           !print*,'1111111111111111',i,zclass(3,i)
           if(zclass(3,i).gt.1e-5)imax=i
          enddo
          icentrality=1+mod((iopcnt-1)/(-jcentrality),imax)
          ffac=zclass(3,icentrality)
          if(ishxxx.eq.1)
     .      write(ifmt,'(a,i3)')'centrality set to',icentrality
        else
          icentrality=ival
          jcentrality=ival
          ffac=1.00000
          imax=1
        endif
      endif

      end

c----------------------------------------------------------------------
      subroutine getMcentr  ! -> mcentrf() -> mmxcentrf()
c----------------------------------------------------------------------
      include "epos.inc"
      common/mcen/mcentr,mmxcentr
      common/nucl3/phi,bimp
      bimevtsave=bimevt
      bimevt=bimp
      call getcentr(3,mcentr,mmxcentr)  ! 3 -> M centrality
      bimevt=bimevtsave
      end
      integer function mcentrf()
      common/mcen/mcentr,mmxcentr
      mcentrf=mcentr
      end
      integer function mmxcentrf()
      common/mcen/mcentr,mmxcentr
      mmxcentrf=mmxcentr
      end

c----------------------------------------------------------------------
      subroutine getcentr(i,icentr,imxcentr) !i: 1 -> I, 2 -> J, 3 -> M
c----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incho"
      include "epos.incems"
      kk=0
      kkx=0
      x=centrVar(jzmode(i))
      !write(ifmt,*)'++++',i,jzmode(i),nprt(1),ikoevt
      if(jtable(i).gt.0)then
        do k=1,nrclass(jtable(i))
        xmin=rclass(jtable(i),1,k)
        xmax=rclass(jtable(i),2,k)
        !write(ifmt,*)'++++',i,xmin,xmax,x
        if(x.ge.xmin.and.x.le.xmax)kk=k
        enddo
        kkx=nrclass(jtable(i))
      endif
      icentr=kk
      imxcentr=kkx
      !write(ifmt,*)'++++',i,icentr
      end

c----------------------------------------------------------------------
      function centrVar(j)
c----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incho"
      include "epos.incems"
      if(j.eq.0)then
        x=0
      elseif(j.eq.1)then !C1
        x=bimevt
      elseif(j.eq.2)then !C2
        !x=nprt(1)
        !x=ikoevt
        x=0
        do i=1,nptl
         call getistptl(i,istx)
         if(istx.eq.30.or.istx.eq.31)x=x+1
        enddo
      elseif(j.eq.3)then
        x=nglevt
      elseif(j.eq.4)then
        x=segevt
      elseif(mod(j,10).eq.5)then
        if(iposi.lt.2)stop'in centrVar: wrong iposi\n\n'
        k=j/10
        nbi=nrclass(k)
        if(nbi.ne.3)stop'in centrVar: wrong nbi\n\n'
        eta1=rclass(k,1,1)
        eta2=rclass(k,2,1)
        pt1=rclass(k,1,2)
        pt2=rclass(k,2,2)
        amult=rclass(k,1,3)
        adivi=rclass(k,2,3)
        mul=0
        do i=maproj+matarg+1,nptl
          if(istptl(i).eq.0)then
            pt=pptl(1,i)**2+pptl(2,i)**2
            pp=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(3,i)**2)
            if(pt.gt.0.)then
              pt=sqrt(pt)
              eta=sign(1.,pptl(3,i))*alog((pp+abs(pptl(3,i)))/pt)
            else
              eta=1000.
            endif
            if(abs(idptl(i)).ge.100
     $     .and.abs(idptl(i)).lt.10000)then
              call idchrg( 3 ,idptl(i),ch)
              if(abs(ch).gt.0.1)then
                if( eta.ge.eta1.and.eta.le.eta2
     *         .and. pt.ge. pt1.and. pt.le. pt2) mul=mul+1
              endif
            endif
          endif
        enddo
        x=mul*amult/adivi
      elseif(mod(j,10).eq.6)then !number of Pomerons
        npom=0
        do n=1,nptl
        if(istptl(n).eq.30.or.istptl(n).eq.31)npom=npom+1
        enddo
        x=npom
      elseif(mod(j,10).eq.7)then !vzero
        x=fmux(2.8, 5.1, 0., 100., 1., 1., 1., -0.468,0.,0.,0.) 
        x=x*0.10  !rescaled!!! to fit into 0-40 range, as Pomerons
      else
        stop'in centrVar: wrong j    '
      endif
      centrVar =x
      end

c----------------------------------------------------------------------
      subroutine getncenthy
c----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incho"
      imax=0
      do i=1,100
        if(zclass(3,i).gt.1e-5)imax=i
      enddo
      ncenthy=imax
      end

c----------------------------------------------------------------------
      subroutine getJcentr
c----------------------------------------------------------------------
      common/jcen/jcentr,jmxcentr
      call getcentr(1,jcentr,jmxcentr)
      end

c----------------------------------------------------------------------
      subroutine getKcentr
c----------------------------------------------------------------------
      common/kcen/kcentr,kmxcentr
      call getcentr(2,kcentr,kmxcentr)
      end

c----------------------------------------------------------------------
      subroutine setcentrVar(x)
c----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incho"
      include "epos.incems"
      common/cirtfile/irtfile
      common/ccc20/icc20  /ciext4/iext4
      
      if(iext4.eq.0)then 
      
      if(izmode.eq.0)then
c        im=mod(irtfile-1,20)+1
        if(icc20.ne.1)stop'ERROR 19042014 '
c        bimevt=(zclass(1,im)+zclass(2,im))/2
      elseif(izmode.eq.1)then
        bimevt=x
      elseif(izmode.eq.2)then
        nprt(1)=x
      elseif(izmode.eq.3)then
        nglevt=x
      elseif(izmode.eq.4)then
        segevt=x
      else
        stop'in setcentrVar: wrong izmode'
      endif
      
      else !read from ordered tree
      
        bimevt=0
        nprt(1)=0
        nglevt=0
        segevt=0
      
      endif
      
      end

c----------------------------------------------------------------------
      subroutine getJKNcentr
c----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incho"
      include "epos.incems"
      common/cen/ncentr
      common/jcen/jcentr,jmxcentr /kcen/kcentr,kmxcentr
      kmax=0
      do k=1,100
        if(zclass(3,k).gt.0.)kmax=kmax+1
      enddo
      ncenthy=kmax
      kk=0
      if(izmode.eq.0.or.izmode.eq.1)then
        do k=1,100
         if(zclass(3,k).gt.0.and.
     .    bimevt.ge.zclass(1,k).and.bimevt.le.zclass(2,k))then
          kk=k
         endif
        enddo
      elseif(izmode.eq.2)then
        zp=nprt(1)
        do k=1,100
         if(zclass(3,k).gt.0.and.
     .    zp.ge.zclass(1,k).and.zp.le.zclass(2,k))then
          kk=k
         endif
        enddo
      elseif(izmode.eq.3)then
        zp=nglevt
        do k=1,5
         if(zclass(3,k).gt.0.and.
     .    zp.ge.zclass(1,k).and.zp.le.zclass(2,k))then
          kk=k
         endif
        enddo
        write(ifmt,'(a,i4)')'centrality class',kk
      else
        stop'in getJKNcentr: wrong jzmode\n\n'
      endif
      ncentr=kk

      if(jzmode(1).gt.0)call getJcentr
      if(jzmode(2).gt.0)call getKcentr

      if(jcentr.gt.ncenthxx)stop'getJKNcentr\n\n'

      end

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      FUNCTION sbessi0(x)
c-----------------------------------------------------------------------
c Returns the modified Bessel function I0(x) for any real x.
c   Numerical Recipes in Fortran 77: The art of scientific computing 
c   (ISBN 0-521-43064-X) Ch 6, page 230-232 
c   T. Pierog, 23/02/2010
c-----------------------------------------------------------------------

      REAL sbessi0,x
      REAL ax

      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,
     *     q8,q9,y !Accumulate polynomials in double precision.
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9

      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     *     1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     *     0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,
     *     0.2635537d-1,-0.1647633d-1,0.392377d-2/

      if (abs(x).lt.3.75) then

        y=(x/3.75)**2

        sbessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))

      else

        ax=abs(x)

        y=3.75/ax

        sbessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4
     *       +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))

      endif

      return

      END

c-----------------------------------------------------------------------
      FUNCTION sbessk0(x)
c-----------------------------------------------------------------------
c Returns the modified Bessel function K0(x) for positive real x.
c Numerical Recipes in Fortran 77: The art of scientific computing 
c (ISBN 0-521-43064-X) Ch 6, page 230-232 
c   T. Pierog, 23/02/2010
c-----------------------------------------------------------------------

      REAL sbessk0,x

C USES bessi0
      REAL sbessi0

      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,
     *     q2,q3,q4,q5,q6,q7,y !Accumulate polynomials in double precision.
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     *     0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     *     -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/

      if (x.le.2.0) then !Polynomial fit.

      y=x*x/4.0

      sbessk0=(-log(x/2.0)*sbessi0(x))+(p1+y*(p2+y*(p3+

     *     y*(p4+y*(p5+y*(p6+y*p7))))))

      else

        y=(2.0/x)

        sbessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+

     *       y*(q4+y*(q5+y*(q6+y*q7))))))

      endif

      return

      END

c-----------------------------------------------------------------------
      FUNCTION sbessi1(x)
c-----------------------------------------------------------------------
c Returns the modified Bessel function I1(x) for any real x.
c Numerical Recipes in Fortran 77: The art of scientific computing 
c (ISBN 0-521-43064-X) Ch 6, page 230-232 
c   T. Pierog, 23/02/2010
c-----------------------------------------------------------------------

      REAL sbessi1,x
      REAL ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,
     *     q8,q9,y          !Accumulate polynomials in double precision.

      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     *     0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     *     -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,
     *     -0.2895312d-1,0.1787654d-1,-0.420059d-2/

      if (abs(x).lt.3.75) then    !Polynomial fit.

        y=(x/3.75)**2

        sbessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))

      else

        ax=abs(x)

        y=3.75/ax

        sbessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+
     *       y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))

        if(x.lt.0.)sbessi1=-sbessi1

      endif

      return

      END

c-----------------------------------------------------------------------
      FUNCTION sbessk1(x)
c-----------------------------------------------------------------------
c Returns the modified Bessel function K1(x) for positive real x.
c Numerical Recipes in Fortran 77: The art of scientific computing 
c (ISBN 0-521-43064-X) Ch 6, page 230-232 
c   T. Pierog, 23/02/2010
c-----------------------------------------------------------------------

      REAL sbessk1,x

C     USES sbessi1
      REAL sbessi1

      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,
     *     q2,q3,q4,q5,q6,q7,y     !Accumulate polynomials in double precision.
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     *     -0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     *     0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/

      if (x.le.2.0) then    !Polynomial fit.

        y=x*x/4.0

        sbessk1=(log(x/2.0)*sbessi1(x))+(1.0/x)*(p1+y*(p2+
     *     y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))

      else

        y=2.0/x

        sbessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+
     *       y*(q4+y*(q5+y*(q6+y*q7))))))

      endif

      return

      END


      
C  (C) Copr. 1986-92 Numerical Recipes Software 4+1$!].


ckwc-----------------------------------------------------------------------
ckw      subroutine testspeed
ckwc-----------------------------------------------------------------------
ckw      double precision a
ckw      common/aaa/nlidtbl(-9990:9900),amtbl(580)
ckw      integer itab(2,580)
ckw      data itab/
ckw     .          120 ,      139570
ckw     .,        -120 ,      139570
ckw     .,         110 ,      134976
ckw     .,         220 ,      547853
ckw     .,         703 ,      800000
ckw     .,         121 ,      775490
ckw     .,        -121 ,      775490
ckw     .,         111 ,      775490
ckw     .,         221 ,      782650
ckw     .,         330 ,      957660
ckw     .,         702 ,      980000
ckw     .,         773 ,      984700
ckw     .,        -773 ,      984700
ckw     .,         782 ,      984700
ckw     .,         331 ,     1019454
ckw     .,         112 ,     1170000
ckw     .,         122 ,     1229500
ckw     .,        -122 ,     1229500
ckw     .,         222 ,     1229500
ckw     .,         123 ,     1230000
ckw     .,        -123 ,     1230000
ckw     .,         223 ,     1230000
ckw     .,         113 ,     1275100
ckw     .,         114 ,     1281800
ckw     .,         729 ,     1294000
ckw     .,         124 ,     1300000
ckw     .,        -124 ,     1300000
ckw     .,         115 ,     1300000
ckw     .,         125 ,     1318300
ckw     .,        -125 ,     1318300
ckw     .,         224 ,     1318300
ckw     .,         332 ,     1350000
ckw     .,         774 ,     1351000
ckw     .,        -774 ,     1351000
ckw     .,         792 ,     1351000
ckw     .,         732 ,     1409800
ckw     .,         334 ,     1426400
ckw     .,         733 ,     1425000
ckw     .,         127 ,     1474000
ckw     .,        -127 ,     1474000
ckw     .,         225 ,     1474000
ckw     .,         126 ,     1465000
ckw     .,        -126 ,     1465000
ckw     .,         117 ,     1465000
ckw     .,         796 ,     1476000
ckw     .,         226 ,     1505000
ckw     .,         337 ,     1525000
ckw     .,         777 ,     1662000
ckw     .,        -777 ,     1662000
ckw     .,         793 ,     1662000
ckw     .,         706 ,     1617000
ckw     .,         707 ,     1670000
ckw     .,         708 ,     1667000
ckw     .,         129 ,     1672400
ckw     .,        -129 ,     1672400
ckw     .,         119 ,     1672400
ckw     .,         797 ,     1680000
ckw     .,         212 ,     1688800
ckw     .,        -212 ,     1688800
ckw     .,         227 ,     1688800
ckw     .,         213 ,     1720000
ckw     .,        -213 ,     1720000
ckw     .,         228 ,     1720000
ckw     .,         339 ,     1724000
ckw     .,         215 ,     1816000
ckw     .,        -215 ,     1816000
ckw     .,         789 ,     1816000
ckw     .,         794 ,     1854000
ckw     .,         714 ,     1944000
ckw     .,         715 ,     2040000
ckw     .,         218 ,     2000999
ckw     .,        -218 ,     2000999
ckw     .,         785 ,     2000999
ckw     .,         717 ,     2017999
ckw     .,         725 ,     2297000
ckw     .,         728 ,     2340000
ckw     .,         130 ,      493677
ckw     .,        -130 ,      493677
ckw     .,         230 ,      497614
ckw     .,        -230 ,      497614
ckw     .,         131 ,      891660
ckw     .,        -131 ,      891660
ckw     .,         231 ,      896000
ckw     .,        -231 ,      896000
ckw     .,         133 ,     1272000
ckw     .,        -133 ,     1272000
ckw     .,         233 ,     1272000
ckw     .,        -233 ,     1272000
ckw     .,         134 ,     1403000
ckw     .,        -134 ,     1403000
ckw     .,         234 ,     1403000
ckw     .,        -234 ,     1403000
ckw     .,         135 ,     1414000
ckw     .,        -135 ,     1414000
ckw     .,         235 ,     1414000
ckw     .,        -235 ,     1414000
ckw     .,         136 ,     1420000
ckw     .,        -136 ,     1420000
ckw     .,         236 ,     1420000
ckw     .,        -236 ,     1420000
ckw     .,         137 ,     1425600
ckw     .,        -137 ,     1425600
ckw     .,         237 ,     1432400
ckw     .,        -237 ,     1432400
ckw     .,         313 ,     1717000
ckw     .,        -313 ,     1717000
ckw     .,         323 ,     1717000
ckw     .,        -323 ,     1717000
ckw     .,         314 ,     1773000
ckw     .,        -314 ,     1773000
ckw     .,         324 ,     1773000
ckw     .,        -324 ,     1773000
ckw     .,         315 ,     1776000
ckw     .,        -315 ,     1776000
ckw     .,         325 ,     1776000
ckw     .,        -325 ,     1776000
ckw     .,         316 ,     1816000
ckw     .,        -316 ,     1816000
ckw     .,         326 ,     1816000
ckw     .,        -326 ,     1816000
ckw     .,         776 ,     2045000
ckw     .,        -776 ,     2045000
ckw     .,         805 ,     2045000
ckw     .,        -805 ,     2045000
ckw     .,         440 ,     2980300
ckw     .,         441 ,     3096916
ckw     .,        -442 ,     3414750
ckw     .,        -443 ,     3510660
ckw     .,        -445 ,     3556200
ckw     .,         446 ,     3637000
ckw     .,         447 ,     3686090
ckw     .,         448 ,     3772920
ckw     .,         752 ,     4039000
ckw     .,         753 ,     4153000
ckw     .,         754 ,     4421000
ckw     .,         551 ,     9460300
ckw     .,         552 ,     9859400
ckw     .,         553 ,     9892800
ckw     .,         555 ,     9912200
ckw     .,        1120 ,      938272
ckw     .,       -1120 ,      938272
ckw     .,        1220 ,      939565
ckw     .,       -1220 ,      939565
ckw     .,        1122 ,     1440000
ckw     .,       -1122 ,     1440000
ckw     .,        1222 ,     1440000
ckw     .,       -1222 ,     1440000
ckw     .,        1123 ,     1520000
ckw     .,       -1123 ,     1520000
ckw     .,        1223 ,     1520000
ckw     .,       -1223 ,     1520000
ckw     .,        1124 ,     1535000
ckw     .,       -1124 ,     1535000
ckw     .,        1224 ,     1535000
ckw     .,       -1224 ,     1535000
ckw     .,        1127 ,     1655000
ckw     .,       -1127 ,     1655000
ckw     .,        1227 ,     1655000
ckw     .,       -1227 ,     1655000
ckw     .,        1128 ,     1675000
ckw     .,       -1128 ,     1675000
ckw     .,        1228 ,     1675000
ckw     .,       -1228 ,     1675000
ckw     .,        1129 ,     1685000
ckw     .,       -1129 ,     1685000
ckw     .,        1229 ,     1685000
ckw     .,       -1229 ,     1685000
ckw     .,        2113 ,     1700000
ckw     .,       -2113 ,     1700000
ckw     .,        2213 ,     1700000
ckw     .,       -2213 ,     1700000
ckw     .,        2114 ,     1710000
ckw     .,       -2114 ,     1710000
ckw     .,        2214 ,     1710000
ckw     .,       -2214 ,     1710000
ckw     .,        2115 ,     1720000
ckw     .,       -2115 ,     1720000
ckw     .,        2215 ,     1720000
ckw     .,       -2215 ,     1720000
ckw     .,        1214 ,     2190000
ckw     .,       -1214 ,     2190000
ckw     .,        2124 ,     2190000
ckw     .,       -2124 ,     2190000
ckw     .,        1111 ,     1232000
ckw     .,       -1111 ,     1232000
ckw     .,        1121 ,     1232000
ckw     .,       -1121 ,     1232000
ckw     .,        1221 ,     1232000
ckw     .,       -1221 ,     1232000
ckw     .,        2221 ,     1232000
ckw     .,       -2221 ,     1232000
ckw     .,        1112 ,     1600000
ckw     .,       -1112 ,     1600000
ckw     .,        1125 ,     1600000
ckw     .,       -1125 ,     1600000
ckw     .,        1225 ,     1600000
ckw     .,       -1225 ,     1600000
ckw     .,        2222 ,     1600000
ckw     .,       -2222 ,     1600000
ckw     .,        1113 ,     1630000
ckw     .,       -1113 ,     1630000
ckw     .,        1126 ,     1630000
ckw     .,       -1126 ,     1630000
ckw     .,        1226 ,     1630000
ckw     .,       -1226 ,     1630000
ckw     .,        2223 ,     1630000
ckw     .,       -2223 ,     1630000
ckw     .,        1114 ,     1700000
ckw     .,       -1114 ,     1700000
ckw     .,        2112 ,     1700000
ckw     .,       -2112 ,     1700000
ckw     .,        2212 ,     1700000
ckw     .,       -2212 ,     1700000
ckw     .,        2224 ,     1700000
ckw     .,       -2224 ,     1700000
ckw     .,        1116 ,     1890000
ckw     .,       -1116 ,     1890000
ckw     .,        2117 ,     1890000
ckw     .,       -2117 ,     1890000
ckw     .,        2217 ,     1890000
ckw     .,       -2217 ,     1890000
ckw     .,        2226 ,     1890000
ckw     .,       -2226 ,     1890000
ckw     .,        1117 ,     1910000
ckw     .,       -1117 ,     1910000
ckw     .,        2118 ,     1910000
ckw     .,       -2118 ,     1910000
ckw     .,        2218 ,     1910000
ckw     .,       -2218 ,     1910000
ckw     .,        2227 ,     1910000
ckw     .,       -2227 ,     1910000
ckw     .,        1118 ,     1920000
ckw     .,       -1118 ,     1920000
ckw     .,        2119 ,     1920000
ckw     .,       -2119 ,     1920000
ckw     .,        2219 ,     1920000
ckw     .,       -2219 ,     1920000
ckw     .,        2228 ,     1920000
ckw     .,       -2228 ,     1920000
ckw     .,        1119 ,     1960000
ckw     .,       -1119 ,     1960000
ckw     .,        1212 ,     1960000
ckw     .,       -1212 ,     1960000
ckw     .,        2122 ,     1960000
ckw     .,       -2122 ,     1960000
ckw     .,        2229 ,     1960000
ckw     .,       -2229 ,     1960000
ckw     .,        7002 ,     1930000
ckw     .,       -7002 ,     1930000
ckw     .,        1213 ,     1930000
ckw     .,       -1213 ,     1930000
ckw     .,        2123 ,     1930000
ckw     .,       -2123 ,     1930000
ckw     .,        7003 ,     1930000
ckw     .,       -7003 ,     1930000
ckw     .,        2130 ,     1115683
ckw     .,       -2130 ,     1115683
ckw     .,        1233 ,     1406000
ckw     .,       -1233 ,     1406000
ckw     .,        1234 ,     1519500
ckw     .,       -1234 ,     1519500
ckw     .,        1235 ,     1600000
ckw     .,       -1235 ,     1600000
ckw     .,        1236 ,     1670000
ckw     .,       -1236 ,     1670000
ckw     .,        1237 ,     1690000
ckw     .,       -1237 ,     1690000
ckw     .,        3124 ,     1800000
ckw     .,       -3124 ,     1800000
ckw     .,        3125 ,     1810000
ckw     .,       -3125 ,     1810000
ckw     .,        3126 ,     1820000
ckw     .,       -3126 ,     1820000
ckw     .,        3127 ,     1830000
ckw     .,       -3127 ,     1830000
ckw     .,        3128 ,     1890000
ckw     .,       -3128 ,     1890000
ckw     .,        3217 ,     2100000
ckw     .,       -3217 ,     2100000
ckw     .,        3214 ,     2110000
ckw     .,       -3214 ,     2110000
ckw     .,        1130 ,     1189370
ckw     .,       -1130 ,     1189370
ckw     .,        1230 ,     1192642
ckw     .,       -1230 ,     1192642
ckw     .,        2230 ,     1197449
ckw     .,       -2230 ,     1197449
ckw     .,        1131 ,     1382800
ckw     .,       -1131 ,     1382800
ckw     .,        1231 ,     1383700
ckw     .,       -1231 ,     1383700
ckw     .,        2231 ,     1387200
ckw     .,       -2231 ,     1387200
ckw     .,        1132 ,     1660000
ckw     .,       -1132 ,     1660000
ckw     .,        1238 ,     1660000
ckw     .,       -1238 ,     1660000
ckw     .,        2232 ,     1660000
ckw     .,       -2232 ,     1660000
ckw     .,        1133 ,     1670000
ckw     .,       -1133 ,     1670000
ckw     .,        1239 ,     1670000
ckw     .,       -1239 ,     1670000
ckw     .,        2233 ,     1670000
ckw     .,       -2233 ,     1670000
ckw     .,        1134 ,     1750000
ckw     .,       -1134 ,     1750000
ckw     .,        3122 ,     1750000
ckw     .,       -3122 ,     1750000
ckw     .,        2234 ,     1750000
ckw     .,       -2234 ,     1750000
ckw     .,        1135 ,     1775000
ckw     .,       -1135 ,     1775000
ckw     .,        3123 ,     1775000
ckw     .,       -3123 ,     1775000
ckw     .,        2235 ,     1775000
ckw     .,       -2235 ,     1775000
ckw     .,        1136 ,     1915000
ckw     .,       -1136 ,     1915000
ckw     .,        3129 ,     1915000
ckw     .,       -3129 ,     1915000
ckw     .,        2236 ,     1915000
ckw     .,       -2236 ,     1915000
ckw     .,        1137 ,     1940000
ckw     .,       -1137 ,     1940000
ckw     .,        3212 ,     1940000
ckw     .,       -3212 ,     1940000
ckw     .,        2237 ,     1940000
ckw     .,       -2237 ,     1940000
ckw     .,        1138 ,     2030000
ckw     .,       -1138 ,     2030000
ckw     .,        3213 ,     2030000
ckw     .,       -3213 ,     2030000
ckw     .,        2238 ,     2030000
ckw     .,       -2238 ,     2030000
ckw     .,        1330 ,     1314860
ckw     .,       -1330 ,     1314860
ckw     .,        2330 ,     1321710
ckw     .,       -2330 ,     1321710
ckw     .,        1331 ,     1531800
ckw     .,       -1331 ,     1531800
ckw     .,        2331 ,     1535000
ckw     .,       -2331 ,     1535000
ckw     .,        1334 ,     1823000
ckw     .,       -1334 ,     1823000
ckw     .,        2334 ,     1823000
ckw     .,       -2334 ,     1823000
ckw     .,        3331 ,     1672450
ckw     .,       -3331 ,     1672450
ckw     .,        1215 ,     2250000
ckw     .,       -1215 ,     2250000
ckw     .,        2125 ,     2250000
ckw     .,       -2125 ,     2250000
ckw     .,        1216 ,     2280000
ckw     .,       -1216 ,     2280000
ckw     .,        2126 ,     2280000
ckw     .,       -2126 ,     2280000
ckw     .,        1218 ,     2600000
ckw     .,       -1218 ,     2600000
ckw     .,        2128 ,     2600000
ckw     .,       -2128 ,     2600000
ckw     .,        7004 ,     2420000
ckw     .,       -7004 ,     2420000
ckw     .,        1217 ,     2420000
ckw     .,       -1217 ,     2420000
ckw     .,        2127 ,     2420000
ckw     .,       -2127 ,     2420000
ckw     .,        7005 ,     2420000
ckw     .,       -7005 ,     2420000
ckw     .,        3216 ,     2350000
ckw     .,       -3216 ,     2350000
ckw     .,         556 ,    10023260
ckw     .,         762 ,    10232500
ckw     .,         763 ,    10255500
ckw     .,         764 ,    10268600
ckw     .,         557 ,    10355200
ckw     .,         558 ,    10579400
ckw     .,         559 ,    10865000
ckw     .,         765 ,    11019000
ckw     .,          10 ,           0
ckw     .,          12 ,         510
ckw     .,         -12 ,         510
ckw     .,          11 ,           0
ckw     .,         -11 ,           0
ckw     .,          14 ,      105658
ckw     .,         -14 ,      105658
ckw     .,          13 ,           0
ckw     .,         -13 ,           0
ckw     .,          16 ,     1776840
ckw     .,         -16 ,     1776840
ckw     .,          15 ,           0
ckw     .,         -15 ,           0
ckw     .,         335 ,     1453000
ckw     .,         336 ,     1518000
ckw     .,         338 ,     1562000
ckw     .,         704 ,     1594000
ckw     .,         128 ,     1647000
ckw     .,        -128 ,     1647000
ckw     .,         118 ,     1647000
ckw     .,         734 ,     1639000
ckw     .,         214 ,     1732000
ckw     .,        -214 ,     1732000
ckw     .,         229 ,     1732000
ckw     .,         709 ,     1756000
ckw     .,         712 ,     1815000
ckw     .,         795 ,     1842000
ckw     .,         216 ,     1909000
ckw     .,        -216 ,     1909000
ckw     .,         783 ,     1909000
ckw     .,         713 ,     1903000
ckw     .,         217 ,     1982000
ckw     .,        -217 ,     1982000
ckw     .,         784 ,     1982000
ckw     .,         716 ,     1992000
ckw     .,         219 ,     2089999
ckw     .,        -219 ,     2089999
ckw     .,         786 ,     2089999
ckw     .,         718 ,     2103000
ckw     .,         719 ,     2156000
ckw     .,         775 ,     2149000
ckw     .,        -775 ,     2149000
ckw     .,         787 ,     2149000
ckw     .,         722 ,     2189000
ckw     .,         723 ,     2231100
ckw     .,         724 ,     2220000
ckw     .,         772 ,     2260000
ckw     .,        -772 ,     2260000
ckw     .,         788 ,     2260000
ckw     .,         726 ,     2320000
ckw     .,         722 ,     2314000
ckw     .,         132 ,      670000
ckw     .,        -132 ,      670000
ckw     .,         232 ,      670000
ckw     .,        -232 ,      670000
ckw     .,         138 ,     1460000
ckw     .,        -138 ,     1460000
ckw     .,         238 ,     1460000
ckw     .,        -238 ,     1460000
ckw     .,         139 ,     1580000
ckw     .,        -139 ,     1580000
ckw     .,         239 ,     1580000
ckw     .,        -239 ,     1580000
ckw     .,         312 ,     1650000
ckw     .,        -312 ,     1650000
ckw     .,         322 ,     1650000
ckw     .,        -322 ,     1650000
ckw     .,         317 ,     1830000
ckw     .,        -317 ,     1830000
ckw     .,         327 ,     1830000
ckw     .,        -327 ,     1830000
ckw     .,         318 ,     1945000
ckw     .,        -318 ,     1945000
ckw     .,         328 ,     1945000
ckw     .,        -328 ,     1945000
ckw     .,         319 ,     1973000
ckw     .,        -319 ,     1973000
ckw     .,         329 ,     1973000
ckw     .,        -329 ,     1973000
ckw     .,         742 ,     2247000
ckw     .,        -742 ,     2247000
ckw     .,         802 ,     2247000
ckw     .,        -802 ,     2247000
ckw     .,         743 ,     2324000
ckw     .,        -743 ,     2324000
ckw     .,         803 ,     2324000
ckw     .,        -803 ,     2324000
ckw     .,         744 ,     2490000
ckw     .,        -744 ,     2490000
ckw     .,         804 ,     2490000
ckw     .,        -804 ,     2490000
ckw     .,        -240 ,     1869620
ckw     .,         240 ,     1869620
ckw     .,        -140 ,     1864840
ckw     .,         140 ,     1864840
ckw     .,        -141 ,     2006969
ckw     .,         141 ,     2006969
ckw     .,        -241 ,     2010270
ckw     .,         241 ,     2010270
ckw     .,        -142 ,     2350000
ckw     .,         142 ,     2350000
ckw     .,        -242 ,     2400000
ckw     .,         242 ,     2400000
ckw     .,        -143 ,     2422300
ckw     .,         143 ,     2422300
ckw     .,        -144 ,     2430000
ckw     .,         144 ,     2430000
ckw     .,        -145 ,     2461100
ckw     .,         145 ,     2461100
ckw     .,        -245 ,     2460100
ckw     .,         245 ,     2460100
ckw     .,        -340 ,     1968490
ckw     .,         340 ,     1968490
ckw     .,        -341 ,     2112300
ckw     .,         341 ,     2112300
ckw     .,        -342 ,     2317800
ckw     .,         342 ,     2317800
ckw     .,        -343 ,     2459600
ckw     .,         343 ,     2459600
ckw     .,        -344 ,     2535400
ckw     .,         344 ,     2535400
ckw     .,         150 ,     5279150
ckw     .,        -150 ,     5279150
ckw     .,         250 ,     5279530
ckw     .,        -250 ,     5279530
ckw     .,         151 ,     5325100
ckw     .,        -151 ,     5325100
ckw     .,         251 ,     5325100
ckw     .,        -251 ,     5325100
ckw     .,         155 ,     5746900
ckw     .,        -155 ,     5746900
ckw     .,         255 ,     5746900
ckw     .,        -255 ,     5746900
ckw     .,         350 ,     5366300
ckw     .,        -350 ,     5366300
ckw     .,         351 ,     5412800
ckw     .,        -351 ,     5412800
ckw     .,         353 ,     5839700
ckw     .,        -353 ,     5839700
ckw     .,         450 ,     6276000
ckw     .,        -450 ,     6276000
ckw     .,         449 ,     3929000
ckw     .,         550 ,     9300000
ckw     .,           1 ,     1900000
ckw     .,           2 ,     1900000
ckw     .,           3 ,     1900000
ckw     .,           4 ,     1900000
ckw     .,           5 ,     1900000
ckw     .,           6 ,     1900000
ckw     .,           7 ,     1900000
ckw     .,           8 ,     1900000
ckw     .,        2140 ,     2286460
ckw     .,       -2140 ,     2286460
ckw     .,        1242 ,     2595400
ckw     .,       -1242 ,     2595400
ckw     .,        1140 ,     2454020
ckw     .,       -1140 ,     2454020
ckw     .,        1240 ,     2452900
ckw     .,       -1240 ,     2452900
ckw     .,        2240 ,     2453760
ckw     .,       -2240 ,     2453760
ckw     .,        1141 ,     2518400
ckw     .,       -1141 ,     2518400
ckw     .,        1241 ,     2517500
ckw     .,       -1241 ,     2517500
ckw     .,        2241 ,     2518000
ckw     .,       -2241 ,     2518000
ckw     .,        3140 ,     2467900
ckw     .,       -3140 ,     2467900
ckw     .,        3240 ,     2471000
ckw     .,       -3240 ,     2471000
ckw     .,        1340 ,     2575700
ckw     .,       -1340 ,     2575700
ckw     .,        2340 ,     2578000
ckw     .,       -2340 ,     2578000
ckw     .,        1341 ,     2646600
ckw     .,       -1341 ,     2646600
ckw     .,        2341 ,     2646100
ckw     .,       -2341 ,     2646100
ckw     .,        3340 ,     2697500
ckw     .,       -3340 ,     2697500
ckw     .,        3341 ,     2768300
ckw     .,       -3341 ,     2768300
ckw     .,        2150 ,     5620200
ckw     .,       -2150 ,     5620200
ckw     .,        1150 ,     5807800
ckw     .,           9 ,     5807800
ckw     .,        1250 ,     5807800
ckw     .,       -1250 ,     5807800
ckw     .,        2250 ,     5815200
ckw     .,          -1 ,     5815200
ckw     .,        1151 ,     5829000
ckw     .,          -2 ,     5829000
ckw     .,        1251 ,     5829000
ckw     .,       -1251 ,     5829000
ckw     .,        2251 ,     5836400
ckw     .,          -3 ,     5836400
ckw     .,        3150 ,     5792400
ckw     .,       -3150 ,     5792400
ckw     .,        3250 ,     5792400
ckw     .,       -3250 ,     5792400/
ckw      !print*,itab(1,1),itab(2,1)
ckw      !print*,itab(1,2),itab(2,2)
ckw      !print*,itab(1,3),itab(2,3)
ckw      do nl=1,580
ckw         id=itab(1,nl)
ckw         nlidtbl(id)=nl
ckw         amtbl(nl)=float( itab(2,nl) ) / 1000000
ckw      enddo
ckw      print*,120, nlidtbl(120),amtbl(nlidtbl(120))
ckw      print*,1120, nlidtbl(1120),amtbl(nlidtbl(1120))
ckw      print*,3331, nlidtbl(3331),amtbl(nlidtbl(3331))
ckw        
ckw      a=0
ckw      print*, 'do the sum ...'
ckw      do i=1,100000000
ckw        call idmssxxx(120,a1)
ckw        call idmssxxx(130,a2)
ckw        call idmssxxx(441,a3)
ckw        call idmssxxx(3331,a4)
ckw        call idmssxxx(140,a5)
ckw        call idmssxxx(1251,a6)
ckw        a=a+a1+a2+a3+a4+a5+a6
ckw      enddo
ckw      print*,'result of sum',a
ckw      end
ckw
ckw      subroutine idmssxxx(id,amass)
ckw      common/aaa/nlidtbl(-9990:9900),amtbl(580)
ckw        nl=nlidtbl(id)
ckw        amass=amtbl(nl) 
ckw      end


