C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
      integer function ideposf(iii,j)
c-----------------------------------------------------------------------
      include "epos.inc"
      if(ityptl(j).eq.61)then
        ideposf=idtrafo('pdg','nxs',idptl(j))
      else
        ideposf=idptl(j)
      endif
      iiii=iii !to avoid warning
ctp VERY time consuming (20% of the TOTAL time) !!!!
c      call idchrg( 998 ,ideposf,ch) !for testing  
c      if(ch.gt.1e29)then
c        print*,'ERROR 11072016a',ideposf,' not in list.'
c     .,ityptl(j),idptl(j),iii
        
c        stop
c      endif
      end

c-----------------------------------------------------------------------
      subroutine istu12(istuse,isu1,isu2)
c-----------------------------------------------------------------------
      include "epos.inc"
      iii=istuse
      if(istfor.gt.-999)
     . iii  = sign(1,istfor)*(abs(istfor)+10*abs(istfor))
      if(iii.gt.-999)then
        isiuse=1
        if(iii.lt.0)isiuse=-1
        isu1x=isiuse*mod(abs(iii),10)
        isu2x=isiuse*abs(iii)/10
        if(istfor.gt.-999)then
          if(isu1x.ne.istfor.or.isu2x.ne.istfor)
     .    stop'##### ERROR 05122015 #####'          
        endif
        isu1=min(isu1x,isu2x)
        isu2=max(isu1x,isu2x)
        if(isu1.eq.6.and.isu2.eq.8)then
          isu1=8  !last generation
          isu2=6
        endif
      else
        isu1=-999
        isu2=-999
      endif
      end
      integer function istu1(istuse)
      call istu12(istuse,isu1,isu2)
      istu1=isu1
      end  
      integer function istu2(istuse)
      call istu12(istuse,isu1,isu2)
      istu2=isu2
      end  

c-----------------------------------------------------------------------
      logical function LongLivPtl(istsoll,i)
c-----------------------------------------------------------------------
c 'Long Lived weakly decaying Ptl'  =   1cm < ctau < 10cm   (Ks,Lda,Sig,Xi,Omg)
c-----------------------------------------------------------------------
      include 'epos.inc'
      LongLivPtl=.false.
      ist=istptl(i)
      idepos=idptl(i)
      if(ist.eq.istsoll)then
        idaepos=abs(idepos)
        if( idaepos.eq.20   .or.idaepos.eq.2130 
     &  .or.idaepos.eq.2230 .or.idaepos.eq.1130 
     &  .or.idaepos.eq.2330 .or.idaepos.eq.1330 
     &  .or.idaepos.eq.3331                    )LongLivPtl=.true.
      endif
      end
c-----------------------------------------------------------------------
      logical function NoLongLivParent(noweak,j)
c-----------------------------------------------------------------------
c 'Long Lived weakly decaying Ptl'  =   1cm < ctau < 10cm   (Ks,Lda,Sig,Xi,Omg)
c-----------------------------------------------------------------------
      include 'epos.inc'
      NoLongLivParent=.true.
      idj=idptl(j)
      idjj=idj
      iorj=iorptl(j) !parent
  77  continue
      if(iorj.gt.0)then
        idorj=idptl(iorj)
        if(idorj.ne.idjj)then !skip ptl copies 0,1 to 8,6
          idora=abs( idorj )
          if(noweak.eq.1.or.noweak.eq.3)then !exclude children from long lived weakly decaying ptls 
            if( idora.eq.20   .or.idora.eq.2130 
     &      .or.idora.eq.2230 .or.idora.eq.1130 
     &      .or.idora.eq.2330 .or.idora.eq.1330 
     &      .or.idora.eq.3331                )NoLongLivParent=.false.
            ! print *, j,n, '   ', idj,idjj,idora,go
          elseif(noweak.eq.2)then
            if( idora.eq.150 .or.idora.eq.250.or.
     &          idora.eq.350 .or.idora.eq.450)NoLongLivParent=.false.
          endif
          istiorj=istptl(iorj)
          if(istiorj.eq.1.or.istiorj.eq.6)then   !parent=hadron
            iorj2=iorptl(iorj)  !grandparent  
            iorj=iorj2
            idjj=idorj!keep parent id
            goto 77
          endif
        endif
      endif
      end

c-----------------------------------------------------------------------
      subroutine setIstXY(n56,istuseN,istxxx,istyyy)
c-----------------------------------------------------------------------
      include 'epos.inc'
      istxxx=0
      istyyy=1
      if(n56.eq.99)return
      if(istu1(istuseN).gt.-999)istxxx=istu1(istuseN)
      if(istu2(istuseN).gt.-999)istyyy=istu2(istuseN)
      if(.not.(   (istxxx.eq.0.and.istyyy.eq.1)  !consider two  \ do noweak
     .        .or.(istxxx.eq.8.and.istyyy.eq.6)  !  generations /  analysis
     .        .or.(istxxx.eq.0.and.istyyy.eq.0)  !only        \  do single
     .        .or.(istxxx.eq.8.and.istyyy.eq.8)  !  single     > generation
     .        .or.(istxxx.eq.6.and.istyyy.eq.6)  ! generation /  analysis
     .        .or.(istxxx.eq.0.and.istyyy.eq.5) 
     .        .or.(istxxx.eq.0.and.istyyy.eq.7) 
     .   ))then
        print*,'ERROR04122021', n56,istuseN,istxxx,istyyy
        stop'ERROR04122021'
      endif  
      if(ihacas.eq.0)then
        if(istxxx.eq.8.and.istyyy.eq.8)then
          istxxx=0
          istyyy=0
        elseif(istxxx.eq.6.and.istyyy.eq.6)then
          istxxx=1
          istyyy=1
        elseif(istxxx.eq.8.and.istyyy.eq.6)then
          istxxx=0
          istyyy=1
        endif
      endif
      end

c-----------------------------------------------------------------------
      subroutine iclass(id,icl)
c-----------------------------------------------------------------------
c      determines hadron class
c-----------------------------------------------------------------------
      ida=iabs(id)
      if(ida.eq.0.or.(ida.ge.17.and.ida.le.19))then
       icl=2
      elseif(ida/10.eq.13.or.ida/10.eq.23.or.ida/10.eq.33
     .   .or.ida.eq.20)then
       icl=3
      elseif(ida/10.eq.14.or.ida/10.eq.24.or.ida/10.eq.34
     .   .or.ida/10.eq.44)then
       icl=4
      elseif(ida.ge.100.and.ida.le.999)then
       icl=1
      elseif(ida.ge.1000.and.ida.le.9900)then
       icl=2
      else
       stop'iclass: id not known'
      endif
      end

c-----------------------------------------------------------------------
      subroutine idspin(id,iso,jspin,istra)
c     computes iso (isospin), jspin and istra (strangeness) of particle
c     with ident code id
c-----------------------------------------------------------------------
      include 'epos.inc'
      dimension ifl(3)
      iso=0
      jspin=0
      istra=0
      idabs=abs(id)
      if(idabs.le.nflav)then
        iso=isospin(idabs)*sign(1,id)
        if(idabs.ge.3)istra=sign(1,id)
        return
      endif
      call idflav(id,ifl(1),ifl(2),ifl(3),jspin,ind)
      iq1=abs(ifl(1))
      iq2=abs(ifl(2))
      iq3=abs(ifl(3))
      if(iq1.ge.3)istra=istra+sign(1,ifl(1))
      if(iq2.ge.3)istra=istra+sign(1,ifl(2))
      if(iq3.ge.3)istra=istra+sign(1,ifl(3))
      if(iq1.ne.0)then         !baryon
        iso=(isospin(iq1)+isospin(iq2)+isospin(iq3))*sign(1,id)
      else
        iso=isospin(iq2)*sign(1,ifl(2))
        iso=iso+isospin(iq3)*sign(1,ifl(3))
      endif
      return
      end

c-----------------------------------------------------------------------
      subroutine idcomk(ic)
c     compactifies ic
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer ic(2),icx(2),jc(nflav,2)
      call idcomp(ic,icx,jc,1)
      ic(1)=icx(1)
      ic(2)=icx(2)
      return
      end

cc-----------------------------------------------------------------------
c      subroutine idcomi(ic,icx)
cc     compactifies ic
cc-----------------------------------------------------------------------
c      parameter (nflav=6)
c      integer ic(2),icx(2),jc(nflav,2)
c      call idcomp(ic,icx,jc,1)
c      return
c      end
c
c-----------------------------------------------------------------------
      subroutine idcomj(jc)
c     compactifies jc
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer ic(2),icx(2),jc(nflav,2)
      call idcomp(ic,icx,jc,2)
      return
      end

c-----------------------------------------------------------------------
      subroutine idcomp(ic,icx,jc,im)
c-----------------------------------------------------------------------
c     compactifies ic,jc
c     input: im (1 or 2)
c            ic (if im=1)
c            jc (if im=2)
c     output: icx (if im=1)
c             jc
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer ic(2),icx(2),jc(nflav,2)
      if(im.eq.1)call iddeco(ic,jc)
      icx(1)=0
      icx(2)=0
           do n=1,nflav
           do j=1,2
      if(jc(n,j).ne.0)goto1
           enddo
           enddo
      return
1     continue
      nq=0
      na=0
           do n=1,nflav
      nq=nq+jc(n,1)
      na=na+jc(n,2)
           enddo
      l=0
           do n=1,nflav
      k=min0(jc(n,1),jc(n,2))
      if(nq.eq.1.and.na.eq.1)k=0
      jc(n,1)=jc(n,1)-k
      jc(n,2)=jc(n,2)-k
      if(jc(n,1).lt.0.or.jc(n,2).lt.0)
     *call utstop('idcomp: jc negative&')
      l=l+jc(n,1)+jc(n,2)
           enddo
           if(l.eq.0)then
      jc(1,1)=1
      jc(1,2)=1
           endif
           if(im.eq.1)then
      call idenco(jc,icx,ireten)
      if(ireten.eq.1)call utstop('idcomp: idenco ret code = 1&')
           endif
      return
      end

c-----------------------------------------------------------------------
      subroutine iddeco(ic,jc)
c     decode particle id
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2),ic(2)
      ici=ic(1)
      jc(6,1)=mod(ici,10)
      jc(5,1)=mod(ici/10,10)
      jc(4,1)=mod(ici/100,10)
      jc(3,1)=mod(ici/1000,10)
      jc(2,1)=mod(ici/10000,10)
      jc(1,1)=mod(ici/100000,10)
      ici=ic(2)
      jc(6,2)=mod(ici,10)
      jc(5,2)=mod(ici/10,10)
      jc(4,2)=mod(ici/100,10)
      jc(3,2)=mod(ici/1000,10)
      jc(2,2)=mod(ici/10000,10)
      jc(1,2)=mod(ici/100000,10)
      return
      end

c-----------------------------------------------------------------------
      subroutine idenco(jc,ic,ireten)
c     encode particle id
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2),ic(2)
      ireten=0
      ic(1)=0
      do 20 i=1,nflav
      if(jc(i,1).ge.10)goto22
20    ic(1)=ic(1)+jc(i,1)*10**(nflav-i)
      ic(2)=0
      do 21 i=1,nflav
      if(jc(i,2).ge.10)goto22
21    ic(2)=ic(2)+jc(i,2)*10**(nflav-i)
      return
22    ireten=1
      ic(1)=0
      ic(2)=0
      return
      end

c-----------------------------------------------------------------------
      subroutine idenct(jc,id,ib1,ib2,ib3,ib4)
c     encode particle id
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2),ic(2)

      id=0
      do 40 nf=1,nflav
      do 40 ij=1,2
        if(jc(nf,ij).ge.10)id=7*10**8
40    continue
           if(id/10**8.ne.7)then
      call idenco(jc,ic,ireten)
      if(ireten.eq.1)call utstop('idenct: idenco ret code = 1&')
      if(mod(ic(1),100).ne.0.or.mod(ic(2),100).ne.0)then
      id=9*10**8
      else
      id=8*10**8+ic(1)*100+ic(2)/100
      endif
           else
      call idtrbi(jc,ib1,ib2,ib3,ib4)
      id=id
     *+mod(jc(1,1)+jc(2,1)+jc(3,1)+jc(4,1),10**4)*10**4
     *+mod(jc(1,2)+jc(2,2)+jc(3,2)+jc(4,2),10**4)
           endif
      return
      end

c-----------------------------------------------------------------------
      subroutine idquacjc(jc,nqu,naq)
c     returns quark content of jc
c        jc(nflav,2) = jc-type particle identification code.
c        nqu = # quarks
c        naq = # antiquarks
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2)
      nqu=0
      naq=0
      do 53 n=1,nflav
        nqu=nqu+jc(n,1)
53      naq=naq+jc(n,2)
      return
      end

c-----------------------------------------------------------------------
      subroutine idqua(i,nq,ns,na)
c-----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incems"
      integer jc(nflav,2)
      call idquac(i,nq,ns,na,jc)
      end

c-----------------------------------------------------------------------
      subroutine idqua6(i,nq,ns,nc,nb,nt,na)
c-----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incems"
      integer jc(nflav,2)
      call idquac6(i,nq,ns,nc,nb,nt,na,jc)
      end

c-----------------------------------------------------------------------
      subroutine idquacic(ic,nqu,naq)
c     returns quark content of ic
c        ic(2) = ic-type particle identification code.
c        nqu = # quarks
c        naq = # antiquarks
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2),ic(2)
      nqu=0
      naq=0
      call iddeco(ic,jc)
      do 53 n=1,nflav
        nqu=nqu+jc(n,1)
53      naq=naq+jc(n,2)
      return
      end

c-----------------------------------------------------------------------
      subroutine idquac6(i,nq,ns,nc,nb,nt,na,jc)
c     returns quark content of ptl i from /cptl/ .
c        nq = # quarks - # antiquarks
c        ns = # strange quarks - # strange antiquarks
c        nc = # charm quarks   - # charm antiquarks
c        nb = # bottom quarks  - # bottom antiquarks
c        nt = # top quarks     - # top antiquarks
c        na = # quarks + # antiquarks
c        jc(nflav,2) = jc-type particle identification code.
c-----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incems"
      integer jc(nflav,2),ic(2)


      if(iabs(idptl(i)).eq.20)then
      idptl(i)=230
      if(rangen().lt..5)idptl(i)=-230
      goto 9999
      endif

      if(iabs(idptl(i)).lt.100)then
      nq=0
      ns=0
      nc=0
      nb=0
      nt=0
      do 1 n=1,nflav
      jc(n,1)=0
1     jc(n,2)=0
      return
      endif

9999  id=idptl(i)
      if(id/10**8.eq.7)then
        call idtrb(ibptl(1,i),ibptl(2,i),ibptl(3,i),ibptl(4,i),jc)
      elseif(id/10**8.eq.8)then
        ic(1)=mod(id,10**8)/10000*100
        ic(2)=mod(id,10**4)*100
        call iddeco(ic,jc)
      else
        call idtr6(idptl(i),ic)
        call iddeco(ic,jc)
      endif
      na=0
      nq=0
      do 53 n=1,nflav
      na=na+jc(n,1)+jc(n,2)
53    nq=nq+jc(n,1)-jc(n,2)
      ns=   jc(3,1)-jc(3,2)
      nc=   jc(4,1)-jc(4,2)
      nb=   jc(5,1)-jc(5,2)
      nt=   jc(6,1)-jc(6,2)
      nsum=0
      do n=1,nflav
        nsum=nsum+jc(n,1)
        nsum=nsum+jc(n,2)
      enddo
      if(nsum.eq.0)then
        write(ifmt,*)'----------------------------------'
        write(ifmt,*)'ERROR in idquac6, zero quarks'
        write(ifmt,*)'i =  ',i
        write(ifmt,*)'id =  ',idptl(i)
        write(ifmt,*)'ist = ',istptl(i)
        write(ifmt,*)'ity = ',ityptl(i)
        write(ifmt,*)'ib =  ',(ibptl(j,i),j=1,4)
        write(ifmt,*)'----------------------------------'
        stop
      endif 
      return
      end


c-----------------------------------------------------------------------
      subroutine idquac(i,nq,ns,na,jc)
c     returns quark content of ptl i from /cptl/ .
c        nq = # quarks - # antiquarks
c        ns = # strange quarks - # strange antiquarks
c        na = # quarks + # antiquarks
c        jc(nflav,2) = jc-type particle identification code.
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      parameter (nflavt=6)
      integer jc(nflavt,2),ic(2)


      if(iabs(idptl(i)).eq.20)then
      idptl(i)=230
      if(rangen().lt..5)idptl(i)=-230
      goto 9999
      endif

      if(iabs(idptl(i)).lt.100)then
      nq=0
      ns=0
      do 1 n=1,nflav
      jc(n,1)=0
1     jc(n,2)=0
      return
      endif
9999  if(idptl(i)/10**8.ne.7)then
      call idtr6(idptl(i),ic)
      call iddeco(ic,jc)
      else
      call idtrb(ibptl(1,i),ibptl(2,i),ibptl(3,i),ibptl(4,i),jc)
      endif
      na=0
      nq=0
      do 53 n=1,nflavt
      na=na+jc(n,1)+jc(n,2)
53    nq=nq+jc(n,1)-jc(n,2)
      ns=   jc(3,1)-jc(3,2)
      return
      end

c
c-----------------------------------------------------------------------
      integer function idraflx(proba,xxx,qqs,icl,jc,jcval,j,iso,c)
c-----------------------------------------------------------------------
c     returns random flavor, according to jc and GRV structure function
c for x(=xxx) and Q2(=qqs) for valence quark (jcval) and sea quarks
c and update jc with quark-antiquark cancellation
c             jc : quark content of remnant
c     j=1 quark, j=2 antiquark,
c     iso : isospin
c     proba : probability of the selected quark (output)
c     c     : "v" is to choose a valence quark, "s" a sea or valence quark
c-----------------------------------------------------------------------
      include 'epos.inc'
      integer jc(nflav,2),jcval(nflav,2)
      double precision s,puv,pdv,psv,pcv,pus,pds,pss,pcs,piso,proba
     &                ,r,drangen
      character c*1

       if(ish.ge.8)then
         write(ifch,10)c,j,xxx,qqs,jc,jcval
       endif
 10    format('entry idraflx, j,x,q: ',a1,i2,2g13.5,/
     &  ,15x,'jc :',2(1x,6i2),/,14x,'jcv :',2(1x,6i2))

       puv=dble(jcval(1,j))
       pdv=dble(jcval(2,j))
       psv=dble(jcval(3,j))
       pcv=dble(jcval(4,j))
       if(c.eq."v")then
         pus=0d0
         pds=0d0
         pss=0d0
         pcs=0d0
       else
         pus=dble(jc(1,j))-puv
         pds=dble(jc(2,j))-pdv
         pss=dble(jc(3,j))-psv
         pcs=dble(jc(4,j))-pcv
       endif

      if(ish.ge.8)then
        write(ifch,'(a,4f6.3)')'idraflx valence:',puv,pdv,psv,pcv
        write(ifch,'(a,4f6.3)')'idraflx sea:',pus,pds,pss,pcs
      endif

       qq=0.
       if(iso.gt.0)then
         if(icl.eq.2)puv=puv*0.5d0   !because GRV already take into account the fact that there is 2 u quark in a proton
         puv=puv*dble(psdfh4(xxx,qqs,qq,icl,1))
         pus=pus*dble(psdfh4(xxx,qqs,qq,icl,-1))
         pdv=pdv*dble(psdfh4(xxx,qqs,qq,icl,2))
         pds=pds*dble(psdfh4(xxx,qqs,qq,icl,-2))
       elseif(iso.lt.0)then
         puv=puv*dble(psdfh4(xxx,qqs,qq,icl,2))
         pus=pus*dble(psdfh4(xxx,qqs,qq,icl,-2))
         if(icl.eq.2)pdv=pdv*0.5d0
         pdv=pdv*dble(psdfh4(xxx,qqs,qq,icl,1))
         pds=pds*dble(psdfh4(xxx,qqs,qq,icl,-1))
       else
         piso=(dble(psdfh4(xxx,qqs,qq,icl,1))
     &        +dble(psdfh4(xxx,qqs,qq,icl,2)))
         if(icl.eq.2)then   !3 quarks
           piso=piso/3d0
         else               !2 quarks
           piso=piso/2d0
         endif
         puv=puv*piso
         pdv=pdv*piso
         piso=0.5d0*(dble(psdfh4(xxx,qqs,qq,icl,-1))
     &              +dble(psdfh4(xxx,qqs,qq,icl,-2)))
         pus=pus*piso
         pds=pds*piso
       endif
       psv=psv*dble(psdfh4(xxx,qqs,qq,icl,3))
       pss=pss*dble(psdfh4(xxx,qqs,qq,icl,-3))
       if(nrflav.ge.4)then
         pcv=pcv*dble(psdfh4(xxx,qqs,qq,icl,4))
         pcs=pcs*dble(psdfh4(xxx,qqs,qq,icl,-4))
       else
         pcv=0d0
         pcs=0d0
       endif

      if(ish.ge.8)then
        write(ifch,'(a,4f6.3)')'idraflx P(valence):',puv,pdv,psv,pcv
        write(ifch,'(a,4f6.3)')'idraflx P(sea):',pus,pds,pss,pcs
      endif

      s=puv+pdv+psv+pcv+pus+pds+pss+pcs
      i=0
      if(s.gt.0.)then
       r=drangen(s)*s
       if(r.gt.(pdv+pus+pds+pss+psv+pcv+pcs).and.puv.gt.0.)then
        i=1
        jcval(i,j)=jcval(i,j)-1
        proba=puv
       elseif(r.gt.(pus+pds+pss+psv+pcv+pcs).and.pdv.gt.0.)then
        i=2
        jcval(i,j)=jcval(i,j)-1
        proba=pdv
       elseif(r.gt.(pds+pss+psv+pcv+pcs).and.pus.gt.0.)then
        i=1
        proba=pus
       elseif(r.gt.(pss+psv+pcv+pcs).and.pds.gt.0.)then
        i=2
        proba=pds
       elseif(r.gt.(psv+pcv+pcs).and.pss.gt.0.)then
        i=3
        proba=pss
       elseif(r.gt.(pcv+pcs).and.psv.gt.0.)then
        i=3
        jcval(i,j)=jcval(i,j)-1
        proba=psv
       elseif(r.gt.pcs.and.pcv.gt.0.)then
        i=4
        jcval(i,j)=jcval(i,j)-1
        proba=pcv
       elseif(pcs.gt.0.)then
        i=4
        proba=pcs
       else
        call utstop("Problem in idraflx, should not be !&")
       endif
      else
        i=idrafl(icl,jc,j,"v",0,iretso)      !no update of jc here
        if(jc(i,j)-jcval(i,j).lt.1)jcval(i,j)=jcval(i,j)-1
        proba=0d0
      endif
      idraflx=i

      if(ish.ge.8)then
        write(ifch,'(a,2(1x,6i2))')'jc before updating:',jc
        write(ifch,20)i,j,jcval,proba
      endif
 20   format('i,j|jcval|P:',2i3,' |',2(1x,6i2),' |',g15.3)

      call idsufl3(i,j,jc)

      if(ish.ge.8)
     & write(ifch,'(a,2(1x,6i2))')'jc after updating:',jc

      return
      end

c-----------------------------------------------------------------------
      integer function idrafl(icl,jc,j,c,imod,iretso)
c-----------------------------------------------------------------------
c     returns random flavor,
c     if : c='v' : according to jc
c          c='s' : from sea
c          c='r' : from sea (always without c quarks)
c          c='d' : from sea for second quark in diquark
c          c='c' : take out c quark first
c             jc : quark content of remnant
c     j=1 quark, j=2 antiquark,
c     imod=0     : returns random flavor of a quark
c     imod=1     : returns random flavor of a quark and update jc
c                 (with quark-antiquark cancellation)
c     imod=2     : returns random flavor of a quark and update jc
c                 (without quark-antiquak cancellation -> accumulate quark)
c     imod=3     : returns random flavor of a quark and update jc with
c                  the corresponding quark-antiquark pair
c                 (without quark-antiquak cancellation -> accumulate quark)
c
c     iretso=0   : ok
c           =1   : more than 9 quarks of same flavor attempted
c-----------------------------------------------------------------------
      include 'epos.inc'
      integer jc(nflav,2),ic(2)
      character c*1

      if(ish.ge.8)
     & write(ifch,*)'entry idrafl, j,imod,c: ',j,imod,' ',c

      pui=1.
      if(c.eq.'s')then
        pu=pui
        pd=pui*exp(-pi*difud/fkappa)
        ps=pui*exp(-pi*difus/fkappa)
        pc=pui*exp(-pi*difuc/fkappa)
        pu=rstrau(icl)*pu
        pd=rstrad(icl)*pd
        ps=rstras(icl)*ps
        pc=rstrac(icl)*pc
      elseif(c.eq.'d')then
        pu=pui*exp(-pi*difuuu/fkappa)
        pd=pui*exp(-pi*difudd/fkappa)
        ps=pui*exp(-pi*difuss/fkappa)
        pc=pui*exp(-pi*difucc/fkappa)
        pu=pu*rstrau(icl)
        pd=pd*rstrad(icl)
        ps=ps*rstras(icl)
        pc=pc*rstrac(icl)
      elseif(c.eq.'v')then
        pu=float(jc(1,j))
        pd=float(jc(2,j))
        ps=float(jc(3,j))
        pc=float(jc(4,j))
      elseif(c.eq.'r')then
        pu=pui
        pd=pui*exp(-pi*difud/fkappa)
        ps=pui*exp(-pi*difus/fkappa)
        pc=0.
        pu=rstrau(icl)*pu
        pd=rstrad(icl)*pd
        ps=rstras(icl)*ps
      elseif(c.eq.'c')then
        pu=0.
        pd=0.
        ps=0.
        pc=1.
      else
        stop'idrafl: dunnowhatodo'
      endif


      s=pu+pd+ps+pc
      if(ish.ge.8)
     & write(ifch,*)'idrafl',s,pu,pd,ps,pc
      if(s.gt.0.)then
       r=rangen()*s
       if(r.gt.(pu+pd+ps).and.pc.gt.0d0)then
        i=4
       elseif(r.gt.(pu+pd).and.ps.gt.0d0)then
        i=3
       elseif(r.gt.pu.and.pd.gt.0d0)then
        i=2
       else
        i=1
       endif
      elseif(iremn.le.1.or.c.ne.'v')then
        i=1+min(2,int((2.+rstras(icl))*rangen()))
      else
        idrafl=0
        return
      endif
      idrafl=i

      if(ish.ge.8)then
        write(ifch,*)'jc before updating',jc
        write(ifch,*)'i,j,jc',i,j
      endif

      if(imod.eq.1)then
        if(iLHC.eq.0.and.iremn.eq.2)then
          call idsufl3(i,j,jc)
c   be sure that jc is not empty
          if(jc(i,j).eq.0)then
            call idenco(jc,ic,iret)
            if(iret.eq.0.and.ic(1).eq.0.and.ic(2).eq.0)then
              jc(i,j)=jc(i,j)+1
              jc(i,3-j)=jc(i,3-j)+1
              iretso=1
            endif
          endif
        elseif(iremn.ge.2)then
          call idsufl3(i,j,jc)
        else
          call idsufl(i,j,jc,iretso)
          if(iretso.ne.0.and.ish.ge.2)then
            call utmsg('idrafl&')
            write(ifmt,*)'iret none 0 in idrafl',iretso
            write(ifch,*)'iret none 0 in idrafl',iretso
            call utmsgf
          endif
        endif
      elseif(imod.eq.2)then
        call idsufl2(i,j,jc)    !do not cancel antiquarks with quarks
      elseif(imod.eq.3)then
        call idsufl2(i,1,jc)    !do not cancel antiquarks with quarks
        call idsufl2(i,2,jc)    !do not cancel antiquarks with quarks
      endif


      if(ish.ge.8)
     & write(ifch,*)'jc after updating',jc

      return
      end


c-----------------------------------------------------------------------
      integer function idraflz(jc,j)
c-----------------------------------------------------------------------
      include 'epos.inc'
      integer jc(nflav,2)

      pu=float(jc(1,j))
      pd=float(jc(2,j))
      ps=float(jc(3,j))
      pc=float(jc(4,j))

      s=pu+pd+ps+pc
      if(s.gt.0.)then
       r=rangen()*s
       if(r.gt.(pu+pd+ps).and.pc.gt.0d0)then
        i=4
       elseif(r.gt.(pu+pd).and.ps.gt.0d0)then
        i=3
       elseif(r.gt.pu.and.pd.gt.0d0)then
        i=2
       else
        i=1
       endif
      else
       stop'in idraflz (1)                      '
      endif
      idraflz=i

      if(jc(i,j).lt.1)stop'in idraflz (2)              '
      jc(i,j)=jc(i,j)-1

      end

c-----------------------------------------------------------------------
      subroutine idsufl(i,j,jc,iretso)
c-----------------------------------------------------------------------
c subtract flavor i, j=1 quark, j=2 antiquark
c add antiflavor if jc(i,j)=0
c iretso=0  ok
c       =1 : more than 9 quarks of same flavor attempted
c-----------------------------------------------------------------------
      integer jc(6,2),ic(2)

      if(jc(i,j).gt.0)then
       jc(i,j)=jc(i,j)-1
       call idenco(jc,ic,iret)
       if(ic(1).eq.0.and.ic(2).eq.0)then
         jc(i,j)=jc(i,j)+1
         if(jc(i,3-j).lt.9.and.iret.eq.0)then
           jc(i,3-j)=jc(i,3-j)+1
         else
           iretso=1
         endif
       endif
      else
        if(j.eq.1)then
          if(jc(i,2).lt.9)then
            jc(i,2)=jc(i,2)+1
          else
            iretso=1
          endif
        else
          if(jc(i,1).lt.9)then
            jc(i,1)=jc(i,1)+1
          else
            iretso=1
          endif
        endif
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine idsufl2(i,j,jc)
c-----------------------------------------------------------------------
c substract flavor i, by adding antiquark i, j=1 quark, j=2 antiquark
c Can replace idsufl if we don't want to cancel quarks and antiquarks
c
c Warning : No protection against jc(i,j)>9 ! should not encode jc without test
c
c-----------------------------------------------------------------------
      parameter(nflav=6)
      integer jc(nflav,2)

      jc(i,3-j)=jc(i,3-j)+1

      return
      end

c-----------------------------------------------------------------------
      subroutine idsufl3(i,j,jc)
c-----------------------------------------------------------------------
c subtract flavor i, j=1 quark, j=2 antiquark
c add antiflavor if jc(i,j)=0
c
c Warning : No protection against jc(i,j)>9 ! should not encode jc without test
c
c-----------------------------------------------------------------------
      parameter(nflav=6)
      integer jc(nflav,2)

      if(jc(i,j).gt.0)then
        jc(i,j)=jc(i,j)-1
      else
        jc(i,3-j)=jc(i,3-j)+1
      endif

      return
      end

cc-----------------------------------------------------------------------
c      subroutine idchfl(jc1,jc2,iret)
cc-----------------------------------------------------------------------
cc checks whether jc1 and jc2 have the same number of quarks and antiquarks
cc if yes: iret=0, if no: iret=1
cc-----------------------------------------------------------------------
c      integer jc1(6,2),jc2(6,2)
c
c      iret=0
c
c      do j=1,2
c       n1=0
c       n2=0
c       do i=1,6
c        n1=n1+jc1(i,j)
c        n2=n2+jc2(i,j)
c       enddo
c       if(n1.ne.n2)then
c        iret=1
c        return
c       endif
c      enddo
c
c      end
c
c
c-----------------------------------------------------------------------
      subroutine idres(idi,am,idr,iadj,iii)
c     returns resonance id idr corresponding to mass am.
c     performs mass adjustment, if necessary (if so iadj=1, 0 else)
c     iii=0 do not mix rho0 and omega
c     iii=-1 mix pions and rhos
c     (only for mesons and baryons, error (stop) otherwise)
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mxindx=9900,mxremax=100,mxmamax=30,mxmxmax=30)
      common/crema/indy(mxindx),rema(mxremax,mxmamax),
     * rewi(mxremax,mxmamax),icre2(mxremax,mxmamax)
     *,idre(mxremax,mxmamax),icre1(mxremax,mxmamax)
     *,mxre,mxma,mxmx 
      character cad*10

      write(cad,'(i10)')idi
      iadj=0
      idr=0
      if(abs(idi).lt.20)then
        idr=idi
        return
      endif
      if(abs(am).lt.1.e-5)am=1e-5
      id=idi
      ami=am
      if(am.lt.0.)then
        call idmass(id,am)
        iadj=1
        if(am.le.0.)then
        write(ifch,*)'*****  warning in idres (0): '
     *,'neg mass returned from idmass'
        write(ifch,*)'id,am(input):',idi,ami
        am=1e-5
        endif
      endif

      if(id.eq.0)goto 9999
      if(abs(id).eq.20)id=sign(230,idi)
      i=indy(abs(id))
      if(i.lt.1.or.i.gt.mxre)then
        write(ifmt,'(a,i8,f10.2,2x,a)')'ERROR idres;',id,am
     .  ,'  particle not in table; see check'
        do n=1,nptl
          write(ifch,'(6i10)')iorptl(n),jorptl(n),n
     .    ,ifrptl(1,n),ifrptl(2,n),idptl(n)
        enddo
        stop
      endif

      do j=1,mxma-1
        remin=rema(i,j)
        remax=rema(i,j+1)
c        if(iii.eq.-1)then
c          if(j.eq.1)then
c            if(i.eq.31)remax=0.3
c            if(i.eq.33)remax=0.4
c          elseif(j.eq.2)then
c            if(i.eq.31)remin=0.3
c            if(i.eq.33)remin=0.4
c          endif
c        endif
c        if(am.ge.rema(i,j).and.am.le.rema(i,j+1))then
        if(am.ge.remin.and.am.le.remax)then
          goto 2
        endif
      enddo
      goto 9999
2     continue

c for NA61 and NA22 data on rho spectra (change pions into rho with 1/2 probability)
c      if(iii.eq.-1)then
c        if(j.eq.1.and.(i.eq.31.or.i.eq.33))then
cc        if(j.eq.1.and.i.eq.31)then     !data favor pi0->rho change only
cc          if(rangen().le.0.5)j=2
c          j=2
c        endif
c      endif

      idr=idre(i,j+1)*id/iabs(id)

      iy=mod(iabs(idr),10)
        if(iy.gt.maxres)then
        iadj=0
        idr=0
        goto 9999
      endif

      if(iy.ne.0.and.iy.ne.1)goto 9999

      call idmass(idr,am)
      if(am.lt.0.)then
      write(*,*)'*****  error in idres: '
     *,'neg mass returned from idmass'
      write(*,*)'id,am(input):',idi,ami
      write(*,*)'idr,am:',idr,am
      call utstop('idres: neg mass returned from idmass&')
      endif
      call idwidth(idr,wi)
c      del=max(1.e-3,2.*wi)    !pb with some resonances
      del=max(1.e-3,wi)
      if(abs(ami-am).gt.del)iadj=1
c      write(ifch,*)'res:',id,idr,ami,am,wi,iadj

9999  continue
      if(iii.eq.1)then
        if(idi.eq.221)stop'\n\n     STOP in idres (1) \n\n'
        if(idr.eq.221)stop'\n\n     STOP in idres (2) \n\n'
        if(idr.eq.111)then
          if(rangen().le.0.5)idr=221
          call idmass(idr,am)
        elseif(idr.eq.330)then
          if(rangen().le.0.5)idr=702
          call idmass(idr,am)
        endif
      endif
      if(.not.(ish.ge.8))return
      write(ifch,*)'return from idres. id,ami,am,idr,iadj:'
      write(ifch,*)idi,ami,am,idr,iadj
      return
      end

cc-----------------------------------------------------------------------
c      integer function idsgl(ic,gen,cmp)
cc     returns 1 for singlets (qqq or qqbar) 0 else.
cc-----------------------------------------------------------------------
c      parameter (nflav=6)
c      integer ic(2),jcx(nflav,2),icx(2)
c      character gen*6,cmp*6
c      idsgl=0
c      if(cmp.eq.'cmp-ys')then
c      call idcomi(ic,icx)
c      else
c      icx(1)=ic(1)
c      icx(2)=ic(2)
c      endif
c      call iddeco(icx,jcx)
c      nq=0
c      na=0
c      do 1 i=1,nflav
c      nq=nq+jcx(i,1)
c1     na=na+jcx(i,2)
c      if(nq.eq.0.and.na.eq.0)return
c      if(gen.eq.'gen-no')then
c      if(nq.eq.3.and.na.eq.0.or.nq.eq.1.and.na.eq.1
c     *.or.nq.eq.0.and.na.eq.3)idsgl=1
c      elseif(gen.eq.'gen-ys')then
c      if(mod(nq-na,3).eq.0)idsgl=1
c      endif
c      return
c      end
c

c-----------------------------------------------------------------------
      subroutine idtaustatus(id,p4,p5,taugm,istatus)
c     returns lifetime(c*tau(fm))*gamma for id with energy p4, mass p5
c-----------------------------------------------------------------------
      include 'epos.inc'
      common/creswi/reswi !for BW
      integer idum
      idum=istatus
c------------------------------------
c      parameter (mxindx=9900,mxremax=100,mxmamax=30,mxmxmax=30)
c      common/crema/indy(mxindx),rema(mxremax,mxmamax),
c     * rewi(mxremax,mxmamax),icre2(mxremax,mxmamax)
c     *,idre(mxremax,mxmamax),icre1(mxremax,mxmamax)
c     *,mxre,mxma,mxmx 
c           if(iabs(id).eq.14)then
c      wi=.197/658.654e15
c           elseif(iabs(id).eq.16)then
c      wi=.197/87.11e9
c           elseif(id.eq.-20)then
c      wi=.197/15.34e15
c           elseif(id.eq.20)then
c      wi=.197/2.675e13
c           elseif(id.eq.221)then
c      wi=0.00849                         !<=width omega
c           elseif(abs(id).lt.100.or.id.gt.1e9)then
c      wi=0
c           elseif(iabs(id).lt.1e8)then
c      ii=indy(abs(id))
c      do 75 ima=2,mxma
c      if(iabs(id).eq.idre(ii,ima))then
c        jj=ima
c        goto 75
c      endif
c75    continue
c      if(ii.lt.1.or.ii.gt.mxre.or.jj.lt.1.or.jj.gt.mxma)then
c      write(ifch,*)'id,ii,jj:',id,'   ',ii,jj
c      call utstop('idtau: ii or jj out of range&')
c      endif
c      wi=rewi(ii,jj)
c           else
c      tauz=1
c-c   tauz=amin1(9./p5**2,tauz)
c-c   tauz=amax1(.2,tauz)
c      wi=.197/tauz
c           endif
c     wi2=wi
c----------------------------
      call idwidth(id,wi)
c----------------------------
c      dif=abs(wi2-wi)
c      if(wi.gt.0)dif=dif/wi
c      if(dif.gt.1e-3)print*,'idtau',id,wi,wi2,'     ',dif 
c----------------------------
      if(wi.eq.0.)then
        if(abs(id).lt.1200)then
          tau=ainfin
        else
          tau=0.
        endif
      else
      tau=.197/wi
      endif
      if(p5.ne.0.)then
      gm=p4/p5
      else
      gm=ainfin
      endif
      if(tau.ge.ainfin.or.gm.ge.ainfin)then
      taugm=ainfin
      else
      taugm=tau*gm
      endif
      reswi=wi  !for BW
      return
      end

c-----------------------------------------------------------------------
      subroutine idtau(id,p4,p5,taugm)
c     wrapper function for idtaustatus that also checks whether IDs exist
c     returns lifetime(c*tau(fm))*gamma for id with energy p4, mass p5
c     fails when index out of range. idtaustatus will return negative
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (indexOutOfRange=-999)
      parameter (iijjOutOfRange=-9999)
      istatus=0

      call idtaustatus(id,p4,p5,taugm,istatus)
      if(istatus.eq.indexOutOfRange) then
        write(ifch,*)"id:",id,p4,p5
c        call alistf('tau&')
        call utstop('idtau: ix out of range.&')
      endif
c      if(istatus.eq.iijjOutOfRange) then
c         write(ifch,*)'id,ii,jj:',id,'   ',ii,jj
c         call utstop('idtau: ii or jj out of range&')
c      endif
      return
      end


c-----------------------------------------------------------------------
      integer function idtra(ic,ier,ires,imix)
c-----------------------------------------------------------------------
c     tranforms from werner-id to paige-id
c         ier .... error message (1) or not (0) in case of problem
c         ires ... dummy variable, not used  any more !!!!
c         imix ... 1 not supported any more
c                  2 010000 010000 -> 110, 001000 000100 -> 110
c                  3 010000 010000 -> 110, 001000 000100 -> 220
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (nidt=59)
      integer idt(3,nidt),ic(2)!,icm(2)
      data idt/
     * 100000,100000, 110   ,100000,010000, 120   ,010000,010000, 220
     *,100000,001000, 130   ,010000,001000, 230   ,001000,001000, 330
     *,100000,000100, 140   ,010000,000100, 240   ,001000,000100, 340
     *,000100,000100, 440
     *,100000,000010, 150   ,010000,000010, 250   ,001000,000010, 350
     *,000100,000010, 450   ,000010,000010, 550
     *,100000,000000,   1   ,010000,000000,   2   ,001000,000000,   3
     *,000100,000000,   4   ,000010,000000,   5   ,000001,000000,   6
     *,330000,000000,  17   ,450000,000000,  18   ,660000,000000,  19
     *,200000,000000,1100   ,110000,000000,1200   ,020000,000000,2200
     *,101000,000000,1300   ,011000,000000,2300   ,002000,000000,3300
     *,100100,000000,1400   ,010100,000000,2400   ,001100,000000,3400
     *,000200,000000,4400   ,100010,000000,1500   ,010010,000000,2500  
     *,001010,000000,3500   ,000110,000000,4500   ,000020,000000,5500 
     *,300000,000000,1111   ,210000,000000,1120   ,120000,000000,1220
     *,030000,000000,2221   ,201000,000000,1130   ,111000,000000,1230
     *,021000,000000,2230   ,102000,000000,1330   ,012000,000000,2330
     *,003000,000000,3331   ,200100,000000,1140   ,110100,000000,1240
     *,020100,000000,2240   ,101100,000000,1340   ,011100,000000,2340
     *,002100,000000,3340   ,100200,000000,1440   ,010200,000000,2440
     *,001200,000000,3440   ,000300,000000,4441/

      idtra=0
      if(ic(1).eq.0.and.ic(2).eq.0)return
      i=1
      do while(i.le.nidt.and.idtra.eq.0)
        if(ic(2).eq.idt(1,i).and.ic(1).eq.idt(2,i))idtra=-idt(3,i)
        if(ic(1).eq.idt(1,i).and.ic(2).eq.idt(2,i))idtra=idt(3,i)
        i=i+1
      enddo
      isi=1
      if(idtra.ne.0)isi=idtra/iabs(idtra)

      jspin=0

      if(imix.eq.1)stop'imix=1 no longer supported'
      if(imix.eq.2)then
      if(idtra.eq.220)idtra=110
      if(idtra.eq.330)idtra=110
      elseif(imix.eq.3)then
      if(idtra.eq.220)idtra=110
      if(idtra.eq.330)idtra=220
      endif

      if(idtra.ne.0)idtra=idtra+jspin*isi

      if(idtra.ne.0)return
      if(ier.ne.1)return
      write(ifch,*)'idtra: ic = ',ic,ires
      call utstop('idtra: unknown code&')

      entry idtrai(num,id,ier)
      idtrai=0
      if(iabs(id).eq.20)then
        j=5
c      elseif(iabs(id).eq.110.or.iabs(id).eq.220)then
c        j=1+2*int(2.*rangen())
      else
        j=0
        do i=1,nidt
          if(iabs(id).eq.idt(3,i))then
            j=i
            goto 2
          endif
        enddo
 2      continue
      endif
      if(j.ne.0)then
        if(id.lt.0)then
          idtrai=idt(3-num,j)
        else
          idtrai=idt(num,j)
        endif
        return
      endif
      if(ier.ne.1)return
      write(ifch,*)'idtrai: id = ',id
      call utstop('idtrai: unknown code&')
      end

c-----------------------------------------------------------------------
      subroutine idtrb(ib1,ib2,ib3,ib4,jc)
c     id transformation ib -> jc
c-----------------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2)
      jc(1,1)=ib1/10**4
      jc(2,1)=ib2/10**4
      jc(3,1)=ib3/10**4
      jc(4,1)=ib4/10**4
      jc(5,1)=0
      jc(6,1)=0
      jc(1,2)=mod(ib1,10**4)
      jc(2,2)=mod(ib2,10**4)
      jc(3,2)=mod(ib3,10**4)
      jc(4,2)=mod(ib4,10**4)
      jc(5,2)=0
      jc(6,2)=0
      return
      end

c-----------------------------------------------------------------------
      subroutine idtrbi(jc,ib1,ib2,ib3,ib4)
c     id transformation jc -> ib
c-----------------------------------------------------------------------
      include 'epos.inc'
      integer jc(nflav,2)
      ib1=jc(1,1)*10**4+jc(1,2)
      ib2=jc(2,1)*10**4+jc(2,2)
      ib3=jc(3,1)*10**4+jc(3,2)
      ib4=jc(4,1)*10**4+jc(4,2)
      ib5=jc(5,1)*10**4+jc(5,2)
      ib6=jc(6,1)*10**4+jc(6,2)
      if(ib5.ne.0.or.ib6.ne.0)then
      write(ifch,*)'***** error in idtrbi: bottom or top quarks'
      write(ifch,*)'jc:'
      write(ifch,*)jc
      call utstop('idtrbi: bottom or top quarks&')
      endif
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cgs following new subroutines based on tables idt.dt and idresi.dt 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c-----------------------------------------------------------------------
      subroutine readidtable 
c     read the table ids.dt which contains all id 
c     information for particles in EPOS
c-----------------------------------------------------------------------

      include 'epos.inc'

      character*130 line

      do n=-999,9900
        nlidtbl(n)=0
      enddo
      do nn=1,nrcode
        do l=1,nidtmax
          idtbl(nn,l)=0
        enddo
      enddo
      do nnn=1,nrcode
        nidtmxx(nnn)=0
      enddo
      do m=1,nidtmax
        ifl1tbl(m)=0
        ifl2tbl(m)=0
        ifl3tbl(m)=0
        jspintbl(m)=0
        chrgtbl(m)=0
        amtbl(m)=0
        mutbl(m)=0
      enddo      

      if(ish.ge.4)write(ifch,'(2a)')'load ',fnnx(1:nfnnx)//'/idt.dt'
      open (95,file= fnnx(1:nfnnx)//'/idt.dt'
     *,STATUS='UNKNOWN')
      nline=0
       
      do
998       read(95,'(a)', end=999)line
          if(line(1:1).eq.'!')goto 998
          jj=0 !count number of words in line
          do jjx=1,120
            if(line(jjx:jjx).eq.' '.and.line(jjx+1:jjx+1).ne.' ')jj=jj+1
          enddo
          if(jj.eq.0)goto 998 

          nline=nline+1
          if(nline.gt.nidtmax)stop'ERROR nidtmax too small'
          read(line,* )idtbl(1,nline), idtbl(2,nline),idtbl(3,nline)
     .    ,idtbl(4,nline),idtbl(5,nline),nametbl(nline),ifl1tbl(nline)
     .    ,ifl2tbl(nline),ifl3tbl(nline),jspintbl(nline)
     .    ,amtbl(nline),chrgtbl(nline),witbl(nline),mutbl(nline)
     .    ,idegtbl(nline),statbl(nline) 

          idepos=idtbl(1,nline)
          if(idepos.ge.-999.and.idepos.le.9900)then
            nlidtbl(idepos)=nline
          endif
c          gs 21/03
          if(mutbl(nline).eq.2)then
            if(idepos.ge.-999.and.idepos.le.9900)then
              nlidtbl(-idepos)=nline
            endif
          endif
c          gs-end

          do i=1,nrcode
             if(idtbl(i,nline).ne.99)nidtmxx(i)=nline
          enddo 

      enddo
      
999   continue
      nidtbl=nline
      close(95)

      !call testNewParticleList

      end

c-----------------------------------------------------------------------
      subroutine getNextParticle
     .           (dname,das,iQ,id,igs,js,iB,iS,last,iflag)
c-----------------------------------------------------------------------
c get the particle properties of the "next" particle in list of hadrons
c  which has been obtaines via readidtable
c-----------------------------------------------------------------------
c in:
c   last = 0 when called for the first time 
c out:
c   dname = nametbl(nline)(1:30)
c   das = mass (double)
c   iQ  = charge
c   id  = epos id code
c   igs = degeneracy
c   js  = excitation counter (former jspin)
c   iB  = baryon number
c   iS  = strangeness 
c   last = -1 for the last particle in the list
c        =  1 otherwise
c   iflag = 0 : no hadron
c         = 1 : basic hadron (liste of 54, the historical hadron list from VENUS) 
c         = 2 : charm=bottom=top=0 
c         = 3 : any hadrons
c-----------------------------------------------------------------------

      include 'epos.inc'

      common/cgetptlidtables/nline,mu
      character dname*30, st*1      
      double precision das
      integer iflav(6)
      do i=1,6
        iflav(i)=0
      enddo
      if(last.eq.0)then
        nline=1
        last=1
        mu=1
      else
        if(mu.eq.1)nline=nline+1
      endif
      isi=3-2*mu
      if(mu.eq.1)dname = nametbl(nline)(1:30)
      if(mu.eq.2)dname = 'anti'//nametbl(nline)(1:26)
      das = amtbl(nline)
      wit = witbl(nline)       !width
      iQ  = isi*nint(chrgtbl(nline))
      id  = isi*idtbl(1,nline) !epos id code
      igs = idegtbl(nline)     !degeneracy
      js  = jspintbl(nline)
      st  = statbl(nline) 
      ifl1=isi*ifl1tbl(nline)
      ifl2=isi*ifl2tbl(nline)
      ifl3=isi*ifl3tbl(nline)
      ia1=abs(ifl1)
      ia2=abs(ifl2)
      ia3=abs(ifl3)
      iB  = ( sign(1,ifl1) + sign(1,ifl2) + sign(1,ifl3) ) / 3
      if(ia1.ne.0) iflav(ia1) = iflav(ia1) + sign(1,ifl1) 
      if(ia2.ne.0) iflav(ia2) = iflav(ia2) + sign(1,ifl2) 
      if(ia3.ne.0) iflav(ia3) = iflav(ia3) + sign(1,ifl3) 
      iS  = -iflav(3)   !strangeness
      jcharm=0
      jbott=0
      jtop=0
      if(ia1.eq.4.or.ia2.eq.4.or.ia3.eq.4)jcharm=1
      if(ia1.eq.5.or.ia2.eq.5.or.ia3.eq.5)jbott=1
      if(ia1.eq.6.or.ia2.eq.6.or.ia3.eq.6)jtop=1
      iflag = 0 !no hadron
      if(ifl2.ne.0.and.ifl3.ne.0)then !hadron
        !database->SetMassRange(0.05, 10.0);SetWidthRange(0., 10.);
        if( das.ge.0.05 !min mass
     .  .and.das.le.10.0 !max mass
     .  .and.wit.ge.0.0  !min width
     .  .and.wit.le.10.0 !max width
     .  .and.(st.eq.'R'.or.st.eq.'D'.or.st.eq.'A'.or.st.eq.'B') !established particles
     .  .and.jbott.eq.0.and.jtop.eq.0  ! b=t=0
     .  )then  
          iflag = 3 !any established b=t=0 hadron within mass and width range
          if(jcharm.eq.0)then
            iflag = 2 !c=b=t=0
            if(js.eq.0.or.js.eq.1)then
              iflag = 1 !basic hadron
            endif
          endif
        endif
      elseif(abs(id).eq.17
     . .and.(st.eq.'R'.or.st.eq.'D'.or.st.eq.'A'.or.st.eq.'B'))then !deuterium (can not do more because of number of quark (as it is))
        iflag=1
      endif
      if(mu.eq.2)then
        mu=1
        if(nline.eq.nidtbl)last=-1 !last particle in list
      else
        mu=mutbl(nline)
        if(mu.eq.1.and.nline.eq.nidtbl)last=-1 !last particle in list
      endif
      end

cc-----------------------------------------------------------------------
c      subroutine testNewParticleList
cc-----------------------------------------------------------------------
c      include 'epos.inc'
c      character dname*30
c      double precision das,das2
c      !load old list
c      write(ifmt,'(a)')'load ptltab, dkytab from '
c      write(ifmt,'(1x,a)')fnnx(1:nfnnx)//'YKt/ptl7.data'
c      write(ifmt,'(1x,a)')fnnx(1:nfnnx)//'YKt/dky7.data'
c
c      stop'comment this stop in order to run the test' 
c 
c      !uncomment the following call in order to run the test
c
cc      call loadptldky(
cc     . fnnx(1:nfnnx)//'YKt/ptl7c.data'//CHAR(0) ,
cc     . fnnx(1:nfnnx)//'YKt/dky7.data'//CHAR(0) ) 
c      write(ifmt,'(80a)')('#',i=1,80)
c      write(ifmt,*)' Nr idEPOS   idPDG     Charge    BaryonNr '
c     .,'Strangeness Degeneracy    Flag' 
c      last=0 
c      n=0 
c      do 
c       call getNextParticle
c     . (dname,das,iQ,id,igs,js,iB,iS,last,iflag)
c       if(iflag.gt.0)then 
c        n=n+1
c        idpdg=idtrafo('nxs','pdg',id)
c        !write(ifmt,*)'NEW',n,id,idpdg
c        !~~~OLD~~~
c        last2=0
c        n2=0
c        do
c        call 
c     .  igetparticle(dname,das2,iQ2,id2,igs2,js2,iB2,iS2,last2,iflag2)
c        if(abs(last2).eq.1)then
c          n2=n2+1
c          !write(ifmt,*),'   OLD',n2,id2
c          if(idpdg.eq.id2)then
c            if(iQ.ne.iQ2.or.iB.ne.iB2.or.iS.ne.iS2
c     .       .or.igs.ne.igs2.or.iflag.ne.iflag2)then
c             write(ifmt,'(i4,i7,i8,5(3x,2i4))')
c     .       n,id,idpdg,iQ,iQ2,iB,iB2,iS,iS2,igs,igs2,iflag,iflag2
c            endif 
c            goto 998
c          endif
c        endif
c        if(last2.lt.0)then
c          write(ifmt,'(2a,i6,i8)')'No particle found corresponding to  '
c     .    ,'idEPOS idPDG =',id,idpdg
c          goto 998
c        endif
c        enddo
c  998   continue 
c       endif
c       if(last.eq.-1) goto 999
c      enddo
c  999 continue
c      write(ifmt,'(a,i5)')'Number of particles chosen:',n  
c      stop'testNewParticleList finished'
c      end
c

c------------------------------------------------------------------------------
      integer function idtrafo(code1,code2,idi)
c------------------------------------------------------------------------------
c wrapper for idtrafo to catch if ID not found
      character*3 code1,code2
      integer istatus=0
      idtrafo = idtrafostatus(code1,code2,idi,istatus)
c      if (istatus.ne.0) then
c         print *,'idtrafo: ',code1,' -> ', code2,idi,' not found.'
c         stop
c      endif
      return
      end

c------------------------------------------------------------------------------
      integer function idtrafostatus(code1,code2,idi,istatus)
     +     result(idtrafo)
c------------------------------------------------------------------------------
c.....tranforms id of code1 (=idi) into id of code2 (=idtrafo)
c.....supported codes:
c.....'nxs' = epos
c.....'pdg' = PDG 1996 (DPMJET)
c.....'qgs' = QGSJet
c.....'III' = QGSJetIII
c.....'ghe' = Gheisha
c.....'sib' = Sibyll
c.....'cor' = Corsika (GEANT)
c.....'flk' = Fluka

c.....returns status 0 if ID could be converted and 1 if it was not found in table

C --- ighenex(I)=EPOS CODE CORRESPONDING TO GHEISHA CODE I ---

      common /ighnx/ ighenex(35)
      data ighenex/
     $               10,   11,   -12,    12,   -14,    14,   120,   110,
     $             -120,  130,    20,   -20,  -130,  1120, -1120,  1220,
     $            -1220, 2130, -2130,  1130,  1230,  2230, -1130, -1230,
     $            -2230, 1330,  2330, -1330, -2330,    17,    18,    19,
     $            3331, -3331,  30/

C --- DATA STMTS. FOR GEANT/GHEISHA PARTICLE CODE CONVERSIONS ---
C --- KIPART(I)=GHEISHA CODE CORRESPONDING TO GEANT   CODE I ---
C --- IKPART(I)=GEANT   CODE CORRESPONDING TO GHEISHA CODE I ---
      DIMENSION        KIPART(48)!,IKPART(35)
      DATA KIPART/
     $               1,   3,   4,   2,   5,   6,   8,   7,
     $               9,  12,  10,  13,  16,  14,  15,  11,
     $              35,  18,  20,  21,  22,  26,  27,  33,
     $              17,  19,  23,  24,  25,  28,  29,  34,
     $              35,  35,  35,  35,  35,  35,  35,  35,
     $              35,  35,  35,  35,  30,  31,  32,  35/

c      DATA IKPART/
c     $               1,   4,   2,   3,   5,   6,   8,   7,
c     $               9,  11,  16,  10,  12,  14,  15,  13,
c     $              25,  18,  26,  19,  20,  21,  27,  28,
c     $              29,  22,  23,  30,  31,  45,  46,  47,
c     $              24,  32,  48/
      INTEGER          ICFTABL(200),IFCTABL(-6:100)
C  ICTABL CONVERTS CORSIKA PARTICLES INTO FLUKA PARTICLES
C  FIRST TABLE ONLY IF CHARMED PARTICLES CAN BE TREATED
      DATA ICFTABL/
     *   7,   4,   3,   0,  10,  11,  23,  13,  14,  12,  ! 10
     *  15,  16,   8,   1,   2,  19,   0,  17,  21,  22,  ! 20
     *  20,  34,  36,  38,   9,  18,  31,  32,  33,  34,  ! 30
     *  37,  39,  24,  25, 6*0,
     *  0,    0,   0,   0,  -3,  -4,  -6,  -5,   0,   0,  ! 50
     *  10*0,
     *   0,   0,   0,   0,   0,   5,   6,  27,  28,   0,  ! 70
     *  10*0,
     *  10*0,
     *  10*0,                                             !100
     *  10*0,
     *   0,   0,   0,   0,   0,  47,  45,  46,  48,  49,  !120
     *  50,   0,   0,   0,   0,   0,   0,   0,   0,   0,  !130
     *  41,  42,  43,  44,   0,   0,  51,  52,  53,   0,  !140
     *   0,   0,  54,  55,  56,   0,   0,   0,  57,  58,  !150
     *  59,   0,   0,   0,  60,  61,  62,   0,   0,   0,  !160
     *  40*0/
C  IFCTABL CONVERTS FLUKA PARTICLES INTO CORSIKA PARTICLES
      DATA IFCTABL/
     *                402, 302, 301, 201,   0,   0,   0,
     *  14,  15,   3,   2,  66,  67,   1,  13,  25,   5,
     *   6,  10,   8,   9,  11,  12,  18,  26,  16,  21,
     *  19,  20,   7,  33,  34,   0,  68,  69,   0,   0,
     *  27,  28,  29,  22,  30,  23,  31,  24,  32,   0,
     * 131, 132, 133, 134, 117, 118, 116, 119, 120, 121,
     * 137, 138, 139, 143, 144, 145, 149, 150, 151, 155,
     * 156, 157,   0,   0,   36*0/
c-------------------------------------------------------------------------------
      include 'epos.inc'

      integer istatus
      character*3 code1,code2
c      double precision drangen,dummy !gs unused after correction with mutbl
      
      istatus=0
      idtrafo=0
      id1=idi
      if(code1.eq.'nxs')then
        i=1
      elseif(code1.eq.'pdg')then
        i=2
      elseif(code1.eq.'qgs')then
        i=3
        if(id1.eq.-10)id1=18
        if(abs(id1).eq.7.or.abs(id1).eq.8)id1=id1*10
      elseif(code1.eq.'III')then
        i=3
      elseif(code1.eq.'cor')then
        i=4
      elseif(code1.eq.'sib')then
        i=5
        if(abs(id1).eq.99)id1=sign(100,id1)
      elseif(code1.eq.'ghe')then
        id1=ighenex(id1)
        i=1
      elseif(code1.eq.'flk')then
        id1=IFCTABL(id1)          !convert to corsika code
        i=4
      else
        stop "unknown code in idtrafo"
      endif
      if(code2.eq.'nxs')then
        j=1
        ji=j
        if(i.eq.2.and.id1.gt.1000000000)then   !nucleus from PDG
          idtrafo=id1
          return
        elseif(i.eq.4.and.id1.gt.402)then               !nucleus from Corsika
          idtrafo=1000000000+mod(id1,100)*10000+(id1/100)*10
          return
        elseif(i.eq.5.and.id1.gt.1004)then               !nucleus from Sibyll
          id1=(id1-1000)
          idtrafo=1000000000+id1/2*10000+id1*10
          return
        elseif(id1.eq.130.and.i.eq.2)then
          idtrafo=-20
          return
        endif
      elseif(code2.eq.'pdg')then
        j=2
        ji=j
        if(i.eq.1.and.id1.gt.1000000000)then !nucleus from NEXUS
          idtrafo=id1
          return
        elseif(i.eq.4.and.id1.gt.402)then               !nucleus from Corsika
          idtrafo=1000000000+mod(id1,100)*10000+(id1/100)*10
          return
        elseif(i.eq.5.and.id1.gt.1004)then               !nucleus from Sibyll
          id1=(id1-1000)
          idtrafo=1000000000+id1/2*10000+id1*10
          return
        elseif(id1.eq.-20.and.i.eq.1)then
          idtrafo=130
          return
        endif
       elseif(code2.eq.'qgs')then
        j=3
        ji=j
       elseif(code2.eq.'III')then
        j=3
        ji=j
      elseif(code2.eq.'cor')then
        j=4
        ji=j
      elseif(code2.eq.'sib')then
        j=5
        ji=j
      elseif(code2.eq.'ghe')then
        j=4
        ji=6
      elseif(code2.eq.'flk')then
        j=4
        ji=7
       else
        stop "unknown code in idtrafo"
      endif
      if(i.eq.4)then !corsika  id always >0 so convert antiparticles
        iadtr=id1
        if(iadtr.eq.25)then
          id1=-13
        elseif(iadtr.eq.15)then
          id1=-14
        elseif(iadtr.ge.26.and.iadtr.le.32)then
          id1=-iadtr+8
        elseif(iadtr.ge.58.and.iadtr.le.61)then
          id1=-iadtr+4
        elseif(iadtr.ge.149.and.iadtr.le.157)then
          id1=-iadtr+12
        elseif(iadtr.ge.171.and.iadtr.le.173)then
          id1=-iadtr+10
        elseif(iadtr.ge.190.and.iadtr.le.195)then
          id1=-iadtr+6
        endif
      endif
      iad1=abs(id1)
      isi=sign(1,id1)

      if(i.ne.j)then
c      nidtmx=max(nidtmxx(i),nidtmxx(j))
       do n=1,nidtmxx(i)
c        write(ifch,*)idtbl(i,n)
c         print *,i,j,n,id1,idtbl(i,n),idtbl(j,n)
        if(id1.eq.idtbl(i,n))then
          idtrafo=idtbl(j,n)
          if(abs(idtrafo).eq.99.and.abs(id1).ne.3340)then
            !print(ifch,*)
            !.      'unknown '//code2//'_id for '//code1//'_id ',idi
            idtrafo=0
          endif
          if(idtrafo.lt.0.and.j.eq.4)then           !corsika  id always >0
            iadtr=abs(idtrafo)
            if(iadtr.eq.13)then
              idtrafo=25
            elseif(iadtr.eq.14)then
              idtrafo=15
            elseif(iadtr.ge.18.and.iadtr.le.24)then
              idtrafo=iadtr+8
            elseif(iadtr.ge.54.and.iadtr.le.57)then
              idtrafo=iadtr+4
            elseif(iadtr.ge.137.and.iadtr.le.145)then
              idtrafo=iadtr+12
            elseif(iadtr.ge.161.and.iadtr.le.163)then
              idtrafo=iadtr+10
            elseif(iadtr.ge.184.and.iadtr.le.189)then
              idtrafo=iadtr+6
            else
              idtrafo=iadtr
            endif
          elseif(abs(idtrafo).eq.100.and.j.eq.5)then
            idtrafo=sign(99,idtrafo)
          elseif(idtrafo.eq.17.and.code2.eq.'qgs')then
            idtrafo=-10
          elseif(abs(idtrafo).eq.80
     .       .or.abs(idtrafo).eq.70.and.code2.eq.'qgs')then
            idtrafo=idtrafo/10
          endif
          if(j.ne.ji)goto 100
          return
        elseif(mutbl(n).eq.2.and.id1.eq.-idtbl(i,n))then
          idtrafo=-idtbl(j,n)
          if(abs(idtrafo).eq.99)then
            !print(ifch,*)
            !.      'unknown '//code2//'_id for '//code1//'_id ',idi
            idtrafo=0
          endif
          if(idtrafo.lt.0.and.j.eq.4)then           !corsika  id always >0
            iadtr=abs(idtrafo)
            if(iadtr.eq.13)then
              idtrafo=25
            elseif(iadtr.eq.14)then
              idtrafo=15
            elseif(iadtr.ge.18.and.iadtr.le.24)then
              idtrafo=iadtr+8
            elseif(iadtr.ge.54.and.iadtr.le.57)then
              idtrafo=iadtr+4
            elseif(iadtr.ge.137.and.iadtr.le.145)then
              idtrafo=iadtr+12
            elseif(iadtr.ge.161.and.iadtr.le.163)then
              idtrafo=iadtr+10
            elseif(iadtr.ge.184.and.iadtr.le.189)then
              idtrafo=iadtr+6
            else
              idtrafo=iadtr
            endif
          elseif(abs(idtrafo).eq.100.and.j.eq.5)then
            idtrafo=sign(99,idtrafo)
          elseif(idtrafo.eq.17.and.code2.eq.'qgs')then
            idtrafo=-10
          elseif(abs(idtrafo).eq.80
     .       .or.abs(idtrafo).eq.70.and.code2.eq.'qgs')then
            idtrafo=idtrafo/10
         endif
          if(j.ne.ji)goto 100
          return
        endif
       enddo
      else
        idtrafo=id1
        if(j.ne.ji)goto 100
        return
      endif

      istatus=1
c      idtrafocx=0
      return

 100  if(j.eq.4)then            !corsika
        if(idtrafo.eq.201)then
          idtrafo=45
        elseif(idtrafo.eq.301)then
          idtrafo=46
        elseif(idtrafo.eq.402)then
          idtrafo=47
        elseif(idtrafo.eq.302)then
          idtrafo=48
        endif
        if(idtrafo.ne.0)then      !air
          if(ji.eq.6)then
            idtrafo=kipart(idtrafo)
          elseif(ji.eq.7)then
            idtrafo=ICFTABL(idtrafo)
          endif
        endif
        return
      else
        call utstop('Should not happen in idtrafo !&')
      endif

      end

c-------------------------------------------------------------------------
      subroutine id_epostopdg(idepos,idpdg)
c-------------------------------------------------------------------------
c tranforms id from EPOS to PDG
c extracted from function idtrafo
c-------------------------------------------------------------------------
      include 'epos.inc'

      if(idepos.gt.1000000000)then !nucleus from NEXUS
        idpdg=idepos
        return
      elseif(idepos.eq.-20)then
        idpdg=130
        return
      endif
      do n=1,nidtmxx(1)
        if(idepos.eq.idtbl(1,n))then
          idpdg=idtbl(2,n)
          if(abs(idpdg).eq.99)then
            !print(ifch,*)
            !.      'unknown '//code2//'_id for '//code1//'_id ',idi
            idpdg=0
          endif
        elseif(mutbl(n).eq.2.and.idepos.eq.-idtbl(1,n))then
          idpdg=-idtbl(2,n)
          if(abs(idpdg).eq.99)then
            !print(ifch,*)
            !.      'unknown '//code2//'_id for '//code1//'_id ',idi
            idpdg=0
          endif
        endif
      enddo
      return

      end

c-----------------------------------------------------------------------
      subroutine idwidth(id,wi)
c     returns the width of the particle with ident code id.
c-----------------------------------------------------------------------

      include 'epos.inc'
                  
      if(id.ge.-9900.and.id.le.9900)then
        nl=0
        if(id.gt.-700)nl=nlidtbl(id)
        if(nl.eq.0)nl=nlidtbl(abs(id))
        wi=witbl(nl)
        return
      endif       
      
      !treat the rare cases outside range -9900 - 9900
      do nl=1,nidtmxx(1)
        if(idtbl(1,nl).eq.id)then
          wi=amtbl(nl)
          return
        endif
      enddo

      !nothing found
      wi=.197/1.0  !tau=1
      return

      end

c-----------------------------------------------------------------------
      subroutine idmass(idi,amass)
c     returns the mass of the particle with ident code id.
c     (deuteron, triton and alpha mass come from Gheisha ???)
c-----------------------------------------------------------------------

      include 'epos.inc'

      id=idi
      if(id.eq.0)then  !?????id????? please change this !!!!!!!!!!!!!
        id=1120                 !for air target
      elseif(abs(id).ge.1000000000)then   !?????id??????????? change this
      !------------------->  LOOK FOR ALL "?????id?????" also in xan
c        print*,id
        goto 600                !nucleus
      endif
            
      if(id.ge.-9900.and.id.le.9900)then
        nl=nlidtbl(abs(id))
        amass=amtbl(nl)  !nl=0 -> should crash with backtrace
        ! print*,'Check mass',id,nl,amass
        if(nl.eq.0)stop'ERROR 07072016' 
        return
      endif       
      
      !treat the rare cases outside range -9900 - 9900
      do nl=1,nidtmxx(1)
        if(idtbl(1,nl).eq.id)then
          amass=amtbl(nl)
          return
        endif
      enddo

      !nothing found
      amass=0
      return
      !stop'####### ERROR 17072016b #######'

600   continue !nuclei
      nbrpro=mod(abs(id/10000),1000)
      nbrneu=mod(abs(id/10),1000)-nbrpro
c        print*,'nuclei',id,nbrpro,nbrneu
      amass=nbrpro*amtbl(nlidtbl(1120))+nbrneu*amtbl(nlidtbl(1220))
      return
      end

c-----------------------------------------------------------------------
      subroutine idflav(id,ifl1,ifl2,ifl3,jspin,idu)
c-----------------------------------------------------------------------
c returns ifl1,ifl2,ifl3,jspin for given id
c  ifl1,ifl2,ifl3: quark flavors
c   baryons: ifl1,ifl2,ifl3
c   mesons: ifl2,ifl3  
c   quarks: ifl3  
c  jspin: excitation counter (0,1,2,3...)
c  idu: dummy (may be used for test purposes)
c-----------------------------------------------------------------------
      
      include 'epos.inc'

      ida=iabs(id)
      idu0=idu !avoid warning
      
      if(ida.le.9900)then
        nl=nlidtbl(abs(id))
        if(nl.eq.0)then
           !~~~~~~~~~~~~~~~~~~~~~~~~ 
           !idu2=iorptl(idu)   !use particle index when calling idflav
           !idu3=iorptl(idu2)
           !print*,'pb idflav',idu,idptl(idu),istptl(idu),ityptl(idu)
           !print*,'         ',idu2,idptl(idu2),istptl(idu2),ityptl(idu2)
           !print*,'         ',idu3,idptl(idu3),istptl(idu3),ityptl(idu3)
           !print*,'             ',ifrptl(1,idu2),ifrptl(2,idu2)
           !~~~~~~~~~~~~~~~~~~~~~~~~
           stop'ERROR 11072016e'
        endif 
        goto 999
      endif  

      !treat the rare cases outside range -9900 - 9900
      do nlx=1,nidtmxx(1)
        if(idtbl(1,nlx).eq.id)then
          nl=nlx
          ida=0 ! to avoid the if... at the end
          goto 999
        endif
      enddo

      !nothing found
      ifl1=0
      ifl2=0
      ifl3=0
      jspin=0
      return

 999  continue

      ifl1=ifl1tbl(nl)
      ifl2=ifl2tbl(nl)
      ifl3=ifl3tbl(nl)
      jspin=jspintbl(nl)
 
      if(id.eq.-ida)then
        ifl1=-ifl1tbl(nl)
        ifl2=-ifl2tbl(nl)
        ifl3=-ifl3tbl(nl)
        jspin=jspintbl(nl)
c        print*,'Check flav 2',id,ifl1,ifl2,ifl3,jspin
      endif       
      
      end

c-------------------------------------------------------------
      subroutine idtr6(id,ic) !id -> ic for via idflav / idt.dt
c-------------------------------------------------------------
      parameter (nflav=6)
      integer jc(nflav,2),ic(2),ifl(3)
      ic(1)=000000
      ic(2)=000000
      if(mod(abs(id),100).eq.99)return !not a particle
      if(iabs(id).gt.16.and.iabs(id).lt.20)then
        if(id.eq.17)then
          ic(1)=330000
          ic(2)=000000
        elseif(id.eq.-17)then
          ic(1)=000000
          ic(2)=330000
        elseif(id.eq.18)then
          ic(1)=450000
          ic(2)=000000
        elseif(id.eq.-18)then
          ic(1)=000000
          ic(2)=450000
        elseif(id.eq.19)then
          ic(1)=660000
          ic(2)=000000
        elseif(id.eq.-19)then
          ic(1)=000000
          ic(2)=660000
        endif
        return
      endif
      if(id.eq.30)then
         ic(1)=222000
         ic(2)=000000
         return
      elseif(id.eq.-30)then
         ic(1)=000000
         ic(2)=222000
         return
      endif
      if(mod(id/10**8,10).eq.8)then
        ic(1)=mod(id,10**8)/10000*100
        ic(2)=mod(id,10**4)*100
        return
      elseif(id/10**9.eq.1.and.mod(id,10).eq.0)then   !nuclei
        nstr=mod(id,10**8)/10000000
        npro=mod(id,10**7)/10000
        nneu=mod(id,10**4)/10
        ic(1)=(2*npro+nneu)*10**5+(2*nneu+npro)*10**4+nstr*10**3
        ic(2)=0
        return
      endif
      do i=1,nflav
        jc(i,1)=0
        jc(i,2)=0
      enddo
      call idflav(id,ifl(1),ifl(2),ifl(3),jspin,idu)
      do j=1,3
         if(ifl(j).gt.0)jc(ifl(j),1)=jc(ifl(j),1)+1
         if(ifl(j).lt.0)jc(-ifl(j),2)=jc(-ifl(j),2)+1
      enddo
      call idenco(jc,ic,iret)
      !print*,'idtr6',idxx,'  ',ifl,'   ',ic
      end

c-----------------------------------------------------------------------
      subroutine idqufl(n,id,nqu,nqd,nqs,nqc)
c     unpacks the ident code of particle (n) and give the number of
c     quarks of each flavour(only u,d,s,c)
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      parameter (nflavt=6)
      integer jc(nflavt,2),ic(2)

      nqu=0
      nqd=0
      nqs=0
      if(iabs(id).ge.7.and.iabs(id).lt.100.and.iabs(id).ne.20)return
      if(iabs(id)/10.eq.11.or.iabs(id)/10.eq.22)return
      if(iabs(id).eq.20)then
        if(iorptl(n).gt.0)then
          if(idptl(iorptl(n)).gt.0)then
            nqd=1
            nqs=-1
          else
            nqd=-1
            nqs=1
          endif
        else
          if(ish.ge.4)write(ifch,*)'Cannot count the number of quark'
        endif
        return
      endif
      if(id.ne.0.and.mod(id,100).eq.0.and.id.le.10**8) goto 300
      if(id/10**8.ne.7)then
        call idtr6(id,ic)
        call iddeco(ic,jc)
      else
        call idtrb(ibptl(1,n),ibptl(2,n),ibptl(3,n),ibptl(4,n),jc)
      endif
      nqu=jc(1,1)-jc(1,2)
      nqd=jc(2,1)-jc(2,2)
      nqs=jc(3,1)-jc(3,2)
      nqc=jc(4,1)-jc(4,2)
      return
 300  i=iabs(id)/1000
      j=mod(iabs(id)/100,10)
      ifl1=isign(i,id)
      ifl2=isign(j,id)
      if(iabs(ifl1).eq.1)nqu=isign(1,ifl1)
      if(iabs(ifl1).eq.2)nqd=isign(1,ifl1)
      if(iabs(ifl1).eq.3)nqs=isign(1,ifl1)
      if(iabs(ifl2).eq.1)nqu=nqu+isign(1,ifl2)
      if(iabs(ifl2).eq.2)nqd=nqd+isign(1,ifl2)
      if(iabs(ifl2).eq.3)nqs=nqs+isign(1,ifl2)
c      write(ifch,*)'id',id,ifl1,ifl2,nqu,nqd,nqs,nqc
      return
      end

c-----------------------------------------------------------------------
      function idlabl(id)
c     returns the character*8 label for the particle id
c-----------------------------------------------------------------------
      include 'epos.inc'
      character*8 idlabl
      real crash(2)
      icrash=3

      nl=0
      ii=0
      if(id.ge.-9900.and.id.le.9900)then
        if(id.ge.-700.and.id.le.9900)then
          nl=nlidtbl(id)
          ii=1
        elseif(id.ge.-9900.and.id.lt.-700)then
          nl=nlidtbl(abs(id))
          ii=2
        endif
        if(nl.eq.0)then
          print*,'id = ',id,nlidtbl(id)
          xx=crash(icrash) !to force crash
        endif
        goto 999
      endif !treat the rare cases outside range -9900 - 9900

      do nlx=1,nidtmxx(1)
        if(idtbl(1,nlx).eq.id)then
          nl=nlx
          ii=1
          goto 999
        endif
      enddo
 
      !nothing found
      idlabl='unknown '
      return

 999  continue

      idlabl=nametbl(nl)(1:8) 
      if(ii.eq.2)idlabl(8:8)='-'//idlabl(1:7)
      length=index(nametbl(nl),' ')-1  !length of string
      if(length.gt.9-ii)idlabl(8:8)='$'  !symbol for truncated

      end

c-----------------------------------------------------------------------
      subroutine idchrg(iii,id,chrg)
c     computes charge of particle with ident code id
c     ichrg must be dimensioned nqlep+12
c-----------------------------------------------------------------------

      include 'epos.inc'
#ifdef CHROMO
Cf2py intent(out) chrg
#endif
      ida=abs(id)
      if(ida.le.9900)then
        nl=nlidtbl(ida)
        if(nl.eq.0)then
          if(iii.ne.998)then
            print*,'ERROR 11072016c id =',id,' not in list.',iii
            stop
          else
            chrg=1e30
            return
          endif
        endif
        chrg=chrgtbl(nl)
        ! print*,'idchrg using table ',id,chrg
        if(id.eq.-ida)chrg=-chrg 
        return
      endif

      !treat the rare cases outside range -9900 - 9900
      do nl=1,nidtmxx(1)
        if(idtbl(1,nl).eq.id)then
          chrg=chrgtbl(nl)
          return
        endif
      enddo

      chrg=0
      !stop'ERROR 12072016c'

      end

c-----------------------------------------------------------------------
      subroutine idresi
c-----------------------------------------------------------------------
c  Initializes /crema/
c  Masses are limits between two states (their average)
c-----------------------------------------------------------------------
c  The former idresi including the data part plus some code to extract
c  The extended table is idresi.dt
c-----------------------------------------------------------------------

      include 'epos.inc'

      real delmrho,delpeta
      integer irasym
      common /cuncertmu/delmrho,delpeta,irasym
      parameter (mxindx=9900,mxremax=100,mxmamax=30,mxmxmax=30)
      common/crema/indy(mxindx),rema(mxremax,mxmamax),
     * rewi(mxremax,mxmamax),icre2(mxremax,mxmamax)
     *,idre(mxremax,mxmamax),icre1(mxremax,mxmamax)
     *,mxre,mxma,mxmx 
      character line*175

      write(ifmt,'(2a)')'load ',fnnx(1:nfnnx)//'/idresi.dt'
      open (98,file= fnnx(1:nfnnx)//'/idresi.dt')

      do i=1,mxindx
        indy(i)=0
      enddo
      do  k=1,mxre
        do  m=1,mxmamax
          rema(k,m)=0
          icre1(k,m)=0
          icre2(k,m)=0
          rewi(k,m)=0
          idre(k,m)=0 
        enddo
      enddo
      
      itab=0
      ntab=0
      mxmx=0
      mxre=0
      mxma=0
      kmx=0

      do
 998    continue
        read(98,'(a)', end=999)line
        if(line(1:7).eq.'!tables')ntab=1
        if(line(1:1).eq.'!')goto 998

        if(ntab.eq.1)then
          
          jj=0 !count number of words in line
          do j=1,150 !up to number of character in one line
            if(line(j+1:j+1).eq.'!')exit
            if(line(j:j).eq.' '.and.line(j+1:j+1).ne.' ')jj=jj+1
          enddo
          if(jj.eq.0)goto 998 
          if(jj.gt.mxmamax-1)stop'ERROR mxmamax too small'           
c          print*,'jj', jj,line
          
          itab=itab+1
          mxre=(itab-1)/3+1
          mtab=mod(itab-1,3)+1
          
          if(mxre.gt.mxremax)stop'ERROR mxremax too small'

          if(mtab.eq.1)then !icre

            read(line,*)ix,ic1,ic2
            if(ix.lt.1.or.ix.gt.mxindx)
     *      stop'idresi: ix out of range.&' 

            do j=2,mxmamax
              icre1(mxre,j)=ic1
              icre2(mxre,j)=ic2
            enddo 
              
            if((ix.eq.111.or.ix.eq.222
     *      .or.ix.eq.333.or.ix.eq.444
     *      .or.ix.eq.555))then !xxx0 does not exist
               icre1(mxre,2)=0
            endif 
          
          elseif (mtab.eq.2)then !idmx

            kmx=0
            if(jj.gt.1)then
              kmx=1
              mxma=max(mxma,jj+1)
              read(line,*)(idre(mxre,j+1),j=1,jj)
              idre(mxre,1)=ix
              do j=1,jj
                if(idre(mxre,j+1).ne.0)then
                  !print*,'idresiA', idre(mxre,j+1),mxre,j+1
                  indy( idre(mxre,j+1) )=mxre
                endif
              enddo
            endif

          elseif (mtab.eq.3)then !rema
            mxma=max(mxma,jj+1)
            read(line,*)(rema(mxre,j+1),j=1,jj)
c           different masses for rho
            if(irasym.eq.0)then
              if(ix.eq.11.or.ix.eq.12)rema(mxre,2)=rema(mxre,2)+delmrho !better to use low mass for rho -> more rhos means less eta in better agreement with e+e- data
            else  !create rho mass asymmetry to get eta, eta' and pion OK without hacas
              if(ix.eq.11)rema(mxre,2)=rema(mxre,2)-delmrho
              if(ix.eq.12)rema(mxre,2)=rema(mxre,2)+delmrho
            endif
            !print*,'rema', ix,(rema(mxre,j+1),mxre,j+1,j=1,jj)
            if(kmx.eq.0)then
              idre(mxre,1)=ix
              do j=1,jj
                idre(mxre,j+1)=idre(mxre,1)*10+j-1
                !print*,'idresiB', idre(mxre,j+1),mxre,j+1
                indy( idre(mxre,j+1) )=mxre
              enddo
            endif

c----------------------------------
c          elseif(mtab.eq.4)then !rewi
c            read(line,*)(rewi(mxre,j+1),j=1,jj)
c---------------------------------
          else
            stop'ERROR mtab too big'
          endif

        endif

      enddo

999   continue
      close(98)

      return
      end


