C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c---------------------------------------------------------------------
      subroutine xiniall
c---------------------------------------------------------------------
      include "epos.inc"
      include "epos.incxan"
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall

      parameter (mxfra=5)
      common/pfra/nfra,ifra(mxfra),ivfra(2,mxhis),itfra(mxtri,mxhis)
     $     ,imofra(3,mxfra),iffra(mxfra),r1fra(3,mxfra),r2fra(3,mxfra)
     $     ,emax(mxfra)
 
      nhis=0
      nfra=0
      istavar=0
      ivnsp=0
      ipairs1=0
      iwww=0
      do n=0,mxxhis
        icorrtrig(n)=0
        ihardevent(n)=0
        ijetfind1(n)=0
        ijetfind2(n)=0
        ifastjet(n)=0
        iphoton(n)=0 
        iCorH(n)=0
        imupurp(n)=0
        ijetevent(n)=0
        icaltrig(n)=0
        imux(n)=0
        iepx(n)=0
        idetphi(n)=0
      enddo
      do n=1,mxhis
        istuse(n)=-999
        do m=1,mxpara
          xpara(m,n)=0
        enddo
        nhisxxx(n)=0
      enddo
      ncontrall=0
      noerrall=0

      end

c---------------------------------------------------------------------
      subroutine xini
c---------------------------------------------------------------------
c  called after beginhisto | bh | beginanalysis
c  returns after endhisto | eh | end 
c---------------------------------------------------------------------
      include "epos.inc"
      include "epos.incho"
      include "epos.incxan"
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall

      parameter (mxfra=5)
      common/pfra/nfra,ifra(mxfra),ivfra(2,mxhis),itfra(mxtri,mxhis)
     $     ,imofra(3,mxfra),iffra(mxfra),r1fra(3,mxfra),r2fra(3,mxfra)
     $     ,emax(mxfra)
      character line*1000,cvar*25
      logical go
      common/nl/noplin  /cnnnhis/nnnhis
      character*1000 cline
      common/cjjj/jjj,cline

      call utpri('xini  ',ish,ishini,5)

      i=1
                                !      iapl=0
                                !      nhis=0
      j=jjj     !-1
      line=cline
                                !      nfra=1
                                !      ifra(1)=iframe
      iapl=0
      if(nfra.eq.0)then
        nfra=1
        ifra(1)=iframe
      endif
      nhis=nhis+1
      nnnhis=nhis
      if(nhis.gt.mxhis)stop'xini: mxhis too small.       '
      noweak(nhis)=0
      ionoerr=0
c      newfra=0
      indfra=1
c      nepfra=0
      inpfra=1
 1    call utword(line,i,j,0)
      if(line(i:j).eq.'application')then !-----------
        call utword(line,i,j,1)
        if(line(i:j).eq.'analysis')then
                                !iapl=0
                                !nhis=nhis+1
                                !newfra=0
                                !indfra=1
                                !nepfra=0
                                !inpfra=1
        else
          iapl=1
        endif
      elseif(line(i:j).eq.'input')then !-----------
        call utword(line,i,j,0)
        if(nopen.ge.0)then
         nopen=nopen+1
         if(nopen.gt.9)stop'too many nested input commands'
         open(unit=20+nopen,file=line(i:j),status='old')
         if(iprmpt.eq.1)iprmpt=-1
        endif
      elseif(line(i:j).eq.'runprogram')then !-----------
        if(iapl.eq.0)then
        else
          goto 9999
        endif
      elseif(line(i:j).eq.'frame'.or.line(i:j).eq.'frame+')then !------
        ifp=1
        if(line(i:j).eq.'frame+')ifp=2
        call utword(line,i,j,1)
        if(line(i:j).eq.'total')then
          nfp=iframe
        elseif(line(i:j).eq.'nucleon-nucleon')then
          nfp=11
        elseif(line(i:j).eq.'target')then
          nfp=12
        elseif(line(i:j).eq.'gamma-nucleon')then
          nfp=21
        elseif(line(i:j).eq.'lab')then
          nfp=22
        elseif(line(i:j).eq.'breit')then
          nfp=23
        elseif(line(i:j).eq.'thrust')then
          nfp=33
        elseif(line(i:j).eq.'sphericity')then
          nfp=32
        else
          nfp=0
          call utstop("Wrong frame in xini !&")
        endif
        go=.true.
        inl=0
        do l=1,nfra
          if(ifra(l).eq.nfp)then
            inl=l
            go=.false.
          endif
        enddo
        if (go) then
          nfra=nfra+1
          inl=nfra
          ifra(nfra)=nfp
        endif
        if(ifp.eq.1)then
          indfra=inl
c          newfra=nfp
          ivfra(1,nhis)=indfra
          ivfra(2,nhis)=indfra
        else
          inpfra=inl
c          nepfra=nfp
        endif
      elseif(line(i:j).eq.'binning')then !-----------
        call utword(line,i,j,1)
        if(line(i:j).eq.'lin')then
          iologb=0
          iocnxb=0
        elseif(line(i:j).eq.'log')then
          iologb=1
          iocnxb=0
        elseif(line(i:j).eq.'clin')then
          iologb=0
          iocnxb=1
        elseif(line(i:j).eq.'clog')then
          iologb=1
          iocnxb=1
        else
          print *, 'what the heck is ',line(i:j),' binning?'
          print *, 'I will use the linear (lin) one'
        endif
      elseif(line(i:j).eq.'setm')then !-----------
        if(iapl.eq.0) then
          print *,"You should use histogram instead of setm, please"
          stop
        endif
      elseif(line(i:j).eq.'set')then !-----------
        call utword(line,i,j,1)
        if(line(i:j).eq.'iologb')then
          call utword(line,i,j,1)
          read(line(i:j),*) iologb
        elseif(line(i:j).eq.'iocnxb')then
          call utword(line,i,j,1)
          read(line(i:j),*) iocnxb
        elseif(line(i:j).eq.'etacut')then
          call utword(line,i,j,1)
          read(line(i:j),*) etacut
        elseif(line(i:j).eq.'nemsi')then
          call utword(line,i,j,1)
          read(line(i:j),*)nemsi
        endif
      elseif(line(i:j).eq.'xpara')then !-----------
        call utword(line,i,j,1)
        read(line(i:j),*)ipara
        if(ipara.gt.mxpara)stop'mxpara too small.         '
        call utword(line,i,j,1)
        read(line(i:j),*)val
        xpara(ipara,nhis)=val
      elseif(line(i:i+5).eq.'xparas'.or.line(i:i+2).eq.'xps')then !-----------
        if(line(i:j).eq.'xparas'.or.line(i:j).eq.'xps')then
          call utword(line,i,j,1)
          ishift=0
          if(line(i:j).eq.'+10')then
           ishift=10 
          elseif(line(i:j).eq.'+15')then
           ishift=15 
          elseif(line(i:j).eq.'+20')then
           ishift=20 
          elseif(line(i:j).eq.'+25')then
           ishift=25 
          elseif(line(i:j).eq.'+30')then
           ishift=30 
          elseif(line(i:j).eq.'+35')then
           ishift=35 
          elseif(line(i:j).eq.'+40')then
           ishift=40 
          elseif(line(i:j).eq.'+45')then
           ishift=45 
          endif
          if(ishift.gt.0) call utword(line,i,j,1)
        elseif(line(i:i+5).eq.'xparas'.and.line(i+6:i+6).eq.'+')then   !case xparas+XX without space
          read(line(i+7:i+8),*)ishift
          call utword(line,i,j,1)
        elseif(line(i:i+2).eq.'xps'.and.line(i+3:i+3).eq.'+')then      !case xps+XX without space
          read(line(i+4:i+5),*)ishift
          call utword(line,i,j,1)
        else
          stop'ERROR xparas: syntax error'
        endif
        read(line(i:j),*,err=99001)ipara
        if(ipara.gt.mxpara)stop'mxpara too small.'
        ii=1
        do while (ii.le.ipara)
         call utword(line,i,j,1)
         if(line(i:i).eq.'>')then
          read(line(i+1:j),*)ii
          do iii=1,ii
            xpara(ishift+iii,nhis)=0
          enddo 
          ii=ii+1
          call utword(line,i,j,1)
         endif 
         read(line(i:j),*)val
         xpara(ishift+ii,nhis)=val
         ii=ii+1
        enddo
      elseif(line(i:j).eq.'echo')then !-----------
        call utword(line,i,j,1)
        if(line(i:j).eq.'on')iecho=1
        if(line(i:j).eq.'off')iecho=0
        if(line(i:j).ne.'on'.and.line(i:j).ne.'off')
     *  stop'invalid option'
      elseif(line(i:j).eq.'weak')then !-----------
        continue !do nothing
      elseif(line(i:j).eq.'noweak')then !-----------
        noweak(nhis)=1
      elseif(line(i:j).eq.'noweak2')then !-----------
        noweak(nhis)=2
      elseif(line(i:j).eq.'noweak3')then !-----------
        noweak(nhis)=3
      elseif(line(i:j).eq.'histogram'
     *       .or.line(i:j).eq.'hi')then !-----------
        nac(nhis)=1
        call utword(line,i,j,1) !xvaria
        cvar='                         '
        cvar(1:j-i+1)=line(i:j)
        call xtrans(cvar,icas,inom,ifrnew,nhis)
        if(inom.eq.-1)then
          if(line(i:i).ge.'0'.and.line(i:i).le.'9')then
            inom=298
            read(line(i:j),*) sval(1,nhis)
          endif
        endif
        icasv(1,nhis)=icas  !from xtrans
        ivar(1,nhis)=inom   !from xtrans
        if(ifrnew.ne.0)then     !check frame for e+e- event
          go=.true.             !shape variables
          do l=1,nfra
            if(ifra(l).eq.ifrnew)then
              indfra=l
              go=.false.        !have it already
            endif
          enddo
          if (go) then
            nfra=nfra+1
            ifra(nfra)=ifrnew
            indfra=nfra
          endif
        endif
        call utword(line,i,j,1) !yvaria
        cvar='                         '
        cvar(1:j-i+1)=line(i:j)
        call xtrans(cvar,icas,inom,ifrnew,nhis)
        icasv(2,nhis)=icas  !from xtrans
        ivar(2,nhis)=inom   !from xtrans
        if(inom.eq.-1)then
          if(line(i:i).ge.'0'.and.line(i:i).le.'9')then
            inom=299
            read(line(i:j),*) sval(2,nhis)
          endif
        endif
        if(inom.eq.-1)then
          ivar(1,nhis)=inom
          icasv(1,nhis)=icas
        endif
        ivfra(1,nhis)=indfra
        ivfra(2,nhis)=indfra

        call utword(line,i,j,1) !normation
        read(line(i:j),*) inorm(nhis)

        call utword(line,i,j,1) !xmin
        if(line(i:i).eq.'z'.and.line(i+2:i+2).eq.'z' !z1z etc
     .   .and.(
     .   line(i+1:i+1).eq.'0'.or.line(i+1:i+1).eq.'1'
     .   .or.line(i+1:i+1).eq.'2'.or.
     .   line(i+1:i+1).eq.'3'.or.line(i+1:i+1).eq.'4'
     .   .or.line(i+1:i+1).eq.'5'.or.
     .   line(i+1:i+1).eq.'6'.or.line(i+1:i+1).eq.'7'
     .   .or.line(i+1:i+1).eq.'8'.or.
     .   line(i+1:i+1).eq.'9'
     .   ))then
         read(line(i+1:i+1),*)ij 
         xmin(nhis)=ij-0.5
        elseif(line(i:i).eq.'z'.and.line(i+3:i+3).eq.'z' !z10z etc
     .   .and.(
     .   line(i+1:i+1).eq.'0'.or.line(i+1:i+1).eq.'1'
     .   .or.line(i+1:i+1).eq.'2'.or.
     .   line(i+1:i+1).eq.'3'.or.line(i+1:i+1).eq.'4'
     .   .or.line(i+1:i+1).eq.'5'.or.
     .   line(i+1:i+1).eq.'6'.or.line(i+1:i+1).eq.'7'
     .   .or.line(i+1:i+1).eq.'8'.or.
     .   line(i+1:i+1).eq.'9')
     .   .and.(
     .   line(i+2:i+2).eq.'0'.or.line(i+2:i+2).eq.'1'
     .   .or.line(i+2:i+2).eq.'2'.or.
     .   line(i+2:i+2).eq.'3'.or.line(i+2:i+2).eq.'4'
     .   .or.line(i+2:i+2).eq.'5'.or.
     .   line(i+2:i+2).eq.'6'.or.line(i+2:i+2).eq.'7'
     .   .or.line(i+2:i+2).eq.'8'.or.
     .   line(i+2:i+2).eq.'9')
     .   )then
         read(line(i+1:i+2),*)ij 
         xmin(nhis)=ij-0.5
        elseif(line(i:j).eq.'egy')then
         if(engy.gt.0)then
          egy=engy
         elseif(ecms.gt.0.)then
          egy=ecms
         elseif(elab.gt.0)then
          call idmass(idproj,apj)
          call idmass(idtarg,atg)
          egy=sqrt( 2*elab*atg+atg**2+apj**2 )
         elseif(ekin.gt.0.)then
          call idmass(idproj,apj)
          call idmass(idtarg,atg)
          egy=sqrt( 2*(ekin+apj)*atg+atg**2+apj**2 )
         elseif(pnll.gt.0.)then
          call idmass(idproj,apj)
          call idmass(idtarg,atg)
          egy=sqrt( 2*sqrt(pnll**2+apj**2)*atg+atg**2+apj**2 )
         else
          stop'pb in xini (1).   '
         endif
         xmin(nhis)=egy-0.5
        elseif(line(i:j).eq.'xp9')then
          xmin(nhis)=xpara(9,nhis)-0.5
        else
         read(line(i:j),*,err=1006) xmin(nhis)
        endif

        call utword(line,i,j,1) !xmax
        if(line(i:i).eq.'z'.and.line(i+2:i+2).eq.'z' !z1z etc
     .   .and.(
     .   line(i+1:i+1).eq.'0'.or.line(i+1:i+1).eq.'1'
     .   .or.line(i+1:i+1).eq.'2'.or.
     .   line(i+1:i+1).eq.'3'.or.line(i+1:i+1).eq.'4'
     .   .or.line(i+1:i+1).eq.'5'.or.
     .   line(i+1:i+1).eq.'6'.or.line(i+1:i+1).eq.'7'
     .   .or.line(i+1:i+1).eq.'8'.or.
     .   line(i+1:i+1).eq.'9'
     .   ))then
         read(line(i+1:i+1),*)ij 
         xmax(nhis)=ij+0.5
        elseif(line(i:i).eq.'z'.and.line(i+3:i+3).eq.'z' !z10z etc
     .   .and.(
     .   line(i+1:i+1).eq.'0'.or.line(i+1:i+1).eq.'1'
     .   .or.line(i+1:i+1).eq.'2'.or.
     .   line(i+1:i+1).eq.'3'.or.line(i+1:i+1).eq.'4'
     .   .or.line(i+1:i+1).eq.'5'.or.
     .   line(i+1:i+1).eq.'6'.or.line(i+1:i+1).eq.'7'
     .   .or.line(i+1:i+1).eq.'8'.or.
     .   line(i+1:i+1).eq.'9')
     .   .and.(
     .   line(i+2:i+2).eq.'0'.or.line(i+2:i+2).eq.'1'
     .   .or.line(i+2:i+2).eq.'2'.or.
     .   line(i+2:i+2).eq.'3'.or.line(i+2:i+2).eq.'4'
     .   .or.line(i+2:i+2).eq.'5'.or.
     .   line(i+2:i+2).eq.'6'.or.line(i+2:i+2).eq.'7'
     .   .or.line(i+2:i+2).eq.'8'.or.
     .   line(i+2:i+2).eq.'9')
     .   )then
         read(line(i+1:i+2),*)ij 
         xmax(nhis)=ij+0.5
        elseif(line(i:j).eq.'egy')then
         if(engy.gt.0)then
          egy=engy
         elseif(ecms.gt.0.)then
          egy=ecms
         elseif(elab.gt.0)then
          call idmass(idproj,apj)
          call idmass(idtarg,atg)
          egy=sqrt( 2*elab*atg+atg**2+apj**2 )
         elseif(ekin.gt.0.)then
          call idmass(idproj,apj)
          call idmass(idtarg,atg)
          egy=sqrt( 2*(ekin+apj)*atg+atg**2+apj**2 )
         elseif(pnll.gt.0.)then
          call idmass(idproj,apj)
          call idmass(idtarg,atg)
          egy=sqrt( 2*sqrt(pnll**2+apj**2)*atg+atg**2+apj**2 )
         else
          stop'pb in xini (2).   '
         endif
         xmax(nhis)=egy+0.5
        elseif(line(i:j).eq.'xp9')then
          xmax(nhis)=xpara(9,nhis)+0.5
        else
         read(line(i:j),*) xmax(nhis)
        endif

        call utword(line,i,j,1) !nbin
        read(line(i:j),*) nbin(nhis)
        do l=1,nbin(nhis)
          bin(l,nac(nhis),nhis)=0.
          zcbin(l,nac(nhis),nhis)=0
        enddo
        lookcontr(nhis)=0
        lookcontrx(nhis)=0
        inoerr(nhis)=0
      elseif(line(i:j).eq.'istuse')then !-----------
        call utword(line,i,j,1) !istuse
        read(line(i:j),*) istuse(nhis)
      elseif(line(i:j).eq.'idcode')then !-----------
        call utword(line,i,j,1) !idcode
        if(line(i:i+2).eq.'995')stop'xini: idcode 995 not supported'
        k1=i
        do while(k1.le.j)
          k2=k1
          do while(line(k2:k2).ne.'/'.and.k2.ne.j+1)
            k2=k2+1
          enddo
          k2=k2-1
          nidcod(nhis)=nidcod(nhis)+1
          read(line(k1:k2),*,err=99002) idcod(nidcod(nhis),nhis)
          idmod(nidcod(nhis),nhis)=.false.
          k1=k2+2
        enddo
      elseif(line(i:j).eq.'idcode+')then !-----------
        stop'xini: idcode+ not supported'
        call utword(line,i,j,1) !idcode
        if(line(i:i+2).eq.'995')stop'xini: idcode 995 not supported'
        nidcod(nhis)=nidcod(nhis)+1
        read(line(i:j),*) idcod(nidcod(nhis),nhis)
        idmod(nidcod(nhis),nhis)=.true.
      elseif(line(i:j).eq.'trigger'.or.line(i:j).eq.'trg')then !-----------
        call utword(line,i,j,1)
        nsingle=0
        if(line(i:j).eq.'-single')then
          nsingle=1
          call utword(line,i,j,1)
        endif
        ntc=1
        imo=1
        ncontr=0
        icontrtyp(nhis)=0
        if(line(i:j).eq.'or'.or.line(i:j).eq.'-or'
     .  .or.line(i:j).eq.'contr')then
          imo=2
          if(line(i:j).eq.'-or')imo=-2
          if(line(i:j).eq.'contr')imo=3
          call utword(line,i,j,1)
          read(line(i:j),*)ztc
          ntc=nint(ztc)
          call utword(line,i,j,1)
          if(imo.eq.3)then
            ncontr=ntc
            ncontrall=ncontrall+ncontr
            if(ncontrall.gt.mxcontr)stop'xini: mxcontr too small.     '
            if(ncontr.gt.mxcnt)stop'xini: mxcnt too small.     '
            lookcontr(nhis)=ncontrall-ncontr+1
            lookcontrx(nhis)=ncontrall
            do nb=1,nbin(nhis)
              do nn=1,ncontr
                bbin(nb,nac(nhis),lookcontr(nhis)-1+nn)=0.d0
                zbbin(nb,nac(nhis),lookcontr(nhis)-1+nn)=0.d0
              enddo
            enddo
            do nn=1,ncontr
                    nccevt(lookcontr(nhis)-1+nn)=0
            enddo
          endif
        endif
        i111=i
        j111=j
        cline=line 
        do n=1,ntc
          if(n.ne.1)then
            if(imo.ne.-2)then ! NOT -or case
              call utword(line,i,j,1) !trigger-name
              i111=i
              j111=j
              cline=line
            endif  
          endif 
          cvar='                         '
          ifp=1
          if(cline(j111:j111).eq.'+')then
            cvar(1:j111-i111+1)=cline(i111:j111-1)
            ifp=2
          else
            cvar(1:j111-i111+1)=cline(i111:j111)
            ifp=1
          endif
          call xtrans(cvar,icas,inom,ifrnew,nhis)
          if(inom.gt.0)then
            ntri(nhis)=ntri(nhis)+1
            if(ntc.eq.1)then
              ntrc(ntri(nhis),nhis)=1
            elseif(n.eq.1)then
              ntrc(ntri(nhis),nhis)=2
            elseif(n.eq.ntc)then
              ntrc(ntri(nhis),nhis)=3
            else
              ntrc(ntri(nhis),nhis)=0
            endif
            if(imo.eq.3)then
              ntrc(ntri(nhis),nhis)=-1
              if(n.eq.1)then
                icontrtyp(nhis)=1+inom/100
              else
                if(icontrtyp(nhis).le.2.or.1+inom/100.le.2)then
                  if(1+inom/100.ne.icontrtyp(nhis))
     &            stop'xini type mismatch'
                endif
              endif
            endif
            icast(ntri(nhis),nhis)=icas  !from xtrans
            itri(ntri(nhis),nhis)=inom   !from xtrans
            if(ifp.eq.1)then
              itfra(ntri(nhis),nhis)=indfra
            else
              itfra(ntri(nhis),nhis)=inpfra
            endif
            xmitrp(ntri(nhis),nhis)=100.
            xmatrp(ntri(nhis),nhis)=100.
            call utword(line,i,j,1) !-----------xmin----------
            if(line(i:j).eq.'inf')then
              xmitri(ntri(nhis),nhis)=1e30
            elseif(line(i:j).eq.'-inf')then
              xmitri(ntri(nhis),nhis)=-1e30
            elseif(line(i:j).eq.'A')then
              xmitri(ntri(nhis),nhis)=maproj
            elseif(line(i:j).eq.'A+1')then
              xmitri(ntri(nhis),nhis)=maproj+1
            elseif(line(i:j).eq.'A+B')then
              xmitri(ntri(nhis),nhis)=maproj+matarg
            elseif(line(i:j).eq.'A+B+1')then
              xmitri(ntri(nhis),nhis)=nptlpt+1
              if(ireadfzo.gt.0)xmitri(ntri(nhis),nhis)=1
            elseif(line(i:j).eq.'lead')then    !leading particle (neads Standard Variable)
              xmitri(ntri(nhis),nhis)=-123456
              istavar=1
            elseif(line(i:j).eq.'jet')then    !jet from fastjet
              xmitri(ntri(nhis),nhis)=nhis*100
              call actifastjet(nhis)
            else
              kk=0
              do k=i+1,j-1
                if(line(k:k).eq.'%')kk=k
              enddo
              if(kk.eq.0)then
                read(line(i:j),*)xmitri(ntri(nhis),nhis)
              else
                read(line(i:kk-1),*)xmitrp(ntri(nhis),nhis)
                read(line(kk+1:j),*)xmitri(ntri(nhis),nhis)
              endif
            endif
            if(nsingle.ne.1)then
             call utword(line,i,j,1) !-----------xmax------------
             if(line(i:j).eq.'inf')then
               xmatri(ntri(nhis),nhis)=1e30
             elseif(line(i:j).eq.'-inf')then
               xmatri(ntri(nhis),nhis)=-1e30
             elseif(line(i:j).eq.'A')then
               xmatri(ntri(nhis),nhis)=maproj
             elseif(line(i:j).eq.'A+1')then
               xmatri(ntri(nhis),nhis)=maproj+1
             elseif(line(i:j).eq.'A+B')then
               xmatri(ntri(nhis),nhis)=maproj+matarg
             elseif(line(i:j).eq.'A+B+1')then
               xmatri(ntri(nhis),nhis)=maproj+matarg+1
             elseif(line(i:j).eq.'lead')then    !leading particle (neads Standard Variable)
               xmatri(ntri(nhis),nhis)=-123456
               istavar=1
             elseif(line(i:j).eq.'jet')then    !jet form fastjet
               xmatri(ntri(nhis),nhis)=nhis*100
             else
               kk=0
               do k=i+1,j-1
                 if(line(k:k).eq.'%')kk=k
               enddo
               if(kk.eq.0)then
                 read(line(i:j),*,err=99003)xmatri(ntri(nhis),nhis)
                 xmatrp(ntri(nhis),nhis)=100.
               else
                 read(line(i:kk-1),*)xmatrp(ntri(nhis),nhis)
                 read(line(kk+1:j),*)xmatri(ntri(nhis),nhis)
               endif
             endif
            else!nsingle=1
              xmatri(ntri(nhis),nhis)=xmitri(ntri(nhis),nhis)*1.00001
              xmitri(ntri(nhis),nhis)=xmitri(ntri(nhis),nhis)*0.99999
            endif
            !---exchange min-max------------------
            if(xmitri(ntri(nhis),nhis).gt.xmatri(ntri(nhis),nhis))then
              xmatri_save=xmatri(ntri(nhis),nhis)
              xmatrp_save=xmatrp(ntri(nhis),nhis)
              xmatri(ntri(nhis),nhis)=xmitri(ntri(nhis),nhis)
              xmatrp(ntri(nhis),nhis)=xmitrp(ntri(nhis),nhis)
              xmitri(ntri(nhis),nhis)=xmatri_save
              xmitrp(ntri(nhis),nhis)=xmatrp_save
            endif
            !-------------------------------------
          else
            ivar(1,nhis)=-1
            call utword(line,i,j,1) !xmin
            call utword(line,i,j,1) !xmax
          endif
        enddo
      elseif(line(i:j).eq.'noerrorbut')then !-----------
        ionoerr=ionoerr+1
        if(ionoerr.gt.2)stop'xini: not more than 2 noerrorbut !   '
        noerrall=noerrall+1
        if(noerrall.gt.mxhis/2)stop'xini: to many noerrorbut     '

        call utword(line,i,j,1) !variable-name
        cvar='                         '
        cvar(1:j-i+1)=line(i:j)
        call xtrans(cvar,icas,inom,ifrnew,nhis)
        if(inom.gt.0)then
          if(inom.gt.100)then
            write(*,*)'xini: noerrorbut can not be used with :',cvar
            stop'xini: error with noerrorbut!'
          endif
          noerrhis(nhis)=noerrall-ionoerr+1
          noerr(noerrhis(nhis),ionoerr)=inom
          do nb=1,nbin(nhis)
             ebin(nb,nac(nhis),ionoerr-1+noerrhis(nhis))=0.d0
                  zebin(nb,nac(nhis),ionoerr-1+noerrhis(nhis))=0.d0
          enddo
        else
          ionoerr=ionoerr-1
          noerrall=noerrall-1
        endif
        inoerr(nhis)=ionoerr
      elseif(line(i:j).eq.'write')then !-----------
        stop'write should not be in bh...eh environment. missing eh?'
      elseif(line(i:j).eq.'writearray'.or.line(i:j).eq.'wa')then !----
        call utword(line,i,j,1)
        iologb=0
        iocnxb=0
      elseif(line(i:j).eq.'writehisto')then !-----------
        call utword(line,i,j,1)
        iologb=0
        iocnxb=0
      elseif(line(i:j).eq.'endhisto'.or.line(i:j).eq.'endanalysis'
     .      .or.line(i:j).eq.'eh')then   !-----------
        ilog(nhis)=.false.
        icnx(nhis)=.false.
        if(iologb.eq.1)ilog(nhis)=.true.
        if(iocnxb.eq.1)icnx(nhis)=.true.
        if(ilog(nhis))then
          xinc(nhis)=1./log(xmax(nhis)/xmin(nhis))*nbin(nhis)
        else
          xinc(nhis)=float(nbin(nhis))/(xmax(nhis)-xmin(nhis))
        endif
        iologb=0
        iocnxb=0
        jjj=j
        cline=line
        goto 9999
      endif
      goto 1

 1006 ierror=1006
      goto 99000
99001 ierror=99001
      goto 99000
99002 ierror=99002
      goto 99000
99003 ierror=99003
      goto 99000
99000 continue
      write(ifmt,'(80a1)')('-',k=1,80)
      write(ifmt,'(a,i7)')'ERROR',ierror
      write(ifmt,'(2a)')'line(i:j): ',line(i:j)
      write(ifmt,'(2a)')'line(1:140): ',line(1:140)
      write(ifmt,'(80a1)')('-',k=1,80)
      stop

 9999 continue
      if(ish.ge.5)then
        do n=1,nhis
          write (ifch,*) n,': ',ivar(1,n),ivar(2,n)
     $         ,'(',ivfra(1,n),ivfra(2,n)
     $         ,')',inorm(n)
     $         ,xmin(n),xmax(n),ilog(n),icnx(n)
     $         ,nbin(n),(idcod(j,n),j=1,nidcod(n))
     $         ,' tri:',ntri(n),(itri(j,n),j=1,ntri(n)),'('
     $         ,(itfra(j,n),j=1,ntri(n)),')'
     $         ,(xmitri(j,n),j=1,ntri(n)) ,(xmatri(j,n),j=1,ntri(n))
        enddo
        write (ifch,*) (ifra(j),j=1,nfra)
      endif
      call utprix('xini  ',ish,ishini,5)
      return
      end


c-----------------------------------------------------------------------
      subroutine actimux(n,ishift)
c-----------------------------------------------------------------------
      include "epos.incxan"
      imux(0)=imux(0)+1
      if(imux(0).gt.mxxhis)stop'mxxhis too small'
      imux(imux(0))=n
      imux(-imux(0))=ishift
      end
c-----------------------------------------------------------------------
      subroutine actiepx(n,ishift)
c-----------------------------------------------------------------------
      include "epos.incxan"
      iepx(0)=iepx(0)+1
      if(iepx(0).gt.mxxhis)stop'mxxhis too small'
      iepx(iepx(0))=n
      iepx(-iepx(0))=ishift
      end

c-----------------------------------------------------------------------
      subroutine actidetphi(n,ishift)
c-----------------------------------------------------------------------
      include "epos.incxan"
      idetphi(0)=idetphi(0)+1
      if(idetphi(0).gt.mxxhis)stop'mxxhis too small'
      iwww=iwww+1
      if(iwww.gt.mxxhis)stop'mxxhis too small'
      idetphi(idetphi(0))=n
      idetphi(-idetphi(0))=ishift
      nhisxxx(n)=iwww
      end

c-----------------------------------------------------------------------
      subroutine actimupurp(n,ishift)
c-----------------------------------------------------------------------
      include "epos.incxan"
      imupurp(0)=imupurp(0)+1
      if(imupurp(0).gt.mxxhis)stop'mxxhis too small'
      iwww=iwww+1
      if(iwww.gt.mxxhis)stop'mxxhis too small'
      imupurp(imupurp(0))=n
      imupurp(-imupurp(0))=ishift
      nhisxxx(n)=iwww
      end


c-----------------------------------------------------------------------
      subroutine actifastjet(n)
c-----------------------------------------------------------------------
      include "epos.incxan"
      iok=0
      do i=1,ifastjet(0)
        if(ifastjet(i).eq.n)iok=1
      enddo
      if(iok.eq.1)return

      ifastjet(0)=ifastjet(0)+1
      if(ifastjet(0).gt.mxxhis)stop'mxxhis too small'
      ifastjet(ifastjet(0))=n
      end              

c-----------------------------------------------------------------------   
      subroutine actijetevent(n)    
c-----------------------------------------------------------------------        
      include "epos.incxan"
      iok=0
      do i=1,ijetevent(0)
        if(ijetevent(i).eq.n)iok=1
      enddo
      if(iok.eq.1)return

      ijetevent(0)=ijetevent(0)+1
      if(ijetevent(0).gt.mxxhis)stop'mxxhis too small'
      ijetevent(ijetevent(0))=n
      end

c---------------------------------------------------------------------
      subroutine xana
c---------------------------------------------------------------------
      include "epos.inc"
      include "epos.incxan"
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall

      parameter (mxfra=5)
      common/pfra/nfra,ifra(mxfra),ivfra(2,mxhis),itfra(mxtri,mxhis)
     $     ,imofra(3,mxfra),iffra(mxfra),r1fra(3,mxfra),r2fra(3,mxfra)
     $     ,emax(mxfra)
      double precision bofra
      common/dfra/bofra(5,mxfra)
      parameter (ntim=1000)
      common/cprt/pprt(5,ntim),q2prt(ntim),idaprt(2,ntim),idprt(ntim)
     &,iorprt(ntim),jorprt(ntim),nprtj

      double precision pgampr,rgampr
      common/cgampr/pgampr(5),rgampr(4)

      common/photrans/phoele(4),ebeam

      dimension ten(4,3)
      logical go,goo(mxcnt),cont,istlo,istlo1,LongLivPtl,NoLongLivParent

      common/cncntje/ncntje

      call utpri('xana  ',ish,ishini,4)

      call xanatest8

      if(ish.ge.2)then
          call alist('fill histograms&',0,0)
      endif

      do n=1,nhis
      do i=1,mypara
        ypara(i,n)=0
      enddo
      enddo

      do n=1,mxxhis
        do m=1,myyarr
          yyarr(n,m)=0
        enddo
      enddo
      
      !if(ish.ge.2)write(ifmt,*)'XANA minfra maxfra :',minfra, maxfra


      if(ish.ge.5)write(ifch,*)'frames ...'

      if(iappl.eq.6)then
        if(mod(iolept/10,10).eq.1) call gakjet(1)
        if(mod(iolept/100,10).eq.1) call gakjet(2)
      endif

      do l=1,nfra
        emax(l)=egyevt/2
        if(ifra(l).eq.12)emax(l)=sqrt(pnll**2+prom**2)
        if(ifra(l).eq.iframe)then
          if(iappl.eq.1.and.iframe.eq.22)emax(l)=ebeam
          imofra(1,l)=0
          imofra(2,l)=0
          imofra(3,l)=0
        elseif(ifra(l).eq.11.or.ifra(l).eq.12)then
          imofra(1,l)=0
          imofra(2,l)=0
          bofra(1,l)=0d0
          bofra(2,l)=0d0
          bofra(3,l)=dsinh(dble(yhaha))
          bofra(4,l)=dcosh(dble(yhaha))
          bofra(5,l)=1d0
          if(ifra(l).eq.11.and.iframe.eq.12)then
            imofra(3,l)=1       ! target -> NN
          elseif(ifra(l).eq.12.and.iframe.eq.11)then
            imofra(3,l)=-1      ! NN -> target
          else
            imofra(3,l)=0       ! not known
          endif
        elseif(ifra(l).eq.21)then
          if(iframe.ne.21)then
            print *, 'invalid frame request'
            print *, 'choose frame gamma-nucleon for event run'
            stop'bye bye'
          endif
        elseif(ifra(l).eq.22)then
          if(iappl.eq.1)emax(l)=ebeam
          if(iframe.eq.21)then
            imofra(1,l)=-1      !'  trafo gN -> lab'
            imofra(1,l)=0
            r1fra(1,l)=rgampr(1)
            r1fra(2,l)=rgampr(2)
            r1fra(3,l)=rgampr(3)
            imofra(2,l)=0
            if(iappl.eq.1)then
              imofra(3,l)=-2
              bofra(1,l)=dsinh(dble(yhaha)) !used for first boost in targ frame
              bofra(2,l)=dcosh(dble(yhaha)) !here : pgampr(1)=pgampr(2)=0.
              bofra(3,l)=pgampr(3)
              bofra(4,l)=pgampr(4)
              bofra(5,l)=pgampr(5)
            else
              imofra(3,l)=-1
              bofra(1,l)=pgampr(1)
              bofra(2,l)=pgampr(2)
              bofra(3,l)=pgampr(3)
              bofra(4,l)=pgampr(4)
              bofra(5,l)=pgampr(5)
            endif
          elseif(iframe.eq.22)then
                                ! nothing to do already gN-frame
          else
            print *, 'invalid frame request'
            print *, 'choose frame gamma-nucleon or lab for event run'
            stop'bye bye'
          endif
        elseif(ifra(l).eq.23)then
          if(iframe.eq.21)then
            imofra(1,l)=0       ! gN -> breit-frame
            r1fra(1,l)=rgampr(1)
            r1fra(2,l)=rgampr(2)
            r1fra(3,l)=rgampr(3)
            imofra(2,l)=0
            imofra(3,l)=1
            bofra(1,l)=0d0
            bofra(2,l)=0d0
            bofra(3,l)=rgampr(4)
            bofra(4,l)=sqrt(rgampr(1)**2+rgampr(2)**2+rgampr(3)**2)
            bofra(5,l)=sqrt( bofra(4,l)**2-rgampr(4)**2)
          elseif(iframe.eq.23)then
                                ! nothing to do already breit-frame
          else
            print *, 'invalid frame request'
            print *, 'choose frame gamma-nucleon or lab for event run'
            stop'bye bye'
          endif
        elseif(ifra(l).eq.33.or.ifra(l).eq.36)then
          if(ifra(l).eq.33)then
            call gakthru(ten,2)
          else
            call gakthru(ten,3)
          endif
          if(ten(4,1).lt.0.)then
            imofra(1,l)=0
            imofra(2,l)=0
            imofra(3,l)=0
          else
            arox=ten(1,1)
            aroy=ten(2,1)
            aroz=ten(3,1)
            brox=ten(1,2)
            broy=ten(2,2)
            broz=ten(3,2)
            call utrota(1,arox,aroy,aroz,brox,broy,broz)
            imofra(1,l)=1
            r1fra(1,l)=arox
            r1fra(2,l)=aroy
            r1fra(3,l)=aroz
            imofra(2,l)=1
            r2fra(1,l)=brox
            r2fra(2,l)=broy
            r2fra(3,l)=broz
            imofra(3,l)=0       !no boost
          endif
          bofra(1,l)=dble(ten(4,1)) !usually this is for boosting
          bofra(2,l)=dble(ten(4,2)) !I abuse it to store the eigenvalues
          bofra(3,l)=dble(ten(4,3)) !
        elseif(ifra(l).eq.32.or.ifra(l).eq.34.or.ifra(l).eq.35)then
          if(ifra(l).eq.32)then
            call gaksphe(ten,2.,2)
          elseif(ifra(l).eq.34)then
            call gaksphe(ten,1.,2)
          else
            call gaksphe(ten,2.,3)
          endif
          if(ten(4,1).lt.0.)then
            imofra(1,l)=0
            imofra(2,l)=0
            imofra(3,l)=0
          else
            arox=ten(1,1)
            aroy=ten(2,1)
            aroz=ten(3,1)
            brox=ten(1,2)
            broy=ten(2,2)
            broz=ten(3,2)
            call utrota(1,arox,aroy,aroz,brox,broy,broz)
            imofra(1,l)=1
            r1fra(1,l)=arox
            r1fra(2,l)=aroy
            r1fra(3,l)=aroz
            imofra(2,l)=1
            r2fra(1,l)=brox
            r2fra(2,l)=broy
            r2fra(3,l)=broz
            imofra(3,l)=0
          endif
          bofra(1,l)=dble(ten(4,1))
          bofra(2,l)=dble(ten(4,2))
          bofra(3,l)=dble(ten(4,3))
        endif
      enddo

      do n=1,nhis
        itrevt(n)=.false.
        if(ivar(1,n).ge.100.and.ivar(1,n).le.199) sval(1,n)=0.
        if(ivar(2,n).ge.100.and.ivar(2,n).le.199) sval(2,n)=0.
        if(ivar(1,n).gt.300.and.ivar(1,n).lt.400)then
          call xval(n,icasv(1,n),ivar(1,n),ivfra(1,n),0,x)!initialization
        endif
        if(ivar(2,n).gt.300.and.ivar(2,n).lt.400)then
          call xval(n,icasv(2,n),ivar(2,n),ivfra(2,n),0,y) !
        endif
        do j=1,ntri(n)
          valtri(j,n)=0.
        enddo
        do j=1,nbin(n)  !copy bins
          bin(j,3-nac(n),n)=bin(j,nac(n),n)
          zcbin(j,3-nac(n),n)=zcbin(j,nac(n),n)
        enddo
        if(lookcontr(n).gt.0)then
          do j=1,nbin(n)
            do loo=lookcontr(n),lookcontrx(n)
                    bbin(j,3-nac(n),loo)=bbin(j,nac(n),loo)
                    zbbin(j,3-nac(n),loo)=zbbin(j,nac(n),loo)
            enddo
          enddo
        endif
        if(inoerr(n).gt.0)then
          do j=1,nbin(n)
            do nn=1,inoerr(n)
              ebin(j,3-nac(n),nn-1+noerrhis(n))=ebin(j,nac(n),
     &                                          nn-1+noerrhis(n))
              zebin(j,3-nac(n),nn-1+noerrhis(n))=zebin(j,nac(n),
     &                                          nn-1+noerrhis(n))
            enddo
          enddo
        endif
      enddo
      if(istavar.eq.1)then
        if(ish.ge.5)write(ifch,*)'Calculate standard variables ...'
        call StandardVariables
      endif
      if(ipairs1.eq.1)then
        if(ish.ge.5)write(ifch,*)'Calculate pair variables ...'
        call PairVariables
      endif
      if(ish.ge.5)write(ifch,*)'Call CorH ...' !bg
      do n=1,iCorH(0)
        call CorH(iCorH(n))
      enddo
      if(ish.ge.5)write(ifch,*)'Call mupurp ...' 
      do n=1,imupurp(0)
        call mupurp(imupurp(n),imupurp(-n))
      enddo
      if(ish.ge.5)write(ifch,*)'Call detphi ...' 
      do n=1,idetphi(0)
        call detphi(idetphi(n),idetphi(-n))
      enddo
      if(ish.ge.5)write(ifch,*)'Call photfill ...' !bg photfill has to be called before corrtrig
      do n=1,iphoton(0)
        call photfill(iphoton(n))
      enddo
      if(ish.ge.5)write(ifch,*)'Call corrtrig ...'
      do n=1,icorrtrig(0)
        call corrtrig(icorrtrig(n))
      enddo
      if(ish.ge.5)write(ifch,*)'Call hardevent ...'
      do n=1,ihardevent(0)
        call hardevent(ihardevent(n))
      enddo
      if(ish.ge.5)write(ifch,*)'Call mux ...'
      do n=1,imux(0)
        call mux(imux(n),imux(-n))
      enddo
      if(ish.ge.5)write(ifch,*)'Call epx ...'
      do n=1,iepx(0)
        call epx(iepx(n),iepx(-n))
      enddo
      if(ish.ge.5)write(ifch,*)'Call caltrig ...'
      do n=1,icaltrig(0)
        call caltrig(icaltrig(n))
      enddo
      if(ish.ge.5)write(ifch,*)'Call jetfind ...'
      do n=1,ijetfind1(0)
        call jetfind(1,ijetfind1(n))
      enddo
      do n=1,ijetfind2(0)
        call jetfind(2,ijetfind2(n))
      enddo
      if(ish.ge.5)write(ifch,*)'Call fastjet ...'
      do n=1,ifastjet(0)
        call fastjet(ifastjet(n))
      enddo
c don't change order here : jetevent should be called always after fastjet !
      if(ish.ge.5)write(ifch,*)'Call jetevent ...'
      ncntje=0
      do n=1,ijetevent(0)
        call jetevent(ijetevent(n))
      enddo

c...........................loop nptl...................................
      ncontr=0
      if(ish.ge.5)write(ifch,*)'Loop nptl ...'

      !call alist('TEST&',1,nptl)
      do iloo=1,max(40,nptl)
      
        j=iloo   

        if(istptl(j).lt.100.and.istptl(j).gt.istmax) goto 8
        if(ish.ge.5)write(ifch,*)'ptl :',j
        do i=1,nfra
          iffra(i)=0            !flag if frame calculated or not
        enddo
        ihad=0
        ibar=0
        jbar=0
        imes=0
        iall4=0
           if(idptl(j).eq.31.or.idptl(j).eq.25)then
             print*,'ERROR 11072016b',iorptl(j),jorptl(j),iloo,j
     .         ,idptl(j),istptl(j),ityptl(j)
           endif
        idepos=ideposf( 1 ,j)
        call idchrg( 10 ,idepos,ch)
        call idflav(idepos,i1,i2,i3,j4,j)
    !~~~~~~~~~~~~~~~~~~
    !  if((abs(i1).eq.4.or.abs(i2).eq.4.or.abs(i3).eq.4)
    ! ..and.(istptl(j).eq.21.or.istptl(j)/2.eq.0))
    ! .print*,'777',j,idepos,istptl(j),'   '
    ! .,idptl(iorptl(j)),istptl(iorptl(j))
    ! .,idptl(iorptl(iorptl(j))),istptl(iorptl(iorptl(j)))
    !~~~~~~~~~~~~~~~~~~
        if(i1.eq.0.and.i2.ne.0.and.i3.ne.0)imes=1
        if(i1.ne.0.and.i2.ne.0.and.i3.ne.0)ibar=1 
        if(i1.gt.0.and.i2.gt.0.and.i3.gt.0)jbar=1 
        if(i1.lt.0.and.i2.lt.0.and.i3.lt.0)jbar=-1 
        if(imes.eq.1.or.ibar.eq.1)ihad=1
        if(abs(idepos).lt.9990)iall4=1 !???????id???? update this !!!!!!

        do n=1,nhis
          if(ivar(1,n).eq.-1.or.ivar(2,n).eq.-1)goto 9
          if(ivar(1,n).ge.200.and.ivar(2,n).ge.200)goto 9 !skip particle loop if event variables
          go=nidcod(n).eq.0

c...........check ist
          ist=istptl(j)
          call setIstXY(0,istuse(n),istxxx,istyyy)
          istlo=ist.eq.istxxx
          istlo1=ist.eq.istxxx.or.ist.eq.istyyy.or.ist.ge.10
          if(noweak(n).eq.1)then  !include long lived weakly decaying ptls 
            if(LongLivPtl(istyyy,j))istlo=.true.
          elseif(noweak(n).eq.3)then  !exclude long lived weakly decaying ptls 
            if(LongLivPtl(istxxx,j))istlo=.false.
          endif
          if(istuse(n).eq.50.or.istuse(n).eq.70)
     .    istlo=ist.eq.istxxx.or.ist.eq.istyyy

c...........check ids
          do i=1,nidcod(n)
            if(ist.eq.0.and.idcod(i,n).eq.10000)then      !all final particle
              go=.true.     
            elseif(istlo.and.idcod(i,n).eq.9995)then      !all particles but nuclei
              if(iall4.eq.1) go=.true.     
            elseif(istlo.and.idcod(i,n).eq.9990)then      !all hadrons
              if(ihad.eq.1) go=.true.      
            elseif(istlo.and.idcod(i,n).eq.9985)then      !neutral particles
              if(abs(ch).lt.0.1.and.iall4.eq.1) go=.true. 
            elseif(istlo.and.idcod(i,n).eq.9980)then      !charged particles
              if(abs(ch).gt.0.1.and.iall4.eq.1) go=.true. 
            elseif(istlo.and.idcod(i,n).eq.9975)then      !neutral hadrons
              if(abs(ch).lt.0.1.and.ihad.eq.1) go=.true.
            elseif(istlo.and.idcod(i,n).eq.9970)then      !charged hadrons
              if(abs(ch).gt.0.1.and.ihad.eq.1) go=.true.
            elseif(istlo.and.idcod(i,n).eq.-9960)then     !negative hadrons
              if(ch.lt.-0.1.and.ihad.eq.1)go=.true.
            elseif(istlo.and.idcod(i,n).eq.9960)then      !positive hadrons
              if(ch.gt.0.1.and.ihad.eq.1)go=.true.
            elseif(istlo.and.idcod(i,n).eq.9940)then      !baryons and antibaryons
              if(ibar.eq.1) go=.true.
            elseif(istlo.and.idcod(i,n).eq.9939)then      !baryons
              if(jbar.eq.1) go=.true.
            elseif(istlo.and.idcod(i,n).eq.-9939)then      !antibaryons
              if(jbar.eq.-1) go=.true.
            elseif(istlo1.and.idcod(i,n).eq.idepos)then
              go=.true.
            elseif(ist.eq.9
     &          .and.idcod(i,n).eq.idepos
     &          .and.(ityptl(j).eq.80.or.ityptl(j).eq.81))then!resonances marked in u.f
              go=.true.
            endif
          enddo
          if(ish.ge.10)write(ifch,*)
c          print *,
     &     n,j,istuse(n),'      ist,id',istptl(j),idepos,'    ',go

c...........check parent
          if(go)then
            if(noweak(n).gt.0)then 
              go=NoLongLivParent(noweak(n),j)
            endif
          endif

c...........check triggers
          if(go)then
            if(ish.ge.7)write(ifch,*)'  check triggers in histogram ',n
            ncontr=0
            do i=1,ntri(n)
              if(ish.ge.7)write(ifch,*)'  trigger variable: ',itri(i,n)
              if(itri(i,n).lt.100)then
                call xval(n,icast(i,n),itri(i,n),itfra(i,n),j,x)
                if(ntrc(i,n).ne.-1)then
                  call triggercondition(i,n,x,go)
                  !if(go)print*,'XAN ptltrg ',i,j,idepos,x,go
                else
                  ncontr=ncontr+1
                  goo(ncontr)=.true.
                  call triggercondition(i,n,x,goo(ncontr))
                  if((ivar(1,n).gt.100.and.ivar(1,n).lt.200)
     .                 .or.(ivar(2,n).gt.100.and.ivar(2,n).lt.200))then
                    print*,'!-----------------------------------------'
                    print*,'!  100-199 event variables can not be used'
                    print*,'! in connection with "trigger contr ..."  '
                    print*,'!-----------------------------------------'
                    stop'in xana (1).                      '
                  endif
                endif
              elseif(itri(i,n).lt.200)then
                if(ntrc(i,n).eq.-1)then
                    print*,'!-----------------------------------------'
                    print*,'!  100-199 event variables can not be used'
                    print*,'! in connection with "trigger contr ..."  '
                    print*,'!-----------------------------------------'
                    stop'in xana (2).                      '
                endif
                call xval(n,icast(i,n),itri(i,n),itfra(i,n),j,x)
                valtri(i,n)=valtri(i,n)+x
              endif
            enddo
          endif
            
c............fill histogram
          if(go)then
            if(ish.ge.6)write(ifch,*)'  fill histogram '
     &            ,n,ivar(1,n),ivar(2,n),ivfra(2,n)
            cont=.false.
            if(ivar(1,n).lt.100.or.ivar(2,n).lt.100)then
              if(ivar(2,n).lt.100)then
                call xval(n,icasv(2,n),ivar(2,n),ivfra(2,n),j,y)
                sval(2,n)=y
                cont=.true.
              endif
              if(ivar(1,n).lt.100)then
                call xval(n,icasv(1,n),ivar(1,n),ivfra(1,n),j,x)
                if(x.ge.xmin(n).and.x.le.xmax(n))then
                  norm3=mod(inorm(n)/100,10)
                  if(norm3.eq.1)then
                    y=y*abs(x)
                  elseif(norm3.eq.2.and.ivar(1,n).eq.63.and.x.ne.0.)then
                    y=y/(x+pptl(5,j))/2/pi
                  elseif(norm3.eq.2.and.ivar(1,n).ne.63.and.x.ne.0.)then
                    y=y/x/2/pi
                  elseif(norm3.eq.4.and.x.ne.0.)then
                    y=y/x**1.5
                  elseif(norm3.eq.5.and.x.ne.0.)then
                    y=y/x
                  elseif(norm3.eq.7.and.x.ne.0.)then
                    y=y/x/sqrt(x-pptl(5,j))
                  endif
                  if(icnx(n))then
                    call fillhistoconex(n,x,y,ivfra(2,n),j)   !for conex
                  else
                    if(ilog(n))then
                      nb=1+int(log(x/xmin(n))*xinc(n))
                    else
                      nb=1+int((x-xmin(n))*xinc(n))
                    endif
                    bin(nb,nac(n),n)=bin(nb,nac(n),n)+y
                   if(ncontr.gt.0)then  !ptl trigger contr
                      do nn=1,ncontr
                        if(goo(nn))then
                          bbin(nb,nac(n),lookcontr(n)-1+nn)=
     &                    bbin(nb,nac(n),lookcontr(n)-1+nn)+y
                          zbbin(nb,nac(n),lookcontr(n)-1+nn)=
     &                    zbbin(nb,nac(n),lookcontr(n)-1+nn)+1
                        endif
                      enddo
                    endif
                    if(inoerr(n).gt.0)then
                    do nn=1,inoerr(n)
                    call xval(n,1,noerr(noerrhis(n),nn),ivfra(2,n),j,y)!?????
                        ebin(nb,nac(n),nn-1+noerrhis(n))=
     &                       ebin(nb,nac(n),nn-1+noerrhis(n))+y
                        zebin(nb,nac(n),nn-1+noerrhis(n))=
     &                       zebin(nb,nac(n),nn-1+noerrhis(n))+1
                    enddo
                    endif
                    zcbin(nb,nac(n),n)=zcbin(nb,nac(n),n)+1
                  endif
                  itrevt(n)=.true.
                endif
              endif
            endif
            if(ivar(1,n).gt.100.and.ivar(1,n).lt.200)then
              call xval(n,icasv(1,n),ivar(1,n),ivfra(1,n),j,x)
              sval(1,n)=sval(1,n)+x
              !print*,'XAN finaladd x',n,j,idptl(j),istptl(j),sval(1,n)
            endif
            if(ivar(2,n).gt.100.and.ivar(2,n).lt.200)then
              call xval(n,icasv(2,n),ivar(2,n),ivfra(2,n),j,y)
              sval(2,n)=sval(2,n)+y
              cont=.true.
              !print*,'XAN finaladd ',y,idepos
            endif
            if(ivar(1,n).gt.300.and.ivar(1,n).lt.400)then
              call xval(n,icasv(1,n),ivar(1,n),ivfra(1,n),j,x)
            endif
            if(ivar(2,n).gt.300.and.ivar(2,n).lt.400)then
              call xval(n,icasv(2,n),ivar(2,n),ivfra(2,n),j,y)
              cont=.true.
            endif
            if(ish.ge.6.and.cont)write (ifch,*)
     *                           '   ---> histo n,x,y:',n,x,y
          endif
   9      continue
        enddo
  8     continue
      
      enddo
c...........................end loop nptl...........................

      do n=1,nhis
      if(ivar(1,n).eq.-1.or.ivar(2,n).eq.-1)goto 99

c........check event triggers

      if(ish.ge.7)write(ifch,*)'  check event triggers in histogram ',n
       go=.true.
        ncontr=0
        do i=1,ntri(n)
          if(itri(i,n).gt.100)then
            if(itri(i,n).lt.200)then
              x=valtri(i,n)
            else
              call xval(n,icast(i,n),itri(i,n),itfra(i,n),0,x)
            endif
            if(ntrc(i,n).ne.-1)then
              call triggercondition(i,n,x,go)
            else
              ncontr=ncontr+1
              goo(ncontr)=.true.
              call triggercondition(i,n,x,goo(ncontr))
            endif
            !print*,'XAN evttrg    HIS TR',n,i,go,x,icast(i,n)
          endif
        enddo

c........event variables > 200

        if(go)then
          if(ivar(1,n).gt.100)then
            if(ivar(1,n).gt.200.and.ivar(1,n).lt.300)then
              call xval(n,icasv(1,n),ivar(1,n),ivfra(1,n),0,x)
            elseif(ivar(1,n).gt.300.and.ivar(1,n).lt.400)then
              call checkcccptl(-99)
              call xval(n,icasv(1,n),ivar(1,n),ivfra(1,n),-99,x)
            elseif(ivar(1,n).gt.100.and.ivar(1,n).lt.200)then
              x=sval(1,n)
            else
              call xval(n,icasv(1,n),ivar(1,n),ivfra(1,n),0,x)
            endif
            if(ivar(2,n).gt.200.and.ivar(2,n).lt.300)then
              call xval(n,icasv(2,n),ivar(2,n),ivfra(2,n),0,y)
            elseif(ivar(2,n).gt.300.and.ivar(2,n).lt.400)then
              call checkcccptl(-99)
              call xval(n,icasv(2,n),ivar(2,n),ivfra(2,n),-99,y)
            elseif(ivar(2,n).gt.0.and.ivar(2,n).lt.200)then
              y=sval(2,n)
            else             !inom>500
              call xval(n,icasv(2,n),ivar(2,n),ivfra(2,n),0,y)
            endif
c The following doesn't work for ivar(2,n)<100, since particle number is not defined !
c            if(ivar(2,n).gt.200.and.ivar(2,n).lt.300)then
c              call xval(n,icasv(2,n),ivar(2,n),ivfra(2,n),0,y)
c            elseif(ivar(2,n).gt.300.and.ivar(2,n).lt.400)then
c              call xval(n,icasv(2,n),ivar(2,n),ivfra(2,n),-99,y)
c            elseif(ivar(2,n).gt.100.and.ivar(2,n).lt.200)then
c              y=sval(2,n)
c            else
c              call xval(n,icasv(2,n),ivar(2,n),ivfra(2,n),0,y)
c            endif
            if(mod(inorm(n)/100,10).eq.1)y=y*x
            if(mod(inorm(n)/100,10).eq.2.and.x.ne.0.)y=y/x/2/pi
            if(mod(inorm(n)/100,10).eq.4.and.x.ne.0.)y=y/x**1.5
            if(mod(inorm(n)/100,10).eq.5.and.x.ne.0.)y=y/x
            sval(1,n)=x
            sval(2,n)=y
            if(ish.ge.6)write (ifch,*) 'histo n,x,y:',n,x,y
            if(x.ge.xmin(n).and.x.le.xmax(n))then
              if(ilog(n))then
                nb=1+int(log(x/xmin(n))*xinc(n))
              else
                nb=1+int((x-xmin(n))*xinc(n))
              endif
              bin(nb,nac(n),n)=bin(nb,nac(n),n)+y
              if(ncontr.gt.0)then
                do nn=1,ncontr
                  if(goo(nn))then
                    bbin(nb,nac(n),lookcontr(n)-1+nn)=
     &              bbin(nb,nac(n),lookcontr(n)-1+nn)+y
                    zbbin(nb,nac(n),lookcontr(n)-1+nn)=
     &              zbbin(nb,nac(n),lookcontr(n)-1+nn)+1
                  endif
                enddo
              endif
              zcbin(nb,nac(n),n)=zcbin(nb,nac(n),n)+1
              itrevt(n)=.true.
            endif
            !if(y.gt.0.)print*,'XAN finalxy :   ',x,y
          endif
        endif

c........particle variables

        if(go)then
          if(ivar(1,n).le.100)then
            if(ncontr.gt.0)then  !event trigger contr
              do nb=1,nbin(n)
                do nn=1,ncontr
                  if(goo(nn))
     &                     bbin(nb,nac(n),lookcontr(n)-1+nn)=
     &              bbin(nb,nac(n),lookcontr(n)-1+nn)
     &              +bin(nb,nac(n),n)-bin(nb,3-nac(n),n)
                enddo
              enddo
            endif
          endif
        endif
 
c............event ok (increase ncevt) or not (take copy)

        if(go)then
          ncevt(n)=ncevt(n)+1
          if(ncontr.gt.0)then
            do nn=1,ncontr
              loo=lookcontr(n)-1+nn
              if(goo(nn))
     &        nccevt(loo)=nccevt(loo)+1
            enddo
          endif
        else
          if(ish.ge.6)write (ifch,*) 'event rejected for histo',n
          nac(n)=3-nac(n)
          itrevt(n)=.false.
        endif

 99   continue
      enddo

      call utprix('xana  ',ish,ishini,4)
      end

c--------------------------------------------------------------------
      subroutine triggercondition(i,n,x,go)
c--------------------------------------------------------------------
c ntrc is used to distinguish the different usage of trigger:
c
c    trigger var xmin xmax
c             ntrc=1
c    trigger or n var1 xmin1 xmax1 var2 xmin2 xmax2 ... varn xminn xmaxn
c          1  ntrc=2
c          2  ntrc=0
c              ...
c         n-1 ntrc=0
c          n  ntrc=3
c    trigger contr n var1 xmin1 xmax1 var2 xmin2 xmax2 ... varn xminn xmaxn
c             ntrc=-1
c--------------------------------------------------------------------
      include "epos.inc"
      include "epos.incxan"
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall
      logical go,gox,ok,goz
                xmn=xmitri(i,n)
                xmx=xmatri(i,n)
                if(xmn.eq.-123456.and.xmx.eq.-123456)then  !for leading part
                  xmn=float(idlead)
                  xmx=float(idlead)
                endif
                pmn=xmitrp(i,n)
                pmx=xmatrp(i,n)
                if(abs(ntrc(i,n)).eq.1)then
                  goz=.true.
                  if(pmn.gt.99.999.and.pmx.gt.99.999)then
                    if(x.lt.xmn.or.x.gt.xmx)goz=.false.
                  else
                    if(x.lt.xmn-0.5.or.x.gt.xmx+0.5)goz=.false.
                    ok=rangen().le.xmitrp(i,n)/100.
                    if(.not.ok.and.x.lt.xmn+0.5)goz=.false.
                    ok=rangen().le.xmatrp(i,n)/100.
                    if(.not.ok.and.x.gt.xmx-0.5)goz=.false.
                  endif
                  if(.not.goz)go=.false.
                else
                  if(ntrc(i,n).eq.2)gox=.false.
                  goz=.true.
                  if(pmn.gt.99.999.and.pmx.gt.99.999)then
                    if(x.lt.xmn.or.x.gt.xmx)goz=.false.
                  else
                    if(x.lt.xmn-0.5.or.x.gt.xmx+0.5)goz=.false.
                    ok=rangen().le.xmitrp(i,n)/100.
                    if(.not.ok.and.x.lt.xmn+0.5)goz=.false.
                    ok=rangen().le.xmatrp(i,n)/100.
                    if(.not.ok.and.x.gt.xmx-0.5)goz=.false.
                  endif
                  if(goz)gox=.true.
                  if(ntrc(i,n).eq.3.and..not.gox)go=.false.
                endif
                if(ish.ge.9)write(ifch,*)'trigger conditions '
     &                                   ,i,n,xmn,x,xmx,go
                end

c-----------------------------------------------------------------------
      subroutine fillhistoconex(n,x,y,lf,iloo)   !for conex
c-----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incxan"
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall

        j=iloo   

           if(.not.(mod(inorm(n),10).ne.4
     &      .and.mod(inorm(n),10).ne.6
     &      .and.mod(inorm(n)/100,10).ne.3))return
                    if(ilog(n))then
                      c=(xmax(n)/xmin(n))**(1./real(nbin(n)))
                      nde=nint(1./log10(c))
                      nb=max(1,1+int(log10(x/xmin(n))*nde))
                      xmb=xmin(n)*c**(nb-0.5)
                      if(x.gt.xmb.and.nb.lt.nbin(n))then
                        if(x.gt.xmax(n))
     &                    write(ifmt,*)'xana max ?',x,xmax(n),nb
                        nbx=1
                        xmx=c*xmb
                      elseif(x.lt.xmb.and.nb.gt.1)then
                        if(x.lt.xmin(n))write(ifmt,*)'xana min ?',x,nb
                        nbx=-1
                        xmx=xmb/c
                      else
                        nbx=0
                        xmx=0.
                      endif
                    else
                      c=(xmax(n)-xmin(n))/real(nbin(n))
                      nb=max(1,1+int((x-xmin(n))/c))
                      xmb=xmin(n)+c*(nb-0.5)
                      if(x.gt.xmb)then
                        nbx=1
                        xmx=c+xmb
                      elseif(x.lt.xmb)then
                        nbx=-1
                        xmx=xmb-c
                      else
                        nbx=0
                        xmx=0.
                      endif
                    endif
                    xc=(x-xmx)/(xmb-xmx)
                    xc=max(0.,min(1.,xc))
                    bin(nb,nac(n),n)=bin(nb,nac(n),n)+xc*y
                    if(nbx.ne.0)bin(nb+nbx,nac(n),n)
     &                  =bin(nb+nbx,nac(n),n)+(1.-xc)*y
                    zcbin(nb,nac(n),n)=zcbin(nb,nac(n),n)+1
                    if(inoerr(n).gt.0)then
                      do nn=1,inoerr(n)
                        call xval(n,1,noerr(noerrhis(n),nn),lf,j,y2)!???????
                        ebin(nb,nac(n),nn-1+noerrhis(n))=
     &                       ebin(nb,nac(n),nn-1+noerrhis(n))+y2
                        zebin(nb,nac(n),nn-1+noerrhis(n))=
     &                       zebin(nb,nac(n),nn-1+noerrhis(n))+1
                      enddo
                    endif
      end

c---------------------------------------------------------------------
      subroutine xhis(n)
c---------------------------------------------------------------------
      include "epos.inc"
      include "epos.incxan"
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall
      dimension xx(mxbin)

      double precision histoweight
      common/chiswei/histoweight
      common/cyield/yield
      common/csigma/sigma
      common/geom/rmproj,rmtarg,bmax,bkmx
      save cnormx

      if(ivar(1,n).eq.-1)then
        nrbins=0
        goto 9999
      endif

c.......here normalization.......................................
c           see also   "..........fill histogram"
c.................................................................
c     the norm ( inorm(n) ) is a number hijk which normalizes to:
c
c  k  0:  * 1
c     1:  / number of events
c     2:  / number of triggered events
c     4:  / bin-counts
c     5:  / bin sum
c     6:  / number of summed bin-counts (yield=1.)
c     7:  uses same normalization as one histo before
c
c  j  0:  * 1
c     1:  / bin-width
c     2:  * sigma_total / bin-width
c     3:  * sigma_diff / bin-width
c
c  i  0:  * 1
c     1:  y => y*|x|
c     2:  y => y/x/2/pi (modified for mt0)
c     3:  kno-scaling
c     4:  y => y/x**1.5
c     5:  y => y/x
c     6:  y => y*xi (for conex, xi=x of the bin)
c     7:  y => y/x/(x-m)
c
c  h  0: normal
c     1: accumulated
c
c.................................................................

      norm1=mod(inorm(n),10)
      norm2=mod(inorm(n)/10,10)
      norm3=mod(inorm(n)/100,10)
      norm4=mod(inorm(n)/1000,10)
      nctbin=0
      sumbin=0
      do l=1,nbin(n)
        nctbin=nctbin+zcbin(l,nac(n),n)
        sumbin=sumbin+bin(l,nac(n),n)
        if(norm1.eq.4.and.zcbin(l,nac(n),n).ne.0d0)then
          bin(l,nac(n),n)=bin(l,nac(n),n)/zcbin(l,nac(n),n)
          if(lookcontr(n).gt.0)then
            do loo=lookcontr(n),lookcontrx(n)
              if(zbbin(l,nac(n),loo).ne.0.)
     &           bbin(l,nac(n),loo)=bbin(l,nac(n),loo)
     &            /zbbin(l,nac(n),loo)
            enddo
          endif
        endif
        if(ilog(n))then
          xx(l)=xmin(n)*(xmax(n)/xmin(n))**((float(l)-.5)/nbin(n))
        else
          xx(l)=(float(l)-0.5)*(xmax(n)-xmin(n))/nbin(n)+xmin(n)
        endif
      enddo
      cnorm=1.
      if(norm1.eq.1.or.norm1.eq.2)cnorm=0 
      if(norm1.eq.5.or.norm1.eq.6)cnorm=0 
      if(norm1.eq.1)then
        if(nevent.ne.0)cnorm=1./float(nevent)
      endif
      if(norm1.eq.2)then
        if(ncevt(n).ne.0)cnorm=1./float(ncevt(n))
      endif
      if(norm1.eq.5.and.sumbin.ne.0.)cnorm=1./sumbin
      if(norm1.eq.6.and.nctbin.ne.0)cnorm=1./float(nctbin)
      if(norm1.eq.7)cnorm=cnormx
      cnormx=cnorm
      sigma=1.
      if(ntevt.ne.0)
     &   sigma=10.*pi*bmax**2.*nevent/ntevt !total (untriggered) sigma
      if(norm2.eq.3)then      !differential (triggered) sigma
        if(ntevt.ne.0)
     &     sigma=10.*pi*bmax**2.*ncevt(n)/ntevt
      endif
c      print *,'---->',sigma,ncevt(n),nevent,ntevt
      if(norm3.eq.3)then      !kno
        first=0.
        secnd=0.
        do l=1,nbin(n)
          if(nctbin.ne.0)first=first+xx(l)*zcbin(l,nac(n),n)/nctbin
          if(nctbin.ne.0)secnd=secnd
     $           +xx(l)**2*zcbin(l,nac(n),n)/nctbin
        enddo
      else
        first=1.
      endif
      if(ilog(n))then
        if(norm2.eq.2.or.norm2.eq.3) cnorm=cnorm*sigma
      else
        if(norm2.ge.1.and.norm2.le.3) cnorm=cnorm*xinc(n)
        if(norm2.eq.2.or.norm2.eq.3) cnorm=cnorm*sigma
      endif
      do l=1,nbin(n)
        bnorm=0.
        if(ilog(n).and.norm2.ge.1.and.norm2.le.3)then
          bnorm=1./(xmin(n)*exp(float(l)/xinc(n))*(1.-exp(-1./xinc(n))))
          bin(l,nac(n),n) =  bin(l,nac(n),n) * bnorm
        endif
        bin(l,nac(n),n) =  bin(l,nac(n),n) * cnorm
        if(lookcontr(n).gt.0)then
          if(ilog(n).and.norm2.ge.1.and.norm2.le.3)then
            do loo=lookcontr(n),lookcontrx(n)
              bbin(l,nac(n),loo)=bbin(l,nac(n),loo) * bnorm
            enddo
          endif
        endif
      enddo
      f=first
      nrbins=nbin(n)
      nctbin=0
      yield=0.
      shft=0
       if(nint(xpara(1,n)).eq.999963)shft=xpara(2,n)
      do ii=1,nbin(n)
        g=1
        if(norm3.eq.1.and.xx(ii).ne.0.)g=1./abs(xx(ii))
        if(norm3.eq.2)g=2*pi*(xx(ii)+shft)
        if(norm3.eq.4)g=xx(ii)**1.5
        if(norm3.eq.5)g=xx(ii)
        if(norm3.eq.7)g=0
        hi=1./xinc(n)
        if(norm2.eq.0)hi=1.
        yield=yield+bin(ii,nac(n),n)*hi*hisfac*f*g
      enddo
      do l=1,nbin(n)
        x=(xx(l)+xshift)      !*xhfact
        ar(l,1)=x/f
        sigbin=0
        if(zcbin(l,nac(n),n).ne.0d0)
     *   sigbin=bin(l,nac(n),n)*hisfac*f/sqrt(zcbin(l,nac(n),n))
        if(norm4.eq.0.or.l.eq.1)then
          ar(l,3)=bin(l,nac(n),n)*hisfac*f
          if(lookcontr(n).gt.0)then
           do loo=lookcontr(n),lookcontrx(n)
             r=1
             if(norm1.eq.2.and.nccevt(loo).ne.0.)
     *          r=float(ncevt(n))/nccevt(loo)
             lo=loo-lookcontr(n)+1
             ary(l,lo)=bbin(l,nac(n),loo)*hisfac*f*cnorm*r
             if(zbbin(l,nac(n),loo).gt.0.)then
               ardy(l,lo)=ary(l,lo)/sqrt(zbbin(l,nac(n),loo))
             else
               ardy(l,lo)=0
             endif
             if(norm1.eq.4)ardy(l,lo)=zbbin(l,nac(n),loo)
           enddo
          endif
          if(norm3.eq.6)then   !conex
           ar(l,3)=ar(l,3)*xx(l)
          endif
          ar(l,4)=sigbin
        else
          ar(l,3)=ar(l-1,3)+bin(l,nac(n),n)*hisfac*f
          ar(l,4)=sqrt(ar(l-1,4)**2+sigbin**2)
        endif
        if(inoerr(n).ge.1)then
          if(zebin(l,nac(n),noerrhis(n)).gt.0.d0)then
         ar(l,4)=ebin(l,nac(n),noerrhis(n))/zebin(l,nac(n),noerrhis(n))
          else
            ar(l,4)=0.
          endif
        endif
        if(inoerr(n).eq.2)then
          if(zebin(l,nac(n),noerrhis(n)+1).gt.0.d0)then
            ar(l,5)=ebin(l,nac(n),noerrhis(n)+1)
     .      /zebin(l,nac(n),noerrhis(n)+1)
          else
            ar(l,5)=0.
          endif
        endif
        if(norm1.eq.4)ar(l,4)=zcbin(l,nac(n),n)
      enddo
      if(lookcontr(n).gt.0)then
        do loo=lookcontr(n),lookcontrx(n)
          lo=loo-lookcontr(n)+1
          histowy(lo)=dble( nccevt(loo) )
          if(norm1.eq.1)histowy(lo)=dble(nevent)
          if(norm1.eq.4)histowy(lo)=0d0
        enddo
      endif
      ionoerr=inoerr(n)
      histoweight=dble(ncevt(n))
      if(norm1.eq.1)histoweight=dble(nevent)
      if(norm1.eq.4)histoweight=0d0
 
 9999 hisfac=1.
      xshift=0
      end

c-----------------------------------------------------------------------
      function underly(nowe,isu,iii)  ! noweak(n) istuse(n)
c-----------------------------------------------------------------------
c Underlying event analysis from ATLAS based on primary charged particles
c (all charged particles with ctau>1cm)
c-----------------------------------------------------------------------
      include "epos.inc"
        logical go,LongLivPtl,NoLongLivParent
        save jj,ptjj

        underly=0
        if(iii.eq.0)then
          jj=0
          ptjj=0.
        endif

        do iloo=1,nptl

        npts=iloo
        i=npts

        go=.false.
        call setIstXY(0,isu,istxxx,istyyy)
        if(istptl(i).eq.istxxx)go=.true.
        if(nowe.eq.1)then
          if(LongLivPtl(istyyy,i))go=.true.
        elseif(nowe.eq.3)then
          if(LongLivPtl(istxxx,i))go=.false.
        endif
        if(go)go=NoLongLivParent(nowe,i)
        if(.not.go)goto 60
        call getidptl(i,id)
        call idchrg( 12 ,id,ch)
        if(abs(ch).lt.0.1)goto 60

        pt=sqrt(pptl(2,npts)**2+pptl(1,npts)**2)

        if(iii.eq.0)then ! iii = 0

          if(jj.eq.0)then
            jj=npts
            ptjj=pt 
          else
            if(pt.gt.ptjj)then
              jj=npts
              ptjj=pt           
            endif
          endif       
          underly=ptjj
c          print*,npts,pt,'    ',jj,ptjj,istxxx

        else ! iii > 0  
 
c          if(npts.eq.jj)goto 60
          pz=pptl(3,npts)
          ppp=sqrt(pz**2+pt**2)
          if(ppp.gt.abs(pz))then
            yyy=.5*log((ppp+pz)/(ppp-pz))
          else
            yyy=sign(100.,pz)
          endif
          if( abs(yyy) .lt. 2.5 .and. pt .gt. 0.5 )then
            px=pptl(1,jj)
            py=pptl(2,jj)
            phi=polar(px,py)
            pxx=pptl(1,npts)
            pyy=pptl(2,npts)
            phii=polar(pxx,pyy)
            delphi=abs(phi-phii)
            if(delphi.gt.pi) delphi=2.*pi-delphi
            !print*,'delphi ',delphi  
            if(delphi.lt.pi/3)then ! ~~~~~~~~~~towards~~~~~~~~~~
              if(iii.eq.2)then
                underly=underly + 1. /(2*pi/3) /5.0
              endif 
              if(iii.eq.5)then
                underly=underly+pt /(2*pi/3) /5.0
              endif 
            elseif(delphi.lt.2*pi/3)then !~~~~~~~~~trans~~~~~~~~~~~~
              if(iii.eq.1)then
                underly=underly + 1. /(2*pi/3) /5.0
              endif 
              if(iii.eq.4)then
                underly=underly+pt /(2*pi/3) /5.0
              endif 
            else !~~~~~~~~~~~~~~~~~~away~~~~~~~~~~~~~~~~~~
              if(iii.eq.3)then
                underly=underly + 1. /(2*pi/3) /5.0
              endif 
              if(iii.eq.6)then
                underly=underly+pt /(2*pi/3) /5.0
              endif 
            endif
          endif

        endif ! iii 

  60    continue
        enddo
        
        !print*,'---------UE------->',iii,underly 
        end

c-----------------------------------------------------------------------
      integer function nsdiff(isu,insdif,nowe)
c-----------------------------------------------------------------------
c  returns  1 if trigger condition for NSD fulfilled and 0 otherwise
c  for  UA1 (insdif=0) or UA5 (insdif=1) or CDF (insdif=2) or STAR (insdif=3,4)
C  or BRAHMS (insdif=5) or NA61 (insdif=6) or CMS (insdif=7)
c  or ATLAS (insdif=8, 17) or ALICE (insdif=9, 10, 11, 16) 
c  or CMS hadron level (insdif=12) or CMS hadron level double sided (insdif=13)
c  or LHCf hadron level (insdif=15)
c  nowe ... noweak(histogram number) 
c-----------------------------------------------------------------------
      include "epos.inc"
      integer ncevt,nsdi(0:20)
      logical cont,go,LongLivPtl,NoLongLivParent
      data ncevt/1/
      save nsdi,ncevt
c initialization for each event
      if(ncevt.eq.nrevt)then
        ncevt=ncevt+1
        do i=0,20
          nsdi(i)=-1
        enddo
      endif
      nsdiff=0

      if(insdif.ge.0)then
ckw         if(nsdi(insdif).lt.0)then
      if(nsdi(insdif).lt.10)then  !ckw  need to always do it
                                   !since there are several isu values possible
                                    !can be done more efficiently of course
        iii1=0
        iii2=0
        iii3=0
        ipos=0
        ineg=0

        do iloo=1,nptl

        npts=iloo
        i=npts

        go=.false.
        call setIstXY(0,isu,istxxx,istyyy)
        if(istptl(i).eq.istxxx)go=.true.
        if(nowe.eq.1)then
          if(LongLivPtl(istyyy,i))go=.true.
        elseif(nowe.eq.3)then
          if(LongLivPtl(istxxx,i))go=.false.
        endif
        if(go)go=NoLongLivParent(nowe,i)
        if(.not.go)goto 60
        call getidptl(i,id)
        call idchrg( 12 ,id,ch)
        if(insdif.ne.7.and.insdif.ne.14.and.insdif.ne.15)then
          if(abs(ch).lt.0.1)goto 60
        endif

        cont=   idptl(npts).ne.120 .and.idptl(npts).ne.-120
     *   .and.idptl(npts).ne.130 .and.idptl(npts).ne.-130
     *   .and.idptl(npts).ne.1120.and.idptl(npts).ne.-1120
     *   .and.idptl(npts).ne.1130.and.idptl(npts).ne.-1130
     *   .and.idptl(npts).ne.2230.and.idptl(npts).ne.-2230
     *   .and.idptl(npts).ne.2330.and.idptl(npts).ne.-2330
     *   .and.idptl(npts).ne.3331.and.idptl(npts).ne.-3331

        pz=pptl(3,npts)
        pt=sqrt(pptl(2,npts)**2+pptl(1,npts)**2)
        ppp=sqrt(pz**2+pt**2)
        Etot=pptl(4,npts)
        if(ppp.gt.abs(pz))then
          yyy=.5*log((ppp+pz)/(ppp-pz))
        else
          yyy=sign(100.,pz)
        endif
        if(insdif.eq.0)then
          if(yyy.gt.1.5  .and. yyy.lt.5.5)iii1=1
          if(yyy.gt.-5.5 .and. yyy.lt.-1.5)iii2=1
        elseif(insdif.eq.1)then
          if(yyy.gt.2.   .and. yyy.lt.5.6)iii1=1
          if(yyy.gt.-5.6 .and. yyy.lt.-2.)iii2=1
        elseif(insdif.eq.2)then
          if(yyy.gt.3.2  .and. yyy.lt.5.9)iii1=1
          if(yyy.gt.-5.9 .and. yyy.lt.-3.2)iii2=1
          if(yyy.gt.0.   .and. yyy.lt.3.0)ipos=ipos+1
          if(yyy.gt.-3.0 .and. yyy.lt.0. )ineg=ineg+1
        elseif(insdif.eq.3)then
          if(yyy.gt.-5.0 .and. yyy.lt.-3.3 )iii1=1
          if(yyy.gt. 3.3 .and. yyy.lt. 5.0 )iii2=1
        elseif(insdif.eq.4)then
          if(yyy.gt.-5.0 .and. yyy.lt.-3.1 )iii1=1
          if(yyy.gt. 3.1 .and. yyy.lt. 5.0 )iii2=1
        elseif(insdif.eq.5)then
          if(yyy.gt.-5.25 .and. yyy.lt.-3.26 )iii1=1
          if(yyy.gt. 3.26 .and. yyy.lt. 5.25 )iii2=1
        elseif(insdif.eq.6)then  !NA61 trigger if NO charged particle with theta<5.26 mrad
          if(pptl(3,npts).gt.0..and.yyy.lt.100.)then
           theta=sqrt(pptl(1,npts)**2+pptl(2,npts)**2)/pptl(3,npts)
           if(theta.lt.5.26e-3)iii1=1
          endif
        elseif(insdif.eq.7)then   !CMS NSD corrected using PYTHIA (2010)
          if(yyy.gt.-5.2 .and. yyy.lt.-2.9 .and. Etot .gt.3.)iii1=1
          if(yyy.gt. 2.9 .and. yyy.lt. 5.2 .and. Etot .gt.3.)iii2=1
          if(yyy.gt. -2.5 .and. yyy.lt. 2.5 .and. pt .gt.0.2 .and. cont)
     &                                                       iii3=1
        elseif(insdif.eq.8)then   !ATLAS
          if(yyy.gt.-2.5 .and. yyy.lt.2.5.and.pt.gt.0.5)iii1=1
          iii2=1
        elseif(insdif.eq.9)then   !ALICE 900 GeV
          yyyy=yyy+rapcms !LAB
          if(yyyy.gt.-3.7 .and. yyyy.lt.-1.7 )iii1=1
          if(yyyy.gt. 2.8 .and. yyyy.lt. 5.1 )iii2=1
        elseif(insdif.eq.10)then   !ALICE 2.36 TeV
          if(yyy.gt.-2 .and. yyy.lt.2 )iii1=1
          iii2=1
        elseif(insdif.eq.11)then   !ALICE Inel>0
          if(yyy.gt.-1 .and. yyy.lt.1 )iii1=1
          iii2=1
        elseif(insdif.eq.12)then   !CMS hadron level NSD trigger (2011)
          if(yyy.gt.-4.4 .and. yyy.lt.-3.9 )iii1=1
          if(yyy.gt. 3.9 .and. yyy.lt. 4.4 )iii2=1
        elseif(insdif.eq.13)then   !CMS hadron level doubl sided trigger (2012)
          if(yyy.gt.-5. .and. yyy.lt.-3. .and. Etot .gt. 3. )iii1=1
          if(yyy.gt. 3. .and. yyy.lt. 5. .and. Etot .gt. 3. )iii2=1
        elseif(insdif.eq.14)then   !CMS hadron level single sided trigger (HF 2012)
          if(yyy.gt.-4.9 .and. yyy.lt.-2.9 .and. Etot .gt. 5. )iii1=1
          if(yyy.gt. 2.9 .and. yyy.lt. 4.9 .and. Etot .gt. 5. )iii2=1
        elseif(insdif.eq.15)then   !LHCf hadron level single sided trigger
          emint=0.
          if(idptl(npts).eq.10.or.abs(idptl(npts)).eq.20)emint=100.
          if(abs(idptl(npts)).eq.1220.or.(abs(idptl(npts)).eq.2130))
     *                                                   emint=300.
          if(yyy.lt.-7.935 .and. Etot .gt. emint )iii1=1
c          if(yyy.gt. 2.9 .and. yyy.lt. 4.9 .and. Etot .gt. 5. )iii2=1
        elseif(insdif.eq.16)then !ALICE inel 
          if(yyy.gt.-3.7 .and. yyy.lt.-1.7 )iii1=1
          if(yyy.gt. 2.8 .and. yyy.lt. 5.1 )iii2=1
          if(yyy.gt. -2  .and. yyy.lt.  2  )iii3=1
        elseif(insdif.eq.17)then !ATLAS 
          if(yyy.gt.-2.5 .and. yyy.lt.2.5.and.pt.gt.0.1)iii1=iii1+1
          iii2=1
        elseif(insdif.eq.18)then !ATLAS 
          if(yyy.gt.-2.5 .and. yyy.lt.2.5.and.pt.gt.0.5)iii1=iii1+1
          iii2=1
        elseif(insdif.eq.19)then !ATLAS 
          if(yyy.gt.-0.8 .and. yyy.lt.0.8.and.pt.gt.0.5)iii1=iii1+1
          iii2=1
        endif
60      continue

        enddo

        if(insdif.le.1)then                                            ! 1
          if(iii1.eq.1 .and. iii2.eq.1) nsdiff=1
        elseif(insdif.eq.2)then                                        ! 2
          if((iii1.eq.1 .and. iii2.eq.1) .and.
     *    ((ipos.ne.0 .and. ineg.ne.0) .and. ipos+ineg.ge.4)) nsdiff=1
        elseif(insdif.eq.3 .or.insdif.eq.4                              ! 3 4 5   8 9 10 11 12 13
     *     .or.insdif.eq.5 .or.insdif.eq.8.or.insdif.eq.9
     *     .or.insdif.eq.10.or.insdif.eq.11.or.insdif.eq.12
     *     .or.insdif.eq.13)then
          if(iii1.eq.1 .and. iii2.eq.1) nsdiff=1
        elseif(insdif.eq.6)then                                        ! 6
          if(iii1.eq.0 .and. iii2.eq.0)then
            nsdiff=1
          endif
        elseif(insdif.eq.7)then                                        ! 7  
          if(iii1.eq.1 .and. iii2.eq.1 .and.iii3.eq.1)then
            nsdiff=1
          endif
        elseif(insdif.eq.14.or.insdif.eq.15)then                       ! 14 15
          if(iii1.eq.1 .or. iii2.eq.1) nsdiff=1
        elseif(insdif.eq.16)then                                       ! 16
          if(iii1.eq.1 .or. iii2.eq.1 .or.iii3.eq.1)then
            nsdiff=1
          endif
        elseif(insdif.eq.17)then                                       ! 17
          if(iii1.ge.2)then
            nsdiff=1
          endif
        elseif(insdif.eq.18)then                                       ! 18
          if(iii1.ge.1)then
            nsdiff=1
          endif
        elseif(insdif.eq.19)then                                       ! 19
          if(iii1.ge.1)then
            nsdiff=1
          endif
        endif

        nsdi(insdif)=nsdiff

      else
        nsdiff=nsdi(insdif)
      endif

      else
        stop'in nsdiff. argument of nsdiff not authorized.        '
      endif

      end

c-----------------------------------------------------------------------
      integer function isdiff(isdif)
c-----------------------------------------------------------------------
c  returns  1 if trigger condition for single diff fulfilled and 0 otherwise
c  for  UA4 Mult distri (isdif=1) or UA4 xsection (isdif=2) 
c  or CDF SD (isdif=3) or CDF DPE (isdif=4) or CDF min bias (for DD) 
c  (isdif=5)
c-----------------------------------------------------------------------
      include "epos.inc"
      isdiff=0
           if(isdif.ge.1)then
      iii0=0
      iii1=0
      iii2=0
      iii3=0
      iii4=0
      Et1=0.
      Et2=0.
      istxxx=0
      if(istfor.gt.-999)istxxx=istfor
      do npts=1,nptl
        if(istptl(npts).ne.istxxx)goto 60
        if(   abs(idptl(npts)).ne.120
     *   .and.abs(idptl(npts)).ne.130
     *   .and.abs(idptl(npts)).ne.1120
     *   .and.abs(idptl(npts)).ne.1130
     *   .and.abs(idptl(npts)).ne.2230
     *   .and.abs(idptl(npts)).ne.2330
     *   .and.abs(idptl(npts)).ne.3331)goto 60
        ppt=pptl(1,npts)**2+pptl(2,npts)**2
        ppp=sqrt(ppt+pptl(3,npts)**2)
        ppt=sqrt(ppt)
        yyy=0.
        if(pptl(3,npts).ne.0..and.ppt.ne.0.)yyy=sign(1.,pptl(3,npts))*
     *   log((ppp+abs(pptl(3,npts)))/ppt)
c        if(ppp.gt.abs(pptl(3,npts)))then
c          yyy=.5*log((ppp+pptl(3,npts))/(ppp-pptl(3,npts)))
c        else
c          yyy=sign(100.,pptl(3,npts))
c        endif
        if(isdif.le.2)yyy=-sign(1.,float(ilprtg))*yyy   !trigger on antiproton (target : ilprt=-1)
c        if(idptl(npts).eq.-1120)then
          if(abs(pptl(3,npts)).gt.0.)then
            theta=sign(1.,float(ilprtg))*ppt/pptl(3,npts)
            if((isdif.le.2.and.theta.gt.2.5e-3.and.theta.lt.4.5e-3).or.
     *         (isdif.gt.2.and.theta.gt.0.2e-3.and.theta.lt.1.2e-3))then
              iii0=iii0+1
c              write(ifch,*)'la',
c              print *,'la',ilprtg,yyy,iii1,iii2
c     &      ,npts,idptl(npts),ppt,pptl(3,npts),theta,ityptl(npts)
            endif
c          endif
        endif
        if(isdif.eq.1)then
          if(yyy.gt.2.5   .and. yyy.lt.5.6)iii1=1
        elseif(isdif.eq.2)then
          if(yyy.gt.3.  .and. yyy.lt.5.6)iii1=1
          if(yyy.gt.-5.6 .and. yyy.lt.-4.4 )iii2=1
        elseif(isdif.eq.3)then
          if(yyy.gt.2.4  .and. yyy.lt.5.9)iii1=1
          if(yyy.gt.-4.2 .and. yyy.lt.1.1)iii3=1
          if(yyy.gt.-5.9 .and. yyy.lt.-2.4)iii2=1
          if(yyy.gt.-1.1 .and. yyy.lt.4.2)iii4=1
        elseif(isdif.eq.4)then
          if(ilprtg.eq.-1)then  !antiproton = target
            if(yyy.gt.2.4  .and. yyy.lt.5.9)iii1=1
            if(yyy.gt.-5.9 .and. yyy.lt.-3.2)iii2=iii2+1
          else                  !antiproton = projectile
            if(yyy.gt.-5.9 .and. yyy.lt.-2.4)iii1=1
            if(yyy.gt. 3.2 .and. yyy.lt.5.9)iii2=iii2+1
          endif
        elseif(isdif.eq.5)then
          if(yyy.gt.3.2  .and. yyy.lt.5.9)iii1=1
          if(yyy.gt.-5.9 .and. yyy.lt.-3.2)iii2=1
          if(abs(yyy).lt.2.4)Et1=Et1+ppt
          if(Et1.gt.0.2)iii3=1
          if(abs(yyy).gt.2.2 .and. abs(yyy).lt.4.2 )Et2=Et2+ppt
          if(Et2.gt.1.)iii4=1
        endif
60      continue
      enddo
      if(isdif.eq.1)then
        if(iii0.eq.1 .and. iii1.eq.1) isdiff=1
      elseif(isdif.eq.2)then
        if(iii0.eq.1 .and. iii1.eq.1 .and. iii2.ne.1) isdiff=1
      elseif(isdif.eq.3)then
        if( (iii1.ne.1 .and. iii3.eq.1) .or.
     &      (iii2.ne.1 .and. iii4.eq.1)       ) isdiff=1
      elseif(isdif.eq.4)then
        if(iii0.eq.1 .and. iii1.ne.1 .and. iii2.le.6) isdiff=1
c        if(isdiff.eq.1)then
c      do npts=1,nptl
c        if(istptl(npts).ne.istxxx)goto 80
c        if(   abs(idptl(npts)).ne.120
c     *   .and.abs(idptl(npts)).ne.130
c     *   .and.abs(idptl(npts)).ne.1120
c     *   .and.abs(idptl(npts)).ne.1130
c     *   .and.abs(idptl(npts)).ne.2230
c     *   .and.abs(idptl(npts)).ne.2330
c     *   .and.abs(idptl(npts)).ne.3331)goto 80
c        ppt=pptl(1,npts)**2+pptl(2,npts)**2
c        ppp=sqrt(ppt+pptl(3,npts)**2)
c        ppt=sqrt(ppt)
c        yyy=0.
c        if(pptl(3,npts).ne.0..and.ppt.ne.0.)yyy=sign(1.,pptl(3,npts))*
c     *   log((ppp+abs(pptl(3,npts)))/ppt)
c        print *,nrevt,yyy,idptl(npts),ityptl(npts)
c80      continue
c        enddo
c          print *,'dpe',iii0
c        endif
      elseif(isdif.eq.5)then
        if(iii1+iii2+iii3+iii4.eq.4) isdiff=1
      else
        stop'in sdiff. argument of sdiff not authorized.        '
      endif
             endif
      end

c----------------------------------------------------------------------
      subroutine xtrans(cvar,icas,inom,ifr,n)
c----------------------------------------------------------------------
      include "epos.incxan"

      character*25 cvar
      character*3 cend
      logical evte,nevte
      nend=index(cvar,' ')-1
      cend=cvar(max(1,nend-2):nend)
      evte=.false.
      if(cend.eq.'evt')evte=.true.
      nevte=.not.evte
      ifr=0
      ishift=0
      do i=1,8
        if(cvar(i:i+2).eq.'+10')ishift=10
        if(cvar(i:i+2).eq.'+15')ishift=15
        if(cvar(i:i+2).eq.'+20')ishift=20
        if(cvar(i:i+2).eq.'+25')ishift=25
        if(cvar(i:i+2).eq.'+30')ishift=30
        if(cvar(i:i+2).eq.'+35')ishift=35
        if(cvar(i:i+2).eq.'+40')ishift=40
        if(cvar(i:i+2).eq.'+45')ishift=45
      enddo
      if(cvar.eq.'numptl')then
        inom=1
      elseif(cvar.eq.'npaptl')then
        inom=2
      elseif(cvar.eq.'npmptl')then
        inom=3
      elseif(cvar.eq.'ispptl')then
        inom=4
      elseif(cvar.eq.'rapx')then
        inom=5
      elseif(cvar.eq.'iptlfr')then
        inom=6
      elseif(cvar.eq.'rinp')then
        inom=7
      elseif(cvar.eq.'eco')then
        inom=8
      elseif(cvar.eq.'tau')then
        inom=9
      elseif(cvar.eq.'ctr')then
        inom=10
      elseif(cvar.eq.'v2np')then
        inom=11
        istavar=1             !to switch on the calculation of "Standard variable"
      elseif(cvar.eq.'absrap')then
        inom=12
      elseif(cvar.eq.'rap')then
        inom=13
      elseif(cvar.eq.'xp')then
        inom=14
      elseif(cvar.eq.'xe')then
        inom=15
      elseif(cvar.eq.'pt')then
        inom=16
      elseif(cvar.eq.'p1a')then
        inom=17
      elseif(cvar.eq.'p2a')then
        inom=18
      elseif(cvar.eq.'xi')then
        inom=19
      elseif(cvar.eq.'xf')then
        inom=20
      elseif(cvar.eq.'t')then
        inom=21
      elseif(cvar.eq.'rapmi')then
        inom=22
      elseif(cvar.eq.'eta')then
        inom=23
      elseif(cvar.eq.'theta')then
        inom=24
      elseif(cvar.eq.'pt2')then
        inom=25
      elseif(cvar.eq.'et')then
        inom=26
      elseif(cvar.eq.'idptl')then
        inom=27
      elseif(cvar.eq.'istptl')then
        inom=28
      elseif(cvar.eq.'mass')then
        inom=29
      elseif(cvar.eq.'idaptl')then
        inom=30
      elseif(cvar.eq.'egy')then
        inom=31
      elseif(cvar.eq.'rapwro')then
        inom=32
      elseif(cvar.eq.'mt')then
        inom=33
      elseif(cvar.eq.'pplus')then
        inom=34
      elseif(cvar.eq.'pminus')then
        inom=35
      elseif(cvar.eq.'p5')then
        inom=36
      elseif(cvar.eq.'pa')then
        inom=37
      elseif(cvar.eq.'sob')then
        inom=38
      elseif(cvar.eq.'idpom')then
        inom=39
      elseif(cvar.eq.'p3a')then
        inom=40
      elseif(cvar.eq.'cmass')then
        inom=41
      elseif(cvar.eq.'arappi')then
        inom=42
      !--------------------------------------------------
      elseif(cvar(1:4).eq.'u2q2'
     .   .or.cvar(1:6).eq.'i-vnsp'.or.       ! vnsp variables
     .   nevte.and.cvar(1:4).eq.'vnsp')then    ! (no evt ending)
      !--------------------------------------------------
        inom=43
        istavar=1   
        ivnsp=1
        if(cvar.eq.'u2q2a')then
          icas=2
        elseif(cvar.eq.'u2q2c')then
          icas=3
        !- - - - - - - - - - - - - - - - - - - - - -
        !the following allows to plot event variables  
        !as a function of a parameter n (n=1-5) using:
        !bh hi iptl ... 0 0.5 5.5 4 trg iptl 1 5  
        !- - - - - - - - - - - - - - - - - - - - - - 
        elseif(cvar.eq.'i-vnsp-qAqB')then
          icas=12
        elseif(cvar.eq.'i-vnsp-qCqA')then
          icas=31
        elseif(cvar.eq.'i-vnsp-qCqB')then
          icas=32
        !- - - - - - - - - - - - - - - - - - - - - -
        !the following are normal ptl variables 
        !- - - - - - - - - - - - - - - - - - - - - -
        elseif(cvar.eq.'vnsp-uBqC-2-')then
          icas=232
        elseif(cvar.eq.'vnsp-uBqC-3-')then
          icas=233
        elseif(cvar.eq.'vnsp-uBqC-4-')then
          icas=234
        endif
      !--------------------------------------------------
      elseif(cvar(1:14).eq.'ncum-evt-drap-')then  ! net-cumlt variables 
      !--------------------------------------------------
        inom=44
        icas=ishift+1
        read(cvar(15:),'(i2)') Ndeta
        xpara(10,n)=Ndeta
        call actimupurp(n,ishift)
      !--------------------------------------------------
      !--------------------------------------------------
      elseif(cvar.eq.'itsptl')then
        inom=50
      elseif(cvar.eq.'ityptl'.or.cvar.eq.'ity25')then
        inom=51
        if(cvar.eq.'ityptl')then
          icas=1
        elseif(cvar.eq.'ity25')then !ityxx
          icas=2
        endif
      elseif(cvar.eq.'idoptl'.or.cvar.eq.'idobmes')then
        inom=52
        if(cvar.eq.'idoptl')then
          icas=1
        elseif(cvar.eq.'idobmes')then 
          icas=2
        endif
      elseif(cvar.eq.'iptl')then
        inom=53
      elseif(cvar.eq.'index')then
        inom=54
      elseif(cvar.eq.'p')then
        inom=55
        icas=10
      elseif(cvar.eq.'p1')then
        inom=55
        icas=1
      elseif(cvar.eq.'p2')then
        inom=56
      elseif(cvar.eq.'p3')then
        inom=57
      elseif(cvar.eq.'p4')then
        inom=58
      elseif(cvar.eq.'xg')then
        inom=59
      elseif(cvar.eq.'ek')then
        inom=60
      elseif(cvar.eq.'beta')then
        inom=61
      elseif(cvar.eq.'mt0')then
        inom=63
      elseif(cvar.eq.'qsqptl')then
        inom=64
      elseif(cvar.eq.'xelab')then
        inom=65
      elseif(cvar.eq.'hgtc05')then
        inom=66
        istavar=1             !to switch on the calculation of "Standard variable"
      elseif(cvar.eq.'hadtyp')then
        inom=67
        istavar=1
      elseif(cvar.eq.'hgtc1')then
        inom=68
        istavar=1
      elseif(cvar.eq.'x4')then
        inom=69
      elseif(cvar.eq.'npn')then
        inom=70
      elseif(cvar.eq.'routp')then
        inom=71
      elseif(cvar.eq.'hgtc3')then
        inom=72
        istavar=1
      elseif(cvar.eq.'mu14')then
        inom=73
        istavar=1
      elseif(cvar.eq.'delphi')then
        inom=74
        iok=0
        !------------------------------------------------------------
        !icorrtrig stores the histogram numbers of those histograms which
        !use the delphi variable (and therfore require a call corrtrig
        !------------------------------------------------------------
        do i=1,icorrtrig(0)
         if(icorrtrig(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          icorrtrig(0)=icorrtrig(0)+1
          if(icorrtrig(0).gt.mxxhis)stop'mxxhis too small'
        icorrtrig(icorrtrig(0))=n
        endif
      elseif(cvar.eq.'v2')then
        inom=75
      elseif(cvar.eq.'pt4')then
        inom=76
      elseif(cvar.eq.'rin')then
        stop'rin not used any more'
        inom=77
      elseif(cvar.eq.'theh1p')then
        inom=78
      elseif(cvar.eq.'theh1t')then
        inom=79
      elseif(cvar.eq.'phi')then
        stop'##### ERROR 12092014 (use: detphi) #####'
      elseif(cvar(1:6).eq.'detphi')then
        inom=80
        icas=nint(xpara(ishift+1,n))
        if(icas.gt.0)then
        call actidetphi(n,ishift)
        endif
      elseif(cvar.eq.'isoft')then
        inom=81
      elseif(cvar(1:3).eq.'mux'.and.cvar(1:6).ne.'muxevt')then
        inom=82
        call actimux(n,ishift)
        icas=ishift+1  
      elseif(cvar.eq.'mupurp')then
        inom=83
        icas=nint(xpara(ishift+1,n))
        call actimupurp(n,ishift)
      elseif(cvar.eq.'x3')then
        inom=84
      elseif(cvar.eq.'jorptl'.or.cvar.eq.'jor25')then
        inom=85
        if(cvar.eq.'jorptl')then
          icas=1
        elseif(cvar.eq.'jor25')then !jorxx
          icas=2
        endif
      elseif(cvar.eq.'ptlead')then
        inom=86
        iok=0
        !------------------------------------------------------------
        !icorrtrig stores the histogram numbers of those histograms which
        !use the ptlead variable (and therfore require a call corrtrig
        !------------------------------------------------------------
        do i=1,icorrtrig(0)
         if(icorrtrig(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          icorrtrig(0)=icorrtrig(0)+1
          if(icorrtrig(0).gt.mxxhis)stop'mxxhis too small'
        icorrtrig(icorrtrig(0))=n
        endif
      elseif(cvar.eq.'mu25')then
        inom=87
        istavar=1
      elseif(cvar.eq.'pai')then
        inom=88
        if(nint(xpara(1,n)).eq.0)stop'##### ERROR 18102014a #####'
        icas=nint(xpara(1,n)) *10000 + nint(xpara(2,n))
        ipairs1=1
      elseif(cvar.eq.'co2')then
        inom=89
        if(nint(xpara(1,n)).eq.0)stop'##### ERROR 18102014b #####'
        icas=nint(xpara(1,n)) *10000 + nint(xpara(2,n))
        ipairs1=1
      elseif(cvar.eq.'co3')then
        inom=90
        if(nint(xpara(1,n)).eq.0)stop'##### ERROR 18102014c #####'
        icas=nint(xpara(1,n)) *10000 + nint(xpara(2,n))
        ipairs1=1
      elseif(cvar.eq.'rad')then
        inom=91
      elseif(cvar.eq.'abseta')then
        inom=92
      elseif(cvar.eq.'phiexp')then
        inom=93
      elseif(cvar(1:3).eq.'epx'.and.cvar(1:6).ne.'epxevt')then
        inom=94
        icas=ishift+1 
        call actiepx(n,ishift)
      elseif(cvar.eq.'mu24')then
        stop'ERROR Use mux rather than mu24!!'  
        !inom=94
        !istavar=1
      elseif(cvar.eq.'rcastp')then
        inom=95
        istavar=1
      elseif(cvar.eq.'rcastt')then
        inom=96
        istavar=1
      elseif(cvar.eq.'corh'.or.cvar.eq.'ncorr')then !ncorr = old name
        inom=97
        iok=0
        do i=1,iCorH(0)
          if(iCorH(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          iCorH(0)=iCorH(0)+1
          iwww=iwww+1
          if(iwww.gt.mxxhis)stop'mxxhis too small'
          iCorH(iCorH(0))=n
          nhisxxx(n)=iwww
        endif
      elseif(cvar.eq.'nmass')then
        inom=98
      elseif(cvar.eq.'ekin')then
        inom=99
      elseif(cvar.eq.'mulevt')then
        inom=101
      elseif(cvar.eq.'etevt')then
        inom=102
      elseif(cvar.eq.'enevt')then
        inom=103
      elseif(cvar.eq.'ev6evt')then
        inom=104
      elseif(cvar.eq.'xenevt')then
        inom=105
      elseif(cvar.eq.'netevt')then
        inom=106
      elseif(cvar.eq.'ptevt')then
        inom=107
      elseif(cvar.eq.'pmxevt')then
        inom=108
      elseif(cvar.eq.'numevt')then
        inom=201
      elseif(cvar.eq.'egyevt')then
        inom=202
      elseif(cvar.eq.'bimevt')then
        inom=203
      elseif(cvar.eq.'xbjevt')then
        inom=204
      elseif(cvar.eq.'qsqevt')then
        inom=205
      elseif(cvar.eq.'yevt')then
        inom=206
      elseif(cvar.eq.'eloevt')then
        inom=207
      elseif(cvar.eq.'nd1evt'.or.cvar(1:3).eq.'nd-')then
        inom=208
        if(cvar.eq.'nd1evt')then
          icas=1
        elseif(cvar.eq.'nd-atlas-pt100-evt')then
          icas=17
        elseif(cvar.eq.'nd-atlas-pt500eta25-evt')then
          icas=18
        elseif(cvar.eq.'nd-atlas-pt500eta08-evt')then
          icas=19
        endif
      elseif(cvar.eq.'nd2evt')then
        inom=209
      elseif(cvar.eq.'theevt')then
        inom=210
      elseif(cvar.eq.'nspevt')then
        inom=211
      elseif(cvar.eq.'nhpevt')then
        inom=212
      elseif(cvar.eq.'sigtot')then
        inom=213
      elseif(cvar.eq.'sigela')then
        inom=214
      elseif(cvar.eq.'sloela')then
        inom=215
      elseif(cvar.eq.'nrgevt')then
        inom=216
      elseif(cvar.eq.'qevt')then
        inom=217
      elseif(cvar.eq.'qtlevt')then
        inom=218
      elseif(cvar.eq.'nd0evt')then
        inom=219
      elseif(cvar.eq.'threvt')then
        inom=220
        ifr=33                  !set thrust-frame
      elseif(cvar.eq.'omtevt')then
        inom=221
        ifr=33                  !set thrust-frame
      elseif(cvar.eq.'tmaevt')then
        inom=222
        ifr=33                  !set thrust-frame
      elseif(cvar.eq.'tmievt')then
        inom=223
        ifr=33                  !set thrust-frame
      elseif(cvar.eq.'oblevt')then
        inom=224
        ifr=33                  !set thrust-frame
      elseif(cvar.eq.'sphevt')then
        inom=230
        ifr=32                  !set sph-frame
      elseif(cvar.eq.'aplevt')then
        inom=231
        ifr=32                  !set sph-frame
      elseif(cvar.eq.'cpaevt')then
        inom=232
        ifr=34                  !set sph2-frame
      elseif(cvar.eq.'dpaevt')then
        inom=233
        ifr=34                  !set sph2-frame
      elseif(cvar.eq.'npoevt')then
        inom=234
      elseif(cvar.eq.'npnevt')then
        inom=235
      elseif(cvar.eq.'ikoevt')then
        inom=236
      elseif(cvar.eq.'iktevt')then
        inom=237
      elseif(cvar.eq.'npxevt')then
        inom=238
      elseif(cvar.eq.'nd6evt')then
        inom=239
      elseif(cvar.eq.'mu1evt')then
        inom=240
        istavar=1
      elseif(cvar.eq.'muievt')then
        inom=241
        istavar=1
      elseif(cvar.eq.'hgtevt')then
        inom=242
        istavar=1
      elseif(cvar.eq.'difevt')then
        inom=243
      elseif(cvar.eq.'dixevt')then
        inom=244
      elseif(cvar.eq.'nd7evt')then
        inom=245
      elseif(cvar.eq.'nd8evt')then
        inom=246
      elseif(cvar.eq.'nd9evt')then
        inom=247
      elseif(cvar.eq.'ndaevt')then
        inom=248
      elseif(cvar.eq.'ndbevt')then
        inom=249
      elseif(cvar.eq.'qinevt')then
        inom=250
      elseif(cvar.eq.'qfievt')then
        inom=251
      elseif(cvar.eq.'einevt')then
        inom=252
      elseif(cvar.eq.'efievt')then
        inom=253
      elseif(cvar.eq.'pinevt')then
        inom=254
      elseif(cvar.eq.'pfievt')then
        inom=255
      elseif(cvar.eq.'pxfevt')then    ! leading proton xf in cms
        inom=256
      elseif(cvar.eq.'pi+xf')then     ! pi+xf: pi+ yield at cms xf>0.01
        inom=257
      elseif(cvar.eq.'pi-xf')then     ! pi-xf: pi- yield at cms xf>0.01
        inom=258
      elseif(cvar.eq.'sigcut')then
        inom=260
      elseif(cvar.eq.'keu')then
        inom=261
      elseif(cvar.eq.'ked')then
        inom=262
      elseif(cvar.eq.'kes')then
        inom=263
      elseif(cvar.eq.'kolevt')then
        inom=265
      elseif(cvar.eq.'sigsd')then
        inom=266
      elseif(cvar.eq.'nglevt')then
        inom=267
      elseif(cvar.eq.'kppevt')then   ! collision numbers per participant
        inom=268
      elseif(cvar.eq.'npievt')then   ! pion + multiplicity per event
        inom=269
      elseif(cvar.eq.'np2evt')then   ! pion + multiplicity per participant
        inom=270
      elseif(cvar.eq.'sigdif'.or.cvar.eq.'sigdifr')then
        inom=271
      elseif(cvar.eq.'koievt')then
        inom=272
      elseif(cvar.eq.'ineevt')then
        inom=273
      elseif(cvar.eq.'elaevt')then
        inom=274
      elseif(cvar.eq.'itgevt')then
        inom=275
        iok=0
        do i=1,icorrtrig(0)
          if(icorrtrig(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          icorrtrig(0)=icorrtrig(0)+1
          if(icorrtrig(0).gt.mxxhis)stop'mxxhis too small'
          icorrtrig(icorrtrig(0))=n
        endif
      elseif(cvar.eq.'hrdevt')then
        inom=276
        iok=0
        do i=1,ihardevent(0)
          if(ihardevent(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          ihardevent(0)=ihardevent(0)+1
          if(ihardevent(0).gt.mxxhis)stop'mxxhis too small'
          ihardevent(ihardevent(0))=n
        endif
      elseif(cvar(2:6).eq.'j1evt'.or.cvar(2:6).eq.'j2evt')then
        iok=0
        do i=1,ijetfind1(0)
          if(ijetfind1(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          ijetfind1(0)=ijetfind1(0)+1
          if(ijetfind1(0).gt.mxxhis)stop'mxxhis too small'
          ijetfind1(ijetfind1(0))=n
        endif
        if(cvar.eq.'ej1evt')inom=277
        if(cvar.eq.'pj1evt')inom=278
        if(cvar(2:6).eq.'j2evt')then
          iok=0
          do i=1,ijetfind2(0)
            if(ijetfind2(i).eq.n)iok=1
          enddo
          if(iok.eq.0)then
            ijetfind2(0)=ijetfind2(0)+1
            if(ijetfind2(0).gt.mxxhis)stop'mxxhis too small'
            ijetfind2(ijetfind2(0))=n
          endif
          if(cvar.eq.'ej2evt')inom=279
          if(cvar.eq.'pj2evt')inom=280
        endif
      elseif(cvar.eq.'zppevt')then
        inom=281
      elseif(cvar.eq.'zptevt')then
        inom=282
      elseif(cvar.eq.'***not used***')then
        inom=283
      elseif(cvar.eq.'nd3evt')then
        inom=284
      elseif(cvar.eq.'nd4evt')then
        inom=285
      elseif(cvar.eq.'mubevt')then
        inom=286
        istavar=1
      elseif(cvar.eq.'nd5evt')then
        inom=287
      elseif(cvar.eq.'ekievt')then
        inom=288
      elseif(cvar.eq.'sd1evt')then
        inom=289
      elseif(cvar.eq.'sd2evt')then
        inom=290
      elseif(cvar.eq.'mdevt')then
        inom=291
        istavar=1     !to switch on the calculation of "Standard variable"
      elseif(cvar.eq.'m2devt')then
        inom=292
        istavar=1
      elseif(cvar.eq.'tdevt')then
        inom=293
        istavar=1
      elseif(cvar.eq.'ndpevt')then
        inom=294
      elseif(cvar.eq.'rapgap')then
        inom=295
        istavar=1
      elseif(cvar.eq.'ng1evt')then
        inom=296
      elseif(cvar.eq.'r21evt')then
        inom=297
      elseif(cvar.eq.'aimevt')then
        inom=301
      elseif(cvar.eq.'wjbevt')then
        inom=302
      elseif(cvar.eq.'njbevt')then
        inom=303
      elseif(cvar.eq.'djbevt')then
        inom=304
      elseif(cvar.eq.'tjbevt')then
        inom=305
      elseif(cvar.eq.'hjmevt')then
        inom=306
      elseif(cvar.eq.'ljmevt')then
        inom=307
      elseif(cvar.eq.'djmevt')then
        inom=308
      elseif(cvar.eq.'ybal')then
        inom=310
      elseif(cvar.eq.'yabal')then
        inom=310
      elseif(cvar.eq.'sigine')then
        inom=312
      elseif(cvar.eq.'sigiaa')then
        inom=313
      elseif(cvar.eq.'alpdsf')then
        inom=314
      elseif(cvar.eq.'alpdsh')then
        inom=315
      elseif(cvar.eq.'betdsf')then
        inom=316
      elseif(cvar.eq.'betdsh')then
        inom=317
      elseif(cvar.eq.'rexdip')then
        inom=318
      elseif(cvar.eq.'rexdit')then
        inom=319
      elseif(cvar.eq.'m14evt')then
        inom=320
        istavar=1
      elseif(cvar.eq.'ht3evt')then
        inom=321
      elseif(cvar.eq.'sigiex')then
        inom=322
      elseif(cvar.eq.'sigdex')then
        inom=323
      elseif(cvar.eq.'sigsex')then
        inom=324
      elseif(cvar.eq.'ekievt')then
        inom=325
      elseif(cvar.eq.'sigcaa')then
        inom=326
      elseif(cvar.eq.'sigtaa')then
        inom=327
      elseif(cvar.eq.'xkappa')then
        inom=328
      elseif(cvar.eq.'gamdsf')then
        inom=329
      elseif(cvar.eq.'gamdsh')then
        inom=330
      elseif(cvar.eq.'deldsf')then
        inom=331
      elseif(cvar.eq.'deldsh')then
        inom=332
      elseif(cvar.eq.'nd6evt')then
        inom=333
      elseif(cvar(1:6).eq.'muxevt')then
        inom=334
        call actimux(n,ishift)
        icas=ishift+1  
      elseif(cvar.eq.'typevt')then
        inom=335
      elseif(cvar.eq.'m25evt')then
        inom=339
        istavar=1
      elseif(cvar.eq.'segevt')then
        inom=340
      elseif(cvar.eq.'ielevt')then
        inom=341
      elseif(cvar.eq.'mc1evt')then
        inom=342
        istavar=1
      elseif(cvar.eq.'sdcdf')then
        inom=343
      elseif(cvar.eq.'dpecdf')then
        inom=344
      elseif(cvar.eq.'ddcdf')then
        inom=345
      elseif(cvar.eq.'phievt')then
        inom=346
      elseif(cvar.eq.'ndcevt')then
        inom=347
      elseif(cvar.eq.'jetevt')then
        inom=348
        call actijetevent(n)
        call actifastjet(n)
      elseif(cvar.eq.'epszer')then
        inom=349
      elseif(cvar.eq.'xsievt')then
        inom=350
        istavar=1
      elseif(cvar.eq.'xsicms')then
        inom=351
      elseif(cvar.eq.'calevt')then
        inom=352
        icaltrig(0)=icaltrig(0)+1
        icaltrig(icaltrig(0))=n
      elseif(cvar.eq.'fgpevt')then
        inom=353
        icaltrig(0)=icaltrig(0)+1
        icaltrig(icaltrig(0))=n
      elseif(cvar.eq.'bgpevt')then
        inom=354
        icaltrig(0)=icaltrig(0)+1
        icaltrig(icaltrig(0))=n
      elseif(cvar.eq.'gapevt')then
        inom=355
        icaltrig(0)=icaltrig(0)+1
        icaltrig(icaltrig(0))=n
      elseif(cvar.eq.'sigdd')then
        inom=356
      elseif(cvar.eq.'ajtevt')then
        inom=357
        call actijetevent(n)
        call actifastjet(n)
      elseif(cvar.eq.'fjtevt')then
        inom=358
        call actijetevent(n)
        call actifastjet(n)
      elseif(cvar.eq.'pjtevt')then
        inom=359
        call actijetevent(n)
        call actifastjet(n)
      elseif(cvar.eq.'styevt')then
        inom=360   !                  moved from 357
      elseif(cvar.eq.'ndsevt')then
        inom=361   !                  moved from 358
      elseif(cvar.eq.'m24evt')then
        inom=362   !                  moved from 359
        istavar=1
      elseif(cvar.eq.'ndhevt')then
        inom=363   !                  moved from 360
      elseif(cvar.eq.'ndfevt')then
        inom=364
      elseif(cvar.eq.'rcpevt')then
        inom=365
        istavar=1
      elseif(cvar.eq.'rctevt')then
        inom=366
        istavar=1
      elseif(cvar.eq.'mdievt')then
        inom=367
        istavar=1
      elseif(cvar.eq.'sig2ex')then
        inom=368
      elseif(cvar.eq.'fapevt')then
        inom=369
      elseif(cvar.eq.'sig2ex')then
        inom=372
      elseif(cvar.eq.'xppevt')then
        inom=494
      elseif(cvar.eq.'xpmevt')then
        inom=495
      elseif(cvar.eq.'xkpevt')then
        inom=496
      elseif(cvar.eq.'xkmevt')then
        inom=497
      elseif(cvar.eq.'pfrevt')then
        inom=498
      elseif(cvar.eq.'tfrevt')then
        inom=499
      elseif(cvar.eq.'ox1evt')then
        inom=501
      elseif(cvar.eq.'ox2evt')then
        inom=502
      elseif(cvar.eq.'ox3evt')then
        inom=503
      elseif(cvar.eq.'ox4evt')then
        inom=504
      elseif(cvar.eq.'eglevt')then  ! eccentricity
        inom=505
      elseif(cvar.eq.'fglevt')then  ! eccentricity_part
        inom=506
      elseif(cvar.eq.'rglevt')then  ! ratio ng2 / ng1
        inom=507
      elseif(cvar.eq.'sglevt')then  ! area S
        inom=508
      elseif(cvar.eq.'ptrevt')then
        inom=509
      elseif(cvar.eq.'rr2evt')then
        inom=510
        istavar=1             !to switch on the calculation of "Standard variable"
      elseif(cvar.eq.'perevt')then
        inom=511
      elseif(cvar.eq.'paievt')then
        inom=512
        if(nint(xpara(1,n)).eq.0)stop'##### ERROR 18102014d #####'
        icas=nint(xpara(1,n)) *10000 
        ipairs1=1
      elseif(cvar.eq.'co2evt')then
        inom=513
        if(nint(xpara(1,n)).eq.0)stop'##### ERROR 18102014e #####'
        icas=nint(xpara(1,n)) *10000 
        ipairs1=1
      elseif(cvar.eq.'co3evt')then
        inom=514
        if(nint(xpara(1,n)).eq.0)stop'##### ERROR 18102014f #####'
        icas=nint(xpara(1,n)) *10000 
        ipairs1=1
      elseif(cvar.eq.'nh1evt')then
        inom=515
      elseif(cvar.eq.'nh2evt')then
        inom=516
      elseif(cvar.eq.'nh3evt')then
        inom=517
      elseif(cvar.eq.'gbyevt')then
        inom=518
      elseif(cvar.eq.'icyevt')then
        inom=519
      elseif(cvar.eq.'xp9evt')then
        inom=520
      elseif(cvar.eq.'hlxevt')then
        inom=521
      elseif(cvar.eq.'icvevt')then
        inom=522
        imux(0)=imux(0)+1
        imux(imux(0))=n
      elseif(cvar.eq.'corhevt'
     .  .or.   cvar.eq.'corhtr'.or.cvar.eq.'trgevt')then ! corhtr,trgevt = old name 
        inom=523
        iok=0
        do i=1,iCorH(0)
          if(iCorH(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          iCorH(0)=iCorH(0)+1
          iwww=iwww+1
          if(iwww.gt.mxxhis)stop'mxxhis too small'
          iCorH(iCorH(0))=n
          nhisxxx(n)=iwww
        endif
      elseif(cvar.eq.'phoevt')then !bg
        inom=524
        iok=0
        do i=1,iphoton(0)
          if(iphoton(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          iphoton(0)=iphoton(0)+1
          if(iphoton(0).gt.mxxhis)stop'mxxhis too small'
          iphoton(iphoton(0))=n
        endif
      !------------------------------------------------
      elseif(cvar.eq.'qacevt'.or.
     .  evte.and.cvar(1:4).eq.'vnsp')then  ! vnsp variables (evt ending)
      !------------------------------------------------
        inom=525
        istavar=1  !call StandardVariables
        ivnsp=1
        if(cvar.eq.'qacevt')then !kept for the moment
          icas=1
        elseif(cvar.eq.'vnsp-qAqB-2-evt')then
          icas=122
        elseif(cvar.eq.'vnsp-qAqB-3-evt')then
          icas=123
        elseif(cvar.eq.'vnsp-qAqB-4-evt')then
          icas=124
        elseif(cvar.eq.'vnsp-qCqA-2-evt')then
          icas=312
        elseif(cvar.eq.'vnsp-qCqA-3-evt')then
          icas=313
        elseif(cvar.eq.'vnsp-qCqA-4-evt')then
          icas=314
        elseif(cvar.eq.'vnsp-qCqB-2-evt')then
          icas=322
        elseif(cvar.eq.'vnsp-qCqB-3-evt')then
          icas=323
        elseif(cvar.eq.'vnsp-qCqB-4-evt')then
          icas=324
        endif
      !--------------------------------------------------
      elseif(evte.and.cvar.eq.'ncum-evt')then  ! net-cumlt variables
      !--------------------------------------------------
        inom=526
        icas=ishift+1
        xpara(10,n)=1
        call actimupurp(n,ishift)
      !------------------------------------------------
      !------------------------------------------------
      elseif(cvar.eq.'ep2evt')then
        inom=530
        icas=1
      elseif(cvar.eq.'ep3evt')then
        inom=530
        icas=2
      elseif(cvar.eq.'ep4evt')then
        inom=530
        icas=3
      elseif(cvar.eq.'ep5evt')then
        inom=530
        icas=4
      elseif(cvar.eq.'rrrevt')then
        inom=530
        icas=5
      !------------------------------------------------
      !------------------------------------------------
      elseif(cvar.eq.'nu1evt')then
        inom=531
      elseif(cvar.eq.'nu2evt')then
        inom=532
      elseif(cvar.eq.'nu3evt')then
        inom=533
      elseif(cvar.eq.'phtevt')then !bg trigger for photon/hadron correlations
        inom=534
        iok=0
        do i=1,icorrtrig(0)
         if(icorrtrig(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          icorrtrig(0)=icorrtrig(0)+1
          if(icorrtrig(0).gt.mxxhis)stop'mxxhis too small'
        icorrtrig(icorrtrig(0))=n
        endif
      elseif(cvar.eq.'mpoevt')then
        inom=535
      elseif(cvar.eq.'ap2evt')then
        inom=536
      elseif(cvar.eq.'aq2evt')then
        inom=537
      elseif(cvar.eq.'ar2evt')then
        inom=538
      elseif(cvar.eq.'corevt')then
        inom=539
       elseif(cvar.eq.'epxevt')then
        inom=540
        icas=ishift+1
        call actiepx(n,ishift)
      elseif(cvar.eq.'ue0evt')then
        inom=541
      elseif(cvar.eq.'ue1evt')then
        inom=542
      elseif(cvar.eq.'ue2evt')then
        inom=543
      elseif(cvar.eq.'ue3evt')then
        inom=544
      elseif(cvar.eq.'ue4evt')then
        inom=545
      elseif(cvar.eq.'ue5evt')then
        inom=546
      elseif(cvar.eq.'ue6evt')then
        inom=547
      elseif(cvar(1:6).eq.'mupevt')then
        inom=548
        icas=ishift+1
        call actimupurp(n,ishift)
      elseif(cvar.eq.'nddevt')then  !ALICE inel 
        inom=549
      elseif(cvar(3:6).eq.'zevt'.or.cvar(4:7).eq.'zevt')then
        inom=551
        if(cvar(3:6).eq.'zevt')read(cvar(2:2),*)ij
        if(cvar(4:7).eq.'zevt')read(cvar(2:3),*)ij
        if(cvar(1:1).eq.'x')icas=-ij
        if(cvar(1:1).eq.'z')icas=ij
      else
        print *,' '
        print *,'              xtrans: unknown variable ',cvar
        print *,' '
c       inom=-1
        stop
      endif
      end

c----------------------------------------------------------------------
      subroutine xval(n,icas,inom,lf,j,x)
c----------------------------------------------------------------------
c   n ...... histogram index
c   inom ... variable index
c              1-100 particle variables
c              101-200 accumulative event variables
c              > 200 other event variables
c   lf ..... frame index
c   j ..... particle index (used for particle variables)
c----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incems"
      include "epos.incho"
      include "epos.incxan"
      common/zeus2/qtl /cgbyjmax/gbyjmax
      parameter (ntim=1000)
      common/cprt/pprt(5,ntim),q2prt(ntim),idaprt(2,ntim),idprt(ntim)
     &,iorprt(ntim),jorprt(ntim),nprtj

      common/cxyzt/xptl(mxptl),yptl(mxptl),zptl(mxptl)
     * ,tptl(mxptl),optl(mxptl),uptl(mxptl),sptl(mxptl)
     *,rptl(mxptl,3)
      parameter(kkkmax=6)
      common/cpairs/paievt(kkkmax),co2evt(kkkmax),co3evt(kkkmax)
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall
      parameter (mxfra=5)
      common/pfra/nfra,ifra(mxfra),ivfra(2,mxhis),itfra(mxtri,mxhis)
     $     ,imofra(3,mxfra),iffra(mxfra),r1fra(3,mxfra),r2fra(3,mxfra)
     $     ,emax(mxfra)
      common/cphi2/phi2pos,phi2neg
      common/cen/ncentr  
      parameter(nbkbin=40)
      common /kfitd/ xkappafit(nbkbin,nclegy,nclha,nclha),xkappa,bkbin
      parameter (mxpaih=10)
      parameter (maxpt=50)
      common/cpairs2/
     .    paih(kkkmax,mxpaih,maxpt)
     .   ,co2h(kkkmax,mxpaih,maxpt),co3h(kkkmax,mxpaih,maxpt)
      parameter (maxcorr=41)        !bg
      common/cpairs3/corr(maxcorr) !bg
      double precision bofra,bofra1,bofra2,bofra3,bofra4,bofra5
      common/dfra/bofra(5,mxfra)
      dimension p(5,mxfra),aimuni(10,mxhis),xor(5,mxfra)
      common/cranphi/ranphi  /cicentrality/icentrality
      double precision xxx
      logical NoLongLivParent
      save p,aimuni,xor,ishift

      istxxx=0
      if(istu1(istuse(n)).gt.-999)istxxx=istu1(istuse(n))

      if(j.gt.0)then
        idepos=ideposf( 3 ,j)
      endif

      call getItyJor25(j,ityxx,jorxx) !ity,jor of 25 parton

      phinll=phievt+ranphi
      if(phinll.lt.-pi)phinll=phinll+2*pi
      if(phinll.gt.pi)phinll=phinll-2*pi
      if(iffra(lf).eq.0.and.j.gt.0)then
        do l=1,5
          p(l,lf)=pptl(l,j)
        enddo
        do l=1,4
          xor(l,lf)=xorptl(l,j)
        enddo
        if(imofra(1,lf).ne.0)then
          call utrota(imofra(1,lf),r1fra(1,lf),r1fra(2,lf),r1fra(3,lf)
     $         ,p(1,lf),p(2,lf),p(3,lf))
          call utrota(imofra(1,lf),r1fra(1,lf),r1fra(2,lf),r1fra(3,lf)
     $         ,xor(1,lf),xor(2,lf),xor(3,lf))
        endif
        if(imofra(2,lf).ne.0)then !the x-z exchanged is ok !!
          call utrota(imofra(2,lf),r2fra(3,lf),r2fra(2,lf),r2fra(1,lf)
     $         ,p(3,lf),p(2,lf),p(1,lf))
          call utrota(imofra(2,lf),r2fra(3,lf),r2fra(2,lf),r2fra(1,lf)
     $         ,xor(3,lf),xor(2,lf),xor(1,lf))
        endif
        if(imofra(3,lf).ne.0)then
          imof3=sign(1,imofra(3,lf))
          if(abs(imofra(3,lf)).gt.1)then
            bofra1=0d0
            bofra2=0d0
            bofra5=1d0
            call utlob5(imof3*yhaha
     $                 ,p(1,lf),p(2,lf),p(3,lf),p(4,lf),p(5,lf))
            bofra3=bofra(1,lf)
            bofra4=bofra(2,lf)
            call utlob4(imof3,bofra1,bofra2,bofra3,bofra4,bofra5
     $         ,xor(1,lf),xor(2,lf),xor(3,lf),xor(4,lf))
            bofra3=bofra(3,lf)
            bofra4=bofra(4,lf)
            bofra5=bofra(5,lf)
          else
            bofra1=bofra(1,lf)
            bofra2=bofra(2,lf)
            bofra3=bofra(3,lf)
            bofra4=bofra(4,lf)
            bofra5=bofra(5,lf)
          endif
          call utlob4(imof3,bofra1,bofra2,bofra3,bofra4,bofra5
     $         ,p(1,lf),p(2,lf),p(3,lf),p(4,lf))
          call utlob4(imof3,bofra1,bofra2,bofra3,bofra4,bofra5
     $         ,xor(1,lf),xor(2,lf),xor(3,lf),xor(4,lf))
        endif
        iffra(lf)=1
      endif

c--------------------------------- 1 - 100 ----------------------------
      if(inom.eq.1)then
        x=1.
        if(idepos.eq.9998 .or.idepos.eq.9997) then !bg weight for photons linked to enhanced production
          x=1./desptl(j)
        endif
      elseif(inom.eq.2)then
        x=isign(1,idepos)
      elseif(inom.eq.3)then
        chrg=0
        if(iabs(idepos).le.9999
     $       .and.mod(iabs(idepos),10).le.1)
     $       call idchrg( 11 ,idepos,chrg)
        if(chrg.eq.0.)then
          x=0
        else
          x=int(sign(1.,chrg))
        endif
      elseif(inom.eq.4)then
        iad=abs(idepos)
        jspin=mod(iad,10)
        x=0.
        if (iad.ge.100.and.iad.lt.1000) x=1./(1.+2.*jspin)
        if (iad.ge.1000.and.iad.lt.9999) x=1./(2.+2*jspin)
      elseif(inom.eq.5)then    !'rapx'  !st-rap for string segments only !!!!!!!!!!
        x=dezptl(j)
      elseif(inom.eq.6)then                                      !'iptlfr'
        x=0
        if(j.ge.minfra.and.j.le.maxfra)x=1
      elseif(inom.eq.7)then                                        !'rinp'
        aa=cos(phinll)
        bb=sin(phinll)
        x=xptl(j)*aa+yptl(j)*bb
      elseif(inom.eq.8)then                        !'eco' !engy in comoving frame
        x=0
        amt=p(5,lf)**2+p(1,lf)**2+p(2,lf)**2
        if(amt.gt.0..and.p(4,lf)+abs(p(3,lf)).gt.0.d0)then
          amt=sqrt(amt)
          rap=sign(1.,p(3,lf))*alog((p(4,lf)+abs(p(3,lf)))/amt)
          rapx=dezptl(j)
          x=amt*cosh(rap-rapx)
        endif
      elseif(inom.eq.9)then                                       !'tau'
        x=-999999
        !if(iorptl(j).ne.0)then
        ! jo=iorptl(j)
         dt=xorptl(4,j)  !-xorptl(4,jo)
         dz=xorptl(3,j)  !-xorptl(3,jo)
         x2=dt**2-dz**2
         if(x2.gt.0.)x=sqrt(x2)
        !endif
      elseif(inom.eq.10)then                                       !'ctr'
        x=ctrevt
      elseif(inom.eq.11)then                                       !'v2np'
        phi=polar( p(1,lf) , p(2,lf) )
        pt=sqrt(p(2,lf)**2+p(1,lf)**2)
        eta=0
        if(p(3,lf).ne.0..and.pt.ne.0.)eta=sign(1.,p(3,lf))*
     *       alog((sqrt(p(3,lf)**2+pt**2)+abs(p(3,lf)))/pt)
        if(eta.gt.0)then
        phi2=phi2neg
        else
        phi2=phi2pos
        endif
        if(maproj.gt.2.and.matarg.le.2)phi2=phi2pos
        if(maproj.le.2.and.matarg.gt.2)phi2=phi2neg
        x=cos(2*(phi-phi2))
      elseif(inom.eq.12)then                                      !'absrap'
        amt=p(5,lf)**2+p(1,lf)**2+p(2,lf)**2
        if(amt.gt.0..and.p(4,lf)+abs(p(3,lf)).gt.0.d0)then
          amt=sqrt(amt)
          x=alog((p(4,lf)+abs(p(3,lf)))/amt)
        else
          x=0.                  !
        endif
      elseif(inom.eq.13)then                                        !'rap'
        amt=p(5,lf)**2+p(1,lf)**2+p(2,lf)**2
        if(amt.gt.0..and.p(4,lf)+abs(p(3,lf)).gt.0.d0)then  
          amt=sqrt(amt)
          x=sign(1.,p(3,lf))*log((p(4,lf)+abs(p(3,lf)))/amt)  
        else
          x=1e10               !
        endif
      elseif(inom.eq.14)then                                         !'xp'
        emx=emax(lf)
        if(iappl.eq.6)emx=sqrt((emx+p(5,lf))*(emx-p(5,lf)))
        x=sqrt(p(3,lf)**2+p(2,lf)**2+p(1,lf)**2)/emx
      elseif(inom.eq.15)then                                         !'xe'
        x=min(1.,p(4,lf)/emax(lf))
      elseif(inom.eq.16)then                                         !'pt'
        x=sqrt(p(2,lf)**2+p(1,lf)**2)
      elseif(inom.eq.17)then
        x=abs(p(1,lf))
      elseif(inom.eq.18)then
        x=abs(p(2,lf))
      elseif(inom.eq.19)then
        x=-log(sqrt(p(3,lf)**2+p(2,lf)**2+p(1,lf)**2)/emax(lf))
      elseif(inom.eq.20)then                                     !'xf'
        m=mod(ifra(lf)/10,10)
        if(m.eq.1.or.noebin.lt.0)then
c          pmax=sqrt((engy/2)**2-prom*2)
          pmax=pnullx               !???????????????????
          if(mod(ifra(lf),10).eq.2)pmax=pnll
          x=p(3,lf)/pmax
        else
          x=p(3,lf)/emax(lf)
        endif
c        if(x.gt.0.95.and.idepos.eq.1220)then
c          write(ifch,'(a,d25.15)')'ici !!!!!!!!!',seedc
c          stop
c        endif
      elseif(inom.eq.21)then
c        pmax=pmxevt
c        pmax=sqrt((engy/2)**2-prom*2)
        pmax=pnullx             !???????????????????
        if(mod(ifra(lf),10).eq.2)pmax=pnll
        x=-(amproj**2-2.*sqrt(amproj**2+pmax**2)*p(4,lf)
     *      +2.*abs(pmax*p(3,lf))+p(5,lf)**2)
      elseif(inom.eq.22)then
        amt=sqrt(p(5,lf)**2+p(1,lf)**2+p(2,lf)**2)
        if(amt.ne.0.)then
          x=-sign(1.,p(3,lf))*alog((p(4,lf)+abs(p(3,lf)))/amt)
        else
          x=0.                  !
        endif
      elseif(inom.eq.23)then                                     !'eta'
        pt=sqrt(p(2,lf)**2+p(1,lf)**2)
        if(p(3,lf).eq.0.)then
          x=0.
        elseif(pt.ne.0.)then
          x=sign(1.,p(3,lf))*
     *       alog((sqrt(p(3,lf)**2+pt**2)+abs(p(3,lf)))/pt)
        else
          x=sign(1000.,p(3,lf))
        endif
      elseif(inom.eq.24)then                                     !'theta (deg)'
        pt=sqrt(p(2,lf)**2+p(1,lf)**2)
        x=90
        if(p(3,lf).ne.0.)x=atan(pt/p(3,lf))/pi*180.
        if(x.lt.0.)x=180.+x
      elseif(inom.eq.25)then                                     !'pt2'
        x=p(2,lf)**2+p(1,lf)**2
      elseif(inom.eq.26)then                                     !'et'
        pt=sqrt(p(2,lf)**2+p(1,lf)**2)
        x=0
        eef=p(4,lf)
c        if(idepos.ge.1000)eef=eef-prom
c        if(idepos.le.-1000)eef=eef+prom
        p2=p(3,lf)**2+p(2,lf)**2+p(1,lf)**2
        if(p2.ne.0.)x=eef*pt/sqrt(p2)
      elseif(inom.eq.27)then                                   !'idptl'
        x=idepos
      elseif(inom.eq.28)then                                   !'istptl'
        x=istptl(j)
      elseif(inom.eq.29)then                                   !'mass'
        x=p(5,lf)
        if(istptl(j).le.1)call idmass(idepos,x)
      elseif(inom.eq.30)then                                   !'idaptl'
        x=abs(idepos)
      elseif(inom.eq.31)then                                   !'egy'
        x=egyevt
      elseif(inom.eq.32)then                                   !'rapwro'
        x=0
        pt2=p(2,lf)**2+p(1,lf)**2
        if(p(3,lf).ne.0.)x=sign(1.,p(3,lf))*
     *       alog((sqrt(p(3,lf)**2+pt2+.13957**2)+abs(p(3,lf)))
     *       /sqrt(pt2+.13957**2))
      elseif(inom.eq.33)then                                  !'mt'
        x=sqrt(p(2,lf)**2+p(1,lf)**2+p(5,lf)**2)
      elseif(inom.eq.34)then                                  !'pplus'
        x=sign(1.,p(3,lf)) * (p(4,lf)+abs(p(3,lf)))
      elseif(inom.eq.35)then                                  !'pminus'
        x=sign(1.,p(3,lf)) * (p(4,lf)-abs(p(3,lf)))
      elseif(inom.eq.36)then                                  !'p5' (mass)
        x=p(5,lf)
      elseif(inom.eq.37)then                                  !pa
        x=sqrt(p(1,lf)**2+p(2,lf)**2+p(3,lf)**2)
      elseif(inom.eq.38)then                                  !'pa'
        if(p(1,lf)**2+p(2,lf)**2+p(3,lf)**2.ne.0)
     *       x=egyevt**2/sqrt(p(1,lf)**2+p(2,lf)**2+p(3,lf)**2)*p(4,lf)
      elseif(inom.eq.39)then                                  !idpom
        x=idepos/1000000
      elseif(inom.eq.40)then                                  !p3a
        x=abs(p(3,lf))
      elseif(inom.eq.41)then
        cm2=p(4,lf)**2-p(3,lf)**2-p(2,lf)**2-p(1,lf)**2         !cmass
        x=sign(sqrt(abs(cm2)),cm2)
      elseif(inom.eq.42)then    !arappi
        x=0
        pt2=p(2,lf)**2+p(1,lf)**2
        if(p(3,lf).ne.0.)
     *       x=alog((sqrt(p(3,lf)**2+pt2+.13957**2)+abs(p(3,lf)))
     *       /sqrt(pt2+.13957**2))
      !------------------------------------------------
      elseif(inom.eq.43)then      ! vnsp variables 
      !------------------------------------------------
        select case (icas)
        case(2) 
          phi=polar(p(1,lf),p(2,lf))
          x = cos(2*phi)*x2aevt +  sin(2*phi)*y2aevt            !'u2q2a'
        case(3)
          phi=polar(p(1,lf),p(2,lf))
          x = cos(2*phi)*x2cevt +  sin(2*phi)*y2cevt            !'u2q2c'
        !- - - - - - - - - - - - - -
        !In the following we fill a 
        !histogram with j=1,2,3,4,5
        !In optns use "trg iptl 1 5" 
        !- - - - - - - - - - - - - - 
        case(12)
          if(j.gt.5)stop'ERROR 19092020x'
          m=j
          x=qvnsp(1,m,1)*qvnsp(2,m,1)+qvnsp(1,m,2)*qvnsp(2,m,2)  !'i-vnsp-qAqB'
        case(31) 
          if(j.gt.5)stop'ERROR 19092020y'
          m=j
          x=qvnsp(3,m,1)*qvnsp(1,m,1)+qvnsp(3,m,2)*qvnsp(1,m,2)  !'i-vnsp-qCqA'
        case(32) 
          if(j.gt.5)stop'ERROR 19092020z'
          m=j 
          x=qvnsp(3,m,1)*qvnsp(2,m,1)+qvnsp(3,m,2)*qvnsp(2,m,2)  !'i-vnsp-qCqB'
        !- - - - - - - - - - - - - - - - - - - - - -
        !the following are normal ptl variables 
        !- - - - - - - - - - - - - - - - - - - - - -
        case(232)
          m=2
          phi=polar(p(1,lf),p(2,lf))
          x = cos(m*phi)*qvnsp(3,m,1) +  sin(m*phi)*qvnsp(3,m,2)   !'vnsp-uBqC-2'
        case(233)
          m=3
          phi=polar(p(1,lf),p(2,lf))
          x = cos(m*phi)*qvnsp(3,m,1) +  sin(m*phi)*qvnsp(3,m,2)   !'vnsp-uBqC-3'
        case(234)
          m=4
          phi=polar(p(1,lf),p(2,lf))
          x = cos(m*phi)*qvnsp(3,m,1) +  sin(m*phi)*qvnsp(3,m,2)   !'vnsp-uBqC-4'
        end select
      !--------------------------------------------------------
      elseif(inom.eq.44)then   ! net-cumlt variables (etagaps) 
      !-------------------------------------------------------- 
        !-----------------------------------------
        !In the following we fill a histogram with 
        !j=1,2,3,...,jmax corresponding to |eta| gaps
        !In optns use "trg iptl 1 jmax" 
        !-----------------------------------------
        if(j.le.xpara(10,n))then
          ishift=icas-1
          x=ypara(ishift+j,n)
        else
          print*," /!\ WARNING : j-value differs between "
     .          ,"'iptl i j' and 'ncum-evt-drap-j' in analyses"
        endif
      !-------------------------------------------------- 
      !-------------------------------------------------- 
      elseif(inom.eq.50)then
        x=itsptl(j)
      elseif(inom.eq.51)then !'ityptl' 'ity25'
        select case (icas)
        case(1) 
          x=ityptl(j)
        case(2) 
          x=ityxx
        end select
      elseif(inom.eq.52)then   !'idoptl' 'idobmes'
        select case (icas)
        case(1) 
          x=0.
          if(iorptl(j).gt.0) x=idptl(iorptl(j))
        case(2) 
          x=0.
          if(.not.NoLongLivParent(2,j))x=1.
        end select
      elseif(inom.eq.53)then
        x=j
      elseif(inom.eq.54)then                       !'sloela'
        call idflav(idepos,ifl1,ifl2,ifl3,jspin,0)
        x=0 !???????????????????????
      elseif(inom.eq.55.and.icas.eq.10)then              !'p'
        x=p(1,lf)**2+p(2,lf)**2+p(3,lf)**2
        x=sqrt(x)
      elseif(inom.eq.55.and.icas.eq.1)then              !'p1'
        x=p(1,lf)
      elseif(inom.eq.56)then                       !'p2'
        x=p(2,lf)
      elseif(inom.eq.57)then                       !'p3'
        x=p(3,lf)
      elseif(inom.eq.58)then                       !'p4'
        x=p(4,lf)
      elseif(inom.eq.59)then                       !E/p_max
c        pmax=sqrt((engy/2)**2-prom*2)
          pmax=pnullx               !???????????????????
        if(mod(ifra(lf),10).eq.2)pmax=pnll
        x=p(4,lf)/pmax
      elseif(inom.eq.60)then                       !'ek'
        x=p(4,lf)-p(5,lf)
      elseif(inom.eq.61)then                       !'beta'
        x=p(3,lf)/p(4,lf)
      elseif(inom.eq.63)then                       !'mt0'
        x=sqrt(p(2,lf)**2+p(1,lf)**2+p(5,lf)**2)-p(5,lf)
      elseif(inom.eq.64)then                       !qsqptl
        x=qsqptl(j)
      elseif(inom.eq.65)then                       !xelab=Elab/Eolab
        x=p(4,lf)/(ecms**2/2/prom-prom)
        if(x.gt.0.9999999) x=.9999999
      elseif(inom.eq.66)then    !'hgtc05' ... charged ptl mult |[c]|<0.5
        x=multc05
      elseif(inom.eq.67)then    !'hadtyp' ... primary (1) or secondary (2) hadron
        if(j.le.nbdky)then
          x=1
        else
          x=2
        endif
      elseif(inom.eq.68)then    !'hgtc1'
        x=multc1
      elseif(inom.eq.69)then                       !'x4'
        x=xor(4,lf)
      elseif(inom.eq.70)then                       !'npn'
        x=ng1evt   !npjevt+ntgevt    !npn
      elseif(inom.eq.71)then                       !'routp'
        cc=-sin(phinll)
        dd=cos(phinll)
        x=xptl(j)*cc+yptl(j)*dd
      elseif(inom.eq.72)then    !'hgtc3' ... charged ptl mult |eta|<3.15  /6.3
        x=multc3/6.3
      elseif(inom.eq.73)then    !'mu14' ... charged ptl mult |eta|<1  pt>.4
        x=multc14
      elseif(inom.eq.74)then    !'delphi' ... azimuthhal correlation or 'xe' xe=-pt_asso*cos(delphi)/pt_trig !bg
        x=10000.
        pt=sqrt(p(1,lf)**2+p(2,lf)**2)
        if(nint(ypara(1,n)).ne.0.and.j.ne.nint(ypara(1,n)).and. !bg j=particle number
     $           pt.gt.0)then
           phi=sign(1.,p(2,lf))*acos(p(1,lf)/pt)
           x=phi-ypara(2,n)
           phiz= ypara(3,n)
           pttrig=ypara(4,n)
           if   (x.lt.(-2+phiz)*pi)then
            x=x+4*pi
           elseif(x.lt.(0+phiz)*pi)then
            x=x+2*pi
           elseif(x.gt.(4+phiz)*pi)then
            x=x-4*pi
           elseif(x.gt.(2+phiz)*pi)then
            x=x-2*pi
           endif
           if(nint(ypara(5,n)).eq.1) then                     !xe without underlying event correction
             if(x.gt.2.*pi/3. .and.x.lt.4.*pi/3.) then
               x=-pt*cos(x)/pttrig
             else
               x=-999.
             endif
           elseif(nint(ypara(5,n)).eq.2) then                 !xe with  underlying event correction
             if(x.gt.4.*pi/3. .and.x.lt.5.*pi/3.) then
               phirand=4.*pi/3.-rangen()*2.*pi/3.       !bg choose phi randomly between [2pi/3,4pi/3]. done by N.Arbor for Alice collaboration
               x=-pt*cos(phirand)/pttrig 
             elseif(x.gt.pi/3. .and.x.lt.2.*pi/3.) then
               phirand=4.*pi/3.-rangen()*2.*pi/3.
               x=-pt*cos(phirand)/pttrig 
             else
               x=-999.
             endif 
           endif
        endif
      elseif(inom.eq.75)then    !'v2'
c        if(iranphi.ne.1)stop'##### ERROR 29062010b #####'
        aa=cos(phinll)
        bb=sin(phinll)
        cc=-sin(phinll)
        dd=cos(phinll)
        px=p(1,lf)*aa+p(2,lf)*bb
        py=p(1,lf)*cc+p(2,lf)*dd
        pt2=p(2,lf)**2+p(1,lf)**2
        x=0
        if(pt2.gt.0.)x=(px**2-py**2)/pt2
      elseif(inom.eq.76)then                                     !'pt4'
        x=(p(2,lf)**2+p(1,lf)**2)**2
      elseif(inom.eq.77)then                                     !'rin'
        stop'rin not used any more'
        x=rinptl(j)
      elseif(inom.eq.78)then              !'theta for H1 (rad) for (proj side)'
        pt=sqrt(p(2,lf)**2+p(1,lf)**2)
        p1=p(1,lf)
        p2=p(2,lf)
        p3=p(3,lf)
        p4=p(4,lf)
        p5=p(5,lf)
        if(abs(p3).gt.1e-5)then
c put the particle in the projectile frame
          call utlob5(yhaha,p1,p2,p3,p4,p5)
c put the particle in a frame where the projectile (proton) has 820 GeV (HERA)
          call utlob4(-1,0d0,0d0,819.99946d0,820.d0,0.938d0,p1,p2,p3,p4)
          x=atan(pt/p3)
          if(x.lt.0.)x=pi+x
        else
          x=0.5*pi
        endif
      elseif(inom.eq.79)then             !'theta for H1 (rad) (for target side)'
        pt=sqrt(p(2,lf)**2+p(1,lf)**2)
        p1=p(1,lf)
        p2=p(2,lf)
        p3=p(3,lf)
        p4=p(4,lf)
        p5=p(5,lf)
        if(abs(p3).gt.1e-5)then
c put the particle in the projectile frame
          call utlob5(-yhaha,p1,p2,p3,p4,p5)
c put the particle in a frame where the projectile (proton) has 820 GeV (HERA)
          call utlob4(-1,0d0,0d0,-819.99946d0,820d0,0.938d0,p1,p2,p3,p4)
          x=atan(pt/p3)
          if(x.gt.0.)x=pi-x
        else
          x=0.5*pi
        endif
      elseif(inom.eq.80)then    !'detphi'
        if(icas.gt.0)then
          x=yyarr( nhisxxx(n) , j )
        else
          x=1000
          pt=sqrt(p(1,lf)**2+p(2,lf)**2)
          if(pt.gt.0.)then
             phi=sign(1.,p(2,lf))*acos(p(1,lf)/pt)
             x=phi-phinll
          endif
          if(x.lt.-pi)x=x+2*pi
          if(x.gt.pi)x=x-2*pi
        endif
      elseif(inom.eq.81)then  !'isoft'
        x=0
        it=ityptl(j)
        if(it.ge.20.and.it.le.29)x=1
        if(it.ge.40.and.it.le.60)x=1
      elseif(inom.eq.82)then  !'mux' ... charged ptl mult 
        x= ypara(icas,n)        
      elseif(inom.eq.83)then    !'mupurp'
cgs-start
        if(icas.eq.1004)then
          phi=polar(p(1,lf),p(2,lf))
          xn=ypara(1,n)
          yn=ypara(2,n)
          iord=ypara(3,n)
          x=cos(iord*phi)*xn+sin(iord*phi)*yn 
        else
          x=yyarr( nhisxxx(n) , j )
        endif
cgs-end
      elseif(inom.eq.84)then                       !'x3'
        x=xor(3,lf)
      elseif(inom.eq.85)then   !'jorptl' 'jor25'
        select case (icas)
        case(1) 
          x=jorptl(j)
        case(2) 
          x=jorxx
        end select
      elseif(inom.eq.86)then    !'ptlead' ... pt of particle with higher pt
        x=0.
        if(nint(ypara(1,n)).ne.0)x=ypara(4,n)
      elseif(inom.eq.87)then    !'mu25' ... charged ptl mult |eta|<2.5  pt>.5
        x=multc25
      elseif(inom.eq.88.or.inom.eq.89.or.inom.eq.90)then !'pai''co2''co3' 
        kkk=icas/10000
        kk1=kkk/10  /8   +1
        kk2=mod(kkk,10)
        kkk=(kk1-1)*(kkkmax/2)+kk2
        id=mod(icas,10000)
        x=0
        iid=0
        if(id.eq.120)iid=1
        if(id.eq.1120)iid=2
        if(id.eq.130)iid=3
        if(id.eq.1130)iid=4
        if(id.eq.20)iid=5
        if(id.eq.2130)iid=6
        if(id.eq.2330)iid=7
        if(id.eq.3331)iid=8
        if(id.eq.331)iid=9
        if(id.eq.9970)iid=10
        if(iid.ne.0.and.inom.eq.88)x=paih(kkk,iid,j)
        if(iid.ne.0.and.inom.eq.89)x=co2h(kkk,iid,j)
        if(iid.ne.0.and.inom.eq.90)x=co3h(kkk,iid,j)
      elseif(inom.eq.91)then                            !'rad'
        x=0.001       !unit is fm !
        x1=xor(1,lf)
        x2=xor(2,lf)
        xd=dble(x1)**2+dble(x2)**2
        if(xd.gt.0.d0.and.xd.eq.xd)x=sqrt(xd)
      elseif(inom.eq.92)then                                     !'abseta'
        pt=sqrt(p(2,lf)**2+p(1,lf)**2)
        if(p(3,lf).eq.0.)then
          x=0.
        elseif(pt.ne.0.)then
          x=sign(1.,p(3,lf))*
     *       alog((sqrt(p(3,lf)**2+pt**2)+abs(p(3,lf)))/pt)
        else
          x=sign(1000.,p(3,lf))
        endif
c        pp=sqrt(p(2,lf)**2+p(1,lf)**2+p(3,lf)**2)
c        x=0
c        if(pp-p(3,lf).gt.0..and.pp+p(3,lf).gt.0.)x=
c     *       0.5*log((pp+p(3,lf))/(pp-p(3,lf)))
         x=abs(x)
      elseif(inom.eq.93)then                                     !'phiexp'
        x=1000
        pt=sqrt(p(1,lf)**2+p(2,lf)**2)
        if(pt.gt.0.)x=sign(1.,p(2,lf))*acos(p(1,lf)/pt)
        if(x.lt.-pi)x=x+2*pi
        if(x.gt.pi)x=x-2*pi
      elseif(inom.eq.94)then  !'epx' ... event plane stuff 
        phi=polar( p(1,lf) , p(2,lf) )
        pt=sqrt(p(2,lf)**2+p(1,lf)**2)
        eta=0
        ivn=ypara(icas+2,n)
        if(p(3,lf).ne.0..and.pt.ne.0.)eta=sign(1.,p(3,lf))*
     *       alog((sqrt(p(3,lf)**2+pt**2)+abs(p(3,lf)))/pt)
        if(eta.gt.0)then
        phi2=ypara(icas,n)   !phi2neg
        else
        phi2=ypara(icas+1,n)   !phi2pos
        endif
        x=cos(ivn*(phi-phi2))/ypara(icas+3,n)
      !---------------------------------
      !The following is obsolete (use mux) !! 
      !elseif(inom.eq.94)then    !'mu24' ... charged ptl mult |eta|<2.4 (CMS)
      !  x=multc24
      !---------------------------------
      elseif(inom.eq.95)then    !'rcastp' ... EM to all energy ratio in Castor (eta>0) (CMS)
        x=rcast(1)
      elseif(inom.eq.96)then    !'rcastt' ... EM to all energy ratio in Castor (eta<0) (CMS)
        x=rcast(2)
      elseif(inom.eq.97)then  ! 'corh'
        x=yyarr( nhisxxx(n) , j )
      elseif(inom.eq.98)then    !'nmass'--------- nuclear mass number
        id=abs(idepos)
        x=0.
        if(id.eq.1220.or.id.eq.1120)then
          x=1.
        elseif(id.eq.17)then
          x=2.
        elseif(id.eq.18)then
          x=3.
        elseif(id.eq.19)then
          x=4.
        elseif(id.gt.1e9)then
          x=mod(id/10,1000)
c          print *,id,x
        endif
      elseif(inom.eq.99)then  ! 'ekin'
        x=ekin
        
c--------------------------------- 101 - 200 ----------------------------

      elseif(inom.eq.101)then           !mulevt
        x=1.
      elseif(inom.eq.102)then                      !'etevt'
        x=0
        if(istptl(j).eq.istxxx)then
         eef=p(4,lf)
         if(maproj.gt.1.or.matarg.gt.1)then
           if(idepos.ge.1000)eef=eef-prom
           if(idepos.le.-1000)eef=eef+prom
         endif
         pp=sqrt(p(1,lf)**2+p(2,lf)**2+p(3,lf)**2)
         if(pp.ne.0.)x=eef*sqrt(p(1,lf)**2+p(2,lf)**2)/pp
         if(x.ne.x)then
           write(ifch,*)x,eef,p(1,lf),p(2,lf),p(3,lf),pp,prom,idepos,j
           call alist('xan&',1,nptl)
           stop 'probleme dans xan'
         endif
        endif
      elseif(inom.eq.103)then
        x=p(4,lf)/1000.
      elseif(inom.eq.104)then                       !'ev6evt'
        x=0
        if(istptl(j).eq.istxxx)then
         pt=sqrt(p(2,lf)**2+p(1,lf)**2)
         eta=0
         if(p(3,lf).ne.0..and.pt.ne.0.)eta=sign(1.,p(3,lf))*
     *   alog((sqrt(p(3,lf)**2+pt**2)+abs(p(3,lf)))/pt)
         if(pt.eq.0.)eta=sign(1e5,p(3,lf))
         if(eta.gt.6.0)then
          eef=p(4,lf)
          if(idepos.ge.1000)eef=eef-prom
          if(idepos.le.-1000)eef=eef+prom
          pp=sqrt(p(1,lf)**2+p(2,lf)**2+p(3,lf)**2)
          if(pp.ne.0.)x=0.001*eef
         endif
        endif
      elseif(inom.eq.105)then
        !etot=maproj*emax(lf)+matarg*0.94  !nur richtig fur target frame!!!!!
        etot=maproj*emax(lf)+matarg*emax(lf)
        x=p(4,lf)/etot
      elseif(inom.eq.106)then
        x=isign(1,idepos)
      elseif(inom.eq.107)then                       !'ptevt'
        x=sqrt(p(2,lf)**2+p(1,lf)**2)
      elseif(inom.eq.108)then                       !'pmxevt'
        x=pmxevt

c--------------------------------- > 200 ----------------------------

      elseif(inom.eq.201)then
        x=1.
      elseif(inom.eq.202)then
        x=egyevt
      elseif(inom.eq.203)then
        x=bimevt
      elseif(inom.eq.204)then                       !'xbjevt'
        x=xbjevt
      elseif(inom.eq.205)then                       !'qsqevt'
        x=qsqevt
      elseif(inom.eq.206)then                       !'yevt'
        x=qsqevt/xbjevt/ecms**2
      elseif(inom.eq.207)then                       !'eloevt'
c        x=qsqevt/4./elepti+elepti*(1.-qsqevt/xbjevt/ecms**2)
        x=elepto
      elseif(inom.eq.208)then                       !'nd1evt' 'nd-...'
        select case (icas)
        case(1)
          x=nsdiff(istuse(n),1,noweak(n))
        case(17)
          x=nsdiff(istuse(n),17,noweak(n))   !'nd-atlas-pt100-evt'
        case(18)
          x=nsdiff(istuse(n),18,noweak(n))   !'nd-atlas-pt500eta25-evt'
        case(19)
          x=nsdiff(istuse(n),19,noweak(n))   !'nd-atlas-pt500eta08-evt'
        case default
          stop'ERROR 211112' 
        end select
      elseif(inom.eq.209)then                       !'nd2evt'
        x=nsdiff(istuse(n),2,noweak(n))
      elseif(inom.eq.210)then                       !'theevt'
c        eloevt=qsqevt/4./elepti+elepti*(1.-qsqevt/xbjevt/ecms**2)
c        x=acos(1-qsqevt/2./elepti/eloevt)/pi*180.
        x=acos(1-qsqevt/2./elepti/elepto)/pi*180.
      elseif(inom.eq.211)then                       !'nspevt'
        x=0
        do i=1,nptl
         call getistptl(i,istx)
         call getidptl(i,id)
         if((istx.eq.30.or.istx.eq.31)
     &      .and.int(id/1000000).eq.1)x=x+1  !1111111111
        enddo
      elseif(inom.eq.212)then                       !'nhpevt'
        x=0
        do i=1,nptl
         call getistptl(i,istx)
         call getidptl(i,id)
         call getpptl(i,p1,p2,p3,p4,p5)
         if((istx.eq.30.or.istx.eq.31)
     &      .and.int(id/1000000).eq.3)x=x+1  !222222222
        enddo
      elseif(inom.eq.213)then                       !'sigtot'
        x=sigtot
      elseif(inom.eq.214)then                       !'sigela'
        x=sigela
      elseif(inom.eq.215)then                       !'sloela'
        x=sloela
      elseif(inom.eq.216)then                       !'nrgevt'
        x=0
        do i=1,nptl
         call getistptl(i,istx)
         call getidptl(i,id)
         call getpptl(i,p1,p2,p3,p4,p5)
          if(istx.eq.31.and.int(id/10000).eq.2)x=x+1
        enddo
      elseif(inom.eq.217)then                       !qevt
        x=sqrt(qsqevt)
      elseif(inom.eq.218)then   !qevt
        if(iappl.eq.8)then
          x=qtl
        else
          x=pprt(1,5)
        endif
      elseif(inom.eq.219)then                       !'nd0evt'  UA1
        x=nsdiff(istuse(n),0,noweak(n))
      elseif(inom.eq.220)then!------------------------------------------
        x=sngl(bofra(1,lf))     !thrust
      elseif(inom.eq.221)then
        x=1.-sngl(bofra(1,lf))  !1-thrust
      elseif(inom.eq.222)then
        x=sngl(bofra(2,lf))     !major
      elseif(inom.eq.223)then
        x=sngl(bofra(3,lf))     !minor
      elseif(inom.eq.224)then
        x=sngl(bofra(2,lf)-bofra(3,lf)) !oblateness
      elseif(inom.eq.230)then!------------------------------------------
        x=1.5*(1.-sngl(bofra(1,lf))) !spherecity
      elseif(inom.eq.231)then
        x=1.5*sngl(bofra(3,lf)) !aplanarity
      elseif(inom.eq.232)then
        x=3.*sngl(bofra(1,lf)*bofra(2,lf)+bofra(1,lf)*bofra(3,lf)
     &       +bofra(2,lf)*bofra(3,lf)) !c-parameter
      elseif(inom.eq.233)then
        x=27.*sngl(bofra(1,lf)*bofra(2,lf)*bofra(3,lf))   !d-parameter
      elseif(inom.eq.234)then                       !npoevt
        x=0
        do i=1,nptl
         call getistptl(i,istx)
         call getidptl(i,id)
         call getpptl(i,p1,p2,p3,p4,p5)
         if(istx.eq.30.or.istx.eq.31)x=x+1
        enddo
      elseif(inom.eq.235)then                       !'npnevt'
        x=ng1evt   !npjevt+ntgevt   !npnevt
      elseif(inom.eq.236)then                       !'ikoevt'
        x=ikoevt
c      elseif(inom.eq.237)then                       !'iktevt'
c        x=zkotest
      elseif(inom.eq.238)then  !npxevt ... nr of pomerons, including absorbed
        x=0
        do i=1,nptl
         call getistptl(i,istx)
         call getidptl(i,id)
         call getpptl(i,p1,p2,p3,p4,p5)
         if(istx.eq.30.or.istx.eq.31)x=x+1
         if(mod(abs(id),100).eq.94)x=x+0.5
        enddo
      elseif(inom.eq.239)then                       !'nd6evt'
        x=nsdiff(istuse(n),6,noweak(n))
      elseif(inom.eq.240)then    !mu1evt ... charged ptl multipl for central rap
        x=multy1
      elseif(inom.eq.241)then    !muievt ... charged ptl multipl
        x=multyi
      elseif(inom.eq.242)then    !hgtevt ... charged ptl multipl for central eta
        x=multc05
      elseif(inom.eq.243)then                       !difevt
        npom=0
        do i=1,nptl
         call getistptl(i,istx)
         call getidptl(i,id)
         call getpptl(i,p1,p2,p3,p4,p5)
         if(istx.eq.30.or.istx.eq.31)npom=npom+1
        enddo
        x=0
        if(npom.eq.0)x=1
      elseif(inom.eq.244)then                       !dixevt
        zpom=0
        do i=1,nptl
         call getistptl(i,istx)
         call getidptl(i,id)
         call getpptl(i,p1,p2,p3,p4,p5)
         if(istx.eq.30.or.istx.eq.31)zpom=zpom+1
         if(mod(abs(id),100).eq.94)zpom=zpom+0.5
        enddo
        x=0
        if(abs(zpom).lt.0.001)x=1
      elseif(inom.eq.245)then                       !'nd7evt' CMS NSD
        x=nsdiff(istuse(n),7,noweak(n))
      elseif(inom.eq.246)then                       !'nd8evt'  ATLAS
        x=nsdiff(istuse(n),8,noweak(n))
      elseif(inom.eq.247)then                       !'nd9evt'  ALICE 900 GeV
        x=nsdiff(istuse(n),9,noweak(n))
      elseif(inom.eq.248)then                       !'ndaevt'  ALICE 2.36 TeV
        x=nsdiff(istuse(n),10,noweak(n))
      elseif(inom.eq.249)then                       !'ndbevt'  ALICE Inel > 0
        x=nsdiff(istuse(n),11,noweak(n))
      elseif(inom.eq.250)then
        if(iappl.eq.8)then      !mass in
          x=-pptl(5,6)
        else
          x=pprt(5,3)
        endif
      elseif(inom.eq.251)then
        if(iappl.eq.8)then      !mass out
          x=pptl(5,7)
        else
          x=pprt(5,2)
        endif
      elseif(inom.eq.252)then
        if(iappl.eq.8)then
          x=-pptl(4,6)
        else
          x=pprt(4,2)
        endif
      elseif(inom.eq.253)then
        if(iappl.eq.8)then
          x=pptl(4,7)
        else
          x=pprt(4,3)
        endif
      elseif(inom.eq.254)then
        if(iappl.eq.8)then
          x=abs(pptl(3,6))
        else
          x=abs(pprt(3,2))
        endif
      elseif(inom.eq.255)then
        if(iappl.eq.8)then
          x=abs(pptl(3,7))
          do l=1,5
            p(l,lf)=pptl(l,7)
          enddo
          if(imofra(1,lf).ne.0)then
            call utrota(imofra(1,lf),r1fra(1,lf),r1fra(2,lf),r1fra(3,lf)
     $           ,p(1,lf),p(2,lf),p(3,lf))
          endif
          if(imofra(2,lf).ne.0)then !the x-z exchanged is ok !!
            call utrota(imofra(2,lf),r2fra(3,lf),r2fra(2,lf),r2fra(1,lf)
     $           ,p(3,lf),p(2,lf),p(1,lf))
          endif
          if(imofra(3,lf).ne.0)then
            call utlob4(imofra(3,lf),bofra(1,lf),bofra(2,lf)
     $           ,bofra(3,lf) ,bofra(4,lf),bofra(5,lf)
     $           ,p(1,lf),p(2,lf),p(3,lf),p(4,lf))
          endif
          x=abs(p(3,lf))
        else
          x=abs(pprt(3,3))
        endif
      elseif(inom.eq.256)then  !pxfevt: leading proton xf in cms
        x=-2
c       pmax=sqrt((ecms/2.)**2-prom**2)
          pmax=pnullx               !???????????????????
        do i=1,nptl
          call getistptl(i,istx)
          call getidptl(i,id)
          call getpptl(i,p1,p2,p3,p4,p5)
          if(id.eq.1120.and.istx.eq.istxxx)then
            if(iframe.eq.11)then
              pz=p3
            else
              amt=sqrt(prom**2+p1**2+p2**2)
              rap=alog((p4+p3)/amt)
     &           -alog((sqrt(pnll**2+ecms**2)+pnll)/ecms)
              pz=amt*sinh(rap)
            endif
            x=max(x,pz/pmax)
          endif
        enddo
      elseif(inom.eq.257)then  !  pi+xf: pi+ yield at cms xf>0.01
        x=0.
c        pmax=sqrt((ecms/2)**2-prom*2)
          pmax=pnullx               !???????????????????
        do i=1,nptl
          call getistptl(i,istx)
          call getidptl(i,id)
          call getpptl(i,p1,p2,p3,p4,p5)
          if(id.eq.120.and.istx.eq.istxxx)then
            if(iframe.eq.11)then
              pz=p3
            else
              amt=sqrt(p5**2+p1**2+p2**2)
              rap=alog((p4+p3)/amt)
     &           -alog((sqrt(pnll**2+ecms**2)+pnll)/ecms)
              pz=amt*sinh(rap)
            endif
            if(pz/pmax.gt.0.01)x=x+1.
          endif
        enddo
      elseif(inom.eq.258)then  !  pi-xf: pi- yield at cms xf>0.01
        x=0.
c        pmax=sqrt((ecms/2)**2-prom*2)
          pmax=pnullx               !???????????????????
        do i=1,nptl
          call getistptl(i,istx)
          call getidptl(i,id)
          call getpptl(i,p1,p2,p3,p4,p5)
          if(id.eq.-120.and.istx.eq.istxxx)then
            if(iframe.eq.11)then
              pz=p3
            else
              amt=sqrt(p5**2+p1**2+p2**2)
              rap=alog((p4+p3)/amt)
     &           -alog((sqrt(pnll**2+ecms**2)+pnll)/ecms)
              pz=amt*sinh(rap)
            endif
            if(pz/pmax.gt.0.01)x=x+1.
          endif
        enddo
      elseif(inom.eq.260)then!------------------------------
        x=sigcut
      elseif(inom.eq.261)then
        x=keu
      elseif(inom.eq.262)then
        x=ked
      elseif(inom.eq.263)then
        x=kes
      elseif(inom.eq.265)then
        x=kolevt
      elseif(inom.eq.266)then
        x=sigsd
      elseif(inom.eq.267)then
        x=nglevt
      elseif(inom.eq.268)then  ! kppevt : collision number per participant
        x=kolevt/float(npjevt+ntgevt)
      elseif(inom.eq.269)then  ! npievt : pion + multi per event
        x=0
        do i=1,nptl
         call getistptl(i,istx)
         call getidptl(i,id)
         call getpptl(i,p1,p2,p3,p4,p5)
         if(id.eq.120)x=x+1
        enddo
      elseif(inom.eq.270)then  ! np2evt : pion + multi per event
        x=0
        do i=1,nptl
         call getistptl(i,istx)
         call getidptl(i,id)
         call getpptl(i,p1,p2,p3,p4,p5)
         if(id.eq.120)x=x+1
        enddo
        x=x/float(npjevt+ntgevt)
      elseif(inom.eq.271)then
        x=sigdif
      elseif(inom.eq.272)then  !number of inelastic collisions per event
        x=koievt
      elseif(inom.eq.273)then  ! inelasticity (energy loss of leading particle)
        x=0.
        do i=maproj+matarg+1,nptl
          call getistptl(i,istx)
          call getidptl(i,id)
          call getpptl(i,p1,p2,p3,p4,p5)
          if(istx.eq.istxxx)then
            if((((abs(id).gt.1000.and.abs(id).lt.10000)
     *           .and.idproj.gt.1000).or.(iabs(id).gt.100
     *           .and.idproj.lt.1000)).and.p4
     *           .gt.x.and.p3.gt.0.)x=p4
          endif
        enddo
        Eini=pptl(4,1)
        if(Eini.gt.0.)x=(Eini-x)/Eini
      elseif(inom.eq.274)then  ! elasticity (energy of leading particle)
        x=0.
        do i=maproj+matarg+1,nptl
          call getistptl(i,istx)
          call getidptl(i,id)
          call getpptl(i,p1,p2,p3,p4,p5)
          if(istx.eq.istxxx)then
            if((((abs(id).gt.1000.and.abs(id).lt.10000)
     *           .and.idproj.gt.1000).or.(iabs(id).gt.100
     *           .and.idproj.lt.1000)).and.p4
     *           .gt.x.and.p3.gt.0.)x=p4
          endif
        enddo
        Eini=pptl(4,1)
        if(Eini.gt.0.)x=x/Eini
      elseif(inom.eq.275)then         !'itgevt'
        x=0
        if(nint(ypara(1,n)).ne.0)x=1
      elseif(inom.eq.276)then         !'hrdevt' ......  1 = hard event
        x=0
        if(nint(ypara(1,n)).ne.0)x=1
      elseif(inom.eq.277)then         !'ej1evt' .... et of jet 1
        x=0
        if(nint(ypara(1,n)).ne.0)
     &  x=ypara(2,n)
      elseif(inom.eq.278)then         !'pj1evt' .... phi of jet 1
        x=1000
        if(nint(ypara(1,n)).ne.0)
     &  x=ypara(4,n)
      elseif(inom.eq.279)then         !'ej2evt' .... et of jet 2
        x=0
        if(nint(ypara(6,n)).ne.0)
     &  x=ypara(7,n)
      elseif(inom.eq.280)then         !'pj2evt' .... delta_phi of jet 2 1
        x=1000
        if(nint(ypara(6,n)).ne.0)then
          x=ypara(9,n)-ypara(4,n)
           if(x.lt.-2.5*pi)then
            x=x+4*pi
           elseif(x.lt.-0.5*pi)then
            x=x+2*pi
          elseif(x.gt.3.5*pi)then
            x=x-4*pi
          elseif(x.gt.1.5*pi)then
            x=x-2*pi
          endif
        endif
      elseif(inom.eq.281)then         !'zppevt'
        x=zppevt
      elseif(inom.eq.282)then         !'zptevt'
        x=zptevt
      elseif(inom.eq.283)then
        stop '**********not used*********'
      elseif(inom.eq.284)then                       !'nd3evt'
        x=nsdiff(istuse(n),3,noweak(n))
      elseif(inom.eq.285)then                       !'nd4evt'
        x=nsdiff(istuse(n),4,noweak(n))
      elseif(inom.eq.286)then         !'mubevt'
        x=multeb
      elseif(inom.eq.287)then                       !'nd5evt'
        x=nsdiff(istuse(n),5,noweak(n))
      elseif(inom.eq.288)then
        x=ekievt
      elseif(inom.eq.289)then                       !'diffmevt'
        x=isdiff(1)
      elseif(inom.eq.290)then                       !'diffxevt'
        x=isdiff(2)
      elseif(inom.eq.291.or.inom.eq.292)then  ! mass of produced system (inelasticity of leading particle )
        x=0.
        i=idlead
        if(i.gt.0)then
          call getistptl(i,istx)
          call getidptl(i,id)
          call getpptl(i,p1,p2,p3,p4,p5)
          pmax=pnullx
          if(mod(ifra(lf),10).eq.2)pmax=pnll
          xxx=dble(abs(pptl(3,i)/pmax))
          xxx=(1.d0-sqrt(xxx))*(1.d0+sqrt(xxx))*dble(engy*engy) !for precision
          if(inom.eq.291.and.xxx.gt.0.d0)then
            x=sngl(sqrt(xxx))
          else
            x=sngl(xxx)
          endif
c         write(ifmt,*)'la',i,idptl(i),x,abs(pptl(3,i)/pmax)
        endif
      elseif(inom.eq.293)then  ! tdevt : -t of leading particle
        x=0.
        i=idlead
        if(i.gt.0)then
          call getistptl(i,istx)
          call getidptl(i,id)
          call getpptl(i,p1,p2,p3,p4,p5)
          pmax=pnullx
          if(mod(ifra(lf),10).eq.2)pmax=pnll
c        xxx=(amproj**2-2.*sqrt(amproj**2+pmax**2)*4
c     *      +2.*abs(pmax*p4)+p5**2)
          ppt=sqrt(p1**2+p2**2)
          if(p3.ne.0.)then
            theta=atan(ppt/p3)
          else
            theta=pi/2.
          endif
          x=abs(p3/pmax)
          if(p3.ne.0.)then
            theta=atan(ppt/p3)
          else
            theta=pi/2.
          endif
c -t definition of UA4 (Phys Lett B136,217)
c          x=p5**2*(1.-x)**2/x+2*x*pmax*pmax*(1.-cos(theta))
c -t definition of TOTEM
          x=(p1**2+p2**2+p3**2)*theta**2
c         write(*,*)'ici',i,idptl(i),theta,x,xxx
        endif
      elseif(inom.eq.294)then          !'ndpevt' pomeron from diffraction
        x=0
        do i=1,nptl
         call getistptl(i,istx)
         call getityptl(i,ity)
         call getidptl(i,id)
         call getpptl(i,p1,p2,p3,p4,p5)
         if((istx.eq.30.or.istx.eq.31)
     &      .and.mod(ity,10).eq.5)x=x+1
        enddo
      elseif(inom.eq.295)then          !'rapgap' rapidity gap
        x=rapgap
      elseif(inom.eq.296)then        !'ng1evt'
        x=ng1evt
      elseif(inom.eq.297)then        !'r21evt'
        x=0
        if(ng1evt.ne.0)x=ng2evt/float(ng1evt)
      elseif(inom.eq.298)then
        x=sval(1,n)
      elseif(inom.eq.299)then
        x=sval(2,n)
      elseif(inom.eq.301)then   !---------------------------------------------
        if(j.eq.0)then          !initialize
          do l=1,4
            aimuni(l,n)=0.
          enddo
        elseif(j.eq.-99)then   !final calculation
          am2=aimuni(4,n)**2-aimuni(3,n)**2
     $         -aimuni(2,n)**2-aimuni(1,n)**2
          x=sign(sqrt(abs(am2)),am2)
c          print *, x
        else                    !routine work
          do l=1,4
            aimuni(l,n)=aimuni(l,n)+p(l,lf)
          enddo
c          print *, j,(p(l,lf),l=1,5)
        endif
      elseif(inom.ge.302.and.inom.le.305)then   !-----------------------
        if(j.eq.0)then          !initialize
          do l=1,4
            aimuni(l,n)=0.
          enddo
        elseif(j.eq.-99)then   !final calculation
          if(inom.eq.302) x=max(aimuni(1,n)/2/(aimuni(2,n)+aimuni(4,n))
     $         ,aimuni(3,n)/2/(aimuni(2,n)+aimuni(4,n)))
          if(inom.eq.303) x=min(aimuni(1,n)/2/(aimuni(2,n)+aimuni(4,n))
     $         ,aimuni(3,n)/2/(aimuni(2,n)+aimuni(4,n)))
          if(inom.eq.304) x=abs(aimuni(1,n)/2/(aimuni(2,n)+aimuni(4,n))
     $         -aimuni(3,n)/2/(aimuni(2,n)+aimuni(4,n)))
          if(inom.eq.305) x=aimuni(1,n)/2/(aimuni(2,n)+aimuni(4,n))
     $         +aimuni(3,n)/2/(aimuni(2,n)+aimuni(4,n))
        else                    !routine work
          l=0
          if(p(3,lf).lt.0.)l=2
          aimuni(1+l,n)=aimuni(1+l,n)+sqrt(p(1,lf)**2+p(2,lf)**2)
          aimuni(2+l,n)=aimuni(2+l,n)
     $         +sqrt(p(1,lf)**2+p(2,lf)**2+p(3,lf)**2)

        endif
      elseif(inom.eq.306.or.inom.eq.307.or.inom.eq.308)then !---------
        if(j.eq.0)then          !initialize
          do ll=1,8
            aimuni(ll,n)=0.
          enddo
        elseif(j.eq.-99)then   !final calculation
          am2a=(aimuni(4,n)-aimuni(3,n))*(aimuni(4,n)+aimuni(3,n))
     $         -aimuni(2,n)**2-aimuni(1,n)**2
          am2b=(aimuni(8,n)-aimuni(7,n))*(aimuni(8,n)+aimuni(7,n))
     $         -aimuni(6,n)**2-aimuni(5,n)**2
          if(inom.eq.306)x=(max(0.,am2a,am2b))/engy**2
          if(inom.eq.307)x=(max(0.,min(am2a,am2b)))/engy**2
          if(inom.eq.308)x=(abs(am2a-am2b))/engy**2
        else                    !routine work
          ll=0
          if(p(3,lf).lt.0.)ll=4
          do l=1,4
            aimuni(l+ll,n)=aimuni(l+ll,n)+p(l,lf)
          enddo
        endif
      elseif (inom.eq.310.or.inom.eq.311) then !---------
        if(j.eq.0)then          !initialize
          aimuni(1,n)=0
          aimuni(2,n)=0
          do i=1,nptl
           call getistptl(i,istx)
           call getidptl(i,id)
           call getpptl(i,p1,p2,p3,p4,p5)
c            charge=0.
             if(istx.eq.istxxx) then
               if (id.eq.idcod(1,n)) aimuni(1,n)=aimuni(1,n)+1.
               if (id.eq.idcod(2,n)) aimuni(2,n)=aimuni(2,n)+1.
             endif
           enddo
        elseif(j.eq.-99)then   !final calculation
          if(aimuni(1,n).eq.0.or.aimuni(2,n).eq.0) then
            ncevt(n)=ncevt(n)-1
          endif
          x=xmin(n)-100.
          do i=1,nbin(n)
            zcbin(i,nac(n),n)=abs(zcbin(i,nac(n),n))
          enddo
        else                    !routine work
          if( istptl(j).eq.istxxx
     $         .and. aimuni(1,n).ne.0. .and. aimuni(2,n).ne.0. ) then
            id1=idepos
            if(id1.eq.idcod(1,n) .or. id1.eq.idcod(2,n)) then
              y1=sign(1.,pptl(3,j))*alog((pptl(4,j)+abs(pptl(3,j)))
     *             /sqrt(pptl(5,j)**2+pptl(1,j)**2+pptl(2,j)**2))

              do iloo=1,nptl
      
                i=iloo   

                if(i.eq.j .or. istptl(i).ne.istxxx) goto 124
                id2=idptl(i)
                if(id2.eq.idcod(1,n) .or. id2.eq.idcod(2,n)) then
                  y2=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))
     *                 /sqrt(pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2))
                  dy=(y2-y1)
                  if(inom.eq.311) dy=abs(dy)
                  ib=1+int((dy-xmin(n))*xinc(n))
                  if(dy.ge.xmin(n).and.dy.le.xmax(n)) then
                    if( id1.eq.idcod(1,n) ) then
                      if( id2.eq.idcod(2,n) ) then
                        bin(ib,nac(n),n)=bin(ib,nac(n),n)+.5/aimuni(2,n)
                        zcbin(ib,nac(n),n)=zcbin(ib,nac(n),n)+1
                      else
                        bin(ib,nac(n),n)=bin(ib,nac(n),n)-.5/aimuni(1,n)
                        zcbin(ib,nac(n),n)=zcbin(ib,nac(n),n)-1
                      endif
                    else        !id1 is idcod(2,n)
                      if(id2.eq.idcod(1,n)) then
                        bin(ib,nac(n),n)=bin(ib,nac(n),n)+.5/aimuni(1,n)
                        zcbin(ib,nac(n),n)=zcbin(ib,nac(n),n)+1
                      else
                        bin(ib,nac(n),n)=bin(ib,nac(n),n)-.5/aimuni(2,n)
                        zcbin(ib,nac(n),n)=zcbin(ib,nac(n),n)-1
                      endif
                    endif
                  endif
                endif
 124            continue 
              
              enddo
          
            endif
          endif
        endif
      elseif (inom.eq.312) then !---------
        x=sigine
      elseif (inom.eq.313) then !---------
        x=sigineaa
      elseif (inom.eq.314) then !--------- 
        x=alpD(idxD0,iclpro,icltar)
      elseif (inom.eq.315) then !---------
        x=alpD(1,iclpro,icltar)
      elseif (inom.eq.316) then !---------) 
        x=0.5*abs(betDp(idxD0,iclpro,icltar)+betDp(idxD0,iclpro,icltar))
      elseif (inom.eq.317) then !--------- 
        x=0.5*abs(betDp(1,iclpro,icltar)+betDpp(1,iclpro,icltar))
      elseif (inom.eq.318) then !---------
        x=rexdif(iclpro)
      elseif (inom.eq.319) then !---------
        x=rexdif(icltar)
      elseif(inom.eq.320)then    !m14evt ... multipl |eta|<1, pt>0.4
        x=multc14
      elseif(inom.eq.321)then    !ht3evt ... height |eta|<3.15
        x=multc3/6.3
      elseif (inom.eq.322) then !---------
        x=sigineex
      elseif (inom.eq.323) then !---------
        x=sigdifex
      elseif (inom.eq.324) then !---------
        x=sigsdex
      elseif (inom.eq.325) then !---------
        x=ekin
      elseif (inom.eq.326) then !---------
        x=sigcutaa
      elseif (inom.eq.327) then !---------
        x=sigtotaa
      elseif (inom.eq.328) then !---------
        x=xkappafit(1,iclegy,iclpro,icltar)
      elseif (inom.eq.329) then !---------
        x=abs(gamD(idxD0,iclpro,icltar))
      elseif (inom.eq.330) then !---------
       x=abs(gamD(1,iclpro,icltar))
      elseif (inom.eq.331) then !---------
       x=delD(idxD0,iclpro,icltar)
      elseif (inom.eq.332) then !---------
        x=delD(1,iclpro,icltar)
      elseif(inom.eq.333)then                       !'nd6evt'
        x=nsdiff(istuse(n),6,noweak(n))
      elseif(inom.eq.334)then !'muxevt' ... multipl 
        x= ypara(icas,n)     
      elseif(inom.eq.335)then                       
        x=abs(typevt)                        !ND(1), DD(2), or SD(3)
      elseif(inom.eq.339)then    !m25evt ... multipl |eta|<2.5, pt>0.5
        x=multc25
      elseif(inom.eq.340)then    !'segevt' ... segment multiplicity
        x=segevt
      elseif(inom.eq.341)then    !'ielevt'
        x=nsdiff(istuse(n),11,noweak(n))
      elseif(inom.eq.342)then                       !'mc1evt' charged particle
        x=multc1                                    ! mult for |eta|<1
      elseif(inom.eq.343)then   !CDF SD trigger 'sdcdf'
        x=isdiff(3)
      elseif(inom.eq.344)then   !CDF DPE trigger 'dpecdf'
        x=isdiff(4)
      elseif(inom.eq.345)then   !CDF DD trigger 'ddcdf'
        x=isdiff(5)
      elseif(inom.eq.346)then   !'phievt'
        x=phievt
      elseif(inom.eq.347)then   !'ndcevt'  CMS hadron level NSD (2011)
        x=nsdiff(istuse(n),12,noweak(n))
      elseif(inom.eq.348)then   !'jetevt' ......  1 = jet event
        x=0
        if(nint(ypara(1,n)).ne.0)x=1
      elseif(inom.eq.349)then   !'epszero (Z for pp)'
        x=epszero
      elseif(inom.eq.350)then  !xsievt:  xsi = (M_X^2/s) 
c (where M_X = sqrt{ (sum E)^2 - (sum vec-p)^2 } with sum E and sum vec-p 
c being resp. the sum of the energies and the sum of the 3-momenta of the 
c generated particles in the event, excluding the proton with the largest 
c laboratory momentum)
        x=xsi
cc xsievt:  xsi = 1-xF_leading
c        x=-2
cc       pmax=sqrt((ecms/2.)**2-prom**2)
c          pmax=pnullx               !???????????????????
c        do i=1,nptl
c          if(abs(idptl(i)).gt.100.and.istptl(i).eq.istxxx)then
c            if(iframe.eq.11)then
c              pz=pptl(3,i)
c            else
c              amt=sqrt(prom**2+pptl(1,i)**2+pptl(2,i)**2)
c              rap=alog((pptl(4,i)+pptl(3,i))/amt)
c     &           -alog((sqrt(pnll**2+ecms**2)+pnll)/ecms)
c              pz=amt*sinh(rap)
c            endif
c            x=max(x,abs(pz/pmax))
c          endif
c        enddo
c        x=max(0.,1.-x)
      elseif(inom.eq.351)then  !xsicms: CMS determination of xsi=M2_X/s
                               !using xsi=Sum [(E+pz)_i/sqrt(s)]      where i runs on every reconstructed particles (= |eta|<4.9 charged + neutral)
        x=-2
c       pmax=sqrt((ecms/2.)**2-prom**2)
        pmax=pnullx               !???????????????????
        Ef=0.
        Eb=0.
        Pf=0.
        Pb=0.
        Esum=0.
        Psum=0.

        do iloo=1,nptl
      
          i=iloo   

          if(istptl(i).eq.istxxx)then
            pt=sqrt(pptl(2,i)**2+pptl(1,i)**2)
            pz=pptl(3,i)
            eta=0.
            if(abs(pz).gt.0..and.pt.gt.0.)eta=sign(1.,pz)*
     *   log((sqrt(pz**2+pt**2)+abs(pz))/pt)
            if(pt.eq.0.)eta=sign(1e5,pz)
            if(eta.ge.4.9)then
              Ef=Ef+pptl(4,i)
              Pf=Pf+pptl(3,i)
            elseif(eta.le.-4.9)then
              Eb=Eb+pptl(4,i)
              Pb=Pb+pptl(3,i)
            else
              Esum=Esum+pptl(4,i)
              Psum=Psum+pptl(3,i)
            endif
c        write(ifch,*)'ici',i,eta,Ef,Eb,Esum,Psum,Ef-Pf,Eb+Pb
          endif

        enddo

        if(Ef.ge.Eb)then
          x=Esum+Psum+Eb+Pb
        else
          x=Esum-Psum+Ef-Pf
        endif
        x=max(0.,x/ecms)
c        write(ifmt,*)'ici',x,min(Esum-Psum,Esum+Psum)/ecms
      elseif(inom.eq.352)then   !'calevt' ......  energy in eta range
        x= ypara(1,n)      
      elseif(inom.eq.353)then   !'fgpevt' ......  max forward rapidity gap in eta range
        x= ypara(2,n)      
      elseif(inom.eq.354)then   !'bgpevt' ......  max backward rapidity gap in eta range
        x= ypara(3,n)      
      elseif(inom.eq.355)then   !'gapevt' ......  max backward rapidity gap in eta range
        x= ypara(4,n)      
      elseif(inom.eq.356)then
        x=sigdd
      elseif(inom.eq.357)then                    !'ajtevt' 
        x=-9999
        if(nint(ypara(1,n)).ne.0)x=ypara(2,n)
      elseif(inom.eq.358)then                    !'fjtevt' 
        x=-9999
        if(nint(ypara(1,n)).ne.0)x=ypara(3,n)
      elseif(inom.eq.359)then                    !'pjtevt' 
        x=-9999
        if(nint(ypara(1,n)).ne.0)x=ypara(4,n)
      elseif(inom.eq.360)then                       
        x=nint(mod(typevt,10.)) !ND(1), DD(2), or SD(3)
      elseif(inom.eq.361)then   !'ndsevt'  CMS hadron level double sided trigger (2012)
        x=nsdiff(istuse(n),13,noweak(n))
      elseif(inom.eq.362)then    !m24evt ... multipl |eta|<2.4 (CMS)
        x=multc24
      elseif(inom.eq.363)then   !'ndhevt'  CMS hadron level single sided trigger (HF 2012)
        x=nsdiff(istuse(n),14,noweak(n))
      elseif(inom.eq.364)then   !'ndfevt'  LHCf trigger
        x=nsdiff(istuse(n),15,noweak(n))
      elseif(inom.eq.365)then    !'rcpevt' ... EM to all energy ratio in Castor (eta>0) (CMS)
        x=rcast(1)
      elseif(inom.eq.366)then    !'rctevt' ... EM to all energy ratio in Castor (eta<0) (CMS)
        x=rcast(2)
      elseif(inom.eq.367)then  !mdievt:  mdiff = M_X^2
c (where M_X = sqrt{ (sum E)^2 - (sum vec-p)^2 } with sum E and sum vec-p 
c being resp. the sum of the energies and the sum of the 3-momenta of the 
c generated particles in the event, excluding the proton with the largest 
c laboratory momentum)
        x=xMdiff
      elseif (inom.eq.368) then !---------
        x=sigddex
      elseif (inom.eq.369) then !---------
        x=0. !facpevt
      elseif(inom.eq.494)then  !weighted spectrum moment for pi+
        x=0. !xmom(1)
      elseif(inom.eq.495)then  !weighted spectrum moment for pi-
        x=0. !xmom(2)
      elseif(inom.eq.496)then  !weighted spectrum moment for K+
        x=0. !xmom(3)
      elseif(inom.eq.497)then  !weighted spectrum moment for K-
        x=0. !xmom(4)
      elseif (inom.eq.498) then !--------- maximum fragment mass projectile
        x=mfragmaxp
      elseif (inom.eq.499) then !--------- maximum fragment mass target
        x=mfragmaxt
      elseif (inom.eq.501) then !---------
        x=sval(1,1)
      elseif (inom.eq.502) then !---------
        x=sval(1,2)
      elseif (inom.eq.503) then !---------
        x=sval(1,3)
      elseif (inom.eq.504) then !---------
        x=sval(1,4)
      elseif(inom.eq.505)then    !'eglevt'  eccentricity
        x=eglevt
      elseif(inom.eq.506)then    !'fglevt'  eccentricity_part
        x=fglevt
      elseif(inom.eq.507)then    !'rglevt'  ratio ng2 / ng1
        x=0
        if(ng1evt.ne.0)
     .  x=ng2evt/float(ng1evt)
      elseif(inom.eq.508)then    !'sglevt'  area S
        x=sglevt
      elseif(inom.eq.509)then                       !'ptrevt'
        x=0
        do i=maproj+matarg+1,minfra
        call getistptl(i,istx)
        call getidptl(i,id)
        call getpptl(i,p1,p2,p3,p4,p5)
        if(istx.eq.25)then
          pt=sqrt(p1**2+p2**2)
          x=max(x,pt)
        endif
        enddo
      elseif(inom.eq.510)then                       !'rr2evt'
        x=cos(2*(phi2neg-phi2pos))
        !write(ifmt,*)'+++++ rr2evt +++++ ',x,phi2neg,phi2pos
      elseif(inom.eq.511)then  !'perevt'
        call getJKNcentr
        x=(ncentr-0.5)*5
      elseif(inom.eq.512)then  !'paievt'
        kkk=icas/10000
        kk1=kkk/10  /8   +1
        kk2=mod(kkk,10)
        kkk=(kk1-1)*(kkkmax/2)+kk2
        x=paievt(kkk)
      elseif(inom.eq.513)then  !'co2evt'
        kkk=icas/10000
        kk1=kkk/10  /8   +1
        kk2=mod(kkk,10)
        kkk=(kk1-1)*(kkkmax/2)+kk2
        x=co2evt(kkk)     
      elseif(inom.eq.514)then  !'co3evt'
        kkk=icas/10000
        kk1=kkk/10  /8   +1
        kk2=mod(kkk,10)
        kkk=(kk1-1)*(kkkmax/2)+kk2
        x=co3evt(kkk)
      elseif(inom.eq.515)then  !'nh1evt'
        x=0
        do i=maproj+matarg+1,minfra
        call getistptl(i,istx)
        if(istx.eq.25)then
          call getpptl(i,p1,p2,p3,p4,p5)
          pt=sqrt(p1**2+p2**2)
          if(pt.gt.10)x=x+1
        endif
        enddo
      elseif(inom.eq.516)then  !'nh2evt'
        x=0
        do i=maproj+matarg+1,minfra
        call getistptl(i,istx)
        if(istx.eq.25)then
          call getpptl(i,p1,p2,p3,p4,p5)
          pt=sqrt(p1**2+p2**2)
          if(pt.gt.20)x=x+1
        endif
        enddo
      elseif(inom.eq.517)then  !'nh3evt'
        x=0
        do i=maproj+matarg+1,minfra
        call getistptl(i,istx)
        if(istx.eq.25)then
          call getpptl(i,p1,p2,p3,p4,p5)
          pt=sqrt(p1**2+p1**2)
          if(pt.gt.50)x=x+1
        endif
        enddo
      elseif(inom.eq.518)then  !'gbyevt'
        x=gbyjmax
      elseif(inom.eq.519)then  !'icyevt'
        call getncenthy
        if(ncenthy.ne.20)stop'##### ERROR 28072011d #####'
        f=(1-(icentrality-0.5)*0.05)
        x=f * 100
      elseif(inom.eq.520)then  !'xp9evt'
        x=xpara(9,n)  !##### do not use xpara in this routine #####
      elseif(inom.eq.521)then  !'hlxevt'
        xi=xpara(9,n)  !##### do not use xpara in this routine #####
        x=helix(xi)
      elseif(inom.eq.522)then  !'icvevt' percentile ALICE
        x=ypara(1,n)
      elseif(inom.eq.523)then  !'corhevt' !corhtr=old
        x=ypara(1,n)
      elseif(inom.eq.524)then  !bg 'phoevt'
        x=1  !dummy
      !-----------------------------------------------------
      elseif(inom.eq.525)then   ! vnsp variables
      !-----------------------------------------------------
        select case (icas)
        case(1)  
          x = x2aevt*x2cevt + y2aevt*y2cevt  !'qacevt'  !kept for the moment
        case(122) 
          m=2 
          x=qvnsp(1,m,1)*qvnsp(2,m,1)+qvnsp(1,m,2)*qvnsp(2,m,2) ! A B 2
        case(123)
          m=3
          x=qvnsp(1,m,1)*qvnsp(2,m,1)+qvnsp(1,m,2)*qvnsp(2,m,2) ! A B 3
        case(124) 
          m=4
          x=qvnsp(1,m,1)*qvnsp(2,m,1)+qvnsp(1,m,2)*qvnsp(2,m,2) ! A B 4
        case(312)
          m=2 
          x=qvnsp(3,m,1)*qvnsp(1,m,1)+qvnsp(3,m,2)*qvnsp(1,m,2) ! C A 2 
        case(313)
          m=3
          x=qvnsp(3,m,1)*qvnsp(1,m,1)+qvnsp(3,m,2)*qvnsp(1,m,2) ! C A 3
        case(314)
          m=4 
          x=qvnsp(3,m,1)*qvnsp(1,m,1)+qvnsp(3,m,2)*qvnsp(1,m,2) ! C A 4
        case(322)
          m=2 
          x=qvnsp(3,m,1)*qvnsp(2,m,1)+qvnsp(3,m,2)*qvnsp(2,m,2) ! C B 2
        case(323)
          m=3
          x=qvnsp(3,m,1)*qvnsp(2,m,1)+qvnsp(3,m,2)*qvnsp(2,m,2) ! C B 3
        case(324)
          m=4 
          x=qvnsp(3,m,1)*qvnsp(2,m,1)+qvnsp(3,m,2)*qvnsp(2,m,2) ! C B 4
        case default
          stop'ERROR 19092020c' 
        end select
      !------------------------------------------------
      elseif(inom.eq.526)then    ! net-cumlt variables
      !------------------------------------------------ 
        x=ypara(ishift+1,n) 
      !-----------------------------------------------------
      !-----------------------------------------------------
      elseif(inom.eq.530)then !'ep2evt''ep3evt''ep4evt' 'ep5evt''rrrevt'
        select case (icas)
        case(1)
          nord=2      
        case(2)
          nord=3
        case(3)
          nord=4
        case(4)
          nord=5
        case(5)
          nord=0
        case default
          stop'ERROR 19092020d' 
        end select
        rrr=0
        avecos=0
        avesin=0
        do i=1,nptl
          call getistptl(i,istx)
          call getpptl(i,p1,p2,p3,p4,p5)
          call getxorptl(i,x1,x2,x3,x4)
          if(istx.eq.7)then
            amt=p5**2+p1**2+p2**2
            if(amt.gt.0..and.p4+abs(p3).gt.0.)then  
              amt=sqrt(amt)
              rap=sign(1.,p3)*log((p4+abs(p3))/amt) 
            else
              rap=1e10                  !
            endif
            if(abs(rap).le.2.5)then
              rr=x1**2+x2**2
              phi=polar(x1,x2)
              rrr=rrr+rr
              avecos=avecos+rr*cos(nord*phi)
              avesin=avesin+rr*sin(nord*phi)
            endif
          endif
        enddo
        if(rrr.gt.0.)then
          avecos=avecos/rrr
          avesin=avesin/rrr
        endif
        if(nord.eq.0)then
          x = rrr
        else
          x = sqrt(avecos**2+avesin**2)
         endif
      !-----------------------------------------------------
      !-----------------------------------------------------
      elseif(inom.eq.531)then  ! 'nu1evt'
        x = 1
      elseif(inom.eq.532)then  ! 'nu2evt'
        x = 2
      elseif(inom.eq.533)then  ! 'nu3evt'
        x = 3
      elseif(inom.eq.534)then  !bg 'phtevt' trigger for photon/hadron correlations
        x=ypara(6,n)
      elseif(inom.eq.535)then ! 'mpoevt' average Pomeron mass
        x=0
        z=0 
        do i=1,nptl
         call getistptl(i,istx)
         call getidptl(i,id)
         call getpptl(i,p1,p2,p3,p4,p5)
         if(istx.eq.30.or.istx.eq.31)then
           z=z+1
           x=x+p5
         endif
        enddo
        if(z.gt.0.)x=x/z
      elseif(inom.eq.536)then ! 'ap2evt' average pt squared particles
        x=0
        z=0 
        do i=1,nptl
         call getistptl(i,istx)
         call getidptl(i,id)
         call getpptl(i,p1,p2,p3,p4,p5)
         if(istx.eq.0)then
           z=z+1
           x=x+p1**2+p2**2
         endif
        enddo
        if(z.gt.0.)x=x/z
      elseif(inom.eq.537)then ! 'aq2evt' Pomeron hardness (mt2)
        x=0
        z=0 
        ii=1
        do while ( ii.lt.nptl )
          call getistptl(ii,istii)
          call getityptl(ii,ityii)
          if(istii.eq.29.and.ityii.ge.20.and.ityii.le.39)then
            call getpptl(ii,p1ii,p2ii,p3ii,p4ii,p5ii)
            p4diff=p4ii
            i=ii+1
            call getiorptl(i,ior)
            do while ( ior.eq.ii )
              call getpptl(i,p1,p2,p3,p4,p5)
              z=z+1
              x=x+p1**2+p2**2+p5**2
              p4diff=p4diff-p4
              i=i+1 
              call getiorptl(i,ior)
            enddo
            !if(p4diff.gt.1.0.and.p4diff.gt.0.05*p4ii)write(ifch,*)
            !.       'WARNING xan p4diff :',ii,nint(p4diff/p4ii*100) , p4ii
            ii=i-1
          endif
          ii=ii+1
        enddo
        if(z.gt.0.)x=x/z
      elseif(inom.eq.538)then ! 'ar2evt' Pomeron hardness (pt2)
        x=0
        z=0 
        ii=1
        do while ( ii.lt.nptl )
          call getistptl(ii,istii)
          call getityptl(ii,ityii)
          if(istii.eq.29.and.ityii.ge.20.and.ityii.le.39)then
            i=ii+1
            call getiorptl(i,ior)
            do while ( ior.eq.ii )
              call getpptl(i,p1,p2,p3,p4,p5)
              z=z+1
              x=x+p1**2+p2**2
              i=i+1 
              call getiorptl(i,ior)
            enddo
            ii=i-1
          endif
          ii=ii+1
        enddo
        if(z.gt.0.)x=x/z
      elseif(inom.eq.539)then                       !'corevt'
        x=0
        do i=minfra,maxfra
        call getistptl(i,istx)
        if(istx.ge.5.and.istx.le.7)x=x+1.
        enddo
      elseif(inom.eq.540)then                       !'epxevt'
        phi2n=ypara(icas,n)  
        phi2p=ypara(icas+1,n)  
        ivn=ypara(icas+2,n)
        x=cos(ivn*(phi2n-phi2p))
      elseif(inom.eq.541)then                       !'ue0evt'
        x=underly(noweak(n),istuse(n),0)
      elseif(inom.eq.542)then                       !'ue1evt'
        x=underly(noweak(n),istuse(n),1)
      elseif(inom.eq.543)then                       !'ue2evt'
        x=underly(noweak(n),istuse(n),2)
      elseif(inom.eq.544)then                       !'ue3evt'
        x=underly(noweak(n),istuse(n),3)
      elseif(inom.eq.545)then                       !'ue4evt'
        x=underly(noweak(n),istuse(n),4)
      elseif(inom.eq.546)then                       !'ue5evt'
        x=underly(noweak(n),istuse(n),5)
      elseif(inom.eq.547)then                       !'ue6evt'
        x=underly(noweak(n),istuse(n),6)
      elseif(inom.eq.548)then                       !'mupevt'
        x=ypara(icas,n)
      elseif(inom.eq.549)then                       !'nddevt'  ALICE inel
        x=nsdiff(istuse(n),16,noweak(n))
      elseif(inom.eq.551)then                       ! ...zevt'  
        if(icas.lt.0)then
          x=-icas
        else
          select case (icas)
          case(1)
            x=nglevt                               !z1zevt = Ncoll
          case(2)
            x=ng1evt                               !z2zevt = Npart
          case(3)
            x=0
            do i=1,nptl
             call getistptl(i,istx)
             if(istx.eq.30.or.istx.eq.31)x=x+1     !z3zevt = Npom
            enddo
          case(4)
            x=ikhevt                 !z4zevt = Nhpom
          case(10)
            x=0.                !z10zevt = Conn1
          case(11)
            x=0.                !z11zevt = Conn2
          case(20)
            x=0.                !z20zevt = BinFac
          case(21)
            x=0.                !z21zevt = BinFacFit
          case default   
            stop'ERROR 210115'
          end select
        endif
      endif 
      end

c----------------------------------------------------------------------
      subroutine epx(n,ishift)  !'epx' 
c----------------------------------------------------------------------
      ! input
      !   n = histogram number
      !   xpara(1,n)-xpara(4,n) ... etamin to etamax with gap
      !   xpara(5,n)-xpara(6,n) ... pmin to pmax 
      !   xpara(7,n) rapidity max
      !   xpara(8,n) choice of weight wi
      !   xpara(9,n) order of anisotropy v_n
      !   xpara(10,n) order of event angle psi_k
      !   xpara(11,n) optimization variable 1 to calc, 0 just result
      !   xpara(12,n) EP resolution correction used (>0) or not (0)
      !                   1=Phobos04, 2=Star15, etc (defined in icl file)
      !   xpara(13,n) only core (1) or all (0) 
      ! output
      !   ypara(1,n) ... phi2neg, 
      !   ypara(2,n) ... phi2pos
      !   ypara(3,n) ... order of anisotropy v_n
      !   ypara(4,n) ... EP resolution correction
      !--------------------------------------------------------------

      include "epos.inc"
      include "epos.incxan"
      common/cepx/avcosnn,avsinnn,avcosnp,avsinnp
      common/cepx2/iwi,ivn,kpsi

      icalc=nint(xpara(ishift+11,n))
      iepreso=nint(xpara(ishift+12,n))
      
      if(icalc.eq.1)then

      eta1=xpara(ishift+1,n)
      eta2=xpara(ishift+2,n)
      eta3=xpara(ishift+3,n)
      eta4=xpara(ishift+4,n)
      pmin=xpara(ishift+5,n)
      pmax=xpara(ishift+6,n)
      rapmax=xpara(ishift+7,n)
      iwi=nint(xpara(ishift+8,n))
      ivn=nint(xpara(ishift+9,n))
      kpsi=nint(xpara(ishift+10,n))
      icore=nint(xpara(ishift+13,n))

      avcosnp=0                 
      avsinnp=0
      avcosnn=0
      avsinnn=0

      do iloo=maproj+matarg+1,nptl
     
        i=iloo   
      
        istxxx=0
c        istx=-2
        if(istu1(istuse(n)).gt.-999)istxxx=istu1(istuse(n))

        if(istptl(i).eq.istxxx)then

         igotype=1
         if(icore.eq.1.and.ityptl(i).ne.60)igotype=0

         if(igotype.eq.1)then

          ida=abs(idptl(i))
          amt=pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2
          px=pptl(1,i)
          py=pptl(2,i)
          pz=pptl(3,i)  
          p4=pptl(4,i)  
          pt=pptl(1,i)**2+pptl(2,i)**2
          pp=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(3,i)**2)
          et=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(5,i)**2)

          if(pt.gt.pmin**2.and.pp.lt.pmax**2)then    
          
            if(amt.gt.0..and.pptl(4,i).gt.0.)then
             amt=sqrt(amt)
             rap=sign(1.,pz)*alog((p4+abs(pz))/amt)
            else
             rap=1000.
            endif
            
            arap=abs(rap)
            if(arap.le.rapmax)then 
            
              if(pt.gt.0.)then
                pt=sqrt(pt)
                p3=pptl(3,i)
                if(p3.ne.0.)then
                  theta=atan(pt/p3)
                else
                  theta=pi/2.
                endif
                if(theta.lt.0.)theta=theta+pi
                eta=sign(1.,pz)*alog((pp+abs(pz))/pt)
              else
                eta=1000.
              endif

              if(iwi.eq.1)then
                if(pt.gt.0.and.pt.lt.2)then !Gev/c
                  weight=pt
                else
                  weight=2
                endif
              elseif(iwi.eq.2)then
                weight=1
              elseif(iwi.eq.3)then
                weight=et
              else
                stop'##### ERROR 18012017 #####' 
              endif

              if(eta.gt.eta3.and.eta.lt.eta4)then
                a=polar(px,py)
                avcosnp=avcosnp+weight*cos(kpsi*a)
                avsinnp=avsinnp+weight*sin(kpsi*a)
              elseif(eta.lt.eta2.and.eta.gt.eta1)then
                a=polar(px,py)
                avcosnn=avcosnn+weight*cos(kpsi*a)
                avsinnn=avsinnn+weight*sin(kpsi*a)
              endif

            endif
 
          endif ! pt...       

         endif ! igotype

        endif ! istptl...

      enddo ! iloo

      endif ! icalc
      
      ypara(ishift+1,n)=polar(avcosnn,avsinnn)/kpsi
      ypara(ishift+2,n)=polar(avcosnp,avsinnp)/kpsi
      ypara(ishift+3,n)=ivn
      ypara(ishift+4,n)=1.0
      if(iepreso.gt.0)then
        do k=1,20
          if(bimevt.le.zlimit(k*5))goto 1
        enddo
  1     k=min(k,20)
        ypara(ishift+4,n)=epreso(iepreso,k)
        if(ypara(ishift+4,n).eq.0.)then
          write(ifmt,*)'ERROR zero in epreso table'
          stop
        endif
      endif
      end

c----------------------------------------------------------------------
      subroutine mux(n,ishift)  !'mux' 'muxevt'
c----------------------------------------------------------------------
      ! input
      !   n = histogram number
      !   xpara(1,n) ... etamin (rapidity for D's)
      !   xpara(2,n) ... etamax (   "    )
      !   xpara(3,n) ... ptmin
      !   xpara(4,n) ... ptmax
      !   xpara(5,n) ... factor    ! 9 ignore 
      !   xpara(6,n) ... divisor   ! 9 istuse
      !   xpara(7,n) ... absolute value of eta (1)
      !   xpara(8,n) ... yboost if >< 0.0001 (boost back)
      !   xpara(9,n) ... etamin2 if >< 0.
      !   xpara(10,n) ... etamax2 if >< 0.
      !   xpara(11,n) ... x0 = normal
      !                   x2 = percentile  
      !                   0x = charged multiplicity
      !                   1x = E (all final hadrons)
      !                   2x = Et (charged hadrons)
      !                   3x = D meson multiplicity
      !                   4x = all particle multiplicity
      !   xpara(12-22,n) ... percentile borders
      ! output
      !   ypara(1or2,n) ... event activity
      !--------------------------------------------------------------
      include "epos.inc"
      include "epos.incxan"
      logical go,LongLivPtl,NoLongLivParent

      sum=0
      io=nint( xpara(ishift+11,n) )
      io1=io/10
      io2=mod(io,10)

      do iloo=1,nptl   !maproj+matarg+1,nptl
      
        i=iloo   

        go=.false.
        n5=nint( xpara(ishift+5,n) )    
        n6=nint( xpara(ishift+6,n) ) 
        n56=n5*10+n6
        call setIstXY(n56,istuse(n),istxxx,istyyy)
        if(io1.eq.1.and.istptl(i).eq.istxxx)then
          go=.true.
        elseif(io1.eq.3.and.
     .   istptl(i).eq.istxxx.or.istptl(i).eq.istyyy)then
          go=.true.
        else 
          if(istptl(i).eq.istxxx)go=.true.
          if(noweak(n).eq.1)then
            if(LongLivPtl(istyyy,i))go=.true.
          elseif(noweak(n).eq.3)then
            if(LongLivPtl(istxxx,i))go=.false.
          endif
          if(go)go=NoLongLivParent(noweak(n),i)
        endif 

        if(go)then
          id=ideposf( 4 ,i)
          ida=abs(id)
          pt2=pptl(1,i)**2+pptl(2,i)**2
          amt2=pt2+pptl(5,i)**2
          yboo=xpara(ishift+8,n)
          p4=pptl(4,i)
          p3=pptl(3,i)
          if(abs(yboo).gt.0.0001)then
          p4= cosh(yboo) * pptl(4,i) + sinh(yboo) *  pptl(3,i)
          p3= sinh(yboo) * pptl(4,i) + cosh(yboo) *  pptl(3,i)
          endif
          pp=sqrt(pt2+p3**2)
          if(pt2.gt.0.)then
            pt=sqrt(pt2)
            eta=sign(1.,p3)*alog((pp+abs(p3))/pt)
          else
            pt=0.
            eta=1000.
          endif
          if(amt2.gt.0.)then
            amt=sqrt(amt2)
            rap=sign(1.,p3)*alog((p4+abs(p3))/amt) 
          else
            amt=0.
            rap=1000.
          endif
          if(ida.eq.240.or.ida.eq.140.or.ida.eq.241)eta=rap
          !if(ida.eq.240.or.ida.eq.140.or.ida.eq.241)
          !.    print*,'XAN mux ',ida,eta,pt
          if(nint(xpara(ishift+7,n)).eq.1.)eta=abs(eta)
          call idflav(id,i1,i2,i3,jdu,i)
          if(i2.ne.0.and.i3.ne.0)then !hadrons
            call idchrg( 12 ,id,ch)
            igo=0
            if(io1.eq.0.and.abs(ch).gt.0.1.or.io1.eq.2)igo=1
            if(io1.eq.3
     .      .and.(ida.eq.240.or.ida.eq.140.or.ida.eq.241))igo=1
            if(io1.eq.4)igo=1
            if(io1.eq.1.and.i.gt.maproj+matarg+1)igo=1
            if(igo.eq.1)then
              val=1
              if(io1.eq.1)then !E 
                val=p4
              elseif(io1.eq.2)then !Et 
                eef=p4
                if(maproj.gt.1.or.matarg.gt.1)then
                  if(i1.gt.0.and.i2.gt.0.and.i3.gt.0)eef=eef-prom !baryons
                  if(i1.lt.0.and.i2.lt.0.and.i3.lt.0)eef=eef+prom !antibar
                endif
                val=0
                if(pp.ne.0.)val=eef*pt/pp
              endif
              if( (eta.ge.xpara(ishift+1,n).and.eta.le.xpara(ishift+2,n)
     *             .or. eta.ge.xpara(ishift+9,n)
     *                .and.eta.le.xpara(ishift+10,n) )
     *           .and.pt .gt.xpara(ishift+3,n)
     *             .and.pt .lt.xpara(ishift+4,n)   )then
                sum=sum+val
c      if(ishift.ne.0)print *,io,i,idptl(i),istptl(i),eta,val,sum,istuse(n)
              endif
            endif
          endif
        endif

      enddo

      z=sum*xpara(ishift+5,n)/xpara(ishift+6,n) 
      
      if(io2.eq.2)then  ! percentile 0-5-10-...-100
        kkx=0
        do kk=1,10
          if( xpara(ishift+11+kk,n).lt.z
     .              .and.z.le. xpara(ishift+12+kk,n) )kkx=kk
        enddo
        z=  100.-(kkx-0.5)*10.
      endif

      ypara(ishift+1,n)=z      

c      if(xpara(ishift+2,n).eq.2.5.and.xpara(ishift+3,n).eq.0.5
c     ..and.z.gt.180)call utstop('test&')
      
      !if(io1.eq.2)print*,'XAN EVA ',z,icas
      !if(io1.eq.3)print*,'XAN NCHARM ',z,icas
      
      end

      function fmux(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11)
      !   xpara(1,n) ... etamin
      !   xpara(2,n) ... etamax
      !   xpara(3,n) ... ptmin
      !   xpara(4,n) ... ptmax
      !   xpara(5,n) ... factor
      !   xpara(6,n) ... divisor
      !   xpara(7,n) ... absolute value of eta (1)
      !   xpara(8,n) ... yboost if >< 0.0001 (boost back)
      !   xpara(9,n) ... etamin2 if >< 0.
      !   xpara(10,n) ... etamax2 if >< 0.
      !   xpara(11,n) ... x0 = normal
      !                   x2 = percentile  
      !                   0x = multiplicity
      !                   2x = Et
      include "epos.incxan"
      n=1
      xpara(1,n)=p1
      xpara(2,n)=p2
      xpara(3,n)=p3
      xpara(4,n)=p4
      xpara(5,n)=p5
      xpara(6,n)=p6
      xpara(7,n)=p7
      xpara(8,n)=p8
      xpara(9,n)=p9
      xpara(10,n)=p10
      xpara(11,n)=p11
      call mux(n,0)
      fmux=ypara(1,n)
      end

c----------------------------------------------------------------------
      function helix(xi)
c----------------------------------------------------------------------
      include "epos.inc"
      parameter (maxpv=100000)
      real phixx(maxpv),etaxx(maxpv),ptxx(maxpv)
      helix=0
      ii=0
      
      do iloo=maproj+matarg+1,nptl
     
        i=iloo   
      
       istxxx=0
       if(istfor.gt.-999)istxxx=istfor
       if(istptl(i).eq.istxxx)then
          px=pptl(1,i)
          py=pptl(2,i)
          pt=pptl(1,i)**2+pptl(2,i)**2
          if(pt.gt.0.)then
            pt=sqrt(pt)
            p3=pptl(3,i)
            if(p3.ne.0.)then
              theta=atan(pt/p3)
            else
              theta=pi/2.
            endif
            if(theta.lt.0.)theta=theta+pi
            sth=sin(theta*0.5)
            cth=cos(theta*0.5)
            if(sth.ne.0..and.cth.ne.0.)then
              tth=sth/cth
              if(tth.gt.0.)then
                eta=-log(tth)
              else
                eta=1000.
              endif
            else
              eta=1000.
            endif
          else
            eta=1000.
          endif
          if(abs(eta).lt.2.5)then
            call idchrg( 13 ,ideposf( 5 ,i),ch)
            if(abs(ch).gt.0.1.and.pt.gt.0.100)then
              ii=ii+1
              if(ii.gt.maxpv)then
                write(ifmt,*)
     .          '***** ERROR: PairVariables: arrays too small'
                 stop'##### PairVariables: arrays too small #####'
              endif
              phixx(ii)=polar(px,py)
              etaxx(ii)=eta
              ptxx(ii)=pt
              !print*,'+++++',ii,pt,idxx(ii)
            endif    
          endif  
        endif

      enddo      

      if(ii.gt.0)then
      do m=1,ii
      do n=1,ii 
      if(n.ne.m)then
        d=xi*(etaxx(m)-etaxx(n))-(phixx(m)-phixx(n))
        helix=helix+cos(d)
      endif
      enddo
      enddo
      helix=helix/ii
      endif
      !print*,'helix',ii,xi,helix
      end
      
c----------------------------------------------------------------------
      subroutine PairVariables !'pai' 'co2' 'co3'
c----------------------------------------------------------------------
      include "epos.inc"
      parameter(kkkmax=6)
      common/cpairs/paievt(kkkmax),co2evt(kkkmax),co3evt(kkkmax)
      parameter (mxpaih=10)
      parameter (maxpt=50)
      common/cpairs2/
     .    paih(kkkmax,mxpaih,maxpt)
     .   ,co2h(kkkmax,mxpaih,maxpt),co3h(kkkmax,mxpaih,maxpt)
      parameter (maxpv=100000)
      integer idxx(maxpv),irxx(maxpv)
      real phixx(maxpv),etaxx(maxpv),ptxx(maxpv)
      real etamax(kkkmax/2),pttmin(kkkmax/2),pttmax(kkkmax/2)
      data etamax / 1.0 ,  2.4 ,  2.4 / 
      data pttmin / 0.  ,  0.  ,  0.3 / 
      data pttmax / 1e30,  1e30,  3.0 /
      yboo=rapcms
      delpt=10./maxpt
      do k=1,kkkmax
      paievt(k)=0
      co2evt(k)=0
      co3evt(k)=0
       do j=1,maxpt
        do i=1,mxpaih
        paih(k,i,j)=0.
        co2h(k,i,j)=0.
        co3h(k,i,j)=0.
        enddo 
       enddo
      enddo

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! kk1 :
      !        1           2
      !      ist=0       ist=8
      ! kk2 :
      !        1            2            3
      !      |eta|<1     |eta|<2.4   |eta|<2.4
      !    |DelEta|>1   |DelEta|>2  |DelEta|>2
      !    pt_ref all   pt_ref all  pt_ref 0.3-3     
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do kk1=1,2
      do kk2=1,kkkmax/2
      kkk=(kk1-1)*(kkkmax/2)+kk2
      
      ii=0
      do iloo=maproj+matarg+1,nptl
     
        i=iloo   

       istxxx=0
       if(istfor.gt.-999)istxxx=istfor
       iok=0
       if( kk1.eq.1 .and. (istptl(i).eq.istxxx
     .   .or.idptl(i).eq.331.and.istptl(i).eq.9.and.ityptl(i).eq.80)
     .                                                          ) iok=1
     
       if( kk1.eq.2 .and. istptl(i).eq.8      ) iok=1
       
       if(iok.eq.1)then
          px=pptl(1,i)
          py=pptl(2,i)
          pt2=pptl(1,i)**2+pptl(2,i)**2
          p4=pptl(4,i)
          p3=pptl(3,i)
          if(abs(yboo).gt.0.0001)then
          p4= cosh(yboo) * pptl(4,i) + sinh(yboo) *  pptl(3,i)
          p3= sinh(yboo) * pptl(4,i) + cosh(yboo) *  pptl(3,i)
          endif
          pp=sqrt(pt2+p3**2)
          if(pt2.gt.0.)then
            pt=sqrt(pt2)
            eta=sign(1.,p3)*alog((pp+abs(p3))/pt)
          else
            pt=0.
            eta=1000.
          endif
           if(abs(eta).lt.etamax(kk2))then
              if(idptl(i).eq.31.or.idptl(i).eq.25)then
                print*,'ERROR 11072016b',iorptl(i),jorptl(i),iloo,i
     .            ,idptl(i),istptl(i),ityptl(i)
              endif
            id=ideposf( 6 ,i)
            call idchrg( 14 ,id,ch)
            if(abs(ch).gt.0.1.or.id.eq.331.or.id.eq.20
     .                .or.abs(id).eq.2130)then
              ii=ii+1
              if(ii.gt.maxpv)then
                write(ifmt,*)
     .          '***** ERROR: PairVariables: arrays too small'
                 stop'##### PairVariables: arrays too small #####'
              endif
              phixx(ii)=polar(px,py)
              etaxx(ii)=eta
              ptxx(ii)=pt
              irxx(ii)=0
              if(pt.ge.pttmin(kk2).and.pt.le.pttmax(kk2)
     .        .and.abs(ch).gt.0.1)irxx(ii)=1
              idxx(ii)=id
              if(idxx(ii).eq.331)then
                if(.not.(istptl(i).eq.9.and.ityptl(i).eq.80))idxx(ii)=0
              endif
            endif    
          endif  
        endif

      enddo      

      do m=1,ii
      if(irxx(m).eq.1)then !reference
      do n=1,ii 
      
        iok=0
        if( kk2.eq.1 .and. 
     .         (  (etaxx(m).lt.-0.5 .and. etaxx(n).gt.0.5)
     .       .or. (etaxx(n).lt.-0.5 .and. etaxx(m).gt.0.5)   )   )iok=1
        if( (kk2.eq.2.or.kk2.eq.3) .and.
     .                    abs(etaxx(m)-etaxx(n)).gt.2            )iok=1
        if(iok.eq.1)then
          if(abs(etaxx(m)-etaxx(n)).lt.0.999)
     .    stop'##### ERROR 04112011 #####' 
          !m.ne.n automatic in this case
          co2mn=cos(2*(phixx(m)-phixx(n)))
          co3mn=cos(3*(phixx(m)-phixx(n)))
          if(irxx(n).eq.1)then
            paievt(kkk)=paievt(kkk)+1
            co2evt(kkk)=co2evt(kkk)+co2mn
            co3evt(kkk)=co3evt(kkk)+co3mn
          endif
          ian=abs(idxx(n))
          idn=    idxx(n)
          in=ptxx(n)/delpt+1
          if(in.ge.1.and.in.le.maxpt)then
            iid=0
            if(ian.eq.120)iid=1
            if(ian.eq.1120)iid=2
            if(ian.eq.130)iid=3
            if(ian.eq.1130)iid=4
            if(idn.eq.20)iid=5
            if(ian.eq.2130)iid=6
            if(ian.eq.2330)iid=7
            if(ian.eq.3331)iid=8
            if(ian.eq.331)iid=9
            if(iid.ne.0)then
              paih(kkk,iid,in)=paih(kkk,iid,in)+1
              co2h(kkk,iid,in)=co2h(kkk,iid,in)+co2mn
              co3h(kkk,iid,in)=co3h(kkk,iid,in)+co3mn
              !print*,'+++hadron++',in, idxx(n),ptxx(n)
            endif
            iid=0
            if(ian.eq.120.or.ian.eq.130.or.ian.eq.1120
     .          .or.ian.eq.1130.or.ian.eq.2230
     .          .or.ian.eq.2330.or.ian.eq.3331)iid=10
            if(iid.ne.0)then
              paih(kkk,iid,in)=paih(kkk,iid,in)+1
              co2h(kkk,iid,in)=co2h(kkk,iid,in)+co2mn
              co3h(kkk,iid,in)=co3h(kkk,iid,in)+co3mn
              !print*,'+++hadron++',in, idxx(n),ptxx(n)
            endif
          endif
        endif
      
      enddo
      endif
      enddo
      
      enddo  !~~~~~kk2~~~~~~~~~~
      enddo  !~~~~~kk1~~~~~~~~~~

      !test removed, see earlier versions
      end
      
c----------------------------------------------------------------------
      subroutine detphi(n,ishift)    ! 'detphi'
c----------------------------------------------------------------------
c  dN_pair / dDeltaPhi per trigger
c----------------------------------------------------------------------
c   xpara(1,n) ... should be NOT 0 !!!
c   xpara(2,n) ... id of particle
c   xpara(3,n) ... ist of particle
c   xpara(4,n) ... min eta
c   xpara(5,n) ... max eta
c   xpara(6,n) ... min pt 
c   xpara(7,n) ... max pt
c   xpara(8,n) ... min deta
c   xpara(9,n) ... max deta
c   xpara(10,n) ... boost
c   xpara(11,n) ... random angle (if 1)
c----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incxan"
      
      parameter (maxpv=100000)
      real phixx(maxpv),etaxx(maxpv),ptxx(maxpv)
      integer idxx(maxpv)

      if(nint(xpara(ishift+1,n)).eq.0)stop'##### ERROR 14102014#####'
      ! 0 not supported any more (referred to former 'phi')

      maxbi=20 
      if(maxbi.gt.myyarr)stop'04042021'
      del=pi/maxbi  
      nxxx=nhisxxx(n)
      
      id= nint(xpara(ishift+2,n))
      call idchrg( 15 ,id,ch)
      if(abs(ch).lt.0.1.and.id.ne.9970)stop'##### ERROR 26102014 #####'
      ist=nint(xpara(ishift+3,n))
      eta1= xpara(ishift+4,n)
      eta2= xpara(ishift+5,n)
      pt1=  xpara(ishift+6,n)
      pt2=  xpara(ishift+7,n)
      deta1=xpara(ishift+8,n)
      deta2=xpara(ishift+9,n)
      yboo= xpara(ishift+10,n)
      ira= nint( xpara(ishift+11,n) )
      
      if(abs(eta1+eta2).gt.1e-4)stop'##### ERROR 12102014a #####'
      !!!!!  mixed event correction done for symmetric case  !!!!

      if(deta2.ge.eta2-eta1)stop'##### ERROR 12102014b #####'
      
      ii=0
      do iloo=maproj+matarg+1,nptl
     
       i=iloo   

       istxxx=0
       if(istu1(istuse(n)).gt.-999)istxxx=istu1(istuse(n))
       if(istptl(i).eq.istxxx)then

        px=pptl(1,i)
        py=pptl(2,i)
        pt2=pptl(1,i)**2+pptl(2,i)**2
        p4=pptl(4,i)
        p3=pptl(3,i)
        if(abs(yboo).gt.0.0001)then
        p4= cosh(yboo) * pptl(4,i) + sinh(yboo) *  pptl(3,i)
        p3= sinh(yboo) * pptl(4,i) + cosh(yboo) *  pptl(3,i)
        endif
        pp=sqrt(pt2+p3**2)
        if(pt2.gt.0.)then
          pt=sqrt(pt2)
          eta=sign(1.,p3)*alog((pp+abs(p3))/pt)
        else
          pt=0.
          eta=1000.
        endif
        if(eta.gt.eta1.and.eta.lt.eta2)then
        if(pt.gt.pt1.and.pt.lt.pt2)then
          call idchrg( 16 ,ideposf( 7 ,i),ch)
          if(abs(ch).gt.0.1)then
            ii=ii+1
            if(ii.gt.maxpv)then
              write(ifmt,*)
     .          '##### ERROR 12092014 : detphi: arrays too small'
               stop'##### ERROR 12092014 #####'
            endif
            phixx(ii)=polar(px,py)
            etaxx(ii)=eta
            ptxx(ii)=pt
            idxx(ii)=idptl(i)
          endif    
        endif  
        endif  

       endif
       
      enddo      

      !------------------------------------------------------------
      ! Approximate treatment mixed event correction via 1/C(deta) factor
      !
      ! The precise teatment is obtained by calling this routine with 
      !   xpara 11 0 for REAL and the  xpara 11 1 for MIXED events
      !    and then REAL / MIXED * MIXED(deta=0)   
      !
      !  CAlculation of C:
      !  C=B(deta)/B(0) ; A=dN/deta ~ const
      !  B(deta) = 1/(eta2-eta1) \int \int A delta(eta-etap-deta) d_etap d_eta
      !  =>  C(deta) = 1 - abs(deta)/(eta2-eta1)
      !------------------------------------------------------------
      
      do m=1,ii
      do k=1,ii 
        
        deta=etaxx(m)-etaxx(k)
        if(deta.gt.deta1.and.deta.lt.deta2 )then
          cmix= 1 - abs(deta)/(eta2-eta1)
          dphi= phixx(m)-phixx(k)
          if(dphi.gt.pi)dphi=dphi-2*pi
          if(dphi.lt.-pi)dphi=dphi+2*pi
          if(ira.eq.1)dphi=(rangen()*2-1)*pi
          if(dphi.ge.0)then
            in=dphi/del+1
            if(in.ge.1.and.in.le.maxbi)then
              ida=abs(idxx(k))
              if(id.eq.ida
     .          .or.id.eq.9970)then
                yyarr(nxxx,in)=yyarr(nxxx,in)+1/cmix
              endif
            endif
          endif
        endif
      
      enddo
      enddo
      
      !-------------------------------------------
      !  / ii    per trigger
      !   double counting taken care of by / ii
      !-------------------------------------------
      do j=1, maxbi  
        if(ii.gt.0) then
          yyarr(nxxx,j)=yyarr(nxxx,j)/real(ii) 
        endif
      enddo

      end

c----------------------------------------------------------------------
      subroutine CorH(n) !'corh' 'corhevt'
c----------------------------------------------------------------------
c correlation function of (mainly) heavy flavor pairs (partons or hadrons)
c----------------------------------------------------------------------
c   xpara(2,n) ... id of the first particle (trigger)
c   xpara(3,n) ... id of the second particle
c   xpara(4,n) ... min rap/pseudorap part 1
c   xpara(5,n) ... max rap/pseudorap part 1
c   xpara(6,n) ... min pt part 1
c   xpara(7,n) ... max pt part 1
c the following ones are not mandatory
c   xpara(1,n) ... 0 (default): normalized by the number of pairs, use rap 
c                  1 : normalized per trigger & use rap for trigger and  eta for associate
c                  2 : take normalization value from preceding call 
c   xpara(8,n) ... min rap/pseudorap part 2
c   xpara(9,n) ... max rap/pseudorap part 2
c   xpara(10,n) ... min pt part 2
c   xpara(11,n) ... max pt part 2
c   xpara(12,n) ... max |delta eta| if > 0.0001, otherwise no condition 
c   xpara(15,n) ... ist part 1 and 2
c   xpara(16,n) ... itymin part 1 and 2
c   xpara(17,n) ... itymax part 1 and 2
c   xpara(18,n) ... 1+jor min part 1 and 2
c   xpara(19,n) ... 1+jor max part 1 and 2
c   xpara(20,n) ... iptcond =1 if > 0.5, otherwise no condition 
c   xpara(21,n) ... use rapidity rather than phi if > 0.5 
c   ypara(1,n) .... normalization
c----------------------------------------------------------------------

      include "epos.inc"
      include "epos.incxan"

      parameter (maxpv=100000)
      real phixx(maxpv),etaxx(maxpv),ptxx(maxpv)
      integer ixx(maxpv),irexx(maxpv)
      logical IOK1,IOK2,IOK3,ptcond
      common/cypara1save/ypara1save 

      inorh=nint(xpara(1,n))
      if(inorh.eq.2)then
        ypara(1,n)=ypara1save
        return
      endif

      yboo=rapcms

      nxxx=nhisxxx(n)
      maxcorr=40 
      if(maxcorr.gt.myyarr)stop'10062020'
      irap=0
      binmin=0.
      binmax=pi
      if(xpara(21,n).gt.0.5)then
       irap=1
       binmin=0.
       binmax=2.
      endif

      idd1=nint(xpara(2,n))
      idd2=nint(xpara(3,n))

      rap1min=xpara(4,n)
      rap1max=xpara(5,n)
      !to keep the old optns format:
      if(xpara(9,n)-xpara(8,n).gt.0.0001)then  
        rap2min=xpara(8,n)
        rap2max=xpara(9,n)
      else
        rap2min=xpara(4,n)
        rap2max=xpara(5,n)
      endif
      ptr1min=xpara(6,n)
      ptr1max=xpara(7,n)
      !to keep the old optns format:
      if(xpara(11,n)-xpara(10,n).gt.0.0001)then  
        ptr2min=xpara(10,n)
        ptr2max=xpara(11,n)
      else
        ptr2min=xpara(6,n)
        ptr2max=xpara(7,n)
      endif

      deta=1e30
      if(xpara(12,n).gt.0.0001)deta=xpara(12,n)

      iptcond=0
      if(xpara(20,n).gt.0.5)iptcond=1

      ist12=0
      if(xpara(15,n).gt.0.0001)ist12=nint(xpara(15,n))
      ity12min=0
      if(xpara(16,n).gt.0.0001)ity12min=nint(xpara(16,n))
      ity12max=0
      if(xpara(17,n).gt.0.0001)ity12max=nint(xpara(17,n))
      jor12min=0 
      if(xpara(18,n).gt.0.0001)jor12min=nint(xpara(18,n))
      jor12max=0 
      if(xpara(19,n).gt.0.0001)jor12max=nint(xpara(19,n))

      ii=0

      do iloo=maproj+matarg+1,nptl
     
        i=iloo   

        j=i 
        call getItyJor25(j,ityxx,jorxx) !ity,jor of 25 parton

        IOK1 = ist12.eq.0 .or. (ist12.gt.0.and.istptl(i).eq.ist12)
        IOK2 = ity12max.eq.0 .or. (ity12max.gt.0
     .      .and.ityxx.ge.ity12min.and.ityxx.le.ity12max) 
        IOK3 = jor12max.eq.0 .or. (jor12max.gt.0
     .      .and.1+jorxx.ge.jor12min.and.1+jorxx.le.jor12max) 
        id=ideposf( 8 ,i)
        idx=9901
        !special code:
        ! D0 or D+ or  D+* (including antiparticles)
        ! excluding   D0 and D+ from D* decay
        idu=9902
        !special code:
        ! electrons from D,B
        idv=9903
        !special code:
        ! electrons from D
        idw=9904
        !special code:
        ! electrons from B
        idy=9970 
        !special code: charged particles
        if(idd1.eq.idx.and.istptl(i)/2.eq.0)then!D-mesons
          if(iorptl(i).gt.0) then
            ifath=iabs(idptl(iorptl(i)))
          else
            ifath=0
          endif
          if(iabs(id).eq.241)id=idx
          if(iabs(id).eq.140.and.ifath.ne.241)id=idx
          if(iabs(id).eq.240.and.ifath.ne.241)id=idx
        endif
        if((idd1.eq.idu.or.idd1.eq.idv.or.idd1.eq.idw)
     .   .and.istptl(i).eq.0)then!electrons from D,B
          if(iorptl(i).gt.0) then
            ifath=iabs(idptl(iorptl(i)))
          else
            ifath=0
          endif
          if(iabs(id).eq.12.and.idd1.eq.idu)then
              if(ifath.eq.240.or.ifath.eq.140.or.ifath.eq.340
     .        .or.ifath.eq.2140)id=idu
              if(ifath.eq.250.or.ifath.eq.150.or.ifath.eq.350
     .        .or.ifath.eq.2150)id=idu
          elseif(iabs(id).eq.12.and.idd1.eq.idv)then
              if(ifath.eq.240.or.ifath.eq.140.or.ifath.eq.340
     .        .or.ifath.eq.2140)id=idv
          elseif(iabs(id).eq.12.and.idd1.eq.idw)then
              if(ifath.eq.250.or.ifath.eq.150.or.ifath.eq.350
     .        .or.ifath.eq.2150)id=idw
          endif
        endif
        if(idd1.eq.idy.and.istptl(i)/2.eq.0)then
          call idchrg( 17 ,id,ch)
          if(abs(ch).gt.0.1)then
              id=idy
          endif
        endif
        if(idd2.eq.idy.and.istptl(i)/2.eq.0)then
          call idchrg( 17 ,id,ch)
          if(abs(ch).gt.0.1)then
              id=idy
          endif
        endif
 
        if(IOK1.and.IOK2.and.IOK3
     .  .and.(id.eq.idd1.or.id.eq.idd2))then

          px=pptl(1,i)
          py=pptl(2,i)
          p3=pptl(3,i)
          p4=pptl(4,i)
          if(abs(yboo).gt.0.0001)then
            p4= cosh(yboo) * pptl(4,i) + sinh(yboo) *  pptl(3,i)
            p3= sinh(yboo) * pptl(4,i) + cosh(yboo) *  pptl(3,i)
          endif
          pt=pptl(1,i)**2+pptl(2,i)**2
         !if(sqrt(pt).gt.3.)
         !.print*,'CorH',iloo,id,istptl(i),ityptl(i),px,py,p3,p4
          amt=pptl(5,i)**2+pt
          if(amt.gt.0..and.pptl(4,i)+abs(pptl(3,i)).gt.0.)then 
            amt=sqrt(amt)
            rap=sign(1.,p3)*log((p4+abs(p3))/amt)
          else
            rap=10000
          endif
          rapxx=rap
          if(pt.gt.0.)then
            pt=sqrt(pt)
            eta=sign(1.,p3)*log((sqrt(p3**2+pt**2)+abs(p3))/pt)
          else
            eta=10000.
          endif

          ireg=0
          if(id.eq.idd1)then
            if(  rap.ge.rap1min.and.rap.le.rap1max
     .      .and. pt.ge.ptr1min.and. pt.le.ptr1max)ireg=1
          endif
          if(inorh.eq.1)rap=eta  
          if(id.eq.idd2)then
            if(  rap.ge.rap2min.and.rap.le.rap2max
     .      .and. pt.ge.ptr2min.and. pt.le.ptr2max)ireg=ireg+2
          endif

          if(ireg.gt.0) then
            ii=ii+1
            if(ii.gt.maxpv)then
              write(ifmt,*)
     .        '***** ERROR: PairVariables: arrays too small'
                stop'##### PairVariables: arrays too small #####'
            endif
            ixx(ii)=iloo
            irexx(ii)=ireg 
            if(irap.eq.0)then
              phixx(ii)=polar(px,py) !bg polar -> phi=[0,2pi]
            else
              phixx(ii)=rapxx
            endif
            etaxx(ii)=eta
            ptxx(ii)=pt
            if((ireg.eq.1.or.ireg.eq.3).and.inorh.eq.1)
     .        ypara(1,n)=ypara(1,n)+1 ! per trigger
c              !~~~~~~~~~~~~~~~~~~~~~~~~~~~
c              if(abs(id).eq.4)
c     .        print*,'TESTc',i,id,istptl(i),ityxx,rap,pt,phixx(ii)
c              !~~~~~~~~~~~~~~~~~~~~~~~~~~~
c          else
c              !~~~~~~~~~~~~~~~~~~~~~~~~~~~
c              if(abs(id).eq.4)
c     .        print*,'-----',i,id,istptl(i),ityxx,rap,pt
c              !~~~~~~~~~~~~~~~~~~~~~~~~~~~
          endif

        endif

      enddo

      do mii=1,ii   
      do nii=1,ii 
        irex1=irexx(mii)
        irex2=irexx(nii) 
        ireg=1
        if(mii.eq.nii)ireg=0
        i=ixx(mii)
        j=ixx(nii)
        call getiorptl(i,iori)
        call getiorptl(j,iorj)
        if(iorj.eq.i)ireg=0 !exclude: father of 2 = 1
        if(iorj.gt.0)then
          call getiorptl(iorj,ioriorj)
          if(ioriorj.eq.i)ireg=0 !exclude: grandfather of 2 = 1
        endif
        if(iori.eq.j)ireg=0 !exclude: father of 1 = 2
        if(ireg.eq.1
     .    .and.(irex1.eq.1.or.irex1.eq.3)
     .    .and.(irex2.eq.2.or.irex2.eq.3))then
          deleta=abs(etaxx(mii)-etaxx(nii))
          ptcond=.true.
          if(iptcond.eq.1)ptcond=ptxx(mii).gt.ptxx(nii)
          if(deleta.le.deta.and.ptcond)then
            delphi=abs(phixx(mii)-phixx(nii))
            if(irap.eq.0.and.delphi.gt.pi) delphi=2.*pi-delphi
            nb=int(delphi/(binmax-binmin)*real(maxcorr)) +1
            if(nb.ge.maxcorr+1)nb=maxcorr
            !print*,'CorH',mii,nii,delphi/pi,nb
            yyarr(nxxx,nb)=yyarr(nxxx,nb)+1
            if(inorh.eq.0)ypara(1,n)=ypara(1,n)+1
            !call getidptl(i,id1)
            !call getidptl(j,id2)
            !print*,'TESTcorh',inorh,ypara(1,n),yyarr(nxxx,nb),nb,id1,id2
          endif
        endif
      enddo
      enddo
      ypara1save=ypara(1,n)
      end

c----------------------------------------------------------------------
      subroutine mupurp_soupai(n)    ! 'mupurp' source function of pairs
c----------------------------------------------------------------------
c   xpara(2,n) ... max relative momentum in pair CMS
c----------------------------------------------------------------------

      include "epos.inc"
      include "epos.incxan"
      logical go,LongLivPtl,NoLongLivParent
      
      parameter (maxsoupai=400)
      real pxx(4,maxsoupai),xxx(4,maxsoupai)
      real p(4),q(4),r(4),x(4),y(4)

      qrmax=xpara(2,n)

      del=0.25   
      nxxx=nhisxxx(n)
            
      ii=0
      do iloo=maproj+matarg+1,nptl
 
       i=iloo   

       go=.false.
       call setIstXY(0,istuse(n),istxxx,istyyy)
       if(istptl(i).eq.istxxx)go=.true.
       if(noweak(n).eq.1)then
         if(LongLivPtl(istyyy,i))go=.true.
       elseif(noweak(n).eq.3)then
         if(LongLivPtl(istxxx,i))go=.false.
       endif
       if(go)go=NoLongLivParent(noweak(n),i)

       if(go)then
         id=ideposf( 26 ,i)
         if(id.eq.1120)then
            ii=ii+1
            if(ii.gt.maxsoupai)then
              write(ifmt,*)
     .        '##### ERROR 15072017 : soupai : arrays too small'
              stop'##### ERROR 15072017 #####'
            endif
            do idim=1,4  
              pxx(idim,ii)=pptl(idim,i)
              xxx(idim,ii)=xorptl(idim,i) 
            enddo
            !print*,ii,xxx(1,ii),xxx(2,ii),xxx(3,ii)
         endif  
       endif
       
      enddo      

      npairs=0
      do m=1,ii-1
      do k=m+1,ii
        do idim=1,4 
          p(idim)= pxx(idim,m) + pxx(idim,k)
          q(idim)= pxx(idim,m)
          r(idim)=               pxx(idim,k)
          x(idim)= xxx(idim,m)
          y(idim)=               xxx(idim,k)
        enddo
        p55=p(4)**2-p(1)**2-p(2)**2-p(3)**2
        if(p55.gt.0)then 
          p5=sqrt(p55) 
          call utlob3(1,p(1),p(2),p(3),p(4),p5,q(1),q(2),q(3),q(4))
          call utlob3(1,p(1),p(2),p(3),p(4),p5,r(1),r(2),r(3),r(4))
          call utlob3(1,p(1),p(2),p(3),p(4),p5,x(1),x(2),x(3),x(4))
          call utlob3(1,p(1),p(2),p(3),p(4),p5,y(1),y(2),y(3),y(4))
          qr2= (q(1)-r(1))**2 + (q(2)-r(2))**2 + (q(3)-r(3))**2
          !print*,sqrt(qr2)
          if(qr2.le.qrmax**2)then
            dd=(x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2 
            d=sqrt(dd)
            !print*,m,k,'   ',x(1),y(1),'   ',x(2),y(2)
            !.                ,'   ',x(3),y(3),'   ',d 
            aa=d/del+1
            if(aa.le.float(myyarr))then
              in=aa
              if(in.le.0)print*,'WARNING soupai ',qrmax,dd,d,del,in
              if(in.le.myyarr)then 
                yyarr(nxxx,in)=yyarr(nxxx,in)+1
              endif
            endif
          endif
        endif
      enddo
      enddo
      
      end

c----------------------------------------------------------------------
      subroutine mupurp(n,ishift)  !'mupurp'
c----------------------------------------------------------------------

      include "epos.inc"
      include "epos.incxan"
      
      if( nint( xpara(1,n) ) .eq. 1001 ) stop'ERROR 11032021' !removed in 3415
      if( nint( xpara(1,n) ) .eq. 1002 ) call mupurp_soupai(n)
      if(nint(xpara(ishift+1,n)).eq. 1003 )call mupurp_cum(n,ishift)
      if(nint(xpara(ishift+1,n)).eq. 1004 )call mupurp_scp(n,ishift)!gs
      if(nint(xpara(ishift+1,n)).eq. 1005 )call mupurp_sc(n,ishift)!gs
      if(nint(xpara(ishift+1,n)).eq. 1006 )call mupurp_susc(n,ishift)

      end

c-----------------------------------------------------------------------------
      subroutine mupurp_susc(n,ishift)  !'mupurp' via net-cumulants N^i  [JJ] 
c-----------------------------------------------------------------------------
      !--------------------------------------------------------------------
      ! Calculation of net-particles cumulants of order i=1-4 
      !  vs. event variables (Nch, Npart...)
      !   - Refs : arxiv:1903.05370v3 / 1312.3572v2
      !   - Inputs :
      !   n = histogram number
      !   xpara(1,n)    1006 
      !   xpara(2,n)    ptmin 
      !   xpara(3,n)    ptmax 
      !   xpara(4,n)    rapidity / pseudorapidity cut 
      !                    0  : rapidity
      !                    1  : pseudorapidity
      !   xpara(5,n)    |y/eta|_max 
      !   xpara(6,n)    feed-down correction
      !                    0  : no feed-down correction (='weak')
      !                    1  : geometrical cut (DCA < 1cm) (see Ref.1)
      !                    2  : full feed-down correction (='noweak')
      !   xpara(7,n)    status of particles considered (ist)
      !                    0  : FS particles+parents (DEFAULT)
      !                    1  : before hadronic cascades
      !   xpara(8,n)    origin of particles considered (ity)
      !                    0  : all particles          (ity=[20-61])
      !                    1  : particles from core    (ity=[60-61]) 
      !                    2  : particles from corona  (ity=[20-59]) 
      !   xpara(9,n)    charges (c) of net-cumulant N^i_(c) considered:
      !                    (c) = a | ab | abc | abcd
      !                     i  = 1 |  2 |   3 |    4
      !
      !                 with a,b,c,d  =  1|2|3|4   |5   |6   |7   |8 |9
      !                 corrspndg to --> B|Q|S|pion|kaon|Prtn|Lbda|Xi|Omg
      !--------------------------------------------------------------------
      include "epos.inc"
      include "epos.incxan"
        logical FDgo,NoLongLivParent
        real pT, pT2, mT, mT2, chrg
        character FDC*4, cut*3, chrg_cum*4
        character(len=3), parameter :: y='rap', eta='eta'
        integer order, dq, dS, sumq
        !--- Highest order calculable by the routine ---
        !-----------------------------------------------
        integer,parameter :: max_order = 4
        !--- Number of flavours considered ---
        !-------------------------------------
        integer,parameter :: dim_fl = 9
        !--- Matrix of net-cumlnt orders / corresponding values ---
        ! cumlts(1,x) : ptl IDs (B=1, Q=2, S=3)
        ! cumlts(2,x) : value
        !    * N_deta (number of (pseudo)rap intervals considered) 
        !----------------------------------------------------------
        integer, dimension(2,max_order,nint(xpara(10,n))) :: cumlts
        !--- List of corresponding IDs ---
        ! chrgID_trans(1,x) = (1, 2, 3, 1120, 120, 130, 2130, 2330, 3331)
        !---------------------------------
        integer :: chrgID_trans(dim_fl)
        !--- List of (pseudo)rapidity intervals to consider ---
        !------------------------------------------------------²
        real :: dcut_list(nint(xpara(10,n)))

        !--- Initialisations ---
        !-----------------------
        chrgID_trans = (/ 1, 2, 3, 120, 130, 1120, 2130, 2330, 3331 /)
        cumlts = 0 ! initialise the whole matrix to 0
        
        !--- Assigning input variables ---
        !---------------------------------
        !-- pT cuts --
        ptmin=xpara(ishift+2,n)
        ptmax=xpara(ishift+3,n)
        !-- y/eta cuts --
        if(int(xpara(ishift+4,n)).eq.0)then
          cut='rap'
          cutmax=xpara(ishift+5,n)
        elseif(int(xpara(ishift+4,n)).eq.1)then
          cut='eta'
          cutmax=xpara(ishift+5,n)
        else
          STOP"ERROR 20042022a (mupurp_susc: xpara(4)=0/1)"
        endif
        !-- Feed-down correction --
        FDC='NO'   !weak
        if(xpara(ishift+6,n).eq.1)then
          FDC='DCA' !DCA < 1cm
        elseif(xpara(ishift+6,n).eq.2)then
          FDC='FULL' !noweak
        elseif(xpara(ishift+6,n).ne.0)then
          STOP"ERROR 20042022b (mupurp_susc: xpara(6)=0/1/2)"
        endif
        !-- Status of particles considered --
        choice=nint(xpara(ishift+7,n))
        !FS particles
        if(choice.eq.0)then 
          !Consider ist=0 + ist=1 too bc of weak dkying ptls
          ist1=0
          ist2=1
        !Before had. cas.
        elseif(choice.eq.1)then
          !Consider ist=3 only -> ptls just after hadronstn/fragmtn
          !(not ist=6/8 bc "no hacas" alternative scenario)
          ist1=3
          ist2=3
        else
          STOP"ERROR 20042022c (mupurp_susc: xpara(7)=0/8)"
        endif
        !-- Origin of particles considered --
        itymin=20
        itymax=61
        if(xpara(ishift+8,n).eq.1)then !only core
          itymin=60
        elseif(xpara(ishift+8,n).eq.2)then !only corona
          itymax=59
        elseif(xpara(ishift+8,n).ne.0)then
          STOP"ERROR 20042022d (mupurp_susc: xpara(8)=0/1/2)"
        endif
        !-- Particle(s)/charge(s) + order of cumulant --
        if(nint(xpara(ishift+9,n)).lt.0
     .     .and.nint(xpara(ishift+9,n)).gt.9999)then
          STOP"ERROR 20042022e (mupurp_susc: xpara(9)=a/ab/abc/abcd>0)"
        endif
          !-Convert xpara(9) as a string
        write(chrg_cum,'(i4)') nint(xpara(ishift+9,n)) 
        chrg_cum=adjustl(chrg_cum) !adjust the string on left
        order=len_trim(chrg_cum)
          !-Fill cumlts(1) with IDs corresponding to numbers in xpara(9)
        do i=1,order
          read(chrg_cum(i:i),'(i1)') cumlts(1,i,1)  !convert chrg_cum(i) into int
          cumlts(1,i,1)=chrgID_trans(cumlts(1,i,1)) !convert number into ID
        enddo
        !-- Calculation of (pseudo)rapidity intervals --
        Ndcut=nint(xpara(10,n))
        do i=1,Ndcut
          dcut_list(i)=(cutmax/Ndcut)*i 
        enddo

        !-------------------------
        !   Loop over particles    
        !...
        do iloo=maproj+matarg+1,nptl 
     
          i=iloo   
          
          !-Testing particle status
          if(istptl(i).eq.ist1.or.istptl(i).eq.ist2)then
            !-Testing particle's origin
            if(ityptl(i).ge.itymin.and.ityptl(i).le.itymax)then
              !-pT cut
              pT2=pptl(1,i)**2+pptl(2,i)**2
              pT=0.
              if(pT2.gt.0.)pT=sqrt(pT2)
              if(pT.ge.ptmin.and.pt.le.ptmax)then
                !-y/eta cut
                rap=1000.
                if(cut.eq.y)then
                  mT2=pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2
                  if(mT2.gt.0..and.pptl(4,i).gt.0.)then
                    mT=sqrt(mT2)
                    rap=sign(1.,pptl(3,i))
     .                    *alog((pptl(4,i)+abs(pptl(3,i)))/mT)
                  endif
                elseif(cut.eq.eta)then
                  theta=PI/2
                  if(pptl(3,i).ne.0.)then
                    theta=atan(pT/pptl(3,i))
                    if(theta.lt.0.)theta=theta+pi
                    if(theta*0.5.ne.0..and.theta.le.3.141592)then
                      rap=-log(tan(theta*0.5))
                    endif
                  endif
                endif
                if(abs(rap).le.cutmax)then
                  !-Loop over cumulants to calculate them
                  do j=1,order
                    !-Net B
                    if(cumlts(1,j,1).eq.1)then
                      dq=0 !#q-#aq
                      if(abs(idptl(i)).ne.20)then
                        call idqua(i,dq,dS,sumq)
                      else !EXCEPTION for K0s/L
                        dq=0
                      endif
                      !-Ensuring to consider only (anti)baryons
                      if(abs(dq).eq.3)then
                        FDgo=.false.
                        call getiorptl(i,iprnt) !get parent's index
                        if(iprnt.lt.1)then
                          !Automatically consider ptl if parent not known
                          FDgo=.true.
                        else
                          !Ignore ptls who's parent has been checked already
                          !to avoid :
                          ! - x2 couting 
                          ! - kinematic contamination
                          if(istptl(iprnt).eq.ist2)then
                            FDgo=.false.
                          else
                            FDgo=.true.
                          endif
                        endif
                        !-Feed-down correction (DCA) if necessary
                        if(FDC.eq.'DCA')then
                          xptl=xorptl(1,i)
                          yptl=xorptl(2,i)
                          zptl=xorptl(3,i)
                          DCA=sqrt(xptl**2.+yptl**2.+zptl**2.)
                          if(DCA.gt.1e13)FDgo=.false.
                        endif
                        if(FDgo)then 
                          !-Loop over (pseudo)rap interval(s)
                          do k=Ndcut,1,-1
                            if(abs(rap).le.dcut_list(k))then
                              cumlts(2,j,k)=cumlts(2,j,k) + dq/3
                            else
                              exit
                            endif
                          enddo
                        endif
                      endif
                    !-Net Q
                    elseif(cumlts(1,j,1).eq.2)then
                      chrg=0.
                      call idchrg(30,idptl(i),chrg)
                      !-Ensuring to consider only charged ptls
                      if(chrg.ne.0.)then
                        FDgo=.false.
                        call getiorptl(i,iprnt) !get parent's index
                        if(iprnt.lt.1)then
                          !Automatically consider ptl if parent not known
                          FDgo=.true.
                        else
                          !Ignore ptls who's parent has been checked already
                          !to avoid :
                          ! - x2 couting 
                          ! - kinematic contamination
                          if(istptl(iprnt).eq.ist2)then
                            FDgo=.false.
                          else
                            FDgo=.true.
                          endif
                        endif
                        !-Feed-down correction (DCA) if necessary
                        if(FDC.eq.'DCA')then
                          xptl=xorptl(1,i)
                          yptl=xorptl(2,i)
                          zptl=xorptl(3,i)
                          DCA=sqrt(xptl**2.+yptl**2.+zptl**2.)
                          if(DCA.gt.1e13)FDgo=.false.
                        endif
                        if(FDgo)then 
                          !-Loop over (pseudo)rap interval(s)
                          do k=Ndcut,1,-1
                            if(abs(rap).le.dcut_list(k))then
                              cumlts(2,j,k)=cumlts(2,j,k) + chrg
                            else
                              exit
                            endif
                          enddo
                        endif
                      endif
                    !-Net S
                    elseif(cumlts(1,j,1).eq.3)then
                      dS=0 !Initialisation of dS=#s-#as from current particle
                      FDgo=.false.
                      if(abs(idptl(i)).ne.20)then
                        call idqua(i,dq,dS,sumq)
                      else !EXCEPTION for K0s/L
                        dq=0
                        !Generate randomly proprties of K0 or K0b
                        if(rangen().lt.0.5)then
                          dS=1
                        else
                          dS=-1
                        endif
                      endif
                      if(dS.ne.0)then
                        call getiorptl(i,iprnt) !get parent's index
                        if(iprnt.lt.1)then
                          !Automatically consider ptls if parent not known
                          FDgo=.true.
                        else
                          !Ignore ptls who's parent has been checked already
                          !to avoid :
                          ! - x2 couting 
                          ! - kinematic contamination
                          ! - feed-down from c/b (inside kinematic region)
                          if(istptl(iprnt).eq.ist2)then
                            FDgo=.false.
                          else
                            FDgo=.true.
                          endif
                        endif
                        if(FDgo)then 
                          !Loop over different (pseudo)rap intervals
                          do k=Ndcut,1,-1
                            if(abs(rap).le.dcut_list(k))then
                              cumlts(2,j,k)=cumlts(2,j,k) - dS
                            else
                              exit
                            endif
                          enddo
                        endif
                      endif
                    !-Net IDfied particles
                    elseif(cumlts(1,j,1).eq.abs(idptl(i)))then
                      !-Feed-down correction
                      FDgo=.true.
                      if(FDC.eq.'DCA')then
                        xptl=xorptl(1,i)
                        yptl=xorptl(2,i)
                        zptl=xorptl(3,i)
                        DCA=sqrt(xptl**2.+yptl**2.+zptl**2.)
                        if(DCA.gt.1e13)FDgo=.false.
                      elseif(FDC.eq.'FULL')then
                        FDgo=NoLongLivParent(1,i)
                      endif
                      if(FDgo)then
                        !Loop over different (pseudo)rap intervals
                        do k=Ndcut,1,-1
                          if(abs(rap).le.dcut_list(k))then
                            netID=idptl(i)/cumlts(1,j,1)
                            cumlts(2,j,k)=cumlts(2,j,k) + netID
                          else
                            exit
                          endif
                        enddo
                      endif
                    endif !ptls studied
                  enddo !cumlts loop
                endif !rap
              endif !pT
            endif !ity
          endif !ist

        !...
        !   End of particles loop    
        !-------------------------
        enddo

        !-------------------------
        !   Final calculations
        !-------------------------
        !Loop over different (pseudo)rap intervals
        do j=1,Ndcut
          cum=1.0 !Initialised to 1
          !-Loop over different net-numbers
          do i=1,order
            cum=cum*cumlts(2,i,j)
          enddo
          !--- Assigning output variables ---
          !----------------------------------
          ypara(ishift+j,n)=cum
        enddo
      end

c------------------------------------------------------------------------
      subroutine mupurp_sc(n,ishift)  !'mupurp' via SC(n,m) !gs
c------------------------------------------------------------------------
      !--------------------------------------------------------------
      !   Calculation of Symmetric Cumulant or Standard Candles
      !   Reference : arxiv:1604.07663 or 1312.3572v2
      !   input
      !   n = histogram number
      !   xpara(2,n)    order n in SC(n,m)
      !   xpara(3,n)    order m in SC(n,m)
      !   xpara(4,n)    id 
      !   xpara(5,n)    etamin 
      !   xpara(6,n)    etamax 
      !   xpara(7,n)    ptmin 
      !   xpara(8,n)    ptmax
      !   xpara(9,n)    action choice
      !                    1  : (W_(2))_i*<2>_i of Q_n (Eq. (2) )
      !                    2  : (W_(2))_i
      !                    3  : (W_(2))_i*<2>_i of Q_m (Eq. (2) )
      !                    4  : (W_(4))_i
      !                    5  : Numerator Eq. (4) in arxiv:1604.07663
      !--------------------------------------------------------------
      include "epos.inc"
      include "epos.incxan"
      logical go,LongLivPtl,NoLongLivParent

      iordn=nint(xpara(ishift+2,n))
      iordm=nint(xpara(ishift+3,n))
      id=nint(xpara(ishift+4,n))
      etamin=xpara(ishift+5,n)
      etamax=xpara(ishift+6,n)
      ptmin=xpara(ishift+7,n)
      ptmax=xpara(ishift+8,n)
      ichoi=nint(xpara(ishift+9,n))

      we2=0
      we4=0
      oevt=0
      oxnevt=0
      oynevt=0
      oxmevt=0
      oymevt=0
      oxnmevtp=0
      oynmevtp=0
      oxnmevtd=0
      oynmevtd=0

      do iloo=maproj+matarg+1,nptl !loop over ptls
     
        i=iloo   

        go=.false.
        call setIstXY(0,istuse(n),istxxx,istyyy)
        if(istptl(i).eq.istxxx)go=.true.
        if(noweak(n).eq.1)then
          if(LongLivPtl(istyyy,i))go=.true.
        elseif(noweak(n).eq.3)then
          if(LongLivPtl(istxxx,i))go=.false.
        endif
        if(go)go=NoLongLivParent(noweak(n),i)

        if(go)then
          amt=pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2
          px=pptl(1,i)
          py=pptl(2,i)
          pz=pptl(3,i)          
          ptt2=pptl(1,i)**2+pptl(2,i)**2
          pp=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(3,i)**2)
          et=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(5,i)**2)
          if(amt.gt.0..and.pptl(4,i).gt.0.)then
            amt=sqrt(amt)
            rap=sign(1.,pptl(3,i))
     .            *alog( ( pptl(4,i)+abs(pptl(3,i)) ) /amt )
          else
           rap=1000.
          endif            
          arap=abs(rap)
          if(ptt2.gt.0)then
            pt=sqrt(ptt2)
            p3=pptl(3,i)
            eta=sign(1.,p3)*log((sqrt(p3**2+pt**2)+abs(p3))/pt)
          else
            pt=0  
            eta=1000.
          endif
          idepos=ideposf(20,i)
          idabs=abs(idepos)
          call idchrg(28,idepos,ch)
          call idflav(idepos,i1,i2,i3,j4,j)

c           Creation of the Q-vector for each order
          if(eta.gt.etamin.and.eta.lt.etamax)then
          if(pt.gt.ptmin.and.pt.lt.ptmax)then
            if(id.eq.idepos
     .     .or.id.eq.9990.and.i1.eq.0.and.i2.ne.0.and.i3.ne.0
     .     .or.id.eq.9990.and.abs(idepos).eq.20
     .     .or.id.eq.9990.and.i1.ne.0.and.i2.ne.0.and.i3.ne.0
     .     .or.id.eq.9970.and.abs(ch).gt.0.1 )then
              phi=polar(px,py)
              oxnevt=oxnevt+cos(iordn*phi)
              oynevt=oynevt+sin(iordn*phi)
              oxmevt=oxmevt+cos(iordm*phi)
              oymevt=oymevt+sin(iordm*phi)
              oxnmevtp=oxnmevtp+cos((iordn+iordm)*phi)
              oynmevtp=oynmevtp+sin((iordn+iordm)*phi)
              oxnmevtd=oxnmevtd+cos((iordn-iordm)*phi)
              oynmevtd=oynmevtd+sin((iordn-iordm)*phi)
              oevt=oevt+1
            endif
          endif
          endif
            
        endif !loop over ist
      enddo

        onref2=oxnevt**2+oynevt**2   !|Qn|_^2
        omref2=oxmevt**2+oymevt**2   !|Qm|_^2
        onmonom=oxmevt*oxnevt*oxnmevtp-oxnmevtp*oymevt*oynevt
     .          +oxnevt*oymevt*oynmevtp+oxmevt*oynevt*oynmevtp  !Re[Q_{n+m}Qn*Qm*]
        ononmom=oxmevt*oxnevt*oxnmevtd+oxnmevtd*oymevt*oynevt
     .          -oxnevt*oymevt*oynmevtd+oxmevt*oynevt*oynmevtd  !Re[QnQ_{n-m}*Qm*]
        onmref2p=oxnmevtp**2+oynmevtp**2                !|Q_{n+m}|^2
        onmref2d=oxnmevtd**2+oynmevtd**2                !|Q_{n-m}|^2
        
        we2=oevt*(oevt-1)
        we4=oevt*(oevt-1)*(oevt-2)*(oevt-3)

        icas=ishift+1
        
        if(ichoi.eq.1)then            !(W_(2))_i*<2>_i of Q_n
          ypara(icas,n)=onref2-oevt
        elseif(ichoi.eq.2)then        !(W_(2))_i
          ypara(icas,n)=we2
        elseif(ichoi.eq.3)then        !(W_(2))_i*<2>_i of Q_m
          ypara(icas,n)=omref2-oevt
        elseif(ichoi.eq.4)then        !(W_(4))_i
          ypara(icas,n)=we4
        elseif(ichoi.eq.5)then        !Numerator Eq. (4) in arxiv:1604.07663
          ypara(icas,n)=omref2*onref2-2*onmonom-2*ononmom
     .    +onmref2p+onmref2d-(oevt-4)*(onref2+omref2)+oevt*(oevt-6)
        else
        write(ifch,*)'wrong choice of ichoi in mupurp_sc',ichoi
        stop'#####  ERROR 22012018 #####'
        endif

      end

c------------------------------------------------------------------------
      subroutine mupurp_scp(n,ishift)  !'mupurp' v_i via SP       xpara1=1004
c------------------------------------------------------------------------
      !--------------------------------------------------------------
      !   Generalization of Scalar product Method in Standard Variables
      !   input
      !   n = histogram number
      !   xpara(2,n)    order of anisotropy
      !   xpara(3,n)    id a
      !   xpara(4,n)    id c
      !   xpara(5,n)    etamin a
      !   xpara(6,n)    etamax a
      !   xpara(7,n)    etamin c
      !   xpara(8,n)    etamax c
      !   xpara(9,n)    action choice
      !                    1  : mupevt : <q_a> . <q_c>
      !                    2  : mupurp : <<u_b.q_a>>
      !                    3  : mupurp : <<u_b.q_c>>
      ! maybe later
      !   xpara(x,n)   ptmin a,c
      !   xpara(x,n)   ptmax a,c
      !--------------------------------------------------------------
      include "epos.inc"
      include "epos.incxan"
      common/scpvar/xnaevt,ynaevt,snaevt,xncevt,yncevt,sncevt,xnevt
     &,ynevt
      logical go,LongLivPtl,NoLongLivParent

      iord=nint(xpara(ishift+2,n))
      ida=nint(xpara(ishift+3,n))
      idc=nint(xpara(ishift+4,n))
      etaamin=xpara(ishift+5,n)
      etaamax=xpara(ishift+6,n)
      etacmin=xpara(ishift+7,n)
      etacmax=xpara(ishift+8,n)
      ichoi=nint(xpara(ishift+9,n))
      
        
      xnaevt=0
      ynaevt=0
      snaevt=0
      xncevt=0
      yncevt=0
      sncevt=0
      xnevt=0
      ynevt=0

      do iloo=maproj+matarg+1,nptl !loop over ptls
     
        i=iloo   

        go=.false.
        call setIstXY(0,istuse(n),istxxx,istyyy)
        if(istptl(i).eq.istxxx)go=.true.
        if(noweak(n).eq.1)then
          if(LongLivPtl(istyyy,i))go=.true.
        elseif(noweak(n).eq.3)then
          if(LongLivPtl(istxxx,i))go=.false.
        endif
        if(go)go=NoLongLivParent(noweak(n),i)

        if(go)then
          amt=pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2
          px=pptl(1,i)
          py=pptl(2,i)
          pz=pptl(3,i)          
          ptt2=pptl(1,i)**2+pptl(2,i)**2
          pp=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(3,i)**2)
          et=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(5,i)**2)
          if(amt.gt.0..and.pptl(4,i).gt.0.)then
            amt=sqrt(amt)
            rap=sign(1.,pptl(3,i))
     .            *alog( ( pptl(4,i)+abs(pptl(3,i)) ) /amt )
          else
           rap=1000.
          endif            
          arap=abs(rap)
          if(ptt2.gt.0)then
            pt=sqrt(ptt2)
            p3=pptl(3,i)
            eta=sign(1.,p3)*log((sqrt(p3**2+pt**2)+abs(p3))/pt)
          else
            pt=0
            eta=1000.
          endif
          idepos=ideposf(19,i)
          idabs=abs(idepos)
          call idchrg(27,idepos,ch)
          call idflav(idepos,i1,i2,i3,j4,j)

          if(idabs.eq.20.or.idabs.gt.100.and.idabs.lt.10000)then
c         Q vector for set A 
          if(eta.gt.etaamin.and.eta.lt.etaamax)then
            if(ida.eq.idepos
     .      .or.ida.eq.9990.and.i1.eq.0.and.i2.ne.0.and.i3.ne.0
     .      .or.ida.eq.9990.and.abs(idepos).eq.20
     .      .or.ida.eq.9990.and.i1.ne.0.and.i2.ne.0.and.i3.ne.0
     .      .or.ida.eq.9970.and.abs(ch).gt.0.1 )then
              phi=polar(px,py)
              xnaevt=xnaevt+cos(iord*phi)
              ynaevt=ynaevt+sin(iord*phi)
              snaevt=snaevt+1
            endif
          endif

c         Q vector for set C 
          if(eta.gt.etacmin.and.eta.lt.etacmax)then
            if(idc.eq.idepos
     .      .or.idc.eq.9990.and.i1.eq.0.and.i2.ne.0.and.i3.ne.0
     .      .or.idc.eq.9990.and.abs(idepos).eq.20
     .      .or.idc.eq.9990.and.i1.ne.0.and.i2.ne.0.and.i3.ne.0
     .      .or.idc.eq.9970.and.abs(ch).gt.0.1 )then
              phi=polar(px,py)
              xncevt=xncevt+cos(iord*phi)
              yncevt=yncevt+sin(iord*phi)
              sncevt=sncevt+1
            endif
          endif
            
          endif ! loop over idabs
        endif !loop over ist
      enddo
        
      if(snaevt.gt.0)xnaevt=xnaevt/snaevt
      if(snaevt.gt.0)ynaevt=ynaevt/snaevt
      if(sncevt.gt.0)xncevt=xncevt/sncevt
      if(sncevt.gt.0)yncevt=yncevt/sncevt

      if(ichoi.eq.1)ypara(1,n)=xnaevt*xncevt+ynaevt*yncevt !mupevt
      if(ichoi.eq.2)then !mupurp <<ubqa>>
        ypara(ishift+1,n)=xnaevt
        ypara(ishift+2,n)=ynaevt
        ypara(ishift+3,n)=iord
      endif
      if(ichoi.eq.3)then !mupurp <<ubqc>>
        ypara(ishift+1,n)=xncevt
        ypara(ishift+2,n)=yncevt
        ypara(ishift+3,n)=iord
      endif
      
      end

c--------------------------------------------------------------------------------
      subroutine mupurp_cum(n,ishift) !'mupurp' v_i via cumulants   xpara1=1003
c--------------------------------------------------------------------------------
      !--------------------------------------------------------------
      !   input
      !   n = histogram number
      !   xpara(2,n)    order of anisotropy
      !   xpara(3,n)    id ref
      !   xpara(4,n)    id poi
      !   xpara(5,n)    etamin ref
      !   xpara(6,n)    etamax ref
      !   xpara(7,n)    etamin poi
      !   xpara(8,n)    etamax poi
      !   xpara(9,n)    ptmin ref
      !   xpara(10,n)   ptmax ref
      !   xpara(11,n)   ptmin poi
      !   xpara(12,n)   ptmax poi
      !   xpara(13,n)   choice of result imp:cum, pair:weights
      !                    1  : (W_(2))_i * <2>_i
      !                    2  : (W_(2))_i
      !                    3  : (W_(4))_i * <4>_i
      !                    4  : (W_(4))_i
      !                    5  : (W_(6))_i * <6>_i
      !                    6  : (W_(6))_i
      !                    7  : (W_(8))_i * <8>_i
      !                    8  : (W_(8))_i
      !                    9  : (W_(2))_i * <2>_i with gap
      !                    10 : (W_(2))_i with gap
      !                    11 : (w_(2'))_i * <2'>_i
      !                    12 : (w_(2'))_i
      !                    13 : (w_(4'))_i * <4'>_i
      !                    14 : (w_(4'))_i
      !   xpara(14,n)   ity choice
      !                    0 all particle with ist=0
      !                    1 corona : itymin 20 itymax 59
      !                    2 core   : itymin 60 itymax 61
      !   xpara(15,n)   1 : calcul - 0 : just results
      !   xpara(16,n)   choice of observables vs diff flow
      !                          1 vn vs pt, lin scale, 0-20, 100 bins  
      !                          2 vn vs eta
      !                          3 vn vs pt, log scale, 0.5-22, 35 bins
      !                          4 vn vs eta, 20 bins
      !                          5 vn vs pt, log scale, 0.5-22, 18 bins
      !   xpara(17,n)   etamin_gap
      !   xpara(18,n)   etamax_gap
      !--------------------------------------------------------------
      include "epos.inc"
      include "epos.incxan"

      common/cumvar/px2evt(myyarr),py2evt(myyarr),qx2evt(myyarr),
     & qy2evt(myyarr),p2o2(myyarr),o2q2(myyarr),
     & q4o4(myyarr),p2o2o2o2(myyarr),q4o2o2(myyarr),
     & p2o2o4(myyarr),q2o2(myyarr),qx4evt(myyarr),qy4evt(myyarr),
     & cor2p(myyarr),cor4p(myyarr),pevt(myyarr),qevt(myyarr),
     & oevt,ox2evt,oy2evt,o4o2o2,we2,we4,wep2(myyarr),
     & ox4evt,oy4evt,o2ref2,o2ref4,o4ref2,wep4(myyarr),
     & we6,we8,o2ref6,o4o2o2o2o2,ox6evt,oy6evt,ox8evt,oy8evt,
     & o6o2o2o2,o6o4o2,o6ref2,o2ref8,o4o2o2o2o2o2o2,
     & o4o4o2o2o2o2,o6o2o2o2o2o2,o6o2o4o2o2,o8o2o2o2o2,
     & o4o4o4o2o2,o8o4o2o2,o6o2o4o4,o8o6o2,o8o4o4,o8ref2,o4ref4
     &,ox2evta,ox2evtb,oy2evta,oy2evtb,we2gap
      logical go,LongLivPtl,NoLongLivParent

      iord=nint(xpara(ishift+2,n))
      id1=nint(xpara(ishift+3,n))
      id2=nint(xpara(ishift+4,n))
      eta1ref=xpara(ishift+5,n)
      eta2ref=xpara(ishift+6,n)
      eta1poi=xpara(ishift+7,n)
      eta2poi=xpara(ishift+8,n)
      pt1ref=xpara(ishift+9,n)
      pt2ref=xpara(ishift+10,n)
      pt1poi=xpara(ishift+11,n)
      pt2poi=xpara(ishift+12,n)
      ichoi=nint(xpara(ishift+13,n))
      iity=xpara(ishift+14,n)
      icalc=nint(xpara(ishift+15,n))
      imode=xpara(ishift+16,n)
      eta3ref=xpara(ishift+17,n)
      eta4ref=xpara(ishift+18,n)

      nxxx=nhisxxx(n)

      kyyarr=100
      if(imode.eq.3)kyyarr=35
      if(imode.eq.4)kyyarr=20
      if(imode.eq.5)kyyarr=18
      if(kyyarr.gt.myyarr)stop'ERROR 20112020'

      if(icalc.eq.1)then

      do j=1,myyarr
        px2evt(j)=0
        py2evt(j)=0
        qx2evt(j)=0
        qy2evt(j)=0
        qx4evt(j)=0
        qy4evt(j)=0
        p2o2(j)=0
        p2o2o2o2(j)=0
        q4o2o2(j)=0
        p2o2o4(j)=0
        q2o2(j)=0
        o2q2(j)=0
        q4o4(j)=0
        pevt(j)=0
        qevt(j)=0
        wep2(j)=0
        wep4(j)=0
      enddo

        we2=0
        we4=0
        we6=0
        we8=0
        we2gap=0
        oevt=0
        oevta=0
        oevtb=0
        ox2evt=0
        oy2evt=0
        ox4evt=0
        oy4evt=0
        ox6evt=0
        oy6evt=0
        ox8evt=0
        oy8evt=0
        ox2evta=0
        oy2evta=0
        ox2evtb=0
        oy2evtb=0
      
        do iloo=maproj+matarg+1,nptl
     
          i=iloo   
      
          go=.false.
          call setIstXY(0,istuse(n),istxxx,istyyy)
          if(istptl(i).eq.istxxx)go=.true.
          if(noweak(n).eq.1)then
            if(LongLivPtl(istyyy,i))go=.true.
          elseif(noweak(n).eq.3)then
            if(LongLivPtl(istxxx,i))go=.false.
          endif
          if(go)go=NoLongLivParent(noweak(n),i)

          if(go)then
            if(iity.eq.0)then
              itymin=0
              itymax=61
            elseif(iity.eq.1)then
              itymin=20
              itymax=59
            elseif(iity.eq.2)then
              itymin=60
              itymax=61
            else
              write(ifch,*)'wrong choice of iity in mupurp_cum',iity
              stop'#####  ERROR 25082017 #####'
            endif
            if(ityptl(i).ge.itymin.and.ityptl(i).le.itymax)then
            irfp=0
            ipoi=0
            amt=pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2
            px=pptl(1,i)
            py=pptl(2,i)
            pz=pptl(3,i)          
            ptt2=pptl(1,i)**2+pptl(2,i)**2
            pp=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(3,i)**2)
            et=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(5,i)**2)
            if(amt.gt.0..and.pptl(4,i).gt.0.)then
              amt=sqrt(amt)
              rap=sign(1.,pptl(3,i))
     .            *alog( ( pptl(4,i)+abs(pptl(3,i)) ) /amt )
            else
             rap=1000.
            endif            
            arap=abs(rap)
            if(ptt2.gt.0)then
              pt=sqrt(ptt2)
              p3=pptl(3,i)
              eta=sign(1.,p3)*log((sqrt(p3**2+pt**2)+abs(p3))/pt)
            else
              pt=0
              eta=1000.
            endif
            idepos=ideposf(18,i)
            ida=abs(idepos)
            call idchrg(26,idepos,ch)
            call idflav(idepos,i1,i2,i3,j4,j)

c           Reference Flow
            if(eta.gt.eta1ref.and.eta.lt.eta2ref)then
              if(pt.gt.pt1ref.and.pt.lt.pt2ref)then
              if(ida.gt.100.and.ida.lt.10000)then
                if(id1.eq.idepos
     .        .or.id1.eq.9990.and.i1.eq.0.and.i2.ne.0.and.i3.ne.0
     .        .or.id1.eq.9990.and.abs(idepos).eq.20
     .        .or.id1.eq.9990.and.i1.ne.0.and.i2.ne.0.and.i3.ne.0
     .        .or.id1.eq.9970.and.abs(ch).gt.0.1 )then
                  phi=polar(px,py)
                  ox2evt=ox2evt+cos(iord*phi)
                  oy2evt=oy2evt+sin(iord*phi)
                  ox4evt=ox4evt+cos(2*iord*phi)
                  oy4evt=oy4evt+sin(2*iord*phi)
                  ox6evt=ox6evt+cos(3*iord*phi)
                  oy6evt=oy6evt+sin(3*iord*phi)
                  ox8evt=ox8evt+cos(4*iord*phi)
                  oy8evt=oy8evt+sin(4*iord*phi)
                  oevt=oevt+1
                  irfp=1
                endif
                endif
              endif
            endif
            
c           Reference Flow with gap eta 
            if(eta.gt.eta1ref.and.eta.lt.eta2ref)then
              if(pt.gt.pt1ref.and.pt.lt.pt2ref)then
              if(ida.gt.100.and.ida.lt.10000)then
                if(id1.eq.idepos
     .        .or.id1.eq.9990.and.i1.eq.0.and.i2.ne.0.and.i3.ne.0
     .        .or.id1.eq.9990.and.abs(idepos).eq.20
     .        .or.id1.eq.9990.and.i1.ne.0.and.i2.ne.0.and.i3.ne.0
     .        .or.id1.eq.9970.and.abs(ch).gt.0.1 )then
                  phi=polar(px,py)
                 if(eta.gt.eta1ref.and.eta.lt.eta3ref)then
                    ox2evta=ox2evta+cos(iord*phi)
                    oy2evta=oy2evta+sin(iord*phi)
                    oevta=oevta+1
                    irfp=1
                  endif
                  if(eta.gt.eta4ref.and.eta.lt.eta2ref)then
                    ox2evtb=ox2evtb+cos(iord*phi)
                    oy2evtb=oy2evtb+sin(iord*phi)
                    oevtb=oevtb+1
                    irfp=1
                  endif
                endif
                endif
              endif
            endif

c           p vector : all Particles of Interest (POI)
            if(eta.gt.eta1poi.and.eta.lt.eta2poi)then
              if(pt.gt.pt1poi.and.pt.lt.pt2poi)then
              if(ida.gt.100.and.ida.lt.10000)then
                if(id2.eq.idepos
     .          .or.id2.eq.9970.and.abs(ch).gt.0.1 )then
                  if(imode.eq.1)then
                    delpt=20./kyyarr   
                    j=pt/delpt +1
                  elseif(imode.eq.2.or.imode.eq.4)then
                    deleta=(eta2poi-eta1poi)/kyyarr 
                    j=(eta-eta1poi)/deleta +1
                  elseif(imode.eq.3.or.imode.eq.5)then !log scale, 0.5-22
                    delpt=log(22./0.5)/kyyarr   
                    j   = log(pt /0.5)/delpt +1
                  else
                    write(ifch,*)'wrong choice imode mupurp_cum',imode
                    stop'#####  ERROR 18102017 #####'
                  endif
                  phi=polar(px,py)
                  if(j.ge.1.and.j.le.kyyarr)then
                    px2evt(j)=px2evt(j)+cos(iord*phi)
                    py2evt(j)=py2evt(j)+sin(iord*phi)
                    pevt(j)=pevt(j)+1
                  endif
                  ipoi=1
                endif
                endif
              endif
            endif

c          q vector : all particles in Reference AND in POI
            if(irfp.eq.1.and.ipoi.eq.1)then
              if(j.ge.1.and.j.le.kyyarr)then
                qx2evt(j)=qx2evt(j)+cos(iord*phi)
                qy2evt(j)=qy2evt(j)+sin(iord*phi)
                qx4evt(j)=qx4evt(j)+cos(2*iord*phi)
                qy4evt(j)=qy4evt(j)+sin(2*iord*phi)
                qevt(j)=qevt(j)+1
              endif
            endif
            
            endif  !if over ity
          endif    !if over ist
        enddo      !loop over particles

        o2ref2=ox2evt**2+oy2evt**2              !|Q2|^2

        o2ref4=ox2evt**4+2*ox2evt**2*oy2evt**2+oy2evt**4         !|Q2|^4
        o4ref2=ox4evt**2+oy4evt**2              !|Q4|^2
        o4o2o2=ox4evt*ox2evt*ox2evt-ox4evt*oy2evt*oy2evt
     .  + 2*oy4evt*ox2evt*oy2evt !Re[Q4Q2*Q2*]

        o2ref6=ox2evt**6+3*ox2evt**4*oy2evt**2
     .   +3*ox2evt**2*oy2evt**4+oy2evt**6   !|Q2|^6
        o4o2o2o2o2=ox2evt**4*ox4evt-ox4evt*oy2evt**4
     .  +2*ox2evt**3*oy2evt*oy4evt+2*ox2evt*oy2evt**3*oy4evt !Re[Q4Q2Q2*Q2*Q2*]
        o6o2o2o2=ox2evt**3*ox6evt-3*ox2evt*ox6evt*oy2evt**2
     .  +3*ox2evt**2*oy2evt*oy6evt-oy2evt**3*oy6evt !Re[Q6Q2*Q2*Q2*]
        o6o4o2=ox2evt*ox4evt*ox6evt-ox6evt*oy2evt*oy4evt
     .  +ox4evt*oy2evt*oy6evt+ox2evt*oy4evt*oy6evt !Re[Q6Q4*Q2*]
        o6ref2=ox6evt**2+oy6evt**2 !|Q6|^6

        o2ref8=ox2evt**8+4*ox2evt**6*oy2evt**2+6*ox2evt**4*oy2evt**4
     .  +4*ox2evt**2*oy2evt**6+oy2evt**8  !|Q8|^2
        o4o2o2o2o2o2o2=ox2evt**6*ox4evt+ox2evt**4*ox4evt*oy2evt**2
     .  -ox2evt**2*ox4evt*oy2evt**4-ox4evt*oy2evt**6
     .  +2*ox2evt**5*oy2evt*oy4evt+4*ox2evt**3*oy2evt**3*oy4evt
     .  +2*ox2evt*oy2evt**5*oy4evt     !Re[Q4Q2Q2Q2*Q2*Q2*Q2*]
        o4o4o2o2o2o2=ox2evt**4*ox4evt**2
     .  -6*ox2evt**2*ox4evt**2*oy2evt**2+ox4evt**2*oy2evt**4
     .  +8*ox2evt**3*ox4evt*oy2evt*oy4evt
     .  -8*ox2evt*ox4evt*oy2evt**3*oy4evt-ox2evt**4*oy4evt**2
     .  +6*ox2evt**2*oy2evt**2*oy4evt**2-oy2evt**4*oy4evt**2  !Re[Q4Q4Q2*Q2*Q2*Q2*]
        o6o2o2o2o2o2=ox2evt**5*ox6evt-2*ox2evt**3*ox6evt*oy2evt**2
     .  -3*ox2evt*ox6evt*oy2evt**4+3*ox2evt**4*oy2evt*oy6evt
     .  +2*ox2evt**2*oy2evt**3*oy6evt-oy2evt**5*oy6evt   !Re[Q6Q2Q2*Q2*Q2*Q2*]
        o6o2o4o2o2=ox2evt**3*ox4evt*ox6evt
     .  +ox2evt*ox4evt*ox6evt*oy2evt**2
     .  -ox2evt**2*ox6evt*oy2evt*oy4evt
     .  -ox6evt*oy2evt**3*oy4evt+ox2evt**2*ox4evt*oy2evt*oy6evt
     .  +ox4evt*oy2evt**3*oy6evt+ox2evt**3*oy4evt*oy6evt
     .  +ox2evt*oy2evt**2*oy4evt*oy6evt !Re[Q6Q2Q4*Q2*Q2*]
        o8o2o2o2o2=ox2evt**4*ox8evt-6*ox2evt**2*ox8evt*oy2evt**2
     .  +ox8evt*oy2evt**4+4*ox2evt**3*oy2evt*oy8evt
     .  -4*ox2evt*oy2evt**3*oy8evt   !Re[Q8Q2*Q2*Q2*Q2*]
        o4o4o4o2o2=ox2evt**2*ox4evt**3-ox4evt**3*oy2evt**2
     .  +2*ox2evt*ox4evt**2*oy2evt*oy4evt+ox2evt**2*ox4evt*oy4evt**2
     .  -ox4evt*oy2evt**2*oy4evt**2+2*ox2evt*oy2evt*oy4evt**3  !Re[Q4Q4Q4*Q2*Q2*]
        o8o4o2o2=ox2evt**2*ox4evt*ox8evt-ox4evt*ox8evt*oy2evt**2
     .  -2*ox2evt*ox8evt*oy2evt*oy4evt+2*ox2evt*ox4evt*oy2evt*oy8evt
     .  +ox2evt**2*oy4evt*oy8evt-oy2evt**2*oy4evt*oy8evt   !Re[Q8Q4*Q2*Q2*]
        o6o2o4o4=ox2evt*ox4evt**2*ox6evt+2*ox4evt*ox6evt*oy2evt*oy4evt
     .  -ox2evt*ox6evt*oy4evt**2-ox4evt**2*oy2evt*oy6evt
     .  +2*ox2evt*ox4evt*oy4evt*oy6evt+oy2evt*oy4evt**2*oy6evt !Re[Q6Q2Q4*Q4*]
        o8o6o2=ox2evt*ox6evt*ox8evt-ox8evt*oy2evt*oy6evt
     .  +ox6evt*oy2evt*oy8evt+ox2evt*oy6evt*oy8evt  !Re[Q8Q6*Q2*]
        o8o4o4=ox4evt**2*ox8evt-ox8evt*oy4evt**2
     .  +2*ox4evt*oy4evt*oy8evt   !Re[Q8Q4*Q4*]
        o8ref2=ox8evt**2+oy8evt**2 !|Q8^2|
        o4ref4=ox4evt**4+2*ox4evt**2*oy4evt**2+oy4evt**4  !|Q4^4|

        we2=oevt*(oevt-1)
        we4=oevt*(oevt-1)*(oevt-2)*(oevt-3)
        we6=oevt*(oevt-1)*(oevt-2)*(oevt-3)*(oevt-4)*(oevt-5)
        we8=oevt*(oevt-1)*(oevt-2)*(oevt-3)*(oevt-4)*(oevt-5)
     .   *(oevt-6)*(oevt-7)
        we2gap=oevta*oevtb

        do j=1,kyyarr 
          p2o2(j)=px2evt(j)*ox2evt+py2evt(j)*oy2evt        !Re[p2Q2*]
          p2o2o2o2(j)=px2evt(j)*ox2evt**3
     .    + px2evt(j)*oy2evt**2*ox2evt
     .    + 2 * (py2evt(j)*ox2evt**2*oy2evt)
     .    - py2evt(j)*ox2evt**2*oy2evt
     .    + py2evt(j)*oy2evt**3              !Re[p2Q2Q2*Q2*]
          q4o2o2(j)=qx4evt(j)*ox2evt**2
     .    - qx4evt(j)*oy2evt**2
     .    + 2 * ox2evt*oy2evt*qy4evt(j) !Re[q4Q2*Q2*]
          p2o2o4(j)=px2evt(j)*ox2evt*ox4evt
     .    + py2evt(j)*ox2evt*oy4evt
     .    - py2evt(j)*oy2evt*ox4evt
     .    + px2evt(j)*oy2evt*oy4evt                   !Re[p2Q2Q4*]
          q2o2(j)=qx2evt(j)*ox2evt+qy2evt(j)*oy2evt    !Re[q2Q2*]
          o2q2(j)=qx2evt(j)*ox2evt+qy2evt(j)*oy2evt    !Re[Q2q2*]
          q4o4(j)=qx4evt(j)*ox4evt+qy4evt(j)*oy4evt    !Re[q4Q4*]
          wep2(j)=pevt(j)*oevt-qevt(j)
          wep4(j)=(pevt(j)*oevt-3*qevt(j))*(oevt-1)*(oevt-2)
        enddo

      endif  ! of icalc
      
      if(icalc.eq.0.or.icalc.eq.1)then !just result

      icas=ishift+1
      
        if(ichoi.eq.1)then
            ypara(icas,n)=o2ref2-oevt                   !weight*<2>
        elseif(ichoi.eq.3)then
            ypara(icas,n)=o2ref4+o4ref2-2*o4o2o2 
     .      -2*(2*(oevt-2)*o2ref2-oevt*(oevt-3))           !weight*<4>
        elseif(ichoi.eq.5)then
          ypara(icas,n)=o2ref6+9*o4ref2*o2ref2-6*o4o2o2o2o2
     .    +4*(o6o2o2o2-3*o6o4o2)+2*(9*(oevt-4)*o4o2o2+2*o6ref2)
     .    -9*(oevt-4)*(o2ref4+o4ref2)+18*(oevt-2)*(oevt-5)*o2ref2
     .    -oevt*(oevt-4)*(oevt-5)*6           !weight*<6>
        elseif(ichoi.eq.7)then
          ypara(icas,n)=o2ref8-12*o4o2o2o2o2o2o2+6*o4o4o2o2o2o2
     .   +16*o6o2o2o2o2o2-96*o6o2o4o2o2-12*o8o2o2o2o2-36*o4o4o4o2o2
     .   +96*(oevt-6)*o4o2o2o2o2+72*o8o4o2o2+48*o6o2o4o4
     .   +64*(oevt-6)*o6o2o2o2+192*(oevt-6)*o6o4o2-96*o8o6o2
     .   -36*o8o4o4-144*(oevt-7)*(oevt-4)*o4o2o2+36*o4ref2
     .   +64*o6ref2*o2ref2-64*(oevt-6)*o6ref2+9*o4ref4
     .   +36*o2ref4*o4ref2-144*(oevt-6)*o4ref2*o2ref2
     .   +72*(oevt-7)*(oevt-4)*(o4ref4+o2ref4)-16*(oevt-6)*o2ref6
     .   -96*(oevt-7)*(oevt-6)*(oevt-2)*o2ref2
     .   +24*oevt*(oevt-7)*(oevt-6)*(oevt-5)  !weight*<8>
        elseif(ichoi.eq.9)then
          ypara(icas,n)= ox2evta*ox2evtb+oy2evta*oy2evtb  !weight*<2> etagap
        endif

        do j=1,kyyarr 
          yyarr(nxxx,j)=0
          if(ichoi.eq.11)then
              yyarr(nxxx,j)=p2o2(j)-qevt(j)  !weight*<2'>  !KW2011: -pevt(j) was wrong !!   
          elseif(ichoi.eq.13)then
              yyarr(nxxx,j)=(p2o2o2o2(j)-q4o2o2(j)-p2o2o4(j)
     .        - 2*oevt*p2o2(j) - 2*qevt(j)*o2ref2+7*q2o2(j)-o2q2(j)
     .        + q4o4(j) + 2*p2o2(j)+ 2*qevt(j)*oevt - 6*qevt(j))  !correlation 4' particles <4'>
          endif
        
        enddo

        if(ichoi.eq.2)then
          ypara(icas,n)=we2
        elseif(ichoi.eq.4)then
          ypara(icas,n)=we4
        elseif(ichoi.eq.6)then
          ypara(icas,n)=we6
        elseif(ichoi.eq.8)then
          ypara(icas,n)=we8
        elseif(ichoi.eq.10)then
          ypara(icas,n)=we2gap
        elseif(ichoi.eq.12)then
          do j=1,kyyarr
            yyarr(nxxx,j)=wep2(j)
          enddo
        elseif(ichoi.eq.14)then
          do j=1,kyyarr
            yyarr(nxxx,j)=wep4(j)
          enddo
        endif

      endif

      if(ichoi.gt.15)then
        write(ifch,*)'wrong choice of ichoi in mupurp_cum',ichoi
        stop'#####  ERROR 23082017 #####'
      endif



      end

c----------------------------------------------------------------------
      subroutine StandardVariables
c----------------------------------------------------------------------
      include "epos.inc"
      include "epos.incxan"
      logical CDF
      common/cphi2/phi2pos,phi2neg
      double precision sumE,sumP(4),sumAl(2),sumEM(2)

      Emax=0.
      Pmax=0.
      multy1=0
      multc05=0
      multc14=0
      multc24=0
      multc25=0
      multc1=0
      multyi=0
      multc3=0
      multeb=0
      multc83=0
      rapgap=0.
      etamn=-1000.
      etamx=1000.
      avcos2p=0
      avsin2p=0
      avcos2n=0
      avsin2n=0
      sumE=0d0
      sumP(1)=0d0
      sumP(2)=0d0
      sumP(3)=0d0
      sumP(4)=0d0
      rcast(1)=-1d0
      rcast(2)=-1d0
      sumEM(1)=0d0
      sumAl(1)=0d0
      sumEM(2)=0d0
      sumAl(2)=0d0
      imax=0
c--- q-vectors for v2 scalar product method              
      x2aevt=0
      y2aevt=0
      s2aevt=0
      x2cevt=0
      y2cevt=0
      s2cevt=0
c--- qvnsp(): q-vectors for vn scalar product method    
      if(ivnsp.eq.1)then
        do k=1,3
        jvnsp(k)=0 
        do m=1,5
        do n=1,2
        qvnsp(k,m,n)=0 
        qvnsp(k,m,n)=0
        qvnsp(k,m,n)=0
        enddo
        enddo
        enddo
      endif

      do iloo=maproj+matarg+1,nptl
     
        i=iloo   
      
        istxxx=0
        if(istfor.gt.-999)istxxx=istfor
        if(istptl(i).eq.istxxx)then
          ida=abs(idptl(i))
          amt=pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2
          px=pptl(1,i)
          py=pptl(2,i)
          pt=pptl(1,i)**2+pptl(2,i)**2
          pp=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(3,i)**2)
          et=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(5,i)**2)
          if(amt.gt.0..and.pptl(4,i).gt.0.)then
            amt=sqrt(amt)
            rap=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))/amt)
          else
            rap=1000.
          endif
          if(pt.gt.0.)then
            pt=sqrt(pt)
            if(pptl(3,i).ne.0.)then
              theta=atan(pt/pptl(3,i))
            else
              theta=pi/2
            endif
            if(theta.lt.0.)theta=theta+pi
            et=sin(theta)*pptl(4,i)
            if(theta*0.5.eq.0.)then
              eta=1000.
            elseif(theta.gt.3.141592)then ! <= pi anyway
              eta=-1000.
            else
              eta=-log(tan(theta*0.5))
            endif 
          else
            eta=1000.
          endif
          if(  (maproj.gt.2.and.matarg.le.2) 
     .      .or. (maproj.le.2.and.matarg.gt.2)   )then
            eta1=2.0
            eta2=5.0
          else     
            eta1=3.2
            eta2=4.8
          endif
          if(eta.gt.eta1.and.eta.lt.eta2)then
            a=polar(px,py)
            avcos2p=avcos2p+cos(2*a)
            avsin2p=avsin2p+sin(2*a)
          elseif(eta.lt.-eta1.and.eta.gt.-eta2)then
            a=polar(px,py)
            avcos2n=avcos2n+cos(2*a)
            avsin2n=avsin2n+sin(2*a)
          endif
          sumE=sumE+dble(pptl(4,i))
          sumP(1)=sumP(1)+dble(pptl(1,i))
          sumP(2)=sumP(2)+dble(pptl(2,i))
          sumP(3)=sumP(3)+dble(pptl(3,i))
          if(eta.gt.5.1.and.eta.lt.6.55)then
            sumAl(1)=sumAl(1)+dble(pptl(4,i))
            if(ida.eq.10.or.ida.eq.12.or.ida.eq.110.or.ida.eq.220)
     $      sumEM(1)=sumEM(1)+dble(pptl(4,i))
          elseif(eta.lt.-5.1.and.eta.gt.-6.55)then
            sumAl(2)=sumAl(2)+dble(pptl(4,i))
            if(ida.eq.10.or.ida.eq.12.or.ida.eq.110.or.ida.eq.220)
     $      sumEM(2)=sumEM(2)+dble(pptl(4,i))
          endif
c          if(idptl(i).eq.idproj.and.pp.gt.Pmax)then
          if(pp.gt.Pmax)then
            imax=i
            Pmax=pp
          endif
          if(ida.ge.100.and.ida.lt.10000)then
            call idchrg( 18 ,ideposf( 9 ,i),ch)
            CDF=.false.
            if(abs(ch).gt.0.1)then
c---multyi---charged ptl multipl
              multyi=multyi+1
c---multy1---charged ptl multipl for central rap
              if(abs(rap).le.1.)multy1=multy1+1
              if(abs(eta).le.0.5)multc05=multc05+1
              if(abs(eta).le.1.and.pt.gt.0.4)multc14=multc14+1
              if(abs(eta).le.0.8.and.pt.gt.0.3
     *        .and.pt.lt.4.)multc83=multc83+1
              if(abs(eta).lt.2.4)multc24=multc24+1
              if(abs(eta).le.2.5.and.pt.gt.0.5)multc25=multc25+1
              if(abs(eta).le.1)multc1=multc1+1
              if(abs(rap).le.3.15)multc3=multc3+1
c---multeb---charged ptl multipl for back rap
              if(eta.gt.-3.8.and.eta.lt.-2.8)multeb=multeb+1
              if(abs(eta).lt.1.2.and.pt.gt.0.3)then
                CDF=.true.      !CDF CTC acceptance
              elseif(abs(eta).gt.3.2.and.abs(eta).lt.5.9)then
                CDF=.true.      !CDF BBC acceptance
              endif
c---x2aevt,y2aevt,x2cevt,y2cevt:v2 scalar product method              
              a=polar(px,py)
              if(eta.gt.-3.7.and.eta.lt.-1.7)then !sub-evt A = VZERO C
                x2aevt=x2aevt+cos(2*a)
                y2aevt=y2aevt+sin(2*a)
                s2aevt=s2aevt+1
              endif  
              if(eta.gt.2.8.and.eta.lt.5.1)then !sub-evt C = VZERO A
                x2cevt=x2cevt+cos(2*a)
                y2cevt=y2cevt+sin(2*a)
                s2cevt=s2cevt+1
              endif  
c--- qvnsp(): q-vectors for vn scalar product method      
              if(ivnsp.eq.1)then
                if(eta.gt.-3.7.and.eta.lt.-1.7)then !sub-evt A = VZERO C
                  do m=1,5
                    qvnsp(1,m,1) = qvnsp(1,m,1) + cos(m*a)
                    qvnsp(1,m,2) = qvnsp(1,m,2) + sin(m*a)
                    jvnsp(1)     = jvnsp(1) + 1
                  enddo
                endif  
                if(eta.gt.-0.8.and.eta.lt.0.8)then !sub-evt B = ITS,TPC
                  do m=1,5
                    qvnsp(2,m,1) = qvnsp(2,m,1) + cos(m*a)
                    qvnsp(2,m,2) = qvnsp(2,m,2) + sin(m*a)
                    jvnsp(2)     = jvnsp(2) + 1
                  enddo
                endif  
                if(eta.gt.2.8.and.eta.lt.5.1)then !sub-evt C = VZERO A
                  do m=1,5
                    qvnsp(3,m,1) = qvnsp(3,m,1) + cos(m*a)
                    qvnsp(3,m,2) = qvnsp(3,m,2) + sin(m*a)
                    jvnsp(3)     = jvnsp(3) + 1
                  enddo
                endif
              endif
c--- endif charged
            endif
            if(abs(eta).lt.2.4.and.et.gt.0.2)then
              CDF=.true.     !CDF central and plug calorimeters acceptance
            elseif(abs(eta).gt.2.2.and.abs(eta).lt.4.2.and.et.gt.1.)then
              CDF=.true.     !CDF forward calorimeters acceptance
            endif
            if(CDF)then
              if(eta.le.0)etamn=max(etamn,eta)
              if(eta.ge.0)etamx=min(etamx,eta)
            endif
          endif
          if(ilprtg.eq.1)then
            if((((ida.gt.1000.and.ida.lt.10000)
     *         .and.abs(idproj).gt.1000).or.(ida.gt.100
     *         .and.abs(idproj).lt.1000)).and.pptl(4,i)
     *         .gt.Emax.and.pptl(3,i).gt.0.)then
              Emax=pptl(4,i)
              idlead=i
            endif
          else
            if(ida.gt.1000.and.ida.lt.10000
     *        .and.pptl(4,i).gt.Emax.and.pptl(3,i).lt.0.)then
              Emax=pptl(4,i)
              idlead=i
            endif
          endif
        endif
     
      enddo
      
      if(imax.gt.0)then
        sumE=sumE-dble(pptl(4,imax))
        sumP(1)=sumP(1)-dble(pptl(1,imax))
        sumP(2)=sumP(2)-dble(pptl(2,imax))
        sumP(3)=sumP(3)-dble(pptl(3,imax))
      endif
      if(sumAl(1).gt.0d0)rcast(1)=sumEM(1)/sumAl(1)
      if(sumAl(2).gt.0d0)rcast(2)=sumEM(2)/sumAl(2)
      sumP(4)=sqrt(sumP(1)**2+sumP(2)**2+sumP(3)**2)
      xsi=sngl((sumE+sumP(4))*(sumE-sumP(4))/dble(engy)**2)
      if((sumE+sumP(4))*(sumE-sumP(4)).gt.0.)then
        xMdiff=sngl(sqrt((sumE+sumP(4))*(sumE-sumP(4))))
      else
        xMdiff=0
      endif
      rapgap=etamx-etamn
      if(rapgap.gt.100)rapgap=-1.    !not defined
      phi2pos=polar(avcos2p,avsin2p)/2.
      phi2neg=polar(avcos2n,avsin2n)/2.
      if(s2aevt.gt.0.)x2aevt=x2aevt/s2aevt
      if(s2aevt.gt.0.)y2aevt=y2aevt/s2aevt
      if(s2cevt.gt.0.)x2cevt=x2cevt/s2cevt
      if(s2cevt.gt.0.)y2cevt=y2cevt/s2cevt
c--- qvnsp(): q-vectors for vn scalar product method / normalization     
      if(ivnsp.eq.1)then
        do k=1,3
          if(jvnsp(k).gt.0.)then
            do m=1,5
              do n=1,2
                qvnsp(3,m,n) = qvnsp(3,m,n) / jvnsp(k)
              enddo
            enddo
          endif
        enddo
      endif
      end

c----------------------------------------------------------------------
      subroutine jetfind(m,n)
c----------------------------------------------------------------------
c   m = 1 ou 2 (two different definitions)
c   n = histogram
c input(jet definition):
c   xpara(1,n) ... output (m=1) (0=et, 1=pt)
c   xpara(2,n) ... etamin (m=1)
c   xpara(3,n) ... etamax (m=1)
c   xpara(4,n) ... rmax   (m=1) (rmax defining the cone)
c   xpara(5,n) ... ichd   (m=1) (1=charged, 0=all)
c   xpara(6,n) ... output (m=2) (0=et, 1=pt)
c   xpara(7,n) ... etamin (m=2)
c   xpara(8,n) ... etamax (m=2)
c   xpara(9,n) ... rmax   (m=2)
c   xpara(10,n) .. ichd   (m=2)
c output (jet properties):
c   ypara(1,n) ... 1 (found) or 0 if not  (m=1)
c   ypara(2,n) ... et or pt               (m=1)
c   ypara(3,n) ... eta of center          (m=1)
c   ypara(4,n) ... phi of center          (m=1)
c   ypara(5,n)
c   ypara(6,n) ... 1 (found) or 0 if not  (m=2)
c   ypara(7,n) ... et or pt               (m=2)
c   ypara(8,n) ... eta of center          (m=2)
c   ypara(9,n) ... phi of center          (m=2)
c   ypara(10,n)
c----------------------------------------------------------------------

      include "epos.inc"
      include "epos.incxan"
      parameter (mxval=5)
      real ptx(mxval),lst(mxval),etax(mxval),phix(mxval)

      if(m.ne.1.and.m.ne.2)stop'jetfind: value of m not valid.      '

      ipt   = nint(xpara(1+5*(m-1),n))
      etamin=           xpara(2+5*(m-1),n)
      etamax=           xpara(3+5*(m-1),n)
      rmax  =           xpara(4+5*(m-1),n)
      ichd  = nint(xpara(5+5*(m-1),n))

      ifound=0
      do l=1,mxval
        ptx(l)=0
        lst(l)=0
        etax(l)=0
        phix(l)=0
      enddo

ctp060829      pp1=0
ctp060829      pp2=0
ctp060829      pp3=0

      istxxx=0
      if(istu1(istuse(n)).gt.-999)istxxx=istu1(istuse(n))

      do iloo=maproj+matarg+1,nptl
    
        i=iloo   
        
        iok=0
        if(istptl(i).eq.istxxx.and.abs(idptl(i)).lt.10000)iok=1
        if(iok.eq.1)call idchrg( 19 ,ideposf( 10 ,i),ch)
        if(ichd.eq.1.and.nint(ch).eq.0)iok=0
        if(iok.eq.1)then
          p1=pptl(1,i)
          p2=pptl(2,i)
          p3=pptl(3,i)
          pt=sqrt(p1**2+p2**2)
                if(pt.gt.0)then
            eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
            phi=sign(1.,p2)*acos(p1/pt)
          else
            eta=10000
            phi=0
          endif
          do k=1,mxval
            iok=1
            if(m.eq.2)then
              dphi=phi-ypara(4,n)
              if(dphi.lt.-pi)dphi=dphi+2*pi
              if(dphi.gt. pi)dphi=dphi-2*pi
              if(abs(dphi).lt.pi/2)iok=0
            endif
            if(iok.eq.1.and.pt.gt.ptx(k)
     &        .and.eta.le.etamax.and.eta.ge.etamin)then
              do l=mxval,k+1,-1
               ptx(l)=ptx(l-1)
               lst(l)=lst(l-1)
               etax(l)=etax(l-1)
               phix(l)=phix(l-1)
              enddo
               ptx(k)=pt
               lst(k)=i
               etax(k)=eta
               phix(k)=phi
              goto2
            endif
          enddo
  2       continue
        endif
      enddo

      kk=0
      etx=0

      do k=1,mxval
       if(lst(k).ne.0)then

        ifound=1
        et=0
        etaxx=etax(k)
        phixx=phix(k)
        do iloo=maproj+matarg+1,nptl
    
          j=iloo   

          iok=0
          if(istptl(j).eq.istxxx.and.abs(idptl(j)).lt.10000)iok=1
          if(iok.eq.1)call idchrg( 20 ,ideposf( 11 ,j),ch)
          if(ichd.eq.1.and.nint(ch).eq.0)iok=0
          if(iok.eq.1)then
            p1=pptl(1,j)
            p2=pptl(2,j)
            p3=pptl(3,j)
            pt=sqrt(p1**2+p2**2)
            am=pptl(5,j)
                  if(pt.gt.0)then
              eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
              phi=sign(1.,p2)*acos(p1/pt)
            else
              eta=-10000
              phi=0
            endif
            if(eta.le.etamax.and.eta.ge.etamin)then
              deta=eta-etaxx
              dphi=phi-phixx
              if(dphi.lt.-pi)dphi=dphi+2*pi
              if(dphi.gt. pi)dphi=dphi-2*pi
              if(deta**2+dphi**2.lt.rmax**2)then
                if(ipt.eq.0)then           !output is et
                  et=et+sqrt(pt**2+am**2)
                else                       !output is pt
                  et=et+pt
                endif
              endif
            endif
          endif

        enddo

        if(et.gt.etx)then
          etx=et
          kk=k
        endif

       endif
      enddo

      ypara(1+5*(m-1),n)=ifound
      ypara(2+5*(m-1),n)=etx
      if(kk.gt.0)then
       ypara(3+5*(m-1),n)=etax(kk)
       ypara(4+5*(m-1),n)=phix(kk)
      endif
      return
      end

c----------------------------------------------------------------------!bg ga
      subroutine photfill(n)
c---------------------------------------------------------------------------
c   n = histogram (to define istptl=100*n of id=9998 isolated)
c                                            id=9997 non isolated
c input:
c   xpara(1,n) ... 1 : sums energy of partons in the cone . 2 : no parton with E>econe inside the cone
c   xpara(2,n) ... econe     energy max (GeV) in the cone of radius rmax
c   xpara(3,n) ... rmax      radius of the cone
c   xpara(4,n) ... ptmin     for selected photons
c output: 
c   new particles (isolated/non isolated photons) in particle list triggered
c          by : trigger istptl pho pho
c          can be used as usual particle to plot pt, phi, etc ...
c   ityptl = 71, 72 : Direct photons
c                73 : Fragmentation photon produced in ISR
c                74 : Fragmentation photon produced in FSR
c----------------------------------------------------------------------!bg ga
      include "epos.inc"
      include "epos.incxan"
      parameter (mxval=10000)
      parameter (mxxeta=120,mxxphi=63)
      common/cjcheck/jcheck

      econe=xpara(2,n)
      rmax=xpara(3,n)
      ptmin=xpara(4,n)
      istpp=nint(xpara(1,n))+10*nint(10.*econe)+1000*nint(10.*rmax)

      istxxx=0
      if(istu1(istuse(n)).gt.-999)istxxx=istu1(istuse(n))

      nsum=nptl
      do i=1,nsum
        call getistptl(i,istx)  
        call getidptl(i,id)    
        call getityptl(i,ity)
        call getpptl(i,p1,p2,p3,p4,p5) 
        call getdesptl(i,des)
        if(id.eq.10 .and.ity.gt.70 .and.
     &                  istx.eq.istxxx)then !Direct and fragmantation photons
          pt=sqrt(p1**2+p2**2)
          if(pt.ge.ptmin) then
            nptl=nptl+1
            call checkcccptl(nptl)
            call setistptl(nptl,istpp)  
            call setityptl(nptl,ity)
            call setpptl(nptl,p1,p2,p3,p4,p5) 
            call setdesptl(nptl,des)
             !bg idividual normalization for photons
            iso=isolation(i,rmax,econe,n)
            if(iso.eq.0) then
              call setidptl(nptl,9997)     !non isolated photon ID 
            else
              call setidptl(nptl,9998)     !isolated photon ID 
            endif
          endif
        elseif(id.eq.110 .and.istx.eq.istxxx) then   ! final pi0
          pt=sqrt(p1**2+p2**2)
          if(pt.ge.ptmin) then
            nptl=nptl+1
            call checkcccptl(nptl)
            call setistptl(nptl,istpp)  
            call setpptl(nptl,p1,p2,p3,p4,p5) 
            iso=isolation(i,rmax,econe,n)
            if(iso.eq.1) then
              call setidptl(nptl,9996)    !isolated pi0 
            else   
              call setidptl(nptl,9995)    !non isolated pi0 
            endif
          endif
        endif
      enddo

      end

c-----------------------------------------------------------------------!bg ga
      function isolation(n,rmax,econe,nn)
c-----------------------------------------------------------------------------
c     n = particle
c     nn = histogram
c---------------------------------------------------------------------------       
c     isolation criteria for photons
c     R=sqrt(delta(phi)**2+delta(eta)**2)=0.4 not always the same value
c     if a particle is in the cone we add its pt to the variable totpt
c     if totpt>x GeV stop photon not isolated. x depends of the experiment 
c-----------------------------------------------------------------------!bg ga
      include "epos.inc"
      include "epos.incsem"
      include "epos.incxan"

      totpt=0.     
      pt=sqrt(pptl(1,n)**2+pptl(2,n)**2)
      ppp=sqrt(pt**2+pptl(3,n)**2)
      eta=0.5*log((ppp+pptl(3,n))/(ppp-pptl(3,n))) !pseudorapidity
      phi=acos(pptl(1,n)/pt)
      if(pptl(2,n).lt.0.) phi=2.*pi-phi

      do iloo=1,nptl
      
        i=iloo   
      
        itrig=0
        istxxx=0
        if(istu1(istuse(nn)).gt.-999)istxxx=istu1(istuse(nn))
        if(i.ne.n .and.istptl(i).eq.istxxx)then    !bg final hadrons, leptons and photons   
          if(idptl(i).eq.1220) itrig=1                     !Neutron aren't detected
          if(itrig.eq.0) then
            pth=sqrt(pptl(1,i)**2+pptl(2,i)**2)
            ppp=sqrt(pth**2+pptl(3,i)**2)
            pz=pptl(3,i)
            pt=pth
            etah=0.
            if(abs(pz).gt.0..and.pt.gt.0.)etah=sign(1.,pz)*
     *       log((ppp+abs(pz))/pt)
            phih=acos(pptl(1,i)/pth)
            if(pptl(2,i).lt.0.) phih=2.*pi-phih
            delphi=abs(phi-phih)
            if(delphi.gt.pi) delphi=2.*pi-delphi
            deleta=abs(eta-etah)
            if(sqrt(delphi**2+deleta**2).lt.rmax)then
              totpt=totpt+sqrt(pptl(1,i)**2+pptl(2,i)**2) !addition of transverse impulsion of particles in the cone
c                if(ityptl(n).eq.74 .and. pt.ge.10.) then !bg enlever
c                  write(*,*)'222',pt,idptl(i),istptl(i),
c     &       sqrt(delphi**2+deleta**2),sqrt(pptl(1,i)**2+pptl(2,i)**2)
c            endif
              if(totpt.gt.econe) then       
                isolation=0  
                return
              endif
            endif
          endif
        endif

      enddo
      
      isolation=1
      end

c----------------------------------------------------------------------
      subroutine fastjet(n)   !  istfor does not work yet ????????????????????
c----------------------------------------------------------------------
c   n = histogram (to define istptl=100*n of id=9999)
c input:
c   xpara(1,n) ... algorithm 1.0=kt, 0.0=Cam/Aachen,  -1.0 = anti-kt
c   xpara(2,n) ... rapmin    for particles used to define jet
c   xpara(3,n) ... rapmax    for particles used to define jet
c   xpara(4,n) ... rmax      rmax defining the cone
c   xpara(5,n) ... ichd      1=charged, 0=all
c   xpara(6,n) ... ptmin     for particles used to define jet
c   xpara(10,n) .. iorap     0 = rap,  >0 = eta
c   xpara(12,n) .. ptcheck   pt threshhold for checks (default 1e30)
c output: 
c   new particles (jets four momentum) in particle list triggered
c          by : trigger istptl jet jet
c          can be used as usual particle to plot pt, phi, etc ...
c----------------------------------------------------------------------

      include "epos.inc"
      parameter (mxval=10000)
      parameter (mxxeta=120,mxxphi=63)
      double precision p(4,mxval), rmax, algo,jets(4,mxval)
      double precision pav(mxval)
      integer jarray1(mxval),jarray2(mxval)
      integer iphimtx(mxval),ietamtx(mxval),npmtx(mxxeta,mxxphi)
      common/cjcheck/jcheck
      real amtx(mxxeta,mxxphi),bmtx(mxxeta,mxxphi)
      include "epos.incxan"
      data ncntfj/0/
      save ncntfj
      ncntfj=ncntfj+1
      ptcheck=1e30
      jcheck=0
      kcheck=0
 
      algo  = dble(xpara(1,n))
      rapmin=      xpara(2,n)
      rapmax=      xpara(3,n)
      rmax  = dble(xpara(4,n))
      ichd  = nint(xpara(5,n))
      ptmin =      xpara(6,n)
      iorap = nint(xpara(10,n))
      if(xpara(12,n).gt.1e-5)ptcheck=xpara(12,n)
      
      mxphi=0
      mxeta=0
      mxphi=63
      if(mxphi.gt.mxxphi)stop'#####  ERROR 04062012a #####'
      phimx=2*3.1415927
      delphi=phimx/mxphi
      mxeta=iorap
      if(mxeta.gt.mxxeta)stop'#####  ERROR 04062012b #####'
      etamin=rapmin
      etamax=rapmax
      deleta=1
      if(mxeta.gt.0)deleta=(etamax-etamin)/mxeta
      do np=1,mxval
        p(4,np)=0
      enddo
      do ieta=1,mxeta
      do iphi=1,mxphi
        amtx(ieta,iphi)=0
      enddo
      enddo
      do ieta=1,mxeta
      do iphi=1,mxphi
        bmtx(ieta,iphi)=0
      enddo
      enddo
       
      do iloo=maproj+matarg+1,nptl
      
        i=iloo   

        if(istptl(i).eq.25)then
          p1=pptl(1,i)
          p2=pptl(2,i)
          pt2=p2**2+p1**2
          if(pt2.gt.ptcheck**2)then
           jcheck=1
           goto 76
          endif
        endif    

      enddo

  76  continue

      nj=0
      npart=0
      nchi=0

      do iloo=maproj+matarg+1,nptl
     
        i=iloo   

        !print*,'fastjet ',iorptl(i),jorptl(i),'   ',i,'   '
        !.      ,idptl(i),istptl(i),ityptl(i),'   '
        !.      ,sqrt(pptl(1,i)**2+pptl(2,i)**2) 
        if(istptl(i).eq.0.and.abs(idptl(i)).lt.10000)then
          call idchrg( 21 ,ideposf( 13 ,i),ch)
          if(nint(ch).ne.0.or.ichd.eq.0)then
            p1=pptl(1,i)
            p2=pptl(2,i)
            p3=pptl(3,i)
            p4=pptl(4,i)
            p5=pptl(5,i)
            amt=p5**2+p1**2+p2**2
            pt=p2**2+p1**2
            rap=0.
            if(iorap.eq.0)then !rap
              if(amt.gt.0..and.p4+abs(p3).gt.0.)then 
                amt=sqrt(amt)
                rap=sign(1.,p3)*log((p4+abs(p3))/amt) !not correct if particles off-shell
c                  rap=0.5*log((p4+p3)/(p4-p3))  !always correct but numerically unstable
              endif
            elseif(iorap.gt.0)then !eta
              if(sqrt(p3**2+pt)+abs(p3).gt.0..and.pt.gt.0.)then
                pt=sqrt(pt)
                rap=sign(1.,p3)*
     .          alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
              endif
            endif
            pt=sqrt(p2**2+p1**2)
            if(rap.ge.rapmin.and.rap.le.rapmax
     .        .and.pt.gt.ptmin)then  !particle used for jet
              if(iorap.gt.1)then
                eta=rap
                ieta=1+(eta-etamin)/deleta
                phi=polar( p1 , p2 )
                iphi=1+phi/delphi
                if(ieta.ge.1.and.ieta.le.mxeta
     .          .and.iphi.ge.1.and.iphi.le.mxphi)then
                  theta=polar(p3,pt)
                  theta2=2*atan(exp(-eta))
                  et=p4*sin(theta)
                  bmtx(ieta,iphi)=bmtx(ieta,iphi)+et
                  !print*,'++++ theta1/2 et1/2',theta,theta2,et,pt
                  !+++++++++++++++++++++++++check+++++++++++++++++++++++++++++++
                  if(jcheck.eq.1)then
                  if(ityptl(i).eq.30)then
                    iori=iorptl(i)
  77                continue               
                    if(iori.gt.0)then
                      if(istptl(iori).eq.1)then
                        !print*,'        ',i,iori,iorptl(iori)
                        iori=iorptl(iori)
                        goto 77
                      endif
                      if(istptl(iori).ne.29)
     .                 stop'##### ERROR 05062012a #####'
                      i1parton=iorptl(iori)
                      i25=i1parton
  78                  i25=i25-1
                      if(i25.eq.1)stop'##### ERROR 05062012b #####'
                      if(istptl(i25).ne.25)goto 78
                      if(istptl(i25-1).ne.25)
     .                 stop'##### ERROR 05062012c #####'
                      if(iorptl(i25-1).ne.iorptl(i1parton)
     .                 .and.iorptl(i25-1).ne.0)
     .                 stop'##### ERROR 05062012d #####'
                      pt1=sqrt(pptl(1,i25)**2+pptl(2,i25)**2)
                      pt2=sqrt(pptl(1,i25-1)**2+pptl(2,i25-1)**2)
                      if(pt1.gt.ptcheck.or.pt2.gt.ptcheck)then
                        nchi=nchi+1
                        !print*,'++++',i,iori,i1parton,iorptl(i1parton)
                        !.       ,i25,pt1,pt2,eta
                        amtx(ieta,iphi)=amtx(ieta,iphi)+et
                      endif
                    endif
                  endif
                  endif
                  !+++++++++++++++++++++++++check+++++++++++++++++++++++++++++++      
                endif          
              else
                npart=npart+1
                if(npart.gt.mxval)then
                  write(ifmt,*)
     .           'Too many particles (mxval) for Fastjet ! skip ...'
                  return
                endif
                p(1,npart)=dble(p1)
                p(2,npart)=dble(p2)
                p(3,npart)=dble(p3)
                p(4,npart)=dble(p4)
              endif
            endif
          endif
        endif

      enddo
             
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if(iorap.gt.1)then
        npart=0
        do iphi=1,mxphi
        do ieta=1,mxeta
          npmtx(ieta,iphi)=0
          et=bmtx(ieta,iphi)
          if(et.gt.0.001)then
            npart=npart+1
            npmtx(ieta,iphi)=npart
            iphimtx(npart)=iphi
            ietamtx(npart)=ieta
            phi=       (iphi-0.5)*delphi
            eta=etamin+(ieta-0.5)*deleta
            p(1,npart)=et*cos(phi)
            p(2,npart)=et*sin(phi)
            p(3,npart)=et*sinh(eta)
            p(4,npart)=et*cosh(eta)
          endif
        enddo  
        enddo  
      endif

      !++++++++++++++++++++++++++++check++++++++++++++++++++++++++++++++++++
      if(jcheck.eq.1)then
        do ieta=1,mxeta
        do iphi=1,mxphi
          amtx(ieta,iphi)=0
          bmtx(ieta,iphi)=0
        enddo
        enddo
        k25p=0
        ihitp=0
        i25=0
        p1p=0.
        p2p=0.
        p3p=0.
        ptp=0.

        do iloo=maproj+matarg+1,nptl
     
          i=iloo   

          ihit=0 
          k25=0
          if(istptl(i).eq.25)then
            k25=1
            p1=pptl(1,i)
            p2=pptl(2,i)
            p3=pptl(3,i)
            pt=p2**2+p1**2
            if(pt.gt.ptcheck**2.or.ihitp.eq.1)then
              if(k25p.eq.1.and.ihitp.eq.0)then
                nj=nj+1 
                ptp=sqrt(ptp)
                phi=polar( p1p , p2p )
                eta
     .           =sign(1.,p3p)*alog((sqrt(p3p**2+ptp**2)+abs(p3p))/ptp)
                print*,'ist25 partons:   i pt phi eta ='
     .           ,i-1,ptp, phi, eta
                ieta=1+(eta-etamin)/deleta
                iphi=1+phi/delphi
                amtx(ieta,iphi)=amtx(ieta,iphi)+ptp
              endif
              nj=nj+1 
              pt=sqrt(pt)
              phi=polar( p1 , p2 )
              eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
              print*,'ist25 partons:   i pt phi eta =',i,pt, phi, eta
              ieta=1+(eta-etamin)/deleta
              iphi=1+phi/delphi
              amtx(ieta,iphi)=amtx(ieta,iphi)+pt
              ihit=1
              if(k25p.eq.1)then
                i25=i
                ior25=iorptl(i25)
                ii=i25
  79            ii=ii+1
                if(istptl(ii).ne.21.or.iorptl(ii).ne.ior25)
     .          stop'#####  ERROR 10062012  #####'           
                !print*,iorptl(ii),ii,istptl(ii)
                p1=pptl(1,ii)
                p2=pptl(2,ii)
                p3=pptl(3,ii)
                pt=p2**2+p1**2
                pt=sqrt(pt)
                phi=polar( p1 , p2 )
                eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
                ieta=1+(eta-etamin)/deleta
                iphi=1+phi/delphi
                if(ieta.ge.1.and.ieta.le.mxeta
     .          .and.iphi.ge.1.and.iphi.le.mxphi)
     .           bmtx(ieta,iphi)=bmtx(ieta,iphi)+pt
                if(istptl(ii+1).eq.21.and.iorptl(ii+1).eq.ior25)goto 79 
                !print*,iorptl(ii+1),ii+1,istptl(ii+1)
              endif
            endif
            p1p=p1
            p2p=p2
            p3p=p3
            ptp=pt
          endif
          k25p=k25
          ihitp=ihit

        enddo
      endif
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if(npart.le.0)goto 1001  !return
      
      do ieta=1,mxeta
      do iphi=1,mxphi
        amtx(ieta,iphi)=0
        bmtx(ieta,iphi)=0
      enddo
      enddo

      !run the clustering with a pp generalised-kt 
      !sequential recombination alg
      ipri=1
      if(ncntfj.gt.1)ipri=0
      if(ipri.eq.1)write(ifmt,'(a,$)')'fastjet:'
      if(ipri.eq.1)write(ifmt,*)npart,' particles'
      !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
          call fastjetppgenkt(p,npart,rmax,algo,jets,njets,
     &          jarray1,jarray2,ipri) 
      !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      if(njets.le.0)goto 1001 !return
      if(.not.(iorap.gt.1))goto 1002
      
      iconp=0
      do i=1,njets
        iprx=jcheck
        if(i.gt.0)iprx=0 !!!!!!!!!!!!!
        p1=sngl(jets(1,i))  !jet px momentum
        p2=sngl(jets(2,i))  !jet py momentum
        p3=sngl(jets(3,i))  !jet pz momentum
        pt=sqrt(p2**2+p1**2)
        phi=polar( p1 , p2 )
        eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
        if(iprx.eq.1)write(ifmt,'(a,i4,3f7.2,$)')
     .  ' +++jet i pt, phi, eta',i,pt, phi, eta
        if(iprx.eq.1)write(ifmt,'(i7,a,$)')jarray1(i)-iconp,'   '
        p1su=0
        p2su=0
        etsu=0
        etmx=0
        nco=0
        do ij=iconp+1,jarray1(i)
          nco=jarray1(i)-iconp
          nij=jarray2(ij)
          et=sqrt(p(1,nij)**2+p(2,nij)**2)
          etsu=etsu+et
          etmx=max(etmx,et)
          !write(ifmt,'(f6.1,$)')et
          p1su=p1su+p(1,nij)
          p2su=p2su+p(2,nij)
        enddo
        ptxx=sqrt(p1su**2+p2su**2)
        if(abs(pt-ptxx).gt.1e-4)then
          print*,'pt ptxx diff',pt,ptxx,abs(pt-ptxx)
          stop'##### ERROR 10062012b #####'
        endif
        dd=etmx/etsu*nco
        if(iprx.eq.1)write(ifmt,*)' ',dd
        if(dd.gt.5)then
          do ij=iconp+1,jarray1(i)
            nij=jarray2(ij)
            et=sqrt(p(1,nij)**2+p(2,nij)**2)
            iphi=iphimtx(nij)
            ieta=ietamtx(nij)
            if(iprx.eq.1)print*,'+++++',i,ieta,iphi,et
            bmtx(ieta,iphi)=1.
          enddo
        endif
        iconp=jarray1(i)
      enddo
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do ieta=1,mxeta
      avpt=0
      do iphi=1,mxphi
        np=npmtx(ieta,iphi)
        if(np.gt.0)then
          p1=p(1,np)
          p2=p(2,np)
          p3=p(3,np)
          p4=p(4,np)
          pt=sqrt(p2**2+p1**2)
          if(nint(bmtx(ieta,iphi)).ne.1)then
            avpt=avpt+pt
          else
            !if(jcheck.eq.1)write(ifmt,*)'excl',ieta,iphi,pt
          endif
        endif
      enddo  
      avpt=avpt/mxphi
      do iphi=1,mxphi
        np=npmtx(ieta,iphi)
        if(np.gt.0)then
          pav(np)=avpt
          amtx(ieta,iphi)=pav(np)
        endif  
      enddo
      enddo  

      iconp=0
      do i=1,njets
        p1nsu=0
        p2nsu=0
        p3nsu=0
        do ij=iconp+1,jarray1(i)
          nij=jarray2(ij)
          p1=p(1,nij)
          p2=p(2,nij)
          p3=p(3,nij)
          pt=sqrt(p2**2+p1**2)
          ptav=pav(nij)
          ptnew=pt-ptav
          iphi=iphimtx(nij)
          ieta=ietamtx(nij)
          phi=       (iphi-0.5)*delphi
          eta=etamin+(ieta-0.5)*deleta
          p1n=ptnew*cos(phi)
          p2n=ptnew*sin(phi)
          p3n=ptnew*sinh(eta)
          p4n=abs(ptnew)*cosh(eta)
          p1nsu=p1nsu+p1n
          p2nsu=p2nsu+p2n
          p3nsu=p3nsu+p3n
        enddo
        p1=p1nsu
        p2=p2nsu
        p3=p3nsu
        p4=sqrt(p1**2+p2**2+p3**2)
        pt=sqrt(p2**2+p1**2)
        phi=polar( p1 , p2 )
        eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
        imark=0
        if(pt.gt.25)then
          do ij=iconp+1,jarray1(i)
            nij=jarray2(ij)
            iphi=iphimtx(nij)
            ieta=ietamtx(nij)
            bmtx(ieta,iphi)=1.
          enddo
          imark=999999
        endif
        !if(jcheck.eq.1.and.i.le.18)then
        !print*,phi,'  ',eta,'  ',pt,'  ',imark
        !endif
        iconp=jarray1(i)
      enddo  

      do ieta=1,mxeta
      do iphi=1,mxphi
        amtx(ieta,iphi)=0
      enddo
      enddo
      
      do ieta=1,mxeta
      avpt=0
      do iphi=1,mxphi
        np=npmtx(ieta,iphi)
        if(np.gt.0)then
          p1=p(1,np)
          p2=p(2,np)
          p3=p(3,np)
          p4=p(4,np)
          pt=sqrt(p2**2+p1**2)
          if(nint(bmtx(ieta,iphi)).ne.1)then
            avpt=avpt+pt
          else
            !if(jcheck.eq.1)write(ifmt,*)'excl',ieta,iphi,pt
          endif
        endif
      enddo  
      avpt=avpt/mxphi
      do iphi=1,mxphi
        np=npmtx(ieta,iphi)
        if(np.gt.0)then
          pav(np)=avpt
          amtx(ieta,iphi)=pav(np)
        endif  
      enddo
      enddo  

      do ieta=1,mxeta
      do iphi=1,mxphi
        amtx(ieta,iphi)=0
      enddo
      enddo
      
      nptl0=nptl+1
      iconp=0
      do i=1,njets
        p1j=sngl(jets(1,i))  !jet px momentum
        p2j=sngl(jets(2,i))  !jet py momentum
        p3j=sngl(jets(3,i))  !jet pz momentum
        p4j=sngl(jets(4,i))  !jet E  momentum
        ptj=sqrt(p2j**2+p1j**2)
        phij=polar( p1j , p2j )
        etaj=sign(1.,p3j)*alog((sqrt(p3j**2+ptj**2)+abs(p3j))/ptj)
        p1su=0
        p2su=0
        p3su=0
        p1nsu=0
        p2nsu=0
        p3nsu=0
        do ij=iconp+1,jarray1(i)
          nij=jarray2(ij)
          p1=p(1,nij)
          p2=p(2,nij)
          p3=p(3,nij)
          p1su=p1su+p1
          p2su=p2su+p2
          p3su=p3su+p3
          pt=sqrt(p2**2+p1**2)
          ptav=pav(nij)
          ptnew=pt-ptav
          iphi=iphimtx(nij)
          ieta=ietamtx(nij)
          phi=       (iphi-0.5)*delphi
          eta=etamin+(ieta-0.5)*deleta
          p1n=ptnew*cos(phi)
          p2n=ptnew*sin(phi)
          p3n=ptnew*sinh(eta)
          p4n=abs(ptnew)*cosh(eta)
          p1nsu=p1nsu+p1n
          p2nsu=p2nsu+p2n
          p3nsu=p3nsu+p3n
          !if(jcheck.eq.1.alnd.i.le.5)print*,phi,eta,pt,ptnew
        enddo
        p1=p1nsu
        p2=p2nsu
        p3=p3nsu
        p4=sqrt(p1**2+p2**2+p3**2)
        pt=sqrt(p2**2+p1**2)
        nptl=nptl+1
        call checkcccptl(nptl)
        call setistptl(nptl,100*n)     !used for trigger
        call setidptl(nptl,9999)       !jet ID
        call setityptl(nptl,30)
        call setpptl(nptl,p1,p2,p3,p4,sqrt(p1**2+p2**2)) 
        phi=polar( p1 , p2 )
        eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
        ieta=1+(eta-etamin)/deleta
        iphi=1+phi/delphi
        if(jcheck.eq.1.and.i.le.10)then
        !print*,i,p1su-p1j,p2su-p2j,p3su-p3j
        !print*,i,p1j,p1nsu,'  ',p2j,p2nsu,'  ',p3j,p3nsu
        !print*,phij,phi,'  ',etaj,eta,'  ',ptj,pt
        !print*,phi,'  ',eta,'  ',pt
        endif
        if(ieta.ge.1.and.ieta.le.mxeta
     .  .and.iphi.ge.1.and.iphi.le.mxphi)
     .   amtx(ieta,iphi)=amtx(ieta,iphi)+pt
        iconp=jarray1(i)
      enddo  
      
      do ieta=1,mxeta
      do iphi=1,mxphi
        amtx(ieta,iphi)=0
      enddo
      enddo
      
      do np=nptl0,nptl-1
        nmax=np
        do mp=np+1,nptl
          call getpptl(mp,p1mp,p2mp,p3mp,p4mp,p5mp)
          call getpptl(nmax,p1nmax,p2nmax,p3nmax,p4nmax,p5nmax)
          if(p5mp.gt.p5nmax)nmax=mp
        enddo
        if(nmax.gt.np)then 
          call getpptl(nmax,p1xx,p2xx,p3xx,p4xx,p5xx)
          do mp=nmax-1,np,-1
            call getpptl(mp,p1,p2,p3,p4,p5)
            call setpptl(mp+1,p1,p2,p3,p4,p5)
          enddo
          call setpptl(np,p1xx,p2xx,p3xx,p4xx,p5xx)
        endif
        if(jcheck.eq.1.and.np-nptl0+1.le.10)then
          call getpptl(np,p1,p2,p3,dum,pt)
          phi=polar( p1 , p2 )
          eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
          if(np-nptl0+1.le.10)
     .     print*,phi,'  ',eta,'  ',pt
          ieta=1+(eta-etamin)/deleta
          iphi=1+phi/delphi
          if(ieta.ge.1.and.ieta.le.mxeta
     .    .and.iphi.ge.1.and.iphi.le.mxphi)
     .     amtx(ieta,iphi)=amtx(ieta,iphi)+pt
        endif
      enddo
      if(ish.ge.5)call alist('list after fastjet&',nptl0,nptl)
      return

 1002 nptl0=nptl+1
      do i=1,njets
        p1=sngl(jets(1,i))  !jet px momentum
        p2=sngl(jets(2,i))  !jet py momentum
        p3=sngl(jets(3,i))  !jet pz momentum
        p4=sngl(jets(4,i))  !jet E  momentum
        pt=sqrt(p2**2+p1**2)
        nptl=nptl+1
        call checkcccptl(nptl)
        call setistptl(nptl,100*n)    !used for trigger
        call setidptl(nptl,9999)      !jet ID
        call setityptl(nptl,30)
        call setpptl(nptl,p1,p2,p3,p4,sqrt(p1**2+p2**2))  
        phi=polar( p1 , p2 )
        eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
        ieta=1+(eta-etamin)/deleta
        iphi=1+phi/delphi
      enddo  
      if(ish.ge.5)call alist('list after fastjet&',nptl0,nptl)
      return

 1001 continue
      return
      end

c----------------------------------------------------------------------
      subroutine jetevent(n)
c----------------------------------------------------------------------
c uses jets found with fastjet (MANDATORY: activate fastjet BEFORE)
c input:
c      n = histogram
c      xpara(6,n) ... pt_min 
c      xpara(7,n) ... njet        (number of jets required)
c      xpara(8,n) ... delta_phi   (0, not used, >0 uses delta_phi-pi (-0 if <0)
c      xpara(9,n) ... pt_min_first  if >0.001
c      xpara(11,n) .. pt_max_first  if >0.001
c output: 
c      ypara(1,n) ... 1 or 0  (valid event or not)
c      ypara(2,n) ... ajt     (Asymmetry parameter A_j)
c      ypara(3,n) ... fjt     (DeltaPhi)
c      ypara(4,n) ... pjt     (ptSubleading/ptLeading)
c----------------------------------------------------------------------

      include "epos.inc"
      include "epos.incxan"
      dimension inumj(10000)
      common/cncntje/ncntje
      common/cjcheck/jcheck
      ncntje=ncntje+1
      
      ypara(1,n)=0
      njet=nint(xpara(7,n))
      numj=0
      ! count the number of jets with Et>Et_min in the event

      do iloo=maproj+matarg+1,nptl
     
        i=iloo   

        if(istptl(i).eq.100*n)then
          if(pptl(5,i).ge.xpara(6,n))then
            if(numj.eq.0.and.xpara(9,n).gt.0.001)then
              if(pptl(5,i).lt.xpara(9,n))goto 999
            endif
            if(numj.eq.0.and.xpara(11,n).gt.0.001)then
              if(pptl(5,i).gt.xpara(11,n))goto 999
            endif
            numj=numj+1
            if(numj.le.10000)then
            !save position of jet in particle list
              inumj(numj)=i
            else
              write(ifmt,*)
     .        "Too many jets in jetevent, last are skipped!"
            endif
          endif
        endif

      enddo

      !write(ifch,*)"jetevent",numj,abs(xpara(8,n)),inumj(1),inumj(2)
      !if enough jets, analyse them
      if(numj.lt.njet)goto 999
      j=inumj(1)  !fastjet provides jet list ordered in pt
      px=pptl(1,j)
      py=pptl(2,j)
      pt1=pptl(5,j)
      phi1=polar(px,py)
      j=inumj(2)
      px=pptl(1,j)
      py=pptl(2,j)
      pt2=pptl(5,j)
      phi2=polar(px,py)
      if(abs(xpara(8,n)).gt.0.001)then
        phi0=pi
        if(xpara(8,n).lt.0.)phi0=0.
        if(abs(abs(phi1-phi2)-phi0).gt.abs(xpara(8,n)))goto 999
      endif  
      ypara(1,n)=1        !valid event
      ajt=0
      if(pt1+pt2.gt.0.)then
        ajt=(pt1-pt2)/(pt1+pt2)
      endif 
      ypara(2,n)=ajt
      fjt=abs(phi1-phi2)
      if(fjt.gt.3.14159)fjt=2*3.14159-fjt
      ypara(3,n)=fjt
      pjt=0
      if(pt1.gt.0.)pjt=pt2/pt1
      ypara(4,n)=pjt
      iprije=0
      if(ncntje.eq.1.and.iprije.eq.1)
     . write(ifmt,'(a,2f8.2,3x,3f7.2)')
     . ' +++pair pt1 pt2 ajt fjt pjt =',pt1,pt2,ajt,fjt,pjt
      return
 999  continue    
      end

c----------------------------------------------------------------------
      subroutine hardevent(n)
c----------------------------------------------------------------------
c   n = histogram
c input(jet event conditions):
c   xpara(2,n) ... pt1
c   xpara(3,n) ... pt2
c   xpara(4,n) ... absetamax
c   xpara(5,n) ... rmax        (r=sqrt(deltaeta**2+deltaphi**2))
c   xpara(6,n) ... Et_min
c output (jet event found or not):
c   ypara(1,n) ... 1 (found) or 0 if not
c----------------------------------------------------------------------

      include "epos.inc"
      include "epos.incxan"

      ypara(1,n)=0
      
      istxxx=0
      if(istu1(istuse(n)).gt.-999)istxxx=istu1(istuse(n))

      do iloo=maproj+matarg+1,nptl
     
        i=iloo   

       if(abs(idptl(i)).ge.100.and.abs(idptl(i)).lt.10000.
     $  and.istptl(i).eq.istxxx)then
        call idchrg( 22 ,ideposf( 14 ,i),ch)
        if(abs(ch).gt.0.1)then
          p1=pptl(1,i)
          p2=pptl(2,i)
          p3=pptl(3,i)
          pt=sqrt(p1**2+p2**2)
          if(pt.gt.0)then
            eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
            phi=sign(1.,p2)*acos(p1/pt)
          else
            eta=10000
            phi=0
          endif
          if(pt.ge.xpara(2,n).and.abs(eta).lt.xpara(4,n))then
            et1=pptl(4,i)*pt/sqrt(p3**2+pt**2)
            do j=maproj+matarg+1,nptl
              if(j.ne.i
     $        .and.abs(idptl(j)).ge.100.and.abs(idptl(j)).lt.10000.
     $        .and.istptl(j).eq.istxxx)then
                call idchrg( 23 ,ideposf( 15 ,j),ch)
                if(abs(ch).gt.0.1.and.abs(idptl(j)).ge.100
     $         .and.abs(idptl(j)).lt.10000.and.istptl(j).eq.istxxx)then
                  p1=pptl(1,j)
                  p2=pptl(2,j)
                  p3=pptl(3,j)
                  pt=sqrt(p1**2+p2**2)
                        if(pt.gt.0)then
                   etax=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
                   phix=sign(1.,p2)*acos(p1/pt)
                  else
                    etax=-10000
                    phix=0
                  endif
                  if(pt.ge.xpara(3,n).and.abs(etax).lt.xpara(4,n))then
                    deta=eta-etax
                    dphi=phi-phix
                    if(dphi.lt.-pi)dphi=dphi+2*pi
                    if(dphi.gt. pi)dphi=dphi-2*pi
                    if(deta**2+dphi**2.lt.xpara(5,n)**2)then
                    et2=pptl(4,j)*pt/sqrt(p3**2+pt**2)
                     if(et1+et2.gt.xpara(6,n))then
                      ypara(1,n)=1
                      goto1
                     endif
                    endif
                  endif
                endif
              endif
            enddo
          endif
        endif
       endif

      enddo

   1  continue
      return
      end

c----------------------------------------------------------------------!bg
      subroutine corrtrig(n)
c----------------------------------------------------------------------
c   n = histogram
c input(trigger conditions):
c   xpara(1,n) ... mode (0,1)
c   xpara(2,n) ... ptmin
c   xpara(3,n) ... ptmax
c   xpara(4,n) ... etamin (ymin for heavy particle)
c   xpara(5,n) ... etamax (ymax for heavy particle)
c   xpara(6,n) ... type of trigger. 0 : charged hadrons ; else idptl of particle of interest
c   xpara(7,n) ... type of correlation (0,1,2). 0=delphi, 1=xegam without corrections, 2=xegam with underlying event correction
c   xpara(8,n)...  istptl for photons
c output (triggered particle (the one with highest pt if there are several)):
c   ypara(1,n) ... iptl or 0 if no particle found
c   ypara(2,n) ... phi of particle
c   ypara(3,n) ... phi_null
c   ypara(4,n) ... pt lead
c   ypara(5,n) ... = xpara(7,n)
c   ypara(6,n) ... N_trigger
c----------------------------------------------------------------------

      include "epos.inc"
      include "epos.incxan"


      ypara(5,n)= xpara(7,n)
      pt0=xpara(2,n)
      idp=nint(xpara(6,n))

      if(nint(xpara(6,n)).eq.0) then !bg trigger particle type : charged hadrons

      do iloo=maproj+matarg+1,nptl
      
        i=iloo   

        istxxx=0
        if(istu1(istuse(n)).gt.-999)istxxx=istu1(istuse(n))

       if(abs(idptl(i)).ge.100.and.abs(idptl(i)).lt.10000.
     $  and.istptl(i).eq.istxxx)then
        call idchrg( 24 ,ideposf( 16 ,i),ch)
        if(abs(ch).gt.0.1)then
          p1=pptl(1,i)
          p2=pptl(2,i)
          p3=pptl(3,i)
          pt=sqrt(p1**2+p2**2)
          pt=max(pt,1e-20)
          eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
          phi=sign(1.,p2)*acos(p1/pt)
          if(pt.ge.pt0.and.pt.le.xpara(3,n).and.eta.gt.xpara(4,n)
     $      .and.eta.lt.xpara(5,n))then
            pt0=pt
            ypara(1,n)=i
            ypara(2,n)=phi
            ypara(4,n)=pt
          endif
        endif
       endif

      enddo

      else

       do iloo=maproj+matarg+1,nptl
     
        i=iloo   

        if(idptl(i).eq.idp.and.istptl(i).eq.nint(xpara(8,n))) then !bg isolated photons
          p1=pptl(1,i)
          p2=pptl(2,i)
          p3=pptl(3,i)
          pt=sqrt(p1**2+p2**2)
          pt=max(pt,1e-20)
          eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
          phi=sign(1.,p2)*acos(p1/pt)
          if(pt.ge.pt0.and.pt.le.xpara(3,n).and.eta.gt.xpara(4,n)
     $     .and.eta.lt.xpara(5,n))then
            pt0=pt            !bg highest pt of the event
            ypara(1,n)=i
            ypara(2,n)=phi
            ypara(4,n)=pt
          endif
        elseif(idptl(i).eq.idp.and.idptl(i).eq.-140) then !bg final identified particle
          p1=pptl(1,i)
          p2=pptl(2,i)
          p3=pptl(3,i)
          pt=sqrt(p1**2+p2**2)
          pt=max(pt,1e-20)
          eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
          if(pptl(5,i).gt.1.) then
            eta=0.5*log((pptl(4,i)+p3)/(pptl(4,i)-p3))  !bg eta = rapidity for heavy particle
          endif
          phi=sign(1.,p2)*acos(p1/pt)
          if(pt.ge.pt0.and.pt.le.xpara(3,n).and.eta.gt.xpara(4,n)
     $     .and.eta.lt.xpara(5,n))then
            pt0=pt            !bg highest pt of the event
            ypara(1,n)=i
            ypara(2,n)=phi
            ypara(4,n)=pt
          endif
        endif

       enddo

      endif  

      if(ypara(1,n).ne.0) ypara(6,n)=1 !bg
      ypara(3,n)=-0.5
      if(nint(xpara(1,n)).eq.1)ypara(3,n)=0.0
      if(xpara(7,n).ge.1)ypara(3,n)=0.0 !bg for xe

      return
      end

c----------------------------------------------------------------------
      subroutine caltrig(n)
c----------------------------------------------------------------------
c   n = histogram
c input(trigger conditions):
c   xpara(1,n) ... mode (0,1,2) (one eta side, 2 eta side independently or, 
c                                2 eta side simultaneously)
c   xpara(2,n) ... etamin
c   xpara(3,n) ... etamax
c   xpara(4,n) ... 0 all or 1 charged or 2 charged + photons or 3 photon or 4 neutron
c   xpara(5,n) ... ptmin
c   xpara(6,n) ... ptmax
c output (triggered energy):
c   ypara(1,n) ... max E for a particle in eta range (per side)
c----------------------------------------------------------------------

      include "epos.inc"
      include "epos.incxan"
      logical cont
      double precision Eforw,Eback

      mode=nint(xpara(1,n))
      etamin=xpara(2,n)
      etamax=xpara(3,n)
      ptmin=xpara(5,n)
      ptmax=xpara(6,n)
      ichrd=nint(xpara(4,n))
      Eforw=0.d0
      Eback=0.d0
      etaf=1000.
      etab=1000.

      istxxx=0
      if(istu1(istuse(n)).gt.-999)istxxx=istu1(istuse(n))

      do iloo=maproj+matarg+1,nptl
      
        i=iloo   
      
       if(istptl(i).eq.istxxx)then
        call idchrg( 25 ,ideposf( 17 ,i),ch)
        if(ichrd.eq.4.and.abs(ideposf(17,i)).eq.1220)then
          cont=.true.
        elseif(ichrd.eq.3.and.abs(ideposf(17,i)).eq.10)then
          cont=.true.
        elseif(ichrd.ge.1)then
          cont=abs(ch).gt.0.1.or.(ideposf(17,i).eq.10.and.ichrd.eq.2)
        else
          cont=.true.
        endif
        if(cont)then
          p1=pptl(1,i)
          p2=pptl(2,i)
          p3=pptl(3,i)
          pt=sqrt(p1**2+p2**2)
          pt=max(pt,1e-20)
          if(pt.gt.ptmin.and.pt.lt.ptmax)then
            eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
c          phi=sign(1.,p2)*acos(p1/pt)
            if(eta.gt.etamin.and.eta.lt.etamax)
     *           Eforw=max(Eforw,dble(pptl(4,i)))
            if(eta.lt.-etamin.and.eta.gt.-etamax)
     *           Eback=max(Eback,dble(pptl(4,i)))
            if(eta.le.etamax)etaf=min(etaf,etamax-eta)
            if(eta.ge.etamin)etab=min(etab,eta-etamin)
          endif
        endif
       endif
      
      enddo
      
      if(mode.eq.0)then
        ypara(1,n)=sngl(Eforw)
      elseif(mode.eq.1)then
        ypara(1,n)=sngl(max(Eforw,Eback))
      elseif(mode.eq.2)then
        ypara(1,n)=sngl(min(Eforw,Eback))
      endif
      ypara(2,n)=etaf
      ypara(3,n)=etab
      ypara(4,n)=max(etab,etaf)
c      if(typevt.eq.3)print *,mode, ypara(1,n),Eforw,Eback,etaf,etab
c      print *,typevt,mode, ypara(1,n),Eforw,Eback,etaf,etab


      return
      end

c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine xanatest8
c----------------------------------------------------------------------
      include "epos.inc"
      return !comment if needed
      x=fmux(-1.0, 1.00, 0., 1.e10, 1., 1., 0., 0., 0.0, 0.0, 20. ) 
      do i=1,nptl
        ida=abs(idptl(i))
        pt2=pptl(1,i)**2+pptl(2,i)**2
        amt2=pptl(5,i)**2+pt2
        p4=pptl(4,i)
        p3=pptl(3,i)
        if(amt2.gt.0..and.p4+abs(p3).gt.0.d0)then      
          amt=sqrt(amt2)
          ay=log( (p4+abs(p3)) / amt ) 
        else
          ay=1e30
        endif 
        if(ay.lt.0.5.and.pt2.gt.8**2.and.pt2.lt.12**2
     .  .and.(ida.eq.240.or.ida.eq.140.or.ida.eq.241))then
          !if(x.gt.50 .and. x.lt.60)then
          write(ifch,*)'8888888888888',ida,ay,sqrt(pt2),x
          !endif
        endif
      enddo 
      end

