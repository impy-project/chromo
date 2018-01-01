c $Id: tabinit.f,v 1.14 2003/05/02 13:14:46 weber Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine loadwtab (io)
c
c     Revision : 1.0
c
coutput   : information in common-block comwid.f
c
c     load the tabulated branching ratios from disk
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'comres.f'
      include 'comwid.f'

      integer ios, nsp, i, io, ver
      character*35 pwdcmd
      character*10 deftab
      character*8 defexe
      logical b

c set the defaultname of the file, containing the table
      parameter (deftab='tables.dat', defexe='uqmd.exe')

      b=io.eq.1

c get the name of the table from the environment variable
      call getenv('URQMD_TAB',tabname)
c if it is empty, use the default name
      if (tabname(1:4).eq.'    ') then
         tabname=deftab
      endif

      if(b)write (6,*) 'Looking for the tabulated decay width...'
c open the table
      open (unit=75,iostat=ios,file=tabname,form='unformatted',
     .      status='old')
c if it fails ...
      if (ios.ne.0) then
c close the file handle
c         close (unit=75, status='delete')
         if(b)write (6,*) 'No file:',tabname,'in this directory'
c get the full path of the executable, ...
         call getenv('_',pwdcmd)
         write (6,*) 'pwd:',pwdcmd
c extract the path 
         i=max(index(pwdcmd,defexe),2)
         write (tabname,*) pwdcmd (1:i-1),deftab
         tabname=tabname(2:)
         if(b)write (6,*) 'Looking for ',tabname,'...'
c and look for a table in the directory of the executable
         open (unit=75,iostat=ios,file=tabname,
     .         form='unformatted',status='old')        
      endif
c if the last 'open' command succeeds read the file
      if (ios.eq.0) then
         if(b)write (6,*) 'O.K.'
         if(b)write (6,*) 'reading...'
c read all tables 
         read (unit=75, iostat=ios) ver, nsp, tabx, fbtaby, pbtaby,
     .          fmtaby, pmtaby, bwbarnorm, bwmesnorm,
     .      tabxnd, frrtaby
c caution! the file is unformatted, therefor it is system dependent!
         if(b)write (6,*) 'version=',ver
c if no errors occur ...
         if (ios.eq.0) then 
            if(b)write (6,*) 'O.K.'
            wtabflg=3
c check, if the version number is correct
            if (ver.eq.tabver) then
               if(b)write (6,*) 'tabver=',ver,'  O.K.'
            else
               write (6,*) 'wrong table!'
               write (6,*) 'tabver should be',tabver,',instead of',ver
               wtabflg=0
            endif
c check, if the table has the correct 'widnsp'
            if (nsp.eq.widnsp) then
               if(b)write (6,*) 'widnsp=',nsp,'  O.K.'
            else
               write (6,*) 'wrong table!'
               write (6,*) 'widnsp should be',widnsp,', instead of',nsp
                  wtabflg=0
               endif
c if table is O.K. close file
            if (wtabflg.eq.3) then
               close (unit=75, status='keep')
c otherwise ...
            else 
c delete the present table
               close (unit=75, status='delete')
               tabname=deftab
c and calculate a new one
               call mkwtab
            endif
c in case of read errors ...
         else
c delete the present table
            close (unit=75, status='delete')
            write (6,*) 'Error while reading ',tabname
            tabname=deftab
c and calculate a new one
            call mkwtab  
         endif
c in any other case ...
      else
         tabname=deftab
c calculate an new table
         call mkwtab
      endif
      
      return
      end





cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine savewtab
c
c     Revision : 1.0
c
c     save the tabulated branching ratios to disk
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'comres.f'
      include 'comwid.f'

      integer ios

      write (6,*) 'Writing new table...'
c try to generate a new file 
      open (unit=75,iostat=ios,file=tabname,form='unformatted',
     .      status='new')
c if it succedds ...
      if (ios.eq.0) then
c write the tables into the file
         write (unit=75, iostat=ios) tabver, widnsp, tabx, fbtaby, 
     .        pbtaby, fmtaby, pmtaby, bwbarnorm, bwmesnorm,
     .      tabxnd, frrtaby
         if (ios.eq.0) write (6,*) 'O.K.'        
c otherwise complain
      else
         write (6,*) 'Error: ',tabname,'exists!'
      endif
c close the file
      close (unit=75, status='keep')
      
      return
      end





cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkwtab
c
c     Revision : 1.0
c
coutput   : information in common-block comwid.f
c
c     tabulate the mass dependent branching ratios 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'comres.f'
      include 'comwid.f'

      real*8 fwidth,m,first,last,delta,abl0,abln,mir,mminit,fbrancx
      real*8 massit,bran,smass,bwnorm,fppfit
      integer i,bchan,itp,isoit,cmin,cmax,i1,i2,i3,i4,ii1

      write (6,*) 'Generating table...'
c this indicates, that all tables are still empty
      wtabflg=0

c high precision splines from mintab to maxtab1
c lower precicision between maxtab1 and maxtab2

c now fill the x-values
c start with 'mintab'
      first=mintab
c 66 % of all fixpoints between mintab and maxtab1
c calculate the steps
      delta=(maxtab1-mintab)/((widnsp-1d0)*2d0/3d0)
      if (delta.le.0d0) then
         write(*,*)'(E) Please allow maxtab1>mintab in comwid'
         stop 137
      endif
c store the values into 'tabx'
      do 10 i=1,int(widnsp*2./3.)
         m=first+(i-1)*delta
         tabx(i)=m
 10   continue
c 33 % of all fixpoints with larger delta between maxtab1 and maxtab2
      delta=(maxtab2-maxtab1)/((widnsp-1d0)*1d0/3d0)
      if (delta.le.0d0) then
         write(*,*)'(E) Please allow maxtab2>maxtab1 in comwid'
         stop 137
      endif
c store the values into 'tabx'
        do 11 i=int(widnsp*2./3.)+1,widnsp
         m=maxtab1+(i-1-int(widnsp*2./3.))*delta
         tabx(i)=m
 11   continue

c now fill the y-values of the full branching ratios

c these are the first derivatives at the first an the last point
c of the interpolating function. a value greater than 1E30 signals the 
c 'spline' routine to set the boundary condition for a natural spline
c with zero second derivative 
      abl0=2D30
      abln=2D30

c loop over all baryons
      do 20 itp=minbar,maxbar
c loop over all x-values
         do 21 i=1,widnsp
c store the values ...
            fbtaby (i,itp,1)=fwidth(itp,isoit(itp),tabx(i))
 21      continue
c calculate the second derivate and store it in 'fbtaby(,,2)'
         call spline (tabx(1),fbtaby(1,itp,1),widnsp,abl0,abln,
     .                fbtaby(1,itp,2))
 20   continue
      write (6,*) '(1/7) ready.'

c loop over all mesons
      do 30 itp=minmes,maxmes
c loop over all x-values
         do 31 i=1,widnsp
c store the values ...
            fmtaby (i,itp,1)=fwidth(itp,isoit(itp),tabx(i))
 31      continue
c calculate the second derivate and store it in 'fmtaby(,,2)'
         call spline (tabx(1),fmtaby(1,itp,1),widnsp,abl0,abln,
     .                fmtaby(1,itp,2))
 30   continue
      write (6,*) '(2/7) ready.'

c the flag indicates, that now all full widths are tabulated 
      wtabflg=1

c now fill the y-values of the partial branching ratios

c loop over all baryons
      do 40 itp=minbar,maxbar
c get the mass of this particle
         mir=massit(itp)
c get the range of possible decay channels
         call brange (itp, cmin, cmax)
c check, if there are any decay channels         
         if (cmax.gt.0) then
c loop over all decay channels
            do 41 bchan=cmin,cmax
c now get the outgoing particles 'i1' and 'i2' for the channel 'j'
c 'bran' is the mass independent branching ratio (tabulated in blockres)
c 'bflag' indicates, if 'i1', 'i2' or both are broad
               call b3type (itp,bchan,bran,i1,i2,i3,i4)
c check, if decay is allowed

               smass=mminit(i2)
               if(i3.ne.0) smass=smass+mminit(i3)
               if(i4.ne.0) smass=smass+mminit(i4)

               if (bran.gt.1d-9.and.mir.gt.mminit(i1)+smass) then
c loop over all x-values               
                  do 42 i=1,widnsp
c store the values
                     pbtaby(i,1,itp,bchan)=
     .                    fbrancx (bchan,itp,isoit(itp),tabx(i),
     .                    bran,i1,i2,i3,i4)
 42               continue
c calculate the second derivate and store it in 'pbtaby(,2,,)'
                  call spline (tabx(1),pbtaby(1,1,itp,bchan),widnsp,
     .                         abl0,abln,pbtaby(1,2,itp,bchan))
               end if
 41         continue
         end if
 40   continue
      write (6,*) '(3/7) ready.'

c loop over all mesons
      do 50 itp=minmes,maxmes
c get the mass of this particle
         mir=massit(itp)
c get the range of possible decay channels 
         call brange (itp, cmin, cmax)
c check, if there are any decay channels         
         if (cmax.gt.0) then
            do 51 bchan=cmin,cmax
c now get the outgoing particles 'i1' and 'i2' for the channel 'j'
c 'bran' is the mass independent branching ratio (tabulated in blockres)
c 'bflag' indicates, if 'i1', 'i2' or both are broad
               call b3type(itp,bchan,bran,i1,i2,i3,i4)
c!!!
               smass=mminit(i2)
               if(i3.ne.0) smass=smass+mminit(i3)
               if(i4.ne.0) smass=smass+mminit(i4)

               if (bran.gt.1d-9.and.mir.gt.mminit(i1)+smass) then
c loop over all x-values               
                  do 52 i=1,widnsp
                     pmtaby(i,1,itp,bchan)=
     .                    fbrancx (bchan,itp,isoit(itp),tabx(i),
     .                    bran,i1,i2,i3,i4)
 52               continue
c calculate the second derivate and store it in 'pmtaby(,2,,)'
                  call spline (tabx(1),pmtaby(1,1,itp,bchan),widnsp,
     .                         abl0,abln,pmtaby(1,2,itp,bchan))
               end if
 51         continue
         end if
 50   continue

      write (6,*) '(4/7) ready.'


c calculate the norm integral of the Breit-Wigner functions
c   with mass dependent widths

c..baryons
        do 60 i=minbar,maxbar
           bwbarnorm(i)=bwnorm(i)
60      continue
      write (6,*) '(5/7) ready.'

c.. mesons
        do 61 i=minmes,maxmes
           bwmesnorm(i)=bwnorm(i)
61      continue
      write (6,*) '(6/7) ready.'

c now all branching ratios and BW-integrals are tabulated 
      wtabflg=2

ce tabulate fppfit
c fill the x-values
c range of tabulated cross sections
      first=2d0*massit(nucleon)+massit(pimeson)
        last=maxtab1
c calculate the steps
c the energies are weighted quadratically
      delta=(last-first)/((widnsp-1)*2./3.)**2
c store the values into 'tabx'
c 66 % of all fixpoints between mintab and maxtab1
      do 69 i=1,int(widnsp*2./3.)
         m=first+(i-1)**2*delta
         tabxnd(i)=m
 69   continue
c 33 % of all fixpoints with larger, constant delta between maxtab1 and maxtab2
        delta=(maxtab2-last)/((widnsp-1)*1./3.)
        do 70 i=int(widnsp*2./3.)+1,widnsp
         m=maxtab1+(i-1-int(widnsp*2./3.))*delta
         tabxnd(i)=m
 70   continue


c.. all pp-exit channels
c loop over first out-particle N & D
        do 81 ii1=1,2
          if(ii1.eq.1)i1=minnuc
          if(ii1.eq.2)i1=mindel
c loop over second out-particle N(1440)..maxdel
          do 82 i2=minnuc+1,maxdel
c loop over all x-values
          do 83 i=1,widnsp
c store the values ...
              frrtaby(i,1,ii1,i2)=fppfit(99,tabxnd(i),i1,i2)
83          continue
c calculate the second derivate and store it in 'frrtaby(,,2)'
          call spline (tabxnd(1),frrtaby(1,1,ii1,i2),widnsp,abl0,abln,
     .                frrtaby(1,2,ii1,i2))
82        continue
81      continue


      write (6,*) '(7/7) ready.'

c pp cross sections are now tabulated
        wtabflg=3

c save the table on disk
      call savewtab

      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function splint (xa,ya,y2a,n,x)
c
c     Unit     : general infrastructure
c     Author   : (C) Copr. 1986-92 Numerical Recipes Software
c     Date     : 03/07/96
c     Revision : 1.1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'comres.f'
      include 'comwid.f'

      integer n
      integer k,khi,klo
      real*8 x,y,xa(n),y2a(n),ya(n)
      real*8 a,b,h

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
         k=(khi+klo)/2d0
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
       write(6,*) 'bad xa input in splint'
      end if
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a*a*a-a)*y2a(klo)+
     .            (b*b*b-b)*y2a(khi))*(h*h)/6d0
      splint=y

      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function splintth (xa,ya,y2a,n,x,th)
c
c     Unit     : general infrastructure
c     Author   : (C) Copr. 1986-92 Numerical Recipes Software
c                modified my H. Weber
c     Date     : 03/07/96
c     Revision : 1.1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     split routine with nice threshold behaviour for cross sections
c

      implicit none

      include 'comres.f'
      include 'comwid.f'

      integer n
      integer k,khi,klo
      real*8 x,y,xa(n),y2a(n),ya(n)
      real*8 a,b,h,th

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
         k=(khi+klo)/2d0
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
       write(6,*) 'bad xa input in splint'
      end if 
      if (xa(khi).lt.(th+2*h)) then
c linear approximation close to threshold (within 2h)
         splintth=ya(khi)*(x-th)/(xa(khi)-th)
      else
         a=(xa(khi)-x)/h
         b=(x-xa(klo))/h
         y=a*ya(klo)+b*ya(khi)+((a*a*a-a)*y2a(klo)+
     .        (b*b*b-b)*y2a(khi))*(h*h)/6d0
         splintth=y
      endif
      
      return
      end



