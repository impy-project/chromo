c $Id: output.f,v 1.24 2007/05/23 14:28:50 bleicher Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine output(iunit)
c
c     Revision : 1.0
c
c     This subroutine writes the event-header to file(iunit)
C
c
cinput iunit  : output-unit
c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      implicit none

      include 'comres.f'
      include 'coms.f'
      include 'options.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'freezeout.f'
      include 'boxinc.f'

c
      integer iunit,i,ttime,iu,app,att,zpp,ztt,l
      integer iiunit,isunit, id, pdgid
      integer timestep,itotcoll,iinelcoll
      real*8 sigmatot,ptsigtot,stot,otime
      common /outco2/sigmatot


      character*3 abox3, abox4, abox5
      character*5 abox6
      character*4 reffram
      character*20 aa,ah,ai,ak,abox1,abox2
      character*36 ae,abt
      character*31 aee
      character*15 ab,aj,al,am
      character*13 ac,ag,pds,tds
      character*12 ad
      character*7 af
      character*9 ag2
      character*1 add
      character*171 apa14,apa15,apav,line
      character*2 apa,aop
      character*35 aboxhead

c file15out
      integer ind,ind1,ind2,nin
      integer istr,ich,ii,iid

      real*8 sqrts, sigpart, colldens, cdens,cdens_
      logical bdum,paulibl

      include 'outcom.f'
     
      integer fchg,strit
      character*1 echar
c     temporary arrays for CTO and CTP when read in from old event
c     (before they are overwritten)
      integer CTOtmp(numcto)
      real*8 CTPtmp(numctp)
c     CTOtc and CTPtc are the temporary fields for CTOdc and CTPdc.
      character CTOtc(numcto)*2
      character CTPtc(numctp)*2
      integer ctoforty, ctofoone
      integer ctolines,ctplines
      
      integer iou(13:20)

chp variables necessary for hydro evolution output in f15
      real*8 thydro_start,thydro
chp hydro flag for visualization output
      logical hydro_flag
chp new variable for vis-out to count correct npart
      integer npart_form


      save 

      data iou/13,14,15,16,17,18,19,20/

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              output formats
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c fileheader
 101  format(a20,3i7,a15,i2)
 301  format(a13,a13,i4,i4,a12,a13,i4,i4,a1)
c 305  format(a36,3f10.7)
 304  format(a36,3f6.2,a31,1f9.2)
 302  format (a7,i9,a13,i12,a9,a20,i7,a20,f11.3)
c 303  format(a20,i3,a15,e10.4,a15,e10.4,a15,e10.4)
 102  format(a2,15(i3,a2))
 103  format(a2,12(e11.4,a2))
 306  format(a171)

 305  format(a36,3f11.7)
 303  format(a20,i3,a15,e11.4,a15,e11.4,a15,e11.4)

c standard particle information vector
 201  format(9e16.8,i11,2i3,i9,i5,i4)
cLHC 201  format(9e24.16,i11,2i3,i9,i5,i4)

c special output for cto40 (restart of old event)
 210  format(9e16.8,i11,2i3,i9,i5,i10,3e16.8,i8)
cLHC 210  format(9e24.16,i11,2i3,i9,i5,i10,3e24.16,i8)

c special output for mmaker
 203  format(9e16.8,i5,2i3,i6,i5,i4,i5,2e16.8)
cLHC 203  format(9e24.16,i5,2i3,i6,i5,i4,i5,2e24.16)

c same with index for file15
 501  format(i5,9e16.8,i11,2i3,i9,i5,i3,i15)
cLHC 501  format(i5,9e24.16,i11,2i3,i9,i5,i3,i15)

c enhanced file16
 503  format(9e15.7,i11,2i3,i9,i5,i4,2i4)
cLHC 503  format(9e24.16,i11,2i3,i9,i5,i4,2i4)

c same including freeze-out coordinates
 213  format(9e16.8,i11,2i3,i9,i5,i4,8e16.8)
cLHC 213  format(9e24.16,i11,2i3,i9,i5,i4,8e24.16)

c collsision stats for file14
 202  format(8i8)
c same with EndOfEvent tag for file16
 602  format(a1,8i8)

c header-line for each collision in file15
 502  format(i8,i8,i4,i7,f8.3,4e12.4)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c header-line for box-mode
 504  format(a20,e14.6,a20,e14.6,a3,i1,a3,i1,a3,i3)
 505  format(a35)
 506  format(a5,2i4,i8,e14.6)

c
      if(iunit.eq.17)return
      if(bf13.and.(iunit.eq.13)) return
      if(bf14.and.(iunit.eq.14)) return
      if(bf15.and.(iunit.eq.15)) return
      if(bf16.and.(iunit.eq.16)) return

c     copy projectile/target info to local vars
      app=ap
      zpp=zp
      att=at
      ztt=zt

      if(iunit.eq.19) return
c
      aa='UQMD   version:     '
      ab='  output_file '
      abt='transformation betas (NN,lab,pro) '
      ac='projectile:  '
      ad='   target: '
      add=' '
      ae='impact_parameter_real/min/max(fm):  '
      aee='  total_cross_section(mbarn):  '
      af='event# '
      ag=' random seed:' 
      ah='equation_of_state: '
      ai=' total_time(fm/c): '
      aj='  E_lab(GeV/u):'
      ak=' Delta(t)_O(fm/c): '
      al='  sqrt(s)(GeV):'
      am='  p_lab(GeV/u):'
      apa='pa'
      aop='op'
      line=''

      abox1='boxmode length(fm): '
      abox2=' tot. energy (GeV): '
      abox3=' s:'
      abox4=' p:'
      abox5=' #:'
      abox6='box: '
      aboxhead='boxh ityp 2i3       N     pmax(GeV)'
      apa14='pvec: '//
     & 'r0              rx              ry              rz          '//
     & '    p0              px              py              pz      '//
     & '        m          ityp 2i3 chg lcl#  ncl or'    
      apa15='pvec:ind   '//
     & 'r0              rx              ry              rz          '//
     & '    p0              px              py              pz      '//
     & '        m          ityp 2i3 chg lcl#  ncl st'    
      if(iunit.eq.15) then
         apav=apa15
      else
         apav=apa14
      endif

      if(fixedseed) then
         ag2=' (fixed) '
      else
         ag2=' (auto)  '
      endif
      if(prspflg.eq.1) then
         pds='(ityp, char) '
         app=spityp(1)
         zpp=fchg(spiso3(1),app)
      else
         pds='(mass, char) '
      endif
      if(trspflg.eq.1) then
         tds='(ityp, char) '
         att=spityp(2)
         ztt=fchg(spiso3(2),att)
      else
         tds='(mass, char) '
      endif

c determine cross section of the projectile-target system
      sigmatot = ptsigtot()
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      otime=outsteps*dtimestep
      ttime=int(nsteps*dtimestep+0.01)

cbb File 15 has the abbreviated event header unless CTO(58) is 1
      if(iunit.eq.15.and.CTOption(58).ne.1)then
       write(iou(15),502)-1,event,Ap,At,bimp,ecm
     ,     ,sigmatot,ebeam,pbeam
      else
cbb How many cto-lines and ctp-lines are to be written:
cbb write all if extended output is on OR if legacy mode CTO(57) is NOT
cbb set.
      if (CTOption(41).ne.0.or.CTOption(57).ne.0) then
        ctolines = 4  ! * 15 = 60 CTOs
        ctplines = 6  ! * 12 = 72 CTPs
      else
        ctolines = 3  ! * 15 = 45 CTOs
        ctplines = 4  ! * 12 = 48 CTPs
      endif
      write(iou(iunit),101) aa,version, sigver, laires, ab,iunit
      write(iou(iunit),301) ac,pds, App, Zpp, ad,tds, Att, Ztt,add
      write(iou(iunit),305) abt,betann,betatar,betapro
      write(iou(iunit),304) ae,bimp,bmin,bdist,aee,sigmatot
      write(iou(iunit),303) ah,eos,aj,ebeam,al,ecm,am,pbeam
      write(iou(iunit),302) af,event,ag,ranseed,ag2,ai,ttime,ak,otime
      do ii = 0,ctolines-1
       write(iou(iunit),102) aop,(CTOption(i),CTOdc(i)
     &                                       ,i=ii*15+1,ii*15+15)
      enddo
      do ii = 0,ctplines-1
       write(iou(iunit),103) apa,(CTParam(i),CTPdc(i)
     &                                      ,i=ii*12+1,ii*12+12)
      enddo
      if(boxflag.eq.1) then
         write(iou(iunit),504) abox1, lbox, abox2, edens, abox3, solid,
     1                         abox4, para, abox5, mbox
         write(iou(iunit),505) aboxhead
         do 507 l=1,mbox
            write(iou(iunit),506) abox6, bptityp(l), bptiso3(l),
     1                            bptpart(l), bptpmax(l)
 507     continue
      end if
      write(iou(iunit),306) apav
      end if

c 
      return
c.....
      entry uounit(iiunit,isunit)
      iou(iiunit)=isunit
      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file14out(timestep)
c
c     Revision : 1.0
c
c     This subroutine writes the standard output-file (unit 14)
c                    
cinput timestep  : timestep of output
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c

c
      if(bf14)return
      ttime=int(timestep*dtimestep+0.01)
      itotcoll=ctag-dectag
      iinelcoll=itotcoll-NBlColl-NElColl
      write(iou(14),*) npart,ttime
      write(iou(14),202) itotcoll,NElColl,iinelcoll,NBlColl,dectag,
     @     NHardRes,NSoftRes,NDecRes

c now write particle-output

c write spectators
      if(CTOption(28).eq.2)then
         if(CTOption(41).eq.0) then
            do 141 i=1,nspec
               write(iou(14),201) r0s(i),rxs(i),rys(i),rzs(i),p0s(i),
     @              pxs(i),pys(i),pzs(i),sfmass(i),
     @              sityp(i),siso3(i),scharge(i),
     @              -1,-1,0
 141        continue
         else
            do 142 i=1,nspec
               write(iou(14),210) r0s(i),rxs(i),rys(i),rzs(i),p0s(i),
     @              pxs(i),pys(i),pzs(i),sfmass(i),
     @              sityp(i),siso3(i),scharge(i),
     @              -1,-1,0,1d34,0d0,1d0,0
 142        continue
         endif
      endif
      if(CTOption(41).eq.0) then
         do 13 i=1,npart
            write(iou(14),201) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           lstcoll(i),ncoll(i),mod(origin(i),100)
 13      continue
      else
         do 31 i=1,npart
            write(iou(14),210) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           lstcoll(i),ncoll(i),origin(i),
     @           dectime(i),tform(i),xtotfac(i),uid(i)
 31      continue
      endif
c 
      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file13out(timestep)
c
c     Revision : 1.0
c
c     This subroutine writes the standard output-file (unit 13),
c     including the freeze-out configuration of the particles
c                    
cinput timestep  : timestep of output
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c

c
      if(bf13)return
      ttime=int(timestep*dtimestep+0.01)
      itotcoll=ctag-dectag
      iinelcoll=itotcoll-NBlColl-NElColl
      write(iou(13),*) npart,ttime
      write(iou(13),202) itotcoll,NElColl,iinelcoll,NBlColl,dectag,
     @     NHardRes,NSoftRes,NDecRes

c now write particle-output

c write spectators
        if(CTOption(28).eq.2)then
          do 191 i=1,nspec
            write(iou(13),213) r0s(i),rxs(i),rys(i),rzs(i),p0s(i),
     @        pxs(i),pys(i),pzs(i),sfmass(i),
     @        sityp(i),siso3(i),scharge(i),
     @        -1,-1,0,r0s(i),rxs(i),rys(i),rzs(i),p0s(i),
     @        pxs(i),pys(i),pzs(i)
 191      continue
        endif


      do 90 i=1,npart
         if(ncoll(i).eq.0) then
            write(iou(13),213) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           lstcoll(i),ncoll(i),mod(origin(i),100),
     @           r0(i),rx(i),ry(i),rz(i),p0(i),px(i)+ffermpx(i),
     @           py(i)+ffermpy(i),pz(i)+ffermpz(i)
         else
            write(iou(13),213) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           lstcoll(i),ncoll(i),mod(origin(i),100),
     @           frr0(i),frrx(i),frry(i),frrz(i),frp0(i),frpx(i),
     @           frpy(i),frpz(i)
         endif
 90   continue
c 
      return

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      entry file15out(ind1,ind2,sqrts,stot,sigpart)
c
c     Revision : 1.0
c
c     This subroutine writes information about the in-channel to file15
c     (the collision statistics file)
c
cinput        ind1    : index of particle 1
cinput        ind2    : index of particle 2 (=0 for decay of {\tt ind1})
cinput        sqrts   : $\sqrt{s}$ of collision
cinput        stot      : total cross section
cinput        sigpart   : partial cross section
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c determine tag for scatter-input or decay-input
c and store entry channel in temporary observables
      bdum=paulibl(ind1,cdens,-1000)
      tsqrts=sqrts
      tstot=stot
      tsigpart=sigpart

      tind(1)=ind1
      tr0(1)=r0(ind1)
      trx(1)=rx(ind1)
      try(1)=ry(ind1)
      trz(1)=rz(ind1)
      tp0(1)=p0(ind1)
      tpx(1)=px(ind1)
      tpy(1)=py(ind1)
      tpz(1)=pz(ind1)
      tm(1)=fmass(ind1)
      tityp(1)=ityp(ind1)
      tiso3(1)=iso3(ind1)
      tstrange(1) = strit(tityp(1))  
      tcoll(1) = ncoll(ind1)
      tlcoll(1)=lstcoll(ind1)
      torigin(1)=origin(ind1)
      tuid(1)=uid(ind1)
      if(ind2.le.0) then
         nin=1
      elseif(ind2.gt.0) then
         bdum=paulibl(ind2,cdens_,-1000)
         cdens=5d-1*(cdens+cdens_)
         nin=2
         tind(2)=ind2
         tr0(2)=r0(ind2)
         trx(2)=rx(ind2)
         try(2)=ry(ind2)
         trz(2)=rz(ind2)
         tp0(2)=p0(ind2)
         tpx(2)=px(ind2)
         tpy(2)=py(ind2)
         tpz(2)=pz(ind2)
         tm(2)=fmass(ind2)
         tityp(2)=ityp(ind2)
         tiso3(2)=iso3(ind2)
         tstrange(2)=strit(tityp(2))  
         tcoll(2) = ncoll(ind2)
         tlcoll(2)=lstcoll(ind2)
         torigin(2)=origin(ind2)
         tuid(2)=uid(ind2)
      endif

      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry f15outch(colldens)

      if (bf15) return
c     This entry writes information about the collision to file15
c     one line for each particle:
c     format: x y z px py pz ityp iso3...

      write(iou(15),502) nin,nexit,iline,ctag,acttime,tsqrts
     ,     ,tstot,tsigpart,cdens
      do 11 i=1,nin
         istr=strit(tityp(i))
         ich = fchg(tiso3(i),tityp(i))
cbb iid is the particle slot id (standard) or the unique particle id (if
c   CTO(56) == 1).
         iid = tind(i)
         if( CTOption(56).eq.1 ) then
             iid = tuid(i)
         endif
         
         write(iou(15),501) iid,tr0(i),trx(i),try(i),trz(i),
     @                   tp0(i),tpx(i),tpy(i),tpz(i),tm(i),
     @                   tityp(i),tiso3(i),ich,tlcoll(i),
     @                   tcoll(i),istr,torigin(i)
 11   continue
      do 20 ii=1,nexit
         i=inew(ii)
         istr=strit(ityp(i))
         iid = i
         if( CTOption(56).eq.1 ) then
             iid = uid(i)
         endif
         write(iou(15),501) iid,r0(i),rx(i),ry(i),rz(i),
     @                   p0(i),px(i),py(i),pz(i),fmass(i),
     @                   ityp(i),iso3(i),charge(i),lstcoll(i),
     @                   ncoll(i),istr,origin(i) 
 20   continue

      return



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry f15hyin(thydro_start)

      if (bf15) return
c     This entry writes information about the collision to file15
c     one line for each particle:
c     format: x y z px py pz ityp iso3...
      nin=npart
      nexit=0

      write(iou(15),502) nin,nexit,91,0,thydro_start,0.d0
     ,     ,0.d0,0.d0,0.d0
      do 250 i=1,nin
         write(iou(15),501) i,r0(i),rx(i),ry(i),rz(i),
     @                   p0(i),px(i),py(i),pz(i),fmass(i),
     @                   ityp(i),iso3(i),charge(i),lstcoll(i),
     @                   ncoll(i),istr,origin(i) 
 250  continue

      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry f15hyout(thydro_start,thydro)

      if (bf15) return
c     This entry writes information about the collision to file15
c     one line for each particle:
c     format: x y z px py pz ityp iso3...
      nin=0
      nexit=npart

      write(iou(15),502) nin,nexit,96,0,thydro_start+thydro,0.d0
     ,     ,0.d0,0.d0,0.d0
      do 251 i=1,nexit
       if(tform(i).lt.1.d-8)then 
         write(iou(15),501) i,r0(i),rx(i),ry(i),rz(i),
     @                   p0(i),px(i),py(i),pz(i),fmass(i),
     @                   ityp(i),iso3(i),charge(i),lstcoll(i),
     @                   ncoll(i),istr,origin(i)
       else   
         write(iou(15),501) i,tform(i),rx(i),ry(i),rz(i),
     @                   p0(i),px(i),py(i),pz(i),fmass(i),
     @                   ityp(i),iso3(i),charge(i),lstcoll(i),
     @                   ncoll(i),istr,origin(i)
       end if 
 251  continue

      return





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry f15outhy(thydro_start,thydro)

cbb Not sure what hydro_flag does and if it cannot be set somewhere
c   else. Nor what it really does. If it can be set somewhere else, I'd
c   prefer doing that; we could get rid of this entry altogether. For
c   the time being, though, we'll leave it in.
      if(thydro.gt.0.d0)then
       hydro_flag=.true.
      end if
      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry f16outch

      if (bf16.or.(CTOption(13).eq.0)) return
      if (nin.eq.1) then
         tityp(2)=0
      endif
      
      do 22 ii=1,nexit
         i=inew(ii)
         write(iou(16),503) r0(i),rx(i),ry(i),rz(i),
     @        p0(i),px(i)+ffermpx(i),py(i)+ffermpy(i),
     @        pz(i)+ffermpz(i),fmass(i),
     @        ityp(i),iso3(i),charge(i),lstcoll(i),
     @        ncoll(i),mod(origin(i),100),tityp(1),tityp(2)
 22   continue

      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file16out

      echar='E'
      itotcoll=ctag-dectag
      iinelcoll=itotcoll-NBlColl-NElColl

c   
      if(bf16) return

c now write particle-output
      if (CToption(13).eq.0) then
c write spectators
         if (CTOption(28).eq.2)then
            do 151 i=1,nspec
               write(iou(16),201)r0s(i),rxs(i),rys(i),rzs(i),p0s(i),
     @              pxs(i),pys(i),pzs(i),sfmass(i),
     @              sityp(i),siso3(i),scharge(i),
     @              -1,-1,0
 151        continue
         endif

         do 12 i=1,npart
            write(iou(16),201) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           dectag+lstcoll(i),ncoll(i),mod(origin(i),100)

 12      continue
      else
         do 14 i=1,npart
            write(iou(16),503) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           dectag+lstcoll(i),ncoll(i),mod(origin(i),100),-99,-99
 14      continue
      endif
c
c write collision counters etc.
       write(iou(16),602) echar,itotcoll,NElColl,iinelcoll,NBlColl,
     @     dectag,NHardRes,NSoftRes,NDecRes
      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file16entry(ind)
c
c     This entry stores one decay for later output (must be done, in case
c     of pauli-blocked decay)
c
      tr0(3)=r0(ind)
      trx(3)=rx(ind)
      try(3)=ry(ind)
      trz(3)=rz(ind)
      tp0(3)=p0(ind)
      tpx(3)=px(ind)
      tpy(3)=py(ind)
      tpz(3)=pz(ind)
      tm(3)=fmass(ind)
      tityp(3)=ityp(ind)
      tind(3)=ind
      tiso3(3)=iso3(ind)
      tcharge(3)=charge(ind)

c     lstcoll is negative to identify decayed particles 


      tlcoll(3)=-(1*lstcoll(ind))
      tcoll(3)=ncoll(ind)
      torigin(3)=origin(ind)
      tuid(3)=uid(ind)

      return
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file16write

c     This entry writes the decay to file
      i=3

      if(bf16)return
      if (CTOption(13).eq.0) then
       
      write(iou(16),201) tr0(i),trx(i),try(i),trz(i),tp0(i),tpx(i),
     @        tpy(i),tpz(i),tm(i),tityp(i),tiso3(i),tcharge(i),
     @        tlcoll(i),tcoll(i),mod(torigin(i),100)
      else
      write(iou(16),503) tr0(i),trx(i),try(i),trz(i),tp0(i),tpx(i),
     @        tpy(i),tpz(i),tm(i),tityp(i),tiso3(i),tcharge(i),
     @        tlcoll(i),tcoll(i),mod(torigin(i),100),-98,-98
      endif

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc_header

      if (bf19) return

chp initialize flag
      hydro_flag=.false.
      
      write (19,901) 'OSC1997A    '
      write (19,901) 'final_id_p_x'

 901  format (a12)

      if (CTOption(27).eq.0) then
         reffram='eqsp'
      elseif (CTOption(27).eq.1) then
         reffram='tar'
      elseif (CTOption(27).eq.2) then
         reffram='pro'
      else
         call error ('osc_header','Unknown Ref-Frame',
     .        dble(CTOption(27)),2)
         reffram='----'
      endif

      write (19,902) 'UrQMD', versiontxt, app, zpp, att, ztt,
     .     reffram, ebeam, 1

 902  format (2(a8,2x),'(',i3,',',i6,')+(',i3,',',i6,')',2x,a4,2x,
     &     e10.4,2x,i8)

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc99_header

c header for OSCAR 99A output format

      if (bf20) return
      
      write (20,991)
      write (20,992)

 991  format ('# OSC1999A')
 992  format ('# full_event_history')

      if (CTOption(27).eq.0) then
         reffram='nncm'
      elseif (CTOption(27).eq.1) then
         reffram='tar'
      elseif (CTOption(27).eq.2) then
         reffram='pro'
      else
         call error ('osc_header','Unknown Ref-Frame',
     .        dble(CTOption(27)),2)
         reffram='----'
      endif

      write (20,993) versiontxt
 993  format ('# UrQMD ',a8)

      write (20,994) app, zpp, att, ztt,reffram, ebeam, 1

 994  format ('# (',i3,',',i6,')+(',i3,',',i6,')',2x,a4,2x,
     &     e10.4,2x,i8)

      return



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc_event

c body for OSCAR 97A format

      if (bf19) return

   
   
      write (19,903) event, npart, bimp, 0D0


 903  format (i10,2x,i10,2x,f8.3,2x,f8.3)


 904  format (i10,2x,i10,2x,9(e12.6,2x))
 

c particles, original

      do 99 i=1,npart
         id = pdgid(ityp(i), iso3(i))
         write(19,904) i, id, 
     .        px(i)+ffermpx(i), py(i)+ffermpy(i), pz(i)+ffermpz(i), 
     .        p0(i), fmass(i),     
     .        frrx(i), frry(i), frrz(i), frr0(i)
 99   continue


      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc_vis(timestep)

c body for OSCAR 97A format adjusted for visualization

      if (bf19) return

      npart_form=0
      
       if(hydro_flag.eqv..true.)then 
        do i=1,npart
         if(tform(i).le.acttime) then
           npart_form=npart_form+1
         end if
        end do 
       end if

       if(hydro_flag.eqv..true.) then   
           write (19,2301) event, npart_form, bimp, 0D0,nsteps,timestep
       else 
           write(19,2301) event, npart, bimp,0D0,nsteps,timestep
       end if  

 2301       FORMAT(I10,2X,I10,2X,F8.3,2X,F8.3,2x,i4,2x,i4)


chp modification for visualization output
      if(hydro_flag.eqv..true.)then
       do 990 i=1,npart
         id = pdgid(ityp(i), iso3(i))
         if(tform(i).le.acttime)then
         write(19,2302) i, id, 
     .       px(i)+ffermpx(i), py(i)+ffermpy(i), pz(i)+ffermpz(i), 
     .       p0(i), fmass(i),     
     .       rx(i), ry(i), rz(i), r0(i),tform(i),
     .       frrx(i), frry(i), frrz(i), frr0(i),ncoll(i)
         end if   
 990    continue
      else 
       do 980 i=1,npart
         id = pdgid(ityp(i), iso3(i))
         write(19,2302) i, id, 
     .       px(i)+ffermpx(i), py(i)+ffermpy(i), pz(i)+ffermpz(i), 
     .       p0(i), fmass(i),     
     .       rx(i), ry(i), rz(i), r0(i),tform(i),
     .       frrx(i), frry(i), frrz(i), frr0(i),ncoll(i) 
 980    continue
      end if
     
 2302  FORMAT(I10,2X,I10,14(2X,E12.6),I10)
 


      return


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc99_event(ind)

c full event info for OSCAR 99A format

      if (bf20) return

      if(ind.eq.-1) then
         write (20,995) 0, npart, event, bimp, 0D0
      elseif(ind.eq.1) then
         write (20,996) npart, 0
      else
         write(6,*) 'fatal error in osc_99_event: wrong tag'
         stop 137
      endif

 995  format (3(i7,2x),2(f8.3,2x))
 996  format (2(i7,2x))

c particles

      do 88 i=1,npart
         id = pdgid(ityp(i), iso3(i))
         write(20,997) uid(i), id, 0,  
     .        px(i)+ffermpx(i), py(i)+ffermpy(i), pz(i)+ffermpz(i), 
     .        p0(i), fmass(i),     
     .        frrx(i), frry(i), frrz(i), frr0(i)
 88   continue

 997  format (3(i10,2x),9(e12.6,2x))

      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc99_coll

      if (bf20) return
c     This entry writes information about the collision to file20
c     one line for each particle:
c     format: x y z px py pz ityp iso3...

      write(iou(20),999) nin,nexit,iline,ctag,acttime,tsqrts
     ,     ,tstot,tsigpart,cdens

      do 911 i=1,nin
         id = pdgid(tityp(i), tiso3(i))
         write(20,997) tuid(i), id, 0,  
     .        tpx(i), tpy(i), tpz(i),tp0(i),tm(i), 
     .        trx(i), try(i), trz(i), tr0(i)
 911   continue
      do 912 ii=1,nexit
         i=inew(ii)
         id = pdgid(ityp(i), iso3(i))
         write(20,997) uid(i), id, 0,  
     .        px(i), py(i), pz(i),p0(i),fmass(i), 
     .        rx(i), ry(i), rz(i), r0(i)
 912  continue


 999  format(3(i7,2x),i7,f8.3,4e12.4)
      return


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc99_eoe

c end of event tag for OSCAR 99A format
      if (bf20) return

      write(20,996) 0,0

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry getoldevent

c     read event header
      read(10,*,end=666) 
     @ aa

      read(10,301) ac,pds, App, Zpp, ad,tds, Att, Ztt,add       
      read(10,305) abt,betann,betatar,betapro
      read(10,304) ae,bimp,bmin,bdist,aee,sigmatot
      read(10,303) ah,eos,aj,ebeam,al,ecm,am,pbeam
      read(10,302) af,event,ag,ranseed,ag2,ai,ttime,ak,otime
C read unspecified number of lines with CTOptions (no more lines than
C possible; since 15 CTOs per line and no more than numcto options,
C there cannot be more than numcto/15+1 lines (note that numcto/15 is a
C truncated integer):
      do ii=0, numcto/15
C get next line:
        read(10,'(A)') line
C if that doesn't start with "op", we have reached the first
C non-CTO-line
        if(line(1:2).ne.'op') then
C ... and need to exit the loop.
          exit
        endif
        ! (else we are where we want to be and read Options)
        read(line,102) aop,(CTOtmp(i),CTOtc(i),i=ii*15+1,ii*15+15)
      enddo
C Same game for CTParams. Here, we inherit one line from the CTO code
C above.
      do ii=0, numctp/12
        if(line(1:2).ne.'pa') then
C doesn't start with 'pa'? Get outta here!
          exit
        endif
        ! (else we are where we want to be and read Parameters)
        read(line,103) apa,(CTPtmp(i),CTPtc(i),i=ii*12+1,ii*12+12)
C read next line for next iteration of loop.
        read(10,'(A)') line
      enddo
C if we got out of the loop, we might have a box header at hand.
      if (line(1:3).eq.'box') then
        boxflag = 1
        ! box header with length, energy, flags and number of particle
        ! specification lines to follow.
        read(line,504) abox1, lbox, abox2, edens, abox3, solid, abox4,
     1                 para, abox5, mbox
        read(10,505) aboxhead
        do l=1, mbox
         ! read from line
         read(10,506) abox6, bptityp(l), bptiso3(l),
     1                       bptpart(l), bptpmax(l)
        enddo
        ! read next line (which we would have if we didn't use up all
        ! our mojo for the box headers):
        read(10,'(A)') line
      endif
C Now, we expect the next line to be the "Particle Vector" line.
      read(line,306) apav

C set CTO only if CTOption(40) is 1 (standard behaviour) or if the
C asterisk '*' has been set (and, implicitly, CTOption(40) > 1).
C CTOption(40) will be overwritten in this loop, so we save the original
C value.
C Also, we set CTOdc(i) only if the above is true: We want to keep
C CTOdc(i) set from the inputfile.
      ctoforty=CTOption(40)
      ctofoone=CTOption(41)
      do i=1, numcto
        if( (ctoforty.eq.1.and.CTOption(i).ne.CTOtmp(i))
     &   .or.CTOtc(i).eq.' *') then
          write(*,*) 'Setting Option ',i,' (default: ',CTOption(i),
     &               ') to ',CTOtmp(i)
          CTOption(i) = CTOtmp(i)
          CTOdc(i) = CTOtc(i)
        endif
      enddo
C same as above, but for CTParam and CTPdc.
      do i=1, numctp
        if( (ctoforty.eq.1.and.CTParam(i).ne.CTPtmp(i))
     &   .or.CTPtc(i).eq.' *') then
          write(*,*) 'Setting Parameter ',i,' (default: ',CTParam(i),
     &               ') to ',CTPtmp(i)
          CTParam(i) = CTPtmp(i)
          CTPdc(i) = CTPtc(i)
        endif
      enddo
c check that CTO(41) is turned on. If not, issue a warning.
      if (CTOption(41).eq.0) then
        write(6,*) 'Warning: Reading in old event which seems to be ',
     +             'generated with CTO(41) = 0!'
      endif
c reset option 40
      CTOption(40)=ctoforty
c also reset option 41: Of course, it is turned on in the output we pipe
c in (or, really, we might try to take a look at that and warn if it
c isn't), but that shouldn't determine if we have the output of the
c continuation still in extended format or not; the current input file
c should do that instead.
      CTOption(41)=ctofoone

c read event body
      read(10,*) npart,ttime
      read(10,202) itotcoll,NElColl,iinelcoll,NBlColl,dectag,
     @     NHardRes,NSoftRes,NDecRes
c      timestep=dble(ttime)/dtimestep
      ctag=itotcoll+dectag
c now read particle-output
      nbar=0
      do 39 i=1,npart
         read(10,210) r0(i),rx(i),ry(i),rz(i),p0(i),
     @        px(i),py(i),pz(i),fmass(i),
     @        ityp(i),iso3(i),charge(i),
     @        lstcoll(i),ncoll(i),origin(i),
     @        dectime(i),tform(i),xtotfac(i)
      if(abs(ityp(i)).le.maxbar)nbar=nbar+1
 39   continue
      nmes=npart-nbar
      acttime=r0(1) 
c     read options-file
      call getparams
      success = .true.
      return
c handle EOF. This is not necessarily an error.
 666  write(*,*) 'No more events to read'
      success = .false.
      return


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry specout(ind,iu)
      i=ind
      if (CTOption(28).lt.0) return 
      if (iu.eq.16.and.bf16) return
      if (iu.eq.14.and.bf14) return
      write(iu,201) r0(i),rx(i),ry(i),rz(i),p0(i),
     @     px(i)+ffermpx(i),py(i)+ffermpy(i),
     @     pz(i)+ffermpz(i),fmass(i),
     @     ityp(i),iso3(i),charge(i),
     @     dectag+lstcoll(i),ncoll(i),mod(origin(i),100)
 
      return

      end



       
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine spectrans(tstep)
c
c  (when cto 28 is set to 2 this subroutine is called
c  to propagate the spectators along straight lines)
c  
cinput tstep : timestep
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 dtime,energ,tstep
      integer j
      include 'coms.f'

      dtime=tstep

      do 1 j=1,nspec
         energ = p0s(j)
         r0s(j) = r0s(j) + dtime
         rxs(j) = rxs(j) + pxs(j)/energ*dtime
         rys(j) = rys(j) + pys(j)/energ*dtime
         rzs(j) = rzs(j) + pzs(j)/energ*dtime
1     continue

      return
      end
       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function ptsigtot()
c
c     Revision : 1.0
c
c     This function caculates the total cross section of the reaction.
c     (Projectile - target total cross section)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

      include 'comres.f'
      include 'coms.f'
      include 'options.f'

      integer indmn,indmx,itypmn,iso3mn,itypmx,iso3mx
      integer isigline,iline,collclass
      real*8 stot,sigel
      real*8 sigtot

c determine total cross section for reaction:
      if(abs(Ap)+abs(At).gt.2) then
         stot=10.d0*pi*(bdist**2-bmin**2)
      elseif(abs(Ap)+abs(At).eq.2) then
         stot=sigtot(1,2,ecm)
cccccccc for CTOption(7)=1 no elastic cross section:
         if(CTOption(7).eq.1) then
c first sort the two itypes for call to collclass and anndec
            if(abs(ityp(1)).lt.abs(ityp(2))) then
               indmn=1
               indmx=2
            else
               indmn=2
               indmx=1
            endif

            itypmn=ityp(indmn)
            iso3mn=iso3(indmn)
            itypmx=ityp(indmx)
            iso3mx=iso3(indmx)
            isigline=collclass(itypmx,iso3mx,itypmn,iso3mn)
c     the elastic cross section is always the first entry (#3)
            iline=SigmaLn(3,1,isigline)
c!!!!DANGER: does not work for unstable particles (-> detailed balance)
            call crossx(iline,ecm,ityp(1),iso3(1),
     &              fmass(1),ityp(2),iso3(2),fmass(2),sigel)
c
            if(stot-sigel.gt.0) then
               stot=stot-sigel
            else
               stot=sigel
            endif
         endif
      else
         stot=0.d0
      endif
c
      ptsigtot=stot
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine urqmdlogo
c
c Displays the UrQMD Logo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'options.f'
c     we need coms.f for versiontxt
      include 'coms.f'

      integer firsttime
      save firsttime
      if (firsttime.eq.1)return
      firsttime=1

      write (*,*)
     $   "#############################################################"
      write (*,*) 
     $   "##                                                         ##"
      write (*,*)
     $   "## UrQMD ",versiontxt,
     $                    "   University of Frankfurt                ##"
      write (*,*)
     $   "##                  http://urqmd.org                       ##"
      write (*,*)             
     $   "##                  bleicher@th.physik.uni-frankfurt.de    ##"
      write (*,*)
     $   "#############################################################"
      write (*,*)
     $   "##                                                         ##"
      write (*,*)
     $   "##     Please cite when using this model:                  ##"
      write (*,*)
     $   "##     S.A.Bass et al., Prog.Part.Nucl.Phys. 41 (1998) 225 ##"
      write (*,*)
     $   "##     M.Bleicher et al., J.Phys. G25  (1999) 1859         ##"
      write (*,*)
     $   "##                                                         ##"
      write (*,*)
     $   "#############################################################"
      write (*,*)
     $   "##     UrQMD uses Pythia6.409 by T. Sjorstrand             ##"
      write (*,*)
     $   "#############################################################"
      write (*,*)
     $   "##                                                         ##"
      write (*,*)
     $   "##     If hydrodynamic evolution is switched on (CTO 45 1) ##"
      write (*,*)
     $   "##     UrQMD uses the SHASTA algorithm by D. Rischke       ##"
      write (*,*)
     $   "##     Please cite when using the hybrid code:             ##"
      write (*,*) 
     $   "##     D. Rischke et al., Nucl.Phys. A 595 (1995) 346      ##"
      write (*,*)
     $   "##     D. Rischke et al., Nucl.Phys. A 595 (1995) 383      ##"
      write (*,*)
     $   "##     H. Petersen et al., Phys.Rev. C78 (2008) 044901     ##"
      write (*,*)        
     $   "##                                                         ##"
      write (*,*)        
     $   "#############################################################"


      return
      end
