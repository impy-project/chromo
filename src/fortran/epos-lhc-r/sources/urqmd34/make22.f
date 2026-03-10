c$Id: make22.f,v 1.42 2007/01/30 14:50:25 bleicher Exp $
C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine make22(iio,e,ii1,iiz1,mm1,xfac1,ii2,iiz2,mm2,xfac2)
c
cinput iio   : label for exit-channel
cinput e     : $\sqrt{s}$ of process
cinput ii1   : ID of incoming particle 1
cinput iiz1  : $2\cdot I_3$ of incoming particle 1
cinput mm1   : mass of incoming particle 1
cinput xfac1 : scaling factor for preformed hadron
cinput ii2   : ID of incoming particle 2
cinput iiz2  : $2\cdot I_3$ of incoming particle 2
cinput mm2   : mass of incoming particle 2
cinput xfac2 : scaling factor for preformed hadron
c
c  output:   exit channel via common-blocks in {\tt newpart.f}
c
c {\tt make22} generates the final state for all scatterings and
c decays. Due to the diverse nature of the interactions handled
c many special cases have to be taken care of. The label {\tt iio}
c matches in most cases the respective label in subroutine {\tt crossx},
c which returns the respective partial cross sections.
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none

      include 'coms.f'
      include 'options.f'
      include 'newpart.f'
      include 'comres.f'

      integer io,iio,i1,i2,i3,i4,i,k,i1p,i1m,i2p,i2m
      integer ifqrk1,ifqrk2,ifdiq1,ifdiq2,ifqrk3,ifqrk4,ifdiq3,ifdiq4
      integer iq1(2),iq2(3),kq1,kq2,kqq1,kqq2,iii,iff(3)
      real*8 sig,e,m1,m2,m3,m4,gam,mm1,mm2,m1m,m2m,xfac1,xfac2

      integer iz1,iz2,iz3,iz4,errflg,icnt,ntry,ib1,ib2,i1old,i2old
      integer jpsi
      logical bit

c...functions
      integer isoit,whichres,strit
      real*8 getmass,fmsr,ranff,pcms,massit,widit,mminit

c...vacuum quantumnumber(s) for special string decay (don't touch)
      real*8 valint(1)
      common /values/ valint

c...string
      real*8 b1,b2,ba1(3),ba2(3)
      integer j,l,ii1,ii2,iiz1,iiz2,iexopt,iddum,ibar,jbar
      logical fboost,switips

c... modified string decay/high mass resonance
      integer idcode,idpars

c.. needed for pythia
      integer pytlist(7),itag0,itag1,itag2,pythflag
      real*8 minsrt,p,ppythia

      common /pythflag/ pythflag

      integer ident(2,mprt)
      real*8 part(9,mprt),ms1,ms2,msmin1,msmin2,tau,esum

c.. pythia particle list
      data pytlist/1,-1,27,40,49,55,101/

c.. Needed for strangeness exchange process
      integer channel
      real*8 dummy
      real*8 rdummy
      real*8 sig1,sig2,sig3,sig4
      real*8 STREXCHANGE
      integer si1,si2,siz1,siz2

      real*8 DIVpcms


      bit=.true.
      switips=.false.
      io=mod(iio,200)
      i1=ii1
      i2=ii2
      iz1=iiz1
      iz2=iiz2
      m1=mm1
      m2=mm2
      icnt=0
      ntry=0
      ibar=0
chp flag for Pythia call because of angular ditribution
      pythflag=0

      if(iabs(i1).lt.minmes)ibar=ibar+isign(1,i1)
      if(iabs(i2).lt.minmes)ibar=ibar+isign(1,i2)

c in case of a MB-reaction, sort particles (but keep track of
c any id-switch with the 'switips'-flag)
      if(iabs(i1).ge.minmes.and.iabs(i1).le.maxmes.and.
     &   iabs(i2).ge.minbar.and.iabs(i2).le.maxbar)then
        call swpizm(i1,iz1,m1,i2,iz2,m2)
        switips=.true.
      endif

      if(i1+i2.eq.0.and.iz1+iz2.eq.0.and.CTOption(20).ne.0.and.
     .    io.gt.20)goto 27 !e+e-

      if(io.lt.0)goto(100,100,100,100,100,100,100,29)-io

c      if(i1+i2.gt.2)write(6,*)'make22:',i1,i2
 1007 continue
      goto(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,9,17,9,17,20,
     , 9,13,23,15,12,26,27,15,14,100,29,100,14,15,14,36,36,13,14,
     , 13,41,42,43,44,45,46,47,48,49,50,51,52,53,54,
     , 13,56,57,58,59,60,61,62,63,64)io

c      write(6,*)'make22: unknown channel requested io:',io
c      write(6,*)'  ',e,i1,iz1,m1,i2,iz2,m2
      stop 137

 1    continue
c...pp->ND
      i3=minnuc
      i4=mindel
      m3=massit(i3)

      call getmas(massit(i4),widit(i4),i4,isoit(i4),mminit(i4),
     .     e-m3,m3,m4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)
      goto 1008

 2    continue
c...pp->pp*
      i3=minnuc
      m3=massit(i3)
      if(bit)call getres(io,e,minnuc+1,maxnuc,i4)
      call getmas(massit(i4),widit(i4),i4,isoit(i4),mminit(i4),
     .     e-m3,m3,m4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)
      goto 1008

 3    continue
c...pp->ND*
      i3=minnuc
      m3=massit(i3)
      if(bit)call getres(io,e,mindel+1,maxdel,i4)
      call getmas(massit(i4),widit(i4),i4,isoit(i4),mminit(i4),
     .     e-m3,m3,m4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)
      goto 1008

 4    continue
c...pp->DD
      i3=mindel
      i4=mindel
      call getmas(massit(i3),widit(i3),i3,isoit(i3),
     .     mminit(i4),e-mminit(i4),mminit(i4),m3)
      call getmas(massit(i4),widit(i4),i4,isoit(i4),mminit(i4),
     .     e-m3,m3,m4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)
      goto 1008

 5    continue
c...pp->DN*, DN* with factor 4/3
      i3=mindel
      if(bit)call getres(io,e,minnuc+1,maxnuc,i4)
      call getmas(massit(i3),widit(i3),i3,isoit(i3),
     .     mminit(i3),e-mminit(i4),mminit(i4),m3)
      call getmas(massit(i4),widit(i4),i4,isoit(i4),mminit(i4),
     .     e-m3,m3,m4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)
      goto 1008

 6    continue
c...pp->DD*, DN* with factor 4/3
      i3=mindel
      if(bit)call getres(io,e,mindel+1,maxdel,i4)
      call getmas(massit(i3),widit(i3),i3,isoit(i3),
     .     mminit(i3),e-mminit(i4),mminit(i4),m3)
      call getmas(massit(i4),widit(i4),i4,isoit(i4),mminit(i4),
     .     e-m3,m3,m4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)
      goto 1008

 7    continue
c...pp -> generate N*N*,N*D*,D*D*
      m3=fmsr(mminit(mindel),e-mresmin)
      m4=fmsr(mminit(mindel),e-m3)
      if(bit)i3=whichres(m3,3)
      if(bit)i4=whichres(m4,3)
      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)
      goto 1008

 8    continue
c...ND->DD
      i3=mindel
      i4=mindel
      call getmas(massit(i3),widit(i3),i3,isoit(i3),
     .     mminit(i4),e-mminit(i4),mminit(i4),m3)
      call getmas(massit(i4),widit(i4),i4,isoit(i4),mminit(i4),
     .     e-m3,m3,m4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)
      goto 1008

 9    continue
      write(6,*)'make22: channel no.',io,'not implemented.'
      stop 137

 10   continue
c...MB->B',MM->M* annihilations
      call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sig,gam)

      goto 2002
c      return

 11   continue
c...MM->M'
      call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sig,gam)

      goto 2002
c      return


 12   continue
      write(6,*)'make22: channel no.',io,'should correspond to',
     ,'total cross section, i.e. make22 should not be called.'
      stop 137

 13   continue
c...elastic scattering
c
      if(switips)then
         call swpizm(i1,iz1,m1,i2,iz2,m2)
         switips=.false.
      endif

      call setizm(i1,iz1,m1,i2,iz2,m2,i3,iz3,m3,i4,iz4,m4)

      if(mminit(i4)+mminit(i3).gt.e)then
c         write(6,*)'make22(el):threshold violated'
c         write(6,*)'m3:',m3,mminit(i3)
c         write(6,*)'m4:',m4,mminit(i4)
      else
c
c
        if(m3.lt.mminit(i3))m3=mminit(i3)
        if(m4.lt.mminit(i4))m4=mminit(i4)
      end if

      goto 2001

 14   continue
c...inelastic scattering (aqm for nonstrange resonances)
c   for arbitrary resonances getinw should be modified
c   arbitrary particles are excited (whithout charge exchange)

c...get minimal masses
c...allow also for delta excitation in nucleon-nucleon reactions
      call getirg(i1,i1m,i1p)
      if (i1p.eq.maxnuc) i1p=maxdel
      m1m=mminit(i1m)
      call getirg(i2,i2m,i2p)
      if (i2p.eq.maxnuc) i2p=maxdel
      m2m=mminit(i2m)

      if(e-m1m-m2m.le.0d0)then
c...elastic scattering
        call setizm(i1,iz1,m1,i2,iz2,m2,i3,iz3,m3,i4,iz4,m4)
        goto 2001
c        return
      end if

c...masses
      if(ranff(0).lt.5d-1)then
         m3=fmsr(m1m,min(mresmax,e-m2m))
         m4=fmsr(m2m,min(mresmax,e-m3))
      else
         m4=fmsr(m2m,min(mresmax,e-m1m))
         m3=fmsr(m1m,min(mresmax,e-m4))
      end if

      if(e-m3-m4.lt.0d0)goto 13

c...itype3
      if(i1m.lt.i1p.and.widit(i1m).lt.1d-3.and.
     .    m3.le.massit(i1m+1)-0.5*widit(i1m+1))then
c...lowest itype of this kind of particles is stable
c & mass .lt. approximate minimal mass of lowest itype + 1
        m3=massit(i1m)
        i3=i1m
      else if(i1m.eq.i1p.and.widit(i1m).lt.1d-3) then
c...class with only one narrow particle
        i3=i1
        m3=massit(i3)
      else
c...itypes of this kind are all unstable
        call whichi(i3,i1m,i1p,m3)
      end if

c...itype4
      if(i2m.lt.i2p.and.widit(i2m).lt.1d-3.and.
     .    m4.le.massit(i2m+1)-0.5*widit(i2m+1))then
c...lowest itype of this kind of particles is stable
        m4=massit(i2m)
        i4=i2m
      else if(i2m.eq.i2p.and.widit(i2m).lt.1d-3) then
c...class with only one narrow particle
        i4=i2
        m4=massit(i4)
      else
c...itypes of this kind are all unstable
        call whichi(i4,i2m,i2p,m4)
      end if

c...no charge transfer
      iz3=iz1
      iz4=iz2
      i3=isign(i3,ii1)
      i4=isign(i4,ii2)
      goto 2001
c      return

 15    continue
c...XX -> 2 strings
      if(CTOption(12).ne.0)then
        write(6,*)' *** error(make22): string section is called ',
     .            'while strings are switched off:CTOption(12).ne.0'
        stop 137
      end if

c
      if(switips)then
         call swpizm(i1,iz1,m1,i2,iz2,m2)
         switips=.false.
      endif


c store old itypes
      i1old=i1
      i2old=i2

 155  continue

c allow for deexcitation: 'excitation' starts from groundstate

      if(iabs(i1).ge.minmes)then
c.. meson resonances may also be deexcitated: assign the particle
c id's of the lowest multiplet (with the same quark content)
        call ityp2id(i1,iz1,iq1(1),iq1(2))
        if(iabs(iq1(1)).gt.iabs(iq1(2)))then
          iddum=iq1(1)
          iq1(1)=iq1(2)
          iq1(2)=iddum
        endif
        iddum=isign(100*iabs(iq1(1))+10*iabs(iq1(2)),iq1(1))
        call id2ityp(iddum,0d0,i1,iz1)
      else
         call getirg(i1,i1m,i1p)
      endif

c same for second particle
      if(iabs(i2).ge.minmes)then
        call ityp2id(i2,iz2,iq1(1),iq1(2))
        if(iabs(iq1(1)).gt.iabs(iq1(2)))then
          iddum=iq1(1)
          iq1(1)=iq1(2)
          iq1(2)=iddum
        endif
        iddum=isign(100*iabs(iq1(1))+10*iabs(iq1(2)),iq1(1))
        call id2ityp(iddum,0d0,i2,iz2)
      else
         call getirg(i2,i2m,i2p)
      endif

c
c BEWARE: now i1 and i2 do not contain anymore the ityps of the ingoing
c particles, but the ityps of the lowest possible states (groundstates).
c This is needed in order for the string-excitation to be able to
c excite ALL states and not only those above the state of the ingoing
c particle.
c If you need the old ityps (i.e. for elastic scattering) then reset them
c via i1old and i2old

c start if/then structure for pythia call here
      minsrt=10.
      itag0=0
      itag1=0
      itag2=0

c.. call pythia only for hard collision (Q>1.5 GeV)
c.. ppythia parametrizes the appropriate probability
      p=ranff(0)
      ppythia=1d0-1.03293d0*exp((-e)/246.771d0)
      if (p.lt.ppythia) itag0=1

c check whether particle type can be handled by PYTHIA
      do 123 i=1,7
         if(i1.eq.pytlist(i)) itag1=1
         if(i2.eq.pytlist(i)) itag2=1
 123  continue

c pythia for high energies

      if(e.gt.minsrt.and.itag0.eq.1.and.itag1.eq.1.and.itag2.eq.1
c... PYTHIA switch
     &  .and.CTOption(44).eq.1) then

         call upyth(i1,iz1,i2,iz2,e)
chp flag is necessary to fix angular distribution
         pythflag=1

      else

c no pythia call, use ordinary string routine
c need to set nstring1 to particle number out of pythia
c set nstring2=0
c initialize inew, pnew, xnew ... arrays from pythia common blocks

c--- end pythia block -----

c if the incoming masses are changed due to excitation, then the
c new masses should not be lower than:
      m1m=mminit(i1)
      m2m=mminit(i2)

c..discriminate between BB and MB collisions:
      ib1=1
      ib2=1
      if(iabs(i1).ge.minmes)ib1=0
      if(iabs(i2).ge.minmes)ib2=0


c set minimum energy for the two strings
      msmin1=m1m+CTParam(2)
      msmin2=m2m+CTParam(2)

c no (meson) string below 1 gev:
      if(msmin1.lt.1d0)msmin1=1d0
      if(msmin2.lt.1d0)msmin2=1d0

c if energy too low: do elastic collision
      if(m1m+m2m+CTParam(34).ge.e)then
c        write(*,*)'make22: not enough energy for string exc. ->elastic/
c     &deexcitation'
c        write(*,*)' i1, i2, m1, m2, e: ',i1,i2,m1,m2,e
         i1=i1old
         i2=i2old
        goto 13
      endif


c convert to quark-IDs
         call ityp2id(i1,iz1,ifdiq1,ifqrk1)
         call ityp2id(i2,iz2,ifdiq2,ifqrk2)

         iexopt=CTOption(22)
 81      continue
c 100 tries for excitation, otherwise elastic scattering
         ntry=ntry+1
         if(ntry.gt.1000)then
c          write(*,*)'make22: too many tries for string exc. ->elastic/
c     &           deexcitation'
c          write(*,*)' i1, i2, m1, m2, e: ',i1,i2,m1,m2,e
         i1=i1old
         i2=i2old
          goto 13
         endif

c string-excitation:
c get string masses ms1,ms2 and the leading quarks
         call STREXCT(IFdiq1,IFqrk1,ib1,M1m,
     &           ifdiq2,ifqrk2,ib2,M2m,E,
     &           iexopt,
     &           ba1,ms1,ba2,ms2,
     &           ifdiq3,ifqrk3,ifdiq4,ifqrk4)
c the boost parameters are now fixed for the masses ms1, ms2. If the
c masses will be changed, set the parameter fboost to "false":
      fboost=.true.

c accept deexcitation of one of the hadrons:
      if(ms1.le.msmin1.and.ms2.ge.msmin2)then
        ms1=massit(i1)
        fboost=.false.
      else if(ms2.le.msmin2.and.ms1.ge.msmin1)then
        ms2=massit(i2)
        fboost=.false.
c don't accept elastic-like (both masses too low):
      else if(ms1.lt.msmin1.and.ms2.lt.msmin2)then
        goto 81
      endif

c in case of deexcitation new masses are necessary
c single diffractive, mass excitation according 1/m
      if(ms1.le.msmin1)then
        ms2=fmsr(msmin2,e-ms1)
        fboost=.false.
      elseif(ms2.le.msmin2)then
        ms1=fmsr(msmin1,e-ms2)
        fboost=.false.
      endif

c quark-quark scattering -> only elastic !
      if(xfac1.lt..999d0.and.xfac2.lt..999d0
     &     .and.ranff(0).lt..25) then
         i1=i1old
         i2=i2old
         goto 13
      endif

c single diffractive if one particle is a quark state,
c mass excitation according 1/m
      if(xfac1.lt..999d0.and.ranff(0).lt.0.5)then
        ms1=massit(i1)
        ms2=fmsr(msmin2,e-ms1)
        fboost=.false.
      elseif(xfac2.lt..999d0.and.ranff(0).lt.0.5)then
        ms2=massit(i2)
        ms1=fmsr(msmin1,e-ms2)
        fboost=.false.
      endif

c take care that the particles will be able to decay lateron:
      if(ms1.lt.m1m)then
        ms1=m1m
        fboost=.false.
      endif
      if(ms2.lt.m2m)then
        ms2=m2m
        fboost=.false.
      endif

c avoid energy conservation violation
      if(ms1+ms2.gt.e)goto 155

       if(CTOption(22).ne.1.or..not.fboost)then
c the boost parameters have to be calculated:
         b1=2.*e*pcms(e,ms1,ms2)/(e**2+ms1**2-ms2**2)
         b2=2.*e*pcms(e,ms1,ms2)/(e**2-ms1**2+ms2**2)
         ba1(3)=+b1
         ba2(3)=-b2
         do 151 j=1,2
            ba1(j)=0d0
 151        ba2(j)=0d0
      end if

c now write information to newpart common-blocks
      mstring(1)=ms1
      mstring(2)=ms2

c string#1
      if(ms1.le.msmin1)then
        nstring1=1
        l=1
        do j=1,3
          part(j,l)=0d0
          part(j+5,l)=0d0
        enddo
        part(4,l)=ms1
        part(4+5,l)=0d0
        part(5,l)=ms1
        ident(1,l)=i1
        ident(2,l)=iz1
      else
       call qstring(ifdiq3,ifqrk3,ms1,part,ident,nstring1)
       if(nstring1.eq.0) goto 155
      end if

      esum=0d0

      jbar=0

      do l=1,nstring1
         pnew(5,l)=part(5,l)
         itypnew(l)=ident(1,l)

         if(iabs(ident(1,l)).lt.minmes)jbar=jbar+isign(1,ident(1,l))

         i3new(l)=ident(2,l)
         do j=1,4
           pnew(j,l)=part(j,l)
           xnew(j,l)=part(j+5,l)
         enddo
         pnew(3,l)=part(3,l)
         xnew(3,l)=part(3+5,l)
         call rotbos(0d0,0d0,ba1(1),ba1(2),ba1(3),
     ,    pnew(1,l),pnew(2,l),pnew(3,l),pnew(4,l))
         call rotbos(0d0,0d0,ba1(1),ba1(2),ba1(3),
     ,    xnew(1,l),xnew(2,l),xnew(3,l),xnew(4,l))
         esum=esum+pnew(4,l)
      enddo
      call leadhad(1,nstring1,1)


c string #2
      if(ms2.le.msmin2) then
        nstring2=1
        l=1
        do j=1,3
          part(j,l)=0d0
          part(j+5,l)=0d0
        enddo
        part(4,l)=ms2
        part(4+5,l)=0d0
        part(5,l)=ms2
        ident(1,l)=i2
        ident(2,l)=iz2
      else
       call qstring(ifdiq4,ifqrk4,ms2,part,ident,nstring2)
       if(nstring2.eq.0) goto 155
      end if

      esum=0d0
       do l=1,nstring2
        pnew(5,nstring1+l)=part(5,l)
        itypnew(nstring1+l)=ident(1,l)

        if(iabs(ident(1,l)).lt.minmes)jbar=jbar+isign(1,ident(1,l))

        i3new(nstring1+l)=ident(2,l)
        do j=1,4
          pnew(j,nstring1+l)=part(j,l)
          xnew(j,nstring1+l)=part(j+5,l)
        enddo
        pnew(3,nstring1+l)=-pnew(3,nstring1+l)
        xnew(3,nstring1+l)=-xnew(3,nstring1+l)
        call rotbos(0d0,0d0,ba2(1),ba2(2),ba2(3),
     ,    pnew(1,nstring1+l),pnew(2,nstring1+l),pnew(3,nstring1+l),
     ,     pnew(4,nstring1+l))
        call rotbos(0d0,0d0,ba2(1),ba2(2),ba2(3),
     ,    xnew(1,nstring1+l),xnew(2,nstring1+l),xnew(3,nstring1+l),
     ,     xnew(4,nstring1+l))
        esum=esum+pnew(4,nstring1+l)
      enddo
      call leadhad(nstring1+1,nstring1+nstring2,1)

c error check
      if(ibar.ne.jbar)then
         write(6,*)' *** (E) no baryon number conservation', ibar,jbar
         write(6,*)'     ',i1,i2,ms1,ms2
         write(6,'(5i4)')(itypnew(l),l=1,nstring1+nstring2)
         end if

cmb,sab
c end of bracket for urqmd vs. pythia string routine
      endif

      return


 718  format(i2,i4,i3,1x,10(f10.4,1x))



 17   continue

      iz3=iz1
      iz4=iz2
      i3=i1
      i4=i2


      m3=m1
      m4=m2

c the following lines MUST be there in order to set nucleons on-shell
c after their first collision
      if(m3.lt.mminit(i3))m3=mminit(i3)
      if(m4.lt.mminit(i4))m4=mminit(i4)


      goto 2001

 20   continue
c...decays
      if(ityptd(1,pslot(1)).eq.0) then
c normal decay
c note: m4,i4 and iz4 are dummies in this call
         call anndec(1,m1,i1,iz1,m4,i4,iz4,e,sig,gam)
      else
c forward time-delay
         pnew(5,1)=fmasstd(1,pslot(1))
         itypnew(1)=ityptd(1,pslot(1))
         i3new(1)=iso3td(1,pslot(1))
         pnew(5,2)=fmasstd(2,pslot(1))
         itypnew(2)=ityptd(2,pslot(1))
         i3new(2)=iso3td(2,pslot(1))
      endif

      if(nexit.eq.2) then
         i3=itypnew(1)
         iz3=i3new(1)
         m3=pnew(5,1)
         i4=itypnew(2)
         iz4=i3new(2)
         m4=pnew(5,2)
         goto 2001
      else
c three or four body decay
         nstring1=1
         nstring2=nexit-1
         do i=1,4
            do j=1,nexit
               pnew(i,j)=0d0
               xnew(i,j)=0d0
            end do
         end do
         mstring(1)=pnew(5,1)
c
         mstring(2)=pnew(5,2)
         do 91 j=3,nexit
            mstring(2)=mstring(2)+pnew(5,j)
 91      continue
c
c now call routine for momentum phase space...
         call nbodydec(e)
         return
      endif

 23   continue
c...annihilation -> string
      if(CTOption(12).ne.0)then
        write(6,*)' *** error(make22): string section is called ',
     .            'while strings are switched off:CTOption(12).ne.0'
        stop 137
      end if

      ms1=e/2.
      ms2=e/2.
      mstring(1)=ms1
      mstring(2)=ms2

c determine flavour content of b-bbar-system
      call ityp2id(i1,iz1,ifdiq1,ifqrk1)
      call ityp2id(i2,iz2,ifdiq2,ifqrk2)
c...create string 1 out of quark-antiquark pair
      call qstring(ifqrk1,ifqrk2,ms1,part,ident,nstring1)
      esum=0d0

      do k=1,nstring1
        l=k
        do j=1,4
          pnew(j,l)=part(j,l)
          xnew(j,l)=part(j+5,l)
        enddo
        esum=esum+pnew(4,l)
        pnew(5,l)=part(5,l)
        itypnew(l)=ident(1,l)
        i3new(l)=ident(2,l)
        tau=part(9,l)/ (part(4,l)/part(5,l))
      enddo
      call leadhad(1,nstring1,0)

c...create string 2 out of diquark-antidiquark-pair
c   use one quark and one antiquark for the string-ends:
       ifqrk1=int(ifdiq1/1000)
       ifqrk2=int(ifdiq2/1000)
c...store remaining flavour quantum numbers in ctp(26), they will
c   be passed to the 'clustr'-routine:
       ifdiq1=mod(ifdiq1/100,10)
       ifdiq2=mod(ifdiq2/100,10)
       valint(1)=dble(((abs(ifdiq1*10.d0)+abs(ifdiq2*1.d0))
     &            *isign(1,ifdiq1)))
       valint(1)=sign(valint(1),dble(ifdiq1))
       call qstring(ifqrk1,ifqrk2,ms2,part,ident,nstring2)
       valint(1)=0.d0
      esum=0d0

      do l=1,nstring2
        do j=1,4
          pnew(j,nstring1+l)=part(j,l)
          xnew(j,nstring1+l)=part(j+5,l)
        enddo
        esum=esum+pnew(4,nstring1+l)
        pnew(5,nstring1+l)=part(5,l)
        itypnew(nstring1+l)=ident(1,l)
        i3new(nstring1+l)=ident(2,l)
        tau=part(9,nstring1+l)/ (part(4,l)/part(5,l))
      enddo
      call leadhad(nstring1+1,nstring1+nstring2,0)

      return

 26   continue
c...elastic MB scattering (the outgoing particle id's must not be
c   assigned randomly like at label 2001)
c
      if(switips)then
         call swpizm(i1,iz1,m1,i2,iz2,m2)
         switips=.false.
      endif

      call setizm(i1,iz1,m1,i2,iz2,m2,i3,iz3,m3,i4,iz4,m4)

      if(mminit(i4)+mminit(i3).gt.e)then
         write(6,*)'make22(el):threshold violated'
         write(6,*)'m3:',m3,mminit(i3)
         write(6,*)'m4:',m4,mminit(i4)
      else
        if(m3.lt.mminit(i3))m3=mminit(i3)
        if(m4.lt.mminit(i4))m4=mminit(i4)
      end if

c... get momenta & fill newpart, 2 particle exit-channel

      nstring1=1
      nstring2=1
      nexit=2
      do i=1,4
         do j=1,2
            pnew(i,j)=0d0
            xnew(i,j)=0d0
         end do
      end do

c...boost to 2-particle cms

      pnew(3,1)=pcms(e,m3,m4)
      pnew(3,2)=-pcms(e,m3,m4)

      pnew(4,1)=sqrt(m3**2+pnew(3,1)**2)
      pnew(4,2)=sqrt(m4**2+pnew(3,2)**2)

      pnew(5,1)=m3
      mstring(1)=m3
      itypnew(1)=i3
      i3new(1)=iz3

      pnew(5,2)=m4
      mstring(2)=m4
      itypnew(2)=i4
      i3new(2)=iz4

      return

 27   continue
c XX-> 1 string : e+e- , MB
      if(CTOption(12).ne.0)then
        write(6,*)' *** error(make22): string section is called ',
     .            'while strings are switched off:CTOption(12).ne.0'
        stop 137
      end if

      ms1=e
      mstring(1)=ms1
      mstring(2)=0d0

c determine flavour content of string
      if(CTOption(20).eq.1)then
c..e+e- annihilation
      if(ranff(0).lt.CTParam(6))then
        ifqrk1=3  ! ssbar
        ifqrk2=-3
      else
        call ityp2id(104,0,ifqrk1,ifqrk2)  ! qqbar
      end if
      else
c...MB annihilation. the quark content must be known:
      call ityp2id(i2,iz2,iq1(1),iq1(2))
      call ityp2id(i1,iz1,ifdiq2,iq2(3))

      if(abs(i1).ge.minmes) then
         iq2(1)=ifdiq2
         iq2(2)=iq2(3)
         iq2(3)=0
      else
         iq2(1)=mod(ifdiq2/100,10)
         iq2(2)=int(ifdiq2/1000)
      endif

      do 312 kq1=1,2
       do 412 kq2=1,3
c..two of the quarks must be able to annihilate:
        if(iq1(kq1)+iq2(kq2).eq.0) then
         kqq1=kq1
         kqq2=kq2
         goto 414
        endif
 412   continue
 312  continue
      goto 26 ! could not create double charged strange baryon string
 414  continue
c.. the 'iff'-quarks constitute the produced (anti-)baron
      iff(1)=iq1(3-kqq1)
      iii=1
      do 512 kq2=1,3
       if (kq2.ne.kqq2) then
        iii=iii+1
        iff(iii)=iq2(kq2)
       endif
 512  continue

      if(abs(i1).lt.minmes) then
         call mquarks(iff,ifqrk1,ifqrk2)
      else
         ifqrk1=iff(1)
         ifqrk2=iff(2)
      endif

      endif

c define energy cut
      if ((ms1.gt.CTParam(60))
     $     .or.(abs(iff(1)+iff(2)+iff(3)).eq.9)) then
c...create string 1 out of quark-antiquark pair
       call qstring(ifqrk1,ifqrk2,ms1,part,ident,nstring1)
      else
        if(ms1.lt.1.67)then
c         normal resonance
          call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sig,gam)

c         some exceptions
          if(pnew(5,1).eq.0.or.
     $        strit(ii1)+strit(ii2).ne.strit(itypnew(1))) then
c          get the quark code first
           idcode=idpars(ifqrk1,ifqrk2,.true.,1)
c          search the meson or baryon lists to find the correct one
           call id2itypnew(idcode,ms1,i1,iz1)
          else
           i1=itypnew(1)
           iz1=i3new(1)
          end if

        else
c          continuous spectra
c          get the quark code first
           idcode=idpars(ifqrk1,ifqrk2,.true.,1)
c          search the meson or baryon lists to find the correct one
           call id2itypnew(idcode,ms1,i1,iz1)

        endif

        nstring1=1
        l=1
        do j=1,3
          part(j,l)=0d0
          part(j+5,l)=0d0
        enddo
        part(4,l)=ms1
        part(4+5,l)=0d0
        part(5,l)=ms1
        ident(1,l)=i1
        ident(2,l)=iz1
c        write(*,*)'mb: ',ms1
      end if


      esum=0d0

c primitive bug-fix:
      if(nstring1.eq.0)then
c        write(6,*)'make22: iline 27 not completed. ->elastic'
        goto 26
      endif

      do 101 k=1,nstring1
        l=k
c no leading hadron in e+e-
        if(CTOption(20).eq.1)then
          leadfac(l)=0.0d0
        endif
        do 102 j=1,4
          pnew(j,l)=part(j,l)
          xnew(j,l)=part(j+5,l)
 102    continue
        esum=esum+pnew(4,l)
        pnew(5,l)=part(5,l)
        itypnew(l)=ident(1,l)
        i3new(l)=ident(2,l)
        tau=part(9,l)/ (part(4,l)/part(5,l))
 101    continue
      if(CTOption(20).eq.1)then
        call leadhad(1,nstring1,3)
      else
        call leadhad(1,nstring1,1)
      endif

      nstring2=0
      nexit=nstring1


      return


 29   continue
c...DD->ND detailed balance
      i3=minnuc
      i4=mindel
      m3=massit(i3)
      m4=getmass(e-m3,0)
      if(iabs(iz1+iz2).gt.isoit(i3)+isoit(i4))then
       iz3=-9
       iz4=-9
      else
         nexit=2
         itot(1)=isoit(i3)
         itot(2)=isoit(i4)
         call isocgk4(isoit(i1),iz1,isoit(i2),iz2,itot,i3new,errflg)
         i3=isign(i3,ii1)
         i4=isign(i4,ii2)
         iz3=i3new(1)
         iz4=i3new(2)
      end if

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)
      goto 2001


 100  continue
c...??->NN detailed balance inverse channels
c iso3 are assigned in detbal
      nexit=2
      i3=minnuc
      i4=minnuc
      m3=massit(minnuc)
      m4=massit(minnuc)
      itot(1)=isoit(i3)
      itot(2)=isoit(i4)

      call isocgk4(isoit(i1),iz1,isoit(i2),iz2,itot,i3new,errflg)
      i3=isign(i3,ii1)
      i4=isign(i4,ii2)
      iz3=i3new(1)
      iz4=i3new(2)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)
      goto 2001


 36   continue
c...MB->B',MM->M* annihilations (forward time delay)
      call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sig,gam)
      goto 2003

 41   continue
c pi D -> rho D*

      if(i1.eq.101)then
       i3=104*isign(1,i1)
       iz3=iz1
       m3=massit(i3)
       i4=134*isign(1,i2)
       iz4=iz2
       m4=massit(i4)
      else
       i3=104*isign(1,i2)
       iz3=iz2
       m3=massit(i3)
       i4=134*isign(1,i1)
       iz4=iz1
       m4=massit(i4)
      endif

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008

 42   continue
c pi D* -> rho D

      if(i1.eq.101)then
       i3=104*isign(1,i1)
       iz3=iz1
       m3=massit(i3)
       i4=133*isign(1,i2)
       iz4=iz2
       m4=massit(i4)
      else
       i3=104*isign(1,i2)
       iz3=iz2
       m3=massit(i3)
       i4=133*isign(1,i1)
       iz4=iz1
       m4=massit(i4)
      endif

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008

 43   continue
c rho D -> pi D*

      if(i1.eq.104)then
       i3=101*isign(1,i1)
       iz3=iz1
       m3=massit(i3)
       i4=134*isign(1,i2)
       iz4=iz2
       m4=massit(i4)
      else
       i3=101*isign(1,i2)
       iz3=iz2
       m3=massit(i3)
       i4=134*isign(1,i1)
       iz4=iz1
       m4=massit(i4)
      endif

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008

 44   continue
c rho D* -> pi D

      if(i1.eq.104)then
       i3=101*isign(1,i1)
       iz3=iz1
       m3=massit(i3)
       i4=133*isign(1,i2)
       iz4=iz2
       m4=massit(i4)
      else
       i3=101*isign(1,i2)
       iz3=iz2
       m3=massit(i3)
       i4=133*isign(1,i1)
       iz4=iz1
       m4=massit(i4)
      endif

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008

 45   continue

c  -> D Dbar
      i3=133
      i4=-133
      m3=massit(i3)
      m4=massit(i4)

      if(iz1.eq.0.and.iz2.eq.0)then
       iz3=-1
       iz4=1
      elseif(iz1.gt.0.or.iz2.gt.0)then
       iz3=1
       iz4=1
      else
       iz3=-1
       iz4=-1
      endif


      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008

 46   continue
c  -> D* Dbar, D*Bar + D


      if(ranff(0).gt.0.5d0)then
       i3=134
       i4=-133
      else
       i3=133
       i4=-134
      endif
      m3=massit(i3)
      m4=massit(i4)

      if(iz1.eq.0.and.iz2.eq.0)then
       iz3=-1
       iz4=1
      elseif(iz1.gt.0.or.iz2.gt.0)then
       iz3=1
       iz4=1
      else
       iz3=-1
       iz4=-1
      endif


      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008


 47   continue
c pi J/psi -> D* D*bar

      i3=134
      i4=-134
      m3=massit(i3)
      m4=massit(i4)

      if(iz1.eq.0.and.iz2.eq.0)then
       iz3=-1
       iz4=1
      elseif(iz1.gt.0.or.iz2.gt.0)then
       iz3=1
       iz4=1
      else
       iz3=-1
       iz4=-1
      endif


      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008

 48   continue
c  -> Ds Dbar, DsBar + D


      if(iz1.eq.0)then
       if(i2.gt.0.and.iz2.gt.0)then
        i3=138
        iz3=0
        i4=-133
        iz4=1
       elseif(i2.lt.0.and.iz2.lt.0)then
        i3=-138
        iz3=0
        i4=133
        iz4=-1
       elseif(i2.gt.0.and.iz2.lt.0)then
        i3=138
        iz3=0
        i4=-133
        iz4=-1
       else
        i3=-138
        iz3=0
        i4=133
        iz4=1
       endif
      else
       if(i1.gt.0.and.iz1.gt.0)then
        i3=138
        iz3=0
        i4=-133
        iz4=1
       elseif(i1.lt.0.and.iz1.lt.0)then
        i3=-138
        iz3=0
        i4=133
        iz4=-1
       elseif(i1.gt.0.and.iz1.lt.0)then
        i3=138
        iz3=0
        i4=-133
        iz4=-1
       else
        i3=-138
        iz3=0
        i4=133
        iz4=1
       endif
      endif

      m3=massit(i3)
      m4=massit(i4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008

 49   continue
c  -> Ds D*bar, DsBar + D*

      if(iz1.eq.0)then
       if(i2.gt.0.and.iz2.gt.0)then
        i3=138
        iz3=0
        i4=-134
        iz4=1
       elseif(i2.lt.0.and.iz2.lt.0)then
        i3=-138
        iz3=0
        i4=134
        iz4=-1
       elseif(i2.gt.0.and.iz2.lt.0)then
        i3=138
        iz3=0
        i4=-134
        iz4=-1
       else
        i3=-138
        iz3=0
        i4=134
        iz4=1
       endif
      else
       if(i1.gt.0.and.iz1.gt.0)then
        i3=138
        iz3=0
        i4=-134
        iz4=1
       elseif(i1.lt.0.and.iz1.lt.0)then
        i3=-138
        iz3=0
        i4=134
        iz4=-1
       elseif(i1.gt.0.and.iz1.lt.0)then
        i3=138
        iz3=0
        i4=-134
        iz4=-1
       else
        i3=-138
        iz3=0
        i4=134
        iz4=1
       endif
      endif

      m3=massit(i3)
      m4=massit(i4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008

 50   continue
c  -> Ds* Dbar, Ds*Bar + D

      if(iz1.eq.0)then
       if(i2.gt.0.and.iz2.gt.0)then
        i3=139
        iz3=0
        i4=-133
        iz4=1
       elseif(i2.lt.0.and.iz2.lt.0)then
        i3=-139
        iz3=0
        i4=133
        iz4=-1
       elseif(i2.gt.0.and.iz2.lt.0)then
        i3=139
        iz3=0
        i4=-133
        iz4=-1
       else
        i3=-139
        iz3=0
        i4=133
        iz4=1
       endif
      else
       if(i1.gt.0.and.iz1.gt.0)then
        i3=139
        iz3=0
        i4=-133
        iz4=1
       elseif(i1.lt.0.and.iz1.lt.0)then
        i3=-139
        iz3=0
        i4=133
        iz4=-1
       elseif(i1.gt.0.and.iz1.lt.0)then
        i3=139
        iz3=0
        i4=-133
        iz4=-1
       else
        i3=-139
        iz3=0
        i4=133
        iz4=1
       endif
      endif

      m3=massit(i3)
      m4=massit(i4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008


 51   continue
c  -> pi + J/Psi

      if(cross(51).eq.1)then
       jpsi=135
      elseif(cross(51).eq.2)then
       jpsi=137
      elseif(cross(51).eq.3)then
       jpsi=136
      else
       write(*,*)'error in make22'
      endif


      if(iz1.gt.0.and.iz2.gt.0)then
       i3=101
       iz3=2
       i4=jpsi
       iz4=0
      elseif(iz1.lt.0.and.iz2.lt.0)then
       i3=101
       iz3=-2
       i4=jpsi
       iz4=0
      else
       i3=101
       iz3=0
       i4=jpsi
       iz4=0
      endif

      m3=massit(i3)
      m4=massit(i4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008

 52   continue
c  -> rho + J/Psi

      if(cross(52).eq.1)then
       jpsi=135
      elseif(cross(52).eq.2)then
       jpsi=137
      elseif(cross(52).eq.3)then
       jpsi=136
      else
       write(*,*)'error in make22'
      endif
c      write(*,*)'JPsityp',jpsi,'52'
      if(iz1.gt.0.and.iz2.gt.0)then
       i3=104
       iz3=2
       i4=jpsi
       iz4=0
      elseif(iz1.lt.0.and.iz2.lt.0)then
       i3=104
       iz3=-2
       i4=jpsi
       iz4=0
      else
       i3=104
       iz3=0
       i4=jpsi
       iz4=0
      endif

      m3=massit(i3)
      m4=massit(i4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008

 53   continue
c  -> K + J/Psi

      if(cross(53).eq.1)then
       jpsi=135
      elseif(cross(53).eq.2)then
       jpsi=137
      elseif(cross(53).eq.3)then
       jpsi=136
      else
       write(*,*)'error in make22'
      endif
      if(iz1.eq.0)then
       if(i1.gt.0)then
        if(iz2.gt.0)then
         i3=106
         iz3=1
         i4=jpsi
         iz4=0
        else
         i3=106
         iz3=-1
         i4=jpsi
         iz4=0
        endif
       else
        if(iz2.gt.0)then
         i3=-106
         iz3=1
         i4=jpsi
         iz4=0
        else
         i3=-106
         iz3=-1
         i4=jpsi
         iz4=0
        endif
       endif
      else
       if(i2.gt.0)then
        if(iz1.gt.0)then
         i3=106
         iz3=1
         i4=jpsi
         iz4=0
        else
         i3=106
         iz3=-1
         i4=jpsi
         iz4=0
        endif
       else
        if(iz1.gt.0)then
         i3=-106
         iz3=1
         i4=jpsi
         iz4=0
        else
         i3=-106
         iz3=-1
         i4=jpsi
         iz4=0
        endif
       endif
      endif

      m3=massit(i3)
      m4=massit(i4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008

 54   continue
c  -> K* + J/Psi

      if(cross(54).eq.1)then
       jpsi=135
      elseif(cross(54).eq.2)then
       jpsi=137
      elseif(cross(54).eq.3)then
       jpsi=136
      else
       write(*,*)'error in make22'
      endif

      if(iz1.eq.0)then
       if(i1.gt.0)then
        if(iz2.gt.0)then
         i3=108
         iz3=1
         i4=jpsi
         iz4=0
        else
         i3=108
         iz3=-1
         i4=jpsi
         iz4=0
        endif
       else
        if(iz2.gt.0)then
         i3=-108
         iz3=1
         i4=jpsi
         iz4=0
        else
         i3=-108
         iz3=-1
         i4=jpsi
         iz4=0
        endif
       endif
      else
       if(i2.gt.0)then
        if(iz1.gt.0)then
         i3=108
         iz3=1
         i4=jpsi
         iz4=0
        else
         i3=108
         iz3=-1
         i4=jpsi
         iz4=0
        endif
       else
        if(iz1.gt.0)then
         i3=-108
         iz3=1
         i4=jpsi
         iz4=0
        else
         i3=-108
         iz3=-1
         i4=jpsi
         iz4=0
        endif
       endif
      endif

      m3=massit(i3)
      m4=massit(i4)

      if(ranff(0).gt.0.5d0)call swpizm(i3,iz3,m3,i4,iz4,m4)

      goto 1008

 55   continue
c Baryon Charmonium -> Baryon Charmonium elastic

      iz3=iz1
      iz4=iz2
      i3=i1
      i4=i2

      m3=m1
      m4=m2
      itot(1)=isoit(i1)
      itot(2)=isoit(i2)

      goto 2001



 56   continue
c Baryon Charmonium -> Baryon

      if(i1.lt.i2)then
      itypnew(1)=i1
      i3new(1)=iz1
      pnew(5,1)=e
      else
      itypnew(1)=i2
      i3new(1)=iz2
      pnew(5,1)=e
      endif
      goto 2003

  57  continue


c Kaon + Baryon strangeness exchange
      nexit=2
c K- + p -> Y + pi
      if(i1.eq.1.and.iz1.eq.1.and.        !proton
     &   i2.eq.-106.and.iz2.eq.-1) then   !K-
         sig1 = STREXCHANGE(e,i1,iz1,i2,iz2,40,2,101,-2) !Sigma+ + Pi-
         sig2 = STREXCHANGE(e,i1,iz1,i2,iz2,40,-2,101,2) !Sigma- + Pi+
         sig3 = STREXCHANGE(e,i1,iz1,i2,iz2,40,0,101,0)  !Lambda0/Sigma0 + Pi0
         call getbran((/sig1,sig2,sig3/),1,3,dummy,1,3,channel)
         i3=40
         i4=101
         select case(channel)
           case(1) !Sigma+ + Pi-
             iz3=2
             iz4=-2
           case(2) !Sigma- + Pi+
             iz3=-2
             iz4=2
           case(3) !Lambda0/Sigma0 + Pi0
             iz3=0
             iz4=0
             if(ranff(0).gt.CTParam(72)) then !Lambda
               i3=27
             endif
           case default
             write(6,*) 'Error: non-existent channel selected'
             write(6,*) i1,iz1,i2,iz2,channel
             stop 137
         end select
c K- + n -> Y + pi
      else if (i1.eq.1.and.iz1.eq.-1.and.  !neutron
     &   i2.eq.-106.and.iz2.eq.-1) then    !K-
        sig1 = STREXCHANGE(e,i1,iz1,i2,iz2,40,-2,101,0) !Sigma- + Pi0
        sig2 = STREXCHANGE(e,i1,iz1,i2,iz2,40,0,101,-2) !Lambda0/Sigma0 + Pi-
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        i3=40
        i4=101
        select case(channel)
          case(1) !Sigma- + Pi0
            iz3=-2
            iz4=0
          case(2) !Lambda0/Sigma0 + Pi-
            iz3=0
            iz4=-2
            if(ranff(0).gt.CTParam(72)) then !Lambda
              i3=27
            endif
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c K+ + ap -> aY + pi
      else if (i1.eq.-1.and.iz1.eq.-1.and. !aProton
     &   i2.eq.106.and.iz2.eq.1) then      !K+
        sig1 = STREXCHANGE(e,i1,iz1,i2,iz2,-40,2,101,-2) !aSigma+ + Pi-
        sig2 = STREXCHANGE(e,i1,iz1,i2,iz2,-40,-2,101,2) !aSigma- + Pi+
        sig3 = STREXCHANGE(e,i1,iz1,i2,iz2,-40,0,101,0)  !aLambda0/aSigma0 + Pi0
        call getbran((/sig1,sig2,sig3/),1,3,dummy,1,3,channel)
        i3=-40
        i4=101
        select case(channel)
          case(1) !aSigma+ + Pi-
            iz3=2
            iz4=-2
          case(2) !aSigma- + Pi+
            iz3=-2
            iz4=2
          case(3) !aLambda0/aSigma0 + Pi0
            iz3=0
            iz4=0
            if(ranff(0).gt.CTParam(72)) then !aLambda
              i3=-27
            endif
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c K+ + an -> aY + pi
      else if (i1.eq.-1.and.iz1.eq.1.and. !aNeutron
     &   i2.eq.106.and.iz2.eq.1) then     !K+
        sig1 = STREXCHANGE(e,i1,iz1,i2,iz2,-40,2,101,0)  !aSigma+ + Pi0
        sig2 = STREXCHANGE(e,i1,iz1,i2,iz2,-40,0,101,2)  !aLambda0/aSigma0 + Pi+
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        i3=-40
        i4=101
        select case(channel)
          case(1) !aSigma+ + Pi0
            iz3=2
            iz4=0
          case(2) !aLambda0/aSigma0 + Pi+
            iz3=0
            iz4=2
            if(ranff(0).gt.CTParam(72)) then !aLambda
              i3=-27
            endif
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c K0 + ap
      else if (i1.eq.-1.and.iz1.eq.-1.and.!aProton
     &   i2.eq.106.and.iz2.eq.-1)  then   !K0
        sig1 = STREXCHANGE(e,i1,iz1,i2,iz2,-40,-2,101,0) !aSigma- + Pi0
        sig2 = STREXCHANGE(e,i1,iz1,i2,iz2,-40,0,101,-2) !aLambda0/aSigma0 + Pi-
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        i3=-40
        i4=101
        select case(channel)
          case(1) !aSigma- + Pi0
            iz3=-2
            iz4=0
          case(2) !aLambda0/aSigma0 + Pi-
            iz3=0
            iz4=-2
            if(ranff(0).gt.CTParam(72)) then !aLambda
              i3=-27
            endif
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c K0 + an
      else if (i1.eq.-1.and.iz1.eq.1.and. !aNeutron
     &   i2.eq.106.and.iz2.eq.-1)  then   !K0
        sig1 = STREXCHANGE(e,i1,iz1,i2,iz2,-40,2,101,-2) !aSigma+ + Pi-
        sig2 = STREXCHANGE(e,i1,iz1,i2,iz2,-40,-2,101,2) !aSigma- + Pi+
        sig3 = STREXCHANGE(e,i1,iz1,i2,iz2,-40,0,101,0)  !aLambda0/aSigma0 + Pi0
        call getbran((/sig1,sig2,sig3/),1,3,dummy,1,3,channel)
        i3=-40
        i4=101
        select case(channel)
          case(1)
            iz3=2
            iz4=-2
          case(2)
            iz3=-2
            iz4=2
          case(3)
            iz3=0
            iz4=0
            if(ranff(0).gt.CTParam(72)) then !aLambda
              i3=-27
            endif
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c aK0 + p
      else if (i1.eq.1.and.iz1.eq.1.and.  !proton
     &   i2.eq.-106.and.iz2.eq.1)  then   !aK0
        sig1 = STREXCHANGE(e,i1,iz1,i2,iz2,40,2,101,0) !Sigma+ + Pi0
        sig2 = STREXCHANGE(e,i1,iz1,i2,iz2,40,0,101,2) !Lambda0/aSigma0 + Pi+
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        i3=40
        i4=101
        select case(channel)
          case(1)
            iz3=2
            iz4=0
          case(2)
            iz3=0
            iz4=2
            if(ranff(0).gt.CTParam(72)) then !Lambda
              i3=27
            endif
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c aK0 + n
      else if (i1.eq.1.and.iz1.eq.-1.and. !neutron
     &   i2.eq.-106.and.iz2.eq.1)  then   !aK0
        sig1 = STREXCHANGE(e,i1,iz1,i2,iz2,40,2,101,-2) !Sigma+ + Pi-
        sig2 = STREXCHANGE(e,i1,iz1,i2,iz2,40,-2,101,2) !Sigma- + Pi+
        sig3 = STREXCHANGE(e,i1,iz1,i2,iz2,40,0,101,0)  !Lambda0/Sigma0 + Pi0
        call getbran((/sig1,sig2,sig3/),1,3,dummy,1,3,channel)
        i3=40
        i4=101
        select case(channel)
          case(1)
            iz3=2
            iz4=-2
          case(2)
            iz3=-2
            iz4=2
          case(3)
            iz3=0
            iz4=0
            if(ranff(0).gt.CTParam(72)) then !Lambda
              i3=27
            endif
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
      else
        write(6,*) 'Error: Strangeness exchange should not be called 57'
        write(6,*) "i1:",i1,"iz1:",iz1,"i2:",i2,"iz2:",iz2
        stop 137
      end if

      m3=massit(i3)
      m4=massit(i4)

      goto 2001

c pi + Y -> K + N strangeness exchange
c all calculated via detailed balance
  58  continue
      nexit=2
cUNCHARGED CHANNELS
c Sigma+ + pi-   -> K-+p  & aK0+n
c and
c Sigma- + pi+   -> K-+p  & aK0+n
      if(i1.eq.40.and.abs(iz1).eq.2.and.       !charged Sigma
     &   i2.eq.101.and.abs(iz2).eq.2.and.      !charged Pi
     &   iz1+iz2.eq.0) then                    !total isospin=0
         i3=1
         i4=-106
         sig1 = STREXCHANGE(e,1,1,-106,-1,i1,iz1,i2,iz2)  !K- + p
     &          *pcms(e,massit(i3),massit(i4))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
         sig2 = STREXCHANGE(e,1,-1,-106,1,i1,iz1,i2,iz2)  !aK0 + n
     &          *pcms(e,massit(i3),massit(i4))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
         call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
         select case(channel)
           case(1)
             iz3=1
             iz4=-1
           case(2)
             iz3=-1
             iz4=1
           case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
         end select
c Sigma0 + pi0   -> K-+p  & aK0+n
      elseif(i1.eq.40.and.iz1.eq.0.and.        !Sigma0
     &       i2.eq.101.and.iz2.eq.0) then      !pi0
        i3=1
        i4=-106
        sig1 = STREXCHANGE(e,1,1,-106,-1,i1,iz1,i2,iz2)  !K- + p
     &         *CTParam(72)
     &         *pcms(e,massit(i3),massit(i4))**2
     &      *DIVpcms(e,massit(i1),massit(i2))**2
        sig2 = STREXCHANGE(e,1,-1,-106,1,i1,iz1,i2,iz2)  !aK0 + n
     &         *CTParam(72)
     &         *pcms(e,massit(i3),massit(i4))**2
     &      *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        select case(channel)
          case(1)
            iz3=1
            iz4=-1
          case(2)
            iz3=-1
            iz4=1
          case default
           write(6,*) 'Error: non-existent channel selected'
           write(6,*) i1,iz1,i2,iz2,channel
           stop 137
        end select
c Lambda0 + pi0  -> K-+p  & aK0+n
      elseif(i1.eq.27.and.iz1.eq.0.and.        !Lambda0
     &       i2.eq.101.and.iz2.eq.0) then      !pi0
        i3=1
        i4=-106
        sig1 = STREXCHANGE(e,1,1,-106,-1,i1,iz1,i2,iz2)  !K- + p
     &         *(1-CTParam(72))
     &         *pcms(e,massit(i3),massit(i4))**2
     &      *DIVpcms(e,massit(i1),massit(i2))**2
        sig2 = STREXCHANGE(e,1,-1,-106,1,i1,iz1,i2,iz2)  !aK0 + n
     &         *(1-CTParam(72))
     &         *pcms(e,massit(i3),massit(i4))**2
     &      *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        select case(channel)
          case(1)
            iz3=1
            iz4=-1
          case(2)
            iz3=-1
            iz4=1
          case default
           write(6,*) 'Error: non-existent channel selected'
           write(6,*) i1,iz1,i2,iz2,channel
           stop 137
        end select
c aSigma+ + pi-  -> K++ap & K0+an
c and
c aSigma- + pi+  -> K++ap & K0+an
      elseif(i1.eq.-40.and.abs(iz1).eq.2.and.      !charged aSigma
     &       i2.eq.101.and.abs(iz2).eq.2.and.      !charged Pi
     &       iz1+iz2.eq.0) then                    !total isospin=0
         i3=-1
         i4=106
         sig1 = STREXCHANGE(e,-1,-1,106,1,i1,iz1,i2,iz2)  !K+ + ap
     &          *pcms(e,massit(i3),massit(i4))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
         sig2 = STREXCHANGE(e,-1,1,106,-1,i1,iz1,i2,iz2)  !K0 + an
     &          *pcms(e,massit(i3),massit(i4))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
         call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
         select case(channel)
           case(1)
             iz3=-1
             iz4=1
           case(2)
             iz3=1
             iz4=-1
           case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
         end select
c aSigma0 + pi0  -> K++ap & K0+an
      elseif(i1.eq.-40.and.iz1.eq.0.and.       !aSigma0
     &       i2.eq.101.and.iz2.eq.0) then      !pi0
        i3=-1
        i4=106
        sig1 = STREXCHANGE(e,-1,-1,106,1,i1,iz1,i2,iz2)  !K+ + ap
     &         *CTParam(72)
     &         *pcms(e,massit(i3),massit(i4))**2
     &      *DIVpcms(e,massit(i1),massit(i2))**2
        sig2 = STREXCHANGE(e,-1,1,106,-1,i1,iz1,i2,iz2)  !K0 + an
     &         *CTParam(72)
     &         *pcms(e,massit(i3),massit(i4))**2
     &      *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        select case(channel)
          case(1)
            iz3=-1
            iz4=1
          case(2)
            iz3=1
            iz4=-1
          case default
           write(6,*) 'Error: non-existent channel selected'
           write(6,*) i1,iz1,i2,iz2,channel
           stop 137
        end select
c aLambda0 + pi0 -> K++ap & K0+an
      elseif(i1.eq.-27.and.iz1.eq.0.and.       !aLambda0
     &       i2.eq.101.and.iz2.eq.0) then      !pi0
        i3=-1
        i4=106
        sig1 = STREXCHANGE(e,-1,-1,106,1,i1,iz1,i2,iz2)  !K+ + ap
     &         *CTParam(72)
     &         *pcms(e,massit(i3),massit(i4))**2
     &      *DIVpcms(e,massit(i1),massit(i2))**2
        sig2 = STREXCHANGE(e,-1,1,106,-1,i1,iz1,i2,iz2)  !K0 + an
     &         *CTParam(72)
     &         *pcms(e,massit(i3),massit(i4))**2
     &      *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        select case(channel)
          case(1)
            iz3=-1
            iz4=1
          case(2)
            iz3=1
            iz4=-1
          case default
           write(6,*) 'Error: non-existent channel selected'
           write(6,*) i1,iz1,i2,iz2,channel
           stop 137
        end select
cCHARGED CHANNELS

c-----positive---------
c Sigma+ + pi0   -> aK0+p
      elseif(i1.eq.40.and.iz1.eq.2.and.      !Sigma+
     &       i2.eq.101.and.iz2.eq.0) then    !pi0
        i3=1
        iz3=1
        i4=-106
        iz4=1
c Sigma0 + pi+   -> aK0+p
      elseif(i1.eq.40.and.iz1.eq.0.and.      !Sigma0
     &       i2.eq.101.and.iz2.eq.2) then    !pi+
        i3=1
        iz3=1
        i4=-106
        iz4=1
c Lambda0 + pi+  -> aK0+p
      elseif(i1.eq.27.and.iz1.eq.0.and.      !Lambda0
     &       i2.eq.101.and.iz2.eq.2) then    !pi+
        i3=1
        iz3=1
        i4=-106
        iz4=1
c aSigma+ + pi0  -> K++an
      elseif(i1.eq.-40.and.iz1.eq.2.and.     !aSigma+
     &       i2.eq.101.and.iz2.eq.0) then    !pi0
        i3=-1
        iz3=1
        i4=106
        iz4=1
c aSigma0 + pi+  -> K++an
      elseif(i1.eq.-40.and.iz1.eq.0.and.     !aSigma0
     &       i2.eq.101.and.iz2.eq.2) then    !pi+
        i3=-1
        iz3=1
        i4=106
        iz4=1
c aLambda0 + pi+ -> K++an
      elseif(i1.eq.-27.and.iz1.eq.0.and.     !aLambda0
     &       i2.eq.101.and.iz2.eq.2) then    !pi+
        i3=-1
        iz3=1
        i4=106
        iz4=1
c-----negative-----------
c Sigma- + pi0   -> K-+n
      elseif(i1.eq.40.and.iz1.eq.-2.and.     !Sigma-
     &       i2.eq.101.and.iz2.eq.0) then    !pi0
        i3  = 1
        iz3 = -1
        i4  = -106
        iz4 = -1
c Sigma0 + pi-   -> K-+n
      elseif(i1.eq.40.and.iz1.eq.0.and.      !Sigma0
     &       i2.eq.101.and.iz2.eq.-2) then   !pi-
        i3  = 1
        iz3 = -1
        i4  = -106
        iz4 = -1
c Lambda0 + pi-  -> K-+n
      elseif(i1.eq.27.and.iz1.eq.0.and.      !Lambda0
     &       i2.eq.101.and.iz2.eq.-2) then   !pi-
        i3  = 1
        iz3 = -1
        i4  = -106
        iz4 = -1
c aSigma- + pi0  -> K0+ap
      elseif(i1.eq.-40.and.iz1.eq.-2.and.    !aSigma-
     &       i2.eq.101.and.iz2.eq.0) then    !pi0
        i3  = -1
        iz3 = -1
        i4  = 106
        iz4 = -1
c aSigma0 + pi-  -> K0+ap
      elseif(i1.eq.-40.and.iz1.eq.0.and.     !aSigma0
     &       i2.eq.101.and.iz2.eq.-2) then   !pi-
        i3  = -1
        iz3 = -1
        i4  = 106
        iz4 = -1
c aLambda0 + pi- -> K0+ap
      elseif(i1.eq.-27.and.iz1.eq.0.and.     !aLambda0
     &       i2.eq.101.and.iz2.eq.-2) then   !pi-
        i3  = -1
        iz3 = -1
        i4  = 106
        iz4 = -1
      else
        write(6,*) 'Error: Strangeness exchange should not be called 58'
        write(6,*) "i1:",i1,"iz1:",iz1,"i2:",i2,"iz2:",iz2
        stop 137
      end if

      m3=massit(i3)
      m4=massit(i4)
      goto 2001

  59  continue
c Kaon + strange Baryon strangeness exchange
      nexit=2
      i4=101
c Sigma+ + K- -> Xi0 + pi0
c             -> Xi- + pi+
      if(i1.eq.40.and.iz1.eq.2.and.        !Sigma+
     &   i2.eq.-106.and.iz2.eq.-1) then    !K-
        sig1=STREXCHANGE(e,i1,iz1,i2,iz2,49,1,101,0)
        sig2=STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,101,2)
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        i3=49
        select case(channel)
          case(1)
            iz3=1
            iz4=0
          case(2)
            iz3=-1
            iz4=2
          case default
           write(6,*) 'Error: non-existent channel selected'
           write(6,*) i1,iz1,i2,iz2,channel
           stop 137
        end select
c Sigma+ + aK0-> Xi0 + pi+
      else if(i1.eq.40.and.iz1.eq.2.and.   !Sigma+
     &   i2.eq.-106.and.iz2.eq.1) then     !aK0
        i3=49
        iz3=1
        iz4=2
c Lambda0/Sigma0 + aK0 -> Xi0 + pi0
c Lambda0/Sigma0 + aK0 -> Xi- + pi+
      else if((i1.eq.40.or.i1.eq.27).and.iz1.eq.0.and.  !Lambda0/Sigma0
     &   i2.eq.-106.and.iz2.eq.1) then                  !aK0
        sig1=STREXCHANGE(e,i1,iz1,i2,iz2,49,1,101,0)
        sig2=STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,101,2)
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        i3=49
        select case(channel)
          case(1)
            iz3=1
            iz4=0
          case(2)
            iz3=-1
            iz4=2
          case default
           write(6,*) 'Error: non-existent channel selected'
           write(6,*) i1,iz1,i2,iz2,channel
           stop 137
        end select
c Lambda0/Sigma0 + K-  -> Xi0 + pi-
c Lambda0/Sigma0 + K-  -> Xi- + pi0
      else if((i1.eq.40.or.i1.eq.27).and.iz1.eq.0.and.  !Lambda0/Sigma0
     &   i2.eq.-106.and.iz2.eq.-1) then                 !aK-
        sig1=STREXCHANGE(e,i1,iz1,i2,iz2,49,1,101,-2)
        sig2=STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,101,0)
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        i3=49
        select case(channel)
          case(1)
            iz3=1
            iz4=-2
          case(2)
            iz3=-1
            iz4=0
          case default
           write(6,*) 'Error: non-existent channel selected'
           write(6,*) i1,iz1,i2,iz2,channel
           stop 137
        end select
c Sigma- + aK0 -> Xi- + pi0
c              -> Xi0 + pi-
      else if(i1.eq.40.and.iz1.eq.-2.and. !Sigma-
     &   i2.eq.-106.and.iz2.eq.1) then    !aK0
        sig1=STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,101,0)
        sig2=STREXCHANGE(e,i1,iz1,i2,iz2,49,1,101,-2)
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        i3=49
        select case(channel)
          case(1)
            iz3=-1
            iz4=0
          case(2)
            iz3=1
            iz4=-2
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c Sigma- + K-  -> Xi- + pi-
      else if(i1.eq.40.and.iz1.eq.-2.and. !Sigma-
     &   i2.eq.-106.and.iz2.eq.-1) then   !aK-
        i3=49
        iz3=-1
        iz4=-2
c aSigma+ + K0 -> aXi+ + pi0
c              -> aXi0 + pi+
      else if(i1.eq.-40.and.iz1.eq.2.and. !aSigma+
     &   i2.eq.106.and.iz2.eq.-1) then    !K0
        sig1=STREXCHANGE(e,i1,iz1,i2,iz2,-49,1,101,0)
        sig2=STREXCHANGE(e,i1,iz1,i2,iz2,-49,-1,101,2)
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        i3=-49
        select case(channel)
          case(1)
            iz3=1
            iz4=0
          case(2)
            iz3=-1
            iz4=2
          case default
           write(6,*) 'Error: non-existent channel selected'
           write(6,*) i1,iz1,i2,iz2,channel
           stop 137
        end select
c aSigma+ + K+ -> aXi+ + pi+
      else if(i1.eq.-40.and.iz1.eq.2.and. !aSigma+
     &   i2.eq.106.and.iz2.eq.1) then     !K+
        i3=-49
        iz3=1
        iz4=2
c aLambda0/aSigma0 + K0 -> aXi0 + pi0
c aLambda0/aSigma0 + K0 -> aXi+ + pi-
      else if((i1.eq.-40.or.i1.eq.-27).and.iz1.eq.0.and.  !aLambda0/aSigma0
     &   i2.eq.106.and.iz2.eq.-1) then                    !K0
        sig1=STREXCHANGE(e,i1,iz1,i2,iz2,-49,-1,101,0)
        sig2=STREXCHANGE(e,i1,iz1,i2,iz2,-49,1,101,-2)
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        i3=-49
        select case(channel)
          case(1)
            iz3=-1
            iz4=0
          case(2)
            iz3=1
            iz4=-2
          case default
           write(6,*) 'Error: non-existent channel selected'
           write(6,*) i1,iz1,i2,iz2,channel
           stop 137
        end select
c aLambda0/aSigma0 + K+ -> aXi0 + pi+
c aLambda0/aSigma0 + K+ -> aXi+ + pi0
      else if((i1.eq.-40.or.i1.eq.-27).and.iz1.eq.0.and.  !aLambda0/aSigma0
     &   i2.eq.106.and.iz2.eq.1) then                     !K+
        sig1=STREXCHANGE(e,i1,iz1,i2,iz2,-49,-1,101,2)
        sig2=STREXCHANGE(e,i1,iz1,i2,iz2,-49,1,101,0)
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        i3=-49
        select case(channel)
          case(1)
            iz3=-1
            iz4=2
          case(2)
            iz3=1
            iz4=0
          case default
           write(6,*) 'Error: non-existent channel selected'
           write(6,*) i1,iz1,i2,iz2,channel
           stop 137
        end select
c aSigma- + K0 -> aXi0 + pi-
      else if(i1.eq.-40.and.iz1.eq.-2.and. !aSigma-
     &   i2.eq.106.and.iz2.eq.-1) then     !K0
        i3=-49
        iz3=-1
        iz4=-2
c aSigma- + K+ -> aXi0 + pi0
c              -> aXi+ + pi-
      else if(i1.eq.-40.and.iz1.eq.-2.and. !aSigma-
     &   i2.eq.106.and.iz2.eq.1) then      !K+
        sig1=STREXCHANGE(e,i1,iz1,i2,iz2,-49,-1,101,0)
        sig2=STREXCHANGE(e,i1,iz1,i2,iz2,-49,1,101,-2)
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        i3=-49
        select case(channel)
          case(1)
            iz3=-1
            iz4=0
          case(2)
            iz3=1
            iz4=-2
          case default
           write(6,*) 'Error: non-existent channel selected'
           write(6,*) i1,iz1,i2,iz2,channel
           stop 137
        end select
      else
        write(6,*) 'Error: Strangeness exchange should not be called 59'
        write(6,*) "i1:",i1,"iz1:",iz1,"i2:",i2,"iz2:",iz2
        stop 137
      end if


      m3=massit(i3)
      m4=massit(i4)
      goto 2001

c Xi + pi -> Sigma/Lambda + K
  60  continue
      nexit=2
c Xi0 + pi0 -> K- + Sigma+
c Xi0 + pi0 -> aK0+ Lambda0/Sigma0
      if(i1.eq.49.and.iz1.eq.1.and.        !Xi0
     &   i2.eq.101.and.iz2.eq.0) then      !pi0
        sig1=STREXCHANGE(e,40,2,-106,-1,i1,iz1,i2,iz2)
     *      *pcms(e,massit(-106),massit(40))**2
     *   *DIVpcms(e,massit(i1),massit(i2))**2
        sig2=STREXCHANGE(e,40,0,-106,1,i1,iz1,i2,iz2)
     *      *pcms(e,massit(-106),massit(40))**2
     *   *DIVpcms(e,massit(i1),massit(i2))**2
        sig3=STREXCHANGE(e,27,0,-106,1,i1,iz1,i2,iz2)
     *      *pcms(e,massit(-106),massit(27))**2
     *   *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2,sig3/),1,3,dummy,1,3,channel)
        select case(channel)
          case(1)
            i3 = 40
            iz3= 2
            i4 = -106
            iz4= -1
          case(2)
            i3 = 40
            iz3= 0
            i4 = -106
            iz4= 1
          case(3)
            i3 = 27
            iz3= 0
            i4 = -106
            iz4= 1
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c Xi0 + pi- -> K- + Lambda0/Sigma0
c Xi0 + pi- -> aK0+ Sigma-
      elseif(i1.eq.49.and.iz1.eq.1.and.    !Xi0
     &   i2.eq.101.and.iz2.eq.-2) then     !pi-
        i4 =-106
        sig1=STREXCHANGE(e,40,0,-106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(-106),massit(40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig2=STREXCHANGE(e,27,0,-106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(-106),massit(27))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig3=STREXCHANGE(e,40,-2,-106,1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(-106),massit(40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2,sig3/),1,3,dummy,1,3,channel)
        select case(channel)
          case(1)
            i3 = 40
            iz3=0
            iz4=-1
          case(2)
            i3 = 27
            iz3=0
            iz4=-1
          case(3)
            i3 = 40
            iz3= -2
            iz4=1
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c Xi0 + pi+ -> aK0+ Sigma+
      elseif(i1.eq.49.and.iz1.eq.1.and.    !Xi0
     &   i2.eq.101.and.iz2.eq.2) then      !pi+
        i3 =40
        iz3=2
        i4 =-106
        iz4=1
c Xi- + pi+ -> aK0 + Lambda0/Sigma0
c Xi- + pi+ -> K-  + Sigma+
      elseif(i1.eq.49.and.iz1.eq.-1.and.   !Xi-
     &   i2.eq.101.and.iz2.eq.2) then      !pi+
        i4 =-106
        sig1=STREXCHANGE(e,40,0,-106,1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(-106),massit(40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig2=STREXCHANGE(e,27,0,-106,1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(-106),massit(27))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig3=STREXCHANGE(e,40,2,-106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(-106),massit(40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2,sig3/),1,3,dummy,1,3,channel)
        select case(channel)
          case(1)
            i3 = 40
            iz3=0
            iz4=1
          case(2)
            i3 = 27
            iz3=0
            iz4=1
          case(3)
            i3 = 40
            iz3= 2
            iz4=-1
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c Xi- + pi0 -> aK0 + Sigma-
c Xi- + pi0 -> K-  + Lambda0/Sigma0
      elseif(i1.eq.49.and.iz1.eq.-1.and.   !Xi-
     &   i2.eq.101.and.iz2.eq.0) then      !pi0
        sig1=STREXCHANGE(e,40,-2,-106,1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(-106),massit(40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig2=STREXCHANGE(e,40,0,-106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(-106),massit(40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig3=STREXCHANGE(e,27,0,-106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(-106),massit(27))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2,sig3/),1,3,dummy,1,3,channel)
        select case(channel)
          case(1)
            i3 = 40
            iz3= -2
            i4 = -106
            iz4= 1
          case(2)
            i3 = 40
            iz3= 0
            i4 = -106
            iz4= -1
          case(3)
            i3 = 27
            iz3= 0
            i4 = -106
            iz4= -1
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c Xi- + pi- -> K-  + Sigma-
      elseif(i1.eq.49.and.iz1.eq.-1.and.   !Xi-
     &   i2.eq.101.and.iz2.eq.-2) then     !pi-
        i3 =40
        iz3=-2
        i4 =-106
        iz4=-1
c aXi0 + pi0 -> K+ + aSigma-
c aXi0 + pi0 -> K0 + aLambda0/aSigma0
      elseif(i1.eq.-49.and.iz1.eq.-1.and.  !Xi0
     &   i2.eq.101.and.iz2.eq.0) then      !pi0
        sig1=STREXCHANGE(e,-40,-2,106,1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig2=STREXCHANGE(e,-40,0,106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig3=STREXCHANGE(e,-27,0,106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-27))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2,sig3/),1,3,dummy,1,3,channel)
        select case(channel)
          case(1)
            i3 = -40
            iz3= -2
            i4 = 106
            iz4= 1
          case(2)
            i3 = -40
            iz3= 0
            i4 = 106
            iz4= -1
          case(3)
            i3 = -27
            iz3= 0
            i4 = 106
            iz4= -1
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c aXi0 + pi+ -> K+ + aLambda0/aSigma0
c aXi0 + pi+ -> K0 + aSigma+
      elseif(i1.eq.-49.and.iz1.eq.-1.and.  !Xi0
     &   i2.eq.101.and.iz2.eq.2) then      !pi+
        i4 =106
        sig1=STREXCHANGE(e,-40,0,106,1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig2=STREXCHANGE(e,-27,0,106,1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-27))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig3=STREXCHANGE(e,-40,2,106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2,sig3/),1,3,dummy,1,3,channel)
        select case(channel)
          case(1)
            i3 = -40
            iz3=0
            iz4=1
          case(2)
            i3 = -27
            iz3=0
            iz4=1
          case(3)
            i3 = -40
            iz3=2
            iz4=-1
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c aXi0 + pi- -> K0 + aSigma-
      elseif(i1.eq.-49.and.iz1.eq.-1.and.  !Xi0
     &   i2.eq.101.and.iz2.eq.-2) then     !pi-
        i3 =-40
        iz3=-2
        i4 =106
        iz4=-1
c aXi+ + pi- -> K0 + aLambda0/aSigma0
c aXi+ + pi- -> K+ + aSigma-
      elseif(i1.eq.-49.and.iz1.eq.1.and.   !Xi+
     &   i2.eq.101.and.iz2.eq.-2) then     !pi-
        i4 =106
        sig1=STREXCHANGE(e,-40,0,106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig2=STREXCHANGE(e,-27,0,106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-27))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig3=STREXCHANGE(e,-40,-2,106,1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2,sig3/),1,3,dummy,1,3,channel)
        select case(channel)
          case(1)
            i3 = -40
            iz3=0
            iz4=-1
          case(2)
            i3 = -27
            iz3=0
            iz4=-1
          case(3)
            i3 = -40
            iz3=-2
            iz4=1
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c aXi+ + pi0 -> K0 + aSigma+
c aXi+ + pi0 -> K+ + aLambda0/aSigma0
      elseif(i1.eq.-49.and.iz1.eq.1.and.   !aXi+
     &   i2.eq.101.and.iz2.eq.0) then      !pi0
        sig1=STREXCHANGE(e,-40,2,106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig2=STREXCHANGE(e,-40,0,106,1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-40))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig3=STREXCHANGE(e,-27,0,106,1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-27))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2,sig3/),1,3,dummy,1,3,channel)
        select case(channel)
          case(1)
            i3 = -40
            iz3= 2
            i4 = 106
            iz4= -1
          case(2)
            i3 = -40
            iz3= 0
            i4 = 106
            iz4= 1
          case(3)
            i3 = -27
            iz3= 0
            i4 = 106
            iz4= 1
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c aXi+ + pi+ -> K+ + aSigma+
      elseif(i1.eq.-49.and.iz1.eq.1.and.   !aXi+
     &   i2.eq.101.and.iz2.eq.2) then      !pi+
        i3 =-40
        iz3=2
        i4 =106
        iz4=1
      endif

      m3=massit(i3)
      m4=massit(i4)
      goto 2001


c Xi + K -> Omega + pi Strangness exchange
61    continue

c Xi0 + aK0 -> Omega- + pi+
      if(i1.eq.49.and.iz1.eq.1.and.
     &   i2.eq.-106.and.iz2.eq.1) then
        i3  = 55
        iz3 = 0
        i4  = 101
        iz4 = 2
c Xi0 + K- -> Omega- + pi0
c Xi- + aK0 -> Omega- + pi0
      elseif((i1.eq.49.and.iz1.eq.1.and.
     &   i2.eq.-106.and.iz2.eq.-1).or.
     &  (i1.eq.49.and.iz1.eq.-1.and.
     &   i2.eq.-106.and.iz2.eq.1)) then
        i3  = 55
        iz3 = 0
        i4  = 101
        iz4 = 0
c Xi- + K- ->  Omega- + pi-
      elseif(i1.eq.49.and.iz1.eq.-1.and.
     &   i2.eq.-106.and.iz2.eq.-1) then
        i3  = 55
        iz3 = 0
        i4  = 101
        iz4 = -2
c aXi0 + K0 -> aOmega+ + pi-
      elseif(i1.eq.-49.and.iz1.eq.-1.and.
     &   i2.eq.106.and.iz2.eq.-1) then
        i3  = -55
        iz3 = 0
        i4  = 101
        iz4 = -2
c aXi0 + K+ -> aOmega+ + pi0
c aXi+ + K0 -> aOmega+ + pi0
      elseif((i1.eq.-49.and.iz1.eq.-1.and.
     &   i2.eq.106.and.iz2.eq.1).or.
     &  (i1.eq.-49.and.iz1.eq.1.and.
     &   i2.eq.106.and.iz2.eq.-1)) then
        i3  = -55
        iz3 = 0
        i4  = 101
        iz4 = 0
c aXi+ + K+ -> aOmega+ + pi+
      elseif(i1.eq.-49.and.iz1.eq.1.and.
     &   i2.eq.106.and.iz2.eq.1) then
        i3  = -55
        iz3 = 0
        i4  = 101
        iz4 = 2
      endif

      m3=massit(i3)
      m4=massit(i4)
      goto 2001

c Omega + pi -> Xi + K Strangeness exchange
62    continue
c Omega- + pi+ -> Xi0 + aK0
      if(i1.eq.55.and.iz1.eq.0.and.
     &   i2.eq.101.and.iz2.eq.2) then
        i3  = 49
        iz3 = 1
        i4  = -106
        iz4 = 1
c Omega- + pi0 -> Xi0 + K-
c Omega- + pi0 -> Xi- + aK0
      elseif(i1.eq.55.and.iz1.eq.0.and.
     &   i2.eq.101.and.iz2.eq.0) then
        i3  = 49
        i4  = -106
        sig1=5d-1*STREXCHANGE(e,49,1,-106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(-106),massit(49))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig2=5d-1*STREXCHANGE(e,49,-1,-106,1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(-106),massit(49))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        select case(channel)
          case(1)
            iz3 = 1
            iz4 = -1
          case(2)
            iz3 = -1
            iz4 = 1
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c Omega- + pi- -> Xi- + K-
      elseif(i1.eq.55.and.iz1.eq.0.and.
     &   i2.eq.101.and.iz2.eq.-2) then
        i3  = 49
        iz3 = -1
        i4  = -106
        iz4 = -1
c aOmega+ + pi- -> aXi0 + K0
      elseif(i1.eq.-55.and.iz1.eq.0.and.
     &   i2.eq.101.and.iz2.eq.-2) then
        i3  = -49
        iz3 = -1
        i4  = 106
        iz4 = -1
c aOmega+ + pi0 -> aXi0 + K+
c aOmega+ + pi0 -> aXi+ + K0
      elseif(i1.eq.-55.and.iz1.eq.0.and.
     &   i2.eq.101.and.iz2.eq.0) then
        i3  = -49
        i4  = 106
        sig1=5d-1*STREXCHANGE(e,-49,-1,106,1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-49))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        sig2=5d-1*STREXCHANGE(e,-49,1,106,-1,i1,iz1,i2,iz2)
     &      *pcms(e,massit(106),massit(-49))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        select case(channel)
          case(1)
            iz3 = -1
            iz4 = 1
          case(2)
            iz3 = 1
            iz4 = -1
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c aOmega+ + pi+ -> aXi+ + K+
      elseif(i1.eq.-55.and.iz1.eq.0.and.
     &   i2.eq.101.and.iz2.eq.2) then
        i3  = -49
        iz3 = 1
        i4  = 106
        iz4 = 1
      end if

      m3=massit(i3)
      m4=massit(i4)
      goto 2001


  63  continue

c Lambda0 + Lambda0 -> Xi0 + n
c                   -> Xi- + p
c Lambda0 + Sigma0  -> Xi0 + n
c                   -> Xi- + p
c Sigma0  + Sigma0  -> Xi0 + n
c                   -> Xi- + p
c Sigma+  + Sigma-  -> Xi0 + n
c                   -> Xi- + p
      if((iz1.eq.0.and.iz2.eq.0.and.  !uncharged channels
     &    (i1.eq.27.or.i1.eq.40).and.
     &    (i2.eq.27.or.i2.eq.40)).or.
     &   ((i1.eq.40.and.i2.eq.40).and.   !Sigma+ + Sigma-
     &     abs(iz1).eq.2.and.(iz1+iz2).eq.0)) then
        i3=49
        i4=1
        if(ranff(0).gt.0.5) then
          iz3=1
          iz4=-1
        else
          iz3=-1
          iz4=1
        endif

c Lambda0 + Sigma+  -> Xi0 + p
c Sigma0  + Sigma+  -> Xi0 + p
      elseif( (iz1.eq.0.and.iz2.eq.2.and.i1.eq.27.and.i2.eq.40).or.
     &        (iz1.eq.2.and.iz2.eq.0.and.i1.eq.40.and.i2.eq.27).or.
     &        (i1.eq.40.and.i2.eq.40.and.((iz1.eq.2.and.iz2.eq.0).or.
     &                                    (iz1.eq.0.and.iz2.eq.2))))then
        i3=49
        iz3=1
        i4=1
        iz4=1


c Lambda0 + Sigma-  -> Xi- + n
c Sigma0  + Sigma-  -> Xi- + n
      elseif( (iz1.eq.0.and.iz2.eq.-2.and.i1.eq.27.and.i2.eq.40).or.
     &        (iz1.eq.-2.and.iz2.eq.0.and.i1.eq.40.and.i2.eq.27).or.
     &        (i1.eq.40.and.i2.eq.40.and.((iz1.eq.-2.and.iz2.eq.0).or.
     &                                   (iz1.eq.0.and.iz2.eq.-2))))then
        i3=49
        iz3=-1
        i4=1
        iz4=-1

c aLambda0 + aLambda0 -> aXi0 + an
c                     -> aXi+ + ap
c aLambda0 + aSigma0  -> aXi0 + an
c                     -> aXi+ + ap
c aSigma0  + aSigma0  -> aXi0 + an
c                     -> aXi+ + ap
c aSigma-  + aSigma+  -> aXi0 + an
c                     -> aXi+ + ap
      elseif((iz1.eq.0.and.iz2.eq.0.and.  !uncharged channels
     &        (i1.eq.-27.or.i1.eq.-40).and.
     &        (i2.eq.-27.or.i2.eq.-40)).or.
     &       ((i1.eq.-40.and.i2.eq.-40).and.   !Sigma+ + Sigma-
     &     abs(iz1).eq.2.and.(iz1+iz2).eq.0)) then
        i3=-49
        i4=-1
        if(ranff(0).gt.0.5) then
          iz3=-1
          iz4=1
        else
          iz3=1
          iz4=-1
        endif

c aLambda0 + aSigma-  -> aXi0 + ap
c aSigma0  + aSigma-  -> aXi0 + ap
      elseif( (iz1.eq.0.and.iz2.eq.-2.and.i1.eq.-27.and.i2.eq.-40).or.
     &        (iz1.eq.-2.and.iz2.eq.0.and.i1.eq.-40.and.i2.eq.-27).or.
     &        (i1.eq.-40.and.i2.eq.-40.and.((iz1.eq.-2.and.iz2.eq.0).or.
     &                                   (iz1.eq.0.and.iz2.eq.-2))))then
        i3=-49
        iz3=-1
        i4=-1
        iz4=-1

c aLambda0 + aSigma+  -> aXi+ + an
c aSigma0  + aSigma+  -> aXi+ + an
      elseif( (iz1.eq.0.and.iz2.eq.2.and.i1.eq.-27.and.i2.eq.-40).or.
     &        (iz1.eq.2.and.iz2.eq.0.and.i1.eq.-40.and.i2.eq.-27).or.
     &        (i1.eq.-40.and.i2.eq.-40.and.((iz1.eq.2.and.iz2.eq.0).or.
     &                                    (iz1.eq.0.and.iz2.eq.2))))then
        i3=-49
        iz3=1
        i4=-1
        iz4=1
      endif

      m3=massit(i3)
      m4=massit(i4)
      goto 2001

  64  continue

c put the lower ityp in the first slot for convenience in
c the later if clauses
      si1=i1
      siz1=iz1
      si2=i2
      siz2=iz2
      if(abs(si1).gt.abs(si2)) then
        call swpizm(si1,siz1,rdummy,si2,siz2,rdummy)
      endif

c Xi0 + n -> Lambda0 + Lambda0
c Xi0 + n -> Lambda0 + Sigma0
c Xi0 + n -> Sigma0  + Sigma0
c Xi0 + n -> Sigma+  + Sigma-
c Xi- + p -> Lambda0 + Lambda0
c Xi- + p -> Lambda0 + Sigma0
c Xi- + p -> Sigma0  + Sigma0
c Xi- + p -> Sigma+  + Sigma-
      if((si1.eq.1.and.siz1.eq.-1.and.
     &    si2.eq.49.and.siz2.eq.1).or.
     &   (si1.eq.1.and.siz1.eq.1.and.
     &    si2.eq.49.and.siz2.eq.-1)) then
        sig1=STREXCHANGE(e,si1,siz1,si2,siz2,27,0,27,0) !Lam + Lam
        sig2=STREXCHANGE(e,si1,siz1,si2,siz2,27,0,40,0) !Lam + Sig0
        sig3=STREXCHANGE(e,si1,siz1,si2,siz2,40,0,40,0) !Sig0 + Sig0
        sig4=STREXCHANGE(e,si1,siz1,si2,siz2,40,2,40,-2)!Sig+ + Sig-
        call getbran((/sig1,sig2,sig3,sig4/),1,4,dummy,1,4,channel)
        select case(channel)
          case(1)
            i3  = 27
            iz3 = 0
            i4  = 27
            iz4 = 0
          case(2)
            i3  = 27
            iz3 = 0
            i4  = 40
            iz4 = 0
          case(3)
            i3  = 40
            iz3 = 0
            i4  = 40
            iz4 = 0
          case(4)
            i3  = 40
            iz3 = 2
            i4  = 40
            iz4 = -2
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select

c Xi0 + p -> Lambda0 + Sigma+
c Xi0 + p -> Sigma0  + Sigma+
      elseif(si1.eq.1.and.siz1.eq.1.and.
     &       si2.eq.49.and.siz2.eq.1) then
        sig1=STREXCHANGE(e,si1,siz1,si2,siz2,27,0,40,2) !Lam + Sig+
        sig2=STREXCHANGE(e,si1,siz1,si2,siz2,40,0,40,2) !Sig0 + Sig+
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        select case(channel)
          case(1)
            i3  = 27
            iz3 = 0
            i4  = 40
            iz4 = 2
          case(2)
            i3  = 40
            iz3 = 0
            i4  = 40
            iz4 = 2
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select

c Xi- + n -> Lambda0 + Sigma-
c Xi- + n -> Sigma0  + Sigma-
      elseif(si1.eq.1.and.siz1.eq.-1.and.
     &       si2.eq.49.and.siz2.eq.-1) then
        sig1=STREXCHANGE(e,si1,siz1,si2,siz2,27,0,40,-2) !Lam + Sig-
        sig2=STREXCHANGE(e,si1,siz1,si2,siz2,40,0,40,2)  !Sig0 + Sig-
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        select case(channel)
          case(1)
            i3  = 27
            iz3 = 0
            i4  = 40
            iz4 = -2
          case(2)
            i3  = 40
            iz3 = 0
            i4  = 40
            iz4 = -2
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
c aXi0 + an -> aLambda0 + aLambda0
c aXi0 + an -> aLambda0 + aSigma0
c aXi0 + an -> aSigma0  + aSigma0
c aXi0 + an -> aSigma-  + aSigma+
c aXi+ + ap -> aLambda0 + aLambda0
c aXi+ + ap -> aLambda0 + aSigma0
c aXi+ + ap -> aSigma0  + aSigma0
c aXi+ + ap -> aSigma-  + aSigma+
      elseif((si1.eq.-1.and.siz1.eq.1.and.
     &        si2.eq.-49.and.siz2.eq.-1).or.
     &       (si1.eq.-1.and.siz1.eq.-1.and.
     &        si2.eq.-49.and.siz2.eq.1)) then
        sig1=STREXCHANGE(e,si1,siz1,si2,siz2,-27,0,-27,0)  !aLam + aLam
        sig2=STREXCHANGE(e,si1,siz1,si2,siz2,-27,0,-40,0) !aLam + aSig0
        sig3=STREXCHANGE(e,si1,siz1,si2,siz2,-40,0,-40,0)  !aSig0 + aSig0
        sig4=STREXCHANGE(e,si1,siz1,si2,siz2,-40,-2,-40,2) !aSig- + aSig+
        call getbran((/sig1,sig2,sig3,sig4/),1,4,dummy,1,4,channel)
        select case(channel)
          case(1)
            i3  = -27
            iz3 = 0
            i4  = -27
            iz4 = 0
          case(2)
            i3  = -27
            iz3 = 0
            i4  = -40
            iz4 = 0
          case(3)
            i3  = -40
            iz3 = 0
            i4  = -40
            iz4 = 0
          case(4)
            i3  = -40
            iz3 = -2
            i4  = -40
            iz4 = 2
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select

c aXi0 + ap -> aLambda0 + aSigma+
c aXi0 + ap -> aSigma0  + aSigma+
      elseif(si1.eq.-1.and.siz1.eq.-1.and.
     &       si2.eq.-49.and.siz2.eq.-1) then
        sig1=STREXCHANGE(e,si1,siz1,si2,siz2,-27,0,-40,-2) !aLam + aSig-
        sig2=STREXCHANGE(e,si1,siz1,si2,siz2,-40,0,-40,-2) !aSig0 + aSig-
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        select case(channel)
          case(1)
            i3  = -27
            iz3 = 0
            i4  = -40
            iz4 = -2
          case(2)
            i3  = -40
            iz3 = 0
            i4  = -40
            iz4 = -2
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select

c aXi- + an -> aLambda0 + aSigma+
c aXi- + an -> aSigma0  + aSigma+
      elseif(si1.eq.-1.and.siz1.eq.1.and.
     &       si2.eq.-49.and.siz2.eq.1) then
        sig1=STREXCHANGE(e,si1,siz1,si2,siz2,-27,0,-40,2) !aLam + aSig+
        sig2=STREXCHANGE(e,si1,siz1,si2,siz2,-40,0,-40,-2)!aSig0 + aSig+
        call getbran((/sig1,sig2/),1,2,dummy,1,2,channel)
        select case(channel)
          case(1)
            i3  = -27
            iz3 = 0
            i4  = -40
            iz4 = 2
          case(2)
            i3  = -40
            iz3 = 0
            i4  = -40
            iz4 = 2
          case default
            write(6,*) 'Error: non-existent channel selected'
            write(6,*) i1,iz1,i2,iz2,channel
            stop 137
        end select
      end if

      m3=massit(i3)
      m4=massit(i4)
      goto 2001

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


 1008 continue
c...get isospin-3 components

      nexit=2
      itot(1)=isoit(i3)
      itot(2)=isoit(i4)

      call isocgk4(isoit(i1),iz1,isoit(i2),iz2,itot,i3new,errflg)
      iz3=i3new(1)
      iz4=i3new(2)


      if(errflg.ne.0)then
        write(6,*)'make22: iso-spin conservation ',
     ,            'not possible in isocgk: error-flag=',errflg
        write(6,*)'      ',isoit(i1),iz1,isoit(i2),iz2,'>',
     ,            isoit(i3),isoit(i4),iz3,iz4,
     ,            ' process:',e,i1,m1,i2,m2,'>',i3,m3,i4,m4,'io=',io
      end if

      if(io.lt.40) then
         i3=isign(i3,ii1)
         i4=isign(i4,ii2)
      endif

 2001 continue


c... get momenta & fill newpart, 2 particle exit-channel

      nstring1=1
      nstring2=1
      nexit=2
      do i=1,4
         do j=1,2
            pnew(i,j)=0d0
            xnew(i,j)=0d0
         end do
      end do

c...boost to 2-particle cms


      if(.not.(ityptd(1,pslot(1)).ne.0.and.CTOption(34).eq.2.and.
     &     iline.eq.20)) then
c normal decay

         pnew(3,1)=pcms(e,m3,m4)
         pnew(3,2)=-pcms(e,m3,m4)

         pnew(4,1)=sqrt(m3**2+pnew(3,1)**2)
         pnew(4,2)=sqrt(m4**2+pnew(3,2)**2)
      else
c forward time delay
         pnew(1,1)=pxtd(1,pslot(1))
         pnew(1,2)=pxtd(2,pslot(1))
         pnew(2,1)=pytd(1,pslot(1))
         pnew(2,2)=pytd(2,pslot(1))
         pnew(3,1)=pztd(1,pslot(1))
         pnew(3,2)=pztd(2,pslot(1))
         pnew(4,1)=p0td(1,pslot(1))
         pnew(4,2)=p0td(2,pslot(1))
      endif

      pnew(5,1)=m3
      mstring(1)=m3
      itypnew(1)=i3
      i3new(1)=iz3

      pnew(5,2)=m4
      mstring(2)=m4
      itypnew(2)=i4
      i3new(2)=iz4

      return

 2002 continue
c...  fill newpart, one particle exit channel
      nstring1=1
      nstring2=0
      nexit=1
      do i=1,4
        pnew(i,1)=0d0
        xnew(i,1)=0d0
      end do
      pnew(4,1)=e
      mstring(1)=e
c     the rest of the relevant new particle data have been
c     filled into the newpart arrays by anndex

      return

 2003 continue
c bookkeeping for forward time-delay
      do 204 j=1,2
         if(pslot(j).lt.1) goto 204
         pold(1,j)=px(pslot(j))
         pold(2,j)=py(pslot(j))
         pold(3,j)=pz(pslot(j))
         pold(4,j)=p0(pslot(j))
         pold(5,j)=fmass(pslot(j))
         itypold(j)=ityp(pslot(j))
         iso3old(j)=iso3(pslot(j))
 204  continue

c...  fill newpart, one particle exit channel
      nstring1=1
      nstring2=0
      nexit=1
      do i=1,4
        pnew(i,1)=0d0
        xnew(i,1)=0d0
      end do
      pnew(4,1)=e
      mstring(1)=e

      return
      end


C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function rnfxf(b,xm,xp)
c
cinput b  : parameter
cinput xm : lower interval boundary
cinput xp : upper interval boundary
c
c {\tt rnfxf} yields a value $x$ between {\tt xm} and {\tt xp}
c distributetd like $(1-x)^b/x$.
c
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit none
      real*8 xm,xp,i,inv,x,fmin,fmax,f,ranff,b

      i(x)=log(x)
      inv(x)=exp(x)

      fmin=i(xm)
      fmax=i(xp)
 3    f=fmin+(fmax-fmin)*ranff(0)
      x=inv(f)
      if(ranff(0).gt.(1d0-x/xp)**b)goto 3
      rnfxf=x

      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine setizm(i1,iz1,m1,i2,iz2,m2,i3,iz3,m3,i4,iz4,m4)
c
c input  : {\tt i1,iz1,m1,i2,iz2,m2}
c output : {\tt i3,iz3,m3,i4,iz4,m4}
c
c This subroutine simply maps {\tt iz1} $\to$ {\tt iz3} $\ldots$
c {\tt m2} $\to$ {\tt m4}
c
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      real*8 m1,m2,m3,m4
      integer i1,iz1,i2,iz2,i3,iz3,i4,iz4
      i3=i1
      i4=i2
      m3=m1
      m4=m2
      iz3=iz1
      iz4=iz2
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine swpizm(i1,iz1,m1,i2,iz2,m2)
c
c
c input   : {\tt i1,iz1,m1,i2,iz2,m2}
c output  : {\tt i1,iz1,m1,i2,iz2,m2}
c
c This subroutine simply swaps {\tt 1} $\to$ {\tt 2} and vice versa.
c
c
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      real*8 m1,m2,ms
      integer i1,iz1,i2,iz2,is,izs
      is=i1
      i1=i2
      ms=m1
      m1=m2
      izs=iz1
      iz1=iz2
      i2=is
      iz2=izs
      m2=ms
      return
      end


C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function rnfpow(n,mmin,mmax)
c
cinput n    : parameter
cinput mmin : lower interval boundary
cinput mmax : upper interval boundary
c
c {\tt rnfpow} yields a value $x$ between {\tt mmin} and {\tt mmax}
c distributetd like $x^n$.
c
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit none
      real*8 mmin,mmax,i,inv,ii,iinv,x,fmin,fmax,f,ranff
      integer n

      i(x)=x**(n+1)/(n+1)
      inv(x)=((n+1)*x)**(1/(n+1))
      ii(x)=log(x)
      iinv(x)=exp(x)

      if(n.eq.-1)then
        fmin=ii(mmin)
        fmax=ii(mmax)
        f=fmin+(fmax-fmin)*ranff(0)
        rnfpow=iinv(f)
      else
        fmin=i(mmin)
        fmax=i(mmax)
        f=fmin+(fmax-fmin)*ranff(0)
        rnfpow=inv(f)
      end if
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function fmsr(mmin,mmax)
c
cinput mmin : minimum mass
cinput mmax : maximum mass
c
c  {\tt fmsr} yields a mass according to the finite mass sum rule (FMSR)
c  A.I.Sanda Phys.~Rev.~{\bf D6}~(1973)~231 and
c  M.B.Einhorn  et al., Phys. Rev. {\bf D5}~(1972)~2063
c
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit none
      real*8 mmin,mmax,i,inv,x,rnfpow
      i(x)=log(x)
      inv(x)=exp(x)

      fmsr=rnfpow(-1,mmin,mmax)
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine getfmsr(io,mmin,mmax,ir,mr)
c
cinput io   : class of resonance
cinput mmin : minimum mass
cinput mmax : maximum mass
coutput ir  : ID of resonance
coutput mr  : mass of resonance
c
c  {\tt getfmsr} is an extension of {\tt fmsr} to simulate
c  a resonance structure at low masses
c
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      implicit none
      integer io,ir,whichres
      real*8 mmin,mmax,mr,massit,widit
      include 'comres.f'
      real*8 fmsr,mmean

      mr=fmsr(mmin,mmax)
c...if strings come into play > modify whichres!
      ir=whichres(mr,io)
      mr=mmean(1,massit(ir),widit(ir),mmin,mmax)
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      SUBROUTINE STREXCT(IFL11,IFL12,IB1,AM1,IFL21,IFL22,IB2,AM2,ECM,
     *IOPT,ba1,ams1,ba2,ams2,ifl31,ifl32,ifl41,ifl42)
c
cinput ifl11  : ID of (di)quark of projectile hadron
cinput ifl12  : ID of remaining quark of projectile hadron
cinput ib1    : baryon number of projectile hadron
cinput am1    : mass of projectile hadron
cinput ifl21  : ID of (di)quark of target hadron
cinput ifl22  : ID of remaining quark of target hadron
cinput ib2    : baryon number of target hadron
cinput am2    : mass of target hadron
cinput ecm    : $\sqrt{s}$ of excitation
cinput iopt   : flag for excitation ansatz (1: Fritiof, 2:QGSM)
coutput ba1   : velocity of 1st string
coutput ams1  : mass of 1st string
coutput ba2   : velocity vector of 2nd string
coutput ams2  : mass of 2nd string
coutput ifl31 : (di)quark content of 1st string
coutput ifl32 : remainig quark content of 1st string
coutput ifl41 : (di)quark content of 2nd string
coutput ifl42 : remaining quark content of  2nd string
c
c output : particle IDs, masses and  momenta via common block {\tt newpart}
c
c This routine performs the excitation of two strings according
c to different ansatzes defined via parameter {\tt iopt}
c
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      include 'options.f'

      real*8 ba1(3),ba2(3),P(2,5),ranff
      integer IFLS(2,2)

      PARAMETER(PI=3.14159)! value of Pi
      PARAMETER(XMAX=1.0)
      LOGICAL SPINT,bt

      sigma=CTParam(31)
      alpha=CTParam(32)
      betav=CTParam(33)
      ampit=CTParam(34)
      betas=CTParam(40)
c
C   CHECK INITIAL ENERGY
      IF(AM1+AM2+ampit.gt.ECM) THEN
         WRITE(6,*) '....STOP! Initial energy is too low for string
     * excitation, ECM=', ECM,'<',am1,'+',am2
         stop 137
      ENDIF

      bt=.false.
      if(ranff(0).lt.5d-1)bt=.true.

C   COMPUTE XMIN
      XMIN=AMPIT/ECM

C   COMPUTE C.M.S. HADRON MOMENTUM AND ENERGIES
      PCMS=SQRT(ALAMB(ECM**2,AM1**2,AM2**2))/(2.0*ECM)
      ECMS1=SQRT(PCMS**2+AM1**2)
      ECMS2=SQRT(PCMS**2+AM2**2)

C   COMPUTE LIGHT-CONE VARIABLES FOR HADRON-PROJECTILE :
C                                           (PPLS01,PMINS01,PTX01,PTY01)
      PPLS01=ECMS1+PCMS
      PMINS01=ECMS1-PCMS
      PTX01=0.
      PTY01=0.

C   COMPUTE LIGHT-CONE VARIABLES FOR HADRON-TARGET :
C                                           (PPLS02,PMINS02,PTX02,PTY02)
      PPLS02=ECMS2-PCMS
      PMINS02=ECMS2+PCMS
      PTX02=0.
      PTY02=0.

C   COMPUTE TRANSFERRED TRANSVERSE MOMENTUM :
C                                           (QX,QY)
100   CALL GAUSPT(QT,SIGMA)
      PHI=2.*PI*ranff(0)
      QX=QT*COS(PHI)
      QY=QT*SIN(PHI)

C   COMPUTE XPLUS AND XMINUS VALUES FOR INTERACTING PARTONS
c      IF(IOPT.EQ.1) THEN
       if(IOPT.eq.1.or.dabs(alpha).lt.1d-4)then
         CALL XSDIS(XPLUS,XMIN,XMAX,BETAS)
         CALL XSDIS(XMINUS,XMIN,XMAX,BETAS)
       else ! IF(IOPT.EQ.2) THEN
         XPLUS=XVDIS(XMIN,ALPHA,BETAV)
         XMINUS=XVDIS(XMIN,ALPHA,BETAV)
c      else
c       write(6,*)'strexct: undefined excitation model'
c       stop 137
      ENDIF
C  COMPUTE LIGHT-CONE PARAMETERS OF INTERACTING PARTONS :
C                                      (PPRPLS1,PPRMINS1,PPRX1,PPRY1)
      PPRPLS1=XPLUS*PPLS01
      PPRMINS1=0.
      PPRX1=0.
      PPRY1=0.
C
C                                      (PPRPLS2,PPRMINS2,PPRX2,PPRY2)
      PPRPLS2=0.
      PPRMINS2=XMINUS*PMINS02
      PPRX2=0.
      PPRY2=0.
C
C  COMPUTE LIGHT-CONE COMPONENT OF TRANSFERRED MOMENTUM :
C                                                       (QPLS,QMINS)
      QPLS=-(QT**2/PPRMINS2)
      QMINS= QT**2/PPRPLS1
C
C  COMPUTE LIGHT-CONE PARAMETERS AND QUARK CONTENTS FOR EXCITED STRINGS
      IF(IOPT.EQ.1) THEN
C  FRITIOF ANSATZ :
C  P1 ---> P1 + Q
C  P2 ---> P2 - Q
         PPLS1=PPLS01+QPLS
         PMINS1=PMINS01+QMINS
         PTX1=PTX01+QX
         PTY1=PTY01+QY
C
         PPLS2=PPLS02-QPLS
         PMINS2=PMINS02-QMINS
         PTX2=PTX02-QX
         PTY2=PTY02-QY
C
C  NO QUARK REARRANGEMENT
         IFLS(1,1)=IFL11
         IFLS(1,2)=IFL12
         IFLS(2,1)=IFL21
         IFLS(2,2)=IFL22
      ELSEIF(IOPT.EQ.2) THEN
C  QGSM ANSATZ :
C  P1 ---> P1 - PPR1 + PPR2 - Q
C  P2 ---> P2 - PPR2 + PPR1 + Q
C
         PPLS1=PPLS01-PPRPLS1+PPRPLS2-QPLS
         PMINS1=PMINS01-PPRMINS1+PPRMINS2-QMINS
         PTX1=PTX01-PPRX1+PPRX2-QX
         PTY1=PTY01-PPRY1+PPRY2-QY
C
         PPLS2=PPLS02-PPRPLS2+PPRPLS1+QPLS
         PMINS2=PMINS02-PPRMINS2+PPRMINS1+QMINS
         PTX2=PTX02-PPRX2+PPRX1+QX
         PTY2=PTY02-PPRY2+PPRY1+QY

C  QUARK REARRANGEMENT IS NEEDED TO CREATE COLOUR NEUTRAL STRINGS
C  (In case of Baryon-Baryon or Meson-Meson or
C  Antibaryon-Antibaryon or Meson-Antibaryon or Antibaryon-Meson
C  interaction)
         IF((IB1.EQ.1.AND.IB2.EQ.1).OR.
     .        (IB1.EQ.-1.AND.IB2.EQ.-1).OR.
     .        (IB1.EQ.0.AND.IB2.EQ.0).OR.
     .        (IB1.EQ.0.AND.IB2.EQ.-1).OR.
     .        (IB1.EQ.-1.AND.IB2.EQ.0)) THEN
            IFLS(1,1)=IFL11
            IFLS(1,2)=IFL22
            IFLS(2,1)=IFL21
            IFLS(2,2)=IFL12
         ENDIF

C  QUARK REARRANGEMENT  (In case of Meson-Baryon or Antibaryon-Baryon or v.v.
C  interaction)
         IF((IB1.EQ.0.AND.IB2.EQ.1).OR.
     .        (IB1.EQ.1.AND.IB2.EQ.0).OR.
     .        (IB1.EQ.-1.AND.IB2.EQ.1).OR.
     .        (IB1.EQ.1.AND.IB2.EQ.-1)) THEN
            IFLS(1,1)=IFL11
            IFLS(1,2)=IFL21
            IFLS(2,1)=IFL22
            IFLS(2,2)=IFL12
         ENDIF
      ELSE
         PPLS1=0
         PMINS1=0
         PTX1=0
         PTY1=0
C
         PPLS2=0
         PMINS2=0
         PTX2=0
         PTY2=0
C
         WRITE(6,*) 'ERROR IN MAKE22: WRONG IOPT'
         stop 137
      ENDIF

C  COMPUTE OUTGOING STRING MASSES. They should be more than stable hadron
C  (with the same quark content) masses. In the contrary case generation
C  should be repeated
      AMS1S=PPLS1*PMINS1-QT**2
      SPINT=.TRUE.
      AMS0=0
C  COMPUTE AMS0
      IF(MOD(IFLS(1,1),100).EQ.0.AND.MOD(IFLS(1,2),100).EQ.0) THEN
C  qq-qqbar string
       IFLU=1  ! add u- quark to construct hadron
C  AMS0 is sum of masses of lowest baryon states
       IFLH1=IDPARS(IFLS(1,1), -ISIGN(IFLU,IFLS(1,1)),SPINT,2)
       IFLH2=IDPARS(IFLS(1,2), -ISIGN(IFLU,IFLS(1,2)),SPINT,2)
       AMS0=amass(IFLH1)+amass(IFLH2)
      ENDIF
      IF(.NOT.(MOD(IFLS(1,1),100).EQ.0
     .    .AND.MOD(IFLS(1,2),100).EQ.0))THEN
C  AMS0 is mass of lowest hadron state
        IKH=IDPARS(IFLS(1,1),IFLS(1,2),SPINT,2)
        AMS0=amass(IKH)
      ENDIF
C
      IF(AMS1S.LT.AMS0**2) GO TO 100
C
      AMS2S=PPLS2*PMINS2-QT**2
      SPINT=.TRUE.
C  COMPUTE AMS0
      IF(MOD(IFLS(2,1),100).EQ.0.AND.MOD(IFLS(2,2),100).EQ.0) THEN
C  qq-qqbar string
       IFLU=1  ! add u- quark to construct hadron
C  AMS0 is sum of masses of lowest baryon states
       IFLH1=IDPARS(IFLS(2,1),-ISIGN(IFLU,IFLS(2,1)),SPINT,2)
       IFLH2=IDPARS(IFLS(2,2),-ISIGN(IFLU,IFLS(2,2)),SPINT,2)
       AMS0=amass(IFLH1)+amass(IFLH2)
      ENDIF
      IF(.NOT.(MOD(IFLS(2,1),100).EQ.0
     .    .AND.MOD(IFLS(2,2),100).EQ.0)) THEN
C  AMS0 is mass of lowest hadron state
        IKH=IDPARS(IFLS(2,1),IFLS(2,2),SPINT,2)
        AMS0=amass(IKH)
      ENDIF
C
      IF(AMS2S.LT.AMS0**2) GO TO 100
C
      AMS1=SQRT(AMS1S)
      AMS2=SQRT(AMS2S)
C
C SUM OF MASSES OF EXCITED STRINGS SHOULD BE LESS THAN INITIAL ENERGY
      IF(AMS1+AMS2.GT.ECM) GO TO 100
C
      P(1,1)=PTX1
      P(1,2)=PTY1
      P(1,3)=0.5*(PPLS1-PMINS1)
      P(1,4)=0.5*(PPLS1+PMINS1)
      P(1,5)=AMS1
      P(2,1)=PTX2
      P(2,2)=PTY2
      P(2,3)=0.5*(PPLS2-PMINS2)
      P(2,4)=0.5*(PPLS2+PMINS2)
      P(2,5)=AMS2
c
      do 152 j=1,3
      ba1(j)=p(1,j)/p(1,4)
 152  ba2(j)=p(2,j)/p(2,4)

      ifl31=ifls(1,1)
      ifl32=ifls(1,2)
      ifl41=ifls(2,1)
      ifl42=ifls(2,2)

C
C  CHECK ENERGY-MOMENTUM CONSERVATION
      ESTR=P(1,4)+P(2,4)
      PSTRX=PTX1+PTX2
      PSTRY=PTY1+PTY2
      PSTRZ=P(1,3)+P(2,3)
C
c      if(abs(estr-ecm).gt.1d-12)then
c           WRITE(6,*) 'ECM=',ECM
c           WRITE(6,*) 'ESTR=',ESTR
c      end if
      if(pstrx.gt.1d-12)WRITE(6,*) 'PSTRX=',PSTRX
      if(pstry.gt.1d-12)WRITE(6,*) 'PSTRY=',PSTRY
      if(pstrz.gt.1d-12)WRITE(6,*) 'PSTRZ=',PSTRZ
C Note! One should also check flavour conservation!
C
      RETURN
      END
C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 FUNCTION ALAMB(X,Y,Z)
c
c input :  {\tt x,y,z}
c
C   THIS ROUTINE COMPUTES KINEMATICAL FUNCTION
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit none
      real*8 x,y,z

      ALAMB=max(0d0,(X-Y-Z)*(X-Y-Z) - 4.D0*Y*Z)

      RETURN
      END

C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 FUNCTION XVDIS(XMIN,ALFA,BETA)
c
cinput xmin : lower boundary
cinput alfa : parameter
cinput beta : parameter
c
C This function returns {\tt xmin}$ < x < 1$
c values distributed according to the Beta-function
C
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit none
      real*8 x,xmin,alfa,beta

 100  CALL SBETA(X,ALFA,BETA)
      IF(X.LE.XMIN) GO TO 100
      XVDIS=X
      RETURN
      END

C####C##1#########2#########3#########4#########5#########6#########7##
      SUBROUTINE SBETA(X,ALFA,BETA)
c
coutput x   : Beta-distributed value
cinput alfa : parameter
cinput beta : parameter

C THIS ROUTINE GENERATES X ACCORDING TO BETA DISTRIBUTION
C   $U(X)=C*X**(ALFA-1)*(1-X)**(BETA-1)$
C IONK,S METHOD IS USED
c
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit none
      real*8 ran1,ran2,alfa,beta,r1a,r2b,r12,x,ranff

 1    RAN1=ranff(0)
      RAN2=ranff(0)
      R1A=RAN1**(1./ALFA)
      R2B=RAN2**(1./BETA)
      R12=R1A+R2B
      IF(R12.GE.1.) GO TO 1
      X=R1A/R12
      RETURN
      END



C####C##1#########2#########3#########4#########5#########6#########7##
      SUBROUTINE XSDIS(X,XMIN,XMAX,BETA)
c
cinput xmin : lower boundary
cinput xmax : upper boundary
cinput beta : parameter
c
C  THIS FUNCTION GENERATES $XMIN < X < XMAX$ ACCORDING TO
C DISTRIBUTION $U(X)= 1./X*(1.-X)**(BETA+1)$ DISTRIBUTION
C PARAMETER $BETA > 0$
c
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      PE1=(1.-XMIN)**(BETA+1.)/(BETA+1.)-(1.-XMAX)**(BETA+1.)/
     *(BETA+1.)
      PE=PE1+DLOG(XMAX/XMIN)
      PE1=PE1/PE
 108  RND=ranff(0)
      RNDMPE1=ranff(0)
      IF(RNDMPE1.GT.PE1) GO TO 200
      X=1.-((1.-XMIN)**(BETA+1.)*(1.-RND)+(1.-XMAX)**(BETA+1.)*
     *RND)**(1./(BETA+1.))
      GO TO 300
200   X=XMIN*(XMAX/XMIN)**RND
300   PPE1=(1.-X)**BETA
      PPE2=1./X
      PPE1=PPE1+PPE2
      PPE2=PPE1*PPE2
      IF(PPE1*ranff(0).GT.PPE2) GO TO 108
      RETURN
      END


C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine getres(io,e,im,ip,ir)
c
cinput io : resonance class
cinput e  : $\sqrt{s}$ of excitation process
cinput im : lower boundary for resonance ID's
cinput ip : upper boundary for resonance ID's
coutput ir : ID of resonance
c
c This subroutine randomly selects a resonance ID for a given
c resonance class ({\tt io}), collision $\sqrt{s}$ and ID possible range.
c The selected resonance choice is weighted according to the respective
c partial cross section for the excitation of such a resonance.
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      include 'comres.f'
      include 'options.f'
      integer io,i,im,ip,ir
      real*8 e,x(minnuc:maxmes),xmax

      do i=im,ip
         call crossf(io,e,i,x(i))
      end do
      call getbran(x,minnuc,maxmes,xmax,im,ip,ir)

      if(ir.gt.ip.or.ir.lt.im) then
         write(6,*)'***(E) getres: no final state selected...'
         write(6,*)io,e,xmax,im,ip,ir
      endif

      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine aqm(i1,i2,xt,xel)
c
cinput :  {\tt i1,i2} : ID's of particle 1 and 2
coutput : {\tt xt,xel}: total and elastic cross section
c
c This subroutine returns cross sections according
c to the additive quark model.
c
C####C##1#########2#########3#########4#########5#########6#########7##
      integer i(2),ij,j,i1,i2,strit
      real*8 s(2),xt,xel,mn(2)
      integer ifa, ifb, isoit
      include 'comres.f'

      i(1)=i1
      i(2)=i2
      do 108 j=1,2
        ij=i(j)
        if(iabs(ij).ge.minmes)then
          s(j)=1d0*abs(strit(ij))
          call ityp2id(ij, isoit(ij), ifa, ifb)
c take care of hidden strangeness
          if (abs(ifa).eq.3 .and. abs(ifb).eq.3) s(j)=2
          mn(j)=1.d0
        else
          s(j)=1d0*abs(strit(ij))
          mn(j)=0.d0
        end if
 108  continue

      xt=max(0d0,40d0*0.666667**(mn(1)+mn(2))
     *  *( 1d0-0.4*s(1)/(3d0-mn(1)) )
     *  *( 1d0-0.4*s(2)/(3d0-mn(2)) ) )
      xel=0.039*xt**1.5d0
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      real* 8 function fppfit(iio,e,i3,i4)
c
c     Version: 1.0
c
cinput iio    : tag for cross section class
cinput e     : $\sqrt{s}$ of collision
cinput i3    : ID of first outgoing particle
cinput i4    : ID of second outgoing particle
c
c
c     {\tt fppfit} returns the isospin-independent part
c     of the production cross section
c     for one or two outgoing resonances ({\tt i3} and {\tt i4})
c     in a proton-proton
c     collision. {\tt io} sets the class of cross section
c     which is returned (i.e. $p p \rightarrow N \Delta$). If
c       {\tt io} is set 99 the class will be determined according
c       to {\tt i3} and {\tt i4}.
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
        include 'comres.f'
      include 'comwid.f'
      include 'options.f'

      real*8 e,tmp
      integer i,iio,io,im,jm,id1,id2,id3,id4,i3,i4
        integer i1m,i1p,i2m,i2p,class1,class2,iim
c funct.
      real*8 pmean,widit,massit,splint
      integer jit,isoit

c...parameters for resonance production in NN
      integer nfit
      parameter (nfit=7)
      real*8 ar(nfit)
      integer rr(2,nfit)
c io =         1      2     3     4       5      6      7
c              nd    nn*    nd*   dd      dn*    dd*    B*B*
      data ar/ 4d4,  6.3d0, 12d0, 2.8d0,   3.5d0, 3.5d0, 0.d0 /
      data rr/ 4,0,   4,1,   4,2,  0,0,   0,1,  0,2,    3,3 /
c rr tells of which particle class the out particles are
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c hard cut: no resonant cross section at highest energies
        if(e.gt.maxtab2)then
          fppfit=0d0
          return
        endif
c

      id1=0
      id2=0
      id3=0
      id4=0
      im=i3
      jm=i4

        if(iio.eq.99)then
c determine iline

c get class of 1 particle
          do 101 i=0,4
            call getrange(i,i1m,i1p)
            if(im.ge.i1m.and.im.le.i1p)then
              class1=i
              goto 102
            endif
101       continue
          write(*,*)'(E) fppfit: no x-section parametrized for ',im
          stop 137
102       continue
c get class of 2 particle
          do 103 i=0,4
            call getrange(i,i2m,i2p)
            if(jm.ge.i2m.and.jm.le.i2p)then
              class2=i
              goto 104
            endif
103       continue
          write(*,*)'(E) fppfit: no x-section parametrized for ',jm
          stop 137
104       continue
c get iline corrsponding to pp->i3,i4
          do 105 i=1,nfit
            if((rr(1,i).eq.class1.and.rr(2,i).eq.class2).or.
     v       (rr(1,i).eq.class2.and.rr(2,i).eq.class1))then
              io=i
              goto 106
            endif
105       continue
c maybe we have an R*R* reaction?
          if(im.ge.minres.and.im.le.maxres
     &     .and.jm.ge.minres.and.jm.le.maxres)then
             io=7
             goto 106
          endif
          write(*,*)'(E) fppfit: no iline found for particles',im,jm
          stop 137
106       continue
        else
          io=iio
        endif

c sort particles acc. to itypes
        if(io.ne.5.and.im.gt.jm)then
          im=i4
          jm=i3
        endif
        if(io.eq.5.and.im.lt.jm)then
          im=i4
          jm=i3
        endif

c consistency checks
        call getrange(rr(1,io),i1m,i1p)
        if(im.lt.i1m.or.im.gt.i1p)then
          write(*,*)'(E) fppfit: wrong iline for outgoing particles',
     &      io,im,jm
          stop 137
        endif
        call getrange(rr(2,io),i2m,i2p)
        if(jm.lt.i2m.or.jm.gt.i2p)then
          write(*,*)'(E) fppfit: wrong iline for outgoing particles',
     &      io,im,jm
          stop 137
        endif


      goto(1,2,2,2,2,2,4) io
      write(6,*) '****(E) wrong x-section ID in fppfit *****',io
      stop 137

c pp ->ND
 1    continue
        if(wtabflg.ge.3.and.CTOption(9).eq.0)then
c table lookup
           if(im.eq.minnuc)then
             iim=1
           elseif(im.eq.mindel)then
             iim=2
           else
             write(*,*)'(E) fppfit: First particle should be N or Delta'
             stop 137
           endif
           fppfit=max(0d0,splint(tabxnd,frrtaby(1,1,iim,jm),
     &       frrtaby(1,2,iim,jm),widnsp,e))
      else
c calculate cross section
         fppfit=pmean(e,im,isoit(im),jm,isoit(jm),id1,id2,id3,id4,1)
     &       /(pmean(e,1,1,1,1,id1,id2,id3,id4,1)
     &       *e*e)*ar(io)
     &       *(massit(jm)**2*widit(jm)**2/((e**2-massit(jm)**2)**2
     &       +widit(jm)**2*massit(jm)**2))
     &       *dble((jit(im)+1)*(jit(jm)+1))
        endif

      return

c pp->NN* pp->ND* pp->DD pp->DN* pp->DD*
 2      continue
        if(wtabflg.ge.3.and.CTOption(9).eq.0)then
c table lookup
           if(im.eq.minnuc)then
             iim=1
           elseif(im.eq.mindel)then
             iim=2
           else
             write(*,*)'(E) fppfit: First particle should be N or Delta'
             stop 137
           endif
           fppfit=max(0d0,splint(tabxnd,frrtaby(1,1,iim,jm),
     .           frrtaby(1,2,iim,jm),widnsp,e))
      else
c calculate cross section
         tmp=pmean(e,im,isoit(im),jm,isoit(jm),id1,id2,id3,id4,1)
     &         /(pmean(e,1,1,1,1,id1,id2,id3,id4,1)
     &         *e*e)*ar(io)
     &         *dble((jit(im)+1)*(jit(jm)+1))
         if(im.ne.jm) then
            tmp=tmp
     &         /((massit(jm)-massit(im))**2*(massit(jm)+massit(im))**2)
         endif
         fppfit=tmp
      endif
      return

c pp->B*B*
 4    continue
c sofar set to zero
      fppfit=0d0
      return

      end

C**********************************************************************
C     Program calculates cross section for reaction m1+m2 <--> m3+m4
C     including spin factors and using detailed balance, [mb]
C     IFL=1 for direct reaction
C     IFL=2 for backward reaction
C     IS=0: no strange quark;  IS=1: single strange quark
C     -----------------------------------------------------------------
      real*8 FUNCTION DMCROSS(SRT,DM1,DM2,DM3,DM4,S1,S2,S3,S4,IFL,IS)
c        INCLUDE'COMMON'
        REAL*8 SRT0(10),MPI,MJP,MRO,MK,MKS,MRHO,MCHI,MPP
        REAL*8 SRT,DM1,DM2,DM3,DM4
        INTEGER S1,S2,S3,S4,IFL,IS

c        write(*,*)SRT,DM1,DM2,DM3,DM4,S1,S2,S3,S4,IFL,IS

C        WRITE(80,*)'DMCROSS started',SRT
        F1=0.0
C----------------------------
C 18.02.03 Enhanced matrix element for chi_c and Psi'
C set 9:
        FMA=0.65
        IF(DM1.GT.3.2 .or. DM2.GT.3.2) FMA=0.65
C-----------------------------
        DMCROSS=0.0
c        PI=4.*ATAN(1.0)
*       Matrix elements^2 |M|^2 in mb/GeV^2
c        FAC=16.*0.13*FMA ! used for comover scenario, "old"

        FAC=16.* 0.13*FMA /1.5 ! used for comover scenario - NEW for B abs. 4.18 mb
c	FAC=0.  ! no meson absorption!
C------------------------------
        IF(IS.EQ.1) FAC=FAC/3.
        IF(IS.EQ.2) FAC=FAC/9.
*
        SOUT=DM3+DM4
        SIN=DM1+DM2
*
        IF(SRT.LE.SOUT)               RETURN
        IF(SRT.LE.SIN)                RETURN
*
        S=SRT**2
        P2OUT=(S-(DM3+DM4)**2)*(S-(DM3-DM4)**2)/(4.*S)
        P2IN= (S-(DM1+DM2)**2)*(S-(DM1-DM2)**2)/(4.*S)
        E1=SQRT(P2IN+DM1**2)
        E2=SQRT(P2IN+DM2**2)
        E3=SQRT(P2OUT+DM3**2)
        E4=SQRT(P2OUT+DM4**2)
        EE=E1*E2*E3*E4
*
        Y1=SOUT/SRT
*
        IF(Y1.GE.1.0)                        RETURN
        SPINOUT=(2.*S3+1.)*(2.*S4+1.)
        IF(SPINOUT.EQ.9) FAC=FAC/3.
        SPININ=(2.*S1+1.)*(2.*S2+1.)
*       constant Matrix element x Y1**6
        F1=EE*SQRT(P2OUT)/SQRT(P2IN)*Y1**6*FAC/S
*
        F1=F1*SPINOUT
*
        IF(IFL.EQ.1)   THEN
        DMCROSS=F1 ! cross section in mb
        RETURN
        ENDIF
*-------------------------------------------------------------------------
*       INVERSE REACTION for IFL=2
*
*
        IF(IFL.EQ.2.AND.P2OUT.GT.0.0)    THEN
        DMCROSS=F1*P2IN/P2OUT*SPININ/SPINOUT
        ENDIF
C       WRITE(80,*)'DMCROSS ended',DMCROSS
*-------------------------------------------------------------------------
        RETURN
        END

C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine crossx(iio,e,ii1,iiz1,mm1,ii2,iiz2,mm2,sig)
c
cinput iio  : pointer to cross section
cinput e    : $\sqrt{s}$ of collision
cinput ii1  : ID of particle 1
cinput iiz1 : $2\cdot I_3$ of particle 1
cinput mm1  : mass of particle 1
cinput ii2  : ID of particle 2
cinput iiz2 : $2\cdot I_3$ of particle 2
cinput mm2  : mass of particle 2
coutput sig : cross section
c
c This routine returns cross sections which are accessed via tags (pointers)
c {\tt iio}. The cross sections can either be partial ones for specific
c exit channels or total cross sections for the two incoming hadrons
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none

      include 'comres.f'
      include 'options.f'
      include 'newpart.f'
      include 'coms.f'

      integer io,iio,i1,i2,iz1,iz2
      integer i,im,ip,j,jm,jp,ii1,iiz1,ii2,iiz2,icharm,ilight
c      integer itmp1,itmp2,itag
      real*8 e,e0,sum,m1,m2,m3,gam,dum,s1,s2,s3
      real*8 DMCROSS,sig,sig1,sig2,sig3
      real*8 t,s,a,b,c,d,f2,g,a3,b3,c3,d3,g3
      real*8 c1,c2,mm1,mm2,aaqm,dbfact,cgkcor,qflg,factor,ggam
      real*8 sighelp,meltpoint2,pfrome,pmelt,pmelt2,plab
      real*8 ma1,ma2,ma3,ma4,pf,pin,cfac
      integer n1,n2,n3,n4,ifl,is
      logical kplflag
c...functions
      real*8 clebsch,siglookup,bcms,massit,widit,mminit,xmelt
      real*8 fppfit,fwidth,fbwnorm,ppiso,pcms
      real*8 STREXWrapper
c...additional functions and variables for ppbar&pp scattering:
      real*8 bbphi,nnphi,sighera,sapptot,
     &       sappela,sappann,sappdiff,dgcgkfct
      integer isoit,iq1(2),iq2(3),ifdiq2,kq1,kq2

c
      real*8 sigheramb,sigmb,meltpoint,sigma1,sigma2,sigma3,sigma4

      logical excl

      real*8 ef,sigf,ee
      integer jf,if

c...parameters for resonance production in NN
      integer nfit
      parameter (nfit=7)
c      real*8 er(nfit)
      integer rr(2,nfit),in(nfit)
      real*8 x,x0,xa,xb,xc,bf
      real*8 f

c     JSFC -
      real*8  decm,aecm
      integer iecm,izcharm,izlight
      integer extra
      real*8 slope
c     JSFC


c...parameters for probability of processes 53,57 & 58
      real*8 ranff


      save kplflag
      f(x,x0,xa,xb,xc)=max(0d0,2.*xb*xa*(x-x0)/(xb**2+(x-x0)**2))*
     *     ((x0+xb)/(x))**xc


c io =         1      2     3     4       5      6      7
c              nd    nn*    nd*   dd      dn*    dd*    B*B*
c      data er/ 4d-2, 0d0,   0d0,  0d0,    0d0,   0d0,   0d0/
      data rr/ 4,0,   4,1,   4,2,  0,0, 0,1,  0,2,    3,3 /
c rr tells of which particle class the out particles are
      data in/ 1,     1,     1,    1,  -1,    1,      1  /
c in tells which out particle has the higher itype


      excl=.false. !flag for exclusive cross section(detailed balance)

      kplflag=.false. ! flag for k+ p total x-section

      goto 107

      entry crossf(if,ef,jf,sigf)
      io=if
      ee=ef
      j=jf
      call getrange(rr(1,io),i,ip) ! call getrange(rr(2,io),jm,jp)

      sigf=fppfit(io,ee,i,j)

      return

      entry crossz(iio,e,ii1,iiz1,mm1,ii2,iiz2,mm2,sig)

      excl=.true.  !flag for exclusive cross section(detailed balance)

 107  continue

c...some settings for all channels enter here
      sig=0d0
      io=iabs(iio)
      aaqm=CTParam(3)      !exponent for bcm scaled AQM-cross sections

      call setizm(ii1,iiz1,mm1,ii2,iiz2,mm2,
     ,i1,iz1,m1,i2,iz2,m2)

c in case of a MB-reaction, sort particles: meson must be second (i2)
      if(iabs(i1).ge.minmes.and.iabs(i1).le.maxmes.and.
     &   iabs(i2).ge.minbar.and.iabs(i2).le.maxbar)then
        call swpizm(i1,iz1,m1,i2,iz2,m2)
      endif

      goto(1,1,1,1,1,1,1,8,9,10,
     ,       11,12,13,14,15,16,17,18,19,9,21,22,23,24,
     ,       25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,
     ,       40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,
     ,       55,56,57,58,59,60,61,62,63,64)io


      write(6,*)'cross[x,z]: ',
     ,' unknown channel requested, wrong io:',io
      stop 137


 1    continue
c...pp-> ND pp->NN* pp->ND* pp->DD pp->DN* pp->DD*
      call setizm(iabs(ii1),iiz1,mm1,iabs(ii2),iiz2,mm2,
     ,i1,iz1,m1,i2,iz2,m2)
c     fist outgoing particle has lower ID
      if(i1.gt.i2)call swpizm(i1,iz1,m1,i2,iz2,m2)

c     determine range for first outgoing particle
      call getrange(rr(1,io),im,ip)


      if(excl)then
cut down the loop to one resonance for inverse exclusive processes
        jm=in(io)*max0(in(io)*i1,in(io)*i2)
        jp=jm
      else
        call getrange(rr(2,io),jm,jp)

      end if
      sum=0.

      do i=im,ip
c...loop over i3 (1st outgoing particle)
        do j=jm,jp
c...loop over i4 (2nd outgoing particle)
c   sum over x-section

c explicit isospin dependence
           cgkcor=ppiso(iio,i1,iz1,i2,iz2,i,j)

           sum=sum+fppfit(io,e,i,j)*cgkcor
        end do
      end do
      sig=sum
c     JS set hard cut for mass threshhold
      if(e.lt.2.015d0)sig=0.0d0

c
c     detailed balance
      if(iio.lt.0.and.sig.gt.1.d-12) then


            call detbal(e,i1,i2,iz1,iz2,max(mminit(i1),m1),
     &           max(mminit(i2),m2),1,1,dbfact)
            sig=sig*dbfact
      endif

c     check energy conservation (cut-off is one 1 MeV)
      if(max(m1,mminit(i1))+max(m2,mminit(i2))+1.d-3.gt.e) sig=0.d0

      return


 8    continue
c...ND->DD  iso-spin summed value
      call setizm(iabs(ii1),iiz1,mm1,iabs(ii2),iiz2,mm2,
     ,i1,iz1,m1,i2,iz2,m2)

      e0=2.*massit(mindel)-widit(mindel)

      c1=clebsch(isoit(i1),isoit(i2),iz1,iz2,1)
      c2=clebsch(isoit(i1),isoit(i2),iz1,iz2,2)

      if(io.lt.0)write(6,*)'crossz(DD->ND):c1,c2=',c1,c2,
     ,           isoit(i1),isoit(i2),iz1,iz2,'itypes:',i1,i2
      sig=f(e,e0,12.0d0,0.02d0,2d0)*(0.66667*c1+4d0/dsqrt(20d0)*c2)

      if(iio.lt.0.and.sig.gt.1.d-12) then
            call detbal(e,i1,i2,iz1,iz2,max(mminit(i1),m1),
     &           max(mminit(i2),m2),minnuc,mindel,dbfact)
            sig=sig*dbfact
      endif
      return

 9    continue
      write(6,*)'crossx: channel no.',io,'not implemented.'
      stop 137

 10   continue
c...MB->B'
      if(e.gt.mresmax)then
         sig=0.0d0
         return
      endif
      call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sig,gam)
      return

 11   continue
c...MM->M'
      call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sig,gam)
      return

 12   continue
c...??->X => additive quark model total cross section
      call aqm(i1,i2,sig,dum)
      if(sig.gt.1d3)then
        write(6,*)'sig=',sig
      end if

c     check energy conservation (cut-off is one 1 MeV)
      if(max(m1,mminit(i1))+max(m2,mminit(i2))+1.d-3.gt.e) sig=0.d0

c     Lambda+Lambda Strangness-Exchange Cross-section
      if(CTOption(59).ne.0) then
        if((abs(i1).eq.27.and.abs(i2).eq.27.and.i1*i2.gt.0).or.
c       Lambda+Sigma Strangness-Exchange Cross-section
     &     (((abs(i1).eq.27.and.abs(i2).eq.40).or.
     &       (abs(i1).eq.40.and.abs(i2).eq.27)).and.
     &             i1*i2.gt.0).or.
c       Sigma+Sigma Strangness-Exchange Cross-section
     &       (abs(i1).eq.40.and.abs(i2).eq.40.and.
     &           i1*i2.gt.0)) then
          sig=sig+STREXWrapper(63,e,i1,iz1,i2,iz2)
c       Xi + Nucleon Strangness-Exchange Cross-section
        elseif(((abs(i1).eq.49.and.abs(i2).eq.1).or.
     &          (abs(i1).eq.1.and.abs(i2).eq.49)).and.
     &             i1*i2.gt.0) then
          sig=sig+STREXWrapper(64,e,i1,iz1,i2,iz2)
        endif
      endif



      return

 13   continue
c...??->X => additive quark model elastic cross section
      call aqm(i1,i2,dum,sig)

c     check energy conservation (cut-off is one 1 MeV)
      if(max(m1,mminit(i1))+max(m2,mminit(i2))+1.d-3.gt.e) sig=0.d0

      return

 14   continue
c...??->X => additive quark model inelastic cross section
      call aqm(i1,i2,s1,s2)
      s3=s1-s2
      s1=s3
c DD->DD is possible only for e<4gev, then DD->Strings (15) is used
      e0=mminit(i1)+mminit(i2)+2d0*CTParam(4)
      bf=2d0*dabs(CTParam(4)-CTParam(2))
      s2=f(e,e0,s3,bf,1d0)
      sig=xmelt(e,s1,s2,max(mminit(i1),m1)+max(mminit(i2),m2),
     @     e0+2d0*CTparam(2))

c     check energy conservation (cut-off is one 1 MeV)
      if(max(m1,mminit(i1))+max(m2,mminit(i2))+1.d-3.gt.e) sig=0.d0

      return

 15   continue
c...??->X => additive quark model string cross section
      if(CTOption(12).ne.0)return
      if(CTOption(12).eq.0)then
        call aqm(i1,i2,sig,dum)
        sig=sig-dum
c dum is the inelastic xsec
c string threshold is at 3.2 GeV by default
        if(e.gt.2d0*(1.08d0+ctparam(2)))then
          sig=sig*bcms(e,1.08d0+ctparam(2),1.08d0+ctparam(2))**(aaqm)
        else
          sig=0.d0
        end if
      end if

c     check energy conservation (cut-off is one 1 MeV)
      if(max(m1,mminit(i1))+max(m2,mminit(i2))+1.d-3.gt.e) sig=0.d0

      return

 16   continue
c...pn-total...low energy tables & high energy fit
      sig=xmelt(e,siglookup(1,e),sighera(.938d0,.938d0,e,3),3.d0,5.d0)
      return

 17   continue
c...pn-elastic...low energy tables & high energy fit (then pn=pp)
      sig=xmelt(e,siglookup(2,e),sighera(.938d0,.938d0,e,2),3.d0,5.d0)
      return

 18   continue
c...pp-total...low energy tables & high energy fit
      sig=xmelt(e,siglookup(3,e),sighera(.938d0,.938d0,e,1),3.d0,5.d0)
      return

 19   continue
c...pp-elastic...low energy tables & high energy fit
      sig=xmelt(e,siglookup(4,e),sighera(.938d0,.938d0,e,2),3.d0,5.d0)
      return

 20   continue
      write(6,*)'crossx: cross section for decay(io=20) requested'
      stop 137
      return

 21   continue
c...bbar total...from ppbar via AQM & phase-space
      call aqm(1,-1,nnphi,dum)
      call aqm(i1,i2,bbphi,dum)
      sig=sapptot(e,m1,m2)*bbphi/nnphi
      return


 22   continue
c...bbar elastic...from ppbar via AQM & phase-space
      call aqm(1,-1,nnphi,dum)
      call aqm(i1,i2,bbphi,dum)
      sig=sappela(e,m1,m2)*bbphi/nnphi
      return

 23   continue
c...bbar annihilation from ppbar via AQM & phase-space
cbb   No cross-section if strings are not allowed.
      if(CTOption(12).ne.0) then
        sig = 0.0
        return
      endif
      call aqm(1,-1,nnphi,dum)
      call aqm(i1,i2,bbphi,dum)
      sig=sappann(e,m1,m2)*bbphi/nnphi
      return

 24   continue
c...bbar string...from ppbar via AQM & phase-space
cbb   No cross-section if strings are not allowed.
      if(CTOption(12).ne.0) then
        sig = 0.0
        return
      endif
      call aqm(1,-1,nnphi,dum)
      call aqm(i1,i2,bbphi,dum)
      sig=xmelt(e,0d0,sappdiff(e,m1,m2)*bbphi/nnphi,3d0,5d0)
      return

 39   continue
c...bbar inelastic...from ppbar via AQM & phase-space
cbb   No cross-section if strings are not allowed.
      if(CTOption(12).ne.0) then
        sig = 0.0
        return
      endif
      call aqm(1,-1,nnphi,dum)
      call aqm(i1,i2,bbphi,dum)
      sig=xmelt(e,sappdiff(e,m1,m2)*bbphi/nnphi,0d0,3d0,5d0)
      return

 25   continue
c...pi^+ +p
      if ((i2.eq.minmes+1).and.(iz2.eq.2)
     &   .and.(i1.eq.1).and.(iz1.eq.1)) then
         sigheramb=sighera(m2,m1,e,7)
         meltpoint=2.00d0
cJS...  pi^- +n
      elseif ((i2.eq.minmes+1).and.(iz2.eq.-2)
     &     .and.(i1.eq.1).and.(iz1.eq.-1)) then
         sigheramb=sighera(m2,m1,e,7)
         meltpoint=2.00d0
c...pi^- +p
      elseif ((i2.eq.minmes+1).and.(iz2.eq.-2)
     &   .and.(i1.eq.1).and.(iz1.eq.1)) then
         sigheramb=sighera(m2,m1,e,9)
         meltpoint=2.18d0
c...K^+ +p
      elseif ((i2.eq.minmes+6).and.(iz2.eq.1)
     &   .and.(i1.eq.1).and.(iz1.eq.1)) then
         sigheramb=sighera(m2,m1,e ,11)
         meltpoint=1.84d0
C        here second meltpoint to enable linear interpolation between
C        meltpoint2 and meltpoint (in plab)
         meltpoint2=1.7d0
         kplflag=.true.
c...K^- +p
      elseif ((i2.eq.-(minmes+6)).and.(iz2.eq.-1)
     &   .and.(i1.eq.1).and.(iz1.eq.1)) then
         sigheramb=sighera(m2,m1,e,14)
         meltpoint=2.12d0
c...K^+ +n
      elseif ((i2.eq.minmes+6).and.(iz2.eq.1)
     &   .and.(i1.eq.1).and.(iz1.eq.-1)) then
         sigheramb=sighera(m2,m1,e,13)
         meltpoint=1.75d0
c...K^- +n
      elseif ((i2.eq.-(minmes+6)).and.(iz2.eq.-1)
     &   .and.(i1.eq.1).and.(iz1.eq.-1)) then
         sigheramb=sighera(m2,m1,e,16)
         meltpoint=1.6d0
c...gamma + p
      elseif ((i2.eq.minmes).and.(iz2.eq.0)
     &   .and.iabs(i1).le.maxbar) then
         sigheramb=sighera(m2,m1,e,6)
         meltpoint=1.75d0
      else
c...MB total (->B*/->Strings/el.)
        call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sigmb,gam)
        if (sigmb.gt.1d-8) then
c meltpoint for resonant meson absorption moved to higher energies
csab: set to 2.18 for bbar-meson scattering
          meltpoint=max(2.18d0, m1+m2+.2d0)
        else
          meltpoint=2.18d0
        endif
        call aqm(pimeson,nucleon,dum,nnphi)
        call aqm(i1,i2,dum,bbphi)
        sigheramb=sighera(m2,m1,e,7)*bbphi/nnphi
C
      endif
C
C
c...breit-wigners...
      call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sigmb,gam)
c now melt the low energy resonance x-sec and the high energy hera-fits:
      if (e.lt.meltpoint) then
C
c cross section for Danielewicz forward delay
c     here the DP-cross section has the same form as the normal one
C
         if(CTOption(34).eq.2) then
            sig=sigmb+CTParam(58)*sigmb
         elseif(CTOption(34).eq.3) then
            m3=e
            call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sigmb,gam)
            if(sig.gt.1d-5) sig=sigmb+CTParam(58)
         elseif(CTOption(34).eq.4) then
            m3=e
            call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sigmb,gam)
            if(sigmb.gt.1d-5)then
               ggam=fwidth(itypnew(1),i3new(1),m3)
               sig=sigmb+CTParam(58)*
     &              sig*dsqrt(2.d0/(ggam*3.1415d0*fbwnorm(m3,
     &              itypnew(1),i3new(1))))
            endif
         else
C
C actually default for e.lt.meltpoint
C    take result from anndec as sigma
C
            if((abs(i1).eq.1).and.
     $           (abs(i2).eq.101.or.abs(i2).eq.106).or.
     $           (i1.eq.minmes).or.(i2.eq.minmes))then

               sig=sigmb
            else
               call aqm(pimeson,nucleon,dum,nnphi)
               call aqm(i1,i2,dum,bbphi)
               sig=sigmb+sighera(m2,m1,e,8)*bbphi/nnphi
            endif
         endif
C
c-- nonres adds s-channel strings for strange meson baryon reactions
c-- resonances are not sufficient to fit data - sig gets modified in nonres !!!
C-- so that k- p cross section describes data reasonable.
C-- nonres call only below meltpoint
C
         if((abs(i1).eq.1).or.(abs(i2).eq.1))then
            call nonres(e,i1,iz1,i2,iz2,sig)
         endif

c Add Strangeness exchange cross-section to total cross-section.
c For Kaon+N reactions this is already done in nonres.
         if((abs(i1).eq.40.or.abs(i1).eq.27).and.
     &       abs(i2).eq.101) then
           sig=sig+STREXWrapper(58,e,i1,iz1,i2,iz2)
         elseif((abs(i1).eq.40.or.abs(i1).eq.27).and.
     &           abs(i2).eq.106) then
           sig=sig+STREXWrapper(59,e,i1,iz1,i2,iz2)
         elseif(abs(i1).eq.49.and.abs(i2).eq.101) then
           sig=sig+STREXWrapper(60,e,i1,iz1,i2,iz2)
         elseif(abs(i1).eq.49.and.abs(i2).eq.106) then
           sig=sig+STREXWrapper(61,e,i1,iz1,i2,iz2)
         elseif(abs(i1).eq.55.and.abs(i2).eq.101) then
           sig=sig+STREXWrapper(62,e,i1,iz1,i2,iz2)
         end if
C
      else
C
C     for e.ge.meltpoint
C
         sig=sigheramb
C
      endif
C


c.. is annihilation possible ? (quark content?)
      call ityp2id(i2,iz2,iq1(1),iq1(2))
      call ityp2id(i1,iz1,ifdiq2,iq2(3))
      if(abs(i1).ge.minmes) then
        iq2(1)=ifdiq2
        iq2(2)=iq2(3)
        iq2(3)=0
      else
        iq2(1)=mod(ifdiq2/100,10)
        iq2(2)=int(ifdiq2/1000)
      endif
      qflg=0.d0
      do 7312 kq1=1,2
       do 7412 kq2=1,3
c.. annihilation ?
        if(iq1(kq1)+iq2(kq2).eq.0) qflg=1.d0
 7412   continue
 7312  continue
c..--> elastic X-section gets minimum val. of 12.5 mb (K+ P)
c..    if annihilation is not possible
c..      (exclude gamma baryon elastic)
      if  ((i2.eq.106).and.(i1.eq.1)) then
          sig=max(sig,12.5d0*(1d0-qflg))

C       for k+ p linear interpolation between meltpoint2 and meltpoint
          if (kplflag.and.
     &        qflg.eq.0.d0.and.e.lt.meltpoint
     &       .and.e.gt.meltpoint2) then
             sighelp=sighera(m2,m1,meltpoint,11)
             pfrome=plab(m2,m1,e)
             pmelt2=plab(m2,m1,meltpoint2)
             pmelt=plab(m2,m1,meltpoint)
             sig=xmelt(pfrome,1.25d1,sighelp,pmelt2,pmelt)
          endif
      endif

c     check energy conservation (cut-off is one 1 MeV)
      if(max(m1,mminit(i1))+max(m2,mminit(i2))+1.d-3.gt.e) sig=0.d0


C
C  end of MB total
C
      return

 26   continue
      if(i2.eq.minmes)then
c..no elastic gamma B scattering:
       sig=0.d0
       return
      endif
c...MB->MB elastic (pi+ p scaled with aqm)
      if(mminit(i1)+mminit(i2).gt.e) then
c     no collision if sqrts is insufficient to put particles on-shell
         sig=0.d0
      else
         call aqm(pimeson,nucleon,dum,nnphi)
         call aqm(i1,i2,dum,bbphi)
         sig=sighera(m2,m1,e,8)*bbphi/nnphi
      endif
      return

 27   continue
c..MB-> 1 String (s-channel) (pi+ p scaled with aqm)
      if(CTOption(12).ne.0)return
      call aqm(pimeson,nucleon,nnphi,dum)
      call aqm(i1,i2,bbphi,dum)
      sig=sighera(m2,m1,e,7)*bbphi/nnphi
c-- minus elastic x-sect. (scaled with el. aqm)
      call aqm(pimeson,nucleon,dum,nnphi)
      call aqm(i1,i2,dum,bbphi)
      sig=sig-sighera(m2,m1,e,8)*bbphi/nnphi
c.. is s-channel possible ? (quark content?)
      call ityp2id(i2,iz2,iq1(1),iq1(2))
      call ityp2id(i1,iz1,ifdiq2,iq2(3))
      if(abs(i1).ge.minmes) then
        iq2(1)=ifdiq2
        iq2(2)=iq2(3)
        iq2(3)=0
      else
        iq2(1)=mod(ifdiq2/100,10)
        iq2(2)=int(ifdiq2/1000)
      endif
      qflg=0.d0
      do 312 kq1=1,2
       do 412 kq2=1,3
c..two of the quarks must be able to annihilate:
        if(iq1(kq1)+iq2(kq2).eq.0) qflg=1.d0
 412   continue
 312  continue
      sig=xmelt(e,sig,0.d0,3d0,6.d0)*qflg
      if(e.le.1.7d0)sig=0.d0
c      call nonres(e,i1,iz1,i2,iz2,sig)

      return

 28   continue
c..MB-> 2 Strings (t-channel) (pi+ p scaled with aqm)
      if(CTOption(12).ne.0)return
      call aqm(pimeson,nucleon,nnphi,dum)
      call aqm(i1,i2,bbphi,dum)
      sig=sighera(m2,m1,e,7)*bbphi/nnphi
c-- minus elastic x-sect. (scaled with el. aqm)
      call aqm(pimeson,nucleon,dum,nnphi)
      call aqm(i1,i2,dum,bbphi)
      sig=sig-sighera(m2,m1,e,8)*bbphi/nnphi
c-- melting
      sig=xmelt(e,0d0,sig,3d0,6.d0)
      return

 29   continue
c...??->X => additive quark model inelastic cross section
c specially reduced to avoid double-counting in case of DN and DD
c    cross sections
      call aqm(i1,i2,sig,dum)
      sig=sig-dum
      sig=sig*bcms(e,max(mminit(i1),m1),max(mminit(i2),m2))**aaqm
      sig=xmelt(e,sig,0d0,2d0*(1.08d0+ctparam(2)),5.0d0)
c here comes the reduction
      sig=xmelt(e,0d0,sig,2.75d0,3.4d0)
      return

 30   continue
c parameterized detailed balance cross sections for ND->NN
      factor=dgcgkfct(i1,i2,iz1,iz2,nucleon,nucleon)
      if(factor.le.1.d-8) then
         sig=0.d0
         return
      endif
      cgkcor=ppiso(-1,i1,iz1,i2,iz2,nucleon,nucleon)

         sig=factor*cgkcor*(1.3d8*(e**(-17.5d0))+3.6d4*(e**(-7d0)))

c         write(6,*)factor,cgkcor,sig

      return

 31   continue
c parameterized detailed balance cross sections for DD->DN
      e0=massit(i1)-0.5*widit(i1)+massit(i2)-0.5*widit(i2)

      factor=dgcgkfct(i1,i2,iz1,iz2,nucleon,mindel)
      if(factor.le.1.d-8.or.e.lt.e0)then
         sig=0.d0
      else
         c1=clebsch(isoit(i1),isoit(i2),iz1,iz2,1)
         c2=clebsch(isoit(i1),isoit(i2),iz1,iz2,2)

c param
         sig=
     &        (2.5d56*exp(-(50.0d0*e))
     &        +4.9d14*exp(-(12.d0*e))
     &        +1.1d6*exp(-(4.50d0*e)))
     &        *factor
c the following factors are from Heinz Sorge's Habilitation
     &        *(0.66667*c1+4d0/dsqrt(20d0)*c2)

      endif
      return

 32   continue
C parameterized detailed balance cross sections for DD->NN
      factor=dgcgkfct(i1,i2,iz1,iz2,nucleon,nucleon)
      if(factor.le.1.d-8.or.e.lt.2.15d0) then
         sig=0.d0
         return
      endif
      cgkcor=ppiso(-4,i1,iz1,i2,iz2,nucleon,nucleon)
c param
      sig=(7.27d0*(e-2.14d0)**(-1.2176d0)+0.05d0*(e-2.14d0)**(-3.257d0))
     &     *factor*cgkcor

      return

 33   continue
c...??->X => additive quark model resonance cross section
       call aqm(i1,i2,sig,dum)
       sig=sig-dum  ! dum is the elastic xsec
      if(e.gt.mminit(i1)+mminit(i2))then
c x_string+x_resonances
          sig=sig*bcms(e,mminit(i1),mminit(i2))**(aaqm*3)
c get x_string
          sig=xmelt(e,sig,0d0,mminit(i1)+mminit(i2)+1.5d0,
     @         mminit(i1)+mminit(i2)+3d0)

      end if

      return

 34   continue
c...??->X => additive quark model string cross section
      if(CTOption(12).ne.0)return
        call aqm(i1,i2,sig,dum)
        sig=sig-dum !  dum is the elastic xsec
       if(e.gt.mminit(i1)+mminit(i2))then
c x_string+x_resonances
          sig=sig*bcms(e,mminit(i1),mminit(i2))**(aaqm*3)
c get x_resonances
          sig=xmelt(e,0d0,sig,mminit(i1)+mminit(i2)+1.5d0,
     @         mminit(i1)+mminit(i2)+3d0)
      end if
      return

 35   continue
      e0=mminit(i1)+mminit(i2)+CTParam(4)*2d0
      call aqm(i1,i2,sig,dum)
      sig=f(e,e0,sig-dum,1d0,1d0)
      return

 36   continue
c cross section for Danielewicz forward delay
c...MB->B'
      if(CTOption(34).eq.2) then
         call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sig,gam)
         sig=CTParam(58)*sig
      elseif(CTOption(34).eq.3) then
         m3=e
         call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sig,gam)
         if(sig.gt.1d-5) sig=CTParam(58)
      elseif(CTOption(34).eq.4) then
         m3=e
         call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sig,gam)
         if(sig.gt.1d-5)then
            ggam=fwidth(itypnew(1),i3new(1),m3)
            sig=CTParam(58)*
     &           sig*dsqrt(2.d0/(ggam*3.1415d0*
     &           fbwnorm(m3,itypnew(1),i3new(1))))
         endif
      else
         sig=0.d0
      endif
      return

 37   continue
c cross section for Danielewicz forward delay
c...MM->M'
      if(CTOption(34).eq.2) then
         call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sig,gam)
         sig=CTParam(58)*sig
       elseif(CTOption(34).eq.3) then
          m3=e
          call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sig,gam)
          if(sig.gt.1d-5) sig=CTParam(58)
       elseif(CTOption(34).eq.4) then
          m3=e
          call anndec(0,m1,i1,iz1,m2,i2,iz2,e,sig,gam)
          if(sig.gt.1d-5)then
             ggam=fwidth(itypnew(1),i3new(1),m3)
             sig=CTParam(58)*
     &            sig*dsqrt(2.d0/(ggam*3.1415d0*
     &           fbwnorm(m3,itypnew(1),i3new(1))))
          endif
      else
         sig=0.d0
      endif
      return
 38   continue
c elastic meson-meson cross section
      sig=5.d0

c     check energy conservation (cut-off is one 1 MeV)
      if(max(m1,mminit(i1))+max(m2,mminit(i2))+1.d-3.gt.e) sig=0.d0

      return

 40   continue
c M_c meson elastic
c Ziwei Lin, C M Ko, J.Phys. G:Nucl.Part.Phys. 27 (2001) 617-623
c Some are assumed to be a constant 5 mb.

      sig=0d0
      icharm=max(iabs(i1),iabs(i2))
      ilight=min(iabs(i1),iabs(i2))

c     JSFC-------------------------------------
      
      if(iabs(i1).ge.133.and.iabs(i2).lt.133)then
c     C + light
         izcharm=iz1
         izlight=iz2
      endif
      if(iabs(i1).lt.133.and.iabs(i2).ge.133)then
c     Light + C
         izcharm=iz2
         izlight=iz1
      endif      
      if(iabs(i1).ge.133.and.iabs(i2).ge.133)then
c     C+C
c     JSF 
         sig=0.0d0
         return
c     sig=5.0d0
c     return      
      endif
      
      if(iabs(icharm).eq.133.and.iabs(ilight).eq.101) then
c     D + pi
         if(e.le.2.005002)then
            sig=0.0d0
            return
         endif
         if(e.gt.2.80d0)then
            iecm=999
            aecm=0.0d0
            decm=0.8d0/999.0d0
         else
            decm=0.8d0/999.0d0
            iecm=int((e-2.005001d0)/decm)+1
            aecm=(e-2.005001d0-iecm*decm)/decm
         endif

         if(izcharm*izlight.eq.2)then 
c     first D+pi+ D-pi-
            sig=aecm*dppip(iecm+1)+(1.0d0-aecm)*dppip(iecm)
            if(e.gt.2.80d0)then
               slope=-dabs((dlog(dppip(999))-dlog(dppip(998)))/decm)
               sig=sig*dexp(slope*(e-2.80d0))
            endif
            return
         else if(izcharm*izlight.eq.-2)then 
c     D+pi- D-pi+
            sig=aecm*dppim(iecm+1)+(1.0d0-aecm)*dppim(iecm)
            if(e.gt.2.80d0)then
               slope=-dabs((dlog(dppim(999))-dlog(dppim(998)))/decm)
               sig=sig*dexp(slope*(e-2.80d0))
            endif
            return
         else if(izlight.eq.0)then 
c     D0pi0 D-pi0 D+pi0
            sig=aecm*dppi0(iecm+1)+(1.0d0-aecm)*dppi0(iecm)
            if(e.gt.2.80d0)then
               slope=-dabs((dlog(dppi0(999))-dlog(dppi0(998)))/decm)
               sig=sig*dexp(slope*(e-2.80d0))
            endif
            return
         else if(izcharm.eq.0.and.i1*i2*izlight.gt.0)then 
c     D0pi+ Dbar0pi- 
            sig=aecm*dzpip(iecm+1)+(1.0d0-aecm)*dzpip(iecm)
            if(e.gt.2.80d0)then
               slope=-dabs((dlog(dzpip(999))-dlog(dzpip(998)))/decm)
               sig=sig*dexp(slope*(e-2.80d0))
            endif
            return
         else if(izcharm.eq.0.and.i1*i2*izlight.lt.0)then 
c     D0pi- Dbar0pi+ 
            sig=aecm*dzpim(iecm+1)+(1.0d0-aecm)*dzpim(iecm)
            if(e.gt.2.80d0)then
               slope=-dabs((dlog(dzpim(999))-dlog(dzpim(998)))/decm)
               sig=sig*dexp(slope*(e-2.80d0))
            endif
            return
         else
            sig=0.0d0
            return
         endif
      endif

      
      if(icharm.eq.138.and.ilight.eq.101) then
c     pi - D_s elastic
c     JS
         sig=0.0d0
         return
cccccFSF
         sigma1= 9.21209 - 1.29031 * exp(-(e-2.06027)**2/1.95275d-3)
         sigma2= 7.4394 + 1985.02 * exp(-3.04699*e)
         if(e.lt.2.3)then
            if(e.lt.2.15)then
               sig=sigma1
            else
               a=(e-2.15)/(2.3-2.15)
               sig=a*sigma2+(1-a)*sigma1
            endif
         else
            sig=sigma2
         endif 
         
         if(sig.lt.0)then
            sig=0
         endif 
         if(sig.gt.30)then
            sig=30
         endif
         return
      endif

      if(icharm.eq.134.and.ilight.eq.101) then
c     D* + pi
         if(e.le.2.147181)then
            sig=0.0d0
            return
         endif
         if(e.gt.2.947)then
            iecm=999
            aecm=0.0d0
            decm=0.8d0/999.0d0
         else
            decm=0.8d0/999.0d0
            iecm=int((e-2.147181d0)/decm)+1
            aecm=(e-2.147181d0-iecm*decm)/decm
         endif

         if(izcharm*izlight.eq.2)then 
c     first D*+pi+ D*-pi-
            sig=aecm*dsppip(iecm+1)+(1.0d0-aecm)*dsppip(iecm)
            if(e.gt.2.947d0)then
               slope=-dabs((dlog(dsppip(999))-dlog(dsppip(998)))/decm)
               sig=sig*dexp(slope*(e-2.947d0))
            endif
            return
         else if(izcharm*izlight.eq.-2)then 
c     D+pi- D-pi+
            sig=aecm*dsppim(iecm+1)+(1.0d0-aecm)*dsppim(iecm)
            if(e.gt.2.947d0)then
               slope=-dabs((dlog(dsppim(999))-dlog(dsppim(998)))/decm)
               sig=sig*dexp(slope*(e-2.947d0))
            endif
            return
         else if(izlight.eq.0)then 
c     D0pi0 D-pi0 D+pi0
            sig=aecm*dsppi0(iecm+1)+(1.0d0-aecm)*dsppi0(iecm)
            if(e.gt.2.947d0)then
               slope=-dabs((dlog(dsppi0(999))-dlog(dsppi0(998)))/decm)
               sig=sig*dexp(slope*(e-2.947d0))
            endif
            return
         else if(izcharm.eq.0.and.i1*i2*izlight.gt.0)then 
c     D0pi+ Dbar0pi- 
            sig=aecm*dszpip(iecm+1)+(1.0d0-aecm)*dszpip(iecm)
            if(e.gt.2.947d0)then
               slope=-dabs((dlog(dszpip(999))-dlog(dszpip(998)))/decm)
               sig=sig*dexp(slope*(e-2.947d0))
            endif
            return
         else if(izcharm.eq.0.and.i1*i2*izlight.lt.0)then 
c     D0pi- Dbar0pi+ 
            sig=aecm*dszpim(iecm+1)+(1.0d0-aecm)*dszpim(iecm)
            if(e.gt.2.947d0)then
               slope=-dabs((dlog(dszpim(999))-dlog(dszpim(998)))/decm)
               sig=sig*dexp(slope*(e-2.947d0))
            endif
            return
         else
            sig=0.0d0
            return
         endif
      endif

      if(icharm.eq.133.and.ilight.eq.106)then
c     D + K
         if(e.le.2.362001)then
            sig=0.0d0
            return
         endif
         if(e.gt.3.162)then
            iecm=999
            aecm=0.0d0
            decm=0.8d0/999.0d0
         else
            decm=0.8d0/999.0d0
            iecm=int((e-2.362001d0)/decm)+1
            aecm=(e-2.362001-iecm*decm)/decm
         endif
         
         if((izcharm*izlight.eq.1).and.(i1*i2.gt.0))then 
c     first D+K+ D-K-
            sig=aecm*dpKp(iecm+1)+(1.0d0-aecm)*dpKp(iecm)
            if(e.gt.3.162d0)then
               slope=-dabs((dlog(dpKp(999))-dlog(dpKp(998)))/decm)
               sig=sig*dexp(slope*(e-3.162d0))
            endif
            return  
         else if((izcharm*izlight.eq.-1).and.(i1*i2.gt.0))then
c     D++K0 D0+K+
            sig=aecm*DpKz(iecm+1)+(1.0d0-aecm)*DpKz(iecm)
            if(e.gt.3.162d0)then
               slope=-dabs((dlog(DpKz(999))-dlog(DpKz(998)))/decm)
               sig=sig*dexp(slope*(e-3.162d0))
            endif
            return
         else if((izcharm*izlight.eq.-1).and.(i1*i2.lt.0))then
c     D++K- D0+Kbar0
            sig=aecm*DpKm(iecm+1)+(1.0d0-aecm)*DpKm(iecm)
            if(e.gt.3.162d0)then
               slope=-dabs((dlog(dpKm(999))-dlog(dpKm(998)))/decm)
               sig=sig*dexp(slope*(e-3.162d0))
            endif
            return
         else if((izcharm*izlight.eq.1).and.(i1*i2.lt.0))then
c     D++aK0 D0+K-
            sig=aecm*DpaKz(iecm+1)+(1.0d0-aecm)*DpaKz(iecm)
            if(e.gt.3.162d0)then
               slope=-dabs((dlog(dpaKz(999))-dlog(dpaKz(998)))/decm)
               sig=sig*dexp(slope*(e-3.162d0))
            endif
            return            
         endif            
      endif

      if(icharm.eq.134.and.ilight.eq.106)then
c     D* + K
         if(e.le.2.504181)then
            sig=0.0d0
            return
         endif
         if(e.gt.3.3041)then
            iecm=999
            aecm=0.0d0
            decm=0.8d0/999.0d0
         else
            decm=0.8d0/999.0d0
            iecm=int((e-2.504181d0)/decm)+1
            aecm=(e-2.504181-iecm*decm)/decm
         endif
         
         if((izcharm*izlight.eq.1).and.(i1*i2.gt.0))then 
c     first D+K+ D-K-
            sig=aecm*dspKp(iecm+1)+(1.0d0-aecm)*dspKp(iecm)
            if(e.gt.3.3041d0)then
               slope=-dabs((dlog(dspKp(999))-dlog(dspKp(998)))/decm)
               sig=sig*dexp(slope*(e-3.3041d0))
            endif
            return  
         else if((izcharm*izlight.eq.-1).and.(i1*i2.gt.0))then
c     D++K0 D0+K+
            sig=aecm*DspKz(iecm+1)+(1.0d0-aecm)*DspKz(iecm)
            if(e.gt.3.3041d0)then
               slope=-dabs((dlog(dspKz(999))-dlog(dspKz(998)))/decm)
               sig=sig*dexp(slope*(e-3.3041d0))
            endif
            return
         else if((izcharm*izlight.eq.-1).and.(i1*i2.lt.0))then
c     D++K- D0+Kbar0
            sig=aecm*DspKm(iecm+1)+(1.0d0-aecm)*DspKm(iecm)
            if(e.gt.3.3041d0)then
               slope=-dabs((dlog(dspKm(999))-dlog(dspKm(998)))/decm)
               sig=sig*dexp(slope*(e-3.3041d0))
            endif
            return
         else if((izcharm*izlight.eq.1).and.(i1*i2.lt.0))then
c     D++aK0 D0+K-
            sig=aecm*DspaKz(iecm+1)+(1.0d0-aecm)*DspaKz(iecm)
            if(e.gt.3.3041d0)then
               slope=-dabs((dlog(dspaKz(999))-dlog(dspaKz(998)))/decm)
               sig=sig*dexp(slope*(e-3.3041d0))
            endif
            return            
         endif            
      endif

      if(icharm.eq.133.and.ilight.eq.102)then
c         D+eta

         if(e.le.2.415001)then
            sig=0.0d0
            return
         endif
         if(e.gt.3.21500)then
            iecm=999
            aecm=0.0d0
            decm=0.8d0/999.0d0
         else
            decm=0.8d0/999.0d0
            iecm=int((e-2.415001d0)/decm)+1
            aecm=(e-2.415001-iecm*decm)/decm
         endif

          sig=aecm*Deta(iecm+1)+(1.0d0-aecm)*Deta(iecm)
          if(e.gt.3.21500d0)then
             slope=-dabs((dlog(deta(999))-dlog(deta(998)))/decm)
             sig=sig*dexp(slope*(e-3.21500d0))
          endif
          return
      endif

      if(icharm.eq.134.and.ilight.eq.102)then
c         D*+eta
         if(e.le.2.557181)then
            sig=0.0d0
            return
         endif
         if(e.gt.3.3571)then
            decm=0.8d0/999.0d0
            iecm=999
            aecm=0.0d0
         else
            decm=0.8d0/999.0d0
            iecm=int((e-2.557181d0)/decm)+1
            aecm=(e-2.557181d0-iecm*decm)/decm
         endif

          sig=aecm*Dseta(iecm+1)+(1.0d0-aecm)*Dseta(iecm)
          if(e.gt.3.3571d0)then
             slope=-dabs((dlog(dseta(999))-dlog(dseta(998)))/decm)
             sig=sig*dexp(slope*(e-3.3571d0))
          endif
          return
      endif

      if((icharm.eq.133.or.icharm.eq.138).and.ilight.eq.104) then
c     rho - D elastic
c     JSF
         sig=0.0d0
         return
cccccccJSF
         
         sigma1=4.87027+1.57437/(e-2.62398)
         sigma2=54.0987-32.1673*e+7.12686*e**2-4.81498d-1*e**3

         if(e.gt.2.5) then
            sig=sigma1
            if(e.gt.3.05) then
               sig=sigma2
            endif
         endif
         
         if(sig.lt.0)then
            sig=0
         endif 
         
         if(sig.gt.30)then
            sig=30
         endif    
         
         return   


      elseif(icharm.eq.139.and.ilight.eq.101) then
c     pi - Ds* elastic

c     JSF
         sig=0.0d0
         return
ccccc JSF
    
         sigma1=8.48692-37.5089*(e-2.45082)**2
         sigma2=7.28792-2.31346/(e+5.96348d-2)
         
         if(e.gt.2.15) then
            if(e.lt.2.74515) then
               a=dsqrt((e-2.15)/(2.74515-2.15))
               sig=sigma2*a+(1-a)*sigma1
            else
               sig=sigma2 
            endif
         endif
         
         if(sig.lt.0)then
            sig=0
         endif 
         
         return   
  
         
      elseif((icharm.eq.134.or.icharm.eq.139).and.ilight.eq.104) then
c     rho - D* elastic

c     JSF
         sig=0.0d0
         return
cccc  JSF

         sigma1=5.17737+1.54699/(e-2.77327)
         sigma2=73.3847-44.8988*e+9.94873*e**2-6.91067d-1*e**3
         
         if(e.gt.2.8) then
            if(e.lt.3.18) then
               sig=sigma1
            else
               sig=sigma2
            endif
         endif

         if(sig.lt.0)then
            sig=0
         endif 
         
         if(sig.gt.30)then
            sig=30
         endif
         
         return   
      endif
      
      if(icharm.ge.133)then
c     JSF
         sig=0.0d0
c     sig=5d0
         return
      endif
      
      return

 41   continue
c pi D -> rho D*
c Z.-w. Lin et al., Nucl. Phys. A 689, 965 (2001) [nucl-th/0006086]
      sig=0d0
c     JSFC-----
      return
c     JSFC-----
      icharm=max(iabs(i1),iabs(i2))
      ilight=min(iabs(i1),iabs(i2))
      if(icharm.ne.abs(133)) return
      if(ilight.ne.101) return

      sigma1= 0.968152-9.34529*(e-3.10002)**2
      sigma2= -0.173629+1.77871/(e-1.66446)

      if(e.lt.3.2) then
        if(e.lt.3.05) then
          sig=sigma1
        else
          a=(e-3.05)/(3.2-3.05)
          sig=a*sigma2+(1-a)*sigma1
        endif
      else
        sig=sigma2
      endif

      if(sig.lt.0)then
        sig=0
      endif


      if(sig.gt.30)then
       sig=30
      endif
      return

 42   continue
c  pi D* -> rho D
c Z.-w. Lin et al., Nucl. Phys. A 689, 965 (2001) [nucl-th/0006086]

      sig=0d0
c     JSFC-----
      return
c     JSFC-----
      ilight=min(iabs(i1),iabs(i2))
      icharm=max(iabs(i1),iabs(i2))
      if(icharm.ne.abs(134)) return
      if(ilight.ne.101) return

      sigma1=2.13419-5.28378d-1*(e-4.79783)**2
      sigma2=4.7814d-2+8.39242d-1/(e-2.38575)

      if(e.le.2.63) then
        sig=sigma1
      elseif(e.lt.3.02) then
        a=dsqrt((e-2.63)/(3.02-2.63))
        sig=a*sigma2+(1-a)*sigma1
      else
        sig=sigma2
      endif

      if(sig.lt.0)then
        sig=0
      endif

      if(sig.gt.30)then
       sig=30
      endif
      return




 43   continue
c rho D -> pi D*
c Z.-w. Lin et al., Nucl. Phys. A 689, 965 (2001) [nucl-th/0006086]

      sig=0d0
c     JSFC-----
      return
c     JSFC-----
      icharm=max(iabs(i1),iabs(i2))
      ilight=min(iabs(i1),iabs(i2))
      if(icharm.ne.abs(133)) return
      if(ilight.ne.104) return

      sigma1=2.14042+0.209326/(e-2.6328)
      sigma2=5.16369d-2+8.19877d-1/(e -2.59231)



      if(e.gt.2.62) then
        if(e.lt.2.87) then
          sig=sigma1
        else
          sig=sigma2
        endif
      endif

      if(sig.lt.0)then
        sig=0
      endif

      if(sig.gt.30)then
       sig=30
      endif
      return

 44   continue
c rho D* -> pi D
c Z.-w. Lin et al., Nucl. Phys. A 689, 965 (2001) [nucl-th/0006086]

      sig=0d0
c     JSFC-----
      return
c     JSFC-----
      icharm=max(iabs(i1),iabs(i2))
      ilight=min(iabs(i1),iabs(i2))
      if(icharm.ne.abs(134)) return
      if(ilight.ne.abs(104)) return

      sigma1=-779.174+385.407/(e-2.30402)
      sigma2=1.20587d-1+3.26789d-2/(e-2.78287)

      if(e.gt.2.76) then
        if(e.lt.2.796) then
          sig=sigma1
        else
          sig=sigma2
        endif
      endif

      if(sig.lt.0)then
        sig=0
      endif

      if(sig.gt.30)then
       sig=30
      endif
      return


 45   continue
c exit channel: D + DBar
c E. Bratkovskaya et al., Phys. Rev. C 67, 054905 (2003) [nucl-th/0301083]


      sig=0d0
c     JSFC-----
      return
c     JSFC-----

      icharm=abs(max(iabs(i1),iabs(i2)))
      ilight=abs(min(iabs(i1),iabs(i2)))
      if(icharm.ne.135.and.icharm.ne.136.and.icharm.ne.137)return
      if(ilight.ne.101.and.ilight.ne.104)return
      cfac=0.18

c initial channel: pi,rho + J/Psi,Psi',Chi_c
       if(icharm.eq.135)then
        ma1=massit(135)
        n1=1
       elseif(icharm.eq.136)then
        ma1=massit(136)
        n1=1
       else
        ma1=massit(137)
        n1=1
       endif
       if(ilight.eq.101)then
        ma2=massit(101)
        n2=0
       else
        ma2=massit(104)
        n2=1
       endif
       ma3=massit(133)
       n3=0
       ma4=massit(133)
       n4=0

       ifl=1
       is=0

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig=0d0
      else
       sig=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif
      return

 46   continue
c exit channel: D* + DBar, D+D*Bar
c E. Bratkovskaya et al., Phys. Rev. C 67, 054905 (2003) [nucl-th/0301083]



      sig=0d0
c     JSFC-----
      return
c     JSFC-----

      icharm=abs(max(iabs(i1),iabs(i2)))
      ilight=abs(min(iabs(i1),iabs(i2)))
      if(icharm.ne.135.and.icharm.ne.136.and.icharm.ne.137)return
      if(ilight.ne.101.and.ilight.ne.104)return
      cfac=0.18

c initial channel: pi,rho + J/Psi,Psi'
       if(icharm.eq.135)then
        ma1=massit(135)
        n1=1
       elseif(icharm.eq.136)then
        ma1=massit(136)
        n1=1
       else
        ma1=massit(137)
        n1=1
       endif
       if(ilight.eq.101)then
        ma2=massit(101)
        n2=0
       else
        ma2=massit(104)
        n2=1
       endif
       ma3=massit(133)
       n3=0
       ma4=massit(134)
       n4=1

       ifl=1
       is=0

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig=0d0
      else
       sig=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif



      return

 47   continue
c exit channel: D* + D*Bar
c E. Bratkovskaya et al., Phys. Rev. C 67, 054905 (2003) [nucl-th/0301083]



      sig=0d0
c     JSFC-----
      return
c     JSFC-----

      icharm=abs(max(iabs(i1),iabs(i2)))
      ilight=abs(min(iabs(i1),iabs(i2)))
      if(icharm.ne.135.and.icharm.ne.136.and.icharm.ne.137)return
      if(ilight.ne.101.and.ilight.ne.104)return
      cfac=0.18

c initial channel: pi,rho + J/Psi,Psi',Chi_c

       if(icharm.eq.135)then
        ma1=massit(135)
        n1=1
       elseif(icharm.eq.136)then
        ma1=massit(136)
        n1=1
       else
        ma1=massit(137)
        n1=1
       endif
       if(ilight.eq.101)then
        ma2=massit(101)
        n2=0
       else
        ma2=massit(104)
        n2=1
       endif
       ma3=massit(134)
       n3=1
       ma4=massit(134)
       n4=1

       ifl=1
       is=0

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig=0d0
      else
       sig=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif

      return

 48   continue
c exit channel: Ds + DBar, DsBar + D
c E. Bratkovskaya et al., Phys. Rev. C 67, 054905 (2003) [nucl-th/0301083]



      sig=0d0
c     JSFC-----
      return
c     JSFC-----
      icharm=abs(max(iabs(i1),iabs(i2)))
      ilight=abs(min(iabs(i1),iabs(i2)))
      if(icharm.ne.135.and.icharm.ne.136.and.icharm.ne.137)return
      if(ilight.ne.106.and.ilight.ne.108)return

      cfac=0.18

c initial channel: K,K* + J/Psi,Psi',Chi_c

       if(icharm.eq.135)then
        ma1=massit(135)
        n1=1
       elseif(icharm.eq.136)then
        ma1=massit(136)
        n1=1
       else
        ma1=massit(137)
        n1=1
       endif
       if(ilight.eq.106)then
        ma2=massit(106)
        n2=0
       else
        ma2=massit(108)
        n2=1
       endif
       ma3=massit(138)
       n3=0
       ma4=massit(133)
       n4=0

       ifl=1
       is=1

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig=0d0
      else
       sig=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif

      return

 49   continue
c exit channel: Ds + D*Bar, DsBar + D*
c E. Bratkovskaya et al., Phys. Rev. C 67, 054905 (2003) [nucl-th/0301083]



      sig=0d0
c     JSFC-----
      return
c     JSFC-----

      icharm=abs(max(iabs(i1),iabs(i2)))
      ilight=abs(min(iabs(i1),iabs(i2)))
      if(icharm.ne.135.and.icharm.ne.136.and.icharm.ne.137)return
      if(ilight.ne.106.and.ilight.ne.108)return
      cfac=0.18

c initial channel: K,K* + J/Psi,Psi',Chi_c

       if(icharm.eq.135)then
        ma1=massit(135)
        n1=1
       elseif(icharm.eq.136)then
        ma1=massit(136)
        n1=1
       else
        ma1=massit(137)
        n1=1
       endif
       if(ilight.eq.106)then
        ma2=massit(106)
        n2=0
       else
        ma2=massit(108)
        n2=1
       endif
       ma3=massit(138)
       n3=0
       ma4=massit(134)
       n4=1

       ifl=1
       is=1

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig=0d0
      else
       sig=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif
      return

 50   continue
c exit channel: Ds* + DBar, Ds*Bar + D
c E. Bratkovskaya et al., Phys. Rev. C 67, 054905 (2003) [nucl-th/0301083]



      sig=0d0
c     JSFC-----
      return
c     JSFC-----

      icharm=abs(max(iabs(i1),iabs(i2)))
      ilight=abs(min(iabs(i1),iabs(i2)))
      if(icharm.ne.135.and.icharm.ne.136.and.icharm.ne.137)return
      if(ilight.ne.106.and.ilight.ne.108)return
      cfac=0.18

c initial channel: K,K* + J/Psi,Psi'

       if(icharm.eq.135)then
        ma1=massit(135)
        n1=1
       elseif(icharm.eq.136)then
        ma1=massit(136)
        n1=1
       else
        ma1=massit(137)
        n1=1
       endif
       if(ilight.eq.106)then
        ma2=massit(106)
        n2=0
       else
        ma2=massit(108)
        n2=1
       endif
       ma3=massit(133)
       n3=0
       ma4=massit(139)
       n4=1

       ifl=1
       is=1

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig=0d0
      else
       sig=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif

      return

 51   continue
c exit channel: pi + J/Psi
c E. Bratkovskaya et al., Phys. Rev. C 67, 054905 (2003) [nucl-th/0301083]


      sig=0d0
c     JSFC-----
      return
c     JSFC-----

      icharm=abs(max(iabs(i1),iabs(i2)))
      ilight=abs(min(iabs(i1),iabs(i2)))
      if(icharm.ne.133.and.icharm.ne.134)return
      if(ilight.ne.133.and.ilight.ne.134)return
      if(i1.gt.0.and.i2.gt.0)return
      if(i1.lt.0.and.i2.lt.0)return

      cfac=0.18


c initial channel: D(Bar),D*(Bar) + D(Bar),D*(Bar)
       if(icharm.eq.133)then
        ma3=massit(133)
        n3=0
       else
        ma3=massit(134)
        n3=1
       endif
       if(ilight.eq.133)then
        n4=0
        ma4=massit(133)
       else
        ma4=massit(134)
        n4=1
       endif

c..J/Psi
       ma1=massit(101)
       n1=0
       ma2=massit(135)
       n2=1

       ifl=2
       is=0

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig1=0d0
      else
       sig1=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif

c..Psi'
       ma1=massit(101)
       n1=0
       ma2=massit(137)
       n2=1

       ifl=2
       is=0

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig2=0d0
      else
       sig2=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif

c..Chi_c
       ma1=massit(101)
       n1=0
       ma2=massit(136)
       n2=1

       ifl=2
       is=0

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig3=0d0
      else
       sig3=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif

      sig=sig1+sig2+sig3

      if(ranff(0).lt.sig1/sig)then
       cross(51)=1
      elseif(ranff(0).lt.sig2/(sig2+sig3))then
       cross(51)=2
      else
       cross(51)=3
      endif

      return

 52   continue
c exit channel: rho + J/Psi
c E. Bratkovskaya et al., Phys. Rev. C 67, 054905 (2003) [nucl-th/0301083]


      sig=0d0
c     JSFC-----
      return
c     JSFC-----

      icharm=abs(max(iabs(i1),iabs(i2)))
      ilight=abs(min(iabs(i1),iabs(i2)))
      if(icharm.ne.133.and.icharm.ne.134)return
      if(ilight.ne.133.and.ilight.ne.134)return
      if(i1.gt.0.and.i2.gt.0)return
      if(i1.lt.0.and.i2.lt.0)return
      cfac=0.18


c initial channel: D(Bar),D*(Bar) + D(Bar),D*(Bar)
       if(icharm.eq.133)then
        ma3=massit(133)
        n3=0
       else
        ma3=massit(134)
        n3=1
       endif
       if(ilight.eq.133)then
        n4=0
        ma4=massit(133)
       else
        n4=1
        ma4=massit(134)
       endif

c..J/Psi
        n1=1
        n2=1
       ma1=massit(104)
       ma2=massit(135)

       ifl=2
       is=0

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig1=0d0
      else
       sig1=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif
c..Psi'
        n1=1
        n2=1
       ma1=massit(104)
       ma2=massit(137)

       ifl=2
       is=0

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig2=0d0
      else
       sig2=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif
c..Chi_c
        n1=1
        n2=1
       ma1=massit(104)
       ma2=massit(136)

       ifl=2
       is=0

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig3=0d0
      else
       sig3=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif

      sig=sig1+sig2+sig3

      if(ranff(0).lt.sig1/sig)then
       cross(52)=1
      elseif(ranff(0).lt.sig2/(sig2+sig3))then
       cross(52)=2
      else
       cross(52)=3
      endif

      return


 53   continue
c exit channel: K + J/Psi
c E. Bratkovskaya et al., Phys. Rev. C 67, 054905 (2003) [nucl-th/0301083]


      sig=0d0
c     JSFC-----
      return
c     JSFC-----


      icharm=abs(max(iabs(i1),iabs(i2)))
      ilight=abs(min(iabs(i1),iabs(i2)))
      if(icharm.ne.138.and.icharm.ne.139)return
      if(ilight.ne.133.and.ilight.ne.134)return
      if(i1.gt.0.and.i2.gt.0)return
      if(i1.lt.0.and.i2.lt.0)return
      cfac=0.18


c initial channel: Ds(Bar),Ds*(Bar) + D(Bar),D*(Bar)
       if(icharm.eq.138)then
        ma3=massit(138)
        n3=0
       else
        ma3=massit(139)
        n3=1
       endif
       if(ilight.eq.133)then
        ma4=massit(133)
        n4=0
       else
        ma4=massit(134)
        n4=1
       endif

c..J/Psi
       ma2=massit(106)
        n2=0
       ma1=massit(135)
        n1=1

       ifl=2
       is=1

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig1=0d0
      else
       sig1=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif


c..Psi'
       ma2=massit(106)
        n2=0
       ma1=massit(137)
        n1=1

       ifl=2
       is=1

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig2=0d0
      else
       sig2=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif

c..Chi_c
       ma2=massit(106)
        n2=0
       ma1=massit(136)
        n1=1

       ifl=2
       is=1

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig3=0d0
      else
       sig3=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif

      sig=sig1+sig2+sig3

      if(ranff(0).lt.sig1/sig)then
       cross(53)=1
      elseif(ranff(0).lt.sig2/(sig2+sig3))then
       cross(53)=2
      else
       cross(53)=3
      endif

      return

 54   continue
c exit channel: K* + J/Psi
c E. Bratkovskaya et al., Phys. Rev. C 67, 054905 (2003) [nucl-th/0301083]


      sig=0d0
c     JSFC-----
      return
c     JSFC-----

      icharm=abs(max(iabs(i1),iabs(i2)))
      ilight=abs(min(iabs(i1),iabs(i2)))
      if(icharm.ne.138.and.icharm.ne.139)return
      if(ilight.ne.133.and.ilight.ne.134)return
      if(i1.gt.0.and.i2.gt.0)return
      if(i1.lt.0.and.i2.lt.0)return
      cfac=0.18

c initial channel: Ds(Bar),Ds*(Bar) + D(Bar),D*(Bar)
       if(icharm.eq.138)then
        ma3=massit(138)
        n3=0
       else
        ma3=massit(139)
        n3=1
       endif
       if(ilight.eq.133)then
        ma4=massit(133)
        n3=0
       else
        ma4=massit(134)
        n3=1
       endif

c..J/Psi
       ma1=massit(108)
        n3=1
       ma2=massit(135)
        n3=1

       ifl=2
       is=1

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig1=0d0
      else
       sig1=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif

c..Psi'
       ma1=massit(108)
        n3=1
       ma2=massit(137)
        n3=1

       ifl=2
       is=1

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig2=0d0
      else
       sig2=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif

c..Chi_c
       ma1=massit(108)
        n3=1
       ma2=massit(136)
        n3=1

       ifl=2
       is=1

      if(e.lt.(ma1+ma2).or.e.lt.(ma3+ma4))then
        sig3=0d0
      else
       sig3=DMCROSS(e,ma1,ma2,ma3,ma4,n1,n2,n3,n4,ifl,is)
      endif

      sig=sig1+sig2+sig3

      if(ranff(0).lt.sig1/sig)then
       cross(54)=1
      elseif(ranff(0).lt.sig2/(sig2+sig3))then
       cross(54)=2
      else
       cross(54)=3
      endif

      return


 55   continue
c Baryon Charmonium -> Baryon elastic


c Baryon C-Meson -> Baryon elastic

      sig=0d0


c     JSFC-------------------------------------
      
      if(iabs(i1).ge.133.and.iabs(i2).lt.133)then
c     C + light
         izcharm=iz1
         izlight=iz2
         ilight=i2
         icharm=i1
      endif
      if(iabs(i1).lt.133.and.iabs(i2).ge.133)then
c     Light + C
         izcharm=iz2
         izlight=iz1
         ilight=i1
         icharm=i2
      endif      
      
      if(iabs(icharm).eq.133.and.iabs(ilight).eq.1) then
c     D + n
         if(icharm*ilight.gt.0)then
c     D scatt
            if(e.le.2.810d0)then
               sig=0.0d0
               return
            endif
            if(e.gt.3.599d0)then
               iecm=7900
               aecm=0.0d0
               decm=0.0001d0
            else
               decm=0.0001d0
               iecm=int((e-2.810d0)/decm)+1
               aecm=(e-2.810d0-iecm*decm)/decm
            endif
            
            if(izcharm*izlight.gt.0)then 
c     first D0+n D++p
               sig=aecm*Dzn(iecm+1)+(1.0d0-aecm)*Dzn(iecm)
               if(e.gt.3.599d0)then
                  slope=-dabs((dlog(Dzn(7900))-dlog(Dzn(7899)))/decm)
                  sig=sig*dexp(slope*(e-3.599d0))
               endif
               return
            else if(izcharm*izlight.lt.0)then 
c     D0+p D++n
               sig=aecm*Dzp(iecm+1)+(1.0d0-aecm)*Dzp(iecm)
               if(e.gt.3.599d0)then
                  slope=-dabs((dlog(Dzp(7900))-dlog(Dzp(7899)))/decm)
                  sig=sig*dexp(slope*(e-3.599d0))
               endif
               return
            endif
         endif

         if(icharm*ilight.lt.0)then
c     Dbar scatt
            if(e.le.2.8092d0)then
               sig=0.0d0
               return
            endif
            if(e.gt.3.599d0)then
               iecm=7908
               aecm=0.0d0
               decm=0.0001d0
            else
               decm=0.0001d0
               iecm=int((e-2.8092d0)/decm)+1
               aecm=(e-2.8092d0-iecm*decm)/decm
            endif
            
            if(izcharm*izlight.gt.0)then 
c     first D-+n Db0+p
               sig=aecm*Dmn(iecm+1)+(1.0d0-aecm)*Dmn(iecm)
               if(e.gt.3.599d0)then
                  slope=-dabs((dlog(Dmn(7908))-dlog(Dmn(7907)))/decm)
                  sig=sig*dexp(slope*(e-3.599d0))
               endif
               return
            else if(izcharm*izlight.lt.0)then 
c     D-+p Db0+n
               sig=aecm*Dmp(iecm+1)+(1.0d0-aecm)*Dmp(iecm)
               if(e.gt.3.599d0)then
                  slope=-dabs((dlog(Dmp(7908))-dlog(Dmp(7907)))/decm)
                  sig=sig*dexp(slope*(e-3.599d0))
               endif
               return
            endif
         endif
      endif

      if(iabs(icharm).eq.133.and.iabs(ilight).eq.17) then
c     D + Delta
         if(icharm*ilight.gt.0)then
c     D scatt
            if(e.le.3.10d0)then
               sig=0.0d0
               return
            endif
            if(e.gt.3.599d0)then
               iecm=4999
               aecm=0.0d0
               decm=0.0001d0
            else
               decm=0.0001d0
               iecm=int((e-3.10d0)/decm)+1
               aecm=(e-3.10d0-iecm*decm)/decm
            endif
            
            if(izcharm*izlight.eq.3)then 
c     first D+ + Del++ 
               sig=aecm*DpDelpp(iecm+1)+(1.0d0-aecm)*DpDelpp(iecm)
               if(e.gt.3.599d0)then
                  slope=-dabs((dlog(Dpdelpp(4999))
     $                 -dlog(Dpdelpp(4998)))/decm)
                  sig=sig*dexp(slope*(e-3.599d0))
               endif
               return
            else if(izcharm*izlight.eq.1)then 
c     D+ + Del+ 
               sig=aecm*DpDelp(iecm+1)+(1.0d0-aecm)*DpDelp(iecm)
               if(e.gt.3.599d0)then
                  slope=-dabs((dlog(Dpdelp(4999))
     $                 -dlog(Dpdelp(4998)))/decm)
                  sig=sig*dexp(slope*(e-3.599d0))
               endif
               return
            else if(izcharm*izlight.eq.-1)then 
c     D+ + Del0 
               sig=aecm*DpDelz(iecm+1)+(1.0d0-aecm)*DpDelz(iecm)
               if(e.gt.3.599d0)then
                  slope=-dabs((dlog(Dpdelz(4999))
     $                 -dlog(Dpdelz(4998)))/decm)
                  sig=sig*dexp(slope*(e-3.599d0))
               endif
               return
            else if(izcharm*izlight.eq.-3)then 
c     D+ + Del- 
               sig=aecm*DpDelm(iecm+1)+(1.0d0-aecm)*DpDelm(iecm)
               if(e.gt.3.599d0)then
                  slope=-dabs((dlog(Dpdelm(4999))
     $                 -dlog(Dpdelm(4998)))/decm)
                  sig=sig*dexp(slope*(e-3.599d0))
               endif
               return
            endif         
         endif

         if(icharm*ilight.lt.0)then
c     Dbar scatt
            if(e.le.3.0797d0)then
               sig=0.0d0
               return
            endif
            if(e.gt.3.599d0)then
               iecm=5203
               aecm=0.0d0
               decm=0.0001d0
            else
               decm=0.0001d0
               iecm=int((e-3.0797d0)/decm)+1
               aecm=(e-3.0797d0-iecm*decm)/decm
            endif
            
            if(izcharm*izlight.eq.3)then 
c     first Db0 + Del++ 
               sig=aecm*DbzDelpp(iecm+1)+(1.0d0-aecm)*DbzDelpp(iecm)
               if(e.gt.3.599d0)then
                  slope=-dabs((dlog(Dbzdelpp(5203))
     $                 -dlog(Dbzdelpp(5202)))/decm)
                  sig=sig*dexp(slope*(e-3.599d0))
               endif
               return
            else if(izcharm*izlight.eq.1)then 
c     Dbz + Del+ 
               sig=aecm*DbzDelp(iecm+1)+(1.0d0-aecm)*DbzDelp(iecm)
               if(e.gt.3.599d0)then
                  slope=-dabs((dlog(Dbzdelp(5203))
     $                 -dlog(Dbzdelp(5202)))/decm)
                  sig=sig*dexp(slope*(e-3.599d0))
               endif
               return
            else if(izcharm*izlight.eq.-1)then 
c     Dbz + Del0 
               sig=aecm*DbzDelz(iecm+1)+(1.0d0-aecm)*DbzDelz(iecm)
               if(e.gt.3.599d0)then
                  slope=-dabs((dlog(Dbzdelz(5203))
     $                 -dlog(Dbzdelz(5202)))/decm)
                  sig=sig*dexp(slope*(e-3.599d0))
               endif
               return
            else if(izcharm*izlight.eq.-3)then 
c     Dbz + Del- 
               sig=aecm*DbzDelm(iecm+1)+(1.0d0-aecm)*DbzDelm(iecm)
               if(e.gt.3.599d0)then
                  slope=-dabs((dlog(Dbzdelm(5203))
     $                 -dlog(Dbzdelm(5202)))/decm)
                  sig=sig*dexp(slope*(e-3.599d0))
               endif
               return
            endif         
         endif         
      endif


c     JSF
      sig=0.0d0
cccc  JSF
      return


 56   continue
c Baryon Charmonium -> Baryon
c B. Alessandro et al. (NA50), Eur. Phys. J. C48, 329 (2006) nucl-ex/061201
c R. Arnaldi et al. (NA60), J. Phys. G 32, S51 (2006) nucl-ex/0609039

        sig=0d0
c     JSFC-----
      return
c     JSFC-----

      if(i1.eq.135.or.i2.eq.135.and.e.gt.0)then
       if(i1.le.maxbar.or.i2.le.maxbar)then
        sig=4.18d0
        return
       endif
      endif
      if(i1.eq.137.or.i2.eq.137.and.e.gt.0)then
       if(i1.le.maxbar.or.i2.le.maxbar)then
        sig=7.6d0
        return
       endif
      endif
      if(i1.eq.136.or.i2.eq.136.and.e.gt.0)then
       if(i1.le.maxbar.or.i2.le.maxbar)then
        sig=4.18d0
       return
       endif
      endif

      return

 57   continue
c Strangeness exchange channel
c Kaon + Nucleon -> Lambda/Sigma + Pi
      sig=STREXWrapper(io,e,i1,iz1,i2,iz2)
      return



  58  continue
c Strangeness exchange channel
c pi + Lambda/Sigam -> K + N
      sig=0.d0
      if(e.le.massit(1)+massit(106)) return !enough energy for N + K
      sig=STREXWrapper(io,e,i1,iz1,i2,iz2)
      return

 59   continue
c Strangeness exchange channel
c K + Sigma/Lambda -> Xi + pi
      sig=0d0
      sig=STREXWrapper(io,e,i1,iz1,i2,iz2)
      return

 60   continue
c Strangeness exchange channel
c Xi + pi -> Sigma/Lambda + K
      sig=0d0
      if(e.le.massit(40)+massit(106)) return !enough energy for Sig/Lam + K
      sig=STREXWrapper(io,e,i1,iz1,i2,iz2)
      return

 61   continue
c Strangeness exchange channel
c Xi + K -> Omega + pi
      sig=0d0
      if(e.le.massit(55)+massit(101)) return !enough energy for Omega + pi
      sig=STREXWrapper(io,e,i1,iz1,i2,iz2)
      return

 62   continue
c Strangeness exchange channel
c Omega + pi -> Xi + K
      sig=0d0
      sig=STREXWrapper(io,e,i1,iz1,i2,iz2)
      return

 63   continue
c Strangeness exchange channel
c Lambda/Sigma + Lambda/Sigma -> Xi + N
      sig=0d0
      if(e.le.massit(49)+massit(1)) return !enough energy for Xi + N
      if(CTOption(59).eq.0) return
      sig=STREXWrapper(io,e,i1,iz1,i2,iz2)
      return

 64   continue
c Strangeness exchange channel
c Xi + N -> Lambda/Sigma + Lambda/Sigma
      sig=0d0
      if(CTOption(59).eq.0) return
      sig=STREXWrapper(io,e,i1,iz1,i2,iz2)
      return

 99   continue
c...single diffr. pp
      sig=.68*(1.+36/e**2)*log(0.6+0.1*e**2)
      return

      end


C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine nonres(e,ii1,iiz1,ii2,iiz2,sig)
c
cinput e    : $\sqrt{s}$ of collision
cinput ii1  : ID of particle 1
cinput iiz1 : $2\cdot I_3$ of particle 1
cinput ii2  : ID of particle 2
cinput iiz2 : $2\cdot I_3$ of particle 2
coutput sig : cross section
c
c {\tt nonres} adds s-channel strings for strange meson baryon reactions,
c since
c resonances are not sufficient to fit data -
c {\tt sig} gets modified in {\tt nonres},
C so that $k^- p$ cross section describes data reasonable.
C {\tt nonres} should be only called below {\tt meltpoint}
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      logical b
      real*8 m1,m2,e,sig
      integer strit,i1,i2,iz1,iz2,ii1,ii2,iiz1,iiz2

      b=.false.
      call setizm(ii1,iiz1,0d0,ii2,iiz2,0d0,
     @      i1,iz1,m1,i2,iz2,m2)
      if(iabs(i1).gt.iabs(i2))
     @     call swpizm(i1,iz1,m1,i2,iz2,m2)

C    check wether sigma should be modified
C        1) strange meson required
C        2) nonstrange baryon required
C        3) combination of strange meson + baryon (k- p)
C         or antistrange meson + antibaryon (k+ pbar) required
C
      if(strit(i2).ne.0.and.i1*i2.lt.0.and.strit(i1).eq.0)b=.true.
C   if this condition is fullfilled modify sigma
C   add exponential underground and corrections
C   (two gaussians,constant,low energy cut off)
C       to describe experimental k- p cross section
C    four  cases for different corrections
C
c   in case of groundstate kaons the constant low energy cut off
c   is replaced by a 1/s dependent rise to account for the large
c   strangeness exchange cross section there
      if(b.and.abs(i2).eq.106.and.e.gt.1.433.and.e.lt.1.4738188) then
        sig=sig+5.93763355d0/(e-1.251377)**2
      elseif ((b).and.e.gt.1.433.and.e.lt.1.4738188) then
        sig=sig+120.
      elseif ((b).and.e.ge.1.4738188.and.e.lt.1.485215) then
        sig=sig-1.296457765d7*(e-1.433)**4+
     &        2.160975431d4*(e-1.433)**2+120.
      elseif ((b).and.e.ge.1.485215.and.e.lt.1.977) then
        sig=sig+1.07769d+06*exp(-(6.44463d0*e))-
     &         10.*exp(-(((e-1.644)**2)/0.004))+
     &         10.*exp(-(((e-1.977)**2)/0.004))
      elseif ((b).and.e.ge.1.977.and.e.lt.2.12) then
C keep maximum value of gaussian above e=1.977 GeV
C                         (e=2.12 GeV is meltpoint)
        sig=sig+1.07769d+06*exp(-(6.44463d0*e))+10.
      endif

C
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function sappann(sroot,m1,m2)
c
cinput sroot : $\sqrt{s}$ of collision
cinput m1    : mass of 1st (anti-)baryon
cinput m2    : mass of 2nd (anti-)baryon
c
c  OUTPUT: {\ttsappann} = annihilation cross section for $\bar N N$
c
c  Taken from: P. Koch, C.B. Dover, Phys. Rev. {\bf C40} (1989) 145
c
C####C##1#########2#########3#########4#########5#########6#########7##

      implicit none

      include 'options.f'

      real*8 sroot,m1,m2,s0,a,b,sig0,s,srootnn,snn
      parameter(a=0.05d0,b=0.6d0,sig0=120.d0,s0=3.52)

      if (CTOption(38).eq.1) then
c evaluate the parametrization at the same relative momentum as in nbar-n
         srootnn=snn(sroot,m1,m2)
         s=srootnn**2
      else
c evaluate parametrization now at the same sqrts
         s=sroot**2
      endif
      sappann=sig0*(s0/s)*(a**2*s0/((s-s0)**2+a**2*s0)+b)

      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function sapptot(sroot,m1,m2)
c
cinput sroot : $\sqrt{s}$ of collision
cinput m1    : mass of 1st (anti-)baryon
cinput m2    : mass of 2nd (anti-)baryon
c
c  OUTPUT: {\ttsapptot} = total cross section for $\bar N N$
c  high energy paramametrization and data taken from PRD 50 (1994)\\
c  $p_{lab} > 5 $GeV:  CERN/HERA parametrization\\
c  0.3 GeV $< p_{lab} <$ 5 GeV: polynomial fit to the data (by C.S.)\\
c  $p_{lab} <0.3$ GeV: another fit, only constrained by sigma annihilation
c
C####C##1#########2#########3#########4#########5#########6#########7##

      implicit none

      include 'options.f'

      real*8 p,sroot,plab,sighera,m1,m2,srootnn,snn

      if (CTOption(38).eq.1) then
c evaluate the parametrization at the same relative momentum as in nbar-n
         srootnn=snn(sroot,m1,m2)
      else
c evaluate parametrization now at the same sqrts
         srootnn=sroot
      endif

      p=plab(0.938d0,0.938d0,srootnn)

      if(p.ge.5.d0)then
         sapptot=sighera(0.938d0,0.938d0,srootnn,4)
        return
      else if(p.ge.0.3d0)then
        sapptot=75.0146d0+43.1276d0/p+2.58298d0/p**2-3.90783d0*p
        return
      else
        sapptot=271.6d0*exp(-(1.1d0*p**2))
        return
      endif

      end


C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function sappela(sroot,m1,m2)
c
cinput sroot : $\sqrt{s}$ of collision
cinput m1    : mass of 1st (anti-)baryon
cinput m2    : mass of 2nd (anti-)baryon
c
c  OUTPUT: {\tt sappela} = elastic cross section for $\bar N N$
c  high energy paramametrization and data taken from PRD 50 (1994)\\
c  $p_{lab} > 5 $GeV:  CERN/HERA parametrization\\
c  0.3 GeV $< p_{lab} <$ 5 GeV: polynomial fit to the data (by C.S.)\\
c  $p_{lab} <0.3$ GeV: no data, set constant
c
C####C##1#########2#########3#########4#########5#########6#########7##

      implicit none

      include 'options.f'

      real*8 p,sroot,plab,sighera,m1,m2,srootnn,snn

      if (CTOption(38).eq.1) then
c evaluate the parametrization at the same relative momentum as in nbar-n
         srootnn=snn(sroot,m1,m2)
      else
c evaluate parametrization now at the same sqrts
         srootnn=sroot
      endif

      p=plab(0.938d0,0.938d0,srootnn)
      if(p.ge.5.d0)then
         sappela=sighera(0.938d0,0.938d0,srootnn,5)
        return
      else if(p.ge.0.3d0)then
        sappela=31.6166d0+18.2842d0/p-1.14896d0/p**2-3.79508d0*p
        return
      else
        sappela=78.6d0
        return
      endif

      end

C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function sappdiff(sroot,m1,m2)
c
cinput sroot : $\sqrt{s}$ of collision
cinput m1    : mass of 1st (anti-)baryon
cinput m2    : mass of 2nd (anti-)baryon
c
c  OUTPUT: {\tt sappdiff} = diffractive cross section for $\bar N N$
c
c  This cross section is totally determined by
c $ \sigma_{diff}=\sigma_{tot}-\sigma_{elast.}-\sigma_{annihil.}$
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      real*8 p,sroot,m1,m2,plab,sappann,sapptot,sappela

      p=plab(0.938d0,0.938d0,sroot)
      if(p.le.0.1d0)then
        sappdiff=0.d0
        return
      else
        sappdiff=max(0.d0,sapptot(sroot,m1,m2)
     &           -sappela(sroot,m1,m2)-sappann(sroot,m1,m2))
        return
      endif
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      blockdata herafits
c  cross section parameters for specific collsion type
c  CERN/HERA fits, taken from PRD 50 (1994)
c  see: function 'sighera'
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      integer nfit
c     integer io
      parameter (nfit=16)
      real*8 a(nfit),b(nfit),n(nfit),c(nfit),d(nfit)
     &       ,p1(nfit),p2(nfit)
      common /HERA/ a,b,n,c,d,p1,p2

c io =          1       2       3       4        5        6
c            pp(tot) pp(ela) pn(tot) app(tot) app(ela) gammap(tot)
c               7       8       9      10       11       12
c         pi+p(tot) pi+p(el) pi-p(tot) pi-p(el) k+p(tot) k+p(el)
c              13      14      15      16
c          k+n(tot) k-p(tot) k-p(ela) k-n(tot)
c
c p1 (p2) give the momentum range (lab momentum) of the fit

      data a/  48.0,   11.9,  47.3,   38.4,     10.2,   0.147,
     @         16.4,     0.,  33.0,   1.76,     18.1,     5.0,
     @         18.7,   32.1,   7.3,   25.2/
      data b/    0.,   26.9,    0.,   77.6,     52.7,      0.,
     @         19.3,   11.4,   14.0,  11.2,       0.,     8.1,
     @           0.,     0.,     0.,    0./
      data n/    0.,  -1.21,    0.,  -0.64,    -1.16,      0.,
     @         -.42,    -.4, -1.36,  -0.64,       0.,    -1.8,
     @           0.,     0.,    0.,     0./
      data c/ 0.522,  0.169, 0.513,   0.26,    0.125,   .0022,
     @          .19,   .079,  .456,   .043,      .26,     .16,
     @          .21,    .66,   .29,    .38/
      data d/ -4.51,  -1.85, -4.27,   -1.2,    -1.28,   -.017,
     @           0.,     0., -4.03,     0.,      -1.,    -1.3,
     @         -.89,   -5.6,  -2.4,   -2.9/
      data p1/   3.,     2.,    3.,     5.,       5.,      3.,
     @           3.,     2.,   1.8,     1.8,       2.,      2.,
     @           2.,     1.75,    3.,    1.8/
c   for (k- p) p1=1.75 GeV/c works reasonable
      data p2/2100.,  2100.,  370.,  1.7d6,    1.7d6,    183.,
     @         340.,   200.,  370.,   360.,     310.,    175.,
     @         310.,   310.,  175.,   310./

      end

C####C##1#########2#########3#########4#########5#########6#########7##
      blockdata reggefits
c  cross section parameters for specific collsion type
c  PDG/regge fits, taken from JPG 33 (2006) p.337
c  see: function 'sighera'
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      integer nfit
c     integer io
      parameter (nfit=16)
      real*8 z(nfit),y1(nfit),y2(nfit),s0,s1,eta1,eta2,bb

      common /REGGE/ z,y1,y2,s0,s1,eta1,eta2,bb

c io =          1       2       3       4        5        6
c            pp(tot) pp(ela) pn(tot) app(tot) app(ela) gammap(tot)
c               7       8       9      10       11       12
c         pi+p(tot) pi+p(el) pi-p(tot) pi-p(el) k+p(tot) k+p(el)
c              13      14      15      16
c          k+n(tot) k-p(tot) k-p(ela) k-n(tot)
c
c p1 (p2) give the momentum range (lab momentum) of the fit

      data z/  35.45,    0.,  35.8,  35.45,       0.,      0.,
     @         20.86,    0., 20.86,     0.,    17.91,      0.,
     @         17.87, 17.91,    0.,  17.87/
      data y1/ 42.53,    0., 40.15,  42.53,       0.,   0.032,
     @         19.24,    0., 19.24,     0.,     7.14,      0.,
     @          5.17,  7.14,    0.,   5.17/
      data y2/ 33.34,    0.,   30., -33.34,       0.,      0.,
     @          6.03,    0., -6.03,     0.,    13.45,      0.,
     @          7.23,-13.45,    0.,  -7.23/
      data s0/ 28.998/
      data s1/ 1./
      data eta1/ 0.458/
      data eta2/ 0.545/
      data bb/ 0.308/

      end



C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function sighera(m1,m2,sroot,io)
c
c
cinput sroot : $\sqrt{s}$ of collision
cinput m1    : mass of 1st (anti-)baryon
cinput m2    : mass of 2nd (anti-)baryon
cinput io    : flag for respective cross section
c
c  OUTPUT: {\tt sighera} = cross section for specific collsion type
c
c {\tt io} can have the following values:
c \begin{tabular}{rl}
c               1 & $pp$ total \\
c               2 & $pp$ elastic \\
c               3 & $pn$ total \\
c               4 & $\bar p p$ total\\
c               5 & $\bar p p$ elastic\\
c               6 & $\gamma p$ (tot) \\
c               7 & $\pi^+ p$ (tot) \\
c               8 & $\pi^+ p$ (el) \\
c               9 & $\pi^- p$ (tot)\\
c              10 & $\ pi- p$ (el) \\
c              11 & $k^+ p$ (tot)\\
c              12 & $k^+ p$ (el)\\
c              13 & $k^+ n$ (tot)\\
c              14 & $k^- p$ (tot)\\
c              15 & $k^- p$ (ela)\\
c              16 & $k^- n$ (tot)\\
c \end{tabular}
c
c This subroutine returns CERN/HERA parametrizations for cross sections
c {\tt p1} and {\tt p2} in the {\tt blockdata} routine
c give the momentum range (lab momentum) of the fit.
c The fits have been taken from
c Phys. Rev. {\bf D50} (1994).
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      include 'coms.f'
      include 'comres.f'
      integer nfit,io
      parameter (nfit=16)
      real*8 a(nfit),b(nfit),n(nfit),c(nfit),d(nfit),sroot,p,plab
     &       ,p1(nfit),p2(nfit),m1,m2,massit
      real*8 z(nfit),y1(nfit),y2(nfit),s0,s1,eta1,eta2,bb,s

      common /REGGE/ z,y1,y2,s0,s1,eta1,eta2,bb
      common /HERA/ a,b,n,c,d,p1,p2

      go to (2,1,2,2,1,2,2,1,2,1,2,1,2,2,1,2),io

c cern/hera cross section parametrization
 1    continue
c      p=plab(m1,m2,sroot)
      if(io.ge.1.and.io.le.5)then
       p=plab(massit(minnuc),massit(minnuc),sroot)
      elseif(io.eq.6)then
       p=plab(massit(minmes),massit(minnuc),sroot)
      elseif(io.ge.7.and.io.le.10)then
       p=plab(massit(pimeson),massit(minnuc),sroot)
      elseif(io.ge.11.and.io.le.16)then
       p=plab(massit(itkaon),massit(minnuc),sroot)
      else
       write(*,*)'#make22 error. sighera called with io=',io
       stop 137
      endif

         if(p.lt.1d-15) then
c if energy conservation is not possible (p.eq.0) then return zero
            sighera=0.d0
            return
         endif
C
      if(p.lt.p1(io))then
        p=p1(io)
C
      elseif (p.gt.p2(io).and.(warn)) then
        write(6,*)'sighera: sroot=',sroot,' high!, io=',io
        write(6,*)'         m1,m2,plab,p2(io)=',m1,m2,p,p2(io)
        write(6,*)'sighera fit used above upper limit (extrapolation)'
      endif
        sighera=a(io)+b(io)*p**n(io)+c(io)*log(p)**2+d(io)*log(p)
      return

c regge cross sections when available (ref. see blockdata reggefit)
 2    continue
      s=sroot**2
      sighera=z(io)+bb*(log(s/s0))**2+y1(io)*(s1/s)**eta1
     &        -y2(io)*(s1/s)**eta2
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function xmelt(e,x1,x2,em,ep)
c
cinput e  :  $\sqrt{s}$ of process
cinput x1 : value at {\tt em}
cinput x2 : value at {\tt ep}
cinput em : lower boundary in $\sqrt{s}$
cinput ep : upper boundary in $\sqrt{s}$
c
c  {\tt xmelt } yields an interpolation beween {\tt x1} and {\tt x2}
c  in the $\sqrt{s}$ range beween {\tt em} and {\tt ep}.
c  For {\tt e < em}  it equals {\tt x2},  if {\tt e > em} it yields
c  a combination of both {\tt x1} and {\tt x2} such that there is a
c  continuous  transition from {\tt x1} to {\tt x2}.
c
c  For parameter b=.false. a linear combination is used,
c  for parameter b=.true. a sin-form is used.
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      real*8 e,x1,x2,em,ep,lam,pi
      logical b
      parameter(b=.false.,pi=.31415927d-1)


      lam=max(0d0,min(1d0,(e-min(em,ep))/abs(ep-em)))
      if(b)lam=5d-1*sin((lam-5d-1)*pi)+5d-1

      xmelt=(1d0-lam)*x1+lam*x2

      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function plab(m1,m2,sroot)
c
cinput sroot : $\sqrt{s}$ of collision
cinput m1    : mass of 1st (anti-)baryon
cinput m2    : mass of 2nd (anti-)baryon
c
c {\tt plab} returns the lab-momentum of particle 1
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      real*8 sroot,m1,m2
      if(sroot-m1-m2.lt.0d0)then
        plab=0d0
      else
        plab=sqrt((sroot**2-(m1+m2)**2)*(sroot**2-(m1-m2)**2))/(2*m2)
      end if
      return
      end
C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function snn(sbb,m1,m2)
c
cinput sbb   : $\sqrt{s}$ of a $BB$ collision
cinput m1    : mass of 1st (anti-)baryon
cinput m2    : mass of 2nd (anti-)baryon
c
c  {\tt snn} returns the  equivalent c.o.m. energy of a $NN$-collision
c                with the same relative momentum
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      real*8 sbb,prel,m1,m2
      prel=sqrt((sbb**2-(m1+m2)**2)*(sbb**2-(m1-m2)**2))/sbb
      snn=sqrt(prel**2+3.52d0)
      return
      end

c New function xsection1(e) used under 42 continue in crossx

        real*8 function xsection1(e)
        implicit none
        real*8 e,xmelt,sig1,sig2
        real*8 a,b,c
        a=12.9723
        b=4.34934
        c=4.02243

        sig2=a*exp(-exp(-b*(e-c)))
        sig1=62.3707-21.8165*e+2.1348*e**2

        xsection1=xmelt(e,sig2,sig1,.41d1,.42d1)
        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  This function returns the partial cross-sections for the
c  strangeness exchange reactions. It was introduced to avoid
c  the introduction of a big number of different processes for
c  the strangeness exchange. At the same time it allows to have
c  all strangeness exchange cross-sections only at one place in
c  the code.
c
c
      recursive function STREXCHANGE(e,i1,iz1,i2,iz2,i3,iz3,i4,iz4)
     &                   result(sig)
      implicit none
      real*8,  intent(in) :: e       !Two particle center of mass energy
      integer, intent(in) :: i1,i2   !Ityp of incoming particles
      integer, intent(in) :: iz1,iz2 !2*Isospin3 of incoming particles
      integer, intent(in) :: i3,i4   !Ityp of outgoing particles
      integer, intent(in) :: iz3,iz4 !2*Isospin3 of outgoing particles
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'options.f'
      real*8 sig,sigaqm1,sigaqm2
      real*8 dummy
      real*8 massit,pcms
      real*8 DIVpcms

      if((abs(i1).ne.1.and.abs(i1).ne.27.and.abs(i1).ne.40
     &                                  .and.abs(i1).ne.49).or.
     &    abs(iz1).gt.2.or.                   !Nucleon/Lambda/Sigma/Xi

     &   (abs(i2).ne.106.and.abs(i2).ne.27.and.abs(i2).ne.40
     &                  .and.abs(i2).ne.1 .and.abs(i2).ne.49).or.
     &   abs(iz2).gt.2.or.                    !Kaons/Nucleon/Lambda/Sigma/Xi
     &  (abs(i3).ne.27.and.abs(i3).ne.40.and.abs(i3).ne.49
     &                 .and.abs(i3).ne.1.and.abs(i3).ne.55).or.
     &   abs(iz3).gt.2.or.                    !Nucleons/Lambda/Sigma/Xi/Omega
     &   (abs(i4).ne.101.and.abs(i4).ne.1.and.abs(i4).ne.49
     &                  .and.abs(i4).ne.27.and.abs(i4).ne.40).or.
     &   abs(iz4).gt.2) then                  !Pions/Nucleons/Lambda/Sigma/Xi
        write(6,*) "Error: Invalid input parameters in STREXCHANGE"
        write(6,*) e,i1,iz1,i2,iz2,i3,iz3,i4,iz4
        stop 137
      endif


      sig=0.d0
c K- + p -> Y + pi
      if(i1.eq.1.and.iz1.eq.1.and.       !proton
     &   i2.eq.-106.and.iz2.eq.-1) then  !K-
        if(i3.eq.40.and.iz3.eq.2.and.                   !Sigma+
     &     i4.eq.101.and.iz4.eq.-2)then                 !Pi-
           sig=0.0788265d0/(e-1.38841)**2
           if(sig.gt.38) sig=38d0
        elseif(i3.eq.40.and.iz3.eq.-2.and.              !Sigma-
     &     i4.eq.101.and.iz4.eq.2) then                 !Pi+
           sig=0.0196741d0/(e-1.42318)**2
           if(sig.gt.100) sig=100d0
        elseif((i3.eq.40.or.i3.eq.27).and.iz3.eq.0.and. !Sigma0/Lambda0
     &     i4.eq.101.and.iz4.eq.0) then                 !Pi0
           sig=0.0508208/(e-1.38837)**2
           if(sig.gt.24) sig=24d0
        end if
c K- + n -> Y + pi
      else if (i1.eq.1.and.iz1.eq.-1.and. !neutron
     &   i2.eq.-106.and.iz2.eq.-1) then   !K-
        if(i3.eq.40.and.iz3.eq.-2.and.
     &     i4.eq.101.and.iz4.eq.0) then
           !same as 0.5*((K-p->S-pi+)+(K-p->S+pi-))
           sig=5d-1*(STREXCHANGE(e,i1,-iz1,i2,iz2,40,-2,101,2)
     &             + STREXCHANGE(e,i1,-iz1,i2,iz2,40,2,101,-2))
        elseif((i3.eq.40.or.i3.eq.27).and.iz3.eq.0.and.
     &     i4.eq.101.and.iz4.eq.-2) then
           sig=STREXCHANGE(e,i1,-iz1,i2,iz2,40,0,101,0) !equal to K-+p->S0/L0+pi0
        endif
c K+ + ap -> aY + pi
      else if (i1.eq.-1.and.iz1.eq.-1.and. !aProton
     &   i2.eq.106.and.iz2.eq.1) then   !K+
        if(i3.eq.-40.and.iz3.eq.2.and.
     &      i4.eq.101.and.iz4.eq.-2) then
          sig=STREXCHANGE(e,-i1,-iz1,-i2,-iz2,40,-2,101,2) !anti reaction
        elseif(i3.eq.-40.and.iz3.eq.-2.and.
     &    i4.eq.101.and.iz4.eq.2) then
          sig=STREXCHANGE(e,-i1,-iz1,-i2,-iz2,40,2,101,-2) !anti reaction
        elseif((i3.eq.-40.or.i3.eq.-27).and.iz3.eq.0.and.
     &    i4.eq.101.and.iz4.eq.0) then
          sig=STREXCHANGE(e,-i1,-iz1,-i2,-iz2,40,0,101,0)  !anti reaction
        endif
c K+ + an -> aY + pi
      else if (i1.eq.-1.and.iz1.eq.1.and. !aNeutron
     &   i2.eq.106.and.iz2.eq.1) then     !K+
        if(i3.eq.-40.and.iz3.eq.2.and.
     &     i4.eq.101.and.iz4.eq.0) then
           sig=STREXCHANGE(e,-i1,-iz1,-i2,-iz2,40,-2,101,0)  !anti reaction
        elseif((i3.eq.-40.or.i3.eq.-27).and.iz3.eq.0.and.
     &     i4.eq.101.and.iz4.eq.2) then
           sig=STREXCHANGE(e,-i1,-iz1,-i2,-iz2,40,0,101,-2)   !anti reaction
        endif
c K0 + ap
      else if (i1.eq.-1.and.iz1.eq.-1.and.!aProton
     &   i2.eq.106.and.iz2.eq.-1)  then   !K0
        if(i3.eq.-40.and.iz3.eq.-2.and.
     &     i4.eq.101.and.iz4.eq.0) then
          sig=STREXCHANGE(e,-i1,-iz1,-i2,-iz2,40,2,101,0)   !anti reaction
        elseif((i3.eq.-40.or.i3.eq.-27).and.iz3.eq.0.and.
     &     i4.eq.101.and.iz4.eq.-2) then
          sig=STREXCHANGE(e,-i1,-iz1,-i2,-iz2,40,0,101,2)   !anti reaction
        endif
c K0 + an
      else if (i1.eq.-1.and.iz1.eq.1.and. !aNeutron
     &   i2.eq.106.and.iz2.eq.-1)  then   !K0
        if(i3.eq.-40.and.iz3.eq.2.and.
     &    i4.eq.101.and.iz4.eq.-2) then
          sig=STREXCHANGE(e,-i1,-iz1,-i2,-iz2,40,-2,101,2) !anti reaction
        elseif(i3.eq.-40.and.iz3.eq.-2.and.
     &    i4.eq.101.and.iz4.eq.2)  then
          sig=STREXCHANGE(e,-i1,-iz1,-i2,-iz2,40,2,101,-2) !anti reaction
        elseif((i3.eq.-40.or.i3.eq.-27).and.iz3.eq.0.and.
     &    i4.eq.101.and.iz4.eq.0) then
          sig=STREXCHANGE(e,-i1,-iz1,-i2,-iz2,40,0,101,0)  !anti reaction
        endif
c aK0 + p
      else if (i1.eq.1.and.iz1.eq.1.and.  !proton
     &   i2.eq.-106.and.iz2.eq.1)  then   !aK0
        if(i3.eq.40.and.iz3.eq.2.and.
     &     i4.eq.101.and.iz4.eq.0) then
          !same as K- + n -> Sigma- + pi0
          sig=STREXCHANGE(e,i1,-iz1,i2,-iz2,40,-2,101,0)
        elseif((i3.eq.40.or.i3.eq.27).and.iz3.eq.0.and.
     &     i4.eq.101.and.iz4.eq.2) then
          !same as K- + n -> Lambda0/Sigma0 + pi-
          sig=STREXCHANGE(e,i1,-iz1,i2,-iz2,40,0,101,-2)
        endif
c aK0 + n
      else if (i1.eq.1.and.iz1.eq.-1.and. !neutron
     &   i2.eq.-106.and.iz2.eq.1)  then   !aK0
        if(i3.eq.40.and.iz3.eq.2.and.
     &    i4.eq.101.and.iz4.eq.-2)then
          sig=STREXCHANGE(e,1,1,-106,-1,40,2,101,-2) !same as K-+p -> Sigma+ + pi-
        elseif(i3.eq.40.and.iz3.eq.-2.and.
     &    i4.eq.101.and.iz4.eq.2) then
          sig=STREXCHANGE(e,1,1,-106,-1,40,-2,101,2) !same as K-+p -> Sigma- + pi+
        elseif((i3.eq.40.or.i3.eq.27).and.iz3.eq.0.and.
     &    i4.eq.101.and.iz4.eq.0) then
          sig=STREXCHANGE(e,1,1,-106,-1,40,0,101,0)  !same as K-+p -> Sigma0 + pi0
        end if



c aK + Sigma
      else if(i1.eq.40.and.i2.eq.-106) then
        if(iz1.eq.2.and.iz2.eq.1.and.                !Sig+ + aK0 -> Xi0 + pi+
     &     i3.eq.49.and.iz3.eq.1.and.
     &     i4.eq.101.and.iz4.eq.2) then
          call aqm(40,-106,sigaqm1,dummy)
          call aqm(1,106,sigaqm2,dummy)
          sig=3.02143713508091d-2/(e-1.66789841701854d0)**2
     &       *sigaqm1/sigaqm2
        elseif(iz1.eq.2.and.iz2.eq.-1.and.           !Sig+ + K-  -> Xi0 + pi0
     &    ((i3.eq.49.and.iz3.eq.1.and.               !           -> Xi- + pi+
     &      i4.eq.101.and.iz4.eq.0).or.
     &     (i3.eq.49.and.iz3.eq.-1.and.
     &      i4.eq.101.and.iz4.eq.2))) then
          sig=STREXCHANGE(e,40,2,-106,1,49,1,101,2)
        elseif(iz1.eq.0.and.iz2.eq.1.and.            !Sig0 + aK0 -> Xi0 + pi0
     &     ((i3.eq.49.and.iz3.eq.1.and.              !           -> Xi- + pi+
     &       i4.eq.101.and.iz4.eq.0).or.
     &      (i3.eq.49.and.iz3.eq.-1.and.
     &       i4.eq.101.and.iz4.eq.2))) then
          sig=STREXCHANGE(e,40,2,-106,1,49,1,101,2)
        elseif(iz1.eq.0.and.iz2.eq.-1.and.           !Sig0 + K-  -> Xi0 + pi-
     &     ((i3.eq.49.and.iz3.eq.1.and.              !           -> Xi- + pi0
     &       i4.eq.101.and.iz4.eq.-2).or.
     &      (i3.eq.49.and.iz3.eq.-1.and.
     &       i4.eq.101.and.iz4.eq.0))) then
          sig=STREXCHANGE(e,40,2,-106,1,49,1,101,2)
        elseif(iz1.eq.-2.and.iz2.eq.1.and.           !Sig- + aK0 -> Xi- + pi0
     &     ((i3.eq.49.and.iz3.eq.-1.and.             !           -> Xi0 + pi-
     &       i4.eq.101.and.iz4.eq.0).or.
     &      (i3.eq.49.and.iz3.eq.1.and.
     &       i4.eq.101.and.iz4.eq.-2))) then
          sig=STREXCHANGE(e,40,2,-106,1,49,1,101,2)
        elseif(iz1.eq.-2.and.iz2.eq.-1.and.          !Sig- + K-  -> Xi- + pi-
     &     i3.eq.49.and.iz3.eq.-1.and.
     &     i4.eq.101.and.iz4.eq.-2) then
          sig=STREXCHANGE(e,40,2,-106,1,49,1,101,2)
        endif
c K + aSigma
      else if(i1.eq.-40.and.i2.eq.106) then
        if(iz1.eq.-2.and.iz2.eq.-1.and.              !aSig- + K0 -> aXi0 + pi-
     &     i3.eq.-49.and.iz3.eq.-1.and.
     &     i4.eq.101.and.iz4.eq.-2) then
          call aqm(-40,106,sigaqm1,dummy)
          call aqm(1,106,sigaqm2,dummy)
          sig=3.02143713508091d-2/(e-1.66789841701854d0)**2
     &       *sigaqm1/sigaqm2
        elseif(iz1.eq.-2.and.iz2.eq.1.and.           !aSig- + K+ -> aXi0 + pi0
     &     ((i3.eq.-49.and.iz3.eq.-1.and.            !           -> aXi+ + pi-
     &       i4.eq.101.and.iz4.eq.0).or.
     &      (i3.eq.-49.and.iz3.eq.1.and.
     &       i4.eq.101.and.iz4.eq.-2))) then
          sig=STREXCHANGE(e,-40,-2,106,-1,-49,-1,101,-2)
        elseif(iz1.eq.0.and.iz2.eq.-1.and.           !aSig0 + K0 -> aXi0 + pi0
     &     ((i3.eq.-49.and.iz3.eq.-1.and.            !           -> aXi+ + pi-
     &       i4.eq.101.and.iz4.eq.0).or.
     &      (i3.eq.-49.and.iz3.eq.1.and.
     &       i4.eq.101.and.iz4.eq.-2))) then
          sig=STREXCHANGE(e,-40,-2,106,-1,-49,-1,101,-2)
        elseif(iz1.eq.0.and.iz2.eq.1.and.            !aSig0 + K+ -> aXi0 + pi+
     &     ((i3.eq.-49.and.iz3.eq.-1.and.            !           -> aXi+ + pi0
     &       i4.eq.101.and.iz4.eq.2).or.
     &      (i3.eq.-49.and.iz3.eq.1.and.
     &       i4.eq.101.and.iz4.eq.0))) then
          sig=STREXCHANGE(e,-40,-2,106,-1,-49,-1,101,-2)
        elseif(iz1.eq.2.and.iz2.eq.-1.and.           !aSig+ + K0 -> aXi+ + pi0
     &     ((i3.eq.-49.and.iz3.eq.1.and.             !           -> aXi0 + pi+
     &       i4.eq.101.and.iz4.eq.0).or.
     &      (i3.eq.-49.and.iz3.eq.-1.and.
     &       i4.eq.101.and.iz4.eq.2))) then
          sig=STREXCHANGE(e,-40,-2,106,-1,-49,-1,101,-2)
        elseif(iz1.eq.2.and.iz2.eq.1.and.            !aSig+ + K+ -> aXi+ + pi+
     &     i3.eq.-49.and.iz3.eq.1.and.
     &     i4.eq.101.and.iz4.eq.2) then
          sig=STREXCHANGE(e,-40,-2,106,-1,-49,-1,101,-2)
        endif
c aK + Lambda
      else if(i1.eq.27.and.i2.eq.-106) then
        if(iz1.eq.0.and.iz2.eq.1.and.                !Lam0 + aK0 -> Xi0 + pi0
     &    ((i3.eq.49.and.iz3.eq.1.and.               !           -> Xi- + pi+
     &      i4.eq.101.and.iz4.eq.0).or.
     &     (i3.eq.49.and.iz3.eq.-1.and.
     &      i4.eq.101.and.iz4.eq.2))) then
          call aqm(-27,-106,sigaqm1,dummy)
          call aqm(1,106,sigaqm2,dummy)
          sig=2.90243238374d-2/(e-1.59226052900406d0)**2
     &       *sigaqm1/sigaqm2
        elseif(iz1.eq.0.and.iz2.eq.-1.and.           !Lam0 + K-  -> Xi0 + pi-
     &     ((i3.eq.49.and.iz3.eq.1.and.              !           -> Xi- + pi0
     &       i4.eq.101.and.iz4.eq.-2).or.
     &      (i3.eq.49.and.iz3.eq.-1.and.
     &       i4.eq.101.and.iz4.eq.0))) then
          sig=STREXCHANGE(e,27,0,-106,1,49,1,101,0)
        endif
c K + aLambda
      else if(i1.eq.-27.and.i2.eq.106) then
        if(iz1.eq.0.and.iz2.eq.1.and.                !aLam0 + K+ -> aXi0 + pi+
     &    ((i3.eq.-49.and.iz3.eq.-1.and.             !           -> aXi+ + pi0
     &      i4.eq.101.and.iz4.eq.2).or.
     &     (i3.eq.-49.and.iz3.eq.1.and.
     &      i4.eq.101.and.iz4.eq.0))) then
          call aqm(-27,-106,sigaqm1,dummy)
          call aqm(1,106,sigaqm2,dummy)
          sig=2.90243238374d-2/(e-1.59226052900406d0)**2
     &       *sigaqm1/sigaqm2
        elseif(iz1.eq.0.and.iz2.eq.-1.and.           !aLam0 + K0 -> aXi0 + pi0
     &     ((i3.eq.-49.and.iz3.eq.-1.and.            !           -> aXi+ + pi-
     &       i4.eq.101.and.iz4.eq.0).or.
     &      (i3.eq.-49.and.iz3.eq.1.and.
     &       i4.eq.101.and.iz4.eq.-2))) then
          sig=STREXCHANGE(e,-27,0,106,1,-49,-1,101,2)
        endif




c Xi0 + aK0 -> Omega- + pi+ (same cross-section for all other combinations)
      elseif(i1.eq.49.and.iz1.eq.1.and.
     &   i2.eq.-106.and.iz2.eq.1.and.
     &   i3.eq.55.and.iz3.eq.0.and.
     &   i4.eq.101.and.iz4.eq.2) then
        call aqm(49,-106,sigaqm1,dummy)
        call aqm(1,106,sigaqm2,dummy)
        sig=3.1975566736027d-2/(e-1.79037415950197d0)**2
     &     *sigaqm1/sigaqm2
c Xi0 + K- -> Omega- + pi0
      elseif(i1.eq.49.and.iz1.eq.1.and.
     &   i2.eq.-106.and.iz2.eq.-1.and.
     &   i3.eq.55.and.iz3.eq.0.and.
     &   i4.eq.101.and.iz4.eq.0) then
        sig=STREXCHANGE(e,49,1,-106,1,55,0,101,2)
c Xi- + aK0 -> Omega- + pi0
      elseif(i1.eq.49.and.iz1.eq.-1.and.
     &   i2.eq.-106.and.iz2.eq.1.and.
     &   i3.eq.55.and.iz3.eq.0.and.
     &   i4.eq.101.and.iz4.eq.0) then
        sig=STREXCHANGE(e,49,1,-106,1,55,0,101,2)
c Xi- + K- ->  Omega- + pi-
      elseif(i1.eq.49.and.iz1.eq.-1.and.
     &   i2.eq.-106.and.iz2.eq.-1.and.
     &   i3.eq.55.and.iz3.eq.0.and.
     &   i4.eq.101.and.iz4.eq.-2) then
        sig=STREXCHANGE(e,49,1,-106,1,55,0,101,2)
c aXi0 + K0 -> aOmega+ + pi-
      elseif(i1.eq.-49.and.iz1.eq.-1.and.
     &   i2.eq.106.and.iz2.eq.-1.and.
     &   i3.eq.-55.and.iz3.eq.0.and.
     &   i4.eq.101.and.iz4.eq.-2) then
        sig=STREXCHANGE(e,49,1,-106,1,55,0,101,2)
c aXi0 + K+ -> aOmega+ + pi0
      elseif(i1.eq.-49.and.iz1.eq.-1.and.
     &   i2.eq.106.and.iz2.eq.1.and.
     &   i3.eq.-55.and.iz3.eq.0.and.
     &   i4.eq.101.and.iz4.eq.0) then
        sig=STREXCHANGE(e,49,1,-106,1,55,0,101,2)
c aXi+ + K0 -> aOmega+ + pi0
      elseif(i1.eq.-49.and.iz1.eq.1.and.
     &   i2.eq.106.and.iz2.eq.-1.and.
     &   i3.eq.-55.and.iz3.eq.0.and.
     &   i4.eq.101.and.iz4.eq.0) then
        sig=STREXCHANGE(e,49,1,-106,1,55,0,101,2)
c aXi+ + K+ -> aOmega+ + pi+
      elseif(i1.eq.-49.and.iz1.eq.1.and.
     &   i2.eq.106.and.iz2.eq.1.and.
     &   i3.eq.-55.and.iz3.eq.0.and.
     &   i4.eq.101.and.iz4.eq.2) then
        sig=STREXCHANGE(e,49,1,-106,1,55,0,101,2)


c Cross-sections from private communications with the
c authors of PRC 85, 064902 (2012)
c Lambda + Lambda
c aLambda + aLambda
      elseif(abs(i1).eq.27.and.abs(i2).eq.27.and.
     &       i1*i2.gt.0) then
        if((e-massit(49)-massit(1)).le.0) then
          sig=0d0
        else
          sig=37.15d0*(e-massit(49)-massit(1))**(-0.16)
     &       *pcms(e,massit(49),massit(1))/pcms(e,massit(27),massit(27))
     &       /2d0
        endif
c Lambda0 + Sigma-  -> Xi- + n
c Lambda0 + Sigma+  -> Xi0 + p
c aLambda0 + aSigma+  -> aXi+ + an
c aLambda0 + aSigma-  -> aXi0 + ap
      elseif( (((abs(i1).eq.27.and.abs(i2).eq.40).and.
     &          (iz1.eq.0.and.abs(iz2).eq.2)).or.
     &         ((abs(i1).eq.40.and.abs(i2).eq.27).and.
     &          (abs(iz1).eq.2.and.iz2.eq.0))).and.
     &         i1*i2.gt.0 ) then
        if((e-massit(40)-massit(27)).le.0) then
          sig=0d0
        else
          sig=24.3781d0*(e-massit(40)-massit(27))**(-0.479d0)
        endif
c Lambda0 + Sigma0  -> Xi0 + n
c                   -> Xi- + p
c aLambda0 + aSigma0  -> aXi0 + an
c                     -> aXi+ + ap
      elseif( ((abs(i1).eq.27.and.abs(i2).eq.40).or.
     &         (abs(i1).eq.40.and.abs(i2).eq.27)).and.
     &         i1*i2.gt.0.and.iz1.eq.0.and.iz2.eq.0) then
        if((e-massit(40)-massit(27)).le.0) then
          sig=0d0
        elseif((e-massit(40)-massit(27)).le.3.336d-2) then
          sig=6.475d0*(e-massit(40)-massit(27))**(-0.4167d0)
        else
          sig=14.5054d0*(e-massit(40)-massit(27))**(-0.1795d0)
        endif
c Sigma0  + Sigma0  -> Xi0 + n
c                   -> Xi- + p
c aSigma0  + aSigma0  -> aXi0 + an
c                     -> aXi+ + ap
      elseif(abs(i1).eq.40.and.abs(i2).eq.40.and.
     &       iz1.eq.0.and.iz2.eq.0.and.i1*i2.gt.0) then
        if((e-2*massit(40)).le.0) then
          sig=0d0
        elseif((e-2*massit(40)).le.9.047d-2) then
          sig=5.625d0*(e-2*massit(40))**(-0.318d0)
        else
          sig=4.174d0*(e-2*massit(40))**(-0.4421d0)
        endif
c Sigma0  + Sigma+  -> Xi0 + p
c Sigma0  + Sigma-  -> Xi- + n
c aSigma0  + aSigma-  -> aXi0 + ap
c aSigma0  + aSigma+  -> aXi+ + an
      elseif(abs(i1).eq.40.and.abs(i2).eq.40.and.
     &       ((iz1.eq.0.and.abs(iz2).eq.2).or.
     &        (abs(iz1).eq.2.and.iz2.eq.0)).and.
     &       i1*i2.gt.0) then
        sig=4d0*STREXCHANGE(e,40,0,40,0,49,-1,1,1)
c Sigma+  + Sigma-  -> Xi0 + n
c                   -> Xi- + p
c aSigma-  + aSigma+  -> aXi0 + an
c                     -> aXi+ + ap
      elseif(abs(i1).eq.40.and.abs(i2).eq.40.and.
     &       abs(iz1).eq.2.and.abs(iz2).eq.2.and.
     &       i1*i2.gt.0.and.(iz1+iz2).eq.0) then
        if((e-2*massit(40)).le.0) then
          sig=0d0
        else
          sig=14.194*(e-2*massit(40))**(-0.442)
        endif
c Xi0 + n -> Lambda0 + Lambda0
c Xi0 + n -> Lambda0 + Sigma0
c Xi0 + n -> Sigma0  + Sigma0
c Xi0 + n -> Sigma+  + Sigma-
c Xi- + p -> Lambda0 + Lambda0
c Xi- + p -> Lambda0 + Sigma0
c Xi- + p -> Sigma0  + Sigma0
c Xi- + p -> Sigma+  + Sigma-
c Xi0 + p -> Lambda0 + Sigma+
c Xi0 + p -> Sigma0  + Sigma+
c Xi- + n -> Lambda0 + Sigma-
c Xi- + n -> Sigma0  + Sigma-
c aXi0 + an -> aLambda0 + aLambda0
c aXi0 + an -> aLambda0 + aSigma0
c aXi0 + an -> aSigma0  + aSigma0
c aXi0 + an -> aSigma+  + aSigma-
c aXi+ + ap -> aLambda0 + aLambda0
c aXi+ + ap -> aLambda0 + aSigma0
c aXi+ + ap -> aSigma0  + aSigma0
c aXi+ + ap -> aSigma+  + aSigma-
c aXi0 + ap -> aLambda0 + aSigma+
c aXi0 + ap -> aSigma0  + aSigma+
c aXi+ + an -> aLambda0 + aSigma+
c aXi+ + an -> aSigma0  + aSigma+
      elseif( ((abs(i1).eq.49.and.abs(i2).eq.1).or.
     &         (abs(i1).eq.1.and.abs(i2).eq.49)).and.
     &         i1*i2.gt.0.and.
     &         abs(iz1).eq.1.and.
     &         abs(iz2).eq.1) then
        sig=STREXCHANGE(e,i3,iz3,i4,iz4,i1,iz1,i2,iz2)
     &      *pcms(e,massit(i3),massit(i4))**2
     &   *DIVpcms(e,massit(i1),massit(i2))**2
c        if(abs(i3).eq.40.and.abs(i4).eq.40.and.iz3.ne.iz4) then
c          sig=2d0*sig !symmetrization for channels Xi + N -> Sigma+ Sigma-
                      !                                   -> Sigma+ Sigma0
                      !                                   -> Sigma- Sigma0
c        endif
      endif

      if(sig.lt.0) sig=0d0
      return
      end function

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function DIVpcms(ecm,m1,m2) !Divde by pcms
      implicit none
      real*8 ecm,m1,m2,X
      real*8 pcms
      X=pcms(ecm,m1,m2)
      if(X.ne.0.)then
       DIVpcms=1.d0/X
      else
       DIVpcms=0  ! to avoid /0 in STREXWrapper,  put sig=0 in that case
      endif
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  This function is a wrapper for the STREXCHANGE function. It
c  is needed to avoid calling the subroutine crossx recursively
c  which might results in unwanted side effects. It returns
c  the cross-section for processes 57 or 58.
      real*8 function STREXWrapper(io,e,i1,iz1,i2,iz2)
      implicit none
      integer,  intent(in) :: io      !process identifier (either 57 or 58)
      real*8 ,  intent(in) :: e       !Two particle center of mass energy
      integer,  intent(in) :: i1,i2   !ityp of incoming particles
      integer,  intent(in) :: iz1,iz2 !2*Isospin3 of incoming particles
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'options.f'
c functions
      real*8 massit,pcms
      real*8 STREXCHANGE 
      real*8 DIVpcms
c variables
      real*8 sig
      real*8 sig_1,sig_2,sig_3,sig_4,sig_5



      sig=0.d0
      if(io.eq.57) then
c K- + p -> Y + pi
        if(i1.eq.1.and.iz1.eq.1.and.        !proton
     &     i2.eq.-106.and.iz2.eq.-1) then   !K-
          sig= STREXCHANGE(e,i1,iz1,i2,iz2,40,2,101,-2) !Sigma+ + Pi-
     &       + STREXCHANGE(e,i1,iz1,i2,iz2,40,-2,101,2) !Sigma- + Pi+
     &       + STREXCHANGE(e,i1,iz1,i2,iz2,40,0,101,0)  !Lambda0/Sigma0 + Pi0
c K- + n -> Y + pi
        else if (i1.eq.1.and.iz1.eq.-1.and.  !neutron
     &     i2.eq.-106.and.iz2.eq.-1) then    !K-
          sig= STREXCHANGE(e,i1,iz1,i2,iz2,40,-2,101,0) !Sigma- + Pi0
     &       + STREXCHANGE(e,i1,iz1,i2,iz2,40,0,101,-2) !Lambda0/Sigma0 + Pi-
c K+ + ap -> aY + pi
        else if (i1.eq.-1.and.iz1.eq.-1.and. !aProton
     &     i2.eq.106.and.iz2.eq.1) then      !K+
          sig= STREXCHANGE(e,i1,iz1,i2,iz2,-40,2,101,-2) !aSigma+ + Pi-
     &       + STREXCHANGE(e,i1,iz1,i2,iz2,-40,-2,101,2) !aSigma- + Pi+
     &       + STREXCHANGE(e,i1,iz1,i2,iz2,-40,0,101,0)  !aLambda0/aSigma0 + Pi0
c K+ + an -> aY + pi
        else if (i1.eq.-1.and.iz1.eq.1.and.  !aNeutron
     &     i2.eq.106.and.iz2.eq.1) then      !K+
          sig= STREXCHANGE(e,i1,iz1,i2,iz2,-40,2,101,0)  !aSigma+ + Pi0
     &       + STREXCHANGE(e,i1,iz1,i2,iz2,-40,0,101,2)  !aLambda0/aSigma0 + Pi+
c K0 + ap
        else if (i1.eq.-1.and.iz1.eq.-1.and. !aProton
     &     i2.eq.106.and.iz2.eq.-1)  then    !K0
          sig= STREXCHANGE(e,i1,iz1,i2,iz2,-40,-2,101,0) !aSigma- + Pi0
     &       + STREXCHANGE(e,i1,iz1,i2,iz2,-40,0,101,-2) !aLambda0/aSigma0 + Pi-
c K0 + an
        else if (i1.eq.-1.and.iz1.eq.1.and.  !aNeutron
     &     i2.eq.106.and.iz2.eq.-1)  then    !K0
          sig= STREXCHANGE(e,i1,iz1,i2,iz2,-40,2,101,-2) !aSigma+ + Pi-
     &       + STREXCHANGE(e,i1,iz1,i2,iz2,-40,-2,101,2) !aSigma- + Pi+
     &       + STREXCHANGE(e,i1,iz1,i2,iz2,-40,0,101,0)  !aLambda0/aSigma0 + Pi0
c aK0 + p
        else if (i1.eq.1.and.iz1.eq.1.and.   !proton
     &     i2.eq.-106.and.iz2.eq.1)  then    !aK0
          sig= STREXCHANGE(e,i1,iz1,i2,iz2,40,2,101,0) !Sigma+ + Pi0
     &       + STREXCHANGE(e,i1,iz1,i2,iz2,40,0,101,2) !Lambda0/aSigma0 + Pi+
c aK0 + n
        else if (i1.eq.1.and.iz1.eq.-1.and.  !neutron
     &     i2.eq.-106.and.iz2.eq.1)  then    !aK0
          sig= STREXCHANGE(e,i1,iz1,i2,iz2,40,2,101,-2) !Sigma+ + Pi-
     &       + STREXCHANGE(e,i1,iz1,i2,iz2,40,-2,101,2) !Sigma- + Pi+
     &       + STREXCHANGE(e,i1,iz1,i2,iz2,40,0,101,0)  !Lambda0/Sigma0 + Pi0
        end if




      else if (io.eq.58) then
cUNCHARGED CHANNELS
c Sigma+ + pi-  -> K- + p
c Sigma+ + pi-  -> aK0 +n
        if(i1.eq.40.and.iz1.eq.2.and.       !Sigma+
     &     i2.eq.101.and.iz2.eq.-2) then    !pi-
           sig= STREXCHANGE(e,1,1,-106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
     &        + STREXCHANGE(e,1,-1,-106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
c Sigma- + pi+  -> K- + p
c Sigma- + pi+  -> aK0 +n
        elseif(i1.eq.40.and.iz1.eq.-2.and.  !Sigma-
     &     i2.eq.101.and.iz2.eq.2)  then    !pi+
           sig= STREXCHANGE(e,1,1,-106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
     &        + STREXCHANGE(e,1,-1,-106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
c Sigma0 + pi0  -> K-  +p
c Sigma0 + pi0  -> aK0 +n
        elseif(i1.eq.40.and.iz1.eq.0.and.  !Sigma0
     &     i2.eq.101.and.iz2.eq.0)  then    !pi0
           sig= (STREXCHANGE(e,1,1,-106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
     &        + STREXCHANGE(e,1,-1,-106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2)
     &        * CTParam(72)
c Lambda0 + pi0 -> K-  +p
c Lambda0 + pi0 -> aK0 +n
        elseif(i1.eq.27.and.iz1.eq.0.and.   !aLambda
     &     i2.eq.101.and.iz2.eq.0)  then    !pi0
           sig= (STREXCHANGE(e,1,1,-106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
     &        + STREXCHANGE(e,1,-1,-106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2)
     &        * (1d0-CTParam(72))
c aSigma+ + pi-  -> K+ + ap
c aSigma+ + pi-  -> K0 + an
        elseif(i1.eq.-40.and.iz1.eq.2.and.  !aSigma+
     &     i2.eq.101.and.iz2.eq.-2) then    !pi-
           sig= STREXCHANGE(e,-1,-1,106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
     &        + STREXCHANGE(e,-1,1,106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
c aSigma- + pi+  -> K+ + ap
c aSigma- + pi+  -> K0 +an
        elseif(i1.eq.-40.and.iz1.eq.-2.and. !aSigma-
     &     i2.eq.101.and.iz2.eq.2)  then    !pi+
           sig= STREXCHANGE(e,-1,-1,106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
     &        + STREXCHANGE(e,-1,1,106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
c aSigma0 + pi0  -> K+  +ap
c aSigma0 + pi0  -> K0  +an
        elseif(i1.eq.-40.and.iz1.eq.0.and. !aSigma0
     &     i2.eq.101.and.iz2.eq.0)  then   !pi0
           sig= (STREXCHANGE(e,-1,-1,106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
     &        + STREXCHANGE(e,-1,1,106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2)
     &        * CTParam(72)
c aLambda0 + pi0 -> K+  +ap
c aLambda0 + pi0 -> K0  +an
        elseif(i1.eq.-27.and.iz1.eq.0.and. !aLambda
     &     i2.eq.101.and.iz2.eq.0)  then   !pi0
           sig= (STREXCHANGE(e,-1,-1,106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
     &        + STREXCHANGE(e,-1,1,106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2)   
     &        * (1d0-CTParam(72))
cCHARGED CHANNELS
c Sigma+ + pi0  -> aK0 + p
        elseif(i1.eq.40.and.iz1.eq.2.and.  !Sigma+
     &     i2.eq.101.and.iz2.eq.0)  then   !pi0
           sig= STREXCHANGE(e,1,1,-106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
c Sigma0 + pi+  -> aK0 + p
        elseif(i1.eq.40.and.iz1.eq.0.and.  !Sigma0
     &     i2.eq.101.and.iz2.eq.2)  then   !pi+
           sig= (STREXCHANGE(e,1,1,-106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2)
     &        * CTParam(72)
c Lambda0 + pi+ -> aK0 + p
        elseif(i1.eq.27.and.iz1.eq.0.and.  !Lambda
     &     i2.eq.101.and.iz2.eq.2)  then   !pi+
           sig= (STREXCHANGE(e,1,1,-106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2)
     &        * (1d0-CTParam(72))
c aSigma++ pi0  -> K+ + an
        elseif(i1.eq.-40.and.iz1.eq.2.and. !aSigma+
     &     i2.eq.101.and.iz2.eq.0)  then   !pi0
           sig= STREXCHANGE(e,-1,1,106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
c aSigma0+ pi+  -> K+ + an
        elseif(i1.eq.-40.and.iz1.eq.0.and. !aSigma0
     &     i2.eq.101.and.iz2.eq.2)  then   !pi+
           sig= (STREXCHANGE(e,-1,1,106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2)
     &        * CTParam(72)
c aLambda0+pi+  -> K+ + an
        elseif(i1.eq.-27.and.iz1.eq.0.and. !aLambda
     &     i2.eq.101.and.iz2.eq.2)  then   !pi+
           sig= (STREXCHANGE(e,-1,1,106,1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2)
     &        * (1d0-CTParam(72))
c Sigma- + pi0  -> K- + n
        elseif(i1.eq.40.and.iz1.eq.-2.and. !Sigma-
     &     i2.eq.101.and.iz2.eq.0)  then   !pi0
           sig= STREXCHANGE(e,1,-1,-106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
c Sigma0 + pi-  -> K- + n
        elseif(i1.eq.40.and.iz1.eq.0.and.  !Sigma0
     &     i2.eq.101.and.iz2.eq.-2)  then  !pi-
           sig= (STREXCHANGE(e,1,-1,-106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2)
     &        * CTParam(72)
c Lambda + pi-  -> K- + n
        elseif(i1.eq.27.and.iz1.eq.0.and.  !Lambda
     &     i2.eq.101.and.iz2.eq.-2)  then  !pi-
           sig= (STREXCHANGE(e,1,-1,-106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(-106),massit(1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2)
     &        * (1d0-CTParam(72))
c aSigma-+ pi0  -> K0 + ap
        elseif(i1.eq.-40.and.iz1.eq.-2.and.!aSigma-
     &     i2.eq.101.and.iz2.eq.0)  then   !pi0
           sig= STREXCHANGE(e,-1,-1,106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2
c aSigma0+ pi-  -> K0 + ap
        elseif(i1.eq.-40.and.iz1.eq.0.and. !aSigma0
     &     i2.eq.101.and.iz2.eq.-2)  then  !pi-
           sig= (STREXCHANGE(e,-1,-1,106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2)
     &        * CTParam(72)
c aLambda0+pi-  -> K0 + ap
        elseif(i1.eq.-27.and.iz1.eq.0.and. !aLambda
     &     i2.eq.101.and.iz2.eq.-2)  then  !pi-
           sig= (STREXCHANGE(e,-1,-1,106,-1,i1,iz1,i2,iz2)
     &          *pcms(e,massit(106),massit(-1))**2
     &       *DIVpcms(e,massit(i1),massit(i2))**2)
     &        * (1d0-CTParam(72))
        endif


      else if(io.eq.59) then
c Sigma+ + K- -> Xi0 + pi0
c Sigma+ + K- -> Xi- + pi+
        if(i1.eq.40.and.iz1.eq.2.and.                     !Sigma+
     &     i2.eq.-106.and.iz2.eq.-1) then                 !K-
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,49,1,101,0)
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,101,2)
c Sigma+ + aK0-> Xi0 + pi+
        else if(i1.eq.40.and.iz1.eq.2.and.                !Sigma+
     &     i2.eq.-106.and.iz2.eq.1) then                  !aK0
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,49,1,101,2)
c Lambda0/Sigma0 + aK0 -> Xi0 + pi0
c Lambda0/Sigma0 + aK0 -> Xi- + pi+
        else if((i1.eq.40.or.i1.eq.27).and.iz1.eq.0.and.  !Lambda0/Sigma0
     &     i2.eq.-106.and.iz2.eq.1) then                  !aK0
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,49,1,101,0)
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,101,2)
c Lambda0/Sigma0 + K-  -> Xi0 + pi-
c Lambda0/Sigma0 + K-  -> Xi- + pi0
      else if((i1.eq.40.or.i1.eq.27).and.iz1.eq.0.and.    !Lambda0/Sigma0
     &     i2.eq.-106.and.iz2.eq.-1) then                 !aK-
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,49,1,101,-2)
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,101,0)
c Sigma- + aK0 -> Xi- + pi0
c Sigma- + aK0 -> Xi0 + pi-
        else if(i1.eq.40.and.iz1.eq.-2.and.               !Sigma-
     &     i2.eq.-106.and.iz2.eq.1) then                  !aK0
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,101,0)
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,49,1,101,-2)
c Sigma- + K-  -> Xi- + pi-
        else if(i1.eq.40.and.iz1.eq.-2.and.               !Sigma-
     &     i2.eq.-106.and.iz2.eq.-1) then                 !K-
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,101,-2)
c aSigma+ + K0 -> aXi+ + pi0
c aSigma+ + K0 -> aXi0 + pi+
        else if(i1.eq.-40.and.iz1.eq.2.and.               !aSigma+
     &     i2.eq.106.and.iz2.eq.-1) then                  !K0
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,-49,1,101,0)
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,-49,-1,101,2)
c aSigma+ + K+ -> aXi+ + pi+
        else if(i1.eq.-40.and.iz1.eq.2.and.               !aSigma+
     &     i2.eq.106.and.iz2.eq.1) then                   !K+
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,-49,1,101,2)
c aLambda0/aSigma0 + K0 -> aXi0 + pi0
c aLambda0/aSigma0 + K0 -> aXi+ + pi-
        else if((i1.eq.-40.or.i1.eq.-27).and.iz1.eq.0.and.!aLambda0/aSigma0
     &     i2.eq.106.and.iz2.eq.-1) then                  !K0
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,-49,-1,101,0)
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,-49,1,101,-2)
c aLambda0/aSigma0 + K+ -> aXi0 + pi+
c aLambda0/aSigma0 + K+ -> aXi+ + pi0
        else if((i1.eq.-40.or.i1.eq.-27).and.iz1.eq.0.and.!aLambda0/aSigma0
     &     i2.eq.106.and.iz2.eq.1) then                   !K+
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,-49,-1,101,2)
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,-49,1,101,0)
c aSigma- + K0 -> aXi0 + pi-
        else if(i1.eq.-40.and.iz1.eq.-2.and.              !aSigma-
     &     i2.eq.106.and.iz2.eq.-1) then                  !K0
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,-49,-1,101,-2)
c aSigma- + K+ -> aXi0 + pi0
c aSigma- + K+ -> aXi+ + pi-
        else if(i1.eq.-40.and.iz1.eq.-2.and.              !aSigma-
     &     i2.eq.106.and.iz2.eq.1) then                   !K+
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,-49,-1,101,0)
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,-49,1,101,-2)
        endif




      else if(io.eq.60) then
c Xi0 + pi0 -> K- + Sigma+
c Xi0 + pi0 -> aK0+ Lambda0/Sigma0
        if(i1.eq.49.and.iz1.eq.1.and.        !Xi0
     &     i2.eq.101.and.iz2.eq.0) then      !pi0
          sig =STREXCHANGE(e,40,2,-106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,40,0,-106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,27,0,-106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(27))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c Xi0 + pi- -> K- + Lambda0/Sigma0
c Xi0 + pi- -> aK0+ Sigma-
        elseif(i1.eq.49.and.iz1.eq.1.and.    !Xi0
     &     i2.eq.101.and.iz2.eq.-2) then     !pi-
          sig =STREXCHANGE(e,40,0,-106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,27,0,-106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(27))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,40,-2,-106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c Xi0 + pi+ -> aK0+ Sigma+
        elseif(i1.eq.49.and.iz1.eq.1.and.    !Xi0
     &     i2.eq.101.and.iz2.eq.2) then      !pi+
          sig =STREXCHANGE(e,40,2,-106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c Xi- + pi+ -> aK0 + Lambda0/Sigma0
c Xi- + pi+ -> K-  + Sigma+
        elseif(i1.eq.49.and.iz1.eq.-1.and.   !Xi-
     &     i2.eq.101.and.iz2.eq.2) then      !pi+
          sig =STREXCHANGE(e,40,0,-106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,27,0,-106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(27))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,40,2,-106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c Xi- + pi0 -> aK0 + Sigma-
c Xi- + pi0 -> K-  + Lambda0/Sigma0
        elseif(i1.eq.49.and.iz1.eq.-1.and.   !Xi-
     &     i2.eq.101.and.iz2.eq.0) then      !pi0
          sig =STREXCHANGE(e,40,-2,-106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,40,0,-106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,27,0,-106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(27))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c Xi- + pi- -> K-  + Sigma-
        elseif(i1.eq.49.and.iz1.eq.-1.and.   !Xi-
     &     i2.eq.101.and.iz2.eq.-2) then     !pi-
          sig =STREXCHANGE(e,40,-2,-106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c aXi0 + pi0 -> K+ + aSigma-
c aXi0 + pi0 -> K0 + aLambda0/aSigma0
        elseif(i1.eq.-49.and.iz1.eq.-1.and.  !Xi0
     &     i2.eq.101.and.iz2.eq.0) then      !pi0
          sig =STREXCHANGE(e,-40,-2,106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,-40,0,106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,-27,0,106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-27))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c aXi0 + pi+ -> K+ + aLambda0/aSigma0
c aXi0 + pi+ -> K0 + aSigma+
        elseif(i1.eq.-49.and.iz1.eq.-1.and.  !Xi0
     &     i2.eq.101.and.iz2.eq.2) then      !pi+
          sig =STREXCHANGE(e,-40,0,106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,-27,0,106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-27))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,-40,2,106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c aXi0 + pi- -> K0 + aSigma-
        elseif(i1.eq.-49.and.iz1.eq.-1.and.  !Xi0
     &     i2.eq.101.and.iz2.eq.-2) then     !pi-
          sig =STREXCHANGE(e,-40,-2,106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c aXi+ + pi- -> K0 + aLambda0/aSigma0
c aXi+ + pi- -> K+ + aSigma-
        elseif(i1.eq.-49.and.iz1.eq.1.and.   !aXi+
     &     i2.eq.101.and.iz2.eq.-2) then     !pi-
          sig =STREXCHANGE(e,-40,0,106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,-27,0,106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-27))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,-40,-2,106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c aXi+ + pi0 -> K0 + aSigma+
c aXi+ + pi0 -> K+ + aLambda0/aSigma0
        elseif(i1.eq.-49.and.iz1.eq.1.and.   !aXi+
     &     i2.eq.101.and.iz2.eq.0) then      !pi0
          sig =STREXCHANGE(e,-40,2,106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,-40,0,106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
     &        +STREXCHANGE(e,-27,0,106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-27))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c aXi+   + pi+ -> K+ + aSigma+
        elseif(i1.eq.-49.and.iz1.eq.1.and.   !aXi+
     &     i2.eq.101.and.iz2.eq.2) then      !pi+
          sig =STREXCHANGE(e,-40,2,106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-40))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
        endif



      else if(io.eq.61) then
c Xi0 + aK0 -> Omega- + pi+
        if(i1.eq.49.and.iz1.eq.1.and.
     &     i2.eq.-106.and.iz2.eq.1) then
          sig =STREXCHANGE(e,i1,iz1,i2,iz2,55,0,101,2)
c Xi0 + K- -> Omega- + pi0
        elseif(i1.eq.49.and.iz1.eq.1.and.
     &     i2.eq.-106.and.iz2.eq.-1) then
          sig =STREXCHANGE(e,i1,iz1,i2,iz2,55,0,101,0)
c Xi- + aK0 -> Omega- + pi0
        elseif(i1.eq.49.and.iz1.eq.-1.and.
     &     i2.eq.-106.and.iz2.eq.1) then
          sig =STREXCHANGE(e,i1,iz1,i2,iz2,55,0,101,0)
c Xi- + K- ->  Omega- + pi-
        elseif(i1.eq.49.and.iz1.eq.-1.and.
     &     i2.eq.-106.and.iz2.eq.-1) then
          sig =STREXCHANGE(e,i1,iz1,i2,iz2,55,0,101,-2)
c aXi0 + K0 -> aOmega+ + pi-
        elseif(i1.eq.-49.and.iz1.eq.-1.and.
     &     i2.eq.106.and.iz2.eq.-1) then
          sig =STREXCHANGE(e,i1,iz1,i2,iz2,-55,0,101,-2)
c aXi0 + K+ -> aOmega+ + pi0
        elseif(i1.eq.-49.and.iz1.eq.-1.and.
     &     i2.eq.106.and.iz2.eq.1) then
          sig =STREXCHANGE(e,i1,iz1,i2,iz2,-55,0,101,0)
c aXi+ + K0 -> aOmega+ + pi0
        elseif(i1.eq.-49.and.iz1.eq.1.and.
     &     i2.eq.106.and.iz2.eq.-1) then
          sig =STREXCHANGE(e,i1,iz1,i2,iz2,-55,0,101,0)
c aXi+ + K+ -> aOmega+ + pi+
        elseif(i1.eq.-49.and.iz1.eq.1.and.
     &     i2.eq.106.and.iz2.eq.1) then
          sig =STREXCHANGE(e,i1,iz1,i2,iz2,-55,0,101,2)
        endif



c all i0=62 channels multiplied by 0.5 due to the effect
c of spin degeneracy of the particles in the detailed balance
c calculation
      else if(io.eq.62) then
c Omega- + pi+ -> Xi0 + aK0
        if(i1.eq.55.and.iz1.eq.0.and.
     &     i2.eq.101.and.iz2.eq.2) then
          sig=5d-1*STREXCHANGE(e,49,1,-106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(49))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c Omega- + pi0 -> Xi0 + K-
c Omega- + pi0 -> Xi- + aK0
        elseif(i1.eq.55.and.iz1.eq.0.and.
     &     i2.eq.101.and.iz2.eq.0) then
          sig=5d-1*(STREXCHANGE(e,49,1,-106,-1,i1,iz1,i2,iz2)
     &             +STREXCHANGE(e,49,-1,-106,1,i1,iz1,i2,iz2))
     &        *pcms(e,massit(-106),massit(49))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c Omega- + pi- -> Xi- + K-
        elseif(i1.eq.55.and.iz1.eq.0.and.
     &     i2.eq.101.and.iz2.eq.-2) then
          sig=5d-1*STREXCHANGE(e,49,-1,-106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(-106),massit(49))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c aOmega+ + pi- -> aXi0 + K0
        elseif(i1.eq.-55.and.iz1.eq.0.and.
     &     i2.eq.101.and.iz2.eq.-2) then
          sig=5d-1*STREXCHANGE(e,-49,-1,106,-1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-49))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c aOmega+ + pi0 -> aXi0 + K+
c aOmega+ + pi0 -> aXi+ + K0
        elseif(i1.eq.-55.and.iz1.eq.0.and.
     &     i2.eq.101.and.iz2.eq.0) then
          sig=5d-1*(STREXCHANGE(e,-49,-1,106,1,i1,iz1,i2,iz2)
     &             +STREXCHANGE(e,-49,1,106,-1,i1,iz1,i2,iz2))
     &        *pcms(e,massit(106),massit(-49))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
c aOmega+ + pi+ -> aXi+ + K+
        elseif(i1.eq.-55.and.iz1.eq.0.and.
     &     i2.eq.101.and.iz2.eq.2) then
          sig=5d-1*STREXCHANGE(e,-49,1,106,1,i1,iz1,i2,iz2)
     &        *pcms(e,massit(106),massit(-49))**2
     &     *DIVpcms(e,massit(i1),massit(i2))**2
        end if


      else if(io.eq.63) then
c Lambda0 + Lambda0   -> Xi0 + n
c                     -> Xi- + p
c aLambda0 + aLambda0 -> aXi0 + an
c                     -> aXi+ + ap
        if(abs(i1).eq.27.and.abs(i2).eq.27.and.
     &     i1*i2.gt.0) then
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,49,1,1,-1)
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,1,1)
c Lambda0 + Sigma-  -> Xi- + n
c Lambda0 + Sigma+  -> Xi0 + p
c aLambda0 + aSigma+  -> aXi+ + an
c aLambda0 + aSigma-  -> aXi0 + ap
        elseif( (((abs(i1).eq.27.and.abs(i2).eq.40).and.
     &            (iz1.eq.0.and.abs(iz2).eq.2)).or.
     &           ((abs(i1).eq.40.and.abs(i2).eq.27).and.
     &            (abs(iz1).eq.2.and.iz2.eq.0))).and.
     &           i1*i2.gt.0 ) then
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,49,1,1,1)
c Lambda0 + Sigma0  -> Xi0 + n
c                   -> Xi- + p
c aLambda0 + aSigma0  -> aXi0 + an
c                     -> aXi+ + ap
        elseif( ((abs(i1).eq.27.and.abs(i2).eq.40).or.
     &           (abs(i1).eq.40.and.abs(i2).eq.27)).and.
     &           i1*i2.gt.0.and.iz1.eq.0.and.iz2.eq.0) then
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,49,1,1,-1)
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,1,1)
c Sigma0  + Sigma0  -> Xi0 + n
c                   -> Xi- + p
c aSigma0  + aSigma0  -> aXi0 + an
c                     -> aXi+ + ap
        elseif(abs(i1).eq.40.and.abs(i2).eq.40.and.
     &         iz1.eq.0.and.iz2.eq.0.and.i1*i2.gt.0) then
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,49,1,1,-1)
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,1,1)
c Sigma0  + Sigma+  -> Xi0 + p
c Sigma0  + Sigma-  -> Xi- + n
c aSigma0  + aSigma-  -> aXi0 + ap
c aSigma0  + aSigma+  -> aXi+ + an
        elseif(abs(i1).eq.40.and.abs(i2).eq.40.and.
     &         ((iz1.eq.0.and.abs(iz2).eq.2).or.
     &          (abs(iz1).eq.2.and.iz2.eq.0)).and.
     &         i1*i2.gt.0) then
          sig=4d0*STREXCHANGE(e,40,0,40,0,49,-1,1,1)
c Sigma+  + Sigma-  -> Xi0 + n
c                   -> Xi- + p
c aSigma-  + aSigma+  -> aXi0 + an
c                     -> aXi+ + ap
        elseif(abs(i1).eq.40.and.abs(i2).eq.40.and.
     &         abs(iz1).eq.2.and.abs(iz2).eq.2.and.
     &         i1*i2.gt.0.and.(iz1+iz2).eq.0) then
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,49,0,1,-1)
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,49,-1,1,1)
        endif






      elseif(io.eq.64) then
c Xi0 + n -> Lambda0 + Lambda0
c Xi0 + n -> Lambda0 + Sigma0
c Xi0 + n -> Sigma0  + Sigma0
c Xi0 + n -> Sigma+  + Sigma-
c Xi- + p -> Lambda0 + Lambda0
c Xi- + p -> Lambda0 + Sigma0
c Xi- + p -> Sigma0  + Sigma0
c Xi- + p -> Sigma+  + Sigma-
        if( (i1.eq.49.and.iz1.eq.1.and.i2.eq.1.and.iz2.eq.-1).or.
     &      (i1.eq.1.and.iz1.eq.-1.and.i2.eq.49.and.iz2.eq.1).or.
     &      (i1.eq.49.and.iz1.eq.-1.and.i2.eq.1.and.iz2.eq.1).or.
     &      (i1.eq.1.and.iz1.eq.1.and.i2.eq.49.and.iz2.eq.-1)) then
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,27,0,27,0) !Lam + Lam
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,27,0,40,0) !Lam + Sig0
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,40,0,40,0) !Sig0 + Sig0
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,40,2,40,-2) !Sig+ + Sig-

c Xi0 + p -> Lambda0 + Sigma+
c Xi0 + p -> Sigma0  + Sigma+
        elseif((i1.eq.49.and.iz1.eq.1.and.i2.eq.1.and.iz2.eq.1).or.
     &         (i1.eq.1.and.iz1.eq.1.and.i2.eq.49.and.iz2.eq.1)) then
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,27,0,40,2) !Lam + Sig+
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,40,0,40,2) !Sig0+ Sig+
c Xi- + n -> Lambda0 + Sigma-
c Xi- + n -> Sigma0  + Sigma-
        elseif((i1.eq.49.and.iz1.eq.-1.and.i2.eq.1.and.iz2.eq.-1).or.
     &         (i1.eq.1.and.iz1.eq.-1.and.i2.eq.49.and.iz2.eq.-1)) then
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,27,0,40,-2) !Lam + Sig-
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,40,0,40,-2) !Sig0+ Sig-
c aXi0 + an -> aLambda0 + aLambda0
c aXi0 + an -> aLambda0 + aSigma0
c aXi0 + an -> aSigma0  + aSigma0
c aXi0 + an -> aSigma+  + aSigma-
c aXi- + ap -> aLambda0 + aLambda0
c aXi- + ap -> aLambda0 + aSigma0
c aXi- + ap -> aSigma0  + aSigma0
c aXi- + ap -> aSigma+  + aSigma-
        elseif((i1.eq.-49.and.iz1.eq.-1.and.i2.eq.-1.and.iz2.eq.1).or.
     &         (i1.eq.-1.and.iz1.eq.1.and.i2.eq.-49.and.iz2.eq.-1).or.
     &         (i1.eq.-49.and.iz1.eq.1.and.i2.eq.-1.and.iz2.eq.-1).or.
     &         (i1.eq.-1.and.iz1.eq.-1.and.i2.eq.-49.and.iz2.eq.1)) then
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,-27,0,-27,0) !aLam + aLam
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,-27,0,-40,0) !aLam + aSig0
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,-40,0,-40,0) !aSig0 + aSig0
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,-40,-2,-40,2) !aSig- + aSig+

c aXi0 + ap -> aLambda0 + aSigma-
c aXi0 + ap -> aSigma0  + aSigma-
        elseif((i1.eq.-49.and.iz1.eq.-1.and.i2.eq.-1.and.iz2.eq.-1).or.
     &         (i1.eq.-1.and.iz1.eq.-1.and.i2.eq.-49.and.iz2.eq.-1))then
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,-27,0,-40,-2) !aLam + aSig-
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,-40,0,-40,-2) !aSig0+ aSig-

c aXi+ + an -> aLambda0 + aSigma+
c aXi+ + an -> aSigma0  + aSigma+
        elseif((i1.eq.-49.and.iz1.eq.1.and.i2.eq.-1.and.iz2.eq.1).or.
     &         (i1.eq.-1.and.iz1.eq.1.and.i2.eq.-49.and.iz2.eq.1)) then
          sig=STREXCHANGE(e,i1,iz1,i2,iz2,-27,0,-40,2) !aLam + aSig+
     &       +STREXCHANGE(e,i1,iz1,i2,iz2,-40,0,-40,2) !aSig0+ aSig+
        endif


      else
        write(*,*) "Error: STREXWrapper called with io: ",io
        stop 137
      endif
      STREXWrapper=sig
      return
      end function
