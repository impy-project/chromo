c $Id: anndec.f,v 1.20 2003/05/02 13:14:47 weber Exp $    
C####C##1#########2#########3#########4#########5#########6#########7##
        subroutine anndec(io,mm1,ii1,iiz1,mm2,ii2,iiz2,sqrts,sig,gam)
c
c
cinput io    : 0: annihilation; 1: decay
cinput mm1   : mass of scattering/decaying particle 1
cinput ii1   : ID of  scattering/decaying particle 1
cinput iiz1  : $2\cdot I_3$ of scattering/decaying particle 1
cinput mm2   : mass of scattering particle 2
cinput ii2   : ID of  scattering particle 2
cinput iiz2  : $2\cdot I_3$ of scattering particle 2
cinput sqrts : $sqrts{s}$ of collision; resonance mass for decay
coutput sig  : resonance scattering cross section
coutput gam  : width of the produced resonance
c
c     {\tt anndec} handles meson-baryon and meson-meson annihilations 
c     as well as all meson and baryon resonance decays. In case of 
c     annihilations it returs the total resonance production cross section 
c     and the decay width of the resonance chosen as final state.
c     The final state itself for both cases, annihilation and decay 
c     is returned via the {\tt newpart} common block. In the case
c     of a decay the final state may consist of up to 4 particles.
c
c     Technically, {\tt anndec} is only an interface to {\tt anndex},
c     which actually handles the summation of the Breit-Wigner formulas
c     in the annihilation case and the final state generation for 
c     annihilation and decay.
c
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      implicit none
      integer io,i1,i2,iz1,iz2,ii1,ii2,iiz1,iiz2,is
      real*8 m1,m2,mm1,mm2,sig,gam,sqrts
      
      include 'comres.f'
      include 'options.f'

      integer strit
C
C    ************************************************************************
C  Case 1 :  Two ingoing Particles --> One outgoing Particle (Resonance,...)
C
C
      i1=ii1
      i2=ii2
      iz1=iiz1
      iz2=iiz2
      m1=mm1
      m2=mm2

      if(io.eq.0)then !annihilation
C
      sig=0d0
      gam=0d0
C
C     Check if (sqrt(s)-masses of ingoing particles) is significant different
C     from zero
C
      if(sqrts-mm1-mm2.le.1d-3)return
C       
C     Check if CTOption(15) is set different from zero-->then skip anndec.f
C
      if(CTOption(15).ne.0)return
C
C
C     Check if itype of particle one is smaller than the one of particle two
C     if so --> interchange particle one and particle two
C      in case of particle one = B and particle two = M --> then
C      new particle one = M and new particle two = B
C
      if(iabs(i1).lt.iabs(i2))call swpizm(i1,iz1,m1,i2,iz2,m2)
C
C     Determination of the amount of the netstrangness
C
      is=iabs(strit(i1)+strit(i2))
C
C       maxbar (Maximum Baryon ityp)
C       if second particle is antibaryon and first particle is strange
C       switch antibaryon to baryon
C

csab: obsolete
c      if(iabs(i2).le.maxbar)then
c         if(i2.lt.0)then 
c           if(strit(i1).ne.0)i1=-i1 ! get corresponding anti-branch
c         end if
c      end if 
C
C
C     Check if both particles are mesons 
C
      if(iabs(i1).ge.minmes.and.iabs(i2).ge.minmes)then
C
C
c... boson+boson sector
C
C
C     Check if amount of netstrangeness is greater than 1
c     currently no resonant processes for |s|>1 are implemented

         if(is.gt.1)return 

         if(is.ne.0)then
           call anndex(0,m1,i1,iz1,m2,i2,iz2,sqrts,
     .        sig,gam,maxbrm,minmes+1,maxmes,bmtype,branmes)
         else
           call anndex(0,m1,i1,iz1,m2,i2,iz2,sqrts,
     .        sig,gam,maxbrm,minmes+1,maxmes,bmtype,branmes)
         endif

C
C        Check if second particle is baryon
C        (with zero amount of netstrangeness) ?
C        e.g. pion-nucleon case
C 
      else if(is.eq.0.and.iabs(i2).le.maxbar)then
c... (anti-)N*,D*
         call anndex(0,m1,i1,iz1,m2,i2,iz2,sqrts,sig,gam,
     .        maxbra,minnuc+1,maxdel,brtype,branres)

C
C       Check if second particle is baryon
C        (with amount one of netstrangeness) ?
C
      else if(is.eq.1.and.iabs(i2).le.maxbar)then
c... (anti-)Y*
         call anndex(0,m1,i1,iz1,m2,i2,iz2,sqrts,sig,gam,
     .        maxbrs1,minlam+1,maxsig,bs1type,branbs1)        
C
C       Check if second particle is baryon
C        (with amount two of netstrangeness) ?
C
      else if(is.eq.2.and.iabs(i2).le.maxbar)then
c... (anti-)X*
         call anndex(0,m1,i1,iz1,m2,i2,iz2,sqrts,sig,gam,
     .        maxbrs2,mincas+1,maxcas,bs2type,branbs2)
C
C
      else 
        sig=0d0
        return
      end if
C
C    *************************************************************
Css  Case 2 : one ingoing particle (resonance,..) --> 2-4 outgoing
Css  particles (decay)
C
      else ! decay !!!!!

         i2=0
         iz2=0
         m2=0.d0
         is=iabs(strit(i1))
c
         if(iabs(i1).ge.minmes)then ! meson dec. 
            call anndex(1,m1,i1,iz1,m2,i2,iz2,sqrts,sig,gam,
     .           maxbrm ,minmes+1,maxmes,bmtype,branmes)


         else if(is.eq.0)then   ! n*,d,d*
            call anndex(1,m1,i1,iz1,m2,i2,iz2,sqrts,sig,gam,
     .           maxbra,minnuc+1,maxdel,brtype,branres)


         else if(is.eq.1)then   ! 
            call anndex(1,m1,i1,iz1,m2,i2,iz2,sqrts,sig,gam,
     .            maxbrs1,minlam+1,maxsig,bs1type,branbs1)


         else if(is.eq.2)then
            call anndex(1,m1,i1,iz1,m2,i2,iz2,sqrts,sig,gam,
     .           maxbrs2,mincas+1,maxcas,bs2type,branbs2)

         else
            write(6,*)'make22(anndex): s=',is,'not included'
            stop 137
         end if
C
C    End of Cases 1 and 2 : annihilation/decay
C
      end if
C    ************************************************************
C     
      return
      end
      


C####C##1#########2#########3#########4#########5#########6#########7##
       subroutine anndex(io,m1,i1,iiz1,m2,i2,iiz2,sqrts,sig,gam,
     &            maxbr,mini,maxi,btype,branch)
c
cinput io     : 0: annihilation; 1: decay
cinput m1     : mass of scattering/decaying particle 1
cinput i1     : ID of  scattering/decaying particle 1
cinput iz1    : $2\cdot I_3$ of scattering/decaying particle 1
cinput m2     : mass of scattering particle 2
cinput i2     : ID of  scattering particle 2
cinput iz2    : $2\cdot I_3$ of scattering particle 2
cinput sqrts  : $sqrts{s}$ of collision; resonance mass for decay
coutput sig   : resonance scattering cross section
coutput gam   : width of the produced resonance
cinput maxbr  : number of decay channels for particle class
cinput mini   : smallest {\tt ityp} of particle class
cinput maxi   : largest {\tt ityp} of particle class
cinput btype  : array with exit channel definitions
cinput branch : array with branching ratios for final state
c
c
c     {\tt anndex} performs meson-baryon and meson-meson annihilations 
c     as well as all meson and baryon resonance decays. In case of 
c     annihilations it returs the total resonance production cross section 
c     and the decay width of the resonance chosen as final state.
c     The final state itself for both cases, annihilation and decay 
c     is returned via the {\tt newpart} common block. In the case
c     of a decay the final state may consist of up to 4 particles.
c
c     In {\tt anndex} the actual summation over Breit-Wigner formulas
c     in the case of annihilations is performed. 
c
c     For decays the branch is choosen according to the mass dependent
c     part of the decay width (call to {\tt fbrancx}); then the final
c     state (which consists of particle {\tt ityp}, $2\cdot I_3$ and
c     mass is generated and transferred to teh {\tt newpart} common
c     block.
c
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'comres.f'
      include 'comwid.f'
      include 'newpart.f'
      include 'options.f'

      real*8 pi,cc,sqrts
      parameter(pi=3.1415927,cc=0.38937966)
      integer maxbr,mini,maxi,btype(4,0:maxbr)
      real*8 branch(0:maxbr,mini:maxi)
      integer icnt,is


      integer io,i,j,i1,i2,iz1,iz2,itag,ii1,ii2
      integer itn1,itnz1,iiz1,iiz2
      real*8 m1,m2,prob(0:100),sum,sig,gam,cgk2
      real*8 sigi(minnuc:maxmes),mmax,mmin,br,mmi1,mmi2,ppcm,gt
      real*8 m,g,mo
        
c      functions
      real*8 fbrwig,pcms,massit
      real*8 mminit,widit,fbrancx,ranf,fcgk,fwidth
      real*8 fprwdt
      integer jit,isoit,strit,chrmit


csab: needed for antibaryon-meson scattering symmetry fix:
      iz1=iiz1
      iz2=iiz2

      if(io.eq.1)then
C
C
C   one ingoing particle --> two,three,four outgoing particles
C
c... decays
  
         do 3 i=0,maxbr
            if(isoit(btype(1,i))+isoit(btype(2,i))+isoit(btype(3,i))+
     &         isoit(btype(4,i)).lt.iabs(iz1).or.
     &           m1.lt.mminit(btype(1,i))+mminit(btype(2,i))
     &                +mminit(btype(3,i))+mminit(btype(4,i)) )then
               prob(i)=0.d0
            else
               prob(i)=fbrancx(i,iabs(i1),iz1,m1,branch(i,iabs(i1)),
     &              btype(1,i),btype(2,i),btype(3,i),btype(4,i))
            endif
 3       continue

         icnt=0
c... find out branch = i
 9       call getbran(prob,0,100,sum,0,maxbr,i)

         if(i.gt.maxbr)then
            write(6,*)'anndex(dec): no final state found for:',i1,m1,iz1
            write(6,*)'please check minimal masses: m1,m1min,m2min'
            write(6,*)'and iso3 of decaying particle'
            write(6,*)(prob(j),j=0,maxbr)
            stop 137
         end if

c... get itypes and set number of outgoing particles, prepare final state
         nexit=2
         itypnew(1)=btype(1,i)
         itot(1)=isoit(itypnew(1))
         itypnew(2)=btype(2,i)
         itot(2)=isoit(itypnew(2))


         itypnew(3)=btype(3,i)
         if(itypnew(3).ne.0) then
            itot(3)=isoit(itypnew(3))
            pnew(5,3)=massit(itypnew(3))
            nexit=nexit+1
         endif
         itypnew(4)=btype(4,i)
         if(itypnew(4).ne.0) then
            itot(4)=isoit(itypnew(4))
            pnew(5,4)=massit(itypnew(4))
            nexit=nexit+1
         endif


c check for some special cases involving decay of antibaryons and 
c strange mesons
         if(iabs(i1).ge.minmes.and.strit(i1).ne.0) then
            do 41 j=1,nexit
               if(strit(itypnew(j)).ne.0)then
c  for anti-K* decays(mesons with one s-quark)
                  itypnew(j)=isign(itypnew(j),i1) 
               end if
 41         continue
c charm mesons
         elseif(iabs(i1).ge.minmes.and.chrmit(i1).ne.0) then
            do 441 j=1,nexit
               if(chrmit(itypnew(j)).ne.0)then
c  for anti-D* decays(mesons with one c-quark)
                  itypnew(j)=isign(itypnew(j),i1) 
               end if
 441         continue
         elseif(iabs(i1).lt.minmes) then
c     the (anti-)baryon MUST always be the first outgoing particle
c     -> conserve baryon-charge
            itypnew(1)=isign(itypnew(1),i1)
            do 42 j=2,nexit
               if(strit(itypnew(j)).ne.0.and.i1.lt.0) then
                  itypnew(j)=(-1)*itypnew(j)
               endif
 42         continue
         endif



c... get isopin-3 components
         itag=-50
         call isonew4(isoit(i1),iz1,itot,i3new,itag)
c     write(6,*)'anndem:',iz3,iz1,iz2,'#',is3,is1,is2

        
c...  get masses
         if(widit(itypnew(1)).ge.1.d-4.and.
     &        widit(itypnew(2)).le.1.d-4)then
c...  i1 is a broad meson
            pnew(5,2)=massit(itypnew(2))
            mmin=mminit(itypnew(1))
            mo = pnew(5,2)
            if(nexit.gt.2) then
               do 39 j=3,nexit
                  mo=mo+pnew(5,j)
 39            continue
            endif

            mmax=sqrts-mo
            call getmas(massit(itypnew(1)),widit(itypnew(1)),itypnew(1)
     &                  ,isoit(itypnew(1)),mmin,mmax,mo,pnew(5,1))
           
         elseif(widit(itypnew(2)).ge.1.d-4
     &           .and.widit(itypnew(1)).le.1.d-4)then
c...  i2 is a broad meson
            pnew(5,1)=massit(itypnew(1))
            mmin=mminit(itypnew(2))
          
            mo = pnew(5,1)
            if(nexit.gt.2) then
               do 49 j=3,nexit
                  mo=mo+pnew(5,j)
 49            continue
            endif
            mmax=sqrts-mo
            call getmas(massit(itypnew(2)),widit(itypnew(2)),itypnew(2)
     &           ,isoit(itypnew(2)),mmin,mmax,mo,pnew(5,2))

         elseif(widit(itypnew(1)).ge.1.d-4
     &           .and.widit(itypnew(2)).ge.1.d-4)then
c...  i1&i2 are both broad 
            if(ranf(0).gt.0.5)then
               mmin=mminit(itypnew(1))
               mo=mminit(itypnew(2))
               if(nexit.gt.2) then
                  do 59 j=3,nexit
                     mo=mo+pnew(5,j)
 59               continue
               endif
               mmax=sqrts-mo

               call getmas(massit(itypnew(1)),widit(itypnew(1)),
     &              itypnew(1),isoit(itypnew(1)),mmin,mmax,mo,pnew(5,1))

               mmin=mminit(itypnew(2))
               mo=pnew(5,1)
               if(nexit.gt.2) then
                  do 69 j=3,nexit
                     mo=mo+pnew(5,j)
 69               continue
               endif
               mmax=sqrts-mo
               call getmas(massit(itypnew(2)),widit(itypnew(2)),
     &              itypnew(2),isoit(itypnew(2)),mmin,mmax,mo,pnew(5,2))

            else ! of ranf.gt.0.5
               mmin=mminit(itypnew(2))
               mo=mminit(itypnew(1))
               if(nexit.gt.2) then
                  do 79 j=3,nexit
                     mo=mo+pnew(5,j)
 79               continue
               endif
               mmax=sqrts-mo
               call getmas(massit(itypnew(2)),widit(itypnew(2)),
     &              itypnew(2),isoit(itypnew(2)),mmin,mmax,mo,pnew(5,2))

               mmin=mminit(itypnew(1))
               mo=pnew(5,2)
               if(nexit.gt.2) then
                  do 89 j=3,nexit
                     mo=mo+pnew(5,j)
 89               continue
               endif
               mmax=sqrts-mo
               call getmas(massit(itypnew(1)),widit(itypnew(1)),
     &              itypnew(1),isoit(itypnew(1)),mmin,mmax,mo,pnew(5,1))
 
            endif
c     none are broad
         else 
            pnew(5,2)=massit(itypnew(2))
            pnew(5,1)=massit(itypnew(1))
         end if

         mmax=0.d0
         do 99 j=1,nexit
            mmax=mmax+pnew(5,j)
 99      continue

         if(sqrts.le.mmax)then
            write(6,*)' *** error(anndex): treshold violated',sqrts-mmax
            stop 137
         end if 
C     
C     
C   two ingoing particles --> one outgoing particle (resonance)
C
C     i.e. (i0=0)       
      else
c.... collisions: find in-branch = j  
         sig=0.0
         gam=0.0         
C
            ii1=i1
            ii2=i2

csab: anti-baryon meson scattering symmetry (baryon always resides in 2nd slot):
         if(iabs(i2).le.99.and.i2.lt.0) then
            iz2=-1*iz2
            iz1=-1*iz1
            if(strit(i1).ne.0) ii1=-1*i1
         endif


c  for strange - nonstrange meson-meson scattering: strip sign
         is=iabs(strit(i1)+strit(i2))
         if(is.ne.0.and.iabs(i1).ge.minmes.and.iabs(i2).ge.minmes) then
            ii1=iabs(i1)
            ii2=iabs(i2)
         endif
c  for meson baryon, strip sign of baryon   
         if(iabs(i2).le.maxbar) then
            ii2=iabs(i2)
         endif

c     
         call getobr(btype,0,maxbr,ii1,ii2,j)

         if(j.eq.-99)return

         mmi1=mminit(i1)
         mmi2=mminit(i2)
c
C   next line outside the loop (compare post-QM-version: inside the loop)
C
C   
         ppcm=pcms(sqrts,m1,m2)
C   
C      Loop over different branches (resonances...)
C  
         do 88 i=mini,maxi
            sigi(i)=0.d0
            br=branch(j,i)
            gt=widit(i)
            if(br*gt.lt.1d-4)goto 88
            
            cgk2=fcgk(i1,i2,iz1,iz2,i)
            if(br*cgk2.gt.0d0.and.sqrts.gt.mmi1+mmi2+1d-2.and.
     &          ppcm.gt.1d-2)then 
C
C
               br=fprwdt(j,i,iz1+iz2,sqrts)/fwidth(i,iz1+iz2,sqrts)
               m=dabs(sqrts)
               g=fwidth(i,iz1+iz2,m)
               sigi(i)=dble(jit(i)+1)
     /                 /dble((jit(i1)+1)*(jit(i2)+1))
     *           *pi/ppcm**2*br
     *           *g*g/((m-massit(i))**2+g*g/4d0)*cgk2*cc
               end if
              if(sigi(i).gt.1e10)then
                write(6,*)' ***error(anndec) cross section too high '
                write(6,*)'anndex(ann):',i,
     ,           br,cgk2,fbrwig(i,iz1+iz2,sqrts,1),
     ,               1/pcms(sqrts,m1,m2),sigi(i)
                write(6,*)m1,m2,sqrts
                write(6,*)i1,i2,i   
                write(6,*)iz1,iz2,iz1+iz2
              end if
C
C
 88      continue
c...  find outgoing resonance      
 108     call getbran(sigi,minnuc,maxmes,sig,mini,maxi,itn1)

csab: iso3 conservation:
            itnz1=iiz1+iiz2

         if(sig.ge.1d-10)then
            gam=fwidth(itn1,itnz1,sqrts)
         end if
         
c     copy created resonance into newpart arrays
         itypnew(1)=itn1
         i3new(1)=itnz1
         pnew(5,1)=sqrts

         if(iabs(i1).ge.minmes.and.iabs(i2).ge.minmes) then
            if(iabs(strit(i1)+strit(i2)).ne.0) then
               itypnew(1)=isign(itypnew(1),i1*i2)
            endif
         else
            itypnew(1)=isign(itypnew(1),i2)
         endif

 3333    continue 


C
C   End of two Cases (annihilation/decay)
C
C
      end if                    !dec/ann

      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function fbrwig(i,iz,mi,bit) 
c
cinput  i  : resonance ID  
cinput  iz : $2\cdot I_3$ of resonance 
cinput  mi : mass of resonance 
cinput bit : sign is used as option to toggle between fixed and m.dep. widths
c
c  {\tt fbrwig} returns a normalized Breit-Wigner Function.
c  Note, that the normalization actually only holds true for
c  fixed decay widths. {\tt fbrwig}, however, uses per default
c  a mass dependent width. You should divide by
c  {\tt bwnorm} to obtain normalized Breit-Wigners for mass dependent 
c  widths also. For {\tt bit} < 0 a fixed width is used.
c
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      implicit none

      integer i,iz,bit
      real*8 pi,mi,g,fwidth,massit,widit
      real*8 f,e,m0,g1,g2
      parameter(pi=3.1415927)

      f(e,m0,g1,g2)=0.5/pi*g1/((e-m0)**2+0.25*g2**2)

      if(bit.lt.0)then
         g=widit(i)
         fbrwig=f(mi,massit(i),g,g)
      else
         g=fwidth(i,iz,mi)
         fbrwig=f(mi,massit(i),g,g)
      end if

      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
        subroutine getbran(x,dmin,dmax,sumx,nmin,nmax,i)
c
c Unit : Infrastructure
c Author :     L.A. Winckelmann, S.A. Bass
c
cinput   x : vector containing weights, dimension is {\tt x(dmin:dmax)}
cinput  dmin : lower dimension of {\tt x}
cinput  dmax : upper dimension of {\tt x}
coutput sumx : sum of elements of {\tt x} from {\tt nmin} to {\tt nmax}
cinput  nmin : lower boundary for {\tt getbran} operation
cinput  nmax : upper boundary for {\tt getbran} operation
coutput i : index of element which has been choosen randomly
c
c     {\tt getbran} takes a vector of weights or probabilities
c     {\tt x(dmin:dmax)} and sums up the elements from
c     {\tt nmin} to {\tt nmax}. It then chooses randomly an element {\tt i}
c     between {\tt nmin} and {\tt nmax}. The probability of
c     choosing {\tt i} depends on the weights contained in {\tt x}.
c
c     \danger{ {\tt i} will be undefined if {\tt sum} is less or
c     equal to zero}
c
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
        implicit none
        integer i,j,nmin,nmax,dmin,dmax,itrial
        real*8 x(dmin:dmax),sumx,ranf,rx,rs,cut
        parameter(cut=1d-20)

        sumx=0.d0
        itrial=0
        do 27 j=nmin,nmax
          sumx=sumx+x(j)
 27     continue
        if(sumx.lt.cut) then
           i=nmax+1
           return
        endif

 28     continue
        itrial=itrial+1
        rs=ranf(0)
        rx=sumx*rs
        do 108 j=nmin,nmax
          if(rx.le.x(j))then
            i=j
            return
          end if
          rx=rx-x(j)
 108   continue

       if(itrial.lt.10) goto 28

       do j=nmin,nmax
          if(x(j).gt.0d0) i=j
       enddo

        write(6,*)'warning: getbran x(',nmin,',',nmax,'):'
        write(6,1008)(x(j),j=nmin,nmax)
        write(6,*)'dmin,dmax,sumx,nmin,nmax,i'
        write(6,*)dmin,dmax,sumx,nmin,nmax,i
        write(6,*) 'rs,sumx,rx ',rs,sumx,rx
        write(6,*) '-> i = max(x(j).gt.0) forced'
        return
c        stop 137
 1008   format(5e10.4)
        end



C####C##1#########2#########3#########4#########5#########6#########7##
        subroutine getobr(x,dmin,dmax,i1,i2,i)
c
c
cinput x : array, either {\tt brtype, bmtype, bs1type} or {\tt bs2type}
cinput dmin : lower dimension of {\tt x(4,dmin:dmax)}
cinput dmin : upper dimension of {\tt x(4,dmin:dmax)}
cinput i1 : ID of first incoming particle
cinput i2 : ID of second incoming particle
coutput i : index of decay branch into {\tt i1} and {\tt i2}
c
c     {\tt getobr} returns the index of the decay branch for the 
c     exit channel $B^* \rightarrow$ {\tt i1} + {\tt i2} 
c     from one of the arrays
c     {\tt brtype, bmtype, bs1type} or {\tt bs2type}. This index
c     is needed for the calculation of the cross section
c     {\tt i1} + {\tt i2} $\rightarrow B^*$. 
c
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
        implicit none
        integer i,j,i1,i2,dmin,dmax,x(4,dmin:dmax)

        do 108 j=dmin,dmax
          if((x(1,j).eq.i1.and.x(2,j).eq.i2.and.x(3,j).eq.0).OR.
     &       (x(1,j).eq.i2.and.x(2,j).eq.i1.and.x(3,j).eq.0))then
            i=j
            return
          end if
  108    continue  
         i=-99        
         return
         end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine normit (sigma,isigline)
c
c     Revision : 1.0
c
cinput sigma : vector with all (partial) cross sections
coutput sigma : unitarized vector with cross sections
cinput isigline : process class of cross sections
c
c     {\tt normit} unitarizes the cross sections contained in the
c     {\tt sigma} array. The total cross section is stored in 
c     {\tt sigma(0)}. The partial cross sections are unitarized
c     (rescaled) such, that their sum adds up to the total cross
c     section. Confidence levels can be assigned to different
c     partial cross sections indicating whether they may be rescaled
c     (i.e. if they are not well known) or whether they must not
c     be rescaled (i.e. because they have been fitted to experimental
c     data).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
        
      include 'options.f'
      include 'comres.f'
      
      real*8 sigma(0:maxpsig)
      integer isigline

      integer i, npsig, restart
      integer uncert(1:maxpsig)
      real*8 diff, sumpart, gsumpart
      real*8 newsig(0:maxpsig)

c get the number of channels
         npsig=sigmaln(1,1,isigline)
c normalize only if sigtot is not given by the sum of sigpart        
      if (sigmaln(2,1,isigline).gt.0) then
c copy array 
         do 10 i=1,npsig
            uncert(i)=sigmaln(i+2,2,isigline)
 10      continue
 100     restart=0
         sumpart=0
         gsumpart=0
c calculate the sum of all sigpart
         do 20 i=1,npsig

            sumpart=sumpart+sigma(i)
            gsumpart=gsumpart+sigma(i)*uncert(i)
 20      continue
c difference between sigtot and the sum of sigpart
         diff=sigma(0)-sumpart
c if all channels are exactly zero, there must be an error in blockres.f!
         if (sumpart.eq.0.0) then
          write (6,*) 'normit: Error! sumpart.eq.0 at sigline ',isigline
c            stop 137
            return
        endif
         if (gsumpart.eq.0.0) then
            do 50 i=1,npsig
c now all channels can be modified
               if (uncert(i).eq.0) then
c                 write (6,*) 'modify channel',i
                  uncert(i)=1
               endif
 50         continue
c restart calculation
            goto 100
         endif
         do 60 i=1,npsig
c rescale channels
            newsig(i)=sigma(i)+uncert(i)*diff*sigma(i)/gsumpart
c if a channel is negative...
            if (newsig(i).lt.0) then
c set it to zero and restart
               sigma(i)=0.0
               restart=1
            endif
 60      continue
         if (restart.eq.1) goto 100
c copy new values to sigma
         do 70 i=1,npsig
            sigma(i)=newsig(i)
 70      continue
      endif

chp cto 7 is now handled in scatter      
c      if (CTOption(7).eq.1.and.sigma(2).gt.1d-10) then
c         sigma(1)=0d0
c      endif
chp 
      if (CTOption(7).eq.-1) then
         do 80 i=2,npsig
            sigma(i)=0d0
 80      continue
      end if
      
      return
      end



C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function fwidth(ir,izr,m)
c
cinput  ir  : resonance ID  
cinput  izr : $2\cdot I_3$ of resonance 
cinput  m   : mass of resonance
c
c     {\tt fwidth} returns the mass-dependent total decay width
c     of the resonance {\tt ir}.
c
c
C####C##1#########2#########3#########4#########5#########6#########7##

      implicit none

      include 'comres.f'
      include 'comwid.f'
      include 'options.f'

      integer i,ir,izr,mm,mp,ires
      real*8 gtot,m,widit,splint
      real*8 minwid, fprwdt

      if (CTOption(1).ne.0) then
         fwidth=widit(ir)
         return
      endif
      if (wtabflg.gt.0.and.CTOption(33).eq.0) then
         ires=iabs(ir)
         minwid=min(widit(ir),1D-8)
         if (ires.ge.minbar.and.ires.le.maxbar) then       !baryons
c widths are continued horicontally outside the spline region
                if(m.le.maxtab2)then
              fwidth=max(splint(tabx(1),fbtaby(1,ires,1),
     .           fbtaby(1,ires,2),widnsp,m),minwid)
                else
              fwidth=max(splint(tabx(1),fbtaby(1,ires,1),
     .           fbtaby(1,ires,2),widnsp,maxtab2),minwid)
                endif
         else if (ires.ge.minmes.and.ires.le.maxmes) then  !mesons
c widths are continued horicontally outside the spline region
                if(m.le.maxtab2)then
              fwidth=max(splint(tabx(1),fmtaby(1,ires,1),
     .           fmtaby(1,ires,2),widnsp,m),minwid)
                else
              fwidth=max(splint(tabx(1),fmtaby(1,ires,1),
     .           fmtaby(1,ires,2),widnsp,maxtab2),minwid)
                endif
         else
            write (6,*) '*** error(fwidth) wrong itype:',ir
            fwidth=0
         endif
      else
         call brange(ir,mm,mp)
         gtot=0d0
         if(mp.gt.0)then
            do 27 i=mm,mp
               gtot=gtot+fprwdt(i,ir,izr,m)
 27         continue
         end if
         fwidth=gtot            !*widit(ir) 
      end if
      
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function fprwdt(i,ir,izr,mi)
c
cinput  i   : decay branch
cinput  ir  : resonance ID  
cinput  izr : $2\cdot I_3$ of resonance 
cinput  mi  : mass of resonance
c
c     {\tt fprwdt} returns the mass dependent partial decay width
c     of the decay channel {\tt i}. 
c
C####C##1#########2#########3#########4#########5#########6#########7##

      implicit none
      real*8 m,br,bi,mir,g,mi
      integer i,ir,izr,i1,i2,i3,i4
      real*8 widit,fbrancx,mminit,massit

      call b3type(ir,i,bi,i1,i2,i3,i4)
      m=dabs(mi)
      g=0d0
      mir=massit(ir)
      if(bi.gt.1d-9.and.mir.gt.mminit(i1)+mminit(i2))then
         br=fbrancx(i,ir,izr,m,bi,i1,i2,i3,i4)
c     write(6,*)'   ',bi,gi,widit(ir)
         g=br*widit(ir) 
      end if
      fprwdt=g
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       real*8 function fbrancx(i,ir,izr,em,bi,b1,b2,b3,b4)
c 
cinput i  : decay branch
cinput ir : ID of resonance
cinput em : actual mass of resonance
cinput bi : branching ration at peak
cinput b1 : itype of 1st outgoing particle
cinput b2 : itype of 2nd outgoing particle
cinput b3 : itype of 3rd outgoing particle
cinput b4 : itype of 4th outgoing particle
c
c     {\tt fbrancx} returns the mass dependent branching ratio for
c     the decay channel {\tt i} of resonance {\tt ir}. This 
c     branching ratio is NOT normalized. To extract the mass dependent
c     decay width, use {\tt fprwdt}.
c    
c {\tt fbrancx} =$
c        \left( \Gamma^{i,j}_{R} \frac{M_{R}}{M}
c        \left( \frac{\langle p_{i,j}(M) \rangle}
c                    {\langle p_{i,j}(M_{R}) \rangle} \right)^{2l+1}
c         \frac{1.2}{1+ 0.2 
c        \left( \frac{\langle p_{i,j}(M) \rangle}
c                    {\langle p_{i,j}(M_{R}) \rangle} \right)^{2l} }
c          \right)  $
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       implicit none
       real*8 kdiv1,kdiv2,em,b,mmin,mn,m1m,m2m
       real*8 bi,minwid
       real*8 fbran,splint,splintth,pmean
       integer i,ires,ir,izr,b1,b2,b3,b4
       include 'comres.f'
       include 'comwid.f'
       include 'options.f'

       real*8 mminit,massit
       integer isoit,flbr,ipwr,ipwr1

       ires=iabs(ir)

       if(iabs(izr).gt.isoit(ires))then 
          fbrancx=0d0
          return
       end if

       if(CTOption(8).ne.0)then
          fbrancx=bi
          return
       end if

       m1m=mminit(b1)
       m2m=mminit(b2)
c in case of three or four particle decays put masses in m2m
       if(b3.ne.0) m2m=m2m+mminit(b3)
       if(b4.ne.0) m2m=m2m+mminit(b4)
       mn=massit(ires)
       mmin= m1m+m2m            ! minimal mass of resonance

       if (wtabflg.ge.2.and.CTOption(33).eq.0) then
          minwid=min(fbran(i,ires),1D-8)
          if (ires.ge.minbar.and.ires.le.maxbar) then !baryons
c branching ratios are continued horicontally outside the spline region
             if(em.le.maxtab2)then
                b=max(splintth(tabx,pbtaby(1,1,ires,i),
     .               pbtaby(1,2,ires,i),widnsp,em,mmin),minwid)
             else
                b=max(splintth(tabx,pbtaby(1,1,ires,i),
     .               pbtaby(1,2,ires,i),widnsp,maxtab2,mmin),minwid)
             endif
          else if (ires.ge.minmes.and.ires.le.maxmes) then !mesons
             if (em.le.maxtab2) then
c branching ratios are continued horicontally outside the spline region
                b=max(splint(tabx,pmtaby(1,1,ires,i),
     .               pmtaby(1,2,ires,i),widnsp,em),minwid)
             else
                b=max(splint(tabx,pmtaby(1,1,ires,i),
     .               pmtaby(1,2,ires,i),widnsp,maxtab2),minwid)
             endif
          else 
             write (6,*) '*** error(fbrancx) wrong id:',ir
             b=0
          endif
       else
         b=0d0
         if (bi.gt.0.and.em.gt.mmin.and.mn.gt.mmin) then

           ipwr=flbr(i,ires)
           ipwr1=ipwr+1

c determine expectation values of outgoing masses  
c call of pmean with -99 instead of iso3 to ensure usage of fixed
c resonance widths: 5% error, but avoids recursion via call
c to fwidth from pmean
           if(CTOption(33).ne.0)then
              kdiv1=pmean(em,b1,-99,b2,-99,b3,-99,b4,-99,ipwr1)/
     &             pmean(mn,b1,-99,b2,-99,b3,-99,b4,-99,ipwr1)
              kdiv2=pmean(em,b1,-99,b2,-99,b3,-99,b4,-99,ipwr)/
     &             pmean(mn,b1,-99,b2,-99,b3,-99,b4,-99,ipwr)
           else
              kdiv1=pmean(em,b1,isoit(b1),b2,isoit(b2),
     &                       b3,isoit(b3),b4,isoit(b4),ipwr1)/
     &              pmean(mn,b1,isoit(b1),b2,isoit(b2),
     &                       b3,isoit(b3),b4,isoit(b4),ipwr1)
              kdiv2=pmean(em,b1,isoit(b1),b2,isoit(b2),
     &                       b3,isoit(b3),b4,isoit(b4),ipwr)/
     &              pmean(mn,b1,isoit(b1),b2,isoit(b2),
     &                       b3,isoit(b3),b4,isoit(b4),ipwr)
           end if
           b=bi*mn/em*kdiv1*1.2/(1.+0.2*kdiv2)
         else
           b=0. 
         end if
       end if
       fbrancx=b
       return
       end

