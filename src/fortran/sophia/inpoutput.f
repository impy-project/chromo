c*****************************************************************************
c**!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***
c**!!              IF YOU USE THIS PROGRAM, PLEASE CITE:                 !!***
c**!! A.M"ucke, Ralph Engel, J.P.Rachen, R.J.Protheroe and Todor Stanev, !!***
c**!!  1999, astro-ph/9903478, to appear in Comp.Phys.Commun.            !!***
c**!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***
c*****************************************************************************
c** Further SOPHIA related papers:                                         ***
c** (1) M"ucke A., et al 1999, astro-ph/9808279, to appear in PASA.        ***
c** (2) M"ucke A., et al 1999, to appear in: Proc. of the                  ***
c**      19th Texas Symposium on Relativistic Astrophysics, Paris, France, ***
c**      Dec. 1998. Eds.: J.~Paul, T.~Montmerle \& E.~Aubourg (CEA Saclay) ***
c** (3) M"ucke A., et al 1999, astro-ph/9905153, to appear in: Proc. of    ***
c**      19th Texas Symposium on Relativistic Astrophysics, Paris, France, ***
c**      Dec. 1998. Eds.: J.~Paul, T.~Montmerle \& E.~Aubourg (CEA Saclay) ***
c** (4) M"ucke A., et al 1999, to appear in: Proc. of 26th Int.Cosmic Ray  ***
c**      Conf. (Salt Lake City, Utah)                                      ***
c*****************************************************************************


c*********************************
c*** Routines related to output: *
c*********************************


       subroutine LISTDISTR(E0,Dg,Dnum,Dnuma,Dnue,Dnuea,Dp,Dn,Dem,
     &                    Dep,nbins,delx)

c*********************************************************************
c** calculates distribution of energy of given particle to incident **
c** proton energy; considered particles are:                        **
c** photons, protons, neutrons, e-neutrinos, nu-neutrinos           **
c** Note: Dg(),Dnum(),Dnue(),Dp(),Dn(),Dem(),Dep(),Dnuea(),Dnuma()  **
c**       gives # of photons per logarithmic bin width: dN/dlog(f)  **
c*********************************************************************
c** Date: 20/01/98   **
c** author: A.Muecke **
c**********************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       SAVE

       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
       DIMENSION Dg(201),Dnum(201),Dnue(201),Dp(201),Dn(201)
       DIMENSION Dem(201),Dep(201),Dnuea(201),Dnuma(201)

        do i=1,201
          Dg(i) = 0.
          Dnum(i) = 0.
          Dnue(i) = 0.
          Dp(i) = 0.
          Dn(i) = 0.
          Dem(i) = 0.
          Dep(i) = 0.
          Dnuma(i) = 0.
          Dnuea(i) = 0.
       enddo

       xini = -nbins*delx
c go through LLIST:
       do 10 i=1,NP
        LA = abs(LLIST(i))
        EI = abs(P(i,4))
        Ep = E0/1.D10
        r = EI/Ep/1.D10
        x = log10(r)
        if (LA.lt.10000) then
c** gamma ray distribution
        if (LA.eq.1) then
        do 20 j=1,nbins
         x1 = xini+delx*(j-1)
         x2 = xini+delx*j
         if (x.ge.x1.and.x.lt.x2) then
          Dg(j) = Dg(j)+1.D0
         endif
         if (x.eq.0.D0) then
          Dg(nbins) = Dg(nbins)+1.D0
         endif
 20     continue
        endif
c** neutron distribution
        if (LA.eq.14) then
        do 21 j=1,nbins
         x1 = xini+delx*(j-1)
         x2 = xini+delx*j
         if (x.ge.x1.and.x.lt.x2) then
          Dn(j) = Dn(j)+1.D0
         endif
         if (x.eq.0.D0) then
          Dn(nbins) = Dn(nbins)+1.D0
         endif
 21     continue
        endif
c** proton distribution
        if (LA.eq.13) then
        do 22 j=1,nbins
         x1 = xini+delx*(j-1)
         x2 = xini+delx*j
         if (x.ge.x1.and.x.lt.x2) then
          Dp(j) = Dp(j)+1.D0
         endif
         if (x.eq.0.D0) then
          Dp(nbins) = Dp(nbins)+1.D0
         endif
 22     continue
        endif
c** neutrino distribution
        if (LA.eq.17) then
        do 23 j=1,nbins
         x1 = xini+delx*(j-1)
         x2 = xini+delx*j
         if (x.ge.x1.and.x.lt.x2) then
          Dnum(j) = Dnum(j)+1.D0
         endif
         if (x.eq.0.D0) then
          Dnum(nbins) = Dnum(nbins)+1.D0
         endif
 23     continue
        endif

        if (LA.eq.18) then
        do 27 j=1,nbins
         x1 = xini+delx*(j-1)
         x2 = xini+delx*j
         if (x.ge.x1.and.x.lt.x2) then
          Dnuma(j) = Dnuma(j)+1.D0
         endif
         if (x.eq.0.D0) then
          Dnuma(nbins) = Dnuma(nbins)+1.D0
         endif
 27     continue
        endif


        if (LA.eq.15) then
        do 24 j=1,nbins
         x1 = xini+delx*(j-1)
         x2 = xini+delx*j
         if (x.ge.x1.and.x.lt.x2) then
          Dnue(j) = Dnue(j)+1.D0
         endif
         if (x.eq.0.D0) then
          Dnue(nbins) = Dnue(nbins)+1.D0
         endif
 24     continue
        endif

        if (LA.eq.16) then
        do 28 j=1,nbins
         x1 = xini+delx*(j-1)
         x2 = xini+delx*j
         if (x.ge.x1.and.x.lt.x2) then
          Dnuea(j) = Dnuea(j)+1.D0
         endif
         if (x.eq.0.D0) then
          Dnuea(nbins) = Dnuea(nbins)+1.D0
         endif
 28     continue
        endif

c** electron distribution
        if (LA.eq.3) then
        do 25 j=1,nbins
         x1 = xini+delx*(j-1)
         x2 = xini+delx*j
         if (x.ge.x1.and.x.lt.x2) then
          Dem(j) = Dem(j)+1.D0
         endif
         if (x.eq.0.D0) then
          Dem(nbins) = Dem(nbins)+1.D0
         endif
 25     continue
        endif

c** positron distribution
        if (LA.eq.2) then
        do 26 j=1,nbins
         x1 = xini+delx*(j-1)
         x2 = xini+delx*j
         if (x.ge.x1.and.x.lt.x2) then
          Dep(j) = Dep(j)+1.D0
         endif
         if (x.eq.0.D0) then
          Dep(nbins) = Dep(nbins)+1.D0
         endif
 26     continue
        endif

        endif
 10     continue

        RETURN

        END

       subroutine output(Dg,Dnum,Dnuma,Dnue,Dnuea,Dp,Dn,Dem,Dep,nbins,
     &   ninc,nameinc,delx,Emin,Emax,E0_arr,epsmin,epsmax)

c********************************************************************
c*** OUTPUT ROUTINE for particle spectra:                     *******
c*** considered particles:                                    *******
c*** photons, protons, neutrons, e-neutrinos, nu-neutrinos,   *******
c*** electron, positron                                       *******
c*** spectra of each particle stored in separate files        *******
c********************************************************************
c** Date: 20/02/98   **
c** author: A.Muecke **
c**********************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-M)
       SAVE

       COMMON/input/ tbb,E0,alpha1,alpha2,
     &           epsm1,epsm2,epsb,L0

       DIMENSION Dg(101,201),Dnum(101,201),Dnue(101,201)
       DIMENSION Dp(101,201),Dn(101,201),Dnuma(101,201),E0_arr(101)
       DIMENSION Dem(101,201),Dep(101,201),Dnuea(101,201)
       character*5 particle
       character*7 spart,fpart
       character*13 filename
       character*6 nameinc
       character mat*20, strnm1*2
       character mat1*20, strnm11*2
       character mat2*20, strnm12*2

 571  format(2(3x,E10.5),3x,I3,5(3x,E10.5),
     &     3x,I3,3x,E10.5,3x,A10,3x,A10)
 572  format(E10.5,3x,2(I3,3x))
 573  format(2x,E10.5)

      print*
      print*,'OUTPUT files:'
c**********************************************
c******** GAMMA spectra: **********************
      particle = 'gamma'
      filename =  nameinc // '.' // particle
      print*,'filename = ',filename
      if (L0.eq.13) spart = 'proton'
      if (L0.eq.14) spart = 'neutron'
      fpart = 'photon'
      if (tbb.gt.0.) then 
        target1 = tbb
        target2 = 0.
      else
        target1 = alpha1
        target2 = alpha2
      endif
      open(1,file=filename)
c... write input parameters:
       write(1,571) Emin,Emax,ninc,target1,target2,
     &     epsmin,epsb,epsmax,nbins,delx,spart,fpart
c... nucleon energy loop:
       do i=1,ninc
c ... determine j-range = range of energy bins not equal zero
        jini = 0
        jfin = 0
c... particle spectrum loop:
        do j=1,nbins
        if (Dg(i,j).gt.0.D0) then
         jfin = j
         if (jini.eq.0) jini = j
        endif
        enddo
      nm = jfin-jini+1
      nm1 = nm+1
      write(strnm1,'(I2)') nm1
      mat = '(' // strnm1 // '(5X,E10.4))'
      nmc = 81
      nmc2 = 82
      write(strnm11,'(I2)') nmc2
      mat1 = '(' // strnm11 // '(5X,E10.4))'
      nmcf = jfin-jini-80+1
      nmcf2 = nmcf+1
      write(strnm12,'(I2)') nmcf2
      mat2 = '(' // strnm12 // '(5X,E10.4))'
        if (jfin.gt.0) then
        write(1,572) E0_arr(i),jini,jfin
c... values written in one line:
         if (jfin-jini.lt.80) then
          write(1,FMT=mat) (Dg(i,jl),jl=jini,jfin)
         else
          jfin0 = jini+80
          write(1,FMT=mat1) (Dg(i,jl),jl=jini,jfin0)
          write(1,FMT=mat2) (Dg(i,jl),jl=jfin0+1,jfin)
         endif
        endif
       enddo
      close(1)

c**********************************************
c******** MU-NEUTRINO spectra: **********************
      particle = 'muneu'
      filename = nameinc // '.' // particle
      print*,'filename = ',filename
      open(1,file=filename)
       write(1,571) Emin,Emax,ninc,target1,target2,
     &     epsmin,epsb,epsmax,nbins,delx,spart,fpart
c... nucleon energy loop:
       do i=1,ninc
c ... determine j-range = range of energy bins not equal zero
        jini = 0
        jfin = 0
c... particle spectrum loop:
        do j=1,nbins
        if (Dnum(i,j).gt.0.D0) then
         jfin = j
         if (jini.eq.0) jini = j
        endif
        enddo
      nm = jfin-jini+1
      nm1 = nm+1
      write(strnm1,'(I2)') nm1
      mat = '(' // strnm1 // '(5X,E10.4))'
      nmc = 81
      nmc2 = 82
      write(strnm11,'(I2)') nmc2
      mat1 = '(' // strnm11 // '(5X,E10.4))'
      nmcf = jfin-jini-80+1
      nmcf2 = nmcf+1
      write(strnm12,'(I2)') nmcf2
      mat2 = '(' // strnm12 // '(5X,E10.4))'
        if (jfin.gt.0) then
        write(1,572) E0_arr(i),jini,jfin
c... values written in one line:
         if (jfin-jini.lt.80) then
          write(1,FMT=mat) (Dnum(i,jl),jl=jini,jfin)
         else
          jfin0 = jini+80
          write(1,FMT=mat1) (Dnum(i,jl),jl=jini,jfin0)
          write(1,FMT=mat2) (Dnum(i,jl),jl=jfin0+1,jfin)
         endif
        endif
       enddo
      close(1)

      particle = 'muane'
      filename = nameinc // '.' // particle
      print*,'filename = ',filename
      open(1,file=filename)
       write(1,571) Emin,Emax,ninc,target1,target2,
     &     epsmin,epsb,epsmax,nbins,delx,spart,fpart
c... nucleon energy loop:
       do i=1,ninc
c ... determine j-range = range of energy bins not equal zero
        jini = 0
        jfin = 0
c... particle spectrum loop:
        do j=1,nbins
        if (Dnuma(i,j).gt.0.D0) then
         jfin = j
         if (jini.eq.0) jini = j
        endif
        enddo
      nm = jfin-jini+1
      nm1 = nm+1
      write(strnm1,'(I2)') nm1
      mat = '(' // strnm1 // '(5X,E10.4))'
      nmc = 81
      nmc2 = 82
      write(strnm11,'(I2)') nmc2
      mat1 = '(' // strnm11 // '(5X,E10.4))'
      nmcf = jfin-jini-80+1
      nmcf2 = nmcf+1
      write(strnm12,'(I2)') nmcf2
      mat2 = '(' // strnm12 // '(5X,E10.4))'
        if (jfin.gt.0) then
        write(1,572) E0_arr(i),jini,jfin
c... values written in one line:
         if (jfin-jini.lt.80) then
          write(1,FMT=mat) (Dnuma(i,jl),jl=jini,jfin)
         else
          jfin0 = jini+80
          write(1,FMT=mat1) (Dnuma(i,jl),jl=jini,jfin0)
          write(1,FMT=mat2) (Dnuma(i,jl),jl=jfin0+1,jfin)
         endif
        endif
       enddo
      close(1)

c**********************************************
c******** ELECTRON NEUTRINO spectra: **********
      particle = 'e_neu'
      filename =  nameinc // '.' // particle
      print*,'filename = ',filename
      open(1,file=filename)
       write(1,571) Emin,Emax,ninc,target1,target2,
     &     epsmin,epsb,epsmax,nbins,delx,spart,fpart
c... nucleon energy loop:
       do i=1,ninc
c ... determine j-range = range of energy bins not equal zero
        jini = 0
        jfin = 0
c... particle spectrum loop:
        do j=1,nbins
        if (Dnue(i,j).gt.0.D0) then
         jfin = j
         if (jini.eq.0) jini = j
        endif
        enddo
      nm = jfin-jini+1
      nm1 = nm+1
      write(strnm1,'(I2)') nm1
      mat = '(' // strnm1 // '(5X,E10.4))'
      nmc = 81
      nmc2 = 82
      write(strnm11,'(I2)') nmc2
      mat1 = '(' // strnm11 // '(5X,E10.4))'
      nmcf = jfin-jini-80+1
      nmcf2 = nmcf+1
      write(strnm12,'(I2)') nmcf2
      mat2 = '(' // strnm12 // '(5X,E10.4))'
        if (jfin.gt.0) then
        write(1,572) E0_arr(i),jini,jfin
c... values written in one line:
         if (jfin-jini.lt.80) then
          write(1,FMT=mat) (Dnue(i,jl),jl=jini,jfin)
         else
          jfin0 = jini+80
          write(1,FMT=mat1) (Dnue(i,jl),jl=jini,jfin0)
          write(1,FMT=mat2) (Dnue(i,jl),jl=jfin0+1,jfin)
         endif
        endif
       enddo
      close(1)

      particle = 'eaneu'
      filename =  nameinc // '.' // particle
      print*,'filename = ',filename
      open(1,file=filename)
       write(1,571) Emin,Emax,ninc,target1,target2,
     &     epsmin,epsb,epsmax,nbins,delx,spart,fpart
c... nucleon energy loop:
       do i=1,ninc
c ... determine j-range = range of energy bins not equal zero
        jini = 0
        jfin = 0
c... particle spectrum loop:
        do j=1,nbins
        if (Dnuea(i,j).gt.0.D0) then
         jfin = j
         if (jini.eq.0) jini = j
        endif
        enddo
      nm = jfin-jini+1
      nm1 = nm+1
      write(strnm1,'(I2)') nm1
      mat = '(' // strnm1 // '(5X,E10.4))'
      nmc = 81
      nmc2 = 82
      write(strnm11,'(I2)') nmc2
      mat1 = '(' // strnm11 // '(5X,E10.4))'
      nmcf = jfin-jini-80+1
      nmcf2 = nmcf+1
      write(strnm12,'(I2)') nmcf2
      mat2 = '(' // strnm12 // '(5X,E10.4))'
        if (jfin.gt.0) then
        write(1,572) E0_arr(i),jini,jfin
c... values written in one line:
         if (jfin-jini.lt.80) then
          write(1,FMT=mat) (Dnuea(i,jl),jl=jini,jfin)
         else
          jfin0 = jini+80
          write(1,FMT=mat1) (Dnuea(i,jl),jl=jini,jfin0)
          write(1,FMT=mat2) (Dnuea(i,jl),jl=jfin0+1,jfin)
         endif
        endif
       enddo
      close(1)

c**********************************************
c******** ELECTRON spectra: **********************
      particle = 'elect'
      filename =  nameinc // '.' // particle
      print*,'filename = ',filename
      open(1,file=filename)
       write(1,571) Emin,Emax,ninc,target1,target2,
     &     epsmin,epsb,epsmax,nbins,delx,spart,fpart
c... nucleon energy loop:
       do i=1,ninc
c ... determine j-range = range of energy bins not equal zero
        jini = 0
        jfin = 0
c... particle spectrum loop:
        do j=1,nbins
        if (Dem(i,j).gt.0.D0) then
         jfin = j
         if (jini.eq.0) jini = j
        endif
        enddo
      nm = jfin-jini+1
      nm1 = nm+1
      write(strnm1,'(I2)') nm1
      mat = '(' // strnm1 // '(5X,E10.4))'
      nmc = 81
      nmc2 = 82
      write(strnm11,'(I2)') nmc2
      mat1 = '(' // strnm11 // '(5X,E10.4))'
      nmcf = jfin-jini-80+1
      nmcf2 = nmcf+1
      write(strnm12,'(I2)') nmcf2
      mat2 = '(' // strnm12 // '(5X,E10.4))'
        if (jfin.gt.0) then
        write(1,572) E0_arr(i),jini,jfin
c... values written in one line:
         if (jfin-jini.lt.80) then
          write(1,FMT=mat) (Dem(i,jl),jl=jini,jfin)
         else
          jfin0 = jini+80
          write(1,FMT=mat1) (Dem(i,jl),jl=jini,jfin0)
          write(1,FMT=mat2) (Dem(i,jl),jl=jfin0+1,jfin)
         endif
        endif
       enddo
      close(1)

c**********************************************
c******** POSITRON spectra: **********************
      particle = 'posit'
      filename =  nameinc // '.' // particle
      print*,'filename = ',filename
      open(1,file=filename)
       write(1,571) Emin,Emax,ninc,target1,target2,
     &     epsmin,epsb,epsmax,nbins,delx,spart,fpart
c... nucleon energy loop:
       do i=1,ninc
c ... determine j-range = range of energy bins not equal zero
        jini = 0
        jfin = 0
c... particle spectrum loop:
        do j=1,nbins
        if (Dep(i,j).gt.0.D0) then
         jfin = j
         if (jini.eq.0) jini = j
        endif
        enddo
      nm = jfin-jini+1
      nm1 = nm+1
      write(strnm1,'(I2)') nm1
      mat = '(' // strnm1 // '(5X,E10.4))'
      nmc = 81
      nmc2 = 82
      write(strnm11,'(I2)') nmc2
      mat1 = '(' // strnm11 // '(5X,E10.4))'
      nmcf = jfin-jini-80+1
      nmcf2 = nmcf+1
      write(strnm12,'(I2)') nmcf2
      mat2 = '(' // strnm12 // '(5X,E10.4))'
        if (jfin.gt.0) then
        write(1,572) E0_arr(i),jini,jfin
c... values written in one line:
         if (jfin-jini.lt.80) then
          write(1,FMT=mat) (Dep(i,jl),jl=jini,jfin)
         else
          jfin0 = jini+80
          write(1,FMT=mat1) (Dep(i,jl),jl=jini,jfin0)
          write(1,FMT=mat2) (Dep(i,jl),jl=jfin0+1,jfin)
         endif
        endif
       enddo
      close(1)

c**********************************************
c******** PROTON spectra: **********************
      particle = 'proto'
      filename =  nameinc // '.' // particle
      print*,'filename = ',filename
      open(1,file=filename)
       write(1,571) Emin,Emax,ninc,target1,target2,
     &     epsmin,epsb,epsmax,nbins,delx,spart,fpart
c... nucleon energy loop:
       do i=1,ninc
c ... determine j-range = range of energy bins not equal zero
        jini = 0
        jfin = 0
c... particle spectrum loop:
        do j=1,nbins
        if (Dp(i,j).gt.0.D0) then
         jfin = j
         if (jini.eq.0) jini = j
        endif
        enddo
      nm = jfin-jini+1
      nm1 = nm+1
      write(strnm1,'(I2)') nm1
      mat = '(' // strnm1 // '(5X,E10.4))'
      nmc = 81
      nmc2 = 82
      write(strnm11,'(I2)') nmc2
      mat1 = '(' // strnm11 // '(5X,E10.4))'
      nmcf = jfin-jini-80+1
      nmcf2 = nmcf+1
      write(strnm12,'(I2)') nmcf2
      mat2 = '(' // strnm12 // '(5X,E10.4))'
        if (jfin.gt.0) then
        write(1,572) E0_arr(i),jini,jfin
c... values written in one line:
         if (jfin-jini.lt.80) then
          write(1,FMT=mat) (Dp(i,jl),jl=jini,jfin)
         else
          jfin0 = jini+80
          write(1,FMT=mat1) (Dp(i,jl),jl=jini,jfin0)
          write(1,FMT=mat2) (Dp(i,jl),jl=jfin0+1,jfin)
         endif
        endif
       enddo
      close(1)

c**********************************************
c******** NEUTRON spectra: **********************
      particle = 'neutr'
      filename =  nameinc // '.' // particle
      print*,'filename = ',filename
      open(1,file=filename)
       write(1,571) Emin,Emax,ninc,target1,target2,
     &     epsmin,epsb,epsmax,nbins,delx,spart,fpart
c... nucleon energy loop:
       do i=1,ninc
c ... determine j-range = range of energy bins not equal zero
        jini = 0
        jfin = 0
c... particle spectrum loop:
        do j=1,nbins
        if (Dn(i,j).gt.0.D0) then
         jfin = j
         if (jini.eq.0) jini = j
        endif
        enddo
      nm = jfin-jini+1
      nm1 = nm+1
      write(strnm1,'(I2)') nm1
      mat = '(' // strnm1 // '(5X,E10.4))'
      nmc = 81
      nmc2 = 82
      write(strnm11,'(I2)') nmc2
      mat1 = '(' // strnm11 // '(5X,E10.4))'
      nmcf = jfin-jini-80+1
      nmcf2 = nmcf+1
      write(strnm12,'(I2)') nmcf2
      mat2 = '(' // strnm12 // '(5X,E10.4))'
        if (jfin.gt.0) then
        write(1,572) E0_arr(i),jini,jfin
c... values written in one line:
         if (jfin-jini.lt.80) then
          write(1,FMT=mat) (Dn(i,jl),jl=jini,jfin)
         else
          jfin0 = jini+80
          write(1,FMT=mat1) (Dn(i,jl),jl=jini,jfin0)
          write(1,FMT=mat2) (Dn(i,jl),jl=jfin0+1,jfin)
         endif
        endif
       enddo
      close(1)

      RETURN
      END
