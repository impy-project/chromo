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



       program SOPHIA20

c************************************************************************
c*** Simulation Of PhotoHadronic Interactions in Astrophysics ***********
c***           VERSION 2.0                                    ***********
c************************************************************************

c************************************************************************
c** Main program for photopion production of relativistic nucleons     **
c**  in a radiation field (blackbody or power law)                     **
c************************************************************************
c** Date: 20/01/98          **
c** correct.:19/02/98       **
c** first release to        ** 
c**  collab.:      25/05/98 **
c** Version 1.3:   31/08/98 **
c** Version 1.4:   12/10/98 **
c** change to DP:  16/11/98 **
c** last corr.:       07/99 **
c** authors: A.G.F. Muecke  **
c**          R.R. Engel     **
c**   in collaboration with **
c**         R.J. Protheroe  **
c**         J.P. Rachen     **
c**         T.S. Stanev     **
c*****************************

c****** INPUT ***********************************************************************
c E0 = energy of incident proton (in lab frame) [in GeV]
c L0 = code number of the incident nucleon (L0=13: proton, L0=14: neutron)
c soft photon field:   for thermal spectrum: temperature tbb [in K]
c                      for PL: (set tbb = 0) alpha = PL index (S \sim nu^{-alpha}), 
c                                            epsmin = low energy cut off  [in eV]
c                                            epsmax = high energy cut off [in eV]
c****** OUTPUT **********************************************************************
c** energy distribution P(x) (logarithmic scale, x=E_particle/E0) for ***
c** photons, protons, neutrons, e-neutrinos, nu-neutrinos             ***
c************************************************************************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       SAVE

       COMMON/input/ tbb,E0,alpha1,alpha2,
     &           epsm1,epsm2,epsb,L0
       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
       COMMON /S_MASS1/ AM(49), AM2(49)
       COMMON /S_CHP/  S_LIFE(49), ICHP(49), ISTR(49), IBAR(49)
       COMMON /S_CSYDEC/ CBR(102), IDB(49), KDEC(612), LBARP(49)

      CHARACTER*6 NAMPRES
      COMMON /RES_PROP/ AMRES(9), SIG0(9),WIDTH(9), 
     +                    NAMPRES(0:9)

      CHARACTER*6 NAMPRESp
      COMMON /RES_PROPp/ AMRESp(9), BGAMMAp(9),WIDTHp(9),  
     +                    RATIOJp(9),NAMPRESp(0:9)

      CHARACTER*6 NAMPRESn
      COMMON /RES_PROPn/ AMRESn(9), BGAMMAn(9),WIDTHn(9),  
     +                    RATIOJn(9),NAMPRESn(0:9)


       DIMENSION DNg(201),DNnum(201),DNnue(201),DNp(201),DNn(201)
       DIMENSION DNnuma(201),DNnuea(201),E0_arr(101)
       DIMENSION Dg(101,201),Dnum(101,201),Dnue(101,201)
       DIMENSION Dp(101,201),Dn(101,201),Dnuma(101,201),Dnuea(101,201)
       DIMENSION Dem(101,201),Dep(101,201),DNem(201),DNep(201)

       dimension c_feps(10000),eps_arr(10000)

       character*6 nameinc
       character*1 ans

       external sample_eps,sample_s,eventgen,
     &  listdistr,output,initial,prob_epskt

       DATA pi /3.141593D0/ 

        do i=1,201
          DNg(i) = 0.D0
          DNnum(i) = 0.D0
          DNnue(i) = 0.D0
          DNp(i) = 0.D0
          DNn(i) = 0.D0
          DNem(i) = 0.D0
          DNep(i) = 0.D0
          DNnuma(i) = 0.D0
          DNnuea(i) = 0.D0

          do jm=1,101
            Dg(jm,i) = 0.D0
            Dnum(jm,i) = 0.D0
            Dnue(jm,i) = 0.D0
            Dp(jm,i) = 0.D0
            Dn(jm,i) = 0.D0
            Dem(jm,i) = 0.D0
            Dep(jm,i) = 0.D0
            Dnuma(jm,i) = 0.D0
            Dnuea(jm,i) = 0.D0
         enddo

       enddo


c****** INPUT **************************************************
       print*
       print*,'Give in code number of incident nucleon: '
       print*,' (13 = proton, 14 = neutron)'
       read(*,*) L0

       delE = 0.D0
       print*
       print*,'Incident nucleon spectrum:'
       print*,'energy grid [s] or single nucleon [n] ? '
       read(*,'(A1)') ans

       if (ans.eq.'n') then
        jm = 1
        ninc = 1
 4     print*,'Give energy [in GeV] of incident nucleon: '
       read(*,*) E0

       if (E0.lt.AM(L0)) then
        print*,'No valid input !'
        goto 4
       endif
       Emin = log10(E0)
       Emax = Emin
       endif
       if (ans.eq.'s') then
  5    print*,'Give in low-energy cut off (in GeV) of nucleon 
     &  energy grid:'
       read(*,*) Emin

       if (Emin.lt.AM(L0)) then
        print*,'No valid input !'
        goto 5
       endif
       Emin = log10(Emin)
       print*,'Give in high-energy cut off (in GeV) of nucleon 
     &  energy grid:'
       read(*,*) Emax

       if (Emax.lt.Emin) then
        print*,'No valid input !'
        goto 5
       endif
       Emax = log10(Emax)
       print*,'Give number of bins (< 101):'
       read(*,*) ninc
       delE = (Emax-Emin)/ninc
       ninc = ninc+1
       endif

       print*
  2    print*,'Give soft photon spectrum: '
       print*,' blackbody spectrum ? (y/n): '
       read(*,'(A1)') ans

       if (ans.eq.'y') then
         rad_den = 1.D0
         print*,'Give temperature [in K]: '
         read(*,*) tbb

       else if (ans.eq.'n') then
       tbb = 0.D0
       print*,' for power law spectrum n(eps) ~ eps^-alpha:'
       print*,'  (alpha = 0 ... 3)'
  8    print*,'straight power law ? (y/n): '

       read(*,'(A1)') ans
       if (ans.eq.'y') then
  3    print*,'   spectral index alpha = '
       read(*,*) alpha2

        if (alpha2.lt.0.D0.or.alpha2.gt.3.D0) then
         print*,'no valid input !'
         goto 3
        endif
       alpha1 = alpha2
       else if (ans.eq.'n') then
        print*,'broken power law:'
        print*,' spectral index alpha1, alpha2 = '
        read(*,*) alpha1,alpha2
        if (alpha1.lt.0.D0.or.alpha1.gt.3.D0.or.
     &  alpha2.lt.0.D0.or.alpha2.gt.3.D0) then
         print*,'no valid input !'
         goto 8
        endif
        print*,'break energy [eV] = '
        read(*,*) epsb        
       else
         print*,'no valid input !'
         goto 8
       endif 
       print*,'low-energy cut off [eV] = '
       read(*,*) epsmin
       print*,'high-energy cut off [eV] = '
       read(*,*) epsmax

       if (ans.eq.'y') epsb = epsmin/100.D0
       else
        print*,'no valid input !'
        goto 2
       endif
      
       print*
       print*,'Give number of trials: '
       read(*,*) ntrial

       print*
       print*,
     &'Give number of bins (<201) for output particle spectra: '
       print*,
     &'(spectra in logarithmic equal bins, with stepsize Delta x'
       print*,'  with x=E_particle/E_initial nucleon )'
       read(*,*) nbins

       print*
       print*,
     &'Give stepsize Delta x for output particle spectra: '
       read(*,*) delx

       print*
       print*,'Give filename (< 7 CH) for output: '
       read(*,'(A6)') nameinc
       print*

       call initial(L0)

c... nucleon energy loop:
       do jm=1,ninc
         Elog = Emin+(jm-1)*delE
         E0 = 10.D0**Elog
         E0_arr(jm) = E0

c... trial loop:
       do nt=1,ntrial 

c*******************************************************
c sample epsilon = energy of photon in lab frame
c*******************************************************
       pm = AM(L0)
  6    call sample_eps(epseV,epsmin,epsmax)
       if (epseV.le.0.) goto 7
c *** eps in GeV:
       eps = epseV/1.D9
       Etot = E0+eps

c*******************************************************
c sample s = total mass of center of momentum frame
c*******************************************************
 18    call sample_s(s,eps)
       gammap = E0/pm
       betap = sqrt(1.D0-1.D0/gammap/gammap)
       theta = ((pm*pm-s)/2.D0/E0/eps+1.D0)/betap
       if (abs(theta).gt.1.D0) STOP
       theta = acos(theta)*180.D0/pi 

c***********************************************************
c*** CALL PHOTOPION EVENT GENERATOR 
c***********************************************************

        call eventgen(L0,E0,eps,theta,Imode)

c*********************************************************
c store decayed particles and their energy distribution: *
c*********************************************************
c... calculate particle distribution dN/dlog(f) with f=E/E0:
c    [note: E*dN/dE = dN/dlog(f)]

        call listdistr(E0,DNg,DNnum,DNnuma,DNnue,DNnuea,
     &   DNp,DNn,DNem,DNep,nbins,delx)

        xini = -nbins*delx
        do 55 i=1,nbins
         Dg(jm,i) = Dg(jm,i)+DNg(i)/delx/ntrial
         Dnum(jm,i) = Dnum(jm,i)+DNnum(i)/delx/ntrial
         Dnue(jm,i) = Dnue(jm,i)+DNnue(i)/delx/ntrial
         Dnuma(jm,i) = Dnuma(jm,i)+DNnuma(i)/delx/ntrial
         Dnuea(jm,i) = Dnuea(jm,i)+DNnuea(i)/delx/ntrial
         Dp(jm,i) = Dp(jm,i)+DNp(i)/delx/ntrial
         Dn(jm,i) = Dn(jm,i)+DNn(i)/delx/ntrial
         Dem(jm,i) = Dem(jm,i)+DNem(i)/delx/ntrial
         Dep(jm,i) = Dep(jm,i)+DNep(i)/delx/ntrial
  55    continue
      
c... trial loop:
       enddo
  7    continue

c... nucleon energy loop:
       enddo

       call output(Dg,Dnum,Dnuma,Dnue,Dnuea,Dp,Dn,Dem,Dep,
     &  nbins,ninc,nameinc,delx,Emin,Emax,E0_arr,epsmin,epsmax)

       STOP
       

       END
