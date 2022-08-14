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


c**********************************************
c** Routines/functions related to sampling  ***
c**  photon energy and squared CMF energy:  ***
c**********************************************



       subroutine sample_s(s,eps)

c***********************************************************************
c samples distribution of s: p(s) = (s-mp^2)sigma_Ngamma
c rejection for s=[sth,s0], analyt.inversion for s=[s0,smax]
c***********************************************************************
c** Date: 20/01/98   **
c** author: A.Muecke **
c**********************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       SAVE

       common/input/ tbb,E0,alpha1,alpha2,
     &           epsm1,epsm2,epsb,L0
       COMMON /S_MASS1/ AM(49), AM2(49)

      external functs,gauss,rndm
      double precision functs,gauss,rndm

c*** calculate smin,smax : ************************
      xmpi = AM(7)
      xmp = AM(L0)
      Pp = sqrt(E0*E0-xmp*xmp)
      smin = 1.1646D0
      smax = max(smin,xmp*xmp+2.D0*eps*(E0+Pp))
      if ((smax-smin).le.1.D-8) then
       s = smin+rndm(0)*1.d-6
       RETURN
      endif
      s0 = 10.D0
c*** determine which method applies: rejection or analyt.inversion: **
      sintegr1 = gauss(functs,smin,s0)
      sintegr2 = gauss(functs,s0,smax)
      if (smax.le.s0) then
c rejection method:
       nmethod=1
       goto 41
      endif
      r1 = rndm(0)
      quo = sintegr1/(sintegr1+sintegr2)
      if (r1.le.quo) then
c rejection method:
       nmethod=1
      else
c analyt. inversion:
       nmethod=2
      endif

  41  continue


c*** rejection method: ************************
      if (nmethod.eq.1) then
  10  continue
c*** sample s random between smin ... s0 **
      r2 = rndm(0)
      s = smin+r2*(smax-smin)
c*** calculate p(s) = pes **********************
      ps = functs(s)
c*** rejection method to sample s *********************
        r3 = rndm(0)
c pmax is roughly p(s) at s=s0
        pmax = 1300.D0/sintegr1
        if (r3*pmax.le.ps/sintegr1) then
          RETURN
        else
         goto 10
        endif
       endif

c*** analyt. inversion method: *******************
      if (nmethod.eq.2) then
       r4 = rndm(0)
       beta = 2.04D0
       betai = 1.D0/beta
       term1 = r4*(smax**beta)
       term2 = (r4-1.D0)*(s0**beta)
       s = (term1-term2)**betai
       RETURN
      endif

       RETURN
       END

       subroutine sample_eps(eps,epsmin,epsmax)

c****************************************************************************
c samples distribution of epsilon p(epsilon) for blackbody spectrum if tbb>0
c  and power law \sim eps^-alpha, epsm1 [eV] < eps [eV] < epsm2 [eV], 
c                       eps in LAB frame if tbb \leq 0
c****************************************************************************
c** Date: 20/01/98   **
c** author: A.Muecke **
c**********************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       SAVE

       common/input/ tbb,E0,alpha1,alpha2,
     &           epsm1,epsm2,epsb,L0
       common/PLindex/ alphaxx
       COMMON /S_MASS1/ AM(49), AM2(49)

      external prob_epskt,prob_epspl,rndm,gauss,functs,probint_pl
      double precision prob_epskt,prob_epspl,rndm,gauss
      double precision functs,probint_pl

      xmpi = AM(7)
      xmp = AM(L0)
      gammap = E0/xmp
      betap = sqrt(1.-1./gammap/gammap)
      Pp = sqrt(E0*E0-xmp*xmp)
  1   continue
      facpmax = 1.6D0
c*** for tbb<=0: power law photon spectrum n(eps) ~ eps^-alpha *********
      if (tbb.gt.0.D0) then
       epsm1 = (1.1646D0-xmp*xmp)/2.D0/(E0+Pp)
       epsm1 = epsm1*1.D9
       epsm2 = 0.007D0*tbb
       epskt = 8.619D-5*tbb
       epspmax = (3.D-3*((E0*epskt*1.D-9)**(-0.97D0))
     &              +0.047D0)/3.9D2*tbb
        if (epsm1.gt.epsm2) then
         print*,
     & 'CMF energy is below threshold for nucleon energy '
     &   ,E0,' GeV !'
         eps = 0.D0
         RETURN
        endif
        cnorm = gauss(prob_epskt,epsm1,epsm2)
         pmaxc = prob_epskt(epspmax)/cnorm
         pmax = facpmax*pmaxc
 
        else
c*** determine distribution:
       epsth = (1.1646D0-xmp*xmp)/2.D0/(E0+Pp)
       epsth = epsth*1.D9
       epsm1 = max(epsmin,epsth)
       epsm2 = epsmax
       if (epsm1.ge.epsm2) then
        eps = 0.
        RETURN
       endif
      endif

      epsmx1 = epsm1
      epsmx2 = epsm2
      epsbx = epsb
      epsdelta = 0.159368/E0*1.d9
      epsxx = 126.D0/E0*1.d9
      alpha12 = alpha2-alpha1
      a1 = 1.D0
  10  continue
c*** sample eps randomly between epsmin ... epsmax **
      r1 = rndm(0)
c*** calculate p(eps) = peps ***************************************
      if (tbb.le.0.D0) then
       rn = rndm(0)
c*******************************************************************
c... sample from straight power law (alpha = alpha2, epsb < epsm1):
      if (alpha12.eq.0.D0.or.epsm1.ge.epsb) then
       if (epsxx.ge.epsm2) then
        alphaxx = alpha2
       else if (epsxx.le.epsm1) then
        alphaxx = (alpha2+2.D0)
       else if (epsm1.lt.epsxx.and.epsxx.lt.epsm2) then
          a2 = epxx*epsxx
          alphaxx = alpha2
          pintegr1 = a1*probint_pl(epsm1,epsxx,alphaxx)
          alphaxx = (alpha2+2.D0)
          pintegr2 = a2*probint_pl(epsxx,epsm2,alphaxx)
          pintegr1 = pintegr1/(pintegr1+pintegr2)
          if (rn.lt.pintegr1) then
            alphaxx = alpha2 
            epsm2 = epsxx 
            ampl = a1          
          else if (pintegr1.le.rn.and.rn.lt.1.D0) then
            alphaxx = alpha2+2.D0
            epsm1  = epsxx
            ampl = a2 
          endif
       endif    
      endif
c... sample from broken power law: input always epsm1 < epsb < epsm2 
      if (epsm1.lt.epsb) then
c... determine where epsb,epsxx lies:
        if (epsm1.lt.epsxx.and.epsxx.lt.epsb) then
          a2 = epxx*epsxx
          a3 = a2*(epsb**(alpha2-alpha1))
          alphaxx = alpha1
          pintegr1 = a1*probint_pl(epsm1,epsxx,alphaxx)
          alphaxx = (alpha1+2.D0)
          pintegr2 = a2*probint_pl(epsxx,epsb,alphaxx)
          alphaxx = (alpha2+2.D0)
          pintegr3 = a3*probint_pl(epsb,epsm2,alphaxx)
          pintegr1 = pintegr1/(pintegr1+pintegr2+pintegr3)
          pintegr2 = (pintegr1+pintegr2)/(pintegr1+pintegr2+pintegr3)
          pintegr3 = 1.D0
          if (rn.lt.pintegr1) then
            alphaxx = alpha1 
            epsm2 = epsxx 
            ampl = a1          
          else if (pintegr1.le.rn.and.rn.lt.pintegr2) then
            alphaxx = alpha1+2.D0
            epsm1  = epsxx
            epsm2 = epsb
            ampl = a2 
          else if (pintegr2.le.rn.and.rn.le.pintegr3) then
            alphaxx = alpha2+2.D0
            epsm1 = epsb
            ampl = a3 
          else
            print*,'error in sampling broken power law: SAMPLE_EPS (1)!'
            STOP
          endif

         else if (epsb.le.epsxx.and.epsxx.lt.epsm2) then
          a2 = epsb**(alpha2-alpha1)
          a3 = a2*epsxx*epsxx
          alphaxx = alpha1
          pintegr1 = a1*probint_pl(epsm1,epsb,alphaxx)
          alphaxx = alpha2
          pintegr2 = a2*probint_pl(epsb,epsxx,alphaxx)
          alphaxx = (alpha2+2.D0)
          pintegr3 = a3*probint_pl(epsxx,epsm2,alphaxx)
          pintegr1 = pintegr1/(pintegr1+pintegr2+pintegr3)
          pintegr2 = (pintegr1+pintegr2)/(pintegr1+pintegr2+pintegr3)
          pintegr3 = 1.D0
          if (rn.lt.pintegr1) then
            alphaxx = alpha1 
            epsm2 = epsb 
            ampl = a1         
          else if (pintegr1.le.rn.and.rn.lt.pintegr2) then
            alphaxx = alpha2
            epsm1  = epsb
            epsm2 = epsxx
            ampl = a2 
          else if (pintegr2.le.rn.and.rn.le.pintegr3) then
            alphaxx = alpha2+2.D0
            epsm1 = epsxx
            ampl = a3 
          else
            print*,'error in sampling broken power law: SAMPLE_EPS (2)!'
            STOP
          endif

         else if (epsxx.ge.epsm2) then
          a2 = epsb**(alpha2-alpha1)
          a3 = 0.D0
          alphaxx = alpha1
          pintegr1 = a1*probint_pl(epsm1,epsb,alphaxx)
          alphaxx = alpha2
          pintegr2 = a2*probint_pl(epsb,epsm2,alphaxx)
          pintegr1 = pintegr1/(pintegr1+pintegr2)
          pintegr2 = 1.D0
          if (rn.lt.pintegr1) then
            alphaxx = alpha1 
            epsm2 = epsb 
            ampl = a1         
          else if (pintegr1.le.rn.and.rn.le.pintegr2) then
            alphaxx = alpha2
            epsm1 = epsb
            ampl = a2 
          else
            print*,'error in sampling broken power law: SAMPLE_EPS (3)!'
            STOP
          endif

         else if (epsxx.le.epsm1) then
          a2 = epsb**(alpha2-alpha1)
          a3 = 0.D0
          alphaxx = (alpha1+2.D0)
          pintegr1 = a1*probint_pl(epsm1,epsb,alphaxx)
          alphaxx = (alpha2+2.D0)
          pintegr2 = a2*probint_pl(epsb,epsm2,alphaxx)
          pintegr1 = pintegr1/(pintegr1+pintegr2)
          pintegr2 = 1.D0
          if (rn.lt.pintegr1) then
            alphaxx = alpha1+2.D0 
            epsm2 = epsb 
            ampl = a1         
          else if (pintegr1.le.rn.and.rn.le.pintegr2) then
            alphaxx = alpha2+2.D0
            epsm1 = epsb
            ampl = a2 
          else
            print*,'error in sampling broken power law: SAMPLE_EPS (4)!'
            STOP
          endif

         endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c... END: sample from broken power law:
       endif
c*****************************************************     
       if (alphaxx.eq.1.D0) then
        term1 = r1*log(epsm2/epsm1)
        eps = epsm1*exp(term1)
       else
        beta = 1.D0-alphaxx
        betai = 1.D0/beta
        term1 = r1*(epsm2**beta)
        term2 = (r1-1.D0)*(epsm1**beta)
        eps = (term1-term2)**betai
       endif


c******************************************************
c*** for thermal spectrum: ***
      else
       eps = epsm1+r1*(epsm2-epsm1)
       peps = prob_epskt(eps)/cnorm
c      endif
c*** rejection method to sample eps *********************
        r2 = rndm(0)
        if (r2*pmax.gt.peps) then
         goto 10
        endif

       endif

       epsm1 = epsmx1
       epsm2 = epsmx2
       epsb = epsbx

c... check maximum of epsilon distribution:
       if (pmax.lt.peps) then
        facpmax = facpmax + 0.1D0
        goto 1
       endif

       RETURN
       END


        DOUBLE PRECISION function prob_epskt(eps)

c*** calculates probability distribution for thermal photon field ***
c*** with temerature tbb (in K), eps (in eV)            *************
c** Date: 20/01/98   **
c** author: A.Muecke **
c**********************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)

       SAVE

       common/input/ tbb,E0,alpha1,alpha2,
     &           epsm1,epsm2,epsb,L0

       external functs,photd,gauss
       double precision functs,photd,gauss

       xmpi = 0.135D0
       xmp = 0.93827D0
       Pp = sqrt(E0*E0-xmp*xmp)
       gammap = E0/xmp
       betap = sqrt(1.D0-1.D0/gammap/gammap)
       deps = photd(eps,tbb)
       if (deps.eq.0.D0) then
        prob_epskt = 0.D0
        RETURN
       else
c*** calculate \int_sth^smax ds (s-mp^2) sigma_pg *******
c*** smin is for head-on collision **********************
        smin = 1.1646D0
        smax = max(smin,xmp*xmp+2.D0*eps/1.D9*E0*(1.D0+betap))
        sintegr = gauss(functs,smin,smax)

        prob_epskt = deps/eps/eps*sintegr/
     &     8.D0/betap/E0/E0*1.D18*1.D6
       endif

        RETURN

        END

        DOUBLE PRECISION function prob_epspl(eps)

c*** calculates probability distribution for power law photon field ***
c*** n = anorm*eps^-alpha, eps=[epsm1,epsm2], eps in eV *************
c** Date: 20/01/98   **
c** author: A.Muecke **
c**********************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)

        SAVE

       common/input/ tbb,E0,alpha1,alpha2,
     &           epsm1,epsm2,epsb,L0

       external functs,gauss
       double precision functs,gauss

       xmpi = 0.135D0
       xmp = 0.93827D0
       Pp = sqrt(E0*E0-xmp*xmp)
       gammap = E0/xmp
       betap = sqrt(1.D0-1.D0/gammap/gammap)
       alpha12 = alpha2-alpha1
       ampl = epsb**alpha12
       if (eps.lt.epsb) then
         deps = eps**(-alpha1)
       else
         deps = ampl*(eps**(-alpha2))
       endif

c*** calculate \int_sth^smax ds (s-mp^2) sigma_pg *******
c*** smin is for head-on collision **********************
        smin = 1.1646D0
        smax = max(smin,xmp*xmp+2.D0*eps/1.D9*(E0+Pp))

         sintegr = gauss(functs,smin,smax)

        prob_epspl = deps/eps/eps*sintegr/
     &     8.D0/betap/E0/E0*1.D18*1.D6

        RETURN

        END

        DOUBLE PRECISION function probint_pl(epsi,epsf,p)

c*** returns \int_epsi^epsf eps^-p   **************
c*** calling program is SAMPLE_EPS ****************
c** Date: 03/03/99   **
c** author: A.Muecke **
c**********************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)

        SAVE

        if (p.eq.1.D0) then
          probint_pl = log(epsf/epsi)
        else
         p1 = 1.D0-p
         probint_pl = ((epsf**p1)-(epsi**p1))/p1
        endif
  
        RETURN

        END


        DOUBLE PRECISION function functs(s)

c*** returns (s-pm^2)*sigma_Ngamma **************
c*** calling program is SAMPLE_S ****************
c** Date: 20/01/98   **
c** author: A.Muecke **
c**********************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)

        SAVE

       common/input/ tbb,E0,alpha1,alpha2,
     &           epsm1,epsm2,epsb,L0

        external crossection
        double precision crossection


        pm = 0.93827D0
        factor = (s-pm*pm)
        epsprime = factor/2.D0/pm
        sigma_pg = crossection(epsprime,3,L0)
        functs = factor*sigma_pg 

        RETURN

        END

	DOUBLE PRECISION FUNCTION PHOTD(EPS,TBB)
C **************************************************************
C    RETURNS THE DENSITY OF BLACKBODY RADIATION OF TEMPERATURE *
C "TBB" DEGREES (DENS1). EPS IN eV, PHOTD IN #/(cm^3.eV)       *
C                                    TSS,  May '92             *
C **************************************************************
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        SAVE
C CONVERT TEMPERATURE TO eV
	EPH = EPS
        EKT = 8.619D-5*TBB
        EPHKT = EPS/EKT
        IF (EPHKT .GT. 80.D0)     GO TO 10
        IF (EPHKT .LT. 1.D-4)   GO TO 11
        FEE = DEXP(EPHKT) - 1.D0
        GO TO 12
   11   FEE = EPHKT
   12   BB = 1.318D13*EPH*EPH/FEE
	GO TO 15
   10   BB = 0.D0
   15	PHOTD = BB
	END


       DOUBLE PRECISION function plinterpol(alpha)

c*** interpolates p(Ep) to give the max. probability p(eps) for ***
c*** a given initial proton energy                              ***
c** Date: 20/01/98   **
c** author: A.Muecke **
c**********************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       SAVE

       DIMENSION AINDEX(4), A(4)
       DATA A / 0.D0,1.D0,2.D0,3.D0/
       DATA AINDEX / 8.D8,5.D8,5.D8,1.D9/

       plinterpol = 0.D0

       do 1 ni=1,3
        p1 = A(ni)
        p2 = A(ni+1)
        if (alpha.le.p2.and.alpha.gt.p1) then
         tang = (log10(AINDEX(ni+1))-log10(AINDEX(ni)))/(p2-p1)
         plinterpol = log10(AINDEX(ni))+(alpha-p1)*tang
         plinterpol = 10.D0**plinterpol
        endif
  1    continue

       if (alpha.eq.0.D0) plinterpol = 5.D8

       if (plinterpol.eq.0.D0) then
        print*,'interpolation not sucessful !'
        pause
       endif

       END

       DOUBLE PRECISION function f_epspl(eps)

c*** gives energy density law of power law photon field    ***
c*** f(epsilon) = eps^-alpha, eps=[epsm1,epsm2], eps in eV *************
c** Date: 14/03/99   **
c** author: A.Muecke **
c**********************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)

       SAVE

       common/input/ tbb,E0,alpha1,alpha2,
     &           epsm1,epsm2,epsb,L0


       alpha12 = alpha2-alpha1
       ampl = epsb**alpha12
       if (eps.lt.epsb) then
         f_epspl = eps*(eps**(-alpha1))
       else
         f_epspl = eps*ampl*(eps**(-alpha2))
       endif

       RETURN

       END


