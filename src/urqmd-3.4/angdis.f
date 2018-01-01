c $Id: angdis.f,v 1.22 2007/05/23 14:28:50 bleicher Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine angdisnew(sqrts,m1,m2,iline,costh,phi)
c
c     Revision : 1.1
c
c input:  sqrts, m1, m2, iline : characteristics of the ingoing channel 
coutput costh   : cos(theta) of theta-angle
coutput phi     : phi-angle
c
c     {\tt angdisnew} delivers phi and cos(theta) scattering angles
c     according to the angular distributions given by Mao et al.
c     {\tt angdisnew} performs numerical inversion of the integral of the 
c     differential cross-section by means of a bisection method.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      real*8 sqrts, m1, m2, costh, costhcoll, phi, ranf, pi, s
      real*8 anginter
      integer iline
      logical symlog(64)

      integer pythflag
      common /pythflag/ pythflag
    
      parameter (pi = 3.14159265358979312d0)

c symmetrize or not angular distribution (depending on iline)
c this data statements may be changed ...
      data symlog /14*.true., 1*.false., 6*.true., 7*.false.,
     &              7*.true., 1*.false., 2*.true., 1*.true.,
     &              25*.true./


      s = sqrts*sqrts
      phi = 2.0d0*pi*ranf(0)
chp fix angular distribution, Pythia rotates strings already
      if(pythflag.eq.1)go to 10

      goto (4, 4, 4, 4, 4, 4, 4, 4, 4, 9,
     &      9, 4, 4, 4, 4, 4, 4, 4, 4, 9,
c ISO-FB interpolation for MB iline 26/27/28
     &      4, 4,10,10,10,11,11,11, 4, 4,
     &      4, 4, 4, 4, 4, 9, 9, 4, 9, 4,
     &      4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
     &      4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
     &      4, 4, 4, 4) abs(iline)

c cross-sections NN (Mao et al.)
   4  continue
c.. quick and dirty fix for very high energies
c.. i.e. at high energies the (elastic) angular distribution is
c.. forward-backward peaked
      if (sqrts.gt.500d0) goto 10
      costh = -costhcoll(s,m1,m2,symlog(abs(iline)))

      return

c isotropic decay
   9  continue
      costh = 1.0d0-2.0d0*ranf(0)
      return

c no deflection at all
  10  continue 
      costh = 1.0d0
      return

c smooth interpolation between iso and f-b:
  11  continue
cbl only for intermediate masses, otherwise no deflection (zero degree scattering)
      if(sqrts.gt.6d0) goto 10
      costh = anginter()
      return

      end
  
      function costhcoll(s,m1,m2,sym)
      implicit none 
      real*8 costhcoll, s, m1, m2, x, dct, ct, dsigma, ranf
      integer j
      logical sym

      x = ranf(0)
      dct = 2.0d0
      costhcoll = -1.0d0
c
c for jmax=12 the accuracy is better than 0.1 degree
c
      do j=1,12 ! accuracy 2**-jmax
        dct = 0.5d0*dct
        ct = costhcoll+dct
        if (dsigma(s,m1,m2,sym,ct).le.x) costhcoll=ct
      enddo
c
c randomize in final interval in order to avoid discrete angles
c
      costhcoll = costhcoll+ranf(0)*dct
      return
      end

      function dsigma(s_in,m1_in,m2_in,sym,costh)
c
cc dsigma(s\_in,chosth) = int_-1^costh dsigma/dOmega(s_in,..) dOmega
cc it is normalized such that dsigma(s\_in,-1) = 0 and
cc                            dsigma(s\_in,1) = 1
      implicit none

      include 'coms.f'

      real*8 dsigma
      real*8 s_in, m1_in, m2_in, costh,
     &       msi, cmsi, gsi, mom, cmom, gom, mpi, cmpi, gpi, m
      real*8 m42, mpi2, cmpi2, d_pi1, d_pi2, cm6gp,
     &       cpi_3, cpi_2, cpi_1, cpi_m, cpi_l, cpi_0,
     &       msi2, cmsi2, cmsi4, cmsi6, d_si1, d_si2, d_si3, cm2gs,
     &       csi_3, csi_2, csi_1, csi_m, csi_l, csi_0,
     &       mom2, cmom2, cmom4, cmom6, d_om1, d_om2, d_om3, s_om1,
     &       cm2go, com_3, com_2, com_1, com_m, com_l,
     &       fac1, d_mx1, d_mx2, d_mx3,
     &       cmx_o1, cmx_s1, cmx_om, cmx_sm, fac2, fac3,
     &       cmx_olc, cmx_ols, cmx_slc, cmx_sls

      real*8 sig, tp_pi, tp_si, tp_om, tm_pi, tm_si, tm_om,
     &       bom_3, bom_2, bom_1, bom_m, bom_0, bom_l,
     &       bmx_o1, bmx_s1, bmx_om, bmx_sm, bmx_ol, bmx_sl

      real*8 s, tmax, tp, to, twos, brak1, norm, 
     &       t1_pi, t2_pi, t1_si, t2_si, t1_om, t2_om,
     &       t3_pi, t4_pi, t3_si, t4_si, t3_om, t4_om

      logical firstlog, sym
      common /nn/ norm


c define masses and coupling constants
      data msi /0.550d0/ cmsi /1.200d0/ gsi  /9.40d0/
      data mom /0.783d0/ cmom /0.808d0/ gom /10.95d0/
      data mpi /0.138d0/ cmpi /0.510d0/ gpi  /7.27d0/
      data m   /0.938d0/        

      save firstlog

      save m42, mpi2, cmpi2, d_pi1, d_pi2, cm6gp,
     &     cpi_3, cpi_2, cpi_1, cpi_m, cpi_l, cpi_0

      save msi2, cmsi2, cmsi4, cmsi6, d_si1, d_si2, d_si3, cm2gs,
     &     csi_3, csi_2, csi_1, csi_m, csi_l, csi_0

      save mom2, cmom2, cmom4, cmom6, d_om1, d_om2, d_om3, s_om1,
     &     cm2go, com_3, com_2, com_1, com_m, com_l 

      save fac1, d_mx1, d_mx2, d_mx3,
     &     cmx_o1, cmx_s1, cmx_om, cmx_sm, fac2, fac3,
     &     cmx_olc, cmx_ols, cmx_slc, cmx_sls
      
      sig(tp_pi,tp_si,tp_om,tm_pi,tm_si,tm_om) =
c pion
     & +((cpi_3*tp_pi + cpi_2)*tp_pi + cpi_1)*tp_pi 
     & + cpi_m/tm_pi + cpi_0 + cpi_l*log(tp_pi*tm_pi)
c sigma
     & +((csi_3*tp_si + csi_2)*tp_si + csi_1)*tp_si 
     & + csi_m/tm_si + csi_0 + csi_l*log(tp_si*tm_si)
c omega 
     & +((bom_3*tp_om + bom_2)*tp_om + bom_1)*tp_om 
     & + bom_m/tm_om + bom_0 + bom_l*log(tp_om*tm_om)
c mix
     & + bmx_o1*(tp_om - 1.0d0)
     & + bmx_s1*(tp_si - 1.0d0)
     & + bmx_om*log(tm_om)
     & + bmx_sm*log(tm_si)
     & + bmx_ol*log(tp_om)
     & + bmx_sl*log(tp_si)

c calculate constants only once!
      if(firstlog) goto 1000
      if (info) write(6,*) 
     $     '(info) dsigma: calculating constants for ang. dist.'

c define constants for pion-Term (no s-dependence)
      m42   = 4.0d0*m*m
      mpi2  = mpi*mpi           
      cmpi2 = cmpi*cmpi
      d_pi1 = cmpi2-mpi2
      d_pi2 = d_pi1*d_pi1
      cm6gp = 1.5d0*cmpi2**3*gpi**4*m42*m42/d_pi2

      cpi_3 = -(cm6gp/3.0d0)
      cpi_2 = -(cm6gp*mpi2/d_pi1)
      cpi_1 = -(cm6gp*mpi2*(2.0d0*cmpi2 + mpi2)/d_pi2)
      cpi_m = -(cm6gp*cmpi2*mpi2/d_pi2)
      cpi_l = -(cm6gp*2.0d0*cmpi2*mpi2*(cmpi2 + mpi2)/d_pi2/d_pi1)
      cpi_0 = -(cpi_3 + cpi_2 + cpi_1 + cpi_m)

c define constants for sigma-Term (no s-dependence)
      msi2  = msi*msi
      cmsi2 = cmsi*cmsi
      cmsi4 = cmsi2*cmsi2
      cmsi6 = cmsi2*cmsi4
      d_si1 = m42-cmsi2
      d_si2 = m42-msi2
      d_si3 = cmsi2-msi2
      cm2gs = 0.5d0*cmsi2*gsi**4/d_si3**2      

      csi_3 = -(cm2gs*d_si1**2/3.0d0)
      csi_2 = -(cm2gs*cmsi2*d_si1*d_si2/d_si3)
      csi_1 = -(cm2gs*cmsi4*(2.0d0*d_si1 + d_si2)*d_si2/d_si3**2)
      csi_m = -(cm2gs*cmsi6*d_si2**2/msi2/d_si3**2)
      csi_l = -(cm2gs*cmsi6*d_si2*(d_si1 + d_si2)*2.0d0/d_si3**3)
      csi_0 = -(csi_3 + csi_2 + csi_1 + csi_m)

c define constants for omega-Term
      mom2  = mom*mom
      cmom2 = cmom*cmom
      cmom4 = cmom2*cmom2
      cmom6 = cmom2*cmom4
      d_om1 = m42-cmom2
      d_om2 = m42-mom2
      d_om3 = cmom2-mom2
      s_om1 = cmom2+mom2
      cm2go = 0.5d0*cmom2*gom**4/d_om3**2

      com_3 =  cm2go/3.0d0 
      com_2 = -(cm2go*cmom2/d_om3)
      com_1 =  cm2go*cmom4/d_om3**2
      com_m =  cm2go*cmom6/(d_om3**2*mom2)
      com_l = -(cm2go*cmom6*4.0d0/d_om3**3)

c define constants for mix-Term
      fac1 = -((gsi*gom*cmsi2*cmom2)**2*m42)  
      d_mx1 = cmom2 - cmsi2
      d_mx2 = cmom2 - msi2
      d_mx3 = cmsi2 - mom2

      cmx_o1 = fac1/(cmom2*d_mx1**2*d_mx2*d_om3)
      cmx_s1 = fac1/(cmsi2*d_mx1**2*d_mx3*d_si3)
      cmx_om = fac1/(d_om3**2*d_mx3**2*(mom2 - msi2))
      cmx_sm = fac1/(d_si3**2*d_mx2**2*(msi2 - mom2)) 
      fac2 = (-fac1)/(d_mx1**3*d_om3**2*d_mx2**2)
      fac3 = (-fac1)/(d_mx1**3*d_mx3**2*d_si3**2) 
      
      cmx_olc =
     &  fac2*(3.0d0*cmom2**3        - cmom2**2*cmsi2 
     &      - 2.0d0*cmom2**2*mom2   - 2.0d0*cmom2**2*msi2 
     &      + cmom2*mom2*msi2       + cmsi2*mom2*msi2 
     &      - 4.0d0*cmom2**2*m42    + 2.0d0*cmom2*cmsi2*m42 
     &      + 3.0d0*cmom2*mom2*m42  - cmsi2*mom2*m42 
     &      + 3.0d0*cmom2*msi2*m42  - cmsi2*msi2*m42 
     &      - 2.0d0*mom2*msi2*m42)

      cmx_ols = 
     &  fac2*(8.0d0*cmom2**2        - 4.0d0*cmom2*cmsi2 
     &      - 6.0d0*cmom2*mom2      + 2.0d0*cmsi2*mom2 
     &      - 6.0d0*cmom2*msi2      + 2.0d0*cmsi2*msi2 
     &      + 4.0d0*mom2*msi2)

      cmx_slc = 
     &  fac3*(cmom2*cmsi2**2        - 3.0d0*cmsi2**3 
     &      + 2.0d0*cmsi2**2*mom2   + 2.0d0*cmsi2**2*msi2 
     &      - cmom2*mom2*msi2       - cmsi2*mom2*msi2 
     &      - 2.0d0*cmom2*cmsi2*m42 + 4.0d0*cmsi2**2*m42 
     &      + cmom2*mom2*m42        - 3.0d0*cmsi2*mom2*m42 
     &      + cmom2*msi2*m42        - 3.0d0*cmsi2*msi2*m42 
     &      + 2.0d0*mom2*msi2*m42) 

      cmx_sls = 
     &  fac3*(4.0d0*cmom2*cmsi2     - 8.0d0*cmsi2**2 
     &      - 2.0d0*cmom2*mom2      + 6.0d0*cmsi2*mom2 
     &      - 2.0d0*cmom2*msi2      + 6.0d0*cmsi2*msi2 
     &      - 4.0d0*mom2*msi2)

      firstlog = .true.
      if (info) write(6,*) '(info) dsigma: calculation finished'

c s-dependence beyond this point

 1000 continue
      s = s_in - (m1_in+m2_in)**2 + m42
      tmax = s-m42
      tp = 0.5d0*(costh+1.0d0)*tmax
      twos = 2.0d0*s

c define s-dependent stuff for omega-Term
      brak1 = (twos-m42)**2
      bom_3 = com_3*(-(2.0d0*cmom2**2) - 2.0d0*cmom2*twos - brak1)
      bom_2 = com_2*(2.0d0*cmom2*mom2 + s_om1*twos + brak1)
      bom_1 = com_1*(-(4.0d0*cmom2*mom2) - 2.0d0*mom2**2 - 
     &               2.0d0*(cmom2+2*mom2)*twos - 3.0d0*brak1)
      bom_m = com_m*(-(2.0d0*mom2**2)- 2.0d0*mom2*twos - brak1)
      bom_l = com_l*(s_om1*mom2 + (cmom2 + 3.0d0*mom2)*s + brak1)
      bom_0 = -(bom_3 + bom_2 + bom_1 + bom_m)

c define s-dependent stuff for mix-Term            
      bmx_o1 = cmx_o1*(d_om1 - twos)
      bmx_s1 = cmx_s1*(d_si1 - twos)
      bmx_om = cmx_om*(d_om2 - twos)
      bmx_sm = cmx_sm*(d_si2 - twos)
      bmx_ol = cmx_olc + cmx_ols*s
      bmx_sl = cmx_slc + cmx_sls*s

      t1_pi = 1.0d0/(1.0d0+tmax/cmpi2)
      t2_pi = 1.0d0+tmax/mpi2
      t1_si = 1.0d0/(1.0d0+tmax/cmsi2)
      t2_si = 1.0d0+tmax/msi2
      t1_om = 1.0d0/(1.0d0+tmax/cmom2)
      t2_om = 1.0d0+tmax/mom2

      norm = sig(t1_pi,t1_si,t1_om,t2_pi,t2_si,t2_om)
      
      t1_pi = 1.0d0/(1.0d0+tp/cmpi2)
      t2_pi = 1.0d0+tp/mpi2
      t1_si = 1.0d0/(1.0d0+tp/cmsi2)
      t2_si = 1.0d0+tp/msi2
      t1_om = 1.0d0/(1.0d0+tp/cmom2)
      t2_om = 1.0d0+tp/mom2

      if (sym) then 
        norm=2.0d0*norm
        to = tmax-tp
        t3_pi = 1.0d0/(1.0d0+to/cmpi2)
        t4_pi = 1.0d0+to/mpi2
        t3_si = 1.0d0/(1.0d0+to/cmsi2)
        t4_si = 1.0d0+to/msi2
        t3_om = 1.0d0/(1.0d0+to/cmom2)
        t4_om = 1.0d0+to/mom2

        dsigma = (sig(t1_pi,t1_si,t1_om,t2_pi,t2_si,t2_om)
     &           -sig(t3_pi,t3_si,t3_om,t4_pi,t4_si,t4_om))/norm+0.5d0
      else
        dsigma = sig(t1_pi,t1_si,t1_om,t2_pi,t2_si,t2_om)/norm
      end if
      return
      end

      function anginter()
      implicit none
      real*8 ranf, anginter, a
c*
      a=8d0
c*
c the costheta distribution p(x) is proportional to exp(a*x)
c x is chosen between -1 and +1
c a=0  corresponds to isotropic distr.
c a=infinity corresponds to exactly forward distr.      
c inverse transform method is used:

      anginter=1d0/a*log( ranf(0)*(exp(a)-exp(-a))+exp(-a) )

      if(anginter.lt.-1d0.or.anginter.gt.1d0)then
       write(*,*)"#angdis# illegal costh value: ",anginter
      endif

      return
      end

