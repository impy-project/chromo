c $Id: comstr.f,v 1.2 1999/01/18 09:56:58 ernst Exp $

      parameter (njspin=8)
  
      real*8  PJSPNS,PMIX1S(3,njspin),PMIX2S(3,njspin),PBARS,
     *PARQLS,PARRS
      real*8  PJSPNC,PMIX1C(3,njspin),PMIX2C(3,njspin),PBARC
 

    
      COMMON/FRGSPA/ PJSPNS,PMIX1S,PMIX2S,PBARS,
     *PARQLS,PARRS
      COMMON/FRGCPA/ PJSPNC,PMIX1C,PMIX2C,PBARC
 
c parm gives the probability for different meson multiplets according
c to spin degeneracy and average mass ratios
c spin-parity 0- : 1- : 0+ : 1+ : 2+ = parm(1):parm(2)...:parm(njspin)  
      real*8 parm(njspin)

      common/coparm/parm

      real*8 pi
      COMMON/CONST/ PI
