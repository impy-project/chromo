!----------------------------------------------------------------------!        
!     Wrappers for Fluka's functions 
!----------------------------------------------------------------------!


      subroutine icode_from_pdg(pdg_id, icode_id)
!----------------------------------------------------------------------!        
!     get internal fluka code of particle type from pdg id  
!----------------------------------------------------------------------!          
      integer pdg_id, icode_id
Cf2py intent(out) icode_id             
      icode_id = MCIHAD(pdg_id)
      end subroutine icode_from_pdg 

      subroutine icode_from_pdg_arr(npart, pdg_id, icode_id)
!----------------------------------------------------------------------!        
!     get internal fluka code of particle type from pdg id  
!----------------------------------------------------------------------!
      integer :: npart                
      integer :: pdg_id(npart), icode_id(npart)
      integer :: i
Cf2py intent(out) icode_id
Cf2py integer intent(hide),depend(pdg_id) :: npart=len(pdg_id)
      do i=1,npart            
            icode_id(i) = MCIHAD(pdg_id(i))
      end do      
      end subroutine icode_from_pdg_arr
      
      subroutine charge_from_pdg_arr(npart, pdg_id, charge)
!----------------------------------------------------------------------!        
!     get internal fluka code of particle type from pdg id  
!----------------------------------------------------------------------!
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(PAPROP)'
      INCLUDE '(PART2)'
      INCLUDE '(NUCFLG)'

      integer :: npart                
      integer :: pdg_id(npart), charge(npart)
      integer :: i
Cf2py intent(out) charge
Cf2py integer intent(hide),depend(pdg_id) :: npart=len(pdg_id)
      do i=1,npart            
            charge(i) = IICH(MCIHAD(pdg_id(i)))
      end do     
      end subroutine charge_from_pdg_arr 

    
      subroutine pdg_from_icode(icode_id, pdg_id)
!----------------------------------------------------------------------!        
!     pdg id from internal particle code 
!----------------------------------------------------------------------!               
      integer icode_id, pdg_id
Cf2py intent(out) pdg_id
      pdg_id = -777
      if ((icode_id.ge.1).and.(icode_id.le.390)) then     
        pdg_id = MPDGHA(icode_id)  
      end if      
      end subroutine pdg_from_icode

      subroutine random_direction(dir_cos)
!----------------------------------------------------------------------! 
!     Wrapper for fluka function SPRNCS
!     returning cosines of direction
!----------------------------------------------------------------------!      
      double precision dir_cos(3)
Cf2py intent(out) dir_cos                  
      CALL SPRNCS (dir_cos(1), dir_cos(2), dir_cos(3))
      end subroutine random_direction
    
      INTEGER FUNCTION ICRVRCK()
!----------------------------------------------------------------------! 
!     Returns a key for correct Fluka version
!     The key is required by some functions
!----------------------------------------------------------------------!          
      INCLUDE '(DBLPRC)'
      ICRVRCK = NINT(AMPRMU * 1.D+12 - 1.0072D+12)
      end function ICRVRCK
      

      subroutine init_rng_state(file_name, logical_unit, 
     &    seed, ntot, ntot2)
!----------------------------------------------------------------------! 
!       Initiate and write a state for "Ranmar" generator. 
!       Wrapper for RNINIT
!
!                    
!        Fluka uses "Ranmar" by Marsaglia and Zaman 
!        of Florida State University. Probably it is the same as
!        "rmmard". 
!         
!        0 <= seed < 2e9 defines a sequence
!        ntot  and ntot2 skeeps (ntot2*m+ntot) numbers
!        In "rmmard" m=1e9, Fluka probably uses the same number
!     
!        BE CAREFULL when setting ntot and ntot2, as the generator
!        just calculates all numbers in a sequence to reach required 
!        state. It can take significant TIME       
!----------------------------------------------------------------------!         
      integer logical_unit, seed, ntot, ntot2
      character(300) file_name

      open(unit=logical_unit, file=trim(file_name))
      call RNINIT(logical_unit, seed, ntot, ntot2)
      close(unit=logical_unit)
      
      end subroutine init_rng_state 


      subroutine load_rng_state(file_name, logical_unit)
!----------------------------------------------------------------------! 
!       Loads state of fluka's random number generator
!----------------------------------------------------------------------!          
        integer logical_unit, seed
        character(300) file_name
        logical success_flag 
        ! seed = -1: reads seed from file
        seed = -1
        open(unit=logical_unit, file=trim(file_name))
        call RNREAD(logical_unit, seed, success_flag)
        close(unit=logical_unit)
      end subroutine load_rng_state
      
      
      subroutine save_rng_state(file_name, logical_unit)
!----------------------------------------------------------------------! 
!       Saves state of fluka's random number generator
!----------------------------------------------------------------------!         
        integer logical_unit
        character(300) file_name

        open(unit=logical_unit, file=trim(file_name))
        call RNWRIT(logical_unit) 
        close(unit=logical_unit)
      end subroutine save_rng_state  


      subroutine fluka_rand(random_number)
!----------------------------------------------------------------------! 
!       Wrapper for fluka's random number generator
!----------------------------------------------------------------------!         
        double precision random_number, FLRNDM
Cf2py intent(out) random_number
        random_number = FLRNDM(0d0)
      end subroutine fluka_rand


      subroutine fluka_particle_scheme
!----------------------------------------------------------------------!
!   Prints particle scheme used by FLUKA
!   Also required to expose some common blocks from FLUKA     
!
!   Inspired by PDGFLK program: 
!   Copyright (C) 2023-2023      by    Alfredo Ferrari & Paola Sala           
!----------------------------------------------------------------------!
        INCLUDE '(DBLPRC)'
        INCLUDE '(DIMPAR)'
        INCLUDE '(PAPROP)'
        INCLUDE '(PART2)'

        WRITE(*, 2500) "Internal", "External", "PDG", 
     &  "NAME", "SHORT", "MASS", "CHARGE", "BARYON"
2500    FORMAT (2A10, A12, A15, 4A10)
2600    FORMAT (I10, I10, I12, A15, A10, F10.4, I10, I10)      
        DO I = -6, 390
           !icode_id = IPTOKP(extcode_id)
           !extcode_id = KPTOIP(icode_id) 
           KPFLK = I
           IPDG  = MPDGHA (KPFLK)
           IPFLK = KPTOIP (KPFLK)
           WRITE (*,2600) KPFLK, IPFLK, IPDG, PRNAME(IPFLK),
     &     ANAME(KPFLK), AAM(KPFLK), IICH(KPFLK), IIBAR(KPFLK)  
        END DO   
      end subroutine fluka_particle_scheme

      DOUBLE PRECISION FUNCTION CHROMO_SGMXYZ ( KPROJ0, MMAT  , EKIN0 ,
     &                                   PPROJ0, IFLXYZ )

*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2022-2023      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*     SiGMa (mb) for X/Y/Z interaction types:                          *
*                                                                      *
*     Authors:                           Alfredo Ferrari & Paola Sala  *
*                                                                      *
*                                                                      *
*     Created on  24 October 2022  by    Alfredo Ferrari & Paola Sala  *
*                                            Private        Private    *
*                                                                      *
*     Last change on  11-May-23    by             Alfredo Ferrari      *
*                                                     Private          *
*                                                                      *
*     Input variables:                                                 *
*                                                                      *
*           Kproj0 = Particle id (Fluka external "Paprop" numbering)   *
*                    A heavy ion projectile can be input using the     *
*                    following (pdg) coding:                           *
*                         Kproj0 = M + A x 10 + Z x 10000              *
*                                + H x 10000000 + 1000000000           *
*           Mmat   = Material number                                   *
*           Ekin0  = Particle kinetic energy (GeV) if > 0, if =< 0 the *
*                    kinetic energy is reconstructed out of Kproj/Pproj*
*           Pproj0 = Particle momentum (GeV/c) if > 0, if =< 0 the     *
*                    momentum is reconstructed out of Kproj/Ekin       *
*                    In case both Ekin0/Pproj0 are > 0, the momentum   *
*                    takes precedence and the kinetic energy is re-    *
*                    constructed from the momentum using the Fluka mass*
*           Iflxyz = Flag defining which processes the cross section   *
*                    is asked for                                      *
*                    Iflxyz =  1 -> only inelastic                     *
*                    Iflxyz = 10 -> only elastic                       *
*                    Iflxyz = 11 -> inelastic + elastic                *
*                    Iflxyz =100 -> only emd                           *
*                    Iflxyz =101 -> inelastic + emd                    *
*                    Iflxyz =110 -> elastic + emd                      *
*                    Iflxyz =111 -> inelastic + elastic + emd          *
*                                                                      *
*     Output variable:                                                 *
*                                                                      *
*           Sgmxyz = Total cross section (mb) for the requested pro-   *
*                    cesses                                            *
*                                                                      *
*----------------------------------------------------------------------*
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
      INCLUDE '(PAPROP)'

      CHROMO_SGMXY = SGMXYZ( KPROJ0, MMAT  , EKIN0 , PPROJ0, IFLXYZ )
      RETURN
      END FUNCTION CHROMO_SGMXYZ

      SUBROUTINE CHROMO_EVTXYZ ( KPROJ0, MMAT  , EKIN0 , PPROJ0, TXX,
     &                    TYY, TZZ, IFLXYZ, CUMSGI, CUMSGE, CUMSGM )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2022-2023      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     EVenT for X/Y/Z interaction types:                               *
*                                                                      *
*     Authors:                           Alfredo Ferrari & Paola Sala  *
*                                                                      *
*                                                                      *
*     Created on  24 October 2022  by    Alfredo Ferrari & Paola Sala  *
*                                            Private        Private    *
*                                                                      *
*     Last change on  18-Jul-23    by             Alfredo Ferrari      *
*                                                     Private          *
*                                                                      *
*           Kproj0 = Particle id (Fluka external "Paprop" numbering)   *
*                    A heavy ion projectile can be input using the     *
*                    following (pdg) coding:                           *
*                         Kproj0 = M + A x 10 + Z x 10000              *
*                                + H x 10000000 + 1000000000           *
*           Mmat   = Material number of the target                     *
*           Ekin0  = Particle kinetic energy (GeV) if > 0, if =< 0 the *
*                    kinetic energy is reconstructed out of Kproj/Pproj*
*           Pproj0 = Particle momentum (GeV/c) if > 0, if =< 0 the     *
*                    momentum is reconstructed out of Kproj/Ekin       *
*                    In case both Ekin0/Pproj0 are > 0, the momentum   *
*                    takes precedence and the kinetic energy is re-    *
*                    constructed from the momentum using the Fluka mass*
*        Txx/yy/zz = Particle direction cosines                        *
*           Iflxyz = Flag defining which processes the cross section   *
*                    is asked for                                      *
*                    Iflxyz =  1 -> only inelastic                     *
*                    Iflxyz = 10 -> only elastic                       *
*                    Iflxyz = 11 -> inelastic + elastic                *
*                    Iflxyz =100 -> only emd                           *
*                    Iflxyz =101 -> inelastic + emd                    *
*                    Iflxyz =110 -> elastic + emd                      *
*                    Iflxyz =111 -> inelastic + elastic + emd          *
*                                                                      *
*     Output variables:                                                *
*                                                                      *
*           Cumsgi(i) = cumulative macroscopic inelastic cross section *
*                       (cm^2/g x rho) for the i_th element of the mmat*
*                       material                                       *
*           Cumsge(i) = cumulative macroscopic elastic   cross section *
*                       (cm^2/g x rho) for the i_th element of the mmat*
*                       material                                       *
*           Cumsgm(i) = cumulative macroscopic emd       cross section *
*                       (cm^2/g x rho) for the i_th element of the mmat*
*                       material                                       *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*

      PARAMETER ( AVOGMB = AVOGAD * 1.D-27 )
*
      INCLUDE '(BALANC)'
      INCLUDE '(CMELDS)'
      INCLUDE '(CRSKCM)'
      INCLUDE '(EVTFLG)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKCMP)'
      INCLUDE '(FLKMAT)'
      INCLUDE '(GENSTK)'
      INCLUDE '(PAPROP)'
      INCLUDE '(PAREVT)'
      INCLUDE '(RESNUC)'
      INCLUDE '(THRSCM)'
      DIMENSION CUMSGI (0:NELEMX), CUMSGE (0:NELEMX), CUMSGM (0:NELEMX)
Cf2py intent(out) CUMSGI, CUMSGE, CUMSGM  
      CALL EVTXYZ ( KPROJ0, MMAT  , EKIN0 , PPROJ0, TXX, TYY, TZZ,
     &                    IFLXYZ, CUMSGI, CUMSGE, CUMSGM )

      RETURN
      END SUBROUTINE CHROMO_EVTXYZ

      SUBROUTINE CHROMO_FLLHEP

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2023-2023      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*     FiLL HEP common:                                                 *
*                                                                      *
*     Authors:                           Alfredo Ferrari & Paola Sala  *
*                                                                      *
*                                                                      *
*     Created on 06 February 2023  by    Alfredo Ferrari & Paola Sala  *
*                                            Private        Private    *
*                                                                      *
*     Last change on  09-Mar-23    by             Alfredo Ferrari      *
*                                                     Private          *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(BALANC)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(GENSTK)'
      INCLUDE '(HEPCMM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(PART2)'
      INCLUDE '(RESNUC)'
      INCLUDE '(THRSCM)'

      CALL FLLHEP
      RETURN
      END SUBROUTINE

      SUBROUTINE CHROMO_STPXYZ ( NMATFL, NELMFL, IZELFL, WFELFL, MXELFL,
     &                    PPTMAX, EF2DP3, DF2DP3, IFLXYZ, LPRINT,
     &                    MTFLKA )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2022-2023      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     SeTuP for X/Y/Z interaction types:                               *
*                                                                      *
*     Authors:                           Alfredo Ferrari & Paola Sala  *
*                                                                      *
*                                                                      *
*     Created on  24 October 2022  by    Alfredo Ferrari & Paola Sala  *
*                                            Private        Private    *
*                                                                      *
*     Last change on  29-Mar-23    by             Alfredo Ferrari      *
*                                                     Private          *
*                                                                      *
*     Input variables:                                                 *
*                                                                      *
*           Nmatfl    = Number of requested Fluka materials            *
*           Nelmfl(m) = Number of elements in the m_th requested Fluka *
*                       material, Nelmfl(m)=1 means simple element, no *
*                       compound                                       *
*           Izelfl(n) = Cumulative array of the Z of the elements con- *
*                       stituting the requested Fluka materials        *
*                      (the Z for the elements of the m_th materials   *
*                       start at n = Sum_i=1^m-1 [Nelmfl(i)] + 1)      *
*           Wfelml(n) = Cumulative array of the weight fractions of the*
*                       elements constituting the requested Fluka      *
*                       materials (the weight fractions for the elem-  *
*                       ents of the m_th materials start at            *
*                        n = Sum_i=1^m-1 [Nelmfl(i)] + 1)              *
*           Mxelfl    = Dimension of the Iz/Wfelml arrays, it must be  *
*                       Mxelfl >= Sum_i=1^nmatfl [Nelmfl(i)]           *
*           Pptmax    = Maximum momentum (GeV/c) to be used in initia- *
*                       lization (optional)                            *
*           Ef2dp3    = Transition energy (GeV) from Fluka (Peanut) to *
*                       Dpmjet3 for hA interactions (optional)         *
*           Df2dp3    = Smearing (+/-Df2dp3) (GeV) of the transition   *
*                       energy from Fluka (Peanut) to Dpmjet3 for hA   *
*                       interactions (optional)                        *
*           Lprint    = Material printout flag                         *
*                                                                      *
*     Output variables:                                                *
*                                                                      *
*           Mtflka(m) = Fluka material number corresponding to the m_th*
*                       requested material                             *
*                                                                      *
*----------------------------------------------------------------------*
*
  
      INCLUDE '(BEAMCM)'
      INCLUDE '(CLSCCM)'
      INCLUDE '(CTITLE)'
      INCLUDE '(CMELDS)'
      INCLUDE '(CRSKCM)'
      INCLUDE '(EVAFLG)'
      INCLUDE '(FLKCMP)'
      INCLUDE '(FLKMAT)'
      INCLUDE '(GENFLG)'
      INCLUDE '(INFLEX)'
      INCLUDE '(NUCDAT)'
      INCLUDE '(NUCGEO)'
      INCLUDE '(PAPROP)'
      INCLUDE '(PAREVT)'
      INCLUDE '(PHNCCM)'
      INCLUDE '(QELCMM)'
      INCLUDE '(RESNUC)'
      INCLUDE '(STNHCM)'
      INCLUDE '(THRSCM)'
      DIMENSION NELMFL (NMATFL), MTFLKA (NMATFL), IZELFL (MXELFL),
     &          WFELFL (MXELFL)
Cf2py intent(out) MTFLKA
      CHARACTER*8 CRVRCK
      INTEGER IKEY

      IKEY = ICRVRCK()
      WRITE (CRVRCK,'(I8)') IKEY
      CALL STPXYZ ( NMATFL, NELMFL, IZELFL, WFELFL, MXELFL,
     &                    PPTMAX, EF2DP3, DF2DP3, IFLXYZ, LPRINT,
     &                    MTFLKA, CRVRCK )
      RETURN

      END SUBROUTINE