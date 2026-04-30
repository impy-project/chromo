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
!     get electric charge of particle from pdg id
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


      subroutine load_rng_state(file_name, logical_unit, status)
!----------------------------------------------------------------------!
!       Loads state of fluka's random number generator.
!       Returns status: 0 = OK, 1 = open failed, 2 = read failed.
!----------------------------------------------------------------------!
        integer logical_unit, seed, status, io_status
Cf2py intent(out) status
        character(300) file_name
        logical success_flag
        seed = -1
        status = 0
        open(unit=logical_unit, file=trim(file_name),
     &       iostat=io_status)
        if (io_status .ne. 0) then
          status = 1
          return
        end if
        call RNREAD(logical_unit, seed, success_flag)
        close(unit=logical_unit)
        if (.not. success_flag) then
          status = 2
        end if
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
!   The GENTHR include pulls in FLUKA's hadronic-generator transition
!   thresholds (Peanut <-> DPMJET for hA, rQMD <-> DPMJET for AA,
!   DPMJET <-> UHE in both) so chromo's Python wrapper can override
!   them via f2py-exposed common-block assignments.
!
!   Inspired by PDGFLK program:
!   Copyright (C) 2023-2023      by    Alfredo Ferrari & Paola Sala
!----------------------------------------------------------------------!
        INCLUDE '(DBLPRC)'
        INCLUDE '(DIMPAR)'
        INCLUDE '(PAPROP)'
        INCLUDE '(PART2)'
        INCLUDE '(GENTHR)'

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

      CHROMO_SGMXYZ = SGMXYZ( KPROJ0, MMAT  , EKIN0 , PPROJ0, IFLXYZ )
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
      SUBROUTINE pdg_to_proj_code(pdg_id, proj_code)
!----------------------------------------------------------------------!
!     Convert a PDG particle id to the FLUKA external "PAPROP" code
!     expected by SGMXYZ. Handles hadrons, the photon, and nuclei.
!
!     Hadron/lepton:   proj_code = KPTOIP(MCIHAD(pdg_id))
!     Photon (pdg=22): proj_code = 7
!     Light nuclei (d, t, 3He, 4He): dedicated native internal codes
!                      -3, -4, -5, -6 respectively (FLUKA's internal
!                      scheme, also valid for EVTXYZ).
!     Heavy ions (A>4): PDG-style extended encoding
!                      A*10 + Z*10000 + L*10000000 + 1e9. SGMXYZ
!                      decodes this internally via PDGION. EVTXYZ
!                      does NOT — callers generating events must
!                      convert via pdg_to_evt_code (which first calls
!                      PDGION and returns the -2 HEAVYION sentinel).
!----------------------------------------------------------------------!
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(PAPROP)'
      INCLUDE '(PART2)'

      INTEGER pdg_id, proj_code
      INTEGER ABSPDG, A, Z, L, KFLK
Cf2py intent(out) proj_code

      IF (ABS(pdg_id) .LT. 1000000000) THEN
         IF (pdg_id .EQ. 22) THEN
            proj_code = 7
         ELSE
            KFLK = MCIHAD(pdg_id)
            IF (KFLK .LT. -6 .OR. KFLK .GT. 390) THEN
               proj_code = 0
            ELSE
               proj_code = KPTOIP(KFLK)
            END IF
         END IF
      ELSE
         ABSPDG = ABS(pdg_id)
         A = MOD(ABSPDG / 10,        1000)
         Z = MOD(ABSPDG / 10000,     1000)
         L = MOD(ABSPDG / 10000000,  10)
         IF (L .EQ. 0 .AND. A .EQ. 2 .AND. Z .EQ. 1) THEN
            proj_code = -3          ! deuteron
         ELSE IF (L .EQ. 0 .AND. A .EQ. 3 .AND. Z .EQ. 1) THEN
            proj_code = -4          ! triton
         ELSE IF (L .EQ. 0 .AND. A .EQ. 3 .AND. Z .EQ. 2) THEN
            proj_code = -5          ! 3-He
         ELSE IF (L .EQ. 0 .AND. A .EQ. 4 .AND. Z .EQ. 2) THEN
            proj_code = -6          ! 4-He (alpha)
         ELSE
!           Heavy ion or hyper-nucleus: PDG-extended encoding (SGMXYZ
!           decodes via PDGION; EVTXYZ requires pre-conversion).
            proj_code = A*10 + Z*10000 + L*10000000 + 1000000000
            IF (pdg_id .LT. 0) proj_code = -proj_code
         END IF
      END IF

      RETURN
      END SUBROUTINE pdg_to_proj_code

      SUBROUTINE pdg_to_evt_code(pdg_id, evt_code)
!----------------------------------------------------------------------!
!     Projectile-code variant for EVTXYZ. Same as pdg_to_proj_code
!     except heavy ions (A > 4) are registered via PDGION into FLUKA's
!     ion common blocks and returned as the -2 HEAVYION sentinel.
!     EVTXYZ rejects the PDG-extended encoding (aborts with
!     KPPRCT=-1/-2), whereas SGMXYZ accepts it directly.
!----------------------------------------------------------------------!
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(PAPROP)'
      INCLUDE '(PART2)'

      INTEGER pdg_id, evt_code
      INTEGER ABSPDG, A, Z, L, KDUMMY
      INTEGER PDGION
Cf2py intent(out) evt_code

!     Start from the xsec encoding, then replace heavy-ion extended
!     codes with the -2 HEAVYION sentinel after PDGION registration.
      CALL pdg_to_proj_code(pdg_id, evt_code)
      IF (ABS(evt_code) .LT. 1000000000) RETURN

      ABSPDG = ABS(pdg_id)
      A = MOD(ABSPDG / 10,        1000)
      Z = MOD(ABSPDG / 10000,     1000)
      L = MOD(ABSPDG / 10000000,  10)
!     Only attempt PDGION for non-strange standard nuclei — the rest
!     falls through to the original extended encoding, which EVTXYZ
!     still rejects but which we at least don't hide.
      IF (L .EQ. 0) THEN
         KDUMMY = PDGION(ABSPDG)
         evt_code = -2
      END IF

      RETURN
      END SUBROUTINE pdg_to_evt_code

      SUBROUTINE fluka_elem_properties(n_materials, mat_idx,
     &                                 z_out, a_out, mass_out)
!----------------------------------------------------------------------!
!     Read FLUKA's FLKMAT common after STPXYZ to verify which
!     materials are registered. Used by Python to build the
!     pdg -> fluka material-index map and by tests.
!
!     Input:
!        n_materials     -- number of materials to read
!        mat_idx(n_mat)  -- FLUKA material indices (from MTFLKA of STPXYZ)
!
!     Output:
!        z_out(n_mat)    -- atomic number Z of each material
!        a_out(n_mat)    -- atomic weight A of each material
!        mass_out(n_mat) -- atomic mass AMSS of each material (g/mol)
!----------------------------------------------------------------------!
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(FLKMAT)'

      INTEGER n_materials
      INTEGER mat_idx(n_materials)
      INTEGER z_out(n_materials)
      INTEGER a_out(n_materials)
      DOUBLE PRECISION mass_out(n_materials)
      INTEGER i, mi
Cf2py intent(out) z_out, a_out, mass_out
Cf2py integer intent(hide),depend(mat_idx) :: n_materials=len(mat_idx)

      DO i = 1, n_materials
         mi = mat_idx(i)
         IF (mi .LT. 1 .OR. mi .GT. MXXMDF) THEN
            z_out(i) = -1
            a_out(i) = -1
            mass_out(i) = -1.0D+00
         ELSE
            z_out(i) = NINT(ZTAR(mi))
            a_out(i) = NINT(AMSS(mi))
            mass_out(i) = AMSS(mi)
         END IF
      END DO

      RETURN
      END SUBROUTINE fluka_elem_properties

      SUBROUTINE fluka_hepevt_summary(nhep_total, n_standard,
     &                                n_heavy, n_residual)
!----------------------------------------------------------------------!
!     Scan HEPEVT (populated by FLLHEP) and count entries by category.
!     Used by tests to assert that nuclear remnants are present.
!
!       n_standard -- standard particles (|pid| < 1e9)
!       n_heavy    -- light nuclei/fragments (d, t, 3He, 4He, ...)
!                     counted as entries with |pid| >= 1e9 AND A <= 4
!       n_residual -- residual nuclei (|pid| >= 1e9 AND A > 4)
!
!     (The A threshold of 4 separates FLUKA's FHEAVY light fragments
!      from RESNUC evaporation residues / projectile remnants.)
!----------------------------------------------------------------------!
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(HEPCMM)'

      INTEGER nhep_total, n_standard, n_heavy, n_residual
      INTEGER i, pid, absp, a
Cf2py intent(out) nhep_total, n_standard, n_heavy, n_residual

      nhep_total = NHEP
      n_standard = 0
      n_heavy    = 0
      n_residual = 0

      DO i = 1, NHEP
         pid = IDHEP(i)
         absp = ABS(pid)
         IF (absp .LT. 1000000000) THEN
            n_standard = n_standard + 1
         ELSE
            a = MOD(absp / 10, 1000)
            IF (a .LE. 4) THEN
               n_heavy = n_heavy + 1
            ELSE
               n_residual = n_residual + 1
            END IF
         END IF
      END DO

      RETURN
      END SUBROUTINE fluka_hepevt_summary

!======================================================================!
!  Radioactive-decay tables: lightweight init (no STPXYZ).             !
!  Mirrors dcytst.f's main-program init sequence.  Idempotent: safe    !
!  to call after STPXYZ (which already loads the same tables).         !
!======================================================================!
      SUBROUTINE chromo_dcy_init()
      INCLUDE '(DBLPRW)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'

      EXTERNAL BDNOPT, BDEVAP, BDPREE, BDINPT, BDTRNS, BDPART, BDPRDC

      INCLUDE '(EMFCMP)'
      INCLUDE '(EMFMAT)'
      INCLUDE '(EMFSWT)'
      INCLUDE '(EMFTHR)'
      INCLUDE '(EVAPRD)'
      INCLUDE '(FLUOXR)'
      INCLUDE '(FLKMAT)'
      INCLUDE '(FRBKCM)'
      INCLUDE '(NUCDAT)'
      INCLUDE '(PAPROP)'
      INCLUDE '(PAREVT)'
      INCLUDE '(TRACKR)'

      INTEGER MMAT, MREG
      LOGICAL LDCYINI
      SAVE    LDCYINI
      DATA    LDCYINI /.FALSE./

      IF (LDCYINI) RETURN

      CALL CMSPPR
      CALL ZEROIN
      LEVPRT = .TRUE.
      LDEEXG = .TRUE.
      LHEAVY = .TRUE.
      LGDHPR = .TRUE.
      LFRMBK = .FALSE.
      CALL NCDTRD
      CALL KPIXSR
      CALL INCINI
      AMUMEV = GEVMEV * AMUAMU

      LFLUKA = .FALSE.
      LEMFON = .TRUE.
      MMAT   = 3
      NMAT   = MMAT
      MREG   = 1
      IPRODC = 2
      MMTRCK = MMAT
      MRTRCK = MREG
      NMDEMF = 1
      METOFL(NMDEMF) = MMAT
      MFLTOE(MMAT)   = NMDEMF
      LXFLUO(NMDEMF) = .TRUE.
      MEDFLK(MREG,1) = MMAT
      MEDFLK(MREG,2) = MMAT
      EPHMIN(NMDEMF) = 1.D-07 * GV2EMF
      EEPMIN(NMDEMF) = (AMELCT + 1.D-07) * GV2EMF
      CALL RDFLUO

      LDCYINI = .TRUE.
      RETURN
      END SUBROUTINE chromo_dcy_init

!======================================================================!
!  ISMRCH wrapper: returns IFOUND=1 if (A,Z,m) has decay data,         !
!  IFOUND=0 otherwise.  Caller is responsible for staying within the   !
!  safe (A,Z) probe band (see chromo_dcy_catalog).                     !
!======================================================================!
      SUBROUTINE chromo_dcy_lookup(ia, iz, im, ifound,
     &                             t12, exm, jsp, jpt)
      INCLUDE '(DBLPRW)'
      INCLUDE '(DIMPAR)'

      INTEGER ia, iz, im, ifound, jsp, jpt
      DOUBLE PRECISION t12, exm
Cf2py intent(in)  ia, iz, im
Cf2py intent(out) ifound, t12, exm, jsp, jpt

      INTEGER kisitp

      kisitp = 0
      ifound = 0
      t12 = 0.0D0
      exm = 0.0D0
      jsp = 0
      jpt = 0

!  Guard against ISMRCH OOB on pathological inputs.  iz>ia is a
!  physical-nuclide constraint (Z<=A always) that also defangs FLUKA's
!  EXMSAZ:KA0-KZ0 abort path; do not remove without retesting (2,50,0).
      IF (ia .LT. 1) RETURN
      IF (iz .LT. 0) RETURN
      IF (ia .EQ. 1 .AND. iz .GT. 1) RETURN
      IF (ia .GT. 1 .AND. iz .LT. 1) RETURN
      IF (iz .GT. ia) RETURN
      IF (im .LT. 0 .OR. im .GT. 4) RETURN

      CALL ISMRCH(ia, iz, im, kisitp, t12, exm, jsp, jpt)

      IF (kisitp .GT. 0) ifound = 1

      RETURN
      END SUBROUTINE chromo_dcy_lookup

!======================================================================!
!  Bulk walk of the FLUKA decay-data catalogue.  Conservative valley-  !
!  of-stability Z-band per A; ground-state probe gates isomer probes.  !
!  Returns N total entries written to caller-allocated arrays.         !
!                                                                      !
!  Arrays are declared (*) (assumed-size) so f2py exposes max_n as a   !
!  regular positional input rather than auto-hiding it from the array  !
!  shape; intent(inplace) keeps writes in the caller's buffers and     !
!  avoids returning resized copies.                                    !
!======================================================================!
      SUBROUTINE chromo_dcy_catalog(max_n, n_out,
     &                              a_out, z_out, m_out,
     &                              t12_out, exm_out,
     &                              jsp_out, jpt_out)
      INCLUDE '(DBLPRW)'
      INCLUDE '(DIMPAR)'

      INTEGER max_n, n_out
      INTEGER a_out(*), z_out(*), m_out(*)
      INTEGER jsp_out(*), jpt_out(*)
      DOUBLE PRECISION t12_out(*), exm_out(*)
Cf2py intent(in) max_n
Cf2py intent(inplace) a_out, z_out, m_out, t12_out, exm_out, jsp_out, jpt_out
Cf2py intent(out) n_out

      INTEGER ia, iz, im, izmin, izmax, kisitp, jsp, jpt
      DOUBLE PRECISION t12, exm

      n_out = 0
!  Empirical ceiling; (ISOTOP)'s NAMSMX=330 but FLUKA's catalogue is
!  empty above A~290.  Probing higher costs ~2700 dead cells with no
!  benefit on the one-shot init path.
      DO ia = 1, 295
         IF (ia .LE. 4) THEN
            izmin = 1
            izmax = MIN(ia, 110)
         ELSE IF (ia .LE. 20) THEN
            izmin = MAX(1, INT(0.30D0*ia) - 4)
            izmax = MIN(ia, 110)
         ELSE
            izmin = MAX(1, INT(0.30D0*ia) - 5)
            izmax = MIN(110, INT(0.55D0*ia) + 8)
         END IF
         DO iz = izmin, izmax
!  Probe ground state first; only try isomers if g.s. exists.
            kisitp = 0
            CALL ISMRCH(ia, iz, 0, kisitp, t12, exm, jsp, jpt)
            IF (kisitp .GT. 0) THEN
               IF (n_out .GE. max_n) RETURN
               n_out = n_out + 1
               a_out(n_out)   = ia
               z_out(n_out)   = iz
               m_out(n_out)   = 0
               t12_out(n_out) = t12
               exm_out(n_out) = exm
               jsp_out(n_out) = jsp
               jpt_out(n_out) = jpt
               IF (ia .GE. 5) THEN
                  DO im = 1, 4
                     kisitp = 0
                     CALL ISMRCH(ia, iz, im, kisitp,
     &                           t12, exm, jsp, jpt)
                     IF (kisitp .GT. 0) THEN
                        IF (n_out .GE. max_n) RETURN
                        n_out = n_out + 1
                        a_out(n_out)   = ia
                        z_out(n_out)   = iz
                        m_out(n_out)   = im
                        t12_out(n_out) = t12
                        exm_out(n_out) = exm
                        jsp_out(n_out) = jsp
                        jpt_out(n_out) = jpt
                     END IF
                  END DO
               END IF
            END IF
         END DO
      END DO

      RETURN
      END SUBROUTINE chromo_dcy_catalog

!======================================================================!
!  Per-isotope decay-channel data: branching ratios, daughter (A,Z,m), !
!  Q value (via QRDDCY).  KIND code: 1=alpha, 2=B-, 3=B+, 4=EC, 5=IT,  !
!  6=SF, 7=B-N, 8=B+P, 9=B-2N, 10=B-3N, 11=B-NA, 12=other.             !
!======================================================================!
      SUBROUTINE chromo_dcy_channels(ia, iz, im, max_ch, n_ch,
     &                               kind_ch, br_ch,
     &                               da_ch, dz_ch, dm_ch, q_ch)
      INCLUDE '(DBLPRW)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IDPPRM)'
      INCLUDE '(ISOTOP)'

      INTEGER ia, iz, im, max_ch, n_ch
      INTEGER kind_ch(*), da_ch(*), dz_ch(*), dm_ch(*)
      DOUBLE PRECISION br_ch(*), q_ch(*)
Cf2py intent(in) ia, iz, im, max_ch
Cf2py intent(inplace) kind_ch, br_ch, da_ch, dz_ch, dm_ch, q_ch
Cf2py intent(out) n_ch

      INTEGER kisitp, jsporu, jptoru, iia, jdcy, kdcy
      DOUBLE PRECISION t12, exm, brdctt
      DOUBLE PRECISION QRDDCY
      EXTERNAL QRDDCY

      n_ch = 0
      kisitp = 0
      CALL ISMRCH(ia, iz, im, kisitp, t12, exm, jsporu, jptoru)
      IF (kisitp .LE. 0) RETURN

      iia = KQMDIN(ia, KMSFVR)
      brdctt = 0.0D0

      DO jdcy = 1, MXDCIS
!  Look up the channel index first; skip empty slots without writing.
         IF (im .LE. 0) THEN
            kdcy = IDCNUC(jdcy, iia, kisitp)
         ELSE
            kdcy = IDCISM(jdcy, kisitp)
         END IF
         IF (kdcy .LE. 0) CYCLE
         IF (n_ch .GE. max_ch) RETURN
         n_ch = n_ch + 1

!  BR write happens AFTER the cap check so an under-sized buffer cannot
!  leave a stale value in the last accepted slot (was a Task-4 bug).
         IF (im .LE. 0) THEN
            IF (jdcy .LT. MXDCIS) THEN
               br_ch(n_ch) = BRDECY(jdcy, iia, kisitp)
            ELSE
               br_ch(n_ch) = 1.0D0 - brdctt
            END IF
         ELSE
            IF (jdcy .LT. MXDCIS) THEN
               br_ch(n_ch) = BRDISM(jdcy, kisitp)
            ELSE
               br_ch(n_ch) = 1.0D0 - brdctt
            END IF
         END IF
         brdctt = brdctt + br_ch(n_ch)

!  Map FLUKA kdcy to a compact KIND code.  See (ISOTOP) for KDCY*
!  parameters: KDCYAL=2 (alpha), KDCYBM=22 (B-), KDCYBP=8 (B+),
!  KDCYEC=15 (EC), KDCYIT=1 (IT), KDCYSF=38 (SF), KDCBMN=23 (B-N),
!  KDCBPP=24 (B+P), KDCB2N=28 (B-2N), KDCB3N=29 (B-3N), KDCBNA=30 (B-NA).
         IF (kdcy .EQ. KDCYAL .OR. kdcy .EQ. KDCALM) THEN
            kind_ch(n_ch) = 1
         ELSE IF (kdcy .EQ. KDCYBM .OR. kdcy .EQ. KDCBMM) THEN
            kind_ch(n_ch) = 2
         ELSE IF (kdcy .EQ. KDCYBP .OR. kdcy .EQ. KDCBPM) THEN
            kind_ch(n_ch) = 3
         ELSE IF (kdcy .EQ. KDCYEC .OR. kdcy .EQ. KDCECM) THEN
            kind_ch(n_ch) = 4
         ELSE IF (kdcy .EQ. KDCYIT .OR. kdcy .EQ. KDCITM) THEN
            kind_ch(n_ch) = 5
         ELSE IF (kdcy .EQ. KDCYSF) THEN
            kind_ch(n_ch) = 6
         ELSE IF (kdcy .EQ. KDCBMN) THEN
            kind_ch(n_ch) = 7
         ELSE IF (kdcy .EQ. KDCBPP) THEN
            kind_ch(n_ch) = 8
         ELSE IF (kdcy .EQ. KDCB2N) THEN
            kind_ch(n_ch) = 9
         ELSE IF (kdcy .EQ. KDCB3N) THEN
            kind_ch(n_ch) = 10
         ELSE IF (kdcy .EQ. KDCBNA) THEN
            kind_ch(n_ch) = 11
         ELSE
            kind_ch(n_ch) = 12
         END IF

         IF (IDCYDA(kdcy) .GT. -100) THEN
            da_ch(n_ch) = ia + IDCYDA(kdcy)
            dz_ch(n_ch) = iz + IDCYDZ(kdcy)
            IF (kdcy .GT. NDCY1M) THEN
               dm_ch(n_ch) = 2
            ELSE IF (kdcy .GT. NDCYGS) THEN
               dm_ch(n_ch) = 1
            ELSE
               dm_ch(n_ch) = 0
            END IF
            q_ch(n_ch) = QRDDCY(ia, iz, im, kdcy, .FALSE.) * GEVMEV
         ELSE
            da_ch(n_ch) = -1
            dz_ch(n_ch) = -1
            dm_ch(n_ch) = -1
            q_ch(n_ch) = 0.0D0
         END IF
      END DO

      RETURN
      END SUBROUTINE chromo_dcy_channels

*$ CREATE QRDDCY.FOR
*COPY QRDDCY
*
*=== Qrddcy ===========================================================*
*
*     QRDDCY is shipped as user-supplied source in FLUKA 2025.x's
*     decay-test harness (dcytst.f) and is NOT exported by libflukahp.a.
*     We replicate it verbatim here so chromo_dcy_channels can resolve
*     the Q-value at link time.
*
*     Source vintage: FLUKA 2025.1 dcytst.f (Last change on 25-Apr-26
*     by Alfredo Ferrari).
*
*     TODO(FLUKA-bump): re-sync this body from $FLUPRO/dcytst.f whenever
*     the FLUKA version is bumped.  See FLUKA_QUESTIONS.md.
*
      DOUBLE PRECISION FUNCTION QRDDCY ( IADCYP, IZDCYP, ISDCYP, IFLDCY,
     &                                   LNCMSS )

      INCLUDE '(DBLPRN)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2024-2026       by   Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*     Q's for RaDioactive DeCaY                                        *
*                                                                      *
*----------------------------------------------------------------------*
*
      LOGICAL          LNCMSS
      INTEGER          IADCYP, IZDCYP, ISDCYP, IFLDCY
*
      INCLUDE '(IDPPRM)'
      INCLUDE '(ISOTOP)'
      INCLUDE '(PAPROP)'

      LOGICAL          LECBPL, LBETAP
      INTEGER          NASUM , IZDUM , JZDUM , JA    , JZ    , JS    ,
     &                 KZDUM , KA    , KZ    , KS    , K
      DOUBLE PRECISION AMDCAY, AMDGHT, AMPRTC
      DOUBLE PRECISION EXMAZM
      EXTERNAL         EXMAZM

      QRDDCY = ZERZER

      IF ( IDCYDA (IFLDCY) .GT. -100 ) THEN
         NASUM  = IADCYP
         AMDCAY = EMVGEV * EXMAZM ( IADCYP, IZDCYP, ISDCYP,
     &                              LNCMSS, IZDUM )
         KA     = IADCYP + IDCYDA (IFLDCY)
         KZ     = IZDCYP + IDCYDZ (IFLDCY)
         IF ( IFLDCY .GT. NDCY1M ) THEN
            KS  = 2
         ELSE IF ( IFLDCY .GT. NDCYGS ) THEN
            KS  = 1
         ELSE
            KS  = 0
         END IF
         NASUM  = NASUM - KA
         AMDGHT = EMVGEV * EXMAZM ( KA, KZ, KS, LNCMSS, KZDUM )
         QRDDCY = AMDCAY - AMDGHT
         LECBPL = .FALSE.
         LBETAP = .FALSE.
         DO 1000 K = 1, NDCYPR (IFLDCY)
            IF ( KDCYPR (K,IFLDCY) .EQ. IJNUEL .OR.
     &           KDCYPR (K,IFLDCY) .EQ. IJANUE ) THEN
               LECBPL = KDCYPR (K,IFLDCY) .EQ. IJNUEL
               GO TO 1000
            ELSE IF ( KDCYPR (K,IFLDCY) .EQ. IJELCT ) THEN
               IF ( LNCMSS ) THEN
                  QRDDCY = QRDDCY - AMELCT
               ELSE
                  GO TO 1000
               END IF
            ELSE IF ( KDCYPR (K,IFLDCY) .EQ. IJPOST ) THEN
               LBETAP = .TRUE.
               IF ( LNCMSS ) THEN
                  QRDDCY = QRDDCY - AMELCT
               ELSE
                  QRDDCY = QRDDCY - TWOTWO * AMELCT
               END IF
            ELSE IF ( KDCYPR (K,IFLDCY) .GE. IJALPH .AND.
     &                KDCYPR (K,IFLDCY) .LE. NALLWP ) THEN
               JA     = IBARCH ( KDCYPR (K,IFLDCY) )
               JZ     = ICHRGE ( KDCYPR (K,IFLDCY) )
               JS     = 0
               NASUM  = NASUM  - JA
               AMPRTC = EMVGEV * EXMAZM ( JA, JZ, JS, LNCMSS, JZDUM )
               QRDDCY = QRDDCY - AMPRTC
            ELSE IF ( KDCYPR (K,IFLDCY) .LT. IJALPH ) THEN
               JA     = MOD ( ABS ( KDCYPR (K,IFLDCY) ), 100000 )
     &                / 100
               JZ     = MOD ( ABS ( KDCYPR (K,IFLDCY) ), 10000000 )
     &                / 100000
               JS     = MOD ( ABS ( KDCYPR (K,IFLDCY) ), 1000000000 )
     &                / 100000000
               NASUM  = NASUM - JA
               AMPRTC = EMVGEV * EXMAZM ( JA, JZ, JS, LNCMSS, JZDUM )
               QRDDCY = QRDDCY - AMPRTC
            ELSE
               CALL FLABRT ( 'QRDDCY', 'INVALID KDCYPR' )
               STOP
            END IF
 1000    CONTINUE
         IF ( NASUM  .NE. 0 ) THEN
            CALL FLABRT ( 'QRDDCY', 'NASUM!=0' )
            STOP
         ELSE IF ( LBETAP .AND. .NOT. LECBPL ) THEN
            CALL FLABRT ( 'QRDDCY', 'LBETAP&&!LECBPL' )
            STOP
         END IF
         IF ( LECBPL .AND. .NOT. LBETAP .AND. LNCMSS ) THEN
            QRDDCY = QRDDCY + AMELCT
            QRDDCY = QRDDCY - AINFNT
         END IF
      END IF

      RETURN
      END

!======================================================================!
!  Per-isotope line lists.  KIND: 1=gamma (BR, E), 2=alpha (BR, E,     !
!  end-level), 3=CE/Auger (BR, E), 4=beta+/- (BR, <E>, end-point E,    !
!  end-level; positron flag encoded in sign of <E>).                   !
!  Energies returned in MeV.                                           !
!  SIGGTT is the REAL*4 view of the FLUKA blank common; it is brought  !
!  in via (DBLPRW) -> (USFLMD) -> USE AABLMD.  Do NOT redeclare it.    !
!======================================================================!
      SUBROUTINE chromo_dcy_lines(ia, iz, im, kind, max_l, n_l,
     &                            br_l, e_l, nlev_l)
      INCLUDE '(DBLPRW)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IDPPRM)'
      INCLUDE '(ISOTOP)'

      INTEGER ia, iz, im, kind, max_l, n_l
      INTEGER nlev_l(*)
      DOUBLE PRECISION br_l(*), e_l(*)
Cf2py intent(in) ia, iz, im, kind, max_l
Cf2py intent(inplace) br_l, e_l, nlev_l
Cf2py intent(out) n_l

      INTEGER kisitp, jsporu, jptoru, iia, kn_lines, kp_lines, ipt
      INTEGER ndata
      DOUBLE PRECISION t12, exm

      n_l = 0
      kisitp = 0
      CALL ISMRCH(ia, iz, im, kisitp, t12, exm, jsporu, jptoru)
      IF (kisitp .LE. 0) RETURN

      iia = KQMDIN(ia, KMSFVR)
      kn_lines = 0
      kp_lines = 0

      IF (im .LE. 0) THEN
         IF (kind .EQ. 1) THEN
            kn_lines = NGMLNS(iia, kisitp)
            kp_lines = KGMLNS(iia, kisitp)
            ndata = 2
         ELSE IF (kind .EQ. 2) THEN
            kn_lines = NALLNS(iia, kisitp)
            kp_lines = KALLNS(iia, kisitp)
            ndata = 3
         ELSE IF (kind .EQ. 3) THEN
            kn_lines = NCELNS(iia, kisitp)
            kp_lines = KCELNS(iia, kisitp)
            ndata = 2
         ELSE IF (kind .EQ. 4) THEN
            kn_lines = NBTSPC(iia, kisitp)
            kp_lines = KBTSPC(iia, kisitp)
            ndata = 4
         ELSE
            RETURN
         END IF
      ELSE
         IF (kind .EQ. 1) THEN
            kn_lines = NGMISM(kisitp)
            kp_lines = KGMISM(kisitp)
            ndata = 2
         ELSE IF (kind .EQ. 2) THEN
            kn_lines = NALISM(kisitp)
            kp_lines = KALISM(kisitp)
            ndata = 3
         ELSE IF (kind .EQ. 3) THEN
            kn_lines = NCEISM(kisitp)
            kp_lines = KCEISM(kisitp)
            ndata = 2
         ELSE IF (kind .EQ. 4) THEN
            kn_lines = NBTISM(kisitp)
            kp_lines = KBTISM(kisitp)
            ndata = 4
         ELSE
            RETURN
         END IF
      END IF

      DO ipt = 1, kn_lines
         IF (n_l .GE. max_l) RETURN
         n_l = n_l + 1
         br_l(n_l) = ABS(DBLE(SIGGTT(kp_lines + ndata*(ipt-1) + 1)))
         e_l(n_l)  = DBLE(SIGGTT(kp_lines + ndata*(ipt-1) + 2)) * GEVMEV
         IF (ndata .EQ. 3) THEN
            nlev_l(n_l) = NINT(SIGGTT(kp_lines + ndata*(ipt-1) + 3))
         ELSE IF (ndata .EQ. 4) THEN
!  Beta: index +2 holds <E> (signed), +3 end-point E, +4 end-level
            e_l(n_l) = DBLE(SIGGTT(kp_lines + ndata*(ipt-1) + 3))
     &                * GEVMEV
            nlev_l(n_l) = NINT(SIGGTT(kp_lines + ndata*(ipt-1) + 4))
         ELSE
            nlev_l(n_l) = 0
         END IF
      END DO

      RETURN
      END SUBROUTINE chromo_dcy_lines
