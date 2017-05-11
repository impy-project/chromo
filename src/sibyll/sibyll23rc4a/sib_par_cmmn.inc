C--------------------------------------------------------------------
C     SIBYLL common blocks containing parameters       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)

      CHARACTER*6 NAMP
      COMMON /S_CNAM/ NAMP (0:99)

      DOUBLE PRECISION STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea

      DOUBLE PRECISION XM2MIN,ALXMIN,SLOP0,ASLOP,BSLOP,XMASS
      COMMON /S_DIFMAss/ XM2MIN(3),ALXMIN(3),SLOP0,ASLOP,BSLOP,XMASS(2)

      DOUBLE PRECISION PPT0,ptflag,PPT02
      COMMON /S_CQDIS/ PPT0(35),ptflag
      COMMON /S_CQDIS2/ PPT02(44)
