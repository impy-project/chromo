*$ CREATE DPMJET.FOR
*COPY DPMJET
*
*===program dpmjet=====================================================*
*
      PROGRAM DPMJET

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

* block data in DPMJET library (uncomment these declarations if library
* option is used)
C     EXTERNAL DT_BDEVAP,DT_BDNOPT,DT_BDPREE,DT_HADPRP,DT_BLKD46,
C    &         DT_BLKD47,DT_RUNTT,DT_NONAME,DT_ZK,DT_BLKD43

C     EXTERNAL PYDATA

* event flag
      COMMON /DTEVNO/ NEVENT,ICASCA

*-----------------------------------------------------------------------
* initialization

*   the following statement provides a call to DT_USRHIS(MODE=1) for
*   histogram initialization etc.
      CALL DT_DTUINI(NEVTS,EPN,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IEMU)
      WRITE (6,*) "****************************"
      WRITE (6,*) NEVTS,EPN,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IEMU
*-----------------------------------------------------------------------
* generation of events

      DO 1 IEVT=1,NEVTS

*   some defaults, do not change!
         NEVENT = IEVT
         KKMAT  = -1
         ELAB   = EPN
*   uncomment if dpmjet3 is linked to particle transport code
C        ICASCA = 1

************************************************************************
* The following lines show how to select the target nucleus for runs
* with composite targets (and fixed projectile and energy!).
*
*   Sampling of the target nucleus (mass number NTMASS, charge NTCHAR)
*   according to the fractions defined with EMULSION input-cards.
*   The different nuclei are numbered as KKMAT = 1,2,3,...  according to
*   their appearance in the input-file.
         IF (IEMU.GT.0) THEN
*   Replace this selection by your own one if needed.
            CALL DT_GETEMU(NTMASS,NTCHAR,KKMAT,0)
*   Kkmat has to be negative for composite targets!
            KKMAT = -KKMAT
         ENDIF
************************************************************************

************************************************************************
* The following lines show how to define projectile, target and energy
* for this event in runs with Glauber-data file pre-initialized for a
* certain range of projectiles, targets and energies. The definitions
* have to be within the pre-initialized parameter range.
*
*   projectile-id (for hadron projectiles)
C        IDP    = 1
*   projectile mass and charge numbers
C        NPMASS = 12
C        NPCHAR = 6
*   target mass and charge numbers
C        NTMASS = 16
C        NTCHAR = 8
*   lab energy
C        ELAB = 200.0D0
************************************************************************

************************************************************************
* If an energy-range has been defined with the ENERGY input-card the
* laboratory energy ELAB can be set to any value within that range. For
* example:
C        ELO  = 10.0D0
C        EHI  = 1000.0D0
C        ELAB = DT_RNDM(ELAB)*(EHI-ELO)+ELO
************************************************************************

*   sampling of one event
         CALL DT_KKINC(NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,ELAB,KKMAT,IREJ)
         IF (IREJ.NE.0) GOTO 1

*   the following statement provides a call to DT_USRHIS(MODE=2) from
*   where the final state particles can be obtained

         CALL PHO_PHIST(2000,DUM)

    1 CONTINUE

*-----------------------------------------------------------------------
* output, statistics etc.

*   the following statement provides a call to DT_USRHIS(MODE=3) in
*   order to calculate histograms etc.
      CALL DT_DTUOUT

      END

*$ CREATE DT_USRHIS.FOR
*COPY DT_USRHIS
*
*===usrhis=============================================================*
*
      SUBROUTINE DT_USRHIS(MODE)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*
* COMMON /DTEVT1/ :
*                   NHKK         number of entries in common block
*                   NEVHKK       number of the event
*                   ISTHKK(i)    status code for entry i
*                   IDHKK(i)     identifier for the entry
*                                (for particles: identifier according
*                                 to the PDG numbering scheme)
*                   JMOHKK(1,i)  pointer to the entry of the first mother
*                                of entry i
*                   JMOHKK(2,i)  pointer to the entry of the second mother
*                                of entry i
*                   JDAHKK(1,i)  pointer to the entry of the first daughter
*                                of entry i
*                   JDAHKK(2,i)  pointer to the entry of the second daughter
*                                of entry i
*                   PHKK(1..3,i) 3-momentum
*                   PHKK(4,i)    energy
*                   PHKK(5,i)    mass
*
* event history
      PARAMETER (NMXHKK=200000)
      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
* extended event history
      COMMON /DTEVT2/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10),
     &                IHIST(2,NMXHKK)

      GOTO (1,2,3) MODE

*------------------------------------------------------------------
*
    1 CONTINUE
*
* initializations
*
*  Called with MODE=1 once at the beginning of the run.
*
      RETURN
*
*------------------------------------------------------------------
*
    2 CONTINUE
*
* scoring of the present event
*
*  Called with MODE=2 every time one event has been finished.
*
*  The final state particles from the actual event (number NEVHKK)
*  can be found in DTEVT1 and identified by their status:
*
*     ISTHKK(i) = 1    final state particle produced in
*                      photon-/hadron-/nucleon-nucleon collisions or
*                      in intranuclear cascade processes
*                -1    nucleons, deuterons, H-3, He-3, He-4 evaporated
*                      from excited nucleus and
*                      photons produced in nuclear deexcitation processes
*                1001  residual nucleus (ground state)
*
*  The types of these particles/nuclei are given in IDHKK as follows
*
*     all final state part. except nuclei :
*       IDHKK(i)=particle identifier according to PDG numbering scheme
*     nuclei (evaporation products, and residual nucleus) :
*       IDHKK(i)=80000, IDRES(i)=mass number, IDXRES(i)=charge number
*
*  The 4-momenta and masses can be found in PHKK (target nucleus rest frame):
*                   PHKK(1..3,i) 3-momentum (p_x,p_y,p_z)
*                   PHKK(4,i)    energy
*                   PHKK(5,i)    mass
*
*
*
*  Pick out the final state particles from DTEVT1 in each event for
*  instance by the following loop (NHKK=number of entries in the present
*  event) and fill your histograms
C     DO 20 I=1,NHKK
C        IF (ABS(ISTHKK(I)).EQ.1) THEN
C        ELSEIF (ABS(ISTHKK(I)).EQ.1001) THEN
C        ENDIF
C  20 CONTINUE

*  At any time during the run a list of the actual entries in DTEVT1 and
*  DTEVT2 can be obtained (output unit 6) by the following statement:
C     CALL DT_EVTOUT(4)

      RETURN
*
*------------------------------------------------------------------
*
    3 CONTINUE
*
* output/statistics/histograms etc.
*
*  Called with MODE=3 once after all events have been sampled.
*
      RETURN

      END
