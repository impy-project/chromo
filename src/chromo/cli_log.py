"""
Logging functionality for the command-line interface.
"""


def get_final_state_particles_info(model):
    """Get information about final state particles."""
    from particle import Particle

    from chromo.constants import long_lived

    particles = model.final_state_particles
    if not particles:
        return ""

    output = "\n\n[Final State Particles]"

    # Check if using default long_lived list
    is_default = set(particles) == set(long_lived)
    n_decayable = len(model._unstable_pids)

    if is_default:
        output += "\nUsing default configuration (particles with lifetime > 30 ps):"
    else:
        output += "\nUsing custom configuration:"

    output += f"\n  Generator can decay: {n_decayable} particle types"
    output += f"\n  Configured as stable: {len(particles)} particle types\n"

    # Sort by absolute value of PDG ID, with negative first, then positive
    sorted_particles = sorted(particles, key=lambda x: (abs(x), x))

    # Table header
    output += f"\n{'Particle':<20} {'PDG ID':>10}  {'Lifetime (s)':>15}"
    output += f"\n{'-'*20} {'-'*10}  {'-'*15}"

    c_mm_per_s = 299792458000.0  # Speed of light in mm/s

    for pdgid in sorted_particles:
        p = Particle.findall(pdgid=pdgid)
        if not p:
            output += f"\n{'Unknown':<20} {pdgid:>10}  {'N/A':>15}"
            continue

        p = p[0]
        name = p.name

        # ctau is in mm, convert to lifetime in seconds
        # lifetime = ctau / c
        if p.ctau is not None and p.ctau != float("inf"):
            lifetime = p.ctau / c_mm_per_s
        else:
            lifetime = float("inf")

        if lifetime == float("inf"):
            lifetime_str = "stable"
        else:
            lifetime_str = f"{lifetime:.6e}"

        output += f"\n{name:<20} {pdgid:>10}  {lifetime_str:>15}"

    return output


def get_pythia8_settings(model, log_level="normal"):
    """Get Pythia8-specific settings."""
    pythia = model._pythia
    output = "\n\n[Pythia8 Settings]"

    output += "\n\nUser configuration:"
    for s in model._final_config:
        output += f"\n  {s}"

    if log_level in ["normal", "verbose"]:
        import io
        from contextlib import redirect_stdout

        buffer = io.StringIO()
        with redirect_stdout(buffer):
            pythia.settings.listChanged()
        settings_output = buffer.getvalue()
        if settings_output.strip():
            output += "\n\nAll changed settings:"
            output += f"\n{settings_output}"

    if log_level == "verbose":
        import io
        from contextlib import redirect_stdout

        buffer = io.StringIO()
        with redirect_stdout(buffer):
            pythia.settings.listAll()
        settings_output = buffer.getvalue()
        if settings_output.strip():
            output += "\n\nAll settings (including defaults):"
            output += f"\n{settings_output}"

    return output


def get_pythia6_settings(model, log_level="normal"):
    """Get Pythia6-specific settings."""
    output = "\n\n[Pythia6 Settings]"

    if log_level in ["normal", "verbose"]:
        # Process selection (MSEL=0 means manual selection via MSUB)
        msel = model._lib.pysubs.msel
        output += f"\nProcess selection (MSEL): {msel}"
        if msel == 0:
            output += " (custom - using MSUB flags)"

        # Enabled subprocesses - only show if custom selection
        if msel == 0:
            subprocess_names = {
                11: "f + f' → f + f' (QCD)",
                12: "f + fbar → f' + fbar' (QCD)",
                13: "f + fbar → g + g (QCD)",
                28: "f + g → f + g (QCD)",
                53: "g + g → f + fbar (QCD)",
                68: "g + g → g + g (QCD)",
                92: "Single diffractive (XB)",
                93: "Single diffractive (AX)",
                94: "Double diffractive",
                95: "Low-pT scattering",
                96: "Semihard QCD 2→2",
            }
            enabled_subs = [i + 1 for i in range(500) if model._lib.pysubs.msub[i] == 1]

            if enabled_subs:
                output += "\n\nEnabled subprocesses:"
                for isub in enabled_subs:
                    name = subprocess_names.get(isub, f"subprocess {isub}")
                    output += f"\n  MSUB({isub}) = 1  ({name})"

        # Multiple interaction model
        mstp82 = model._lib.pypars.mstp[81]  # MSTP(82): MPI model
        mpi_models = {
            1: "old model (default)",
            2: "intermediate model",
            3: "new model",
            4: "newest model (recommended)",
        }
        output += f"\n\nMultiple interactions: MSTP(82) = {mstp82}"
        if mstp82 in mpi_models:
            output += f" ({mpi_models[mstp82]})"

        # PDF set
        mstp51 = model._lib.pypars.mstp[50]  # MSTP(51): PDF set
        output += f"\nPDF set: MSTP(51) = {mstp51}"

        # Fragmentation model
        mstj11 = model._lib.pydat1.mstj[10]  # MSTJ(11): fragmentation function
        frag_models = {
            1: "Lund symmetric",
            2: "Lund asymmetric",
            3: "independent fragmentation",
            4: "Field-Feynman",
            5: "Lund symmetric + independent",
        }
        output += f"\nFragmentation: MSTJ(11) = {mstj11}"
        if mstj11 in frag_models:
            output += f" ({frag_models[mstj11]})"

    if log_level == "verbose":
        # Hadronization parameters
        parj1 = model._lib.pydat1.parj[0]  # PARJ(1): diquark suppression
        parj2 = model._lib.pydat1.parj[1]  # PARJ(2): s-quark suppression
        parj3 = model._lib.pydat1.parj[2]  # PARJ(3): strange diquark suppression
        parj21 = model._lib.pydat1.parj[20]  # PARJ(21): fragmentation width

        output += "\n\nFragmentation parameters:"
        output += f"\n  PARJ(1) = {parj1:.3f}  (diquark suppression)"
        output += f"\n  PARJ(2) = {parj2:.3f}  (s-quark suppression)"
        output += f"\n  PARJ(3) = {parj3:.3f}  (strange diquark suppression)"
        output += f"\n  PARJ(21) = {parj21:.3f}  (width of fragmentation pT)"

        # MPI parameters
        parp82 = model._lib.pypars.parp[81]  # PARP(82): MPI cutoff
        parp89 = model._lib.pypars.parp[88]  # PARP(89): reference energy
        parp90 = model._lib.pypars.parp[89]  # PARP(90): energy scaling

        output += "\n\nMultiple interaction parameters:"
        output += f"\n  PARP(82) = {parp82:.3f}  (cutoff scale)"
        output += f"\n  PARP(89) = {parp89:.0f}  (reference energy)"
        output += f"\n  PARP(90) = {parp90:.3f}  (energy dependence)"

        output += "\n\nNote: Settings can be modified via model._lib common blocks:"
        output += "\n  - model._lib.pypars (MSTP/PARP: process switches & parameters)"
        output += "\n  - model._lib.pysubs (MSEL/MSUB: subprocess selection)"
        output += (
            "\n  - model._lib.pydat1 (MSTJ/PARJ: hadronization switches & parameters)"
        )

    return output


def get_sibyll_settings(model, log_level="normal"):
    """Get Sibyll-specific settings."""
    output = f"\n\n[{model.label} Settings]"

    if log_level in ["normal", "verbose"]:
        # Debug level from S_DEBUG common block (user-settable via chromo.debug_level)
        output += f"\nDebug level: {model._lib.s_debug.ndebug}"
        output += "\n  (Set via chromo.debug_level before model initialization)"

        # S* parameter for StarMixed models from S_STAR common block
        model_type = type(model)
        from chromo import models

        if model_type in (
            models.Sibyll23dStarMixed,
            models.Sibyll23eStarMixed,
        ):
            imod_labels = {
                0: "no enhancement",
                1: "rho enhancement",
                2: "baryon enhancement",
                3: "strange enhancement",
                4: "mixed enhancement",
            }
            imod = int(model._lib.s_star.imod)
            imod_desc = imod_labels.get(imod, "unknown")
            output += f"\nS* production mode: {imod} ({imod_desc})"
            output += "\n  (Set via model._lib.s_star.imod before event generation)"

        # Critical model switches (IPAR array) - exposed via s_cflafr
        try:
            # IPAR(15) - Charm production
            charm_enabled = int(model._lib.s_cflafr.ipar[14])
            charm_status = "ENABLED" if charm_enabled else "DISABLED"
            output += f"\n\nCharm production: {charm_status}"
            if not charm_enabled:
                output += (
                    "\n  ⚠️  Must be enabled for neutrino physics "
                    "and prompt lepton studies"
                )
            output += "\n  (Enable via model._lib.s_cflafr.ipar[14] = 1)"

            # IPAR(11) - Leading rho0 mode
            rho0_mode = int(model._lib.s_cflafr.ipar[10])
            rho0_labels = {
                0: "disabled",
                1: "standard",
                2: "enhanced (Sibyll* +10-15% muons)",
            }
            rho0_desc = rho0_labels.get(rho0_mode, f"unknown mode {rho0_mode}")
            output += f"\nLeading ρ⁰ production: {rho0_desc}"
            output += "\n  (Set via model._lib.s_cflafr.ipar[10] = 0/1/2)"

            # IPAR(12) - Inelastic screening in p-Air
            screening = int(model._lib.s_cflafr.ipar[11])
            screening_status = "ON" if screening else "OFF"
            output += f"\nInelastic screening (p-Air): {screening_status}"
            output += "\n  (Set via model._lib.s_cflafr.ipar[11] = 0/1)"
        except (AttributeError, IndexError):
            # Sibyll 2.1 or other versions may not have these
            pass

        # Nuclear fragmentation mode (for nucleus-nucleus collisions)
        try:
            kodfrag = model._lib.ckfrag.kodfrag
            output += f"\nNuclear fragmentation mode: {kodfrag}"
            output += "\n  (Set via model._lib.ckfrag.kodfrag = value)"
        except AttributeError:
            pass

    if log_level == "verbose":
        # String mass cutoffs from S_CUTOFF (settable via PAR array)
        try:
            str_val = model._lib.s_cutoff.str_mass_val
            str_sea = model._lib.s_cutoff.str_mass_sea
            str_hyp = model._lib.s_cutoff.str_mass_val_hyp
            output += "\n\nString fragmentation masses:"
            output += "\n  Controls minimum invariant mass for string fragmentation"
            output += f"\n  Valence: {str_val:.3f} GeV"
            output += "\n    (model._lib.s_cutoff.str_mass_val = value)"
            output += f"\n  Sea: {str_sea:.3f} GeV"
            output += "\n    (model._lib.s_cutoff.str_mass_sea = value)"
            output += f"\n  Hyperon: {str_hyp:.3f} GeV"
            output += "\n    (model._lib.s_cutoff.str_mass_val_hyp = value)"
        except AttributeError:
            # Sibyll 2.1 has different structure
            pass

        # Fragmentation function parameters from S_CZDIS (settable via PAR array)
        try:
            fa = model._lib.s_czdis.fain
            fb0 = model._lib.s_czdis.fb0in
            output += "\n\nFragmentation function:"
            output += "\n  Controls z-distribution in string fragmentation"
            output += f"\n  FA: {fa:.3f}"
            output += "\n    (model._lib.s_czdis.fain = value)"
            output += f"\n  FB0: {fb0:.3f}"
            output += "\n    (model._lib.s_czdis.fb0in = value)"
        except AttributeError:
            pass

        # Leading baryon parameter from S_CZLEAD (settable via PAR array)
        try:
            clead = model._lib.s_czlead.clead
            output += "\n\nLeading baryon mixing:"
            output += "\n  Controls fraction of hard vs soft leading protons"
            output += f"\n  CLEAD: {clead:.3f}"
            output += "\n    (model._lib.s_czlead.clead = value)"
        except AttributeError:
            pass

    return output


def get_epos_settings(model, log_level="normal"):
    """Get EPOS-specific settings."""
    output = f"\n\n[{model.label} Settings]"

    if log_level in ["normal", "verbose"]:
        # EPOS specific parameters from initialization
        output += f"\nHadronic rescattering: {model._hadronic_rescattering}"

    return output


def get_qgsjet_settings(model, log_level="normal"):
    """Get QGSJet-specific settings."""
    output = f"\n\n[{model.label} Settings]"

    # QGSJet has no user-configurable Fortran settings exposed
    if log_level == "verbose":
        output += "\nNo exposed Fortran settings"

    return output


def get_dpmjet_settings(model, log_level="normal"):
    """Get DPMJet-specific settings."""
    output = f"\n\n[{model.label} Settings]"

    if log_level in ["normal", "verbose"]:
        output += f"\nParameter file: {model._param_file_name}"

    return output


def get_phojet_settings(model, log_level="normal"):
    """Get Phojet-specific settings."""
    output = f"\n\n[{model.label} Settings]"

    if log_level in ["normal", "verbose"]:
        output += f"\nParameter file: {model._param_file_name}"

    return output


def get_urqmd_settings(model, log_level="normal"):
    """Get UrQMD-specific settings."""
    output = f"\n\n[{model.label} Settings]"

    # UrQMD has no user-configurable Fortran settings exposed
    if log_level == "verbose":
        output += "\nNo exposed Fortran settings"

    return output


def get_sophia_settings(model, log_level="normal"):
    """Get Sophia-specific settings."""
    output = f"\n\n[{model.label} Settings]"

    # Sophia has no user-configurable Fortran settings exposed
    if log_level == "verbose":
        output += "\nNo exposed Fortran settings"

    return output


def get_model_settings(model, log_level="normal"):
    """Get model-specific settings as a string."""
    from chromo import models

    # Check model type and dispatch to appropriate function
    model_type = type(model)

    # Pythia8
    if model_type is models.Pythia8:
        return get_pythia8_settings(model, log_level)

    # Pythia6 (different implementation)
    if model_type is models.Pythia6:
        return get_pythia6_settings(model, log_level)

    # Sibyll models
    if model_type in (
        models.Sibyll21,
        models.Sibyll23c,
        models.Sibyll23d,
        models.Sibyll23dStarMixed,
        models.Sibyll23e,
        models.Sibyll23eStarMixed,
    ):
        return get_sibyll_settings(model, log_level)

    # EPOS models
    if model_type in (
        models.EposLHC,
        models.EposLHCR,
        models.EposLHCRHadrRescattering,
    ):
        return get_epos_settings(model, log_level)

    # QGSJet models
    if model_type in (
        models.QGSJet01d,
        models.QGSJetII03,
        models.QGSJetII04,
        models.QGSJetIII,
    ):
        return get_qgsjet_settings(model, log_level)

    # DPMJet models
    if model_type in (
        models.DpmjetIII191,
        models.DpmjetIII193,
        models.DpmjetIII307,
    ):
        return get_dpmjet_settings(model, log_level)

    # Phojet models
    if model_type in (models.Phojet112, models.Phojet191, models.Phojet193):
        return get_phojet_settings(model, log_level)

    # UrQMD models
    if model_type is models.UrQMD34:
        return get_urqmd_settings(model, log_level)

    # Sophia models
    if model_type is models.Sophia20:
        return get_sophia_settings(model, log_level)

    # Fallback for unknown models
    output = f"\n\n[{model.label} Settings]"
    output += f"\nSeed: {model.seed}"
    return output


def write_log(args, p1, p2, model):
    """Generate and write log file with all model information.

    Args:
        args: Command line arguments
        p1: Projectile particle
        p2: Target particle
        model: Model instance (after run completion)
    """
    import sys

    from chromo import __version__

    cmdline = " ".join(sys.argv)
    pr = args.projectile_momentum
    ta = args.target_momentum

    msg = f"""Command:
{cmdline}

[General Information]
Chromo version: {__version__}
Model: {model.label}

[Collision Setup]
Projectile: {p1.name} ({int(args.projectile_id)})
Target: {p2.name} ({int(args.target_id)})
"""
    if pr != 0 or ta != 0:
        msg += f"""Projectile momentum: {pr:g} GeV/c
Target momentum: {ta:g} GeV/c
"""
    msg += f"""sqrt(s): {args.sqrts:g} GeV

[Run Configuration]
Collisions: {args.number}
Seed: {args.seed}
Format: {args.output}
"""
    if args.config:
        msg += f"Configuration file: {args.config}\n"

    msg += get_final_state_particles_info(model)
    msg += get_model_settings(model, args.log_level)

    if args.log_file:
        with open(args.log_file, "w") as f:
            f.write(msg)
