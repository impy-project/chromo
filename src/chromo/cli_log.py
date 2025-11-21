"""
Logging functionality for the command-line interface.
"""


def get_final_state_particles_info(model):
    """Get information about final state particles."""
    from particle import Particle

    if not hasattr(model, "final_state_particles"):
        return ""

    particles = model.final_state_particles
    if not particles:
        return ""

    output = "\n\n[Final State Particles]"
    output += "\nParticles considered stable in this run:\n"

    # Sort by absolute value of PDG ID, with negative first, then positive
    sorted_particles = sorted(particles, key=lambda x: (abs(x), x))

    # Table header
    output += f"\n{'Particle':<20} {'PDG ID':>10}  {'Lifetime (s)':>15}"
    output += f"\n{'-'*20} {'-'*10}  {'-'*15}"

    # Convert to list for error handling
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


def get_model_settings(model, log_level="normal"):
    """Get model-specific settings as a string."""
    output = ""

    if hasattr(model, "_pythia"):
        pythia = model._pythia
        output += "\n\n[Pythia8 Settings]"
        if hasattr(model, "_final_config"):
            output += "\n\nUser configuration:"
            for s in model._final_config:
                output += f"\n  {s}"

        if log_level in ["normal", "verbose"]:
            try:
                import io
                from contextlib import redirect_stdout

                buffer = io.StringIO()
                with redirect_stdout(buffer):
                    pythia.settings.listChanged()
                settings_output = buffer.getvalue()
                if settings_output.strip():
                    output += "\n\nAll changed settings:"
                    output += f"\n{settings_output}"
            except Exception as e:
                output += f"\n\nError getting changed settings: {e}"

        if log_level == "verbose":
            try:
                import io
                from contextlib import redirect_stdout

                buffer = io.StringIO()
                with redirect_stdout(buffer):
                    pythia.settings.listAll()
                settings_output = buffer.getvalue()
                if settings_output.strip():
                    output += "\n\nAll settings (including defaults):"
                    output += f"\n{settings_output}"
            except Exception as e:
                output += f"\n\nError getting all settings: {e}"

    else:
        # For other models, extract available information
        output += f"\n\n[{model.label} Settings]"

        # Check for common model attributes
        if hasattr(model, "seed"):
            output += f"\nSeed: {model.seed}"

        if log_level in ["normal", "verbose"]:
            # Sibyll-specific settings
            if hasattr(model, "_sstar_param"):
                output += f"\nS* parameter (imod): {model._sstar_param}"
            if hasattr(model, "_lib") and hasattr(model._lib, "s_debug"):
                try:
                    debug_level = model._lib.s_debug.ndebug
                    output += f"\nDebug level: {debug_level}"
                except Exception:
                    pass

            # EPOS-specific settings
            if hasattr(model, "_lib") and hasattr(model._lib, "ceposi"):
                try:
                    eposi = model._lib.ceposi
                    if hasattr(eposi, "model"):
                        output += f"\nEPOS model: {eposi.model}"
                    if hasattr(eposi, "nfull"):
                        output += f"\nFull event: {eposi.nfull}"
                except Exception:
                    pass

            # Check for configuration file/parameters
            if hasattr(model, "_param_file"):
                output += f"\nParameter file: {model._param_file}"

            # Add general model info
            if hasattr(model, "_restartable") and model._restartable:
                output += "\nRestartable: Yes"
            else:
                output += "\nRestartable: No"

            if hasattr(model, "_frame"):
                output += f"\nFrame: {model._frame.name}"

    return output


def get_log_text(args, p1, p2, model, include_final_state=True):
    """Generate a log string from the given arguments.

    Args:
        args: Command line arguments
        p1: Projectile particle
        p2: Target particle
        model: Model instance
        include_final_state: If True, include final state particles info.
                           Set to False before model run, True after.
    """
    # Get command line
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

    if include_final_state:
        msg += get_final_state_particles_info(model)
    msg += get_model_settings(model, args.log_level)
    return msg


def write_log_file(log_file, log_text):
    """Write log text to a file."""
    if log_file:
        with open(log_file, "w") as f:
            f.write(log_text)


def log_settings(args, p1, p2, model):
    """Log settings to file and console."""
    from rich.console import Console
    from rich.panel import Panel

    log_text = get_log_text(args, p1, p2, model)
    if args.log_file:
        with open(args.log_file, "w") as f:
            f.write(log_text)

    # convert to rich markup
    msg = log_text.replace(":", "[repr.str]:[/repr.str]")
    for x in (
        "Model",
        "Projectile",
        "Target",
        "Projectile momentum",
        "Target momentum",
        "sqrt(s)",
        "Collisions",
        "Seed",
        "Format",
        "Configuration",
    ):
        msg = msg.replace(f"[repr.str]{x}[/repr.str]", f"[repr.str]{x}[/repr.str]\t")

    from chromo import __version__

    console = Console()
    console.print(
        Panel(
            msg,
            title=f"[bold]chromo [green]{__version__}[/green][/bold]",
            width=78,
            highlight=True,
        )
    )
    return log_text
