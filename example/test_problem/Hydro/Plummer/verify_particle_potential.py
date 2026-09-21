"""
Verify GAMER's interpolated particle potential (OPT__OUTPUT_PAR_POT) against
the analytic Plummer potential:

    Phi(r) = -G*M_tot / sqrt(r**2 + R0**2)

Usage:

    python verify_particle_potential.py [Particle_XXXXXX.txt] [-o out.png]

The input file is a text-format particle dump produced with
OPT__OUTPUT_PAR_MODE=1 (default: Particle_000000.txt in the current
directory). Column layout is assumed to be:

    ParMass ParPosX ParPosY ParPosZ ParVelX ParVelY ParVelZ
    ParAccX ParAccY ParAccZ ParPot ParTime ParType ParPUID ParFlag

Assumes G=1 (GAMER internal units) and a single, non-colliding Plummer
cloud centered on the box center. Adjust R0/RHO0/BOX_SIZE below (or via
CLI flags) to match Input__TestProb / Input__Parameter if they differ.

Note: particles are sampled only out to Plummer_MaxR, so beyond that radius
the *particle* component's enclosed mass is capped while the untruncated
analytic curve -GM_tot/sqrt(r**2+R0**2) keeps adding the (missing) outer
mass. That mismatch is a truncation artifact, not Poisson-solver/
interpolation error, so residuals are only reported/fit for r < max-r.
"""

import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("particle_file", nargs="?", default="Particle_000000.txt",
                         help="text-format particle dump [Particle_000000.txt]")
    parser.add_argument("-o", "--output", default="pot_verification.png",
                         help="output plot filename [pot_verification.png]")
    parser.add_argument("--box-size", type=float, default=3.0,
                         help="BOX_SIZE from Input__Parameter [3.0]")
    parser.add_argument("--r0", type=float, default=0.1,
                         help="Plummer_R0 scale radius [0.1]")
    parser.add_argument("--rho0", type=float, default=1.0,
                         help="Plummer_Rho0 peak density [1.0]")
    parser.add_argument("--fit-radius", type=float, default=0.2,
                         help="only use r < fit-radius to fit the constant "
                              "potential offset [0.2]")
    parser.add_argument("--max-r", type=float, default=0.375,
                         help="Plummer_MaxR: particles beyond this radius don't "
                              "exist, so the untruncated analytic profile is not "
                              "a valid reference there and such particles are "
                              "excluded from the residual statistics [0.375]")
    args = parser.parse_args()

    data = np.loadtxt(args.particle_file, skiprows=2)
    mass = data[:, 0]
    pos = data[:, 1:4]
    pot = data[:, 10]

    box_center = np.full(3, 0.5 * args.box_size)
    r = np.linalg.norm(pos - box_center, axis=1)

    print("N particles:", len(mass))
    print("total particle mass:", mass.sum())
    print("r range:", r.min(), r.max())
    print("Pot range:", pot.min(), pot.max())

    # analytic total mass of an (untruncated) Plummer sphere, G=1
    tot_m_inf = 4.0 / 3.0 * np.pi * args.r0 ** 3 * args.rho0
    print("TotM_Inf (analytic, G=1):", tot_m_inf)
    print("particle mass sum vs TotM_Inf ratio:", mass.sum() / tot_m_inf)

    phi_analytic = -tot_m_inf / np.sqrt(r ** 2 + args.r0 ** 2)

    # GAMER's potential has an arbitrary additive constant (from the Poisson
    # solver's boundary conditions); fit it out using particles well inside
    # the truncation radius, where the analytic profile is most reliable
    inner = r < args.fit_radius
    offset = np.median(pot[inner] - phi_analytic[inner])
    print("fitted constant offset:", offset)

    # particles beyond Plummer_MaxR don't exist, so the untruncated analytic
    # profile systematically overestimates |Phi| there (missing outer mass);
    # exclude them so that truncation error doesn't contaminate the residual
    # statistics for the Poisson-solver/interpolation error we actually want
    valid = r < args.max_r
    n_excluded = np.size(r) - np.count_nonzero(valid)
    if n_excluded > 0:
        print(f"excluding {n_excluded} particle(s) with r >= max-r ({args.max_r}) "
              "from residual statistics (truncation artifact)")

    resid = pot - (phi_analytic + offset)
    rms_inner = np.sqrt(np.mean(resid[inner] ** 2))
    rms_valid = np.sqrt(np.mean(resid[valid] ** 2))
    print("RMS residual (inner r < %.2f): %.3e" % (args.fit_radius, rms_inner))
    print("RMS residual (valid, r < max-r): %.3e" % rms_valid)
    print("relative RMS residual (inner, relative to |phi|): %.2f%%" %
          (100 * rms_inner / np.mean(np.abs(phi_analytic[inner]))))

    order = np.argsort(r)
    r_s = r[order]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7.5), sharex=True,
                                    gridspec_kw={"height_ratios": [2.5, 1]})

    valid_s = r_s < args.max_r

    ax1.scatter(r, pot, s=4, alpha=0.15, color="#1f77b4",
                label="GAMER: interpolated $\\Phi$ at each particle\n"
                      "(OPT__OUTPUT_PAR_POT)")
    ax1.plot(r_s[valid_s], (phi_analytic[order] + offset)[valid_s], color="#d62728", lw=2,
              label=r"Analytic Plummer: $-GM_{\rm tot}/\sqrt{r^2+a^2}$ + const")
    ax1.axvline(args.max_r, color="gray", ls="--", lw=1, label="Plummer_MaxR (truncation)")
    ax1.set_ylabel(r"Gravitational potential $\Phi$ (code units, $G=1$)")
    ax1.set_title("Plummer sphere: GAMER particle-potential verification")
    ax1.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax1.grid(alpha=0.3)

    ax2.scatter(r[valid], resid[valid], s=4, alpha=0.15, color="#2ca02c")
    ax2.scatter(r[~valid], resid[~valid], s=4, alpha=0.15, color="gray")
    ax2.axhline(0.0, color="#d62728", lw=1.5)
    ax2.axvline(args.max_r, color="gray", ls="--", lw=1)
    ax2.set_ylabel(r"Residual $\Phi_{\rm GAMER}-\Phi_{\rm analytic}$")
    ax2.set_xlabel(r"Particle radius $r$ from cluster center (code units)")
    ax2.grid(alpha=0.3)
    ax2.text(0.02, 0.90,
              f"RMS residual (r < max-r) = {rms_valid:.2e}\n"
              f"(relative: {rms_valid / np.mean(np.abs(phi_analytic[valid])):.2%}; "
              f"gray points beyond Plummer_MaxR excluded)",
              transform=ax2.transAxes, va="top", fontsize=9,
              bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print("saved", args.output)


if __name__ == "__main__":
    main()
