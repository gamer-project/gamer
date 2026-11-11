"""Generate core-collapse supernova (CCSN) feedback catalogs.

Units
-----
Stellar masses and gas/metal yields are in Msun, lifetimes in Myr,
energies in erg, and the metallicity Z is a dimensionless mass fraction.
The number IMF xi(m) is normalized per unit formed stellar mass by
integral[m*xi(m) dm] = 1 over [0.08, 120] Msun.

Sampling and caveats
--------------------
The randomly sampled catalog draws independent Poisson-distributed star
counts in logarithmic mass bins, and random masses within each bin.
The total initial mass is therefore random: its expectation equals
target_particle_mass, but individual realizations differ.

The fixed-count SN catalogs are deterministic.
For k events, the conditional SN number IMF is cut into k slices of
equal expected number, and each slice is represented by its mean mass.
They are reproducible, but they are not random mass realizations, and
the explosion-time distribution is resolved only at finite resolution.
By default, the gas and metal yields are evaluated at the slice-mean mass,
which is biased high because the yields are nonlinear in mass.
The option slice_averaged_yields (--slice-averaged-yields) stores the
slice-averaged yields instead, which removes this bias.
The catalog with the largest count also represents the remaining upper
tail of the Poisson distribution.

After tabulation, the Poisson-weighted totals are compared with the
IMF-integrated expectations, and a warning is issued if any of them
differs by more than 1 percent.

The fixed-count grid stores k events for every k up to the maximum count, so
its memory use and run time grow quadratically with the expected SN count.
If more than 1000 SNe are expected, the fixed-count catalogs are skipped
with a warning.
"""

import argparse
from functools import lru_cache
import os
import sys
import warnings
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
from scipy.special import erf, erfinv
from scipy.stats import poisson

# Scalar or array argument of the vectorized IMF and yield functions.
FloatOrArray = Union[float, np.ndarray]

# Metallicity mass fractions over which the stellar-lifetime fit is used.
Z_MIN = 0.0004
Z_MAX = 0.03

# Conservative guard on the expected all-star sample size. A float64 mass
# array needs 8 bytes per star, and bin arrays, temporary copies and plots
# increase peak memory substantially. This is not a strict limit on the
# stochastic realized count; users needing larger populations should use
# an aggregate/streaming sampler rather than materializing every star.
MAX_SAMPLED_STARS = 2_000_000

# Smallest value drawn on a logarithmic plot axis. Empty histogram bins are
# raised to it, because zero cannot be shown on a log scale.
LOG_PLOT_FLOOR = 1.0e-30

class StarParticle:
    """Represent a coeval stellar population and its CCSN feedback.

    Parameters
    ----------
    target_particle_mass : float, optional
        Expected mass of the population in Msun, not a realized hard budget.
    random_seed : int or None, optional
        Seed for the per-instance numpy.random.Generator.
    imf_type : {'Kroupa', 'Chabrier'}, optional
        Initial mass function used for mass sampling and expectations.
    stellar_metallicity : float, optional
        Stellar metallicity for the lifetime formula; default Z=0.02.
        Only stellar lifetimes depend on metallicity; yields do not.
        The accepted range is [0.0004, 0.03].
    slice_averaged_yields : bool, optional
        If True, each event of a fixed-count catalog carries the explosion
        time and the gas and metal yields averaged over its equal-number IMF
        slice instead of the values evaluated at the slice-mean mass.
        Default False.

    Attributes
    ----------
    target_particle_mass : float
        Nominal formed initial mass (Msun).
    slice_averaged_yields : bool
        Whether the fixed-count catalogs use slice-averaged yields.
    random_seed : int or None
        Recorded seed; it does not change deterministic fixed-count nodes.
    imf_type : str
        Selected number IMF.
    stellar_metallicity : float
        Metallicity passed into the stellar lifetime model.
    rng : numpy.random.Generator
        Independent pseudorandom generator for Poisson and uniform draws.
    M_MIN, M_MAX : float
        IMF mass range in Msun, 0.08 to 120.
    KROUPA_M_BREAK, KROUPA_ALPHA_LOWER, KROUPA_ALPHA_UPPER : float
        Kroupa break mass (Msun) and low/high number-IMF power-law slopes.
    CHABRIER_M_C, CHABRIER_LOG10_SIGMA : float
        Lognormal location (Msun) and log10-mass width.
    CHABRIER_M_BREAK, CHABRIER_ALPHA_UPPER : float
        Chabrier transition mass (Msun) and high-mass power-law exponent.
    M_SN_MIN, M_SN_YIELD_BREAK, M_SN_MAX : float
        Eligible progenitor bounds (8 to 120 Msun) and 40-Msun yield switch.
    KROUPA_NORM_LOWER, KROUPA_NORM_UPPER : float
        Mass-normalized continuous Kroupa branch coefficients.
    CHABRIER_NORM_LOWER, CHABRIER_NORM_UPPER : float
        Mass-normalized continuous Chabrier branch coefficients.
    sampled_star_masses : numpy.ndarray
        Realized star masses (Msun), ordered by their logarithmic source bins.
    sampled_stellar_mass : float
        Sum of realized masses (Msun), which fluctuates about nominal mass.
    sampled_star_count : int
        Realized count of all stars, including low-mass ones.
    expected_star_count, expected_sn_count : float
        Number-IMF expectations for the full IMF and CCSN mass range.
    expected_stellar_mass, expected_sn_progenitor_mass : float
        Mean total initial mass and mean initial SN-progenitor mass (Msun).
    expected_sn_energy : float
        Mass-dependent IMF-integrated SN energy expectation (erg).
    expected_ejected_gas_mass, expected_ejected_metal_mass : float
        Mass-dependent IMF-integrated gas and metal expectations (Msun).
    expected_sn_mean_time : float
        Number-IMF-averaged explosion time of one SN (Myr);
        independent of the particle mass.
    agora_expected_sn_count, agora_expected_sn_energy : float
        Independent AGORA-style SN count and energy expectations.
    agora_expected_ejected_gas_mass, agora_expected_ejected_metal_mass : float
        AGORA-style expected ejected gas and metal masses (Msun).
    sampled_sn_events : list of dict
        Time-sorted sampled SN events with five physical fields.
    sampled_sn_count : int
        Realized SN count in the all-star catalog.
    sampled_first_sn_time : float
        First sampled explosion time (Myr), or NaN for zero SNe.
    sampled_last_sn_time : float
        Last sampled explosion time (Myr), or NaN for zero SNe.
    sampled_sn_progenitor_mass : float
        Initial mass of all sampled SN progenitors (Msun).
    sampled_sn_energy : float
        Total sampled SN explosion energy (erg).
    sampled_ejected_gas_mass : float
        Total sampled gas ejecta (Msun).
    sampled_ejected_metal_mass : float
        Total sampled metal ejecta (Msun).
    max_sn_count : int
        Largest deterministic count k,
        max(2, ceil(lambda + 5*sqrt(lambda))),
        with lambda = expected_sn_count.
        Its catalog also stands for every larger count
        (the residual Poisson upper tail).
        It is 0 if the fixed-count grid was skipped.
    fixed_count_sn_events : list of list of dict
        k-th list holds k slice-mean progenitors (with slice-averaged yields
        and explosion times if slice_averaged_yields is True,
        else yields at the mean masses); list zero is empty.
        If the fixed-count grid was skipped (expected_sn_count above 1000),
        only list zero exists.
    fixed_count_sn_counts : numpy.ndarray
        Fixed event counts, including zero, indexed by count index k.
    fixed_count_first_sn_times : numpy.ndarray
        Earliest explosion time by count (Myr); zero-count entry is NaN.
    fixed_count_last_sn_times : numpy.ndarray
        Latest explosion time by count (Myr); zero-count entry is NaN.
    fixed_count_progenitor_masses : numpy.ndarray
        Total initial SN-progenitor mass by count (Msun).
    fixed_count_energies : numpy.ndarray
        Summed explosion energy by count (erg).
    fixed_count_gas_masses : numpy.ndarray
        Summed gas ejecta by count (Msun).
    fixed_count_metal_masses : numpy.ndarray
        Summed metal ejecta by count (Msun).
    fixed_count_poisson_weights : numpy.ndarray
        Count probabilities; final bin includes the residual upper tail.
        If the grid was skipped, the array holds only entry zero, P(N=0).
    poisson_weighted_sn_count : float
        Poisson-weighted tabulated mean count.
    poisson_weighted_progenitor_mass : float
        Poisson-weighted initial progenitor mass (Msun).
    poisson_weighted_sn_energy : float
        Poisson-weighted ejected energy (erg).
    poisson_weighted_ejected_gas_mass : float
        Poisson-weighted ejected gas (Msun).
    poisson_weighted_ejected_metal_mass : float
        Poisson-weighted ejected metals (Msun).
    poisson_weighted_mean_sn_time : float
        Poisson-weighted mean explosion time per SN in the fixed-count
        catalogs (Myr); zero if the grid was skipped.

    Notes
    -----
    The randomly sampled catalog and the deterministic fixed-count catalogs
    serve different purposes.
    A single sampled population does not contain exactly target_particle_mass.
    Every eligible massive star is assumed to explode successfully as a CCSN
    with the prescribed yields.

    If the expected SN count exceeds 1000, the fixed-count catalogs are not
    built (a warning is issued): max_sn_count is 0,
    all poisson_weighted_* attributes remain zero, and only the randomly
    sampled catalog carries information.
    The summary, table, and plotting methods still run, but the fixed-count
    tables and curves are empty.
    """

    def __init__(self,
                 target_particle_mass: float = 1000.0,
                 random_seed: Optional[int] = 123,
                 imf_type: str = 'Kroupa',
                 stellar_metallicity: float = 0.02,
                 slice_averaged_yields: bool = False
                ) -> None:
        """Initialize a coeval star population and its feedback catalogs.

        Parameters
        ----------
        target_particle_mass : float, optional
            Expected initial stellar mass, in solar masses (Msun).
            The sampled population is not constrained to contain exactly this mass.
        random_seed : int or None, optional
            Seed for this object's independent NumPy Generator.
            Identical seeds and settings reproduce a sample in the same
            software environment.
        imf_type : {'Kroupa', 'Chabrier'}, optional
            IMF selected for sampling and integrating stellar populations.
        stellar_metallicity : float, optional
            Stellar metallicity mass fraction Z used in event lifetimes.
            The default is 0.02.
            It does not modify the chosen yield prescription.
        slice_averaged_yields : bool, optional
            Use slice-averaged explosion times and gas and metal yields
            in the fixed-count catalogs
            (see _fixed_count_slice_yields and _fixed_count_slice_times).
            The default, False, evaluates the yields at the slice-mean
            masses.
            The sampled catalog and all IMF-integrated expectations
            are unaffected.

        Raises
        ------
        ValueError
            If target_particle_mass is not finite and positive, if
            stellar_metallicity is outside [0.0004, 0.03], if imf_type is not
            'Kroupa' or 'Chabrier', or if slice_averaged_yields is not a bool.

        Notes
        -----
        The constructor normalizes both IMFs, samples all stars, computes the
        expectations, and builds the sampled and fixed-count event lists.
        The fixed-count lists are skipped, with a warning, if the expected SN
        count exceeds 1000.
        The inputs are validated before any calculation, except the seed,
        which NumPy validates itself.
        A very small population can contain no sampled stars or SNe.
        """

        if not np.isfinite(target_particle_mass) or target_particle_mass <= 0:
            raise ValueError(f"target_particle_mass must be a finite positive float, got {target_particle_mass}")
        if not np.isfinite(stellar_metallicity) or not Z_MIN <= stellar_metallicity <= Z_MAX:
            raise ValueError(f"stellar_metallicity must lie in [{Z_MIN:g}, {Z_MAX:g}], got {stellar_metallicity}")
        if imf_type not in ("Kroupa", "Chabrier"):
            raise ValueError(f"Unknown imf_type: {imf_type}")
        if not isinstance(slice_averaged_yields, (bool, np.bool_)):
            raise ValueError(f"slice_averaged_yields must be a bool, got {slice_averaged_yields!r}")

        self.target_particle_mass = target_particle_mass
        self.random_seed = random_seed
        self.imf_type = imf_type
        self.stellar_metallicity = stellar_metallicity
        self.slice_averaged_yields = bool(slice_averaged_yields)

        self.rng = np.random.default_rng(self.random_seed)

        # IMF mass limits.
        # Lower mass limit of the IMF [Msun].
        self.M_MIN = 0.08
        # Upper mass limit of the IMF [Msun].
        self.M_MAX = 120.0

        # Kroupa IMF parameters (Kroupa 2001).
        # Break mass between the two power-law segments [Msun].
        self.KROUPA_M_BREAK = 0.5
        # Number-IMF power-law index below the break mass.
        self.KROUPA_ALPHA_LOWER = 1.3
        # Number-IMF power-law index above the break mass.
        self.KROUPA_ALPHA_UPPER = 2.3

        # Chabrier IMF parameters (Chabrier 2003).
        # Characteristic mass of the log-normal segment [Msun].
        self.CHABRIER_M_C = 0.08
        # Standard deviation of log10(m) in the log-normal segment.
        self.CHABRIER_LOG10_SIGMA = 0.69
        # Transition mass from the log-normal to the power-law segment [Msun].
        self.CHABRIER_M_BREAK = 1.0
        # Number-IMF power-law index above the transition mass.
        self.CHABRIER_ALPHA_UPPER = 2.3

        # CCSN progenitor limits and the transition in the yield model.
        # Minimum initial mass of a CCSN progenitor [Msun].
        self.M_SN_MIN = 8.0
        # Initial mass above which the adopted gas and metal yields are constant [Msun].
        self.M_SN_YIELD_BREAK = 40.0
        # Maximum initial mass of a CCSN progenitor [Msun].
        self.M_SN_MAX = self.M_MAX
        # Reference compact remnant mass below the yield break [Msun].
        # Above the break, the gas yield saturates and the unreturned mass
        # rises with progenitor mass; it is not a fixed 1.4-Msun remnant.
        self.COMPACT_REMNANT_MASS = 1.4

        self._normalize_imfs()
        self._sample_population()
        self._compute_expectations()
        self._build_sampled_sn_events()
        self._build_fixed_count_sn_events()

    def _integral(self, f: Callable, a: FloatOrArray, b: FloatOrArray, breaks: Sequence[float] = ()) -> FloatOrArray:
        """Integrate a scalar function over scalar or broadcastable bounds.

        Parameters
        ----------
        f : callable
            Scalar integrand f(m), evaluated by scipy.integrate.quad.
        a : float or numpy.ndarray
            Lower mass bound; may be broadcast with b.
        b : float or numpy.ndarray
            Upper mass bound; may be broadcast with a.
        breaks : sequence of float, optional
            Masses (Msun) where the integrand is not smooth, such as an IMF
            branch boundary or the yield break.
            Those lying strictly inside (a, b) are passed to quad as
            interior points.
            Default: none.

        Returns
        -------
        float or numpy.ndarray
            Integral for each pair of bounds. Identical bounds give zero.

        Notes
        -----
        Adaptive quadrature uses epsabs=0, epsrel=1e-10, and limit=200.
        The caller supplies the breakpoints of its own integrand, so
        this method does not depend on the selected IMF.
        A warning is issued if the absolute error estimated by quad is
        large compared with the integral.
        For array bounds, np.vectorize calls quad serially.
        """

        breaks = tuple(breaks)

        # Quadrature tolerances: relative error only, with no absolute floor.
        epsabs, epsrel = 0.0, 1.0e-10

        def integrate_bounds(lower: float, upper: float) -> float:
            """Integrate one pair of scalar bounds, splitting at IMF breaks.

            Parameters
            ----------
            lower : float
                Lower integration bound.
            upper : float
                Upper integration bound.

            Returns
            -------
            float
                Adaptive-quadrature estimate between lower and upper.

            Notes
            -----
            Reversed bounds are passed to quad without a separate check.
            """

            if lower == upper:
                return 0.0

            # Breakpoints strictly inside the interval, for quad to subdivide at.
            pts = [p for p in breaks if lower < p < upper and not np.isclose(p, lower) and not np.isclose(p, upper)]

            val, err = scipy.integrate.quad(f, lower, upper,
                                            points=pts if pts else None,
                                            epsabs=epsabs,
                                            epsrel=epsrel,
                                            limit=200)
            if err > 1.0e-6 * abs(val) and err > 1.0e-14:
                warnings.warn(f"quad error {err:.2e} vs value {val:.2e} on [{lower}, {upper}]")

            return val

        if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
            return np.vectorize(integrate_bounds, otypes=[float])(a, b)

        return integrate_bounds(a, b)

    def _integral_power_law(self, alpha: float, a: FloatOrArray, b: FloatOrArray) -> FloatOrArray:
        """Calculate integral_{a}^{b} m**(-alpha) dm analytically.

        Parameters
        ----------
        alpha : float
            Power-law index of the integrand.
        a : float or numpy.ndarray
            Lower mass bound; may be broadcast with b.
        b : float or numpy.ndarray
            Upper mass bound; may be broadcast with a.

        Returns
        -------
        float or numpy.ndarray
            Integral for each pair of bounds.
            Identical bounds give zero.
            Bounds must be positive; close endpoints are evaluated
            with cancellation-resistant functions.
        """

        # The antiderivative is m**p / p, or ln(m) for alpha = 1.
        p = 1.0 - alpha

        # Compute differences near a == b with log1p/expm1 rather than
        # subtracting two nearly equal powers (catastrophic cancellation).
        # All IMF mass bounds are strictly positive.
        a_arr, b_arr = np.broadcast_arrays(np.asarray(a, dtype=float),
                                           np.asarray(b, dtype=float))
        log_ratio = np.log1p((b_arr - a_arr) / a_arr)

        if np.isclose(p, 0.0):
            integral = log_ratio
        else:
            integral = np.power(a_arr, p) * np.expm1(p * log_ratio) / p

        integral = np.where(a_arr == b_arr, 0.0, integral)
        return float(integral) if integral.ndim == 0 else integral

    def _raw_kroupa_lower_imf(self, m: FloatOrArray) -> FloatOrArray:
        """Evaluate unnormalized low-mass Kroupa number-IMF power law.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun, possibly an array.

        Returns
        -------
        float or numpy.ndarray
            m to the negative low-mass IMF exponent.

        Raises
        ------
        ValueError
            If any mass is not strictly positive.

        Notes
        -----
        Used for masses below the Kroupa break.
        """

        if not np.all(np.asarray(m) > 0.0):
            raise ValueError("Stellar masses must be positive")

        return m**(-self.KROUPA_ALPHA_LOWER)

    def _raw_kroupa_upper_imf(self, m: FloatOrArray) -> FloatOrArray:
        """Evaluate unnormalized high-mass Kroupa number-IMF power law.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun, possibly an array.

        Returns
        -------
        float or numpy.ndarray
            m to the negative high-mass IMF exponent.

        Raises
        ------
        ValueError
            If any mass is not strictly positive.

        Notes
        -----
        Used for masses above the Kroupa break.
        """

        if not np.all(np.asarray(m) > 0.0):
            raise ValueError("Stellar masses must be positive")

        return m**(-self.KROUPA_ALPHA_UPPER)

    def _raw_chabrier_lower_imf(self, m: FloatOrArray) -> FloatOrArray:
        """Evaluate unnormalized low-mass Chabrier number IMF per unit mass.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun, possibly an array.

        Returns
        -------
        float or numpy.ndarray
            Lognormal-in-log10(m) kernel divided by m.

        Raises
        ------
        ValueError
            If any mass is not strictly positive.

        Notes
        -----
        m must be positive; division by m converts from logarithmic mass to
        linear mass.
        """

        if not np.all(np.asarray(m) > 0.0):
            raise ValueError("Stellar masses must be positive")

        return np.exp(-(np.log10(m / self.CHABRIER_M_C)**2) / (2 * self.CHABRIER_LOG10_SIGMA**2)) / m

    def _raw_chabrier_upper_imf(self, m: FloatOrArray) -> FloatOrArray:
        """Evaluate unnormalized high-mass Chabrier number-IMF power law.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun, possibly an array.

        Returns
        -------
        float or numpy.ndarray
            m raised to minus CHABRIER_ALPHA_UPPER.

        Raises
        ------
        ValueError
            If any mass is not strictly positive.

        Notes
        -----
        Used for masses above the Chabrier break.
        """

        if not np.all(np.asarray(m) > 0.0):
            raise ValueError("Stellar masses must be positive")

        return m**(-self.CHABRIER_ALPHA_UPPER)

    def _normalize_kroupa_imf(self) -> None:
        """Normalize both Kroupa IMF branches per unit formed stellar mass.

        Returns
        -------
        None
            Sets KROUPA_NORM_LOWER and KROUPA_NORM_UPPER.

        Notes
        -----
        For number IMF xi(m), integral(m * xi(m) dm) = 1 Msun.
        """

        continuity_factor = self._raw_kroupa_lower_imf(self.KROUPA_M_BREAK) / self._raw_kroupa_upper_imf(self.KROUPA_M_BREAK)

        # Unnormalized mass contained in each branch; their sum is rescaled to 1 Msun.
        raw_total_mass_lower = self._integral_power_law(self.KROUPA_ALPHA_LOWER - 1.0, self.M_MIN, self.KROUPA_M_BREAK)
        raw_total_mass_upper = continuity_factor * self._integral_power_law(self.KROUPA_ALPHA_UPPER - 1.0, self.KROUPA_M_BREAK, self.M_MAX)
        # or, raw_total_mass_lower = self._integral(lambda m: self._raw_kroupa_lower_imf(m)*m, self.M_MIN, self.KROUPA_M_BREAK,
        #                                           breaks=(self.KROUPA_M_BREAK,))
        #     raw_total_mass_upper = self._integral(lambda m: continuity_factor*self._raw_kroupa_upper_imf(m)*m, self.KROUPA_M_BREAK, self.M_MAX,
        #                                           breaks=(self.KROUPA_M_BREAK,))

        self.KROUPA_NORM_LOWER = 1.0 / (raw_total_mass_lower + raw_total_mass_upper)
        self.KROUPA_NORM_UPPER = continuity_factor * self.KROUPA_NORM_LOWER

        return

    def _normalize_chabrier_imf(self) -> None:
        """Normalize both Chabrier IMF branches per unit formed stellar mass.

        Returns
        -------
        None
            Sets CHABRIER_NORM_LOWER and CHABRIER_NORM_UPPER.

        Notes
        -----
        For number IMF xi(m), integral(m * xi(m) dm) = 1 Msun.
        The low-mass number-density includes an explicit 1/m Jacobian.
        """

        # Scale factor that makes the power-law branch continuous with the
        # log-normal branch at the transition mass.
        continuity_factor = self._raw_chabrier_lower_imf(self.CHABRIER_M_BREAK) / self._raw_chabrier_upper_imf(self.CHABRIER_M_BREAK)

        # Unnormalized mass contained in each branch; their sum is rescaled to 1 Msun.
        raw_total_mass_lower = self._integral(lambda m: self._raw_chabrier_lower_imf(m) * m, self.M_MIN, self.CHABRIER_M_BREAK,
                                              breaks=(self.CHABRIER_M_BREAK,))
        raw_total_mass_upper = continuity_factor * self._integral_power_law(self.CHABRIER_ALPHA_UPPER - 1.0, self.CHABRIER_M_BREAK, self.M_MAX)
        # or, raw_total_mass_upper = self._integral(lambda m: continuity_factor*self._raw_chabrier_upper_imf(m)*m, self.CHABRIER_M_BREAK, self.M_MAX,
        #                                           breaks=(self.CHABRIER_M_BREAK,))

        self.CHABRIER_NORM_LOWER = 1.0 / (raw_total_mass_lower + raw_total_mass_upper)
        self.CHABRIER_NORM_UPPER = continuity_factor * self.CHABRIER_NORM_LOWER

        return

    def _normalize_imfs(self) -> None:
        """Initialize coefficients for both supported initial mass functions.

        Returns
        -------
        None
            Sets both branch normalization coefficients for each IMF.

        Notes
        -----
        Both IMFs are normalized even when only one is selected.
        """

        self._normalize_kroupa_imf()
        self._normalize_chabrier_imf()

        return

    def kroupa_number_imf(self, m: FloatOrArray) -> FloatOrArray:
        """Evaluate the Kroupa IMF normalized to 1 Msun.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun, possibly an array.

        Returns
        -------
        xi(m), a stellar number density per Msun formed.

        Notes
        -----
        The integral of m*xi over M_MIN to M_MAX is one; xi itself does not
        integrate to one.
        Input masses outside the IMF interval return zero.
        A clipped evaluation mass prevents invalid intermediate
        power-law values.
        """

        # Fast path that avoids NumPy overhead in scalar quadrature calls.
        if np.isscalar(m):
            if not (self.M_MIN <= m <= self.M_MAX):
                return 0.0
            return self.KROUPA_NORM_LOWER * self._raw_kroupa_lower_imf(m) if m < self.KROUPA_M_BREAK else self.KROUPA_NORM_UPPER * self._raw_kroupa_upper_imf(m)

        # np.where evaluates both branches, so clip the mass to keep the
        # power laws in their domain; xi is zero outside [M_MIN, M_MAX].
        m_eval = np.maximum(m, self.M_MIN)
        valid_mask = (m >= self.M_MIN) & (m <= self.M_MAX)
        xi = np.where(valid_mask,
                      np.where((m < self.KROUPA_M_BREAK),
                               self.KROUPA_NORM_LOWER * self._raw_kroupa_lower_imf(m_eval),
                               self.KROUPA_NORM_UPPER * self._raw_kroupa_upper_imf(m_eval)),
                      0.0)

        return xi

    def chabrier_number_imf(self, m: FloatOrArray) -> FloatOrArray:
        """Evaluate the Chabrier IMF normalized to 1 Msun.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun, possibly an array.

        Returns
        -------
        xi(m), a stellar number density per Msun formed.

        Notes
        -----
        The integral of m*xi over M_MIN to M_MAX is one; xi itself does not
        integrate to one.
        Input masses outside the IMF interval return zero.
        A clipped evaluation mass prevents log10 of nonpositive values.
        """

        # Fast path that avoids NumPy overhead in scalar quadrature calls.
        if np.isscalar(m):
            if not (self.M_MIN <= m <= self.M_MAX):
                return 0.0
            return self.CHABRIER_NORM_LOWER * self._raw_chabrier_lower_imf(m) if m < self.CHABRIER_M_BREAK else self.CHABRIER_NORM_UPPER * self._raw_chabrier_upper_imf(m)

        # np.where evaluates both branches, so clip the mass to avoid
        # log10(m) at m <= 0; xi is zero outside [M_MIN, M_MAX].
        m_eval = np.maximum(m, self.M_MIN)
        valid_mask = (m >= self.M_MIN) & (m <= self.M_MAX)
        xi = np.where(valid_mask,
                      np.where((m < self.CHABRIER_M_BREAK),
                               self.CHABRIER_NORM_LOWER * self._raw_chabrier_lower_imf(m_eval),
                               self.CHABRIER_NORM_UPPER * self._raw_chabrier_upper_imf(m_eval)),
                      0.0)

        return xi

    def kroupa_mass_weighted_imf(self, m: FloatOrArray) -> FloatOrArray:
        """Evaluate mass density m*xi(m) of the normalized Kroupa IMF.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun, possibly an array.

        Returns
        -------
        Fraction of formed mass per unit stellar-mass interval.

        Notes
        -----
        Integrates to one across the full modeled IMF range.
        """

        return m * self.kroupa_number_imf(m)

    def chabrier_mass_weighted_imf(self, m: FloatOrArray) -> FloatOrArray:
        """Evaluate mass density m*xi(m) of the normalized Chabrier IMF.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun, possibly an array.

        Returns
        -------
        Fraction of formed mass per unit stellar-mass interval.

        Notes
        -----
        Integrates to one across the full modeled IMF range.
        """

        return m * self.chabrier_number_imf(m)

    def _kroupa_closed_form_cumulative_number(self, m: FloatOrArray) -> FloatOrArray:
        """Evaluate the closed-form cumulative stellar number of the Kroupa IMF.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun; clipped to [M_MIN, M_MAX].

        Returns
        -------
        float or numpy.ndarray
            Number of stars per Msun of formed stellar mass with initial mass
            below m, i.e., the integral of xi from M_MIN to m (not multiplied
            by target_particle_mass).

        Notes
        -----
        Adds the analytic power-law integrals of the two Kroupa branches.
        Used by _expected_star_count.
        """

        m = np.clip(np.asarray(m, dtype=np.float64), self.M_MIN, self.M_MAX)
        number_segment_lower = self.KROUPA_NORM_LOWER * self._integral_power_law(self.KROUPA_ALPHA_LOWER, self.M_MIN, np.minimum(m, self.KROUPA_M_BREAK))
        number_segment_upper = self.KROUPA_NORM_UPPER * self._integral_power_law(self.KROUPA_ALPHA_UPPER, self.KROUPA_M_BREAK, np.maximum(m, self.KROUPA_M_BREAK))

        return number_segment_lower + number_segment_upper

    def _kroupa_closed_form_cumulative_mass(self, m: FloatOrArray) -> FloatOrArray:
        """Evaluate the closed-form cumulative stellar mass of the Kroupa IMF.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun; clipped to [M_MIN, M_MAX].

        Returns
        -------
        float or numpy.ndarray
            Mass of stars per Msun of formed stellar mass with initial mass
            below m, i.e., the integral of m*xi from M_MIN to m (not multiplied
            by target_particle_mass).

        Notes
        -----
        Adds the analytic power-law integrals of the two Kroupa branches.
        Used by _expected_stellar_mass.
        """

        m = np.clip(np.asarray(m, dtype=np.float64), self.M_MIN, self.M_MAX)
        mass_segment_lower = self.KROUPA_NORM_LOWER * self._integral_power_law(self.KROUPA_ALPHA_LOWER - 1.0, self.M_MIN, np.minimum(m, self.KROUPA_M_BREAK))
        mass_segment_upper = self.KROUPA_NORM_UPPER * self._integral_power_law(self.KROUPA_ALPHA_UPPER - 1.0, self.KROUPA_M_BREAK, np.maximum(m, self.KROUPA_M_BREAK))

        return mass_segment_lower + mass_segment_upper

    def _chabrier_closed_form_cumulative_number(self, m: FloatOrArray) -> FloatOrArray:
        """Evaluate the closed-form cumulative stellar number of the Chabrier IMF.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun; clipped to [M_MIN, M_MAX].

        Returns
        -------
        float or numpy.ndarray
            Number of stars per Msun of formed stellar mass with initial mass
            below m, i.e., the integral of xi from M_MIN to m (not multiplied
            by target_particle_mass).

        Notes
        -----
        The lognormal branch is integrated in log10(m) through the error
        function, and the power-law branch analytically.
        Used by _expected_star_count.
        """

        m = np.clip(np.asarray(m, dtype=np.float64), self.M_MIN, self.M_MAX)
        c_loc = np.log10(self.CHABRIER_M_C)
        c_scale = np.sqrt(2.0) * self.CHABRIER_LOG10_SIGMA

        def erf_of_log10_mass(x: np.ndarray) -> np.ndarray:
            """Evaluate the error function of the standardized log10(m / Msun).

            Parameters
            ----------
            x : numpy.ndarray
                Base-10 logarithm of the stellar mass in Msun.

            Returns
            -------
            numpy.ndarray
                erf((x - log10(M_C)) / (sqrt(2) * sigma)).
            """

            return erf((x - c_loc) / c_scale)

        number_segment_lower = (self.CHABRIER_NORM_LOWER * np.log(10.0) * self.CHABRIER_LOG10_SIGMA * np.sqrt(np.pi / 2.0)
                                * (erf_of_log10_mass(np.log10(np.minimum(m, self.CHABRIER_M_BREAK))) - erf_of_log10_mass(np.log10(self.M_MIN))))
        number_segment_upper = self.CHABRIER_NORM_UPPER * self._integral_power_law(self.CHABRIER_ALPHA_UPPER, self.CHABRIER_M_BREAK, np.maximum(m, self.CHABRIER_M_BREAK))

        return number_segment_lower + number_segment_upper

    def _chabrier_closed_form_cumulative_mass(self, m: FloatOrArray) -> FloatOrArray:
        """Evaluate the closed-form cumulative stellar mass of the Chabrier IMF.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun; clipped to [M_MIN, M_MAX].

        Returns
        -------
        float or numpy.ndarray
            Mass of stars per Msun of formed stellar mass with initial mass
            below m, i.e., the integral of m*xi from M_MIN to m (not multiplied
            by target_particle_mass).

        Notes
        -----
        With x = log10(m / M_C), a = ln(10) and sigma = CHABRIER_LOG10_SIGMA,
        the lognormal branch is
        m*xi dm = NORM_LOWER * exp(-x^2 / (2 sigma^2)) * M_C * a * 10^x dx.
        Completing the square in the exponent, a*x - x^2 / (2 sigma^2)
        = -(x - a sigma^2)^2 / (2 sigma^2) + a^2 sigma^2 / 2,
        reduces it to an error function.
        The power-law branch is analytic.
        Used by _expected_stellar_mass.
        """

        m = np.clip(np.asarray(m, dtype=np.float64), self.M_MIN, self.M_MAX)
        a = np.log(10.0)
        sigma = self.CHABRIER_LOG10_SIGMA
        c_loc = np.log10(self.CHABRIER_M_C)

        def erf_shifted(x: np.ndarray) -> np.ndarray:
            """Evaluate erf((x - a sigma^2) / (sqrt(2) sigma)) for x = log10(m / M_C).

            Parameters
            ----------
            x : numpy.ndarray
                Base-10 logarithm of the stellar mass in Msun, minus log10(M_C).

            Returns
            -------
            numpy.ndarray
                Error function of the shifted, standardized argument.
            """

            return erf((x - a * sigma**2) / (np.sqrt(2.0) * sigma))

        prefactor = self.CHABRIER_NORM_LOWER * self.CHABRIER_M_C * a * np.exp(0.5 * (a * sigma)**2) * sigma * np.sqrt(np.pi / 2.0)
        mass_segment_lower = prefactor * (erf_shifted(np.log10(np.minimum(m, self.CHABRIER_M_BREAK)) - c_loc)
                                          - erf_shifted(np.log10(self.M_MIN) - c_loc))
        mass_segment_upper = self.CHABRIER_NORM_UPPER * self._integral_power_law(self.CHABRIER_ALPHA_UPPER - 1.0, self.CHABRIER_M_BREAK, np.maximum(m, self.CHABRIER_M_BREAK))

        return mass_segment_lower + mass_segment_upper

    def _inverse_closed_form_cumulative_number(self, cumulative_number: np.ndarray) -> np.ndarray:
        """Invert the closed-form cumulative stellar number of the selected IMF.

        Parameters
        ----------
        cumulative_number : numpy.ndarray
            Values of the cumulative number of stars per Msun of formed mass,
            as returned by _kroupa_closed_form_cumulative_number or
            _chabrier_closed_form_cumulative_number, within their range.

        Returns
        -------
        numpy.ndarray
            Stellar masses in Msun, clipped to [M_MIN, M_MAX], with the shape
            of cumulative_number.

        Notes
        -----
        The inverse is exact and piecewise. Each power-law branch is a
        power of m, and the Chabrier lognormal branch is an error function
        whose inverse is scipy.special.erfinv:
        log10(m) = log10(M_C) + sqrt(2) sigma erfinv(erf(z_min) + N / A), with
        A = NORM_LOWER ln(10) sigma sqrt(pi/2) the amplitude used in
        _chabrier_closed_form_cumulative_number.
        """

        n = np.asarray(cumulative_number, dtype=np.float64)

        if self.imf_type == 'Kroupa':
            n_break = self._kroupa_closed_form_cumulative_number(self.KROUPA_M_BREAK)
            p_lower = 1.0 - self.KROUPA_ALPHA_LOWER
            p_upper = 1.0 - self.KROUPA_ALPHA_UPPER
            m_lower = (self.M_MIN**p_lower + p_lower * np.minimum(n, n_break) / self.KROUPA_NORM_LOWER)**(1.0 / p_lower)
            m_upper = (self.KROUPA_M_BREAK**p_upper + p_upper * np.maximum(n - n_break, 0.0) / self.KROUPA_NORM_UPPER)**(1.0 / p_upper)

        else:
            n_break = self._chabrier_closed_form_cumulative_number(self.CHABRIER_M_BREAK)
            p_upper = 1.0 - self.CHABRIER_ALPHA_UPPER
            sigma = self.CHABRIER_LOG10_SIGMA
            c_loc = np.log10(self.CHABRIER_M_C)
            amplitude = self.CHABRIER_NORM_LOWER * np.log(10.0) * sigma * np.sqrt(np.pi / 2.0)
            erf_min = erf((np.log10(self.M_MIN) - c_loc) / (np.sqrt(2.0) * sigma))
            # Keep erfinv inside its domain when n lies a rounding error beyond the break.
            erf_value = np.clip(erf_min + np.minimum(n, n_break) / amplitude, -1.0, 1.0)
            m_lower = 10.0**(c_loc + np.sqrt(2.0) * sigma * erfinv(erf_value))
            m_upper = (self.CHABRIER_M_BREAK**p_upper
                       + p_upper * np.maximum(n - n_break, 0.0) / self.CHABRIER_NORM_UPPER)**(1.0 / p_upper)

        return np.clip(np.where(n <= n_break, m_lower, m_upper), self.M_MIN, self.M_MAX)

    def _number_imf_integral(self, m_l: FloatOrArray, m_r: FloatOrArray) -> FloatOrArray:
        """Calculate the unscaled mathematical integral of the number IMF.

        Parameters
        ----------
        m_l : float or numpy.ndarray
            Lower stellar mass bound in Msun.
        m_r : float or numpy.ndarray
            Upper stellar mass bound in Msun.

        Returns
        -------
        float or numpy.ndarray
            Integral of the number IMF in the interval; generally fractional.

        Notes
        -----
        The number IMF is integrated exactly as the difference of the
        closed-form cumulative numbers,
        and clamped to zero to prevent negative expectations.
        Bounds outside [M_MIN, M_MAX] are clipped by the helpers,
        which equals integrating the IMF as zero there.
        Inverted bounds raise ValueError.
        """

        # Reject inverted intervals.
        if not np.all(np.asarray(m_l) <= np.asarray(m_r)):
            raise ValueError("Lower mass bounds must not exceed upper bounds")

        # Number of stars per Msun of formed mass (not yet scaled by the particle mass).
        if self.imf_type == 'Kroupa':
            # or, num_stars_per_Msun = self._integral( lambda m: self.kroupa_number_imf(m), m_l, m_r )
            num_stars_per_Msun = (self._kroupa_closed_form_cumulative_number(m_r) -
                                  self._kroupa_closed_form_cumulative_number(m_l))

        elif self.imf_type == 'Chabrier':
            # or, num_stars_per_Msun = self._integral( lambda m: self.chabrier_number_imf(m), m_l, m_r )
            num_stars_per_Msun = (self._chabrier_closed_form_cumulative_number(m_r) -
                                  self._chabrier_closed_form_cumulative_number(m_l))

        else:
            raise ValueError(f"Unknown imf_type: {self.imf_type}")

        return np.maximum(0.0, num_stars_per_Msun)

    def _expected_star_count(self, m_l: FloatOrArray, m_r: FloatOrArray) -> FloatOrArray:
        """Scale the selected number IMF integral by the particle mass.

        Parameters
        ----------
        m_l : float or numpy.ndarray
            Lower stellar mass bound in Msun.
        m_r : float or numpy.ndarray
            Upper stellar mass bound in Msun.

        Returns
        -------
        float or numpy.ndarray
            Expected number of stars in the interval; generally fractional.

        Notes
        -----
        The integral of the number IMF is scaled by target_particle_mass
        and clamped to zero to prevent negative expectations.
        """

        num_stars_per_Msun = self._number_imf_integral(m_l, m_r)
        expected = np.maximum(0.0, num_stars_per_Msun * self.target_particle_mass)

        return float(expected) if np.ndim(expected) == 0 else expected

    def _expected_stellar_mass(self, m_l: FloatOrArray, m_r: FloatOrArray) -> FloatOrArray:
        """Integrate selected mass-weighted IMF over a mass interval.

        Parameters
        ----------
        m_l : float or numpy.ndarray
            Lower stellar mass bound in Msun.
        m_r : float or numpy.ndarray
            Upper stellar mass bound in Msun.

        Returns
        -------
        float or numpy.ndarray
            Expected summed initial stellar mass in the interval, in Msun.

        Notes
        -----
        Integral is in units of formed-mass fraction before multiplying by
        target_particle_mass.
        Clamped to zero to prevent negative mass expectations.
        """

        # Reject inverted intervals.
        if not np.all(np.asarray(m_l) <= np.asarray(m_r)):
            raise ValueError("Lower mass bounds must not exceed upper bounds")

        # Fraction of the formed mass in the requested interval.
        if self.imf_type == 'Kroupa':
            # or, mass_stars_per_Msun = self._integral(lambda m: self.kroupa_mass_weighted_imf(m), m_l, m_r,
            #                                          breaks=(self.KROUPA_M_BREAK,))
            mass_stars_per_Msun = (self._kroupa_closed_form_cumulative_mass(m_r) -
                                   self._kroupa_closed_form_cumulative_mass(m_l))

        elif self.imf_type == 'Chabrier':
            # or, mass_stars_per_Msun = self._integral(lambda m: self.chabrier_mass_weighted_imf(m), m_l, m_r,
            #                                          breaks=(self.CHABRIER_M_BREAK,))
            mass_stars_per_Msun = (self._chabrier_closed_form_cumulative_mass(m_r) -
                                   self._chabrier_closed_form_cumulative_mass(m_l))

        else:
            raise ValueError(f"Unknown imf_type: {self.imf_type}")
        expected = np.maximum(0.0, mass_stars_per_Msun * self.target_particle_mass)

        return float(expected) if np.ndim(expected) == 0 else expected

    def _single_power_law_alpha(self, m_l: float, m_r: float) -> Optional[float]:
        """Return the number-IMF slope if [m_l, m_r] lies in one power law.

        Parameters
        ----------
        m_l : float
            Lower stellar mass bound in Msun.
        m_r : float
            Upper stellar mass bound in Msun.

        Returns
        -------
        float or None
            The exponent alpha of xi(m) proportional to m**(-alpha) if the
            interval lies entirely above or below the selected IMF break
            (Kroupa: either side; Chabrier: above its break only), else
            None (interval straddles a break or touches the lognormal part).
        """

        if self.imf_type == 'Kroupa':
            if m_l >= self.KROUPA_M_BREAK:
                return self.KROUPA_ALPHA_UPPER
            if m_r <= self.KROUPA_M_BREAK:
                return self.KROUPA_ALPHA_LOWER
        elif self.imf_type == 'Chabrier':
            if m_l >= self.CHABRIER_M_BREAK:
                return self.CHABRIER_ALPHA_UPPER

        return None

    def _number_cdf(self, m: FloatOrArray, m_l: float, m_r: float, total_weight: Optional[float] = None) -> FloatOrArray:
        """Evaluate the conditional stellar *number* CDF on an interval.

        Parameters
        ----------
        m : float or numpy.ndarray
            Query mass or masses, in Msun.
        m_l : float
            Conditional interval's lower mass bound, in Msun.
        m_r : float
            Conditional interval's upper mass bound, in Msun.
        total_weight : float, optional
            Pre-computed expected star count on [m_l, m_r] to avoid recalculation.

        Returns
        -------
        float or numpy.ndarray
            Conditional fraction of stellar numbers at masses <= m.

        Notes
        -----
        Returns zero for a degenerate or inverted interval.
        Clipping query masses to [m_l, m_r] gives zero below the interval
        and unity above it, without integrating over reversed bounds.
        The conditional CDF requires a nonzero IMF number integral over
        the chosen interval.
        It is used for the CDF value at the yield break in
        _fixed_count_slice_yields.
        """

        if m_l >= m_r:
            return 0.0 if np.ndim(m) == 0 else np.zeros_like(m)

        # Clip the query masses so that the CDF is 0 below and 1 above the interval.
        m_clipped = np.clip(m, m_l, m_r)
        total = self._number_imf_integral(m_l, m_r) if total_weight is None else total_weight
        if not np.isfinite(total) or total <= 0:
            raise ValueError("Number CDF requires positive finite IMF weight")
        cdf = self._number_imf_integral(m_l, m_clipped) / total

        return cdf

    def _mass_cdf(self, m: FloatOrArray, m_l: float, m_r: float) -> FloatOrArray:
        """Evaluate the conditional cumulative *stellar mass* fraction.

        Parameters
        ----------
        m : float or numpy.ndarray
            Query mass or masses, in Msun.
        m_l : float
            Lower interval bound, in Msun.
        m_r : float
            Upper interval bound, in Msun.

        Returns
        -------
        float or numpy.ndarray
            Cumulative fraction of initial stellar mass, not a number CDF.

        Notes
        -----
        The denominator is the mass integral over [m_l, m_r].
        Queries are clipped before integrating, so out-of-range masses
        give zero or one.
        Degenerate or inverted intervals return zero.
        A zero-mass IMF interval has no defined conditional CDF.
        """

        if m_l >= m_r:
            return 0.0 if np.ndim(m) == 0 else np.zeros_like(m)

        m_clipped = np.clip(m, m_l, m_r)
        total = self._expected_stellar_mass(m_l, m_r)
        if not np.isfinite(total) or total <= 0:
            raise ValueError("Mass CDF requires positive finite IMF weight")
        cdf = self._expected_stellar_mass(m_l, m_clipped) / total

        return cdf

    def _inverse_number_cdf(self, quantiles: FloatOrArray, m_l: float, m_r: float) -> np.ndarray:
        """Invert the conditional number CDF (exact, closed form).

        Parameters
        ----------
        quantiles : array_like
            Scalar or array of finite quantiles in the closed interval [0, 1].
        m_l : float
            Lower stellar mass bound, in Msun.
        m_r : float
            Upper stellar mass bound, in Msun.

        Returns
        -------
        numpy.ndarray
            Inverted masses in Msun, retaining array shape; scalar inputs
            produce a one-element array.

        Notes
        -----
        If [m_l, m_r] lies inside a single power-law segment (Kroupa above or
        below its break, Chabrier above its break), the inverse is analytic,
        m = [m_l**p + u*(m_r**p - m_l**p)]**(1/p) with p = 1 - alpha.
        The result is clipped to [m_l, m_r], and the endpoints are reproduced to
        floating-point precision.
        This covers the SN progenitor range [M_SN_MIN, M_SN_MAX] for both IMFs.
        Otherwise (an interval straddling a break, or Chabrier's lognormal
        segment) the target cumulative number N(m_l) + u (N(m_r) - N(m_l)) is
        inverted exactly and piecewise by _inverse_closed_form_cumulative_number
        (powers and erfinv), and the endpoint masses are forced for u=0 and 1.
        The result is accurate to floating-point precision and needs no
        iteration.
        NumPy inputs are converted to finite floating-point arrays.
        Empty arrays return empty arrays without evaluating the CDF.
        """

        if not np.isfinite(m_l) or not np.isfinite(m_r) or m_r < m_l:
            raise ValueError("Inverse-CDF bounds must be finite and ordered")
        u = np.atleast_1d(np.asarray(quantiles, dtype=np.float64))
        if not np.all(np.isfinite(u) & (u >= 0.0) & (u <= 1.0)):
            raise ValueError("CDF quantiles must be finite and within [0, 1]")
        if u.size == 0:
            return np.empty_like(u)

        # Closed-form inverse for an interval inside one power-law segment.
        # The number CDF is proportional to m^p - m_l^p with p = 1 - alpha, so
        # m = [m_l^p + u (m_r^p - m_l^p)]^(1/p).
        alpha = self._single_power_law_alpha(m_l, m_r)
        if alpha is not None:
            p = 1.0 - alpha
            masses = (m_l**p + u * (m_r**p - m_l**p))**(1.0 / p)
            return np.clip(masses, m_l, m_r)

        if m_l == m_r:
            return np.full_like(u, m_l)

        # Otherwise (interval straddles a break or touches the log-normal
        # segment), invert the closed-form cumulative number exactly.
        if self.imf_type == 'Kroupa':
            cumulative_number = self._kroupa_closed_form_cumulative_number
        else:
            cumulative_number = self._chabrier_closed_form_cumulative_number
        n_l, n_r = float(cumulative_number(m_l)), float(cumulative_number(m_r))
        masses = np.clip(self._inverse_closed_form_cumulative_number(n_l + u * (n_r - n_l)), m_l, m_r)

        # Or, bisect on the number CDF.
        # Thirty bisections reduce the
        # bracket width to (m_r - m_l) / 2**30.
        # num_iterations = 30

        # # Lower and upper brackets of the inverted mass for each quantile [Msun].
        # m_low = np.full(u.shape, m_l, dtype=np.float64)
        # # m_high stores the upper bracket of each quantile-inverted stellar
        # # mass (Msun).
        # m_high = np.full(u.shape, m_r, dtype=np.float64)
        # total_weight = self._number_imf_integral(m_l, m_r)
        # for _ in range(num_iterations):
        #     m_mid = (m_low + m_high) * 0.5
        #     cdf_midpoint = self._number_cdf(m_mid, m_l, m_r, total_weight=total_weight)
        #     m_high = np.where((cdf_midpoint >= u), m_mid, m_high)
        #     m_low = np.where((cdf_midpoint < u), m_mid, m_low)

        # masses = 0.5 * (m_low + m_high)

        masses = np.where(u == 0.0, m_l, masses)
        masses = np.where(u == 1.0, m_r, masses)

        return masses

    def _sample_fixed_count_stars(self, m_l: float, m_r: float, star_count: int, is_stochastic: bool = False) -> np.ndarray:
        """Construct a fixed number of stellar masses within a mass interval.

        Parameters
        ----------
        m_l : float
            Lower stellar mass bound in Msun (scalar).
        m_r : float
            Upper stellar mass bound in Msun (scalar).
        star_count : int
            Strictly positive integer count of stars to generate.
        is_stochastic : bool
            Use random uniform quantiles if True.
            If False, use deterministic equal-number-slice
            mean masses (evenly spaced midpoint quantiles
            on intervals that are not a single power law).

        Returns
        -------
        Array of progenitor or ordinary stellar masses, in Msun.

        Raises
        ------
        ValueError
            If star_count is not positive, or if the bounds are not finite
            and ordered (raised by _inverse_number_cdf).

        Notes
        -----
        When is_stochastic=False the interval is cut into star_count slices
        of equal expected number.
        On a single power-law interval (always the case for the SN range)
        each slice is represented by its mean mass, in closed form,
        so the summed mass equals star_count times the IMF-mean mass exactly.
        Intervals that are not a single power law fall back to evenly spaced
        midpoint quantiles, which truncate the high-mass tail.
        These nodes are not random realizations.
        Because the yields are nonlinear in mass, evaluating
        them at these node masses overestimates gas and metals (10-25 percent
        for N of order one); the slice_averaged_yields option avoids this by
        using _fixed_count_slice_yields instead.
        """

        if star_count <= 0:
            raise ValueError("star_count must be positive")

        if is_stochastic:
            # Independent uniform quantiles for random mass draws.
            u = self.rng.random(star_count)
        else:
            # Evenly spaced midpoint quantiles.
            # They are used only if the
            # interval is not a single power law, where slice means have no
            # closed form.
            u = (np.arange(star_count) + 0.5) / star_count

            # On a single power-law interval, use the mean mass of each
            # equal-number slice.
            # The midpoint quantile would underestimate the
            # mean mass because it truncates the high-mass tail.
            alpha = self._single_power_law_alpha(m_l, m_r)
            if alpha is not None and m_l < m_r:
                edges = self._inverse_number_cdf(np.arange(star_count + 1) / star_count, m_l, m_r)
                slice_number = self._integral_power_law(alpha, edges[:-1], edges[1:])
                slice_mass = self._integral_power_law(alpha - 1.0, edges[:-1], edges[1:])
                return np.clip(np.divide(slice_mass, slice_number,
                                         out=np.full_like(slice_mass, m_l),
                                         where=slice_number != 0),
                               m_l, m_r)

        # Masses from the inverse number CDF.
        sampled_stellar_masses = self._inverse_number_cdf(u, m_l, m_r)

        return sampled_stellar_masses

    @staticmethod
    @lru_cache(maxsize=8)
    def _legendre_quadrature_rule(order: int) -> Tuple[np.ndarray, np.ndarray]:
        """Cache Gauss-Legendre nodes and weights for repeated fixed-count slices.

        Parameters
        ----------
        order : int
            Positive quadrature order.

        Notes
        -----
        Returned arrays are shared and must not be modified by callers.
        """

        if order < 1:
            raise ValueError("Quadrature order must be positive")

        return np.polynomial.legendre.leggauss(order)

    def _slice_quadrature_masses(self, m_l: float, m_r: float, star_count: int, quadrature_order: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Place Gauss-Legendre nodes inside each equal-number IMF slice.

        Parameters
        ----------
        m_l : float
            Lower progenitor mass bound in Msun (scalar).
        m_r : float
            Upper progenitor mass bound in Msun (scalar).
        star_count : int
            Strictly positive number of equal-number slices.
        quadrature_order : int
            Number of Gauss-Legendre nodes inside each slice.

        Returns
        -------
        masses : numpy.ndarray
            Node masses in Msun, shape (star_count, quadrature_order); row j
            holds the nodes of slice j in increasing quantile order.
        nodes, weights : numpy.ndarray
            Gauss-Legendre nodes and weights on [-1, 1].
            A slice mean of a function y is np.sum(y * 0.5 * weights, axis=1),
            because the number-CDF variable u is uniform within a slice.

        Raises
        ------
        ValueError
            If star_count is not positive.

        Notes
        -----
        Shared by _fixed_count_slice_yields and _fixed_count_slice_times,
        which differ only in the function averaged over the slice.
        """

        if star_count <= 0:
            raise ValueError("star_count must be positive")

        nodes, weights = self._legendre_quadrature_rule(quadrature_order)
        # u[j, i] is node i inside slice j, mapped from [-1, 1] to the slice.
        u = (np.arange(star_count)[:, None] + 0.5 * (nodes[None, :] + 1.0)) / star_count
        masses = self._inverse_number_cdf(u.ravel(), m_l, m_r).reshape(u.shape)

        return masses, nodes, weights

    def _fixed_count_slice_yields(self, m_l: float, m_r: float, star_count: int, quadrature_order: int = 16) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Average the SN yields over each equal-number IMF slice.

        Parameters
        ----------
        m_l : float
            Lower progenitor mass bound in Msun (scalar).
        m_r : float
            Upper progenitor mass bound in Msun (scalar).
        star_count : int
            Strictly positive number of equal-number slices.
        quadrature_order : int, optional
            Number of Gauss-Legendre nodes used inside each slice; the
            underlying rule is cached across fixed-count catalogs.

        Returns
        -------
        tuple of numpy.ndarray
            Slice-mean (energy [erg], gas [Msun], metals [Msun]), each of
            length star_count and in the same slice order as
            _sample_fixed_count_stars(..., is_stochastic=False).

        Notes
        -----
        Used only when slice_averaged_yields is True.
        A star drawn uniformly from one of the k equal-number slices has
        the expected yield <y> = k * integral y(m(u)) du over the slice, with
        u the conditional number-CDF variable. Summing these k values
        therefore reproduces k times the IMF-averaged yield, whereas summing
        y evaluated at the slice-mean masses is biased high for the concave
        gas yield and the convex-then-saturating metal yield (about +14 and
        +25 percent at one expected SN).
        The same holds for the explosion times: the stellar lifetime is
        convex in mass, so the lifetime at the slice-mean mass is earlier
        than the slice-mean lifetime.
        The integral is done with Gauss-Legendre quadrature in u.
        The one slice that straddles M_SN_YIELD_BREAK is split at the exact
        number-CDF breakpoint before quadrature, avoiding integration across
        the yield-law kink.
        """

        masses, nodes, weights = self._slice_quadrature_masses(m_l, m_r, star_count, quadrature_order)

        # Weights of the slice mean sum to one over the nodes.
        means = [np.sum(y * 0.5 * weights[None, :], axis=1)
                 for y in self._sn_yields_array(masses)]

        # Split only the slice crossing the yield-law kink at 40 Msun.
        # This removes Gauss-Legendre error from placing nodes on both branches.
        if m_l < self.M_SN_YIELD_BREAK < m_r:
            u_break = float(self._number_cdf(self.M_SN_YIELD_BREAK, m_l, m_r))
            crossing = min(int(u_break * star_count), star_count - 1)
            u_left, u_right = crossing / star_count, (crossing + 1) / star_count
            if u_left < u_break < u_right:
                corrected = np.zeros(3)
                for lower, upper in ((u_left, u_break), (u_break, u_right)):
                    sample_u = 0.5 * (lower + upper) + 0.5 * (upper - lower) * nodes
                    sample_m = self._inverse_number_cdf(sample_u, m_l, m_r)
                    for i, component in enumerate(self._sn_yields_array(sample_m)):
                        corrected[i] += star_count * (upper - lower) * 0.5 * np.dot(weights, component)
                for i in range(3):
                    means[i][crossing] = corrected[i]

        return tuple(means)

    def _fixed_count_slice_times(self, m_l: float, m_r: float, star_count: int, quadrature_order: int = 16) -> np.ndarray:
        """Average the explosion time over each equal-number IMF slice.

        Parameters
        ----------
        m_l : float
            Lower progenitor mass bound in Msun (scalar).
        m_r : float
            Upper progenitor mass bound in Msun (scalar).
        star_count : int
            Strictly positive number of equal-number slices.
        quadrature_order : int, optional
            Number of Gauss-Legendre nodes used inside each slice.

        Returns
        -------
        numpy.ndarray
            Slice-mean stellar lifetime [Myr] of length star_count, in the
            same slice order as _sample_fixed_count_stars(...,
            is_stochastic=False).

        Notes
        -----
        Used only when slice_averaged_yields is True.
        The lifetime is convex in mass, so its value at the slice-mean mass
        is earlier than the slice mean (the Poisson-weighted mean explosion
        time is 22 percent too early at one expected SN). Unlike the yields,
        the lifetime has no kink inside the SN range, so plain Gauss-Legendre
        quadrature in the number-CDF variable u is sufficient.
        """

        masses, _, weights = self._slice_quadrature_masses(m_l, m_r, star_count, quadrature_order)
        lifetimes = self._stellar_lifetime(masses, Z=self.stellar_metallicity)

        return np.sum(lifetimes * 0.5 * weights[None, :], axis=1)

    def _sample_poisson_stars(self, m_l: float, m_r: float) -> np.ndarray:
        """Draw a Poisson number of stars and randomly sample their masses.

        Parameters
        ----------
        m_l : float
            Lower stellar mass bound in Msun (scalar).
        m_r : float
            Upper stellar mass bound in Msun (scalar).

        Returns
        -------
        Possibly empty one-dimensional array of stellar masses (Msun).

        Notes
        -----
        Independent Poisson counts give correct expected mass, not exact
        particle-mass conservation.
        """

        # Poisson rate (expected star count) of this mass bin.
        expected_number_of_stars = self._expected_star_count(m_l, m_r)

        # Realized star count of this bin; it can be zero.
        sampled_bin_count = self.rng.poisson(expected_number_of_stars)

        if sampled_bin_count > 0:
            sampled_stellar_masses = self._sample_fixed_count_stars(m_l, m_r,
                                                                    sampled_bin_count,
                                                                    is_stochastic=True,
                                                                   )

        else:
            sampled_stellar_masses = np.array([])

        return sampled_stellar_masses

    def _sample_stars_in_bins(self, m_l: float, m_r: float, number_of_bins: int = 32) -> np.ndarray:
        """Sample a full interval by independent logarithmic-mass bins.

        Parameters
        ----------
        m_l : float
            Lower stellar mass bound in Msun (scalar, positive).
        m_r : float
            Upper stellar mass bound in Msun (scalar).
        number_of_bins : int
            Number of logarithmic mass bins for independent Poisson sampling.

        Returns
        -------
        Concatenated masses from all logarithmic bins, in bin order.

        Notes
        -----
        The bin counts are independent Poisson variates.
        Different bin choices are intended to describe the same underlying
        IMF statistically but change the exact seeded realization.
        """

        # Logarithmically spaced bin edges.
        mass_bin_edges = np.geomspace(m_l, m_r, number_of_bins + 1)
        # Enforce the exact end points despite floating-point rounding.
        mass_bin_edges[0], mass_bin_edges[-1] = m_l, m_r

        # Sample each mass bin independently.
        sampled_stellar_masses_list = [self._sample_poisson_stars(mass_bin_edges[i], mass_bin_edges[i + 1]) for i in range(number_of_bins)]

        return np.concatenate(sampled_stellar_masses_list)

    def _sample_population(self) -> None:
        """Build and summarize the full sampled stellar population.

        Returns
        -------
        None
            Sets sampled_star_masses, sampled_stellar_mass, and
            sampled_star_count.

        Raises
        ------
        MemoryError
            If the expected number of all stars exceeds MAX_SAMPLED_STARS.

        Notes
        -----
        The sampled total is a random variable; it should not be forced equal
        to target_particle_mass unless a mass-constrained model is intended.
        """

        # Check before drawing millions of samples.
        preliminary_star_count = self._expected_star_count(self.M_MIN, self.M_MAX)
        if preliminary_star_count > MAX_SAMPLED_STARS:
            raise MemoryError(f"Target particle mass implies {preliminary_star_count:.1e} sampled stars; "
                              f"the safe materialization limit is {MAX_SAMPLED_STARS:,.0f}. "
                              "Reduce --mass or use an aggregate population sampler.")

        # Independent bin counts make the total star number and mass random;
        # the mass is not constrained to equal target_particle_mass.
        self.sampled_star_masses = self._sample_stars_in_bins(self.M_MIN, self.M_MAX)
        self.sampled_stellar_mass = np.sum(self.sampled_star_masses)
        self.sampled_star_count = len(self.sampled_star_masses)

        return

    def _compute_expectations(self) -> None:
        """Store IMF-integrated means and separate AGORA comparison values.

        Notes
        -----
        expected_sn_energy, expected_ejected_gas_mass and
        expected_ejected_metal_mass are integrals of the mass-dependent
        _sn_yields(m) over the selected IMF; they are not fixed per-SN
        AGORA yields.
        The agora_expected_* attributes use their own 0.011 SNe/Msun event
        frequency and fixed ejecta per SN.
        Thus, these two sets of mean feedback values intentionally differ.
        The mass-dependent integrations rely on adaptive quadrature across
        the M_SN_YIELD_BREAK yield transition.
        expected_sn_mean_time is the number-IMF average of the stellar
        lifetime over the SN mass range; it lets the fixed-count catalogs
        be checked in time as well.
        """

        self.expected_star_count = self._expected_star_count(self.M_MIN, self.M_MAX)
        self.expected_stellar_mass = self._expected_stellar_mass(self.M_MIN, self.M_MAX)
        self.expected_sn_progenitor_mass = self._expected_stellar_mass(self.M_SN_MIN, self.M_SN_MAX)
        self.expected_sn_count, self.expected_sn_energy, self.expected_ejected_gas_mass, self.expected_ejected_metal_mass = self._integrate_sn_yields(self.M_SN_MIN, self.M_SN_MAX)
        self.agora_expected_sn_count, self.agora_expected_sn_energy, self.agora_expected_ejected_gas_mass, self.agora_expected_ejected_metal_mass = self._agora_expectations()
        self.expected_sn_mean_time = self._expected_sn_mean_time()

        return

    def _expected_sn_mean_time(self) -> float:
        """Average the explosion time of one SN over the number IMF.

        Returns
        -------
        float
            Integral of xi(m) * lifetime(m) over the SN range divided by the
            integral of xi(m), in Myr. The particle mass cancels.

        Notes
        -----
        The adaptive quadrature breaks at the selected IMF break mass, which
        lies outside the SN range and is harmless there.
        """

        if self.imf_type == 'Kroupa':
            number_imf, imf_breaks = self.kroupa_number_imf, (self.KROUPA_M_BREAK,)
        else:
            number_imf, imf_breaks = self.chabrier_number_imf, (self.CHABRIER_M_BREAK,)

        weighted_time = self._integral(lambda m: number_imf(m) * self._stellar_lifetime(m, Z=self.stellar_metallicity),
                                       self.M_SN_MIN, self.M_SN_MAX, breaks=imf_breaks)

        return weighted_time / self._number_imf_integral(self.M_SN_MIN, self.M_SN_MAX)

    def _stellar_lifetime(self, m: FloatOrArray, Z: float = 0.02) -> FloatOrArray:
        """Return the Raiteri-style stellar lifetime in Myr.

        Parameters
        ----------
        m : float or numpy.ndarray
            Stellar mass in Msun; must be positive.
        Z : float, optional
            Metallicity mass fraction for the fit (default 0.02). The event
            catalog passes self.stellar_metallicity explicitly.

        Returns
        -------
        float or numpy.ndarray
            Lifetime in Myr.

        Notes
        -----
        The metallicity-dependent coefficients a_0, a_1, a_2 define a
        quadratic in log10(m); its result is converted from years to Myr.
        Masses must be finite and strictly positive, and Z must lie in the
        closed interval [0.0004, 0.03]; this validation is explicit, so it
        remains enabled when Python runs with -O.
        The fit is not guaranteed monotonic for arbitrary Z and mass limits.
        At Z=0.02 it is indeed not monotonic: its minimum, 3.275 Myr, is at
        106 Msun, so a 120-Msun star (3.293 Myr) outlives a 106-Msun one.
        """

        if not np.all(np.isfinite(m) & (np.asarray(m) > 0.0)):
            raise ValueError("Stellar masses must be finite and positive")
        if not np.isfinite(Z) or not Z_MIN <= Z <= Z_MAX:
            raise ValueError(f"Metallicity must lie in [{Z_MIN:g}, {Z_MAX:g}]")

        log10_Z = np.log10(Z)

        # Metallicity-dependent coefficients of the quadratic in log10(m).
        a_0 = 10.13  + 0.07547 * log10_Z - 0.008084 * log10_Z**2
        a_1 = -4.424 - 0.7939  * log10_Z - 0.1187   * log10_Z**2
        a_2 =  1.262 + 0.3385  * log10_Z + 0.05417  * log10_Z**2

        log10_m = np.log10(m)
        lifetime_in_Myr = 10**(a_0 + a_1 * log10_m + a_2 * log10_m**2) / 1.0e6  # Convert yr to Myr.

        return lifetime_in_Myr

    def _sn_yields(self, m: float) -> Tuple[float, float, float]:
        """Calculate the adopted ejecta of one eligible SN progenitor.

        Parameters
        ----------
        m : float
            Scalar initial progenitor mass in Msun. Arrays are not supported.

        Returns
        -------
        tuple of float
            (ejected_energy [erg], ejected_gas_mass [Msun], ejected_metal_mass [Msun]).

        Notes
        -----
        Non-finite progenitor masses are rejected.
        Outside [M_SN_MIN, M_SN_MAX], all yields are zero.
        Each eligible progenitor releases 1e51 erg.
        At m <= M_SN_YIELD_BREAK, gas is m - COMPACT_REMNANT_MASS Msun
        and metals follow the oxygen/iron fitting formulas.
        Above M_SN_YIELD_BREAK, gas is held at the constant
        M_SN_YIELD_BREAK - COMPACT_REMNANT_MASS Msun and metals at 16.084 Msun,
        independent of m.
        These prescriptions are modeling assumptions, not numerical
        constraints.
        The gas yield is continuous at M_SN_YIELD_BREAK; the adopted
        constant metal yield approximates the fitted value there.
        Above the break, m - ejected_gas_mass is retained/untracked mass;
        do not interpret COMPACT_REMNANT_MASS as the final remnant there.
        This scalar wrapper delegates to _sn_yields_array, which holds the
        only copy of the yield formulas.
        """

        ejected_energy, ejected_gas_mass, ejected_metal_mass = self._sn_yields_array(m)

        return float(ejected_energy), float(ejected_gas_mass), float(ejected_metal_mass)

    def _sn_yields_array(self, m: FloatOrArray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Vectorized adopted ejecta of SN progenitors (the single yield model).

        Parameters
        ----------
        m : float or numpy.ndarray
            Initial progenitor mass(es) in Msun; every value must be finite.

        Returns
        -------
        tuple of numpy.ndarray
            (ejected_energy [erg], ejected_gas_mass [Msun], ejected_metal_mass [Msun]),
            each with the shape of m. Outside [M_SN_MIN, M_SN_MAX] all are zero.

        Notes
        -----
        _sn_yields wraps this method, so the yield prescription lives in
        exactly one place.
        Ineligible masses are replaced by M_SN_MIN inside the exponentials,
        which keeps every intermediate value finite, and are zeroed afterwards.
        The output is checked for nonnegative, mass-conserving gas and
        nonnegative metal yields no larger than the total gas ejecta.
        """

        m = np.asarray(m, dtype=np.float64)
        if not np.all(np.isfinite(m)):
            raise ValueError("Progenitor mass must be finite")

        eligible = (m >= self.M_SN_MIN) & (m <= self.M_SN_MAX)
        m_safe = np.where(eligible, m, self.M_SN_MIN)
        below_break = m_safe <= self.M_SN_YIELD_BREAK

        # Fixed explosion energy [erg] for every eligible progenitor.
        # All progenitors are assumed to explode successfully, including those above
        # M_SN_YIELD_BREAK; failed explosions, fallback-dependent energies, and
        # pair-instability SNe are ignored.
        ejected_energy = np.where(eligible, 1e51, 0.0)

        # Gas ejecta [Msun] after removing the remnant.
        # Above the yield break, they stay at their value at M_SN_YIELD_BREAK.
        ejected_gas_mass = np.where(eligible,
                                    np.minimum(m_safe, self.M_SN_YIELD_BREAK) - self.COMPACT_REMNANT_MASS,
                                    0.0)

        # Metal ejecta [Msun] from oxygen and iron fits (Kim et al. 2014,
        # Figure 2) to the solar-metallicity models of Woosley & Heger (2007).
        # The fits are used only for progenitors of 8-40 Msun. For 0.011 CCSNe
        # per Msun of formed stellar mass, they correspond to an average of
        # 0.0133 Msun of O and 0.0011 Msun of Fe per Msun of formed mass.
        ejected_O_mass = 27.66 * np.exp(-51.81 / m_safe)  # Equation (6) in Kim et al. 2014
        ejected_Fe_mass = 0.375 * np.exp(-17.94 / m_safe)  # Equation (5) in Kim et al. 2014
        # Total metals from the scaled O and Fe masses (Equation (4) in
        # Kim et al. 2014); above the yield break, a constant is adopted.
        ejected_metal_mass = np.where(eligible,
                                      np.where(below_break, 2.09 * ejected_O_mass + 1.06 * ejected_Fe_mass, 16.084),
                                      0.0)

        # Catch unphysical edits to the empirical yield prescription early.
        if np.any(eligible & ((ejected_gas_mass < 0)
                              | (ejected_gas_mass > m)
                              | (ejected_metal_mass < 0)
                              | (ejected_metal_mass > ejected_gas_mass))):
            raise ValueError("SN yield model violates ejecta mass conservation")

        return ejected_energy, ejected_gas_mass, ejected_metal_mass

    def _agora_expectations(self) -> Tuple[float, float, float, float]:
        """Compute independent AGORA-style benchmark feedback expectations.

        Returns
        -------
        tuple of float
            Expected SN count, explosion energy (erg), gas ejecta (Msun),
            and metal ejecta (Msun).

        Notes
        -----
        Uses 0.011 SNe/Msun formed, 1e51 erg/SN, 14.8 Msun gas/SN,
        and 2.63 Msun metals/SN.
        These values are not calculated from the selected IMF or the
        event-level _sn_yields() prescription.
        The variable names agora_* refer only to this comparison model.
        """

        agora_sn_rate_per_msun = 0.011
        agora_sn_count = self.target_particle_mass * agora_sn_rate_per_msun

        # Apply the same fixed energy to every benchmark event.
        agora_energy_per_sn = 1e51
        agora_sn_energy = agora_sn_count * agora_energy_per_sn

        # Gas ejected per SN, in Msun (Kim et al. 2016, Section 3.2).
        agora_gas_mass_per_sn = 14.8
        agora_gas_mass = agora_sn_count * agora_gas_mass_per_sn

        # Metals ejected per SN, in Msun.
        # This equals (2.09*0.0133 + 1.06*0.0011) / 0.011;
        # see Section 3.2 of Kim et al. (2016)
        # and Equation (4) of Kim et al. (2014).
        agora_metal_mass_per_sn = 2.63
        agora_metal_mass = agora_sn_count * agora_metal_mass_per_sn

        return agora_sn_count, agora_sn_energy, agora_gas_mass, agora_metal_mass

    def _integrate_sn_yields(self, m_l: float, m_r: float) -> Tuple[float, float, float, float]:
        """Integrate event-level SN feedback over a chosen mass interval.

        Parameters
        ----------
        m_l : float
            Requested lower stellar mass limit (Msun).
        m_r : float
            Requested upper stellar mass limit (Msun).

        Returns
        -------
        tuple of float
            (expected number of SNe, energy [erg], gas [Msun], metals [Msun]).

        Notes
        -----
        Clips the interval to the configured CCSN progenitor range and
        integrates xi(m)*_sn_yields(m)[component] independently for each
        quantity.
        All integrals are multiplied by target_particle_mass.
        The integral is split at M_SN_YIELD_BREAK when the threshold lies
        inside the requested interval; the quadrature helper also considers
        the selected IMF break.
        Returns zeros for an empty intersection.
        Non-finite mass limits or inverted input intervals are rejected.
        """

        if not np.isfinite(m_l) or not np.isfinite(m_r) or m_l > m_r:
            raise ValueError("SN mass limits must be finite and ordered")

        # Clamp the requested interval to SN progenitor limits.
        m_upper = min(m_r, self.M_SN_MAX)
        m_lower = max(m_l, self.M_SN_MIN)

        if m_lower >= m_upper:
            return (0.0,) * 4

        if self.imf_type == 'Kroupa':
            number_imf = self.kroupa_number_imf
            imf_breaks = (self.KROUPA_M_BREAK,)
        elif self.imf_type == 'Chabrier':
            number_imf = self.chabrier_number_imf
            imf_breaks = (self.CHABRIER_M_BREAK,)
        else:
            raise ValueError(f"Unknown imf_type: {self.imf_type}")

        # Expected (non-integer) number of CCSNe.
        expected_number_of_sne = self._expected_star_count(m_lower, m_upper)

        # Split at the yield-prescription transition so each interval uses
        # a single branch of the mass-dependent model.
        bounds = ([m_lower, self.M_SN_YIELD_BREAK, m_upper]
                  if m_lower < self.M_SN_YIELD_BREAK < m_upper else [m_lower, m_upper])

        def integrate_component(index: int) -> float:
            """Integrate xi(m) * _sn_yields(m)[index] over the SN range.

            Parameters
            ----------
            index : int
                0 for the energy (erg), 1 for gas (Msun), 2 for metals
                (Msun).

            Returns
            -------
            float
                The expected total of that component for the particle,
                summed over the sub-intervals in bounds.
            """

            return self.target_particle_mass * sum(self._integral(lambda m: number_imf(m) * self._sn_yields(m)[index], a, b, breaks=imf_breaks) for a, b in zip(bounds[:-1], bounds[1:]))

        expected_ejected_energy = integrate_component(0)
        expected_ejected_gas_mass = integrate_component(1)
        expected_ejected_metal_mass = integrate_component(2)

        return expected_number_of_sne, expected_ejected_energy, expected_ejected_gas_mass, expected_ejected_metal_mass

    def _append_sn_events(self, sn_events: List[Dict[str, Any]], progenitor_masses: Union[Sequence[float], np.ndarray],
                          yields: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None,
                          explosion_times: Optional[np.ndarray] = None) -> None:
        """Append one event dictionary per progenitor, computed in one batch.

        Parameters
        ----------
        sn_events : list of dict
            Mutable event list to append to.
        progenitor_masses : array_like
            Initial stellar masses of the exploding progenitors in Msun.
        yields : tuple of numpy.ndarray, optional
            (energy [erg], gas [Msun], metals [Msun]) arrays, one entry per
            progenitor, to store instead of _sn_yields_array(progenitor_masses).
            The fixed-count catalogs pass the slice-averaged yields here when
            slice_averaged_yields is True (see _fixed_count_slice_yields).
            All three arrays must match progenitor_masses in shape and contain
            finite, nonnegative values; otherwise ValueError is raised.
        explosion_times : numpy.ndarray, optional
            Explosion times [Myr], one per progenitor, to store instead of the
            lifetimes of progenitor_masses.
            The fixed-count catalogs pass the slice-averaged lifetimes here
            when slice_averaged_yields is True (see _fixed_count_slice_times).
            Same shape and finiteness rules as yields,
            and the times must be positive.

        Returns
        -------
        None
            Appends the events to sn_events; an empty input appends nothing.

        Notes
        -----
        Each dictionary stores progenitor_mass (Msun), explosion_time (Myr),
        ejected_energy (erg), ejected_gas_mass and ejected_metal_mass (Msun).
        Assumes coeval star formation; explosion time equals lifetime.
        The lifetime and yield models are evaluated on whole arrays.
        The list is not re-sorted here; callers sort it by explosion_time once,
        after all events have been appended.
        Shape checks prevent zip from silently discarding events
        when caller-supplied yield arrays have fewer elements than masses.
        """

        masses = np.asarray(progenitor_masses, dtype=np.float64)
        if masses.size == 0:
            return

        # The population is coeval, so the explosion time is the stellar lifetime.
        if explosion_times is None:
            explosion_times = self._stellar_lifetime(masses, Z=self.stellar_metallicity)
        else:
            explosion_times = np.asarray(explosion_times, dtype=np.float64)
            if explosion_times.shape != masses.shape or not np.all(np.isfinite(explosion_times) & (explosion_times > 0)):
                raise ValueError("explosion_times must be positive, finite, and match the progenitor-mass shape")
        if yields is None:
            ejected_energy, ejected_gas_mass, ejected_metal_mass = self._sn_yields_array(masses)
        else:
            if len(yields) != 3:
                raise ValueError("yields must contain energy, gas, and metal arrays")
            components = tuple(np.asarray(component, dtype=np.float64) for component in yields)
            if any(component.shape != masses.shape for component in components):
                raise ValueError("Each yield array must match the progenitor-mass array shape")
            if any(not np.all(np.isfinite(component) & (component >= 0)) for component in components):
                raise ValueError("Yield arrays must contain finite, nonnegative values")
            ejected_energy, ejected_gas_mass, ejected_metal_mass = components

        sn_events.extend({'progenitor_mass': m,
                          'explosion_time': t,
                          'ejected_energy': e,
                          'ejected_gas_mass': g,
                          'ejected_metal_mass': z}
                         for m, t, e, g, z in zip(masses.tolist(), explosion_times.tolist(),
                                                  np.asarray(ejected_energy).tolist(),
                                                  np.asarray(ejected_gas_mass).tolist(),
                                                  np.asarray(ejected_metal_mass).tolist()))

        return

    def _build_sampled_sn_events(self) -> None:
        """Collect and summarize sampled progenitors eligible for CCSNe.

        Notes
        -----
        Populates sampled_sn_events and the sampled_sn_count,
        sampled_first_sn_time, sampled_last_sn_time,
        sampled_sn_progenitor_mass, sampled_sn_energy,
        sampled_ejected_gas_mass and sampled_ejected_metal_mass attributes.
        The first and last explosion times are NaN if there are no SNe;
        all summed counts, progenitor masses, and yields then equal zero.
        Events are sorted by explosion time once, after all have been
        collected.
        """

        self.sampled_sn_events = []
        # Select the stars that are eligible CCSN progenitors.
        valid_mask = (self.sampled_star_masses >= self.M_SN_MIN) & (self.sampled_star_masses <= self.M_SN_MAX)
        sn_progenitor_masses = self.sampled_star_masses[valid_mask]
        self._append_sn_events(self.sampled_sn_events, sn_progenitor_masses)
        self.sampled_sn_events.sort(key=lambda x: x['explosion_time'])

        self.sampled_sn_count = len(self.sampled_sn_events)
        has_sampled = self.sampled_sn_count > 0
        self.sampled_first_sn_time = self.sampled_sn_events[0]['explosion_time'] if has_sampled else np.nan
        self.sampled_last_sn_time = self.sampled_sn_events[-1]['explosion_time'] if has_sampled else np.nan

        self.sampled_sn_progenitor_mass = sum(sn['progenitor_mass'] for sn in self.sampled_sn_events)
        self.sampled_sn_energy = sum(sn['ejected_energy'] for sn in self.sampled_sn_events)
        self.sampled_ejected_gas_mass = sum(sn['ejected_gas_mass'] for sn in self.sampled_sn_events)
        self.sampled_ejected_metal_mass = sum(sn['ejected_metal_mass'] for sn in self.sampled_sn_events)

        return

    def _build_fixed_count_sn_events(self) -> None:
        """Build deterministic catalogs and Poisson-weighted feedback averages.

        Notes
        -----
        The tabulated maximum is max(2, ceil(lambda+5*sqrt(lambda))),
        or 0 if the grid is skipped (see below),
        with lambda=expected_sn_count; the floor of 2 keeps at least
        two nonempty catalogs for very small populations.
        Catalog k has k fixed slice-mean masses of the conditional
        SN-number IMF, which set the explosion times and progenitor masses.
        Their gas and metal ejecta are _sn_yields at those masses by default,
        which is biased high for small populations; with
        slice_averaged_yields=True they are the slice-averaged values
        from _fixed_count_slice_times and _fixed_count_slice_yields.
        Catalog zero is empty.
        P(0) is retained explicitly.
        For k < max, weights equal the Poisson PMF; the final weight uses
        scipy.stats.poisson.sf to retain precision when the upper tail is
        very small.
        The last catalog represents the entire remaining Poisson tail with
        only max events.
        Its tabulated averages are thus approximate, and
        finite-count quantile integration also creates quadrature error.
        The five weighted totals track event count, progenitor mass, energy,
        gas ejecta, and metal ejecta.
        Afterwards each Poisson-weighted total (and the SN count and the mean
        explosion time per SN) is compared with the IMF-integrated expectation,
        and a warning is issued if they differ by more than 1 percent (usually
        due to upper Poisson tail truncation limits and slice-quadrature error).
        A warning is also issued if the Poisson weights do not sum to one.
        Because catalog k stores k events, the grid grows quadratically.
        When lambda exceeds 1000, a warning is issued and max_sn_count is 0.
        No nonempty catalog is built: fixed_count_sn_events is
        [[]], the per-count arrays have the single entry zero, the weight
        array holds only P(N=0), and all poisson_weighted_* totals stay zero.
        Every later loop over catalogs 1..max_sn_count is then empty, so
        print_summary, write_feedback_tables, plot_imf and
        plot_sn_feedback_rates still run and produce empty fixed-count
        output.
        """

        # Cover the Poisson distribution up to 5 standard deviations above the
        # mean; keep at least two nonempty catalogs for small expected counts.
        poisson_sigma_cap = 5.0
        # Above this expected count, skip the (quadratically growing) catalogs.
        max_expected_sn_count = 1000.0
        if self.expected_sn_count > max_expected_sn_count:
            warnings.warn(f"expected_sn_count={self.expected_sn_count:.0f} > {max_expected_sn_count:.0f}: "
                           "fixed-count catalogs are impractical; use a smaller particle mass "
                           "or the randomly sampled catalog only")
            self.max_sn_count = 0
        else:
            self.max_sn_count = max(2,
                                    int(np.ceil(self.expected_sn_count +
                                                poisson_sigma_cap * np.sqrt(self.expected_sn_count))))

        # Catalogs and summary arrays are indexed by the event count k.
        self.fixed_count_sn_events = [[] for _ in range(self.max_sn_count + 1)]
        self.fixed_count_sn_counts = np.zeros(self.max_sn_count + 1, dtype=int)
        self.fixed_count_first_sn_times = np.full(self.max_sn_count + 1, np.nan)
        self.fixed_count_last_sn_times = np.full(self.max_sn_count + 1, np.nan)
        self.fixed_count_progenitor_masses = np.zeros(self.max_sn_count + 1)
        self.fixed_count_energies = np.zeros(self.max_sn_count + 1)
        self.fixed_count_gas_masses = np.zeros(self.max_sn_count + 1)
        self.fixed_count_metal_masses = np.zeros(self.max_sn_count + 1)
        self.fixed_count_poisson_weights = np.zeros(self.max_sn_count + 1)

        if self.max_sn_count > 0:
            # The last catalog represents N >= max_sn_count.
            # The survival function keeps the small upper-tail probability accurate.
            self.fixed_count_poisson_weights[:-1] = poisson.pmf(np.arange(self.max_sn_count), self.expected_sn_count)
            self.fixed_count_poisson_weights[-1] = poisson.sf(self.max_sn_count - 1, self.expected_sn_count)
        else:
            self.fixed_count_poisson_weights[0] = poisson.pmf(0, self.expected_sn_count)

        self.poisson_weighted_sn_count = 0.0
        self.poisson_weighted_progenitor_mass = 0.0
        self.poisson_weighted_sn_energy = 0.0
        self.poisson_weighted_ejected_gas_mass = 0.0
        self.poisson_weighted_ejected_metal_mass = 0.0
        self.poisson_weighted_mean_sn_time = 0.0
        poisson_weighted_time_sum = 0.0

        for sn_count in range(1, self.max_sn_count + 1):
            # Each equal-number IMF slice is represented by its mean mass,
            # and optionally by its mean yields.
            sn_progenitor_masses = self._sample_fixed_count_stars(self.M_SN_MIN, self.M_SN_MAX, sn_count, is_stochastic=False)
            slice_yields = (self._fixed_count_slice_yields(self.M_SN_MIN, self.M_SN_MAX, sn_count)
                            if self.slice_averaged_yields else None)
            slice_times = (self._fixed_count_slice_times(self.M_SN_MIN, self.M_SN_MAX, sn_count)
                           if self.slice_averaged_yields else None)
            self._append_sn_events(self.fixed_count_sn_events[sn_count], sn_progenitor_masses,
                                   yields=slice_yields, explosion_times=slice_times)
            # Sort chronologically; the lifetime is not monotonic in mass.
            self.fixed_count_sn_events[sn_count].sort(key=lambda x: x['explosion_time'])

            self.fixed_count_sn_counts[sn_count] = len(self.fixed_count_sn_events[sn_count])
            self.fixed_count_first_sn_times[sn_count] = self.fixed_count_sn_events[sn_count][0]['explosion_time']
            self.fixed_count_last_sn_times[sn_count] = self.fixed_count_sn_events[sn_count][-1]['explosion_time']

            for sn in self.fixed_count_sn_events[sn_count]:
                self.fixed_count_progenitor_masses[sn_count] += sn['progenitor_mass']
                self.fixed_count_energies[sn_count] += sn['ejected_energy']
                self.fixed_count_gas_masses[sn_count] += sn['ejected_gas_mass']
                self.fixed_count_metal_masses[sn_count] += sn['ejected_metal_mass']

            probability = self.fixed_count_poisson_weights[sn_count]
            self.poisson_weighted_sn_count += probability * self.fixed_count_sn_counts[sn_count]
            self.poisson_weighted_progenitor_mass += probability * self.fixed_count_progenitor_masses[sn_count]
            self.poisson_weighted_sn_energy += probability * self.fixed_count_energies[sn_count]
            self.poisson_weighted_ejected_gas_mass += probability * self.fixed_count_gas_masses[sn_count]
            self.poisson_weighted_ejected_metal_mass += probability * self.fixed_count_metal_masses[sn_count]
            poisson_weighted_time_sum += probability * sum(sn['explosion_time'] for sn in self.fixed_count_sn_events[sn_count])

        if self.poisson_weighted_sn_count > 0:
            self.poisson_weighted_mean_sn_time = poisson_weighted_time_sum / self.poisson_weighted_sn_count

        # Compare the Poisson-weighted totals with the IMF-integrated
        # expectations, skipping grids that were not built.
        for label, expected, weighted in (
                                             ("SN count",            self.expected_sn_count,           self.poisson_weighted_sn_count),
                                             ("progenitor mass",     self.expected_sn_progenitor_mass, self.poisson_weighted_progenitor_mass),
                                             ("energy",              self.expected_sn_energy,          self.poisson_weighted_sn_energy),
                                             ("gas mass",            self.expected_ejected_gas_mass,   self.poisson_weighted_ejected_gas_mass),
                                             ("metal mass",          self.expected_ejected_metal_mass, self.poisson_weighted_ejected_metal_mass),
                                             ("mean explosion time", self.expected_sn_mean_time,       self.poisson_weighted_mean_sn_time),
                                         ):
            if self.max_sn_count > 0 and expected > 0 and abs(weighted / expected - 1.0) > 0.01:
                hint = (" (the bias of mean-mass nodes; try slice_averaged_yields=True / --slice-averaged-yields)"
                        if label in ("gas mass", "metal mass", "mean explosion time") and not self.slice_averaged_yields else "")
                warnings.warn(f"Poisson-weighted {label} differs from the IMF-integrated "
                              f"expectation by {100*(weighted/expected - 1.0):+.1f}%{hint}")

        # The weights P(0), ..., P(K-1), P(N >= K) must sum to one.
        if self.max_sn_count > 0 and not np.isclose(self.fixed_count_poisson_weights.sum(), 1.0, rtol=0.0, atol=1.0e-9):
            warnings.warn(f"Poisson weights sum to {self.fixed_count_poisson_weights.sum():.12f}, not 1")

        return

    def plot_imf(self, output_suffix: str = '', output_dir: Optional[str] = None) -> None:
        """Plot analytic, randomly sampled, and deterministic SN mass spectra.

        Parameters
        ----------
        output_suffix : str, optional
            Suffix incorporated into the output PNG filename.
        output_dir : str or None, optional
            Directory receiving the PNG; created if missing. The default is
            the current working directory.

        Notes
        -----
        The horizontal axis is log10(m/Msun).
        The analytic number IMF xi(m)=dN/dm is converted
        into dN/dlog10(m) using m*ln(10).
        The sampled histogram uses density=False; raw counts are divided by
        log10-mass bin widths, preserving the realized star count.
        An empty sampled catalog produces zero bin counts.
        Every fixed-count histogram is normalized and then multiplied by its
        SN count and Poisson weight.
        A separate analytic CCSN expectation makes the finite-node
        histogram approximation visible.
        The resulting curves have units of number of stars per
        log10-mass interval.
        Logarithmic y-limits can hide bins containing no stars or small
        positive densities.
        If the fixed-count catalogs were skipped (expected_sn_count above
        1000), the averaged-SN histogram is identically zero and does not
        appear on the logarithmic axis.
        """

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

        imf_plot_points = 256
        m = np.geomspace(self.M_MIN, self.M_MAX, imf_plot_points)
        m[0], m[-1] = self.M_MIN, self.M_MAX
        # Convert the IMF from xi(m)=dN/dm into dN/dlog10(m).
        jacobian = m * np.log(10)
        if self.imf_type == 'Kroupa':
            ax.plot(np.log10(m),
                    self.target_particle_mass * self.chabrier_number_imf(m) * jacobian,
                    color='darkorange',
                    ls='--',
                    label='Chabrier',
                    zorder=1,
                   )
            ax.plot(np.log10(m),
                    self.target_particle_mass * self.kroupa_number_imf(m) * jacobian,
                    color='b',
                    label=f"Kroupa, $N_\\star$={int(round(self.expected_star_count))}",
                    zorder=1,
                   )

        elif self.imf_type == 'Chabrier':
            ax.plot(np.log10(m),
                    self.target_particle_mass * self.kroupa_number_imf(m) * jacobian,
                    color='b',
                    ls='--',
                    label='Kroupa',
                    zorder=1,
                   )
            ax.plot(np.log10(m),
                    self.target_particle_mass * self.chabrier_number_imf(m) * jacobian,
                    color='darkorange',
                    label=f"Chabrier, $N_\\star$={int(round(self.expected_star_count))}",
                    zorder=1,
                   )

        number_of_bins = 32
        sampled_star_density, sampled_stars_bin_edges = np.histogram(np.log10(self.sampled_star_masses),
                                                                     bins=number_of_bins,
                                                                     range=(np.log10(self.M_MIN), np.log10(self.M_MAX)),
                                                                     density=False,
                                                                    )
        # Convert raw counts to counts per dex before plotting.
        sampled_star_density = sampled_star_density / np.diff(sampled_stars_bin_edges)
        ax.stairs(np.maximum(sampled_star_density, LOG_PLOT_FLOOR), sampled_stars_bin_edges,
                  color='r',
                  label=f"Rand. Sampled, $N_\\star$={self.sampled_star_count}, $N_{{\\rm SNe}}$={self.sampled_sn_count}",
                 )

        # Sum the fixed-count progenitor distributions, weighted by the event
        # count and the Poisson probability. This stays zero if the grid was skipped.
        averaged_sne_bin_edges = np.linspace(np.log10(self.M_MIN), np.log10(self.M_MAX), number_of_bins + 1)
        averaged_sne_histogram = np.zeros(number_of_bins)
        for sn_count in range(1, self.max_sn_count + 1):
            # Fixed bins are already stored in averaged_sne_bin_edges;
            # the per-catalog histogram edges are not needed.
            sne_histogram, _ = np.histogram(np.log10([sn['progenitor_mass'] for sn in self.fixed_count_sn_events[sn_count]]),
                                            bins=number_of_bins,
                                            range=(np.log10(self.M_MIN), np.log10(self.M_MAX)),
                                            density=True,
                                           )
            averaged_sne_histogram += sne_histogram * self.fixed_count_sn_counts[sn_count] * self.fixed_count_poisson_weights[sn_count]

        if self.max_sn_count > 0:
            ax.stairs(np.maximum(averaged_sne_histogram, LOG_PLOT_FLOOR), averaged_sne_bin_edges,
                      color='g',
                      label=f"Det. Averaged, $N_{{\\rm SNe}}$={self.poisson_weighted_sn_count:.2f}",
                     )

        # Integrate the conditional SN number IMF exactly in each displayed
        # log-mass bin.
        # This reference is independent of representative nodes.
        bin_cdf = self._number_cdf(10.0**averaged_sne_bin_edges,
                                   self.M_SN_MIN, self.M_SN_MAX)
        analytic_sne_histogram = (self.expected_sn_count * np.diff(bin_cdf)
                                  / np.diff(averaged_sne_bin_edges))
        ax.stairs(np.maximum(analytic_sne_histogram, LOG_PLOT_FLOOR),
                  averaged_sne_bin_edges, color='purple', linestyle=':',
                  label='Analytic expected CCSN IMF')

        ax.set_xlim(np.log10(self.M_MIN)-0.1, np.log10(self.M_MAX)+0.2)
        ax.set_ylim(5.0e-4*self.target_particle_mass, 5.0e+0*self.target_particle_mass)

        ax.set_xlabel(r'$\log\left(\frac{m}{\rm M_\odot}\right)$')
        ax.set_ylabel(r'$\frac{{\rm d}N}{{\rm d\,}\log\left(\frac{m}{\rm M_\odot}\right)}$')

        ax.legend(loc='upper right')
        ax.set_yscale('log')

        ax.minorticks_on()
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.tick_params(which='both', direction='in')

        def log10_to_mass(x: FloatOrArray) -> FloatOrArray:
            """Convert base-10 logarithmic stellar mass to linear mass.

            Parameters
            ----------
            x : float or numpy.ndarray
                Input value or NumPy array for the logarithmic axis conversion.

            Returns
            -------
            Mass equal to 10**x.

            Notes
            -----
            Used solely for the secondary plotting axis.
            """

            return 10**x

        def mass_to_log10(x: FloatOrArray) -> FloatOrArray:
            """Convert positive linear stellar mass to base-10 logarithmic
            mass.

            Parameters
            ----------
            x : float or numpy.ndarray
                Input value or NumPy array for the logarithmic axis conversion.

            Returns
            -------
            log10(x); zero/negative inputs can yield -inf or nan.

            Notes
            -----
            Used solely for the secondary plotting axis.
            """

            return np.log10(np.clip(x, 1e-30, None))

        # Secondary x-axis on top, showing the linear mass.
        secax = ax.secondary_xaxis('top', functions=(log10_to_mass, mass_to_log10))
        ticks = [self.M_MIN, self.KROUPA_M_BREAK, self.CHABRIER_M_BREAK, self.M_SN_MIN, self.M_SN_YIELD_BREAK, 80.0, self.M_MAX]
        secax.set_xticks(np.array(ticks))
        secax.set_xticklabels([f"{t:g}" for t in ticks])
        secax.set_xlabel(r'$m\quad {\rm[M_\odot]}$')

        # Save the figure to the output directory.
        plt.tight_layout()
        if output_dir is None:
            output_dir = '.'
        os.makedirs(output_dir, exist_ok=True)
        fig.savefig(os.path.join(output_dir, f'fig_IMF__{output_suffix}.png'))
        plt.close(fig)

        return

    def plot_sn_feedback_rates(self, output_suffix: str = '', output_dir: Optional[str] = None) -> None:
        """Plot sampled and Poisson-averaged SN feedback *rates* over time.

        Parameters
        ----------
        output_suffix : str, optional
            Filename suffix for the three-panel PNG.
        output_dir : str or None, optional
            Directory receiving the PNG; created if missing. The default is
            the current working directory.

        Notes
        -----
        Histogram counts are divided by the bin width to obtain energy
        (erg/Myr), gas ejecta (Msun/Myr), and metal ejecta (Msun/Myr) rates.
        Integrating each rate over all plotted bins recovers the
        corresponding catalog total if all events fall within the plotted
        time interval.
        Every stored explosion time is verified against the plot limits,
        which properly accounts for non-monotonic lifetime relations.
        Curves for k >= 1 are combined using their stored Poisson weights.
        Logarithmic y-limits adapt to the positive plotted rates,
        spanning at most three decades below the peak with headroom above it
        for the legend, and use a fallback range for all-zero feedback.
        Time bins are one Myr wide.
        The upper time bound is at least 45 Myr and is extended to the next
        integer Myr when any stored event explodes later.
        If the fixed-count catalogs were skipped (expected_sn_count above
        1000), the Poisson-averaged curves are identically zero; the y-limits
        then come from the sampled curves alone (or the fallback range).
        """

        # Energy, gas, and metal injection rates are shown in separate panels.
        fig, (ax_energy, ax_gas, ax_metal) = plt.subplots(1, 3, figsize=(16, 6))

        # Time range [Myr]: coeval formation begins at t=0.
        time_min = 0.0

        # Include the lifetime of the least massive eligible progenitor.
        last_lifetime = self._stellar_lifetime(self.M_SN_MIN, Z=self.stellar_metallicity)
        time_max = max(45.0, np.ceil(last_lifetime))

        # Check all stored events, since the lifetime is not monotonic in mass.
        for catalog in [self.sampled_sn_events, *self.fixed_count_sn_events]:
            for sn in catalog:
                t = sn['explosion_time']
                if np.isfinite(t) and t > time_max:
                    time_max = np.ceil(t)  # Extend the window to include this event.
                elif not np.isfinite(t) or t < time_min:
                    warnings.warn(f"Explosion at {t} Myr is outside the plot window")


        # An integer bin count gives exactly 1-Myr-wide bins.
        time_bin_count = int(time_max - time_min)

        def format_rate_axes(ax: plt.Axes, sampled_rate: np.ndarray, averaged_rate: np.ndarray) -> None:
            """Format an SN feedback-rate axis with logarithmic y-scaling.

            Parameters
            ----------
            ax : matplotlib.axes.Axes
                Axis receiving the time range, legend, ticks, and log y-scale.
            sampled_rate : numpy.ndarray
                Rates of the random catalog, used only for the y limits.
            averaged_rate : numpy.ndarray
                Poisson-weighted rates, used only for the y limits.

            Notes
            -----
            This helper determines logarithmic y-limits from finite positive
            rates: the upper limit is 10 times the peak (legend headroom) and
            the lower limit is the larger of half the smallest rate and
            1e-3 of the peak, so the plotted range never exceeds about four
            decades.
            Empty or zero-only distributions use fallback limits.
            """

            ax.set_xlim(time_min, time_max)
            ax.set_xlabel(r'$t\quad {\rm [Myr]}$')
            ax.legend(loc='upper right')
            ax.set_yscale('log')

            # Set the log-scale y limits from the finite positive rates.
            rates = np.concatenate((sampled_rate, averaged_rate))
            positive = rates[np.isfinite(rates) & (rates > 0)]
            if positive.size:
                max_dynamic_range = 1.0e3  # ratio peak/floor: negligible Poisson-tail bins
                legend_headroom = 10.0  # keeps the legend clear of the tallest bars
                ax.set_ylim(max(positive.min() / 2, positive.max() / max_dynamic_range),
                            positive.max() * legend_headroom)
            else:
                ax.set_ylim(1.0e-3, 1.0)

            ax.minorticks_on()
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.tick_params(which='both', direction='in')

            return

        def feedback_rate_histogram(sn_events: List[Dict[str, Any]], ejected_quantity: str) -> Tuple[np.ndarray, np.ndarray]:
            """Histogram one feedback quantity and convert totals to rates.

            Parameters
            ----------
            sn_events : list of dict
                Sorted or unsorted SN events with explosion_time and yield
                keys.
            ejected_quantity : {'ejected_energy', 'ejected_gas_mass',
            'ejected_metal_mass'}
                Dictionary key identifying the weight contributed by each
                event.

            Returns
            -------
            feedback_rate : numpy.ndarray
                Weighted quantity per Myr in each time bin (units follow
                weights).
            bin_edges : numpy.ndarray
                Shared time-bin boundaries in Myr, including both endpoints.

            Notes
            -----
            The histogram uses density=False and np.diff(bin_edges) for
            normalization, so empty SN lists give all-zero rates.
            Events outside the chosen time range are omitted by np.histogram;
            the caller checks all stored event times before histogramming.
            """

            hist, bin_edges = np.histogram([sn['explosion_time'] for sn in sn_events],
                                           weights=[sn[ejected_quantity] for sn in sn_events],
                                           bins=time_bin_count,
                                           range=(time_min, time_max),
                                           density=False,
                                          )

            bin_widths = np.diff(bin_edges)

            feedback_rate = hist / bin_widths

            return feedback_rate, bin_edges

        # One panel per feedback quantity:
        # ax_energy shows the energy rate dE/dt in erg/Myr,
        # ax_gas the gas return rate dM_gas/dt in Msun/Myr,
        # and ax_metal the metal ejection rate dM_metal/dt in Msun/Myr.
        # Each tuple holds the axis, the event key, the y label, the legend tag,
        # the sampled and Poisson-averaged totals, and the legend unit.
        # All panels share the same time bins.
        panels = (
                     (ax_energy,
                      'ejected_energy',
                      r'$\frac{{\rm d}E}{{\rm d\,} t}\quad [{\rm erg/Myr}]$',
                      'E',
                      self.sampled_sn_energy,
                      self.poisson_weighted_sn_energy,
                      ' erg'),
                     (ax_gas,
                      'ejected_gas_mass',
                      r'$\frac{{\rm d}M_{\rm gas}}{{\rm d\,} t}\quad [{\rm M_\odot/Myr}]$',
                      'M_gas',
                      self.sampled_ejected_gas_mass,
                      self.poisson_weighted_ejected_gas_mass,
                      r' M$_\odot$'),
                     (ax_metal,
                      'ejected_metal_mass',
                      r'$\frac{{\rm d}M_{\rm metal}}{{\rm d\,} t}\quad [{\rm M_\odot/Myr}]$',
                      'M_metal',
                      self.sampled_ejected_metal_mass,
                      self.poisson_weighted_ejected_metal_mass,
                      r' M$_\odot$'),
                 )

        for ax, key, ylabel, tag, sampled_total, averaged_total, unit in panels:
            # sampled_rate is the rate of the random catalog; time_bin_edges
            # are the time-bin limits shared by all three panels.
            sampled_rate, time_bin_edges = feedback_rate_histogram(self.sampled_sn_events, key)
            ax.stairs(np.maximum(sampled_rate, LOG_PLOT_FLOOR), time_bin_edges,
                      color='r',
                      label=f"Rand. Sampled, {tag}={sampled_total:.2e}{unit}",
                     )

            # Ensemble-mean rate per one-Myr bin: the conditional rate of each
            # fixed-count catalog k weighted by its Poisson probability.
            averaged_rate = np.zeros(time_bin_count)
            for sn_count in range(1, self.max_sn_count + 1):
                # The ignored second return (_) is the common set of bin edges.
                this_rate, _ = feedback_rate_histogram(self.fixed_count_sn_events[sn_count], key)
                averaged_rate += this_rate * self.fixed_count_poisson_weights[sn_count]

            if self.max_sn_count > 0:
                ax.stairs(np.maximum(averaged_rate, LOG_PLOT_FLOOR), time_bin_edges,
                          color='g',
                          label=f"Det. Averaged, {tag}={averaged_total:.2e}{unit}",
                         )
            ax.set_ylabel(ylabel)
            format_rate_axes(ax, sampled_rate, averaged_rate)

        # Save the figure to the output directory.
        plt.tight_layout()
        if output_dir is None:
            output_dir = '.'
        os.makedirs(output_dir, exist_ok=True)
        fig.savefig(os.path.join(output_dir, f'fig_CCSNe_feedback_time_evolution__{output_suffix}.png'))
        plt.close(fig)

        return

    def print_summary(self, show_fixed_counts: bool = True) -> None:
        """Print input parameters and expected, sampled, and fixed-count CCSN results.

        Parameters
        ----------
        show_fixed_counts : bool, optional
            Print one block per deterministic fixed count if True (default).
            With about 1000 expected SNe that is over a thousand blocks, so
            False prints only the Poisson-averaged totals.

        Notes
        -----
        Expected feedback comes from _sn_yields() integrals, while AGORA
        numbers are separate fixed-yield comparison values.
        If there are no sampled SNe, the rate is zero for positive
        sampled mass but its reciprocal is undefined.
        An empty sampled population has an undefined number of SNe per
        unit formed mass.
        These printed ratios are frequencies per mass, not time rates;
        time-dependent rates appear in the PNG figure.
        Random seeds may be integers or None.
        If the fixed-count catalogs were skipped (expected_sn_count above
        1000), no per-count block is printed and a "Skipped" line replaces
        the averaged quantities.
        """

        separation_line = "-" * 70

        print("\nInput Parameters")
        print(separation_line)
        print(f"Target particle mass                    :  {self.target_particle_mass:12.4e} Msun")
        print(f"IMF                                     :  {self.imf_type:>12s}")
        print(f"Stellar metallicity                     :  {self.stellar_metallicity:12.4e}")
        print(f"Random seed                             :  {str(self.random_seed):>12s}")
        print(separation_line)

        print("\nConstants")
        print(separation_line)
        print(f"M_MIN                                   :  {self.M_MIN:12.4e} Msun")
        print(f"M_MAX                                   :  {self.M_MAX:12.4e} Msun")
        print(f"KROUPA_M_BREAK                          :  {self.KROUPA_M_BREAK:12.4e} Msun")
        print(f"KROUPA_ALPHA_LOWER                      :  {self.KROUPA_ALPHA_LOWER:12.4e}")
        print(f"KROUPA_ALPHA_UPPER                      :  {self.KROUPA_ALPHA_UPPER:12.4e}")
        print(f"CHABRIER_M_C                            :  {self.CHABRIER_M_C:12.4e} Msun")
        print(f"CHABRIER_LOG10_SIGMA                    :  {self.CHABRIER_LOG10_SIGMA:12.4e}")
        print(f"CHABRIER_M_BREAK                        :  {self.CHABRIER_M_BREAK:12.4e} Msun")
        print(f"CHABRIER_ALPHA_UPPER                    :  {self.CHABRIER_ALPHA_UPPER:12.4e}")
        print(f"M_SN_MIN                                :  {self.M_SN_MIN:12.4e} Msun")
        print(f"M_SN_YIELD_BREAK                        :  {self.M_SN_YIELD_BREAK:12.4e} Msun")
        print(f"M_SN_MAX                                :  {self.M_SN_MAX:12.4e} Msun")
        print(f"COMPACT_REMNANT_MASS                    :  {self.COMPACT_REMNANT_MASS:12.4e} Msun")
        print(separation_line)

        print("\nExpected Quantities")
        print(separation_line)
        print(f"Expected number of stars                :  {self.expected_star_count:12.2f}")
        print(f"Expected number of supernovae           :  {self.expected_sn_count:12.2f}")
        print(f"Expected total stellar mass             :  {self.expected_stellar_mass:12.4e} Msun")
        if self.expected_stellar_mass > 0 and self.expected_sn_count > 0:
            print(f"Expected SNe per formed mass            :  {self.expected_sn_count/self.expected_stellar_mass:12.4e} Msun^-1 = 1 SN per {self.expected_stellar_mass/self.expected_sn_count:12.4e} Msun")
        elif self.expected_stellar_mass > 0:
            print(f"Expected SNe per formed mass            :  {0.0:12.4e} Msun^-1 (no SNe; reciprocal undefined)")
        else:
            print(f"Expected SNe per formed mass            :  undefined (zero stellar mass)")
        print(f"Expected total SN progenitor mass       :  {self.expected_sn_progenitor_mass:12.4e} Msun")
        print(f"Expected mean SN explosion time         :  {self.expected_sn_mean_time:12.4e} Myr")
        print(f"Expected total SN ejected energy        :  {self.expected_sn_energy:12.4e} erg")
        print(f"Expected total SN ejected gas mass      :  {self.expected_ejected_gas_mass:12.4e} Msun")
        print(f"Expected total SN ejected metal mass    :  {self.expected_ejected_metal_mass:12.4e} Msun")
        print(f"AGORA expected number of SNe            :  {self.agora_expected_sn_count:12.2f}")
        print(f"AGORA total SN ejected energy           :  {self.agora_expected_sn_energy:12.4e} erg")
        print(f"AGORA total SN ejected gas mass         :  {self.agora_expected_ejected_gas_mass:12.4e} Msun")
        print(f"AGORA total SN ejected metal mass       :  {self.agora_expected_ejected_metal_mass:12.4e} Msun")
        print(separation_line)

        print("\nSampled Quantities")
        print(separation_line)
        print(f"Sampled number of stars                 :  {self.sampled_star_count:12d}")
        print(f"Sampled number of supernovae            :  {self.sampled_sn_count:12d}")
        print(f"Sampled total stellar mass              :  {self.sampled_stellar_mass:12.4e} Msun")
        if self.sampled_stellar_mass > 0 and self.sampled_sn_count > 0:
            print(f"Sampled SNe per formed mass             :  {self.sampled_sn_count/self.sampled_stellar_mass:12.4e} Msun^-1 = 1 SN per {self.sampled_stellar_mass/self.sampled_sn_count:12.4e} Msun")
        elif self.sampled_stellar_mass > 0:
            print(f"Sampled SNe per formed mass             :  {0.0:12.4e} Msun^-1 (no SNe; reciprocal undefined)")
        else:
            print(f"Sampled SNe per formed mass             :  undefined (zero sampled stellar mass)")
        print(f"Sampled total SN progenitor mass        :  {self.sampled_sn_progenitor_mass:12.4e} Msun")
        print(f"Sampled total SN ejected energy         :  {self.sampled_sn_energy:12.4e} erg")
        print(f"Sampled total SN ejected gas mass       :  {self.sampled_ejected_gas_mass:12.4e} Msun")
        print(f"Sampled total SN ejected metal mass     :  {self.sampled_ejected_metal_mass:12.4e} Msun")
        print(separation_line)

        for sn_count in (range(1, self.max_sn_count + 1) if show_fixed_counts else ()):
            print(f"\nQuantities of Supernovae for Deterministic Fixed Count = {sn_count:02d}")
            print(separation_line)
            print(f"        Number of supernovae            :  {self.fixed_count_sn_counts[sn_count]:12d}")
            print(f"        Total SN progenitor mass        :  {self.fixed_count_progenitor_masses[sn_count]:12.4e} Msun")
            print(f"        Total SN ejected energy         :  {self.fixed_count_energies[sn_count]:12.4e} erg")
            print(f"        Total SN ejected gas mass       :  {self.fixed_count_gas_masses[sn_count]:12.4e} Msun")
            print(f"        Total SN ejected metal mass     :  {self.fixed_count_metal_masses[sn_count]:12.4e} Msun")
            print(separation_line)

        print("\nAveraged Quantities of All Supernovae with a Deterministic Fixed Count")
        print(separation_line)
        if self.max_sn_count > 0:
            print(f"Number of supernovae                    :  {self.poisson_weighted_sn_count:12.2f}")
            print(f"Total SN progenitor mass                :  {self.poisson_weighted_progenitor_mass:12.4e} Msun")
            print(f"Mean SN explosion time                  :  {self.poisson_weighted_mean_sn_time:12.4e} Myr")
            print(f"Total SN ejected energy                 :  {self.poisson_weighted_sn_energy:12.4e} erg")
            print(f"Total SN ejected gas mass               :  {self.poisson_weighted_ejected_gas_mass:12.4e} Msun")
            print(f"Total SN ejected metal mass             :  {self.poisson_weighted_ejected_metal_mass:12.4e} Msun")
        else:
            print("Skipped (expected SN count exceeded fixed-count evaluation threshold)")
        print(separation_line)

        return

    def write_feedback_tables(self, output_suffix: str, output_dir: Optional[str] = None) -> None:
        """Write sampled and deterministic fixed-count SN event tables.

        Parameters
        ----------
        output_suffix : str
            Text suffix included in every generated catalog filename.
        output_dir : str or None, optional
            Directory receiving all tables; created if missing.
            The default is "Tables__<output_suffix>" in the current
            working directory.

        Notes
        -----
        Each file stores five columns:
        explosion time [Myr], progenitor initial mass [Msun],
        returned gas [Msun], metals [Msun], and energy [erg].
        Zero-count tables contain only the header.
        Rows are sorted by
        time. Files with matching names inside output_dir are overwritten.
        Every file begins with commented header lines: two provenance lines
        (round-trip-safe input mass/metallicity, IMF, seed, SN limits, and a short
        yield description), one catalog line, and the line of column names,
        which carry their units.
        np.loadtxt skips these lines as comments.
        The catalog line of a fixed-count table gives its count k,
        a slice_averaged=yields+times flag, present only if
        slice_averaged_yields is True (explosion time, gas and metals are then
        slice means, not evaluated at the progenitor mass in the same row), its
        Poisson weight P(k) (for k = max_sn_count > 0 this is the upper tail
        P(N >= k)), the expected SN count lambda and max_sn_count, so the
        tables can be recombined into the Poisson average without the
        original object.
        One extra manifest file,
        Table_CCSNe_fixed_count_manifest__<output_suffix>.txt, lists the same
        weights with the per-count totals in one place (columns: k, Poisson
        weight, event count, progenitor mass [Msun], gas [Msun], metals [Msun],
        energy [erg], first and last explosion time [Myr]; the times are NaN
        for k = 0), so the Poisson average needs no parsing of the individual
        headers.
        The catalog line of the sampled table gives the realized star count,
        stellar mass and SN count.
        The yield description is hard-coded text and must be updated
        if _sn_yields changes.
        Values are written with 17 significant digits (%27.16e).
        If the fixed-count catalogs were skipped (expected_sn_count above
        1000), only the sampled-star table and the empty zero-count
        table, and the one-row manifest are written.
        """

        # Keep the table provenance consistent with the event-yield model.
        provenance = (
            f"IMF={self.imf_type} stellar_metallicity={self.stellar_metallicity:.17g} "
            f"target_particle_mass={self.target_particle_mass:.17g}Msun random_seed={self.random_seed}\n"
            f"M_SN_MIN={self.M_SN_MIN:g} M_SN_YIELD_BREAK={self.M_SN_YIELD_BREAK:g} "
            f"M_SN_MAX={self.M_SN_MAX:g} yields=Kim+2014 fits below break, constants above")

        def write_sn_table(filename: str, sn_events: List[Dict[str, Any]], catalog_header: str) -> None:
            """Write one SN catalog as a five-column ASCII text file.

            Parameters
            ----------
            filename : str
                File path to create or overwrite.
            sn_events : list of dict
                Event dictionaries containing five quantity keys.
            catalog_header : str
                One-line description of this catalog (count, weight, ...)

            Notes
            -----
            np.array assembles time, initial mass, gas, metals, and energy
            from event dictionaries;
            np.savetxt writes a header-only table for N=0.
            The provenance and catalog lines are prepended to
            the column-name header.
            """

            if not sn_events:
                data_matrix = np.empty((0, 5))
            else:
                data_matrix = np.array([[sn['explosion_time'],
                                         sn['progenitor_mass'],
                                         sn['ejected_gas_mass'],
                                         sn['ejected_metal_mass'],
                                         sn['ejected_energy']] for sn in sn_events])

            np.savetxt(filename,
                       data_matrix,
                       fmt='%27.16e',
                       delimiter=' ',
                       header=f"{provenance}\n{catalog_header}\n{'explosion_time_in_Myr':>25}{'progenitor_mass_in_Msun':>28}{'ejected_gas_mass_in_Msun':>28}{'ejected_metal_mass_in_Msun':>28}{'ejected_energy_in_erg':>28}",
                      )

        # Default to one directory per run.
        if output_dir is None:
            output_dir = f"Tables__{output_suffix}"
        os.makedirs(output_dir, exist_ok=True)

        write_sn_table(os.path.join(output_dir, f"Table_CCSNe_in_randomly_sampled_all_stars__{output_suffix}.txt"),
                       self.sampled_sn_events,
                       f"catalog=randomly_sampled_all_stars sampled_star_count={self.sampled_star_count} "
                       f"sampled_stellar_mass={self.sampled_stellar_mass:.17e}Msun sampled_sn_count={self.sampled_sn_count}")

        # Flag tables whose yields are not evaluated at the progenitor mass of the row.
        yield_flag = " slice_averaged=yields+times" if self.slice_averaged_yields else ""

        # One-file index of the fixed-count grid: Poisson weights and per-count totals.
        manifest_columns = np.column_stack((np.arange(self.max_sn_count + 1),
                                            self.fixed_count_poisson_weights,
                                            self.fixed_count_sn_counts,
                                            self.fixed_count_progenitor_masses,
                                            self.fixed_count_gas_masses,
                                            self.fixed_count_metal_masses,
                                            self.fixed_count_energies,
                                            self.fixed_count_first_sn_times,
                                            self.fixed_count_last_sn_times))
        np.savetxt(os.path.join(output_dir, f"Table_CCSNe_fixed_count_manifest__{output_suffix}.txt"),
                   manifest_columns,
                   fmt=['%13d', '%27.16e', '%9d'] + ['%27.16e'] * 6,
                   header=(f"{provenance}\n"
                           f"catalog=fixed_count_manifest expected_sn_count={self.expected_sn_count:.17e} "
                           f"max_sn_count={self.max_sn_count} (the last weight is the upper tail P(N >= k) if k > 0)\n"
                           "fixed_count              poisson_weight  sn_count     progenitor_mass_in_Msun    ejected_gas_mass_in_Msun  "
                           "ejected_metal_mass_in_Msun       ejected_energy_in_erg           first_time_in_Myr            last_time_in_Myr"))

        # sn_count selects the fixed-count table filename and catalog.
        for sn_count in range(0, self.max_sn_count + 1):
            write_sn_table(os.path.join(output_dir, f"Table_CCSNe_with_deterministic_fixed_count_{sn_count:04d}__{output_suffix}.txt"),
                           self.fixed_count_sn_events[sn_count],
                           f"catalog=deterministic_fixed_count fixed_count={sn_count}{yield_flag} "
                           f"poisson_weight={self.fixed_count_poisson_weights[sn_count]:.17e} "
                           f"weight_is_upper_tail={bool(sn_count == self.max_sn_count and sn_count > 0)} "
                           f"expected_sn_count={self.expected_sn_count:.17e} max_sn_count={self.max_sn_count}")

        return


def format_roundtrip(value: float) -> str:
    """Format a float with the shortest text that reads back as the same float.

    Parameters
    ----------
    value : float
        Number to format.

    Returns
    -------
    str
        repr of the float without a trailing ".0", e.g. "1000", "0.1", "1e-05".
    """

    text = repr(float(value))

    return text[:-2] if text.endswith('.0') else text


def main(argv: Optional[List[str]] = None) -> int:
    """Run the calculation from the command line; write catalogs and figures.

    Parameters
    ----------
    argv : list of str or None, optional
        Command-line arguments; None means sys.argv[1:].
        Options: --mass (target particle mass in Msun, default 1000),
        --imf (Kroupa or Chabrier, default Kroupa),
        --seed (integer, default 123),
        --metallicity (mass fraction Z in [0.0004, 0.03], default 0.02),
        --slice-averaged-yields (unbiased gas/metal yields in the
        fixed-count catalogs; off by default), and
        --output-dir (directory for tables and figures; default: tables in
        Tables__<suffix>, figures in the current directory), and
        --brief (omit the per-count blocks of the printed summary),
        --no-plots (write tables without PNG figures), and
        --allow-no-fixed-count (accept a run without fixed-count catalogs).

    Returns
    -------
    int
        Process exit status: 0 on success; 3 if the fixed-count catalogs
        were skipped (expected SN count above 1000) and
        --allow-no-fixed-count was not given.
        All outputs that could be produced
        (sampled table, zero-count table, figures) are still
        written in that case.

    Notes
    -----
    The command-line arguments select population mass, IMF, random seed,
    metallicity, and output directory.
    The output suffix includes all four values to distinguish metallicity
    choices (plus _SliceAveragedYields when that option is on), but the
    numeric values use shortest round-trip-safe text (format_roundtrip)
    so distinct float inputs cannot collide through filename rounding;
    model-code revisions are not encoded in filenames.
    Outputs are produced in this order: tables, printed summary, then
    optional IMF/rate figures; --no-plots avoids figure generation in
    batch production and does not change the numerical catalogs.
    Tables go to the directory Tables__<suffix>, figures to the
    current working directory unless --output-dir is given, in which case
    both go to that directory.
    Re-running identical settings overwrites matching output files.
    If the expected SN count exceeds 1000, the fixed-count catalogs are skipped
    (with a warning): only the sampled-star and zero-count tables are written,
    and the fixed-count curves in the figures are empty
    (see StarParticle._build_fixed_count_sn_events).
    Since the fixed-count tables are the main product, this is reported on stderr
    and through exit status 3 so that scripts do not mistake the incomplete
    output for a success.
    main() switches matplotlib to the non-interactive Agg backend; merely
    importing this module leaves the caller's backend untouched.
    """

    parser = argparse.ArgumentParser(description="Generate CCSN feedback catalogs for a coeval population.")
    parser.add_argument('--mass', type=float, default=1000.0, help='Target particle mass [Msun]')
    parser.add_argument('--imf', choices=['Kroupa', 'Chabrier'], default='Kroupa')
    parser.add_argument('--seed', type=int, default=123)
    parser.add_argument('--metallicity', type=float, default=0.02, help='mass fraction Z in [0.0004, 0.03]')
    parser.add_argument('--slice-averaged-yields', action='store_true',
                        help='store slice-averaged gas/metal yields in the fixed-count catalogs instead of '
                             'yields at the slice-mean masses')
    parser.add_argument('--output-dir', default=None,
                        help='directory for tables and figures (default: Tables__<suffix> for tables; current directory for figures)')
    parser.add_argument('--brief', action='store_true',
                        help='print only the Poisson-averaged fixed-count totals, not one block per count')
    parser.add_argument('--no-plots', action='store_true',
                        help='skip PNG figures; generate only feedback tables and the text summary')
    parser.add_argument('--allow-no-fixed-count', action='store_true',
                        help='exit with status 0 even if the expected SN count is too large for the '
                             'fixed-count catalogs, which are then not written')
    args = parser.parse_args(argv)

    # The command line only writes files, so use the non-interactive backend
    # here rather than at import time, which would override the backend of any
    # session (a notebook, say) that imports StarParticle as a library.
    plt.switch_backend('Agg')

    input_target_particle_mass  = args.mass
    input_imf_type              = args.imf
    input_random_seed           = args.seed
    input_stellar_metallicity   = args.metallicity
    input_slice_averaged_yields = args.slice_averaged_yields

    # Use shortest round-trip-safe float formatting
    output_suffix = (f"{input_imf_type}IMF_ParticleMass{format_roundtrip(input_target_particle_mass)}"
                     f"_RandomSeed{input_random_seed:05d}_Z{format_roundtrip(input_stellar_metallicity)}")
    if input_slice_averaged_yields:
        # Keep the two yield modes from overwriting each other's files.
        output_suffix += "_SliceAveragedYields"

    # Build the population before writing any output.
    try:
        population = StarParticle(target_particle_mass=input_target_particle_mass,
                                  random_seed=input_random_seed,
                                  imf_type=input_imf_type,
                                  stellar_metallicity=input_stellar_metallicity,
                                  slice_averaged_yields=input_slice_averaged_yields,
                                 )
    except (ValueError, MemoryError) as err:
        parser.error(str(err))

    population.write_feedback_tables(output_suffix, args.output_dir)
    population.print_summary(show_fixed_counts=not args.brief)
    if not args.no_plots:
        population.plot_imf(output_suffix, args.output_dir)
        population.plot_sn_feedback_rates(output_suffix, args.output_dir)

    if population.max_sn_count == 0 and not args.allow_no_fixed_count:
        print(f"ERROR: no fixed-count catalogs were written (expected SN count "
              f"{population.expected_sn_count:.0f} exceeds the limit); only the sampled "
              f"catalog exists. Use a smaller --mass or pass --allow-no-fixed-count.",
              file=sys.stderr)
        return 3

    return 0


if __name__ == '__main__':
    sys.exit(main())
