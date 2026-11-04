"""ICRU 90 tabulated electronic mass stopping powers in liquid water.

All values are embedded directly in this Python module (no external text files).
Energies are kinetic energies in MeV; stopping powers are in MeV cm2/g.
The tables are transcribed from the supplied electron/proton/alpha TXT files.

At tabulated energies the exact table value is returned; between entries,
log-log interpolation is applied. No extrapolation outside the tables.

CSDA ranges (g/cm2) are numerically integrated from these stopping powers,
with the original low-energy approximation retained. They are not tabulated
reference ranges.
"""
import bisect
import math

# Tables embedded from the three supplied ICRU90 files.
ELECTRON_ENERGY = [
    0.0010, 0.0015, 0.0020, 0.0030, 0.0040, 0.0050, 0.0060, 0.0080,
    0.0100, 0.0150, 0.0200, 0.0300, 0.0400, 0.0500, 0.0600, 0.0800,
    0.1000, 0.1500, 0.2000, 0.3000, 0.4000, 0.5000, 0.6000, 0.8000,
    1.0000, 1.5000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 8.0000,
    10.0000, 15.0000, 20.0000, 30.0000, 40.0000, 50.0000, 60.0000, 80.0000,
    100.0000, 150.0000, 200.0000, 300.0000, 400.0000, 500.0000, 600.0000, 800.0000,
    1000.0000,
]

ELECTRON_STOPPING = [
    1.181E+02, 9.063E+01, 7.436E+01, 5.563E+01, 4.500E+01, 3.806E+01, 3.315E+01, 2.659E+01,
    2.239E+01, 1.635E+01, 1.308E+01, 9.591E+00, 7.729E+00, 6.564E+00, 5.763E+00, 4.731E+00,
    4.093E+00, 3.221E+00, 2.779E+00, 2.344E+00, 2.138E+00, 2.025E+00, 1.956E+00, 1.880E+00,
    1.845E+00, 1.819E+00, 1.821E+00, 1.844E+00, 1.869E+00, 1.891E+00, 1.910E+00, 1.942E+00,
    1.967E+00, 2.013E+00, 2.045E+00, 2.088E+00, 2.117E+00, 2.139E+00, 2.156E+00, 2.182E+00,
    2.202E+00, 2.238E+00, 2.263E+00, 2.298E+00, 2.322E+00, 2.341E+00, 2.357E+00, 2.382E+00,
    2.401E+00,
]

PROTON_ENERGY = [
    0.0010, 0.0015, 0.0020, 0.0030, 0.0040, 0.0050, 0.0060, 0.0080,
    0.0100, 0.0150, 0.0200, 0.0300, 0.0400, 0.0500, 0.0600, 0.0800,
    0.1000, 0.1500, 0.2000, 0.3000, 0.4000, 0.5000, 0.6000, 0.8000,
    1.0000, 1.5000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 8.0000,
    10.0000, 15.0000, 20.0000, 30.0000, 40.0000, 50.0000, 60.0000, 80.0000,
    100.0000, 150.0000, 200.0000, 300.0000, 400.0000, 500.0000, 600.0000, 800.0000,
    1000.0000, 1500.0000, 2000.0000, 3000.0000, 4000.0000, 5000.0000, 6000.0000, 8000.0000,
    10000.0000,
]

PROTON_STOPPING = [
    1.337E+02, 1.638E+02, 1.891E+02, 2.316E+02, 2.675E+02, 2.990E+02, 3.276E+02, 3.782E+02,
    4.229E+02, 5.036E+02, 5.673E+02, 6.628E+02, 7.290E+02, 7.740E+02, 8.026E+02, 8.241E+02,
    8.145E+02, 7.360E+02, 6.585E+02, 5.435E+02, 4.643E+02, 4.065E+02, 3.624E+02, 2.997E+02,
    2.574E+02, 1.934E+02, 1.569E+02, 1.160E+02, 9.319E+01, 7.842E+01, 6.801E+01, 5.417E+01,
    4.532E+01, 3.269E+01, 2.589E+01, 1.864E+01, 1.479E+01, 1.238E+01, 1.072E+01, 8.578E+00,
    7.250E+00, 5.417E+00, 4.470E+00, 3.504E+00, 3.018E+00, 2.731E+00, 2.544E+00, 2.323E+00,
    2.203E+00, 2.065E+00, 2.017E+00, 1.998E+00, 2.010E+00, 2.029E+00, 2.050E+00, 2.090E+00,
    2.124E+00,
]

ALPHA_ENERGY = [
    0.0010, 0.0015, 0.0020, 0.0030, 0.0040, 0.0050, 0.0060, 0.0080,
    0.0100, 0.0150, 0.0200, 0.0300, 0.0400, 0.0500, 0.0600, 0.0800,
    0.1000, 0.1500, 0.2000, 0.3000, 0.4000, 0.5000, 0.6000, 0.8000,
    1.0000, 1.5000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 8.0000,
    10.0000, 15.0000, 20.0000, 30.0000, 40.0000, 50.0000, 60.0000, 80.0000,
    100.0000, 150.0000, 200.0000, 300.0000, 400.0000, 500.0000, 600.0000, 800.0000,
    1000.0000,
]

ALPHA_STOPPING = [
    9.891E+01, 1.228E+02, 1.431E+02, 1.776E+02, 2.069E+02, 2.330E+02, 2.568E+02, 2.993E+02,
    3.370E+02, 4.181E+02, 4.871E+02, 6.039E+02, 7.030E+02, 7.905E+02, 8.696E+02, 1.009E+03,
    1.131E+03, 1.383E+03, 1.582E+03, 1.873E+03, 2.062E+03, 2.178E+03, 2.240E+03, 2.256E+03,
    2.190E+03, 1.877E+03, 1.599E+03, 1.239E+03, 1.021E+03, 8.747E+02, 7.687E+02, 6.239E+02,
    5.290E+02, 3.894E+02, 3.118E+02, 2.267E+02, 1.802E+02, 1.506E+02, 1.300E+02, 1.030E+02,
    8.593E+01, 6.194E+01, 4.923E+01, 3.586E+01, 2.885E+01, 2.452E+01, 2.156E+01, 1.780E+01,
    1.550E+01,
]

T0 = 0.001               # MeV, lowest energy of every table
RANGE_MIN_MEV = {"e-": 0.01, "proton": 0.1, "alpha": 0.5}  # CSDA ranges reliable from here up
STEPS_PER_DECADE = 400   # integration grid of the CSDA range

PARTICLES = ("e-", "proton", "alpha")


def _interpolate(x, xs, ys, log=True):
    """Interpolation in a table, log-log (stopping powers) or linear."""
    i = min(max(bisect.bisect_right(xs, x) - 1, 0), len(xs) - 2)
    x0, x1, y0, y1 = xs[i], xs[i + 1], ys[i], ys[i + 1]
    if x == x0:
        return y0
    if x == x1:
        return y1
    if log:
        f = math.log(x / x0) / math.log(x1 / x0)
        return y0 * (y1 / y0) ** f
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def energy_limits(particle):
    """Energy range (MeV) of the reference for this particle."""
    if particle == "e-":
        return ELECTRON_ENERGY[0], ELECTRON_ENERGY[-1]
    if particle == "proton":
        return PROTON_ENERGY[0], PROTON_ENERGY[-1]
    if particle == "alpha":
        return ALPHA_ENERGY[0], ALPHA_ENERGY[-1]
    raise ValueError("No ICRU 90 data for particle %r (only %s)" % (particle, ", ".join(PARTICLES)))


def covers(particle, energy):
    low, high = energy_limits(particle)
    return low <= energy <= high


def stopping_power(particle, energy):
    """ICRU 90 mass stopping power (MeV cm2/g) of liquid water at energy (MeV)."""
    if not covers(particle, energy):
        raise ValueError("%g MeV is outside the ICRU 90 table for %s" % (energy, particle))
    if particle == "e-":
        return _interpolate(energy, ELECTRON_ENERGY, ELECTRON_STOPPING)
    if particle == "proton":
        return _interpolate(energy, PROTON_ENERGY, PROTON_STOPPING)
    return _interpolate(energy, ALPHA_ENERGY, ALPHA_STOPPING)


def csda_range(particle, energy):
    """CSDA range (g/cm2) in liquid water from the ICRU 90 stopping powers."""
    if not covers(particle, energy):
        raise ValueError("%g MeV is outside the ICRU 90 table for %s" % (energy, particle))
    s0 = stopping_power(particle, T0)
    result = (T0 if particle == "e-" else 2 * T0) / s0
    # Simpson's rule in u = ln(T): dR = T/S(T) du.
    n = max(2, 2 * math.ceil(STEPS_PER_DECADE * math.log10(energy / T0) / 2))
    if energy > T0:
        h = math.log(energy / T0) / n
        total = 0.0
        for i in range(n + 1):
            t = T0 * math.exp(i * h) if i < n else energy
            weight = 1 if i in (0, n) else (4 if i % 2 else 2)
            total += weight * t / stopping_power(particle, t)
        result += total * h / 3
    return result
