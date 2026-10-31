"""ICRU Report 90 stopping powers and CSDA ranges of e-, protons and alphas in liquid water.

Shared by the dna-range and dna-spower parsers (keep the two copies identical).
Everything is pure Python: energies in MeV, mass stopping powers in MeV cm2/g,
CSDA ranges in g/cm2.

Protons and alphas: electronic stopping power, G4_WATER rows (index 1) of
G4ICRU90StoppingData, the ICRU 90 table embedded in Geant4 (Geant4 v11.3.2,
source/materials/src/G4ICRU90StoppingData.cc). The alpha table is in alpha
kinetic energy. Nuclear stopping is not included: it is below about 1% of the
electronic stopping above 10 keV for protons and above 100 keV for alphas.

Electrons: collision stopping power from the Bethe formula of ICRU Report 37,
eq. 2.1, with the ICRU 90 mean excitation energy of liquid water, I = 78 eV, and
the ICRU 90 density-effect correction (EGSnrc water_icru90.density, computed
with ESTAR.f; zero below 0.5 MeV). This is how ESTAR and ICRU 90 compute it; with
I = 75 eV the same code reproduces the ESTAR table to 4 digits. Geant4-DNA has no
bremsstrahlung, so the collision stopping power, not the total, is the right
reference. The Bethe formula has no shell corrections and becomes less reliable
at a few keV.

CSDA range: R(T) = R(T0) + integral from T0 to T of dT'/S(T'), with T0 = 1 keV
(lowest tabulated energy). Below T0 the stopping power is taken proportional to
the velocity for protons and alphas, R(T0) = 2 T0/S(T0), and constant for
electrons, R(T0) = T0/S(T0). This estimate, and for ions the missing nuclear
stopping at low energy, change R by a few percent near the lowest energies only:
use the ranges from RANGE_MIN_MEV[particle] up.
"""
import bisect
import math

# Kinetic energy (MeV) and electronic mass stopping power (MeV cm2/g).
PROTON_ENERGY = [
    0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008, 0.01, 0.015, 0.02, 0.03, 0.04,
    0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0,
    6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 150.0, 200.0, 300.0,
    400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0,
    8000.0, 10000.0]
PROTON_STOPPING = [
    133.7, 163.8, 189.1, 231.6, 267.5, 299.0, 327.6, 378.2, 422.9, 503.6, 567.3, 662.8, 729.0,
    774.0, 802.6, 824.1, 814.5, 736.0, 658.5, 543.5, 464.3, 406.5, 362.4, 299.7, 257.4, 193.4,
    156.9, 116.0, 93.19, 78.42, 68.01, 54.17, 45.32, 32.69, 25.89, 18.64, 14.79, 12.38, 10.72,
    8.578, 7.25, 5.417, 4.47, 3.504, 3.018, 2.731, 2.544, 2.323, 2.203, 2.065, 2.017, 1.998,
    2.01, 2.029, 2.05, 2.09, 2.124]
ALPHA_ENERGY = PROTON_ENERGY[:49]
ALPHA_STOPPING = [
    98.91, 122.8, 143.1, 177.6, 206.9, 233.0, 256.8, 299.3, 337.0, 418.1, 487.1, 603.9, 703.0,
    790.5, 869.6, 1009.0, 1131.0, 1383.0, 1582.0, 1873.0, 2062.0, 2178.0, 2240.0, 2256.0,
    2190.0, 1877.0, 1599.0, 1239.0, 1021.0, 874.7, 768.7, 623.9, 529.0, 389.4, 311.8, 226.7,
    180.2, 150.6, 130.0, 103.0, 85.93, 61.94, 49.23, 35.86, 28.85, 24.52, 21.56, 17.8, 15.5]

# Electrons: Bethe formula constants and density-effect correction delta(T).
K = 0.307075             # 4 pi N_A r_e^2 m_e c^2, MeV cm2/mol
ELECTRON_MASS = 0.51099895  # MeV
Z_OVER_A = 10 / 18.0153  # mol/g, H2O
I_WATER = 78e-6          # MeV, ICRU 90
ELECTRON_MAX = 10.0      # MeV, end of the density-effect table below
DELTA_ENERGY = [0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0,
                4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0]
DELTA = [0.0, 6.126e-4, 1.501e-2, 5.338e-2, 1.005e-1, 1.529e-1, 2.086e-1, 3.535e-1, 4.982e-1,
         6.377e-1, 7.704e-1, 1.014, 1.231, 1.426, 1.601, 1.760, 1.906, 2.041, 2.165, 2.390,
         2.589, 2.767, 2.928]

T0 = 0.001               # MeV, lowest energy of every table
RANGE_MIN_MEV = {"e-": 0.01, "proton": 0.1, "alpha": 0.5}  # CSDA ranges reliable from here up
STEPS_PER_DECADE = 400   # integration grid of the CSDA range

PARTICLES = ("e-", "proton", "alpha")


def _interpolate(x, xs, ys, log=True):
    """Interpolation in a table, log-log (stopping powers) or linear."""
    i = min(max(bisect.bisect_right(xs, x) - 1, 0), len(xs) - 2)
    x0, x1, y0, y1 = xs[i], xs[i + 1], ys[i], ys[i + 1]
    if log:
        f = math.log(x / x0) / math.log(x1 / x0)
        return y0 * (y1 / y0) ** f
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def _electron_stopping(energy):
    tau = energy / ELECTRON_MASS
    gamma = tau + 1
    beta2 = 1 - 1 / gamma ** 2
    f_minus = 1 - beta2 + (tau ** 2 / 8 - (2 * tau + 1) * math.log(2)) / gamma ** 2
    delta = _interpolate(energy, DELTA_ENERGY, DELTA, log=False) if energy > DELTA_ENERGY[0] else 0.0
    return (K / 2 * Z_OVER_A / beta2
            * (math.log(tau ** 2 * (tau + 2) / (2 * (I_WATER / ELECTRON_MASS) ** 2)) + f_minus - delta))


def energy_limits(particle):
    """Energy range (MeV) of the reference for this particle."""
    if particle == "e-":
        return T0, ELECTRON_MAX
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
        return _electron_stopping(energy)
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
