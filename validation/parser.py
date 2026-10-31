#!/usr/bin/env python

# dna-spower portal parser: metadata() reads the macro, parse() the spower.txt
# table that Run::EndOfRun appends after every /run/beamOn (one macro = one
# particle, one energy per run): energy (eV), mean stopping power (keV/um) and
# its standard deviation over the events (keV/um).
#
# The stopping power is computed in stationary mode (G4EmDNAPhysics_stationary*):
# the primary does not lose energy, and SteppingAction aborts each event after
# 1000 steps (e-, protons) or 10000 steps (alphas) and scores
# energy deposited / path length. It is compared with the ICRU 90 stopping power
# of liquid water (validation/icru90.py) wherever the ICRU 90 table covers the
# beam energy.
import math
import os

import icru90
from geantval import energy_mev, getJSON, one_command

TEST = "dna-spower"
PHYSLISTS = {"dna_stat", "dna_stat_opt2", "dna_stat_opt4", "dna_stat_opt6"}
DENSITY = 1.0  # g/cm3, G4_WATER


def energy_scan(commands):
    """Energy (MeV) and number of events of every /run/beamOn, in macro order."""
    energies, events = [], []
    energy = None
    for command, value in commands:
        if command == "/gun/energy":
            energy = energy_mev(value)
        elif command == "/run/beamOn":
            n = int(value)
            if energy is None or n < 2:
                raise ValueError("Each /run/beamOn needs an explicit /gun/energy before it "
                                 "and at least 2 events")
            energies.append(energy)
            events.append(n)
    if not energies:
        raise ValueError("No /run/beamOn in the macro")
    return energies, events


def metadata(commands):
    physlist = one_command(commands, "/spower/phys/addPhysics")
    if physlist not in PHYSLISTS:
        # The application silently keeps its default for an unknown name.
        raise ValueError("Unknown physics constructor %r (expected one of %s)"
                         % (physlist, ", ".join(sorted(PHYSLISTS))))
    if one_command(commands, "/spower/det/setMat") != "G4_WATER":
        raise ValueError("The ICRU 90 comparison needs /spower/det/setMat G4_WATER")
    particle = one_command(commands, "/gun/particle")
    if particle not in icru90.PARTICLES:
        raise ValueError("Unsupported particle %r (expected one of %s)"
                         % (particle, ", ".join(icru90.PARTICLES)))
    energies, events = energy_scan(commands)
    return {"TEST": TEST, "PHYSLIST": physlist, "PARTICLE": particle,
            "ENERGIES": energies, "NEVENTS": events}


def read_table(path, expected_energies):
    """Rows of spower.txt, checked against the energies of the macro."""
    rows = []
    with open(path) as table:
        for line in table:
            if line.strip():
                rows.append([float(v) for v in line.split()])
    if len(rows) != len(expected_energies):
        raise ValueError("spower.txt has %d rows, expected %d (one per /run/beamOn)"
                         % (len(rows), len(expected_energies)))
    for row, energy in zip(rows, expected_energies):
        if len(row) != 3:
            raise ValueError("Expected 3 columns in spower.txt, got %r" % row)
        if not math.isclose(row[0] * 1e-6, energy, rel_tol=1e-5):
            raise ValueError("spower.txt energy %g eV does not match the macro (%g MeV)"
                             % (row[0], energy))
    return rows


def parse(job):
    rows = read_table(os.path.join(job["path"], "spower.txt"), job["ENERGIES"])
    particle = job["PARTICLE"]
    energies = job["ENERGIES"]
    stopping = [row[1] for row in rows]
    errors = [row[2] / math.sqrt(n) for row, n in zip(rows, job["NEVENTS"])]

    common = dict(observableName="stopping power",
                  targetName="water",
                  beamParticle=particle,
                  secondaryParticle="None",
                  title="Stopping power of %s in liquid water" % particle,
                  xAxisName="E, MeV",
                  yAxisName="S, keV/um")

    yield getJSON(job, "chart",
                  mctool_name="GEANT4",
                  mctool_model=job["PHYSLIST"],
                  beamEnergies=energies,
                  xValues=energies,
                  yValues=stopping,
                  yStatErrorsPlus=errors,
                  yStatErrorsMinus=errors,
                  **common)

    covered = [e for e in energies if icru90.covers(particle, e)]
    if covered:
        # MeV cm2/g * g/cm3 = MeV/cm = 0.1 keV/um
        reference = [icru90.stopping_power(particle, e) * DENSITY * 0.1 for e in covered]
        yield getJSON(job, "chart",
                      mctool_name="ICRU",
                      mctool_model="ICRU90",
                      beamEnergies=covered,
                      xValues=covered,
                      yValues=reference,
                      **common)
