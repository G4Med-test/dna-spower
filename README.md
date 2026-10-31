# dna-spower

**Stopping power of electrons, protons and alpha particles in liquid water**
with the Geant4-DNA physics in stationary mode, compared with ICRU Report 90.
Built on the Geant4 example
[`extended/medical/dna/spower`](https://github.com/Geant4/geant4/tree/v11.3.2/examples/extended/medical/dna/spower).

---

## 🔬 What this test does

Mono-energetic particles start at the centre of a 1 m radius sphere of
`G4_WATER`, in a random direction. The `G4EmDNAPhysics_stationary*` constructors
keep the primary energy constant; `SteppingAction` kills secondary electrons and
photons, and after 1000 steps (e-, protons and H) or 10000 steps (alphas and
their charge states) scores energy deposited / path length and ends the event.
For each energy the application writes to `spower.txt` the mean stopping power
(keV/µm) and its standard deviation over the events.

| Macro | Particle | Physics | Energies | Events per energy |
| --- | --- | --- | --- | --- |
| `e-_dna_stat_opt2.mac` | e- | `dna_stat_opt2` | 50 eV – 1 MeV | 2000 |
| `proton_dna_stat_opt2.mac` | proton | `dna_stat_opt2` | 1 keV – 50 MeV | 2000 |
| `alpha_dna_stat_opt2.mac` | alpha | `dna_stat_opt2` | 1 keV – 100 MeV | 2000 |

`validation/parser.py` publishes one plot per particle (x: energy in MeV, y:
stopping power in keV/µm, error: standard deviation / √N) with the **ICRU 90**
stopping power at the simulated energies it covers (from 1 keV), converted with
the `G4_WATER` density of 1 g/cm³.

## 📚 Reference: ICRU 90

`validation/icru90.py` (identical in `dna-range`):

- **protons, alphas**: electronic stopping power of `G4_WATER` from
  `G4ICRU90StoppingData` (the ICRU 90 table embedded in Geant4,
  `source/materials/src/G4ICRU90StoppingData.cc`, 1 keV – 10 GeV for protons,
  1 keV – 1 GeV for alphas), log-log interpolated. Nuclear stopping is not
  included: it exceeds about 1% of the electronic stopping below 10 keV
  (protons) and 100 keV (alphas).
- **electrons**: collision stopping power from the Bethe formula (ICRU 37
  eq. 2.1) with the ICRU 90 mean excitation energy I = 78 eV and density-effect
  correction (zero below 0.5 MeV), as ESTAR computes it (with I = 75 eV it
  reproduces ESTAR to 4 digits). Geant4-DNA has no bremsstrahlung, so the
  collision, not the total, stopping power is the right reference. The Bethe
  formula has no shell corrections and becomes less reliable at a few keV.

## Import notes

- `main.cc`, `include/`, `src/` are `spower.cc`, `include/` and `src/` of the
  Geant4 v11.3.2 example (the version of the default base image), unchanged.
- `macro/spower.in` is the example's macro, kept for reference; the validation
  macros are the `*_dna_stat_opt2.mac` files listed in `validation/config.json`.
  `docs/geant4-example-README.txt` is the example's README.
- `macro/unit.mac` runs one 1 keV electron on every push.

For the CI, the validation flow and running everything locally, see the
[template README](https://github.com/G4Med-test/template#readme).

References for this example: Nucl. Instrum. Meth. B 397 (2017) 45–50, and the
Geant4-DNA publications Med. Phys. 51 (2024) 5873–5889; Med. Phys. 45 (2018)
e722–e739; Phys. Med. 31 (2015) 861–874; Med. Phys. 37 (2010) 4692–4708;
Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178.
