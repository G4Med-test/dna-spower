# dna-spower

**Stopping power of electrons, protons and alpha particles in liquid water**
with the Geant4-DNA physics in stationary mode, compared with ICRU Report 90.
Built on the Geant4 example
[`extended/medical/dna/spower`](https://github.com/Geant4/geant4/tree/v11.3.2/examples/extended/medical/dna/spower).

---

## 🔬 What this test does

Mono-energetic particles start at the centre of a 1 m radius sphere of
`G4_WATER`, in a random direction. The `G4EmDNAPhysics*` constructors are used in stationary mode
keep the primary energy constant; `SteppingAction` kills secondary electrons and
photons, and after 1000 steps (e-, protons and H) or 10000 steps (alphas and
their charge states) scores energy deposited / path length and ends the event.
For each energy the application writes to `spower.txt` the mean stopping power
(keV/µm) and its standard deviation over the events.

| Macro | Particle | Physics | Energies | Events per energy |
| --- | --- | --- | --- | --- |
| `e-_dna_opt2.mac` | e- | `dna_opt2` | 10 eV – 1 MeV | 10000 |
| `e-_dna_opt4.mac` | e- | `dna_opt4` | 10 eV – 10 MeV | 10000 |
| `e-_dna_opt6.mac` | e- | `dna_opt6` | 10 eV – 1 MeV | 10000 |
| `e-_dna_opt8.mac` | e- | `dna_opt8` | 10 eV – 1 MeV | 10000 |
| `proton_dna_opt2.mac` | proton | `dna_opt2` | 100 eV – 300 MeV | 1000 |
| `proton_dna_opt4.mac` | proton | `dna_opt4` | 100 eV – 300 MeV | 1000 |
| `proton_dna_opt6.mac` | proton | `dna_opt6` | 100 eV – 300 MeV | 1000 |
| `proton_dna_opt8.mac` | proton | `dna_opt8` | 100 eV – 300 MeV | 1000 |
| `alpha_dna_opt2.mac` | alpha | `dna_opt2` | 100 eV – 400 MeV | 100 |
| `alpha_dna_opt4.mac` | alpha | `dna_opt4` | 100 eV – 400 MeV | 100 |
| `alpha_dna_opt6.mac` | alpha | `dna_opt6` | 100 eV – 400 MeV | 100 |
| `alpha_dna_opt8.mac` | alpha | `dna_opt8` | 100 eV – 400 MeV | 100 |

`validation/parser.py` publishes one plot per particle (x: energy in MeV, y:
stopping power in keV/µm, error: standard deviation / √N) with the **ICRU 90**
stopping power at the simulated energies it covers (from 1 keV), converted with
the `G4_WATER` density of 1 g/cm³.

## 📚 Reference: ICRU 90

`validation/icru90.py` returns the tabulated ICRU90 range values.

## Import notes

- `main.cc`, `include/`, `src/` are `spower.cc`, `include/` and `src/` of the
  Geant4 v11.3.2 example (the version of the default base image), unchanged.
- `macro/spower.in` is the example's macro, kept for reference; the validation
  macros are the `*_dna_opt2.mac` files listed in `validation/config.json`.
  `docs/geant4-example-README.txt` is the example's README.
- `macro/unit.mac` runs one 1 keV electron on every push.

For the CI, the validation flow and running everything locally, see the
[template README](https://github.com/G4Med-test/template#readme).

References for this example: Nucl. Instrum. Meth. B 397 (2017) 45–50, and the
Geant4-DNA publications Med. Phys. 51 (2024) 5873–5889; Med. Phys. 45 (2018)
e722–e739; Phys. Med. 31 (2015) 861–874; Med. Phys. 37 (2010) 4692–4708;
Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178.
