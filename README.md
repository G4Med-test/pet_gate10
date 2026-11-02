# pet_gate10

**In-vivo PET validation test built on GATE 10** ([opengate](https://github.com/OpenGATE/opengate)),
compiled from source against the Geant4 of the G4Med
[`geant4-alma9`](https://github.com/G4Med-test/geant4-alma9) image.

The simulation reproduces [`pet_imaging`](https://github.com/G4Med-test/pet_imaging),
the macro-based GATE 9 version of this test: same scanner, phantom, sources,
digitizer, acquisition and observables. This allows the two GATE versions to be compared on the
same Geant4.

---

## How the container is built

`Apptainer.def` starts from `geant4-alma9` and builds:

1. **ITK** (`ITK_TAG`, default `v5.4.4`), required by `opengate_core`. Only the
   modules it uses (`ITKCommon`, `ITKImageFilterBase`, `ITKImageIntensity`) are built,
   as static, position-independent libraries.
2. **`opengate_core`** (the C++/pybind11 part of GATE 10), from the
   `OPENGATE_TAG` sources (default `10.1.1`). It is compiled by its own `setup.py`
   against `/opt/geant4` and ITK. The `opengate-core` wheels on PyPI are **not**
   used, because they bundle their own Geant4.
3. The pure-Python **`opengate`** package and its dependencies (itk, uproot,
   pandas, ...) from PyPI. None of them contains Geant4, and the `opengate-core`
   requirement is already satisfied by step 2.

The `%test` section checks that `opengate_core` links only the Geant4 libraries of
`/opt/geant4` and reports the same Geant4 version as `geant4-config`.

The base Geant4 already has the options of the official opengate "novis" build:
multithreaded, `GEANT4_BUILD_TLS_MODEL=global-dynamic` (needed to load Geant4 from a
Python module), no GDML, no Qt or OpenGL.

### Geant4 version

opengate 10.1.1 targets **Geant4 11.4**: its CI and dataset list are those of 11.4.
That is why the default `TAG` (the base image's Geant4) is `v11.4.3`, unlike the
`v11.3.2` of most tests. Change `OPENGATE_TAG` and `TAG` together:

```bash
apptainer build --build-arg TAG=v11.4.3 --build-arg OPENGATE_TAG=10.1.1 pet_gate10.sif Apptainer.def
```

### Geant4 datasets

When imported, `opengate_core` normally sets `G4LEDATA` & co. to a `geant4_data`
folder inside the package and downloads the datasets there. The image sets
`GATEONRTD=1`, the switch opengate uses for its ReadTheDocs build, so this step is
skipped and Geant4 finds the datasets mounted at `/g4data` through
`GEANT4_DATA_DIR`, like every other test. Check that this switch still exists when
changing `OPENGATE_TAG`.

## The simulation

[`simulation/pet_gate10.py`](simulation/pet_gate10.py) translates the macros of
pet_imaging (`macro/main.mac` and `macro/mac/*.mac` there) to opengate. Each
setting quotes the GATE 9 command it comes from.

- **Scanner** (`pet_head.mac`, "NOT a true PET device", built to check several
  features): a cylindrical PET of 8 heads (rsectors) at 44 cm. Each head has
  4×4 modules of 5×5 blocks of 5×5 crystals. Every crystal is an LSO/BGO phoswich
  of two 1.5 cm layers. Each head also has lead end shielding and tungsten septa.
  The gantry orbits at 0.1875 deg/s around the axis.
- **Phantom**: a water cylinder of 10 cm radius and 70 cm length.
- **Sources** (`sources.mac`): ¹⁸F and ¹⁵O line sources of 100 kBq each,
  0.5 mm radius and 68 cm length, at y = −2 cm and +2 cm. They emit positrons with
  the GATE analytic spectra and decay with their half-lives.
- **Physics** (`physics.mac`): the physics list and production cut from the
  macro. The cut applies in the phantom and the crystals; it is 1 km elsewhere.
  Materials come from pet_imaging's `GateMaterials.db`.
- **Digitizer** (`pet_digitizer.mac`):
  1. adder in each layer;
  2. readout at the head;
  3. energy resolution of 26 % at 511 keV;
  4. 350–650 keV window, giving the `Singles` tree;
  5. prompt (`Coincidences`) and delayed (`delay`, 500 ns offset) coincidence
     sorters, with a 120 ns window.
- **Acquisition**: time slices of `/pet/timeSlice` up to `/pet/timeStop`. GATE 9
  moves the geometry at the start of each slice, so there is one run per slice.
  `main.mac` runs 2 × 120 s, like pet_imaging.
- **Output** in `output/`: `pet.root` (`Singles`, `Coincidences`, `delay`) and
  `stat.json`.

### Differences from GATE 9

| GATE 9 (pet_imaging) | GATE 10 (here) |
| --- | --- |
| merger of `Singles_LSO` into `Singles_BGO` | one hits collection on both layers |
| `sourceID` | `EventPosition` of the decay: the two sources are 4 cm apart |
| `comptonPhantom`, `RayleighPhantom` (PhantomSD on phantom, end shielding and septa) | `ProcessDefinedStepInVolumeAttribute` counters on the same volumes, propagated to secondaries and summed by the analysis |
| `minSectorDifference 2` (default): no pairs from the same or adjacent heads | the sorters reject pairs from the same head (`group_volume`); the analysis removes the pairs from adjacent heads |
| `delay` sorter policy: default `takeWinnerIfAllAreGoods` | `TakeWinnerIfAllAreGoods` |
| time in s, a single thread | time in ns; with several threads, the activity is split among them. Event numbers restart in every thread, so the analysis compares (`EventID`, `ThreadID`) pairs |

The one approximation is the adjacent-head cut. Applying it after the
`TakeWinnerOfGoods` policy, instead of inside it, can make a different choice
only when more than two singles share the 120 ns window. At 200 kBq this is
rare.

## Running

The pipeline starts the container as `<image> <macro>`. The runscript runs
[`simulation/pet_gate10.py`](simulation/pet_gate10.py) with that file, which is a
**parameter macro**: one `/command value` per line, in the syntax of a Geant4 macro,
so that `metadata()` in `validation/parser.py` reads it like any other macro. Every
command must appear exactly once. An unknown command stops the run with
`COMMAND NOT FOUND`.

| Command | Example | Meaning |
| --- | --- | --- |
| `/pet/physicsList` | `G4EmStandardPhysics_option3` | Geant4 physics list or EM constructor (GATE 9 `emstandard_opt3`) |
| `/pet/productionCut` | `1 mm` | production cut for γ, e⁻ and e⁺ in the phantom and crystals |
| `/pet/randomSeed` | `auto` | MersenneTwister seed, or `auto` |
| `/pet/numberOfThreads` | `1` | Geant4 threads |
| `/pet/timeSlice` | `120 s` | length of a time slice (one run each) |
| `/pet/timeStop` | `240 s` | end of the acquisition |

`macro/main.mac` is the validation run (`validation/config.json`), as in pet_imaging.
`macro/unit.mac` is the unit test run on every push. It is `main.mac` with two
1 s slices, about 400 000 decays. Its log ends with the number of entries of
the output trees.

## Validation

`validation/parser.py` runs [`validation/scripts/analyse_pet.py`](validation/scripts/analyse_pet.py)
on `output/pet.root`. The script is a port of pet_imaging's analysis, and the
portal plots keep the same names:

- transaxial detection positions, for both runs;
- axial distribution of the detected coincidences;
- axial sensitivity (unscattered coincidences);
- axial scatter fraction;
- ¹⁵O and ¹⁸F decay time distributions, with an exponential fit of the ¹⁵O one;
- randoms and delayed coincidences in time.

The script header gives the branch mapping between the two GATE versions.

```bash
apptainer build pet_gate10.sif Apptainer.def
apptainer run -B /path/to/geant4/data:/g4data pet_gate10.sif macro/unit.mac
```
