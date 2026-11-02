# pet_gate10

**In-vivo PET validation test built on GATE 10** ([opengate](https://github.com/OpenGATE/opengate)),
compiled from source against the Geant4 of the G4Med
[`geant4-alma9`](https://github.com/G4Med-test/geant4-alma9) image.

> **Status: skeleton.** The container and the unit test work, but the simulation is
> only a placeholder: a water cylinder with a positron (¹⁸F spectrum) line source and
> simulation statistics. The PET geometry, digitizer, observables and reference data
> are still to be defined with the test referents. Until then there is no
> `validation/config.json`, so the pipeline skips the validation.
> The macro-based GATE 9 version of this test is
> [`pet_imaging`](https://github.com/G4Med-test/pet_imaging).

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

## Running

The pipeline starts the container as `<image> <macro>`. The runscript runs
[`simulation/pet_gate10.py`](simulation/pet_gate10.py) with that file, which is a
**parameter macro**: one `/command value` per line, in the syntax of a Geant4 macro,
so that `metadata()` in `validation/parser.py` reads it like any other macro. Every
command must appear exactly once. An unknown command stops the run with
`COMMAND NOT FOUND`.

| Command | Example | Meaning |
| --- | --- | --- |
| `/pet/physicsList` | `G4EmStandardPhysics_option3` | Geant4 physics list or EM constructor |
| `/pet/productionCut` | `1 mm` | production cut for γ, e⁻, e⁺ and p |
| `/pet/randomSeed` | `123456789` | MersenneTwister seed |
| `/pet/numberOfThreads` | `1` | Geant4 threads |
| `/pet/phantom/material` | `G4_WATER` | phantom material |
| `/pet/phantom/radius` | `10 cm` | phantom cylinder radius |
| `/pet/phantom/halfLength` | `10 cm` | phantom cylinder half-length |
| `/pet/source/particle` | `e+` | source particle |
| `/pet/source/spectrum` | `F18` | opengate energy spectrum type |
| `/pet/source/radius` | `0.5 mm` | line source radius |
| `/pet/source/halfLength` | `5 cm` | line source half-length |
| `/pet/source/events` | `100` | number of primaries |
| `/pet/output/statistics` | `statistics.txt` | simulation statistics file (in the working directory) |

`macro/unit.mac` (100 positrons) is the unit test run on every push.

```bash
apptainer build pet_gate10.sif Apptainer.def
apptainer run -B /path/to/geant4/data:/g4data pet_gate10.sif macro/unit.mac
```
