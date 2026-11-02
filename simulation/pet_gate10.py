#!/usr/bin/env python3
"""GATE 10 (opengate) simulation of the pet_gate10 test.

Usage: pet_gate10.py <macro>

The setup is the one of the GATE 9 test pet_imaging (macro/main.mac and
macro/mac/*.mac there), translated to opengate: a cylindrical PET of 8 heads with
LSO/BGO phoswich crystals, lead end shielding and tungsten septa, orbiting around
the axis; a water cylinder with an 18F and an 15O positron line source; the same
digitizer and coincidence sorters. The GATE 9 command is quoted next to each
setting; README.md lists where the two versions differ.

The pipeline starts the container as `<image> <macro>`, so every input parameter
is in the macro: a parameter file with one `/command value` per line, in the
syntax of a Geant4 macro (comments start with #). Every command of COMMANDS
must appear exactly once; an unknown command stops the run with
"COMMAND NOT FOUND", which the pipeline treats as a failed run.

Output, in ./output like in pet_imaging: pet.root (trees Singles, Coincidences
and delay) and stat.json (simulation statistics).
"""
import sys
from importlib.metadata import version
from pathlib import Path

MATERIALS = Path(__file__).resolve().parent / "GateMaterials.db"
OUTPUT_DIR = "output"

# Interactions in the volumes that pet_imaging attaches to the phantom sensitive
# detector (GATE 9 comptonPhantom and RayleighPhantom): one counter attribute
# per process and volume, summed by the analysis.
SCATTER_PROCESSES = {"Compton": "compt", "Rayleigh": "Rayl"}
SCATTER_VOLUMES = ("phantom", "endshielding", "septa")


def text(value):
    return value


def integer(value):
    return int(value)


def seed(value):
    return value if value == "auto" else int(value)


def quantity(value):
    """'10 cm' -> 10 * cm in Geant4 internal units."""
    from opengate.utility import g4_units

    number, unit = value.split()
    if unit not in g4_units:
        raise ValueError("Unknown unit: " + unit)
    return float(number) * g4_units[unit]


COMMANDS = {
    "/pet/physicsList": text,
    "/pet/productionCut": quantity,
    "/pet/randomSeed": seed,
    "/pet/numberOfThreads": integer,
    "/pet/timeSlice": quantity,
    "/pet/timeStop": quantity,
}


def read_macro(path):
    """Return {command: raw value}, checking every command of COMMANDS once."""
    values = {}
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        command, _, value = line.partition(" ")
        if command not in COMMANDS:
            sys.exit("***** COMMAND NOT FOUND <%s> *****" % line)
        if command in values:
            sys.exit("Repeated command in the macro: " + command)
        values[command] = value.strip()
    missing = [c for c in COMMANDS if c not in values]
    if missing:
        sys.exit("Missing commands in the macro: " + ", ".join(missing))
    return values


def grid(numbers, steps, centre=(0, 0, 0)):
    """Translations of a GATE cubicArray repeater (centred on the original position)."""
    from opengate.geometry.utility import get_grid_repetition

    return [[t[i] + centre[i] for i in range(3)]
            for t in get_grid_repetition(numbers, steps)]


def run_intervals(time_slice, time_stop):
    """GATE 9 setTimeStart 0, setTimeSlice, setTimeStop: one run per slice."""
    intervals = []
    start = 0
    while start < time_stop * (1 - 1e-9):
        intervals.append([start, min(start + time_slice, time_stop)])
        start += time_slice
    return intervals


def add_scanner(sim, intervals):
    """mac/pet_head.mac: "NOT a true PET device", designed to check several features."""
    from opengate.geometry.utility import get_circular_repetition
    from opengate.utility import g4_units
    from scipy.spatial.transform import Rotation

    mm, cm, deg, s = g4_units.mm, g4_units.cm, g4_units.deg, g4_units.s

    pet = sim.add_volume("Tubs", "cylindricalPET")
    pet.rmin = 39.9 * cm
    pet.rmax = 52.0 * cm
    pet.dz = 40.2 * cm / 2
    pet.material = "Air"
    # /gate/cylindricalPET/orbiting/setSpeed .1875 deg/s around the z axis: GATE 9
    # moves the geometry at the start of every time slice.
    speed = 0.1875 * deg / s
    pet.add_dynamic_parametrisation(
        translation=[[0, 0, 0] for _ in intervals],
        rotation=[Rotation.from_rotvec([0, 0, speed * start]).as_matrix()
                  for start, _ in intervals])

    # rsector: ring repeater of 8
    head = sim.add_volume("Box", "head")
    head.mother = pet.name
    head.size = [8.0 * cm, 32.0 * cm, 40.2 * cm]
    head.material = "Air"
    head.translation, head.rotation = get_circular_repetition(8, [44.0 * cm, 0, 0])

    endshielding = sim.add_volume("Box", "endshielding")
    endshielding.mother = head.name
    endshielding.size = [8.0 * cm, 32.0 * cm, 1.0 * mm]
    endshielding.material = "Lead"
    endshielding.translation = grid([1, 1, 2], [0, 0, 40.1 * cm])

    septa = sim.add_volume("Box", "septa")
    septa.mother = head.name
    septa.size = [5.0 * cm, 32.0 * cm, 0.5 * mm]
    septa.material = "Tungsten"
    septa.translation = grid([1, 1, 3], [0, 0, 10.0 * cm], [-1.5 * cm, 0, 0])

    module = sim.add_volume("Box", "module")
    module.mother = head.name
    module.size = [3.0 * cm, 8.0 * cm, 10.0 * cm]
    module.material = "Air"
    module.translation = grid([1, 4, 4], [0, 8.0 * cm, 10.0 * cm], [2.5 * cm, 0, 0])

    block = sim.add_volume("Box", "block")
    block.mother = module.name
    block.size = [30 * mm, 15.9 * mm, 19.9 * mm]
    block.material = "Air"
    block.translation = grid([1, 5, 5], [0, 1.6 * cm, 2.0 * cm])

    crystal = sim.add_volume("Box", "crystal")
    crystal.mother = block.name
    crystal.size = [3.0 * cm, 3.0 * mm, 3.8 * mm]
    crystal.material = "Air"
    crystal.translation = grid([1, 5, 5], [0, 3.2 * mm, 4.0 * mm])

    layers = []
    for name, x in (("LSO", -0.75 * cm), ("BGO", 0.75 * cm)):
        layer = sim.add_volume("Box", name)
        layer.mother = crystal.name
        layer.size = [1.5 * cm, 3.0 * mm, 3.8 * mm]
        layer.material = name
        layer.translation = [x, 0, 0]
        layers.append(layer.name)
    return head, layers


def add_digitizer(sim, head, layers):
    """mac/pet_digitizer.mac (GATE >= 9.3 syntax) and mac/output.mac."""
    from opengate.utility import g4_units

    keV, ns = g4_units.keV, g4_units.ns

    def digitizer(kind, name, source):
        actor = sim.add_actor(kind, name)
        actor.attached_to = layers
        actor.authorize_repeated_volumes = True
        actor.input_digi_collection = source.name
        actor.root_output.write_to_disk = False
        return actor

    scatter = []
    for label, process in SCATTER_PROCESSES.items():
        for volume in SCATTER_VOLUMES:
            name = label + volume.capitalize()
            attribute = sim.activate_auxiliary_attribute(
                "ProcessDefinedStepInVolumeAttribute", name)
            attribute.process_name = process
            attribute.volume_name = volume
            attribute.propagate_from_parent_track = True
            scatter.append(name)

    # attachCrystalSD on both layers. EventPosition (the decay point) tells the
    # two sources apart, like the GATE 9 sourceID.
    hits = sim.add_actor("DigitizerHitsCollectionActor", "Hits")
    hits.attached_to = layers
    hits.authorize_repeated_volumes = True
    hits.attributes = [
        "EventID", "ThreadID", "RunID", "GlobalTime", "TotalEnergyDeposit",
        "PostPosition", "PreStepUniqueVolumeID", "EventPosition",
    ] + scatter
    hits.root_output.write_to_disk = False

    # adder in each layer (default policy energyWeightedCentroid); one collection
    # for both layers replaces the GATE 9 merger of Singles_LSO into Singles_BGO.
    adder = digitizer("DigitizerAdderActor", "Adder", hits)
    adder.policy = "EnergyWeightedCentroidPosition"

    # readout setDepth 1: the rsector, i.e. the head (default policy TakeEnergyWinner)
    readout = digitizer("DigitizerAdderActor", "Readout", adder)
    readout.policy = "EnergyWinnerPosition"
    readout.group_volume = head.name

    # energyResolution fwhm 0.26 at 511 keV (inverse square law)
    blurring = digitizer("DigitizerBlurringActor", "EnergyResolution", readout)
    blurring.blur_attribute = "TotalEnergyDeposit"
    blurring.blur_method = "InverseSquare"
    blurring.blur_resolution = 0.26
    blurring.blur_reference_value = 511 * keV

    # energyFraming 350-650 keV: its only channel is the Singles tree
    window = digitizer("DigitizerEnergyWindowsActor", "EnergyFraming", blurring)
    window.channels = [{"name": "Singles", "min": 350 * keV, "max": 650 * keV}]
    window.output_filename = "pet.root"

    # Prompt and delayed coincidence sorters on the singles. As in GATE 9, a single
    # cannot open a window while another one is open (multi_window False), and two
    # singles of the same head never pair (group_volume). The GATE 9 default
    # minSectorDifference 2 (no pairs of adjacent heads either) has no opengate
    # option: the analysis applies it.
    for name, offset, policy in (("Coincidences", 0, "TakeWinnerOfGoods"),
                                 ("delay", 500 * ns, "TakeWinnerIfAllAreGoods")):
        sorter = sim.add_actor("CoincidenceSorterActor", name)
        sorter.input_digi_collection = "Singles"
        sorter.window = 120 * ns
        sorter.offset = offset
        sorter.multiples_policy = policy
        sorter.multi_window = False
        sorter.group_volume = head.name
        sorter.output_filename = "pet.root"


def build_simulation(p):
    from opengate.managers import Simulation
    from opengate.utility import g4_units

    mm, cm, m, s, Bq = (g4_units.mm, g4_units.cm, g4_units.m, g4_units.s,
                        g4_units.Bq)
    threads = p["/pet/numberOfThreads"]

    sim = Simulation()
    sim.g4_verbose = False
    sim.visu = False
    sim.random_engine = "MersenneTwister"
    sim.random_seed = p["/pet/randomSeed"]
    sim.number_of_threads = threads
    sim.output_dir = Path.cwd() / OUTPUT_DIR
    intervals = run_intervals(p["/pet/timeSlice"], p["/pet/timeStop"])
    sim.run_timing_intervals = intervals

    # mac/world.mac
    sim.volume_manager.add_material_database(str(MATERIALS))
    sim.world.size = [1.5 * m] * 3
    sim.world.material = "Air"

    head, layers = add_scanner(sim, intervals)

    # mac/cylindrical_phantom.mac
    phantom = sim.add_volume("Tubs", "phantom")
    phantom.rmin = 0
    phantom.rmax = 10.0 * cm
    phantom.dz = 70.0 * cm / 2
    phantom.material = "Water"

    # mac/physics.mac: 1 km everywhere, the macro's cut in phantom and crystals
    physics = sim.physics_manager
    physics.physics_list_name = p["/pet/physicsList"]
    particles = ("gamma", "electron", "positron")
    for particle in particles:
        physics.global_production_cuts[particle] = 1 * g4_units.km
    for volume in [phantom.name] + layers:
        for particle in particles:
            physics.set_production_cut(volume, particle, p["/pet/productionCut"])

    # mac/sources.mac: two line sources of 100 kBq at t = 0 with forced half-lives.
    # Every thread runs its own copy of a source, hence the activity per thread.
    for name, spectrum, half_life, y in (
            ("F18LineSource", "F18_analytic", 6586.2 * s, -2.0 * cm),
            ("O15LineSource", "O15_analytic", 122.24 * s, 2.0 * cm)):
        source = sim.add_source("GenericSource", name)
        source.particle = "e+"
        source.energy.type = spectrum
        source.activity = 100000.0 * Bq / threads
        source.half_life = half_life
        source.position.type = "cylinder"
        source.position.radius = 0.5 * mm
        source.position.dz = 34.0 * cm
        source.position.translation = [0, y, 0]
        source.direction.type = "iso"

    add_digitizer(sim, head, layers)

    stats = sim.add_actor("SimulationStatisticsActor", "stat")
    stats.track_types_flag = True
    stats.output_filename = "stat.json"
    return sim, stats


def summary():
    """Number of entries of the output trees, as a check in the log."""
    import uproot

    with uproot.open(Path.cwd() / OUTPUT_DIR / "pet.root") as f:
        for tree in ("Singles", "Coincidences", "delay"):
            print("%-13s %d" % (tree, f[tree].num_entries if tree in f else 0))
        if "Coincidences" in f and f["Coincidences"].num_entries:
            ids = f["Coincidences"]["PreStepUniqueVolumeID1"].array(
                library="np", entry_stop=3)
            print("PreStepUniqueVolumeID1 of the first coincidences:", list(ids))


def main(argv):
    if len(argv) != 2:
        sys.exit("Usage: %s <macro>" % argv[0])
    raw = read_macro(argv[1])

    import opengate_core as g4

    params = {command: COMMANDS[command](value) for command, value in raw.items()}
    print("opengate", version("opengate"))
    print("Geant4  ", g4.GateInfo.get_G4Version().replace("$Name:", "").strip(" $"))
    print("ITK     ", g4.GateInfo.get_ITKVersion())
    for command, value in raw.items():
        print("%-22s %s" % (command, value))

    sim, stats = build_simulation(params)
    sim.run()
    print(stats)
    summary()


if __name__ == "__main__":
    main(sys.argv)
