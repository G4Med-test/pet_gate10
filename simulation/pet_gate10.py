#!/usr/bin/env python3
"""GATE 10 (opengate) simulation of the pet_gate10 test.

Usage: pet_gate10.py <macro>

The pipeline starts the container as `<image> <macro>`, so every input parameter
is in the macro: a parameter file with one `/command value` per line, in the
syntax of a Geant4 macro (comments start with #). Every command of COMMANDS
must appear exactly once; an unknown command stops the run with
"COMMAND NOT FOUND", which the pipeline treats as a failed run.

The geometry is a first skeleton (a cylindrical phantom with a positron line
source and simulation statistics), to be replaced by the PET setup defined with
the test referents.
"""
import sys
from importlib.metadata import version
from pathlib import Path


def text(value):
    return value


def integer(value):
    return int(value)


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
    "/pet/randomSeed": integer,
    "/pet/numberOfThreads": integer,
    "/pet/phantom/material": text,
    "/pet/phantom/radius": quantity,
    "/pet/phantom/halfLength": quantity,
    "/pet/source/particle": text,
    "/pet/source/spectrum": text,
    "/pet/source/radius": quantity,
    "/pet/source/halfLength": quantity,
    "/pet/source/events": integer,
    "/pet/output/statistics": text,
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


def build_simulation(p):
    from opengate.managers import Simulation
    from opengate.utility import g4_units

    sim = Simulation()
    sim.g4_verbose = False
    sim.visu = False
    sim.random_engine = "MersenneTwister"
    sim.random_seed = p["/pet/randomSeed"]
    sim.number_of_threads = p["/pet/numberOfThreads"]
    sim.output_dir = Path.cwd()

    sim.world.size = [1 * g4_units.m] * 3
    sim.world.material = "G4_AIR"

    phantom = sim.add_volume("Tubs", "phantom")
    phantom.rmin = 0
    phantom.rmax = p["/pet/phantom/radius"]
    phantom.dz = p["/pet/phantom/halfLength"]
    phantom.material = p["/pet/phantom/material"]

    sim.physics_manager.physics_list_name = p["/pet/physicsList"]
    cuts = sim.physics_manager.global_production_cuts
    cuts.gamma = cuts.electron = cuts.positron = cuts.proton = p["/pet/productionCut"]

    source = sim.add_source("GenericSource", "source")
    source.particle = p["/pet/source/particle"]
    source.energy.type = p["/pet/source/spectrum"]
    source.position.type = "cylinder"
    source.position.radius = p["/pet/source/radius"]
    source.position.dz = p["/pet/source/halfLength"]
    source.direction.type = "iso"
    source.n = p["/pet/source/events"]

    stats = sim.add_actor("SimulationStatisticsActor", "statistics")
    stats.track_types_flag = True
    stats.output_filename = p["/pet/output/statistics"]
    return sim, stats


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
        print("%-26s %s" % (command, value))

    sim, stats = build_simulation(params)
    sim.run()
    print(stats)


if __name__ == "__main__":
    main(sys.argv)
