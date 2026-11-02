"""Turns the output of one validation run into plots for the Geant Validation Portal.

metadata() reads the parameter macro of simulation/pet_gate10.py; parse() runs
validation/scripts/analyse_pet.py on the run's output/pet.root and converts each
.dat file it writes into a scatter or histogram record. The observables, their
.dat format and the record fields are those of pet_imaging, the GATE 9 version of
this test, so that the two can be compared on the portal.
"""
import os
import subprocess
import sys

from geantval import getJSON, one_command

PARTICLE = "Fluor18"
ENERGY_GEV = 0.512
DETECTOR = "PET"


def metadata(commands):
    """Describe one run, reading every setting from its macro."""
    # Like GATE 9, the acquisition is set by time, not by a number of events.
    one_command(commands, "/pet/timeStop")
    return {"TEST": "pet_gate10",
            "PHYSICS_LIST": one_command(commands, "/pet/physicsList")}


def extract_observable(filename):
    x, y, z, name = [], [], [], []
    with open(filename) as f:
        for line in f:
            if "#" in line:
                label = line.split(":", 1)[1].strip()
                label = label.split(" = ")[0].strip()
                name.append(label)
                continue
            fields = line.split(";")
            x.append(fields[0].strip())
            y.append(fields[1].strip())
            if len(fields) == 3:
                z.append(fields[2].strip())
    return x, y, z, name


def parse(job):
    """Read the output of the run in job["path"] and yield one plot per `yield`."""
    output_dir = os.path.join(job["path"], "output")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    subprocess.run([sys.executable, os.path.join(script_dir, "scripts", "analyse_pet.py"), "pet.root"],
                   cwd=output_dir, check=True)

    for filename in sorted(os.listdir(output_dir)):
        if not filename.endswith(".dat"):
            continue
        x, y, z, name = extract_observable(os.path.join(output_dir, filename))
        if len(z) == 0:
            yield getJSON(job, "chart",
                          mctool_name="GEANT4",
                          mctool_model=job["PHYSICS_LIST"],
                          plotType="SCATTER2D",
                          observableName=name[0],
                          targetName=DETECTOR,
                          beamParticle=PARTICLE,
                          beamEnergies=[ENERGY_GEV],
                          title=name[0],
                          secondaryParticle="None",
                          xAxisName=name[1],
                          yAxisName=name[2],
                          xValues=[float(v) for v in x],
                          yValues=[float(v) for v in y])
        else:
            yield getJSON(job, "histogram",
                          mctool_name="GEANT4",
                          mctool_model=job["PHYSICS_LIST"],
                          observableName=name[0],
                          targetName=DETECTOR,
                          beamParticle=PARTICLE,
                          beamEnergies=[ENERGY_GEV],
                          secondaryParticle="None",
                          title=name[0],
                          xAxisName=name[1],
                          yAxisName=name[2],
                          binEdgeLow=[float(v) for v in x],
                          binEdgeHigh=[float(v) for v in y],
                          binContent=[float(v) for v in z])
