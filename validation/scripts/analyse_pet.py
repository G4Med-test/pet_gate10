#!/usr/bin/env python3
"""PET analysis of pet_gate10: analyse_pet.py <pet.root>

Port of pet_imaging's validation/scripts/analyse_pet.py to the opengate output:
it writes the same observables, in the same .dat files (r00a ... r12b) with the
same titles, so that the plots of the two tests can be compared. The GATE 9
branches map to the opengate ones as follows:

    runID, eventID1           RunID1, EventID1 (+ ThreadID1: event numbers
                              restart in every thread)
    globalPosX1 ...           PostPosition1_X ... (mm in both)
    time1                     GlobalTime1, in ns instead of s
    sourceID1                 EventPosition1_Y: the F18 line source (sourceID 0)
                              is at y = -2 cm, the O15 one (sourceID 1) at +2 cm
    comptonPhantom1           ComptonPhantom1 + ComptonEndshielding1 + ComptonSepta1
    RayleighPhantom1          the same with Rayleigh

The coincidence sorters of GATE 9 also reject pairs of adjacent heads (default
minSectorDifference 2), which opengate cannot do: those pairs are removed here,
from both the prompt and the delayed coincidences, before anything else.

stat.json (opengate SimulationStatisticsActor) next to the ROOT file gives the
number of decays for the sensitivity. The plots of the GATE 9 script are not
drawn: the portal only uses the .dat files.
"""
import json
import os
import sys

import numpy as np
import scipy.optimize
import uproot

HEADS = 8
MIN_SECTOR_DIFFERENCE = 2
NS_PER_S = 1e9
N_TRANSAXIAL = 1000  # restrict the transaxial plots to the first n coincidences


def head_index(volume_ids):
    """Head copy number from PreStepUniqueVolumeID, "<volume>-<copy numbers>".

    The copy numbers go from the world down: world, cylindricalPET, head, ...
    """
    heads = np.array([int(v.rsplit("-", 1)[1].split("_")[2]) for v in volume_ids],
                     dtype=int)
    if len(heads) and (heads.min() < 0 or heads.max() >= HEADS):
        sys.exit("Unexpected head copy numbers in PreStepUniqueVolumeID: %s"
                 % sorted(set(heads)))
    return heads


def read_coincidences(tree):
    """Branches of a coincidence tree, without the pairs of adjacent heads."""
    c = tree.arrays(library="np")
    if len(c["GlobalTime1"]) == 0:
        return c
    difference = np.abs(head_index(c["PreStepUniqueVolumeID1"])
                        - head_index(c["PreStepUniqueVolumeID2"])) % HEADS
    difference = np.minimum(difference, HEADS - difference)
    keep = difference >= MIN_SECTOR_DIFFERENCE
    print("%s: %d pairs, %d of adjacent heads removed"
          % (tree.name, len(keep), np.count_nonzero(~keep)))
    return {k: v[keep] for k, v in c.items()}


def scatter_count(c, n):
    return sum(c["%s%s%d" % (process, volume, n)]
               for process in ("Compton", "Rayleigh")
               for volume in ("Phantom", "Endshielding", "Septa"))


def write_scatter(filename, title, x_axis, y_axis, x, y):
    with open(filename, "w") as f:
        print("# Title: " + title, file=f)
        print("# X axis: " + x_axis, file=f)
        print("# Y axis: " + y_axis, file=f)
        for a, b in zip(x, y):
            print("{0};{1}".format(a, b), file=f)


def write_histogram(filename, title, x_axis, y_axis, counts, edges):
    with open(filename, "w") as f:
        print("# Title: " + title, file=f)
        print("# X axis: " + x_axis, file=f)
        print("# Y axis: " + y_axis, file=f)
        for low, high, n in zip(edges[:-1], edges[1:], counts):
            print("{0};{1};{2}".format(low, high, n), file=f)


def number_of_events(directory):
    try:
        with open(os.path.join(directory, "stat.json")) as f:
            return json.load(f)["events"]["value"]
    except (OSError, KeyError, ValueError) as e:
        print("No number of events from stat.json:", e)
        return None


def analyse_pet(filename):
    print("Filename", filename)
    f = uproot.open(filename)
    n_singles = f["Singles"].num_entries
    coinc = read_coincidences(f["Coincidences"])
    delays = read_coincidences(f["delay"])
    print("nb of singles", n_singles)
    print("nb of coincidences", len(coinc["GlobalTime1"]))
    print("nb of delays", len(delays["GlobalTime1"]))

    # Transaxial detection positions of the first coincidences of each run. The
    # x1, x2, x1 sequence is the one of the GATE 9 script.
    run = coinc["RunID1"]
    for run_id, filename_dat in ((0, "r00a.dat"), (1, "r00b.dat")):
        mask = run == run_id
        x1 = coinc["PostPosition1_X"][mask][:N_TRANSAXIAL]
        x2 = coinc["PostPosition2_X"][mask][:N_TRANSAXIAL]
        y1 = coinc["PostPosition1_Y"][mask][:N_TRANSAXIAL]
        y2 = coinc["PostPosition2_Y"][mask][:N_TRANSAXIAL]
        write_scatter(filename_dat,
                      "transaxial detection position (runid=={})".format(run_id),
                      "X", "Y", np.concatenate((x1, x2, x1)),
                      np.concatenate((y1, y2, y1)))

    # Axial detection
    z1 = coinc["PostPosition1_Z"]
    z2 = coinc["PostPosition2_Z"]
    counts, edges = np.histogram(np.concatenate((z1, z2)), bins=100)
    write_histogram("r01.dat", "axial coincidences detection position",
                    "Z position", "N", counts, edges)

    # True unscattered (tuc) and scattered (tsc) coincidences
    z = (z1 + z2) / 2
    unscattered = (scatter_count(coinc, 1) == 0) & (scatter_count(coinc, 2) == 0)
    tuc = z[unscattered]
    tsc = z[~unscattered]
    print("\tscattered", len(tsc))
    print("\tunscattered", len(tuc))
    counts, edges = np.histogram(tuc, bins=100)
    write_histogram("r02.dat", "axial sensitivity detection",
                    "Z position", "N", counts, edges)

    # Same computation as the GATE 9 script (two independent binnings)
    counts_scattered, bins_scattered = np.histogram(tsc, bins=100)
    counts_all, _ = np.histogram(z, bins=100)
    with np.errstate(divide="ignore", invalid="ignore"):
        fraction = counts_scattered.astype(float) / counts_all
    counts, edges = np.histogram(bins_scattered[:-1], bins=100, weights=fraction)
    write_histogram("r10.dat", "axial scatter fraction",
                    "Z position", "N", counts, edges)

    # Radioactive decays
    time = coinc["GlobalTime1"] / NS_PER_S
    y_source1 = coinc["EventPosition1_Y"]
    y_source2 = coinc["EventPosition2_Y"]
    decay_f18 = time[(y_source1 < 0) & (y_source2 < 0)]
    decay_o15 = time[(y_source1 > 0) & (y_source2 > 0)]

    heights, borders = np.histogram(decay_o15, bins="auto", density=True)
    centres = borders[:-1] + np.diff(borders) / 2

    def exponential(x, a, b):
        return a * np.exp(-b * x)

    popt, _ = scipy.optimize.curve_fit(exponential, centres, heights,
                                       p0=(heights[0], np.log(2) / 122.24))
    xx = np.linspace(0, 240, 240)
    yy = exponential(xx, *popt)
    half_life = np.log(2) / popt[1]

    counts, edges = np.histogram(decay_o15, bins=100, density=True)
    write_histogram("r11a.dat", "rad decays, O15 HL = 122.24 sec",
                    "time", "N", counts, edges)
    counts, edges = np.histogram(decay_f18, bins=100, density=True)
    write_histogram("r11b.dat", "rad decays, F18 HL = 6586.2 sec",
                    "time", "N", counts, edges)
    write_scatter("r11c.dat", "rad decays, O15 fit HL = {:.2f} sec".format(half_life),
                  "time", "N", xx, yy)

    # Randoms: the two singles come from different decays
    different = ((coinc["EventID1"] != coinc["EventID2"])
                 | (coinc["ThreadID1"] != coinc["ThreadID2"]))
    randoms = time[different]
    delay_time = delays["GlobalTime1"] / NS_PER_S
    print("nb of randoms", len(randoms))
    print("nb of delays", len(delay_time))
    counts, edges = np.histogram(randoms, bins=100)
    write_histogram("r12a.dat", "randoms, random = {}".format(len(randoms)),
                    "time", "N", counts, edges)
    counts, edges = np.histogram(delay_time, bins=100)
    write_histogram("r12b.dat",
                    "randoms, delays with coinc sorter = {}".format(len(delay_time)),
                    "time", "N", counts, edges)

    n_events = number_of_events(os.path.dirname(os.path.abspath(filename)))
    print("Number of events", n_events)
    print("Number of singles", n_singles)
    print("Number of coincidences", len(time))
    print("Number of true", len(tuc))
    print("Number of randoms", len(randoms))
    print("Number of scatter", len(tsc))
    if n_events:
        print("Absolute sensitivity {:.2f} %".format(len(tuc) / n_events * 100.0))
    print("O15 fit half-life {:.2f} s".format(half_life))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: analyse_pet.py <root file>")
    analyse_pet(sys.argv[1])
