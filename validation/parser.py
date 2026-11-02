"""Turns the output of one validation run into plots for the Geant Validation Portal.

Not used yet: there is no validation/config.json, so the pipeline skips the
validation until the observables and reference data of this test are defined
with its referents. metadata() already reads the parameter macro of
simulation/pet_gate10.py.
"""
from geantval import one_command


def metadata(commands):
    """Describe one run, reading every setting from its macro."""
    return {
        "TEST": "pet_gate10",
        "PHYSICS_LIST": one_command(commands, "/pet/physicsList"),
        "PARTICLE": one_command(commands, "/pet/source/particle"),
        "SPECTRUM": one_command(commands, "/pet/source/spectrum"),
        "EVENTS": int(one_command(commands, "/pet/source/events")),
        "STATISTICS": one_command(commands, "/pet/output/statistics"),
    }


def parse(job):
    """Read the output of the run in job["path"] and yield one plot per `yield`."""
    raise NotImplementedError("pet_gate10: observables still to be defined")
    yield  # pragma: no cover
