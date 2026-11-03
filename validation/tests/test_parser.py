"""parse() on .dat files in the format of scripts/analyse_pet.py, without running it.

Run by the self-test of ci-workflows, with its validation/ directory (geantval,
export) in PYTHONPATH and in G4MED_TOOLS.
"""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import parser as pet_parser
import export

SCATTER = """# Title: transaxial detection position (runid==0)
# X axis: X
# Y axis: Y
1.5;-2.0
3.0;4.5
"""
HISTOGRAM = """# Title: axial coincidences detection position
# X axis: Z position
# Y axis: N
-10.0;0.0;3
0.0;10.0;5
"""


class ParseTests(unittest.TestCase):
    def test_records_pass_the_checks_of_the_exporter(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "output"
            output.mkdir()
            (output / "r00a.dat").write_text(SCATTER)
            (output / "r01.dat").write_text(HISTOGRAM)
            job = {"TEST": "pet_gate10", "PHYSICS_LIST": "G4EmStandardPhysics_option3",
                   "path": tmp, "VERSION": "11.4.3"}
            with mock.patch.object(pet_parser.subprocess, "run") as run:
                records = list(pet_parser.parse(job))
        self.assertEqual(run.call_args.kwargs["cwd"], str(output))
        self.assertEqual([r["plotType"] for r in records], ["SCATTER2D", "TH1"])
        for record in records:
            export.validate_record(record)
        self.assertEqual(records[0]["chart"]["xValues"], [1.5, 3.0])
        self.assertEqual(records[1]["histogram"]["binContent"], [3.0, 5.0])
        self.assertEqual(records[1]["metadata"]["observableName"],
                         "axial coincidences detection position")


if __name__ == "__main__":
    unittest.main()
