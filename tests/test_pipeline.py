"""Small CPU-only regression tests; no AMBER execution or production MD."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import MDAnalysis as mda
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pnsp.provenance import execution_provenance
from pnsp.run_md import protocol
from pnsp.analyse_occupancy import events


class PipelineTests(unittest.TestCase):
    def test_provenance(self):
        record = execution_provenance("test", "phenol", "rep1", [__file__], ROOT)
        self.assertEqual(record["contaminant_id"], "phenol")
        self.assertEqual(record["replica_id"], "rep1")
        self.assertTrue(record["recorded_at_utc"].endswith("+00:00"))
        self.assertEqual(record["input_paths"], [str(Path(__file__).resolve())])

    def test_protocol(self):
        stages = protocol(15, 303, 1)
        self.assertEqual(len(stages), 11)
        self.assertIn("nstlim=7500000, dt=0.002", stages[-1].mdin)
        self.assertIn("ntwx=5000", stages[-1].mdin)
        self.assertNotIn("restraintmask", stages[-1].mdin)

    def test_events(self):
        self.assertEqual(events(np.array([1, 1, 0, 1], dtype=bool)), [(0, 2), (3, 4)])
        self.assertEqual(events(np.zeros(4, dtype=bool)), [])

    def test_strict_cutoff_and_copy_aggregation(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            u = mda.Universe.empty(3, n_residues=3, atom_resindex=[0, 1, 2], trajectory=True)
            u.add_TopologyAttr("names", ["C", "C", "C"])
            u.add_TopologyAttr("resnames", ["F01", "LIG", "LIG"])
            u.add_TopologyAttr("resids", [1, 2, 3])
            u.add_TopologyAttr("masses", [12., 12., 12.])
            u.dimensions = [32, 32, 32, 90, 90, 90]
            u.atoms.positions = [[0, 0, 0], [3, 0, 0], [10, 0, 0]]
            u.atoms.write(str(work / "system.pdb"))
            with mda.Writer(str(work / "prod.nc"), n_atoms=3) as writer:
                for distance in (3, 3, 4, 3):
                    u.atoms.positions = [[0, 0, 0], [distance, 0, 0], [10, 0, 0]]
                    writer.write(u.atoms)
            subprocess.run([sys.executable, str(ROOT / "pnsp/analyse_occupancy.py"),
                            "--topology", str(work / "system.pdb"),
                            "--trajectory", str(work / "prod.nc"), "--output", str(work / "out"),
                            "--system", "test", "--replica", "rep1", "--expected-frames", "4"], check=True)
            copies = pd.read_csv(work / "out/copy_metrics.csv")
            self.assertEqual(copies.occupancy.tolist(), [.75, 0])
            self.assertEqual(copies.n_events.tolist(), [2, 0])
            self.assertAlmostEqual(copies.mean_event_duration_ns.iloc[0], .015)
            self.assertTrue(np.isnan(copies.mean_event_duration_ns.iloc[1]))
            manifest = json.loads((work / "out/analysis_manifest.json").read_text())
            self.assertEqual(manifest["provenance"]["replica_id"], "rep1")


if __name__ == "__main__":
    unittest.main()
