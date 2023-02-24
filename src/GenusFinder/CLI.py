import logging
from pathlib import Path
import subprocess as sp
import sys


class CLI:
    """
    A template wrapper class for CLI interactions
    """

    def __init__(self) -> None:
        self.args = []

    def _call(self):
        try:
            sp.run(self.args, check=True)
        except sp.CalledProcessError as e:
            logging.error(f"{' '.join(e.cmd)} returned code {e.returncode}")
            sys.exit()
        except sp.TimeoutExpired as e:
            logging.error(f"{' '.join(e.cmd)} timed out (timeout: {e.timeout})")
            sys.exit()


class MuscleAligner(CLI):
    def __init__(self):
        super().__init__()
        self.args = ["muscle"]
    
    def call_simple(self, align: Path, output: Path):
        self.args += [
            "-align",
            align,
            "-output",
            output
        ]
        self._call()

    def call_profile(self, profile: bool, in1: Path, in2: Path, out: Path):
        if profile:
            self.args.append("-profile")
        self.args += ["-in1", in1, "-in2", in2, "-out", out]

        self._call()


class RAxMLTreeBuilder(CLI):
    def __init__(self) -> None:
        super().__init__()
        self.args = ["raxmlHPC"]

    def call(self, b: int, f: str, N: int, m: str, n: str, p: int, s: Path, t: Path, w: Path, z: Path):
        self.args += ["-b", str(b)] if b else None
        self.args += ["-f", f] if f else None
        self.args += ["-N", str(N)] if N else None
        self.args += ["-m", m] if m else None
        self.args += ["-n", n] if n else None
        self.args += ["-p", p] if p else None
        self.args += ["-s", s] if s else None
        self.args += ["-t", t] if t else None
        self.args += ["-w", w.resolve()] if w else None
        self.args += ["-z", z] if z else None

        self._call()


class VsearchSearcher(CLI):
    def __init__(self) -> None:
        super().__init__()
        self.args = ["vsearch"]
    
    def call(self, u: Path, db: Path, id: float, fp: Path):
        self.args += [
            "--usearch_global",
            u,
            "--db",
            db,
            "--id",
            str(id),
            "--fastapairs",
            fp,
        ]
        self._call()