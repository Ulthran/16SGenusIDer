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

    def call(self, profile: bool, in1: Path, in2: Path, out: Path):
        if profile:
            self.args.append("-profile")
        self.args += ["-in1", in1, "-in2", in2, "-out", out]

        self._call()


class RAxMLTreeBuilder(CLI):
    def __init__(self) -> None:
        super().__init__()
        self.args = ["raxmlHPC"]

    def call(self, m: str, n: str, p: int, f: str, s: Path, t: Path, w: Path):
        self.args += [
            "-m",
            m,
            "-n",
            n,
            "-p",
            str(p),
            "-f",
            f,
            "-s",
            s,
            "-t",
            t,
            "-w",
            w.resolve(),
        ]
        self._call()
