import logging
from pathlib import Path
import subprocess as sp
import sys

# A template wrapper class for aligner CLI interactions
class Aligner:
    def __init__(self):
        pass

    @staticmethod
    def _call(args: list):
        try:
            sp.run(args, check=True)
        except sp.CalledProcessError as e:
            logging.error(f"{e.cmd} returned code {e.returncode}")
            sys.exit()
        except sp.TimeoutExpired as e:
            logging.error(f"{e.cmd} timed out (timeout: {e.timeout})")
            sys.exit()


class MuscleAligner(Aligner):
    def __init__(self):
        super().__init__()

    def call(self, profile: bool, in1: Path, in2: Path, out: Path):
        args = ["muscle"]
        if profile:
            args.append("-profile")
        args += [in1, in2, out]

        return self._call(args)
