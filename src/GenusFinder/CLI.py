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
        logging.info(f"Calling: {self.args}")
        try:
            sp.run(self.args, check=True)
            logging.info(f"Completed process: {' '.join(self.args)}")
        except sp.CalledProcessError as e:
            logging.error(f"{' '.join(e.cmd)} returned code {e.returncode}")
            sys.exit()
        except sp.TimeoutExpired as e:
            logging.error(f"{' '.join(e.cmd)} timed out (timeout: {e.timeout})")
            sys.exit()
        
        self.__init__()


class MuscleAligner(CLI):
    """
    v3
    """

    def call_simple(self, align: Path, output: Path):
        self.args += ["muscle", "-in", str(align), "-out", str(output)]
        self._call()

    def call_profile(self, profile: bool, in1: Path, in2: Path, out: Path):
        self.args.append("muscle")
        if profile:
            self.args.append("-profile")
        self.args += ["-in1", str(in1), "-in2", str(in2), "-out", str(out)]

        self._call()


class RAxMLTreeBuilder(CLI):
    def call(
        self,
        b: int = None,
        f: str = None,
        N: int = None,
        m: str = None,
        n: str = None,
        p: int = None,
        s: Path = None,
        t: Path = None,
        w: Path = None,
        z: Path = None,
    ):
        if n and w:
            if (w / f"RAxML_info.{n}").exists():
                logging.warning(f"{str(w / f'RAxML_info.{n}')} exists, skipping this tree-building step...")
                return None
        self.args += ["raxmlHPC"]
        self.args += ["-b", str(b)] if b else []
        self.args += ["-f", f] if f else []
        self.args += ["-N", str(N)] if N else []
        self.args += ["-m", m] if m else []
        self.args += ["-n", n] if n else []
        self.args += ["-p", str(p)] if p else []
        self.args += ["-s", str(s)] if s else []
        self.args += ["-t", str(t)] if t else []
        self.args += ["-w", str(w.resolve())] if w else []
        self.args += ["-z", str(z)] if z else []

        self._call()


class VsearchSearcher(CLI):
    def call(self, u: Path, db: Path, id: float, fp: Path):
        self.args += [
            "vsearch",
            "--usearch_global",
            str(u),
            "--db",
            str(db),
            "--id",
            str(id),
            "--fastapairs",
            str(fp),
        ]
        self._call()
