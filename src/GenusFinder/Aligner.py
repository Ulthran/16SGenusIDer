import subprocess as sp

# A template wrapper class for aligner CLI interactions
class Aligner:
    def __init__(self, logger) -> None:
        self.l = logger

    def _call(args: list) -> int:
        sp.check_output(args)


class MuscleAligner(Aligner):
    def __init__(self, logger) -> None:
        super().__init__(logger)
