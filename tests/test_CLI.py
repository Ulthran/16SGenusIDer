import os
import pytest
import sys
from pathlib import Path

print(sys.path)
sys.path.append(Path(os.path.dirname(__file__)).parent)
print(sys.path)
from src.GenusFinder.CLI import CLI, MuscleAligner, RAxMLTreeBuilder


@pytest.fixture
def muscle_fixture():
    yield MuscleAligner()


@pytest.fixture
def raxml_fixture():
    yield RAxMLTreeBuilder()


def test_muscle_call(muscle_fixture):
    cli: MuscleAligner = muscle_fixture
    print(f"PYTHONPATH: {os.environ.get('PYTHONPATH')}")
    print(f"sys.path: {sys.path}")
    cli.call()
