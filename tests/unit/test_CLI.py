import pytest
from src.GenusFinder.CLI import CLI, MuscleAligner, RAxMLTreeBuilder

@pytest.fixture
def muscle_fixture():
    yield MuscleAligner()

@pytest.fixture
def raxml_fixture():
    yield RAxMLTreeBuilder()

def test_muscle_call(muscle_fixture):
    cli: MuscleAligner = muscle_fixture
    
    cli.call()