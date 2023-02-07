import os
import pytest
import shutil
import tempfile
from src.GenusFinder.OutputDir import OutputDir
from pathlib import Path


@pytest.fixture
def temp_dir():
    yield Path(tempfile.mkdtemp())

@pytest.fixture
def out_fixture(temp_dir):
    temp_dir = temp_dir
    seq_fp = temp_dir / "query.txt"
    with open(seq_fp, "w") as f:
        f.write("ACTGCTGACTGACTGACGATCGTACGTGCAGTCAGTGTCAGTCAGTACGGTACTGAC")

    yield OutputDir(temp_dir, seq_fp, False)
    
    shutil.rmtree(temp_dir)

def test_get_query(out_fixture):
    out: OutputDir = out_fixture
    assert out.get_query() == out.root_fp / "query.txt"

def test_write_probs(out_fixture):
    out: OutputDir = out_fixture
    out.write_probs({"GENUS1": 0.95, "GENUS2": 0.4})
    assert out.probs_fp.exists()

@pytest.fixture
def out_overwrite_seq_fixture(temp_dir):
    temp_dir = temp_dir

    seq = "ACTGCTGACTGACTGACGATCGTACGTGCAGTCAGTGTCAGTCAGTACGGTACTGAC"
    
    yield OutputDir(temp_dir, seq, True)
    
    shutil.rmtree(temp_dir)

def test_query_as_seq(out_overwrite_seq_fixture):
    out: OutputDir = out_overwrite_seq_fixture
    assert out.get_query().exists()

def test_overwrite(out_overwrite_seq_fixture):
    out: OutputDir = out_overwrite_seq_fixture
    assert len(os.listdir(out.root_fp)) == 1


