import os
import pytest
import shutil
import tempfile
from .. import INC
from src.GenusFinder.DBDir import DBDir
from pathlib import Path


@pytest.fixture
def db_fixture():
    temp_dir = tempfile.mkdtemp()
    api_key = os.environ.get("NCBI_API_KEY", "")
    yield DBDir(temp_dir, api_key)
    shutil.rmtree(temp_dir)


def test_get_16S_db(db_fixture):
    db = db_fixture
    with open(db._16S_db, "w") as f:
        f.write(" ")
    assert db.get_16S_db().exists()
    assert db.get_16S_db().stat().st_size > 0


def test_get_LTP_aligned(db_fixture):
    db = db_fixture
    with open(db.LTP_aligned_fp, "w") as f:
        f.write(" ")
    assert db.get_LTP_aligned().exists()
    assert db.get_LTP_aligned().stat().st_size > 0


def test_get_LTP_blastdb(db_fixture):
    db = db_fixture
    with open(db.LTP_blastdb_fp, "w") as f:
        f.write(" ")
    assert db.get_LTP_blastdb().exists()
    assert db.get_LTP_blastdb().stat().st_size > 0


def test_get_LTP_tree(db_fixture):
    db = db_fixture
    with open(db.LTP_tree_fp, "w") as f:
        f.write(" ")
    assert db.get_LTP_tree().exists()
    assert db.get_LTP_tree().stat().st_size > 0


def test_get_LTP_csv(db_fixture):
    db = db_fixture
    with open(db.LTP_csv_fp, "w") as f:
        f.write(" ")
    assert db.get_LTP_csv().exists()
    assert db.get_LTP_csv().stat().st_size > 0
