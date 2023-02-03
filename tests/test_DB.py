import os
import pytest
import shutil
import tempfile
from GenusFinder.DB import DB
from pathlib import Path

@pytest.fixture
def db_fixture():
    temp_dir = tempfile.mkdtemp()
    api_key = os.environ.get("NCBI_API_KEY", "")
    yield DB(temp_dir, api_key), Path(temp_dir)
    shutil.rmtree(temp_dir)

def test_get_16S_db(db_fixture):
    db, root_fp = db_fixture
    assert db.get_16S_db().exists()
    assert db.get_16S_db().stat().st_size > 0