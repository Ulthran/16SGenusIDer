import os
import shutil
import tempfile
import pytest
from .. import INC
from src.GenusFinder.command import main


@pytest.fixture
def temp_dir():
    temp_dir = tempfile.mkdtemp()
    db_fp = os.path.join(temp_dir, "16S.db")
    db_content = "> genus1\nacgtacgg\n"
    with open(db_fp, "w") as f:
        f.write(db_content)

    yield temp_dir, db_fp

    shutil.rmtree(temp_dir)


def test_seq_id(temp_dir):
    temp_dir, db_fp = temp_dir

    # main(["--seq", "acgtacgg", "--output", temp_dir])

    assert 1 == 1
