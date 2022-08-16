import os
import shutil
import tempfile
import unittest

from GenusFinder.command import main

class CommandTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.db_fp = os.path.join(self.temp_dir, "16S.db")
        self.db_content = ("> genus1\n"
            "acgtacgg\n")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_collect_genomes(self):
        main([
            "acgtacgt",
            "--db", self.db_fp
        ])

        self.assertEqual(1, 1)