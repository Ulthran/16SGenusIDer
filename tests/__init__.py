import os
import sys
from pathlib import Path

# Necessary on CI builds for some reason
INC = sys.path.insert(0, str(Path(os.path.dirname(__file__)).parent))