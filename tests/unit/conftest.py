import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LIB = ROOT / "workflow" / "lib"

if str(LIB) not in sys.path:
    sys.path.insert(0, str(LIB))
