"""Configuration shared by the two prospectively specified cohorts."""
import os
POWER = os.environ.get("GUIDE_STUDY","initial") == "power"
SEED = 20262720 if POWER else 20261120
BOUND = 262139 if POWER else 4093
REFS = 8192 if POWER else 256
FIT_LIMIT = 400000 if POWER else 160000
PER_FAMILY = 32 if POWER else 20
PERM_PER_GEOMETRY = 4096 if POWER else 1024
GRID = [60,68,76,84,92,98,104,110]
BITS = [GRID, [64,78,94,106], [66,82,96,108], GRID*2] if POWER else [GRID,[74,98],[70,102],GRID]
