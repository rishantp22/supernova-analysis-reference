# Read-only extracted notebook code: /scratch/rpraka15/analyze.ipynb
from pathlib import Path
import os

file_path = Path("/scratch/svierlin/prism-v02/gsvance/supernova_data/jet3b/jet3b/run3g_50Am_jet3b_sph.2000")

print("Exists:", file_path.exists())
print("File:", file_path)

if file_path.exists():
    print("Size:", file_path.stat().st_size / 1024 / 1024, "MB")

# Look at the first few bytes to guess whether it is text or binary

with open(file_path, "rb") as f:
    first_bytes = f.read(200)

print(first_bytes[:200])

# Try reading it as plain text

try:
    with open(file_path, "r", errors="replace") as f:
        for i in range(20):
            print(f.readline().rstrip())
except Exception as e:
    print("Could not read as text:")
    print(e)

from pathlib import Path

file_path = Path("/scratch/svierlin/prism-v02/gsvance/supernova_data/jet3b/jet3b/run3g_50Am_jet3b_sph.2000")

with open(file_path, "rb") as f:
    header_lines = []

    while True:
        line = f.readline()

        # stop if empty file/end
        if not line:
            break

        # SDF header is readable ASCII; binary starts when decoding gets messy
        try:
            decoded = line.decode("ascii")
        except UnicodeDecodeError:
            break

        header_lines.append(decoded.rstrip())

        # safety limit
        if len(header_lines) > 500:
            break

for i, line in enumerate(header_lines):
    print(f"{i:03d}: {line}")

from pathlib import Path
import numpy as np
import pandas as pd

file_path = Path("/scratch/svierlin/prism-v02/gsvance/supernova_data/jet3b/jet3b/run3g_50Am_jet3b_sph.2000")

from pathlib import Path
import numpy as np
import pandas as pd

# ------------------------------------------------------------
# File path
# ------------------------------------------------------------

file_path = Path(
    "/scratch/svierlin/prism-v02/gsvance/supernova_data/jet3b/jet3b/"
    "run3g_50Am_jet3b_sph.2000"
)

print("File exists:", file_path.exists())
print("File path:", file_path)

if not file_path.exists():
    raise FileNotFoundError(file_path)

# ------------------------------------------------------------
# Find where the SDF header ends and binary data begins
# ------------------------------------------------------------

def find_sdf_data_offset(path):
    with open(path, "rb") as f:
        while True:
            line = f.readline()

            if not line:
                raise RuntimeError("Could not find '# SDF-EOH' marker.")

            if line.strip() == b"# SDF-EOH":
                return f.tell()

data_offset = find_sdf_data_offset(file_path)

print("Binary data starts at byte:", data_offset)

# ------------------------------------------------------------
# Header information from your SDF file
# ------------------------------------------------------------

npart = 945628

# byteorder = 0x78563412 means little-endian
# <f8 = little-endian double
# <f4 = little-endian float
# <u4 = little-endian unsigned int
# <i4 = little-endian int

particle_dtype = np.dtype([
    ("x", "<f8"),
    ("y", "<f8"),
    ("z", "<f8"),

    ("mass", "<f4"),

    ("vx", "<f4"),
    ("vy", "<f4"),
    ("vz", "<f4"),

    ("u", "<f4"),
    ("h", "<f4"),
    ("rho", "<f4"),
    ("drho_dt", "<f4"),
    ("udot", "<f4"),

    ("ax", "<f4"),
    ("ay", "<f4"),
    ("az", "<f4"),

    ("lax", "<f4"),
    ("lay", "<f4"),
    ("laz", "<f4"),

    ("phi", "<f4"),
    ("idt", "<f4"),

    ("nbrs", "<u4"),
    ("ident", "<u4"),
    ("windid", "<u4"),

    ("temp", "<f4"),
    ("Y_el", "<f4"),

    *[(f"f{i}", "<f4") for i in range(1, 23)],
    *[(f"p{i}", "<i4") for i in range(1, 23)],
    *[(f"m{i}", "<i4") for i in range(1, 23)],
])

print("Bytes per particle:", particle_dtype.itemsize)
print("Expected particle data size:", particle_dtype.itemsize * npart / 1024**2, "MB")

# ------------------------------------------------------------
# Read particle data using memory map
# ------------------------------------------------------------

particles = np.memmap(
    file_path,
    dtype=particle_dtype,
    mode="r",
    offset=data_offset,
    shape=(npart,)
)

print("Particle array shape:", particles.shape)
print("Available fields:")
print(particles.dtype.names)

# ------------------------------------------------------------
# Quick inspection
# ------------------------------------------------------------

print("\nFirst 5 raw particles:")
print(particles[:5])

# ------------------------------------------------------------
# Convert selected useful fields to pandas DataFrame
# ------------------------------------------------------------

selected_columns = [
    "ident",
    "x", "y", "z",
    "vx", "vy", "vz",
    "mass",
    "rho",
    "temp",
    "u",
    "h",
    "nbrs",
    "windid"
]

df = pd.DataFrame({
    col: particles[col] for col in selected_columns
})

print("\nDataFrame shape:", df.shape)
display(df.head())

# ------------------------------------------------------------
# Summary statistics
# ------------------------------------------------------------

display(df.describe())

# ------------------------------------------------------------
# Example analyses
# ------------------------------------------------------------

print("\n20 hottest particles:")
display(df.sort_values("temp", ascending=False).head(20))

print("\n20 densest particles:")
display(df.sort_values("rho", ascending=False).head(20))

# ------------------------------------------------------------
# Optional: include all f1-f22, p1-p22, m1-m22 columns
# ------------------------------------------------------------

all_columns = list(particles.dtype.names)

df_all = pd.DataFrame({
    col: particles[col] for col in all_columns
})

print("\nFull DataFrame shape:", df_all.shape)
display(df_all.head())

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

selected_out_csv = "run3g_50Am_jet3b_sph_2000_selected_fields.csv"
full_out_parquet = "run3g_50Am_jet3b_sph_2000_all_fields.parquet"

df.to_csv(selected_out_csv, index=False)
df_all.to_parquet(full_out_parquet, index=False)

print("\nSaved selected CSV:", selected_out_csv)
print("Saved full Parquet:", full_out_parquet)

# Build isotope/species mapping from the first particle

element_symbols = {
    0: "n",
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
}

species = []

for i in range(1, 23):
    Z = int(particles[f"p{i}"][0])
    N = int(particles[f"m{i}"][0])
    A = Z + N

    if Z == 0 and N == 1:
        name = "neutron"
    elif Z == 1 and N == 0:
        name = "proton"
    else:
        symbol = element_symbols.get(Z, f"Z{Z}")
        name = f"{symbol}-{A}"

    species.append({
        "field": f"f{i}",
        "Z_protons": Z,
        "N_neutrons": N,
        "A_mass_number": A,
        "species": name,
    })

species_df = pd.DataFrame(species)
display(species_df)