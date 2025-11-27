"""
Batch-build RSI 4D lookup tables (STRICT PARITY, MPI-enabled) for multiple sandy layers.

- Reuses the fast route: first builds a small forward MFP table per soil, then builds the 4D RSI table.
- Uses the same grids/dtype/physics as the baseline (strict parity).
- Appends a suffix "_mpi" to BOTH output filenames (MFP + 4D table) when running under MPI.

Usage examples:
  # MPI (8 ranks)
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  mpiexec -n 8 --bind-to core python3 fast_table_sand_batch.py

  # Serial fallback
  python3 fast_table_sand_batch.py

Outputs per layer name L:
  results/L[_mpi].npz              → 4D RSI interface table (shape (150,150,100,30), float64)
  results/L_mfp[_mpi].npz          → forward MFP table bundle (h_tab, mfp_tab, soil)
"""
from __future__ import annotations
import sys, os
sys.path.append("../../..")
sys.path.append("../../../src/")

import numpy as np

# --- Core deps (same as fast_table.py) ---
import functional.van_genuchten as vg
from functional.Perirhizal import (
    create_lookup_mpi_interp,   # strict-parity MPI builder
    create_lookup_interp,       # strict-parity serial builder
)

# --- Try MPI early so we can coordinate side effects (file writes) safely ---
try:
    from mpi4py import MPI  # type: ignore
    COMM = MPI.COMM_WORLD
    RANK = COMM.Get_rank()
    SIZE = COMM.Get_size()
    USE_MPI = True
    MPI_SUFFIX = "_mpi"
except Exception:
    COMM = None
    RANK = 0
    SIZE = 1
    USE_MPI = False
    MPI_SUFFIX = ""

# --- Where to write results ---
BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Soil3_look_ups")
if RANK == 0:
    os.makedirs(BASE_DIR, exist_ok=True)
    
if USE_MPI:
    COMM.Barrier()

# --- Define all sandy layers (name, [theta_r, theta_s, alpha, n, Ks]) ---
SOILS = [
    ("layer0_af", [0.056, 0.372, 0.0162, 4.48, 218.0]),
    ("layer1_af", [0.056, 0.372, 0.0162, 4.48, 218.0]),
    ("layer2_af", [0.047, 0.386, 0.0191, 3.82, 529.0]),
    ("layer3_af", [0.001, 0.707, 0.0452, 1.41, 538.0]),
    ("layer4_af", [0.033, 0.340, 0.0127, 2.046, 96.0]),
    ("layer5_af", [0.019, 0.253, 0.0120, 1.922, 36.67]),

    ("layer0_nf", [0.048, 0.351, 0.0190, 4.887, 421.67]),
    ("layer1_nf", [0.048, 0.351, 0.0190, 4.887, 421.67]),
    ("layer2_nf", [0.062, 0.371, 0.0160, 3.573, 373.3]),
    ("layer3_nf", [0.046, 0.368, 0.0152, 2.736, 105.67]),
    ("layer4_nf", [0.021, 0.304, 0.01529, 1.89, 127.67]),
    ("layer5_nf", [0.028, 0.243, 0.0124, 1.999, 7.89]),
]

# Optional: limit to a subset by name
# SUBSET = {"layer0_af", "layer1_af"}
SUBSET = None  # set to a set[str] to filter
if SUBSET is not None:
    SOILS = [item for item in SOILS if item[0] in SUBSET]

# --- Helper: build a single layer ---
def build_one(name: str, params: list[float]) -> None:
    # Filenames with optional MPI suffix
    fourd_path_no_ext = os.path.join(BASE_DIR, f"{name}{MPI_SUFFIX}")
    mfp_npz_path      = os.path.join(BASE_DIR, f"{name}_mfp{MPI_SUFFIX}.npz")

    # 1) Forward MFP (rank 0 creates & saves; others wait)
    if RANK == 0:
        sp = vg.Parameters(params)
        h_tab, mfp_tab = vg.make_forward_mfp_table(sp)   # default h in [-20000, 0]
        vg.save_mfp_tables(mfp_npz_path, h_tab=h_tab, mfp_tab=mfp_tab, sp=sp)
    if USE_MPI:
        COMM.Barrier()

    # 2) Build strict-parity 4D RSI lookup using that MFP bundle
    sp = vg.Parameters(params)
    if USE_MPI:
        create_lookup_mpi_interp(fourd_path_no_ext, sp, tables_npz=mfp_npz_path)
    else:
        create_lookup_interp(fourd_path_no_ext, sp, tables_npz=mfp_npz_path)

    # 3) Sanity print from rank 0
    if USE_MPI:
        COMM.Barrier()
    if RANK == 0:
        z = np.load(f"{fourd_path_no_ext}.npz", allow_pickle=True)
        iface = z["interface"]
        print(f"\n[{name}] interface shape/dtype:", iface.shape, iface.dtype)
        print(f"[{name}] axes:",
              "rx_",  z["rx_"].shape,
              "sx_",  z["sx_"].shape,
              "akr_", z["akr_"].shape,
              "rho_", z["rho_"].shape)
        vmin, vmax = float(np.nanmin(iface)), float(np.nanmax(iface))
        print(f"[{name}] value range: [{vmin:.3f}, {vmax:.3f}]  (expected ~[-16000, 0])")
        print(f"[{name}] wrote 4D: {fourd_path_no_ext}.npz\n[{name}] wrote MFP: {mfp_npz_path}")

# --- Batch loop ---
for idx, (name, params) in enumerate(SOILS):
    if RANK == 0:
        print(f"\n=== ({idx+1}/{len(SOILS)}) Building {name} ===")
    if USE_MPI:
        COMM.Barrier()
    build_one(name, params)
    if USE_MPI:
        COMM.Barrier()

if RANK == 0:
    print("\nAll layers done.")
