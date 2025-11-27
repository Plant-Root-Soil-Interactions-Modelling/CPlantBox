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
BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Soil3_scenario_look_ups")
if RANK == 0:
    os.makedirs(BASE_DIR, exist_ok=True)
    
if USE_MPI:
    COMM.Barrier()

# Scenario VG parameters for layer3_af
# Source: VG_params_af3_from_scenario_analysis.json (Fokker–Stokes script)
SOILS = [
    ("layer3_af_scenario1",
     [0.00122, 0.64477, 0.02042, 1.60506, 332.31264]),   # BD = 0.88929

    ("layer3_af_scenario2",
     [0.00132, 0.64481, 0.02275, 1.58101, 407.63021]),   # BD = 0.88926

    ("layer3_af_scenario3",
     [0.00382, 0.52254, 0.01191, 1.55696, 104.92684]),   # BD = 0.98611

    ("layer3_af_scenario4",
     [0.00431, 0.52296, 0.01646, 1.50886, 213.73309]),   # BD = 0.98577

    ("layer3_af_scenario5",
     [0.01307, 0.39421, 0.00624, 1.43671, 19.11943]),    # BD = 1.08775

    ("layer3_af_scenario6",
     [0.01321, 0.39431, 0.01069, 1.41266, 85.32898]),    # BD = 1.08768
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
