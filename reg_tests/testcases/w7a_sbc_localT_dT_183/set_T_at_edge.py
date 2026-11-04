"""
Set edge (boundary-node) T to a single hardcoded value, DC harmonic, value-slot only.
"""
import h5py, shutil, os
from pathlib import Path

base_dir = Path(os.environ.get("DATA_DIR", "."))
IN  = base_dir / "jorek000001.h5"
OUT = base_dir / "jorek_restart.h5"

T_EDGE    = 1.0e-3     # hardcoded target edge T (comfortably above vpar_sbc_T_floor=1e-4)
VAR_T_F90 = 6          # 1-based

shutil.copy(IN, OUT)
with h5py.File(OUT, "r+") as f:
    v = f["values"][:]                       # (n_var, n_degrees, n_tor, n_nodes)
    b = f["boundary"][:].ravel()              # (n_nodes,)
    assert v.shape[0] == int(f["n_var"][()])
    assert v.shape[-1] == b.shape[0]
    v[VAR_T_F90 - 1, 0, 0, b != 0] = T_EDGE   # (var, degree=value, tor=DC, node)
    f["values"][...] = v
    print(f"\n\nMANUALLY set T={T_EDGE:.3e} at {(b!=0).sum()} boundary nodes-> {OUT}\n\n")