========================
t-Se(10-10) surface slab
========================

The (10-10) surface of t-Se contains whole helical chains lying parallel to the surface.
It is the natural cleavage plane (only weak inter-chain interactions are broken).

Note : default ase.build.surface cuts THROUGH the helices (Se atoms with 0 or 1 bonds).
The script finds the chains as connected Se-Se bonds and keeps whole chains only.
Check in the logfile : all Se coordination numbers must be 2.

Slab : 4 chain layers, 3x3 surface supercell, 108 Se atoms, 10 A vacuum on each side,
cell 12.63 x 15.82 A (y = chain axis).

run
---
python tSe_10-10_slab_build_01.py > tSe_10-10_slab_build_01.py_logfileSAVED

python tSe_slab_relax_MACE_D3_02.py > tSe_slab_relax_MACE_D3_02.py_logfileSAVED

Relaxation : 2 bottom chain layers fixed (bulk), 2 top chain layers free, BFGS fmax = 0.01 eV/A.

Result : relaxation energy -0.053 eV, max displacement 0.07 A only - the chain surface is
practically unreconstructed, as expected for a van der Waals type surface.

E(slab) = -383.593299 eV (MACE-MP-0+D3)

Pictures : slab_side.png (view along the chains), slab_top.png
