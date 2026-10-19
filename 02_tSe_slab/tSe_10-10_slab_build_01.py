from ase.io import read, write
from ase.build import surface
from ase.neighborlist import neighbor_list
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
import numpy as np

# Builds the t-Se(10-10) surface : the helical Se chains lie PARALLEL to the surface,
# so the surface is made of whole chains held together only by weak inter-chain (vdW) forces.
# This is the natural cleavage plane of trigonal selenium.
#
# Default ase.build.surface cuts through the helices (Se atoms with 0 or 1 bonds),
# therefore we keep whole chains only (every Se must have 2 Se-Se bonds).

N_CHAIN_LAYERS = 4   # chain layers in the slab
REPEAT = (3, 3, 1)   # 3 x 3 surface supercell
VACUUM = 10.0        # Å on each side
BOND_CUT = 2.7       # Å, Se-Se bond cutoff (bulk bond = 2.38 Å)

print("=" * 70)
print("t-Se(10-10) SLAB FROM MACE-MP-0+D3 RELAXED BULK")
print("=" * 70)

bulk = read('../01_tSe_bulk/tSe_bulk_relaxed_MACE-MP-0+D3.vasp')
a, b, c = bulk.cell.cellpar()[:3]
print(f"\nBulk lattice used: a = {a:.4f} Å, c = {c:.4f} Å  (MACE-MP-0+D3)")

# periodic slab (no vacuum yet) - we can safely find the chains
slab = surface(bulk, (1, 0, 0), N_CHAIN_LAYERS, vacuum=0.0, periodic=True)

# find chains = connected components of the Se-Se bond graph
i, j = neighbor_list('ij', slab, BOND_CUT)
graph = coo_matrix((np.ones(len(i)), (i, j)), shape=(len(slab), len(slab)))
n_chains, chain_id = connected_components(graph)
print(f"Chains found in the periodic slab: {n_chains}")

# unwrap every chain along z so that its atoms sit together, then sort chains by height
pos = slab.get_positions()
Lz = slab.cell[2, 2]
for k in range(n_chains):
    idx = np.where(chain_id == k)[0]
    z0 = pos[idx[0], 2]
    dz = pos[idx, 2] - z0
    pos[idx, 2] -= Lz * np.round(dz / Lz)
slab.set_positions(pos)

# remove periodicity in z and add the vacuum
slab.pbc = (True, True, False)
slab.center(vacuum=VACUUM, axis=2)
slab.pbc = (True, True, True)
slab = slab.repeat(REPEAT)

# sanity check : every Se has exactly 2 Se neighbours
i, j = neighbor_list('ij', slab, BOND_CUT)
cn = np.bincount(i, minlength=len(slab))
print(f"\nSlab: {len(slab)} Se atoms, {N_CHAIN_LAYERS} chain layers, {REPEAT[0]}x{REPEAT[1]} supercell")
print(f"Cell: {slab.cell.cellpar()[:3].round(3)} Å")
print(f"Se coordination numbers found: {sorted(set(cn))}  (must be [2] : whole chains only)")
assert set(cn) == {2}, "broken chains in the slab !"

zs = slab.positions[:, 2]
print(f"Slab thickness: {zs.max() - zs.min():.3f} Å")

write('tSe_10-10_slab.vasp', slab, format='vasp', direct=True, sort=True)
print("\nSlab saved as tSe_10-10_slab.vasp (open in VESTA)")
