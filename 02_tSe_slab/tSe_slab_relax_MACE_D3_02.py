from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.neighborlist import neighbor_list
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
from mace.calculators import mace_mp
import numpy as np
import time

# Relaxation of the clean t-Se(10-10) slab with MACE-MP-0+D3.
# The two bottom chain layers are kept fixed at the bulk positions,
# the two top chain layers are free.

N_FIXED_LAYERS = 2
FMAX = 0.01  # eV/Å

print("=" * 70)
print("CLEAN t-Se(10-10) SLAB RELAXATION - MACE-MP-0+D3")
print("=" * 70)

slab = read('tSe_10-10_slab.vasp')

# find the chains (connected Se-Se bonds), then order the chain layers by their mean height
i, j = neighbor_list('ij', slab, 2.7)
graph = coo_matrix((np.ones(len(i)), (i, j)), shape=(len(slab), len(slab)))
n_chains, chain_id = connected_components(graph)
chain_z = np.array([slab.positions[chain_id == k, 2].mean() for k in range(n_chains)])
layer_of_chain = np.unique(np.round(chain_z, 1), return_inverse=True)[1]
layer = layer_of_chain[chain_id]   # 0 = bottom chain layer
print(f"\nAtoms per chain layer (bottom -> top): {np.bincount(layer)}")

fixed = [k for k in range(len(slab)) if layer[k] < N_FIXED_LAYERS]
slab.set_constraint(FixAtoms(indices=fixed))
print(f"Fixed atoms: {len(fixed)} of {len(slab)}")

calc = mace_mp(model="medium", device="cpu", default_dtype="float64", dispersion=True)
slab.calc = calc

t0 = time.time()
e_start = slab.get_potential_energy()
print(f"\nSingle point energy: {e_start:.6f} eV  ({time.time()-t0:.1f} s per evaluation)")

opt = BFGS(slab, trajectory='tSe_slab_relax.traj', logfile='tSe_slab_relax.log')
opt.run(fmax=FMAX, steps=300)
e_slab = slab.get_potential_energy()

print(f"\nRelaxed slab energy : {e_slab:.6f} eV")
print(f"Relaxation energy   : {e_slab - e_start:.6f} eV  ({opt.nsteps} BFGS steps)")
print(f"Max displacement of free atoms: "
      f"{np.abs(slab.positions - read('tSe_10-10_slab.vasp').positions).max():.3f} Å")

write('tSe_slab_relaxed.vasp', slab, format='vasp', direct=True)
print("\nRelaxed slab saved as tSe_slab_relaxed.vasp")
