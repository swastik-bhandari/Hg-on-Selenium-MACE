from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.neighborlist import neighbor_list
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
from mace.calculators import mace_mp
import numpy as np

# Why does Hg stick much more strongly to RED AMORPHOUS Se than to TRIGONAL Se ?
#   experiment : -dH_ads(t-Se) < 60 kJ/mol , -dH_ads(red a-Se) > 85 kJ/mol
#
# Red amorphous Se is a mixture of Se8 rings and disordered chains, with many
# under-coordinated Se atoms (chain ends). We test two simple local models :
#
#   A) Hg on a Se8 crown ring (building block of red Se), in a big box
#   B) Hg at a chain-end defect of the t-Se(10-10) surface :
#      one Se removed from a top chain -> two 1-coordinated Se atoms
#
# Both are MODELS of local sites, not of the real amorphous solid.

FMAX = 0.01
EV2KJ = 96.485
BOX = 22.0

print("=" * 70)
print("Hg ON MODELS OF RED AMORPHOUS Se - MACE-MP-0+D3")
print("=" * 70)

calc = mace_mp(model="medium", device="cpu", default_dtype="float64", dispersion=True)

hg = Atoms('Hg', cell=[BOX] * 3, pbc=True)
hg.calc = calc
e_hg = hg.get_potential_energy()
print(f"\nE(Hg atom) = {e_hg:.6f} eV")


def relax(atoms, name, steps=400):
    atoms.calc = calc
    opt = BFGS(atoms, trajectory=f'{name}.traj', logfile=f'{name}.log')
    opt.run(fmax=FMAX, steps=steps)
    return atoms.get_potential_energy(), opt.nsteps


def hg_se_distances(atoms):
    ihg = [a.index for a in atoms if a.symbol == 'Hg'][0]
    ise = [a.index for a in atoms if a.symbol == 'Se']
    return np.sort(atoms.get_distances(ihg, ise, mic=True))


results = []

# ============================================
# A) Se8 crown ring
# ============================================
print("\n" + "-" * 35)
print("A) Se8 CROWN RING")
print("-" * 35)

# crown D4d ring : 8 atoms on a circle, alternately up/down
R, h = 2.643, 0.588   # Å, gives d(Se-Se) = 2.34 Å and angle ~106 deg (starting guess)
ring = Atoms('Se8', [(R * np.cos(k * np.pi / 4), R * np.sin(k * np.pi / 4), h * (-1) ** k)
                     for k in range(8)], cell=[BOX] * 3, pbc=True)
ring.center()
e_ring, n = relax(ring, 'Se8_ring')
d = ring.get_distances(0, [1, 7])
ang = ring.get_angle(7, 0, 1)
print(f"Relaxed Se8: d(Se-Se) = {d.mean():.3f} Å (exp. 2.34), angle = {ang:.1f} deg (exp. 105.7), {n} steps")
write('Se8_ring_relaxed.vasp', ring, format='vasp', direct=True)

centre = ring.get_center_of_mass()
starts = {
    'Se8_axis': centre + [0, 0, 2.8],                                  # over the crown cavity
    'Se8_atom': ring.positions[0] + 2.7 * (ring.positions[0] - centre) / np.linalg.norm(ring.positions[0] - centre),
    'Se8_bond': (ring.positions[0] + ring.positions[1]) / 2 + np.array([0, 0, 2.8]) * np.sign(ring.positions[0, 2] - centre[2]),
}
for name, p in starts.items():
    s = ring.copy()
    s.append('Hg')
    s.positions[-1] = p
    e, n = relax(s, f'Hg_{name}')
    e_ads = e - e_ring - e_hg
    dd = hg_se_distances(s)
    print(f"  {name:<10} E_ads = {e_ads:8.4f} eV = {e_ads*EV2KJ:8.2f} kJ/mol ; "
          f"Hg-Se = {dd[0]:.3f} {dd[1]:.3f} Å ; {n} steps")
    write(f'Hg_{name}_relaxed.vasp', s, format='vasp', direct=True)
    results.append(('A: ' + name, e_ads, dd[0], dd[1]))

# ============================================
# B) chain-end defect on t-Se(10-10)
# ============================================
print("\n" + "-" * 35)
print("B) CHAIN-END DEFECT ON t-Se(10-10)")
print("-" * 35)

slab = read('../02_tSe_slab/tSe_slab_relaxed.vasp')
n_se = len(slab)
i, j = neighbor_list('ij', slab, 2.7)
graph = coo_matrix((np.ones(len(i)), (i, j)), shape=(n_se, n_se))
n_chains, chain_id = connected_components(graph)
chain_z = np.array([slab.positions[chain_id == k, 2].mean() for k in range(n_chains)])
layer = np.unique(np.round(chain_z, 1), return_inverse=True)[1][chain_id]

# remove the topmost Se atom of the slab (it belongs to a top-layer chain)
k_remove = int(np.argmax(slab.positions[:, 2]))
vac_pos = slab.positions[k_remove].copy()
del slab[k_remove]
layer = np.delete(layer, k_remove)
fixed = [k for k in range(len(slab)) if layer[k] < 2]
slab.set_constraint(FixAtoms(indices=fixed))

i, j = neighbor_list('ij', slab, 2.7)
cn = np.bincount(i, minlength=len(slab))
ends = np.where(cn == 1)[0]
print(f"Removed Se atom {k_remove} ; 1-coordinated chain-end Se atoms: {list(ends)}")

e_def, n = relax(slab, 'tSe_chainend_slab')
print(f"Defect slab relaxed in {n} steps, E = {e_def:.6f} eV")
write('tSe_chainend_slab_relaxed.vasp', slab, format='vasp', direct=True)
ends_pos = slab.positions[ends]

starts = {
    'chainend_gap': vac_pos + [0, 0, 0.5],                 # Hg into the gap between the two chain ends
    'chainend_top': ends_pos[0] + [0, 0, 2.6],             # Hg on top of one chain end
}
for name, p in starts.items():
    s = slab.copy()
    s.append('Hg')
    s.positions[-1] = p
    s.set_constraint(FixAtoms(indices=fixed))
    e, n = relax(s, f'Hg_{name}')
    e_ads = e - e_def - e_hg
    dd = hg_se_distances(s)
    print(f"  {name:<13} E_ads = {e_ads:8.4f} eV = {e_ads*EV2KJ:8.2f} kJ/mol ; "
          f"Hg-Se = {dd[0]:.3f} {dd[1]:.3f} Å ; {n} steps")
    write(f'Hg_{name}_relaxed.vasp', s, format='vasp', direct=True)
    results.append(('B: ' + name, e_ads, dd[0], dd[1]))

# ============================================
# Summary
# ============================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"{'Model / site':<18} {'E_ads (eV)':>11} {'E_ads (kJ/mol)':>15} {'Hg-Se 1st (Å)':>14} {'Hg-Se 2nd (Å)':>14}")
print("-" * 76)
for name, e, d1, d2 in results:
    print(f"{name:<18} {e:>11.4f} {e*EV2KJ:>15.2f} {d1:>14.3f} {d2:>14.3f}")
print("-" * 76)
best = min([r for r in results if r[0].startswith('B')], key=lambda x: x[1])
write('Hg_aSe_best_relaxed.vasp', read(f"Hg_{best[0][3:]}_relaxed.vasp"), format='vasp', direct=True)
print(f"Most stable chain-end site: {best[0]} -> saved as Hg_aSe_best_relaxed.vasp (used in ../05_thermochemistry)")
print("Reference: Hg-Se bond in zinc-blende HgSe crystal = 2.63 Å")
print("Experiment: -dH_ads(red a-Se) > 85 kJ/mol , -dH_ads(t-Se) < 60 kJ/mol")
