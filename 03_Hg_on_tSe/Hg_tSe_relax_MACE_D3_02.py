from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.neighborlist import neighbor_list
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
from mace.calculators import mace_mp
import numpy as np

# Full relaxation of Hg on t-Se(10-10), starting from the best points of the site scan
# (Hg_tSe_site_scan.dat). Hg and the two top chain layers are free,
# the two bottom chain layers are fixed (as for the clean slab).
#
#   E_ads = E(Hg/slab) - E(slab) - E(Hg atom)

N_START = 3        # how many distinct starting points from the scan
MIN_SEP = 1.5      # Å, starting points must be at least this far apart (in x,y)
FMAX = 0.01        # eV/Å
EV2KJ = 96.485

# Experiment : Hg on trigonal Se, thermochromatography (Aksenov/Chiera et al.)
EXP_LIMIT_TSE = 60.0   # -dH_ads < 60 kJ/mol
EXP_LIMIT_ASE = 85.0   # -dH_ads > 85 kJ/mol on red amorphous Se

print("=" * 70)
print("Hg ON t-Se(10-10) : FULL RELAXATION OF THE BEST SITES, MACE-MP-0+D3")
print("=" * 70)

calc = mace_mp(model="medium", device="cpu", default_dtype="float64", dispersion=True)

slab = read('../02_tSe_slab/tSe_slab_relaxed.vasp')
slab.calc = calc
e_slab = slab.get_potential_energy()
hg = Atoms('Hg', cell=[20, 20, 20], pbc=True)
hg.calc = calc
e_hg = hg.get_potential_energy()
n_se = len(slab)
z_top = slab.positions[:, 2].max()

# same fixed atoms as in the clean slab relaxation (two bottom chain layers)
i, j = neighbor_list('ij', slab, 2.7)
graph = coo_matrix((np.ones(len(i)), (i, j)), shape=(n_se, n_se))
n_chains, chain_id = connected_components(graph)
chain_z = np.array([slab.positions[chain_id == k, 2].mean() for k in range(n_chains)])
layer = np.unique(np.round(chain_z, 1), return_inverse=True)[1][chain_id]
fixed = [k for k in range(n_se) if layer[k] < 2]

# distinct starting points from the scan
scan = np.loadtxt('Hg_tSe_site_scan.dat')
scan = scan[np.argsort(scan[:, 5])]
starts = []
for r in scan:
    if all(np.hypot(r[2] - s[2], r[3] - s[3]) > MIN_SEP for s in starts):
        starts.append(r)
    if len(starts) == N_START:
        break

results = []
for n, r in enumerate(starts, 1):
    name = f"site{n}"
    print("\n" + "-" * 35)
    print(f"{name}: start x={r[2]:.3f} y={r[3]:.3f} h={r[4]:.3f} Å  (scan E_ads = {r[5]:.4f} eV)")
    print("-" * 35)
    system = slab.copy()
    system.append('Hg')
    system.positions[-1] = [r[2], r[3], z_top + r[4]]
    system.set_constraint(FixAtoms(indices=fixed))
    system.calc = calc
    opt = BFGS(system, trajectory=f'Hg_tSe_{name}.traj', logfile=f'Hg_tSe_{name}.log')
    opt.run(fmax=FMAX, steps=400)

    e_ads = system.get_potential_energy() - e_slab - e_hg
    d = system.get_distances(n_se, range(n_se), mic=True)
    nn = np.sort(d)[:3]
    h = system.positions[-1, 2] - system.positions[:n_se, 2].max()
    dxy = np.hypot(*(system.positions[-1, :2] - [r[2], r[3]]))
    print(f"  E_ads              = {e_ads:.4f} eV = {e_ads*EV2KJ:.2f} kJ/mol")
    print(f"  3 nearest Hg-Se    = {nn[0]:.3f} {nn[1]:.3f} {nn[2]:.3f} Å")
    print(f"  Hg height over top Se = {h:.3f} Å ; lateral move during relax = {dxy:.3f} Å")
    print(f"  BFGS steps = {opt.nsteps}")
    write(f'Hg_tSe_{name}_relaxed.vasp', system, format='vasp', direct=True)
    results.append((name, e_ads, nn[0], nn[1], h))

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"{'Site':<8} {'E_ads (eV)':>11} {'E_ads (kJ/mol)':>15} {'Hg-Se 1st (Å)':>14} {'Hg-Se 2nd (Å)':>14}")
print("-" * 66)
for name, e, d1, d2, h in sorted(results, key=lambda x: x[1]):
    print(f"{name:<8} {e:>11.4f} {e*EV2KJ:>15.2f} {d1:>14.3f} {d2:>14.3f}")
print("-" * 66)
best = min(results, key=lambda x: x[1])
print(f"\nMost stable: {best[0]}, -E_ads = {-best[1]*EV2KJ:.1f} kJ/mol")
write('Hg_tSe_best_relaxed.vasp', read(f'Hg_tSe_{best[0]}_relaxed.vasp'), format='vasp', direct=True)
print(f"Saved as Hg_tSe_best_relaxed.vasp (used in ../05_thermochemistry)")
print(f"Experiment (thermochromatography): t-Se  -dH_ads < {EXP_LIMIT_TSE:.0f} kJ/mol")
print(f"                                   red a-Se -dH_ads > {EXP_LIMIT_ASE:.0f} kJ/mol")
print("NOTE: E_ads is a 0 K electronic energy ; the thermal correction to dH is in ../05_thermochemistry")
