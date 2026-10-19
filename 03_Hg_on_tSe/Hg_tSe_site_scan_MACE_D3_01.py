from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms, FixCartesian
from ase.optimize import BFGS
from mace.calculators import mace_mp
import numpy as np

# Screening of Hg adsorption sites on t-Se(10-10).
#
# The Hg atom is placed over a 6 x 6 grid covering one surface unit cell
# (a along x, c = chain axis along y). At every grid point the slab is frozen,
# Hg keeps its (x,y) and only its height z is optimized.
# The best grid points are then fully relaxed in Hg_tSe_relax_MACE_D3_02.py
#
#   E_ads = E(Hg/slab) - E(slab) - E(Hg atom)

NX, NY = 6, 6
H0 = 2.8      # Å, starting height of Hg above the topmost Se atom
FMAX = 0.02   # eV/Å

print("=" * 70)
print("Hg ON t-Se(10-10) : SITE SCREENING, MACE-MP-0+D3")
print("=" * 70)

calc = mace_mp(model="medium", device="cpu", default_dtype="float64", dispersion=True)

# reference energies
slab = read('../02_tSe_slab/tSe_slab_relaxed.vasp')
slab.calc = calc
e_slab = slab.get_potential_energy()

hg = Atoms('Hg', cell=[20, 20, 20], pbc=True)
hg.calc = calc
e_hg = hg.get_potential_energy()

print(f"\nE(slab)    = {e_slab:.6f} eV")
print(f"E(Hg atom) = {e_hg:.6f} eV")

# surface unit cell (the 3x3 supercell has 3 x 3 primitive surface cells)
a1 = slab.cell[0] / 3.0
a2 = slab.cell[1] / 3.0
z_top = slab.positions[:, 2].max()
n_se = len(slab)

print(f"Surface unit cell: |a1| = {np.linalg.norm(a1):.3f} Å, |a2| = {np.linalg.norm(a2):.3f} Å (chain axis)")
print(f"Topmost Se at z = {z_top:.3f} Å\n")

# start from the centre of the supercell, away from the cell boundary
origin = slab.cell[0] + slab.cell[1]
origin = origin / 2.0

rows = []
print(f"{'i':>3} {'j':>3} {'x (Å)':>8} {'y (Å)':>8} {'h (Å)':>8} {'E_ads (eV)':>11} {'E_ads (kJ/mol)':>15}")
print("-" * 62)
for i in range(NX):
    for j in range(NY):
        xy = origin + (i / NX) * a1 + (j / NY) * a2
        system = slab.copy()
        system.append('Hg')
        system.positions[-1] = [xy[0], xy[1], z_top + H0]
        system.set_constraint([FixAtoms(indices=range(n_se)),
                               FixCartesian(n_se, mask=(True, True, False))])
        system.calc = calc
        BFGS(system, logfile=None).run(fmax=FMAX, steps=100)
        e_ads = system.get_potential_energy() - e_slab - e_hg
        h = system.positions[-1, 2] - z_top
        rows.append((i, j, xy[0], xy[1], h, e_ads))
        print(f"{i:>3} {j:>3} {xy[0]:>8.3f} {xy[1]:>8.3f} {h:>8.3f} {e_ads:>11.4f} {e_ads*96.485:>15.2f}")

rows = np.array(rows)
np.savetxt('Hg_tSe_site_scan.dat', rows,
           header='i j x(A) y(A) height_above_top_Se(A) E_ads(eV)', fmt='%4d %4d %10.4f %10.4f %8.4f %10.5f')

best = rows[np.argsort(rows[:, 5])]
print("\nBest grid points (most negative E_ads):")
for r in best[:5]:
    print(f"  i={int(r[0])} j={int(r[1])}  x={r[2]:.3f}  y={r[3]:.3f}  E_ads = {r[5]:.4f} eV = {r[5]*96.485:.2f} kJ/mol")
print(f"\nCorrugation (max - min E_ads over the grid): {(rows[:,5].max()-rows[:,5].min())*1000:.1f} meV")
print("Scan saved in Hg_tSe_site_scan.dat")
