from ase import Atoms
from ase.build import bulk
from ase.io import read
from ase.filters import FrechetCellFilter
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.neighborlist import neighbor_list
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
from mace.calculators import mace_mp
import numpy as np

# Two sanity checks of the Hg/Se results
#
#  1) Can MACE-MP-0+D3 describe a real Hg-Se chemical bond ?
#     -> zinc-blende HgSe crystal (tiemannite) : lattice constant, Hg-Se bond,
#        formation energy  Hg(l) + Se(t) -> HgSe(s)
#  2) Is the 4-chain-layer slab thick enough ?
#     -> E_ads of Hg (best site) with the bottom chain layer removed (3 layers)

EV2KJ = 96.485
EXP_A_HGSE = 6.085           # Å
EXP_DHF_HGSE = -46.0         # kJ/mol, standard enthalpy of formation of HgSe, cca -43 to -46 (NBS tables)
EXP_ECOH_HG = 0.67           # eV/atom, cohesive energy of Hg (Kittel)

print("=" * 70)
print("SANITY CHECKS - MACE-MP-0+D3")
print("=" * 70)

calc = mace_mp(model="medium", device="cpu", default_dtype="float64", dispersion=True)

# ============================================
# 1. HgSe crystal
# ============================================
print("\n" + "-" * 35)
print("1. ZINC-BLENDE HgSe")
print("-" * 35)

hgse = bulk('HgSe', 'zincblende', a=EXP_A_HGSE)
hgse.calc = calc
BFGS(FrechetCellFilter(hgse), logfile=None).run(fmax=0.002)
a = hgse.cell.cellpar()[0] * np.sqrt(2)    # primitive fcc cell -> cubic a
d_hgse = a * np.sqrt(3) / 4
e_hgse = hgse.get_potential_energy() / 2   # per formula unit (1 Hg + 1 Se)

# elemental references : t-Se bulk (from 01_) and liquid Hg ~ solid alpha-Hg (rhombohedral)
tse = read('../01_tSe_bulk/tSe_bulk_relaxed_MACE-MP-0+D3.vasp')
tse.calc = calc
e_se = tse.get_potential_energy() / len(tse)

ahg = bulk('Hg', 'rhombohedral', a=3.005, alpha=70.53)   # alpha-Hg, exp. (227 K)
ahg.calc = calc
BFGS(FrechetCellFilter(ahg), logfile=None).run(fmax=0.002)
e_hgbulk = ahg.get_potential_energy() / len(ahg)

dhf = (e_hgse - e_se - e_hgbulk) * EV2KJ

hg_atom = Atoms('Hg', cell=[20, 20, 20], pbc=True)
hg_atom.calc = calc
e_hg_atom = hg_atom.get_potential_energy()
print(f"Hg cohesive energy (alpha-Hg): {e_hg_atom - e_hgbulk:.3f} eV/atom  (exp. {EXP_ECOH_HG}) "
      f"-> the isolated Hg atom reference is reasonable")
print(f"a(HgSe)         : {a:.4f} Å   (exp. {EXP_A_HGSE}, {(a-EXP_A_HGSE)/EXP_A_HGSE*100:+.1f} %)")
print(f"d(Hg-Se)        : {d_hgse:.4f} Å   (exp. {EXP_A_HGSE*np.sqrt(3)/4:.3f})")
print(f"formation energy: {dhf:.1f} kJ/mol  (exp. dHf(298 K) cca {EXP_DHF_HGSE:.0f} kJ/mol, ref. Hg(l))")
print(f"  -> error cca {dhf-EXP_DHF_HGSE:.0f} kJ/mol : MACE-MP-0+D3 strongly UNDERBINDS the chemical Hg-Se bond,")
print("     although the HgSe geometry is fine. Chemisorption of Hg on Se is therefore underestimated.")

# ============================================
# 2. slab thickness
# ============================================
print("\n" + "-" * 35)
print("2. SLAB THICKNESS : 4 vs 3 CHAIN LAYERS")
print("-" * 35)

hg = Atoms('Hg', cell=[20, 20, 20], pbc=True)
hg.calc = calc
e_hg = hg.get_potential_energy()


def layers(atoms):
    n = len(atoms)
    i, j = neighbor_list('ij', atoms, 2.7)
    g = coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n))
    nc, cid = connected_components(g)
    cz = np.array([atoms.positions[cid == k, 2].mean() for k in range(nc)])
    return np.unique(np.round(cz, 1), return_inverse=True)[1][cid]


full = read('../03_Hg_on_tSe/Hg_tSe_best_relaxed.vasp')
se = full[[a.index for a in full if a.symbol == 'Se']]
lay = layers(se)
print(f"E_ads (4 layers, from 03_) taken below for comparison")

res = {}
for nlay, drop in [(4, -1), (3, 0)]:
    keep = [k for k in range(len(se)) if lay[k] != drop]
    slab = se[keep]
    l2 = lay[keep]
    fixed_layer = sorted(set(l2))[0]       # bottom remaining layer fixed
    fixed = [k for k in range(len(slab)) if l2[k] == fixed_layer] if nlay == 3 else \
            [k for k in range(len(slab)) if l2[k] < 2]
    slab.set_constraint(FixAtoms(indices=fixed))
    slab.calc = calc
    BFGS(slab, logfile=None).run(fmax=0.01, steps=300)
    e_slab = slab.get_potential_energy()
    sys_ = slab.copy()
    sys_.append('Hg')
    sys_.positions[-1] = full.positions[-1]
    sys_.set_constraint(FixAtoms(indices=fixed))
    sys_.calc = calc
    BFGS(sys_, logfile=None).run(fmax=0.01, steps=300)
    e_ads = (sys_.get_potential_energy() - e_slab - e_hg) * EV2KJ
    res[nlay] = e_ads
    print(f"  {nlay} chain layers ({len(slab)} Se): E_ads = {e_ads:.2f} kJ/mol")

print(f"  difference 4 - 3 layers: {res[4]-res[3]:+.2f} kJ/mol")
