from ase import Atoms
from ase.spacegroup import crystal
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from ase.io import write
from mace.calculators import mace_mp
import numpy as np

# Experimental reference values for trigonal selenium (t-Se), space group P3_121 (152)
# Se on Wyckoff 3a site (x, 0, 1/3)
EXP_A = 4.3662   # Å
EXP_C = 4.9536   # Å
EXP_X = 0.2254
EXP_INTRA_CHAIN = 2.373  # Å, Se-Se bond inside the helical chain
EXP_INTER_CHAIN = 3.436  # Å, shortest Se...Se distance between neighbouring chains
EXP_COHESIVE_ENERGY = 2.46  # eV/atom (Kittel)

print("=" * 70)
print("TRIGONAL SELENIUM (t-Se) BULK - MACE-MP-0 vs MACE-MP-0+D3")
print("=" * 70)


def build_tse():
    """Trigonal Se from the experimental crystal structure"""
    return crystal('Se', basis=[(EXP_X, 0.0, 1.0 / 3.0)], spacegroup=152,
                   cellpar=[EXP_A, EXP_A, EXP_C, 90, 90, 120])


def se_distances(atoms):
    """Return the shortest (intra-chain) and second shell (inter-chain) Se-Se distances"""
    big = atoms.repeat((3, 3, 3))
    d = big.get_all_distances(mic=True)
    d = np.unique(np.round(d[d > 0.1], 3))
    return d[0], d[1]


results = {}

for label, dispersion in [("MACE-MP-0", False), ("MACE-MP-0+D3", True)]:
    print("\n" + "-" * 35)
    print(f"{label} CALCULATIONS")
    print("-" * 35)

    calc = mace_mp(model="medium", device="cpu", default_dtype="float64",
                   dispersion=dispersion)

    # isolated Se atom in a big box, for the cohesive energy
    atom = Atoms('Se', cell=[15, 15, 15], pbc=True)
    atom.calc = calc
    e_atom = atom.get_potential_energy()

    tse = build_tse()
    tse.calc = calc
    print(f"\nRelaxing atoms + cell of t-Se ({len(tse)} atoms)")
    opt = BFGS(FrechetCellFilter(tse), trajectory=f"tSe_bulk_{label}.traj")
    opt.run(fmax=0.005, steps=300)

    a, b, c, alpha, beta, gamma = tse.cell.cellpar()
    d_intra, d_inter = se_distances(tse)
    e_coh = e_atom - tse.get_potential_energy() / len(tse)

    print(f'\n{label} Results:')
    print(f"  a                       : {a:.4f} Å")
    print(f"  c                       : {c:.4f} Å")
    print(f"  Se-Se intra-chain       : {d_intra:.4f} Å")
    print(f"  Se...Se inter-chain     : {d_inter:.4f} Å")
    print(f"  Cohesive energy         : {e_coh:.4f} eV/atom")

    write(f"tSe_bulk_relaxed_{label}.vasp", tse, format="vasp", direct=True, sort=True)
    results[label] = (a, c, d_intra, d_inter, e_coh)

# ============================================
# Comparison with Experiment
# ============================================
print("\n" + "=" * 70)
print("COMPARISON WITH EXPERIMENT")
print("=" * 70)

exp = (EXP_A, EXP_C, EXP_INTRA_CHAIN, EXP_INTER_CHAIN, EXP_COHESIVE_ENERGY)
names = ["a (Å)", "c (Å)", "d intra (Å)", "d inter (Å)", "E_coh (eV/at)"]

print(f"{'Property':<16} {'Exp.':<10} {'MACE':<18} {'MACE+D3':<18}")
print("-" * 62)
for i, n in enumerate(names):
    m = results["MACE-MP-0"][i]
    md = results["MACE-MP-0+D3"][i]
    print(f"{n:<16} {exp[i]:<10.4f} {m:<8.4f}({(m-exp[i])/exp[i]*100:+5.1f}%)  "
          f"{md:<8.4f}({(md-exp[i])/exp[i]*100:+5.1f}%)")
print("-" * 62)

print("\nRelaxed structures saved as tSe_bulk_relaxed_*.vasp (open in VESTA)")
