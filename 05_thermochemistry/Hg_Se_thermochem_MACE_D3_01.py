from ase import Atoms
from ase.io import read
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from mace.calculators import mace_mp
import numpy as np
import sys

# Adsorption ENTHALPY of Hg on selenium surfaces, to be compared with thermochromatography.
#
# Experiment measures -dH_ads at finite temperature, while 03_ and 04_ give
# the 0 K electronic adsorption energy E_ads. Here we add the thermal part :
#
#   adsorbed Hg : 3 vibrational modes (frustrated translations), harmonic,
#                 only the Hg atom is displaced, the surface is kept frozen
#                 -> HarmonicThermo : U_ads(T) = E_ads + ZPE + vibrational energy
#   gas Hg      : ideal monatomic gas -> IdealGasThermo : H_gas(T) = 5/2 kT
#
#   dH_ads(T) = U_ads(T) - H_gas(T)        (pV of the adsorbed atom neglected)
#   dS_ads(T) = S_ads(T) - S_gas(T, 1 bar)   (pressure in Pa ; the ASE default units.bar is in eV/A^3 and gives a wrong S, see readme)
#   dG_ads(T) = dH_ads(T) - T dS_ads(T)
#
# usage : python Hg_Se_thermochem_MACE_D3_01.py
#         (structures are taken from 03_Hg_on_tSe and 04_Hg_on_aSe_models)

EV2KJ = 96.485
TEMPS = [200, 298.15, 400, 500, 600]
T_REF = 298.15
P_GAS = 1.0e5    # Pa (1 bar) - give it explicitly !

SYSTEMS = [
    # label, adsorbate structure, clean reference structure
    ("t-Se(10-10)", "../03_Hg_on_tSe/Hg_tSe_best_relaxed.vasp", "../02_tSe_slab/tSe_slab_relaxed.vasp"),
    ("t-Se chain end", "../04_Hg_on_aSe_models/Hg_aSe_best_relaxed.vasp",
     "../04_Hg_on_aSe_models/tSe_chainend_slab_relaxed.vasp"),
]

print("=" * 70)
print("THERMOCHEMISTRY OF Hg ADSORPTION ON SELENIUM - MACE-MP-0+D3")
print("=" * 70)

calc = mace_mp(model="medium", device="cpu", default_dtype="float64", dispersion=True)

# ============================================
# 1. Gas phase Hg atom
# ============================================
hg = Atoms('Hg', cell=[20, 20, 20], pbc=True)
hg.calc = calc
e_hg = hg.get_potential_energy()
gas = IdealGasThermo(vib_energies=[], potentialenergy=0.0, atoms=Atoms('Hg'),
                     geometry='monatomic', symmetrynumber=1, spin=0)
print(f"\nHg atom: H(298 K) - E = {gas.get_enthalpy(T_REF, verbose=False):.4f} eV (= 5/2 kT), "
      f"S(298 K, 1 bar) = {gas.get_entropy(T_REF, P_GAS, verbose=False)*1000:.4f} meV/K")
print("  experiment: S(Hg, gas, 298 K) = 174.97 J/(mol K) = 1.8135 meV/K (NIST)")

summary = []

for label, f_ads, f_ref in SYSTEMS:
    print("\n" + "=" * 70)
    print(f"Hg on {label}")
    print("=" * 70)

    ref = read(f_ref)
    ref.calc = calc
    e_ref = ref.get_potential_energy()

    ads = read(f_ads)
    ads.calc = calc
    e_ads = ads.get_potential_energy() - e_ref - e_hg
    ihg = [a.index for a in ads if a.symbol == 'Hg']
    print(f"E_ads (0 K, electronic) = {e_ads:.4f} eV = {e_ads*EV2KJ:.2f} kJ/mol")

    # ============================================
    # 2. Vibrations of the adsorbed Hg atom
    # ============================================
    name = 'vib_' + label.replace(' ', '_').replace('(', '').replace(')', '').replace('-', '')
    vib = Vibrations(ads, indices=ihg, name=name, delta=0.01, nfree=4)
    vib.run()
    freqs = vib.get_frequencies()
    print("\nVibrational frequencies of adsorbed Hg (cm-1):")
    for k, fr in enumerate(freqs):
        tag = "  (imaginary !)" if abs(fr.imag) > 1e-6 else ""
        print(f"  Mode {k+1}: {abs(fr):8.2f}{tag}")
    vib_energies = np.array([abs(e) for e in vib.get_energies()])   # eV, |.| for soft modes
    zpe = 0.5 * vib_energies.sum()
    print(f"ZPE = {zpe*1000:.2f} meV = {zpe*EV2KJ:.3f} kJ/mol")

    ads_thermo = HarmonicThermo(vib_energies=vib_energies, potentialenergy=e_ads)

    # ============================================
    # 3. dH, dS, dG of adsorption
    # ============================================
    print(f"\n{'T (K)':>8} {'dH_ads (kJ/mol)':>16} {'dS_ads (J/mol/K)':>17} {'dG_ads (kJ/mol)':>16}")
    print("-" * 61)
    for T in TEMPS:
        dH = ads_thermo.get_internal_energy(T, verbose=False) - gas.get_enthalpy(T, verbose=False)
        dS = ads_thermo.get_entropy(T, verbose=False) - gas.get_entropy(T, P_GAS, verbose=False)
        dG = dH - T * dS
        print(f"{T:>8.2f} {dH*EV2KJ:>16.2f} {dS*EV2KJ*1000:>17.2f} {dG*EV2KJ:>16.2f}")
        if T == T_REF:
            summary.append((label, e_ads, zpe, dH, dS, dG))
    print("-" * 61)

# ============================================
# 4. Comparison with experiment
# ============================================
print("\n" + "=" * 70)
print("COMPARISON WITH EXPERIMENT (298 K)")
print("=" * 70)
print(f"{'Surface':<16} {'-E_ads':>8} {'-dH_ads':>9} {'dH-E':>7} {'experiment -dH_ads':>22}")
print(f"{'':<16} {'kJ/mol':>8} {'kJ/mol':>9} {'kJ/mol':>7}")
print("-" * 66)
exp = {"t-Se(10-10)": "< 60 (t-Se)", "t-Se chain end": "> 85 (red a-Se)"}
for label, e_ads, zpe, dH, dS, dG in summary:
    print(f"{label:<16} {-e_ads*EV2KJ:>8.1f} {-dH*EV2KJ:>9.1f} {(dH-e_ads)*EV2KJ:>7.2f} {exp[label]:>22}")
print("-" * 66)
print("Experiment: N.M. Chiera et al., J. Radioanal. Nucl. Chem. 311 (2017) 99")
