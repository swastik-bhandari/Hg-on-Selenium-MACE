=============================
Trigonal selenium (t-Se) bulk
=============================

t-Se is the stable form of selenium : helical Se chains along the c axis (3 atoms per turn),
space group P3_121 (152), Se on Wyckoff site 3a (x, 0, 1/3), x = 0.2254.

Inside the chain the Se-Se bonds are covalent (2.37 A), between the chains there are only
weak (van der Waals) interactions (3.44 A). So dispersion is important - this is the test
if MACE-MP-0 alone, or MACE-MP-0 + D3, describes the chains packing.

Experimental data : a = 4.3662 A, c = 4.9536 A, cohesive energy cca 2.46 eV/atom (Kittel).

run
---
python tSe_bulk_relax_MACE_D3_01.py > tSe_bulk_relax_MACE_D3_01.py_logfileSAVED

Atoms and cell are relaxed (FrechetCellFilter + BFGS, fmax = 0.005 eV/A).

=============  =======  ================  ================
Property       Exp.     MACE-MP-0         MACE-MP-0+D3
=============  =======  ================  ================
a (A)          4.3662   4.6893 (+7.4%)    4.2095 (-3.6%)
c (A)          4.9536   5.1631 (+4.2%)    5.2720 (+6.4%)
d intra (A)    2.373    2.369  (-0.2%)    2.382  (+0.4%)
d inter (A)    3.436    3.792  (+10.4%)   3.416  (-0.6%)
E_coh (eV/at)  2.46     2.49   (+1.2%)    2.78   (+13%)
=============  =======  ================  ================

Comments
~~~~~~~~
- the covalent Se-Se bond is fine with both
- without D3 the chains drift apart (inter-chain distance +10 %) - MACE-MP-0 has no dispersion,
  it was trained on PBE data
- with D3 the inter-chain distance is within 1 %, so MACE-MP-0+D3 is used further
- D3 overbinds the cohesive energy (+13 %); the chains are a bit stretched along c (+6 %)

The MACE-MP-0+D3 relaxed bulk (tSe_bulk_relaxed_MACE-MP-0+D3.vasp) is used to build the surface,
so that the slab is not strained. Open the .vasp files in VESTA.
