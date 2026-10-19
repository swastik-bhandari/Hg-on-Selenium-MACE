===========================
Hg on t-Se(10-10) with MACE
===========================

Single Hg atom adsorbed on the trigonal selenium surface (chains parallel to the surface),
MACE-MP-0 medium + D3 dispersion (torch-dftd), float64, ASE BFGS.

   E_ads = E(Hg/slab) - E(slab) - E(Hg atom)

Experiment (thermochromatography, Aksenov/Chiera et al. 2016/2017) : -dH_ads(Hg, t-Se) < 60 kJ/mol

run
---
1) site screening : Hg over a 6x6 grid of the surface unit cell, slab frozen, only Hg height optimized

python Hg_tSe_site_scan_MACE_D3_01.py > Hg_tSe_site_scan_MACE_D3_01.py_logfileSAVED

python plot_site_scan_01.py   ->  Hg_tSe_site_scan.png

2) full relaxation of the 3 best distinct grid points (Hg + 2 top chain layers free)

python Hg_tSe_relax_MACE_D3_02.py > Hg_tSe_relax_MACE_D3_02.py_logfileSAVED

Results
-------
Site screening : the surface is corrugated by cca 100 meV. Hg avoids the on-top position over
a Se atom and prefers a hollow between several outer Se atoms of the chains (see the png map).

=====  ==========  ==============  =============  =============
Site   E_ads (eV)  E_ads (kJ/mol)  Hg-Se 1st (A)  Hg-Se 2nd (A)
=====  ==========  ==============  =============  =============
site1  -0.2525     -24.37          3.669          3.860
site2  -0.2477     -23.90          3.633          3.918
site3  -0.2459     -23.73          3.641          3.829
=====  ==========  ==============  =============  =============

All 3 starting points end in equivalent hollow sites (within 0.6 kJ/mol).

- Hg is PHYSISORBED : shortest Hg-Se distance 3.67 A, much longer than the Hg-Se bond in HgSe (2.63 A)
- -E_ads = 24 kJ/mol, well inside the experimental limit -dH_ads < 60 kJ/mol for t-Se
- slab thickness checked in ../06_checks (4 vs 3 chain layers : 0.6 kJ/mol)

Hg_tSe_best_relaxed.vasp is used for the thermochemistry in ../05_thermochemistry
