=============
Sanity checks
=============

run
---
python checks_MACE_D3_01.py > checks_MACE_D3_01.py_logfileSAVED

1) Can MACE-MP-0+D3 describe a real Hg-Se chemical bond ?
---------------------------------------------------------

Test : zinc-blende HgSe (tiemannite) crystal, cell relaxed.

=======================  ===========  ===================
Property                 MACE+D3      Experiment
=======================  ===========  ===================
a (A)                    6.201        6.085 (+1.9 %)
d(Hg-Se) (A)             2.685        2.635
dHf(HgSe) (kJ/mol)       +197.5       cca -46 (NBS)
Ecoh(Hg) (eV/atom)       0.765        0.67 (Kittel)
=======================  ===========  ===================

- the geometry is fine, the Hg reference is fine (Hg cohesive energy close to experiment)
- BUT HgSe is predicted UNSTABLE by cca 240 kJ/mol : the foundation model strongly underbinds
  the chemical Hg-Se bond
- consequence : physisorption of Hg (dispersion, D3) is OK, chemisorption on under-coordinated
  Se (amorphous Se) is underestimated. This explains why the chain-end model gives only
  -40 kJ/mol instead of the experimental > 85 kJ/mol.
- remedy : fine-tune MACE on DFT (e.g. Quantum ESPRESSO PBE-D3) energies and forces of
  Hg/Se configurations - HgSe bulk, Hg on t-Se, Hg at Se defects.

2) Slab thickness
-----------------

Hg at the best t-Se site, slab with 4 and with 3 chain layers (re-relaxed from the same start) :

E_ads(4 layers) = -24.66 kJ/mol , E_ads(3 layers) = -24.06 kJ/mol , difference 0.6 kJ/mol

4 chain layers are enough. (03 folder gives -24.37 kJ/mol for 4 layers; the 0.3 kJ/mol difference
comes from a slightly different relaxed local minimum of the soft surface.)
