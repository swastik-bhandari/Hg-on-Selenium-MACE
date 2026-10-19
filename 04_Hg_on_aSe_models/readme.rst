================================
Hg on models of red amorphous Se
================================

Experiment says Hg sticks much more strongly to red amorphous Se than to trigonal Se :

  -dH_ads(red a-Se) > 85 kJ/mol ,  -dH_ads(t-Se) < 60 kJ/mol   (Aksenov/Chiera et al.)

Why ? Red amorphous Se is a mixture of Se8 rings and disordered chains, with under-coordinated
Se atoms (chain ends). The real amorphous solid is too big for a quick test, so two LOCAL models :

A) Hg on a Se8 crown ring (building block of red Se) in a big box

B) Hg at a chain-end defect on t-Se(10-10) : the topmost Se atom removed -> two 1-coordinated Se

run
---
python Hg_aSe_models_MACE_D3_01.py > Hg_aSe_models_MACE_D3_01.py_logfileSAVED

Results (MACE-MP-0+D3)
----------------------

=================  ==============  =============
Model / site       E_ads (kJ/mol)  Hg-Se 1st (A)
=================  ==============  =============
A: Se8 axis        -21.00          3.600
A: Se8 atom        -18.50          3.641
A: Se8 bond        -21.27          3.853
B: chain-end gap   -39.67          2.968
B: chain-end top   -20.49          3.625
t-Se (03 folder)   -24.37          3.669
=================  ==============  =============

Comments
~~~~~~~~
- the isolated Se8 ring binds Hg like t-Se (physisorption, cca -20 kJ/mol) - rings alone do not
  explain the experiment
- at the chain-end defect Hg binds 1.6 x more strongly (-40 kJ/mol) and comes closer to Se (2.97 A)
- during relaxation of the defect slab one chain end bonds sideways to a neighbouring chain
  (3-coordinated Se, d = 2.42 A) while the other end stays 1-coordinated. This over/under
  coordinated pair is the known "valence alternation pair" defect of amorphous Se.
  Hg binds next to the remaining dangling Se.
- so the TREND is right : under-coordinated Se binds Hg more strongly. But -40 kJ/mol is still
  far from > 85 kJ/mol.
- reason, see ../06_checks : MACE-MP-0+D3 strongly underbinds the chemical Hg-Se bond (HgSe
  formation energy wrong by cca 240 kJ/mol). Physisorption (dispersion) is fine, chemisorption is not.

Caveats : MACE-MP-0 is spin-unpolarized in practice, a dangling Se bond is a radical-like site;
one defect geometry only; no relativistic spin-orbit effects beyond what is in the PBE training data.
