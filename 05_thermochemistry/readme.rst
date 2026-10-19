================================
Thermochemistry of Hg adsorption
================================

Thermochromatography gives the adsorption ENTHALPY at finite temperature, -dH_ads,
not the 0 K electronic energy E_ads. Here the thermal part is added :

- adsorbed Hg : 3 vibrational modes (only Hg displaced, surface frozen), ASE HarmonicThermo
- gas Hg : ideal monatomic gas, ASE IdealGasThermo, H = E + 5/2 kT

   dH_ads(T) = U_ads(T) - H_gas(T)

run
---
python Hg_Se_thermochem_MACE_D3_01.py > Hg_Se_thermochem_MACE_D3_01.py_logfileSAVED

Results at 298 K (kJ/mol)
-------------------------

==============  =======  ========  ======  ============  ===============
Surface         -E_ads   -dH_ads   dH-E    dS (J/mol/K)  experiment -dH
==============  =======  ========  ======  ============  ===============
t-Se(10-10)     24.4     23.1      +1.24   -77           < 60 (t-Se)
t-Se chain end  39.7     38.4      +1.25   -97           > 85 (red a-Se)
==============  =======  ========  ======  ============  ===============

Hg frequencies : 10-15 cm-1 on t-Se, 19-39 cm-1 at the chain end (very soft, "frustrated translations").

Comments
~~~~~~~~
- the thermal correction to the ENTHALPY is tiny, cca +1.2 kJ/mol. Check : 3 soft classical modes give 3kT,
  the gas atom has 5/2 kT, so dH - E = 1/2 kT = 1.24 kJ/mol at 298 K - exactly what we get.
- so for Hg the comparison with experiment is decided by E_ads itself, not by the thermal part
- dG_ads(298 K) is cca 0 on t-Se : Hg is only weakly held, consistent with the low limit for t-Se
- the modes are so soft that the harmonic ENTROPY is not reliable (a hindered/free translator
  model would be better); the enthalpy is not sensitive to this

Note on ASE (version 3.29)
~~~~~~~~~~~~~~~~~~~~~~~~~~
IdealGasThermo.get_entropy(T) with the DEFAULT pressure gives S(Hg) = 389 J/mol/K, wrong !
The default is units.bar (in eV/A^3) but the function expects Pa. With pressure = 1.0e5 (Pa) given
explicitly we get S(Hg, 298 K) = 174.97 J/mol/K, exactly the NIST value 174.97 J/mol/K.
