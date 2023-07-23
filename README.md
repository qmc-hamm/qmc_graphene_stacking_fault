Supporting data for the publication:
    [https://arxiv.org/abs/2307.07210](https://arxiv.org/abs/2307.07210).

### Links to recommended files: 

1. Interlayer energy of bilayer graphene from QMC
[`0_interlayer_energy/data/qmc.csv`](https://github.com/WagnerGroup/tblg_corrugation_qmc/blob/main/0_interlayer_energy/data/qmc.csv)

2. Kolmogorov-Crespi potential parameters for LAMMPS
  - Near equilibrium interlayer spacing (for relaxation calculations):
    [`1_kc_fitting/fit/QMC_kT0.004/CH_taper.KC`](https://github.com/WagnerGroup/tblg_corrugation_qmc/blob/main/1_kc_fitting/fit/QMC_kT0.004/CH_taper.KC)
  - Large interlayer spacing (for binding energy calculations):
    [`1_kc_fitting/fit/QMC_kTinf/CH_taper.KC`](https://github.com/WagnerGroup/tblg_corrugation_qmc/blob/main/1_kc_fitting/fit/QMC_kTinf/CH_taper.KC)
  
3. Relaxed twisted bilayer graphene structured from QMC-trained potential: [`2_optimized_geometry/kc_qmc/raw/simulations/0-99/dump_final.txt`](https://github.com/WagnerGroup/tblg_corrugation_qmc/blob/main/2_optimized_geometry/kc_qmc/raw/simulations/0-99/dump_final.txt)



### Explanation of directories:

`0_interlayer_energy` : Workflows and raw data from quantum Monte Carlo calculations of the interlayer energy.

`1_kc_fitting`: scripts and plots from the fits of the raw data to a Kolmogorov-Crespi model.

`2_optimized_geometry`: Minimum energy structures computed using the Kolmogorov-Crespi potential.

`3_letb`: Tight-binding calculations of the band structure of the above geometries.
