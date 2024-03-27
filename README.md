# Mackintosh_et_al_2023_rearrangement_fixation
Mathematica notebooks and python scripts associated with Mackintosh et al. 2023

## Mathematica notebooks
The notebook brenthis_sweeps_chromosome_scan.nb was used for fitting sweep models to blockwise data. This requires functions from two other notebooks (sit.likelihood_inference_code.nb and sweepsInTime_S1.nb) which are also available here but were written by Bisschop et al. (2021). The notebook finite_island_model.nb was used for fitting models of intraspecific population structure.

## Python scripts
The python script get_3D_SFS.py was used to generate an unfolded 3D-SFS from a vcf (without missingness) and tsv containing which population a RG belongs to. The script six_lineage_bSFS.py was used to generate a folded bSFS. This requires a vcf and also a bed file of callable sites. The raw bSFS was formatted with format_blocks.py. Coalescent simulations were performed with anc_sweep_simulator.py. The script six_lineage_bSFS_from_sims.py was used to generate a bSFS from simulated data.

All python scripts require `docopt` and anc_sweep_simulator.py requires `msprime`. Both can be installed via conda.
