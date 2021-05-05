Code for manuscript "Machine learning approach for higher-order interactions detection to ecological communities management" by María Evarista Arellano-García, José Ariel Camacho-Gutiérrez, and Selene Solorza-Calderón.

Notes:
- Make sure "/data" folder contains experimental data in appropriate csv files.
- Make sure to have Julia lang installed, and change line 15 of "hoi-master-V4.r" to indicate Julia location.
- Run "hoi-master-V4.r" in the same folder as data folder and the remaining R files.
- Required structure for experimental data is described in hoi-master-V4.r (lines 31, 113, 212).

TODO:
- Automatize number of species, number of replicates, and time interval.
- Incorporate hyperparameters in master file instead of individual scripts.

# Summary

In files get_single_V4.r, get_pairwise_V4.r, get_three_V4.r, change estimated parameter intervals as required.

In file gen_sample_V4.r, modify "TOL" (lines 452, 1444, 2080), "flags" (lines 495, 1487, 2125) and "jMul" (line 2183) variables as required.

In the hoi-master-V4.r file:
1. Put csv files in /data folder. Current version is hard coded for three species, three replicates and time interval of 7 days. Files should be 1-species, 2-species and 3-species experiment replicates.
2. To estimate $a_{ij}$ parameters, modify lines 123-124 as required.
3. To estimate $b_{ijk}$ parameters, modify lines 223-224 as required.
4. To generate testing samples to mimic experimental data, modify line 332 to change number of required samples. Modify noise levels as required.
5. To generate training samples for ML models, modify line 647. Modify noise levels as required.
