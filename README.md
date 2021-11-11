Code for manuscript "Machine learning approach for higher-order interactions detection to ecological communities management" by María Evarista Arellano-García, José Ariel Camacho-Gutiérrez, and Selene Solorza-Calderón.

# Notes
- Make sure "/data" folder contains experimental data in appropriate csv files.
- Make sure to have Julia lang installed, and change line 15 of "hoi-master-V4.r" to indicate Julia location.
- Run "hoi-master-V4.r" in the same folder as data folder and the remaining R files.
- Required structure for experimental data is described in hoi-master-V4.r (lines 31, 113, 212).

TODO:
- Automatize number of species, number of replicates, and time interval.
- Incorporate hyperparameters in master file instead of individual scripts.

# Custom variables

In files get_single_V4.r, get_pairwise_V4.r, get_three_V4.r, change estimated parameter intervals as required.

In file gen_sample_V4.r, modify "TOL" (lines 452, 1444, 2080), "flags" (lines 495, 1487, 2125) and "jMul" (line 2183) variables as required.

In the hoi-master-V4.r script:
1. Put csv files in /data folder. Current version is hard coded for three species, three replicates and time interval of 7 days. Files should be 1-species, 2-species and 3-species experiment replicates.
2. To estimate $a_{ij}$ parameters, modify lines 123-124 as required.
3. To estimate $b_{ijk}$ parameters, modify lines 223-224 as required.
4. To generate testing samples to mimic experimental data, modify line 332 to change number of required samples. Modify noise levels as required.
5. To generate training samples for ML models, modify line 647. Modify noise levels as required.

# hoi-master-V4 code steps

In the hoi-master-V4.r script, the following steps are performed:
- **STEP 1.** Read data and estimate r, K from 1-species data.
- **STEP 2.** Read data and estimate a_ij from 2-species data.
- **STEP 3.** Read data and estimate b_ijk from 3-species data.
- **STEP 4.** Generate testing samples.
- **STEP 5.** Perform classical tests (Bender-Case, Wootton) on synthetic samples.
- **STEP 6.** Generate synthetic training samples for ML models.
- **STEP 7.** Challenge ML models on synthetic samples.
- **STEP 8.** Predict HOI/NON-HOI label for experimental data.

# References for R Packages

- Soetaert, K., & Petzoldt, T. (2010). Inverse modelling, sensitivity and Monte Carlo analysis in R using package FME. Journal of statistical software, 33, 1-28.
- Kuhn, M. (2015). Caret: classification and regression training. Astrophysics Source Code Library, ascl-1505.
- Xavier Robin, Natacha Turck, Alexandre Hainard, Natalia Tiberti, Frédérique Lisacek, Jean-Charles Sanchez and Markus Müller (2011). “pROC: an open-source package for R and S+ to analyze and compare ROC curves”. BMC Bioinformatics, 12, p. 77. DOI: doi: 10.1186/147121051277
- Meyer, D., Dimitriadou, E., Hornik, K., Weingessel, A., Leisch, F., Chang, C. C., & Lin, C. C.  (2019). Package ‘e1071’. The R Journal.
- Pinheiro, J., Bates, D., DebRoy, S., Sarkar, D., Heisterkamp, S., Van Willigen, B., & Ranke, J. (2017). Package ‘nlme’. Linear and nonlinear mixed effects models, version, 3(1).
- Kassambara, A. (2020). Package ‘ggpubr’.
- Meyer, D., Hornik, K., & Buchta, C. (2017). Package ‘sets’.
- Petzoldt, T. (2019). Estimation of Growth Rates with Package Growthrates.
- Bache, S. M., & Wickham, H. (2020). Package ‘magrittr’.
- Li, C. (2019). JuliaCall: an R package for seamless integration between R and Julia. Journal of Open Source Software, 4(35), 1284.
- Calaway, R., Weston, S., & Tenenbaum, D. (2015). Package ‘doParallel’.
