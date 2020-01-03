Scripts used in the manuscript Thomaz and Knowles 2020 - Molecular Ecology

Two folders:
(1) Stacks - one script to use after running the POPULATIONS module with loose parameters. This script excludes loci with high theta values and excludes the final positions of each locus due to an increase in variable sites. Outputs a whitelist to run POPULATION in STACKS again.
(2) Fastsimcoal - several scripts used to generate the SFS and in Fastsimcoal analyses
(2.1) sampleDownGeno2SFS_Dea_mod.py - modified from He and Knowles 2016 to calculate the SFS (NOTE: input needs to be a vcf file from Stacks 1.41)
(2.2) HEATMAP_vcf_SFS.r - generates a heatmap for each SFS
(2.3) Fastsimcoal_selecting_PointEstimate.r - after running several fastsimcoal replicates, it generates a table summarizing all parameter estimations and a table with the best point estimate for each scenario tested, this best point estimate is used to generate the simulated SFS for parametric bootstrap
(2.4) writeMSFS.py - script used to generate a SFS from each line in the file with all simulated SFS to run parametric bootstrap (by Qixin He)
(2.5) ScriptToOrganizeCIResults.r - script that selects the best point estimate from each parametric bootstrap run (out of 40 runs)
(2.6) Cal95CI.py - script to calculate the 95% confidence interval based on 100 point estimates parametric bootstrap from Fastsimcoal (by Qixin He)
(2.7) Fig3Plot_FastsimcoalANDsea.r - script to plot Figure 3 in Thomaz and Knowles (2020) - inputs are available on Dryad
