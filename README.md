# VE303 phase 1 study in healthy volunteers

## Reference

Dsouza M, Menon R, et al. "Colonization of the Live Biotherapeutic Product VE303 in Healthy Volunteers and Associations with the Resident Microbiota and Metabolites" (submitted)

Dsouza et al. describe VE303, a live biotherapeutic product developed for the treatment of CDI. A phase 1 study in healthy volunteers found that VE303 was safe and well tolerated. Dosing VE303 daily for up to 14 days after vancomycin led to rapid and durable VE303 strain engraftment. 

**Highlights**
- The VE303 LBP is effective at treating C. difficile infection (CDI) in mice
- In a phase 1, dose escalation study, VE303 was safe and well-tolerated
- Microbiota depletion and multi-days of dosing lead to VE303 engraftment to 1 year
- Higher VE303 doses were associated with microbiota and metabolite changes


## Analysis summary

### R Scripts

###### AUC and VE303 spike in analysis
- VE303_logisticAUC.Rmd - fit VE303 strain colonization to logistic model; calculate AUC
- VE303_spike_in_recovery.Rmd - colonization assay validation
###### Plots
- scfa_plots.R - Displays SCFAs profile over time
- marker_panel_plots.R - Generates plots for VE303 abundances over time across different cohorts
######  Microbiome recovery post VE303
- taxa_recov_comp_Vanco.R - Use linear mixed effect (LME) model to find taxa that recover during early recovery period (post VE303) compared to Vanco only cohort
- taxa_recov_comp_Vanco_heatmap.R - Creates heatmap of all significant taxon over time displaying taxa recovery
###### Microbiome diversity analysis
- diversity_comp.R -  Generate PCA/diversity plots and compare diversities across different time points within cohorts.
###### Models to predict metabolite(SCFA/BA) abundances
- SCFA_lme.R - LME analysis for SCFA and BA
- rf_metabolites.R - Creates RF models to predict metabolites (SCFA and BA) abundances as a function of microbiome.
- ALE_RF.R - Generates ALE plots/tables for each metabolites using all predictors
- ALE_vimp_plots.R - Creates heatmaps displaying important taxa associated with SCFA/BA
###### Clostridia displacement post VE303
- displace_clostridia.R - Use LME to find firmicutes that are displaced by VE303.
###### Taxa associated with VE303 colonization
- make_input_matrix.R - Creates input matrix (LFC and Baseline abundances) for RF models at different taxonomic levels. 
- run_rf_colonization.R - Runs RF models using LFC and baseline abudances of taxa to predict VE303 colonization.
- viz_colonization_rf.R - Generates heatmaps for significant taxa associated with VE303 colonization
- lme_ve303_colonization.R - Uses LME models to find taxa associated with total VE303 colonization
###### Assess potential of Secondary Bile Acid(SBA) production using VE303 genomes 
- make_heatmap_RBH.R - Creates heatmap of SBA encoding genes in VE303 using reciprocal blast hits.
