library(ggplot2)
library(DESeq2)
library(phyloseq)
library(ggthemes)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape2)
library(nlme)
library(lme4)
library(ggpubr)
library(gtools)
library(data.table)
library(here)
library(stringr)
library(tidyverse)


########Create result directory#############
mainDir <- "../results"
subDir <- "Input_matrix"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

###############Read input files########################
col_day_start <- c("Baseline"= "#7570b3","Vanco" = "#a6761d","Early recovery" = "#1b9e77",
                   "Late recovery" = "#d95f02","Early no vanco" = "#e7298a","Late no vanco" = "#66a61e")

# Import metadata 
meta <-  read.csv("../Data/2019_07_24_Ph1_metadata_all.csv")
#Filter out the Cohort 8 samples and week 24 - These are excluded from the publication
# Filter out the samples collected during Vanco treatment - these were excluded from publication
meta <- meta %>% 
  filter(!cohort_id_long %in% c("Cohort 8")) %>% 
  # filter(Day.in.Treatment.original != "Week 24") %>% 
  filter(!Timepoint.Calc.Vanco %in% c("On vanco")) 

# Set factor levels for plotting
meta$Day.in.Treatment.original <- factor(meta$Day.in.Treatment.original, levels = c("Screening", "Day -1", "Day 0", "Day 02", "Day 03", "Day 04", "Day 05", "Day 06", "Day 08", "Day 10", "Day 13", "Day 14", "Day 15", "Day 19", "Day 20", "Day 22", "Day 23", "Day 24", "Day 26", "Day 28", "Day 33", "Day 35", "Week 4", "Week 6", "Week 8", "Week 12", "Week 24", "Week 52"))
meta$Timepoint.Calc <- factor(meta$Timepoint.Calc, levels = c("Baseline", "Day -1", "Day 0", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7", "Day 8", "Day 10", "Day 13", "Day 14", "Day 15", "Week 3", "Week 4", "Week 6", "Week 8", "Week 12", "Week 24", "Week 52"))
meta$Timepoint.Calc.Vanco <- factor(meta$Timepoint.Calc.Vanco, levels = c("Baseline", "On vanco", "Vanco", "Day 2", "Day 4", "Day 6", "Day 8", "Day 10", "Week 2", "Week 3", "Week 4", "Week 5", "Week 7", "Week 11", "Week 24", "Week 52"))
meta$Day.from.Start <- factor(meta$Day.from.Start, levels = c("Baseline", "On vanco", "Vanco", "Early recovery", "Late recovery", "Early no vanco", "Late no vanco"))
meta$cohort_id_long <- factor(meta$cohort_id_long, levels = c("Vanco", "Cohort 1", "Cohort 2", "Cohort 3", "Cohort 4", "Cohort 5", "Cohort 6"))
# Create a new plot_order variable and set levels for plotting
# Also remove duplicate samples collected on the same date for each subject
library(tidyr)
sample.order <- meta %>% 
  arrange(cohort_id_long, Collection.Day.Norm, Subject.ID) %>% 
  distinct(sample_name, .keep_all = TRUE) %>% 
  select(sample_name, Subject.Cohort, Collection.Day.Norm, Subject.ID) %>% 
  unite("desc.ID", c("Subject.Cohort", "Subject.ID", "Collection.Day.Norm"), sep = "_") %>% 
  distinct(desc.ID, .keep_all = TRUE) %>% 
  mutate(plot_order = seq(1, nrow(.)))

sample.order$plot_order <- factor(sample.order$plot_order, levels = c(seq(1, nrow(sample.order))))
sample.order$desc.ID <- factor(sample.order$desc.ID, levels = c(sample.order$desc.ID))

meta <- meta %>% 
  select(-desc.ID, -plot_order) %>% 
  left_join(., sample.order, by = "sample_name") %>% 
  group_by(Subject.ID) %>% 
  distinct(Collection.Date, .keep_all = TRUE) %>% 
  ungroup()
meta <-  as.data.frame(meta)
rownames(meta) <-  meta$sample_name


# Read different tax level data
all_abun <-  readRDS("../Data/oc_tax_level_abun.rds")

one_codex_data <- all_abun[["Species"]]

# RElative abundance
rel_dt <-  one_codex_data[,c("sample_name","Species.corrected","abundance")]
rel_dt <-  rel_dt[rel_dt$abundance != 0,]
unique(rel_dt$Species.corrected)
str(rel_dt)
rel_dt$sample_name <- as.character(rel_dt$sample_name)

library(tidyr)
#tmp <- rel_dt[16740:16780,]
rel_dt <-  rel_dt %>%
  pivot_wider(names_from = sample_name,values_from = abundance,values_fill = 0) %>%
  as.data.frame()
# Remove the VE303 bugs
grep("VE303",rel_dt$Species.corrected,value = T)
colSums(rel_dt[2:ncol(rel_dt)])  
rel_dt <- rel_dt[-grep("VE303",rel_dt$Species.corrected),]  
# Taxonomy:
tax_dt <- one_codex_data[,c("superkingdom_name" ,   "phylum_name"   ,       "class_name"  ,         "order_name",           "family_name",         
                            "genus_name" ,"Species.corrected")]
names(tax_dt) <- c("Kingdom","Phylum","Class","Order","Family", "Genus","Species" )
tax_dt <- unique(tax_dt)
tax_dt <- tax_dt[tax_dt$Species %in% rel_dt$Species.corrected,]
rownames(tax_dt) <- tax_dt$Species
tax_dt$Family <- NULL
# Clean up the names of taxa with unknowns 

# Now Clean up the missing taxonomy assignment with last known
tax_dt_clean <- tax_dt
# Replace the missing NCBI_species with  species name
# # Modify taxa table
modifyTax_table <- function( tax_matrix ){
  tax_matrix <-  as.matrix(tax_matrix)
  for(i in 2:ncol(tax_matrix)) {
    index_of_emptystring <- as.numeric(which(is.na(tax_matrix[,i]), arr.ind=TRUE))
    #Loop over each NA in a column i
    if(length(index_of_emptystring) > 0){
      for(j in 1:length(index_of_emptystring) ) {
        # Loop over remaining columns
        for(k in i:(ncol(tax_matrix))) {
          # Check if the current index is NA
          if(is.na(tax_matrix[index_of_emptystring[j],k])){
            tax_matrix[index_of_emptystring[j],k] =
              paste(c("uncl_",colnames(tax_matrix)[k],"_of_",tax_matrix[index_of_emptystring[j],i-1]),collapse='')
          }
        }
      }
    }
  }
  tax_matrix <-  as.data.frame(tax_matrix)
  return(tax_matrix)
}
tax_ncbi <- modifyTax_table(tax_dt_clean)

rownames(rel_dt) <- rel_dt$Species.corrected
rel_dt$Species.corrected <-  NULL

# Phyloseq 
phy <-  phyloseq(otu_table(rel_dt,taxa_are_rows = T),
                 tax_table(as.matrix(tax_ncbi)),
                 sample_data(meta))
# Only Keep the bacterial reads
phy <- subset_taxa(phy,Kingdom == "Bacteria")

# Save input data matrix for RF models
tax_levels <- c("Class","Order", "Genus","Species" )
data_list <- list()
abun_list <- list()
phy_list <- list()
abun_coh6_list <- list()
lev_tax <- tax_levels[4]
for(lev_tax in tax_levels){
  print(paste0("Glomming at ",lev_tax))
  phy_sel <- phy
  if(lev_tax != "Species"){
    phy_sel  <- tax_glom(phy, lev_tax)
    taxa_names(phy_sel) <-  as.character(tax_table(phy_sel)[,c(lev_tax)])
  }
  
  phy_list[[lev_tax]] <- phy_sel
  
  # Only want Cohort 1-5:
  phy_sel <-  subset_samples(phy_sel, !cohort_id_long %in% c("Vanco","Cohort 6") )
  phy_sel <-  prune_taxa(taxa_sums(phy_sel)>0, phy_sel)
  # Baseline samples / Pre
  phy_pre <- subset_samples(phy_sel,Day.from.Start == "Baseline")
  phy_vanco <- subset_samples(phy_sel,Day.from.Start == "Vanco")
  sm_dt <-  data.frame(sample_data(phy_vanco))
  # Post Vanco samples for early recovery only:
  phy_post <- subset_samples(phy_vanco,Collection.Day.Norm %in% c(5,6))
  # Merge baseline samples and pre-VE303 samples 
  phy_mer <-  merge_phyloseq(phy_pre,phy_post)
  phy_mer <-  prune_taxa(taxa_sums(phy_mer)>0, phy_mer)
  
  # For RF model we increase the time points post Vanco
  # Subset samples less than 40 days
  phy_post_rf <- subset_samples(phy_sel,Day.from.Start %in% c("Early recovery","Late recovery"))
  phy_post_rf <- subset_samples(phy_post_rf,Collection.Day.Norm <= 40)
  
  phy_rf <-   merge_phyloseq(phy_pre,phy_post_rf)
  phy_rf <-  prune_taxa(taxa_sums(phy_rf)>0, phy_rf)
  
  # Keep only prevalence species:
  ### filter by prevalence
  prevdf = apply(X = otu_table(phy_mer),
                 MARGIN = ifelse(taxa_are_rows(phy_mer), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(phy_mer),
                      tax_table(phy_mer))
  prevdf$OTU <-  rownames(prevdf)
  
  prevdf <- prevdf[order(prevdf$TotalAbundance,decreasing = T), ]
  prevdf$OTU <- factor(as.character(prevdf$OTU) , levels = as.character(prevdf$OTU))
  prevdf$Prop <- prevdf$TotalAbundance/sum(prevdf$TotalAbundance)
  
  prev_frac<-0.01
  prev_cutoff <- prev_frac*nsamples(phy_mer) # Cut-off
  abun_cutoff <- 0.00001
  prevdf_fil <- prevdf[prevdf$Prevalence >= prev_cutoff & prevdf$Prop >= abun_cutoff, ]
  
  p1 <-  ggplot(prevdf,aes(x = Prevalence, y = Prop)) +  geom_jitter(size = 2)+
    xlab("Prevalence")+
    ylab("Proportion")+
    ggtitle(paste0("Prevalence_vs_Abundance"))+
    scale_y_log10()+
    geom_hline(yintercept=0.00001, linetype="dashed",
               color = "red", size=1)+
    theme_base()+
    geom_vline(xintercept=prev_cutoff, linetype="dashed",
               color = "red", size=1)+
    theme(axis.text.x=element_text(angle = -90, hjust = 0),panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(colour = "grey50"),
          panel.ontop = TRUE)
  print(p1)
  most_prevalent_taxa<-unique(rownames(prevdf_fil))
  phy_mer = prune_taxa(most_prevalent_taxa, phy_mer)
  
  # Import marker panel data:
  m_panel_dt = read.csv("../Data/VE303_Ph2_Extended_Marker_Data_20210607.csv")
  # Three replicates of this sample. Remove one
  m_panel_dt <- m_panel_dt[-grep("S04310504100010-166-S28_S28_L001_001.fastq.gz",m_panel_dt$filename),]
  
  #Generate a list of -subsampled files and  remove the non-subsetted samples
  # This was done because some files were sequenced much deeper than others due to pooling
  subsampled <- m_panel_dt %>%
    filter(grepl(pattern = "-subsampled", filename)) %>%
    distinct(filename)
  subsampled$filename <- str_replace_all(subsampled$filename, "-subsampled", "")
  m_panel_dt <- m_panel_dt %>%
    filter(!filename %in% subsampled$filename)
  m_panel_dt$filename <- str_replace_all(m_panel_dt$filename, "-subsampled", "")
  m_panel_dt$SampleID <-  gsub("-.*","",m_panel_dt$filename)
  m_panel_dt$SampleID[-grep("S",m_panel_dt$SampleID)] <-paste0("S0",m_panel_dt$SampleID[-grep("S",m_panel_dt$SampleID)])
  #Edit the sample_name to ensure consistency with the meta file. Some have the leading zero dropped.
  m_panel_dt$SampleID <- str_replace_all(m_panel_dt$SampleID, "S3", "S03")
  m_panel_dt$SampleID <- str_replace_all(m_panel_dt$SampleID, "S4", "S04")
  
  colnames(m_panel_dt)
  unique(m_panel_dt$organism)
  names(m_panel_dt)
  unique(m_panel_dt$SampleID)
  
  # Merge two tables 
  mp_dt <-  m_panel_dt %>% inner_join(meta, by =  c("SampleID" = "sample_name"))
  mp_dt$N_det <-  ifelse(mp_dt$detection_status == "Detected", 1,0)
  unique(mp_dt$SampleID)
  # Compute absolute bacterial abundance
  #Calculate the Absolute Abundance at each taxonomic depth
  ### note: 
  # f = relative abundance
  # M = mass of stool 
  # V = volume suspension (2 ml)
  # Vh = volume homogeneizer (0.25)
  # DNA = mass of DNA from the homogeneized stool
  # absolute = f * DNA / M * V / Vh
  mp_dt <- mp_dt %>% 
    mutate(dna.norm = DNA.Yield...U.00B5.g. / Weight.of.collected.sample..mg.) %>% 
    mutate(absolute_abund = targeted_panel_relative_abundance * dna.norm * (2/0.25))
  
  library(gtools)
  mp_dt$absolute_abund_adj <- mp_dt$absolute_abund
  mp_dt$absolute_abund_adj[mp_dt$N_det == 0] <- 0
  
  mp_dt$rel_abund_adj <- mp_dt$targeted_panel_relative_abundance
  mp_dt$rel_abund_adj[mp_dt$N_det == 0] <- 0
  mp_dt <-  mp_dt[,c("SampleID","Subject.ID","Subject.Cohort","cohort_id_long" ,                                               
                     "desc.ID" ,"Day.in.Treatment.original","Collection.Day.Norm",
                     "Day.from.vanco","Timepoint.Calc","Timepoint.Calc.Vanco","Day.from.Start","Collection.Date",
                     "Collection.Time","Vanco.Start.Date","organism","detection_status",
                     "DNA.Yield...U.00B5.g.",
                     "Weight.of.collected.sample..mg.",
                     "dna.norm",
                     "targeted_panel_relative_abundance","absolute_abund","absolute_abund_adj","rel_abund_adj")]
  
  
  # Total VE303 abundance 
  sum_ve303_dt <- mp_dt %>%
    group_by(Subject.ID,cohort_id_long,SampleID,Timepoint.Calc,Day.from.Start) %>% 
    summarise(sum_abs_abun = sum(absolute_abund_adj),
              sum_rel_abun =  sum(rel_abund_adj))
  names(sum_ve303_dt)[3] <- "sample_name"
  
  # For Cohort 6
  sum_ve303_coh6 <-  sum_ve303_dt %>% 
    filter(cohort_id_long == "Cohort 6")
  abun_coh6_list[[lev_tax]] <- sum_ve303_coh6
  sum_ve303_dt <- sum_ve303_dt[sum_ve303_dt$sample_name %in% sample_names(phy_rf),]
  
  # Pre-Post of Sig taxa that vary between pre and post Vanco LFC
  melt_dt <- psmelt(phy_mer)
  melt_dt$Treatment <-  ifelse(melt_dt$Day.from.Start == "Baseline","Pre","Post")
  mean_dt <-  melt_dt %>%
    select(OTU,Abundance,Treatment,Subject.ID) %>%
    group_by(OTU,Treatment,Subject.ID) %>%
    summarize(Mean = mean(Abundance)) %>% 
    pivot_wider(id_cols = c(OTU,Subject.ID),names_from = Treatment,values_from = Mean,
                values_fill = NA) 
  mean_dt <-  mean_dt[complete.cases(mean_dt),]
  
  # Baseline abundance 
  baseline_abun <-  mean_dt %>%
    select(Subject.ID,OTU,Pre)%>%
    pivot_wider(id_cols = c(Subject.ID),
                names_from = OTU,values_from = Pre,values_fill = NA)
  
  # Remove columns with sum zero
  baseline_mat <-  baseline_abun[,2:ncol(baseline_abun)]
  baseline_mat <- baseline_mat[,colSums(baseline_mat)>0]
  
  post_vanco_abun <- mean_dt %>%
    select(Subject.ID,OTU,Post)%>%
    pivot_wider(id_cols = c(Subject.ID),
                names_from = OTU,values_from = Post,values_fill = NA)
  # Remove columns with sum zero
  post_vanco_mat <-  post_vanco_abun[,2:ncol(post_vanco_abun)]
  post_vanco_mat <- post_vanco_mat[,colSums(post_vanco_mat)>0]
  
  # Pseudo abundance for samples to computed the LFC
  
  mean_dt[mean_dt$Post == 0,"Post"] <- 1e-5
  mean_dt[mean_dt$Pre == 0,"Pre"] <- 1e-5
  
  mean_dt$FC <- mean_dt$Post/mean_dt$Pre        
  mean_dt$logFC <-  log2(mean_dt$FC)
  
  #sm_dt <-  data.frame(sample_data(ps_sig))
  sm_dt <-  data.frame(sample_data(phy_mer))
  #Use only subjects present in mean_LFC
  mean_dt <- mean_dt[mean_dt$Subject.ID %in% sm_dt$Subject.ID,]
  mean_dt_LFC <-  mean_dt[,c("OTU","logFC","Subject.ID")] %>% data.frame()
  mean_LFC <- reshape(mean_dt_LFC, idvar = c("Subject.ID"), timevar = "OTU", direction = "wide")
  names(mean_LFC) <- make.names(names(mean_LFC))
  mean_LFC <-mean_LFC[complete.cases(mean_LFC), ]
  
  
  # Import Diversity:
  alpha_dt <-  read.csv("../data/data_Mel/2021-04-26 VE303_Ph1_alpha_diversity.csv")
  # Remove Inf div samples
  alpha_dt <-  alpha_dt[is.finite(alpha_dt$invsimpson),]
  alpha_dt <-  alpha_dt %>%
    select(sample_name,shannon,invsimpson)%>%
    left_join(meta)%>%
    select(sample_name,Subject.ID,shannon,invsimpson,Day.from.Start) 
  
  # Now, select  subject IDs that are present in 
  alpha_dt <- alpha_dt[alpha_dt$Subject.ID %in% mean_LFC$Subject.ID,]
  alpha_dt$Treatment <-  ifelse(alpha_dt$Day.from.Start == "Baseline","Pre","Post")
  # Average of diversity per subject pre and post
  div_mean_dt <-  alpha_dt %>%
    select(shannon,invsimpson,Treatment,Subject.ID) %>%
    group_by(Treatment,Subject.ID) %>%
    summarize(Mean_Shannon = mean(shannon),
              Mean_invsimpson = mean(invsimpson)) %>% as.data.frame()
  div_mean_dt_w <- div_mean_dt %>%
    pivot_wider(id_cols = c(Subject.ID),names_from = Treatment,values_from = c(Mean_Shannon,Mean_invsimpson),
                values_fill = NA)
  div_mean_dt_w <-div_mean_dt_w[complete.cases(div_mean_dt_w), ]
  div_mean_dt_w$FC_Shannon <-log2(div_mean_dt_w$Mean_Shannon_Post/div_mean_dt_w$Mean_Shannon_Pre)        
  div_mean_dt_w$FC_InvSimpson <- log2(div_mean_dt_w$Mean_invsimpson_Post/div_mean_dt_w$Mean_invsimpson_Pre)
  
  
  ## Baseline biomass:
  biomass_dt <-  mp_dt
  biomass_dt$Treatment <- ifelse(biomass_dt$Day.from.Start == "Baseline","Pre","Post")
  # Average of diversity per subject pre and post
  bio_dt <-  biomass_dt %>%
    select(dna.norm,Treatment,Subject.ID,SampleID) %>%
    unique()%>%
    group_by(Treatment,Subject.ID) %>%
    summarize(Mean_dna = mean(dna.norm))%>%
    pivot_wider(id_cols = c(Subject.ID),names_from = Treatment,values_from = c(Mean_dna),
                values_fill = NA)
  bio_dt <-bio_dt[complete.cases(bio_dt), ]%>%  #filter(Treatment == "Pre") %>%
    as.data.frame()
  
  bio_dt$FC_biomass <-log2(bio_dt$Post/bio_dt$Pre)        
  names(bio_dt)[2:3] <- paste0("biomass_",names(bio_dt)[2:3])  
  
  bio_dt <- bio_dt[bio_dt$Subject.ID %in% div_mean_dt_w$Subject.ID,]
  
  # Dosing information:
  # Import Dose duration 
  dose_dt <-  read.csv("../data/Microbiome/VE303_cohorts.csv")
  dose_dt <-  dose_dt[,c("Subject.Cohort","CFU_loading","loading_dose_duration","dose_duration")]
  coh_dt <-  unique(meta[,c("Subject.Cohort","Subject.ID")]) 
  coh_dt <-  coh_dt[coh_dt$Subject.ID %in% bio_dt$Subject.ID,]
  coh_dt$Subject.Cohort[coh_dt$Subject.Cohort == "Sentinel"] <- "1"
  
  CFU_dt <-  coh_dt %>%
    left_join(dose_dt)
  CFU_dt$CFU_loading <- log10(CFU_dt$CFU_loading + 1)
  CFU_dt$dose_duration <-  factor(CFU_dt$dose_duration)
  
  # Import Vanco concentration
  vanco_dt <-  read.csv("../data/Microbiome/2019-9-18 VE303_vanco.csv")
  vanco_dt <- vanco_dt[vanco_dt$Subject.ID %in% CFU_dt$Subject.ID,]
  vanco_dt$Result <- as.numeric(vanco_dt$Result)
  vanco_dt$Result[is.na(vanco_dt$Result)] <- 2
  vanco_dt$Treatment <- ifelse(vanco_dt$Day.in.Treatment.original %in% c("Day 08","Day 10"),"Post","Pre")
  vanco_dt <-  vanco_dt %>%
    select(Subject.ID,Treatment,Result) %>%
    group_by(Treatment,Subject.ID) %>%
    summarize(Mean_Vanco = mean(Result)) %>% as.data.frame()
  vanco_dt_w <- vanco_dt %>%
    pivot_wider(id_cols = c(Subject.ID),names_from = Treatment,values_from = c(Mean_Vanco),
                values_fill = NA)
  vanco_dt_w$FC_Vanco_con <-log2(vanco_dt_w$Post/vanco_dt_w$Pre)        
  
  
  # Combine all the predictors 
  set.seed(30)
  # Make a dataframe including LFC and baseline diversity/biomass
  # Baseline abundance
  names(baseline_mat)[1:ncol(baseline_mat)] <- paste0("baseline_", names(baseline_mat)) 
  base_abun_X_dt <- baseline_mat
  
  # Post-Vanco abun
  names(post_vanco_mat) <- paste0("post_vanco_", names(post_vanco_mat)) 
  post_vanco_abun_X_dt <- post_vanco_mat
  
  # Fold change in abundance of taxa with sig diff
  fc_X_tax_dt <-  mean_LFC[,2:ncol(mean_LFC)]
  
  # Baseline diversity and FC Diversity
  div_X_dt <-  div_mean_dt_w[,c("Mean_Shannon_Pre","Mean_invsimpson_Pre","FC_Shannon","FC_InvSimpson" )]
  
  # Baseline biomass and FC biomass
  biomass_X_dt <- bio_dt[,c("biomass_Pre" ,"FC_biomass")]
  
  # Dose CFUs
  dose_X_dt <-  CFU_dt[,c("CFU_loading","loading_dose_duration","dose_duration")]
  
  # Vanco Conc
  vanco_X_dt <-  vanco_dt_w
  
  
  # Total Input matrix
  tot_X_dt  <- data.frame(cbind(Subject.ID = mean_LFC$Subject.ID,base_abun_X_dt,post_vanco_abun_X_dt,
                                fc_X_tax_dt,
                                div_X_dt,
                                biomass_X_dt,
                                dose_X_dt,
                                FC_Vanco_con = vanco_dt_w$FC_Vanco_con))
  
  
  # Dependent variable : Total VE303 abundance
  abun_data <-  sum_ve303_dt %>%
    filter(Day.from.Start != "Baseline") 
  
  
  data_list[[lev_tax]] <- tot_X_dt
  abun_list[[lev_tax]] <-  abun_data
  
  
}
saveRDS(data_list,paste0(results_folder,"/input_data_rf.RDS"))
saveRDS(abun_list,paste0(results_folder,"/ve303_abun_rf.RDS"))
saveRDS(abun_coh6_list,paste0(results_folder,"/ve303_abun_coh6.RDS"))
saveRDS(phy_list,paste0(results_folder,"/phy_list.RDS"))


