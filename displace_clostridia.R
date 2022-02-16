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
subDir <- paste0("Displaced_clostridia/",gsub("-","_",Sys.Date()))
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

###############Read input files########################
col_day_start <- c("Baseline"= "#7570b3","Vanco" = "#a6761d","Early recovery" = "#1b9e77",
                   "Late recovery" = "#d95f02","Early no vanco" = "#e7298a","Late no vanco" = "#66a61e")
# Import metadata 
meta <-  read.csv("../data/Microbiome/2019_07_24_Ph1_metadata_all.csv")
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
abun_dt <-  read.csv("../data/Microbiome/2019-07-24 all_one_codex_abundance.csv")
one_codex_data <- abun_dt[abun_dt$sample_name %in% meta$sample_name,]
sum(unique(one_codex_data$sample_name) %in% meta$sample_name)
# Remove rows with NA is Abundance
one_codex_data <- one_codex_data[ !is.na(one_codex_data$abundance),]
# Only keep rows with non zero abundance
one_codex_data <-  one_codex_data[one_codex_data$abundance != 0,]

# Select only species
one_codex_data <- one_codex_data[ one_codex_data$tax_rank %in% c("species"),]
head(one_codex_data)
one_codex_data <- unique(one_codex_data)

tail(one_codex_data$abundance)

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


# Clean up the names of taxawith unknowns 
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

phy_sel <-  phy
# Only want Pre and post-Vanco samples:
# Remove Cohort Vanco
#sm_dt <-  data.frame(sample_data(phy))
phy_sel <-  subset_samples(phy_sel, !cohort_id_long %in% c("Vanco","Cohort 6") )
# Only select Cohort 4 and Cohort 5
#phy_sel <-  subset_samples(phy_sel, cohort_id_long %in% c("Cohort 4","Cohort 5") )

phy_sel <-  prune_taxa(taxa_sums(phy_sel)>0, phy_sel)

# Baseline samples / Pre
phy_pre <- subset_samples(phy_sel,Day.from.Start == "Baseline")

# Keep only prevalence species:
### filter by prevalence
prevdf = apply(X = otu_table(phy_pre),
               MARGIN = ifelse(taxa_are_rows(phy_pre), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phy_pre),
                    tax_table(phy_pre))
prevdf$OTU <-  rownames(prevdf)

prevdf <- prevdf[order(prevdf$TotalAbundance,decreasing = T), ]
prevdf$OTU <- factor(as.character(prevdf$OTU) , levels = as.character(prevdf$OTU))
prevdf$Prop <- prevdf$TotalAbundance/sum(prevdf$TotalAbundance)

prev_frac<-0.1
prev_cutoff <- prev_frac*nsamples(phy_pre) # Cut-off
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

phy_pre <- prune_taxa(as.character(prevdf_fil$OTU), phy_pre)



# Step 1: Subset samples around final week of dosing
# Subset samples less than 40 days
phy_post_dose <- subset_samples(phy_sel,Day.from.Start %in% c("Late recovery"))
sm_dt <- data.frame(sample_data(phy_post_dose))
# Post VE303 samples Day 14- Day 30 (Week 3 and Week 4)
phy_post_dose <- prune_taxa(taxa_names(phy_pre),phy_post_dose)
phy_post_dose <- subset_samples(phy_post_dose,Collection.Day.Norm >=14 & Collection.Day.Norm <=30)
sm_dt <- data.frame(sample_data(phy_post_dose))
# Subset samples less than 30 days
phy_mer <-  merge_phyloseq(phy_pre,phy_post_dose)
phy_mer <-  prune_taxa(taxa_sums(phy_mer)>0, phy_mer)


# Post VE303 samples
phy_post <- subset_samples(phy_sel,Day.from.Start %in% c("Early recovery","Late recovery"))
sm_dt <- data.frame(sample_data(phy_post))
# Post VE303 samples limit to Day 50
phy_post  <- subset_samples(phy_post,Collection.Day.Norm >=14 & Collection.Day.Norm <=50)


# Only select Firmicutes
phy_mer_clos <- subset_taxa(phy_mer, Phylum == "Firmicutes")
# Perform LME to identify clostridia that are subdued compared to baseline

# Determine taxa that are significantly different between pre and post-Vanco:
result_lme <-  list()
melt_dt <-  psmelt(phy_mer_clos)

taxa <- taxa_names(phy_mer_clos)[1]
for (taxa in taxa_names(phy_mer_clos)){
  phy_m <-  melt_dt[melt_dt$OTU %in% taxa,]
  phy_m$cohort_id_long <- factor(phy_m$cohort_id_long)
  #phy_m$total_Abundance <- phy_m$Abundance * phy_m$dna_conc
  phy_m$Treatment <- ifelse(phy_m$Day.from.Start=="Baseline","Pre","Post")
  phy_m$Treatment <- factor(phy_m$Treatment,levels = c("Pre","Post"))
  phy_m$t_Abundance  <- asin(sqrt(phy_m$Abundance ))
  library("nlme")
  mod_bac <- lme(t_Abundance ~ Treatment  ,random = ~1 |cohort_id_long/Subject.ID, phy_m)
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Bac <- taxa
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[taxa]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
# Remove intercept
result_lme_dt <- result_lme_dt[result_lme_dt$Var != "(Intercept)",]
result_lme_dt$p_adj <- p.adjust(result_lme_dt$pval, method = "BH")
result_lme_dt <- result_lme_dt[order(result_lme_dt$p_adj,decreasing = F),]

p_cutoff <-  0.05
result_sig <- result_lme_dt[result_lme_dt$p_adj < p_cutoff ,]

# Filter by negative coefficient (Lower in Post compared to baseline)
result_sig <- result_sig[result_sig$Value < 0 ,]
result_sig <- result_sig[abs(result_sig$Value) >= 0.005,]

ps_sig <- prune_taxa(result_sig$Bac,phy_mer_clos)

# Visualize these species in heatmap
final_sig_dt <-  psmelt(ps_sig)
sig_tax_w <- data.frame(otu_table(ps_sig))
names(sig_tax_w)

metadata <- data.frame(sample_data(ps_sig))
metadata <- metadata  %>% arrange(plot_order) %>% data.frame()

sig_tax_w <-  sig_tax_w[,match(metadata$sample_name,colnames(sig_tax_w))]

mat_ht <- as.matrix(sig_tax_w)

row_ann_dt <-  result_sig
row_ann_dt$taxon <- row_ann_dt$Bac
row_ann_dt$Value <- ifelse(row_ann_dt$Value < 0, -1,1)
row_ann_dt <-  row_ann_dt[match(rownames(mat_ht),row_ann_dt$taxon),]

mat2 <- mat_ht
library(gtools)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
col_bar = colorRamp2(c(0, 10), c("white", "blue"))
# Taxa annotations
library(yingtools2)
library(data.table)
library(hues)
set.seed(1057)
dist_colors <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6",
                 "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3","#808000","#ffd8b1",
                 "#000080")
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
library(scico)
col_cfu <- scico(5, palette = 'lajolla')
col_time <- scico(2,palette = "devon")

ant_col <-     data.frame(Cohort= metadata$cohort_id_long,
                          Time = metadata$Day.from.Start)


col_func_time <- c("Baseline"= "#7570b3","Early recovery" = "#1b9e77","Late recovery" = "#e7298a")

ha_column = HeatmapAnnotation(Time = ant_col$Time,
                              col=list(Time = col_func_time))
splitcols <-ant_col$Cohort

library(scico)
col_mat <- scico(15, palette = 'lajolla')

mat4 <- log10(mat2 + 1e-5)
rownames(mat4) <-  gsub(".*__","",rownames(mat4))


#row_ann_dt_w <- row_ann_dt_w[match(row_ann_dt_w$taxon,rownames(mat4)),]
col_b <- c("0" = "grey","1"= "blue","-1"= "red")
ha_left = HeatmapAnnotation(Coeff_lme = row_ann_dt$Value,
                            col = list(Coeff_lme = col_b), which = "row")


library(calecopal)
col_mat <- cal_palette(name = "desert", n = 15, type = "continuous")
library(wacolors)
col_mat <-  wacolors$vantage[1:10]
library(wesanderson)

col_mat <- wes_palette("Zissou1", 15, type = "continuous")
col_mat <-  rev(rainbow(10)[1:7])
col_mat <-  viridis::viridis(10)


library(ComplexHeatmap)
ht  =  Heatmap(mat4,name = "Log10(Rel_Abun)",
               # row_split = rsplit[r_order],
               gap = unit(2, "mm"),
               column_split = splitcols,
               top_annotation = ha_column,
               left_annotation = ha_left,
               #col = jet.colors(5),
               col = col_mat,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = F,
               cluster_row_slices = F,
               cluster_columns = F,
               border = T,
               show_column_names = F,
               show_row_dend = F,
               row_names_max_width = max_text_width(rownames(mat4),
                                                    gp = gpar(fontsize = 14)),
               column_names_max_height = max_text_width(colnames(mat4), gp = gpar(fontsize = 12)),
               heatmap_legend_param = list(title = "RA", at = c(0,-1,-2,-3,-4,-5), 
                                           labels = c(expression("10"^0), expression("10"^-1), 
                                                      expression("10"^-2),expression("10"^-3),expression("10"^-4),expression("10"^-5))) )


pdf(paste0(results_folder,"/Heatmap_Sig_Firm_Baseline_Post_VE303.pdf"),height = 25, width =30)
draw(ht)
dev.off()


# Now post VE303 abundance:
# Import marker panel data for total ve303 abundance
m_panel_dt = read.csv("../data/Microbiome/Updated_Marker_panel/VE303_Ph2_Extended_Marker_Data_20210607.csv")

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
#mp_dt <- mp_dt[mp_dt$SampleID %in% sample_names(phy_rf),]
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

sum_ve303_dt <- sum_ve303_dt[sum_ve303_dt$sample_name %in% sample_names(phy_post),]

head(sum_ve303_dt)


# Now perform lme where post Clostridia abundace ~ total VE303 + time
result_lme <-  list()

phy_sel_post <-  prune_taxa(result_sig$Bac,phy_post)

melt_dt <-  psmelt(phy_sel_post)
melt_dt <- melt_dt[,c("Subject.ID","sample_name","OTU","Abundance","Collection.Day.Norm")]

melt_dt <-  melt_dt %>%
  inner_join(sum_ve303_dt)


result_lme <-  list()

taxa <- taxa_names(phy_sel_post)[4]

for (taxa in  taxa_names(phy_sel_post)){
  print(taxa)
  sum_mod_dt <- tryCatch({
    #phy_tax <- prune_taxa(taxa,phy_mer)
    phy_m <-  melt_dt[melt_dt$OTU %in% taxa,]
    phy_m$cohort_id_long <- factor(phy_m$cohort_id_long)
    phy_m$t_Abundance  <- asin(sqrt(phy_m$Abundance ))
    phy_m$t_sum_rel_abun  <- asin(sqrt(phy_m$sum_rel_abun ))
    phy_m$time <-  as.numeric(scale(phy_m$Collection.Day.Norm))
    library("nlme")
    mod_bac <- lme(t_Abundance ~ t_sum_rel_abun + time  ,random = ~1 |cohort_id_long/Subject.ID, phy_m)
    sum_mod <-  summary(mod_bac)
    sum_mod_dt <- data.frame(sum_mod$tTable)
    sum_mod_dt$Bac <- taxa
    sum_mod_dt$Var <-  rownames(sum_mod_dt)
    result_lme[[taxa]] <-  sum_mod_dt
    
  }, error = function(err){
    sum_mod_dt <- NA
    sum_mod_dt
  }
  )
  result_lme[[taxa]] <-  sum_mod_dt
}


result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
result_lme_dt <- result_lme_dt[!is.na(result_lme_dt$Var),]

result_lme_dt$pval_adj <-  p.adjust(result_lme_dt$pval,method = "BH")
result_lme_dt <- result_lme_dt[!result_lme_dt$Var == "(Intercept)",]
result_lme_dt <-  result_lme_dt[result_lme_dt$Var == "t_sum_rel_abun",]
result_lme_dt <-  result_lme_dt[result_lme_dt$pval < 0.1,]
result_lme_dt <-  result_lme_dt[result_lme_dt$Value < 0,]
phy_sig_post <-  prune_taxa(result_lme_dt$Bac,phy_sel_post)

# Make XY plot for these bugs

sel_melt_dt <- melt_dt[melt_dt$OTU %in% result_lme_dt$Bac,]


p_line <- ggplot(sel_melt_dt,aes(x = Abundance+1e-5, y = sum_rel_abun+1e-5, group = Subject.ID)) +
  geom_point(size = 2, color = "black", shape = 21, fill = "blue", alpha = 0.3)+
  facet_wrap(~OTU,scales = "free")+
  # geom_line()+
  theme_clean()+
  scale_y_log10()+
  scale_x_log10()+
  ylab("TOtal VE303 abundance")

print(p_line)

pdf(paste0(results_folder,paste0("/Sig_Firmicutes_Coh1_5.pdf")),height = 6, width = 10)
print(p_line)
dev.off()  



mat_ht <- data.frame(otu_table(phy_sig_post))
mat_ht <-  log10(mat_ht + 1e-6)
# Add ve303 abundance
#Arrange VE303 abundance per subject
mean_ve_abun <-  sum_ve303_dt[order(sum_ve303_dt$sum_rel_abun),]
mat_ht <- mat_ht[,match(mean_ve_abun$sample_name,colnames(mat_ht))]

library(gtools)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(data.table)
library(hues)
library(scales)
col_b <- c("0" = "grey","-1"= muted("blue"),"1"= muted("red"))
library(scico)
wacolors$sound_sunset
col_ve303 <- colorRamp2(rev(c(seq(-0.5,-5,-0.5))),wes_palette("Zissou1", 10, type = "continuous"))

coh.cols <- c("Vanco" = "#e41a1c", "Cohort 1" = "#377eb8", "Cohort 2" = "#4daf4a", "Cohort 3" = "#ff7f00",
              "Cohort 4" = "#a65628", "Cohort 5" = "#984ea3", "Cohort 6" = "#999999")

ha_column = HeatmapAnnotation(Total_VE303 = log10(mean_ve_abun$sum_rel_abun + 1e-5),
                              #Cohort = factor(mean_ve_abun$cohort_id_long),
                              col=list(Total_VE303 = col_ve303,
                                       Cohort =  coh.cols),annotation_name_side = "left",
                              annotation_legend_param = list (Total_VE303 = list(title = "VE303", at = c(0,-1,-2,-3,-4,-5), 
                                                                                 labels = c(expression("10"^0), expression("10"^-1), 
                                                                                            expression("10"^-2),expression("10"^-3),expression("10"^-4),expression("10"^-5)))))

col_mat <- viridis::viridis(10)



tax_dt_sel <- tax_dt[rownames(tax_dt_clean) %in% rownames(mat_ht),]  
tax_dt_sel <-  tax_dt_sel[match( rownames(mat_ht),rownames(tax_dt_sel)),]

phy.cols <- c("Actinobacteria" = "#ff3017", "Bacteroidetes" = "#ffd73a", "Firmicutes" = "#73c347", "Fusobacteria" = "#6bc77e", "Proteobacteria" = "#3c8fcd", "Spirochaetes" = "#a65628", "Synergistetes" = "#f781bf", "Tenericutes" = "#999999", "Verrucomicrobia" = "#d04ea6", "Euryarchaeota" = "black")

class.cols <-  c("Clostridia" = "#45752a","Erysipelotrichia" =  "#17270e")

library(hues)
ha_left <- HeatmapAnnotation(Class = tax_dt_sel$Class,
                             which = "row",
                             col = list(Class = class.cols
                             ))
splitcols <-  mean_ve_abun$cohort_id_long
#col_mat <-  scico(15, palette = 'lajolla')
library(ComplexHeatmap)
ht  =  Heatmap(as.matrix(mat_ht),name = "RA",
               # row_split = rsplit[r_order],
               gap = unit(2, "mm"),
               column_split = splitcols,
               #              #col = col_mat,
               col =  col_mat,
               top_annotation = ha_column,
               left_annotation = ha_left,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_title_rot = 0,
               #column_title = "Cohort 1-5",
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = T,
               cluster_row_slices = F,
               cluster_columns = F,
               show_column_names = F,
               border = T,
               show_row_dend = F,
               row_names_max_width = max_text_width(rownames(mat_ht),
                                                    gp = gpar(fontsize = 14)),
               column_names_max_height = max_text_width(colnames(mat_ht),
                                                        gp = gpar(fontsize = 12)),
               heatmap_legend_param = list(title = "RA", at = c(0,-1,-2,-3,-4,-5), 
                                           labels = c(expression("10"^0), expression("10"^-1), 
                                                      expression("10"^-2),expression("10"^-3),expression("10"^-4),expression("10"^-5))) )

pdf(paste0(results_folder,paste0("/Heatmap_Sig_Coh1_5.pdf")),height = 4, width = 15)
draw(ht)
dev.off()  

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
result_lme_dt <- result_lme_dt[!is.na(result_lme_dt$Var),]

result_lme_dt$pval_adj <-  p.adjust(result_lme_dt$pval,method = "BH")
result_lme_dt <- result_lme_dt[!result_lme_dt$Var == "(Intercept)",]

result_lme_dt <-  result_lme_dt[result_lme_dt$Var == "t_sum_rel_abun",]
write.csv(result_lme_dt,paste0(results_folder,paste0("/LME_displace_Clostridia.csv")))



