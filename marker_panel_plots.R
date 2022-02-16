# Load libraries
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
library(lmerTest)
library(ggpubr)
library(gtools)
library(data.table)
library(here)
library(stringr)

# Create a directory to store analysis
mainDir <- "../results"
subDir <- "Marker_panel_plots"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

library(stringr)

###############Read marker panel data ########################
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

dup_dt  <- stack(table(m_panel_dt$SampleID))
meta<-  read.csv("../data/Microbiome/2019_07_24_Ph1_metadata_all.csv")
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

meta <- meta %>% 
  select(sample_name, everything())


# Remove duplicate Timepoint.Calc for Marker panel plots
meta_unique <- meta %>% 
  group_by(Subject.ID) %>% 
  distinct(Timepoint.Calc, .keep_all = TRUE) %>% 
  ungroup()




dup_dt  <- stack(table(meta$sample_name))
meta$sample_name %in% m_panel_dt$SampleID
# Merge two tables 
panel_dt <-  m_panel_dt %>% inner_join(meta_unique, by =  c("SampleID" = "sample_name"))

# Make names
names(panel_dt) <- make.names(names(panel_dt))
unique(panel_dt$SampleID)
names(panel_dt)

# Subset 
panel_dt <- unique(panel_dt[,c("SampleID","Subject.ID","Subject.Cohort","cohort_id_long" ,                                               
                               "desc.ID" ,"Day.in.Treatment.original","Collection.Day.Norm",
                               "Day.from.vanco","Timepoint.Calc","Timepoint.Calc.Vanco","Day.from.Start","Collection.Date",
                               "Collection.Time","Vanco.Start.Date","organism","detection_status",
                               "targeted_panel_relative_abundance","normalized_marker_depth","n_reads")])

panel_dt$N_det <-  ifelse(panel_dt$detection_status == "Detected", 1,0)
# Lets compute total number of detected strains in each sample
num_dt <- panel_dt

library(hues)
set.seed(100)
iwanthue(11, plot=F)
# Color Schemes
det.cols <- c("Detected" = "#4daf4a", "Not detected" = "#e41a1c", "Insufficient data" = "#984ea3", "Probable" = "#377eb8")
strain.cols <- c("VE303-01" = "#1b9e77", "VE303-02" = "#d95f02", "VE303-03" = "#7570b3", "VE303-04" = "#e7298a", "VE303-05" = "#66a61e", "VE303-06" = "#e6ab02", "VE303-07" = "#a6761d", "VE303-08" = "#666666")
coh.cols <- c("Vanco" = "#0487e3", "Cohort 1" = "#dc2d00", "Cohort 2" = "#5b0076", "Cohort 3" = "#629449", "Cohort 4" = "#561600", "Cohort 5" = "#004631", "Cohort 6" = "#b87800", "Cohort 8" = "#b97696")


library(gtools)

num_dt$est_rel_abundance_adj <- 100*num_dt$targeted_panel_relative_abundance
num_dt$est_rel_abundance_adj[num_dt$N_det == 0] <- 0

num_dt$normalized_marker_depth_adj <- num_dt$normalized_marker_depth
num_dt$normalized_marker_depth_adj[num_dt$N_det == 0] <- 0



# Total number of strains per time point
sum_num_dt <- num_dt %>% 
  group_by(Subject.ID,cohort_id_long,SampleID,Timepoint.Calc) %>% 
  summarise(N_sum_det = sum(N_det),
            sum_abun = sum(est_rel_abundance_adj))


pdf(paste0(results_folder,"/Total_Num_Strain_det_MP.pdf"),width = 20,height = 4)
p <- ggplot(sum_num_dt, aes(x = Timepoint.Calc,y = N_sum_det,fill = cohort_id_long,color = cohort_id_long))+
  geom_point(size = 2,alpha = 0.3)+
  geom_boxplot(alpha = 0.5,outlier.colour = NA)+
  scale_fill_manual(name = "Cohort",values = coh.cols)+
  scale_color_manual(name = "Cohort",values = coh.cols)+
  facet_wrap(~cohort_id_long,scales = "free_x",nrow = 1)+
  scale_y_continuous(breaks = c(0:8),labels = c(0:8))+
  theme_bw()+
  ylab("Total VE303 Strains Detected (N)")+
  xlab("Time")+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"))
print(p)  



# Median of the detected strains
med_det_strain <-  sum_num_dt %>%
  group_by(cohort_id_long,Timepoint.Calc)%>%
  summarise(median_det = median(N_sum_det),Total_samples = n())


p <- ggplot(sum_num_dt, aes(x = Timepoint.Calc,y = N_sum_det,fill = cohort_id_long,group = cohort_id_long,color = cohort_id_long))+
  geom_point(size = 2,alpha = 0.3)+
  stat_summary(fun=median, geom="line",size = 1, alpha = 0.8)+
  geom_line(aes(group = Subject.ID), alpha = .2 )+
  geom_vline(xintercept = "0",color = "red",linetype = 2)+
  scale_color_manual(name = "Cohort",values = coh.cols) +
  scale_fill_manual(name = "Cohort",values = coh.cols) +
  facet_wrap(~cohort_id_long,scales = "free_x",nrow = 1)+
  scale_y_continuous(breaks = c(0:11),labels = c(0:11))+
  theme_bw()+
  ylab("Total VE303 Strains Detected (N)")+
  xlab("Time")+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"))
print(p)  
dev.off()

# Bubble plot for Number of Strains detected per subject

pdf(paste0(results_folder,"/Total_Num_Strain_det_Subject_level.pdf"),width = 20,height = 4)
p <- ggplot(sum_num_dt, aes(x = Timepoint.Calc,y = factor(Subject.ID),
                            size = N_sum_det,
                            fill = cohort_id_long,
                            group = cohort_id_long,
                            color = cohort_id_long))+
  geom_point(alpha = 0.5)+
  scale_color_manual(name = "Cohort",values = coh.cols) +
  scale_fill_manual(name = "Cohort",values = coh.cols) +
  facet_wrap(~cohort_id_long,scales = "free",nrow = 1)+
  scale_size(range = c(0,8))+
  theme_bw()+
  ylab("Subjects")+
  xlab("Time")+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"))
print(p)  
dev.off()


coh_list <- c("Vanco","Cohort 1", "Cohort 2", "Cohort 3", "Cohort 4", "Cohort 5" ,"Cohort 6") 
coh <-  coh_list[1]
plot_list <-  list()

for(coh in coh_list){
  
  
  sub_dt <-  num_dt[num_dt$cohort_id_long == coh,]
  library(dplyr)
  p <-   ggplot(sub_dt,
                aes(x = Timepoint.Calc,y = organism, 
                    size = est_rel_abundance_adj+0.01, color=organism)) +
    scale_y_discrete() +
    geom_point(alpha = 0.6, show.legend = T) + 
    scale_size(breaks = c(0.01,.1,1,10),labels = c(0.01,.1,1,10),
               range = c(0.01,10),
               limits = c(0.01,12),
               trans = "log10", name="Strain Abundance") +
    theme_bw() + 
    labs(title= element_blank(), x="Day", y="Strains",shape=1) + 
    theme(axis.title = element_text(size=16), 
          axis.text.y = element_text(size = 10,face = "bold"),
          axis.text.x = element_text(size=8, angle = 90,hjust = 1,vjust = 0.5,face = "bold"), 
          panel.grid.major.y = element_line(size = 1, color = "grey"),
          panel.grid.major.x = element_line(size = .2, color = "grey"), 
          strip.text.y = element_text(angle = 0, size = 12)) + 
    facet_grid(~Subject.ID, scales = "free_x", space = "free_y") + 
    scale_color_manual("Species",values = strain.cols) +
    guides(colour = guide_legend(override.aes = list(size=5)))+
    ggtitle(coh)
  
  plot_list[[coh]] <- p
  #print(p)
}

p0 <- plot_list[[1]]
pdf(paste0(results_folder,"/VE303_Strain_Abun_Subject_level_Vanco.pdf"),width = 12,height = 4)
p0+theme(legend.position="none")
dev.off()    
p1 <- plot_list[[2]]
pdf(paste0(results_folder,"/VE303_Strain_Abun_Subject_level_Cohort_1.pdf"),width = 8,height = 4)
p1+theme(legend.position="none")
dev.off()    
p2 <- plot_list[[3]]
pdf(paste0(results_folder,"/VE303_Strain_Abun_Subject_level_Cohort_2.pdf"),width = 8,height = 4)
p2+theme(legend.position="none")
dev.off()    
p3 <- plot_list[[4]]
pdf(paste0(results_folder,"/VE303_Strain_Abun_Subject_level_Cohort_3.pdf"),width = 8,height = 4)
p3+theme(legend.position="none")
dev.off()    
p4 <- plot_list[[5]]
pdf(paste0(results_folder,"/VE303_Strain_Abun_Subject_level_Cohort_4.pdf"),width = 15,height = 4)
p4+theme(legend.position="none")
dev.off()    
p5 <- plot_list[[6]]
pdf(paste0(results_folder,"/VE303_Strain_Abun_Subject_level_Cohort_5.pdf"),width = 18,height = 4)
p5+theme(legend.position="none")
dev.off()    
p6 <- plot_list[[7]]
pdf(paste0(results_folder,"/VE303_Strain_Abun_Subject_level_Cohort_6.pdf"),width = 12,height = 4)
p6+theme(legend.position="none")
dev.off()    

library(qpdf)
f_names <- (list.files(path = "../results/Marker_panel_plots",pattern = "VE303_Strain",full.names = T))
f_names <-  c(f_names[1:6],f_names[7])

qpdf::pdf_combine(input = f_names,
                  output = "../results/Marker_panel_plots/Combined_pdf.pdf")

# Draw legend:
# Using the cowplot package
legend <- cowplot::get_legend(p1)

pdf(paste0(results_folder,"/Legend.pdf"),width = 2,height = 5)
grid.newpage()
grid.draw(legend)
dev.off()

p2 <- plot_list[[3]]
p3 <- plot_list[[4]]
p4 <- plot_list[[5]]
p5 <- plot_list[[6]]
p6 <- plot_list[[7]]

pdf(paste0(results_folder,"/Subject_line_Strain_VE303_Log_RA_MP_Mean.pdf"),width = 20,height = 12)

p_ra <- ggplot(num_dt, aes(x = Timepoint.Calc,y = est_rel_abundance_adj+0.01,fill = cohort_id_long,group = cohort_id_long,color = cohort_id_long))+
  geom_point(size = 1,alpha = 0.5)+
  stat_summary(fun=mean, geom="line",size = 1.5, alpha = 0.5)+
  geom_line(aes(group = Subject.ID), alpha = .2 )+
  geom_hline(yintercept = .1,color = "red",linetype = 2, alpha = 0.5)+
  scale_color_manual(name = "Cohort",values = coh.cols) +
  scale_fill_manual(name = "Cohort",values = coh.cols) +
  facet_grid(organism~cohort_id_long,scales = "free_x",space = "free_x")+
  scale_y_log10(breaks = c(0.01,0.1,1,10,100),labels = c(0.01,0.1,1,10,100))+
  theme_bw()+
  ylab("% of Total VE303 Abundance")+
  xlab("Time")+
  coord_cartesian(ylim = c(0.01, 100))+
  theme(axis.text.x=element_text(size=8,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=8,face="bold"),
        axis.title.y = element_text(size=8,face="bold"),
        axis.title.x = element_text(size=8,face="bold"))
print(p_ra)  
dev.off()


pdf(paste0(results_folder,"/Total_VE303_Col_RA_MP.pdf"),width = 20,height = 4)
p <- ggplot(sum_num_dt, aes(x = Timepoint.Calc,y = sum_abun,fill = cohort_id_long,color = cohort_id_long))+
  geom_point(size = 2,alpha = 0.3)+
  geom_boxplot(alpha = 0.5,outlier.colour = NA)+
  scale_fill_manual(name = "Treatment",values = coh.cols)+
  scale_color_manual(name = "Treatment",values = coh.cols)+
  facet_wrap(~cohort_id_long,scales = "free_x",nrow = 1)+
  theme_bw()+
  ylab("% of Total VE303 Abundance")+
  xlab("Time")+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"))
print(p)


p <- ggplot(sum_num_dt, aes(x = Timepoint.Calc,y = sum_abun+0.01,fill = cohort_id_long,color = cohort_id_long))+
  geom_point(size = 2,alpha = 0.3)+
  geom_boxplot(alpha = 0.5,outlier.colour = NA)+
  geom_hline(yintercept = .1,color = "red",linetype = 2)+
  scale_fill_manual(name = "Treatment",values = coh.cols)+
  scale_color_manual(name = "Treatment",values = coh.cols)+
  facet_wrap(~cohort_id_long,scales = "free_x",nrow = 1)+
  scale_y_log10(breaks = c(0.01,0.1,1,5,10,100),labels = c(0.01,0.1,1,5,10,100))+
  theme_bw()+
  ylab("% of Total VE303 Abundance")+
  xlab("Time")+
  coord_cartesian(ylim = c(0.01, 100))+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"))
print(p)  

dev.off()


# Plot VE303 absolute abundances:
# Merge two tables 
mp_dt <-  m_panel_dt %>% inner_join(meta_unique, by =  c("SampleID" = "sample_name"))
mp_dt$N_det <-  ifelse(mp_dt$detection_status == "Detected", 1,0)

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


#tmp <- mp_dt[mp_dt$organism == "VE303-01",]

mp_dt <-  mp_dt[,c("SampleID","Subject.ID","Subject.Cohort","cohort_id_long" ,                                               
                   "desc.ID" ,"Day.in.Treatment.original","Collection.Day.Norm",
                   "Day.from.vanco","Timepoint.Calc","Timepoint.Calc.Vanco","Day.from.Start","Collection.Date",
                   "Collection.Time","Vanco.Start.Date","organism","detection_status",
                   "DNA.Yield...U.00B5.g.",
                   "Weight.of.collected.sample..mg.",
                   #"dna_norm",
                   "targeted_panel_relative_abundance","absolute_abund","absolute_abund_adj")]

sum_ve303_dt <- mp_dt %>%
  group_by(Subject.ID,cohort_id_long,SampleID,Timepoint.Calc) %>% 
  summarise(sum_abun = sum(absolute_abund_adj))

sum_ve303_dt$Timepoint.Calc <- factor(sum_ve303_dt$Timepoint.Calc, levels = c("Baseline", "Day -1", "Day 0", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7", "Day 8", "Day 10", "Day 13", "Day 14", "Day 15", "Week 3", "Week 4", "Week 6", "Week 8", "Week 12", "Week 24", "Week 52"))


phy_p <- ggplot(sum_ve303_dt, aes( y= sum_abun +1e-8 , x= Timepoint.Calc,
                                   fill = cohort_id_long,color = cohort_id_long)) + 
  geom_point(size = 2,alpha = 0.3)+
  geom_boxplot(alpha = 0.5,outlier.colour = NA)+
  scale_fill_manual(name = "Treatment",values = coh.cols)+
  scale_color_manual(name = "Treatment",values = coh.cols)+
  facet_wrap(~cohort_id_long,scales = "free_x",nrow = 1)+
  scale_y_log10()+
  theme_bw()+
  ylab("VE303 DNA(ug/Fecal(mg))")+
  xlab("Time")+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"))

pdf(paste0(results_folder,"/Absolute_VE303_abundance.pdf"),width = 20,height = 4)
print(phy_p)      
dev.off()


# Total VE303 abundances:
meta_abs <- meta_unique %>% 
  mutate(dna.norm = DNA.Yield...U.00B5.g. / Weight.of.collected.sample..mg.) %>% 
  mutate(absolute_abund = dna.norm * (2/0.25))%>%
  filter(absolute_abund>0)

meta_abs$Timepoint.Calc <- factor(meta_abs$Timepoint.Calc, levels = c("Baseline", "Day -1", "Day 0", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7", "Day 8", "Day 10", "Day 13", "Day 14", "Day 15", "Week 3", "Week 4", "Week 6", "Week 8", "Week 12", "Week 24", "Week 52"))


phy_p <- ggplot(meta_abs, aes( y= absolute_abund, x= Timepoint.Calc,
                               fill = cohort_id_long,color = cohort_id_long)) + 
  geom_point(size = 2,alpha = 0.3)+
  geom_boxplot(alpha = 0.5,outlier.colour = NA)+
  scale_fill_manual(name = "Treatment",values = coh.cols)+
  scale_color_manual(name = "Treatment",values = coh.cols)+
  facet_wrap(~cohort_id_long,scales = "free_x",nrow = 1)+
  scale_y_log10()+
  theme_bw()+
  ylab("Absolute DNA(ug/Fecal(mg))")+
  xlab("Time")+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"))

pdf(paste0(results_folder,"/Absolute_total_abundance.pdf"),width = 20,height = 4)
print(phy_p)      
dev.off()


# Major Phylum abundance 
# Read phylum level 
phyla_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_phlyum_RA.csv")
unique(phyla_dt$sample_name)
phyla_dt <- phyla_dt %>%
  filter(absolute_abund > 0)
phyla_dt <- phyla_dt[phyla_dt$sample_name %in% meta_unique$sample_name,]
phyla_dt$cohort_id_long <-  factor(phyla_dt$cohort_id_long)
phyla_dt$cohort_id_long <- relevel(phyla_dt$cohort_id_long,ref = "Vanco") 
phyla_dt$Timepoint.Calc <- factor(phyla_dt$Timepoint.Calc, levels = c("Baseline", "Day -1", "Day 0", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7", "Day 8", "Day 10", "Day 13", "Day 14", "Day 15", "Week 3", "Week 4", "Week 6", "Week 8", "Week 12", "Week 24", "Week 52"))
sub_phy_dt <-  phyla_dt[phyla_dt$phylum_name %in% c("Proteobacteria","Bacteroidetes","Firmicutes"),]

# Compute median at each time point for each cohort 
# Median of the detected strains
med_dt <-  phyla_dt %>%
  group_by(cohort_id_long,Timepoint.Calc,phylum_name)%>%
  summarise(median_det = median(absolute_abund),Total_samples = n()) 

write.csv(med_dt,paste0(results_folder,"/Median_Absolute_abundance_Phylum_level.csv"))

sub_med_dt <-  med_dt[med_dt$phylum_name %in% c("Proteobacteria","Bacteroidetes","Firmicutes"),]

# Phylum cols
phy.cols <- c("Actinobacteria" = "#e41a1c", "Bacteroidetes" = "#ffd73a", "Firmicutes" = "#4daf4a",
              "Fusobacteria" = "#ff7f00", "Proteobacteria" = "#377eb8", "Spirochaetes" = "#a65628",
              "Synergistetes" = "#f781bf", "Tenericutes" = "#999999", "Verrucomicrobia" = "#984ea3", "Euryarchaeota" = "black")

unique(phyla_dt$phylum_name)

library(ggplot2)
library(dplyr)
library(tidyr)

pdf(paste0(results_folder,"/Absolute_abundance_Phylum_level.pdf"),width = 20,height = 4)

phy_p <- ggplot(sub_phy_dt, aes( y= absolute_abund, x= Timepoint.Calc,
                                 fill = phylum_name,color = phylum_name)) + 
  geom_point(aes( group=interaction(Subject.ID, phylum_name)),size = 2,alpha = 0.15)+
  geom_line(aes( group=interaction(Subject.ID, phylum_name)),alpha = 0.15)+
  geom_point(data = sub_med_dt, aes( y= median_det, x= Timepoint.Calc,
                                     fill = phylum_name,color = phylum_name),size = 3,alpha = 0.7)+
  geom_line(data = sub_med_dt, aes( y= median_det, x= Timepoint.Calc,
                                    color = phylum_name,group = phylum_name),alpha = 0.7)+
  scale_fill_manual(name = "Phylum",values = phy.cols)+
  scale_color_manual(name = "Phylum",values = phy.cols)+
  facet_wrap(~cohort_id_long,scales = "free_x",nrow = 1)+
  scale_y_log10()+
  theme_bw()+
  ylab("Absolute DNA(ug/Fecal(mg))")+
  xlab("Time")+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"))
print(phy_p)      
#dev.off()

phy_p <- ggplot(phyla_dt, aes( y= absolute_abund, x= Timepoint.Calc,
                               fill = phylum_name,color = phylum_name)) + 
  geom_point(aes( group=interaction(Subject.ID, phylum_name)),size = 2,alpha = 0.15)+
  geom_line(aes( group=interaction(Subject.ID, phylum_name)),alpha = 0.15)+
  geom_point(data = med_dt, aes( y= median_det, x= Timepoint.Calc,
                                 fill = phylum_name,color = phylum_name),size = 3,alpha = 0.7)+
  geom_line(data = med_dt, aes( y= median_det, x= Timepoint.Calc,
                                color = phylum_name,group = phylum_name),alpha = 0.7)+
  #geom_boxplot(alpha = 0.5,outlier.colour = NA)+
  scale_fill_manual(name = "Phylum",values = phy.cols)+
  scale_color_manual(name = "Phylum",values = phy.cols)+
  facet_wrap(~cohort_id_long,scales = "free_x",nrow = 1)+
  scale_y_log10()+
  theme_bw()+
  ylab("Absolute DNA(ug/Fecal(mg))")+
  xlab("Time")+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"))

#pdf(paste0(results_folder,"/Absolute_abundance_All_Phylum_level.pdf"),width = 20,height = 4)
print(phy_p)      
dev.off()




# Do similar to Class Clostridia
# Read CLass level 
class_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_class_RA.csv")
unique(class_dt$sample_name)

class_dt <- class_dt %>% 
  mutate(dna.norm = DNA.Yield...U.00B5.g. / Weight.of.collected.sample..mg.) %>% 
  mutate(absolute_abund = rel_abund * dna.norm * (2/0.25))

class_dt <- class_dt %>%
  filter(absolute_abund > 0)

class_dt <- class_dt[class_dt$sample_name %in% meta_unique$sample_name,]
class_dt$cohort_id_long <- factor(class_dt$cohort_id_long)
class_dt$cohort_id_long <- relevel(class_dt$cohort_id_long,ref = "Vanco") 
class_dt$Timepoint.Calc <- factor(class_dt$Timepoint.Calc, levels = c("Baseline", "Day -1", "Day 0", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7", "Day 8", "Day 10", "Day 13", "Day 14", "Day 15", "Week 3", "Week 4", "Week 6", "Week 8", "Week 12", "Week 24", "Week 52"))

class_dt$GTDB_class <-  gsub("c__","",class_dt$GTDB_class)

sub_class_dt <-  class_dt[class_dt$GTDB_class %in% c("Clostridia"),]
unique(sub_class_dt$GTDB_class)

# Compute median at each time point for each cohort 
# Median of the detected strains
med_dt <-  class_dt %>%
  group_by(cohort_id_long,Timepoint.Calc,GTDB_class)%>%
  summarise(median_det = median(absolute_abund),Total_samples = n()) 
write.csv(med_dt,paste0(results_folder,"/Median_Absolute_abundance_Class_level.csv"))


sub_med_dt <-  med_dt[med_dt$GTDB_class %in% c("Clostridia"),]

# Phylum cols
class.cols <- c("Clostridia" = "#4daf4a")


unique(class_dt$GTDB_class)


library(ggplot2)
library(dplyr)
library(tidyr)



pdf(paste0(results_folder,"/Absolute_abundance_CLass_Clostridia.pdf"),width = 20,height = 4)


phy_p <- ggplot(sub_class_dt, aes( y= absolute_abund, x= Timepoint.Calc,
                                   fill = GTDB_class,color = GTDB_class)) + 
  geom_point(aes( group=interaction(Subject.ID, GTDB_class)),size = 2,alpha = 0.15)+
  geom_line(aes( group=interaction(Subject.ID, GTDB_class)),alpha = 0.15)+
  geom_point(data = sub_med_dt, aes( y= median_det, x= Timepoint.Calc,
                                     fill = GTDB_class,color = GTDB_class),size = 3,alpha = 0.7)+
  geom_line(data = sub_med_dt, aes( y= median_det, x= Timepoint.Calc,
                                    color = GTDB_class,group = GTDB_class),alpha = 0.7)+
  scale_fill_manual(name = "Class",values = class.cols)+
  scale_color_manual(name = "Class",values = class.cols)+
  facet_wrap(~cohort_id_long,scales = "free_x",nrow = 1)+
  scale_y_log10()+
  theme_bw()+
  ylab("Absolute DNA(ug/Fecal(mg))")+
  xlab("Time")+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"))
print(phy_p)      
dev.off()




