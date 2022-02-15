# Copyright (c) <2019>, <Shakti K Bhattarai & Vanni Bucci>
#   
#   Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the"Software"),
# to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#   
#   The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(dplyr)
library(phyloseq)
library(gtools)

mainDir <- "../results"
subDir <- paste0("diversity/",gsub("-","_",Sys.Date()))
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

# Read Alpha diversity data
alpha_dt <-  read.csv("../data/data_Mel/2021-04-26 VE303_Ph1_alpha_diversity.csv")
col_day_start <- c("Baseline"= "#7570b3","Vanco" = "#a6761d","Early recovery" = "#1b9e77",
                   "Late recovery" = "#d95f02","Early no vanco" = "#e7298a","Late no vanco" = "#66a61e")
# Remove Inf div samples
alpha_dt <-  alpha_dt[is.finite(alpha_dt$invsimpson),]
alpha_dt$Treatment <- alpha_dt$Day.from.Start  
alpha_dt$Treatment <-   factor(alpha_dt$Treatment , levels = names(col_day_start))

# Order by time and cohort
alpha_dt <-  alpha_dt[order(alpha_dt$cohort_id_long, alpha_dt$Treatment),]
alpha_dt$sample_name <- as.character(alpha_dt$sample_name)
alpha_dt$sample_name <- factor(alpha_dt$sample_name, levels = unique(alpha_dt$sample_name))

alpha_dt$cohort_id_long <-  factor(alpha_dt$cohort_id_long, levels = unique(alpha_dt$cohort_id_long))
alpha_dt$cohort_id_long <- relevel(alpha_dt$cohort_id_long,ref = "Vanco")

alpha_dt$Timepoint.Calc <-  as.character(alpha_dt$Timepoint.Calc)
unique(alpha_dt$Timepoint.Calc)

alpha_dt$Timepoint.Calc <-  factor(alpha_dt$Timepoint.Calc, levels = c("Baseline","Day -1","Day 0","Day 2","Day 4", "Day 5","Day 6","Day 7","Day 8","Day 10","Day 13","Day 14",
                                                                       "Week 3","Week 4","Week 6","Week 8","Week 12","Week 24","Week 52" ))
coh.cols <- c("Vanco" = "#0487e3", "Cohort 1" = "#dc2d00", "Cohort 2" = "#5b0076", "Cohort 3" = "#629449", "Cohort 4" = "#561600", "Cohort 5" = "#004631", "Cohort 6" = "#b87800", "Cohort 8" = "#b97696")
library(ggplot2)
library(ggthemes)
# Time plot 
div_index <-  c( "shannon","invsimpson","richness" )
pdf(paste(results_folder,paste0("Diversity_timeplot",'.pdf'),sep="/"),height =5 , width = 40, useDingbats = FALSE)
var <- div_index[1]
for(var in div_index){
  
  # Plot all time points  
  p_div <- ggplot(alpha_dt) +
    geom_bar(aes_string(x = "factor(sample_name)",y= var, fill = "Treatment"),size=1,stat = "identity")+
    scale_fill_manual(name = "Day.from.Start", values = col_day_start )+
    theme_base()+
    facet_grid(~cohort_id_long,scales = "free", space = "free")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    xlab("")+
    ylab(var)+
    ggtitle(var)
  
  print(p_div)
}

dev.off()

# Choose unique samples per subject per timepoint
# There are few samples per subject that are at the same timepoint.Calc
alpha_dt_unique <- alpha_dt %>% 
  group_by(Subject.ID) %>% 
  distinct(Timepoint.Calc, .keep_all = TRUE) %>% 
  ungroup()
# Make a boxplot at all time points:
pdf(paste(results_folder,paste0("Diversity_boxplot",'.pdf'),sep="/"),height =4 , width = 20, useDingbats = FALSE)

var <- div_index[1]
for(var in div_index){
  p_box_div <- ggplot(alpha_dt_unique, aes_string(x="Timepoint.Calc", y=var ,fill = "cohort_id_long",
                                                  color = "cohort_id_long")) + 
    geom_point(size = 2,alpha = 0.3)+
    geom_boxplot(alpha = 0.5,outlier.colour = NA)+
    scale_color_manual(name = "Cohort",values = coh.cols) +
    scale_fill_manual(name = "Cohort",values = coh.cols) +
    facet_wrap(~cohort_id_long,scales = "free_x",nrow = 1)+
    theme_bw()+
    xlab("Time")+
    theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
          axis.text.y=element_text(size=10,face="bold"),
          axis.title.y = element_text(size=10,face="bold"),
          axis.title.x = element_text(size=10,face="bold"))
  print(p_box_div)
}
dev.off()





# LME 
var <- div_index[3]
result_var_lme <-  list()
for(var in div_index){
  print(var)
  # LME comparing MHI pre vanco, during and early recovery
  result_lme <- list()
  coh <-  unique(alpha_dt$cohort_id_long)[1]
  for(coh in unique(alpha_dt$cohort_id_long)){
    
    lme_dt <- alpha_dt[alpha_dt$cohort_id_long == coh,]
    names(lme_dt)
    library("nlme")
    if(var == "richness"){
      fml <- as.formula( paste( "log","(",var,")", "~", paste(c("Timepoint.Calc"), collapse="+") ) )
    }else{
      fml <- as.formula( paste( var, "~", paste(c("Timepoint.Calc"), collapse="+") ) )
    }
    mod_bac <- lme(fml,random = ~1|Subject.ID, lme_dt)
    
    sum_mod <-  summary(mod_bac)
    sum_mod_dt <- data.frame(sum_mod$tTable)
    sum_mod_dt$coh <-  coh
    sum_mod_dt$Index <-  var
    sum_mod_dt$Var <-  rownames(sum_mod_dt)
    result_lme[[coh]] <- sum_mod_dt
  }
  
  lme_coh_dt <- do.call("rbind",result_lme)  
  lme_coh_dt$padj <- p.adjust(lme_coh_dt$p.value,method = "BH") 
  
  result_var_lme[[var]] <-  lme_coh_dt
}




final_lme_dt <-  do.call("rbind",result_var_lme)
require(openxlsx)
write.xlsx(result_var_lme, file = paste(results_folder,paste0("Diversity_lme_results_comp_baseline.xlsx"),sep="/"))

# Now, Beta diversity:
library(CoDaSeq)
# Metadata
meta <-  read.csv("../data/data_Mel/2021-05-14 VE303_species_table_for_beta_diversity_metadata.csv")
rownames(meta) <-  meta$sample_name
mat <-  read.csv("../data/data_Mel/2021-05-14 VE303_species_table_for_beta_diversity.csv")
rownames(mat) <-  mat$sample_name
mat$sample_name <-  NULL

fil_prop_dt <- codaSeq.filter(t(mat), min.prop=0.001, samples.by.row=FALSE)
# replace 0 values with an estimate based on Bayesian-multiplicative replacement 
imp_prop_dt <- cmultRepl(t(fil_prop_dt), method="CZM", label=0)
#Center Log-Ration Transform
clr_mat <- codaSeq.clr(imp_prop_dt)
#Now use phyloseq object to perform PCA analysis using CLR transformed data
otu_dat <- otu_table(t(clr_mat), taxa_are_rows = T)
phy_clr <-  phyloseq(otu_dat, sample_data(meta))
library(vegan)
#Well first perform a PCoA using euclidean distance from CLR transformed data
# Aitchison distance
set.seed(1057)
out.pcoa<- ordinate(phy_clr,  method = "PCoA", distance = "euclidean")
evals <- out.pcoa$values[,1]
phyloseq::plot_scree(out.pcoa) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

p <- plot_ordination(phy_clr, out.pcoa,axes = c(1,2)) 

data_pca <- p$data
data_pca$Day.from.Start
col_day_start <- c("Baseline"= "#7570b3","Vanco" = "#a6761d","Early recovery" = "#1b9e77",
                   "Late recovery" = "#d95f02","Early no vanco" = "#e7298a","Late no vanco" = "#66a61e")
data_pca[1,]
shape_time <-  c("Baseline"= 8,"Vanco" = 21,"Early recovery" = 22,
                 "Late recovery" = 24,"Early no vanco" = 10,"Late no vanco" = 13)

data_pca$Day.from.Start <-  factor(as.character(data_pca$Day.from.Start),levels = names(col_day_start))
data_pca$cohort_id_long <-  factor(data_pca$cohort_id_long)
data_pca$cohort_id_long <-  relevel(data_pca$cohort_id_long,ref = "Vanco")
# Fig 2 D
library(ggrepel)
library(lemon)
p_pca <- ggplot(data_pca, aes(x =Axis.1, y = Axis.2
))+ 
  geom_point(aes(fill = Day.from.Start),color = "black", pch = 21,stroke = 1,
             alpha = .7,size=3) + 
  theme_bw()+
  facet_rep_wrap(~cohort_id_long,repeat.tick.labels = T,nrow = 1)+
  theme(legend.text = element_text(size=15),
        legend.title=element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_fill_manual(name = "Time",values=col_day_start) +
  xlab(p$labels$x)+
  ylab (p$labels$y)


pdf(paste(results_folder,paste0("PCoA_plot",'.pdf'),sep="/"),height =4 , width = 25, useDingbats = FALSE)
print(p_pca)
dev.off()


library(ggrepel)
p_pca <- ggplot(data_pca, aes(x =Axis.1, y = Axis.2
))+ 
  geom_point(aes(color = cohort_id_long, 
                 fill = cohort_id_long,group=Subject.ID, shape = Day.from.Start),
             alpha = .5,size=2) + 
  theme_classic()+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_color_manual(name = "Cohort",values=coh.cols) +
  scale_fill_manual(name = "Cohort",values=coh.cols) +
  scale_shape_manual(name = "Time",values = shape_time)+
  xlab(p$labels$x)+
  ylab (p$labels$y)


pdf(paste(results_folder,paste0("PCoA_all_cohorts",'.pdf'),sep="/"),height =6 , width = 8, useDingbats = FALSE)
print(p_pca)
dev.off()


library(vegan)
library(pairwiseAdonis)
# For all cohorts 
# LME comparing MHI pre vanco, during and early recovery
result_lme <- list()
coh <-  unique(alpha_dt$cohort_id_long)[1]
p_manova_res_list <-  list()
for(coh in unique(alpha_dt$cohort_id_long)){
  
  phy_coh <-  subset_samples(phy_clr, cohort_id_long == coh)
  metadata <- as(sample_data(phy_coh), "data.frame")
  p_manova  <- pairwise.adonis(phyloseq::distance(phy_coh, method="euclidean") ,factors =  metadata$Day.from.Start)
  # Filter p_manova results and only include comparisons with Baseline:
  p_manova <-  data.frame(p_manova)
  p_manova$p.adjusted <- NULL
  p_manova$sig <-  NULL
  p_manova <- p_manova[grepl("Baseline",p_manova$pairs),]
  p_manova$Cohort <- coh
  p_manova_res_list[[coh]] <-  p_manova
  
}

p_manova_coh <-  do.call("rbind",p_manova_res_list)
p_manova_coh$padj <- p.adjust(p_manova_coh$p.value,method = "BH") 

write.csv(p_manova_coh,file = paste(results_folder,paste0("Beta_Div_pair_manova_Time_Cohorts.csv"),sep="/"))

