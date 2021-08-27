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
subDir <- paste0("VE303_recov/all_Cohorts")
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")
  
col_day_start <- c("Baseline"= "#7570b3","Vanco" = "#a6761d","Early recovery" = "#1b9e77",
                  "Late recovery" = "#d95f02","Early no vanco" = "#e7298a","Late no vanco" = "#66a61e")

# Now do the same thing but for all the cohorts
# Read different tax level data
tax_list <- list()
tax_lev <- c("phylum","class","family","genus")[2]
for(tax_lev in c("phylum","class","family","genus")){
  
  if(tax_lev == "phylum"){
    tax_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_phlyum_RA.csv")
  }else if(tax_lev == "class"){
    tax_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_class_RA.csv")
  }else if(tax_lev == "family"){
    tax_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_family_RA.csv")
  }else if(tax_lev == "genus"){
    tax_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_genus_RA.csv")
  }else{
    
  }
  
  tax_dt$Treatment <- tax_dt$Day.from.Start  
  tax_dt$cohort_id_long <-  factor(tax_dt$cohort_id_long, levels = unique(tax_dt$cohort_id_long))
  tax_dt$cohort_id_long <- relevel(tax_dt$cohort_id_long,ref = "Vanco")
  tax_dt$tax_name <- tax_dt[,2]
  
  
  library(ggplot2)
  library(ggthemes)
  # Filter low prevalent taxa from the analysis 
  sum_tax_dt <- tax_dt[tax_dt$Treatment %in% c("Baseline","Early recovery"),]
  sum_tax_dt <-  sum_tax_dt[sum_tax_dt$cohort_id_long != "Cohort 6",]
  
  tax.summary <- sum_tax_dt %>%
    group_by(tax_name) %>%
    summarise(
      sd = sd(rel_abund)/sqrt(n()),
      Mean_abun = 100*mean(rel_abund),
      Prevalence = sum(rel_abund> 10^-4)
    )%>% data.frame()
  tax.summary
  # Remove taxa that are less than 1 %
  tax.summary$prev_frac <-  tax.summary$Prevalence/length(unique(sum_tax_dt$sample_name))  
  fil_tax <-  tax.summary$tax_name[tax.summary$prev_frac>= 0.2]
  fil_tax_dt <-  sum_tax_dt[sum_tax_dt$tax_name %in% as.character(fil_tax),]  
  
  
  
  result_lme <- list()
  
  pdf(paste(results_folder,paste0("Rel_Abun_",tax_lev,"_early_recovery",'.pdf'),sep="/"),height = 4, width = 6, useDingbats = FALSE)
  
  taxon <-  unique(fil_tax_dt$tax_name)[5]
  
  for(taxon in  unique(fil_tax_dt$tax_name)){
    
    early_bac_dt <-  fil_tax_dt[fil_tax_dt$tax_name ==  taxon,]    
    
    
    
    box_dt <-  tax_dt[tax_dt$Treatment %in% c("Baseline","Early recovery","Early no vanco"),]
    box_dt <-  box_dt[box_dt$tax_name == taxon,]
    
    set.seed(100)
    p_box_tax <- ggplot(box_dt, aes(x=cohort_id_long, y = rel_abund*100 + 0.01,fill = Treatment )) + 
      geom_point(size = 2,position=position_jitterdodge(dodge.width=0.9),
                 color = "black", pch = 21, stroke = 0.5, alpha = 0.8) +
      geom_boxplot( outlier.colour=NA, position=position_dodge(width=0.9),alpha = 0.2)+
      
      scale_fill_manual(name = "Time", values = col_day_start )+
      scale_color_manual(name = "Time", values = col_day_start )+
      #geom_line( color = "#1b9e77")+
      theme_bw()+
      scale_y_continuous(trans='log10',breaks = c(0.01,0.1,1,10,100),
                         labels = c(0.01,0.1,1,10,100),
                         limits = c(0.01,100))+
      theme(axis.text.x = element_text(color = "grey20", size = 10, face = "bold",
                                       angle = 45, vjust = 1, hjust=1),
            axis.text.y = element_text(color = "grey20", size = 10, face = "bold"),  
            axis.title.x = element_text(color = "grey20", size = 14, face = "bold"),
            axis.title.y = element_text(color = "grey20", size = 14,face = "bold"))+
      xlab("Cohort")+
      ylab("100*Rel_abund + 0.01")+
      ggtitle(gsub(".*__","",taxon))
    
    # Use geom_line()+geom_pointrange()
    print(p_box_tax)
    
    
    
    # Loop it across different baseline
    
    basline_lme <- list()
    baseline_cohorts <- c("Vanco","Cohort 1","Cohort 2","Cohort 3")
    baseline <- baseline_cohorts[1]
    for(baseline in baseline_cohorts[1]){
      
      print(baseline)
      # Linear mixed effect model
      lme_dt <-  early_bac_dt
      #names(lme_dt)[8] <-  "VE_303"
      #names(lme_dt) <-  gsub("Abundance.","",names(lme_dt))
      lme_dt$t_Abundance  <- asin(sqrt(lme_dt$rel_abund ))
      lme_dt$cohort_id_long <- relevel(lme_dt$cohort_id_long,baseline)
      
      # Use geom_line()+geom_pointrange()
      #print(p_box_tax)
      
      #lme_dt <-  lme_dt[lme_dt$cohort_id_long %in% c("Vanco","Cohort 4","Cohort 5"),]
      #lme_dt <-  lme_dt[lme_dt$Treatment %in% c("Early recovery"),]
      
      lme_dt$Treatment <-  relevel(lme_dt$Treatment,ref = "Early recovery")
      library("nlme")
      fml <- as.formula( paste( "t_Abundance", "~", paste(c("Treatment*cohort_id_long "), collapse="+") ) )
      #fml <- as.formula( paste( "t_Abundance", "~", paste(c("cohort_id_long "), collapse="+") ) )
      mod_bac <- lme(fml,random = ~1|Subject.ID, lme_dt)
      sum_mod <-  summary(mod_bac)
      sum_mod_dt <- data.frame(sum_mod$tTable)
      sum_mod_dt$Var <-  rownames(sum_mod_dt)
      sum_mod_dt$tax_lev <- tax_lev
      sum_mod_dt$taxon <-  taxon
      sum_mod_dt$baseline <-  baseline
      # library(emmeans)
      # em_sum <- emmeans(mod_bac,specs = pairwise ~ cohort_id_long,adjust="none")
      # contrast_dt  <-  data.frame(em_sum$contrasts)
      # rownames(contrast_dt) <- paste0("Contrast ",contrast_dt$contrast)
      # contrast_dt$Phylum <-  phyla
      # contrast_dt$padj <- p.adjust(contrast_dt$p.value,method = "BH")
      basline_lme[[baseline]] <- sum_mod_dt
    }  
    
    lme_base_dt <-   do.call("rbind",basline_lme)
    #lme_base_dt$padj <- p.adjust(lme_base_dt$p.value,method = "BH")
    
    result_lme[[taxon]] <-  lme_base_dt
    
    # result_lme_dt <-  result_lme_dt[result_lme_dt$Var != "(Intercept)",]
    #result_lme_dt$padj <- p.adjust(result_lme_dt$p.value,method = "BH")
    # 
  }  
  dev.off() 
  result_lme_dt <-  do.call("rbind",result_lme)
  
  # Print out results for different baseline
  result_lme_dt <- result_lme_dt[result_lme_dt$baseline == "Vanco",]
  
  result_lme_dt$padj <- p.adjust(result_lme_dt$p.value,method = "BH")
  
  write.csv(result_lme_dt,paste(results_folder,paste0(tax_lev,"_lme_results_baseline_vanco",'.csv'),sep="/"))
  # Remove intercept
  result_lme_dt <-  result_lme_dt[result_lme_dt$Var != "(Intercept)",]
  # Also remove TreatmentBaseline (Detects the different between Baseline and Early recovery in Vanco)
  result_lme_dt <-  result_lme_dt[result_lme_dt$Var != "TreatmentBaseline",]
  
  # Remove interaction
  result_lme_dt <-  result_lme_dt[!grepl(":",result_lme_dt$Var),]
  
  # Filter 
  result_lme_dt_sig <- result_lme_dt[result_lme_dt$padj < 0.05,]
  
  
  # Plot the significant taxa 
  pdf(paste(results_folder,paste0("Rel_Abun_",tax_lev,"_early_recovery_sig",'.pdf'),sep="/"),height = 4, width = 6, useDingbats = FALSE)
  for(taxon in  unique(result_lme_dt_sig$taxon)){
    
    #sub_bac_dt <-  fil_tax_dt[fil_tax_dt$tax_name ==  taxon,]    
    
    sub_bac_dt <-  tax_dt[tax_dt$Treatment %in% c("Baseline","Early recovery","Early no vanco"),]
    sub_bac_dt <-  sub_bac_dt[sub_bac_dt$tax_name == taxon,]
    
    p_box_tax <- ggplot(sub_bac_dt,  aes(x=cohort_id_long, y = rel_abund*100 + 0.01,fill = Treatment )) + 
      geom_point(size = 2,position=position_jitterdodge(dodge.width=0.9),
                 color = "black", pch = 21, stroke = 0.5, alpha = 0.8) +
      geom_boxplot( outlier.colour=NA, position=position_dodge(width=0.9),alpha = 0.2)+
      
      scale_fill_manual(name = "Time", values = col_day_start )+
      scale_color_manual(name = "Time", values = col_day_start )+
      #geom_line( color = "#1b9e77")+
      theme_bw()+
      scale_y_continuous(trans='log10',breaks = c(0.01,0.1,1,10,100),
                         labels = c(0.01,0.1,1,10,100),
                         limits = c(0.01,100))+
      theme(axis.text.x = element_text(color = "grey20", size = 10, face = "bold",
                                       angle = 45, vjust = 1, hjust=1),
            axis.text.y = element_text(color = "grey20", size = 10, face = "bold"),  
            axis.title.x = element_text(color = "grey20", size = 14, face = "bold"),
            axis.title.y = element_text(color = "grey20", size = 14,face = "bold"))+
      xlab("Cohort")+
      ylab("100*Rel_abund + 0.01")+
      ggtitle(gsub(".*__","",taxon))
    
    # Use geom_line()+geom_pointrange()
    print(p_box_tax)
    
    
  }
  dev.off()
  
  
}





# Use MHI to do the similar analysis
mhi_dt <-  read.csv("../data/data_Mel/2021-04-26 VE303_Ph1_MHI_values.csv")

# MHI index
col_day_start <- c("Baseline"= "#7570b3","Vanco" = "#a6761d","Early recovery" = "#1b9e77",
                   "Late recovery" = "#d95f02","Early no vanco" = "#e7298a","Late no vanco" = "#66a61e")
mhi_dt$Treatment <- mhi_dt$Day.from.Start  
mhi_dt$Treatment <-   factor(mhi_dt$Treatment , levels = names(col_day_start))

# Plot FC for Bacteroidetes/Proteobacteria

# Order by time and cohort
mhi_dt <-  mhi_dt[order(mhi_dt$Subject.Cohort, mhi_dt$Treatment),]

mhi_dt$sample_name <- as.character(mhi_dt$sample_name)
mhi_dt$sample_name <- factor(mhi_dt$sample_name, levels = unique(mhi_dt$sample_name))


mhi_dt$cohort_id_long <-  factor(mhi_dt$cohort_id_long, levels = unique(mhi_dt$cohort_id_long))
mhi_dt$cohort_id_long <- relevel(mhi_dt$cohort_id_long,ref = "Vanco")


p_MHI <- ggplot(mhi_dt) +
  geom_bar(aes(x = factor(sample_name),y=MHI, fill = Treatment),size=1,stat = "identity")+
  scale_fill_manual(name = "Day.from.Start", values = col_day_start )+
  theme_base()+
  facet_grid(~cohort_id_long,scales = "free", space = "free")+
  scale_y_continuous(trans='log10')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  xlab("")+
  ylab("MHI")

print(p_MHI)

#print(p_FC_bai)
pdf(paste(results_folder,paste0("MHI_timeplot",'.pdf'),sep="/"),height =6 , width = 40, useDingbats = FALSE)
print(p_MHI)
dev.off()


library(dplyr)
# Group by plot just for  Eary Recovery
early_mhi_dt <-  mhi_dt[mhi_dt$Treatment %in% c("Baseline","Early recovery","Early no vanco"),]
#early_mhi_dt <-  early_mhi_dt[early_mhi_dt$cohort_id_long !=  c("Cohort 6"),]


mhi.summary <- early_mhi_dt %>%
  group_by(cohort_id_long,Treatment) %>%
  summarise(
    sd = sd(MHI)/sqrt(n()),
    MHI = mean(MHI)
  )%>% data.frame()
mhi.summary

pdf(paste(results_folder,paste0("MHI_box_early_recovery",'.pdf'),sep="/"),height =4 , width = 5, useDingbats = FALSE)

# Use geom_pointrange
set.seed(100)
p_box_mhi <- ggplot(early_mhi_dt, aes(x=cohort_id_long, y=MHI ,fill = Treatment)) + 
  #geom_boxplot()+
  #geom_point(size = 2,aes(color = Treatment))+
  geom_point(position=position_jitterdodge(dodge.width=0.9),
             aes(fill = Treatment),color = "black", pch = 21, stroke = 0.5, alpha = 0.8) +
  geom_boxplot( outlier.colour=NA, position=position_dodge(width=0.9),alpha = 0.2)+
  
  scale_fill_manual(name = "Time", values = col_day_start )+
  scale_color_manual(name = "Time", values = col_day_start )+
  #geom_line( color = "#1b9e77")+
  scale_y_continuous(trans='log10')+
  theme_bw()+
  #theme_base()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, face = "bold",
                                   angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(color = "grey20", size = 10, face = "bold"),  
        axis.title.x = element_text(color = "grey20", size = 14, face = "bold"),
        axis.title.y = element_text(color = "grey20", size = 14,face = "bold"))+
  xlab("Cohort")

# Use geom_line()+geom_pointrange()
print(p_box_mhi)

dev.off()

# LME comparing MHI pre vanco, during and early recovery
# lme_dt <- mhi_dt[mhi_dt$Treatment %in% c("Baseline","Early recovery"),]
# lme_dt <- lme_dt[lme_dt$cohort_id_long != "Cohort 6",]

lme_dt <-   early_mhi_dt[early_mhi_dt$cohort_id_long !=  c("Cohort 6"),]
names(lme_dt)
lme_dt$Treatment <-  relevel(lme_dt$Treatment,ref = "Early recovery")

library("nlme")
fml <- as.formula( paste( "MHI_log10", "~", paste(c("Treatment*cohort_id_long "), collapse="+") ) )
mod_bac <- lme(fml,random = ~1|Subject.ID, lme_dt)
sum_mod <-  summary(mod_bac)
sum_mod_dt <- data.frame(sum_mod$tTable)
write.csv(sum_mod_dt,paste(results_folder,paste0("MHI_lme_results",'.csv'),sep="/"))


