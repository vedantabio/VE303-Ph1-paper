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
subDir <- paste0("LME_ve303_colonization_model/",gsub("-","_",Sys.Date()))
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

# Import dataset
# Add appropriated dated folder 
input_data <-  readRDS("../results/Input_matrix/2021_12_20/input_data_rf.RDS")[["Species"]]
abun_data <- readRDS("../results/Input_matrix/2021_12_20/ve303_abun_rf.RDS")[["Species"]]
phy_orig <-  readRDS("../results/Input_matrix/2021_12_20/phy_list.RDS")[["Species"]]

phy_list <-  subset_samples(phy_orig, !cohort_id_long %in% c("Vanco","Cohort 6") )
phy_list <-  prune_taxa(taxa_sums(phy_list)>0, phy_list)

phy_pre <- subset_samples(phy_list,Day.from.Start == "Baseline")
phy_vanco <- subset_samples(phy_list,Day.from.Start == "Vanco")
sm_dt <-  data.frame(sample_data(phy_vanco))
# Post Vanco samples for early recovery only:
phy_post <- subset_samples(phy_vanco,Collection.Day.Norm %in% c(5,6))
# Merge baseline samples and pre-VE303 samples 
phy_mer <-  merge_phyloseq(phy_pre,phy_post)
phy_mer <-  prune_taxa(taxa_sums(phy_mer)>0, phy_mer)

# For post ve303 we increase the time points post Vanco
# Subset samples less than 40 days
phy_post_ve303 <- subset_samples(phy_list,Day.from.Start %in% c("Early recovery","Late recovery"))
phy_post_ve303 <- subset_samples(phy_post_ve303,Collection.Day.Norm <= 40)

phy_pre_post_ve303 <-   merge_phyloseq(phy_pre,phy_post_ve303)
phy_pre_post_ve303 <-  prune_taxa(taxa_sums(phy_pre_post_ve303)>0, phy_pre_post_ve303)

# Keep only prevalent species in baseline:
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

prev_frac<-0.05
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
most_prevalent_taxa<-unique(rownames(prevdf_fil))
phy_mer = prune_taxa(most_prevalent_taxa, phy_mer)


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


mean_dt$Post <-  mean_dt$Post + 1e-5 
mean_dt$Pre <-  mean_dt$Pre + 1e-5 
# mean_dt[mean_dt$Post == 0,"Post"] <- 1e-5
# mean_dt[mean_dt$Pre == 0,"Pre"] <- 1e-5

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


# Select only samples that are post VE303
phy_sel <- prune_samples(abun_data$sample_name,phy_orig)
taxa_names(phy_sel) <-  make.names(taxa_names(phy_sel))

# Log fold change 
lfc_dt <- mean_LFC[,grep("logFC",names(mean_LFC))]
names(lfc_dt) <-  gsub("logFC.","",names(lfc_dt))
lfc_dt <-  cbind(Subject.ID = mean_LFC$Subject.ID,lfc_dt)
lfc_dt_long <- lfc_dt %>%
  pivot_longer(!Subject.ID,names_to = "Taxa",values_to = "LFC")
unique(lfc_dt_long$Taxa)
# Select bugs that have LogFC
phy_sel <-  prune_taxa(names(lfc_dt),phy_sel) 

# Run LME
# Perform LME to identify species that relate to VE303 colonization
result_lme <-  list()
melt_dt <-  psmelt(phy_sel)
# Now merge the data with ve303 abundance 
comb_dt <-  melt_dt %>%
  left_join(abun_data)
comb_dt <-  comb_dt[,c("Subject.ID","OTU","Sample","Abundance","Collection.Day.Norm","Day.from.Start","cohort_id_long","sum_rel_abun")]
taxa <- taxa_names(phy_sel)[4]

for (taxa in  taxa_names(phy_sel)){
  print(taxa)
  sum_mod_dt <- tryCatch({
    #phy_tax <- prune_taxa(taxa,phy_mer)
    phy_m <-  comb_dt[comb_dt$OTU %in% taxa,]
    lfc_dt_sub <-  lfc_dt_long[lfc_dt_long$Taxa %in% taxa,]
    # Select only important covariates
    phy_m <- phy_m %>%
      left_join(lfc_dt_sub)
    phy_m <-  phy_m[!is.na(phy_m$LFC),]
    
    phy_m$cohort_id_long <- factor(phy_m$cohort_id_long)
    phy_m$t_Abundance  <- asin(sqrt(phy_m$Abundance ))
    phy_m$t_sum_rel_abun  <- asin(sqrt(phy_m$sum_rel_abun ))
    phy_m$time <-  as.numeric(scale(phy_m$Collection.Day.Norm))
    library("nlme")
    mod_bac <- lme(t_sum_rel_abun ~ t_Abundance +  LFC + time,random = ~1 |cohort_id_long/Subject.ID, phy_m)
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


result_lme_w <-  result_lme_dt %>%
  select(Value,pval_adj,Bac,Var)%>%
  pivot_wider(names_from = Var, values_from = c(Value,pval_adj))

write.csv(result_lme_w,paste0(results_folder,"/LME_tab.csv"))
### LFC
result_lme_lfc <-  result_lme_w[result_lme_w$pval_adj_LFC <= 0.1,]

### Abundances
result_lme_abun <-  result_lme_w[result_lme_w$pval_adj_t_Abundance <= 0.05,]

# Make a heatmap of species with logFC corresponding to  VE303-total abundance:
abun_data <- abun_data[abun_data$Subject.ID %in% lfc_dt$Subject.ID,]

mat_logFC <-  lfc_dt[,-1]
mat_logFC <-  mat_logFC[,colnames(mat_logFC) %in% result_lme_lfc$Bac]
#Summarise VE303 abundance per subject
mean_ve_abun <- abun_data %>%
  group_by(Subject.ID)%>%
  summarise(Mean_abun = mean(sum_rel_abun))
mean_ve_abun <-  mean_ve_abun[order(mean_ve_abun$Mean_abun),]


mat_ht <-  t(mat_logFC)
mat_ht <- mat_ht[,match(mean_ve_abun$Subject.ID,lfc_dt$Subject.ID)]


# Row annotation:
row_ann_dt <-  result_lme_lfc
row_ann_dt <-  row_ann_dt[match(rownames(mat_ht),row_ann_dt$Bac  ),]


library(gtools)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(data.table)
library(hues)
library(scales)
#splitcols <- c(rep("Baseline",ncol(mat_baseline)), rep("LogFC",ncol(mat_logFC)),"VE303")
col_b <-  colorRamp2( seq(-0.04,0.04,0.02), colorRampPalette(c(muted("blue"),"white",muted("red")))(5) )
ha_left = HeatmapAnnotation(Coeff = row_ann_dt$Value_LFC,
                            col = list(Coeff = col_b), which = "row")

library(scico)
library(wacolors)
col_ve303 <- colorRamp2(rev(c(seq(-0.5,-6,-0.5))),viridis::viridis(12))

ha_column = HeatmapAnnotation(Total_VE303 = log10(mean_ve_abun$Mean_abun + 1e-5),
                              col=list(Total_VE303 = col_ve303),
                              annotation_legend_param = list(title = "Total_VE303",
                                                             col_fun = col_ve303,
                                                             at = c(0,-1:-6),
                                                             legend_height = unit(7, "cm")),
                              annotation_name_side = "left")


min = -12.7
max = 6.27
gap = 4

make_seq <- function(min, max, gap){
  # Make it even
  min_r <- ifelse(round(min)%%2 == 0, round(min),round(min) -1)
  max_r <- ifelse(round(max)%%2 == 0, round(max),round(max) +1)
  max_val <-  max(abs(min_r),abs(max_r))
  
  inc <-max_val/gap
  seq_tot <-  seq(min_r-inc, max_r+inc,by = inc)
  
  return(seq_tot)
}

vec_mat <-  make_seq(min(mat_ht),max(mat_ht),4)
col_mat <- colorRamp2(c(min(mat_ht), 0, max(mat_ht)), c("blue", "white", "red"))
#col_mat <- colorRamp2(seq(-4,4,by = 2),colorRampPalette(c("#7f0000", "white","#00007f"))(6))
library(ComplexHeatmap)
ht  =  Heatmap(mat_ht,name = "LogFC",
               # row_split = rsplit[r_order],
               gap = unit(2, "mm"),
               #column_split = splitcols,
               #col = col_mat,
               top_annotation = ha_column,
               left_annotation = ha_left,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_title_rot = 0,
               column_title = "Cohort 1-5",
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = T,
               cluster_row_slices = F,
               cluster_columns = F,
               border = T,
               show_row_dend = F,
               heatmap_legend_param = list(title = "log2FC",
                                           col_fun = col_mat,
                                           at = vec_mat,
                                           legend_height = unit(6, "cm")),
               
               row_names_max_width = max_text_width(rownames(mat_ht),
                                                    gp = gpar(fontsize = 14)),
               column_names_max_height = max_text_width(colnames(mat_ht),
                                                        gp = gpar(fontsize = 12)))


pdf(paste0(results_folder,paste0("/Heatmap_Sig_LogFC_Coh1_5.pdf")),height =4, width = 9)
draw(ht)
dev.off()


# Lets just look at the abundance of the bugs post VE303
# Do the similar with baseline abundance:
# Make a heatmap of species with logFC corresponding to  VE303-total abundance:
phy_abun_sel <- phy_post_ve303
taxa_names(phy_abun_sel) <-  make.names(taxa_names(phy_abun_sel))  
phy_abun_sel <-   prune_taxa(result_lme_abun$Bac,phy_abun_sel)



mat_ht <- data.frame(otu_table(phy_abun_sel))
mat_ht <-  log10(mat_ht + 1e-5)

# Add ve303 abundance
#Arrange VE303 abundance per subject
mean_ve_abun <-  abun_data[order(abun_data$sum_rel_abun),]


mat_ht <- mat_ht[,match(mean_ve_abun$sample_name,colnames(mat_ht))]

library(gtools)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(data.table)
library(hues)
library(scales)
#splitcols <- c(rep("Baseline",ncol(mat_baseline)), rep("LogFC",ncol(mat_logFC)),"VE303")
col_b <- c("0" = "grey","-1"= muted("blue"),"1"= muted("red"))
# ha_left = HeatmapAnnotation(Relation = row_ann_dt$Dir,
#                             col = list(Relation = col_b), which = "row")
# 
library(scico)
library(wesanderson)
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


tax_dt_sel <-  data.frame(tax_table(phy_abun_sel)@.Data)
tax_dt_sel <-  tax_dt_sel[match( rownames(mat_ht),rownames(tax_dt_sel)),]

# Row annotation:
row_ann_dt <-  result_lme_abun
row_ann_dt <-  row_ann_dt[match(rownames(mat_ht),row_ann_dt$Bac  ),]


col_coef <-  colorRampPalette(c(muted("blue"),"white",muted("red")))(7)
col_b <-  colorRamp2(c(seq(-6,6,2)),col_coef  )
ha_left = HeatmapAnnotation(Coeff = row_ann_dt$Value_t_Abundance,
                            col = list(Coeff = col_b),
                            annotation_legend_param = list (Coeff = list(title = "Coeff_LME", at = seq(-6,6,2))),
                            which = "row")


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


pdf(paste0(results_folder,paste0("/Heatmap_RA_Coh1_5.pdf")),height = 6, width = 15)
draw(ht)
dev.off()  
