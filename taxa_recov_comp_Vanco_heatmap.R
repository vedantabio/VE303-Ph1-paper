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

# Includes all cohorts:
mainDir <- "../results"
subDir <- paste0("VE303_recov/Heatmap_sig")
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")


ve303_abun <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 ve303_total_absolute_abundance.csv")
ve303_abun <-  unique(ve303_abun[,c("sample_name","absolute_abund","cohort_id_long","Day.from.Start")])
names(ve303_abun)[2] <- "VE303_tot_abs_abund" 

col_day_start <- c("Baseline"= "#7570b3","Vanco" = "#a6761d","Early recovery" = "#1b9e77",
                   "Late recovery" = "#d95f02","Early no vanco" = "#e7298a","Late no vanco" = "#66a61e")

# Read different tax level data
tax_list <- list()
tax_lev <- c("phylum","class","family","genus")[3]

tax_list <-  list()
sig_coh_list <- list()

for(tax_lev in c("phylum","class","family","genus")){
  
  # if(tax_lev == "phylum"){
  #   tax_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_phlyum_RA.csv")
  #   result_lme_dt <-  read.csv("../results/VE303_recov/all_Cohorts/2021_05_13/phylum_lme_results_baseline_vanco.csv")
  # }else if(tax_lev == "class"){
  #   tax_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_class_RA.csv")
  #   result_lme_dt <-  read.csv("../results/VE303_recov/all_Cohorts/2021_05_13/class_lme_results_baseline_vanco.csv")
  # }else if(tax_lev == "family"){
  #   tax_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_family_RA.csv")
  #   result_lme_dt <-  read.csv("../results/VE303_recov/all_Cohorts/2021_05_13/family_lme_results_baseline_vanco.csv")
  # }else if(tax_lev == "genus"){
  #   tax_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_genus_RA.csv")
  #   result_lme_dt <-  read.csv("../results/VE303_recov/all_Cohorts/2021_05_13/genus_lme_results_baseline_vanco.csv")
  # }else{
  #   
  # }
  if(tax_lev == "phylum"){
    tax_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_phlyum_RA.csv")
    result_lme_dt <-  read.csv("../results/VE303_recov/all_Cohorts/phylum_lme_results_baseline_vanco.csv")
  }else if(tax_lev == "class"){
    tax_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_class_RA.csv")
    result_lme_dt <-  read.csv("../results/VE303_recov/all_Cohorts/class_lme_results_baseline_vanco.csv")
  }else if(tax_lev == "family"){
    tax_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_family_RA.csv")
    result_lme_dt <-  read.csv("../results/VE303_recov/all_Cohorts/family_lme_results_baseline_vanco.csv")
  }else if(tax_lev == "genus"){
    tax_dt <-  read.csv("../data/data_Mel/Taxon_data_5_12/2021-05-12 VE303_Ph1_genus_RA.csv")
    result_lme_dt <-  read.csv("../results/VE303_recov/all_Cohorts/genus_lme_results_baseline_vanco.csv")
  }else{
    
  }
  
  
  #tax_dt <- tax_dt[tax_dt$cohort_id_long %in% c("Vanco","Cohort 4","Cohort 5"),]
  
  tax_dt$Treatment <- tax_dt$Day.from.Start  
  tax_dt$cohort_id_long <-  factor(tax_dt$cohort_id_long, levels = unique(tax_dt$cohort_id_long))
  tax_dt$cohort_id_long <- relevel(tax_dt$cohort_id_long,ref = "Vanco")
  tax_dt$tax_name <- tax_dt[,2]
  
  library(ggplot2)
  library(ggthemes)
  # Filter significant taxa
  # Remove intercept
  result_lme_dt <-  result_lme_dt[result_lme_dt$Var != "(Intercept)",]
  # Also remove TreatmentBaseline (Detects the different between Baseline and Early recovery in Vanco)
  result_lme_dt <-  result_lme_dt[result_lme_dt$Var != "TreatmentBaseline",]
  
  # Remove interaction
  result_lme_dt <-  result_lme_dt[!grepl(":",result_lme_dt$Var),]
  
  # Filter 
  result_lme_dt_sig <- result_lme_dt[result_lme_dt$padj < 0.05,]
  
  
  sig_coh_dt <-  unique(result_lme_dt_sig[,c("taxon","Var","Value")])
  sig_coh_list[[tax_lev]] <- sig_coh_dt
  
  # Filter low prevalent taxa from the analysis 
  sum_tax_dt <- tax_dt[tax_dt$Treatment %in% c("Baseline","Early recovery","Early no vanco"),]
  #sum_tax_dt <-  sum_tax_dt[sum_tax_dt$cohort_id_long != "Cohort 6",]
  
  sig_tax_dt <-  sum_tax_dt[sum_tax_dt$tax_name %in% unique(result_lme_dt_sig$taxon),]  
  
  # Only select necessary columns for heatmap:
  sig_tax_dt <- sig_tax_dt[, c("sample_name","tax_name","plot_order","rel_abund","cohort_id_long","Treatment")]
  
  #sig_tax_dt$level <-  tax_lev
  tax_list[[tax_lev]] <-  sig_tax_dt
  
}  

final_sig_dt <- do.call("rbind",tax_list) 


row_ann_dt <-  do.call("rbind",sig_coh_list)
names(row_ann_dt)[2] <- "Cohort"
row_ann_dt$Cohort <-  gsub("cohort_id_long","",row_ann_dt$Cohort)
row_ann_dt$taxon <- gsub(".*__","",row_ann_dt$taxon)
row_ann_dt$Value <- ifelse(row_ann_dt$Value < 0, -1,1)

# Widen the dataframe into cohort 
library(tidyr)

row_ann_dt_w <- row_ann_dt %>%
  pivot_wider(names_from = Cohort,values_from = Value) %>%
   data.frame()
names(row_ann_dt_w) <- gsub("\\."," ",names(row_ann_dt_w))

row_ann_dt_w <-  row_ann_dt_w[,c("taxon","Cohort 1", "Cohort 2", "Cohort 3", "Cohort 4", "Cohort 5")]



sig_tax_w <- final_sig_dt %>%
  pivot_wider(
    names_from = tax_name,
    values_from = c(rel_abund)
  )
library(dplyr)
sig_tax_w <-  sig_tax_w %>% arrange(plot_order) %>% data.frame()

mat_ht <-  t(sig_tax_w[,5:ncol(sig_tax_w)]) 

metadata <- sig_tax_w[,1:4]

rsplit  <-  rownames(mat_ht)
rsplit <-  rep("Phylum",length(rsplit) )
rsplit[grep("f_",rownames(mat_ht))] <- "Family"
rsplit[grep("g_",rownames(mat_ht))] <- "Genus"
rsplit[grep("c_",rownames(mat_ht))] <- "Class"

rsplit <- factor(rsplit,levels = c("Genus","Family","Class","Phylum")) 


mat2 <- mat_ht
#rownames(mat2) <-  rownames(dt_mat)
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
# Order matrix 
#mat2 <- mat2[match(match_tax_dt$otu,rownames(mat2)),]
# Order columns based on CFU
dist_colors <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6",
                 "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3","#808000","#ffd8b1",
                 "#000080")
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))



library(scico)
col_cfu <- scico(5, palette = 'lajolla')
col_time <- scico(2,palette = "devon")

ant_col <-     data.frame(Cohort= metadata$cohort_id_long,
                          Time = metadata$Treatment)


col_func_time <- c("Baseline"= "#7570b3","Early recovery" = "#1b9e77","Early no vanco" = "#e7298a")

ha_column = HeatmapAnnotation(Time = ant_col$Time,
                              col=list(Time = col_func_time))
splitcols <-ant_col$Cohort

library(scico)
col_mat <- scico(15, palette = 'lajolla')

#mat4 <-sqrt(sqrt(mat2))
mat4 <- log10(mat2 + 1e-5)
rownames(mat4) <-  gsub(".*__","",rownames(mat4))


row_ann_dt_w <- row_ann_dt_w[match(row_ann_dt_w$taxon,rownames(mat4)),]
col_b <- c("0" = "grey","1"= "blue","-1"= "red")
ha_left = HeatmapAnnotation(Cohort_1 = row_ann_dt_w$`Cohort 1`,
                            Cohort_2 = row_ann_dt_w$`Cohort 2`,
                            Cohort_3 = row_ann_dt_w$`Cohort 3`,
                            Cohort_4 = row_ann_dt_w$`Cohort 4`,
                            Cohort_5 = row_ann_dt_w$`Cohort 5`,
                            col = list(Cohort_1 = col_b,
                                       Cohort_2 = col_b,
                                       Cohort_3 = col_b,
                                       Cohort_4 = col_b,
                                       Cohort_5 = col_b), which = "row")



library(ComplexHeatmap)
# ht  =  Heatmap(mat4,name = "sqrt(sqrt(Rel_Abun))",
#                row_split = rsplit,
#                gap = unit(2, "mm"),
#                column_split = splitcols,
#                top_annotation = ha_column,
#                left_annotation = ha_left,
#                #col = jet.colors(5),
#                col = col_mat,
#                row_names_side = "left",
#                row_gap = unit(2, "mm"),
#                row_title_gp = gpar(fontsize = 10,fontface = "bold"),
#                row_title_rot = 0,
#                column_title_gp = gpar(fontsize = 10,fontface = "bold"),
#                column_title_rot = 0,
#                cluster_rows = T,
#                cluster_columns = F,
#                border = T,
#                show_row_dend = F,
#                row_names_max_width = max_text_width(rownames(mat4),
#                                                     gp = gpar(fontsize = 12)),
#                column_names_max_height = max_text_width(colnames(mat4),
#                                                         gp = gpar(fontsize = 12)))
# pdf(paste0(results_folder,"/Heatmap_Sig_tax_Early_recov_All_Cohorts.pdf"),height = 6, width = 25)
# draw(ht)
# dev.off()

library(ComplexHeatmap)
ht  =  Heatmap(mat4,name = "Log10(Rel_Abun)",
               row_split = rsplit,
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
               cluster_rows = T,
               #cluster_row_slices = T,
               cluster_columns = F,
               border = T,
               show_row_dend = F,
               row_names_max_width = max_text_width(rownames(mat4),
                                                    gp = gpar(fontsize = 12)),
               column_names_max_height = max_text_width(colnames(mat4),
                                                        gp = gpar(fontsize = 12)))

r_order <- as.numeric(unlist(row_order(ht)))

mat4 <-  mat4[r_order,]
row_ann_dt_w <- row_ann_dt_w[match(rownames(mat4),row_ann_dt_w$taxon),]
col_b <- c("0" = "grey","1"= "blue","-1"= "red")
ha_left = HeatmapAnnotation(Cohort_1 = row_ann_dt_w$`Cohort 1`,
                            Cohort_2 = row_ann_dt_w$`Cohort 2`,
                            Cohort_3 = row_ann_dt_w$`Cohort 3`,
                            Cohort_4 = row_ann_dt_w$`Cohort 4`,
                            Cohort_5 = row_ann_dt_w$`Cohort 5`,
                            col = list(Cohort_1 = col_b,
                                       Cohort_2 = col_b,
                                       Cohort_3 = col_b,
                                       Cohort_4 = col_b,
                                       Cohort_5 = col_b), which = "row")




ht  =  Heatmap(mat4,name = "Log10(Rel_Abun)",
               row_split = rsplit[r_order],
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
               show_row_dend = F,
               row_names_max_width = max_text_width(rownames(mat4),
                                                    gp = gpar(fontsize = 12)),
               column_names_max_height = max_text_width(colnames(mat4),
                                                        gp = gpar(fontsize = 12)))


pdf(paste0(results_folder,"/Heatmap_Sig_tax_Early_recov_All_Cohorts.pdf"),height = 6, width = 25)
draw(ht)
dev.off()



# Now only display species that are Significant in Cohort 4 and 5

final_sig_dt <- do.call("rbind",tax_list) 

row_ann_dt <-  do.call("rbind",sig_coh_list)
names(row_ann_dt)[2] <- "Cohort"
row_ann_dt$Cohort <-  gsub("cohort_id_long","",row_ann_dt$Cohort)
row_ann_dt$taxon <- gsub(".*__","",row_ann_dt$taxon)
row_ann_dt$Value <- ifelse(row_ann_dt$Value < 0, -1,1)
row_ann_dt <-  row_ann_dt[row_ann_dt$Cohort %in% c("Cohort 4", "Cohort 5"),]

# Only select Vanco Cohort 4 and 5

# Widen the dataframe into cohort 
#library(tidyr)

row_ann_dt_w <- row_ann_dt %>%
  pivot_wider(names_from = Cohort,values_from = Value) %>%
  data.frame()
names(row_ann_dt_w) <- gsub("\\."," ",names(row_ann_dt_w))

row_ann_dt_w <-  row_ann_dt_w[,c("taxon", "Cohort 4", "Cohort 5")]


# Select only Vanco Cohort 4, 5 and 6
final_sig_dt <-  final_sig_dt[final_sig_dt$cohort_id_long %in% c("Vanco","Cohort 4", "Cohort 5","Cohort 6"),]
#final_sig_dt$tax_name_clean <- gsub(".*__","",final_sig_dt$tax_name)
final_sig_dt <-  final_sig_dt[gsub(".*__","",final_sig_dt$tax_name) %in% row_ann_dt_w$taxon,]

sig_tax_w <- final_sig_dt %>%
           pivot_wider(
                names_from = tax_name,
               values_from = c(rel_abund))

sig_tax_w <-  sig_tax_w %>% arrange(plot_order) %>% data.frame()

mat_ht <-  t(sig_tax_w[,5:ncol(sig_tax_w)]) 

metadata <- sig_tax_w[,1:4]

rsplit  <-  rownames(mat_ht)
rsplit <-  rep("Phylum",length(rsplit) )
rsplit[grep("f_",rownames(mat_ht))] <- "Family"
rsplit[grep("g_",rownames(mat_ht))] <- "Genus"
rsplit[grep("c_",rownames(mat_ht))] <- "Class"

rsplit <- factor(rsplit,levels = c("Genus","Family","Class","Phylum")) 


mat2 <- mat_ht
#rownames(mat2) <-  rownames(dt_mat)
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
# Order matrix 
#mat2 <- mat2[match(match_tax_dt$otu,rownames(mat2)),]
# Order columns based on CFU
dist_colors <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6",
                 "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3","#808000","#ffd8b1",
                 "#000080")
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))



library(scico)
col_cfu <- scico(5, palette = 'lajolla')
col_time <- scico(2,palette = "devon")

ant_col <-     data.frame(Cohort= metadata$cohort_id_long,
                          Time = metadata$Treatment)


col_func_time <- c("Baseline"= "#7570b3","Early recovery" = "#1b9e77","Early no vanco" = "#e7298a")

ha_column = HeatmapAnnotation(Time = ant_col$Time,
                              col=list(Time = col_func_time))
splitcols <-ant_col$Cohort

library(scico)
col_mat <- scico(15, palette = 'lajolla')

#mat4 <-sqrt(sqrt(mat2))
mat4 <- log10(mat2 + 1e-5)
rownames(mat4) <-  gsub(".*__","",rownames(mat4))


row_ann_dt_w <- row_ann_dt_w[match(row_ann_dt_w$taxon,rownames(mat4)),]
col_b <- c("0" = "grey","1"= "blue","-1"= "red")
ha_left = HeatmapAnnotation(Cohort_1 = row_ann_dt_w$`Cohort 1`,
                            Cohort_2 = row_ann_dt_w$`Cohort 2`,
                            Cohort_3 = row_ann_dt_w$`Cohort 3`,
                            Cohort_4 = row_ann_dt_w$`Cohort 4`,
                            Cohort_5 = row_ann_dt_w$`Cohort 5`,
                            col = list(Cohort_1 = col_b,
                                       Cohort_2 = col_b,
                                       Cohort_3 = col_b,
                                       Cohort_4 = col_b,
                                       Cohort_5 = col_b), which = "row")



library(ComplexHeatmap)
# ht  =  Heatmap(mat4,name = "sqrt(sqrt(Rel_Abun))",
#                row_split = rsplit,
#                gap = unit(2, "mm"),
#                column_split = splitcols,
#                top_annotation = ha_column,
#                left_annotation = ha_left,
#                #col = jet.colors(5),
#                col = col_mat,
#                row_names_side = "left",
#                row_gap = unit(2, "mm"),
#                row_title_gp = gpar(fontsize = 10,fontface = "bold"),
#                row_title_rot = 0,
#                column_title_gp = gpar(fontsize = 10,fontface = "bold"),
#                column_title_rot = 0,
#                cluster_rows = T,
#                cluster_columns = F,
#                border = T,
#                show_row_dend = F,
#                row_names_max_width = max_text_width(rownames(mat4),
#                                                     gp = gpar(fontsize = 12)),
#                column_names_max_height = max_text_width(colnames(mat4),
#                                                         gp = gpar(fontsize = 12)))
# pdf(paste0(results_folder,"/Heatmap_Sig_tax_Early_recov_All_Cohorts.pdf"),height = 6, width = 25)
# draw(ht)
# dev.off()

library(ComplexHeatmap)
ht  =  Heatmap(mat4,name = "Log10(Rel_Abun)",
               row_split = rsplit,
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
               cluster_rows = T,
               #cluster_row_slices = T,
               cluster_columns = F,
               border = T,
               show_row_dend = F,
               row_names_max_width = max_text_width(rownames(mat4),
                                                    gp = gpar(fontsize = 12)),
               column_names_max_height = max_text_width(colnames(mat4),
                                                        gp = gpar(fontsize = 12)))

r_order <- as.numeric(unlist(row_order(ht)))

mat4 <-  mat4[r_order,]
row_ann_dt_w <- row_ann_dt_w[match(rownames(mat4),row_ann_dt_w$taxon),]
col_b <- c("0" = "grey","1"= "blue","-1"= "red")
ha_left = HeatmapAnnotation(Cohort_1 = row_ann_dt_w$`Cohort 1`,
                            Cohort_2 = row_ann_dt_w$`Cohort 2`,
                            Cohort_3 = row_ann_dt_w$`Cohort 3`,
                            Cohort_4 = row_ann_dt_w$`Cohort 4`,
                            Cohort_5 = row_ann_dt_w$`Cohort 5`,
                            col = list(Cohort_1 = col_b,
                                       Cohort_2 = col_b,
                                       Cohort_3 = col_b,
                                       Cohort_4 = col_b,
                                       Cohort_5 = col_b), which = "row")




ht  =  Heatmap(mat4,name = "Log10(Rel_Abun)",
               row_split = rsplit[r_order],
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
               show_row_dend = F,
               row_names_max_width = max_text_width(rownames(mat4),
                                                    gp = gpar(fontsize = 12)),
               column_names_max_height = max_text_width(colnames(mat4),
                                                        gp = gpar(fontsize = 12)))


pdf(paste0(results_folder,"/Heatmap_Sig_tax_Early_recov_Cohorts_4_5.pdf"),height = 5, width = 15)
draw(ht)
dev.off()




