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
library(tidyverse)


det.cols <- c("Detected" = "#4daf4a", "Not detected" = "#e41a1c", "Insufficient data" = "#984ea3", "Probable" = "#377eb8")
strain.cols <- c("VE303-01" = "#1b9e77", "VE303-02" = "#d95f02", "VE303-03" = "#7570b3", "VE303-04" = "#e7298a", "VE303-05" = "#66a61e", "VE303-06" = "#e6ab02", "VE303-07" = "#a6761d", "VE303-08" = "#666666")
coh.cols <- c("Vanco" = "#e41a1c", "Sentinel/Cohort 1" = "#377eb8", "Cohort 2" = "#4daf4a", "Cohort 3" = "#ff7f00", "Cohort 4" = "#a65628", "Cohort 5" = "#984ea3", "Cohort 6" = "#999999")
phy.cols <- c("Actinobacteria" = "#ff3017", "Bacteroidetes" = "#ffd73a", "Firmicutes" = "#73c347", "Fusobacteria" = "#6bc77e", "Proteobacteria" = "#3c8fcd", "Spirochaetes" = "#a65628", "Synergistetes" = "#f781bf", "Tenericutes" = "#999999", "Verrucomicrobia" = "#d04ea6", "Euryarchaeota" = "black")

########Create result directory#############
mainDir <- "../results"
subDir <-"RF_Col_Viz"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

###############Read input files########################
col_day_start <- c("Baseline"= "#7570b3","Vanco" = "#a6761d","Early recovery" = "#1b9e77",
                   "Late recovery" = "#d95f02","Early no vanco" = "#e7298a","Late no vanco" = "#66a61e")

# Import dataset
# Add appropriated dated folder 
input_data <-  readRDS("../results/Input_matrix/input_data_rf.RDS")
abun_data <- readRDS("../results/Input_matrix/ve303_abun_rf.RDS")[["Species"]]
abun_data_coh6 <- readRDS("../results/Input_matrix/ve303_abun_coh6.RDS")[["Species"]]

# Visualize Importance
lev_tax <- "Class"
for(lev_tax in c("Class","Genus","Species")){
#  print(lev_tax)
  # First visualize the importance along with frequency:
  # Remove negative importance and 0 importance and p-val < 0.05
  final_imp_dt <-  read.csv(paste0("../results/RF_colonization_model/","VarImp_",lev_tax,".csv"))
  ntrials <-  as.numeric(table(final_imp_dt$Predictors)[1])
  final_imp_dt <- final_imp_dt[final_imp_dt$VarImp > 0,]
  final_imp_dt <- final_imp_dt[final_imp_dt$p.value < 0.1,]
  
  # Summarize the variable importance across the runs:
  rf_imp_sum <- final_imp_dt %>%
    group_by(Predictors) %>%
    summarize(Mean_VarImp = mean(VarImp),
              Freq = n() )%>% arrange(desc(Mean_VarImp))%>%
    as.data.frame()
  
  # Filter Predictors that are atleast signficant in half of the total reps
  rf_imp_sum <-  rf_imp_sum[rf_imp_sum$Freq >= 0.4*ntrials,]
  rf_imp_sum <- rf_imp_sum %>% arrange(desc(Mean_VarImp))
  final_imp_dt <- final_imp_dt[final_imp_dt$Predictors %in% rf_imp_sum$Predictors,]
  final_imp_dt$Predictors <- factor(final_imp_dt$Predictors,levels = rev(as.character(rf_imp_sum$Predictors)))
  rf_imp_sum$Predictors <- factor(rf_imp_sum$Predictors,levels = rev(as.character(rf_imp_sum$Predictors)))
  
  
  library(ggprism)
  freq_p <- ggplot(rf_imp_sum, aes(x = Predictors, y = Freq)) + geom_bar(stat = "identity",fill = "#008080",color = "black") +
    coord_flip()+
    theme_classic()+
    scale_y_reverse()+
    theme_prism()+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x =  element_text(size = 25))
  
  # Mean variable importance
  vimp_p <- ggplot(rf_imp_sum, aes(y = Mean_VarImp, x = Predictors)) +
    geom_linerange(
      aes(x = Predictors, ymin = 0, ymax = Mean_VarImp), 
      color = "lightgray", size = 1.5
    )+
    geom_point(color = "blue", size = 3)+
    ggpubr::color_palette("jco")+
    # theme_pubclean()+
    #geom_boxplot(outlier.colour = NA, fill = "lightblue") +
    coord_flip()+
    theme_prism()+
    theme(axis.text.y = element_text(size = 25),
          axis.title.x =  element_text(size = 25),
          axis.title.y =  element_text(size = 25))  
  
  
  imp_bugs <- as.character(rf_imp_sum$Predictors)
  # Import ALE
  tot_slope_dt <- read.csv(paste0("../results/RF_colonization_model/","ALE_Dir_",lev_tax,".csv"))
  tot_slope_dt$Dir <- ifelse(tot_slope_dt$Slope > 0,1,-1)
  tot_slope_dt <-  tot_slope_dt[tot_slope_dt$pred %in% imp_bugs,]
  
  
  library(tidyverse)
  sum_slope_dt <-  tot_slope_dt %>%
    group_by(pred)%>%
    summarise(Freq = n(),
              Sum = sum(Dir))%>%
    data.frame()
  sum_slope_dt$Dir <- as.character(ifelse(sum_slope_dt$Sum > 0,1,-1))
  
  sum_slope_dt  <- sum_slope_dt[sum_slope_dt$pred %in% imp_bugs,]
  sum_slope_dt$pred <-  factor(sum_slope_dt$pred, levels = rev(imp_bugs))
  
  # library(ggplot2)
  # gg_slope <- ggplot(sum_slope_dt, aes(x = 1, y = pred )) +
  #   geom_tile(color = "white", size = 0.1, aes(fill = Dir))+
  #   # scale_fill_gradient2( high = "#e34a33", mid = "#2c7fb8")+
  #   #scale_fill_gradient( low = scales::muted("blue"), high = scales::muted("red"))+
  #   scale_fill_manual(name = "",values = c("-1" = scales::muted('blue'), "1" = scales::muted('red')))+
  #   theme_classic()+
  #   theme_prism()+
  #   theme(
  #     axis.title.x =element_text(size = 25),
  #     legend.position = "none",
  #     axis.ticks.x = element_blank(),
  #     axis.text.y = element_blank(),
  #     axis.ticks.y = element_blank(),
  #     axis.text.x=element_blank(),
  #     axis.title.y = element_blank(),
  #     axis.text.y.left= element_blank())+
  #   xlab("Relation")
  
  library(ggplot2)
  gg_slope <- ggplot(sum_slope_dt, aes(x = factor(1), y = pred )) +
    geom_raster(aes(fill = Dir),vjust = 1)+
    # scale_fill_gradient2( high = "#e34a33", mid = "#2c7fb8")+
    #scale_fill_gradient( low = scales::muted("blue"), high = scales::muted("red"))+
    scale_fill_manual(name = "",values = c("-1" = scales::muted('blue'), "1" = scales::muted('red')))+
    theme_classic()+
    theme_prism()+
    # scale_x_discrete(expand = c(0,0))+
    #  scale_y_discrete(expand = c(0,0)) +
    coord_fixed(ratio=1)+
    theme(axis.title.x =element_text(size = 20),
          legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x=element_blank(),
          axis.title.y = element_blank(),
          axis.text.y.left= element_blank())+
    xlab("Relation")

  #gg_slope
  
  library(cowplot)
  pdf(paste(results_folder,paste0('VIMP_Freq_',lev_tax,'.pdf'),sep="/"),height = 8, width = 20,useDingbats = F)
  print(plot_grid(freq_p,gg_slope, vimp_p, labels = c('', '',''),nrow = 1,rel_widths = c(0.5,0.1,1),label_size = 10))
  dev.off()
  
} 
  # XY Plots  for baseline
  # Subset from input data
  sp_abun_dt <-  input_data[["Species"]]
  sp_abun_dt_sel  <- sp_abun_dt[,colnames(sp_abun_dt) %in% imp_bugs]
  sp_abun_dt_sel <- sp_abun_dt_sel[,!grepl("logFC",names(sp_abun_dt_sel))]
  sp_abun_dt_sel <-  cbind(Subject.ID = sp_abun_dt$Subject.ID,sp_abun_dt_sel)
  #sp_abun_dt_sel <- sp_abun_dt_sel[,!grepl("logFC",names(sp_abun_dt_sel))]
  
  #Summarise VE303 abundance per subject
  mean_ve_abun <- abun_data %>%
    group_by(Subject.ID)%>%
    summarise(Mean_abun = mean(sum_rel_abun))
  mean_ve_abun <-  mean_ve_abun[order(mean_ve_abun$Mean_abun),]
  
  sp_abun_dt_sel <-  sp_abun_dt_sel %>%
    inner_join(mean_ve_abun)
  
  sel_melt_dt <- sp_abun_dt_sel %>%
    pivot_longer(cols = baseline_Alistipes.putredinis.CAG.67:post_vanco_Veillonella.sp..DORA_B_18_19_23,
                 names_to = "Species",
                 values_to = "RA")
  # Make XY plot for these bugs
  
  p_line <- ggplot(sel_melt_dt,aes(x = RA+1e-5, y = Mean_abun+ 1e-5)) +
    geom_point(size = 2, color = "black", shape = 21, fill = "blue", alpha = 0.3)+
    facet_wrap(~Species,scales = "free")+
    # geom_line()+
    theme_clean()+
    scale_y_log10()+
    scale_x_log10()+
    ylab("TOtal VE303 abundance")
  
  print(p_line)
  
  pdf(paste0(results_folder,paste0("/XY_Plot_Species_abundances_Coh1_5.pdf")),height = 6, width = 10)
  #pdf(paste0(results_folder,paste0("/Sig_Firmicutes_Coh1_5_wo_time.pdf")),height = 6, width = 10)
  print(p_line)
  dev.off()  
  
  
  # Make a heatmap of species with logFC corresponding to  VE303-total abundance:
  input_mat <-  input_data[["Species"]]
  tot_X_dt <- input_mat
  tot_X_dt$Subject.ID <-  NULL
  abun_data <- abun_data[abun_data$Subject.ID %in% input_mat$Subject.ID,]
  
  mat_ht <- tot_X_dt[,c(as.character(rf_imp_sum$Predictors))]
  
  # Separate into two matrices and combine
  mat_baseline <-  mat_ht[,grep("baseline",colnames(mat_ht))]
  mat_baseline <-  log10(mat_baseline + 1e-5)
  mat_logFC <-  mat_ht[,grep("logFC",colnames(mat_ht))]
  
  #Summarise VE303 abundance per subject
  mean_ve_abun <- abun_data %>%
    group_by(Subject.ID)%>%
    summarise(Mean_abun = mean(sum_rel_abun))
  mean_ve_abun <-  mean_ve_abun[order(mean_ve_abun$Mean_abun),]
  
  
  tot_mat_ht <- cbind(Subject.ID = input_mat$Subject.ID,mat_baseline,mat_logFC,VE303 = log10(mean_ve_abun$Mean_abun+1e-5))
  mat_ht <-  t(mat_logFC)
  mat_ht <- mat_ht[,match(mean_ve_abun$Subject.ID,tot_mat_ht$Subject.ID)]
  
  
  # Row annotation:
  row_ann_dt <-  sum_slope_dt
  row_ann_dt <-  row_ann_dt[match(rownames(mat_ht),row_ann_dt$pred  ),]
  
  
  library(gtools)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(data.table)
  library(hues)
  library(scales)
  #splitcols <- c(rep("Baseline",ncol(mat_baseline)), rep("LogFC",ncol(mat_logFC)),"VE303")
  col_b <- c("0" = "grey","-1"= muted("blue"),"1"= muted("red"))
  ha_left = HeatmapAnnotation(Relation = row_ann_dt$Dir,
                              col = list(Relation = col_b), which = "row")
  
  library(scico)
  library(wacolors)
  col_ve303 <- colorRamp2(rev(c(seq(-0.5,-6,-0.5))),viridis::viridis(12))
  
  ha_column = HeatmapAnnotation(Total_VE303 = log10(mean_ve_abun$Mean_abun + 1e-5),
                                col=list(Total_VE303 = col_ve303),annotation_name_side = "left")
  col_mat <- colorRamp2(seq(-15,15,by = 5),colorRampPalette(c("#7f0000", "white","#00007f"))(7))
  library(ComplexHeatmap)
  ht  =  Heatmap(mat_ht,name = "LogFC",
                 # row_split = rsplit[r_order],
                 gap = unit(2, "mm"),
                 #column_split = splitcols,
                 col = col_mat,
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
                 row_names_max_width = max_text_width(rownames(mat_ht),
                                                      gp = gpar(fontsize = 14)),
                 column_names_max_height = max_text_width(colnames(mat_ht),
                                                          gp = gpar(fontsize = 12)))
  
  
  pdf(paste0(results_folder,paste0("/Heatmap_Sig_",lev_tax,"_Coh1_5.pdf")),height =4, width = 10)
  draw(ht)
  dev.off()
  
  
  # Now Display the prediction error:
  library(caret)
  error_func <-  function(pred_dt){
    
    
    r2 <-  R2(pred_dt$Pred, pred_dt$True)
    error_dt <- data.frame(
      R2 = 1 - sum(((pred_dt$True-pred_dt$Pred)^2)/sum((pred_dt$True-mean(pred_dt$True))^2)),
      #R2adj = 1- (1-r2)*(21-1)/(21-) 
      RMSE = RMSE(pred_dt$Pred, pred_dt$True),
      MAE = MAE(pred_dt$Pred, pred_dt$True)
    )
    return(error_dt)
  }
  
  # Visualize Importance
  lev_tax <- "Species"
  pred_list <-  list()
  err_list <- list()
  for(lev_tax in c("Class","Order","Genus","Species","null_data")[1:5]){
    
    pred_dt <-  read.csv(paste0("../results/RF_colonization_model/","Prediction_error_",lev_tax,".csv"))
    pred_list[[lev_tax]] <- pred_dt
    
    i = 1
    err_sub_list <-  list()
    for(i in unique(pred_dt$trial)){
      pred_sub  <- pred_dt %>%
        filter(trial == i)
      error  <- error_func(pred_sub)
      err_sub_list[[i]] <- error
      
    }
    error_dt <- do.call("rbind",err_sub_list) 
    error_dt$Taxa <-  lev_tax
    err_list[[lev_tax]] <-  error_dt
  }  
  
  
  err_final_dt <-  do.call("rbind",err_list)
  err_final_dt_m  <- melt(err_final_dt)
  
  library(ggthemes)
  library(cowplot)
  library(wesanderson)
  library(wacolors)
  library(ggprism)
  err_p <- ggplot(err_final_dt_m, aes(x = variable,y = value ,fill = Taxa))+
    geom_boxplot(alpha = 0.8)+
    facet_wrap(~variable,scales = "free")+
    scale_fill_wa_d(wacolors$rainier, reverse=TRUE)+
    #scale_fill_manual(values = wes_palette("Rushmore1"))+
    theme_prism()+
    xlab("Error Type")+
    ylab("Value")
  
  pdf(paste0(results_folder,paste0("/Error_boxplot.pdf")),height =5, width = 10)
  print(err_p)
  dev.off()
  
  