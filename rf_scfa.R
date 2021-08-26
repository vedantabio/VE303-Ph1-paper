
# Libraries 

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
library(data.table)
library(here)
library(latex2exp)
library(randomForest)
library(parallel)

rm(list=ls())




theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.6, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

# Color schemes 

det.cols <- c("Detected" = "#4daf4a", "Not detected" = "#e41a1c", "Insufficient data" = "#984ea3", "Probable" = "#377eb8")
strain.cols <- c("VE303-01" = "#1b9e77", "VE303-02" = "#d95f02", "VE303-03" = "#7570b3", "VE303-04" = "#e7298a", "VE303-05" = "#66a61e", "VE303-06" = "#e6ab02", "VE303-07" = "#a6761d", "VE303-08" = "#666666")
coh.cols <- c("Vanco" = "#e41a1c", "Sentinel/Cohort 1" = "#377eb8", "Cohort 2" = "#4daf4a", "Cohort 3" = "#ff7f00", "Cohort 4" = "#a65628", "Cohort 5" = "#984ea3", "Cohort 6" = "#999999")
phy.cols <- c("Actinobacteria" = "#ff3017", "Bacteroidetes" = "#ffd73a", "Firmicutes" = "#73c347", "Fusobacteria" = "#6bc77e", "Proteobacteria" = "#3c8fcd", "Spirochaetes" = "#a65628", "Synergistetes" = "#f781bf", "Tenericutes" = "#999999", "Verrucomicrobia" = "#d04ea6", "Euryarchaeota" = "black")
cluster.cols <- c("I" = "#fb9a99", "II" = "#cab2d6", "IV" = "#a6cee3", "IX" = "#e31a1c", "XIVa" = "#1f78b4", "XV" = "#fdbf6f", "XVI" = "#ff7f00", "XVII" = "#b2df8a", "XVIII" = "#33a02c")

########Create result directory#############
mainDir <- "../results"
subDir <- paste0("metabolomics/",gsub("-","_",Sys.Date()),'_SCFA_RF_results')
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

# RF data directory
main_data_Dir <- "../RF_data"
sub_data_Dir <- paste0("data/",gsub("-","_",Sys.Date()))
dir.create(file.path(main_data_Dir, sub_data_Dir), showWarnings = TRUE,recursive = TRUE)
data_folder <- paste(main_data_Dir,sub_data_Dir,sep="/")

###############Read input files########################

library("phyloseq")
phy_mic <- readRDS("../data/phy_mic.rds")


#### SCFA dynamics ####

phy_met <-  readRDS("../data/phy_met.rds")
phy_sel <- prune_samples(sample_names(phy_met), phy_mic)


# Get the OTU table 
met_dt <-  data.frame(t(otu_table(phy_met)))
met_dt$Sample <-  rownames(met_dt)

sam_dt <- data.frame(sample_data(phy_sel))

merg_sam_dt <- merge(sam_dt, met_dt, by.x = "Ph_id",by.y =  "Sample")
rownames(merg_sam_dt) <- merg_sam_dt$Ph_id

sample_data(phy_sel ) <-  sample_data(merg_sam_dt)

# Remove Cohort 6 from the analysis
phy_sel <-  subset_samples(phy_sel, cohort_id_long !="Cohort 6")
phy_sel <- prune_taxa(taxa_sums(phy_sel)>0, phy_sel)


met_names <-  names(met_dt)[1:(length(names(met_dt)) -2)]




# # Melting the phyloseq object (Longer format from wider format)
# phy_m_scfa<-psmelt(phy_met)
# # SCFA and BA names
# SCFA_array<-names(phy_m_scfa)[13:36]
# 
# 
# 
# phy_m_scfa<-phy_m_scfa[,c("oc_time","oc_subject","oc_cohort",SCFA_array)]
# phy_m_scfa<-unique(phy_m_scfa)
# nrow(phy_m_scfa)
# 


phy_post<-phy_sel
phy_post<-subset_samples(phy_post,Collection.Day.Norm > 6)
phy_post<-subset_samples(phy_post,Collection.Day.Norm <= 60)

phy_post <- prune_taxa(taxa_sums(phy_post) > 0, phy_post)

phy_post <- prune_samples(sample_sums(phy_post) > 0, phy_post)





# prevalence/abundance analysis first
prevdf = apply(X = otu_table(phy_post),
               MARGIN = ifelse(taxa_are_rows(phy_post), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phy_post))
prevdf$OTU <-  rownames(prevdf)

prevdf <- prevdf[order(prevdf$TotalAbundance,decreasing = T), ]
prevdf$OTU <- factor(as.character(prevdf$OTU) , levels = as.character(prevdf$OTU))
prevdf$Prop <- prevdf$TotalAbundance/sum(prevdf$TotalAbundance)

prev_frac<-0.05
prev_cutoff <- prev_frac*nsamples(phy_post) # Cut-off
ab_cutoff <- 1e-6 # Cut-off
# Prevalence
prevdf_fil <- prevdf[prevdf$Prevalence >= prev_cutoff, ]
# Abundance
prevdf_fil <- prevdf_fil[prevdf_fil$Prop >= ab_cutoff, ]

library("ggplot2")

p1 <-  ggplot(prevdf,aes(x = Prevalence, y = Prop)) +  geom_jitter(size = 2 ,alpha=0.5, shape=21, fill='purple')+
  xlab("Prevalence")+
  ylab("Proportion")+
  ggtitle(paste0("Prevalence vs. Abundance"))+
  scale_y_log10()+
  theme_Publication()+
  geom_vline(xintercept=prev_cutoff, linetype="dashed",
             color = "red", size=1)+
  geom_hline(yintercept=ab_cutoff, linetype="dashed",
             color = "red", size=1)+
  theme(axis.text.x=element_text(angle = -90, hjust = 0),panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey50"),
        panel.ontop = TRUE)
print(p1)

pdf(paste(results_folder,paste0('prev_ab','_xy_plot_v1.pdf'),sep="/"),height = 4, width = 5)
print(p1)
dev.off()

# Filter prevalent and abundant taxa
phy_mer_f = prune_taxa(rownames(prevdf_fil), phy_post)

grep("VE303",taxa_names(phy_mer_f),value = T)

taxa_names(phy_mer_f)


sample_sums(phy_mer_f)



#### Run random forest by choosing random samples equally from each subject ####
ntrials = 50
nRep_rf = 20
ncores <- 20
nTopFeatures <- 30

library(randomForest)
#library("plyr")
library(reshape2)



library(doParallel)
cl <- makeCluster(ncores)
registerDoParallel(cl)


is<-"Ursodeoxycholic.Acid"

# Removing Vancomycin
SCFA_array <-  met_names

# We look at the time window between 6 and 20 days (Effect of VE303)
phy_rf<-phy_mer_f
# phy_rf<-subset_samples(phy_rf,oc_time>6)
# phy_rf<-subset_samples(phy_rf,oc_time<=20)



tp_pred_union <- list()
#slope_dt_list <- list()
#cor_dt_list <- list()
for (is in SCFA_array){
  
  # Create an empty list to save the randomForest object
  rf_list<-list()
  itrial <- 1
  # Remove samples where SCFA/BA are NA
  good_samples <- sample_names(phy_rf)[!is.na(get_variable(phy_rf,is))]
  # Print SCFA/BA and total number of samples
  print(is)
  print(length(good_samples))
  if(length(good_samples)> 20){
    phy_day_sel <- prune_samples(good_samples,phy_rf)
    # Remove taxa that have zero contribution 
    phy_day_sel <- prune_taxa(taxa_sums(phy_day_sel)>0,phy_day_sel)
    phy_day_sel_sam_sub <- sample_data(phy_day_sel)[,c("Subject.ID","cohort_id_long","Ph_id")]
    
    # list of training data for random forest in each trial
    sam_train_data <-  list()
    # Loop over ntrials for RF
    for (itrial in seq(1,ntrials)){
      
      set.seed(itrial + 1000)
      #Pick one random sample(irrespective of time) from each subjects
      sam_sub <- phy_day_sel_sam_sub %>% group_by(Subject.ID) %>% sample_n(1)
      # Subset samples that are sampled
      phy_sel_s <- subset_samples(phy_day_sel,Ph_id %in% sam_sub$Ph_id)
      # Remove taxa that have zero abundance
      phy_sel_s <-  prune_taxa(taxa_sums(phy_sel_s)> 0, phy_sel_s)
      # Variable /SCFA to be predicted
      Y <-  sample_data(phy_sel_s)[,is] 
      # Predictors  
      X <- data.frame(otu_table(phy_sel_s))
      tdata <- t(as.matrix(X))
      # Total X and Y combined
      tdata_comb <-  cbind(Y=Y,tdata)
      train <- tdata_comb
      names(train)[1]<-"Y"
      rf_model <- foreach(nrep = seq(1,nRep_rf), .packages=c('randomForest')) %dopar%{
        set.seed(nrep + 2000)
        rf <- randomForest(Y ~.,train,ntree = 1001, importance = T,keep.inbag = T)
      }  
      rf_list[[itrial]]<-rf_model
      
      # Save the training set for each random forest in a list for Partial dependence
      sam_train_data[[itrial]] <- train
    }
    
    # Save all RFR objects  
    saveRDS(rf_list,paste(data_folder,paste0('RF_list_',is,'.RDS'),sep="/"))
    
    # Save all training data for Partial Dependence plots
    saveRDS(sam_train_data,paste(data_folder,paste0('train_RF_',is,'.RDS'),sep="/"))
    
    
    # Compute variable importance dataframe across all ntrials and nreps and 
    # Error
    # il <-1
    # jl <- 1
    
    p_val_pimp <-  0.1
    library(vita)
    rf_error_list <- list()
    rf_sig_list <-  list()
    for (il in seq(1,length(rf_list))){
      
      rf_list_rep<-list()
      rf_oob_error_rep <-list()
      
      for (jl in seq(1,nRep_rf)){
        t_dt <- sam_train_data[[il]]
        pimp <- PIMP(t_dt[,-1], t_dt$Y,rf_list[[il]][[jl]],parallel=TRUE, ncores=20 ,seed = 340)
        p_test <- PimpTest(pimp)
        #pimp_all <- data.frame( summary(p_test,pless = 0.05))
        pimp_all <- data.frame(orig =  p_test$VarImp, p_test$pvalue)
        
        #pimp_all <- pimp_all[pimp_all$VarImp > 0 & pimp_all$p.value <= 0.1,]
        #pimp_all <-  pimp$PerVarImp
        #pimp_av <- data.frame(av = rowMeans(pimp_all), orig =  pimp$VarImp)
        
        #pimp.t.reg = PimpTest(pimp)
        #dt_p <- data.frame(summary(pimp.t.reg,pless = 0.1))
        #pimp <-compVarImp(t_dt[,-1], t_dt$Y,rf_list[[il]][[jl]] )
        #pimp_dt <- data.frame(imp = pimp$importance, imp_sd =  pimp$importanceSD)  
        #pimp_dt$orig_imp <-  importance(rf_list[[il]][[jl]],type=1, scale=FALSE)
        
        
        imp_dt <- pimp_all     
        #names(imp_dt)
        #imp_dt <- as.data.frame(importance(rf_list[[il]][[jl]]))
        imp_dt$Predictors <- rownames(imp_dt)
        rownames(imp_dt) <- NULL
        rf_list_rep[[jl]]<-imp_dt
        
        # Error
        err_dt <- as.data.frame(rf_list[[il]][[jl]]$mse)
        err_dt$n_tree <- 1:nrow(err_dt)
        err_dt$n_rep <- jl
        rf_oob_error_rep[[jl]]<-err_dt
      }
      
      rf_list_dt<-do.call('rbind',rf_list_rep)
      # rf_list_dt_summary<- rf_list_dt %>% group_by(Predictors) %>% summarise_all(funs(mean))
      # rf_list_dt_summary <-  rf_list_dt_summary[order(rf_list_dt_summary$X.IncMSE, decreasing = T),]
      # rf_sig<- rf_list_dt_summary
      #[1:nTopFeatures,]
      rf_list_dt$ntrial <- il
      #rf_sig$ntrial <- il
      rf_sig_list[[il]] <-  rf_list_dt
      
      # Error
      rf_err_dt<-do.call('rbind',rf_oob_error_rep)
      names(rf_err_dt)[1] <- "MSE"
      rf_err_dt$ntrial <-  il
      rf_error_list[[il]] <-  rf_err_dt
    }
    #  stopCluster(cl)
  #  library(dplyr)
    # Variable Importance plot for ntrials
    rf_list_dt_sel <- do.call("rbind",rf_sig_list) 
    
    rf_list_dt_freq <-  rf_list_dt_sel %>% 
      filter(p.value < p_val_pimp & VarImp > 0)%>%
      add_count(Predictors) %>%
      select(Predictors,n)%>% unique
    
    
    
    
    #rf_list_dt_sel <- rf_list_dt_sel[rf_list_dt_sel$p.value < 0.1 & rf_list_dt_sel$VarImp > 0,]
    
    rf_mean_dt_sel_summary <-  rf_list_dt_sel %>%
      group_by(Predictors )%>%
      summarise(av_vimp=mean(VarImp),sd_mse =sd(VarImp)/sqrt(length(VarImp))) 
    
    rf_mean_dt_sel_summary <- rf_mean_dt_sel_summary[rf_mean_dt_sel_summary$Predictors %in% rf_list_dt_freq$Predictors,]
    
    #rf_mean_dt_all <- rf_mean_dt_sel_summary[order(rf_mean_dt_sel_summary$pimp,decreasing = T),]
    
    rf_mean_dt_all <- rf_mean_dt_sel_summary[match(rf_list_dt_freq$Predictors,rf_mean_dt_sel_summary$Predictors),]
    rf_mean_dt_all$freq <-  rf_list_dt_freq$n
    
    rf_mean_dt_all <-  rf_mean_dt_all[order(rf_mean_dt_all$av_vimp,decreasing = T),]
    
    
    # Save the dataframe for variable importance for all predictors
    write.csv(rf_mean_dt_all,paste(results_folder,paste0(Sys.Date(),'-RF_regression-',is,'-VIMP.csv'),sep="/"))
    if(nrow(rf_mean_dt_all) >= nTopFeatures){
      topPred<-rf_mean_dt_all[1:nTopFeatures,]
    }else {
      topPred <- rf_mean_dt_all
    }
    topPred$Predictors <- factor(topPred$Predictors , levels = rev(topPred$Predictors))
    topPred$freq <-  topPred$freq/(nRep_rf*ntrials)
    
    #topPred <-  topPred[topPred$freq >= 0.1,]                               
    vimp_p <- ggplot() +
      geom_point(data=topPred, aes(x = Predictors, y = av_vimp,fill=freq),shape=21, size = 3)+
      #geom_errorbar(data=topPred, aes(x = Predictors, ymin = av_vimp-sd_mse,ymax=av_vimp+sd_mse,color=freq))+
      scale_fill_gradient(low = "white", high = "darkred")+
      #scale_color_gradient(low = "white", high = "darkred")+
      coord_flip()+
      # theme_base()+
      theme_Publication()+
      xlab("")+
      ylab("Increase in MSE if permuted")+
      ggtitle(gsub("\\.|X"," ",is))+
      theme(legend.position="top",legend.text = element_text(size = 10))
    #vimp_p
    
    pdf(paste(results_folder,paste0(Sys.Date(),'-RF_regression-',is,'-VIMP.pdf'),sep="/"),height = 8, width = 8)
    print(vimp_p)
    dev.off()
    
    # Save the dataframe for variable importance for top predictors
    topPred$Analyte <- is
    write.csv(topPred,paste(data_folder,paste0('RFR_VIMP_',is,'.csv'),sep="/"))

  }
}


# Compute the union of the top predictors
tp_pred <- Reduce(union, tp_pred_union)
unique(unlist(tp_pred_union))

# Heatmap using the coefficient of the linear model Y(Strain Abundance) ~ Bacterial Abundance   
slope_dt <-  bind_rows(slope_dt_list, .id = 'strain')
slope_dt$Intercept <-  NULL
slope_dt_w <- dcast(slope_dt, variable ~ strain)
names(slope_dt_w)
rownames(slope_dt_w) <- slope_dt_w$variable
slope_dt_w$variable <-  NULL 


SBA <-  c("Deoxycholic.Acid","Glycoursodeoxycholic.Acid","Lithocholic.Acid")
PBA <-  setdiff(grep("Acid",names(slope_dt_w), value =T),SBA)
SCFA <- setdiff(names(slope_dt_w),c(SBA,PBA))
met_analyte <-  data.frame(Type = c(rep("Primary BA", length(PBA)),rep("Secon BA", length(SBA)),
                                    rep("SCFA",length(SCFA))),
                           Analyte =  c(PBA,SBA,SCFA))           
# Bile Acids (BA) 
slope_dt_w_BA <- slope_dt_w[,names(slope_dt_w) %in% c(SBA,PBA)]

slope_dt_w_BA <- slope_dt_w_BA[rowSums(is.na(slope_dt_w_BA)) != ncol(slope_dt_w_BA), ]

met_analyte_ba <- met_analyte[met_analyte$Type != "SCFA",]
met_analyte_ba <- met_analyte_ba[match(names(slope_dt_w_BA), met_analyte_ba$Analyte),]

library("ComplexHeatmap")

annotations <- data.frame(Class =  met_analyte_ba$Type)
ha_column = HeatmapAnnotation(annotations,col=list(Class=c("Primary BA" = "#1b9e77",
                                                           "Secon BA"= "#d95f02")))
library("circlize")
mat<-as.matrix(slope_dt_w_BA)     
mat2<-scale(t(mat), scale = F, center = F)
mat2<-t(mat2)
# Setting NA to zero
mat2[is.na(mat2)] <- 0


jet.colors <-   colorRamp2(c(-0.3, 0, 0.3), c("blue", "white", "red"))
library("circlize")
ht1 = Heatmap(mat2, name = "Slope", column_title = NA, top_annotation = ha_column,
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "complete",row_names_side = "left", km=1, color_space = "LAB",
              row_dend_side="right",
              col=jet.colors,
              clustering_method_columns = "ward.D",
              width=2, row_names_max_width = unit(8, "cm"),show_column_names= T,
              row_names_gp = gpar(fontsize = 9), cluster_columns = T,cluster_rows = T,na_col="white")
#ht1

pdf(paste(results_folder,'/',Sys.Date(),'-Heatmap_slope_str_BA_pred.pdf',sep=""),height = 25, width = 10)
draw(ht1)
dev.off()


# SCFA
slope_dt_w_SCFA <- slope_dt_w[,names(slope_dt_w) %in% SCFA]
slope_dt_w_SCFA <- slope_dt_w_SCFA[rowSums(is.na(slope_dt_w_SCFA)) != ncol(slope_dt_w_SCFA), ]
met_analyte_scfa <- met_analyte[met_analyte$Type != "SCFA",]
met_analyte_scfa <- met_analyte_scfa[match(names(slope_dt_w_SCFA), met_analyte_scfa$Analyte),]


mat<-as.matrix(slope_dt_w_SCFA)     
mat2<-scale(t(mat), scale = F, center = F)
mat2<-t(mat2)
# Setting NA to zero
mat2[is.na(mat2)] <- 0
library(circlize)
jet.colors <-   colorRamp2(c(-.8, 0, .8), c("blue", "white", "red"))
library("circlize")
ht1 = Heatmap(mat2, name = "Correlation", column_title = NA, 
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "complete",row_names_side = "left", km=1, color_space = "LAB",
              row_dend_side="right",
              col=jet.colors,
              clustering_method_columns = "ward.D",
              width=2, row_names_max_width = unit(8, "cm"),show_column_names= T,
              row_names_gp = gpar(fontsize = 9), cluster_columns = T,cluster_rows = T,na_col="white")
#ht1

pdf(paste(results_folder,'/',Sys.Date(),'_Heatmap_Slope_SCFA_str_scfa_pred.pdf',sep=""),height = 10, width = 10)
draw(ht1)
dev.off()


# Heatmap using the correlation coefficient using rmcorr   

is <- SCFA_array[1]
cor_res_list <-  list()
for (is in names(cor_dt_list)){
  corr_dt <-  cor_dt_list[[is]]  
  corr_dt$logY <- log10(corr_dt[,c(is)])
  corr_dt$oc_cohort <- factor(corr_dt$oc_subject)
  other_spp <- names(corr_dt)[34:(ncol(corr_dt) -1)]
  
  cor_res_dt <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(cor_res_dt) <- c("Species","Corr","p_value","conf_int1","conf_int2")
  
  fit.values <- list()
  library("rmcorr")
  other <- other_spp[1]
  for(other in other_spp){
    
    my.rmc <- rmcorr(participant = oc_cohort, measure1 = "logY", measure2 = other, dataset = corr_dt)
    cor_res_dt <- rbind(cor_res_dt, data.frame(Species = other,Corr = my.rmc$r, p_value = my.rmc$p, conf_int1 = my.rmc$CI[1], conf_int2 = my.rmc$CI[1]))
    fit.values[[other]] <- my.rmc$model$fitted.values
  }
  
  #fit_dt <- data.frame(do.call(cbind, fit.values))
  #fit_dt  <-  cbind(sample_dt,data.frame(do.call(cbind, fit.values)))
  
  #cor_res_dt$Species <- other_spp
  cor_res_dt <- cor_res_dt[ order(cor_res_dt$p_value),]
  cor_res_dt$p_adjust <- p.adjust(cor_res_dt$p_value, method = "BH")
  cor_res_list[[is]] <- cor_res_dt
}
#write.csv(cor_res_dt,paste(results_folder,paste0('Corr_result_MDRO_CFU',treatment,'.csv'),sep="/"))

# Save the correlation data in a list
rmcorr_dt <-  bind_rows(cor_res_list, .id = 'strain')
# P-value cutoff <= 0.05
rmcorr_dt <- rmcorr_dt[rmcorr_dt$p_value <=  0.05,]
cor_dt_w <- dcast(rmcorr_dt, Species ~ strain,value.var = "Corr")
names(cor_dt_w)
rownames(cor_dt_w) <- cor_dt_w$Species
cor_dt_w$Species <-  NULL 


SBA <-  c("Deoxycholic.Acid","Glycoursodeoxycholic.Acid","Lithocholic.Acid")
PBA <-  setdiff(grep("Acid",names(cor_dt_w), value =T),SBA)
SCFA <- setdiff(names(cor_dt_w),c(SBA,PBA))
met_analyte <-  data.frame(Type = c(rep("Primary BA", length(PBA)),rep("Secon BA", length(SBA)),
                                    rep("SCFA",length(SCFA))),
                           Analyte =  c(PBA,SBA,SCFA))           

# Bile Acids (BA) 
cor_dt_w_BA <- cor_dt_w[,names(cor_dt_w) %in% c(SBA,PBA)]

cor_dt_w_BA <- cor_dt_w_BA[rowSums(is.na(cor_dt_w_BA)) != ncol(cor_dt_w_BA), ]

met_analyte_ba <- met_analyte[met_analyte$Type != "SCFA",]
met_analyte_ba <- met_analyte_ba[match(names(cor_dt_w_BA), met_analyte_ba$Analyte),]

library("ComplexHeatmap")

annotations <- data.frame(Class =  met_analyte_ba$Type)
ha_column = HeatmapAnnotation(annotations,col=list(Class=c("Primary BA" = "#1b9e77",
                                                           "Secon BA"= "#d95f02")))
library("circlize")

mat<-as.matrix(cor_dt_w_BA)     
mat2<-scale(t(mat), scale = F, center = F)
mat2<-t(mat2)
# Setting NA to zero
mat2[is.na(mat2)] <- 0
jet.colors <-   colorRamp2(c(-.8, 0, .8), c("blue", "white", "red"))

ht1 = Heatmap(mat2, name = "Correlation", column_title = NA, top_annotation = ha_column,
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "complete",row_names_side = "left", km=1, color_space = "LAB",
              row_dend_side="right",
              col=jet.colors,
              clustering_method_columns = "ward.D",
              width=2, row_names_max_width = unit(8, "cm"),show_column_names= T,
              row_names_gp = gpar(fontsize = 9), cluster_columns = T,cluster_rows = T,na_col="white")
#ht1

pdf(paste(results_folder,'/',Sys.Date(),'_Heatmap_Correlation_BA_str_scfa_pred.pdf',sep=""),height = 12, width = 10)
draw(ht1)
dev.off()

# SCFA
cor_dt_w_SCFA <- cor_dt_w[,names(cor_dt_w) %in% SCFA]
cor_dt_w_SCFA <- cor_dt_w_SCFA[rowSums(is.na(cor_dt_w_SCFA)) != ncol(cor_dt_w_SCFA), ]
met_analyte_scfa <- met_analyte[met_analyte$Type != "SCFA",]
met_analyte_scfa <- met_analyte_scfa[match(names(cor_dt_w_SCFA), met_analyte_scfa$Analyte),]


mat<-as.matrix(cor_dt_w_SCFA)     
mat2<-scale(t(mat), scale = F, center = F)
mat2<-t(mat2)
# Setting NA to zero
mat2[is.na(mat2)] <- 0
library(circlize)
jet.colors <-   colorRamp2(c(-.8, 0, .8), c("blue", "white", "red"))
library("circlize")
ht1 = Heatmap(mat2, name = "Correlation", column_title = NA, 
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "complete",row_names_side = "left", km=1, color_space = "LAB",
              row_dend_side="right",
              col=jet.colors,
              clustering_method_columns = "ward.D",
              width=2, row_names_max_width = unit(8, "cm"),show_column_names= T,
              row_names_gp = gpar(fontsize = 9), cluster_columns = T,cluster_rows = T,na_col="white")
#ht1

pdf(paste(results_folder,'/',Sys.Date(),'_Heatmap_Correlation_SCFA_str_scfa_pred.pdf',sep=""),height = 10, width = 10)
draw(ht1)
dev.off()

