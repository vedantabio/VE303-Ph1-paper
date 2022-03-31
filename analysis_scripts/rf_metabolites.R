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
phy_mic <- readRDS("../Data/phy_mic_abs.rds")
#### SCFA dynamics ####
phy_met <-  readRDS("../Data/phy_met.rds")
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

met_names <-  names(met_dt)[1:(length(names(met_dt)) -1)]

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
SCFA_array <-  met_names[!met_names %in% c("Vancomycin")]
phy_rf<-phy_mer_f
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
    phy_day_sel_sam_sub <- sample_data(phy_day_sel)[,c("Subject.ID","cohort_id_long","Ph_id")] %>% data.frame()
    
    # list of training data for random forest in each trial
    sam_train_data <-  list()
    # Loop over ntrials for RF
    for (itrial in seq(1,ntrials)){
      
      set.seed(itrial + 1000)
      #Pick one random sample(irrespective of time) from each subjects
      sam_sub <- phy_day_sel_sam_sub %>% group_by(Subject.ID) %>% sample_n(1)%>% data.frame()
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
        pimp_all <- data.frame(orig =  p_test$VarImp, p_test$pvalue)
        imp_dt <- pimp_all     
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
      rf_list_dt$ntrial <- il
      rf_sig_list[[il]] <-  rf_list_dt
      
      # Error
      rf_err_dt<-do.call('rbind',rf_oob_error_rep)
      names(rf_err_dt)[1] <- "MSE"
      rf_err_dt$ntrial <-  il
      rf_error_list[[il]] <-  rf_err_dt
    }
    # Variable Importance plot for ntrials
    rf_list_dt_sel <- do.call("rbind",rf_sig_list) 
    
    rf_list_dt_freq <-  rf_list_dt_sel %>% 
      filter(p.value < p_val_pimp & VarImp > 0)%>%
      add_count(Predictors) %>%
      select(Predictors,n)%>% unique
    
    rf_mean_dt_sel_summary <-  rf_list_dt_sel %>%
      group_by(Predictors )%>%
      summarise(av_vimp=mean(VarImp),sd_mse =sd(VarImp)/sqrt(length(VarImp))) 
    
    rf_mean_dt_sel_summary <- rf_mean_dt_sel_summary[rf_mean_dt_sel_summary$Predictors %in% rf_list_dt_freq$Predictors,]
    
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



