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
subDir <- "RF_colonization_model"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

df <- file.info(list.files("../results/Input_matrix", full.names = T))
rownames(df)[which.max(df$mtime)]


# Import dataset
# Add the dated folder for the input data
input_data <-  readRDS("../results/Input_matrix/input_data_rf.RDS")
abun_data <- readRDS("../results/Input_matrix/ve303_abun_rf.RDS")[["Species"]]

names(input_data)

# Repeats on RF 
ntrials = 100
# Number of cores to be used
ncores<-15
# Top 30 important features
nTopFeatures <- 30

# Run feature selection before running RF:
sig_boruta_list <- list()
lev_tax <- "Species"
for(lev_tax in names(input_data)[1:4]){
  input_mat <-  input_data[[lev_tax]]
  tot_X_dt <- input_mat
  tot_X_dt$Subject.ID <-  NULL
  abun_data <- abun_data[abun_data$Subject.ID %in% input_mat$Subject.ID,]
  
  tot_X_dt <-  tot_X_dt[,c(grep("logFC",names(tot_X_dt)) ,
                           grep("baseline",names(tot_X_dt)),
                           grep("post",names(tot_X_dt)))]
  library("randomForest")
  
  itrial <- 1
  dt_list <- list() 
  boruta_sp <-  c()
  for (itrial in seq(1,ntrials)){
    set.seed(itrial + 1000)
    # Pick one random sample(irrespective of time) from each subjects
    abun_data_sampled <- abun_data %>% group_by(Subject.ID) %>% sample_n(1)
    abun_data_sampled <- abun_data_sampled[match(input_mat$Subject.ID,abun_data_sampled$Subject.ID),]
    tdata_comb <-  cbind(tot_X_dt ,Y = abun_data_sampled$sum_rel_abun)
    # Boruta to select important species.
    # First step to filter out species that relates to response
    library(Boruta)
    set.seed(3000)
    print("Selecting important species using Boruta")
    boruta_model <- Boruta(Y~., data=tdata_comb, doTrace=2,mcAdj = T, maxRuns = 200)
    #Significant predictors from Boruta are also the Top 10 predictors from RF
    boruta_signif <- names(boruta_model$finalDecision[boruta_model$finalDecision %in% c("Confirmed","Tentative")])  # collect Confirmed and Tentative variables
    boruta_sp <- c(boruta_sp,boruta_signif)
  }
  
  sig_boruta_list[[lev_tax]] <- unique(boruta_sp)
  
}  


saveRDS(sig_boruta_list,paste0(results_folder,"/sig_boruta_list.RDS"))

# Make a random matrix in R to test the null dataset

row_size = nrow(input_data[["Species"]])
col_size = length(sig_boruta_list[["Species"]])
set.seed(100)
mat_null = matrix(rnorm(row_size*col_size),nrow=row_size)
mat_null <- data.frame(Subject.ID = input_data[["Species"]]$Subject.ID,mat_null)

input_data[["null_data"]] <- mat_null
# Add null data in input_data
# Perform leave on out cross validation
lev_tax <- "Species"
for(lev_tax in c(names(input_data)[1:4],"null_data")){
  
  input_mat <-  input_data[[lev_tax]]
  tot_X_dt <- input_mat
  tot_X_dt$Subject.ID <-  NULL
  abun_data <- abun_data[abun_data$Subject.ID %in% input_mat$Subject.ID,]
  
  sig_boruta <-  NULL
  if(lev_tax != "null_data"){
    sig_boruta <- sig_boruta_list[[lev_tax]]
  }else{
    tot_X_dt <- tot_X_dt
    sig_boruta <-  colnames(tot_X_dt)
  }
  
  # Test on LOOCV
  library("randomForest")
  library("rsample")
  pred_list<-list()
  itrial <- 1
  
  for (itrial in seq(1,ntrials)){
    
    print(paste0("Running LOOCV in: ",itrial ))
    set.seed(itrial + 1000)
    # Pick one random sample(irrespective of time) from each subjects
    abun_data_sampled <- abun_data %>% group_by(Subject.ID) %>% sample_n(1)
    abun_data_sampled <- abun_data_sampled[match(input_mat$Subject.ID,abun_data_sampled$Subject.ID),]
    tdata_comb <-  cbind(tot_X_dt[,c(sig_boruta)] ,Y = abun_data_sampled$sum_rel_abun) 
    
    # Use LOOCV and gather the prediction on test cases
    set.seed(199) 
    folds<- vfold_cv(tdata_comb, v=nrow(tdata_comb),repeats=1) 
    len<-1
    
    pred_vec <- c()
    true_vec <- c()
    for(len in 1:nrow(folds)){
      split<-folds$splits[[len]]
      # Training set
      train_dt <- analysis(split)
      # Only Predictors
      train_X <-  train_dt[,!colnames(train_dt) %in% "Y"]
      # Test set
      test_dt <- assessment(split)
      # Only Predictors
      test_X <-  test_dt[,!colnames(test_dt) %in% "Y"]
      
      set.seed(itrial + 1000)
      rf <- randomForest(Y ~.,train_dt,ntree = 3000, importance = T)
      RFestimated =  predict(rf, newdata= test_X)
      
      # Prediction
      pred_vec <- c(pred_vec,as.numeric(RFestimated))
      true_vec <- c(true_vec,test_dt$Y)
      
    }
    
    pred_dt <-  data.frame(Pred = pred_vec,True = true_vec)
    pred_dt$trial <- itrial 
    pred_dt$Taxa_level <-  lev_tax
    pred_list[[itrial]] <- pred_dt
    
  }
  
  pred_final_dt <-  do.call("rbind",pred_list)
  write.csv(pred_final_dt,paste(results_folder,paste0('Prediction_error_',lev_tax,'.csv'),sep="/"))
  # pdf(paste(results_folder,paste0(Sys.Date(),'-VIMP.pdf'),sep="/"),height = 6, width = 8)
  # print( vimp_p)
  # dev.off()
  
  
  library(caret)
  error_dt <- data.frame(
    R2 = R2(pred_final_dt$Pred, pred_final_dt$True),
    RMSE = RMSE(pred_final_dt$Pred, pred_final_dt$True),
    MAE = MAE(pred_final_dt$Pred, pred_final_dt$True)
  )
  
  # ggplot(data = pred_final_dt, aes(x= True, y = Pred, group = trial))+
  #    geom_point(size = 3,shape = 21, color = "black")+
  #   geom_smooth(se = F,method = "lm")+
  #   geom_abline(slope = 1,intercept = 0, color = "red")+
  #   theme_classic()+
  #   coord_fixed()
  
  
  # Run randomforest using all Boruta selected predictors
  library("randomForest")
  rf_list<-list()
  itrial <- 1
  dt_list <- list() 
  for (itrial in seq(1,ntrials)){
    
    set.seed(itrial + 1000)
    # Pick one random sample(irrespective of time) from each subjects
    abun_data_sampled <- abun_data %>% group_by(Subject.ID) %>% sample_n(1)
    abun_data_sampled <- abun_data_sampled[match(input_mat$Subject.ID,abun_data_sampled$Subject.ID),]
    tdata_comb <-  cbind(tot_X_dt[,c(sig_boruta)] ,Y = abun_data_sampled$sum_rel_abun) 
    
    train <- tdata_comb
    dt_list[[itrial]] <- train 
    rf <- randomForest(Y ~.,train,ntree = 3000, importance = T)
    rf_list[[itrial]]<-rf
    
    
  }
  
  library(vita)
  # Permutated importance:
  
  # A combined list of random forest object.
  rf_imp_list<- list()
  
  itrial <-  1
  # Now run PIMP on each model and save it 
  for (itrial in seq(1,ntrials)){
    model <- rf_list[[itrial]]
    train_dt <- dt_list[[itrial]]
    train_X <- train_dt[,1:(ncol(train_dt)-1)] 
    print("Running PIMP")
    pimp <- PIMP(train_X, train_dt$Y,model, parallel=TRUE, ncores = 20, seed = 1)
    pimp_test <- PimpTest(pimp)
    pimp_all <- data.frame(orig =  pimp_test$VarImp, pimp_test$pvalue)
    imp_df<- pimp_all     
    imp_df$Predictors <- rownames(imp_df)
    rownames(imp_df) <- NULL
    #imp_df$nrep <-  rep
    rf_imp_list[[itrial]]<-imp_df
    # rf_comb_model_list <- list.append(rf_comb_model_list,rf_list[[il]][[jl]])
    # rf_comb_t_data_list <-list.append(rf_comb_t_data_list,dt_list[[il]]) 
  }
  
  # Process the importance dataframe for each repetition:
  final_imp_dt  <-   do.call("rbind",rf_imp_list)
  write.csv(final_imp_dt,paste(results_folder,paste0('VarImp_',lev_tax,'.csv'),sep="/"))
  
  # Summarize the variable importance across the runs:
  rf_imp_sum <- final_imp_dt %>%
    group_by(Predictors) %>%
    summarize(Mean_VarImp = mean(VarImp),
              Freq = n() )%>% arrange(desc(Mean_VarImp))%>%
    as.data.frame()
  
  # Filter Predictors that are atleast signficant in half of the total reps
  rf_imp_sum <-  rf_imp_sum[rf_imp_sum$Freq >= 0.3*ntrials,]
  rf_imp_sum <- rf_imp_sum %>% arrange(desc(Mean_VarImp))
  final_imp_dt <- final_imp_dt[final_imp_dt$Predictors %in% rf_imp_sum$Predictors,]
  final_imp_dt$Predictors <- factor(final_imp_dt$Predictors,levels = rev(as.character(rf_imp_sum$Predictors)))
  rf_imp_sum$Predictors <- factor(rf_imp_sum$Predictors,levels = rev(as.character(rf_imp_sum$Predictors)))
  
  
  # ALE plots for the sign:
  library(ALEPlot)
  library(parallel)
  yhat <- function(X.model, newdata) as.numeric(randomForest:::predict.randomForest(X.model, newdata))
  ALE_summary_list <- list()
  itrial <-  1
  # Now run PIMP on each model and save it 
  for (itrial in seq(1,ntrials)){
    model <- rf_list[[itrial]]
    train_dt <- dt_list[[itrial]]
    train_X <- train_dt[,1:(ncol(train_dt)-1)] 
    
    imp_bugs <-  as.character(rf_imp_sum$Predictors)
    
    train_data <-  train_dt[,!colnames(train_dt) %in% "Y"]
    
    ale_func <-  function(rep){
      
      if(imp_bugs[len] %in% colnames(train_data)){
        idx <- which(colnames(train_data) == imp_bugs[len]) 
        ale_dt  <- data.frame(ALEPlot(X = train_data,X.model =  rf_list[[itrial]], 
                                      pred.fun = yhat, J=idx, K = 40, NA.plot =F))
        
        ale_dt$pred <-  imp_bugs[len]
        ale_dt$nrep <- rep
        res <-  ale_dt
        # rm(ale_dt)
      }
    }
    
    pp_bug_list <- list()
    for(len in seq_along(imp_bugs)){
      print(len)
      ale_dt <- ale_func(itrial)
      pp_bug_list[[len]] <- ale_dt
    }
    
    pp_dt_final <- do.call("rbind",pp_bug_list)
    
    pp_dt_final$f.values <-  as.numeric(pp_dt_final$f.values)
    pp_dt_final$x.values <-  as.numeric(pp_dt_final$x.values)
    pp_dt_final <-  na.omit(pp_dt_final)
    
    
    pp_dt_final$pred <- factor(pp_dt_final$pred,levels = imp_bugs)
    
    library(tidyr)
    library(dplyr)
    slope_dt <- pp_dt_final %>%
      group_by(pred) %>%
      do({
        mod <-  data.frame(summary(lm( f.values ~ x.values , data = .))$coefficients)[,c(1,4)]
        names(mod) <- c("val","p_val")
        mod$coeff <- rownames(mod)
        
        data.frame(Intercept = mod$val[1],
                   Slope = mod$val[2],
                   p_val_int =  mod$p_val[1],
                   p_val_slope =  mod$p_val[2] )
      }) %>% as.data.frame()
    
    head(slope_dt)
    
    
    slope_dt <- slope_dt[,c("pred","Slope")]
    slope_dt$ntrial <- itrial
    
    ALE_summary_list[[itrial]] <-  slope_dt
    
    
  }
  
  
  tot_slope_dt <- do.call("rbind",ALE_summary_list)
  write.csv(tot_slope_dt,paste(results_folder,paste0('ALE_Dir_',lev_tax,'.csv'),sep="/"))
  
  
}



