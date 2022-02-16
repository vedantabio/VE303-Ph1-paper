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
library(tidyverse)
library(randomForest)
library(parallel)
library(ggplot2)
library(ggthemes)
library(ALEPlot)
library(dplyr)


# Function to compute ale plots
ale_func <-  function(rep){
  if(imp_bugs[len] %in% names(train_list[[rep]])){
    idx <- which(names(train_list[[rep]]) == imp_bugs[len]) -1
    ale_dt  <- data.frame(ALEPlot(X = train_list[[rep]][,-1],X.model =  rf_list_lin[[rep]], 
                                  pred.fun = yhat, J=idx, K = 40, NA.plot =F))
    
    ale_dt$pred <-  imp_bugs[len]
    ale_dt$nrep <- rep
    res <-  ale_dt
    # rm(ale_dt)
  }
}

########Create result directory#############
mainDir <- "../results"
subDir <- paste0("metabolomics/",gsub("-","_",Sys.Date()),'_SCFA_ALE_results')
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")


# RF data directory
main_data_Dir <- "../RF_data"
# Make a sub directory of ALE under RF_data
sub_data_Dir <- paste0("data/ALE/",gsub("-","_",Sys.Date()))
dir.create(file.path(main_data_Dir, sub_data_Dir), showWarnings = TRUE,recursive = TRUE)
data_folder <- paste(main_data_Dir,sub_data_Dir,sep="/")


# Add appropriate dated folder from RF analysis
vimp_dt <- list.files(path = "../RF_data/data/2019_08_04/",pattern = ".csv",full.names = T) %>%
  map_dfr(read_csv)

scfa <-  unique(vimp_dt$Analyte)
# Make the list of random forest linear to parallelize
# Loop over analyte
analyte <- scfa[1]
for(analyte in scfa){
  
  
  print("Reading RF list ")
  rf_list <- readRDS(paste0("../RF_data/data/2019_08_04/RF_list_",analyte,".RDS"))
  nrep <-  length(rf_list[[1]])
  rf_list_lin <-  do.call("c",rf_list)
  # Remove rf_list 
  rm(rf_list)
  # Read training data 
  print("Reading training data for each RF")
  #train_list <- readRDS(paste0("../RF_data/data/2019_05_01/train_RF_",analyte,".RDS"))
  train_list <- readRDS(paste0("../RF_data/data/2019_08_04/train_RF_",analyte,".RDS"))
  train_list <-  rep(train_list, each = nrep)
  # VIMP
  vimp_scfa <-  vimp_dt[vimp_dt$Analyte == analyte,]
  vimp_scfa <- vimp_scfa[vimp_scfa$freq > 0.2,]
  # Loop over bugs 
  imp_bugs <- as.character(vimp_scfa$Predictors)
  print("Computing ALE for important predictors")
  yhat <- function(X.model, newdata) as.numeric(randomForest:::predict.randomForest(X.model, newdata))
  pp_bug_list <- list()
  for(len in seq_along(imp_bugs)){
    ale_list <- mclapply(1:length(rf_list_lin),ale_func,mc.cores = 20)
    pp_bug_list[[len]] <- do.call("rbind",ale_list)
  }
  
  # Clean up snapshot files
  # It fills up the home directory quick
  # find . -name "*.snapshot" -delete
  #system('find /home/sbhattarai/.rstudio/sessions/active/session-1f8c8f24/graphics-r3/ -name "*.snapshot" -delete', intern = T)
  system('find /home/sbhattarai/.rstudio/sessions/active/*/graphics-r3/ -name "*.snapshot" -delete', intern = T)
  rm(rf_list_lin)
  
  # ALE  plot for ntrials* nrep
  pp_dt_final <- do.call("rbind",pp_bug_list)
  pp_dt_final$Analyte <-  analyte
  # Save the ALE data
  write.csv(pp_dt_final,paste(data_folder,paste0('RF_ALE_',analyte,'.csv'),sep="/"))
  
  pp_dt_final$rep<-  pp_dt_final$nrep %% nrep
  pp_dt_final$rep[pp_dt_final$rep == 0 ] <-  nrep
  pp_dt_final$trial<- pp_dt_final$nrep %/% nrep + 1
  pp_dt_final$trial[pp_dt_final$rep == nrep ] <-  pp_dt_final$trial[pp_dt_final$rep == nrep] + 1
  
  pp_dt_final$id <- paste0(pp_dt_final$rep,pp_dt_final$trial)
  pp_dt_final$pred <- factor(pp_dt_final$pred,levels = rev(imp_bugs))
  p_p <- ggplot(data = pp_dt_final ,aes(x = x.values, y = f.values)) +
    geom_line(aes(group = nrep), alpha = 0.2)+
    #geom_line(data = pp_dt_final_median ,aes(x = x, y = y_med, color = trial))+
    #geom_smooth()+
    #geom_line(data = pp_dt_final_median ,aes(x = x, y = y_med, color = "red"))+
    theme_base()+
    facet_wrap(~pred,scales = "free")+
    xlab("Microbes")+
    ylab(analyte)+
    ggtitle(analyte)+
    theme(legend.position="none")
  
  pdf(paste(results_folder,paste0(Sys.Date(),'_RF_ALE_',analyte,'.pdf'),sep="/"),height = 20, width = 25)
  print(p_p)
  dev.off()
  rm(pp_dt_final)
  rm(pp_bug_list)
}






