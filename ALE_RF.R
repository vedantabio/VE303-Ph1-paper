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
#Function to compute partial dependence from random forest

# pp_func <- function(rep){
#   if(imp_bugs[len] %in% names(train_list[[rep]])){
#     pp_dt <- data.frame(partialPlot(x = rf_list_lin[[rep]],
#                                     pred.data = train_list[[rep]],
#                                     imp_bugs[len],
#                                     plot = F))
#     pp_dt$pred <-  imp_bugs[len]
#     pp_dt$nrep <- rep
#     res <-  pp_dt
#   }
# }


rep <-  1
len <- 1
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
sub_data_Dir <- paste0("data/ALE/",gsub("-","_",Sys.Date()))
dir.create(file.path(main_data_Dir, sub_data_Dir), showWarnings = TRUE,recursive = TRUE)
data_folder <- paste(main_data_Dir,sub_data_Dir,sep="/")



# Get variable importance 
# vimp_dt <- list.files(path = "../RF_data/data/2019_05_01/",pattern = ".csv",full.names = T) %>%
#   map_dfr(read_csv)

vimp_dt <- list.files(path = "../RF_data/data/2019_08_04/",pattern = ".csv",full.names = T) %>%
  map_dfr(read_csv)

scfa <-  unique(vimp_dt$Analyte)

# Make the list of random forest linear to parallelize
# Loop over analyte
analyte <- scfa[1]
for(analyte in scfa){
  
  
  print("Reading RF list ")
  #rf_list <- readRDS(paste0("../RF_data/data/2019_05_01/RF_list_",analyte,".RDS"))
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
  
  # Partial dependence (Half parallel)
  print("Computing ALE for important predictors")
  
  yhat <- function(X.model, newdata) as.numeric(randomForest:::predict.randomForest(X.model, newdata))
  #yhat(rf_list_lin[[1]],train_list[[1]][,-1])

  #pred <- randomForest:::predict.randomForest(rf_list_lin[[1]],train_list[[1]][,-1])

  #ALEPlot(X = train_list[[1]][,-1],X.model =  rf_list_lin[[1]], pred.fun = yhat, J=idx, K = 40, NA.plot = TRUE)
  # RMSE <- function(x,y){
  #   a <- sqrt(sum((x-y)^2)/length(y))
  #   return(a)
  # }
  # RMSE1 <- RMSE(as.numeric(pred), rf_list_lin[[1]]$y)/mean(rf_list_lin[[1]]$y)
  # 
  pp_bug_list <- list()
  for(len in seq_along(imp_bugs)){
     ale_list <- mclapply(1:length(rf_list_lin),ale_func,mc.cores = 20)
    pp_bug_list[[len]] <- do.call("rbind",ale_list)
  }
  
  # Clean up snapshot files
  # It fills up the home directory quick
  #library(lime)
 # find . -name "*.snapshot" -delete

  
  
  
  
  #system('find /home/sbhattarai/.rstudio/sessions/active/session-1f8c8f24/graphics-r3/ -name "*.snapshot" -delete', intern = T)
  
  # 
  # pp_bug_list <- list()
  # for(len in seq_along(imp_bugs)){
  #   pp_bug_list[[len]] <- do.call("rbind",mclapply(1:length(rf_list_lin),pp_func,mc.cores = 20))
  # }
  
  rm(rf_list_lin)
  
  # Partial Dependence  plot for ntrials* nrep
  pp_dt_final <- do.call("rbind",pp_bug_list)
  pp_dt_final$Analyte <-  analyte
  # Save the partial dependence data
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
#  print(p_p)
  pdf(paste(results_folder,paste0(Sys.Date(),'_RF_ALE_',analyte,'.pdf'),sep="/"),height = 20, width = 25)
  print(p_p)
  dev.off()
  #ggsave(paste(results_folder,paste0(Sys.Date(),'_RF_ALE_',analyte,'.pdf'),sep="/"),height = 20, width = 25)
  #rm(p_p)
  rm(pp_dt_final)
  rm(pp_bug_list)
  #   geom_line()+
  # #  theme_Publication()+
  #   facet_wrap(~pred,scales = "free")+
  #   xlab("Microbes")+
  #   ylab(analyte)+
  #   ggtitle(analyte)+
  #   theme(legend.position="none")
  # #vimp_p
  # 
  # pdf(paste(results_folder,paste0(Sys.Date(),'_RF_ALE_',analyte,'.pdf'),sep="/"),height = 20, width = 25)
  # print(p_p)
  # dev.off()
  
  # # Plot the log of the SCFA
  # p_p <- ggplot(data = pp_dt_final ,aes(x = log10(x), y = log10(y),color = id)) +
  #   geom_line()+
  #   theme_Publication()+
  #   facet_wrap(~pred,scales = "free")+
  #   xlab("Microbes")+
  #   ylab(is)+
  #   ggtitle(is)+
  #   theme(legend.position="none")
  # #vimp_p
  # 
  # pdf(paste(results_folder,paste0(Sys.Date(),'_Log_RF_partial_dependence_',is,'.pdf'),sep="/"),height = 20, width = 25)
  # print(p_p)
  # dev.off()
  # 
  
  # # Use the lm
  # # Slope calculation SCFA Vs Microbes to identify positive/ negative relation
  # pp_dt_final_slope <- pp_dt_final
  # pp_dt_final$logx <- log10(pp_dt_final$x + 10^(-8))
  # pp_dt_final$logy <- log10(pp_dt_final$y)
  # 
  # 
  # pp_dt_final <- pp_dt_final[complete.cases(pp_dt_final),]
  # #comb_dt <- comb_dt[comb_dt$variable =="base_bm",]
  # slope_dt_list <- pp_dt_final %>%
  #   group_by(pred,nrep) %>%
  #   do({
  #     mod <-  data.frame(summary(lm( logy ~ logx , data = .))$coefficients)[,c(1,4)]
  #     names(mod) <- c("val","p_val")
  #     mod$coeff <- rownames(mod)
  #     
  #     data.frame(Intercept = mod$val[1],
  #                Slope = mod$val[2],
  #                p_val_int =  mod$p_val[1],
  #                p_val_slope =  mod$p_val[2] )
  #   }) %>% as.data.frame()
  # 
  # slope_dt <- data.frame(slope_dt_list)
  # 
  # slope_dt$pred <-  factor(slope_dt$pred, levels = rev(imp_bugs))
  # slope_dt <- slope_dt[complete.cases(slope_dt),]
  # # Boxplot for slopes ordered by the VIMP
  # # Plot the log of the SCFA
  # #slope_dt <- slope_dt[slope_dt$p_val_slope <= 0.05,]
  # 
  # slope_dt$dir <-  ifelse(slope_dt$Slope >0 ,"pos","neg")
  # 
  # slope_box_p <- ggplot(data = slope_dt ,aes(y = Slope, x = pred, color = dir)) +
  #   geom_boxplot()+
  #   geom_jitter(alpha = 0.4)+
  #   theme_base()+
  #   scale_color_manual(name = "Slope",values = c("neg" = "red","pos" = "blue"))+
  #   xlab("Microbes")+
  #   ylab("Slope")+
  #   ggtitle(analyte)+
  #   theme(axis.text.x=element_text(angle=90, hjust=1))
  # #vimp_p
  # 
  # pdf(paste(results_folder,paste0(Sys.Date(),'_Box_plot_slope_ALE_',analyte,'.pdf'),sep="/"),height = 10, width = 20)
  # print(slope_box_p)
  # dev.off()
  # 
  
  
}








# Partial dependence (Half parallel)
pp_bug_list <- list()
for(len in seq_along(imp_bugs)){
    pp_bug_list[[len]] <- do.call("rbind",mclapply(1:length(rf_list_lin),pp_func,mc.cores = 20))
}







imp_bugs <- as.character(rf_list_dt_sel_summary$Predictors)


#Function to compute partial dependence from random forest

pp_func <- function(il, jl){
  if(imp_bugs[len] %in% names(train)){
    pp_dt <- data.frame(partialPlot(x = rf_list[[il]][[jl]],
                                    pred.data = sam_train_data[[il]],
                                    imp_bugs[len],
                                    plot = F))
    pp_dt$pred <-  imp_bugs[len]
    pp_dt$nrep <- jl
    pp_dt$ntrial <-  il
    res <-  pp_dt
  }
}

# Partial dependence (Half parallel)
pp_bug_list <- list()
for(len in seq_along(imp_bugs)){
  pp_rep_list <- list()
  for(jl in 1:nRep_rf){
    pp_rep_list[[jl]]  <-  do.call("rbind",mclapply(1:length(rf_list),pp_func, jl = jl))
  }
  pp_bug_list[[len]] <-   do.call("rbind",pp_rep_list)
}


# Partial Dependence  plot for ntrials* nrep
pp_dt_final <- do.call("rbind",pp_bug_list)

# Save the partial dependence data
write.csv(pp_dt_final,paste(results_folder,paste0(Sys.Date(),'_RF_partial_dependence_',is,'.csv'),sep="/"))



pp_dt_final$id <- paste0(pp_dt_final$nrep,pp_dt_final$ntrial)
pp_dt_final$pred <- factor(pp_dt_final$pred,levels = rev(levels(rf_list_dt_sel_summary$Predictors)))
p_p <- ggplot(data = pp_dt_final ,aes(x = x, y = y,color = id)) +
  geom_line()+
  theme_Publication()+
  facet_wrap(~pred,scales = "free")+
  xlab("Microbes")+
  ylab(is)+
  ggtitle(is)+
  theme(legend.position="none")
#vimp_p

pdf(paste(results_folder,paste0(Sys.Date(),'_RF_partial_dependence_',is,'.pdf'),sep="/"),height = 20, width = 25)
print(p_p)
dev.off()

# Plot the log of the SCFA
p_p <- ggplot(data = pp_dt_final ,aes(x = log10(x), y = log10(y),color = id)) +
  geom_line()+
  theme_Publication()+
  facet_wrap(~pred,scales = "free")+
  xlab("Microbes")+
  ylab(is)+
  ggtitle(is)+
  theme(legend.position="none")
#vimp_p

pdf(paste(results_folder,paste0(Sys.Date(),'_Log_RF_partial_dependence_',is,'.pdf'),sep="/"),height = 20, width = 25)
print(p_p)
dev.off()


# Use the lm
# Slope calculation SCFA Vs Microbes to identify positive/ negative relation
pp_dt_final <- do.call("rbind",pp_bug_list)
pp_dt_final$id <- paste0(pp_dt_final$nrep,pp_dt_final$ntrial)
pp_dt_final$logx <- log10(pp_dt_final$x + 10^(-8))
pp_dt_final$logy <- log10(pp_dt_final$y)


pp_dt_final <- pp_dt_final[complete.cases(pp_dt_final),]
#comb_dt <- comb_dt[comb_dt$variable =="base_bm",]
slope_dt_list[[is]] <- pp_dt_final %>%
  group_by(pred,id) %>%
  do({
    mod <-  data.frame(summary(lm( logy ~ logx , data = .))$coefficients)[,c(1,4)]
    names(mod) <- c("val","p_val")
    mod$coeff <- rownames(mod)

    data.frame(Intercept = mod$val[1],
               Slope = mod$val[2],
               p_val_int =  mod$p_val[1],
               p_val_slope =  mod$p_val[2] )
  }) %>% as.data.frame()

slope_dt <- data.frame(slope_dt_list[[is]])

slope_dt$pred <-  factor(slope_dt$pred, levels = rev(imp_bugs))
slope_dt <- slope_dt[complete.cases(slope_dt),]
# Boxplot for slopes ordered by the VIMP
# Plot the log of the SCFA
#slope_dt <- slope_dt[slope_dt$p_val_slope <= 0.05,]

slope_dt$dir <-  ifelse(slope_dt$Slope >0 ,"pos","neg")

slope_box_p <- ggplot(data = slope_dt ,aes(y = Slope, x = pred, color = dir)) +
  geom_boxplot()+
  geom_jitter(alpha = 0.5)+
  theme_base()+
  scale_color_manual(name = "Slope",values = c("neg" = "red","pos" = "blue"))+
  xlab("Microbes")+
  ylab("Slope")+
  ggtitle(is)+
  theme(axis.text.x=element_text(angle=90, hjust=1))
#vimp_p

pdf(paste(results_folder,paste0(Sys.Date(),'_Box_plot_slope_partial_dependence_',is,'.pdf'),sep="/"),height = 10, width = 20)
print(slope_box_p)
dev.off()
