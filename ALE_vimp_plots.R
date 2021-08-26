# Make the VIMP 
library(tidyverse)
library(randomForest)
library(parallel)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)
library(phyloseq)

#  scale y in the range of -1 to 1
norm_range <- function(x , a , b){
  xnorm = (b- a)*((x - min(x))/(max(x)- min(x))) + a  
}

assignCols <- function(names,selection) {
  mymap<-list()
  availableCols <-c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                    "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
                    "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
                    "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
                    "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58")
  mymap[[1]]<-rev(availableCols[match(intersect(names,selection),names)])
  mymap[[2]]<-intersect(names,selection)
  mymap
};

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
            legend.key.size= unit(0.8, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


########Create result directory#############
mainDir <- "../results"
subDir <- paste0("metabolomics/",gsub("-","_",Sys.Date()),'_ALE_Vimp_plots')
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

# Read the VIMP files
# Get variable importance 
# vimp_dt <- list.files(path = "../RF_data/data/2019_05_01/",pattern = ".csv",full.names = T) %>%
#   map_dfr(read_csv)
# 
# ale_dt <-  list.files(path = "../RF_data/data/ALE/2019_05_03/",pattern = ".csv",full.names = T) %>%
#   map_dfr(read_csv)
# 
# vimp_dt <- list.files(path = "../RF_data/data/2019_05_29/",pattern = ".csv",full.names = T) %>%
#   map_dfr(read_csv)
# 
# ale_dt <-  list.files(path = "../RF_data/data/ALE/2019_05_30/",pattern = ".csv",full.names = T) %>%
#   map_dfr(read_csv)

vimp_dt <- list.files(path = "../RF_data/data/2019_08_04/",pattern = ".csv",full.names = T) %>%
  map_dfr(read_csv)

ale_dt <-  list.files(path = "../RF_data/data/ALE/2019_08_08/",pattern = ".csv",full.names = T) %>%
  map_dfr(read_csv)


all_analytes <-  unique(vimp_dt$Analyte)


phy_codex <- readRDS("../data/phy_mic.rds")


phy_sel <- prune_taxa(unique(vimp_dt$Predictors), phy_codex)

tax_dt <-  tax_table(phy_sel)
tax_dt <- data.frame(tax_dt@.Data)

tax_dt[] <- lapply(tax_dt, as.character)

tax_dt$Species <- rownames(tax_dt)


gtdb_tax <-  tax_dt[, grep("GTDB",names(tax_dt), value = T)]

incomplete_gtdb <-  gtdb_tax[!complete.cases(gtdb_tax),]
complete_gtdb <-  gtdb_tax[complete.cases(gtdb_tax),]

ncbi_tax <-  tax_dt[tax_dt$Species %in% rownames(incomplete_gtdb),][,1:7]
names(ncbi_tax) <-  names(gtdb_tax)

tot_gtdb_dt <- rbind(complete_gtdb, ncbi_tax)
tot_gtdb_dt$Species <-  rownames(tot_gtdb_dt)

match_tax_dt <-  tot_gtdb_dt
match_tax_dt$GTDB_phylum <- gsub("p__","",match_tax_dt$GTDB_phylum)
match_tax_dt$GTDB_phylum <- gsub("_.*","",match_tax_dt$GTDB_phylum)


match_tax_dt$GTDB_order <-  gsub("o__","",match_tax_dt$GTDB_order)

match_tax_dt[match_tax_dt$Species == "Odoribacter_sp._UNK.MGS_12","GTDB_phylum"] <-  "Bacteroidetes"
match_tax_dt[match_tax_dt$Species == "Odoribacter_sp._UNK.MGS_12","GTDB_order"] <-  "Bacteroidales"


match_tax_dt <- match_tax_dt[order(match_tax_dt$GTDB_phylum,match_tax_dt$GTDB_order),]

# Phylum and Order colors
unique(match_tax_dt$GTDB_phylum)

phylum_col <-  c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
names(phylum_col) <- unique(match_tax_dt$GTDB_phylum)
dt_col_phylum <- data.frame(phylum_col)
dt_col_phylum$Phylum <-  rownames(dt_col_phylum)
shades_num <- match_tax_dt %>%
  group_by(GTDB_phylum) %>%
  mutate(unique_types = n_distinct(GTDB_order))%>% select(GTDB_phylum,unique_types) %>% unique %>% as.data.frame()

shades_num <- merge(shades_num, dt_col_phylum, by.y = "Phylum", by.x = "GTDB_phylum")
shades_num$phylum_col <-  as.character(shades_num$phylum_col )
tax_sel <-  unique(match_tax_dt[,c("GTDB_phylum","GTDB_order")])

library(yingtools2)
#shades_num$phy_col <- as.vector(phylum_col[order(names(phylum_col))])
order_col <- mapply(function(x, y){ shades(color = x, ncolor = y, variation = 1)}, x= shades_num$phylum_col ,y = shades_num$unique_types,
                    SIMPLIFY = F)

order_col <- as.vector(unlist(order_col))
names(order_col) <- tax_sel$GTDB_order


# Make the list of random forest linear to parallelize
# Loop over analyte
analyte <- all_analytes[1]
slope_all_list <- list()
for(analyte in all_analytes){
  
  # VIMP 
  vimp_sel <-  vimp_dt %>% filter(Analyte == analyte & freq >= 0.2) %>% as.data.frame()
   # Select top 10
  vimp_sel <-  vimp_sel %>% top_n(10, freq) 
  
  
  
  vimp_sel$Predictors <- factor(vimp_sel$Predictors,levels = rev(vimp_sel$Predictors))
  
  # vimp_p <- ggplot() +
  #   geom_point(data=vimp_sel, aes(x = Predictors, y = av_vimp,fill=freq),shape=21, size = 3)+
  #   #geom_errorbar(data=topPred, aes(x = Predictors, ymin = av_vimp-sd_mse,ymax=av_vimp+sd_mse,color=freq))+
  #   scale_fill_gradientn(name = "Norm Freq    ",colours = c("white","darkred"),limits=c(0,1))+
  #   #scale_color_gradient(low = "white", high = "darkred")+
  #   coord_flip()+
  #   # theme_base()+
  #   theme_Publication()+
  #   xlab("")+
  #   ylab("Increase in MSE if permuted")+
  #   ggtitle(gsub("\\.|X"," ",analyte))+
  #   theme(legend.position="top",legend.text = element_text(size = 10))
  # #vimp_p
  # 
  # pdf(paste(results_folder,paste0(Sys.Date(),'-RF_regression-',analyte,'-VIMP.pdf'),sep="/"),height = 8, width = 8)
  # print(vimp_p)
  # dev.off()
  
  # ALE Plots 
  ale_sel_dt <- ale_dt %>% filter(Analyte == analyte) 
  ale_sel_dt <-  ale_sel_dt %>% filter(pred %in% levels(vimp_sel$Predictors)) 
  ale_sel_dt$pred <- factor(ale_sel_dt$pred,levels = rev(levels(vimp_sel$Predictors)))
 
  #ale_sel_dt <-  ale_sel_dt[ale_sel_dt$pred == "VE303_03",]
  
  # p_p <- ggplot(data = ale_sel_dt ,aes(x = x.values, y = f.values)) +
  #   geom_line(aes(group = nrep), alpha = 0.2)+
  #   #geom_line(data = pp_dt_final_median ,aes(x = x, y = y_med, color = trial))+
  #   #geom_smooth()+
  #   #geom_line(data = pp_dt_final_median ,aes(x = x, y = y_med, color = "red"))+
  #   theme_base()+
  #   facet_wrap(~pred,scales = "free")+
  #   xlab("Microbes")+
  #   ylab(analyte)+
  #   ggtitle(analyte)+
  #   theme(legend.position="none")
  # #  print(p_p)
  # pdf(paste(results_folder,paste0(Sys.Date(),'_RF_ALE_',analyte,'.pdf'),sep="/"),height = 20, width = 25)
  # print(p_p)
  # dev.off()
  
  
  # Fit a logistic growth curve and measure the  growth factor
  ale_sel_dt$logx <- log10(ale_sel_dt$x.values + 10^(-8))
  ale_sel_dt$norm_y <- norm_range(ale_sel_dt$f.values, 0, 1)
  
  # p_p <- ggplot(data = ale_sel_dt ,aes(x = logx, y = norm_y)) +
  #   geom_point()+
  #   #geom_line(aes(group = nrep), alpha = 0.2)+
  #   #geom_line(data = pp_dt_final_median ,aes(x = x, y = y_med, color = trial))+
  #   geom_smooth(method = "lm")+
  #   #geom_line(data = pp_dt_final_median ,aes(x = x, y = y_med, color = "red"))+
  #   theme_base()+
  #   facet_wrap(~pred,scales = "free")+
  #   xlab("Microbes")+
  #   ylab(analyte)+
  #   ggtitle(analyte)+
  #   theme(legend.position="none")
  # #p_p
  # pdf(paste(results_folder,paste0(Sys.Date(),'_RF_ALE_fit_',analyte,'.pdf'),sep="/"),height = 20, width = 25)
  # print(p_p)
  # dev.off()
  # 
  
  slope_dt_list <- ale_sel_dt %>%
    group_by(pred) %>%
    do({
      mod <-  data.frame(summary(lm( norm_y ~ logx , data = .))$coefficients)[,c(1,4)]
      names(mod) <- c("val","p_val")
      mod$coeff <- rownames(mod)
      
      data.frame(Intercept = mod$val[1],
                 Slope = mod$val[2],
                 p_val_int =  mod$p_val[1],
                 p_val_slope =  mod$p_val[2] )
    }) %>% as.data.frame()
  
  slope_dt <- data.frame(slope_dt_list)
  slope_dt$analyte <-  analyte
  slope_dt <- slope_dt[complete.cases(slope_dt),]
  
  slope_all_list[[analyte]] <- slope_dt
  
}  


# Create a heatmap with slope
ht_slope_dt <- do.call('rbind', slope_all_list)

ht_slope_sub <-  ht_slope_dt[,c("pred","analyte","Slope")]
library("reshape2")
slope_dt_w <- dcast(ht_slope_sub, pred ~ analyte)
names(slope_dt_w)
rownames(slope_dt_w) <- slope_dt_w$pred
slope_dt_w$pred <-  NULL 


# Bile Acids (BA) 

SBA <-  c("Deoxycholic.Acid","Lithocholic.Acid","Ursodeoxycholic.Acid",
          "Glycoursodeoxycholic.Acid","Tauroursodeoxycholic.Acid","Taurolithocholic.Acid",
          "Glycolithocholic.Acid","Taurodeoxycholic.Acid","Glycodeoxycholic.Acid")
PBA <-  setdiff(grep("Acid",names(slope_dt_w), value =T),SBA)
SCFA <- setdiff(names(slope_dt_w),c(SBA,PBA))
met_analyte <-  data.frame(Type = c(rep("Primary BA", length(PBA)),rep("Secon BA", length(SBA)),
                                    rep("SCFA",length(SCFA))),
                           Analyte =  c(PBA,SBA,SCFA))            

cor_dt_w_BA <- slope_dt_w[,grep("Acid",names(slope_dt_w)) ]

cor_dt_w_BA <- cor_dt_w_BA[rowSums(is.na(cor_dt_w_BA)) != ncol(cor_dt_w_BA), ]



met_analyte_ba <- met_analyte[met_analyte$Type != "SCFA",]
met_analyte_ba <- met_analyte_ba[match(names(cor_dt_w_BA), met_analyte_ba$Analyte),]

met_analyte_ba$Mechanism <- "bai Dehydroxylation"
met_analyte_ba$Mechanism[met_analyte_ba$Analyte %in% c("Glycoursodeoxycholic.Acid","Ursodeoxycholic.Acid","Tauroursodeoxycholic.Acid")]<- "Epimerization"
met_analyte_ba$Mechanism[met_analyte_ba$Type == "Primary BA"] <- NA

# Categorize the BA
met_analyte_ba$Cat <- "Unconjugated BA"
met_analyte_ba$Cat[met_analyte_ba$Analyte %in% 
                     c("Glycocholic.Acid","Taurocholic.Acid",
                       "Glycochenodeoxycholic.Acid","Taurochenodeoxycholic.Acid")]<- "Primary Conjugated BA"
met_analyte_ba$Cat[met_analyte_ba$Analyte %in% 
                     c("Glycodeoxycholic.Acid","Taurodeoxycholic.Acid",
                       "Glycolithocholic.Acid","Taurolithocholic.Acid",
                       "Glycoursodeoxycholic.Acid","Tauroursodeoxycholic.Acid")]<- "Secondary Conjugated BA"

met_analyte_ba$Cat[met_analyte_ba$Analyte %in% 
                     c("Deoxycholic.Acid",
                       "Lithocholic.Acid","Ursodeoxycholic.Acid")]<- "Secondary DeConjugated BA"



library("ComplexHeatmap")

annotations <- data.frame(Cat = met_analyte_ba$Cat,Class =  met_analyte_ba$Type, Mechanism = met_analyte_ba$Mechanism)



ha_column = HeatmapAnnotation(Cat = met_analyte_ba$Cat,
                              Class =  met_analyte_ba$Type,
                              Mechanism = met_analyte_ba$Mechanism,
                              col=list(Cat = c("Unconjugated BA" = "#b2182b","Primary Conjugated BA" = "#ef8a62",
                                                           "Secondary Conjugated BA"= "#67a9cf","Secondary DeConjugated BA" = "#2166ac" ),
                                                   Class=c("Primary BA" = "#1b9e77","Secon BA"= "#d95f02"),
                                                   Mechanism = c("bai Dehydroxylation" = "#f26391","Epimerization"= "#6e3999")
))
library("circlize")

mat<-as.matrix(cor_dt_w_BA)     
mat2<-scale(t(mat), scale = F, center = F)
mat2<-t(mat2)
# Setting NA to zero
mat2[is.na(mat2)] <- 0
jet.colors <-   colorRamp2(c(-0.01, 0, 0.05), c("blue", "white", "red"))

tax_dt <- match_tax_dt

tax_ba <- tax_dt[tax_dt$Species %in% rownames(mat2),]
tax_ba[] <- lapply(tax_ba, as.character)
tax_ba <- tax_ba[order(tax_ba$GTDB_phylum,tax_ba$GTDB_order),]


#tax_ba <- tax_ba[match( rownames(mat2),tax_ba$X),]


# # Get unique Phylum
# 
# cart_col <- paste0(c("blue","orange","red","brown","green", "purple","pink","grey",
#                      "turquoise","sand","taupe","kaki","harmo")
#                    ,".pal")
# 
# 
# phylum_col <- cart_col[1:length(unique(tax_ba$Phylum))]
# 
# unique(tax_ba$Phylum)
# unique(tax_ba$Class)
# 
# shades_num <- tax_ba %>%
#   group_by(Phylum) %>%
#   mutate(unique_types = n_distinct(Order))%>% select(Phylum,unique_types) %>% unique
# 
# library(paletteer)
# order_col <- mapply(function(x ,y){ paletteer_dynamic("cartography", !!x,y)},x = phylum_col ,y = shades_num$unique_types,
#                     SIMPLIFY = F)
# unlist(order_col)
# 
# tax_ba_sel <-  unique(tax_ba[,c("Phylum","Order")])
# 
# order_col <- as.vector(unlist(order_col))
# names(order_col) <- tax_ba_sel$Order
# 
# phylum_col <- gsub(".pal","",phylum_col)
# names(phylum_col) <- unique(tax_ba_sel$Phylum)



# tax_ba_sel$Order_col <- as.vector(unlist(order_col))
# tax_ba_sel$Phylum_col <- rep(gsub(".pal","",phylum_col),shades_num$unique_types)
# 
# order = as.character(tax_ba$Order)
# Phylum = as.character(tax_ba$Phylum)






# col_class <- assignCols(unique(Class),unique(Class))
# class_col <- col_class[[1]]
# names(class_col) <- col_class[[2]] 
# 
# library(RColorBrewer)
#    phy.col <-gsub(".pal","",rep(phylum_col, as.numeric(table(Phylum))))
# 
#    order.col <-  rep(as.vector(unlist(order_col)), as.numeric(table(order)))
#    
# phylum_col <- brewer.pal(length(unique(Phylum)),"Set1")
# names(phylum_col) <- unique(Phylum) 



#tax_ba <- tax_ba[match( rownames(mat2),tax_ba$X),]

#tmp_dt <- data.frame(mat2)

mat2 <- mat2[match(tax_ba$Species,rownames(mat2)),]


rownames(mat2) <- gsub("_"," ",rownames(mat2))
colnames(mat2) <- gsub("\\."," ",colnames(mat2))

ha1 = HeatmapAnnotation(Order = tax_ba$GTDB_order,
                       col = list(Order = order_col), which = "row")
ha2 = HeatmapAnnotation(Phylum = tax_ba$GTDB_phylum,
                        col = list(Phylum = phylum_col), which = "row")

split_rows <-  tax_ba$GTDB_phylum
split_rows <- factor(split_rows, levels= unique(split_rows))

#"Legends"
# lgd = Legend(labels = month.name[1:10], legend_gp = gpar(fill = 1:10), 
#              title = "foo", ncol = 3)
# draw(lgd)
# 
# tax_ba$Order
# 
# lgd <- 

ht1 = Heatmap(mat2, name = "Slope", column_title = NA, top_annotation = ha_column,
              clustering_distance_rows = "euclidean",
              split = split_rows, row_gap = unit(2, "mm"),border = TRUE,
              row_title_gp = gpar(fontsize = 10,fontface = "bold"),
              row_title_rot = 0,
              left_annotation = ha1,
              clustering_method_rows = "complete",row_names_side = "left", km=1, color_space = "LAB",
              row_dend_side="right",
              col=jet.colors,
              #heatmap_legend_param = list(),
              clustering_method_columns = "ward.D",
              width=2, row_names_max_width = unit(8, "cm"),show_column_names= T,
              row_names_gp = gpar(fontsize = 9), cluster_columns = T,cluster_rows = T,na_col="white")
#ht1
#ha1+ht1+ha2

#ht2 <- Heatmap(Class,col = class_col, name = "Class", width = unit(5, "mm"),show_heatmap_legend = T)
#ht3 <- Heatmap(Phylum,col = phylum_col, name = "Phylum", width = unit(5, "mm"),show_heatmap_legend = T)

#ht_ba <- ht1 + ht2 + ht3

pdf(paste(results_folder,'/',Sys.Date(),'_Heatmap_BA.pdf',sep=""),height = 13, width = 13)
draw(ht1)
dev.off()


# Presence absence
mat_logic  <- mat2 
mat_logic[mat_logic > 0] <- 1
mat_logic[mat_logic < 0] <- -1

mat_logic[mat_logic == 1] <- "Positive"
mat_logic[mat_logic == -1] <- "Negative"
mat_logic[mat_logic == 0] <- "NR"

# Also group by categories

id_order <- order(annotations$Class,annotations$Mechanism,annotations$Cat)

order_annot <- annotations[id_order,]

mat_logic <- mat_logic[,id_order]

colnames(mat_logic)
ha_column = HeatmapAnnotation(Class =  order_annot$Class,
                              Mechanism = order_annot$Mechanism,
                              Cat = order_annot$Cat,
                              col=list(Cat = c("Unconjugated BA" = "#b2182b","Primary Conjugated BA" = "#ef8a62",
                                               "Secondary Conjugated BA"= "#67a9cf","Secondary DeConjugated BA" = "#2166ac" ),
                                       Class=c("Primary BA" = "#1b9e77","Secon BA"= "#d95f02"),
                                       Mechanism = c("bai Dehydroxylation" = "#f26391","Epimerization"= "#6e3999")
                              ))


cols <- structure(c("#3b5998","#990000", "white"), names = c("Positive","Negative","NR"))
ht1 = Heatmap(mat_logic, name = "Slope", column_title = NA, top_annotation = ha_column,
              #clustering_distance_rows = "euclidean",
              split = split_rows, row_gap = unit(2, "mm"),border = TRUE,
              row_title_gp = gpar(fontsize = 10,fontface = "bold"),
              row_title_rot = 0,
              left_annotation = ha1,
              clustering_method_rows = "complete",row_names_side = "left", km=1, color_space = "LAB",
              #row_dend_side="right",
              col=cols,
              #heatmap_legend_param = list(),
              #clustering_method_columns = "ward.D",
              width=2, row_names_max_width = unit(8, "cm"),show_column_names= T,
              row_names_gp = gpar(fontsize = 9),
              #cluster_columns = T,cluster_rows = T,
              heatmap_legend_param = list(title = "Relation",
                                          at = c("Positive", "Negative"), labels = c("Positive", "Negative")),
              na_col="white")

pdf(paste(results_folder,'/',Sys.Date(),'_Heatmap_logical_ba.pdf',sep=""),height = 13, width = 13)
draw(ht1)
dev.off()




# SCFA
slope_dt_w_SCFA <- slope_dt_w[,names(slope_dt_w) %in% SCFA]

slope_dt_w_SCFA <- slope_dt_w_SCFA[rowSums(is.na(slope_dt_w_SCFA)) != ncol(slope_dt_w_SCFA), ]

mat<-as.matrix(slope_dt_w_SCFA)     
mat2<-scale(t(mat), scale = F, center = F)
mat2<-t(mat2)
# Setting NA to zero
mat2[is.na(mat2)] <- 0
library(circlize)
jet.colors <-   colorRamp2(c(-.05, 0, .05), c("blue", "white", "red"))

tax_dt <- match_tax_dt

tax_scfa <- tax_dt[tax_dt$Species %in% rownames(mat2),]
tax_scfa[] <- lapply(tax_scfa, as.character)
tax_scfa <- tax_scfa[order(tax_scfa$GTDB_phylum,tax_scfa$GTDB_order),]


#tax_scfa <- tax_scfa[match(rownames(mat2),tax_scfa$X),]


# Get unique Phylum

# cart_col <- paste0(c("blue","orange","red","brown","green", "purple","pink","wine","grey",
#                      "turquoise","sand","taupe","kaki","harmo")
#                    ,".pal")
# 
# 
# phylum_col <- cart_col[1:length(unique(tax_scfa$Phylum))]
# 
# unique(tax_scfa$Phylum)
# unique(tax_scfa$Class)
# 
# shades_num <- tax_scfa %>%
#   group_by(Phylum) %>%
#   mutate(unique_types = n_distinct(Order))%>% select(Phylum,unique_types) %>% unique
# 
# library(paletteer)
# order_col <- mapply(function(x ,y){ paletteer_dynamic("cartography", !!x,y)},x = phylum_col ,y = shades_num$unique_types,
#                     SIMPLIFY = F)
# unlist(order_col)
# 
# tax_scfa_sel <-  unique(tax_scfa[,c("Phylum","Order")])
# 
# order_col <- as.vector(unlist(order_col))
# names(order_col) <- tax_scfa_sel$Order
# 
# phylum_col <- gsub(".pal","",phylum_col)
# names(phylum_col) <- unique(tax_scfa_sel$Phylum)

mat2 <- mat2[match(tax_scfa$Species,rownames(mat2)),]


rownames(mat2) <- gsub("_"," ",rownames(mat2))
colnames(mat2) <- gsub("\\."," ",colnames(mat2))



colnames(mat2) <- gsub("X2.","2-",colnames(mat2))

ha1 = HeatmapAnnotation(Order = tax_scfa$GTDB_order,
                        col = list(Order = order_col), which = "row")
ha2 = HeatmapAnnotation(Phylum = tax_scfa$GTDB_phylum,
                        col = list(Phylum = phylum_col), which = "row")

split_rows <-  tax_scfa$GTDB_phylum
split_rows <- factor(split_rows, levels= unique(split_rows))


ht1 = Heatmap(mat2, name = "Slope", column_title = NA, 
              clustering_distance_rows = "euclidean",
              split = split_rows, row_gap = unit(2, "mm"),border = TRUE,
              row_title_gp = gpar(fontsize = 10,fontface = "bold"),
              row_title_rot = 0,
              left_annotation = ha1,
              clustering_method_rows = "complete",row_names_side = "left", km=1, color_space = "LAB",
              row_dend_side="right",
              col=jet.colors,
              #heatmap_legend_param = list(),
              clustering_method_columns = "ward.D",
              width=2, row_names_max_width = unit(8, "cm"),show_column_names= T,
              row_names_gp = gpar(fontsize = 9), cluster_columns = T,cluster_rows = T,na_col="white")
#ht1
#ha1+ht1+ha2

#ht2 <- Heatmap(Class,col = class_col, name = "Class", width = unit(5, "mm"),show_heatmap_legend = T)
#ht3 <- Heatmap(Phylum,col = phylum_col, name = "Phylum", width = unit(5, "mm"),show_heatmap_legend = T)

#ht_ba <- ht1 + ht2 + ht3

pdf(paste(results_folder,'/',Sys.Date(),'_Heatmap_SCFA.pdf',sep=""),height = 13, width = 10)
draw(ht1)
dev.off()


# Presence absence
mat_logic  <- mat2 
mat_logic[mat_logic > 0] <- 1
mat_logic[mat_logic < 0] <- -1

mat_logic[mat_logic == 1] <- "Positive"
mat_logic[mat_logic == -1] <- "Negative"
mat_logic[mat_logic == 0] <- "NR"

cols <- structure(c("#3b5998","#990000", "white"), names = c("Positive","Negative","NR"))
ht1 = Heatmap(mat_logic, name = "Slope", column_title = NA, 
              #clustering_distance_rows = "euclidean",
              split = split_rows, row_gap = unit(2, "mm"),border = TRUE,
              row_title_gp = gpar(fontsize = 10,fontface = "bold"),
              row_title_rot = 0,
              left_annotation = ha1,
              clustering_method_rows = "complete",row_names_side = "left", km=1, color_space = "LAB",
              #row_dend_side="right",
              col=cols,
              #heatmap_legend_param = list(),
              #clustering_method_columns = "ward.D",
              width=2, row_names_max_width = unit(8, "cm"),show_column_names= T,
              row_names_gp = gpar(fontsize = 9),
              #cluster_columns = T,cluster_rows = T,
              heatmap_legend_param = list(title = "Relation",
                                          at = c("Positive", "Negative"), labels = c("Positive", "Negative")),
              na_col="white")



pdf(paste(results_folder,'/',Sys.Date(),'_Heatmap_SCFA_logical.pdf',sep=""),height = 13, width = 10)
draw(ht1)
dev.off()
















