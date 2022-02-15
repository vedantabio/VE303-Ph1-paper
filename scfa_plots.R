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

# Create a directory to store analysis
mainDir <- "../results"
subDir <- "SCFA_plots"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")


# Read SCFA 
scfa_dt <-  read.csv("../data/Metabolomics/Pharmaron_VE303-01_Metabolon_30NOV2018.csv")
scfa_dt$Pharmaron.Sample.Barcode <- as.character(scfa_dt$Pharmaron.Sample.Barcode)
names(scfa_dt)

scfa_dt$Cohort.ID
# Remove Subject 160
scfa_dt <- scfa_dt[scfa_dt$Subject.ID != 160,]
# Remove Cohort 8
scfa_dt <- scfa_dt[scfa_dt$Cohort.ID != "8",]

# Fix the Timepoint
scfa_dt$Time <- as.character(scfa_dt$Timepoint...day.of.treatment.)
unique(scfa_dt$Time)
scfa_dt$Time <- gsub("DAY","Day",scfa_dt$Time)
scfa_dt$Time <- gsub("SCREENING","Screening",scfa_dt$Time)
# Remove Unscheduled
scfa_dt <- scfa_dt[scfa_dt$Time != "Unscheduled",]
scfa_dt$Time <-  factor( scfa_dt$Time, c("Screening",mixedsort(unique(scfa_dt$Time)[-1])))

# Fix Cohort ID
unique(scfa_dt$Cohort.ID)
scfa_dt$Cohort.ID <- gsub("Sentinel","1",scfa_dt$Cohort.ID)
scfa_dt$Cohort.ID[-grep("Vanco",scfa_dt$Cohort.ID)] <-  paste0("Cohort ",scfa_dt$Cohort.ID[-grep("Vanco",scfa_dt$Cohort.ID)])
unique(scfa_dt$Cohort.ID )
scfa_dt$Cohort.ID <- factor(scfa_dt$Cohort.ID, levels = c("Vanco",paste0("Cohort ",1:6)))

unique(scfa_dt$Time )
unique(scfa_dt$Analyte )

# Filter Analyte
scfa_dt$Analyte <-  as.character(scfa_dt$Analyte)
scfa_names  <- c("2-Methylbutyrate","Acetate","Butyrate","Hexanoate","Isobutyrate","Isovalerate","Propionate","Valerate")  
scfa_dt <-  scfa_dt[scfa_dt$Analyte %in% scfa_names,]


# Replace BLOQ with NA
unique(scfa_dt$Result)
scfa_dt$Abundance <-  as.numeric(as.character(scfa_dt$Result))
scfa_dt$Abundance[grep("BLOQ",scfa_dt$Result)] <-  scfa_dt$LLOQ[grep("BLOQ",scfa_dt$Result)]
names(scfa_dt)
# Sum the abundances per Sample 
sum_abun_dt <- scfa_dt %>% 
  group_by(Pharmaron.Sample.Barcode,Subject.ID,Cohort.ID ,Time) %>% 
  summarise(sum_abun = sum(Abundance))

# Now compute median of the sum:
med_abun_det <-  sum_abun_dt %>%
  group_by(Cohort.ID ,Time)%>%
  summarise(median_abun = median(sum_abun),
            mean_abun = mean(sum_abun),
            Total_samples = n())

coh.cols <- c("Vanco" = "#0487e3", "Cohort 1" = "#dc2d00", "Cohort 2" = "#5b0076", "Cohort 3" = "#629449", "Cohort 4" = "#561600", "Cohort 5" = "#004631", "Cohort 6" = "#b87800", "Cohort 8" = "#b97696")

pdf(paste0(results_folder,"/Total_Conc_SCFA.pdf"),width = 20,height = 4)

plot_abun <- ggplot(sum_abun_dt, aes( y= log10(sum_abun), x= Time,
                                      fill = Cohort.ID,color = Cohort.ID)) +
  geom_point(aes( group=Subject.ID),size = 2,alpha = 0.15)+
  geom_line(aes( group=Subject.ID),alpha = 0.15)+
  geom_point(data = med_abun_det, aes( y= log10(median_abun), x= Time,
                                       color = Cohort.ID,group = Cohort.ID),size = 3,alpha = 0.7)+
  geom_line(data = med_abun_det, aes( y= log10(median_abun), x= Time,
                                      color = Cohort.ID,group = Cohort.ID),size = 1,alpha = 0.7)+
  #geom_boxplot(alpha = 0.5,outlier.colour = NA)+
  scale_fill_manual(name = "Cohort",values = coh.cols)+
  scale_color_manual(name = "Cohort",values = coh.cols)+
  facet_wrap(~Cohort.ID,scales = "free_x",nrow = 1)+
  #scale_y_log10()+
  theme_bw()+
  ylab("Log10 (Total Conc ug/mg)")+
  xlab("Time")+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"))
print(plot_abun)      

plot_abun <- ggplot(sum_abun_dt, aes( y= sum_abun, x= Time,
                                      fill = Cohort.ID,color = Cohort.ID)) +
  geom_point(aes( group=Subject.ID),size = 2,alpha = 0.15)+
  geom_line(aes( group=Subject.ID),alpha = 0.15)+
  geom_point(data = med_abun_det, aes( y= median_abun, x= Time,
                                       color = Cohort.ID,group = Cohort.ID),size = 3,alpha = 0.7)+
  geom_line(data = med_abun_det, aes( y= median_abun, x= Time,
                                      color = Cohort.ID,group = Cohort.ID),size = 1,alpha = 0.7)+
  #geom_boxplot(alpha = 0.5,outlier.colour = NA)+
  scale_fill_manual(name = "Cohort",values = coh.cols)+
  scale_color_manual(name = "Cohort",values = coh.cols)+
  facet_wrap(~Cohort.ID,scales = "free_x",nrow = 1)+
  scale_y_log10()+
  theme_bw()+
  ylab("Total Conc (ug/mg)")+
  xlab("Time")+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"))
print(plot_abun)      
dev.off()
