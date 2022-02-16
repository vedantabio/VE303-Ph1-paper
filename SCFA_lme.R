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

mainDir <- "../results"
subDir <- paste0("SCFA_BA_lme/",gsub("-","_",Sys.Date()))
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")



# Read SCFA 
phy_scfa <-  readRDS("../data/processed_data/phy_met.rds")
phy_scfa_m <- psmelt(phy_scfa)
head(phy_scfa_m)
phy_scfa_m <- phy_scfa_m[,c("Sample","OTU","cohort_id_long","Subject.ID","desc.ID","Collection.Day.Norm","Timepoint.Calc","Abundance")]

# Plot VE303 Total Abundance on the x and scfa on the y from each cohort 
phy_mic <-  readRDS("../data/processed_data/phy_mic_abs.rds")
phy_sel <- prune_samples(sample_names(phy_scfa), phy_mic)
phy_sel <- prune_taxa(taxa_sums(phy_sel)>0, phy_sel)

# Select only VE303 genes
ve303_st <- grep("VE303",taxa_names(phy_sel), value = T)
phy_ve303  <- prune_taxa(ve303_st, phy_sel)  

tot_ve303 <- tax_glom(phy_ve303, "superkingdom_name")  
taxa_names(tot_ve303) <- "VE-303"
ve303_dt  <-  t(data.frame(otu_table(tot_ve303)))
ve303_dt    <-  data.frame(Sample = rownames(ve303_dt), tot_abun =  as.vector(ve303_dt))
ve303_dt$Sample <-  gsub("X","",ve303_dt$Sample)

scfa_dt <-  merge( x =phy_scfa_m, y = ve303_dt, by =  "Sample" )
unique(scfa_dt$OTU)
scfa <- unique(scfa_dt$OTU)
names(scfa_dt)

# LME on SCFA to find if the SCFA recover after addition of VE303
# Remove time points lower than 0 
scfa_dt$Treatment <- "Pre_Vanco"
scfa_dt$Treatment[scfa_dt$Collection.Day.Norm < 0] <- "Pre_Vanco"
scfa_dt$Treatment[scfa_dt$Collection.Day.Norm >= 0 & scfa_dt$Collection.Day.Norm <= 6] <- "Vanco"
scfa_dt$Treatment[scfa_dt$Collection.Day.Norm >6 ] <- "Post_Vanco"
# Remove Cohort 6 for Pre and Post differential analysis
# Vanco was not administered in Cohort 6 
scfa_dt <- scfa_dt[- which(scfa_dt$cohort_id_long== "Cohort 6"), ] 

# Which SCFA goes up in Vanco Treatment?
# scfa_dt <- scfa_dt[-which(scfa_dt$Treatment == "Post_Vanco"),]
result_lme <-  list()
scfa <- unique(scfa_dt$OTU)
scfa_var <- scfa[15]
for (scfa_var in scfa){
  
  scfa_sel <- scfa_dt[scfa_dt$OTU %in% scfa_var,]
  
  # Remove the species with NA absolute abundance
  scfa_sel <- scfa_sel[complete.cases(scfa_sel[, "Abundance"]),]
  scfa_sel$Treatment <- factor(scfa_sel$Treatment, levels = c("Pre_Vanco","Vanco","Post_Vanco"))
  scfa_sel$tot_abun <- log10(scfa_sel$tot_abun+ 10^-7) 
  #scfa_sel$tot_abun <-  sqrt(scfa_sel$tot_abun) 
  scfa_sel$Subject.ID <- factor(scfa_sel$Subject.ID)
  scfa_sel$Abundance <- log10(scfa_sel$Abundance)
  
  library("nlme")
  fml <- as.formula( paste( "Abundance", "~", paste(c("Treatment + tot_abun"), collapse="+") ) )
  mod_bac <- lme(fml,random = ~1|cohort_id_long/Subject.ID, scfa_sel)
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$scfa_var <- scfa_var
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[scfa_var]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
# Remove intercept
result_lme_dt <- result_lme_dt[result_lme_dt$Var != "(Intercept)",]
result_lme_dt$p_adj <- p.adjust(result_lme_dt$pval, method = "BH")
result_lme_dt <- result_lme_dt[order(result_lme_dt$p_adj,decreasing = F),]

result_lme_dt <- result_lme_dt[,c("Value","scfa_var","Var","p_adj" )]
dt_w <- reshape(result_lme_dt, idvar = "scfa_var", timevar = "Var", direction = "wide")

dt_w <- dt_w[,mixedsort(names(dt_w))]
rownames(dt_w) <- dt_w$scfa_var

dt_w_BA <-  dt_w[grep("Acid",dt_w$scfa_var),]
# Categorize the BA
dt_w_BA$Cat <- "Unconjugated BA"
dt_w_BA$Cat[dt_w_BA$scfa_var %in% 
              c("Glycocholic.Acid","Taurocholic.Acid",
                "Glycochenodeoxycholic.Acid","Taurochenodeoxycholic.Acid")]<- "Primary Conjugated BA"
dt_w_BA$Cat[dt_w_BA$scfa_var %in% 
              c("Glycodeoxycholic.Acid","Taurodeoxycholic.Acid",
                "Glycolithocholic.Acid","Taurolithocholic.Acid",
                "Glycoursodeoxycholic.Acid","Tauroursodeoxycholic.Acid")]<- "Secondary Conjugated BA"

dt_w_BA$Cat[dt_w_BA$scfa_var %in% 
              c("Deoxycholic.Acid",
                "Lithocholic.Acid","Ursodeoxycholic.Acid")]<- "Secondary DeConjugated BA"





write.csv(dt_w_BA,paste0(results_folder,"/LME_results_ve303_no_cohort_6_add_vanco_BA.csv"))

dt_w_scfa <-  dt_w[!dt_w$scfa_var %in% rownames(dt_w_BA),]
dt_w_scfa$Type <-  "SCFA"

write.csv(dt_w_scfa,paste0(results_folder,"/LME_results_ve303_no_cohort_6_add_vanco_SCFA.csv"))


