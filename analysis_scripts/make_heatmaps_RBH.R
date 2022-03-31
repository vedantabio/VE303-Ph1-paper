library(Biostrings)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(gtools)

# Read all the RBH files:
rbh_dt <- read.csv("../Data/RBH_results_VE303_Control.csv")

# Now read fasta file to extract the full names:
fasta_db <- readAAStringSet("../Data/SBA_Genes_dataset_v3_cdhit.fasta")
full_id <-  names(fasta_db)
fasta_id <-  gsub(" .*","",full_id)
id_dt <-  data.frame(fasta_id, full_id)

# Now combine the fasta_id from rbh to its full name 
library(tidyverse)
rbh_dt <-  rbh_dt %>%
             left_join(id_dt, by = c("B_id" = "fasta_id"))

ba_Genes_new <-rbh_dt
ba_Genes_new$Genome<- gsub("GCA_000156055", "C. hiranonis DSM 13275", ba_Genes_new$Genome)
ba_Genes_new$Genome <- gsub("GCA_000156515", "C. hylemonae DSM 15053", ba_Genes_new$Genome)
ba_Genes_new$Genome <- gsub("GCA_004558675", "C. scindens W0P25", ba_Genes_new$Genome)
ba_Genes_new$Genome <- gsub("GCA_000024265", "E. lenta DSM 2243", ba_Genes_new$Genome)
pos_control <-  c("C. hiranonis DSM 13275","C. hylemonae DSM 15053","C. scindens W0P25","E. lenta DSM 2243")

to_avoid <- c(pos_control ,grep("VE303", ba_Genes_new$Genome,value = T))
idx_to_avoid <- which(ba_Genes_new$Genome %in% to_avoid)

ba_Genes_new$ba <- gsub("OS.*","",sub(".*? (.+)", "\\1", ba_Genes_new$full_id))
ba_Genes_new$ba <- gsub("Dehydrogenases with different specificities.*", "Dehydrogenases with different specificities", ba_Genes_new$ba)

mat_new <-  ba_Genes_new[,c("Genome","pident","ba")]
mat1_new <- mat_new[!duplicated(mat_new),]

mat2_new <- mat1_new %>%
  #  select(pident,ba,ID)%>%
  group_by(Genome,ba) %>%
  filter(pident == max(pident))%>% data.frame()

mat_w_new <-reshape(mat2_new, idvar = c("Genome"), timevar = "ba", direction = "wide")
names(mat_w_new) <- gsub("pident.","",names(mat_w_new))
mat_w_new[is.na(mat_w_new)] <-  0
rownames(mat_w_new) <- mat_w_new$Genome
mat_w_new$Genome <- NULL

Group <- rep("Control", nrow(mat_w_new))
Group[grep('VE3', rownames(mat_w_new))] <- "VE303"
group_col <- c("red","darkblue")
names(group_col)<-c("Control","VE303")
group_col

library(ComplexHeatmap)
ha3 = HeatmapAnnotation(Group = Group,
                        col = list(Group = group_col), which = "row")
split_rows <-  Group
split_rows <- factor(split_rows, levels= unique(split_rows))

color.palette <- colorRampPalette(c("white", "blue", "#007FFF", "cyan",
                                    "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(20)


split_cols <-  colnames(mat_w_new)
split_cols <-  rep("Other",length(split_cols))
split_cols[grep("Bai|Bile acid",colnames(mat_w_new))] <- "Bai_genes"
split_cols[grep("3a|3-alpha|3alpha",colnames(mat_w_new))] <- "3-a-HSDH"
split_cols[grep("3b|3-beta|3beta",colnames(mat_w_new))] <- "3-b-HSDH"
split_cols[grep("7a|7-alpha|7alpha",colnames(mat_w_new))] <- "7-a-HSDH"
split_cols[grep("7b|7-beta|7beta",colnames(mat_w_new))] <- "7-b-HSDH"
split_cols[grep("Choloylglycine hydrolase|choloylglycine hydrolase",colnames(mat_w_new))] <- "BSH"

library(gtools)
mixedsort(unique(split_cols))

split_cols <- factor(split_cols, mixedsort(unique(split_cols)))

ht2 = Heatmap(as.matrix(mat_w_new), name = "Identity", column_title = NA,
              #top_annotation = ha,
              show_parent_dend_line = F,
              show_row_dend = F,
              show_column_dend = F,
              #clustering_distance_rows = "euclidean",
              split = split_rows,
              column_split = split_cols,
              column_title_side = "top",
              column_title_rot = 90,
              cluster_column_slices = F,
              row_gap = unit(2, "mm"),border = TRUE,
              row_title_gp = gpar(fontsize = 10,fontface = "bold"),
              row_title_rot = 0,
              left_annotation = ha3,
              row_names_side = "left", km=1, color_space = "LAB",
              #row_dend_side="right",
              col=color.palette ,
              #heatmap_legend_param = list(),
              #clustering_method_columns = "ward.D",
              width=2, show_column_names= T,
              row_names_gp = gpar(fontsize = 9),
              row_names_max_width = max_text_width(rownames(mat_w_new), gp = gpar(fontsize = 12)),
              column_names_max_height = max_text_width(colnames(mat_w_new), gp = gpar(fontsize = 12)),
              column_names_side = "bottom",
              na_col="white")

# pdf('Heatmap_RBH_blast_ba_genes_v3.pdf',height = 10, width = 12)
# draw(ht2,padding = unit(c(2, 2,15, 2), "mm"))
# dev.off()


pdf('Heatmap_RBH_blast_ba_genes.pdf',height = 10, width = 14)
draw(ht2)
dev.off()

