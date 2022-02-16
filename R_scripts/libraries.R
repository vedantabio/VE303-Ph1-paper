library(patchwork)
library(tidyverse)
library(openxlsx)
library(knitr)
library(viridis)
library(here)

# Custom Plotting parameters
pretty_plot = theme_classic() + 
  theme(
    text = element_text(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    axis.line.y.left = element_line(color = "black"),
    # panel.border = element_blank(),
    # strip.background = element_blank(),
    strip.text = element_text(size = 12),
    # panel.grid.major = element_blank(), 
    # panel.grid.minor = element_blank(),
    plot.title = element_text(size=15, face="bold"),
    axis.title=element_text(size=12,face="bold"),
    axis.text.y = element_text(size=12, color="#000000"),
    axis.text.x = element_text(size=12, color="#000000")
  )