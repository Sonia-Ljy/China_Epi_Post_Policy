library(ggtree)
library(treeio)
library(ggplot2)
library(tidyverse)

setwd('/home/soniali/Desktop/02_china_recom_renew/3_recom/snipit/')
tree <- read.tree("left/FR.1.1_EG.5.1.1_masked_left.fasta.nwk")# read.newick
tree <- drop.tip(tree,'Reference')

p <-ggtree(tree,color='#C0C0C0',linetype=1,size=1.5,ladderize=T) + xlim(0,0.003) + #  0.03  0.004
  geom_tiplab(hjust = -0.05,size=3,fontface="italic",align = F,angle = 0,layout_dendrogram())+#
  geom_treescale()+geom_nodelab(hjust=-0.05,size=2)
p
ggsave(filename = paste0("FE.1_recom_FR.1.1_left_v3",".pdf"),plot = p,device = "pdf",width = 8, height = 5)
