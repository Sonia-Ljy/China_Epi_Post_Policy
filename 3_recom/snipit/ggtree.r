
library(ggtree)
library(treeio)
library(ggplot2)
library(tidyverse)


setwd('/home/soniali/Desktop/02_china_recom/2_XCN/snipit/left/')
tree <- read.newick("FR.1.1_XCN_EG.5.1.1_masked_left.fasta.treefile",node.label = 'support')
tree <- drop.tip(tree,'Reference') #删除游离的参考基因组样本

# MRCA(tree, "EPI_ISL_18101315_FE.1", "EPI_ISL_18101313_FE.1")
# MRCA(tree, "C_AA030871.1_FR.1.1","C_AA037154.1_FR.1.1","C_AA037256.1_FR.1.1","C_AA037234.1_FR.1.1","C_AA037204.1_FR.1.1","C_AA037540.1_FR.1.1","C_AA038331.1_FR.1.1","C_AA038701.1_FR.1.1","C_AA036799.1_FR.1.1","C_AA036729.1_FR.1.1","EPI_ISL_17994535_FR.1.1","EPI_ISL_18060821_FR.1.1","EPI_ISL_18060751_FR.1.1","EPI_ISL_18074541_FR.1.1","EPI_ISL_18074511_FR.1.1","EPI_ISL_18074780_FR.1.1","EPI_ISL_18074563_FR.1.1","EPI_ISL_18074461_FR.1.1","EPI_ISL_18105577_FR.1.1","EPI_ISL_18108456_FR.1.1")
# MRCA(tree, "C_AA038410.1_recom","C_AA044601.1_recom","C_AA047643.1_recom","EPI_ISL_18284946_recom","EPI_ISL_18289774_recom","EPI_ISL_18289734_recom","EPI_ISL_18289815_recom","EPI_ISL_18401744_recom","EPI_ISL_18105656_recom")

p <-ggtree(tree,color='#C0C0C0',linetype=1,size=1.5,ladderize=T) + xlim(0,0.003) + #  0.03  0.004
  geom_tiplab(hjust = -0.05,size=3,fontface="italic",align = F,angle = 0,layout_dendrogram()) +#
  # geom_point2(aes(subset=!isTip),fill="red",shape=21,size=4) + #A62B19
  # geom_text2(aes(subset=!isTip,label=support,color=support, hjust=-0.6),size=4) +
  # scale_color_gradient(high='red', low='Navy') +
  # geom_cladelabel(node=61, label="zygote", offset = 10,barsize = 1,color=c("black","red2"), align=TRUE) +
  theme(legend.position='left') +
  geom_treescale()#+
  # geom_highlight(node = 61,fill = "#E4C4CB")+
  # geom_highlight(node = 33,fill = "#9FC6C4")
  
p

ggsave(filename = paste0("FE.1_recom_FR.1.1_left_v1",".pdf"),plot = p,device = "pdf",width = 8, height = 5)
#######

setwd('/home/soniali/Desktop/02_china_recom/2_XCN/snipit/right/')
tree <- read.newick("FR.1.1_XCN_EG.5.1.1_masked_right.fasta.treefile",node.label = 'support')
tree <- drop.tip(tree,'Reference') #删除游离的参考基因组样本

# MRCA(tree, "EPI_ISL_18101315_FE.1", "EPI_ISL_18101313_FE.1")
# MRCA(tree, "C_AA030871.1_FR.1.1","C_AA037154.1_FR.1.1","C_AA037256.1_FR.1.1","C_AA037234.1_FR.1.1","C_AA037204.1_FR.1.1","C_AA037540.1_FR.1.1","C_AA038331.1_FR.1.1","C_AA038701.1_FR.1.1","C_AA036799.1_FR.1.1","C_AA036729.1_FR.1.1","EPI_ISL_17994535_FR.1.1","EPI_ISL_18060821_FR.1.1","EPI_ISL_18060751_FR.1.1","EPI_ISL_18074541_FR.1.1","EPI_ISL_18074511_FR.1.1","EPI_ISL_18074780_FR.1.1","EPI_ISL_18074563_FR.1.1","EPI_ISL_18074461_FR.1.1","EPI_ISL_18105577_FR.1.1","EPI_ISL_18108456_FR.1.1")
# MRCA(tree, "C_AA038410.1_recom","C_AA044601.1_recom","C_AA047643.1_recom","EPI_ISL_18284946_recom","EPI_ISL_18289774_recom","EPI_ISL_18289734_recom","EPI_ISL_18289815_recom","EPI_ISL_18401744_recom","EPI_ISL_18105656_recom")

p <-ggtree(tree,color='#C0C0C0',linetype=1,size=1.5,ladderize=T) + xlim(0,0.0045) + #  0.03  0.004
  geom_tiplab(hjust = -0.05,size=3,fontface="italic",align = F,angle = 0,layout_dendrogram()) +#
  # geom_point2(aes(subset=!isTip),fill="red",shape=21,size=4) + #A62B19
  # geom_text2(aes(subset=!isTip,label=support,color=support, hjust=-0.6),size=4) +
  # scale_color_gradient(high='red', low='Navy') +
  # geom_cladelabel(node=61, label="zygote", offset = 10,barsize = 1,color=c("black","red2"), align=TRUE) +
  theme(legend.position='top') +
  geom_treescale() #+
  # geom_highlight(node = 61,fill = "#E4C4CB")+
  # geom_highlight(node = 33,fill = "#9FC6C4")
p

ggsave(filename = paste0("FE.1_recom_FR.1.1_right_v1",".pdf"),plot = p,device = "pdf",width = 8, height = 5)

