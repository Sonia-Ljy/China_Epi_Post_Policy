setwd("/home/soniali/Desktop/02_china_recom_github/0_raw_data/covsampler_300/")
library(ggtree)
library(ggplot2)
library(treeio)
library(stringr)
library(randomcoloR)
require(graphics)

# 读取树文件
tree <- read.tree("sub300ref_iqtree.fasta.treefile")
tree <- drop.tip(tree,'Reference') #删除游离的参考基因组样本
# 读取分组信息
group_file <- read.csv("meta_300.csv",header = T,row.names = 1)
##将地区替换为颜色
palette <- distinctColorPalette(length(unique(group_file$merged_lineage))+1)
palette = c("#82ADCF","#4F72B1","#F0B26D","#9C9A99","#5D9488","#9372A9","#F4D586","#C5432E","#B7D7E9")

class(palette)
# 按类分组
head(group_file)
groupInfo <- split(row.names(group_file), group_file$merged_lineage)
# 将分组信息添加到树中
tree <- groupOTU(tree, groupInfo)
# 输出样本的节点和分组信息
data=fortify(tree)

# 绘制进化树
p1 <- ggtree(tree,linetype=1,size=1.2,ladderize=T,color = "#3F3C3C")+ #xlim(0,0.0015) +
  geom_tippoint(aes(color = group), show.legend = T,size = 3) +
  xlim(0,0.0021)+geom_treescale()+
  scale_color_manual( values = palette)
p1

ggsave(filename = paste0("china_covsamper_compr_300_GTR_lin",".pdf"),plot = p1,device = "pdf",width = 5, height = 8)
