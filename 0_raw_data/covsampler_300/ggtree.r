setwd("/home/soniali/Desktop/02_china_recom_renew/0_raw_data/covsampler_300/")
# 加载R包
library("ggtree")
library("ggplot2")
library("treeio")
library(stringr)
library(randomcoloR)
require(graphics)


# 读取树文件（ggtree针对不同工具生成的树文件有不同的函数进行读取，如没有对应的函数，可以使用通用的办法读取）
tree <- read.tree("sub300ref_iqtree.fasta.nwk")
# tree<-read.newick("sub300ref_iqtree.fasta.nwk",node.label = "support")
# tree
as_tibble(tree)
tree <- drop.tip(tree,'Reference') #删除游离的参考基因组样本
# 读取分组信息
group_file <- read.csv("meta_300.csv",header = T,row.names = 1)
##将地区替换为颜色
palette <- distinctColorPalette(length(unique(group_file$merged_lineage))+1)
# palette=c("#d62728","#1f77b4",  "#2ca02c", "#ff7f0e","#9467bd", "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f","#800080")
palette = c("#82ADCF","#4F72B1","#F0B26D","#9C9A99","#5D9488","#9372A9","#F4D586","#C5432E","#B7D7E9")

class(palette)
# 按类分组
head(group_file)
groupInfo <- split(row.names(group_file), group_file$merged_lineage)
# 将分组信息添加到树中
tree <- groupOTU(tree, groupInfo)
# days_name <- unique(group_file$merged_lineage)
# 输出样本的节点和分组信息
data=fortify(tree)

# tree@data$support1<-ifelse(tree@data$support<75,NA,tree@data$support)

# 绘制进化树
p1 <- ggtree(tree,linetype=1,size=1.2,ladderize=T,color = "#3F3C3C")+ #xlim(0,0.0015) +
  geom_tippoint(aes(color = group), show.legend = T,size = 3) +
  xlim(0,0.0021)+geom_treescale()+
  # theme(legend.position='none')+
  scale_color_manual( values = palette)+
  # geom_text(aes(subset=!isTip,label=support,color=support, hjust=-0.5),size=4)
  # geom_nodelab(hjust=-0.05,size=4) 
  # geom_tiplab(size=1.5)+
  geom_nodelab(hjust=-0.05,size=2)#+
  # geom_nodepoint(size = color="#b5e521", alpha = 0.5, size = 7)
  # geom_nodepoint(aes( x=x-0.5))
  
p1

ggsave(filename = paste0("china_covsamper_compr_300_GTR_lin2",".pdf"),plot = p1,device = "pdf",width = 5, height = 8)
###############

################################################################################################################
################################################################################################################

palette <- distinctColorPalette(length(unique(group_file$province))+1)
# palette=c("#d62728","#1f77b4",  "#2ca02c", "#ff7f0e","#9467bd", "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f","#800080")
class(palette)
# 按类分组
head(group_file)
groupInfo <- split(row.names(group_file), group_file$province)
# 将分组信息添加到树中
tree <- groupOTU(tree, groupInfo)
# days_name <- unique(group_file$merged_lineage)
# 输出样本的节点和分组信息
data=fortify(tree)
# write.table(data, file = "epi_node_lineages.txt", append = FALSE, quote = FALSE, sep = "\t",
#             na = "NA", row.names = FALSE,
#             col.names = TRUE)
# 绘制进化树
p1 <- ggtree(tree,linetype=1,size=1.2,ladderize=T,color = "#3F3C3C")+ #xlim(0,0.0015) +
  geom_tippoint(aes(color = group), show.legend = T,size = 3) +
  xlim(0,0.0021)+geom_treescale()+
  # theme(legend.position='none')+
  scale_color_manual( values = palette)

p1

ggsave(filename = paste0("china_covsamper_compr_300_GTR",".pdf"),plot = p1,device = "pdf",width = 5, height = 8)

###############################################################################################################
################################################################################################################
library(RColorBrewer)
display.brewer.all()
require(viridis)
palette <- viridis(length(unique(group_file$date)))
# palette <- distinctColorPalette(length(unique(group_file$date))+1)
# palette=c("#d62728","#1f77b4",  "#2ca02c", "#ff7f0e","#9467bd", "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f","#800080")

class(palette)
# 按类分组
head(group_file)
groupInfo <- split(row.names(group_file), group_file$date)
# 将分组信息添加到树中
tree <- groupOTU(tree, groupInfo)
# days_name <- unique(group_file$merged_lineage)
# 输出样本的节点和分组信息
data=fortify(tree)
# write.table(data, file = "epi_node_lineages.txt", append = FALSE, quote = FALSE, sep = "\t",
#             na = "NA", row.names = FALSE,
#             col.names = TRUE)
# 绘制进化树
p1 <- ggtree(tree,linetype=1,size=1.2,ladderize=T,color = "#3F3C3C")+ #xlim(0,0.0015) +
  geom_tippoint(aes(color = group), show.legend = T,size = 3) +
  xlim(0,0.0021)+geom_treescale()+
  theme(legend.position='none')+
  scale_color_manual( values = palette)

p1
ggsave(filename = paste0("china_covsamper_compr_300_GTR_date",".pdf"),plot = p1,device = "pdf",width = 5, height = 8)

p1 <- ggtree(tree,linetype=1,size=1.2,ladderize=T,color = "#3F3C3C")+ #xlim(0,0.0015) +
  geom_tippoint(aes(color = group), show.legend = T,size = 3) +
  xlim(0,0.0021)+geom_treescale()+
  # theme(legend.position='none')+
  scale_color_manual( values = palette)

p1

ggsave(filename = paste0("china_covsamper_compr_300_GTR_date2",".pdf"),plot = p1,device = "pdf",width = 12, height = 8)
