library(ggplot2)
library(ggsci)
library(forcats)

setwd("/home/soniali/Desktop/02_china_recom_github/1_epi/")
df <- read.csv("china_lineage_num_draw.csv", header = TRUE)
df$Nextstrain_clade <- factor(df$Nextstrain_clade,levels=c("22B",'22D',"22F","23A","23B","23E","Others"))
df1 = aggregate(df$count, by = list(df$Pangolin_lineage), FUN = sum)
df1$per1 = df1$x / sum(df1$x)
for (i in seq(nrow(df1), 1)) {
  if (i == nrow(df1)) {
    df1$per.y1[i] = df1$per1[i] / 2
  }else{
    df1$per.y1[i] = sum(df1$per1[(i + 1):nrow(df1)]) + df1$per1[i] / 2
  }
}

df1$label1 = paste(df1$Group.1,'(',round(df1$per1*100, 2),'%',')', sep = '')
df = merge(df, df1[,c(1,3,4,5)], by.x = 'Pangolin_lineage', by.y = 'Group.1')
df2 = aggregate(df$count, by = list(df$Nextstrain_clade), FUN = sum)
df2$per2 = df2$x / sum(df2$x)

for (i in seq(nrow(df2), 1)) {
  if (i == nrow(df2)) {
    df2$per.y2[i] = df2$per2[i] / 2
  }else{
    df2$per.y2[i] = sum(df2$per2[(i + 1):nrow(df2)]) + df2$per2[i] / 2
  }
}

df2$label2 = paste(df2$Group.1,'(',round(df2$per2*100, 2),'%',')', sep = '')
df = merge(df, df2[,c(1,3,4,5)], by.x = 'Nextstrain_clade', by.y = 'Group.1')

df$y = c("a","a","a","a","a","a","a","a","a")
df$x = c("b","b","b","b","b","b","b","b","b")

p1 <-ggplot(df) +
  # 绘制柱状图
  geom_bar(aes(y,
               ifelse(Nextstrain_clade %in% c('22D',"22F","23B","23E","Others"), per2, per2/2), #"22B""22D""22F" "23A" "23B""23E" "Others"
               fill = Nextstrain_clade),
           stat = 'identity', width = 1.3) +
  # 添加标签
  geom_text(aes(1.25, as.numeric(per.y2), 
                label = label2),
            size =2.5, color = 'black') +
  # 绘制柱状图
  geom_bar(aes(x, per1, fill = Pangolin_lineage), 
           stat = 'identity', width = .6, color = 'white') +
  # 添加标签
  geom_text(aes(2, as.numeric(per.y1),label = label1),
            size = 2.5, color = 'black') +
  # 设置Y轴刻度
  scale_y_continuous(labels = scales::percent) +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_igv(alpha = 0.6) +
  theme(legend.position = 'none')

p1
ggsave(filename = "china_lineage_num_plot.pdf",width = 8,height = 8)
