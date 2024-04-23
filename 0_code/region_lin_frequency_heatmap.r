library(pheatmap)
library(RColorBrewer)
library(circlize)
library(ggplot2)
dirpath = "/home/soniali/Desktop/02_china_recom_github/1_epi/"
setwd(dirpath)

data <- read.csv(paste0(dirpath,"region_lin_frequency_heatmap.csv"), header = TRUE, row.names = 1)
bk <- c(seq(0,0.2,by=0.01),seq(0.21,1,by=0.01))
province_group_color = c("A" = "#C59D83", "B" = "#D8BFAB","C" = "#FBEBDA","D" = "#B9B8A5","E" = "#A5A48E","F" = "#6C705E")
ann_colors = list(Group = province_group_color)


data1 <- t(data)
p1 <- pheatmap(data1,
               scale = "none",cluster_cols=TRUE,clustering_distance_cols = "manhattan",display_numbers = matrix(ifelse(data1 > 0.5, "*", ""), nrow(data1)),
               cellwidth = 12, cellheight = 15,cutree_cols=8,
               # annotation_col = groufile,annotation_colors = ann_colors,
               color = c(colorRampPalette(colors = c("#5774E1","white"))(length(seq(0,0.2,by=0.01))),colorRampPalette(colors = c("white","#B70E29"))(length(seq(0.21,1,by=0.01)))),
               legend_breaks=seq(0,1,0.1),angle_col = 45,fontsize_row = 8,fontsize_col = 7,
               breaks=bk)
p1
# ggsave(filename = paste0(paste0(dirpath,"region_lin_frequency_heatmap0"),'.pdf'),plot = p1,device = "pdf",width = 9, height = 5)

groufile <- read.csv(paste0(dirpath,"region_lin_freq_group.csv"), header = TRUE, row.names = 1)
groufile <- data.frame(groufile)
bk <- c(seq(0,0.2,by=0.01),seq(0.21,1,by=0.01))
province_group_color = c("A" = "#C59D83", "B" = "#D8BFAB","C" = "#FBEBDA","D" = "#B9B8A5","E" = "#A5A48E","F" = "#6C705E")
ann_colors = list(Group = province_group_color)
data1 <- t(data)
p11 <- pheatmap(data1,
               scale = "none",cluster_cols=TRUE,clustering_distance_cols = "manhattan",display_numbers = matrix(ifelse(data1 > 0.5, "*", ""), nrow(data1)),
               cellwidth = 12, cellheight = 15,cutree_cols=8,
               annotation_col = groufile,annotation_colors = ann_colors,
               color = c(colorRampPalette(colors = c("#5774E1","white"))(length(seq(0,0.2,by=0.01))),colorRampPalette(colors = c("white","#B70E29"))(length(seq(0.21,1,by=0.01)))),
               legend_breaks=seq(0,1,0.1),angle_col = 45,fontsize_row = 8,fontsize_col = 7,
               breaks=bk)
p11
ggsave(filename = paste0(paste0(dirpath,"region_lin_frequency_heatmap"),'.pdf'),plot = p11,device = "pdf",width = 9, height = 5)

library(maps)
library(ggplot2)
library(rmapshaper)
library(hchinamap)

Rdata<-read.csv("/home/soniali/Desktop/02_china_recom_github/1_epi/China_except_sh.csv")
p_china_except_sh <- hchinamap(name = Rdata$name, 
          value = Rdata$value,
          width = "100%",
          height = "400px",
          minColor = "#F0F0F0", 
          maxColor = "#B8142E",
          region = "China"
)
# ggsave(filename = paste0(paste0(dirpath,"China_except_sh.jpeg"),'.jpeg'),plot = p_china_except_sh,device = "pdf")
ggsave("China_except_sh.jpeg", plot = p_china_except_sh, width = 10, height = 8, units = "in")


library(mapchina)
library(sysfonts)
library(showtextdb)
library(showtext)
library(tidyverse)
library(sf)
library(ggrepel)

code <- china$Code_Province
pro_eng <- china$Name_Province
unique(pro_eng)
values <- c("安徽省","北京市","重庆市","福建省","甘肃省","广东省","广西壮族自治区","贵州省","海南省","河北省","黑龙江省","河南省","香港省","湖北省","湖南省","内蒙古自治区","江苏省","江西省","吉林省","辽宁省","宁夏回族自治区","青海省","山东省","上海市","山西省","陕西省","四川省","台湾省","天津市","新疆维吾尔自治区","西藏自治区","云南省","浙江省")
names <- c("Anhui","Beijing","Chongqing","Fujian","Gansu","Guangdong","Guangxi","Guizhou","Hainan","Hebei","Heilongjiang","Henan","Hong Kong","Hubei","Hunan","Inner Mongolia","Jiangsu","Jiangxi","Jilin","Liaoning","Ningxia","Qinghai","Shandong","Shanghai","Shanxi","Shaanxi","Sichuan","Taiwan","Tianjin","Xinjiang","Xizang","Yunnan","Zhejiang")
my_dict_EC <- setNames( values,names)
my_dict_CE <- setNames( names,values)

df <- data.frame(code = code,pro_eng = pro_eng)

if (!require("showtext")) {
  install.packages("showtext")}
showtext::showtext_auto()
sf_use_s2(TRUE)

groufile_g <- subset(groufile, Group == "A")
all_province = rownames(groufile_g)
all_code = c()
for (i in seq(from=1, to=length(all_province))){
  temp_code = as.integer(unique(subset(china, Name_Province == my_dict_EC[all_province[i]])$Code_Province))
  all_code = c(all_code,temp_code)
}
df3 <- china %>%
  filter(Code_Province %in% as.character(all_code))
df3 <- df3 %>%
  group_by(Name_Province) %>%
  summarise(geometry = st_union(geometry))

select_province_engA = c()
for (j in df3$Name_Province){
  select_province_engA = c(select_province_engA,my_dict_CE[j])
}

pA <- ggplot(data = df3) +
  geom_sf(fill = province_group_color["A"]) +
  theme_light()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text=element_text(size=12),
        plot.title = element_text(colour = "black", face = "bold", size = 18, vjust = 0.5,hjust = 0.5))+
  theme(legend.position = "none") +
  ggtitle(paste0("Group ","A"))+
  # geom_sf_text(aes(label = select_province_engA),size=1.5,color = "black")
  geom_sf_label(aes(label = select_province_engA),label.size = 0.2,size = 2)
pA

groufile_g <- subset(groufile, Group == "B")
all_province = rownames(groufile_g)
all_code = c()
for (i in seq(from=1, to=length(all_province))){
  temp_code = as.integer(unique(subset(china, Name_Province == my_dict_EC[all_province[i]])$Code_Province))
  all_code = c(all_code,temp_code)
}
df3 <- china %>%
  filter(Code_Province %in% as.character(all_code))
df3 <- df3 %>%
  group_by(Name_Province) %>%
  summarise(geometry = st_union(geometry))
select_province_engB = c()
for (j in df3$Name_Province){
  select_province_engB = c(select_province_engB,my_dict_CE[j])
}

pB <- ggplot(data = df3) +
  geom_sf(fill = province_group_color["B"]) +
  theme_light()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text=element_text(size=12),
        plot.title = element_text(colour = "black", face = "bold", size = 18, vjust = 0.5,hjust = 0.5))+
  theme(legend.position = "none") +
  ggtitle(paste0("Group ","B"))+
  # geom_sf_text(aes(label = select_province_engB),size=1.5,color = "black")
  geom_sf_label(aes(label = select_province_engB),label.size = 0.2,size = 2)
pB

groufile_g <- subset(groufile, Group == "C")
all_province = rownames(groufile_g)
all_code = c()
for (i in seq(from=1, to=length(all_province))){
  temp_code = as.integer(unique(subset(china, Name_Province == my_dict_EC[all_province[i]])$Code_Province))
  all_code = c(all_code,temp_code)
}
df3 <- china %>%
  filter(Code_Province %in% as.character(all_code))
df3 <- df3 %>%
  group_by(Name_Province) %>%
  summarise(geometry = st_union(geometry))
select_province_engC = c()
for (j in df3$Name_Province){
  select_province_engC = c(select_province_engC,my_dict_CE[j])
}
pC <- ggplot(data = df3) +
  geom_sf(fill = province_group_color["C"]) +
  theme_light()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text=element_text(size=12),
        plot.title = element_text(colour = "black", face = "bold", size = 18, vjust = 0.5,hjust = 0.5))+
  theme(legend.position = "none") +
  ggtitle(paste0("Group ","C"))+
  # geom_sf_text(aes(label = select_province_engC),size=1.5,color = "black")
  geom_sf_label(aes(label = select_province_engC),label.size = 0.2,size = 2)
pC

groufile_g <- subset(groufile, Group == "D")
all_province = rownames(groufile_g)
all_code = c()
for (i in seq(from=1, to=length(all_province))){
  temp_code = as.integer(unique(subset(china, Name_Province == my_dict_EC[all_province[i]])$Code_Province))
  all_code = c(all_code,temp_code)
}
df3 <- china %>%
  filter(Code_Province %in% as.character(all_code))
df3 <- df3 %>%
  group_by(Name_Province) %>%
  summarise(geometry = st_union(geometry))
select_province_engD = c()
for (j in df3$Name_Province){
  select_province_engD = c(select_province_engD,my_dict_CE[j])
}
pD <- ggplot(data = df3) +
  geom_sf(fill = province_group_color["D"]) +
  theme_light()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text=element_text(size=12),
        plot.title = element_text(colour = "black", face = "bold", size = 18, vjust = 0.5,hjust = 0.5))+
  theme(legend.position = "none") +
  ggtitle(paste0("Group ","D"))+
  geom_sf_label(aes(label = select_province_engD),label.size = 0.2,size = 2)
pD

groufile_g <- subset(groufile, Group ==  "E")
all_province = rownames(groufile_g)
all_code = c()
for (i in seq(from=1, to=length(all_province))){
  temp_code = as.integer(unique(subset(china, Name_Province == my_dict_EC[all_province[i]])$Code_Province))
  all_code = c(all_code,temp_code)
}

df3 <- china %>%
  filter(Code_Province %in% as.character(all_code))
sf_use_s2(FALSE)
df3 <- df3 %>%
  group_by(Name_Province) %>%
  summarise(geometry = st_union(geometry))

select_province_engE = c()
for (j in df3$Name_Province){
  select_province_engE = c(select_province_engE,my_dict_CE[j])
}

pE <- ggplot(data = df3) +
  geom_sf(fill = province_group_color[ "E"]) +
  theme_light()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text=element_text(size=12),
        plot.title = element_text(colour = "black", face = "bold", size = 18, vjust = 0.5,hjust = 0.5))+
  theme(legend.position = "none") +
  ggtitle(paste0("Group ", "E"))+
  geom_sf_label(aes(label = select_province_engE),label.size = 0.2,size = 2)
pE

groufile_g <- subset(groufile, Group == "F")
all_province = rownames(groufile_g)
all_code = c()
for (i in seq(from=1, to=length(all_province))){
  temp_code = as.integer(unique(subset(china, Name_Province == my_dict_EC[all_province[i]])$Code_Province))
  all_code = c(all_code,temp_code)
}
df3 <- china %>%
  filter(Code_Province %in% as.character(all_code))
df3 <- df3 %>%
  group_by(Name_Province) %>%
  summarise(geometry = st_union(geometry))
select_province_engF = c()
for (j in df3$Name_Province){
  select_province_engF = c(select_province_engF,my_dict_CE[j])
}

pF <- ggplot(data = df3) +
  geom_sf(fill = province_group_color["F"]) +
  theme_light()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text=element_text(size=12),
        plot.title = element_text(colour = "black", face = "bold", size = 18, vjust = 0.5,hjust = 0.5))+
  theme(legend.position = "none") +
  ggtitle(paste0("Group ","F"))+
  geom_sf_label(aes(label = select_province_engF),label.size = 0.2,size = 2)
pF

library(patchwork)
pAll = pA+pB+pC+pD+pE+pF+plot_layout(ncol = 3) & scale_y_continuous(limits = c(18, 55)) & scale_x_continuous(limits = c(73, 135))
pAll
ggsave(filename = paste0(paste0(dirpath,"province_groups"),'.pdf'),plot = pAll,device = "pdf",width = 12, height = 8)
