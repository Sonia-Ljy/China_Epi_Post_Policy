library(ggplot2)
library(reshape2)
library(dplyr)
library(ggthemes)
library(IDPmisc)
library(tidyr)

out_dir = raw_data <- '/home/soniali/Desktop/02_china_recom_github/'
file <- "0_raw_data/Qualified_china_meta_merged.txt"
setwd(raw_data)

file_name <- print(strsplit(file, "\\.txt")[[1]][1])
data <- read.table(file = paste0(raw_data,file),sep=',',header = T,quote = "")
data <-as.data.frame(data)
head(data)

total_sample = dim(data)[1]
title <- "Source of samples in China"
test1 <- table(data$province)

names(test1) 
as.numeric(test1) 
test2 <- as.data.frame(test1)
test3 <-arrange(test2,desc(Freq))

dat <-as.data.frame(test3)
dat$fraction = dat$Freq / sum(dat$Freq)
dat$ymax = cumsum(dat$fraction)
dat$ymin = c(0, head(dat$ymax, n=-1))

levels<- as.character(dat$Var1)
dat$category <- factor(dat$Var1, levels = levels)
#############
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold",hjust = 0.5)
  )

p1 = ggplot(dat, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect() +
  coord_polar(theta="y") +
  xlim(c(1, 4))
pie<-p1 + scale_fill_manual(values=alpha(c("#3a7f5e","#c4921a","#8e4d93","#257eb2","#f04765","#6b8c42","#a13f7d","#5392a8","#e96d2f","#4c6157","#d78b46","#725ec8","#a8352e","#619b84","#c24d8f","#3f9fca","#8a6e29","#d14779","#53845f","#b54f21","#6d3b8e","#2083b5","#f2734c","#5e7943","#c1749e","#4e7db3","#c96f3a","#587d68","#d05481","#368ba5","#a95835","#7b976c","#e4425d","#B22222","#32CD32","#20B2AA"),0.6))#+scale_fill_brewer("City")
pie_2<-pie+blank_theme +
  theme(axis.text.x=element_blank()) + theme(legend.position=c(.5, .5)) + ggtitle(title) +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  #theme(legend.title = element_text(size=16, face="bold")) +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(size = 8, face = "bold"))

pie_3 <-  pie_2 + geom_text(aes(label = paste0(Freq,"(",round(fraction,2) * 100,"%)"),x = 3.5,y =(ymin + ymax)/ 2),inherit.aes = TRUE,show.legend = FALSE,size = 3.5)
pie_3

ggsave(filename = paste0(paste0(out_dir,"1_epi/China_province_merge"),'.pdf'),plot = pie_3,device = "pdf",width = 8, height = 8)
# 
# 
# out_dir = raw_data <- '/home/soniali/Desktop/02_china_recom/0_raw_data/4_treetime_3945/'
# file <- "Qualified_china_meta-comp-3945_others.csv"
# setwd(raw_data)
# 
# file_name <- print(strsplit(file, "\\.tsv")[[1]][1])
# data <- read.table(file = paste0(raw_data,file),sep=',',header = T,quote = "")
# data <-as.data.frame(data)
# head(data)
# 
# total_sample = dim(data)[1]
# title <- "Source of samples in China"
# test1 <- table(data$province)
# 
# names(test1)
# as.numeric(test1)
# test2 <- as.data.frame(test1)
# test3 <-arrange(test2,desc(Freq))
# 
# dat <-as.data.frame(test3)
# dat$fraction = dat$Freq / sum(dat$Freq)
# dat$ymax = cumsum(dat$fraction)
# dat$ymin = c(0, head(dat$ymax, n=-1))
# 
# levels<- as.character(dat$Var1)
# dat$category <- factor(dat$Var1, levels = levels)
# #############
# blank_theme <- theme_minimal()+
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.border = element_blank(),
#     panel.grid=element_blank(),
#     axis.ticks = element_blank(),
#     plot.title=element_text(size=14, face="bold",hjust = 0.5)
#   )
# 
# 
# nature_colors <- c("#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#E69F00",
#                    "#0072B2", "#F0E442", "#666666", "#000000", "#D8B365",
#                    "#5DA5DA", "#FAA43A", "#B276B2", "#DECF3F", "#F15854",
#                    "#4D4D4D", "#7F7F7F", "#A65628", "#984EA3", "#FF7F00")
# 
# 
# p1 = ggplot(dat, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
#   geom_rect() +
#   coord_polar(theta="y") +
#   xlim(c(1, 4))
# pie<-p1 + scale_fill_manual(values=alpha(nature_colors,0.8))#+scale_fill_brewer("City")
# pie_2<-pie+blank_theme +
#   theme(axis.text.x=element_blank()) + theme(legend.position=c(.5, .5)) + ggtitle(title) +
#   theme(panel.grid=element_blank()) +
#   theme(axis.text=element_blank()) +
#   theme(axis.ticks=element_blank()) +
#   #theme(legend.title = element_text(size=16, face="bold")) +
#   theme(legend.title=element_blank()) +
#   theme(legend.text = element_text(size = 8, face = "bold"))
# 
# pie_3 <-  pie_2 + geom_text(aes(label = paste0(round(fraction,2) ),x = 3.5,y =(ymin + ymax)/ 2),inherit.aes = TRUE,show.legend = FALSE,size = 3.5)
# pie_3
# 
# ggsave(filename = paste0(paste0(out_dir,"China_subsampled+3945"),'.pdf'),plot = pie_3,device = "pdf",width = 8, height = 8)
# 
