library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(forcats)
library(dplyr)
library(reshape2)
library(patchwork) 

setwd("/home/soniali/Desktop/02_china_recom_renew/1_epi/world_china_compare_GISAID/")
df <- read.table("Qualified_world_china_meta_merged_202212_202310.txt",sep = ",",header = TRUE)

unique(df$merged_lineage)
df <- subset(df, df$merged_lineage != "Others")
new_df <- subset(df, df$country != "South Africa")
df_time = as.Date(new_df$date,"%Y-%m-%d")
df_time2 <- as.numeric(df_time)
new_df$date_num <- df_time2
new_df$merged_lineage <- factor(
  new_df$merged_lineage,
 )
head(new_df)
unique(new_df$merged_lineage)

df_date_country  <- data.frame()
colors_list = c("#82ADCF","#4F72B1","#F0B26D","#F4D586","#B7D7E9","#5D9488","#9372A9","#C5432E") 
lineages = c('BA.5.2.48*', 'BF.7.14*', 'BA.2.75*', 'XBB.1.5*', 'XBB.1.9.1*', 'XBB.1.22.1*', 'XBB.1.16*', 'EG.5.1*')
new_df$merged_lineage <- factor(new_df$merged_lineage,level = lineages) 

p1 <- ggplot(new_df, aes(x=date_num, y=country, color=merged_lineage, point_color=merged_lineage, fill=merged_lineage)) +
  geom_density_ridges(
    jittered_points=FALSE, scale = .95, rel_min_height = .01,alpha = 0.6,
    position = position_points_jitter(height = 0)
  ) + facet_wrap(~merged_lineage,nrow = 2)+
  scale_y_discrete(expand = c(.01, 0)) +
  scale_x_continuous(expand = c(0, 0), name = "Collection Date") +
  scale_fill_manual(values = colors_list) + 
  scale_color_manual(values = colors_list, guide = "none") +
  scale_discrete_manual("point_color", values = colors_list, guide = "none") +
  theme_ridges(center = TRUE)+
  theme(legend.position = "none")
p1

ggsave(filename = paste0("plot_top8_lineages",".pdf"),plot = p1,device = "pdf",width = 12, height = 5)

new_df_china <- subset(new_df, new_df$country == "China")
write.csv(new_df_china, "plot_top8_lineages_china.csv", row.names = FALSE)

p2 <- ggplot(new_df, aes(x=date_num, y=country, color=merged_lineage, point_color=merged_lineage, fill=merged_lineage)) +
  geom_density_ridges(
    jittered_points=FALSE, scale = .95, rel_min_height = .01,alpha = 0.3,
    position = position_points_jitter(height = 0)
  ) + 
  scale_y_discrete(expand = c(.01, 0)) +
  scale_x_continuous(expand = c(0, 0), name = "Collection Date") +
  scale_fill_manual(values = colors_list) +
  scale_color_manual(values = colors_list, guide = "none") +
  scale_discrete_manual("point_color", values = colors_list, guide = "none") +
  theme_ridges(center = TRUE)+
  theme(legend.position = "bottom")
p2

ggsave(filename = paste0("plot_top8_lineages_merge",".pdf"),plot = p2,device = "pdf",width = 8, height = 10)

as.Date(c(19400,19500,19600),origin = "1970-01-01")
# [1] "2023-02-12" "2023-05-23" "2023-08-31"
