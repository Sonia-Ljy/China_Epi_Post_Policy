library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)
library(gridExtra)
library(data.table)
library(scales)
library(lubridate)
library(patchwork)
library(sp)
library(rworldmap)
library(RColorBrewer)
library(pheatmap)

country2lat = function(x)
{  
  country_data = subset(world_data, NAME==x)
  return (as.numeric(country_data[['LAT']]))  
}

country2long = function(x)
{  
  country_data = subset(world_data, NAME==x)
  return (as.numeric(country_data[['LON']]))  
}


export2type = function(x,y)
{  
  if (x==y) {
    return (as.character("Regional"))
  }
  else {
    return (as.character("Global"))
  }
}

dirpath = "/home/soniali/Desktop/02_china_recom_renew/2_treetime/"
setwd(dirpath)
worldmap <- getMap()
world_map_data<-worldmap@data
world_map_data_short<-world_map_data  %>% dplyr::select("NAME","POP_EST")

theUrl <-"provinve_raw_latitude.csv"
world_data <-read.table(file=theUrl, header=TRUE, sep=",")

wave_number = "1"
replicate1<-read.table(file=paste0("mugration_wave",wave_number,".csv"), sep = ',', header = TRUE)
replicate1$replicate <- "1"
head(replicate1)
alpha_all_replicates<-rbind(replicate1)

alpha_all_replicates$date<-date_decimal(alpha_all_replicates$EventTime)
alpha_all_replicates$days<-as.Date(cut(alpha_all_replicates$date,breaks = "day",start.on.monday = FALSE))
alpha_all_replicates$weeks<-as.Date(cut(alpha_all_replicates$date,breaks = "week",start.on.monday = FALSE))
alpha_all_replicates$twoweeks<-as.Date(cut(alpha_all_replicates$date,breaks = "2 week",start.on.monday = FALSE))
alpha_all_replicates$months<-as.Date(cut(alpha_all_replicates$date,breaks = "1 month",start.on.monday = FALSE))

alpha_all_replicates$destination_lat<- lapply(alpha_all_replicates$Destination,country2lat)
alpha_all_replicates$destination_long<- lapply(alpha_all_replicates$Destination,country2long)
alpha_all_replicates$origin_lat<- lapply(alpha_all_replicates$Origin,country2lat)
alpha_all_replicates$origin_long<- lapply(alpha_all_replicates$Origin,country2long)

world1 <- map_data("world", region = c("China"))

world2<- world1 %>% 
  left_join(alpha_all_replicates, by = c("region" = "Destination"))

alpha_all_replicates$Destination<-unlist(alpha_all_replicates$Destination)
alpha_all_replicates<-alpha_all_replicates %>% 
  group_by(Origin) %>% 
  mutate(mean_origin_lat=mean(unlist(origin_lat)), mean_origin_long=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination) %>% 
  mutate(mean_destination_lat=mean(unlist(destination_lat)), mean_destination_long=mean(unlist(destination_long))) %>% 
  ungroup()
alpha_all_replicates<-subset(subset(subset(alpha_all_replicates, Origin!="character(0)"),mean_origin_lat!='numeric(0)'),mean_origin_long!='numeric(0)')
alpha_all_replicates<-alpha_all_replicates %>% 
  group_by(Origin, Destination) %>% 
  mutate(mean_date=mean(days)) %>% 
  ungroup() 

alpha_all_replicates<-alpha_all_replicates %>% 
  group_by(Origin) %>% 
  mutate(mean_origin_lat2=mean(unlist(origin_lat)), mean_origin_long2=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination) %>% 
  mutate(mean_destination_lat2=mean(unlist(destination_lat)), mean_destination_long2=mean(unlist(destination_long))) %>% 
  ungroup()

alpha_origin_continental<-alpha_all_replicates %>%
  group_by(replicate)  %>%
  count(Origin)
colnames(alpha_origin_continental)<-c("replicate","Province","Exports")

# alpha_source_sink <- replace(alpha_source_sink, is.na(alpha_source_sink), 0)

alpha_destination_continental<-alpha_all_replicates %>%
  group_by(replicate)  %>%
  count(Destination)
colnames(alpha_destination_continental)<-c("replicate","Province","Imports")

alpha_source_sink<-left_join(subset(alpha_origin_continental,Province!='character(0)'),subset(alpha_destination_continental,Province!='character(0)'))
alpha_source_sink <- replace(alpha_source_sink, is.na(alpha_source_sink), 0)

alpha_source_sink$net_movements<-alpha_source_sink$Exports-alpha_source_sink$Imports
colnames(alpha_source_sink)<-c("replicate","Origin_Continent","Exports","Imports",'net_movements')


alpha_source_sink$rate = alpha_source_sink$Exports / (alpha_source_sink$Imports +alpha_source_sink$Exports)
alpha_source_sink <- alpha_source_sink %>%
  mutate(name = fct_reorder(Origin_Continent, desc(rate))) #%>%

alpha_source_sink2 <-   mutate(alpha_source_sink, name = fct_reorder(Origin_Continent, desc(rate)))
sorted_data <- alpha_source_sink2[order(alpha_source_sink2$rate, decreasing = TRUE), ]
levels(sorted_data$name)

max(alpha_source_sink$Exports)
max(alpha_source_sink$Imports)

alpha_absolutes <- ggplot(alpha_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('dodgerblue4','red4'),name="",labels=c("Exports","Imports"))+
  geom_bar(stat='summary',aes(x = reorder(Origin_Continent,rate),y=Exports, fill='Exports'))+ #
  geom_bar(stat='summary',aes(x = reorder(Origin_Continent,rate),y=-Imports, fill='Imports'))+
  theme(legend.justification = "left")+
  ylab('No. Viral Movements')+
  ggtitle("The first wave")+
  theme(axis.text=element_text(size=7))+
  scale_y_continuous(limits=c(-100,300),breaks=seq(-100,300,50),labels=abs(seq(-100,300,50)))+
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+ #scale_y_reverse()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))+
  theme(legend.position = "top")

alpha_absolutes
setwd(paste0("/home/soniali/Desktop/02_china_recom_renew/2_treetime/mugration_3wave/wave",wave_number))
ggsave(filename = paste0("source_sink_wave",wave_number,".pdf"),plot = alpha_absolutes,device = "pdf",width = 5, height = 8)

colors <- rainbow(31)
colors
regional_cols<-c("Zhejiang"="#FF0000","Yunnan"="#FF3100","Xinjiang"="#FF6300","Tianjin"="#FF9400","Taiwan"="#FFC500","Sichuan"="#FFF700","Shanxi"="#D6FF00","Shanghai"="#A5FF00","Shandong"="#73FF00","Shaanxi"="#42FF00","Qinghai"="#10FF00","Ningxia"="#00FF21","Liaoning"="#00FF52","Jilin"="#00FF84","Jiangxi"="#00FFB5","Jiangsu"="#00FFE6","Inner Mongolia"="#00E6FF","Hunan"="#00B5FF","Hubei"="#0084FF","Hong Kong"="#0052FF","Henan"="#0021FF","Heilongjiang"="#1000FF","Hebei"="#4200FF","Hainan"="#7300FF","Guizhou"="#A500FF","Guangxi"="#D600FF","Guangdong"="#FF00F7","Gansu"="#FF00C5","Fujian"="#FF0094","Beijing"="#FF0063","Anhui"="#FF0031")

alpha_all_replicates_global_regional<-subset(subset(alpha_all_replicates, Origin!='character(0)'),Destination!='character(0)')
alpha_all_replicates_global_regional$Export_Type<- mapply(export2type,alpha_all_replicates_global_regional$Origin,alpha_all_replicates_global_regional$Destination)
alpha_all_replicates_global_regional_count<-alpha_all_replicates_global_regional %>%
  group_by(Origin,Export_Type)  %>%
  count(Export_Type)  %>%
  ungroup()

colnames(alpha_all_replicates_global_regional_count)<-c("Origin","Export_Type","Counts")

alpha_all_replicates_regional_count<-subset(alpha_all_replicates_global_regional_count,Export_Type=='Regional')
total_regional_counts<-sum(alpha_all_replicates_regional_count$Counts)
alpha_all_replicates_regional_count$ProportionsOfRegional<-alpha_all_replicates_regional_count$Counts/total_regional_counts
alpha_all_replicates_regional_count<-alpha_all_replicates_regional_count%>%select("Origin","ProportionsOfRegional")

mycolors<-brewer.pal(6, "RdBu") 
# ---------------------------------------------
d_alpha_links<-alpha_all_replicates %>%
  group_by(Origin,Destination,mean_origin_lat,mean_origin_long,mean_destination_lat,mean_destination_long, mean_date)  %>%
  count()
d_alpha_links2 <- d_alpha_links[order(d_alpha_links$n, decreasing = TRUE), ]
threshold = "t20"
d_alpha_links <- head(d_alpha_links2, n = 20)
# d_alpha_links<-subset(subset(subset(d_alpha_links,n>=threshold),Origin!='character(0)'),Destination!='character(0)')
write.table(d_alpha_links,file = "d_alpha_links_wave1.csv",append = FALSE,sep = ",")

d_alpha_origins<-alpha_all_replicates %>%
  group_by(Origin,mean_origin_lat,mean_origin_long)  %>%
  count()
# d_alpha_origins <-subset(subset(subset(d_alpha_origins,n>=threshold)))
d_alpha_links

lab_dates <- pretty(d_alpha_links$mean_date)

alpha_map <- ggplot() +
  theme_void()+
  geom_map(data=world1,map=world1, aes(long, lat,map_id=region),fill='snow3',size=2)+ #color="white",
  geom_curve(data = d_alpha_links,
             aes(x = as.double(mean_origin_long),
                 y = as.double(mean_origin_lat),
                 xend = as.double(mean_destination_long),
                 yend = as.double(mean_destination_lat), colour=as.Date(mean_date),size=(n/10)),alpha=0.8,
             arrow = arrow(length = unit(0.03, "npc")))+
  geom_point(data = d_alpha_links,
             aes(x = as.double(mean_origin_long), y = as.double(mean_origin_lat), size=12),
             color='#9C9A99', shape=21)+
  scale_color_gradientn(colours =mycolors,breaks = lab_dates, 
                        labels = lab_dates, name='Inferred Dispersal Dates (Mean)')+
  scale_fill_manual(values=regional_cols, name='Regions of Origin')+
  theme(legend.position = 'bottom')+
  theme(legend.direction = 'horizontal')+
  coord_fixed()+
  guides(colour = guide_colourbar(barwidth = 20, barheight = 0.5,title.position = 'top',title.hjust = 1,ticks.colour = "white",
                                  ticks.linewidth = 2), size='none')+
  scale_y_continuous(limits=c(0,60))+
  scale_x_continuous(limits=c(70,140))+
  theme(plot.title = element_text(size = rel(1.2), hjust=0.5,lineheight = 0.5,family = "Arial", colour = "black"))

alpha_map
# ggsave(filename = paste0("wave",wave_number,"_",threshold,"_lines.pdf"),plot = alpha_map,device = "pdf",width = 6, height = 6)
