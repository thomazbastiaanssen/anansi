# make map showing origins of P. trichocarpa genotypes

library(tidyverse)
library(magrittr)
library(ggmap)
library(ggrepel)
library(cowplot)
library(magick)

points <- read.csv("data/GenotypeOrigins.csv",header=T,as.is=T) %>%
  transmute(Latitude=Latitude,
            Longitude=Longitude,
            Label=paste0("textstyle(`",Genotype,"`)"),
            Region=Region)
fungi <- read.csv("data/FungiOrigins.csv",as.is=T) %>%
  transmute(Latitude=Latitude,
            Longitude=Longitude,
            Label=paste0("italic(",Genus,")"),
            Region=NA)
merged <-  bind_rows(points,fungi)


ptri <- png::readPNG("data/500px-Populus_trichocarpa_range_map.svg.png")

map <- make_bbox(points$Longitude,points$Latitud,f=0.35) %>%
  get_map(zoom=7, maptype = "watercolor",source="osm",force=F) %>%
  ggmap()+
  geom_label_repel(data=merged,aes(x=Longitude,y=Latitude,label=Label),
                   point.padding = 0.3,box.padding = 0.5, parse=T,
                   alpha=1, segment.color="black",
                   fill=c(rep("white",nrow(points)),rep("grey85",nrow(fungi))),
                   color=c(rep("black",nrow(points)),rep("black",nrow(fungi))))+
  geom_point(data=points,aes(x=Longitude,y=Latitude,fill=Region,shape=Region),size=3)+
  geom_point(data=fungi,aes(x=Longitude,y=Latitude),shape=21,size=2,fill="white")+
  scale_shape_manual("Ecotype",values=c(23,21))+
  scale_fill_manual("Ecotype",values = c("#fb832d","#016392"))+

  labs(x="Longitude",y="Latitude")

ggdraw()+
  draw_plot(map)+
  draw_image("data/500px-Populus_trichocarpa_range_map.svg.png",scale=1,width=0.32,x=0.52,y=0.24)
ggsave("output/figs/Fig.S1.jpg",width=20,height=16,units="cm")
  
