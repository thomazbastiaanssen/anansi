# Analyze uninoculated control plants to determine baseline rust susceptibility.

library(tidyverse)
library(magrittr)
library(ggthemes)
library(betareg)

source("code/Rfunctions.R")
source("code/colors.R")

# load rust data and subset to uninoculated plants
rust.negs <- loadRust() %>% filter(Treatment=="Control") 

# Estimate weighted means and weighted standard errors
rust.negs.means <- rust.negs %>% 
  group_by(Genotype,Region) %>%
  summarise(pctLesion.se=diagis::weighted_se(pctLesion,leafArea),
            pctLesion=weighted.mean(pctLesion,leafArea),
            pctRust.se=diagis::weighted_se(pctRust,leafArea),
            pctRust=weighted.mean(pctRust,leafArea))

# Plot rust severity in negative control plants with point size proportional to leaf area. Vertical lines indicate the weighted means. This figure confirms our expectation of greater rust resistance in western genotypes. However, one western genotype is very susceptible
rust.negs %>% 
  ggplot(aes(x=Genotype,y=pctLesion)) +
  geom_point(aes(size=leafArea),shape=19,alpha=0.3)+
  geom_point(data=rust.negs.means,shape="|",size=5)+
  facet_wrap(~Region,scales="free_y",ncol=1) +
  labs(x="",y="Proportion chlorotic lesion") +
  coord_flip()+
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid.minor.y=element_line(color="black",size=1),
        legend.position = "none")

# Similar results for uredinia. 
rust.negs %>% 
  ggplot(aes(x=Genotype,y=pctRust)) +
  geom_point(aes(size=leafArea),shape=19,alpha=0.3)+
  geom_point(data=rust.negs.means,shape="|",size=5)+
  facet_wrap(~Region,scales="free_y",ncol=1) +
  labs(x="",y="Proportion uredinia") +
  coord_flip()+
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid.minor.y=element_line(color="black",size=1),
        legend.position = "none")

# Test the overall trend in rust severity among populations east and west of Cascades.  
# Test for difference in chlorotic lesions
betareg(pctLesion  ~ Region, data=rust.negs.means) %>% summary
# Test for difference in uredinia
betareg(pctRust  ~ Region, data=rust.negs.means) %>% summary

# Use k-means clustering to partition genotypes into 2 groups (resistant and susceptible) based on both % lesion and % uredinia.  
# Log transformation is applied prior to clustering to emphasize lower values. This is consistent with emphasizing a split on with/without major gene resistance, as opposed to identifying clusters of varying quantative resitance.
rust.negs.means$Cluster <- kmeans(log(rust.negs.means[,3:4]),2)$cluster
susceptible <- rust.negs.means$Cluster[order(rust.negs.means$pctLesion, decreasing = T)[1]]
rust.negs.means$Cluster <-  ifelse(rust.negs.means$Cluster==susceptible,"susceptible","resistant")

# Scatter plot of clustering results confirms that one Western genotype can be considered susceptible.
ggplot(rust.negs.means,aes(x=pctLesion,y=pctRust,fill=Region,shape=Cluster))+
  geom_errorbar(aes(ymin=pctRust-pctRust.se,ymax=pctRust+pctRust.se),width=0)+
  geom_errorbarh(aes(xmin=pctLesion-pctLesion.se,xmax=pctLesion+pctLesion.se),height=0)+
  geom_point(size=4)+
  scale_shape_manual(values = c(21,23),name="K-means\nclustering")+
  scale_fill_manual("Ecotype",values = pal.region)+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="Proportion chlorotic lesion", y="Proportion uredinia")+
  guides(fill=guide_legend(order=1,override.aes=list(shape=21)))+
  theme_few()+
  theme(legend.position = c(0.05,0.98),
        legend.justification = c(0, 1))
ggsave("output/figs/Fig.S3.jpg",width=6,height=4)

# Save susceptibility data
write.csv(rust.negs.means,"output/tabs/susceptibility.csv")



