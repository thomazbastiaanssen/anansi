# Analysis of the relative advantage of preemptive colonization.

library(tidyverse)
library(magrittr)
library(phyloseq)
library(foreach)
library(doMC)
library(ggthemes)

source("code/Rfunctions.R")
source("code/colors.R")

# register cores for parallel processing
registerDoMC(parallel::detectCores())

# load phyloseq data
(phy <- loadPhyloseq())

# identify focal taxa for plotting 
focalTaxa <- unique(sample_data(phy)$Treatment)

# Extract OTU data, convert to proportions and long format
bias <- read.csv("output/tabs/bias.csv")
df <- phy %>% unbias(.,bias) %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  otu_table %>% data.frame %>%
  dplyr::select(all_of(focalTaxa)) %>%
  bind_cols(sample_data(phy) %>% data.frame %>% 
              dplyr::select(Region,Genotype,Treatment)) %>%
  pivot_longer(all_of(focalTaxa),names_to="Taxa",values_to="proportion") %>%
  mutate(Focal=ifelse(Taxa==Treatment,T,F),
         Region=ifelse(Region=="East","E","W") %>%
           factor(levels=c("W","E")))

# Define a function to calculate the estimated priority effect strength for each species on each host genotype as the log-ratio of the proportional abundance when arriving early vs not arriving early
getPEs <- function(df){
  df %>% group_by(Focal,Region,Taxa,Genotype,Treatment) %>%
    summarize_all(gm_mean) %>% ungroup %>%
    group_by(Region,Genotype,Taxa) %>%
    summarise(PE=log(mean(proportion[Focal])/mean(proportion[!Focal])))
}
meanPEs <- getPEs(df) 
meanPEs$Genotype %<>% factor(.,levels=unique(.))

#####################
### Bootstrap CIs ###
#####################

# Get bootstrap confidence intervals on point estimates of priority effects for each species on each genotype
nboots <- 10000
bootPEs <- foreach(i=1:nboots, .combine=bind_rows) %dopar% {
  df %>% group_by(Region, Taxa, Genotype,Focal,Treatment) %>%
    sample_frac(replace=T) %>%
    getPEs() %>% mutate(bootID=i)
}
bootPEs$Genotype %<>% factor(.,levels=unique(bootPEs$Genotype))

# Get bias corrected and accelerated confidence intervals
bootPE.ci <- bootPEs %>% 
  group_by(Genotype,Taxa) %>%
  summarize(LCI=coxed::bca(PE)[1],
            UCI=coxed::bca(PE)[2]) %>% full_join(meanPEs)

###################
### Region test ###
###################

# First get the mean priority effect for each species on eastern and western genotypes
regionPEs <- meanPEs %>% group_by(Region,Taxa) %>%
  summarize(meanPE=mean(PE),
            tstat=t.test(PE,mu=0,alternative="greater") %>% .$statistic)

# Run bootstrapped t-test using 100 bootstraps at the region level for each of the previously generated genotype-level bootstraps
bootRegion <- foreach(i=1:100,.combine=bind_rows) %dopar% {
  bootPEs %>% 
    group_by(Region,Taxa,bootID) %>%
    sample_frac(replace=T) %>%
    left_join(regionPEs,by=c("Region","Taxa")) %>%
    mutate(center=PE-meanPE) %>%
    summarize(mu.boot=mean(PE),
              tstat.obs=mean(tstat),
              tstat.boot=ifelse(var(center)==0,ifelse(mean(PE)>mean(meanPE),Inf,-Inf),
                                t.test(center,mu=0) %>% .$statistic))}
# calculate bootstrapped p-values (one-tailed test for significant positive priority effects)
bootRegionPvals <- bootRegion %>% drop_na %>%
  mutate(test = tstat.boot+1>tstat.obs+1) %>% 
  summarize(pval=mean(as.numeric(test))) %>%
  mutate(stars=gtools:::stars.pval(pval) %>% gsub(".","+",fixed = T,.)) 
regionSig <- bootRegion %>% summarize(mu.boot.max=max(mu.boot)) %>%
  left_join(bootRegionPvals)

############
### PLOT ###
############

#' Define a function to plot the regional priority effects bootstrap results 
get_inset <- function(df,inset.ymin,inset.ymax){
  ggplot(df, 
         aes(x=Region,y=mu.boot)) +
    geom_violin(scale="width",fill="grey75",color="grey75") + 
    geom_text(data=df %>% dplyr::select(Region,Taxa,stars,mu.boot.max) %>%
                slice(1),
              aes(label=stars,y=mu.boot.max),size=5)+
    geom_hline(yintercept = 0,alpha=0.4,linetype="dotted")+
    coord_cartesian(clip='off')+
    ylim(c(inset.ymin,inset.ymax))+
    theme_few()+
    theme(strip.text = element_blank(),
          axis.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color="black",size=0.2))
}

# create inset plots as a list
insets <- bootRegion %>% 
  left_join(dplyr::select(regionSig,Region,Taxa,stars,mu.boot.max)) %>%
  split(f = .$Taxa) %>%
  purrr::map(~annotation_custom2(
    grob = ggplotGrob(get_inset(.,min(bootRegion$mu.boot),max(bootRegion$mu.boot))), 
    data = data.frame(Taxa=unique(.$Taxa)),
    ymin = 2.4, ymax=Inf, xmin=4, xmax=Inf)
  )

# make full plot of bootstrapped priority effects results
ggplot(bootPE.ci,aes(x=Genotype,y=PE))+
  geom_pointrange(aes(ymin=LCI,ymax=UCI,fill=Genotype,shape=Genotype),
                  size=0.5,stroke=0.35,fatten=5.5)+
  scale_shape_manual("Host genotype", values=c(rep(21,7),rep(23,5)))+
  scale_fill_manual("Host genotype",values=pal.genotype)+
  geom_hline(yintercept = 0,alpha=0.6,linetype="dotted")+
  scale_x_discrete(breaks = levels(bootPEs$Genotype),
                   limits = c(levels(bootPEs$Genotype)[1:7], "skip1",
                              levels(bootPEs$Genotype)[8:12]),
                   expand=expansion(0.075))+
  ylab("Strength of priority effect")+
  labs(tag="Eastern              Western")+
  coord_cartesian(clip='off')+
  ylim(c(-1,3.5))+
  facet_wrap(~Taxa,scales="free_x",nrow=1)+
  guides(fill = guide_legend(override.aes = list(size=0.8)))+
  theme_few()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size=14,face="italic"),
        axis.title.y = element_text(size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.margin=margin(6,6,6,12),
        plot.tag = element_text(angle=90,hjust=0),
        plot.tag.position = c(0.848,0.115))+
  insets
ggsave("output/figs/Fig.3.pdf",width=24,height=10,units="cm")

