# This script will estimate taxon specific bias following [McLaren et al (2019)](https://doi.org/10.7554/eLife.46923.001).  
# Function from the R-package associated with the manuscript [metacal](https://github.com/elifesciences-publications/metacal) are used.

library(tidyverse)
library(magrittr)
library(phyloseq)
library(ggthemes)
library(ggbeeswarm)
library(metacal)
library(patchwork)

# Read in phyloseq object and subset experimental samples
# Removing singletons and doubletons leaves only the expected species
phy <- readRDS("output/compiled/phy.rds") %>% 
  prune_samples(grepl("Mock",sample_names(.)),.) %>% 
  filter_taxa(., function(x) {sum(x>0) > 2}, TRUE)

# Extract sequence counts and convert to long format
dat <- phy %>% otu_table %>% data.frame %>%
  rownames_to_column("Sample") %>%
  gather("Taxon","Observed",-1)

# Read in known concentrations of the mock communities.
# Listed values are the dilution factors used when making 
# the mock communities and need to be inverted to represent 
# known relative abundance
mock <- read.csv("data/MockCommunities.csv",row.names = 1) %>% 
  dplyr::select(-Sym4) %>%
  rownames_to_column("Sample") %>%
  gather("Taxon","Actual",-1) %>%
  mutate(Actual=1/Actual) 

# Merge known and observed abundances
dat %<>% full_join(mock,by=c("Sample","Taxon"))

# Preliminary analysis showed sample 'Mock.5' to be an outlier, so we will remove it 
# (though this has little effect on resulting bias estimates)
dat %<>% filter(Sample!="Mock.5")

# Estimate error 
dat %<>% mutate(Error=Observed/Actual)
error_mat <- build_matrix(dat, Sample, Taxon, Error)

# Estimate bias
(bias <- center(error_mat, enframe = TRUE) %>%
  dplyr::rename(Bhat = Center))

# Use bootstrapping to estimate error associated with bias estimate
bootreps <- bootrep_center(error_mat) %>%
  dplyr::rename(Bhat = Center)
bootreps.summary <- bootreps %>%
  group_by(Taxon) %>%
  summarize(Gm_mean = gm_mean(Bhat), Gm_se = gm_sd(Bhat))
(bias0 <- left_join(bias, bootreps.summary, by = "Taxon"))

# Write bias estimates to csv file
write.csv(bias0,"output/tabs/bias.csv",row.names = F)

# Figure of bias estimates
fig.a <- bias0 %>%
  dplyr::rename(Error=Bhat) %>%
  mutate(Taxon=factor(Taxon,levels=bias0$Taxon[order(bias0$Bhat)])) %>%
  ggplot(aes(x=Taxon,y=Error,fill=Taxon))+
  geom_hline(yintercept = 1, color = "grey") +
  geom_pointrange(aes(ymin = Error / Gm_se^2, ymax = Error * Gm_se^2),
                  shape=21) +
  geom_point(shape=21,size=2.5)+
  scale_fill_tableau("Color Blind")+
  scale_y_log10() +
  labs(x="",y="Bias estimate") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size=12,face="italic"),
        axis.title = element_text(size=12))

# Figure showing effect of bias correction
fig.b <- dat %>% 
  left_join(bias0)  %>%
  mutate_by("Sample", Corrected=logit(close_elts(Observed/Bhat)),
            Actual=logit(close_elts(Actual)),
            Observed=logit(close_elts(Observed))) %>%
  gather("Results","Value",c(3,9)) %>%
  mutate(Results=ifelse(Results=="Observed",
                        "Raw data",
                        "Bias-corrected data") %>%
           factor(.,levels=c("Raw data","Bias-corrected data"))) %>%
  ggplot(aes(x=Actual,y=Value,fill=Taxon)) +
  geom_quasirandom(width=0.2,shape=21,size=2.5)+
  scale_fill_tableau("Color Blind")+
  labs(x="logit (known proportions)",
       y="logit (observed proportions)")+
  facet_wrap(Results~.,ncol=1)+
  theme_few()+
  theme(legend.position = "none",
        axis.title = element_text(size=10),
        strip.text = element_text(size=12,hjus=0))

fig.a + fig.b +
  plot_layout(ncol=2,widths=c(1,1.2)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size=14,face="bold"))
ggsave("output/figs/Fig.S2.jpg",width=20,height=10,units="cm")



