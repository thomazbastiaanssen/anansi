# Explore possibility of correlations between rust severity and composition of foliar microbiome

library(tidyverse)
library(magrittr)
library(phyloseq)
library(foreach)
library(betareg)
library(vegan)
library(cowplot)
library(ggthemes)
library(ggtext)

source("code/Rfunctions.R")
source("code/colors.R")

### Rust data ###
rust <- loadRust() %>% filter(Treatment!="Control") 

# subset susceptible and resistant genotypes
susceptible <- read.csv("output/tabs/susceptibility.csv", stringsAsFactors = F) 
susceptible$Genotype %<>% factor(.,levels=.[order(susceptible$pctLesion)])

rust.s <-  rust %>%
  filter(Genotype %in% filter(susceptible,Cluster=="susceptible")$Genotype) %>%
  mutate(Genotype=factor(Genotype, levels=levels(susceptible$Genotype)),
         weights=nleaf*length(nleaf)/sum(nleaf)) #proportionality weight for betagreg model

# fit betaregression model of rust variation by genotype and get residuals
resid <- betareg(pctLesion ~ Genotype,weights=weights,data=rust.s) %>% residuals %>%
  data.frame(SampID=rust.s$SampID,
             resid=.)

### Species data ###
phy <- loadPhyloseq() %>% unbias(read.csv("output/tabs/bias.csv"))

# dbRDA to get vectors representing community-level correlation with species arrival treatment, conditioned on host genotype  
# Extract Jensen-Shannon Distance
phy.jsd <- phy %>% distance("jsd") %>% sqrt
# dbRDA
dbRDA.treat <- capscale(phy.jsd~Condition(Genotype)+Treatment, data=data.frame(sample_data(phy)))
# How many significant axes (3)
# anova(dbRDA.treat, by="axis", parallel=15)
# Extract significant axes
dbRDA.scores <- scores(dbRDA.treat, choices=1:3, display="sites") %>% 
  data.frame %>% rownames_to_column("SampID") %>%
  rename("dbRDA axis 1"="CAP1","dbRDA axis 2"="CAP2","dbRDA axis 3"="CAP3") %>%
  mutate(SampID=gsub(".TP1","",SampID,fixed=T))

# extract table of species relative abundance and convert to log odds
otuTab <- phy %>% phy_logit() %>%
  otu_table() %>% data.frame() %>%
  rownames_to_column("SampID") %>%
  mutate(SampID=gsub(".TP1","",SampID,fixed=T))

### Merge data ###
dat <- full_join(resid,otuTab) %>% 
  full_join(dbRDA.scores) %>% drop_na() %>% 
  dplyr::select(-SampID) %>% pivot_longer(-1,"variable") %>%
  mutate(variable=factor(variable,levels=c(colnames(dbRDA.scores)[-1],taxa_names(phy))))

#########################
### Correlation tests ###
#########################

correlation <- foreach(i=unique(dat$variable), .combine=bind_rows) %do% {
  foo <- dat %>% filter(variable==i)
  p <- cor.test(foo$resid,foo$value, method="kendall")$p.value %>% round(2)
  names(p) <- "pval"
  r <- cor.test(foo$resid,foo$value, method="kendall")$estimate %>% round(2)
  c(p,r)
} %>% mutate(variable=unique(dat$variable))

################
### Make fig ###
################

# Markdown formatted facet labels
dat$variable_lab <- ifelse(dat$variable %in% taxa_names(phy),
                           paste0("*",dat$variable,"*"),
                           paste0(dat$variable))
dat$variable_lab %<>% factor(.,levels=unique(.)[c(grep("*",.,fixed=T,invert = T),grep("*",.,fixed=T))])
# Merge correlation test results and format for markdown
cor.dat <- dat %>% select(variable,variable_lab) %>% unique() %>%
  left_join(correlation) %>%
  mutate(stats=paste0("*p*=",pval,"; ","τ=",tau))
#make figure
ggplot(dat, aes(x=value,y=resid))+
  geom_point()+
  facet_wrap(~variable_lab,scales = "free",nrow=4,as.table = T)+
  geom_richtext(data=cor.dat,aes(x=-Inf,y=Inf,label=stats),
                hjust=0,vjust=1,
                fill = NA, label.color = NA)+
  ylim(c(-3.5,4))+
  labs(x="Predictor value",
       y="Residual rust severity")+
  theme_few()+
  panel_border("black",size=0.2)+
  theme(strip.background = element_blank(), 
        strip.text = element_markdown(size=12),
        axis.line = element_blank())
ggsave("output/figs/Fig.S4.jpg",width=7,height=8)
