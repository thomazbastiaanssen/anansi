# Make multipanel Figure 1.
# Fig 1a = heatmap of standardized relative abundance
# Fig 1b = heatmap summarizing univariate model results from Join-species models
# Fig 1c = db-RDA ordination showing initial colonist treatment effects, conditioned on host genotype
# Fig 1d = db-RDA ordination showing host genotype / ecotype effects, conditioned on initial colonist

# Load packages
library(tidyverse)
library(magrittr)
library(mvabund)
library(phyloseq)
library(foreach)
library(MASS)
library(vegan)
library(ggtext)
library(ggnewscale)
library(ggvegan)
library(cowplot)
library(scico)
library(patchwork)

source("code/Rfunctions.R")
source("code/colors.R")

# load phyloseq data
(phy <- loadPhyloseq())
#Taxon-specific sequencing bias estimates
bias <- read.csv("output/tabs/bias.csv")

# Extract OTU data and convert to proportional and tidy format
df.tidy <- phy %>% unbias(.,bias) %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  otu_table %>% data.frame %>% 
  bind_cols(sample_data(phy) %>% data.frame %>% 
              dplyr::select(Region,Genotype,Treatment)) %>%
  pivot_longer(all_of(taxa_names(phy)),names_to="Taxa",values_to="abundance") %>%
  mutate(ID=paste(Treatment,Genotype),
         Region=ifelse(Region=="East","E","W") %>%
           factor(levels=c("W","E")),
         Focal=ifelse(Taxa==Treatment,T,F)) 

#################################################
### Fig 1c - dbRDA of initial colonist effect ###
#################################################

# Summarize and spread to wide format
df.1c <- df.tidy %>% 
  group_by(Region,Genotype,Treatment,Taxa) %>%
  summarize(abundance=median(abundance)) %>%
  ungroup %>%
  pivot_wider(names_from = Taxa, values_from=abundance)

# Calculate Jensen-Shannon distance on standardized matrix
JSDdist <- df.1c %>% select_if(is.numeric) %>%
  wisconsin %>% as.matrix %>% 
  philentropy::JSD(est.prob = "empirical") %>% sqrt

# fit dbRDA constrained by treatment, conditioned on genotype
(dbrda.treatment <- capscale(JSDdist~Condition(Genotype)+Treatment, data=df.1c))
# get variance explained by the ordination axes
dbrda.treatment.r2 <- (eigenvals(dbrda.treatment)/sum(eigenvals(dbrda.treatment))) %>% multiply_by(100) %>% round(1)
# plot ordination
(ordination.treatment <- fortify(dbrda.treatment,display="wa") %>%
    filter(Score=="sites") %>% 
    bind_cols(df.1c) %>%
    group_by(Treatment) %>%
    summarize(n=n(),x=mean(CAP1),y=mean(CAP2),sd1=sd(CAP1),sd2=sd(CAP2)) %>%
    ggplot(aes(x=x,y=y,fill=Treatment)) +
    geom_errorbar(aes(ymin=y-sd2,ymax=y+sd2),width=0)+
    geom_errorbarh(aes(xmin=x-sd1,xmax=x+sd1),height=0)+
    geom_point(size=3.5,shape=21)+
    scale_fill_manual("Initial\ncolonist",values=pal.treatment)+
    guides(fill=guide_legend(nrow=3))+
    coord_fixed(eigenvals(dbrda.treatment)[2]/eigenvals(dbrda.treatment)[1],
                clip="off")+
    labs(x=paste0("dbRDA1 [",dbrda.treatment.r2[1],"%]"),
         y=paste0("dbRDA2 [",dbrda.treatment.r2[2],"%]"))+
    ggthemes::theme_few()+
    theme(legend.justification = "left",
          legend.position = "none",
          axis.title = element_text(size=12),
          axis.text = element_text(size=10),
          panel.border = element_rect(fill=NA, colour = "black", size=0.6)))

########################################
### Fig 1d - dbRDA of ecotype effect ###
########################################

# fit dbRDA constrained by genotype, conditioned on treatment
(dbrda.genotype <- capscale(JSDdist~Condition(Treatment)+Genotype, data=df.1c))
# get variance explained by the ordination axes
dbrda.genotype.r2 <- (eigenvals(dbrda.genotype)/sum(eigenvals(dbrda.genotype))) %>% multiply_by(100) %>% round(1)
#' plot ordination
(ordination.genotype <- fortify(dbrda.genotype,display="wa") %>%
    filter(Score=="sites") %>% 
    bind_cols(df.1c) %>%
    group_by(Genotype) %>%
    summarize(n=n(),x=mean(CAP1),y=mean(CAP2),sd1=sd(CAP1),sd2=sd(CAP2)) %>%
    left_join(sample_data(phy) %>%
                data.frame %>%
                dplyr::select(Genotype,Region) %>%
                unique) %>%
    ggplot(aes(x=x,y=y,fill=Region)) +
    geom_errorbar(aes(ymin=y-sd2,ymax=y+sd2),width=0)+
    geom_errorbarh(aes(xmin=x-sd1,xmax=x+sd1),height=0)+
    geom_point(size=3.5,shape=21)+
    scale_fill_manual("Host\necotype", values=pal.region)+
    coord_fixed(eigenvals(dbrda.genotype)[2]/eigenvals(dbrda.genotype)[1])+
    labs(x=paste0("dbRDA1 [",dbrda.genotype.r2[1],"%]"),
         y=paste0("dbRDA2 [",dbrda.genotype.r2[2],"%]"))+
    guides(fill=guide_legend(nrow=2))+
    theme_few()+
    theme(legend.justification = "left",
          legend.position = "none",
          axis.title = element_text(size=12),
          axis.text = element_text(size=10),
          panel.border = element_rect(fill=NA, colour = "black", size=1)))

###########################################
### FIG 1A - Relative abundance heatmap ###
###########################################

# Summarize replicates and Standardize by taxa
df.1a <- df.tidy %>% 
  group_by(ID,Region,Genotype,Treatment,Taxa) %>%
  summarize(abundance=median(abundance)) %>%
  ungroup() %>%
  mutate_by(Taxa,
            normed=scale(abundance),
            xlabs=ifelse(grepl("East",ID),
                         paste0("E<br>",stringr::str_sub(ID,-1)),
                         paste0("W<br>",stringr::str_sub(ID,-1))))

# Set order of taxa based on relative abundance for platting
#taxa.order <- tapply(df.1a$abundance,df.1a$Taxa,median) %>% sort(decreasing = F) %>% names
taxa.order <- c("Trichoderma", "Penicillium", "Epicoccum", "Fusarium",
                "Dioszegia", "Cladosporium", "Aureobasidium", "Alternaria")
df.1a$Taxa %<>% factor(levels=taxa.order)

# Set order of genotypes for plotting 
gt_order <- fortify(dbrda.genotype,display="wa") %>%
  filter(Score=="sites") %>% 
  bind_cols(df.1c) %>%
  group_by(Genotype) %>%
  summarize(n=n(),x=mean(CAP1),y=mean(CAP2),sd1=sd(CAP1),sd2=sd(CAP2)) %>%
  arrange(x) %$% Genotype %>% paste(rep(unique(df.1a$Treatment),each=length(.)),.)
df.1a$ID %<>% factor(levels=gt_order)  

# Coordinates of focal taxa annotations
focal.rect <- data.frame(Treatment=unique(df.1a$Treatment),
                             xmin=0.5,xmax=12.5,
                             ymin=c(7.5,6.5,5.5,4.5,3.5),
                             ymax=c(8.5,7.5,6.5,5.5,4.5))

# Hack to add legend for panel d to x-axis without changing y-axis scale
ggplot(data.frame(lab=c("western","eastern")),
       aes(x=c(1,1.13),y=c(1,1),fill=lab)) +
  lims(x=c(0.985,1.22))+
  geom_point(shape=21,size=3) +
  scale_fill_manual(values=pal.region)+
  geom_text(aes(label=lab),hjust=-0.2)+
  theme_nothing()
ggsave("output/figs/ecotypeHack.jpg",width=4.5,height=0.4,units = "cm")

# Plot Fig 1a heatmap
(fig1a <- ggplot(df.1a,aes(x=ID,y=Taxa))+
    geom_tile(color="grey50",size=0.2,aes(fill=normed)) +
    scale_fill_scico("Standardized\nrelative abundance   ",
                     palette = "vik", rescaler = mid_rescaler(),
                     limits=c(-2.5,4.5),breaks=c(-2,0,2,4)) +
    geom_rect(data=focal.rect,fill="NA",inherit.aes = F,color="black",size=0.65,
              aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))+
    geom_richtext(data=filter(df.1a,Taxa=="Alternaria"),
                  aes(label=xlabs,y=0.4), family = "mono",
                  fill = ifelse(filter(df.1a,Taxa=="Alternaria")$Region=="E",
                                pal.region[1],pal.region[2]),
                  label.padding = grid::unit(c(3,1.5,1,1.5), "pt"),label.color = "grey50",
                  label.r = unit(0, "lines"), size = 3, vjust=1) +
    geom_richtext(data=filter(df.1a,Taxa=="Alternaria" & Genotype=="East-1"),
              aes(label=paste0("*",Treatment,"*"),y=8.95,x=2),
              hjust=0,size=4.25,
              fill = NA, label.color = NA,
              label.padding = grid::unit(rep(0, 4), "pt"))+
    geom_point(data=filter(df.1a,Taxa=="Alternaria" & Genotype=="East-1"),
              aes(y=9.1,x=1),
              size=3.25,shape=19,show.legend=F, color=pal.treatment)+
    geom_point(data=filter(df.1a,Taxa=="Alternaria" & Genotype=="East-1"),
               aes(y=9.1,x=1),
               size=3.25,shape=21,show.legend=F)+
    coord_cartesian(clip="off") +
    guides(fill=guide_colourbar(ticks=F,
                                direction="horizontal",
                                title.position="left"))+
    facet_wrap(~Treatment,scales="free_x",nrow=1) +
    labs(title="Initial colonist:\n",
         x="Host ecotype: <img src='output/figs/ecotypeHack.jpg' width='150' />")+
    theme_minimal_grid()+
    theme(panel.grid.major = element_blank(),
          legend.position = "bottom",
          legend.margin=margin(0,0,0,0),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          axis.title.x = element_markdown(size=14,margin=margin(26,0,0,0)),
          plot.title = element_text(size=14,margin=margin(0,0,0,0), 
                                    hjust=0.5, face="plain"),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=12,face="italic"),
          axis.text.x =  element_blank(),
          strip.text = element_blank(),
          legend.text = element_text(size=10),
          legend.title = element_text(size=12),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt")))

##########################################################
### Figure 1b - heatmap of species-level model results ###
##########################################################

# Create offset variable for glms (effort)
effort <- (sample_sums(unbias(phy, bias)) %*% t(bias$Bhat[match(taxa_names(phy),bias$Taxon)])) %>% log 
colnames(effort) <- taxa_names(phy)
effort.long <- effort %>% data.frame %>%
  mutate(SampID=sample_names(phy)) %>%
  pivot_longer(1:8,"Taxon",values_to = "effort")

# Fit individual nb.glm models and extract pseudo R2s
uni.r2 <- otu_table(phy) %>% data.frame %>%
  rownames_to_column("SampID") %>%
  mutate(Genotype=sample_data(phy)$Genotype,
         Treatment=sample_data(phy)$Treatment,
         Region=sample_data(phy)$Region) %>%
  pivot_longer(2:9,"Taxon",values_to="abund") %>%
  left_join(effort.long) %>%
  split(.$Taxon) %>%
  map_df( function(x){
    full.genotype <- glm.nb(abund~Treatment*Genotype+offset(effort),data=x)
    full.region <- glm.nb(abund~Treatment*Region+offset(effort),data=x)
    list(Treatment=rsq::rsq.partial(full.genotype,update(full.genotype, .~. -Treatment/Genotype),type='n',adj = T)$partial.rsq,
         Genotype=rsq::rsq.partial(full.genotype,update(full.genotype, .~. -Genotype/Treatment),type='n',adj = T)$partial.rsq,
         'Genotype:Treatment'=rsq::rsq.partial(full.genotype,update(full.genotype, .~. -Genotype:Treatment),type='n',adj = T)$partial.rsq,
         Region=rsq::rsq.partial(full.region,update(full.region, .~. -Region/Treatment),type='n',adj = T)$partial.rsq,
         'Region:Treatment'=rsq::rsq.partial(full.region,update(full.region, .~. -Region:Treatment),type='n',adj = T)$partial.rsq) %>% return}, .id="Taxon") %>%
  pivot_longer(., 2:6,"predictor",values_to = "partialR2")

# Get univariate p-values from mvabund model 
mv.genotype <- readRDS("output/rds/mv.genotype.rds")
mv.region <- readRDS("output/rds/mv.region.rds")
uni.pvals <- mv.genotype$uni.p %>% data.frame %>% .[-1,] %>%
  rbind(mv.region$uni.p %>% data.frame %>% .[c(-1,-3),]) %>%
  rownames_to_column("predictor") %>%
  pivot_longer(2:9,"Taxon",values_to = "pval") %>%
  mutate(stars=gtools::stars.pval(pval) %>% gsub(".","+",.,fixed=T))

# Merge univariate data for plotting
uni.dat <- full_join(uni.r2,uni.pvals) 
uni.dat$predictor %<>% factor(levels=c("Region","Genotype","Treatment","Genotype:Treatment","Region:Treatment"))
uni.dat$Taxon %<>% factor(levels=taxa.order)

uni.dat %<>% mutate(xlabs=gsub("Region","Ecotype",predictor) %>% 
                      gsub("Genotype:Treatment","Treatment x<br>Genotype",.) %>% 
                      gsub("Ecotype:Treatment","Treatment x<br>Ecotype",.))

# make plot
(fig1b <- ggplot(uni.dat,aes(x=predictor,y=Taxon,fill=partialR2))+
  geom_tile(color="grey30",size=0.2)+
    scale_fill_scico("Partial<br>*pseudo*-*r*^2",
                     palette="vik",begin=0.5, end=0.9,
                     breaks=c(0,0.25,0.5),labels=c("0.0","0.25","0.5"))+
    geom_text(aes(label=stars),size=5)+
    geom_richtext(data=filter(uni.dat,Taxon=="Alternaria"),
                  aes(label=xlabs,y=0.4), size=3.2,
                  fill = NA, angle=90, hjust=1,vjust=0.5,
                  label.color = NA, lineheight=0)+
    geom_point(data=filter(uni.dat,Taxon=="Alternaria"),
               aes(y=9.1,x=1),
               size=3.25,shape=21,show.legend=F, color="white", fill="white")+
    guides(fill=guide_colourbar(ticks=F))+
    coord_cartesian(clip="off") +
    guides(fill=guide_colourbar(ticks=F,
                                direction="horizontal",
                                title.position="left"))+
    theme_minimal_grid()+
    theme(panel.grid.major = element_blank(),
          legend.position = "bottom",
          legend.margin=margin(0,50,0,10),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x =  element_blank(),
          legend.text = element_text(size=10),
          legend.title = element_markdown(size=12),
          plot.margin = unit(c(0,5.5,0,5.5), "pt")))

########################
### Compile Figure 1 ###
########################

((fig1a + fig1b + plot_layout(widths = c(5,1), ncol=2)) / 
(ordination.treatment + ordination.genotype + plot_layout(ncol=2))) + 
  plot_layout(heights=c(5,4)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size=14,face="bold"))
ggsave("output/figs/Fig.2.pdf",units="cm",width=28,height=18)

#Remove tmp file
unlink("output/figs/ecotypeHack.jpg")
