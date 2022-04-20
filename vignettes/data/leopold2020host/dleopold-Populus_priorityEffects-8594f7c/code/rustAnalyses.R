# Effects of genotype x arrival order interactions on rust severity
# Analyses conducted separately for resistant and susceptible genotypes (likely representing variation in major gene resistance)

library(tidyverse)
library(magrittr)
library(betareg)
library(emmeans)
library(lmtest)
library(ggthemes)
library(patchwork)
library(ggtext)
library(gt)

source("code/Rfunctions.R")
source("code/colors.R")

# load rust data and subset to inoculated plants
rust <- loadRust() %>% filter(Treatment!="Control") 

# subset susceptible and resistant genotypes
susceptible <- read.csv("output/tabs/susceptibility.csv", stringsAsFactors = F) 
susceptible$Genotype %<>% factor(.,levels=.[order(susceptible$pctLesion)])

rust.s <-  rust %>%
  filter(Genotype %in% filter(susceptible,Cluster=="susceptible")$Genotype) %>%
  mutate(Genotype=factor(Genotype, levels=levels(susceptible$Genotype)),
         weights=nleaf*length(nleaf)/sum(nleaf)) #proportionality weight for betagreg model

rust.r <-  rust %>%
  filter(Genotype %in% filter(susceptible,Cluster=="resistant")$Genotype) %>%
  mutate(Genotype=factor(Genotype, levels=levels(susceptible$Genotype)),
         weights=nleaf*length(nleaf)/sum(nleaf)) #proportionality weight for betagreg model

#############################
### Betaregression models ###
#############################

### Lesions ###

# Susceptible genotypes
lesion.s <- betareg(pctLesion ~ Genotype*Treatment, weights=weights, data=rust.s)
lesion.s.noInt <- betareg(pctLesion ~ Genotype+Treatment, weights=weights, data=rust.s)
susceptible.interaction <- lrtest(lesion.s, lesion.s.noInt) #test interation
susceptible.treatment <- lrtest(lesion.s.noInt, .~. -Treatment) #test treatment 
susceptible.genotype <- lrtest(lesion.s.noInt, .~. -Genotype) #test genotype

# Resistant genotypes - precision parameter modeled as function of genotype to account for heteroskedasticity
lesion.r <- betareg(pctLesion ~ Genotype*Treatment|Genotype, weights=weights, data=rust.r)
lesion.r.noInt <- betareg(pctLesion ~ Genotype+Treatment|Genotype, weights=weights, data=rust.r)
resistant.interaction <- lrtest(lesion.r, lesion.r.noInt) #test interation
resistant.treatment<- lrtest(lesion.r.noInt, .~. -Treatment) #test treatment 
resistant.genotype <- lrtest(lesion.r.noInt, .~. -Genotype) #test genotype

results.susceptible <- data.frame(group="Susceptible",
                          Predictor=c("Genotype","Treatment","Genotype:Treatment"),
                          'ΔDf'=abs(c(susceptible.genotype$Df[2],
                                      susceptible.treatment$Df[2],
                                      susceptible.interaction$Df[2])),
                          LR=round(c(susceptible.genotype$Chisq[2],
                               susceptible.treatment$Chisq[2],
                               susceptible.interaction$Chisq[2]),1),
                          P=round(c(susceptible.genotype$`Pr(>Chisq)`[2],
                                     susceptible.treatment$`Pr(>Chisq)`[2],
                                     susceptible.interaction$`Pr(>Chisq)`[2]),3))
results.resistant <- data.frame(group="Resistant",
                          Predictor=c("Genotype","Treatment","Genotype:Treatment"),
                          'ΔDf'=abs(c(resistant.genotype$Df[2],
                                      resistant.treatment$Df[2],
                                      resistant.interaction$Df[2])),
                          LR=round(c(resistant.genotype$Chisq[2],
                                     resistant.treatment$Chisq[2],
                                     resistant.interaction$Chisq[2]),1),
                          P=round(c(resistant.genotype$`Pr(>Chisq)`[2],
                                    resistant.treatment$`Pr(>Chisq)`[2],
                                    resistant.interaction$`Pr(>Chisq)`[2]),3))
# Summarize results
bind_cols(results.susceptible,results.resistant) %>% 
  select(-group,-group1, -Predictor1) %>%
  gt(rowname_col = "Predictor") %>%
  tab_spanner(
    label = "Susceptible",
    columns = vars(ΔDf, LR, P)
  ) %>%
  tab_spanner(
    label = "Resistant",
    columns = vars(ΔDf1, LR1, P1)
  ) %>%
  cols_label(ΔDf1='ΔDf',
              P=md("*P*-value"),
              P1=md("*P*-value"),
              LR=md("LRT-χ<sup>2<sup>"),
              LR1=md("LRT-χ<sup>2<sup>")) %>%
  fmt(
    columns = vars(P),
    fns = function(x) {
      ifelse(x>=0.001,
             as.character(x),
             "< 0.001")
    }
  ) %>%
  cols_align("center") %>%
  gtsave("output/figs/betaReg.lesion.png")
    

### Uridinia ###

# Susceptible genotypes
#uridinia.s <- betareg(pctRust ~ Genotype*Treatment, weights=weights, data=rust.s)
#uridinia.s.noInt <- betareg(pctRust ~ Genotype+Treatment, weights=weights, data=rust.s)
#lrtest(uridinia.s, uridinia.s.noInt) #test interation
#lrtest(uridinia.s.noInt, .~. -Treatment) #test treatment 
#lrtest(uridinia.s.noInt, .~. -Genotype) #test genotype

# Resistant genotypes 
#uridinia.r <- betareg(pctRust ~ Genotype*Treatment|Genotype, weights=weights, data=rust.r)
#uridinia.r.noInt <- betareg(pctRust ~ Genotype+Treatment|Genotype, weights=weights, data=rust.r)
#lrtest(uridinia.r, uridinia.r.noInt) #test interation
#lrtest(uridinia.r.noInt, .~. -Treatment) #test treatment 
#lrtest(uridinia.r.noInt, .~. -Genotype) #test genotype

####################
### Make figures ###
####################

# Susceptible genotypes
(fig.a <- emmeans(lesion.s, ~Treatment|Genotype) %>% data.frame %>%
    mutate(Genotype=factor(Genotype,levels=levels(susceptible$Genotype)),
           upper.CL=asymp.UCL,
           lower.CL=asymp.LCL) %>%
    ggplot(aes(x=Treatment,y=emmean,fill=Treatment)) +
    geom_hline(data=filter(susceptible,Cluster=="susceptible"),
               aes(yintercept=pctLesion),linetype="dotted")+
    geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0)+
    geom_point(shape=21,size=2)+
    geom_point(data=rust.s,
               aes(y=pctLesion,alpha=nleaf),shape=21,size=1,
               show.legend = F)+
    scale_alpha(range=c(0.4,1))+
    labs(y="Rust lesion (%)")+
    scale_fill_manual("Initial colonist:",values=pal.treatment)+
    scale_y_continuous(labels = function(x) x*100)+
   #coord_equal(20)+
   ggtitle("Rust susceptible genotypes")+
    facet_wrap(~Genotype,nrow=1)+
    theme_few()+
    theme(strip.background = element_blank(), 
          axis.text.x = element_blank(),
          strip.text = element_text(size=9),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_markdown(size=9),
          legend.text = element_text(size=9,face="italic")))

# Resistant genotypes
(fig.b <- emmeans(lesion.r.noInt, ~Treatment|Genotype) %>% data.frame %>%
  mutate(Genotype=factor(Genotype,levels=levels(susceptible$Genotype)),
         upper.CL=asymp.UCL,
         lower.CL=asymp.LCL) %>%
  ggplot(aes(x=Treatment,y=emmean,fill=Treatment)) +
  geom_hline(data=filter(susceptible,Cluster=="resistant"),
             aes(yintercept=pctLesion),linetype="dotted")+
  geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0)+
  geom_point(shape=21,size=2)+
  geom_point(data=rust.r,
             aes(y=pctLesion,alpha=nleaf),shape=21,size=1,
             show.legend = F)+
  scale_alpha(range=c(0.4,1))+
  labs(y="Rust lesion (%)")+
  scale_y_continuous(limits=c(0,0.14),
                     breaks=c(0,0.1),
                     labels = function(x) x*100)+
  scale_fill_manual("Initial colonist:",values=pal.treatment)+
    ggtitle("Rust resistant genotypes")+
  facet_wrap(~Genotype,nrow=1)+
  theme_few()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=9),
        legend.text = element_text(size=9,face="italic")))

# Combine into multipanel figure
fig.a / fig.b +
  plot_layout(guides = 'collect',heights = c(4,1.35)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size=10,face="bold"),
        legend.position = 'bottom', legend.title = element_text(size=9),
        plot.title = element_text(size=10,hjust=0.5))
ggsave("output/figs/Fig.4.pdf",width=174,height=110,units="mm")



