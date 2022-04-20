# Custom R functions

#__________________________________________________#
### load and subset phyloseq object for analyses ###
loadPhyloseq <- function(){
  readRDS("output/compiled/phy.rds") %>% 
    #Remove control and mock community samps
    subset_samples(Samp_type=="Experiment") %>%  
    #Remove uninoculated plants
    subset_samples(Treatment!="Negative") %>%
    #Subset to time point 1 (time point 2 was durring rust sampling and community data was overwhelmed with rust reads)
    subset_samples(Timepoint==1) %>%
    #Remove non-target taxa
    prune_taxa(taxa_names(.)!="Melampsora" & !grepl("OTU",taxa_names(.)),.) %>%
    #Remove two outlier samples that cause problems for model fitting
    #prune_samples(!(grepl("G6.T3.R4.TP1",sample_names(.))),.) #%>%
    #This sample has an order of magnitude greater Trichoderma which causes instability in the mvabund model
    prune_samples(!(grepl("G4.T2.R5.TP1",sample_names(.))),.)
}

#_______________________________________________________________#
### Load rust severity data and aggregate to plant-level data ###
loadRust <- function(){
  dat <-   read.csv("data/rust_measurements.csv",as.is=T) %>%
        group_by(SampID,Genotype,Region,Treatment,Replicate) %>% 
        summarise(pctRust=sum(Rust_cm2)/sum(Leaf_cm2),
                  pctLesion=sum(Lesion_cm2)/sum(Leaf_cm2),
                  leafArea=sum(Leaf_cm2),
                  nleaf=n()) %>% ungroup
  # Transform uridinia data to deal with zeros, following Smithson & Verkuilen (2006)
  dat$pctRust <- (dat$pctRust*(nrow(dat)-1)+0.5)/nrow(dat)
  return(dat)
    }

#______________________________________________#
### apply bias correction to phyloseq object ###
unbias <- function(phy.in,bias){
  otuTab <- phy.in %>% otu_table %>% data.frame
  #Impute values to replace zeros using a method that retains the correlations structure of the compositional community data
  #this is necessary so that the bias correction can account for variation in the probability of zero counts being real due to 
  #unequal sequencing depth
  if(sum(otuTab==0)>0){otuTab %<>% zCompositions::cmultRepl(method="GBM",output="p-counts",suppress.print = T)}
  #Construct bias correction vector 
  bias.adjust <- rep(1,ntaxa(phy.in)) 
  for(i in 1:ntaxa(phy.in)){
    txx <- taxa_names(phy.in)[i]
    if(txx %in% bias$Taxon){
      bias.adjust[i] <- bias$Bhat[match(txx,bias$Taxon)]
    }
  }
  #Apply Bias correction to OTU table
  otuTab %<>% sweep(.,2,bias.adjust,FUN=`/`) 
  #output phyloseq object with bias corrected OTU table
  otu_table(phy.in) <- otu_table(otuTab, taxa_are_rows = F)
  return(phy.in)
}

#_____________________________________________#
### logit transform data in phyloseq object ###
phy_logit <- function(phy.in){
  #Impute values for zeros
  otuTab <- otu_table(phy.in) %>% data.frame
  if(sum(otuTab==0)>0){otuTab %<>% zCompositions::cmultRepl(method="GBM",output="p-counts",suppress.print = T)
    otu_table(phy) <- otu_table(otu_tab, taxa_are_rows = F)}
  # apply logit transformation
  phy.in %<>% transform_sample_counts(function(x){log(x/sum(x) - log(1-(x/sum(x))))})
  #return updated phylseq object
  return(phy.in)
}

#_________________________________________________________________#
### custom annotation function to add plots as insets in facets ###
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = F, params = list(grob = grob, 
                                       xmin = xmin, xmax = xmax, 
                                       ymin = ymin, ymax = ymax))
}

#______________________________#
### mutate over ragged array ###
mutate_by <- function(.data, group_vars, ...) {
  gvs <- rlang::enquos(group_vars)
  .data %>%
    group_by_at(vars(!!!gvs)) %>%
    mutate(...) %>%
    ungroup
}

#_____________________________#
### calculate geometric mean ###
gm_mean <- function(x, na_rm = FALSE) {
  exp(mean(log(x), na.rm = na_rm))
}

#__________________________________________________#
### rescaling assymetric diverging color palette ###
mid_rescaler <- function(mid = 0) {
  function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
    scales::rescale_mid(x, to, from, mid)
  }
}


