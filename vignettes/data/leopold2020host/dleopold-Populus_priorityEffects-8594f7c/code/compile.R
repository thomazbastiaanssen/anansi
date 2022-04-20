#Compile MiSeq library and sample data into single phyloseq object 
#Also add taxonomy and remove non-target sequences
 
library(tidyverse)
library(magrittr)
library(foreach)
library(dada2)
library(Biostrings)
library(phyloseq)
library(DECIPHER)

#read in denoised sequence table
seqTab <- readRDS("output/dada/seqTab.rds")
#read in sample data
meta <- read.csv("data/Sample_data.csv",as.is=T)

#merge into phyloseq object and give OTUs arbitrary names
seqs <- getSequences(seqTab) %>% DNAStringSet()
colnames(seqTab) <- paste0("OTU.",1:ncol(seqTab))
names(seqs) <- colnames(seqTab)
meta %<>% column_to_rownames("SampID")
phy <- phyloseq(otu_table(seqTab, taxa_are_rows = F),
                sample_data(meta),
                refseq(seqs))

#remove host contamination by aligning against genome with Bowtie2
#Requires a precompiled bowtie database file created from the Populus trichocarpa v3 genome
dir.create("output/scratch")
refseq(phy) %>% DNAStringSet %>% writeXStringSet("output/scratch/otus.fa")
system("/home/harpua/miniconda3/bin/bowtie2 --threads 10 -x data/referenceDBs/Ptri.v.3.db -f output/scratch/otus.fa --un output/scratch/noHost.fa --al output/scratch/host.fa --very-sensitive")
hostFiltered <- readDNAStringSet("output/scratch/noHost.fa") %>% names
phy %<>% prune_taxa(hostFiltered, .)
unlink("output/scratch",recursive=T)

#collapse sequences with >= 99% similarity
cluster <- function(phy.in,method="single",dissimilarity=0.01){
  clusters <- DECIPHER::DistanceMatrix(refseq(phy.in), includeTerminalGaps = T, processors=NULL) %>%
    DECIPHER::IdClusters(method=method, cutoff=dissimilarity, processors=NULL) 
  clusters[,1] %<>% as.character()
  tax_table(phy.in) <- clusters %>% as.matrix %>% tax_table
  phy.in %<>% speedyseq::tax_glom("cluster")
  tax_table(phy.in) <- NULL
  return(phy.in)
}  
phy %<>% cluster  

#Identify expected taxa (within 99% of expected species)
focalSeqs <- readDNAStringSet("data/SangerSeqs/dadaPriors.fasta")
seqDist <- append(focalSeqs,refseq(phy)) %>% DistanceMatrix(verbose=F, processors=NULL) %>%
  .[!(rownames(.) %in% names(focalSeqs)) , (rownames(.) %in% names(focalSeqs))]
toMerge <- foreach(i=1:dim(seqDist)[1], .combine=rbind) %do% {
  minDist <- min(seqDist[i,])
  if(minDist > 0.02){return(c(rownames(seqDist)[i],NA))}
  if(minDist <= 0.02){return(c(rownames(seqDist)[i],colnames(seqDist)[seqDist[i,]==minDist]))}
} %>% data.frame(stringsAsFactors = F) %>% drop_na()
toMerge$X2[toMerge$X2=="MelampsoraB"] <- "Melampsora"
for(i in unique(toMerge$X2)){
  otus <- toMerge$X1[toMerge$X2==i]
  taxa_names(phy)[taxa_names(phy)==otus[1]] <- i; otus[1] <- i
  phy %<>% merge_taxa(otus)}

#assign taxonomy with database including expected taxa
tax <- assignTaxonomy(refseq(phy),"data/referenceDBs/sh_general_release_dynamic_s_04.02.2020.fasta.gz",multithread=T)
rownames(tax) <- taxa_names(phy)
tax_table(phy) <- tax_table(tax)

#remove OTUs not assigned to Order
phy %<>% subset_taxa(!is.na(Order))

#Remove putative contam (taxa more preveent or proportionally abundanct in negative control samples)
#After this filtering a number of unexpected OTUs remain. These appeaar to be primarily associated with uninoculated control plants. Their 
contam.dat <- data.frame(OTU=taxa_names(phy),
                         prev.neg = phy %>% prune_samples(grepl("NEG",sample_names(.)),.) %>% 
                           otu_table %>% data.frame %>% apply(2,function(x){sum(x>0)/length(x)}),
                         prev.samps = phy %>% prune_samples(grepl("^G",sample_names(.)),.) %>% 
                           otu_table %>% data.frame %>% apply(2,function(x){sum(x>0)/length(x)}),
                         mean.neg = phy %>% prune_samples(grepl("NEG",sample_names(.)),.) %>% 
                           otu_table %>% data.frame %>% apply(2,function(x){ifelse(sum(x)==0,0, sum(x)/sum(x>0))}),
                         mean.samps = phy %>%  prune_samples(grepl("^G",sample_names(.)),.) %>% 
                           otu_table %>% data.frame %>% apply(2,function(x){ifelse(sum(x)==0,0, sum(x)/sum(x>0))}))
contam.dat %<>% mutate(contam=(prev.neg>prev.samps|mean.neg>mean.samps))
phy %<>% prune_taxa(!contam.dat$contam,.)

#remove negative control samps
phy %<>% prune_samples(!grepl("NEG",sample_names(.)),.)

#' ### Look at sequencing depth
minDepth <- 5000
data.frame(SeqDepth=sample_sums(phy), TP=factor(sample_data(phy)$Timepoint),Type=sample_data(phy)$Treatment ) %>%
  mutate(cutoff=SeqDepth>minDepth,neg=ifelse(Type=="Negative","neg","pos"),TP=paste(TP,neg)) %>%
  ggplot(aes(x=TP, y=SeqDepth)) +
  geom_violin() +
  geom_point(aes(color=cutoff),position=position_jitter(width=0.1)) +
  theme_classic()
#' ### Remove samples below sequencing depth cutoff
(phy %<>% prune_samples(sample_sums(.)>minDepth,.))

###############
### OUTPUTs ###
###############

#Save phyloseq object
saveRDS(phy,"output/compiled/phy.rds")
#Save components for manual inspection
otu_table(phy) %>% write.csv("output/compiled/OTU.table.csv")
tax_table(phy) %>% write.csv("output/compiled/taxonomy.table.csv")
