#denoise MiSeq data using dada2 

#load packages
library(tidyverse)
library(magrittr)
library(dada2)
library(ShortRead)
library(Biostrings)

#path for output
out.path <- "output/dada"

#Idenify trimmed fwd and rev reads
in.path <- file.path("output/trim")
path.fwd <- in.path %>% list.files(pattern="R1.fq.gz",full.names = T) %>% sort 
path.rev <- in.path %>% list.files(pattern="R2.fq.gz",full.names = T) %>% sort 
  
#make file names / paths for trimmed and filtered files
filt.path <- file.path(out.path, "filt")
filt.fwd <- path.fwd %>% gsub(in.path,filt.path,.)
filt.rev <- path.rev %>% gsub(in.path,filt.path,.)
  
#preview read quality of 12 randomly selected samples before trimming
qual.samps <- sample(1:length(path.fwd),12)
qual.fwd <- plotQualityProfile(path.fwd[qual.samps]) + ggtitle("Qual profiles fwd")
file.path(out.path, "qual.fwd.pdf") %>% ggsave(width=12,height=9)
qual.rev <- plotQualityProfile(path.rev[qual.samps]) + ggtitle("Qual profiles rev")
file.path(out.path, "qual.rev.pdf") %>% ggsave(width=12,height=9)

#Quality filter 
trim <- filterAndTrim(path.fwd, filt.fwd,
                      path.rev, filt.rev,
                      maxN=0, maxEE=2, 
                      multithread=TRUE)
  
#Update list of trimmed file paths to exclude samples with no reads passing filters
filt.fwd <- list.files(filt.path, pattern = "R1.fq.gz",full.names = T)
filt.rev <- list.files(filt.path, pattern = "R2.fq.gz",full.names = T)
  
#Check quality of trimmed and filtered reads
qual.samps <- sample(1:length(filt.fwd),12)
qual.filt.fwd <- plotQualityProfile(filt.fwd[qual.samps]) + ggtitle("Qual profiles fwd filtered")
file.path(out.path, "qual.fwd.filtered.pdf") %>% ggsave(width=12,height=9)
qual.filt.rev <- plotQualityProfile(filt.rev[qual.samps]) + ggtitle("Qual profiles rev filtered")
file.path(out.path, "qual.rev.filtered.pdf") %>% ggsave(width=12,height=9)
  
#dereplicate
derep.fwd <- derepFastq(filt.fwd, verbose=F)
derep.rev <- derepFastq(filt.rev, verbose=F) 
  
#Trim names of derep objects to just sample names 
names(derep.fwd) %<>% gsub(".R1.fq.gz","",.)
names(derep.rev) %<>% gsub(".R2.fq.gz","",.)
  
#Learn errors 
err.fwd <- learnErrors(filt.fwd, multithread=TRUE, nbases = 5e08) 
err.plot.fwd <- plotErrors(err.fwd, nominalQ=TRUE) + ggtitle("Forward reads error model")
file.path(out.path, "errMod.fwd.pdf") %>% ggsave(width=5,height=5)
err.rev <- learnErrors(filt.rev, multithread=TRUE, nbases = 5e08)
err.plot.rev <- plotErrors(err.rev, nominalQ=TRUE) + ggtitle("Reverse reads error model")
file.path(out.path, "errMod.rev.pdf") %>% ggsave(width=5,height=5)
    
#Import known sequences for the expected taxa to use as priors during denoising.
priors <- readDNAStringSet("data/SangerSeqs/dadaPriors.fasta")
#Trim to same length as truncation applied in the last step of the trim.sh script
priors.fwd.200bp <- priors %>% substr(1,200) %>% DNAStringSet()
priors.rev.175bp <- priors %>% dada2:::rc() %>% substr(1,175) %>% DNAStringSet()
 
#Denoise 
dada.fwd <- dada(derep.fwd, err=err.fwd, multithread=TRUE, priors=as.character(priors.fwd.200bp))
dada.rev <- dada(derep.rev, err=err.rev, multithread=TRUE, priors=as.character(priors.rev.175bp))
  
#Merge reads
merged <- mergePairs(dada.fwd, derep.fwd, dada.rev, derep.rev, trimOverhang = T)
    
#Make sequence table
seqtab <- merged %>% makeSequenceTable %>% collapseNoMismatch
    
#Remove chimeras
seqtab.nonChimeras <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE,verbose=T)
    
#Make summary report of what happened
getN <- function(x) sum(getUniques(x))
trim.summary <- trim %>% data.frame %>% rownames_to_column("Sample") 
trim.summary$Sample %<>% gsub(".R1.fq.gz","",.)
track <- cbind(sapply(dada.fwd, getN), 
               sapply(dada.rev, getN), 
               sapply(merged, getN),
               rowSums(seqtab.nonChimeras)) %>% 
  data.frame %>% rownames_to_column("Sample") 
track %<>% left_join(trim.summary,.)
colnames(track) <- c("sample", "input", "filtered", "denoised.fwd", "denoised.rev", "merged", "ChimeraFiltered")
file.path(out.path,"dadaSummary.csv") %>% write.csv(track,.,row.names = F)
  
# Save output
file.path(out.path,"seqTab.rds") %>% saveRDS(seqtab.nonChimeras,.)
