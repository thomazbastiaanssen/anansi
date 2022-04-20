makeNotes: 
	$(info make clean = delete all output to start over)
	$(info make demux = demultiplex raw MiSeq data with Pheniqs)
	$(info make trim = remove primers and sequencing adapters) 
	$(info make denoise = denoise with DADA2)
	$(info make compile = merge libraries, add metadata and taxonomy)
	$(info make analysis = run analyses in R and make figures for publication)
	
clean:
	rm -r output
	$(info you have a clean slate)
	
##############################
### demultiplex MiSeq data ###
##############################

demux: output/demux/

output/demux/: code/demux.config.json 
	mkdir -p output/demux
	pheniqs mux -R output/demux/demux.report.txt -c $< --base-input data/MiSeq_raw --base-output output/demux/
	touch output/demux/

###############################
### trim primers / adapters ###
###############################

trim: output/trim/

output/trim/: code/trim.sh output/demux/
	mkdir -p output/trim
	$<
	touch output/trim/

##############
## denoise ###
##############

denoise: output/dada/

output/dada/: code/denoise.R output/trim/
	mkdir -p output/dada
	Rscript $< 
	touch output/dada/

##############################
### compile processed data ###
##############################

compile: output/compiled/

output/compiled/: code/compile.R output/dada/ data/Sample_data.csv
	mkdir -p output/compiled/
	Rscript $<
	touch output/compiled/

################
### Analysis ###
################

tb := output/tabs/
fg := output/figs/
rds := output/rds/

analysis: analysisOut \
	${tb}bias.csv ${fg}Fig.S2.jpg \
	${rds}mv.genotype.rds ${rds}mv.region.rds \
	${fg}Fig.2.pdf \
	${fg}Fig.3.pdf \
	${fg}Fig.4.pdf \
	${fg}Fig.S3.jpg ${tb}susceptibility.csv \
	${fg}Fig.S4.jpg \
	${fg}Fig.S1.jpg
	
analysisOut:
	mkdir -p output/figs
	mkdir -p output/tabs
	mkdir -p output/rds

# Estimate sequencing bias from mock community data
${tb}bias.csv ${fg}Fig.S2.jpg: code/biasEstimates.R output/compiled/ | analysisOut
	Rscript $<

### Community analyses ###
# Joint-species distribution models
${rds}mv.genotype.rds ${rds}mv.region.rds: code/jsdModels.R ${tb}bias.csv | analysisOut
	Rscript $<
# Multipanel figure showing variation in community composition
${fg}Fig.2.pdf: code/communityFigure.R ${rds}mv.genotype.rds ${rds}mv.region.rds | analysisOut
	Rscript $<
# Single-species priority effects
${fg}Fig.3.pdf: code/priorityEffects.R ${tb}bias.csv | analysisOut
	Rscript $<

### Rust analyses ###
${fg}Fig.S3.jpg ${tb}susceptibility.csv: code/rustSusceptibility.R | analysisOut
	Rscript $<
${fg}Fig.4.pdf: code/rustAnalyses.R code/rustSusceptibility.R | analysisOut
	Rscript $<
${fg}Fig.S4.jpg: code/rustCor.R | analysisOut
	Rscript $<
	
### Make map of geontype origins for supplement
${fg}Fig.S1.jpg: code/mapFigS1.R | analysisOut
	Rscript $<


