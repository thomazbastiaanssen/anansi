#!/bin/bash

#This script will remove gene primers from both ends of paired-end read prior to denoising
#requires cutadapt (v2.7+) https://cutadapt.readthedocs.io/en/stable/
#requires SeqPurge https://github.com/imgag/ngs-bits

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#Steps marked with alternating hash/dash require adjustment for each data set

#make directory for temporary files
scratch="output/scratch/"
mkdir -p $scratch

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#set gene primers for trimming
fwd_primer_fungi="CTTGGTCATTTAGAGGAAGTAA" #ITS1F
rc_fwd_primer_fungi="TTACTTCCTCTAAATGACCAAG"
rev_primer_fungi="CAAAGAYTCGATGATTCAC" #glITS7r
rc_rev_primer_fungi="GTGAATCATCGARTCTTTG"

#Make empty file for appending summary data
log="output/trim/summary.tab"
echo -e "Sample\tRaw\tFungi.trim1\tFungi.trim2" > ${log}
 
#Loop over samples
for FWD in $( find output/demux/ -name "*R1.fastq.gz" | grep -v "undetermined"); do
    REV=`echo $FWD | sed 's/R1.fastq.gz/R2.fastq.gz/'`
    
    #identify current sample 
    SAMP=`echo $FWD | sed 's/.R1.fastq.gz//' | sed 's/output\/demux\///'`  

    #trim R1 gene primer
    FWDtrim1="${scratch}${SAMP}.R1.trim1.fq.gz"
    REVtrim1="${scratch}${SAMP}.R2.trim1.fq.gz"
    cutadapt --quiet -g $fwd_primer_fungi -G $rev_primer_fungi --discard-untrimmed --overlap 20 -e 0.17 --minimum-length 223:226 --maximum-length 226:229 -o $FWDtrim1 -p $REVtrim1 $FWD $REV

    #Count seqs after first trim step 
    if [[ -f $FWDtrim1 && -f $REVtrim1 ]]; then
        trim1Ct=`gzip -cd $FWDtrim1 | grep -c '^@M01498'`
    else
        trim1Ct=0
    fi

    #trim read through contamination
    if [ $trim1Ct -gt 0 ]; then
        FWDtrim2="${scratch}${SAMP}.R1.trim2.fq.gz"
        REVtrim2="${scratch}${SAMP}.R2.trim2.fq.gz"
        SeqPurge -in1 $FWDtrim1 -in2 $REVtrim1 -out1 $FWDtrim2 -out2 $REVtrim2 -a1 $rc_rev_primer_fungi -a2 $rc_fwd_primer_fungi -qcut 0 -ncut 0 -min_len 100 -summary ${scratch}${SAMP}.log.txt
    fi

    #Count seqs after final step 
    if [[ -f $FWDtrim2 && -f $REVtrim2 ]]; then
        trim2Ct=`gzip -cd $FWDtrim2 | grep -c '^@M01498'`
    else
        trim2Ct=0
    fi

    #Keep only samples that made it through the processing (also truncate reads to remove low quality tails)
    if [ $trim2Ct -gt 0 ]; then
        FWDout="output/trim/${SAMP}.R1.fq.gz"
        REVout="output/trim/${SAMP}.R2.fq.gz"
        #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
        #set an appropriate --length setting here depending on your data 
        cutadapt --quiet -g XX --length 200 -o $FWDout $FWDtrim2
        cutadapt --quiet -g XX --length 175 -o $REVout $REVtrim2
    fi

    #make summary and append to file
    RAW=`gzip -cd $FWD | grep -c '^@M01498'`
    echo -e $SAMP"\t"$RAW"\t"$trim1Ct"\t"$trim2Ct | cat >> ${log}

    #remove tmp files
    rm ${scratch}${SAMP}.*
done

#delete scratch folder
rm -r ${scratch}


