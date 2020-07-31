#!/bin/bash
indexedloc=$1
bamfile1=$2
bamfile2=$3
outputloc=$4
## Run this from folder michel/home/*
sudo -v

## Compute Psi values for control sample
##Spliced Exon
sudo miso --run ${indexedloc}SE_indexmap ${bamfile1} --output-dir ${outputloc}/SE/CEMWT --read-len 100

##Retained Intron
sudo miso --run ${indexedloc}RI_indexmap ${bamfile1} --output-dir ${outputloc}/RI/CEMWT --read-len 100

##A5SS
sudo miso --run ${indexedloc}A5SS_indexmap ${bamfile1} --output-dir ${outputloc}/A5SS/CEMWT --read-len 100

##A3SS
sudo miso --run ${indexedloc}A3SS_indexmap ${bamfile1} --output-dir ${outputloc}/A3SS/CEMWT --read-len 100

##MXE
sudo miso --run ${indexedloc}MXE_indexmap ${bamfile1} --output-dir ${outputloc}/MXE/CEMWT --read-len 100


## Compute Psi values for knockdown sample
##SE
sudo miso --run ${indexedloc}SE_indexmap ${bamfile2} --output-dir ${outputloc}/SE/R30DM --read-len 100 

##RI
sudo miso --run ${indexedloc}RI_indexmap ${bamfile2} --output-dir ${outputloc}/RI/R30DM --read-len 100 

##A5SS
sudo miso --run ${indexedloc}A5SS_indexmap ${bamfile2} --output-dir ${outputloc}/A5SS/R30DM --read-len 100 

##A3SS
sudo miso --run ${indexedloc}A3SS_indexmap ${bamfile2} --output-dir ${outputloc}/A3SS/R30DM --read-len 100 

##MXE
sudo miso --run ${indexedloc}MXE_indexmap ${bamfile2} --output-dir ${outputloc}/MXE/R30DM --read-len 100 

## Summarize the output (only run this once --run finished!)
## This will create a "summary" directory in SE/control/ and in SE/knockdown/

##Spliced Exon
#sudo summarize_miso --summarize-samples ${outputloc}/SE/CEMWT/ ${outputloc}/SE_summ/CEMWT  
#sudo summarize_miso --summarize-samples ${outputloc}/SE/R30dm/ ${outputloc}/SE_summ/R30dm

##Retained Intron
#sudo summarize_miso --summarize-samples MISO_output/hg19/RI/CEMWT/ MISO_output/hg19/RI_summ/CEMWT
#sudo summarize_miso --summarize-samples MISO_output/hg19/RI/R30dm/ MISO_output/hg19/RI_summ/R30dm

##A5SS
#sudo summarize_miso --summarize-samples MISO_output/hg19/A5SS/CEMWT/ MISO_output/hg19/A5SS_summ/CEMWT
#sudo summarize_miso --summarize-samples MISO_output/hg19/A5SS/R30dm/ MISO_output/hg19/A5SS_summ/R30dm

##A3SS
#sudo summarize_miso --summarize-samples MISO_output/hg19/A3SS/CEMWT/ MISO_output/hg19/A3SS_summ/CEMWT
#sudo summarize_miso --summarize-samples MISO_output/hg19/A3SS/R30dm/ MISO_output/hg19/A3SS_summ/R30dm

##MXE
#sudo summarize_miso --summarize-samples MISO_output/hg19/MXE/CEMWT/ MISO_output/hg19/MXE_summ/CEMWT
#sudo summarize_miso --summarize-samples MISO_output/hg19/MXE/R30dm/ MISO_output/hg19/MXE_summ/R30dm


## Detect differentially expressed isoforms between "control" and "knockdown"
## This will compute Bayes factors and delta Psi values between the samples
## and place the results in the directory SE/comparisons/control_vs_knockdown

##Spliced Exon
sudo compare_miso --compare-samples ${outputloc}/SE/CEMWT ${outputloc}/SE/R30DM ${outputloc}/comparison_CEMWTvsR30DM_SE/

##Retained intron
sudo compare_miso --compare-samples ${outputloc}/RI/CEMWT ${outputloc}/RI/R30DM ${outputloc}/comparison_CEMWTvsR30DM_RI/

##A5SS
sudo compare_miso --compare-samples ${outputloc}/A5SS/CEMWT ${outputloc}/A5SS/R30DM ${outputloc}/comparison_CEMWTvsR30DM_A5SS/

##A3SS
sudo compare_miso --compare-samples ${outputloc}/A3SS/CEMWT ${outputloc}/A3SS/R30DM ${outputloc}/comparison_CEMWTvsR30DM_A3SS/

##MXE
sudo compare_miso --compare-samples ${outputloc}/MXE/CEMWT ${outputloc}/MXE/R30DM ${outputloc}/comparison_CEMWTvsR30DM_MXE/
