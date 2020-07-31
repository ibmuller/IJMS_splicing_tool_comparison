#library(Rcpp)
#library(plyr)
#library(dplyr)
#library(tidyr)
#library(ggplot2)
#library(tibble)
#library(splitstackshape)
#library(stringr)

#unfiltered
MATS_SE <- read.delim("SE.MATS.JCEC.txt", header=TRUE)
MISO_SE <- read.delim("CEMWT_vs_R30DM_SE.miso_bf")
SUPPA_SE <- read.delim("CEMWT_r30dm_diffSplice.dpsi.temp.0", header=TRUE)

subset_MATS <- subset(MATS_SE, MATS_SE$FDR<0.05)

##change MATS coords
MATS_SE$exonStart_0base <- as.numeric(as.character(MATS_SE[,6])) + 1
MATS_SE$upstreamES <- as.numeric(as.character(MATS_SE[,8])) + 1
MATS_SE$downstreamES <- as.numeric(as.character(MATS_SE[,10])) + 1

#create event for MATS that matches MISO
MATS_SE <- cbind(MATS_SE, as.data.frame(paste0(MATS_SE$chr,":",MATS_SE$downstreamES,":",MATS_SE$downstreamEE,":-@",MATS_SE$chr,":",MATS_SE$exonStart_0base,":",MATS_SE$exonEnd,":-@",MATS_SE$chr,":",MATS_SE$upstreamES,":",MATS_SE$upstreamEE,":", MATS_SE$strand)))
colnames(MATS_SE)[24] <- "event_name" 


#merge on event_name
miso_mats <- merge(MISO_SE, MATS_SE, by = "event_name", all=FALSE)

miso_mats_sig <- subset(miso_mats, miso_mats$bayes_factor > 10 & miso_mats$FDR < 0.05)

#miso_mats_sig <- subset(miso_mats, miso_mats$bayes_factor > 10)
#miso_mats_sig <- subset(miso_mats_sig, miso_mats_sig$FDR < 0.05)


#miso_mats_sig <- subset(miso_mats, miso_mats$bayes < 1000 & miso_mats$bayes >= 10 )
#miso_mats_sig <- subset(miso_mats_sig, miso_mats_sig$PValue < 0.05)

##MISO
#split coordinte by ":" and transpose
miso_split <- t(as.data.frame(strsplit(as.character(MISO_SE$event_name), "-|\\+|\\:")))

#change column names and bind to original MISO file
colnames(miso_split) <-as.character(c("chrom", "upstreamES", "upstreamEE","V4", "chr.x", "exonStart_0base", "exonEnd", "V8", "chr.y", "downstreamES", "downstreamEE", "V12"))
new_miso <- cbind(MISO_SE, miso_split)
colnames(new_miso) <- make.unique(names(new_miso))
#new_miso <- separate(data = new_miso, col= "upstreamEE", into = c("riExonStart_0base", "upstreamEE"), sep="\\|", fill = "right")

##SUPPA
#import table with SUPPA2 dpsi values and change column names
colnames(SUPPA_SE) <- c("event", "dPSI","p_value")

#split columns and separate Ensembl ID 
SUPPA_SE <- cbind.data.frame(str_split_fixed(SUPPA_SE$event, ';', 2), SUPPA_SE$dPSI, SUPPA_SE$p_value)
colnames(SUPPA_SE) <- c("event", "coord", "dPSI","p_value")

#split columns and separate coordinate)
SUPPA_SE <- cbind.data.frame(SUPPA_SE$event, str_split_fixed(SUPPA_SE$coord, ':', 2), SUPPA_SE$dPSI, SUPPA_SE$p_value)
colnames(SUPPA_SE) <- c("ensembl_gene_id", "event_type", "event_name_suppa", "dPSI","p_value")

#Remove all rows with NA
SUPPA_SE_omit <- na.omit(SUPPA_SE)
SUPPA_SE_omit <- SUPPA_SE_omit[grep("SE", SUPPA_SE_omit$event_type),]

#put chr before event    
SUPPA_SE_omit$event_name_suppa <- paste0("chr", SUPPA_SE_omit$event_name_suppa)

#Get Gene list from BioMart and match genesymbols with ensembl IDs
#library('biomaRt')
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#gene_list_mart <- getBM(filters= "ensembl_gene_id",attributes= c("hgnc_symbol", "ensembl_gene_id"),values= SUPPA_SE$ensembl_gene_id, mart= mart, uniqueRows = TRUE)
#SUPPA_SE <- merge(dpsi,gene_list_mart)

#colnames(new_miso)[15] <- "chrom.x"

## SUPPA MISO ## merge on event_name
#remove first column name and create suppa event coords SE
rownames(new_miso) <- NULL
new_miso <- cbind(new_miso, as.data.frame(paste0(new_miso$chrom,":",new_miso$upstreamEE,"-",new_miso$exonStart_0base,":",new_miso$exonEnd,"-",new_miso$downstreamES,":",new_miso$strand)))
colnames(new_miso)[31] <- "event_name_suppa"

#colnames(new_miso)[15] <- "chrom.y"
#colnames(new_miso)[16] <- "chrom.x"
miso_suppa <- merge(new_miso, SUPPA_SE_omit, by = "event_name_suppa", all=FALSE)
miso_suppa_sig <- subset(miso_suppa, miso_suppa$bayes_factor > 10)
miso_suppa_sig <- subset(miso_suppa_sig, miso_suppa_sig$p_value < 0.05)


################ MATS SUPPA
#mats_suppa <- cbind(MATS_SE, as.data.frame(paste0(MATS_SE$chr,":",MATS_SE$longExonEnd,"-",MATS_SE$flankingES,":",MATS_SE$shortEE,"-",MATS_SE$flankingES,":",MATS_SE$strand)))
mats_suppa <- cbind(MATS_SE, as.data.frame(paste0(MATS_SE$chr,":",MATS_SE$upstreamEE,"-",MATS_SE$exonStart_0base,":",MATS_SE$exonEnd,"-",MATS_SE$downstreamES,":",MATS_SE$strand)))
colnames(mats_suppa)[25] <- "event_name_suppa"
mats_suppa <- merge(mats_suppa, SUPPA_SE_omit, by = "event_name_suppa", all=FALSE)
mats_suppa_sig <- subset(mats_suppa, mats_suppa$FDR < 0.05)
mats_suppa_sig <- subset(mats_suppa_sig, mats_suppa_sig$p_value < 0.05)

#significant hits
miso_se_sig <- subset(MISO_SE, MISO_SE$bayes_factor > 10)
mats_se_sig <- subset(MATS_SE, MATS_SE$FDR < 0.05)
suppa_se_sig <- subset(SUPPA_SE_omit, SUPPA_SE_omit$p_value < 0.05)


######STATISTICS
sink('output_SE.txt')

#print events
print("MATS SE + FDR < 0.05")
dim(MATS_SE)
dim(mats_se_sig)
print("MISO SE")
dim(MISO_SE)
dim(miso_se_sig)
print("SUPPA SE")
dim(SUPPA_SE_omit)
dim(suppa_se_sig)
print("MISO MATS")
dim(miso_mats)
dim(miso_mats_sig)
print("MISO SUPPA")
dim(miso_suppa)
dim(miso_suppa_sig)
print("MATS SUPPA")
dim(mats_suppa)
dim(mats_suppa_sig)

print("CORRELATION MISO MATS")
lm_miso_mats <-lm(miso_mats$IncLevelDifference~miso_mats$diff, data = miso_mats)
summary(lm_miso_mats)

lm_miso_mats_sig <-lm(miso_mats_sig$IncLevelDifference~miso_mats_sig$diff, data = miso_mats_sig)
summary(lm_miso_mats_sig)

print("CORRELATION MISO SUPPA")
lm_miso_suppa <-lm(miso_suppa$diff~miso_suppa$dPSI, data = miso_suppa)
summary(lm_miso_suppa)

lm_miso_suppa_sig <-lm(miso_suppa_sig$diff~miso_suppa_sig$dPSI, data = miso_suppa_sig)
summary(lm_miso_suppa_sig)

print("CORRELATION MATS SUPPA")
lm_mats_suppa <-lm(mats_suppa$dPSI~mats_suppa$IncLevelDifference, data = mats_suppa)
summary(lm_mats_suppa)

lm_mats_suppa_sig <-lm(mats_suppa_sig$dPSI~mats_suppa_sig$IncLevelDifference, data = mats_suppa_sig)
summary(lm_mats_suppa_sig)

sink()

rm(list = ls())

#unfiltered
MATS_RI <- read.delim("RI.MATS.JCEC.txt", header=TRUE)
MISO_RI <- read.delim("CEMWT_vs_R30DM_RI.miso_bf")
SUPPA_RI <- read.delim("CEMWT_r30dm_diffSplice.dpsi.temp.0", header=TRUE)

##change MATS coords
MATS_RI$riExonStart_0base <- as.numeric(as.character(MATS_RI[,6])) + 1
MATS_RI$upstreamES <- as.numeric(as.character(MATS_RI[,8])) + 1
MATS_RI$downstreamES <- as.numeric(as.character(MATS_RI[,10])) + 1

#create event for MATS that matches MISO
MATS_RI <- cbind(MATS_RI, as.data.frame(paste0(MATS_RI$chr,":",MATS_RI$riExonStart_0base,"-",MATS_RI$upstreamEE,":",MATS_RI$strand,"@",MATS_RI$chr,":",MATS_RI$downstreamES,"-",MATS_RI$downstreamEE,":",MATS_RI$strand)))
colnames(MATS_RI)[24] <- "event_name" 


#merge on event_name
miso_mats <- merge(MISO_RI, MATS_RI, by = "event_name", all=FALSE)
miso_mats_sig <- subset(miso_mats, miso_mats$bayes_factor > 10 & miso_mats$FDR < 0.05)

##MISO
#split coordinte by ":" and transpose
miso_split <- t(as.data.frame(strsplit(as.character(MISO_RI$event_name), "-|\\+|\\:")))

#change column names and bind to original MISO file
colnames(miso_split) <-as.character(c("chrom", "riExonStart_0base", "upstreamEE","V4", "@chr", "downstreamES", "downstreamEE", "V7"))
new_miso <- cbind(MISO_RI, miso_split)
colnames(new_miso) <- make.unique(names(new_miso))
#new_miso <- separate(data = new_miso, col= "upstreamEE", into = c("riExonStart_0base", "upstreamEE"), sep="\\|", fill = "right")

##SUPPA
#import table with SUPPA2 dpsi values and change column names
colnames(SUPPA_RI) <- c("event", "dPSI","p_value")

#split columns and separate Ensembl ID 
SUPPA_RI <- cbind.data.frame(str_split_fixed(SUPPA_RI$event, ';', 2), SUPPA_RI$dPSI, SUPPA_RI$p_value)
colnames(SUPPA_RI) <- c("event", "coord", "dPSI","p_value")

#split columns and separate coordinate)
SUPPA_RI <- cbind.data.frame(SUPPA_RI$event, str_split_fixed(SUPPA_RI$coord, ':', 2), SUPPA_RI$dPSI, SUPPA_RI$p_value)
colnames(SUPPA_RI) <- c("ensembl_gene_id", "event_type", "event_name_suppa", "dPSI","p_value")

#Remove all rows with NA
SUPPA_RI_omit <- na.omit(SUPPA_RI)
SUPPA_RI_omit <- SUPPA_RI_omit[grep("RI", SUPPA_RI_omit$event_type),]

#put chr before event    
SUPPA_RI_omit$event_name_suppa <- paste0("chr", SUPPA_RI_omit$event_name_suppa)

#Get Gene list from BioMart and match genesymbols with ensembl IDs
#library('biomaRt')
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#gene_list_mart <- getBM(filters= "ensembl_gene_id",attributes= c("hgnc_symbol", "ensembl_gene_id"),values= SUPPA_RI$ensembl_gene_id, mart= mart, uniqueRows = TRUE)
#SUPPA_RI <- merge(dpsi,gene_list_mart)

#colnames(new_miso)[15] <- "chrom.x"

## SUPPA MISO ## merge on event_name
#remove first column name and create suppa event coords RI
rownames(new_miso) <- NULL
new_miso <- cbind(new_miso, as.data.frame(paste0(new_miso$chrom,":",new_miso$riExonStart_0base,":",new_miso$upstreamEE,"-",new_miso$downstreamES,":",new_miso$downstreamEE,":",new_miso$strand)))
colnames(new_miso)[27] <- "event_name_suppa"

#colnames(new_miso)[15] <- "chrom.y"
#colnames(new_miso)[16] <- "chrom.x"
miso_suppa <- merge(new_miso, SUPPA_RI_omit, by = "event_name_suppa", all=FALSE)
miso_suppa_sig <- subset(miso_suppa, miso_suppa$bayes_factor > 10 & miso_suppa$p_value < 0.05)

################ MATS SUPPA
#mats_suppa <- cbind(MATS_RI, as.data.frame(paste0(MATS_RI$chr,":",MATS_RI$longExonEnd,"-",MATS_RI$flankingES,":",MATS_RI$shortEE,"-",MATS_RI$flankingES,":",MATS_RI$strand)))
mats_suppa <- cbind(MATS_RI, as.data.frame(paste0(MATS_RI$chr,":",MATS_RI$riExonStart_0base,":",MATS_RI$upstreamEE,"-",MATS_RI$downstreamES,":",MATS_RI$downstreamEE,":",MATS_RI$strand)))
colnames(mats_suppa)[25] <- "event_name_suppa"
mats_suppa <- merge(mats_suppa, SUPPA_RI_omit, by = "event_name_suppa", all=FALSE)
mats_suppa_sig <- subset(mats_suppa, mats_suppa$FDR < 0.05)
mats_suppa_sig <- subset(mats_suppa_sig, mats_suppa_sig$p_value < 0.05)

#significant hits
miso_ri_sig <- subset(MISO_RI, MISO_RI$bayes_factor > 10)
mats_ri_sig <- subset(MATS_RI, MATS_RI$FDR < 0.05)
suppa_ri_sig <- subset(SUPPA_RI_omit, SUPPA_RI_omit$p_value < 0.05)



######STATISTICS
sink('output_RI.txt')

#print events
print("MATS RI + FDR < 0.05")
dim(MATS_RI)
dim(mats_ri_sig)
print("MISO RI")
dim(MISO_RI)
dim(miso_ri_sig)
print("SUPPA RI")
dim(SUPPA_RI_omit)
dim(suppa_ri_sig)
print("MISO MATS")
dim(miso_mats)
dim(miso_mats_sig)
print("MISO SUPPA")
dim(miso_suppa)
dim(miso_suppa_sig)
print("MATS SUPPA")
dim(mats_suppa)
dim(mats_suppa_sig)

print("CORRELATION MISO MATS")
lm_miso_mats <-lm(miso_mats$IncLevelDifference~miso_mats$diff, data = miso_mats)
summary(lm_miso_mats)

lm_miso_mats_sig <-lm(miso_mats_sig$IncLevelDifference~miso_mats_sig$diff, data = miso_mats_sig)
summary(lm_miso_mats_sig)

print("CORRELATION MISO SUPPA")
lm_miso_suppa <-lm(miso_suppa$diff~miso_suppa$dPSI, data = miso_suppa)
summary(lm_miso_suppa)

lm_miso_suppa_sig <-lm(miso_suppa_sig$diff~miso_suppa_sig$dPSI, data = miso_suppa_sig)
summary(lm_miso_suppa_sig)

print("CORRELATION MATS SUPPA")
lm_mats_suppa <-lm(mats_suppa$dPSI~mats_suppa$IncLevelDifference, data = mats_suppa)
summary(lm_mats_suppa)

lm_mats_suppa_sig <-lm(mats_suppa_sig$dPSI~mats_suppa_sig$IncLevelDifference, data = mats_suppa_sig)
summary(lm_mats_suppa_sig)

sink()

rm(list = ls())

#unfiltered
MATS_A5 <- read.delim("A5SS.MATS.JCEC.txt", header=TRUE)
MISO_A5 <- read.delim("CEMWT_vs_R30DM_A5.miso_bf")
SUPPA_A5 <- read.delim("CEMWT_r30dm_diffSplice.dpsi.temp.0", header=TRUE)

##MATS
##change MATS coords
MATS_A5$longExonStart_0base <- as.numeric(as.character(MATS_A5[,6])) + 1
MATS_A5$flankingES <- as.numeric(as.character(MATS_A5[,10])) + 1

#create event for MATS that matches MISO
MATS_A5 <- cbind(MATS_A5, as.data.frame(paste0 (MATS_A5$chr,":",MATS_A5$longExonStart_0base, ":", MATS_A5$shortEE,"|",MATS_A5$longExonEnd,":", MATS_A5$strand, "@", MATS_A5$chr,":",MATS_A5$flankingES,":",MATS_A5$flankingEE,":",MATS_A5$strand)))
colnames(MATS_A5)[24] <- "event_name" 

#merge on event_name
miso_mats <- merge(MATS_A5, MISO_A5, by = "event_name", all=FALSE)
miso_mats_sig <- subset(miso_mats, miso_mats$bayes_factor > 10)
miso_mats_sig <- subset(miso_mats_sig, miso_mats_sig$PValue < 0.05)

##MISO
#split coordinte by ":" and transpose
miso_split <- t(as.data.frame(strsplit(as.character(MISO_A5$event_name), "-|\\+|\\:")))

#change column names and bind to original MISO file
colnames(miso_split) <-as.character(c("chrom", "riExonStart_0base", "upstreamEE","V4", "@chr", "downstreamES", "downstreamEE", "V7"))
new_miso <- cbind(MISO_A5, miso_split)
colnames(new_miso) <- make.unique(names(new_miso))
new_miso <- separate(data = new_miso, col= "upstreamEE", into = c("riExonStart_0base", "upstreamEE"), sep="\\|", fill = "right")

##SUPPA
#import table with SUPPA2 dpsi values and change column names
colnames(SUPPA_A5) <- c("event", "dPSI","p_value")

#split columns and separate Ensembl ID 
SUPPA_A5 <- cbind.data.frame(str_split_fixed(SUPPA_A5$event, ';', 2), SUPPA_A5$dPSI, SUPPA_A5$p_value)
colnames(SUPPA_A5) <- c("event", "coord", "dPSI","p_value")

#split columns and separate coordinate)
SUPPA_A5 <- cbind.data.frame(SUPPA_A5$event, str_split_fixed(SUPPA_A5$coord, ':', 2), SUPPA_A5$dPSI, SUPPA_A5$p_value)
colnames(SUPPA_A5) <- c("ensembl_gene_id", "event_type", "event_name_suppa", "dPSI","p_value")

#Remove all rows with NA
SUPPA_A5_omit <- na.omit(SUPPA_A5)
SUPPA_A5_omit <- SUPPA_A5_omit[grep("A5", SUPPA_A5_omit$event_type),]

#put chr before event    
SUPPA_A5_omit$event_name_suppa <- paste0("chr", SUPPA_A5_omit$event_name_suppa)

#Get Gene list from BioMart and match genesymbols with ensembl IDs
#library('biomaRt')
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#gene_list_mart <- getBM(filters= "ensembl_gene_id",attributes= c("hgnc_symbol", "ensembl_gene_id"),values= SUPPA_A5$ensembl_gene_id, mart= mart, uniqueRows = TRUE)
#SUPPA_A5 <- merge(dpsi,gene_list_mart)

#colnames(new_miso)[15] <- "chrom.x"

## SUPPA MISO ## merge on event_name
#remove first column name and create suppa event coords A5
rownames(new_miso) <- NULL
new_miso <- cbind(new_miso, as.data.frame(paste0(new_miso$chrom,":",new_miso$upstreamEE,"-",new_miso$downstreamES, ":", new_miso$riExonStart_0base,"-",new_miso$downstreamES,":",new_miso$strand)))
colnames(new_miso)[27] <- "event_name_suppa"


#colnames(new_miso)[15] <- "chrom.y"
#colnames(new_miso)[16] <- "chrom.x"
miso_suppa <- merge(new_miso, SUPPA_A5_omit, by = "event_name_suppa", all=FALSE)
miso_suppa_sig <- subset(miso_suppa, miso_suppa$bayes_factor > 10)
miso_suppa_sig <- subset(miso_suppa_sig, miso_suppa_sig$p_value < 0.05)



################ MATS SUPPA
mats_suppa <- cbind(MATS_A5, as.data.frame(paste0(MATS_A5$chr,":",MATS_A5$longExonEnd,"-",MATS_A5$flankingES,":",MATS_A5$shortEE,"-",MATS_A5$flankingES,":",MATS_A5$strand)))
colnames(mats_suppa)[25] <- "event_name_suppa"
mats_suppa <- merge(mats_suppa, SUPPA_A5_omit, by = "event_name_suppa", all=FALSE)
mats_suppa_sig <- subset(mats_suppa, mats_suppa$FDR < 0.05)
mats_suppa_sig <- subset(mats_suppa_sig, mats_suppa_sig$p_value < 0.05)

#significant hits
miso_a5_sig <- subset(MISO_A5, MISO_A5$bayes_factor > 10)
mats_a5_sig <- subset(MATS_A5, MATS_A5$FDR < 0.05)
suppa_a5_sig <- subset(SUPPA_A5_omit, SUPPA_A5_omit$p_value < 0.05)

######STATISTICS
sink('output_A5.txt')

#print events
print("MATS A5")
dim(MATS_A5)
dim(mats_a5_sig)
print("MISO A5")
dim(MISO_A5)
dim(miso_a5_sig)
print("SUPPA A5")
dim(SUPPA_A5)
dim(suppa_a5_sig)
print("MISO MATS")
dim(miso_mats)
dim(miso_mats_sig)
print("MISO SUPPA")
dim(miso_suppa)
dim(miso_suppa_sig)
print("MATS SUPPA")
dim(mats_suppa)
dim(mats_suppa_sig)

print("CORRELATION MISO MATS")
lm_miso_mats <-lm(miso_mats$IncLevelDifference~miso_mats$diff, data = miso_mats)
summary(lm_miso_mats)

lm_miso_mats_sig <-lm(miso_mats_sig$IncLevelDifference~miso_mats_sig$diff, data = miso_mats_sig)
summary(lm_miso_mats_sig)

print("CORRELATION MISO SUPPA")
lm_miso_suppa <-lm(miso_suppa$diff~miso_suppa$dPSI, data = miso_suppa)
summary(lm_miso_suppa)

lm_miso_suppa_sig <-lm(miso_suppa_sig$diff~miso_suppa_sig$dPSI, data = miso_suppa_sig)
summary(lm_miso_suppa_sig)

print("CORRELATION MATS SUPPA")
lm_mats_suppa <-lm(mats_suppa$dPSI~mats_suppa$IncLevelDifference, data = mats_suppa)
summary(lm_mats_suppa)

lm_mats_suppa_sig <-lm(mats_suppa_sig$dPSI~mats_suppa_sig$IncLevelDifference, data = mats_suppa_sig)
summary(lm_mats_suppa_sig)

sink()

rm(list = ls())

#unfiltered
MATS_A3 <- read.delim("A3SS.MATS.JCEC.txt", header=TRUE)
MISO_A3 <- read.delim("CEMWT_vs_R30DM_A3.miso_bf")
SUPPA_A3 <- read.delim("CEMWT_r30dm_diffSplice.dpsi.temp.0", header=TRUE)

##MATS
##change MATS coords
MATS_A3$longExonStart_0base <- as.numeric(as.character(MATS_A3[,6])) + 1
MATS_A3$flankingES <- as.numeric(as.character(MATS_A3[,10])) + 1

#create event for MATS that matches MISO
MATS_A3 <- cbind(MATS_A3, as.data.frame(paste0 (MATS_A3$chr,":",MATS_A3$flankingES, ":", MATS_A3$flankingEE,":-@", MATS_A3$chr, ":", MATS_A3$shortEE, "|", MATS_A3$longExonEnd,":", MATS_A3$longExonStart_0base, ":",MATS_A3$strand)))
colnames(MATS_A3)[24] <- "event_name" 

#merge on event_name
miso_mats <- merge(MATS_A3, MISO_A3, by = "event_name", all=FALSE)
miso_mats_sig <- subset(miso_mats, miso_mats$bayes_factor > 10)
miso_mats_sig <- subset(miso_mats_sig, miso_mats_sig$PValue < 0.05)

##MISO
#split coordinte by ":" and transpose
miso_split <- t(as.data.frame(strsplit(as.character(MISO_A3$event_name), "-|\\+|\\:")))

#change column names and bind to original MISO file
colnames(miso_split) <-as.character(c("chrom", "riExonStart_0base", "upstreamEE","V4", "@chr", "downstreamES", "downstreamEE", "V7"))
miso_suppa <- cbind(MISO_A3, miso_split)
colnames(miso_suppa) <- make.unique(names(miso_suppa))
miso_suppa <- separate(data = miso_suppa, col= "downstreamES", into = c("downstreamES.x", "downstreamES"), sep="\\|", fill = "right")

##SUPPA
#import table with SUPPA2 dpsi values and change column names
colnames(SUPPA_A3) <- c("event", "dPSI","p_value")

#split columns and separate Ensembl ID 
SUPPA_A3 <- cbind.data.frame(str_split_fixed(SUPPA_A3$event, ';', 2), SUPPA_A3$dPSI, SUPPA_A3$p_value)
colnames(SUPPA_A3) <- c("event", "coord", "dPSI","p_value")

#split columns and separate coordinate)
SUPPA_A3 <- cbind.data.frame(SUPPA_A3$event, str_split_fixed(SUPPA_A3$coord, ':', 2), SUPPA_A3$dPSI, SUPPA_A3$p_value)
colnames(SUPPA_A3) <- c("ensembl_gene_id", "event_type", "event_name_suppa", "dPSI","p_value")

#Remove all rows with NA
SUPPA_A3_omit <- na.omit(SUPPA_A3)
SUPPA_A3_omit <- SUPPA_A3_omit[grep("A3", SUPPA_A3_omit$event_type),]

#put chr before event    
SUPPA_A3_omit$event_name_suppa <- paste0("chr", SUPPA_A3_omit$event_name_suppa)

#Get Gene list from BioMart and match genesymbols with ensembl IDs
#library('biomaRt')
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#gene_list_mart <- getBM(filters= "ensembl_gene_id",attributes= c("hgnc_symbol", "ensembl_gene_id"),values= SUPPA_A3$ensembl_gene_id, mart= mart, uniqueRows = TRUE)
#SUPPA_A3 <- merge(dpsi,gene_list_mart)

#colnames(miso_suppa)[15] <- "chrom.x"

## SUPPA MISO ## merge on event_name
#remove first column name and create suppa event coords A3
rownames(miso_suppa) <- NULL
miso_suppa <- cbind(miso_suppa, as.data.frame(paste0(miso_suppa$chrom,":",miso_suppa$downstreamES, "-", miso_suppa$riExonStart_0base,":",miso_suppa$downstreamES.x,"-",miso_suppa$riExonStart_0base,":", miso_suppa$strand)))
colnames(miso_suppa)[28] <- "event_name_suppa"


#colnames(miso_suppa)[15] <- "chrom.y"
#colnames(miso_suppa)[16] <- "chrom.x"
miso_suppa <- merge(miso_suppa, SUPPA_A3_omit, by = "event_name_suppa", all=FALSE)
miso_suppa_sig <- subset(miso_suppa, miso_suppa$bayes_factor > 10)
miso_suppa_sig <- subset(miso_suppa_sig, miso_suppa_sig$p_value < 0.05)



################ MATS SUPPA
mats_suppa <- cbind(MATS_A3, as.data.frame(paste0(MATS_A3$chr,":",MATS_A3$longExonEnd,"-",MATS_A3$flankingES,":",MATS_A3$shortEE,"-",MATS_A3$flankingES,":",MATS_A3$strand)))
colnames(mats_suppa)[25] <- "event_name_suppa"
mats_suppa <- merge(mats_suppa, SUPPA_A3_omit, by = "event_name_suppa", all=FALSE)
mats_suppa_sig <- subset(mats_suppa, mats_suppa$FDR < 0.05)
mats_suppa_sig <- subset(mats_suppa_sig, mats_suppa_sig$p_value < 0.05)

#significant hits
miso_a3_sig <- subset(MISO_A3, MISO_A3$bayes_factor > 10)
mats_a3_sig <- subset(MATS_A3, MATS_A3$FDR < 0.05)
suppa_a3_sig <- subset(SUPPA_A3_omit, SUPPA_A3_omit$p_value < 0.05)

######STATISTICS
sink('output_A3.txt')

#print events
print("MATS A3")
dim(MATS_A3)
dim(mats_a3_sig)
print("MISO A3")
dim(MISO_A3)
dim(miso_a3_sig)
print("SUPPA A3")
dim(SUPPA_A3)
dim(suppa_a3_sig)
print("MISO MATS")
dim(miso_mats)
dim(miso_mats_sig)
print("MISO SUPPA")
dim(miso_suppa)
dim(miso_suppa_sig)
print("MATS SUPPA")
dim(mats_suppa)
dim(mats_suppa_sig)

print("CORRELATION MISO MATS")
lm_miso_mats <-lm(miso_mats$IncLevelDifference~miso_mats$diff, data = miso_mats)
summary(lm_miso_mats)

lm_miso_mats_sig <-lm(miso_mats_sig$IncLevelDifference~miso_mats_sig$diff, data = miso_mats_sig)
summary(lm_miso_mats_sig)

print("CORRELATION MISO SUPPA")
lm_miso_suppa <-lm(miso_suppa$diff~miso_suppa$dPSI, data = miso_suppa)
summary(lm_miso_suppa)

lm_miso_suppa_sig <-lm(miso_suppa_sig$diff~miso_suppa_sig$dPSI, data = miso_suppa_sig)
summary(lm_miso_suppa_sig)

print("CORRELATION MATS SUPPA")
lm_mats_suppa <-lm(mats_suppa$dPSI~mats_suppa$IncLevelDifference, data = mats_suppa)
summary(lm_mats_suppa)

lm_mats_suppa_sig <-lm(mats_suppa_sig$dPSI~mats_suppa_sig$IncLevelDifference, data = mats_suppa_sig)
summary(lm_mats_suppa_sig)

sink()

rm(list = ls())

SE_output <- read.csv("output_SE.txt", sep="", header=FALSE)
SE_output_dim <- SE_output[c(2,3,5,6,8,9,11,12,14,15,17,18),]
SE_output_dim$V3 <- NULL
rownames(SE_output_dim) <- c("All MATS events", "Significant MATS events", "All MISO Events", "significant MISO events", "All SUPPA events", "Significant SUPPA events", "MISO MATS overlap", "MISO MATS sig", "MISO SUPPA overlap", "MISO SUPPA sig", "MATS SUPPA overlap", "MATS SUPPA sig")
colnames(SE_output_dim) <- c("Events")
rownames(SE_output) <- make.names(names=SE_output$V1, unique=TRUE)
SE_output_Rsquare <- SE_output[c("Adjusted.1","Adjusted.3","Adjusted.5"),]
SE_output_Rsquare$V1 <- NULL
SE_output_Rsquare$V2 <- NULL 
rownames(SE_output_Rsquare) <- c("Adj. Rsq. MATS MISO", "Adj. Rsq MISO SUPPA", "Adj. Rsq. MATS SUPPA")

RI_output <- read.csv("output_RI.txt", sep="", header=FALSE)
RI_output_dim <- RI_output[c(2,3,5,6,8,9,11,12,14,15,17,18),]
RI_output_dim$V3 <- NULL
rownames(RI_output_dim) <- c("All MATS events", "Significant MATS events", "All MISO Events", "significant MISO events", "All SUPPA events", "Significant SUPPA events", "MISO MATS overlap", "MISO MATS sig", "MISO SUPPA overlap", "MISO SUPPA sig", "MATS SUPPA overlap", "MATS SUPPA sig")
colnames(RI_output_dim) <- c("Events")
rownames(RI_output) <- make.names(names=RI_output$V1, unique=TRUE)
RI_output_Rsquare <- RI_output[c("Adjusted.1","Adjusted.3","Adjusted.5"),]                    
RI_output_Rsquare$V1 <- NULL
RI_output_Rsquare$V2 <- NULL 
rownames(RI_output_Rsquare) <- c("Adj. Rsq. MATS MISO", "Adj. Rsq MISO SUPPA", "Adj. Rsq. MATS SUPPA")

A5_output <- read.csv("output_A5.txt", sep ="", header=FALSE)
A5_output_dim <- A5_output[c(2,3,5,6,8,9,11,12,14,15,17,18),]
A5_output_dim$V3 <- NULL
rownames(A5_output_dim) <- c("All MATS events", "Significant MATS events", "All MISO Events", "significant MISO events", "All SUPPA events", "Significant SUPPA events", "MISO MATS overlap", "MISO MATS sig", "MISO SUPPA overlap", "MISO SUPPA sig", "MATS SUPPA overlap", "MATS SUPPA sig")
colnames(A5_output_dim) <- c("Events")
rownames(A5_output) <- make.names(names=A5_output$V1, unique=TRUE)
A5_output_Rsquare <- A5_output[c("Adjusted.1","Adjusted.3","Adjusted.5"),]
A5_output_Rsquare$V1 <- NULL
A5_output_Rsquare$V2 <- NULL 
rownames(A5_output_Rsquare) <- c("Adj. Rsq. MATS MISO", "Adj. Rsq MISO SUPPA", "Adj. Rsq. MATS SUPPA")

A3_output <- read.csv("output_A3.txt", sep="", header=FALSE)
A3_output_dim <- A3_output[c(2,3,5,6,8,9,11,12,14,15,17,18),]
A3_output_dim$V3 <- NULL
rownames(A3_output_dim) <- c("All MATS events", "Significant MATS events", "All MISO Events", "significant MISO events", "All SUPPA events", "Significant SUPPA events", "MISO MATS overlap", "MISO MATS sig", "MISO SUPPA overlap", "MISO SUPPA sig", "MATS SUPPA overlap", "MATS SUPPA sig")
colnames(A3_output_dim) <- c("Events")
rownames(A3_output) <- make.names(names=A3_output$V1, unique=TRUE)
A3_output_Rsquare <- A3_output[c("Adjusted.1","Adjusted.3","Adjusted.5"),]
A3_output_Rsquare$V1 <- NULL
A3_output_Rsquare$V2 <- NULL 
rownames(A3_output_Rsquare) <- c("Adj. Rsq. MATS MISO", "Adj. Rsq MISO SUPPA", "Adj. Rsq. MATS SUPPA")

sink(file= "clean_output.csv")

print("SE")
SE_output_dim
SE_output_Rsquare

print("RI")
RI_output_dim
RI_output_Rsquare

print("A5")
A5_output_dim
A5_output_Rsquare

print("A3")
A3_output_dim
A3_output_Rsquare

sink()