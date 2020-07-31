

ibrary(dplyr)
library(stringr)
#command line arguments.
#First command is MATS, second is MISO, third is SUPPA2 file 
args = commandArgs(trailingOnly = TRUE)
df <- read.delim(args[1], header=TRUE, row.names=NULL)

#df <- read.delim("D2_100m_2v2_MATS.log", row.names=NULL)

df <- df[,c("row.names","X","Time")]
colnames(df) <- c("Mem","CPU","TIME")
df$CPU <- str_remove(df$CPU, "[%]")

max(df$Mem)
max(df$CPU)

