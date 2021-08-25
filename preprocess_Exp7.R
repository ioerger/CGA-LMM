#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

COUNTS_FILE=args[1] # "Exp7_counts.txt"
DRUG = args[2]

data = read.table(COUNTS_FILE,sep='\t',head=T)
data = data[data$note=='sample',] # keep only rows with note==sample

sample_info = data[,c("Drug","xMIC")]
colnames(sample_info) = c("Drug","Conc")
sample_info$ConcLabl = sample_info$Conc

NO_DRUG_CONC = '0.000xMIC'
NO_DRUG_LOG_CONC = -4 # Exp7, lowest / 2
sample_info[sample_info$Drug=='Input','Conc'] = NO_DRUG_CONC
sample_info$logConc[sample_info$Conc==NO_DRUG_CONC] = NO_DRUG_LOG_CONC

sample_info$logConc[sample_info$Conc=="8xMIC"] = 3
sample_info$logConc[sample_info$Conc=="4xMIC"] = 2
sample_info$logConc[sample_info$Conc=="2xMIC"] = 1
sample_info$logConc[sample_info$Conc=="1xMIC"] = 0
sample_info$logConc[sample_info$Conc=="0.5xMIC"] = -1
sample_info$logConc[sample_info$Conc=="0.25xMIC"] = -2
sample_info$logConc[sample_info$Conc=="0.125xMIC"] = -3

counts = data[,20:190] # exclude controls at the end of each, row starting with H37Rv
counts = counts[sample_info$logConc <= 0,] # keep MICs only up to 1xMIC (486 samples)
sample_info = sample_info[sample_info$logConc <= 0,] 


###########################

source("CGA_functions.R")

melted = get_filtered_and_melted_data(DRUG,sample_info,counts)

write.table(melted,sprintf("%s_melted.txt",DRUG),sep='\t',row.names=F,quote=F)
