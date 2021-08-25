#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

COUNTS_FILE=args[1] # "poscon_lib1_matrix.txt"
DRUG = args[2]

data = read.table(COUNTS_FILE,sep='\t',head=T)
data = data[data$conc <= 1,] # remove concs above 1uM for TMP and MTX in Broad data

sample_info = data[,c("drug","conc")]
colnames(sample_info) = c("Drug","Conc")
sample_info$ConcLabl = sprintf("%guM",data$conc)

sample_info$logConc[data$conc!=0] = log(data[data$conc!=0,"conc"],2) # calc logConc for non-zeros

# for each drug, set 0uM conc to half of lowest non-zero conc

for (drug in unique(sample_info$Drug))
{
  m = min(data[data$drug==drug & data$conc!=0,"conc"])
  NO_DRUG_LOG_CONC = log(m/2,2)
  sample_info$logConc[data$drug==drug & data$conc==0] = NO_DRUG_LOG_CONC
  cat(sprintf("%s, min_conc=%s (log2=%s), NO_DRUG_LOG_CONC=%s (log2=%s)\n",drug,m,round(log(m,2),5),m/2,round(NO_DRUG_LOG_CONC,5)))
}

counts = data[,8:162] # 528x155 in Broad data, library1, 4 drugs

###########################

source("CGA_functions.R")

melted = get_filtered_and_melted_data(DRUG,sample_info,counts,30) # lower MIN_COUNT_PER_CELL

write.table(melted,sprintf("%s_melted.txt",DRUG),sep='\t',row.names=F,quote=F)

