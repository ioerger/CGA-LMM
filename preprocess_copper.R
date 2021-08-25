#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

COUNTS_FILE=args[1] "Supplemental_Table_T3_copper_barcode_counts.txt"
CARBON = args[2]
SSB = args[3]
DRUG = sprintf("Cu_%s_%s",CARBON,SSB)

data = read.table(COUNTS_FILE,sep='\t',head=T)
data[is.na(data)] = 0

# note: I changed the header of the column from "Concentration (uM)" to "Concentration_uM" in Supplemental_Table_T3_copper_barcode_counts.txt
sample_info = data[,c("Concentration_uM","Carbon_source","SSB_strength")]
colnames(sample_info) = c("Conc","Carbon","SSB")
sample_info$Drug = sprintf("Cu_%s_%s",sample_info$Carbon,sample_info$SSB)
sample_info$ConcLabl = sprintf("%suM",sample_info$Conc)


sample_info$logConc[sample_info$Conc!=0] = log(sample_info[sample_info$Conc!=0,"Conc"],2) # calc logConc for non-zeros
# for 0uM, set logConc to log of half of min conc
# in the case of copper, all drugs have conc range [1,2,4,8], so "0uM" maps to 0.5, and log("0uM",2) maps to -1 for all drugs
for (drug in unique(sample_info$Drug))
{
  m = min(sample_info[sample_info$Drug==drug & sample_info$Conc!=0,"Conc"])
  NO_DRUG_LOG_CONC = log(m/2,2)
  sample_info$logConc[sample_info$Drug==drug & sample_info$Conc==0] = NO_DRUG_LOG_CONC
  #cat(sprintf("%s, min_conc=%s (log2=%s), NO_DRUG_LOG_CONC=%s (log2=%s)\n",drug,m,round(log(m,2),5),m/2,round(NO_DRUG_LOG_CONC,5)))
}

counts = data[,6:470]

###########################

source("CGA_functions.R")

melted = get_filtered_and_melted_data(DRUG,sample_info,counts)

write.table(melted,sprintf("%s_melted.txt",DRUG),sep='\t',row.names=F,quote=F)
