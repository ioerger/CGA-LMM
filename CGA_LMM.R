#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

meltedFilename = args[1]
DRUG = args[2]

# melted should have columns: logConc, Conc, gene, Abund

melted = read.table(meltedFilename,sep='\t',head=T)

####################################
# run linear mixed model

source("CGA_functions.R")

model_data = melted

tuple = run_LMM(model_data)
results = tuple$results # data frame with slope, Zrobust, etc for each gene
model = tuple$model

allSlopes = results["LMM_slope"] # for histogram
predicted = as.data.frame(predict(model,model_data)) # predicted relAbund for every data point in melted, for regression lines in dot plots

################################
# generate plots

Genes = unique(model_data$gene)
Genes = sort(as.character(Genes))

GenesToHighlight = NULL
if (DRUG=="Sulfa") { GenesToHighlight=c("thyA.Rv2764c") }
if (DRUG=="INH") { GenesToHighlight=c("ino1.Rv0046c","kasB.Rv2246") } # removed "fas.Rv2524c"
if (DRUG=="Levo") { GenesToHighlight=c("gyrA.Rv0006") }
if (DRUG=="Moxi") { GenesToHighlight=c("gyrA.Rv0006") }
if (DRUG=="Fida") { GenesToHighlight=c("rpoB.Rv0667"); }
if (DRUG=="BDQ") { GenesToHighlight=c("atpB.Rv1304","atpH.Rv1307","atpF.Rv1306","atpG.Rv1309"); }
if (DRUG=="trimethoprim" | DRUG=="methotrexate") { GenesToHighlight = c("trpG"); }
if (DRUG=="rifampin") { GenesToHighlight = c("rpoB"); }
if (substr(DRUG,1,3)=="Cu_") { GenesToHighlight = c("Rv2158c.murE","Rv1315.murA","Rv2157c.murF","Rv2155c.murD","Rv2152c.murC","Rv2156c.mraY") }


slope_histogram(DRUG,allSlopes,GenesToHighlight)

make_fan_plot(DRUG,model_data,Genes,GenesToHighlight)

plot_genes_PDF(DRUG,Genes,model_data,predicted)

combined_dot_plot(DRUG,model_data,Genes,GenesToHighlight)
