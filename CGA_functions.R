##################################
# melt data, and filter out genes or drugs with too low or high abundance

# sample_info and counts are parallel data frames; 
#   sample_info has columns: Drug, Conc, ConcLabl, logConc, and other stuff...
#   counts is a matrix of counts for sample (row) and gene (col)
# this does filtering of genes and/or drugs with too little counts or too high abundance
# the last 2 args are optional, default is 100 and 0.1, but can change it

get_filtered_and_melted_data = function(DRUG,sample_info,counts,MIN_COUNT_PER_CELL=100,MAX_ABUNDANCE=0.1)
{
  # drop genes (columns) with less than 100 counts per cell on average

  coltot = colSums(counts)
  print("total counts for each gene over all samples:")
  print(sort(coltot))
  cat("dropping these columns because not enough counts:\n")
  print(colnames(counts[,coltot<MIN_COUNT_PER_CELL*nrow(counts)]))
  counts1 = counts[,coltot>=MIN_COUNT_PER_CELL*nrow(counts)] 
  
  # filtering rows and cols with low counts
  # drop samples with rowtot<20000 (45 samples), or less than 100 per cell on average
  
  rowtot = rowSums(counts1)
  cat("dropping these samples (rows) because not enough counts:\n")
  print(sample_info[rowtot<MIN_COUNT_PER_CELL*ncol(counts1),c("Drug","ConcLabl")])
  data2 = sample_info[rowtot>=MIN_COUNT_PER_CELL*ncol(counts1),]
  
  print("filtered:"); print(dim(counts1)) # this is after lower-bound filters, but filtering high-abundance genes is below...

  # compute relative abundances
  
  relAbund1 = sweep(counts1,1,rowtot,'/')
  relAbund2 = relAbund1[rowtot>=MIN_COUNT_PER_CELL*ncol(counts1),] 
  # note: data2 and relAbund2 are parallel arrays, with metadata in data2 and abundance for each drug in relAbund2
  
  require(reshape)
  
  drug_data = relAbund2[data2$Drug==DRUG,]
  drug_data$ConcLabl = data2[data2$Drug==DRUG,"ConcLabl"]
  drug_data$logConc = data2[data2$Drug==DRUG,"logConc"]
  melted = melt(drug_data,id=c("ConcLabl","logConc") )
  colnames(melted) = c("ConcLabl","logConc","gene","Abund")
  
  # filter out genes with >10% abundance at any conc (avg'd over replicates) for a given drug
  
  agg1 = aggregate(Abund~gene+logConc,data=melted,FUN=mean)
  agg2 = aggregate(Abund~gene,data=agg1,FUN=max)
  print("max relative abundance at any concentration:")
  print(agg2[order(agg2$Abund),])

  for (gene in unique(agg2$gene))
  {
    if (agg2$Abund[agg2$gene==gene]>MAX_ABUNDANCE)
    {
      cat(sprintf("eliminating %s due to excessive abundance\n",gene))
      print(agg1[agg1$gene==gene,])
      melted = melted[melted$gene!=gene,] 
    }
  }
  print("melted:"); print(dim(melted))

  return(melted)
}

####################################
# run linear mixed model

library(reshape)
library(lme4)

run_LMM = function(model_data)
{
  mod1 = lmer(Abund~1+logConc+(1|gene)+(0+logConc|gene),data=model_data)
  print(summary(mod1))
  predicted = as.data.frame(predict(mod1,model_data))

  coeffs = ranef(mod1)$gene
  colnames(coeffs)[2] = 'slope'


  # calculate significance of outlier slopes as Zrobust

  allSlopes = coeffs$slope
  med = median(allSlopes)
  MAD = median(abs(allSlopes-med))
  Zrobust = 0.6745*(allSlopes-med)/MAD


  # do simple linear regressions to test significance of coefficient for each gene

  lm_results = NULL
  for (gene in rownames(coeffs))
  {
    temp = melted[melted$gene==gene,]
    lmod = lm(Abund ~ logConc,data=temp)
    lmcoeffs = summary(lmod)$coefficients
    lmcoeffs = as.data.frame(lmcoeffs)
    lmslope = lmcoeffs["logConc",1]
    pval = lmcoeffs["logConc",4]
    #cat(sprintf("%s lm_slope=%s, pval=%s\n",gene,lmslope,pval))
    lm_results = rbind(lm_results,c(gene=gene,lm_slope=lmslope,pval=pval))
  }
  lm_results = as.data.frame(lm_results)
  lm_results$lm_slope = as.numeric(as.character(lm_results$lm_slope))
  lm_results$pval = as.numeric(as.character(lm_results$pval))
  lm_results$padj = p.adjust(lm_results$pval,method="BH")
  cat(sprintf("num of genes with signif neg slope in lm: %s (pval<0.05), %s (padj<0.05)\n",sum(lm_results$pval<0.05 & lm_results$lm_slope<0),sum(lm_results$padj<0.05 & lm_results$lm_slope<0)))


  # write out results

  #res = cbind(LM_slope=lm_results$lm_slope,Padj=lm_results$padj,meanAbundNoDrug=ctrlMean$Abund,LMM_intercept=coeffs[,1],LMM_slope=coeffs[,2],Zrobust)
  res = cbind(LM_slope=lm_results$lm_slope,Padj=lm_results$padj,LMM_intercept=coeffs[,1],LMM_slope=coeffs[,2],Zrobust)
  res = as.data.frame(res)
  rownames(res) = rownames(coeffs)
  res = res[order(res$Zrobust),] 

  write.table(res,sprintf("%s_coeffs.txt",DRUG),sep='\t',quote=F)

  return(list(results=res,model=mod1))
}

##############################
# plot histogram, with significance cutoffs

# allSlopes is a data frame with Genes as rownames and LMM_slope as column
# GenesToHighlight is a list of gene names; can be NULL

slope_histogram = function(DRUG,allSlopes,GenesToHighlight)
{
  fname = sprintf("%s_hist.png",DRUG)
  cat(sprintf("generating %s\n",fname))
  png(fname)

  slopes = allSlopes$LMM_slope
  y1 = 20; y2=10
  med = median(slopes)
  MAD = median(abs(slopes-med))
  cutoff = 3.5*MAD/0.6745
  a=med-cutoff; b=med+cutoff
  CEX = 1.5
  hist(slopes,breaks=30,main=sprintf("Slope Distribution for %s",DRUG),cex=CEX,cex.lab=CEX,cex.axis=CEX)
  abline(v=med,lty=2)
  abline(v=med-cutoff,lty=2,col=2)
  abline(v=med+cutoff,lty=2,col=2)

  count = 0
  for (Gene in GenesToHighlight)
  {
    x = allSlopes[Gene,1]
    arrows(x,y1+5*count,x,y2,col=4,length=0.1)
    text(x,y1+1+5*count,label=Gene,cex=CEX)
    count = count+1
  }
  dev.off()
}

#######################################3
# make fan plot

# extract a mapping between logConcs (numeric) and Concs (as text labels)

get_concs = function(model_data)
{
  temp = melted[,c("ConcLabl","logConc")]
  mapping = unique(temp) 
  mapping = mapping[order(mapping$logConc),]
  return(mapping)
}

make_fan_plot = function(Drug,Melted,Genes,GenesToHighlight)
{
  mapping = get_concs(model_data)
  LOG_CONC_RANGE = mapping$logConc
  CONC_LABELS = mapping$ConcLabl
  NO_DRUG_LOG_CONC = min(LOG_CONC_RANGE)

  fname = sprintf("%s_fan_plot.png",DRUG)
  cat(sprintf("generating %s\n",fname))
  png(fname)
  first = T
  for (Gene in c(Genes,GenesToHighlight))
  {
   if (!(Gene %in% Melted$Gene)) { next }
   cat(sprintf("fan: %s\n",Gene))
   subset = Melted[Melted$gene==Gene,]
   pred = predicted[Melted$gene==Gene,]
   X = subset$logConc
   Y = subset$Abund

   pred0 = mean(pred[subset$logConc==NO_DRUG_LOG_CONC])
   m = mean(subset$logConc)
   m0 = mean(Y[subset$logConc==NO_DRUG_LOG_CONC])

   col = 1 
   if (Gene %in% GenesToHighlight) { col = 2 }

   YLIM = c(-2,2)
   if (first==T) 
   {
     plot(X,100*(pred-pred0),main=sprintf("drug=%s",Drug),xlab="",ylab="normalized abundance",ylim=YLIM,xaxt="n",type='l',col=col)
     axis(1, at=LOG_CONC_RANGE, labels=CONC_LABELS,las=2)   
     first = F
   } else 
   { 
     plot(X,100*(pred-pred0),xlab="",ylab="",main="",ylim=YLIM,xaxt="n",yaxt="no",type='l',col=col) # intersect at 1 for conc=0
     if (Gene %in% GenesToHighlight) { a = max(X); b = pred[X==a]; text(a,100*(b-pred0),label=Gene,col=col) } 
   } 
   par(new=T)
 }
 dev.off()
}

########################
# generate dot plots

# dot plot for individual gene
# model_data is melted

dot_plot_gene = function(Drug,Gene,model_data,predicted,LOG_CONC_RANGE,CONC_LABELS)
{
  NO_DRUG_LOG_CONC = min(LOG_CONC_RANGE)

  subset = model_data[model_data$gene==Gene,]
  pred = predicted[model_data$gene==Gene,]
  X = subset$logConc
  Y = subset$Abund
  m = mean(Y[subset$logConc==NO_DRUG_LOG_CONC])

  plot(X,Y/m,main=sprintf("drug=%s ; gene=%s",Drug,Gene),xlab="",ylab="normalized abundance",ylim=c(0,5),xaxt="n")
  axis(1, at=LOG_CONC_RANGE, labels=CONC_LABELS, las=2)
  par(new=T)
  lines(X,pred/m,col=2) 
}

# this assumes model_data is melted and already been selected for DRUG (just used for filename and plot label)

plot_genes_PDF = function(DRUG,Genes,model_data,predicted)
{
  mapping = get_concs(model_data)
  LOG_CONC_RANGE = mapping$logConc
  CONC_LABELS = mapping$ConcLabl

  fname = sprintf("%s_dot_plots.pdf",DRUG)
  cat(sprintf("generating %s\n",fname))
  pdf(fname)
  par(mfrow=c(3,2))
  for (gene in Genes) { dot_plot_gene(DRUG,gene,model_data,predicted,LOG_CONC_RANGE,CONC_LABELS) }
  dev.off()
}

# make combined dot plot as .png

combined_dot_plot = function(DRUG,model_data,Genes,GenesToHighlight)
{
  mapping = get_concs(model_data)
  LOG_CONC_RANGE = mapping$logConc
  CONC_LABELS = mapping$ConcLabl
  NO_DRUG_LOG_CONC = min(LOG_CONC_RANGE)

  fname = sprintf("%s_dot_plots_combined.png",DRUG)
  cat(sprintf("generating %s\n",fname))
  png(fname)
  par(mar=c(7,4.1,4.1,2.1))
  firstpass = T
  for (Gene in c(Genes,GenesToHighlight))
  {
    col = "gray"
    if (!(Gene %in% model_data$gene)) { next }
    if (Gene %in% GenesToHighlight) { col = 2}
    subset = model_data[model_data$gene==Gene,]
    means = aggregate(Abund ~ logConc,data=subset,FUN=mean)
    m0 = means[means$logConc==NO_DRUG_LOG_CONC,"Abund"]
    X = means$logConc
    Y = means$Abund

    CEX = 1.5
    YLIM = c(-2,2)
    if (firstpass==T) { 
      firstpass = F
      plot(X,100*(Y-m0),main=sprintf("drug=%s",DRUG),xlab="",ylab="change in percent abundance",ylim=YLIM,xaxt="n",col=col,pch=".",cex=CEX,cex.lab=CEX,cex.axis=CEX)
      axis(1, at=LOG_CONC_RANGE, labels=CONC_LABELS,las=2,cex.axis=CEX) 
    }
    else { plot(X,100*(Y-m0),ylim=YLIM,xaxt="n",yaxt="no",xlab="",ylab="",col=col,pch=".") }
    lines(X,100*(Y-m0),col=col)
    if (Gene %in% GenesToHighlight) { a = max(X); b = Y[X==a]; text(a-0.05,100*(b-m0),label=Gene,cex=CEX) }
    par(new=T)
  }
  dev.off()
}

