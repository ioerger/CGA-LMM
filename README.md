CGA_LMM: Chemical-Genomics Analysis with Linear Mixed Models
============================================================

This github repository holds scripts for statistical analysis of
chemical-genomics data using linear mixed models, as described in:


>  Esha Dutta, Michael A. DeJesus, Nadine Ruecker, Anisha Zaveri, Eun-Ik Koh, Christopher M. Sassetti, Dirk Schnappinger, and Thomas R. Ioerger (2021).  "An Improved Statistical Method to Identify Chemical-Genetic Interactions by Exploiting Concentration-Dependence", submitted to PLOS ONE.


Overview of the code
--------------------

The scripts are written in R, and were developed and tested on v3.6.3 of R.

Some R packages are required to be installed: lme4, reshape

The two main files to be used for analyzing any dataset are:
CGA_functions.R and CGA_LMM.R.  CGA_functions contains utility
functions to be used in the analysis, and is *sourced* by the other
scripts.  CGA_LMM.R is called to perform the LMM analysis and generate
some plots and output files.  In addition, there is a custom
preprocessing script for each of the 3 example datasets.



Workflow
--------

In a typical C-G experiment, a hypomorph library of bacterial
depletion mutants will be treated with a growth inhibitor, the DNA
will be extracted and sequenced, and nucleotide barcode counts for
each strain will be tabulated and saved in a matrix (spreadsheet,
saved as a tab-separated text file) where the rows are samples
(treatments with different drugs at different concentrations, columns
are genes, and the counts represent abundance of each mutant in each
sample.

Step 1. Preprocessing

Preprocess the input file of 'counts' to generated a file with 'melted' data.
This step can be customized for each dataset, to accommodate anything
unique about the input file (such as extra columns or rows that can be discarded,
or column headers that have to be renamed, etc.)
This step also performs filtering (to remove samples or drugs with too few counts,
or drugs with excessive abundance).
It calculates relative abundance (normalizing by total counts for each sample),
and also calculated the log2 of drug concentrations.  Finally it melts the 
data into a column format suitable for input to the next step.
All of these operations are performed in a standardized way by calling the
get_filtered_and_metled_data() function definde in CGA_functions.R
See any of the preprocess_DATASET.R scripts as examples.

> usage: Rscript preprocess_DATASET.txt DRUG


Step 2. Run the LMM and generate output files and plots.  

This is done by running the CGA_LMM.R, which reads the
DRUG_melted.txt file as input and generates DRUG_coeffs.txt file
as output.  For interpretation, see Output section below.

> Rscript CGA_LMM.R DRUG_melted.txt DRUG


Example 
-------

Here is an example of running the analysis on levofloxcin in the Exp7 dataset:

1. preprocessing input counts file to generate melted data

  > Rscript preprocess_Exp7.R Exp7_counts.txt Levo

  generates Levo_melted.txt

2. running the analysis and generating output files and plots

  > Rscript CGA_LMM.R Levo_melted.txt Levo

  * reads Levo_melted.txt
  * generates: 
    * Levo_coeffs.txt - a spreadsheet with analysis of each gene (slope, significance, Zrobust)
    * Levo_hist.png - histogram of distribution of slopes
    * Levo_dot_plots.pdf - plot showing the raw data for abundance of each gene and the regression line fit by the LMM
    * Levo_dot_plot_combined.png - plot showing the mean abundance for each gene superimposed
    * Levo_fan_plot.png - plot showing the data reduced to linear trends for each gene


Interpetation of Output (coeffs.txt file)
------------------------

The main output file is DRUG_coeffs.txt, which is a tab-separated
text file that can be opened as a spreadsheet.  The primary column to
focus on is Zrobust, which gives the robust Z-score based on slope
(trend of relative abundance with increasing drug concentration) for
each gene.  Genes that are outliers or candidate interactions are
those with |Zrobust|>3.5.  However, the most interesting interactions
are usually the ones with negative interactions, Zrobust<-3.5.

Furthermore, to assess which genes have a regression sloped that is
significantly different than 0, a regular linear model was fit for
each gene, and the slope coefficient was tested for significance
using a t-statistic; an adjusted P-value is derived from this.
Note that many genes in the test datasets have negative slopes that
are significantly different than 0, but there are usually only a
small number that are outliers, and hence identified as candidate interactions.

The output file has the following columns:

  * LM_slope:      regression coefficient of relative abundance against log2 of drug concentration
  * Padj:          P-value for where coefficient is significantly different than 0, adjusted for multiple comparisons by the Benjamini-Hochberg procedure; genes with slopes with Padj>=0.05 could be disregarded
  * LMM_intercept: random effect intercept for each gene in LMM
  * LMM_slope:     random effect slope for each gene in LMM: trend of relative abundance against log2 of drug concentration
  * Zrobust:       robust Z-score (Iglewicz and Hoaglin, 1993) of LMM_slope with respect to distribution of slopes for all genes; outliers (interations with drug) are those genes with |Zrobust|>3.5

> Iglewicz, B., &amp; Hoaglin, D. (1993). How to Detect and Handle Outliers (E. F. Mykytka Ed. Vol. 16): ASCQ Quality Press.


Functions in CGA_functions.R (API)
----------------------------

* **get_filtered_and_melted_data**(DRUG,sample_info,counts,MIN_COUNT_PER_CELL=100,MAX_ABUNDANCE=0.1):

  * This function takes 2 data frames as input (and the name of a DRUG to analyze; the last 2 args are optional and are used in filtering).

  * sample_info has columns: Drug, Conc, ConcLabl, logConc (multiple replicates are stored on separate rows) (ConcLabl contains text symbols for labeling X-axis of plots)
  * counts is a matrix of counts for each sample (row) and gene (col); it is intended to be parallel to sample_info (i.e. corresponding rows)


* **run_LMM**(model_data):

  * The model_data is a melted data frame with ConcLabl, logConc, gene, and Abund, from the preprocessing step.
  * This runs lmer and performs statistical analysis 


The following functions are for generating the plot files.  They are fairly
self-explanatory.  If one wants to highlight certain genes in the
plot, one can pass in a list of gene names (or alternatively, leave it
as NULL).  Some examples of GenesToHighlight are defined in CGA_LMM.R


* **slope_histogram**(DRUG,allSlopes,GenesToHighlight) 

  * generates DRUG_hist.png
  * histogram of distribution of slopes

* **make_fan_plot**(Drug,Melted,Genes,GenesToHighlight)

  * generate DRUG_fan_plot.png
  * plot showing the data reduced to linear trends for each gene

* **plot_genes_PDF**(DRUG,Genes,Melted,predicted) 

  * generates DRUG_dot_plots.pdf
  * plot showing the raw data for abundance of each gene and the regression line fit by the LMM

* **combined_dot_plot**(DRUG,Melted,Genes,GenesToHighlight):

  * generates DRUG_dot_plots_combined.png
  * plot showing the mean abundance for each gene superimposed


Datasets 
-------- 

Datasets for 3 hypomorph libraries of M. tuberculosis H37Rv treated
with 11 drugs (combined) are provided in the data/ sub-directory.
Here the drugs and abbreviations tested in each library:

* Broad (input file: poscon_lib1_matrix.txt, downloaded from https://www.broadinstitute.org/chemical-biology/initiative-chemical-genetics ; library1, positive controls)

  * trimethoprim
  * methotrexate
  * rifampin
  * BRD-4592

> Rscript preprocess_Broad.R poscon_lib1_matrix.txt trimethoprim \
> Rscript CGA_LMM.R trimethoprim_melted.txt trimethoprim
>
> Rscript preprocess_Broad.R poscon_lib1_matrix.txt methotrexate \
> Rscript CGA_LMM.R methotrexate_melted.txt methotrexate
>
> Rscript preprocess_Broad.R poscon_lib1_matrix.txt rifampin \
> Rscript CGA_LMM.R rifampin_melted.txt rifampin
>
> Rscript preprocess_Broad.R poscon_lib1_matrix.txt BRD-4592 \
> Rscript CGA_LMM.R BRD-4592_melted.txt BRD-4592


* Exp7 (input file: Exp7_counts.txt)

  * Levo - levofloxacin
  * Moxi - moxifloxacin
  * INH - isoniazid
  * Fida - fidaxomycin
  * Sulfa - sulfamethoxazole
  * BDQ - bedaquiline

> Rscript preprocess_Exp7.R Exp7_counts.txt Levo \
> Rscript CGA_LMM.R Levo_melted.txt Levo
>
> Rscript preprocess_Exp7.R Exp7_counts.txt Moxi \
> Rscript CGA_LMM.R Moxi_melted.txt Moxi
>
> Rscript preprocess_Exp7.R Exp7_counts.txt INH \
> Rscript CGA_LMM.R INH_melted.txt INH
>
> Rscript preprocess_Exp7.R Exp7_counts.txt Fida \
> Rscript CGA_LMM.R Fida_melted.txt Fida
>
> Rscript preprocess_Exp7.R Exp7_counts.txt Sulfa \
> Rscript CGA_LMM.R Sulfa_melted.txt Sulfa
>
> Rscript preprocess_Exp7.R Exp7_counts.txt BDQ \
> Rscript CGA_LMM.R BDQ_melted.txt BDQ


* Cu (input file: Supplemental_Table_T3_copper_barcode_counts.txt; 
treated with copper, with 3 different carbon sources in the medium:)

  * cholesterol
  * acetate
  * glycerol

For the copper dataset, the carbon source is specified as the 
argument on the command line for the preprocessing step, instead of the drug 
(since there is only 1 drug: Cu).
The melted file it creates is called 'Cu_CARBON_SOURCE_ssb1'
The suffix 'ssb1' refers to the specific Tet promoter used for controlling sspB expression.


> Rscript preprocess_copper.R copper_barcode_counts.txt cholesterol \
> Rscript CGA_LMM.R Cu_cholesterol_ssb1_melted.txt Cu_cholesterol_ssb1
>
> Rscript preprocess_copper.R copper_barcode_counts.txt acetate \
> Rscript CGA_LMM.R Cu_acetate_ssb1_melted.txt Cu_acetate_ssb1
>
> Rscript preprocess_copper.R copper_barcode_counts.txt glycerol \
> Rscript CGA_LMM.R Cu_glycerol_ssb1_melted.txt Cu_glycerol_ssb1 



