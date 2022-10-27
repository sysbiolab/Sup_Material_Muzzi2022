Supporting material for Muzzi et al., 2022
================
João C. D. Muzzi, Jessica M. Magno, Jean S. S. Resende, Larissa M. Alvarenga, Juliana F. de Moura, Bonald C. Figueiredo, and Mauro A. A. Castro <br>

August 22, 2022


# Context
----
Adrenocortical carcinoma (ACC) is a rare and aggressive malignancy with few treatment options and a low survival rate. Our previous work (Muzzi et al., 2021) showed that the Steroid Phenotype defined by Zheng et al. (2016) distinguishes the immune activation in the ACC microenvironment. Here, we used The Cancer Genome Atlas (TCGA) ACC cohort to assess the regulatory network. We focused on regulons (a transcription factor and its targets) related to overall survival. We also evaluated how these regulons correlated with immune variables and molecular phenotypes. Finally, we assessed its activity and prognostic association in the ENSAT cohort. This script reproduces all results published in Muzzi et al. (2022) and serves as complementary material.
In Preprocessing.R file, we show how to obtain and preprocess the data from TCGA-ACC and ENSAT cohorts used in the study. At the end of the section, the data sets are saved as RData and can be imported for further analysis. Here we use the preprocessed data and show the steps for data analysis and how to obtain the results and plots shown in <a href="https://www.mdpi.com/2072-6694/14/21/5279/htm">Muzzi et al. (2022)<a>. 

# Loading packages
----
```r 
## For GDC Data Portal ACC RNAseq import and data manipulation
library(SummarizedExperiment)
## For filtering protein coding genes
library(AnnotationHub)
## Differential expression
library(DESeq2)
## Data manipulation and plots
library(tidyverse)
library(ggfortify)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(ComplexHeatmap)
## Transcription network inference
library(RTN)
library(RTNduals)
library(RTNsurvival)
## Clustering and dendrogram ordering
library(ConsensusClusterPlus)
library(dendextend)
library(reshape2)
source("./getConsensusClusterCalls.R")
## MSigDB gene sets
library(msigdbr)

dir.create("./ACC_Cancers2022")
setwd("./ACC_Cancers2022")
```

# TCGA-ACC Cohort
## Code snippet for Supplementary Figure 1 (targets distribution) and Supplementary Table 2 (Regulons' size and number of targets)

```r
# Importing rtni object created in the Preprocessing Step
load("./rtni.RData")

# Checking distribution of targets for the inferred regulons
tni.plot.checks(rtni_tcgaACC, option = "barplot")
# export pdf 6.27 x 5.38"
tni.plot.checks(rtni_tcgaACC, option = "edf")
tni.plot.checks(rtni_tcgaACC, option = "points")

# Retrieving regulons table with size and the number of targets
regulons <- tni.get(rtni_tcgaACC, what = "regulonSize", idkey = "SYMBOL") 
regulons <- regulons[regulons$Size>0,]
regulons$Targets.Dif <- regulons$Positive-regulons$Negative
regulons <- regulons[order(regulons$Targets.Dif),]

# Exporting table - Supplementary Table 2
writexl::write_xlsx(as.data.frame(cbind(
  ENSEMBL=rownames(regulons), regulons)),
  path = "./Regulons.xlsx")

```

## Code snippet for Supplementary table 1 (regulatory network - regulons and target)

```r
# Extracting tnet (table with regulons and targets) from rtni object
tnet <- tni.get(rtni_tcgaACC, what="tnet", idkey = "SYMBOL" )
# Exporting table - Supplementary table 1
write_csv(as.data.frame(cbind(Genes=rownames(tnet), tnet)),
          col_names = T,
          file="./tnet.csv")

```

## Computing Regulon activity - Code snippet for Supplementary table 5
```r
# Computing regulon activity with Two-tailed gene set enrichment analysis
rtni_tcgaACC <- tni.gsea2(rtni_tcgaACC)

# Extracting regulon activity matrix
regact <- tni.get(rtni_tcgaACC, what="regulonActivity", idkey = "SYMBOL")
regact <- as.data.frame(t(regact$differential))

# Exporting regulon activity for TCGA-ACC
writexl::write_xlsx(regact, path = "./Regulon_Activity_TCGA.xlsx")

# Retrieving clinical data for TCGA-ACC samples
coldat <- tni.get(rtni_tcgaACC, what="colAnnotation")

```

## RTNSurvival pipeline - computing survival curves for TCGA-ACC samples

```r
## Importing ACC survival data from Xena Browser
# https://xenabrowser.net/datapages/?dataset=survival%2FACC_survival.txt&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#### Download data
repo_link <- "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
name.file <- "download/survival%2FACC_survival.txt"
download.file(
  url = paste0(repo_link,name.file),
  destfile = "ACC_survival.txt")
ACC_survival<- 
  read.delim("ACC_survival.txt",
             header = T,
             sep="\t")
names(ACC_survival) <- c("sample", "patient","OS","OS_time",
                         "DSS","DSS_time","DFI","DFI_time",
                         "PFI","PFI_time","Redaction")
## Filtering for ACC patients with mRNA data
ACC_survival<- dplyr::filter(ACC_survival, patient %in% coldat$patient)
## Defines ACC-survival data rownames as patient identification
rownames(ACC_survival)<- ACC_survival$patient

## Preparing annotations data.frame
cols<- c("barcode","patient", "age_at_diagnosis", "tumor_stage",
         "days_to_death", "vital_status", 
         "paper_C1A.C1B","paper_mRNA_K4","paper_MethyLevel", 
         "paper_COC", "paper_SCNA.cluster",
         "paper_OncoSign","paper_purity","paper_ploidy",
         "paper_genome_doublings",
         "paper_ADS", "Steroid", "Proliferation_mRNA")
annot<- coldat[,cols]
annot <- annot %>%
  mutate(Death = ifelse(annot$vital_status == "Dead", 1, 0)) %>%
  mutate(age_at_diagnosis = age_at_diagnosis/365)
## Order annot by patient identification
annot <- annot[order(annot$patient),]
identical(annot$patient, ACC_survival$patient) #TRUE
annot <- cbind(annot,
               ACC_survival[,-c(2,11)]) # excluding patient and Redaction
## Changing survival times to months
annot[,c("OS_time", "DSS_time", "DFI_time","PFI_time")] <- 
  lapply(annot[,c("OS_time", "DSS_time", "DFI_time","PFI_time")],
         function(a){a/30})
## Adding dummy variables for Molecular Classifications and tumor stage
annot <- annot %>%
  mutate("Steroid_High" = ifelse(
    annot$Steroid == "Steroid_High", 1, 0)) %>%
  mutate("Steroid_Low" = ifelse(
    annot$Steroid == "Steroid_Low", 1, 0)) %>%
  mutate("Proliferation" = ifelse(
    annot$Proliferation_mRNA == "Proliferation_High", 1, 0)) %>%
  mutate(tumor_stage = if_else(
    annot$tumor_stage == "stage iv", 4, 
    if_else(annot$tumor_stage == "stage iii", 3,
            if_else(annot$tumor_stage == "stage ii", 2,
                    if_else(annot$tumor_stage == "stage i",1, 
                            0))))) %>%
  mutate(C1A = if_else(
    annot$paper_C1A.C1B== "C1A", 1, 0)) %>%
  mutate(C1B = if_else(
    annot$paper_C1A.C1B== "C1B", 1, 0)) %>%
  mutate(CIMP_low = if_else(
    annot$paper_MethyLevel== "CIMP-low", 1, 0 )) %>%
  mutate(CIMP_interm = if_else(
    annot$paper_MethyLevel== "CIMP-intermediate", 1, 0 )) %>%
  mutate(CIMP_high = if_else(
    annot$paper_MethyLevel== "CIMP-high", 1, 0 ))

## Setting patient barcode as rownames
rownames(annot) <- annot[,1]

cols <- c("patient", "age_at_diagnosis", "tumor_stage",
          "OS", "OS_time", "DSS", "DSS_time",
          "DFI", "DFI_time", "PFI", "PFI_time",
          "Steroid_High", "Steroid_Low", "Proliferation",
          "C1A", "C1B", 
          "CIMP_low", "CIMP_interm", "CIMP_high")
# Date.frame for Cox
annot_cox <- annot[,cols] # tumor stage, survival data and binary molecular classification
## No NA values on OS and PFI data
any(is.na(annot_cox[,c("OS","OS_time","PFI","PFI_time")])) #FALSE

names(annot_cox)[2] <- "Age"
names(annot_cox)[3] <- "Tumor.Stage"
annot_cox$Stage1 <- ifelse(annot_cox$Tumor.Stage==1,1,0)
annot_cox$Stage2 <- ifelse(annot_cox$Tumor.Stage==2,1,0)
annot_cox$Stage3 <- ifelse(annot_cox$Tumor.Stage==3,1,0)
annot_cox$Stage4 <- ifelse(annot_cox$Tumor.Stage==4,1,0)

# Selecting regulons with more than 15 positive and negative targets
regs <- tni.get(rtni_tcgaACC, what="regulonSize", idkey = "SYMBOL" )
regs <- regs[regs$Positive>=15 & regs$Negative>=15,] #611 regulons

# RTNS pipeline
rtns <- tni2tnsPreprocess(rtni_tcgaACC, 
                          regulatoryElements = rownames(regs),
                          survivalData = annot_cox,
                          time = 5, #11 for PFI
                          event = 4, #10 for PFI
                          keycovar = c("Age", "Tumor.Stage"),
                          excludeMid = F)
rtns <- tnsGSEA2(rtns)
```

### Cox multivariate analysis for overall survival - Code snippet for Figure 7A (Cox) and Supplementary table 3

```r
### Cox multivariate analysis
rtns <- tnsCox(rtns)
# Retrieving table with Cox results
reg.cox <- rtns@results[["Cox"]][["Table"]]
# Saving the results
writexl::write_xlsx(reg.cox, path="./TCGA_cox.xlsx")

# Reordering Cox results table to merge with Master Table
reg.cox <- reg.cox[regs$SYMBOL,]
identical(rownames(reg.cox), regs$SYMBOL) # T

regs <- cbind(regs,
              Cox.HR = reg.cox$HR,
              Cox.adj.pvalue = reg.cox$Adjusted.Pvalue)

# Creating data.frame for the prognostic regulons
regs.risk <- regs
regs.risk <- regs.risk[regs.risk$Cox.adj.pvalue < 0.05,]
nrow(regs.risk) #369
nrow(regs.risk[regs.risk$Cox.HR < 1,]) #188
nrow(regs.risk[regs.risk$Cox.HR > 1,]) #181

# Creating categorical variable for Cox HR result
regs.risk$HR <- ifelse(regs.risk$Cox.adj.pvalue > 0.05,
                       0,
                       ifelse(regs.risk$Cox.HR > 1,
                              ">1",
                              "<1"))

# Selecting activity only for the 369 prognostic regulons
regact <- regact[regs.risk$SYMBOL,]

# Plotting results
tnsPlotCox(rtns, regs = c("NR5A1", "CENPA"), sortregs = F,
           xlim = c(0.1,15),
           plotpdf = F, #T, 
           fname= "Cox_TCGA", 
           fpath = "./", width = 4, height = 2)
```

### Kaplan-Meier analysis for overall survival (NR5A1 and CENPA) - Code snippet for Figure 7A (KM)

```r
# Kaplan Meier analysis
rtns <- tnsKM(rtns, regs = c("NR5A1", "CENPA"))
                #regs.risk$SYMBOL)
reg.KM <- rtns@results[["KM"]][["Table"]]

# Plotting Kaplan Meier curves
cols <- c("Steroid_High", "Steroid_Low","Proliferation",
          "C1A", "C1B", "Stage1", "Stage2", "Stage3", "Stage4")

atribs <- names(annot_cox[,cols])

tnsPlotKM(rtns, regs=c("NR5A1"), attribs = list(atribs),
          #plotpdf = T, 
          plotpdf = F,
          fname = "KM_TCGA_stage", fpath = "./")

tnsPlotKM(rtns, regs=c("CENPA"), attribs = list(atribs),
          plotpdf = F, #T,
          fname = "KM_TCGA_stage", fpath = "./")

# Cleaning environment
rm(reg.cox, reg.KM, ACC_survival, cols, atribs, res.df)
```

### Cox multivariate analysis for progression-free interval - Code snippet for Supplementary table 3

```r
# rtns for PFI
rtns.PFI <- tni2tnsPreprocess(rtni_tcgaACC, 
                          regulatoryElements = rownames(regs.risk),
                          survivalData = annot_cox,
                          time = 11, #11 for PFI
                          event = 10, #10 for PFI
                          keycovar = c("Age", "Tumor.Stage"),
                          excludeMid = F)
rtns.PFI <- tnsGSEA2(rtns.PFI)

### Cox multivariate analysis
rtns.PFI <- tnsCox(rtns.PFI)
# Generating table to export results
reg.cox.PFI <- rtns.PFI@results[["Cox"]][["Table"]]
writexl::write_xlsx(reg.cox.PFI, path="./TCGA_cox_PFI_v3.xlsx")

# Reorgering results to merge with master table
reg.cox.PFI <- reg.cox.PFI[regs.risk$SYMBOL,]
identical(rownames(reg.cox.PFI), regs.risk$SYMBOL) # T
nrow(reg.cox.PFI[reg.cox.PFI$Adjusted.Pvalue < 0.05,]) #331

regs.risk <- cbind(ENSEMBL = rownames(regs.risk),
                   regs.risk,
                   Cox.HR.PFI = reg.cox.PFI$HR,
                   Cox.adj.pvalue.PFI = reg.cox.PFI$Adjusted.Pvalue)

rm(reg.cox.PFI, reg.KM.PFI, ACC_survival, cols, atribs, res.df)
```
## TNA pipeline - Assessing correlation between Steroid and Proliferation Phenotypes with regulon activity

### Steroid Phenotype - Code snippet for Supplementary table 6
```r
# Transcriptional Network Analysis (TNA)
load("./tcgaACC_Steroid.RData")
######## Differential expression analysis with DESeq2
dds <- DESeqDataSet(tcgaACC, 
                    design = ~ Steroid + Proliferation_mRNA) # for information on design argument, see DESeq2 vignette
dds <- DESeq(dds)

# Presenting D.E. results
res <- results(dds, contrast = c("Steroid",
                                 "Steroid_High",
                                 "Steroid_Low"))

## Generating data.frame with DEGs
identical(rownames(res),rownames(rowData(tcgaACC))) #TRUE
res.df <- cbind(as.data.frame(rowData(tcgaACC)[,c("ensembl_gene_id", "external_gene_name")]),
                as.data.frame(res)) #25476 genes

res.df <- res.df[res.df$baseMean>0,] #24229 genes
res.df <- res.df[order(res.df$padj),]

## Hits for the phenotype - DEGs with adj p < 0.05
hits <- rownames(res.df[res.df$padj<0.05,])
## Log2 Fold Change for definying the phenotype
phenotype <- res.df$log2FoldChange
names(phenotype) <- rownames(res.df)
## Information for the DEGs
phenoID <- cbind(rownames(res.df),
                 res.df[,c("external_gene_name",
                           "ensembl_gene_id")])
# RTNA pipeline - creating rtna object for Steroid Phenotype
rtna <- tni2tna.preprocess(object = rtni_tcgaACC, 
                           phenotype = phenotype,
                           hits = hits,
                           phenoIDs = phenoID)

# Run the MRA method
rtna <- tna.mra(rtna, tfs = regs.risk$SYMBOL)

# Get MRA results;
#..setting 'ntop = -1' will return all results, regardless of a threshold
mra <- tna.get(rtna, what="mra", ntop = -1)

# Exporting MRA results
writexl::write_xlsx(mra, path = "./MRA_Steroid.xlsx")

# Run the GSEA-2T method
rtna <- tna.gsea2(rtna, nPermutations = 1000, tfs = regs.risk$SYMBOL)
gsea2 <- tna.get(rtna, what = "gsea2", ntop = -1)
# Retrieving results
gsea2.dif <- as.data.frame(gsea2$differential)
# Exporting Enrichment analysis result
writexl::write_xlsx(gsea2.dif, path="./GSEA2_Steroid.xlsx")

# Merging GSEA results with master table
identical(regs.risk$SYMBOL, rownames(gsea2.dif)) # F
gsea2.dif <- gsea2.dif[regs.risk$SYMBOL,]
identical(regs.risk$SYMBOL, rownames(gsea2.dif)) # T
regs.risk$GSEA2.Ster <- gsea2.dif$Observed.Score
regs.risk$GSEA2.Ster.Padj <- gsea2.dif$Adjusted.Pvalue

```

### Proliferation phenotype - Code snippet for Supplementary table 6
```r
# Transcriptional Network Analysis (TNA)
######## Results from Differential expression analysis with DESeq2
res.prol <- results(dds, 
                    contrast = c("Proliferation_mRNA",
                                 "Proliferation_High",
                                 "Proliferation_Low"))

## Generating data.frame
identical(rownames(res.prol),rownames(rowData(tcgaACC))) #TRUE
res.prol.df <- cbind(as.data.frame(rowData(tcgaACC)[,2:3]),
                     as.data.frame(res.prol)) #25476 genes

res.prol.df <- res.prol.df[res.prol.df$baseMean>0,] #24229 genes
res.prol.df <- res.prol.df[order(res.prol.df$padj),]

## Hits for the phenotype - DEGs with adj p < 0.05
hits <- rownames(res.prol.df[res.prol.df$padj<0.05,])
## Log2 Fold Change for definying the phenotype
phenotype <- res.prol.df$log2FoldChange
names(phenotype) <- rownames(res.prol.df)
## Information for the DEGs
phenoID <- cbind(rownames(res.df),
                 res.df[,c("external_gene_name",
                           "ensembl_gene_id")])

# RTNA pipeline - creating rtna object for Proliferation phenotype
rtna.prol <- tni2tna.preprocess(object = rtni_tcgaACC, 
                           phenotype = phenotype,
                           hits = hits,
                           phenoIDs = phenoID)

# Run the MRA method
rtna.prol <- tna.mra(rtna.prol, tfs = regs.risk$SYMBOL)

# Get MRA results;
#..setting 'ntop = -1' will return all results, regardless of a threshold
mra <- tna.get(rtna.prol, what="mra", ntop = -1)

# Exporting MRA results
writexl::write_xlsx(mra, path = "./MRA_Proliferation.xlsx")

# Run the GSEA-2T method
rtna.prol <- tna.gsea2(rtna.prol, nPermutations = 1000, 
                       tfs = regs.risk$SYMBOL)
gsea2 <- tna.get(rtna.prol, what = "gsea2", ntop = -1)
# Retrieving results
gsea2.dif.prol <- as.data.frame(gsea2$differential)
# Exporting Enrichment analysis result
writexl::write_xlsx(gsea2.dif.prol, 
                    path="./GSEA2_Proliferation.xlsx")

# Merging GSEA results with master table
identical(regs.risk$SYMBOL, rownames(gsea2.dif.prol)) # F
gsea2.dif.prol <- gsea2.dif.prol[regs.risk$SYMBOL,]
identical(regs.risk$SYMBOL, rownames(gsea2.dif.prol)) # T
regs.risk$GSEA2.Prol <- gsea2.dif.prol$Observed.Score
regs.risk$GSEA2.Prol.Padj <- gsea2.dif.prol$Adjusted.Pvalue

# Cleaning environment
rm(dds, gsea2, mra, phenoID, res, res.df, res.prol, tcgaACC,
   res.prol.df, rtna, rtna.prol, hits, phenotype, rtns.PFI)
```

## Assessing Thorsson master table for ACC samples
```r
###### Thorsson et al. (2018) data-sets for TCGA samples
# Download Immune Features and Cibersort scores from Thorsson et al 2018
repo_link <- "https://www.cell.com/cms/10.1016/j.immuni.2018.03.023/attachment/"
name.file <- "1b63d4bc-af31-4a23-99bb-7ca23c7b4e0a/mmc2.xlsx"
download.file(
  url= paste0(repo_link,name.file),
  destfile = "thorsson2018.xlsx")
thorsson <- 
  data.frame(readxl::read_xlsx(
    path= "./thorsson2018.xlsx"))

## Select patients with Leukocyte Fraction information
thorsson <- thorsson[thorsson$TCGA.Study == "ACC",]
thorsson[, 5:ncol(thorsson)] <- apply(thorsson[,5:ncol(thorsson)], 2, as.numeric)
thorsson <- thorsson[complete.cases(thorsson$Leukocyte.Fraction),]
## Filter for TCGA-ACC cohort and merge with clinical information
thorsson <- subset(thorsson,
                   thorsson$TCGA.Participant.Barcode %in% coldat$patient)
rownames(thorsson) <- thorsson$TCGA.Participant.Barcode
thorsson <- thorsson[coldat$patient,]
rownames(thorsson) <- rownames(coldat)
identical(rownames(annot), rownames(thorsson)) #F
thorsson <- thorsson[rownames(annot),]

identical(rownames(annot), rownames(thorsson)) #T
annot$Leukocyte.Fraction <- thorsson$Leukocyte.Fraction
annot$BCR.Shannon <- thorsson$BCR.Shannon
annot$TCR.Shannon <- thorsson$TCR.Shannon
```

### Correlation between regulon activity and Leukocyte Fraction

```r
identical(rownames(thorsson), colnames(regact)) #F
thorsson <- thorsson[colnames(regact),]
identical(rownames(thorsson), colnames(regact)) #T

LF <- thorsson$Leukocyte.Fraction
names(LF) <- rownames(thorsson)

regact_t <- as.data.frame(t(regact))
# calculating the spearman correlation between regulon activity and leukocyte fraction
LF.cor <- (apply(regact_t, 2, 
                       function(a){cor.test(LF, a, exact = F,
                                            method="spearman")}))
# Creating data frame to merge with master table (Supplementary Table 4)
LF.cor.df <- data.frame(LF.rho = rep(NA, length(LF.cor)),
                        LF.pval = rep(NA, length(LF.cor)), 
                        row.names = names(LF.cor))

for(i in 1:nrow(LF.cor.df)){
  LF.cor.df[i,"LF.rho"] <- 
    LF.cor[[rownames(LF.cor.df)[i]]][["estimate"]][["rho"]]
  LF.cor.df[i,"LF.pval"] <- LF.cor[[rownames(LF.cor.df)[i]]][["p.value"]]
}

LF.cor.df$LF.rho.padj <- p.adjust(LF.cor.df$LF.pval, method = "BH")
identical(rownames(LF.cor.df), regs.risk$SYMBOL) #T

# Mergind with master table
regs.risk <- cbind(regs.risk, 
                   LF.cor.df[,c(1,3)])

rm(LF.cor, LF.cor.df, LF)
```

### Correlation between immune signatures and regulon activity

```r
# Extracting immune signatures values for correlation calculation
dat <- as.data.frame(t(regact))
identical(rownames(thorsson), rownames(dat)) #TRUE
dat <- cbind(dat, 
             LF = thorsson$Leukocyte.Fraction, 
             thorsson[,c(9:14, 26:28)])
cols.immune <- names(dat)[370:379]

regs.risk.cor <- regs.risk[,c("SYMBOL", "HR",
                              "LF.rho", "LF.rho.padj")]
# Creating loop to correlate with all signatures
for(i in c(9:14,26:28)){
  print(names(thorsson)[i])
  dat_imun <- as.numeric(thorsson[,i])
  names(dat_imun) <- rownames(thorsson)
  
  dat_imun.cor <- apply(regact_t, 2, 
                         function(a){cor.test(dat_imun, a, exact=F,
                                              method="spearman" )})
  dat_imun.cor.df <- data.frame(rho = rep(NA, length(dat_imun.cor)),
                        pval = rep(NA, length(dat_imun.cor)), 
                        row.names = names(dat_imun.cor))

  for(k in 1:nrow(dat_imun.cor.df)){
    dat_imun.cor.df[k,"rho"] <- 
      dat_imun.cor[[rownames(dat_imun.cor.df)[k]]][["estimate"]][["rho"]]
    dat_imun.cor.df[k,"pval"] <-
      dat_imun.cor[[rownames(dat_imun.cor.df)[k]]][["p.value"]]
  }

  dat_imun.cor.df$rho.padj <- 
    p.adjust(dat_imun.cor.df$pval, method = "BH")
  
  colnames(dat_imun.cor.df) <- c(paste0(colnames(thorsson)[i], ".rho"),
                                 paste0(colnames(thorsson)[i], ".pvalue"),
                                 paste0(colnames(thorsson)[i], ".padj"))
  
  identical(rownames(dat_imun.cor.df), regs.risk.cor$SYMBOL) #T
  # Mergind with data frame
  regs.risk.cor <- cbind(regs.risk.cor, 
                         dat_imun.cor.df[,c(1,3)])
}
```

## Regulon activity clustering - Code snippet for Supplementary Figure 2

```r

#--- run a simple example, same available from the documentation of the
#--- ConsensusClusterPlus function
dc <- as.matrix(regact_t)
# run consensus cluster, with standard options
met = "maximum"
rcc = ConsensusClusterPlus(dc, maxK=8, reps=1000, pItem=.8, pFeature=1,
                           innerLinkage="ward.D2", finalLinkage="ward.D2",
                           distance=met, clusterAlg="hc", plot="png",
                           title = "maximum_80_100_hc_def_v2", seed = 123
                          )

#--- extract consensus calls from rcc for k=4
res <- getConsensusClusterCalls(rcc, k= 4)

cluster_calls <- res$clusterCalls
hclust_obj <- res$consensusTree

# Assessing dendrogram
dend <- as.dendrogram(hclust_obj)
colors = c("#a6cee3", "#1f78b4", 
           "#fb9a99", "#e31a1c" )
# Definying branches colors for each cluster
dend <- color_branches(dend, 
                       k=4, 
                       col = colors,
                       groupLabels = T)
# View dendrogram
plot(dend)
labs.regs <- labels(dend)

# Assessing 4 branches
dend[[1]][[1]] # 62 lsp
dend[[1]][[2]] # 126 low prol
dend[[2]][[1]] # 113 hsp
dend[[2]][[2]] # 68 high prol

# Retrieving regulons for each cluster
c1 <- labs.regs[1:62]
c2 <- labs.regs[63:188]
c3 <- labs.regs[189:301]
c4 <- labs.regs[302:369]

#--- check clusterCalls
head(cluster_calls[order(cluster_calls$Silhouette),])
#        Item TreeOrder Silhouette ConsensusClass
# NR1I3 NR1I3       191  0.5447255              3
# ZBTB5 ZBTB5       195  0.5516564              3
# NR3C1 NR3C1        16  0.5531460              1
# FOXO6 FOXO6         9  0.5539559              1
# PURB   PURB       194  0.5555581              3
# CREM   CREM       203  0.5639320              3

cluster_calls %>% 
  group_by(ConsensusClass) %>% 
  summarise(mean=mean(Silhouette))
# A tibble: 4 × 2
#   ConsensusClass  mean
#            <int> <dbl>
# 1              1 0.763
# 2              2 0.982
# 3              3 0.751
# 4              4 0.892

# All clusters have silhouette > 0.75
# Definying cluster colors
cluster_calls$col <- ifelse(cluster_calls$ConsensusClass == 1,
                            "#a6cee3",
                            ifelse(cluster_calls$ConsensusClass==2,
                                   "#1f78b4",
                                   ifelse(cluster_calls$ConsensusClass==3,
                                          "#fb9a99", "#e31a1c")))
# Exporting cluster results
 writexl::write_xlsx(cluster_calls,
                     path="./Cluster_calls.xlsx")

#--- plot consensus silhouette
barplot(cluster_calls$Silhouette, col = cluster_calls$col, 
        beside = T, border = cluster_calls$col, las=2)
cluster_calls <- cluster_calls[order(cluster_calls$Item),]
cluster_calls$Item <- factor(cluster_calls$Item, 
                             levels = cluster_calls$Item)
cluster_calls$ConsensusClass <- as.factor(cluster_calls$ConsensusClass)

# Cluster matrix with silhouette
cluster_matrix <- as.data.frame(rcc[[4]][["consensusMatrix"]])

cluster_calls$ConsensusClass <- as.factor(cluster_calls$ConsensusClass)

cluster_calls <- cluster_calls[order(cluster_calls$Item),]

# Plotting results for consensus clustering
h1 <-Heatmap(
  cluster_matrix,
  name = "Silhouette",
  column_title = "Consensus Matrix k=4",
  column_title_gp = gpar(fontsize=8, fontface="bold"),
  show_row_dend = F,
  show_column_names = F,
  cluster_rows=as.dendrogram(hclust_obj),
  cluster_columns=as.dendrogram(hclust_obj),
  heatmap_legend_param = list(
    direction = "horizontal",
    at=c(0, .5, 1),
    title="Consensus\nmatrix",
    title_gp = gpar(fontsize=8,
                    fontface="bold"),
    labels_gp=gpar(fontsize=8),
    border="black",
    title_position="topcenter"),
  col = circlize::colorRamp2(breaks = c(0,.5,1),
                             colors=c("white", "blue","blue3")),
  top_annotation = columnAnnotation(
    df = cluster_calls[,
                       "ConsensusClass", drop=F],
    simple_anno_size=unit(3,"mm"),
    height=unit(1,"cm"),
    annotation_name_gp = gpar(fontsize=0),
    annotation_legend_param= list(
      ConsensusClass=list(
        title="Regulon\nCluster", 
        title_position="topcenter",
        title_gp = gpar(fontsize=8, fontface="bold"),
        labels = c("RC1", "RC2", "RC3", "RC4"),
        labels_gp=gpar(fontsize=8),
        nrow=2,
        border="black",
        direction="vertical")),
    col= list(
      ConsensusClass=c("1"="#a6cee3",
                       "2"="#1f78b4",
                       "3"="#fb9a99",
                       "4"="#e31a1c"))),
  bottom_annotation = columnAnnotation(
    Silhouette = anno_barplot(
      cluster_calls[,
                    "Silhouette", drop=F],
      gp = gpar(fill = cluster_calls[order(cluster_calls$Item),"col"],
                col = cluster_calls[order(cluster_calls$Item),"col"]),
      annotation_label_gp = gpar(fontsize=6))))
draw(h1, heatmap_legend_side="bottom", 
     annotation_legend_side="bottom", merge_legends=T)

# export pdf 4.98 x 5.36 in

rownames(cluster_calls) <- cluster_calls$Item
cluster_calls <- cluster_calls[regs.risk$SYMBOL,]

identical(rownames(cluster_calls), regs.risk$SYMBOL) #T

# Merging clustering results with master table
regs.risk$Cluster <- 
  ifelse(cluster_calls$ConsensusClass==1, "RC1",
         ifelse(cluster_calls$ConsensusClass==2, "RC2",
                ifelse(cluster_calls$ConsensusClass==3, "RC3","RC4")))
regs.risk$Silhouette <- cluster_calls$Silhouette
regs.risk$Tree.order <- cluster_calls$TreeOrder

```

## Code snippet for Figure 2A

```r
# Main heatmap
# Top (columns) annotation
cols <- c("tumor_stage", "Steroid", "Proliferation",
          "paper_C1A.C1B", "paper_purity",
          "Leukocyte.Fraction")
colannot <- annot[,cols]

colnames(colannot) <- c("Tumor.Stage", "Steroid", "Proliferation",
                        "C1A.C1B", "Purity",
                        "Leukocyte.Fraction")

# Left (rows) annotation
rowannot <- regs.risk
rownames(rowannot) <- rowannot$SYMBOL
rowannot$Leukocyte.Fraction.rho <- ifelse(rowannot$LF.rho.padj < 0.05,
                                          rowannot$LF.rho,
                                          0)

identical(rownames(rowannot), rownames(regact)) #T

identical(rownames(rowannot), rownames(gsea2.dif)) # T
rowannot$GSEA2.Ster <- ifelse(gsea2.dif$Adjusted.Pvalue > 0.01,
                              0,
                              gsea2.dif$Observed.Score)

identical(rownames(rowannot), rownames(gsea2.dif.prol)) # T
rowannot$GSEA2.Prol <- ifelse(gsea2.dif.prol$Adjusted.Pvalue > 0.01,
                              0,
                              gsea2.dif.prol$Observed.Score)

identical(colnames(regact), rownames(colannot)) # F
colannot <- colannot[colnames(regact),]
identical(colnames(regact), rownames(colannot)) # T
colannot$Tumor.Stage[colannot$Tumor.Stage == 0] <- "NA"
colannot$Tumor.Stage <- as.factor(colannot$Tumor.Stage)

colannot$Purity <- as.numeric(gsub(",", ".", colannot$Purity))

# Plotting main heatmap
h2<- Heatmap(
  matrix = regact,
  name = "Regulon_Activity",
  show_column_names = F,
  cluster_rows = dend, 
  row_dend_side = "left",
  show_row_names = F, 
  row_split = 4,
  row_gap = unit(0.3,"mm"),
  use_raster = F,
  show_column_dend = T,
  clustering_distance_columns = "maximum",
  clustering_method_columns = "ward.D2",
  column_title = "TCGA-ACC cohort (n=78)",
  column_title_rot = 0,
  column_title_gp = gpar(fontsize=8, fontface="bold"),
  column_title_side = "top",
  show_row_dend = T,  
  row_title_gp = gpar(fontsize=8, fontface="bold"),
  row_title_rot = 90,
  row_title = c("369 prognostic regulons"),
  column_names_gp = gpar(fontsize=8), 
  column_names_side = "top", 
  column_names_centered = T,
  column_names_rot = 90,
  col = circlize::colorRamp2(
    breaks = seq(-2,2, length.out=7),
    colors=rev(RColorBrewer::brewer.pal(7,"RdBu"))),
  heatmap_legend_param = list(
    direction = "horizontal",
    at=c(-2, 0, 2),
    title="Regulon\nActivity",
    title_gp = gpar(fontsize=6, fontface="bold"),
    labels_gp=gpar(fontsize=6),
    title_position="topcenter"),
  top_annotation = columnAnnotation(
    df=colannot,
    simple_anno_size=unit(3,"mm"),
    height=unit(1,"cm"),
    annotation_name_side="left", 
    annotation_name_gp=gpar(fontsize=8),
    annotation_label=c("Tumor Stage", "Steroid Phenotype",
                       "Proliferation", "C1A/C1B", "Purity",
                       "Leukocyte Fraction"),
    show_legend=c(T, T, T, T, T, F),
    na_col = "white",
    annotation_legend_param= list(
      Steroid=list(
        title="Steroid", 
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        border="black",
        direction="vertical"),
      Proliferation=list(
        title="Proliferation", 
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        at = c(0,1),
        labels = c("Low Proliferation", "High Proliferation"),
        labels_gp=gpar(fontsize=6),
        border="black",
        direction="vertical"),
      C1A.C1B=list(
        title="C1A.C1B", 
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        border="black",
        direction="vertical"),
      Leukocyte.Fraction=list(
        title="Leukocyte\nFraction", 
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        border="black",
        direction="horizontal",
        legend_width= unit(1.5,"cm"),
        at=c(0,.2,.4)),
      Purity=list(
        title="Purity", 
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        border="black",
        direction="horizontal",
        legend_width= unit(1.5,"cm"),
        at=c(0,.5,1)),
      Tumor.Stage=list(
        title="Tumor Stage", 
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        nrow=2,
        border="black",
        direction="vertical")),
    col= list(
      Tumor.Stage= c("NA"="white",
                     "1"="gray90",
                     "2"="gray66",
                     "3"="gray33",
                     "4"="black"),
      Leukocyte.Fraction= 
        circlize::colorRamp2(c(0,.4),
                             c("gray90","darkgreen")),
      Purity= 
        circlize::colorRamp2(c(0, .5, 1),
                             c("gray90","gray50", "navy")),
      Steroid= c("Steroid_High" ="black",
                 "Steroid_Low"= "gray90"),
      C1A.C1B= c("C1A" ="black",
                 "C1B"= "gray90"),
      Proliferation=c("1" = "black",
                      "0" = "gray90"))),
  left_annotation = rowAnnotation(
    df=rowannot[c("HR", "GSEA2.Ster",
                  "GSEA2.Prol", "LF.rho")],
    annotation_label = c("HR", "Steroid IP",
                         "Prolif. IS", "LF rho"),
    annotation_name_gp=gpar(fontsize=6),
    annotation_name_rot=45,
    gap = unit(1,"mm"),
    simple_anno_size=unit(2,"mm"),
    show_legend=T,
    annotation_legend_param= list(
      HR=list(
        title="Hazard\nRatio",
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        direction="horizontal"),
      LF.rho=list(
        title="Leuk. Fraction\nCorrelation",
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        direction="horizontal",
        at=c(-1, 0, 1),
        legend_width= unit(1.5,"cm")),
      GSEA2.Ster=list(
        title="Steroid\nIP Score", 
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        direction="horizontal",
        at=c(-2, 0, 2),
        legend_width= unit(1.5,"cm")),
      GSEA2.Prol=list(
        title="Proliferation\nIS Score",
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        direction="horizontal",
        at=c(-2, 0, 2),
        legend_width= unit(1.5,"cm"))),
    col= list(
      HR = c("<1" = "#1B7837",
             ">1" = "#BF812D"),
      LF.rho = 
        circlize::colorRamp2(c(-.75,0,.75),
                             c("darkorange", "white", "blue")),
      GSEA2.Ster = 
        circlize::colorRamp2(c(-2,0,2),
                             c("#00BFC4","white","#F8766D")),
      GSEA2.Prol = 
        circlize::colorRamp2(c(-2,0,2),
                             c( "chartreuse4", "white", "red")))))
# Plotting right annotation (with Regulon cluster)
clusters <- rowAnnotation(
  df=rowannot[,c("Cluster"), drop=F], 
  annotation_label = "Cluster",
  annotation_name_gp=gpar(fontsize=6),
  annotation_name_rot=45,
  simple_anno_size=unit(3,"mm"),
  show_legend=T,
  annotation_legend_param= list(
    Cluster=list(
      title="Regulon\nCluster",
      title_position="topcenter",
      title_gp = gpar(fontsize=6, fontface="bold"),
      labels_gp=gpar(fontsize=6),
      direction="horizontal",
      nrow=2)),
  col= list(
    Cluster = c("RC1" = "#a6cee3",
                "RC2" = "#1f78b4",
                "RC3" = "#fb9a99",
                "RC4" = "#e31a1c")))
# Merging main heatmap with right annotation
ht_list = draw(h2+clusters, 
               heatmap_legend_side="bottom",
               annotation_legend_side="bottom")
# Exporting results
g0 <- grid.grabExpr(draw(ht_list,
                         heatmap_legend_side="bottom", 
                         annotation_legend_side = "bottom",
                         merge_legends=T),
                   height = 7.5, width = 6)
plot_grid(g0)
ggsave(g0, filename = "./2A Heatmap_regact_TCGA.pdf", 
        device="pdf", width = 6, height = 7.5, units="in")

```
## Theme for plots

```r
theme_plots <- theme(
        plot.title=element_text(hjust=0),
        legend.position="none",
        legend.margin = margin(-5, -5, 5, -5),
        legend.spacing.y = unit(0, 'cm'), 
        title = element_text(size=10), 
        axis.title.y = element_text(size=8, face = "bold"),
        axis.title.x = element_text(size=8, face = "bold"),
        axis.text.y=element_text(size=8.5), 
        axis.text.x = element_text(size=8, face = "bold"),
        axis.ticks = element_line(linetype=1, color="grey"),
        legend.text=element_text(size=8),
        legend.title= element_text(face="bold", size=8),
        legend.key=element_rect(size=10, color=NA, fill=NA), 
        legend.key.size=unit(5,"mm"),
        legend.direction = "vertical",
        legend.box = "vertical",
        axis.line = element_blank(),
        panel.grid.major.x = element_line(linetype=111, 
                                          color="grey80", 
                                          size = 0.4),
        panel.grid.major = element_line(linetype=3, 
                                        color="grey", 
                                        size = 0.2),
        panel.background = element_rect(fill = "grey98",
                                        colour = "grey50"),
        panel.border = element_rect(colour = "grey", 
                                    fill=NA, size=1))

blankPlot <- ggplot()+
  geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(), 
   panel.border = element_blank(),
   panel.background = element_blank(),
   axis.title.x = element_blank(),
   axis.title.y = element_blank(),
   axis.text.x = element_blank(), 
   axis.text.y = element_blank(),
   axis.ticks = element_blank())
```

## Code snippet for Figure 2B, C, D

```r
rm(cluster_calls, cluster_matrix, dc, g0,  h1, h2,  
   ht_list,  rcc, res,  met)

regs.risk$Cluster <- as.factor(regs.risk$Cluster)
summary(regs.risk$Cluster)
#  C1  C2  C3  C4 
# 62  126 113  68 

# Selecting features to show in boxplot
dat <- regs.risk[,c("HR", "LF.rho", 
                    "GSEA2.Ster", "GSEA2.Prol")]
cols <- c("Hazard Ratio", "Leukocyte Fraction\nSpearman Correlation",
          "Steroid IP\nScore", "Proliferation IS\nScore")
# Create empty list to save plots
g <- list()
# Loop to generate a plot for each feature analyzed
for(i in 2:ncol(dat)){
  lab <- cols[i]
  ymax <- round(max(dat[,i]))
  g1<- eval(substitute(
    ggplot(dat,
           aes(x=HR, y= dat[,i], fill=HR))+
    geom_boxplot(width=c(0.5), lwd=0.2, outlier.size = 0.5)+  
    geom_violin(
      aes(fill=NULL, col=HR),
      alpha=0.4, 
      show.legend = F)+
    geom_jitter(size=0.2, width = 0.1)+
    geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
    scale_x_discrete(labels=element_blank())+
    scale_y_continuous(expand = expansion(mult = c(0.2,0.2)),
                       breaks=c(-ymax,-ymax/2,0,ymax/2,ymax))+
    scale_color_manual(values=c("<1" = "#1B7837",
                                ">1" = "#BF812D"))+
    scale_fill_manual(values=c("<1" = "#1B7837",
                               ">1" = "#BF812D"),
                      labels = c("Regulons with HR <1 (n=188)",
                                 "Regulons with HR >1 (n=181)"))+ 
    labs(y= lab, 
         x= element_blank())+
    theme_plots+
    stat_compare_means(comparisons = list(c("<1",">1")),
                       label = "p.signif",
                       method = "wilcox",
                       label.x = 1.5)
    ,list(i=i)))
  g[[i-1]] <- g1
}

# Definying x axis labels
g1 <- g1+scale_x_discrete(labels = c("Low\nRisk\nRegulons",
                                     "High\nRisk\nRegulons"))
g[[3]] <- g1
# Assessing legend
leg <- get_legend(
  g1+
    theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow=2, 
                             title = "Prognostic regulons",
                             title.hjust = .5)))
# generating panel with 3 boxplots
g <- plot_grid(g[[1]], g[[2]], g[[3]], leg, ncol = 1, 
               rel_heights = c(1,1,1,.25),
               labels = c("B","C","D",""))
g
# Saving panel
ggsave(g, filename = "./1_boxplot_risk_regs.pdf",
       device="pdf", height = 7.5, width = 2)
```

## Code snippet for Figure 3A, B, C, D, E, and F

```r
# Selecting variables
dat <- regs.risk[,c("Cluster","Cox.HR", "LF.rho", 
                    "GSEA2.Ster", "GSEA2.Prol")]
cols <- c("Cluster", "Hazard Ratio", 
          "Leukocyte Fraction Correlation",
          "Steroid IP\nscore", "Proliferation IS\nscore")
# Loop through the variables to generate and store the plots in a list
g <- list()
for(i in 3:ncol(dat)){
  form <- as.formula(paste0(names(dat[i]), "~", "Cluster"))
  stat.test <- rstatix::dunn_test(data=dat, 
                                  formula=form, 
                                  p.adjust.method = "BH")
  stat.test$p.adj <- format(stat.test$p.adj, 
                            scientific = T, 
                            digits = 3)
  stat.test$size <- ifelse(stat.test$p.adj.signif== "ns", 2, 4)
  size <- c()
  for(j in 1:nrow(stat.test)){
    size <- c(size, rep(stat.test$size[j], 3))
  }
  lab <- cols[i]
  ymax <- round(max(dat[,i]))
  ypos <- c(ymax, ymax*1.2, ymax*1.6, ymax, ymax*1.4, ymax)
  g1<- eval(substitute(
    ggplot(dat,
          aes(x=Cluster, y= dat[,i], fill=Cluster))+
    geom_boxplot(width=c(0.5), lwd=0.2, outlier.size = 0.5)+  
    geom_violin(
      aes(fill=NULL, col=Cluster),
      alpha=0.4, 
      show.legend = F)+
    geom_jitter(size=0.2, width = 0.1)+
    geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
    scale_y_continuous(expand = expansion(mult = c(0.1,0.1)),
                       breaks=c(-ymax,-ymax/2,0,ymax/2,ymax))+
    scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                                "RC3"="#fb9a99","RC4"="#e31a1c"))+
    scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                               "RC3"="#fb9a99","RC4"="#e31a1c"),
                      labels = c("RC1 (n=62)", "RC2 (n=126)",
                                 "RC3 (n=113)", "RC4 (n=68)"))+ 
    labs(y= lab, 
         x= element_blank())+
    theme_plots+
    stat_compare_means(method = "kruskal",
                       label.y = ymax*1.8,
                       size = 2.5)+
    stat_pvalue_manual(stat.test, 
                       inherit.aes=F,
                       label="p.adj.signif", 
                       bracket.shorten = .1,
                       tip.length = .02,
                       y.position = ypos, 
                       label.size = size)
    ,list(i=i)))
  g[[i-2]] <- g1
}

# Boxplot for Hazard Ratio (log10 transform in y axis)
stat.test <- rstatix::dunn_test(data=dat, 
                                formula=Cox.HR~Cluster, 
                                p.adjust.method = "BH")
stat.test$p.adj <- format(stat.test$p.adj, 
                          scientific = T, 
                          digits = 3)
stat.test$size <- ifelse(stat.test$p.adj.signif== "ns", 2, 4)
size <- c()
for(j in 1:nrow(stat.test)){
  size <- c(size, rep(stat.test$size[j], 3))
}
ymax <- round(max(dat$Cox.HR))
ypos <- c(ymax*1.3, ymax*1.7, ymax*3, ymax*1.3, ymax*2.3, ymax*1.3)
g4<- 
  ggplot(dat,
         aes(x=Cluster, y= Cox.HR, fill=Cluster))+
  geom_boxplot(width=c(0.5), lwd=0.2, outlier.size = 0.5)+  
  geom_violin(
    aes(fill=NULL, col=Cluster),
    alpha=0.4, 
    show.legend = F)+
  geom_jitter(size=0.2, width = 0.1)+
  geom_hline(yintercept=1, linetype=2, color="#800026", size=.25)+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)),
                     breaks=c(0, 0.25, 0.5, 1,2,3,5),
                     trans="log10")+
  scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c"))+
  scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                             "RC3"="#fb9a99","RC4"="#e31a1c"),
                    labels = c("RC1 (n=62)", "RC2 (n=126)",
                               "RC3 (n=113)", "RC4 (n=68)"))+ 
  labs(y= "Hazard Ratio", 
       x= element_blank())+
  theme_plots+
  stat_compare_means(method = "kruskal",
                     label.y = log10(ymax*4.2),
                     size = 2.5)+
  stat_pvalue_manual(stat.test, 
                     inherit.aes=F,
                     label="p.adj.signif", 
                     bracket.shorten = .1,
                     tip.length = .02,
                     y.position = log10(ypos), 
                     label.size = size)
# Assessing legend
leg <- get_legend(
  g4+
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow=1,
                             title= "Regulon cluster", 
                             title.hjust = .5)))
# Generating scatter plots (Fig 3E and F)
## Fig 3E
scatterPlot <- 
  ggplot(regs.risk, aes(x=LF.rho, y=GSEA2.Ster ))+
  geom_point(size=.5, aes(color=Cluster))+
  labs(y = "Steroid IP\nscore", 
       x = "Leukocyte Fraction correlation",
       color="Regulon Cluster")+
  stat_cor(method="spearman", 
           size = 2.5, 
           label.x = .2,
           label.y = 1.5)+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)),
                     breaks=c(-2,-1,0,1,2))+
  geom_rug(aes(color=Cluster), show.legend = F)+
  scale_color_manual(
      values=c("RC1"="#a6cee3","RC2"="#1f78b4",
               "RC3"="#fb9a99","RC4"="#e31a1c"),
      labels = c("RC1 (n=62)", "RC2 (n=126)",
                 "RC3 (n=113)", "RC4 (n=68)"),
      na.value="white")+
  theme_plots
## Getting legend
leg2 <- get_legend(
  scatterPlot + 
  theme(legend.position = "bottom")+
  guides(color = guide_legend(
    nrow = 1,
    title.position = "top", 
    title.hjust = .5,
    override.aes=list(size=5))))
## Fig 3F
scatterPlot2 <- 
  ggplot(regs.risk, aes(x=LF.rho, y=GSEA2.Prol))+
  geom_point(size=.5, aes(color=Cluster))+
  labs(y = "Proliferation IS\nscore", 
       x = "Leukocyte Fraction correlation",
       color="Regulon Cluster")+
  stat_cor(method="spearman", 
           size = 2.5, 
           label.x = .2,
           label.y = 1.5)+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.15)),
                     breaks=c(-2,-1,0,1,2))+
  geom_rug(aes(color=Cluster), show.legend = F)+
  scale_color_manual(
      values=c("RC1"="#a6cee3","RC2"="#1f78b4",
               "RC3"="#fb9a99","RC4"="#e31a1c"),
      na.value="white")+
  theme_plots+
  guides(color=guide_legend(nrow= 1, 
                            title.position = "top"))
# Creating panel with the 6 plots
g0 <- plot_grid(g4, g[[1]], g[[2]], g[[3]],
              scatterPlot, scatterPlot2,
              leg, leg2,
              ncol=2, align = "hv",
              rel_heights = c(1,1,1,.25),
              labels = c("A","B","C","D","E","F", "", ""))
g0
# Saving plots
ggsave(g0, filename = "./3_clusters_boxplot_scatter.pdf", device="pdf", 
       width = 7, height = 9)
```

## Code snippet for Supplementary Figure 4

```r
# Selecting variables (immune signature rho)
identical(rownames(regs.risk.cor), rownames(regs.risk)) #T
regs.risk.cor$Cluster <- regs.risk$Cluster
dat <- regs.risk.cor[,c(seq(from=5, to=23, by=2))]
names(dat)
# Looping through variables and storing plots in a list
g <- list()
for(i in 1:9){
  print(names(dat)[i])
  form <- as.formula(paste0(names(dat)[i], "~", "Cluster"))
  stat.test <- rstatix::dunn_test(data=dat, 
                                  formula=form, 
                                  p.adjust.method = "BH")
  stat.test$p.adj <- format(stat.test$p.adj, 
                            scientific = T, 
                            digits = 3)
  print(stat.test)
  stat.test$size <- ifelse(stat.test$p.adj.signif== "ns", 2, 4)
  size <- c()
  for(j in 1:nrow(stat.test)){
    size <- c(size, rep(stat.test$size[j], 3))
  }
  lab <- names(dat)[i]
  if(i==4){
    lab <-"Lymphocyte.Inf.Sig.Score.rho"
  }
  ymax <- round(max(dat[,i]))
  if(ymax == 0){ymax = 1}
  ypos <- c(ymax, ymax*1.2, ymax*1.6, ymax, ymax*1.4, ymax)
  g1<- eval(substitute(
    ggplot(dat,
          aes(x=Cluster, y= dat[,i], fill=Cluster))+
    geom_boxplot(width=c(0.5), lwd=0.2, outlier.size = 0.5)+  
    geom_violin(
      aes(fill=NULL, col=Cluster),
      alpha=0.4, 
      show.legend = F)+
    geom_jitter(size=0.2, width = 0.1)+
    geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
    scale_y_continuous(expand = expansion(mult = c(0.05,0.1)),
                       breaks=c(-1,0,1), limits = c(-1,1.8))+
    scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                                "RC3"="#fb9a99","RC4"="#e31a1c"))+
    scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                               "RC3"="#fb9a99","RC4"="#e31a1c"),
                      labels = c("RC1 (n=62)", "RC2 (n=126)",
                                 "RC3 (n=113)", "RC4 (n=68)"))+ 
    labs(y= lab, 
    #"Spearman correlation with regulon activity", 
         x= element_blank())+
    theme_plots+
    stat_compare_means(method = "kruskal",
                       label.y = ymax*1.8,
                       size = 2.5)+
    stat_pvalue_manual(stat.test, 
                       inherit.aes=F,
                       label="p.adj.signif", 
                       bracket.shorten = .1,
                       tip.length = .02,
                       y.position = ypos, 
                       label.size = size)
    ,list(i=i)))
  g[[i]] <- g1
}
# Getting legend
leg <- get_legend(
  g1+
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow=1,
                             title= "Regulon cluster", 
                             title.hjust = .5)))
# Creating panel
## For 6 immune signatures
row1 <- plot_grid(g[[1]], g[[2]], g[[3]],
                  g[[4]], g[[5]], g[[6]],
                  ncol=3)
## For 3 TCR metrics
row2 <- plot_grid(g[[7]], g[[8]], g[[9]],
                  ncol=3)
## adding legend
leg <- plot_grid(blankPlot, leg, blankPlot,
                 ncol = 3)
# Plotting panel
g0 <- plot_grid(row1, row2, leg, 
                ncol=1,labels = c("A","B",""),
                rel_heights = c(2,1,.25))
g0
# Saving
ggsave(g0, 
       filename = "./Sup_4_immune_signatures_rho.pdf", 
       device="pdf",
       width = 7, height=8)

```
## Code snippet for Figure 5A
```r
# TCGA C1A C1B
# Calculating median activity within regulon clusters
dat <- as.data.frame(regact)

dat.rc1 <- dat[c1[c1 %in% rownames(dat)],]
dat.rc1 <- apply(dat.rc1, 2, median)

dat.rc2 <- dat[c2[c2 %in% rownames(dat)],]
dat.rc2 <- apply(dat.rc2, 2, median)

dat.rc3 <- dat[c3[c3 %in% rownames(dat)],]
dat.rc3 <- apply(dat.rc3, 2, median)

dat.rc4 <- dat[c4[c4 %in% rownames(dat)],]
dat.rc4 <- apply(dat.rc4, 2, median)

identical(rownames(dat.rc1), rownames(dat.rc2)) #T
identical(rownames(dat.rc1), rownames(dat.rc3)) #T
identical(rownames(dat.rc3), rownames(dat.rc4)) #T

dat <- as.data.frame(
  cbind(RC1 = dat.rc1,
        RC2 = dat.rc2,
        RC3 = dat.rc3,
        RC4 = dat.rc4))

identical(rownames(dat), rownames(coldat))#T

dat <- cbind(dat,
             C1A.C1B = coldat$paper_C1A.C1B,
             Vital.Status = coldat$vital_status)

dat <- dat[complete.cases(dat$C1A.C1B),]
dat$C1A.C1B <- as.factor(dat$C1A.C1B)

dat_bp <- pivot_longer(dat,
                       cols = -contains(c("C1A.C1B",
                                          "Vital.Status")),
                       values_to = "value",
                       names_to = "Cluster")

# C1A
stat.test <- rstatix::dunn_test(
    data=dat_bp[dat_bp$C1A.C1B=="C1A",], 
    formula= value ~ Cluster, 
    p.adjust.method = "BH")
stat.test$p.adj <- format(stat.test$p.adj, 
                          scientific = T, 
                          digits = 3)
print(stat.test)
stat.test$size <- ifelse(stat.test$p.adj.signif== "ns", 2, 4)
size <- c()
for(j in 1:nrow(stat.test)){
  size <- c(size, rep(stat.test$size[j], 3))
}
ymax <- round(max(dat_bp$value))
ymax <- 1.8
ypos <- c(ymax, ymax*1.2, ymax*1.6, ymax, ymax*1.4, ymax)

g1.C1A.tcga<-  
  ggplot(dat_bp[dat_bp$C1A.C1B=="C1A",],
         aes(x=Cluster, y=value, fill=Cluster))+
  labs(title= "C1A\nTCGA (n=43)",
       y= "Median Regulon\nCluster Activity", 
       x= element_blank()) +
  geom_violin(
    aes(fill=NULL, col=Cluster),
    alpha=0.4, 
    show.legend = c(T))+
  ggbeeswarm::geom_beeswarm(
    aes(color=Vital.Status),
    alpha =.5, 
    cex = 2, size=1, show.legend = T)+
  geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)),
                     breaks=c(-2,-1,0,1,2))+
  scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c",
                              "Alive"="gray50", "Dead"="red"))+
  scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c"),
                    labels = c("RC1 (n=62)", "RC2 (n=126)",
                               "RC3 (n=113)", "RC4 (n=68)"))+ 
  theme_plots +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10, 
                                  face = "bold"),
          axis.title.y = element_text(size=8))+
  stat_compare_means(method = "kruskal",
                     label.y = ymax*1.8,
                     size = 2.5)+
  stat_pvalue_manual(stat.test, 
                     inherit.aes=F,
                     label="p.adj.signif", 
                     bracket.shorten = .1,
                     tip.length = .02,
                     y.position = ypos, 
                     label.size = size)

## Tumor Stages 3 and 4
stat.test <- rstatix::dunn_test(
    data=dat_bp[dat_bp$C1A.C1B=="C1B",], 
    formula= value ~ Cluster, 
    p.adjust.method = "BH")
stat.test$p.adj <- format(stat.test$p.adj, 
                          scientific = T, 
                          digits = 3)
print(stat.test)
stat.test$size <- ifelse(stat.test$p.adj.signif== "ns", 2, 4)
size <- c()
for(j in 1:nrow(stat.test)){
  size <- c(size, rep(stat.test$size[j], 3))
}

g2.C1B.TCGA <-  
  ggplot(dat_bp[dat_bp$C1A.C1B=="C1B",],
         aes(x=Cluster, y=value, fill=Cluster))+
  labs(title= "C1B\nTCGA (n=35)",
       y= "Median Regulon\nCluster Activity", 
       x= element_blank())+
  geom_violin(
    aes(fill=NULL, col=Cluster),
    alpha=0.4, 
    show.legend = c(T))+
  ggbeeswarm::geom_beeswarm(
    aes(color=Vital.Status), 
    alpha =.5, cex = 2, size=1,
    show.legend = T)+
  geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.15)),
                     breaks=c(-2,-1,0,1,2))+
  scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c",
                              "Alive"="gray50", "Dead"="red"))+
  scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                             "RC3"="#fb9a99","RC4"="#e31a1c"),
                    labels = c("RC1 (n=62)", "RC2 (n=126)",
                               "RC3 (n=113)", "RC4 (n=68)"))+ 
  theme_plots +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10, 
                                  face = "bold"),
        axis.title.y = element_text(size=8))+
  stat_compare_means(method = "kruskal",
                     label.y = ymax*1.8,
                     size = 2.5)+
  stat_pvalue_manual(stat.test, 
                     inherit.aes=F,
                     label="p.adj.signif", 
                     bracket.shorten = .1,
                     tip.length = .02,
                     y.position = ypos, 
                     label.size = size)
# Rearranging labels
g1.C1A.tcga <- g1.C1A.tcga +
  labs(title= "C1A TCGA-ACC\n(n=43 samples)",
       y= "Regulon activity (median of the cluster)", 
       x= element_blank())+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 8, 
                                  face = "bold"),
        axis.title.y = element_text(size=6))  

g2.C1B.TCGA <-  g2.C1B.TCGA+
  labs(title= "C1B TCGA-ACC\n(n=35 samples)",
       y= element_blank(), 
       x= element_blank())+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 8, 
                                  face = "bold"))    

```

## Code snippet for Figure 5C

```r
identical(rownames(colannot), rownames(regact_t)) # T
identical(rownames(colannot), colnames(regact)) #T
dat <- regact

# Calculating median of regulon activity for each cluster
dat.rc1 <- dat[c1,]
dat.rc1 <- apply(dat.rc1, 2, median)

dat.rc2 <- dat[c2,]
dat.rc2 <- apply(dat.rc2, 2, median)

dat.rc3 <- dat[c3,]
dat.rc3 <- apply(dat.rc3, 2, median)

dat.rc4 <- dat[c4,]
dat.rc4 <- apply(dat.rc4, 2, median)

identical(rownames(dat.rc1), rownames(dat.rc2)) #T
identical(rownames(dat.rc1), rownames(dat.rc3)) #T
identical(rownames(dat.rc3), rownames(dat.rc4)) #T

dat <- as.data.frame(
  cbind(RC1 = dat.rc1,
        RC2 = dat.rc2,
        RC3 = dat.rc3,
        RC4 = dat.rc4))

identical(rownames(dat), rownames(annot))#F
dat <- dat[rownames(annot),]
identical(rownames(dat), rownames(annot))#T

dat <- cbind(dat,
             Tumor.Stage = annot$tumor_stage,
             Vital.Status = annot$vital_status)
# Filtering samples without tumor stage annotation
dat <- dat[dat$Tumor.Stage != 0,]
# Merging stages low and high stages
dat$Tumor.Stage <- as.factor(ifelse(dat$Tumor.Stage %in% c(1,2),
                                    "Stage 1 and 2",
                                    "Stage 3 and 4"))
# Adapting data frame to boxplot
dat_bp <- pivot_longer(dat,
                       cols = -contains(c("Tumor.Stage",
                                          "Vital.Status")),
                       values_to = "value",
                       names_to = "Cluster")

# Plots for Tumor stages 1 and 2
## Statistical test
stat.test <- rstatix::dunn_test(
    data=dat_bp[dat_bp$Tumor.Stage=="Stage 1 and 2",], 
    formula= value ~ Cluster, 
    p.adjust.method = "BH")
stat.test$p.adj <- format(stat.test$p.adj, 
                          scientific = T, 
                          digits = 3)
print(stat.test)
stat.test$size <- ifelse(stat.test$p.adj.signif== "ns", 2, 4)
size <- c()
for(j in 1:nrow(stat.test)){
  size <- c(size, rep(stat.test$size[j], 3))
}
ymax <- round(max(dat_bp$value))
ypos <- c(ymax, ymax*1.2, ymax*1.6, ymax, ymax*1.4, ymax)
## Plot
g1<-  
  ggplot(dat_bp[dat_bp$Tumor.Stage=="Stage 1 and 2",],
         aes(x=Cluster, y=value, fill=Cluster))+
  labs(title= "Tumor Stage 1 and 2\nTCGA-ACC (n=46)",
       y= "Median Regulon\nCluster Activity", 
       x= element_blank()) +
  geom_violin(
    aes(fill=NULL, col=Cluster),
    alpha=0.4, 
    show.legend = c(T))+
  ggbeeswarm::geom_beeswarm(
    aes(color=Vital.Status),
    alpha =.5, 
    cex = 2, size=1, show.legend = T)+
  geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.1)),
                     breaks=c(-2,-1,0,1,2))+
  scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c",
                              "Alive"="gray50", "Dead"="red"))+
  scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c"),
                    labels = c("RC1 (n=62)", "RC2 (n=126)",
                               "RC3 (n=113)", "RC4 (n=68)"))+ 
  theme_plots +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10, 
                                  face = "bold"),
          axis.title.y = element_text(size=8))+
  stat_compare_means(method = "kruskal",
                     label.y = ymax*1.8,
                     size = 2.5)+
  stat_pvalue_manual(stat.test, 
                     inherit.aes=F,
                     label="p.adj.signif", 
                     bracket.shorten = .1,
                     tip.length = .02,
                     y.position = ypos, 
                     label.size = size)

# Plots for Tumor Stages 3 and 4
## Statistical test
stat.test <- rstatix::dunn_test(
    data=dat_bp[dat_bp$Tumor.Stage!="Stage 1 and 2",], 
    formula= value ~ Cluster, 
    p.adjust.method = "BH")
stat.test$p.adj <- format(stat.test$p.adj, 
                          scientific = T, 
                          digits = 3)
print(stat.test)
stat.test$size <- ifelse(stat.test$p.adj.signif== "ns", 2, 4)
size <- c()
for(j in 1:nrow(stat.test)){
  size <- c(size, rep(stat.test$size[j], 3))
}
## Plots
g2 <-  
  ggplot(dat_bp[dat_bp$Tumor.Stage!="Stage 1 and 2",],
         aes(x=Cluster, y=value, fill=Cluster))+
  labs(title= "Tumor Stage 3 and 4\nTCGA-ACC (n=30)",
       y= "Median Regulon\nCluster Activity", 
       x= element_blank())+
  geom_violin(
    aes(fill=NULL, col=Cluster),
    alpha=0.4, 
    show.legend = c(T))+
  ggbeeswarm::geom_beeswarm(
    aes(color=Vital.Status), 
    alpha =.5, cex = 2, size=1,
    show.legend = T)+
  geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)),
                     breaks=c(-2,-1,0,1,2))+
  scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c",
                              "Alive"="gray50", "Dead"="red"))+
  scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                             "RC3"="#fb9a99","RC4"="#e31a1c"),
                    labels = c("RC1 (n=62)", "RC2 (n=126)",
                               "RC3 (n=113)", "RC4 (n=68)"))+ 
  theme_plots +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10, 
                                  face = "bold"),
        axis.title.y = element_text(size=8))+
  stat_compare_means(method = "kruskal",
                     label.y = ymax*1.8,
                     size = 2.5)+
  stat_pvalue_manual(stat.test, 
                     inherit.aes=F,
                     label="p.adj.signif", 
                     bracket.shorten = .1,
                     tip.length = .02,
                     y.position = ypos, 
                     label.size = size)
# Assessing legend
leg <- get_legend(
  g1+
    theme(legend.position = "bottom")+
  guides(fill = "none",
         color = guide_legend(nrow=2, title=element_blank(),
                              override.aes = list(alpha=1,
                                                  size=2))))

# Creating panel
row1 <- plot_grid(g1, g2, ncol=2) 
row2 <- plot_grid(leg)

g0 <- plot_grid(row1, row2, nrow=2, rel_heights = c(1,.25))
g0

# Storing in new object to merge with ENSAT results
g10 <- g1
g20 <- g2
leg20 <- leg
g0.stagestcga <- g0

# Rearranging title and labels
g10 <- g10 +
  labs(title= "Tumor Stages 1 and 2\nTCGA-ACC (n=46 samples)",
       y= "Regulon activity (median of the cluster)", 
       x= element_blank())+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 8, 
                                  face = "bold"),
        axis.title.y = element_text(size=6))

g20 <- g20+
  labs(title= "Tumor Stages 3 and 4\nTCGA-ACC (n=30 samples)",
       y= element_blank(), 
       x= element_blank())+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 8, 
                                  face = "bold")) 

```

## MSigDb Hallmarks enrichment

### Creating Hallmarks data frame with categories
```r
# Cleaning environment
rm(gsea2.dif, g, gsea2.dif.prol, rtns, 
   g, g0, g1, g2, g3, g4, leg, leg2, dat_bp, dat_imun,
   dat_imun.cor, dat_imun.cor.df, row1, row2, scatterPlot,
   scatterPlot2, stat.test, cols, cols.immune, dat.rc1,
   dat.rc2, dat.rc3, dat.rc4, form, i, j, k ,lab, size, ymax, ypos)

hallmark <- data.frame(
  stringsAsFactors = FALSE,
  Number = c(1L,2L,3L,4L,5L,6L,7L,8L,
             9L,10L,11L,12L,13L,14L,15L,16L,
             17L,18L,19L,20L,21L,22L,23L,
             24L,25L,26L,27L,28L,29L,30L,31L,
             32L,33L,34L,35L,36L,37L,38L,
             39L,40L,41L,42L,43L,44L,45L,46L,
             47L,48L,49L,50L),
  Hallmark.Name = c("APICAL_JUNCTION",
                    "APICAL_SURFACE", "PEROXISOME",
                    "ADIPOGENESIS", "ANGIOGENESIS",
                    "EPITHELIAL_MESENCHYMAL_TRANSITION",
                    "MYOGENESIS",
                    "SPERMATOGENESIS","PANCREAS_BETA_CELLS",
                    "DNA_REPAIR","UV_RESPONSE_DN",
                    "UV_RESPONSE_UP","ALLOGRAFT_REJECTION",
                    "COAGULATION","COMPLEMENT",
                    "INTERFERON_ALPHA_RESPONSE",
                    "INTERFERON_GAMMA_RESPONSE",
                    "IL6_JAK_STAT3_SIGNALING",
                    "INFLAMMATORY_RESPONSE",
                    "BILE_ACID_METABOLISM",
                    "CHOLESTEROL_HOMEOSTASIS",
                    "FATTY_ACID_METABOLISM", "GLYCOLYSIS",
                    "HEME_METABOLISM", "OXIDATIVE_PHOSPHORYLATION",
                    "XENOBIOTIC_METABOLISM","APOPTOSIS",
                    "HYPOXIA","PROTEIN_SECRETION",
                    "UNFOLDED_PROTEIN_RESPONSE", 
                    "REACTIVE_OXYGEN_SPECIES_PATHWAY",
                    "E2F_TARGETS", "G2M_CHECKPOINT",
                    "MYC_TARGETS_V1", "MYC_TARGETS_V2",
                    "P53_PATHWAY", "MITOTIC_SPINDLE",
                    "ANDROGEN_RESPONSE", "ESTROGEN_RESPONSE_EARLY",
                    "ESTROGEN_RESPONSE_LATE","IL2_STAT5_SIGNALING",
                    "KRAS_SIGNALING_UP","KRAS_SIGNALING_DN",
                    "MTORC1_SIGNALING","NOTCH_SIGNALING",
                    "PI3K_AKT_MTOR_SIGNALING",
                    "HEDGEHOG_SIGNALING","TGF_BETA_SIGNALING",
                    "TNFA_SIGNALING_VIA_NFKB",
                    "WNT_BETA_CATENIN_SIGNALING"),
  Process.Category = c("cellular component",
                       "cellular component","cellular component",
                       "development","development",
                       "development","development","development",
                       "development","DNA damage","DNA damage",
                       "DNA damage","immune","immune",
                       "immune","immune","immune","immune",
                       "immune","metabolic","metabolic",
                       "metabolic","metabolic","metabolic",
                       "metabolic","metabolic","pathway",
                       "pathway","pathway","pathway",
                       "pathway","proliferation", "proliferation",
                       "proliferation","proliferation",
                       "proliferation","proliferation",
                       "signaling","signaling","signaling",
                       "signaling","signaling","signaling",
                       "signaling","signaling","signaling",
                       "signaling","signaling","signaling",
                       "signaling"),
  Description = c(
    "apical junction complex consisting of adherens and tight junctions",
    "membrane proteins in the apical domain","peroxisomes",
    "adipocyte development", "blood vessel formation",
    "epithelial mesenchymal transition",
    "muscle differentiation",
    "sperm development and male fertility",
    "genes specific to pancreatic beta cells",
    "DNA repair","UV response: downregulated genes",
    "UV response: upregulated genes",
    "allograft rejection",
    "blood coagulation cascade","complement cascade",
    "interferon alpha response",
    "interferon gamma response",
    "IL6 STAT3 signaling during acute phase response",                  "inflammation","biosynthesis of bile acids",
    "cholesterol homeostasis",
    "fatty acid metabolism",
    "glycolysis and gluconeogenesis",
    "heme metabolism and erythroid lineage",
    "oxidative phosphorylation and citric acid cycle",
    "metabolism of xenobiotics",
    "programmed cell death; caspase pathway",
    "response to hypoxia; HIF1A targets","protein secretion",  
    "unfolded protein response; ER stress",
    "reactive oxygen species pathway",
    "cell cycle progression: E2F targets",
    "cell cycle progression: G2/M checkpoint",
    "MYC targets, variant 1",
    "MYC targets, variant 2","p53 pathway",
    "cell cycle progression: mitotic spindle assembly",
    "androgen response", "early estrogen response",
    "late estrogen response","IL2 STAT5 signaling",
    "KRAS signaling, upregulated genes",
    "KRAS signaling, downregulated genes",
    "mTORC1 signaling","Notch signaling",
    "PI3K signaling via AKT to mTORC1",
    "Hedgehog signaling","TGF beta signaling",
    "TNFA signaling via NFκB",
    "canonical beta catenin pathway"),
  Number.of.Founder.Sets = c(37L,12L,28L,36L,14L,107L,
                             64L,24L,24L,44L,17L,16L,190L,
                             71L,71L,82L,82L,24L,120L,28L,28L,
                             53L,87L,36L,93L,124L,80L,87L,
                             74L,22L,13L,420L,420L,404L,6L,
                             85L,108L,8L,61L,61L,13L,14L,16L,
                             487L,49L,591L,79L,29L,132L,49L),
  Number.of.Genes = c(200L,44L,107L,200L,36L,
                      200L,200L,135L,40L,150L,144L,158L,
                      200L,138L,200L,97L,200L,87L,200L,
                      112L,74L,158L,200L,200L,200L,
                      200L,161L,200L,96L,113L,49L,200L,
                      200L,200L,58L,200L,200L,117L,
                      200L,200L,200L,200L,200L,200L,32L,
                      105L,36L,54L,200L,42L))

```

### Enrichment analysis for regulon activity - Code snippet for Figure 4A and Supplementary Table 7 (TCGA)

```r

# For hallmarks with RTN pipeline
# Assessing Hallmarks gene sets
h_gene_sets = msigdbr(species = "Homo sapiens", 
                      category = "H")

h_gene_sets = split(x = h_gene_sets$gene_symbol, 
                    f = h_gene_sets$gs_name)
# Annotate regulons
tni.annotate.Hallmark <- 
  tni.annotate.regulons(rtni_tcgaACC,
                        regulatoryElements = regs.risk$SYMBOL,
                        geneSetList = h_gene_sets)
# Retrieving results as data frame
tni.annotate.Hallmark <- as.data.frame(t(tni.annotate.Hallmark))
rownames(tni.annotate.Hallmark) <- 
  gsub("HALLMARK_", "",
       rownames(tni.annotate.Hallmark))
# Exporting results 
writexl::write_xlsx(cbind(
   Hallmark = rownames(tni.annotate.Hallmark),
   tni.annotate.Hallmark),
   path="./TNI_Annotate_Hallmark.xlsx")

# Preparing annotations for Heatmap
colannot <- rowannot[colnames(tni.annotate.Hallmark),]
rownames(hallmark) <- hallmark$Hallmark.Name
hallmark <- 
  hallmark[order(hallmark$Process.Category,
                 hallmark$Hallmark.Name),]

tni.annotate.Hallmark <- tni.annotate.Hallmark[rownames(hallmark),]

rowannot <- hallmark[rownames(tni.annotate.Hallmark),
                     "Process.Category", drop=F]
rowannot$Process.Category <- as.factor(rowannot$Process.Category)

# Color pallete
pal <- c("#003C30", "#01665E", "#35978F", "#80CDC1", "#C7EAE5", 
         "#ffffff",
         "#F6E8C3", "#DFC27D", "#BF812D", "#8C510A", "#543005")

# Right annotation - Hallmarks category
hallmark.annot <- rowannot$Process.Category
names(hallmark.annot) <- rownames(rowannot)
# Top annotation (regulon characteristics)
reg.annot <- colannot
reg.annot <- reg.annot[colnames(tni.annotate.Hallmark),]
reg.annot$LF.rho <- ifelse(reg.annot$LF.rho.padj < 0.05,
                           reg.annot$LF.rho, 0)
# Splitting by Hallmarks category
split <- as.factor(rowannot$Process.Category)

tni.annotate.Hallmark <- tni.annotate.Hallmark[,rownames(reg.annot)]
identical(colnames(tni.annotate.Hallmark), rownames(reg.annot)) #T
# Clustering enrichment results within regulon clusters
RC1 <- dist(t(
  tni.annotate.Hallmark[, rownames(
    reg.annot[reg.annot$Cluster=="RC1",])]),
  method = "euclidean")
RC1 <- hclust(RC1, method = "ward.D2")
RC1 <- RC1$labels[RC1$order]

RC2 <- dist(t(
  tni.annotate.Hallmark[,rownames(
    reg.annot[reg.annot$Cluster=="RC2",])]),
  method = "euclidean")
RC2 <- hclust(RC2, method = "ward.D2")
RC2 <- rev(RC2$labels[RC2$order])

RC3 <- dist(t(
  tni.annotate.Hallmark[,rownames(
    reg.annot[reg.annot$Cluster=="RC3",])]),
  method = "euclidean")
RC3 <- hclust(RC3, method = "ward.D2")
RC3 <- RC3$labels[RC3$order]

RC4 <- dist(t(
  tni.annotate.Hallmark[,rownames(
    reg.annot[reg.annot$Cluster=="RC4",])]),
  method = "euclidean")
RC4 <- hclust(RC4, method = "ward.D2")
RC4 <- RC4$labels[RC4$order]
# Ordering data frame
tni.annotate.Hallmark <- tni.annotate.Hallmark[,c(RC1,RC2,RC3,RC4)]
reg.annot <- reg.annot[colnames(tni.annotate.Hallmark),]
# Plotting heatmap
h1<- Heatmap(
  matrix = tni.annotate.Hallmark,
  name = "Hallmark GSEA",
  show_column_names = F,
  cluster_rows = T,
  show_row_names = T, 
  use_raster = F,
  split = split,
  show_column_dend = T,
  cluster_columns = F,
  column_split = as.factor(reg.annot$Cluster),
  row_title_gp = gpar(fontsize=8),
  column_title = "369 prognostic regulons",
  column_title_rot = 0,
  column_title_gp = gpar(fontsize=8, fontface="bold"),
  column_title_side = "top",
  show_row_dend = F,  
  row_names_gp = gpar(fontsize=6), 
  row_title_rot = 0,
  column_names_side = "bottom", 
  column_names_centered = T,
  column_names_rot = 90,
  heatmap_legend_param = list(
    direction = "horizontal", 
    title="Hallmark\nGSEA",
    title_gp = gpar(fontsize=6, fontface="bold"),
    labels_gp=gpar(fontsize=6),
    at=c(-2, 0, 2),
    title_position="topcenter"),
  top_annotation = columnAnnotation(
    df=reg.annot[c("Cluster", "HR", "GSEA2.Ster",
                  "GSEA2.Prol", "LF.rho")],
    annotation_name_side="left", 
    annotation_name_gp=gpar(fontsize=8),
    annotation_label=c("Regulon Cluster", 
                       "Hazard Ratio", "Steroid IP score",
                       "Proliferation IS score",
                       "Leuk.Fraction rho"),
    simple_anno_size=unit(3,"mm"),
    show_legend=T,
    annotation_legend_param= list(
      Cluster=list(
        title="Regulon\nCluster",
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        nrow=2,
        direction="horizontal"),
      HR=list(
        title="Hazard\nRatio",
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        direction="horizontal"),
      LF.rho=list(
        title="Leuk. Fraction\nCorrelation",
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        direction="horizontal",
        at=c(-1, 0, 1),
        legend_width= unit(1.5,"cm")),
      GSEA2.Ster=list(
        title="Steroid\nIP Score", 
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        direction="horizontal",
        at=c(-2, 0, 2),
        legend_width= unit(1.5,"cm")),
      GSEA2.Prol=list(
        title="Proliferation\nIS Score",
        title_position="topcenter",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        direction="horizontal",
        at=c(-2, 0, 2),
        legend_width= unit(1.5,"cm"))),
    col= list(
      Cluster = c("RC1" = "#a6cee3",
                  "RC2" = "#1f78b4",
                  "RC3" = "#fb9a99",
                  "RC4" = "#e31a1c"),
      HR = c("<1" = "#1B7837",
             ">1" = "#BF812D"),
      LF.rho = 
        circlize::colorRamp2(c(-.75,0,.75),
                             c("darkorange", "white", "blue")),
      GSEA2.Ster = 
        circlize::colorRamp2(c(-2,0,2),
                             c("#00BFC4","white","#F8766D")),
      GSEA2.Prol = 
        circlize::colorRamp2(c(-2,0,2),
                             c( "chartreuse4", "white", "red")))),
  right_annotation = rowAnnotation(
    df=hallmark.annot,
    simple_anno_size=unit(3,"mm"), 
    annotation_name_gp=gpar(fontsize=0),
    show_legend=T,
    annotation_legend_param= list(
      df=list(
        title="Hallmark Category",
        title_gp = gpar(fontsize=6, fontface="bold"),
        labels_gp=gpar(fontsize=6),
        nrow = 4,
        direction="horizontal")),
     col= list(
       df=c(immune="#7FC97F", development= "#BEAED4",
            signaling="#FDC086", 'cellular component' ="#FFFF99",
            metabolic="#BF5B17", 'DNA damage'= "#F0027F",
            proliferation="#386CB0", pathway= "#666666"))))
    
ht <- draw(h1, heatmap_legend_side="bottom", merge_legends=T)

g0 <- grid.grabExpr(draw(ht,
                         heatmap_legend_side="bottom", 
                         annotation_legend_side = "bottom",
                         merge_legends=T),
                   height = 7, width = 5.5)
plot_grid(g0)
# Saving
ggsave(g0, filename = "./6 Heatmap_Hallmark_TCGA.pdf",
        device="pdf",
        width = 6, height = 8, units="in")
```

## Code snippet for Figure 4B, C, and D
```r
# Hallmarks selected
H <- c("INFLAMMATORY_RESPONSE", 
       "WNT_BETA_CATENIN_SIGNALING",
       "PI3K_AKT_MTOR_SIGNALING")

## Preparing data.frame for boxplot
identical(rownames(reg.annot), colnames(tni.annotate.Hallmark)) #T
dat <- as.data.frame(cbind(
  Cluster=reg.annot$Cluster, 
  t(tni.annotate.Hallmark[H,])))
dat[,2:ncol(dat)] <- lapply(dat[,2:ncol(dat)] , as.numeric)
cols <- c("Inflammatory Response",
          "WNT Beta Catenin Signaling",
          "PI3K Akt mTOR Signaling")
# Looping through the 3 hallmarks and storing plots
g <- list()
for(i in 2:ncol(dat)){
  form <- as.formula(paste0(names(dat)[i], "~", "Cluster"))
  stat.test <- rstatix::dunn_test(data=dat, 
                                  formula=form, 
                                  p.adjust.method = "BH")
  stat.test$p.adj <- format(stat.test$p.adj, 
                            scientific = T, 
                            digits = 3)
  stat.test$size <- ifelse(stat.test$p.adj.signif== "ns", 2, 4)
  size <- c()
  for(j in 1:nrow(stat.test)){
    size <- c(size, rep(stat.test$size[j], 3))
  }
  lab <- cols[i-1]
  ymax <- round(max(dat[,i]))
  ypos <- c(ymax, ymax*1.2, ymax*1.6, ymax, ymax*1.4, ymax)
  g1<- eval(substitute(
    ggplot(dat,
          aes(x=Cluster, y= dat[,i], fill=Cluster))+
    geom_boxplot(width=c(0.5), lwd=0.2, outlier.size = 0.5)+  
    geom_violin(
      aes(fill=NULL, col=Cluster),
      alpha=0.4, 
      show.legend = F)+
    geom_jitter(size=0.2, width = 0.1)+
    geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
    scale_y_continuous(expand = expansion(mult = c(0.1,0.1)),
                       breaks=c(-ymax,-ymax/2,0,ymax/2,ymax))+
    scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                                "RC3"="#fb9a99","RC4"="#e31a1c"))+
    scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                               "RC3"="#fb9a99","RC4"="#e31a1c"),
                      labels = c("RC1 (n=62)", "RC2 (n=126)",
                                 "RC3 (n=113)", "RC4 (n=68)"))+ 
    labs(y= lab, 
         x= element_blank())+
    theme_plots+
    stat_compare_means(method = "kruskal",
                       label.y = ymax*1.8,
                       size = 2.5)+
    stat_pvalue_manual(stat.test, 
                       inherit.aes=F,
                       label="p.adj.signif", 
                       bracket.shorten = .1,
                       tip.length = .02,
                       y.position = ypos, 
                       label.size = size)
    ,list(i=i)))
  g[[i-1]] <- g1
}
# Assessing the legend
leg <- get_legend(
  g1+
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow=2,
                             title= "Regulon cluster", 
                             title.hjust = .5)))
# Mounting panel
g0 <- plot_grid(g[[2]], g[[1]], g[[3]], leg,
                ncol = 1, rel_heights = c(1,1,1,.25),
                labels = c("B", "C", "D", ""))
g0
# Saving
ggsave(g0, 
        filename = "./6_Main_Boxplot_hallmarks_v4.pdf", 
        device="pdf", 
        width = 2.5, height = 7)
```

## Code snippet for Supplementary Figure 3
```r
# Creating boxplot for all Hallmarks
identical(colnames(tni.annotate.Hallmark), rownames(reg.annot)) #T
# Separating Hallmarks by category
levs <- levels(split)
levs <- c("immune", "proliferation", "signaling",
          "pathway", "metabolic", "DNA damage", 
          "cellular component", "development")
# Function to make first letter uppercase
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Looping through Hallmark categories
for(i in 1:length(levs)){
  dat_bp <- cbind(Cluster=reg.annot$Cluster, 
                  t(tni.annotate.Hallmark[split %in% levs[i],]))
  dat_bp <- as.data.frame(dat_bp)
  dat_bp[,2:ncol(dat_bp)] <- lapply(dat_bp[,2:ncol(dat_bp)] , as.numeric)
  ## Preparing data.frame for boxplot
  dat_bp <- pivot_longer(dat_bp,
                         cols = -contains(c("Cluster")),
                         values_to = "value",
                         names_to = "variable")
  dat_bp$variable <- factor(dat_bp$variable)
  dat_bp$Cluster <- factor(dat_bp$Cluster)
  dat_bp$value <- as.numeric(dat_bp$value)
  ## Plotting
  g1<- ggplot(dat_bp,
              aes(x=Cluster, y= value, fill=Cluster))+
    facet_wrap(.~variable, ncol=5)+
    geom_boxplot(width=c(0.4), lwd=0.2, outlier.size = 0.5)+  
    geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
    geom_violin(
      aes(fill=NULL, col=Cluster),
      alpha=0.4, 
      show.legend = F)+
    scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                                "RC3"="#fb9a99","RC4"="#e31a1c"),
                       labels = c("RC1 (n=62)", "RC2 (n=126)",
                                  "RC3 (n=113)", "RC4 (n=68)"))+
    labs(y="Regulon Hallmark Score", x= element_blank(), 
         title = paste0(firstup(levs[i]), " Hallmarks")) +
    theme_plots+
    theme(strip.text = element_text(size=6))+
    stat_compare_means(
      label.x = 1.5, 
      label.y = 2.5,
      size=2)+
    scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                               "RC3"="#fb9a99","RC4"="#e31a1c"),
                      labels = c("RC1 (n=62)", "RC2 (n=126)",
                                 "RC3 (n=113)", "RC4 (n=68)"))+
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
    guides(fill=guide_legend(nrow = 2, title = "Regulon Cluster"))
  # definying width and height of panel depending on the number of hallmarks in each category
  if(length(unique(dat_bp$variable)) < 5 ){ 
    width=5
    height=3
  }else{
      width=8
      if(length(unique(dat_bp$variable)) < 10 ){
        height=4.8
        if(length(unique(dat_bp$variable)) == 5 ){height=3}
      }else{
        height=7
      }
  }
  g1
  # saving
  ggsave(plot_grid(g1, labels = LETTERS[i]), 
          filename = paste0("./Sup_3",LETTERS[i],
                            "_boxplot_hallmarks_", levs[i],".pdf"),
          device="pdf",
          width = width, height=height)
}

# cleaning environment
rm(g, h_gene_sets, h1, tni.annotate.Hallmark, hallmark.annot, pal, split, i, H, ht, g1,  clusters, g0, hallmark, leg, rowannot)

```

## RTNDuals pipeline - Code snippet for Supplementary Table 8

```r
# RTNDuals
# Inferring targets sharing between regulons
rmbr <- tni2mbrPreprocess(rtni_tcgaACC, 
                          regulatoryElements = regs.risk$SYMBOL)
rmbr <- mbrAssociation(rmbr, pValueCutoff = 0.05)
# Saving data
save(rmbr, file="./rmbr.RData")
mbrGet(rmbr, what="summary")
overlap <- mbrGet(rmbr, what="dualsOverlap")
correlation <- mbrGet(rmbr, what="dualsCorrelation") #150 duals
# Checking regulons with higher number of associations
summary(as.factor(
  correlation$Regulon1))[order(summary(as.factor(correlation$Regulon1)), 
                               decreasing = T)]
# CENPA   LIN54   DEAF1    DPF1    ETV4   FOXK1   DOT1L    E4F1    MITF  SETDB1 CREB3L2 
#       9       7       5       5       5       5       4       4       4       4       3 
# Big group centered in CENPA (HR > 1) and LIN54  (HR < 1)

# Checking regulons with higher number of associations
summary(as.factor(
  correlation$Regulon2))[order(summary(as.factor(correlation$Regulon2)), 
                               decreasing = T)]
# SNAPC4  THAP10    MITF   MYBL2    TCF3  ZNF496  ZNF770  ZNF787    KLF6    MXD3   PROX2 
#       8       6       4       4       4       4       4       4       3       3       3 

# Exporting correlation table
writexl::write_xlsx(correlation, path="./Duals_correlation.xlsx")

# Cleaning environment
rm(correlation,  overlap, rmbr, colannot)

```

#  ENSAT (GSE49280) 

## Computing regulon activity - Code snippet for Supplementary table 5
```r
# Load ENSAT data prepared in Preprocessing.R
## Download data

## Loading data
load("./ENSAT.RData")
gexp_ENSAT <- ENSAT[[1]]
ENSAT_coldat <- ENSAT[[2]]
rm(ENSAT)

# Replace samples
rtni_ENSAT <- tni.replace.samples(rtni_tcgaACC, gexp_ENSAT)
#Warning in tni.replace.samples(rtni_tcgaACC, gexp_ENSAT) :
#  NOTE: 25.7046176701193% of the TNI is not represented in new 'rowAnnotation'!
# TFs found in TCGA-ACC that were not present in the ENSAT gene expression matrix

# Compute regulon activity for the new samples
rtni_ENSAT <- tni.gsea2(rtni_ENSAT)

regact_ENSAT <- tni.get(rtni_ENSAT, 
                        what = "regulonActivity",
                        idkey = "SYMBOL")
regact_ENSAT <- as.data.frame(regact_ENSAT$differential)

# Exporting regulon activity results
writexl::write_xlsx(regact_ENSAT, path = "./Regact_ENSAT.xlsx")
```

## RTNSurvival pipeline - Computing survival curves for ENSAT cohort - Code snippet for Figure 7B and Supplementary table 3 (ENSAT)

```r
# Dummy variables for survival analysis
ENSAT_coldat$C1A <- ifelse(ENSAT_coldat$Consensus.clustering.K2=="C1A", 1, 0)
ENSAT_coldat$C1B <- ifelse(ENSAT_coldat$Consensus.clustering.K2=="C1B",1,0)
ENSAT_coldat$CIMP.High <- ifelse(ENSAT_coldat$CIMP=="CIMP.high", 1, 0)
ENSAT_coldat$CIMP.Low <- ifelse(ENSAT_coldat$CIMP=="CIMP.low", 1, 0)
ENSAT_coldat$CIMP.High[is.na(ENSAT_coldat$CIMP.High)] <- 0
ENSAT_coldat$CIMP.Low[is.na(ENSAT_coldat$CIMP.Low)] <- 0

ENSAT_annot_cox <- ENSAT_coldat
# Covariates for Multivariate Cox
colnames(ENSAT_annot_cox)[8] <- "Age"
colnames(ENSAT_annot_cox)[13] <- "Tumor.Stage"
# Filtering for samples with Tumor Stage annotation
ENSAT_annot_cox$Tumor.Stage <- ifelse(
  complete.cases(ENSAT_annot_cox$Tumor.Stage), ENSAT_annot_cox$Tumor.Stage, 0)
# Dummy variables
ENSAT_annot_cox$Stage1 <- ifelse(ENSAT_annot_cox$Tumor.Stage==1,1,0)
ENSAT_annot_cox$Stage2 <- ifelse(ENSAT_annot_cox$Tumor.Stage==2,1,0)
ENSAT_annot_cox$Stage3 <- ifelse(ENSAT_annot_cox$Tumor.Stage==3,1,0)
ENSAT_annot_cox$Stage4 <- ifelse(ENSAT_annot_cox$Tumor.Stage==4,1,0)

# RTNSurvival pipeline
rtns <- tni2tnsPreprocess(
  rtni_ENSAT,
  regulatoryElements = 
    colnames(regact_ENSAT)[colnames(regact_ENSAT) %in% regs.risk$SYMBOL],
  survivalData = ENSAT_annot_cox,
  time = 22,
  event = 21,
  keycovar = c("Age", "Tumor.Stage"),
  excludeMid = F)
rtns <- tnsGSEA2(rtns)

### selecting regulons for Cox multivariate analysis
rtns <- tnsCox(rtns)
reg.cox <- rtns@results[["Cox"]][["Table"]]
# Exporting results
writexl::write_xlsx(reg.cox, path="./ENSAT/Cox_ENSAT_v4.xlsx")

# Forest plot
tnsPlotCox(rtns, regs = c("NR5A1", "CENPA"), sortregs = F,
           xlim = c(0.1,15),
           plotpdf = F,#T, 
           fname= "Cox_ENSAT", 
           fpath = "./", width = 4, height = 2)

# Kaplan Meier Analysis
rtns <- tnsKM(rtns, 
              regs = c("NR5A1", "CENPA"))
reg.KM <- rtns@results[["KM"]][["Table"]]

# Filtering the covariates, maintaining only regulons
reg.cox <- reg.cox[!rownames(reg.cox) %in% c("Age", "Tumor.Stage"),]
dim(reg.cox[reg.cox$Adjusted.Pvalue<0.05,]) # 328 of 361

# Kaplan Meier plot
atribs <- c("C1A", "C1B", "Stage1", "Stage2", "Stage3", "Stage4")
tnsPlotKM(rtns, regs=c("NR5A1"), attribs = list(atribs),
          plotpdf = F,#T,
          fname = "ENSAT_KM")
tnsPlotKM(rtns, regs=c("CENPA"), attribs = list(atribs),
          plotpdf = F,#T, 
          fname = "ENSAT_KM")


```
## Code snippet for Supplementary Figure 5

```r
# Selecting activity for regulons in the TCGA prognostic list
regact_ENSAT <- regact_ENSAT[,colnames(regact_ENSAT) %in% regs.risk$SYMBOL] 
reg.cox <- reg.cox[rownames(reg.cox) %in% c(regs.risk$SYMBOL),]
# Preparing right annotation (from TCGA)
rowannot <- reg.annot[reg.annot$SYMBOL %in% colnames(regact_ENSAT),]
rownames(rowannot) <- rowannot$SYMBOL
rowannot <- rowannot[colnames(regact_ENSAT),]

identical(rownames(rowannot), colnames(regact_ENSAT)) #T
# Preparing top (columns) annotation
colannot <- ENSAT_coldat
identical(rownames(colannot), rownames(regact_ENSAT)) #T

rowannot.ensat <- reg.cox[rownames(reg.cox) %in% colnames(regact_ENSAT),]
rowannot.ensat <- rowannot.ensat[colnames(regact_ENSAT),]

identical(colnames(regact_ENSAT), rownames(rowannot.ensat)) #T

rowannot.ensat$Cox.HR <- ifelse(rowannot.ensat$Adjusted.Pvalue < 0.05,
                                ifelse(rowannot.ensat$HR < 1, "<1", ">1"), 0)

colannot$ENSAT.staging[is.na(colannot$ENSAT.staging)] <- "NA"

identical(colnames(regact_ENSAT), rownames(rowannot)) #T

# Semi supervised clustering within regulon clusters
RC1 <- dist(t(
  regact_ENSAT[, rownames(
    rowannot[rowannot$Cluster=="RC1",])]),
  method = "euclidean")
RC1 <- hclust(RC1, method = "ward.D2")
RC1 <- RC1$labels[RC1$order]

RC2 <- dist(t(
  regact_ENSAT[,rownames(
    rowannot[rowannot$Cluster=="RC2",])]),
  method = "euclidean")
RC2 <- hclust(RC2, method = "ward.D2")
RC2 <- rev(RC2$labels[RC2$order])

RC3 <- dist(t(
  regact_ENSAT[,rownames(
    rowannot[rowannot$Cluster=="RC3",])]),
  method = "euclidean")
RC3 <- hclust(RC3, method = "ward.D2")
RC3 <- RC3$labels[RC3$order]

RC4 <- dist(t(
  regact_ENSAT[,rownames(
    rowannot[rowannot$Cluster=="RC4",])]),
  method = "euclidean")
RC4 <- hclust(RC4, method = "ward.D2")
RC4 <- RC4$labels[RC4$order]
# Ordering data frame and annotations
regact_ENSAT <- regact_ENSAT[,c(RC1,RC2,RC3,RC4)]
rowannot <- rowannot[colnames(regact_ENSAT),]
rowannot.ensat <- rowannot.ensat[rownames(rowannot),]

identical(colnames(regact_ENSAT), rownames(rowannot)) #T
identical(colnames(regact_ENSAT), rownames(rowannot.ensat)) #T

# Plotting heatmap
h1<- Heatmap(
  matrix = t(regact_ENSAT),
  name = "Regulon Activity",
  show_column_names = F,
  row_split = as.factor(rowannot$Cluster),
  row_gap = unit(0.3,"mm"),
  cluster_rows = F,
  show_row_names = F, 
  use_raster = F,
  show_column_dend = T,
  clustering_distance_columns = "maximum",
  clustering_method_columns = "ward.D2",
  column_title = "ENSAT cohort (n=44)",
  column_title_rot = 0,
  column_title_gp = gpar(fontsize=8, fontface="bold"),
  column_title_side = "top",
  show_row_dend = T,  
  row_title_gp = gpar(fontsize=8, fontface="bold"),
  row_title_rot = 90,
  row_title = "361 Prognostic Regulons",
  column_names_gp = gpar(fontsize=8), 
  column_names_side = "top", 
  column_names_centered = T,
  column_names_rot = 90,
  na_col = "white",
  col = circlize::colorRamp2(
    breaks = seq(-2,2, length.out=7),
    colors=rev(RColorBrewer::brewer.pal(7,"RdBu"))),
  heatmap_legend_param = list(
    direction = "horizontal", 
    at=c(-2, 0, 2),
    title = "Regulon\nActivity",
    title_gp = gpar(fontsize=6, fontface="bold"),
    labels_gp=gpar(fontsize=6),
    title_position="topcenter"),
  top_annotation = columnAnnotation(
    df=colannot[,c("Consensus.clustering.K2", "ENSAT.staging")],
    annotation_name_side="left", 
    annotation_name_gp=gpar(fontsize=8),
    annotation_label=c("C1A/C1B", "Tumor Stage"),
    simple_anno_size=unit(3,"mm"),
    show_legend=c(T),
    na_col = "white",
    annotation_legend_param= list(
      Consensus.clustering.K2=list(
       title="C1A/C1B", 
       title_position="topcenter",
       title_gp = gpar(fontsize=6, fontface="bold"),
       labels_gp=gpar(fontsize=6),
       border="black",
       direction="vertical"),
      ENSAT.staging=list(
       title="Tumor Stage", 
       title_position="topcenter",
       title_gp = gpar(fontsize=6, fontface="bold"),
       labels_gp=gpar(fontsize=6),
       border="black",
       nrow = 2,
       direction="vertical")),
    col= list(
      ENSAT.staging=c("1"="gray90",
                      "2"="gray66",
                      "3"="gray33",
                      "4"="black",
                      "NA" = "white"),
      Consensus.clustering.K2= c("C1A" ="black" ,
                                 "C1B"= "gray90"))),
  left_annotation = rowAnnotation(
    df=rowannot.ensat[,c("Cox.HR"), drop=F],
    annotation_label="ENSAT HR", 
    annotation_name_gp=gpar(fontsize=6),
    annotation_name_rot=45,
    simple_anno_size=unit(2,"mm"),
    show_legend=T,
    na_col = "white",
    annotation_legend_param= list(
        Cox.HR=list(
          title="Hazard\nRatio",
          title_position="topcenter",
          title_gp = gpar(fontsize=6, fontface="bold"),
          labels = c("<1", ">1", "NS"),
          labels_gp=gpar(fontsize=6),
          border = "black",
          direction="horizontal")),
    col= list(
      Cox.HR = c("<1"= "#35978F",
                 "0" = "white",
                 ">1"="#BF812D"))),
  right_annotation = rowAnnotation(
    df=rowannot[,c("HR", "GSEA2.Ster", "GSEA2.Prol",
                   "Leukocyte.Fraction.rho", "Cluster"), drop=F],
    annotation_label = c("TCGA HR", "Steroid IP",
                         "Prolif. IS", "LF rho", "Cluster"),
    annotation_name_gp=gpar(fontsize=6),
    annotation_name_rot=45,
    gap = unit(1,"mm"),
    simple_anno_size=unit(2,"mm"),
    show_legend=c(F,T,T,T,T),
    annotation_legend_param= list(
        HR=list(
          title="Hazard\nRatio",
          title_position="topcenter",
          title_gp = gpar(fontsize=6, fontface="bold"),
          labels_gp=gpar(fontsize=6),
          direction="horizontal"),
        Leukocyte.Fraction.rho=list(
          title="Leuk. Fraction\nCorrelation",
          title_position="topcenter",
          title_gp = gpar(fontsize=6, fontface="bold"),
          labels_gp=gpar(fontsize=6),
          direction="horizontal",
          at=c(-1, 0, 1),
          legend_width= unit(1.5,"cm")),
        Cluster=list(
          title="Regulon\nCluster",
          title_position="topcenter",
          nrow=2,
          title_gp = gpar(fontsize=6, fontface="bold"),
          labels_gp=gpar(fontsize=6),
          direction="horizontal"),
        GSEA2.Ster=list(
          title="Steroid\nIP Score", 
          title_position="topcenter",
          title_gp = gpar(fontsize=6, fontface="bold"),
          labels_gp=gpar(fontsize=6),
          direction="horizontal",
          at=c(-2, 0, 2),
          legend_width= unit(1.5,"cm")),
        GSEA2.Prol=list(
          title="Proliferation\nIS Score",
          title_position="topcenter",
          title_gp = gpar(fontsize=6, fontface="bold"),
          labels_gp=gpar(fontsize=6),
          direction="horizontal",
          at=c(-2, 0, 2),
          legend_width= unit(1.5,"cm"))),
    col= list(
      HR = c("<1" = "#1B7837",
             ">1" = "#BF812D"),
      Leukocyte.Fraction.rho = 
        circlize::colorRamp2(c(-.75,0,.75),
                             c("darkorange", "white", "blue")),
      Cluster = c("RC1" = "#a6cee3",
                  "RC2" = "#1f78b4",
                  "RC3" = "#fb9a99",
                  "RC4" = "#e31a1c"),
      GSEA2.Ster = 
        circlize::colorRamp2(c(-2,0,2),
                             c("#00BFC4","white","#F8766D")),
      GSEA2.Prol = 
        circlize::colorRamp2(c(-2,0,2),
                             c( "chartreuse4", "white", "red")))))

ht <- draw(h1, heatmap_legend_side="bottom", merge_legends=T)

g0 <- grid.grabExpr(draw(ht,
                         heatmap_legend_side="bottom", 
                         annotation_legend_side = "bottom",
                         merge_legends=T),
                   height = 7, width = 5)
plot_grid(g0)
# Saving
ggsave(g0, 
      filename = "./ENSAT/Heatmap_regact_ENSAT.pdf", 
      device="pdf", width = 5, height = 8.5, units="in")

# Cleaning environment
rm( h1, dup, rowdat.up_ENSAT, atribs)
```

## Code snippet for Figure 5D

```r
# Calculating median acitivity for each regulon cluster
dat <- as.data.frame(t(regact_ENSAT))

dat.rc1 <- dat[c1[c1 %in% rownames(dat)],]
dat.rc1 <- apply(dat.rc1, 2, median)

dat.rc2 <- dat[c2[c2 %in% rownames(dat)],]
dat.rc2 <- apply(dat.rc2, 2, median)

dat.rc3 <- dat[c3[c3 %in% rownames(dat)],]
dat.rc3 <- apply(dat.rc3, 2, median)

dat.rc4 <- dat[c4[c4 %in% rownames(dat)],]
dat.rc4 <- apply(dat.rc4, 2, median)

identical(rownames(dat.rc1), rownames(dat.rc2)) #T
identical(rownames(dat.rc1), rownames(dat.rc3)) #T
identical(rownames(dat.rc3), rownames(dat.rc4)) #T

dat <- as.data.frame(
  cbind(RC1 = dat.rc1,
        RC2 = dat.rc2,
        RC3 = dat.rc3,
        RC4 = dat.rc4))

identical(rownames(dat), rownames(ENSAT_coldat))#T
# Merging with clinical data
dat <- cbind(dat,
             Tumor.Stage = ENSAT_coldat$ENSAT.staging,
             Vital.Status = ENSAT_coldat$OS)

dat <- dat[complete.cases(dat$Tumor.Stage),]
# Grouping low and high tumor stages
dat$Tumor.Stage <- as.factor(ifelse(dat$Tumor.Stage %in% c(1,2),
                                    "Stage 1 and 2",
                                    "Stage 3 and 4"))
dat$Vital.Status <- as.factor((ifelse(dat$Vital.Status == 0,
                                      "Alive", "Dead")))
# Preparing data frame for boxplot
dat_bp <- pivot_longer(dat,
                       cols = -contains(c("Tumor.Stage",
                                          "Vital.Status")),
                       values_to = "value",
                       names_to = "Cluster")

# Tumor stages 1 and 2
## Statistical test
stat.test <- rstatix::dunn_test(
    data=dat_bp[dat_bp$Tumor.Stage=="Stage 1 and 2",], 
    formula= value ~ Cluster, 
    p.adjust.method = "BH")
stat.test$p.adj <- format(stat.test$p.adj, 
                          scientific = T, 
                          digits = 3)
print(stat.test)
stat.test$size <- ifelse(stat.test$p.adj.signif== "ns", 2, 4)
size <- c()
for(j in 1:nrow(stat.test)){
  size <- c(size, rep(stat.test$size[j], 3))
}
ymax <- round(max(dat_bp$value))
ymax <- 1.8
ypos <- c(ymax, ymax*1.2, ymax*1.6, ymax, ymax*1.4, ymax)
## Plotting
g1<-  
  ggplot(dat_bp[dat_bp$Tumor.Stage=="Stage 1 and 2",],
         aes(x=Cluster, y=value, fill=Cluster))+
  labs(title= "Tumor Stage 1 and 2\nENSAT (n=28)",
       y= "Median Regulon\nCluster Activity", 
       x= element_blank()) +
  geom_violin(
    aes(fill=NULL, col=Cluster),
    alpha=0.4, 
    show.legend = c(T))+
  ggbeeswarm::geom_beeswarm(
    aes(color=Vital.Status),
    alpha =.5, 
    cex = 2, size=1, show.legend = T)+
  geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)),
                     breaks=c(-2,-1,0,1,2))+
  scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c",
                              "Alive"="gray50", "Dead"="red"))+
  scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c"),
                    labels = c("RC1 (n=62)", "RC2 (n=126)",
                               "RC3 (n=113)", "RC4 (n=68)"))+ 
  theme_plots +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10, 
                                  face = "bold"),
          axis.title.y = element_text(size=8))+
  stat_compare_means(method = "kruskal",
                     label.y = ymax*1.8,
                     size = 2.5)+
  stat_pvalue_manual(stat.test, 
                     inherit.aes=F,
                     label="p.adj.signif", 
                     bracket.shorten = .1,
                     tip.length = .02,
                     y.position = ypos, 
                     label.size = size)

# Tumor Stages 3 and 4
## Statistical test
stat.test <- rstatix::dunn_test(
    data=dat_bp[dat_bp$Tumor.Stage!="Stage 1 and 2",], 
    formula= value ~ Cluster, 
    p.adjust.method = "BH")
stat.test$p.adj <- format(stat.test$p.adj, 
                          scientific = T, 
                          digits = 3)
print(stat.test)
stat.test$size <- ifelse(stat.test$p.adj.signif== "ns", 2, 4)
size <- c()
for(j in 1:nrow(stat.test)){
  size <- c(size, rep(stat.test$size[j], 3))
}
## Plotting
g2 <-  
  ggplot(dat_bp[dat_bp$Tumor.Stage!="Stage 1 and 2",],
         aes(x=Cluster, y=value, fill=Cluster))+
  labs(title= "Tumor Stage 3 and 4\nENSAT (n=15)",
       y= "Median Regulon\nCluster Activity", 
       x= element_blank())+
  geom_violin(
    aes(fill=NULL, col=Cluster),
    alpha=0.4, 
    show.legend = c(T))+
  ggbeeswarm::geom_beeswarm(
    aes(color=Vital.Status), 
    alpha =.5, cex = 2, size=1,
    show.legend = T)+
  geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.15)),
                     breaks=c(-2,-1,0,1,2))+
  scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c",
                              "Alive"="gray50", "Dead"="red"))+
  scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                             "RC3"="#fb9a99","RC4"="#e31a1c"),
                    labels = c("RC1 (n=62)", "RC2 (n=126)",
                               "RC3 (n=113)", "RC4 (n=68)"))+ 
  theme_plots +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10, 
                                  face = "bold"),
        axis.title.y = element_text(size=8))+
  stat_compare_means(method = "kruskal",
                     label.y = ymax*1.8,
                     size = 2.5)+
  stat_pvalue_manual(stat.test, 
                     inherit.aes=F,
                     label="p.adj.signif", 
                     bracket.shorten = .1,
                     tip.length = .02,
                     y.position = ypos, 
                     label.size = size)
# Assessing legend
leg <- get_legend(
  g1+
    theme(legend.position = "bottom")+
  guides(fill = "none",
         color = guide_legend(nrow=2, title=element_blank(),
                              override.aes = list(alpha=1,
                                                  size=2))))
# Storing
g10.ensat <- g1
g20.ensat <- g2
leg.tumorstage <- leg

# Changing labels
g10.ensat <- g10.ensat +
  labs(title= "Tumor Stages 1 and 2\nENSAT (n=28 samples)",
       y= element_blank(), 
       x= element_blank())+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 8, 
                                  face = "bold"))     

g20.ensat <- g20.ensat +
  labs(title= "Tumor Stages 3 and 4\nENSAT (n=15 samples)",
       y= element_blank(), 
       x= element_blank())+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 8, 
                                  face = "bold"))  

 # Mounting panel with TCGA-ACC and ENSAT boxplots
row1 <- plot_grid(g10, g10.ensat, ncol=2)
row2 <- plot_grid(g20, g20.ensat, ncol=2)
row3<- plot_grid(blankPlot, leg20, blankPlot, ncol=3, rel_widths = c(.5,1,.5))
g0 <- plot_grid(row1, row2, row3,
                ncol=1, rel_heights = c(1,1,.25))
g0
# Saving
ggsave(g0, 
      filename = "./Tumor_stage_TCGA_ENSAT.pdf", 
      device="pdf",
      width = 7, height=7)
```

## Code snippet for Figure 5B

```r

# Calculating median acitivity within regulon clusters
dat <- as.data.frame(t(regact_ENSAT))

dat.rc1 <- dat[c1[c1 %in% rownames(dat)],]
dat.rc1 <- apply(dat.rc1, 2, median)

dat.rc2 <- dat[c2[c2 %in% rownames(dat)],]
dat.rc2 <- apply(dat.rc2, 2, median)

dat.rc3 <- dat[c3[c3 %in% rownames(dat)],]
dat.rc3 <- apply(dat.rc3, 2, median)

dat.rc4 <- dat[c4[c4 %in% rownames(dat)],]
dat.rc4 <- apply(dat.rc4, 2, median)

identical(rownames(dat.rc1), rownames(dat.rc2)) #T
identical(rownames(dat.rc1), rownames(dat.rc3)) #T
identical(rownames(dat.rc3), rownames(dat.rc4)) #T

dat <- as.data.frame(
  cbind(RC1 = dat.rc1,
        RC2 = dat.rc2,
        RC3 = dat.rc3,
        RC4 = dat.rc4))

identical(rownames(dat), rownames(ENSAT_coldat))#T
# Mergind with clinical data
dat <- cbind(dat,
             C1A.C1B = ENSAT_coldat$Consensus.clustering.K2,
             Vital.Status = ENSAT_coldat$OS)
# Selecting samples with molecular classification (C1A, C1B)
dat <- dat[complete.cases(dat$C1A.C1B),]
dat$C1A.C1B <- as.factor(dat$C1A.C1B)
                                    
dat$Vital.Status <- as.factor((ifelse(dat$Vital.Status == 0,
                                      "Alive", "Dead")))
# Preparing data frame
dat_bp <- pivot_longer(dat,
                       cols = -contains(c("C1A.C1B",
                                          "Vital.Status")),
                       values_to = "value",
                       names_to = "Cluster")

# Tumor stages 1 and 2
## Statistical tests
stat.test <- rstatix::dunn_test(
    data=dat_bp[dat_bp$C1A.C1B=="C1A",], 
    formula= value ~ Cluster, 
    p.adjust.method = "BH")
stat.test$p.adj <- format(stat.test$p.adj, 
                          scientific = T, 
                          digits = 3)
print(stat.test)
stat.test$size <- ifelse(stat.test$p.adj.signif== "ns", 2, 4)
size <- c()
for(j in 1:nrow(stat.test)){
  size <- c(size, rep(stat.test$size[j], 3))
}
ymax <- round(max(dat_bp$value))
ymax <- 1.8
ypos <- c(ymax, ymax*1.2, ymax*1.6, ymax, ymax*1.4, ymax)
## Plotting
g1.C1A.ensat<-  
  ggplot(dat_bp[dat_bp$C1A.C1B=="C1A",],
         aes(x=Cluster, y=value, fill=Cluster))+
  labs(title= "C1A\nENSAT (n=26)",
       y= "Median Regulon\nCluster Activity", 
       x= element_blank()) +
  geom_violin(
    aes(fill=NULL, col=Cluster),
    alpha=0.4, 
    show.legend = c(T))+
  ggbeeswarm::geom_beeswarm(
    aes(color=Vital.Status),
    alpha =.5, 
    cex = 2, size=1, show.legend = T)+
  geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)),
                     breaks=c(-2,-1,0,1,2))+
  scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c",
                              "Alive"="gray50", "Dead"="red"))+
  scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c"),
                    labels = c("RC1 (n=62)", "RC2 (n=126)",
                               "RC3 (n=113)", "RC4 (n=68)"))+ 
  theme_plots +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10, 
                                  face = "bold"),
          axis.title.y = element_text(size=8))+
  stat_compare_means(method = "kruskal",
                     label.y = ymax*1.8,
                     size = 2.5)+
  stat_pvalue_manual(stat.test, 
                     inherit.aes=F,
                     label="p.adj.signif", 
                     bracket.shorten = .1,
                     tip.length = .02,
                     y.position = ypos, 
                     label.size = size)

# Tumor Stages 3 and 4
## Statistical tests
stat.test <- rstatix::dunn_test(
    data=dat_bp[dat_bp$C1A.C1B=="C1B",], 
    formula= value ~ Cluster, 
    p.adjust.method = "BH")
stat.test$p.adj <- format(stat.test$p.adj, 
                          scientific = T, 
                          digits = 3)
print(stat.test)
stat.test$size <- ifelse(stat.test$p.adj.signif== "ns", 2, 4)
size <- c()
for(j in 1:nrow(stat.test)){
  size <- c(size, rep(stat.test$size[j], 3))
}
## Plotting
g2.C1B.Ensat <-  
  ggplot(dat_bp[dat_bp$C1A.C1B=="C1B",],
         aes(x=Cluster, y=value, fill=Cluster))+
  labs(title= "C1B\nENSAT (n=18)",
       y= "Median Regulon\nCluster Activity", 
       x= element_blank())+
  geom_violin(
    aes(fill=NULL, col=Cluster),
    alpha=0.4, 
    show.legend = c(T))+
  ggbeeswarm::geom_beeswarm(
    aes(color=Vital.Status), 
    alpha =.5, cex = 2, size=1,
    show.legend = T)+
  geom_hline(yintercept=0, linetype=2, color="#800026", size=.25)+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.15)),
                     breaks=c(-2,-1,0,1,2))+
  scale_color_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                              "RC3"="#fb9a99","RC4"="#e31a1c",
                              "Alive"="gray50", "Dead"="red"))+
  scale_fill_manual(values=c("RC1"="#a6cee3","RC2"="#1f78b4",
                             "RC3"="#fb9a99","RC4"="#e31a1c"),
                    labels = c("RC1 (n=62)", "RC2 (n=126)",
                               "RC3 (n=113)", "RC4 (n=68)"))+ 
  theme_plots +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10, 
                                  face = "bold"),
        axis.title.y = element_text(size=8))+
  stat_compare_means(method = "kruskal",
                     label.y = ymax*1.8,
                     size = 2.5)+
  stat_pvalue_manual(stat.test, 
                     inherit.aes=F,
                     label="p.adj.signif", 
                     bracket.shorten = .1,
                     tip.length = .02,
                     y.position = ypos, 
                     label.size = size)

# Changing labels
g1.C1A.ensat <- g1.C1A.ensat +
  labs(title= "C1A ENSAT\n(n=26 samples)",
       y= element_blank(), 
       x= element_blank())+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 8, 
                                  face = "bold"))  

g2.C1B.Ensat <- g2.C1B.Ensat +
  labs(title= "C1B ENSAT\n(n=18 samples)",
       y= element_blank(), 
       x= element_blank())+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 8, 
                                  face = "bold"))  

row1 <- plot_grid(g10, g10.ensat, g20, g20.ensat, 
                  ncol=4, align = "hv",
                  labels = c("A","B","C","D"))
row2 <- plot_grid(g1.C1A.tcga, g1.C1A.ensat, g2.C1B.TCGA, g2.C1B.Ensat,
                  ncol=4, align = "hv",
                  labels = c("E","F","G","H"))
row3 <- plot_grid(blankPlot, leg.tumorstage, blankPlot, 
                  ncol=3, rel_widths = c(.25,1,.25))

gp <- plot_grid(row1, row2, row3, ncol = 1,
                rel_heights = c(1,1,.25))

 ggsave(gp, 
        filename = "./XX_C1_Stages.pdf",
        device="pdf", width = 8, height = 5, units="in")

 
 # Mounting panel with TCGA-ACC and ENSAT boxplots
row1 <- plot_grid(g1.C1A.tcga, g1.C1A.ensat, ncol=2)
row2 <- plot_grid(g2.C1B.TCGA, g2.C1B.Ensat, ncol=2)
row3<- plot_grid(blankPlot, leg20, blankPlot, ncol=3, rel_widths = c(.5,1,.5))
g0 <- plot_grid(row1, row2, row3,
                ncol=1, rel_heights = c(1,1,.25))
g0
# Saving
ggsave(g0, 
      filename = "./C1A_C1B_TCGA_ENSAT.pdf", 
      device="pdf",
      width = 7, height=7)
 
```

## Code snippet for Supplementary table 4 - Master Table
```r
## Preparing Master table
reg.cox <- reg.cox[colnames(regact_ENSAT),]
reg.cox <- cbind(reg.cox[,c("Regulons", "HR", "Adjusted.Pvalue")])

colnames(reg.cox) <- c("Regulons", "ENSAT.HR", 
                       "ENSAT.Cox.adjp")
# Arranging same number of rows
df <- regs.risk$SYMBOL[!regs.risk$SYMBOL %in% rownames(reg.cox)]
# [1] "ASCL5"    "C11orf95" "CENPBD1"  "DLX2"     "NME2"     "SRCAP"    "ZNF177"  
# [8] "ZNF595"  

df <- data.frame(Regulons = df, 
                 "ENSAT.HR"=rep(NA,8), 
                 "ENSAT.Cox.adjp" = rep(NA,8),
                 row.names = df)
# Merging with ENSAT Survival analysis
reg.cox <- rbind(reg.cox, df)
reg.cox <- reg.cox[regs.risk$SYMBOL,]
identical(rownames(reg.cox), regs.risk$SYMBOL) #T
regs.risk <- cbind(regs.risk, reg.cox)

# Merging with immune correlation
identical(rownames(regs.risk), rownames(regs.risk.cor)) #T
regs.risk <- cbind(regs.risk, regs.risk.cor[,3:22])

# Saving and exporting
writexl::write_xlsx(regs.risk, path = "./Risk_regulons_ENSAcorrigido.xlsx")
rm(colannot, df, rowannot , reg.cox, reg.KM)
```

# Code snippet for Supplementary Figure 6
## NR5A1 Network - Supplementary Figure 6A
```r
# Creating graph with NR5A1 targets
g <- tni.graph(rtni_tcgaACC, regulatoryElements = "NR5A1")
library(RedeR)
rdp <- RedPort()
calld(rdp, checkcalls=TRUE)
resetd(rdp)
addGraph(rdp, g, layout=NULL)

addLegend.color(rdp, g, type="edge")
addLegend.shape(rdp, g)
relax(rdp, p1=10, p2=10, p3=50, p4=10, p5=5, ps = TRUE)
rm(g, rdp)
```

## NR5A1 Karyogram - Supplementary Figure 6B
```r
# Retrieving targets for NR5A1 regulon
gen <- tnet[,"NR5A1"]
names(gen) <- rownames(tnet)
# Filtering and ordering for significant associations
gen <- gen[gen!=0]
gen <- gen[order(gen)]
# Listing gene names
gen <- names(gen)

# Retrieving row data for TCGA-ACC genes
genes <- data.frame(rowData(tcgaACC))
# Filtering for NR5A1 targets
genes <- genes[genes$gene_name %in% gen,]
# create a summarizedExperiment with the 248 NR5A1 targets
genes <- tcgaACC[rownames(tcgaACC) %in% genes$gene_id,]
# Extract GRanges from the SummarizedExperiment
genes.ranges <- rowRanges(genes)
# Retrieving metadata of Granges
metadata <- as.data.frame(values(genes.ranges))

# Creating data frame with mutual information of NR5A1 targets
gen <- tnet[,"NR5A1", drop=F]
rownames(gen) <- rownames(tnet)
gen <- as.data.frame(gen[gen!=0,,drop=F])
gen$Direction <- ifelse(gen$NR5A1>0, "Positive","Negative")

# Reordering data frame to merge with metadata adding direction of association
gen <- gen[metadata$gene_name,]
metadata <- cbind(metadata, rownames(gen), gen)
# Inserting new metadata information in the Granges object
values(genes.ranges) <- metadata

# Arranging seqlevels
genes.ranges = keepSeqlevels(genes.ranges, paste0("chr", c(1:22, "X", "Y")))
# Plotting karyogram
g<-ggbio::autoplot(genes.ranges, 
            layout = "karyogram", 
            aes(color=Direction, fill=Direction), 
            main="NR5A1 targets Karyogram",
            show.legend = c(T))+
  scale_color_manual(values = c("#00BFC4","#F8766D"),
                     labels = c("Negative (n=176)",
                               "Positive (n=72)"))+
  scale_fill_manual(values = c("#00BFC4","#F8766D"),
                    labels = c("Negative (n=176)",
                               "Positive (n=72)"))+
  theme(legend.position = "bottom", 
        legend.title.align = .5,
        legend.title = element_text(size=8), 
        legend.text = element_text(size=8))+
  guides(fill= guide_legend(title = "Association Direction",
                             title.position = "top"), 
         color="none")
g
# export 3.5x5.5"

```

# Code snippet for Figure 6
## CENPA Network - Figure 6A
```r
# Creating graph with CENPA targets
g <- tni.graph(rtni_tcgaACC, regulatoryElements = "CENPA")
library(RedeR)
rdp <- RedPort()
calld(rdp, checkcalls=TRUE)
resetd(rdp)
addGraph(rdp, g, layout=NULL)

addLegend.color(rdp, g, type="edge")
addLegend.shape(rdp, g)
relax(rdp, p1=10, p2=10, p3=50, p4=10, p5=5, ps = TRUE)
rm(g, rdp)
```

## CENPA Karyogram - Figure 6B
```r
# Retrieving targets for CENPA regulon
gen <- tnet[,"CENPA"]
names(gen) <- rownames(tnet)
# Filtering and ordering for significant associations
gen <- gen[gen!=0]
gen <- gen[order(gen)]
# Listing gene names
gen <- names(gen)

# Retrieving row data for TCGA-ACC genes
genes <- data.frame(rowData(tcgaACC))
# Filtering for NR5A1 targets
genes <- genes[genes$gene_name %in% gen,]
# create a summarizedExperiment with the 248 NR5A1 targets
genes <- tcgaACC[rownames(tcgaACC) %in% genes$gene_id,]
# Extract GRanges from the SummarizedExperiment
genes.ranges <- rowRanges(genes)
# Retrieving metadata of Granges
metadata <- as.data.frame(values(genes.ranges))

# Creating data frame with mutual information of CENPA targets
gen <- tnet[,"CENPA", drop=F]
rownames(gen) <- rownames(tnet)
gen <- as.data.frame(gen[gen!=0,,drop=F])
gen$Direction <- ifelse(gen$CENPA>0, "Positive","Negative")

# Reordering data frame to merge with metadata adding direction of association
gen <- gen[metadata$gene_name,]
metadata <- cbind(metadata, rownames(gen), gen)
# Inserting new metadata information in the Granges object
values(genes.ranges) <- metadata

# Arranging seqlevels
genes.ranges = keepSeqlevels(genes.ranges, paste0("chr", c(1:22, "X", "Y")))

g<-ggbio::autoplot(genes.ranges, 
            layout = "karyogram", 
            aes(color=Direction, fill=Direction), 
            main="CENPA targets Karyogram",
            show.legend = c(T))+
  scale_color_manual(values = c("#00BFC4","#F8766D"),
                     labels = c("Negative (n=44)",
                               "Positive (n=120)"))+
  scale_fill_manual(values = c("#00BFC4","#F8766D"),
                    labels = c("Negative (n=44)",
                               "Positive (n=120)"))+
  theme(legend.position = "bottom", 
        legend.title.align = .5,
        legend.title = element_text(size=8), 
        legend.text = element_text(size=8))+
  guides(fill= guide_legend(title = "Association Direction",
                             title.position = "top"), 
         color="none")
g
# export 3.5x5.5"

```

```r
# sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
# 
# locale:
#  [1] LC_CTYPE=pt_BR.UTF-8       LC_NUMERIC=C               LC_TIME=pt_BR.UTF-8       
#  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=pt_BR.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=pt_BR.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=pt_BR.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
#  [1] RedeR_2.0.0                          ensembldb_2.20.1                    
#  [3] AnnotationFilter_1.20.0              GenomicFeatures_1.48.0              
#  [5] hugene20sttranscriptcluster.db_8.8.0 org.Hs.eg.db_3.15.0                 
#  [7] AnnotationDbi_1.58.0                 GEOquery_2.64.0                     
#  [9] msigdbr_7.5.1                        reshape2_1.4.4                      
# [11] dendextend_1.16.0                    ConsensusClusterPlus_1.60.0         
# [13] RTNsurvival_1.20.0                   RTNduals_1.20.0                     
# [15] RTN_2.20.0                           ComplexHeatmap_2.12.0               
# [17] cowplot_1.1.1                        ggpubr_0.4.0                        
# [19] ggrepel_0.9.1                        ggfortify_0.4.14                    
# [21] forcats_0.5.1                        stringr_1.4.0                       
# [23] dplyr_1.0.9                          purrr_0.3.4                         
# [25] readr_2.1.2                          tidyr_1.2.0                         
# [27] tibble_3.1.8                         ggplot2_3.3.6                       
# [29] tidyverse_1.3.2                      DESeq2_1.36.0                       
# [31] AnnotationHub_3.4.0                  BiocFileCache_2.4.0                 
# [33] dbplyr_2.1.1                         SummarizedExperiment_1.26.1         
# [35] Biobase_2.56.0                       GenomicRanges_1.48.0                
# [37] GenomeInfoDb_1.32.1                  IRanges_2.30.0                      
# [39] S4Vectors_0.34.0                     BiocGenerics_0.42.0                 
# [41] MatrixGenerics_1.8.0                 matrixStats_0.62.0                  
# 
# loaded via a namespace (and not attached):
#   [1] rappdirs_0.3.3                rtracklayer_1.56.0           
#   [3] GGally_2.1.2                  minet_3.54.0                 
#   [5] R.methodsS3_1.8.2             ragg_1.2.2                   
#   [7] bit64_4.0.5                   knitr_1.39                   
#   [9] R.utils_2.12.0                DelayedArray_0.22.0          
#  [11] data.table_1.14.2             rpart_4.1.16                 
#  [13] KEGGREST_1.36.0               RCurl_1.98-1.8               
#  [15] doParallel_1.0.17             generics_0.1.3               
#  [17] snow_0.4-4                    RSQLite_2.2.15               
#  [19] proxy_0.4-27                  bit_4.0.4                    
#  [21] tzdb_0.3.0                    xml2_1.3.3                   
#  [23] lubridate_1.8.0               httpuv_1.6.5                 
#  [25] assertthat_0.2.1              viridis_0.6.2                
#  [27] gargle_1.2.0                  xfun_0.31                    
#  [29] hms_1.1.1                     babelgene_22.3               
#  [31] evaluate_0.15                 promises_1.2.0.1             
#  [33] fansi_1.0.3                   restfulr_0.0.15              
#  [35] progress_1.2.2                readxl_1.4.0                 
#  [37] igraph_1.3.4                  DBI_1.1.3                    
#  [39] geneplotter_1.74.0            htmlwidgets_1.5.4            
#  [41] reshape_0.8.9                 googledrive_2.0.0            
#  [43] ellipsis_0.3.2                backports_1.4.1              
#  [45] annotate_1.74.0               biomaRt_2.52.0               
#  [47] deldir_1.0-6                  pwr_1.3-0                    
#  [49] vctrs_0.4.1                   Cairo_1.6-0                  
#  [51] abind_1.4-5                   cachem_1.0.6                 
#  [53] withr_2.5.0                   BSgenome_1.64.0              
#  [55] checkmate_2.1.0               GenomicAlignments_1.32.0     
#  [57] prettyunits_1.1.1             cluster_2.1.4                
#  [59] segmented_1.6-0               lazyeval_0.2.2               
#  [61] crayon_1.5.1                  genefilter_1.78.0            
#  [63] labeling_0.4.2                pkgconfig_2.0.3              
#  [65] nlme_3.1-159                  vipor_0.4.5                  
#  [67] ProtGenerics_1.28.0           nnet_7.3-17                  
#  [69] rlang_1.0.4                   lifecycle_1.0.1              
#  [71] filelock_1.0.2                modelr_0.1.8                 
#  [73] dichromat_2.0-0.1             cellranger_1.1.0             
#  [75] graph_1.74.0                  Matrix_1.4-1                 
#  [77] carData_3.0-5                 reprex_2.0.1                 
#  [79] base64enc_0.1-3               beeswarm_0.4.0               
#  [81] GlobalOptions_0.1.2           googlesheets4_1.0.0          
#  [83] pheatmap_1.0.12               viridisLite_0.4.0            
#  [85] png_0.1-7                     rjson_0.2.21                 
#  [87] dunn.test_1.3.5               bitops_1.0-7                 
#  [89] R.oo_1.25.0                   KernSmooth_2.23-20           
#  [91] Biostrings_2.64.0             blob_1.2.3                   
#  [93] shape_1.4.6                   jpeg_0.1-9                   
#  [95] rstatix_0.7.0                 ggsignif_0.6.3               
#  [97] scales_1.2.0                  memoise_2.0.1                
#  [99] viper_1.30.0                  magrittr_2.0.3               
# [101] plyr_1.8.7                    zlibbioc_1.42.0              
# [103] compiler_4.2.1                BiocIO_1.6.0                 
# [105] RColorBrewer_1.1-3            clue_0.3-61                  
# [107] Rsamtools_2.12.0              cli_3.3.0                    
# [109] XVector_0.36.0                htmlTable_2.4.1              
# [111] Formula_1.2-4                 MASS_7.3-58                  
# [113] tidyselect_1.1.2              stringi_1.7.8                
# [115] textshaping_0.3.6             yaml_2.3.5                   
# [117] locfit_1.5-9.6                latticeExtra_0.6-30          
# [119] VariantAnnotation_1.42.0      tools_4.2.1                  
# [121] parallel_4.2.1                circlize_0.4.15              
# [123] rstudioapi_0.13               foreach_1.5.2                
# [125] foreign_0.8-82                gridExtra_2.3                
# [127] farver_2.1.1                  digest_0.6.29                
# [129] BiocManager_1.30.18           shiny_1.7.2                  
# [131] Rcpp_1.0.9                    egg_0.4.5                    
# [133] car_3.1-0                     broom_1.0.0                  
# [135] BiocVersion_3.15.2            later_1.3.0                  
# [137] writexl_1.4.0                 OrganismDbi_1.38.0           
# [139] httr_1.4.3                    ggbio_1.44.0                 
# [141] biovizBase_1.44.0             kernlab_0.9-31               
# [143] colorspace_2.0-3              rvest_1.0.2                  
# [145] XML_3.99-0.10                 fs_1.5.2                     
# [147] splines_4.2.1                 RBGL_1.72.0                  
# [149] systemfonts_1.0.4             xtable_1.8-4                 
# [151] jsonlite_1.8.0                R6_2.5.1                     
# [153] Hmisc_4.7-0                   pillar_1.8.0                 
# [155] htmltools_0.5.3               mime_0.12                    
# [157] glue_1.6.2                    fastmap_1.1.0                
# [159] BiocParallel_1.30.0           class_7.3-20                 
# [161] interactiveDisplayBase_1.34.0 codetools_0.2-18             
# [163] utf8_1.2.2                    lattice_0.20-45              
# [165] mixtools_1.2.0                curl_4.3.2                   
# [167] ggbeeswarm_0.6.0              magick_2.7.3                 
# [169] interp_1.1-3                  survival_3.4-0               
# [171] limma_3.52.0                  rmarkdown_2.14               
# [173] munsell_0.5.0                 e1071_1.7-11                 
# [175] GetoptLong_1.0.5              GenomeInfoDbData_1.2.8       
# [177] iterators_1.0.14              haven_2.5.0                  
# [179] gtable_0.3.0    
```



