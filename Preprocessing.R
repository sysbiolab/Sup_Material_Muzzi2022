

## For GDC Data Portal ACC RNAseq import and data manipulation
library(TCGAbiolinks)
library(SummarizedExperiment)
## For filtering protein coding genes
library(AnnotationHub)
## Differential expression
library(DESeq2)
## Data manipulation and plots
library(tidyverse)
## Transcription network inference
library(RTN)
library(RTNduals)
library(RTNsurvival)
library(snow)

dir.create("~/ACC_Cancers2022")
setwd("~/ACC_Cancers2022")

# Download TCGA-ACC cohort data (clinical and RNA-seq) from GDC Portal
query <- GDCquery(project = "TCGA-ACC", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts",
                  sample.type = c("Primary Tumor"))

GDCdownload(query)
# Setting as SummarizedExperiment object 
tcgaACC <- GDCprepare(query, summarizedExperiment = TRUE) # 56457 genes

######################################################### Beginning - Prof Mauro Castro script
##--- Set a 'rowAnnotation' data frame for gene annotation 

ah <- AnnotationHub()
query(ah, "EnsDb")
edb <- query(ah, pattern = c("Homo sapiens", "EnsDb", 101))[[1]]
gns <- genes(edb, return.type= "GRanges")
gns <- gns[gns$gene_biotype=="protein_coding",]
rowAnnotation <- as.data.frame(gns)
######################################################## End - Prof Mauro Castro script

### Changing seqlevels style from AnnotationHub for subsetting
gns <- keepStandardChromosomes(gns, pruning.mode="coarse")
gns@seqnames #Levels(25): 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9 MT X Y

seqlevelsStyle(gns) #NCBI
seqlevelsStyle(tcgaACC) #"UCSC"

ns <- mapSeqlevels(seqlevels(gns), "UCSC")
gns <- renameSeqlevels(gns,ns)
seqlevelsStyle(gns) #[1] "NCBI"

# Filtering tcgaACC gene expression matrix for protein coding genes
tcgaACC <- subsetByOverlaps(tcgaACC, gns) #25476 genes
## Select patients with molecular classification on Steroid Phenotype
tcgaACC <- subset(tcgaACC, 
                  select = complete.cases(colData(tcgaACC)$paper_mRNA_K4)) 
# From 79 to 78 patients

# Cleaning environment
rm(ah,edb,gns,query,rowAnnotation, ns)

# Separating mRNA K4 classification in Steroid and Proliferation classification
colData(tcgaACC)$Steroid <- 
  as.factor(ifelse(colData(tcgaACC)$paper_mRNA_K4=="steroid-phenotype-low " |
                     colData(tcgaACC)$paper_mRNA_K4=="steroid-phenotype-low+proliferation",
                   "Steroid_Low",
                   "Steroid_High"))
colData(tcgaACC)$Proliferation_mRNA <- 
  as.factor(ifelse(colData(tcgaACC)$paper_mRNA_K4=="steroid-phenotype-high+proliferation" |  
                     colData(tcgaACC)$paper_mRNA_K4=="steroid-phenotype-low+proliferation",
                   "Proliferation_High",
                   "Proliferation_Low"))

coldat <- data.frame(colData(tcgaACC))

# Saving R.Data for further analysis
save(tcgaACC, file = "./tcgaACC_Steroid.RData")

######## Differential expression analysis with DESeq2 using the Steroid Phenotype as design
###### Normalizing counts - vst = variance stabilizing transformation from DESeq2
dds <- DESeqDataSet(tcgaACC, design = ~ Steroid)
dds <- DESeq(dds)
vsd <- vst(dds, blind = F) # for information on "blind" argument, see DESeq2 vignette

# Extracting normalized gene expression matrix
gen_exp <- as.data.frame(assay(vsd))

# Changing column name to Symbol for RTNI construction using Symbol gene names
colnames(rowData(tcgaACC)) #[1] "ensembl_gene_id" "external_gene_name" "original_ensembl_gene_id"
colnames(rowData(tcgaACC)) <- c("ENSEMBL", "SYMBOL", "OG_ENSEMBL")
geneannot <- rowData(tcgaACC)
# Retrieving TFs from Lambert et al 2018 present in Gene Expression Matrix
data("tfsData")
regulatoryElements <- intersect(tfsData$Lambert2018$SYMBOL, geneannot$SYMBOL)
# Constructing RNTI object with TFs from Lambert et al 2018 as Regulatory Elements
identical(colnames(tcgaACC), colnames(gen_exp)) #TRUE
identical(rownames(tcgaACC), rownames(gen_exp)) #TRUE
rtni_tcgaACC <- tni.constructor(expData = as.matrix(gen_exp), 
                                regulatoryElements = regulatoryElements,
                                rowAnnotation = as.data.frame(rowData(tcgaACC)),
                                colAnnotation = as.data.frame(colData(tcgaACC)))
# --Removing inconsistent data: standard deviation is zero for 1243 gene(s)! 
rm(geneannot, tcgaACC, tfsData, regulatoryElements)
# Allocating memory for permutation analysis (memory consuming)
options(cluster=snow::makeCluster(spec=4, "SOCK"))

# Defining p-value for RTN permutation analysis
# Robertson et al. (2017) used p-value = 1e-5 for n= 405
tni.alpha.adjust(nB=78, nA=405, alphaA = 1e-5) # 0.15
# From RTN vignette:
tni.alpha.adjust(nB = 78, nA = 300, alphaA = 1e-5, betaA = 0.2) # 0.067
# Permutation analysis for Regulons inference, p-value<0.05 for small cohorts
rtni_tcgaACC <- tni.permutation(rtni_tcgaACC, 
                                pValueCutoff = 0.05,
                                boxcox = T)
# Bootstrap analysis
rtni_tcgaACC <- tni.bootstrap(rtni_tcgaACC)
# DPI filter
rtni_tcgaACC <- tni.dpi.filter(rtni_tcgaACC)
stopCluster(getOption("cluster"))
# Saving RTNI object
save(rtni_tcgaACC, file = "rtni.RData")

#############################################################################
# ENSAT cohort
#############################################################################
library(GEOquery)
library(Biobase)

#Download GEO processed data with GEOquery
#Sys.setenv("VROOM_CONNECTION_SIZE"=131072*4)
gds <- getGEO("GSE49278")
e <- gds[[1]]
dim(e) 
#Features  Samples 
#   53617       44 
ENSAT_pdat <- pData(e)

# Retrieve Genes annotation data
rowdat_ENSAT <- fData(e)
gexp_ENSAT <- exprs(e)
rm(e, gds)

# Getting gene names from Human Gene 2.0 ST Array
library(hugene20sttranscriptcluster.db)

gene.names <- select(hugene20sttranscriptcluster.db,
                     rownames(gexp_ENSAT), 
                     c("ENTREZID","SYMBOL","GENENAME"))

gene.names <- gene.names[complete.cases(gene.names$SYMBOL),]

# Getting updated gene annotation for filtering for protein coding genes
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, "EnsDb")
edb <- query(ah, pattern = c("Homo sapiens", "EnsDb", 103))[[1]]
gns <- genes(edb, return.type= "data.frame")
gns <- gns[gns$gene_biotype=="protein_coding",]
gns <- gns[,c("gene_id","gene_name","symbol","entrezid","gene_biotype")]

# Selecting ENTREZ names, between duplicated entries, we will select the ones with higher coeffient of variation
rowdat.up_ENSAT <- gene.names[gene.names$SYMBOL %in% gns$symbol,]

gexp2 <- as.data.frame(gexp_ENSAT)
gexp2 <- gexp2[rowdat.up_ENSAT$PROBEID,]
sd <- apply(gexp2, 1, sd)
mean <- apply(gexp2, 1, mean)
cv <- as.vector(sd/mean)
gexp2$cv <- cv

rowdat.up_ENSAT <- rowdat.up_ENSAT[rowdat.up_ENSAT$PROBEID %in% rownames(gexp2),]
gexp2 <- gexp2[rownames(gexp2) %in% rowdat.up_ENSAT$PROBEID,] #20793 genes

# Merging gexp data with Rowdata
rowdat.up_ENSAT <- rowdat.up_ENSAT[!duplicated(rowdat.up_ENSAT$PROBEID),]
identical(rowdat.up_ENSAT$PROBEID, rownames(gexp2)) # T
rowdat.up_ENSAT$CV <- gexp2$cv
# Finding duplicated entries
dup <- rowdat.up_ENSAT[duplicated(rowdat.up_ENSAT$ENTREZID) | duplicated(rowdat.up_ENSAT$ENTREZID, fromLast = T),]

# Ordering by the coefficient of variation to select duplicated genes with higher cv
dup <- dup[order(dup$ENTREZID, dup$CV, decreasing = T),]
# Discover the duplicated entries to remove from gexp
# This will retain only the maximum value of each duplicated gene
dup <- dup[duplicated(dup$ENTREZID),]

rowdat.up_ENSAT <- rowdat.up_ENSAT[!rowdat.up_ENSAT$PROBEID %in% dup$PROBEID,]

# Find if there is still any duplicated entries in PROBE ID, Entrez or Symbol
any(duplicated(rowdat.up_ENSAT$PROBEID)) #F
any(duplicated(rowdat.up_ENSAT$ENTREZID)) #F
any(duplicated(rowdat.up_ENSAT$SYMBOL)) #F

# Filtering gene expression matrix for genes selected
gexp_ENSAT <- gexp_ENSAT[rowdat.up_ENSAT$PROBEID,] # 18624 genes
identical(rownames(gexp_ENSAT), rowdat.up_ENSAT$PROBEID) #T
# Changing rownames to symbols
rownames(gexp_ENSAT) <- rowdat.up_ENSAT$SYMBOL

# In Assié et al., 2014 this data is already normalized by RMA
rm(gexp2, gns, ah, edb, gene.names, rowdat_ENSAT, cv, mean, sd)

# Clinical data from Sup tab 1b (Assié et al. 2014)
ENSAT_clinicaldat <-data.frame(
  stringsAsFactors = FALSE,
  check.names = FALSE,
  Sample = c("ACC1","ACC2","ACC3","ACC4","ACC5","ACC6",
             "ACC7","ACC8","ACC9","ACC10","ACC11",
             "ACC12","ACC13","ACC14","ACC15","ACC16",
             "ACC17","ACC18","ACC19","ACC20","ACC21",
             "ACC22","ACC23","ACC24","ACC25","ACC26",
             "ACC27","ACC28","ACC29","ACC30","ACC31",
             "ACC32","ACC33","ACC34","ACC35","ACC36",
             "ACC37","ACC38","ACC39","ACC40",
             "ACC41","ACC42","ACC43","ACC44","ACC45",
             "ACC46","ACC47","ACC48","ACC49","ACC50",
             "ACC51","ACC52","ACC55"),
  Sex = c("F",
          "F","M","F","M","F","F","F","M","F",
          "M","F","F","F","F","F","F","F","F",
          "M","F","F","F","F","F","F","F",
          "F","F","F","M","F","F","F","F","M",
          "F","F","F","F","F","M","M","F","F",
          "F","F","M","F","F","M","F","M"),
  Age.at.diagnosis = c(70.3,
                       25.6,40,53.3,72.9,18.3,77.5,50.7,
                       63.9,27,29.6,79.3,46.2,43,53.9,45,41,
                       37.2,81.6,67.5,42.3,39.7,25.2,41.7,
                       37.9,23.9,59.5,75.5,37.6,34.1,26.1,26.5,
                       48.4,58.8,49.6,54.3,79.6,29,44.5,
                       28.5,68.9,28.9,52.4,30,46.3,59.4,18.6,
                       39.6,40.2,53.8,30.6,44.6,52.9),
  Tumor.size = c(55,
                 50,70,80,120,45,270,80,80,
                 130,160,140,110,160,100,100,50,
                 240,90,50,55,85,110,150,90,
                 95,110,170,90,90,140,80,150,
                 110,70,160,90,70,220,65,120,
                 150,200,75,100,160,40,110,60,
                 150,250,NA,NA),
  Tumor.side = c("left","left","right","right","right",
                 "right","right","right","right","left",
                 "right","right","left","left","right","right",
                 "left","left","right","right","left",
                 "left","left","left","left","left",
                 "left","right","left","right","right",
                 "left","left","right","right","left","right",
                 "right","right","left","right","left",
                 "right","left","right","left","right",
                 "right","right","left","right","Abdo",
                 "left"),
  Hormonal.Secretion = c("no",
                         "yes","no","yes","yes","yes","yes",
                         "yes","yes","yes","yes","yes","yes",
                         "yes","no","yes","no","yes","yes","yes",
                         "yes","no","yes","yes","no","yes","yes",
                         "yes","yes","yes","yes","yes","yes",
                         "yes","yes","yes","yes","yes","yes",
                         "yes","no","yes","yes","yes","yes","yes",
                         "yes","yes","yes","yes","no","no",
                         "yes"),
  Weiss.score = c(4,
                  3,7,4,8,2,8,5,8,7,8,7,
                  8,6,6,7,3,9,2,6,4,6,5,
                  8,3,9,5,5,6,7,7,9,6,3,
                  2,7,6,2,4,3,6,4,8,5,3,
                  6,6,3,2,6,5,8,5),
  ENSAT.staging = c(2,
                    1,4,2,4,1,2,2,4,2,4,2,
                    4,3,2,NA,1,4,2,1,2,2,4,
                    2,2,4,2,2,2,3,4,4,4,2,
                    2,3,4,2,2,2,3,4,4,2,2,
                    2,2,2,2,4,2,NA,3),
  Recurrence = c("no",
                 "no","yes","no","yes","no","yes",
                 "yes","yes","no","yes","no","yes","yes",
                 "no","yes","no","yes","no","yes","no",
                 "no","yes","yes","yes","yes","no",
                 "yes","yes","yes","yes","yes","yes","yes",
                 "no","yes","yes","no","no","no","yes",
                 "yes","yes","no","no","yes","yes",
                 "yes","no","yes","no","yes","yes"),
  `Time.to.recurrence.(months)` = c(NA,
                                    NA,0,NA,0,NA,0,7.3,0,NA,0,NA,0,
                                    6.2,NA,21.8,NA,0,NA,19.3,NA,NA,0,18,
                                    36,0,NA,10,56.8,7.1,0,0,11.8,31,
                                    NA,29.3,0,NA,NA,NA,7,0,0,NA,NA,11,
                                    18.4,17.3,NA,0,NA,240,4.5),
  Specific.Death = c("no",
                     "no","yes","no","yes","no","yes",
                     "yes","yes","no","yes","no","yes","yes",
                     "no","yes","no","yes","no","yes","no",
                     "no","yes","no","no","yes","no","no",
                     "no","yes","yes","yes","yes","yes",
                     "no","yes","yes","no","no","no","yes",
                     "yes","yes","no","no","no","yes","yes",
                     "no","no","no","no","yes"),
  `Follow.up.(months)` = c(151.8,
                           131.1,23,147.8,0.5,142.8,5.2,36,0.6,
                           73.6,74.4,12,12.7,9.5,115.1,85.1,
                           59.6,11.3,40.2,40.2,41.6,57.2,48.3,35.9,
                           99.4,9.4,40.9,28,81.8,11.7,30.2,2,
                           19.8,34.3,129.8,57.9,7.6,118,119.5,
                           144.3,19,21,24.3,111.5,97.6,12.9,55.4,
                           26,158.8,154.2,85.1,351.7,11.4))

#C1A C1B classification Sup table 11
ENSAT_c1a.c1b <- data.frame(
  stringsAsFactors = FALSE,
  check.names = FALSE,
  Sample = c("ACC1","ACC10","ACC11","ACC12", "ACC13","ACC14",
             "ACC15","ACC16","ACC17","ACC18","ACC19",
             "ACC2","ACC20","ACC21","ACC22","ACC23","ACC24",
             "ACC25","ACC26","ACC27","ACC28","ACC29",
             "ACC31","ACC32","ACC33","ACC35",
             "ACC36","ACC37","ACC38","ACC39","ACC4",
             "ACC40","ACC42","ACC43","ACC44",
             "ACC45","ACC46","ACC47","ACC48","ACC49",
             "ACC5","ACC50","ACC52","ACC55",
             "ACC6","ACC8","ACC9"),
  `Consensus.clustering.K2` = c(
    "C1B","C1B","C1B","C1A","C1A", "C1A","C1B","C1A",
    "C1B","C1A","C1A","C1B","C1A","C1B","C1B","C1A",
    "C1A","C1B","C1A","C1A","C1B","C1B","C1A","C1A","C1A",
    "C1B","C1A","C1A","C1A","C1A","C1B","C1B","C1A","C1A",
    "C1A","C1B","C1A","C1A","C1A","C1B","C1A","C1B",
    "C1B","C1A","C1A","C1A","C1B"),
  CIMP = c(
    "CIMP.low", "CIMP.low", "CIMP.low", "CIMP.low", 
    "CIMP.low", "CIMP.low", "CIMP.low", "CIMP.low", 
    "CIMP.low", "CIMP.high", "non.CIMP", "non.CIMP", 
    NA, NA, "CIMP.high", NA, "CIMP.high", "CIMP.high",
    "CIMP.high", "CIMP.low", "CIMP.high", "CIMP.high",
    "non.CIMP", "non.CIMP", "non.CIMP", "non.CIMP",
    "non.CIMP", "non.CIMP", "CIMP.low", "non.CIMP", 
    "non.CIMP", "non.CIMP", "non.CIMP", "CIMP.low",
    "CIMP.low", "non.CIMP", "non.CIMP", "non.CIMP", 
    "non.CIMP", "non.CIMP", "non.CIMP", "non.CIMP", 
    "non.CIMP", "CIMP.low", "CIMP.low", "non.CIMP", "non.CIMP"))

# Merging clinical data

ENSAT_pdat[!ENSAT_pdat$title %in% ENSAT_clinicaldat$Sample, "title"] # 0

ENSAT_pdat[!ENSAT_pdat$title %in% ENSAT_c1a.c1b$Sample, "title"] #0

ENSAT_clinicaldat[!ENSAT_clinicaldat$Sample %in%  ENSAT_pdat$title, "Sample"] 
# "ACC3"  "ACC7"  "ACC30" "ACC31" "ACC34" "ACC41" "ACC51" "ACC52" "ACC55"

ENSAT_c1a.c1b[!ENSAT_c1a.c1b$Sample %in%  ENSAT_pdat$title, "Sample"] 
#"ACC31" "ACC52" "ACC55"

identical(colnames(gexp_ENSAT), ENSAT_pdat$geo_accession) #T

ENSAT_clinicaldat <- 
  ENSAT_clinicaldat[ENSAT_clinicaldat$Sample %in%  ENSAT_pdat$title,]
rownames(ENSAT_clinicaldat) <- ENSAT_clinicaldat$Sample

ENSAT_c1a.c1b <- ENSAT_c1a.c1b[ENSAT_c1a.c1b$Sample %in%  ENSAT_pdat$title,]
rownames(ENSAT_c1a.c1b) <- ENSAT_c1a.c1b$Sample

identical(rownames(ENSAT_c1a.c1b), rownames(ENSAT_clinicaldat)) #F
ENSAT_c1a.c1b <- ENSAT_c1a.c1b[ENSAT_pdat$title,]
ENSAT_clinicaldat <- ENSAT_clinicaldat[ENSAT_pdat$title,]
identical(rownames(ENSAT_c1a.c1b), rownames(ENSAT_clinicaldat)) #T
identical(ENSAT_pdat$title, rownames(ENSAT_clinicaldat)) #T

# Creating column annotation for ENSAT cohort
ENSAT_coldat <- cbind(ENSAT_pdat[,c(1,2,34,35,36)], 
                      ENSAT_clinicaldat, 
                      ENSAT_c1a.c1b)

ENSAT_coldat$OS <- ifelse(ENSAT_coldat$Specific.Death == "yes", 1, 0)
ENSAT_coldat$OS.time <- as.numeric(ENSAT_coldat$`Follow.up.(months)`)

# Saving and exporting objects
ENSAT <- list(gexp_ENSAT, ENSAT_coldat)
save(ENSAT, file = "./ENSAT.RData")



