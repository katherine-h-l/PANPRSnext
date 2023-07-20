# ----- Packages -----#
### Use GenomicSEM to obtain munged summary statistics
# devtools::install_github("GenomicSEM/GenomicSEM")
library(GenomicSEM)

# install.packages("~/Downloads/PANPRSnext_0.1.0.tar.gz", repos = NULL, type="source")
library(PANPRSnext)
library(PANPRS)

# ---- 0. Download and preprocess data ---- #
# I followed the tutorials:
# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
# https://gist.github.com/bulik/c194b929a340ed3ff610

{
# Download data (in Terminal)
# cd ~/Your/Download/Folder
# wget www.med.unc.edu/pgc/files/resultfiles/pgc.cross.bip.zip
# wget www.med.unc.edu/pgc/files/resultfiles/pgc.cross.scz.zip
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
## Extract
# tar -jxvf eur_w_ld_chr.tar.bz2
# bunzip2 w_hm3.snplist.bz2
# unzip -o pgc.cross.bip.zip
# unzip -o pgc.cross.scz.zip
}


# ---- 1. Pre-process info ---- #

setwd("~/Your/Download/Folder")

### Load summary statistics
BIP <- read.table("pgc.cross.BIP11.2013-05.txt", header=T)
BIP$CEUaf <- as.numeric(BIP$CEUaf)
summary(BIP)
SCZ <- read.table("pgc.cross.SCZ17.2013-05.txt", header=T)
summary(SCZ)

### Load LD scores and num of SNPs
ldscore0 <- do.call(rbind, lapply(1:22, function(chr){
  ldsc_chr <- read.table(gzfile(paste0("eur_w_ld_chr/",chr,".l2.ldscore.gz")), header=T)
  return(ldsc_chr)}))
dim(ldscore0[ldscore0$MAF>=0.05,])

M_5_50 <- do.call(rbind, lapply(1:22, function(chr){
  ldsc_chr <- read.table(paste0("eur_w_ld_chr/",chr,".l2.M_5_50"), header=F)
  return(ldsc_chr)}))
sum(M_5_50$V1)

### Load hm3 SNPs
hm3 <- read.table("w_hm3.snplist", header=T)

### Following
# https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects

### Munge summary stats
#create vector of the summary statistics files
files <- c("pgc.cross.BIP11.2013-05.txt", "pgc.cross.SCZ17.2013-05.txt")

#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3
hm3 <- "eur_w_ld_chr/w_hm3.snplist"

#name the traits
trait.names <- c("BIP","SCZ")

#list the sample sizes. All but PTSD have SNP-specific sum of effective sample sizes so only its
#sample size is listed here
N = c(11810,17115)

#definte the imputation quality filter
info.filter = 0.9

#define the MAF filter
maf.filter = 0.01

#run munge function to generate munged summary statistics files: "BIP.sumstats.gz" and "SCZ.sumstats.gz"
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)


# ---- 2. Load needed files ---- #

BIPsumm <- read.table(gzfile("~/Desktop/PANPRS-next-Katya/NewDataSet/BIP.sumstats.gz"), header=T)
SCZsumm <- read.table(gzfile("~/Desktop/PANPRS-next-Katya/NewDataSet/SCZ.sumstats.gz"), header=T)
### Ensure A1 is the same between pairs and flip the sign of one if not
summaryZ1 = merge(BIPsumm, SCZsumm, by="SNP", all.x=T, all.y=T)
# all(summaryZ1$A1.x==summaryZ1$A1.y, na.rm=T)
rownames(summaryZ1) = summaryZ1$SNP
summaryZ1 <- summaryZ1[,c("Z.x","Z.y")]
names(summaryZ1) <- c("ZBIP","ZSCZ")

### One could subset data using:
set.seed(1)
summaryZ1sub <- summaryZ1[sample(1:nrow(summaryZ1),8e4),]

### Load precomputed data.frame "plinkLD" (on OneDrive)
load(file="~/Desktop/PANPRS-next-Katya/NewDataSet/1KG_QC_EUR_pruned_plinkLD.RData")

# ---- 3. Run PANPRS ---- #
output0a <- PANPRS::gsPEN(na.omit(summaryZ1sub), Nvec = c(11810, 17115), plinkLD = plinkLD, warmStart = 1)

## Mock funcIndex (if needed)

Nrws <- nrow(na.omit(summaryZ1sub))
funcIndex_mock <- data.frame(matrix(0, nrow=Nrws, ncol=2))
names(funcIndex_mock) <- c("V1","V2")
set.seed(1)
funcIndex_mock$V1[sample(Nrws,ceiling(0.25*Nrws))] <- 1
funcIndex_mock$V2[sample(Nrws,ceiling(0.55*Nrws))] <- 1


output0a <- PANPRS::gsfPEN(na.omit(summaryZ1sub), Nvec = c(11810, 17115), plinkLD = plinkLD, funcIndex = funcIndex_mock, numfunc=2, warmStart = 1)

output0b <- PANPRSnext::gsfPEN_R(summary_z=na.omit(summaryZ1sub), n_vec = c(11810, 17115), plinkLD = plinkLD, func_index = funcIndex_mock)

