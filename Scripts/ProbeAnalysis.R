library(limma)
library(sva)

targets <- read.table("Targets.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rawdata <- read.ilmn(files="Sample_probe_Profile.txt",ctrlfiles="Control_probe_Profile.txt")
#Order arrays in targets table as the expression data 
edata <- rawdata$E
batch <- targets$Batch
mod <- model.matrix(~as.factor(Covariate3), data=targets)
mod0 <- model.matrix(~1,data=targets)
num.sv(edata,mod,method="leek")
combat_edata = ComBat(dat=edata, batch=batch,mod=mod,numCovs=NULL, par.prior=TRUE,prior.plots=T)
rawdata$E <- combat_edata
bkcor_data<- neqc(rawdata)
pe <- propexpr(rawdata)
expressed <- rowSums(bkcor_data$other$Detection < 0.05) >= 3
data <- bkcor_data[expressed,]

#########################################################################################################
################################## Differential expression ##############################################
#########################################################################################################

control_critical <- factor(targets$Severity)
# > control_critical
  # [1] Moderate_0h  Control_0h   Control_0h   Critical_0h  Moderate_0h  Critical_0h  ...........
# Levels: Control_0h Control_24h Control_72h Critical_0h Critical_24h Critical_72h Moderate_0h Moderate_24h Moderate_72h
design <- model.matrix(~0+control_critical)
colnames(design) <- levels(control_critical)
dupcor<- duplicateCorrelation(data,design,block=targets$Subject,ndups=2)
fit <- lmFit(data,design,block=targets$Subject,correlation=dupcor$consensus.correlation)

## Critical Control at time points ##
contrasts <- makeContrasts(Critical_0h-Control_0h, Critical_24h-Control_24h,Critical_72h-Control_72h,levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2,trend=TRUE)
Critical_Control_TP_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")

## Criticals Time points  ##
contrasts <- makeContrasts(Critical_24h-Critical_0h,Critical_72h-Critical_24h,Critical_72h-Critical_0h,levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
Critical_TP_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")

## MODS vs NO_MODS ##
mods_tags <- factor(targets$MODS)
levels(mods_tags)
# [1] "MODS_0hr"     "MODS_24hr"    "MODS_72hr"    "No"           "NO_MODS_0hr"  "NO_MODS_24hr" "NO_MODS_72hr"
design <- model.matrix(~0+mods_tags)
colnames(design) <- levels(mods_tags)
dupcor<- duplicateCorrelation(data,design,block=targets$Subject,ndups=2)
fit <- lmFit(data,design,block=targets$Subject,correlation=dupcor$consensus.correlation)
contrasts <- makeContrasts(Mods0hr=(MODS_0hr - NO_MODS_0hr), MODS24hr=(MODS_24hr - NO_MODS_24hr), MODS72hr=(MODS_72hr - NO_MODS_72hr),levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
MODs_TP_tb <- topTable(fit2, number=dim(fit2)[1], adjust="BH")