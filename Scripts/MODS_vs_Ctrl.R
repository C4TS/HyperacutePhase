# Load libraries, register cores, set seed
library(data.table)
library(limma)
library(qusage)
library(dplyr)
library(doMC)
registerDoMC(4)
set.seed(123)

### GENES ###

# Import data
clin <- fread('./Data/Clinical.csv') 
eset <- read.csv('./Data/NormData.csv', 
                 check.names = FALSE, row.names = 1) %>%
  as.matrix()
colnames(eset) <- gsub(' ', '_', colnames(eset))  
eset <- eset[, match(clin$Array, colnames(eset))]

# Design
des <- model.matrix(~ 0 + Time + Time:Group + Time:MODS, data = clin) 
des <- des[, !grepl('Critical', colnames(des))]
coefs <- c('noMODS_vs_Ctrl_0h', 'noMODS_vs_Ctrl_24h', 'noMODS_vs_Ctrl_72h', 
           'MODS_vs_Ctrl_0h', 'MODS_vs_Ctrl_24h', 'MODS_vs_Ctrl_72h')
colnames(des)[grepl('MODS', colnames(des))] <- coefs

# Fit model
icc <- duplicateCorrelation(eset, des, block = clin$Subject) 
fit <- lmFit(eset, des, correlation = icc$cor, block = clin$Subject)  
fit <- eBayes(fit)

# Results
res <- function(coef) {
  topTable(fit, coef = coef, number = Inf, sort.by = 'none') %>%
  mutate(Probe = rownames(fit)) %>%
  rename(AvgExpr = AveExpr,
         p.value = P.Value, 
             FDR = adj.P.Val) %>%
  arrange(p.value) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, Probe, AvgExpr, logFC, p.value, FDR) %>%
  fwrite(paste0('./Results/', coef, '.Genes.csv'))
}
foreach(coef = coefs) %dopar% res(coef)

### MODULES ###

# Import data
anno <- fread('./Data/ProbeMap.csv')
eset <- read.csv('./Data/norm_matrix_probes.csv', 
                 check.names = FALSE, row.names = 1) %>%
  as.matrix()
colnames(eset) <- gsub(' ', '_', colnames(eset))  
eset <- eset[, match(clin$Array, colnames(eset))]
df <- data.frame(ArrayAddress = as.numeric(rownames(eset))) %>%
  inner_join(anno, by = 'ArrayAddress')
rownames(eset) <- df$nuID 
mods <- fread('./Data/ChaussabelModules.csv')
mod_list <- lapply(unique(mods$Module), function(m) mods[Module == m, nuID])
names(mod_list) <- unique(mods$Module)

# Refit model
icc <- duplicateCorrelation(eset, des, block = clin$Subject) 
fit <- lmFit(eset, des, correlation = icc$cor, block = clin$Subject)  
fit <- eBayes(fit)

# Results 
res <- function(coef) {

  resid_mat <- residuals(fit, eset)
  mean <- fit$coefficients[, coef] 
  SD <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
  sd.alpha <- SD / (fit$sigma * fit$stdev.unscaled[, coef])
  sd.alpha[is.infinite(sd.alpha)] <- 1L
  dof <- fit$df.total
  
  # Run QuSAGE functions
  res <- newQSarray(mean = mean, SD = SD, sd.alpha = sd.alpha, dof = dof,
                    labels = rep('resid', ncol(eset)))  # Create QSarray obj
  res <- aggregateGeneSet(res, mod_list, 2L^14L)        # PDF per gene set
  res <- calcVIF(resid_mat, res, useCAMERA = FALSE)     # VIF on resid_mat
  
  # Export
  qsTable(res, number = Inf, sort.by = 'p') %>%
    mutate(Rank = row_number()) %>%
    rename(Module = pathway.name,
            logFC = log.fold.change,
          p.value = p.Value) %>%
    select(Rank, Module, logFC, p.value, FDR) %>%
    fwrite(paste0('./Results/', coef, '.Modules.csv'))
  
}
foreach(coef = coefs) %dopar% res(coef)
