# Load libraries
library(data.table)
library(limma)
library(qvalue)
library(qusage)
library(dplyr)
set.seed(123)

# Import, format data
clin <- fread('./Data/Clinical.csv')
eset <- read.ilmn(files = './Data/Sample_probe_Profile_5_12_13.txt',
                  ctrlfiles = './Data/Control_probe_Profile_5_12_13.txt')
mods <- fread('./Data/ChaussabelModules.csv')
mod_list <- lapply(unique(mods$Module), function(m) mods[Module == m, nuID])
names(mod_list) <- unique(mods$Module)
anno <- fread('./Data/ProbeMap.csv')

# Background correct, normalize, filter, reset dimnames
eset <- neqc(eset)
keep <- rowSums(eset$other$Detection < 0.05) >= 3
colnames(eset) <- gsub(' ', '_', colnames(eset))
eset <- eset[keep, match(clin$Array, colnames(eset))]
rownames(eset) <- anno$nuID

# Big ol' loop
loop <- function(variable) {
  
  # Fit limma models
  if (variable == 'Time') {
    des <- model.matrix(~ 0 + Subject + Batch + Time:Group, data = clin)
    des <- des[, !grepl('0h', colnames(des))]
    coefs <- c('0h_vs_24h_Crit', '0h_vs_72h_Crit')
    colnames(des)[grepl('Critical', colnames(des))] <- coefs
    fit <- eBayes(lmFit(eset, des), trend = TRUE, robust = TRUE)
  } else {
    if (variable == 'Crit') {
      des <- model.matrix(~ 0 + Time + Batch + Sex + Age + Time:Group, data = clin)
      coefs <- c('Crit_vs_Ctrl_0h', 'Crit_vs_Ctrl_24h', 'Crit_vs_Ctrl_72h')
      colnames(des)[grepl('Critical', colnames(des))] <- coefs
      icc <- duplicateCorrelation(eset, des, block = clin$Subject)
    } else if (variable == 'MODS') {
      des <- model.matrix(~ 0 + Time + Batch + Sex + Age + Time:MODS, data = clin)
      coefs <- c('MODS_vs_noMODS_0h', 'MODS_vs_noMODS_24h', 'MODS_vs_noMODS_72h')
      colnames(des)[grepl('Yes', colnames(des))] <- coefs
      icc <- duplicateCorrelation(eset, des, block = clin$Subject)
    }
    fit <- lmFit(eset, des, correlation = icc$cor, block = clin$Subject)
    fit <- eBayes(fit, robust = TRUE, trend = TRUE)
  }
  
  # Run QuSAGE pipeline
  for (coef in coefs) {
    
    # Prep data
    j <- which(colnames(des) == coef)
    se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, j]
    sd_a <- se / (fit$sigma * fit$stdev.unscaled[, j])
    sd_a[is.infinite(sd_a)] <- 1
    resid_mat <- residuals.MArrayLM(fit, eset)
    overlap <- sapply(mod_list, function(x) sum(x %in% rownames(fit)))
    mod_list <- mod_list[overlap > 0]
    
    # QuSAGE functions
    res <- newQSarray(mean = fit$coefficients[, j],         # Create QSarray obj
                      SD = se,
                      sd.alpha = sd_a,
                      dof = fit$df.total,
                      labels = rep('Resid', ncol(eset)))
    res <- aggregateGeneSet(res, mod_list, n.points = 2^14)  # PDF per gene set
    res <- calcVIF(resid_mat, res, useCAMERA = FALSE)        # VIF on resid_mat
    
    # Export
    qsTable(res, number = Inf, sort.by = 'p') %>%
      rename(p.value = p.Value,
             Module = pathway.name,
             logFC = log.fold.change) %>%
      mutate(q.value = qvalue(p.value)$qvalues) %>%
      select(Module:p.value, q.value) %>%
    fwrite(paste0('./Results/Modules/', coef, '.Modules.csv'))
    
  }
  
}

# Compute in parallel
library(doParallel)
registerDoParallel(3)
foreach(v = c('Time', 'Crit', 'MODS')) %dopar% loop(v)


