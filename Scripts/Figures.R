### BOX PLOTS ###

# Load libraries
library(data.table)
library(tidyverse)

# Import data
clin <- fread('./Data/Clinical.csv')[Group != 'Moderate']
dat <- read.csv('./Data/NormData.csv', row.names = 1, check.names = FALSE) 

### Crit vs. Ctrl ###
genes_1 <- c('CRISP3', 'EPHB4', 'FAM131A', 'GFI1', 'GJB6', 'HLX', 'LOC388514', 
             'LOC441582', 'LOC646575', 'LOC728835', 'OLR1', 'PARVG', 'PDXK', 
             'PTGS2', 'PTX3', 'SP100')

# Tidy
mat <- dat[genes_1, match(clin$Array, colnames(dat))]
df <- gather(mat, Sample, Expression) %>%
  mutate(Group = rep(clin$Group, each = nrow(mat)),
         Time = rep(clin$Time, each = nrow(mat)),
         Gene = rep(genes_1, times = ncol(mat)))

# Plot
ggplot(df, aes(Time, Expression, fill = Group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c('#00BA38', '#F8766D')) + 
  labs(title = 'Gene Expression by Group, Time', 
       x = 'Time', 
       y = 'Normalized Expression') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Gene, nrow = 4, ncol = 4)


### Ctrl, MODS, noMODS ###
genes_2 <- c('AXUD1', 'FRMD4B', 'GJB6', 'JMJD6', 'JOSD1', 'JUNB', 'LOC100127984', 
             'LOC728877', 'OLR1', 'OSM', 'PTX3', 'RNASE3', 'RPS6KA2', 'SH2D2A', 
             'SOCS3', 'TFF3')

# Tidy
clin_ctrl <- clin[Group == 'Control']
clin_mods <- clin[Group == 'Critical'
  ][, MODS := ifelse(MODS == 'Yes', 'MODS', 'NO_MODS')
  ][, MODS := relevel(as.factor(MODS), ref = 'NO_MODS')]
mat_ctrl <- dat[genes_2, match(clin_ctrl$Array, colnames(dat))]
df_ctrl <- gather(mat_ctrl, Sample, Expression) %>%
  mutate(Group = 'Control',
         Time = rep(clin_ctrl$Time, each = nrow(mat_ctrl)),
         Gene = rep(genes_2, times = ncol(mat_ctrl)))
mat_mods <- dat[genes_2, match(clin_mods$Array, colnames(dat))]
df_mods <- gather(mat_mods, Sample, Expression) %>%
  mutate(Group = rep(clin_mods$MODS, each = nrow(mat_mods)),
         Time = rep(clin_mods$Time, each = nrow(mat_mods)),
         Gene = rep(genes_2, times = ncol(mat_mods)))
df <- rbind(df_ctrl, df_mods)

# Plot
ggplot(df, aes(Time, Expression, fill = Group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c('#00BA38', '#619CFF', '#F8766D')) + 
  labs(title = 'Gene Expression by Group, Time', 
       x = 'Time', 
       y = 'Normalized Expression') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Gene, nrow = 4, ncol = 4)



### MODULE HEATMAPS ###

# Load libraries
library(data.table)
library(NMF)
library(RColorBrewer)
library(dplyr)

# Import, organise data
res1 <- fread('./Results/MODS_vs_noMODS_0h.Modules.csv') %>%
  rename(lfc0h = logFC) %>%
  filter(q.value <= 0.1)
res2 <- fread('./Results/MODS_vs_noMODS_24h.Modules.csv') %>%
  rename(lfc24h = logFC)
res <- fread('./Results/MODS_vs_noMODS_72h.Modules.csv') %>%
  rename(lfc72h = logFC) %>%
  inner_join(res1, by = 'Module') %>%
  inner_join(res2, by = 'Module') %>%
  select(Module, lfc0h, lfc24h, lfc72h)
mat <- res %>% 
  select(-Module) %>% 
  as.matrix()
dimnames(mat) <- list(res$Module, c('0h', '24h', '72h'))

# Plot
rg <- colorRampPalette(c('green', 'black', 'red'))(256)
aheatmap(mat, distfun = 'pearson', scale = 'row', col = rg, border_color = 'black', 
         cellwidth = 20, cexCol = 1, main = 'MODS vs NO_MODS')



