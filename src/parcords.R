#parcords

library(ProjectTemplate)
load.project()

cord <- X2021.11.26.FINAL.cryoem.energy.bond.calc %>% 
  filter(interact != 'MRAS-PP1C') %>% 
  select(1:4) 

cord[cord == 0] <- NA

colnames(cord)

library(GGally)

GGally::ggparcoord(cord, columns = 2:4, scale = 'globalminmax')
