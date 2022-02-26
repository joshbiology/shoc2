library(ProjectTemplate)
load.project()

x1 <- SHOC2_mutation_sequence %>% 
  rename(index = mutation,
         x = sequence)
x2 <- per_variant_LFC_SHOC2_DMS %>% 
  rename(index = variant.by.aa,
         y = LFC_ct_frctn_adj.z) 
  
inner_join(x1, x2, by="index") %>% 
  select(x, y) %>% 
  write_csv('./output/input_for_mavenn.csv')
