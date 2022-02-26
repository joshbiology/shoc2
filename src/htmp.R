library(ProjectTemplate)
load.project()


# Matrix ------------------------------------------------------------------


SHOC2_DMS <- read.csv(here::here("data","File1_reformatted_SHOC2_callapsedToAA.csv"))
xpheatmap <- read.csv(here::here('data','xping.heatmap.csv'))
REAL_LRR_POSITIONS <- read.csv(here::here('data', 'Correct_LRRnum_CIRCpos.csv'))
SHOC2_DMS<- xpheatmap[,c(5,7)] %>% left_join(SHOC2_DMS, by = "variant.by.aa" )

ggplot(SHOC2_DMS, aes(x=LFC_ct_frctn_adj.z, color = Type)) + geom_histogram()

SHOC2_DMS <- SHOC2_DMS %>% 
  dplyr::select(c(POS, Vt_aa, LFC_ct_frctn_adj.z)) %>% 
  pivot_wider(names_from = POS, values_from = LFC_ct_frctn_adj.z) %>% 
  column_to_rownames(var = "Vt_aa") %>% 
  as.matrix() 


silent_effect <- SHOC2_DMS['B',] %>% mean(na.rm=T)


tmp <- SHOC2_DMS - silent_effect

nonsense_effect <- tmp['X',] %>% mean(na.rm=T)

rescaled <- tmp/abs(nonsense_effect)

rescaled['X',] %>% mean(na.rm=T)
rescaled['B',] %>% mean(na.rm=T)


mat <- rescaled[setdiff(rownames(rescaled), c('B', 'X')),]


stability <- tibble(ddg = avg_ddg %>% column_to_rownames("res") %>% as.matrix() %>% colMeans(),
              POS = 1:length(ddg)) %>% 
  filter(ddg > 2) %>% 
  mutate(Destabilizing = "destabilized") %>% 
  select(POS, Destabilizing)


evo <- Aminode_accessed_2021.11.23[,c(1,3)]

# Annotations -------------------------------------------------------------


tmp <- list("Small Hydrophobic" = c("A", "I", "L", "M", "V"), "Special Cases" = c("C","G","P"), "Polar Uncharged" = c("S","T","N","Q"), "Positive" = c("K","R","H"), "Negative" = c("D","E"), "Large Hydrophobic" = c("F", "W", "Y"), "Nonsense" = c("X"), "Silent" = c("B")) %>% 
  enframe() %>% 
  unnest()

contact_types <- split(contact_v_nocontact, contact_v_nocontact$Core_resi)

tmp2 <- tibble(POS = 1:582) %>% #left_join(stability) %>% 
  left_join(REAL_LRR_POSITIONS %>% 
              select(POS, LRR_xtal)) %>% 
  left_join(contact_v_nocontact %>% 
              filter(Core_resi %in% c('PP1C_interact', 'MRASandPP1C_interact', 'MRAS_interact')) %>% 
              rename(PPI = Core_resi)) %>% 
  #left_join(contact_types$Core_Residue %>% rename(Core = Core_resi)) %>% 
  #left_join(contact_types$Non_Contact_Surface_Residue %>% rename(Surface = Core_resi)) %>% 
  left_join(evo) %>% 
  column_to_rownames(var = "POS") 


gaps <- REAL_LRR_POSITIONS %>% 
  group_by(LRR_xtal) %>% 
  summarize(POS = min(POS)) %>% 
  pull(POS)

# Plotting ----------------------------------------------------------------

ann_colors = list(
  Substitution.Score = c("white", "black"),
  LRR_xtal = rep("#00B18B",20) %>% set_names(REAL_LRR_POSITIONS$LRR_xtal %>% unique()),
  PPI = c("MRAS_interact" = '#ED2386', 
          "PP1C_interact" = '#FAA41A', 
          "MRASandPP1C_interact" = '#675EA9')
)

library(pheatmap)
pheatmap::pheatmap(mat[tmp$value[1:20],],
                   border_color = "black",
                   cluster_rows=FALSE,
                   breaks = seq(-3, 1, length.out=1000),
                   color=colorRampPalette(c("navy",
                                            "navy",
                                            "navy",
                                            "white", 
                                            "red"))(1000), 
                   cluster_cols=FALSE, 
                   annotation_row=tmp %>% column_to_rownames("value"), 
                   annotation_col = tmp2,
                   na_col = 'grey', show_colnames=F,
                   annotation_colors = ann_colors,
                   cellwidth = 1, cellheight = 10,
                   gaps_col = which(colnames(mat) %in% c(gaps-1, 553))
                   
                   ,
                   filename = "shoc2-heatmap.pdf")




#With numbers

jason_labels <- c("223", "177", "131", "203", "200", "178", "63", "64", "65", "66", "411", "434", "154", "316", "244", "288", "133", "156", "316", "293", "129")
mask <- which(!(colnames(mat) %in% jason_labels))

new_colnames <- colnames(mat)

new_colnames[mask] <- ""

library(pheatmap)
pheatmap::pheatmap(mat[tmp$value[1:20],] %>% magrittr::set_colnames(new_colnames),
                   border_color = "black",
                   cluster_rows=FALSE,breaks = seq(-3, 3, length.out=1000),
                   color=colorRampPalette(c("navy",
                                            "white", 
                                            "red"))(1000), 
                   cluster_cols=FALSE, 
                   annotation_row=tmp %>% column_to_rownames("value"), 
                   annotation_col = tmp2,
                   na_col = 'grey', show_colnames=T,fontsize_col = 4,
                   annotation_colors = ann_colors,
                   cellwidth = 1, cellheight = 10,
                   gaps_col = which(colnames(mat) %in% c(gaps-1, 553)),
                   filename = "shoc2_heatmap_w_numbers.pdf")




# Mega heatmap ------------------------------------------------------------

pheatmap::pheatmap(mat[tmp$value[1:20],],
                   border_color = "black",
                   cluster_rows=FALSE,breaks = seq(-3, 3, length.out=1000),
                   color=colorRampPalette(c("navy",
                                            "white", 
                                            "red"))(1000), 
                   cluster_cols=FALSE, 
                   annotation_row=tmp %>% column_to_rownames("value"), 
                   annotation_col = tmp2,
                   na_col = 'grey', show_colnames=F,
                   annotation_colors = ann_colors,
                   cellwidth = 10, cellheight = 10,
                   gaps_col = which(colnames(mat) %in% c(gaps-1, 553))
                   
                   ,
                   filename = "shoc2-heatmap_large.pdf")


# Amino acid level features -----------------------------------------------
aa_groups = list("Neg/Acidic"= c("D", "E"),
     "Pos/Basic"= c("K", "R"),
     "Polar_uncharged"= c("S", "T", "C", "Y", "N", "Q"),
     "Non-polar_non-aromatic"= c("G", "A", "V", "L", "I", "M"),
     #"Proline"= c("P"),
     "Non-polar_large-aromatic"= c("F", "W"),
     "Helix_Breaker"= c("G", "P"),
     "Aromatic"= c("F","W", "Y", "H"),
     "Charged"= c("R", "K", "D", "E", "H"),
     "Aliphatic"= c("P", "A", "L", "V", "I"))

map_df(aa_groups, function(x) colMeans(mat[x,],na.rm = T), .id = "group") %>% 
  pivot_longer(names_to = "Index", values_to = "Score", -group) %>% 
  mutate(Index= as.numeric(Index)) %>% 
  ggplot(aes(Index, Score, group = group, color = group)) + 
  geom_line()+
  facet_wrap(~group, ncol=1)


map_df(aa_groups, function(x) colMeans(mat[x,],na.rm = T), .id = "group") %>% 
  pivot_longer(names_to = "Index", values_to = "Score", -group) %>% 
  mutate(Index= as.numeric(Index)) %>% 
  write_tsv("./output/grouped_aa_mean.tsv")



map(aa_groups, function(x) colMeans(mat[x,],na.rm = T)) %>% 
  do.call(rbind, .) %>% 
  pheatmap(
           border_color = "black",
           cluster_rows=FALSE,breaks = seq(-3, 3, length.out=1000),
           color=colorRampPalette(c("navy",
                                    "white", 
                                    "red"))(1000), 
           cluster_cols=FALSE)
           

# Export for Jason --------------------------------------------------------
output <- File1_reformatted_SHOC2_callapsedToAA %>% left_join(
  
  rescaled %>% 
    as_tibble(rownames = "Vt_aa") %>% 
    pivot_longer(names_to = "POS", values_to = "Josh_P_Rescaled", -Vt_aa) %>% 
    mutate(POS = as.numeric(POS))
)

output %>% 
  filter(Vt_aa == "B") %>% 
  pull("Josh_P_Rescaled") %>% 
  mean()
#result = 0

output %>% 
  filter(Vt_aa == "X") %>% 
  pull("Josh_P_Rescaled") %>% 
  mean()
#result: -1

write_tsv(output,"./output/rescaled-data-for_Jason.tsv")

