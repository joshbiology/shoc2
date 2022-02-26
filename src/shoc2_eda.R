library(ProjectTemplate)
library(tidyverse)
load.project()

eve_SHOC2_HUMAN <- eve_SHOC2_HUMAN %>% 
  unite(variant.by.aa, c(wt_aa, position, mt_aa), sep = "")


dplyr::inner_join(per_variant_LFC_SHOC2_DMS, eve_SHOC2_HUMAN)  %>% 
  ggplot(aes(abs(LFC_ct_frctn_adj.z), EVE_scores_ASM)) + 
  geom_point() + 
  geom_smooth(method = "lm")

#Variant table

var_table <- tribble(
  ~Class, ~variant.by.aa,
  "known",   "S2G",
  "known",   "M173I",
  "personal",   "Q61K",
  "personal",   "G63R",
  "personal",   "T411A",
  "personal",   "M173V",
  "personal",   "Q249K"
)

jk_guess <-left_join(var_table, combo)

jk_guess %>% 
  ggplot(aes(LFC_ct_frctn_adj.z, EVE_scores_ASM, label = variant.by.aa, color = Class)) +
  geom_point() + 
  geom_text()

write_tsv(jk_guess, "./output/jason_annotated_vars.tsv")



#plddt

af2 <- read_table('./data/SHOC2Human_d1190_unrelaxed_model_1_rank_1.pdb', skip = 1) %>% 
  select(6, 11) %>% 
  magrittr::set_colnames(c("POS", "pLDDT")) %>% 
  group_by(POS) %>% 
  summarize(pLDDT = median(pLDDT))


#Data annotation
#Correct_LRRnum_CIRCpos where LRR_xtal column has this info. It would be cool to also layer in all the other data that's available in this folder as well! Thanks again for all your help!



combo <- dplyr::inner_join(per_variant_LFC_SHOC2_DMS, eve_SHOC2_HUMAN) %>% 
  group_by(POS) %>% 
  summarise(mean_z = median(LFC_ct_frctn_adj.z),
            mean_eve = median(EVE_scores_ASM, na.rm =T)) %>% 
  left_join(Correct_LRRnum_CIRCpos) %>% 
  left_join(Aminode_accessed_2021.11.23 %>% select(POS, sub_score = Relative.Substitution.Score)) %>% 
  left_join(af2)


#scatter
combo %>% 
  ggplot(aes(mean_z, mean_eve)) + 
  geom_point()


#colplot
#            #data_domain = as.factor(LRR_xtal),

combo %>% 
  select(POS, mean_z, mean_eve, sub_score, pLDDT) %>% 
  pivot_longer(!POS, names_to = "Dataset", values_to = "Value") %>% 
  ggplot(aes(POS, Value)) +
  geom_col() + 
  facet_wrap(~Dataset,ncol = 1, scales = 'free_y')


#Amino acid features



