library(ggraph)
library(tidygraph)


# Begin -------------------------------------------------------------------

library(ProjectTemplate)
load.project()

#Jason - once your file is correct, save it here.
#par_cord_complex_example <- X21.11.27CORRECT_FINAL.cryoem.energy.bond.calc

par_cord_complex_example
#The key variables are as follows:
#interaction - indicates which protein interactions the  are
#POS - is the amino acid position of SHOC2. I left POS so that its easier to merge other relevant data such as the LFC from DMS data.
#MRAS - indicates position of MRAS residue in interaction
#PP1C - indicates position of PP1C residue in interaction
#Type - indicates the type of interaction; I was hoping we could indicate this by having the edges dashed for D (D is likely just hydrophobic interaction) and solid for all others (polar contacts) etc...
#Energy - this is the energy of each pairwise interaction. I was hoping we could have the thickness of the edge be based on this maybe?
  
  
#Reshape Jason's data

#Split table into three parts
tmp0 <- per_variant_LFC_SHOC2_DMS %>% 
  group_by(POS) %>% 
  dplyr::summarize(mean_z = mean(LFC_ct_frctn_adj.z))

tmp <- par_cord_complex_example

split_tmp <- split(tmp, tmp$interaction)


#Rewrite the from to columns based on the splits
shoc2_mras <- split_tmp$`SHOC2-MRAS` %>% 
  left_join(tmp0) %>% 
  mutate(from = paste("SHOC2", POS, sep = '_'),
         to = paste("MRAS", MRAS, sep = '_')) %>% 
  select(from, to, energy = Energy, type = Type, mean_z)

shoc2_pp1c <- split_tmp$`SHOC2-PP1C` %>% 
  left_join(tmp0) %>% 
  mutate(from = paste("SHOC2", POS, sep = '_'),
         to = paste("PP1C", PP1C, sep = '_')) %>% 
  select(from, to, energy = Energy, type = Type, mean_z)


mras_pp1c <- split_tmp$`MRAS-PP1C` %>% 
  mutate(from = paste("MRAS", MRAS, sep = '_'),
         to = paste("PP1C", PP1C, sep = '_'),
         mean_z = NA) %>% 
  select(from, to, energy = Energy, type = Type, mean_z)


edge_table <- rbind(shoc2_mras, shoc2_pp1c) %>% 
  mutate(energy = exp(energy),
         type = type == "D")

#Edges
#Type (hydrophobic or polar) - dotted or dashed
#Energy of each pairwise interaction - thickness
#Z score from fitness screen - color
#from: protein_#
#to: protein_#

#Nodes
#Synthesize
#Name: SHOC2_#
#Degree: equiv to POS, 1-length for each protein
#Protein: SHOC2.

shoc2_len <- 582
mras_len <- 208
pp1c_len <- 323


generate_node_table <- function(prot_name, prot_len){
  
  tibble(prot = rep(prot_name, prot_len),
         y = seq_len(prot_len),
         node_key = paste(prot, y, sep = '_')) %>% 
    select(node_key, everything())
}


#Idea: custom layout 
shoc2_x <- 0
mras_x <- -100
pp1c_x <- 100


node_table <- rbind(
  generate_node_table("SHOC2", shoc2_len) %>% mutate(x = shoc2_x),
  generate_node_table("MRAS", mras_len) %>% mutate(x = mras_x),
  generate_node_table("PP1C", pp1c_len) %>% mutate(x = pp1c_x))


#Flip y

node_table <- node_table %>% mutate(y=-y)



# Construct table graph ---------------------------------------------------




g <- tbl_graph(nodes = node_table, edges = edge_table)

ggraph(g, x=x, y=y)+
  geom_edge_link(aes(colour = mean_z, linetype = type), width =1) + #alpha = energy
  geom_node_point(aes(colour = prot), size = 0.5) + 
  coord_fixed() +
  scale_edge_color_gradient2(low = 'blue', high = 'red', mid='white',midpoint = 0) +
  theme_void() + 
  scale_color_manual(values = c("SHOC2" = "#F47721",
                                "MRAS" = "#B91D7B",
                                "PP1C" = "#00A783"))









