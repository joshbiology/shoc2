library(drawProteins)
library(ProjectTemplate)
load.project()
# accession numbers of rel A

  

#shoc2_data <- drawProteins::feature_to_dataframe(drawProteins::get_features("Q9UQ13"))

#write_tsv(shoc2_data, "./output/uniprot_shoc2_features.tsv")

#add extra features

shoc2_data <- read_tsv('./output/uniprot_shoc2_features.tsv')
p <- draw_canvas(shoc2_data) 
p <- draw_chains(p, shoc2_data)
p <- draw_domains(p, shoc2_data)
p <- p +ggplot2::geom_rect(data=  shoc2_data[shoc2_data$type == "REPEAT",],
                           aes(xmin=begin,
                               xmax=end,
                               ymin=order-0.25,
                               ymax=order+0.25,
                               colour = description,
                               fill = description)) 


p


p <- p + theme_bw(base_size = 20) + # white backgnd & change text size
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank())

p



p <- draw_repeat(p, rel_data, fill = , label_repeats = F)
p

p <- p + theme_bw(base_size = 20) + # white background
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank())
p


draw_regions(p, rel_data) # adds activation domain


draw_repeat(p, rel_data, fill = desecription) # doesn't add anything in this case


draw_motif(p, rel_data) # adds 9aa Transactivation domain & NLS


draw_phospho(p, rel_data, size = 8)


draw_canvas(rel_data) -> p
p <- draw_chains(p, rel_data)
p <- draw_domains(p, rel_data)
p <- draw_regions(p, rel_data)
p <- draw_motif(p, rel_data)
p <- draw_phospho(p, rel_data, size = 8) 

p <- p + theme_bw(base_size = 20) + # white backgnd & change text size
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank())
p
