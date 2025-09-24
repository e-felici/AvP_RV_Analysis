# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(ggpubr) 
library(patchwork)
library(ggtext)
library(ggstatsplot)
library(taylor)


#Define paths
results_path <- "~/Busqueda_antigenos/All_Final_results/AllStrains_AgProtect_Final_results.tsv"
output_path <- "~/Desktop/Graficos"

Results <- read_tsv(results_path)

#Remove underscore
Results$Strain <- str_replace_all(Results$Strain,
                                  "Experimental_Antigens",
                                  "Experimental Antigens")


Results$Host_Homologue_Result_All <- str_replace_all(Results$Host_Homologue_Result_All,
                                                  "Non Host Homologue",
                                                  "Non-Host Homologue")

# Find total number of proteins per strain
Strain_number = group_by(Results, Strain) %>% 
  summarise("Total_Number_of_proteins_per_Strain"=n()) 

# Select relevant columns and group by Strain and COG, then summarize
COG <- Results %>%
  select(ID, COG, Strain) %>%
  group_by(Strain, COG) %>%
  summarise(Total_Number_of_proteins_COG = n()) %>%
  full_join(Strain_number, by = "Strain")

##### Figure 2A
# Create and polish Host Homologue tibble
Host <- Results %>%
  group_by(Host_Homologue_Result_All, Strain) %>%
  summarise(Host_Homologue_or_NonHH_Proteins = n()) %>%
  full_join(Strain_number, by = "Strain") %>%
  mutate(Percentage_Host = Host_Homologue_or_NonHH_Proteins * 100 / Total_Number_of_proteins_per_Strain)

write_tsv(Host, paste0(output_path,"/Host.tsv"))

# Calculate the mean percentage of non-host homologue proteins
Mean_Percentage_Host <- Host %>%
  filter(Host_Homologue_Result_All == "Non-Host Homologue", 
         Strain != "Experimental Antigens") %>%
  summarise(mean_percentage = round(mean(Percentage_Host), digits = 2)) %>%
  pull(mean_percentage)


#######Figure 2B
#Create Host and COG tibble
Host_COG <- Results %>%
  group_by(Host_Homologue_Result_All, Strain, COG_category_description) %>%
  summarise(Total_Number_of_proteins_Host_COG = n()) %>%
  full_join(Strain_number, by = c("Strain")) %>%
  mutate(Percentage_Host_COG = Total_Number_of_proteins_Host_COG * 100 / Total_Number_of_proteins_per_Strain) %>%
  filter(!str_starts(Strain, "Exp"))

# Calculate the mean percentage of host homologue proteins per COG category
means_Host_COG <- Host_COG %>%
  group_by(COG_category_description, Host_Homologue_Result_All) %>%
  summarize(mean_value = mean(Percentage_Host_COG))

means_Host_COG$COG_category_description <- str_replace_all(means_Host_COG$COG_category_description,
                                  "_",
                                  " ")

means_Host_COG$COG_category_description <- str_replace_all(means_Host_COG$COG_category_description,
                                                           "-",
                                                           "Proteins not assigned to any COG")

means_Host_COG$COG_category_description <- str_replace_all(means_Host_COG$COG_category_description,
                                                           "Posttranslational",
                                                           "Post-translational")

COG_category_description_order <- c("Translation, ribosomal structure and biogenesis", 
                                    "RNA processing and modification", 
                                    "Transcription", 
                                    "Replication, recombination and repair", 
                                    "Chromatin structure and dynamics", 
                                    "Cell cycle control, cell division, chromosome partitioning", 
                                    "Nuclear structure", 
                                    "Defense mechanisms", 
                                    "Signal transduction mechanisms", 
                                    "Cell wall/membrane/envelope biogenesis", 
                                    "Cell motility",
                                    "Cytoskeleton", 
                                    "Extracellular structures", 
                                    "Intracellular trafficking, secretion, and vesicular transport", 
                                    "Post-translational modification, protein turnover, chaperones", 
                                    "Mobilome: prophages, transposons", 
                                    "Energy production and conversion", 
                                    "Carbohydrate transport and metabolism", 
                                    "Amino acid transport and metabolism", 
                                    "Nucleotide transport and metabolism", 
                                    "Coenzyme transport and metabolism", 
                                    "Lipid transport and metabolism", 
                                    "Inorganic ion transport and metabolism", 
                                    "Secondary metabolites biosynthesis, transport and catabolism", 
                                    "General function prediction only", 
                                    "Function unknown", 
                                    "Proteins not assigned to any COG")

all_levels <- c("Host Homologue", 
                "Host Homologue or Non-Host Homologue, depending on the Host",
                "Non-Host Homologue")

means_Host_COG <- means_Host_COG %>%
  bind_rows(tibble(
    Host_Homologue_Result_All = "Host Homologue or Non-Host Homologue, depending on the Host",
    mean_value = 0,
    COG_category_description = "Proteins not assigned to any COG"
  ))

means_Host_COG <- means_Host_COG %>%
  mutate(Host_Homologue_Result_All = factor(Host_Homologue_Result_All, 
                                            levels = all_levels)
         )

# Create the second plot (Figure 2B)
p2B <- ggplot(means_Host_COG, aes(fill=Host_Homologue_Result_All, 
                                  y=mean_value, 
                                  x=fct_rev(factor(COG_category_description, 
                                                   levels=COG_category_description_order))
                                    )) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values =c("Host Homologue"="#9e2f28", 
                              "Host Homologue or Non-Host Homologue, depending on the Host"="sienna2",
                              "Non-Host Homologue"="olivedrab"),
                    drop = FALSE) +
  theme(axis.title = element_markdown(face="bold"), 
        axis.text = element_text(family = "Times New Roman", size = 10, color = "black"),
        title = element_text(family = "Times New Roman"), 
        legend.position= "bottom",
  #      legend.justification = c("left", "bottom"),
 #       legend.box.just = "left",
#        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        text = element_text(family = "Times New Roman", size = 14), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "Mean % of Total Proteins in<br>*Av. paragallinarum* Strains", 
       x = "COG category",
       tag = "B") +
  coord_flip() + 
  ylim(0, 25)

p2B


Host <- Host %>% 
  filter(Host_Homologue_Result_All == "Non-Host Homologue") %>%
  mutate(
    HH_Percent = Host_Homologue_or_NonHH_Proteins * 100 / Total_Number_of_proteins_per_Strain,
    Group = ifelse(Strain == "Experimental Antigens",
                   "Experimental Antigens",
                   "*Av. paragallinarum* strains")
  )

AgProt <- Host %>% 
  ungroup() %>%
  filter(Strain== "Experimental Antigens") %>%
  select(HH_Percent) %>%
  as.double()

Host <- Host %>%
  filter(Strain != "Experimental Antigens")

p2A <- ggbetweenstats(
  data = Host,
  x = Group,
  y = HH_Percent,
  centrality.label.args = list(color = "transparent", size=0, box.padding=200),  
  centrality.point.args = list(color = "#25482f", size = 6),
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.6),
    alpha = 0.7, 
    size = 3, 
    stroke = 0
  ),
  ggplot.component = list(
  #  ggplot2::scale_y_continuous(
   #   limits = c(70,85),
    #  breaks = seq(from = 70, to = 85, by = 5)
  #  ),
    scale_color_manual(values = c(
      "*Av. paragallinarum* strains" = "olivedrab", 
      "Experimental Antigens" = "#0073C2FF"
    )),
#    scale_y_continuous(breaks = seq(70, 85, 5), limits = c(70, 85)),
    theme(
      text = element_text(family = "Times New Roman", size = 14, color = "black"), 
      axis.text = element_text(family = "Times New Roman", size = 14, color = "black"),
      legend.position = "bottom",
      #legend.position.inside = c(0.5, 0.005),
      legend.box = "vertical",
      legend.direction = "vertical",
      legend.text = ggtext::element_markdown(),
      legend.title = element_blank(),
      axis.text.x = element_blank(),
      panel.background = element_rect(fill = "white"), 
      panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                      linewidth = 0.3),
      panel.grid.minor = element_blank()
    ),
    labs(
      y = "Percentage of Non-Host Homologue Proteins",
      tag = "A",
      x = ""
    )
  )
) + 
  geom_point(aes(y = AgProt, color = "Experimental Antigens"), size = 5, shape = 18)

p2A

#######Figure 2B
#Create Host and COG tibble
Host_COG_Ag <- Results %>%
  group_by(Host_Homologue_Result_All, Strain, COG_category_description) %>%
  summarise(Total_Number_of_proteins_Host_COG = n()) %>%
  full_join(Strain_number, by = c("Strain")) %>%
  mutate(Percentage_Host_COG = Total_Number_of_proteins_Host_COG * 100 / Total_Number_of_proteins_per_Strain) %>%
  filter(str_starts(Strain, "Exp"))

# Calculate the mean percentage of host homologue proteins per COG category
means_Host_COG <- Host_COG_Ag %>%
  group_by(COG_category_description, Host_Homologue_Result_All) %>%
  summarize(mean_value = mean(Percentage_Host_COG))

means_Host_COG$COG_category_description <- str_replace_all(means_Host_COG$COG_category_description,
                                                           "_",
                                                           " ")

means_Host_COG$COG_category_description <- str_replace_all(means_Host_COG$COG_category_description,
                                                           "-",
                                                           "Proteins not assigned to any COG")

means_Host_COG$COG_category_description <- str_replace_all(means_Host_COG$COG_category_description,
                                                           "Posttranslational",
                                                           "Post-translational")

COG_category_description_order <- c("Translation, ribosomal structure and biogenesis", 
                                    "RNA processing and modification", 
                                    "Transcription", 
                                    "Replication, recombination and repair", 
                                    "Chromatin structure and dynamics", 
                                    "Cell cycle control, cell division, chromosome partitioning", 
                                    "Nuclear structure", 
                                    "Defense mechanisms", 
                                    "Signal transduction mechanisms", 
                                    "Cell wall/membrane/envelope biogenesis", 
                                    "Cell motility",
                                    "Cytoskeleton", 
                                    "Extracellular structures", 
                                    "Intracellular trafficking, secretion, and vesicular transport", 
                                    "Post-translational modification, protein turnover, chaperones", 
                                    "Mobilome: prophages, transposons", 
                                    "Energy production and conversion", 
                                    "Carbohydrate transport and metabolism", 
                                    "Amino acid transport and metabolism", 
                                    "Nucleotide transport and metabolism", 
                                    "Coenzyme transport and metabolism", 
                                    "Lipid transport and metabolism", 
                                    "Inorganic ion transport and metabolism", 
                                    "Secondary metabolites biosynthesis, transport and catabolism", 
                                    "General function prediction only", 
                                    "Function unknown", 
                                    "Proteins not assigned to any COG")

means_Host_COG$Host_Homologue_Result_All <- str_replace_all(means_Host_COG$Host_Homologue_Result_All,
                                                           "and",
                                                           "or")

means_Host_COG <- means_Host_COG %>%
  mutate(Host_Homologue_Result_All = factor(Host_Homologue_Result_All, 
                                            levels = all_levels)
  )

p2C <- ggplot(means_Host_COG, aes(fill=Host_Homologue_Result_All, 
                                  y=mean_value, 
                                  x=fct_rev(factor(COG_category_description, 
                                                   levels=COG_category_description_order))
)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values =c("Host Homologue"="#9e2f28", 
                              "Host Homologue or Non-Host Homologue, depending on the Host"="sienna2",
                              "Non-Host Homologue"="olivedrab"),
                    drop = FALSE) +
  theme(axis.title.x = element_markdown(face="bold"), 
        axis.title.y = element_blank(),
        axis.text = element_text(family = "Times New Roman", size = 10, color = "black"),
        axis.text.y = element_blank(),
        title = element_text(family = "Times New Roman"), 
        legend.position = "bottom",
        #legend.position.inside = c(0.75, 0.9),
        legend.box = "vertical",
        legend.box.just = "bottom",
        legend.direction = "vertical",
        legend.title = element_blank(),
        text = element_text(family = "Times New Roman", size = 14), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "% of Total<br>Experimental Antigens", 
       x = "COG category",
       tag = "C") +
  coord_flip() + 
  ylim(0, 25)

p2C

layout <- '
AAABBBBBBBB
'

patch <- p2A + plot_layout(guides = 'keep') 
patch2 <- (p2B | p2C) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "bottom",
        legend.direction = "vertical",
        legend.title = element_blank())

wrap_plots(A = patch, B = patch2, design = layout)

# Save the combined plot as a PNG file
ggsave("2A&B.png", device = "png", path = output_path, 
       width =3700, height = 2200, units="px")


