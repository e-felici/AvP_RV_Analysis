# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(ggpubr) 
library(ggtext)

#Define paths
results_path <- "~/Busqueda_antigenos/All_Final_results/AllStrains_AgProtect_Final_results.tsv"
output_path <- "~/Desktop/Graficos"

Results <- read_tsv(results_path, col_names = T)

#Remove underscore
Results$Strain <- str_replace_all(Results$Strain,
                                  "Experimental_Antigens",
                                  "Experimental Antigens")

# Find total number of proteins per strain
Strain_number = group_by(Results, Strain) %>% 
  summarise("Total_Number_of_proteins_per_Strain"=n()) 

COG <- Results %>%
  filter(Host_Homologue_Result_All == "Non Host Homologue") %>%
  select(ID, COG_category_description, Strain) %>%
  group_by(Strain, COG_category_description) %>%
  summarise(Total_Number_of_proteins_COG = n()) %>%
  full_join(Strain_number, by = "Strain")

Ag_Total <- Results %>%
  group_by(AntigenicityResult, Strain) %>%
  summarise(Ag_or_NonAg_proteins_total = n())

Ag_Total <- full_join(Ag_Total, Strain_number)

Ag_Total <- Ag_Total %>%
  mutate(Percentage_Ag_Total = Ag_or_NonAg_proteins_total * 100 / Total_Number_of_proteins_per_Strain )

AgProt_Tot <- Ag_Total %>%
  filter(AntigenicityResult == "ANTIGEN", 
         Strain == "Experimental Antigens") %>%
  summarise(mean_percentage = round(mean(Percentage_Ag_Total), digits = 2)) %>%
  pull(mean_percentage)

# Filter non-host homologue proteins
Antigenics <- filter(Results, Host_Homologue_Result_All == "Non Host Homologue")

Ag_COG <- Antigenics %>%
  group_by(AntigenicityResult, Strain, COG_category_description) %>%
  summarise(Ag_or_NonAg_proteins_COG = n())

Antigenics <- Antigenics %>%
  group_by(AntigenicityResult, Strain) %>%
  summarise(Ag_or_NonAg_proteins = n())

# Add the total number of non-host homologue proteins per strain
Host <- read_tsv(paste0(output_path,"/Host.tsv"))
Host <- filter(Host, Host_Homologue_Result_All == "Non Host Homologue")
Antigenics <- full_join(Antigenics, Host, by = "Strain")
rm(Host)

# Rename columns for clarity
Antigenics <- Antigenics %>%
  rename(Total_Non_HH_proteins = Host_Homologue_or_NonHH_Proteins )

# Remove unnecessary column
Antigenics <- select(Antigenics, -Host_Homologue_Result_All, 
                     -Total_Number_of_proteins_per_Strain, -Percentage_Host)

# Add percentage columns
Antigenics <- Antigenics %>%
  mutate(
    Percentage_Antigenic_or_Not = Ag_or_NonAg_proteins * 100 / Total_Non_HH_proteins,
    Percentage_antigenic_Total = Ag_or_NonAg_proteins * 100 / Strain_number$Total_Number_of_proteins_per_Strain
  )

write_tsv(Antigenics, paste0(output_path,"/Antigenics.tsv"))

# Add the total number of proteins per strain and COG category
Ag_COG <- full_join(Ag_COG, Antigenics, by = c("Strain", "AntigenicityResult"))

Ag_COG <- Ag_COG %>% select(AntigenicityResult, Strain, COG_category_description,
                            Ag_or_NonAg_proteins_COG, Total_Non_HH_proteins)

# Add percentage column
Ag_COG_AgProt <- Ag_COG %>%
  mutate(Ag_Percentage = Ag_or_NonAg_proteins_COG * 100 / Total_Non_HH_proteins) %>%
  filter(str_starts(Strain, "Exp"))

Ag_COG <- Ag_COG %>%
  mutate(Ag_Percentage = Ag_or_NonAg_proteins_COG * 100 / Total_Non_HH_proteins) %>%
  filter(!str_starts(Strain, "Exp"))

# Calculate the mean percentage of antigenic proteins per COG category
means_ag <- Ag_COG %>%
  group_by(COG_category_description, AntigenicityResult) %>%
  summarize(mean_value = mean(Ag_Percentage))

means_ag$COG_category_description <- str_replace_all(means_ag$COG_category_description,
                                                     "_",
                                                     " ")

means_ag$COG_category_description <- str_replace_all(means_ag$COG_category_description,
                                                     "-",
                                                     "Proteins not assigned to any COG")

means_ag$COG_category_description <- str_replace_all(means_ag$COG_category_description,
                                                     "Posttranslational",
                                                     "Post-translational")

means_ag$AntigenicityResult <- str_replace_all(means_ag$AntigenicityResult,
                                               "NON-ANTIGEN",
                                               "Non antigenic")
means_ag$AntigenicityResult <- str_replace_all(means_ag$AntigenicityResult,
                                               "ANTIGEN",
                                               "Antigenic")

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

# Create the second plot (Figure 2B)
p3B <- ggplot(means_ag, aes(fill=AntigenicityResult, 
                            y=mean_value, 
                            x=fct_rev(factor(COG_category_description, 
                                             levels=COG_category_description_order))
)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values =c("olivedrab","#9e2f28")) +
  theme(axis.title = element_markdown(face="bold"), 
        axis.text = element_text(family = "Times New Roman", size = 10, color = "black"),
        title = element_text(family = "Times New Roman"), 
        legend.position= "none",
        legend.title = element_blank(),
        text = element_text(family = "Times New Roman", size = 14), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "Mean % of Non Homologous<br>to Host Proteins in *Av.<br>paragallinarum* Strains", 
       x = "COG category",
       tag = "B") +
  coord_flip() + 
  ylim(0, 25)

p3B

Antigenics <- Antigenics %>% 
  filter(AntigenicityResult == "ANTIGEN") %>%
  mutate(
    Ag_Percent = Ag_or_NonAg_proteins * 100 / Total_Non_HH_proteins,
    Group = ifelse(Strain == "Experimental Antigens",
                   "Experimental Antigens",
                   "*Av. paragallinarum* strains")
  )

AgProt <- Antigenics %>% 
  ungroup() %>%
  filter(Strain== "Experimental Antigens") %>%
  select(Ag_Percent) %>%
  as.double()

Antigenics <- Antigenics %>%
  filter(Strain != "Experimental Antigens")

p3A <- ggbetweenstats(
  data = Antigenics,
  x = Group,
  y = Ag_Percent,
  centrality.label.args = list(color = "transparent", size=0, box.padding=200),  
  centrality.point.args = list(color = "#25482f", size = 6),
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.6),
    alpha = 0.7, 
    size = 3, 
    stroke = 0
  ),
  ggplot.component = list(
    scale_color_manual(values = c(
      "*Av. paragallinarum* strains" = "olivedrab", 
      "Experimental Antigens" = "#0073C2FF"
    )),
    theme(
      axis.text = ggtext::element_markdown(
        family = "Times New Roman", size = 14, color = "black",
        fill = "white", box.colour = "white", linetype = "solid",
        linewidth = 1, halign = 0.5, valign = 0.5, 
        padding = margin(2, 4, 2, 4), r = unit(5, "pt"),
        align_widths = TRUE, align_heights = TRUE, 
        rotate_margins = TRUE
      ),
      axis.title = ggtext::element_markdown(family = "Times New Roman",size = 14),
      legend.text = ggtext::element_markdown(family = "Times New Roman",size = 11),
      legend.title = element_blank(),
      axis.text.x = element_blank(),
      panel.background = element_rect(fill = "white"), 
      panel.grid.major = element_line(colour = "grey", linetype = "dotted", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.direction = "vertical"
    ),
    labs(
      y = "Percentage of Antigenic Proteins<br>(Non-Homologous to Host)",
      tag = "A",
      x = ""
    )
  )
) + 
  geom_point(aes(y = AgProt, color = "Experimental Antigens"), size = 5, shape = 18)


p3A

means_ag <- Ag_COG_AgProt %>%
  group_by(COG_category_description, AntigenicityResult) %>%
  summarize(mean_value = mean(Ag_Percentage))

means_ag$COG_category_description <- str_replace_all(means_ag$COG_category_description,
                                                     "_",
                                                     " ")

means_ag$COG_category_description <- str_replace_all(means_ag$COG_category_description,
                                                     "-",
                                                     "Proteins not assigned to any COG")

means_ag$COG_category_description <- str_replace_all(means_ag$COG_category_description,
                                                     "Posttranslational",
                                                     "Post-translational")

means_ag$AntigenicityResult <- str_replace_all(means_ag$AntigenicityResult,
                                               "NON-ANTIGEN",
                                               "Non antigenic")
means_ag$AntigenicityResult <- str_replace_all(means_ag$AntigenicityResult,
                                               "ANTIGEN",
                                               "Antigenic")

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

# Create the second plot (Figure 2B)
p3C <- ggplot(means_ag, aes(fill=AntigenicityResult, 
                            y=mean_value, 
                            x=fct_rev(factor(COG_category_description, 
                                             levels=COG_category_description_order))
))  + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values =c("olivedrab","#9e2f28")) +
  theme(axis.title.x = element_markdown(face="bold"), 
        axis.title.y = element_blank(),
        axis.text = element_text(family = "Times New Roman", size = 10, color = "black"),
        axis.text.y = element_blank(),
        title = element_text(family = "Times New Roman"), 
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "bottom",
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(family = "Times New Roman"),
        text = element_text(family = "Times New Roman", size = 14), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "Mean % of Non<br>Homologous to Host<br>Experimental Antigens", 
       x = "COG category",
       tag = "C") +
  coord_flip() + 
  ylim(0, 25)
p3C

layout <- '
AAABBBBCCCC
AAABBBBCCCC
AAABBBBCCCC
AAABBBBCCCC
DDDBBBBCCCC
'
wrap_plots(A = p3A, B = p3B, C= p3C, D = guide_area(), design = layout)  +
  plot_layout(guides = 'collect')

ggsave("3A&B.png", device = "png", path = output_path, 
       width =3500, height = 2000, units="px")
