library(tidyverse)
library(ggplot2)
library(ggpubr) 
library(ggtext)
library(patchwork)
library(ggstatsplot)

#Define paths
results_path <- "~/Busqueda_antigenos/All_Final_results/AllStrains_AgProtect_Final_results.tsv"
output_path <- "~/Desktop/Graficos"

Results <- read_tsv(results_path)

Results$Host_Homologue_Result_All <- str_replace_all(Results$Host_Homologue_Result_All,
                                                     "Host Homologue and Non Host Homologue, depending on the Host",
                                                     "Proteins with **sequence similarity** to at least one protein in **some hosts**")
Results$Host_Homologue_Result_All <- str_replace_all(Results$Host_Homologue_Result_All,
                                                     "Non Host Homologue",
                                                     "Proteins with **no sequence similarity** to any protein in any host")
Results$Host_Homologue_Result_All <- str_replace_all(Results$Host_Homologue_Result_All,
                                                     "Host Homologue",
                                                     "Proteins with **sequence similarity** to at least one protein in **all hosts**")

Results$Strain <- str_replace_all(Results$Strain,
                                  "Experimental_Antigens",
                                  "Experimental antigens")

# Find total number of proteins per strain
Strain_number <- group_by(Results, Strain) %>% 
  summarise("Total_Number_of_proteins_per_Strain"=n()) 

COG <- Results %>%
  filter(Host_Homologue_Result_All == "Non Host Homologue",
         AntigenicityResult == "ANTIGEN") %>%
  select(ID, COG_category_description, Strain) %>%
  group_by(Strain, COG_category_description) %>%
  summarise(Total_Number_of_proteins_COG = n()) %>%
  full_join(Strain_number, by = "Strain")


Exp_Total <- Results %>%
  group_by(SubcellularLocalization, Strain) %>%
  summarise(Exp_or_NonExp_proteins_total = n())

Exp_Total <- full_join(Exp_Total, Strain_number)

Exp_Total <- Exp_Total %>%
  mutate(Percentage_Exp_Total = Exp_or_NonExp_proteins_total * 100 / Total_Number_of_proteins_per_Strain )

AgProt_Tot <- Exp_Total %>%
  filter(Strain == "Experimental antigens") 

Mean_Percentage_Exp_Total <- Exp_Total %>%
  filter(Strain != "Experimental antigens") %>%
  summarise(mean_percentage = round(mean(Percentage_Exp_Total), digits = 2)) 

##### Grafico Exposition ####
#filtro
Exposed <- filter(Results, AntigenicityResult ==  "ANTIGEN" & 
                   Host_Homologue_Result_All ==  "Proteins with **no sequence similarity** to any protein in any host")

SL_COG <- group_by(Exposed, SubcellularLocalization, Strain, COG_category_description) %>%
  summarise(Exposed_or_Not_proteins_COG=n())

#agrupo
Exposed <- group_by(Exposed, Exposition, SubcellularLocalization, Strain) %>%
  summarise(Exposed_or_Not_proteins=n()) 

#Agrego el total de cada grupo
# Add the total number of non-host homologue proteins per strain
Antigenics <- read_tsv(paste0(output_path,"/Antigenics.tsv"))
Antigenics <- filter(Antigenics, AntigenicityResult == "ANTIGEN")
Exposed <- full_join(Exposed, Antigenics, by = "Strain")
rm(Antigenics)

#acomodo columnas para menos confusion
Exposed <- select(Exposed, SubcellularLocalization, Strain, 
                  Exposed_or_Not_proteins, Ag_or_NonAg_proteins, Exposition)

Exposed <- Exposed %>%
  rename(Total_Ag_proteins=Ag_or_NonAg_proteins)

Exposed<- left_join(Exposed,Strain_number, by = "Strain")

#agrego columna porcentaje
Exposed <- Exposed %>%
  mutate(
Percentage_Exposed_or_Not = Exposed_or_Not_proteins * 100 / Total_Ag_proteins,
Percentage_Total = Exposed_or_Not_proteins * 100 / Total_Number_of_proteins_per_Strain
)

#cambio columna de character a factor
Exposed$SubcellularLocalization = as.factor(Exposed$SubcellularLocalization )

write_tsv(Exposed, paste0(output_path,"/Exposed.tsv"))


# Calculate the mean percentage of non-antigenic proteins
Mean_Percentage_Exposed <- Exposed %>%
  filter(SubcellularLocalization == "OuterMembrane" | SubcellularLocalization == "Extracellular",
         Strain != "Experimental antigens") %>%
  summarise(mean_percentage = round(mean(Percentage_Exposed_or_Not), digits = 2)) 

Mean_Percentage_Exposed <- sum(Mean_Percentage_Exposed$mean_percentage)

#Agrego el total de cada grupo
SL_COG = full_join(SL_COG, Exposed, by=c("Strain", "SubcellularLocalization"))

SL_COG <- SL_COG %>% select(SubcellularLocalization, Strain, COG_category_description,
                            Exposed_or_Not_proteins_COG, Total_Ag_proteins)

#agrego columna porcentaje
SL_COG_Ag <- SL_COG %>% 
  mutate(Exposed_Percentage = Exposed_or_Not_proteins_COG *100 / Total_Ag_proteins) %>% 
  filter(substr(Strain, 1, 3) == "Exp")

SL_COG <- SL_COG %>% 
  mutate(Exposed_Percentage = Exposed_or_Not_proteins_COG *100 / Total_Ag_proteins) %>% 
  filter(substr(Strain, 1, 3) != "Exp")


means_SL <- SL_COG %>%
  group_by(COG_category_description, SubcellularLocalization) %>%
  summarize(mean_value = mean(Exposed_Percentage)) %>%
  na.omit(means_SL)

means_SL <- bind_rows(means_SL, tibble(COG_category_description = "-", 
                                       SubcellularLocalization = "Cell wall",
                                       mean_value = 0))

means_SL$COG_category_description <- str_replace_all(means_SL$COG_category_description,
                                                     "_",
                                                     " ")

means_SL$COG_category_description <- str_replace_all(means_SL$COG_category_description,
                                                     "-",
                                                     "Proteins not assigned to any COG")

means_SL$COG_category_description <- str_replace_all(means_SL$COG_category_description,
                                                     "Posttranslational",
                                                     "Post-translational")

means_SL$SubcellularLocalization <- str_replace_all(means_SL$SubcellularLocalization,
                                                     "CytoplasmicMembrane",
                                                     "Cytoplasmic membrane")

means_SL$SubcellularLocalization <- str_replace_all(means_SL$SubcellularLocalization,
                                                    "OuterMembrane",
                                                    "Outer membrane")

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

#remove NAs, since no protein fit into Cell Wall category
means_SL <- filter(means_SL, !is.na(COG_category_description))

means_SL$SubcellularLocalization <- as.factor(means_SL$SubcellularLocalization)

p4B <- ggplot(means_SL, aes(fill=factor(SubcellularLocalization, 
                                       levels=c("Unknown", "Cytoplasmic", "Cytoplasmic membrane", 
                                                "Periplasmic","Outer membrane",
                                                "Extracellular","Cell wall")), 
                           y=mean_value, x=fct_rev(factor(COG_category_description, 
                                                          levels=COG_category_description_order)))) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values =c("#e1ccd1","#9e2f28","sienna2" , "#deb867",
                             "#25482f","olivedrab","#4f2f4a","#4d0f0b")) +
  theme(axis.title = element_markdown(face="bold"), 
                 axis.text = element_text(family = "Times New Roman", size = 10, color = "black"),
                 title = element_text(family = "Times New Roman"), 
                 legend.position= "bottom",
                 legend.location = "plot",
                 axis.title.x.bottom = element_markdown(face="bold"), 
                 legend.byrow = T,
                 legend.title = element_blank(),
                 text = element_text(family = "Times New Roman", size = 14), 
                 panel.background =  element_rect(fill = "white"), 
                 panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                                 linewidth = 0.3)) + 
           labs(y = "Mean % of *Av. paragallinarum*<br> antigenic proteins with no sequence<br>similarity to host proteins", 
                x = "COG category",
                tag = "B") +
           coord_flip() + 
  ylim(0, 40)
         

p4B

Exposed <- Exposed %>% 
  mutate(
    SL_Percent = Exposed_or_Not_proteins * 100 / Total_Ag_proteins,
    Group = ifelse(Strain == "Experimental antigens",
                   "Experimental antigens",
                   "*Av. paragallinarum* strains")
  )

Exposed <- Exposed %>%
  mutate(SubcellularLocalization = recode(SubcellularLocalization,
                                          "Cellwall" = "Cell wall", 
                                          "Cytoplasmic" = 'Cytoplasmic',
                                          "CytoplasmicMembrane" = 'Cytoplasmic\nmembrane', 
                                          "Extracellular" = 'Extracellular',
                                          "OuterMembrane" = 'Outer\nmembrane',
                                          "Periplasmic" = 'Periplasmic', 
                                          "Unknown" = 'Unknown'))

Exposed$SubcellularLocalization <- factor(Exposed$SubcellularLocalization, 
                                          levels=c("Unknown", "Cytoplasmic", "Cytoplasmic\nmembrane", 
                                                   "Periplasmic","Outer\nmembrane",
                                                   "Extracellular","Cell wall"))

AgProt <- Exposed %>%
  filter(Strain == "Experimental antigens")

Exposed <- Exposed %>%
  filter(Strain != "Experimental antigens")


p4A <- ggplot(Exposed, aes(x = SubcellularLocalization, y = SL_Percent)) +
  geom_violin(width = 1.4, aes(fill = SubcellularLocalization), alpha= 0.4) +
  scale_fill_manual(values =c("#e1ccd1","#9e2f28","sienna2" , "#deb867",
                              "#25482f","olivedrab","#4f2f4a","#4d0f0b")) +
  geom_point(data=AgProt, aes(x=SubcellularLocalization, y=SL_Percent, 
                              color=SubcellularLocalization), 
             size = 7, shape = 18) +
  scale_color_manual(values =c("#e1ccd1","#9e2f28","sienna2" , "#deb867",
                               "#25482f","olivedrab","#4f2f4a","#4d0f0b")) +
  theme(
    text = element_text(family = "Times New Roman", size = 14, color = "black"), 
    axis.title = element_markdown(),
    legend.position = "none",
    plot.caption = element_markdown(halign = 0, hjust = 0.95),
    legend.title = element_blank(),
    plot.tag.location = "margin",
    axis.title.y.left = element_markdown(face="bold"), 
    axis.text.y = element_text(family = "Times New Roman", size = 14, color = "black"),
    axis.title.y = element_text(family = "Times New Roman", size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(family = "Times New Roman", size = 14, color = "black",face = "bold"),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                    linewidth = 0.3)
  ) + 
  labs(
    y = "% of antigenic proteins with no<br>sequence similarity to host<br>in each localisation",
    tag = "A",
    x = "",
    caption = "Violin plots: *Av. paragallinarum* strains<br>◆: Experimental Antigens  "
  ) + 
  ylim(0, 60)

p4A

means_SL <- SL_COG_Ag %>%
  group_by(COG_category_description, SubcellularLocalization) %>%
  summarize(mean_value = mean(Exposed_Percentage)) %>%
  na.omit(means_SL)

means_SL$COG_category_description <- str_replace_all(means_SL$COG_category_description,
                                                     "_",
                                                     " ")

means_SL$COG_category_description <- str_replace_all(means_SL$COG_category_description,
                                                     "-",
                                                     "Proteins not assigned to any COG")

means_SL$COG_category_description <- str_replace_all(means_SL$COG_category_description,
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

means_SL <- means_SL %>%
  mutate(SubcellularLocalization = recode(SubcellularLocalization,
                                          "Cellwall" = "Cell wall", 
                                          "Cytoplasmic" = 'Cytoplasmic',
                                          "CytoplasmicMembrane" = 'Cytoplasmic\nmembrane', 
                                          "Extracellular" = 'Extracellular',
                                          "OuterMembrane" = 'Outer\nmembrane',
                                          "Periplasmic" = 'Periplasmic', 
                                          "Unknown" = 'Unknown'))

means_SL$SubcellularLocalization <-  factor(means_SL$SubcellularLocalization, 
                                          levels=c("Unknown", "Cytoplasmic", "Cytoplasmic\nmembrane", 
                                                   "Periplasmic","Outer\nmembrane",
                                                   "Extracellular","Cell wall"))

p4C <- ggplot(means_SL, aes(fill=SubcellularLocalization,
                            y=mean_value, x=fct_rev(factor(COG_category_description, 
                                                           levels=COG_category_description_order)))) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values =c("#e1ccd1","#9e2f28","sienna2" , "#deb867",
                              "#25482f","olivedrab","#4f2f4a","#4d0f0b")) +
  theme(axis.title = element_markdown(face="bold"), 
        axis.text = element_text(family = "Times New Roman", size = 10, color = "black"),
        title = element_text(family = "Times New Roman"), 
        legend.position= "none",      
        axis.title.x.bottom = element_markdown(face="bold"), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        text = element_text(family = "Times New Roman", size = 14), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "% of experimental antigens<br>(antigenic, with no sequence<br>similarity to host proteins)", 
       x = "COG category",
       tag = "C") +
  coord_flip() + 
  ylim(0, 40)

p4C


layout <- '
AAAA
AAAA
AAAA
AAAA
AAAA
BBCC
BBCC
BBCC
BBCC
BBCC
DD##
'
wrap_plots(A = free(p4A), B = p4B, C= p4C, D = guide_area(), design = layout)  +
  plot_layout(guides = 'collect')

ggsave("4A&B.png", device = "png", path = output_path, 
       width =2900, height = 3000, units="px")

