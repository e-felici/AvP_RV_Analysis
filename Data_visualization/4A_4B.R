library(tidyverse)
library(ggplot2)
library(ggpubr) 
library(ggtext)

#Define paths
results_path <- "~/Busqueda_antigenos/All_Final_results/AllStrains_AgProtect_Final_results.tsv"
output_path <- "~/Desktop/Graficos"

Results <- read_tsv(results_path)

#Remove underscore
Results$Strain <- str_replace_all(Results$Strain,
                                  "Experimental_Antigens",
                                  "Experimental Antigens")

# Find total number of proteins per strain
Strain_number <- group_by(Results, Strain) %>% 
  summarise("Total_Number_of_proteins_per_Strain"=n()) 

# Select relevant columns and group by Strain and COG_category_u, then summarize
COG <- Results %>%
  select(ID, COG_category_u, Strain) %>%
  group_by(Strain, COG_category_u) %>%
  summarise(Total_Number_of_proteins_COG = n()) %>%
  full_join(Strain_number, by = "Strain")


Exp_Total <- Results %>%
  group_by(SubcellularLocalization, Strain) %>%
  summarise(Exp_or_NonExp_proteins_total = n())

Exp_Total <- full_join(Exp_Total, Strain_number)

Exp_Total <- Exp_Total %>%
  mutate(Percentage_Exp_Total = Exp_or_NonExp_proteins_total * 100 / Total_Number_of_proteins_per_Strain )

AgProt_Tot <- Exp_Total %>%
  filter(Strain == "Experimental Antigens") 

Mean_Percentage_Exp_Total <- Exp_Total %>%
  filter(Strain != "Experimental Antigens") %>%
  summarise(mean_percentage = round(mean(Percentage_Exp_Total), digits = 2)) 

##### Grafico Exposition ####
#filtro
Exposed <- filter(Results, AntigenicityResult ==  "ANTIGEN" & 
                   Host_Homologue_Result_All ==  "Non Host Homologue")

SL_COG <- group_by(Exposed, SubcellularLocalization, Strain, COG_category_u) %>%
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
         Strain != "Experimental Antigens") %>%
  summarise(mean_percentage = round(mean(Percentage_Exposed_or_Not), digits = 2)) 

Mean_Percentage_Exposed <- sum(Mean_Percentage_Exposed$mean_percentage)

p4C <- ggplot(Exposed, aes(
  fill= factor(SubcellularLocalization, 
               levels=c( "Cellwall", "Cytoplasmic",  "Unknown",
                         "CytoplasmicMembrane","Periplasmic", "Extracellular", "OuterMembrane")), 
  y=Percentage_Exposed_or_Not, 
  x=Strain)) + 
  geom_bar(position="stack", 
           stat="identity",
           show.legend = FALSE) +
  scale_fill_manual(values =c("sienna2", "#9e2f28","#e1ccd1", "#deb867","#4f2f4a", "olivedrab", "#25482f"),
                    labels=c("Cell Wall", "Cytoplasmic", "Unknown Localization", 
                             "Cytoplasmic Membrane","Periplasmic Space", "Extracellular Space", "Outer Membrane")
  ) +
  theme(axis.title = element_text(family = "Times New Roman", size=18, color = "black"), 
        text = element_text(family = "Times New Roman", size = 16, color = "black"),
        axis.text = element_text(family = "Times New Roman", size = 11, color = "black"),title = element_text(family = "Times New Roman"),
        legend.title=element_blank(), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) +
  coord_flip() + 
  labs(y = "Percentage of proteins",
       tag = "4A")  +
  geom_hline(yintercept = Mean_Percentage_Exposed, color = "blue")

p4C

#Agrego el total de cada grupo
SL_COG = full_join(SL_COG, COG, by=c("Strain", "COG_category_u"))

#agrego columna porcentaje
SL_COG <- SL_COG %>% 
  mutate(Exposed_Percentage = Exposed_or_Not_proteins_COG *100 / Total_Number_of_proteins_per_Strain) %>% 
  filter(substr(Strain, 1, 3) != "Exp")


means_SL <- SL_COG %>%
  group_by(COG_category_u,SubcellularLocalization) %>%
  summarize(mean_value = mean(Exposed_Percentage)) %>%
  na.omit(means_SL)

means_SL <- bind_rows(means_SL, tibble(COG_category_u = "-", 
                                       SubcellularLocalization = "Cell Wall",
                                       mean_value = 0))


p4B <- ggplot(means_SL, aes(fill=factor(SubcellularLocalization, 
                                       levels=c("Cell Wall","Unknown",
                                                "Cytoplasmic", "OuterMembrane",
                                                "CytoplasmicMembrane","Periplasmic",
                                                "Extracellular")), 
                           y=mean_value, x=COG_category_u)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values =c("sienna2", "#e1ccd1","#9e2f28","#25482f",
                              "#deb867","#4f2f4a","olivedrab"),
                    labels=c("Cell Wall",'Unknown Localization', 'Cytoplasm', 
                             'Outer Membrane', 'Citoplasmic Membrane', 
                             'Periplasmic Space', 'Extracellular Space')) +
  theme(axis.title = element_text(family = "Times New Roman", size = 16), 
        axis.text = element_text(family = "Times New Roman", size = 16, color = "black"),
        legend.title=element_blank(),
        text = element_text(family = "Times New Roman", size = 16), 
        legend.text = element_text(family = "Times New Roman", size = 16),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.8),
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.location = "plot",
        legend.key.size = unit(0.5, 'cm'),
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "Mean Percentage of Proteins", 
       x = "COG category",
       tag = "B")

p4B

Exposed <- Exposed %>% 
  mutate(
    SL_Percent = Exposed_or_Not_proteins * 100 / Total_Ag_proteins,
    Group = ifelse(Strain == "Experimental Antigens",
                   "Experimental Antigens",
                   "*Av. paragallinarum* strains")
  )

Exposed <- Exposed %>%
  mutate(SubcellularLocalization = recode(SubcellularLocalization,
                                          "Cellwall" = "Cell Wall", 
                                          "Cytoplasmic" = 'Cytoplasm',
                                          "CytoplasmicMembrane" = 'Citoplasmic\nMembrane', 
                                          "Extracellular" = 'Extracellular\nSpace',
                                          "OuterMembrane" = 'Outer\nMembrane',
                                          "Periplasmic" = 'Periplasmic\nSpace', 
                                          "Unknown" = 'Unknown\nLocalization'))
AgProt <- Exposed %>%
  filter(Strain == "Experimental Antigens")

Exposed <- Exposed %>%
  filter(Strain != "Experimental Antigens")

p4A <- ggplot(Exposed, aes(x = SubcellularLocalization, y = SL_Percent)) +
  geom_violin(width = 1.4, aes(fill = SubcellularLocalization), alpha= 0.4) +
  scale_fill_manual(values =c("#9e2f28","#deb867","olivedrab","#25482f",
                              "#4f2f4a","#e1ccd1", "sienna2")) +
  geom_point(data=AgProt, aes(x=SubcellularLocalization, y=SL_Percent, 
                              color=SubcellularLocalization), 
             size = 7, shape = 18) +
  scale_color_manual(values =c("sienna2","#9e2f28","#deb867",
                               "olivedrab","#25482f",
                               "#4f2f4a","#e1ccd1" )) +
  theme(
    text = element_text(family = "Times New Roman", size = 16, color = "black"), 
    axis.title = element_markdown(),
    legend.position = "none",
    plot.caption = element_markdown(halign = 0, hjust = 0.95),
    legend.title = element_blank(),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 13, color = "black"),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                    linewidth = 0.3)
  ) + 
  labs(
    y = "Percentage of Proteins in each subcellular<br>localization (Antigenic, Non-Homologous to Host)",
    tag = "A",
    x = "",
    caption = "٥ *Av. paragallinarum* strains<br>◆ Experimental Antigens  "
  ) + 
  ylim(0, 60)

p4A


layout <- '
AAA
BBB
'
wrap_plots(A = p4A, B = p4B, design = layout)

# Save the combined plot as a PNG file
#ggsave("4A&B.png", device = "png", path = output_path, 
 #      width =3500, height = 2000, units="px")



