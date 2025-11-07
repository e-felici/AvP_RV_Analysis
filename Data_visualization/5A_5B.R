library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(patchwork)

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
Strain_number = group_by(Results, Strain) %>% 
  summarise("Total_Number_of_proteins_per_Strain"=n()) 

COG <- Results %>%
  filter(Host_Homologue_Result_All == "Proteins with **no sequence similarity** to any protein in any host",
         AntigenicityResult == "ANTIGEN", 
         Exposition == "Exposed") %>%
  select(ID, COG_category_description, Strain) %>%
  group_by(Strain, COG_category_description) %>%
  summarise(Total_Number_of_proteins_COG = n()) %>%
  full_join(Strain_number, by = "Strain")

Cons_Total <- Results %>%
  group_by(Conservation_Results, Strain) %>%
  summarise(Cons_or_NonCons_proteins_total = n())

Cons_Total <- full_join(Cons_Total, Strain_number)

Cons_Total <- Cons_Total %>%
  mutate(Percentage_Cons_Total = Cons_or_NonCons_proteins_total * 100 / Total_Number_of_proteins_per_Strain )

Mean_Percentage_Cons_Total <- Cons_Total %>%
  filter(Strain != "Experimental antigens") %>%
  summarise(mean_percentage = round(mean(Percentage_Cons_Total), digits = 2)) 


#### Grafico Conservation ####
#filtro

Conserved <- filter(Results, AntigenicityResult ==  "ANTIGEN" & 
                      Host_Homologue_Result_All ==  "Proteins with **no sequence similarity** to any protein in any host" &
                      Exposition == "Exposed")

Conserved <- Conserved %>%
  mutate(Conservation_Results = recode(Conservation_Results,
                                       "Conservation Score < 0.80, less than 90% of the strains" = "**Sequence not conserved** +<br>**Not prevalent**",
                                       "Conservation Score > 0.80, less than 90% of the strains" = "Conserved sequence (CS > 0.80)<br>but **not prevalent**",
                                       "Conservation Score > 0.95, less than 90% of the strains" = "Conserved sequence (CS > 0.80)<br>but **not prevalent**",
                                       "Conservation Score < 0.80, more than 90% of the strains" = "Prevalent in > 90 % of the strains<br>but **sequence not conserved**",
                                       "Conservation Score > 0.80, more than 90% of the strains" = "**Conserved sequence** (CS > 0.80) +<br>**prevalent** in > 90 % of the strains",
                                       "Conservation Score > 0.95, more than 90% of the strains" = "**Conserved sequence** (CS > 0.80) +<br>**prevalent** in > 90 % of the strains"
  ))

Cons_COG <- group_by(Conserved, Conservation_Results, Strain, COG_category_description) %>%
  summarise(Conserved_or_Not_proteins_COG=n())

Conserv1 <- group_by(Conserved, Conservation_Results, Strain) %>% 
  summarise(Conserved_or_Not_proteins=n()) 
Conserv2= group_by(Conserved, Strain) %>% summarise("Total_Exposed_proteins"=n()) 
Conserved= full_join(Conserv1, Conserv2, by="Strain")
rm(Conserv1, Conserv2)

#agrego columna porcentaje
Conserved <- Conserved %>%
  mutate(Percentage_Conserved_or_Not = Conserved_or_Not_proteins *100 / Total_Exposed_proteins)

write_tsv(Conserved, paste0(output_path,"/Conserved.tsv"))

Conserved <- Conserved%>% filter(substr(Strain, 1, 3) != "Exp")

Mean_Percentage_Conserved <- Conserved%>%
  filter(Conservation_Results == "**Conserved sequence** (CS > 0.80) +<br>**prevalent** in > 90 % of the strains") %>%
  summarise(mean_percentage = round(mean(Percentage_Conserved_or_Not), digits = 2)) 

Mean_Percentage_Conserved <- sum(Mean_Percentage_Conserved$mean_percentage)

#Agrego el total de cada grupo
Cons_COG = full_join(Cons_COG, Conserved, by=c("Strain", "Conservation_Results"))

Cons_COG <- Cons_COG %>% select(Strain, Conservation_Results, COG_category_description,
                                Conserved_or_Not_proteins_COG, Total_Exposed_proteins)

#agrego columna porcentaje
Cons_COG <- Cons_COG %>% 
  mutate(Conserved_Percentage = Conserved_or_Not_proteins_COG * 100 / Total_Exposed_proteins)

Cons_COG <- Cons_COG%>% filter(substr(Strain, 1, 3) != "Exp")


means_Cons <- Cons_COG %>%
  group_by(COG_category_description, Conservation_Results) %>%
  summarize(mean_value = mean(Conserved_Percentage))

means_Cons$COG_category_description <- str_replace_all(means_Cons$COG_category_description,
                                                     "_",
                                                     " ")

means_Cons$COG_category_description <- str_replace_all(means_Cons$COG_category_description,
                                                     "-",
                                                     "Proteins not assigned to any COG")

means_Cons$COG_category_description <- str_replace_all(means_Cons$COG_category_description,
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

means_Cons= na.omit(means_Cons)

means_Cons$Conservation_Results <- factor(means_Cons$Conservation_Results, 
                                          levels = c("**Sequence not conserved** +<br>**Not prevalent**",
                                                     "Conserved sequence (CS > 0.80)<br>but **not prevalent**",
                                                     "Prevalent in > 90 % of the strains<br>but **sequence not conserved**",
                                                     "**Conserved sequence** (CS > 0.80) +<br>**prevalent** in > 90 % of the strains"))

means_Cons$COG_category_description <- factor(means_Cons$COG_category_description,
                                              levels = COG_category_description_order)

p5B <- ggplot(means_Cons, aes( fill= Conservation_Results, 
                               y = mean_value, x = fct_rev(COG_category_description))) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values=c("#4f2f4a","#deb867","#9e2f28","olivedrab")) +
  #scale_fill_viridis(discrete=T, begin = 0.3, end=0.65)  +
  theme(axis.title = element_markdown(face="bold"), 
        axis.text = element_text(family = "Times New Roman", size = 10, color = "black"),
        title = element_text(family = "Times New Roman"), 
        legend.position= "bottom",
        legend.location = "plot",
        legend.byrow = T,
        legend.title = element_blank(),
        legend.text = element_markdown(),
        axis.title.x.bottom = element_markdown(face="bold"), 
        text = element_text(family = "Times New Roman", size = 14), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "Mean % of non-similar to host<br>sequences, surface-exposed,<br>antigenic proteins in *Av.<br>paragallinarum* strains", 
       x = "COG category",
       tag = "B") +
  coord_flip()

p5B

p5A <- ggplot(Conserved, aes(x = Conservation_Results, y = Percentage_Conserved_or_Not)) +
  geom_violin(width = 1.4, aes(fill = Conservation_Results)) +
  scale_fill_manual(values =c("olivedrab","#4f2f4a","#deb867","#9e2f28")) +
  theme(
    text = element_text(family = "Times New Roman", size = 14, color = "black"), 
    axis.title = element_markdown(face="bold"),
    legend.position = "none",
    plot.caption = element_markdown(halign = 0, hjust = 0.95),
    legend.title = element_blank(),
    axis.title.x.bottom = element_markdown(face="bold"), 
    axis.text.y = element_markdown(family = "Times New Roman", size = 12, 
                                   color = "black"),
    axis.text.x = element_markdown(family = "Times New Roman", size = 14, 
                                   color = "black", face = "bold"),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                    linewidth = 0.3)
  ) + 
  labs(
    y = "Conservation across *Av.paragallinarum*<br>proteins (non-similar to host sequences,<br>surface-exposed and antigenic)",
    tag = "A",
    x = ""
  ) + 
  coord_flip()

p5A

layout <- '
AAABBBB
AAABBBB
AAABBBB
AAABBBB
DDDD###
'
wrap_plots(A = p5A, B = p5B, D = guide_area(), design = layout) +
  plot_layout(guides = 'collect')

ggsave("5A&B.png", device = "png", path = output_path, 
       width =3500, height = 1800, units="px")




