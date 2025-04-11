# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(ggpubr) 

#Define paths
results_path <- "~/Busqueda_antigenos/All_Final_results/AllStrains_AgProtect_Final_results.tsv"
output_path <- "~/Desktop/Graficos"

Results <- read_tsv(results_path)

#Remove underscore
Results$Strain <- str_replace_all(Results$Strain,
                                  "Experimental_Antigens",
                                  "Experimental Antigens")

# Find total number of proteins per strain
Strain_number = group_by(Results, Strain) %>% 
  summarise("Total_Number_of_proteins_per_Strain"=n()) 

# Select relevant columns and group by Strain and COG_category_u, then summarize
COG <- Results %>%
  select(ID, COG_category_u, Strain) %>%
  group_by(Strain, COG_category_u) %>%
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
  
Mean_Percentage_Ag_Total <- Ag_Total %>%
  filter(AntigenicityResult == "ANTIGEN", 
         Strain != "Experimental Antigens") %>%
  summarise(mean_percentage = round(mean(Percentage_Ag_Total), digits = 2)) %>%
  pull(mean_percentage)

# Filter non-host homologue proteins
Antigenics <- filter(Results, Host_Homologue_Result_All == "Non Host Homologue")

# We assume that proteins that could not be processed by VaxiJen (because they 
#had rare amino acids) are non antigenic
Antigenics$AntigenicityResult <- ifelse(Antigenics$AntigenicityResult == "-",
                                        "NON-ANTIGEN",
                                        Antigenics$AntigenicityResult)

Ag_COG <- Antigenics %>%
  group_by(AntigenicityResult, Strain, COG_category_u) %>%
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

# Calculate the mean percentage of non-antigenic proteins
Mean_Percentage_Antigenic <- Antigenics %>%
  filter(AntigenicityResult == "NON-ANTIGEN", 
         Strain != "Experimental Antigens") %>%
  summarise(mean_percentage = round(mean( Percentage_Antigenic_or_Not), digits = 2)) %>%
  pull(mean_percentage)

# Create the first plot (Grafico 2)
p3C <- ggplot(Antigenics, aes(fill = AntigenicityResult, y =  Percentage_Antigenic_or_Not, 
                             x = Strain)) + 
  geom_bar(position = "stack", stat = "identity", show.legend = FALSE) + 
  scale_fill_manual(values = c("olivedrab", "#9e2f28"), labels = c('Antigenic', 'Non Antigenic')) +
  theme(
    axis.title = element_text(family = "Times New Roman", size = 16, color = "black"), 
    title = element_text(family = "Times New Roman", size = 16, color = "black"), 
    text = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 11, color = "black"),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_line(colour = "grey", linetype = "dotted", linewidth = 0.3)
  ) +
  coord_flip() + 
  labs(y = "Percentage of proteins",
       tag = "3A") +
  geom_hline(yintercept = Mean_Percentage_Antigenic, color = "blue")

p3C

# Add the total number of proteins per strain and COG category
Ag_COG <- full_join(Ag_COG, COG, by = c("Strain", "COG_category_u"))

# Add percentage column
Ag_COG <- Ag_COG %>%
  mutate(Ag_Percentage = Ag_or_NonAg_proteins_COG * 100 / Total_Number_of_proteins_per_Strain) %>%
  na.omit() %>%
  filter(!str_starts(Strain, "Exp"))

# Calculate the mean percentage of antigenic proteins per COG category
means_ag <- Ag_COG %>%
  group_by(COG_category_u, AntigenicityResult) %>%
  summarize(mean_value = mean(Ag_Percentage))

# Create the second plot (Grafico 2)
p3B <- ggplot(means_ag, aes(fill = AntigenicityResult, y = mean_value, x = COG_category_u)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = c("olivedrab", "#9e2f28"), labels = c('Antigenic', 'Non Antigenic')) +
  theme(axis.title = element_markdown(), 
    axis.text = element_text(family = "Times New Roman", size = 16, color = "black"),
    title = element_text(family = "Times New Roman"),
    legend.title = element_blank(),  
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.9),
    legend.box = "vertical",
    legend.direction = "vertical",
    text = element_text(family = "Times New Roman", size = 16), 
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_line(colour = "grey", linetype = "dotted", linewidth = 0.3)
  ) + 
  labs(y = "Mean Percentage of Proteins", 
       x = "COG category",
       tag = "B")

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

p3A <- ggplot(Antigenics, aes(x = Group, y = Ag_Percent)) +
  geom_violin(width = 1.4, aes(fill = Group)) +
  scale_fill_manual(values = "olivedrab") +
  geom_point(aes(y = AgProt, color = "Experimental Antigens"), size = 5, shape = 18) +  
  scale_color_manual(values = c("Experimental Antigens" = "#0073C2FF")) +  
  theme(
    text = element_text(family = "Times New Roman", size = 16, color = "black"), 
    axis.text = element_text(family = "Times New Roman", size = 16, color = "black"),
    legend.position = "inside",
    legend.position.inside = c(0.5, 0),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.text = element_markdown(),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                    linewidth = 0.3)
  ) + 
  labs(
    y = "Percentage of Antigenic Proteins (Non-Homologous to Host)",
    tag = "A",
    x = ""
  ) + 
  ylim(0, 83)

p3A


layout <- '
ABBBB
'
wrap_plots(A = p3A, B = p3B, design = layout)

# Save the combined plot as a PNG file
ggsave("3A&B.png", device = "png", path = output_path, 
       width =3500, height = 2000, units="px")
