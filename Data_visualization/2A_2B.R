# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(ggpubr) 
library(patchwork)
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
Strain_number = group_by(Results, Strain) %>% 
  summarise("Total_Number_of_proteins_per_Strain"=n()) 

# Select relevant columns and group by Strain and COG_category_u, then summarize
COG <- Results %>%
  select(ID, COG_category_u, Strain) %>%
  group_by(Strain, COG_category_u) %>%
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
  filter(Host_Homologue_Result_All == "Non Host Homologue", 
         Strain != "Experimental Antigens") %>%
  summarise(mean_percentage = round(mean(Percentage_Host), digits = 2)) %>%
  pull(mean_percentage)

# Create the first plot (Figure 2A)
p2C <- ggplot(Host, aes(fill=Host_Homologue_Result_All, y=Percentage_Host, x=Strain)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values =c("#9e2f28","#e1ccd1","olivedrab")) +
  theme(axis.title = element_text(family = "Times New Roman", size=16, color = "black"), 
        title = element_text(family = "Times New Roman", size=16, color = "black"), 
        text = element_text(family = "Times New Roman", size = 16, color = "black"),
        axis.text = element_text(family = "Times New Roman", size = 11, color = "black"),
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3),
        legend.title = element_blank(),
        legend.position = "bottom"
  ) +
  coord_flip() + 
  labs(y = "Percentage of proteins",
       tag = "2C") +
  geom_hline(yintercept = Mean_Percentage_Host, color = "blue") 

p2C


#######Figure 2B
#Create Host and COG tibble
Host_COG <- Results %>%
  group_by(Host_Homologue_Result_All, Strain, COG_category_u) %>%
  summarise(Total_Number_of_proteins_Host_COG = n()) %>%
  full_join(COG, by = c("Strain", "COG_category_u")) %>%
  mutate(Percentage_Host_COG = Total_Number_of_proteins_Host_COG * 100 / Total_Number_of_proteins_per_Strain) %>%
  filter(!str_starts(Strain, "Exp"))


# Calculate the mean percentage of host homologue proteins per COG category
means_Host_COG <- Host_COG %>%
  group_by(COG_category_u,Host_Homologue_Result_All) %>%
  summarize(mean_value = mean(Percentage_Host_COG))

# Create the second plot (Figure 2B)
p2B <- ggplot(means_Host_COG, aes(fill=Host_Homologue_Result_All, 
                                  y=mean_value, x=COG_category_u)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values =c("#9e2f28","olivedrab")) +
  theme(axis.title = element_markdown(), 
        axis.text = element_text(family = "Times New Roman", size = 16, color = "black"),
        title = element_text(family = "Times New Roman"), 
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.9),
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.title = element_blank(),
        text = element_text(family = "Times New Roman", size = 16), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "Mean Percentage of Proteins", 
       x = "COG category",
       tag = "B")

p2B


Host <- Host %>% 
  filter(Host_Homologue_Result_All == "Non Host Homologue") %>%
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

p2A <- ggplot(Host, aes(x = Group, y = HH_Percent)) +
  geom_violin(width = 1.4, aes(fill = Group)) +
  scale_fill_manual(values = "olivedrab") +
  geom_point(aes(y = AgProt, color = "Experimental Antigens"), size = 5, shape = 18) +  
  scale_color_manual(values = c("Experimental Antigens" = "#0073C2FF")) +  
  theme(
    text = element_text(family = "Times New Roman", size = 16, color = "black"), 
    axis.text = element_text(family = "Times New Roman", size = 16, color = "black"),
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.05),
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
    y = "Percentage of Proteins Non-Homologous to Host",
    tag = "A",
    x = ""
  ) + 
  ylim(69, 83)

p2A


layout <- '
ABBBB
'
wrap_plots(A = p2A, B = p2B, design = layout)

# Save the combined plot as a PNG file
ggsave("2A&B.png", device = "png", path = output_path, 
       width =3500, height = 2000, units="px")
