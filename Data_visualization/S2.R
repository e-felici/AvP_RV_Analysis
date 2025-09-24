library(tidyverse)
library(ggtext)

results_path <- "~/Busqueda_antigenos/All_Final_results"
output_path <- "~/Desktop/Graficos"

Results <- read_tsv(paste0(results_path, "/AllStrains_AgProtect_Final_results.tsv"))

Results <- select(Results, ID, Strain, Cluster_Number, Description)

fil_haema_hxua <- filter(Results, Cluster_Number == 49)

hmtp210 <- filter(Results, Cluster_Number == 22 | 
                    Cluster_Number == 19 | 
                    Cluster_Number == 20)

hagA <- filter(Results, Cluster_Number == 1580)

fil_haema <- Results   %>%
  filter(Cluster_Number != 49) %>%
  filter(if_any(everything(), ~ str_detect(as.character(.), regex("filamentous_hemagglutinin", ignore_case = TRUE))))

haema_rep <- Results   %>%
  filter(if_any(everything(), ~ str_detect(as.character(.), regex("hemagglutinin_repeat", ignore_case = TRUE))))

hagA$Haemagglut <- "HagA (Porin OmpA)"
hmtp210$Haemagglut <- "HMTp210 (YadA-like Family Protein)"
haema_rep$Haemagglut <- "Haemagglutinin Repeat Containing Protein"
fil_haema$Haemagglut <- "Another Filamentous Haemagglutinin"
fil_haema_hxua$Haemagglut <- "Filamentous Haemagglutinin (Homologue to HxuA from *H. influenzae*)"

haema <- full_join(hmtp210, hagA)
haema <- full_join(haema, fil_haema)  
haema <- full_join(haema, haema_rep)
haema <- full_join(haema, fil_haema_hxua)
rm(haema_rep, fil_haema, fil_haema_hxua, hagA, hmtp210)

haema <- haema %>%
  filter(Strain != "Experimental_Antigens")

#Preliminary plot
ggplot(haema, aes(as.factor(Cluster_Number), Strain, 
                  fill= Haemagglut)) +
  geom_tile() + 
  scale_fill_manual(values = c("olivedrab", "sienna2","#9e2f28",
                               "#25482f","#deb867","white")) +
  theme(panel.background = element_rect(fill = "white"),
        legend.text = element_markdown(),
        text = element_text(family = "Times New Roman", size = 14)) + 
  labs(y = "Strain", 
       x = "Clusters Containing Haemagglutinin Proteins",
       tag = "S2") 

pattern <- unique(haema$Cluster_Number)

Results <- Results %>%
  filter(Cluster_Number %in% pattern)

Results <- anti_join(
  Results,
  haema,
  by = c("ID", "Cluster_Number", "Strain", "Description")
)

Results$Haemagglut <- "Protein Not Designated as Haemagglutinin"

haema2 <- filter(haema, Cluster_Number == "355" | Cluster_Number == "192" | 
                   Cluster_Number == "324" )

haema <- anti_join(haema, haema2)

haema <- full_join(haema, Results)

haema <- full_join(haema, haema2)

rm(Results, haema2)

haema <- haema %>% 
  select(ID, Cluster_Number, Strain, Haemagglut)

haema$Cluster_Number <- as.factor(haema$Cluster_Number)

haema <- haema %>%
  mutate(Haemagglut = factor(Haemagglut,
                             levels = c("Filamentous Haemagglutinin (Homologue to HxuA from *H. influenzae*)",
                                        "Another Filamentous Haemagglutinin",
                                        "Haemagglutinin Repeat Containing Protein",
                                        "HMTp210 (YadA-like Family Protein)",
                                        "HagA (Porin OmpA)",
                                        "Protein Not Designated as Haemagglutinin"
                                        )))

ggplot(haema, aes(Cluster_Number, Strain, 
                  fill= Haemagglut)) +
  geom_tile() + 
  scale_fill_manual(values = c("olivedrab", "sienna2","#9e2f28",
                               "#deb867","#25482f","#4f2f4a")) +
  theme(panel.background = element_rect(fill = "white"),
        legend.text = element_markdown(),
        legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(family = "Times New Roman", size = 14)) + 
  labs(y = "Strain", 
       x = "Clusters Containing Haemagglutinin Proteins",
       tag = "S2") 

ggsave("S2.png", device = "png", path = output_path, 
       width =4000, height = 3500, units="px")
