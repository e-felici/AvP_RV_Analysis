library(tidyverse)
library(ggplot2)
library(ggimage)
library(patchwork)
library(png)

results_path <- "~/Busqueda_antigenos/All_Final_results"
output_path <- "~/Desktop/Graficos"
WorkDir <- "~/Busqueda_antigenos/ALL"

Results <- read_tsv(paste0(results_path, "/AllStrains_AgProtect_Final_results.tsv"))

Results$Clasif = ifelse(Results$Host_Homologue_Result_All ==  "Host Homologue",
                        "Host Homologue",
                        ifelse(Results$AntigenicityResult == "NON-ANTIGEN",
                               "Non Antigen",
                               ifelse(Results$Exposition == "Unknown" | Results$Exposition == "Inside the cell",
                                      "Non Exposed or Unknown",
                                      ifelse(Results$Conservation_Results == "Conservation Score < 0.80, less than 90% of the strains" |
                                               Results$Conservation_Results == "Conservation Score < 0.80, more than 90% of the strains" |
                                               Results$Conservation_Results == "Conservation Score > 0.80, less than 90% of the strains" |
                                               Results$Conservation_Results == "Conservation Score > 0.95, less than 90% of the strains",
                                             "Non Conserved",
                                             "Candidate protein"))))


Candidates <- Results %>% filter(
  Host_Homologue_Result_All == "Non Host Homologue" &
    AntigenicityResult == "ANTIGEN" &
    Exposition == "Exposed"
)

Candidates <- Candidates %>% filter(
  Conservation_Results == "Conservation Score > 0.95, more than 90% of the strains" |
    Conservation_Results == "Conservation Score > 0.80, more than 90% of the strains" )

#Check again that the candidate antigens are effectively present in more than 90% of the strains,
#and also share the apropiate characteristics
subfolders <- list.dirs(WorkDir, full.names = TRUE, 
                        recursive = FALSE)
num_subfolders <- length(subfolders)
num_subfolders <- round(90 * num_subfolders / 100, digits = 0)

Candidates <- Candidates %>%
  group_by(Cluster_Number) %>%
  filter(n() >= num_subfolders) %>%
  ungroup()

cluster_mapping <- c("3841" = "Peptidoglycan-associated lipoprotein Pal", 
                     "3299" = "OmpH family outer membrane protein",
                     "2960" = "OmpW family outer membrane protein",
                     "1580" = "Porin OmpA",
                     "1579" = "Outer membrane beta-barrel protein",
                     "1528" = "Peptidoglycan DD-metalloendopeptidase family protein",
                     "876" = "Outer membrane protein transport protein",
                     "864" = "TolC family protein",
                     "522" = "BamA/TamA family outer membrane protein", 
                     "267" = "TonB-dependent hemoglobin/transferrin/lactoferrin family receptor",
                     "209" = "Outer membrane protein assembly factor BamA",
                     "191" = "TonB-dependent receptor plug domain-containing protein",
                     "49" = "Filamentous hemagglutinin N-terminal domain-containing protein")

# Replace the values in Cluster_Number with letters
Candidates <- Candidates %>%
  mutate(Protein = recode(Cluster_Number, !!!cluster_mapping))

Candidates <- Candidates %>%
  mutate(VF_or_Adhe = ifelse(VirulenceFactor == "Probable Virulence Factor",
                             ifelse(AdhesinProbability2 == "Adhesin",
                                    "Virulence Factor (Probable Adhesin)",
                                    "Virulence Factor"),
                             "Non Virulence Factor")
         )

write_tsv(Candidates, paste0(output_path,"/Antigen_Candidates.tsv"))

temp_pi <- Candidates %>%
  group_by(Protein) %>%
  summarise(Isoelectric_Point = mean(Isoelectric_Point))

temp_mw <- Candidates %>%
  group_by(Protein) %>%
  summarise(MW = mean(MW))

temp_essen <- Candidates %>%
  group_by(Protein) %>%
  count(EssentialProtein) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

temp_type <- Candidates %>%
  group_by(Protein) %>%
  count(Type_of_Protein2) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

temp_vf <- Candidates %>%
  group_by(Protein) %>%
  count(VF_or_Adhe) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

temp_a <- inner_join(temp_essen, temp_mw, by = "Protein")
temp_b <- inner_join(temp_pi, temp_type, by = "Protein")
temp_c <- inner_join(temp_vf, temp_a, by = "Protein")
All <- inner_join(temp_c, temp_b, by = "Protein")

rm(list = ls(pattern = '^temp_'))

All <- select(All, !starts_with("n"))

All <- All %>% 
  mutate(MW = round(MW/1000, digits = 2), 
         Isoelectric_Point = round(Isoelectric_Point, digits = 2))

VF_mapping <- c("Virulence Factor" = "~/Desktop/Graficos/VF.png",
                "Non Virulence Factor" =  "~/Desktop/Graficos/NonVF.png")

All <- All %>%
  mutate(VF_or_Adhe2 = recode(VF_or_Adhe, !!!VF_mapping))


Type_mapping <- c("Beta Barrel" = "~/Desktop/Graficos/beta.png",
                  "Secreted" = "~/Desktop/Graficos/glob.png",
                  "Lipoprotein" = "~/Desktop/Graficos/lipo.png")
All <- All %>%
  mutate(Type_of_Protein = recode(Type_of_Protein2, !!!Type_mapping))

Essen_mapping <- c("Essential" = "~/Desktop/Graficos/esencial.png",
                  "Non Essential" = "~/Desktop/Graficos/nonesen.png")
All <- All %>%
  mutate(EssentialProtein2 = recode(EssentialProtein, !!!Essen_mapping))

All <- All  %>%
  mutate(category = rep(c(0, 1), length.out = n()))

a <- ifelse(All$category == 0, "black", "grey35")

p1 <- ggplot(All, aes(x = Protein)) + 
  geom_point(aes(y = 50, fill = Type_of_Protein2), size = 5) +
  geom_image(aes(image = Type_of_Protein), y = 50, size = 0.8) +  theme(
    axis.title = element_text(family = "Times New Roman", size = 12, color = "black"), 
    text = element_text(family = "Times New Roman", size = 12, color = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks = element_blank(),
    legend.position="right", 
    legend.justification = c(0, 1), 
    legend.key.size =unit(1, 'cm'),
    legend.text = element_text(family = "Times New Roman", size = 11, color = "black"),
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major.y = element_line(colour = "white"),
    panel.grid.major.x = element_line(colour = "grey", linetype = "dotted", linewidth = 0.3)
  ) + 
  labs(y = "", x = "") +
  scale_x_discrete(labels = scales::label_wrap(25))


p2 <- ggplot(All, aes(x = Protein, y = 50)) + 
  geom_point(aes(size = MW, fill = MW), alpha = 1, shape = 21) +
  scale_fill_gradientn(colours = c("#25482f","olivedrab","#deb867","sienna2", "#9e2f28"),
                       limits=c(10, 150), breaks=seq(10, 150, by=20),
                       labels = scales::comma) +
  guides(fill= guide_legend(), size=guide_legend()) + 
  scale_size_continuous(limits=c(10, 150), 
                        breaks=seq(10, 150, by=20), 
                        labels = scales::comma,
                        range = c(1, 10)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(family = "Times New Roman"),
        legend.title=element_blank(), 
        legend.justification = c(0, 1), 
        legend.position="right",
        legend.byrow = T,
        legend.location = "plot", 
        legend.direction = "horizontal",
        legend.key.size =unit(1, 'cm'),
        legend.text = element_text(family = "Times New Roman", size = 11, color = "black"),
        text = element_text(family = "Times New Roman", size = 12), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major.y = element_line(colour = "white"),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted", linewidth = 0.3)
  ) + 
  labs(y = "", x = "") +
  scale_x_discrete(labels = scales::label_wrap(25))

p2

p3 <- ggplot(All, aes(x = Protein, y = 50)) + 
  geom_point(aes(size = Isoelectric_Point, fill = Isoelectric_Point),
             alpha = 1, shape = 21) +
  scale_fill_gradientn(colours = c("#4f2f4a","#25482f","olivedrab","#deb867","#9e2f28"),
                       limits=c(5.5, 10.5), breaks=seq(5.5, 10.5, by=1),
                       labels = scales::comma) +
  guides(fill= guide_legend(), size=guide_legend()) + 
  scale_size_continuous(limits=c(5.5, 10.5), 
                        breaks=seq(5.5, 10.5, by=1), 
                        labels = scales::comma,
                        range = c(1, 10)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        title = element_text(family = "Times New Roman"),
        legend.title=element_blank(), 
        legend.position="right",
        legend.justification = c(0, 1), 
        legend.byrow = T,
        legend.location = "plot", 
        legend.direction = "horizontal", 
        legend.key.size =unit(1, 'cm'),
        legend.text = element_text(family = "Times New Roman", size = 11, color = "black"),
        axis.ticks = element_blank(),
        text = element_text(family = "Times New Roman", size = 12), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major.y = element_line(colour = "white"),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted", linewidth = 0.3)
  ) + 
  labs(y = "", x = "") +
  scale_x_discrete(labels = scales::label_wrap(25))

p3

p4 <- ggplot(All, aes(x = Protein)) + 
  geom_point(aes(y = 50, fill = VF_or_Adhe), size = 5) +
  geom_image(aes(image = VF_or_Adhe2), y = 50, size = 0.65) +  theme(
    axis.title = element_text(family = "Times New Roman", size = 12, color = "black"), 
    text = element_text(family = "Times New Roman", size = 12, color = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks = element_blank(), 
    legend.justification = c(0, 1), 
    legend.position="right", 
    legend.key.size =unit(1, 'cm'),
    legend.text = element_text(family = "Times New Roman", size = 11, color = "black"),
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major.y = element_line(colour = "white"),
    panel.grid.major.x = element_line(colour = "grey", linetype = "dotted", linewidth = 0.3)
      ) + 
  labs(y = "", x = "") +
  scale_x_discrete(labels = scales::label_wrap(25))

p4

p5 <- ggplot(All, aes(x = Protein)) + 
  geom_point(aes(y = 50, fill = EssentialProtein), size = 5) +
  geom_image(aes(image = EssentialProtein2), y = 50, size = 0.6) +  theme(
    axis.title = element_text(family = "Times New Roman", size = 10, color = "black"), 
    text = element_text(family = "Times New Roman", size = 10, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", angle = 90, 
                               size = 10, color = a, hjust = 1, vjust = 0.25),
    axis.title.y = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    legend.position="right",  
    legend.justification = c(0, 1), 
    legend.key.size =unit(1, 'cm'),
    legend.text = element_text(family = "Times New Roman", size = 11, color = "black"),
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major.y = element_line(colour = "white"),
    panel.grid.major.x = element_line(colour = "grey", linetype = "dotted", linewidth = 0.3)
      ) + 
  labs(y = "", x = "") +
  scale_x_discrete(labels = scales::label_wrap(25))

p5

p1 / p2 / p3 / p4 / p5

ggsave("7_2.png", device = "png", path = output_path, 
       width =3300, height = 2500, units="px")




