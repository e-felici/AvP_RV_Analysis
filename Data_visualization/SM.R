library(tidyverse)
library(ggplot2)
library(patchwork)

output_path <- "~/Desktop/Graficos"

Candidates <- read_tsv(paste0(output_path,"/Antigen_Candidates.tsv"))

Candidates <- Candidates %>% complete(Strain, Protein,
                                      fill = list(MW = 0, 
                                                  Isoelectric_Point = 0,
                                                  AdhesinProbability2 = "Removed by the RV pipeline or non present in this strain",
                                                  Type_of_Protein2 = "Removed by the RV pipeline or non present in this strain",
                                                  VF_or_Adhe = "Removed by the RV pipeline or non present in this strain",
                                                  EssentialProtein = "Removed by the RV pipeline or non present in this strain",
                                                  COG_category_description = "Removed by the RV pipeline or non present in this strain"))


#####MW
MW <- select(Candidates, Protein, Strain, MW, Isoelectric_Point)

MW <- MW %>%
  mutate(category = rep(c(0, 1), length.out = n()))

a <- ifelse(MW$category == 0, "black", "grey35")
#bubble plot
p1 <- ggplot(MW, aes(x = Protein, y = Strain)) + 
  geom_point(aes(size = MW, fill = MW), alpha = 1, shape = 21) +
  scale_fill_gradientn(colours = c("#4f2f4a","#25482f","olivedrab","#deb867","#9e2f28"),
                       limits=c(0, 150000), breaks=seq(0, 150000, by=10000),
                       labels = scales::comma) +
  guides(fill= guide_legend(nrow = 2), size=guide_legend(nrow = 2)) + 
  scale_size_continuous(limits=c(0, 150000), breaks=seq(0, 150000, by=10000), 
                        labels = scales::comma) +
  theme(axis.title = element_text(family = "Times New Roman", color = "black"),
        axis.text.x = element_text(size = 11, angle = 90, hjust = 1, 
                                   vjust = 0.5, color = a),
        axis.text.y = element_text(size = 11, color = "black"),
        title = element_text(family = "Times New Roman"),
        legend.title=element_blank(), 
        legend.position="bottom",
        legend.byrow = T,
        legend.location = "plot", 
        legend.direction = "horizontal", 
        text = element_text(family = "Times New Roman", size = 11), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "Strain",
       x = "") +
  scale_x_discrete(labels = scales::label_wrap(25))

ggsave("Mw.png", device = "png", path = output_path, 
       width =3200, height = 4200, units="px")


pI <- select(Candidates, Protein, Strain, Isoelectric_Point)

pI <- pI %>%
  mutate(category = rep(c(0, 1), length.out = n()))

a <- ifelse(pI$category == 0, "black", "grey35")

#bubble plot
p2 <- ggplot(pI, aes(x = Protein, y = Strain)) + 
  geom_point(aes(size = Isoelectric_Point, fill = Isoelectric_Point), alpha =1, shape = 21) +
  scale_fill_gradientn(colours = c("#4f2f4a","#25482f","olivedrab","#deb867","#9e2f28"),
                       limits=c(0, 12), breaks=seq(0, 12, by=2)) +
  guides(fill= guide_legend(nrow = 1), size=guide_legend(nrow = 1)) + 
  scale_size_continuous(limits=c(0, 12), breaks=seq(0, 12, by=2)) + 
  theme(axis.title = element_text(family = "Times New Roman"),
        axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5,
                                   color = a),
        axis.text.y = element_text(size = 11, color = a),
        title = element_text(family = "Times New Roman"),
        legend.title=element_blank(), 
        legend.position="bottom",
        legend.location = "plot", 
        text = element_text(family = "Times New Roman", size = 11), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "", x = "") +
  scale_x_discrete(labels = scales::label_wrap(25))

ggsave("pi.png", device = "png", path = output_path, 
       width =3200, height = 4200, units="px")


VF <- select(Candidates, Protein, Strain, VF_or_Adhe)

VF_mapping <- c("Virulence Factor" = 10,
                  "Non Virulence Factor" =  5,
                "Virulence Factor (Probable Adhesin)"= 15,
                "Removed by the RV pipeline or non present in this strain" = 0)

VF <- VF %>%
  mutate(VF_or_Adhe2 = recode(VF_or_Adhe, !!!VF_mapping))

VF <- VF %>%
  mutate(category = rep(c(0, 1), length.out = n()))

a <- ifelse(VF$category == 0, "black", "grey35")

p3 <- ggplot(VF, aes(x = Protein, y = Strain)) + 
  geom_point(aes(size = VF_or_Adhe2, fill = VF_or_Adhe2), alpha = 1, shape = 21) +
  scale_fill_gradientn(colours = c("#9e2f28","#deb867","olivedrab","#25482f"),
                        guide = guide_legend(nrow = 1),
                       labels= c("Removed by the RV pipeline or non present in this strain",
                                 "Non Virulence Factor",
                                 "Virulence Factor",
                                 "Virulence Factor (Probable Adhesin)")
                       ) +
  scale_size_continuous(guide = guide_legend(nrow = 1),
                        labels= c("Removed by the RV pipeline or non present in this strain",
                                  "Non Virulence Factor",
                                  "Virulence Factor",
                                  "Virulence Factor (Probable Adhesin)")
                        ) +
  theme(axis.title = element_text(family = "Times New Roman"),
        axis.text.x = element_text(size = 11, angle = 90, hjust = 1, 
                                   vjust = 0.5, color = a),
        axis.text.y = element_text(size = 11, color = a),
        title = element_text(family = "Times New Roman"),
        legend.title=element_blank(),
        legend.position="bottom", 
        legend.location = "plot", 
        text = element_text(family = "Times New Roman", size = 11), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "", x = "") +
  scale_x_discrete(labels = scales::label_wrap(25))

ggsave("vfad.png", device = "png", path = output_path, 
       width =3200, height = 4200, units="px")



esen <- select(Candidates, Protein, Strain, EssentialProtein)

esen_mapping <- c("Essential" = 20,
                "Non Essential" =  10,
                "Removed by the RV pipeline or non present in this strain" = 0)

esen <- esen %>%
  mutate(EssentialProtein2 = recode(EssentialProtein, !!!esen_mapping))

esen <- esen %>%
  mutate(category = rep(c(0, 1), length.out = n()))

a <- ifelse(esen$category == 0, "black", "grey35")

p4 <- ggplot(esen, aes(x = Protein, y = Strain)) + 
  geom_point(aes(size = as.factor(EssentialProtein2), fill = as.factor(EssentialProtein2)), alpha = 1, shape = 21) +
  scale_fill_manual(values = c("#9e2f28","#deb867", "olivedrab"),
                    labels= c("Removed by the RV pipeline or non present in this strain",
                              "Non Essential",
                              "Essential")
                    ) +
  scale_size_discrete(labels= c("Removed by the RV pipeline or non present in this strain",
                                "Non Essential",
                                "Essential")
                      ) +
  theme(axis.title = element_text(family = "Times New Roman"),
        axis.text.x = element_text(size = 11, angle = 90, hjust = 1,
                                   vjust = 0.5, color = a),
        axis.text.y = element_text(size = 11, color = a),
        title = element_text(family = "Times New Roman"),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.location = "plot", 
        text = element_text(family = "Times New Roman", size = 11), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "", x = "") +
  scale_x_discrete(labels = scales::label_wrap(25)) +
  guides(fill = guide_legend(), size = guide_legend())

ggsave("essen.png", device = "png", path = output_path, 
       width =3200, height = 4200, units="px")


types <- select(Candidates, Protein, Strain, Type_of_Protein2)

types <- types %>%
  mutate(category = rep(c(0, 1), length.out = n()))

a <- ifelse(types$category == 0, "black", "grey35")

types$Type_of_Protein2 <- factor(types$Type_of_Protein2, levels = c("Removed by the RV pipeline or non present in this strain",
                                                                    "Beta Barrel",
                                                                    "Lipoprotein",
                                                                    "Secreted"))

p5 <- ggplot(types, aes(x = Protein, y = Strain)) + 
  geom_point(aes(size = Type_of_Protein2, 
                 fill = Type_of_Protein2), 
             alpha = 1, shape = 21) + 
  scale_fill_manual(
    values = c("#4f2f4a", "#9e2f28","olivedrab", "#25482f")) +
  scale_size_discrete() +
  theme(
    axis.title = element_text(family = "Times New Roman", size=11, color = "black"), 
    text = element_text(family = "Times New Roman", size = 11, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", angle = 90, 
                               size = 11, color = a, hjust = 1),
    axis.text.y = element_text(size = 11, color = a),
    legend.title = element_blank(), 
    legend.position="bottom",
    legend.location = "plot", 
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                    linewidth = 0.3)) + 
  labs(y = "", x = "") +
  scale_x_discrete(labels = scales::label_wrap(25))

ggsave("types.png", device = "png", path = output_path, 
       width =3200, height = 4200, units="px")


(p1 | p2 | p4) /
  ( p3 | p5 )

ggsave("7.png", device = "png", path = output_path, 
       width =7000, height = 7000, units="px")

