library(tidyverse)
library(ggplot2)
library(grid)
library(GGally)
library(ggtext)

#Define paths
output_path <- "~/Desktop/Graficos"

Host <- read_tsv(paste0(output_path,"/Host.tsv"))
Antigenics <- read_tsv(paste0(output_path,"/Antigenics.tsv"))
Exposed <- read_tsv(paste0(output_path,"/Exposed.tsv"))
Conserved <- read_tsv(paste0(output_path,"/Conserved.tsv"))

Host <- Host %>%
  filter(Host_Homologue_Result_All == "Non Host Homologue") %>%
  select(Strain, Host_Homologue_or_NonHH_Proteins, 
         Total_Number_of_proteins_per_Strain)

Antigenics <- Antigenics %>%
  filter(AntigenicityResult == "ANTIGEN") %>%
  select(Strain, Ag_or_NonAg_proteins)
  
Exposed <- Exposed %>%
  filter(SubcellularLocalization %in% c("Extracellular", "OuterMembrane")) %>%  
  group_by(Strain) %>%
  summarise(Exposed_proteins = sum(Exposed_or_Not_proteins))

Conserved <- Conserved %>%
  filter(Conservation_Results %in% c("Conservation Score > 0.95, more than 90% of the strains",
                                     "Conservation Score > 0.80, more than 90% of the strains")) %>%  
  group_by(Strain) %>%
  summarise(Conserved_proteins = sum(Conserved_or_Not_proteins))

All <- full_join(Host, Antigenics, by = "Strain")
All <- full_join(All, Exposed, by = "Strain")
All <- full_join(All, Conserved, by = "Strain")

rm(Antigenics,Exposed,Conserved,Host)

All <- All %>%
  rename("Total\nProteins" = Total_Number_of_proteins_per_Strain, 
         "Homology\nwith Host" = Host_Homologue_or_NonHH_Proteins, 
         "Antigenicity" = Ag_or_NonAg_proteins,
         "Subcellular\nLocalization" = Exposed_proteins,
         "Conservation" = Conserved_proteins)

All <- All %>%
  bind_rows(
    All %>%  
      filter(Strain != "Experimental Antigens") %>%    
      select(-Strain) %>%
      summarise(across(everything(), ~ round(mean(.x, na.rm = TRUE), 0))) %>%   
      mutate(Strain = "Mean")             
  )

#fictitiously assume al experimental antigens are conserved
All <- All %>%
  mutate(Conservation = if_else(Strain == "Experimental Antigens", `Subcellular
Localization`, 
Conservation))

All_Ag <- All %>%
  filter(Strain=="Experimental Antigens")%>%
  mutate(Group= "Experimental Antigens")

All_mean <- All %>%
  filter(Strain=="Mean")%>%
  mutate(Group= "Mean of Av. paragallinarum strains")

All_AvP <- All %>%
  filter(Strain!="Experimental Antigens") %>%
  filter(Strain!="Mean") %>%
  mutate(Group= "Av. paragallinarum individual strains")

# Combine the two datasets into one
All <- rbind(All_AvP, All_mean, All_Ag)

rm(All_AvP, All_mean, All_Ag)

All$Group <- as.factor(All$Group)

All$alphaLevel <- c("Av. paragallinarum individual strains" = 0.15,
                       "Experimental Antigens" = 1, 
                       "Mean of Av. paragallinarum strains" = 1)[All$Group]

grob <- grobTree(textGrob("*", x=0.89,  y=0.27, hjust=0,
                          gp=gpar(col="orangered3", fontsize=18, fontface="bold")))

All$labelss <- c("Av. paragallinarum individual strains" = "*Av. paragallinarum* individual strains",
                    "Experimental Antigens" = "Experimental Antigens", 
                    "Mean of Av. paragallinarum strains" = "Mean of *Av. paragallinarum* strains")[All$Group]


# Plot combined data using ggparcoord
ggparcoord(All,
           columns = c(3,2,4,5,6),  
           groupColumn = "labelss", 
           scale = "globalminmax",
           alphaLines = "alphaLevel",
           showPoints = TRUE,
           mapping = ggplot2::aes(linewidth = 1.5)
) +
  ggplot2::scale_linewidth_identity() +
  scale_color_manual(values = c("deepskyblue", "orangered3","navyblue"))  +
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        text = element_text(family = "Times New Roman", size = 16), 
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3),
        legend.position = "inside",
        legend.position.inside = c(0.70,0.88),
        legend.text = element_markdown(size = 16),
        legend.title = element_blank()) + 
  labs(y = "Number of proteins",
       x = "") + 
  annotation_custom(grob) +
  guides(alpha="none")

ggsave("6.png", device = "png", path = output_path, 
       width =4000, height = 2500, units="px", bg='transparent')

