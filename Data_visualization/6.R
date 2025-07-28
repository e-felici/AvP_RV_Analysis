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
  filter(Conservation_Results == "**Conserved sequence** (CS > 0.80) +<br>**prevalent** in > 90 % of the strains") %>%  
  group_by(Strain) %>%
  summarise(Conserved_proteins = sum(Conserved_or_Not_proteins))

All <- full_join(Host, Antigenics, by = "Strain")
All <- full_join(All, Exposed, by = "Strain")
All <- full_join(All, Conserved, by = "Strain")

rm(Antigenics,Exposed,Conserved,Host)

All <- All %>%
  rename("Total\nProteins" = Total_Number_of_proteins_per_Strain, 
         "No Homology\nwith Host" = Host_Homologue_or_NonHH_Proteins, 
         "Antigenic" = Ag_or_NonAg_proteins,
         "Exposed on\nCell Surface" = Exposed_proteins,
         "Conserved and\nPrevalent" = Conserved_proteins)

All <- All %>%
  bind_rows(
    All %>%  
      filter(Strain != "Experimental Antigens") %>%    
      select(-Strain) %>%
      summarise(across(everything(), ~ round(mean(.x, na.rm = TRUE), 0))) %>%   
      mutate(Strain = "Mean")             
  )

All <- All %>% 
  mutate(Percentage_retained = All$`Conserved and
Prevalent` * 100 / All$`No Homology
with Host`)

All <- All %>%
  mutate(Percentage_retained = if_else(Strain == "Experimental Antigens", `Exposed on
Cell Surface` * 100 / `No Homology
with Host`, 
Percentage_retained))

All_Ag <- All %>%
  filter(Strain=="Experimental Antigens")%>%
  mutate(Group= "Experimental Antigens")

All_mean <- All %>%
  filter(Strain=="Mean")%>%
  mutate(Group= "Mean of *Av. paragallinarum* strains")

All_AvP <- All %>%
  filter(Strain!="Experimental Antigens") %>%
  filter(Strain!="Mean") %>%
  mutate(Group= "*Av. paragallinarum* individual strains")

All_AvP$Group <- as.factor(All_AvP$Group)
All_mean$Group <- as.factor(All_mean$Group)
All_Ag$Group <- as.factor(All_Ag$Group)

All_AvP$alphaLevel <- 0.15
All_Ag$alphaLevel <- 1
All_mean$alphaLevel <-1

#NA placeholder
All_Ag <- All_Ag %>% mutate(`Conserved and
Prevalent` = `Exposed on
Cell Surface`)

All <- full_join(All_AvP, All_Ag)
All <- full_join(All, All_mean)

rm(All_AvP, All_Ag, All_mean)

grob <- grobTree(textGrob("*", x=0.70,  y=0.25, hjust=0,
                          gp=gpar(col="#0073C2FF", fontsize=18, fontface="bold")))

# Plot combined data using ggparcoord
ggparcoord(All,
           columns = c(3,2,4,5,6),  
           groupColumn = "Group", 
           scale = "globalminmax",
           alphaLines = "alphaLevel",
           showPoints = TRUE,
           mapping = ggplot2::aes(linewidth = 1.5)
) +
  ggplot2::scale_linewidth_identity() +
  scale_color_manual(values = c("olivedrab","#0073C2FF", "orangered3"))  +
  theme(axis.text.x=element_text(colour="black", face = "bold"),
        axis.text.y=element_text(colour="black"),
        axis.title.y = element_text(face = "bold"),
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
  ylim(0,2500) +
  annotation_custom(grob) +
  guides(alpha="none")


ggsave("6.png", device = "png", path = output_path, 
       width =2500, height = 2500, units="px", bg='transparent')

###Note: manually fixed the experimental antigens column