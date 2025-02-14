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



#### Grafico Conservation ####
#filtro

Conserved <- filter(Results, AntigenicityResult ==  "ANTIGEN" & 
                      Host_Homologue_Result_All ==  "Non Host Homologue" &
                      Exposition == "Exposed")

Cons_COG <- group_by(Conserved, Conservation_Results, Strain, COG_category_u) %>%
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


#calculo porcentaje de homologas
Mean_Percentage_Conserved <- Conserved%>%
  filter(Conservation_Results == paste0("Conservation Score > 0.95, more than 90% of the strains")
         | Conservation_Results == paste0("Conservation Score > 0.80, more than 90% of the strains") ,
         Strain != "Experimental Antigens") %>%
  summarise(mean_percentage = round(mean(Percentage_Conserved_or_Not), digits = 2)) 

Mean_Percentage_Exposed <- sum(Mean_Percentage_Conserved$mean_percentage)


##grafico
p5A <- ggplot(Conserved, aes(
  fill= factor(Conservation_Results, 
               levels=c("Conservation Score < 0.80, less than 90% of the strains",
                        "Conservation Score > 0.80, less than 90% of the strains",
                        "Conservation Score > 0.95, less than 90% of the strains",
                        "Conservation Score < 0.80, more than 90% of the strains",
                        "Conservation Score > 0.80, more than 90% of the strains",
                        "Conservation Score > 0.95, more than 90% of the strains")),
  y=Percentage_Conserved_or_Not, 
  x=Strain)) +  
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(name = " ", values=c("#9e2f28","#e1ccd1","#4f2f4a","#deb867","olivedrab","#25482f")) +
  theme(axis.title = element_text(family = "Times New Roman", size=16, color = "black"), 
        title = element_text(family = "Times New Roman", size=16, color = "black"), 
        text = element_text(family = "Times New Roman", size = 16, color = "black"),
        axis.text = element_text(family = "Times New Roman", size = 11, color = "black"),
        legend.title=element_blank(), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) +
  coord_flip() + 
  labs(y = "Percentage of proteins", 
       tag = "5A") +
  geom_hline(yintercept = Mean_Percentage_Exposed, color = "blue")

p5A


#Agrego el total de cada grupo
Cons_COG = full_join(Cons_COG, COG, by=c("Strain", "COG_category_u"))

#agrego columna porcentaje
Cons_COG <- Cons_COG %>% 
  mutate(Conserved_Percentage = Conserved_or_Not_proteins_COG * 100 / Total_Number_of_proteins_per_Strain)

Cons_COG <- Cons_COG%>% filter(substr(Strain, 1, 3) != "Exp")


means_Cons <- Cons_COG %>%
  group_by(COG_category_u, Conservation_Results) %>%
  summarize(mean_value = mean(Conserved_Percentage))

means_Cons= na.omit(means_Cons)

p5B <- ggplot(means_Cons, aes( fill= factor(Conservation_Results, 
                                            levels=c("Conservation Score < 0.80, less than 90% of the strains",
                                                     "Conservation Score > 0.80, less than 90% of the strains",
                                                     "Conservation Score > 0.95, less than 90% of the strains",
                                                     "Conservation Score < 0.80, more than 90% of the strains",
                                                     "Conservation Score > 0.80, more than 90% of the strains",
                                                     "Conservation Score > 0.95, more than 90% of the strains")), 
                               y = mean_value, x = COG_category_u)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values=c("#9e2f28","#e1ccd1","#4f2f4a","#deb867","olivedrab","#25482f")) +
  #scale_fill_viridis(discrete=T, begin = 0.3, end=0.65)  +
  theme(axis.title = element_text(family = "Times New Roman", size=18, color = "black"), 
        text = element_text(family = "Times New Roman", size = 16, color = "black"),
        axis.text = element_text(family = "Times New Roman", size = 16, color = "black"),
        legend.title=element_blank(), 
        panel.background =  element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", 
                                        linewidth = 0.3)) + 
  labs(y = "Mean for all strains", 
       x = "COG category",
       tag = "5B")

p5B

ggarrange(p5A, p5B, 
          ncol = 2, nrow = 1,
          legend = "bottom",
          common.legend = T)

ggsave("5A&B.png", device = "png", path = output_path, 
       width =4200, height = 3500, units="px")




