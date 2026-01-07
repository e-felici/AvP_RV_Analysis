library(tidyverse)

output_path <- "~/Desktop/Graficos"

AgProt <- read_table("~/Ablacion_ag.txt", col_names = F)
colnames(AgProt)[4] <- "AgProt"
AgProt <- AgProt %>% select(-X5)

AvP <- read_table("~/Ablacion.txt", col_names = F)
AvP <- AvP %>% group_by(X1,X2,X3) %>%
  summarise(AvP = mean(X4, na.rm = TRUE), .groups = "drop")

All <- full_join(AvP, AgProt)
rm(AgProt, AvP)

All$AgProt2 <- 1 - All$AgProt

colnames(All) <- c("Identity_percent", "e_value", "bitscore", 
                   "Av._paragallinarum_proteome_HOST_HOMOLOGUE",
                   "Experimental_Antigens_HOST_HOMOLOGUE",
                   "Experimental_Antigens_PASS")

ggplot(All, aes(x=Experimental_Antigens_PASS,
                y=Av._paragallinarum_proteome_HOST_HOMOLOGUE)) +
  geom_point(color="olivedrab") +
  geom_abline(slope=-1,intercept=1, color="#25482f") +
  theme(text = element_text(family = "Times New Roman", size = 9, color = "black"),
        axis.title.y = element_markdown(family = "Times New Roman", 
                                      color = "black")) +
  coord_cartesian(xlim=c(0.75,1), ylim = c(0,0.25)) + 
    labs(y = "Fraction of discarded *Av. paragallinarum* proteins",
         x = "Fraction of experimentally validated antigens passing the filters"
         ) 

ggsave("Figure_S1.jpeg", device = "jpeg", path = output_path, 
       width =110, height = 110, units="mm", dpi = 500, bg = "white")

