library(tidyverse)
library(rPref)


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
  geom_point() +
  geom_abline(slope=-1,intercept=1) +
  coord_cartesian(xlim=c(0.75,1), ylim = c(0,0.25))

All$fraction <- All$Av._paragallinarum_proteome_HOST_HOMOLOGUE / All$Experimental_Antigens_PASS



# Si x e y son valores a minimizar, entonces estás buscando los puntos en tu 
#conjunto que no están dominados por ningún otro:  es decir, la frontera de 
#Pareto de minimización.
#Formalmente:
#Una observación A “domina” a otra B si A tiene ambos valores (x, y) menores o 
#iguales, y al menos uno estrictamente menor.
#Las observaciones no dominadas son las soluciones óptimas de Pareto.
#p <- low(x) * low(y)   -->  ambas se minimizan
#frontera <- psel(df, p)

p <- high(All$Av._paragallinarum_proteome_HOST_HOMOLOGUE) * low(All$Experimental_Antigens_HOST_HOMOLOGUE) 
frontera <- psel(All, p)


