library(tidyverse)


output_path <- "~/Desktop/Graficos"
AgProtDir <- "~/Busqueda_antigenos/AgProtect"

AgProt_Homology <- read_tsv(paste0(AgProtDir,
                                   "/Homology_Analysis_results/AgProtect-Homology-all.tsv"))

HostList <- read_tsv(paste0(AgProtDir,"/HostList.tsv"), 
                     col_names = TRUE)

Ag_IDs <- read_tsv(paste0(AgProtDir,"/AllProteinIds-AgProtect.txt"), 
                   col_names = TRUE)



AgProt_Homology <- AgProt_Homology %>% mutate(Homology = ifelse(Host.Homologue=="-",
                                                                "Non-similar to host proteins",
                                                                "Similar"))

#Then, we filter only the hosts that are relevant for that antigen
#Split the Host column into separate rows
HostList2 <- HostList %>%
  separate_rows(Host, sep = ";")

AgProt_Homology <- AgProt_Homology %>%
  mutate(Host = str_remove(Host, "\\.out$"))

# Keep distinct strings in 'ID' that maximize Bit Score
AgProt_Homology <- AgProt_Homology %>%
  group_by(ID, Host) %>% 
  arrange(desc(Host.Bitscore), desc(Host.IdentityPercent)) %>% # Sort by BitScore and then idpercent (both descending)
  slice(1) %>%                                # Keep the first row in each group
  ungroup()       

AgProt_Homology <- inner_join(AgProt_Homology, HostList2, by= c("ID","Host"))

rm(HostList2, HostList)

#Add Missing IDs  
AgProt_Homology <- left_join(Ag_IDs, AgProt_Homology, by = "ID")

AgProt_Homology <- AgProt_Homology %>%
  select(ID, Host, Homology) 

AgProt<- AgProt_Homology %>%
  group_by(Host, Homology) %>%
  count(Host, Homology) 

AgProt<- AgProt %>%
  ungroup() %>%
  complete(Host, Homology, fill = list(n = 0))


AgProt$Host <- str_replace_all(AgProt$Host,
                                "_",
                                " ")

ggplot(AgProt, aes(fill=Homology, y=n, x=Host)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c( "olivedrab","#9e2f28")) +
  theme(text = element_text(family = "Times New Roman"),
        axis.text.y = element_text(face = "italic"),
        legend.title = element_blank())+ 
  labs(x = "Host", 
       y = "Proportion of experimental antigens\ncorresponding to each host")  +
  coord_flip()

ggsave("Figure_S4.jpeg", device = "jpeg", path = output_path, 
       width =190, height = 110, units="mm", dpi = 500, bg = "white")

  
