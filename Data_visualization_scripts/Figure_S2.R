library(tidyverse)
library(patchwork)
library(ggthemes)


output_path <-"~/Desktop/Graficos"
Results <- read_tsv("~/Busqueda_antigenos/All_Final_results/AllStrains_AgProtect_Final_results.tsv")

#Total proteome:
Total_proteins <- Results %>%
  group_by(Cluster_Number) %>%
  n_groups()

#Filtered proteins
Results$Host_Homologue_Result_All <- str_replace_all(Results$Host_Homologue_Result_All,
                                                     "Host Homologue and Non Host Homologue, depending on the Host",
                                                     "Proteins with **sequence similarity** to at least one protein in **some hosts**")
Results$Host_Homologue_Result_All <- str_replace_all(Results$Host_Homologue_Result_All,
                                                     "Non Host Homologue",
                                                     "Proteins with **no sequence similarity** to any protein in any host")
Results$Host_Homologue_Result_All <- str_replace_all(Results$Host_Homologue_Result_All,
                                                     "Host Homologue",
                                                     "Proteins with **sequence similarity** to at least one protein in **all hosts**")
Filtered <- Results %>%
  filter(AntigenicityResult ==  "ANTIGEN" & 
           Host_Homologue_Result_All ==  "Proteins with **no sequence similarity** to any protein in any host" &
           Exposition == "Exposed")

Total_filtered_proteins <- Filtered %>%
  group_by(Cluster_Number) %>%
  n_groups()

Results <- select(Results, Cluster_Number, Conservation, Strain_count)
Filtered <- select(Filtered, Cluster_Number, Conservation, Strain_count)

#If the cluster has more than 95 members, it is that we have 2 proteins per strain that 
#fit inside that cluster. These proteins should not pass the filter
Results <- filter(Results, Strain_count < 95 )
Filtered  <- filter(Filtered, Strain_count < 95 )

cons <- c(0.75, 0.8, 0.85, 0.90, 0.95)
prev <- c(75, 80, 85, 90, 95)

#Generate all possible combinations of thresholds
thresholds <- crossing(cons, prev)

# Apply filters and summarize how many rows pass each combination
summary_results <- thresholds %>%
  pmap_dfr(\(cons, prev) {
    Results %>%
      filter(Strain_count >= prev, Conservation > cons) %>%
      distinct %>%
      summarise(
        Conservation_score = cons,
        Strain_prevalence = prev,
        n_passed = n(),
        fraction_passed = n()/Total_proteins
              )
  })

# Apply filters and summarize how many rows pass each combination
summary_filtered_results <- thresholds %>%
  pmap_dfr(\(cons, prev) {
    Filtered %>%
      filter(Strain_count >= prev, Conservation > cons) %>%
      distinct %>%
      summarise(
        Conservation_score = cons,
        Strain_prevalence = prev,
        n_passed = n(),
        fraction_passed = n()/Total_proteins
      )
  })

fill_limits <- range(
  c(summary_results$n_passed,
    summary_filtered_results$n_passed),
  na.rm = TRUE
)

scale_fill <-  scale_fill_gradient2_tableau(palette = "Temperature Diverging",
                                            limits = fill_limits,
                                            trans = scales::pseudo_log_trans(base = 10))


S4A <- ggplot(summary_results, aes(x = Conservation_score, y = Strain_prevalence, fill = n_passed)) +
  geom_tile(color = "white") +
  scale_fill +
  labs(
    title = "Total proteome",
    x = "Conservation score threshold",
    y = "Prevalence threshold",
    fill = "N° of proteins\npassed (pseudo-log10)") +
  geom_text(aes(label = round(fraction_passed, 4)), size = 1.5) +
  scale_x_continuous(limits = c(0.70, 1)) +
  scale_y_continuous(limits = c(70, 100)) +
  theme( axis.title = element_text(size = 8),
         axis.text = element_text(size = 7),
         legend.text = element_text(size = 7),
         panel.background = element_rect(fill = "white"), 
         panel.grid = element_line(colour = "grey", linetype = "dotted", 
                                         linewidth = 0.3)) + 
  labs(tag = "A")

S4A


S4B <- ggplot(summary_filtered_results, aes(x = Conservation_score, y = Strain_prevalence, fill = n_passed)) +
  geom_tile(color = "white") +
 scale_fill +
  labs(
    title = "Filtered proteome",
    x = "Conservation score threshold",
    y = "Prevalence threshold",
    fill =  "N° of proteins\npassed (pseudo-log10)") +
  geom_text(aes(label = round(fraction_passed, 4)), size = 1.5) +
  scale_x_continuous(limits = c(0.70, 1)) +
  scale_y_continuous(limits = c(70, 100)) +
  theme( axis.title = element_text(size = 8),
         axis.text = element_text(size = 7),
         legend.text = element_text(size = 7),
         panel.background = element_rect(fill = "white"), 
         panel.grid = element_line(colour = "grey", linetype = "dotted", 
                                   linewidth = 0.3))+ 
  labs(tag = "B") 

S4B



combined <- S4A + S4B & theme(text = element_text(family = "Times New Roman", size = 8, 
                                                  color = "black"),
                              legend.key.height = unit(1.7, "cm"))
combined + plot_layout(guides = "collect")


ggsave("Figure_S2.jpeg", device = "jpeg", path = output_path, 
       width =190, height = 110, units="mm", dpi = 500, bg = "white")
