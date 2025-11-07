library(tidyverse)

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

ggplot(summary_results, aes(x = Conservation_score, y = Strain_prevalence, fill = n_passed)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "olivedrab", high = "firebrick4") +
  labs(
    title = "Number of proteins passing each threshold combination",
    subtitle = "Total proteome reduction",
    x = "Prevalence threshold",
    y = "Conservation threshold",
    fill = "n passed") +
  geom_text(aes(label = round(fraction_passed, 4))) +
  scale_x_continuous(limits = c(0.70, 1)) +
  scale_y_continuous(limits = c(70, 100)) +
  theme( panel.background = element_rect(fill = "white"), 
         panel.grid = element_line(colour = "grey", linetype = "dotted", 
                                         linewidth = 0.3))

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

ggplot(summary_filtered_results, aes(x = Conservation_score, y = Strain_prevalence, fill = n_passed)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "olivedrab", high = "firebrick4") +
  labs(
    title = "Number of proteins passing each threshold combination",
    subtitle = "Filtered proteome (antigenic, exposed in the surface non-similar to chicken proteins) reduction",
    x = "Prevalence threshold",
    y = "Conservation threshold",
    fill = "n passed") +
  geom_text(aes(label = round(fraction_passed, 4))) +
  scale_x_continuous(limits = c(0.70, 1)) +
  scale_y_continuous(limits = c(70, 100)) +
  theme( panel.background = element_rect(fill = "white"), 
         panel.grid = element_line(colour = "grey", linetype = "dotted", 
                                   linewidth = 0.3))

