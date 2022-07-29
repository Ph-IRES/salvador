
library(purrr)
library(gridExtra)

#### Faceted Species Richness by Shallow vs. Deep Reef ####
habitatcolors <- c("#6FAFC6","#F08080")
habitat(habitatcolors) <- c("Shallow Reef", "Deep Reef")
family_chao_s <- function(
    family = "Serranidae", 
    data = data_vegan,
    data.env = data_vegan.env){
  
  # attach(data.env)
  
  pool <- 
    estimateR(x = data %>%
                dplyr::select(contains(family))) %>% 
    t() %>%
    as_tibble()
  
  
  pool %>%
    clean_names() %>%
    bind_cols(data.env) %>%
    pivot_longer(cols = s_chao1:se_ace) %>%
    mutate(se_value = case_when(str_detect(name,
                                           "se_") ~ "se",
                                TRUE ~ "value"),
           name = str_remove(name,
                             "se*_"),
           study_locations = case_when(
             site_code == "CAG" ~ "CAGAYANCILLO",
             site_code == "TUB" ~ "TRNP"
           )) %>% 
    mutate(study_locations = factor(study_locations,
                                    levels = c("TRNP",
                                               "CAGAYANCILLO"))) %>%
    pivot_wider(names_from = se_value) %>% 
    dplyr::rename(sp_richness_est = value,
           estimator = name) %>% 
    filter(estimator != "ace") %>% 
    group_by(study_locations, habitat) %>%
    dplyr::summarise(mean_chao_s = mean(sp_richness_est),
              se_chao_s = sd(sp_richness_est)/sqrt(n()),
              mean_s = mean(s_obs)) %>%
    ggplot(aes(x= study_locations,
               y= mean_chao_s,
               fill = habitat)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = mean_chao_s - se_chao_s,
                      ymax = mean_chao_s + se_chao_s), position ="dodge")+
    geom_point(aes(y = mean_s),
               color = "red3",
               position = position_dodge(width = .9)) +
    theme_classic() +
    # theme(legend.position = "none") +
    scale_fill_manual(values = habitatcolors,
                      labels =
                        c("Mesophotic Reef",
                          "Shallow Reef")) +
    labs(title = family) +
    xlab("Study Locations") +
    ylab("Mean Chao Estimate of Species Richness")
}


families <-
  colnames(data_vegan) %>%
  str_remove("_.*$") %>%
  unique()

families %>%
  purrr::map(~family_chao_s(family = .x, 
                             data = data_vegan,
                             data.env = data_vegan.env))

as.list(families) %>%
  purrr::pmap(~family_chao_s(family = ..1, 
                            data = data_vegan,
                            data.env = data_vegan.env))

plots <- families %>%
  purrr::map(~family_chao_s(family = .x))

ppl <- list(p1 = arrangeGrob(grobs=plots[1:7]))

class(ppl) <- c("arrangelist", class(plots))

ggsave("FacetedSpeciesRichness.pdf", ppl, height = 11, width = 8.5, units = "in")

#### Faceted Species Richness with Bait Type ####
family_chao_s <- function(
    family = "Serranidae", 
    data = data_vegan,
    data.env = data_vegan.env){
  
  # attach(data.env)
  
  pool <- 
    estimateR(x = data %>%
                dplyr::select(contains(family))) %>% 
    t() %>%
    as_tibble()
  
  
  pool %>%
    clean_names() %>%
    bind_cols(data.env) %>%
    pivot_longer(cols = s_chao1:se_ace) %>%
    mutate(se_value = case_when(str_detect(name,
                                           "se_") ~ "se",
                                TRUE ~ "value"),
           name = str_remove(name,
                             "se*_"),
           study_locations = case_when(
             site_code == "CAG" ~ "CAGAYANCILLO",
             site_code == "TUB" ~ "TRNP"
           )) %>% 
    mutate(study_locations = factor(study_locations,
                                    levels = c("TRNP",
                                               "CAGAYANCILLO"))) %>%
    pivot_wider(names_from = se_value) %>% 
    dplyr::rename(sp_richness_est = value,
                  estimator = name) %>% 
    filter(estimator != "ace") %>% 
    group_by(study_locations, bait_type) %>%
    dplyr::summarise(mean_chao_s = mean(sp_richness_est),
                     se_chao_s = sd(sp_richness_est)/sqrt(n()),
                     mean_s = mean(s_obs)) %>%
    ggplot(aes(x= study_locations,
               y= mean_chao_s,
               fill = bait_type)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = mean_chao_s - se_chao_s,
                      ymax = mean_chao_s + se_chao_s), position ="dodge")+
    geom_point(aes(y = mean_s),
               color = "red3",
               position = position_dodge(width = .9)) +
    theme_classic() +
    theme(legend.position = "none") +
    # scale_fill_manual(values = habitatcolors,
    #                   labels = 
    #                     c("Mesophotic Reef",
    #                       "Shallow Reef")) +
    labs(title = family) +
    xlab("Study Locations") +
    ylab("Mean Chao Estimate of Species Richness")
}


families <-
  colnames(data_vegan) %>%
  str_remove("_.*$") %>%
  unique()

families %>%
  purrr::map(~family_chao_s(family = .x, 
                            data = data_vegan,
                            data.env = data_vegan.env))

as.list(families) %>%
  purrr::pmap(~family_chao_s(family = ..1, 
                             data = data_vegan,
                             data.env = data_vegan.env))

plots_bait <- families %>%
  purrr::map(~family_chao_s(family = .x))

ppl_bait<- list(p1 = arrangeGrob(grobs=plots_bait[1:7]))

class(ppl_bait) <- c("arrangelist", class(plots_bait))

ggsave("FacetedSpeciesRichnesswithBaitType.pdf", ppl_bait, height = 11, width = 8.5, units = "in")

#### Faceted Species Richness with Groupings ####
habitatcolors <- c("#6FAFC6","#F08080")
habitat(habitatcolors) <- c("Shallow Reef", "Deep Reef")
grouping_chao_s <- function(
    groupings = "Serranidae", 
    data = data_vegan_grouping,
    data.env = data_vegan.env){
  
  # attach(data.env)
  
  pool <- 
    estimateR(x = data %>%
                dplyr::select(contains(groupings))) %>% 
    t() %>%
    as_tibble()
  
  
  pool %>%
    clean_names() %>%
    bind_cols(data.env) %>%
    pivot_longer(cols = s_chao1:se_ace) %>%
    mutate(se_value = case_when(str_detect(name,
                                           "se_") ~ "se",
                                TRUE ~ "value"),
           name = str_remove(name,
                             "se*_"),
           study_locations = case_when(
             site_code == "CAG" ~ "CAGAYANCILLO",
             site_code == "TUB" ~ "TRNP"
           )) %>% 
    mutate(study_locations = factor(study_locations,
                                    levels = c("TRNP",
                                               "CAGAYANCILLO"))) %>%
    mutate(habitat = factor(habitat,
                            levels = c("Shallow Reef",
                                       "Deep Reef"))) %>%
    pivot_wider(names_from = se_value) %>% 
    dplyr::rename(sp_richness_est = value,
                  estimator = name) %>% 
    filter(estimator != "ace") %>% 
    group_by(study_locations, habitat) %>%
    dplyr::summarise(mean_chao_s = mean(sp_richness_est),
                     se_chao_s = sd(sp_richness_est)/sqrt(n()),
                     mean_s = mean(s_obs)) %>%
    ggplot(aes(x= study_locations,
               y= mean_chao_s,
               fill = habitat)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = mean_chao_s - se_chao_s,
                      ymax = mean_chao_s + se_chao_s), position ="dodge")+
    geom_point(aes(y = mean_s),
               color = "red3",
               position = position_dodge(width = .9)) +
    theme_classic() +
    # theme(legend.position = "none") +
    scale_fill_manual(values = habitatcolors,
                      labels =
                        c("Shallow Reef",
                          "Mesophotic Reef")) +
    labs(title = groupings, fill = "Habitat") +
    xlab("Study Locations") +
    ylab("Mean Chao Estimate of Species Richness")
}


groupings1 <-
  unique(data_all_removed_sp$groupings)

groupings1 %>%
  purrr::map(~grouping_chao_s(groupings = .x, 
                            data = data_vegan_grouping,
                            data.env = data_vegan.env))

# as.list(families) %>%
#   purrr::pmap(~grouping_chao_s(groupings = ..1, 
#                              data = data_vegan_grouping,
#                              data.env = data_vegan.env))

plots_grouping <- groupings1 %>%
  purrr::map(~grouping_chao_s(groupings = .x, 
                              data = data_vegan_grouping,
                              data.env = data_vegan.env))

ppl_grouping <- list(p1 = arrangeGrob(grobs=plots_grouping[1:6]))

class(ppl_grouping) <- c("arrangelist", class(plots_grouping))

ggsave("FacetedSpeciesRichnesswithGroupings.pdf", ppl_grouping, height = 11, width = 8.5, units = "in")
