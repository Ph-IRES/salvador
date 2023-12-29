#### INITIALIZATION ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#install.packages("remotes")
#remotes::install_github("gavinsimpson/ggvegan")
library(ggvegan)
library(tidyverse)
library(janitor)
library(readxl)
library(sjPlot)
# install.packages("maps")
# install.packages("viridis")
# require(maps)
# require(viridis)
# theme_set(
#   theme_void()
# )
# install.packages("vegan")
library(vegan)
library(devtools)
# devtools::install_github("jfq3/ggordiplots", dependencies = TRUE)
library(ggordiplots)
# devtools::install_github("vqv/ggbiplot")
library(ggbiplot)
# devtools::install_github("gavinsimpson/ggvegan")
library(ggvegan)
library(ggpubr)

library(ggrepel)
# install.packages("indicspecies")
library(indicspecies)


## install iNEXT from github
# install.packages('devtools')
# library(devtools)
# devtools::install_github('AnneChao/iNEXT')

## import packages
library(iNEXT)

# install.packages("patchwork")
library(patchwork)

# remove.packages("gtable")
# install.packages("gtable")
library(gtable)

#### USER DEFINED VARIABLES ####

inFilePath1 = "./WorkingData_CLEANED_TUB,CAG.xlsx"
inFilePath2 = "./PHIRES_MetaData.xlsx"

#### READ IN DATA & CURATE ####

data <-
  read_excel(inFilePath1,
             na="NA") %>%
  clean_names() %>%
  # there is a different depth for CAG_024 in data vs metadata
  # solution: go with metadata depth, confirm with Gene & Rene
  mutate(family = str_to_title(family),
         genus = str_to_title(genus),
         species = str_to_lower(species),
         family_clean = case_when(
           family == "Epinephelidae" ~ "Serranidae",
           TRUE ~ family),
         taxon = str_c(family_clean,
                       genus,
                       species,
                       sep = "_")) %>%
  mutate(depth_m = case_when(op_code == "CAG_024" ~ 8,
                             TRUE ~ depth_m)) %>%
  # keep only max_n
  group_by(op_code,
           taxon) %>%
  filter(max_n == max(max_n)) %>%
  ungroup() %>%
  # remove duplicated rows
  distinct(op_code,
           taxon,
           .keep_all = TRUE) %>% 
  filter(family_clean !="Alopiidae",
         family_clean != "Sphyrnidae",
         family_clean != "Carcharhinidae")

metadata <-
  read_excel(inFilePath2,
             na="NA") %>%
  clean_names() %>%
  dplyr::rename(bait_weight_grams = weight_grams) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Deep Reef")))

#### SUBSET AND MODIFY DATA: remove taxa not identified to species level for maxN analysis ####

data_removed_sp <- 
  data %>%
  filter(species != "sp") %>% 
  # mutate(family = str_to_title(family),
  #        genus = str_to_title(genus),
  #        species = str_to_lower(species),
  #        family_clean = case_when(
  #          family == "Epinephelidae" ~ "Serranidae",
  #          TRUE ~ family)) %>%
  mutate(groupings = case_when(
    family == "Labridae" ~ "Cheilinus undulatus",
    family == "Sphyrnidae" ~ "Galeomorphii",
    family == "Carcharhinidae" ~ "Galeomorphii",
    family == "Alopiidae" ~ "Galeomorphii",
    family == "Epinephelidae" ~ "Serranidae",
    TRUE ~ family))%>%
  mutate(taxon = str_c(groupings,
                       genus,
                       species,
                       sep = "_")) %>%
  # mutate(depth_m = case_when(op_code == "CAG_024" ~ 8,
  #                            TRUE ~ depth_m)) %>%
  # keep only max_n
  group_by(op_code,
           taxon) %>%
  filter(max_n == max(max_n)) %>%
  ungroup() %>%
  # remove duplicated rows
  distinct(op_code,
           taxon,
           .keep_all = TRUE)

# View(data_removed_sp)


#### COMBINE DATA & METADATA ####

data_all <-
  data %>%
  left_join(metadata,
            by = c("op_code" = "opcode",
                   "depth_m" = "depth_m")) %>%
  # rearrange order of columns, metadata then data
  dplyr::select(op_code,
                site:long_e,
                depth_m,
                time_in:bait_weight_grams,
                everything()) %>%
  mutate(study_locations = case_when(
    site == "Cawili" ~ "CAGAYANCILLO",
    site == "Calusa" ~ "CAGAYANCILLO",
    site == "Cagayancillo" ~ "CAGAYANCILLO",
    site == "TUBBATAHA" ~ "TRNP")) %>%
  mutate(study_locations = factor(study_locations,
                                  levels = c("TRNP", 
                                             "CAGAYANCILLO")))

data_all_only_sp <- data_all %>%
  filter(species == "sp")

data_all_removed_sp <- 
  data_removed_sp %>%
  left_join(metadata,
            by = c("op_code" = "opcode",
                   "depth_m" = "depth_m")) %>%
  # rearrange order of columns, metadata then data
  dplyr::select(op_code,
                site:long_e,
                depth_m,
                time_in:bait_weight_grams,
                everything()) %>%
  mutate(study_locations = case_when(
    site == "Cawili" ~ "CAGAYANCILLO",
    site == "Calusa" ~ "CAGAYANCILLO",
    site == "Cagayancillo" ~ "CAGAYANCILLO",
    site == "TUBBATAHA" ~ "TRNP")) %>%
  mutate(study_locations = factor(study_locations, 
                                  levels = c("TRNP", 
                                             "CAGAYANCILLO")))

#### PREP DATA FOR VEGAN ####

data_vegan <-
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) 


# convert metadata into tibble for vegan ingestion
# each row is a site
# each column is a dimension of site, such as depth, lat, long, region, etc

data_vegan.env <-
  data_all_removed_sp %>%
  # sum all max_n counts for a taxon and op_code
  dplyr::select(taxon,
                op_code,
                site,
                study_locations,
                survey_area,
                habitat,
                lat_n,
                long_e,
                depth_m,
                survey_length_hrs,
                bait_type,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  dplyr::select(op_code:bait_type) %>%
  mutate(site_code = str_remove(op_code,
                                "_.*$"),
         site_code = factor(site_code),
         study_locations = factor(study_locations),
         habitat = factor(habitat),
         bait_type = factor(bait_type),
         site = factor(site),
         survey_area = factor(survey_area),
         habitat_mpa = str_c(habitat,
                             study_locations,
                             sep = " ")) %>%
  mutate(studylocation_habitat = str_c(site_code,
                                       habitat,
                                       sep = "_"))

# and now we "attach" the metadata to the data

# Upon reviewing the function of `attach` we should involk it just prior to running code that uses it, then involk `detach` when we no longer need the data to be attached.  this will prevent interference between different versions of the vegan data tibbles
# https://statisticsglobe.com/r-warning-the-following-objects-are-masked
#attach(data_vegan.env)

#### PREP DATA FOR VEGAN: TRNP & Cagayancillo ####

##Create data_vegan and data_vegan.env for TRNP

data_vegan_TRNP.env <- data_vegan.env %>%
  filter(study_locations == "TRNP")

data_vegan_TRNP <- bind_cols(data_vegan, data_vegan.env) %>%
  filter(study_locations == "TRNP") %>%
  dplyr::select(-op_code:-studylocation_habitat) %>%
  dplyr::select_if(colSums(.) != 0)


attach(data_vegan_TRNP.env)

##Create data_vegan and data_vegan.env for Cagayancillo
data_vegan_CAG.env <- data_vegan.env %>%
  filter(study_locations == "CAGAYANCILLO")

data_vegan_CAG <- bind_cols(data_vegan, data_vegan.env) %>%
  filter(study_locations == "CAGAYANCILLO") %>%
  dplyr::select(-op_code:-studylocation_habitat) %>%
  dplyr::select_if(colSums(.) != 0)


attach(data_vegan_CAG.env)

#### PREP DATA FOR VEGAN: shallow & deep reefs ####

##Create data_vegan for shallow reefs
data_vegan_shallow.env <- data_vegan.env %>%
  filter(habitat == "Shallow Reef")

data_vegan_shallow <- bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat == "Shallow Reef") %>% 
  select(-op_code:-studylocation_habitat) %>%
  select_if(colSums(.) != 0)

#attach(data_vegan_shallow.env)

##Create data_vegan for deep reefs
data_vegan_deep.env <- data_vegan.env %>%
  filter(habitat == "Deep Reef")

data_vegan_deep <- bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat == "Deep Reef") %>%
  select(-op_code:-studylocation_habitat) %>%
  select_if(colSums(.) != 0)

#attach(data_vegan_deep.env)

#### vegan::specpool - Extrapolated Species Richness in a Species Pool Based on Incidence (Presence Absence) ####

# view(data_vegan.env)
pool <- 
  with(data_vegan.env, 
       specpool(x = data_vegan, 
                pool = habitat_mpa,
                smallsample = TRUE))

# View(pool)

pool %>%
  clean_names() %>%
  # mutate(depth_cat = rownames(.),
  #        depth_cat = factor(depth_cat,
  #                           levels = c("<2m",
  #                                      "2-15m",
  #                                      ">15m"))) %>%
  pivot_longer(cols = chao:boot_se) %>%
  mutate(se_value = case_when(str_detect(name,
                                         "_se") ~ "se",
                              TRUE ~ "value"),
         name = str_remove(name,
                           "_se")) %>% 
  pivot_wider(names_from = se_value) %>%
  dplyr::rename(sp_richness_est = value,
                estimator = name) %>% 
  ggplot(aes(x= habitat_mpa,
             y= sp_richness_est)) +
  geom_col() +
  geom_point(aes(y = species),
             color = "red3") +
  geom_errorbar(aes(ymin = sp_richness_est - se,
                    ymax = sp_richness_est + se)) +
  theme_classic() +
  labs(title = "Extrapolated Species Richness - Presence/Absence") +
  facet_grid(. ~ estimator)


#### vegan::estaccumR Extrapolated Species Richness Curve in a Species Pool Based on Abundance ####

# see `specpool` in https://cloud.r-project.org/web/packages/vegan/vegan.pdf

# abundance based richness rarefaction curve
# creates 1 curve per data frame, so if you want multiple curves, have to make them separately then combine into 1 tibble to plot
# increase permutations to 999 if you use this for your project

#Check the run time
# system.time({
#   p <- vegan::estaccumR(data_vegan, permutations = 9999)
#   })

p <- vegan::estaccumR(data_vegan, permutations = 9999)
# filter(estimator != "ace")

data_estaccumR_plot <-
  p$chao %>% 
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(chao_mean = mean(value),
                   chao_ci_lower = quantile(value,
                                            probs = 0.025),
                   chao_ci_upper = quantile(value,
                                            probs = 0.975))

str(data_estaccumR_plot)

data_estaccumR_plots <- 
  data_estaccumR_plot %>%
  ggplot(aes(x=N,
             y=chao_mean)) +
  geom_ribbon(aes(ymin=chao_ci_lower,
                  ymax=chao_ci_upper),
              fill = "#ebe8f3") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Chao1 Species Richness") +
  theme_classic() +
  labs(title =  "Overall Species Richness")

data_estaccumR_plots
ggsave("SpeciesRarefactionCurveAbundance.png",
       data_estaccumR_plots)

##Old Plot
# p.plot<-plot(p, 
#              display = c("chao"))
# p.plot
# save_plot("SpeciesRichnessRarefactionCurveAbundance.png", 
#           p.plot)

## Species Rarefaction Curve for Shallow Reef at Cagayancillo ##
data_vegan_CAG_shallow <- 
  bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat_mpa == "Shallow Reef CAGAYANCILLO") %>%
  dplyr::select(colnames(data_vegan))

# data_vegan_CAG_shallow.env <- 
#   bind_cols(data_vegan, data_vegan.env) %>%
#   filter(habitat_mpa == "Shallow Reef CAGAYANCILLO") %>%
#   dplyr::select(-colnames(data_vegan))

p_CAG_shallow <- estaccumR(data_vegan_CAG_shallow, permutations = 9999)

p_CAG_shallow_plot <-
  p_CAG_shallow$chao %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(chao_mean = mean(value),
                   chao_ci_lower = quantile(value,
                                            probs = 0.025),
                   chao_ci_upper = quantile(value,
                                            probs = 0.975))
p_CAG_shallow_plots <- p_CAG_shallow_plot %>%
  ggplot(aes(x=N,
             y=chao_mean)) +
  geom_ribbon(aes(ymin=chao_ci_lower,
                  ymax=chao_ci_upper),
              fill = "#f8bac4") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Chao1 Species Richness") +
  labs(title ="Shallow Reef at Cagayancillo") +
  ylim(0,100) +
  theme_classic() 

p_CAG_shallow_plots

##Old Plot
# p_CAG_shallow_plot <- plot(p_CAG_shallow,
#                            display = c("chao"),
#                            main = "Shallow Reef at Cagayancillo")
# p_CAG_shallow_plot

## Species Rarefaction Curve for Deep Reef at Cagayancillo ##
data_vegan_CAG_deep <- 
  bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat_mpa == "Deep Reef CAGAYANCILLO") %>%
  dplyr::select(colnames(data_vegan))

p_CAG_deep <- estaccumR(data_vegan_CAG_deep, permutations = 9999)

p_CAG_deep_plot <-
  p_CAG_deep$chao %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(chao_mean = mean(value),
                   chao_ci_lower = quantile(value,
                                            probs = 0.025),
                   chao_ci_upper = quantile(value,
                                            probs = 0.975))
p_CAG_deep_plots <- p_CAG_deep_plot %>%
  ggplot(aes(x=N,
             y=chao_mean)) +
  geom_ribbon(aes(ymin=chao_ci_lower,
                  ymax=chao_ci_upper),
              fill = "#c5d8ea") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Chao1 Species Richness") +
  labs(title ="Mesophotic Reef at Cagayancillo") +
  ylim(0,100) +
  theme_classic()

p_CAG_deep_plots

##Old Plot
# p_CAG_deep_plot <- plot(p_CAG_deep,
#                            display = c("chao"),
#                         main = "Mesophotic Reef at Cagayancillo")
# 
# p_CAG_deep_plot

## Species Rarefaction Curve for Shallow Reef at TRNP ##
data_vegan_TUB_shallow <- 
  bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat_mpa == "Shallow Reef TRNP") %>%
  dplyr::select(colnames(data_vegan))

p_TUB_shallow <- estaccumR(data_vegan_TUB_shallow, permutations = 9999)

p_TUB_shallow_plot <-
  p_TUB_shallow$chao %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(chao_mean = mean(value),
                   chao_ci_lower = quantile(value,
                                            probs = 0.025),
                   chao_ci_upper = quantile(value,
                                            probs = 0.975))
p_TUB_shallow_plots <- p_TUB_shallow_plot %>%
  ggplot(aes(x=N,
             y=chao_mean)) +
  geom_ribbon(aes(ymin=chao_ci_lower,
                  ymax=chao_ci_upper),
              fill = "#f8bac4") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Chao1 Species Richness") +
  labs(title ="Shallow Reef at TRNP") +
  ylim(0,100) +
  theme_classic()

p_TUB_shallow_plots

##Old Plot

# p_TUB_shallow_plot <- plot(p_TUB_shallow,
#                         display = c("chao"),
#                         main = "Shallow Reef at TRNP")
# p_TUB_shallow_plot

## Species Rarefaction Curve for Deep Reef at TRNP ##
data_vegan_TUB_deep <- 
  bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat_mpa == "Deep Reef TRNP") %>%
  dplyr::select(colnames(data_vegan))

p_TUB_deep <- estaccumR(data_vegan_TUB_deep, permutations = 9999)
p_TUB_deep_plot <-
  p_TUB_deep$chao %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(chao_mean = mean(value),
                   chao_ci_lower = quantile(value,
                                            probs = 0.025),
                   chao_ci_upper = quantile(value,
                                            probs = 0.975))
p_TUB_deep_plots <- p_TUB_deep_plot %>%
  ggplot(aes(x=N,
             y=chao_mean)) +
  geom_ribbon(aes(ymin=chao_ci_lower,
                  ymax=chao_ci_upper),
              fill = "#c5d8ea") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Chao1 Species Richness") +
  labs(title ="Mesophotic Reef at TRNP") +
  ylim(0,100) +
  theme_classic()

p_TUB_deep_plots

##Old Plot
# p_TUB_deep_plot <- plot(p_TUB_deep,
#                            display = c("chao"),
#                         main = "Mesophotic Reef at TRNP")
# 
# p_TUB_deep_plot



##All Habitat and Study Locations Chao Plot ##
p_habitat_locations_plot <- ggarrange(p_TUB_shallow_plots,
                                      p_TUB_deep_plots,
                                      p_CAG_shallow_plots,
                                      p_CAG_deep_plots,
                                      ncol = 2,
                                      nrow = 2)

# dev.off()
ggsave("SpeciesRichnessRarefactionCurveHabitatandStudyLocations.png", 
       p_habitat_locations_plot,
       height = 11,
       width = 8.5,
       units = "in")

## Species Rarefaction Curve for Serranidae ##
data_vegan_Serranidae <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Serranidae"))


# data_vegan_Serranidae_nozeros <- data_vegan_Serranidae %>%
#                             select_if(colSums(.) != 0)
# 
# p_Serranidae_nz <- estaccumR(data_vegan_Serranidae_nozeros, permutations = 999)
# p_Serranidae_plot_nz <- plot(p_Serranidae_nz,
#                           display = c("chao"),
#                           main = "Serranidae")
# p_Serranidae_plot_nz



p_Serranidae <- estaccumR(data_vegan_Serranidae, permutations = 9999)
View(p_Serranidae)

p_Serranidae_plot <-
  p_Serranidae$S %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>% 
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))
p_Serranidae_plots <- p_Serranidae_plot %>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#ebe8f3") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Serranidae") +
  ylim(0,20) +
  theme_classic()

p_Serranidae_plots

##Old Plots

# p_Serranidae_plot <- plot(p_Serranidae,
#                         display = c("chao"),
#                         main = "Serranidae")
# 
# p_Serranidae_plot

## Species Rarefaction Curve for Lutjanidae ##
data_vegan_Lutjanidae <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Lutjanidae"))

p_Lutjanidae <- estaccumR(data_vegan_Lutjanidae, permutations = 9999)

p_Lutjanidae_plot <-
  p_Lutjanidae$S %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))
p_Lutjanidae_plots <- p_Lutjanidae_plot %>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#ebe8f3") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Lutjanidae") +
  ylim(0,20) +
  theme_classic()

p_Lutjanidae_plots

##Old Plots
# p_Lutjanidae_plot <- plot(p_Lutjanidae,
#                           display = c("chao"),
#                           main = "Lutjanidae")
# 
# p_Lutjanidae_plot

## Species Rarefaction Curve for Lethrinidae ##
data_vegan_Lethrinidae <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Lethrinidae"))

p_Lethrinidae <- estaccumR(data_vegan_Lethrinidae, permutations = 9999)
p_Lethrinidae_plot <-
  p_Lethrinidae$S %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))
p_Lethrinidae_plots <- p_Lethrinidae_plot %>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#ebe8f3") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Lethrinidae") +
  ylim(0,20) +
  theme_classic()

p_Lethrinidae_plots

##Old Plots
# p_Lethrinidae_plot <- plot(p_Lethrinidae,
#                           display = c("chao"),
#                           main = "Lethrinidae")
# 
# p_Lethrinidae_plot

## Species Rarefaction Curve for Carangidae ##
data_vegan_Carangidae <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Carangidae"))

p_Carangidae <- estaccumR(data_vegan_Carangidae, permutations = 9999)
p_Carangidae_plot <-
  p_Carangidae$S %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))
p_Carangidae_plots <- p_Carangidae_plot %>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#ebe8f3") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Carangidae") +
  ylim(0,20) +
  theme_classic()

p_Carangidae_plots

#Old Plots
# p_Carangidae_plot <- plot(p_Carangidae,
#                            display = c("chao"),
#                           main = "Carangidae")
# 
# p_Carangidae_plot


## Species Rarefaction Curve for Galeomorphii ##
data_vegan_Galeomorphii <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Galeomorphii"))

view(data_vegan_Galeomorphii)

p_Galeomorphii <- estaccumR(data_vegan_Galeomorphii, permutations = 999)

p_Galeomorphii_plot <-
  p_Galeomorphii$S %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(chao_mean = mean(value),
                   chao_ci_lower = quantile(value,
                                            probs = 0.025),
                   chao_ci_upper = quantile(value,
                                            probs = 0.975))
p_Galeomorphii_plots <- p_Galeomorphii_plot %>%
  ggplot(aes(x=N,
             y=chao_mean)) +
  geom_ribbon(aes(ymin=chao_ci_lower,
                  ymax=chao_ci_upper),
              fill = "#ebe8f3") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Galeomorphii") +
  ylim(0,20) +
  theme_classic()

p_Galeomorphii_plots

##Old Plots
# p_Galeomorphii_plot <- plot(p_Galeomorphii,
#                           display = c("chao"),
#                           main = "Galeomorphii")
# 
# p_Galeomorphii_plot

##Species Rarefaction Curve for Cheilinus undulatus ##
data_vegan_Cheilinus_undulatus <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Cheilinus undulatus"))

View(data_vegan_Cheilinus_undulatus)

p_Cheilinus_undulatus <- estaccumR(data_vegan_Cheilinus_undulatus, permutations = 999)

p_Cheilinus_undulatus_plot <- plot(p_Cheilinus_undulatus,
                                   display = c("chao"),
                                   main = "Cheilinus undulatus")

p_Cheilinus_undulatus_plot


## Species Rarefaction Curves for All Groupings ##
p_groupings_plot <- ggpubr::ggarrange(data_estaccumR_plots,
                                      p_Serranidae_plots,
                                      p_Lutjanidae_plots,
                                      p_Lethrinidae_plots,
                                      p_Carangidae_plots,
                                      ncol = 2,
                                      nrow = 3)
p_groupings_plot

ggsave("SpeciesRichnessRarefactionCurveGroupings.png", 
       p_groupings_plot, 
       height = 11,
       width = 8.5,
       units = "in")

## Species 
# to make plot w ggplot, see https://stackoverflow.com/questions/52652195/convert-rarefaction-plots-from-vegan-to-ggplot2-in-r

# #### iNEXT::iNEXT Extrapolated Species Richness Curve in a Species Pool Based on Abundance ####
# 
# # note, iNext can only provide one curve per video using abundance data.
# # it has no equivalent function to estaccumR
# # the iNEXT abundance-based rarefaction curves have number of fish on the x axis, whereas
# # estaccumR has number of videos. They are fundamentally different.
# # where the vegan functions do overlap with iNEXT is the chao1 ests by video
# 
# inext_output <-
#   t(data_vegan) %>%
#   iNEXT(
#     q=0, 
#     datatype="abundance",
#     endpoint = 100
#   )
# 
# # ggiNEXT(inext_output, type=1)
# 
# inext_output$AsyEst %>%
#   clean_names() %>%
#   filter(diversity == "Species richness") %>%
#   mutate(
#     # zero pad the assemblage names so they can be sorted to match data_vegan.env
#     assemblage = str_replace(
#       assemblage, "\\d+", 
#       str_pad(
#         as.numeric(
#           str_extract(
#             assemblage, 
#             "\\d+"
#           )
#         ), 
#         width = 2, 
#         pad = "0"
#       )
#     ) 
#   ) %>%
#   arrange(assemblage) %>%
#   bind_cols(data_vegan.env) %>% #pull("depth_m")
#   ggplot() +
#   aes(
#     x=depth_m,
#     y=estimator,
#     color = site_code
#   ) +
#   geom_point() +
#   
#   geom_errorbar(aes(ymin = lcl,
#                     ymax = ucl)) +
#   theme_classic() 
# 
# 
# # plot size_based rarefaction and extrapolation
# # http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/iNEXT_Introduction.pdf
# 
# # Sample-size-based R/E sampling curves: iNEXT computes diversity estimates for
# # rarefied and extrapolated samples up to double the reference sample size (by
# # default) or a user-specified size. This type of sampling curve plots the diversity
# # estimates with respect to sample size. Sample size refers to the number of
# # individuals in a sample for abundance data, whereas it refers to the number of
# # sampling units for incidence data.
# 
# inext_output$iNextEst$size_based %>%
#   clean_names() %>%
#   mutate(
#     # zero pad the assemblage names so they can be sorted to match data_vegan.env
#     assemblage = str_replace(
#       assemblage, "\\d+", 
#       str_pad(
#         as.numeric(
#           str_extract(
#             assemblage, 
#             "\\d+"
#           )
#         ), 
#         width = 2, 
#         pad = "0"
#       )
#     ),
#     method = factor(
#       method,
#       levels = c(
#         "Rarefaction",
#         "Extrapolation",
#         "Observed"
#       )
#     )
#   ) %>%
#   filter(
#     order_q == 0,
#     method != "Observed") %>%
#   left_join(
#     data_vegan.env %>%
#       mutate(
#         assemblage = sprintf(
#           "assemblage%02d", 
#           1:n()
#         )
#       ) 
#   ) %>% #filter(method == "Rarefaction") %>% pull(m) %>% max()
#   ggplot() +
#   aes(
#     x=m,
#     y=q_d,
#     color = site_code,
#     linetype = method,
#     group = interaction(assemblage,method)
#   ) +
#   geom_ribbon(aes(ymin = q_d_lcl,
#                   ymax = q_d_ucl),
#               color = "grey90",
#               alpha = 0.2) +
#   geom_line(size = 1) +
#   theme_classic() +
#   facet_grid(site_code ~ habitat)  
# 
# 
# # plot coverage_based rarefaction and extrapolation
# # http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/iNEXT_Introduction.pdf
# 
# # Coverage-based R/E sampling curves: iNEXT computes diversity estimates for
# # rarefied and extrapolated samples with sample completeness (as measured by
# # sample coverage) up to the coverage value of double the reference sample size (by
# # default) or a user-specified coverage. This type of sampling curve plots the
# # diversity estimates with respect to sample coverage.
# 
# inext_output$iNextEst$size_based %>%
#   clean_names() %>%
#   mutate(
#     # zero pad the assemblage names so they can be sorted to match data_vegan.env
#     assemblage = str_replace(
#       assemblage, "\\d+", 
#       str_pad(
#         as.numeric(
#           str_extract(
#             assemblage, 
#             "\\d+"
#           )
#         ), 
#         width = 2, 
#         pad = "0"
#       )
#     ),
#     method = factor(
#       method,
#       levels = c(
#         "Rarefaction",
#         "Extrapolation",
#         "Observed"
#       )
#     )
#   ) %>%
#   filter(
#     order_q == 0,
#     method != "Observed") %>%
#   left_join(
#     data_vegan.env %>%
#       mutate(
#         assemblage = sprintf(
#           "assemblage%02d", 
#           1:n()
#         )
#       ) 
#   ) %>%
#   ggplot() +
#   aes(
#     x=m,
#     y=sc,
#     color = site_code,
#     linetype = method,
#     group = interaction(assemblage,method)
#   ) +
#   geom_ribbon(aes(ymin = sc_lcl,
#                   ymax = sc_ucl),
#               color = "grey90",
#               alpha = 0.5) +
#   geom_line(size = 1) +
#   theme_classic() +
#   facet_grid(site_code ~ habitat) 
# 
# inext_output$iNextEst$coverage_based %>%
#   clean_names() %>%
#   mutate(
#     # zero pad the assemblage names so they can be sorted to match data_vegan.env
#     assemblage = str_replace(
#       assemblage, "\\d+", 
#       str_pad(
#         as.numeric(
#           str_extract(
#             assemblage, 
#             "\\d+"
#           )
#         ), 
#         width = 2, 
#         pad = "0"
#       )
#     ),
#     method = factor(
#       method,
#       levels = c(
#         "Rarefaction",
#         "Extrapolation",
#         "Observed"
#       )
#     )
#   ) %>%
#   filter(
#     order_q == 0,
#     method != "Observed") %>%
#   left_join(
#     data_vegan.env %>%
#       mutate(
#         assemblage = sprintf(
#           "assemblage%02d", 
#           1:n()
#         )
#       ) 
#   ) %>%
#   ggplot() +
#   aes(
#     x=sc,
#     y=q_d,
#     color = site_code,
#     linetype = method,
#     group = interaction(assemblage,method)
#   ) +
#   geom_ribbon(aes(ymin = q_d_lcl,
#                   ymax = q_d_ucl),
#               color = "grey90",
#               alpha = 0.5) +
#   geom_line(size = 1) +
#   theme_classic() +
#   facet_grid(site_code ~ habitat) 
# 
# 
# 
#### iNEXT::iNEXT Extrapolated Species Richness Curve in a Species Pool Based on Incidence ####
# https://chat.openai.com/share/1573861d-142f-4e7f-ba40-dc34bbd3e2d1

# all data
inext_output <- 
  data_vegan %>%
  #convert to incidence data, 0 or 1, rather than abundance
  mutate(
    across(
      everything(),
      ~if_else(. > 0, 1, 0)
    )
  ) %>%
  # convert to format expected by iNEXT
  # Incidence-raw data (datatype="incidence_raw"): for each assemblage, input
  # data for a reference sample consist of a species-by-sampling-unit matrix; when
  # there are N assemblages, input data consist of N lists of matrices, and each matrix
  # is a species-by-sampling-unit matrix.
  t() %>%
  list() %>%
  iNEXT(
    q=0, 
    datatype="incidence_raw"
  )

ggiNEXT(inext_output, 
        type=1) +
  theme_classic() +
  labs(y = "Species Richness",
       x = "Number of BRUVs") +
  theme(legend.position = "none")
# guides(color = FALSE)

ggiNEXT(inext_output,
        type=2)
ggiNEXT(inext_output,
        type=3)


# plot size_based rarefaction and extrapolation
# http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/iNEXT_Introduction.pdf

# Sample-size-based R/E sampling curves: iNEXT computes diversity estimates for
# rarefied and extrapolated samples up to double the reference sample size (by
# default) or a user-specified size. This type of sampling curve plots the diversity
# estimates with respect to sample size. Sample size refers to the number of
# individuals in a sample for abundance data, whereas it refers to the number of
# sampling units for incidence data.

##Chris's Original Code

inextoutput_overall <- inext_output$iNextEst$size_based %>%
  clean_names() %>%
  mutate(
    # zero pad the assemblage names so they can be sorted to match data_vegan.env
    assemblage = str_replace(
      assemblage, "\\d+", 
      str_pad(
        as.numeric(
          str_extract(
            assemblage, 
            "\\d+"
          )
        ), 
        width = 2, 
        pad = "0"
      )
    ),
    method = factor(
      method,
      levels = c(
        "Rarefaction",
        "Extrapolation",
        "Observed"
      )
    )
  ) %>%
  filter(
    order_q == 0,
    method != "Observed"
  ) %>%
  left_join(
    data_vegan.env %>%
      mutate(
        assemblage = sprintf(
          "assemblage%02d", 
          1:n()
        )
      ) 
  ) %>% 
  ggplot() +
  aes(
    x=t,
    y=q_d,
    # color = site_code,
    linetype = method,
    group = interaction(assemblage,method)
  ) +
  geom_ribbon(aes(ymin = q_d_lcl,
                  ymax = q_d_ucl,
                  fill = method),
               # color = "black",
              alpha = 0.2) +
  geom_line(size = 1) +
  geom_point(
    data = inext_output$iNextEst$size_based %>%
      clean_names() %>%
      filter(method == "Observed"),
    size = 3
  ) +
  scale_fill_manual(values = c("#9C8EC4", "#9C8EC4"))+
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Species Richness",
       x = "Number of BRUVs")
# facet_grid(site_code ~ habitat)  

inextoutput_overall
##Code I took from Later in Taxon but applied here bc the fill is missing around the Observed Point in Chris's original code

inext_overall_u <- 
  inext_output$iNextEst$size_based %>%
  clean_names() %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Rarefaction",
        "Extrapolation",
        "Observed"
      )
    )
  ) %>%
  filter(
    order_q == 0,
    method != "Observed"
  ) %>% 
  ggplot() +
  aes(
    x=t,
    y=q_d,
    # fill = taxon,
  ) +
  geom_ribbon(
    aes(ymin = q_d_lcl,
        ymax = q_d_ucl),
    alpha = 0.5,
    fill = "#DBD6EA"
  ) +
  geom_line(size = 1,
            aes(linetype = method)) +
  geom_point(
    data = inext_output$iNextEst$size_based %>%
      clean_names() %>%
      filter(method == "Observed"),
    size = 3
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Species Richness",
       x = "Number of BRUVs")

inext_overall_u

# # Check data for observed
# Print data for the "Observed" method
# observed_data <- inext_output$iNextEst$size_based %>%
#   clean_names() %>%
#   filter(method == "Observed")
# print(observed_data)
# 
# summary_stats <- inext_output$iNextEst$size_based %>%
#   clean_names() %>%
#   filter(method == "Observed") %>%
#   summarise(
#     min_q_d_lcl = min(q_d_lcl),
#     max_q_d_lcl = max(q_d_lcl),
#     min_q_d_ucl = min(q_d_ucl),
#     max_q_d_ucl = max(q_d_ucl)
#   )
# print(summary_stats)

ggsave("iNEXT Incidence Overall Species Richness Curves.png",
       inext_overall_u,
       height = 6,
       width = 8,
       units = "in")


# plot coverage_based rarefaction and extrapolation
# http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/iNEXT_Introduction.pdf

# Coverage-based R/E sampling curves: iNEXT computes diversity estimates for
# rarefied and extrapolated samples with sample completeness (as measured by
# sample coverage) up to the coverage value of double the reference sample size (by
# default) or a user-specified coverage. This type of sampling curve plots the
# diversity estimates with respect to sample coverage.

inext_output$iNextEst$size_based %>%
  clean_names() %>%
  mutate(
    # zero pad the assemblage names so they can be sorted to match data_vegan.env
    assemblage = str_replace(
      assemblage, "\\d+", 
      str_pad(
        as.numeric(
          str_extract(
            assemblage, 
            "\\d+"
          )
        ), 
        width = 2, 
        pad = "0"
      )
    ),
    method = factor(
      method,
      levels = c(
        "Rarefaction",
        "Extrapolation",
        "Observed"
      )
    )
  ) %>%
  filter(
    order_q == 0,
    method != "Observed") %>%
  left_join(
    data_vegan.env %>%
      mutate(
        assemblage = sprintf(
          "assemblage%02d", 
          1:n()
        )
      ) 
  ) %>%
  ggplot() +
  aes(
    x=t,
    y=sc,
    # color = site_code,
    linetype = method,
    group = interaction(assemblage,method)
  ) +
  geom_ribbon(aes(ymin = sc_lcl,
                  ymax = sc_ucl),
              color = "grey90",
              alpha = 0.2) +
  geom_line(size = 1) +
  geom_point(
    data = inext_output$iNextEst$size_based %>%
      clean_names() %>%
      filter(method == "Observed"),
    size = 3
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Sample Coverage",
       x = "Number of Individuals")
# facet_grid(site_code ~ habitat) 

inext_output$iNextEst$coverage_based %>%
  clean_names() %>%
  mutate(
    # zero pad the assemblage names so they can be sorted to match data_vegan.env
    assemblage = str_replace(
      assemblage, "\\d+", 
      str_pad(
        as.numeric(
          str_extract(
            assemblage, 
            "\\d+"
          )
        ), 
        width = 2, 
        pad = "0"
      )
    ),
    method = factor(
      method,
      levels = c(
        "Rarefaction",
        "Extrapolation",
        "Observed"
      )
    )
  ) %>%
  filter(
    order_q == 0,
    method != "Observed") %>%
  left_join(
    data_vegan.env %>%
      mutate(
        assemblage = sprintf(
          "assemblage%02d", 
          1:n()
        )
      ) 
  ) %>%
  ggplot() +
  aes(
    x=sc,
    y=q_d,
    # color = site_code,
    linetype = method,
    group = interaction(assemblage,method)
  ) +
  geom_ribbon(aes(ymin = q_d_lcl,
                  ymax = q_d_ucl),
              color = "grey90",
              alpha = 0.2) +
  geom_line(size = 1) +
  geom_point(
    data = inext_output$iNextEst$size_based %>%
      clean_names() %>%
      filter(method == "Observed"),
    size = 3
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Species Richness",
       x = "Sample Coverage")
# facet_grid(site_code ~ habitat) 





#### iNEXT::iNEXT Extrapolated Species Richness Curve in a Species Pool Based on Incidence: Location x Depth ####
# https://chat.openai.com/share/1573861d-142f-4e7f-ba40-dc34bbd3e2d1

# location x depth
abundance2incidence <-
  function(data){
    data %>%
      mutate(
        across(
          everything(),
          ~if_else(. > 0, 1, 0)
        )
      ) %>%
      t()
  }

inext_output <- 
  list(
    data_vegan_CAG_deep %>%
      abundance2incidence(),
    data_vegan_CAG_shallow %>%
      abundance2incidence(),
    data_vegan_TUB_deep %>%
      abundance2incidence(),
    data_vegan_TUB_shallow %>%
      abundance2incidence()
  ) %>%
  iNEXT(
    q=0, 
    datatype="incidence_raw",
    endpoint = 40
  )

ggiNEXT(
  inext_output, 
  type=1,
  facet.var = "Assemblage"
) +
  theme_classic() +
  labs(y = "Species Richness",
       x = "Number of BRUVs") +
  theme(legend.position = "none")

ggiNEXT(inext_output,
        type=2,
        facet.var = "Assemblage")
ggiNEXT(inext_output,
        type=3,
        facet.var = "Assemblage")

# edit iNEXT output to add our treatments, location and depth_cat
add_treatment <-
  function(data){
    data %>%
      mutate(location = case_when(str_detect(Assemblage,
                                             "[12]") ~
                                    "cag",
                                  str_detect(Assemblage,
                                             "[34]") ~
                                    "tub",
                                  TRUE ~ NA_character_),
             depth_cat = case_when(str_detect(Assemblage,
                                              "[13]") ~
                                     "mesophotic",
                                   str_detect(Assemblage,
                                              "[24]") ~
                                     "shallow",
                                   TRUE ~ NA_character_),
             location = factor(location,
                               levels = c("tub",
                                          "cag")),
             depth_cat = factor(depth_cat,
                                levels = c("shallow",
                                           "mesophotic")))
  }

inext_output$DataInfo <-
  inext_output$DataInfo %>%
  add_treatment()

inext_output$iNextEst$size_based <-
  inext_output$iNextEst$size_based %>%
  add_treatment()

inext_output$iNextEst$coverage_based <-
  inext_output$iNextEst$coverage_based %>%
  add_treatment()

inext_output$AsyEst <-
  inext_output$AsyEst %>%
  add_treatment()

# ggiNEXT has too much hardcoding so even if I add variables, it can't handle them
ggiNEXT(
  inext_output,
  type=1,
  color.var = "depth_cat",
  facet.var = "Assemblage"
) +
  theme_classic() +
  labs(y = "Species Richness",
       x = "Number of BRUVs") +
  theme(legend.position = "none")


# plot size_based rarefaction and extrapolation
# http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/iNEXT_Introduction.pdf

# Sample-size-based R/E sampling curves: iNEXT computes diversity estimates for
# rarefied and extrapolated samples up to double the reference sample size (by
# default) or a user-specified size. This type of sampling curve plots the diversity
# estimates with respect to sample size. Sample size refers to the number of
# individuals in a sample for abundance data, whereas it refers to the number of
# sampling units for incidence data.

inextoutput_studylocation_depth <- 
  inext_output$iNextEst$size_based %>%
  clean_names() %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Rarefaction",
        "Extrapolation",
        "Observed"
      )
    )
  ) %>%
  filter(
    order_q == 0,
    method != "Observed"
  ) %>% 
  ggplot() +
  aes(
    x=t,
    y=q_d,
    fill = depth_cat,
  ) +
  geom_ribbon(aes(ymin = q_d_lcl,
                  ymax = q_d_ucl),
              alpha = 0.5) +
  geom_line(size = 1,
            aes(linetype = method)) +
  geom_point(
    data = inext_output$iNextEst$size_based %>%
      clean_names() %>%
      filter(method == "Observed"),
    size = 3
  ) +
  scale_fill_manual(values = c("#F8BAC4",
                              "#C5D8EA",
                              "#F8BAC4",
                              "#C5D8EA"
                              )) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Species Richness",
       x = "Number of BRUVs") +
  # facet_grid(location ~ depth_cat,
  #            scales = "free")
  facet_wrap(
    location ~ depth_cat,
    scales = "free"  # Set scales to "free" for both x and y axes
  ) +
  ylim(0, 80)

inextoutput_studylocation_depth

ggsave("iNEXTrarefactioncurves_location_depth.png",
       inextoutput_studylocation_depth,
       height = 8,
       width = 8,
       units = "in")

# plot coverage_based rarefaction and extrapolation
# http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/iNEXT_Introduction.pdf

# Coverage-based R/E sampling curves: iNEXT computes diversity estimates for
# rarefied and extrapolated samples with sample completeness (as measured by
# sample coverage) up to the coverage value of double the reference sample size (by
# default) or a user-specified coverage. This type of sampling curve plots the
# diversity estimates with respect to sample coverage.

inext_output$iNextEst$size_based %>%
  clean_names() %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Rarefaction",
        "Extrapolation",
        "Observed"
      )
    )
  ) %>%
  filter(
    order_q == 0,
    method != "Observed") %>%
  ggplot() +
  aes(
    x=t,
    y=sc,
    fill = depth_cat,
  ) +
  geom_ribbon(aes(ymin = sc_lcl,
                  ymax = sc_ucl),
              color = "grey90",
              alpha = 0.5) +
  geom_line(aes(linetype = method),
            size = 1) +
  geom_point(
    data = inext_output$iNextEst$size_based %>%
      clean_names() %>%
      filter(method == "Observed"),
    size = 3
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Sample Coverage",
       x = "Number of Individuals") +
  facet_grid(location ~ depth_cat)

inext_output$iNextEst$coverage_based %>%
  clean_names() %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Rarefaction",
        "Extrapolation",
        "Observed"
      )
    )
  ) %>%
  filter(
    order_q == 0,
    method != "Observed"
  ) %>%
  ggplot() +
  aes(
    x=sc,
    y=q_d,
    fill = depth_cat,
  ) +
  geom_ribbon(aes(ymin = q_d_lcl,
                  ymax = q_d_ucl),
              color = "grey90",
              alpha = 0.2) +
  geom_line(aes(linetype = method),
            size = 1) +
  geom_point(
    data = inext_output$iNextEst$size_based %>%
      clean_names() %>%
      filter(method == "Observed"),
    size = 3
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Species Richness",
       x = "Sample Coverage") +
  facet_grid(location ~ depth_cat)


#### iNEXT::iNEXT Extrapolated Species Richness Curve in a Species Pool Based on Incidence: Taxon ####
# https://chat.openai.com/share/1573861d-142f-4e7f-ba40-dc34bbd3e2d1

# location x depth
abundance2incidence <-
  function(data){
    data %>%
      mutate(
        across(
          everything(),
          ~if_else(. > 0, 1, 0)
        )
      ) %>%
      t()
  }

inext_output <- 
  list(
    data_vegan_Serranidae %>%
      abundance2incidence(),
    data_vegan_Lutjanidae %>%
      abundance2incidence(),
    data_vegan_Lethrinidae %>%
      abundance2incidence(),
    data_vegan_Carangidae %>%
      abundance2incidence()
  ) %>%
  iNEXT(
    q=0, 
    datatype="incidence_raw"
  )

ggiNEXT(
  inext_output, 
  type=1,
  facet.var = "Assemblage"
) +
  theme_classic() +
  labs(y = "Species Richness",
       x = "Number of BRUVs") +
  theme(legend.position = "none")

ggiNEXT(inext_output,
        type=2,
        facet.var = "Assemblage")
ggiNEXT(inext_output,
        type=3,
        facet.var = "Assemblage")

# edit iNEXT output to add our treatments, location and depth_cat
add_taxon <-
  function(data){
    data %>%
      mutate(taxon = case_when(str_detect(Assemblage,
                                          "[1]") ~
                                 "Serranidae",
                               str_detect(Assemblage,
                                          "[2]") ~
                                 "Lutjanidae",
                               str_detect(Assemblage,
                                          "[3]") ~
                                 "Lethrinidae",
                               str_detect(Assemblage,
                                          "[4]") ~
                                 "Carangidae",
                               TRUE ~ NA_character_))
  }

inext_output$DataInfo <-
  inext_output$DataInfo %>%
  add_taxon()

inext_output$iNextEst$size_based <-
  inext_output$iNextEst$size_based %>%
  add_taxon()

inext_output$iNextEst$coverage_based <-
  inext_output$iNextEst$coverage_based %>%
  add_taxon()

inext_output$AsyEst <-
  inext_output$AsyEst %>%
  add_taxon()

# ggiNEXT has too much hardcoding so even if I add variables, it can't handle them
ggiNEXT(
  inext_output,
  type=1,
  color.var = "taxon",
  facet.var = "Assemblage"
) +
  theme_classic() +
  labs(y = "Species Richness",
       x = "Number of BRUVs") +
  theme(legend.position = "none")


# plot size_based rarefaction and extrapolation
# http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/iNEXT_Introduction.pdf

# Sample-size-based R/E sampling curves: iNEXT computes diversity estimates for
# rarefied and extrapolated samples up to double the reference sample size (by
# default) or a user-specified size. This type of sampling curve plots the diversity
# estimates with respect to sample size. Sample size refers to the number of
# individuals in a sample for abundance data, whereas it refers to the number of
# sampling units for incidence data.


inext_taxon <- 
  inext_output$iNextEst$size_based %>%
  clean_names() %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Rarefaction",
        "Extrapolation",
        "Observed"
      )
    )
  ) %>%
  filter(
    order_q == 0,
    method != "Observed"
  ) %>% 
  ggplot() +
  aes(
    x=t,
    y=q_d,
    # fill = taxon,
  ) +
  geom_ribbon(
    aes(ymin = q_d_lcl,
        ymax = q_d_ucl),
    alpha = 0.5,
    fill = "#DBD6EA"
  ) +
  geom_line(size = 1,
            aes(linetype = method)) +
  geom_point(
    data = inext_output$iNextEst$size_based %>%
      clean_names() %>%
      filter(method == "Observed"),
    size = 3
  ) +
  theme_classic() +
  # theme(legend.position = "none") +
  labs(y = "Species Richness",
       x = "Number of BRUVs") +
  # facet_wrap(~factor(taxon, levels = c("Serranidae", "Lutjanidae", "Lethrinidae", "Carangidae")))
  facet_wrap(
    ~factor(taxon, levels = c("Serranidae", "Lutjanidae", "Lethrinidae", "Carangidae")),
    scales = "free", # Set scales to "free" to have individual x and y-axes
    ncol = 2 # Adjust the number of columns as needed
  )

inext_taxon

ggsave("inext_RarefactionCurves_taxonGroupings.png", 
       inext_taxon,
       height = 8,
       width = 11.5,
       units = "in")

# plot coverage_based rarefaction and extrapolation
# http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/iNEXT_Introduction.pdf

# Coverage-based R/E sampling curves: iNEXT computes diversity estimates for
# rarefied and extrapolated samples with sample completeness (as measured by
# sample coverage) up to the coverage value of double the reference sample size (by
# default) or a user-specified coverage. This type of sampling curve plots the
# diversity estimates with respect to sample coverage.

inext_output$iNextEst$size_based %>%
  clean_names() %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Rarefaction",
        "Extrapolation",
        "Observed"
      )
    )
  ) %>%
  filter(
    order_q == 0,
    method != "Observed") %>%
  ggplot() +
  aes(
    x=t,
    y=sc,
    fill = depth_cat,
  ) +
  geom_ribbon(aes(ymin = sc_lcl,
                  ymax = sc_ucl),
              color = "grey90",
              alpha = 0.5) +
  geom_line(aes(linetype = method),
            size = 1) +
  geom_point(
    data = inext_output$iNextEst$size_based %>%
      clean_names() %>%
      filter(method == "Observed"),
    size = 3
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Sample Coverage",
       x = "Number of Individuals") +
  facet_grid(location ~ depth_cat)

inext_output$iNextEst$coverage_based %>%
  clean_names() %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Rarefaction",
        "Extrapolation",
        "Observed"
      )
    )
  ) %>%
  filter(
    order_q == 0,
    method != "Observed"
  ) %>%
  ggplot() +
  aes(
    x=sc,
    y=q_d,
    fill = depth_cat,
  ) +
  geom_ribbon(aes(ymin = q_d_lcl,
                  ymax = q_d_ucl),
              color = "grey90",
              alpha = 0.2) +
  geom_line(aes(linetype = method),
            size = 1) +
  geom_point(
    data = inext_output$iNextEst$size_based %>%
      clean_names() %>%
      filter(method == "Observed"),
    size = 3
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Species Richness",
       x = "Sample Coverage") +
  facet_grid(location ~ depth_cat)




#### vegan::estaccumR using species observations instead of Chao ####
p_obs<- vegan::estaccumR(data_vegan, permutations = 999)
# filter(estimator != "ace")

data_estaccumR_obs <-
  p$S %>% 
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))


data_estaccumR_obs_plot <- data_estaccumR_obs %>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#ebe8f3") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Chao1 Species Richness") +
  theme_classic() +
  labs(title =  "Overall Species Richness") +
  ylim(0,65)

data_estaccumR_obs_plot
ggsave("SpeciesRarefactionCurveAbundance_spobs.png",
       data_estaccumR_obs_plot)



## Species Rarefaction Curve for Shallow Reef at Cagayancillo ##
data_vegan_CAG_shallow <- 
  bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat_mpa == "Shallow Reef CAGAYANCILLO") %>%
  dplyr::select(colnames(data_vegan))

# data_vegan_CAG_shallow.env <- 
#   bind_cols(data_vegan, data_vegan.env) %>%
#   filter(habitat_mpa == "Shallow Reef CAGAYANCILLO") %>%
#   dplyr::select(-colnames(data_vegan))

p_CAG_shallow <- estaccumR(data_vegan_CAG_shallow, permutations = 999)

p_CAG_shallow_obs <-
  p_CAG_shallow$S %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))
p_CAG_shallow_plot <- p_CAG_shallow_obs %>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#f8bac4") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Shallow Reef at Cagayancillo") +
  ylim(0,50) +
  theme_classic() 

p_CAG_shallow_plot

##Old Plot
# p_CAG_shallow_plot <- plot(p_CAG_shallow,
#                            display = c("chao"),
#                            main = "Shallow Reef at Cagayancillo")
# p_CAG_shallow_plot

## Species Rarefaction Curve for Deep Reef at Cagayancillo ##
data_vegan_CAG_deep <- 
  bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat_mpa == "Deep Reef CAGAYANCILLO") %>%
  dplyr::select(colnames(data_vegan))

p_CAG_deep <- estaccumR(data_vegan_CAG_deep, permutations = 999)

p_CAG_deep_obs <-
  p_CAG_deep$S %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))
p_CAG_deep_plot <- p_CAG_deep_obs %>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#c5d8ea") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Mesophotic Reef at Cagayancillo") +
  ylim(0,50) +
  xlim(0,15)+
  theme_classic()

p_CAG_deep_plot

##Old Plot
# p_CAG_deep_plot <- plot(p_CAG_deep,
#                            display = c("chao"),
#                         main = "Mesophotic Reef at Cagayancillo")
# 
# p_CAG_deep_plot

## Species Rarefaction Curve for Shallow Reef at TRNP ##
data_vegan_TUB_shallow <- 
  bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat_mpa == "Shallow Reef TRNP") %>%
  dplyr::select(colnames(data_vegan))

p_TUB_shallow <- estaccumR(data_vegan_TUB_shallow, permutations = 999)

p_TUB_shallow_obs <-
  p_TUB_shallow$S %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))
p_TUB_shallow_plot <- p_TUB_shallow_obs %>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#f8bac4") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Shallow Reef at TRNP") +
  ylim(0,50) +
  xlim(0,15)+
  theme_classic()

p_TUB_shallow_plot

##Old Plot

# p_TUB_shallow_plot <- plot(p_TUB_shallow,
#                         display = c("chao"),
#                         main = "Shallow Reef at TRNP")
# p_TUB_shallow_plot

## Species Rarefaction Curve for Deep Reef at TRNP ##
data_vegan_TUB_deep <- 
  bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat_mpa == "Deep Reef TRNP") %>%
  dplyr::select(colnames(data_vegan))

p_TUB_deep <- estaccumR(data_vegan_TUB_deep, permutations = 999)
p_TUB_deep_obs <-
  p_TUB_deep$S%>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))
p_TUB_deep_plot <- p_TUB_deep_obs %>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#c5d8ea") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Mesophotic Reef at TRNP") +
  ylim(0,50) +
  xlim(0,20) +
  theme_classic()

p_TUB_deep_plot

##Old Plot
# p_TUB_deep_plot <- plot(p_TUB_deep,
#                            display = c("chao"),
#                         main = "Mesophotic Reef at TRNP")
# 
# p_TUB_deep_plot



##All Habitat and Study Locations Chao Plot ##
p_locations_depth_plot <- ggarrange(p_TUB_shallow_plot,
                                    p_TUB_deep_plot,
                                    p_CAG_shallow_plot,
                                    p_CAG_deep_plot,
                                    ncol = 2,
                                    nrow = 2)

ggsave("SpeciesRarefactionCurves_Study_Location_DepthCategory.png",
       p_locations_depth_plot,
       height = 14,
       width = 12,
       units = "in")

## Species Rarefaction Curve for Serranidae ##
data_vegan_Serranidae <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Serranidae"))


# data_vegan_Serranidae_nozeros <- data_vegan_Serranidae %>%
#                             select_if(colSums(.) != 0)
# 
# p_Serranidae_nz <- estaccumR(data_vegan_Serranidae_nozeros, permutations = 999)
# p_Serranidae_plot_nz <- plot(p_Serranidae_nz,
#                           display = c("chao"),
#                           main = "Serranidae")
# p_Serranidae_plot_nz



p_Serranidae <- estaccumR(data_vegan_Serranidae, permutations = 999)
View(p_Serranidae)

p_Serranidae_obs <-
  p_Serranidae$S %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>% 
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))
p_Serranidae_plot <- p_Serranidae_obs %>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#ebe8f3") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Serranidae") +
  ylim(0,20) +
  theme_classic()

p_Serranidae_plot

##Old Plots

# p_Serranidae_plot <- plot(p_Serranidae,
#                         display = c("chao"),
#                         main = "Serranidae")
# 
# p_Serranidae_plot

## Species Rarefaction Curve for Lutjanidae ##
data_vegan_Lutjanidae <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Lutjanidae"))

p_Lutjanidae <- estaccumR(data_vegan_Lutjanidae, permutations = 999)

p_Lutjanidae_obs <-
  p_Lutjanidae$S %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))
p_Lutjanidae_plot <- p_Lutjanidae_obs%>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#ebe8f3") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Lutjanidae") +
  ylim(0,20) +
  theme_classic()

p_Lutjanidae_plot


## Species Rarefaction Curve for Lethrinidae ##
data_vegan_Lethrinidae <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Lethrinidae"))

p_Lethrinidae <- estaccumR(data_vegan_Lethrinidae, permutations = 999)
p_Lethrinidae_obs <-
  p_Lethrinidae$S %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))
p_Lethrinidae_plot <- p_Lethrinidae_obs %>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#ebe8f3") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Lethrinidae") +
  ylim(0,20) +
  theme_classic()

p_Lethrinidae_plot

##Old Plots
# p_Lethrinidae_plot <- plot(p_Lethrinidae,
#                           display = c("chao"),
#                           main = "Lethrinidae")
# 
# p_Lethrinidae_plot

## Species Rarefaction Curve for Carangidae ##
data_vegan_Carangidae <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Carangidae"))

p_Carangidae <- estaccumR(data_vegan_Carangidae, permutations = 999)
p_Carangidae_obs <-
  p_Carangidae$S %>%
  # t() %>%
  as_tibble() %>%
  dplyr::mutate(N = row_number()) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "permutation") %>%
  group_by(N) %>%
  dplyr::summarize(S_mean = mean(value),
                   S_ci_lower = quantile(value,
                                         probs = 0.025),
                   S_ci_upper = quantile(value,
                                         probs = 0.975))
p_Carangidae_plot <- p_Carangidae_obs %>%
  ggplot(aes(x=N,
             y=S_mean)) +
  geom_ribbon(aes(ymin=S_ci_lower,
                  ymax=S_ci_upper),
              fill = "#ebe8f3") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Species Richness") +
  labs(title ="Carangidae") +
  ylim(0,20) +
  theme_classic()

p_Carangidae_plot

## Species Rarefaction Curves for All Groupings ##
p_groupings_obs_plot <- ggpubr::ggarrange(data_estaccumR_obs_plot,
                                          p_Serranidae_plot,
                                          p_Lutjanidae_plot,
                                          p_Lethrinidae_plot,
                                          p_Carangidae_plot,
                                          ncol = 2,
                                          nrow = 3)


ggsave("SpeciesRichnessRarefactionCurveGroupings_speciesobs.png", 
       p_groupings_obs_plot, 
       height = 14,
       width = 12,
       units = "in")


#### vegan::poolaccum Extrapolated Species Richness Curve in a Species Pool Based on Presence Absence ####

# incidence based (presence / absence) richness rarefaction curve
# creates 1 curve per data frame, so if you want multiple curves, have to make them separately then combine into 1 tibble to plot
# increase permutations to 999 if you use this for your project
p_presenceabsence<-poolaccum(data_vegan, permutations = 999)
p.plot_presenceabsence <-plot(p, display = c("chao"))
p.plot_presenceabsence
# autoplot(p_presenceabsence)
save_plot("SpeciesRichnessRarefactionCurvePresenceAbsence.png", 
          p.plot)

## Species Rarefaction Curve for Shallow Reef at Cagayancillo ##
data_vegan_CAG_shallow <- 
  bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat_mpa == "Shallow Reef CAGAYANCILLO") %>%
  dplyr::select(colnames(data_vegan))

# data_vegan_CAG_shallow.env <- 
#   bind_cols(data_vegan, data_vegan.env) %>%
#   filter(habitat_mpa == "Shallow Reef CAGAYANCILLO") %>%
#   dplyr::select(-colnames(data_vegan))

p_CAG_shallow_pa <- poolaccum(data_vegan_CAG_shallow, permutations = 999)

p_CAG_shallow_plot_pa <- plot(p_CAG_shallow_pa,
                              display = c("chao"),
                              main = "Shallow Reef at Cagayancillo")
p_CAG_shallow_plot_pa

## Species Rarefaction Curve for Deep Reef at Cagayancillo ##
data_vegan_CAG_deep <- 
  bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat_mpa == "Deep Reef CAGAYANCILLO") %>%
  dplyr::select(colnames(data_vegan))

p_CAG_deep_pa <- poolaccum(data_vegan_CAG_deep, permutations = 999)

p_CAG_deep_plot_pa <- plot(p_CAG_deep_pa,
                           display = c("chao"),
                           main = "Mesophotic Reef at Cagayancillo")

p_CAG_deep_plot_pa

## Species Rarefaction Curve for Shallow Reef at TRNP ##
data_vegan_TUB_shallow <- 
  bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat_mpa == "Shallow Reef TRNP") %>%
  dplyr::select(colnames(data_vegan))

p_TUB_shallow_pa <- poolaccum(data_vegan_TUB_shallow, permutations = 999)

p_TUB_shallow_plot_pa <- plot(p_TUB_shallow_pa,
                              display = c("chao"),
                              main = "Shallow Reef at TRNP")
p_TUB_shallow_plot_pa

## Species Rarefaction Curve for Deep Reef at TRNP ##
data_vegan_TUB_deep <- 
  bind_cols(data_vegan, data_vegan.env) %>%
  filter(habitat_mpa == "Deep Reef TRNP") %>%
  dplyr::select(colnames(data_vegan))

p_TUB_deep_pa <- poolaccum(data_vegan_TUB_deep, permutations = 999)

p_TUB_deep_plot_pa <- plot(p_TUB_deep_pa,
                           display = c("chao"),
                           main = "Mesophotic Reef at TRNP")

p_TUB_deep_plot_pa

##All Habitat and Study Locations Chao Plot Presence Absence ##
p_habitat_locations_plot_pa <- ggarrange(p_CAG_shallow_plot_pa,
                                         p_CAG_deep_plot_pa,
                                         p_TUB_shallow_plot_pa,
                                         p_TUB_deep_plot_pa,
                                         ncol = 2,
                                         nrow = 2)
ggsave("SpeciesRichnessRarefactionCurveHabitatandStudyLocationsPresenceAbsence.pdf", 
       p_habitat_locations_plot_pa, 
       height = 11,
       width = 8.5,
       units = "in")

## Species Rarefaction Curve for Serranidae ##
data_vegan_Serranidae <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Serranidae"))

p_Serranidae_pa <- poolaccum(data_vegan_Serranidae, permutations = 999)

p_Serranidae_plot_pa <- plot(p_Serranidae_pa,
                             display = c("chao"),
                             main = "Serranidae")

p_Serranidae_plot_pa

## Species Rarefaction Curve for Lutjanidae ##
data_vegan_Lutjanidae <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Lutjanidae"))

p_Lutjanidae_pa <- poolaccum(data_vegan_Lutjanidae, permutations = 999)

p_Lutjanidae_plot_pa <- plot(p_Lutjanidae_pa,
                             display = c("chao"),
                             main = "Lutjanidae")

p_Lutjanidae_plot_pa

## Species Rarefaction Curve for Lethrinidae ##
data_vegan_Lethrinidae <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Lethrinidae"))

p_Lethrinidae_pa <- poolaccum(data_vegan_Lethrinidae, permutations = 999)

p_Lethrinidae_plot_pa <- plot(p_Lethrinidae_pa,
                              display = c("chao"),
                              main = "Lethrinidae")

p_Lethrinidae_plot_pa

## Species Rarefaction Curve for Carangidae ##
data_vegan_Carangidae <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Carangidae"))


p_Carangidae_pa <- poolaccum(data_vegan_Carangidae, permutations = 999)

p_Carangidae_plot_pa <- plot(p_Carangidae_pa,
                             display = c("chao"),
                             main = "Carangidae")

p_Carangidae_plot_pa


# ## Species Rarefaction Curve for Galeomorphii ##
# data_vegan_Galeomorphii <- 
#   data_removed_sp %>%
#   dplyr::select(op_code,
#                 taxon,
#                 max_n) %>%
#   # convert tibble from long to wide format
#   pivot_wider(names_from = taxon,
#               values_from = max_n,
#               values_fill = 0) %>%
#   # sort by op_code
#   arrange(op_code) %>%
#   # remove the op_code column for vegan
#   dplyr::select(-op_code) %>%
#   dplyr::select(contains("Galeomorphii"))

# view(data_vegan_Galeomorphii)

# p_Galeomorphii_pa <- poolaccum(data_vegan_Galeomorphii, permutations = 999)
# 
# p_Galeomorphii_plot_pa <- plot(p_Galeomorphii_pa,
#                                display = c("chao"),
#                                main = "Galeomorphii")
# 
# p_Galeomorphii_plot_pa

##Species Rarefaction Curve for Cheilinus undulatus ##
data_vegan_Cheilinus_undulatus <- 
  data_removed_sp %>%
  dplyr::select(op_code,
                taxon,
                max_n) %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code) %>%
  dplyr::select(contains("Cheilinus undulatus"))
# select_if(colSums(.) != 0)


View(data_vegan_Cheilinus_undulatus)

p_Cheilinus_undulatus_pa <- poolaccum(data_vegan_Cheilinus_undulatus, permutations = 999) 


p_Cheilinus_undulatus_plot_pa <- plot(p_Cheilinus_undulatus_pa,
                                      display = c("chao"),
                                      main = "Cheilinus undulatus")

p_Cheilinus_undulatus_plot_pa

## Species Rarefaction Curves for All Groupings ##
p_groupings_plot_pa <- ggarrange(p_Serranidae_plot_pa,
                                 p_Lutjanidae_plot_pa,
                                 p_Lethrinidae_plot_pa,
                                 p_Carangidae_plot_pa,
                                 p_Galeomorphii_plot_pa,
                                 # p_Cheilinus_undulatus_plot,
                                 ncol = 2,
                                 nrow = 3)
ggsave("SpeciesRichnessRarefactionCurveGroupingsPresenceAbsence.pdf", 
       p_groupings_plot_pa, 
       height = 11,
       width = 8.5,
       units = "in")

# data_vegan_CAG_shallow <- 
# bind_cols(data_vegan, data_vegan.env) %>%
#   filter(habitat_mpa == "Shallow Reef CAGAYANCILLO") %>%
#   dplyr::select(colnames(data_vegan))
# 
# data_vegan_CAG_shallow.env <- 
#   bind_cols(data_vegan, data_vegan.env) %>%
#   filter(habitat_mpa == "Shallow Reef CAGAYANCILLO") %>%
#   dplyr::select(-colnames(data_vegan))
# 
# p_CAG_shallow <- poolaccum(data_vegan_CAG_shallow, permutations = 9999)
# autoplot(p_CAG_shallow)


