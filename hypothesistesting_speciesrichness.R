#### INITIALIZATION ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(janitor)
library(magrittr)
library(gridExtra)
library(readxl)

library(devtools)

library(vegan)

library(fitdistrplus)
library(emmeans)
library(multcomp)
library(multcompView)
library(ggeffects)

library(rlang)
library(afex)
library(ggbeeswarm)
library(performance)
library(optimx)
library(effects)
library(prediction)
library(ggforce)

library(readr)
library(ggpubr)
library(sjPlot)

#### PACKAGES ####
packages_used <- 
  c("tidyverse", # have to install from github, messes up tidyverse code, have to add dplyr:: to several commands
    "janitor",
    "magrittr",
    "gridExtra",
    "readxl",
    "devtools",
    "vegan",
    "fitdistrplus",
    "emmeans",
    "multcomp",
    "multcompView",
    "rlang",
    "afex",
    "ggbeeswarm",
    "performance",
    "optimx",
    "effects",
    "prediction",
    "ggforce",
    "readr",
    "ggpubr",
    "sjPlot")

# NOTE: after loading these packages, you may find that tidyverse commands are affected 
#       the solution is to add the appropriate package name before commands that break code
#       such as `dplyr::select` if `select` doesn't work correctly anymore
#       this happens when multiple packages have the same command names. 
#       The last package loaded takes precidence, and your tidyverse commands break.
#       you could load tidyverse last, but invariably, you will load a package after tidyverse
#       so it's impossible to avoid this

packages_to_install <- 
  packages_used[!packages_used %in% installed.packages()[,1]]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install, 
                   Ncpus = Sys.getenv("NUMBER_OF_PROCESSORS") - 1)
}

lapply(packages_used, 
       require, 
       character.only = TRUE)


#### USER DEFINED VARIABLES ####
inFilePath2 = "./meso_euphotic_carniv_fish_videobaitstations_all.rds"
functionPath = "./model_fitting_functions.R"
inFilePath3 = "./WorkingData_CLEANED_TUB,CAG.xlsx"
inFilePath4 = "./PHIRES_MetaData.xlsx"

#### READ IN DATA ####

# read in data and remove rows with missing data or multiple bands
data_all <-
  read_rds(inFilePath2) 


data <-
  read_excel(inFilePath3,
             na="NA") %>%
  janitor::clean_names() %>%
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

data_removed_sp <- data %>%
  filter(species != "sp") %>%
  mutate(family = str_to_title(family),
         genus = str_to_title(genus),
         species = str_to_lower(species),
         family_clean = case_when(
           family == "Epinephelidae" ~ "Serranidae",
           TRUE ~ family)) %>%
  mutate(groupings = case_when(
    family == "Labridae" ~ "Cheilinus undulatus",
    family == "Epinephelidae" ~ "Serranidae",
    TRUE ~ family))%>%
  mutate(taxon = str_c(groupings,
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
           .keep_all = TRUE) %>% view()


metadata <-
  read_excel(inFilePath4,
             na="NA") %>%
  clean_names() %>%
  dplyr::rename(bait_weight_grams = weight_grams)




####COMBINE DATA ####
data_all_removed_sp <- 
  data_removed_sp %>%
  left_join(metadata,
            by = c("op_code" = "opcode",
                   "depth_m" = "depth_m")) %>%
  # rearrange order of columns, metadata then data
  # select(op_code,
  #        site:long_e,
  #        depth_m,
  #        time_in:bait_weight_grams,
  #        everything()) %>%
  mutate(study_locations = case_when(
    site == "Cawili" ~ "CAGAYANCILLO",
    site == "Calusa" ~ "CAGAYANCILLO",
    site == "Cagayancillo" ~ "CAGAYANCILLO",
    site == "TUBBATAHA" ~ "TRNP")) %>%
  mutate(study_locations = factor(study_locations,
                                  levels = c("TRNP", 
                                             "CAGAYANCILLO"))) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Deep Reef"))) %>%
  mutate(habitat = case_when(
    habitat == "Shallow Reef" ~ "Shallow Reef",
    habitat == "Deep Reef" ~ "Mesophotic Reef"
  )) %>%
  mutate(groupings = case_when(
    family == "Labridae" ~ "Cheilinus undulatus",
    family == "Epinephelidae" ~ "Serranidae",
    TRUE ~ family))

#### Prep Data for Vegan ####

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
         study_locations = factor(study_locations,
                                  levels = c("TRNP",
                                             "CAGAYANCILLO")),
         habitat = factor(habitat),
         bait_type = factor(bait_type),
         site = factor(site),
         survey_area = factor(survey_area),
         habitat_mpa = str_c(habitat,
                             study_locations,
                             sep = " "))

attach(data_vegan.env)
#### Overall mean_chao_s: Make Visualization of Data ####
pool <- 
  estimateR(x = data_vegan) %>%
  t() %>%
  as_tibble()

# vis_dists(pool,
#           "S.chao1")
# par(mar = c(1,1,1,1))

pool %>%
  clean_names() %>%
  bind_cols(data_vegan.env) %>%
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
                                             "CAGAYANCILLO")))%>%
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
  # scale_fill_manual(values = habitatcolors,
  #                   labels = c("Shallow Reef",
  #                              "Mesophotic Reef")) +
  labs(title = "Species Richness at Cagayancillo vs. Tubbataha", fill = "Habitat") +
  xlab("Study Locations") +
  ylab("Mean Chao Estimate of Species Richness")
#### Overall mean_chao_s: Mixed Effects Hypothesis Test ####
pool <- 
  estimateR(x = data_vegan) %>%
  t() %>%
  as_tibble()

data_chao_s <- 
  pool %>%
  clean_names() %>%
  bind_cols(data_vegan.env) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef")))



## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(s_chao1) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "Gamma"


alpha_sig = 0.05


# we start with the loci subjected to 11 primer concentrations (we removed loci with no sum_max_n to simplify)


sampling_design = "s_chao1 ~  habitat * study_locations + (1|study_locations:bait_type)"


# # fit mixed model
model <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_chao_s)

model
anova(model)

# visualize summary(model)
emmip(model, 
      study_locations ~ habitat,    # type = "response" for back transformed values
      cov.reduce = range) +
  # geom_vline(xintercept=mean(data_all_summaxn_$primer_x),
  #            color = "grey",
  #            linetype = "dashed") +
  # geom_text(aes(x = mean(data_all_summaxn_$primer_x),
  #               y = -2,
  #               label = "mean primer_x"),
  #           color = "grey") +
  theme_classic() +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")
### mean_chao_s: Conduct A priori contrast tests for differences among sites ####
emmeans_model_sr <<-
  emmeans(model,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_sr,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_sr), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")
#### mean_chao_s: Group Sites Based on Model Results ####
groupings_model_sr <<-
  multcomp::cld(emmeans_model_sr, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_sr             # these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_sr <<-
  summary(emmeans_model_sr,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_sr %>%
              dplyr::select(-response:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_sr <- groupings_model_fixed_sr %>%
  mutate(habitat = factor(habitat,
                          levels = c(
                            "Shallow Reef",
                            "Mesophotic Reef")))

habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Mesophotic Reef")

#### mean_chao_s: Visualize Estimated Marginal Means Output with Group Categories ####
p_sr <- 
  groupings_model_fixed_sr %>%
  ggplot(aes(x=study_locations,
             y=response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  # geom_point(data = data_chao_s,
  #            aes(x = study_locations,
  #                y = !!response_var,
  #                shape = habitat
  #            ),
  #            position = position_jitterdodge(),
  #            color = "grey50",
  #            # shape = 1,
  #            size = 3) +
geom_point(data = data_chao_s,
           aes(x = study_locations,
               y = !!response_var,
               shape = habitat
           ),
           position = position_jitterdodge(
             jitter.width = 0.6,
             dodge.width = 0.7
             # jitter.height = 0.05
           ),
           color = "grey50",
           # shape = 1,
           size = 3) +
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "black",
                size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  # geom_text(aes(label=group),
  #           position = position_dodge(width=0.9),
  #           vjust = -0.5,
  #           hjust = -0.15,
  #           size = 8 / (14/5)) +  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  theme_classic() +
  labs(x = "Study Locations",
       y = "Species Richness (Chao1)",
       title = "Overall Species Richness") +
  # ylim(ymin, 
  #      ymax) +
  # labs(title = "Species Richness at TRNP vs. Cagayancillo",
  #      subtitle = "Distribution Family = Gamma",
  #       x = "Study Locations",
  #      y = "EM Means of Chao Estimate of Species Richness") +
  theme(legend.position=c(0.5,0.8),  
        legend.title=element_blank()) +
  ylim(0,40)+
  scale_fill_manual(values = habitatcolors,
                    labels = c("Shallow",
                               "Mesophotic"))

p_sr
ggsave("EMMeansofSpeciesRichnessGamma.png",
       p_sr)


#### Overall means_chao_s with Species Observations instead of Chao estimate ####
pool <- 
  estimateR(x = data_vegan) %>%
  t() %>%
  as_tibble()

data_chao_s <- 
  pool %>%
  clean_names() %>%
  bind_cols(data_vegan.env)


## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(s_obs) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"


alpha_sig = 0.05


# we start with the loci subjected to 11 primer concentrations (we removed loci with no sum_max_n to simplify)

sampling_design = "s_obs ~  habitat * study_locations"
#fit glm model
model <<- 
  glm(formula = sampling_design, 
      family = distribution_family,
      data = data_chao_s)
# sampling_design = "s_chao1 ~  habitat * study_locations + (1|study_locations:bait_type)"


# # fit mixed model
# model <<-
#   afex::mixed(formula = sampling_design,
#               family = distribution_family,
#               method = "LRT",
#               sig_symbols = rep("", 4),
#               # all_fit = TRUE,
#               data = data_chao_s)

model
anova(model)

# visualize summary(model)
emmip(model, 
      study_locations ~ habitat,    # type = "response" for back transformed values
      cov.reduce = range) +
  # geom_vline(xintercept=mean(data_all_summaxn_$primer_x),
  #            color = "grey",
  #            linetype = "dashed") +
  # geom_text(aes(x = mean(data_all_summaxn_$primer_x),
  #               y = -2,
  #               label = "mean primer_x"),
  #           color = "grey") +
  theme_classic() +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")
## mean_chao_s: Conduct A priori contrast tests for differences among sites ##
emmeans_model_sr <<-
  emmeans(model,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_sr,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_sr), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")
## mean_chao_s: Group Sites Based on Model Results ##
groupings_model_sr <<-
  multcomp::cld(emmeans_model_sr, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_sr             # these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_sr <<-
  summary(emmeans_model_sr,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_sr %>%
              dplyr::select(-response:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_sr <- groupings_model_fixed_sr %>%
  mutate(habitat = factor(habitat,
                          levels = c(
                            "Shallow Reef",
                            "Mesophotic Reef")))

habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Mesophotic Reef")

## mean_chao_s: Visualize Estimated Marginal Means Output with Group Categories ##
p_sr <- 
  groupings_model_fixed_sr %>%
  ggplot(aes(x=study_locations,
             y=response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  geom_point(data = data_chao_s,
             aes(x = study_locations,
                 y = !!response_var
             ),
             position = position_jitterdodge(),
             # color = "grey70",
             # shape = 1,
             size = 1) +
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "grey50",
                # size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  geom_text(aes(label=group),
            position = position_dodge(width=0.9),
            vjust = -0.5,
            hjust = -0.15,
            size = 8 / (14/5)) +  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  theme_classic() +
  # ylim(ymin, 
  #      ymax) +
  labs(title = "Species Observations at TRNP vs. Cagayancillo",
       subtitle = "Distribution Family = Poisson",
       x = "Study Locations",
       y = "Estimated Marginal Means of Species Richness") +
  theme(legend.position=c(0.6,0.8),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors)

p_sr
save_plot("EMMeansofSpeciesRichnessPoisson.png")


#### mean_chao_s: Serranidae ####
## Make New Data Vegan for Serranidae
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
         study_locations = factor(study_locations,
                                  levels = c("TRNP",
                                             "CAGAYANCILLO")),
         habitat = factor(habitat),
         bait_type = factor(bait_type),
         site = factor(site),
         survey_area = factor(survey_area),
         habitat_mpa = str_c(habitat,
                             study_locations,
                             sep = " "))

attach(data_vegan.env)

pool_Serranidae <- 
  estimateR(x = data_vegan_Serranidae) %>%
  t() %>%
  as_tibble()

#Visualize statistical distribution
# source(functionPath)
# vis_dists(pool_Serranidae,
#           "S.chao1")
# vis_dists(data_chao_s_Serranidae,
#           "s_chao1")

data_chao_s_Serranidae <- pool_Serranidae %>%
  clean_names() %>%
  bind_cols(data_vegan.env) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef")))

# data_chao_s_Serranidae %>%
#   ggplot(aes(x = s_chao1)) +
#   geom_histogram()

## Enter Information About Your Data for A Hypothesis Test ##


# define your response variable, here it is binomial
response_var = quo(s_obs) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"

alpha_sig = 0.05


## Histogram and Visualizing Distance Matrix
# data_chao_s_Serranidae %>%
#   ggplot(aes(x = s_obs)) + 
#   geom_histogram() +
#   facet_grid(habitat ~ study_locations)
# 
# data_chao_s_Serranidae %>%
#   ggplot(aes(x = s_chao1)) + 
#   geom_histogram() +
#   facet_grid(habitat ~ study_locations)
# 
# vis_dists(data_chao_s_Serranidae,
#           "s_chao1")
# 
# vis_dists(data_chao_s_Serranidae,
#           "s_obs")


sampling_design = "s_obs ~  habitat * study_locations"
#fit glm model
model_Serranidae <<- 
  glm(formula = sampling_design, 
      family = distribution_family,
      data = data_chao_s_Serranidae)


# sampling_design = "s_obs ~  habitat * study_locations + (1|study_locations:bait_type)"


# # # fit mixed model
# model_Serranidae <<-
#   afex::mixed(formula = sampling_design,
#               family = distribution_family,
#               method = "LRT",
#               sig_symbols = rep("", 4),
#               # all_fit = TRUE,
#               data = data_chao_s_Serranidae)

model_Serranidae
anova(model_Serranidae)

# visualize summary(model)
emmip(model_Serranidae, 
      study_locations ~ habitat,    # type = "response" for back transformed values
      cov.reduce = range) +
  # geom_vline(xintercept=mean(data_all_summaxn_$primer_x),
  #            color = "grey",
  #            linetype = "dashed") +
  # geom_text(aes(x = mean(data_all_summaxn_$primer_x),
  #               y = -2,
  #               label = "mean primer_x"),
  #           color = "grey") +
  theme_classic() +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")

## mean_chao_s Serranidae: Conduct A priori contrast tests for differences among sites ##
emmeans_model_sr_Serranidae <<-
  emmeans(model_Serranidae,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_sr_Serranidae,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_sr_Serranidae), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")

## mean_chao_s: Group Sites Based on Model Results ##
groupings_model_sr_Serranidae <<-
  multcomp::cld(emmeans_model_sr_Serranidae, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_sr_Serranidae  
# these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_sr_Serranidae <<-
  summary(emmeans_model_sr_Serranidae,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_sr_Serranidae %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_sr_Serranidae <- groupings_model_fixed_sr_Serranidae %>%
  mutate(habitat = factor(habitat,
                          levels = c(
                            "Shallow Reef",
                            "Mesophotic Reef")))

habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Mesophotic Reef")


## mean_chao_s Serranidae: Visualize Estimated Marginal Means Output with Group Categories ##
p_sr_Serranidae <- 
  groupings_model_fixed_sr_Serranidae %>%
  ggplot(aes(x=study_locations,
             y= response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  geom_point(data = data_chao_s_Serranidae,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat
             ),
             position = position_jitterdodge(
               jitter.width = 0.651,
               dodge.width = 0.7
               # jitter.height = 0.05
             ),
             color = "grey50",
             # shape = 1,
             size = 3) +
  # geom_beeswarm(data = data_chao_s_Serranidae,
  #            aes(x = study_locations,
  #                y = !!response_var,
  #                shape = habitat
  #            ),
  #            # position = position_jitterdodge(
  #            #   jitter.width = 0.3,
  #            #   dodge.width = 0.9
  #            #   # jitter.height = 0.05
  #            # ),
  #            color = "grey50",
  #            # shape = 1,
  #            size = 3) +
# geom_quasirandom(data = data_chao_s_Serranidae,
#            aes(x = study_locations,
#                y = !!response_var,
#                shape = habitat
#            ),
#            position = position_jitterdodge(
#              jitter.width = 0.3,
#              dodge.width = 0.9
#              # jitter.height = 0.05
#            ),
#            color = "grey50",
#            # shape = 1,
#            size = 3) +
  
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "black",
                size = 1,
                position = position_dodge(width=.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  # geom_text(aes(label=group),
  #           position = position_dodge(width=0.9),
  #           vjust = -0.5,
  #           hjust = -0.15,
  #           size = 8 / (14/5)) +  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  theme_classic() +
  labs(x = "Study Locations",
       y = "Observed Species Richness",
       title = "Serranidae") +
  ylim(0, 5) +
  # labs(title = "Serranidae",
  #      subtitle = "Distribution Family = Poisson",
  #     x = "Study Locations",
  #      y = "EM Means of Species Richness") +
  theme(legend.position=c(0.5,0.8),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors,
                    labels = c("Shallow",
                               "Mesophotic"))

p_sr_Serranidae
  

#### mean_chao_s of Lutjanidae ####
## Make New Data Vegan for Lutjanidae
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
         study_locations = factor(study_locations,
                                  levels = c("TRNP",
                                             "CAGAYANCILLO")),
         habitat = factor(habitat),
         bait_type = factor(bait_type),
         site = factor(site),
         survey_area = factor(survey_area),
         habitat_mpa = str_c(habitat,
                             study_locations,
                             sep = " "))

attach(data_vegan.env)

pool_Lutjanidae <- 
  estimateR(x = data_vegan_Lutjanidae) %>%
  t() %>%
  as_tibble()


data_chao_s_Lutjanidae <- pool_Lutjanidae %>%
  clean_names() %>%
  bind_cols(data_vegan.env) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef")))

# data_chao_s_Lutjanidae %>%
#   ggplot(aes(x = s_chao1)) +
#   geom_histogram()

## Enter Information About Your Data for A Hypothesis Test ##


# define your response variable, here it is binomial
response_var = quo(s_obs) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"

alpha_sig = 0.05


## Histogram and Visualizing Distance Matrix
# data_chao_s_Lutjanidae %>%
#   ggplot(aes(x = s_obs)) + 
#   geom_histogram() +
#   facet_grid(habitat ~ study_locations)
# 
# data_chao_s_Lutjanidae %>%
#   ggplot(aes(x = s_chao1)) + 
#   geom_histogram() +
#   facet_grid(habitat ~ study_locations)
# 
# vis_dists(data_chao_s_Lutjanidae,
#           "s_chao1")
# 
# vis_dists(data_chao_s_Lutjanidae,
#           "s_obs")


sampling_design = "s_obs ~  habitat * study_locations"
#fit glm model
model_Lutjanidae <<- 
  glm(formula = sampling_design, 
      family = distribution_family,
      data = data_chao_s_Lutjanidae)


# sampling_design = "s_obs ~  habitat * study_locations + (1|study_locations:bait_type)"


# # # fit mixed model
# model_Lutjanidae <<-
#   afex::mixed(formula = sampling_design,
#               family = distribution_family,
#               method = "LRT",
#               sig_symbols = rep("", 4),
#               # all_fit = TRUE,
#               data = data_chao_s_Lutjanidae)

model_Lutjanidae
anova(model_Lutjanidae)

# visualize summary(model)
emmip(model_Lutjanidae, 
      study_locations ~ habitat,    # type = "response" for back transformed values
      cov.reduce = range) +
  # geom_vline(xintercept=mean(data_all_summaxn_$primer_x),
  #            color = "grey",
  #            linetype = "dashed") +
  # geom_text(aes(x = mean(data_all_summaxn_$primer_x),
  #               y = -2,
  #               label = "mean primer_x"),
  #           color = "grey") +
  theme_classic() +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")

## mean_chao_s Lutjanidae: Conduct A priori contrast tests for differences among sites ##
emmeans_model_sr_Lutjanidae <<-
  emmeans(model_Lutjanidae,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_sr_Lutjanidae,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_sr_Lutjanidae), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")

## mean_chao_s: Group Sites Based on Model Results ##
groupings_model_sr_Lutjanidae <<-
  multcomp::cld(emmeans_model_sr_Lutjanidae, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_sr_Lutjanidae  
# these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_sr_Lutjanidae <<-
  summary(emmeans_model_sr_Lutjanidae,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_sr_Lutjanidae %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_sr_Lutjanidae <- groupings_model_fixed_sr_Lutjanidae %>%
  mutate(habitat = factor(habitat,
                          levels = c(
                            "Shallow Reef",
                            "Mesophotic Reef")))

habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Mesophotic Reef")


## mean_chao_s Lutjanidae: Visualize Estimated Marginal Means Output with Group Categories ##
p_sr_Lutjanidae <- 
  groupings_model_fixed_sr_Lutjanidae %>%
  ggplot(aes(x=study_locations,
             y= response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  # geom_point(data = data_chao_s_Lutjanidae,
  #            aes(x = study_locations,
  #                y = !!response_var,
  #                shape = habitat
  #            ),
  #            position = position_jitterdodge(),
  #            color = "grey50",
  #            # shape = 1,
  #            size = 3) +
  geom_point(data = data_chao_s_Lutjanidae,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat
             ),
             position = position_jitterdodge(
               jitter.width = 0.65,
               dodge.width = 0.7
               # jitter.height = 0.05
             ),
             color = "grey50",
             # shape = 1,
             size = 3) +
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "black",
                size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  # geom_text(aes(label=group),
  # position = position_dodge(width=0.9),
  # vjust = -0.5,
  # hjust = -0.15,
  # size = 8 / (14/5)) +  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  theme_classic() +
  labs(x = "Study Locations",
       y = "Observed Species Richness",
       title = "Lutjanidae") +
  ylim(0,5) +
  # labs(title = "Lutjanidae",
  #      subtitle = "Distribution Family = Poisson",
  #      x = "Study Locations",
  #      y = "EM Means of Species Richness") +
  theme(legend.position=c(0.5,0.8),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors,
                    labels = c("Shallow",
                               "Mesophotic"))

p_sr_Lutjanidae
#### mean_chao_s of Lethrinidae ####
## Make New Data Vegan for Lethrinidae
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
         study_locations = factor(study_locations,
                                  levels = c("TRNP",
                                             "CAGAYANCILLO")),
         habitat = factor(habitat),
         bait_type = factor(bait_type),
         site = factor(site),
         survey_area = factor(survey_area),
         habitat_mpa = str_c(habitat,
                             study_locations,
                             sep = " "))

attach(data_vegan.env)

pool_Lethrinidae <- 
  estimateR(x = data_vegan_Lethrinidae) %>%
  t() %>%
  as_tibble()


data_chao_s_Lethrinidae <- pool_Lethrinidae %>%
  clean_names() %>%
  bind_cols(data_vegan.env) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef")))

# data_chao_s_Lethrinidae %>%
#   ggplot(aes(x = s_chao1)) +
#   geom_histogram()

## Enter Information About Your Data for A Hypothesis Test ##


# define your response variable, here it is binomial
response_var = quo(s_obs) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"

alpha_sig = 0.05


## Histogram and Visualizing Distance Matrix
# data_chao_s_Lethrinidae %>%
#   ggplot(aes(x = s_obs)) + 
#   geom_histogram() +
#   facet_grid(habitat ~ study_locations)
# 
# data_chao_s_Lethrinidae %>%
#   ggplot(aes(x = s_chao1)) + 
#   geom_histogram() +
#   facet_grid(habitat ~ study_locations)
# 
# vis_dists(data_chao_s_Lethrinidae,
#           "s_chao1")
# 
# vis_dists(data_chao_s_Lethrinidae,
#           "s_obs")


sampling_design = "s_obs ~  habitat * study_locations"
#fit glm model
model_Lethrinidae <<- 
  glm(formula = sampling_design, 
      family = distribution_family,
      data = data_chao_s_Lethrinidae)


# sampling_design = "s_obs ~  habitat * study_locations + (1|study_locations:bait_type)"


# # # fit mixed model
# model_Lethrinidae <<-
#   afex::mixed(formula = sampling_design,
#               family = distribution_family,
#               method = "LRT",
#               sig_symbols = rep("", 4),
#               # all_fit = TRUE,
#               data = data_chao_s_Lethrinidae)

model_Lethrinidae
anova(model_Lethrinidae)

# visualize summary(model)
emmip(model_Lethrinidae, 
      study_locations ~ habitat,    # type = "response" for back transformed values
      cov.reduce = range) +
  # geom_vline(xintercept=mean(data_all_summaxn_$primer_x),
  #            color = "grey",
  #            linetype = "dashed") +
  # geom_text(aes(x = mean(data_all_summaxn_$primer_x),
  #               y = -2,
  #               label = "mean primer_x"),
  #           color = "grey") +
  theme_classic() +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")

## mean_chao_s Lethrinidae: Conduct A priori contrast tests for differences among sites ##
emmeans_model_sr_Lethrinidae <<-
  emmeans(model_Lethrinidae,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_sr_Lethrinidae,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_sr_Lethrinidae), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")

## mean_chao_s: Group Sites Based on Model Results ##
groupings_model_sr_Lethrinidae <<-
  multcomp::cld(emmeans_model_sr_Lethrinidae, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_sr_Lethrinidae  
# these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_sr_Lethrinidae <<-
  summary(emmeans_model_sr_Lethrinidae,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_sr_Lethrinidae %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_sr_Lethrinidae <- groupings_model_fixed_sr_Lethrinidae %>%
  mutate(habitat = factor(habitat,
                          levels = c(
                            "Shallow Reef",
                            "Mesophotic Reef")))

habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Mesophotic Reef")


## mean_chao_s Lethrinidae: Visualize Estimated Marginal Means Output with Group Categories ##
p_sr_Lethrinidae <- 
  groupings_model_fixed_sr_Lethrinidae %>%
  ggplot(aes(x=study_locations,
             y= response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  # geom_point(data = data_chao_s_Lethrinidae,
  #            aes(x = study_locations,
  #                y = !!response_var,
  #                shape = habitat
  #            ),
  #            position = position_jitterdodge(),
  #            color = "grey50",
  #            # shape = 1,
  #            size = 3) +
geom_point(data = data_chao_s_Lethrinidae,
           aes(x = study_locations,
               y = !!response_var,
               shape = habitat
           ),
           position = position_jitterdodge(
             jitter.width = 0.6,
             dodge.width = 0.7
             # jitter.height = 0.05
           ),
           color = "grey50",
           # shape = 1,
           size = 3) +
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "black",
                size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  # geom_text(aes(label=group),
  #           position = position_dodge(width=0.9),
  #           vjust = -0.5,
  #           hjust = -0.15,
  #           size = 8 / (14/5)) +  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  theme_classic() +
  labs(x = "Study Locations",
       y = "Observed Species Richness",
       title = "Lethrinidae") +
  ylim(0,5) +
  # labs(title = "Lethrinidae",
  #      subtitle = "Distribution Family = Poisson",
  #      x = "Study Locations",
  #      y = "EM Means of Species Richness") +
  theme(legend.position=c(0.5,0.8),
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors,
                    labels = c("Shallow",
                               "Mesophotic"))

p_sr_Lethrinidae
#### mean_chao_s of Carangidae ####
## Make New Data Vegan for Carangidae
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
         study_locations = factor(study_locations,
                                  levels = c("TRNP",
                                             "CAGAYANCILLO")),
         habitat = factor(habitat),
         bait_type = factor(bait_type),
         site = factor(site),
         survey_area = factor(survey_area),
         habitat_mpa = str_c(habitat,
                             study_locations,
                             sep = " "))

attach(data_vegan.env)

pool_Carangidae <- 
  estimateR(x = data_vegan_Carangidae) %>%
  t() %>%
  as_tibble()


data_chao_s_Carangidae <- pool_Carangidae %>%
  clean_names() %>%
  bind_cols(data_vegan.env) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef")))

# data_chao_s_Carangidae %>%
#   ggplot(aes(x = s_chao1)) +
#   geom_histogram()

## Enter Information About Your Data for A Hypothesis Test ##


# define your response variable, here it is binomial
response_var = quo(s_obs) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"

alpha_sig = 0.05


## Histogram and Visualizing Distance Matrix
# data_chao_s_Carangidae %>%
#   ggplot(aes(x = s_obs)) + 
#   geom_histogram() +
#   facet_grid(habitat ~ study_locations)
# 
# data_chao_s_Carangidae %>%
#   ggplot(aes(x = s_chao1)) + 
#   geom_histogram() +
#   facet_grid(habitat ~ study_locations)
# 
# vis_dists(data_chao_s_Carangidae,
#           "s_chao1")
# 
# vis_dists(data_chao_s_Carangidae,
#           "s_obs")


sampling_design = "s_obs ~  habitat * study_locations"
#fit glm model
model_Carangidae <<- 
  glm(formula = sampling_design, 
      family = distribution_family,
      data = data_chao_s_Carangidae)


# sampling_design = "s_obs ~  habitat * study_locations + (1|study_locations:bait_type)"


# # # fit mixed model
# model_Carangidae <<-
#   afex::mixed(formula = sampling_design,
#               family = distribution_family,
#               method = "LRT",
#               sig_symbols = rep("", 4),
#               # all_fit = TRUE,
#               data = data_chao_s_Carangidae)

model_Carangidae
anova(model_Carangidae)

# visualize summary(model)
emmip(model_Carangidae, 
      study_locations ~ habitat,    # type = "response" for back transformed values
      cov.reduce = range) +
  # geom_vline(xintercept=mean(data_all_summaxn_$primer_x),
  #            color = "grey",
  #            linetype = "dashed") +
  # geom_text(aes(x = mean(data_all_summaxn_$primer_x),
  #               y = -2,
  #               label = "mean primer_x"),
  #           color = "grey") +
  theme_classic() +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")

## mean_chao_s Carangidae: Conduct A priori contrast tests for differences among sites ##
emmeans_model_sr_Carangidae <<-
  emmeans(model_Carangidae,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_sr_Carangidae,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_sr_Carangidae), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")

## mean_chao_s: Group Sites Based on Model Results ##
groupings_model_sr_Carangidae <<-
  multcomp::cld(emmeans_model_sr_Carangidae, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_sr_Carangidae  
# these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_sr_Carangidae <<-
  summary(emmeans_model_sr_Carangidae,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_sr_Carangidae %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_sr_Carangidae <- groupings_model_fixed_sr_Carangidae %>%
  mutate(habitat = factor(habitat,
                          levels = c(
                            "Shallow Reef",
                            "Mesophotic Reef")))

habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Mesophotic Reef")


## mean_chao_s Carangidae: Visualize Estimated Marginal Means Output with Group Categories ##
p_sr_Carangidae <- 
  groupings_model_fixed_sr_Carangidae %>%
  ggplot(aes(x=study_locations,
             y= response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  # geom_point(data = data_chao_s_Carangidae,
  #            aes(x = study_locations,
  #                y = !!response_var,
  #                shape = habitat
  #            ),
  #            position = position_jitterdodge(),
  #            color = "grey50",
  #            # shape = 1,
  #            size = 3) +
geom_point(data = data_chao_s_Carangidae,
           aes(x = study_locations,
               y = !!response_var,
               shape = habitat
           ),
           position = position_jitterdodge(
             jitter.width = 0.648,
             dodge.width = 0.7
             # jitter.height = 0.05
           ),
           color = "grey50",
           # shape = 1,
           size = 3) +
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "black",
                size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  # geom_text(aes(label=group),
  #           position = position_dodge(width=0.9),
  #           vjust = -0.5,
  #           hjust = -0.15,
  #           size = 8 / (14/5)) +  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  theme_classic() +
  labs(x = "Study Locations",
       y = "Observed Species Richness",
       title = "Carangidae") +
  ylim(0,5) +
  # labs(title = "Carangidae",
  #      subtitle = "Distribution Family = Poisson",
  #      x = "Study Locations",
  #      y = "EM Means of Species Richness") +
  theme(legend.position=c(0.5,0.8),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors,
                    labels = c("Shallow",
                               "Mesophotic"))

p_sr_Carangidae

#### mean_chao_s of Galeomorphii ####
## Make New Data Vegan for Galeomorphii
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
         study_locations = factor(study_locations,
                                  levels = c("TRNP",
                                             "CAGAYANCILLO")),
         habitat = factor(habitat),
         bait_type = factor(bait_type),
         site = factor(site),
         survey_area = factor(survey_area),
         habitat_mpa = str_c(habitat,
                             study_locations,
                             sep = " "))

attach(data_vegan.env)

pool_Galeomorphii <- 
  estimateR(x = data_vegan_Galeomorphii) %>%
  t() %>%
  as_tibble()


data_chao_s_Galeomorphii <- pool_Galeomorphii %>%
  clean_names() %>%
  bind_cols(data_vegan.env)

# data_chao_s_Galeomorphii %>%
#   ggplot(aes(x = s_chao1)) +
#   geom_histogram()

## Enter Information About Your Data for A Hypothesis Test ##


# define your response variable, here it is binomial
response_var = quo(s_obs) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"

alpha_sig = 0.05


## Histogram and Visualizing Distance Matrix
# data_chao_s_Galeomorphii %>%
#   ggplot(aes(x = s_obs)) + 
#   geom_histogram() +
#   facet_grid(habitat ~ study_locations)
# 
# data_chao_s_Galeomorphii %>%
#   ggplot(aes(x = s_chao1)) + 
#   geom_histogram() +
#   facet_grid(habitat ~ study_locations)
# 
# vis_dists(data_chao_s_Galeomorphii,
#           "s_chao1")
# 
# vis_dists(data_chao_s_Galeomorphii,
#           "s_obs")


sampling_design = "s_obs ~  habitat * study_locations"
#fit glm model
model_Galeomorphii <<- 
  glm(formula = sampling_design, 
      family = distribution_family,
      data = data_chao_s_Galeomorphii)


# sampling_design = "s_obs ~  habitat * study_locations + (1|study_locations:bait_type)"


# # # fit mixed model
# model_Galeomorphii <<-
#   afex::mixed(formula = sampling_design,
#               family = distribution_family,
#               method = "LRT",
#               sig_symbols = rep("", 4),
#               # all_fit = TRUE,
#               data = data_chao_s_Galeomorphii)

model_Galeomorphii
anova(model_Galeomorphii)

# visualize summary(model)
emmip(model_Galeomorphii, 
      study_locations ~ habitat,    # type = "response" for back transformed values
      cov.reduce = range) +
  # geom_vline(xintercept=mean(data_all_summaxn_$primer_x),
  #            color = "grey",
  #            linetype = "dashed") +
  # geom_text(aes(x = mean(data_all_summaxn_$primer_x),
  #               y = -2,
  #               label = "mean primer_x"),
  #           color = "grey") +
  theme_classic() +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")

## mean_chao_s Galeomorphii: Conduct A priori contrast tests for differences among sites ##
emmeans_model_sr_Galeomorphii <<-
  emmeans(model_Galeomorphii,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_sr_Galeomorphii,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_sr_Galeomorphii), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")

## mean_chao_s: Group Sites Based on Model Results ##
groupings_model_sr_Galeomorphii <<-
  multcomp::cld(emmeans_model_sr_Galeomorphii, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_sr_Galeomorphii  
# these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_sr_Galeomorphii <<-
  summary(emmeans_model_sr_Galeomorphii,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_sr_Galeomorphii %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_sr_Galeomorphii <- groupings_model_fixed_sr_Galeomorphii %>%
  mutate(habitat = factor(habitat,
                          levels = c(
                            "Shallow Reef",
                            "Mesophotic Reef")))

habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Mesophotic Reef")


## mean_chao_s Galeomorphii: Visualize Estimated Marginal Means Output with Group Categories ##
p_sr_Galeomorphii <- 
  groupings_model_fixed_sr_Galeomorphii %>%
  ggplot(aes(x=study_locations,
             y= response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  geom_point(data = data_chao_s_Galeomorphii,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat
             ),
             position = position_jitterdodge(),
             # color = "grey70",
             # shape = 1,
             size = 3) +
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "grey50",
                # size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  # geom_text(aes(label=group),
  #           position = position_dodge(width=0.9),
  #           vjust = -0.5,
  #           hjust = -0.15,
  #           size = 8 / (14/5)) +  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  theme_classic() +
  labs(x = "Study Locations",
       y = "Species Richness (Chao1)",
       title = "Galeomorphii") +
  ylim(0,5) +
  # labs(title = "Galeomorphii",
  #      subtitle = "Distribution Family = Poisson",
  #      x = "Study Locations",
  #      y = "EM Means of Species Richness") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors,
                    labels = c("Shallow",
                               "Mesophotic"))

p_sr_Galeomorphii

#### mean_chao_s of Cheilinus undulatus ####
## Make New Data Vegan for Cheilinus_undulatus
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
         study_locations = factor(study_locations,
                                  levels = c("TRNP",
                                             "CAGAYANCILLO")),
         habitat = factor(habitat),
         bait_type = factor(bait_type),
         site = factor(site),
         survey_area = factor(survey_area),
         habitat_mpa = str_c(habitat,
                             study_locations,
                             sep = " "))

attach(data_vegan.env)

pool_Cheilinus_undulatus <- 
  estimateR(x = data_vegan_Cheilinus_undulatus) %>%
  t() %>%
  as_tibble()


data_chao_s_Cheilinus_undulatus <- pool_Cheilinus_undulatus %>%
  clean_names() %>%
  bind_cols(data_vegan.env)

# data_chao_s_Cheilinus_undulatus %>%
#   ggplot(aes(x = s_chao1)) +
#   geom_histogram()

## Enter Information About Your Data for A Hypothesis Test ##


# define your response variable, here it is binomial
response_var = quo(s_obs) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"

alpha_sig = 0.05


## Histogram and Visualizing Distance Matrix
# data_chao_s_Cheilinus_undulatus %>%
#   ggplot(aes(x = s_obs)) + 
#   geom_histogram() +
#   facet_grid(habitat ~ study_locations)
# 
# data_chao_s_Cheilinus_undulatus %>%
#   ggplot(aes(x = s_chao1)) + 
#   geom_histogram() +
#   facet_grid(habitat ~ study_locations)
# 
# vis_dists(data_chao_s_Cheilinus_undulatus,
#           "s_chao1")
# 
# vis_dists(data_chao_s_Cheilinus_undulatus,
#           "s_obs")


sampling_design = "s_obs ~  habitat * study_locations"
#fit glm model
model_Cheilinus_undulatus <<- 
  glm(formula = sampling_design, 
      family = distribution_family,
      data = data_chao_s_Cheilinus_undulatus)


# sampling_design = "s_obs ~  habitat * study_locations + (1|study_locations:bait_type)"


# # # fit mixed model
# model_Cheilinus_undulatus <<-
#   afex::mixed(formula = sampling_design,
#               family = distribution_family,
#               method = "LRT",
#               sig_symbols = rep("", 4),
#               # all_fit = TRUE,
#               data = data_chao_s_Cheilinus_undulatus)

model_Cheilinus_undulatus
anova(model_Cheilinus_undulatus)

# visualize summary(model)
emmip(model_Cheilinus_undulatus, 
      study_locations ~ habitat,    # type = "response" for back transformed values
      cov.reduce = range) +
  # geom_vline(xintercept=mean(data_all_summaxn_$primer_x),
  #            color = "grey",
  #            linetype = "dashed") +
  # geom_text(aes(x = mean(data_all_summaxn_$primer_x),
  #               y = -2,
  #               label = "mean primer_x"),
  #           color = "grey") +
  theme_classic() +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")

## mean_chao_s Cheilinus_undulatus: Conduct A priori contrast tests for differences among sites ##
emmeans_model_sr_Cheilinus_undulatus <<-
  emmeans(model_Cheilinus_undulatus,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_sr_Cheilinus_undulatus,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_sr_Cheilinus_undulatus), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")

## mean_chao_s: Group Sites Based on Model Results ##
groupings_model_sr_Cheilinus_undulatus <<-
  multcomp::cld(emmeans_model_sr_Cheilinus_undulatus, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_sr_Cheilinus_undulatus  
# these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_sr_Cheilinus_undulatus <<-
  summary(emmeans_model_sr_Cheilinus_undulatus,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_sr_Cheilinus_undulatus %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_sr_Cheilinus_undulatus <- groupings_model_fixed_sr_Cheilinus_undulatus %>%
  mutate(habitat = factor(habitat,
                          levels = c(
                            "Shallow Reef",
                            "Mesophotic Reef")))

habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Mesophotic Reef")


## mean_chao_s Cheilinus_undulatus: Visualize Estimated Marginal Means Output with Group Categories ##
p_sr_Cheilinus_undulatus <- 
  groupings_model_fixed_sr_Cheilinus_undulatus %>%
  ggplot(aes(x=study_locations,
             y= response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  # geom_point(data = data_chao_s_Cheilinus_undulatus,
  #            aes(x = study_locations,
  #                y = !!response_var
  #            ),
  #            position = position_jitterdodge(),
  #            # color = "grey70",
  #            # shape = 1,
  #            size = 1) +
geom_point(data = data_chao_s_Cheilinus_undulatus,
           aes(x = study_locations,
               y = !!response_var,
               shape = habitat
           ),
           position = position_jitterdodge(
             jitter.width = 0.6,
             dodge.width = 0.7
             # jitter.height = 0.05
           ),
           color = "grey50",
           # shape = 1,
           size = 3) +
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "grey50",
                # size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  geom_text(aes(label=group),
            position = position_dodge(width=0.9),
            vjust = -0.5,
            hjust = -0.15,
            size = 8 / (14/5)) +  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  theme_classic() +
  # ylim(ymin, 
  #      ymax) +
  labs(title = "Cheilinus undulatus",
       x = "Study Locations",
       y = "EM Means of Species Richness") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors)

p_sr_Cheilinus_undulatus

#### Histograms ####
## For Overall Species Richness
histogram_overall <- 
  data_chao_s %>%
  ggplot() +
  aes(x = s_chao1,
      fill = habitat) +
  geom_histogram(alpha = 0.5,
                 position = "dodge",
                 binwidth = 10) +
  theme_classic() +
  facet_wrap(~study_locations) +
labs(subtitle = "Overall Species Richness",
     x = "Chao1 Species Richness",
     y = "Count",
     fill = "Depth Category") +
  scale_fill_manual(values = habitatcolors)
# facet_grid(habitat~study_locations)

histogram_overall

histogram_Serranidae <- 
  data_chao_s_Serranidae %>%
  bind_cols() %>%
  filter() %>%
  ggplot() +
  aes(x = s_obs,
      fill = habitat) +
  geom_histogram(alpha = 0.5,
                 position = "dodge",
                 binwidth = 1) +
  theme_classic() +
  facet_wrap(~study_locations) +
  guides(fill = "none") +
  labs(subtitle = "Serranidae",
       x = "Observed Species Richness",
       y = "Count") +
  scale_fill_manual(values = habitatcolors)
# facet_grid(habitat~study_locations)

histogram_Serranidae

histogram_Lutjanidae <- 
  data_chao_s_Lutjanidae %>%
  bind_cols() %>%
  filter() %>%
  ggplot() +
  aes(x = s_obs,
      fill = habitat) +
  geom_histogram(alpha = 0.5,
                 position = "dodge",
                 binwidth = .75) +
  theme_classic() +
  facet_wrap(~study_locations) +
  guides(fill = "none") +
  labs(subtitle = "Lutjanidae",
       x = "Observed Species Richness",
       y = "Count") +
  scale_fill_manual(values = habitatcolors)
# facet_grid(habitat~study_locations)

histogram_Lutjanidae

histogram_Lethrinidae <- 
  data_chao_s_Lethrinidae %>%
  bind_cols() %>%
  filter() %>%
  ggplot() +
  aes(x = s_obs,
      fill = habitat) +
  geom_histogram(alpha = 0.5,
                 position = "dodge",
                 binwidth = .75) +
  theme_classic() +
  facet_wrap(~study_locations) +
  guides(fill = "none") +
  labs(subtitle = "Lethrinidae",
       x = "Observed Species Richness",
       y = "Count") +
  scale_fill_manual(values = habitatcolors)
# facet_grid(habitat~study_locations)

histogram_Lethrinidae

histogram_Carangidae <- 
  data_chao_s_Carangidae %>%
  bind_cols() %>%
  filter() %>%
  ggplot() +
  aes(x = s_obs,
      fill = habitat) +
  geom_histogram(alpha = 0.5,
                 position = "dodge",
                 binwidth = 1) +
  theme_classic() +
  facet_wrap(~study_locations) +
  guides (fill = "none") +
  labs(subtitle = "Carangidae",
       x = "Observed Species Richness",
       y = "Count") +
  scale_fill_manual(values = habitatcolors)
# facet_grid(habitat~study_locations)

histogram_Carangidae

## Compile and save all histograms into one plot

histograms_sr <- ggarrange(histogram_overall,
                           histogram_Serranidae,
                           histogram_Lutjanidae,
                           histogram_Lethrinidae,
                           histogram_Carangidae,
                           ncol = 2,
                           nrow = 3)

ggsave("FacetedHistogram_SR.png",
       histograms_sr, 
       height = 14,
       width = 12,
       units = "in")


#### Save Overall Species Richness Plot ####
emmeans_sr <- ggarrange(p_sr,
                        p_sr_Serranidae,
                        p_sr_Lutjanidae,
                        p_sr_Lethrinidae, 
                        p_sr_Carangidae,
                        ncol = 2,
                        nrow = 3)
ggsave("FacetedEmMeansSpeciesRichness.pdf", 
       emmeans_sr, height = 14, width = 12, units = "in")
ggsave("FacetedEmMeansSpeciesRichness.png", 
       emmeans_sr, height = 14, width = 12, units = "in")

emmeans_sr_reorganized <- ggarrange(p_sr,
                                    p_sr_Serranidae,
                                    p_sr_Lutjanidae,
                                    p_sr_Lethrinidae, 
                                    p_sr_Carangidae,
                                    ncol = 2,
                                    nrow = 3)
ggsave("FacetedSpeciesRichnesswOverall.png",
       emmeans_sr_reorganized,
       height = 11,
       width = 8.5,
       units = "in")

#### Testing for Bait Type on Species Richness ####
pool <- 
  estimateR(x = data_vegan) %>%
  t() %>%
  as_tibble()

data_chao_s <- 
  pool %>%
  clean_names() %>%
  bind_cols(data_vegan.env) 

##Bait Type Effect on TRNP 
data_chao_s_TRNP <- 
  pool %>%
  clean_names() %>%
  bind_cols(data_vegan.env) %>%
  filter(study_locations == "TRNP")

##Visualize using boxplot##
ggplot(data_chao_s_TRNP, aes(bait_type, s_chao1, fill=habitat)) +
  geom_boxplot() 

## Enter Information About Your Data for A Hypothesis Test ##

#set sampling design
sampling_design = "s_chao1 ~ bait_type * habitat"

distribution_family = "Gamma"

#glm model (afex::mixed requires a random factor, which I'm not sure is appropriate here)
model_bait_TRNP <<-
  glm(formula=sampling_design,
      family = distribution_family, 
      data = data_chao_s_TRNP)

summary(model_bait_TRNP)


##Bait Type Effect on Cag
data_chao_s_Cag <- 
  pool %>%
  clean_names() %>%
  bind_cols(data_vegan.env) %>%
  filter(study_locations == "CAGAYANCILLO") %>%
  filter(bait_type != "Black Jack/Bluefin Trevally")

##Visualize using boxplot##
ggplot(data_chao_s_Cag, aes(bait_type, s_chao1, fill=habitat)) +
  geom_boxplot() 

## Enter Information About Your Data for A Hypothesis Test ##

#set sampling design
sampling_design = "s_chao1 ~ bait_type * habitat"

distribution_family = "Gamma"

#glm model (afex::mixed requires a random factor, which I'm not sure is appropriate here)
model_bait_Cag <<-
  glm(formula=sampling_design,
      family = distribution_family, 
      data = data_chao_s_Cag)

summary(model_bait_Cag)

