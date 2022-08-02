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
           .keep_all = TRUE)

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
           family == "Sphyrnidae" ~ "Galeomorphii",
           family == "Carcharhinidae" ~ "Galeomorphii",
           family == "Alopiidae" ~ "Galeomorphii",
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
           .keep_all = TRUE) 


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
    family == "Sphyrnidae" ~ "Galeomorphii",
    family == "Carcharhinidae" ~ "Galeomorphii",
    family == "Alopiidae" ~ "Galeomorphii",
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

vis_dists(pool,
          "S.chao1")
par(mar = c(1,1,1,1))

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
  scale_fill_manual(values = habitatcolors,
                    labels = c("Shallow Reef",
                               "Mesophotic Reef")) +
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
  bind_cols(data_vegan.env)



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
  labs(x = "Study Locations",
       y = "Estimated Marginal Means of Chao Estimate of Species Richness") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors)

p_sr
save_plot("EMMeansofSpeciesRichness.png")

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

view(data_vegan_Serranidae)

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

vis_dists(pool_Serranidae,
          "S.chao1")

data_chao_s_Serranidae <- pool_Serranidae %>%
  clean_names() %>%
  bind_cols(data_vegan.env)

## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(s_chao1) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "Gamma"


alpha_sig = 0.05


# we start with the loci subjected to 11 primer concentrations (we removed loci with no sum_max_n to simplify)


sampling_design = "s_chao1 ~  habitat * study_locations + (1|study_locations:bait_type)"


# # fit mixed model
model_Serranidae <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_chao_s_Serranidae)

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

## mean_chao_s Serranidae: Conduct A priori contrast tests for differences among sites ####
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
