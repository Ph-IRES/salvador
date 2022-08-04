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


metadata <-
  read_excel(inFilePath4,
             na="NA") %>%
  clean_names() %>%
  dplyr::rename(bait_weight_grams = weight_grams)

habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Deep Reef")

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
    family == "Sphyrnidae" ~ "Galeoidea",
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
#### sum_max_n: Make Visualization of Hypothesis Test ####
response_var = quo(sum_max_n)

habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Mesophotic Reef")

data_all_summaxn <- 
  data_all %>%
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
  group_by(op_code,
           study_locations,
           habitat,
           bait_type) %>%
  dplyr::summarize(sum_max_n = sum(max_n))
View(data_all_summaxn)

data_all %>%
  group_by(study_locations, 
           habitat,
           op_code) %>% 
  dplyr::summarize(sum_max_n = sum(max_n)) %>%
  dplyr::summarize(mean_sum_max_n = mean(sum_max_n),
                   se_sum_max_n = sd(sum_max_n)/sqrt(n())) %>% 
  ggplot(aes(x = study_locations,
             y = mean_sum_max_n,
             fill = habitat))+
  geom_bar(position = "dodge", 
           stat = "identity") +
  geom_point(data = data_all_summaxn,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat),
             position=position_jitterdodge(),
             size = 3,
             color = "darkgrey",
             inherit.aes = FALSE) +
  xlab("Study Locations") +
  ylab("Mean MaxN per BRUV Deployment") +
  labs(title = "Mean Sum of MaxN at TRNP vs. Cagayancillo",
       fill = "Habitat") +
  theme_classic() +
  scale_fill_manual(values = habitatcolors, 
                    labels = c("Shallow Reef",
                               "Mesophotic Reef")) +
  geom_errorbar(aes(ymax = mean_sum_max_n + se_sum_max_n,
                    ymin = mean_sum_max_n - se_sum_max_n), 
                position = "dodge") +
  scale_y_continuous(trans='log10') 


#### sum_max_n: Mixed Effects Hypothesis Test ####

data_all_summaxn <- 
  data_all %>%
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
  group_by(op_code,
           study_locations,
           habitat,
           bait_type) %>%
  dplyr::summarize(sum_max_n = sum(max_n))

## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(sum_max_n) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"


alpha_sig = 0.05

# we start with the loci subjected to 11 primer concentrations (we removed loci with no sum_max_n to simplify)


sampling_design = "sum_max_n ~  habitat * study_locations + (1|study_locations:bait_type)"

# # fit mixed model
model <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all_summaxn)

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

# sum_max_n: Conduct A priori contrast tests for differences among sites ####

# now we move on to finish the hypothesis testing.  Are there differences between the sites?
# estimated marginal means 

emmeans_model <<-
  emmeans(model,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")


#### sum_max_n: Group Sites Based on Model Results ####

groupings_model <<-
  multcomp::cld(emmeans_model, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model             # these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed <<-
  summary(emmeans_model,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed  <- groupings_model_fixed %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef")))     # cld messes up back transformation, this takes values from emmeans and groupings from cld


#### sum_max_n: Visualize Estimated Marginal Means Output With Group Categories ####

p <- 
  groupings_model_fixed %>%
  ggplot(aes(x=study_locations,
             y=response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "grey50",
                # size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  geom_point(data = data_all_summaxn,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat),
             position=position_jitterdodge(),
             size = 3,
             color = "darkgrey",
             inherit.aes = FALSE) +
  geom_text(aes(label=group),
            position = position_dodge(width=0.9),
            vjust = -0.5,
            hjust = -0.15,
            size = 8 / (14/5)) +
  scale_y_continuous(trans='log10') +
  theme_classic() +
  # ylim(ymin, 
  #      ymax) +
  labs(title = sampling_design,
       subtitle = "Distribution Family = Poisson",
       x = "Depth",
       y = "Mean Sum_Max_N") +
  theme(legend.position=c(0.33,0.9),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors)

p
save_plot("EstimatedMarginalMeansMaxN.png")
#### sum_max_n: Serranidae ####
data_all_summaxn_Serranidae <- 
  data_all %>%
  filter(groupings == "Serranidae") %>%
  mutate(study_locations = factor(study_locations,
                                  levels = c("TRNP", 
                                             "CAGAYANCILLO"))) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Deep Reef"))) %>%
  mutate(habitat = case_when(
    habitat == "Shallow Reef" ~ "Shallow Reef",
    habitat == "Deep Reef" ~ "Mesophotic Reef")) %>%
  
  group_by(op_code,
           habitat,
           study_locations) %>%
  summarize(sum_max_n = sum(max_n)) %>%

  mutate(opcode_habitat = str_c(op_code,
                                habitat,
                                sep = ",")) %>% 
  ungroup() %>%
  dplyr::select(-op_code,
                -habitat) %>% 
  complete(study_locations,
           opcode_habitat,
           fill = list(sum_max_n = 0)) %>%
  separate(opcode_habitat,
           into = c("op_code",
                    "habitat"),
           sep = ",") %>% 
  left_join(data_all %>%
              dplyr::select(op_code,
                            bait_type) %>%
              distinct()) 

View(data_all_summaxn_Serranidae)

## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(sum_max_n) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"


alpha_sig = 0.05

# we start with the loci subjected to 11 primer concentrations (we removed loci with no sum_max_n to simplify)


sampling_design = "sum_max_n ~  habitat * study_locations + (1|study_locations:bait_type)"

# # fit mixed model
model_Serranidae <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all_summaxn_Serranidae)

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
       subtitle = "Serranidae",
       y = "Linear Prediciton",
       x = "MPA")

## sum_max_n Serranidae: Conduct A priori contrast tests for differences among sites ##

# now we move on to finish the hypothesis testing.  Are there differences between the sites?
# estimated marginal means 

emmeans_model_Serranidae <<-
  emmeans(model_Serranidae,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_Serranidae,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_Serranidae), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")

## sum_max_n Serranidae: Group Sites Based on Model Results ##

groupings_model_Serranidae <<-
  multcomp::cld(emmeans_model_Serranidae, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_Serranidae             # these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_Serranidae <<-
  summary(emmeans_model_Serranidae,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_Serranidae %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_Serranidae  <- groupings_model_fixed_Serranidae %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef")))# cld messes up back transformation, this takes values from emmeans and groupings from cld

## sum_max_n Serranidae: Visualize Estimated Marginal Means Output With Group Categories 

p_Serranidae <- 
  groupings_model_fixed_Serranidae %>%
  ggplot(aes(x=study_locations,
             y=response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "grey50",
                # size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  geom_point(data = data_all_summaxn_Serranidae,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat),
             position=position_jitterdodge(),
             size = 3,
             color = "darkgrey",
             inherit.aes = FALSE) +
  geom_text(aes(label=group),
            position = position_dodge(width=0.9),
            vjust = -0.5,
            hjust = -0.15,
            size = 8 / (14/5)) +
  scale_y_continuous(trans='log10') +
  theme_classic() +
  # ylim(ymin, 
  #      ymax) +
  labs(title = "Serranidae",
       subtitle = "Distribution Family = Poisson",
       x = "Depth",
       y = "Estimated Marginal Means of Sum of MaxN") +
  theme(legend.position=c(0.33,0.85),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors)

p_Serranidae



#### sum_max_n: Lutjanidae ####
data_all_summaxn_Lutjanidae <- 
  data_all %>%
  filter(groupings == "Lutjanidae") %>%
  mutate(study_locations = factor(study_locations,
                                  levels = c("TRNP", 
                                             "CAGAYANCILLO"))) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Deep Reef"))) %>%
  mutate(habitat = case_when(
    habitat == "Shallow Reef" ~ "Shallow Reef",
    habitat == "Deep Reef" ~ "Mesophotic Reef")) %>%
  
  group_by(op_code,
           habitat,
           study_locations) %>%
  summarize(sum_max_n = sum(max_n)) %>%
  
  mutate(opcode_habitat = str_c(op_code,
                                habitat,
                                sep = ",")) %>% 
  ungroup() %>%
  dplyr::select(-op_code,
                -habitat) %>% 
  complete(study_locations,
           opcode_habitat,
           fill = list(sum_max_n = 0)) %>%
  separate(opcode_habitat,
           into = c("op_code",
                    "habitat"),
           sep = ",") %>% 
  left_join(data_all %>%
              dplyr::select(op_code,
                            bait_type) %>%
              distinct()) 


## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(sum_max_n) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"


alpha_sig = 0.05


sampling_design = "sum_max_n ~  habitat * study_locations + (1|study_locations:bait_type)"


# # fit mixed model
model_Lutjanidae <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all_summaxn_Lutjanidae)

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
       subtitle = "Lutjanidae",
       y = "Linear Prediciton",
       x = "MPA")

## sum_max_n Lutjanidae: Conduct A priori contrast tests for differences among sites ##

# now we move on to finish the hypothesis testing.  Are there differences between the sites?
# estimated marginal means 

emmeans_model_Lutjanidae <<-
  emmeans(model_Lutjanidae,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_Lutjanidae,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_Lutjanidae), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")

## sum_max_n Lutjanidae: Group Sites Based on Model Results ##

groupings_model_Lutjanidae <<-
  multcomp::cld(emmeans_model_Lutjanidae, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_Lutjanidae             # these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_Lutjanidae <<-
  summary(emmeans_model_Lutjanidae,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_Lutjanidae %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_Lutjanidae    <- groupings_model_fixed_Lutjanidae %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef")))   # cld messes up back transformation, this takes values from emmeans and groupings from cld

## sum_max_n Lutjanidae: Visualize Estimated Marginal Means Output With Group Categories 
habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Mesophotic Reef")
p_Lutjanidae <- 
  groupings_model_fixed_Lutjanidae %>%
  ggplot(aes(x=study_locations,
             y=response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "grey50",
                # size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
   # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  geom_point(data = data_all_summaxn_Lutjanidae,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat),
             position=position_jitterdodge(),
             size = 3,
             color = "darkgrey",
             inherit.aes = FALSE) +
  geom_text(aes(label=group),
            position = position_dodge(width=0.9),
            vjust = -0.5,
            hjust = -0.15,
            size = 8 / (14/5)) + 
  scale_y_continuous(trans='log10') +
  theme_classic() +
  # ylim(ymin, 
  #      ymax) +
  labs(title = "Lutjanidae",
       subtitle = "Distribution Family = Poisson",
       x = "Depth",
       y = "Estimated Marginal Means of Sum of MaxN") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors)

p_Lutjanidae

#### sum_max_n: Lethrinidae ####
data_all_summaxn_Lethrinidae <- 
  data_all %>%
  filter(groupings == "Lethrinidae") %>%
  mutate(study_locations = factor(study_locations,
                                  levels = c("TRNP", 
                                             "CAGAYANCILLO"))) %>%
  mutate(habitat = case_when(
    habitat == "Shallow Reef" ~ "Shallow Reef",
    habitat == "Deep Reef" ~ "Mesophotic Reef")) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef"))) %>%
  group_by(op_code,
           habitat,
           study_locations) %>%
  summarize(sum_max_n = sum(max_n)) %>%
  
  mutate(opcode_habitat = str_c(op_code,
                                habitat,
                                sep = ",")) %>% 
  ungroup() %>%
  dplyr::select(-op_code,
                -habitat) %>% 
  complete(study_locations,
           opcode_habitat,
           fill = list(sum_max_n = 0)) %>%
  separate(opcode_habitat,
           into = c("op_code",
                    "habitat"),
           sep = ",") %>% 
  left_join(data_all %>%
              dplyr::select(op_code,
                            bait_type) %>%
              distinct()) 


## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(sum_max_n) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"


alpha_sig = 0.05


sampling_design = "sum_max_n ~  habitat * study_locations + (1|study_locations:bait_type)"

# # fit mixed model
model_Lethrinidae <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all_summaxn_Lethrinidae)

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
       subtitle = "Lethrinidae",
       y = "Linear Prediciton",
       x = "MPA")

## sum_max_n Lethrinidae: Conduct A priori contrast tests for differences among sites ##

# now we move on to finish the hypothesis testing.  Are there differences between the sites?
# estimated marginal means 

emmeans_model_Lethrinidae <<-
  emmeans(model_Lethrinidae,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_Lethrinidae,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_Lethrinidae), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")

## sum_max_n Lethrinidae: Group Sites Based on Model Results ##

groupings_model_Lethrinidae <<-
  multcomp::cld(emmeans_model_Lethrinidae, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_Lethrinidae             # these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_Lethrinidae <<-
  summary(emmeans_model_Lethrinidae,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_Lethrinidae %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_Lethrinidae <-  groupings_model_fixed_Lethrinidae %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef")))    # cld messes up back transformation, this takes values from emmeans and groupings from cld

## sum_max_n Lethrinidae: Visualize Estimated Marginal Means Output With Group Categories 

p_Lethrinidae <- 
  groupings_model_fixed_Lethrinidae %>%
  ggplot(aes(x=study_locations,
             y=response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "grey50",
                # size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
    # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  geom_point(data = data_all_summaxn_Lethrinidae,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat),
             position=position_jitterdodge(),
             size = 3,
             color = "darkgrey",
             inherit.aes = FALSE) +
  geom_text(aes(label=group),
            position = position_dodge(width=0.9),
            vjust = -0.5,
            hjust = -0.15,
            size = 8 / (14/5)) +
  scale_y_continuous(trans='log10') +
  theme_classic() +
  # ylim(ymin, 
  #      ymax) +
  labs(title = "Lethrinidae",
       subtitle = "Distribution Family = Poisson",
       x = "Depth",
       y = "Estimated Marginal Means of Sum of MaxN") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors)

p_Lethrinidae

#### sum_max_n: Carangidae ####
data_all_summaxn_Carangidae <- 
  data_all_summaxn_Lethrinidae <- 
  data_all %>%
  filter(groupings == "Carangidae") %>%
  mutate(study_locations = factor(study_locations,
                                  levels = c("TRNP", 
                                             "CAGAYANCILLO"))) %>%
  mutate(habitat = case_when(
    habitat == "Shallow Reef" ~ "Shallow Reef",
    habitat == "Deep Reef" ~ "Mesophotic Reef")) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef"))) %>%
  group_by(op_code,
           habitat,
           study_locations) %>%
  summarize(sum_max_n = sum(max_n)) %>%
  
  mutate(opcode_habitat = str_c(op_code,
                                habitat,
                                sep = ",")) %>% 
  ungroup() %>%
  dplyr::select(-op_code,
                -habitat) %>% 
  complete(study_locations,
           opcode_habitat,
           fill = list(sum_max_n = 0)) %>%
  separate(opcode_habitat,
           into = c("op_code",
                    "habitat"),
           sep = ",") %>% 
  left_join(data_all %>%
              dplyr::select(op_code,
                            bait_type) %>%
              distinct()) 


## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(sum_max_n) # quo() allows column names to be put into variables 


# enter the distribution family for your response variable
distribution_family = "poisson"


alpha_sig = 0.05


sampling_design = "sum_max_n ~  habitat * study_locations + (1|study_locations:bait_type)"

# # fit mixed model
model_Carangidae <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all_summaxn_Carangidae)

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
       subtitle = "Carangidae",
       y = "Linear Prediciton",
       x = "MPA")


## sum_max_n Carangidae: Conduct A priori contrast tests for differences among sites ##

# now we move on to finish the hypothesis testing.  Are there differences between the sites?
# estimated marginal means 

emmeans_model_Carangidae <<-
  emmeans(model_Carangidae,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_Carangidae,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_Carangidae), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")

## sum_max_n Carangidae: Group Sites Based on Model Results ##

groupings_model_Carangidae <<-
  multcomp::cld(emmeans_model_Carangidae, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_Carangidae             # these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_Carangidae <<-
  summary(emmeans_model_Carangidae,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_Carangidae %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_Carangidae <- groupings_model_fixed_Carangidae %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef")))  # cld messes up back transformation, this takes values from emmeans and groupings from cld

## sum_max_n Carangidae: Visualize Estimated Marginal Means Output With Group Categories 

p_Carangidae <- 
  groupings_model_fixed_Carangidae %>%
  ggplot(aes(x=study_locations,
             y=response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "grey50",
                # size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  geom_point(data = data_all_summaxn_Carangidae,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat),
             position=position_jitterdodge(),
             size = 3,
             color = "darkgrey",
             inherit.aes = FALSE) +
  geom_text(aes(label=group),
            position = position_dodge(width=0.9),
            vjust = -0.5,
            hjust = -0.15,
            size = 8 / (14/5)) +  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  scale_y_continuous(trans='log10') +
  theme_classic() +
  # ylim(ymin, 
  #      ymax) +
  labs(title = "Carangidae",
       subtitle = "Distribution Family = Poisson",
       x = "Depth",
       y = "Estimated Marginal Means of Sum of MaxN") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) + 
  scale_fill_manual(values = habitatcolors)

p_Carangidae

#### sum_max_n: Galeomorphii ####
data_all_summaxn_Galeomorphii <- 
  data_all %>%
  mutate(groupings = case_when(
    groupings == "Galeoidea" ~ "Galeomorphii",
    TRUE ~ groupings
  )) %>%
  filter(groupings == "Galeomorphii") %>%
  mutate(study_locations = factor(study_locations,
                                  levels = c("TRNP", 
                                             "CAGAYANCILLO"))) %>%
  mutate(habitat = case_when(
    habitat == "Shallow Reef" ~ "Shallow Reef",
    habitat == "Deep Reef" ~ "Mesophotic Reef")) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef"))) %>%
  group_by(op_code,
           habitat,
           study_locations) %>%
  summarize(sum_max_n = sum(max_n)) %>%
  
  mutate(opcode_habitat = str_c(op_code,
                                habitat,
                                sep = ",")) %>% 
  ungroup() %>%
  dplyr::select(-op_code,
                -habitat) %>% 
  complete(study_locations,
           opcode_habitat,
           fill = list(sum_max_n = 0)) %>%
  separate(opcode_habitat,
           into = c("op_code",
                    "habitat"),
           sep = ",") %>% 
  left_join(data_all %>%
              dplyr::select(op_code,
                            bait_type) %>%
              distinct()) 


## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(sum_max_n) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"


alpha_sig = 0.05


sampling_design = "sum_max_n ~  habitat * study_locations + (1|study_locations:bait_type)"

# # fit mixed model
model_Galeomorphii <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all_summaxn_Galeomorphii)

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
       subtitle = "Galeomorphii",
       y = "Linear Prediciton",
       x = "MPA")

## sum_max_n Galeomorphii: Conduct A priori contrast tests for differences among sites ##

# now we move on to finish the hypothesis testing.  Are there differences between the sites?
# estimated marginal means 

emmeans_model_Galeomorphii <<-
  emmeans(model_Galeomorphii,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_Galeomorphii,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_Galeomorphii), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")

## sum_max_n Galeomorphii: Group Sites Based on Model Results ##

groupings_model_Galeomorphii <<-
  multcomp::cld(emmeans_model_Galeomorphii, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_Galeomorphii             # these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_Galeomorphii <<-
  summary(emmeans_model_Galeomorphii,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_Galeomorphii %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_Galeomorphii   <- groupings_model_fixed_Galeomorphii %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef")))     # cld messes up back transformation, this takes values from emmeans and groupings from cld

## sum_max_n Galeomorphii: Visualize Estimated Marginal Means Output With Group Categories 

p_Galeomorphii <- 
  groupings_model_fixed_Galeomorphii %>%
  ggplot(aes(x=study_locations,
             y=response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "grey50",
                # size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  geom_point(data = data_all_summaxn_Galeomorphii,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat),
             position=position_jitterdodge(),
             size = 3,
             color = "darkgrey",
             inherit.aes = FALSE) +
  geom_text(aes(label=group),
            position = position_dodge(width=0.9),
            vjust = -0.5,
            hjust = -0.15,
            size = 8 / (14/5)) +  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  scale_y_continuous(trans='log10') +
  theme_classic() +
  # ylim(ymin, 
  #      ymax) +
  labs(title = "Galeomorphii",
       subtitle = "Distribution Family = Poisson",
       x = "Depth",
       y = "Estimated Marginal Means of Sum of MaxN") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors)

p_Galeomorphii

#### sum_max_n: Cheilinus undulatus ####
data_all_summaxn_Cheilinus_undulatus <- 
  data_all %>%
  filter(groupings == "Cheilinus undulatus") %>%
  mutate(study_locations = factor(study_locations,
                                  levels = c("TRNP", 
                                             "CAGAYANCILLO"))) %>%
  mutate(habitat = case_when(
    habitat == "Shallow Reef" ~ "Shallow Reef",
    habitat == "Deep Reef" ~ "Mesophotic Reef")) %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef"))) %>%
  group_by(op_code,
           habitat,
           study_locations) %>%
  summarize(sum_max_n = sum(max_n)) %>%
  
  mutate(opcode_habitat = str_c(op_code,
                                habitat,
                                sep = ",")) %>% 
  ungroup() %>%
  dplyr::select(-op_code,
                -habitat) %>% 
  complete(study_locations,
           opcode_habitat,
           fill = list(sum_max_n = 0)) %>%
  separate(opcode_habitat,
           into = c("op_code",
                    "habitat"),
           sep = ",") %>% 
  left_join(data_all %>%
              dplyr::select(op_code,
                            bait_type) %>%
              distinct()) 

## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(sum_max_n) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"


alpha_sig = 0.05


sampling_design = "sum_max_n ~  habitat * study_locations + (1|study_locations:bait_type)"

# # fit mixed model
model_Cheilinus_undulatus <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all_summaxn_Cheilinus_undulatus)

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
       subtitle = "Cheilinus_undulatus",
       y = "Linear Prediciton",
       x = "MPA")

## sum_max_n Cheilinus_undulatus: Conduct A priori contrast tests for differences among sites ##

# now we move on to finish the hypothesis testing.  Are there differences between the sites?
# estimated marginal means 

emmeans_model_Cheilinus_undulatus <<-
  emmeans(model_Cheilinus_undulatus,
          ~ habitat * study_locations,
          alpha = alpha_sig)

# emmeans back transformed to the original units of response var
summary(emmeans_model_Cheilinus_undulatus,      
        type="response")

# contrasts between sites
contrast(regrid(emmeans_model_Cheilinus_undulatus), # emmeans back transformed to the original units of response var
         method = 'pairwise', 
         simple = 'each', 
         combine = FALSE, 
         adjust = "bh")

## sum_max_n Cheilinus_undulatus: Group Sites Based on Model Results ##

groupings_model_Cheilinus_undulatus <<-
  multcomp::cld(emmeans_model_Cheilinus_undulatus, 
                alpha = alpha_sig,
                Letters = letters,
                type="response",
                adjust = "bh") %>%
  as.data.frame %>%
  mutate(group = str_remove_all(.group," "),
         group = str_replace_all(group,
                                 "(.)(.)",
                                 "\\1,\\2")) 

groupings_model_Cheilinus_undulatus             # these values are back transformed, groupings based on transformed


# i noticed that the emmeans from groupings don't match those from emmeans so this is the table to use for making the figure
# the emmeans means and conf intervals match those produced by afex_plot, so I think those are what we want
groupings_model_fixed_Cheilinus_undulatus <<-
  summary(emmeans_model_Cheilinus_undulatus,      # emmeans back transformed to the original units of response var
          type="response") %>%
  tibble() %>%
  left_join(groupings_model_Cheilinus_undulatus %>%
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_Cheilinus_undulatus <- groupings_model_fixed_Cheilinus_undulatus %>%
  mutate(habitat = factor(habitat,
                          levels = c("Shallow Reef",
                                     "Mesophotic Reef")))      # cld messes up back transformation, this takes values from emmeans and groupings from cld

## sum_max_n Cheilinus_undulatus: Visualize Estimated Marginal Means Output With Group Categories 

p_Cheilinus_undulatus <- 
  groupings_model_fixed_Cheilinus_undulatus %>%
  ggplot(aes(x=study_locations,
             y=response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  
  geom_errorbar(aes(ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width = 0.2,
                color = "grey50",
                # size = 1,
                position = position_dodge(width=0.9)) +
  guides(color = "none",
         shape = "none") +   #remove color legend
  geom_point(data = data_all_summaxn_Cheilinus_undulatus,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat),
             position=position_jitterdodge(),
             size = 3,
             color = "darkgrey",
             inherit.aes = FALSE) +
  geom_text(aes(label=group),
            position = position_dodge(width=0.9),
            vjust = -0.5,
            hjust = -0.15,
            size = 8 / (14/5)) +  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  scale_y_continuous(trans='log10') +
  theme_classic() +
  # ylim(ymin, 
  #      ymax) +
  labs(title = "Cheilinus_undulatus",
       subtitle = "Distribution Family = Poisson",
       x = "Depth",
       y = "Estimated Marginal Means of Sum of MaxN") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) +
  scale_fill_manual(values = habitatcolors)

p_Cheilinus_undulatus

emmeans_maxn <- ggarrange(p_Serranidae, p_Lutjanidae, p_Lethrinidae, p_Carangidae, p_Galeomorphii, p_Cheilinus_undulatus, 
          ncol = 2,
          nrow = 3)
ggsave("FacetedEmMeansMaxN.pdf", 
       emmeans_maxn, height = 11, width = 8.5, units = "in")
