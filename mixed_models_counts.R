#### INITIALIZATION ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(janitor)
library(magrittr)
library(gridExtra)
library(readxl)

# install.packages("rlang")
# install.packages("emmeans")
library(devtools)
library(vegan)
# devtools::install_github("singmann/afex@master")
# devtools::install_github("eclarke/ggbeeswarm")
# install.packages("multcomp")
# install.packages("multcompView")
# install.packages("performance")
# install.packages("fitdistrplus")
# install.packages("optimx")
# install.packages("effects")
# install.packages("prediction")
# install.packages("ggforce")

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

# NOTE: after loading these packages, you may find that tidyverse commands are affected 
#       the solution is to add the appropriate package name before commands that break code
#       such as `dplyr::select` if `select` doesn't work correctly anymore
#       this happens when multiple packages have the same command names. 
#       The last package loaded takes precidence, and your tidyverse commands break.
#       you could load tidyverse last, but invariably, you will load a package after tidyverse
#       so it's impossible to avoid this

#### USER DEFINED VARIABLES ####

# path to data set
inFilePath2 = "./meso_euphotic_carniv_fish_videobaitstations_all.rds"
functionPath = "./model_fitting_functions.R"
inFilePath3 = "./WorkingData_CLEANED_TUB,CAG.xlsx"
inFilePath4 = "./PHIRES_MetaData.xlsx"

# you can make a default theme for your publication's figures.  This makes things easier for you. 
# feel free to customize as necessary
theme_myfigs <- 
  theme_classic() +
  theme(panel.background = element_rect(fill = 'white', 
                                        color = 'white'),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color="grey95", 
                                          size=0.25),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 9, 
                                   color = 'black'),
        axis.text.x = element_text(size = 9, 
                                   color = 'black'),
        # axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, 
                                    color = 'black'),
        plot.title = element_text(size = 10, 
                                  color = 'black'),
        plot.subtitle = element_text(size = 9, 
                                     color = 'black'),
        plot.caption = element_text(size = 9, 
                                    color = 'black', 
                                    hjust = 0),
        legend.text = element_text(size = 9, 
                                   color = 'black'),
        legend.title = element_text(size = 9, 
                                    color = 'black'),
        legend.background = element_blank(),
        legend.position="right"
  )



#### READ IN DATA ####

# read in data and remove rows with missing data or multiple bands
data_all <-
  read_rds(inFilePath2) 
View(data_all)

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
  mutate(groupings = case_when(
    family == "Labridae" ~ "Cheilinus undulatus",
    family == "Sphyrnidae" ~ "Galeoidea",
    family == "Carcharhinidae" ~ "Galeoidea",
    family == "Alopiidae" ~ "Galeoidea",
    family == "Epinephelidae" ~ "Serranidae",
    TRUE ~ family))

#### FUNCTIONS ####
# visualize statistical distributions (see fitdistrplus: An R Package for Fitting Distributions, 2020)
source(functionPath)

#### Explore Your Data For the Hypothesis Tests ####

# it is important to understand the nature of your data
# histograms can help you make decisions on how your data must be treated to conform with the assumptions of statistical models

#histogram raw data
data_all %>%
  ggplot(aes(x = max_n)) +
  geom_histogram()+
  facet_grid(study_locations ~ habitat) +
  theme_myfigs +
  labs(x = "Max N per Species Per Video")

#histogram for the sum of max_n by op_code
data_all %>%
  group_by(op_code,
           study_locations,
           habitat) %>%
  dplyr::summarize(sum_max_n = sum(max_n)) %>%
  ggplot(aes(x = sum_max_n)) +
  geom_histogram()+
  theme_myfigs +
  facet_grid(study_locations ~ habitat) +
  labs(x= "Sum of Max N per Video")

#histogram for the sum of max_n by op_code
data_all %>%
  group_by(op_code,
           study_locations,
           habitat) %>%
  dplyr::summarize(sum_max_n = sum(max_n)) %>%
  ggplot(aes(x = op_code,
             y = sum_max_n)) +
  geom_col()+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(. ~ study_locations + habitat,
             scales = "free")


# visualize statistical distributions (see fitdistrplus: An R Package for Fitting Distributions, 2020)
#  vis_dists() is a function that I made above in the FUNCTIONS section.  It accepts the tibble and column name to visualize.
#  vis_dists() creates three figures

# results in error making third plot because some values are zero and some of the distibutions are incompatible with zeros in data
par(mar = c(1,1,1,1)) #to fix plot margins
vis_dists(data_all,
          "max_n")

vis_dists(data_all %>%
            mutate(log10_max_n = log10(max_n)),
          "log10_max_n")

vis_dists(data_all %>%
            group_by(op_code,
                     study_locations,
                     habitat) %>%
            dplyr::summarize(sum_max_n = sum(max_n)),
          "sum_max_n")


#### max_n: Make Visualization of Hypothesis Test ####
habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Deep Reef")

response_var = quo(max_n)

data_all %>%
  group_by(study_locations, 
           habitat,
           op_code) %>%
  dplyr::summarize(mean_max_n = mean(max_n)) %>%
  dplyr::summarize(mean_mean_max_n = mean(mean_max_n),
                   se_mean_max_n = sd(mean_max_n)/sqrt(n())) %>%
  ggplot(aes(x = study_locations,
             y = mean_mean_max_n,
             fill = habitat))+
  geom_bar(position = "dodge", 
           stat = "identity") +
  xlab("Study Locations") +
  ylab("Mean Sum of MaxN per BRUV Deployment") +
  labs(title = "Mean MaxN at TRNP vs. Cagayancillo",
       fill = "Habitat") +
  theme_classic() +
  scale_fill_manual(values = habitatcolors, 
                    labels = c("Shallow Reef",
                               "Mesophotic Reef")) +
  geom_errorbar(aes(ymax = mean_mean_max_n + se_mean_max_n,
                    ymin = mean_mean_max_n - se_mean_max_n), 
                position = "dodge") +
  geom_point(data = data_all,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat),
             position=position_jitterdodge(),
             size = 3,
             color = "darkgrey",
             inherit.aes = FALSE) +
  scale_y_continuous(trans='log10')


#### max_n: Mixed Effects Hypothesis Test ####

## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(max_n) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"


alpha_sig = 0.05

# we start with the full model, then reduce complexity
sampling_design = "max_n ~  habitat * study_locations + (1 | study_locations:bait_type) + (1|study_locations:bait_type:op_code)"
sampling_design2 = "max_n ~  habitat * study_locations + (1 | study_locations:bait_type)"
sampling_design3 = "max_n ~  habitat * study_locations + (1 | study_locations:op_code)"
sampling_design4 = "max_n ~  habitat * study_locations + (1|study_locations:bait_type:op_code) "

# # fit mixed model
model <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all)

model2 <<-
  afex::mixed(formula = sampling_design2,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all)

model3 <<-
  afex::mixed(formula = sampling_design3,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all)

model4 <<-
  afex::mixed(formula = sampling_design4,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all)

# show
model
model2
model3
model4

model = model4

# view anova table
anova(model)

# visualize anova(model)
emmip(model, 
      study_locations ~ habitat,    # type = "response" for back transformed values
      cov.reduce = range) +
  # geom_vline(xintercept=mean(data_all_$primer_x),
  #            color = "grey",
  #            linetype = "dashed") +
  # geom_text(aes(x = mean(data_all_$primer_x),
  #               y = -2,
  #               label = "mean primer_x"),
  #           color = "grey") +
  theme_myfigs +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")

#### max_n: Conduct A priori contrast tests for differences among sites ####

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


#### max_n: Group Sites Based on Model Results ####

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

groupings_model_fixed       # cld messes up back transformation, this takes values from emmeans and groupings from cld

#### max_n: Visualize Estimated Marginal Means Output With Group Categories ####

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
  # geom_point(data = data,
  #            aes(x = location,
  #                y = !!response_var,
  #                color = location
  #            ),
  #            position = position_dodge(width = 0.9),
  #            # color = "grey70",
#            # shape = 1,
#            size = 1)
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
  theme_myfigs +
  # ylim(ymin, 
  #      ymax) +
  labs(x = "Depth",
       y = "Mean Max_N") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) 

p






#### max_n bait: Make Visualization of Hypothesis Test ####
data_all %>%
  group_by(bait_type, 
           habitat,
           op_code) %>%
  dplyr::summarize(mean_max_n = mean(max_n)) %>%
  dplyr::summarize(mean_mean_max_n = mean(mean_max_n),
                   se_mean_max_n = sd(mean_max_n)/sqrt(n())) %>%
  ggplot(aes(x = bait_type,
             y = mean_mean_max_n,
             fill = habitat))+
  geom_bar(position = "dodge", 
           stat = "identity") +
  xlab("Study Locations") +
  ylab("Mean MaxN per BRUV Deployment") +
  labs(title = "Mean MaxN at TRNP vs. Cagayancillo",
       fill = "Habitat") +
  theme_classic() +
  scale_fill_manual(values = habitatcolors, 
                    labels = c("Shallow Reef",
                               "Mesophotic Reef")) +
  geom_errorbar(aes(ymax = mean_mean_max_n + se_mean_max_n,
                    ymin = mean_mean_max_n - se_mean_max_n), 
                position = "dodge") 

data_all_bait <-
  data_all %>%
  filter(bait_type != "Black Jack/Bluefin Trevally",
         bait_type != "Skipjack Tuna")

#### max_n bait: Mixed Effects Hypothesis Test ####

#Here we use the visayan deer data set to demonstrate a mxed model with both fixed and randome factor.scope(

## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(max_n) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"


alpha_sig = 0.05

# we start with the loci subjected to 11 primer concentrations (we removed loci with no max_n to simplify)
sampling_design = "max_n ~  habitat * study_locations + (1 | study_locations:bait_type) + (1|study_locations:bait_type:op_code)"

# # fit mixed model
model <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all_bait)

model
anova(model)
summary(model)

# visualize summary(model)
emmip(model, 
      study_locations ~ habitat,    # type = "response" for back transformed values
      cov.reduce = range) +
  # geom_vline(xintercept=mean(data_all_$primer_x),
  #            color = "grey",
  #            linetype = "dashed") +
  # geom_text(aes(x = mean(data_all_$primer_x),
  #               y = -2,
  #               label = "mean primer_x"),
  #           color = "grey") +
  theme_myfigs +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")

#### max_n bait: Conduct A priori contrast tests for differences among sites ####

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


#### max_n bait: Group Sites Based on Model Results ####

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

groupings_model_fixed       # cld messes up back transformation, this takes values from emmeans and groupings from cld

#### max_n bait: Visualize Estimated Marginal Means Output With Group Categories ####

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
  # geom_point(data = data,
  #            aes(x = location,
  #                y = !!response_var,
  #                color = location
  #            ),
  #            position = position_dodge(width = 0.9),
  #            # color = "grey70",
#            # shape = 1,
#            size = 1)
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
  geom_point(data = data_all,
             aes(x = study_locations,
                 y = !!response_var,
                 shape = habitat),
             position=position_jitterdodge(),
             size = 3,
             color = "darkgrey",
             inherit.aes = FALSE) +
  theme_myfigs +
  # ylim(ymin, 
  #      ymax) +
  labs(x = "Depth",
       y = "Mean Max_N") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) +
  scale_y_continuous(trans='log10')


p


#### max_n bait 2: Make Visualization of Hypothesis Test ####
data_all %>%
  group_by(bait_type, 
           habitat,
           op_code) %>%
  dplyr::summarize(mean_max_n = mean(max_n)) %>%
  dplyr::summarize(mean_mean_max_n = mean(mean_max_n),
                   se_mean_max_n = sd(mean_max_n)/sqrt(n())) %>%
  ggplot(aes(x = bait_type,
             y = mean_mean_max_n,
             fill = habitat))+
  geom_bar(position = "dodge", 
           stat = "identity") +
  xlab("Study Locations") +
  ylab("Mean MaxN per BRUV Deployment") +
  labs(title = "Mean MaxN at TRNP vs. Cagayancillo",
       fill = "Habitat") +
  theme_classic() +
  scale_fill_manual(values = habitatcolors, 
                    labels = c("Shallow Reef",
                               "Mesophotic Reef")) +
  geom_errorbar(aes(ymax = mean_mean_max_n + se_mean_max_n,
                    ymin = mean_mean_max_n - se_mean_max_n), 
                position = "dodge") 

data_all_bait <-
  data_all %>%
  filter(bait_type != "Black Jack/Bluefin Trevally",
         bait_type != "Skipjack Tuna")

#### max_n bait 2: Mixed Effects Hypothesis Test ####

#Here we use the visayan deer data set to demonstrate a mxed model with both fixed and randome factor.scope(

## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(max_n) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"


alpha_sig = 0.05

# we start with the loci subjected to 11 primer concentrations (we removed loci with no max_n to simplify)
sampling_design = "max_n ~  habitat * study_locations + (1 | study_locations:bait_type) + (1|study_locations:bait_type:op_code)"

sampling_design = "max_n ~  habitat + bait_type + (1|bait_type:op_code)"
sampling_design = "max_n ~  bait_type * habitat + (1|bait_type:op_code)"

# # fit mixed model
model <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all_bait)

model
anova(model)
summary(model)

# visualize summary(model)
emmip(model, 
      bait_type ~ habitat,    # type = "response" for back transformed values
      cov.reduce = range) +
  # geom_vline(xintercept=mean(data_all_$primer_x),
  #            color = "grey",
  #            linetype = "dashed") +
  # geom_text(aes(x = mean(data_all_$primer_x),
  #               y = -2,
  #               label = "mean primer_x"),
  #           color = "grey") +
  theme_myfigs +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")

#### max_n bait 2: Conduct A priori contrast tests for differences among sites ####

# now we move on to finish the hypothesis testing.  Are there differences between the sites?
# estimated marginal means 

emmeans_model <<-
  emmeans(model,
          ~ habitat * bait_type,
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


#### max_n bait 2: Group Sites Based on Model Results ####

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
                   "bait_type")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed       # cld messes up back transformation, this takes values from emmeans and groupings from cld

#### max_n bait 2: Visualize Estimated Marginal Means Output With Group Categories ####

p <- 
  groupings_model_fixed %>%
  ggplot(aes(x=bait_type,
             y=response,
             fill = habitat)) +
  geom_col(position = "dodge",
           color = "black") +
  # scale_fill_manual(values = c("lightgrey",
  #                              "white"),
  #                   labels = c('Pre-Screen', 
  #                              'Post-Screen')) +
  # geom_point(data = data,
  #            aes(x = location,
  #                y = !!response_var,
  #                color = location
  #            ),
  #            position = position_dodge(width = 0.9),
  #            # color = "grey70",
#            # shape = 1,
#            size = 1)
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
  theme_myfigs +
  # ylim(ymin, 
  #      ymax) +
  labs(x = "Depth",
       y = "Mean Max_N") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) 

p







#### sum_max_n: Make Visualization of Hypothesis Test ####
response_var = quo(sum_max_n)

habitatcolors <- c("#F08080","#6FAFC6")
habitat(habitatcolors) <- c("Shallow Reef", "Deep Reef")
data_all_summaxn <- 
  data_all %>%
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
  theme_myfigs +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")

#### sum_max_n: Conduct A priori contrast tests for differences among sites ####

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

groupings_model_fixed       # cld messes up back transformation, this takes values from emmeans and groupings from cld

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
  geom_text(aes(label=group),
            position = position_dodge(width=0.9),
            vjust = -0.5,
            hjust = -0.15,
            size = 8 / (14/5)) +  # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  geom_point(data = data_all_summaxn,
              aes(x = study_locations,
                  y = !!response_var,
                  shape = habitat),
              position=position_jitterdodge(),
              size = 3,
              color = "darkgrey",
              inherit.aes = FALSE) +
  scale_y_continuous(trans='log10') +
  theme_myfigs +
  # ylim(ymin, 
  #      ymax) +
  labs(title = sampling_design,
       subtitle = "Distibution Family = Poisson",
       x = "Depth",
       y = "Mean Sum_Max_N") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) 

p





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
                                             "Cagayancillo")),
         habitat = factor(habitat),
         bait_type = factor(bait_type),
         site = factor(site),
         survey_area = factor(survey_area),
         habitat_mpa = str_c(habitat,
                             study_locations,
                             sep = " "))

attach(data_vegan.env)
#### mean_chao_s: Make Visualization of Hypothesis Test ####
pool <- 
  estimateR(x = data_vegan) %>%
  t() %>%
  as_tibble()

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

#### mean_chao_s: Mixed Effects Hypothesis Test ####
data_all_chao_s <- 
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
                   mean_s = mean(s_obs))

## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(mean_chao_s) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "gamma"


alpha_sig = 0.05

# we start with the loci subjected to 11 primer concentrations (we removed loci with no sum_max_n to simplify)


sampling_design = "mean_chao_s ~  habitat * study_locations + (1|study_locations:bait_type)"

# # fit mixed model
model <<-
  afex::mixed(formula = sampling_design,
              family = distribution_family,
              method = "LRT",
              sig_symbols = rep("", 4),
              # all_fit = TRUE,
              data = data_all_chao_s)

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
  theme_myfigs +
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
groupings_model_Sr <<-
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
              dplyr::select(-rate:-asymp.UCL),
            # by = c(str_replace(fixed_vars,
            #                    "[\\+\\*]",
            #                    '" , "'))) %>%
            by = c("habitat",
                   "study_locations")) %>%
  dplyr::rename(response = 3)

groupings_model_fixed_sr 

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
  # geom_point(data = data,
  #            aes(x = location,
  #                y = !!response_var,
  #                color = location
  #            ),
  #            position = position_dodge(width = 0.9),
  #            # color = "grey70",
#            # shape = 1,
#            size = 1)
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
  theme_myfigs +
  # ylim(ymin, 
  #      ymax) +
  labs(x = "Depth",
       y = "Mean Sum_Max_N") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) 

p_sr
#### sum_max_n bait: Make Visualization of Hypothesis Test ####
data_all %>%
  group_by(bait_type, 
           habitat,
           study_locations,
           op_code) %>%
  dplyr::summarize(sum_max_n = sum(max_n)) %>%
  dplyr::summarize(mean_sum_max_n = mean(sum_max_n),
                   se_sum_max_n = sd(sum_max_n)/sqrt(n())) %>%
  ggplot(aes(x = bait_type,
             y = mean_sum_max_n,
             fill = habitat))+
  geom_bar(position = "dodge", 
           stat = "identity") +
  xlab("Study Locations") +
  ylab("Mean MaxN per BRUV Deployment") +
  labs(title = "Mean MaxN at TRNP vs. Cagayancillo",
       fill = "Habitat") +
  theme_classic() +
  scale_fill_manual(values = habitatcolors, 
                    labels = c("Shallow Reef",
                               "Mesophotic Reef")) +
  geom_errorbar(aes(ymax = mean_sum_max_n + se_sum_max_n,
                    ymin = mean_sum_max_n - se_sum_max_n), 
                position = "dodge") +
  facet_grid(. ~ study_locations,
             scales = "free_x")



#### sum_max_n bait: Mixed Effects Hypothesis Test ####

data_all_summaxn_bait <- 
  data_all %>%
  group_by(op_code,
           study_locations,
           habitat,
           bait_type) %>%
  dplyr::summarize(sum_max_n = sum(max_n)) %>%
  filter(bait_type != "Black Jack/Bluefin Trevally",
         bait_type != "Skipjack Tuna")

## Enter Information About Your Data for A Hypothesis Test ##

# define your response variable, here it is binomial
response_var = quo(sum_max_n) # quo() allows column names to be put into variables 

# enter the distribution family for your response variable
distribution_family = "poisson"


alpha_sig = 0.05

# we start with the loci subjected to 11 primer concentrations (we removed loci with no sum_max_n to simplify)

sampling_design = "sum_max_n ~  habitat * study_locations + study_locations:bait_type "

# # fit mixed model
model <<- 
  glm(formula = sampling_design, 
      family = distribution_family,
      data = data_all_summaxn)

model
anova(model)
summary(model)
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
  theme_myfigs +
  labs(title = "Visualization of `summary(model)`",
       subtitle = "",
       y = "Linear Prediciton",
       x = "MPA")

#### sum_max_n bait: Conduct A priori contrast tests for differences among sites ####

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


#### sum_max_n bait: Group Sites Based on Model Results ####

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

groupings_model_fixed       # cld messes up back transformation, this takes values from emmeans and groupings from cld

#### sum_max_n bait: Visualize Estimated Marginal Means Output With Group Categories ####

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
  # geom_point(data = data,
  #            aes(x = location,
  #                y = !!response_var,
  #                color = location
  #            ),
  #            position = position_dodge(width = 0.9),
  #            # color = "grey70",
#            # shape = 1,
#            size = 1)
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
  theme_myfigs +
  # ylim(ymin, 
  #      ymax) +
  labs(x = "Depth",
       y = "Mean Sum_Max_N") +
  theme(legend.position=c(0.33,0.8),  
        legend.title=element_blank()) 

p



