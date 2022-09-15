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
require(maps)
require(viridis)
theme_set(
  theme_void()
)
# install.packages("vegan")
library(vegan)
library(devtools)
#devtools::install_github("jfq3/ggordiplots", dependencies = TRUE)
library(ggordiplots)
library(ggbiplot)
library(ggvegan)
library(ggpubr)

library(ggrepel)

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
           .keep_all = TRUE)

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
  select(-op_code:-studylocation_habitat) %>%
  select_if(colSums(.) != 0)


#attach(data_vegan_TRNP.env)

##Create data_vegan and data_vegan.env for Cagayancillo
data_vegan_CAG.env <- data_vegan.env %>%
  filter(study_locations == "CAGAYANCILLO")

data_vegan_CAG <- bind_cols(data_vegan, data_vegan.env) %>%
  filter(study_locations == "CAGAYANCILLO") %>%
  select(-op_code:-studylocation_habitat) %>%
  select_if(colSums(.) != 0)


#attach(data_vegan_CAG.env)

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

#### Visualize Counts of Top 10 Most Abundant Species by MPA and Depth ####
abundant_species <- data_vegan %>%
                    clean_names() %>%
                    bind_cols(data_vegan.env) %>% 
                    pivot_longer(carangidae_carangoides_coeruleopinnatus:galeomorphii_sphyrna_lewini,
                                 names_to = "taxon") %>%
                    group_by(taxon,
                            habitat,
                            study_locations) %>%
                    dplyr::summarize(taxon_sum = sum(value)) %>% 
                    filter(taxon_sum >= 3) %>% 
                    ##visualize
                    ggplot(aes(x = reorder(taxon, - taxon_sum, sum),
                               y = taxon_sum)) +
                    geom_bar(stat = "identity") +
                    theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Species") +
  ylab("Species Count") +
  scale_y_continuous(trans='log10') +
  facet_grid(habitat ~ study_locations)

                    
abundant_species
ggsave("MostAbundantSpecies.png",
       abundant_species,
       height = 7, 
       width = 7, 
       units = "in")
                  
#### ORDINATION: Non-metric multidimensional scaling (NMDS) ####

ord <- 
  metaMDS(data_vegan %>%
            filter(row_number() != 17), # CAG_017 very different from all other bruvs
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500) 
# View(ord)
ord

# use envfit to generate species loading vectors
# https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2
ord_species_vectors <- 
  envfit(ord$points, 
         data_vegan %>%
           filter(row_number() != 17), 
         perm=1000)

ord_species_vectors

# convert envfit output to ggplot tibble
ggord_species_vectors <- 
  as.data.frame(scores(ord_species_vectors, 
                       display = "vectors")) %>%
  clean_names() %>%
  dplyr::rename(nmds1 = mds1,
                nmds2 = mds2) %>%
  mutate(taxon = rownames(.),
         # convert taxon to plotmath format for figure text/labels
         # https://stackoverflow.com/questions/41528953/how-do-i-include-italic-text-in-geom-text-repel-or-geom-text-labels-for-ggplot
         # https://stackoverflow.com/questions/18237134/line-break-in-expression
         taxon = str_replace(taxon,
                             "^",
                             "atop("),
         taxon = str_replace(taxon,
                             "_",
                             ",~italic('"),
         taxon = str_replace(taxon,
                             "_",
                             " "),
         taxon = str_replace(taxon,
                             "$",
                             "'))"),
          pvals = ord_species_vectors$vectors$pvals) %>%
  # separate(taxon,
  #          into=c("family",
  #                 "species"),
  #          sep = "-") %>%
  filter(pvals <= 0.05)

ggord_species_vectors

# convert metaMDS output to ggplot tibble
ggord <- 
  fortify(ord) %>% 
  tibble() %>% 
  clean_names() %>%
  filter(score == "sites") %>% 
  bind_cols(data_vegan.env %>%
              filter(op_code != "CAG_017") # very different from all other bruvs
            ) %>%
  clean_names()

ggord

#### nMDS Plot of Fish Assemblage at Different Shallow and Deep Reefs at Cagayancillo and Tubbataha ####

habitatcolors <- c("#F08080",
                   "#6FAFC6")
habitatlabels <- c("Shallow Reef",
                   "Mesophotic Reef")

studylocationcolors <- c("#C97CD5",
                         "#79CE7A")
studylocationlabels <- c("CAGAYANCILLO", 
                         "TRNP")
# View(ord)
# View(ggord)
vector_scale_factor = 3.5  #will multiply the vector length by this value to improve clarity of plot
ggord_plot <- 
  ggord %>%
  # filter(label != 17) %>%
  ggplot(aes(x = nmds1,
             y = nmds2,
             color = habitat,
             shape = study_locations)) +
  # scale_x_continuous(limits = c(NA,2)) +
  scale_y_continuous(limits = c(NA,2)) +
    coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = ggord_species_vectors,
               aes(x = 0, 
                   xend = nmds1*vector_scale_factor, 
                   y = 0, 
                   yend = nmds2*vector_scale_factor),
               arrow = arrow(length = unit(.25,
                                           "cm")),
               color = "grey",
               inherit.aes = FALSE) +
  
  geom_point(size = 5) +
  scale_color_manual(values = habitatcolors,
                     labels = habitatlabels)+
  scale_shape_manual(values = c(16,2)) +
  stat_ellipse(aes(group = studylocation_habitat,
                   lty=factor(study_locations))) +
  scale_linetype_manual(values=c(1,2,1,2)) +
  # geom_label(data = ggord_species_vectors,
  #            aes(x = nmds1*vector_scale_factor,
  #                y = nmds2*vector_scale_factor,
  #                label = taxon),
  #            # label.padding = 0,
  #            size = 3,
  #            color = "grey20",
  #            inherit.aes = FALSE) +
  geom_text_repel(data = ggord_species_vectors,
                  parse = TRUE,
                  aes(x = nmds1*vector_scale_factor,
                      y = nmds2*vector_scale_factor,
                      label = taxon),
                  # label.padding = 0,
                  size = 3,
                  color = "grey30",
                  inherit.aes = FALSE) +
  
  theme_classic() +
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  labs(color = "Habitat", 
       shape = "Study Locations", 
       linetype = "Study Locations",
       title = "NMDS: Fish Assemblages") 

ggord_plot

ggsave("NMDSfishassemblageversion3.png",
       ggord_plot,
       height = 5,
       width = 7,
       units = "in")
save_plot("NMDSfishassemblageversion2groupings.png")

ggord_plot_2_3 <- 
  ggord %>%
  # filter(label != 17) %>%
  ggplot(aes(x = nmds2,
             y = nmds3,
             color = habitat,
             shape = study_locations)) +
  # scale_x_continuous(limits = c(-3,3)) +
  geom_point(size = 5) +
  stat_ellipse(aes(group = studylocation_habitat)) +
  theme_classic() +
  xlab("NMDS 2") +
  ylab("NMDS 3") +
  labs(color = "Habitat", 
       shape = "Study Locations", 
       title = "NMDS Plots of Fish Assemblage") + 
  scale_color_manual(values = habitatcolors,
                     labels = habitatlabels)


# ggord %>%
#   ggbiplot(ellipse = TRUE)

# ggord %>%
#   ggbiplot()
# 
# gg_ordiplot(ord, groups = data_vegan.env$habitat_mpa, pt.size = 3)

#### SIMPER ####
##Overall Influential Species 
data.simper <- data_vegan %>%
  simper()

summary(data.simper)

##Overall Influential Species of Shallow vs. Deep Reefs
data.simper.habitat <- data_vegan %>%
  simper(group = habitat)
summary(data.simper.habitat)


##Overall Influential Species of TRNP vs. Cagayancillo
data.simper.study.locations <- data_vegan %>%
  simper(group = study_locations)
summary(data.simper.study.locations)


##Influential Species of TRNP

data.simper.TRNP <- data_vegan_TRNP %>%
  simper()

summary(data.simper.TRNP)


##Influential Species of TRNP shallow and deep reefs
data.simper.TRNP.habitat <- data_vegan_TRNP %>%
  simper(group = habitat)

summary(data.simper.TRNP.habitat)


##Influential Species of CAG
data.simper.CAG <- data_vegan_CAG %>%
  simper()

summary(data.simper.CAG)

##Influential Species of Cagayancillo shallow and deep reefs
data.simper.CAG.habitat <- data_vegan_CAG %>%
  simper(group = habitat)

summary(data.simper.CAG.habitat)

##Influential Species of Shallow Reefs
data.simper.shallow.reefs <- data_vegan_shallow %>%
  simper()

summary(data.simper.shallow.reefs)

##Influential Species of Shallow Reefs between TRNP and Cagayancillo
data.simper.shallow.reefs.study.locations <- data_vegan_shallow %>%
  simper(group = study_locations)

summary(data.simper.shallow.reefs.study.locations)

##Influential Species of Deep Reefs
data.simper.deep.reefs <- data_vegan_deep %>%
  simper()

summary(data.simper.deep.reefs)

##Influential Species of Deep Reefs between TRNP and Cagayancillo
data.simper.deep.reefs.study.locations <- data_vegan_deep %>%
  simper(group = study_locations)

summary(data.simper.deep.reefs.study.locations)

#### Varpart ####
varpart.p <- varpart(X = vegdist(data_vegan), 
        ~ habitat + study_locations,
        data = data_vegan.env)
       # distance = "bray")

#### nMDS of Fish Assemblage at Cagayancillo and Tubbataha w/ Bait Type ####
ggord %>%
  filter (label != 17) %>%
  ggplot(aes(x = nmds1,
             y= nmds2,
             color = bait_type,
             shape = site_code)) +
  geom_point(size = 5) +
  scale_x_continuous(limits = c(-3,3)) +
  # geom_text(aes(label = label), color = "black") +
  theme_classic() +
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  labs(color = "Bait Type", shape = "Study Locations", title = "NMDS Plots of Fish Assemblage") 
save_plot("NMDSfishassemblagebaittype.png")

gg_ordiplot(ord, groups = data_vegan.env$bait_type, pt.size = 3)
# scale_color_manual(values = habitatcolors)
#### vegan::specpool - Extrapolated Species Richness in a Species Pool Based on Incidence (Presence Absence) ####

view(data_vegan.env)
pool <- with(data_vegan.env, specpool(x = data_vegan, 
                                      pool = habitat_mpa,
                                      smallsample = TRUE))

View(pool)

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

p<- vegan::estaccumR(data_vegan, permutations = 999)
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

data_estaccumR_plot

data_estaccumR_plots <- data_estaccumR_plot %>%
  ggplot(aes(x=N,
             y=chao_mean)) +
  geom_ribbon(aes(ymin=chao_ci_lower,
                  ymax=chao_ci_upper),
              fill = "#ebe8f3") +
  geom_line() +
  xlab("Sample Size") +
  ylab("Mean Chao1 Species Richness") +
  theme_classic()

data_estaccumR_plots
ggsave("SpeciesRarefactionCurveAbundance.png",
          data_estaccumR_plots)

##Old Plot
# p.plot<-plot(p, 
#              display = c("chao"))
# p.plot
save_plot("SpeciesRichnessRarefactionCurveAbundance.png", 
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

p_CAG_shallow <- estaccumR(data_vegan_CAG_shallow, permutations = 999)

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

p_CAG_deep <- estaccumR(data_vegan_CAG_deep, permutations = 999)

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

p_TUB_shallow <- estaccumR(data_vegan_TUB_shallow, permutations = 999)

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

p_TUB_deep <- estaccumR(data_vegan_TUB_deep, permutations = 999)
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

dev.off()
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

  

p_Serranidae <- estaccumR(data_vegan_Serranidae, permutations = 999)
View(p_Serranidae)

p_Serranidae_plot <-
  p_Serranidae$S %>%
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
p_Serranidae_plots <- p_Serranidae_plot %>%
  ggplot(aes(x=N,
             y=chao_mean)) +
  geom_ribbon(aes(ymin=chao_ci_lower,
                  ymax=chao_ci_upper),
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

p_Lutjanidae <- estaccumR(data_vegan_Lutjanidae, permutations = 999)

p_Lutjanidae_plot <-
  p_Lutjanidae$S %>%
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
p_Lutjanidae_plots <- p_Lutjanidae_plot %>%
  ggplot(aes(x=N,
             y=chao_mean)) +
  geom_ribbon(aes(ymin=chao_ci_lower,
                  ymax=chao_ci_upper),
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

p_Lethrinidae <- estaccumR(data_vegan_Lethrinidae, permutations = 999)
p_Lethrinidae_plot <-
  p_Lethrinidae$S %>%
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
p_Lethrinidae_plots <- p_Lethrinidae_plot %>%
  ggplot(aes(x=N,
             y=chao_mean)) +
  geom_ribbon(aes(ymin=chao_ci_lower,
                  ymax=chao_ci_upper),
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

p_Carangidae <- estaccumR(data_vegan_Carangidae, permutations = 999)
p_Carangidae_plot <-
  p_Carangidae$S %>%
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
p_Carangidae_plots <- p_Carangidae_plot %>%
  ggplot(aes(x=N,
             y=chao_mean)) +
  geom_ribbon(aes(ymin=chao_ci_lower,
                  ymax=chao_ci_upper),
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
p_groupings_plot <- ggarrange(p_Serranidae_plots,
                              p_Lutjanidae_plots,
                              p_Lethrinidae_plots,
                              p_Carangidae_plots,
                              p_Galeomorphii_plots,
                                      ncol = 2,
                                      nrow = 3)
ggsave("SpeciesRichnessRarefactionCurveGroupings.png", 
       p_groupings_plot, 
       height = 11,
       width = 8.5,
       units = "in")

## Species 
# to make plot w ggplot, see https://stackoverflow.com/questions/52652195/convert-rarefaction-plots-from-vegan-to-ggplot2-in-r

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

p_Galeomorphii_pa <- poolaccum(data_vegan_Galeomorphii, permutations = 999)

p_Galeomorphii_plot_pa <- plot(p_Galeomorphii_pa,
                            display = c("chao"),
                            main = "Galeomorphii")

p_Galeomorphii_plot_pa

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

#### ORDINATION: Fitting Environmental Variables ####

# # Let us test for an effect of site and depth on the NMDS
# 
# ord.fit <- 
#   envfit(ord ~ depth_m + site + bait_type, 
#          data=data_vegan.env, 
#          perm=999,
#          na.rm = TRUE)
# ord.fit
# plot(ord, dis="site")
# ordiellipse(ord, site, col=1:4, kind = "ehull", lwd=3)
# plot(ord.fit)
# 
# # Plotting fitted surface of continuous variables (depth_m) on ordination plot
# ordisurf(ord, depth_m, add=TRUE)
# 
# 
# #### CONSTRAINED ORDINATION ####
# 
# # constrained or “canonical” correspondence analysis (function cca)
# ord <- cca(data_vegan ~ depth_m + site, 
#            data=data_vegan.env)
# ord
# plot(ord, dis="site")
# points(ord, disp="site", pch=21, col=1:2, bg="yellow", cex=1.3)
# ordiellipse(ord, site, col=1:4, kind = "ehull", lwd=3)
# 
# # Significance tests of constraints
# anova(ord)
# anova(ord, by="term", permutations=999)
# anova(ord, by="mar", permutations=999)
# anova(ord, by="axis", permutations=999)
# 
# 
# #### Conditioned or partial ordination ####
# 
# ord <- cca(data_vegan ~ depth_m + site + Condition(bait_type), 
#            data=data_vegan.env)
# anova(ord, by="term", permutations=999)


