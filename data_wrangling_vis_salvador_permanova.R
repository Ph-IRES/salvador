#### INITIALIZATION ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(janitor)
library(readxl)
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
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

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
  mutate(depth_m = case_when(op_code == "CAG_024" ~ 8,
                             TRUE ~ depth_m),
         family = str_to_title(family),
         genus = str_to_title(genus),
         species = str_to_lower(species),
         trophic_groups = str_to_title(trophic_groups),
         family_clean = case_when(
           family == "Epinephelidae" ~ "Serranidae",
           TRUE ~ family),
         taxon = str_c(family_clean,
                       genus,
                       species,
                       sep = "_"))  %>%
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
         family_clean != "Carcharhinidae") %>% view()

metadata <-
  read_excel(inFilePath2,
             na="NA") %>%
  clean_names() %>%
  dplyr::rename(bait_weight_grams = weight_grams) %>%
  mutate(site = str_to_title(site),
         survey_area = str_to_title(survey_area),
         habitat = str_to_title(habitat),
         bait_type = str_to_title(bait_type))

#### COMBINE DATA ####

data_all <-
  data %>%
  left_join(metadata,
            by = c("op_code" = "opcode",
                   "depth_m" = "depth_m")) %>%
  # rearrange order of columns, metadata then data
  select(op_code,
         site:long_e,
         depth_m,
         time_in:bait_weight_grams,
         everything()) 

data_all <- data %>%
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


#### PREP DATA FOR VEGAN ####

# convert species count data into tibble for vegan ingestion
# each row is a site
# each column is a taxon
# data are counts

# data_vegan <-
#   data %>%
#   # make unique taxa
#   mutate(taxon = str_c(family,
#                        genus,
#                        species,
#                        sep = "_")) %>%
#   # sum all max_n counts for a taxon and op_code
#   group_by(taxon,
#            op_code) %>%
#   summarize(sum_max_n = sum(max_n)) %>%
#   ungroup() %>%
#   # convert tibble from long to wide format
#   pivot_wider(names_from = taxon,
#               values_from = sum_max_n,
#               values_fill = 0) %>%
#   # sort by op_code
#   arrange(op_code) %>%
#   # remove the op_code column for vegan
#   dplyr::select(-op_code)

data_vegan <- 
  data %>%
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

# data_vegan.env <-
#   data_all %>%
#   # make unique taxa
#   mutate(taxon = str_c(family,
#                        genus,
#                        species,
#                        sep = "_")) %>%
#   # sum all max_n counts for a taxon and op_code
#   group_by(taxon,
#            op_code,
#            site,
#            survey_area,
#            habitat,
#            lat_n,
#            long_e,
#            depth_m,
#            bait_type) %>%
#   summarize(sum_max_n = sum(max_n)) %>%
#   ungroup() %>%
#   # convert tibble from long to wide format
#   pivot_wider(names_from = taxon,
#               values_from = sum_max_n,
#               values_fill = 0) %>%
#   # sort by op_code
#   arrange(op_code) %>%
#   # remove the op_code column for vegan
#   dplyr::select(op_code:bait_type) %>%
#   mutate(site_code = str_remove(op_code,
#                                 "_.*$"),
#          site_code = factor(site_code),
#          habitat = factor(habitat),
#          bait_type = factor(bait_type),
#          site = factor(site),
#          survey_area = factor(survey_area))

data_vegan.env <-
  data_all %>%
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
         habitat = factor(habitat),
         bait_type = factor(bait_type),
         site = factor(site),
         survey_area = factor(survey_area)) %>%
  mutate(study_locations = case_when(
    site == "Calusa" ~ "Cagayancillo",
    site == "Cawili" ~ "Cagayancillo",
    site == "Cagayancillo" ~ "Cagayancillo",
    site == "Tubbataha" ~ "TRNP"
  )) %>% 
  mutate(study_locations = factor(study_locations,
                                  levels = c("TRNP", 
                                             "Cagayancillo"))) %>% 
  mutate(habitat_mpa = str_c(habitat,
                             study_locations,
                             sep = " ")) %>% view()

# and now we "attach" the metadata to the data

attach(data_vegan.env)


#### PERMANOVA W ADONIS2 ####

# vegan manual - https://cloud.r-project.org/web/packages/vegan/vegan.pdf

# global test of model, differences in species composition with depth and site
#The global test of the whole model is the most powerful test of your hypothesis that you can perform. The result of this test should be the first reported in your results for the test of your model.  If the global test of the model is not significant, then there is no reason to test the individual terms of the model.  In the example here, the global test is significant (see the `Pr(>F)` column in the PERMANOVA table.)
adonis2(data_vegan ~ depth_m*site,
        data = data_vegan.env,
        by = NULL)


# test for differences in species composition with depth and site by each predictor, this is the default behavior, so `by` is not necessary
# Once we have found the model to be significant, we can move on to testing whether each term in the statistical model is non-randomly related to the response variables.

adonis2(data_vegan ~ depth_m*site,
        data = data_vegan.env,
        by = "terms")

# the rest of these examples demonstrate additional functionality. your data and sampling design dicate how your parameterize `adonis2`

# test for differences in species composition with depth and site by each predictor, marginal effects
adonis2(data_vegan ~ depth_m*site,
        data = data_vegan.env,
        by = "margin")

# dissimilarity indices can be selected, see `vegdist` in the vegan manual for options
adonis2(data_vegan ~ depth_m*site,
        data = data_vegan.env,
        method = "bray")

adonis2(data_vegan ~ depth_m*site,
        data = data_vegan.env,
        method = "euclidean")

# by default, missing data will cause adonis2 to fail, but there are other alternatives
# only non-missing site scores, remove all rows with missing data
adonis2(data_vegan ~ depth_m*site,
        data = data_vegan.env,
        na.action = na.omit)

# do not remove rows with missing data, but give NA for scores of missing observations or results that cannot be calculated
adonis2(data_vegan ~ depth_m*site,
        data = data_vegan.env,
        na.action = na.exclude)

# constrain permutations within sites, if site is a "block" factor, then this is correct and including site as a factor is incorrect
adonis2(data_vegan ~ bait_type*habitat,
        data = data_vegan.env,
        permutations = 10000)

##Effects of Study Location (site) and depth category 
adonis2(data_vegan ~ site*habitat,
        data = data_vegan.env,
        na.action = na.exclude)

#### PERMANOVA W ADONIS2 - manuscript version ####

adonis2(data_vegan ~ study_locations * habitat,
        data = data_vegan.env,
        permutations = 10000)

# post hoc test
# installed pairwiseAdonis from instructions on https://github.com/pmartinezarbizu/pairwiseAdonis
pairwise_results <- 
  pairwise.adonis(data_vegan, 
                  interaction(data_vegan.env$study_locations,
                              data_vegan.env$habitat), 
                  p.adjust.m = "BH",
                  perm = 10000)
print(pairwise_results)

# alternative test: this confirms the result with pairwiseAdonis, all sig diff
adonis2(data_vegan ~ interaction(data_vegan.env$study_locations,
                                 data_vegan.env$habitat),
        data = data_vegan.env,
        permutations = 10000,
        by = "onedf")

#### PeRMANOVA w bait type ####

adonis2(data_vegan ~ study_locations * habitat * bait_type,
        data = data_vegan.env,
        permutations = 10000)


# constrain permutations within sites, if site is a "block" factor, then this is correct and including site as a factor is incorrect
adonis2(data_vegan ~ bait_type*habitat,
        data = data_vegan.env,
        permutations = 10000)

pairwise_results <- 
  pairwise.adonis(data_vegan, 
                  interaction(data_vegan.env$bait_type,
                              data_vegan.env$habitat), 
                  p.adjust.m = "BH",
                  perm = 10000)
print(pairwise_results)

write.csv(pairwise_results, file = "posthoctests_baittype_depthcategory.csv")
