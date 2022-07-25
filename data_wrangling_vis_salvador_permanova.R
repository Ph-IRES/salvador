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
         trophic_groups = str_to_title(trophic_groups))

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

#### PREP DATA FOR VEGAN ####

# convert species count data into tibble for vegan ingestion
  # each row is a site
  # each column is a taxon
  # data are counts

data_vegan <-
  data %>%
  # make unique taxa
  mutate(taxon = str_c(family,
                       genus,
                       species,
                       sep = "_")) %>%
  # sum all max_n counts for a taxon and op_code
  group_by(taxon,
           op_code) %>%
  summarize(sum_max_n = sum(max_n)) %>%
  ungroup() %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = sum_max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(-op_code)

# convert metadata into tibble for vegan ingestion
# each row is a site
# each column is a dimension of site, such as depth, lat, long, region, etc

data_vegan.env <-
  data_all %>%
  # make unique taxa
  mutate(taxon = str_c(family,
                       genus,
                       species,
                       sep = "_")) %>%
  # sum all max_n counts for a taxon and op_code
  group_by(taxon,
           op_code,
           site,
           survey_area,
           habitat,
           lat_n,
           long_e,
           depth_m,
           bait_type) %>%
  summarize(sum_max_n = sum(max_n)) %>%
  ungroup() %>%
  # convert tibble from long to wide format
  pivot_wider(names_from = taxon,
              values_from = sum_max_n,
              values_fill = 0) %>%
  # sort by op_code
  arrange(op_code) %>%
  # remove the op_code column for vegan
  dplyr::select(op_code:bait_type) %>%
  mutate(site_code = str_remove(op_code,
                                "_.*$"),
         site_code = factor(site_code),
         habitat = factor(habitat),
         bait_type = factor(bait_type),
         site = factor(site),
         survey_area = factor(survey_area))

# and now we "attach" the metadata to the data

attach(data_vegan.env)


#### PERMANOVA W ADONIS2 ####

# vegan manual - https://cloud.r-project.org/web/packages/vegan/vegan.pdf

# global test of model, differences in species composition with depth and site
adonis2(data_vegan ~ depth_m*site,
        data = data_vegan.env,
        by = NULL)

# test for differences in species composition with depth and site by each predictor, this is the default behavior, so `by` is not necessary
adonis2(data_vegan ~ depth_m*site,
        data = data_vegan.env,
        by = "terms")

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
        strata = site)

