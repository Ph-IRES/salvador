#### INITIALIZATION ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(janitor)
library(readxl)
# install.packages("maps")
# install.packages("viridis")
require(maps)
require(viridis)
library(sjPlot)
library(vegan)
theme_set(
  theme_void()
)
library(purrr)

#### USER DEFINED VARIABLES ####

inFilePath1 = "./WorkingData_CLEANED_TUB,CAG.xlsx"
inFilePath2 = "./PHIRES_MetaData.xlsx"

# outFilePath = "./data_combined.tsv"

#### READ IN DATA & CURATE ####

View(data)

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
  read_excel(inFilePath2,
             na="NA") %>%
  clean_names() %>%
  dplyr::rename(bait_weight_grams = weight_grams)

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
         everything()) %>%
  mutate(study_locations = case_when(
    site == "Cawili" ~ "CAGAYANCILLO",
    site == "Calusa" ~ "CAGAYANCILLO",
    site == "Cagayancillo" ~ "CAGAYANCILLO",
    site == "TUBBATAHA" ~ "TRNP"))

data_all_removed_sp <- 
  data_removed_sp %>%
  left_join(metadata,
            by = c("op_code" = "opcode",
                   "depth_m" = "depth_m")) %>%
  # rearrange order of columns, metadata then data
  select(op_code,
         site:long_e,
         depth_m,
         time_in:bait_weight_grams,
         everything()) %>%
  mutate(study_locations = case_when(
    site == "Cawili" ~ "CAGAYANCILLO",
    site == "Calusa" ~ "CAGAYANCILLO",
    site == "Cagayancillo" ~ "CAGAYANCILLO",
    site == "TUBBATAHA" ~ "TRNP"))

#### CHECK DATA INTEGRITY ####
# isolate duplicated rows, if no dups, then tibbles should be empty
data_all %>%
  group_by(op_code,
           taxon) %>%
  filter(n()>1) 


data %>%
  group_by(op_code,
           taxon) %>%
  filter(n()>1) 

# check unique taxon names
data$taxon %>% unique() %>% sort()


#### VISuALIZE METADATA ####

metadata %>%
  ggplot(aes(x=depth_m,
             fill = habitat)) +
  geom_histogram() +
  theme_classic() 
ggsave("histogram_depth-x-habitat.png")

metadata %>%
  ggplot(aes(x=depth_m,
             fill = habitat)) +
  geom_histogram() +
  theme_classic() +
  facet_grid(bait_type ~ bait_weight_grams)
ggsave("histogram_depth-x-habitat-x-bait.png")

metadata %>%
  ggplot(aes(x=habitat,
             y=survey_length_hrs,
             fill = habitat)) +
  # geom_violin() +
  geom_boxplot() +
  theme_classic() 

metadata %>%
  ggplot(aes(y=lat_n,
             x=long_e,
             color = habitat)) +
  geom_point(size = 5) +
  theme_classic() 

#### VISUALIZE ALL DATA ####

data_all %>%
  ggplot(aes(x=habitat,
             y=max_n,
             color = habitat)) +
  geom_boxplot() +
  theme_classic() +
  facet_grid(trophic_groups ~ .)

data_all %>%
  ggplot(aes(x=max_n,
             fill = habitat)) +
  geom_histogram(position="identity") +
  theme_classic() +
  facet_grid(. ~ trophic_groups,
             scales = "free_x")

#### MAP DATA ####
#https://www.datanovia.com/en/blog/how-to-create-a-map-using-ggplot2/

plot1 <- 
  metadata %>%
    ggplot(aes(y=lat_n,
               x=long_e,
               color = habitat)) +
    geom_point(size = 5) +
    theme_classic()

world_map <- map_data("world")

world_map %>%
ggplot(aes(x = long, 
           y = lat, 
           group = group)) +
  geom_polygon(fill="lightgray", 
               color = "red")

map_data("world",
         region = "Philippines") %>%
  ggplot(aes(x = long, 
             y = lat,
             group = group)) +
  geom_polygon(fill="lightgray",
               colour = "black") 

subregion_label_data <- 
  map_data("world",
           region = "Philippines") %>%
  dplyr::group_by(subregion,
                  group) %>%
  dplyr::summarize(long = mean(long), 
                   lat = mean(lat)) %>%
  filter(subregion == "Negros" |
           subregion == "Cebu")

region_label_data <- 
  map_data("world",
           region = "Philippines") %>%
  dplyr::group_by(region) %>%
  dplyr::summarize(long = mean(long), 
                   lat = mean(lat))

studylocationcolors <- c("#C97CD5","#79CE7A")
#labels(studylocationcolors) <- ("CAGAYANCILLO", "TRNP")

# insert in your ggplot commands to control colors of MPAs
# scale_fill_manual(values = c('CAGAYANCILLO' = "lightblue",
#                              'TRNP' = "palegreen")) +
# scale_color_manual(values = c('CAGAYANCILLO' = "lightblue",
#                               'TRNP' = "#00BA38"))
minLat = 5
minLong = 116
maxLat = 20
maxLong = 130
subregions_keep <- 
  map_data("world",
         region = "Philippines") %>%
  filter(long > minLong,
         long < maxLong,
         lat > minLat,
         lat< maxLat) %>%
  distinct(subregion) %>%
  pull()
  

subregions_keep %>%
  purrr::map(~ map_data("world",
                        region = "Philippines") %>%
               filter(subregion == .x)) %>%
mutate(lat = case_when(lat < minLat ~ minLat,
                       lat > maxLat ~ maxLat, TRUE ~ lat),
       long = case_when(long < minLong ~ minLong,
                       long > maxLong ~ maxLong,
                       TRUE ~ long))

map_data("world",
         region = "Philippines") %>%
  ggplot(aes(long,
             lat,
             group=group)) +
  geom_polygon(fill="lightgray",
               color = "black") +
  # geom_text(data = subregion_label_data,
  #           aes(label = subregion),
  #           size = 6,
  #           hjust = 0.5) +
  geom_text(data = region_label_data,
            aes(x = long,
                y= lat,
                label = ""),
            size = 10,
            hjust = 0.5,
            inherit.aes = FALSE) +
  geom_point(data = data_all,
             aes(x = long_e,
                 y = lat_n,
                 color = study_locations),
             inherit.aes = FALSE) +
  theme_classic() + 
  xlab("Longitude East (in degrees)") +
  ylab("Latitude North (in degrees)") +
  labs(title = "Map of BRUV Deployments", 
       color = "Study Locations") +
  scale_color_manual(values = studylocationcolors)
  save_plot("MapofBRUVDeployments.png")
  
#### Map Data ####
 
  minLat = 7
  minLong = 119
  maxLat = 10
  maxLong = 122.5
  
  studylocationcolors <- c("#C97CD5","#79CE7A")
  study_locations(studylocationcolors) <- c("CAGAYANCILLO", "TRNP")

  region_label_data <- 
    map_data("world",
             region = "Philippines") %>%
    dplyr::group_by(region) %>%
    dplyr::summarize(long = mean(long), 
                     lat = mean(lat))
  #Zoomed Out Map
  map_data("world",
           region = "Philippines") %>%
    ggplot(aes(long,
               lat,
               group=group)) +
    geom_polygon(fill="lightgray",
                 color = "black") +
    # geom_text(data = subregion_label_data,
    #           aes(label = subregion),
    #           size = 6,
    #           hjust = 0.5) +
    geom_text(data = region_label_data,
              aes(x = long,
                  y= lat,
                  label = ""),
              size = 10,
              hjust = 0.5,
              inherit.aes = FALSE) +
    geom_point(data = data_all,
               aes(x = long_e,
                   y = lat_n,
                   color = study_locations),
               inherit.aes = FALSE) +
    theme_classic() + 
    xlab("Longitude East (in degrees)") +
    ylab("Latitude North (in degrees)") +
    labs(title = "Map of BRUV Deployments", 
         color = "Study Locations") +
    scale_color_manual(values = studylocationcolors)
  save_plot("MapofBRUVDeployments.png")
  
  #Zoomed in Maps
  
  subregions_keep <-
    map_data("world",
             region = "Philippines") %>%
    filter(long > minLong,
           long < maxLong,
           lat > minLat,
           lat < maxLat) %>%
    distinct(subregion) %>%
    pull()
  
  subregions_keep %>%
    purrr::map_df(~ map_data("world",
                             region = "Philippines") %>%
            filter(subregion == .x)) %>%
 mutate(lat = case_when(lat < minLat ~ minLat,
                        lat > maxLat ~ maxLat,
                        TRUE ~ lat),
        long = case_when(long < minLong ~ minLong,
                            long > maxLong ~ maxLong,
                            TRUE ~ long)) %>%
    ggplot(aes(long,
               lat,
               group=group)) +
    geom_polygon(fill="lightgray") +
                 #color = "black") +
    # geom_text(data = subregion_label_data,
    #           aes(label = subregion),
    #           size = 6,
    #           hjust = 0.5) +
    geom_text(data = region_label_data,
              aes(x = long,
                  y= lat,
                  label = "Sulu Sea"),
              size = 5,
              hjust = 1,
              vjust = 10,
              inherit.aes = FALSE) +
    geom_point(data = data_all,
               aes(x = long_e,
                   y = lat_n,
                   color = study_locations),
               inherit.aes = FALSE) +
    scale_color_manual(values = studylocationcolors) +
    theme_classic() +
    xlab("Longitude East (in degrees)") + 
    ylab("Latitude North (in degrees)") +
    labs(color = "Study Locations", title = "Map of Study Locations")
  save_plot("MapofStudyLocationsZoomedin.png")
  
  data_all %>%
    ggplot(aes(y=lat_n,
               x=long_e,
               color = study_locations)) +
    geom_point(size = 5) +
    theme_classic() +
    scale_color_manual(values = studylocationcolors)+
    xlab("Latitude North (in degrees)") +
    ylab("Longitude East (in degrees)") +
    labs(title = "Map of Study Locations", color = "Study Locations") 
  save_plot("closemapstudylocations.png")
  

# #### Mikaela's Data CleanUp and Modifications ####
# #View(data_all)
# data_all <- data_all %>%
#   mutate(study_locations = case_when(
#     site == "Cawili" ~ "CAGAYANCILLO",
#     site == "Calusa" ~ "CAGAYANCILLO",
#     site == "Cagayancillo" ~ "CAGAYANCILLO",
#     site == "TUBBATAHA" ~ "TRNP")) %>%
#   mutate(family_clean = case_when(
#     family == "Epinephelidae" ~ "Serranidae",
#     TRUE ~ family))
# 
# View(data_all)

#### Mikaela's Data Visualization ####
habitatcolors <- c("#6FAFC6","#F08080")
habitat(habitatcolors) <- c("Shallow Reef", "Deep Reef") # habitat() is not a function
#histogram
data_all %>%
  ggplot(aes(x = max_n)) +
  geom_histogram()+
  facet_grid(study_locations ~ habitat)
save_plot("DistributionofMaxNatdepth.png")

#histogram for the sum of max_n by op_code
data_all %>%
  group_by(op_code,
           study_locations,
           habitat) %>%
  summarize(sum_max_n = sum(max_n)) %>%
  ggplot(aes(x = sum_max_n)) +
  geom_histogram()+
  theme_classic() +
  facet_grid(study_locations ~ habitat)
ggsave("DistributionofMaxNatdepth_summaxn.png")

#Barplot of MaxN per BRUV Station at TRNP and Cagayancillo
View(data_all)
data_compiled <- 
  data_all %>%
  group_by(study_locations, 
           habitat,
           op_code) %>%
  summarize(sum_max_n = sum(max_n)) %>%
  summarize(mean_sum_max_n = mean(sum_max_n),
            se_sum_max_n = sd(sum_max_n)/sqrt(n())) %>%
  ggplot(aes(x = study_locations,
             y = mean_sum_max_n,
             fill = habitat))+
  geom_bar(position = "dodge", 
           stat = "identity") +
  xlab("Study Locations") +
  ylab("Mean MaxN per BRUV Deployment") +
  labs(title = "Mean MaxN at TRNP vs. Cagayancillo",
       fill = "Habitat") +
  theme_classic() +
  scale_fill_manual(values = habitatcolors) +
  geom_errorbar(aes(ymax = mean_sum_max_n + se_sum_max_n,
                    ymin = mean_sum_max_n - se_sum_max_n), 
                position = "dodge")
data_compiled  
save_plot("MeanMaxNatTRNPvs.Cagayancillo.png")

#Barplot of MaxN at TRNP and Cagayancillo, filled by bait_type instead of shallow vs. deep reef
View(data_all)
data_compiled_bait <- 
  data_all %>%
  group_by(study_locations, 
           op_code,
           bait_type) %>%
  summarize(sum_max_n = sum(max_n)) %>%
  summarize(mean_sum_max_n = mean(sum_max_n),
            se_sum_max_n = sd(sum_max_n)/sqrt(n())) %>%
  ggplot(aes(x = study_locations,
             y = mean_sum_max_n,
             fill = bait_type)) +
  geom_bar(position = "dodge", 
           stat = "identity") +
  xlab("Study Locations") +
  ylab("Mean MaxN per BRUV Deployment") +
  labs(title = "Mean MaxN at TRNP vs. Cagayancillo",
       fill = "Bait Type") +
  theme_classic() +
  geom_errorbar(aes(ymax = mean_sum_max_n + se_sum_max_n,
                    ymin = mean_sum_max_n - se_sum_max_n), 
                position = "dodge")
data_compiled_bait  
save_plot("MeanMaxNatTRNPvs.Cagayancilloincludingbait.png")

#Faceted Barplot of MaxN at TRNP and Cagayancillo faceted by family 
data_compiled_faceted <- 
  data_all %>%
  group_by(study_locations, 
           habitat,
           op_code,
           family_clean) %>%
  summarize(sum_max_n = sum(max_n)) %>% 
  group_by(study_locations, 
           habitat,
           family_clean) %>%
  summarize(mean_sum_max_n = mean(sum_max_n),
            se_sum_max_n = sd(sum_max_n)/sqrt(n())) %>% 
  
  ggplot(aes(x = study_locations,
             y = mean_sum_max_n,
             fill = habitat))+
  geom_bar(position = "dodge", 
           stat = "identity") +
  xlab("Study Locations") +
  ylab("Mean MaxN per BRUV Deployment") +
  labs(title = "Mean MaxN at TRNP vs. Cagayancillo",
       fill = "Habitat") +
  theme_classic() +
  scale_fill_manual(values = habitatcolors) +
  geom_errorbar(aes(ymax = mean_sum_max_n + se_sum_max_n,
                    ymin = mean_sum_max_n - se_sum_max_n), 
                position = "dodge") +
  facet_wrap(family_clean ~ .,
             scales = "free_y") +
  theme(strip.text.y.right = element_text(angle = 0)) 
data_compiled_faceted  
ggsave("FacetedMaxN.pdf", 
       data_compiled_faceted, 
       height = 11, 
       width = 8.5, 
       units = "in")

#Faceted Barplot of MaxN at TRNP and Cagayancillo faceted by family, INCLUDING ZEROS
data_compiled_faceted_zeros <- 
  data_all %>%
  group_by(study_locations, 
           habitat,
           op_code,
           family_clean) %>%
  summarize(sum_max_n = sum(max_n)) %>%
  ungroup() %>%
  arrange(op_code) %>%

  #add observations of zero individuals in a family
  mutate(opcode_habitat = str_c(op_code,
                                habitat,
                                sep = ",")) %>%
  select(-op_code,
         -habitat) %>%
  complete(study_locations,
           opcode_habitat,
           family_clean,
           fill = list(sum_max_n =0)) %>%
  separate(opcode_habitat,
           into = c("op_code",
                    "habitat"),
           sep = ",") %>%
  
  group_by(study_locations, 
           habitat,
           family_clean) %>%
  summarize(mean_sum_max_n = mean(sum_max_n),
            se_sum_max_n = sd(sum_max_n)/sqrt(n())) %>%
  ggplot(aes(x = study_locations,
             y = mean_sum_max_n,
             fill = habitat))+
  geom_bar(position = "dodge", 
           stat = "identity") +
  xlab("Study Locations") +
  ylab("Mean MaxN per BRUV Deployment") +
  labs(title = "Mean MaxN at TRNP vs. Cagayancillo w/ Zeros",
       fill = "Habitat") +
  theme_classic() +
  scale_fill_manual(values = habitatcolors) +
  geom_errorbar(aes(ymax = mean_sum_max_n + se_sum_max_n,
                    ymin = mean_sum_max_n - se_sum_max_n), 
                position = "dodge") +
  facet_wrap(family_clean ~ .,
             scales = "free_y") +
  theme(strip.text.y.right = element_text(angle = 0)) 
data_compiled_faceted_zeros  
ggsave("FacetedMaxNwZeros.pdf", 
       data_compiled_faceted_zeros, 
       height = 11, 
       width = 8.5, 
       units = "in")

#Faceted Barplot of MaxN at TRNP and Cagayancillo faceted by bait_type, INCLUDING ZEROS
data_compiled_faceted_zeros_bait <- 
  data_all %>%
  group_by(study_locations, 
           bait_type,
           op_code,
           family_clean) %>%
  summarize(sum_max_n = sum(max_n)) %>%
  ungroup() %>%
  arrange(op_code) %>%
  
  #add observations of zero individuals in a family
  mutate(opcode_bait_type = str_c(op_code,
                                bait_type,
                                sep = ",")) %>%
  select(-op_code,
         -bait_type) %>%
  complete(study_locations,
           opcode_bait_type,
           family_clean,
           fill = list(sum_max_n =0)) %>%
  separate(opcode_bait_type,
           into = c("op_code",
                    "bait_type"),
           sep = ",") %>%
  
  group_by(study_locations, 
           bait_type,
           family_clean) %>%
  summarize(mean_sum_max_n = mean(sum_max_n),
            se_sum_max_n = sd(sum_max_n)/sqrt(n())) %>%
  ggplot(aes(x = study_locations,
             y = mean_sum_max_n,
             fill = bait_type))+
  geom_bar(position = "dodge", 
           stat = "identity") +
  xlab("Study Locations") +
  ylab("Mean MaxN per BRUV Deployment") +
  labs(title = "Mean MaxN at TRNP vs. Cagayancillo w/ Zeros",
       fill = "Bait Type") +
  theme_classic() +
  geom_errorbar(aes(ymax = mean_sum_max_n + se_sum_max_n,
                    ymin = mean_sum_max_n - se_sum_max_n), 
                position = "dodge") +
  facet_wrap(family_clean ~ .,
             scales = "free_y") +
  theme(strip.text.y.right = element_text(angle = 0)) 
data_compiled_faceted_zeros_bait  
dev.off()
ggsave("FacetedMaxNwZeroswBait.pdf", 
       data_compiled_faceted_zeros_bait, 
       height = 11, 
       width = 8.5, 
       units = "in")

#### PREP DATA FOR VEGAN ####

# convert species count data into tibble for vegan ingestion
# each row is a site
# each column is a taxon
# data are counts



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
  dplyr::select(-op_code) %>% view()
 

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
                             sep = " "))

View(data_vegan.env)

# and now we "attach" the metadata to the data

attach(data_vegan.env)

#### PERMANOVA W ADONIS2 ####

# vegan manual - https://cloud.r-project.org/web/packages/vegan/vegan.pdf

# global test of model, differences in species composition with depth and site
adonis2(data_vegan ~ habitat * study_locations,
        data = data_vegan.env,
        by = NULL)

# test for differences in species composition with depth and site by each predictor, this is the default behavior, so `by` is not necessary
adonis2(data_vegan ~ habitat*study_locations,
        data = data_vegan.env,
        by = "terms")

# test for differences in species composition with depth and site by each predictor, marginal effects
adonis2(data_vegan ~ depth_m*study_locations,
        data = data_vegan.env,
        by = "margin")

# dissimilarity indices can be selected, see `vegdist` in the vegan manual for options
adonis2(data_vegan ~ depth_m*study_locations,
        data = data_vegan.env,
        method = "bray")

adonis2(data_vegan ~ depth_m*study_locations,
        data = data_vegan.env,
        method = "euclidean")

# by default, missing data will cause adonis2 to fail, but there are other alternatives
# only non-missing site scores, remove all rows with missing data
adonis2(data_vegan ~ depth_m*study_locations,
        data = data_vegan.env,
        na.action = na.omit)

adonis2(data_vegan ~ depth_m *bait_type,
        data = data_vegan.env,
        na.action = na.omit)

adonis2(data_vegan ~ study_locations*bait_type,
        data = data_vegan.env,
        na.action = na.omit)

# do not remove rows with missing data, but give NA for scores of missing observations or results that cannot be calculated
adonis2(data_vegan ~ depth_m*study_locations,
        data = data_vegan.env,
        na.action = na.exclude)

dis<-vegdist(data_vegan, method="bray")
#creates distance matrix
beta<-betadisper(dis, group = habitat_mpa)
anova(beta)
permutest(beta, pairwise=TRUE, permutations=999)
TukeyHSD(beta)
boxplot(beta)


# constrain permutations within sites, if site is a "block" factor, then this is correct and including site as a factor is incorrect
 adonis2(data_vegan ~ bait_type*habitat,
         data = data_vegan.env,
         strata = site)

####Rarefaction Curve Species Richness
#### vegan::estimateR - Extrapolated Species Richness in a Species Pool Based on Incidence (Abundance) ####

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
           site_code == "TUB" ~ "TUBBATAHA"
         )) %>% 
  pivot_wider(names_from = se_value) %>% 
  rename(sp_richness_est = value,
         estimator = name) %>% 
  filter(estimator != "ace") %>% 
  group_by(study_locations, habitat) %>%
  summarise(mean_chao_s = mean(sp_richness_est),
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
  scale_fill_manual(values = habitatcolors) +
  labs(title = "Species Richness at Cagayancillo vs. Tubbataha", fill = "Habitat") +
  xlab("Study Locations") +
  ylab("Mean Chao Estimate of Species Richness")

save_plot("ChaoEstimateSpeciesRichnessBarplot.png")

#Faceted Barplot of Species Richness
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
           site_code == "TUB" ~ "TUBBATAHA"
         )) %>% 
  pivot_wider(names_from = se_value) %>% 
  rename(sp_richness_est = value,
         estimator = name) %>% 
  filter(estimator != "ace") %>% 
  group_by(study_locations, habitat) %>%
  summarise(mean_chao_s = mean(sp_richness_est),
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
  scale_fill_manual(values = habitatcolors) +
  labs(title = "Species Richness at Cagayancillo vs. Tubbataha", fill = "Habitat") +
  xlab("Study Locations") +
  ylab("Mean Chao Estimate of Species Richness") 

