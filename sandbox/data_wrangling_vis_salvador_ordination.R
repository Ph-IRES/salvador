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
                             sep = " "))

View(data_vegan.env)

# and now we "attach" the metadata to the data

attach(data_vegan.env)

#### ORDINATION: Detrended correspondence analysis ####

# ord <- decorana(data_vegan)
# ord
# summary(ord)
# #boring plot
# plot(ord)
# #fancier plot
# plot(ord, type = "n")
# points(ord, display = "sites", cex = 0.8, pch=21, col="black", bg="yellow")
# text(ord, display = "spec", cex=0.7, col="red")
# #fanciest plot
# plot(ord, disp="sites", type="n")
# ordihull(ord, habitat, col=1:2, lwd=3)
# ordiellipse(ord, habitat, col=1:2, kind = "ehull", lwd=3)
# ordiellipse(ord, habitat, col=1:2, draw="polygon")
# points(ord, disp="sites", pch=21, col=1:2, bg="yellow", cex=1.3)
# ordispider(ord, habitat, col=1:2, label = TRUE)


#### ORDINATION: Non-metric multidimensional scaling (NMDS) ####

# ord
# summary(ord)
# #fanciest plot
# plot(ord, disp="sites", type="n")
# ordihull(ord, habitat, col=1:2, lwd=3)
# ordiellipse(ord, habitat, col=1:2, kind = "ehull", lwd=3)
# ordiellipse(ord, habitat, col=1:2, draw="polygon")
# points(ord, disp="sites", pch=21, col=1:2, bg="yellow", cex=1.3)
# ordispider(ord, habitat, col=1:2, label = TRUE)
ord <- metaMDS(data_vegan)
ggord <- 
  fortify(ord) %>% 
  tibble() %>% 
  clean_names() %>%
  filter(score == "sites") %>% 
  bind_cols(tibble(data_vegan.env)) %>% 
  clean_names()

habitatcolors <- c("#6FAFC6", "#F08080")
habitat(habitatcolors) <- ("Shallow Reef", "Deep Reef")

ggord %>%
  ggplot(aes(x = nmds1,
             y= nmds2,
             color = habitat,
             shape = site_code)) +
  geom_point(size = 5) +
  theme_classic() +
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  labs(color = "Habitat", shape = "Study Locations", title = "NMDS Plots of Fish Assemblage") +
  scale_color_manual(values = habitatcolors)
save_plot("NMDSfishassemblage.png")



#### ORDINATION: Fitting Environmental Variables ####

# Let us test for an effect of site and depth on the NMDS

ord.fit <- 
  envfit(ord ~ depth_m + site + bait_type, 
         data=data_vegan.env, 
         perm=999,
         na.rm = TRUE)
ord.fit
plot(ord, dis="site")
ordiellipse(ord, site, col=1:4, kind = "ehull", lwd=3)
plot(ord.fit)

# Plotting fitted surface of continuous variables (depth_m) on ordination plot
ordisurf(ord, depth_m, add=TRUE)


#### CONSTRAINED ORDINATION ####

# constrained or “canonical” correspondence analysis (function cca)
ord <- cca(data_vegan ~ depth_m + site, 
           data=data_vegan.env)
ord
plot(ord, dis="site")
points(ord, disp="site", pch=21, col=1:2, bg="yellow", cex=1.3)
ordiellipse(ord, site, col=1:4, kind = "ehull", lwd=3)

# Significance tests of constraints
anova(ord)
anova(ord, by="term", permutations=999)
anova(ord, by="mar", permutations=999)
anova(ord, by="axis", permutations=999)


#### Conditioned or partial ordination ####

ord <- cca(data_vegan ~ depth_m + site + Condition(bait_type), 
           data=data_vegan.env)
anova(ord, by="term", permutations=999)


