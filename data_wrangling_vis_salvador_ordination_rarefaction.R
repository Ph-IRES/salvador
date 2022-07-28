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
  select(op_code,
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

View(data_vegan.env)

# and now we "attach" the metadata to the data

attach(data_vegan.env)

#### ORDINATION: Detrended correspondence analysis ####

#ord <- decorana(data_vegan)
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
#plot(ord, disp="sites", type="n")
#ordihull(ord, habitat, col=1:2, lwd=3)
# ordiellipse(ord, habitat, col=1:2, kind = "ehull", lwd=3)
#ordiellipse(ord, habitat, col=1:2, draw="polygon")
#points(ord, disp="sites", pch=21, col=1:2, bg="yellow", cex=1.3)
# ordispider(ord, habitat, col=1:2, label = TRUE)
ord <- metaMDS(data_vegan %>%
                 filter(op_code != "CAG_017"),
               distance = "bray",
               k = 3,
               maxit = 999, 
               trymax = 500,
               wascores = TRUE) 
View(ord)

ggord <- 
  fortify(ord) %>% 
  tibble() %>% 
  clean_names() %>%
  filter(score == "sites") %>% 
  bind_cols(tibble(data_vegan.env)%>%
              filter(op_code != "CAG_017")) %>% 
  clean_names()

habitatcolors <- c("#6FAFC6", "#F08080")
habitat(habitatcolors) <- c("Shallow Reef", "Deep Reef")

#### nMDS of Fish Assemblage at Different Shallow and Deep Reefs at Cagayancillo and Tubbataha ####
studylocationcolors <- c("#C97CD5","#79CE7A")
study_locations(studylocationcolors) <- c("CAGAYANCILLO", "TRNP")
View(ord)
View(ggord)
ggord %>%
  # filter(label != 17) %>%
  ggplot(aes(x = nmds1,
             y= nmds2,
             color = studylocation_habitat,
             shape = site_code)) +
  scale_x_continuous(limits = c(-3,3)) +
  geom_point(size = 5) +
  stat_ellipse(aes(group = studylocation_habitat)) +
  theme_classic() +
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  labs(color = "Habitat", shape = "Study Locations", title = "NMDS Plots of Fish Assemblage") 
  # scale_color_manual(values = habitatcolors) 
save_plot("NMDSfishassemblageversion2groupings.png")

# ggord %>%
#   ggbiplot(ellipse = TRUE)

# ggord %>%
#   ggbiplot()
# 
# gg_ordiplot(ord, groups = data_vegan.env$habitat_mpa, pt.size = 3)

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

pool <- with(data_vegan.env, specpool(x = data_vegan, 
                                      pool = depth_cat,
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
  rename(sp_richness_est = value,
         estimator = name) %>%
  ggplot(aes(x= depth_cat,
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
p<-estaccumR(data_vegan, permutations = 50) %>%
  filter(estimator != "ace")
View(p)
p.plot<-plot(p, display = c("chao","ace"))
p.plot
# to make plot w ggplot, see https://stackoverflow.com/questions/52652195/convert-rarefaction-plots-from-vegan-to-ggplot2-in-r

#### vegan::poolaccum Extrapolated Species Richness Curve in a Species Pool Based on Presence Absence ####

# incidence based (presence / absence) richness rarefaction curve
# creates 1 curve per data frame, so if you want multiple curves, have to make them separately then combine into 1 tibble to plot
# increase permutations to 999 if you use this for your project
p<-poolaccum(data_vegan, permutations = 50)
p.plot<-plot(p, display = c("jack1","jack2", "chao", "boot"))
p.plot

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


