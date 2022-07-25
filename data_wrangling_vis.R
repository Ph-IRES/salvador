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
theme_set(
  theme_void()
)
# install.packages("sf")


#### USER DEFINED VARIABLES ####

inFilePath1 = "./WorkingData_CLEANED_TUB,CAG.xlsx"
inFilePath2 = "./PHIRES_MetaData.xlsx"

# outFilePath = "./data_combined.tsv"

#### READ IN DATA & CURATE ####

data <-
  read_excel(inFilePath1,
             na="NA") %>%
  clean_names() %>%
  # there is a different depth for CAG_024 in data vs metadata
  # solution: go with metadata depth, confirm with Gene & Rene
  mutate(depth_m = case_when(op_code == "CAG_024" ~ 8,
                             TRUE ~ depth_m))

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
         everything())

#### WRITE LONG FORMAT FILE ####

# data_all %>%
#   write_tsv(outFilePath)

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
study_locations(studylocationcolors) <- ("CAGAYANCILLO", "TRNP")
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
  
install.packages("purrr")
library(purrr)

subregions_keep %>%
  purrr::map(~ map_data("world",
                        region = "Philippines") %>%
               filter(subregion == .x))%>%
mutate(lat = case_when(lat < minLat ~ minLat,
                       lat > maxLat ~ maxLat, TRUE ~ lat),
       long = case_when(long < minLong ~ minLong,
                       long > maxLong ~ maxLong,
                       TRUE ~ long)) %>%
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
                  label = region),
              size = 10,
              hjust = 0.5,
              inherit.aes = FALSE) +
    geom_point(data = metadata,
               aes(x = long_e,
                   y = lat_n,
                   color = habitat),
               inherit.aes = FALSE) +
    theme_classic()
  
  metadata %>%
    ggplot(aes(y=lat_n,
               x=long_e,
               color = habitat)) +
    geom_point(size = 5) +
    theme_classic()
  
  

#### Mikaela's Data CleanUp and Modifications ####
View(data_all)
data_all <- data_all %>%
  mutate(study_locations = case_when(
    site == "Cawili" ~ "CAGAYANCILLO",
    site == "Calusa" ~ "CAGAYANCILLO",
    site == "Cagayancillo" ~ "CAGAYANCILLO",
    site == "TUBBATAHA" ~ "TRNP")) %>%
  mutate(family_clean = case_when(
    family == "Epinephelidae" ~ "Serranidae",
    TRUE ~ family))

View(data_all)

#### Mikaela's Data Visualization ####
habitatcolors <- c("#6FAFC6","#F08080")
habitat(habitatcolors) <- ("Shallow Reef", "Deep Reef")

#Barplot of MaxN per BRUV Station at TRNP and Cagayancillo
data_compiled <- data_all %>%
  group_by(study_locations, 
           habitat, 
           family_clean) %>%
  summarise(maxn_per_opcode = mean(max_n),
            sd = sd(max_n),
            n = n(),
            se_per_opcode = sd/sqrt(n)) %>%
  ggplot(aes(x = study_locations,
             y = maxn_per_opcode,
             fill = habitat))+
  geom_bar(position = "dodge", 
           stat = "identity") +
  xlab("Study Locations") +
  ylab("Mean MaxN per BRUV Deployment") +
  labs(title = "Mean MaxN at TRNP vs. Cagayancillo",
       fill = "Habitat") +
  theme_classic() +
  scale_fill_manual(values = habitatcolors) +
  geom_errorbar(aes(ymax = maxn_per_opcode + se_per_opcode,
                    ymin = maxn_per_opcode -
                      se_per_opcode), 
                position = "dodge")
data_compiled  
save_plot("MeanMaxNatTRNPvs.Cagayancillo.png")

#Faceted Barplot of MaxN at TRNP and Cagayancillo faceted by family 
data_compiled_faceted <- data_all %>%
  group_by(study_locations, 
           habitat, 
           family_clean) %>%
  summarise(maxn_per_opcode = mean(max_n),
            sd = sd(max_n),
            n = n(),
            se_per_opcode = sd/sqrt(n)) %>%
  ggplot(aes(x = study_locations,
             y = maxn_per_opcode,
             fill = habitat))+
  geom_bar(position = "dodge", 
           stat = "identity") +
  xlab("Study Locations") +
  ylab("Mean MaxN per BRUV Deployment") +
  labs(title = "Mean MaxN at TRNP vs. Cagayancillo",
       fill = "Habitat") +
  theme_classic() +
  scale_fill_manual(values = habitatcolors) +
  geom_errorbar(aes(ymax = maxn_per_opcode + se_per_opcode,
                    ymin = maxn_per_opcode -
                      se_per_opcode), 
                position = "dodge") +
  facet_grid(family_clean ~ .) +
  theme(strip.text.y.right = element_text(angle = 0))
data_compiled_faceted  
save_plot("FacetedMeanMaxNatTRNPvs.Cagayancillo.png")
