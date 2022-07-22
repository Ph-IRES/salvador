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
# install.packages("sf")
library(sf)

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


metadata %>%
  ggplot(aes(y=lat_n,
             x=long_e,
             color = habitat)) +
  geom_point(size = 5) +
  theme_classic()

world_map <- map_data("world")

ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="lightgray", colour = "white")

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
                   lat = mean(lat))

region_label_data <- 
  map_data("world",
           region = "Philippines") %>%
  dplyr::group_by(region) %>%
  dplyr::summarize(long = mean(long), 
                   lat = mean(lat))

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
                label = region),
            size = 10,
            hjust = 0.5,
            inherit.aes = FALSE) +
  geom_point(data = metadata,
             aes(x = long_e,
                 y = lat_n),
             inherit.aes = FALSE,
             color = "red")

