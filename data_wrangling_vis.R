#### INITIALIZATION ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(janitor)
library(readxl)

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
  rename(bait_weight_grams = weight_grams)

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
  ggplot(aes(x=lat_n,
             y=long_e,
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

# 
# data %>%
#   group_by(pipettor,
#            channel,
#            trial) %>%
#   summarize(mean_mass_g = mean(mass_g),
#             sd_mass_g = sd(mass_g)) %>%
#   ggplot(aes(x=channel,
#              y=mean_mass_g,
#              color = pipettor)) +
#   geom_point() +
#   geom_errorbar(aes(ymin=mean_mass_g - sd_mass_g,
#                     ymax = mean_mass_g + sd_mass_g)) +
#   geom_hline(yintercept = 0.013,
#              color = "grey",
#              linetype = "dashed") +
#   theme_classic() +
#   facet_grid(trial ~ pipettor,
#              scales = "free_x")
# ggsave("mean-mass_vs_channel_x_pipettor.png")
# 
# data %>%
#   ggplot(aes(x=order,
#              y=mass_g,
#              color = pipettor)) +
#   geom_point() +
#   geom_smooth() +
#   geom_hline(yintercept = 0.013,
#              color = "grey",
#              linetype = "dashed") +
#   theme_classic() +
#   facet_grid(. ~ trial,
#              scales = "free_x")
# ggsave("mass_vs_order_x_pipettor.png")
# 
# data %>%
#   group_by(pipettor,
#            channel,
#            trial) %>%
#   summarize(mean_mass_g = mean(mass_g),
#             sd_mass_g = sd(mass_g),
#             order = min(order)) %>%
#   ggplot(aes(x=order,
#              y=sd_mass_g,
#              color = pipettor)) +
#   geom_point() +
#   geom_smooth() +
#   geom_hline(yintercept = 0,
#              color = "grey",
#              linetype = "dashed") +
#   theme_classic() +
#   facet_grid(. ~ trial,
#              scales = "free_x")
# ggsave("sd-mass_vs_order_x_pipettor.png")
