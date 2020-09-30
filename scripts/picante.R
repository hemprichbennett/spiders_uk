library(picante)
library(readr)
library(dplyr)
library(tidyr)

prey_tree <- read.tree('data/co1_sequences/RAxML_bestTree.temp11')

phydist <- cophenetic(prey_tree)


edgelist <- read_csv('data/processed_data/edges_lab_and_field.csv')


# make a summary table of BINs per site
sample_table <- edgelist %>%
  select(SiteLampType, 
         #month_name, 
         consumer_BIN, consumer_id) %>%
  distinct() %>%
  select(-consumer_id) %>%
  group_by_all() %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = SiteLampType, values_from = n, 
              values_fill = 0)

sample_table

# the sample table shows that only two bins actually have enough 
# samples accross treatments to be viable for the analysis. This'll hopefully
# change when more samples are sequenced, and we're using binomials instead of 
# BINs. The code below also only retains them if they're present in ALL treatments,
# when perhaps being present in 2/3 should be sufficient

good_or_bad_df <- sample_table %>%
  mutate(good = ifelse(HPS > 0 & LED > 0 | 
                      unlit > 0 & LED >0 |
                        HPS > 0 & unlit >0,
                       T, F)) %>%
  select(consumer_BIN, good)

# make an adjacency matrix to use for calculating mpd
adj_mat <- edgelist %>%
  # remove any rows of BINs whicher were excluded above
  left_join(good_or_bad_df) %>%
  filter(good == T) %>%
  # replace spaces with underscores in species names
  mutate(diet_Species = gsub(' ', '_', diet_Species)) %>%
  # for the time being we'll remove all spiders
  filter(diet_Order != 'Araneae') %>%
  filter(!is.na(diet_Species)) %>%
  filter(diet_Species %in% rownames(phydist)) %>%
  select(consumer_id, diet_Species, SiteLampType) %>%
  group_by(consumer_id, diet_Species) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = diet_Species, values_from = n,
              values_fill = 0) 

consumer_ids <- adj_mat$consumer_id

full_mat <- adj_mat %>% 
  ungroup() %>% 
  select(-consumer_id) %>%
  as.matrix()

rownames(full_mat) <- consumer_ids



mpd_df <- ses.mpd(full_mat, phydist, null.model="taxa.labels",
        abundance.weighted=FALSE, runs=999) %>%
  tibble::rownames_to_column(var = 'consumer_id')

  

# add the field data etc to it

field_and_mpd <- edgelist %>%
  select(consumer_id, consumer_BIN, SiteName, Treatment, SiteLampType, 
         Date, year, month, month_name, day, VisitNum, location_type) %>%
  # get rid of non-unique rows
  distinct() %>%
  left_join(mpd_df) %>%
  # as a lot of the BINs aren't identified yet
  filter(!is.na(consumer_BIN)) %>%
  # as we've not calculated mpd on all samples, remove them
  filter(!is.na(mpd.obs)) 




temp_subset <- filter(field_and_mpd, consumer_BIN == 'BOLD:AAE4245')

trial_lm <- lm(mpd.obs ~ factor(consumer_BIN) + 
             factor(SiteLampType) 
             + factor(month_name), data = field_and_mpd)
summary(trial_lm)

MASS::stepAIC(trial_lm)

library(ggplot2)


ggplot(temp_subset, aes(x = SiteLampType, y = mpd.obs, colour = month_name))+
  geom_jitter()+
  facet_wrap(. ~ consumer_BIN) +
  theme_classic()+
  scale_colour_viridis_d()


# the same for faiths pd, out of curiosity


pd_df <- pd(full_mat, prey_tree) %>%
  tibble::rownames_to_column(var = 'consumer_id') 

pd_all <- edgelist %>%
  select(consumer_id, consumer_BIN, SiteName, Treatment, SiteLampType, 
         Date, year, month, month_name, day, VisitNum, location_type) %>%
  # get rid of non-unique rows
  distinct() %>%
  left_join(pd_df) %>%
  # as a lot of the BINs aren't identified yet
  filter(!is.na(consumer_BIN)) %>%
  # as we've not calculated mpd on all samples, remove them
  filter(!is.na(PD)) 

pd_subset <- filter(pd_all, consumer_BIN == 'BOLD:AAE4245')

pd_lm <- lm(PD ~ #factor(consumer_BIN) + 
                 factor(SiteLampType) 
               + factor(month_name), data = pd_all)
summary(pd_lm)
