#### Header ####
## Project: spiders_uk
## Script purpose: analysing the taxonomic composition of the ASVs
## Date: 2020-03-24
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

#### Setup ####

# Prevent partial-match errors
options(warnPartialMatchArgs = TRUE)

# enable god-mode debugging
# options(error = recover)

library(tidyverse)
library(gridExtra)


# Format data ---------------------------------------------------------

# Edgelist

meta_edgelist <- read_csv("data/processed_data/meta_edgelist.csv")

# the filenames were slightly different from expected, so the consumer_id
# field is wrong. Thankfully the information is still in the file field, so
# we can get it from there. Also, interaction_type and run are redundant for
# this analysis, so lets get rid of them

meta_edgelist <- meta_edgelist %>%
  # get rid of redundant columns
  select(-run, -interaction_type, -consumer_id) %>%
  # split the file column into three columns
  separate(file, into = c("primer", "consumer_id", "replicate"), sep = "-") %>%
  # get rid of superfluous text from 'replicate' field
  mutate(replicate = gsub("_.+", "", replicate))



# Field data

field_data <- read_csv("data/to_sequence.csv")

# left-pad the ID column, to match the corresponding column in meta_edgelist

field_data <- field_data %>%
  mutate(Sample.Unique.ID = str_pad(
    string = Sample.Unique.ID, width = 3,
    side = "left", pad = "0"
  ))


# Now join both tibbles together, to make one object to rule them all

all_data <- left_join(meta_edgelist, field_data, 
                      by = c("consumer_id" = "Sample.Unique.ID")) %>%
  # Make the taxonomic Order a factor, then reverse it as ggplot by default
  # plots it the wrong way round ('A' at bottom, 'Z' at top)
  mutate(diet_Order = as.factor(diet_Order),
         diet_Order = fct_rev(diet_Order))


# quick summary of the spider/non-spider reads and ASVs

all_data %>%
  mutate(araneae = ifelse(diet_Order == 'Araneae', T, F)) %>%
  group_by(primer, araneae) %>%
  summarise(n_reads = sum(n_copies),
            n_asvs = length(unique(diet_asv)))


# Now we need to do 4 things:

# Compare the taxonomic composition of the overall dataset for each primer -----

no_spiders <- all_data %>%
  # Get rid of all spider ASVs for now, as we assume they're consumer
  filter(diet_Order != "Araneae") %>%
  # get rid of superfluous columns
  select(
    primer, consumer_id, replicate, diet_asv, n_copies, SiteLampType,
    SiteName, diet_Order
  )

for_taxa_plot <- no_spiders %>%
  group_by(primer, diet_Order) %>%
  # make summary columns to give the values for plotting
  summarise(
    # The number of ASVs, recounting for each sample
    n_ASVs = n(), 
    # The number of reads overall
    n_copies = sum(n_copies), 
    # The number of unique ASVs overall
    n_unique_ASVs = length(unique(diet_asv))) 


# Make a basic theme for use in the three tileplots below
tileplot_theme <- theme_classic() +
  theme(legend.position = "bottom",
  text = element_text(size=15))

# Make the three plots
copies_plot <- ggplot(for_taxa_plot, aes(x = primer, y = diet_Order)) +
  geom_tile(aes(fill = n_copies)) +
  scale_fill_gradient(
    name = "Number of reads overall assigned to Order",
    breaks = seq(0,max(for_taxa_plot$n_copies), 4000)
  ) +
  tileplot_theme +
  labs(x = "Primer", y = "Taxonomic Order")

copies_plot

overall_asvs_plot <- ggplot(for_taxa_plot, aes(x = primer, y = diet_Order)) +
  geom_tile(aes(fill = n_ASVs)) +
  # Make the fill go between white and blue, depending on number of ASVs
  scale_fill_gradient(
    name = "Number of ASVs assigned overall\n
    (e.g. ASVs counted twice if found in two samples)"
  ) + tileplot_theme+
  labs(x = "Primer", y = "Taxonomic Order")

overall_asvs_plot

unique_asvs_plot <- ggplot(for_taxa_plot, aes(x = primer, y = diet_Order)) +
  geom_tile(aes(fill = n_unique_ASVs)) +
  # Make the fill go between white and blue, depending on number of ASVs
  scale_fill_gradient(
    name = "Number of unique ASVs detected\n
    (i.e. ASV only counted once, even if found in 100 samples)"
  ) + tileplot_theme+
  labs(x = "Primer", y = "Taxonomic Order")

unique_asvs_plot

# put plots into a grid
taxa_grid <- grid.arrange(copies_plot, overall_asvs_plot, unique_asvs_plot, ncol = 3)

# save grid
ggsave("figures/primer_biases/overall_primer_biases.pdf", taxa_grid, width = 24)
ggsave("figures/primer_biases/overall_primer_biases.jpg", taxa_grid, width = 24)

# Compare the results within-sample for the two primer pairs
# Compare the replicats
# Look at differences between ecological treatments
