#### Header ####
## Project:
## Script purpose:
## Date:
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(tidyverse)
library(here)

# Read in data ------------------------------------------------------------



basic_df <- read.csv(
  here("data", "sequencing_planning", "basic_sequence_layout.csv"),
  header = T, stringsAsFactors = F
)

lab_data <- read.csv(
  here("data", "raw_data", "lab_data", "DNA Extractions Master.csv"),
  stringsAsFactors = F
)


# Add columns to the sequencing df ----------------------------------------

# First, remove all lab info other than sample ID and conc

lab_data <- lab_data %>%
  select(Sample.Unique.ID, BDNA.Conc.ng.uL) %>%
  # rename the column to match the one from the basic_df, so we can merge the
  # two dfs
  rename(Sample.ID = Sample.Unique.ID)

extra_data <- basic_df %>%
  # first, add the concentrations
  left_join(lab_data) %>%
  # then add a load of new columns
  mutate(
    # Isolate the primer name from the plate name
    primer = gsub(" .+", "", Plate.ID),
    # Isolate the run number from the plate name
    run = gsub(".+ |-.+", "", Plate.ID),
    # Pad the Sample.ID column so all numeric IDs are 3 digits
    Sample.ID = str_pad(Sample.ID, width = 3, side = "left", pad = "0"),
    # Get the row and column
    `Plate col` = gsub("[A-Z]", "", Well),
    `Plate row` = gsub("[0-9]", "", Well)
  ) %>%
  # Make a new column where primer and Sample.ID are combined
  unite(sample_UID, primer, Sample.ID,
    sep = "_",
    remove = F
  ) 
extra_data <- extra_data%>%
  # Make a column of if the sample is a replicate or not
  mutate(replicate = ifelse(duplicated(extra_data$sample_UID), "02", "01")) %>%
  # Add the replicate info to the Sample UID column
  unite(sample_UID, sample_UID, replicate,
    remove = F
  ) %>%
  # Finally, reorder the columns
  select(Plate.ID, Well, sample_UID, everything())


# Save the data
write_csv(
  extra_data,
  here("data", "sequencing_planning", "sequencing_metadata.csv")
)
