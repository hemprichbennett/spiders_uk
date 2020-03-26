#### Header ####
## Project:
## Script purpose: script isolates the good species scores, for use in a later
## python script to assign them to long sequences from BOLD
## Date: 2019-04-26
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################


library(dplyr)
args <- as.character(commandArgs(trailingOnly = TRUE))

if (length(args) > 0) {
  infile <- args[1]
  outfile <- args[2]
} else {
  infile <- "data/processed_data/classifications.tsv"
  outfile <- "temp_output_files/species_names.txt"
}


# Read in and format the data
assignment_df <- read.table(infile,
  sep = "\t", row.names = NULL,
  stringsAsFactors = F
)

colnames(assignment_df) <- c(
  NA, NA, NA, "Domain_name", "Domain_score",
  "Superkingdom_name", NA, "Superkingdom_score",
  "Kingdom_name", NA, "Kingdom_score",
  "Phylum_name", NA, "Phylum_score",
  "Class_name", NA, "Class_score",
  "Order_name", NA, "Order_score",
  "Family_name", NA, "Family_score",
  "Genus_name", NA, "Genus_score",
  "Species_name", NA, "Species_score"
)

# kill the unneccesary colnames
assignment_df <- assignment_df[, -which(
  is.na(
    colnames(assignment_df)
  )
)]

# Filter by confidence score, retain only species name
goodsp <- assignment_df %>%
  filter(Species_score >= 0.6) %>%
  select(Species_name)

goodsp <- gsub("_", " ", goodsp$Species_name) # Replace underscores in spnames

goodsp <- unique(goodsp)

# save file
write(goodsp, outfile)
