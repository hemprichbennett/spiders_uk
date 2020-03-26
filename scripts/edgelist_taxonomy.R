#### Header ####
## Project:
## Script purpose: Adding taxonomy to an edgelist
## Date: 2019-04-25
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(tidyr)
library(dplyr)

args <- as.character(commandArgs(trailingOnly = TRUE))

if (length(args) > 0) {
  in_matrix_file <- args[1]
  in_assignments <- args[2]
  basedir <- args[3]
  args_provided <- T
} else {
  in_matrix_file <- "temp_output_files/out_edgelist.csv"
  in_assignments <- "temp_output_files/classified_sequences.tsv"
  basedir <- getwd()
  args_provided <- F
}

##### Read in primer info #####

primers <- read.csv(paste0(basedir,'/data/primers/primertypes.csv'), stringsAsFactors = F, 
                    row.names = NULL)

##### Read in and format interaction data #####
if(args_provided == T){
  edgelist <- read.csv(in_matrix_file,
                       stringsAsFactors = F,
                       row.names = NULL
  )
}else{
  edgedir <- 'data/processed_data/interactions/'
  edgefiles <- list.files(edgedir, pattern = 'simple_edgelist', recursive = T, include.dirs = T)
  edgefiles <- paste0(edgedir, edgefiles)
  edgefiles <- lapply(edgefiles, function(x) read.csv(x, stringsAsFactors = F, row.names = NULL))
  edgelist <- do.call(rbind, edgefiles)
}

# Get primer names
edgelist$primer <- gsub('-.+$', '', edgelist$file)
# get filenames
edgelist$consumer_id <- gsub('^.+-', '', edgelist$file)
edgelist$consumer_id <- gsub('_L001', '', edgelist$consumer_id)

colnames(edgelist) <- gsub('^asv$', 'diet_asv', colnames(edgelist))

# Add interaction types

edgelist <- edgelist %>%
  left_join(primers, by = 'primer')

#### Setup the assignment data from RDP classifier ####
assignment_df <- read.table(in_assignments,
  sep = "\t",
  row.names = NULL,
  stringsAsFactors = F
)

# Get rid of redundant columns
assignment_df <- assignment_df[, -c(2, 3, 7, 10, 13, 16, 19, 22, 25, 28)]

colnames(assignment_df) <- c(
  "diet_asv", "diet_Domain_name", "diet_Domain_score",
  "diet_Superkingdom_name", "diet_Superkingdom_score",
  "diet_Kingdom_name", "diet_Kingdom_score",
  "diet_Phylum_name", "diet_Phylum_score",
  "diet_Class_name", "diet_Class_score",
  "diet_Order_name", "diet_Order_score",
  "diet_Family_name", "diet_Family_score",
  "diet_Genus_name", "diet_Genus_score",
  "diet_Species_name", "diet_Species_score"
)


# get rif of any asvs which aren't arthropoda
assignment_df <- filter(assignment_df, diet_Phylum_name == "Arthropoda")

# get rid of any asvs which couldnt confidently be assigned to at least family
assignment_df <- filter(assignment_df, diet_Family_score >= 0.6)

##### Find best assignments ####

# Make simple vectors to be used in the function 'best_assignment'
scorecols <- grep("score", colnames(assignment_df))
namecols <- grep("name", colnames(assignment_df))

taxa_vec <- colnames(assignment_df)[scorecols] # The taxonomic levels, in order
taxa_vec <- gsub("_score", "", taxa_vec)



# Function giving the best taxonomic assignment possible for each asv which
# met the above criteria
best_assignment <- function(desired_row) {
  bestpos <- max(
    which(
      desired_row[, scorecols] >= 0.6
    )
  ) # This is the position in both score and names of the best match

  # The name of the taxonomic level of the best match
  best_diet_taxa <- taxa_vec[bestpos]
  best_diet_taxa <- gsub('diet_', '', best_diet_taxa) # Remove the diet_ bit from the column name

  best_confidence <- desired_row[, scorecols[bestpos]] # The score
  best_diet_ID <- desired_row[, namecols[bestpos]] # The actual taxonomic label
  best_diet_ID <- gsub("_", " ", best_diet_ID) # remove any underscores (usually in species names)

  outdf <- data.frame(diet_asv = desired_row$diet_asv, best_diet_taxa, best_confidence, best_diet_ID)
  return(outdf)
}

# Use the function to find the best match for each asv
bestmatches <- do.call(rbind, lapply(
  seq(1, nrow(assignment_df)),
  function(x) best_assignment(assignment_df[x, ])
))

##### Make a dataframe where there are only NAs for taxa with a score < 0.6 ####

filtered_df <- assignment_df

for (i in 1:length(scorecols)) {
  badrows <- which(assignment_df[, scorecols[i]] < 0.6) # the rows which had
  # too-low a confidence level
  if (length(badrows) > 0) { # If there is a single value below 0.6 in that column

    # Make NA the taxa which can't confidently be assigned
    filtered_df[badrows, scorecols[i] - 1] <- NA
  }
}

# Now get rid of the scorecols and rename the columns,
# for better format post-join
filtered_df <- filtered_df[, -scorecols]
colnames(filtered_df) <- gsub("_name", "", colnames(filtered_df))

##### Join datasets ####
uber_edgelist <- edgelist %>%
  inner_join(bestmatches, by = "diet_asv") %>%
  inner_join(filtered_df, by = "diet_asv")
# This removes any asvs from the edgelist which aren't matched in the
# assignment_df. Any missing values will be because we removed them in the
# above QC

uber_edgelist$diet_Species <- gsub('_', ' ', uber_edgelist$diet_Species)

out_df_name <- paste0(basedir, "/data/processed_data/meta_edgelist.csv")

write.csv(uber_edgelist, out_df_name, row.names = F)
