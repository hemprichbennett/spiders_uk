#### Header ####
## Project: CRAGinformatics
## Script purpose:
## Date:
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

start_time <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(rgbif)
library(purrr)

# Make a custom function function to use in the map call below
gbif_querying <- function(desired_country, desired_taxa) {
  cat(as.character(Sys.time()), desired_country, desired_taxa, "starting\n")
  output <- occ_search(
    scientificName = desired_taxa, limit = 10000,
    country = desired_country,
    hasCoordinate = T, return = "data"
  )
  # if there are no matches, occ_search returns a gbif object rather than a list
  # as gbif objects break the select call and subsequent dataframe making,
  # check the type of 'output'
   if (typeof(output) == "list") {
     # Some records don't have species names etc. These records break the select
     # call, are no use to us and should be discarded
     required_cols <- c("species", "name", "decimalLongitude",
     "decimalLatitude", "year",
     "individualCount", "country")
     # if statement translates to 'if all of the required_cols are in 
     # names(output)
     if(!FALSE %in% (is.element(required_cols, names(output)))){
       
       # Simplify occurrence dataframe, get rid of superfluous columns. This
       # is necessary as GBIF won't return empty columns, and this can mess 
       # with any attempts at combining dataframes later
       output <- dplyr::select(output,
                               species, name, decimalLongitude,
                               decimalLatitude, year,
                               individualCount, country
       ) %>%
         mutate(taxon = desired_taxa)
       
       cat(as.character(Sys.time()), desired_country, desired_taxa, "done\n")
       return(output) 
     }
     
  }
}

# Make a vector of shortnames for all African countries (GBIF doesn't accept
# the full name of a country, but instead their 2-letter ISO-3166-1 code, see
# http://en.wikipedia.org/wiki/ISO_3166-1_alpha-2 for details)

african_countries <- c("AO", "BF", "BI", "BJ", "BW", "CD", "CF", "CG", "CI", 
                       "CM", "CV", "DJ", "DZ", "EG", "EH", "ER", "ET", "GA", 
                       "GH", "GM", "GN", "GQ", "GW", "KE", "KM", "LR", "LS", 
                       "LY", "MA", "MG", "ML", "MR", "MU", "MW", "MZ", "NA", 
                       "NE", "NG", "RE", "RW", "SC", "SD", "SH", "SL", "SN", 
                       "SO", "SS", "ST", "SZ", "TD", "TG", "TN", "TZ", "UG", 
                       "YT", "ZA", "ZM", "ZW")

# Make a vector of desired taxonomic groups

taxa_vec <- c("Chiroptera", "Aves", "Squamata", "Anura")

# Make a dataframe with all possible combinations of the two desired variables.
# This makes it easier to use purrr later
combinations <- crossing(african_countries, taxa_vec)

# Use our custom function on the above dataframe to query for every combination
# of country and taxonomic group
occurrence_df <- purrr::map2_df(
  combinations$african_countries,
  combinations$taxa_vec,
  gbif_querying
)

# Save the output
write.csv(occurrence_df,
          "data/processed_data/taxa_occurences.csv",
          row.names = F)

end_time <- Sys.time()


# temp <- gbif_querying(desired_country = 'CV', desired_taxa = 'Anura')
# str(temp)
