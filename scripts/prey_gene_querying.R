#### Header ####
## Project: 
## Script purpose: Finding which genes are available for the species matched
##                 by the RDP classifier
## Date: 2020-08-19
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

#### Setup ####

# Prevent partial-match errors 
options(warnPartialMatchArgs = TRUE)

# enable god-mode debugging
options(error=recover)

library(rentrez)
library(tidyverse)

prey_sp <- read_lines('results/diet_species.txt')

set_entrez_key("d688720c8f7d7a03bc9eeb2c3e33986e8508")

# dbs <- entrez_dbs()
# 
#   for(i in 1:length(dbs)){
#     print(entrez_db_summary(dbs[i]))}
# 

tax_rec <- entrez_search(
  db = "nuccore", 'Tibellus oblongus[ORGN]',
  rettype = "xml", parsed = TRUE, retmax = 500
)

summ <- entrez_summary(id = tax_rec$ids[3], db = 'nucleotide'#,
                     #rettype = 'fasta'
                 )


for(i in 1:length(summ)){
  print(i)
  print(names(summ)[i])
  print(summ[i])
}

outlist <- list()
it <- 1
for(i in 1:length(prey_sp)){
  # select the prey species of interest
  prey_species <- prey_sp[i]
  
  # query NCBI for what sequences exist for that species
  prey_seqs <- entrez_search(
    db = "nuccore", paste0(prey_species,'[ORGN]'),
    rettype = "xml", parsed = TRUE, retmax = 5000
  )
  
  if(length(prey_seqs$ids) >0){
    
    # query NCBI for summary info on each of the returned seqs
    
    for(s in 1:length(prey_seqs$ids)){
      
      # get the summary info
      summ <- entrez_summary(id = prey_seqs$ids[s], db = 'nucleotide')
  
      # store the sequence length
      sequence_length <- summ$slen
      # store the species name
      sp_name <- summ$organism
      gene_name <- 
      it <- it + 1
    }
    
    
  }else{
    cat('something went wrong with', prey_species, '\n')
  }
}


str_split(summ$extra, pattern = '\\|',simplify = T)[1,5] %>% gsub('.+\\.', '', .)
# for(i in 2:length(prey_sp)){
#   print(prey_sp[i])
#   
#   print(
#     entrez_search(
#       db = "taxonomy", term = paste0(prey_sp[i],'[ORGN]'),
#       retmax = 500
#     )
#   )
# }
