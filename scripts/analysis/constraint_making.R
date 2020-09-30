# making constraint file to use when building a phylogeny in RAXML

library(readr)
library(dplyr)

df <- read_csv('data/processed_data/meta_edgelist.csv')

summary_df <- df %>%
  select(diet_Species, diet_Family, diet_Order) %>%
  distinct() %>%
  filter(!is.na(diet_Species)) %>%
  arrange(diet_Species)


fasta_files <- list.files('data/co1_sequences/fastas/')
fasta_files <- gsub('.fasta', '', fasta_files)


outstr <- '('
for(i in 1:length(unique(summary_df$diet_Order))){
  ord <- unique(summary_df$diet_Order)[i]
  sp <- summary_df %>%
    filter(diet_Order == ord) %>%
    filter(diet_Species %in% fasta_files) %>%
    pull(diet_Species) %>%
    gsub(' ', '_', .)
  if(length(sp) > 1){
    sp_string <- paste0('(', paste(sp, collapse = ','), ')')
  }else{sp_string <- sp}
  
  
  outstr <- paste0(outstr, sp_string)
  if(i != length(unique(summary_df$diet_Order))){
    outstr <- paste0(outstr, ',')
  }else{
    outstr <- paste0(outstr ,');')
  }
}
outstr

write_file(outstr, 'data/co1_sequences/constraints.txt')


