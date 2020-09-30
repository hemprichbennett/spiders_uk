#### Header ####
## Project: 
## Script purpose: 
## Date: 
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

#### Setup ####

# Prevent partial-match errors 
options(warnPartialMatchArgs = TRUE)

# enable god-mode debugging
#options(error=recover)

library(tidyverse)
library(magrittr)
library(igraph)
library(bipartite)

input_df <- read_csv('data/processed_data/edges_lab_and_field.csv')


# make a list of tibbles, with each containing an edgelist of consumer BINs
# and their prey
nets <- input_df %>%
  # filter out any na values
  filter(!is.na(consumer_BIN) & !is.na(SiteLampType) & !is.na(diet_Species)) %>%
  # make a list for each of the SiteLampType variables
  group_split(SiteLampType) %>%
  # Set the list element names
  setNames(sort(unique(input_df$SiteLampType))) 
  

# now make a list of networks
copy_number_nets <- list()
for(i in 1:length(nets)){
  
  
  
  # make the tibble into a adjacency matrix
 network <- nets[[i]] %>%
   # make a summary column for the number of copies per consumer BIN and 
   # dietary species
  group_by(consumer_BIN, diet_Species) %>%
  summarise(n_copies = sum(n_copies), .groups = 'drop') %>%
   # now pivot the columns to make the edgelist into an adjacency matrix
   pivot_wider(names_from = diet_Species, values_from = n_copies,
               # fill the empty cells with 0
              values_fill = 0)
   
 # make a vector of all the consumer BINs, to use as rownames for the matrix
 bins_vec <- network$consumer_BIN
 
 network <- network %>%
   select(-consumer_BIN) %>%
   as.matrix()

 rownames(network) <- bins_vec
    
 copy_number_nets[[i]] <- network
 names(copy_number_nets)[i] <- names(nets)[i]
}


# calculate the real values for isa, and null-model values

isa_vals_list <- list() # list for the output to go into
it <- 1 # an iterator

for(i in 1:length(copy_number_nets)){
  
  netname <- names(copy_number_nets)[i] # the name of the network in question
  isa_vals_list[[it]] <- c( netname,
    networklevel(copy_number_nets[[i]], index = 'ISA'),
    'real')
  it <- it + 1
  
  
  for(j in 1:1000){
      null_val <- permatswap(copy_number_nets[[i]], fixedmar='both',mtype="count",times=1, method="quasiswap")$perm[[1]] %>%
        networklevel(index = 'ISA')
      
      isa_vals_list[[it]] <- c(netname, null_val, 'random')
      
      it <- it + 1
  }
  
  }

isa_tib <- bind_rows(isa_vals_list)

colnames(isa_tib) <- c('network', 'isa', 'origin')

isa_tib$isa <- as.numeric(isa_tib$isa)


isa_quartiles <- isa_tib %>% 
  filter(origin == 'random') %>%
  group_by(network) %>% 
  summarise(x = quantile(isa, c(0.05, 0.95)), q = c(0.05, 0.95)) %>%
  rename(value = x, quant = q) %>%
  pivot_wider(names_from = quant, values_from = value)

isa_values_final <- isa_tib %>%
  filter(origin == 'real') %>%
  select(-origin) %>%
  left_join(isa_quartiles)

isa_plot <- ggplot(isa_values_final, aes(x = network, y = isa)) +
  geom_point() +
  geom_errorbar(aes(ymin = `0.05`, ymax = `0.95`))+ 
  ylab('Interaction strength asymmetry') +
#  ylim(c(0,1))+
  theme_bw()

isa_plot  
ggsave('figures/isa.pdf', isa_plot)
