
# Use the working directory to identify if we're working on ARCUS or locally
working_dir <- getwd()
basedir <- strsplit(working_dir, split = "/")[[1]][2]



if (basedir == "Users") {
  local_job <- T

  library(rentrez)
  library(lubridate)
  library(purrr)
} else {
  local_job <- F

  library(rentrez, lib.loc = "/home/zool2291/r_packages")
  library(lubridate, lib.loc = "/home/zool2291/r_packages")
  library(purrr, lib.loc = "/home/zool2291/r_packages")
}




# This is Dave messing around, trying to identify how best to get all available
# arthropod CO1 sequences from genbank, and their attached phylogenies. In theory
# once this has been done, I can then use RDPtools to turn this into a reference
# database (a la Porter and Hajibabaeim 'Automated high throughput animal CO1
# metabarcode classification') to query with our fasta file of ASVs

# This code was written using the guide available at
# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html


##### Function defining

# Now what we want to do is make a list where we have the species name,
# the ids of all of the sequences, and the taxonomy of the species

ncbi_querying <- function(sp) {
  bad <- F # to be changed below if this is a bad query
  sp_taxonomic_id <- sp$uid
  sp_name <- sp$scientificname
  cat("starting looking at", sp_name, "
")


  ## First, lets get the taxonomy

  tax_rec <- entrez_fetch(
    db = "taxonomy", sp_taxonomic_id,
    rettype = "xml", parsed = TRUE
  )
  print("obtained tax_rec")
  print(str(tax_rec))
  taxlist <- XML::xmlToList(tax_rec)$Taxon
  print("obtained taxlist")
  print(str(taxlist))
  lineagelist <- taxlist$LineageEx
  print("obtained lineaglist")
  print(str(lineagelist))
  print("obtained taxonomy")

  # Lots of the returned lineage items per species are useless as they have no
  # rank. Lets find out which and dispose of them
  ranknames <- purrr::map_chr(lineagelist, "Rank")
  badtaxons <- grep("no rank", ranknames)
  lineagelist[badtaxons] <- NULL
  if (length(lineagelist) == 0) {
    bad <- T
    # For some reason lineagelist was empty, so we don't want to retain this
  } else {
    ## Then lets get the sequences

    seq_search_term <- paste0(
      "cox1[gene] OR coxI[gene] OR CO1[gene] OR COI[gene] AND ",
      sp_name,
      "[organism]"
    )
    cat("starting looking for sequences for", sp_name, "
")
    seq_matches <- entrez_search(
      db = "nucleotide", term = seq_search_term,
      use_history = T
    )
    cat("got sequence ids for", sp_name, "
")
    # We need to store the web history to reduce the number of objects passing
    # between our session and NCBI (see
    # https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#web_history)
    #
    nseqs <- seq_matches$count # The number of sequences
    cat("this iteration has", nseqs, "sequences", "
")

    if (nseqs == 0) {
      bad <- T
      # Here there were no sequences matching the species name,
      # so we don't want it
    } else {
      print("starting sequence download")
      # Now get all of the sequences that are stored in our web history
      seqs <- entrez_fetch(
        db = "nucleotide", web_history = seq_matches$web_history,
        rettype = "fasta"
      )

      print("sequence download successful")
      seqs_split <- strsplit(seqs, split = ">")
      # The first item of the first vector in the list is now empty, due to the
      # splitting process. Delete it
      seqs_split <- seqs_split[[1]][-1]



      # Use the format_seqs function, as defined elsewhere, to format the list

      all_cat <- purrr::map(seqs_split, format_seqs)




      # The guide online states 'If you really wanted to you could also use
      # web_history objects to download all those thousands of COI sequences.
      # When downloading large sets of data, it is a good idea to take advantage of
      # the arguments retmax and restart to split the request up into smaller chunks.
      # For instance, we could get the first 200 sequences in 50-sequence chunks'.
      # Its probably worth trying to incorporate this, but I haven't yet...

      # Save each seq file here, they may be useful later when we want to reuse the
      # reference sequences for any sort of phylogenetic stuff
      outfastapath <- paste0(
        seqs_dir,
        gsub(" ", "_", sp_name),
        ".fasta"
      )
      print(outfastapath)
      # sink(outfastapath)
      cat(unlist(all_cat), file = outfastapath)
      # sink()
      cat("saving", sp_name, "successful\n")
      ncbi_objects <- list(lineagelist, all_cat)

      return(ncbi_objects)
    }
  }
}

# Function to tidy up each sequence list item, used in the ncbi_querying function
format_seqs <- function(desired_sequence) {
  # split it by newlines
  desired_split <- strsplit(desired_sequence, split = "\n")[[1]]
  seqname <- desired_split[1] # Get the sequence name

  # Because of using it in the strsplit call above, only the first sequence name
  # will start with a '>'. Replace the missing ones

  if (grepl(">", seqname) == F) {
    seqname <- paste0(">", seqname)
  }

  # Get the sequence itself
  sequence <- paste(desired_split[2:length(desired_split)], collapse = "")
  out <- c(seqname, "\n", sequence, "\n")

  return(out)
}




##### Setup

# For the non-interactive function, its important to store when the project
# started and the desired taxa to be worked on


desired_taxa <- "arthropoda" # in the final form of the script this should be specified when starting
start_time <- now()
start_string <- gsub(" ", "_", as.character(start_time))
start_string <- gsub(":", "-", start_string)
out_base_dir <- paste0("data/my_reference_files/", desired_taxa, "/", start_string)
seqs_dir <- paste0(out_base_dir, "/sequences/")
taxa_dir <- paste0(out_base_dir, "/taxonomy/")
dir.create(seqs_dir, recursive = T)
dir.create(taxa_dir, recursive = T)


# Set API key. This is Dave's, not to be used by anyone else as its paired to
# my email address!
set_entrez_key("d688720c8f7d7a03bc9eeb2c3e33986e8508")

# Lets quickly query the taxonomy database to find out how many matches there are
# for our desired taxa.
# This still needs to be automated by using the desired_taxa variable above!

full_search <- entrez_search(
  db = "taxonomy",
  term = '"Arthropoda"[ORGN] AND "species"[RANK]',
  retmax = 0
)
n_taxonomies <- full_search$count


cat("There are", n_taxonomies, "matches to our search query")


if (local_job == T) {
  # an arbitrary number to use if on a local iteration
  big_n <- 50000
} else {
  # if on arcus run the full thing
  big_n <- n_taxonomies
}




# the use of a loop, with seq_start provides a limiter, letting us only query
# NCBI for that many species records at a time. Whilst working in batches in a loop
# is computationally inefficient, if we ask for too much at once NCBI will kill
# the connection, therefore crashing everything.

search_max <- 50 # The maximum to use at a time

big_query_list <- list() # Where everything will ultimately be saved
for (seq_start in seq(1, big_n, search_max)) {
  cat("seq_start is", seq_start, "
")

  print("starting r_search")
  r_search <- entrez_search(
    db = "taxonomy",
    term = '"Arthropoda"[ORGN] AND "species"[RANK] ',
    retmax = search_max, retstart = seq_start, use_history = T
  )
  print("r_search complete")

  # Ideally in future we'll get all the arthropod CO1 genes, plus some outgroup
  # sequences


  print("starting multi_summs")
  # Get multiple summary items (just the first hundred for now though)
  multi_summs <- entrez_summary(db = "taxonomy", id = r_search$ids)

  sp_names <- extract_from_esummary(multi_summs, "scientificname")
  print("multi_summs finished")


  # Filter the returned species we now have taxonomy for by species name,
  # as we only want species with proper names

  # head(sp_names)
  # We only want true, named species
  bad_locs <- grep("sp\\.|nr\\.|aff\\.|cf\\.|\\/", sp_names)
  # head(bad_locs)
  multi_summs <- multi_summs[-bad_locs]
  cat("(length(multi_summs)) is", length(multi_summs), "
")


  ##### Use the functions on the queried subset
  if (length(multi_summs) == 0) {
    next()
  }

  for (m in 1:length(multi_summs)) {
    id <- multi_summs[[m]]$uid
    big_query_list[[id]] <- ncbi_querying(multi_summs[[m]])
  }

  # big_query_list[[seq_start]] <- map(multi_summs, ncbi_querying)
  # The names of this are the uids of each item
}
