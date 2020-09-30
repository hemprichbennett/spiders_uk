#### Header ####
## Project: CRAGinformatics
## Script purpose: using DADA2 to learn errors, merge sequences to ASVS, then
## spit out useful outputs
## Date: 2019-05-01
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################



##### Set up: read in files, sort their names ####

library(dada2)

args <- as.character(commandArgs(trailingOnly = TRUE))
if (length(args) > 0) {
  filtpathF <- args[1] # the directory containing your filtered forward fastqs
  filtpathR <- args[2] # the directory containing your filtered reverse fastqs
  run_path <- args[3]
  basepath <- args[3]
} else {
  filtpathF <- "temp_output_files/f_dir/"
  filtpathR <- "temp_output_files/r_dir/"
  run_path <- "data/raw_data/interactions/run1/"
  basepath <- getwd()
}

cat("filtpathF is", filtpathF, "\n")
cat("filtpathR is", filtpathR, "\n")
cat("run_path is", run_path, "\n")
cat("basepath is", basepath, "\n")

# Isolate the run name
run_name <- strsplit(run_path, split = "\\/")[[1]] # make a vector of split path
run_name <- run_name[length(run_name)] # select its final element


# Where all the files will all be saved
outdir <- paste0(basepath, "/data/processed_data/interactions/", run_name)
cat("outdir is", outdir)
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

# where the fastas will be saved
individual_fasta_dir <- paste0(outdir, "/sample_asv_fastas/")
if (!dir.exists(individual_fasta_dir)) {
  dir.create(individual_fasta_dir)
}



# File parsing

filtFs <- list.files(filtpathF, full.names = TRUE)
filtRs <- list.files(filtpathR, full.names = TRUE)

sample.names <- gsub(".+\\/", "", filtFs)
sample.names <- gsub("_f$|_r$", "", sample.names)
names(filtFs) <- sample.names
names(filtRs) <- sample.names


set.seed(100)

##### Learn error rates for sequence files #####

# Learn forward error rates
errF <- learnErrors(filtFs, nbases = 1e8, multithread = TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases = 1e8, multithread = TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names


##### Merge the files ####

for (sam in sample.names) {
  cat("Processing:", sam, "\n")
  # r_name <- gsub('f', 'r', sam)
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err = errF, multithread = TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err = errR, multithread = TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  #  print(merger)
  mergers[[sam]] <- merger
}
rm(derepF)
rm(derepR)


# Save a fasta for each sequencing file
for (nm in names(mergers)) {
  mrg <- mergers[[nm]]
  if (nrow(mrg) > 0) {
    uniquesToFasta(mrg,
                   paste0(outdir, "/sample_asv_fastas/", nm, ".fasta"))
  }
}



#### Construct sequence table and remove chimeras #####
seqtab <- makeSequenceTable(mergers)
seqtab <- removeBimeraDenovo(seqtab,
                             method = "consensus",
                             multithread = TRUE,
                             verbose = TRUE)


# save a fasta for the whole run

uniquesToFasta(seqtab,
               paste0(outdir, "/all_asvs.fasta"),
               ids = colnames(seqtab))


# Convert seqtab to an edgelist
seqtab <- as.data.frame(seqtab)
seqtab <- cbind(rownames(seqtab), seqtab)
colnames(seqtab)[1] <- "file"

edgelist <- reshape2::melt(seqtab,
  id.vars = "file", variable.name = "asv",
  value.name = "n_copies"
)

# add run name
edgelist$run <- run_name

# Get rid of null-edges to save space
edgelist <- dplyr::filter(edgelist, n_copies != 0)

write.csv(edgelist,
  paste0(outdir, "/simple_edgelist.csv"),
  row.names = F
)

write.csv(seqtab, 
          paste0(outdir, "/seqtab.csv"),
          row.names = F)