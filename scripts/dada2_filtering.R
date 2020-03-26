args <- as.character(commandArgs(trailingOnly = TRUE))

fwd_in <- args[1]

rev_in <- args[2]

fwd_out <- args[3]

rev_out <- args[4]

library(dada2)

out <- filterAndTrim(fwd_in, fwd_out, rev_in, rev_out,
  maxN = 0, maxEE = c(2, 2), rm.phix = TRUE,
  compress = TRUE, multithread = TRUE
) # On Windows set multithread=FALSE
