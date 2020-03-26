##### setup ####

library(ggplot2)
library(reshape2)
library(ggbeeswarm)
library(dplyr)
library(rgbif)

args <- as.character(commandArgs(trailingOnly = TRUE))

if(length(args) >0){
  infile <- args[1]
  basedir <- args[2]
}else{
  infile <- 'temp_output_files/classified_sequences.tsv'
  basedir <- getwd()
}

# Where figures will be plotted
plotdir <- paste0(basedir, "/figures/")
if(!dir.exists(plotdir)){dir.create(plotdir)}


assignment_df <- read.table(infile, sep = "\t", row.names = NULL,
                            stringsAsFactors = F)

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


only_scores <- assignment_df[, -grep("_name", colnames(assignment_df))]


##### Plot scores for each level #####


melted_scores <- melt(only_scores)
colnames(melted_scores) <- c("Taxa", "Confidence")
melted_scores$Taxa <- gsub("_score", "", melted_scores$Taxa)
melted_scores$Taxa <- factor(melted_scores$Taxa,
  levels = unique(melted_scores$Taxa)
)

bee_plot <- ggplot(melted_scores, aes(x = Taxa, y = Confidence)) + geom_quasirandom() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  ) + # Kill annoying formatting
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # rotate axis labels
  geom_hline(
    yintercept = 0.95, linetype = "dashed",
    color = "red", size = 0.5
  )
pdf(paste0(basedir, "/figures/assignment_scores.pdf"))
bee_plot
dev.off()

##### show which sequences met the required scores #####

score_retrieval <- function(in_vector) {
  if (length(which(in_vector >= 0.95)) == 0) {
    outval <- NA
  } else {
    matches <- which(in_vector >= 0.95)
    outval <- colnames(only_scores)[
      max(matches)
    ]
  }
  return(outval)
}
bestmatches <- apply(only_scores, 1, score_retrieval)
best_df <- as.data.frame(table(bestmatches))
colnames(best_df) <- c("Taxa", "n_matches")
best_df$Taxa <- gsub("_score", "", best_df$Taxa)

best_df$Taxa <- factor(best_df$Taxa,
  levels = c(
    "Domain", "Superkingdom", "Kingdom", "Phylum",
    "Class", "Order", "Family", "Genus", "Species"
  )
)

best_plot <- ggplot(best_df, aes(x = Taxa, y = n_matches)) + geom_bar(stat = "identity") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  ) + # Kill annoying formatting
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # rotate axis labels

pdf(paste0(basedir, "/figures/taxa_matches.pdf"))
best_plot
dev.off()


##### Taxonomy querying ####

validated_df <- assignment_df
score_columns <- grep('score', colnames(validated_df))

for(i in 1:length(score_columns)){ #Kill any taxonomic names which have a confidence lower than 0.6
  score_col <- score_columns[i]
  to_kill <- which(validated_df[,score_col] < 0.6)
  validated_df[to_kill, score_col-1] <- "NO MATCH"
}

validated_df <- validated_df[,-grep('score', colnames(validated_df))] #Get rid of the score columns
validated_df$seqname <- rownames(validated_df)

melted_validated <- melt(validated_df, id.vars = 'seqname')
colnames(melted_validated)[c(2,3)] <- c('Level', 'Name')
melted_validated$Level <- gsub('_name', '', melted_validated$Level)

melted_validated <- melted_validated %>%
  group_by(Level, Name) %>%
  summarise(count = n())

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

ggplot(filter(melted_validated, Level == 'Order'), aes(x = '', y = Name, fill = Name))+
  geom_bar(width = 1, stat = 'identity')+
  coord_polar('y', start = 0)+  # make it a pie chart
  scale_fill_viridis_d()+
  blank_theme




#Gbif stuff (unfinished)
goodsp <- filter(assignment_df, Species_score >= 0.6)

goodsp$Species_name <- gsub('_', ' ', goodsp$Species_name)

z <- name_backbone(name = goodsp$Species_name[1])
