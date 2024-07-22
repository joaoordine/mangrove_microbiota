# Genome Mining - RAST web server - processing resulting outputs in R

library(dplyr)
library(tidyverse)

## Load the two TSV files into R (the category table and the feature table for each isolate)
features_SA1 <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/features_SA1.tsv")  
features_SA6 <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/features_SA6.tsv")
features_SB1 <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/features_SB1.tsv")
features_SB8 <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/features_SB8.tsv")
features_WA5 <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/features_WA5.tsv")
features_WB3 <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/features_WB3.tsv")
features_WB8 <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/features_WB8.tsv")
SA1_df <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/SA1_df.tsv")
SA6_df <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/SA6_df.tsv")
SB1_df <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/SB1_df.tsv")
SB8_df <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/SB8_df.tsv")
WA5_df <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/WA5_df.tsv")
WB3_df <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/WB3_df.tsv")
WB8_df <- read_tsv("/home/strawberry/Documents/Collaborations/mangrove/RAST-genome_mining/WB8_df.tsv")

## Remove NA columns 
clean_na_columns <- function(df) {
  df %>%
    select_if(~ !all(is.na(.)))
} # function to remove NA columns

features_SA1 <- clean_na_columns(features_SA1)
features_SA6 <- clean_na_columns(features_SA6)
features_SB1 <- clean_na_columns(features_SB1)
features_SB8 <- clean_na_columns(features_SB8)
features_WA5 <- clean_na_columns(features_WA5)
features_WB3 <- clean_na_columns(features_WB3)
features_WB8 <- clean_na_columns(features_WB8)
SA1_df <- clean_na_columns(SA1_df)
SA6_df <- clean_na_columns(SA6_df)
SB1_df <- clean_na_columns(SB1_df)
SB8_df <- clean_na_columns(SB8_df)
WA5_df <- clean_na_columns(WA5_df)
WB3_df <- clean_na_columns(WB3_df)
WB8_df <- clean_na_columns(WB8_df) # apply the function to each dataframe 

## Make 'Features' column in long format, separating each feature identified using commas
make_long_format <- function(df) {
  df %>%
    separate_rows(Features, sep = ",")
} # Function to make 'Features' column in long format

SA1_df <- make_long_format(SA1_df)
SA6_df <- make_long_format(SA6_df)
SB1_df <- make_long_format(SB1_df)
SB8_df <- make_long_format(SB8_df)
WA5_df <- make_long_format(WA5_df)
WB3_df <- make_long_format(WB3_df)
WB8_df <- make_long_format(WB8_df)

## Rename features column so it has the same name on both dataframes 
rename_features <- function(df) {
  df %>%
    rename(Features = `Feature ID`)
}

## Merge both dataframes in rows where there's a match between features in both dataframes, filtering out columns I'm not interested in 
columns_to_keep <- c("Category", "Subsystem", "Role",  "Contig", "Start", "Stop")  # Define the columns to keep

merged_SA1 <- inner_join(SA1_df, features_SA1, by = "Features") %>% select(all_of(columns_to_keep))
merged_SA6 <- inner_join(SA6_df, features_SA6, by = "Features") %>% select(all_of(columns_to_keep))
merged_SB1 <- inner_join(SB1_df, features_SB1, by = "Features") %>% select(all_of(columns_to_keep))
merged_SB8 <- inner_join(SB8_df, features_SB8, by = "Features") %>% select(all_of(columns_to_keep))
merged_WA5 <- inner_join(WA5_df, features_WA5, by = "Features") %>% select(all_of(columns_to_keep))
merged_WB3 <- inner_join(WB3_df, features_WB3, by = "Features") %>% select(all_of(columns_to_keep))
merged_WB8 <- inner_join(WB8_df, features_WB8, by = "Features") %>% select(all_of(columns_to_keep))

## Save them as tsv files 
write_tsv(merged_SA1, "genome_mining_SA1.tsv")
write_tsv(merged_SA6, "genome_mining_SA6.tsv")
write_tsv(merged_SB1, "genome_mining_SB1.tsv")
write_tsv(merged_SB8, "genome_mining_SB8.tsv")
write_tsv(merged_WA5, "genome_mining_WA5.tsv")
write_tsv(merged_WB3, "genome_mining_WB3.tsv")
write_tsv(merged_WB8, "genome_mining_WB8.tsv")











