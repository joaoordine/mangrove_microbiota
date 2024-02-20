# Preprocessing 

## Installing required R packages and activating them
install.packages("tidyverse")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ShortRead")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")

library(ShortRead)
library(dada2)
library(tidyverse)

## Appending `Sample´ in all fastq files to get a pattern
for file in ID*fastq; do mv "$file" "Sample-$file"; done

## Quality control - before trimming/clipping 
mkdir fastqc_before
fastqc -o fastqc_before *.fastq
multiqc fastqc_before
mv multiqc_report.html multiqc_report_before.html

## Defining parameters for trimminging and filtering
minLength <- 1200 #Removes reads shorter than this length. Minimum length is enforced AFTER trimming
maxLength <- 1800 #Removes reads longer than this length. Maximum length is enforced BEFORE trimming 
trimLeft <- 100 #The number of nucleotides to remove from the start of each read. Must cover the primer length
trimRight <- 100 #The number of nucleotides to remove from the end of each read. Must cover the primer length

## Handling the files
path <- getwd()
files <- list.files(pattern = "Sample") # create object for files based on patter 
sample.names <- tools::file_path_sans_ext(basename(files))  # remove extension
rawFiles <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
filtFiles <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq"))

## Sending filtered files to output directory
### I previously created the directory called 'filtered' and saved the table 
out <- filterAndTrim(rawFiles, filtFiles, trimLeft = trimLeft, 
                     trimRight = trimRight, maxLen = maxLength, 
                     minLen = minLength,  truncQ = 0, 
                     compress = FALSE)
out
write.table(out, file = "read_count_QC.txt", sep = "\t", row.names = TRUE)

## Quality control - after trimming/clipping 
mkdir fastqc_after # I created this directory inside the `filtered`one 
fastqc -o fastqc_after *.fastq
multiqc fastqc_after
mv multiqc_report.html multiqc_report_after.html
