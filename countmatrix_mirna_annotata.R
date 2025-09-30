# Load libraries
library(dplyr)
library(readr)

# Input paths
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
mirdeep_file <- args[2]
output_file <- args[3]

# Read the counts matrix
counts <- read_tsv(counts_file)

# Read miRDeep2 file as plain text
all_lines <- readLines(mirdeep_file)

# Find line index where "tag id" header starts
header_line <- grep("^tag id\\t", all_lines)[1]

# Subset to actual data block (header + data lines)
mirdeep_data_lines <- all_lines[header_line:length(all_lines)]

# Write to temp file and read with read_tsv
tmp_file <- tempfile(fileext = ".tsv")
writeLines(mirdeep_data_lines, tmp_file)
mirdeep <- read_tsv(tmp_file, show_col_types = FALSE)

# Create annotation mapping
mir_map <- mirdeep %>%
  select(`tag id`, `mature miRBase miRNA`) %>%
  mutate(Geneid = paste0("known:", `tag id`),
         new_id = paste0(Geneid, "_", `mature miRBase miRNA`))

# Annotate the counts
annotated_counts <- counts %>%
  left_join(mir_map, by = "Geneid") %>%
  mutate(Geneid = ifelse(!is.na(new_id), new_id, Geneid)) %>%
  select(-new_id)

# Write output
#write_tsv(annotated_counts, output_file)
# Write output (remove extra columns from final file)
write_tsv(annotated_counts %>% select(-`tag id`, -`mature miRBase miRNA`), output_file)

