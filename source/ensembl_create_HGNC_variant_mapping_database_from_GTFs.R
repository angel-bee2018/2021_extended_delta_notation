
script_description <- "# CREATE HGNC VARIANT MAPPING FILE ###
This script should be run in the folder that contains all historical GTFs.
Input: multiple GTFs
Output: Mapping between ENST IDs and HGNC variant numbers and retirement status"

# print the arguments received by the R script
cat("Arguments input:", commandArgs(), sep = "\n")
args = 
  commandArgs(trailingOnly = TRUE)
cat(args)
cat("number of arguments specified:", length(args))

# SET ENVIRONMENT ##########
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(c("seqinr", "tidyverse", "purrr", "dplyr", "rtracklayer", "data.table", "furrr", "RhpcBLASctl", "optparse"))

library(tidyverse)
library(furrr)
options(future.globals.maxSize = 30000000000, future.fork.enable = TRUE)
library(rtracklayer)
library(data.table)
library(optparse)
library(gtools)

library(tictoc)
# start counting execution time of the whole script
tictoc::tic("Overall execution time")

# manage arguments
list_input_arg_info = list(
  "1" = make_option(c("-D", "--input_GTF_dir"), type = "character", default = NULL, 
                    help = "Compulsory. Directory that contains all the historical GTFs.", metavar = "character"),
  "2" = make_option(c("-O", "--output_dir"), type = "character", default = NULL, 
                    help = "Compulsory. Directory that the mapping info gets outputted into.", metavar = "character"),
  "3" = make_option(c("-C", "--ncores"), type = "character", default = 0, 
                    help = "Optional. Number of cores to use. possible inputs: numbers 1 to any integer. By default, uses all cores (ncores = 0). If a single number is specified, it will just tell future to loop thru chromosomes in parallel using the specified core count. If numberxnumber for example 7x4 then 28 cores will be used. 7 for chromosomes and 4 for inside each chromosome.", metavar = "character")
)

input_arg_info <- OptionParser(option_list = list_input_arg_info, description = script_description)
input_args <- input_arg_info %>% parse_args

# check if the input arguments are O.K
if (list(input_args$input_GTF_dir, input_args$output_dir) %>% lapply(is.null) %>% unlist %>% any == TRUE) {
  
  print_help(input_arg_info)
  
  stop("Make sure you entered all arguments correctly", call. = FALSE)
  
}

input_GTF_dir <- input_args$input_GTF_dir
output_dir <- input_args$output_dir
ncores <- input_args$ncores

# DEBUG ########
#input_GTF_dir <- "/mnt/LTS/projects/2020_isoform_nomenclature/ensembl_GTF_dump/"
#output_dir <- "/mnt/LTS/projects/2020_isoform_nomenclature/ENST_to_HGNC_stable_variant_mapping/"
#ncores <- "8x4"
################

cat("input_GTF_dir:", input_GTF_dir, "\n")
cat("output_dir:", output_dir, "\n")

# manage parrallellisation rrlllRll

if (grepl(x = ncores, pattern = "x") == FALSE) {
  
  if (ncores != 0) {
    number_of_workers <- ncores
    cat(future::availableCores(), "cores will be used\n")
  } else {
    number_of_workers <- future::availableCores()
    cat(future::availableCores(), "cores will be used\n")
  } 
  
} else if (grepl(x = ncores, pattern = "x") == TRUE) {
  
  plan(list(tweak(multiprocess, workers = ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert, gc = TRUE), 
            tweak(multiprocess, workers = ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert, gc = TRUE))
  )
  
  cat((ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert) * (ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert), "cores will be used in total\n")
  cat("first layer:", ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert, "cores\n")
  cat("second layer:", ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert, "cores\n")
  
}

# BEGIN EXECUTION ###

plan(multiprocess)
options(mc.cores = 72)

# import all GTFs
list_all_GTFs <- furrr::future_map(
  .x = list.files(input_GTF_dir),
  .f = ~rtracklayer::import(paste(input_GTF_dir, .x, sep = "")) %>% as_tibble,
  .progress = TRUE )

# get vector of pure release numbers
vector_release_numbers <- gsub(x = list.files(input_GTF_dir), pattern = ".*\\.(\\d+)\\.gtf", replacement = "\\1") %>% type.convert

names(list_all_GTFs) <- vector_release_numbers

# get vector of assembly numbers
vector_assembly_numbers <- gsub(x = list.files(input_GTF_dir), pattern = ".*\\.(.*)\\.(\\d+)\\.gtf", replacement = "\\1")

setDTthreads(66)

# add column to each GTF indicating assembly and release version and tibblise
list_all_GTFs <- purrr::pmap(
  .l = list(
    "a1" = list_all_GTFs,
    "a2" = vector_release_numbers,
    "a3" = vector_assembly_numbers
  ),
  .f = function(a1, a2, a3) {
    
    return(a1 %>% add_column("ensembl_release_version" = a2, "genome_assembly" = a3))
    
  } )

tibble_all_GTFs <- list_all_GTFs %>% rbindlist(fill = TRUE, use.names = TRUE) %>% as_tibble

# tibble_all_GTFs[tibble_all_GTFs$transcript_version %>% is.na == FALSE, ] %>% .$ensembl_release_version %>% min

# to record the historical development of ENST IDs, retirement status of ENST ID as well as the latest ENST ID it occurs in, split by ENST ID

## detect latest release
latest_release_number <- tibble_all_GTFs$ensembl_release_version %>% max

## ENST TRACKING TABLE: filter table for ENST ID and version 
tibble_ENST_ID_version_tracked_by_release <- tibble_all_GTFs[, c("gene_name", "transcript_id", "ensembl_release_version", "genome_assembly")] %>% dplyr::distinct()

## ENST RETIREMENT TABLE: last release + retirement status
tibble_ENST_retirement_status <- tibble_ENST_ID_version_tracked_by_release %>% 
  dplyr::distinct() %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarise("release_last_seen" = max(ensembl_release_version), 
                   "retirement_status" = if (max(ensembl_release_version) == latest_release_number) {"active"} else if (max(ensembl_release_version) < latest_release_number) {"retired"} )

# to enumerate the HGNC stable variant ID, split by ENSG ID
list_tibble_all_ENST_IDs_split_by_ENSG_ID <- tibble_all_GTFs %>% 
  dplyr::distinct(gene_id, transcript_id, gene_name, transcript_version) %>% 
  .[!is.na(.$gene_id) & !is.na(.$transcript_id), ] %>%
  dplyr::group_split(gene_id) %>%
  set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$gene_id %>% unique) %>% unlist)

options(mc.cores = 10)

# HGNC STABLE VARIANT ID MAPPING TABLE: create HGNC stable variant IDs
tibble_HGNC_stable_variant_IDs <- furrr::future_map(
  .x = list_tibble_all_ENST_IDs_split_by_ENSG_ID,
  .f = function(a1) {
    
    # DEBUG ###
    # a1 <- list_tibble_all_ENST_IDs_split_by_ENSG_ID$ENSG00000229425
    ###########
    
    # sort by ENST ID
    sorted_tibble <- a1[mixedorder(a1$transcript_id), ]
    
    # IN THE CASE OF BLANK gene_ids - USE THE ENSG. jeez.
    sorted_tibble[is.na(sorted_tibble$gene_name), "gene_name"] <- sorted_tibble[is.na(sorted_tibble$gene_name), "gene_id"]
    
    # add HGNC_stable_variant_ID
    output_tibble <- sorted_tibble %>%
      dplyr::group_by(gene_id, transcript_id) %>% 
      tibble::add_column("stable_variant_number" = group_indices(.)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate("hgnc_stable_variant_ID" = paste(sorted_tibble$gene_name, "-enst", stable_variant_number, ".", sorted_tibble$transcript_version, sep = ""))
    
  }, .progress = TRUE ) %>% rbindlist %>% as_tibble

# get rid of NA transcript_versions
tibble_HGNC_stable_variant_IDs <- tibble_HGNC_stable_variant_IDs %>% 
  dplyr::mutate("hgnc_stable_variant_ID" = gsub(x = `hgnc_stable_variant_ID`, pattern = ".NA$", replacement = ""))

# TOTAL MAPPING TABLE
## we are going to combine the information from all 3 tables into one
tibble_total_mapping_table <- dplyr::left_join(tibble_HGNC_stable_variant_IDs, tibble_ENST_retirement_status, by = c("transcript_id"))

# ANNOTATE ALL RELEASES
tibble_all_annotated_releases <- dplyr::left_join(tibble_all_GTFs, tibble_total_mapping_table)

# list-ify and save
list_all_annotated_releases <- tibble_all_annotated_releases %>%
  dplyr::group_split(ensembl_release_version) %>%
  set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$ensembl_release_version %>% unique) %>% unlist)

# save ALL annotated releases
furrr::future_map2(
  .x = list_all_annotated_releases,
  .y = names(list_all_annotated_releases),
  .f = ~write.table(x = .x, file = paste(output_dir, "annotated_ensembl_gtf_release_", .y, ".txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE), 
  .progress = TRUE )

# save mapping table
write.table(x = tibble_total_mapping_table, file = paste(output_dir, "latest_mapping_table.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# add information about latest update
cat(paste("Date of last update: ", date(), "\nLatest Ensembl release version: ", latest_release_number, "\n\nSystem info of machine that retrieved last update:\n\n", sep = ""), file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""))

write.table(Sys.info() %>% tibble::enframe(name = "id", value = "value"), file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""), sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

write("\n", file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""), append = TRUE)

write.table(R.version %>% unlist(use.names = TRUE) %>% tibble::enframe(name = "id", value = "value"), file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""), sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

# finish counting
tictoc::toc()

q()

