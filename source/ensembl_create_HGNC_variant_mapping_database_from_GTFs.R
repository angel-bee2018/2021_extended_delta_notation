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
# input_GTF_dir <- "/mnt/LTS/projects/2020_isoform_nomenclature/ensembl_GTF_dump/"
# output_dir <- "/mnt/LTS/projects/2020_isoform_nomenclature/ENST_to_HGNC_stable_variant_mapping/"
# ncores <- "6x4"
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
options(mc.cores = 66)

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

list_all_GTFs <- list_all_GTFs[names(list_all_GTFs) %>% mixedsort]

# ADD IN STABLE UNIVERSAL ensg NUMBER
## SOMETIMES A GENE NAME CAN HAVE MORE THAN ONE ENSG ID
## we MUST know how to group ENSG and gene names together in order to create the HGNC stable variant IDs.
## in other words, we have to know which genes are the same
## go up the releases by ENSG IDs. match each ENSG ID to the next release. for each release, assign the same universal ENSG number to the ENSG IDs that are attributed to the same gene name. if there is no gene name, use the ENSG ID.

# Return values: 1. each release annotated with the universal_ensg_number, 2. the tibble that maps ALL ensg_gene_name_combined IDs to the gene_name and ENSG ID (can be done later)

## add a dummy column specifying the dummy gene name = HGNC gene name + ENSG ID if gene name is not available.
list_all_GTFs <- purrr::map(
  .x = list_all_GTFs,
  .f = function(a1) {
    
    a1$ensg_gene_name_combined <- a1$gene_name
    
    a1[is.na(a1$ensg_gene_name_combined), "ensg_gene_name_combined"] <- a1[is.na(a1$ensg_gene_name_combined), "gene_id"]
    
    return(a1)
    
  } )

tibble_all_GTFs <- list_all_GTFs %>% rbindlist(fill = TRUE, use.names = TRUE) %>% as_tibble

# join by combined name and ensembl_release_version because that's the only known commonality between apparently unrelated ENSG ids
tibble_joined_by_combined_name_and_ensembl_release_version <- dplyr::full_join(
  tibble_all_GTFs[, c("ensg_gene_name_combined", "gene_id", "ensembl_release_version")] %>% unique,
  tibble_all_GTFs[, c("ensg_gene_name_combined", "gene_id", "ensembl_release_version")] %>% unique,
  by = c("ensg_gene_name_combined", "ensembl_release_version")
)

# this freezes interactions between apparently unrelated ENSG IDs, as their only point of contact is same combined names in the same assembly version.
# all we have to do now is pull out each ENSG ID from the list and collapse. 
# we should get a list of vectors - each element corresponds to all the ENSG IDs that describes the same underlying gene.

# equivalent to hard thresholding cluster detection problem
## do not unique-ify because you'll miss out on the well-behaved ENSG IDs
tibble_edge_table <- tibble_joined_by_combined_name_and_ensembl_release_version[, c("gene_id.x", "gene_id.y")] %>% unique %>% setNames(nm = c("a", "b"))

# iterative grow from each node 
vector_all_available_nodes <- c(tibble_edge_table$a, tibble_edge_table$b) %>% unique

# create logical signifying whether all nodes have been considered
logical_all_nodes_considered <- FALSE

# work through the list of all availbble nodes
vector_nodes_remaining <- vector_all_available_nodes

# pre-allocate a tibble with seed value and members of association
tibble_association_table <- tibble("seed" = character(0),
                                   "members_of_association" = character(0),
                                   "association_size" = integer(0))

temp_tibble_edge_table <- tibble_edge_table

while (logical_all_nodes_considered == FALSE) {
  
  current_node <- vector_nodes_remaining[1]
  
  cat("current_node: ", current_node, "\n")
  
  # iteratively grow an association out from the current node until it can't get any bigger
  iteration_difference <- 1
  # the list of association members grows as the tree moves outwards
  vector_current_association_members <- current_node
  while (iteration_difference > 0 & nrow(temp_tibble_edge_table) > 0) {
    
    previous_number_of_members_in_association <- vector_current_association_members %>% length
    
    # check source nodes
    vector_logical_temp_tibble_indices_source_node_hit <- temp_tibble_edge_table$a %in% vector_current_association_members
    
    # check target nodes
    vector_logical_temp_tibble_indices_target_node_hit <- temp_tibble_edge_table$b %in% vector_current_association_members
    
    # find all interactors of current association members
    vector_all_interactors_of_current_association_members <- c(
      temp_tibble_edge_table[vector_logical_temp_tibble_indices_source_node_hit, "b"] %>% unlist,
      temp_tibble_edge_table[vector_logical_temp_tibble_indices_target_node_hit, "a"] %>% unlist)
    
    # remove already traversed edges from the temp edge table
    temp_tibble_edge_table <- temp_tibble_edge_table[!vector_logical_temp_tibble_indices_source_node_hit & !vector_logical_temp_tibble_indices_target_node_hit, ]
    
    cat("temp tibble rows: ", nrow(temp_tibble_edge_table), "\n")
    
    vector_current_association_members <- c(vector_current_association_members, vector_all_interactors_of_current_association_members) %>% unique
    
    current_number_of_members_in_association <- vector_current_association_members %>% length
    
    iteration_difference <- current_number_of_members_in_association - previous_number_of_members_in_association
    
    cat("iteration_difference: ", iteration_difference, "\n")
    
  }
  
  tibble_association_table <- tibble_association_table %>% 
    tibble::add_row(tibble("seed" = current_node %>% as.character, "members_of_association" = vector_current_association_members %>% paste(collapse = ",") %>% as.character, "association_size" = current_number_of_members_in_association))
  
  vector_nodes_remaining <- vector_nodes_remaining[!vector_nodes_remaining %in% vector_current_association_members]
  
  cat("number of nodes remaining: ", length(vector_nodes_remaining), "\n")
  
  logical_all_nodes_considered <- length(vector_nodes_remaining) == 0
  
  flush.console() 
  
}

# retrieve unique gene mapping 
tibble_tracked_ENSG_id_trees <- tibble_association_table %>% 
  tibble::add_column("universal_ensg_number" = 1:nrow(.))

tibble_tracked_ENSG_id_trees <- tibble_tracked_ENSG_id_trees %>% dplyr::mutate("members_of_association" = strsplit(x = members_of_association, split = ","))

## unnest
tibble_tracked_ENSG_id_trees <- tibble_tracked_ENSG_id_trees %>% unnest(cols = members_of_association)

# tack on the universal ENSG number to the master tibble
tibble_all_GTFs_universal_ensg_id <- dplyr::left_join(tibble_all_GTFs, tibble_tracked_ENSG_id_trees %>% dplyr::rename("gene_id" = "members_of_association") %>% dplyr::select(-seed, -association_size))

# Major analysis point no.1: Record the historical development of ENST IDs, retirement status of ENST ID as well as the latest ENST ID it occurs in, split by ENST ID
## detect latest release
latest_release_number <- tibble_all_GTFs_universal_ensg_id$ensembl_release_version %>% max

## ENST TRACKING TABLE: filter table for ENST ID and version 
tibble_ENST_ID_version_tracked_by_release <- tibble_all_GTFs_universal_ensg_id[, c("gene_name", "transcript_id", "ensembl_release_version", "genome_assembly", "universal_ensg_number")] %>% dplyr::distinct()

## ENST RETIREMENT TABLE: last release + retirement status
tibble_ENST_retirement_status <- tibble_ENST_ID_version_tracked_by_release %>% 
  dplyr::distinct() %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarise("release_last_seen" = max(ensembl_release_version), 
                   "retirement_status" = if (max(ensembl_release_version) == latest_release_number) {"active"} else if (max(ensembl_release_version) < latest_release_number) {"retired"} )

# Major analysis point no.2: Write down the ENSG tree, as well as the associated gene_name, combined gene_name and the release version
tibble_historical_ENSG_tree_and_gene_names <- tibble_all_GTFs_universal_ensg_id[, c("universal_ensg_number", "gene_id", "gene_name", "ensg_gene_name_combined", "ensembl_release_version", "genome_assembly")] %>% dplyr::distinct() %>% dplyr::arrange(universal_ensg_number)

# Major analysis point no.3: HGNC STABLE VARIANT ID MAPPING.
## to enumerate the HGNC stable variant ID, split by UNIVERSAL ensg number
## if multiple gene names have been assigned to the gene tree in the past, then we take the most recent gene name used (that's why we need the ensembl_release_version.)
## IMPORTANT NOTE: as of ensembl_104, there does not exist ENST IDs which have been attributed to more than one HGNC gene name in the same version, which is very good news for us. in the future, we still have to check:
# test <- tibble_all_GTFs_universal_ensg_id %>% .[!is.na(.$universal_ensg_number) & !is.na(.$transcript_id), ] %>%
#   dplyr::distinct(transcript_id, gene_name, ensembl_release_version)
# test2 <- test %>% dplyr::group_by(transcript_id, ensembl_release_version) %>% dplyr::summarise("tally" = n()) %>% dplyr::arrange(desc(tally))

# unfortunately, *different* ENST IDs in the same gene tree can have different HGNC gene names in the same ensembl release.
# WE CAN LITERALLY JUST PUT IT IN.
# the HGNC_stable_variant_ID will still be completely correct for the gene tree - just different ENST IDs will have their own HGNC gene symbol.
list_tibble_all_ENST_IDs_split_by_ENSG_ID <- tibble_all_GTFs_universal_ensg_id %>% 
  .[!is.na(.$universal_ensg_number) & !is.na(.$transcript_id), ] %>%
  dplyr::distinct(universal_ensg_number, transcript_id, ensg_gene_name_combined, transcript_version, ensembl_release_version) %>% 
  dplyr::group_split(universal_ensg_number) %>%
  set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$universal_ensg_number %>% unique) %>% unlist)

options(mc.cores = 16)

# create HGNC stable variant IDs
tibble_HGNC_stable_variant_IDs <- furrr::future_imap(
  .x = list_tibble_all_ENST_IDs_split_by_ENSG_ID,
  .f = function(a1, a2) {
    
    # DEBUG ###
    # a1 <- list_tibble_all_ENST_IDs_split_by_ENSG_ID$`1`
    ###########
    
    cat(a2, "\n")
    
    # sort by ENST ID
    sorted_tibble <- a1[mixedorder(a1$transcript_id), ]
    
    # add HGNC_stable_variant_ID
    ## first group by ENST ID
    ## then add a column for the ensg_gene_name_combined that is found in the latest release
    output_tibble <- sorted_tibble %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::mutate("release_last_seen" = max(ensembl_release_version))
    
    output_tibble <- dplyr::left_join(output_tibble,
                                      output_tibble[output_tibble$ensembl_release_version == output_tibble$release_last_seen, c("transcript_id", "ensg_gene_name_combined")] %>% dplyr::rename("latest_ensg_gene_name_combined" = "ensg_gene_name_combined")
                                      )
    
    output_tibble <- output_tibble[mixedorder(output_tibble$transcript_id), ]
    
    # finally add in the stable variant number
    output_tibble <- output_tibble %>% 
      tibble::add_column("stable_variant_number" = group_indices(.)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate("hgnc_stable_variant_ID" = paste(latest_ensg_gene_name_combined, "-enst", stable_variant_number, ".", transcript_version, sep = ""))
    
    return(output_tibble)
    
  }, .progress = TRUE ) %>% rbindlist %>% as_tibble

# get rid of NA transcript_versions
tibble_HGNC_stable_variant_IDs <- tibble_HGNC_stable_variant_IDs %>% 
  dplyr::mutate("hgnc_stable_variant_ID" = gsub(x = `hgnc_stable_variant_ID`, pattern = ".NA$", replacement = ""))

# TOTAL MAPPING TABLE
## we are going to combine the information from all 3 tables into one
tibble_total_mapping_table <- dplyr::left_join(tibble_HGNC_stable_variant_IDs, tibble_ENST_retirement_status)

if (nrow(tibble_total_mapping_table) != nrow(tibble_HGNC_stable_variant_IDs)) {
  stop("imperfect join. check whether there may be some ENST IDs that have been assigned to more than one ENSG ID in the past.")
}

# ANNOTATE ALL RELEASES
tibble_all_annotated_releases <- dplyr::left_join(tibble_all_GTFs_universal_ensg_id, tibble_total_mapping_table)

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

# 1. save ENST retirement tracking
write.table(x = tibble_ENST_retirement_status, file = paste(output_dir, "latest_ENST_retirement_status.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 2. save gene tree tracking
write.table(x = tibble_historical_ENSG_tree_and_gene_names, file = paste(output_dir, "latest_gene_tree_tracking.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 3. save HGNC stable variant mapping table
write.table(x = tibble_total_mapping_table, file = paste(output_dir, "latest_mapping_table.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# add information about latest update
cat(paste("Date of last update: ", date(), "\nLatest Ensembl release version: ", latest_release_number, "\n\nSystem info of machine that retrieved last update:\n\n", sep = ""), file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""))

write.table(Sys.info() %>% tibble::enframe(name = "id", value = "value"), file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""), sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

write("\n", file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""), append = TRUE)

write.table(R.version %>% unlist(use.names = TRUE) %>% tibble::enframe(name = "id", value = "value"), file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""), sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

# finish counting
tictoc::toc()

q()

