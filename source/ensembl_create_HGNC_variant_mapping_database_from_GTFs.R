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

# Major analysis point no.1a: Record the historical development of ENST IDs, retirement status of ENST ID as well as the latest ENST ID it occurs in, split by ENST ID
## detect latest release
latest_release_number <- tibble_all_GTFs_universal_ensg_id$ensembl_release_version %>% max

## ENST TRACKING TABLE: filter table for ENST ID and version 
tibble_ENST_ID_version_tracked_by_release <- tibble_all_GTFs_universal_ensg_id[, c("gene_name", "transcript_id", "ensembl_release_version", "genome_assembly", "universal_ensg_number")] %>% dplyr::distinct()

## ENST RETIREMENT TABLE: last release + retirement status
tibble_ENST_retirement_status <- tibble_ENST_ID_version_tracked_by_release %>% 
  dplyr::distinct() %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarise("release_last_seen_transcript_id" = max(ensembl_release_version), 
                   "retirement_status_transcript_id" = if (max(ensembl_release_version) == latest_release_number) {"active"} else if (max(ensembl_release_version) < latest_release_number) {"retired"} )

# Major analysis point no.1b: Record the historical development of ENSP IDs, retirement status of ENSP ID as well as the latest ENSP ID it occurs in, split by ENSP ID

## ENSP TRACKING TABLE: filter table for ENSP ID and version 
tibble_ENSP_ID_version_tracked_by_release <- tibble_all_GTFs_universal_ensg_id[, c("gene_name", "protein_id", "ensembl_release_version", "genome_assembly", "universal_ensg_number")] %>% dplyr::distinct() %>% .[!is.na(.$protein_id), ]

## ENSP RETIREMENT TABLE: last release + retirement status
tibble_ENSP_retirement_status <- tibble_ENSP_ID_version_tracked_by_release %>% 
  dplyr::distinct() %>%
  dplyr::group_by(protein_id) %>%
  dplyr::summarise("release_last_seen_protein_id" = max(ensembl_release_version), 
                   "retirement_status_protein_id" = if (max(ensembl_release_version) == latest_release_number) {"active"} else if (max(ensembl_release_version) < latest_release_number) {"retired"} )

# Major analysis point no.2: Write down the ENSG tree, as well as the associated gene_name, combined gene_name and the release version
tibble_historical_ENSG_tree_and_gene_names <- tibble_all_GTFs_universal_ensg_id[, c("universal_ensg_number", "gene_id", "gene_name", "ensg_gene_name_combined", "ensembl_release_version", "genome_assembly")] %>% dplyr::distinct() %>% dplyr::arrange(universal_ensg_number)

# Major analysis point no.3a: HGNC STABLE TRANSCRIPT ID MAPPING.
## to enumerate the HGNC stable transcript ID, split by UNIVERSAL ensg number
## if multiple gene names have been assigned to the gene tree in the past, then we take the most recent gene name used (that's why we need the ensembl_release_version.)
## IMPORTANT NOTE: as of ensembl_104, there does not exist ENST IDs which have been attributed to more than one HGNC gene name in the same version, which is very good news for us. in the future, we still have to check:
# test <- tibble_all_GTFs_universal_ensg_id %>% .[!is.na(.$universal_ensg_number) & !is.na(.$transcript_id), ] %>%
#   dplyr::distinct(transcript_id, gene_name, ensembl_release_version)
# test2 <- test %>% dplyr::group_by(transcript_id, ensembl_release_version) %>% dplyr::summarise("tally" = n()) %>% dplyr::arrange(desc(tally))

# unfortunately, *different* ENST IDs in the same gene tree can have different HGNC gene names in the same ensembl release.
## there's no guarantee that two ENSG IDs with unrelated histories may be merged into the same gene tree sometime in the future, thus destroying the ascending ENST order
## therefore we have to enumerate ENST IDs by ENSG IDs, *AND ADD A NUMBER BEFORE "enst" THAT INDICATES THE HISTORICAL ORDER IN WHICH THE ENSG ID HAS BEEN ADDED TO THE GENE TREE.*
list_tibble_all_ENST_IDs_split_by_gene_tree <- tibble_all_GTFs_universal_ensg_id %>% 
  .[!is.na(.$universal_ensg_number) & !is.na(.$transcript_id), ] %>%
  dplyr::distinct(gene_id, universal_ensg_number, transcript_id, ensg_gene_name_combined, transcript_version, ensembl_release_version, .keep_all = FALSE) %>% 
  dplyr::group_split(universal_ensg_number) %>%
  set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$universal_ensg_number %>% unique) %>% unlist)

plan(list(tweak(multiprocess, workers = 24),
          tweak(multiprocess, workers = 1)))

# fetch the table join data for each gene tree
# list_ENSG_ID_table_join_per_gene_tree <- purrr::imap(
#   .x = list_tibble_all_ENST_IDs_split_by_gene_tree %>% purrr::map(~.x$gene_id %>% unique),
#   .f = function(a1) {
#     
#     # DEBUG ###
#     # a1 <- list_tibble_all_ENST_IDs_split_by_gene_tree %>% purrr::map(~.x$gene_id %>% unique) %>% .[[1]]
#     ###########
#     
#     tibble_joined_by_combined_name_and_ensembl_release_version[tibble_joined_by_combined_name_and_ensembl_release_version$gene_id.x %in% a1 | tibble_joined_by_combined_name_and_ensembl_release_version$gene_id.y %in% a1, ]
#     
#   } )

# create HGNC stable transcript IDs
tibble_HGNC_stable_transcript_IDs <- furrr::future_imap(
  .x = list_tibble_all_ENST_IDs_split_by_gene_tree,
  .f = function(a1, a2) {
    
    # DEBUG ###
    # a2 <- list_tibble_all_ENST_IDs_split_by_gene_tree[[17819]]
    # a1 <- list_tibble_all_ENST_IDs_split_by_gene_tree[[76616]]
    # a1 <- list_tibble_all_ENST_IDs_split_by_gene_tree[[24]]
    # a1 <- list_tibble_all_ENST_IDs_split_by_gene_tree[[32]]
    # a1 <- list_tibble_all_ENST_IDs_split_by_gene_tree[[1]]
    ###########
    
    cat(a2, "\n")
    
    # generate historical order of ENSG IDs joining the gene tree.
    ## rejoin each ENSG ID, ensembl_release_version by the ensg_gene_name_combined. 
    ## This will help us determine the historical order in which each ENSG ID joined the gene tree.
    tibble_gene_tree_join_order0 <- dplyr::inner_join(
      a1[, c("gene_id", "ensembl_release_version", "ensg_gene_name_combined")],
      a1[, c("gene_id", "ensembl_release_version", "ensg_gene_name_combined")],
      by = c("ensg_gene_name_combined", "ensembl_release_version")
    ) %>% .[.$`gene_id.x` != .$`gene_id.y`, ]
    
    # for each ENSG ID, get the earliest release version where it joined the gene tree (due to common gene symbol)
    if (nrow(tibble_gene_tree_join_order0) > 0) {
      tibble_gene_tree_join_order1 <- tibble_gene_tree_join_order0 %>% dplyr::group_by(ensembl_release_version) %>% dplyr::summarise("gene_id" = c(gene_id.x, gene_id.y) %>% unique) %>% dplyr::group_by(gene_id) %>% dplyr::summarise("first_release_ENSG_joined" = min(ensembl_release_version)) %>% .[mixedorder(.$gene_id), ] %>% dplyr::arrange(first_release_ENSG_joined) %>% tibble::add_column("ENSG_historical_tree_order" = 1:nrow(.))
    } else if (nrow(tibble_gene_tree_join_order0) == 0) {
      tibble_gene_tree_join_order1 <- a1[, c("gene_id", "ensembl_release_version", "ensg_gene_name_combined")] %>% dplyr::group_by(gene_id) %>% dplyr::summarise("first_release_ENSG_joined" = min(ensembl_release_version)) %>% tibble::add_column("ENSG_historical_tree_order" = 1:nrow(.))
    }
    
    output_tibble0 <- dplyr::left_join(
      a1, 
      tibble_gene_tree_join_order1
    )
    
    # add HGNC_stable_transcript_ID
    ## first group by ENST ID + ENSG_historical_tree_order
    ## then add a column for the ensg_gene_name_combined that is found in the latest release
    ## this makes it resistant to changes in ENST to ENSG mapping, i.e. in the event that and ENST ID gets moved to another ENSG ID.
    output_tibble1 <- output_tibble0 %>%
      dplyr::group_by(transcript_id, ENSG_historical_tree_order) %>%
      dplyr::mutate("last_release_seen_gene_transcript_id" = max(ensembl_release_version))
    
    output_tibble0 <- dplyr::left_join(
      output_tibble1,
      output_tibble1[output_tibble1$ensembl_release_version == output_tibble1$last_release_seen_gene_transcript_id, c("gene_id", "transcript_id", "ensg_gene_name_combined")] %>% dplyr::rename("latest_ensg_gene_name_combined" = "ensg_gene_name_combined")
    )
    
    # finally add in the stable transcript number
    output_tibble1 <- output_tibble0 %>% 
      ungroup() %>% 
      group_split(gene_id) %>%
      purrr::imap(
        .x = .,
        .f = function(b1, b2) {
          
          # DEBUG ###
          # b1 <- output_tibble1 %>% ungroup() %>% group_split(gene_id) %>% .[[1]]
          ###########
          
          print(b2)
          
          output <- b1[mixedorder(b1$transcript_id), ]
          
          output %>% 
            dplyr::group_by(transcript_id) %>% 
            tibble::add_column("stable_transcript_number" = group_indices(.)) %>% 
            return
          
        } ) %>% dplyr::bind_rows() %>%
      dplyr::ungroup() %>%
      dplyr::mutate("hgnc_stable_transcript_ID" = paste(latest_ensg_gene_name_combined, "-", ENSG_historical_tree_order, "enst", stable_transcript_number, ".", transcript_version, sep = "") %>% gsub(pattern = "\\-1enst", replacement = "-enst"))
      
    return(output_tibble1)
    
  }, .progress = TRUE ) %>% rbindlist %>% as_tibble

# get rid of NA transcript_versions
tibble_HGNC_stable_transcript_IDs <- tibble_HGNC_stable_transcript_IDs %>% 
  dplyr::mutate("hgnc_stable_transcript_ID" = gsub(x = `hgnc_stable_transcript_ID`, pattern = ".NA$", replacement = ""))

# Major analysis point no.3b: HGNC STABLE PROTEIN ID MAPPING.
## to enumerate the HGNC stable protein ID, split by UNIVERSAL ensg number
## if multiple gene names have been assigned to the gene tree in the past, then we take the most recent gene name used (that's why we need the ensembl_release_version.)
## IMPORTANT NOTE: as of ensembl_104, there does not exist ENST IDs which have been attributed to more than one HGNC gene name in the same version, which is very good news for us. in the future, we still have to check:
# test <- tibble_all_GTFs_universal_ensg_id %>% .[!is.na(.$universal_ensg_number) & !is.na(.$protein_id), ] %>%
#   dplyr::distinct(protein_id, gene_name, ensembl_release_version)
# test2 <- test %>% dplyr::group_by(protein_id, ensembl_release_version) %>% dplyr::summarise("tally" = n()) %>% dplyr::arrange(desc(tally))

# unfortunately, *different* ENST IDs in the same gene tree can have different HGNC gene names in the same ensembl release.
## there's no guarantee that two ENSG IDs with unrelated histories may be merged into the same gene tree sometime in the future, thus destroying the ascending ENST order
## therefore we have to enumerate ENST IDs by ENSG IDs, *AND ADD A NUMBER BEFORE "enst" THAT INDICATES THE HISTORICAL ORDER IN WHICH THE ENSG ID HAS BEEN ADDED TO THE GENE TREE.*
list_tibble_all_ENSP_IDs_split_by_gene_tree <- tibble_all_GTFs_universal_ensg_id %>% 
  .[!is.na(.$universal_ensg_number) & !is.na(.$protein_id), ] %>%
  dplyr::distinct(gene_id, universal_ensg_number, protein_id, ensg_gene_name_combined, protein_version, ensembl_release_version, .keep_all = FALSE) %>% 
  dplyr::group_split(universal_ensg_number) %>%
  set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$universal_ensg_number %>% unique) %>% unlist)

plan(list(tweak(multiprocess, workers = 24),
          tweak(multiprocess, workers = 1)))

# fetch the table join data for each gene tree
# list_ENSG_ID_table_join_per_gene_tree <- purrr::imap(
#   .x = list_tibble_all_ENSP_IDs_split_by_gene_tree %>% purrr::map(~.x$gene_id %>% unique),
#   .f = function(a1) {
#     
#     # DEBUG ###
#     # a1 <- list_tibble_all_ENSP_IDs_split_by_gene_tree %>% purrr::map(~.x$gene_id %>% unique) %>% .[[1]]
#     ###########
#     
#     tibble_joined_by_combined_name_and_ensembl_release_version[tibble_joined_by_combined_name_and_ensembl_release_version$gene_id.x %in% a1 | tibble_joined_by_combined_name_and_ensembl_release_version$gene_id.y %in% a1, ]
#     
#   } )

# create HGNC stable protein IDs
tibble_HGNC_stable_protein_IDs <- furrr::future_imap(
  .x = list_tibble_all_ENSP_IDs_split_by_gene_tree,
  .f = function(a1, a2) {
    
    # DEBUG ###
    # a2 <- list_tibble_all_ENSP_IDs_split_by_gene_tree[[17819]]
    # a1 <- list_tibble_all_ENSP_IDs_split_by_gene_tree[[76616]]
    # a1 <- list_tibble_all_ENSP_IDs_split_by_gene_tree[[24]]
    # a1 <- list_tibble_all_ENSP_IDs_split_by_gene_tree[[32]]
    # a1 <- list_tibble_all_ENSP_IDs_split_by_gene_tree[[1]]
    ###########
    
    cat(a2, "\n")
    
    # generate historical order of ENSG IDs joining the gene tree.
    ## rejoin each ENSG ID, ensembl_release_version by the ensg_gene_name_combined. 
    ## This will help us determine the historical order in which each ENSG ID joined the gene tree.
    tibble_gene_tree_join_order0 <- dplyr::inner_join(
      a1[, c("gene_id", "ensembl_release_version", "ensg_gene_name_combined")],
      a1[, c("gene_id", "ensembl_release_version", "ensg_gene_name_combined")],
      by = c("ensg_gene_name_combined", "ensembl_release_version")
    ) %>% .[.$`gene_id.x` != .$`gene_id.y`, ]
    
    # for each ENSG ID, get the earliest release version where it joined the gene tree (due to common gene symbol)
    if (nrow(tibble_gene_tree_join_order0) > 0) {
      tibble_gene_tree_join_order1 <- tibble_gene_tree_join_order0 %>% dplyr::group_by(ensembl_release_version) %>% dplyr::summarise("gene_id" = c(gene_id.x, gene_id.y) %>% unique) %>% dplyr::group_by(gene_id) %>% dplyr::summarise("first_release_ENSG_joined" = min(ensembl_release_version)) %>% .[mixedorder(.$gene_id), ] %>% dplyr::arrange(first_release_ENSG_joined) %>% tibble::add_column("ENSG_historical_tree_order" = 1:nrow(.))
    } else if (nrow(tibble_gene_tree_join_order0) == 0) {
      tibble_gene_tree_join_order1 <- a1[, c("gene_id", "ensembl_release_version", "ensg_gene_name_combined")] %>% dplyr::group_by(gene_id) %>% dplyr::summarise("first_release_ENSG_joined" = min(ensembl_release_version)) %>% tibble::add_column("ENSG_historical_tree_order" = 1:nrow(.))
    }
    
    output_tibble0 <- dplyr::left_join(
      a1, 
      tibble_gene_tree_join_order1
    )
    
    # add HGNC_stable_protein_ID
    ## first group by ENST ID + ENSG_historical_tree_order
    ## then add a column for the ensg_gene_name_combined that is found in the latest release
    ## this makes it resistant to changes in ENST to ENSG mapping, i.e. in the event that and ENST ID gets moved to another ENSG ID.
    output_tibble1 <- output_tibble0 %>%
      dplyr::group_by(protein_id, ENSG_historical_tree_order) %>%
      dplyr::mutate("last_release_seen_gene_protein_id" = max(ensembl_release_version))
    
    output_tibble0 <- dplyr::left_join(
      output_tibble1,
      output_tibble1[output_tibble1$ensembl_release_version == output_tibble1$last_release_seen_gene_protein_id, c("gene_id", "protein_id", "ensg_gene_name_combined")] %>% dplyr::rename("latest_ensg_gene_name_combined" = "ensg_gene_name_combined")
    )
    
    # finally add in the stable protein number
    output_tibble1 <- output_tibble0 %>% 
      ungroup() %>% 
      group_split(gene_id) %>%
      purrr::imap(
        .x = .,
        .f = function(b1, b2) {
          
          # DEBUG ###
          # b1 <- output_tibble1 %>% ungroup() %>% group_split(gene_id) %>% .[[1]]
          ###########
          
          print(b2)
          
          output <- b1[mixedorder(b1$protein_id), ]
          
          output %>% 
            dplyr::group_by(protein_id) %>% 
            tibble::add_column("stable_protein_number" = group_indices(.)) %>% 
            return
          
        } ) %>% dplyr::bind_rows() %>%
      dplyr::ungroup() %>%
      dplyr::mutate("hgnc_stable_protein_ID" = paste(latest_ensg_gene_name_combined, "-", ENSG_historical_tree_order, "ensp", stable_protein_number, ".", protein_version, sep = "") %>% gsub(pattern = "\\-1ensp", replacement = "-ensp"))
    
    return(output_tibble1)
    
  }, .progress = TRUE ) %>% rbindlist %>% as_tibble

# get rid of NA protein_versions
tibble_HGNC_stable_protein_IDs <- tibble_HGNC_stable_protein_IDs %>% 
  dplyr::mutate("hgnc_stable_protein_ID" = gsub(x = `hgnc_stable_protein_ID`, pattern = ".NA$", replacement = ""))

# TOTAL MAPPING TABLE
## we are going to combine the information from all 1 + 1 + 2 tables into one
tibble_total_mapping_table_transcripts <- dplyr::left_join(tibble_HGNC_stable_transcript_IDs, tibble_ENST_retirement_status)
tibble_total_mapping_table_proteins <- dplyr::left_join(tibble_HGNC_stable_protein_IDs, tibble_ENSP_retirement_status %>% dplyr::rename("release_last_seen_protein_id" = "release_last_seen", "retirement_status_protein_id" = "retirement_status"))

if (nrow(tibble_total_mapping_table_transcripts) != nrow(tibble_HGNC_stable_transcript_IDs)) {
  stop("imperfect join. check whether there may be some ENST IDs that have been assigned to more than one ENSG ID in the past.")
}

if (nrow(tibble_total_mapping_table_proteins) != nrow(tibble_HGNC_stable_protein_IDs)) {
  stop("imperfect join. check whether there may be some ENSP IDs that have been assigned to more than one ENSG ID in the past.")
}

# ANNOTATE ALL RELEASES
# tibble_all_annotated_releases <- dplyr::left_join(tibble_all_GTFs_universal_ensg_id %>% dplyr::mutate("id" = 1:nrow(.), .before = 1), tibble_total_mapping_table)
tibble_all_annotated_releases <- dplyr::left_join(tibble_all_GTFs_universal_ensg_id, tibble_total_mapping_table_transcripts) %>%
  dplyr::left_join(., tibble_total_mapping_table_proteins)

# list-ify and save
list_all_annotated_releases <- tibble_all_annotated_releases %>%
  dplyr::group_split(ensembl_release_version) %>%
  set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$ensembl_release_version %>% unique) %>% unlist)

plan(list(tweak(multiprocess, workers = 16),
          tweak(multiprocess, workers = 4)))

# for each table in the list, if there is no "transcript" in the column type, we add it in
list_all_annotated_releases_fixed <- furrr::future_map(
  .x = list_all_annotated_releases,
  .f = function(a1) {
    
    # DEBUG ###
    # a1 <- list_all_annotated_releases$Homo_sapiens.GRCh38.76.gtf
    ###########
    
    if ("transcript" %in% (a1$type %>% unique)) {
      
      return(a1)
      
    } else {
      
      # list-ify the tibble by transcript_id
      list_a1_split_by_transcript_id <- a1 %>% 
        dplyr::group_split(transcript_id) %>%
        set_names(nm = purrr::map(.x = ., .f = ~.x$transcript_id %>% unique) %>% unlist)
      
      # create entries for all transcript_ids
      tibble_transcript_entries <- furrr::future_map2(
        .x = list_a1_split_by_transcript_id,
        .y = names(list_a1_split_by_transcript_id),
        .f = function(b1, b2) {
          
          # DEBUG ###
          # b1 <- list_a1_split_by_transcript_id[[1]]
          # b2 <- names(list_a1_split_by_transcript_id)[[1]]
          ###########
          
          return(
            
            dplyr::bind_rows(
              b1,
              b1[1, ] %>% dplyr::mutate(
                "start" = min(`start`),
                "end" = max(`end`),
                "width" = max(`end`) - min(`start`) + 1,
                "type" = "transcript",
                "exon_number" = NA,
                "protein_id" = NA
              )
            )
            
          )
          
        } )
      
      return(tibble_transcript_entries %>% rbindlist(use.names = TRUE, fill = TRUE) %>% as_tibble)
      
    }
    
  }, .progress = TRUE )

# save ALL annotated releases
furrr::future_map2(
  .x = list_all_annotated_releases_fixed,
  .y = names(list_all_annotated_releases_fixed),
  .f = ~write.table(x = .x, file = paste(output_dir, "annotated_ensembl_gtf_release_", .y, ".txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE), 
  .progress = TRUE )

# 1. save ENST retirement tracking
write.table(x = tibble_ENST_retirement_status, file = paste(output_dir, "latest_ENST_retirement_status.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 2. save gene tree tracking
write.table(x = tibble_historical_ENSG_tree_and_gene_names, file = paste(output_dir, "latest_gene_tree_tracking.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 3. save HGNC stable transcript mapping table
write.table(x = tibble_total_mapping_table, file = paste(output_dir, "latest_mapping_table.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# add information about latest update
cat(paste("Date of last update: ", date(), "\nLatest Ensembl release version: ", latest_release_number, "\n\nSystem info of machine that retrieved last update:\n\n", sep = ""), file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""))

write.table(Sys.info() %>% tibble::enframe(name = "id", value = "value"), file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""), sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

write("\n", file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""), append = TRUE)

write.table(R.version %>% unlist(use.names = TRUE) %>% tibble::enframe(name = "id", value = "value"), file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""), sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

# finish counting
tictoc::toc()

q()

