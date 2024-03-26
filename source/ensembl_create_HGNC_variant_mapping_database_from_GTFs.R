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
                    help = "Optional. Number of cores to use. possible inputs: numbers 1 to any integer. By default, uses all cores (ncores = 0). If a single number is specified, it will just tell future to loop thru chromosomes in parallel using the specified core count. If numberxnumber for example 7x4 then 28 cores will be used. 7 for chromosomes and 4 for inside each chromosome.", metavar = "character"),
  "4" = make_option(c("-T", "--tempdir"), type = "character", default = 0, 
                    help = "Optional. Specify the temp directory for storing intermediate files during parallel processing. This will default to the output dir if not explicitly specified.", metavar = "character")
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
tempdir <- input_args$tempdir

# DEBUG ########
# input_GTF_dir <- "/mnt/LTS/projects/2020_isoform_nomenclature/ensembl_GTF_dump/"
# output_dir <- "/mnt/LTS/projects/2020_isoform_nomenclature/ENST_to_HGNC_stable_variant_mapping/"
# ncores <- "8x8"
# tempdir <- 
################

cat("input_GTF_dir:", input_GTF_dir, "\n")
cat("output_dir:", output_dir, "\n")
cat("tempdir:", tempdir, "\n")

# manage parrallellisation rrlllRll

max_number_of_cores <- parallel::detectCores()

if (grepl(x = ncores, pattern = "x") == FALSE) {
  
  if (ncores != 0) {
    number_of_workers <- ncores
    cat(future::availableCores(), "cores will be used\n")
  } else {
    number_of_workers <- future::availableCores()
    cat(future::availableCores(), "cores will be used\n")
  } 
  
} else if (grepl(x = ncores, pattern = "x") == TRUE) {
  
  ncores_level_1 <- ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert
  ncores_level_2 <- ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert
  
  plan(list(tweak(multicore, workers = ncores_level_1, gc = TRUE), 
            tweak(multicore, workers = ncores_level_2, gc = TRUE))
  )
  
  cat((ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert) * (ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert), "cores will be used in total\n")
  cat("first layer:", ncores_level_1, "cores\n")
  cat("second layer:", ncores_level_2, "cores\n")
  
}

options(future.globals.maxSize = 30000000000, future.fork.enable = TRUE)

# BEGIN EXECUTION ###

plan(list(tweak(multicore, workers = ncores_level_1 * ncores_level_2, gc = TRUE), 
          tweak(multicore, workers = 1, gc = TRUE))
)

# import all GTFs
list_all_GTFs <- furrr::future_map(
  .x = list.files(input_GTF_dir, pattern = "gtf", include.dirs = FALSE),
  .f = ~rtracklayer::import(paste(input_GTF_dir, .x, sep = "")) %>% as_tibble,
  .progress = TRUE )

# get vector of pure release numbers
vector_release_numbers <- gsub(x = list.files(input_GTF_dir, pattern = "gtf", include.dirs = FALSE), pattern = "^([^\\.]+)\\.([^\\.]+)\\.(\\d+)\\.(.*)", replacement = "\\3") %>% type.convert

names(list_all_GTFs) <- vector_release_numbers

# get vector of assembly numbers
vector_assembly_numbers <- gsub(x = list.files(input_GTF_dir, pattern = "gtf", include.dirs = FALSE), pattern = "^([^\\.]+)\\.([^\\.]+)\\.(\\d+)\\.(.*)", replacement = "\\2") %>% type.convert

setDTthreads(ncores_level_1 * ncores_level_2)

# add column to each GTF indicating assembly and release version
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

tibble_all_GTFs <- list_all_GTFs %>% rbindlist(fill = TRUE, use.names = TRUE) %>% tibble::as_tibble() %>% readr::type_convert()

# Major analysis point no.1a: Record the historical development of ENST IDs, retirement status of ENST ID as well as the latest ENST ID it occurs in, split by ENST ID
## detect latest release
latest_release_number <- tibble_all_GTFs$ensembl_release_version %>% max

## ENST TRACKING TABLE: filter table for ENST ID and version 
tibble_ENST_ID_version_tracked_by_release <- tibble_all_GTFs[, c("gene_name", "transcript_id", "ensembl_release_version", "genome_assembly")] %>% dplyr::distinct()

## ENST RETIREMENT TABLE: last release + retirement status
tibble_ENST_retirement_status <- tibble_ENST_ID_version_tracked_by_release %>% 
  dplyr::distinct() %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarise("release_last_seen_transcript_id" = max(ensembl_release_version), 
                   "retirement_status_transcript_id" = if (max(ensembl_release_version) == latest_release_number) {"active"} else if (max(ensembl_release_version) < latest_release_number) {"retired"} )

# Major analysis point no.1b: Record the historical development of ENSP IDs, retirement status of ENSP ID as well as the latest ENSP ID it occurs in, split by ENSP ID
## ENSP TRACKING TABLE: filter table for ENSP ID and version 
tibble_ENSP_ID_version_tracked_by_release <- tibble_all_GTFs[, c("gene_name", "protein_id", "ensembl_release_version", "genome_assembly")] %>% dplyr::distinct() %>% .[!is.na(.$protein_id), ]

## ENSP RETIREMENT TABLE: last release + retirement status
tibble_ENSP_retirement_status <- tibble_ENSP_ID_version_tracked_by_release %>% 
  dplyr::distinct() %>%
  dplyr::group_by(protein_id) %>%
  dplyr::summarise("release_last_seen_protein_id" = max(ensembl_release_version), 
                   "retirement_status_protein_id" = if (max(ensembl_release_version) == latest_release_number) {"active"} else if (max(ensembl_release_version) < latest_release_number) {"retired"} )

# Major analysis point no.2: Write down the ENSG tree, as well as the associated gene_name, combined gene_name and the release version
# tibble_historical_ENSG_tree_and_gene_names <- tibble_all_GTFs[, c("gene_id", "gene_name", "ensg_gene_name_combined", "ensembl_release_version", "genome_assembly")] %>% dplyr::distinct()

# Major analysis point no.3a: create HGNC stable IDs
## according to our HGNC gene symbol centric view of the database, each unique HGNC gene symbol constitutes its own family
## gene symbol families are associated with ENSG, ENST and ENSP IDs
## over time, ENSG, ENST and ENSP IDs appear to come and go from the perspective of the gene symbol
## thus, we assign a number to each unique ENSG, ENST and ENSP ID that has ever been assigned to the gene symbol, according to the chronological order in which it made its first ever appearance
## If multiple IDs have been associated with a gene symbol simultaneously in any given release, we distinguish by the ascending order of ID.
list_tibble_all_ENST_IDs_split_by_gene_name <- tibble_all_GTFs %>% 
  .[!is.na(.$ensg_gene_name_combined) & !is.na(.$transcript_id), ] %>%
  dplyr::distinct(gene_id, transcript_id, ensg_gene_name_combined, transcript_version, ensembl_release_version, .keep_all = FALSE) %>% 
  dplyr::group_split(ensg_gene_name_combined) %>%
  set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$ensg_gene_name_combined %>% unique) %>% unlist)

# create HGNC stable transcript IDs
round_robin_pmap_callr(
  .l = list(
    "a1" = purrr::map(.x = parallel::splitIndices(nx = length(list_tibble_all_ENST_IDs_split_by_gene_name), ncl = ceiling(length(list_tibble_all_ENST_IDs_split_by_gene_name)/1000)), .f = ~list_tibble_all_ENST_IDs_split_by_gene_name[.x]),
    "a2" = 1:(purrr::map(.x = parallel::splitIndices(nx = length(list_tibble_all_ENST_IDs_split_by_gene_name), ncl = ceiling(length(list_tibble_all_ENST_IDs_split_by_gene_name)/1000)), .f = ~list_tibble_all_ENST_IDs_split_by_gene_name[.x]) %>% length)
  ),
  .num_workers = ncores_level_1 * ncores_level_2, 
  .env_flag = "user",
  .objects = character(0),
  .temp_path = paste(tempdir, "/tibble_HGNC_stable_transcript_IDs_env.Rdata", sep = ""),
  .status_messages_dir = tempdir, 
  .job_name = "tibble_HGNC_stable_transcript_IDs",
  .f = function(a1, a2) {
    
    # DEBUG ###
    # a1 <- purrr::map(.x = parallel::splitIndices(nx = length(list_tibble_all_ENST_IDs_split_by_gene_name), ncl = ceiling(length(list_tibble_all_ENST_IDs_split_by_gene_name)/500)), .f = ~list_tibble_all_ENST_IDs_split_by_gene_name[.x]) %>% .[[32]]
    # a1 <- list_tibble_all_ENST_IDs_split_by_gene_name$TBCE
    # a1 <- list_tibble_all_ENST_IDs_split_by_gene_name$U2
    # a1 <- list_tibble_all_ENST_IDs_split_by_gene_name$RUNX2
    ###########
    
    cat(a2, "\n")
    
    # CALLR ###
    library(tidyverse)
    library(furrr)
    library(rtracklayer)
    library(data.table)
    library(optparse)
    library(gtools)
    
    library(tictoc)
    ###########
    
    L1_result <- purrr::map(
      .x = a1,
      .f = function(b1) { 
        
        # DEBUG ###
        # b1 <- a1[[1]]
        ###########
        
        # fetch the earliest ensembl release where each ENSG ID joined the gene_name umbrella
        # order by earliest release for each ENSG ID
        tibble_ENSG_earliest_release_joined <- b1 %>%
          dplyr::group_by(gene_id) %>%
          dplyr::summarise("gene_id_earliest_release_joined" = min(ensembl_release_version)) %>%
          dplyr::arrange(gene_id) %>%
          dplyr::arrange(gene_id_earliest_release_joined) %>%
          tibble::add_column("ENSG_historical_join_order" = 1:nrow(.))
        
        tibble_ENST_earliest_release_joined <- b1 %>%
          dplyr::group_by(transcript_id) %>%
          dplyr::summarise("transcript_id_earliest_release_joined" = min(ensembl_release_version)) %>%
          dplyr::arrange(transcript_id) %>%
          dplyr::arrange(transcript_id_earliest_release_joined) %>%
          tibble::add_column("ENST_historical_join_order" = 1:nrow(.))
        
        # table join
        tibble_with_join_orders <- dplyr::left_join(b1, tibble_ENSG_earliest_release_joined) %>%
          dplyr::left_join(., tibble_ENST_earliest_release_joined)
        
        # fix NA transcript_version
        tibble_with_join_orders[is.na(tibble_with_join_orders$transcript_version), "transcript_version"] <- 1
        
        # finally add in the stable transcript number
        tibble_hgnc_stable_transcript_id <- tibble_with_join_orders %>%
          dplyr::mutate("hgnc_stable_transcript_ID" = paste(ensg_gene_name_combined, "-", ENSG_historical_join_order, "enst", ENST_historical_join_order, ".", transcript_version, sep = "") %>% gsub(pattern = "\\-1enst", replacement = "-enst"))
        
        return(tibble_hgnc_stable_transcript_id)
        
        } )
    
    saveRDS(object = L1_result, file = paste(tempdir, "/tibble_HGNC_stable_transcript_IDs_L1_", a2, ".rds", sep = ""), compress = FALSE)
    
    return(NULL)
    
  } ) 

tibble_HGNC_stable_transcript_IDs <- purrr::reduce2(
  .x = purrr::splice(list(tibble()), list.files(path = tempdir, pattern = "tibble_HGNC_stable_transcript_IDs_L1_.*\\.rds", full.names = TRUE) %>% as.list),
  .y = 1:length(list.files(path = tempdir, pattern = "tibble_HGNC_stable_transcript_IDs_L1_.*\\.rds", full.names = TRUE)),
  .f = function(a1, a2, a3) {
    
    # DEBUG ###
    # a1 <- tibble()
    # a2 <- list.files(path = tempdir, pattern = "tibble_HGNC_stable_transcript_IDs_L1_.*\\.rds", full.names = TRUE) %>% .[[1]]
    ###########
    
    message(a3)
    
    tibble_right <- readRDS(file = a2) %>% data.table::rbindlist() %>% tibble::as_tibble()
    
    tibble_rbind <- dplyr::bind_rows(a1, tibble_right)
    
    return(tibble_rbind)
    
  } )

# remove intermediate parallel files
file.remove(list.files(path = tempdir, pattern = "tibble_HGNC_stable_transcript_IDs.*", full.names = TRUE))

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

# create HGNC stable protein IDs
list_tibble_all_ENSP_IDs_split_by_gene_name <- tibble_all_GTFs %>% 
  .[!is.na(.$ensg_gene_name_combined) & !is.na(.$protein_id), ] %>%
  dplyr::distinct(gene_id, protein_id, ensg_gene_name_combined, protein_version, ensembl_release_version, .keep_all = FALSE) %>% 
  dplyr::group_split(ensg_gene_name_combined) %>%
  set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$ensg_gene_name_combined %>% unique) %>% unlist)

round_robin_pmap_callr(
  .l = list(
    "a1" = purrr::map(.x = parallel::splitIndices(nx = length(list_tibble_all_ENSP_IDs_split_by_gene_name), ncl = ceiling(length(list_tibble_all_ENSP_IDs_split_by_gene_name)/1000)), .f = ~list_tibble_all_ENSP_IDs_split_by_gene_name[.x]),
    "a2" = 1:(purrr::map(.x = parallel::splitIndices(nx = length(list_tibble_all_ENSP_IDs_split_by_gene_name), ncl = ceiling(length(list_tibble_all_ENSP_IDs_split_by_gene_name)/1000)), .f = ~list_tibble_all_ENSP_IDs_split_by_gene_name[.x]) %>% length)
  ),
  .num_workers = ncores_level_1 * ncores_level_2, 
  .env_flag = "user",
  .objects = character(0),
  .temp_path = paste(tempdir, "/tibble_HGNC_stable_protein_IDs_env.Rdata", sep = ""),
  .status_messages_dir = tempdir, 
  .job_name = "tibble_HGNC_stable_protein_IDs",
  .f = function(a1, a2) {
    
    # DEBUG ###
    # a1 <- purrr::map(.x = parallel::splitIndices(nx = length(list_tibble_all_ENSP_IDs_split_by_gene_name), ncl = ceiling(length(list_tibble_all_ENSP_IDs_split_by_gene_name)/500)), .f = ~list_tibble_all_ENSP_IDs_split_by_gene_name[.x]) %>% .[[32]]
    # a1 <- list_tibble_all_ENSP_IDs_split_by_gene_name$TBCE
    # a1 <- list_tibble_all_ENSP_IDs_split_by_gene_name$U2
    # a1 <- list_tibble_all_ENSP_IDs_split_by_gene_name$RUNX2
    ###########
    
    cat(a2, "\n")
    
    # CALLR ###
    library(tidyverse)
    library(furrr)
    library(rtracklayer)
    library(data.table)
    library(optparse)
    library(gtools)
    
    library(tictoc)
    ###########
    
    L1_result <- purrr::map(
      .x = a1,
      .f = function(b1) { 
        
        # DEBUG ###
        # b1 <- a1[[1]]
        ###########
        
        # fetch the earliest ensembl release where each ENSG ID joined the gene_name umbrella
        # order by earliest release for each ENSG ID
        tibble_ENSG_earliest_release_joined <- b1 %>%
          dplyr::group_by(gene_id) %>%
          dplyr::summarise("gene_id_earliest_release_joined" = min(ensembl_release_version)) %>%
          dplyr::arrange(gene_id) %>%
          dplyr::arrange(gene_id_earliest_release_joined) %>%
          tibble::add_column("ENSG_historical_join_order" = 1:nrow(.))
        
        tibble_ENSP_earliest_release_joined <- b1 %>%
          dplyr::group_by(protein_id) %>%
          dplyr::summarise("protein_id_earliest_release_joined" = min(ensembl_release_version)) %>%
          dplyr::arrange(protein_id) %>%
          dplyr::arrange(protein_id_earliest_release_joined) %>%
          tibble::add_column("ENSP_historical_join_order" = 1:nrow(.))
        
        # table join
        tibble_with_join_orders <- dplyr::left_join(b1, tibble_ENSG_earliest_release_joined) %>%
          dplyr::left_join(., tibble_ENSP_earliest_release_joined)
        
        # fix NA protein_version
        tibble_with_join_orders[is.na(tibble_with_join_orders$protein_version), "protein_version"] <- 1
        
        # finally add in the stable protein number
        tibble_hgnc_stable_protein_id <- tibble_with_join_orders %>%
          dplyr::mutate("hgnc_stable_protein_ID" = paste(ensg_gene_name_combined, "-", ENSG_historical_join_order, "enst", ENSP_historical_join_order, ".", protein_version, sep = "") %>% gsub(pattern = "\\-1ensp", replacement = "-ensp"))
        
        return(tibble_hgnc_stable_protein_id)
        
      } )
    
    saveRDS(object = L1_result, file = paste(tempdir, "/tibble_HGNC_stable_protein_IDs_L1_", a2, ".rds", sep = ""), compress = FALSE)
    
    return(NULL)
    
  } ) 

tibble_HGNC_stable_protein_IDs <- purrr::reduce2(
  .x = purrr::splice(list(tibble()), list.files(path = tempdir, pattern = "tibble_HGNC_stable_protein_IDs_L1_.*\\.rds", full.names = TRUE) %>% as.list),
  .y = 1:length(list.files(path = tempdir, pattern = "tibble_HGNC_stable_protein_IDs_L1_.*\\.rds", full.names = TRUE)),
  .f = function(a1, a2, a3) {
    
    # DEBUG ###
    # a1 <- tibble()
    # a2 <- list.files(path = tempdir, pattern = "tibble_HGNC_stable_protein_IDs_L1_.*\\.rds", full.names = TRUE) %>% .[[1]]
    ###########
    
    message(a3)
    
    tibble_right <- readRDS(file = a2) %>% data.table::rbindlist() %>% tibble::as_tibble()
    
    tibble_rbind <- dplyr::bind_rows(a1, tibble_right)
    
    return(tibble_rbind)
    
  } )

# remove intermediate parallel files
file.remove(list.files(path = tempdir, pattern = "tibble_HGNC_stable_protein_IDs.*", full.names = TRUE))

# get rid of NA protein_versions
tibble_HGNC_stable_protein_IDs <- tibble_HGNC_stable_protein_IDs %>% 
  dplyr::mutate("hgnc_stable_protein_ID" = gsub(x = `hgnc_stable_protein_ID`, pattern = ".NA$", replacement = ""))

# TOTAL MAPPING TABLE
## we are going to combine the information from all 1 + 1 + 2 tables into one
tibble_total_mapping_table_transcripts <- dplyr::left_join(tibble_HGNC_stable_transcript_IDs, tibble_ENST_retirement_status)
tibble_total_mapping_table_proteins <- dplyr::left_join(tibble_HGNC_stable_protein_IDs, tibble_ENSP_retirement_status)

if (nrow(tibble_total_mapping_table_transcripts) != nrow(tibble_HGNC_stable_transcript_IDs)) {
  stop("imperfect join. check whether there may be some ENST IDs that have been assigned to more than one ENSG ID in the past.")
}

if (nrow(tibble_total_mapping_table_proteins) != nrow(tibble_HGNC_stable_protein_IDs)) {
  stop("imperfect join. check whether there may be some ENSP IDs that have been assigned to more than one ENSG ID in the past.")
}

# ANNOTATE ALL RELEASES
tibble_all_annotated_releases <- dplyr::left_join(tibble_all_GTFs, tibble_total_mapping_table_transcripts) %>%
  dplyr::left_join(., tibble_total_mapping_table_proteins)

# list-ify and save
list_all_annotated_releases <- tibble_all_annotated_releases %>%
  dplyr::group_split(ensembl_release_version) %>%
  set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$ensembl_release_version %>% unique) %>% unlist)

plan(list(tweak(multicore, workers = ncores_level_1),
          tweak(multicore, workers = ncores_level_2)))

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

plan(list(tweak(multicore, workers = ncores_level_1 * ncores_level_2, gc = TRUE), 
          tweak(multicore, workers = 1, gc = TRUE))
)

# save ALL annotated releases
furrr::future_map2(
  .x = list_all_annotated_releases_fixed,
  .y = names(list_all_annotated_releases_fixed),
  .f = ~write.table(x = .x, file = paste(output_dir, "annotated_ensembl_gtf_release_", .y, ".txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE), 
  .progress = TRUE )

# 1. save ENST retirement tracking
write.table(x = tibble_ENST_retirement_status, file = paste(output_dir, "ENST_retirement_status.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(x = tibble_ENSP_retirement_status, file = paste(output_dir, "ENSP_retirement_status.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 2. save gene tree tracking
# write.table(x = tibble_historical_ENSG_tree_and_gene_names, file = paste(output_dir, "latest_gene_tree_tracking.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 3. save HGNC stable transcript mapping table
write.table(x = tibble_total_mapping_table_transcripts, file = paste(output_dir, "total_mapping_table_transcripts.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(x = tibble_total_mapping_table_proteins, file = paste(output_dir, "total_mapping_table_proteins.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# add information about latest update
cat(paste("Date of last update: ", date(), "\nLatest Ensembl release version: ", latest_release_number, "\n\nSystem info of machine that retrieved last update:\n\n", sep = ""), file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""))

write.table(Sys.info() %>% tibble::enframe(name = "id", value = "value"), file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""), sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

write("\n", file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""), append = TRUE)

write.table(R.version %>% unlist(use.names = TRUE) %>% tibble::enframe(name = "id", value = "value"), file = paste(output_dir, "latest_mapping_metadata.txt", sep = ""), sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

# finish counting
tictoc::toc()

q()

