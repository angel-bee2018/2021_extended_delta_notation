# The EDN Suite, a collection of 3 tools to aid the usage of Extended Delta Notation
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

user_instructions <- c("

- Make sure ALL the exons are contained within the VSR!

                       ")

# DEFINE ENVIRONMENT
library(shiny)
library(shinyWidgets)

library(tidyverse)
library(ggh4x)

library(rtracklayer)
library(gtools)
library(data.table)
data.table::setDTthreads(12)

library(furrr)

library(gridSVG)

options(future.globals.maxSize = 30000000000, future.fork.enable = TRUE)

options(shiny.maxRequestSize = 1500*1024^2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# tibble_ref_gtf <- rtracklayer::import(con = "data/latest_ensembl_GTF.gtf", format = "gtf") %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)

# write.table(tibble_ref_gtf, file = "data/latest_ensembl_GTF.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# tibble_ref_gtf_hg18 <- vroom::vroom(file = "data/annotated_ensembl_gtf_release_54.txt", delim = "\t") %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
# tibble_ref_gtf_hg19 <- vroom::vroom(file = "data/annotated_ensembl_gtf_release_75.txt", delim = "\t") %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
# tibble_ref_gtf_hg38 <- vroom::vroom(file = "data/annotated_ensembl_gtf_release_102.txt", delim = "\t") %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)

# DEFINE FUNCTIONS ####

## FUNCTION TO EXTRACT REFERENCE EXONS WHICH OVERLAP EXACTLY WITH QUERY EXONS
## NOTE: to be used with purrr
## input: spliceregion_list: a list containing details of ONE junction: $chr, $diff_exon_start, $diff_exon_end
## tibble_gtf_table: rtracklayer::import(...) %>% as_tibble. can be any GTF. reconstructed or reference.
## index: loop progress marker to be used with imap

extract_matching.exons <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 0, right_tolerance = 0, return_type = "exon") {
    
    # DEBUG ###################
    # index <- 1
    # spliceregion_list <- wide_tibble_of_all_unique_VSR_and_exon_coords_array.tree_not_IR[[index]]
    # # tibble_gtf_table <- tibble_ref_gtf
    # tibble_gtf_table <- tibble_recon_gtf
    # stranded = FALSE
    ###########################
    
    # print(query_chr, "\n")
    # print(query_start, "\n")
    # print(query_end, "\n")
    # print(query_strand, "\n")
    # print(tibble_gtf_table, "\n")
    
    # print(paste("now processing junction number", index))
    
    # left_query_shift <<- left_query_shift
    # right_query_shift <<- right_query_shift
    # left_tolerance <<- left_tolerance
    # right_tolerance <<- right_tolerance
    
    if (!(query_strand == "+" | query_strand == "-")) {
        
        # cat("a\n")
        
        # +/- 1 nt tolerance
        tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% 
            .[which(.$start > ((query_start %>% as.numeric) - 1 + left_query_shift - left_tolerance) & .$end < ((query_end %>% as.numeric) + 1 + right_query_shift + right_tolerance)), ] %>% 
            .[which((.$start < ((query_start %>% as.numeric) + 1 + left_query_shift + left_tolerance) & .$end > ((query_end %>% as.numeric) - 1 + right_query_shift - right_tolerance))), ] %>% 
            .[which(.$type %in% return_type), ]
        
    } else if (query_strand == "+" | query_strand == "-") {
        
        # cat("b\n")
        
        # +/- 1 nt tolerance
        tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws &
                                                                 tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
            .[which(.$start > ((query_start %>% as.numeric) - 1 + left_query_shift - left_tolerance) & .$end < ((query_end %>% as.numeric) + 1 + right_query_shift + right_tolerance)), ] %>% 
            .[which((.$start < ((query_start %>% as.numeric) + 1 + left_query_shift + left_tolerance) & .$end > ((query_end %>% as.numeric) - 1 + right_query_shift - right_tolerance))), ] %>% 
            .[which(.$type %in% return_type), ]
        
    }
    
    return(tibble_gtf_subset_matching_exons)
    
}

## END extract_matching.exons() ###

## FUNCTION TO EXTRACT TRANSCRIPTS WITH JUNCTION-FLANKING EXONS.
## NOTE: to be used with purrr
## details of ONE junction: $chr, $start, $end, $strand
## tibble_gtf_table: rtracklayer::import(...) %>% as_tibble. can be any GTF. reconstructed or reference.
## index: loop progress marker to be used with imap

extract_junction.flanking.exons <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 0, right_tolerance = 0, match_consecutive = TRUE, match_inside_same_transcript = TRUE, return_type = "exon") {
    
    # DEBUG ###################
    
    # query_chr = query_chr
    # query_start = a1$event_region_start %>% type.convert
    # query_end = a1$event_region_end %>% type.convert
    # query_strand = "*"
    # tibble_gtf_table = tibble_ref_gtf
    # tolerance_left = 0
    # tolerance_right = 0
    # tolerance_inside = 0
    # tolerance_outside = 0
    # match_consecutive = FALSE
    # return_type = "exon"
    
    ###########################
    
    # print(paste("now processing junction number", index))
    
    # left_query_shift <<- left_query_shift
    # right_query_shift <<- right_query_shift
    # left_tolerance <<- left_tolerance
    # right_tolerance <<- right_tolerance
    
    if (query_strand == "." | query_strand == 0 | query_strand == "*") {
        
        tibble_gtf_subset_flanking_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% .[.$start <= ((query_end %>% as.numeric) + 1 + right_query_shift + right_tolerance) & .$end >= ((query_start %>% as.numeric) - 1 + left_query_shift - left_tolerance), ] %>% .[!(.$start <= ((query_end %>% as.numeric) + right_query_shift - right_tolerance) & .$end >= ((query_start %>% as.numeric) + left_query_shift + left_tolerance)), ] %>% .[.$type %in% return_type, ]
        
    } else if (query_strand == "+" | query_strand == "-") {
        
        tibble_gtf_subset_flanking_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% .[.$strand == query_strand %>% trimws, ] %>% .[.$start <= ((query_end %>% as.numeric) + 1 + right_query_shift + right_tolerance) & .$end >= ((query_start %>% as.numeric) - 1 + left_query_shift - left_tolerance), ] %>% .[!(.$start <= ((query_end %>% as.numeric) + right_query_shift - right_tolerance) & .$end >= ((query_start %>% as.numeric) + left_query_shift + left_tolerance)), ] %>% .[.$type %in% return_type, ]
        
    }
    
    list_of_junction_associated_transcripts <- tibble_gtf_subset_flanking_exons$transcript_id %>% unique %>% array_tree %>% flatten
    
    # make a list for each transcript that directly flanks a junction.
    # then filter so that there are only a) exon PAIRS which b) are directly connected in the mature (spliced) transcript
    
    if (match_inside_same_transcript == FALSE) {
        
        list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~tibble_gtf_subset_flanking_exons[tibble_gtf_subset_flanking_exons$transcript_id == .x, ] %>% dplyr::arrange(exon_number %>% as.numeric)) %>% set_names(list_of_junction_associated_transcripts)
        
    } else if (match_inside_same_transcript == TRUE & match_consecutive == TRUE) {
        
        list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~tibble_gtf_subset_flanking_exons[tibble_gtf_subset_flanking_exons$transcript_id == .x, ] %>% dplyr::arrange(exon_number %>% as.numeric)) %>% set_names(list_of_junction_associated_transcripts) %>% keep(.x = ., .p = ~nrow(.x) == 2) %>% keep(.x = ., .p = ~abs((.x[2, "exon_number"] %>% paste %>% as.numeric) - (.x[1, "exon_number"] %>% paste %>% as.numeric)) == 1)
        
    } else if (match_inside_same_transcript == TRUE & match_consecutive == FALSE) {
        
        list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~tibble_gtf_subset_flanking_exons[tibble_gtf_subset_flanking_exons$transcript_id == .x, ] %>% dplyr::arrange(exon_number %>% as.numeric)) %>% set_names(list_of_junction_associated_transcripts) %>% keep(.x = ., .p = ~nrow(.x) == 2)
        
    }
    
    return(list_of_tibbles_flanking_exon_gtf.entries_per_transcript)
    
}

## END extract_junction.flanking.exons_JUM() ###

## FUNCTION TO EXTRACT REFERENCE EXONS WHICH OVERLAP EXACTLY WITH QUERY EXONS
## NOTE: to be used with purrr
## input: spliceregion_list: a list containing details of ONE junction: $chr, $diff_exon_start, $diff_exon_end
## tibble_gtf_table: rtracklayer::import(...) %>% as_tibble. can be any GTF. reconstructed or reference.
## index: loop progress marker to be used with imap

extract_overlapping_features <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 0, right_tolerance = 0, return_type = NULL) {
    
    # DEBUG ###################
    # index <- 1
    # spliceregion_list <- wide_tibble_of_all_unique_VSR_and_exon_coords_array.tree_not_IR[[index]]
    # # tibble_gtf_table <- tibble_ref_gtf
    # tibble_gtf_table <- tibble_recon_gtf
    # stranded = FALSE
    ###########################
    
    # cat(query_chr, "\n")
    # cat(query_start, "\n")
    # cat(query_end, "\n")
    # cat(query_strand, "\n")
    # print(tibble_gtf_table, "\n")
    # cat(left_query_shift, "\n")
    # cat(right_query_shift, "\n")
    # cat(left_tolerance, "\n")
    # cat(right_tolerance, "\n")
    
    # print(paste("now processing junction number", index))
    
    # left_query_shift <<- left_query_shift
    # right_query_shift <<- right_query_shift
    # left_tolerance <<- left_tolerance
    # right_tolerance <<- right_tolerance
    
    if (!(query_strand == "+" | query_strand == "-")) {
        
        # +/- 1 nt tolerance
        tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% 
            .[which(.$end >= ((query_start %>% as.numeric) + left_query_shift - left_tolerance) & .$start <= ((query_end %>% as.numeric) + right_query_shift + right_tolerance)), ]
        
    } else if (query_strand == "+" | query_strand == "-") {
        
        # +/- 1 nt tolerance
        tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws &
                                                                 tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
            .[which(.$end >= ((query_start %>% as.numeric) + left_query_shift - left_tolerance) & .$start <= ((query_end %>% as.numeric) + right_query_shift + right_tolerance)), ]
        
    }
    
    if (is.null(return_type) == FALSE) {
        return(tibble_gtf_subset_matching_exons[tibble_gtf_subset_matching_exons$type %in% return_type, ])
    } else if (is.null(return_type) == TRUE) {
        return(tibble_gtf_subset_matching_exons)
    }
    
    return(tibble_gtf_subset_matching_exons)
    
}

## END extract_overlapping_features() ###

## FUNCTIONS TO MAGNETISE A USER VALUE TO REFERENCE START COORDS
### ref_start/end_shift: use +/-1 to get intronic starts/ends
### match to reference starts
magnetise_genome_position_to_ref_starts <- function(query_chr, query_coord, query_strand, tibble_gtf_table, query_shift = 0, query_tolerance = 0, ref_start_shift = 0, return_type) {
    
    # DEBUG ###################
    # query_chr <- query_chr
    # query_coord <- 153411582
    # query_strand <- query_strand
    # query_shift <- right_query_shift
    # query_tolerance <- right_tolerance
    # ref_start_shift <- -1
    # return_type <- "exon"
    ###########################
    
    if (!(query_strand == "+" | query_strand == "-")) {
        
        # +/- 1 nt tolerance
        tibble_ref_starts_matched_to_query_coord <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% 
            .[which((.$start + ref_start_shift) >= (query_coord + query_shift - query_tolerance) & (.$start + ref_start_shift) <= (query_coord + query_shift + query_tolerance)), ] %>% 
            .[which(.$type %in% return_type), ]
    } else if (query_strand == "+" | query_strand == "-") {
        
        # +/- 1 nt tolerance
        tibble_ref_starts_matched_to_query_coord <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws & tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
            .[which((.$start + ref_start_shift) >= (query_coord + query_shift - query_tolerance) & (.$start + ref_start_shift) <= (query_coord + query_shift + query_tolerance)), ] %>% 
            .[which(.$type %in% return_type), ]
        
    }
    
    magnetised_coord <- (tibble_ref_starts_matched_to_query_coord$start + ref_start_shift) %>%
        .[which(abs(tibble_ref_starts_matched_to_query_coord$start - query_coord) == min(abs(tibble_ref_starts_matched_to_query_coord$start - query_coord)))] %>% 
        .[1]
    
    return(list(
        "tibble_ref_starts_matched_to_query_coord" = tibble_ref_starts_matched_to_query_coord,
        "magnetised_coord" = magnetised_coord
    ) )
    
}

## END magnetise_genome_position_to_ref_starts() ###

### match to reference ends
magnetise_genome_position_to_ref_end <- function(query_chr, query_coord, query_strand, tibble_gtf_table, query_shift = 0, query_tolerance = 0, ref_end_shift = 0, return_type) {
    
    # DEBUG ###################
    # query_chr <- query_chr
    # query_coord <- query_VSR_end
    # query_strand <- query_strand
    # query_shift <- right_query_shift
    # query_tolerance <- right_tolerance
    # ref_end_shift <- -1
    ###########################
    
    if (!(query_strand == "+" | query_strand == "-")) {
        
        # +/- 1 nt tolerance
        tibble_ref_ends_matched_to_query_coord <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% 
            .[which((.$end + ref_end_shift) >= (query_coord + query_shift - query_tolerance) & (.$end + ref_end_shift) <= (query_coord + query_shift + query_tolerance)), ] %>% 
            .[which(.$type %in% return_type), ]
    } else if (query_strand == "+" | query_strand == "-") {
        
        # +/- 1 nt tolerance
        tibble_ref_ends_matched_to_query_coord <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws & tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
            .[which((.$end + ref_end_shift) >= (query_coord + query_shift - query_tolerance) & (.$end + ref_end_shift) <= (query_coord + query_shift + query_tolerance)), ] %>% 
            .[which(.$type %in% return_type), ]
        
    }
    
    magnetised_coord <- (tibble_ref_ends_matched_to_query_coord$end + ref_end_shift) %>%
        .[which(abs(tibble_ref_ends_matched_to_query_coord$end - query_coord) == min(abs(tibble_ref_ends_matched_to_query_coord$end - query_coord)))] %>% 
        .[1]
    
    return(list(
        "tibble_ref_ends_matched_to_query_coord" = tibble_ref_ends_matched_to_query_coord,
        "magnetised_coord" = magnetised_coord
    ) )
    
}

## END magnetise_genome_position_to_ref_starts() ###

## VSRs/LISs: FUNCTION TO FIND THE STABLE VARIANT NUMBER ###
### Behaviour: accepts the start/end coords of the VSR and a tibble of alternative exon start/ends
### find the transcript with the greatest number of vertices in common.
### as PSI-Sigma says the VSR region is the full intron for IR events, we will not use the VSR for IR matching.

VSR_select_reference_transcript_variant <- function(VSR_coordinates, tibble_VSR_exon_start_end, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 0, right_tolerance = 0) {
    
    # DEBUG ###
    
    # VSR_coordinates <- "16:89737976-89740903"
    # 
    # tibble_VSR_exon_start_end <- tribble(
    #     ~start, ~end,
    #     89740100, 89740803
    # )
    
    # matches exactly
    # VSR_coordinates <- "16:2757946-2768996:+"
    # doesn't match
    # VSR_coordinates <- "16:2758130-2768358:+"
    # tibble_VSR_exon_start_end <- tribble(
    # ~start, ~end,
    # boundary 5' extension
    # 2757946, 2757975,
    # IR
    # 2758548, 2758984,
    # SE
    # 2759140, 2759172,
    # internal A3SS
    # 2760350, 2760499,
    # completely novel exon
    # 2760593, 2760958,
    # completely novel exon
    # 2761023, 2761661
    # )
    
    VSR_coordinates <- VSR_coordinates
    tibble_VSR_exon_start_end <- list_tibble_exon_start_end_per_LIS[[1]]
    left_query_shift <- 0
    right_query_shift <- 0
    left_tolerance <- 1
    right_tolerance <- 1
    tibble_gtf_table <- tibble_ref_gtf
    
    ###########
    
    # left_query_shift <- left_query_shift %>% paste %>% type.convert
    # right_query_shift <- right_query_shift %>% paste %>% type.convert
    # left_tolerance <- left_tolerance %>% paste %>% type.convert
    # right_tolerance <- right_tolerance %>% paste %>% type.convert
    # 
    # print(left_query_shift)
    # print(right_query_shift)
    # print(left_tolerance)
    # print(right_tolerance)
    
    # demultiplex the coordinates
    query_chr <- gsub(x = VSR_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
    query_VSR_start <- gsub(x = VSR_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2")
    query_VSR_end <- gsub(x = VSR_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3")
    query_strand <- gsub(x = VSR_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
    
    # EXON TIBBLE
    # sort the exon tibble in order of increasing co-ordinates
    tibble_VSR_exon_start_end <- tibble_VSR_exon_start_end %>% dplyr::arrange(start)
    
    # detect A3SS/A5SS
    ## if there are any alternative exons which have a co-ordinate in common with the VSR, then it's most likely a A3SS/A5SS exon.
    ## however, there is the extremely rare chance that it's actually an exon in its own right.
    ## in that case, we flag it and find whether it has a match in the reference. If it indeed doesn't have a match in the reference (i.e. the provided exon marks an extension), then we will only use the inner vertex.
    tibble_VSR_exon_start_end_flagged_extensions <- tibble_VSR_exon_start_end %>% 
        dplyr::mutate("left_end_of_VSR" = tibble_VSR_exon_start_end$start == query_VSR_start,
                      "right_end_of_VSR" = tibble_VSR_exon_start_end$end == query_VSR_end)
    
    # detect novel/known exons
    ## match exons with the reference to determine which are novel.
    ### preallocate
    tibble_VSR_exon_start_end_flagged_extensions <- tibble_VSR_exon_start_end_flagged_extensions %>% add_column("exon_matches_to_any_reference" = purrr::map2(.x = tibble_VSR_exon_start_end_flagged_extensions$start, .y = tibble_VSR_exon_start_end_flagged_extensions$end, .f = ~extract_matching.exons(query_chr = query_chr, query_start = .x, query_end = .y, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon") %>% nrow > 0) %>% unlist)
    
    # detect intron retention
    ## we will match every exon in the tibble to check for whether they are junction-spanning.
    tibble_VSR_exon_start_end_flagged_IR <- tibble_VSR_exon_start_end_flagged_extensions %>%
        dplyr::mutate("IR" = purrr::map2(.x = tibble_VSR_exon_start_end_flagged_extensions$start, .y = tibble_VSR_exon_start_end_flagged_extensions$end, .f = ~extract_junction.flanking.exons(query_chr = query_chr, query_start = .x, query_end = .y, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, match_consecutive = FALSE, return_type = "exon") %>% length > 0) %>% unlist)
    
    # generate A3SS/A5SS flags
    flag_boundary_A5SS_detected <- any(tibble_VSR_exon_start_end_flagged_IR$left_end_of_VSR == TRUE & tibble_VSR_exon_start_end_flagged_IR$exon_matches_to_any_reference == FALSE & tibble_VSR_exon_start_end_flagged_IR$IR == FALSE)
    flag_boundary_A3SS_detected <- any(tibble_VSR_exon_start_end_flagged_IR$right_end_of_VSR == TRUE & tibble_VSR_exon_start_end_flagged_IR$exon_matches_to_any_reference == FALSE & tibble_VSR_exon_start_end_flagged_IR$IR == FALSE)
    # also generate a flag indicating whether all exons are IR.
    flag_is_all_IR <- all(tibble_VSR_exon_start_end_flagged_IR$IR == TRUE)
    
    # MATCH VSR TO THE REFERENCE GTF
    ## the final transcript variant used must absolutely contain the VSR.
    ## however, if we had a bona-fide A3/A5SS extension at the boundary, then we have to replace the covered end with the exon interior.
    
    # this is priority. matched transcripts from VSR comes first.
    VSR_selected_transcripts <- NULL
    
    effective_VSR_start <- query_VSR_start
    effective_VSR_end <- query_VSR_end
    
    ### modify the effective VSR ends if there were boundary A3/5SS exon extensions.
    if (flag_boundary_A5SS_detected == TRUE) {
        
        effective_VSR_start <- tibble_VSR_exon_start_end_flagged_IR[tibble_VSR_exon_start_end_flagged_IR$left_end_of_VSR == TRUE & tibble_VSR_exon_start_end_flagged_IR$exon_matches_to_any_reference == FALSE & tibble_VSR_exon_start_end_flagged_IR$IR == FALSE, "end"] %>% unlist + 1
        
    }
    
    if (flag_boundary_A3SS_detected == TRUE) {
        
        effective_VSR_end <- tibble_VSR_exon_start_end_flagged_IR[tibble_VSR_exon_start_end_flagged_IR$right_end_of_VSR == TRUE & tibble_VSR_exon_start_end_flagged_IR$exon_matches_to_any_reference == FALSE & tibble_VSR_exon_start_end_flagged_IR$IR == FALSE, "start"] %>% unlist - 1
        
    }
    
    # match effective VSR to the reference GTF 
    list_VSR_matched_transcripts <- extract_junction.flanking.exons(query_chr = query_chr, query_start = effective_VSR_start, query_end = effective_VSR_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, match_consecutive = FALSE, return_type = "exon")
    
    ## test for transcript overlap
    tibble_trancripts_overlapping_VSR <- extract_overlapping_features(query_chr = query_chr, query_start = query_VSR_start, query_end = query_VSR_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "transcript")
    
    # find vertices that overlap with the VSR ends
    ## LEFT VSR END
    tibble_reference_entries_with_exon_ends_overlapping_VSR_left <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr, ] %>% .[.$type == "exon", ] %>% .[.$end <= (effective_VSR_start %>% type.convert - 1 + left_query_shift + left_tolerance) & .$end >= (effective_VSR_start %>% type.convert - 1 + left_query_shift - left_tolerance), ] 
    effective_magnetised_VSR_start <- tibble_reference_entries_with_exon_ends_overlapping_VSR_left %>% .$end %>% unique + 1
    ## RIGHT VSR END
    tibble_reference_entries_with_exon_starts_overlapping_VSR_right <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr, ] %>% .[.$type == "exon", ] %>% .[.$start <= (effective_VSR_end %>% type.convert + 1 + right_query_shift + right_tolerance) & .$start >= (effective_VSR_end %>% type.convert + 1 + right_query_shift - right_tolerance), ]
    effective_magnetised_VSR_end <- tibble_reference_entries_with_exon_starts_overlapping_VSR_right %>% .$start %>% unique - 1
    
    # VSR MATCH TO ONLY ONE TRANSCRIPT
    ## if there is only one possible VSR match, then we can only use that transcript.
    if (list_VSR_matched_transcripts %>% length == 1) {
        
        VSR_selected_transcripts <- names(list_VSR_matched_transcripts) %>% paste
        
        flag_is_intergenic <- FALSE
        flag_VSR_is_exactly_matched <- TRUE
        
        # VSR MATCH TO MORE THAN ONE TRANSCRIPT
    } else if (list_VSR_matched_transcripts %>% length > 1) {
        
        flag_is_intergenic <- FALSE
        flag_VSR_is_exactly_matched <- TRUE
        
        # VSR NO MATCH TO ANY TRANSCRIPTS
        # condition to account for VSRs which didn't match to any transcript
    } else if (list_VSR_matched_transcripts %>% length == 0) {
        
        flag_VSR_is_exactly_matched <- FALSE
        
        # test for intergenic VSR
        ## if not intergenic, then existing exons can be modified.
        
        flag_is_intergenic <- tibble_trancripts_overlapping_VSR %>% nrow == 0
        
    }
    
    # EXON MATCHING
    # if not intergenic, then we can proceed as normal.
    if (flag_is_intergenic == FALSE) {
        
        # consolidate vertices
        ## these consist of exon start and end vertices
        ## bona fide A3/5SS extensions at the ends (common end with VSR + not boundary skipped exon) will have only the interior vertex kept.
        ## IR exons will have their start/ends swapped and expanded 
        vector_consolidated_exon_starts <- c(
            tibble_VSR_exon_start_end_flagged_IR[tibble_VSR_exon_start_end_flagged_IR$left_end_of_VSR == FALSE & tibble_VSR_exon_start_end_flagged_IR$IR == FALSE, "start"] %>% unlist, 
            tibble_VSR_exon_start_end_flagged_IR[tibble_VSR_exon_start_end_flagged_IR$left_end_of_VSR == TRUE & tibble_VSR_exon_start_end_flagged_IR$exon_matches_to_any_reference == TRUE & tibble_VSR_exon_start_end_flagged_IR$IR == FALSE, "start"] %>% unlist, 
            tibble_VSR_exon_start_end_flagged_IR[tibble_VSR_exon_start_end_flagged_IR$IR == TRUE, "end"] + 1 %>% unlist)
        
        vector_consolidated_exon_ends <- c(
            tibble_VSR_exon_start_end_flagged_IR[tibble_VSR_exon_start_end_flagged_IR$right_end_of_VSR == FALSE & tibble_VSR_exon_start_end_flagged_IR$IR == FALSE, "end"] %>% unlist, 
            tibble_VSR_exon_start_end_flagged_IR[tibble_VSR_exon_start_end_flagged_IR$right_end_of_VSR == TRUE & tibble_VSR_exon_start_end_flagged_IR$exon_matches_to_any_reference == TRUE & tibble_VSR_exon_start_end_flagged_IR$IR == FALSE, "end"] %>% unlist, 
            tibble_VSR_exon_start_end_flagged_IR[tibble_VSR_exon_start_end_flagged_IR$IR == TRUE, "start"] - 1 %>% unlist)
        
        # if all events are IR OR if the VSR didn't match to anything, we can't base everything on the VSR anymore.
        if (flag_is_all_IR == FALSE & flag_VSR_is_exactly_matched == TRUE) {
            
            # subset reference tibble for only exonic entries which contain the consolidated vertices
            if (query_strand %in% c("+", "-")) {
                
                tibble_ref_gtf_exons_with_common_vertices <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr, ] %>% .[.$strand == query_strand, ] %>% .[.$start %in% vector_consolidated_exon_starts | .$end %in% vector_consolidated_exon_ends, ] %>% .[.$type == "exon" & .$transcript_id %in% names(list_VSR_matched_transcripts), ]
                
            } else if (!query_strand %in% c("+", "-")) {
                
                tibble_ref_gtf_exons_with_common_vertices <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr, ] %>% .[.$start %in% vector_consolidated_exon_starts | .$end %in% vector_consolidated_exon_ends, ] %>% .[.$type == "exon" & .$transcript_id %in% names(list_VSR_matched_transcripts), ]
                
            }
            
        } else if (flag_is_all_IR == TRUE | flag_VSR_is_exactly_matched == FALSE) {
            
            # subset reference tibble for only exonic entries which contain the consolidated vertices
            if (query_strand %in% c("+", "-")) {
                
                tibble_ref_gtf_exons_with_common_vertices <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr, ] %>% .[.$strand == query_strand, ] %>% .[.$start %in% vector_consolidated_exon_starts | .$end %in% vector_consolidated_exon_ends, ] %>% .[.$type == "exon", ]
                
            } else if (!query_strand %in% c("+", "-")) {
                
                tibble_ref_gtf_exons_with_common_vertices <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr, ] %>% .[.$start %in% vector_consolidated_exon_starts | .$end %in% vector_consolidated_exon_ends, ] %>% .[.$type == "exon", ]
                
            }
            
        }
        
        # ACCOUNT FOR THE POSSIBILITY OF NO INTERIOR EXONS IN COMMON
        ## IF: 
        ## - exons are completely novel BUT the VSR is matched: take the lowest matched transcript number. The internal exons will be defined in terms of the VSR boundary. It will look like (SRRM2-v1 E2) E(2344_2389) E(545986_546433) (SRRM2-v1 E3) OR (SRRM-v1 E2) J(2343nt) E(46nt) etc...
        ## - exons are completely novel AND the VSR is not matched: take the transcript which has exons closest to the VSR boundaries. The VsR boundary will be defined in terms of this. The internal exons will be defined in terms of the VSR boundary.
        ## these will be dealt with later.
        if (tibble_ref_gtf_exons_with_common_vertices %>% nrow == 0) {
            
            flag_no_common_vertices <- TRUE
            
            if (flag_VSR_is_exactly_matched == TRUE) {
                
                selected_transcript_id <- mixedsort(tibble_trancripts_overlapping_VSR$transcript_id) %>% .[1]
                
            } else if (flag_VSR_is_exactly_matched == FALSE) {
                
                # find the transcript with exon end closest to the VSR start, and exon start closest to the VSR end.
                ## retrieve exon-only entries
                tibble_exonic_entries_transcripts_overlapping_VSR <- tibble_gtf_table[tibble_gtf_table$transcript_id %in% (tibble_trancripts_overlapping_VSR$transcript_id %>% unique) &
                                                                                          tibble_gtf_table$type == "exon", ]
                ## list-ify by transcript_id
                list_exonic_entries_transcripts_overlapping_VSR_by_transcript_id <- tibble_exonic_entries_transcripts_overlapping_VSR %>% 
                    dplyr::group_split(transcript_id) %>% 
                    set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$transcript_id %>% unique) %>% unlist)
                
                ## for each transcript_id, for each row, calculate the vertex differences, hence the sum of minimum absolute differences
                list_rowwise_vertex_differences_by_transcript_id <- purrr::map(
                    .x = list_exonic_entries_transcripts_overlapping_VSR_by_transcript_id,
                    .f = function(a1) {
                        
                        a1 %>% 
                            dplyr::mutate("VSR_start_minus_exon_ends" = (query_VSR_start %>% type.convert) - a1$end,
                                          "exon_starts_minus_VSR_end" = a1$start - (query_VSR_end %>% type.convert)) %>%
                            return
                        
                    } )
                
                vector_sum_of_vertex_differences_by_transcript_id <- purrr::map(
                    .x = list_rowwise_vertex_differences_by_transcript_id,
                    .f = ~min(abs(.x$VSR_start_minus_exon_ends)) + min(abs(.x$exon_starts_minus_VSR_end))) %>% unlist
                
                # retrieve the items with the lowest sum.
                vector_minimum_sum <- vector_sum_of_vertex_differences_by_transcript_id[which(vector_sum_of_vertex_differences_by_transcript_id == min(vector_sum_of_vertex_differences_by_transcript_id))]
                # extract only the lowest transcript_id out of these
                VSR_selected_transcripts <- names(vector_minimum_sum) %>% mixedsort %>% .[1]
                
                # # find the exons overlapping with the VSR, choose the transcript with the least number of exons.
                # ## find reference exons in the VSR
                # tibble_exons_overlapping_VSR <- extract_overlapping_features(query_chr = query_chr, query_start = query_VSR_start, query_end = query_VSR_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1, return_type = "exon")
                # ## create a tally of number of exons for each transcript
                # tibble_tally_number_of_exons_in_VSR_by_transcript_id <- tibble_exons_overlapping_VSR %>% 
                #     dplyr::group_by(transcript_id) %>%
                #     dplyr::summarise("tally" = n()) %>% 
                #     dplyr::arrange(tally)
                # # get selected transcript id
                # selected_transcript_id <- tibble_tally_number_of_exons_in_VSR_by_transcript_id[1, "transcript_id"] %>% paste
                
            }
            
        } else if (tibble_ref_gtf_exons_with_common_vertices %>% nrow > 0) {
            
            flag_no_common_vertices <- FALSE
            
            # list-ify by transcript_id and test for the number of vertices in common
            ## group_split
            list_ref_gtf_exons_with_common_vertices_by_transcript_id <- tibble_ref_gtf_exons_with_common_vertices %>% 
                dplyr::group_split(transcript_id) %>% 
                set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$transcript_id %>% unique) %>% unlist)
            ## test for number of vertices in common
            list_number_of_reference_vertices_in_common <- purrr::map(
                .x = list_ref_gtf_exons_with_common_vertices_by_transcript_id,
                .f = ~intersect(c(.x$start, .x$end) %>% unique, c(vector_consolidated_exon_starts, vector_consolidated_exon_ends) %>% unique) %>% length
            )
            ## percolate and tibblise
            tibble_common_vertices_count <- list_number_of_reference_vertices_in_common %>% as_tibble(.name_repair = "unique") %>% t %>% as_tibble(rownames = "transcript_id") %>% setNames(c("transcript_id", "number_vertices_in_common"))
            ## extract the lowest transcript_id with the highest common vertices
            selected_transcript_id <- tibble_common_vertices_count[tibble_common_vertices_count$number_vertices_in_common == max(tibble_common_vertices_count$number_vertices_in_common), "transcript_id"] %>% unlist %>% mixedsort %>% .[1]
            
        }
        
        # if VSR authority exists, then override body exon matches.
        if (VSR_selected_transcripts %>% is.null == FALSE) {
            
            selected_transcript_id <- VSR_selected_transcripts
            
        }
        
        ## extract all the entries of the selected variant
        tibble_selected_transcript_entries <- tibble_gtf_table[which(tibble_gtf_table$transcript_id == selected_transcript_id), ]
        
        ## extract stable HGNC variant name
        selected_hgnc_variant_name <- tibble_selected_transcript_entries$hgnc_stable_variant_ID %>% unique
        
        ## extract parent gene entries
        tibble_parent_gene_entries <- tibble_gtf_table[which(tibble_gtf_table$gene_id %in% (tibble_selected_transcript_entries$gene_id %>% unique)), ]
        
        # determine the lowest HGNC variant ID what contains an exon matching/if not then flanking/if not then overlapping with each of the query exons 
        tibble_VSR_exon_start_end_flagged_overlaps <- tibble_VSR_exon_start_end_flagged_IR %>% add_column("lowest_stable_variant_ID_matched_to_exon" = purrr::map2(
            .x = tibble_VSR_exon_start_end_flagged_IR$start, 
            .y = tibble_VSR_exon_start_end_flagged_IR$end, 
            .f = function(a1, a2) {
                
                tibble_matching <- extract_matching.exons(query_chr = query_chr, query_start = a1, query_end = a2, query_strand = query_strand, tibble_gtf_table = tibble_parent_gene_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon")
                
                tibble_flanking <- extract_junction.flanking.exons(query_chr = query_chr, query_start = a1, query_end = a2, query_strand = query_strand, tibble_gtf_table = tibble_parent_gene_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, match_consecutive = FALSE, return_type = "exon") %>% rbindlist
                
                tibble_overlapping <- extract_overlapping_features(query_chr = query_chr, query_start = a1, query_end = a2, query_strand = query_strand, tibble_gtf_table = tibble_parent_gene_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon")
                
                if (tibble_matching %>% nrow != 0) {
                    return(tibble_matching %>% .$hgnc_stable_variant_ID %>% mixedsort %>% .[1])
                } else if (tibble_flanking %>% nrow != 0) {
                    return(tibble_flanking %>% .$hgnc_stable_variant_ID %>% mixedsort %>% .[1])
                } else if (tibble_overlapping %>% nrow != 0) {
                    return(tibble_overlapping %>% .$hgnc_stable_variant_ID %>% mixedsort %>% .[1])
                }
                
            } ) %>% unlist )
        
        # update the strand information if it was not known before
        ## update based on the selected transcript strand. if that information is not available, then use the intersection of the magnetised VSR vertex.
        if (!query_strand %in% c("+", "-")) {
            
            if (exists("tibble_selected_transcript_entries") == TRUE) {
                
                query_strand <- tibble_selected_transcript_entries$strand %>% unique
                
            } else if (exists("tibble_selected_transcript_entries") == FALSE & exists("tibble_reference_entries_with_exon_ends_overlapping_VSR_left") == TRUE) {
                
                query_strand <- tibble_reference_entries_with_exon_ends_overlapping_VSR_left$strand %>% unique %>% sort(decreasing = TRUE) %>% .[1]
                
            } else if (exists("tibble_selected_transcript_entries") == FALSE & exists("tibble_reference_entries_with_exon_starts_overlapping_VSR_right") == TRUE) {
                
                query_strand <- tibble_reference_entries_with_exon_starts_overlapping_VSR_right$strand %>% unique %>% sort(decreasing = TRUE) %>% .[1]
                
            }
            
        }
        
    } else if (flag_is_intergenic == TRUE) {
        
        selected_hgnc_variant_name <- "Intergenic"
        
        tibble_VSR_exon_start_end_flagged_overlaps <- tibble_VSR_exon_start_end_flagged_IR %>% 
            add_column("lowest_stable_variant_ID_matched_to_exon" = "")
        
    }
    
    return(list(
        "query_chr" = query_chr,
        "query_VSR_start" = query_VSR_start %>% type.convert,
        "query_VSR_end" = query_VSR_end %>% type.convert,
        "query_strand" = query_strand,
        "effective_VSR_start" = effective_VSR_start %>% paste %>% type.convert,
        "effective_VSR_end" = effective_VSR_end %>% paste %>% type.convert,
        "effective_magnetised_VSR_start" = effective_magnetised_VSR_start,
        "effective_magnetised_VSR_end" = effective_magnetised_VSR_end,
        "selected_hgnc_variant_name" = selected_hgnc_variant_name,
        "flag_boundary_A5SS_detected" = flag_boundary_A5SS_detected,
        "flag_boundary_A3SS_detected" = flag_boundary_A3SS_detected,
        "flag_is_all_IR" = flag_is_all_IR,
        "flag_is_intergenic" = flag_is_intergenic,
        "flag_VSR_is_exactly_matched" = flag_VSR_is_exactly_matched,
        "flag_no_common_vertices" = flag_no_common_vertices,
        "tibble_body_exon_info" = tibble_VSR_exon_start_end_flagged_overlaps,
        "tibble_selected_transcript_entries" = tibble_selected_transcript_entries,
        "tibble_parent_gene_entries" = if (tibble_parent_gene_entries %>% is.null != TRUE) {tibble_parent_gene_entries} else {NA}
    ))
    
}

# END VSR_select_reference_transcript_variant() ###

## VSRs/LISs: FUNCTION TO NAME THE BODY EXONS ###

### takes the list input from VSR_select_reference_transcript_variant()
### we write this as a separate modular function because multiple LISs can be called, and we have to reconcile the selected variant names.
### so that means, we'll be doing a first-pass naming scheme for each individual LIS, then determining which one we want to roll with, then we refresh each LIS with the final choice.

### there can only be 5 scenarios: 
### 1. the exon matches exactly
### 2. there is exon truncation/extension
### 3. retained intron
### 4. the exon is completely unannotated
### 5. (unrelated to event exons) the reference transcript contains extra exons. this is not a problem for VSRs. however for FLI's we have to extensively use the Delta notation.

### if exact match, we don't need to do anything.
### otherwise, we need to find any mutual vertex overlap
### also test event vertex -> reference exon overlap
### then find the originating start/end reference vertex
### if they're from different exons, then it's an exon-joining event. 
### assign truncations/extensions accordingly (if vertex in exon, then use next further exon vertex (truncation). if vertex not in exon, use next closest exon vertex (extension))
### we need to pay attention to the flags. if it's a boundary A5/3SS then we have to check: 1. if there is an extended exon which extends past the VSR 2. If not, then it's a novel extension, so the VSR will change
### if the exon is IR, then we skip everything and directly call a join between the two flanking exons.

### for the rest, we can think of the nomenclature as slots. For example: -32(E6)+46 OR -32(E6)_(E8)+42
### slots1: -, 32, E6, +, 46
### slots2: -, 32, E6, +, 46, _, +, 56, E8, +, 42
### therefore we only need left/right truncation/extensions and the main exon which can only be a single exon or two joined.
### therefore, 1. test whether the query exon overlaps more than one exon. Call the middle slot. 2. for BOTH ends, if it's inside an exon, then it's a truncation. if it's outside an exon then it's an extension.

VSR_name_body_exons_and_VSR <- function(list_output_from_VSR_select_reference_transcript_variant, left_query_shift = 0, right_query_shift = 0, left_tolerance = 0, right_tolerance = 0) {
    
    # DEBUG ###
    
    # list_output_from_VSR_select_reference_transcript_variant <- VSR_select_reference_transcript_variant(VSR_coordinates = "16:89737976-89740903:*", tibble_VSR_exon_start_end = tribble(
    #     ~start, ~end,
    #     89740100, 89740803
    # ), left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1) 
    
    # list_output_from_VSR_select_reference_transcript_variant <- VSR_select_reference_transcript_variant(VSR_coordinates = "16:2757946-2768996:+", tibble_VSR_exon_start_end = tribble(
    #     ~start, ~end,
    #     2757946, 2757975,
    #     2758548, 2758984,
    #     2759140, 2759172,
    #     2760350, 2760499,
    #     2760593, 2760958,
    #     2761023, 2761661
    # ), left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1)
    
    # matches
    # "3:197051669-197065274:-"
    # "1:15152785-15214893:+"
    # no match
    # "3:197029997-197300102:-"
    
    # - strand    
    # 197161122, 197297255,
    # 197065155, 197065448,
    # 197059888, 197059998
    # list_output_from_VSR_select_reference_transcript_variant <- VSR_select_reference_transcript_variant(VSR_coordinates = "1:15152785-15214893:+", tibble_VSR_exon_start_end = tribble(
    #     ~start, ~end,
    #     15210489, 15210562
    # 
    # ), left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1)
    # 
    
    # list_output_from_VSR_select_reference_transcript_variant <- list_initial_lowest_variant_detection_force_lowest_variant_ID[[3]]
    # 
    # left_query_end_shift <- 0
    # right_query_end_shift <- 0
    # left_match_tolerance <- 1
    # right_match_tolerance <- 1
    
    # list_output_from_VSR_select_reference_transcript_variant <- list_selected_transcript_variant
    
    ###########
    
    # left_query_shift <<- left_query_shift
    # right_query_shift <<- right_query_shift
    # left_tolerance <<- left_tolerance
    # right_tolerance <<- right_tolerance
    
    tibble_body_exon_info <- list_output_from_VSR_select_reference_transcript_variant$tibble_body_exon_info
    
    if (list_output_from_VSR_select_reference_transcript_variant$flag_is_intergenic == FALSE) {
        
        # pool together all reference transcript vertices, starts and ends
        vector_selected_transcript_vertices <- c(list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries %>% dplyr::filter(type == "exon") %>% .$start, list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries %>% dplyr::filter(type == "exon") %>% .$end) %>% unlist %>% unique %>% sort
        
        # vector_selected_transcript_starts <- c(list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries %>% dplyr::filter(type == "exon") %>% .$start) %>% unlist %>% unique %>% sort
        
        # vector_selected_transcript_ends <- c(list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries %>% dplyr::filter(type == "exon") %>% .$end) %>% unlist %>% unique %>% sort
        
        # MATCH VSR TO SELECTED TRANSCRIPT
        ## we only need the equal to or one further than VSR ends for the exon number.
        ## for exon vertices, we have to consider all possibilities for it being an extension/truncation
        
        ## VSR VERTICES -> SELECTED EXON NUMBERS
        ## VSR VERTICES -> SELECTED VERTICES
        # VSR start
        if (list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_start %>% length > 0) {
            closest_LHS_exon_number_to_VSR_start <- list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries %>% .[.$start <= max(list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_start, .$start %>% min) & .$type == "exon", ] %>% .[.$start == max(.$start), "exon_number"] %>% paste %>% type.convert
            
            closest_LHS_reference_vertex_to_VSR_start <- max(vector_selected_transcript_vertices[vector_selected_transcript_vertices <= list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_start] %>% max,
                                                             list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries$end %>% min
            )
            
            closest_RHS_reference_vertex_to_VSR_start <- min(vector_selected_transcript_vertices[vector_selected_transcript_vertices >= list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_start] %>% min,
                                                             list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries$start %>% max
            )
            
            ## check whether each VSR vertex overlaps an exon
            ## note: if the VSR was matched then obviously it doesn't lie in an exon lool
            flag_VSR_start_is_exonic <- FALSE
        } else {
            closest_LHS_exon_number_to_VSR_start <- list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries %>% .[.$start <= max(list_output_from_VSR_select_reference_transcript_variant$effective_VSR_start, .$start %>% min) & .$type == "exon", ] %>% .[.$start == max(.$start), "exon_number"] %>% paste %>% type.convert
            
            closest_LHS_reference_vertex_to_VSR_start <- max(vector_selected_transcript_vertices[vector_selected_transcript_vertices <= list_output_from_VSR_select_reference_transcript_variant$effective_VSR_start] %>% max,
                                                             list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries$end %>% min
            )
            
            closest_RHS_reference_vertex_to_VSR_start <- min(vector_selected_transcript_vertices[vector_selected_transcript_vertices >= list_output_from_VSR_select_reference_transcript_variant$effective_VSR_start] %>% min,
                                                             list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries$start %>% max
            )
            
            ## check whether each VSR vertex overlaps an exon
            flag_VSR_start_is_exonic <- extract_overlapping_features(query_chr = list_output_from_VSR_select_reference_transcript_variant$query_chr, query_start = list_output_from_VSR_select_reference_transcript_variant$effective_VSR_start, query_end = list_output_from_VSR_select_reference_transcript_variant$effective_VSR_start, query_strand = "*", tibble_gtf_table = list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon") %>% nrow > 0
        }
        
        # VSR end
        if (list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_end %>% length > 0) {
            closest_RHS_exon_number_to_VSR_end <- list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries %>% .[.$end >= min(list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_end, .$end %>% max) & .$type == "exon", ] %>% .[.$end == min(.$end), "exon_number"] %>% paste %>% type.convert
            
            closest_LHS_reference_vertex_to_VSR_end <- max(vector_selected_transcript_vertices[vector_selected_transcript_vertices <= list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_end] %>% max,
                                                           list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries$end %>% min
            )
            
            closest_RHS_reference_vertex_to_VSR_end <- min(vector_selected_transcript_vertices[vector_selected_transcript_vertices >= list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_end] %>% min,
                                                           list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries$start %>% max
            )
            
            flag_VSR_end_is_exonic <- FALSE
        } else {
            closest_RHS_exon_number_to_VSR_end <- list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries %>% .[.$end >= min(list_output_from_VSR_select_reference_transcript_variant$effective_VSR_end, .$end %>% max) & .$type == "exon", ] %>% .[.$end == min(.$end), "exon_number"] %>% paste %>% type.convert
            
            closest_LHS_reference_vertex_to_VSR_end <- max(vector_selected_transcript_vertices[vector_selected_transcript_vertices <= list_output_from_VSR_select_reference_transcript_variant$effective_VSR_end] %>% max,
                                                           list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries$end %>% min
            )
            
            closest_RHS_reference_vertex_to_VSR_end <- min(vector_selected_transcript_vertices[vector_selected_transcript_vertices >= list_output_from_VSR_select_reference_transcript_variant$effective_VSR_end] %>% min,
                                                           list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries$start %>% max
            )
            
            ## check whether each VSR vertex overlaps an exon
            flag_VSR_end_is_exonic <- extract_overlapping_features(query_chr = list_output_from_VSR_select_reference_transcript_variant$query_chr, query_start = list_output_from_VSR_select_reference_transcript_variant$effective_VSR_end, query_end = list_output_from_VSR_select_reference_transcript_variant$effective_VSR_end, query_strand = "*", tibble_gtf_table = list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon") %>% nrow > 0
            
        }
        
        # MATCH BODY EXONS TO SELECTED TRANSCRIPT
        # match to check if IR (we have to do this because this function accepts transcript overrides)
        tibble_body_exon_info$IR <- purrr::map2(.x = tibble_body_exon_info$start, .y = tibble_body_exon_info$end, .f = ~extract_junction.flanking.exons(query_chr = list_output_from_VSR_select_reference_transcript_variant$query_chr, query_start = .x, query_end = .y, query_strand = "*", tibble_gtf_table = list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, match_consecutive = FALSE, return_type = "exon") %>% length > 0) %>% unlist
        
        list_output_from_VSR_select_reference_transcript_variant$flag_is_all_IR <- all(tibble_body_exon_info$IR == TRUE)
        
        ## test for exact exon matches in our selected transcript
        tibble_body_exon_info_matched_to_selected_transcript <- tibble_body_exon_info %>%
            dplyr::mutate("matches_to_selected_exons" = purrr::map2(.x = tibble_body_exon_info$start, .y = tibble_body_exon_info$end, .f = ~extract_matching.exons(query_chr = list_output_from_VSR_select_reference_transcript_variant$query_chr, query_start = .x, query_end = .y, query_strand = "*", tibble_gtf_table = list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon") %>% nrow > 0) %>% unlist)
        ## test for exon overlaps with our selected transcript
        tibble_body_exon_info_matched_to_selected_transcript <- tibble_body_exon_info_matched_to_selected_transcript %>%
            dplyr::mutate("overlaps_selected_exons" = purrr::map2(.x = tibble_body_exon_info_matched_to_selected_transcript$start, .y = tibble_body_exon_info_matched_to_selected_transcript$end, .f = ~extract_overlapping_features(query_chr = list_output_from_VSR_select_reference_transcript_variant$query_chr, query_start = .x, query_end = .y, query_strand = "*", tibble_gtf_table = list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon") %>% nrow > 0) %>% unlist)
        ## find the exon numbers overlapped by query exons
        tibble_body_exon_info_matched_to_selected_transcript <- tibble_body_exon_info_matched_to_selected_transcript %>%
            dplyr::mutate("exon_numbers_overlapped_by_query_exon" = purrr::map2(.x = tibble_body_exon_info_matched_to_selected_transcript$start, .y = tibble_body_exon_info_matched_to_selected_transcript$end, .f = ~extract_overlapping_features(query_chr = list_output_from_VSR_select_reference_transcript_variant$query_chr, query_start = .x, query_end = .y, query_strand = "*", tibble_gtf_table = list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon") %>% .$exon_number) )
        ## find the selected transcript vertices overlapped by query vertices
        tibble_body_exon_info_matched_to_selected_transcript <- tibble_body_exon_info_matched_to_selected_transcript %>%
            dplyr::mutate("selected_vertex_matched_to_query_start" = purrr::map(.x = tibble_body_exon_info_matched_to_selected_transcript$start, .f = ~vector_selected_transcript_vertices[vector_selected_transcript_vertices <= (.x + left_query_shift + left_tolerance) & vector_selected_transcript_vertices >= (.x + left_query_shift - left_tolerance)] ) %>% purrr::map_if(.p = ~.x %>% length == 0, .f = ~NA) %>% unlist %>% type.convert,
                          "selected_vertex_matched_to_query_end" = purrr::map(.x = tibble_body_exon_info_matched_to_selected_transcript$end, .f = ~vector_selected_transcript_vertices[vector_selected_transcript_vertices <= (.x + right_query_shift + right_tolerance) & vector_selected_transcript_vertices >= (.x + right_query_shift - right_tolerance)] ) %>% purrr::map_if(.p = ~.x %>% length == 0, .f = ~NA) %>% unlist %>% type.convert)
        ## find the exon numbers overlapped by query vertices
        tibble_body_exon_info_matched_to_selected_transcript <- tibble_body_exon_info_matched_to_selected_transcript %>%
            dplyr::mutate("exon_numbers_overlapped_by_query_start" = purrr::map(.x = tibble_body_exon_info_matched_to_selected_transcript$start, .f = ~extract_overlapping_features(query_chr = list_output_from_VSR_select_reference_transcript_variant$query_chr, query_start = .x, query_end = .x, query_strand = "*", tibble_gtf_table = list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon") %>% .$exon_number) %>% purrr::map_if(.p = ~.x %>% length == 0, .f = ~NA) %>% unlist )
        #
        tibble_body_exon_info_matched_to_selected_transcript <- tibble_body_exon_info_matched_to_selected_transcript %>%
            dplyr::mutate("exon_numbers_overlapped_by_query_end" = purrr::map(.x = tibble_body_exon_info_matched_to_selected_transcript$end, .f = ~extract_overlapping_features(query_chr = list_output_from_VSR_select_reference_transcript_variant$query_chr, query_start = .x, query_end = .x, query_strand = "*", tibble_gtf_table = list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon") %>% .$exon_number) %>% purrr::map_if(.p = ~.x %>% length == 0, .f = ~NA) %>% unlist )
        ## find the nearest exon numbers upstream and downstream of the query vertices
        tibble_body_exon_info_matched_to_selected_transcript <- tibble_body_exon_info_matched_to_selected_transcript %>%
            dplyr::mutate("closest_upstream_exon_number_to_query_start" = purrr::map(.x = tibble_body_exon_info_matched_to_selected_transcript$start, .f = ~list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries %>% .[.$start <= .x & .$type == "exon", ] %>% .[.$start == max(.$start), "exon_number"] %>% paste) %>% unlist %>% type.convert,
                          "closest_downstream_exon_number_to_query_start" = purrr::map(.x = tibble_body_exon_info_matched_to_selected_transcript$start, .f = ~list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries %>% .[.$end >= .x & .$type == "exon", ] %>% .[.$end == min(.$end), "exon_number"] %>% paste) %>% unlist %>% type.convert,
                          
                          "closest_upstream_exon_number_to_query_end" = purrr::map(.x = tibble_body_exon_info_matched_to_selected_transcript$end, .f = ~list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries %>% .[.$start <= .x & .$type == "exon", ] %>% .[.$start == max(.$start), "exon_number"] %>% paste) %>% unlist %>% type.convert,
                          "closest_downstream_exon_number_to_query_end" = purrr::map(.x = tibble_body_exon_info_matched_to_selected_transcript$end, .f = ~list_output_from_VSR_select_reference_transcript_variant$tibble_selected_transcript_entries %>% .[.$end >= .x & .$type == "exon", ] %>% .[.$end == min(.$end), "exon_number"] %>% paste) %>% unlist %>% type.convert)
        ## find the first reference exon vertices upstream and downstream of the query start/ends
        tibble_body_exon_info_matched_to_selected_transcript <- tibble_body_exon_info_matched_to_selected_transcript %>%
            dplyr::mutate("closest_LHS_reference_vertex_to_query_start" = purrr::map(.x = tibble_body_exon_info_matched_to_selected_transcript$start, .f = ~vector_selected_transcript_vertices[vector_selected_transcript_vertices <= .x] %>% max) %>% unlist,
                          "closest_RHS_reference_vertex_to_query_start" = purrr::map(.x = tibble_body_exon_info_matched_to_selected_transcript$start, .f = ~vector_selected_transcript_vertices[vector_selected_transcript_vertices >= .x] %>% min) %>% unlist,
                          
                          "closest_LHS_reference_vertex_to_query_end" = purrr::map(.x = tibble_body_exon_info_matched_to_selected_transcript$end, .f = ~vector_selected_transcript_vertices[vector_selected_transcript_vertices <= .x] %>% max) %>% unlist,
                          "closest_RHS_reference_vertex_to_query_end" = purrr::map(.x = tibble_body_exon_info_matched_to_selected_transcript$end, .f = ~vector_selected_transcript_vertices[vector_selected_transcript_vertices >= .x] %>% min) %>% unlist)
        
        # account for strand. if - strand, then we gona be flipping all those upstream/downstream column names around
        if (list_output_from_VSR_select_reference_transcript_variant$query_strand == "-") {
            
            tibble_body_exon_info_matched_to_selected_transcript <- tibble_body_exon_info_matched_to_selected_transcript %>% 
                dplyr::rename("closest_downstream_exon_number_to_query_start" = "closest_upstream_exon_number_to_query_start",
                              "closest_upstream_exon_number_to_query_start" = "closest_downstream_exon_number_to_query_start",
                              
                              "closest_downstream_exon_number_to_query_end" = "closest_upstream_exon_number_to_query_end",
                              "closest_upstream_exon_number_to_query_end" = "closest_downstream_exon_number_to_query_end")
            
        }
        
    } else if (list_output_from_VSR_select_reference_transcript_variant$flag_is_intergenic == TRUE) {
        
        tibble_body_exon_info_matched_to_selected_transcript <- tibble_body_exon_info
        
    }
    
    # GENERATE BODY EXON NOMENCLATURE BLOCKS
    # if not intergenic...
    if (list_output_from_VSR_select_reference_transcript_variant$flag_is_intergenic == FALSE) {
        
        # loop thru the whole table and generate nomenclature
        list_body_exon_nomenclature <- purrr::map(
            .x = tibble_body_exon_info_matched_to_selected_transcript %>% purrr::array_tree(),
            .f = function(a1) { 
                
                # DEBUG ###
                # a1 <- tibble_body_exon_info_matched_to_selected_transcript %>% purrr::array_tree() %>% .[[1]]
                ###########
                
                # account for boundary A3/A5SS events
                # return NA immediately
                if ((a1$left_end_of_VSR == TRUE | a1$right_end_of_VSR == TRUE) & a1$exon_matches_to_any_reference == FALSE) {
                    
                    return(NA)
                    
                    # account for IR exons
                } else if (a1$IR == TRUE) {
                    
                    return(paste("E", min(a1$closest_upstream_exon_number_to_query_start, a1$closest_downstream_exon_number_to_query_end), "_E", max(a1$closest_upstream_exon_number_to_query_start, a1$closest_downstream_exon_number_to_query_end), sep = ""))
                    
                    # account for completely new exons in intronic regions
                } else if (a1$exon_numbers_overlapped_by_query_exon %>% length == 0) {
                    
                    if (list_output_from_VSR_select_reference_transcript_variant$query_strand == "-") {
                        
                        return(paste("(", list_output_from_VSR_select_reference_transcript_variant$effective_VSR_end - a1$end + 1, "_", list_output_from_VSR_select_reference_transcript_variant$effective_VSR_end - a1$start + 1, ")", sep = ""))
                        
                    } else {
                        
                        return(paste("(", a1$start - list_output_from_VSR_select_reference_transcript_variant$effective_VSR_start + 1, "_", a1$end - list_output_from_VSR_select_reference_transcript_variant$effective_VSR_start + 1, ")", sep = ""))
                        
                    }
                    
                    # finally, we deal with everything else that remains.   
                } else {
                    
                    # MIDDLE SLOT
                    ## check how many selected exons are covered by the query exon
                    if (a1$exon_numbers_overlapped_by_query_exon %>% length == 1) {
                        middle_slot <- paste("E", a1$exon_numbers_overlapped_by_query_exon, sep = "")
                    } else if (a1$exon_numbers_overlapped_by_query_exon %>% length > 1) {
                        middle_slot <- paste("E", a1$exon_numbers_overlapped_by_query_exon %>% min, "_E", a1$exon_numbers_overlapped_by_query_exon %>% max, sep = "")
                    }
                    
                    # LEFT SLOT (QUERY EXON START)
                    ## check if it matched any of the selected vertices
                    ## if it matched to any of the selected vertices, then the slot is empty.
                    ## if not, then check whether it's in an intronic or exonic region and call truncation/extension accordingly
                    if (is.na(a1$selected_vertex_matched_to_query_start) == FALSE) {
                        
                        left_slot <- ""
                        
                    } else if (is.na(a1$selected_vertex_matched_to_query_start) == TRUE) {
                        
                        # INTRONIC QUERY VERTEX - LEFT EXTENSION
                        if (is.na(a1$exon_numbers_overlapped_by_query_start) == TRUE) {
                            
                            left_slot <- paste("+", a1$closest_RHS_reference_vertex_to_query_start - a1$start, sep = "")
                            
                            # EXONIC QUERY VERTEX - LEFT TRUNCATION   
                        } else if (is.na(a1$exon_numbers_overlapped_by_query_start) == FALSE) {
                            
                            left_slot <- paste("–", a1$start - a1$closest_LHS_reference_vertex_to_query_start, sep = "")
                            
                        }
                        
                    }
                    
                    # RIGHT SLOT (QUERY EXON END)
                    ## check if it matched any of the selected vertices
                    ## if it matched to any of the selected vertices, then the slot is empty.
                    ## if not, then check whether it's in an intronic or exonic region and call truncation/extension accordingly
                    if (is.na(a1$selected_vertex_matched_to_query_end) == FALSE) {
                        
                        right_slot <- ""
                        
                    } else if (is.na(a1$selected_vertex_matched_to_query_end) == TRUE) {
                        
                        # INTRONIC QUERY VERTEX - RIGHT EXTENSION
                        if (is.na(a1$exon_numbers_overlapped_by_query_end) == TRUE) {
                            
                            right_slot <- paste("+", a1$end - a1$closest_LHS_reference_vertex_to_query_end, sep = "")
                            
                            # EXONIC QUERY VERTEX - RIGHT TRUNCATION   
                        } else if (is.na(a1$exon_numbers_overlapped_by_query_end) == FALSE) {
                            
                            right_slot <- paste("–", a1$closest_RHS_reference_vertex_to_query_end - a1$end, sep = "")
                            
                        }
                        
                    }
                    
                    # JOIN SLOTS FOR THE FINAL NAME MODULE
                    
                    if (list_output_from_VSR_select_reference_transcript_variant$query_strand == "-") {
                        
                        return(paste(right_slot, middle_slot, left_slot, sep = ""))
                        
                    } else {
                        
                        return(paste(left_slot, middle_slot, right_slot, sep = ""))
                        
                    }
                    
                }
                
            } )
        
        # intergenic exons - output genomic coords
    } else if (list_output_from_VSR_select_reference_transcript_variant$flag_is_intergenic == TRUE) {
        
        if (list_output_from_VSR_select_reference_transcript_variant$query_strand == "-") {
            
            list_body_exon_nomenclature <- paste("(", list_output_from_VSR_select_reference_transcript_variant$effective_VSR_end - tibble_body_exon_info_matched_to_selected_transcript$end + 1, "_", list_output_from_VSR_select_reference_transcript_variant$effective_VSR_end - tibble_body_exon_info_matched_to_selected_transcript$start + 1, ")", sep = "") %>% array_tree
            
        } else {
            
            list_body_exon_nomenclature <- paste("(", tibble_body_exon_info_matched_to_selected_transcript$start - list_output_from_VSR_select_reference_transcript_variant$effective_VSR_start + 1, "_", tibble_body_exon_info_matched_to_selected_transcript$end - list_output_from_VSR_select_reference_transcript_variant$effective_VSR_start + 1, ")", sep = "") %>% array_tree
            
        }
        
    }
    
    # NAME THE VSR
    ## collect the VSR start and end into the magnetised values
    if (list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_start %>% length == 0) {
        effective_magnetised_VSR_start <- list_output_from_VSR_select_reference_transcript_variant$effective_VSR_start
    } else {
        effective_magnetised_VSR_start <- list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_start
    }
    
    if (list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_end %>% length == 0) {
        effective_magnetised_VSR_end <- list_output_from_VSR_select_reference_transcript_variant$effective_VSR_end
    } else {
        effective_magnetised_VSR_end <- list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_end
    }
    
    ## rematch the VSR if a match was originally found.
    ## CHECK FOR ALL IR FIRST. If all IR, then we remove the VSR because it doesn't tell us anything.
    ## if match was not found, then the exons would have determined the selected transcript. we have to modify the exon ends.
    ## Otherwise, if the splice event was intergenic, then the VSR will just be coordinates.
    # match effective VSR to the parent gene
    if (list_output_from_VSR_select_reference_transcript_variant$flag_is_all_IR == TRUE) {
        
        upstream_VSR_slot <- ""
        downstream_VSR_slot <- ""
        
        # upstream_VSR_slot <- paste("E", min(closest_LHS_exon_number_to_VSR_start, closest_RHS_exon_number_to_VSR_end), sep = "")
        # 
        # downstream_VSR_slot <- paste("E", max(closest_LHS_exon_number_to_VSR_start, closest_RHS_exon_number_to_VSR_end), sep = "")
        
    } else if (list_output_from_VSR_select_reference_transcript_variant$flag_is_intergenic == FALSE) {
        
        # if no exact VSR match, look for the nearest upstream/downstream exon end/start to the VSR start/end and modify the exons.
        if (list_output_from_VSR_select_reference_transcript_variant$flag_VSR_is_exactly_matched == FALSE) {
            
            # retrieve proximal exon number
            left_VSR_exon_slot <- paste("E", closest_LHS_exon_number_to_VSR_start, sep = "")
            right_VSR_exon_slot <- paste("E", closest_RHS_exon_number_to_VSR_end, sep = "")
            
            # test exonic for truncation/extension
            ## left VSR slot
            if (flag_VSR_start_is_exonic == TRUE) {
                left_VSR_modification_slot <- paste("–", closest_RHS_reference_vertex_to_VSR_start - effective_magnetised_VSR_start, sep = "")
            } else if (flag_VSR_start_is_exonic == FALSE) {
                left_VSR_modification_slot <- paste("+", effective_magnetised_VSR_start - closest_LHS_reference_vertex_to_VSR_start, sep = "")
            }
            
            ## right VSR slot
            if (flag_VSR_end_is_exonic == TRUE) {
                right_VSR_modification_slot <- paste("–", effective_magnetised_VSR_end - closest_LHS_reference_vertex_to_VSR_end, sep = "")
            } else if (flag_VSR_end_is_exonic == FALSE) {
                right_VSR_modification_slot <- paste("+", closest_RHS_reference_vertex_to_VSR_end - effective_magnetised_VSR_end, sep = "")
            }
            
            ## create the VSR slots
            ### MODIFIERS MUST ALWAYS BE IN THE MIDDLE (because they're meeting the VSR)
            if (list_output_from_VSR_select_reference_transcript_variant$query_strand == "-") {
                upstream_VSR_slot <- paste(right_VSR_exon_slot, right_VSR_modification_slot, sep = "")
                downstream_VSR_slot <- paste(left_VSR_modification_slot, left_VSR_exon_slot, sep = "")
            } else {
                upstream_VSR_slot <- paste(left_VSR_exon_slot, left_VSR_modification_slot, sep = "")
                downstream_VSR_slot <- paste(right_VSR_modification_slot, right_VSR_exon_slot, sep = "")
            }
            
            # now account for if the left VSR was exactly on the start of an exon
            if (closest_RHS_reference_vertex_to_VSR_start - effective_magnetised_VSR_start == 0) {
                if (list_output_from_VSR_select_reference_transcript_variant$query_strand == "-") {
                    downstream_VSR_slot <- paste("(", left_VSR_exon_slot, ")3'", sep = "")
                } else {
                    upstream_VSR_slot <- paste("5'(", left_VSR_exon_slot, ")", sep = "")
                }
            }
            # same thing for the VSR right
            # now account for if the left VSR was exactly on the start of an exon
            if (effective_magnetised_VSR_end - closest_LHS_reference_vertex_to_VSR_end == 0) {
                if (list_output_from_VSR_select_reference_transcript_variant$query_strand == "-") {
                    upstream_VSR_slot <- paste("5'(", right_VSR_exon_slot, ")", sep = "")
                } else {
                    downstream_VSR_slot <- paste("(", right_VSR_exon_slot, ")3'", sep = "")
                }
            }
            
        }
        
        # override exact exon onto the matched start if a match was found on the start
        if (list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_start %>% length > 0) {
            
            if (list_output_from_VSR_select_reference_transcript_variant$query_strand == "-") {
                downstream_VSR_slot <- paste("E", closest_LHS_exon_number_to_VSR_start, sep = "")
            } else {
                upstream_VSR_slot <- paste("E", closest_LHS_exon_number_to_VSR_start, sep = "")
            }
            
        }
        # override exact exon onto the matched end if a match was found on the end
        if (list_output_from_VSR_select_reference_transcript_variant$effective_magnetised_VSR_end %>% length > 0) {
            
            if (list_output_from_VSR_select_reference_transcript_variant$query_strand == "-") {
                upstream_VSR_slot <- paste("E", closest_RHS_exon_number_to_VSR_end, sep = "")
            } else {
                downstream_VSR_slot <- paste("E", closest_RHS_exon_number_to_VSR_end, sep = "")
            }
            
        }
        
        # deal with VSRs extending past the transcript boundaries
        upstream_VSR_slot <- gsub(x = upstream_VSR_slot, pattern = "(.*)\\+\\–(.*)", replacement = "(\\1–\\2)3'")
        downstream_VSR_slot <- gsub(x = downstream_VSR_slot, pattern = "\\+\\–(.*)", replacement = "5'(–\\1)")
        
        # intergenic condition    
    } else if (list_output_from_VSR_select_reference_transcript_variant$flag_is_intergenic == TRUE) {
        
        if (list_output_from_VSR_select_reference_transcript_variant$query_strand == "-") {
            upstream_VSR_slot <- paste(effective_magnetised_VSR_end, "(g)", sep = "")
            downstream_VSR_slot <- paste(effective_magnetised_VSR_start, "(g)", sep = "")
        } else {
            upstream_VSR_slot <- paste(effective_magnetised_VSR_start, "(g)", sep = "")
            downstream_VSR_slot <- paste(effective_magnetised_VSR_end, "(g)", sep = "")
        }
        
    }
    
    return(list_output_from_VSR_select_reference_transcript_variant %>% 
               purrr::modify_at(.at = "tibble_body_exon_info", .f = ~tibble_body_exon_info_matched_to_selected_transcript) %>%
               purrr::splice(
                   "list_body_exon_nomenclature" = if (list_output_from_VSR_select_reference_transcript_variant$query_strand == "-") {
                       list_body_exon_nomenclature %>% rev %>% list
                   } else {
                       list_body_exon_nomenclature %>% list
                   },
                   "upstream_VSR_slot" = upstream_VSR_slot,
                   "downstream_VSR_slot" = downstream_VSR_slot
               ))
    
}

# END VSR_name_body_exons_and_VSR() ###

## FLIs: FUNCTION TO ORGANISE THE TRANSCRIPT VARIANT CALLING AND NAMING ###
### Behaviour: accepts the chr/start/end/strand coords of the FLI in a tibble
### find the transcript with the greatest number of vertices in common.
### name using the VSR function
### this is essentially the same as VSR naming but without the complexities of VSRs
### this function automates the process for a single reassembled transcript. Parallelisation is called from the server-side Shiny code.
### Behaviour: 1. find transcripts overlapping with the isoform, 2. find transcripts with the most number of common vertices, if not then those with the transcript start/end closest to the query isoform start/end will be chosen. 3. pick lowest variant number 4. pass through to the VSR naming function. 5. Return name.

FLI_organise_matching <- function(tibble_FLI_chr_start_end_strand, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 0, right_tolerance = 0) {
    
    # DEBUG ###
    
    # tibble_FLI_chr_start_end_strand <- rtracklayer::import(con = "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/BM_MSC_to_OB_12d_r1_denovo_reconstructed_transcriptome.gtf", format = "gtf") %>% as_tibble %>% .[.$transcript_id == "BM_MSC_to_OB_12d_r1_Aligned.100.1" & .$type == "exon", ] %>% dplyr::rename("chr" = "seqnames") %>% dplyr::mutate_if(is.factor, as.character)
    # left_query_shift <- 0
    # right_query_shift <- 0
    # left_tolerance <- 1
    # right_tolerance <- 1
    # tibble_gtf_table <- tibble_ref_gtf
    
    ###########
    
    # left_query_shift <- left_query_shift %>% paste %>% type.convert
    # right_query_shift <- right_query_shift %>% paste %>% type.convert
    # left_tolerance <- left_tolerance %>% paste %>% type.convert
    # right_tolerance <- right_tolerance %>% paste %>% type.convert
    # 
    # print(left_query_shift)
    # print(right_query_shift)
    # print(left_tolerance)
    # print(right_tolerance)
    
    # DETERMINE THE EXTENT OF THE TRANSCRIPT
    isoform_chr <- tibble_FLI_chr_start_end_strand$chr %>% unique
    isoform_strand <- tibble_FLI_chr_start_end_strand$strand %>% unique
    # extract isoform start/end
    isoform_start <- tibble_FLI_chr_start_end_strand$start %>% min
    isoform_end <- tibble_FLI_chr_start_end_strand$end %>% max
    
    # EXON TIBBLE
    # sort the exon tibble in order of increasing co-ordinates
    tibble_FLI_chr_start_end_strand <- tibble_FLI_chr_start_end_strand %>% dplyr::arrange(start)
    
    # detect intron retention
    ## we will match every exon in the tibble to check for whether they are junction-spanning.
    tibble_FLI_chr_start_end_strand_flagged_IR <- tibble_FLI_chr_start_end_strand %>%
        dplyr::mutate("IR" = purrr::map2(.x = tibble_FLI_chr_start_end_strand$start, .y = tibble_FLI_chr_start_end_strand$end, .f = ~extract_junction.flanking.exons(query_chr = isoform_chr, query_start = .x, query_end = .y, query_strand = isoform_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, match_consecutive = FALSE, return_type = "exon") %>% length > 0) %>% unlist,
                      "left_end_of_VSR" = FALSE,
                      "right_end_of_VSR" = FALSE,
                      "exon_matches_to_any_reference" = FALSE)
    
    ## find reference transcripts which overlap the query isoform
    tibble_trancripts_overlapping_query_isoform <- extract_overlapping_features(query_chr = isoform_chr, query_start = isoform_start, query_end = isoform_end, query_strand = isoform_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "transcript")
    
    vector_trancript_ids_overlapping_query_isoform <- tibble_trancripts_overlapping_query_isoform$transcript_id %>% unique
    
    ## extract parent entries of overlapped transcripts
    tibble_isoform_overlapped_transcripts_parent_entries <- tibble_gtf_table[tibble_gtf_table$transcript_id %in% vector_trancript_ids_overlapping_query_isoform, ]
    
    # generate isoform overlap flags
    ## if there was only one overlapping transcript, then we can only use that transcript.
    if (vector_trancript_ids_overlapping_query_isoform %>% length == 1) {
        
        isoform_selected_transcript_id <- tibble_trancripts_overlapping_query_isoform$transcript_id
        
        flag_isoform_overlaps_ref_transcripts <- TRUE
        
        flag_is_intergenic <- FALSE
        
        # VSR MATCH TO MORE THAN ONE TRANSCRIPT
    } else if (vector_trancript_ids_overlapping_query_isoform %>% length > 1) {
        
        flag_isoform_overlaps_ref_transcripts <- TRUE
        
        flag_is_intergenic <- FALSE
        
        # VSR NO MATCH TO ANY TRANSCRIPTS
        # condition to account for VSRs which didn't match to any transcript
    } else if (vector_trancript_ids_overlapping_query_isoform %>% length == 0) {
        
        # test for intergenic VSR
        ## if not intergenic, then existing exons can be modified.
        
        flag_isoform_overlaps_ref_transcripts <- FALSE
        
        flag_is_intergenic <- TRUE
        
    }
    
    
    
    # EXON MATCHING
    # if not intergenic, then we can proceed as normal.
    if (flag_is_intergenic == FALSE) {
        
        # consolidate vertices
        ## these consist of exon start and end vertices
        ## bona fide A3/5SS extensions at the ends (common end with VSR + not boundary skipped exon) will have only the interior vertex kept.
        ## IR exons will have their start/ends swapped and expanded 
        vector_consolidated_exon_starts <- c(
            tibble_FLI_chr_start_end_strand_flagged_IR[tibble_FLI_chr_start_end_strand_flagged_IR$IR == FALSE, "start"], 
            tibble_FLI_chr_start_end_strand_flagged_IR[tibble_FLI_chr_start_end_strand_flagged_IR$IR == TRUE, "end"] + 1) %>% unlist
        
        vector_consolidated_exon_ends <- c(
            tibble_FLI_chr_start_end_strand_flagged_IR[tibble_FLI_chr_start_end_strand_flagged_IR$IR == FALSE, "end"], 
            tibble_FLI_chr_start_end_strand_flagged_IR[tibble_FLI_chr_start_end_strand_flagged_IR$IR == TRUE, "start"] - 1) %>% unlist
        
        # subset reference tibble for only exonic entries which contain the consolidated vertices
        if (isoform_strand %in% c("+", "-")) {
            
            tibble_ref_gtf_exons_with_common_vertices <- tibble_gtf_table[tibble_gtf_table$seqnames == isoform_chr, ] %>% .[.$strand == isoform_strand, ] %>% .[.$start %in% vector_consolidated_exon_starts | .$end %in% vector_consolidated_exon_ends, ] %>% .[.$type == "exon" & .$transcript_id %in% vector_trancript_ids_overlapping_query_isoform, ]
            
        } else if (!isoform_strand %in% c("+", "-")) {
            
            tibble_ref_gtf_exons_with_common_vertices <- tibble_gtf_table[tibble_gtf_table$seqnames == isoform_chr, ] %>% .[.$start %in% vector_consolidated_exon_starts | .$end %in% vector_consolidated_exon_ends, ] %>% .[.$type == "exon" & .$transcript_id %in% vector_trancript_ids_overlapping_query_isoform, ]
            
        }
        
        # ACCOUNT FOR THE POSSIBILITY OF NO EXONS IN COMMON
        ## IF: 
        ## - exons are completely novel BUT the isoform overlaps with a transcript: take the lowest matched transcript number. The isoform exons will be defined in terms of the isoform boundary. It will look like (SRRM2-v1 E2) E(2344_2389) E(545986_546433) (SRRM2-v1 E3) OR (SRRM-v1 E2) J(2343nt) E(46nt) etc...
        ## - exons are completely novel AND the VSR is not matched: take the transcript which has exons closest to the VSR boundaries. The VsR boundary will be defined in terms of this. The internal exons will be defined in terms of the VSR boundary.
        ## these will be dealt with later.
        if (tibble_ref_gtf_exons_with_common_vertices %>% nrow == 0) {
            
            flag_no_common_vertices <- TRUE
            
            selected_transcript_id <- mixedsort(vector_trancript_ids_overlapping_query_isoform) %>% .[1]
            
        } else if (tibble_ref_gtf_exons_with_common_vertices %>% nrow > 0) {
            
            flag_no_common_vertices <- FALSE
            
            # list-ify by transcript_id and test for the number of vertices in common
            ## group_split
            list_ref_gtf_exons_with_common_vertices_by_transcript_id <- tibble_ref_gtf_exons_with_common_vertices %>% 
                dplyr::group_split(transcript_id) %>% 
                set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$transcript_id %>% unique) %>% unlist)
            ## test for number of vertices in common
            list_number_of_reference_vertices_in_common <- purrr::map(
                .x = list_ref_gtf_exons_with_common_vertices_by_transcript_id,
                .f = ~intersect(c(.x$start, .x$end) %>% unique, c(vector_consolidated_exon_starts, vector_consolidated_exon_ends) %>% unique) %>% length
            )
            ## percolate and tibblise
            tibble_common_vertices_count <- list_number_of_reference_vertices_in_common %>% as_tibble(.name_repair = "unique") %>% t %>% as_tibble(rownames = "transcript_id") %>% setNames(c("transcript_id", "number_vertices_in_common"))
            ## extract the lowest transcript_id with the highest common vertices
            selected_transcript_id <- tibble_common_vertices_count[tibble_common_vertices_count$number_vertices_in_common == max(tibble_common_vertices_count$number_vertices_in_common), "transcript_id"] %>% unlist %>% mixedsort %>% .[1]
            
        }
        
        # if isoform range authority exists, then override body exon matches.
        if (exists("isoform_selected_transcript_id") == TRUE) {
            
            selected_transcript_id <- isoform_selected_transcript_id
            
        }
        
        ## extract all the entries of the selected variant
        tibble_selected_transcript_entries <- tibble_gtf_table[which(tibble_gtf_table$transcript_id == selected_transcript_id), ]
        
        ## extract stable HGNC variant name
        selected_hgnc_variant_name <- tibble_selected_transcript_entries$hgnc_stable_variant_ID %>% unique
        
        ## extract parent gene entries
        tibble_parent_gene_entries <- tibble_gtf_table[which(tibble_gtf_table$gene_id %in% (tibble_selected_transcript_entries$gene_id %>% unique)), ]
        
    } else if (flag_is_intergenic == TRUE) {
        
        flag_no_common_vertices <- TRUE
        
        selected_hgnc_variant_name <- "Intergenic"
        
    }
    
    # update the strand information if it was not known before
    ## update based on the selected transcript strand. if that information is not available, then use the intersection of the magnetised VSR vertex.
    if (!isoform_strand %in% c("+", "-")) {
        
        if (selected_hgnc_variant_name != "Intergenic") {
            
            isoform_strand <- tibble_selected_transcript_entries$strand %>% unique
            
        } else if (selected_hgnc_variant_name == "Intergenic" & exists("tibble_reference_entries_with_exon_ends_overlapping_VSR_left") == TRUE) {
            
            isoform_strand <- tibble_reference_entries_with_exon_ends_overlapping_VSR_left$strand %>% unique %>% .[1]
            
        } else if (selected_hgnc_variant_name == "Intergenic" & exists("tibble_reference_entries_with_exon_starts_overlapping_VSR_right") == TRUE) {
            
            isoform_strand <- tibble_reference_entries_with_exon_starts_overlapping_VSR_right$strand %>% unique %>% .[1]
            
        }
        
    }
    
    list_selected_transcript_variant <- list(
        "query_chr" = isoform_chr,
        "query_VSR_start" = isoform_start %>% type.convert,
        "query_VSR_end" = isoform_end %>% type.convert,
        "query_strand" = isoform_strand,
        "effective_VSR_start" = isoform_start %>% type.convert,
        "effective_VSR_end" = isoform_end %>% type.convert,
        "effective_magnetised_VSR_start" = isoform_start %>% type.convert,
        "effective_magnetised_VSR_end" = isoform_end %>% type.convert,
        "selected_hgnc_variant_name" = selected_hgnc_variant_name,
        "flag_boundary_A5SS_detected" = FALSE,
        "flag_boundary_A3SS_detected" = FALSE,
        "flag_is_all_IR" = FALSE,
        "flag_is_intergenic" = flag_is_intergenic,
        "flag_VSR_is_exactly_matched" = FALSE,
        "flag_isoform_overlaps_ref_transcripts" = flag_isoform_overlaps_ref_transcripts,
        "flag_no_common_vertices" = flag_no_common_vertices,
        "tibble_body_exon_info" = tibble_FLI_chr_start_end_strand_flagged_IR,
        "tibble_selected_transcript_entries" = if (exists("tibble_selected_transcript_entries") == TRUE) {tibble_selected_transcript_entries} else {NA},
        "tibble_parent_gene_entries" = if (exists("tibble_parent_gene_entries") == TRUE) {tibble_parent_gene_entries} else {NA}
    )
    
    # first pass naming using the lowest variant ID for ALL exons.
    list_first_pass_naming <- VSR_name_body_exons_and_VSR(list_output_from_VSR_select_reference_transcript_variant = list_selected_transcript_variant, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance)
    
    if (flag_is_intergenic == FALSE) {
        
        # RETURN NOMENCLATURE
        ## collapse by spaces between exons of the same LIS
        ## NOTE: ONLY RETURN THE EXONS WHICH ARE DIFFERENT AND DIDN'T MATCH TO THE SELECTED TRANSCRIPT.
        ## ALSO PUT A DELTA ON THE EXONS WHICH 
        list_isoform_only_nomenclature_event_raw <- list_first_pass_naming$list_body_exon_nomenclature %>% unlist
        # create a vector of logicals indicating which isoform exons didnt have an exact match. we only use these in the naming. the rest are omitted.
        if (isoform_strand == "-") {
            vector_isoform_exons_with_exact_match <- which(list_first_pass_naming$tibble_body_exon_info$matches_to_selected_exons %>% rev == FALSE)
        } else {
            vector_isoform_exons_with_exact_match <- which(list_first_pass_naming$tibble_body_exon_info$matches_to_selected_exons == FALSE)
        }
        
        # get a vector of exon numbers in the transcripts so we can see which is missing from the isoform
        vector_selected_exon_numbers <- list_first_pass_naming$tibble_selected_transcript_entries$exon_number %>% na.omit %>% unique
        # get a vector of selected exon numbers overlapped by the isoform
        vector_selected_exons_overlapped_by_isoform <- list_first_pass_naming$tibble_body_exon_info$exon_numbers_overlapped_by_query_exon %>% unlist %>% na.omit %>% unique
        # take the set difference. these are the deltas.
        vector_selected_exon_number_not_in_event <- setdiff(vector_selected_exon_numbers, vector_selected_exons_overlapped_by_isoform)
        
        # create the delta slots and final event-only exons
        vector_isoform_only_nomenclature_simplified <- list_isoform_only_nomenclature_event_raw[vector_isoform_exons_with_exact_match]
        if (vector_selected_exon_number_not_in_event %>% length > 0) {
            vector_delta_slots <- paste("ΔE", vector_selected_exon_number_not_in_event, sep = "")
        } else {
            vector_delta_slots <- ""
        }
        
        # mix and re-sort in proper exon order
        vector_combined_nomenclature <- c(vector_isoform_only_nomenclature_simplified, vector_delta_slots)
        # gsub for the exon number
        vector_combined_exon_numbers_only <- gsub(x = vector_combined_nomenclature, pattern = ".*E([0-9]+)(.*){0,1}", replacement = "\\1") %>% type.convert
        
        # sort by exon number
        tibble_sorted_combined_nomenclature <- tibble("slots" = vector_combined_nomenclature, "exon_numbers" = vector_combined_exon_numbers_only) %>% dplyr::arrange(exon_numbers)
        
        # finally, extract the transcript version
        variant_slot <- paste(list_first_pass_naming$selected_hgnc_variant_name %>% na.omit %>% unique, ".", tibble_selected_transcript_entries[tibble_selected_transcript_entries$hgnc_stable_variant_ID == list_first_pass_naming$selected_hgnc_variant_name %>% na.omit %>% unique, "transcript_version"] %>% unlist %>% na.omit %>% unique %>% .[1], sep = "")
        
        final_VSR_nomenclature <- paste(variant_slot, "t ", tibble_sorted_combined_nomenclature$slots %>% paste(collapse = " "), sep = "") %>% 
            trimws
        
    } else if (flag_is_intergenic == TRUE) {
        
        if (isoform_strand == "-") {
            intergenic_slot <- paste(isoform_strand, isoform_chr, ":", isoform_end, sep = "")
        } else {
            intergenic_slot <- paste(isoform_strand, isoform_chr, ":", isoform_start, sep = "")
        }
        
        final_VSR_nomenclature <- paste("Intergenic ", intergenic_slot, " ", list_first_pass_naming$list_body_exon_nomenclature %>% paste(collapse = " "), sep = "") %>% 
            trimws
        
    }
    
    return(final_VSR_nomenclature)
    
}

# END FLI_organise_matching() ###

# VSRs: FUNCTION TO ORGANISE BODY EXON REMATCHING TO SECONDARY SPLICE VARIANTS
## after the VSR is named, it is set. Nothing changes. 
## but because some body exons missing from the select variant may in fact have an equivalent in another splice variant, we will have to systematically match every body exon according to their lowest overlapped variant, and repeatedly re-call the VSR_select_reference_transcript_variant() function.
## NOTE: we are only re-matching body exons which were apparently intronic in the selected transcript.
## Behaviour: Select lowest variant for all LISs -> rename using the common lowest variant -> for the remaining un-named intronic exons, rename if they overlapped with a separate variant.
## inputs: VSR coords (pre-checked) and a list of start/end info for each LIS. 
## mode: VSR or LIS
VSR_LIS_organise_exon_naming <- function(VSR_coordinates, list_tibble_exon_start_end_per_LIS, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1, mode = NULL) {
    
    # DEBUG ###
    
    # list_first_pass_naming <- VSR_select_reference_transcript_variant(VSR_coordinates = "1:15152785-15214893:+", tibble_VSR_exon_start_end = tribble(
    #     ~start, ~end,
    #     15210489, 15210562
    #     
    # ), left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1) %>% VSR_name_body_exons_and_VSR(list_output_from_VSR_select_reference_transcript_variant = ., left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1)
    
    # match
    # VSR_coordinates <- "16:89737976-89740903:+"
    # no match
    # VSR_coordinates <- "16:89700400-89740923:+"
    
    # list_tibble_exon_start_end_per_LIS <- list(
    #     "LIS_1" = tribble(
    #         ~start, ~end,
    #         89740100, 89740803
    #     ),
    #     "LIS_2" = tribble(
    #         ~start, ~end,
    #         89739554, 89739993,
    #         89738975, 89739477,
    #         89738709, 89738881
    #     ),
    #     "LIS_3" = tribble(
    #         ~start, ~end,
    #         89720400, 89720582,
    #         89737976, 89738589,
    #         89739267, 89739993
    #     ),
    #     "LIS_4" = tribble(
    #         ~start, ~end,
    #         89739290, 89739477,
    #         89740100, 89740803
    #     )
    # ) 
    
    # VSR_coordinates <- automator_input_alternative_event_region
    # NOTE: we will deliberately choose to exclude rev() sequences which means that it won't mean anything to input reversed start/end coords.
    # This is because we decided that this VSR matching program has to automate strand matching (often users aren't given strand info for input). Also it's because we can fix the user input + start/end has meaning.
    # also, most DS tools output the LIS coords out of strand order. Therefore, we actually can't assume the user actually knows the sequence of matches. Have to consult the matched strand for that.
    # So that means we have to report the most common strand - **all elements which are on the oppoosite strand are made to be OS(), meaning the reverse complement.**
    # ALSO: We also require the user to know the exact exon connectivity of each LIS. This allows us to be able to give a more general description of heterogeneous LISs, such as those derived from multiple transcripts - ordering by exon number is meaningless between transcripts.
    # Also, all LIS coords must be in between the VSR coords. or it wont work. Triage will take care of this.
    # list_tibble_exon_start_end_per_LIS <- list_of_exon_start_end_tibbles
    # tibble_gtf_table <- tibble_ref_gtf

    # left_query_shift <- 0
    # right_query_shift <- 0
    # left_tolerance <- 1
    # right_tolerance <- 1
    
    # left_query_shift <- left_query_shift %>% paste %>% type.convert
    # right_query_shift <- right_query_shift %>% paste %>% type.convert
    # left_tolerance <- left_tolerance %>% paste %>% type.convert
    # right_tolerance <- right_tolerance %>% paste %>% type.convert
    # 
    # print(left_query_shift)
    # print(right_query_shift)
    # print(left_tolerance)
    # print(right_tolerance)
    
    ###########
    
    # demultiplex the VSR coordinates - currently we only use the chr and strand info from the VSR.
    query_chr <- gsub(x = VSR_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
    query_VSR_start <- gsub(x = VSR_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2") %>% type.convert
    query_VSR_end <- gsub(x = VSR_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3") %>% type.convert
    query_strand <- gsub(x = VSR_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
    
    # detect A3SS/A5SS - they are special cases because by definition, LISs which are flush with one side of the VSR are A3/5SS events.
    ## if there are any alternative exons which have a co-ordinate in common with the VSR, then it's BY DEEFINITION a A3SS/A5SS exon.
    ## modify the effective VSR ends if there were boundary A3/5SS exon extensions.
    
    ## in addition, if the VSR == alternative exon, then we assume it's IR.
    list_tibble_exon_start_end_per_LIS_flagged_extensions <- purrr::map(
        .x = list_tibble_exon_start_end_per_LIS,
        .f = ~.x %>% 
            dplyr::mutate("left_end_of_VSR" = .x$start + left_tolerance >= query_VSR_start & .x$start - left_tolerance <= query_VSR_start,
                          "right_end_of_VSR" = .x$end + left_tolerance >= query_VSR_end & .x$end - left_tolerance <= query_VSR_end)
    )
    
    list_tibble_exon_start_end_per_LIS_A35SS_corrected <- purrr::map(
        .x = list_tibble_exon_start_end_per_LIS_flagged_extensions,
        .f = function(a1) {
            
            # DEBUG ###
            # a1 <- list_tibble_exon_start_end_per_LIS_flagged_extensions[[1]]
            ###########
            
            output_tibble <- a1 %>% dplyr::mutate(
                "effective_VSR_start" = query_VSR_start,
                "effective_VSR_end" = query_VSR_end,
                "logical_is_IR" = FALSE
            )
            
            if (any(output_tibble$left_end_of_VSR == TRUE & output_tibble$right_end_of_VSR == FALSE)) {
                output_tibble$effective_VSR_start <- output_tibble[output_tibble$left_end_of_VSR == TRUE & output_tibble$right_end_of_VSR == FALSE, "end"] %>% unlist + 1
            }
            
            if (any(output_tibble$left_end_of_VSR == FALSE & output_tibble$right_end_of_VSR == TRUE)) {
                output_tibble$effective_VSR_end <- output_tibble[output_tibble$left_end_of_VSR == FALSE & output_tibble$right_end_of_VSR == TRUE, "start"] %>% unlist - 1
            }
            
            output_tibble[output_tibble$left_end_of_VSR == TRUE & output_tibble$right_end_of_VSR == TRUE, "logical_is_IR"] <- TRUE
            
            output_tibble[output_tibble$logical_is_IR == TRUE, "left_end_of_VSR"] <- FALSE
            output_tibble[output_tibble$logical_is_IR == TRUE, "right_end_of_VSR"] <- FALSE
            
            return(output_tibble)
            
        } )
    
    # magnetise all coords by matching query ends to ref vertices
    # in the process, count the number of common vertices for each ref transcript. this will be the main selected variant.
    # starts: magnetise to the nearest vertex in range and convert to exonic/intronic start. ends: magnetise to the nearest vertex in range and convert to exonic/intronic end.
    
    ## MAGNETISE QUERY VSR START
    list_magnetised_ref_start_to_query_VSR_start <- magnetise_genome_position_to_ref_starts(query_chr = query_chr, query_coord = query_VSR_start, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, query_shift = left_query_shift, query_tolerance = left_tolerance, ref_start_shift = 0, return_type = "exon")
    
    list_magnetised_ref_end_to_query_VSR_start <- magnetise_genome_position_to_ref_end(query_chr = query_chr, query_coord = query_VSR_start, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, query_shift = left_query_shift, query_tolerance = left_tolerance, ref_end_shift = 1, return_type = "exon")
    
    query_VSR_start_magnetised <- c(list_magnetised_ref_start_to_query_VSR_start$magnetised_coord, list_magnetised_ref_end_to_query_VSR_start$magnetised_coord) %>% na.omit %>% .[abs(. - query_VSR_start) == min(abs(. - query_VSR_start))] %>% .[1]
    
    tibble_ref_entries_containing_magnetised_query_VSR_start <- dplyr::bind_rows(list_magnetised_ref_start_to_query_VSR_start$tibble_ref_starts_matched_to_query_coord %>% .[which(.$start == query_VSR_start_magnetised), ], list_magnetised_ref_end_to_query_VSR_start$tibble_ref_ends_matched_to_query_coord %>% .[which((.$end + 1) == query_VSR_start_magnetised), ])
    
    if (length(query_VSR_start_magnetised) == 0 | is.na(query_VSR_start_magnetised) == TRUE) {
        query_VSR_start_magnetised <- query_VSR_start
    }
    
    ## MAGNETISE QUERY VSR END
    list_magnetised_ref_start_to_query_VSR_end <- magnetise_genome_position_to_ref_starts(query_chr = query_chr, query_coord = query_VSR_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, query_shift = left_query_shift, query_tolerance = left_tolerance, ref_start_shift = -1, return_type = "exon")
    
    list_magnetised_ref_end_to_query_VSR_end <- magnetise_genome_position_to_ref_end(query_chr = query_chr, query_coord = query_VSR_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, query_shift = left_query_shift, query_tolerance = left_tolerance, ref_end_shift = 0, return_type = "exon")
    
    query_VSR_end_magnetised <- c(list_magnetised_ref_start_to_query_VSR_end$magnetised_coord, list_magnetised_ref_end_to_query_VSR_end$magnetised_coord) %>% na.omit %>% .[abs(. - query_VSR_end) == min(abs(. - query_VSR_end))] %>% .[1]
    
    tibble_ref_entries_containing_magnetised_query_VSR_end <- dplyr::bind_rows(list_magnetised_ref_start_to_query_VSR_end$tibble_ref_starts_matched_to_query_coord %>% .[which((.$start - 1) == query_VSR_end_magnetised), ], list_magnetised_ref_end_to_query_VSR_end$tibble_ref_ends_matched_to_query_coord %>% .[which(.$end == query_VSR_end_magnetised), ])
    
    if (length(query_VSR_end_magnetised) == 0 | is.na(query_VSR_end_magnetised) == TRUE) {
        query_VSR_end_magnetised <- query_VSR_end
    }
    
    vector_VSR_matched_hgnc_variant_names <- c(tibble_ref_entries_containing_magnetised_query_VSR_start$hgnc_stable_variant_ID, tibble_ref_entries_containing_magnetised_query_VSR_end$hgnc_stable_variant_ID) %>% na.omit %>% unique
    
    tibble_global_VSR_possible_names <- name_a_single_junction(query_chr = query_chr, query_start = query_VSR_start_magnetised, query_end = query_VSR_end_magnetised, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, return_all_possibilities = TRUE, premagnetised = TRUE, left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1)
    
    ## LIS starts, ends and effective VSR starts and ends
    list_tibble_exon_start_end_per_LIS_magnetised <- purrr::map(
        .x = list_tibble_exon_start_end_per_LIS_A35SS_corrected,
        .f = function(a1) {
            
            # DEBUG ###
            # a1 <- list_tibble_exon_start_end_per_LIS_A35SS_corrected[[1]]
            ###########
            
            purrr::map(
                .x = a1 %>% rowwise() %>% dplyr::group_split(),
                .f = function(b1) {
                    
                    # DEBUG ###
                    # b1 <- a1 %>% rowwise() %>% dplyr::group_split() %>% .[[1]]
                    ###########
                    
                    query_start <- b1$start
                    query_end <- b1$end
                    
                    effective_VSR_start <- b1$effective_VSR_start
                    effective_VSR_end <- b1$effective_VSR_end
                    
                    ## MAGNETISE QUERY START COORD
                    list_magnetised_ref_start_to_query_start <- magnetise_genome_position_to_ref_starts(query_chr = query_chr, query_coord = query_start, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, query_shift = left_query_shift, query_tolerance = left_tolerance, ref_start_shift = 0, return_type = "exon")
                    
                    list_magnetised_ref_end_to_query_start <- magnetise_genome_position_to_ref_end(query_chr = query_chr, query_coord = query_start, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, query_shift = left_query_shift, query_tolerance = left_tolerance, ref_end_shift = 1, return_type = "exon")
                    
                    query_start_magnetised <- c(list_magnetised_ref_start_to_query_start$magnetised_coord, list_magnetised_ref_end_to_query_start$magnetised_coord) %>% na.omit %>% .[abs(. - query_start) == min(abs(. - query_start))] %>% .[1]
                    
                    tibble_ref_entries_containing_magnetised_query_start <- dplyr::bind_rows(list_magnetised_ref_start_to_query_start$tibble_ref_starts_matched_to_query_coord %>% .[which(.$start == query_start_magnetised), ], list_magnetised_ref_end_to_query_start$tibble_ref_ends_matched_to_query_coord %>% .[which((.$end + 1) == query_start_magnetised), ])
                    
                    if (length(query_start_magnetised) == 0 | is.na(query_start_magnetised) == TRUE) {
                        query_start_magnetised <- query_start
                    }
                    
                    ## MAGNETISE QUERY END COORD
                    list_magnetised_ref_start_to_query_end <- magnetise_genome_position_to_ref_starts(query_chr = query_chr, query_coord = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, query_shift = left_query_shift, query_tolerance = left_tolerance, ref_start_shift = -1, return_type = "exon")
                    
                    list_magnetised_ref_end_to_query_end <- magnetise_genome_position_to_ref_end(query_chr = query_chr, query_coord = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, query_shift = left_query_shift, query_tolerance = left_tolerance, ref_end_shift = 0, return_type = "exon")
                    
                    query_end_magnetised <- c(list_magnetised_ref_start_to_query_end$magnetised_coord, list_magnetised_ref_end_to_query_end$magnetised_coord) %>% na.omit %>% .[abs(. - query_end) == min(abs(. - query_end))] %>% .[1]
                    
                    tibble_ref_entries_containing_magnetised_query_end <- dplyr::bind_rows(list_magnetised_ref_start_to_query_end$tibble_ref_starts_matched_to_query_coord %>% .[which((.$start - 1) == query_end_magnetised), ], list_magnetised_ref_end_to_query_end$tibble_ref_ends_matched_to_query_coord %>% .[which(.$end == query_end_magnetised), ])
                    
                    if (length(query_end_magnetised) == 0 | is.na(query_end_magnetised) == TRUE) {
                        query_end_magnetised <- query_end
                    }
                    
                    ## MAGNETISE EFFECTIVE VSR START
                    list_magnetised_ref_start_to_effective_VSR_start <- magnetise_genome_position_to_ref_starts(query_chr = query_chr, query_coord = effective_VSR_start, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, query_shift = left_query_shift, query_tolerance = left_tolerance, ref_start_shift = 0, return_type = "exon")
                    
                    list_magnetised_ref_end_to_effective_VSR_start <- magnetise_genome_position_to_ref_end(query_chr = query_chr, query_coord = effective_VSR_start, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, query_shift = left_query_shift, query_tolerance = left_tolerance, ref_end_shift = 1, return_type = "exon")
                    
                    effective_VSR_start_magnetised <- c(list_magnetised_ref_start_to_effective_VSR_start$magnetised_coord, list_magnetised_ref_end_to_effective_VSR_start$magnetised_coord) %>% na.omit %>% .[abs(. - effective_VSR_start) == min(abs(. - effective_VSR_start))] %>% .[1]
                    
                    tibble_ref_entries_containing_magnetised_effective_VSR_start <- dplyr::bind_rows(list_magnetised_ref_start_to_effective_VSR_start$tibble_ref_starts_matched_to_query_coord %>% .[which(.$start == effective_VSR_start_magnetised), ], list_magnetised_ref_end_to_effective_VSR_start$tibble_ref_ends_matched_to_query_coord %>% .[which((.$end + 1) == effective_VSR_start_magnetised), ])
                    
                    if (length(effective_VSR_start_magnetised) == 0 | is.na(effective_VSR_start_magnetised) == TRUE) {
                        effective_VSR_start_magnetised <- effective_VSR_start
                    }
                    
                    ## MAGNETISE EFFECTIVE VSR END
                    list_magnetised_ref_start_to_effective_VSR_end <- magnetise_genome_position_to_ref_starts(query_chr = query_chr, query_coord = effective_VSR_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, query_shift = left_query_shift, query_tolerance = left_tolerance, ref_start_shift = -1, return_type = "exon")
                    
                    list_magnetised_ref_end_to_effective_VSR_end <- magnetise_genome_position_to_ref_end(query_chr = query_chr, query_coord = effective_VSR_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, query_shift = left_query_shift, query_tolerance = left_tolerance, ref_end_shift = 0, return_type = "exon")
                    
                    effective_VSR_end_magnetised <- c(list_magnetised_ref_start_to_effective_VSR_end$magnetised_coord, list_magnetised_ref_end_to_effective_VSR_end$magnetised_coord) %>% na.omit %>% .[abs(. - effective_VSR_end) == min(abs(. - effective_VSR_end))] %>% .[1]
                    
                    tibble_ref_entries_containing_magnetised_effective_VSR_end <- dplyr::bind_rows(list_magnetised_ref_start_to_effective_VSR_end$tibble_ref_starts_matched_to_query_coord %>% .[which((.$start - 1) == effective_VSR_end_magnetised), ], list_magnetised_ref_end_to_effective_VSR_end$tibble_ref_ends_matched_to_query_coord %>% .[which(.$end == effective_VSR_end_magnetised), ])
                    
                    if (length(effective_VSR_end_magnetised) == 0 | is.na(effective_VSR_end_magnetised) == TRUE) {
                        effective_VSR_end_magnetised <- effective_VSR_end
                    }
                    
                    ## retrieve matched HGNC stable variant IDs
                    ### for each *LIS*, collect all the UNIQUE matched HGNC stable variant IDs.
                    ### if there are non-A3/5SS events that make up an LIS, we will also take into account the VSR-matched ref entries too.
                    
                    return(b1 %>% 
                               dplyr::mutate(
                                   "query_start_magnetised" = query_start_magnetised,
                                   "query_end_magnetised" = query_end_magnetised,
                                   "effective_VSR_start_magnetised" = effective_VSR_start_magnetised,
                                   "effective_VSR_end_magnetised" = effective_VSR_end_magnetised
                               ) %>% as.list %>%
                               purrr::splice(
                                   "tibble_ref_entries_containing_magnetised_query_start" = tibble_ref_entries_containing_magnetised_query_start,
                                   "tibble_ref_entries_containing_magnetised_query_end" = tibble_ref_entries_containing_magnetised_query_end,
                                   "tibble_ref_entries_containing_magnetised_effective_VSR_start" = tibble_ref_entries_containing_magnetised_effective_VSR_start,
                                   "tibble_ref_entries_containing_magnetised_effective_VSR_end" = tibble_ref_entries_containing_magnetised_effective_VSR_end
                               )
                    )
                    
                } ) %>% 
                return
            
        } )
    
    # retrieve the HGNC stable variant IDs matched to LIS vertices
    # list_tibble_exon_start_end_per_LIS_magnetised <- purrr::map(
    #     .x = list_tibble_exon_start_end_per_LIS_magnetised,
    #     .f = function(a1) {
    #         
    #         # DEBUG ###
    #         # a1 <- list_tibble_exon_start_end_per_LIS_magnetised[[1]]
    #         ###########
    #         
    #         a1 %>% 
    #             purrr::splice(
    #                 "vector_hgnc_stable_variant_ids_matched_to_LIS_vertices" = list(purrr::map(
    #                     .x = a1,
    #                     .f = function(b1) {
    #                         
    #                         # DEBUG ###
    #                         # b1 <- a1[[1]]
    #                         ###########
    #                         
    #                         if (b1$left_end_of_VSR == TRUE | b1$right_end_of_VSR == TRUE) {
    #                             return(b1[c("tibble_ref_entries_containing_magnetised_query_start", "tibble_ref_entries_containing_magnetised_query_end")] %>% purrr::map(~.x$hgnc_stable_variant_ID) %>% unlist %>% unique)
    #                         } else {
    #                             return(b1[c("tibble_ref_entries_containing_magnetised_query_start", "tibble_ref_entries_containing_magnetised_query_end", "tibble_ref_entries_containing_magnetised_effective_VSR_start", "tibble_ref_entries_containing_magnetised_effective_VSR_end")] %>% purrr::map(~.x$hgnc_stable_variant_ID) %>% unlist %>% unique)
    #                         }
    #                         
    #                         
    #                         
    #                     } ) %>% unlist %>% unique)
    #             ) %>%
    #             return
    #         
    #     } )
    
    ## get list of the LIS-matched HGNC stable variant IDs
    # list_LIS_matched_hgnc_stable_variant_IDs <- purrr::map(
    #     .x = list_tibble_exon_start_end_per_LIS_magnetised,
    #     .f = ~.x$vector_hgnc_stable_variant_ids_matched_to_LIS_vertices
    # )
    # 
    # tally up the HGNC stable variant IDs
    # tibble_tally_hgnc_stable_variant_ids_matched_to_VSR_LIS_vertices <- c(list_LIS_matched_hgnc_stable_variant_IDs %>% unlist, vector_VSR_matched_hgnc_variant_names) %>% 
    #     table %>%
    #     as_tibble %>%
    #     set_names(nm = c("hgnc_stable_variant_ID", "tally")) %>%
    #     dplyr::arrange(desc(tally))
    
    # roll with the most commonly matched transcript IDs for now
    # vector_vertex_matched_hgnc_stable_variant_IDs <- tibble_tally_hgnc_stable_variant_ids_matched_to_VSR_LIS_vertices[tibble_tally_hgnc_stable_variant_ids_matched_to_VSR_LIS_vertices$tally == max(tibble_tally_hgnc_stable_variant_ids_matched_to_VSR_LIS_vertices$tally), ] %>% .$hgnc_stable_variant_ID
    
    # NAME THE BODY EXONS ###
    ## loop through each LIS. in each LIS, loop through each exon and call `name_a_single_exon()`
    ## variant override with the hgnc_stable_variant_ID tally, if any matches present.
    list_LIS_exons_named <- purrr::map2(
        .x = list_tibble_exon_start_end_per_LIS_magnetised,
        .y = 1:length(list_tibble_exon_start_end_per_LIS_magnetised),
        .f = function(a1, a2) {
            
            # print(a2)
            
            # DEBUG ###
            # a1 <- list_tibble_exon_start_end_per_LIS_magnetised[[1]]
            ###########
            
            purrr::map2(
                .x = a1,
                .y = 1:length(a1),
                .f = function(b1, b2) {
                    
                    # print(b2)
                    
                    # DEBUG ###
                    # b1 <- a1[[2]]
                    ###########
                    
                    if (b1$left_end_of_VSR == TRUE | b1$right_end_of_VSR == TRUE) {
                        # find free end
                        free_A3_5SS_vertex <- setdiff(c(b1$query_start_magnetised, b1$query_end_magnetised), c(b1$effective_VSR_start, b1$effective_VSR_end, b1$effective_VSR_start_magnetised, b1$effective_VSR_end_magnetised) %>% .[duplicated(.)])
                    }
                    
                    b1 %>% 
                        purrr::splice(
                            "exonic_matches" = if (b1$left_end_of_VSR == FALSE & b1$right_end_of_VSR == FALSE) {
                                name_a_single_exon(query_chr = b1$chr, query_start = b1$query_start_magnetised, query_end = b1$query_end_magnetised, query_strand = b1$strand, tibble_gtf_table = tibble_gtf_table, return_all_possibilities = TRUE, premagnetised = TRUE, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance) %>% 
                                    tibble::add_column(
                                        "query_start" = b1$query_start_magnetised,
                                        "query_end" = b1$query_end_magnetised
                                    )
                            },
                            "effective_VSR_matches" = if (b1$left_end_of_VSR == TRUE | b1$right_end_of_VSR == TRUE) {
                                name_a_single_junction(query_chr = b1$chr, query_start = b1$effective_VSR_start_magnetised, query_end = b1$effective_VSR_end_magnetised, query_strand = query_strand, tibble_gtf_table, return_all_possibilities = TRUE, premagnetised = TRUE, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = tibble_gtf_table, right_tolerance = right_tolerance) %>% 
                                    tibble::add_column(
                                        "query_start" = b1$effective_VSR_start,
                                        "query_end" = b1$effective_VSR_end
                                    )
                            },
                            "ref_exons_matched_to_free_A3_5SS_vertex" = if (b1$left_end_of_VSR == TRUE | b1$right_end_of_VSR == TRUE) {
                                dplyr::bind_rows(
                                    magnetise_genome_position_to_ref_end(query_chr = b1$chr, query_coord = free_A3_5SS_vertex, query_strand = b1$strand, tibble_gtf_table = tibble_gtf_table, query_shift = 0, query_tolerance = 0, ref_end_shift = 0, return_type = "exon") %>% .$tibble_ref_ends_matched_to_query_coord,
                                    magnetise_genome_position_to_ref_starts(query_chr = b1$chr, query_coord = free_A3_5SS_vertex, query_strand = b1$strand, tibble_gtf_table = tibble_gtf_table, query_shift = 0, query_tolerance = 0, ref_start_shift = 0, return_type = "exon") %>% .$tibble_ref_ends_matched_to_query_coord
                                )
                            }
                        ) %>% return
                    
                } ) %>% return
            
        } )
    
    # NAME EACH LIS and combine with VSR data. ###
    list_LIS_VSR_named_records <- purrr::imap(
        .x = list_LIS_exons_named,
        .f = function(a1, a2) {
            
            print(a2)
            
            # DEBUG ###
            # a1 <- list_LIS_exons_named[[1]]
            ###########
            
            if ( a1[purrr::map(.x = a1, .f = ~is.null(.x$exonic_matches) == FALSE) %>% unlist] %>% length > 0 ) {
                
                # deal with exonic matches.
                tibble_LIS_exonic_matches <- purrr::map2(
                    .x = a1[purrr::map(.x = a1, .f = ~is.null(.x$exonic_matches) == FALSE) %>% unlist],
                    .y = 1:length(a1[purrr::map(.x = a1, .f = ~is.null(.x$exonic_matches) == FALSE) %>% unlist]),
                    .f = function(b1, b2) {
                        
                        # DEBUG ###
                        # b1 <- a1[[1]]
                        # b2 <- 1:length(a1) %>% .[[1]]
                        ###########
                        
                        # add the exon numbers within the LIS and rbind into long tibble
                        if ( is.null(b1$exonic_matches) == FALSE ) {
                            long_tibble_named_exons_in_LIS <- b1$exonic_matches %>% 
                                dplyr::mutate("number_in_LIS" = b2)
                        } else if ( is.null(b1$exonic_matches) == TRUE ) {
                            long_tibble_named_exons_in_LIS <- tibble()
                        }
                        
                        
                    } ) %>% dplyr::bind_rows() %>% 
                    dplyr::group_by(variant_ID_slot) %>%
                    dplyr::mutate("number_of_LIS_exons_matched_to_variant" = n(),
                                  "match_type" = "exonic")
                
            } else {
                tibble_LIS_exonic_matches <- tibble()
            }
            
            # deal with A3/5SS matches.
            # the A3/5SS junction is the effective VSR and everything splits up.
            # the effective VSR is the same for the ENTIRE LIS, because the A3/5SS event has already defined the bounds of the VSR.
            # if there are any A3/5SS events, then we match them together with the other exons in the LIS.
            # if not, then continue with the global VSR as usual.
            if ( any(purrr::map(.x = a1, .f = ~c(.x$left_end_of_VSR, .x$right_end_of_VSR)) %>% unlist == TRUE) ) {
                
                tibble_A35SS_effective_VSR_matches <- a1[[1]]$effective_VSR_matches %>%
                    dplyr::mutate("number_in_LIS" = if ( nrow(tibble_LIS_exonic_matches) > 0 ) {max(tibble_LIS_exonic_matches$number_in_LIS) + 1} else {1},
                                  "vector_vertex_differences" = `query_start_match_distance` + `query_end_match_distance`,
                                  "number_of_ref_elements_to_describe_exon" = 0,
                                  "match_type" = "effective_VSR")
                
                tibble_LIS_match_entries <- dplyr::bind_rows(tibble_LIS_exonic_matches, tibble_A35SS_effective_VSR_matches)
                
                # feed into global VSR if there are no A3/5SS events
            } else if ( all(purrr::map(.x = a1, .f = ~c(.x$left_end_of_VSR, .x$right_end_of_VSR)) %>% unlist == FALSE) ) {
                
                tibble_LIS_match_entries <- tibble_LIS_exonic_matches
                
            }
            
            # for each LIS match entry, subset by matched HGNC stable variant ID 
            list_whole_LIS_named_all_matched_exons_only_split_by_hgnc_stable_variant_ID <- tibble_LIS_match_entries %>% ungroup() %>% dplyr::group_split(variant_ID_slot) %>% set_names(nm = purrr::map(.x = ., .f = ~.x$variant_ID_slot %>% unique) %>% unlist)
            
            # add in the VSR matches for each HGNC stable variant ID. 
            # we only use a common HGNC stable variant ID for VSRs ONLY if there is at least one LIS which has all exons matched to the HGNC stable variant ID
            # we take the first match with the lowest delta
            list_VSR_possible_names_split_by_hgnc_stable_variant_ID <- tibble_global_VSR_possible_names %>% dplyr::group_split(variant_ID_slot) %>% set_names(nm = purrr::map(.x = ., .f = ~.x$variant_ID_slot %>% unique) %>% unlist)
            
            # commonise the lists 
            vector_hgnc_stable_variant_IDs_in_common <- intersect(list_whole_LIS_named_all_matched_exons_only_split_by_hgnc_stable_variant_ID %>% names, list_VSR_possible_names_split_by_hgnc_stable_variant_ID %>% names)
            
            if (length(vector_hgnc_stable_variant_IDs_in_common) > 0) {
                
                list_named_LIS_and_VSR_split_by_hgnc_stable_variant_ID <- purrr::map2(
                    .x = list_whole_LIS_named_all_matched_exons_only_split_by_hgnc_stable_variant_ID[vector_hgnc_stable_variant_IDs_in_common],
                    .y = list_VSR_possible_names_split_by_hgnc_stable_variant_ID[vector_hgnc_stable_variant_IDs_in_common],
                    .f = ~list(
                        "LIS" = .x,
                        "VSR" = .y)
                )
                
                # find the transcript variant with the lowest delta
                list_LIS_VSR_record <- list_named_LIS_and_VSR_split_by_hgnc_stable_variant_ID %>% 
                    (function(x) {
                        
                        # DEBUG ###
                        # x <- list_named_LIS_and_VSR_split_by_hgnc_stable_variant_ID
                        ###########
                        
                        # whole LIS matches take priority
                        vector_metric <- purrr::map(.x = x, .f = ~nrow(.x$LIS)) %>% unlist
                        list_output <- x[which(vector_metric == max(vector_metric))]
                        
                        # full matches take priority
                        vector_metric <- purrr::map(.x = list_output, .f = ~which(c(.x$LIS$flag_is_exact_match, .x$VSR$flag_is_exact_match) == "FULL") %>% length) %>% unlist
                        list_output <- list_output[which(vector_metric == max(vector_metric))]
                        
                        # followed by half matches
                        vector_metric <- purrr::map(.x = list_output, .f = ~which(c(.x$LIS$flag_is_exact_match, .x$VSR$flag_is_exact_match) == "HALF") %>% length) %>% unlist
                        list_output <- list_output[which(vector_metric == max(vector_metric))]
                        
                        # followed by deltas
                        vector_metric <- purrr::map(.x = list_output, .f = ~c(.x$LIS$vector_vertex_differences, .x$VSR$query_start_match_distance, .x$VSR$query_end_match_distance) %>% sum) %>% unlist
                        list_output <- list_output[which(vector_metric == min(vector_metric))]
                        
                        # followed by number of ref exons needed to describe
                        vector_metric <- purrr::map(.x = list_output, .f = ~c(.x$LIS$number_of_ref_elements_to_describe_exon) %>% sum) %>% unlist
                        list_output <- list_output[which(vector_metric == min(vector_metric))]
                        
                        # and finally lowest hgnc_stable+_variant_ID
                        list_output <- list_output[names(list_output) == (mixedsort(names(list_output)) %>% .[1])] %>% .[[1]]
                        
                        return(list_output)
                        
                    } )
                
            } else {
                list_LIS_VSR_record <- list()
            }
            
            # organise LIS matching without taking into account the global VSR. 
            # whole LIS matches take priority
            if ( length(list_whole_LIS_named_all_matched_exons_only_split_by_hgnc_stable_variant_ID) > 0 ) {
                
                # find the transcript variant with the lowest delta
                tibble_LIS_record <- list_whole_LIS_named_all_matched_exons_only_split_by_hgnc_stable_variant_ID %>% 
                    (function(x) {
                        
                        # DEBUG ###
                        # x <- list_whole_LIS_named_all_matched_exons_only_split_by_hgnc_stable_variant_ID
                        ###########
                        
                        # whole LIS matches take priority
                        vector_metric <- purrr::map(.x = x, .f = ~nrow(.x)) %>% unlist
                        list_output <- x[which(vector_metric == max(vector_metric))]
                        
                        # full matches take priority
                        vector_metric <- purrr::map(.x = list_output, .f = ~which(c(.x$flag_is_exact_match) == "FULL") %>% length) %>% unlist
                        list_output <- list_output[which(vector_metric == max(vector_metric))]
                        
                        # followed by half matches
                        vector_metric <- purrr::map(.x = list_output, .f = ~which(c(.x$flag_is_exact_match) == "HALF") %>% length) %>% unlist
                        list_output <- list_output[which(vector_metric == max(vector_metric))]
                        
                        # followed by deltas
                        vector_metric <- purrr::map(.x = list_output, .f = ~.x$vector_vertex_differences %>% sum) %>% unlist
                        list_output <- list_output[which(vector_metric == min(vector_metric))]
                        
                        # followed by number of ref exons needed to describe
                        vector_metric <- purrr::map(.x = list_output, .f = ~.x$number_of_ref_elements_to_describe_exon %>% sum) %>% unlist
                        list_output <- list_output[which(vector_metric == min(vector_metric))]
                        
                        # and finally lowest hgnc_stable+_variant_ID
                        list_output <- list_output[names(list_output) == (mixedsort(names(list_output)) %>% .[1])]
                        
                        return(list_output)
                        
                    } ) %>% .[[1]]
                
            # if no ref matches for the entire LIS, then match each LIS exon individually 
            } else {
                # generate a tibble that records the HGNC stable variant ID, the exon entry and the LIS number.
                tibble_LIS_record <- tibble(
                    "variant_ID_slot" = character(),
                    "exon_slot" = character(),
                    "number_in_LIS" = numeric()
                )
                
            }
            
            if (length(tibble_LIS_record$number_in_LIS) < max(tibble_LIS_match_entries$number_in_LIS)) {
                
                tibble_LIS_record <- tibble_LIS_match_entries %>% 
                    dplyr::ungroup() %>%
                    dplyr::group_split(number_in_LIS) %>% 
                    purrr::map(
                        .f = ~.x %>%
                            # full matches take priority
                            .[length(which(.$flag_is_exact_match == "FULL")) == max(length(which(.$flag_is_exact_match == "FULL"))), ] %>%
                            # followed by half matches
                            .[length(which(.$flag_is_exact_match == "HALF")) == max(length(which(.$flag_is_exact_match == "HALF"))), ] %>%
                            # followed by deltas
                            .[.$vector_vertex_differences == min(.$vector_vertex_differences), ] %>%
                            # followed by number of ref exons needed to describe
                            .[.$number_of_ref_elements_to_describe_exon == min(.$number_of_ref_elements_to_describe_exon), ] %>%
                            # and finally lowest hgnc_stable+_variant_ID
                            .[.$variant_ID_slot == (.$variant_ID_slot %>% mixedsort %>% .[1]), ]
                    ) %>% dplyr::bind_rows() %>% dplyr::bind_rows(.[!.$number_in_LIS %in% .$number_in_LIS, ])
                
            }
            
            return(list(
                "list_LIS_VSR_record" = list_LIS_VSR_record,
                "tibble_LIS_record" = tibble_LIS_record
            ) )
            
        } )
    
    # finalise matching and names
    vector_segments_with_global_VSR_authority <- which(purrr::map(
        .x = list_LIS_VSR_named_records,
        .f = ~.x$list_LIS_VSR_record %>% length
    ) %>% unlist > 0)
    
    # if there is VSR authority, then we extract the LIS/VSR combinations and go thru the list once again.
    if (length(vector_segments_with_global_VSR_authority) > 0) {
        
        tibble_global_VSR_authority <- purrr::map(
            .x = list_LIS_VSR_named_records[vector_segments_with_global_VSR_authority],
            .f = ~tibble(
                "variant_ID_slot" = .x$list_LIS_VSR_record$VSR$variant_ID_slot,
                "full_match_count" = length(which(.x$list_LIS_VSR_record$VSR$flag_is_exact_match == "FULL")) + length(which(.x$list_LIS_VSR_record$LIS$flag_is_exact_match == "FULL")),
                "half_match_count" = length(which(.x$list_LIS_VSR_record$VSR$flag_is_exact_match == "HALF")) + length(which(.x$list_LIS_VSR_record$LIS$flag_is_exact_match == "HALF")),
                "delta_sum" = c(.x$list_LIS_VSR_record$VSR$query_start_match_distance, .x$list_LIS_VSR_record$VSR$query_end_match_distance, .x$list_LIS_VSR_record$LIS$vector_vertex_differences) %>% sum,
                "number_of_ref_elements_sum" = .x$list_LIS_VSR_record$LIS$number_of_ref_elements_to_describe_exon %>% sum
            )
        ) %>% dplyr::bind_rows()
        
        global_VSR_variant_ID_slot <- tibble_global_VSR_authority %>% 
            .[.$full_match_count == max(.$full_match_count), ] %>% 
            .[.$half_match_count == max(.$half_match_count), ] %>% 
            .[.$delta_sum == min(.$delta_sum), ] %>% 
            .[.$number_of_ref_elements_sum == min(.$number_of_ref_elements_sum), ] %>% 
            .[.$variant_ID_slot == (mixedsort(.$variant_ID_slot %>% .[1])), ] %>% 
            .$variant_ID_slot %>% .[1]
        
        # retrieve list indices of LISs which have VSR override as a result
        # these are the LISs which have been matched to the same hgnc_stable_variant_ID as the VSR variant ID
        logical_indices_LIS_with_global_VSR_override <- (1:length(list_LIS_VSR_named_records)) %in% (vector_segments_with_global_VSR_authority[purrr::map(.x = list_LIS_VSR_named_records[vector_segments_with_global_VSR_authority], .f = ~.x$list_LIS_VSR_record$VSR$variant_ID_slot == global_VSR_variant_ID_slot) %>% unlist])
        
        # retrieve VSR tibble from the first matching VSR table in the list
        tibble_global_VSR_final_naming <- list_LIS_VSR_named_records[logical_indices_LIS_with_global_VSR_override][[1]]$list_LIS_VSR_record$VSR
        
        if (tibble_global_VSR_final_naming$matched_strand == "+") {
            global_VSR_left_slot <- tibble_global_VSR_final_naming$exon_slot_query_start
            global_VSR_right_slot <- tibble_global_VSR_final_naming$exon_slot_query_end
        } else if (tibble_global_VSR_final_naming$matched_strand == "-") {
            global_VSR_left_slot <- tibble_global_VSR_final_naming$exon_slot_query_end
            global_VSR_right_slot <- tibble_global_VSR_final_naming$exon_slot_query_start
        }
        
    # if no VSR authority, then VSR and LISs will be independently named 
    } else if (length(vector_segments_with_global_VSR_authority) == 0) {
        
        logical_indices_LIS_with_global_VSR_override <- FALSE %>% rep(times = length(list_LIS_VSR_named_records))
        
        tibble_global_VSR_final_naming <- tibble_global_VSR_possible_names %>% 
            .[if (any(.$flag_is_exact_match == "FULL")) {.$flag_is_exact_match == "FULL"} else {TRUE}, ] %>% 
            .[if (any(.$flag_is_exact_match == "HALF")) {.$flag_is_exact_match == "HALF"} else {TRUE}, ] %>% 
            .[(.$query_start_match_distance + .$query_end_match_distance) == min(.$query_start_match_distance + .$query_end_match_distance), ] %>% 
            .[.$variant_ID_slot == (mixedsort(.$variant_ID_slot %>% .[1])), ] 
            
        global_VSR_variant_ID_slot <- tibble_global_VSR_final_naming$variant_ID_slot %>% .[1]
        
        if (tibble_global_VSR_final_naming$matched_strand == "+") {
            global_VSR_left_slot <- tibble_global_VSR_final_naming$exon_slot_query_start
            global_VSR_right_slot <- tibble_global_VSR_final_naming$exon_slot_query_end
        } else if (tibble_global_VSR_final_naming$matched_strand == "-") {
            global_VSR_left_slot <- tibble_global_VSR_final_naming$exon_slot_query_end
            global_VSR_right_slot <- tibble_global_VSR_final_naming$exon_slot_query_start
        }
        
    }
    
    # name the LIS
    ## deal with VSR override 
    ## VSR override: modify the list of LIS/VSR records by preserving only the LIS consistent with VSR match
    list_LIS_VSR_named_records[logical_indices_LIS_with_global_VSR_override] <- list_LIS_VSR_named_records[logical_indices_LIS_with_global_VSR_override] %>% purrr::map(~.x$list_LIS_VSR_record$LIS)
    ## no VSR override: modify by keeping only the independent LIS match
    list_LIS_VSR_named_records[!logical_indices_LIS_with_global_VSR_override] <- list_LIS_VSR_named_records[!logical_indices_LIS_with_global_VSR_override] %>% purrr::map(~.x$tibble_LIS_record)
    
    # find where A3/5SS events are. we will have to split them up and refactorise.
    logical_LIS_containing_A35SS <- purrr::map(.x = list_tibble_exon_start_end_per_LIS_A35SS_corrected, .f = ~any(c(.x$left_end_of_VSR, .x$right_end_of_VSR) == TRUE)) %>% unlist
    
    # finalise LIS naming
    ## IF VSR override present: don't explicitly specify the variant ID.
    ## NOTE: there will never be mixed variant IDs whenever there is VSR authority because we only considered VSR authority if there was an exact LIS match to the variant ID.
    
    list_final_LIS_name <- purrr::pmap(
        .l = list(
            "a1" = list_LIS_VSR_named_records,
            "a2" = logical_indices_LIS_with_global_VSR_override,
            "a3" = logical_LIS_containing_A35SS,
            "a4" = 1:length(logical_LIS_containing_A35SS)
        ),
        .f = function(a1, a2, a3, a4) {
        
            print(a4)
                
            # DEBUG ###
            # a1 <- list_LIS_VSR_named_records[[1]]
            # a2 <- logical_indices_LIS_with_global_VSR_override[[1]]
            # a3 <- logical_LIS_containing_A35SS[[1]]
            ###########
            
            tibble_exonic_records <- a1[a1$match_type == "exonic", ]
            tibble_effective_VSR_records <- a1[a1$match_type == "effective_VSR", ]
            
            # if global VSR override is active here, then we don't specify the effective VSR variant ID
            if (a2 == TRUE | a3 == FALSE) {
                effective_VSR_variant_ID_slot <- ""
            } else if (a2 == FALSE & a3 == TRUE) {
                effective_VSR_variant_ID_slot <- tibble_effective_VSR_records$variant_ID_slot
            }
            # do VSR slots
            if (a3 == FALSE) {
                effective_VSR_left_slot <- ""
                effective_VSR_right_slot <- ""
            } else if (tibble_effective_VSR_records$matched_strand == "+") {
                effective_VSR_left_slot <- tibble_effective_VSR_records$exon_slot_query_start
                effective_VSR_right_slot <- tibble_effective_VSR_records$exon_slot_query_end
            } else if (tibble_effective_VSR_records$matched_strand == "-") {
                effective_VSR_left_slot <- tibble_effective_VSR_records$exon_slot_query_end
                effective_VSR_right_slot <- tibble_effective_VSR_records$exon_slot_query_start
            }
            
            if ( nrow(tibble_exonic_records) > 0 ) {
                
                # detect majority strand and sort tibble by strand
                ## + strand
                if ( length(which(tibble_exonic_records$matched_strand == "+")) > length(which(tibble_exonic_records$matched_strand == "-")) ) {
                    flag_main_strand <- "+"
                    tibble_exonic_records <- tibble_exonic_records %>% dplyr::arrange(query_start)
                    ## - strand
                } else if ( length(which(tibble_exonic_records$matched_strand == "+")) < length(which(tibble_exonic_records$matched_strand == "-")) ) {
                    flag_main_strand <- "-"
                    tibble_exonic_records <- tibble_exonic_records %>% dplyr::arrange(desc(query_start))
                }
                
                tibble_exonic_records[tibble_exonic_records$matched_strand != flag_main_strand, "exon_slot"] <- paste("OS(", tibble_exonic_records[tibble_exonic_records$matched_strand != flag_main_strand, "exon_slot"], ")", sep = "")
                
                # test whether the HGNC stable variant ID is the same for the effective VSR and exonic slots. if they are, then we can factorise.
                ## we have to show the exonic variant ID no matter what if the effective VSR (A35SS) exists AND the variant matches are not consistent with the exonic match.
                ## also if there is no override.
                if ( (length(setdiff(c(tibble_exonic_records$variant_ID_slot, tibble_effective_VSR_records$variant_ID_slot), intersect(tibble_exonic_records$variant_ID_slot, tibble_effective_VSR_records$variant_ID_slot))) > 0 & nrow(tibble_effective_VSR_records) > 0) | (nrow(tibble_effective_VSR_records) == 0 & a2 == FALSE) ) {
                    exonic_slot <- paste("(", tibble_exonic_records$variant_ID_slot %>% unique, " ", paste(tibble_exonic_records$exon_slot, collapse = " "), ")", sep = "")
                    ## don't have to show the exonic variant only if either the effective VSR variant is consistent with exonic, or if there are no A3/5SS events and there is override.
                } else {
                    exonic_slot <- paste(tibble_exonic_records$exon_slot, collapse = " ")
                }
                
            } else {
                exonic_slot <- ""
            }
            
            # if no global VSR override AND there is A3/5SS event, then we bracket off the effective VSR variant ID.
            if (a2 == FALSE & nrow(tibble_effective_VSR_records) > 0) {
                segment_name <- paste("(", effective_VSR_variant_ID_slot, " ", effective_VSR_left_slot, " ", exonic_slot, " ", effective_VSR_right_slot, ")", sep = "")
            # global VSR override
            } else {
                segment_name <- paste(effective_VSR_variant_ID_slot, " ", effective_VSR_left_slot, " ", exonic_slot, " ", effective_VSR_right_slot, sep = "")
            }
            
            return(segment_name %>% trimws() %>% stringr::str_squish())
            
        } ) 
    
    # deal with A3/5SS events which will cause everything to split up
    # in those cases, we will need to add the global VSR in every LIS except the A3/5SS-containing ones
    if ( any(logical_LIS_containing_A35SS == TRUE) ) {
        
        list_final_LIS_name[!logical_LIS_containing_A35SS] <- purrr::map(
            .x = list_final_LIS_name[!logical_LIS_containing_A35SS],
            .f = ~paste("(", global_VSR_variant_ID_slot, " ", global_VSR_left_slot, " ", .x, " ", global_VSR_right_slot, ")", sep = "")
        )
        
        final_nomenclature <- list_final_LIS_name %>% paste(collapse = "/")
        
    # if no A3/5SS events at all, then simply apply the global VSR as usual   
    } else if ( all(logical_LIS_containing_A35SS == FALSE) ) {
        
        # we do the splitting up thing if we have just a few alternative events that have a mix of variant override and no override
        if ( any(logical_indices_LIS_with_global_VSR_override == FALSE) ) {
            
            list_final_LIS_name[logical_indices_LIS_with_global_VSR_override] <- purrr::map(
                .x = list_final_LIS_name[logical_indices_LIS_with_global_VSR_override],
                .f = ~paste("(", global_VSR_variant_ID_slot, " ", global_VSR_left_slot, " ", .x, " ", global_VSR_right_slot, ")", sep = "")
            )
            
            final_nomenclature <- list_final_LIS_name %>% paste(collapse = "/")
            
        } else {
            final_nomenclature <- paste(global_VSR_variant_ID_slot, " ", global_VSR_left_slot, " ", if ( length(list_final_LIS_name) > 1 ) { paste("(", paste(list_final_LIS_name, collapse = "/"), ")", sep = "") } else { paste(list_final_LIS_name, collapse = "/") }, " ", global_VSR_right_slot, sep = "")
        }
        
    }
    
    return(final_nomenclature)
    
}

# END VSR_LIS_organise_exon_naming() ###

# LSVs, AJ: FUNCTION TO GENERALLY MATCH AN UNLIMITED NUMBER OF JUNCTIONS IN A REGION OF VARIABLE SPLICING
## accepts a tibble of chr/start/end/strand of junctions, assumed from the same LSV.
## finds the transcript with most vertex matches and labels junctions with the caret notation from AS Code.
LSV_AJ_organise_junction_matching <- function(tibble_LSV_coords, tibble_gtf_table, single_junction = NULL, left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1) {
    
    # DEBUG ###
    # tibble_LSV_coords <- tibble(chr = c("1", "1", "1", "1", "1", "1"), start = c(150267478, 150267789, 150267991, 150268127, 150267478, 150267789), end = c(150268696, 150268696, 150268696, 150268696, 150267714, 150267955), strand = c("-", "-", "-", "-", "-", "-"))
    
    # left_query_shift <- 0
    # right_query_shift <- 0
    # left_tolerance <- 1
    # right_tolerance <- 1
    # tibble_gtf_table <- tibble_ref_gtf
    
    # tibble_LSV_coords <- tibble_LSV_chr_start_end_strand
    ###########
    
    left_query_shift <- left_query_shift %>% paste %>% type.convert
    right_query_shift <- right_query_shift %>% paste %>% type.convert
    left_tolerance <- left_tolerance %>% paste %>% type.convert
    right_tolerance <- right_tolerance %>% paste %>% type.convert
    
    # check the extents of the LSV
    LSV_chr <- tibble_LSV_coords$chr %>% unique
    LSV_start <- min(tibble_LSV_coords$start, tibble_LSV_coords$end)
    LSV_end <- max(tibble_LSV_coords$start, tibble_LSV_coords$end)
    LSV_strand <- tibble_LSV_coords$strand %>% unique
    
    # sort in increasing order of start co-ords
    tibble_LSV_coords <- tibble_LSV_coords %>% dplyr::arrange(start)
    
    # CHECK INTERGENIC
    tibble_trancripts_overlapping_query_LSV <- extract_overlapping_features(query_chr = LSV_chr, query_start = LSV_start - 1, query_end = LSV_end + 1, query_strand = LSV_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "transcript")
    
    flag_is_intergenic <- tibble_trancripts_overlapping_query_LSV %>% nrow == 0
    
    # SELECT REFERENCE VARIANT
    ## if not intergenic, pick the lowest variant ID with the most amount of common vertices. 
    ## since we're talking about LSVs here, the most common vertices should bias the count. 
    ## If intergenic, we're not going to even bother with naming. It would get way too big.
    if (flag_is_intergenic == FALSE) {
        
        # extract only the exonic entries from the overlapping transcripts 
        tibble_exonic_overlapped_transcripts_parent_entries <- tibble_gtf_table[which(tibble_gtf_table$transcript_id %in% (tibble_trancripts_overlapping_query_LSV$transcript_id %>% unique) & tibble_gtf_table$type == "exon"), ]
        # extract all the overlapping transcript vertices
        # vector_overlapping_transcript_vertices <- c(tibble_exonic_overlapped_transcripts_parent_entries$start, tibble_exonic_overlapped_transcripts_parent_entries$end) %>% unlist %>% unique
        
        # pool together all LSV vertices
        vector_LSV_starts_minus_one <- tibble_LSV_coords$start - 1
        vector_LSV_ends_plus_one <- tibble_LSV_coords$end + 1
        
        # find which exons have a common vertex with the LSV starts/ends
        list_of_tibbles_exons_with_common_end_with_LSV_starts <- purrr::map(.x = tibble_LSV_coords$start, .f = ~tibble_exonic_overlapped_transcripts_parent_entries[tibble_exonic_overlapped_transcripts_parent_entries$end < .x + left_query_shift + left_tolerance & tibble_exonic_overlapped_transcripts_parent_entries$end > .x - 2 + left_query_shift - left_tolerance, ])
        
        list_of_tibbles_exons_with_common_start_with_LSV_ends <- purrr::map(.x = tibble_LSV_coords$end, .f = ~tibble_exonic_overlapped_transcripts_parent_entries[tibble_exonic_overlapped_transcripts_parent_entries$start > .x + right_query_shift - right_tolerance & tibble_exonic_overlapped_transcripts_parent_entries$start < .x + 2 + right_query_shift + right_tolerance, ])
        
        # rbind and tibblise
        tibble_exons_vertex_matched <- dplyr::bind_rows(list_of_tibbles_exons_with_common_end_with_LSV_starts %>% rbindlist, list_of_tibbles_exons_with_common_start_with_LSV_ends %>% rbindlist) %>% as_tibble
        
        # in the event of no vertex matches at all, select the lowest variant_ID from the overlapped transcripts
        if (tibble_exons_vertex_matched %>% nrow == 0) {
            
            selected_transcript_id <- tibble_trancripts_overlapping_query_LSV$transcript_id %>% unique %>% mixedsort %>% .[1]
            
        } else {
            
            # tally up the number of times each transcript_id appeared
            suppressMessages(suppressWarnings( tibble_tally_transcript_ids_flanking_LSV_junctions <- tibble_exons_vertex_matched %>% dplyr::group_by(transcript_id) %>% dplyr::summarise("tally" = n()) %>% dplyr::arrange(desc(tally)) ))
            
            selected_transcript_id <- tibble_tally_transcript_ids_flanking_LSV_junctions[tibble_tally_transcript_ids_flanking_LSV_junctions$tally == max(tibble_tally_transcript_ids_flanking_LSV_junctions$tally), "transcript_id"] %>% unlist %>% mixedsort %>% .[1]
            
        }
        
        # MAIN SELECTED HGNC VARIANT NAME
        # retrieve stable HGNC variant ID and transcript version
        selected_hgnc_variant_name <- tibble_exonic_overlapped_transcripts_parent_entries[tibble_exonic_overlapped_transcripts_parent_entries$transcript_id == selected_transcript_id, "hgnc_stable_variant_ID"] %>% unlist %>% unique 
        
        # NAME THE CONSTITUENT JUNCTIONS
        # retrieve the gene name from selected transcript_id
        selected_gene_id <- tibble_exonic_overlapped_transcripts_parent_entries[tibble_exonic_overlapped_transcripts_parent_entries$transcript_id == selected_transcript_id, "gene_id"] %>% unlist %>% unique 
        
        # extract the parent exonic entries of the selected transcript
        tibble_parent_selected_transcript_exonic_entries <- tibble_exonic_overlapped_transcripts_parent_entries[tibble_exonic_overlapped_transcripts_parent_entries$transcript_id == selected_transcript_id & tibble_exonic_overlapped_transcripts_parent_entries$type == "exon", ]
        
        # update strand if it wasn't previously provided.
        if (!LSV_strand %in% c("+", "-")) {
            LSV_strand <- tibble_parent_selected_transcript_exonic_entries$strand %>% unique
        }
        
        # do the naming of junctions
        # loop thru each junction to find all the junction-flanking exons for our selected transcript name
        # if there's not exact flanking matches, then we look for any other vertices automatically.
        # automatically generate nomenclature blocks
        list_junction_nomenclature_slots <- purrr::map2(
            .x = tibble_LSV_coords$start, 
            .y = tibble_LSV_coords$end, 
            .f = function(a1, a2) {
                
                # DEBUG ###
                # a1 <- tibble_LSV_coords$start %>% .[3] + 5
                # a2 <- tibble_LSV_coords$end %>% .[3] - 23
                ###########
                
                # get exact flanking matches
                list_tibble_exact_flanking_matches <- extract_junction.flanking.exons(query_chr = LSV_chr, query_start = a1, query_end = a2, query_strand = "*", tibble_gtf_table = tibble_parent_selected_transcript_exonic_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, match_consecutive = FALSE, return_type = "exon")
                
                # if no exactly flanking matches, then consider matches in another transcript (but still overlapping the VSR)
                if (list_tibble_exact_flanking_matches %>% length == 1) {
                    
                    tibble_exact_flanking_matches <- list_tibble_exact_flanking_matches %>% .[[1]]
                    
                    # generate junction nomenclature block
                    return(paste("E", tibble_exact_flanking_matches$exon_number %>% min, "^E", tibble_exact_flanking_matches$exon_number %>% max, sep = ""))
                    
                } else if (list_tibble_exact_flanking_matches %>% length == 0) {
                    
                    list_flanking_matches_in_another_transcript <- extract_junction.flanking.exons(query_chr = LSV_chr, query_start = a1, query_end = a2, query_strand = "*", tibble_gtf_table = tibble_exonic_overlapped_transcripts_parent_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, match_consecutive = FALSE, match_inside_same_transcript = TRUE, return_type = "exon")
                    
                    if (list_flanking_matches_in_another_transcript %>% length > 0) {
                        
                        # choose the lowest transcript_id
                        tibble_flanking_matches_in_another_transcript <- list_flanking_matches_in_another_transcript[[names(list_flanking_matches_in_another_transcript) %>% mixedsort %>% .[1]]]
                        
                        # generate junction nomenclature block
                        if (tibble_flanking_matches_in_another_transcript$strand %>% unique == LSV_strand) {
                            
                            return(paste("(", tibble_flanking_matches_in_another_transcript$hgnc_stable_variant_ID %>% unique, " E", tibble_flanking_matches_in_another_transcript$exon_number %>% min, "^E", tibble_flanking_matches_in_another_transcript$exon_number %>% max, ")", sep = ""))
                            
                        } else {
                            
                            return(paste("(", tibble_flanking_matches_in_another_transcript$hgnc_stable_variant_ID %>% unique, " E", tibble_flanking_matches_in_another_transcript$exon_number %>% max, "^E", tibble_flanking_matches_in_another_transcript$exon_number %>% min, ")", sep = ""))
                            
                        }
                        
                        
                        # if the junction doesn't have a match in another transcript, then try higgedly piggeldy matching between different transcripts
                    } else if (list_flanking_matches_in_another_transcript %>% length == 0) {
                        
                        tibble_flanking_matches_between_different_transcripts <- extract_junction.flanking.exons(query_chr = LSV_chr, query_start = a1, query_end = a2, query_strand = "*", tibble_gtf_table = tibble_exonic_overlapped_transcripts_parent_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, match_consecutive = FALSE, match_inside_same_transcript = FALSE, return_type = "exon") %>% rbindlist %>% as_tibble
                        
                        # check which sides matched 
                        flag_junction_start_matched <- any(tibble_flanking_matches_between_different_transcripts$end < a1 + left_query_shift + left_tolerance & tibble_flanking_matches_between_different_transcripts$end > a1 - 2 + left_query_shift - left_tolerance)
                        flag_junction_end_matched <- any(tibble_flanking_matches_between_different_transcripts$start > a2 + left_query_shift - left_tolerance & tibble_flanking_matches_between_different_transcripts$start < a2 + 2 + left_query_shift + left_tolerance)
                        
                        # if both sides matched, then create nomenclature.
                        # use the lowest transcript id, using the same gene if possible
                        if (all(c(flag_junction_start_matched, flag_junction_end_matched) == TRUE)) {
                            
                            # JUNCTION START
                            tibble_junction_start_overlapping_exons <- tibble_flanking_matches_between_different_transcripts[tibble_flanking_matches_between_different_transcripts$end < a1 + left_query_shift + left_tolerance & tibble_flanking_matches_between_different_transcripts$end > a1 - 2 + left_query_shift - left_tolerance, ]
                            
                            if (selected_gene_id %in% tibble_junction_start_overlapping_exons$gene_id) {
                                junction_start_selected_transcript_id <- tibble_junction_start_overlapping_exons[tibble_junction_start_overlapping_exons$gene_id == selected_gene_id, ] %>% .[.$transcript_id == (.$transcript_id %>% mixedsort %>% .[1]), "transcript_id"] %>% unlist %>% unique
                            } else {
                                junction_start_selected_transcript_id <- tibble_junction_start_overlapping_exons[tibble_junction_start_overlapping_exons$transcript_id == (tibble_junction_start_overlapping_exons$transcript_id %>% mixedsort %>% .[1]), "transcript_id"] %>% unlist %>% unique
                            }
                            # extract exon entry
                            tibble_junction_start_selected_exon <- tibble_junction_start_overlapping_exons[tibble_junction_start_overlapping_exons$transcript_id == junction_start_selected_transcript_id, ]
                            
                            start_slot <- paste("(", tibble_junction_start_selected_exon$hgnc_stable_variant_ID, " E", tibble_junction_start_selected_exon$exon_number, ")", sep = "")
                            if (tibble_junction_start_selected_exon$strand != LSV_strand) {
                                start_slot <- paste(start_slot, "revcomp", sep = "")
                            }
                            
                            # JUNCTION END
                            tibble_junction_end_overlapping_exons <- tibble_flanking_matches_between_different_transcripts[tibble_flanking_matches_between_different_transcripts$start > a2 + left_query_shift - left_tolerance & tibble_flanking_matches_between_different_transcripts$start < a2 + 2 + left_query_shift + left_tolerance, ]
                            
                            if (selected_gene_id %in% tibble_junction_end_overlapping_exons$gene_id) {
                                junction_start_selected_transcript_id <- tibble_junction_end_overlapping_exons[tibble_junction_end_overlapping_exons$gene_id == selected_gene_id, ] %>% .[.$transcript_id == (.$transcript_id %>% mixedsort %>% .[1]), "transcript_id"] %>% unlist %>% unique
                            } else {
                                junction_start_selected_transcript_id <- tibble_junction_end_overlapping_exons[tibble_junction_end_overlapping_exons$transcript_id == (tibble_junction_end_overlapping_exons$transcript_id %>% mixedsort %>% .[1]), "transcript_id"] %>% unlist %>% unique
                            }
                            # extract exon entry
                            tibble_junction_end_selected_exon <- tibble_junction_end_overlapping_exons[tibble_junction_end_overlapping_exons$transcript_id == junction_start_selected_transcript_id, ]
                            
                            end_slot <- paste("(", tibble_junction_end_selected_exon$hgnc_stable_variant_ID, " E", tibble_junction_end_selected_exon$exon_number, ")", sep = "")
                            if (tibble_junction_end_selected_exon$strand != LSV_strand) {
                                end_slot <- paste(end_slot, "revcomp", sep = "")
                            }
                            
                            if (LSV_strand == "-") {
                                return(paste(end_slot, "^", start_slot, sep = ""))
                            } else {
                                return(paste(start_slot, "^", end_slot, sep = ""))
                            }
                            
                            # if all else fails then we will use the nearest exon extensions/truncations of our selected transcript ONLY.
                        } else {
                            
                            # pool all the selected transcript junctions
                            vector_selected_transcript_vertices <- tibble_parent_selected_transcript_exonic_entries[, c("start", "end")] %>% unlist %>% unique %>% sort
                            
                            # JUNCTION START
                            flag_junction_start_is_exonic <- extract_overlapping_features(query_chr = LSV_chr, query_start = a1, query_end = a1, query_strand = "*", tibble_gtf_table = tibble_parent_selected_transcript_exonic_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon") %>% nrow > 0
                            
                            closest_LHS_exon_number_to_junction_start <- tibble_parent_selected_transcript_exonic_entries[tibble_parent_selected_transcript_exonic_entries$start <= a1, ] %>% .[.$start == max(.$start), "exon_number"] %>% paste %>% type.convert
                            
                            left_exon_slot <- paste("E", closest_LHS_exon_number_to_junction_start, sep = "")
                            
                            closest_LHS_reference_vertex_to_junction_start <- max(vector_selected_transcript_vertices[vector_selected_transcript_vertices <= a1] %>% max,
                                                                                  tibble_parent_selected_transcript_exonic_entries$end %>% min
                            )
                            
                            closest_RHS_reference_vertex_to_junction_start <- min(vector_selected_transcript_vertices[vector_selected_transcript_vertices >= a1] %>% min,
                                                                                  tibble_parent_selected_transcript_exonic_entries$start %>% max
                            )
                            
                            if (flag_junction_start_is_exonic == TRUE) {
                                left_modifier_slot <- paste("–", closest_RHS_reference_vertex_to_junction_start - a1, sep = "")
                            } else if (flag_junction_start_is_exonic == FALSE) {
                                left_modifier_slot <- paste("+", a1 - closest_LHS_reference_vertex_to_junction_start, sep = "")
                            }
                            
                            # JUNCTION END
                            flag_junction_end_is_exonic <- extract_overlapping_features(query_chr = LSV_chr, query_start = a2, query_end = a2, query_strand = "*", tibble_gtf_table = tibble_parent_selected_transcript_exonic_entries, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon") %>% nrow > 0
                            
                            closest_RHS_exon_number_to_junction_end <- tibble_parent_selected_transcript_exonic_entries[tibble_parent_selected_transcript_exonic_entries$end >= a2, ] %>% .[.$end == min(.$end), "exon_number"] %>% paste %>% type.convert
                            
                            right_exon_slot <- paste("E", closest_RHS_exon_number_to_junction_end, sep = "")
                            
                            closest_LHS_reference_vertex_to_junction_end <- max(vector_selected_transcript_vertices[vector_selected_transcript_vertices <= a2] %>% max,
                                                                                tibble_parent_selected_transcript_exonic_entries$end %>% min
                            )
                            
                            closest_RHS_reference_vertex_to_junction_end <- min(vector_selected_transcript_vertices[vector_selected_transcript_vertices >= a2] %>% min,
                                                                                tibble_parent_selected_transcript_exonic_entries$start %>% max
                            )
                            
                            if (flag_junction_end_is_exonic == TRUE) {
                                right_modifier_slot <- paste("–", a2 - closest_LHS_reference_vertex_to_junction_end, sep = "")
                            } else if (flag_junction_end_is_exonic == FALSE) {
                                right_modifier_slot <- paste("+", closest_RHS_reference_vertex_to_junction_end - a2, sep = "")
                            }
                            
                            if (LSV_strand == "-") {
                                return(paste("(", right_exon_slot, right_modifier_slot, ")^(", left_modifier_slot, left_exon_slot, ")", sep = ""))
                            } else {
                                return(paste("(", left_exon_slot, left_modifier_slot, ")^(", right_modifier_slot, right_exon_slot, ")", sep = ""))
                            }
                            
                        }
                        
                    }
                    
                } 
                
            } )
        
        
        if (single_junction == TRUE) {
            
            if (LSV_strand == "-") {
                return(paste(selected_hgnc_variant_name, " ", list_junction_nomenclature_slots %>% unlist %>% rev %>% paste(collapse = " / "), sep = ""))
            } else {
                return(paste(selected_hgnc_variant_name, " ", list_junction_nomenclature_slots %>% unlist %>% paste(collapse = " / "), sep = ""))
            }
            
        } else if (single_junction == FALSE) {
            
            if (LSV_strand == "-") {
                return(paste(selected_hgnc_variant_name, " (", list_junction_nomenclature_slots %>% unlist %>% rev %>% paste(collapse = " / "), ")", sep = ""))
            } else {
                return(paste(selected_hgnc_variant_name, " (", list_junction_nomenclature_slots %>% unlist %>% paste(collapse = " / "), ")", sep = ""))
            }
            
        }
        
    } else if (flag_is_intergenic == TRUE) {
        
        return("Intergenic")
        
    }
    
}

# END LSV_AJ_organise_junction_matching() ###

# AE: FUNCTION TO NAME A SINGLE EXON
# we will take an exon's coordinates and find overlapping counterparts in the reference. Then add exon modifiers as required.
# return_all_possibilities: return not only the lowest HGNC stable variant ID, but all matches. THis is for doing LIS or VSR.
# premagnetised: if TRUE, then this will do magnetisation for edge tolerance cases. if FALSE, such as calls from VSR naming, then treat query starts and ends as is and don't use tolerance.
name_a_single_exon <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, return_all_possibilities = NULL, premagnetised = NULL, left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1) {
    
    # DEBUG ###
    # query_chr <- AE_query_chr
    # query_start <- AE_query_start %>% type.convert
    # query_end <- AE_query_end %>% type.convert
    # query_strand <- AE_query_strand
    # tibble_gtf_table <- tibble_ref_gtf
    # variant_ID_override <- vector_vertex_matched_hgnc_stable_variant_IDs
    # left_query_shift <- 0
    # right_query_shift <- 0
    # left_tolerance <- 1
    # right_tolerance <- 1
    
    # query_chr <- b1$chr
    # query_start <- b1$query_start_magnetised
    # query_end <- b1$query_end_magnetised
    # query_strand <- b1$strand
    # return_all_possibilities <- TRUE
    # premagnetised <- TRUE
    
    # print(query_chr, "\n")
    # print(query_start, "\n")
    # print(query_end, "\n")
    # print(query_strand, "\n")
    
    # data.class(left_query_shift)
    # print(right_query_end_shift, "\n")
    # print(left_match_tolerance, "\n")
    # print(right_match_tolerance, "\n")
    
    # print(tibble_gtf_table)
    ###########
    
    # left_query_shift <<- left_query_shift
    # right_query_shift <<- right_query_shift
    # left_tolerance <<- left_tolerance
    # right_tolerance <<- right_tolerance
    
    if (premagnetised == TRUE) {
        left_query_shift <- 0 
        right_query_shift <- 0
        left_tolerance <- 0
        right_tolerance <- 0
    }
    
    ## look for exact exon match in the reference
    # tibble_matched_reference_exons <- extract_matching.exons(query_chr = query_chr, query_start = query_start, query_end = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon")
    ## look for exact junction match in the reference
    # list_matched_reference_junctions <- extract_junction.flanking.exons(query_chr = query_chr, query_start = query_start, query_end = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, match_consecutive = FALSE, return_type = "exon")
    
    # look for reference matches overlapping the query
    tibble_overlapping_reference_transcripts <- extract_overlapping_features(query_chr = query_chr, query_start = query_start, query_end = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "transcript")
    
    # catch intergenic situation
    if (nrow(tibble_overlapping_reference_transcripts) == 0) {
        
        return(
            tibble(
                "variant_ID_slot" = "Intergenic",
                "exon_slot" = paste(query_chr, ":", query_start, "-", query_end, ":", query_strand, sep = ""),
                "flag_is_exact_match" = "NO",
                "intergenic" = TRUE,
                "vector_distance_between_query_and_overlapped_ref_exons" = NA
            )
        )
        
    } else if (nrow(tibble_overlapping_reference_transcripts) > 0) {
        
        # return ALL elements belonging to the overlapping reference transcripts. this will be the only table we use from now on.
        tibble_gtf_table <- tibble_gtf_table[which(tibble_gtf_table$transcript_id %in% tibble_overlapping_reference_transcripts$transcript_id), ]
        
        # check whether query exon overlaps with reference exons, as the transcript with exon with the highest overlap takes priority.
        tibble_overlapping_reference_exons <- extract_overlapping_features(query_chr = query_chr, query_start = query_start, query_end = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon")
        # calculate euclidean distances between the query exon and each overlapped ref exon.
        vector_distance_between_query_and_overlapped_ref_exons <- (abs(tibble_overlapping_reference_exons$start - query_start - left_query_shift) - left_tolerance) %>% (function(x) {x[x < 0] <- 0; return(x)} ) + (abs(tibble_overlapping_reference_exons$end - query_end - right_query_shift) - right_tolerance) %>% (function(x) {x[x < 0] <- 0; return(x)} )
        
        # if we don't use return_all_possibilities, then take the lowest ref. transcript that has the highest exonic overlap.
        if (return_all_possibilities == FALSE) {
            best_match_hgnc_variant_name <- tibble_overlapping_reference_exons[vector_distance_between_query_and_overlapped_ref_exons == min(vector_distance_between_query_and_overlapped_ref_exons), ] %>% .[.$hgnc_stable_variant_ID == (.$hgnc_stable_variant_ID %>% mixedsort %>% .[1]), ] %>% .$hgnc_stable_variant_ID
            # go back and select only the best matched overlapped reference transcript for feeding into the next steps
            tibble_overlapping_reference_transcripts <- tibble_overlapping_reference_transcripts[tibble_overlapping_reference_transcripts$hgnc_stable_variant_ID == best_match_hgnc_variant_name, ]
        }
        
        # loop thru each overlapped reference transcript and get the nomenclature
        list_tibble_parent_exons_introns_of_overlapped_ref_transcripts_split_by_hgnc_stable_variant_ID <- purrr::imap(
            .x = tibble_overlapping_reference_transcripts$hgnc_stable_variant_ID %>% mixedsort,
            .f = function(a1, a2) {
                
                print(a2)
                
                # DEBUG ###
                # a1 <- tibble_overlapping_reference_transcripts$hgnc_stable_variant_ID %>% mixedsort %>% .[[32]]
                ###########
                
                tibble_subset_ref_exons <- tibble_gtf_table[which(tibble_gtf_table$hgnc_stable_variant_ID == a1 & tibble_gtf_table$type == "exon"), ] %>% .[mixedorder(.$exon_number), ]
                
                ## add in intronic entries
                if (tibble_subset_ref_exons$strand %>% .[1] == "+") {
                    
                    for (i in (tibble_subset_ref_exons$exon_number %>% .[1:(length(.) - 1)] %>% rev)) {
                        
                        tibble_subset_ref_exons <- tibble_subset_ref_exons %>%
                            tibble::add_row(tibble_subset_ref_exons[1, ] %>% dplyr::mutate(
                                "start" = (tibble_subset_ref_exons[i, ] %>% .$end) + 1 ,
                                "end" = (tibble_subset_ref_exons[i + 1, ] %>% .$start) - 1,
                                "type" = "intron",
                                "exon_number" = i + 0.5
                            ) %>% 
                                dplyr::mutate(
                                    "width" = `end` - `start` + 1
                                ), .after = i )
                        
                    }
                    
                } else if (tibble_subset_ref_exons$strand %>% .[1] == "-") {
                    
                    for (i in (tibble_subset_ref_exons$exon_number %>% .[1:(length(.) - 1)] %>% rev)) {
                        
                        tibble_subset_ref_exons <- tibble_subset_ref_exons %>%
                            tibble::add_row(tibble_subset_ref_exons[1, ] %>% dplyr::mutate(
                                "start" = (tibble_subset_ref_exons[i + 1, ] %>% .$end) + 1 ,
                                "end" = (tibble_subset_ref_exons[i, ] %>% .$start) - 1,
                                "type" = "intron",
                                "exon_number" = i + 0.5
                            ) %>% 
                                dplyr::mutate(
                                    "width" = `end` - `start` + 1
                                ), .after = i )
                        
                    }
                    
                } # elif
                
                ## test query ends to see which exon/intron number the query spans
                ### test query start 
                tibble_entry_overlapping_query_start <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$start <= (query_start + left_query_shift + left_tolerance) & tibble_subset_ref_exons$end >= (query_start + left_query_shift - left_tolerance)), ]
                ### test query end 
                tibble_entry_overlapping_query_end <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$start <= (query_end + right_query_shift + right_tolerance) & tibble_subset_ref_exons$end >= (query_end + right_query_shift - right_tolerance)), ]
                
                # if left/right end fell within a transcript, measure whether it crosses more than 50% of the exon or junction
                # if less than 50%, then previous element extension.
                # if more than 50%, then current element truncation.
                
                # if the tolerance made any query cross an exon/intron boundary then we will select the leftmost START vertex or the rightmost END vertex.
                if (nrow(tibble_entry_overlapping_query_start) > 1) {
                    tibble_entry_overlapping_query_start <- tibble_entry_overlapping_query_start[tibble_entry_overlapping_query_start$start == min(tibble_entry_overlapping_query_start$start %>% .[. <= (query_start + left_query_shift + left_tolerance) & . >= (query_start + left_query_shift - left_tolerance)]), ]
                } # elif
                
                # if the tolerance made any query cross an exon/intron boundary then we will select the leftmost START vertex or the rightmost END vertex.
                if (nrow(tibble_entry_overlapping_query_end) > 1) {
                    tibble_entry_overlapping_query_end <- tibble_entry_overlapping_query_end[tibble_entry_overlapping_query_end$end == max(tibble_entry_overlapping_query_end$end %>% .[. <= (query_end + right_query_shift + right_tolerance) & . >= (query_end + right_query_shift - right_tolerance)]), ]
                } # elif
                
                # query lies within the confines of the transcript
                if (nrow(tibble_entry_overlapping_query_start) > 0 & nrow(tibble_entry_overlapping_query_end) > 0) {
                    
                    # QUERY START ###
                    logical_query_start_is_further_than_50_pct <- query_start <= 0.5*(tibble_entry_overlapping_query_start$start + tibble_entry_overlapping_query_start$end)
                    
                    # make special case for if query exon is confined within a single element of the transcript. then it will always "cross 50%" for truncation.
                    if (tibble_entry_overlapping_query_start$exon_number == tibble_entry_overlapping_query_end$exon_number) {
                        logical_query_start_is_further_than_50_pct <- TRUE
                    }
                    
                    # use the crossing logicals to determine the element bounds
                    if (logical_query_start_is_further_than_50_pct == FALSE) {
                        
                        # find the lowest start coord greater than the left query end
                        lowest_start_coord_greater_than_query_start <- tibble_subset_ref_exons$start %>% .[which(. >= query_start)] %>% min
                        
                        # find out the exon that has this coord
                        query_start_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$start == lowest_start_coord_greater_than_query_start), ] %>% .$exon_number
                        
                        query_start_distance_to_vertex <- lowest_start_coord_greater_than_query_start - query_start
                        
                        # get distance modifier
                        query_start_distance_modifier <- paste("+", query_start_distance_to_vertex, sep = "")
                        
                    } else if (logical_query_start_is_further_than_50_pct == TRUE) {
                        
                        # find the highest start coord smaller than the left query end
                        highest_start_coord_less_than_query_start <- tibble_subset_ref_exons$start %>% .[which(. <= query_start)] %>% max
                        
                        # find out the exon that has this coord
                        query_start_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$start == highest_start_coord_less_than_query_start), ] %>% .$exon_number
                        
                        query_start_distance_to_vertex <- query_start - highest_start_coord_less_than_query_start
                        
                        # get distance modifier
                        query_start_distance_modifier <- paste("–", query_start_distance_to_vertex, sep = "")
                        
                    }
                    
                    # QUERY END ###
                    logical_query_end_is_further_than_50_pct <- query_end >= 0.5*(tibble_entry_overlapping_query_end$start + tibble_entry_overlapping_query_end$end)
                    
                    # make special case for if query exon is confined within a single element of the transcript. then it will always "cross 50%" for truncation.
                    if (tibble_entry_overlapping_query_start$exon_number == tibble_entry_overlapping_query_end$exon_number) {
                        logical_query_end_is_further_than_50_pct <- TRUE
                    }
                    
                    # use the crossing logicals to determine the element bounds
                    if (logical_query_end_is_further_than_50_pct == FALSE) {
                        
                        # find the highest end coord less than the query end
                        highest_end_coord_less_than_query_end <- tibble_subset_ref_exons$end %>% .[which(. <= query_end)] %>% max
                        
                        # find out the exon that has this coord
                        query_end_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$end == highest_end_coord_less_than_query_end), ] %>% .$exon_number
                        
                        query_end_distance_to_vertex <- query_end - highest_end_coord_less_than_query_end
                        
                        # get distance modifier
                        query_end_distance_modifier <- paste("+", query_end_distance_to_vertex, sep = "")
                        
                    } else if (logical_query_end_is_further_than_50_pct == TRUE) {
                        
                        # find the lowest end coord larger than the left query end
                        lowest_end_coord_greater_than_query_end <- tibble_subset_ref_exons$end %>% .[which(. >= query_end)] %>% min
                        
                        # find out the exon that has this coord
                        query_end_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$end == lowest_end_coord_greater_than_query_end), ] %>% .$exon_number
                        
                        query_end_distance_to_vertex <- lowest_end_coord_greater_than_query_end - query_end
                        
                        # get distance modifier
                        query_end_distance_modifier <- paste("–", query_end_distance_to_vertex, sep = "")
                        
                    }
                    
                    # if the query start doesn't lie within the transcript, then it's an extension of the leftmost exon.    
                    # in that case, we have to be careful to make sure that the end's exon is always a truncation if it lies to the left of the leftmost end
                } else if (nrow(tibble_entry_overlapping_query_start) == 0 & nrow(tibble_entry_overlapping_query_end) > 0) {
                    
                    # QUERY END ###
                    logical_query_end_is_further_than_50_pct <- query_end >= 0.5*(tibble_entry_overlapping_query_end$start + tibble_entry_overlapping_query_end$end)
                    # test for leftmost end
                    if (length(which(tibble_subset_ref_exons$end < query_end)) == 0) {
                        logical_query_end_is_further_than_50_pct <- TRUE
                    }
                    
                    # use the crossing logicals to determine the element bounds
                    if (logical_query_end_is_further_than_50_pct == FALSE) {
                        
                        # find the highest end coord less than the query end
                        highest_end_coord_less_than_query_end <- tibble_subset_ref_exons$end %>% .[which(. <= query_end)] %>% max
                        
                        # find out the exon that has this coord
                        query_end_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$end == highest_end_coord_less_than_query_end), ] %>% .$exon_number
                        
                        query_end_distance_to_vertex <- query_end - highest_end_coord_less_than_query_end
                        
                        # get distance modifier
                        query_end_distance_modifier <- paste("+", query_end_distance_to_vertex, sep = "")
                        
                    } else if (logical_query_end_is_further_than_50_pct == TRUE) {
                        
                        # find the lowest end coord larger than the left query end
                        lowest_end_coord_greater_than_query_end <- tibble_subset_ref_exons$end %>% .[which(. >= query_end)] %>% min
                        
                        # find out the exon that has this coord
                        query_end_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$end == lowest_end_coord_greater_than_query_end), ] %>% .$exon_number
                        
                        query_end_distance_to_vertex <- lowest_end_coord_greater_than_query_end - query_end
                        
                        # get distance modifier
                        query_end_distance_modifier <- paste("–", query_end_distance_to_vertex, sep = "")
                        
                    }
                    
                    # if the query end doesn't lie within the transcript, then it's an extension of the rightmost exon.
                    # in that case, we have to be careful to make sure that the start's exon is always a truncation if it lies to the right of the rightmost start.
                } else if (nrow(tibble_entry_overlapping_query_start) > 0 & nrow(tibble_entry_overlapping_query_end) == 0) {
                    
                    # QUERY START ###
                    logical_query_start_is_further_than_50_pct <- query_start <= 0.5*(tibble_entry_overlapping_query_start$start + tibble_entry_overlapping_query_start$end)
                    # test for rightmost start
                    if (length(which(tibble_subset_ref_exons$start > query_start)) == 0) {
                        logical_query_start_is_further_than_50_pct <- TRUE
                    }
                    
                    # use the crossing logicals to determine the element bounds
                    if (logical_query_start_is_further_than_50_pct == FALSE) {
                        
                        # find the lowest start coord greater than the left query end
                        lowest_start_coord_greater_than_query_start <- tibble_subset_ref_exons$start %>% .[which(. >= query_start)] %>% min
                        
                        # find out the exon that has this coord
                        query_start_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$start == lowest_start_coord_greater_than_query_start), ] %>% .$exon_number
                        
                        query_start_distance_to_vertex <- lowest_start_coord_greater_than_query_start - query_start
                        
                        # get distance modifier
                        query_start_distance_modifier <- paste("+", query_start_distance_to_vertex, sep = "")
                        
                    } else if (logical_query_start_is_further_than_50_pct == TRUE) {
                        
                        # find the highest start coord smaller than the left query end
                        highest_start_coord_less_than_query_start <- tibble_subset_ref_exons$start %>% .[which(. <= query_start)] %>% max
                        
                        # find out the exon that has this coord
                        query_start_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$start == highest_start_coord_less_than_query_start), ] %>% .$exon_number
                        
                        query_start_distance_to_vertex <- query_start - highest_start_coord_less_than_query_start
                        
                        # get distance modifier
                        query_start_distance_modifier <- paste("–", query_start_distance_to_vertex, sep = "")
                        
                    }
                    
                } # elif
                
                if (nrow(tibble_entry_overlapping_query_start) == 0) {
                    
                    # QUERY START ###
                    # find the lowest end coord larger than the left query end
                    leftmost_exon_start <- tibble_subset_ref_exons$start %>% min
                    
                    query_start_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$start == leftmost_exon_start), ] %>% .$exon_number
                    
                    query_start_distance_to_vertex <- leftmost_exon_start - query_start
                    
                    # get distance modifier
                    query_start_distance_modifier <- paste("+", query_start_distance_to_vertex, sep = "")
                    
                }
                
                if (nrow(tibble_entry_overlapping_query_end) == 0) {
                    
                    # QUERY END ###
                    # find the lowest end coord larger than the left query end
                    rightmost_exon_end <- tibble_subset_ref_exons$end %>% max
                    
                    query_end_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$end == rightmost_exon_end), ] %>% .$exon_number
                    
                    query_end_distance_to_vertex <- query_end - rightmost_exon_end
                    
                    # get distance modifier
                    query_end_distance_modifier <- paste("+", query_end_distance_to_vertex, sep = "")
                    
                }
                
                list_exon_slot <- purrr::pmap(
                    .l = list(
                        "exon_intron_number" = seq(from = min(c(query_start_exon_number, query_end_exon_number)), to = max(c(query_start_exon_number, query_end_exon_number)), by = 0.5) %/% 1,
                        "logical_is_intron" = seq(from = min(c(query_start_exon_number, query_end_exon_number)), to = max(c(query_start_exon_number, query_end_exon_number)), by = 0.5) %% 1
                    ),
                    .f = function(exon_intron_number, logical_is_intron) {
                        
                        if (logical_is_intron == 0.0) {
                            exon_intron <- "E"
                        } else if (logical_is_intron == 0.5) {
                            exon_intron <- "I"
                        }
                        
                        return(paste(exon_intron, exon_intron_number, sep = ""))
                        
                    } )
                
                number_of_ref_elements_to_describe_exon <- length(list_exon_slot)
                
                if (length(list_exon_slot) > 3) {
                    list_exon_slot <- list(list_exon_slot[[1]], "++", list_exon_slot[[length(list_exon_slot)]])
                }
                
                if (query_start_distance_to_vertex == 0 & query_end_distance_to_vertex == 0) {
                    flag_is_exact_match <- "FULL"
                    exon_slot <- paste(list_exon_slot, collapse = " ")
                } else if (query_start_distance_to_vertex == 0 | query_end_distance_to_vertex == 0) {
                    flag_is_exact_match <- "HALF"
                    if (tibble_subset_ref_exons$strand %>% .[1] == "+") {
                        exon_slot <- paste(query_start_distance_modifier, paste(list_exon_slot, collapse = " "), query_end_distance_modifier, sep = "") %>% gsub(pattern = "\\+0", replacement = "") %>% gsub(pattern = "\\–0", replacement = "")
                    } else if (tibble_subset_ref_exons$strand %>% .[1] == "-") {
                        exon_slot <- paste(query_end_distance_modifier, paste(list_exon_slot, collapse = " "), query_start_distance_modifier, sep = "") %>% gsub(pattern = "\\+0", replacement = "") %>% gsub(pattern = "\\–0", replacement = "")
                    }
                } else {
                    flag_is_exact_match <- "NO"
                    if (tibble_subset_ref_exons$strand %>% .[1] == "+") {
                        exon_slot <- paste(query_start_distance_modifier, paste(list_exon_slot, collapse = " "), query_end_distance_modifier, sep = "")
                    } else if (tibble_subset_ref_exons$strand %>% .[1] == "-") {
                        exon_slot <- paste(query_end_distance_modifier, paste(list_exon_slot, collapse = " "), query_start_distance_modifier, sep = "")
                    }
                }
                
                exon_slot <- gsub(x = exon_slot, pattern = " \\+\\+ ", replacement = "++")
                
                return(
                    tibble(
                        "variant_ID_slot" = a1,
                        "exon_slot" = exon_slot,
                        "vector_vertex_differences" = (query_start_distance_to_vertex + query_end_distance_to_vertex),
                        "number_of_ref_elements_to_describe_exon" = number_of_ref_elements_to_describe_exon,
                        "matched_strand" = tibble_subset_ref_exons$strand %>% .[1],
                        "flag_is_exact_match" = flag_is_exact_match,
                        "intergenic" = FALSE
                    )
                )
                
            } ) # purrr::map
        
        list_tibble_parent_exons_introns_of_overlapped_ref_transcripts_split_by_hgnc_stable_variant_ID %>% 
            rbindlist %>% as_tibble %>% 
            return
        
    } # elif intergenic test
    
}

# END name_a_single_exon() ###

# AJ: FUNCTION TO NAME A SINGLE JUNCTIOn
# very similar to alternative exon naming, except we use a caret and also the extension/trunction rules are swapped around
# we will take a junction's coordinates and find overlapping counterparts in the reference. Then add modifiers as required. 
# we will not use intron symbols at this stage.
# return all possible matches for EACH junction end. This is because a junction can match two different transcripts at the same time
# return_all_possibilities: return not only the lowest HGNC stable variant ID, but all matches. THis is for doing LIS or VSR.
# premagnetised: if TRUE, then this will do magnetisation for edge tolerance cases. if FALSE, such as calls from VSR naming, then treat query starts and ends as is and don't use tolerance.
name_a_single_junction <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, return_all_possibilities = NULL, premagnetised = NULL, left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1) {
    
    # DEBUG ###
    # query_chr <- AE_query_chr
    # query_start <- AE_query_start %>% type.convert
    # query_end <- AE_query_end %>% type.convert
    # query_strand <- AE_query_strand
    # tibble_gtf_table <- tibble_ref_gtf
    # variant_ID_override <- vector_vertex_matched_hgnc_stable_variant_IDs
    # left_query_shift <- 0
    # right_query_shift <- 0
    # left_tolerance <- 1
    # right_tolerance <- 1
    
    # query_chr <- "11"
    # query_start <- 12864838
    # query_end <- 12879706
    # query_strand <- "*"
    # return_all_possibilities <- TRUE
    # premagnetised <- TRUE
    
    # print(query_chr, "\n")
    # print(query_start, "\n")
    # print(query_end, "\n")
    # print(query_strand, "\n")
    
    # data.class(left_query_shift)
    # print(right_query_end_shift, "\n")
    # print(left_match_tolerance, "\n")
    # print(right_match_tolerance, "\n")
    
    # print(tibble_gtf_table)
    ###########
    
    # left_query_shift <<- left_query_shift
    # right_query_shift <<- right_query_shift
    # left_tolerance <<- left_tolerance
    # right_tolerance <<- right_tolerance
    
    if (premagnetised == TRUE) {
        left_query_shift <- 0 
        right_query_shift <- 0
        left_tolerance <- 0
        right_tolerance <- 0
    }
    
    ## look for exact exon match in the reference
    # tibble_matched_reference_exons <- extract_matching.exons(query_chr = query_chr, query_start = query_start, query_end = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon")
    ## look for exact junction match in the reference
    # list_matched_reference_junctions <- extract_junction.flanking.exons(query_chr = query_chr, query_start = query_start, query_end = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, match_consecutive = FALSE, return_type = "exon")
    
    # look for reference matches overlapping the query
    tibble_overlapping_reference_transcripts <- extract_overlapping_features(query_chr = query_chr, query_start = query_start, query_end = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "transcript")
    
    # catch intergenic situation
    if (nrow(tibble_overlapping_reference_transcripts) == 0) {
        
        return(
            tibble(
                "variant_ID_slot" = "Intergenic",
                "exon_slot" = paste(query_chr, ":", query_start, "-", query_end, ":", query_strand, sep = ""),
                "flag_is_exact_match" = "NO",
                "intergenic" = TRUE,
                "vector_distance_between_query_and_overlapped_ref_exons" = NA
            )
        )
        
    } else if (nrow(tibble_overlapping_reference_transcripts) > 0) {
        
        # return ALL elements belonging to the overlapping reference transcripts. this will be the only table we use from now on.
        tibble_gtf_table <- tibble_gtf_table[which(tibble_gtf_table$transcript_id %in% tibble_overlapping_reference_transcripts$transcript_id), ]
        
        # check whether query exon overlaps with reference exons, as the transcript with exon with the highest overlap takes priority.
        # tibble_overlapping_reference_exons <- extract_overlapping_features(query_chr = query_chr, query_start = query_start, query_end = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon")
        # calculate euclidean distances between the query exon and each overlapped ref exon.
        # vector_distance_between_query_and_overlapped_ref_exons <- (abs(tibble_overlapping_reference_exons$start - query_start - left_query_shift) - left_tolerance) %>% (function(x) {x[x < 0] <- 0; return(x)} ) + (abs(tibble_overlapping_reference_exons$end - query_end - right_query_shift) - right_tolerance) %>% (function(x) {x[x < 0] <- 0; return(x)} )
        
        # if we don't use return_all_possibilities, then take the lowest ref. transcript that has the highest exonic overlap.
        # if (return_all_possibilities == FALSE) {
        #     best_match_hgnc_variant_name <- tibble_overlapping_reference_exons[vector_distance_between_query_and_overlapped_ref_exons == min(vector_distance_between_query_and_overlapped_ref_exons), ] %>% .[.$hgnc_stable_variant_ID == (.$hgnc_stable_variant_ID %>% mixedsort %>% .[1]), ] %>% .$hgnc_stable_variant_ID
        #     # go back and select only the best matched overlapped reference transcript for feeding into the next steps
        #     tibble_overlapping_reference_transcripts <- tibble_overlapping_reference_transcripts[tibble_overlapping_reference_transcripts$hgnc_stable_variant_ID == best_match_hgnc_variant_name, ]
        # }
        
        # loop thru each overlapped reference transcript and get the nomenclature
        list_tibble_parent_exons_introns_of_overlapped_ref_transcripts_split_by_hgnc_stable_variant_ID <- purrr::imap(
            .x = tibble_overlapping_reference_transcripts$hgnc_stable_variant_ID %>% mixedsort,
            .f = function(a1, a2) {
                
                print(a2)
                
                # DEBUG ###
                # a1 <- tibble_overlapping_reference_transcripts$hgnc_stable_variant_ID %>% mixedsort %>% .[[1]]
                ###########
                
                tibble_subset_ref_exons <- tibble_gtf_table[which(tibble_gtf_table$hgnc_stable_variant_ID == a1 & tibble_gtf_table$type == "exon"), ] %>% .[mixedorder(.$exon_number), ]
                
                ## add in intronic entries
                # if (tibble_subset_ref_exons$strand %>% .[1] == "+") {
                #     
                #     for (i in (tibble_subset_ref_exons$exon_number %>% .[1:(length(.) - 1)] %>% rev)) {
                #         
                #         tibble_subset_ref_exons <- tibble_subset_ref_exons %>%
                #             tibble::add_row(tibble_subset_ref_exons[1, ] %>% dplyr::mutate(
                #                 "start" = (tibble_subset_ref_exons[i, ] %>% .$end) + 1 ,
                #                 "end" = (tibble_subset_ref_exons[i + 1, ] %>% .$start) - 1,
                #                 "type" = "intron",
                #                 "exon_number" = i + 0.5
                #             ) %>% 
                #                 dplyr::mutate(
                #                     "width" = `end` - `start` + 1
                #                 ), .after = i )
                #         
                #     }
                #     
                # } else if (tibble_subset_ref_exons$strand %>% .[1] == "-") {
                #     
                #     for (i in (tibble_subset_ref_exons$exon_number %>% .[1:(length(.) - 1)] %>% rev)) {
                #         
                #         tibble_subset_ref_exons <- tibble_subset_ref_exons %>%
                #             tibble::add_row(tibble_subset_ref_exons[1, ] %>% dplyr::mutate(
                #                 "start" = (tibble_subset_ref_exons[i + 1, ] %>% .$end) + 1 ,
                #                 "end" = (tibble_subset_ref_exons[i, ] %>% .$start) - 1,
                #                 "type" = "intron",
                #                 "exon_number" = i + 0.5
                #             ) %>% 
                #                 dplyr::mutate(
                #                     "width" = `end` - `start` + 1
                #                 ), .after = i )
                #         
                #     }
                #     
                # } # elif
                
                # test query ends for the closest start/end coords, and where they lie relative to each other.
                # QUERY START
                closest_ref_vertex_to_the_left_of_query_start <- c(tibble_subset_ref_exons$start, tibble_subset_ref_exons$end) %>% .[. <= (query_start + left_query_shift + left_tolerance)] %>% max 
                closest_ref_vertex_to_the_right_of_query_start <- c(tibble_subset_ref_exons$start, tibble_subset_ref_exons$end) %>% .[. >= (query_start + left_query_shift - left_tolerance)] %>% min
                # if one distance is negative, then the vertex lies within the transcript. if positive, then it lies outside
                logical_query_start_is_in_transcript <- length((!is.infinite(c(closest_ref_vertex_to_the_left_of_query_start, closest_ref_vertex_to_the_right_of_query_start))) %>% which) == 2
                
                # deal with the condition that query start is outside the transcript
                if (logical_query_start_is_in_transcript == FALSE) {
                    
                    # find the lowest end coord larger than the left query end
                    leftmost_exon_start <- tibble_subset_ref_exons$start %>% min
                    
                    query_start_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$start == leftmost_exon_start), ] %>% .$exon_number
                    
                    query_start_distance_to_vertex <- leftmost_exon_start - query_start
                    
                    # get distance modifier
                    query_start_distance_modifier <- paste("+", query_start_distance_to_vertex, sep = "")
                    
                    if (tibble_subset_ref_exons$strand %>% .[1] == "+") {
                        exon_slot_query_start <- paste("5'(", query_start_distance_modifier,"E", query_start_exon_number, ")", sep = "")
                    } else if (tibble_subset_ref_exons$strand %>% .[1] == "-") {
                        exon_slot_query_start <- paste("(", "E", query_start_exon_number, query_start_distance_modifier, ")3'", sep = "")
                    }
                    
                    # if closest_ref_vertex_to_the_left... is a start, then we're in an exon (truncation). otherwise we're in an intron if it's an end (extension).
                    ## exonic
                } else if ((any(closest_ref_vertex_to_the_left_of_query_start %in% tibble_subset_ref_exons$start) & any(closest_ref_vertex_to_the_right_of_query_start %in% tibble_subset_ref_exons$end)) | (closest_ref_vertex_to_the_left_of_query_start == closest_ref_vertex_to_the_right_of_query_start)) {
                    
                    if (closest_ref_vertex_to_the_left_of_query_start == closest_ref_vertex_to_the_right_of_query_start) {
                        query_start_exon_number <- tibble_subset_ref_exons[which((tibble_subset_ref_exons$end == closest_ref_vertex_to_the_right_of_query_start) | (tibble_subset_ref_exons$start == closest_ref_vertex_to_the_right_of_query_start)), ] %>% .$exon_number
                    } else {
                        query_start_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$end == closest_ref_vertex_to_the_right_of_query_start), ] %>% .$exon_number
                    }
                    
                    query_start_distance_to_vertex <- (tibble_subset_ref_exons[tibble_subset_ref_exons$exon_number == query_start_exon_number, ] %>% .$end) - query_start + 1
                    
                    query_start_distance_modifier <- paste("–", query_start_distance_to_vertex, sep = "")
                    
                    if (tibble_subset_ref_exons$strand %>% .[1] == "+") {
                        exon_slot_query_start <- paste("E", query_start_exon_number, query_start_distance_modifier, sep = "")
                    } else if (tibble_subset_ref_exons$strand %>% .[1] == "-") {
                        exon_slot_query_start <- paste(query_start_distance_modifier, "E", query_start_exon_number, sep = "")
                    }
                    
                    ## intronic
                } else if (any(closest_ref_vertex_to_the_left_of_query_start %in% tibble_subset_ref_exons$end) & any(closest_ref_vertex_to_the_right_of_query_start %in% tibble_subset_ref_exons$start)) {
                    
                    query_start_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$end == closest_ref_vertex_to_the_left_of_query_start), ] %>% .$exon_number
                    
                    query_start_distance_to_vertex <- query_start - closest_ref_vertex_to_the_left_of_query_start - 1
                    
                    query_start_distance_modifier <- paste("+", query_start_distance_to_vertex, sep = "")
                    
                    if (tibble_subset_ref_exons$strand %>% .[1] == "+") {
                        exon_slot_query_start <- paste("E", query_start_exon_number, query_start_distance_modifier, sep = "")
                    } else if (tibble_subset_ref_exons$strand %>% .[1] == "-") {
                        exon_slot_query_start <- paste(query_start_distance_modifier, "E", query_start_exon_number, sep = "")
                    }
                    
                }
                
                # QUERY END
                closest_ref_vertex_to_the_left_of_query_end <- c(tibble_subset_ref_exons$start, tibble_subset_ref_exons$end) %>% .[. <= (query_end + right_query_shift + right_tolerance)] %>% max
                closest_ref_vertex_to_the_right_of_query_end <- c(tibble_subset_ref_exons$start, tibble_subset_ref_exons$end) %>% .[. >= (query_end + right_query_shift - right_tolerance)] %>% min
                # if one distance is negative, then the vertex lies within the transcript. if positive, then it lies outside
                logical_query_end_is_in_transcript <- length((!is.infinite(c(closest_ref_vertex_to_the_left_of_query_end, closest_ref_vertex_to_the_right_of_query_end))) %>% which) == 2
                
                # deal with the condition that query start is outside the transcript
                if (logical_query_end_is_in_transcript == FALSE) {
                    
                    # find the lowest end coord larger than the left query end
                    rightmost_exon_end <- tibble_subset_ref_exons$end %>% min
                    
                    query_end_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$end == rightmost_exon_end), ] %>% .$exon_number
                    
                    query_end_distance_to_vertex <- query_end - rightmost_exon_end
                    
                    # get distance modifier
                    query_end_distance_modifier <- paste("+", query_end_distance_to_vertex, sep = "")
                    
                    if (tibble_subset_ref_exons$strand %>% .[1] == "+") {
                        exon_slot_query_end <- paste("(", "E", query_end_exon_number, query_end_distance_modifier, ")3'", sep = "")
                    } else if (tibble_subset_ref_exons$strand %>% .[1] == "-") {
                        exon_slot_query_end <- paste("5'(", query_end_distance_modifier,"E", query_end_exon_number, ")", sep = "")
                    }
                    
                    # if closest_ref_vertex_to_the_left... is a start, then we're in an exon (truncation). otherwise we're in an intron if it's an end (extension).
                    ## exonic
                } else if ((any(closest_ref_vertex_to_the_left_of_query_end %in% tibble_subset_ref_exons$start) & any(closest_ref_vertex_to_the_right_of_query_end %in% tibble_subset_ref_exons$end)) | (closest_ref_vertex_to_the_left_of_query_end == closest_ref_vertex_to_the_right_of_query_end)) {
                    
                    if (closest_ref_vertex_to_the_left_of_query_end == closest_ref_vertex_to_the_right_of_query_end) {
                        query_end_exon_number <- tibble_subset_ref_exons[which((tibble_subset_ref_exons$end == closest_ref_vertex_to_the_left_of_query_end) | (tibble_subset_ref_exons$start == closest_ref_vertex_to_the_left_of_query_end)), ] %>% .$exon_number
                    } else {
                        query_end_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$start == closest_ref_vertex_to_the_left_of_query_end), ] %>% .$exon_number
                    }
                    
                    query_end_distance_to_vertex <- query_end - (tibble_subset_ref_exons[tibble_subset_ref_exons$exon_number == query_end_exon_number, ] %>% .$start) + 1
                    
                    query_end_distance_modifier <- paste("–", query_end_distance_to_vertex, sep = "")
                    
                    if (tibble_subset_ref_exons$strand %>% .[1] == "+") {
                        exon_slot_query_end <- paste(query_end_distance_modifier, "E", query_end_exon_number, sep = "")
                    } else if (tibble_subset_ref_exons$strand %>% .[1] == "-") {
                        exon_slot_query_end <- paste("E", query_end_exon_number, query_end_distance_modifier, sep = "")
                    }
                    
                    ## intronic
                } else if (any(closest_ref_vertex_to_the_left_of_query_end %in% tibble_subset_ref_exons$end) & any(closest_ref_vertex_to_the_right_of_query_end %in% tibble_subset_ref_exons$start)) {
                    
                    query_end_exon_number <- tibble_subset_ref_exons[which(tibble_subset_ref_exons$start == closest_ref_vertex_to_the_right_of_query_end), ] %>% .$exon_number
                    
                    query_end_distance_to_vertex <- closest_ref_vertex_to_the_right_of_query_end - query_end - 1
                    
                    query_end_distance_modifier <- paste("+", query_end_distance_to_vertex, sep = "")
                    
                    if (tibble_subset_ref_exons$strand %>% .[1] == "+") {
                        exon_slot_query_end <- paste(query_end_distance_modifier, "E", query_end_exon_number, sep = "")
                    } else if (tibble_subset_ref_exons$strand %>% .[1] == "-") {
                        exon_slot_query_end <- paste("E", query_end_exon_number, query_end_distance_modifier, sep = "")
                    }
                    
                }
                
                return(
                    tibble(
                        "variant_ID_slot" = a1,
                        "exon_slot_query_start" = exon_slot_query_start %>% gsub(pattern = "\\+0", replacement = "") %>% gsub(pattern = "\\–0", replacement = ""),
                        "exon_slot_query_end" = exon_slot_query_end %>% gsub(pattern = "\\+0", replacement = "") %>% gsub(pattern = "\\–0", replacement = ""),
                        "flag_is_exact_match" = if ((query_start_distance_to_vertex + query_end_distance_to_vertex) == 0) {"FULL"} else if (query_start_distance_to_vertex == 0 | query_end_distance_to_vertex == 0) {"HALF"} else {"NO"},
                        "query_start_match_distance" = query_start_distance_to_vertex,
                        "query_end_match_distance" = query_end_distance_to_vertex,
                        "matched_strand" = tibble_subset_ref_exons$strand %>% .[1],
                        "intergenic" = FALSE
                    )
                )
                
            } ) # purrr::map
        
        list_tibble_parent_exons_introns_of_overlapped_ref_transcripts_split_by_hgnc_stable_variant_ID %>% 
            rbindlist %>% as_tibble %>% 
            return
        
    } # elif intergenic test
    
}

# END name_a_single_junction() ###

# SHINY: FUNCTION TO SCREEN INPUT COORDINATES
## input: a vector of co-ordinates to be tested
triage_input_coordinates <- function(vector_input_coordinates, vector_of_expected_chromosomes, expect_stranded = NULL) {
    
    # DEBUG ###
    # vector_input_coordinates <- c(
    #     "16:89740100-89740803:8",
    #     "16:89738975-89739477:+",
    #     "16:89739267-89739993:*",
    #     "16:89739290-89739477:0",
    #     "16:89739554-89739993:.",
    #     "8:89738975-89739137",
    #     "89738975-89739132",
    #     "16:89738975-89739128:-",
    #     "16:89738709:89738881"
    # )
    
    # vector_input_coordinates <- c(
    #     "1:2-3:+"
    # )
    # expect_stranded <- TRUE
    # tibble_gtf_table <- tibble_ref_gtf
    
    ###########
    
    if (expect_stranded == TRUE) {
        
        # print(vector_input_coordinates)
        # print(tibble_gtf_table)
        
        # check for basic format
        if (grep(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)") %>% length != vector_input_coordinates %>% length) {
            warning("Error: Please check the format of your co-ordinates and make sure you have specified the strand.")
            return("triage fail")
        }
        
        vector_query_chr <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
        vector_query_VSR_start <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2")
        vector_query_VSR_end <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3")
        vector_query_strand <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
        
    } else if (expect_stranded == FALSE) {
        
        # check for basic format
        if (grep(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)") %>% length != vector_input_coordinates %>% length) {
            warning("Error: Please check the format of your co-ordinates and make sure you have specified the strand.")
            return("triage fail")
        }
        
        vector_query_chr <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)", replacement = "\\1")
        vector_query_VSR_start <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)", replacement = "\\2")
        vector_query_VSR_end <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)", replacement = "\\3")
        
    }
    
    # check chromosomes
    if (vector_query_chr %in% (vector_of_expected_chromosomes) %>% all == FALSE) {
        warning("Error: Please check the format of your chromosomes")
        return("triage fail")
    }
    
    # check start/end for numeric
    if (vector_query_VSR_start %>% type.convert %>% is.numeric == FALSE | vector_query_VSR_end %>% type.convert %>% is.numeric == FALSE) {
        warning("Error: Please check the format of your start/end co-ordinates")
        return("triage fail")
    }
    
    # check start/end for order. end must be > start
    if (any((vector_query_VSR_start %>% type.convert) > (vector_query_VSR_end %>% type.convert))) {
        warning("Error: Please ensure that all the start co-ordinates are not greater than your end co-ordinates")
        return("triage fail")
    }
    
    message("triage successful")
    
    return("triage successful")
    
}

# SHINY: ALTERNATIVE FUNCTION TO SCREEN INPUT COORDINATES
## based on IGV behaviour. 
parse_input_coordinates <- function(input_coordinates, vector_of_expected_chromosomes) {
    
    # DEBUG ###
    # vector_input_coordinates <- c(
    #     "16:89740100-89740803:8",
    #     "16:89738975-89739477:+",
    #     "16:89739267-89739993:*",
    #     "16:89739290-89739477:0",
    #     "16:89739554-89739993:.",
    #     "8:89738975-89739137",
    #     "89738975-89739132",
    #     "16:89738975-89739128:-",
    #     "16:89738709:89738881"
    # )
    
    # vector_input_coordinates <- c(
    #     "1:2-3:+"
    # )
    # expect_stranded <- TRUE
    # tibble_gtf_table <- tibble_ref_gtf
    
    ###########
    
    # split into vector of chr start end strand.
    vector_chr_start_end_strand <- input_coordinates %>% trimws() %>% gsub(pattern = "[^0-9]", replacement = " ") %>% stringr::str_squish() %>% strsplit(split = " ") %>% unlist %>% type.convert
    
    # if length is 2, then it's a single nt. do a +/- 50nt window
    if (vector_chr_start_end_strand %>% length == 2) {
        output <- c(vector_chr_start_end_strand[1], vector_chr_start_end_strand[2] - 50, vector_chr_start_end_strand[2] + 50)
    # if length is 3, then it's a geniune range.
    } else if (vector_chr_start_end_strand %>% length == 3) {
        output <- c(vector_chr_start_end_strand[1], min(vector_chr_start_end_strand[2], vector_chr_start_end_strand[3]), max(vector_chr_start_end_strand[2], vector_chr_start_end_strand[3]))
    } else if (vector_chr_start_end_strand %>% length == 4) {
        output <- c(vector_chr_start_end_strand[1], min(vector_chr_start_end_strand[2], vector_chr_start_end_strand[3]), max(vector_chr_start_end_strand[2], vector_chr_start_end_strand[3]))
    } else {
        output <- "non_coord"
    }
    
    return(output)
    
}


### SHINY ####

ui <- fluidPage(
    
    navbarPage(title = "EDN Suite",               
               tabPanel(icon("info"),
                        
                        fluidRow(
                            
                            column(width = 12, ######
                                   
                                   br(),
                                   p("This is a collection of 4 tools to help generate, visualise, interpret and compare splicing events using the Extended Delta Notation.",
                                     br(),
                                     br(),
                                     strong("EDN Automator: "), "Automatically generates a publication-ready shorthand based on genomic co-ordinates.",
                                     br(),
                                     strong("EDN Workshop: "), "From user-given ranges provided in genomic co-ordinates, draws a schematic of matches to reference transcripts. HGNC stable variant IDs and exon numbering is explicitly shown for each transcript and distances to the nearest exon vertices are automatically calculated. This is useful for contextualising a given co-ordinate range in terms of reference transcripts, as well as for manually building the shorthand.",
                                     br(),
                                     strong("EDN Reverse Translate: "), "From any user-given EDN, produces genomic co-ordinates, draws a schematic defined by the shorthand and the associated reference transcripts. It also checks whether a pair of shorthands are equivalent.", style = "text-align:justify; color:black; background-color:lavender;padding:15px; border-radius:10px",
                                     br(),
                                     strong("Library of EDN: "), "A collection of curated literature re-annotations using EDN. Any researcher is free to submit a re-annotation of a splicing event previously described in literature using EDN. Multiple pieces of evidence are required in the submission. Once submitted, articles are curated and checked for duplication, before being added to a public collection. In this way, splicing events in literature are all in one place and searchable using EDN with a single click. There is also a leaderboard for high scores - this ranks researchers according to the number of submissions they have made.", style = "text-align:justify; color:black; background-color:lavender;padding:15px; border-radius:10px"),
                                   br(),
                                   p("All reference data is retrieved from Ensembl and RefSeq.", style="text-align:justify; color:black; background-color:papayawhip; padding:15px ;border-radius:10px")
                            ),
                            
                            # save this for publication link shenanigans.
                            # column(
                            #     br(),
                            #     tags$img(src="Gobernacion.png",width="200px",height="130px"),
                            #     br(),
                            #     br(),
                            #     p("For more information please check the",em("Anuario Estadístico de Antioquia's"),"page clicking",
                            #       br(),
                            #       a(href="http://www.antioquiadatos.gov.co/index.php/anuario-estadistico-de-antioquia-2016", "Here",target="_blank"),style="text-align:center;color:black"),
                            #     
                            #     width = 6))
                        ) # fluidRow ######
                        
               ), # tabPanel
               
               tabPanel("EDN Automator",
                        
                        fluidRow(
                            
                            column(width = 5,
                                   
                                   titlePanel("EDN Automator"),
                                   
                                   selectInput("automator_structure_type", 
                                               label = "Select the type of splice structure to describe", 
                                               choices = list("Full-length isoform" = c("Upload assembled transcriptome (GTF)", "Manually enter exon co-ordinates"), 
                                                              "Alternative splicing" = c("Alternative exon", "Alternative junction", "Local Isoform Segment (LIS)"), 
                                                              "Alternatively spliced regions" = c("Variable Splice Region (VSR)", "Local Splice Variation (LSV) (junctions only)")),
                                               width = "300px"),
                                   
                                   selectInput("automator_genome_assembly", 
                                               label = "Select the genome assembly to use (Ensembl release)", 
                                               choices = list("GRCh38 (hg38)" = list.files("data") %>% gsub(pattern = "annotated_ensembl_gtf_release_(.*).txt", replacement = "\\1") %>% type.convert %>% max %>% as.character,
                                                              "GRCh37 (hg19)" = "75",
                                                              "NCBI36 (hg18)" = "54"
                                               ),
                                               width = "300px"
                                   ),
                                   
                                   actionButton("automator_import_GTF", "Import genome assembly", icon = icon("file-import"),
                                                class = "btn btn-primary", width = "300px"),
                                   
                                   br(),
                                   br(),
                                   
                                   # magnetisation options required
                                   div(style = "align: left; text-align: center; width: 300px;",  
                                       
                                       shinyWidgets::dropdownButton(
                                           
                                           textInput("automator_left_query_end_shift", label = "Left (5') query end shift (default = 0)", value = 0, placeholder = "e.g. 0"),
                                           
                                           textInput("automator_right_query_end_shift", label = "Right (3') query end shift (default = 0)", value = 0, placeholder = "e.g. 0"),
                                           
                                           textInput("automator_left_match_tolerance", label = "Left (5') match tolerance (default = 1)", value = 1, placeholder = "e.g. 0"),
                                           
                                           textInput("automator_right_match_tolerance", label = "Right (3') match tolerance (default = 1)", value = 1, placeholder = "e.g. 0"),
                                           
                                           circle = FALSE, status = "info", icon = icon("wrench"),
                                           label = "Tweak match tolerances",
                                           inline = FALSE
                                       )
                                       
                                   ),
                                   
                                   br(),
                                   
                                   # STRATEGY: 
                                   # FLI & LIS: vertex matching
                                   # VSRs and AEs: VSR and exon coord matching, splicemode autodetect
                                   # LSVs and AJs: Junction matching. Forward-slashes everywhere.
                                   
                                   # generate an upload box for users to upload full-length reconstructed isoforms
                                   conditionalPanel(
                                       condition = "input.automator_structure_type == 'Upload assembled transcriptome (GTF)'",
                                       fileInput("automator_path_full_length_recon_gtf", "Choose GTF file",
                                                 accept = c(".gtf"), 
                                                 width = "300px")
                                   ),
                                   
                                   # generate a text box for users to state the number of furrr multilayer cores
                                   conditionalPanel(
                                       condition = "input.automator_structure_type == 'Upload assembled transcriptome (GTF)'",
                                       textInput("automator_ncores", label = "Specify the number of parallel processing cores (default: 1 (serial)) (example: 8x4, 4)", value = 1, width = "300px")
                                   ),
                                   
                                   # generate a text box for users to state the number of exons in the full length isoform
                                   conditionalPanel(
                                       condition = "input.automator_structure_type == 'Manually enter exon co-ordinates'",
                                       textInput("automator_full_length_number_of_exons", label = "Enter the number of exons in the transcript", value = 1, width = "300px"),
                                   ),
                                   
                                   # generate a text box for users to enter the genome-relative coords of alternative exon
                                   conditionalPanel(
                                       condition = "input.automator_structure_type == 'Alternative exon'",
                                       textInput("automator_alternative_exon_coords", label = "Enter the genome-relative co-ordinates of the exon", placeholder = "e.g. 16:2756334-2756606", width = "300px")
                                   ),
                                   
                                   # generate a text box for users to state the number of exons in the full length isoform
                                   conditionalPanel(
                                       condition = "input.automator_structure_type == 'Alternative junction'",
                                       textInput("automator_alternative_junc_coords", label = "Enter the genome-relative co-ordinates of the junction", placeholder = "e.g. 16:2756607-2757471", width = "300px")
                                   ),
                                   
                                   # generate a text box for users to state the number of exons in the LIS
                                   conditionalPanel(
                                       condition = "input.automator_structure_type == 'Local Isoform Segment (LIS)'",
                                       textInput("automator_LIS_number_of_exons", label = "Enter the number of exons inside the LIS (not including the constitutive exons at the ends)", value = 1, width = "300px")
                                   ),
                                   
                                   # generate a text box for users to state the number of exons in the LIS
                                   conditionalPanel(
                                       condition = "input.automator_structure_type == 'Local Isoform Segment (LIS)' || input.automator_structure_type == 'Variable Splice Region (VSR)'",
                                       textInput("automator_alternative_event_region", label = "Enter the co-ordinates of the alternative event region (bounded by constitutive exons)", placeholder = "e.g. 16:2756607-2757471", width = "300px")
                                   ),
                                   
                                   # generate a text box for users to state the number of independent splicing events to be included in the VSR
                                   conditionalPanel(
                                       condition = "input.automator_structure_type == 'Variable Splice Region (VSR)' || input.automator_structure_type == 'Local Splice Variation (LSV) (junctions only)'",
                                       textInput("automator_alternative_region_number_of_independent_events", label = "State the number of independent events in this region (an independent region is an exon or LIS or one junction only)", value = 1, width = "300px")
                                   ),
                                   
                                   div(style = "padding-left: 50px; width: 300px;", 
                                       uiOutput('automator_reactive_UI_1')
                                   ),
                                   
                                   div(style = "padding-left: 50px; width: 300px;",
                                       uiOutput('automator_reactive_UI_2')
                                   ),
                                   
                            ),
                            
                            column(width = 7,
                                   
                                   div(style = "display:inline-block;
        box-shadow: 0px 0px 0px #888888;
        width:400px;
        height:0px;
        padding-top:50px;
        position:relative;",
                                       actionButton("automator_button_execute", "Generate shorthand", icon = icon("searchengin"),
                                                    class = "btn btn-primary", width = "300px")
                                   ),
                                   
                                   br(),
                                   
                                   tags$style("#automator_nomenclature_output {
    font-family:'Helvetica';
    font-size:15px;
               color:blue;
               display:block;
               width: 300px; 
               max-width: 200%;
               white-space: pre-wrap; }"),
                                   
                                   div(style = "text-align:center;
        box-shadow: 0px 0px 0px #888888;
        width:300px;
        height:200px;
        padding-top:40px;
        position:relative;",
                                       verbatimTextOutput("automator_nomenclature_output", placeholder = TRUE)
                                   ),
                                   
                                   br(),
                                   
                                   div(style = "display:inline-block;
        box-shadow: 0px 0px 0px #888888;
        width:400px;
        height:0px;
        padding-top:50px;
        position:relative;",
                                       # generate an upload box for users to upload full-length reconstructed isoforms
                                       conditionalPanel(
                                           condition = "input.automator_structure_type == 'Upload assembled transcriptome (GTF)'",
                                           downloadButton("automator_download_GTF_naming", "Download Results")
                                       ),
                                   ),
                                   
                                   div(style = "display:inline-block;
        box-shadow: 0px 0px 0px #888888;
        width:400px;
        height:0px;
        padding-top:50px;
        position:relative;",
                                       tableOutput("automator_named_GTF_entries")
                                   )
                                   
                            )
                            
                        )
                        
               ), # tabPanel
               
               tabPanel("EDN Workshop",
                        
                        fluidRow(
                            
                            column(width = 3,
                                   
                                   # Application title
                                   titlePanel("EDN Workshop"),
                                   
                                   # Let user import a custom GTF
                                   selectInput(inputId = "workshop_import_file_type_selection", 
                                               label = "Choose the annotation files to import", 
                                               choices = list("Reference GTF" = c("Ensembl"),
                                                              "Custom GTF" = c("Upload custom GTF")),
                                               width = "200px"),
                                   
                                   conditionalPanel(
                                       condition = "input.workshop_import_file_type_selection == 'Ensembl'",
                                       selectInput("workshop_genome_assembly", 
                                                   label = "Select the genome assembly to use (Ensembl release)", 
                                                   choices = list("GRCh38 (hg38)" = list.files("data") %>% gsub(pattern = "annotated_ensembl_gtf_release_(.*).txt", replacement = "\\1") %>% type.convert %>% max %>% as.character,
                                                                  "GRCh37 (hg19)" = "75",
                                                                  "NCBI36 (hg18)" = "54"
                                                   ), 
                                                   width = "200px")
                                   ),
                                   
                                   conditionalPanel(
                                       condition = "input.workshop_import_file_type_selection == 'Upload custom GTF'",
                                       fileInput("workshop_path_to_custom_gtf", "Choose GTF file",
                                                 accept = c(".gtf"), 
                                                 width = "200px")
                                   ),
                                   
                                   conditionalPanel(
                                       condition = "input.workshop_import_file_type_selection == 'Upload custom GTF'",
                                       textInput("workshop_custom_file_import_name", label = "Enter a name for this file", placeholder = "e.g. Fibroblast nanopore assembly", width = "200px")
                                   ),
                                   
                                   actionButton("workshop_button_import_annotation_file", "Import", icon = icon("file-import"),
                                                class = "btn btn-primary", width = "200px"),
                                   
                                   br(),
                                   br(),
                                   
                                   # magnetisation options required
                                   div(style = "align: left; text-align: center; width: 200px;", 
                                       
                                       shinyWidgets::dropdownButton(
                                           
                                           textInput("workshop_left_query_end_shift", label = "Left (5') query end shift (default = 0)", value = 0, placeholder = "e.g. 0"),
                                           
                                           textInput("workshop_right_query_end_shift", label = "Right (3') query end shift (default = 0)", value = 0, placeholder = "e.g. 0"),
                                           
                                           textInput("workshop_left_match_tolerance", label = "Left (5') match tolerance (default = 1)", value = 1, placeholder = "e.g. 0"),
                                           
                                           textInput("workshop_right_match_tolerance", label = "Right (3') match tolerance (default = 1)", value = 1, placeholder = "e.g. 0"),
                                           
                                           circle = FALSE, status = "info", icon = icon("wrench"),
                                           label = "Tweak match tolerances",
                                           inline = FALSE
                                       )
                                       
                                   ),
                                   
                                   br(),
                                   
                                   selectInput("workshop_range_type", 
                                               label = "Select the type of range to describe", 
                                               choices = list("Exon", 
                                                              "Junction"), 
                                               width = "200px", 
                                   ),
                                   
                                   textInput("workshop_input_range", label = "Enter genome-relative co-ordinates", placeholder = "e.g. 16:2756607-2757471:+", width = "200px"),
                                   
                                   actionButton("workshop_add_user_range", "Add range", icon = icon("plus"),
                                                class = "btn btn-primary", width = "200px"),
                                   
                                   #        div(style = "text-align:center;
                                   # box-shadow: 0px 0px 0px #888888;
                                   # padding-top:40px;
                                   # position:relative;",
                                   #            verbatimTextOutput("console_output", placeholder = TRUE)
                                   #        ),
                                   
                                   div(style = "text-align: center;
        box-shadow: 0px 0px 0px #888888;
        width: 200px;
        padding-top: 40px;
        position: relative;",
                                       verbatimTextOutput("workshop_nomenclature_output", placeholder = TRUE)
                                   ),
                                   
                                   div(style = "text-align: center;
        box-shadow: 0px 0px 0px #888888;
        width: 200px;
        height: 400px;
        padding-top: 40px;
        position: relative;",
                                       verbatimTextOutput("workshop_graph_click_info", placeholder = TRUE)
                                   )
                                   
                            ),
                            
                            column(width = 9,
                                   
                                   shinyWidgets::materialSwitch(inputId = "workshop_plot_is_active", label = "Interactive", status = "success", value = FALSE, right = FALSE, inline =TRUE),
                                   
                                   div(style = "display: inline-block; vertical-align: middle; align: left;",
                                       br()
                                   ),
                                   
                                   # manage user ranges
                                   div(style = "display: inline-block; vertical-align: middle; align: left;",
                                       uiOutput("workshop_reactive_UI_1")
                                   ),
                                   
                                   div(style = "display: inline-block; vertical-align: middle; align: left;",
                                       br()
                                   ),
                                   
                                   shinyWidgets::dropdownButton(
                                       
                                       sliderInput("workshop_slider_plot_width", "Plot width:",
                                                   min = 100, max = 4000, step = 100,
                                                   value = 800),
                                       sliderInput("workshop_slider_plot_height", "Plot height:",
                                                   min = 100, max = 4000, step = 100,
                                                   value = 1500),
                                       
                                       sliderInput("workshop_slider_plot_x_scale", "Base zoom x-axis:",
                                                   min = -5, max = 1.5, step = 0.01,
                                                   value = 1),
                                       sliderInput("workshop_slider_plot_y_scale", "Base zoom y-axis:",
                                                   min = -5, max = 5, step = 0.01,
                                                   value = 1.45),
                                       
                                       sliderInput("workshop_slider_plot_x_offset", "Base offset left/right:",
                                                   min = -5, max = 5, step = 0.01,
                                                   value = 0),
                                       sliderInput("workshop_slider_plot_y_offset", "Base offset up/down:",
                                                   min = -5, max = 5, step = 0.01,
                                                   value = 0),
                                       
                                       sliderInput("workshop_slider_table_height", "Height of table:",
                                                   min = 100, max = 3000, step = 50,
                                                   value = 250),
                                       sliderInput("workshop_slider_table_width", "Width of table:",
                                                   min = 100, max = 3000, step = 50,
                                                   value = 800),
                                       sliderInput("workshop_slider_table_font_size", "Table font size (%):",
                                                   min = 10, max = 800, step = 5,
                                                   value = 80),
                                       
                                       actionButton("workshop_reset_sliders", "Reset sliders",
                                                    class = "btn btn-primary", width = "200px"),
                                       
                                       circle = FALSE, status = "info", icon = icon("gear"),
                                       label = "Layout settings",
                                       inline = TRUE
                                   ),
                                   
                                   br(),
                                   
                                   div(style = "display: inline-block; vertical-align: middle; align: left;",
                                       h3("Reference exon table"),
                                   ),
                                   
                                   # generate an download box for users to save the plot
                                   div(style = "display: inline-block; position: relative; top: 15px; float: right;",
                                       downloadButton("workshop_download_table", "Save Table", style = "padding: 7.5px;")
                                   ),
                                   
                                   # exon table
                                   uiOutput("workshop_reactive_exon_table"),
                                   
                                   div(style = "display: inline-block; vertical-align: middle; align: left;",
                                       h3("Schematic")
                                   ),
                                   
                                   div(style = "display: inline-block; vertical-align: middle; align: left;",
                                       br()
                                   ),
                                   div(style = "display: inline-block; vertical-align: middle; align: left;",
                                       br()
                                   ),
                                   
                                   # manage plot panels
                                   div(style = "display: inline-block; position: relative; top: 5.5px; align: left;",
                                       uiOutput("workshop_reactive_UI_2")
                                   ),
                                   
                                   div(style = "display: inline-block; vertical-align: middle; align: left;",
                                       br()
                                   ),
                                   div(style = "display: inline-block; vertical-align: middle; align: left;",
                                       br()
                                   ),
                                   
                                   div(style = "display: inline-block; position: relative; top: 6px; width: 200px;",
                                       textInput("workshop_jump_to_coords", label = NULL, placeholder = "e.g. 16:2756607-2757471")
                                   ),
                                   
                                   div(style = "display: inline-block; position: relative; top: 5px;",
                                       actionButton("workshop_coord_jump_button", "Jump to", icon = icon("fast-forward"),
                                                    class = "btn btn-primary", width = "100px")
                                   ),
                                   
                                   # generate an download box for users to save the plot
                                   div(style = "display: inline-block; position: relative; top: 15px; float: right;",
                                       downloadButton("workshop_download_plot", "Save Plot", style = "padding: 7.5px;")
                                   ),
                                   
                                   h5("Click and drag + double-click to zoom. Double-click to reset."),
                                   
                                   plotOutput("workshop_plot_output", height = 300,
                                              click = "workshop_plot_output_sglclick",
                                              dblclick = "workshop_plot_output_dblclick",
                                              brush = brushOpts(
                                                  id = "workshop_plot_output_brush",
                                                  resetOnNew = TRUE
                                              )),
                                   
                            )
                            
                        )
                        
               ), # tabPanel
               
               tabPanel("EDN Reverse Translate",
                        
                        fluidRow(column(width = 12,
                                        br(),
                                        p("Coming soon...",)
                        )
                        )
                        
               ), # tabPanel
               
               tabPanel("Library of EDN",
                        
                        fluidRow(column(width = 12,
                                        br(),
                                        p("Coming soon...",)
                        )
                        )
                        
               ) # tabPanel
               
    ) # navbarPage
    
) # fluidPage

# Define server logic required to generate the nomenclature
server <- function(input, output, session) {
    
    # AUTOMATOR ###
    
    ## import GTF
    automator_reactive_tibble_ref_gtf <- eventReactive(input$automator_import_GTF, {
        
        # tibble_ref_gtf <- data.table::fread(file = "/mnt/LTS/projects/2020_isoform_nomenclature/nomenclature_app/app_native/EDN_suite/data/annotated_ensembl_gtf_release_104.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
        
        import_tibble_ref_gtf <- data.table::fread(file = paste("data/annotated_ensembl_gtf_release_", input$automator_genome_assembly, ".txt", sep = ""), sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
        
        
        return(import_tibble_ref_gtf)
        
    } )
    
    observeEvent(input$automator_import_GTF, {
        
        showModal(modalDialog(paste("Importing release ", input$automator_genome_assembly, ". Please wait...\n", sep = ""), footer = NULL))
        
        observeEvent(automator_reactive_tibble_ref_gtf(), {
            
            output$automator_nomenclature_output <- renderPrint( { cat("Importing done.\n") })
            
            removeModal()
            
        } )
        
    } )
    
    # reactive ui
    # if the user wants to manualy enter exon co-ordinates for full-length isoform or LISs, then open up multiple text boxes for entering individual co-ords.
    output$automator_reactive_UI_1 <- renderUI( {
        
        if (input$automator_structure_type == "Manually enter exon co-ordinates") {
            
            purrr::map(.x = 1:input$automator_full_length_number_of_exons, .f = ~textInput(paste("full_length_exon_genome_relative_coordinate_", .x, sep = ""), label = paste("Enter the genome-relative co-ordinates of exon #", .x, sep = ""), placeholder = "e.g. 16:2756334-2756606", width = "300px"))
            
        } else if (input$automator_structure_type == "Local Isoform Segment (LIS)") {
            
            purrr::map(.x = 1:input$automator_LIS_number_of_exons, .f = ~textInput(paste("LIS_exon_genome_relative_coordinate_", .x, sep = ""), label = paste("Enter the genome-relative co-ordinates of alternative exon #", .x, sep = ""), placeholder = "e.g. 16:2756334-2756606", width = "300px"))
            
        } else if (input$automator_structure_type == "Variable Splice Region (VSR)") {
            
            # VSRs
            # after the user has entered the number of independent events, for EACH independent event in the region, we have to ask the user how many exons there are, since alternative regions are a collection of LISs. we're basically asking for multiple LISs.
            purrr::map(.x = 1:input$automator_alternative_region_number_of_independent_events, .f = ~textInput(paste("VSR_number_of_exons_for_LIS_", .x, sep = ""), label = paste("Enter the number of exons for LIS #", .x, sep = ""), value = 1, width = "300px"))
            
        } else if (input$automator_structure_type == "Local Splice Variation (LSV) (junctions only)") {
            
            # VSRs
            # after the user has entered the number of independent events, for EACH independent event in the region, we have to ask the user how many exons there are, since alternative regions are a collection of LISs. we're basically asking for multiple LISs.
            purrr::map(.x = 1:input$automator_alternative_region_number_of_independent_events, .f = ~textInput(paste("LSV_junction_genome_relative_coordinate_", .x, sep = ""), label = paste("Enter the genome-relative co-ordinates of constitutent junction #", .x, sep = ""), placeholder = "e.g. 16:2756607-2757471", width = "300px"))
            
        } # else if
        
    } ) # renderUI
    
    output$automator_reactive_UI_2 <- renderUI( {
        
        
        if (input$automator_structure_type == "Variable Splice Region (VSR)") {
            
            reactive_automator_VSR_number_of_exons_for_each_LIS <- reactive({
                
                purrr::map(
                    .x = 1:input$automator_alternative_region_number_of_independent_events,
                    .f = ~input[[paste("VSR_number_of_exons_for_LIS_", .x, sep = "")]])
                
            })
            
            # finishing off the VSR options...
            # after use has entered the number of exons for each LIS,
            # then map across each LIS, creating textbox options for co-ordinates of exons
            
            purrr::imap(
                .x = reactive_automator_VSR_number_of_exons_for_each_LIS(),
                .f = function(a1, a2) {
                    
                    number_of_exons <- a1 %>% paste %>% type.convert
                    
                    purrr::imap(
                        .x = 1:number_of_exons,
                        .f = function(b1, b2) {
                            
                            textInput(paste("VSR_exon_genome_relative_coordinate_exon_number_", b2, "_LIS_number_", a2, sep = ""),
                                      label = paste("Enter the genome-relative co-ordinates of alternative exon #", b2, " in LIS #", a2, sep = ""),
                                      placeholder = "16:2756334-2756606", 
                                      width = "300px")
                            
                            # print(b1)
                            # print(b2)
                            
                            # print(a1)
                            # print(a2)
                            
                        } )
                    
                } )
            
        }
        
    } )
    
    observeEvent(input$automator_button_execute, {
        
        showModal(modalDialog(paste("Calculating nomenclature... \n", sep = ""), footer = NULL))
        
        tibble_ref_gtf <- automator_reactive_tibble_ref_gtf()
        
        # PROCESS INPUTS ####
        automator_input_left_query_end_shift <- input$automator_left_query_end_shift %>% type.convert
        automator_input_right_query_end_shift <- input$automator_right_query_end_shift %>% type.convert
        automator_input_left_match_tolerance <- input$automator_left_match_tolerance %>% type.convert
        automator_input_right_match_tolerance <- input$automator_right_match_tolerance %>% type.convert
        
        automator_input_full_length_recon_gtf <- input$automator_full_length_recon_gtf
        automator_input_structure_type <- input$automator_structure_type
        automator_input_alternative_exon_coords <- input$automator_alternative_exon_coords
        automator_input_alternative_junc_coords <- input$automator_alternative_junc_coords
        
        automator_input_alternative_event_region <- input$automator_alternative_event_region
        
        automator_input_structure_type <- input$automator_structure_type
        
        # invoke FLI/GTF pipeline
        if (automator_input_structure_type == "Upload assembled transcriptome (GTF)") {
            
            # DEBUG ###
            # path_recon_GTF <- "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/BM_MSC_to_OB_12d_r1_denovo_reconstructed_transcriptome.gtf"
            
            # tibble_FLI_chr_start_end_strand <- rtracklayer::import(con = "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/BM_MSC_to_OB_12d_r1_denovo_reconstructed_transcriptome.gtf", format = "gtf") %>% as_tibble %>% .[.$transcript_id == "BM_MSC_to_OB_12d_r1_Aligned.2.1" & .$type == "exon", ] %>% dplyr::rename("chr" = "seqnames") %>% dplyr::mutate_if(is.factor, as.character)
            #
            # automator_input_left_query_end_shift <- 0
            # automator_input_right_query_end_shift <- 0
            # automator_input_left_match_tolerance <- 1
            # automator_input_right_match_tolerance <- 1
            # automator_ncores <- "2x32"
            ###########
            
            automator_path_recon_GTF <- input$automator_path_full_length_recon_gtf$datapath
            
            output$automator_nomenclature_output <- renderPrint( {
                # manage parrallellisation
                automator_ncores <- input$automator_ncores
                
                if (automator_ncores != 1) {
                    
                    if (grepl(x = automator_ncores, pattern = "x") == FALSE) {
                        
                        if (automator_ncores != 0) {
                            automator_number_of_workers <- automator_ncores
                            cat(future::availableCores(), "cores will be used\n")
                        } else {
                            automator_number_of_workers <- future::availableCores()
                            cat(future::availableCores(), "cores will be used\n")
                        }
                        
                    } else if (grepl(x = automator_ncores, pattern = "x") == TRUE) {
                        
                        plan(list(tweak(multiprocess, workers = min(automator_ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert, Inf), gc = TRUE),
                                  tweak(multiprocess, workers = min(automator_ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert, Inf), gc = TRUE))
                        )
                        
                        cat((automator_ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert) * (automator_ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert), "cores will be used in total\n")
                        cat("first layer:", automator_ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert, "cores\n")
                        cat("second layer:", automator_ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert, "cores\n")
                        
                    }
                    
                }
                
            } )
            
            tibble_FLI_chr_start_end_strand <- rtracklayer::import(con = automator_path_recon_GTF, format = "gtf") %>% .[.$type == "exon", ] %>% as_tibble %>% dplyr::rename("chr" = "seqnames") %>% dplyr::mutate_if(is.factor, as.character) %>% head(n = 50000)
            # remove chrX... etc if required
            tibble_FLI_chr_start_end_strand$chr <- gsub(x = tibble_FLI_chr_start_end_strand$chr, pattern = "chr(.*)", replacement = "\\1")
            # change chromosome M to MT.
            tibble_FLI_chr_start_end_strand[tibble_FLI_chr_start_end_strand$chr == "M", "chr"] <- "MT"
            
            # list-ify reference GTF by choromosome
            list_tibble_ref_gtf_by_chr <- tibble_ref_gtf %>%
                dplyr::group_split(seqnames) %>%
                set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$seqnames %>% unique) %>% unlist)
            # list-ify FLI table by chromosome
            list_FLI_table_by_chr <- tibble_FLI_chr_start_end_strand %>%
                dplyr::group_split(chr) %>%
                set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$chr %>% unique) %>% unlist)
            # list-ify again by transcript_id
            list_FLI_table_by_chr_by_transcript_id <- purrr::map(
                .x = list_FLI_table_by_chr,
                .f = function(a1) {
                    
                    return(a1 %>%
                               dplyr::group_split(transcript_id) %>%
                               set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$transcript_id %>% unique) %>% unlist))
                    
                } )
            
            output$automator_nomenclature_output <- renderPrint( { cat("naming each transcript in the GTF...\n") })
            
            reactive_tibble_nomenclature <- reactive({
                
                suppressMessages(suppressWarnings(    # execute naming in parallel
                    list_nomenclature <- furrr::future_map2(
                        .x = list_FLI_table_by_chr_by_transcript_id,
                        .y = list_tibble_ref_gtf_by_chr[names(list_FLI_table_by_chr_by_transcript_id)],
                        .f = function(a1, a2){
                            
                            list_nomenclature_by_chr <- furrr::future_map2(
                                .x = a1,
                                .y = names(a1),
                                .f = function(b1, b2) {
                                    
                                    # cat(b2, "\n")
                                    
                                    FLI_organise_matching(b1, tibble_gtf_table = a2, left_query_shift = automator_input_left_query_end_shift, right_query_shift = automator_input_right_query_end_shift, left_tolerance = automator_input_left_match_tolerance, right_tolerance = automator_input_right_match_tolerance) %>%
                                        return
                                    
                                } )
                            
                        } )
                    
                ))
                
                tibble_nomenclature <- list_nomenclature %>% flatten %>% flatten %>% as_tibble(.name_repair = "unique") %>% t %>% as_tibble(rownames = "transcript_id") %>% setNames(c("transcript_id", "shorthand"))
                
                return(tibble_nomenclature)
                
            }) 
            
            # percolate and tibblise
            output$automator_named_GTF_entries <- renderTable(reactive_tibble_nomenclature())
            
            output$automator_nomenclature_output <- renderPrint( { cat("Naming done. Click below to download results.\n") } )
            
            # make tibble available for download
            output$automator_download_GTF_naming <- downloadHandler(
                filename = function() {
                    paste(Sys.time() %>% as.numeric, ".tab", sep = "")
                },
                content = function(file) {
                    write.table(reactive_tibble_nomenclature(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
                } )
            
            # invoke FLI/manual pipeline
        } else if (automator_input_structure_type == "Manually enter exon co-ordinates") {
            
            reactive_FLI_exon_genome_relative_coordinates <- reactive({
                
                purrr::map(.x = 1:input$automator_full_length_number_of_exons, .f = ~input[[paste("full_length_exon_genome_relative_coordinate_", .x, sep = "")]])
                
            })
            
            vector_FLI_exon_genome_relative_coordinates <- reactive_FLI_exon_genome_relative_coordinates() %>% unlist
            
            # cat(vector_FLI_exon_genome_relative_coordinates)
            
            # DEBUG ###
            # SRRM2-v1 as is
            # vector_FLI_exon_genome_relative_coordinates <- c("16:2752638-2752846:+", "16:2756334-2756606:+", "16:2757472-2757579:+", "16:2757781-2757945:+", "16:2758470-2758547:+", "16:2758985-2759047:+", "16:2759140-2759172:+", "16:2759352-2759402:+", "16:2759569-2759661:+", "16:2760301-2760499:+", "16:2761561-2768261:+", "16:2768997-2769284:+", "16:2770352-2770465:+", "16:2770604-2770717:+", "16:2770858-2771412:+")
            # with stuff changed
            # vector_FLI_exon_genome_relative_coordinates <- c("16:2752638-2752846:+", "16:2756304-2756606:+", "16:2757472-2757584:+", "16:2757801-2757850:+", "16:2758470-2758547:+", "16:2758985-2759047:+", "16:2759140-2759172:+", "16:2759352-2759402:+", "16:2760301-2760499:+", "16:2761561-2768261:+", "16:2768997-2770465:+", "16:2770604-2770717:+", "16:2770858-2771412:+")
            # vector_FLI_exon_genome_relative_coordinates <- "1:39340502-39340717"
            #
            # automator_input_left_query_end_shift <- 0
            # automator_input_right_query_end_shift <- 0
            # automator_input_left_match_tolerance <- 0
            # automator_input_right_match_tolerance <- 0
            ###########
            
            output$automator_nomenclature_output <- renderText( {
                
                triage_input_coordinates(vector_input_coordinates = vector_FLI_exon_genome_relative_coordinates %>% unlist, vector_of_expected_chromosomes = tibble_ref_gtf$seqnames %>% unique, expect_stranded = TRUE)
                
                # create tibble of chr start end strand
                vector_isoform_chr <- gsub(x = vector_FLI_exon_genome_relative_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
                vector_isoform_start <- gsub(x = vector_FLI_exon_genome_relative_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2")
                vector_isoform_end <- gsub(x = vector_FLI_exon_genome_relative_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3")
                vector_isoform_strand <- gsub(x = vector_FLI_exon_genome_relative_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
                
                tibble_FLI_chr_start_end_strand <- tibble("chr" = vector_isoform_chr,
                                                          "start" = vector_isoform_start,
                                                          "end" = vector_isoform_end,
                                                          "strand" = vector_isoform_strand) %>% type_convert
                
                # observeEvent(input$automator_button_execute, {
                
                
                paste("Suggested shorthand notation: \n", FLI_organise_matching(tibble_FLI_chr_start_end_strand, tibble_gtf_table = tibble_ref_gtf, left_query_shift = automator_input_left_query_end_shift, right_query_shift = automator_input_right_query_end_shift, left_tolerance = automator_input_left_match_tolerance, right_tolerance = automator_input_right_match_tolerance), "\n", sep = "") 
                
            } )
            
            # })
            
            # invoke VSR pipeline
        } else if (automator_input_structure_type == "Variable Splice Region (VSR)") {
            
            reactive_VSR_number_of_exons_for_each_LIS <- reactive({
                
                purrr::map(
                    .x = 1:input$automator_alternative_region_number_of_independent_events,
                    .f = ~input[[paste("VSR_number_of_exons_for_LIS_", .x, sep = "")]])
                
            } )
            
            reactive_VSR_exon_genome_relative_coordinates <- reactive({
                
                purrr::imap(
                    .x = reactive_VSR_number_of_exons_for_each_LIS(),
                    .f = function(a1, a2) {
                        number_of_exons <- a1 %>% paste %>% type.convert
                        
                        purrr::imap(
                            .x = 1:number_of_exons,
                            .f = function(b1, b2) {
                                
                                input[[paste("VSR_exon_genome_relative_coordinate_exon_number_", b2, "_LIS_number_", a2, sep = "")]]
                                
                            } )
                        
                    } )
                
            } )
            
            list_of_VSR_exon_genome_relative_coordinates <- reactive_VSR_exon_genome_relative_coordinates()
            
            # DEBUG ###
            # automator_input_alternative_event_region <- "4:82425668-82426036:*"
            # list_of_VSR_exon_genome_relative_coordinates <- list("4:82425668-82426036:*", "4:82424884-82426036:*", "4:82424884-82425562:*")
            
            # automator_input_alternative_event_region <- "7:153405599-153413612:*"
            # list_of_VSR_exon_genome_relative_coordinates <- list(c("7:153406960-153407090:*", "7:153399920-153405492:*"), "7:153400547-153401340:*", c("7:153405599-153406000:*", "7:153410612-153413612:*", "7:153406964-153407086:*"))
            
            # automator_input_alternative_event_region <- "X:11759402-11759875:*"
            # list_of_VSR_exon_genome_relative_coordinates <- list("X:11759713-11759792:*", "X:11759402-11759875:*")
            
            # automator_input_left_query_end_shift <- 0
            # automator_input_right_query_end_shift <- 0
            # automator_input_left_match_tolerance <- 1
            # automator_input_right_match_tolerance <- 1
            ###########
            
            triage_input_coordinates(vector_input_coordinates = automator_input_alternative_event_region, vector_of_expected_chromosomes = tibble_ref_gtf$seqnames %>% unique, expect_stranded = TRUE) 
            triage_input_coordinates(vector_input_coordinates = list_of_VSR_exon_genome_relative_coordinates %>% unlist, vector_of_expected_chromosomes = tibble_ref_gtf$seqnames %>% unique, expect_stranded = TRUE)
            
            output$automator_nomenclature_output <- renderText( {
                
                # create list of tibble of exon start/ends
                list_of_exon_start_end_tibbles <- purrr::map(
                    .x = list_of_VSR_exon_genome_relative_coordinates,
                    .f = function(a1) {
                        
                        vector_exon_coords <- a1 %>% unlist
                        
                        # vector_query_chr <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
                        # vector_query_VSR_start <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2")
                        # vector_query_VSR_end <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3")
                        # vector_query_strand <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
                        return(tibble("chr" = gsub(x = vector_exon_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1"),
                                      "start" = gsub(x = vector_exon_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2"),
                                      "end" = gsub(x = vector_exon_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3"),
                                      "strand" = gsub(x = vector_exon_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")) %>% type_convert)
                        
                    } )
                
                paste("Suggested shorthand notation: \n", VSR_LIS_organise_exon_naming(VSR_coordinates = automator_input_alternative_event_region, list_tibble_exon_start_end_per_LIS = list_of_exon_start_end_tibbles, tibble_gtf_table = tibble_ref_gtf, left_query_shift = automator_input_left_query_end_shift, right_query_shift = automator_input_right_query_end_shift, left_tolerance = automator_input_left_match_tolerance, right_tolerance = automator_input_right_match_tolerance), "\n", sep = "") 
                
            })
            
            # invoke VSR pipeline
        } else if (input$automator_structure_type == "Local Isoform Segment (LIS)") {
            
            reactive_LIS_exon_genome_relative_coordinates <- reactive({
                
                purrr::map(.x = 1:input$automator_LIS_number_of_exons, .f = ~input[[paste("LIS_exon_genome_relative_coordinate_", .x, sep = "")]])
                
            })
            
            vector_LIS_exon_genome_relative_coordinates <- reactive_LIS_exon_genome_relative_coordinates() %>% unlist
            
            # DEBUG ###
            # automator_input_alternative_event_region <- "9:137613615-137614211:*"
            # list_of_LIS_exon_genome_relative_coordinates <- list(c("9:137613770-137614031:*", "9:137613770-137614139:*", "9:137613770-137613980:*"))
            #
            # automator_input_left_query_end_shift <- 0
            # automator_input_right_query_end_shift <- 0
            # automator_input_left_match_tolerance <- 0
            # automator_input_right_match_tolerance <- 0
            ###########
            
            output$automator_nomenclature_output <- renderText( { 
                
                triage_input_coordinates(vector_input_coordinates = automator_input_alternative_event_region, vector_of_expected_chromosomes = tibble_ref_gtf$seqnames %>% unique, expect_stranded = TRUE)
                triage_input_coordinates(vector_input_coordinates = vector_LIS_exon_genome_relative_coordinates, vector_of_expected_chromosomes = tibble_ref_gtf$seqnames %>% unique, expect_stranded = TRUE) 
                
                # observeEvent(input$automator_button_execute, {
                
                # create list of tibble of exon start/ends
                list_of_exon_start_end_tibbles <- purrr::map(
                    .x = list(vector_LIS_exon_genome_relative_coordinates),
                    .f = function(a1) {
                        
                        # vector_query_chr <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
                        # vector_query_VSR_start <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2")
                        # vector_query_VSR_end <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3")
                        # vector_query_strand <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
                        return(tibble("start" = gsub(x = a1, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2"),
                                      "end" = gsub(x = a1, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3")) %>% type_convert)
                        
                    } )
                
                paste("Suggested shorthand notation: \n", VSR_LIS_organise_exon_naming(VSR_coordinates = automator_input_alternative_event_region, list_tibble_exon_start_end_per_LIS = list_of_exon_start_end_tibbles, tibble_gtf_table = tibble_ref_gtf, left_query_shift = automator_input_left_query_end_shift, right_query_shift = automator_input_right_query_end_shift, left_tolerance = automator_input_left_match_tolerance, right_tolerance = automator_input_right_match_tolerance), "\n", sep = "")
                
            })
            
            # invoke single exon pipeline
        } else if (automator_input_structure_type == "Alternative exon") {
            
            # DEBUG ###
            # automator_input_alternative_exon_coords <- "7:7157427-7234302:*"
            # 
            # automator_input_left_query_end_shift <- 0
            # automator_input_right_query_end_shift <- 0
            # automator_input_left_match_tolerance <- 1
            # automator_input_right_match_tolerance <- 1
            ###########
            
            output$automator_nomenclature_output <- renderText( {
                
                triage_input_coordinates(vector_input_coordinates = automator_input_alternative_exon_coords, vector_of_expected_chromosomes = tibble_ref_gtf$seqnames %>% unique, expect_stranded = TRUE)
                
                AE_query_chr <- gsub(x = automator_input_alternative_exon_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
                AE_query_start <- gsub(x = automator_input_alternative_exon_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2")
                AE_query_end <- gsub(x = automator_input_alternative_exon_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3")
                AE_query_strand <- gsub(x = automator_input_alternative_exon_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
                
                match_result <- name_a_single_exon(query_chr = AE_query_chr, query_start = AE_query_start %>% type.convert, query_end = AE_query_end %>% type.convert, query_strand = AE_query_strand, tibble_gtf_table = tibble_ref_gtf, variant_ID_override = NULL, left_query_shift = automator_input_left_query_end_shift, right_query_shift = automator_input_right_query_end_shift, left_tolerance = automator_input_left_match_tolerance, right_tolerance = automator_input_right_match_tolerance)
                
                paste("Suggested shorthand notation: \n", match_result$variant_ID_slot, " ", match_result$exon_slot, "\n", sep = "") 
                
            } )
            
        } else if (automator_input_structure_type == "Local Splice Variation (LSV) (junctions only)") {
            
            reactive_LSV_junction_genome_relative_coordinates <- reactive({
                
                purrr::map(.x = 1:input$automator_alternative_region_number_of_independent_events, .f = ~input[[paste("LSV_junction_genome_relative_coordinate_", .x, sep = "")]])
                
            } )
            
            vector_LSV_junction_genome_relative_coordinates <- reactive_LSV_junction_genome_relative_coordinates() %>% unlist
            
            # DEBUG ###
            # tibble_LSV_chr_start_end_strand <- tibble(chr = c("1", "1", "1", "1", "1", "1"), start = c(150267478, 150267789, 150267991, 150268127, 150267478, 150267789), end = c(150268696, 150268696, 150268696, 150268696, 150267714, 150267955), strand = c("-", "-", "-", "-", "-", "-"))
            # automator_left_query_shift <- 0
            # automator_right_query_shift <- 0
            # automator_left_tolerance <- 1
            # automator_right_tolerance <- 1
            ###########
            
            output$automator_nomenclature_output <- renderText( {
                
                triage_input_coordinates(vector_input_coordinates = vector_LSV_junction_genome_relative_coordinates, vector_of_expected_chromosomes = tibble_ref_gtf$seqnames %>% unique, expect_stranded = TRUE)
                
                LSV_query_chr <- gsub(x = vector_LSV_junction_genome_relative_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
                LSV_query_start <- gsub(x = vector_LSV_junction_genome_relative_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2")
                LSV_query_end <- gsub(x = vector_LSV_junction_genome_relative_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3")
                LSV_query_strand <- gsub(x = vector_LSV_junction_genome_relative_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
                
                # create list of tibble of junction start/ends
                tibble_LSV_chr_start_end_strand <- tibble("chr" = LSV_query_chr, "start" = LSV_query_start, "end" = LSV_query_end, "strand" = LSV_query_strand) %>% type_convert
                
                paste("Suggested shorthand notation: \n", LSV_AJ_organise_junction_matching(tibble_LSV_coords = tibble_LSV_chr_start_end_strand, tibble_gtf_table = tibble_ref_gtf, single_junction = FALSE, left_query_shift = automator_input_left_query_end_shift, right_query_shift = automator_input_right_query_end_shift, left_tolerance = automator_input_left_match_tolerance, right_tolerance = automator_input_right_match_tolerance), "\n", sep = "")
                
            })
            
        } else if (automator_input_structure_type == "Alternative junction") {
            
            vector_LSV_junction_genome_relative_coordinates <- automator_input_alternative_junc_coords
            
            output$automator_nomenclature_output <- renderText( {
                
                triage_input_coordinates(vector_input_coordinates = vector_LSV_junction_genome_relative_coordinates, vector_of_expected_chromosomes = tibble_ref_gtf$seqnames %>% unique, expect_stranded = TRUE)
                
                # DEBUG ###
                # vector_LSV_junction_genome_relative_coordinates <- "1:180178884-180182172:+"
                # automator_input_left_query_end_shift <- 0
                # automator_input_right_query_end_shift <- 0
                # automator_input_left_match_tolerance <- 1
                # automator_input_right_match_tolerance <- 1
                ###########
                
                LSV_query_chr <- gsub(x = vector_LSV_junction_genome_relative_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
                LSV_query_start <- gsub(x = vector_LSV_junction_genome_relative_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2")
                LSV_query_end <- gsub(x = vector_LSV_junction_genome_relative_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3")
                LSV_query_strand <- gsub(x = vector_LSV_junction_genome_relative_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
                
                # create list of tibble of junction start/ends
                tibble_LSV_chr_start_end_strand <- tibble("chr" = LSV_query_chr, "start" = LSV_query_start, "end" = LSV_query_end, "strand" = LSV_query_strand) %>% type_convert
                
                paste("Suggested shorthand notation: \n", LSV_AJ_organise_junction_matching(tibble_LSV_coords = tibble_LSV_chr_start_end_strand, tibble_gtf_table = tibble_ref_gtf, single_junction = TRUE, left_query_shift = automator_input_left_query_end_shift, right_query_shift = automator_input_right_query_end_shift, left_tolerance = automator_input_left_match_tolerance, right_tolerance = automator_input_right_match_tolerance), "\n", sep = "") 
                
            })
            
        }
        
        removeModal()
        
    }, ignoreNULL = FALSE, ignoreInit = TRUE)
    
    # END AUTOMATOR ###
    
    # WORKSHOP ###
    
    # import files
    
    workshop_reactive_button_import_annotation_file <- reactive({input$workshop_button_import_annotation_file}) %>% debounce(1000) %>% throttle(1000)
    
    ## deal with the annotation tracks
    ## create metadata
    workshop_reactiveValues_annotation_files <- reactiveValues( 
        "annotation_files" = list(
            "reference_gtf" = list(),
            "custom_gtf" = list()
        )
    )
    
    workshop_reactiveValues_custom_file_import_name <- reactiveValues(
        "workshop_custom_file_import_name" = character()
    )
    
    observe( {
        workshop_reactiveValues_custom_file_import_name$workshop_custom_file_import_name <- input$workshop_custom_file_import_name
    } )
    
    ## import reference/custom GTF
    observeEvent(workshop_reactive_button_import_annotation_file(), {
        
        if (input$workshop_import_file_type_selection == "Ensembl") {
            
            showModal(modalDialog(paste("Importing release ", input$workshop_genome_assembly, ". Please wait...\n", sep = ""), footer = NULL))
            
        } else if (input$workshop_import_file_type_selection == "Upload custom GTF" & is.null(input$workshop_path_to_custom_gtf$datapath) == FALSE) {
            
            showModal(modalDialog(paste("Importing custom GTF. Please wait...\n", sep = ""), footer = NULL))
            
        }
        
    }, ignoreNULL = FALSE, ignoreInit = TRUE )
    
    workshop_reactive_temp_imported_annotation_tibble <- eventReactive(workshop_reactive_button_import_annotation_file(), {
        
        # tibble_ref_gtf <- data.table::fread(file = "/mnt/LTS/projects/2020_isoform_nomenclature/nomenclature_app/app_native/EDN_suite/data/annotated_ensembl_gtf_release_104.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
        
        if (input$workshop_import_file_type_selection == "Ensembl") {
            
            import_tibble_ref_gtf <- data.table::fread(file = paste("data/annotated_ensembl_gtf_release_", input$workshop_genome_assembly, ".txt", sep = ""), sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
            
            workshop_reactiveValues_custom_file_import_name$workshop_custom_file_import_name <- paste("ensembl_", input$workshop_genome_assembly, sep = "")
            
            return(import_tibble_ref_gtf)
            
        } else if (input$workshop_import_file_type_selection == "Upload custom GTF" & is.null(input$workshop_path_to_custom_gtf$datapath) == FALSE) {
            
            # import_tibble_custom_gtf <- rtracklayer::import(con = "/mnt/LTS/projects/2019_hmsc_spliceome/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/BM_MSC_to_OB_3d_denovo_reconstructed_stringtiemerged.gtf", format = "gtf") %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character) %>% readr::type_convert() 
            
            workshop_path_recon_GTF <- input$workshop_path_to_custom_gtf$datapath
            
            import_tibble_custom_gtf <- rtracklayer::import(con = workshop_path_recon_GTF, format = "gtf") %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character) %>% readr::type_convert() %>% dplyr::rename("exon_number" = "exon_id")
            
            # remove chrX... etc if required
            import_tibble_custom_gtf$seqnames <- gsub(x = import_tibble_custom_gtf$seqnames, pattern = "chr(.*)", replacement = "\\1")
            # change chromosome M to MT.
            import_tibble_custom_gtf[import_tibble_custom_gtf$seqnames == "M", "seqnames"] <- "MT"
            
            return(import_tibble_custom_gtf)
            
        }
        
    }, ignoreNULL = FALSE, ignoreInit = TRUE )
    
    observeEvent(workshop_reactive_temp_imported_annotation_tibble(), {
        
        # check if file has indeed been imported
        if ((input$workshop_import_file_type_selection == "Ensembl") | 
            (input$workshop_import_file_type_selection == "Upload custom GTF" & is.null(input$workshop_path_to_custom_gtf$datapath) == FALSE)) {
            
            input_annotation_tibble <- workshop_reactive_temp_imported_annotation_tibble()
            
            if (input$workshop_import_file_type_selection == "Ensembl" & !workshop_reactiveValues_custom_file_import_name$workshop_custom_file_import_name %in% names(workshop_reactiveValues_annotation_files$annotation_files$reference_gtf)) {
                workshop_reactiveValues_annotation_files$annotation_files$reference_gtf <- workshop_reactiveValues_annotation_files$annotation_files$reference_gtf %>% purrr::splice(
                    input_annotation_tibble %>% tibble::add_column("panel" = paste("Reference GTF: ", workshop_reactiveValues_custom_file_import_name$workshop_custom_file_import_name, sep = ""))
                )
                names(workshop_reactiveValues_annotation_files$annotation_files$reference_gtf)[length(workshop_reactiveValues_annotation_files$annotation_files$reference_gtf)] <- workshop_reactiveValues_custom_file_import_name$workshop_custom_file_import_name
                # automatically checkbox the added GTF
                workshop_reactiveValues_annotation_files_selected$vector_checked_items <- c(workshop_reactiveValues_annotation_files_selected$vector_checked_items, paste("reference_gtf", "|", workshop_reactiveValues_custom_file_import_name$workshop_custom_file_import_name, sep = ""))
            } else if (input$workshop_import_file_type_selection == "Upload custom GTF" & is.null(input$workshop_path_to_custom_gtf$datapath) == FALSE) {
                workshop_reactiveValues_annotation_files$annotation_files$custom_gtf <- workshop_reactiveValues_annotation_files$annotation_files$custom_gtf %>% purrr::splice(
                    input_annotation_tibble %>% tibble::add_column("panel" = paste("Custom GTF: ", workshop_reactiveValues_custom_file_import_name$workshop_custom_file_import_name, sep = ""))
                )
                names(workshop_reactiveValues_annotation_files$annotation_files$custom_gtf)[length(workshop_reactiveValues_annotation_files$annotation_files$custom_gtf)] <- workshop_reactiveValues_custom_file_import_name$workshop_custom_file_import_name
                # automatically checkbox the added GTF
                workshop_reactiveValues_annotation_files_selected$vector_checked_items <- c(workshop_reactiveValues_annotation_files_selected$vector_checked_items, paste("custom_gtf", "|", workshop_reactiveValues_custom_file_import_name$workshop_custom_file_import_name, sep = ""))
            }
            
            print(workshop_reactiveValues_annotation_files$annotation_files$reference_gtf)
            print(workshop_reactiveValues_annotation_files$annotation_files$custom_gtf)
            
            output$workshop_nomenclature_output <- renderPrint( { cat("Importing done.\n") })
            
            removeModal()
            
        }
        
    } )
    
    # reactive ui ####
    ## this is for managing the user ranges
    output$workshop_reactive_UI_1 <- renderUI( {
        
        shinyWidgets::dropdownButton(
            
            label = "Manage ranges",
            
            div(style = "display: inline-block; float: left",
                actionButton("workshop_select_user_range", "Select",
                             class = "btn btn-primary", inline = TRUE)
            ),
            
            div(style = "display: inline-block; float: right",
                actionButton("workshop_reset_user_table", "Delete all ranges", icon = icon("eraser"),
                             class = "btn btn-primary", style = "color: #fff; background-color: #ff0000", inline = TRUE)
            ),
            
            div(style = "display: inline-block; float: right",
                actionButton("workshop_delete_user_range", "Delete", icon = icon("eraser"),
                             class = "btn btn-primary", style = "color: #fff; background-color: #ff0000", inline = TRUE)
            ),
            
            br(),
            br(),
            
            if (length(workshop_reactiveValues_user_ranges$id) > 1) {
                
                div(style = "display: block; text-align: left; font-size: 150%",
                    c("Input range history")
                )
                
            },
            
            if (length(workshop_reactiveValues_user_ranges$id) > 1) {
                
                radioButtons(
                    inputId = "workshop_user_range_id_selection", 
                    label = NULL, 
                    choiceValues = workshop_reactiveValues_user_ranges$id %>% .[2:length(workshop_reactiveValues_user_ranges$id)],
                    choiceNames = purrr::map(.x = 2:length(workshop_reactiveValues_user_ranges$id), .f = function(a1) {
                        
                        range_id <- workshop_reactiveValues_user_ranges$id %>% .[a1]
                        
                        range_chr <- workshop_reactiveValues_user_ranges$chr %>% .[a1]
                        range_start <- workshop_reactiveValues_user_ranges$start %>% .[a1]
                        range_end <- workshop_reactiveValues_user_ranges$end %>% .[a1]
                        range_strand <- workshop_reactiveValues_user_ranges$strand %>% .[a1]
                        range_range_type <- workshop_reactiveValues_user_ranges$range_type %>% .[a1]
                        
                        paste("#", range_id, "    ", range_chr, ":", range_start, "-", range_end, ":", range_strand, "  -  ", range_range_type, sep = "") %>% return
                        
                    } ) %>% unlist,
                    selected = 1)
                
            },
            circle = FALSE, status = "info", icon = icon("file-alt"),
            width = 500
        )
        
    } ) # renderUI
    
    # reactive ui ####
    ## this is for managing the plot panels
    output$workshop_reactive_UI_2 <- renderUI( {
        
        shinyWidgets::dropdownButton(
            
            label = "Show/hide panels",
            
            checkboxGroupInput(
                inputId = "workshop_select_plot_panels",
                label = NULL,
                inline = FALSE,
                choices = NULL,
                choiceNames = purrr::map2(
                    .x = workshop_reactiveValues_annotation_files$annotation_files,
                    .y = names(workshop_reactiveValues_annotation_files$annotation_files),
                    .f = function(a1, a2) {
                        
                        purrr::map2(
                            .x = a1,
                            .y = names(a1),
                            .f = function(b1, b2) {
                                
                                paste(a2 %>% stringr::str_to_sentence() %>% gsub(pattern = "gtf", replacement = "GTF") %>% gsub(pattern = "_", replacement = " "), ": ", b2, sep = "")  %>%
                                    return
                                
                            } )
                        
                    } ) %>% unlist(use.names = FALSE),
                choiceValues = purrr::map2(
                    .x = workshop_reactiveValues_annotation_files$annotation_files,
                    .y = names(workshop_reactiveValues_annotation_files$annotation_files),
                    .f = function(a1, a2) {
                        
                        purrr::map2(
                            .x = a1,
                            .y = names(a1),
                            .f = function(b1, b2) {
                                
                                paste(a2, "|", b2, sep = "") %>%
                                    return
                                
                            } )
                        
                    } ) %>% unlist(use.names = FALSE),
                selected = workshop_reactiveValues_annotation_files_selected$vector_checked_items %>% .[. %in% (purrr::map2(
                    .x = workshop_reactiveValues_annotation_files$annotation_files,
                    .y = names(workshop_reactiveValues_annotation_files$annotation_files),
                    .f = function(a1, a2) {
                        
                        purrr::map2(
                            .x = a1,
                            .y = names(a1),
                            .f = function(b1, b2) {
                                
                                paste(a2, "|", b2, sep = "") %>%
                                    return
                                
                            } )
                        
                    } ) %>% unlist(use.names = FALSE))]
            ),
            
            circle = FALSE, status = "info", icon = icon("fa-eye"),
            width = 500
        )
        
    } ) # renderUI
    
    workshop_reactive_slider_table_height <- reactive({input$workshop_slider_table_height})
    workshop_reactive_slider_table_width <- reactive({input$workshop_slider_table_width})
    workshop_reactive_slider_table_font_size <- reactive({input$workshop_slider_table_font_size})
    
    output$workshop_reactive_exon_table <- renderUI( {
        
        div(style = paste("height: ", as.numeric(workshop_reactive_slider_table_height()), "px; width: ", as.numeric(workshop_reactive_slider_table_width()), "px; overflow-y: scroll; overflow-x: scroll; font-size: ", as.numeric(workshop_reactive_slider_table_font_size()), "%", sep = ""), 
            dataTableOutput("workshop_ref_table_output")
        )
        
    } ) # renderUI
    
    # handle plot
    workshop_reactive_plot_width <- reactive({input$workshop_slider_plot_width}) %>% debounce(1000) %>% throttle(1000)
    workshop_reactive_plot_height <- reactive({input$workshop_slider_plot_height}) %>% debounce(1000) %>% throttle(1000)
    
    workshop_reactive_plot_x_scale <- reactive({input$workshop_slider_plot_x_scale}) %>% debounce(1000) %>% throttle(1000)
    workshop_reactive_plot_y_scale <- reactive({input$workshop_slider_plot_y_scale}) %>% debounce(1000) %>% throttle(1000)
    
    workshop_reactive_plot_x_offset <- reactive({input$workshop_slider_plot_x_offset}) %>% debounce(1000) %>% throttle(1000)
    workshop_reactive_plot_y_offset <- reactive({input$workshop_slider_plot_y_offset}) %>% debounce(1000) %>% throttle(1000)
    
    observeEvent(input$workshop_reset_sliders, {
        updateSliderInput(session, "workshop_slider_plot_width", value = 800)
        updateSliderInput(session, "workshop_slider_plot_height", value = 1500)
        
        updateSliderInput(session, "workshop_slider_plot_x_scale", value = 1)
        updateSliderInput(session, "workshop_slider_plot_y_scale", value = 0)
        
        updateSliderInput(session, "workshop_slider_plot_x_offset", value = 0)
        updateSliderInput(session, "workshop_slider_plot_y_offset", value = 0)
        
        updateSliderInput(session, "workshop_slider_table_height", value = 250)
        updateSliderInput(session, "workshop_slider_table_width", value = 800)
        updateSliderInput(session, "workshop_slider_table_font_size", value = 80)
    } )
    
    # handle user selection
    workshop_reactive_user_range_id_selection <- reactive({input$workshop_user_range_id_selection})
    
    # observe( {
    
    # initialise user table
    workshop_reactiveValues_user_ranges <- reactiveValues( 
        id = 0,
        chr = "3",
        start = 11807833,
        end = 11809682,
        strand = "*",
        range_type = "exon",
        panel = "user_ranges"
    )
    
    workshop_reactiveValues_selected_user_range <- reactiveValues( 
        id = 0,
        chr = "3",
        start = 11807833,
        end = 11809682,
        strand = "*",
        range_type = "exon",
        panel = "user_ranges"
    )
    
    # deal with user editing plot panels
    # we also have to keep a record of what has been selected in the past, so the the checkbox doesn't automatically reset itself when files are loaded of removed.
    workshop_reactiveValues_annotation_files_selected <- reactiveValues(
        "vector_checked_items" = c(),
        "annotation_files" = list()
    )
    
    observeEvent(input$workshop_select_plot_panels, {
        
        if (input$workshop_select_plot_panels %>% length > 0 | workshop_reactiveValues_annotation_files_selected$annotation_files %>% length > 0) {
            
            workshop_reactiveValues_annotation_files_selected$vector_checked_items <- input$workshop_select_plot_panels
            
            if (input$workshop_select_plot_panels %>% length > 0) {
                
                # list-ify
                list_selected_plot_panels <- input$workshop_select_plot_panels %>% 
                    strsplit(split = "\\|") %>%
                    purrr::map(~.x %>% unlist)
                
                names(list_selected_plot_panels) <- list_selected_plot_panels %>% purrr::map(~.x[1]) %>% unlist
                
                list_selected_plot_panels <- list_selected_plot_panels %>%
                    purrr::map(~.x[2])
                
                list_selected_plot_panels <- purrr::map(
                    .x = names(list_selected_plot_panels) %>% unique,
                    .f = function(a1) {
                        
                        list_selected_plot_panels[names(list_selected_plot_panels) == a1] %>% unlist %>%
                            return
                        
                    } ) %>% 
                    set_names(nm = names(list_selected_plot_panels) %>% unique)
                
                # subset the annotation files according to selection
                workshop_reactiveValues_annotation_files_selected$annotation_files <- purrr::map2(
                    .x = workshop_reactiveValues_annotation_files$annotation_files %>% .[names(list_selected_plot_panels)],
                    .y = list_selected_plot_panels,
                    .f = ~.x[.y]
                )
                
            } else {
                
                workshop_reactiveValues_annotation_files_selected$annotation_files <- list()
                
            }
            
        }
        
        print(workshop_reactiveValues_annotation_files_selected$annotation_files)
        
    }, ignoreNULL = FALSE, ignoreInit = TRUE)
    
    # deal with user adding range values
    observeEvent(input$workshop_add_user_range, {
        
        # parse input coords
        parse_result <- parse_input_coordinates(input_coordinates = input$workshop_input_range, vector_of_expected_chromosomes = c(1:22, "X", "Y", "M"))
        
        if (parse_result == "non_coord") {
            
            output$workshop_nomenclature_output <- renderText({"triage fail"})
            
        } else {
                
            input_chr <- parse_result[1]
            input_start <- parse_result[2]
            input_end <- parse_result[3]
            input_strand <- "*"
            
            workshop_reactiveValues_user_ranges$id <- c(workshop_reactiveValues_user_ranges$id, if (length(workshop_reactiveValues_user_ranges$id) == 0) {"1"} else {as.character(max(workshop_reactiveValues_user_ranges$id %>% as.numeric) + 1)} )
            workshop_reactiveValues_user_ranges$chr <- c(workshop_reactiveValues_user_ranges$chr, input_chr)
            workshop_reactiveValues_user_ranges$start <- c(workshop_reactiveValues_user_ranges$start, input_start)
            workshop_reactiveValues_user_ranges$end <- c(workshop_reactiveValues_user_ranges$end, input_end)
            workshop_reactiveValues_user_ranges$strand <- c(workshop_reactiveValues_user_ranges$strand, input_strand)
            workshop_reactiveValues_user_ranges$range_type <- c(workshop_reactiveValues_user_ranges$range_type, input$workshop_range_type)
            workshop_reactiveValues_user_ranges$panel <- c(workshop_reactiveValues_user_ranges$panel, "user_ranges")
            
            output$workshop_nomenclature_output <- renderText({"triage successful"})
            
            # automatically select the one that was added
            workshop_reactiveValues_selected_user_range$id <- workshop_reactiveValues_user_ranges$id %>% .[length(workshop_reactiveValues_user_ranges$id)]
            workshop_reactiveValues_selected_user_range$chr <- workshop_reactiveValues_user_ranges$chr %>% .[length(workshop_reactiveValues_user_ranges$id)]
            workshop_reactiveValues_selected_user_range$start <- workshop_reactiveValues_user_ranges$start %>% .[length(workshop_reactiveValues_user_ranges$id)]
            workshop_reactiveValues_selected_user_range$end <- workshop_reactiveValues_user_ranges$end %>% .[length(workshop_reactiveValues_user_ranges$id)]
            workshop_reactiveValues_selected_user_range$strand <- workshop_reactiveValues_user_ranges$strand %>% .[length(workshop_reactiveValues_user_ranges$id)]
            workshop_reactiveValues_selected_user_range$range_type <- workshop_reactiveValues_user_ranges$range_type %>% .[length(workshop_reactiveValues_user_ranges$id)]
            workshop_reactiveValues_selected_user_range$panel <- workshop_reactiveValues_user_ranges$panel %>% .[length(workshop_reactiveValues_user_ranges$id)]
            
            updateRadioButtons(session, inputId = "workshop_user_range_id_selection", selected = workshop_reactiveValues_user_ranges$id %>% .[length(workshop_reactiveValues_user_ranges$id)])
            
            print("workshop_reactiveValues_selected_user_range$id outside")
            print(workshop_reactiveValues_selected_user_range$id)
            
            # autojump to selected range
            workshop_reactiveValues_current_plot_range$chr <- workshop_reactiveValues_selected_user_range$chr
            workshop_reactiveValues_current_plot_range$start <- workshop_reactiveValues_selected_user_range$start
            workshop_reactiveValues_current_plot_range$end <- workshop_reactiveValues_selected_user_range$end
            
        }
        
    } )
    
    # deal with user selecting range from the list
    observeEvent(input$workshop_select_user_range, {
        
        workshop_user_range_id_selection <- workshop_reactive_user_range_id_selection() %>% type.convert
        
        vector_subsetting_condition <- workshop_reactiveValues_user_ranges$id == workshop_user_range_id_selection
        
        # automatically select the one that was added
        workshop_reactiveValues_selected_user_range$id <- workshop_reactiveValues_user_ranges$id %>% .[vector_subsetting_condition]
        workshop_reactiveValues_selected_user_range$chr <- workshop_reactiveValues_user_ranges$chr %>% .[vector_subsetting_condition]
        workshop_reactiveValues_selected_user_range$start <- workshop_reactiveValues_user_ranges$start %>% .[vector_subsetting_condition]
        workshop_reactiveValues_selected_user_range$end <- workshop_reactiveValues_user_ranges$end %>% .[vector_subsetting_condition]
        workshop_reactiveValues_selected_user_range$strand <- workshop_reactiveValues_user_ranges$strand %>% .[vector_subsetting_condition]
        workshop_reactiveValues_selected_user_range$range_type <- workshop_reactiveValues_user_ranges$range_type %>% .[vector_subsetting_condition]
        workshop_reactiveValues_selected_user_range$panel <- workshop_reactiveValues_user_ranges$panel %>% .[vector_subsetting_condition]
        
        print("workshop_reactiveValues_selected_user_range$id outside")
        print(workshop_reactiveValues_selected_user_range$id)
        
        # autojump to selected range
        workshop_reactiveValues_current_plot_range$chr <- workshop_reactiveValues_selected_user_range$chr
        workshop_reactiveValues_current_plot_range$start <- workshop_reactiveValues_selected_user_range$start
        workshop_reactiveValues_current_plot_range$end <- workshop_reactiveValues_selected_user_range$end
        
    } )
    
    # deal with user deleting range from the list
    observeEvent(input$workshop_delete_user_range, {
        
        workshop_user_range_id_selection <- workshop_reactive_user_range_id_selection() %>% type.convert
        
        vector_subsetting_condition <- (workshop_reactiveValues_user_ranges$id != workshop_user_range_id_selection)
        
        workshop_reactiveValues_user_ranges$id <- workshop_reactiveValues_user_ranges$id %>% .[vector_subsetting_condition]
        workshop_reactiveValues_user_ranges$chr <- workshop_reactiveValues_user_ranges$chr %>% .[vector_subsetting_condition]
        workshop_reactiveValues_user_ranges$start <- workshop_reactiveValues_user_ranges$start %>% .[vector_subsetting_condition]
        workshop_reactiveValues_user_ranges$end <- workshop_reactiveValues_user_ranges$end %>% .[vector_subsetting_condition]
        workshop_reactiveValues_user_ranges$strand <- workshop_reactiveValues_user_ranges$strand %>% .[vector_subsetting_condition]
        workshop_reactiveValues_user_ranges$range_type <- workshop_reactiveValues_user_ranges$range_type %>% .[vector_subsetting_condition]
        workshop_reactiveValues_user_ranges$panel <- workshop_reactiveValues_user_ranges$panel %>% .[vector_subsetting_condition]
        
        print("workshop_reactiveValues_user_ranges")
        print(reactiveValuesToList(workshop_reactiveValues_user_ranges))
        
        print("workshop_reactiveValues_selected_user_range")
        print(reactiveValuesToList(workshop_reactiveValues_selected_user_range))
        
        # check for out of bounds and no items left after delete
        if (!workshop_user_range_id_selection %in% workshop_reactiveValues_user_ranges$id & length(workshop_reactiveValues_user_ranges$id) > 1) {
            
            vector_subsetting_condition <- length(workshop_reactiveValues_user_ranges$id)
            
            workshop_reactiveValues_selected_user_range$id <- workshop_reactiveValues_user_ranges$id %>% .[vector_subsetting_condition]
            workshop_reactiveValues_selected_user_range$chr <- workshop_reactiveValues_user_ranges$chr %>% .[vector_subsetting_condition]
            workshop_reactiveValues_selected_user_range$start <- workshop_reactiveValues_user_ranges$start %>% .[vector_subsetting_condition]
            workshop_reactiveValues_selected_user_range$end <- workshop_reactiveValues_user_ranges$end %>% .[vector_subsetting_condition]
            workshop_reactiveValues_selected_user_range$strand <- workshop_reactiveValues_user_ranges$strand %>% .[vector_subsetting_condition]
            workshop_reactiveValues_selected_user_range$range_type <- workshop_reactiveValues_user_ranges$range_type %>% .[vector_subsetting_condition]
            workshop_reactiveValues_selected_user_range$panel <- workshop_reactiveValues_user_ranges$panel %>% .[vector_subsetting_condition]
            
            # autojump to selected range
            workshop_reactiveValues_current_plot_range$chr <- workshop_reactiveValues_selected_user_range$chr
            workshop_reactiveValues_current_plot_range$start <- workshop_reactiveValues_selected_user_range$start
            workshop_reactiveValues_current_plot_range$end <- workshop_reactiveValues_selected_user_range$end
            
        } else if (length(workshop_reactiveValues_user_ranges$id) == 1) {
            
            workshop_reactiveValues_selected_user_range$id <- workshop_reactiveValues_user_ranges$id %>% .[1]
            workshop_reactiveValues_selected_user_range$chr <- workshop_reactiveValues_user_ranges$chr %>% .[1]
            workshop_reactiveValues_selected_user_range$start <- workshop_reactiveValues_user_ranges$start %>% .[1]
            workshop_reactiveValues_selected_user_range$end <- workshop_reactiveValues_user_ranges$end %>% .[1]
            workshop_reactiveValues_selected_user_range$strand <- workshop_reactiveValues_user_ranges$strand %>% .[1]
            workshop_reactiveValues_selected_user_range$range_type <- workshop_reactiveValues_user_ranges$range_type %>% .[1]
            workshop_reactiveValues_selected_user_range$panel <- workshop_reactiveValues_user_ranges$panel %>% .[1]
            
            # autojump to selected range
            workshop_reactiveValues_current_plot_range$chr <- workshop_reactiveValues_selected_user_range$chr
            workshop_reactiveValues_current_plot_range$start <- workshop_reactiveValues_selected_user_range$start
            workshop_reactiveValues_current_plot_range$end <- workshop_reactiveValues_selected_user_range$end
            
        }
        
        if (vector_subsetting_condition > 1) {
            
            updateRadioButtons(session, inputId = "workshop_user_range_id_selection", selected = workshop_reactiveValues_user_ranges$id %>% .[vector_subsetting_condition])
            
        }
        
    } )
    
    # deal with user resetting user ranges
    observeEvent(input$workshop_reset_user_table, {
        
        workshop_reactiveValues_user_ranges$id <- workshop_reactiveValues_user_ranges$id %>% .[1]
        workshop_reactiveValues_user_ranges$chr <- workshop_reactiveValues_user_ranges$chr %>% .[1]
        workshop_reactiveValues_user_ranges$start <- workshop_reactiveValues_user_ranges$start %>% .[1]
        workshop_reactiveValues_user_ranges$end <- workshop_reactiveValues_user_ranges$end %>% .[1]
        workshop_reactiveValues_user_ranges$strand <- workshop_reactiveValues_user_ranges$strand %>% .[1]
        workshop_reactiveValues_user_ranges$range_type <- workshop_reactiveValues_user_ranges$range_type %>% .[1]
        workshop_reactiveValues_user_ranges$panel <- workshop_reactiveValues_user_ranges$panel %>% .[1]
        
        # select default range
        workshop_reactiveValues_selected_user_range$id <- workshop_reactiveValues_user_ranges$id %>% .[1]
        workshop_reactiveValues_selected_user_range$chr <- workshop_reactiveValues_user_ranges$chr %>% .[1]
        workshop_reactiveValues_selected_user_range$start <- workshop_reactiveValues_user_ranges$start %>% .[1]
        workshop_reactiveValues_selected_user_range$end <- workshop_reactiveValues_user_ranges$end %>% .[1]
        workshop_reactiveValues_selected_user_range$strand <- workshop_reactiveValues_user_ranges$strand %>% .[1]
        workshop_reactiveValues_selected_user_range$range_type <- workshop_reactiveValues_user_ranges$range_type %>% .[1]
        workshop_reactiveValues_selected_user_range$panel <- workshop_reactiveValues_user_ranges$panel %>% .[1]
        
        # autojump to default range
        # autojump to selected range
        workshop_reactiveValues_current_plot_range$chr <- workshop_reactiveValues_selected_user_range$chr
        workshop_reactiveValues_current_plot_range$start <- workshop_reactiveValues_selected_user_range$start
        workshop_reactiveValues_current_plot_range$end <- workshop_reactiveValues_selected_user_range$end
        
    } )
    
    workshop_reactiveValues_current_plot_range <- reactiveValues(
        "chr" = "3",
        "start" = 11807833,
        "end" = 11809682
    )
    
    workshop_plot_brush_ranges <- reactiveValues(
        x = c(), 
        y = c()
    )
    
    # process the viewable co-ordinate range
    observeEvent(input$workshop_coord_jump_button, {
        
        # enable interactivity automatically
        updateMaterialSwitch(session, inputId = "workshop_plot_is_active", value = TRUE)
        
        # parse input coords
        parse_result <- parse_input_coordinates(input_coordinates = input$workshop_jump_to_coords, vector_of_expected_chromosomes = c(1:22, "X", "Y", "M"))
        
        print("input$workshop_jump_to_coords")
        print(input$workshop_jump_to_coords)
        
        print("parse_result")
        print(parse_result)
        
        # this is the case when the jump-to coords aren't coords - we have to search the HGNC variant ID, HGNC gene symbol, gene_id, transcript_id and protein_id.
        if (parse_result == "non_coord" & length(workshop_reactiveValues_annotation_files_selected$annotation_files %>% purrr::flatten()) > 0) {
            
            # rbind all the annotation files together
            long_tibble_all_annotation_files <- workshop_reactiveValues_annotation_files_selected$annotation_files %>% flatten %>% rbindlist(use.names = TRUE, fill = TRUE) %>% as_tibble
            
            # list-ify by column
            long_tibble_all_annotation_files <- long_tibble_all_annotation_files %>% dplyr::select(seqnames, start, end, contains("gene_name"), contains("hgnc_stable_variant_ID"), contains("transcript_id"), contains("gene_id"), contains("protein_id"), panel )
            list_all_annotation_files_by_column <- long_tibble_all_annotation_files %>% purrr::array_tree(margin = 2)
            
            # grep each column
            list_grep_result_per_column <- list_all_annotation_files_by_column %>% purrr::map(~which(.x == input$workshop_jump_to_coords)) %>% purrr::discard(.p = ~.x %>% length == 0)
            
            if (length(list_grep_result_per_column) > 0) {
                
                # get the first match and only consider the annotation table associated with the first result
                tibble_matches <- long_tibble_all_annotation_files[list_grep_result_per_column[[1]], ] %>% .[.$panel == (long_tibble_all_annotation_files[list_grep_result_per_column[[1]][1], ] %>% .$panel), ]
                
                workshop_reactiveValues_current_plot_range$chr <- tibble_matches$seqnames %>% unique %>% .[1]
                workshop_reactiveValues_current_plot_range$start <- tibble_matches$start %>% min
                workshop_reactiveValues_current_plot_range$end <- tibble_matches$start %>% max
                
            }
            
        } else {
            
            input_chr <- parse_result[1]
            input_start <- parse_result[2]
            input_end <- parse_result[3]
            input_strand <- "*"
            
            workshop_reactiveValues_current_plot_range$chr <- input_chr
            workshop_reactiveValues_current_plot_range$start <- input_start
            workshop_reactiveValues_current_plot_range$end <- input_end
        }
        
        # output$workshop_nomenclature_output <- renderText({triage_result})
        
    } )
    
    # render ggplots having been given the coords of the range to be viewed
    workshop_reactive_final_plot <- reactive( {
        
        # TEST ###
        # vector_input_range <- "7:7157427-7234302:*"
        # selected_user_range_chr <- "7"
        # selected_user_range_start <- 7157427
        # selected_user_range_end <- 7234302
        # selected_user_range_strand <- "*"
        
        # plot_x_scale <- 0
        # plot_y_scale <- 0
        # 
        # plot_x_offset <- 0
        # plot_y_offset <- 0
        
        # update jump-to slider
        updateTextInput(
            session = session,
            inputId = "workshop_jump_to_coords",
            placeholder = paste(workshop_reactiveValues_current_plot_range$chr, ":", workshop_reactiveValues_current_plot_range$start, "-", workshop_reactiveValues_current_plot_range$end, sep = "")
        )
        
        # global_workshop_reactiveValues_current_plot_range <<- reactiveValuesToList(workshop_reactiveValues_current_plot_range)
        
        # smuggle and set up variables
        plot_x_scale <- workshop_reactive_plot_x_scale() %>% as.numeric
        plot_y_scale <- workshop_reactive_plot_y_scale() %>% as.numeric
        
        plot_x_offset <- workshop_reactive_plot_x_offset() %>% as.numeric
        plot_y_offset <- workshop_reactive_plot_y_offset() %>% as.numeric
        
        # test <- mtcars %>% as_tibble(rownames = "model") %>% dplyr::group_split(model) %>% purrr::map2(.x = ., .y = 1:length(.), .f = ~.x %>% tibble::add_column("panel" = .y))
        # 
        # list(ggplot(), ggplot2::facet_grid(panel ~ .)) %>% purrr::splice(purrr::map(.x = test, .f = function(a1) {
        # list(geom_point(data = a1, aes(x = mpg, y = cyl)),
        #     geom_label(data = a1, aes(x = mpg, y = cyl, label = cyl)))
        #     }) %>% purrr::flatten()) %>% purrr::reduce(ggplot2:::`+.gg`)
        
        # set up ALL user ranges
        tibble_all_user_ranges <- tibble(
            "id" = workshop_reactiveValues_user_ranges$id %>% as.character,
            "chr" = workshop_reactiveValues_user_ranges$chr,
            "start" = workshop_reactiveValues_user_ranges$start,
            "end" = workshop_reactiveValues_user_ranges$end,
            "strand" = workshop_reactiveValues_user_ranges$strand,
            "range_type" = workshop_reactiveValues_user_ranges$range_type,
            "panel" = "user_ranges"
        )
        
        # global_tibble_all_user_ranges <<- tibble_all_user_ranges
        
        print("tibble_all_user_ranges")
        print(tibble_all_user_ranges)
        
        # global_workshop_reactiveValues_selected_user_range <<- reactiveValuesToList(workshop_reactiveValues_selected_user_range)
        
        # set up the ranges that the user has selected to highlight and calculate distances for
        selected_user_range_chr <- workshop_reactiveValues_selected_user_range$chr
        selected_user_range_start <- workshop_reactiveValues_selected_user_range$start
        selected_user_range_end <- workshop_reactiveValues_selected_user_range$end
        selected_user_range_strand <- workshop_reactiveValues_selected_user_range$strand
        
        global_selected_user_range_chr <<- selected_user_range_chr
        global_selected_user_range_start <<- selected_user_range_start
        global_selected_user_range_end <<- selected_user_range_end
        global_selected_user_range_strand <<- selected_user_range_strand
        
        if (selected_user_range_strand == "*") {
            selected_user_range_strand <- c("+", "-")
        }
        
        # print("workshop_reactiveValues_current_plot_range$chr")
        # print(workshop_reactiveValues_current_plot_range$chr)
        # print("workshop_reactiveValues_current_plot_range$start")
        # print(workshop_reactiveValues_current_plot_range$start)
        # print("workshop_reactiveValues_current_plot_range$end")
        # print(workshop_reactiveValues_current_plot_range$end)
        
        print("plot_view_initial_x_start beginning")
        print(workshop_reactiveValues_current_plot_range$start)
        
        print("plot_view_initial_x_end beginning")
        print(workshop_reactiveValues_current_plot_range$end)
        
        # subset all tables to be plotted, for user-specified range (1.5x jump to/user range selection)
        plot_view_initial_x_start0 <- workshop_reactiveValues_current_plot_range$start - 1.5*(workshop_reactiveValues_current_plot_range$end - workshop_reactiveValues_current_plot_range$start) + (plot_x_scale^3)*(workshop_reactiveValues_current_plot_range$end - workshop_reactiveValues_current_plot_range$start)
        plot_view_initial_x_end0 <- workshop_reactiveValues_current_plot_range$end + 1.5*(workshop_reactiveValues_current_plot_range$end - workshop_reactiveValues_current_plot_range$start) - (plot_x_scale^3)*(workshop_reactiveValues_current_plot_range$end - workshop_reactiveValues_current_plot_range$start)
        
        print("plot_view_initial_x_start middle")
        print(plot_view_initial_x_start0)
        
        print("plot_view_initial_x_end middle")
        print(plot_view_initial_x_end0)
        
        print("plot_x_offset")
        print(plot_x_offset)
        
        # make sure the start coord is to the left of the end
        if (plot_view_initial_x_start0 < plot_view_initial_x_end0) {
            
            # apply offset
            plot_view_initial_x_start <- plot_view_initial_x_start0 + (plot_view_initial_x_end0 - plot_view_initial_x_end0)*plot_x_offset
            plot_view_initial_x_end <- plot_view_initial_x_end0 + (plot_view_initial_x_end0 - plot_view_initial_x_start0)*plot_x_offset
            
        } else {
            
            plot_view_initial_x_start <- workshop_reactiveValues_current_plot_range$start - 1.5*(workshop_reactiveValues_current_plot_range$end - workshop_reactiveValues_current_plot_range$start)
            plot_view_initial_x_end <- workshop_reactiveValues_current_plot_range$end + 1.5*(workshop_reactiveValues_current_plot_range$end - workshop_reactiveValues_current_plot_range$start)
            
        }
        
        # print("workshop_reactiveValues_annotation_files_selected$annotation_files")
        # print(workshop_reactiveValues_annotation_files_selected$annotation_files)
        # 
        # print("workshop_reactiveValues_current_plot_range$chr")
        # print(workshop_reactiveValues_current_plot_range$chr)
        # 
        print("plot_view_initial_x_start final")
        print(plot_view_initial_x_start)
        
        print("plot_view_initial_x_end final")
        print(plot_view_initial_x_end)
        
        # global_workshop_reactiveValues_annotation_files_selected <<- reactiveValuesToList(workshop_reactiveValues_annotation_files_selected)
        
        # plot calculation
        ## find the reference elements visible in the viewing range
        list_tibbles_track_features_visible <- purrr::map2(
            .x = workshop_reactiveValues_annotation_files_selected$annotation_files,
            .y = names(workshop_reactiveValues_annotation_files_selected$annotation_files),
            .f = function(a1, a2) {
                
                purrr::map2(
                    .x = a1,
                    .y = names(a1),
                    .f = function(b1, b2) {
                        
                        tibble_captured_in_range <- b1[which(b1$seqnames == workshop_reactiveValues_current_plot_range$chr & b1$start <= plot_view_initial_x_end & b1$end >= plot_view_initial_x_start), ]
                        
                        return(tibble_captured_in_range)
                        
                    } ) %>% 
                    purrr::keep(.p = ~.x %>% nrow > 0) %>%
                    return
                
            } )
        
        # print("list_tibbles_track_features_visible")
        # print(list_tibbles_track_features_visible)
        
        ## flatten into single list with names like: "category | name"
        list_tibbles_track_features_visible_flattened <- list_tibbles_track_features_visible %>% flatten
        
        # names(list_tibbles_track_features_visible_flattened) <- purrr::map2(
        #     .x = list_tibbles_track_features_visible, 
        #     .y = names(list_tibbles_track_features_visible), 
        #     .f = function(a1, a2) {
        #         
        #         purrr::map2(
        #             .x = a1,
        #             .y = names(a1),
        #             .f = function(b1, b2) {
        #                 
        #                 # return(paste(a2 %>% stringr::str_to_sentence() %>% gsub(pattern = "gtf", replacement = "GTF") %>% gsub(pattern = "\\_", replacement = " "), ": ", b2, sep = ""))
        #                 
        #                 
        #             } ) %>% unlist
        #         
        #     } ) %>% unlist
        
        names(list_tibbles_track_features_visible_flattened) <- purrr::map(.x = list_tibbles_track_features_visible_flattened, .f = ~.x$panel %>% unique)
        
        print("list_tibbles_track_features_visible_flattened")
        print(list_tibbles_track_features_visible_flattened)
        
        # list_tibbles_track_features_visible_flattened <- purrr::map2(
        #     .x = list_tibbles_track_features_visible_flattened,
        #     .y = names(list_tibbles_track_features_visible_flattened),
        #     .f = ~.x %>% tibble::add_column("panel" = .y)
        # )
        
        ## find the user ranges visible in the viewing range
        tibble_user_ranges_visible <- tibble_all_user_ranges[which(tibble_all_user_ranges$chr == workshop_reactiveValues_current_plot_range$chr & tibble_all_user_ranges$start <= plot_view_initial_x_end & tibble_all_user_ranges$end >= plot_view_initial_x_start), ]
        
        # global_1_list_tibbles_track_features_visible_flattened <<- list_tibbles_track_features_visible_flattened
        
        # distance shenanigans
        ## calculate distances for selected range
        ## only do this if the selected range is visible
        if (workshop_reactiveValues_selected_user_range$id %in% tibble_user_ranges_visible$id & nrow(tibble_user_ranges_visible) > 0 & length(list_tibbles_track_features_visible_flattened) > 0) {
            
            # flatten the original full GTF table
            list_tibbles_track_features_all_flattened <- workshop_reactiveValues_annotation_files_selected$annotation_files %>% flatten
            
            # DEBUG ###
            global_list_tibbles_track_features_visible_flattened_1 <<- list_tibbles_track_features_visible_flattened
            global_list_tibbles_track_features_all_flattened_1 <<- list_tibbles_track_features_all_flattened
            ###########
            
            ## having determined the plot window, calculate distances to every exon in the plot range
            list_distance_annotation_data_flattened <- purrr::pmap(
                .l = list(
                    "a1" = list_tibbles_track_features_visible_flattened,
                    "a2" = names(list_tibbles_track_features_visible_flattened),
                    "a3" = list_tibbles_track_features_all_flattened
                ),
                .f = function(a1, a2, a3) {
                    
                    # DEBUG ###
                    # a1 <- global_list_tibbles_track_features_visible_flattened_1[[1]]
                    # a2 <- names(global_list_tibbles_track_features_visible_flattened_1)[[1]]
                    # a3 <- global_list_tibbles_track_features_all_flattened_1[[1]]
                    # selected_user_range_chr <- global_selected_user_range_chr
                    # selected_user_range_start <- global_selected_user_range_start
                    # selected_user_range_end <- global_selected_user_range_end
                    ###########
                    
                    # determine visible transcript ids overlapped by user range
                    tibble_ref_transcripts_overlapped_by_user_query <- extract_overlapping_features(query_chr = selected_user_range_chr, query_start = selected_user_range_start, query_end = selected_user_range_end, query_strand = "*", tibble_gtf_table = a1, left_query_shift = input$workshop_left_query_end_shift %>% type.convert, right_query_shift = input$workshop_right_query_end_shift %>% type.convert, left_tolerance = input$workshop_left_match_tolerance %>% type.convert, right_tolerance = input$workshop_right_match_tolerance %>% type.convert, return_type = "transcript")
                    
                    print("tibble_ref_transcripts_overlapped_by_user_query")
                    print(tibble_ref_transcripts_overlapped_by_user_query)
                    
                    # check if the selected transcript overlaps a ref transcript
                    if (tibble_ref_transcripts_overlapped_by_user_query %>% nrow > 0) {
                        
                        tibble_all_exons_of_overlapped_parent_transcript <- a3[which(a3$type == "exon" & a3$transcript_id %in% (tibble_ref_transcripts_overlapped_by_user_query$transcript_id %>% unique)), ]
                        
                        print("tibble_all_exons_of_overlapped_parent_transcript")
                        print(tibble_all_exons_of_overlapped_parent_transcript)
                        
                        tibble_distance_annotations_based_on_user_query <- purrr::map2(
                            # overlapping transcript entries
                            .x = tibble_ref_transcripts_overlapped_by_user_query %>% dplyr::group_split(transcript_id),
                            # exons belonging to the transcript
                            .y = tibble_all_exons_of_overlapped_parent_transcript %>% dplyr::group_split(transcript_id) %>% set_names(nm = purrr::map(.x = ., .f = ~.x$transcript_id %>% unique) %>% unlist) %>% .[tibble_ref_transcripts_overlapped_by_user_query %>% dplyr::group_split(transcript_id) %>% purrr::map(~.x$transcript_id %>% unique) %>% unlist],
                            .f = function(a1, a2) {
                                
                                # DEBUG ###
                                # a1 <- tibble_ref_transcripts_overlapped_by_user_query %>% dplyr::group_split(transcript_id) %>% .[[2]]
                                # a2 <- tibble_all_exons_of_overlapped_parent_transcript %>% dplyr::group_split(transcript_id) %>% set_names(nm = purrr::map(.x = ., .f = ~.x$transcript_id %>% unique) %>% unlist) %>% .[tibble_ref_transcripts_overlapped_by_user_query %>% dplyr::group_split(transcript_id) %>% purrr::map(~.x$transcript_id %>% unique) %>% unlist] %>% .[[2]]
                                ###########
                                
                                vector_all_ref_vertices <- a2[, c("start", "end")] %>% unlist %>% sort
                                
                                # strategy: grow left and right ends of the user vertices until it touches a vertex. 
                                left_ref_vertex_grown_from_user_query_start <- vector_all_ref_vertices[vector_all_ref_vertices <= selected_user_range_start] %>% max
                                right_ref_vertex_grown_from_user_query_start <- vector_all_ref_vertices[vector_all_ref_vertices >= selected_user_range_start] %>% min
                                left_ref_vertex_grown_from_user_query_end <- vector_all_ref_vertices[vector_all_ref_vertices <= selected_user_range_end] %>% max
                                right_ref_vertex_grown_from_user_query_end <- vector_all_ref_vertices[vector_all_ref_vertices >= selected_user_range_end] %>% min
                                
                                tibble_vertices_with_distances <- tibble(
                                    "transcript_id" = a1$transcript_id %>% unique,
                                    "ref_vertex" = c(left_ref_vertex_grown_from_user_query_start, right_ref_vertex_grown_from_user_query_start, left_ref_vertex_grown_from_user_query_end, right_ref_vertex_grown_from_user_query_end),
                                    "query_vertex" = c(selected_user_range_start, selected_user_range_start, selected_user_range_end, selected_user_range_end)
                                ) %>% dplyr::mutate("ref_vertex_minus_query_vertex" = `ref_vertex` - `query_vertex`)
                                
                                # test for redundant overlapping distances. this happens when 1. both query ends find a ref transcript and 2. distance to left overlaps and/or distance to right overlaps.
                                if (left_ref_vertex_grown_from_user_query_end < selected_user_range_start) {
                                    tibble_vertices_with_distances <- tibble_vertices_with_distances[!(tibble_vertices_with_distances$ref_vertex == left_ref_vertex_grown_from_user_query_end & tibble_vertices_with_distances$query_vertex == selected_user_range_end), ]
                                }
                                
                                if (right_ref_vertex_grown_from_user_query_start > selected_user_range_end) {
                                    tibble_vertices_with_distances <- tibble_vertices_with_distances[!(tibble_vertices_with_distances$ref_vertex == right_ref_vertex_grown_from_user_query_start & tibble_vertices_with_distances$query_vertex == selected_user_range_start), ]
                                }
                                
                                return(tibble_vertices_with_distances)
                                
                            } ) %>% dplyr::bind_rows() %>% unique
                        
                        # make the distances directional
                        tibble_distance_annotations_based_on_user_query <- tibble_distance_annotations_based_on_user_query %>% 
                            dplyr::mutate("ref_vertex_minus_query_vertex" = purrr::map(.x = `ref_vertex_minus_query_vertex`, .f = function(x) {
                                
                                if (x < 0) {
                                    return(paste(abs(x), ">", sep = ""))
                                } else if (x > 0) {
                                    return(paste("<", abs(x), sep = ""))
                                } else {
                                    return(x)
                                }
                                
                            } ) %>% unlist  )
                        
                        tibble_distance_annotations_based_on_user_query$panel <- a2
                        
                    } else {
                        
                        tibble_distance_annotations_based_on_user_query <- tibble()
                        tibble_all_exons_of_overlapped_parent_transcript <- tibble()
                        
                    }
                    
                    # invoke tolerances
                    # tibble_selected_range_measurement_against_reference[abs(tibble_selected_range_measurement_against_reference$distance_to_ref) <= 1, "distance_to_ref"] <- 0
                    
                    ##########
                    
                    return(
                        list(
                            "tibble_distance_annotations_based_on_user_query" = tibble_distance_annotations_based_on_user_query,
                            "tibble_all_exons_of_overlapped_parent_transcript" = tibble_all_exons_of_overlapped_parent_transcript
                        )
                    )
                    
                } )
            
            list_distances_between_user_ranges_and_reference_annotations <- list_distance_annotation_data_flattened %>% purrr::map(~.x$tibble_distance_annotations_based_on_user_query) %>% purrr::keep(.p = ~.x %>% length > 0)
            
            # add back in the exons touched by the distance range
            list_tibbles_track_features_visible_flattened <- purrr::pmap(
                .l = list(
                    "a1" = list_tibbles_track_features_visible_flattened,
                    "a2" = list_distance_annotation_data_flattened %>% purrr::map(~.x$tibble_all_exons_of_overlapped_parent_transcript),
                    "a3" = names(list_tibbles_track_features_visible_flattened)
                ), .f = function(a1, a2, a3) {
                    
                    dplyr::bind_rows(a1, a2) %>% dplyr::mutate("panel" = a3) %>% 
                        return
                    
                } 
            ) %>% set_names(nm = names(list_tibbles_track_features_visible_flattened))
            
        } else {
            
            list_distances_between_user_ranges_and_reference_annotations <- list()
            
            if (length(list_tibbles_track_features_visible_flattened) == 0) {
                list_tibbles_track_features_visible_flattened <- list()
            }
            
        }
        
        number_of_transcripts_captured <- list_tibbles_track_features_visible_flattened %>% purrr::map(~.x$transcript_id %>% unique %>% length) %>% unlist %>% max
        
        # ONLY NOW we can set the y-viewing range
        plot_view_initial_y_start <- 1 + plot_y_scale*(number_of_transcripts_captured - 1)
        #  - 1.5*(number_of_transcripts_captured - 1)
        plot_view_initial_y_end <- number_of_transcripts_captured - plot_y_scale*(number_of_transcripts_captured - 1)
        #  + 1.5*(number_of_transcripts_captured - 1)
        
        # apply offset
        ## minuses because the y axis is in reversed order.
        plot_view_initial_y_start <- plot_view_initial_y_start - (plot_view_initial_y_end - plot_view_initial_y_start)*plot_y_offset
        plot_view_initial_y_end <- plot_view_initial_y_end - (plot_view_initial_y_end - plot_view_initial_y_start)*plot_y_offset
        
        return(list(
            "plot_view_initial_x_start" = plot_view_initial_x_start,
            "plot_view_initial_x_end" = plot_view_initial_x_end,
            "plot_view_initial_y_start" = plot_view_initial_y_start,
            "plot_view_initial_y_end" = plot_view_initial_y_end,
            "list_tibbles_track_features_visible_flattened" = list_tibbles_track_features_visible_flattened,
            "tibble_user_ranges_visible" = tibble_user_ranges_visible,
            "list_distances_between_user_ranges_and_reference_annotations" = list_distances_between_user_ranges_and_reference_annotations
        ) )
        
    } )
    
    workshop_plot_brush_ranges <- reactiveValues(
        x = c(), 
        y = c(),
        logical_brush_zoom_on = FALSE
    )
    
    observe( {
        
        workshop_plot_brush_ranges$x <- c(workshop_reactive_final_plot() %>% .$plot_view_initial_x_start, 
                                          workshop_reactive_final_plot() %>% .$plot_view_initial_x_end)
        
        # workshop_plot_brush_ranges$y <- c(workshop_reactive_final_plot() %>% .$plot_view_initial_y_start, 
        #                                   workshop_reactive_final_plot() %>% .$plot_view_initial_y_end)
        
    } )
    
    workshop_reactiveValues_plot_metadata <- reactiveValues(
        "list_track_data_in_viewing_range" = list(),
        "list_y_axis_scale" = list(),
        "list_y_axis_scale_initial" = list(),
        "vector_number_of_features_per_track" = list(),
        "list_distances_between_user_ranges_and_reference_annotations" = list()
    )
    
    
    
    
    # Zoomable plot + ref table output
    observe( {
        
        if (input$workshop_plot_is_active == TRUE) {
            
            # smuggle in variables
            plot_height <- workshop_reactive_plot_height()
            plot_width <- workshop_reactive_plot_width()
            
            plot_view_initial_y_start <- workshop_reactive_final_plot() %>% .$plot_view_initial_y_start
            plot_view_initial_y_end <- workshop_reactive_final_plot() %>% .$plot_view_initial_y_end
            
            # print("plot_height")
            # print(plot_height)
            # print("plot_width")
            # print(plot_width)
            
            slider_table_height <- workshop_reactive_slider_table_height()
            
            
            list_tibbles_track_features_visible_flattened <- workshop_reactive_final_plot() %>% .$list_tibbles_track_features_visible_flattened
            
            tibble_user_ranges_visible <- workshop_reactive_final_plot() %>% .$tibble_user_ranges_visible
            # print("tibble_user_ranges_visible 2")
            # print(tibble_user_ranges_visible)
            
            list_distances_between_user_ranges_and_reference_annotations <- workshop_reactive_final_plot() %>% .$list_distances_between_user_ranges_and_reference_annotations
            
            # update the metadata record of items in the viewing range
            ## track data
            workshop_reactiveValues_plot_metadata$list_track_data_in_viewing_range <- list(
                
                if (length(list_tibbles_track_features_visible_flattened) > 0){
                    list_tibbles_track_features_visible_flattened
                },
                if (nrow(tibble_user_ranges_visible) > 0){
                    list("tibble_user_ranges_visible" = tibble_user_ranges_visible)
                }
                
            ) %>% flatten
            
            if (workshop_plot_brush_ranges$logical_brush_zoom_on == FALSE) {
                
                ## get the y-axis values 
                workshop_reactiveValues_plot_metadata$list_y_axis_scale <- list(
                    
                    if (length(list_tibbles_track_features_visible_flattened[grep(x = names(list_tibbles_track_features_visible_flattened), pattern = "^Reference GTF")]) > 0) {
                        list_tibbles_track_features_visible_flattened[grep(x = names(list_tibbles_track_features_visible_flattened), pattern = "^Reference GTF")] %>% purrr::map(~.x[mixedorder(.x$hgnc_stable_variant_ID), ] %>% .$transcript_id %>% unique %>% na.omit %>% rev)
                    },
                    if (length(list_tibbles_track_features_visible_flattened[grep(x = names(list_tibbles_track_features_visible_flattened), pattern = "^Custom GTF")]) > 0) {
                        list_tibbles_track_features_visible_flattened[grep(x = names(list_tibbles_track_features_visible_flattened), pattern = "^Custom GTF")] %>% purrr::map(~.x$transcript_id %>% unique %>% na.omit %>% rev)
                    },
                    if (nrow(tibble_user_ranges_visible) > 0) {
                        list("user_ranges" = tibble_user_ranges_visible$id %>% as.character())
                    }
                    
                ) %>% flatten
                
                workshop_reactiveValues_plot_metadata$list_y_axis_scale_initial <- workshop_reactiveValues_plot_metadata$list_y_axis_scale
                
            }
            
            ## number of features per track
            workshop_reactiveValues_plot_metadata$vector_number_of_features_per_track <- workshop_reactiveValues_plot_metadata$list_y_axis_scale %>% 
                purrr::map(~.x %>% length) %>% unlist
            ## distances
            workshop_reactiveValues_plot_metadata$list_distances_between_user_ranges_and_reference_annotations <- list_distances_between_user_ranges_and_reference_annotations
            
            # DEBUG ###
            print("plot debug")
            global_workshop_reactiveValues_selected_user_range <<- reactiveValuesToList(workshop_reactiveValues_selected_user_range)
            global_list_distances_between_user_ranges_and_reference_annotations <<- list_distances_between_user_ranges_and_reference_annotations
            global_workshop_reactiveValues_current_plot_range <<- reactiveValuesToList(workshop_reactiveValues_current_plot_range)
            global_tibble_user_ranges_visible <<- tibble_user_ranges_visible
            global_workshop_reactiveValues_plot_metadata <<- reactiveValuesToList(workshop_reactiveValues_plot_metadata)
            global_2_list_tibbles_track_features_visible_flattened <<- list_tibbles_track_features_visible_flattened
            
            print("tibble_user_ranges_visible")
            print(tibble_user_ranges_visible)
            
            print("list_distances_between_user_ranges_and_reference_annotations")
            print(list_distances_between_user_ranges_and_reference_annotations)
            ###########
            
            # CREATE GGPLOT
            ggplot_final_plot <- list(
                ggplot(),
                # these mark the original viewing window
                geom_vline(colour = "green", lty = 2, xintercept = workshop_reactiveValues_current_plot_range$start),
                geom_vline(colour = "green", lty = 2, xintercept = workshop_reactiveValues_current_plot_range$end),
                theme_bw(),
                theme(text = element_text(family = "Helvetica")),
                ggplot2::xlab(paste("chr", workshop_reactiveValues_current_plot_range$chr, sep = "")),
                # brush resizing (x)
                # NOTE we CANNOT use coord_cartesion for facet-specific y. MUST use the scales.
                coord_cartesian(xlim = workshop_plot_brush_ranges$x
                                # ,
                                # ylim = c(plot_view_initial_y_start, plot_view_initial_y_end)
                )
            ) %>% 
                
                purrr::splice(
                    
                    # reference_gtf
                    if (length(list_tibbles_track_features_visible_flattened[grep(x = names(list_tibbles_track_features_visible_flattened), pattern = "^Reference GTF")]) > 0) {
                        
                        purrr::map(
                            .x = list_tibbles_track_features_visible_flattened[grep(x = names(list_tibbles_track_features_visible_flattened), pattern = "^Reference GTF")], 
                            .f = function(a1) {
                                
                                list(
                                    
                                    geom_segment(data = a1 %>% dplyr::filter(type == "transcript"), colour = "slateblue1", mapping = aes(x = start, xend = end, y = transcript_id, yend = transcript_id)),
                                    geom_text(data = a1 %>% dplyr::filter(type == "transcript"), nudge_y = 0.25, fontface = "italic", mapping = aes(x = mean(workshop_plot_brush_ranges$x), y = transcript_id, label = purrr::pmap(.l = list("b1" = strand, "b2" = hgnc_stable_variant_ID, "b3" = transcript_version), .f = function(b1, b2, b3) {if (b1 == "+") {paste("> > > > > > ", b2, " > > > > > >", sep = "")} else if (b1 == "-") {paste("< < < < < < ", b2, " < < < < < <", sep = "")} else {b2} } ) %>% unlist)),
                                    geom_segment(data = a1 %>% dplyr::filter(type == "exon"), colour = "slateblue1", mapping = aes(x = start, xend = end, y = transcript_id, yend = transcript_id), size = 10),
                                    geom_label(data = a1 %>% dplyr::filter(type == "exon"), colour = "black", nudge_y = 0.15, fontface = "bold.italic", mapping = aes(x = purrr::map2(.x = start, .y = end, .f = ~c(.x, .y) %>% mean) %>% unlist, y = transcript_id, label = paste("E", exon_number, sep = "")))
                                    
                                )
                                
                            } ) %>% purrr::flatten()
                        
                    },
                    
                    # custom_gtf
                    if (length(list_tibbles_track_features_visible_flattened[grep(x = names(list_tibbles_track_features_visible_flattened), pattern = "^Custom GTF")]) > 0) {
                        
                        purrr::map(
                            .x = list_tibbles_track_features_visible_flattened[grep(x = names(list_tibbles_track_features_visible_flattened), pattern = "^Custom GTF")],
                            .f = function(a2) {
                                
                                list(
                                    
                                    geom_segment(data = a2 %>% dplyr::filter(type == "transcript"), colour = "slateblue1", mapping = aes(x = start, xend = end, y = transcript_id, yend = transcript_id)),
                                    geom_text(data = a2 %>% dplyr::filter(type == "transcript"), nudge_y = 0.25, mapping = aes(x = mean(workshop_plot_brush_ranges$x), y = transcript_id, label = purrr::pmap(.l = list("b1" = strand, "b2" = transcript_id), .f = function(b1, b2) {if (b1 == "+") {paste("> > > > > > ", b2, " > > > > > >", sep = "")} else if (b1 == "-") {paste("< < < < < < ", b2, " < < < < < <", sep = "")} else {b2} } ) %>% unlist)),
                                    geom_segment(data = a2 %>% dplyr::filter(type == "exon"), colour = "slateblue1", mapping = aes(x = start, xend = end, y = transcript_id, yend = transcript_id), size = 10),
                                    geom_label(data = a2 %>% dplyr::filter(type == "exon"), colour = "black", nudge_y = 0.15, fontface = "bold.italic", mapping = aes(x = purrr::map2(.x = start, .y = end, .f = ~c(.x, .y) %>% mean) %>% unlist, y = transcript_id, label = paste("E", exon_number, sep = "")))
                                    
                                )
                                
                            } ) %>% purrr::flatten()
                        
                    },
                    
                    # user ranges
                    if (nrow(tibble_user_ranges_visible[tibble_user_ranges_visible$range_type == "Junction" & tibble_user_ranges_visible$id != "0", ]) > 0) {
                        geom_curve(data = tibble_user_ranges_visible[tibble_user_ranges_visible$range_type == "Junction" & tibble_user_ranges_visible$id != "0", ], colour = "grey50", size = 2, curvature = -0.25, mapping = aes(x = start, xend = end, y = id, yend = id))
                    },
                    if (nrow(tibble_user_ranges_visible[tibble_user_ranges_visible$range_type == "Exon" & tibble_user_ranges_visible$id != "0", ]) > 0) {
                        geom_segment(data = tibble_user_ranges_visible[tibble_user_ranges_visible$range_type == "Exon" & tibble_user_ranges_visible$id != "0", ], colour = "grey50", size = 5, mapping = aes(x = start, xend = end, y = id, yend = id))
                    },
                    if (nrow(tibble_user_ranges_visible) > 1) {
                        geom_text(data = tibble_user_ranges_visible[tibble_user_ranges_visible$id != "0", ], colour = "black", nudge_y = 0.25, fontface = "bold", mapping = aes(x = purrr::map2(.x = start, .y = end, .f = ~c(.x, .y) %>% mean) %>% unlist, y = id, label = id))
                    },
                    if (workshop_reactiveValues_selected_user_range$id %in% tibble_user_ranges_visible$id) {
                        geom_vline(colour = "red", lty = 2, xintercept = workshop_reactiveValues_selected_user_range$start)
                    },
                    if (workshop_reactiveValues_selected_user_range$id %in% tibble_user_ranges_visible$id) {
                        geom_vline(colour = "red", lty = 2, xintercept = workshop_reactiveValues_selected_user_range$end)
                    },
                    
                    # distance annotation
                    if (length(list_distances_between_user_ranges_and_reference_annotations) > 0) {
                        
                        purrr::map(
                            .x = list_distances_between_user_ranges_and_reference_annotations,
                            .f = function(a3) {
                                
                                list(
                                    
                                    geom_segment(data = a3[a3$ref_vertex_minus_query_vertex != 0, ], colour = "red", arrow = arrow(angle = 30), mapping = aes(x = ref_vertex, xend = query_vertex, y = transcript_id, yend = transcript_id)),
                                    geom_label(data = a3, colour = "red", nudge_y = -0.25, mapping = aes(x = purrr::map2(.x = ref_vertex, .y = query_vertex, .f = ~c(.x, .y) %>% mean) %>% unlist, y = transcript_id, label = ref_vertex_minus_query_vertex))
                                    
                                )
                                
                            } ) %>% purrr::flatten()
                        
                    },
                    
                    # FACETS 
                    if (length(list_tibbles_track_features_visible_flattened %>% flatten) > 0 | nrow(tibble_user_ranges_visible) > 0) {
                        ggplot2::facet_grid(factor(panel, level = workshop_reactiveValues_plot_metadata$list_y_axis_scale %>% names) ~ ., scales = "free_y")
                    },
                    
                    # adaptive facet aspect ratio
                    if (length(list_tibbles_track_features_visible_flattened %>% flatten) > 0 | nrow(tibble_user_ranges_visible) > 0) {
                        ggh4x::force_panelsizes(rows = workshop_reactiveValues_plot_metadata$vector_number_of_features_per_track %>%
                                                    (function(x) {
                                                        if (sum(x) != 0) {
                                                            return(x/sum(x))
                                                        } else {
                                                            return(0)}
                                                    } ) )
                    },
                    
                    # facet-specific brush resizing (y)
                    if (length(list_tibbles_track_features_visible_flattened %>% flatten) > 0 | nrow(tibble_user_ranges_visible) > 0) {
                        ggh4x::facetted_pos_scales(
                            y = workshop_reactiveValues_plot_metadata$list_y_axis_scale %>% purrr::map(~scale_y_discrete(limits = .x, breaks = .x, labels = .x))
                        )
                    }
                    
                )  %>% purrr::reduce(ggplot2:::`+.gg`)
            
            output$workshop_plot_output <- renderPlot( {
                
                ggplot_final_plot
                
            }, height = plot_height, width = plot_width )
            
            output$workshop_ref_table_output <- renderDataTable(
                {workshop_reactive_final_plot() %>% .$list_tibbles_track_features_visible_flattened %>% rbindlist(use.names = TRUE, fill = TRUE) %>%
                        dplyr::select(contains("hgnc_stable_variant_ID"), contains("transcript_version"), contains("type"), contains("exon_number"), contains("seqnames"), contains("start"), contains("end"), contains("width"), contains("strand"), contains("gene_id"), contains("transcript_id"), contains("protein_id"), contains("gene_biotype"), contains("transcript_biotype"), contains("panel"), contains("retirement_status"), contains("release_last_seen")) %>%
                        .[mixedorder(.$hgnc_stable_variant_ID), ] %>%
                        dplyr::rename_all(function(x) {x %>% stringr::str_to_sentence() %>% gsub(pattern = "\\_", replacement = " ") %>% return}) %>%
                        dplyr::mutate("id" = 1:nrow(.), .before = 1) %>%
                        return},
                options = list(fixedHeader = TRUE, lengthMenu = list(c(25, 50, 100, -1), c("25", "50", "100", "All")))
            )
            
            # Create the button to download the scatterplot as PDF
            output$workshop_download_plot <- downloadHandler(
                filename = function() {
                    paste('EDN_workshop_', Sys.Date(), ".", "A" %>% (function(x) {options(digits.secs = 9); Sys.time() %>% as.numeric %>% return}), '.pdf', sep = "")
                },
                content = function(file) {
                    
                    ggplot_final_plot <- ggplot_final_plot
                    
                    ggsave(file, ggplot_final_plot, width = 20, height = 20*(plot_height/plot_width), dpi = 600, units = "cm")
                }
            )
            
        } else {
            
        } 
        
    } )
    
    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observeEvent(input$workshop_plot_output_dblclick, {
        
        # print("input$workshop_plot_output_dblclick")
        # print(input$workshop_plot_output_dblclick)
        # test_workshop_plot_output_dblclick <<- input$workshop_plot_output_dblclick
        
        # print("input$workshop_plot_output_brush")
        # print(input$workshop_plot_output_brush %>% head)
        brush <- input$workshop_plot_output_brush
        
        global_brush <<- brush
        
        if (!is.null(brush)) {
            workshop_plot_brush_ranges$x <- c(brush$xmin, brush$xmax)
            workshop_reactiveValues_plot_metadata$list_y_axis_scale[[brush$panelvar1]] <- workshop_reactiveValues_plot_metadata$list_y_axis_scale[[brush$panelvar]]  %>% .[(brush$ymin %>% round(0)):(brush$ymax %>% round(0))] %>% na.omit
            
            workshop_plot_brush_ranges$logical_brush_zoom_on <- TRUE
            # workshop_plot_brush_ranges$y <- c(workshop_reactive_final_plot() %>% .$plot_view_initial_y_start,
            #                                   workshop_reactive_final_plot() %>% .$plot_view_initial_y_end)
            # workshop_plot_brush_ranges$y <- c(brush$ymin, brush$ymax)          
        } else {
            workshop_plot_brush_ranges$x <- c(workshop_reactive_final_plot() %>% .$plot_view_initial_x_start, 
                                              workshop_reactive_final_plot() %>% .$plot_view_initial_x_end)
            workshop_reactiveValues_plot_metadata$list_y_axis_scale <- workshop_reactiveValues_plot_metadata$list_y_axis_scale_initial
            
            workshop_plot_brush_ranges$logical_brush_zoom_on <- FALSE
            
            # workshop_plot_brush_ranges$y <- c(workshop_reactive_final_plot() %>% .$plot_view_initial_y_start,
            #                                   workshop_reactive_final_plot() %>% .$plot_view_initial_y_end)
        }
        
        print("workshop_plot_brush_ranges$logical_brush_zoom_on")
        print(workshop_plot_brush_ranges$logical_brush_zoom_on)
        
        
        # print("workshop_reactiveValues_plot_metadata$list_y_axis_scale_initial")
        # print(workshop_reactiveValues_plot_metadata$list_y_axis_scale_initial)
        # print("workshop_reactiveValues_plot_metadata$list_y_axis_scale")
        # print(workshop_reactiveValues_plot_metadata$list_y_axis_scale)
        
    }, ignoreNULL = TRUE, ignoreInit = TRUE )
    
    # When a single-click happens, output the coords which were clicked
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observeEvent(input$workshop_plot_output_sglclick, {
        
        output$workshop_graph_click_info <- renderText( {
            
            paste("x = chr", workshop_reactiveValues_current_plot_range$chr, ":", input$workshop_plot_output_sglclick$x %>% round(digits = 0), "\ny = ", input$workshop_plot_output_sglclick$domain$discrete_limits$y %>% .[[input$workshop_plot_output_sglclick$y %>% round(digits = 0)]], 
                  sep = "")
            
        } )
        
        # global_workshop_plot_output_sglclick <<- input$workshop_plot_output_sglclick
        
    }, ignoreNULL = TRUE, ignoreInit = TRUE )
    
    # test <- ggplot() +
    #     ggplot2::facet_grid(panel ~ ., scales = "free_y") +
    #     geom_bar(data = mtcars %>% as_tibble %>% .[, "mpg"] %>% tibble::add_column("panel" = "one"), aes(x = mpg)) +
    #     geom_bar(data = mtcars %>% as_tibble %>% .[, "disp"] %>% tibble::add_column("panel" = "one"), aes(x = disp)) +
    #     geom_bar(data = mtcars %>% as_tibble %>% .[, "drat"] %>% tibble::add_column("panel" = "two"), aes(x = drat)) +
    #     ggh4x::facetted_pos_scales(
    #         y = list(NULL, NULL)
    #     )
    
    
    
    # END WORKSHOP ###
    
    outputOptions(output, "automator_reactive_UI_1", suspendWhenHidden = TRUE)
    outputOptions(output, "automator_reactive_UI_2", suspendWhenHidden = TRUE)
    outputOptions(output, "workshop_reactive_UI_1", suspendWhenHidden = TRUE)
    outputOptions(output, "workshop_reactive_UI_2", suspendWhenHidden = TRUE)
    outputOptions(output, "workshop_reactive_exon_table", suspendWhenHidden = TRUE)
    
} # server

# Run the application 
shinyApp(ui = ui, server = server)
