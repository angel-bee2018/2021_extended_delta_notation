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
            .[which(.$type == return_type), ]
        
    } else if (query_strand == "+" | query_strand == "-") {
        
        # cat("b\n")
        
        # +/- 1 nt tolerance
        tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws &
                                                                 tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
            .[which(.$start > ((query_start %>% as.numeric) - 1 + left_query_shift - left_tolerance) & .$end < ((query_end %>% as.numeric) + 1 + right_query_shift + right_tolerance)), ] %>% 
            .[which((.$start < ((query_start %>% as.numeric) + 1 + left_query_shift + left_tolerance) & .$end > ((query_end %>% as.numeric) - 1 + right_query_shift - right_tolerance))), ] %>% 
            .[which(.$type == return_type), ]
        
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

extract_overlapping_features <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 0, right_tolerance = 0, return_type) {
    
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
            .[which(.$end >= ((query_start %>% as.numeric) + left_query_shift - left_tolerance) & .$start <= ((query_end %>% as.numeric) + right_query_shift + right_tolerance)), ] %>% 
            .[which(.$type %in% return_type), ]
        
    } else if (query_strand == "+" | query_strand == "-") {
        
        # +/- 1 nt tolerance
        tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws &
                                                                 tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
            .[which(.$end >= ((query_start %>% as.numeric) + left_query_shift - left_tolerance) & .$start <= ((query_end %>% as.numeric) + right_query_shift + right_tolerance)), ] %>% 
            .[which(.$type %in% return_type), ]
        
    }
    
    return(tibble_gtf_subset_matching_exons)
    
}

## END extract_overlapping_features() ###

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
    
    # VSR_coordinates <- VSR_coordinates
    # tibble_VSR_exon_start_end <- list_tibble_exon_start_end_per_LIS[[1]]
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
        
        final_VSR_nomenclature <- paste(variant_slot, " ", tibble_sorted_combined_nomenclature$slots %>% paste(collapse = " "), sep = "") %>% 
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
VSR_LIS_organise_exon_rematching <- function(VSR_coordinates, list_tibble_exon_start_end_per_LIS, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1) {
    
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
    # list_tibble_exon_start_end_per_LIS <- list_of_exon_start_end_tibbles
    # tibble_gtf_table <- tibble_ref_gtf
    # 
    # left_query_shift <- 0
    # right_query_shift <- 0
    # left_tolerance <- 1
    # right_tolerance <- 1
    
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
    
    # DETECT VSR OR LIS
    ## if the list has more than one element, then it's a VSR.
    flag_is_LIS <- list_tibble_exon_start_end_per_LIS %>% length == 1
    
    # LOWEST VARIANT DETECTION FOR EACH LIS
    list_initial_lowest_variant_detection <- list_tibble_exon_start_end_per_LIS %>% 
        purrr::map(
            .f = ~VSR_select_reference_transcript_variant(VSR_coordinates = VSR_coordinates,
                                                          tibble_VSR_exon_start_end = .x,
                                                          tibble_gtf_table = tibble_gtf_table,
                                                          left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance))
    
    # get lowest common variant ID
    lowest_common_variant_ID <- purrr::map(.x = list_initial_lowest_variant_detection, .f = ~.x$selected_hgnc_variant_name) %>% unlist %>% mixedsort %>% .[1]
    
    # edit the initial detection lists to force the lowest common variant ID.
    list_initial_lowest_variant_detection_force_lowest_variant_ID <- purrr::map(
        .x = list_initial_lowest_variant_detection,
        .f = function(a1) {
            
            list_edited <- a1
            tibble_parent_gene_entries <- list_edited$tibble_parent_gene_entries
            
            # rename the variant ID flag
            list_edited$selected_hgnc_variant_name <- lowest_common_variant_ID
            # refresh the selected parent entries for the new variant ID
            list_edited$tibble_selected_transcript_entries <- tibble_parent_gene_entries[which(tibble_parent_gene_entries$hgnc_stable_variant_ID == lowest_common_variant_ID), ]
            
            # mark the tibble of body exon info to indicate which entries overlap with a different variant ID and need to be renamed
            list_edited$tibble_body_exon_info <- list_edited$tibble_body_exon_info %>% 
                add_column("still_needs_renaming" = list_edited$tibble_body_exon_info$lowest_stable_variant_ID_matched_to_exon != lowest_common_variant_ID)
            
            return(list_edited)
            
        } )
    
    # first pass naming using the lowest variant ID for ALL exons.
    list_first_pass_naming <- list_initial_lowest_variant_detection_force_lowest_variant_ID %>% purrr::map(.f = ~VSR_name_body_exons_and_VSR(list_output_from_VSR_select_reference_transcript_variant = .x, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance))
    
    # check if we need second-pass naming
    if (any(list_first_pass_naming %>% purrr::map(.f = ~.x$tibble_body_exon_info$still_needs_renaming) %>% unlist %>% na.omit == TRUE)) {
        
        # second-pass naming. inexact over-ride naming of each exon using the lowest overlapped variant ID.
        list_second_pass_naming <- purrr::imap(
            .x = list_first_pass_naming,
            .f = function(a1, a2) {
                
                # DEBUG ###
                # a1 <- list_first_pass_naming[[3]]
                ###########
                
                # cat(a2, "\n")
                
                # retrieve the row indices of the exons that still need naming.
                ## note that we do not match for A3/5SS exons, as the information is already captured in the VSR.
                row.indices_exons_still_need_renaming <- which(a1$tibble_body_exon_info$still_needs_renaming == TRUE &
                                                                   ((a1$tibble_body_exon_info$left_end_of_VSR == TRUE | a1$tibble_body_exon_info$right_end_of_VSR == TRUE) & a1$tibble_body_exon_info$exon_matches_to_any_reference == FALSE) == FALSE)
                
                # map thru all the row indices of exons that still need renaming
                list_renamed_exons <- purrr::map(.x = a1$tibble_body_exon_info %>% .[row.indices_exons_still_need_renaming, ] %>% purrr::array_tree(),
                                                 .f = ~name_a_single_exon(query_chr = a1$query_chr, query_start = .x$start, query_end = .x$end, query_strand = "*", tibble_gtf_table = a1$tibble_parent_gene_entries, variant_ID_override = .x$lowest_stable_variant_ID_matched_to_exon, override_mode = "inexact", left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance))
                
                list_edited <- a1
                
                list_edited[["list_body_exon_nomenclature"]][row.indices_exons_still_need_renaming] <- list_renamed_exons %>% purrr::map(~paste("(", .x$variant_ID_slot, ".", a1$tibble_parent_gene_entries %>% .[.$hgnc_stable_variant_ID == .x$variant_ID_slot, ] %>% .$transcript_version %>% na.omit %>% unique, " ", .x$exon_slot, ")", sep = ""))
                
                # don't need renaming anymore! :)
                list_edited$tibble_body_exon_info$still_needs_renaming <- FALSE
                
                # done
                
                return(list_edited)
                
            } )
        
    } else {
        
        list_second_pass_naming <- list_first_pass_naming
        
    }
    
    # RETURN NOMENCLATURE
    ## collapse by spaces between exons of the same LIS
    list_final_LIS_nomenclature <- purrr::map(
        .x = list_second_pass_naming,
        .f = function(a1) {
            
            # DEBUG ###
            # a1 <- list_second_pass_naming[[1]]
            ###########
            
            collapsed_nomenclature_for_one_LIS <- a1$list_body_exon_nomenclature %>% unlist %>% na.omit %>% paste(collapse = " ")
            
            return(collapsed_nomenclature_for_one_LIS)
            
        } )
    
    # finally, extract the transcript version
    variant_slot <- paste(lowest_common_variant_ID, ".", tibble_gtf_table[tibble_gtf_table$hgnc_stable_variant_ID == lowest_common_variant_ID, "transcript_version"] %>% unlist %>% na.omit %>% unique %>% .[1], sep = "")
    
    ## collapse by "/" between different LISs and sandwich between the VSR upstream and downstream slots
    final_VSR_nomenclature <- paste(variant_slot, " ", list_first_pass_naming[[1]]$upstream_VSR_slot, " (", list_final_LIS_nomenclature %>% paste(collapse = " / "), ") ", list_first_pass_naming[[1]]$downstream_VSR_slot, sep = "") %>% 
        trimws %>%
        gsub(pattern = "\\(\\((.*)\\)\\)", replacement = "(\\1)") %>% 
        gsub(pattern = "\\(\\)", replacement = "") 
    
    return(final_VSR_nomenclature)
    
}

# END VSR_LIS_organise_exon_rematching() ###

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
        
        selected_transcript_version <- tibble_exonic_overlapped_transcripts_parent_entries[tibble_exonic_overlapped_transcripts_parent_entries$transcript_id == selected_transcript_id, "transcript_version"] %>% unlist %>% unique
        
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
                            
                            return(paste("(", tibble_flanking_matches_in_another_transcript$hgnc_stable_variant_ID %>% unique, ".", tibble_flanking_matches_in_another_transcript$transcript_version %>% unique, " E", tibble_flanking_matches_in_another_transcript$exon_number %>% min, "^E", tibble_flanking_matches_in_another_transcript$exon_number %>% max, ")", sep = ""))
                            
                        } else {
                            
                            return(paste("(", tibble_flanking_matches_in_another_transcript$hgnc_stable_variant_ID %>% unique, ".", tibble_flanking_matches_in_another_transcript$transcript_version %>% unique, " E", tibble_flanking_matches_in_another_transcript$exon_number %>% max, "^E", tibble_flanking_matches_in_another_transcript$exon_number %>% min, ")", sep = ""))
                            
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
                            
                            start_slot <- paste("(", tibble_junction_start_selected_exon$hgnc_stable_variant_ID, ".", tibble_junction_start_selected_exon$transcript_version, " E", tibble_junction_start_selected_exon$exon_number, ")", sep = "")
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
                            
                            end_slot <- paste("(", tibble_junction_end_selected_exon$hgnc_stable_variant_ID, ".", tibble_junction_end_selected_exon$transcript_version, " E", tibble_junction_end_selected_exon$exon_number, ")", sep = "")
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
                return(paste(selected_hgnc_variant_name, ".", selected_transcript_version, " ", list_junction_nomenclature_slots %>% unlist %>% rev %>% paste(collapse = " / "), sep = ""))
            } else {
                return(paste(selected_hgnc_variant_name, ".", selected_transcript_version, " ", list_junction_nomenclature_slots %>% unlist %>% paste(collapse = " / "), sep = ""))
            }
            
        } else if (single_junction == FALSE) {
            
            if (LSV_strand == "-") {
                return(paste(selected_hgnc_variant_name, ".", selected_transcript_version, " (", list_junction_nomenclature_slots %>% unlist %>% rev %>% paste(collapse = " / "), ")", sep = ""))
            } else {
                return(paste(selected_hgnc_variant_name, ".", selected_transcript_version, " (", list_junction_nomenclature_slots %>% unlist %>% paste(collapse = " / "), ")", sep = ""))
            }
            
        }
        
    } else if (flag_is_intergenic == TRUE) {
        
        return("Intergenic")
        
    }
    
}

# END LSV_AJ_organise_junction_matching() ###

# AE: FUNCTION TO NAME A SINGLE EXON
# we will take an exon's coordinates and find overlapping counterparts in the reference. Then add exon modifiers as required.
# alternatively, a splice variant over-ride can be given. This is used as part of the VSR pipeline.
# override_mode - "exact": specifically look for desired variant ID for all matching. "inexact": only look for the desired variant ID if exact matching failed.
name_a_single_exon <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, variant_ID_override = NULL, override_mode = "inexact", left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1) {
    
    # DEBUG ###
    # query_chr <- AE_query_chr
    # query_start <- AE_query_start %>% type.convert
    # query_end <- AE_query_end %>% type.convert
    # query_strand <- AE_query_strand
    # tibble_gtf_table <- tibble_ref_gtf
    # variant_ID_override <- NULL
    # override_mode <- "inexact"
    # 
    # left_query_shift <- 0
    # right_query_shift <- 0
    # left_tolerance <- 0
    # right_tolerance <- 0
    
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
    
    # EXACT MATCH - direct naming of the reference exon numbers
    ## look for exact match in the reference
    tibble_matched_reference_exons <- extract_matching.exons(query_chr = query_chr, query_start = query_start, query_end = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon")
    
    # check if the override variant ID was in the exact match list. if so, then use it.
    if (is.null(variant_ID_override) == FALSE & override_mode == "exact") {
        
        tibble_matched_reference_exons <- tibble_matched_reference_exons[tibble_matched_reference_exons$hgnc_stable_variant_ID == variant_ID_override, ]
        
        
    }
    
    if (tibble_matched_reference_exons %>% nrow > 0) {
        
        ## select lowest variant ID
        selected_hgnc_variant_name <- tibble_matched_reference_exons$hgnc_stable_variant_ID %>% mixedsort %>% .[1]
        ## return matched exon number
        matched_exon_number <- tibble_matched_reference_exons[tibble_matched_reference_exons$hgnc_stable_variant_ID == selected_hgnc_variant_name, "exon_number"] %>% paste %>% type.convert
        
        return(list(
            "variant_ID_slot" = selected_hgnc_variant_name,
            "exon_slot" = paste("E", matched_exon_number, sep = "")
        ))
        
        # INEXACT MATCH - exon modifiers required
    } else {
        
        tibble_overlapping_reference_exon <- extract_overlapping_features(query_chr = query_chr, query_start = query_start, query_end = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon")
        
        # check if the override variant ID was in the inexact match list. if so, then use it.
        if (is.null(variant_ID_override) == FALSE) {
            
            tibble_overlapping_reference_exon <- tibble_overlapping_reference_exon[tibble_overlapping_reference_exon$hgnc_stable_variant_ID == variant_ID_override, ]
            
            
        }
        
        # if exons dont overlap then last resort is to overlap with reference transcripts.
        # this can occur if the exon is entirely intronic.
        if (tibble_overlapping_reference_exon %>% nrow == 0) {
            
            tibble_overlapping_reference_exon <- extract_overlapping_features(query_chr = query_chr, query_start = query_start, query_end = query_end, query_strand = query_strand, tibble_gtf_table = tibble_gtf_table, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "transcript")
            
            # check if the override variant ID was in the inexact match list. if so, then use it.
            if (is.null(variant_ID_override) == FALSE) {
                
                # extract all the exonic entries from the transcript_id.
                tibble_overlapping_reference_exon <- tibble_gtf_table[tibble_gtf_table$hgnc_stable_variant_ID == variant_ID_override, ]
                
            }
            
        }
        
        # if still no match then this is intergenic. report genomic coordinates.
        if (tibble_overlapping_reference_exon %>% nrow == 0) {
            
            return(list(
                "variant_ID_slot" = "Intergenic",
                "exon_slot" = paste(query_chr, ":", query_start, "-", query_end, sep = "")))
            
        } else {
            
            ## select lowest variant ID
            selected_hgnc_variant_name <- tibble_overlapping_reference_exon$hgnc_stable_variant_ID %>% mixedsort %>% .[1]
            
            ## filter the overlapped reference exon table
            tibble_overlapping_reference_exon_filtered_for_selected <- tibble_overlapping_reference_exon[tibble_overlapping_reference_exon$hgnc_stable_variant_ID == selected_hgnc_variant_name, ]
            
            exon_numbers_overlapped_by_query_start <- extract_overlapping_features(query_chr = query_chr, query_start = query_start, query_end = query_start, query_strand = "*", tibble_gtf_table = tibble_overlapping_reference_exon_filtered_for_selected, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon") %>% .$exon_number
            
            exon_numbers_overlapped_by_query_end <- extract_overlapping_features(query_chr = query_chr, query_start = query_end, query_end = query_end, query_strand = "*", tibble_gtf_table = tibble_overlapping_reference_exon_filtered_for_selected, left_query_shift = left_query_shift, right_query_shift = right_query_shift, left_tolerance = left_tolerance, right_tolerance = right_tolerance, return_type = "exon") %>% .$exon_number
            
            ## retrieve all exonic entries for the selected transcript
            tibble_selected_exonic_entries <- tibble_gtf_table[tibble_gtf_table$hgnc_stable_variant_ID == selected_hgnc_variant_name & tibble_gtf_table$type == "exon", ]
            
            ## pool together all selected vertices
            vector_selected_transcript_vertices <- tibble_selected_exonic_entries[, c("start", "end")] %>% unlist %>% unique
            
            ## find any overlap between query and the selected vertices
            selected_vertex_matched_to_query_start <- vector_selected_transcript_vertices[vector_selected_transcript_vertices <= (query_start + left_query_shift + left_tolerance) & vector_selected_transcript_vertices >= (query_start + left_query_shift - left_tolerance)] 
            selected_vertex_matched_to_query_end <- vector_selected_transcript_vertices[vector_selected_transcript_vertices <= (query_end + left_query_shift + left_tolerance) & vector_selected_transcript_vertices >= (query_end + left_query_shift - left_tolerance)] 
            
            ## find the closest upstream/downstream exon vertices
            closest_LHS_reference_vertex_to_query_start <- vector_selected_transcript_vertices[vector_selected_transcript_vertices <= query_start] %>% max
            closest_RHS_reference_vertex_to_query_start = vector_selected_transcript_vertices[vector_selected_transcript_vertices >= query_start] %>% min
            
            closest_LHS_reference_vertex_to_query_end = vector_selected_transcript_vertices[vector_selected_transcript_vertices <= query_end] %>% max
            closest_RHS_reference_vertex_to_query_end = vector_selected_transcript_vertices[vector_selected_transcript_vertices >= query_end] %>% min
            
            # CREATE EXON MODIFIERS
            # MIDDLE SLOT
            ## check how many selected exons are covered by the query exon
            if (tibble_overlapping_reference_exon_filtered_for_selected %>% nrow == 1) {
                middle_slot <- paste("E", tibble_overlapping_reference_exon_filtered_for_selected$exon_number, sep = "")
            } else if (tibble_overlapping_reference_exon_filtered_for_selected %>% nrow > 1) {
                middle_slot <- paste("E", tibble_overlapping_reference_exon_filtered_for_selected$exon_number %>% min, "_E", tibble_overlapping_reference_exon_filtered_for_selected$exon_number %>% max, sep = "")
            }
            
            # LEFT SLOT (QUERY EXON START)
            ## check if it matched any of the selected vertices
            ## if it matched to any of the selected vertices, then the slot is empty.
            ## if not, then check whether it's in an intronic or exonic region and call truncation/extension accordingly
            if (selected_vertex_matched_to_query_start %>% length > 0) {
                
                left_slot <- ""
                
            } else if (selected_vertex_matched_to_query_start %>% length == 0) {
                
                # INTRONIC QUERY VERTEX - LEFT EXTENSION
                if (exon_numbers_overlapped_by_query_start %>% length == 0) {
                    
                    left_slot <- paste("+", closest_RHS_reference_vertex_to_query_start - query_start, sep = "")
                    
                    # EXONIC QUERY VERTEX - LEFT TRUNCATION   
                } else if (exon_numbers_overlapped_by_query_start %>% length > 0) {
                    
                    left_slot <- paste("–", query_start - closest_LHS_reference_vertex_to_query_start, sep = "")
                    
                }
                
            }
            
            # RIGHT SLOT (QUERY EXON END)
            ## check if it matched any of the selected vertices
            ## if it matched to any of the selected vertices, then the slot is empty.
            ## if not, then check whether it's in an intronic or exonic region and call truncation/extension accordingly
            if (selected_vertex_matched_to_query_end %>% length > 0) {
                
                right_slot <- ""
                
            } else if (selected_vertex_matched_to_query_end %>% length == 0) {
                
                # INTRONIC QUERY VERTEX - RIGHT EXTENSION
                if (exon_numbers_overlapped_by_query_end %>% length == 0) {
                    
                    right_slot <- paste("+", query_end - closest_LHS_reference_vertex_to_query_end, sep = "")
                    
                    # EXONIC QUERY VERTEX - RIGHT TRUNCATION   
                } else if (exon_numbers_overlapped_by_query_end %>% length > 0) {
                    
                    right_slot <- paste("–", closest_RHS_reference_vertex_to_query_end - query_end, sep = "")
                    
                }
                
            }
            
            # JOIN SLOTS FOR THE FINAL NAME MODULE
            
            if (query_strand == "-") {
                
                return(list(
                    "variant_ID_slot" = selected_hgnc_variant_name,
                    "exon_slot" = paste(right_slot, middle_slot, left_slot, sep = "")))
                
            } else {
                
                return(list(
                    "variant_ID_slot" = selected_hgnc_variant_name,
                    "exon_slot" = paste(left_slot, middle_slot, right_slot, sep = "")))
                
            }
            
        }
        
    }
    
}

# END name_a_single_exon() ###

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
        
        # tibble_ref_gtf <- data.table::fread(file = "/mnt/LTS/projects/2020_isoform_nomenclature/nomenclature_app/app_native/EDN_workshop/data/annotated_ensembl_gtf_release_102.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
        
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
            # 
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
                        return(tibble("start" = gsub(x = vector_exon_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2"),
                                      "end" = gsub(x = vector_exon_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3")) %>% type_convert)
                        
                    } )
                
                paste("Suggested shorthand notation: \n", VSR_LIS_organise_exon_rematching(VSR_coordinates = automator_input_alternative_event_region, list_tibble_exon_start_end_per_LIS = list_of_exon_start_end_tibbles, tibble_gtf_table = tibble_ref_gtf, left_query_shift = automator_input_left_query_end_shift, right_query_shift = automator_input_right_query_end_shift, left_tolerance = automator_input_left_match_tolerance, right_tolerance = automator_input_right_match_tolerance), "\n", sep = "") 
                
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
                
                paste("Suggested shorthand notation: \n", VSR_LIS_organise_exon_rematching(VSR_coordinates = automator_input_alternative_event_region, list_tibble_exon_start_end_per_LIS = list_of_exon_start_end_tibbles, tibble_gtf_table = tibble_ref_gtf, left_query_shift = automator_input_left_query_end_shift, right_query_shift = automator_input_right_query_end_shift, left_tolerance = automator_input_left_match_tolerance, right_tolerance = automator_input_right_match_tolerance), "\n", sep = "")
                
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
        
        # tibble_ref_gtf <- data.table::fread(file = "/mnt/LTS/projects/2020_isoform_nomenclature/nomenclature_app/app_native/EDN_workshop/data/annotated_ensembl_gtf_release_102.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
        
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
                    input_annotation_tibble
                )
                names(workshop_reactiveValues_annotation_files$annotation_files$reference_gtf)[length(workshop_reactiveValues_annotation_files$annotation_files$reference_gtf)] <- workshop_reactiveValues_custom_file_import_name$workshop_custom_file_import_name
                # automatically checkbox the added GTF
                workshop_reactiveValues_annotation_files_selected$vector_checked_items <- c(workshop_reactiveValues_annotation_files_selected$vector_checked_items, paste("reference_gtf", "|", workshop_reactiveValues_custom_file_import_name$workshop_custom_file_import_name, sep = ""))
            } else if (input$workshop_import_file_type_selection == "Upload custom GTF" & is.null(input$workshop_path_to_custom_gtf$datapath) == FALSE) {
                workshop_reactiveValues_annotation_files$annotation_files$custom_gtf <- workshop_reactiveValues_annotation_files$annotation_files$custom_gtf %>% purrr::splice(
                    input_annotation_tibble
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
        
        triage_result <- triage_input_coordinates(vector_input_coordinates = input$workshop_input_range, vector_of_expected_chromosomes = c(1:22, "X", "Y", "M"), expect_stranded = TRUE)
        
        if (triage_result == "triage fail") {
            
            triage_result_2 <- triage_input_coordinates(vector_input_coordinates = input$workshop_input_range, vector_of_expected_chromosomes = c(1:22, "X", "Y", "M"), expect_stranded = FALSE)
            
            if (triage_result_2 == "triage successful") {
                
                input_chr <- gsub(x = input$workshop_input_range, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)", replacement = "\\1")
                input_start <- gsub(x = input$workshop_input_range, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)", replacement = "\\2") %>% type.convert
                input_end <- gsub(x = input$workshop_input_range, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)", replacement = "\\3") %>% type.convert
                input_strand <- gsub(x = input$workshop_input_range, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)([\\:])?", replacement = "\\4") %>% gsub(pattern = "\\:", replacement = "")
                
                workshop_reactiveValues_user_ranges$id <- c(workshop_reactiveValues_user_ranges$id, if (length(workshop_reactiveValues_user_ranges$id) == 0) {"1"} else {as.character(max(workshop_reactiveValues_user_ranges$id %>% as.numeric) + 1)} )
                workshop_reactiveValues_user_ranges$chr <- c(workshop_reactiveValues_user_ranges$chr, input_chr)
                workshop_reactiveValues_user_ranges$start <- c(workshop_reactiveValues_user_ranges$start, input_start)
                workshop_reactiveValues_user_ranges$end <- c(workshop_reactiveValues_user_ranges$end, input_end)
                workshop_reactiveValues_user_ranges$strand <- c(workshop_reactiveValues_user_ranges$strand, input_strand)
                workshop_reactiveValues_user_ranges$range_type <- c(workshop_reactiveValues_user_ranges$range_type, input$workshop_range_type)
                workshop_reactiveValues_user_ranges$panel <- c(workshop_reactiveValues_user_ranges$panel, "user_ranges")
                
                triage_result <- triage_result_2
                
            }
            
        } else if (triage_result == "triage successful") {
            
            input_chr <- gsub(x = input$workshop_input_range, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
            input_start <- gsub(x = input$workshop_input_range, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2") %>% type.convert
            input_end <- gsub(x = input$workshop_input_range, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3") %>% type.convert
            input_strand <- gsub(x = input$workshop_input_range, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
            
            workshop_reactiveValues_user_ranges$id <- c(workshop_reactiveValues_user_ranges$id, if (length(workshop_reactiveValues_user_ranges$id) == 0) {"1"} else {as.character(max(workshop_reactiveValues_user_ranges$id %>% as.numeric) + 1)} )
            workshop_reactiveValues_user_ranges$chr <- c(workshop_reactiveValues_user_ranges$chr, input_chr)
            workshop_reactiveValues_user_ranges$start <- c(workshop_reactiveValues_user_ranges$start, input_start)
            workshop_reactiveValues_user_ranges$end <- c(workshop_reactiveValues_user_ranges$end, input_end)
            workshop_reactiveValues_user_ranges$strand <- c(workshop_reactiveValues_user_ranges$strand, input_strand)
            workshop_reactiveValues_user_ranges$range_type <- c(workshop_reactiveValues_user_ranges$range_type, input$workshop_range_type)
            workshop_reactiveValues_user_ranges$panel <- c(workshop_reactiveValues_user_ranges$panel, "user_ranges")
            
        }
        
        output$workshop_nomenclature_output <- renderText({triage_result})
        
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
        
        # check input coords, double triage
        triage_result <- triage_input_coordinates(vector_input_coordinates = input$workshop_jump_to_coords, vector_of_expected_chromosomes = c(1:22, "X", "Y", "M"), expect_stranded = TRUE)
        
        if (triage_result == "triage fail") {
            
            triage_result_2 <- triage_input_coordinates(vector_input_coordinates = input$workshop_jump_to_coords, vector_of_expected_chromosomes = c(1:22, "X", "Y", "M"), expect_stranded = FALSE)
            
            if (triage_result_2 == "triage successful") {
                
                input_chr <- gsub(x = input$workshop_jump_to_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)", replacement = "\\1")
                input_start <- gsub(x = input$workshop_jump_to_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)", replacement = "\\2") %>% type.convert
                input_end <- gsub(x = input$workshop_jump_to_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)", replacement = "\\3") %>% type.convert
                input_strand <- "*"
                
                workshop_reactiveValues_current_plot_range$chr <- input_chr
                workshop_reactiveValues_current_plot_range$start <- input_start
                workshop_reactiveValues_current_plot_range$end <- input_end
                
                triage_result <- triage_result_2
                
            }
            
        } else if (triage_result == "triage successful") {
            
            input_chr <- gsub(x = input$workshop_jump_to_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
            input_start <- gsub(x = input$workshop_jump_to_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2") %>% type.convert
            input_end <- gsub(x = input$workshop_jump_to_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3") %>% type.convert
            input_strand <- gsub(x = input$workshop_jump_to_coords, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
            
            workshop_reactiveValues_current_plot_range$chr <- input_chr
            workshop_reactiveValues_current_plot_range$start <- input_start
            workshop_reactiveValues_current_plot_range$end <- input_end
            
        }
        
        output$workshop_nomenclature_output <- renderText({triage_result})
        
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
        
        global_workshop_reactiveValues_current_plot_range <<- reactiveValuesToList(workshop_reactiveValues_current_plot_range)
        
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
        
        print("tibble_all_user_ranges")
        print(tibble_all_user_ranges)
        
        global_workshop_reactiveValues_selected_user_range <<- reactiveValuesToList(workshop_reactiveValues_selected_user_range)
        
        # set up the ranges that the user has selected to highlight and calculate distances for
        selected_user_range_chr <- workshop_reactiveValues_selected_user_range$chr
        selected_user_range_start <- workshop_reactiveValues_selected_user_range$start
        selected_user_range_end <- workshop_reactiveValues_selected_user_range$end
        selected_user_range_strand <- workshop_reactiveValues_selected_user_range$strand

        if (selected_user_range_strand == "*") {
            selected_user_range_strand <- c("+", "-")
        }
        
        # print("workshop_reactiveValues_current_plot_range$chr")
        # print(workshop_reactiveValues_current_plot_range$chr)
        # print("workshop_reactiveValues_current_plot_range$start")
        # print(workshop_reactiveValues_current_plot_range$start)
        # print("workshop_reactiveValues_current_plot_range$end")
        # print(workshop_reactiveValues_current_plot_range$end)
        
        print("plot_view_initial_x_end beginning")
        print(workshop_reactiveValues_current_plot_range$end)
        
        print("plot_view_initial_x_start beginning")
        print(workshop_reactiveValues_current_plot_range$start)
        
        # subset all tables to be plotted, for user-specified range (1.5x jump to/user range selection)
        plot_view_initial_x_start0 <- workshop_reactiveValues_current_plot_range$start - 1.5*(workshop_reactiveValues_current_plot_range$end - workshop_reactiveValues_current_plot_range$start) + (plot_x_scale^3)*(workshop_reactiveValues_current_plot_range$end - workshop_reactiveValues_current_plot_range$start)
        plot_view_initial_x_end0 <- workshop_reactiveValues_current_plot_range$end + 1.5*(workshop_reactiveValues_current_plot_range$end - workshop_reactiveValues_current_plot_range$start) - (plot_x_scale^3)*(workshop_reactiveValues_current_plot_range$end - workshop_reactiveValues_current_plot_range$start)
        
        print("plot_view_initial_x_end middle")
        print(plot_view_initial_x_end0)
        
        print("plot_view_initial_x_start middle")
        print(plot_view_initial_x_start0)
        
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
        print("plot_view_initial_x_end final")
        print(plot_view_initial_x_end)

        print("plot_view_initial_x_start final")
        print(plot_view_initial_x_start)
        
        global_workshop_reactiveValues_annotation_files_selected <<- reactiveValuesToList(workshop_reactiveValues_annotation_files_selected)
        
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
        
        names(list_tibbles_track_features_visible_flattened) <- purrr::map2(
            .x = list_tibbles_track_features_visible, 
            .y = names(list_tibbles_track_features_visible), 
            .f = function(a1, a2) {
                
                purrr::map2(
                    .x = a1,
                    .y = names(a1),
                    .f = function(b1, b2) {
                        
                        return(paste(a2 %>% stringr::str_to_sentence() %>% gsub(pattern = "gtf", replacement = "GTF") %>% gsub(pattern = "\\_", replacement = " "), ": ", b2, sep = ""))
                        
                    } ) %>% unlist
                
            } ) %>% unlist
        
        list_tibbles_track_features_visible_flattened <- purrr::map2(
            .x = list_tibbles_track_features_visible_flattened,
            .y = names(list_tibbles_track_features_visible_flattened),
            .f = ~.x %>% tibble::add_column("panel" = .y)
        )
        
        ## find the user ranges visible in the viewing range
        tibble_user_ranges_visible <- tibble_all_user_ranges[which(tibble_all_user_ranges$chr == workshop_reactiveValues_current_plot_range$chr & tibble_all_user_ranges$start <= plot_view_initial_x_end & tibble_all_user_ranges$end >= plot_view_initial_x_start), ]
        
        global_list_tibbles_track_features_visible_flattened <<- list_tibbles_track_features_visible_flattened
        
        # distance shenanigans
        ## calculate distances for selected range
        ## only do this if the selected range is visible
        if (workshop_reactiveValues_selected_user_range$id %in% tibble_user_ranges_visible$id & nrow(tibble_user_ranges_visible) > 0 & length(list_tibbles_track_features_visible_flattened) > 0) {
            
            # flatten the original full GTF table
            list_tibbles_track_features_all_flattened <- workshop_reactiveValues_annotation_files_selected$annotation_files %>% flatten
            
            global_list_tibbles_track_features_all_flattened <<- list_tibbles_track_features_all_flattened
            
            ## having determined the plot window, calculate distances to every exon in the plot range
            list_distance_annotation_data_flattened <- purrr::pmap(
                .l = list(
                    "a1" = list_tibbles_track_features_visible_flattened,
                    "a2" = names(list_tibbles_track_features_visible_flattened),
                    "a3" = list_tibbles_track_features_all_flattened
                ),
                .f = function(a1, a2, a3) {
                    
                    # DEBUG ###
                    # a1 <- list_tibbles_track_features_visible_flattened[[1]]
                    # a2 <- names(list_tibbles_track_features_visible_flattened)[[1]]
                    # a3 <- list_tibbles_track_features_all_flattened[[1]]
                    # selected_user_range_chr <- "7"
                    # selected_user_range_start <- 7157427
                    # selected_user_range_end <- 7234302
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
                                left_ref_vertex_grown_from_user_query_start <- vector_all_ref_vertices[vector_all_ref_vertices < selected_user_range_start] %>% max
                                right_ref_vertex_grown_from_user_query_start <- vector_all_ref_vertices[vector_all_ref_vertices > selected_user_range_start] %>% min
                                left_ref_vertex_grown_from_user_query_end <- vector_all_ref_vertices[vector_all_ref_vertices < selected_user_range_end] %>% max
                                right_ref_vertex_grown_from_user_query_end <- vector_all_ref_vertices[vector_all_ref_vertices > selected_user_range_end] %>% min
                                
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
                                
                            } ) %>% dplyr::bind_rows()
                        
                        # make the distances directional
                        tibble_distance_annotations_based_on_user_query <- tibble_distance_annotations_based_on_user_query %>% 
                            dplyr::mutate("ref_vertex_minus_query_vertex" = purrr::map(.x = `ref_vertex_minus_query_vertex`, .f = function(x) {
                                
                                if (x < 0) {
                                    return(paste(abs(x), ">", sep = ""))
                                } else if (x > 0) {
                                    return(paste("<", abs(x), sep = ""))
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
            
            list_tibbles_track_features_visible_flattened <- list()
            
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
            global_list_tibbles_track_features_visible_flattened <<- list_tibbles_track_features_visible_flattened
            
            print("tibble_user_ranges_visible")
            print(tibble_user_ranges_visible)
            
            print("list_distances_between_user_ranges_and_reference_annotations")
            print(list_distances_between_user_ranges_and_reference_annotations)
            ###########
            
            # CREATE GGPLOT
            ggplot_final_plot <- list(
                ggplot(),
                ggplot2::facet_grid(factor(panel, level = workshop_reactiveValues_plot_metadata$list_y_axis_scale %>% names) ~ ., scales = "free_y"),
                # these mark the original viewing window
                geom_vline(colour = "green", lty = 2, xintercept = workshop_reactiveValues_current_plot_range$start),
                geom_vline(colour = "green", lty = 2, xintercept = workshop_reactiveValues_current_plot_range$end)) %>% 
                
                purrr::splice(
                    
                    # reference_gtf
                    if (length(list_tibbles_track_features_visible_flattened[grep(x = names(list_tibbles_track_features_visible_flattened), pattern = "^Reference GTF")]) > 0) {
                        
                        purrr::map(
                            .x = list_tibbles_track_features_visible_flattened[grep(x = names(list_tibbles_track_features_visible_flattened), pattern = "^Reference GTF")], 
                            .f = function(a1) {
                                
                                list(
                                    
                                    geom_segment(data = a1 %>% dplyr::filter(type == "transcript"), colour = "slateblue1", mapping = aes(x = start, xend = end, y = transcript_id, yend = transcript_id)),
                                    geom_text(data = a1 %>% dplyr::filter(type == "transcript"), nudge_y = 0.25, fontface = "italic", mapping = aes(x = mean(workshop_plot_brush_ranges$x), y = transcript_id, label = purrr::pmap(.l = list("b1" = strand, "b2" = hgnc_stable_variant_ID, "b3" = transcript_version), .f = function(b1, b2, b3) {if (b1 == "+") {paste("> > > > > > ", b2, if (is.na(b3) == FALSE) {"."}, if (is.na(b3) == FALSE) {b3}, " > > > > > >", sep = "")} else if (b1 == "-") {paste("< < < < < < ", b2, if (is.na(b3) == FALSE) {"."}, if (is.na(b3) == FALSE) {b3}, " < < < < < <", sep = "")} else {paste(b2, if (is.na(b3) == FALSE) {"."}, if (is.na(b3)) {b3}, sep = "")} } ) %>% unlist)),
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

                                    geom_segment(data = a3, colour = "red", arrow = arrow(angle = 45), mapping = aes(x = ref_vertex, xend = query_vertex, y = transcript_id, yend = transcript_id)),
                                    geom_label(data = a3, colour = "red", nudge_y = -0.25, mapping = aes(x = purrr::map2(.x = ref_vertex, .y = query_vertex, .f = ~c(.x, .y) %>% mean) %>% unlist, y = transcript_id, label = ref_vertex_minus_query_vertex))

                                )
                                
                            } ) %>% purrr::flatten()

                    },

                    # adaptive facet aspect ratio
                    ggh4x::force_panelsizes(rows = workshop_reactiveValues_plot_metadata$vector_number_of_features_per_track %>%
                                                (function(x) {
                                                    if (sum(x) != 0) {
                                                        return(x/sum(x))
                                                        } else {
                                                            return(0)}
                                                    } ) ),

                    # facet-specific brush resizing (y)
                    ggh4x::facetted_pos_scales(
                        y = workshop_reactiveValues_plot_metadata$list_y_axis_scale %>% purrr::map(~scale_y_discrete(limits = .x, breaks = .x, labels = .x))
                    ),
                    # brush resizing (x)
                    # NOTE we CANNOT use coord_cartesion for facet-specific y. MUST use the scales.
                    coord_cartesian(xlim = workshop_plot_brush_ranges$x
                                    # ,
                                    # ylim = c(plot_view_initial_y_start, plot_view_initial_y_end)
                                    ),
                    theme_bw(),
                    theme(text = element_text(family = "Helvetica"))
                )  %>% purrr::reduce(ggplot2:::`+.gg`)
            
            output$workshop_plot_output <- renderPlot( {
                
                ggplot_final_plot
                
            }, height = plot_height, width = plot_width )
            
            # output$workshop_ref_table_output <- renderDataTable(
            #     {workshop_reactive_final_plot() %>% .$list_tibbles_track_features_visible_flattened %>% rbindlist(use.names = TRUE, fill = TRUE) %>%
            #             dplyr::select(contains("hgnc_stable_variant_ID"), contains("transcript_version"), contains("type"), contains("exon_number"), contains("seqnames"), contains("start"), contains("end"), contains("width"), contains("strand"), contains("gene_id"), contains("transcript_id"), contains("protein_id"), contains("gene_biotype"), contains("transcript_biotype"), contains("panel"), contains("retirement_status"), contains("release_last_seen")) %>% 
            #             .[mixedorder(.$hgnc_stable_variant_ID), ] %>%
            #             dplyr::rename_all(function(x) {x %>% stringr::str_to_sentence() %>% gsub(pattern = "\\_", replacement = " ") %>% return}) %>%
            #             dplyr::mutate("id" = 1:nrow(.), .before = 1) %>% 
            #             return}, 
            #     options = list(fixedHeader = TRUE, lengthMenu = list(c(25, 50, 100, -1), c("25", "50", "100", "All")))
            # )
            
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
        # test_workshop_plot_output_brush <<- input$workshop_plot_output_brush
        
        brush <- input$workshop_plot_output_brush
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
