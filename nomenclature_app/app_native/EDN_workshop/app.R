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

library(vroom)
library(tidyverse)
library(rtracklayer)
library(gtools)
library(data.table)

library(furrr)

library(svgPanZoom)
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
            .[which(.$type == return_type), ]
        
    } else if (query_strand == "+" | query_strand == "-") {
        
        # +/- 1 nt tolerance
        tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws &
                                                                 tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
            .[which(.$end >= ((query_start %>% as.numeric) + left_query_shift - left_tolerance) & .$start <= ((query_end %>% as.numeric) + right_query_shift + right_tolerance)), ] %>% 
            .[which(.$type == return_type), ]
        
    }
    
    return(tibble_gtf_subset_matching_exons)
    
}

## END extract_overlapping_features() ###

## VSRs/LIVs: FUNCTION TO FIND THE STABLE VARIANT NUMBER ###
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
    # tibble_VSR_exon_start_end <- list_tibble_exon_start_end_per_LIV[[1]]
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

## VSRs/LIVs: FUNCTION TO NAME THE BODY EXONS ###

### takes the list input from VSR_select_reference_transcript_variant()
### we write this as a separate modular function because multiple LIVs can be called, and we have to reconcile the selected variant names.
### so that means, we'll be doing a first-pass naming scheme for each individual LIV, then determining which one we want to roll with, then we refresh each LIV with the final choice.

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
        ## collapse by spaces between exons of the same LIV
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
## Behaviour: Select lowest variant for all LIVs -> rename using the common lowest variant -> for the remaining un-named intronic exons, rename if they overlapped with a separate variant.
## inputs: VSR coords (pre-checked) and a list of start/end info for each LIV. 
VSR_LIV_organise_exon_rematching <- function(VSR_coordinates, list_tibble_exon_start_end_per_LIV, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1) {
    
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
    
    # list_tibble_exon_start_end_per_LIV <- list(
    #     "LIV_1" = tribble(
    #         ~start, ~end,
    #         89740100, 89740803
    #     ),
    #     "LIV_2" = tribble(
    #         ~start, ~end,
    #         89739554, 89739993,
    #         89738975, 89739477,
    #         89738709, 89738881
    #     ),
    #     "LIV_3" = tribble(
    #         ~start, ~end,
    #         89720400, 89720582,
    #         89737976, 89738589,
    #         89739267, 89739993
    #     ),
    #     "LIV_4" = tribble(
    #         ~start, ~end,
    #         89739290, 89739477,
    #         89740100, 89740803
    #     )
    # ) 
    
    # VSR_coordinates <- input_alternative_event_region
    # list_tibble_exon_start_end_per_LIV <- list_of_exon_start_end_tibbles
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
    
    # DETECT VSR OR LIV
    ## if the list has more than one element, then it's a VSR.
    flag_is_LIV <- list_tibble_exon_start_end_per_LIV %>% length == 1
    
    # LOWEST VARIANT DETECTION FOR EACH LIV
    list_initial_lowest_variant_detection <- list_tibble_exon_start_end_per_LIV %>% 
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
    ## collapse by spaces between exons of the same LIV
    list_final_LIV_nomenclature <- purrr::map(
        .x = list_second_pass_naming,
        .f = function(a1) {
            
            # DEBUG ###
            # a1 <- list_second_pass_naming[[1]]
            ###########
            
            collapsed_nomenclature_for_one_LIV <- a1$list_body_exon_nomenclature %>% unlist %>% na.omit %>% paste(collapse = " ")
            
            return(collapsed_nomenclature_for_one_LIV)
            
        } )
    
    # finally, extract the transcript version
    variant_slot <- paste(lowest_common_variant_ID, ".", tibble_gtf_table[tibble_gtf_table$hgnc_stable_variant_ID == lowest_common_variant_ID, "transcript_version"] %>% unlist %>% na.omit %>% unique %>% .[1], sep = "")
    
    ## collapse by "/" between different LIVs and sandwich between the VSR upstream and downstream slots
    final_VSR_nomenclature <- paste(variant_slot, " ", list_first_pass_naming[[1]]$upstream_VSR_slot, " (", list_final_LIV_nomenclature %>% paste(collapse = " / "), ") ", list_first_pass_naming[[1]]$downstream_VSR_slot, sep = "") %>% 
        trimws %>%
        gsub(pattern = "\\(\\((.*)\\)\\)", replacement = "(\\1)") %>% 
        gsub(pattern = "\\(\\)", replacement = "") 
    
    return(final_VSR_nomenclature)
    
}

# END VSR_LIV_organise_exon_rematching() ###

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
triage_input_coordinates <- function(vector_input_coordinates, tibble_gtf_table, expect_stranded = NULL) {
    
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
            stop("Please check the format of your co-ordinates and make sure you have specified the strand.")
        }
        
        vector_query_chr <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
        vector_query_VSR_start <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2")
        vector_query_VSR_end <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3")
        vector_query_strand <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
        
    } else if (expect_stranded == FALSE) {
        
        # check for basic format
        if (grep(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)") %>% length != vector_input_coordinates %>% length) {
            stop("Please check the format of your co-ordinates and make sure you have specified the strand.")
        }
        
        vector_query_chr <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)", replacement = "\\1")
        vector_query_VSR_start <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)", replacement = "\\2")
        vector_query_VSR_end <- gsub(x = vector_input_coordinates, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)", replacement = "\\3")
        
    }
    
    # check chromosomes
    if (vector_query_chr %in% (tibble_gtf_table$seqnames %>% unique) %>% all == FALSE) {
        stop("Please check the format of your chromosomes")
    }
    
    # check start/end for numeric
    if (vector_query_VSR_start %>% type.convert %>% is.numeric == FALSE | vector_query_VSR_end %>% type.convert %>% is.numeric == FALSE) {
        stop("Please check the format of your start/end co-ordinates")
    }
    
    # check start/end for order. end must be > start
    if (any((vector_query_VSR_start %>% type.convert) > (vector_query_VSR_end %>% type.convert))) {
        stop("Please ensure that all the start co-ordinates are not greater than your end co-ordinates")
    }
    
}


### SHINY ####

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    fluidRow(
        
        column(3,
               
               # Application title
               titlePanel("EDN workshop"),
               
               selectInput("range_type", 
                           label = "Select the type of range to describe", 
                           choices = list("Junction", 
                                          "Exon"), 
               ),
               
               # selectInput("genome_assembly", 
               #             label = "Select the genome assembly to use (Ensembl release)", 
               #             choices = list("NCBI36 (hg18)" = c(48:54) %>% as.character, 
               #                            "GRCh37 (hg19)" = c(55:75) %>% as.character, 
               #                            "GRCh38 (hg38)" = c(76:(list.files("data") %>% gsub(pattern = "annotated_ensembl_gtf_release_(.*).txt", replacement = "\\1") %>% type.convert %>% max)) %>% as.character)),
               
               selectInput("genome_assembly", 
                           label = "Select the genome assembly to use (Ensembl release)", 
                           choices = list("NCBI36 (hg18)" = "54",
                                          "GRCh37 (hg19)" = "75",
                                          "GRCh38 (hg38)" = list.files("data") %>% gsub(pattern = "annotated_ensembl_gtf_release_(.*).txt", replacement = "\\1") %>% type.convert %>% max %>% as.character)
               ),
               
               actionButton("import_GTF", "Import genome assembly", icon = icon("file-import"),
                            class = "btn btn-primary", width = "200px"),
               
               # magnetisation options required
               radioButtons("set_tolerances", label = "Tweak match tolerances",
                            choices = list("no" = "no", "yes" = "yes"), 
                            selected = 1),
               conditionalPanel(
                   condition = "input.set_tolerances == 'yes'",
                   textInput("left_query_end_shift", label = "Left (5') query end shift (default = 0)", value = 0, placeholder = "0")
               ),
               conditionalPanel(
                   condition = "input.set_tolerances == 'yes'",
                   textInput("right_query_end_shift", label = "Right (3') query end shift (default = 0)", value = 0, placeholder = "0")
               ),
               conditionalPanel(
                   condition = "input.set_tolerances == 'yes'",
                   textInput("left_match_tolerance", label = "Left (5') match tolerance (default = 0)", value = 1, placeholder = "1")
               ),
               conditionalPanel(
                   condition = "input.set_tolerances == 'yes'",
                   textInput("right_match_tolerance", label = "Right (3') match tolerance (default = 0)", value = 1, placeholder = "1")
               ),
               
               # STRATEGY: 
               # FLI & LIV: vertex matching
               # VSRs and AEs: VSR and exon coord matching, splicemode autodetect
               # LSVs and AJs: Junction matching. Forward-slashes everywhere.
               
               textInput("input_range", label = "Enter genome-relative co-ordinates", placeholder = "e.g. 16:2756607-2757471:+"),
               
               sliderInput("slider_plot_width", "Plot width:",
                           min = 100, max = 4000, step = 100,
                           value = 800),
               sliderInput("slider_plot_height", "Plot height:",
                           min = 100, max = 4000, step = 100,
                           value = 1500),
               
               sliderInput("slider_plot_x_scale", "Base zoom x-axis:",
                           min = -5, max = 5, step = 0.01,
                           value = 1),
               sliderInput("slider_plot_y_scale", "Base zoom y-axis:",
                           min = -5, max = 5, step = 0.01,
                           value = 1.3),
               
               sliderInput("slider_plot_x_offset", "Base offset left/right:",
                           min = -5, max = 5, step = 0.01,
                           value = 0),
               sliderInput("slider_plot_y_offset", "Base offset up/down:",
                           min = -5, max = 5, step = 0.01,
                           value = 0),
               
               mainPanel(
                   uiOutput('reactive_UI_1')
                   # ,
                   # uiOutput('reactive_UI_2')
               ),
               
               shinyWidgets::materialSwitch(inputId = "plot_is_active", label = "Activate plot", status = "success", value = FALSE, right = FALSE, inline = TRUE),
              
               actionButton("reset_user_table", "Reset history",
                            class = "btn btn-primary", width = "200px"),
                
               actionButton("add_user_range", "Add range", icon = icon("searchengin"),
                            class = "btn btn-primary", width = "200px"),
               
               br(),
               br(),
               
               verbatimTextOutput("nomenclature_output", placeholder = TRUE),
               
               br(),
               
               # generate an upload box for users to upload full-length reconstructed isoforms
               downloadButton("download_plot", "Download Results"),
               
        ),
        
        column(9,
               
               h5("Click and drag + double-click to zoom. Double-click to reset."),
               
               plotOutput("plot1", height = 300,
                          dblclick = "plot1_dblclick",
                          brush = brushOpts(
                              id = "plot1_brush",
                              resetOnNew = TRUE
                          )),
               
               # svgPanZoomOutput(outputId = "plot1", width = "100%", height = "400px")
               
        )
        
    )
    
)

# Define server logic required to generate the nomenclature
server <- function(input, output) {
    
    # REACTIVE UI ####
    
    # if the user wants to manualy enter exon co-ordinates for full-length isoform or LIVs, then open up multiple text boxes for entering individual co-ords.
    output$reactive_UI_1 <- renderUI( {
        
        if (input$range_type == "Manually enter exon co-ordinates") {
            
            purrr::map(.x = 1:input$full_length_number_of_exons, .f = ~textInput(paste("full_length_exon_genome_relative_coordinate_", .x, sep = ""), label = paste("Enter the genome-relative co-ordinates of exon #", .x, sep = ""), placeholder = "e.g. 16:2756334-2756606"))
            
        } else if (input$range_type == "Local Isoform Variant (LIV)") {
            
            purrr::map(.x = 1:input$LIV_number_of_exons, .f = ~textInput(paste("LIV_exon_genome_relative_coordinate_", .x, sep = ""), label = paste("Enter the genome-relative co-ordinates of alternative exon #", .x, sep = ""), placeholder = "e.g. 16:2756334-2756606"))
            
        } else if (input$range_type == "Variable Splice Region (VSR)") {
            
            # VSRs
            # after the user has entered the number of independent events, for EACH independent event in the region, we have to ask the user how many exons there are, since alternative regions are a collection of LIVs. we're basically asking for multiple LIVs.
            purrr::map(.x = 1:input$alternative_region_number_of_independent_events, .f = ~textInput(paste("VSR_number_of_exons_for_LIV_", .x, sep = ""), label = paste("Enter the number of exons for LIV #", .x, sep = ""), value = 1))
            
        } else if (input$range_type == "Local Splice Variation (LSV) (junctions only)") {
            
            # VSRs
            # after the user has entered the number of independent events, for EACH independent event in the region, we have to ask the user how many exons there are, since alternative regions are a collection of LIVs. we're basically asking for multiple LIVs.
            purrr::map(.x = 1:input$alternative_region_number_of_independent_events, .f = ~textInput(paste("LSV_junction_genome_relative_coordinate_", .x, sep = ""), label = paste("Enter the genome-relative co-ordinates of constitutent junction #", .x, sep = ""), placeholder = "e.g. 16:2756607-2757471"))
            
        } # else if
        
    } ) # renderUI
    
    # output$reactive_UI_2 <- renderUI( {
    #     
    #     
    #     if (input$structure_type == "Variable Splice Region (VSR)") {
    #         
    #         reactive_VSR_number_of_exons_for_each_LIV <- reactive({
    #             
    #             purrr::map(
    #                 .x = 1:input$alternative_region_number_of_independent_events,
    #                 .f = ~input[[paste("VSR_number_of_exons_for_LIV_", .x, sep = "")]])
    #             
    #         })
    #         
    #         # finishing off the VSR options...
    #         # after use has entered the number of exons for each LIV,
    #         # then map across each LIV, creating textbox options for co-ordinates of exons
    #         
    #         purrr::imap(
    #             .x = reactive_VSR_number_of_exons_for_each_LIV(),
    #             .f = function(a1, a2) {
    #                 
    #                 number_of_exons <- a1 %>% paste %>% type.convert
    #                 
    #                 purrr::imap(
    #                     .x = 1:number_of_exons,
    #                     .f = function(b1, b2) {
    #                         
    #                         textInput(paste("VSR_exon_genome_relative_coordinate_exon_number_", b2, "_LIV_number_", a2, sep = ""),
    #                                   label = paste("Enter the genome-relative co-ordinates of alternative exon #", b2, " in LIV #", a2, sep = ""),
    #                                   placeholder = "16:2756334-2756606")
    #                         
    #                         # print(b1)
    #                         # print(b2)
    #                         
    #                         # print(a1)
    #                         # print(a2)
    #                         
    #                     } )
    #                 
    #             } )
    #         
    #     }
    #     
    # } )
    
    observeEvent(input$import_GTF, {
        
        showModal(modalDialog(paste("Importing release ", input$genome_assembly, ". Please wait...\n", sep = ""), footer = NULL))
        
        observeEvent(reactive_tibble_ref_gtf(), {
            
            output$nomenclature_output <- renderPrint( { cat("Importing done.\n") })
            
            removeModal()
            
        } )
        
    } )
    
    reactive_tibble_ref_gtf <- eventReactive(input$import_GTF, {
        
        # tibble_ref_gtf_hg18 <- vroom::vroom(file = "data/annotated_ensembl_gtf_release_54.txt", delim = "\t") %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
        
        # import_tibble_ref_gtf <- read.delim(file = paste("data/annotated_ensembl_gtf_release_", input$genome_assembly, ".txt", sep = ""), sep = "\t", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
        
        import_tibble_ref_gtf <- data.table::fread(file = paste("data/annotated_ensembl_gtf_release_", input$genome_assembly, ".txt", sep = ""), sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
        
        # output$nomenclature_output <- renderPrint( { cat("Importing done.\n") })
        # 
        # removeModal()
        
        return(import_tibble_ref_gtf)
        
    } )
    
    # handle plot
    
    # reactive_active_plot <- reactive({input$plot_is_active})
    
    reactive_plot_width <- reactive({input$slider_plot_width})
    reactive_plot_height <- reactive({input$slider_plot_height})
    
    reactive_plot_x_scale <- reactive({input$slider_plot_x_scale})
    reactive_plot_y_scale <- reactive({input$slider_plot_y_scale})
    
    reactive_plot_x_offset <- reactive({input$slider_plot_x_offset})
    reactive_plot_y_offset <- reactive({input$slider_plot_y_offset})
    
    observeEvent(input$reset_sliders,{
        updateSliderInput('Var1',value = 0)
        updateSliderInput('Var2',value = 0)
        updateSliderInput('Var3',value = 0)
    })
    
    # initialise user table
    reactive_tibble_user_ranges <- reactive({
        
        tibble(
            "id" = integer(0),
            "chr" = character(0),
            "start" = integer(0),
            "end" = integer(0),
            "strand" = character(0),
            "range_type" = character(0)
        ) %>% return
        
    } ) 
    
    reactive_final_plot <- reactive( {
        
        # TEST ###
        # vector_input_range <- "3:11807833-11809682:*"
        # input_chr <- "3"
        # input_start <- 11807833
        # input_end <- 11809682
        # input_strand <- "*"
        
        tibble_ref_gtf <- reactive_tibble_ref_gtf()
        
        triage_input_coordinates(vector_input_coordinates = input$input_range, tibble_gtf_table = tibble_ref_gtf, expect_stranded = TRUE)
        
        input_chr <- gsub(x = input$input_range, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\1")
        input_start <- gsub(x = input$input_range, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\2") %>% type.convert
        input_end <- gsub(x = input$input_range, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\3") %>% type.convert
        input_strand <- gsub(x = input$input_range, pattern = "^([^\\:]+)\\:([^\\-]+)\\-([^\\:]+)\\:(.*)", replacement = "\\4")
        
        plot_x_scale <- reactive_plot_x_scale() %>% as.numeric
        plot_y_scale <- reactive_plot_y_scale() %>% as.numeric
        
        plot_x_offset <- reactive_plot_x_offset() %>% as.numeric
        plot_y_offset <- reactive_plot_y_offset() %>% as.numeric
        
        # plot_view_initial_x_start <- input_start - plot_x_scale*(input_end - input_start)
        # plot_view_initial_x_end <- input_end + plot_x_scale*(input_end - input_start)
        # 
        # tibble_captured_in_range <- tibble_ref_gtf[which(tibble_ref_gtf$seqnames == input_chr & tibble_ref_gtf$start <= plot_view_initial_x_end & tibble_ref_gtf$end >= plot_view_initial_x_start & tibble_ref_gtf$strand %in% input_strand), ]
        # 
        # tibble_captured_in_range$transcript_id <- factor(tibble_captured_in_range$transcript_id, levels = tibble_captured_in_range$transcript_id %>% unique %>% na.omit %>% mixedsort(decreasing = TRUE))
        # 
        # tibble_captured_in_range$panel <- "transcripts"
        
        tibble_user_ranges <- reactive_tibble_user_ranges()
        
        tibble_user_ranges <- tibble_user_ranges %>%
            dplyr::add_row("id" = nrow(.) + 1,
                           "chr" = input_chr,
                           "start" = input_start,
                           "end" = input_end,
                           "strand" = input_strand,
                           "range_type" = input$range_type)
        
        tibble_user_ranges$panel <- "user_ranges"
        
        # create tibble that draws vertical lines matching user range to reference features ###
        tibble_selected_range <- tibble_user_ranges[tibble_user_ranges$id == nrow(tibble_user_ranges), ]
        
        selected_chr <- tibble_selected_range$chr
        selected_start <- tibble_selected_range$start
        selected_end <- tibble_selected_range$end
        selected_strand <- tibble_selected_range$strand
        
        if (selected_strand == "*") {
            selected_strand <- c("+", "-")
        }
        
        # snap range to the selected item
        plot_view_initial_x_start <- selected_start - 1.5*(selected_end - selected_start) + plot_x_scale*(selected_end - selected_start)
        plot_view_initial_x_end <- selected_end + 1.5*(selected_end - selected_start) - plot_x_scale*(selected_end - selected_start)
        
        # apply offset
        plot_view_initial_x_start <- plot_view_initial_x_start + (plot_view_initial_x_end - plot_view_initial_x_start)*plot_x_offset
        plot_view_initial_x_end <- plot_view_initial_x_end + (plot_view_initial_x_end - plot_view_initial_x_start)*plot_x_offset
        
        # plot shenanigans
        tibble_captured_in_range <- tibble_ref_gtf[which(tibble_ref_gtf$seqnames == input_chr & tibble_ref_gtf$start <= plot_view_initial_x_end & tibble_ref_gtf$end >= plot_view_initial_x_start & tibble_ref_gtf$strand %in% selected_strand), ]
        
        tibble_captured_in_range$transcript_id <- factor(tibble_captured_in_range$transcript_id, levels = tibble_captured_in_range$transcript_id %>% unique %>% na.omit %>% mixedsort(decreasing = TRUE))
        
        tibble_captured_in_range$panel <- "transcripts"
        
        tibble_ref_transcripts_overlapped_by_user_query <- extract_overlapping_features(query_chr = selected_chr, query_start = selected_start, query_end = selected_end, query_strand = tibble_selected_range$strand, tibble_gtf_table = tibble_ref_gtf, left_query_shift = 0, right_query_shift = 0, left_tolerance = 1, right_tolerance = 1, return_type = "transcript")
        
        tibble_all_exons_of_overlapped_parent_transcript <- tibble_ref_gtf[which(tibble_ref_gtf$type == "exon" & tibble_ref_gtf$transcript_id %in% (tibble_ref_transcripts_overlapped_by_user_query$transcript_id %>% unique)), ]
        
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
                left_ref_vertex_grown_from_user_query_start <- vector_all_ref_vertices[vector_all_ref_vertices < selected_start] %>% max
                right_ref_vertex_grown_from_user_query_start <- vector_all_ref_vertices[vector_all_ref_vertices > selected_start] %>% min
                left_ref_vertex_grown_from_user_query_end <- vector_all_ref_vertices[vector_all_ref_vertices < selected_end] %>% max
                right_ref_vertex_grown_from_user_query_end <- vector_all_ref_vertices[vector_all_ref_vertices > selected_end] %>% min
                
                tibble_vertices_with_distances <- tibble(
                    "transcript_id" = a1$transcript_id %>% unique,
                    "ref_vertex" = c(left_ref_vertex_grown_from_user_query_start, right_ref_vertex_grown_from_user_query_start, left_ref_vertex_grown_from_user_query_end, right_ref_vertex_grown_from_user_query_end),
                    "query_vertex" = c(selected_start, selected_start, selected_end, selected_end)
                ) %>% dplyr::mutate("ref_vertex_minus_query_vertex" = `ref_vertex` - `query_vertex`)
                
                # test for redundant overlapping distances. this happens when 1. both query ends find a ref transcript and 2. distance to left overlaps and/or distance to right overlaps.
                if (left_ref_vertex_grown_from_user_query_end < selected_start) {
                    tibble_vertices_with_distances <- tibble_vertices_with_distances[!(tibble_vertices_with_distances$ref_vertex == left_ref_vertex_grown_from_user_query_end & tibble_vertices_with_distances$query_vertex == selected_end), ]
                }
                
                if (right_ref_vertex_grown_from_user_query_start > selected_end) {
                    tibble_vertices_with_distances <- tibble_vertices_with_distances[!(tibble_vertices_with_distances$ref_vertex == right_ref_vertex_grown_from_user_query_start & tibble_vertices_with_distances$query_vertex == selected_start), ]
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
        
        tibble_distance_annotations_based_on_user_query$panel <- "transcripts"
        
        # make the distance table common with the reference plotting table
        tibble_distance_annotations_based_on_user_query <- tibble_distance_annotations_based_on_user_query[tibble_distance_annotations_based_on_user_query$transcript_id %in% tibble_captured_in_range$transcript_id, ]
        
        # invoke tolerances
        # tibble_selected_range_measurement_against_reference[abs(tibble_selected_range_measurement_against_reference$distance_to_ref) <= 1, "distance_to_ref"] <- 0
        
        ##########
        
        number_of_transcripts_captured <- tibble_captured_in_range$transcript_id %>% unique %>% length
        
        plot_view_initial_y_start <- 1 - 1.5*(number_of_transcripts_captured - 1) + plot_y_scale*(number_of_transcripts_captured - 1)
        plot_view_initial_y_end <- number_of_transcripts_captured + 1.5*(number_of_transcripts_captured - 1) - plot_y_scale*(number_of_transcripts_captured - 1)
        
        # apply offset
        ## minuses because the y axis is in reversed order.
        plot_view_initial_y_start <- plot_view_initial_y_start - (plot_view_initial_y_end - plot_view_initial_y_start)*plot_y_offset
        plot_view_initial_y_end <- plot_view_initial_y_end - (plot_view_initial_y_end - plot_view_initial_y_start)*plot_y_offset
        
        # # # Zoomable plot
        # plot_brush_ranges <- reactiveValues(
        #     x = c(plot_view_initial_x_start,
        #           plot_view_initial_x_end),
        #     y = c(plot_view_initial_y_start,
        #           plot_view_initial_y_end)
        # )
        
        # plot_brush_ranges <- reactiveValues(x = NULL, y = NULL)
        
        tibble_combined <- dplyr::full_join(tibble_captured_in_range, tibble_distance_annotations_based_on_user_query)
        
        final_plot <- ggplot() +
            facet_grid(panel ~ ., scales = "free_y") +
            
            # transcripts
            geom_segment(data = tibble_combined %>% dplyr::filter(type == "transcript"), colour = "slateblue1", arrow = arrow(angle = 30), mapping = aes(x = start, xend = end, y = transcript_id, yend = transcript_id)) +
            geom_text(data = tibble_combined %>% dplyr::filter(type == "transcript"), nudge_y = 0.25, mapping = aes(x = mean(c(plot_view_initial_x_start, plot_view_initial_x_end)), y = transcript_id, label = purrr::map2(.x = strand, .y = hgnc_stable_variant_ID, .f = function(.x, .y) {if (.x == "+") {paste("> > > > > > ", .y, " > > > > > >", sep = "")} else if (.x == "-") {paste("< < < < < < ", .y, " < < < < < <", sep = "")} else {.y} } ) %>% unlist)) +
            geom_segment(data = tibble_combined %>% dplyr::filter(type == "exon"), colour = "slateblue1", mapping = aes(x = start, xend = end, y = transcript_id, yend = transcript_id), size = 10) +
            geom_label(data = tibble_combined %>% dplyr::filter(type == "exon"), colour = "black", nudge_y = 0.15, fontface = "bold.italic", mapping = aes(x = purrr::map2(.x = start, .y = end, .f = ~c(.x, .y) %>% mean) %>% unlist, y = transcript_id, label = paste("E", exon_number, sep = ""))) +
            
            # user ranges
            geom_curve(data = tibble_user_ranges[tibble_user_ranges$range_type == "Junction", ], colour = "grey50", size = 2, curvature = -0.25, mapping = aes(x = start, xend = end, y = id, yend = id)) +
            geom_segment(data = tibble_user_ranges[tibble_user_ranges$range_type == "Exon", ], colour = "grey50", size = 5, mapping = aes(x = start, xend = end, y = id, yend = id)) +
            geom_text(data = tibble_user_ranges, colour = "black", nudge_y = 0.25, fontface = "bold", mapping = aes(x = purrr::map2(.x = start, .y = end, .f = ~c(.x, .y) %>% mean) %>% unlist, y = id, label = id)) +
            geom_vline(colour = "red", lty = 2, xintercept = selected_start) +
            geom_vline(colour = "red", lty = 2, xintercept = selected_end) +
            geom_segment(data = tibble_combined, colour = "red", arrow = arrow(angle = 45), mapping = aes(x = ref_vertex, xend = query_vertex, y = transcript_id, yend = transcript_id)) +
            geom_label(data = tibble_combined, colour = "red", nudge_y = -0.25, mapping = aes(x = purrr::map2(.x = ref_vertex, .y = query_vertex, .f = ~c(.x, .y) %>% mean) %>% unlist, y = transcript_id, label = ref_vertex_minus_query_vertex)) +
            
            ggh4x::force_panelsizes(rows = c(1, 0.3)) +
            theme_bw() +
            theme(text = element_text(family = "Helvetica"))
        
        # svgPanZoom(p, controlIconsEnabled = T, width = 800, height = 1500)
        
        return(list(
            "final_plot" = final_plot,
            "plot_view_initial_x_start" = plot_view_initial_x_start,
            "plot_view_initial_x_end" = plot_view_initial_x_end,
            "plot_view_initial_y_start" = plot_view_initial_y_start,
            "plot_view_initial_y_end" = plot_view_initial_y_end
        ) )
        
    } )
    
    # Zoomable plot
    observe( {
        
        if (input$plot_is_active == TRUE) {
            
            # smuggle in variables
            plot_height <- reactive_plot_height()
            plot_width <- reactive_plot_width()
            
            plot_brush_ranges <- reactiveValues(
                x = c(reactive_final_plot() %>% .$plot_view_initial_x_start, 
                      reactive_final_plot() %>% .$plot_view_initial_x_end), 
                y = c(reactive_final_plot() %>% .$plot_view_initial_y_start, 
                      reactive_final_plot() %>% .$plot_view_initial_y_end)
            )
            
            final_plot <- reactive_final_plot() %>% .$final_plot
            
            output$plot1 <- renderPlot( {
                
                final_plot +
                    coord_cartesian(xlim = plot_brush_ranges$x, ylim = plot_brush_ranges$y)
                
            }, height = plot_height, width = plot_width )
            
            # When a double-click happens, check if there's a brush on the plot.
            # If so, zoom to the brush bounds; if not, reset the zoom.
            observeEvent(input$plot1_dblclick, {
                
                brush <- input$plot1_brush
                if (!is.null(brush)) {
                    plot_brush_ranges$x <- c(brush$xmin, brush$xmax)
                    plot_brush_ranges$y <- c(brush$ymin, brush$ymax)          
                } else {
                    plot_brush_ranges$x <- c(reactive_final_plot() %>% .$plot_view_initial_x_start, 
                                             reactive_final_plot() %>% .$plot_view_initial_x_end)
                    plot_brush_ranges$y <- c(reactive_final_plot() %>% .$plot_view_initial_y_start, 
                                             reactive_final_plot() %>% .$plot_view_initial_y_end)
                }
                
            } )
            
        } else {
            
        } 
        
    } )
    
    
    
    
    removeModal()
    
    # Create the button to download the scatterplot as PDF
    # output$download_plot <- downloadHandler(
    #     filename = function() {
    #         paste('scatterPlot_', Sys.Date(), '.pdf', sep='')
    #     },
    #     content = function(file) {
    #         ggsave(file, makeScatPlot(), width = 11, height = 4, dpi = 300, units = "in")
    #     }
    # )
    
    
    # }, ignoreNULL = FALSE, ignoreInit = TRUE)
    
    
    
    outputOptions(output, "reactive_UI_1", suspendWhenHidden = FALSE)
    # outputOptions(output, "reactive_UI_2", suspendWhenHidden = FALSE)
    # outputOptions(output, "plot1", suspendWhenHidden = FALSE)
    # outputOptions(output, "nomenclature_output", suspendWhenHidden = FALSE)
    
} # server

# Run the application 
shinyApp(ui = ui, server = server)
