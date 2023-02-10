
#' Transcripts considered
#' @return GRanges, one range per transcript 
#' 
get_txs <- function(){
#=====================
    return(future::value(._EA_txs_gr))
}

#' Transcript mask
#' @return GRanges, one range for exon - to - transcript mapping 
#' 
get_txs_mask <- function(){
    return(future::value(._EA_txs_mask_gr))
}

#' Splice mask
#' @return GRanges splice regions
#' 
get_spl_mask <- function(){
    return(future::value(._EA_spl_mask_gr))
}

#' Exons, stratified by transcript
#' @param txname String. Name of transcript 
#' @return GRanges. Exons, sorted 5' to 3' 
#' 
get_exns_by_tx <- function(txname){
    return(get0(txname, future::value(._EA_exn_env)))
}

#' Coding sequence, stratified by transcript
#' @param txname String. Name of transcript 
#' @return DNAString. Coding sequence 
#' 
get_cds_by_tx <- function(txname){
    return(get0(txname, future::value(._EA_cds_env)))
}

#' Query for single exon transcripts
#' @param txname String. Name of transcript 
#' @return Logical. TRUE if txname is a single exon transcript.
#' @details Caveat: This function returns "FALSE" for single exon transcripts that are not in the transcript set.
#' 
is_single_exn_tx <- function(txname){
        return(exists(txname, future::value(._EA_sel_env)))
}

#' Query stop-making SNVs by transcript
#' @param keys Character vector. SNVs to be checked.
#' @param txname String. Transcript name. 
#' @return Logical. For each key, whether it is PTC-generating in txname.
is_ptc_snv <- function(keys, txname){
    triebeard::longest_match(future::value(._EA_snv_tri) ,keys) |>
                  stringr::str_detect(pattern=txname)
}






