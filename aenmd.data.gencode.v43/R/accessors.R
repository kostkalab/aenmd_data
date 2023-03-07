
#- Functions that provide access to annotation data. 
#- ad = annotation data
#
#- aenmd package then uses this function to access data.
#
#  At startup, or when switching annotation packages, 
#  aenmd then uses
#
#  FNAME <- get(FNAME, envir = asNamespace(._EA_dataPackage_name))
#  
#  for actually getting access to annotation data.
#  This works b/c get returns function and environment, which is the
#  specific annotation data package's namespace.
#  ._EA_dataPackage_name is a string variable aenmd maintains that
#  contains the annotation package it is currently using.


#' Genome info
#' @param details Logical. If TRUE, return seqinfo object. Default: FALSE.
#' @return GRanges, one range per transcript 
#' 
ad_get_genome <- function(details = FALSE){
#=====================
    if(details) return( GenomeInfoDb::seqinfo(ad_get_txs()))
    return("GRCh38")
}

#' Transcripts considered
#' @return GRanges, one range per transcript 
#' 
ad_get_txs <- function(){
#=====================
    return(future::value(._EA_txs_gr))
}

#' Transcript mask
#' @return GRanges, one range for exon - to - transcript mapping 
#' 
ad_get_txs_mask <- function(){
    return(future::value(._EA_txs_mask_gr))
}

#' Splice mask
#' @return GRanges splice regions
#' 
ad_get_spl_mask <- function(){
    return(future::value(._EA_spl_mask_gr))
}

#' Exons, stratified by transcript
#' @param txname String. Name of transcript 
#' @return GRanges. Exons, sorted 5' to 3' 
#' 
ad_get_exns_by_tx <- function(txname){
    return(get0(txname, future::value(._EA_exn_env)))
}

#' Coding sequence, stratified by transcript
#' @param txname String. Name of transcript 
#' @return DNAString. Coding sequence 
#' 
ad_get_cds_by_tx <- function(txname){
    return(get0(txname, future::value(._EA_cds_env)))
}

#' Query for single exon transcripts
#' @param txname String. Name of transcript 
#' @return Logical. TRUE if txname is a single exon transcript.
#' @details Caveat: This function returns "FALSE" for single exon transcripts that are not in the transcript set.
#' 
ad_is_single_exn_tx <- function(txname){
        return(exists(txname, future::value(._EA_set_env)))
}

#' Query stop-making SNVs by transcript
#' @param keys Character vector. SNVs to be checked.
#' @param txname String. Transcript name. 
#' @return Logical. For each key, whether it is PTC-generating in txname.
ad_is_ptc_snv <- function(keys, txname = NULL){
    tmp <- triebeard::longest_match(future::value(._EA_snv_tri) ,keys)
    ind <- !is.na(tmp) #- these make PTCs

    if(!is.null(txname)){
        #- these make PTCs in a specific transcript
        ind2     <- stringr::str_detect(tmp[ind], pattern = txname)
        ind[ind] <- ind2
    }

    return(ind)
}





