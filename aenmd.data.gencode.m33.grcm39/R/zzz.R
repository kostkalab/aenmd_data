
#- Load processed ANNOTATION INFORMATION
#---------------------------------------

#- we keep the info a dedicated environment, in the package name

._EA_exn_env      <- NULL #- exons we use
._EA_cds_env      <- NULL #- cds sequences we use
._EA_txs_gr       <- NULL #- mask for all splice regions
._EA_txs_mask_gr  <- NULL #- mask for all splice regions
._EA_spl_mask_gr  <- NULL #- mask for all splice regions
._EA_set_env      <- NULL #- dictionary of single exon transcripts 
._EA_snv_tri      <- NULL #- SNVs that cause stop-gains in transcripts we use


.onLoad <- function(libname, pkgname) {

    #- load envs containing data we need
    #-----------------------------------
    prefix       <- system.file("extdata", package = "aenmd.data.gencode.v43")

    ._EA_exn_env     <<- future::future({readRDS(paste0(prefix,'/','env_gencode_v43_exns_byTx.rds'))}, lazy = TRUE)
    ._EA_cds_env     <<- future::future({readRDS(paste0(prefix,'/','env_gencode_v43_seqs_byTx.rds'))}, lazy = TRUE)
    ._EA_txs_gr      <<- future::future({readRDS(paste0(prefix,'/','gr_gencode_v43_txs.rds'))}, lazy = TRUE)
    ._EA_txs_mask_gr <<- future::future({readRDS(paste0(prefix,'/','gr_gencode_v43_txs-mask.rds'))}, lazy = TRUE)
    ._EA_set_env     <<- future::future({readRDS(paste0(prefix,'/','env_gencode_v43_set.rds'))}, lazy = TRUE)
    ._EA_spl_mask_gr <<- future::future({readRDS(paste0(prefix,'/','gr_gencode_v43_spl-mask.rds'))}, lazy = TRUE)
    ._EA_snv_tri     <<- future::future({
        m_keys <- readRDS(paste0(prefix,'/','tri-keys_gencode_v43_all-stop-making-snvs.rds'))
        m_vals <- readRDS(paste0(prefix,'/','tri-vals_gencode_v43_all-stop-making-snvs.rds'))
        m_trie <- triebeard::trie(keys = m_keys, values = m_vals)},
        lazy = TRUE, seed = TRUE)
}

#' Exons (coding DNA) of all the transcripts in this package's transcript set
#'
#' @name EA_exn_env
#' @docType data
#' @author Dennis Kostka
#' @keywords data
#' @details Environment, with transcript names as keys and exons (as \code{GenomicRanges::GRanges}) as values.
NULL

#' DNA sequence (coding) of all the transcripts in this package's transcript set
#'
#' @name EA_cds_env
#' @docType data
#' @author Dennis Kostka
#' @keywords data
#' @details Envirnonment, with transcript names as keys, and DNA sequence (as \code{Biostrings::DNAString}) as values.
NULL

#' Splice regions of all the transcripts in this package's transcript set, reduced to non-overlapping ranges.
#'
#' @name EA_spl_mask_gr
#' @docType data
#' @author Dennis Kostka
#' @keywords data
#' @details GRanges object containing splice regions. Transcripts where splice regions are overlapping are in the metadata.
#' Splice regions are 3bp of the exon and 8bp of the intron.
NULL

#' Single exon transcripts contained in this package's transcript set
#'
#' @name EA_set_env
#' @docType data
#' @author Dennis Kostka
#' @keywords data
#' @details Envirnonment, with transcript names *of single exon transcripts only* as keys, and \code{TRUE} as values.
NULL

#' Splice regions of all the transcripts in this package's transcript set
#'
#' @name EA_spl_grl
#' @docType data
#' @author Dennis Kostka
#' @keywords data
#' @details Splice regions of all transcripts, as \code{GenomicRanges::GRangesList}.
NULL

#' Exons (coding DNA only) of all the transcripts in this package's transcript set
#'
#' @name EA_exn_grl
#' @docType data
#' @author Dennis Kostka
#' @keywords data
#' @details Exons of all transcripts, as \code{GenomicRanges::GRangesList}.
NULL

#' Transcripts in package's transcript set
#'
#' @name EA_txs_grl
#' @docType data
#' @author Dennis Kostka
#' @keywords data
#' @details All transcripts in this package's transcript list as \code{GenomicRanges::GRangesList}.
NULL

#' Keys and values containing all STOP-generating SNVs in this package's transcript  set
#'
#' @name EA_snv_tri
#' @docType data
#' @author Dennis Kostka
#' @keywords data
#' @details A \code{triebeard::trie}. Contains the keys for all SNVs that generate 
#' stop codons. The key format is \code{chr:start}, where start is a zero-padded 9-digit number.
#' The values are strings, each string contains all ensembl transcript ids for which the SNV 
#' causes a stop. Transcript ids in the string are separated by the "|" character.
NULL
