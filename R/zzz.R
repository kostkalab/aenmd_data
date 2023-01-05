#- Load processed ENSEMBL annotation
#-----------------------------------

._EA_exn_env <- NULL #- exons we use
._EA_cds_env <- NULL #- cds sequences we use
._EA_spl_env <- NULL #- splice regions we use
._EA_spl_grl <- NULL #- splice regions we use as grl
._EA_exn_grl <- NULL #- exons we use as grl
._EA_txs_grl <- NULL #- transcripts we use as grl
._EA_snv_tri <- NULL #- SNVs that caust stop-gains in transcripts we use
._EA_set_env <- NULL #- list of single exon transcripts ; TODO: make me a tri


.onLoad <- function(libname, pkgname) {

    #- load envs containing data we need
    #-----------------------------------
    prefix       <- system.file("extdata", package = "aenmd.data.ensdb.v105")

    ._EA_exn_env <<- future::future({readRDS(paste0(prefix,'/','env_ensdb_v105_exns_byTx_fil.rds'))}, lazy = TRUE)
    #environment(._EA_exn_env) <- asNamespace('aenmd')

    ._EA_cds_env <<- future::future({readRDS(paste0(prefix,'/','env_ensdb_v105_seqs_byTx_fil.rds'))}, lazy = TRUE)
    #environment(._EA_cds_env) <- asNamespace('aenmd')

    ._EA_spl_env <<- future::future({readRDS(paste0(prefix,'/','env_ensdb_v105_splc_byTx_fil.rds'))}, lazy = TRUE)
    #environment(._EA_spl_env) <- asNamespace('aenmd')

    ._EA_set_env <<- future::future({readRDS(paste0(prefix,'/','env_ensdb_v105_setx_byTx_fil.rds'))}, lazy = TRUE)
    #environment(._EA_set_env) <- asNamespace('aenmd')

    ._EA_spl_grl <<- future::future({readRDS(paste0(prefix,'/','grl_ensdb_v105_splc_byTx_fil.rds'))}, lazy = TRUE)
    #environment(._EA_spl_grl) <- asNamespace('aenmd')

    ._EA_exn_grl <<- future::future({readRDS(paste0(prefix,'/','grl_ensdb_v105_exns_byTx_fil.rds'))}, lazy = TRUE)
    #environment(._EA_exn_grl) <- asNamespace('aenmd')

    ._EA_txs_grl <<- future::future({readRDS(paste0(prefix,'/','grl_ensdb_v105_trnscrpts_fil.rds'))}, lazy = TRUE)
    #environment(._EA_txl_grl) <- asNamespace('aenmd')

    ._EA_snv_tri <<- future::future({
        m_keys <- readRDS(paste0(prefix,'/','tri-keys_ensdb_v105_fil_all-stop-making-snvs.rds'))
        m_vals <- readRDS(paste0(prefix,'/','tri-vals_ensdb_v105_fil_all-stop-making-snvs.rds'))
        m_trie <- triebeard::trie(keys = m_keys, values = m_vals)},
        lazy = TRUE, seed = TRUE)
    #environment(._EA_snv_tri) <- asNamespace('aenmd')
}

#' Exons (coding DNA) of all the transcripts in this package's transcript set
#'
#' @name EA_exn_env
#' @docType data
#' @author Dennis Kostka
#' @keywords data
#' @details Environment, with transcript names as keys and exons (as \code{GenomicRanges::GRanges}) as values.
NULL

#' DNA sequence (codong) of all the transcripts in this package's transcript set
#'
#' @name EA_cds_env
#' @docType data
#' @author Dennis Kostka
#' @keywords data
#' @details Envirnonment, with transcript names as keys, and DNA sequence (as \code{Biostrings::DNAString}) as values.
NULL

#' Splice regions of all the transcripts in this package's transcript set
#'
#' @name EA_spl_env
#' @docType data
#' @author Dennis Kostka
#' @keywords data
#' @details Envirnonment, with transcript names as keys, and splice regions (as \code{GenomicRanges::GRangesList}) as values.
#' Splice regions extend 3bp into the exon and 8bp into the intron.
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
