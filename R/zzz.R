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
        lazy = TRUE)
    #environment(._EA_snv_tri) <- asNamespace('aenmd')
}


