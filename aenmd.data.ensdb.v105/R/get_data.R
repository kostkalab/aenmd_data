
#- THIS FILE MAKES DATA aenmd needs
#- Currently using:
#
#  AnnotationHub / Ensembl v. 105 ("AH98047")

FORCE_CREATE = FALSE #- set to TRUE to overwrite existing data

prefix <- here::here("inst","extdata")

datfiles <- c(  "env_ensdb_v105_exns_byTx_fil.rds",
                "env_ensdb_v105_seqs_byTx_fil.rds",
                "env_ensdb_v105_splc_byTx_fil.rds",
                "env_ensdb_v105_setx_byTx_fil.rds",
                "grl_ensdb_v105_splc_byTx_fil.rds",
                "grl_ensdb_v105_exns_byTx_fil.rds",
                "grl_ensdb_v105_trnscrpts_fil.rds")

if((file.exists(paste(prefix, datfiles, sep="/")) |> all()) && (!FORCE_CREATE)){
    
    #- nothing to do

} else { #- generate data.

Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens #- BSgenome.Hsapiens.UCSC.hg38_1.4.4
seqlevelsStyle(Hsapiens) <- 'NCBI'
genome(Hsapiens) <- 'GRCh38' #- ugly, but otherwise the ensdbmerge down there does not work.
                              #- also pretty sure this is corredt in that the ensembl annoations are  are also a newer patch 

AH_VERSION <- 'AH98047'

#- Get the data
ah   <- AnnotationHub::AnnotationHub()
edb  <- ah[[AH_VERSION]]

#- COMPILE TRANSCRIPT SET
#------------------------

#- all protein-coding txs on standard chromosomes
flt_chr <- AnnotationFilter::SeqNameFilter(as.character(c(1:22,"X","Y","MT")))
flt_tbt <- AnnotationFilter::TxBiotypeFilter("protein_coding")
flt_lst <- AnnotationFilter::AnnotationFilterList(flt_chr, flt_tbt, logicOp = "&")

txs <- ensembldb::transcripts(edb, filter = flt_lst,
                               columns = c( "tx_id_version","gene_id", "tx_support_level", "tx_is_canonical"))

#- corresponding exon ranges
all_exn_rng        <- ensembldb::cdsBy(edb, "tx", filter = AnnotationFilter::TxIdFilter(names(txs)))   #- cds only
all_tx_exn_rng     <- ensembldb::exonsBy(edb, "tx", filter = AnnotationFilter::TxIdFilter(names(txs))) #- whole transcript
#- these are exactly the same transcripts
if(! all(names(all_exn_rng) == names(all_tx_exn_rng)) ) stop("Transcripts naming issue")
txs                <- txs[names(all_exn_rng),]
txs$num_exons_cds  <- S4Vectors::elementNROWS(all_exn_rng)
txs$num_exons_txs  <- S4Vectors::elementNROWS(all_tx_exn_rng)
#- sanity check
if(! all(txs$num_exons_cds <= txs$num_exons_txs) ) stop("Transcripts error counting issue")

#- tsl 1 transcripts plus
#- recover 1-exon transcripts withou transcript support level info
ind <- txs$tx_support_level == 1
ind <- ind | ((txs$num_exons_txs == 1) & (is.na(txs$tx_support_level)))

txs_used     <- txs[sort(which(ind))]
exn_rng_used <- all_exn_rng[names(txs_used)]
cds_used     <- BSgenome::getSeq(Hsapiens, exn_rng_used)

#- FILTER CDS must be divisible by three
#- some of these sequences are sort of incomplete
#  in this case the cds is not divisible by three.
#  we will exclude these but keep track.
out <- which(width(exn_rng_used) |> lapply(sum) |> unlist() %% 3 != 0) |> sort()
length(out)
# 859
write.table(names(out) |> sort(), file="../inst/extdata/enst_length_excluded.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
txs_used     <- txs_used[ -out  ]
exn_rng_used <- exn_rng_used[ -out  ]
cds_used     <- cds_used[ -out  ]

#- TODO:
#- FILTER CDS must have at least 3 codons (start, content, stop)

#- annotate splice regions based on:
#  3bp into the exon and 8bp into the intron;

#- exons (with utrs) ; but for ANY protein-coding transcript on standard chromosomes
ex   <- ensembldb::exonsBy(edb, "tx", filter = flt_lst, columns = c("exon_id", "tx_name"))
ex_s <- resize(ex, width=3, fix = "start")
ex_e <- resize(ex, width=3, fix = "end")
sr_s <- resize(ex_s, width = 11, fix = "end")
sr_e <- resize(ex_e, width = 11, fix = "start")
exn_spl_rng <-  pc(sr_s,sr_e) |> sort()

#- aggregate for all txs_used ALL overlapping splice regions
#  This is not strictly necessary, could actually just filter 
#  based on all the splice regions above...
ov <- GenomicRanges::findOverlaps(txs_used, exn_spl_rng)

afu <- function(ind){
	#tx_nme                <- names(txs_used)[ind]
    sh                    <- S4Vectors::subjectHits(ov)[which(S4Vectors::queryHits(ov)==ind)]
	res                   <- exn_spl_rng[sh] |> unlist() |> sort()
	mcols(res)$tx_biotype <- NULL
	return(res)
}

#- FIXME: slow, there needs to be a better way.
spl_rng_used        <- sapply(seq_len(length(txs_used)), afu) |> GenomicRanges::GRangesList()
names(spl_rng_used) <- names(txs_used)

#- NOTE
#======
#  It appears several orders of magnitude faster to
#  look up GRanges in an environment compared with a GRangesList.
#  Therefore we make these.
#  Also, we pre-make all the reference sequences.
#  Takes long to make, though.

#- Empty envs
exon_env   <- new.env(hash=TRUE, parent = emptyenv())
cds_env    <- new.env(hash=TRUE, parent = emptyenv())
splice_env <- new.env(hash=TRUE, parent = emptyenv())
set_env    <- new.env(hash=TRUE, parent = emptyenv()) #- list of single exon transcripts

#- Fill envs with content
for( i in seq_len(length(txs_used))){
        if(i %% 10 == 0) message( i )
        txn               <- names(txs_used)[i]
        exns              <- exn_rng_used[[txn]]
        exon_env[[txn]]   <- exns
        seq_ref_exns      <- cds_used[[txn]]
        seq_ref_c         <- unlist(seq_ref_exns)
        cds_env[[txn]]    <- seq_ref_c
    	spl               <- spl_rng_used[[txn]]
	splice_env[[txn]] <- spl
	if(txs_used$num_exons_txs[i] == 1){
		set_env[[txn]] <- TRUE #- true 1-exon transcripts
	}
}

#- we put the environments into inst/extdata
if(TRUE){
    saveRDS(exon_env,     file = "../inst/extdata/env_ensdb_v105_exns_byTx_fil.rds")
    saveRDS(cds_env,      file = "../inst/extdata/env_ensdb_v105_seqs_byTx_fil.rds")
    saveRDS(splice_env,   file = "../inst/extdata/env_ensdb_v105_splc_byTx_fil.rds")
    saveRDS(set_env,      file = "../inst/extdata/env_ensdb_v105_setx_byTx_fil.rds")
    saveRDS(spl_rng_used, file = "../inst/extdata/grl_ensdb_v105_splc_byTx_fil.rds")
    saveRDS(exn_rng_used, file = "../inst/extdata/grl_ensdb_v105_exns_byTx_fil.rds")
    saveRDS(txs_used,     file = "../inst/extdata/grl_ensdb_v105_trnscrpts_fil.rds")
}

} #- end else