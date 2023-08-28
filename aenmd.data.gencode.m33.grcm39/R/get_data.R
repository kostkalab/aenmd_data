#library(BSgenome.Mmusculus.UCSC.mm39)  

FORCE_CREATE = FALSE #- set to TRUE to overwrite existing data

PREFIX <- paste0(rprojroot::find_root("DESCRIPTION", path="."),"/inst/extdata/")

datfiles <- c(  "env_gencode_v43_exns_byTx.rds",
                "env_gencode_v43_set.rds",
                "env_gencode_v43_seqs_byTx.rds",
                "gr_gencode_v43_spl-mask.rds",
                "gr_gencode_v43_txs-mask.rds",
                "gr_gencode_v43_txs.rds")

if((file.exists(paste(PREFIX, datfiles, sep="/")) |> all()) && (!FORCE_CREATE)){

    #- nothing to do

} else { #- generate data.

#========================================
#- need ENSEMBL for finding canonical txs : https://annotationhub.bioconductor.org/package2/AHEnsDbs
#========================================
AH_VERSION    <- 'AH104895' #- version 107
ah            <- AnnotationHub::AnnotationHub()
edb           <- ah[[AH_VERSION]]
CANONICAL_TXS <- ensembldb::genes(edb)$canonical_transcript |> sort() |> unique()

#- ALL other info is based on the GENCODE GTF
#=============================================
#- GENCODE 33 Mice on GRCm39
gtff <- 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf.gz'
#- this keeps all info, makeTxdbFromGFF does not
ALLgtf <- rtracklayer::import(BiocFileCache::bfcrpath(BiocFileCache::BiocFileCache(), gtff))

#- info we keep for each tx, see https://www.gencodegenes.org/pages/data_format.html
#====================================================================================
TXinfoCols <- c('source',
			"transcript_id",
		   	'transcript_support_level',
			'gene_id',
			'gene_name',
			'protein_id',
			'ccdsid')

#- make a txdb
#=============
txdb <- GenomicFeatures::makeTxDbFromGFF( file = BiocFileCache::bfcrpath(BiocFileCache::BiocFileCache(), gtff),
					  format = 'gtf',
					  dataSource = 'GENCODE_m33',
					  organism = 'Mus musculus',
					  circ_seqs = 'chrM',
					  chrominfo = BSgenome::seqinfo(BSgenome.Mmusculus.UCSC.mm39::Mmusculus)) 
#Nb of transcripts: 149488

#- change seqlevelstyle
txdb <- GenomeInfoDb::keepStandardChromosomes(txdb)
GenomeInfoDb::genome(txdb) <- ""
GenomeInfoDb::seqlevelsStyle(txdb) <- 'NCBI'
GenomeInfoDb::genome(txdb) <- "GRCm39"

GenomeInfoDb::genome(ALLgtf) <- ""
GenomeInfoDb::seqlevelsStyle(ALLgtf) <- 'NCBI'
GenomeInfoDb::genome(ALLgtf) <- "GRCm39"

#- collect transcripts with coding sequence
#==========================================
txs <- GenomicFeatures::transcripts(txdb,
                               columns = c( "TXNAME","GENEID", "TXTYPE", "CDSID", "CDSNAME", "CDSSTART", "CDSEND"))

cds_tx_ind <- sapply(txs$CDSID, \(x) any(!is.na(x)))
txs_cds <- txs[cds_tx_ind]
txs_cds$GENEID |> unlist() |> unique() |> length()
#22127 -> reasonable ?

tx2gn           <- txs$GENEID |> unlist()
names(tx2gn)    <- txs$TXNAME |> unlist()
txs_exn_rng     <- GenomicFeatures::cdsBy(txdb, "tx", use.names=TRUE) #- whole transcript
txs_cds_exn_rng <- txs_exn_rng[ ( txs_cds$TXNAME |> unlist() |> unique() |> sort() ) ]

#- filter out txs with cds not divisible by 3
#--------------------------------------------
out <- which(GenomicRanges::width(txs_cds_exn_rng) |> lapply(sum) |> unlist() %% 3L != 0) |> sort()
length(out) # 10568
txs_cds_exn_rng <- txs_cds_exn_rng[-out]
length(txs_cds_exn_rng )
#[1] 56168 number of transcripts

#- only keep transcripts with at least 3 codons (start, content, end)
#--------------------------------------------------------------------
out <- which((GenomicRanges::width(txs_cds_exn_rng) |> lapply(sum) |> unlist()) <9 ) |> sort()
txs_cds_exn_rng <- txs_cds_exn_rng[-out]
length(txs_cds_exn_rng )
# [1] 56158 number of transcripts
tx2gn[txs_cds_exn_rng |> names()] |> unique() |> length()
#21686 20k genes, reasonable
CDS_EXN_RNG_BY_TX <- txs_cds_exn_rng

#- collect info for those transcripts (from ALLgtf)
#--------------------------------------------------
if(length(CDS_EXN_RNG_BY_TX |> names())!=  length(CDS_EXN_RNG_BY_TX |> unique() |> names())) stop()
txnms <- CDS_EXN_RNG_BY_TX |> names()
tmp <- ALLgtf[ALLgtf$type =='transcript']
tmp <- tmp[tmp$transcript_id %in% txnms]
names(tmp) <- tmp$transcript_id

#- remove IG genes
#-----------------
ind <- tmp$gene_type == "protein_coding"
tmp <- tmp[ind]
CDS_EXN_RNG_BY_TX <- CDS_EXN_RNG_BY_TX[names(tmp)]
TX_INFO <- tmp[,TXinfoCols]

TX_INFO$gene_id |> unique() |> length()
#- [1] 21534  about 20k protein-coding genes.

#- look up which of the transcripts are canonical txs 
#----------------------------------------------------
TX_INFO$is_canonical <- FALSE
TX_INFO$is_canonical[ (TX_INFO$transcript_id |> stringr::str_remove("\\..*")) %in% CANONICAL_TXS] <- TRUE
 
#- number of (coding & non-coding) exons
#---------------------------------------
#- need all exons (not just coding) for finding 1-exon TXs
exns_all     <- GenomicFeatures::exonsBy(txdb, "tx", use.names=TRUE) #- whole transcript
exns_all     <- exns_all[names(CDS_EXN_RNG_BY_TX)]
num_exns <- S4Vectors::elementNROWS(exns_all)
TX_INFO$num_exons <- num_exns[names(TX_INFO)] |> as.integer()

#- convert transcript support level to integer
#- still has NAs here 
#---------------------------------------------
TX_INFO$transcript_support_level <- TX_INFO$transcript_support_level |> as.integer() 

#- add genome info
GenomeInfoDb::seqinfo(TX_INFO) <- GenomeInfoDb::seqinfo(txdb)

#- SAVE RESULT
#-------------
saveRDS(CDS_EXN_RNG_BY_TX,     file = "../inst/extdata/grl_gencode_v43_exns_byTx.rds")
################################################################################
saveRDS(TX_INFO,    		   file = paste0(PREFIX,"/gr_gencode_v43_txs.rds"))
################################################################################

#- for each transcript, assemble its coding sequence
#===================================================

Mmusculus <- BSgenome.Mmusculus.UCSC.mm39::Mmusculus
si <- seqinfo(Mmusculus)
genome(si) <- ""
GenomeInfoDb::seqlevelsStyle(si) <- 'NCBI'
genome(si) <- 'GRCm39'
seqinfo(Mmusculus) <- si
cds_used     <- BSgenome::getSeq(Mmusculus, CDS_EXN_RNG_BY_TX)

exon_env   <- new.env(hash=TRUE, parent = emptyenv())
cds_env    <- new.env(hash=TRUE, parent = emptyenv())
set_env    <- new.env(hash=TRUE, parent = emptyenv()) #- list of single exon transcripts

#- Fill envs with content... slow....
for( i in seq_len(length(TX_INFO))){
        if(i %% 10 == 0) message( i )
        txn               <- names(TX_INFO)[i]
		nex               <- TX_INFO$num_exons[i]
		#- put in env containing exons
        exns              <- CDS_EXN_RNG_BY_TX[[txn]]
        exon_env[[txn]]   <- exns

		#- keep reference coding sequence for tx
        seq_ref           <- cds_used[[txn]] |> unlist() #- unlisting concatenates
        cds_env[[txn]]    <- seq_ref

		#- keep track with transcripts are "true" 1-exons (not just coding 1-exon)
		if(nex==1) set_env[[txn]] <- TRUE #- true 1-exon transcripts

}
#- exons we need fast, so we put them into an environment

###########################################################################
saveRDS(exon_env,  file = paste0(PREFIX, "/env_gencode_v43_exns_byTx.rds"))
saveRDS(cds_env,   file = paste0(PREFIX, "/env_gencode_v43_seqs_byTx.rds"))
saveRDS(set_env,   file = paste0(PREFIX, "/env_gencode_v43_set.rds"))
###########################################################################


#==============================================
#- MASK of all the coding sequence we consider
#- with map to overlapping transcripts
#==============================================

tmp               <- CDS_EXN_RNG_BY_TX |> unlist()
strand(tmp)       <- '*' #- don't need strand info
tmp$transcript_id <- names(tmp)
tst               <- reduce(tmp, with.revmap = TRUE)
tst_tx            <- lapply(tst$revmap, \(x) tmp$transcript_id[x])
tst$transcript_id <- tst_tx
tst$revmap        <- NULL #- don't need that column

CDS_TO_TX <- tst
#CDS_TO_TX 
#GRanges object with 204560 ranges and 1 metadata column:
#           seqnames        ranges strand |                       transcript_id

###########################################################################
saveRDS(CDS_TO_TX,   file = "../inst/extdata/gr_gencode_v43_txs-mask.rds")
###########################################################################


#- need a splice mask for efficient filtering
#============================================


#- exons (with utrs) ; but for ANY protein-coding transcript on standard chromosomes
ex   <- exns_all #- see above for exns_all
ex_s <- resize(ex, width=3, fix = "start")
ex_e <- resize(ex, width=3, fix = "end")
sr_s <- resize(ex_s, width = 11, fix = "end")
sr_e <- resize(ex_e, width = 11, fix = "start")
exn_spl_rng <-  pc(sr_s,sr_e) |> sort()

tmp               <- exn_spl_rng |> unlist()
strand(tmp)       <- '*'
tmp$transcript_id <- names(tmp)
tst               <- reduce(tmp, with.revmap = TRUE)
tst_tx            <- lapply(tst$revmap, \(x) tmp$transcript_id[x])
tst$transcript_id <- tst_tx
tst$revmap        <- NULL #- don't need that column

SPL_MSK <- tst

########################################################################
saveRDS(SPL_MSK,   file = paste0(PREFIX, "/gr_gencode_v43_spl-mask.rds"))
########################################################################

}