## aenmd_data

### Data packages for the `aenmd` R package

This repository contains R data packages for the `aenmd` package; `aenmd` annotates transcripts with variants that cause premature termination codons with regard to predicted escape nonsense-mediated decay (NMD). `aenmd` uses data objects derived from transcript models, and here we provide them for transcript sets for the GRCh37 and GRCh38 assemblies. This makes it straight-forward to use `aenmd` with variant annotations from either assembly.

Currently, the repository contains three data packages:

- `aenmd.data.ensdb.v105` High-confidence transcript set based on Ensembl version 105. High-confidence means protein-coding transcripts with transcript support level 1 (or single exon transcripts).

- `aenmd.data.gencode.v43` Protein-coding transcripts from [GENCODE version 43](https://www.gencodegenes.org/human/release_43.html).

- `aenmd.data.gencode.v43.grch37` Protein-coding transcripts from [GENCODE version 43, mapped to GRCh37](https://www.gencodegenes.org/human/release_43lift37.html).

#### Installation

Install packages by selecting the appropriate subdirectory from this `aenmd_data` repository.

```R
#- install ENSEMBL v105 transcripts on GRCh38
remotes::install_github(repo = "kostkalab/aenmd_data",
                          subdir = "aenmd.data.ensdb.v105")

#- install GENCODE v43 transcripts on GRCh38
remotes::install_github(repo = "kostkalab/aenmd_data",
                          subdir = 'aenmd.data.gencode.v43')

#- install GENCODE v43 transcripts on GRCh37
remotes::install_github(repo = "kostkalab/aenmd_data",
                          subdir = 'aenmd.data.gencode.v43.grch37')

```

####Annotating variants from GRCh37 or GRCh38 assemblies

```R
library(aenmd)

#- load variants (two assemblies)
vars_grch38 <- system.file('extdata/clinvar_20221211_noinfo_sample1k.vcf.gz', package = 'aenmd') |> 
               aenmd:::parse_vcf_vcfR()

vars_grch37 <- system.file('extdata/clinvar_grch37_20221211_noinfo_sample1k.vcf.gz', package = 'aenmd') |> 
               aenmd:::parse_vcf_vcfR()

#- ANALYZE GRCH38 VARIANTS
#-------------------------

#- per default we have
ad_get_packagename()
#[1] "aenmd.data.ensdb.v105"
ad_get_genome()
#[1] "GRCh38"

vars_rng_38 <- process_variants(vars_grch38$vcf_rng)
res_38      <- annotate_nmd(vars_rng_38)

#- ANALYZE GRCH37 VARIANTS
#-------------------------

#- change to a GRCh37 annotation package
ad_swap_annotation('aenmd.data.gencode.v43.grch37')

ad_get_packagename()
#[[1] "aenmd.data.gencode.v43.grch37"
ad_get_genome()
#[1] "GRCh37"

#- process variants
vars_rng_37 <- process_variants(vars_grch37$vcf_rng)

#- filter out out variants where we don't know the alternative 100%
n_N         <- Biostrings::alphabetFrequency(vars_rng_37$alt) [,-(1:4)] |> rowSums()
ind_out     <- which(n_N > 0)
vars_rng_37 <- vars_rng_37[-ind_out]

#- finally, annotate on GRCH37
res_37  <- annotate_nmd(vars_rng_37)

```

***Additional information***

Each annotation package provide transcript information, which in turn is available via the `aenmd` package to the user. Typically the functions below are used by `aenmd` internally, but users can use them to get information about `aenmd`s transcript set, and about PTC-generating SNVs. 
Below, transcripts  are identified by their GENCODE/ENSEMBL transcript ID, e.g. `txname <- ENSG00000187634.13`; SNVs are identified by their location and alleles, e.g. `snvname <- 10:000047074|C|A`. 

*Information about transcripts:*

* *Transcripts in transcript set.* `ad_get_txs()` List of transcripts in transcript set.
* *Transcript coding exons.* `ad_get_exns_by_tx(txname)` Coding exons for a ranscript.
* *Transcript coding sequence.* `ad_get_cds_by_tx(txname)` Coding (reference) sequence for a transcript. 
* *Transcript mask.* `ad_get_txs_mask()` Reduced genomic ranges for coding sequence considered in that package.
* *Single exon transcripts.* `ad_is_single_exn_tx(txname)` Query whether a transcript only has a single exon.
* *Splice regions.* `ad_get_splice_mask()` Reduced genomic ranges for for splice regions.

*Information about PTC-generating SNVs:*

* *PTC-generating SNVs.* `ad_is_ptc_snv(snvname)` Does a SNV generate a PTC in a transcript? 




