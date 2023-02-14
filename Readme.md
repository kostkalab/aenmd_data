## aenmd_data

This repository contains R data packages for the `aenmd` package that annotates transcripts with variants that cause premature termination codons from nonsense-mediated decay. These annotations depend data structures derived from transcript models, and here we provide them for different transcript annotations for the GRCh37 and GRCh38 assemblies. This makes it straight-forward to use `aenmd` with variant annotations from either assembly:

***Annotating variants from GRCh37 or GRCh38 assemblies*** 

```
> library(aenmd)
> vars_grch38 <- 
> vars_grch37 <- 

#- annotating variants from GRCh38
> aenmd::ad_get_genome()
> vars_grch38_annotated <- 

#- annotating variants from GRCh37
#- load GRCh37 transcripts
> library(aenmd.data.gencode.v43.grch37)
#- ask aenmd for the genome version we are using
> ad_get_genome()
#- change to GRCh37annotation
> ad_change_annotation('aenmd.data.gencode.v43.grch37')
#- ask aenmd for the genome version we are using
> ad_get_genome()
#- annotate variants using GRCh37 annotation
> vars_grch37_annotated <- 
```

***More information***

Each annotation package provide transcript information, which in turn is available via the `aenmd` package to the user. Data packages implement the functions listed below. When changing data packages, `aenmd` replaces said functions. Below, transcripts  are identified by their GENCODE/ENSEMBL transcript ID, e.g. `txname <- ENSG00000187634.13`; SNVs are identified by their locations and alleles, e.g. `snvname <- 10:000047074|C|A`. 

*Information about transcripts:*

* *Transcripts in transcript set.* `ad_get_txs()` List of transcripts in transcript set.
* *Transcript coding exons.* `ad_get_exns_by_tx(txname)` Coding exons for a ranscript.
* *Transcript coding sequence.* `ad_get_cds_by_tx(txname)` Coding (reference) sequence for a transcript. 
* *Transcript mask.* `ad_get_txs_mask()` Reduced genomic ranges for coding sequence considered in that package.
* *Single exon transcripts.* `ad_is_single_exn_tx(txname)` Query whether a transcript only has a single exon.
* *Splice regions.* `ad_get_splice_mask()` Reduced genomic ranges for for splice regions.

*Information about PTC-generating SNVs:*

* *PTC-generating SNVs.* `ad_is_ptc_snv(snvname)` Query whether a SNV generates a PTC in a transcript. Example:




