## aenmd_data

This repository contains R data packages for the `aenmd` package that annotates transcripts with variants that cause premature termination codons from nonsense-mediated decay. These annotations depend data structures derived from transcript models, and here we provide them for different transcript annotations for the GRCh37 and GRCh38 assemblies.

Each annotation package provide transcript information via accessor functions that are used by `aenmd`:

**1) Information about transcripts:**

Transcripts are identified by their GENCODE/ENSEMBL transcript ID, e.g. `ENSG00000187634.13`. 

* *Transcripts in transcript set.* `get_txs()` List of transcripts, with information:
    - Its annoation source
    - Transcript ID
    - Transcript support level
    - Gene ID
    - Protein ID
    - CCDS ID
    - Whether the transcript is canonical as per ENSEMBL
* *Transcript coding exons.* `get_exns_by_tx()` Coding exons for all coding transcripts.
* *Transcript coding sequence.* `get_cds_by_tx()` Coding (reference) sequence for a transcript. 
* *Transcript mask.* `get_txs_mask()` Reduced mask for coding sequence considered in that package.
* *Single exon transcripts.* `is_single_exn_tx()` Query whether a transcript only has a single exon.
* *Splice regions.* `get_splice_mask()` Reduced mask for splice regions.

**2) Information about PTC-generating SNVs:**

SNVs are identified by their locations and alleles, e.g. `10:000047074|C|A`. 

* *All possible PTC-generating SNVs.* `is_ptc_snv()` Query whether a SNV generates a PTC in a transcript. Example:
  ```
  > is_ptc_snv(c("10:000047074|C|A","10:000047074|C|A"),'ENST00000568584.6')
  [1] TRUE TRUE
  ```




