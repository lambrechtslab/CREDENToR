This repo is moved from Jieyi-DiLaKULeuven/CREDENToR.
# CREDENToR 
_(**CR**yptic **E**lements’ **D**ifferential **E**xpression by de **N**ovo **T**ranscript**o**me **R**econstruction)_

**A bioinformatics pipeline to detect and quantify cryptic transcripts associated with retrotranspon repeats.**

This pipeline is coded in Julia language. Users need to map the RNA-seq data first (hg38 or mm10, usage of STAR is recommended). REDENToR subsequently assembles transcriptomes using StringTie, under guidance of transcript annotation (Ensembl is recommended). Then, all de novo assembled transcription annotations from the given samples are merged using “StringTie --merge”. HTSeq-counts are used to count the number of reads in known and novel genes. Non-coding transcripts (transcripts not overlapping annotated coding genes) in the merged transcription annotations are assigned as cryptic transcripts when any of their exons overlaps with a retrotransposon repeat annotation (LTR, LINE, or SINE, based on RepeatMasker annotation from UCSC). If a transcript overlaps with >1 annotated repeat, the retrotransposon with the highest overlap is assigned to this cryptic transcripts.

## Dependencies
### Software
The softwares below are required:
- [Julia](https://julialang.org/) v1.2.0+
- [StringTie](https://ccb.jhu.edu/software/stringtie/)
- [HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/)

The Julia package below is required:
- [SortingAlgorithms.jl](https://github.com/JuliaCollections/SortingAlgorithms.jl)

### Data
- Mapped bam file(s), sorted by coordinates;
- ENSEMBL GTF file.

## Usage
`julia `_\[_`-p `_n \]_` credentor.jl `_\<parameters\>_

### Input
- `--bam` (or `-b`) _\<file1.bam,file2.bam,...\>_ Input mapped bam file (sorted by coordinate). Multiple input files should be seperated with a comma.
- `--id`  (or `-i`) _\<id1,id2,...\>_ ID of each sample. Should be the same length as bam files.
- `--output` (or `-o`) _\<dir\>_ The directory for the output. Should be pre-created.
- `--genome` (or `-g`) `hg38`|`mm10` The sequenced genome type. Only hg38 and mm10 are supported so far.
- `--annotation` (or `-a`) _\<file.gtf\>_ The transcription annotation file. Tested on Ensembl annotation only.
- `--seqstrand` (or `-s`) `first`|`second`|`none`  The sequenced strand. Optional (default="none").
  - `first` : Assumes a stranded library fr-firststrand.
  - `second`: Assumes a stranded library fr-secondstrand.
- `--help` (or `-h`) Print help document.
- `--version` (or `-v`) Print version.
- `-p` _\<thread_num\>_ Assign thread number. A Julia parameter that should be put between `julia` and `redentor.jl`. Optional (default=1).

### Output
`RAT_info.tsv`, `RAT_readcount.tsv`, and `RAT_RPKM.tsv` are the information, read count, and read per kilobase of transcript, per million mapped reads (RPKM) file of all detected cryptic transcripts, respectively. For the exon structure of cryptic transcripts, please see `work/StringTieMerged.gtf`.

### Tips
- When the pipeline is interrupted and rerun in the same folder, the previously finished steps will be skipped.

### Contacts
Jieyi Xiong (jieyi.xiong[at]kuleuven.be); Bernard Thienpont (bernard.thienpont[at]kuleuven.be); Diether Lambrechts (diether.lambrechts[at]kuleuven.be).
