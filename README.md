# REDENToR
A pipeline to quantify retrotransposon-associated transcriptions (RAT).

This pipeline was coded in Julia language. User need map the RNA-seq data first (hg38 or mm10, using STAR is recommended).In REDENToR, transcriptomes were subsequently assembled using StringTie, under guidance of the transcript annotation (Ensembl is recommended). All de novo assembled transcription annotations from the given samples were merged using “StringTie --merge”. HTSeq-counts were used to count the read numbers of known and novel genes. Non-coding transcripts (transcripts not overlapping annotated coding genes) in the merged transcription annotations were assigned as RATs when any of its exons overlapped with a retrotransposon repeat annotation (LTR, LINE, or SINE, based on RepeatMasker annotation from UCSC). If a transcript overlapped with >1 annotated repeat, the retrotransposon with the highest overlap was assigned to this RAT.

## Prerequirement
### Software
Below softwares should be preinstalled and be in $PATH:
- [Julia](https://julialang.org/) v1.2.0 above
- [StringTie](https://ccb.jhu.edu/software/stringtie/)
- [HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/)

Below Julia packages should be preinstalled:
- [SortingAlgorithms.jl](https://github.com/JuliaCollections/SortingAlgorithms.jl)

### Data
- Mapped bam files, sorted by coordinates;
- ENSEMBL GTF file;

## Usage
`julia [-p n] redentor.jl <parameters>`

### Parameters:
 `--bam (or -b) file1.bam,file2.bam,...` Input mapped bam file (sorted by coordinate). For multiple files input, seperate them with comma.
 `--id  (or -i) id1,id2,...` ID of each sample. Should be the same number as bam files.
 `--output (or -o)` The directory for the output. Should be pre-created.
 `--genome (or -g)` hg38|mm10 The sequenced genome type. Only support hg38 and mm10 so far.
 `--annotation (or -a)` file.gtf The transcription annotation file. Tested on Ensembl annotation.
 `--seqstrand (or -s)` first|second|none  The sequenced strand. Optional (default="none").
 `-p thread_n` Assign thread number. It is a Julia parameter and should be putted between `julia` and `redentor.jl`. Optional (default=1).

## Output
RAT_info.tsv, RAT_readcount.tsv, and RAT_RPKM.tsv are the information, read count, and read per kilobase of transcript, per million mapped reads of all detected RATs. For the RAT gene structure, see work/StringTieMerged.gtf.