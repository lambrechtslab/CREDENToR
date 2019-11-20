using Distributed
@everywhere include((@__DIR__)*"/libs/my.jl")

@everywhere module CREDENToR #A supportive module for biodata handling.
using Distributed
using DelimitedFiles #for writedlm
using Printf
using ..Youtiao

THREAD=nprocs()
BAMs=AbstractString[]
IDs=String[]
SEQSTRAND=:none
GTF=""
OUTPUT_PATH="."
REPEAT_FILE=""

#{{ main
function main()
    parse_arguments() && return nothing
    software_check()
    run_assemble_and_count()
    prepare_output()
    println("[CREDNEToR] All finished.")
end
#}}
#{{ parse_arguments
function parse_arguments()
    global THREAD, BAMs, IDs, SEQSTRAND, GTF, OUTPUT_PATH, REPEAT_FILE
    i=1
    while i<=length(ARGS)
        if ARGS[i][1]=='-'
            if ARGS[i]=="-b" || ARGS[i]=="--bam"
                BAMs=split(ARGS[i+1], ',')
                i+=2
            elseif ARGS[i]=="-i" || ARGS[i]=="--id"
                IDs=split(ARGS[i+1], ',')
                i+=2
            elseif ARGS[i]=="-o" || ARGS[i]=="--output"
                OUTPUT_PATH=ARGS[i+1]
                i+=2
            elseif ARGS[i]=="-g" || ARGS[i]=="--genome"
                REPEAT_FILE=(@__DIR__)*"/repeat_ann/"*ARGS[i+1]*"_rmsk.tsv.gz"
                i+=2
            elseif ARGS[i]=="-a" || ARGS[i]=="--annotation"
                GTF=ARGS[i+1]
                i+=2
            elseif ARGS[i]=="-s" || ARGS[i]=="--seqstrand"
                SEQSTRAND=Symbol(ARGS[i+1])
                i+=2
            elseif ARGS[i]=="-h" || ARGS[i]=="--help"
                print_help()
                return true
            elseif ARGS[i]=="-v" || ARGS[i]=="--version"
                print_version()
                return true
            else
                error("Invalid argument: "*ARGS[i])
            end
        else
            error("Invalid argument: "*ARGS[i])
            # i+=1
        end
    end
    length(BAMs)==length(IDs) || error("Bam files and IDs have not the same number, or any of them is not assigned.")
    if !isdir(OUTPUT_PATH)
        mkpath(OUTPUT_PATH)
        println("[CREDENToR] Cannot find output path $OUTPUT_PATH. Created.")
    end
    isfile(GTF) || error("Cannot find GTF file at "*GTF)
    in(SEQSTRAND, [:none, :first, :second]) || error("Invalid --seqstrand "*SEQSTRAND)
    for bam in BAMs
        isfile(bam) || error("Cannot find bam file "*bam)
    end
    isfile(GTF) || error("Cannot find GTF file "*GTF)
    isfile(REPEAT_FILE) || error("Cannot find repeat file "*REPEAT_FILE)
    return false
end
#}}
#{{ software_check
function software_check()
    success(`stringtie --version`) || error("Cannot find Samtools.")
    success(`htseq-count -h`) || error("Cannot find THSeq.")
end
#}}
#{{ run_assemble_and_count
function run_assemble_and_count()
    strandtag=SEQSTRAND==:second ? "--fr" : SEQSTRAND==:first ? "--rf" : ""
    mkpath("$OUTPUT_PATH/work/")
    if isfile("$OUTPUT_PATH/work/StringTieMerged.gtf")
        println("[CREDENToR] StringTieMerged.gtf is preexisted. Skip....")
    else
        for (bam, lb) in zip(BAMs, IDs)
            if isfile("$OUTPUT_PATH/work/$lb.StringTie.gtf")
                println("[CREDENToR] $lb.StringTie.gtf is preexisted. Skip....")
            else
                println("[CREDENToR] Assembling $lb annotation....")
                runstr("stringtie $bam -p $THREAD -G $GTF -o $OUTPUT_PATH/work/$lb.StringTie.gtf.tmp $strandtag")
                runstr("mv $OUTPUT_PATH/work/$lb.StringTie.gtf.tmp $OUTPUT_PATH/work/$lb.StringTie.gtf")
            end
        end
        println("[CREDENToR] Merging all assembled annotation....")
        writedlm("$OUTPUT_PATH/work/stringtie_output_list.txt", f"$OUTPUT_PATH/work/$1.StringTie.gtf".(IDs))
        runstr("stringtie --merge -p $THREAD -o $OUTPUT_PATH/work/StringTieMerged.gtf.tmp -G $GTF $OUTPUT_PATH/work/stringtie_output_list.txt")
        runstr("mv $OUTPUT_PATH/work/StringTieMerged.gtf.tmp $OUTPUT_PATH/work/StringTieMerged.gtf")
    end
    println("[CREDENToR] Counting reads ....")
    l=isfile.(f"$OUTPUT_PATH/work/$1-htseq.tsv".(IDs))
    for lb in IDs[l]
        println("[CREDENToR] $lb-htseq.tsv is preexisted. Skip....")
    end
    strandarg=SEQSTRAND==:second ? "yes" : SEQSTRAND==:first ? "reverse" : "no"
    (THREAD>1 ? pmapwd : foreach)(IDs[.!l], BAMs[.!l], env=(strandarg, OUTPUT_PATH)) do lb, bam, strandarg, OUTPUT_PATH
        runstr("htseq-count -f bam -s $strandarg -r pos $bam $OUTPUT_PATH/work/StringTieMerged.gtf > $OUTPUT_PATH/work/$lb-htseq.tsv.tmp")
        runstr("mv $OUTPUT_PATH/work/$lb-htseq.tsv.tmp $OUTPUT_PATH/work/$lb-htseq.tsv")
        nothing
    end
end
#}}
#{{ prepare_output
function prepare_output()
    println("[CREDNEToR] Preparing output....")
    T=parseGTF("$OUTPUT_PATH/work/StringTieMerged.gtf")
    rename!(T, "geneid"=>"gene_STRid")
    exn=parseGTF(GTF)
    D=Dict(zip(exn["transid"], exn["geneid"]))
    cdTransSet=Set(exn["transid"][exn["gene_typ"].=="protein_coding"])
    T["geneid"]=grpfunexp(T["gene_STRid"], T["transid"], T["gene_STRid"]) do x, y
        # l=startswith.(x, "ENST")
        l=.!startswith.(x, "MSTRG.")
        if any(l)
            ll=ismbr(x[l], cdTransSet)
            if any(ll)
                l=falsesbut(length(l), findall(l)[ll])
            end
            modebyr(map(xx->get(D, xx, y), x[l]))
        else
            y
        end
    end
    T["gene_typ"]=dtshift(T["geneid"], exn["geneid"], exn["gene_typ"], "novel_detected")

    rep=readtb(`zcat $REPEAT_FILE`)
    rep["chrno"]=chr2no(rep["chr"])

    ncT=rec(T, .!ismbr(T["gene_typ"], c"protein_coding"))
    ri, ti, ovlplen=genomemap(d"chrno, pos"rep, d"chrno, exn_pos"ncT, touch=true) do si, ti, sp, tp
        (si, ti, seglen([max(sp[1], tp[1]) min(sp[2], tp[2])]))
    end
    t=grpfun(ncT["gene_STRid"][ti], rep["rep_recNo"][ri], ovlplen) do RT, cov
        sumcov, uRT=grpfun(sum, RT, cov)
        uRT[argmax(sumcov)]
    end
    rt_ge=tb(c"rep_recNo, gene_STRid", t)
    d"RT_name, Fam, subFam, RT_strand"rt_ge=dtshift(rt_ge["rep_recNo"], rep["rep_recNo"], d"rep_name, rep_Fam, rep_subFam, strand"rep)
    T=tbx(T, rt_ge, "gene_STRid", fillempty=:B, order=:A, buniq=true)

    T["gene_group"]=map(T["gene_typ"].=="protein_coding", T["Fam"]) do x, y
        if x
            "cdgene"
        elseif isempty(y)
            "nonRT"
        else
            "RT_"*y
        end
    end
    
    ge=tbuniq(T, "gene_STRid", trim=true)
    d"gene_pos, gene_len"ge=grpfun(fastgrp(T["gene_STRid"], ge["gene_STRid"]), T["exn_pos"]) do x
        ([minimum(x[:, 1]) maximum(x[:, 2])],
         sum(seglen(segunion(x)[1])))
    end

    #Read HTSeq output
    C=map(IDs) do lb
        T=readtb("$OUTPUT_PATH/work/$lb-htseq.tsv", head=c"gene_STRid, N::readnum")
        rec(T, .!startswith.(T["gene_STRid"], "__"))
    end
    @assert iselequal(x->x["gene_STRid"], C)
    ge0=tb(gene_STRid=C[1]["gene_STRid"],
           readnum=hcat(map(x->x["readnum"], C)...))
    ge=tbx(ge0, ge, "gene_STRid")

    ge["RPKM"]=(1e+9)*ge["readnum"]./sum(ge["readnum"], dims=1)./ge["gene_len"]
    rt=rec(ge, startswith.(ge["gene_group"], "RT_"))
    
    writedlm("$OUTPUT_PATH/RAT_RPKM.tsv", asmx("RAT_ID"=>rt["gene_STRid"], IDs=>rt["RPKM"]))
    writedlm("$OUTPUT_PATH/RAT_readcount.tsv", asmx("RAT_ID"=>rt["gene_STRid"], IDs=>rt["readnum"]))
    writedlm("$OUTPUT_PATH/RAT_info.tsv",
             asmx("RAT_ID"=>rt["gene_STRid"],
                  "dominant_retrotransposon"=>rt["RT_name"],
                  "subfamily"=>rt["subFam"],
                  "family"=>rt["gene_group"],
                  "transcription_length"=>rt["gene_len"],
                  "associated_Ensembl_gene_id"=>map(x->startswith(x, "MSTRG") ? x : "", rt["geneid"]))
             )

end
#}}
#{{ print_help
function print_help()
    helpdoc="""
CREDENToR v0.1.0 usage:
julia [-p n] redentor.jl <parameters...>
A bioinformatics pipeline to detect and quantify cryptic transcripts associated with retrotranspon repeats.
Parameters:
--bam (or -b) file1.bam,file2.bam,... Input mapped bam file (sorted by coordinate). Multiple input files should be seperated with a comma.
--id  (or -i) id1,id2,... ID of each sample. Should be the same length as bam files.
--output (or -o) The directory for the output. Should be pre-created.
--genome (or -g) hg38|mm10 The sequenced genome type. Only hg38 and mm10 are supported so far.
--annotation (or -a) file.gtf The transcription annotation file. Tested on Ensembl annotation only.
--seqstrand (or -s) first|second|none  The sequenced strand. Optinal (default="none").
--version (or -v) Print version information.
--help (or -h) Print help.
-p thread_n Assign thread number. A Julia parameter that should be put between `julia' and `redentor.jl'. Optional (default=1).
"""
    println(helpdoc)
end
#}}
#{{ print_version
function print_version()
    println("CRETENDoR pipeline v0.1.0")
end
#}}
end

CREDENToR.main()
