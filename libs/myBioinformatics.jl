#{{ genomemap coordinatemap

#[[ genomemap coordinatemap ]]
#(SiteIdx, TagIdx)=genomemap(Site, Tag; fun=(si,ti,sp,tp), touch=false*, revstrand=false,
#                     minoverlap=0~1, no_match_output=(Int[],Int[]))
#Or ... = genomemap(fun, ...; ...)
#(SiteIdx, TagIdx)=coordinatemap(SitePos, TagPos; SiteId=index, TagId=index,
#                     sitesorted=false, tagsorted=false, sortcheck=true,
#                     minoverlap=0~1, no_match_output=(Int[],Int[]) )
#Or ... = coordinatemap(fun, ... ; ... )
# Assignment of `minoverlap=` implys `touch=true`.
#Site and Tag are two tuples with (Chromosome, pos[, strand]). `pos` is a n x 2 matrix. Chromosoem could be other group IDs.
#When `revstrand=ture`, function will reverse Site's strand in the first.
#If the Tags is sorted by (chr[, strand], pos) order, set tagsorted as true to save time.
#In non-touch mode, any tag fall in site completely will be reported. In touch mode, the overlapped_region / min(site_length, tag_length) >= minoverlap will be reported.
#The input of `fun` is (SiteIdx::Int, TagIdx::Int, SitePos::[1x2 Matrix], TagPos::[1x2 Matrix]). The function output will be pipeup by addrows!, except it be nothing or empty tuple ().
#In the case of no any match, function return `no_match_output`.
#`coordinatemap`: Mapping 1D coordinates without consider chromosome and strand.
#See also: coordinatemap
#Xiong Jieyi, 8 Jul 2014 >December 8, 2014>February 8, 2015>Dec 22, 2015>Jan 4, 2016>8 Mar 2017

# using my
export genomemap, coordinatemap
genomemap(Site::Tuple{Array,Matrix},Tag::Tuple{Array,Matrix}; kw...)=
    genomemap(tuple(Site...,Char[]),tuple(Tag...,Char[]); kw...)
function genomemap(Site::Tuple{Array{Tchr,N1},Matrix{Tpos},AbstractArray{Tstrand,N2}},
                   Tag::Tuple{Array{Tchr,N1},Matrix{Tpos},AbstractArray{Tstrand,N2}};
                   revstrand::Bool=false,no_match_output=(Int[],Int[]), kw...) where {Tchr,Tpos<:Real,Tstrand<:Char,N1,N2}
    
    @assert(length(Site)==length(Tag),"Input should both have or haven't strand information.")
    if isempty(Site[3]) && isempty(Tag[3])
        Site=Site[1:2]
        Tag=Tag[1:2]
    elseif revstrand
        Site=(Site[1],Site[2],dtshift(Site[3],['+','-'],['-','+']))
    end

    #Remove NaN
    SiteNo=collect(1:size(Site[2],1))
    TagNo=collect(1:size(Tag[2],1))

    l=.!isnan.(Site[2][:,1])
    Site=getrow(Site,l)
    SiteNo=SiteNo[l]

    l=.!isnan.(Tag[2][:,1])
    Tag=getrow(Tag,l)
    TagNo=TagNo[l]
    
    if length(Site)==2 #no strand
        tfg=fastgrp(Tag[1])
        sfg=fastgrp(Site[1],tfg.grpid)
    else
        tfg=fastgrp((Tag[1],Tag[3]))
        sfg=fastgrp((Site[1],Site[3]),tfg.grpid)
    end
    rlt=rowpile()
    for gi=1:sfg.grpnum
        ls=want(sfg,gi)
        if isempty(ls)#Add in 28 Aug 2014
            continue
        end
        lt=want(tfg,gi)
        crlt=coordinatemap(Site[2][ls,:], Tag[2][lt,:], SiteId=SiteNo[ls], TagId=TagNo[lt];
                           no_match_output=nothing, kw...)
        if crlt!=nothing
            addrows!(rlt, crlt)
        end
    end

    isvirgin(rlt) ? no_match_output : value(rlt)
end
genomemap(Site::Tuple{Array,Matrix,AbstractArray},Tag::Tuple{Array,Matrix,AbstractArray};args...)=
    genomemap(promote_tuple(Site,Tag)...;args...)
genomemap(fun::Function,args...; kw...)=genomemap(args...; kw..., fun=fun)

function coordinatemap(Site::Matrix{Ts},Tag::Matrix{Tt};
                       SiteId::Group=collect(1:size(Site,1)),TagId::Group=collect(1:size(Tag,1)),
                       sitesorted::Bool=false,tagsorted::Bool=false,sortcheck::Bool=true,
                       fun::Function=(si,ti,sp,tp)->(si,ti),
                       minoverlap::AbstractFloat=-2.0,touch::Bool=minoverlap>-1.0,
                       no_match_output=(Int[],Int[])) where {Ts<:Real,Tt<:Real}
    if sitesorted
        if sortcheck && !issortedr(Site)
            error("Site is not sorted.")
        end
    else
        (Site,t)=sortr(Site)
        SiteId=getrow(SiteId,t)
    end
    if tagsorted
        if sortcheck && !issortedr(Tag)
            error("Tag is not sorted.")
        end
    else
        (Tag,t)=sortr(Tag)
        TagId=getrow(TagId,t)
    end

    rlt=rowpile()
    ps=1
    pt=1
   if touch
        #Touch model
        while ps<=size(Site,1) && pt<=size(Tag,1)
            if Tag[pt,2]<Site[ps,1]
                pt+=1
            elseif Site[ps,2]<Tag[pt,1]
                ps+=1
            else
                oldpt=pt
                while true
                    if Site[ps,1]<=Tag[pt,2] #Hitted
                        sitepos=Site[ps,:]
                        tagpos=Tag[pt,:]
                        if (minoverlap<=0.0) || ( (min(sitepos[2],tagpos[2])-max(sitepos[1],tagpos[1])+1) / (min(sitepos[2]-sitepos[1],tagpos[2]-tagpos[1])+1) )>=minoverlap
                            crlt=fun(getrow(SiteId,ps),getrow(TagId,pt),sitepos,tagpos)
                            if crlt!=nothing && crlt!=()
                                addrows!(rlt,crlt)
                            end
                        end
                    end
                    pt+=1
                    if pt>size(Tag,1) || Site[ps,2]<Tag[pt,1]
                        break
                    end
                end
                pt=oldpt
                ps+=1
            end
        end
    else
        #Match model
        while ps<=size(Site,1) && pt<=size(Tag,1)
            if Tag[pt,1]<Site[ps,1]
                pt+=1
            elseif Site[ps,2]<Tag[pt,1]
                ps+=1
            else
                oldpt=pt
                while true
                    if Tag[pt,2]<=Site[ps,2]#Hitted
                        crlt=fun(getrow(SiteId,ps),getrow(TagId,pt),Site[ps,:],Tag[pt,:])
                        # if !isa(crlt,@compat(Tuple{}))
                        #     # addrow!(rlt,crlt...)
                        #     addrows!(rlt,crlt...) #Changed in February 8, 2015
                        # end
                        if crlt!=nothing && crlt!=()
                            addrows!(rlt,crlt)
                        end
                    end
                    pt+=1
                    if pt>size(Tag,1) || Site[ps,2]<Tag[pt,1]
                        break
                    end
                end
                pt=oldpt
                ps+=1
            end
        end
    end
    # out=value(rlt)
    # return out==nothing?no_match_output:out
    isvirgin(rlt) ? no_match_output : value(rlt)
 end
coordinatemap(fun::Function, args...;kw...)=coordinatemap(args...; kw..., fun=fun)

#}}

#{{ chr2no, no2chr, no2chr_celeg
#[[ chr2no no2chr no2chr_celeg ]]
# chrNO|chrNO = chr2no(chrID|chrID_vector)
# chrId = no2chr(chrNO; chr=TF|prefix="chr") | no2chr_celeg(chrNo; chr=TF|prefix="CHROMOSOME_") #C_elegans style: CHROMOSOME_III
#Convert chrID to chrNO or NO. NN=chrNN, 31=chrX, 32=chrY, 41=chrM( or 41=MT when chr=false or prefix is empty), 51=chr2a, 51=chr2b, 0=all others.
#See also: genomemap
#Xiong Jieyi, 13 Jul 2014>February 25, 2015>31 Aug 2019

export chr2no, no2chr, no2chr_celeg
function chr2no(chr::AbstractString)
    chr=replace(chr, r"^CHROMOSOME_|^Chr|^chr|^CHR" => "")
    Tb=Dict{AbstractString,Int}("I"=>1,"II"=>2,"III"=>3,"IV"=>4,"V"=>5,
        "X"=>31,"Y"=>32,
        "M"=>41,"MT"=>41,"Mt"=>41,"mt"=>41,"mtDNA"=>41,"MTDNA"=>41,"MtDNA"=>41,
        "2a"=>51,"2A"=>51,"2b"=>52,"2B"=>52)
    chrno=get(Tb,chr,0)
    if chrno==0
        try
            chrno=int(chr)
        catch
            chrno=0
        end    
    end
    return chrno
end
chr2no(chr::Vector{T}) where {T<:Union{AbstractString,Symbol}} = map(chr2no,chr)
chr2no(chr::Symbol)=chr2no(string(chr))
function no2chr(no::Real; chr::Bool=true, prefix::AbstractString=chr ? "chr" : "")
    Tb=Dict{Int,AbstractString}(31=>"X", 32=>"Y", 41=>isempty(prefix) ? "MT" : "M" ,
                                51=>"2a",52=>"2b")
    chr=get(Tb,no,"")
    if chr==""
        chr=string(int(no))
    end
    return prefix*chr
end
no2chr(chr::Vector{T};args...) where {T<:Real} = map(x->no2chr(x;args...),chr)
function no2chr_celeg(no::Real; chr::Bool=true, prefix::AbstractString=chr ? "CHROMOSOME_" : "")
    Tb=Dict{Int,AbstractString}(31=>"X",32=>"Y",41=>isempty(prefix) ? "MT" : "M",
                                1=>"I",2=>"II",3=>"III",4=>"IV",5=>"V")
    chr=get(Tb,no,"")
    if chr==""
        error("Unknown chrNO $no.")
    end
    return prefix*chr
end
no2chr_celeg(chr::Vector{T};args...) where {T<:Real} = map(x->no2chr_celeg(x;args...),chr)
#}}

#{{ seglen
#[[ seglen ]]
# lengths=seglen(pos)
#Equal to pos[:,2]-pos[:,1]+1, but seglen([0 0])=0.
#A error will be reported if any segment is less than 1nt.
#See also: chr2no, no2chr
#Xiong Jieyi, 11 Sep 2014>5 Oct 2014

#<v0.6# function seglen{T<:Real}(pos::Array{T})
function seglen(pos::Array{T}) where {T<:Real}
    l=.!isnan.(pos[:,1])
    all(pos[.!l,1].<=pos[.!l,2]) || error("Invalid position number.")
    len=pos[:,2]-pos[:,1]+1
    len[(pos[:,1].==0) .& (pos[:,2].==0)].=0
    return len
end
export seglen
#}}

#{{ readGeneralSamFile, readGeneralSamLines
#[[ readGeneralSamFile readGeneralSamLines ]]
# Table = readGeneralSamLines(iterator; pairend=false, strand=true, readidfun=ascii, chrfun=x->(chr2no(x)|>x->x==0?nothing:Int16(x)), recordno_ref::Vector{Int}=[0], fields::Dict=tb(...), fieldfun::Dict=tb(cs=parsefun))
# ... = readGeneralSamFile(alignment.sam|alignment.bam|`pipe`; ...)
# Input:
#   iterator: Any iteratiable data, each item is a line of sam file (no newline in the end).
#   pairend: Paired reads?
#   strand:  Read has strand information?
#   readidfun and chrfun: function to handle readid and chr. Input is a SubString and output must be a row or nothing. if function output nothing, this record will be ignored.
#   recordno_ref: A 1-element Int vector to assign and access the start and end recordno.
#   fields: a one row table, keys are the needed optional fields in sam file, and values are the default values for the missing fields.
#   fieldfun: a dict, keys are the fields and values are parsing functions. Note that the input of this funciton are always a string.
#
# Output field:
#  record_no, readid, flag, chrno, pos, juncno, read_len, ...
#
#See also: readStarSamFile
#Xiong Jieyi, 12 Apr 2018 > 16 May 2018 >30 May 2018

export readGeneralSamFile, readGeneralSamLines
function readGeneralSamLines(samlines;
                             pairend::Bool=false,
                             strand::Bool=true,
                             readidfun::Function=ascii,
                             chrfun::Function=x->(chr2no(x)|>x->x==0 ? nothing : Int16(x)),
                             recordno_ref::Vector{Int}=[0],
                             fields::Dict=tb(),
                             fieldfun::Dict=tb())
    RECORD_NO=recordno_ref[1]
    rlt=rowpile()
    rlt_fields=rowpile()
    
    for ctxt in samlines
        tc=split(chomp(ctxt),'\t')
        if length(tc)<3
            @warn("Ignore invalid line (record NO: $RECORD_NO): \"$ctxt\"")
            continue
        end
        if tc[3]=="*"
            continue
        end
        CHRNO=chrfun(tc[3])
        if CHRNO==nothing
            continue
        end
        
        READID=readidfun(tc[1])
        if READID==nothing
            continue
        end
        
        RECORD_NO+=1
        FLAG=parse(Int16,tc[2])
        PAIRNO=pairend ? (FLAG&128>0)+1 : 0
        READLEN=Int32(length(tc[10])) #Add in 23 Sep 2014

        #Parse Appendences
        ADDFDS=deepcopy(fields)
        
        for ctc=tc[12:end]
            fdnm, fdty, fdval = split(ctc,':', limit=3)
            fdfun=get(fieldfun, fdnm, nothing)
            if fdfun==nothing
                fdfun=if fdty=="i"
                    int
                elseif fdty=="f"
                    float
                elseif fdty=="Z"
                    ascii
                elseif fdty=="A" #Char
                    x->x[1]
                else
                    error("Unknown field type: $fdty")
                end
            end
            
            if haskey(ADDFDS, fdnm)
                ADDFDS[fdnm]=fdfun(fdval)
            end
        end
        
        #Parse CIGAR
        STR=int(tc[4])
        POS=zeros(Int,1,2)
        LOC=zeros(Int,1,2)
        cigar=tc[6];p=1;rdp=1;rfp=STR;
        aL=findall(['A'<=x<='Z' for x in cigar])
        LOC_read_start=[];
        cjuncno=Int16(0);
        for cL=aL
            NO=int(cigar[p:cL-1]);p=cL+1;
            ccL=cigar[cL]
            if ccL=='S' || ccL=='H'
                rdp+=NO
                #S is not be calculated in the mapping start position.
                ## Remove this: ## rfp+=NO
                #Changed in Dec 14, 2015
            elseif ccL=='M'
                if !(POS[1]>0)
                    LOC[1]=rdp
		    if isempty(LOC_read_start)
		        LOC_read_start=rdp
		    end
                    POS[1]=rfp
                end
                rfp+=NO
                rdp+=NO
                LOC[2]=rdp-1
                POS[2]=rfp-1
            elseif ccL=='D'
                rfp+=NO;
            elseif ccL=='I'
                rdp+=NO
            elseif ccL=='N'
	        cjuncno+=1
                addrow!(rlt,RECORD_NO,READID,FLAG,CHRNO,POS,cjuncno,READLEN,PAIRNO)
                addrow!(rlt_fields, ADDFDS)
                # SEQ="" #For junction reads, only the first one assigned with sequence.
                # CIGAR4save=""
                rfp+=NO
                POS=zeros(1,2)
            else
                error("Unknown cigar char in read: $ccL")
            end
        end
        if cjuncno>0
            cjuncno+=1
        end
        addrow!(rlt,RECORD_NO,READID,FLAG,CHRNO,POS,cjuncno,READLEN,PAIRNO)
        addrow!(rlt_fields, ADDFDS)
    end
    rltvalue=value(rlt)
    if rltvalue==nothing
        return nothing
    end
    read=tb(value(rlt_fields), c"record_no,readid,flag,chrno,pos,juncno,read_len,pairno"=>rltvalue)
    if !pairend
        delete!(read, "pairno")
    end
    if strand
        read["strand"]=map(x->x&16==0 ? '+' : '-',read["flag"])
        if pairend
            #Reverse second half-read in pair-end mode
            read["strand"][read["pairno"].==2]=map(x->x=='-' ? '+' : '-',read["strand"][read["pairno"].==2])
        end
    end
    recordno_ref[1]=RECORD_NO
    return read
end
function readGeneralSamFile(file::Union{AbstractString, Cmd}, args...; wargs...)
    if isa(file,AbstractString) && splitext(file)[2]==".bam"
        readGeneralSamLines(eachline(sfc"samtools view $1"(file)), args...; wargs...)
    else
        readGeneralSamLines(eachline(file), args...; wargs...)
    end
end
#}}

#{{ readStarSamFile readStarSamLine
#[[ readStarSamFile readStarSamLine ]]
# read=readStarSamFile(filename|IOStream; savefile="", onlyuniq=true,pairend=false, strand=true, chrfun=int16(chr2no), keepHI=false, onlyjunc=false,
#                      [readprefixlen=auto, readiddm=auto] | readidfun=fun(Field1_Str),
#                      region=(chr|chrno|vec, [posstart posend])) or (chr|chrno|vec,),
#                      addfields=(c"name1,name2,...", [N1, N2, ...])
#Read Star/Tophat output sam/bam file (depended by its ext name.).
#
# read|nothing=readStarSamLine(String/String_Vector; recordno_ref=[0], pairend=false, strand=true, readprefixlen=auto | readidfun=fun(Field1_Str), chrfun=int16(chr2no), addfields=...)
#Read lines of sam file.
#
#Output field:
#  record_no, readid, flag, chrno, pos, juncno, read_len,
#  align_score(AS), mapped_loci_num(NH), alignment_index(HI, Only for STAR), mismatch_num(nM in STAR or NM in Tophat).
#
#In readStarSamFile, you can use addfields to add other fields in SAM file. for example, readStarSamFile(..., addfields=(c"sequence,CIGAR", [9,6])).
#
#If onlyuniq is true, function will remove low-score multiple mapping reads. Note that the multiple mapped reads with equal highest score will still be keeped. You can pick up these reads by is_uniq_map field.
#In pair-end mode, paired reads just be regarded as two indipendent reads, with field pairno=1 or 2. If strand is needed, the strand of the second half pair will be reversed to make sure both half pairs have consistent strand information.
#If onlyjunc is ture, function will only keep junction reads (e.g. for alternative splicing analysis). Note that as the non-junction reads will be ignored before calculate uniquity, the duplicated reads selection (onlyuniq=true) will only considered among junction reads, and the is_uniq_map field will just be replaced by mapped_loci_num.==1.
#readprefixlen assigned the id of reads which will be removed. It will be detected automaticly if omitted. readiddm is the delimiter. If ignored, function will check if readid has ':', '.' or '_' (and error if none is found.)
#Read will be sorted as (chrno,strand,pos).
#chrfun: The function convert chromosome field (#3) to a number. If the number is zero, this record will be discarded.
#
#region only work in indexed bam file. region can be single one or mulitple regions. In the case of multiple regions, function will segunion them first to avoid duplicated output records (For the reasion, see samtools manual.). If input is chrno rather than chr-name, function will check the first record of bam file first, to determind which prefix ("chr", "CHR" or "") should be used for no2chr. If only chromosome is given, function will read all the reads in given chromosome(s).
#
#For readStarSamFile, recordno_ref should be a reference (one element vector) for record number counter. Use [0] for the start value. If it is missed, the read in all output will start with 1.
#
#When no any valid read was found, readStarSamFile throw an error, while readStarSamLine return nothing.
#
#See also: readGeneralSamFile, readGeneralSamLines, seglen, segconc, genomemap, parseGTF, parseDemulReport, parseBamMutation
#Xiong Jieyi, 13 Jul 2014>29 Jul 2014>25 Aug 2014>23 Sep 2014>January 12, 2015>May 17, 2015>Jun 12, 2015>Aug 4, 2015>Aug 21, 2015>Aug 23, 2015>Nov 23, 2015>Dec 10, 2015>Dec 14, 2015>Apr 19, 2016>Apr 25, 2016>May 17, 2016>Sep 9, 2016

function readStarSamFile(fp::IO;
                         pairend::Bool=false,
                         savefile::AbstractString="",
                         onlyuniq::Bool=true,
                         strand::Bool=true,
                         readprefixlen::Int=-1,readiddm='\0',
                         readidfun::Union{Function,Void}=nothing,
                         chrfun::Function=x->Int16(chr2no(x)),
                         keepHI::Bool=false,
                         addfields::Tuple{Vector,Vector}=([],[]),
                         onlyjunc::Bool=false )#Add in Sep 9, 2016

    isaddfields=!isempty(addfields[1])
    rlt=rowpile()
    ctxt=chomp(readline(fp))
    if isempty(ctxt)
        @warn("File is empty.")
        return nothing
    end
    while ctxt[1]=='@'
        ctxt=chomp(readline(fp))
    end
    RECORD_NO=0
    while true
        tc=split(ctxt,'\t')
        if length(tc)<3
            @warn("Ignore invalid line (record NO: $RECORD_NO): \"$ctxt\"")
            if eof(fp)
                break
            end
            ctxt=chomp(readline(fp))
            continue
        end
        if tc[3]=="*" || (CHRNO=chrfun(tc[3]))==0
            if eof(fp)
                break
            end
            ctxt=chomp(readline(fp))
            continue
        end

        #Parse READID
        READID=nothing
        try
            if isa(readidfun,Void)
                if readiddm=='\0'
                    for t in [':','.','_']
                        if in(t,tc[1])
                            readiddm=t
                            break
                        end
                    end
                end
                in(readiddm,tc[1])||error("readiddm is not exist in readid.")
                if readprefixlen<0
                    t=match(Regex("\\$readiddm[\\$readiddm\\d]+\$"),tc[1])
                    # t=match(r"\:[\:\d]+$",tc[1])
                    # if t==nothing
                    #     readprefixlen=findonly(collect(tc[1]).=='.')
                    #     readidfun=X->parse(Int32,X[readprefixlen+1:end])
                    # else
                    readprefixlen=length(tc[1])-length(t.match)+1
                    readidfun=X->map(x->parse(Int32,x),split(X[readprefixlen+1:end],readiddm))
                    # end
                else
                    @assert(tc[1]==0 || tc[1][readprefixlen]==readiddm || tc[1][readprefixlen]=='.', "Inproperate readprefixlen. Try to set readidfun and readprefixlen for the specific SAM file.")
                    if tc[1][readprefixlen]==readiddm
                        readidfun=X->map(x->parse(Int32,x),split(X[readprefixlen+1:end],':'))
                    else
                        # readidfun=X->parse(Int32,X[readprefixlen+1:end])
                        error("readprefixlen($readprefixlen) is not point to a readiddm($readiddm)")
                    end
                end
            end
            READID=readidfun(tc[1])
        catch err
            println("Error in parsing read-ID with prefix length $readprefixlen: "*tc[1])
            rethrow(err)
        end
        
        RECORD_NO+=1
        FLAG=parse(Int16,tc[2])
        if pairend
            PAIRNO=Int8((FLAG&128>0)+1)
        else
            PAIRNO=Int8(0)
            # READID=[PAIRNO;READID]
        end
        READLEN=Int32(length(tc[10])) #Add in 23 Sep 2014

        if isaddfields
            ADDFDS=map(x->ascii(tc[x]),addfields[2])
        else
            ADDFDS=ASCIIString[]
        end
        
        #Parse Appendences
        AS=Int32(-1)
        NH=Int16(-1);HI=Int16(-1);nM=Int16(-1)
        for ctc=tc[12:end]
            if ctc[1:5]=="AS:i:"
                AS=parse(Int32,ctc[6:end])
            elseif ctc[1:5]=="NH:i:"
                NH=parse(Int16,ctc[6:end])
            elseif ctc[1:5]=="HI:i:" #Tophat have no HI.
                HI=parse(Int16,ctc[6:end])
            elseif ctc[1:5]=="nM:i:" || ctc[1:5]=="NM:i:" #STAR using nM for mismatch number, but Tophat use NM.
                nM=parse(Int16,ctc[6:end])
            end
        end
        
        #Parse CIGAR
        STR=int(tc[4])
        POS=zeros(Int,1,2)
        LOC=zeros(Int,1,2)
        # t=zeros(1,MAXMIS)
        # INS=copy(t);DEL=copy(t);INSLEN=copy(t);DELLEN=copy(t)
        # Ip=1;Dp=1;
        cigar=tc[6];p=1;rdp=1;rfp=STR;
        # aL=find(['A'.<=x.<='Z' for x in cigar])
        aL=findall(['A'<=x<='Z' for x in cigar])
        LOC_read_start=[];
        cjuncno=Int16(0);
        for cL=aL
            NO=int(cigar[p:cL-1]);p=cL+1;
            ccL=cigar[cL]
            if ccL=='S' || ccL=='H' #Added H in Dec 14, 2015. H used in BWA.
                rdp+=NO
                #S is not be calculated in the mapping start position.
                ## Remove this: ## rfp+=NO
                #Changed in Dec 14, 2015
            elseif ccL=='M'
                if !(POS[1]>0)
                    LOC[1]=rdp
		    if isempty(LOC_read_start)
			LOC_read_start=rdp
		    end
                    POS[1]=rfp
                end
                rfp+=NO
                rdp+=NO
                LOC[2]=rdp-1
                POS[2]=rfp-1
            elseif ccL=='D'
                # if Dp<=MAXMIS
                #     DEL[Dp]=rdp
                #     DELLEN[Dp]=NO
                # end
                # Dp+=1
                rfp+=NO;
            elseif ccL=='I'
                # if Ip<=MAXMIS
                #     INS[Ip]=rdp
                #     INSLEN[Ip]=NO
                # end
                # Ip+=1
                rdp+=NO
            elseif ccL=='N'
		cjuncno+=1
                addrow!(rlt,RECORD_NO,READID,PAIRNO,FLAG,CHRNO,POS,cjuncno,AS,NH,HI,nM,READLEN,ADDFDS...)
                rfp+=NO
                POS=zeros(1,2)
            else
                error("Unknown cigar char in read: $ccL")
            end
        end
        if cjuncno>0
            cjuncno+=1
        end
        if (!onlyjunc) || cjuncno>0
            addrow!(rlt,RECORD_NO,READID,PAIRNO,FLAG,CHRNO,POS,cjuncno,AS,NH,HI,nM,READLEN,ADDFDS...)
        end
        if eof(fp)
            break
        else
            ctxt=chomp(readline(fp))
        end
        
        if mod(RECORD_NO,1000000)==0
            println("Reading record: ", RECORD_NO);
        end
    end

    rltvalue=value(rlt)
    @assert(rltvalue!=nothing,"No any valid records. Please check, e.g., if the chromosome fields (#3) is valid.")
    read=tb([c"record_no,readid,pairno,flag,chrno,pos,juncno,align_score,mapped_loci_num,alignment_index,mismatch_num,read_len";addfields[1]],rltvalue)

    # Tophat and Bowtie are not support HI field.
    if !keepHI && all(read["alignment_index"]==-1)
        delete!(read,"alignment_index")
    end
    
    if onlyuniq
        println("Finding unique-mapped reads....")
        read=rec(read,sortri((read["readid"],read["pairno"],read["record_no"],-read["align_score"])))
        read["is_uniq_map"]=fill(true,rownum(read))
        lkeep=fill(true,rownum(read))

        preadid=pairend ? d"readid,pairno"read : read["readid"]
        cscore=read["align_score"][1]
        flag=true
        recordstrp=1
        for i=2:rownum(read)
            # if read["readid"][i,:]==read["readid"][i-1,:]
            if getrow(preadid,i)==getrow(preadid,i-1)
                if read["record_no"][i,:]==read["record_no"][i-1,:]
                    read["is_uniq_map"][i]=read["is_uniq_map"][i-1]
                    lkeep[i]=lkeep[i-1]
                else
                    if flag && read["align_score"][i]==cscore
                        read["is_uniq_map"][recordstrp:i].=false
                    else
                        flag=false
                        lkeep[i]=false
                    end
                    recordstrp=i
                end
            else
                flag=true
                recordstrp=i
                cscore=read["align_score"][i]
            end
        end
        read=rec(read,lkeep)
        if onlyjunc
            read["is_uniq_map"]=read["mapped_loci_num"].==1
        end
    end


    if !pairend
        delete!(read,"pairno")
    end

    if strand
        read["strand"]=map(x->x&16==0 ? '+' : '-',read["flag"])
        if pairend
            #Reverse second half-read in pair-end mode
            read["strand"][read["pairno"].==2]=map(x->x=='-' ? '+' : '-',read["strand"][read["pairno"].==2])
        end
        read=rec(read,sortri(d"chrno,strand,pos"read))
    else
        read=rec(read,sortri(d"chrno,pos"read))        
    end
    if !isempty(savefile)
        println("Saving file to $savefile ....")
        jldsave(savefile,read)
    end
    return read
end
function readStarSamFile(filename::Union{AbstractString,Cmd};region::Tuple=(),param...)
    if isa(filename,AbstractString) && ismatch(r"\.bam$",filename)
        if isempty(region)
            file=`samtools view $filename`
            println("Reading BAM file ....")
        else
            #Check if bam file is indexed
            if !isfile(filename[1:end-1]*"i") && !isfile(filename*".bai")
                println("Indexing bam file $filename")
                run(`samtools index $filename`)
                println("Indexing done.")
            end
            if length(region)==1
                #Read all reads by chromosomes
                chrs=isa(region[1],Group) ? unival(region[1]) : [region[1]]
                if eltype(chrs)<:Real
                    #Check if chromosome in bam file is started with chr
                    S=readchomp(pipeline(`samtools view -H $filename`,`grep '^@SQ'`))
                    M=match(r"\bSN\:(chr|CHR|Chr|chromosome|CHROMOSOME|Chromosome)",S)
                    if M!=nothing
                        chrs=no2chr(chrs,prefix=M.captures[1])
                    else
                        chrs=no2chr(chrs,prefix="")
                    end
                end
                T,=rowfunwith(chrs) do chr
                    file=`samtools view $filename $chr`
                    println("Reading BAM file chromosome $chr ....")
                    open((x)->readStarSamFile(x;keepHI=true,param...),file)
                end
                return T
            else
                #Read reads in some regions
                if isa(region[1],Group)
                    region,=segunion(region)
                else
                    region=([region[1]],region[2])
                end
                if eltype(region[1])<:Real
                    #Check if chromosome in bam file is started with chr
                    # fstbamli=open(readline,pipeline(`samtools view $filename`,`head -n1`))
                    # S=split(fstbamli,'\t')[3]|>x->x[1:min(3,length(x))]
                    # if lowercase(S)=="chr"
                    #     region=(no2chr(region[1],prefix=S), Array{Int}(region[2]))
                    # else
                    #     region=(no2chr(region[1],prefix=""), Array{Int}(region[2]))
                    # end
                    #Check if chromosome in bam file is started with chr
                    S=readchomp(pipeline(`samtools view -H $filename`,`grep '^@SQ'`))
                    M=match(r"\bSN\:(chr|CHR|Chr|chromosome|CHROMOSOME|Chromosome)",S)
                    if M!=nothing
                        region=(no2chr(region[1],prefix=M.captures[1]), Array{Int}(region[2]))
                    else
                        region=(no2chr(region[1],prefix=""), Array{Int}(region[2]))
                    end
                end
                T,=rowfunwith(region) do crg
                    rgstr="$(crg[1]):$(crg[2][1])-$(crg[2][2])"
                    file=`samtools view $filename $rgstr`
                    println("Reading BAM file region $rgstr ....")
                    open((x)->readStarSamFile(x;keepHI=true,param...),file)
                end
                return T
            end
        end
    else
        @assert(isempty(region) || isempty(region[1]),"region parameter only support bam file so far.")
        file=filename
        println("Reading SAM file ....")
    end
    open((x)->readStarSamFile(x;param...),file)
end
export readStarSamFile

#readStarSamLine
export readStarSamLine
readStarSamLine(samline::AbstractString,args...;wargs...)=readStarSamLine([samline],args...;wargs...)
function readStarSamLine(samlines::Vector{T};
                         pairend::Bool=false,
                         strand::Bool=true,
                         readprefixlen::Int=-2,
                         readidfun::Union{Function,Void}=nothing,
                         chrfun::Function=x->Int16(chr2no(x)),
                         recordno_ref::Vector{Int}=[0],
                         addfields::Tuple{Vector,Vector}=([],[])) where {T<:AbstractString}
    RECORD_NO=recordno_ref[1]
    isaddfields=!isempty(addfields[1])
    rlt=rowpile()
    # SEQ=""
    # CIGAR4save=""
    for ctxt in samlines
        tc=split(chomp(ctxt),'\t')
        if length(tc)<3
            @warn("Ignore invalid line (record NO: $RECORD_NO): \"$ctxt\"")
            continue
        end
        if tc[3]=="*" || (CHRNO=chrfun(tc[3]))==0
            continue
        end
        
        #Parse READID
        if isa(readidfun,Void)
            if readprefixlen<0
                t=match(r"\:[\:\d]+$",tc[1])
                if t==nothing
                    readprefixlen=findonly(collect(tc[1]).=='.')
                    readidfun=X->parse(Int32,X[readprefixlen+1:end])
                else
                    readprefixlen=length(tc[1])-length(t.match)+1
                    readidfun=X->map(x->parse(Int32,x),split(X[readprefixlen+1:end],':'))
                end
            else
                @assert(tc[1]==0 || tc[1][readprefixlen]==':' || tc[1][readprefixlen]=='.', "Inproperate readprefixlen. Try to set readidfun and readprefixlen for the specific SAM file.")
                if tc[1][readprefixlen]==':'
                    readidfun=X->map(x->parse(Int32,x),split(X[readprefixlen+1:end],':'))
                else
                    readidfun=X->parse(Int32,X[readprefixlen+1:end])
                end
            end
        end
        READID=readidfun(tc[1])
        
        # if isa(readidfun,Void)
        #     if readprefixlen<0
        #         if readprefixlen<-1
        #             error("readprefixlen must be set in readStarSamLine. If you want it be signed automatically, set it as -1.")
        #         end
        #         readprefixlen=length(tc[1])-length(match(r"\:[\:\d]+$",tc[1]).match)+1
        #     end
        #     # READID=readidfun(tc[1],readprefixlen,readsplitchar)
        #     @assert(tc[1]==0 || tc[1][readprefixlen]==':', "Inproperate readprefixlen. Try to set readidfun and readprefixlen for the specific SAM file.")
        #     READID=map(x->parse(Int32,x),split(tc[1][readprefixlen+1:end],':'))
        # else
        #     if readprefixlen>0
        #         tc1=tc[1][readprefixlen+1:end]
        #     else
        #         tc1=tc[1]
        #     end
        #     READID=readidfun(tc1)
        # end
        
        RECORD_NO+=1
        # if getseq
        #     SEQ=ascii(tc[10])
        # end
        # if getCIGAR
        #     CIGAR4save=ascii(tc[6])
        # end
        FLAG=parse(Int16,tc[2])
        if pairend
            PAIRNO=(FLAG&128>0)+1
            READID=[PAIRNO;READID]
        end
        READLEN=Int32(length(tc[10])) #Add in 23 Sep 2014

        if isaddfields
            ADDFDS=map(x->ascii(tc[x]),addfields[2])
        else
            ADDFDS=ASCIIString[]
        end
        
        #Parse Appendences
        AS=Int32(-1)
        NH=Int16(-1);HI=Int16(-1);nM=Int16(-1)
        for ctc=tc[12:end]
            if ctc[1:5]=="AS:i:"
                AS=parse(Int32,ctc[6:end])
            elseif ctc[1:5]=="NH:i:"
                NH=parse(Int16,ctc[6:end])
            elseif ctc[1:5]=="HI:i:" #Tophat have no HI.
                HI=parse(Int16,ctc[6:end])
            elseif ctc[1:5]=="nM:i:" || ctc[1:5]=="NM:i:" #STAR using nM for mismatch number, but Tophat use NM.
                nM=parse(Int16,ctc[6:end])
            end
        end
        
        #Parse CIGAR
        STR=int(tc[4])
        POS=zeros(Int,1,2)
        LOC=zeros(Int,1,2)
        cigar=tc[6];p=1;rdp=1;rfp=STR;
        aL=findall(['A'<=x<='Z' for x in cigar])
        LOC_read_start=[];
        cjuncno=Int16(0);
        for cL=aL
            NO=int(cigar[p:cL-1]);p=cL+1;
            ccL=cigar[cL]
            if ccL=='S' || ccL=='H'
                rdp+=NO
                #S is not be calculated in the mapping start position.
                ## Remove this: ## rfp+=NO
                #Changed in Dec 14, 2015
            elseif ccL=='M'
                if !(POS[1]>0)
                    LOC[1]=rdp
		    if isempty(LOC_read_start)
		        LOC_read_start=rdp
		    end
                    POS[1]=rfp
                end
                rfp+=NO
                rdp+=NO
                LOC[2]=rdp-1
                POS[2]=rfp-1
            elseif ccL=='D'
                rfp+=NO;
            elseif ccL=='I'
                rdp+=NO
            elseif ccL=='N'
	        cjuncno+=1
                addrow!(rlt,RECORD_NO,READID,FLAG,CHRNO,POS,cjuncno,AS,NH,HI,nM,READLEN,SEQ,CIGAR4save)
                # SEQ="" #For junction reads, only the first one assigned with sequence.
                # CIGAR4save=""
                rfp+=NO
                POS=zeros(1,2)
            else
                error("Unknown cigar char in read: $ccL")
            end
        end
        if cjuncno>0
            cjuncno+=1
        end
        addrow!(rlt,RECORD_NO,READID,FLAG,CHRNO,POS,cjuncno,AS,NH,HI,nM,READLEN,ADDFDS...)
        # addrow!(rlt,RECORD_NO,READID,FLAG,CHRNO,POS,cjuncno,AS,NH,HI,nM,READLEN,SEQ,CIGAR4save)
        # SEQ=""
        # CIGAR4save=""
    end
    rltvalue=value(rlt)
    if rltvalue==nothing
        return nothing
    end
    read=tb([c"record_no,readid,flag,chrno,pos,juncno,align_score,mapped_loci_num,alignment_index,mismatch_num,read_len";addfields[1]],rltvalue)
    # read=tb(c"record_no,readid,flag,chrno,pos,juncno,align_score,mapped_loci_num,alignment_index,mismatch_num,read_len,sequence,CIGAR",rltvalue)
    # if !getseq
    #     delete!(read,"sequence")
    # end
    # if !getCIGAR
    #     delete!(read,"CIGAR")
    # end
    if pairend
        read["pairno"]=Array{Int8}(read["readid"][:,1])
        read["readid"]=read["readid"][:,2:end]
    end
    if strand
        read["strand"]=map(x->x&16==0 ? '+' : '-',read["flag"])
        if pairend
            #Reverse second half-read in pair-end mode
            read["strand"][read["pairno"].==2]=map(x->x=='-' ? '+' : '-',read["strand"][read["pairno"].==2])
        end
    end
    recordno_ref[1]=RECORD_NO
    return read
end

#}}

#{{ parseBamMutation
#[[ parseBamMutation ]]
# T = parseBamMutation(Cigar, MD, read_seq; vcf_indel=false)
# Parse the mutation information according to the Cigar and MD fields.
# Output a table with fields: mut_typ, ref_loc, read_loc, ref_allele, alt_allele
# Where mut_typ is a charactor vector indicating the SNP ('S'), Insert ('I') or Delete ('D') mutations. ref_loc and read_loc are N x 2 matrix showing the mutation happened at the position on reference sequences (based on the left alignment point N0) and on read respectively. The absolute position of mutation can be calculated by N0+ref_loc-1.
# If vcf_indel is true, the left one nt adjacent to indel region will be included in the record.
#See also: readStarSamFile, readStarSamLine
#Xiong Jieyi, 7 Aug 2017

export parseBamMutation
function parseBamMutation(cigar::String, MD0::String, read::String; vcf_indel::Bool=false)
    rpl=rowpile()

    cgN,cgC=mapr(eachmatch(r"(\d+)([A-Z])",cigar)) do M
        (int(M.captures[1]), M.captures[2])
    end
    cg_readloc=zeros(Int,length(cgN),2)
    cg_refloc=zeros(Int,length(cgN),2)
    readIdxNoIns=collect(1:length(read)) #The index in MD is not include inserts and the "S" part of CIGAR
    rdp=1
    rfp=1
    for i=1:length(cgN)
        if cgC[i]=="S"
            cg_readloc[i,:]=[rdp,rdp+cgN[i]-1]
            cg_refloc[i,:]=[rfp,rfp-1]
            rdp+=cgN[i]
            deleteat!(readIdxNoIns, cg_readloc[i,1].<=readIdxNoIns.<=cg_readloc[i,2])
        elseif cgC[i]=="M"
            cg_readloc[i,:]=[rdp,rdp+cgN[i]-1]
            cg_refloc[i,:]=[rfp,rfp+cgN[i]-1]
            rdp+=cgN[i]
            rfp+=cgN[i]
        elseif cgC[i]=="I"
            cg_readloc[i,:]=[rdp,rdp+cgN[i]-1]
            cg_refloc[i,:]=[rfp,rfp-1]
            rdp+=cgN[i]
            addrow!(rpl, 'I', cg_refloc[i,:], cg_readloc[i,:], "", read[cg_readloc[i,:]|>x->x[1]:x[2]])
            deleteat!(readIdxNoIns, cg_readloc[i,1].<=readIdxNoIns.<=cg_readloc[i,2])
        elseif cgC[i]=="D"
            cg_refloc[i,:]=[rfp,rfp+cgN[i]-1]
            cg_readloc[i,:]=[rdp,rdp-1]
            rfp+=cgN[i]
        else
            error("Unknown GIGAR $cgC")
        end
    end
    @assert sum(seglen(cg_readloc))==length(read)
    
    MD=replace(MD0, r"(\d+[A-Z])0([A-Z])" => s"\1\2")
    mts=eachmatch(r"(\d+)(\^?[A-Z]+)",MD)
    if !isempty(mts)
        mdN,mdC=mapr(mts) do M
            (int(M.captures[1]), ascii(M.captures[2]))
        end
        cgi=1
        rdp_noins=1
        for i=1:length(mdN)
            rdp_noins+=mdN[i]
            rdp=readIdxNoIns[rdp_noins]
            
            if mdC[i][1]=='^'
                while cg_readloc[cgi,1]<rdp
                    cgi+=1
                end
                @assert cg_readloc[cgi,2]==rdp-1
                addrow!(rpl, 'D', cg_refloc[cgi,:], cg_readloc[cgi,:], mdC[i][2:end], "")
            else
                while cg_readloc[cgi,2]<rdp
                    cgi+=1
                end
                @assert cgC[cgi]=="M" && cg_readloc[cgi,1]<=rdp
                len=length(mdC[i])
                offset=cg_refloc[cgi,1]-cg_readloc[cgi,1]
                addrow!(rpl, 'S', [rdp rdp+len-1]+offset, [rdp rdp+len-1], mdC[i], read[rdp:rdp+len-1])
                rdp_noins+=len
            end
        end
        @assert int(match(r"\d+$",MD).match)+rdp_noins-1==length(readIdxNoIns)
    end
    
    T=tb(c"mut_typ, ref_loc, read_loc, ref_allele, alt_allele", value(rpl))

    if vcf_indel
        #Add the very left NT for INDEL
        for i=1:rownum(T)
            if T["mut_typ"][i]=='I' || T["mut_typ"][i]=='D'
                C=read[[T["read_loc"][i,1]-1]]
                T["ref_allele"][i]=C*T["ref_allele"][i]
                T["alt_allele"][i]=C*T["alt_allele"][i]
                T["ref_loc"][i,1]-=1
                T["read_loc"][i,1]-=1
            end
        end
    end
    
    T
end

#}}

#{{ parseCIGARpos
#[[ parseCIGARpos ]]
# pos_mx = parseCIGARpos(cigar::String, start_position::Int)
# Parse CIGAR code in bam/sam file, and output the position segments. If N exists in the CIGAR, the output will have multi-rows for each segment of junction read. This funciton is useful along with the BioAlignments.jl.
# See also:
# Xiong Jieyi, 17 Apr 2019

export parseCIGARpos
function parseCIGARpos(cigar::String, STR::Int)
    POS=zeros(Int,1,2)
    LOC=zeros(Int,1,2)
    p=1;rdp=1;rfp=STR;
    aL=findall(['A'<=x<='Z' for x in cigar])
    LOC_read_start=[];
    cjuncno=Int(0);
    outpos=rowmx(Int, 2)
    for cL=aL
        NO=int(cigar[p:cL-1]);p=cL+1;
        ccL=cigar[cL]
        if ccL=='S' || ccL=='H'
            rdp+=NO
            #S is not be calculated in the mapping start position.
            ## Remove this: ## rfp+=NO
            #Changed in Dec 14, 2015
        elseif ccL=='M'
            if !(POS[1]>0)
                LOC[1]=rdp
	        if isempty(LOC_read_start)
		    LOC_read_start=rdp
	        end
                POS[1]=rfp
            end
            rfp+=NO
            rdp+=NO
            LOC[2]=rdp-1
            POS[2]=rfp-1
        elseif ccL=='D'
            rfp+=NO;
        elseif ccL=='I'
            rdp+=NO
        elseif ccL=='N'
	    cjuncno+=1
            addrow!(outpos, POS)
            rfp+=NO
            POS=zeros(1,2)
        else
            error("Unknown cigar char in read: $ccL")
        end
    end
    addrow!(outpos, POS)
    value(outpos)
end
#}}

#{{ bamReadCount
#[[ bamReadCount ]]
# using BioAlignments #or using Bio.Align #Required.
# readnum, site_len = bamReadCount(indexed_bam_file, (SiteChr, SitePos[, SiteStrand]); revstrand=false, UMI="", rm_overlap_sites=true)
# readnum, site_len, site_grp = bamReadCount(...; sitegrp=Group, ...)
# readnum, site_len, site_index::Int, smp_barcode = bamReadCount(...; smpBC="BC_field", ...)
# readnum, site_len, site_grp, smp_barcode = bamReadCount(...; smpBC="BC_field",  sitegrp=Group, ...)
# readnum, site_len, site_grp|index, smp_barcode = bamReadCount(...; is10x=true, ...)
# Counting overlapped unique-mapped read number in each site or site group. A read overlapped with more than one site in a site group will only be counted once.
# `UMI` is the field of UMI barcode. e.g., for 10x data, set UMI="UB". When UMI is setted, for the reads with the same left positions and UMIs, only the first read was used.
# `smpBC` is the sample barcode. e.g. for 10x data, set smpBC="CB" to get the gene read number of every cell x every gene.
# is10x=true equals to UMI="UB", smpBC="CB".
# When rm_overlap_sites=true, the site segments overlapped by different site groups will be removed in advance. These parts will also not be counted for `site_len`.
# Xiong Jieyi, 27 Sep 2019

export bamReadCount
function bamReadCount(bamfn::AbstractString,
                      Site::Tuple{AbstractVector{String}, AbstractMatrix{Int}, AbstractVector{Char}};
                      revstrand::Bool=false, rm_overlap_sites::Bool=true,
                      UMI::String="", smpBC::String="", is10x::Bool=false, sitegrp::Union{Group, Nothing}=nothing)
    if is10x
        UMI="UB"
        smpBC="CB"
    end
    BAM=importpkg(:BioAlignments, preloaded=true).BAM::Module
    eachoverlap=importpkg(:BioAlignments, preloaded=true).eachoverlap::Function
    isBC=!isempty(smpBC)
    isstrand=length(Site)==3
    if isstrand
        Site::Tuple{AbstractVector{String}, AbstractMatrix{Int}, AbstractVector{Char}}
    else
        Site::Tuple{AbstractVector{String}, AbstractMatrix{Int}}
    end
    
    # Segment proprecess
    if isnothing(sitegrp)
        isSiteGrp=false
        sitegrp=1:rnum(Site)
        ositeidx=1:rnum(Site)
    else
        isSiteGrp=true
        Site, sitegrp=grpfunwith(sitegrp, Site) do csite
            pos=segunion(csite[2])[1]
            pos_len=size(pos, 1)
            if isstrand
                (fill(csite[1][1], pos_len), pos, fill(csite[3][1], pos_len))
            else
                (fill(csite[1][1], pos_len), pos)
            end
        end
    end
    if rm_overlap_sites
        segSt, t=segconc(Site)
        dupSt=getrow(segSt, t.>1)
        Site, t=segmentdiff(Site, dupSt, BnoOverlap=true)
        sitegrp=getrow(sitegrp, t)
    end

    if !isfile(bamfn*".bai")
        println("Bam file $bamfn haven't be indexed yet. Indexing....")
        sf"samtools index $1"(bamfn)
    end
    bamO=open(BAM.Reader, bamfn, index=bamfn*".bai")
    l=ismbr(Site[1], bamO.refseqnames)
    Site=getrow(Site, l)
    sitegrp=getrow(sitegrp, l)
    siteLenTb=tb(c"sitelen, sitegrp"=>grpfun(x->sum(seglen(x)), sitegrp, Site[2]))
    println("Site proprecessing has done.")
    
    siteGrpNum=rownum(siteLenTb)

    isUMI=!isempty(UMI)
    if isBC
        bcpool=Dict{String, Set{String}}()
        grpfun1=grpfunwith
    else
        rdpool=Set{String}()
        grpfun1=grpfun
    end
    # if isUMI
    #     umipool=Set{String}()
    #     prdpos=0
    # end
    C, stgrp=grpfunwith(sitegrp, Site, no=true, id=true) do sitegrpno::Int, sitegrpid, csite::Tuple
        @timetip println(f"$1 out of $2 ($(3=.3g)%) site groups are done."(sitegrpno, siteGrpNum, 100*sitegrpno/siteGrpNum))
        empty!(isBC ? bcpool : rdpool)
        cstrand=isstrand ? csite[3][1] : '.'
        in(cstrand, Set(['+', '-', '.'])) || error("Invalid strand '$cstrand'.")
        chr=csite[1][1]::String
        for pos in eachrow(csite[2])
            # empty!(umipool)
            for rd in eachoverlap(bamO, chr, pos[1]:pos[2])
                if cstrand=='.' || xor(BAM.ispositivestrand(rd)==(cstrand=='+'), revstrand) && rd["NH"]::UInt8==0x01 && (!isUMI || haskey(rd, UMI))
                    # if isUMI
                    #     if prdpos==(crdpos=BAM.position(rd);)
                    #         if in((cumi=rd[UMI]::String;), umipool)
                    #             continue
                    #         else
                    #             push!(umipool, cumi)
                    #         end
                    #     else
                    #         empty!(umipool)
                    #         prdpos=crdpos
                    #     end
                    # end
                    touched=true
                    if in('N', (cigar=BAM.cigar(rd);))
                        tagpos=parseCIGARpos(cigar, BAM.position(rd))
                        for i=1:size(tagpos, 2)-1
                            if tagpos[i, 2]<pos[1] && pos[2]<tagpos[i+1, 1]
                                touched=false
                                break
                            end
                        end
                    end
                    if touched
                        rdid=if isUMI
                            rd[UMI]::String
                        else
                            BAM.tempname(rd)
                        end
                        if isBC
                            if haskey(rd, smpBC)
                                cbc=rd[smpBC]
                                if haskey(bcpool, cbc)
                                    push!(bcpool[cbc], rdid)
                                else
                                    bcpool[cbc]=Set{String}([rdid])
                                end
                            end
                        else
                            push!(rdpool, rdid)
                        end
                    end
                end
            end
        end
        if isBC
            (collect(keys(bcpool)), map(length, values(bcpool)))
        else
            length(rdpool)
        end
    end
    close(bamO)
    println("All finished.")
    if isBC
        smp_bc, rdnum=C
    else
        rdnum=C
    end
    if isBC
        rdnum, dtshift(stgrp, siteLenTb["sitegrp"], siteLenTb["sitelen"]), stgrp, smp_bc
    elseif isSiteGrp
        rdnum, dtshift(stgrp, siteLenTb["sitegrp"], siteLenTb["sitelen"]), stgrp
    else
        dtshift(ositeidx, stgrp, rdnum, sitelen, 0), dtshift(ositeidx, siteLenTb["sitegrp"], siteLenTb["sitelen"], 0)
    end
end
#}}

#{{ parseGTF
#[[ parseGTF ]]
# exon_table = parseGTF( GTF_filename[, additional_field1, additional_field2, ... ])
# Parse GTF file. Output will contain below fields:
#       chrno, exn_pos, cds_pos, strand,
#       geneid="",
#       transid="",
#       exn_rank=0,
#       gene_name="",
#       gene_typ="",
#       trans_name="",
#       protid="",
#       exnid=""
#If additional fields are needed, the 2-N parameters should be (OutFieldName::AbstractString, GTFFieldName::AbstractString, Function, DefaultValue::Any) format. Only "exon" rows in GTF files are considered. For the rows missed this item, default values will be filled.
#See also: readStarSamFile, genome_fasta2jld, parseDemulReport
#Xiong Jieyi, May 21, 2015 >Jul 30, 2015 >Aug 8, 2015

export parseGTF
function parseGTF(fn, fields::Tuple{AbstractString,AbstractString,Function,Any}...)
    function phasefd(X)
        m=match(r"^\s*(\S+)\s+\"([^\"]*)\"\s*$",X)
        if m==nothing
            m=match(r"^\s*(\S+)\s+(\S*)\s*$",X)
        end
        if m==nothing
            @warn("Invalid item: $X")
            ("","")
        else
            (m.captures[1],m.captures[2])
        end
    end
    fielddict=Dict{AbstractString, Tuple{AbstractString,Function}}("gene_id"=>("geneid",ascii),
               "transcript_id"=>("transid",ascii),
               "exon_number"=>("exn_rank",int),
               "gene_name"=>("gene_name",ascii),
               "gene_biotype"=>("gene_typ",ascii),
               "gene_type"=>("gene_typ",ascii),
               "transcript_name"=>("trans_name",ascii),
               "protein_id"=>("protid",ascii),
               "exon_id"=>("exnid",ascii))
    empexn=tb(geneid="",
              transid="",
              exn_rank=0,
              gene_name="",
              gene_typ="",
              trans_name="",
              protid="",
              exnid="")
    for c in fields
        # @assert(isa(c,(AbstractString,AbstractString,Function,Any)),"Invalid additional fields. The input should be (OutFieldName::AbstractString, GTFFieldName::AbstractString, Function, DefaultValue::Any).")
        fielddict[c[2]]=(c[1],c[3])
        empexn[c[1]]=c[4]
    end
    exn=rowpile()
    cds=rowpile()
    fileloop(fn, ignore='#')do li
        C=split(rstrip(li),'\t')
        if C[3]!="exon" && C[3]!="CDS"
            return nothing
        end
        chrno=chr2no(C[1])
        if chrno==0
            return nothing
        end
        if length(C)>9
            Ced=join(C[9:end]," ")
        else
            Ced=C[9]
        end
        D=split(rstrip(Ced,';'),';')
        cexn=copy(empexn)
        d"chrno,exn_pos,strand"cexn=(chrno,int(C[4:5])',C[7][1])
        for cfd in D
            k,v=phasefd(cfd)
            if isempty(k)
                return nothing #When parsed an invalid item, give up this record.
            end
            if haskey(fielddict,k)
                kk,fun=fielddict[k]
                cexn[kk]=fun(v)
            end
        end
        if C[3]=="exon"
            addrow!(exn,cexn)
        else
            # C[3]=="CDS"
            addrow!(cds,cexn)
        end
     end
    exn=value(exn)
    cds=value(cds)
    if cds!=nothing
        (ie,ic)=genomemap(d"transid,exn_pos"exn,d"transid,exn_pos"cds)
        exn["cds_pos"]=zeros(Int,rownum(exn),2)
        exn["cds_pos"][ie,:]=cds["exn_pos"][ic,:]
        exn["protid"][ie,:]=cds["protid"][ic,:]
    end
    return exn
end
#}}

#{{ writeGTF
#[[ writeGTF ]]
#writeGTF(filename, chr|chrno, pos, geneid; strand='.', prefix="chr")
#Write a simplest GTF file. This GTF file can be used by HTSeq-count. prefix is used in no2chr.
#See also:no2chr
#Xiong Jieyi, May 15, 2016

export writeGTF
#<v0.6# function writeGTF{T<:Real}(filename::AbstractString,chr::Vector,pos::Matrix{T},geneid::Vector;strand::Union{Void,Vector{Char}}=nothing,prefix="chr")
function writeGTF(filename::AbstractString,chr::Vector,pos::Matrix{T},geneid::Vector;strand::Union{Void,Vector{Char}}=nothing,prefix="chr") where {T<:Real}
    if eltype(chr)<:Real
        chr=no2chr(chr,prefix=prefix)
    end
    if strand==nothing
        strand=fill('.',length(chr))
    end
    open(filename,"w") do fid
        rowloop(chr,pos,geneid,strand) do cchr,cpos,cgeid,cstrand
            println(fid,"$cchr\t.\texon\t$(cpos[1])\t$(cpos[2])\t.\t$cstrand\t.\tgene_id \"$cgeid\";")
        end
    end
    nothing
end
#}}

#{{ segunion

#[[ segunion ]]
# (seg_pos, score)=segunion(pos; score=1s, weld=true)
# ((seg_chr, seg_pos), score)=segunion(chr, pos;...)
# ((seg_chr, seg_pos, seg_strand), score)=segunion(chr, pos, strand;...)|segunion((chr,pos,strand);...)
#Union the segments. The output score is the sum of inputed scores for each union segment.
#If weld is true, [1 3] will be merged with [4 6].
#See also: segconc, seg2cvg
#Xiong Jieyi, 1 Sep 2014 >October 21, 2014 >2 Mar 2017>14 May 2018

function segunion(pos::Matrix{T};score::AbstractVector{T2}=ones(Int,size(pos,1)), weld::Bool=true) where {T<:Real,T2<:Real}
    if isempty(pos)
        return zeros(T,0,2)
    end
    len=size(pos,1)
    W=weld ? 1 : 0
    t,idx=sortr([pos[:,1]-W ones(T,len); pos[:,2] 2*ones(T,len)])
    st=t[:,1]
    yy=t[:,2]
    yy[yy.==2].=-1
    dscore=[score;zeros(typein(score),length(score))][idx]
    
    cur=0
    count=0
    bd=st[1]+W
    rlt=rowpile()
    for i=1:length(st)
        cur+=yy[i]
        if cur==0
            addrow!(rlt,[bd, st[i]],count)
            count=0
            maxcover=0
            if i<length(st)
                bd=st[i+1]+W
            end
        elseif yy[i]==1
            count+=dscore[i]                
        end
    end
    return value(rlt)
end
#<v0.6# function segunion{T<:Real}(chr::Group,pos::Array{T};args...)
function segunion(chr::Group,pos::Array{T};args...) where {T<:Real}
    (seg_pos,seg_score),seg_chr=grpfunwith(x->segunion(x;args...),chr,pos)
    return ((seg_chr,seg_pos),seg_score)
    # O,t=grpfun(x->segunion(x;args...),chr,pos)
    # seg_chr,t=vcatr_with(t,O...)
    # return ((seg_chr,t[1]),t[2])
end
#<v0.6# function segunion{T<:Real}(chr::Array,pos::Array{T},strand::Array;args...)
function segunion(chr::Array,pos::Array{T},strand::Array;args...) where {T<:Real}
    (seg_pos,seg_score),(seg_chr,seg_strand)=grpfunwith(x->segunion(x;args...),(chr,strand),pos)
    return ((seg_chr,seg_pos,seg_strand),seg_score)
    # O,t=grpfun(x->segunion(x;args...),(chr,strand),pos)
    # tt,t=vcatr_with(t,O...)
    # return ((tt[1],t[1],tt[2]),t[2])
end
segunion(X::Tuple;args...)=segunion(X...;args...)
export segunion

#}}

#{{ findintron

#[[ findintron ]]
#intron_pos = findintron(exon_pos)
#(transid, intron_pos)=findintron(transid,exon_pos)
#Looking for intron using given exon position.
#See also: segunion, segconc, findorf
#Xiong Jieyi, 16 Sep 2014

#<v0.6# function findintron{T<:Real}(pos::Matrix{T})
function findintron(pos::Matrix{T}) where {T<:Real}
    allpos=[minimum(pos[:,1]) maximum(pos[:,2])]
    segconc(allpos,1,pos,-1)[1]
end
#<v0.6# function findintron{T<:Real}(transid::Group,pos::Matrix{T})
function findintron(transid::Group,pos::Matrix{T}) where {T<:Real}
    grpfunwith(x->findintron(x),transid,pos)[[2,1]]
end
export findintron

#}}

#{{ getseq

#[[ getseq ]]
#Seqs = getseq((chr, pos[, strand]);path="path_of_sequence_lib",fun=x->x)
#Or Seqs = getseq(fun, ...)
#Seqs = getseq(pos[, strand];path="path_of_sequence_lib",fun=x->x)
#Fetch the sequences from reference genomes. If the strand is given, the '-' sequence will be reverse-complement transferred. You can build the path_of_sequence_lib by genome_fasta2jld.
#See also: catseq, segunion, genome_fasta2jld
#Xiong Jieyi, 5 Sep 2014>February 26, 2015>Jan 1, 2016

# function needBioSeq()
#     isdefined(:PyCall)||require("PyCall")
#     if !isdefined(:BioSeq)
#         @eval Main.PyCall.@pyimport Bio.Seq as BioSeq
#         @eval Main.PyCall.@pyimport Bio.SeqIO as BioSeqIO
#     end
# end
function getseq(chrpos::Tuple{Group,Array};
                path::ASCIIString="",fun::Function=x->typeof(x)[x])
    chr,pos=chrpos
    grpfunexp(chr,pos;input_grpid=true) do cchr, cpos
        if isa(cchr,Real)
            cchr=no2chr(cchr)
        end
        file=joinpath(path,cchr)
        getseq(cpos,file;fun=fun)
    end
end
function getseq(chrposstrand::Tuple{Group,Array,Vector{Char}};
                path::ASCIIString="",fun::Function=x->typeof(x)[x])
    chr,pos,strand=chrposstrand
    grpfunexp((chr,strand),pos;input_grpid=true) do cchr_strand, cpos
        cchr,cstrand=cchr_strand
        if isa(cchr,Real)
            cchr=no2chr(cchr)
        end
        file=joinpath(path,cchr)
        getseq(cpos,file;strand=cstrand,fun=fun)
    end
end
function getseq(pos::Array{T},file::AbstractString; strand::Char='+',fun::Function=x->typeof(x)[x]) where {T<:Real}
    # needBioSeq()
    BioSeq=importpy("Bio.Seq")
    seq=jldload(file, "seq") #Data still saved in JLD format, since JLD2 do not support long string. 15 Jun 2019
    rlt=rowpile()
    for i=1:size(pos,1)
        csq=seq[pos[i,1]:pos[i,2]]
        if strand=='-'
            csq=BioSeq.reverse_complement(csq)
        end
        addrow!(rlt,fun(csq))
    end
    return value(rlt)
end
getseq(fun::Function,args...;wargs...)=getseq(args...;fun=fun,wargs...)
# getseq(X::Tuple;args...)=getseq(X...;args...)
export getseq

#}}

#{{ genome_fasta2jld
#[[ genome_fasta2jld ]]
# genome_fasta2jld(fasta_filename, output_path)
#Convert genome fasta file to binrary sequence library. This path can be used by getseq.
#The reason this function still saves data using JLD.jl rather than JLD2.jl is JLD2.jl do not support storage long string so far (15 Jun 2019).
#See also: getseq, parseGTF, readStarSamFile
#Xiong Jieyi,February 26, 2015 > 14 May 2018

export genome_fasta2jld
function genome_fasta2jld(fn::AbstractString,pth::AbstractString)
    if !isdir(pth)
        println("Create path $pth")
        mkdir(pth)
    end
    open(fn)do fp
        chr=""
        iob=IOBuffer()
        while !eof(fp)
            ln=strip(chomp(readline(fp)))
            if ln[1]=='>'
                if !isempty(chr)
                    # seq=takebuf_string(iob)
                    seq=String(take!(iob))
                    jldsave(joinpath(pth,chr), "seq"=>seq, "chr"=>chr, "length"=>length(seq))
                end
                chr=ascii(split(strip(ln[2:end]))[1])
                iob=IOBuffer()
            else
                print(iob,ln)
            end
        end
        # seq=takebuf_string(iob)
        seq=String(take!(iob))
        jldsave(joinpath(pth,chr), "seq"=>seq, "chr"=>chr, "length"=>length(seq))
    end
end
#}}

#{{ catseq

#[[ catseq ]]
#(trans_seq, exn_loc_in_trans)=catseq(exn_seqs,exn_pos,strand='+')
#((trans_seq, transid), exn_loc_in_trans)=catseq(exn_seq,(transid,pos,strand))
#Caternate exon sequence together, and also return the location of each exon in the transcrition.
#exn_loc_in_trans is a N_exn x 2 matrix with the corresponding line as input.
#Supported sequence could be AbstractString or any vector.
#The coverlaps of exn_pos will be checked automatically, and a error will occur if happened.
#See also: getseq, segmask
#Xiong Jieyi, 4 Oct 2014>May 28, 2015>Aug 19, 2015

function catseq(seq::Vector{T},pos::Matrix,strand::Char='+') where {T}
    if size(pos,2)>1
        if strand=='+'
            idx=sortri(pos)
            spos=pos[idx,:]
            @assert(all(spos[1:end-1,2].<spos[2:end,1]),"Segments have overlap.")
        elseif strand=='-'
            idx=sortri(pos,rev=true)
            spos=pos[idx,:]
            @assert(all(spos[2:end,2].<spos[1:end-1,1]),"Segments have overlap.")
        else
            error("Strand could neither be + or -.")
        end
    else
        spos=pos
    end
    sseq=seq[idx]
    # spos=pos[idx,:]
    if T<:AbstractString
        trseq=*(sseq... )
    else
        trseq=[sseq...]
    end
    t=cumsum(spos[:,2]-spos[:,1]+1)
    trloc=[[1;t[1:end-1]+1] t][rvorder(idx),:]
    return (trseq,trloc)
end
function catseq(seq::Vector{T},trid_pos_strand::Tuple{Group,Matrix,Vector{Char}}) where {T}
    (transid,pos,strand)=trid_pos_strand
    trT=rowpile()
    trans_loc=grpfunexp((transid,strand),seq,pos;input_grpid=true) do trid_strand, cseq, cpos
        (trseq,trloc)=catseq(cseq,cpos,trid_strand[2])
        if T<:AbstractString
            addrow!(trT,trseq,trid_strand[1])
        else
            addrow!(trT,Any[trseq],trid_strand[1])
        end
        trloc
    end
    return (value(trT),trans_loc)
end
export catseq

#}}

#{{ findorf
#[[ findorf ]]
#pos_mx = findorf(Seq)
#Greedy search all open reading frame (ATG-TAA|TGA|TAG) in the sequence.
#See also: getseq, segmask
#Xiong Jieyi, 7 Sep 2014>7 Oct 2014

function findorf(Seq::ASCIIString)
    stri=1
    pos=rowmx()
    while (m=match(r"ATG(?:[ATGC]{3})*(?:TAA|TGA|TAG)"i,Seq,stri))!=nothing
        addrow!(pos,[m.offset,m.offset+length(m.match)-1])
        stri=m.offset+1
    end
    return value(pos)
end
export findorf
#}}

#{{ segmask
#[[ segmask ]]
#(mask_pos, mask_loc)=segmask(seg_pos, mask_start, mask_len; strand='+')
#(mask_pos, mask_loc)=segmask(seg_pos, [mask_start mask_end]; strand='+')
#Mask the given segments, and output the masked absolute-position and relative-location. 0 in output means this segment is not in the mask region. When strand='-', segment will be masked from right. The inputed seg_pos isn't needed to be sorted, and the output is always consistent with inputed seg_pos rowly-by-rowly.
#e.g., cds_pos=segmask(exn_pos,orf_loc_in_trans; strand=strand)
#See also: findorf, getseq, catseq
#Xiong Jieyi, 5 Oct 2014

#<v0.6# function segmask{T}(pos::Matrix{T},mask_str::Real,mask_len::Real;strand::Char='+')
function segmask(pos::Matrix{T},mask_str::Real,mask_len::Real;strand::Char='+') where {T}
    if strand=='-'
        (mask_pos, mask_loc)=segmask(-fliplr(pos),mask_str,mask_len;strand='+')
        return (-fliplr(mask_pos), mask_loc)
    end

    idx=sortri(pos)
    pos=pos[idx,:]
    
    mask_skp=mask_str-1
    mask_pos=zeros(T,size(pos))
    mask_loc=zeros(T,size(pos))
    seg_len=seglen(pos)
    p=1
    while p<=size(pos,1) && seg_len[p]<=mask_skp
        mask_skp-=seg_len[p]
        p+=1
    end
    next_mask_loc=1
    while p<=size(pos,1)
        mask_loc[p,1]=next_mask_loc
        mask_pos[p,1]=pos[p,1]+mask_skp
        if mask_len<=seg_len[p]-mask_skp
            mask_pos[p,2]=mask_pos[p,1]+mask_len-1
            mask_loc[p,2]=mask_loc[p,1]+mask_len-1
            break
        else
            mask_pos[p,2]=pos[p,2]
            mask_loc[p,2]=mask_loc[p,1]+seg_len[p]-mask_skp-1
            mask_len-=seg_len[p]-mask_skp
            next_mask_loc=mask_loc[p,2]+1
        end
        p+=1
        mask_skp=0
    end
    vidx=rvorder(idx)
    (mask_pos[vidx,:], mask_loc[vidx,:])
end
segmask(pos::Matrix,mask_pos::Array;wargs...)=segmask(pos,mask_pos[1],mask_pos[2]-mask_pos[1]+1;wargs...)
export segmask
#}}

#{{ segconc

#[[ segconc ]]
# (pos,score)=segconc(pos,score=1;cutoff=1)
# ((chr,pos),score)=segconc((chr,pos),score=1;cutoff=1)
# ((chr,pos,strand),score)=segconc((chr,pos,strand),score=1;cutoff=1)
# ...=segconc(A,score_A,B,score_B;cutoff)
#Sort, merge segments and tide scores. The result bed are sorted and unoverlapped.
#See also: segunion, seg2cvg, segintersect segmentdiff, segmask
#Xiong Jieyi, 8 Sep 2014>22 Sep 2014

function segconc(pos::Matrix{T1},
               score::Union{AbstractVector,Real}=1;cutoff::Real=1) where {T1<:Real}
    if isa(score,Real)
        score=fill(score,size(pos,1))
    end
    P0,t=sortr([pos[:,1];pos[:,2]+1])
    tt=[score;-score]
    S0=tt[t]
    ldel=fill(false,length(S0))
    if S0[1]==0 && cutoff>0
        ldel[1]=true
    end
    for i=2:length(S0)
        if P0[i]==P0[i-1]
            ldel[i-1]=true
            S0[i]+=S0[i-1]
        end
        if S0[i]==0 && (cutoff>0 || i<length(S0))
            ldel[i]=true
        end
    end
    P=P0[.!ldel]
    S=S0[.!ldel]

    rlt=rowpile()
    if isempty(S)
        return (zeros(typein(pos),0,2),ones(typein(score),0))
    else
        V=S[1]
        for i=2:length(S)
            if V>=cutoff
                addrow!(rlt,[P[i-1] P[i]-1],V)
            end
            V+=S[i]
        end
        ret=value(rlt)
        if ret==nothing
            return (zeros(typein(pos),0,2),ones(typein(score),0))
        else
            return ret
        end
    end
end
function segconc(chrpos::Tuple{Group,Matrix,AbstractVector},
                 score::Union{AbstractVector,Real}=1;cutoff::Real=1,args...)
    (chr,pos,strand)=chrpos
    if isa(score,Real)
        score=fill(score,size(pos,1))
    end
    # (pos_score,chr_strand)=grpfun((chr,strand),pos,score) do cpos,cscore
    #     segconc(cpos,cscore;args...)
    # end
    # ((ochr,ostrand),(opos,oscore))=vcatr_with(chr_strand,pos_score...)
    (opos,oscore),(ochr,ostrand)=grpfunwith((chr,strand),pos,score) do cpos,cscore
        segconc(cpos,cscore;args...)
    end
    return ((ochr,opos,ostrand),oscore)
end
function segconc(chrpos::Tuple{Group,Matrix},score::Union{AbstractVector,Real}=1;args...)
    (chr,pos)=chrpos
    if isa(score,Real)
        score=fill(score,size(pos,1))
    end
    # (pos_score,gchr)=grpfun(chr,pos,score) do cpos,cscore
    #     segconc(cpos,cscore;args...)
    # end
    # (gchr,(opos,oscore))=vcatr_with(gchr,pos_score...)
    (opos,oscore),ochr=grpfunwith(chr,pos,score) do cpos,cscore
        segconc(cpos,cscore;args...)
    end
    return ((ochr,opos),oscore)
end
function segconc(Ap::Union{Matrix,Tuple}, As::Union{AbstractVector,Real},
                 Bp::Union{Matrix,Tuple}, Bs::Union{AbstractVector,Real};args...)
    if !isa(As,AbstractVector)
        As=fill(As,rownum(Ap))
    end
    if !isa(Bs,AbstractVector)
        Bs=fill(Bs,rownum(Bp))
    end
    segconc(vcatr(Ap,Bp),vcatr(As,Bs);args...)
end
export segconc

#}}

#{{ segintersect segmentdiff
#[[ segintersect segmentdiff ]]
#(pos_mx|(chrno,pos[,strand]), Idx_segA)=segintersect|segmentdiff(segA,segB; BnoOverlap=false)
# pos_mx, Idx_segA = segintersect(A_pos_mx, B_pos_mx; BnoOverlap=false, Asorted=false, Bsorted=false, sortcheck=true)
#Get the seg parts of segA in/not in region of segB. Input could be pos_matrix or (chr,pos[, strand]) tuple.
#If BnoOverlap is ture, segintersect will not union B first. Note that when B have overlap, the results are wrong. At the moment, when BnoOverlap=false, the segunion() step will not be optimized even if Bsorted=true (5 Sep 2019).
#See also: segconc, segunion, seg2cvg, segmask
#Xiong Jieyi, February 8, 2015 > 5 Sep 2019

export segintersect, segmentdiff
function segintersect(segA::Tuple,segB::Tuple; BnoOverlap::Bool=false)
    if !BnoOverlap && rownum(segB)>1
        segB,=segunion(segB)
    end
    difseg,Ai=genomemap(segA,segB,touch=true,
              fun=(ai,bi,ap,bp)->([max(ap[1],bp[1]) min(ap[2],bp[2])], ai))
    if length(segA)==2
        difseg=tuple(getrow(segA[1],Ai),difseg)
    else
        @assert(length(segA)==3)
        difseg=tuple(getrow(segA[1],Ai),difseg,getrow(segA[3],Ai))
    end
    return (difseg, Ai)
end
function segintersect(segA::Matrix{T1},segB::Matrix{T2}; BnoOverlap::Bool=false, Asorted::Bool=false, Bsorted::Bool=false, sortcheck::Bool=true) where {T1<:Real,T2<:Real}
    if !BnoOverlap && size(segB,1)>1
        segB,=segunion(segB)
    end
    coordinatemap(segA,segB,touch=true, no_match_output=(zeros(T1, 0, 2), Int[]),
                  fun=(ai,bi,ap,bp)->([max(ap[1],bp[1]) min(ap[2],bp[2])], ai),
                  sitesorted=Asorted, tagsorted=Bsorted, sortcheck=sortcheck)
end
function segmentdiff(segA,segB; BnoOverlap::Bool=false)
    if !BnoOverlap && size(segB,1)>1
        segB,=segunion(segB)
    end
    segD,=segconc(segunion(segA)[1],1,segB,-1)
    segintersect(segA,segD;BnoOverlap=true)
end
#}}

#{{ seg2cvg
#[[ seg2cvg ]]
#(coverage::Vector{Int}, X::Range) = seg2cvg(pos[, start, end]; max_pixel=0)
#coverage::Vector{Int} = seg2cvg(pos, x_range::Range)
#Calculate coverage using segments.
#If max_pixel >0 and end-start+1>max_pixel, coverage will be conpressed, and the length of output will be between 1/2 max_pixel and max_pixel.
#See also: segconc, segunion
#Xiong Jieyi, 10 Sep 2014>December 12, 2014>Mar 3, 2016

function seg2cvg(pos::Array{T},
                 limsta::Real=minimum(pos[:,1]),limend::Real=maximum(pos[:,2]); max_pixel::Integer=0) where {T<:Real}
    len=limend-limsta+1
    if 0<max_pixel<len
        N=ceil(Int,len/max_pixel)
    else
        N=1
    end
    X=limsta:N:limend
    cvg=seg2cvg(pos,X)
    return (cvg,X)
end
function seg2cvg(pos::Array{T}, xrang::Range) where {T<:Real}
    len=length(xrang)
    cvg=zeros(Int,len)
    rpos=round.(Int,(pos-first(xrang))/step(xrang)+1)
    for i=1:size(rpos,1)
        cvg[max(rpos[i,1],1):min(rpos[i,2],len)]+=1
    end
    return cvg
end
export seg2cvg
#}}

#{{ writefasta
#[[ writefasta ]]
#writefasta("filename",seqs,heads;append=false)
#Write sequence to fasta file. head could be string vector or numerical vector.
#See also: getseq, readfasta
#Xiong Jieyi, 16 Sep 2014

#<v0.6# function writefasta{T1<:ASCIIString,T2<:AbstractString}(filename::AbstractString,seq::Vector{T1},head::Vector{T2};append::Bool=false)
function writefasta(filename::AbstractString,seq::Vector{T1},head::Vector{T2};append::Bool=false) where {T1<:ASCIIString,T2<:AbstractString}
    open(filename,append ? "w+" : "w") do fp
        rowloop(seq,head) do cseq,chead
            write(fp,">"*chead*"\n")
            write(fp,cseq*"\n\n")
        end
    end
end
writefasta(filename::AbstractString,seq::Vector{T1},head::AbstractVector{T2};args...) where {T1<:ASCIIString,T2<:Real} =
    writefasta(filename,seq,num2str(head);args...)
export writefasta
#}}

#{{ readfasta
#[[ readfasta ]]
#(seqeuence, head)=readfasta(filename|pipe_cmd; skipline::Int=0, ignore::Char='\0')
#Read fasta file. The line breaks and side blanks in the sequence will be cleaned.
#See also: writefasta, getseq
#Xiong Jieyi, 16 Sep 2014

# function readfasta(filename)
#     # needBioSeq()
#     importpkg(:PyCall, preloaded=true)
#     BioSeqIO=importpy("Bio.SeqIO")
#     oo=BioSeqIO.parse(filename,"fasta")
#     str(x)=Main.PyCall.pycall(Main.PyCall.pybuiltin("str"),Main.PyCall.PyAny,x)
#     t=map(oo) do x
#         (str(x.seq), x.id)
#     end
#     return vcatr(t...)
# end
export readfasta
function readfasta(fn; skipline::Int=0, ignore::Char='\0')
    head=String[]
    seq=String[]
    buf=IOBuffer()
    fileloop(fn, no=true, skipline=skipline, ignore=ignore) do no, li
        isempty(strip(li)) && return nothing
        if li[1]=='>'
            no>1 && push!(seq, String(take!(buf)))
            push!(head, strip(li[2:end]))
        else
            write(buf, strip(li))
        end
    end
    push!(seq, String(take!(buf)))
    (seq, head)
end
#}}

#{{ balancemask balancemask2
#[[ balancemask balancemask2 ]]
#lA, lB = balancemask(A, B; N=10)
#lB, times_B/A = balancemask2(A_small, B_big; N=1~10)
#Divide [A,B] into N group, and mask in order to keep the number of A and B balance in the two groups. A[lA] and B[lB] is balanced. (lA and lB are BitVector).
#When using balancemask2, #B could no more than #A in any group (A error will occur if this assumption is failed). Function will try to keep as many as B while keep A and B in balance.
#see also: drawer, GOenrichFB, GOenrichGrp
#Xiong Jieyi, October 19, 2014>March 9, 2015>Nov 12, 2015>8 Sep 2018

function balancemask(A::Vector,B::Vector;N::Int=0)
    Idx=sortri([A;B])
    if N==0
        N=min(max(floor(Int,length(Idx)./10),1),10)
        println("Group by $N bins.")
    end
    sL=[trues(length(A)); falses(length(B))][Idx]
    oL=trues(length(sL))
    grpelnum=drawer(length(sL),N)
    p=1
    for gi=1:N
        l=p:p+grpelnum[gi]-1
        cL=sL[l]
        nA=sum(cL)
        nB=length(cL)-nA
        if nA>nB
            oL[l[find(cL)[randperm(nA)[1:(nA-nB)]]]].=false
        elseif nA<nB
            oL[l[find(.!cL)[randperm(nB)[1:(nB-nA)]]]].=false
        end
        p=p+grpelnum[gi]
    end
    ooL=oL[rvorder(Idx)]
    return ooL[1:length(A)],ooL[length(A)+1:end]
end
function balancemask2(A::Vector,B::Vector;N::Int=0)
    Idx=sortri([A;B])
    if N==0
        N=min(max(floor(Int,length(Idx)./10),1),10)
        println("Group by $N bins.")
    end
    sL=[trues(length(A)); falses(length(B))][Idx]
    oL=trues(length(sL))
    grpelnum=drawer(length(sL),N)
    p=1
    minratio=Inf
    for gi=1:N
        l=p:p+grpelnum[gi]-1
        cL=sL[l]
        nA=sum(cL)
        nB=length(cL)-nA
        cratio=nB/nA
        if cratio<minratio
            minratio=cratio
        end
        p=p+grpelnum[gi]
    end
    if minratio<1
        error("B in some group is even less than A. Try to decrease N or using balancemask instead.")
    end
    p=1
    for gi=1:N
        l=p:p+grpelnum[gi]-1
        cL=sL[l]
        nA=sum(cL)
        nB=length(cL)-nA
        oL[l[find(.!cL)[randperm(nB)[1:(nB-round(Int,nA*minratio))]]]].=false
        p=p+grpelnum[gi]
    end
    ooL=oL[rvorder(Idx)]
    return ooL[length(A)+1:end],minratio
end
export balancemask, balancemask2
#}}

#{{ GOenrichFB GOenrichGrp
#[[ GOenrichFB GOenrichGrp ]]
#Table|nothing = GOenrichFB(foreground, background; balance=(fg,bg), ann="org.Hs.eg.db", ID="ensembl", outfile="filename.tsv", sigP=0.05|0.01, topGO_algorithm=:classic, padj_method=:fdr|:none)
#Table|nothing, GOdata, GOtest = GOenrichFB(...; outRobj=true) #Output topR objects also.
#Table|nothing = GOenrichGrp(gene, group; balance=(fg,[ bg=fg]), bg=gene, ann="org.Hs.eg.db", ID="ensembl", outfile="filename.tsv", sigFDR=0.05)
#Check the GO enrichment using Fisher test in topGO. If outfile is given, a text file will be saved. If any GO terms are enriched, a table will be return, otherwise will return nothing instead.
#If the background is need to be balanced, set balance value( e.g. as gene expression). Function will use balancemask2, i.e., only subsample background but not subsample forground. In the case that background is not enough, an error is occur. The shared item in background, and the duplicated item will be remove, so does their balance values, before apply to balancemask function.
#The supported ID could be "entrez", "genbank", "alias", "ensembl", "symbol", "genename".
#padj_method: the method of p.adjust in R. Can be :holm, :hochberg, :hommel, :bonferroni, :BH, :BY, :fdr, :none. The default is :fdr if method is :classic, or :none otherwise.
#sigP: The default is 0.05 if p-value adjustment is used, or 0.01 otherwise.
#topGO_algorithm: can be :classic|:weight01|...
#* (i) 'elim' : this method processes the GO terms by traversing the GO hierarchy from bottom to top, ie. it first assesses the most specific (bottom-most) GO terms, and proceeds later to more general (higher) GO terms. When it assesses a higher (more general) GO term, it discards any genes that are annotated with significantly enriched descendant GO terms (considered significant using a pre-defined P-value threshold). This method does tend to miss some true positives at higher (more general) levels of the GO hierarchy.
#* (ii) 'weight': this is designed to miss less true positives than the 'elim' method. The significance scores of connected nodes (a parent and its child in the GO hierarchy) are compared in order to detect the locally most significant terms in the GO hierarchy. 
#* (iii) 'classic': each GO term is tested independently, not taking the GO hierarchy into account.
#* (iv) 'weight01': this is the default method used by TopGO, and is a mixture of the 'elim' and 'weight' methods.
#* (v) 'parentchild': described by Grossman et al (2007) : when assessing a GO term, it takes into accoount the annotation of terms to the current term's parents, and so reduces false positives due to the inheritance problem. You can also use this algorithm in a very nice tool called Ontologizer , which has a lovely GUI (I used it in a collaboration on gene duplicates with Shane Woods, described in his paper Woods et al (2013)). 
# The 'elim' and 'weight' methods are designed to be more conservative (have less false positives) than the 'classic' method, so should have larger p-values in general. Simulation results reported by Alexa et al (2006) show that the 'weight' algorithm has less false positives than 'classic' and misses few true positives, while 'elim' has even less false positives than 'weight' but misses more true positives.  Grossman et al (2007) compare the 'elim', 'weight' and 'parentchild' algorithms, and say that each has advantages in certain scenarios. 
#See also: balancemask, balancemask2, GSVA
#Xiong Jieyi, October 27, 2014>Aug 28, 2015>Oct 2, 2015>23 Apr 2017>16 Oct 2019

export GOenrichFB,GOenrichGrp
function _GOenrichCore(val::BitVector,genes::AbstractVector{T},ann::AbstractString,ID::AbstractString; topGO_algorithm::Symbol=:classic, padj_method::Symbol=topGO_algorithm==:classic ? (:fdr) : (:none), sigP::Real=padj_method==:none ? 0.01 : 0.05) where {T<:AbstractString}
    myr=importpkg(:myR)
    # DF=importpkg(:DataFrames)
    @pkgfun(myR, eR, callr, relem, j2r, r2j, callrj, callr1)
    # @pkgfun(DataFrames,DataFrame,vcat=>DF_vcat,hcat=>DF_hcat)
    # rlt=DataFrame()
    rlt=rowpile()
    GOdatas=Dict{String,myr.RObject}()
    GOtests=Dict{String,myr.RObject}()
    ontology_ls=c"BP,MF,CC"
    for ontology in ontology_ls
        println("*** TEST $ontology ONTOLOGY ***")
        GOdata=callr("new","topGOdata",
                     ontology=ontology,
                     allGenes=j2r(int(val),names=genes),
                     geneSel=eR("function(x){x==1}"),
                     annot=eR("annFUN.org"), mapping=ann, ID=ID,
                     nodeSize=5)
        GOdatas[ontology]=GOdata
        GOtest=callr("runTest",GOdata,algorithm=topGO_algorithm,statistic="fisher")
        GOtests[ontology]=GOtest
        allGOlen=callr("length", callr("usedGO",object = GOdata))
        GenTb=callr("GenTable",GOdata,classicFisher = GOtest, topNodes = allGOlen, numChar=9999)
        FDR=callr("p.adjust",relem(GenTb,"classicFisher"),method=padj_method)

        #For unknown reason, classicFisher and Term occationally become invalid StringVector after transferred to Julia DataFrame.
        GenTb=callr("function(y,z,w){y['classicFisher']<-as.vector(data.matrix(y['classicFisher']));subset(cbind(y, Padj_$padj_method=z), Padj_$padj_method<=w)}",GenTb,FDR,sigP)
        # end

        N=callrj("nrow", GenTb)
        if N>0
            println("$N GO item(s) is/are enriched.")
            # rlt=DF_vcat(rlt,DF_hcat(DataFrame(Ontology=fill(ontology,N)), r2j(GenTb))) #Will be error in DataFrames v0.11.2
            # rlt=DF_hcat(DataFrame(Ontology=fill(ontology,N)), r2j(GenTb))|>x->isempty(rlt) ? x : DF_vcat(rlt,x)
            addrows!(rlt, tb(r2j(GenTb), Ontology=fill(ontology, N)))
         else
            minFDR=callr1("min",FDR)
            println("No enriched item was found. The lowest p_adj($padj_method) is (p_adj=$minFDR): ")
        end
    end
    return (value(rlt), GOdatas, GOtests)
end
function GOenrichFB(fg::AbstractVector{T},bg::AbstractVector{T};
                    ann::AbstractString="org.Hs.eg.db",
                    balance::Tuple{Vector,Vector}=([],[]),
                    ID::AbstractString="ensembl",
                    outfile::AbstractString="",
                    outRobj::Bool=false,
                    wargs...) where {T<:AbstractString}
    if !isempty(balance[1])
        fgv,bgv=balance
        @assert(length(fg)==length(fgv),"Foreground items and their balance values have inconsistent lengths.")
        @assert(length(bg)==length(bgv),"Background items and their balance values have inconsistent lengths.")
        #First remove shared items in background.
        l=.!ismbr(bg,fg)
        bg=bg[l]
        bgv=bgv[l]

        #Remove duplicated items.
        l=uniqr(fg)[1]
        fg=fg[l]
        fgv=fgv[l]
        
        l=uniqr(bg)[1]
        bg=bg[l]
        bgv=bgv[l]
        
        l,t=balancemask2(fgv,bgv)
        bg=bg[l]
        println(f"Balance FG($1):BG($2) 1:$3"(length(fg),length(bg),t))
    end
    # mr=importpkg(:myR)
    @pkgfun(myR,eR,callr)
    eR("library(topGO)")
    callr("library",ann)
    ug=unival(vcatr(fg,bg))
    rlt, GOdata, GOtest=_GOenrichCore(ismemberr(ug,fg),ug,ann,ID;wargs...)
    
    if !isempty(outfile)
        # @pkgfun(DataFrames,writetable)
        # writetable(outfile, rlt)
        if isnothing(rlt)
            writeall(outfile, "No any enriched GO item was found.")
        else
            writetb(outfile, rlt, c"Ontology, GO.ID, Annotated, Significant, Expected, classicFisher, Padj_fdr, Term")
        end
    end
    if outRobj
        return (rlt, GOdata, GOtest)
    else
        return rlt
    end
end
function GOenrichGrp(gene::Vector{T},grp::Group;
                     bg::Vector{T}=gene,outfile::AbstractString="",
                     balance::Tuple{Vector,Vector}=([],[]),
                     wargs...) where {T<:AbstractString}
    # DF=importpkg(:DataFrames)
    # @pkgfun(DataFrames,DataFrame,writetable,hcat=>DF_hcat,vcat=>DF_vcat,size)

    if !isempty(balance[1])
        if length(balance)==1
            allfgv=balance[1]
            bgv=balance[1]
        else
            allfgv,bgv=balance
        end
        rlt,_=grpfunwith(grp,gene,allfgv;input_grpid=true)do grp,fg,fgv
            println("--- Group $grp ---")
            cgo=GOenrichFB(fg,bg;balance=(fgv,bgv),wargs...)
            # Any[DF_hcat(DataFrame(Group=fill(grp,size(cgo,1))),cgo)]
            cgo["Group"]=fill(grp,size(cgo,1))
            cgo
        end
    else
        rlt,_=grpfunwith(grp,gene;input_grpid=true)do grp,fg
            println("--- Group $grp ---")
            cgo=GOenrichFB(fg,bg;wargs...)
            # Any[DF_hcat(DataFrame(Group=fill(grp,size(cgo,1))),cgo)]
            cgo["Group"]=fill(grp,size(cgo,1))
            cgo
        end
    end
    # rlt=DF_vcat(C[.!isempty.(C)]...)
    if !isempty(outfile)
        if isnothing(rlt)
            writeall(outfile, "No any enriched GO item was found.")
        else
            writetb(outfile, rlt, c"Group, Ontology, GO.ID, Annotated, Significant, Expected, classicFisher, Padj_fdr, Term")
        end
    end
    return rlt
end

#}}

#{{ drawann
#[[ drawann drawannlines ]]
#drawann(exn_pos; trans::Group=..., strand='+/-'vec, cds_pos=...;
#        trans_label=(trans, string_vec), txtplot_kw=(fontsize=8), xlim=(left, right))
#drawannlines(pos, styleNO::int_vec, (inp(style1...), inp(style2...), ...); mingap=0, textplot_kw=kw(fontsize=8), transid::Group=1:size(pos,1), showtrans::Bool=false, transplot_inp=inp("k:"), trans_label=([trans=1:N, ]string_vec), label=string_vec(Same rownumber as pos))
#Draw annotation using PyPlot. txtplot_inp is the parameters of function txtplot except the first three arguments.
#In drawannlines, when showtrans is true, function will draw a line for the whole transcriptions at very background. label is used to label every lines while trans_label is used to label every transcriptions.
#The xlim can use to mark the directions when transcrition are out of window.
#See also: inp, segpileup
#Xiong Jieyi, December 9, 2014 > Nov 11, 2015 >Dec 14, 2015>Sep 26, 2016

function drawann(exn_pos::Matrix;trans::Group=[ones(Int,size(exn_pos,1))],
                 strand::Vector{Char}=Char[],cds_pos::Matrix=zeros(Int,size(exn_pos)),
                 mingap::Real=0.05*(maximum(exn_pos[:,2])-minimum(exn_pos[:,1])+1),
                 xlim=nothing,
                 wargs...)
    cds_trans=trans
    l=cds_pos[:,1].>0
    cds_pos=cds_pos[l,:]
    cds_trans=getrow(cds_trans,l)
    int_trans,int_pos=findintron(trans,exn_pos)
    empty_strand_pos=zeros(Int,0,2)
    empty_strand_trans=zerorow(trans)
    if isempty(strand)
        p_strand_pos=empty_strand_pos
        p_strand_trans=empty_strand_trans
        n_strand_pos=empty_strand_pos
        n_strand_trans=empty_strand_trans
    else
        p_rl=rowpile()
        n_rl=rowpile()
        grploop(trans, strand,exn_pos, input_grpid=true) do ctrans, cstrand, cexn_pos
            if cstrand[1]=='+'
                t=maximum(cexn_pos[:,2])
                if xlim!=nothing && t>xlim[2]
                    markpos=xlim[2]-(xlim[2]-xlim[1])*0.02
                    if any(cexn_pos[:,1].<markpos)
                        t=isa(t,Integer) ? round(typeof(t),markpos) : markpos
                    end
                end
                addrow!(p_rl, fill(t,1,2),ctrans)
            elseif cstrand[1]=='-'
                # @assert(cstrand[1]=='-',"Invalid strand.")
                t=minimum(cexn_pos[:,1])
                if xlim!=nothing && t<xlim[1]
                    markpos=xlim[1]+(xlim[2]-xlim[1])*0.02
                    if any(cexn_pos[:,2].>markpos)
                        t=isa(t,Integer) ? round(typeof(t),markpos) : markpos
                    end
                end
                addrow!(n_rl, fill(t,1,2),ctrans)
            else
                if cstrand[1]!='.'
                    @warn f"Unknown annotation strand '$1'. Considered it as no-strand."(cstrand[1])
                end
            end
        end
        (p_strand_pos, p_strand_trans)=isvirgin(p_rl) ? (empty_strand_pos,empty_strand_trans) : value(p_rl)
        (n_strand_pos, n_strand_trans)=isvirgin(n_rl) ? (empty_strand_pos,empty_strand_trans) : value(n_rl)
    end
    style=[1*ones(Int,size(p_strand_pos,1));
           2*ones(Int,size(n_strand_pos,1));
           3*ones(Int,size(int_pos,1));
           4*ones(Int,size(exn_pos,1));
           5*ones(Int,size(cds_pos,1))]
    styledetail=(inp("r>",markersize=4),inp("r<",markersize=4),inp("-",lw=0.5,color="0.5"),inp("b-",linewidth=2),inp("g-",linewidth=4))
    drawannlines([p_strand_pos;n_strand_pos;int_pos;exn_pos;cds_pos],
                 style,styledetail;
                 transid=vcatr(p_strand_trans,n_strand_trans,int_trans,trans,cds_trans),
                 showtrans=false,mingap=mingap,xlim=xlim,wargs...)
    nothing
end

function drawannlines(pos::Matrix{Tn1},style::Vector{Tn2},styledetail::Union{Vector,Tuple};mingap::Real=0,trans_label::Tuple=(),label::Vector=ASCIIString[],txtplot_inp::Tuple{Tuple,Any}=inp(fontsize=8),txtplot_kw=txtplot_inp[2], transid::Group=1:size(pos,1), showtrans::Bool=false, transplot_inp::Tuple{Tuple,Any}=inp("k:"), xlim=nothing) where {Tn1<:Real,Tn2<:Real}
    PP=importpkg(:PyPlot)
    
    (Tpos,Ttrans)=grpfun(x->[minimum(x[:,1]), maximum(x[:,2])],transid,pos)
    Tlayer=segpileup([Tpos[:,1]-(mingap/2) Tpos[:,2]+(mingap/2)])
    layer=dtshift(transid,Ttrans,Tlayer)
    # holdstat=PP.ishold()
    # PP.hold(true)
    if showtrans
        for i=1:length(Tlayer)
            invokelatest(PP.plot, vec(Tpos[i,:]),[Tlayer[i],Tlayer[i]],transplot_inp[1]...;solid_capstyle="butt", transplot_inp[2]...)
        end
    end
    grploop(style,pos,layer,input_grpid=true) do cstyle, cpos, clayer
        for i=1:size(cpos,1)
            invokelatest(PP.plot, vec(cpos[i,:]),fill(clayer[i],2),styledetail[cstyle][1]...;solid_capstyle="butt", styledetail[cstyle][2]...)
        end
    end

    if !isempty(label)
        xx=meanh(pos)
        l=if xlim!=nothing
            xlim[1].<=xx.<=xlim[2]
        else
            trues(length(xx))
        end
        txtplot(xx[l], layer[l], label[l]; verticalalignment="bottom", txtplot_kw...)
    end
    
    if !isempty(trans_label)
        tl=tb()
        if length(trans_label)>1
            tl["transid"]=trans_label[1]
            tl["label"]=trans_label[2]
        else
            tl["label"]=trans_label[1]
            tl["transid"]=collect(1:length(trans_label[1]))
        end
        d"pos,layer"tl=dtshift(tl["transid"],Ttrans,(Tpos,Tlayer))
        tl["midpos"]=meanh(tl["pos"])
        if xlim!=nothing
            tl=rec(tl, xlim[1].<=tl["midpos"].<=xlim[2])
        end
        txtplot(tl["midpos"],tl["layer"],tl["label"];verticalalignment="bottom",txtplot_kw...)
    end
    
    invokelatest(PP.ylim, 0.5,maximum(layer)+0.5)

    if xlim!=nothing
        invokelatest(PP.xlim, xlim)
    end
    nothing
end
export drawann, drawannlines

#}}

#{{ drawcvg
#[[ drawcvg ]]
#drawcvg(read_pos; readid=..., strand=..., lim=(min,max), log=false, is_uniq_map=BitVector, juncno::Vector=[],max_pixel::Integer=10000,ylimby=(bottom, top)|(),minyhight=10,normP::Float=NaN,mixstrand=false,ylim_by_uni_reads=false,fancy_junc=true, rasterized=true)
#Draw coverage figure. readid is used to draw junctions.
#See also: drawann, drawannlines, drawcvgfromfile, drawcvg_core
#Xiong Jieyi, December 10, 2014 >Jun 7, 2015 >Aug 26, 2015 >Mar 3, 2016>Apr 24, 2016>Jun 1, 2016

function drawcvg(pos::Matrix; strand::Vector{Char}=fill('+',size(pos,1)),lim::Array=[minimum(pos[:,1]),maximum(pos[:,2])],readid::Group=[],log::Bool=false,is_uniq_map::Union{BitVector,Vector{Bool}}=trues(size(pos,1)),juncno::Vector=[],max_pixel::Integer=10000,ylimby::Tuple=(),minyhight::Real=6,normP::AbstractFloat=NaN,mixstrand::Bool=false, shadow::Matrix=[], ylim_by_uni_reads::Bool=false, fancy_junc::Bool=true, rasterized::Bool=true)
    PP=importpkg(:PyPlot)
    # @pkgfun(PyPlot, getindex)
    if mixstrand
        strand=fill('+', size(pos,1))
    end
    
    ppos=pos[strand.=='+',:]
    pcvg,xrang=seg2cvg(ppos,lim...; max_pixel=max_pixel)
    pupos=pos[(strand.=='+') .& is_uniq_map,:]
    pucvg=seg2cvg(pupos, xrang)
    
    npos=pos[strand.=='-',:]
    ncvg=seg2cvg(npos, xrang)
    nupos=pos[(strand.=='-') .& is_uniq_map,:]
    nucvg=seg2cvg(nupos, xrang)
    
    xx=[minimum(xrang);xrang;maximum(xrang)]

    if log
        pcvg=log10.(pcvg.+1)
       pucvg=log10.(pucvg.+1)
        ncvg=log10.(ncvg.+1)
        nucvg=log10.(nucvg.+1)
    end
    
    if !isempty(ylimby) || ylim_by_uni_reads
        tpcvg=max.(pcvg,pucvg)
        t=ylim_by_uni_reads ? pucvg : tpcvg
        t=isempty(ylimby) ? t : t[ylimby[1].<=xrang.<=ylimby[2]]
        ylimmax_p=maximum(t)*1.3

        tncvg=max.(ncvg,nucvg)
        t=ylim_by_uni_reads ? nucvg : tncvg
        t=isempty(ylimby) ? t : t[ylimby[1].<=xrang.<=ylimby[2]]
        ylimmax_n=maximum(t)*1.3

        ylimmax=ceil(Int,max(ylimmax_p,ylimmax_n,minyhight))
        isoverflow_p=tpcvg.>ylimmax
        if any(isoverflow_p)
            pcvg[pcvg.>ylimmax].=ylimmax
            pucvg[pucvg.>ylimmax].=ylimmax
        end
        
        isoverflow_n=tncvg.>ylimmax
        if any(isoverflow_n)
            ncvg[ncvg.>ylimmax].=ylimmax
            nucvg[nucvg.>ylimmax].=ylimmax
        end
    else
        ylimmax=Inf
    end
    invokelatest(PP.fill, xx,[0;pcvg;0],"y",
                 xx,[0;-ncvg;0],"y",
                 xx,[0;pucvg;0],"#af8dc3",
                 xx,[0;-nucvg;0],"#7fbf7b",
                 rasterized=rasterized)
    # holdstat=PP.ishold()
    # PP.hold(true)
    if !isempty(readid)
        if !isempty(juncno)
            l=juncno.>0
            readid=getrow(readid,l)
            pos=pos[l,:]
            strand=strand[l]
        end
        if !isempty(pos)
            junc_readid,junc_pos=findintron(readid,pos)
            if !isempty(junc_readid)
                junc_strand=dtshift(junc_readid,readid,strand)
                (uj_freq,(uj_pos,uj_strand))=freq((junc_pos,junc_strand))
                if log
                    uj_freq=log10.(uj_freq.+1)
                end
                uj_freq[uj_strand.=='-']*=-1
                uj_len=seglen(uj_pos)
                # uj_color=fancy_junc ? invokelatest( invokelatest( getindex, invokelatest( getindex, PP.plt, :cm), :jet), (uj_len-minimum(uj_len))./(maximum(uj_len)-minimum(uj_len)))[:,1:3] : []
                if fancy_junc
                    uj_color=f"C$1".(mod.(uj_pos[:, 1]+uj_pos[:, 2], 10))
                end
                for i=sortri(uj_len)
                    yy=uj_freq[i]
                    yy=if yy>0
                        min(yy,ylimmax)
                    else
                        -min(-yy,ylimmax)
                    end
                    invokelatest(PP.plot, vec(uj_pos[i,:]),fill(yy,2),".",markersize=4,lw=0.5,color=fancy_junc ? uj_color[i] : "r")
                end
                for i=sortri(uj_len, rev=true)
                    yy=uj_freq[i]
                    yy=if yy>0
                        min(yy,ylimmax)
                    else
                        -min(-yy,ylimmax)
                    end
                    invokelatest(PP.plot, vec(uj_pos[i,:]),fill(yy,2),"--",markersize=4,lw=0.5,color=fancy_junc ? uj_color[i] : "r")
                end
            end
        end
    end
    if !isempty(ylimby) || ylim_by_uni_reads
        if any(isoverflow_p)
            invokelatest(PP.plot, xrang[isoverflow_p],fill(ylimmax,sum(isoverflow_p)),"c.",markersize=4)
            invokelatest(PP.ylim, (invokelatest(PP.ylim)[1],ylimmax)) #To make the figure tight.
        end
        if any(isoverflow_n)
            invokelatest(PP.plot, xrang[isoverflow_n],fill(-ylimmax,sum(isoverflow_n)),"c.",markersize=4)
            invokelatest(PP.ylim, (-ylimmax,invokelatest(PP.ylim)[2])) #To make the figure tight.
        end
    end
    #Adjust y lim when it is too small
    t=[invokelatest(PP.ylim)...]
    @assert(t[1]<=0 && t[2]>=0)
    ot=copy(t)
    if t[1]<0 && -t[1]<minyhight
        t[1]=mixstrand ? -0.5 : -minyhight
    end
    if t[2]>0 && t[2]<minyhight
        t[2]=minyhight
    end
    t!=ot && invokelatest(PP.ylim, t)
    
    #Add right y-ticks for normalized coverage
    if !isnan(normP)
        t=invokelatest(PP.ylim)
        ca=invokelatest(PP.gca)
        invokelatest(getidx(invokelatest(getidx(ca, :twinx)), :set_ylim), t[1]*normP, t[2]*normP)
        invokelatest(PP.yticks, color="r")
        invokelatest(PP.sca, ca)
    end
    
    # PP.hold(holdstat)
    #PP.ylim(-maximum(ncvg)-1,maximum(pcvg)+1)
    if !isempty(shadow)
        ca=invokelatest(PP.gca)
        yy=invokelatest(PP.ylim)
        for i in 1:size(shadow,1)
            invokelatest(getidx(ca, :add_patch), invokelatest(getidx(PP.plt, :matplotlib, :patches, :Rectangle), (shadow[i,1], yy[1]), shadow[i,2]-shadow[i,1], yy[2]-yy[1], facecolor="#00ffff", edgecolor="none", zorder=0))
        end
    end
    invokelatest(PP.plot, lim, [0, 0], "-", color="0.7", linewidth=0.5, zorder=-1) #0-base line
    invokelatest(PP.xlim, lim...)
end
export drawcvg

#}}

#{{ drawcvgsbypos drawcvgsbygene drawcvgs_core
#[[ drawcvgsbypos drawcvgsbygene drawcvgs_core ]]
# read_vec = drawcvgsbypos(chrno_vec, pos_mx, annfile|anntable_vec, readfile_vec|readtable_vec;
#  samplelabel=String_Vec|fastgrp(...), figfile="file.pdf", readprop=c"...", readfilter=fun(read), title=vec, annlabel=vec, log=false, expand=1, expandnt=..., parallel=false, subfigsize=(8,2), readsaminp=inp(readStarSamFile_params, default: pairend=true), tight=true, ylimfit=false, totalreadnum=vector, mixstrand=false, ylim_by_uni_reads=false, fancy_junc=true, shadow=Any[(chrno_vec, pos_mx), (...), ...], rasterized=true,
# as10x=( Vector{Dict(barcode.=>cell_grp)}, Vector{Vector{cell_grp}} ) ) #Need using BioAlignments
# ... = drawcvgsbygene(geneid_vec, annfile.jld|.mat|anntable_vec, readfile_vec|readtable_vec; ...)
# Draw annotation and coverage figures.
# chrno_vec and pos_mx should have the same row number.
# drawcvgs_core(chrno,pos,ann::Dict,reads_vec::Vector{Dict};samplelabel=vec,vline=...)
# Draw annotation and coverage file by data.
# samplelabel could be a string vector or a fastgrp, with the same length of readfile_vec. The read sets with the same samplelabels will be merged and draw in one panel.
# If expand is not need, set expand=0.
# mixstrand=true: function will ignore the read strand information and just regard all the reads as in positive strand.
# parallel: loading file parallelly. Only useful when input read_filenames rather than read tables.
# vline: A vector of the x axis, where vertical lines will draw.
# ylimfit: ylim assigned by only the center (true. e.g. the candidate gene itself) or the whole coverage (false. i.e. include expanded region).
# ylim_by_uni_reads: The hight of ylim is decided by unique-mapped reads.
# subfigsize: size of each subfigures. The total size will be (subfigsize[1],subfigsize[2]*fignum).
# element of reads should have fields: readid,chrno,pos,strand,is_uniq_map(option)
# ann should have fields: transid, chrno, exn_pos, cds_pos(option), strand(option), label(option)*
# annlalel will be showed at the left of each annotation lane.
# readfilter: given a function, input: read_table, output: BitVector
# readprop: Additional fields needed in readfilter, except readid,chrno,pos,strand,is_uniq_map.
# Function will output the read_vec, which contains all the read used in the drawing regions. To draw another figure in the same regions, please first draw a wider scale coverage then draw narrower scale coverage.
#When readfile is indexed bam files, readsaminp is used to pass parameter to readStarSamFile function.
#tight: Whether show xlabel in every subfigures (but occupy more figure space) or only the first and last subfigures. If tight is false, xaxis in all subfigures will be shalled, which will be benefit for interactive GUI browsing.
#totalreadnum: a vector for the total read number of each read file. If it is given, the RPKM axis will show as the right-y-axis.
#fancy_junc: When it is true, draw junction lines in rainbow colors for different lengths.
#shadow: A vector of tuple with the length of coverage panel number. Each tuple composited a chrno list and n x 2 pos matrix.
# The label in ann table will label each transcription in the figure. A transcription with any of its exon having a label will be labeled. Labels for other exons can be "".
# as10x: Treat bam file as 10x data, and draw coverage for cell groups. The first and second parameter is a vector of (Dict and Vector, seperately), for each bam files. When as10x is used, the samplelabel should equal to the length of all cell groups number. Note that in as10x mode, using BioAlignments is reqired.
#See also: readStarSamFile
#Xiong Jieyi, December 11, 2014>March 18, 2015>May 31, 2015>Aug 22, 2015>Aug 26, 2015>Nov 5, 2015 >Nov 11, 2015>Dec 22, 2015>Apr 24, 2016>Jul 11, 2016>Feb 10, 2017>22 Jan 2018>23 Apr 2019>25 Sep 2019

export drawcvgsbypos, drawcvgsbygene, drawcvgs_core
  #{{ drawcvgsbypos

function drawcvgsbypos(
    chrno::Group,pos::Matrix{T1},
    ann::Vector,
    readfiles::Vector{rdType};
    parallel::Bool=false,
    readsaminp=inp(pairend=true),
    readprop::Array=AbstractString[],
    readfilter::Function=x->trues(rownum(x)),
    expand::Real=1,
    expandnt::Real=-1,
    as10x::Union{Tuple, Nothing}=nothing,
    wargs...) where {T1<:Real,rdType<:AbstractString}
    @pkgfun(PyCall, getindex, getproperty)

    if !isnothing(as10x)
        return _drawcvgsbypos_10x(
            chrno, pos, ann, readfiles, as10x...; readsaminp=readsaminp, readprop=readprop,
            readfilter=readfilter, expand=expand, expandnt=expandnt, wargs...)
    end
    
    opos=deepcopy(pos)
    if expandnt>=0
        pos=[pos[:,1]-expandnt pos[:,2]+expandnt]
    elseif expand>0
        t=(pos[:,2]-pos[:,1]+1)*expand/2
        pos=[pos[:,1]-t pos[:,2]+t]
    end
    
    if isempty(readfiles)
        rds=Dict[]
    else
        if occursin(r"\.(bam|sam)$",readfiles[1])
            loadfun=x->begin
                T=fd(readStarSamFile(x;region=(chrno,[floor.(Int,pos[:,1]-100000) ceil.(Int,pos[:,2]+100000)]),strand=true,onlyuniq=false,readidfun=string,readsaminp[2]...),[c"readid,chrno,pos,strand,juncno,mapped_loci_num";readprop])
                T["is_uniq_map"]=.!(T["mapped_loci_num"].>1)
                T
            end
        elseif endswith(readfiles[1], ".jld")
            loadfun=x->jldload(x,[c"readid,chrno,pos,strand,is_uniq_map,juncno";readprop])
        else
            loadfun=x->jload(x,[c"readid,chrno,pos,strand,is_uniq_map,juncno";readprop])
        end
        if parallel
            println("Loading $(length(readfiles)) files parallelly...")
            rds=pmapwd(readfiles,env=(chrno,pos)) do rfn,chrno,pos
                crd=loadfun(rfn)
                crd=rec(crd,readfilter(crd))
                (ignore,t)=genomemap((chrno,pos),d"chrno,pos"crd,touch=true)
                rec(crd, ismbr(crd["readid"], Set(eachrow(getrow(crd["readid"], unival(t))))))
            end
        else
            rds=map(readfiles) do rfn
                println("Loading read file $rfn")
                crd=loadfun(rfn)
                crd=rec(crd,readfilter(crd))
                (ignore,t)=genomemap((chrno,pos),d"chrno,pos"crd,touch=true)
                rec(crd, ismbr(crd["readid"], Set(eachrow(getrow(crd["readid"], unival(t))))))
            end
        end
    end
    drawcvgsbypos(chrno,opos,ann,rds;expand=expand,expandnt=expandnt,wargs...)
end
function _drawcvgsbypos_10x(
    chrno::Group,pos::Matrix{T1},
    ann::Vector,
    readfiles::AbstractVector{T2},
    barcode::AbstractVector{Dict{T3, T4}},
    cellgrpid::AbstractVector{T5};
    UMI::String="UB",
    smpBC::String="CB",
    readsaminp=inp(pairend=true),
    readprop::Array=AbstractString[],
    readfilter::Function=x->trues(rownum(x)),
    expand::Real=1,
    expandnt::Real=-1,
    samplelabel::Union{AbstractVector, fastgrp}=map(x->"#$x", 1:length(readfiles)),
    wargs...) where {T1<:Integer, T2<:AbstractString, T3<:AbstractString, T4<:AbstractString, T5<:AbstractVector}
    # 30 Sep 2019: using @il has a problem here, as the for ... in ... calls iterate() implicitly, which cannot be invokelatest-ed by @il.
    
    iselequal((length(readfiles), length(barcode), length(cellgrpid))) || error("In 10x mode, inputed bam file number, barcode table number and cell_group_id list number does not match.")
    
    if isa(samplelabel, fastgrp)
        error("In 10x mode, samplelabel=fastgrp(...) is not supported yet.")
    end
    if length(samplelabel)==length(readfiles)
        C=map(samplelabel, cellgrpid) do x1, x2
            # @BK x1 x2 samplelabel cellgrpid
            x1*":".*x2
        end
        samplelabel=[C...;]
    else
        if length(samplelabel)!=sum(length.(cellgrpid))
            error("The length of samplelabel should either equals to the bam file number or equals to the file_number X cell_group_number.")
        end
    end
    
    opos=deepcopy(pos)
    if expandnt>=0
        pos=[pos[:,1]-expandnt pos[:,2]+expandnt]
    elseif expand>0
        t=(pos[:,2]-pos[:,1]+1)*expand/2
        pos=[pos[:,1]-t pos[:,2]+t]
    end

    peudoRdP=Ref(1)
    if isempty(readfiles)
        rds=Dict[]
    else
        rds=mapr(readfiles, barcode, cellgrpid, mrows=true) do rfn, cbarcode, ccellgrpid
            println("Loading read file $rfn")
            _read10xBam(rfn, cbarcode, ccellgrpid, chrno, pos, smpBC, UMI, peudoRdP)
        end
    end
    drawcvgsbypos(chrno,opos,ann,rds;expand=expand,expandnt=expandnt, samplelabel=samplelabel, wargs...)
end
function _read10xBam(bamfn::AbstractString, cellBC::Dict, cellgrpid::AbstractVector, chrno::AbstractVector, pos::AbstractMatrix, smpBC::String, UMI::String, peudoRdP::Base.RefValue{Int64})
    BA=importpkg(:BioAlignments, preloaded=true)
    rddict=Dict{String, rowpile}()
    rdtmp=Dict{String, rowpile}()
    for x in cellgrpid
        rddict[x]=rowpile()
    end
    cellBC=filter(((k, v),)->haskey(rddict, v), cellBC)
    Pjuncno=1
    # umipool=Set{String}()
    # prdpos=-1
    if !isfile(bamfn*".bai")
        println("Bam file $bamfn haven't be indexed yet. Indexing....")
        sf"samtools index $1"(bamfn)
    end
    readno=0
    open(BA.BAM.Reader, bamfn, index=bamfn*".bai") do bamO
        bamrefset=Set(bamO.refseqnames) #@il
        rowloop(chrno, pos) do cchrno, cpos
            chr=no2chr(cchrno, chr=false)
            in(chr, bamrefset) || return nothing
            println(f"$1:$2-$3"(chr, floor.(Int, cpos[1]-100000), ceil.(Int, cpos[2]+100000)))
            for rd in BA.eachoverlap(bamO, chr, floor.(Int, cpos[1]-100000):ceil.(Int, cpos[2]+100000))
                BA.haskey(rd, UMI) || continue
                BA.haskey(rd, smpBC) || continue
                rd_cellBC=cut"-1"(rd[smpBC]::String)
                cellgrp=get(cellBC, rd_cellBC, nothing)
                isnothing(cellgrp) && continue
                readno+=1
                
                rdpos=parseCIGARpos(BA.BAM.cigar(rd), BA.BAM.position(rd))
                t=tb(readid=readno,
                     chrno=cchrno,
                     strand=BA.BAM.ispositivestrand(rd) ? '+' : '-',
                     is_uniq_map=rd["NH"]::UInt8==0x01,
                     UMI=rd[smpBC]::String)
                T=reprow(t, size(rdpos, 1))
                T["pos"]=rdpos
                T["juncno"]=(t=rownum(T);)==1 ? 0 : (1:t)
                if !haskey(rdtmp, cellgrp)
                    rdtmp[cellgrp]=rowpile()
                end
                addrows!(rdtmp[cellgrp], T)
            end
            
            for (k, v) in rdtmp
                if !haskey(rddict, k)
                    rddict[k]=rowpile()
                end
                vT=value(v)
                grploop(vT["UMI"], fd_i(vT, "UMI"), id=true) do cumi, T
                    if any(T["is_uniq_map"])
                        T=getrow(T, T["is_uniq_map"])
                    end
                    cstrand=freq(T["strand"], revsort=true)[2][1]
                    T=getrow(T, T["strand"].==cstrand)
                    T=getrow(T, sortri(T["pos"]))
                    isJcPoint::Matrix{Bool}=grpfunexp(T["readid"], T["juncno"]) do x
                        t=[x x].>0
                        t[1, 1]=false
                        t[end, 2]=false
                        t
                    end
                    L=Int[1]
                    for i=[2:rownum(T); 0]
                        if i>0 && T["pos"][i, 1]<=T["pos"][i-1, 2]
                            push!(L, i)
                        else
                            pos2=T["pos"][L, 2]
                            P2=if any((t=isJcPoint[L, 2];))
                                pos2=pos2[t]
                                pos2[sortri((idennum(pos2), pos2))[end]]
                            else
                                maximum(pos2)
                            end
                            pos1=T["pos"][L, 1]
                            P1=if any((t=isJcPoint[L, 1];))
                                pos1=pos1[t]
                                pos1[sortri((idennum(pos1), -pos1))[end]]
                            else
                                peudoRdP.x+=1
                                minimum(pos1)
                            end
                            crd=tb(
                                readid=f"$1_$2"(cumi, peudoRdP.x),
                                chrno=cchrno,
                                pos=[P1 P2],
                                is_uniq_map=any(T["is_uniq_map"][L]),
                                juncno=maximum(T["juncno"][L]),
                                strand=cstrand)
                            #Note here the juncno may not be in prefect order. Anyhow, it will not disturb the final figure since the drawcvg() only check whether juncno==0. 27 Sep 2019
                            addrow!(rddict[k], crd)
                            L=Int[i]
                        end
                    end
                end
            end
            empty!(rdtmp)
        end #rowloop
    end #open bam file
    map(cellgrpid) do x
        if isnothing((t=value(rddict[x]);))
            tb(readid=zerorow(String),
               chrno=zerorow(Int),
               pos=zerorow(Int, 2),
               strand=zerorow(Char),
               juncno=zerorow(Int),
               is_uniq_map=zerorow(Bool))
        else
            t
        end
    end
end

function drawcvgsbypos(chrno::Group,pos::Matrix{T1},ann::Vector{annType},readfiles::Vector{rdType};
                       samplelabel::Union{AbstractVector, fastgrp}=map(x->"#$x", 1:length(readfiles)),
                       annlabel::Vector{T2}=convert(Vector{ASCIIString},map(x->"ann$x",1:length(ann))),
                       readsaminp=inp(), #No use here.
                       figfile::AbstractString="",dpi::Int=0,title::Vector=[],expand::Real=1,expandnt::Real=-1,log::Bool=false,subfigsize::Tuple{Real,Real}=(8,2),ylimfit::Bool=false,totalreadnum::AbstractVector{Ttr}=fill(0,length(readfiles)),mixstrand::Bool=false,
                       shadow::Vector=Any[],
                       ylim_by_uni_reads::Bool=false,
                       fancy_junc::Bool=true,
                       rasterized::Bool=true) where {T1<:Real,T2<:AbstractString,annType,rdType<:Union{Dict, Nothing},Ttr<:Real}
    
    vlinepos=deepcopy(pos)
    if expandnt>=0
        pos=[pos[:,1]-expandnt pos[:,2]+expandnt]
    elseif expand>0
        expandnt=(pos[:,2]-pos[:,1]+1)*expand/2
        pos=[pos[:,1]-expandnt pos[:,2]+expandnt]
    end
    
    length(totalreadnum)==length(readfiles) || error("totalreadnum should have the same length as readfiles.")
    
    needmerge=isa(samplelabel, fastgrp) || !isuniqr(samplelabel)
    if needmerge
        fg=isa(samplelabel, fastgrp) ? samplelabel : fastgrp(samplelabel, stable=true)
        rds=Dict[]
        rdp=0
        grploop(fg, readfiles) do creadfiles
            rd=nothing
            for fi=1:length(creadfiles)
                crd=creadfiles[fi]
                t1,t2=uniqr(crd["readid"])
                crd["readid"]=t2+rdp
                rdp=rdp+length(t1)
                
                rd=vcatr(rd,crd)
            end
            push!(rds,rd)
        end
        samplelabel=fg.grpid
        readnumsum=grpfun(sum,fg,totalreadnum)
    else
        length(samplelabel)==length(readfiles) || error("sample and samplelabel should be in the same length.")
        rds=readfiles
        readnumsum=totalreadnum
    end
    
    if isempty(ann)
        ann=Dict[]
    elseif annType<:AbstractString
        ann=map(ann) do fn
            if endswith(fn, ".jld")
                jldload(fn, c"geneid,transid,chrno,exn_pos,cds_pos,strand")
            else
                jload(fn, c"geneid,transid,chrno,exn_pos,cds_pos,strand")
            end
        end
    end
    
    println("Drawing...")
    PP=importpkg(:PyPlot)
    
    if !isempty(figfile)
        invokelatest(PP.pygui, false)
        py_backend_pdf=importpy("matplotlib.backends.backend_pdf")
        if !ismatch(r"\.pdf$",figfile)
            figfile=figfile*".pdf"
        end
        pdf=invokelatest(invokelatest(getproperty, py_backend_pdf, :PdfPages), figfile)
    end
    for i=1:size(pos,1)
        drawcvgs_core(getrow(chrno,i),pos[i,:],ann,rds,samplelabel=samplelabel,log=log,annlabel=annlabel,vline=vlinepos[i,:],subfigsize=subfigsize,title=isempty(title) ? "" : title[i],ylimby=ylimfit ? (vlinepos[i,1],vlinepos[i,2]) : (),totalreadnum=readnumsum,mixstrand=mixstrand, shadow=shadow, ylim_by_uni_reads=ylim_by_uni_reads, fancy_junc=fancy_junc, rasterized=rasterized)
        if !isempty(figfile)
            invokelatest(PP.rc,"pdf",fonttype=42)
            invokelatest(getidx(pdf, :savefig), dpi=300)
            invokelatest(PP.close)
        end
        println("Finished $i / $(size(pos,1)).")
    end
    if !isempty(figfile)
        invokelatest(getidx(pdf, :close))
        if isijulia()
            invokelatest(PP.pygui, true)
        end
        
        #Below codes may considered be removed.
        if dpi>0
            println("Rasterizing vector graphic to $dpi dpi...")
            run(`convert -density $dpi -compress Zip $figfile $figfile`)
        end
        #end
        println("Figure is saved to $figfile")
    end
    return readfiles
end
drawcvgsbypos(chrno::Group,pos::Matrix{T1},ann::Union{Dict,AbstractString},args...;wargs...) where {T1<:Real} = drawcvgsbypos(chrno,pos,[ann],args...;wargs...)

  #}}
  #{{ drawcvgsbygene
function drawcvgsbygene(genes::Vector{T1},anns::Vector{T2},args...;wargs...) where {T1<:AbstractString,T2<:Dict}
    fg=fastgrp(anns[1]["geneid"],genes)
    chrnos,poss=grpfun(fg,anns[1],multi_output=true)do cann
        (cann["chrno"][1],cann["exn_pos"]|>x->[minimum(x[:,1]) maximum(x[:,2])])
    end 
    drawcvgsbypos(chrnos,poss,anns,args...;title=genes,wargs...)
end
drawcvgsbygene(genes::Vector{T1},anns::Dict,args...;wargs...) where {T1<:AbstractString} =drawcvgsbygene(genes,[anns],args...;wargs...)
function drawcvgsbygene(genes::Vector{T},annfile::AbstractString,args...;wargs...) where {T<:AbstractString}
    if splitext(annfile)[2]==".mat"
        drawcvgsbygene(genes::Vector{T},[matload(annfile,c"geneid,transid,chrno,exn_pos,cds_pos,strand")],args...;wargs...)
    else
        drawcvgsbygene(genes::Vector{T},[jldload(annfile,c"geneid,transid,chrno,exn_pos,cds_pos,strand")],args...;wargs...)
    end
end
  #}}
  #{{ drawcvgs_core

function drawcvgs_core(chrno::Real,pos::Array{T1},anns::Vector{T2},rds::Vector{T3};
                       samplelabel::Vector=map(x->"Sample$x",1:length(rds)),
                       log::Bool=false,annlabel::Vector{T4}=map(x->"Ann$x",1:length(anns)),
                       vline::Array=[],subfigsize::Tuple{Real,Real}=(8,2),
                       title::AbstractString="",
                       xlabel::AbstractString=f"$1:$2-$3"(no2chr(chrno), round(Int, pos[1]), round(Int, pos[2]))*( isempty(vline) ? "" : " vline-length: "*num2str(vline[2]-vline[1]+1) ),
                       tight::Bool=true,
                       max_pixel::Integer=10000,
                       ylimby::Tuple=(),
                       totalreadnum::AbstractVector{Ttr}=fill(0,length(rds)),
                       mixstrand::Bool=false,
                       shadow::Vector=Any[],
                       ylim_by_uni_reads::Bool=false,
                       fancy_junc::Bool=true,
                       rasterized::Bool=true
                       ) where {T1<:Real,T2<:Dict,T3<:Union{Dict, Nothing},T4<:AbstractString,Ttr<:Real}
    # rd should have fields: readid,chrno,pos,strand,is_uniq_map(option)
    # ann should have fields: transid, chrno, exn_pos, cds_pos, strand
    
    # @pkgfun(PyPlot, getindex)
    if all(totalreadnum.<=0)
        normP=fill(NaN,length(rds))
    else
        normP=10^9 ./ totalreadnum
    end
    
    isdefined(Main,:PyPlot) || require("PyPlot")
    PP=Main.PyPlot
    
    fignum=length(anns)+length(rds)
    invokelatest(PP.figure,figsize=(subfigsize[1],subfigsize[2]*fignum))

    ax=nothing #The first subfigure
    axs=Any[]  #Other subfigures
    
    #Draw annotation
    for i=1:length(anns)
        if isa(ax,Void)
            ax=invokelatest(PP.subplot,fignum,1,i)
        elseif tight
            push!(axs,invokelatest(PP.subplot,fignum,1,i))
        else
            push!(axs,invokelatest(PP.subplot,fignum,1,i,sharex=ax))
        end
        t=rec(anns[i],(anns[i]["chrno"].==chrno) .& (anns[i]["exn_pos"][:,2].>pos[1]) .& (anns[i]["exn_pos"][:,1].<pos[2]))
        cann=rec(anns[i],ismemberr(anns[i]["transid"],t["transid"]))
        mingap=0.05*(pos[2]-pos[1]+1)
        drawannkw=Any[]
        if haskey(cann,"strand")
            push!(drawannkw, :strand=>cann["strand"])
        end
        if haskey(cann,"cds_pos")
            push!(drawannkw, :cds_pos=>cann["cds_pos"])
        end
        if haskey(cann,"label")
            l=uniqr(getrow((cann["transid"],cann["label"]),cann["label"].!=""))[1]
            push!(drawannkw,:trans_label=>(getrow(cann["transid"],l),cann["label"][l]))
        end
        if rownum(cann)>0
            drawann(cann["exn_pos"];trans=cann["transid"],mingap=mingap,xlim=pos,drawannkw...)
        end
        if !isempty(vline)
            t=invokelatest(PP.ylim)
            invokelatest(PP.vlines,vec(vline),t...,colors="y")
            invokelatest(PP.ylim,t)
        end
        invokelatest(PP.ylabel,annlabel[i])

        #Retick coordinate
        if i==1
            xtk=invokelatest(PP.xticks)[1]|>x->x[pos[1].<=x.<=pos[2]]
            xtks=map(x->num2str(round(Int, x)), xtk)
            if uninum(length.(xtks))==1
                N=findfirst(i->uninum(map(x->x[i], xtks))>1, 1:length(xtks[1]))
                
                if N!=nothing && N>1
                    xtks[1]=xtks[1][1:N-1]*","*xtks[1][N:end]
                    for i=2:length(xtks)
                        xtks[i]=xtks[i][N:end]
                    end
                end
            end
            invokelatest(PP.xticks, xtk, xtks)
        end

        invokelatest(PP.xlim,pos...)
    end
    
    #Draw coverage
    for i=1:length(rds)
        if isa(ax,Void)
            ax=invokelatest(PP.subplot,fignum,1,i+length(anns))
        elseif tight
            push!(axs,invokelatest(PP.subplot,fignum,1,i+length(anns)))
        else
            push!(axs,invokelatest(PP.subplot,fignum,1,i+length(anns),sharex=ax))
        end
        rd=rds[i]
        if !isnothing(rd)
            rd=rec(rd,(rd["chrno"].==chrno) .& (rd["pos"][:,2].>pos[1]-100000) .& (rd["pos"][:,1].<pos[2]+100000)) #To reduce memory requirement of the next step.
            rd=rec(rd, ismbr(rd["readid"], Set(eachrow(getrow(rd["readid"], (rd["chrno"].==chrno) .& (rd["pos"][:,2].>pos[1]) .& (rd["pos"][:,1].<pos[2])))))) #To show complete junction lines
            if !isempty(shadow) && !isempty(shadow[i])
                l=(shadow[i][1].==chrno) & (shadow[i][2]|>x->(x[:,2].>=pos[1]) & (x[:,1].<=pos[2]))
                cshadow=shadow[i][2][l,:]
            else
                cshadow=zerorow(eltype(pos),2)
            end
            if haskey(rd,"is_uniq_map")
                drawcvg(rd["pos"];readid=rd["readid"],strand=rd["strand"],is_uniq_map=rd["is_uniq_map"],lim=pos,log=log,juncno=rd["juncno"],max_pixel=max_pixel,ylimby=ylimby,normP=normP[i],mixstrand=mixstrand,shadow=cshadow,ylim_by_uni_reads=ylim_by_uni_reads,fancy_junc=fancy_junc, rasterized=rasterized)
            else
                drawcvg(rd["pos"];readid=rd["readid"],strand=rd["strand"],lim=pos,log=log,juncno=rd["juncno"],max_pixel=max_pixel,ylimby=ylimby,normP=normP[i],mixstrand=mixstrand,shadow=cshadow,ylim_by_uni_reads=ylim_by_uni_reads,fancy_junc=fancy_junc, rasterized=rasterized)
            end
        end

        #Draw vertical lines
        if !isempty(vline)
            # PP.hold(true)
            t=invokelatest(PP.ylim)
            invokelatest(PP.vlines,vec(vline),t..., linestyles=":", colors="b", linewidth=0.5)
            invokelatest(PP.ylim,t)
        end
        
        invokelatest(PP.ylabel,samplelabel[i])
        invokelatest(PP.xlim,pos...)
    end
    
    #Retick coordinate
    xtk=invokelatest(PP.xticks)[1]|>x->x[pos[1].<=x.<=pos[2]]
    xtks=map(x->num2str(round(Int, x)), xtk)
    if uninum(length.(xtks))==1
        N=findfirst(i->uninum(map(x->x[i], xtks))>1, 1:length(xtks[1]))
        
        if N!=nothing && N>1
            xtks[1]=xtks[1][1:N-1]*","*xtks[1][N:end]
            for i=2:length(xtks)
                xtks[i]=xtks[i][N:end]
            end
        end
    end
    invokelatest(PP.xticks, xtk, xtks)
    invokelatest(PP.xlim,pos...)
    
    if !isempty(title)
        # PP.suptitle(title)
        invokelatest(getidx(ax, :set_title), title)
    end
    if !isempty(xlabel)
        t=isempty(axs) ? ax : axs[end]
        invokelatest(getidx(t, :set_xlabel), xlabel)
    end

    if tight
        if !isempty(axs)
            if isempty(title)
                invokelatest(getidx(ax, :xaxis, :tick_top))
            else
                #If title is given, set top-ticklabel as empty.
                invokelatest(getidx(ax, :xaxis, :set_ticklabels), [])
            end
            for h in axs[1:end-1]
                invokelatest(getidx(h, :xaxis, :set_ticklabels), [])
            end
        end
        invokelatest(PP.tight_layout)
        invokelatest(PP.subplots_adjust,hspace = .05)
    else
        invokelatest(PP.tight_layout)
    end
end

 #}}
#}}

#{{ parseDemulReport
#[[ parseDemulReport ]]
#T = parseDemulReport(demultiplex_path)
#Parse demultiplex_path/Reports and get information of %>=Q30, MeanQual and PASSchastity. If a sample were sequenced in multiple lanes, this funciton will output the mean value.
#Output fields: lane,PF_clusters,%ofTheLane,%PrefectBarcode,%OneMismatchBarcode,Yield_Mbases,%PF_clusters,%>=Q30bases,meanQualityScore
#See also: parseGTF
#Xiong Jieyi, Feb 5, 2016

export parseDemulReport
function parseDemulReport(reportpth)
    pth=getfiles(joinpath(reportpth,"Reports/html/"),typ=:dir)[1][1]
    t1,t2=getfiles(pth,typ=:dir)
    pth=t1[findonly(!ismemberr(t2,c"all,default"))]
    t1,t2=getfiles(pth,typ=:dir)
    l=t2.!="all"
    smp=tb(ID=t2[l],
           lane_report_html = map(f"$1/all/lane.html",t1[l]))

    t=rowfunwith(smp) do csmp
        txt=readchomp(csmp["lane_report_html"])
        rpl=rowpile()
        pso=ParseStr(txt)
        next(pso,"Lane Summary","Score</th>","</tr>")
        while next(pso,(:check,"<tr>"))[1]
            C=tb(c"lane,PF_clusters,%ofTheLane,%PrefectBarcode,%OneMismatchBarcode,Yield_Mbases,%PF_clusters,%>=Q30bases,meanQualityScore",
                 next(pso,
                     "<td>",("</td>",int), #Lane
                     # "<td>",("</td>",ascii), #Barcode
                     "<td>",("</td>",x->int(replace(x, "," => ""))), #PF Clusters
                     "<td>",("</td>",x->float(x)/100), # % of the lane
                     "<td>",("</td>",x->float(x)/100), # % Perfect barcode
                     "<td>",("</td>",x->float(x)/100), # % One mismatch barcode
                     "<td>",("</td>",x->int(replace(x, "," => ""))), #Yield (Mbases)
                     "<td>",("</td>",x->float(x)/100), #% PF Clusters
                     "<td>",("</td>",x->float(x)/100), #% >= Q30 bases
                     "<td>",("</td>",float), #Mean Quality Score
                     "</tr>"))
            C["sample"]=csmp["ID"]
            addrow!(rpl,C)
        end
        value(rpl)
    end
    t[1]
end
#}}

#{{ DESeq
#[[ DESeq DESeqMF ]]
# Table = DESeq(A_counts::Matrix|Vector, B_counts::Matrix|Vector; norep=false, sizeFactors=(vecA, vecB))
#     --Run DESeq to test different gene expression.
# (FDRs, raw_p-values, log2_coef::Dict) = DESeqMF( dataMatrix, keyFactor, otherFactor1, otherFactor2, ...; sizeFactors=vec)
#     --Run DESeq to test if the key factor model is better than non-key-factor model, with other factors under consideration.
# If no any repeats, set norep=true to use "blind" method in estimateDispersions.
# Each row of counts is a gene, while each column of counts is a sample. nbinomTest will do between A and B.
# DESeq output table contains fields:
# "id"            
# "baseMean"      
# "baseMeanA"     
# "baseMeanB"     
# "foldChange"    
# "log2FoldChange"
# "pval"          
# "padj"
# In the outpus, the positive direction of fold change is B-A in DESeq().
# For the log2_coef output of DESeqMF(), the alphabetical smallest key factor will not be given.
#
#See also: edgeR_anova, GOenrichFB, GOenrichGrp, myR, GSVA
#Xiong Jieyi, Mar 10, 2016 > 23 Jul 2017 >24 Jul 2017 >31 May 2018 >16 Aug 2018

export DESeq, DESeqMF
function DESeq(A::AbstractArray{T1}, B::AbstractArray{T2}; norep::Bool=false, sizeFactors::Tuple=()) where {T1<:Real, T2<:Real}
    @pkgfun(myR,callr,r2j,eR,rlibrary)
    X=hcat(A,B)
    L=[fill(:A,size(A,2));fill(:B,size(B,2))]
    rlibrary("DESeq")
    cds=callr("newCountDataSet",X,L)
    if isempty(sizeFactors)
        cds=callr("estimateSizeFactors", cds)
    else
        # callr("`sizeFactors<-`", cds, vcat(sizeFactors...))
        cds=callr("function(x, y){sizeFactors(x)<-y;x}", cds, vcat(sizeFactors...)|>x->x./mean(x))
    end
    if norep
        cds=callr("estimateDispersions",cds,method="blind",sharingMode="fit-only")
    else
        cds=callr("estimateDispersions",cds)
    end
    r2j(callr("nbinomTest",cds,"A","B"), NA=NaN, keepvec=true)
end
function DESeqMF(M::AbstractMatrix{T1}, keyF::AbstractVector{T2}, Fs::AbstractVector...; sizeFactors::AbstractVector=[]) where {T1<:Real, T2<:Union{String, Symbol}}
    @pkgfun(myR,callr,callrj,eR,r2j,rnames,relemj,rlibrary)
    if isempty(Fs)
        fml1=eR("count~keyF")
        fml0=eR("count~1")
        L=callr("data.frame", keyF=keyF)
    else
        Flb=f"F$1".(1:length(Fs))
        t=join(Flb,"+")
        fml1=eR("count~$t+keyF")
        fml0=eR("count~$t")
        L=callr("data.frame", tb(tb(Flb,Fs),keyF=keyF))
    end
    
    rlibrary("DESeq")
    cds=callr("newCountDataSet",M,L)
    if isempty(sizeFactors)
        cds=callr("estimateSizeFactors",cds)
    else
        # callr("`sizeFactors<-`", cds, sizeFactors)
        cds=callr("function(x, y){sizeFactors(x)<-y;x}", cds, sizeFactors./mean(sizeFactors))
    end
    cds=callr("estimateDispersions",cds, method="pooled-CR", modelFormula=fml1)
    fit1=callr("fitNbinomGLMs", cds, fml1)
    fit0=callr("fitNbinomGLMs", cds, fml0)
    pval=callr("nbinomGLMTest",fit1,fit0)
    kys=filter(r"^keyF", rnames(fit1))
    T=rename!(relemj(fit1, kys, NA=NaN, keepvec=true), r"^keyF", "")
    (r2j(callr("p.adjust",pval,method="fdr"), keepvec=true), r2j(pval, keepvec=true), T)
end
#}}

#{{ edgeR_anova
#[[ edgeR_anova ]]
#Syntax: Table, comparison_group = edgeR_anova(Readnum_Matrix, Group; lib!size=vec, norm!factors=vec, remove!zeros=false)
#Run edgeR ANOVA mode to test different gene expression. 
#Each row is a gene and each column is a sample. Group length should equal to matrix column number. Group could be a vector of string, symbol or integer.
# lib!size: numeric vector giving the total count (sequence depth) for each library.
# norm!factors: numeric vector of normalization factors that modify the library sizes.
# remove!zeros: logical, whether to remove rows that have 0 total count.
# See ?DGEList in edgeR.
#Returns fields:  "logFC[_group...]", "LR", "PValue", "FDR"
# Updates in 5 Aug 2019:
# The 2nd output "comparison_group" is the "comparison" field of edgeR result. When Group has only two factors, it is the "higher" factor, in String or in Symbol depending on the input, showing the direction of logFCis this field - other filed. (I have no idea about this factor so far when Group has three or more factors.)
#See also: DESeq, DESeqMF, GOenrichFB, GOenrichGrp, GSVA
#Xiong Jieyi, Jun 7, 2016 > 2 Jun 2018 > 5 Aug 2019

export edgeR_anova
function edgeR_anova(M::AbstractMatrix{T1}, G::AbstractVector{T2}; wargs...) where {T1<:Real,T2<:Union{AbstractString,Symbol,Integer}}
    rnum(G)==size(M,2) || error("Group length should equal to matrix column number.")
    @pkgfun(myR, callr, eR, relem, relem1, j2r, r2j, rlibrary)
    rlibrary("edgeR")
    y=callr(:DGEList,counts=j2r(M,row!names=1:size(M,1)),
            group=callr(:factor,G); wargs...)
    design=callr(:model!matrix,eR("~group"),data=relem(y,"samples"))
    y=callr(:estimateGLMCommonDisp,y,design)
    y=callr(:estimateGLMTrendedDisp,y,design)
    y=callr(:estimateGLMTagwiseDisp,y,design)
    fit=callr(:glmFit,y,design)
    lrt=callr(:glmLRT,fit,coef=2:uninum(G))
    # tag=callr(:topTags,lrt,n=size(M,1))["table"]
    toptags=callr(:topTags,lrt,n=size(M,1))
    tag=relem(toptags, "table")
    T=r2j(tag, NA=NaN, keepvec=true) #Changed on 6 Nov 2018
    
    rowidx=int(r2j(callr(:row!names,tag), keepvec=true))

    compgrp=relem1(toptags, "comparison") #Now I have no idea what compgrp will be when facor number >=3. 5 Aug 2019
    if isa(compgrp, String)
        compgrp=replace(compgrp, r"^group"=>"")
        if eltype(G)<:Symbol
            compgrp=Symbol(compgrp)
        end
    end
    (dtshift(1:size(M,1),rowidx,T), compgrp)
end
#}}

#{{ GSVA
#[[ GSVA ]]
# sig1, sig2, ... = GSVA(Matrix, geneid, gene_list1[, gene_list2, ...]; counts=false, pn=1, method="gsva|ssgsea|zscore|plage")
# ... = GSVA(Matrix, Bool_vec1=ismbr(geneid, gene_list)[, Bool_vec2...]; ... )
# Calculate siguature value using R-GSVA pacakge.
# The matrix of first input: gene for row, sample for column. The output is one or more vectors with the length of sample number. If the input is read numbers rather than expression values, set counts=true.
# pn is the parallel computation core number.
# method="ssgsea" will be much faster than the default method "gsva".
# See the manual on https://www.rdocumentation.org/packages/GSVA/versions/1.20.0/topics/gsva
# `method`: Method to employ in the estimation of gene-set enrichment scores per sample. By default this is set to `gsva` (Hänzelmann et al, 2013) and other options are `ssgsea` (Barbie et al, 2009), `zscore` (Lee et al, 2008) or `plage` (Tomfohr et al, 2005). The latter two standardize first expression profiles into z-scores over the samples and, in the case of zscore, it combines them together as their sum divided by the square-root of the size of the gene set, while in the case of plage they are used to calculate the singular value decomposition (SVD) over the genes in the gene set and use the coefficients of the first right-singular vector as pathway activity profile.
# See also: DESeq, edgeR_anova, GOenrichFB, GOenrichGrp
# Jieyi Xiong, 13 Sep 2018 > 27 Feb 2019

export GSVA
function GSVA(count::Matrix, geneid::Vector{String}, gene_list1::Vector{String}, gene_lists::Vector{String}...; counts::Bool=false, pn::Int=1, kwargs...)
    if counts && !(eltype(count)<:Integer)
        error("Counts should be integers when counts=true.")
    end
    if eltype(count)<:Integer
        counts || @warn("Input count matrix is integer, but counts=false.")
    else
        counts && @warn("Input count matrix is not integer, but counts=true.")
    end
    @pkgfun(myR, rlibrary, callr, j2r, r2j)
    rlibrary("GSVA")
    if pn>1
        rlibrary("snow")
    end
    M=r2j(callr("gsva", j2r(count, rownames=geneid), callr("list", gene_list1, gene_lists...), parallel!sz=pn, kcdf=counts ? "Poisson" : "Gaussian", kwargs...), NA=NaN, keepvec=true)
    if isempty(gene_lists)
        vec(M)
    else
        tuple(eachrow(M, as_vec=true)...)
    end
end
function GSVA(count::Matrix, isfg1::AbstractVector{Bool}, isfgs::AbstractVector{Bool}...; wargs...)
    geneid=string.(1:length(isfg1))
    gene_lists=map(x->geneid[x], [isfg1, isfgs...])
    GSVA(count, geneid, gene_lists...; wargs...)
end
#}}

#{{ GEOdownload
#[[ GEOdownload ]]
#smp_info = GEOdownload(c"GSMXXXXX, GSMXXXXX, ..."; label=GSM_IDs, path="GEO_downloaded", trymax=10, paral=T|F, directly=true)
#Downloaded files from GEO, convert to fastq file, and merge them if needed. Multi-threads is supported.
#The finial files will be $path/$(label)_1.fastq.gz. Besides, under $path, SRA_file/ is the sra files; fastq/ is the fastq.gz files before merging; GEO_sample_info.jld is the sample info table the same as function output.
#If directly=true, function will download and convert fastq file using fastq-dump directly. Otherwise, a .sra file will be download by wget first, then converted to fastq using fastq-dump.
#If FTP download is failed, function will retry until get trymax times.
#If path exists and some finished fastq files are in it, these files will not be re-downloaded.
#Entrezdirect (ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/) and SRA Toolkit (https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) need be installed.
#See also:
#Xiong Jieyi, 16 May 2017 > 7 Jun 2018 > 30 Oct 2019

export GEOdownload
function GEOdownload(GSMs::Vector{String}; label::Vector{String}=GSMs, path::String="GEO_downloaded", trymax::Int=10, paral::Bool=nprocs()>1, directly::Bool=true)
    
    if !(isuniqr(GSMs) && isuniqr(label))
        error("GSM ID or labels is not unique.")
    end
    
    asmp=smp=tb(GSM=GSMs, label=label)
    if isdir(path)
        println("Path $path has already exists. Check unfinished file...")
        l=cd(path) do
            map(x->isfile("$(x)_1.fastq.gz") ,label)
        end
        println(f"$1 finished fastq file(s) has(ve) found. Skip them."(sum(l)))
        smp=rec(smp, .!l)
    else
        mkdir(path)
    end
    
    println("Fetching sample info")
    smp["SRX"],smp["SRRs"]=rowfun(smp) do T
        println(T["GSM"])
        S=sfc"esearch -db sra -query \"$1\"| efetch -format docsum"(T["GSM"])|>readchomp
        isempty(S) && error("GEO ID is not found: "*T["GSM"])
        O=ParseStr(S)
        SRX=next(O,"<Experiment acc=\"",("\"",),(:no,"SRX"))[1] #SRX is no more used.
        # SRP=next(O,"<Study acc=\"",("\"",),(:no,"SRP"))[1]
        SRRs=String[]
        while next(O, (:test,"<Run acc=\""))[1]
            push!(SRRs, next(O,("\"",))[1])
        end
        (SRX,join(SRRs,','))
    end
    
    cd(path) do
        # jsave("GEO_sample_info", asmp)
        writetb("GEO_sample_info.tsv", smp)
        smp["isPairedEnd"]=falses(rnum(smp))
        # for (fi,(lb,SRX,t)) in enumerate(eachrow((smp["label"],smp["SRX"],smp["SRRs"])))
        mapfun=paral ? pmapwd : map
        mapfun(1:rnum(smp), smp["label"], smp["SRX"], smp["SRRs"]) do fi, lb, SRX, t
            SRRs=split(t,',')
            println(lb*" -- Downloading")
            multiSRR=length(SRRs)>1
            if directly
                for SRR in SRRs
                    for tryi=1:trymax
                        try
                            sf"fastq-dump --outdir fastq --gzip --split-files $1"(SRR)
                            break
                        catch err
                            if tryi==trymax
                                rethrow(err)
                            else
                                println(f"Download failed $1 times. Wait $2 second and try again...."(tryi, tryi*5))
                                sleep(tryi*5)
                            end
                        end
                    end
                end
            else
                for SRR in SRRs
                    for tryi=1:trymax
                        try
                            # sf"wget -c -P SRA_file ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/$1/$2/$3/$3.sra"(SRX[1:6],SRX,SRR)
                            sf"wget -c -P SRA_file ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/$1/$2/$2.sra"(SRR[1:6],SRR)
                            break
                        catch err
                            if tryi==trymax
                                rethrow(err)
                            else
                                println(f"Download failed $1 times. Wait $2 second and try again...."(tryi, tryi*5))
                                sleep(tryi*5)
                            end
                        end
                    end
                    println("Converting to Fastq")
                    # sf"sra2fastqgz SRA_file/$1.sra"(SRR) #fastq file is at fastq/
                    sf"fastq-dump --outdir fastq --gzip --skip-technical  --readids --dumpbase --split-files --clip SRA_file/$1.sra"(SRR)
                end
            end
            if multiSRR
                println("Merging")
                s1=f"fastq/$(1)_1.fastq.gz".(SRRs)
                sf"zcat $1|gzip -c > $(2)_1.fastq.gz"(join(s1, ' '), lb)
                if isfile(f"fastq/$(1)_2.fastq.gz"(SRRs[1]))
                    smp["isPairedEnd"][fi]=true
                    s2=f"fastq/$(1)_2.fastq.gz".(SRRs)
                    sf"zcat $1|gzip -c > $(2)_2.fastq.gz"(join(s2, ' '), lb)
                    rm.(s2)
                end
                rm.(s1)
            else
                if isfile(f"fastq/$(1)_2.fastq.gz"(SRRs[1]))
                    smp["isPairedEnd"][fi]=true
                    sf"mv fastq/$(1)_2.fastq.gz $(2)_2.fastq.gz"(SRRs[1], lb)
                end
                sf"mv fastq/$(1)_1.fastq.gz $(2)_1.fastq.gz"(SRRs[1], lb)
            end
            if !directly
                for SRR in SRRs
                    rm("SRA_file/$SRR.sra")
                end
            end
            println("Done.")
        end
        directly || rm("SRA_file")
        rm("fastq")
    end
    jldsave("GEO_sample_info", smp)
    println("Finished.")
    smp
end
#}}

#{{ survdiff, survcoxph
#[[ survdiff survcoxph ]]
# pvalue = survdiff(surv_time, status, group; vb=0)
# pvalue, coef = survcoxph(surv_time, status, mainFactor[, otherFactor1, otherfactor2, ...]; vb=0)
# Test survival differences among catalogus or to a numeric variable.
# `status` could be a vector of 0/1 or 1/2 to represent censored vs “dead”.
# verbose: Whether show model detail on screen.
# When mainFactor is numeric, the coefficients in a Cox regression relate to hazard; a positive coefficient indicates a worse prognosis and a negative coefficient indicates a protective effect of the variable with which it is associated.
# When mainFactor is factorial, the coeff is a Pair vector for each factor => coefficient.
# If calculation is failed, function return NaN.
# It needs myR package to be preloaded, and survival package has installed (see http://bioconnector.org/workshops/r-survival.html).
# See also: survivalplot, survfit
#Xiong Jieyi, 15 Jun 2018

export survdiff, survcoxph
function survdiff(time::AbstractVector{Tt}, status::AbstractVector{Ts}, grp::AbstractVector; vb::Int=0) where {Tt<:Real, Ts<:Integer}
    0 < length(time) == length(status) == length(grp) || error("Inputs have different lengths or zero lengths.")
    @pkgfun(myR, rlibrary, callr, r2j)
    rlibrary("survival")
    rfun=k"""function(x, y, z, vb){
        sd<-survdiff(Surv(x, y)~z)
        if(vb==1 || vb==3){show(sd)}
        if(vb==2 || vb==3){show(summary(sd))}
        pchisq(sd$chisq, length(sd$n)-1, lower.tail=FALSE)
        }    
        """
    r2j(callr(rfun, time, status, grp, vb), keepvec=false, NA=NaN)
end
function survcoxph(time::AbstractVector{Tt}, status::AbstractVector{Ts}, X::AbstractVector, Fs::AbstractVector...; vb::Int=0) where {Tt<:Real, Ts<:Integer}
    0 < length(time) == length(status) == length(X) || error("Inputs have different lengths or zero lengths.")
    @pkgfun(myR, rlibrary, callr, r2j)
    rlibrary("survival")
    if isempty(Fs)
        f1=""
        f2=""
        f3="1"
    else
        f1=join(f", F$1".(1:length(Fs)))
        f2=join(f"F$1+".(1:length(Fs)))
        f3=join(f"F$1".(1:length(Fs)), "+")
    end
    rfun="""function(x, y, z, vb$(f1)){
        m1<-coxph(Surv(x, y)~$(f2)z)
        m0<-coxph(Surv(x, y)~$(f3))
        anv<-anova(m1, m0)
        if(vb==1 || vb==3){show(m1); show(anv)}
        if(vb==2 || vb==3){show(summary(m1)); show(summary(anv))}
        list(tail(anv[["P(>|Chi|)"]], n=1), m1\$coefficients, names(m1\$coefficients))
        }    
        """
    t=callr(rfun, time, status, X, vb, Fs...)
    pval=r2j(t[1], keepvec=false, NA=NaN)
    coef, coefname =( r2j(t[2]), r2j(t[3]) )
    l=startswith.(coefname, "z")
    if sum(l)==1 && coefname[l][1]=="z"
        (pval, coef[l][1])
    else
        (pval, replace.(coefname[l], [r"^z"=>""]).=>coef[l])
    end
end
#}}

#{{ survfit
#[[ survfit ]]
# nothing = survfit(time, status[, group])
# Show mean and median survival time on screen. No output value.
# Equal to print(survfit(Surv(time, status) ~ group), print.rmean=TRUE) in R.
# See also: survdiff, survcoxph
# Jieyi Xiong, 2 Apr 2019

export survfit
function survfit(time::AbstractVector{Tt}, status::AbstractVector{Ts}, grp::AbstractVector) where {Tt<:Real, Ts<:Integer}
    @pkgfun(myR, rlibrary, callr)
    rlibrary("survival")
    rfun="function(time, status, group){print(survfit(Surv(time, status) ~ group), print.rmean=TRUE)}"
    callr(rfun, time, status, grp)
    nothing         
end
function survfit(time::AbstractVector{Tt}, status::AbstractVector{Ts}) where {Tt<:Real, Ts<:Integer}
    @pkgfun(myR, rlibrary, callr)
    rlibrary("survival")
    rfun="function(time, status){print(survfit(Surv(time, status) ~ 1), print.rmean=TRUE)}"
    callr(rfun, time, status)
    nothing         
end
#}}

#{{ survivalplot
#[[ survivalplot ]]
# handles = survivalplot(sfit; legend_inp=nothing|inp([labels;] ...) | legend_kw=(...))
# Or ...   = survivalplot(time, status[, group]; vb=0, ...) #Fit the model automatically.
# Draw the survival curves fitted by R package "survival".
# For R:survival package, see http://bioconnector.org/workshops/r-survival.html
# `status` could be a vector of 0/1 or 1/2 to represent censored vs “dead”.
# When legend_inp is nothing, no legend will be show.
# See also: survdiff, survcoxph
#Xiong Jieyi, 1 Jul 2017 > 3 Jan 2018 >15 Jun 2018

export survivalplot
function survivalplot(sfit; legend_inp::Union{Tuple, Void}=inp(), legend_kw=legend_inp[2])
    # myr=importpkg(:myR, preloaded=true)
    @pkgfun(myR, rlibrary, callr, callr1, relem, relemj, rnames)
    plt=importpkg(:PyPlot, preloaded=true)
    xx=relemj(sfit,"time", keepvec=true)
    yy=relemj(sfit,"surv", keepvec=true)
    ss=relemj(sfit,"n.censor", keepvec=true).>0
    num=relemj(sfit,"n", keepvec=true)
    
    is_strata=callr1("function(x,y) y %in% names(x)", sfit, "strata")
    if is_strata
        strata=[0;cumsum(relemj(sfit,"strata", keepvec=true))]
        labels=map(x->replace(x, r"^(z|group)\=" => ""), rnames(relem(sfit, "strata")))
        labels=f"$1 (N=$2)".(labels,num)
    else
        strata=[0,length(xx)]
    end
    
    hds=[]
    for i=1:length(strata)-1
        l=strata[i]+1:strata[i+1]
        x=[0;xx[l]]
        y=[1;yy[l]]
        s=[false;ss[l]]
        h=plt.step(x,y,where="post")[1] #Here draw step curves.
        plt.plot(x[s],y[s],"w+",mec=h.get_color()) #Here draw cencor dots.
        push!(hds,h)
    end
    plt.xlabel("Time")
    plt.ylabel("Survival probability")
    if is_strata && legend_inp!=nothing
        if !isempty(legend_inp[1])
            labels=legend_inp[1]
        end
        plt.legend(hds, labels; legend_kw...)
    end
    
    hds
end

function survivalplot(time::Vector{Tt}, status::Vector{Ts}, grp::Vector=[]; legend_inp::Union{Tuple, Void}=inp(), legend_kw=legend_inp[2], vb::Int=0) where {Tt<:Real, Ts<:Integer}
    0 < length(time) == length(status) == length(grp) || error("Inputs have different lengths or zero lengths.")
    @pkgfun(myR, rlibrary, callr, callrw, eR, rshow)
    rlibrary("survival")
    if isempty(grp)
        data=callr("data.frame",time=time,status=status)
        sfit=callrw("survfit",eR("Surv(time, status)~1"),data)
    else
        data=callr("data.frame", group=grp,time=time,status=status)
        sfit=callrw("survfit",eR("Surv(time, status)~group"),data)
    end
    rshow(sfit, vb=vb)
    survivalplot(sfit; legend_inp=legend_inp, legend_kw=legend_kw)
end
#}}

#{{ liftover
#[[ liftover ]]
# chr|chrno, pos[, strand] = liftover(chr|chrno, pos, chainfile; strand=..., liftover_param="liftOver")
# Or ... = liftover((chr|chrno, pos[, strand]), chainfile; ...)
# Warped liftOver in Julia. The output always has the same row number as the input. The unliftoverable coordinates will be filled with empty values. The liftOver program and other parameters can be assigned in liftover_param.
# See also:
# Jieyi Xiong, 25 Apr 2018

export liftover
function liftover(chr::Vector{T1}, pos::AbstractMatrix{T2}, chainfile::AbstractString; strand::AbstractVector{T3}=Char[], liftover_param::AbstractString="liftOver") where T1<:AbstractString where T2<:Number where T3<:Char
    fn=tempname()
    fnout=tempname()
    open(fn,"w") do fid
        for i=1:length(chr)
            println(fid, join((chr[i],pos[i,1]-1,pos[i,2],i,255,isempty(strand) ? '+' : strand[i]), '\t'))
        end
    end
    liftover_param=isempty(liftover_param) ? "" : " "*liftover_param
    sf"$1 $2 $3 $4 /dev/null"(liftover_param, fn, chainfile, fnout)
    T=readtb(fnout, head=c"chr, N::pos, <, N::NO, X::, C::strand")
    rm(fn)
    rm(fnout)
    T["pos"][:,1]+=1
    fT=dtshift(1:length(chr), T["NO"], T, nothing)
    if isempty(strand)
        d"chr, pos"fT
    else
        d"chr, pos, strand"fT
    end
end
function liftover(chr::Vector{T1}, arg...; warg...) where T1<:Number
    C=liftover(no2chr(chr), arg...; warg...)
    tuple(chr2no(C[1]), C[2:end]...)
end
function liftover(chrpos::Tuple, arg...; warg...)
    if length(chrpos)==2
        liftover(chrpos..., arg...; warg...)
    else
        length(chrpos)==3 || error("Invalid tuple input.")
        liftover(chrpos[1], chrpos[2], arg...; strand=chrpos[3], warg...)
    end
end
#}}

# Code needs be further developed
#=
 #{{ segcloest
# #[[ segcloest ]]
#15 May 2017: Problem for this code: It's beheave is hard to be understand. The repeated segments is not be considered. It costs more coding time to use this function than do it manually by genomemap.
#
# Host_idx, Guest_idx, distance=segcloest(pos|(chr, pos); grp=group_ID, lim=(-XX, XX)|farest=XX, strand='+'s, uniqhost=false, dup=false) #A/B is the item with smaller/larger group ID in each pair.
# ...             =segcloest(Host_pos, Guest_pos; lim=(-XX, XX)|farest=XX, strand=CharVector_of_host, uniqhost=true)
#Choose the cloest pairs in any different groups. Their index are Ai and Bi. grp[Ai]!=grp[Bi]. The Ai is the group with smaller ID (host-group) and Bi is the group with larger ID (guest-group). The strand input are only for the host-group. If the guest group is in the downstream/upstream of host-group, distance is positive/negative. If host group are fully or partly overlapped with guest group, distance is 0. Only lim[1]<=distance<=lim[2] will be reported. If farest is assigned, lim=(-farest, farest). Note that if two segments are joined but not overlapped, the distance is 1/-1 rather than 0.
#In the pos-grp model, although all the strands are inputed, only the strand of host group is considered. In the pos1-pos2 model, only strand of host group is required.
#If grp is missing or be [], function will consider every row of pos are different, and the pos with smaller number are always be the host segment in each output pairs, irrelavent with the order of input pos.
#If uniqhost=true, function will try to uniquefy Ai at last and only report the pairs with the minimum abs(distance) for each Ai, unless there is more than one minimum value. In pos-grp module, this argument only work when grp is given. Note that the default values in pair mode and pos-grp mode are different.
#If dup is true, function will return ([Ai,Bi],[Bi,Ai],[D,+/-D]) rather than (Ai,Bi,D).
#See also: genomemap
#Xiong Jieyi, Sep 3, 2015

export segcloest
#<v0.6# @compat function segcloest{T<:Real}(pos::Matrix{T}; grp::Group=[], lim::Tuple{Real,Real}=(-Inf,Inf), farest::Real=NaN, strand::Vector{Char}=fill('+',size(pos,1)), uniqhost::Bool=false, dup::Bool=false)
@compat function segcloest(pos::Matrix{T}; grp::Group=[], lim::Tuple{Real,Real}=(-Inf,Inf), farest::Real=NaN, strand::Vector{Char}=fill('+',size(pos,1)), uniqhost::Bool=false, dup::Bool=false) where {T<:Real}
    if !isnan(farest)
        lim=(-farest,farest)
    end
    
    spos,idx=sortr(pos)
    sS=strand[idx]
    if isempty(grp)
        sG=1:size(pos,1)
        uniqhost=false #Not need to unique host.
    else
        sG=getrow(grp,idx)
    end

    rlt=rowpile()
    #Remove overlapped items
    len=size(spos,1)
    l=trues(len)
    op=1
    np=2
    while np<=len
        if spos[op,2]>spos[np,2]
            if getrow(sG,op)!=getrow(sG,np)
                if getrow(sG,op)>getrow(sG,np)
                    cAi=np
                    cBi=op
                else
                    cAi=op
                    cBi=np
                end
                addrow!(rlt,idx[cAi],idx[cBi],0)
            end
            l[np]=false
            np+=1
        else
            op=np
            np+=1
        end
    end
    sG=sG[l]
    sS=sS[l]
    spos=spos[l,:]
    idx=idx[l]
    
    for i=2:size(spos,1)
        Gp=getrow(sG,i-1)
        Gi=getrow(sG,i)
        if Gp==Gi
            continue
        end
        if spos[i-1,2]>=spos[i,1]
            D=0
        else
            D=spos[i,1]-spos[i-1,2]
        end
        if Gp<Gi
            cAi=i-1
            cBi=i
            rev=sS[i-1]=='-'
        else
            cAi=i
            cBi=i-1
            rev=sS[i]!='-'
        end
        if rev
            D=-D
        end
        if lim[1]<=D<=lim[2]
            addrow!(rlt,idx[cAi],idx[cBi],D)
        end
    end
    if isvirgin(rlt)
        return (Int[],Int[],T[])
    end
    Ai,Bi,D=value(rlt)
    if dup
        @assert(all((strand.=='+')|(strand.=='-')),"When dup is using, strand can only be + or -.")
        nD=D
        for i=1:length(D)
            if strand[Ai[i]]==strand[Bi[i]]
                nD[i]=-nD[i]
            end
        end
        (Ai,Bi,D)=([Ai,Bi],[Bi,Ai],[D,nD])
    end
    if uniqhost
        l=grpfunexp(x->x.==minimum(x),Ai,abs(D))
        Ai=Ai[l]
        Bi=Bi[l]
        D=D[l]
    # else
    #     if any((pos[Ai,1].==pos[Bi,1]) & (pos[Ai,2].==pos[Bi,2]))
    #         warn("Some host segment and guest segment are identical. In uniqhost=false, one of two pairs in these identical segments will lost in the report due to the algorithm limitation.")
    #     end
    end
    Ai,Bi,D
end
@compat function segcloest(chrpos::Tuple{Group,Matrix};grp=1:rownum(chrpos),strand::Vector{Char}=fill('+',rownum(chrpos)),wargs...)
    rlt,=grpfunwith(chrpos[1],chrpos[2], grp, strand, 1:rownum(chrpos), multi_output=true) do pos, cgrp, cstrand, l
        Ai,Bi,D=segcloest(pos;grp=cgrp,strand=cstrand,wargs...)
        (l[Ai],l[Bi],D)
    end
    rlt
end
@compat segcloest(chrpos::Tuple{Group,Matrix,Vector{Char}};wargs...)=segcloest(chrpos[1:2];strand=chrpos[3],wargs...)
function segcloest(pos1,pos2;grp=nothing,strand::Vector{Char}=fill('+',rownum(pos1)),uniqhost::Bool=true,dup::Bool=false,wargs...)
    @assert(grp==nothing,"If grp is given, function only accept one pos argument.")
    @assert(!dup,"dup=true is only supported in group mode so far.")
    l=[1:rownum(pos1);1:rownum(pos2)]
    grp,(pos,astrand)=vcatr_with(1:2,(pos1,strand),(pos2,fill('.',rownum(pos2))))
    Ai,Bi,D=segcloest(pos;grp=grp,strand=astrand,uniqhost=uniqhost,dup=false,wargs...)
    (l[Ai],l[Bi],D)
end
# function segcloestself(pos;grp::Void=nothing,uniqhost::Bool=false,dup::Bool=false,wargs...)
#     @assert(uniqhost==false,"In segcloestself, uniqhost can only be false.")
#     if isa(pos,Tuple) && length(pos)==3
#         strand=pos[3]
#         wargs=push!(wargs,(:strand,pos[3]))
#         pos=tuple(pos[1],pos[2])
#     end
#     Ai,Bi,D=segcloest(pos,pos;uniqhost=false,strand=strand,wargs...)
#     l=Ai.!=Bi
#     Ai=Ai[l]
#     Bi=Bi[l]
#     D=D[l]
#     if dup
#         nD=D
#         for i=1:length(D)
#             if strand[Ai[i]]==strand[Bi[i]]
#                 nD[i]=-nD[i]
#             end
#         end
#         ([Ai,Bi],[Bi,Ai],[D,nD])
#     else
#         (Ai,Bi,D)
#     end
# end
#}}
=#
