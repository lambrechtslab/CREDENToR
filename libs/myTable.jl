#{{ tb ds

#[[ tb tb! ds ds! ]]
#Dict{String,Any}(...)=tb|ds([Dict, ]"key1"=>val1, "key2"=>val2; key3=val3, key4=val4, ...)
# ... = tb|ds([Dict, ]c"key1,key2,..." => (val1, val2, ...), ...; ...)
# ... = tb|ds([Dict, ]("key1", ("key2", "key3")) => (value1, (value2, value3)), ...; ...)
# ... = tb|ds(A, B; ...) = tb|ds(A => B; ...) #Only support 2-3 input.
# ... = tb!|ds!(Dict, ...; ...) #The inputed Dict will be modified.
#Create a table-style Dict (Dict{String,Any}()). In tb, the rownum will be checked but not in ds.
#See also: dictfun, tb2DataFrame, pandasDF, rec, fd*
#Xiong Jieyi, February 9, 2015 >May 27, 2015>Jun 3, 2015>Nov 4, 2015>Nov 4, 2015>23 Jul 2017>18 Dec 2018

export tb, ds, tb!, ds!

function __deepzip__(P::Pair)
    if isa(P[1], AbstractString)
        return [convert(Pair{String, Any}, P)]
    end
    isa(P[1], Union{AbstractVector, Tuple}) || error("Key of a Table should always be String.")
    length(P[1])==length(P[2]) || error("Two inputs cannot be paired due to different lengths or strcutures.")
    O=Pair{String, Any}[]
    for (a, b) in zip(P[1], P[2])
        append!(O, __deepzip__(a=>b))
    end
    O
end

function tb!(D::Dict{T, Any}, args::Pair...; wargs...) where T<:AbstractString
    if isempty(D)
        len=-1
    else
        len=rnum(D)
    end
    for carg in args
        for (k, v) in __deepzip__(carg)
            if len==-1
                len=rnum(v)
            else
                len==rnum(v) || error("Inputs $k doesn't have the same row number $len.")
            end
            D[k]=v
        end
    end
    for (k, v) in wargs
        if len==-1
            len=rnum(v)
        else
            len==rnum(v) || error("Inputs $k doesn't have the same row number $len.")
        end
        D[string(k)]=v
    end
    D
end

# function tb!(D::Dict=Dict{ASCIIString,Any}(); wargs...)
#     if isempty(D)
#         len=-1
#     else
#         # @assert(istable(D),"First input is not a table.")
#         len=rownumsf(D)
#     end
#     for (k,v) in wargs
#         if len==-1
#             len=rownum(v)
#         else
#             @assert(len==rownum(v),"Inputs doesn't have the same row number.")
#         end
#         D[string(k)]=v
#     end
#     D
# end

# tb!(T::Dict, fds::Union{Vector,Tuple}, X::Union{Tuple,Vector}; warg...)=tb!(T, fds => X; warg...)
# function tb!(T::Dict,fds::Union{Vector,Tuple},X::Union{Tuple,Vector}, arg...; warg...)
#     @assert(length(fds)==length(X),"Key and value number is not identical.")
#     if length(X)>1
#         len=rownum(X[1])
#         for i=2:length(X)
#             @assert(rownum(X[i])==len,"Inputs doesn't have the same row number.")
#         end
#     end
#     for (fd, val) in zip(fds, X)
#         T[fd]=val
#     end
#     if isempty(arg) && isempty(warg)
#         T
#     else
#         tb!(T, arg...; warg...)
#     end
# end

# tb()=Dict{String, Any}()
tb(fds::Union{AbstractVector,Tuple}, X::Union{Tuple, AbstractVector}; warg...)=tb!(tb(), fds => X; warg...)
tb(K::Pair, args::Pair...; wargs...)=tb!(tb(), K, args...; wargs...)
tb(; wargs...)=tb!(Dict{String, Any}(); wargs...)
tb(D::Dict, args...; wargs...)=tb!(deepcopy(D), args...; wargs...)

function ds!(D::Dict{T, Any}, args::Pair...; wargs...) where T<:AbstractString
    for carg in args
        for (k, v) in __deepzip__(carg)
            D[k]=v
        end
    end
    for (k, v) in wargs
        D[string(k)]=v
    end
    D
end
# function ds!(D::Dict=Dict{ASCIIString,Any}();wargs...)
#     for (k, v) in wargs
#         D[string(k)]=v
#     end
#     D
# end
# ds()=Dict{String, Any}()
ds(K::Pair, args...; wargs...)=ds!(ds(), K, args...; wargs...)
ds(; wargs...)=ds!(Dict{String, Any}(); wargs...)
ds(D::Dict,args...;wargs...)=ds!(deepcopy(D),args...;wargs...)
# ds!(T::Dict, fds::Union{Vector,Tuple}, X::Union{Tuple,Vector}; wargs...)=ds!(T, fds => X; wargs...)
# function ds!(T::Dict,fds::Union{Vector,Tuple},X::Union{Tuple,Vector}, arg...; warg...)
#     for (fd, val) in zip(fds, X)
#         T[fd]=val
#     end
#     if isempty(arg) && isempty(warg)
#         T
#     else
#         ds!(T, arg...; warg...)
#     end
# end
ds(fds::Union{AbstractVector,Tuple}, X::Union{Tuple, AbstractVector}; warg...)=ds!(ds(), fds => X; warg...)

#}}

#{{ dictfun dictfun!
#[[ dictfun dictfun! ]]
#outD =  dictfun(func, D::Dict...; key=keys(D1)|:shared, copy=false)
#  D1 = dictfun!(func, D::Dict...; key=keys(D1)|:shared, replace=true)
#  DF = dataframefun(func, DF::DataFrame...)
#Run function on each value of Dict under the given keys. Note that both dictfun and dictfun! will change the value of input Dict if the given function do that.
#The only difference of these two function is dictfun collect the output of function as a new dict and return it, while dictfun! just simplely return the first argument.
#If copy=true, the missing D-fields in key will be deepcopy to outD.
#If replace=true, D1[key] will be assigned to function output.
#If key=:shared, only the shared keys among all inputed Dict are calculated. Note that Dict number should be at least two.
#Change in 16 May 2019: changed the default value of dectfun!( replace=false) to true.
#Xiong Jieyi, 13 May 2014>11 Sep 2014>October 14, 2014>Nov 25, 2015>Jul 14, 2016>28 Jan 2019>16 May 2019

function dictfun(f::Function,D::Dict{Tk,}; key=keys(D), copy::Bool=false) where {Tk}
    if isa(key, Symbol)
        if key==:shared
            error("key=:shared is not supported when only one Dict inputted.")
        else
            error("Invalid key=.... It should be a container or :shared.")
        end
    end
    O=Dict{Tk,Any}()
    for k in key
        O[k]=f(D[k])
    end
    if copy
        for k in keys(D)
            if !in(k,key)
                O[k]=deepcopy(D[k])
            end
        end
    end
    return O
end
function dictfun(f::Function,args::Dict...;key=keys(args[1]),copy::Bool=false)
    if isa(key, Symbol)
        if key==:shared
            key=intersect(keys.(args)...)
        else
            error("Invalid key=.... It should be a container or :shared.")
        end
    end
    O=Dict{keytype(args[1]),Any}()
    for k in key
        O[k]=f(map(x->x[k],args)...)
    end
    if copy
        for k in keys(args[1])
            if !in(k,key)
                O[k]=deepcopy(args[1][k])
            end
        end
    end
    return O
end
function dataframefun(f::Function, args::AbstractDataFrame...)
    O=DataFrame()
    nms=names(args[1])
    for k in nms
        O[k]=f(map(x->x[k],args)...)
    end
    return O
end

function dictfun!(f::Function, D::Dict; key=keys(D), replace::Bool=true)
    if isa(key, Symbol)
        if key==:shared
            error("key=:shared is not supported when only one Dict inputted.")
        else
            error("Invalid key=.... It should be a container or :shared.")
        end
    end
    if replace
        for k in key
            D[k]=f(D[k])
        end
    else
        for k in key
            f(D[k])
        end
    end
    return D
end
function dictfun!(f::Function,args::Dict...;key=keys(args[1]),replace::Bool=true)
    if isa(key, Symbol)
        if key==:shared
            key=intersect(keys.(args)...)
        else
            error("Invalid key=.... It should be a container or :shared.")
        end
    end
    if replace
        for k in key
            args[1][k]=f(map(x->x[k],args)...)
        end
    else
        for k in key
            #If below error occurs, it may because of the confliction of Rif.keys against Base.keys.
            # ERROR: `keys` has no method matching keys(::Dict{ASCIIString,Array{Int64,1}})
            # you may have intended to import Base.keys
            f(map(x->x[k],args)...)
        end
    end
    return args[1]
end
export dictfun, dictfun!, dataframefun
#}}

#{{ istable
#[[ istable ]]
#istable(Dict)
#If the Dict is a table, i.e., have the same row number for each elements, return true. Otherwise, return false.
#See also: vcatr, rec, fd, fd_i, ft, tbx, tbuniq, tbuniq_i, dictconc, dictconc!, tbcomp
#Xiong Jieyi, 13 May 2014

function istable(D::Dict)
    
    rn=-1;
    for (k,v) in D
        if rn==-1
            rn=rownum(v)
        elseif rn!=rownum(v)
            return false
        end
    end
    return true
end
export istable
#}}

#{{ fd fd_i fd! fd_i!
#[[ fd fd_i fd! fd_i! ]]
#sub_Dict = fd|fd!(Dict, fields::AnyIterateable|Regex|AbstractString)
#Or ...= fd_i|fd_i!(Dict, fields::AnyIterateable|Regex|AbstractString)
#Or fd|fd_i|fd_i!(Dict, X1, X2...) #Union or substract all given field|Regex. [*]
#Get structure with subfields. fd_i only reture the Dict with the feilds not included or matched the Regex.
#The ! function will also change the input Dict.
#The difference between delete!(D::Dict,f::AbstractString) and fd_i!(D::Dict,f::AbstractString): fd_i! will check if given field exists in the Dict, but delete! will not.
#[*] Try to avoid fd(D,c"a",c"b",...) but use fd(D, c"a,b,...") instead for speed reason. You can use fd(D,c"field", r"Regex1", r"Regex2"...).
#See also: sk, rename, dictfun
#Xiong Jieyi, 3 Sep 2014 > February 10, 2015 >Oct 9, 2015 >Aug 12, 2016

import Base.fd #fd is used in Base but for other useage.
export fd, fd_i, fd!, fd_i!
function fd(D::Dict{Tk,Tv}, F) where {Tk,Tv}
    O=Dict{Tk,Tv}()
    for f in F
        O[f]=D[f]
    end
    return O
end
fd(D::Dict,F::AbstractString)=fd(D,[F])
function fd(D::Dict, R::Regex)
    # filter((k,v)->ismatch(R,k),D)
    filter(((k, _),)->occursin(R, k), D)
end
function fd_i(D::Dict, F)
    K=collect(keys(D))
    fd(D,setdiff(K,F))
end
fd_i(D::Dict,F::AbstractString)=fd_i(D,[F])
function fd_i(D::Dict, R::Regex)
    filter(((k, _),)->!occursin(R, k), D)
end
function fd_i!(D::Dict, F::Vector{Tf}) where {Tf<:AbstractString}
    for f in F
        @assert(haskey(D,f),"Field $f is not exist.")
        delete!(D,f)
    end
    D
end
function fd!(D::Dict, F)
    K=collect(keys(D))
    for f in setdiff(K,F)
        delete!(D,f)
    end
    D
end
fd!(D::Dict,F::AbstractString)=fd(D,[F])
function fd!(D::Dict, R::Regex)
    filter!(((k, _),)->occursin(R, k), D)
end
function fd_i!(D::Dict, R::Regex)
    filter!(((k, _),)->!occursin(R, k), D)
end
function fd_i!(D::Dict,f::AbstractString)
    haskey(D,f) || error("Field $f is not exist.")
    delete!(D,f)
end
fd_i(D::Dict,X,P...)=fd_i(fd_i(D,X),P...)
fd_i!(D::Dict,X,P...)=fd_i!(fd_i!(D,X),P...)
fd(D::Dict,P...)=dictconc(map(x->fd(D,x),P)...)
#}}

#{{ @d_str

#[[ @d_str  ]]
#d"field1[, field2, ...]"dict or d"""..."""dict
#return (dict["field1"], dict["field2"],...)
#field will be strip.
#d"field"dict equal to dict["field"]
#For the simplicity to form a tuple with the table fields.
#See also: @c_str, @dd_str, dict2tuple, @with
#Xiong Jieyi, 9 Sep 2014>January 6, 2015>Oct 19, 2015>May 13, 2016

macro d_str(S,P...)
    fields=map(strip,split(S,','))
    # for i=1:length(fields)
    #     fields[i]=replace(fields[i],r"^[\ \t]*\n","")
    # end
    cmd=if isempty(P)
        "D->("*join(["""D["$f"]""" for f in fields],',')*")"
    else
        dct=onlyone(P)
        join(["""$dct["$f"]""" for f in fields],',')
    end
    esc(Meta.parse(cmd))
end
export @d_str

# if VERSION<v"0.4"
#     macro d_mstr(fields,dct)
#         fields=map(lstrip,split(fields,','))
#         esc(parse(join(["""$dct["$f"]""" for f in fields],',')))
#     end
#     export @d_mstr
# end
#}}

#{{ @dd_str
#[[ @dd_str ]]
# SubDict = dd" field1[ => rename1], field2[ => prefix_# ],..."Dict
# Or ...  = dd"..."(Dict)
# fd and rename! Dict. If # is in new field name but not in the old one, it will be replaced by the old field name. All the blanks at the flanks of names will be stripped.
#See also: @d_str, fd*, rename!, @with
#Xiong Jieyi, May 13, 2016

export @dd_str
macro dd_str(S,P...)
    fds1=AbstractString[]
    fds2=AbstractString[]
    fds3=AbstractString[]
    for C in split(S,',')
        C=strip(C)
        tC=split(C,"=>")
        if length(tC)>1
            tC[1]=strip(tC[1])
            tC[2]=strip(tC[2])
            push!(fds1,tC[1])
            push!(fds2,tC[1])
            if !in('#',tC[1])
                tC[2]=replace(tC[2], "#" => tC[1])
            end
            push!(fds3,tC[2])
        else
            push!(fds1,C)
        end
    end
    cmd=if isempty(P)
        if isempty(fds2)
            "D->fd(D,$fds1)"
        else
            "D->rename!(fd(D,$fds1),$fds2,$fds3)"
        end
    else
        D=onlyone(P)
        if isempty(fds2)
            "fd($D,$fds1)"
        else
            "rename!(fd($D,$fds1),$fds2,$fds3)"
        end
    end
    esc(Meta.parse(cmd))
end
#}}

#{{ @with
#[[ @with ]]
#@with Dict any-cmd-with-$field-or-$"field"
#A lasy way to use Dict["field"]. In the second parameter, any $field or $"field" will be replaced as Dict["field"].
#See also: @d_str, @dd_str
#Xiong Jieyi, 27 Oct 2017

macro with(T, C)
    function recur(C)
        if C.head==:$ && length(C.args)==1 && (isa(C.args[1], Symbol) || isa(C.args[1], String))
            C.head=:ref
            C.args=Any[T,string(C.args[1])]
        else
            for (i, v) in enumerate(C.args)
                if isa(v, Expr)
                    recur(v)
                end
            end
        end
    end
    isa(C, Expr) && recur(C)
    esc(C)
end
export @with
#}}

#{{ rec rec_i
#[[ rec rec_i @rec @rec_i ]]
#Table = rec(Table, index)
#Table = rec_i(Table, index)
# ...  = @rec[_i](Table, $field.>0)
#Get the rows of each element of Dict (rec), or remove the rows of Dict (rec_i)
#The difference between rec and getrow: rec will check if the input is a table firstly, but getrow will not.
#See also: getrow, fd, fd_i, ft, sk, dictconc, tbx, dict2tuple
#Xiong Jieyi, 8 Jul 2014>28 Sep 2014>October 16, 2014>14 Nov 2017

function rec(D::Dict,I::Union{Int,AbstractVector})
    istable(D) || error("Input is not a table.")
    getrow(D,I)
end
rec_i(D::Dict,I::BitVector)=rec(D, .!I)
rec_i(D::Dict,I::AbstractVector{Bool})=rec(D, .!I)
rec_i(D::Dict,I::AbstractVector)=rec(D,fill(true,rownum(D))|>x->(x[I].=false;x))
rec_i(D::Dict,i::Integer)=rec(D,fill(true,rownum(D))|>x->(x[i]=false;x))

macro rec(T, C)
    function recur(C)
        if C.head==:$ && length(C.args)==1 && (isa(C.args[1], Symbol) || isa(C.args[1], String))
            C.head=:ref
            C.args=Any[T,string(C.args[1])]
        else
            for (i, v) in enumerate(C.args)
                if isa(v, Expr)
                    recur(v)
                end
            end
        end
    end
    isa(C, Expr) && recur(C)
    esc(:(rec($T, $C)))
end

macro rec_i(T, C)
    function recur(C)
        if C.head==:$ && length(C.args)==1 && (isa(C.args[1], Symbol) || isa(C.args[1], String))
            C.head=:ref
            C.args=Any[T,string(C.args[1])]
        else
            for (i, v) in enumerate(C.args)
                if isa(v, Expr)
                    recur(v)
                end
            end
        end
    end
    isa(C, Expr) && recur(C)
    esc(:(rec_i($T, $C)))
end

export rec, rec_i, @rec, @rec_i
#}}

#{{ dictconc dictconc!
#[[ dictconc dictconc! ]]
#CombindDict = dictconc(A::Dict, B::Dict, ...; delconflict=false | aprefix="Pre-"|function, bprefix=...)
#          A = dictconc!(A::Dict, B::Dict, ...; delconflict=false )
#Combind all the fields of two dictionaries. If two field's values with the same fieldnames are inconsistent, a error will occur in default, or this field will be deleted (delconflict=true) or renamed (a/bprefix=...). aprefix/bprefix can also be a function, similar as rename!(function, ...). aprefix/bprefix can only be used when input Dict number = 2.
#In dictconc!, A will be added B's fields and be outputed, but B is ALWAYS intact (even when bprefix=... is assigned).
#See also: vcatr, vcatr_with, rec, fd, fd_i, ft, tbx, tbuniq, tbuniq_i, dict2tuple, tbcomp, rename!
#Xiong Jieyi, 16 May 2014 >May 22, 2015 >Oct 9, 2015>Apr 27, 2016>13 Jan 2018>3 May 2019

function dictconc!(O::Dict,B::Dict; delconflict::Bool=false,
                   aprefix::Union{AbstractString, Function}="",
                   bprefix::Union{AbstractString, Function}="")
    if !(isempty(aprefix) && isempty(bprefix))
        delconflict && error("rename_conflict=true conflicts with field renaming.")
        confkys=collect(filter(ky->!isequal(O[ky], B[ky]), intersect(keys(O), keys(B))))
        if !isempty(aprefix)
            O=rename!(isa(aprefix, AbstractString) ? x->aprefix*x : aprefix, O, confkys)
        end
        # if !isempty(bprefix)
        #     B=rename!(isa(bprefix, AbstractString) ? x->bprefix*x : bprefix, deepcopy(B), confkys)
        # end
    end
    
    for (k, v) in B
        if haskey(O, k)
            if !isequal(O[k], v)
                if delconflict
                    delete!(O, k)
                elseif !isempty(bprefix)
                    nk=isa(bprefix, AbstractString) ? bprefix*k : bprefix(k)
                    haskey(O, nk) && error("New field $k => $nk is still conflicted.")
                    O[nk]=B[k]
                else
                    error("The coexisted field \"$k\" is heterogeneous.")
                end
            end
        else
            O[k]=B[k]
        end
    end
    return O
end
dictconc!(O::Dict,B::Dict,P::Dict...; delconflict::Bool=false)=dictconc!(dictconc!(O::Dict, B::Dict; delconflict=delconflict), P...; delconflict=delconflict)
dictconc(A::Dict,P::Dict...; warg...)=dictconc!(deepcopy(A),P...; warg...)
export dictconc,dictconc!
#}}

#{{ tbuniq, tbuniq_i
#[[ tbuniq tbuniq_i ]]
# UniqTb=tbuniq(Tb)
# UniqTb=tbuniq(Tb, fields|nothing*[, "field"=>function, ...]; trim=false, check=false)
#  ...  =tbuniq(function, Tb, fields, "field"=>:, ...; ...) # ":" will be replaced by the first funciton input.
# UniqTb=tbuniq_i(Tb, ignore_fields)
# Get the unique table by the given fields. If the fields is omitted in tbuniq, function using all the fields.
# If trim or check is true, function will check each of non-key-fields to make sure it can be unique following the key-fields. If failed, function will remove this field (trim=true) or throw an error (check=true).
# When "field"=>function is given, these fields will be firstly excluded in the unique step, and then be added as To["field"]=grpfun(fun, ..., Ti["field"]). 
# * Arg#2=nothing means all fields except the ones specificly handled in 3+ args. Only allowed with 3+ args.
#See also: coo2mx, tbcomp, grpfun
#Xiong Jieyi, 30 Aug 2014>7 Oct 2014>Sep 8, 2017>9 May 2018>26 Aug 2019

tbuniq(T::Dict{Tk,}) where {Tk<:AbstractString} =
    rec(T,uniqr(tuple(values(T)...))[1]) 

function tbuniq(T::Dict{Tk,}, fld, Pr::Pair...; trim::Bool=false, check::Bool=false) where {Tk<:AbstractString}
    fg=nothing
    if isempty(Pr)
        isnothing(fld) && error("Arg#2 = nothing is not allowed when no \"field\"=>function behand. Just use tbuniq(T) instead of tbuniq(T, nothing)")
    else
        T0=T
        T=fd_i(T, map(x->x.first, Pr))
        fg=isnothing(fld) ? fastgrp(dict2tuple(T)) : fastgrp(dict2tuple(T, fld))
    end
    T=if trim || check
        isnothing(fld) && error("Arg#2 = nothing is not allowed when trim=true or check=true.")
        otherfd=setdiff(collect(keys(T)), fld)
        if isnothing(fg)
            fg=fastgrp(dict2tuple(T, fld))
        end
        isuni=map(otherfd) do ofd
            all(isrowsame, eachgrp(fg, T[ofd]))
        end
        if check
            if !all(isuni)
                t=otherfd[.!isuni]
                if length(t)>50
                    t=[t[1:50];"..."]
                end
                error("Below fields are not intra-group unique: "*join(t, ", "))
            end
            rec(T,uniqr(dict2tuple(T,fld))[1])
        else
            rec(fd(T, [fld; otherfd[isuni]]), fg.grpfirsti)
        end
    else
        if isnothing(fg)
            rec(T, uniqr(dict2tuple(T, fld))[1])
        else
            rec(T, fg.grpfirsti)
        end
    end
    for (ky, fun) in Pr
        T[ky]=grpfun(fun, fg, T0[ky])
    end
    T
end
tbuniq(T::Dict{Tk,}, fld::AbstractString, Pr::Pair...; kw...) where {Tk<:AbstractString} = tbuniq(T, [fld], Pr...; kw...)
tbuniq(fun::Function, T::Dict{Tk,}, fld, Pr1::Pair{Tk2, Colon}, Pr::Pair...; kw...)  where {Tk<:AbstractString, Tk2<:AbstractString} =tbuniq(T, fld, Pr1.first=>fun, Pr...; kw...)

tbuniq_i(T::Dict{Tk,},fld) where {Tk<:AbstractString} =
    rec(T,uniqr(tuple(values(fd_i(T,fld))...))[1])
tbuniq(T::Dict{Tk,}, fld::AbstractString; warg...) where {Tk<:AbstractString} =tbuniq(T, [fld]; warg...)
tbuniq_i(T::Dict{Tk,}, fld::AbstractString; warg...) where {Tk<:AbstractString} =tbuniq_i(T, [fld]; warg...)
export tbuniq, tbuniq_i
#}}

#{{ coo2mx
#[[ coo2mx ]]
# row_id, col_id, Matrix1, ... = coo2mx(row_id::Group, col_id::Group, Vector1, ...;
#                     default=...|nothing, defaults=(d1, d2, ...), row_order=..., col_order=...)
# Change 3-column-table data to a matrix. When default is given, all the inputs will using the same default value. When the default value is nothing, the emptyval() will be used. When default value is ignored, unassinged matrix element will trigger an error.
# See also: pysparse2coo, tbuniq, mxexpand, disent
# Xiong Jieyi, 26 Sep 2019

export coo2mx
function coo2mx(R::Group, C::Group, Vs::AbstractVector...; default=undef, defaults::Tuple=(), row_order::Union{Group, Nothing}=nothing, col_order::Union{Group, Nothing}=nothing)
    siR, iR=uniqr(R)
    sR=getrow(R, siR)
    siC, iC=uniqr(C)
    sC=getrow(C, siC)
    nrow=length(siR)
    ncol=length(siC)
    Ms=map(1:length(Vs)) do i
        cdef=isempty(defaults) ? default : defaults[i]
        if isnothing(cdef)
            cdef=emptyval(eltype(Vs[i]))
        end
        if cdef===undef
            Matrix{eltype(Vs[i])}(undef, nrow, ncol)
        else
            fill(cdef, nrow, ncol)
        end
    end
    ML=falses(nrow, ncol)
    for (i, (ir, ic)) in enumerate(zip(iR, iC))
        if ML[ir, ic]
            error("Duplicated row-column was detected.")
        end
        ML[ir, ic]=true
        for mi=1:length(Vs)
            Ms[mi][ir, ic]=Vs[mi][i]
        end
    end
    if default===undef && !all(ML)
        error(f"$1 out of $2 matrix elements are unassigned. Try set default value."(sum(.!ML), length(ML)))
    end
    if !isnothing(row_order)
        nrow=rnum(row_order)
        OMs=map(1:length(Vs)) do i
            cdef=isempty(defaults) ? default : defaults[i]
            if isnothing(cdef)
                cdef=emptyval(eltype(Vs[i]))
            end
            if cdef===undef
                Matrix{eltype(Vs[i])}(undef, nrow, ncol)
            else
                fill(cdef, nrow, ncol)
            end
        end
        Rfg=fastgrp(sR, row_order, Asorted=true, sortcheck=false)
        for gi=1:Rfg.grpnum
            if !isempty((rl=want(Rfg, gi);))
                for mi=1:length(Vs)
                    OMs[mi][gi, :]=Ms[mi][rl, :]
                end
            end
        end
        Ms=OMs
        sR=row_order
    end
    if !isnothing(col_order)
        ncol=rnum(col_order)
        OMs=map(1:length(Vs)) do i
            cdef=isempty(defaults) ? default : defaults[i]
            if isnothing(cdef)
                cdef=emptyval(eltype(Vs[i]))
            end
            if cdef===undef
                Matrix{eltype(Vs[i])}(undef, nrow, ncol)
            else
                fill(cdef, nrow, ncol)
            end
        end
        Cfg=fastgrp(sC, col_order, Asorted=true, sortcheck=false)
        for gi=1:Cfg.grpnum
            if !isempty((cl=want(Cfg, gi);))
                for mi=1:length(Vs)
                    OMs[mi][:, gi]=Ms[mi][:, cl]
                end
            end
        end
        Ms=OMs
        sC=col_order
    end
    (sR, sC, Ms...)
end
#}}

#{{ tbx
#[[ tbx ]]
#merged = tbx( A::Dict, B::Dict, key;
#             auniq=false, buniq=false, afull=false, bfull=false, akeep=false, bkeep=false,
#             exp=:AB(Default)|:A|:none, order=:none(Default)|:A|:B,
#             aprefix="Pre-"|function, bprefix=..., rename_conflict=true, delconflict=false,
#             fillempty=:none(Default)|:A|:B|:AB, emptyA=ds(...)|missing, emptyB=... )
#outA, outB = tbx(...; nomerge=true, ...)
#key could be "fieldname", ("fieldname_A", "fieldname_B"), (group1, group2), (group1, "fieldname_B"). If the key is omitted, the only shared key will be used.
#Or ... = tbx (A, B;< key=c"key1, key2, ..." >|< akey=c"...", bkey=c"..." >, ...)
#Merge two tables together. When key=(keyA, keyB), the keyB will be not in the output.
#The aprefix and bprefix will be added to the names of conflict fields (rename_conflict=true) or all non-key-fields (rename_conflict=false) in output table.
#delconflict: For the field with conflict values, throw an error (false) or delete this filed (true).
#aprefix/bprefix can also be a function, similar as rename!(function, ...). 
#emptyA and emptyB is Dict for the empty values in A and B. If any fields in emptyA/B are missing, emptyval(...) will be used. Note the field names in emptyA/B should be as the names in the input A/B, rather than the modified names in the output.
#When nomerge=true, tbx() will output two tables with marched rows. The key field(s) in B will be kept, and `delconflict` should only be false.
#See also: vcatr, rec, fd, fd_i, ft, tbx, tbuniq, tbuniq_i, dictconc, dictconc!, tbcomp
#Xiong Jieyi, 16 May 2014>11 Sep 2014>December 7, 2014>March 26, 2015>Jun 2, 2015>14 Oct 2017>13 Jan 2018>13 Mar 2018>7 Feb 2019>17 Sep 2019

export tbx
function __tbx__(A::Dict, B::Dict, key::Tuple{Group,Group},
             keynameA::AbstractVector{T1}, keynameB::AbstractVector{T2};
             exp::Symbol=:AB, fillempty::Symbol=:none,
             auniq::Bool=false, buniq::Bool=false, afull::Bool=false, bfull::Bool=false, akeep::Bool=false, bkeep::Bool=false, order::Symbol=:none, aprefix::Union{AbstractString, Function}="", bprefix::Union{AbstractString, Function}="", delconflict::Bool=false, rename_conflict::Bool=true, emptyA::Union{Dict, Missing}=ds(), emptyB::Union{Dict, Missing}=ds(), nomerge::Bool=false) where {T1<:AbstractString,T2<:AbstractString}

    # Input check
    in(fillempty, [:A,:B,:AB,:none]) || error("'fillempty' should only be :A, :B, :AB or :none.")
    in(order, [:A,:B,:none]) || error("'order' should only be :A, :B, :AB or :none.")
    if akeep
        afull=buniq=true
    end
    if bkeep
        bfull=auniq=true
    end
    if fillempty!=:none
        # delconflict && error("delconflict=true is not supported in fillempty model so far.")
        afull && (fillempty==:A || fillempty==:AB) && error("afull/akeep=true conflicts with fillempty=:A/:AB.")
        bfull && (fillempty==:B || fillempty==:AB) && error("bfull/bkeep=true conflicts with fillempty=:B/:AB.")
    end

    # Index match calculation
    if exp==:AB
        (uA,uB)=setxri(key...)
    elseif exp==:A
        (fillempty==:B || fillempty==:AB) && error("When exp=:A, fillempty must be :none or :A.")
        (uA,uB)=setxri(key...;uniqB=true)
    elseif exp==:none
        fillempty==:none || error("When exp=:none, fillempty must be :none.")
        (uA,uB)=intersectri(key...)
    else
        error("exp should only be :AB,:A or :none.")
    end

    auniq && (isuniqr(uB) || error("The A key is not row-unique"))
    buniq && (isuniqr(uA) || error("The B key is not row-unique"))

    afull && (length(uniqr(uA)[1])==rownum(A) || error("A is not full."))
    bfull && (length(uniqr(uB)[1])==rownum(B) || error("B is not full."))

    if order==:A && (fillempty==:none || fillempty==:A)
        (uA,idx)=sortr(uA)
        uB=uB[idx]
    elseif order==:B && (fillempty==:none || fillempty==:B)
        (uB,idx)=sortr(uB)
        uA=uA[idx]
    end
    oA=rec(A,uA)
    oB=rec(B,uB)
    
    #Fill empty
    if fillempty!=:none
        sharedfields=if isempty(aprefix) && isempty(bprefix)
            # setdiff(intersect(keys(A), keys(B)), keyname)
            intersect(setdiff(keys(A), keynameA), setdiff(keys(B), keynameB))
        elseif rename_conflict
            filter(x->isequal(oA[x], oB[x]), intersect(setdiff(keys(A), keynameA), setdiff(keys(B), keynameB)))
        else
            String[]
        end
        if fillempty==:B || fillempty==:AB
            cAi=setdiff(1:rownum(key[1]), uA)
            cA=getrow(A, cAi)
            oA=vcatr(oA, cA)
            
            # cB=emptyrows(B, length(cAi))
            cB=tb()
            for (k, v) in B
                cB[k]=if ismissing(emptyB)
                    missing
                elseif haskey(emptyB, k)
                    reprow(emptyB[k], length(cAi))
                else
                    emptyrows(v, length(cAi))
                end
            end
            
            for x in sharedfields
                cB[x]=cA[x]
            end
            oB=vcatr(oB, cB)

            if order==:A
                idx=sortri([uA; cAi])
                oA=getrow(oA, idx)
                oB=getrow(oB, idx)
            elseif order==:B
                uB=[uB; rownum(key[2])+(1:length(cAi))]
            end
        end
        if fillempty==:A || fillempty==:AB
            cBi=setdiff(1:rownum(key[2]), uB)
            cB=getrow(B, cBi)
            oB=vcatr(oB, cB)
            
            # cA=emptyrows(A, length(cBi))
            cA=tb()
            for (k, v) in A
                cA[k]=if ismissing(emptyA)
                    missing
                elseif haskey(emptyA, k)
                    reprow(emptyA[k], length(cBi))
                else
                    emptyrows(v, length(cBi))
                end
            end
            
            if !isempty(keynameA)
                for (cka, ckb) in zip(keynameA, keynameB)
                    cA[cka]=cB[ckb]
                end
                # cA[keyname[1]]=cB[keyname[2]]
            end
            for x in sharedfields
                cA[x]=cB[x]
            end
            oA=vcatr(oA, cA)

            if order==:B
                idx=sortri([uB; cBi])
                oA=getrow(oA, idx)
                oB=getrow(oB, idx)
            end
        end
    end
    
    #Note that I cannot operate delete! or rename! on A or B directly, since A and B are not copied.
    #Remove B key field
    if !nomerge && !isempty(keynameB)
        for ckb in keynameB
            delete!(oB, ckb)
        end
    end
    # Field rename
    if !(isempty(aprefix) && isempty(bprefix))
        if rename_conflict
            delconflict && error("rename_conflict=true conflicts with delconflict=true")
            confkys=collect(filter(ky->!isequal(oA[ky], oB[ky]), intersect(keys(oA), keys(oB))))
            if !isempty(aprefix)
                oA=rename!(isa(aprefix, AbstractString) ? x->aprefix*x : aprefix, oA, confkys)
            end
            if !isempty(bprefix)
                oB=rename!(isa(bprefix, AbstractString) ? x->bprefix*x : bprefix, oB, confkys)
            end
        else
            if !isempty(aprefix)
                if isempty(keynameA)
                    oA=rename!(isa(aprefix, AbstractString) ? x->aprefix*x : aprefix, oA)
                else
                    oA=rename!(isa(aprefix, AbstractString) ? x->aprefix*x : aprefix, oA, exclude=keynameA)
                end
            end
            if !isempty(bprefix)
                oB=rename!(isa(bprefix, AbstractString) ? x->bprefix*x : bprefix, oB)
            end    
        end
    end
    if nomerge
        delconflict && error("nomerge=true conflicts with delconflict=true.")
        (oA, oB)
    else
        dictconc!(oA,oB; delconflict=delconflict)
    end
end
tbx(A::Dict, B::Dict, key::Tuple{T1,T2}; P...) where {T1<:AbstractString,T2<:AbstractString} =
    __tbx__(A, B, (A[key[1]], B[key[2]]), [key[1]], [key[2]]; P...)
tbx(A::Dict, B::Dict, key::Tuple{T,Group}; P...) where {T<:AbstractString} =
    __tbx__(A, B, (A[key[1]], key[2]), [key[1]], String[]; P...)
tbx(A::Dict, B::Dict, key::Tuple{Group, T}; P...) where {T<:AbstractString} =
    __tbx__(A, B, (key[1], B[key[2]]), String[], [key[2]]; P...)
tbx(A::Dict, B::Dict, key::Tuple{Group, Group}; P...)=__tbx__(A, B, key, String[], String[]; P...)
tbx(A::Dict, B::Dict, key::AbstractString; P...)=tbx(A, B, (key, key); P...)
function tbx(A::Dict,B::Dict; key=String[], akey::AbstractVector{T1}=key, bkey::AbstractVector{T2}=key, wargs...) where {T1<:AbstractString,T2<:AbstractString}
    if isempty(akey)
        sharekey=intersect(keys(A),keys(B))
        length(sharekey)==1 || error("More than one key is shared, or no key is shared. Please assign the linking key.")
        tbx(A,B,first(sharekey); wargs...)
    else
        length(akey)==length(bkey) || error("akey=... and bkey=... should be in the same lengths.")
        __tbx__(A, B, (dict2tuple(A, akey), dict2tuple(B, bkey)), akey, bkey; wargs...)
    end
end

#}}

#{{ tbcomp
#[[ tbcomp ]]
#T|F=tbcomp(Tb_Small, Tb_Big[, c"fields"]; eq=false)
#Test if the combinations of given fields in first table is exist in second table (eq=fasle), or if their uniqued combinations are the same (eq=true). The repeated combinations will be uniqued.
#See also: vcatr, rec, fd, fd_i, ft, tbx, tbuniq, tbuniq_i, dictconc, dictconc!
#Xiong Jieyi, January 8, 2015 >Apr 19, 2016

export tbcomp
function tbcomp(A::Dict, B::Dict, fds::AbstractVector; eq::Bool=false)
    @assert(istable(A),"A is not a table.")
    @assert(istable(B),"B is not a table.")
    GA=dict2tuple(A,fds)
    GB=dict2tuple(B,fds)
    @assert(typeof(GA)==typeof(GB),"Type of two tables is inconsistent. A: $(typeof(GA)), B: $(typeof(GB)).")
    if eq
        unival(GA)==unival(GB)
    else
        all(ismemberr(GA,GB))
    end
end
tbcomp(A::Dict, B::Dict; wargs...)=tbcomp(A::Dict, B::Dict, collect(intersect(keys(A), keys(B))); wargs...)
#}}

#{{ readtb writetb

#[[ readtb writetb ]]
# table = readtb(filename; dm='\t', head=c"data1,N::data2,...", skipline=0, autotype=false)
# writetb(filename, Table, "key1[=>rename1], key2[=>rename1], ..."; dm='\t', with_type=true)
# writetb(filename, Table, c"key1, key2, ..."; dm='\t', with_type=true) #not parse '=>'.
#Read tsv file, using the first line as table field name (or use the head=c"..." instead). The #-started lines will be regarded as comments, but not for #::[Any space]. So it is recommended (but not required) to start the head line like this so other software will regarded it as comment line.
#The elements of first line could be started with below symbol. Otherwise, data will be readed as string:
#N:: Int number *
#n:: Int number, empty will be converted to 0.
#F:: Float number, only NaN will be converted to NaN. *
#f:: Float number, empty and NaN, NA, na, Na, n.a., "NaN", "NA", "na", "Na", "n.a.", "NaN" will be converted to NaN.
#B:: Bool number (lowercase of values should be any of "t,true,f,false,y,yes,n,no,1,0") *
#S:: Symbol *
#T:: String, keeping quote mark (") also. In default, the "..." will be striped.
#C:: Char (Only support one charactor) *
#X:: Ignore this column
#The head could also be < , which means this column will be concaterate with its left column and output a matrix. *
# --Items with * are supported by writetb.
#In the case header is in a shorter length than column number, set the last header as "..." will ignore the right columns, set the last header as <<< will merge all the right columns (equal to '<' for all right columns), Otherwise an error will be throwed.
#autotype: transfer all string fields to BitVector, Int or Float64 once possible ("NaN" is accepted as Float64 but "nan", et al. is not accepted; BitVector only accept TRUE, True, true, FALSE, False, false).
#with_type: weather write with the type mark, e.g. N::xxx. If with_type is false, 'field < ...' will be replaced to 'field.1 field.2 ...'.
#See also: filefun
#Xiong Jieyi, February 22, 2015 >Sep 4, 2015 >Aug 17, 2016 >15 May 2017>31 May 2017>22 Nov 2018>17 Dec 2018

export readtb, writetb

function __parse_tb_head__(head::Vector{T}) where {T<:AbstractString}
    
    #In the case field name is quited by `"':
    for i=1:length(head)
        if length(head[i])>=2 && head[i][1]=='"' && head[i][end]=='"'
            head[i]=head[i][2:end-1]
        end
    end
    
    if head[end]=="..."
        rightCols=-1
        head=head[1:end-1]
    elseif head[end]=="<<<"
        rightCols=1
        head=head[1:end-1]
    else
        rightCols=0
    end
    dtfun=Function[]
    fdnm=ASCIIString[]
    colidx=Array{Int}[]
    for (i,cf) in enumerate(head)
        if startswith(cf,"X::")
            continue
        elseif cf=="<"
            push!(colidx[end],i)
            continue
        else
            push!(colidx,[i])
        end
        if startswith(cf,"N::")
            push!(dtfun,int)
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"n::")
            push!(dtfun,x->isempty(x) ? 0 : int(x))
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"F::")
            push!(dtfun,float)
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"f::")
            push!(dtfun,x->(isempty(x)||in(x,["NA","na","Na","n.a.","\"NA\"","\"na\"","\"Na\"","\"n.a.\"","\"NaN\""])) ? NaN : float(x))
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"B::")
            push!(dtfun,x::AbstractString->begin
                  if lowercase(x) in c"t,true,1,y,yes"
                  true
                  else
                  @assert(lowercase(x) in c"f,false,0,n,no",
                          "Invalid logical symbol $x")
                  false
                  end
                  end)
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"C::")
            push!(dtfun,x::AbstractString->begin
                  @assert(length(x)==1,
                          "Char value should only be one charactors.")
                  x[1]
                  end)
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"S::")
            push!(dtfun,x->Symbol(x))
            push!(fdnm,cf[4:end])
        elseif startswith(cf,"T::")
            push!(dtfun,ascii)
            push!(fdnm,cf[4:end])
        else
            push!(dtfun,x::AbstractString->(length(x)>=2 && x[1]=='"' && x[end]=='"') ? ascii(x[2:end-1]) : ascii(x))
            push!(fdnm,cf)
        end
    end
    (dtfun,fdnm,colidx,rightCols)
end

function readtb(fn; dm='\t', head::Vector{Ts}=ASCIIString[], skipline=0, autotype::Bool=false) where {Ts<:AbstractString}
    column=Ref(-1)
    rightCols=Ref(0)
    if isempty(head)
        dtfun=Function[]
        fdnm=ASCIIString[]
        colidx=Array{Int}[]
    else
        (dtfun,fdnm,colidx,rightCols.x)=__parse_tb_head__(head)
        column.x=length(head)-(rightCols.x!=0)
    end
    T=filefun(fn; input_rowno=true, skipline=skipline)do no,ln
        if (length(ln)>1) && ln[1]=='#'
            m=match(r"^#\:\:\s+",ln)
            if m==nothing
                return nothing
            else
                ln=ln[length(m.match)+1:end]
            end
        end
        C=split(ln,dm)
        if no==1 && isempty(head)
            (dtfun,fdnm,colidx,rightCols.x)=__parse_tb_head__(C)
            column.x=length(C)-(rightCols.x!=0)
            nothing
        else
            if rightCols.x==0 ? length(C)!=column.x : length(C)<column.x
                error("Invalid column number $(length(C)) in row $no.")
            end
            T=ds()
            for i=1:length(colidx)
                cidx=colidx[i]
                if i==length(colidx) && rightCols.x==1
                    T[fdnm[i]]=map(dtfun[i],C[cidx[1]:end])
                elseif length(cidx)==1
                    T[fdnm[i]]=dtfun[i](C[cidx[1]])
                else
                    T[fdnm[i]]=map(dtfun[i],C[cidx])
                end
            end
            T
        end
    end
    if autotype
        for (k, v) in T
            if eltype(v)<:AbstractString
                v=strip.(v)
                if all(x->in(x, c"TRUE, True, true, FALSE, False, false"), v)
                    T[k]=BitVector(parse.(Bool, lowercase.(v)))
                elseif all(x->occursin(r"^\-?\d+$", x), v)
                    T[k]=parse.(Int, v)
                elseif all(x->occursin(r"^[e\d\-\+\.]+$|^NaN$", x), v)
                    try
                        T[k]=parse.(Float64, v)
                    catch
                    end
                end
            end
        end
    end
    T
end

function writetb(filename::AbstractString, T::Dict{Ts1, Any}, keyls::Union{AbstractString, AbstractVector{Ts2}}=sort(collect(keys(T))); dm='\t', with_type::Bool=true) where {Ts1<:AbstractString, Ts2<:AbstractString}
    istable(T) || error("Input is not a table.")
    if isa(keyls, AbstractString)
        keyls=strip.(split(keyls, ','))
        parseflag=true
    else
        parseflag=false
    end
    cols=map(enumerate(keyls)) do (i, ky)
        if parseflag
            t=split(ky, "=>")
            if length(t)==2
                ky, wky=strip.(t)
                keyls[i]=ky
            else
                wky=ky
            end
        else
            wky=ky
        end
        occursin(dm, wky) && error("Delim char (dm=$dm) is not allowed in key $wky.")
        vv=T[ky]
        if with_type
            vt=eltype(vv)
            pf=if vt<:Bool
                "B::"
            elseif vt<:Integer
                "N::"
            elseif vt<:AbstractFloat
                "F::"
            elseif vt<:Symbol
                "S::"
            elseif vt<:Char
                "C::"
            elseif vt<:AbstractString
                ""
            else
                error("Unsupported type $vt in field $ky.")
            end
            ["$pf$wky", fill("<", size(vv, 2)-1)...]
        else
            if size(vv, 2)>1
                f"$wky.$1".(1:size(vv, 2))
            else
                [wky]
            end
        end
    end
    writedlm(filename, asmx(Pair.(cols, dict2tuple(T, keyls))...), dm)
end

#}}

#{{ rename!
#[[ rename! ]]
#Dict=rename!(Dict, Key1=>NewKey1, Key2=>NewKey2, ...) or rename!(Dict, Key_vec, NewKey_vec)
#... =rename!(Dict, Regex=>Replace_string[, c"field1,field2,..."|Regex];
#            exclude=Vector|Tuple|Set|AbstractString|Regex) #use replace() for fields.
#... =rename!(function, Dict[, key_vec|Regex]; exclude=...) #replace key to function(key).
#Change the key in Dict to the new name, and return the changed Dict.
#e.g. rename!(D, regex1=>s"...", regex2; exclude=regex3), the keys of D matches regex2 and not match regex3 will be renamed as replace(key, regex1 => s"...").
#For compat reason, rename!(Dict, "key"|Regex, newkey, ...) is still supported but it is not recommended.
#See also: fd, fd_i, delete!
#Xiong Jieyi, 10 Sep 2014>December 12, 2014>March 25, 2015>Jun 2, 2015>Nov 12, 2015>29 Apr 2019

export rename!
function rename!(D::Dict{Tk,},A,B) where {Tk<:AbstractString}
    if A!=B
        D[B]=D[A]
        delete!(D,A)
    else
        D
    end
end
function rename!(D::Dict{Tk,},A::T,B::T) where {Tk<:AbstractString, T<:Union{Vector,Tuple}}
    for i=1:length(A)
        rename!(D,A[i],B[i])
    end
    D
end
function rename!(D::Dict{Tk,},A::Regex,B::AbstractString,fds::Union{AbstractVector,Tuple}=collect(keys(D)); exclude::Union{AbstractVector,AbstractSet,Tuple,AbstractString,Regex}=()) where {Tk<:AbstractString}
    if !isempty(exclude)
        if isa(exclude,AbstractString)
            fds=fds[fds.!=exclude]
        elseif isa(exclude, Regex)
            fds=filter(x->!occursin(exclude, x), fds)
        else
            fds=setdiff(fds,exclude)
        end
    end
    nfds=map(x->replace(x, A => B),fds)
    l=fds.!=nfds
    fds=fds[l]
    nfds=nfds[l]
    isuniqr(nfds) || error("Duplicated fieldname occurred after replacement.")
    rename!(D,fds,nfds)
end
function rename!(F::Function, D::Dict{Tk,}, fds::Union{AbstractVector,Tuple}=collect(keys(D)); exclude::Union{AbstractVector,AbstractSet,Tuple,AbstractString,Regex}=())  where {Tk<:AbstractString}
    if !isempty(exclude)
        if isa(exclude,AbstractString)
            fds=fds[fds.!=exclude]
        elseif isa(exclude, Regex)
            fds=filter(x->!occursin(exclude, x), fds)
        else
            fds=setdiff(fds,exclude)
        end
    end
    nfds=map(F, fds)
    l=fds.!=nfds
    fds=fds[l]
    nfds=nfds[l]
    isuniqr(nfds) || error("Duplicated fieldname occurred after replacement.")
    rename!(D,fds,nfds)
end
rename!(F::Function, D::Dict{Tk,}, ft::Regex; kw...) where {Tk<:AbstractString}=rename!(F, D, collect(filter(ft, keys(D))); kw...)
rename!(D::Dict{Tk,}, P::Pair, ft::Regex; kw...) where {Tk<:AbstractString} = rename!(D, P, collect(filter(ft, keys(D))); kw...)
rename!(D::Dict{Tk,}, P::Pair{Regex, T}, arg...; kw...) where {Tk<:AbstractString, T<:Union{AbstractString, Function}} = rename!(D, P.first, P.second, arg...; kw...)
function rename!(D::Dict{Tk,}, pr1::Pair, pr::Pair...) where {Tk<:AbstractString}
    for (a,b) in (pr1, pr...)
        rename!(D,a,b)
    end
    D
end
#}}

#{{ tb2DataFrame, DataFrame2tb
#[[ tb2DataFrame ]]
# dataframe = tb2DataFrame( Table[, fields ])
#Convert table to DataFrame
#If field "F" is a matrix, its column will be splitted into :F_1, :F_2, ....
#See also: readtb, sk, sk_s, pandasDF, DataFrame2tb
#Xiong Jieyi, Jul 10, 2015 > 2 May 2017

export tb2DataFrame
#<v0.6# function tb2DataFrame{Tk<:AbstractString}(T::Dict{Tk,}, fields=keys(T))
function tb2DataFrame(T::Dict{Tk,}, fields=keys(T)) where {Tk<:AbstractString}
    D=importpkg(:DataFrames, preloaded=true).DataFrame()
    function foo(k,v::AbstractMatrix)
        for i=1:size(v,2)
            D[Symbol("$(k)_$i")]=v[:,i]
        end
    end
    function foo(k,v)
        D[Symbol(k)]=v
    end
    for k in fields
        foo(k,T[k])
    end
    D
end

#[[ DataFrame2tb ]]
# dataframe = DataFrame2tb( DataFrame )
#Convert DataFrame to table.
#NA in data will be converted to empty value, and a warning will occur.
#See also: readtb, sk, sk_s, pandasDF, tb2DataFrame
#Xiong Jieyi, Oct 6, 2015

export DataFrame2tb
function DataFrame2tb(D::DataFrame)
    T=tb()
    # isna(x)=invokelatest(importpkg(:DataFrames).isna, x)
    # @pkgfun(DataFrames,DataFrame,isna,getindex,names,collect)
    # DF=importpkg(:DataFrames)
    for k in names(D)
        F=D[k]
        l=ismissing.(F)
        if any(l)
            emp=emptyval(F[findfirst(.!l)])
            V=fill(emp,length(F))
            V[.!l]=F[.!l]
            @warn(f"$1/$2 are NA in field $3. Replaced by empty value `$4'($5)."(sum(l),length(l),k,emp,typeof(emp)))
        else
            # V=[F...]
            V=collect(F)
        end
        T[string(k)]=V
    end
    T
end
#}}

#{{ pandasDF pandas2tb
#[[ pandasDF pandas2tb ]]
# dataframe = pandasDF( Table )
# dataframe = pandasDF(tb_style_params...) # =pandasDF(tb(...))
# dataframe = pandasDF(data; index=..., colomns=...) # =pandas.DaraFrame(...)
# table = pandas2tb( pandas_dataframe[, "index_name" ]) # Requiring PyCall.
#Convert table to python pandas DataFrame or vise versa.
#For pandasDF, if field "F" is a matrix, its column will be splitted into "F:1", "F:2", ....
#For pandas2tb, the index will also be converted to a column if exists.
#See also: readtb, sk, sk_s, tb, ds, tb2DataFrame
#Xiong Jieyi, Aug 16, 2015 > Jan 12, 2016 >26 Jun 2019

export pandasDF, pandas2tb
#<v0.6# function pandasDF{Tk<:AbstractString}(T::Dict{Tk,})
function pandasDF(T::Dict{Tk,}) where {Tk<:AbstractString}
    D=tb()
    function foo(k,v::AbstractMatrix)
        for i=1:size(v,2)
            D["$k:$i"]=v[:,i]
        end
    end
    function foo(k,v)
        D[k]=v
    end
    for (k,v) in T
        foo(k,v)
    end
    importpy(:pandas).DataFrame(D)
end
pandasDF(args...;wargs...)=pandasDF(tb(args...;wargs...))
pandasDF(D::AbstractArray;wargs...)=importpy(:pandas).DataFrame(D;wargs...)

function pandas2tb(pd, index_name::Union{AbstractString, Nothing}=pd.index.name)
    # T=Main.PyCall.pycall(pd.to_dict, Dict{String, Any}, "list")
    T=tb()
    for ky in pd
        T[ky]=getproperty(pd, ky).to_list()
    end
    if !isnothing(index_name)
        T[index_name]=pd.index.to_list()
    end
    T
end
#}}

#{{ asmx
#[[ asmx ]]
# Matrix = asmx(Matrix[, ...]; row=row_header_vec, col=col_header_vec, rowexp=N)
# Or ... = asmx(c"header1, header2, ..."=>Matrix|Tuple_of_columns, "header3"=>Vector, ...; row=...)
#Composite a matrix, and ddd the row header and column header to a matrix.
#In the pairs input mode, the length of header and column number will be checked in each pair. It could avoid programmer's mistake.
#If rowexp>0, the single-row input will be expand to N, and other input should has N row.
#See also: writedlm, tocb, coo2mx
#Xiong Jieyi, Jul 6, 2017 > 28 Oct 2017 > 21 Oct 2019

export asmx
function asmx(M::Union{AbstractVector, AbstractMatrix}; row::Union{AbstractVector,Void}=nothing, col::Union{AbstractVector,Void}=nothing,corner="")
    if row==nothing && col==nothing
        M
    elseif row==nothing
        size(M,2)==length(col) || error(f"Column header vector length should be $1 but it is $2."(size(M,2), length(col)))
        [trsp(col); M]
    elseif col==nothing
        size(M,1)==length(row) || error(f"Row header vector length should be $1 but it is $2."(size(M,1), length(row)))
        [row M]
    else
        size(M,1)==length(row) || error(f"Row header vector length should be $1 but it is $2."(size(M,1), length(row)))
        size(M,2)==length(col) || error(f"Column header vector length should be $1 but it is $2."(size(M,2), length(col)))
        [trsp([corner;col]); [row M]]
    end
end

function asmx(M1, Ms...; rowexp::Integer=0, wargs...)
    C=Any[M1, Ms...]
    if rowexp>0
        for (i, cM) in enumerate(C)
            rn=rownum(cM)
            if rn==1
                C[i]=reprow(cM, rowexp)
            else
                rn==rowexp || error("The $(i+1)-th input ($rn rows) is neither single row nor match the row number $rowexp.")
            end
        end
    end
    asmx(hcat(C...); wargs...)
end

function asmx(Ps::Pair... ; wargs...)
    for (p1, p2) in Ps
        if isa(p1, AbstractArray)
            length(p1)==(isa(p2, Tuple) ? sum(size.(p2, 2)) : size(p2,2)) || error(f"Header \"$1\"...\"$2\" has $3 elements, but matrix has $4 columns."(p1[1],p1[end],length(p1),size(p2,2)))
        else
            isa(p2, AbstractVector) || (size(p2,2)==1) || error("Single header \"$p1\" does not match a vector column.")
        end
    end
    M=hcat(map(x->isa(x.second, Tuple) ? hcat(x.second...) : x.second, Ps)...)
    header=vcat(map(x->x.first, Ps)...)
    asmx(M; col=header, wargs...)
end
#}}
