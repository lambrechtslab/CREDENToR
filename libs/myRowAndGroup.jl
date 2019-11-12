export names
#eachcol and rename! is conflicting with DataFrames.
import Base: keys
keys(X::DataFrame)=names(X) #For the autocompletion.

# const Group=Union{AbstractVector,AbstractMatrix,Tuple,Dict,Any} #Changed in 13 Mar 2018
const Group=Union{AbstractVector,AbstractMatrix,Tuple{Vararg{Union{AbstractVector,AbstractMatrix,Tuple,AbstractDataFrame}}},Dict,AbstractDataFrame}
export Group

#{{ islessr2 isequalr2

#[[ islessr2 isequalr2 ]]
# islessr2(A, row_a, B, row_b) only check the row of a and b.
# isequalr2(A, row_a, B, row_b) check if the row of a and b are equal. NaN is equal to NaN.
#Low-level function to compare two rows from two groups. It is faster than isequal(getrow(A, row_a),getrow(B, row_b))
#Note that for the speed reason, if size(A,2)<size(B,2), function will only compare the shared columns and not throw a error.
# See also: sortr, uniqr, isuniqr
#Xiong Jieyi,10 Jul 2014 > Sep 23, 2016 > 26 Apr 2017

function islessr2(A::AbstractMatrix,ia::Integer,B::AbstractMatrix,ib::Integer)
    lA=size(A,2)
    lB=size(B,2)
    N=lA<lB ? lA : lB
    for i=1:N
        if A[ia,i]<B[ib,i]
            return true
        elseif B[ib,i]<A[ia,i]
            return false
        end
    end
    return lA<lB
end
function islessr2(A::T,ia::Integer,B::T,ib::Integer) where {T<:Tuple}
    lA=length(A)
    lB=length(B)
    N=lA<lB ? lA : lB
    for i=1:N
        if islessr2(A[i],ia,B[i],ib)
            return true
        elseif islessr2(B[i],ib,A[i],ia)
            return false
        end
    end
    return lA<lB
end
islessr2(A::AbstractVector,ia::Integer,B::AbstractVector,ib::Integer)=A[ia]<B[ib]
function islessr2(A::T,ia::Integer,B::T,ib::Integer) where {T<:Dict} #Add in Sep 23, 2016
    tA=tuple(values(A)...)
    tB=tuple(map(x->B[x],keys(A))...)
    islessr2(tA,ia,tB,ib)
end
function islessr2(A::T,ia::Integer,B::T,ib::Integer) where {T<:AbstractDataFrame}
    N=size(A, 2)
    size(B, 2)==N || error("Two DataFrames are not in the same lengths.")
    for i=1:N
        if islessr2(A[i], ia, B[i], ib)
            return true
        elseif islessr2(B[i], ib, A[i], ia)
            return false
        end
    end
    return false
end

function isequalr2(A::AbstractMatrix,ia::Integer,B::AbstractMatrix,ib::Integer)
    for i=1:size(A,2)
        if !isequal(A[ia,i],B[ib,i])
            return false
        end
    end
    return true
end
isequalr2(A::AbstractVector,ia::Integer,B::AbstractVector,ib::Integer)=isequal(A[ia],B[ib])
function isequalr2(A::Tuple,ia::Integer,B::Tuple,ib::Integer) #Changed in 8 Dec 2017
    lenA=length(A)
    lenA==length(B) || error("Two tuples are not in the same lengths.")
    for i=1:lenA
        if !isequalr2(A[i],ia,B[i],ib)
            return false
        end
    end
    return true
end
function isequalr2(A::T,ia::Integer,B::T,ib::Integer) where {T<:Dict} #Add in Sep 23, 2016
    tA=tuple(values(A)...)
    
    # tB=tuple(map(x->B[x],keys(A))...) #Modified in 31 Oct 2018 for Julia1.0. map doesn't support keys(...) in Julia1.0
    tB=tuple(broadcast(x->B[x], keys(A))...)
    
    isequalr2(tA,ia,tB,ib)
end
function isequalr2(A::T,ia::Integer,B::T,ib::Integer) where {T<:AbstractDataFrame}
    N=size(A, 2)
    size(B, 2)==N || error("Two DataFrames are not in the same lengths.")
    for i=1:N
        if !isequalr2(A[i], ia, B[i], ib)
            return false
        end
    end
    return true
end
export islessr2,isequalr2

#}}

#{{ isgroup isrowdt

#[[ isgroup ]]
#T|F=isgroup(X)
#Return if input is a group. If input is a Tuple or Dict{String,}, type and rownum of each element will be checked recurly.
#Any non-group input is false. e.g. isgroup(1)=false
#See also: istable, rownum, iseqrow
#Xiong Jieyi, Jun 16, 2015 > Sep 23, 2016

export isgroup
function isgroup(X::Tuple)
    if isempty(X) || !isgroup(X[1])
        return true
    end
    rN=rownum(X[1])
    for v in X[2:end]
        if !(isgroup(v) && (rN==rownum(v)))
            return false
        end
    end
    return true
end
isgroup(X::Dict{T,}) where T = T<:AbstractString && isgroup(dict2tuple(X))
isgroup(::AbstractVector{Any})=false
isgroup(::AbstractVector{Tuple})=false
isgroup(::AbstractVector{Dict})=false
isgroup(::AbstractVector{T}) where {T<:AbstractArray} = false
isgroup(::AbstractMatrix{Any})=false
isgroup(::AbstractMatrix{Tuple})=false
isgroup(::AbstractMatrix{Dict})=false
isgroup(::AbstractMatrix{T}) where {T<:AbstractArray} = false
isgroup(::AbstractVector)=true
isgroup(::AbstractMatrix)=true
isgroup(::AbstractDataFrame)=true
isgroup(::Any)=false

#[[ iseqrownum ]]
#T|F=iseqrownum(X1, X2, ...)
#Return if all inputs is the data with identical row number. If input is a Tuple or Dict, rownum of each element will be checked recurly.
#Any non-Tuple or non-Dict input is true. e.g. iseqrow(1)=true
#See also: isgroup, istable
#Xiong Jieyi, Jun 16, 2015

export iseqrownum
function iseqrownum(X::Tuple)
    if isempty(X)
        return true
    end
    rN=rownum(X)
    for v in X
        if !(iseqrownum(v) && (rN==rownum(v)))
            return false
        end
    end
    return true
end
function iseqrownum(X::Dict)
    if isempty(X)
        return true
    end
    rN=-1
    for v in values(X)
        if !iseqrownum(v)
            return false
        end
        if rN<0
            rN=rownum(v)
        elseif rN!=rownum(v)
            return false
        end
    end
    return true
end
iseqrownum(::Any)=true
iseqrownum(X1,X...)=iseqrownum((X1,X))

#}}

#{{ getrow setrow! setrowexp! rownum

#[[ getrow ]]
#getrow(Group, RowN)
#Return the rows. If input is a Tuple, return a Tuple also.
#getrow(Dict,I) is equal to rec(Dict,I)
#When Group is a Any-Vector, Array-Vector, Tuple-Vector or Dict-Vector, the output is one of the same container even if RowN is a number.
#See also: unslice, rownum, setrow!, setrow2!, vcatr, colmx2vec, iseqrow, isgroup
#Xiong Jiyei, 14 May 2014 > May 23, 2015>Jun 16, 2015>Dec 20, 2015

getrow(X::AbstractVector{T},I::Integer) where {T<:AbstractArray} = X[[I]]
getrow(X::AbstractVector{Any},I::Integer)=X[[I]]
getrow(X::AbstractVector{T},I::Integer) where {T<:Tuple} = X[[I]]
getrow(X::AbstractVector{T},I::Integer) where {T<:Dict} = X[[I]]
getrow(X::AbstractDataFrame, I)=X[I, :]
getrow(X::AbstractVector,I)=X[I]
getrow(X::AbstractMatrix,I::AbstractVector)=X[I,:]
getrow(X::AbstractMatrix,I::Real)=X[[I],:] #Add in Oct 19, 2016 for Julia 0.5
getrow(X::AbstractVector,I::AbstractMatrix)=invoke(getrow,(Group,AbstractMatrix),X,I)
getrow(X::AbstractMatrix,I::AbstractMatrix)=invoke(getrow,(Group,AbstractMatrix),X,I)
function getrow(X::Union{Group,Dict},I::AbstractMatrix)
    @assert(size(I,2)==1,"Only support vector or one-colunm matrix.")
    getrow(X,vec(I))
end
getrow(X::Tuple,I::AbstractMatrix)=getrow(X,I[:])
getrow(X::Tuple,I)=map(x->getrow(x,I),X)
getrow(X::Dict,I::Union{Integer,AbstractVector})=dictfun(x->getrow(x,I),X)
getrow(X::Union{AbstractString,Number},I)=getrow([X],I)
export getrow

#[[ unslice unslice! ]]
# getrow(X, [i]) = unslice|unslice!(getrow(X, i))
#Only when working at a Dict, unslice and unslice! have different.
#See also: getrow
#Xiong Jieyi, Nov 25, 2015 > 16 May 2019

function unslice(X::AbstractArray)
    @assert(rownum(X)==1,"Invalid slice.")
    X
end
unslice(X::Tuple)=map(unslice,X)
unslice(X::Dict)=dictfun(unslice,X)
unslice(X)=[X]

unslice!(X::Tuple)=map(unslice!,X)
unslice!(X::Dict)=dictfun!(unslice!, X, replace=true)
unslice!(X)=unslice(X)
export unslice, unslice!


#[[ setrow! setrowexp!]]
#setrow!(X, I, V) or setrowexp!(X, I, V)
#X[I,:]=V. If X is a Tuple, function will do in every elements of X.
#setrowexp! will expand the row of V (V should be a one-row matrix).
#When X is a Any vector or array-vector, and I is a number, V must be a singleton homogenious vector.
#X and V could be Group or Table.
#See also: getrow, rownum, vcatr
#Xiong Jiyei, 14 May 2014>11 Sep 2014>Jun 16, 2015

function setrow!(X::AbstractVector,I::Integer,V)
    X[I]=V
    X
end
function setrow!(X::AbstractVector,I::AbstractVector,V)
    X[I].=V
    X
end
function setrow!(X::AbstractVector,I::AbstractVector,V::AbstractVector)
    X[I].=V
    X
end
function setrow!(X::AbstractMatrix,I,V)
    X[I, :]=V
    X
end
function setrow!(X::AbstractDataFrame, I, V::AbstractDataFrame)
    X[I, :]=V
    X
end
function setrow!(X::Tuple,I,V::Tuple)
    for i=1:length(X)
        setrow!(X[i],I,V[i])
    end
    X
end
setrow!(X::Dict,I,V::Dict)=dictfun!((x,v)->setrow!(x,I,v),X,V, replace=false)
setrow!(X::AbstractVector{T},I::Integer,V::AbstractVector{T}) where {T<:Union{Tuple,Dict,Any}} = setrow!(X,[I],onlyone(V))
setrow!(X::AbstractVector{T},I::Integer,V::AbstractVector{TV}) where {T<:AbstractArray,TV<:AbstractArray} = setrow!(X,[I],onlyone(V))
setrow!(X::AbstractVector{T},I::Integer,V::T) where {T<:Union{Tuple,Dict,Any}} = setrow!(X,[I],V)
setrow!(X::AbstractVector{T},I::Integer,V::TV) where {T<:AbstractArray,TV<:AbstractArray} = setrow!(X,[I],V)
function setrow!(X::AbstractVector{T},I,V::AbstractVector{T}) where {T<:Union{AbstractArray,Tuple,Dict,Any}}
    X[I]=rownum(V)==1 ? onlyone(V) : V
    X
end

setrowexp!(X,I,V)=setrow!(X,I,V)
setrowexp!(X::AbstractMatrix,I::AbstractVector{T},V::AbstractMatrix) where {T<:Integer} = (length(I)>1 && size(V,1)==1) ? setrow!(X,I,repeat(V,length(I))) : setrow!(X,I,V)
setrowexp!(X::AbstractMatrix,I::AbstractVector{Bool},V::AbstractMatrix)=(size(V, 1)==1 && (N=sum(I);)>0) ? setrow!(X, I, repeat(V, N)) : setrow!(X, I, V)
function setrowexp!(X::Tuple,I,V::Tuple)
    for i=1:length(X)
        setrowexp!(X[i],I,V[i])
    end
    X
end
setrowexp!(X::Dict,I,V::Dict)=dictfun!((x,v)->setrowexp!(x,I,v),X,V, replace=false)
function setrowexp!(X::AbstractVector{T},I::AbstractVector{TI},V::AbstractVector{TV}) where {T<:AbstractArray,TI<:Integer,TV<:AbstractArray}
    X[I]=(length(V)==1 && length(I)>1) ? repeat(V, length(I)) : V
    X
end
function setrowexp!(X::AbstractVector{T},I::AbstractVector{TI},V::AbstractVector{TV}) where {T<:AbstractArray,TI<:Bool,TV<:AbstractArray}
    X[I]=(length(V)==1 && (N=sum(I);)>0) ? repeat(V, N) : V
    X
end
function setrowexp!(X::AbstractVector{T},I::AbstractVector{TI},V::TV) where {T<:AbstractArray,TI<:Integer,TV<:AbstractArray}
    X[I]=fill(V, length(I))
    X
end
function setrowexp!(X::AbstractVector{T},I::AbstractVector{TI},V::TV) where {T<:AbstractArray,TI<:Bool,TV<:AbstractArray}
    X[I]=fill(V, sum(I))
    X
end

# setrowexp!(X::AbstractVector,I,V)=setrow!(X,I,V)
# setrowexp!(X::AbstractVector{Any},I,V)=setrow!(X,I,repmat(V,rownum(I)))
# setrowexp!{T<:AbstractArray}(X::AbstractVector{T},I,V)=setrow!(X,I,repmat(V,rownum(I)))
# setrowexp!{T<:Tuple}(X::AbstractVector{T},I,V)=setrow!(X,I,repmat(V,rownum(I)))
# setrowexp!{T<:Dict}(X::AbstractVector{T},I,V)=setrow!(X,I,repmat(V,rownum(I)))

# function setrowexp!(X::AbstractMatrix,I::@compat(Union{AbstractVector,Range}),V::AbstractMatrix)
#     @assert(size(V,1)==1,"The third input should be a one-row matrix.");
#     setrow!(X,I,repmat(V,length(I),1))
# end
# setrowexp!(X,I,V)=setrow!(X,I,V)

export setrow!, setrowexp!

#[[ rownum rownumsf rnum ]]
#N = rownum[sf](Group)
#Return row number of any input. rnum=rownumsf.
#The rownum of any other type excepting AbstractArray, Tuple and Dict is always 1.
#The rownum of empty Tuple and Dict is always 0.
#When input is Tuple or Dict, rownumsf will check if all elements in the tuple or Dict have the same row number (otherwise a error will occur), but rownum only return the row number of the first element without check if other elements have the same row nubmer.
#See also: getrow, istable, iseqrow, isgroup
#Xiong Jieyi, 14 May 2014> January 25, 2015>May 20, 2015>Jun 16, 2015

export rownum, rownumsf, rnum
rownum(X::AbstractArray)=size(X,1)
rownum(X::Tuple)=isempty(X) ? 0 : rownum(X[1])
rownum(X::Dict)=isempty(X) ? 0 : rownum(first(values(X)))
rownum(X::AbstractDataFrame)=size(X, 1)
rownum(::Any)=1 #Is that OK?
# rownum(X::AbstractString)=1
# rownum(X::Char)=1
# rownum(X::Number)=1
# rownum(X::Symbol)=1
# rownum(X::VersionNumber)=1
function rownumsf(X::Tuple)
    if isempty(X)
        return 0
    end
    rN=rownum(X[1])
    for v in X[2:end]
        @assert(rN==rownumsf(v),"Tuple elements have different row number.")
    end
    rN
end
function rownumsf(X::Dict)
    if isempty(X)
        return 0
    end
    rN=-1
    for v in values(X)
        if rN<0
            rN=rownumsf(v)
        else
            @assert(rN==rownumsf(v),"Dict elements have different row number.")
        end
    end
    rN
end
rownumsf(X)=rownum(X)
const rnum=rownumsf

#}}

#{{ eachr
#[[ eachr asgrp ]]
# Iterators = eachr(Group; as_vec=false)
# Group = asgrp(Iteratiors)
#Convert a group to an iterators, for its each row.
#If group is matrix and as_vec is true, function will vecterize row. This convertion is not recurly.
#For historical reason, my.eachrow(...) is equal to eachr(...). So far my module is still export eachrow() in order to mask the Base.eachrow(). It will be removed in the future.
#See also: rowfun, eachcol
#Xiong Jieyi, May 5, 2015>Mar 10, 2016>30 Mar 2016>8 Apr 2019

mutable struct eachrow
    value::Union{Group,Dict}
    len::Int
    as_vec::Bool
    eachrow(X::Union{Group,Dict};as_vec::Bool=false)=new(X,rownum(X),as_vec)
end
import Base: length, keys
# export eachr, eachrow, length, asgrp #Stopping export eachrow() in the future (to free Base.eachrow().)
export eachr, length, asgrp, keys

eachr(arg...; kw...)=eachrow(arg...; kw...)

start(::eachrow)=1
function next(X::eachrow,state::Int)
    O=getrow(X.value,state)
    (X.as_vec&&isa(O,AbstractMatrix) ? vec(O) : O,state+1)
end
done(X::eachrow,state::Int)=X.len<state
length(X::eachrow)=X.len
keys(X::eachrow)=LinearIndices(1:X.len) #Added in 16 Apr 2019
asgrp(X::eachrow)=vcatr(X...)
#}}

#{{ eachcomb eachcombr
#[[ eachcomb eachcombr ]]
# Iterator = eachcomb(Iterator1, Iterator2, ...)
# Iterator = eachcombr(Group1, Group2, ...) # -- *
#Iterate combination of each element. Each output will be a iterator of tuple.
#* -- Iterate combination of rows. Each iterate is a tuple of rows. Group1 has the inerest cycle. You can use vcatr(eachcombr(...)...) to get a group-structure of combinations.
# eachcomb() has a similar function as Base.Iterators.product(), except collect(eachcomb(...)) always return a vector. For both eachcomb(), eachcombr() and Iterates.product(), the first iterator changes the fastest.
#See also: Iterates.product, eachr, splitdim, splitdimv, vcatr
#Xiong Jieyi, May 14, 2016

if VERSION<v"0.7.0"
    import Base: start, next, done, length
    export eachcombr, start, next, done, length
else
    import Base: length
    export eachcombr, length
end

mutable struct eachcombr
    V::Tuple
    N::Tuple{Vararg{Int}}
    len::Int
    function eachcombr(X::Group...)
        t=map(rnum,X)
        new(X,t,prod(t))
    end
end
start(::eachcombr)=1
next(X::eachcombr,state::Int)=(map(getrow,X.V,ind2sub(X.N,state)), state+1)
done(X::eachcombr,state::Int)=X.len<state
length(X::eachcombr)=X.len

if VERSION<v"0.7.0"
    import Base: start, next, done, length
    export eachcomb, start, next, done, length
else
    import Base: length
    export eachcomb, length
end
mutable struct eachcomb
    V::Tuple
    eachcomb(X...)=new(X)
end
if VERSION<v"0.7.0"
    start(X::eachcomb)=Any[start(v) for v in X.V]
    function next(X::eachcomb,state::Vector{Any})
        t=map(next,X.V,state)
        rlt=tuple(map(x->x[1],t)...)
        ns=map(x->x[2],t)
        p=1
        while p<length(X.V)
            if done(X.V[p],ns[p])
                state[p]=start(X.V[p])
                p+=1
            else
                state[p]=ns[p]
                break
            end
        end
        if p>=length(X.V)
            state[p]=ns[p]
        end
        (rlt,state)
    end
    done(X::eachcomb,state::Vector{Any})=done(X.V[end],state[end])
else
    function iterate(X::eachcomb)
        t=map(iterate, X.V)
        if any(x->x==nothing, t)
            return nothing
        end
        rlt=map(x->x[1], t)
        (rlt, (rlt, map(x->x[2], t)))
    end
    function iterate(X::eachcomb, st::Tuple{Tuple, Tuple})
        rlt=Any[st[1]...]
        state=Any[st[2]...]
        p=1
        while p<=length(X.V)
            t=iterate(X.V[p], state[p])
            if t==nothing
                rlt[p], state[p]=iterate(X.V[p])
                p+=1
            else
                rlt[p], state[p]=t
                break
            end
        end
        if p>length(X.V)
            nothing
        else
            trlt=tuple(rlt...)
            (trlt, (trlt, tuple(state...)))
        end
    end
end
length(X::eachcomb)=prod(map(length,X.V))
#}}

#{{ eachcol [ not be exported ]
# Iterator = eachcol(Matrix)
# Iterate each column of a matrix. Each column will be as a vector in the outputed iterator.
# Since Julia1.0, the Base.eachcol does the similar job as this function, therefore this function is no more exported to avoid confliction. If you need this function, use my.eachcol(...).
#See also: eachrow, splitdim, splitdimv
#Xiong Jieyi, Apr 25, 2016

if VERSION<v"0.7.0"
    import Base: start, next, done, length
    export eachcol, start, next, done, length
else
    import Base: length
    export length
end
mutable struct eachcol
    value::AbstractMatrix
    eachcol(X::AbstractMatrix)=new(X)
end
start(::eachcol)=1
next(X::eachcol,state::Int)=(X.value[:,state], state+1)
done(X::eachcol,state::Int)=size(X.value,2)<state
length(X::eachcol)=size(X.value,2)
#}}

#{{ tuplefun

#[[ tuplefun ]]
#Since the bug in julia map is fixed, this function could be
#replaced by map.
#tuplefun(fun, (X1, X2, ...), (Y1, Y2, ...),...)
#Return (fun(X1,Y1,...), fun(X2,Y2,...), ...)
#See also: map
#Xiong Jieyi, 16 May 2014

tuplefun(fun::Function,T::Tuple)=tuple([fun(x) for x in T]...)
function tuplefun(fun::Function, Ts::Tuple...)
    i_num=length(Ts[1])
    j_num=length(Ts)
    Pc=cell(j_num)
    O=cell(length(Ts[1]))
    for i=1:i_num
        for j=1:j_num
            Pc[j]=Ts[j][i]
        end
        O[i]=fun(Pc...)
    end
    return tuple(O...)
end
export tuplefun

#}}

#{{ vcatr vcatr_with

#[[ vcatr vcatr_with ]]
#Groups = vcatr(Group1, Group2, ...)
#Or Table = vcatr(Table1, Table2, ...; fields=N(default=1)|:intersect|:union|fields_vec, checktable=true)
#vcatr(nothing, else...)=vcatr(else...)
#(labels, Groups) = vcatr_with(Rep, Group1, Group2, ...;...)
#Or Table = vcatr_with("label_field"=>Rep::Group, Table1, Table2, ...;...)
#Catenate the groups. Rep are labels of each group. In the outputs, the labels or label_field has the same row number as catenated group. i.g. vcatr_with(1:3=>"NO", T1, T2, T3) will return a table with new field "NO", which composited by number 1,2,3 assigning the source of each row.
#When input is Dict, fields=:intersect will only combine shared fields, while fields=:union will fill the missing fields as empty value (Unfinished function). fields=N will using the fields of TableN. Or you can also assign the fields as a vector. fields=1 have the farest speed since dictfun only considering fields of the first input by default.
#The nothing or () at the beginning of input will be ignored (but not tolerate the nothing or () in the middle of input).
#
#Note that the tuple input will vcatr recurly. For example:
#vcatr((1,(2,3)), (4,(5,6)))=([1,4],([2,5],[3,6]))
#If you don't want the tuple vcatr recurly, try:
#vcatr((1,[(2,3)]), (4,[(5,6)])) = ([1,4],[(2,3),(5,6)])
#
#See also: rownum, rnum, unslice, reprow, getrow, segrow!, segrow2!, setrowexp, rownum, dictfun, 
#Xiong Jieyi, 8 Jul 2014>15 Sep 2014>December 12, 2014>February 5, 2015>February 8, 2015>May 13, 2016>Jun 5, 2016

export vcatr,vcatr_with
vcatr()=()
vcatr(Xs::Tuple...)=map(vcatr,Xs...)
vcatr(::Tuple{},Xs::Tuple...)=vcatr(Xs...)
vcatr(::Union{Void,Tuple{}},Xs...)=vcatr(Xs...)
vcatr(X...)=vcat(X...)
function vcatr(X1::Dict,Xs::Dict...; fields::Union{Integer,Symbol,Vector}=1)
    Xs=tuple(X1,Xs...)
    if isa(fields,Integer)
        #Since dictfun only considering the fields of the first input, nothing need to do when field=1.
        fds=keys(Xs[fields])
        # Xs=[i==fields?Xs[i]:fd(Xs[i],fds) for i=1:length(Xs)]
    elseif isa(fields,AbstractVector)
        fds=fields
        # Xs=map(x->fd(x,fields),Xs)
    elseif isa(fields,Symbol)
        if fields==:union
            #to be code here... February 5, 2015
            error("Unfinished function.")
        elseif fields==:intersect
            fds=intersect(map(keys,Xs)...)
            # Xs=map(x->fd(x,fds),Xs)
        else
            error("Invalid field $fields.")
        end
    end
    # if checktable
    #     @assert(all(map(istable,Xs)),"Not all input are table.")
    # end
    dictfun(vcatr,Xs...;key=fds)
end
function vcatr_with(rep::Group,X...;wargs...)
    @assert(rownum(rep)==length(X),"Label number ($(length(rep))) is not consistent with group number ($(length(X))).")
    grpnumarr=map(rownum,X)
    lb=vcatr([reprow(getrow(rep,i),grpnumarr[i]) for i=1:length(X)]...)
    return (lb,vcatr(X...;wargs...))
end
function vcatr_with(rep::Pair{T1,T2},args::Dict...;wargs...) where {T1<:AbstractString,T2<:Group}
    t,T=vcatr_with(rep.second,args...;wargs...)
    T[rep.first]=t
    T
end
function vcatr_with(rep::Pair{T1,T2},args::Dict...;wargs...) where {T1<:Group,T2<:AbstractString}
    @warn("vcatr_with(Rep=>\"field\",...) will be replaced by vcatr_with(\"field\"=>Rep,...).")
    vcatr_with(rep.second=>rep.first,args...;wargs...)
end

#}}

#{{ vcatdo vcatdov
#[[ vcatdo vcatdov ]]
# (O1,O2,...) = vcatdo(function, I1, I2, ...)
# T[O1,O2,...] = vcatdov(function, I1, I2, ...) #Faster
# Do as vcatr(O1,O2,...)=function( vcatr(I1,I2,...))
# Each input and output have the same row number. vcatdo return a tuple while vcatdov return a group-container vector.
#See also: asvecdo, splitdim, splitdimv
#Xiong Jieyi, Sep 20, 2015

export vcatdo, vcatdov
function vcatdov(fun::Function,X...)
    L=1:length(X)
    idx,V=vcatr_with(L,X...)
    W=fun(V)
    grpfun(x->typeof(W)[x],fastgrp(idx,L),W)
end
vcatdo(arg...;warg...)=tuple(vcatdov(arg...;warg...)...)
#}}

#{{ sortr sortri uniqr isuniqr unival uninum
#[[ sortr sortri issortedr ]]
#Syntax: Index=sortri(Array{Num|AbstractString,1~2}|Tuple_Group; rev=false)
#(SortedX, Index)=sortr(Array{Num|AbstractString,1~2}|Tuple_Group; rev=false)
#Where getrow(X,Index) is SortedX. For multiple column sort, this algorithm is faster than Base.sortrows.
#See also: uniqr, uniq_sorted_r, issegregatedr
#Xiong Jieyi, 7 May 2014.>February 8, 2015

function sortri(X::AbstractVector{T};rev::Bool=false) where {T<:Union{Number,AbstractString,Symbol,Char,Any}}
    return sortperm(X,rev=rev,alg=length(X)>500 ? TimSort : MergeSort)
end
function sortri(X::AbstractDataFrame; rev::Bool=false)
    return sortperm(X,rev=rev,alg=length(X)>500 ? TimSort : MergeSort)
end

function sortri(X::AbstractMatrix{T};rev::Bool=false) where {T<:Union{Number,AbstractString,Symbol,Char,Any}}
    # S=sortrows([X 1:size(X,1)],rev=rev)
    function sort_core(Ci::Integer,Rg::Union{Range,Vector})
        #idx=Rg[sortperm(X[Rg,Ci],rev=rev,alg=TimSort)]
        idx=Rg[sortri(X[Rg,Ci],rev=rev)]#Changed in February 8, 2015
        if Ci<size(X,2)
            p=1
            for i=2:length(idx)
                if X[idx[p],Ci]!=X[idx[i],Ci]
                    if p+1<i
                        idx[p:i-1]=sort_core(Ci+1,idx[p:i-1])
                    end
                    p=i
                end
            end
            if p<length(idx)
                idx[p:end]=sort_core(Ci+1,idx[p:end])
            end
        end
        return idx
    end
    return sort_core(1,1:size(X,1))
end
function sortri(X::Tuple;rev::Bool=false)
#Solution 1
    # idx=sort([1:rownum(X)],
    #          lt=(ia,ib)->islessr2(X,ia,X,ib))
    # (getrow(X,idx), idx)

#Solution 2
#But seens no obvious speed difference. Both sorting is correct, but the rank is different.
    function sort_core(Ci::Integer,Rg::Union{Range,Vector})
        idx=Rg[sortri(getrow(X[Ci],Rg),rev=rev)]
        if Ci<length(X)
            p=1
            for i=2:length(idx)
                if X[Ci][idx[p],:]!=X[Ci][idx[i],:]
                    if p+1<i
                        idx[p:i-1]=sort_core(Ci+1,idx[p:i-1])
                    end
                    p=i
                end
            end
            if p<length(idx)
                idx[p:end]=sort_core(Ci+1,idx[p:end])
            end
        end
        return idx
    end
    return sort_core(1,1:rownum(X))
end
sortri(X::Dict;rev::Bool=false)=sortri(dict2tuple(X);rev=rev)

function sortr(X::Group;rev::Bool=false)
    I=sortri(X,rev=rev)
    (getrow(X,I), I)
end

#[[ issortedr ]]
#TF=issortedr(X)
#Check if X is row-sorted.
#See also: sortr, sortri, uniqr
#Xiong Jieyi, 15 Jul 2014 > 11 Jun 2019

issortedr(X::AbstractVector)=issorted(X) #To speed up.
function issortedr(X::Group)
    for i=2:rownum(X)
        if islessr2(X,i,X,i-1)
            return false
        end
    end
    return true
end
export sortri, sortr, issortedr

#[[ uniqr unival ]]
#Syntax: (ia, ib)=uniqr(X; stable=false, sorted=false, sortcheck=true)
#        uniqued_X = unival(X; stable=false, sorted=false, sortcheck=true)
#ia is a vector with the length of unique number, and X[ia,:] is the row unique of X, which is also the output of unival(X; stable=false). Note that X[ia,:] is sorted.
#ib is a vector with the length of X_rownum. X[i,:] = uniqued_X[ib[i],:].
#If stable is true, the uniqued_X will in the order of X; otherwise, uniqued_X is sorted.
#If X is sorted, set sorted=true will speeds up. When X is actually not sorted, an error will occur if sortcheck=true, or output unpredictable result if sortcheck=false (DANGEROUS!!).
#See also: isuniqr, rownum, getrow, setrow!
#Xiong Jieyi, Mar 12, 2014>December 11, 2014>Dec 22, 2016>Feb 10, 2017

function uniq_sorted_r(X::Group) #Not be exported since 11 Jun 2019
    if rownum(X)==0  #Bug fixed at 24 Aug 2017
        return (Int[],Int[])
    end
    ia=ones(Int,rownum(X))
    ib=[1]
    p=1
    c=1
    for ri=2:rownum(X)
        if !isequalr2(X,p,X,ri)
            p=ri
            c+=1
            push!(ib,ri)
        end
        ia[ri]=c
    end
    (ib, ia)
end

function uniqr(X::Group; stable::Bool=false, sorted::Bool=false, sortcheck::Bool=true)
    S,Idx=if sorted
        if sortcheck && !issortedr(X)
            error("A is not sorted.")
        end
        (X, 1:rownum(X))
    else
        sortr(X)
    end
    ia,ib=uniq_sorted_r(S)
    if stable
        Oa, Oai=sortr(Idx[ia])
        Ob=rvorder(Oai)[ib[rvorder(Idx)]]
        (Oa, Ob)
    else
        (Idx[ia], ib[rvorder(Idx)])
    end
end

function unival(X::Group; stable::Bool=false, kw...)
    if stable
        getrow(X,falsesbut(rownum(X),uniqr(X; kw...)[1]))
    else
        getrow(X,uniqr(X; kw...)[1])
    end
end
export uniqr, unival

#[[ uninum ]]
# Uniq_row_num=uninum(X; sorted=false, sortcheck=true)
#See also: uniqr, unival, isuniqr
#Xiong Jieyi, 11 Sep 2014
uninum(X::Group; kw...)=length(uniqr(X; kw...)[1])
export uninum

#[[ isuniqr ]]
#T|F = isuniqr(Group; sorted=false, sortcheck=true)
#Return if the Group is unique for all its rows, i.e., no repeat.
#See also: uniqr, uni_sorted_r, rownum, getrow, setrow!, uninum
#Xiong Jieyi, 11 Sep 2014
isuniqr(G::Group; kw...)=rownum(G)==uninum(G; kw...)
export isuniqr

#}}

#{{ issegregatedr
#[[ issegregatedr ]]
# T|F = issegregatedr(Group[, GroupNum])
#Test if a group is segregated, i.e., if the elements of each group are in a continuesrange in the list. If the GroupNum is given, function will not sort and calculated it so it will be faster.
#See also: isuniqr
#Xiong Jieyi, 5 May 2017

export issegregatedr
function issegregatedr(G::Group, GN::Int=uninum(G))
    N=0
    for i=1:rownum(G)-1
        N+=!isequalr2(G,i,G,i+1)
        if N>=GN
            return false
        end
    end
    return true
end
#}}

#{{ intersectri
#[[ intersectri ]]
#(uAi,uBi)=intersectri(A_group, B_group; Asorted=false, Bsorted=false, sortcheck=true)
#All output are vectors of indices. getrow(A_group,uAi) = getrow(B_group,uBi) are the intersect result of A_group and B_group. They are uniqued and in-paired, i.e., only the first of repeated items will be in the output indices. To exhaust all combination between A and B, try setxri.
#See also: uniqr, sortr, vcatr, setxri
#Xiong Jieyi, 4 Sep 2014 >May 31, 2015
export intersectri
function intersectri(A::Group,B::Group; Asorted::Bool=false, Bsorted::Bool=false, sortcheck::Bool=true)
    (A,B)=promote_group(A,B)
    (sA,siA)=if Asorted
        if sortcheck && !issortedr(A)
            error("A is not sorted.")
        end
        (A, 1:rownum(A))
    else
        sortr(A)
    end
    (sB,siB)=if Bsorted
        if sortcheck && !issortedr(B)
            error("B is not sorted.")
        end
        (B, 1:rownum(B))
    else
        sortr(B)
    end
    lenA=length(siA)
    lenB=length(siB)
    if lenA==0
        return (Int[],Int[],Int[],[1:lenB])
    elseif lenB==0
        return (Int[],Int[],[1:lenA],Int[])
    end
    pA=1
    pB=1
    uA=Int[];uB=Int[];
    while true
        if isequalr2(sA,pA,sB,pB)
            push!(uA,pA)
            push!(uB,pB)
            pA+=1
            while pA<=lenA && isequalr2(sA,pA,sA,pA-1)
                pA+=1
            end
            pA>lenA && break
            pB+=1
            pB>lenB && break
        elseif islessr2(sA,pA,sB,pB)
            pA+=1
            pA>lenA && break
        else
            pB+=1
            pB>lenB && break
        end
    end
    return (siA[uA],siB[uB])
end

#}}

#{{ setxri
#[[ setxri ]]
#(uAi,uBi)=setxri(A_group, B_group; uniqB=false, Asorted=false, Bsorted=false, sortcheck=true)
#All output are vectors. uAi and uBi have the same lengths showing the index (row number) pairs of each A and B. A[uAi, :]==B[uBi, :]. Both A[uAi, :] and B[uBi, :] will be sorted.
#If uniqB is ture, only the first item of repeated elements in B will be reported.
#See also: uniqr, sortr, vcatr
#Xiong Jieyi, 4 Sep 2014 > 11 Jun 2019
function setxri(A::Group,B::Group; uniqB::Bool=false, Asorted::Bool=false, Bsorted::Bool=false, sortcheck::Bool=true)
    A,B=promote_group(A,B)
    (sA,siA)=if Asorted
        if sortcheck && !issortedr(A)
            error("A is not sorted.")
        end
        (A, 1:rownum(A))
    else
        sortr(A)
    end
    (sB,siB)=if Bsorted
        if sortcheck && !issortedr(B)
            error("B is not sorted.")
        end
        (B, 1:rownum(B))
    else
        sortr(B)
    end
    lenA=length(siA)
    lenB=length(siB)
    if lenA==0
        return (Int[],Int[],Int[],collect(1:lenB))
    elseif lenB==0
        return (Int[],Int[],collect(1:lenA),Int[])
    end
    pA=1
    opA=1
    npA=2
    pB=1
    islasteq=false
    uA=Int[];uB=Int[];rA=Int[];rB=Int[]
    while true
        if islasteq
            if pA<=lenA && isequalr2(sA,pA,sB,pB)
                push!(uA,pA)
                push!(uB,pB)
                pA+=1
            else
                pB+=1
                if pB>lenB
                    break
                end
                if !uniqB
                    npA=pA
                    pA=opA
                else
                    if pA>lenA
                        break
                    end
                    npA=pA+1
                end
                islasteq=false
            end
        else
            if isequalr2(sA,pA,sB,pB)
                push!(uA,pA)
                push!(uB,pB)
                opA=pA
                islasteq=true
                pA+=1
            elseif islessr2(sA,pA,sB,pB)
                if npA>lenA
                    break
                end
                pA=npA
                npA+=1
            else
                pB+=1
                if pB>lenB
                    break
                end
            end            
        end
    end
    return (siA[uA],siB[uB])
end
# setxri(A::Group,B::Group;args...)=setxri(promote_group(A,B)...;args...)
export setxri

#}}

#{{ memberr ismbr
#[[ memberr ismbr ismemberr ]]
#I = memberr(X, SetGroup; Asorted=false, Bsorted=false, sortcheck=true)
#T|F= ismbr(X, SetGroup; Asorted=false, Bsorted=false, sortcheck=true)
#T|F= ismbr(X, Set(eachrow(SetGroup)), as_vec=false) #Run in hash rather than in sorting algorithm.
#I is Int vector with the length of X, showing the elements index in Set.
#The elements absent in Set will represent as 0.
#ismemberr only return true or false.
#ismbr is an alias of ismemberr.
#See also: setxr, intersectri, setcmpr, eachrow
#Xiong Jieyi, 14 May 2014>3 Sep 2014>23 Jan 2018

function memberr(A::Group, B::Group; Asorted::Bool=false, Bsorted::Bool=false, sortcheck::Bool=true)
    (uA,uB)=setxri(A,B;uniqB=true, Asorted=Asorted, Bsorted=Bsorted, sortcheck=sortcheck)
    O=zeros(Int,rownum(A))
    O[uA]=uB
    return O
end
function ismemberr(X...; kw...)
    # warn("ismemberr will be replaced by ismbr in the furture.")
    memberr(X...; kw...).>0
end
ismbr(X...; kw...)=memberr(X...; kw...).>0
ismbr(A::Group, B::Set; as_vec::Bool=false)=BitVector(map(x->in(x, B), eachrow(A; as_vec=as_vec)))
export memberr, ismemberr, ismbr
#}}

#{{ dtshift, dtshift!
#[[ dtshift dtshift! ]]
#targetV = dtshift(targetI, sourceI, sourceV[, DefaultVal|nothing=emptyval][; safe=false])
#(tV1,tV2,...) = dtshift(targetI, sourceI, (sV1, sV2, ...)[, (defV1, devV2,...)]; ...)
#Dict  = dtshift(targetI, sourceI, Dict[, Dict]; ...)
# ...  = dtshift!(..., defaultVal; ...) #The defaultVal may changed after run, but it is not guaranteed due to the promoting reason. Try replacer!() for in-situ replacement.
#DefaultVal could be:
#  a scalar when sourceV is a vector or matrix (1)
#  a single row as sourceV
#  a matrix with row number as targetI and size as sourceV (2)
#When input is Dict, (1) is not supported.
#In dtshift! function, DefaultVal must be like (2). 
#When defaultVal is nothing, defaultVal will be given by emptyval(getrow(sourceV,1)).
#When safe=true (moderately slower), the ambiguious assignment will trigger an error.
#See also: replacer!, emptyexp, memberr, tbx, grpfun, grpfunexp, reprow
#Xiong Jieyi, 15 May 2014>12 Sep 2014>10 Oct 2014>May 21, 2015>Feb 8, 2017>22 Apr 2017>7 Jun 2019

function dtshift_old(I::Group,VI::Group,V::Union{Group,Dict}) #Only be temperately kepted for checking. Remove it in the future. 11 Jun 2019.
    idx=memberr(I, VI)
    try
        return getrow(V, idx)
    catch err
        if isa(err, BoundsError)
            if rnum(VI)!=rnum(V)
                error("source_index(arg#2) and source_value(arg#3) have different row number.") #Added on 22 Apr 2017
            else
                error("Some target ID is not exist in source ID. Try set default value (nothing = emptyrows).")
            end
        else
            rethrow()
        end
    end
end
function dtshift(I::Group,VI::Group,V::Union{Group,Dict}; safe::Bool=false, Asorted::Bool=false, Bsorted::Bool=false, sortcheck::Bool=true)
    (sVI, siVI)=if Bsorted
        if sortcheck && !issortedr(VI)
            error("B is not sorted.")
        end
        (VI, 1:rownum(VI))
    else
        sortr(VI)
    end
    if Asorted && sortcheck && !issortedr(I)
        error("A is not sorted.")
    end
    idx=memberr(I, sVI; Asorted=Asorted, Bsorted=true, sortcheck=false)
    if safe
        t=try
            getrow(sVI, unique(idx))
        catch err
            if isa(err, BoundsError) && any(isequal(0), idx)
                error("Some target ID is not exist in source ID. Try set default value (nothing = emptyrows).")
            else
                rethrow(err)
            end
        end
        ul=ismbr(sVI, t; Asorted=true, Bsorted=true, sortcheck=false)
        for x in eachgrp(fastgrp(getrow(sVI, ul); sorted=true, sortcheck=false), getrow(V, siVI[ul]))
            isrowsame(x) || error("Ambiguious assignment detected.")
        end
    end
    try
        return getrow(V, siVI[idx])
    catch err
        if isa(err, BoundsError)
            if rnum(VI)!=rnum(V)
                error("source_index(arg#2) and source_value(arg#3) have different row number.") #Added on 22 Apr 2017
            elseif any(isequal(0), idx)
                error("Some target ID is not exist in source ID. Try set default value (nothing = emptyrows).")
            else
                rethrow()
            end
        else
            rethrow()
        end
    end
end
function dtshift!(I::Group,VI::Group,V::Union{Group,Dict},default::Union{Group,Dict},rowcheck::Bool=true; safe::Bool=false, Asorted::Bool=false, Bsorted::Bool=false, sortcheck::Bool=true)
    default,V=promote_anyrow!(default,V)
    # idx=memberr(I,VI)
    # idxl=idx.>0
    # if rowcheck
    #     @assert(rownumsf(default)==rownum(I),"Inproporate row number of default groups.")
    # end
    # if safe
    #     ul=ismbr(VI, getrow(VI, unique(idx[idxl])))
    #     for x in eachgrp(getrow(VI, ul), getrow(V, ul))
    #         isrowsame(x) || error("Ambiguious assignment detected.")
    #     end
    # end
    (sVI, siVI)=if Bsorted
        if sortcheck && !issortedr(VI)
            error("B is not sorted.")
        end
        (VI, 1:rownum(VI))
    else
        sortr(VI)
    end
    if Asorted && sortcheck && !issortedr(I)
        error("A is not sorted.")
    end
    idx=memberr(I, sVI; Asorted=Asorted, Bsorted=true, sortcheck=false)
    idxl=idx.>0
    if safe
        ul=ismbr(sVI, getrow(sVI, unique(idx[idxl])); Asorted=true, Bsorted=true, sortcheck=false)
        for x in eachgrp(fastgrp(getrow(sVI, ul); sorted=true, sortcheck=false), getrow(V, siVI[ul]))
            isrowsame(x) || error("Ambiguious assignment detected.")
        end
    end
    setrow!(default, idxl, getrow(V, siVI[idx[idxl]]))
    return default
end
function _dtshift_core(I::Group,VI::Group,V::Tuple,default::Tuple; kw...) #Add in Feb 8, 2017
    length(V)==length(default) || error(f"Wrong default value tuple length. It should be $1."(length(V)))
    O=Vector{Any}(undef, length(default))
    for i=1:length(default)
        D=default[i]
        if rownum(D)==1
            t=if isa(V[i],Matrix) && !(isa(D,Matrix)) &&
                (eltype(V[i])==typeof(D) || (eltype(V[i])<:Real && isa(D,Real)))
                fill(D,1,size(V[i],2))
            else
                D
            end
            O[i]=reprow(t,rownum(I))
        else
            defrn==rownum(I) || error("Inproporate row number of default groups.")
            O[i]=deepcopy(D)
        end
    end
    return dtshift!(I,VI,V,tuple(O...),false; kw...)
end
function _dtshift_core(I::Group,VI::Group,V::Union{Group,Dict},default::Union{Group,Dict}; kw...)
    defrn=rownum(default)
    if defrn==1
        O=reprow(default,rownum(I))
    else
        @assert(defrn==rownum(I),"Inproporate row number of default groups.")
        O=deepcopy(default)
    end
    return dtshift!(I,VI,V,O,false; kw...)
end
dtshift(I::Group,VI::Group,V::Group,default; kw...)=
    _dtshift_core(I,VI,V,default; kw...)
function dtshift(I::Group,VI::Group,V::Dict,default::Dict; kw...)
    @assert(istable(V) && istable(default),"Input #3 or #4 is not a table.")
    _dtshift_core(I,VI,V,default; kw...)
end
dtshift(I::Group,VI::Group,V::AbstractVector,default::AbstractVector; kw...)=
    _dtshift_core(I,VI,V,default; kw...)
dtshift(I::Group,VI::Group,V::AbstractVector,default; kw...)=
    dtshift!(I,VI,promote_group(V,fill(default,rownum(I)))...,false; kw...)
dtshift(I::Group,VI::Group,V::AbstractMatrix,default::AbstractMatrix; kw...)=
    _dtshift_core(I,VI,V,default; kw...)
dtshift(I::Group,VI::Group,V::AbstractMatrix,default; kw...)=
    dtshift!(I,VI,V,fill(default,rownum(I),size(V,2)),false; kw...)
dtshift(I::Group,VI::Group,V::AbstractVector,::Void; kw...)=
    dtshift!(I,VI,V,emptyrows(V,rownum(I)),false; kw...)
dtshift(I::Group,VI::Group,V::AbstractMatrix,::Void; kw...)=
    dtshift!(I,VI,V,emptyrows(V,rownum(I)),false; kw...)
dtshift(I::Group,VI::Group,V::Group,::Void; kw...)=
    dtshift!(I,VI,V,emptyrows(V,rownum(I)),false; kw...)
function dtshift(I::Group,VI::Group,V::Dict,::Void; kw...)
    @assert(istable(V),"Input #3 is not a table.")
    dtshift!(I,VI,V,emptyrows(V,rownum(I)),false; kw...)
end
dtshift(::Group,::Group,::AbstractVector,::Tuple; kw...)=error("The third and fourth inputs should be both Tuples.")
dtshift(::Group,::Group,::AbstractMatrix,::Tuple; kw...)=error("The third and fourth inputs should be both Tuples.")
dtshift(::Group,::Group,::AbstractVector,::AbstractMatrix; kw...)=error("The third and fourth inputs should be both Matrix or Vector.")
dtshift(::Group,::Group,::AbstractMatrix,::AbstractVector; kw...)=error("The third and fourth inputs should be both Matrix or Vector.")
dtshift(::Group,::Group,::AbstractArray,::AbstractArray; kw...)=error("The third and fourth inputs should be both Matrix or Vector.")
export dtshift, dtshift!
#}}

#{{ replacer!
#[[ replacer! ]]
#    X = replacer!(X::Vector|Matrix|Tuple|Dict, A, B; others=default_value)
#  ... = replacer!(X, A1=>B1, A2=>B2, ...; ...)
# In-situ replace the A-like row in X to B. X and A are homogerious groups. B could be the same row number of A or be one-row.
#See also: dtshift, fastgrp, reprow
#Xiong Jieyi, 28 Aug 2014 >Sep 16, 2015 >Feb 14, 2017 >12 Jun 2017>27 Oct 2019

# # output = my.replacer( X, A, B; new=false|others=default_value) #Unrecommended*.
# # replacer() was abandoned as its function covered by dtshift(). 27 Oct 2019
# # In replacer but not in replacer!, when new=true, the unmatched row in X will be replaced by "others"(or empty value if missed) rather than the value in X. B can be in different type from X or A. Setting others=value implying new=true.
# function replacer(X::Group, A::Group, B::Group; others=nothing, new::Bool=others!=nothing)
#     if rownum(B)==1
#         getB=(x)->getrow(B,1)
#     else
#         getB=(x)->getrow(B,x)
#     end
#     if new
#         if others==nothing
#             others=emptyval(getB(1))
#         end
#         O=reprow(others,rownum(X))
#     else
#         O=copy(X)
#     end
#     fg=fastgrp(X,A)
#     for i=1:fg.grpnum
#         l=want(fg,i)
#         if !isempty(l)
#             setrowexp!(O,l,getB(i));
#         end
#     end
#     return O
# end

function replacer!(X::Union{Array, Tuple, Dict}, A::Group, B::Group; others=nothing)
    if rownum(B)==1
        getB=(x)->getrow(B,1)
    else
        getB=(x)->getrow(B,x)
    end
    fg=fastgrp(X,A)
    if others==nothing
        for i=1:fg.grpnum
            l=want(fg,i)
            if !isempty(l)
                setrowexp!(X,l,getB(i));
            end
        end
    else
        nl=trues(rownum(X))
        for i=1:fg.grpnum
            l=want(fg,i)
            nl[l].=false
            if !isempty(l)
                setrowexp!(X,l,getB(i));
            end
        end
        setrowexp!(X,nl, rownum(others)==1 ? others : getrow(others,nl))
    end
    return X
end
function replacer!(X, P1::Pair, Ps::Pair...; kw...)
    A=vcatr(P1.first, map(x->x.first, Ps)...)
    B=vcatr(P1.second, map(x->x.second, Ps)...)
    replacer!(X, A, B; kw...)
end
export replacer!
#}}

#{{ mapr
#[[ mapr ]]
# rows = mapr(fun, ...; mrows=false, ignore=mrows)
# Run function and collect output in a rowpile by addrow!(mrows=false) or addrows!(mrows=true). When ignore is true, nothing output will be ignored.
# See also: eachrow, eachgrp, eachline
#Xiong Jieyi, Feb 23, 2016>Sep 12, 2016

export mapr
function mapr(fun::Function, Xs...; mrows::Bool=false, ignore::Bool=mrows)
    rlt=rowpile()
    for x in zip(Xs...)
        O=fun(x...)
        if !ignore || O!=nothing
            if mrows
                addrows!(rlt,O)
            else
                addrow!(rlt,O)
            end
        end
    end
    value(rlt)
end
#}}

#{{ fastgrp
#[[ fastgrp ]]
#Obj=fastgrp(X; sorted=false, sortcheck=true)
#Obj=fastgrp(X, Groupid; Asorted=false, Bsorted=false, sortcheck=true)
#Obj.grpid: Get uniqued group ID;
#Obj.grpnum: Get group number;
#Obj.itemnum: Get item number;
#Obj.grpfirsti: The index of the first group index;
#want(Obj,#i): Get the index for the ith group.
#When sorted is ture, function will regarded X as sorted and saving the time of sorting.
#fastgrp is also a Iterator.
#See also: grpfun, grpfunexp, dtshift, hashgrp, mapr
#Xiong Jieyi, 7 May 2014 >Feb 17, 2016>3 Oct 2019

struct fastgrp
    grpid::Group
    grpfirsti::Vector{Int}
    grpnum::Int
    itemnum::Int
    _index::Vector{Int}
    _boundary::Vector{Int}
    _groupmap::Union{Vector{Int},Range}
    function fastgrp(X::Group; sorted::Bool=false, sortcheck::Bool=true, stable::Bool=false)
        if stable && !sorted
            #This code could be optimized to avoid sorting twice. 3 Oct 2019
            return fastgrp(X, unival(X, stable=true))
        end
        if sorted
            if sortcheck && !issortedr(X)
                error("A is not sorted.")
            end
            sX=X
            si=1:rownum(X)
        else
            (sX,si)=sortr(X)
        end
        (iu,id)=uniq_sorted_r(sX)
        grpid=getrow(sX,iu)
        new(grpid,si[iu],length(iu),length(id),si,
            [iu;length(id)+1],1:length(id))
    end
    function fastgrp(X::Group, L::Group; Asorted::Bool=false, Bsorted::Bool=false, sortcheck::Bool=true)
        if Asorted
            if sortcheck && !issortedr(X)
                error("A is not sorted.")
            end
            sX=X
            si=1:rownum(X)
        else
            (sX,si)=sortr(X)
        end
        (iu,id)=uniq_sorted_r(sX)
        grpid=getrow(sX,iu)
        if Bsorted && sortcheck && !issortedr(L)
            error("B is not sorted.")
        end
        mp=memberr(L, grpid; Asorted=Bsorted, Bsorted=true, sortcheck=false)
        #In 26 Sep 2019, I found and fixed a serious bug in above sentense. It was:
        # mp=memberr(L, grpid; Asorted=Asorted, Bsorted=true, sortcheck=false)
        #Which will output wrong result without any error. So far this bug haven't be used in all other functions of my module.

        grpnum=rownum(L)
        fsti=zeros(Int,grpnum)
        l=mp.>0
        fsti[l]=si[iu][mp[l]]
        
        new(L,fsti,grpnum,length(id),si,
            [iu;length(id)+1],mp)
    end
end
function want(O::fastgrp,X::Int)
    I=O._groupmap[X]
    if I==0
        return Int[]
    else
        O._index[O._boundary[I]:O._boundary[I+1]-1]
    end
end
start(::fastgrp)=1
next(o::fastgrp,i::Int)=(want(o,i),i+1)
done(o::fastgrp,i::Int)=o.grpnum<i
length(o::fastgrp)=o.grpnum

import Base: length
export fastgrp, want, length
#}}

#{{ freq freqmulti freqas freqmultias idennum

#[[ freq freqmulti freqas freqmultias idennum ]]
# (Frequence, outGroupId) = freq(GroupId; count=vector, revsort=false,...)
# (Freq_Matrix, outGroupId) = freqmulti(Grp1,Grp2,...; count=(c1,c2,...), ...)
# Freq_Matrix = freqmulti(Grp1,Grp2,...; order=outGroupId, count=(c1,c2,...), ...)
# Freq = idennum(GroupId) #Freq is expanded as the same size of GroupId.
# Frequence|Freq_Matrix = freq|freqmulti(...; ..., order=outGroupId)
#    Equal to freqas|freqmultias(outGroupId,...; ...)
#Calculate the frequece of group.
#freqas should be faster than freq(..., order=...)
# ...multi... functions always return the frequence in matrix, even only one group input.
# revsort=true: Output will be reversely sorted by frequence.
#See also: fastgrp, grpfun, isgrpident, twofactormx
#Xiong Jieyi, 27 Aug 2014>11 Sep 2014>30 Sep 2014>March 27, 2015

export freq, freqmulti, idennum, freqas, freqmultias
function freqas(order::Group,X::Group;count::Union{Void,AbstractVector}=nothing)
    if isa(count,Void)
        (N,id)=__freq_sorted_r(sortr(X)[1])
    else
        (N,id,ignore)=grpfun(sum,X,count)
    end
    dtshift(order,id,N,0)
end
function freqmultias(order::Group,Xs::Group...;count::Union{Void,AbstractVector}=nothing)
    if isempty(Xs)
        return zeros(Int,length(order),0)
    end
    if isa(count,Void)
        N=hcat(map(x->freq(x;order=order),Xs)...)
    else
        N=hcat(map((x,y)->freq(x;order=order,count=y),Xs,count)...)
    end
    N
end
function freq(X::Group;order::Union{Void,Group}=nothing,
              count::Union{Void,AbstractVector}=nothing,
              revsort::Bool=false)
    if isa(count,Void)
        (N,id)=__freq_sorted_r(sortr(X)[1])
    else
        (N,id,ignore)=grpfun(sum,X,count)
    end
    if isa(order,Void)
        if revsort
            (N,l)=sortr(N,rev=true)
            (N, getrow(id,l))
        else
            (N, id)
        end
    else
        dtshift(order,id,N,0)
    end
end
function freqmulti(Xs::Group...;count::Union{Void,Tuple}=nothing,order::Union{Void,Group}=nothing)
    noorder=isa(order,Void)
    if noorder
        order=unival(vcatr(Xs...))
    end
    if isempty(Xs)
        if noorder
            return (zeros(Int,length(order),0), order)
        else
            return zeros(Int,length(order),0)
        end
    end
    if isa(count,Void)
        N=hcat(map(x->freq(x;order=order),Xs)...)
    else
        N=hcat(map((x,y)->freq(x;order=order,count=y),Xs,count)...)
    end
    if noorder
        (N, order)
    else
        N
    end
end
idennum(X::Group)=grpfunexp(length,X,1:rownum(X))

function __freq_sorted_r(X::Group)
    if rownum(X)==0
        return (Int[], zerorow(X))
    end
    grpi=[1]
    frq=Int[]
    c=1
    for ri=2:rownum(X)
        if isequalr2(X,ri-1,X,ri)
            c+=1
        else
            push!(frq,c) #For the pervious group
            push!(grpi,ri) #For the next group
            c=1
        end
    end
    push!(frq,c)
    (frq, getrow(X,grpi))
end

#}}

#{{ grpfun, grpfunwith, grpfunexp, grploop, eachgrp, grpvec

#[[ grpfun grpfunwith grploop eachgrp grpvec ]]
#(rlt, grpid) = grpfun(fun, group, param1, param2,...;
#                          input_grpno|no=false, input_grpid|id=false)
#rlt=grpfun(fun, fastgrp(grp, grp_order),...; default=..., ...) #Group by given order. Note that if function input a fastgrp object rather than a group ID, only one output.
#grpvec(...) = grpfun(x->typeof(x)[x], ...)
#(rlt, grpid)=grpfunwith(fun, group|fastgrp, param1, param2,...;
#                          input_grpno|no=false, input_grpid|id=false, default=...)
#grpid = grploop(fun, group|fastgrp, param1, param2,...;
#                          input_grpno|no=false, input_grpid|id=false)
#Or:
# for ([grpno][, grpid] contents) in eachgrp(group|fastgrp, ...; ...)
#        #do something with contents
# end
#Run function for each group of the inputs.
#If input_grpno is true, the first argument for fun is the group_NO.
#If input_grpid is true, the first argument for fun is the group_ID.
#If both input_grpno and input_grpid is ture, the first and second argument for fun is group_NO and group_ID respectively.
#The differences between grpfun and grpfunwith: grpfun require the output should be one row, while grpfunwith can have any number line of output. Note that the vector output in grpfun is regarded as a row while in grpfunwith, be regarded as a column. grpfunwith is also supported empty output. 
#Only for grpfunwith but not for grpfun, if the given function return nothing, this record will be ignored.
#grploop doesn't collect any output.
#Change in May 29, 2017: The third output (fastgrp-object) in grpfun is cancelled.
#See also: grpfunexp, dtshift, rowfun, dictfun, vcatr, vcatr_with
#Xiong Jieyi, 26 May 2014 >5 Sep 2014>10 Sep 2014>1 Oct 2014>10 Oct 2014>December 12, 2014>Dec 19, 2015>Feb 17, 2016>May 29, 2017

#Discarded parameter: multi_output=false, row_exp=false (only for grpfunwith)
#When multi_output=true, function should return a tuple. Each elements in the tuple will be pileed-up by row. rowexp is only useful while multi_output is true.
#If group is a fastgrp object, the behave of empty group is depend on whether default is given. In the case the defalut is missed, a zerorow group will be inputed to the fun. Otherwise, fun will not be called but just regard the default as output. Therefore, grpfunwith(...;default=nothing) means ignore the empty group.
function grpfun(fun::Function, G::fastgrp, V...; no::Bool=false, id::Bool=false,
                input_grpno::Bool=no, input_grpid::Bool=id, multi_output::Bool=false, default=())
    (isempty(V) || G.itemnum==rownumsf(V[1])) || error("Group number is not equal to the row number of the first input.")
    if G.grpnum==0
        @warn("Group number is zero, return nothing instead.")
        return nothing
    end
    isemptyval=default!=()
    rlt=rowpile()
    for i=1:G.grpnum
        l=want(G,i)
        if isemptyval && isempty(l)
            crlt=default
        else
            crlt=nothing
            try
                if input_grpno && input_grpid
                    crlt=fun(i,getrow(G.grpid,i),getrow(V,l)...)
                elseif input_grpno
                    crlt=fun(i,getrow(V,l)...)
                elseif input_grpid
                    crlt=fun(getrow(G.grpid,i),getrow(V,l)...)
                else
                    crlt=fun(getrow(V,l)...)
                end
            catch err
                println("Error occured in group $i-$(getrow(G.grpid,i))")
                rethrow(err)
            end
        end
        try
            if multi_output && isa(crlt,Tuple)
                addrow!(rlt,crlt...)
            else
                addrow!(rlt,crlt)
            end
        catch err
            println("Error occured in group $i-$(getrow(G.grpid,i)), with function returns:")
            sk(crlt)
            rethrow(err)
        end
    end
    return value(rlt, keep_tuple=multi_output)
end
function grpfun(fun::Function, G::Group, V...; P...)
    fg=fastgrp(G)
    # return (grpfun(fun, fg, V...; P...),fg.grpid,fg)
    return (grpfun(fun, fg, V...; P...),fg.grpid)
end

function grpfunwith(fun::Function, G::fastgrp, V...; no::Bool=false, id::Bool=false,
                 input_grpno::Bool=no, input_grpid::Bool=id, multi_output::Bool=false, rowexp::Bool=false, default=())
    (isempty(V) || G.itemnum==rownumsf(V[1])) || error("Group number is not equal to the row number of the first input.")
    isemptyval=default!=()
    rlt=rowpile()
    for i=1:G.grpnum
        l=want(G,i)
        cid=getrow(G.grpid,i)
        if isemptyval && isempty(l)
            crlt=default
        else
            crlt=nothing
            try
                if input_grpno && input_grpid
                    crlt=fun(i,cid,getrow(V,l)...)
                elseif input_grpno
                    crlt=fun(i,getrow(V,l)...)
                elseif input_grpid
                    crlt=fun(cid,getrow(V,l)...)
                else
                    crlt=fun(getrow(V,l)...)
                end
            catch err
                println("Error occured in group $i-$(getrow(G.grpid,i))")
                rethrow(err)
            end
        end
        try
            if !isa(crlt,Void)
                cids=reprow(cid,rownum(crlt))
                if multi_output && isa(crlt,Tuple)
                    addrows!(rlt,cids,crlt...;rowexp=rowexp)
                else
                    addrows!(rlt,cids,crlt)
                end
            end
        catch err
            println("Error occured in group $i-$(getrow(G.grpid,i)), with function returns:")
            sk(crlt)
            rethrow(err)
        end
        
    end
    if isvirgin(rlt)
        @warn("No any output collected, return nothing instead.")
        return (nothing, zerorow(G.grpid))
    else
        t=value(rlt)
        if multi_output
            return (t[2:end],t[1])
        else
            return (t[2],t[1])
        end
    end
end
function grpfunwith(fun::Function, G::Group, V...; P...)
    fg=fastgrp(G)
    grpfunwith(fun, fg, V...; P...)
end
function grploop(fun::Function, G::fastgrp, V...; no::Bool=false, id::Bool=false,
                 input_grpno::Bool=no, input_grpid::Bool=id)
    (isempty(V) || G.itemnum==rownumsf(V[1])) || error("Group number is not equal to the row number of the first input.")
    for i=1:G.grpnum
        l=want(G,i)
        try
            if input_grpno && input_grpid
                fun(i,getrow(G.grpid,i),getrow(V,l)...)
            elseif input_grpno
                fun(i,getrow(V,l)...)
            elseif input_grpid
                fun(getrow(G.grpid,i),getrow(V,l)...)
            else
                fun(getrow(V,l)...)
            end
        catch err
            println("Error occured in group $i-$(getrow(G.grpid,i))")
            rethrow(err)
        end
    end
    return G.grpid
end
function grploop(fun::Function, G::Group, V...; P...)
    fg=fastgrp(G)
    return grploop(fun, fg, V...; P...)
end
grpvec(args...;wargs...)=grpfun(x->typeof(x)[x],args...;wargs...)
export grpfun, grpfunwith, grploop, grpvec

struct eachgrp
    grp::fastgrp
    contents::Tuple
    input_grpno::Bool
    input_grpid::Bool
    eachgrp(G::fastgrp, V...; no::Bool=false, id::Bool=false,
            input_grpno::Bool=no, input_grpid::Bool=id)=
                new(G,V,input_grpno,input_grpid)
end
eachgrp(G::Group,args...;wargs...)=eachgrp(fastgrp(G),args...;wargs...)
start(::eachgrp)=1
function next(o::eachgrp,i::Int)
    l=want(o.grp,i)
    rlt=if o.input_grpno && o.input_grpid
        (i,getrow(o.grp.grpid,i),getrow(o.contents,l)...)
    elseif o.input_grpno
        (i,getrow(o.contents,l)...)
    elseif o.input_grpid
        (getrow(o.grp.grpid,i),getrow(o.contents,l)...)
    elseif length(o.contents)>1
        getrow(o.contents,l)
    else
        getrow(o.contents,l)[1]
    end
    (rlt,i+1)
end
done(o::eachgrp,i::Int)=o.grp.grpnum<i
length(o::eachgrp)=o.grp.grpnum
if VERSION<v"0.7.0"
    import Base: start, next, done, length
    export eachgrp, start, next, done, length
else
    import Base: length
    export eachgrp, length
end

#}}

#{{ grpfunexp
#[[ grpfunexp ]]
# (output1, output2,...) = grpfunexp(fun, group|fastgrp, v1[, v2...];
#                          input_grpno|no=false, input_grpid|id=false)
#Run function for each group, and expand the result to the same row number of the input.
#If input_grpno is true, the first argument for fun is the group_NO.
#If input_grpid is true, the first argument for fun is the group_ID.
#If both input_grpno and input_grpid is ture, the first and second argument for fun is group_NO and group_ID respectively.
#
#See also: grpfun, grpfunwith, rowfun, dictfun
#Xiong Jieyi, 7 Jul 2014 >5 Sep 2014 > Nov 7, 2016 > 19 Feb 2019 > 16 Sep 2019

#Discarded parameter: multi_output, always_tuple.
#When multi_output=true, function should return a tuple. Each elements in the tuple will be pileed-up by row.
#The differences between grpfun and grpfunwith: grpfun require the output should be one row, while grpfunwith can have any number line of output. Note that the vector output in grpfun is regarded as a row while in grpfunwith, be regarded as a column. grpfunwith is also supported empty output. rowexp is only useful while multi_output is true.
function grpfunexp(fun::Function, G::fastgrp, V...; no::Bool=false, id::Bool=false,
                   input_grpno::Bool=no, input_grpid::Bool=id)
    @assert(G.itemnum==rownumsf(V[1]),
            "Group number is not equal to the row number of first input.");
    rlt=undef
    # isrlttuple=false
    for i=1:G.grpnum
        l=want(G,i)
        crlt=nothing
        try
            if input_grpno && input_grpid
                crlt=fun(i,getrow(G.grpid,i),getrow(V,l)...)
            elseif input_grpno
                crlt=fun(i,getrow(V,l)...)
            elseif input_grpid
                crlt=fun(getrow(G.grpid,i),getrow(V,l)...)
            else
                crlt=fun(getrow(V,l)...)
            end
        catch err
            println("Error occured in group $i-$(getrow(G.grpid,i))")
            rethrow(err)
        end
        try
            if rlt==undef
                rlt=emptyrows(crlt, G.itemnum)
            end
            setrowexp!(rlt,l,crlt)
        catch err
            println("Error occured in group $i-$(getrow(G.grpid,i)), with function returns:")
            display(crlt)
            rethrow(err)
        end
    end
    return rlt
end
function grpfunexp(fun::Function, G, V...; P...)
    fg=fastgrp(G)
    grpfunexp(fun, fg, V...; P...)
end
export grpfunexp
#}}

#{{ hashgrp
#[[ hashgrp ]]
# Dict(Gid => Index) = hashgrp(G::Group)
# [index_of_Gi...] = D[getrow(G,i)]
# Dict(Gid => ValueByGroup) = hashgrp(G::Group, Value)
#Transfer group index as a Dict.
#See also: fastgrp
#Xiong Jieyi, February 26, 2015 > Jan 12, 2017 >Apr 24, 2017

export hashgrp
function hashgrp(G::Group)
    if rownum(G)==0
        return Dict{None,Vector{Int}}()
    end
    D=Dict{typeof(getrow(G,1)),Vector{Int}}()
    for i=1:rownum(G)
        ckey=getrow(G,i)
        if haskey(D,ckey)
            push!(D[ckey],i)
        else
            D[ckey]=[i]
        end
    end
    D
end
function hashgrp(G::Group, X)
    if rownum(G)==0
        return nothing
    end
    D=Dict{typeof(getrow(G,1)), typeof(X)}()
    grploop(G,X,input_grpid=true)do gid, cX
        D[gid]=copy(cX) #copy is VERY IMPORTANT here!! Apr 24, 2017
    end
    D
end
#}}

#{{ rowmx

#[[ rowmx addrow! addrows! ]]
#Using a vector to represent a column-number-fixed matrix.
# obj=rowmx([Type, colnum_num]) #Create a row-matrix objects.
# obj=addrow!(obj, new_row) #Add a new row to the row-matrix. The vector input is also regarded as a row.
# Special case: For addrow!(::Vector{T1}, X::Vector{T2}), length(X) must be one and T1 must be equal to T2 (otherwise an error will occur). Function actually do push!(..., X[1]).
# obj=addrows!(obj, matrix) #Add multiple rows. The vector input is regarded as N x 1 matrix. matrix can also be a empty row.
# obj=addrows!(V::Vector,X::Vector) push!(V,each_of_X).
# N=rownum(obj) #Get row number.
# obj1=pileup!(obj1, obj2, ...) #pile up all matrix to obj1.
# obj0 = pileup(obj1, obj2, ...) #pile up all matrix to obj0.
# mx = value(obj) #Return the matrix.
# Update in 17 Sep 2018: Missing is supported.
# See also: rowpile
#Xiong Jieyi, 7 Jul 2014>December 9, 2014>May 27, 2015>Jun 7, 2015>Jun 16, 2015>Dec 3, 2015>2 Nov 2017

mutable struct rowmx
    vec::Vector
    col_len::Int
    rowmx(T,collen=-1)=new(T[],collen)
    rowmx()=new(Void[],-1)
end
function addrow!(X::rowmx,V::AbstractArray{T,N}) where {T,N}
    if X.col_len==-1
        isempty(V) && error("addrows! can add a empty row, but addrow! can not.")
        if N==1
            X.col_len=length(V)
        else
            #Added in Dec 3, 2015
            size(V,1)==1 || error("Input is neither a vector nor an one-row matrix.")
            X.col_len=size(V,2)
        end
        if isa(X.vec,AbstractArray{Void,1})
            X.vec=T[]
        end
    else
        X.col_len==length(V) || error("Wrong new row length. It should be $(X.col_len).");
    end
    if ( eltype(X.vec)==Missing && !(Missing<:T) ) || ( Missing<:T && !( Missing<:eltype(X.vec) ) )
        X.vec=[X.vec..., V...]
    else
        append!(X.vec, V)
    end
    X
end
function addrow!(X::Dict, V::Dict)
    for (kx, x) in X
        v=V[kx]
        if ( isa(x, AbstractVector{Missing}) && !isa(v, Missing) ) || ( isa(v, Missing) && isa(x, AbstractVector) && !(Missing<:eltype(x)) )
            X[kx]=[x..., v]
        else
            addrow!(X[kx], v)
        end
    end
    X
end
function addrow!(X::AbstractVector{Any},V::AbstractVector{Any})
    length(V)==1 || error("addrow! to a Any-vector, input should be single-element Any-vector. For more general purpose, try push!.")
    push!(X,V[1])
end
#<v0.6# function addrow!{T<:AbstractArray}(X::AbstractVector,V::AbstractVector{T})
function addrow!(X::AbstractVector,V::AbstractVector{T}) where {T<:AbstractArray}
    length(V)==1 || error("addrow! to a vector-of-Array, input should be single-element vector-of-array. For more general purpose, try push!.")
    push!(X,V[1])
end
#<v0.6# function addrow!{T<:Tuple}(X::AbstractVector,V::AbstractVector{T})
function addrow!(X::AbstractVector,V::AbstractVector{T}) where {T<:Tuple}
    length(V)==1 || error("addrow! to a vector-of-Tuple, input should be single-element vector-of-Tuple. For more general purpose, try push!.")
    push!(X,V[1])
end
#<v0.6# function addrow!{T<:Dict}(X::AbstractVector,V::AbstractVector{T})
function addrow!(X::AbstractVector,V::AbstractVector{T}) where {T<:Dict}
    length(V)==1 || error("addrow! to a vector-of-Dict, input should be single-element vector-of-Dict. For more general purpose, try push!.")
    push!(X,V[1])
end
function addrow!(X::AbstractVector,V::AbstractVector)
    (eltype(V)<:eltype(X) && length(V)==1) || error("Invalid type or input is not a one-element vector")
    push!(X,V[1])
end
addrow!(X::AbstractVector,V)=push!(X,V)
function addrow!(X::Tuple,V::Tuple)
    for (x,v) in zip(X,V)
        addrow!(x,v)
    end
    X
end
function addrow!(X::AbstractDataFrame,V::AbstractDataFrame)
    size(V, 1)==1 || error("Second input is not a one-row DataFrame. Try addrows! instead.")
    append!(X, V)
end

function addrows!(X::rowmx,V::AbstractArray{T}) where {T}
    if X.col_len==-1#Added in December 9, 2014
        X.vec=T[]
        X.col_len=size(V,2)
    end
    for ri=1:size(V,1)
        addrow!(X,getrow(V,ri))
    end
    X
end
addrows!(X::rowmx,V)=addrow!(X,V)
# function addrows!(X::AbstractVector,V::AbstractVector)
#     for ri=1:length(V)
#         push!(X,V[ri])
#     end
#     X
# end
addrows!(X::AbstractVector,V::AbstractVector)=append!(X,V)
function addrows!(X::Tuple,V::Tuple)
    for (x,v) in zip(X,V)
        addrows!(x,v)
    end
    X
end
# addrows!(X::Dict,V::Dict)=dictfun!(addrows!,X,V)
function addrows!(X::Dict, V::Dict)
    for (kx, x) in X
        v=V[kx]
        if ( isa(x, AbstractVector{Missing}) && eltype(v)!=Missing ) || ( isa(x, AbstractVector) && !(Missing<:eltype(x)) && Missing<:eltype(v) )
            X[kx]=isa(v, Missing) ? [x..., v] : [x..., v...]
        else
            addrows!(X[kx], v)
        end
    end
    X
end
addrows!(X::AbstractVector,V)=push!(X,V)
addrows!(X::AbstractDataFrame,V::AbstractDataFrame)=append!(X, V)
function pileup!(X::rowmx,Y::rowmx)
    X.col_len==Y.col_len || error("Column number is not identical.")
    X.vec=[X.vec,Y.vec]
    X
end
function pileup!(X::rowmx,Y::rowmx,args::rowmx...)
    pileup!(X,Y)
    pileup!(X,args...)
end
function pileup(X::rowmx,args::rowmx...)
    O=deepcopy(X)
    pileup!(O,args...)
    return O
end
function value(X::rowmx)
    if X.col_len==-1
        return []
    end
    typ=eltype(X.vec)
    cN=X.col_len
    rN=int(length(X.vec)/X.col_len)
    
    #For a bug in Julia0.7.0 and Julia 1.0.0
    if (VERSION==v"0.7.0" || VERSION==v"1.0.0") && Missing<:eltype(X.vec)
        X.vec=copy(X.vec)
    end
    
    if typ==Any || typ<:AbstractArray  #To avoid transpose recurly.
        if cN==1
            return reshape(X.vec,:,1)
        else
            # O=typ==Any?cell(rN,cN):fill(typ[],rN,cN)
            O=Array{typ}(undef, rN, cN) #Changed in Feb 20, 2017 > 17 Dec 2018
            for ri=1:rN, ci=1:cN
                O[ri,ci]=X.vec[(ri-1)*cN+ci]
            end
            return O
        end
    else
        # return reshape(X.vec,(cN,rN))' #Avoid deprecated warning in Julia v0.5 when X.vec is String.
        return permutedims(reshape(X.vec,(cN,rN)),[2,1])
    end
end
rownum(X::rowmx)=length(X.vec)/X.col_len
export rowmx,addrow!,addrows!,value,rownum,pileup!,pileup

#}}

#{{ rowpile

#[[ rowpile add! value rownum elem_num isvirgin pileup pileup! ]]
# rp=rowpile() #Create a rowpile
# rp=addrow!(rp, V1, V2, ...) #Add value to rowpile. Input could be either scalar or row-matrix. If you hope to add a any value or any array to rowpile, make V as single-element Any or Array vector.
# rp=addrows!(rp, V1, V2, ...[, rowexp=false]) #Add multiple lines or empty lines. When rowexp is true, all the single row input will be expanded to the same row number of the first input.
# (V1,V2,...)=value(rp; keep_tuple=false, default=nothing) #Get the value out. When the 2nd argument is false and elem_num(rp)==1, value output the value directly rather than a tuple. When rp is empty, return default value. For compatible reason, syntax value(rp, keep_tuple) is still allowed but not recommended.
# rownum(rp) #Get the row number. Return 0 for the virginal rowpile.
# elem_num(rp) #Get the element numbar. Return 0 for the vinginal rowpile.
# T|F=isvirgin(rp) #Test if rowpile is vinginal.
# pileup!(obj1, obj2, ...) #pile up all matrix to obj1.
# obj0 = pileup(obj1, obj2, ...) #pile up all matrix to obj0.
#When X is not a vector, addrow!(rp,(...)[X]) is equal to addrow!(rp,X).
#Update in Jun 13, 2015: Bool input will be saved as BitArray.
#Update in 17 Sep 2018: Missing is supported.
#See also: rowmx
#Xiong Jieyi, 8 Jul 2014>9 Sep 2014>November 23, 2014>Jun 7, 2015>Jun 16, 2015>Dec 6, 2016>29 Sep 2017

export rowpile,add!,value,rownum,elem_num,isvirgin,pileup,pileup!
mutable struct rowpile
    _container::Vector{Any}
    rowpile()=new(Any[])
end
function __createrp__(X::AbstractMatrix)
    Y=rowmx()
    addrow!(Y,X)
end
function __createrp__(X::AbstractVector)
    if length(X)==1
        deepcopy(X)
    else
        isempty(X) && error("addrow! is not support empty input.")
        __createrp__(trsp(X))
    end
end
__createrp__(X::Range)=__createrp__(collect(X))
__createrp__(X::AbstractArray)=error("addrow! only support 1 or 2 dimension array.")
__createrp__(X::Dict)=dictfun(__createrp__,X)
__createrp__(X::Tuple)=map(__createrp__,X)
__createrp__(X::Bool)=BitArray([X])
__createrp__(X::AbstractDataFrame)=deepcopy(X)
__createrp__(X)=[deepcopy(X)]
function addrow!(R::rowpile,V...)
    if isempty(R._container)
        R._container=Any[__createrp__(v) for v in V]
    else
        length(R._container)==length(V) || error("Inproporate input number.")
        for i=1:length(V)
            if ( isa(R._container[i], AbstractVector{Missing}) && !isa(V[i], Missing) ) || ( isa(V[i], Missing) && isa(R._container[i], AbstractVector) && !(Missing<:eltype(R._container[i])) )
                !isa(V[i], AbstractArray) || length(V[i])==1 || error("Input can only be a element or an one-element-vector.")
                R._container[i]=[R._container[i]..., V[i]]
            else
                addrow!(R._container[i], V[i])
            end
        end
    end
    R
end
function __createrprows__(X::AbstractMatrix)
    Y=rowmx()
    addrows!(Y,X)
end
__createrprows__(X::AbstractVector)=deepcopy(X)
__createrprows__(X::AbstractArray)=error("addrow! only support 1 or 2 dimension array.")
__createrprows__(X::Dict)=dictfun(__createrprows__,X)
__createrprows__(X::Tuple)=map(__createrprows__,X)
__createrprows__(X::Bool)=BitArray([X])
__createrprows__(X::Range)=collect(X)
__createrprows__(X::AbstractDataFrame)=deepcopy(X)
__createrprows__(X)=[deepcopy(X)]

function addrows!(R::rowpile,V...;rowexp::Bool=false)
    rowlen=rnum(V[1])
    if rowexp
        V=collect(V)
        for i=2:length(V)
            vN=rnum(V[i])
            if vN==1
                V[i]=reprow(V[i],rowlen)
            else
                vN==rowlen || error("Inconsistent row lengths.")
            end
        end
    else
        for cV in V[2:end]
            rnum(cV)==rowlen || error("Inconsistent row lengths.")
        end
    end
    if isempty(R._container)
        R._container=Any[__createrprows__(v) for v in V]
    else
        length(R._container)==length(V) || error("Inproporate input number.")
        for i=1:length(V)
            if ( isa(R._container[i], AbstractVector{Missing}) && eltype(V[i])!=Missing ) || ( isa(R._container[i], AbstractVector) && !(Missing<:eltype(R._container[i])) && Missing<:eltype(V[i]) )
                R._container[i]=[R._container[i]..., V[i]...]
            else
                addrows!(R._container[i],V[i])
            end
        end
    end
    R
end

function value(R::rowpile; keep_tuple::Bool=false, default=nothing)
    if isempty(R._container)
        return default
    end
    rlt=tuple(map(__value__, R._container)...)
    if !keep_tuple && length(R._container)==1
        return rlt[1]
    else
        return rlt
    end
end
value(R::rowpile, keep_tuple::Bool)=value(R, keep_tuple=keep_tuple) #Only for compatible. 29 Sep 2017
__value__(X::Dict)=dictfun(__value__,X)
__value__(X::Tuple)=map(__value__,X)
__value__(X::rowmx)=value(X)
__value__(X)=X

rownum(R::rowpile)=isempty(R._container) ? 0 : rownum(R._container[1])
isvirgin(R::rowpile)=isempty(R._container)
elem_num(R::rowpile)=length(R._container)
function pileup!(A::rowpile,B::rowpile)
    length(A._container)==length(B._container) || error("Input number is not identical.")
    for i=1:length(A._container)
        if isa(A._container[i],rowmx)
            pileup!(A._container[i],B._container[i])
        else
            A._container[i]=[A._container[i],B._container[i]]
        end
    end
    A
end
function pileup!(A::rowpile,B::rowpile,args::rowpile...)
    pileup!(A,B)
    pileup!(A,args...)
end
function pileup(X::rowpile,args::rowpile...)
    O=deepcopy(X)
    pileup!(O,args...)
    return O
end

#}}

#{{ emptyval zerorow emptyrows
#[[ emptyval zerorow emptyrows ]]
#empty_val = emptyval|zerorow(Data|Type|(Type1, Type2, ...))
#Or     ...= emptyval|zerorow(Vector|Matrix|Tuple|Dict)
#Or     zerorow(Type[, col_num]) Return a 0-row vector (col_num omitted) or matrix.
#empty_rows = emptyrows(Type|Group, row_num)
#Return empty value for given type.
#emptyval(x::Tuple) will output emptyval every element of the tuple. The tuple could be nested.
#emptyval(x::Dict) will output emptyval every field of the Dict. The Dict could be nested.
#zeorrow and emptyrows are restrict row-functions. The first input of emptyrows must be a Group or Type. e.g. emptyrows([1], 4) is fine but emptyrows(1, 4) is not accepted.
#See also: grpfun, grpfunexp, reprow, rowpile
#Xiong Jieyi, 16 May 2014>January 14, 2015>May 23, 2015>Jan 1, 2016>Nov 7, 2016>13 Mar 2018>19 Feb 2019

export emptyval, zerorow, emptyrows
emptyval(::Type{T}) where {T<:Integer} =zero(T) #Include Bool=false
emptyval(::Type{T}) where {T<:Real} =convert(T,NaN)
emptyval(::Type{Char})='\0'
emptyval(::Type{T}) where {T<:AbstractString} =""
emptyval(::Type{Symbol})=Symbol("")
emptyval(::Tuple{})=()
emptyval(T::Tuple{Vararg{Type}})=map(emptyval,T) #Allow nested group
emptyval(T::Tuple)=map(emptyval,T)
# emptyval(::Type{T}) where {T>:Missing} = missing
emptyval(X::Type)=error("Not support empty value $X.")
emptyval(X::Array{T}) where {T} =fill(emptyval(T), size(X))
emptyval(X::Array{T}) where {T>:Missing} = convert(typeof(X), fill(missing, size(X)))
emptyval(X::BitArray)=falses(size(X))
emptyval(::Range{T}) where {T} =zero(T):zero(T)
function emptyval(X::Array{Any})
    X=Array{Any}(undef, size(X))
    X[:].=nothing
    X
end
emptyval(X::Array{T}) where {T<:Tuple} = convert(Array{Tuple}, fill(tuple(), size(X)))
emptyval(X::Array{T}) where {T<:Union{AbstractArray,Dict}} =fill(T(),size(X))
emptyval(X)=emptyval(typeof(X)) #If input is not a type.
emptyval(D::Dict)=dictfun(emptyval,D)
emptyval(D::AbstractDataFrame)=dataframefun(emptyval, D)

zerorow(X::Type)=fill(emptyval(X),0)
zerorow(X::Type,N::Real)=fill(emptyval(X),0,N)
zerorow(X,N::Real)=fill(X,0,N)
zerorow(X::Array{T}) where {T} =fill(emptyval(T),0,size(X,2))
zerorow(X::Array{T}) where {T>:Missing} = convert(typeof(X), fill(missing, 0, size(X,2)))
zerorow(::BitVector)=falses(0)
zerorow(X::BitMatrix)=falses(0,size(X,2))
zerorow(X::Range{T}) where {T} =zero(T):(zero(T)-1)
zerorow(T::Tuple{Vararg{Type}})=map(zerorow,T) #Allow nested group
zerorow(X::Tuple)=map(zerorow,X) #Allow nested group
zerorow(D::Dict)=dictfun(zerorow,D)
zerorow(D::AbstractDataFrame)=dataframefun(zerorow, D)
zerorow(::Tuple{Vararg{Type}},::Real)=error("Column number can not be assigned for a non-singleton value.")
zerorow(X::Group,N::Real)=error("Column number can not be assigned for a non-singleton value.")
zerorow(X)=fill(X,0)

emptyrows(::Type{T}, rn::Int) where {T<:Integer} =zeros(T, rn)
emptyrows(::Type{T}, rn::Int) where {T<:Real} =reprow(convert(T,NaN), rn)
emptyrows(::Type{Char}, rn::Int)=fill('\0', rn)
emptyrows(::Type{T}, rn::Int) where {T<:AbstractString} =fill("", rn)
emptyrows(::Type{Symbol}, rn::Int)=fill(Symbol(""), rn)
emptyrows(::Tuple{}, rn::Int)=()
emptyrows(T::Tuple{Vararg{Type}}, rn::Int)=map(x->emptyrows(x, rn), T) #Allow nested group
emptyrows(T::Tuple, rn::Int)=map(x->emptyrows(x, rn), T)
emptyrows(X::Type, rn::Int)=error("Not support empty value $X.")
emptyrows(X::AbstractVector{T}, rn::Int) where {T} =fill(emptyval(T), rn)
emptyrows(X::AbstractVector{T}, rn::Int) where {T>:Missing} =convert(typeof(X), fill(missing, rn))
emptyrows(X::AbstractVector{T}, rn::Int) where {T<:AbstractArray} =fill(eltype(T)[], rn)
emptyrows(X::AbstractVector{T}, rn::Int) where {T<:Tuple} =convert(Vector{Tuple}, fill(tuple(), rn))
emptyrows(X::AbstractVector{Any}, rn::Int)=Vector{Any}(undef, rn)
emptyrows(X::AbstractMatrix{T}, rn::Int) where {T} =fill(emptyval(T), (rn, size(X,2)))
emptyrows(X::AbstractMatrix{T}, rn::Int) where {T>:Missing} =convert(typeof(X), fill(missing, (rn, size(X,2))))
emptyrows(X::AbstractMatrix{T}, rn::Int) where {T<:AbstractArray} =fill(eltype(T)[], (rn, size(X,2)))
emptyrows(X::AbstractMatrix{T}, rn::Int) where {T<:Tuple} =convert(Matrix{Tuple}, fill(tuple(), (rn, size(X,2))))
emptyrows(X::AbstractMatrix{Any}, rn::Int)=Matrix{Any}(undef, rn, size(X,2))
emptyrows(X::BitVector, rn::Int)=falses(rn)
emptyrows(X::BitMatrix, rn::Int)=falses(rn, size(X, 2))
emptyrows(D::Dict, rn::Int)=dictfun(x->emptyrows(x, rn), D)
emptyrows(D::AbstractDataFrame, rn::Int)=dataframefun(x->emptyrows(x, rn), D)
emptyrows(X, rn::Int)=fill(emptyval(getrow(X, 1)), rn)
#}}

#{{ reprow
#[[ reprow ]]
# row_pile=reprow(Scalar|Row_Mx|Table, row_num)
# Expand the row to matrix or vector. Input should be a row.
#See also: vcatr, vcatr_with, emptyval, dtshift
#Xiong Jieyi, 2 Sep 2014>January 14, 2015>17 Sep 2019

export reprow
#reprow(X::AbstractVector,N::Real)=repmat(X',N)
function reprow(X::AbstractArray,N::Real)
    @assert(size(X,1)==1,"Input is not a single-row matrix.")
    # repmat(X,N)
    repeat(X, outer=(N, 1))
end
reprow(X::Tuple,N::Real)=map(x->reprow(x,N),X)
reprow(X::Dict,N::Real)=dictfun(x->reprow(x,N),X)
reprow(X,N::Real)=fill(X,N)
#}}

#{{ col2vec
#[[ col2vec ]]
# ... = col2vec(Vector|Matrix|Dict|Tuple)
#Convert single-column matrix to a column vector.
#See also: tidytype, getrow, rownum, setrow!, setrow2!, vcatr
#Xiong Jieyi, 27 Aug 2014

col2vec(X::AbstractMatrix)=size(X,2)==1 ? X[:] : X;
col2vec(X::Dict)=dictfun(col2vec,X)
col2vec(X::Tuple)=map(col2vec,X)
col2vec(X)=X
export col2vec

#}}

#{{ rowfun rowfunwith rowloop
#[[ rowfun rowfunwith rowloop ]]
# (Array1,Array2,...)=rowfun(function, arg1, arg2,...; as_vec=false, input_rowno|no=false)
# ((Array1,Array2,...), row_NO)=rowfunwith(function, arg1, arg2,...; as_vec=false, input_rowno|no=false)
# rownum=rowloop(function, arg1, arg2,...; as_vec=false, input_rowno|no=false)
# Run function(arg1_rowi,arg2_rowi,...) for each row of input and pile-up the results.
# If input is a vector, funciton will recieve the element of the vector. If input is a matrix, function will receive a 1xN matrix(as_vec=false) or a vector(as_vec=true).
# If as_vec is true, the row will be convert matrix to vector before input to function. This function is not recurly.
# If the input_rowno is true, the first argument for given function is row number.
# The differences between rowfun and rowfunwith: rowfun require the output should be one row, while rowfunwith can have any number line of output. Note that the vector output in rowfun is regarded as a row, while in rowfunwith, be regarded as a column. rowfunwith is also supported empty output.
# Only for rowfunwith but not for rowfun, if the given function return nothing, this record will be ignored.
# rowloop() just run the function but do not record any result.
#
#See also: grpfun, grpfunwith, grpfunexp, filefun, filefunwith, fileloop, eachrow, mapr
#Xiong Jieyi, 5 Sep 2014 >October 24, 2014>December 12, 2014>Feb 10, 2016

#Discarded parameter: rowfun(...; ..., multi_output=false), rowfunwith(...; ..., multi_output=false)
#When multi_output=true, function should return a tuple. Each elements in the tuple will be pileed-up by row.
#Discarded parameter: rowfunwith(...; ..., rowexp=false)
#When rowexp=true and multi_output=true, the other outputs of given function will be row-expanded to the row number of the first output before being piped up. Note it is only work when multi_output=true.

function rowfun(fun::Function, args::Union{Group,Dict}...;
                as_vec::Bool=false, no::Bool=false, input_rowno::Bool=no, multi_output::Bool=false)
    if rownum(args)==0
        @warn("Row number is zero, return nothing instead.")
        return nothing
    end
    rlt=rowpile()
    if as_vec
        ff=(x,ri)->isa(x,Matrix) ? vec(x[ri,:]) : getrow(x,ri)
    else
        ff=getrow
    end
    for ri=1:rownumsf(args)
        crlt=nothing
        try
            if input_rowno
                crlt=fun(ri,map(x->ff(x,ri),args)...)
            else
                crlt=fun(map(x->ff(x,ri),args)...)
            end
        catch err
            println("Error occured in raw $ri.")
            rethrow(err)
        end
        try
            if multi_output
                addrow!(rlt,crlt...)
            else
                addrow!(rlt,crlt)
            end
        catch err
            println("Error occured in raw $ri, with function returns:")
            display(crlt)
            rethrow(err)
        end
    end
    return value(rlt, keep_tuple=multi_output)
end
function rowfunwith(fun::Function, args::Union{Group,Dict}...; as_vec::Bool=false,
                no::Bool=false, input_rowno::Bool=no, multi_output::Bool=false, rowexp::Bool=false)
    rlt=rowpile()
    # ff=as_vec?vec:(x->x)
    ff(x)=x
    if as_vec
        ff(x::Matrix)=vec(x)
    end
    for ri=1:rownumsf(args)
        crlt=nothing
        try
            if input_rowno
                crlt=fun(ri,map(x->ff(getrow(x,ri)),args)...)
            else
                crlt=fun(map(x->ff(getrow(x,ri)),args)...)
            end
        catch err
            println("Error occured in row $ri.")
            rethrow(err)
        end
        try
            if !isa(crlt,Void)
                rowno=fill(ri,rownum(crlt))
                if multi_output
                    addrows!(rlt,rowno,crlt...;rowexp=rowexp)
                else
                    addrows!(rlt,rowno,crlt)
                end
            end
        catch err
            println("Error occured in row $ri, with function returns:")
            display(crlt)
            rethrow(err)
        end
    end
    if isvirgin(rlt)
        @warn("No any output collected, return nothing instead.")
        return (nothing,Int[])
    else
        t=value(rlt)
        if multi_output
            return (t[2:end],t[1])
        else
            return (t[2],t[1])
        end
    end
end
function rowloop(fun::Function, args::Union{Group,Dict}...;
                 as_vec::Bool=false, no::Bool=false, input_rowno::Bool=no)
    rn=rownumsf(args)
    ff(x)=x
    if as_vec
        ff(x::Matrix)=vec(x)
    end
    for ri=1:rn
        try
            if input_rowno
                fun(ri,map(x->ff(getrow(x,ri)),args)...)
            else
                fun(map(x->ff(getrow(x,ri)),args)...)
            end
        catch err
            println("Error occured in row $ri.")
            rethrow(err)
        end
    end
    return rn
end
export rowfun, rowfunwith, rowloop

#}}

#{{ isrowsame
#[[ isrowsame ]]
# T|F = isrowsame(X::Group|Dict)
#Check if all the rows are in the same values. NaNs are considered as the same, and -0.0 is not equal to 0.0.
#For Dict input, each field will be checked. Table-style is not required.
#See also: isuniqr, isgrpident
#Xiong Jieyi, May 28, 2015 > Apr 26, 2017

export isrowsame
function isrowsame(X::Group)
    # rN=rownum(X)
    # if rN<=1
    #     return true
    # end
    # row1=getrow(X,1)
    # for crow in drop(eachrow(X),1)
    #     # if row1!=crow
    #     if !isequal(row1,crow)
    #         return false
    #     end
    # end
    # return true
    for i=2:rownum(X)
        if !isequalr2(X,1,X,i)
            return false
        end
    end
    return true
end
function isrowsame(X::Dict)
    for (k,v) in X
        if !isrowsame(v)
            return false
        end
    end
    return true
end

#}}

#{{ isgrpident
#[[ isgrpident ]]
# T|F = isgrpident(G::Group, V::Union{Group, Dict}; eq=false)
#Check if all the V in each group are in the same values.
#When eq=true, function return T|F = isgrpident(A, B) && isgrpident(B, A).
#See also: isuniqr, idennum, grpfun
#Xiong Jieyi, May 28, 2015 > 24 Jan 2018

export isgrpident
function isgrpident(G::Group,X::Union{Group,Dict}; eq::Bool=false)
    flag=true
    grploop(G,X) do cX
        if !isrowsame(cX)
            flag=false
            return nothing
        end
    end
    flag && (!eq || isgrpident(X, G))
end
#}}

#{{ setcmpr
#[[ setcmpr ]]
# setcmpr(A::Group, B::Group; venn=false, label=c"Set A,Set B")
#Print the shared relationship of Row-sets A and B. It has no output value.
#If venn is true, function will also draw a venn figure using matplotlib_venn python package.
#See also: ls, quantilels, eachrow, intersectri, ismemberr
#Xiong Jieyi, May 31, 2015 > Jan 14, 2016

export setcmpr
function setcmpr(A::Group, B::Group; venn::Bool=false, label=c"Set A,Set B")

    #Added in 22 Sep 2018, only for user error prevention.
    if xor(typeof(A)==Any, typeof(B)==Any)
        @warn("Only one of the inputs is Any[].")
    end

    A=Set(eachrow(A))
    Al=length(A)
    @printf("A uni-num: %d\n", Al)
    B=Set(eachrow(B))
    Ao=length(setdiff(A,B))
    @printf("A only: %d (%f%%)\n", Ao, 100*Ao/Al)
    AB=length(intersect(A,B))
    Bl=length(B)
    @printf("AB shared: %d (%f%% in A, %f%% in B)\n", AB, 100*AB/Al, 100*AB/Bl)
    Bo=length(setdiff(B,A))
    @printf("B only: %d (%f%%)\n", Bo, 100*Bo/Bl)
    @printf("B uni-num: %d\n", length(B))

    if venn
        ob=importpy(:matplotlib_venn)
        ob.venn2(subsets=ds("10"=>Ao, "01"=>Bo, "11"=>AB), set_labels=label)
    else
        nothing
    end
end
#}}

#{{ asvecdo
#[[ asvecdo ]]
# Matrix = asvecdo(fun, Matrix...)
# The inputed matrix will be convert as vector then as arguments of function. The result will convert back to Matrix. If the function output is a Tuple, every vector in the tuple will be converted.
#See also: dtshift, vcatdo
#Xiong Jieyi, Sep 18, 2015

export asvecdo
function asvecdo(fun::Function,Xs::AbstractMatrix...)
    Xv=map(vec,Xs)
    O=fun(Xv...)
    sz=size(Xs[1])
    if isa(O,Tuple)
        map(O) do x
            @assert(isa(x,AbstractVector),"Function output is not a vector.")
            reshape(x,sz)
        end
    else
        @assert(isa(O,AbstractVector),"Function output is not a vector.")
        reshape(O,sz)
    end
end
#}}
