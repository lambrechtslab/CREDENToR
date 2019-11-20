# Youtiao module is designed to facilitate biology data analysis in general purpose. It needs to be cleaned specifically for CREDENToR pipeline in the future.

module Youtiao
using SortingAlgorithms

# Loading DataFrames is too slow.
#using DataFrames
#export DataFrame
struct DataFrame; end
struct AbstractDataFrame; end

if VERSION==v"1.1.0"
    #Fix a bug of Julia v1.1.0
    import Base.string
    string(x::Vector{String})=replace(string(convert(Vector{Any}, x)), r"^Any"=>"")
    export string
end

using Printf
using Distributed
using DelimitedFiles
using Statistics
using Random
using Dates

export randperm

#Recover some abolished functions:
const Void = Nothing
const Range = AbstractRange
const find = Base.findall
eval(m, x) = Core.eval(m, x)
ismatch(R::Regex, S::AbstractString) = Base.occursin(R, S)
contains(A::AbstractString, B::AbstractString) = occursin(B, A)
search(s::AbstractString, t::AbstractString) = something(findfirst(t, s), 0:-1)
search(s::AbstractString, t::AbstractString, i::Integer) = something(findnext(t, s, i), 0:-1)
rsearch(s::AbstractString, t::AbstractString) = something(findlast(t, s), 0:-1)
rsearch(s::AbstractString, t::AbstractString, i::Integer) = something(findprev(t, s, i), 0:-1)

mean(X::AbstractArray, n::Int) = Statistics.mean(X, dims=n)
mean(A...; B...) = Statistics.mean(A...; B...)

sum(X::AbstractArray, n::Int) = Base.sum(X, dims=n)
sum(A...; B...) = Base.sum(A...; B...)

median(X::AbstractArray, n::Int) = Statistics.median(X, dims=n)
median(A...; B...) = Statistics.median(A...; B...)

var(X::AbstractArray, n::Int) = Statistics.var(X, dims=n)
var(A...; B...) = Statistics.var(A...; B...)

std(X::AbstractArray, n::Int) = Statistics.std(X, dims=n)
std(A...; B...) = Statistics.std(A...; B...)

maximum(X::AbstractArray, n::Int) = Base.maximum(X, dims=n)
maximum(A...; B...) = Base.maximum(A...; B...)

minimum(X::AbstractArray, n::Int) = Base.minimum(X, dims=n)
minimum(A...; B...) = Base.minimum(A...; B...)

any(X::AbstractArray, n::Int) = Base.any(X, dims=n)
any(A...; B...) = Base.any(A...; B...)

all(X::AbstractArray, n::Int) = Base.all(X, dims=n)
all(A...; B...) = Base.all(A...; B...)

(+)(A::AbstractArray, b::Real) = A .+ b
(+)(a::Real, B::AbstractArray) = a .+ B
(+)(A...; B...) = Base.:+(A...; B...)
(-)(A::AbstractArray, b::Real) = A .- b
(-)(a::Real, B::AbstractArray) = a .- B
(-)(A...; B...) = Base.:-(A...; B...)

findfirst(A::AbstractString, v::AbstractChar) = Base.findfirst(isequal(v), A)
findfirst(A...; B...) = Base.findfirst(A...; B...)

ind2sub(dims, ind) = Tuple(CartesianIndices(dims)[ind])

import Base: float, iterate, length
#Support float("1.23")
float(X::AbstractString) = Base.parse(Float64, X)
float(X::AbstractArray{T}) where {T<:AbstractString} = float.(X)

#Let Regex callable (so filter(r"...", vec) will work).
(R::Regex)(X::AbstractString) = occursin(R, X)

#Let Regex work as a scalar so it can be Broadcast (It's default in Julia 0.6 and also will be a default feature in Julia 1.1.0)
length(::Regex)=1
iterate(R::Regex)=(R, nothing)
iterate(R::Regex, ::Nothing)=nothing
export float, iterate, length

#Support ("String")' = "String"
import Base: adjoint
adjoint(S::Union{AbstractString, Symbol, AbstractChar})=S
export adjoint

const ASCIIString = String
const UTF8String = String
const ByteString = String

function getidx(obj, fds...)
    for fd in fds
        obj=Base.invokelatest(getindex, obj, fd)
    end
    obj
end

#[[ my.@il ]]
# my.@il(any_scripts) will automatically add Base.invokelatest to every function calling.
# 7 Aug 2019
_ilrecur!(X)=X
function _ilrecur!(X::Expr)
    if X.head==:function #In the case @il used before a function defination
        foreach(_ilrecur!, X.args[1].args[2:end])
        foreach(_ilrecur!, X.args[2:end])
    elseif X.head==:(=) && isa(X.args[1], Expr)
        if X.args[1].head==:ref # like a[b, c...]=d
            X.head=:call
            X.args=[:(Base.invokelatest), :setindex!, X.args[1].args[1], X.args[2], X.args[1].args[2:end]...]
            foreach(_ilrecur!, X.args[3:end])
        elseif X.args[1].head==:. && (!isa(X.args[2], Expr) || X.args[1].args[2].head!=:tuple) # like a.b=c
            X.head=:call
            X.args=[:(Base.invokelatest), :setproperty!, X.args[1].args[1], X.args[1].args[2], X.args[2]]
            foreach(_ilrecur!, X.args[3:end])
        else
            foreach(_ilrecur!, X.args)
        end
    elseif X.head==:call && X.args[1]!=:(Base.invokelatest) && X.args[1]!=:invokelatest # a(b)
        X.args=if length(X.args)>=2 && isa(X.args[2], Expr) && X.args[2].head==:parameters
            Any[:(Base.invokelatest), X.args[2], X.args[1], X.args[3:end]...]
        else
            Any[:(Base.invokelatest), X.args...]
        end
        foreach(_ilrecur!, X.args[2:end])
    elseif X.head==:. && (!isa(X.args[2], Expr) || X.args[2].head!=:tuple) # a.b
        X.head=:call
        X.args=Any[:(Base.invokelatest), :getproperty, X.args...]
        foreach(_ilrecur!, X.args[3:end])
    elseif X.head==:ref # a[b, c...]
        X.head=:call
        X.args=Any[:(Base.invokelatest), :getindex, X.args...]
        foreach(_ilrecur!, X.args[3:end])
    else
        foreach(_ilrecur!, X.args)
    end
    X
end
macro il(X)
    esc(_ilrecur!(X))
end

includelist=["myEnhanceSyntax.jl"
             "myRowAndGroup.jl"
             "myTable.jl"
             "myCalculation.jl"
             "myBioinformatics.jl"]
for fn in includelist
    include(fn)
end

function iterate(X::Union{eachrow, eachcol, eachcombr, eachgrp, fastgrp}, S=start(X))
    if done(X, S)
        nothing
    else
        next(X, S)
    end
end

end#module my

