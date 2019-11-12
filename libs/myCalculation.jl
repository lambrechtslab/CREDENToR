#{{ quantilels
#[[ quantilels ]]
#     quantilels(Data1[, Data2, ...]; q=0:0.1:1)
# Or: quantilels(Data; grp=Group|fastgrp, q=...)
#List the quantile and value on screen.
#So far this function output a Float or Any Matrix, but it is not guarantied in the future version.
#See also: ls, showprop, twofactormx
#Xiong Jieyi, 5 Sep 2014 > 25 May 2018 > 1 Nov 2019

function quantilels(Xs::AbstractVector...;q::AbstractVector{T}=Float64[], grp::Union{Group, fastgrp, Nothing}=nothing) where {T<:AbstractFloat}
    if !isnothing(grp)
        vX, vgrp=if isa(grp, fastgrp)
            (grpvec(grp, onlyone(Xs)), grp.grpid)
        else
            grpvec(grp, onlyone(Xs))
        end
        out=quantilels(vX...; q=q)
        return ["q" vgrp'; out]
    end
    if isempty(q)
        N=maximum(map(rownum,Xs))
        if N>199
            q=[0.0,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99,1.0]
        elseif N>49
            q=0.0:0.1:1
        else
            q=[0.0,0.25,0.5,0.75,1.0]
        end
    end
    O=map(Xs,1:length(Xs)) do X, i
        l=isnan.(X)
        if any(l)
            @warn(f"$1 ($2%) NaNs in data #$i are ignored."(sum(l),mean(l)*100))
            X=X[.!l]
        end
        X
    end
    [q map(x->quantile(vec(x),q),O)...]
end
quantilels(X::AbstractMatrix;args...)=quantilels(splitdim(X,1)...;args...)
quantilels(Xs::AbstractVector{Any}...;wargs...)=quantilels(map(x->retype(x,Real),Xs)...;wargs...) #For the convience of julia -e 'quantilels(SIC[:,1])'.
export quantilels
#}}

#{{ twofactormx
#[[ twofactormx ]]
# AnyMx = twofactormx(factorA, factorB; showsum=false)
#Or FreqMx, row, col = twofactormx(...; out1mx=false, ...)
#Calculated frequenes for each factor combinations, and output the result as a matrix. Input can be groups or fastgrp objects. factorA and factorB will be as row and columns respectively.
#If showsum=true, the last row and column in the outputed matrix is the sum values. Note that when out1mx=false, the outputed row and col will not change no matter showsum=true or false.
#See also: freq
#Xiong Jieyi, 8 Nov 2017

export twofactormx
function twofactormx(A::fastgrp, B::fastgrp; out1mx::Bool=true, showsum::Bool=false)
    M=fill(0, A.grpnum, B.grpnum)
    for ia=1:A.grpnum
        la=want(A,ia)
        for ib=1:B.grpnum
            lb=want(B,ib)
            M[ia, ib]=length(intersect(la,lb))
        end
    end
    if showsum
        M=[M sum(M,2);sum(M,1) sum(M)]
        if out1mx
            asmx(M, row=[A.grpid;"SUM"], col=[B.grpid;"SUM"])
        else
            (M, A.grpid, B.grpid)
        end
    else
        if out1mx
            asmx(M, row=A.grpid, col=B.grpid)
        else
            (M, A.grpid, B.grpid)
        end
    end
end
twofactormx(A, args...; wargs...)=twofactormx(fastgrp(A), args...; wargs...)
twofactormx(A::fastgrp, B; wargs...)=twofactormx(A, fastgrp(B); wargs...)
#}}

#{{ rankgrp

#[[ rankgrp ]]
# GroupID=rankgrp(V,N)
#Group by the rank of V based on drawer rule. Return a int vector of group NO.
#Note that considering the repeat elements in the vector, different groups could be have same values.
#See also: drawer, balancemask
#Xiong Jieyi, October 22, 2014

#<v0.6# function rankgrp{T<:Real}(V::AbstractVector{T},N::Real)
function rankgrp(V::AbstractVector{T},N::Real) where {T<:Real}
    gn=drawer(length(V),N)
    gl=vcat([fill(i,gn[i]) for i=1:N]...)
    return gl[rvorder(sortri(V))]
end
export rankgrp

#}}

#{{ kmean, kmeancost
#[[ kmean kmeancost ]]
# using Clustering
# ClusterNO=kmean(X, k; kmeans_param...)
# totalcost=kmeancost(X, k|k_vec; kmeans_param...)
#Calculate cluster of X' using kmeans, and return the assignments(group NO) or totalcost.
#Note that X are clusted by row, which is not like the Clustering.kmeans.
#See also: kmeancostplot, grpboxplot, fastcluster, fcluster, dendrogram.
#Xiong Jieyi, 3 Oct 2014 > 13 May 2019

# function needClustering()
#     isdefined(:Clustering) || require("Clustering")
#     nothing
# end
kmean(X::Matrix,k::Real;wargs...)=invokelatest(importpkg(:Clustering).kmeans, X', k; wargs...).assignments
kmeancost(X::Matrix,k::Real;wargs...)=invokelatest(importpkg(:Clustering).kmeans, X', k; wargs...).totalcost
function kmeancost(X::Matrix,ks::AbstractVector{T};wargs...) where {T<:Real} 
    @pkgfun(Clustering, kmeans)
    map(k->kmeans(X', k; wargs...).totalcost, ks)
end
export kmean, kmeancost
#}}

#{{ fastcluster,fcluster,dendrogram,leaves
#[[ fastcluster linkage fcluster dendrogram leaves ]]
# Z=fastcluster(Matrix; method="ward", metric="euclidean") #Make tree by data
# Z=linkage(DistancesMatrix|DistancesVector; method="ward") #Make tree by distance matrix.
# cluster_NO=fcluster(Z, cutoff; criterion="distance(D)|maxclust", ...)   #Cut tree
# dendrogram(Z; labels=None, color_threshold=height_value, ...)     #Draw tree(Only used in Qt so far.)
# leaves_index = leaves(Z)  #The return corresponds to the observation vector index as it appears in the tree from left to right. e.g. leaves_list(Z).+1 in Python.
#Use scipy.cluster.hierarchy and fastcluster to create a hierarchy tree.
# method could be single, complete, average, mcquitty, centroid, median, and ward.
# metric could be euclidean, maximum, manhattan, canberra, binary, minkowski, correlation, and spearman. The correlation distance is actually calculated by scipy.spatial.distance, and the spearman distance calculated by 1-corspearman(Matrix'). These two metrics are not originally supported by fastcluster.
# criterion in fcluster could be inconsistent, distance(default), maxclust, moncrit and maxclust_monocrit. For the detail, see: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.cluster.hierarchy.fcluster.html#scipy.cluster.hierarchy.fcluster
# fastcluster is a stable python package (https://pypi.org/project/fastcluster/). It needs to be preinstalled in PYTHONPATH, like  PYTHONPATH=/home/jxiong/programs/fastcluster/fastcluster-1.1.23/lib/python2.7/site-packages:$PATYONPATH .
#See also: kmean, kmeancost, distmx2vec
#Xiong Jieyi, 3 Oct 2014>October 31, 2014>Aug 14, 2015>11 Jan 2019

function fastcluster(X::AbstractMatrix;
                     method::Union{String,Symbol}="ward",
                     metric::Union{String,Symbol}="euclidean", wargs...)
    @pypkgfun(fastcluster, linkage)
    if metric=="correlation"
        @pkgfun(Distances,pairwise,CorrDist)
        D=distmx2vec(pairwise(CorrDist(),X'))
        linkage(D; method=method, wargs...)
    elseif metric=="spearman"
        D=max(distmx2vec(1-corspearman(X')),0)
        linkage(D; method=method, wargs...)
    else
        linkage(X; method=method, metric=metric, wargs...)
    end
end
linkage(X::AbstractMatrix; method::Union{String,Symbol}="ward", wargs...)=invokelatest(importpy("fastcluster").linkage, distmx2vec(X); method=method, wargs...)
linkage(X::AbstractVector; method::Union{String,Symbol}="ward", wargs...)=invokelatest(importpy("fastcluster").linkage, X;method=method,wargs...)
fcluster(Z::Matrix,t::Real; criterion::Union{String,Symbol}="distance",wargs...)=invokelatest(importpy("scipy.cluster.hierarchy").fcluster, Z,t;criterion=criterion,wargs...)
dendrogram(args...;wargs...)=invokelatest(importpy("scipy.cluster.hierarchy").dendrogram,args...;wargs...)
leaves(Z::Matrix)=invokelatest(importpy("scipy.cluster.hierarchy").leaves_list,Z).+1
export fastcluster, linkage, fcluster, dendrogram, leaves
#}}

#{{ distmx2vec
#[[ distmx2vec distvec2mx ]]
# distance_vector = dismx2vec(distance_matrix)
# distance_matrix = disvec2mx(distance_vector)
#Convert between distance matrix and consense distance vector.
#The distances are arranged in the order (2,1), (3,1), ..., (m,1), (3,2), ..., (m,2), ..., (m,m–1)).
#See also: fastcluster, distvec
#Xiong Jieyi, 3 Oct 2014 > 3 Jun 2019

export distmx2vec, distvec2mx
function distmx2vec(M::AbstractMatrix{T}) where {T<:Real}
    n=size(M, 1)
    size(M, 2)==n || error("Inputed matrix is not a squire.")
    V=T[]
    for i=1:n
        for j=i+1:n
            M[i, j]==M[j, i] || error("Inputed matrix is not symmatric.")
            push!(V, M[j, i])
        end
    end
    V
end
function distvec2mx(V::AbstractVector{T}) where {T<:Real}
    t=sqrt(8*length(V)+1)/2+0.5
    if !isinteger(t)
        error("Invalid distance vector length.")
    end
    n=Int(t)
    
    M=fill(zero(T), n, n)
    p=1
    for i=1:n
        for j=i+1:n
            M[j, i]=V[p]
            M[i, j]=V[p]
            p+=1
        end
    end
    M
end
#}}

#{{ distancevec
#[[ distancevec ]]
# vector = distvec(fun(row_i, col_i)::Float64, n::Int)
# Calculate consense distance vector. Note that the input of given function is a pair of scaler indice.
# The distances are arranged in the order (2,1), (3,1), ..., (n,1), (3,2), ..., (n,2), ..., (n,n–1)).
# See also: fastcluster, distmx2vec, distvec2mx
# Xiong Jieyi, 3 Jun 2019

export distancevec
function distancevec(fun::Function, N::Int)
    V=Float64[]
    for i=1:N
        for j=i+1:N
            push!(V, fun(j,i))
        end
    end
    V
end
#}}

#{{ StatsBase
#[[ corspearman ]]
# pho=corspearman(A[, B])
#Calculate Spearman correlation. Input could be vector or matrix. If input is matrix, function calculated column-wisely.
#See also: cor
#StatsBase Pkg, 3 Oct 2014

@pkgfun(StatsBase, corspearman)
# function corspearman(args...;wargs...)
#     isdefined(:StatsBase)||require("StatsBase")
#     Main.StatsBase.corspearman(args...;wargs...)
# end
export corspearman

#[[ nrank ]]
#rand=nrank(x[, dim]; rev=false)
#Compute tied ranking (also known as fractional ranking or 1 2.5 2.5 4 ranking) for x, using StatsBase Pkg. rev=true will output the reversed rank (the maximum is 1 and the minimum is N).
#The difference between nrank and tiedrank in StatsBase.pkg: nrank will ignore NaNs in calculation (with a warning) but refill NaNs in the rank.
#See also: rankgrp
#Xiong Jieyi, November 2, 2014 >Mar 3, 2016>May 9, 2016
function nrank(X::AbstractVector; rev::Bool=false)
    # sb=importpkg(:StatsBase)
    @pkgfun(StatsBase, tiedrank)
    l=isnan(X)
    any(l) && @warn(f"$1 of $2 data are NaNs. They are also NaNs in the ranks."(sum(l),length(l)))
    oR=tiedrank(X[~l])
    R=fill(NaN,size(X))
    R[~l]=oR
    rev ? length(oR)-R+1 : R
end
function nrank(X::AbstractMatrix,dim::Real;wargs...)
    if dim==1
        rowfun(x->nrank(x;wargs...),X',as_vec=true)'
    elseif dim==2
        rowfun(x->nrank(x;wargs...),X,as_vec=true)
    else
        error("dim only support 1 and 2 so far.")
    end
end
export nrank
#}}

#{{ HypothesisTest [ obsoleted ]
#pval=pvalue(fisherexacttest(a,b,c,d))
#pval=pvalue(fisherexacttest([a,b,c,d]))
#pval=pvalue(fisherexacttest(BitVector1,BitVector2))
#Calculate two-tailed Fisher exact test.
#See also: pvalue, fishertest
#Xiong Jieyi, 3 Oct 2014

# @pkgfun(HypothesisTests,FisherExactTest=>fisherexacttest,pvalue)
# function fisherexacttest(M::Matrix)
#     @assert(size(M)==(2,2),"FisherExactTest should be 2x2 matrix.")
#     fisherexacttest(M[1,1],M[1,2],M[2,1],M[2,2])
# end
# fisherexacttest(A::BitVector,B::BitVector)=fisherexacttest(sum(A),sum(B),sum(.!A),sum(.!B))
# export fisherexacttest, pvalue
#}}

#{{ glm12test

#[[ glm12test ]]
# pval=glm12test(Y,X;family="",test="F|Chisq",outall=false)
# Test if data fit linear or quadrotic model, using lm (when family is omitted) or glm. Return the minimal F-test p-value of these two models, or the both two p-values and the p-value between two models if outall is ture.
# family could be:
# binomial(link = "logit")
# gaussian(link = "identity")
# Gamma(link = "inverse")
# inverse.gaussian(link = "1/mu^2")
# poisson(link = "log")
# quasi(link = "identity", variance = "constant")
# quasibinomial(link = "logit")
# quasipoisson(link = "log")
# nb --Need run eR("library(MASS)") in advance.
#See also: callr, rdata, relem, rshow, callr("p.adjust",pval,method="BH")
#Xiong Jieyi, October 15, 2014>October 23, 2014

function glm12test(Y::Array,X::Vector;family::AbstractString="",test::AbstractString="F",outall::Bool=false)
    @assert(size(Y,1)==length(X),"Y and X have different length.")
    # Rif.isinitialized()||Rif.initr()
    # import needR
    @pkgfun(myR, callr, eR)
    if isa(Y,Matrix) && size(Y,2)==2 && size(Y,1)>1
        df=callr("data.frame",Y1=Y[:,1],Y2=Y[:,2],X=X)
        sY="cbind(Y1,Y2)"
    else
        df=callr("data.frame",Y=Y,X=X)
        sY="Y"
    end
    if isempty(family)
        md0=callr("lm",eR("$(sY)~1"),data=df)
        # md1=callr("lm",R("$(sY)~X"),data=df)
        # md2=callr("lm",R("$(sY)~I(X^2)+X+1"),data=df)
    elseif family=="nb"
        md0=callr("glm.nb",eR("$(sY)~1"),data=df)
        # md1=callr("glm.nb",R("$(sY)~X"),data=df)
        # md2=callr("glm.nb",R("$(sY)~I(X^2)+X+1"),data=df)
    else
        md0=callr("glm",eR("$(sY)~1"),data=df,family=eR(family))
        # md1=callr("glm",R("$(sY)~X"),data=df,family=R(family))
        # md2=callr("glm",R("$(sY)~I(X^2)+X+1"),data=df,family=R(family))
    end
    md1=callr("update",md0,eR(".~.+X"))
    md2=callr("update",md1,eR(".~.+I(X^2)"))

    pval1=callr("anova",md1,md0,test=test)|>x->relem(x,length(x))|>x->[x...]|>x->x[~isnan(x)]|>x->(@assert(length(x)==1);x[1])
    pval2=callr("anova",md2,md0,test=test)|>x->relem(x,length(x))|>x->[x...]|>x->x[~isnan(x)]|>x->(@assert(length(x)==1);x[1])
    if outall
        pval12=callr("anova",md2,md1,test=test)|>x->relem(x,length(x))|>x->[x...]|>x->x[~isnan(x)]|>x->(@assert(length(x)==1);x[1])
        return pval1,pval2,pval12
    else
        return min(pval1,pval2)
    end
end
export glm12test

#}}

#{{ nthmax nthmin nthmaxh nthminh
#[[ nthmin nthmax nthmaxh nthminh ]]
# N-th-max|min = nthmax|nthmin(N, Vector)
# N-th-max|min = nthmax|nthmin(N, Matrix, dim)
# Vector = nthmaxh|nthminh(N, Matrix) # equal to vec(nthmax|nthmin(N,M,2))
#Calculate the Nth minimum or maximum values. Running time is linearly increased with N.
#See also:
#Xiong Jieyi, March 7, 2015>Jun 6, 2016

export nthmin, nthmax, nthminh, nthmaxh
function nthmin(n::Integer,X::AbstractVector)
    if n>1
        X=deepcopy(X)
        for i=2:n
            deleteat!(X,argmin(X))
        end
    end
    minimum(X)
end
#<v0.6# function nthmin{T}(n::Integer,X::AbstractMatrix{T},dim::Int)
function nthmin(n::Integer,X::AbstractMatrix{T},dim::Int) where {T}
    if dim==1
        T[nthmin(n,X[:,i]) for i=1:size(X,2)]'
    else
        @assert(dim==2,"dim can only be 1 or 2.")
        hcat(T[nthmin(n,vec(X[i,:])) for i=1:size(X,1)])
    end
end
function nthmax(n::Integer,X::AbstractVector)
    if n>1
        X=deepcopy(X)
        for i=2:n
            deleteat!(X,argmax(X))
        end
    end
    maximum(X)
end
#<v0.6# function nthmax{T}(n::Integer,X::AbstractMatrix{T},dim::Int)
function nthmax(n::Integer,X::AbstractMatrix{T},dim::Int) where {T}
    if dim==1
        T[nthmax(n,X[:,i]) for i=1:size(X,2)]'
    else
        @assert(dim==2,"dim can only be 1 or 2.")
        hcat(T[nthmax(n,vec(X[i,:])) for i=1:size(X,1)])
    end
end
nthminh(n::Integer,X::AbstractMatrix)=vec(nthmin(n,X,2))
nthmaxh(n::Integer,X::AbstractMatrix)=vec(nthmax(n,X,2))
#}}

#{{ slidwin
#[[ slidwin ]]
# (fun_out, out_x) = slidwin(fun, Y::Group; x::Group=1:N, winlen= ~N/10(oddnum), x_sorted=false)
#Slid-window calculation and output fun(Y_in_window) and the mid-x of window. For input Y and x in row number N, the output will have N-winlen+1 rows. Function output will be pile up by mapr.
#If x is sorted, set x_sorted=true to save time. Note that output will be wrong without any warning if x is not in sort.
#See also: mapr
#Xiong Jieyi, Feb 29, 2016

export slidwin
function slidwin(fun::Function, Y::Group; x::Group=1:length(Y), winlen::Integer=floor(Int,length(Y)/10)|>x->mod(x,2)==1 ? x : x+1, x_sorted::Bool=x==1:length(Y))
    @assert(mod(winlen,2)==1 && winlen>0, "winlen must be a positive odd integer.")
    if !x_sorted
        x,t=sortr(x)
        Y=getrow(Y,t)
    end
    wing=int((winlen-1)/2)
    mapr(wing+1:length(Y)-wing) do i
        (fun(Y[i-wing:i+wing]),x[i])
    end
end
#=
#<v0.6# function slidwin{T<:Real}(fun::Function,V::Vector{T};lim::(Real,Real)=(minimum(V),maximum(V)),winlen::Real=0.1*(lim[2]-lim[1]),stepn::Integer=0,steplen::Real=NaN)
function slidwin(fun::Function,V::Vector{T};lim::(Real,Real)=(minimum(V),maximum(V)),winlen::Real=0.1*(lim[2]-lim[1]),stepn::Integer=0,steplen::Real=NaN) where {T<:Real}
    if stepn==0
        if isnan(steplen)
            stepn=100
        else
            stepn=int((lim[2]-winlen-lim[1])/steplen)+1
        end
    end
    if isnan(steplen)
        steplen=(lim[2]-winlen-lim[1])/(stepn-1)
    end
    rp=rowpile()
    for i=1:stepn
        lim[1]+steplen*(i-1)
    end
    end
    =#
#}}

#{{ sumh,medianh,meanh,varh,stdh,maximumh,minimumh,anyh,allh
#[[ sumh medianh meanh varh stdh maximumh maxh minimumh minh anyh allh ]]
# sumh(X::AbstractMatrix)=vec(sum(X,2)) #The same for other functions.
# maxh and minh is equal to maximumh and minimumh.
#Just a syntax suger.
#see also: truesbut, falsesbut
#Xiong Jieyi, Jun 24, 2015 >Jan 19, 2016

export sumh,medianh,meanh,varh,stdh,maximumh,maxh,minimumh,minh,anyh,allh
sumh(X::AbstractMatrix)=vec(sum(X,2))
medianh(X::AbstractMatrix)=vec(median(X,2))
meanh(X::AbstractMatrix)=vec(mean(X,2))
varh(X::AbstractMatrix)=vec(var(X,2))
stdh(X::AbstractMatrix)=vec(std(X,2))
maximumh(X::AbstractMatrix)=vec(maximum(X,2))
maxh(X::AbstractMatrix)=vec(maximum(X,2))
minimumh(X::AbstractMatrix)=vec(minimum(X,2))
minh(X::AbstractMatrix)=vec(minimum(X,2))
anyh(X::AbstractMatrix)=vec(any(X,2))
allh(X::AbstractMatrix)=vec(all(X,2))
#}}

#{{ mx2DataFrame

#[[ mx2DataFrame ]]
# DF =  mx2DataFrame(Matrix; col=vector, row=vector)
#Create a DataFrame from matrix. This funcion is useful to output data tidily on screen.
#See also: pandasDF, tb2DataFrame
#Xiong Jieyi, Feb 13, 2016

export mx2DataFrame
function mx2DataFrame(M::Matrix; col::Vector=map(f"C$1",1:size(M,2)), row::Vector=[])
    df=importpkg(:DataFrames)
    C=splitdimv(M,1)
    if !isempty(row)
        C=Vector[row, C...]
        col=[""; col]
    end
    df.DataFrame(;map((x,y)->Pair(symbol(x),y),col,C)...)
end

#}}

#{{ quantilenorm
#[[ quantilenorm ]]
# M=quantilenorm(Matrix)
# V1,V2,...=quantilenorm(V1,V2,...)
#Quantile normalization for each column of matrix or each inputed vector.
#See also:
#Xiong Jieyi, Aug 25, 2016

export quantilenorm
#<v0.6# function quantilenorm{T<:Real}(M::AbstractMatrix{T})
function quantilenorm(M::AbstractMatrix{T}) where {T<:Real}
    if size(M,2)==1
        return M
    end
    ln=size(M,1)
    IM=hcat(map(sortperm,eachcol(M))...)
    SM=hcat(map((x,l)->x[l],eachcol(M),eachcol(IM))...)
    Sme=rowfun(mean,SM,as_vec=true)
    reshape(Sme[vec(IM)],size(IM))
end
quantilenorm(Vs::AbstractVector...)=splitdim(quantilenorm(hcat(Vs...)),1)
#}}

#{{ showprop
#[[ showprop ]]
# showprop(Bool_Array)
#Print sum, length and mean of the Bool vector.
#See also: quantilels, saynum
#Nov 10, 2016 > 9 Nov 2018

export showprop
function showprop(X::AbstractArray{Bool})
    N=length(X)
    TN=sum(X)
    FN=N-TN
    n, tn, fn=num2str([N, TN, FN], fixdigit=true, leftfill=' ')
    hn=saynum(N)
    htn=saynum(TN)
    hfn=saynum(FN)
    TP=TN/N
    FP=FN/N
    println("T: $tn\t$htn\t$TP")
    println("F: $fn\t$hfn\t$FP")
    println("S: $n\t$hn")
end
#}}

#{{ princomp
#[[ princomp ]]
#(score_matrix, contribution_precentiages) = princomp(M; method=:pca|:pca_R|:tSNE_py|:mds, tsne_init_pc=0 or 50)
#Calcuate PCA score matrix and contribution precentages using MultivariateStats package (method=:pca) or R"prcomp" (method=:pca_R).
#This function also support tSNE (method=:tSNE_py) and MDS (method=:mds). In these modes, the 2nd output is always nothing. For tSNE, tsne_init_pc assign the initial PC number. By default, tsne_init_pc=50 if size(M, 2)>100, or no pre-PCA otherwise.
#See also: pcaplot
#Xiong Jieyi, Nov 21, 2016 > 10 Oct 2019

export princomp
function princomp(M::AbstractMatrix{TT}; method::Symbol=:pca, tsne_init_pc::Int=size(M, 2)>100 ? 50 : 0) where {TT<:Real}
    if method==:pca
        @pkgfun(MultivariateStats, fit=>__mvs_fit, transform=>__mvs_transform, principalvars=>__mvs_principalvars, principalratio=>__mvs_principalratio)
        pca_obj=__mvs_fit(importpkg(:MultivariateStats).PCA, M')
        pcaM=trsp(__mvs_transform(pca_obj,M'))
        pc_contri=__mvs_principalvars(pca_obj)|>x->__mvs_principalratio(pca_obj)*x./sum(x)
    elseif method==:pca_R
        @pkgfun(myR,callr,relemj,relem,r2j)
        pca_Robj=callr(:prcomp,M,center=true,scale=false)
        pcaM=relemj(pca_Robj,"x")
        t=relem(callr(:summary,pca_Robj),"importance")
        rn=r2j(callr(:row!names,t), keepvec=true)
        pc_contri=vec(r2j(t, keepvec=true)[findonly(rn.=="Proportion of Variance"),:])
    elseif method==:tSNE_py
        #According to http://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html#examples-using-sklearn-manifold-tsne
        if tsne_init_pc>0
            M=princomp(M, method=:pca)[1]
            if size(M, 2)>tsne_init_pc
                M=M[:, 1:tsne_init_pc]
            end
        end
        pcaM = @il importpy("sklearn.manifold").TSNE().fit_transform(M)
        pc_contri=nothing
    elseif method==:mds
        @pkgfun(MultivariateStats, classical_mds)
        pcaM=classical_mds(M, 2)'
        pc_contri=nothing
    else
        error("Unsupported method $method.")
    end
    (pcaM, pc_contri)
end
#}}

#{{ fdrcutoff
#[[ fdrcutoff ]]
# pval_cutoff = fdrcutoff(pval_test, pval_H0; fdr=0.05, max_fdr=min(3*fdr, 1))
#Calculate the p value cutoff, based on test sets and H0-control sets. The calculation will stop when fdr screening reached the max_fdr. If chosen a cutoff is unfeasible, function return zero or the maximum inputed value instead, along with a warning.
#See also: glm12test
#Xiong Jieyi, Jun 27, 2017

export fdrcutoff
#<v0.6# function fdrcutoff{T1<:Real, T2<:Real}(P1::Vector{T1}, P0::Vector{T2}; fdr::AbstractFloat=0.05, max_fdr::AbstractFloat=min(3*fdr,1.0))
function fdrcutoff(P1::Vector{T1}, P0::Vector{T2}; fdr::AbstractFloat=0.05, max_fdr::AbstractFloat=min(3*fdr,1.0)) where {T1<:Real, T2<:Real}
    (Ps, si)=sortr([P1;P0])
    Ls=[trues(length(P1));falses(length(P0))][si]
    N1=0
    N0=0
    Pcutoff=zero(Ps[1])
    for i=1:length(Ps)
        if Ls[i]
            N1+=1
        else
            N0+=1
        end
        if i<length(Ps) && Ps[i+1]==Ps[i]
            continue
        end
        cP=N0/N1
        if cP<=fdr
            Pcutoff=Ps[i]
        elseif N1>=100 && cP>max_fdr
            break
        end
    end
    if Pcutoff==0
        @warn("No any cutoff could satisfy the FDR. fdrcutoff output zero instead.")
    elseif Pcutoff==Ps[end]
        @warn("Any cutoff could satisfy the FDR. fdrcutoff output the maximum inputed value instead.")
    end
    Pcutoff
end
#}}

#{{ zscore
#[[ zscore ]]
# score_mx = zscore( Matrix; dims=0, baseon=X )
# Normalize data to mean=0, std=1. If baseon is assigned, the mean and std will calculated on the assigned data.
# See also:
#Xiong Jieyi, 27 Oct 2017 > 25 Oct 2018

export zscore
function zscore(X::Union{AbstractArray, Base.SkipMissing}, D::Integer=0; baseon::Union{AbstractArray, Base.SkipMissing}=X, dims::Integer=D)
    if dims==0
        (X.-mean(baseon))./std(baseon)
    else
        (X.-mean(baseon, dims=dims))./std(baseon, dims=dims)
    end
end
#}}

#{{ modebyr
#[[ modebyr ]]
# mode_element = modebyr(Group; big=false, sorted=false, sortcheck=true)
# Find the mode, i.e. the most occurred element (by row) in a group-data. If the input is sorted, set sorted=true to speed up.
# If more than one category have the same maximum counts, the one with the smallest group ID (big=false), or the one with the biggest group ID (big=true), will be returned.
# See also: freq, idennum, uninum
# Xiong Jieyi, 5 Nov 2019

export modebyr
function modebyr(X::Group; sorted::Bool=false, sortcheck::Bool=true, big::Bool=false)
    idx=if sorted
        if sortcheck && !issortedr(X)
            error("Input is not sorted.")
        end
        1:rownum(X)
    else
        sortri(X)
    end
    P=1
    maxct=0
    maxI=1
    for i=2:length(idx)
        if !isequalr2(X, idx[i], X, idx[i-1])
            cct=i-P
            P=i
            if big ? (maxct<=cct) : (maxct<cct)
                maxct=cct
                maxI=i-1
            end
        end
    end
    cct=length(idx)-P+1
    if big ? (maxct<=cct) : (maxct<cct)
        maxI=length(idx)
    end
    getrow(X, idx[maxI])
end
#}}

#{{ iselequal
#[[ iselequal ]]
# T|F = iselequal([function, ]Iternation)
#Calculate if all the elements, or all the outputs of fun(element), are equal (isequal() return true). Empty or one-element iternation is always return true.
#See also:
#Xiong Jieyi, 9 Jul 2018

export iselequal
if VERSION<v"0.7.0"
    function iselequal(X)
        st=start(X)
        done(X, st) && return true
        x1, st=next(X, st)
        while !done(X, st)
            cx, st=next(X, st)
            isequal(cx, x1) || return false
        end
        return true
    end
else
    function iselequal(X)
        out=iterate(X)
        if out==nothing
            return true
        else
            x1, st=out
        end
        while true
            out=iterate(X, st)
            out==nothing && return true
            cx, st=out
            isequal(cx, x1) || return false
        end
    end
end
iselequal(fun::Function, X)=iselequal(fun(cx) for cx in X)
#}}

#{{ valid
#[[ valid ]]
# Y1, Y2, ... = valid(X1, X2, ...; delif=x->ismissing(x) || isnan(x) || isinf(x), follower=Set|indics|TF_vec)
# Remove the invalid values (missing, NaN and Inf by default, or any others where delif return true) in float vectors. Inputs should be vectors in the same lengths. The element "invalid" in any of input vectors will be removed. The outputs are always in the same lengths.
# The #i input where i in `follower` will not be checked by `delif`, but only be filtered at last.
# It is useful such as cor(valid(X, Y)...). Avoid to use it in the preformance-critical codes.
#See also:
#Xiong Jieyi, 17 Aug 2018 > 30 Jan 2019>21 Feb 2019

export valid
function valid(X1::AbstractVector, Xs::AbstractVector...; delif::Function=x->ismissing(x) || isnan(x) || isinf(x), follower=Set())
    if isa(follower, AbstractVector{Bool})
        follower=Set(findall(follower))
    elseif !isa(follower, AbstractSet)
        follower=Set(follower)
    end
    l=in(1, follower) ? falses(length(X1)) : delif.(X1)
    if isempty(Xs)
        return X1[.!l]
    else
        for i in length(Xs)
            if !in(i+1, follower)
                l=l .| delif.(Xs[i])
            end
        end
        return X1[.!l], map(x->x[.!l], Xs)...
    end
end
#}}

#{{ grp2mx
#[[ grp2mx ]]
# RowId, ColId, Mx1, Mx2, ... = grp2mx(RowGrp, ColGrp, V1[=>default], V2[=>default], ...;
#                               (option I)    stable|row_stable|col_stable=false,
#                               (option II)   row_order=..., col_order=...) #Should be row-unique.
# Convert a table like vector to a matrix.
# If the '=>default' is not given, the default value will be the cooresponding emptyval. Note that input as fastgrp object is not supported. default value should be in the same type as V, or be missing.
# Note that this function will not check 1) whether (RowGrp, ColGrp) are row-unique; 2) whether row_order/col_order are row-unique. For 2, the value will NOT be duplicated but just assigned with default values in all-but-one row/columns. Utilize this feature is not recommended (may change in the future).
# See also: grpvec
# Xiong Jieyi, 3 Jun 2019

export grp2mx
function grp2mx(R::Group, C::Group, Vs...;
                stable::Bool=false, row_stable::Bool=stable, col_stable::Bool=stable,
                row_order::Union{Group, Nothing}=nothing,
                col_order::Union{Group, Nothing}=nothing)
    if row_order==nothing
        ri, rl=uniqr(R; stable=row_stable)
        uR=getrow(R, ri)
        rlen=length(ri)
    else
        uR=row_order
        rl=memberr(R, uR)
        rlen=rownum(uR)
    end
    if col_order==nothing
        ci, cl=uniqr(C; stable=col_stable)
        uC=getrow(C, ci)
        clen=length(ci)
    else
        uC=col_order
        cl=memberr(C, uC)
        clen=rownum(uC)
    end
    Ms=map(Vs) do V
        if isa(V, Pair)
            df=V[2]
            V=V[1]
        else
            df=emptyval(eltype(V))
        end
        M=fill(df, rlen, clen)
        if ismissing(df)
            M=convert(Matrix{Union{eltype(V), Missing}}, M)
        end
        foreach(rl, cl, V) do r, c, v
            if r>0 && c>0
                M[r, c]=v
            end
        end
        M
    end
    (uR, uC, Ms...)
end
#}}

#{{ pysparse2coo pysparse2jl
#[[ pysparse2coo pysparse2jl ]]
# rowIdx, colIdx, value = pysparse2coo( python_sparseMatrix_obj ) # indice are 1-based.
# julia_sparse_mx = pysparse2jl( python_sparseMatrix_obj ) # `using SparseArrays` required.
# it_self = pysparse2jl( julia_mx ) # Just for convenience.
# Convert python sparse matrix to three vectors. The output indice are 1 based. PyCall should be preloaded.
# See also: grp2mx
# Jieyi Xiong 17 Jun 2019

export pysparse2coo, pysparse2jl
function pysparse2coo(X)
    # ri, ci, v=importpy("scipy.sparse").find(X::Main.PyCall.PyObject)
    # (ri.+1, ci.+1, v)
    X=(X::Main.PyCall.PyObject).tocoo()
    (X.row.+1, X.col.+1, X.data)
end
pysparse2jl(X)=invokelatest(importpkg(:SparseArrays).sparse, pysparse2coo(X)..., X.shape...)
pysparse2jl(X::AbstractArray)=X
#}}

#{{ Some common R functions
export wilcoxtest, studenttest, padjust, cortest, fishertest, binomtest
#[[ wilcoxtest ]]
# p.value = wilcoxtest(x, y = NULL,
#     alternative = c("two.sided", "less", "greater"),
#     mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
#     conf.int = FALSE, conf.level = 0.95, ...; vb=0)
# Performs one- and two-sample Wilcoxon tests on vectors of data; the latter is also known as 'Mann-Whitney' test.
# x: numeric vector of data values.  Non-finite (e.g., infinite or missing) values will be omitted.
# y: an optional numeric vector of data values: as with 'x' non-finite values will be omitted.
# alternative: a character string specifying the alternative hypothesis, must be one of '"two.sided"' (default), '"greater"' or '"less"'.  You can specify just the initial letter.
# paired: a logical indicating whether you want a paired test.
#For detail, see rhelp("wilcox.test")
#See also: studenttest, padjust, cortest, fishertest, binomtest
#Xiong Jieyi, 27 Oct 2017

function wilcoxtest(arg...; vb::Int=0, warg...)
    @pkgfun(myR, callrw, relem, r2j, rshow)
    o=callrw("wilcox.test", arg...; warg...)
    rshow(o; vb=vb)
    r2j(relem(o, "p.value"), NA=NaN, keepvec=false)
end

#[[ studenttest ]]
# p.value, meanX - meanY, [CI_low, CI_high] = studenttest(x, y = NULL,
#        alternative = c("two.sided", "less", "greater"),
#        mu = 0, paired = FALSE, var.equal = FALSE,
#        conf.level = 0.95, ...; vb=0)
# Performs one and two sample t-tests on vectors of data.
# x: a (non-empty) numeric vector of data values.
# y: an optional (non-empty) numeric vector of data values.
# paired: a logical indicating whether you want a paired t-test.
# For detail, see rhelp("t.test")
# See also: studenttest, padjust, cortest, fishertest, binomtest
# Xiong Jieyi, 27 Oct 2017

function studenttest(arg...; vb::Int=0, warg...)
    @pkgfun(myR, callrw, relem, r2j, rshow)
    o=callrw("t.test", arg...; warg...)
    rshow(o; vb=vb)
    pval=r2j(relem(o, "p.value"), NA=NaN, keepvec=false)
    mn=r2j(relem(o, "estimate"), NA=NaN, keepvec=true)
    dmn=if length(mn)==1 # if paired=true
        mn[1]
    else
        @assert length(mn)==2
        mn[1]-mn[2]
    end
    cfi=r2j(relem(o, "conf.int"), NA=NaN, keepvec=true)
    (pval, dmn, cfi)
end
#[[ padjust ]]
#  padj = padjust(p, method = p.adjust.methods)
#  Given a set of p-values, returns p-values adjusted using one of several methods.
#  p: numeric vector of p-values (possibly with 'NA's).  Any other R is coerced by 'as.numeric'.
#  p.adjust.methods: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
# For detail, see rhelp("p.adjust")
# See also: wilcoxtest, studenttest, cortest, fishertest, binomtest
# Xiong Jieyi, 27 Oct 2017

function padjust(arg...; warg...)
    @pkgfun(myR, callr, r2j)
    r2j(callr("p.adjust", arg...; warg...), NA=NaN, keepvec=false)
end

#[[ cortest ]]
# r, p.value = cortest(x, y,
#          alternative = c("two.sided", "less", "greater"),
#          method = c("pearson", "kendall", "spearman"),
#          exact = NULL, conf.level = 0.95, continuity = FALSE, ...; vb=0)
# Test for association between paired samples, using one of Pearson's product moment correlation coefficient, Kendall's tau or Spearman's rho.
# x, y: numeric vectors of data values.  'x' and 'y' must have the same length.
# method: a character string indicating which correlation coefficient is to be used for the test.  One of '"pearson"', '"kendall"', or '"spearman"', can be abbreviated.
# For detail, see rhelp("cor.test")
# See also: wilcoxtest, studenttest, padjust, fishertest, binomtest
# Xiong Jieyi, 27 Oct 2017

function cortest(arg...; vb::Int=0, warg...)
    @pkgfun(myR, callrw, relem, r2j, rshow)
    o=callrw("cor.test", arg...; warg...)
    rshow(o; vb=vb)
    r2j.(relem(o, "estimate", "p.value"), NA=NaN, keepvec=false)
end

#[[ fishertest ]]
# p.value = fishertest(x, y = NULL, workspace = 200000, hybrid = FALSE,
#   control = list(), or = 1, alternative = "two.sided(default)|greater|less"",
#   conf.int = TRUE, conf.level = 0.95,
#   simulate.p.value = FALSE, B = 2000; vb=0)
# Performs Fisher's exact test for testing the null of independence of rows and columns in a contingency table with fixed marginals.
# x: either a two-dimensional contingency table in matrix form, or a factor object.
# y: a factor object; ignored if 'x' is a matrix.
# alternative: indicates the alternative hypothesis and must be one of '"two.sided"', '"greater"' or '"less"'.  You can specify just the initial letter.  Only used in the 2 by 2 case.
# For detail, see rhelp("fisher.test")
# See also: wilcoxtest, studenttest, padjust, cortest, fisherexacttest, binomtest
# Xiong Jieyi, 27 Oct 2017

function fishertest(arg...; vb::Int=0, warg...)
    @pkgfun(myR, callrw, relem, r2j, rshow)
    o=callrw("fisher.test", arg...; warg...)
    rshow(o; vb=vb)
    r2j(relem(o, "p.value"), NA=NaN, keepvec=false)
end

#[[ binomtest ]]
# p.value, conf.int = binomtest(x, n, p = 0.5,
#            alternative = c("two.sided", "less", "greater"),
#            conf.level = 0.95; vb=0)
# binomtest(682, 682 + 243, p = 3/4) is equal to binomtest([682, 243], p = 3/4).
# Performs an exact test of a simple null hypothesis about the probability of success in a Bernoulli experiment.
# x: number of successes, or a vector of length 2 giving the numbers of successes and failures, respectively.
# n: number of trials; ignored if x has length 2.
# p: hypothesized probability of success.
# alternative: indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter.
# conf.level: confidence level for the returned confidence interval.
# For detail, see rhelp("binom.test")
# See also: wilcoxtest, studenttest, padjust, cortest, fishertest
# Xiong Jieyi, 27 Oct 2017

function binomtest(arg...; vb::Int=0, warg...)
    @pkgfun(myR, callrw, relem, r2j, rshow)
    ro=callrw("binom.test", arg...; warg...)
    rshow(ro; vb=vb)
    (r2j(relem(ro, "p.value"), NA=NaN, keepvec=false), r2j(relem(ro, "conf.int"), NA=NaN, keepvec=true))
end
#}}
