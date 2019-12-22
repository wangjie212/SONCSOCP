using JuMP
using Mosek
using MosekTools
using LinearAlgebra
using TypedPolynomials
using MultivariatePolynomials

n=2
@polyvar x[1:2]
f=1+x[1]^4*x[2]^4+x[1]^4+x[2]^4-x[1]*x[2]^2-x[1]^2*x[2]+x[1]*x[2]
mon=monomials(f)
coe=coefficients(f)
lm=length(mon)
supp=zeros(Int8,n,lm)
for i=1:lm
    for j=1:n
        supp[j,i]=degree(mon[i],x[j])
    end
end

n=2
ms=MSession()
mat"% x = sym('x',[1 $n]);
poly=1+x(1)^4*x(2)^4+x(1)^4+x(2)^4-x(1)*x(2)^2-x(1)^2*x(2)+x(1)*x(2);
% poly=1+3*x(1)^2*x(2)^6+2*x(1)^6*x(2)^2+6*x(1)^2*x(2)^2-x(1)*x(2)^2-2*x(1)^2*x(2)-3*x(1)^3*x(2)^3;
% poly=x(1)^4*x(2)^2+x(1)^2*x(2)^4+1-3*x(1)^2*x(2)^2;
% syms x0 x1 x2 x3 x4 x5 x6 x7 x8 x9;
% x=[x0 x1 x2 x3 x4 x5 x6 x7 x8 x9];
% load E:\\Programs\\sonc\\example_4_2.mat;
[coe, terms] = coeffs(poly,x);
lt=length(terms);
supp=zeros($n,lt);
for i=1:lt
    for j=1:$n
        supp(j,i)=feval(symengine,'degree',terms(i),x(j));
    end
end
coe=double(coe)"
coe=jarray(get_mvariable(ms,:coe))
supp=jarray(get_mvariable(ms,:supp))
supp=convert(Array{Int8},supp)

@time begin
opt,sol=soncsocp(n,supp,coe)
end

function soncsocp(n,supp,coe)
lsupp=size(supp,2)
pos,neg=posneg(n,supp,coe)
if size(pos,2)==n+1
    l=size(neg,2)
    A=vcat(pos,ones(Int,1,n+1))
    B=zeros(Rational,n,3,1)
    trellis=[k for k=1:n+1]
    for i=1:l
        b=vcat(neg[:,i],Int(1))
        lambuda=LinearSolve(A,b)
        new=MedSet(neg[:,i],trellis,lambuda,pos,n)
        B=cat(B,new,dims=3)
    end
    B=B[:,:,2:end]
else
    pos=sortslices(pos,dims=2)
    inset,trellis,lambuda,l=simplexcover(pos,neg,n)
    B=zeros(Rational,n,3,1)
    for i=1:l
        new=MedSet(inset[:,i],trellis[i],lambuda[i],pos,n)
        B=cat(B,new,dims=3)
    end
    B=B[:,:,2:end]
end
num=size(B,3)
model=Model(with_optimizer(Mosek.Optimizer, QUIET=false))
socp=Array{Any}(undef, num)
supp1=[pos B[:,1,1:end]]
for i=1:num
    socp[i]=@variable(model, [1:3])
    @constraint(model, socp[i] in RotatedSecondOrderCone())
end
supp1=sortslices(supp1,dims=2)
supp1=unique(supp1,dims=2)
lsupp1=size(supp1,2)
cons=Array{Any}(undef, lsupp1)
cons.=AffExpr(0)
for i=1:num
    Locb=bfind(supp1,lsupp1,B[:,2,i],n)
    cons[Locb]+=2*socp[i][1]
    Locb=bfind(supp1,lsupp1,B[:,3,i],n)
    cons[Locb]+=socp[i][2]
    Locb=bfind(supp1,lsupp1,B[:,1,i],n)
    cons[Locb]-=2*socp[i][3]
end
bc=zeros(1,lsupp1)
for i=1:lsupp
    Locb=bfind(supp1,lsupp1,supp[:,i],n)
    bc[Locb]=coe[i]
end
@constraint(model, cons[2:end].==bc[2:end])
@variable(model, lower)
@constraint(model, cons[1]+lower==bc[1])
@objective(model, Max, lower)
optimize!(model)
status=termination_status(model)
if status == MOI.OPTIMAL
   objv = objective_value(model)
   println("optimum = $objv")
   sol=zeros(3,num)
   for i=1:num
       sol[:,i]=value.(socp[i])
   end
else
   objv = objective_value(model)
   println("$status")
   println("optimum = $objv")
   sol=zeros(3,num)
   for i=1:num
       sol[:,i]=value.(socp[i])
   end
end
return objv,sol
end

function posneg(n,supp,coe)
lo=size(supp,2)
posb=[0]
negb=[0]
for i=1:lo
    bi=supp[:,i]
    if sum(Int[iseven(bi[j]) for j=1:n])==n&&coe[i]>0
       posb=[posb i]
    else
       negb=[negb i]
    end
end
posb=posb[2:end]
negb=negb[2:end]
return supp[:,posb],supp[:,negb]
end

function simsel(inner,pos,loc)
    lpos=size(pos,2)
    model=Model(with_optimizer(Mosek.Optimizer, QUIET=true))
    @variable(model, x[1:lpos]>=0)
    @constraint(model, pos*x.==inner)
    @constraint(model, sum(x)==1)
    @objective(model, Max, x[loc])
    optimize!(model)
    return value.(x)
end

function simplexcover(pos,neg,n)
U=pos
lpos=size(pos,2)
V=neg
k=0
lneg=size(neg,2)
treb=Array{Any}(undef,lpos+lneg)
lambuda=Array{Any}(undef,lpos+lneg)
inset=zeros(UInt8,n,1)
while size(U,2)>0&&size(V,2)>0
    lU=size(U,2)
    k=k+1
    inset=[inset V[:,1]]
    loc=bfind(pos,lpos,U[:,1],n)
    bary=simsel(V[:,1],pos,loc)
    treb[k]=[0]
    for i=1:lpos
        if bary[i]>0.000001
           treb[k]=[treb[k] i]
        end
    end
    treb[k]=treb[k][2:end]
    ltreb=length(treb[k])
    A=vcat(pos[:,treb[k]],ones(Int,1,ltreb))
    b=vcat(V[:,1],Int(1))
    lambuda[k]=LinearSolve(A,b)
    del=[0]
    for i=1:ltreb
        bi=pos[:,treb[k][i]]
        loc=bfind(U,lU,bi,n)
        if loc!=0
           del=[del loc]
        end
    end
    del=del[2:end]
    ldel=length(del)
    Ub=[i for i=1:lU]
    for i=1:ldel
        loc=lbfind(Ub,lU,del[i])
        deleteat!(Ub,loc)
        lU=lU-1
    end
    U=U[:,Ub]
    V=V[:,2:end]
end
if size(V,2)>0
    while size(V,2)>0
        if size(U,2)==0
            U=pos
        end
        lU=size(U,2)
        k=k+1
        inset=[inset V[:,1]]
        loc=bfind(pos,lpos,U[:,1],n)
        bary=simsel(V[:,1],pos,loc)
        treb[k]=[0]
        for i=1:lpos
            if bary[i]>0.000001
               treb[k]=[treb[k] i]
            end
        end
        treb[k]=treb[k][2:end]
        ltreb=length(treb[k])
        A=vcat(pos[:,treb[k]],ones(Int,1,ltreb))
        b=vcat(V[:,1],Int(1))
        lambuda[k]=LinearSolve(A,b)
        del=[0]
        for i=1:ltreb
            bi=pos[:,treb[k][i]]
            loc=bfind(U,lU,bi,n)
            if loc!=0
               del=[del loc]
            end
        end
        del=del[2:end]
        ldel=length(del)
        Ub=[i for i=1:lU]
        for i=1:ldel
            loc=lbfind(Ub,lU,del[i])
            deleteat!(Ub,loc)
            lU=lU-1
        end
        U=U[:,Ub]
        V=V[:,2:end]
    end
else
    while size(U,2)>0
        if size(V,2)==0
            V=neg
        end
        lU=size(U,2)
        k=k+1
        inset=[inset V[:,1]]
        loc=bfind(pos,lpos,U[:,1],n)
        bary=simsel(V[:,1],pos,loc)
        treb[k]=[0]
        for i=1:lpos
            if bary[i]>0.000001
               treb[k]=[treb[k] i]
            end
        end
        treb[k]=treb[k][2:end]
        ltreb=length(treb[k])
        A=vcat(pos[:,treb[k]],ones(Int,1,ltreb))
        b=vcat(V[:,1],Int(1))
        lambuda[k]=LinearSolve(A,b)
        del=[0]
        for i=1:ltreb
            bi=pos[:,treb[k][i]]
            loc=bfind(U,lU,bi,n)
            if loc!=0
               del=[del loc]
            end
        end
        del=del[2:end]
        ldel=length(del)
        Ub=[i for i=1:lU]
        for i=1:ldel
            loc=lbfind(Ub,lU,del[i])
            deleteat!(Ub,loc)
            lU=lU-1
        end
        U=U[:,Ub]
        V=V[:,2:end]
    end
end
return inset[:,2:end],treb[1:k],lambuda[1:k],k
end

function simplexcover2(pos,neg,n)
U=pos
lpos=size(pos,2)
V=neg
k=0
lneg=size(neg,2)
treb=Array{Any}(undef,lpos+lneg)
lambuda=Array{Any}(undef,lpos+lneg)
inset=zeros(UInt8,n,1)
while size(V,2)>0
    if size(U,2)==0
        U=pos
    end
    lU=size(U,2)
    k=k+1
    inset=[inset V[:,1]]
    loc=bfind(pos,lpos,U[:,1],n)
    bary=simsel(V[:,1],pos,loc)
    treb[k]=[0]
    for i=1:lpos
        if bary[i]>0.000001
           treb[k]=[treb[k] i]
        end
    end
    treb[k]=treb[k][2:end]
    ltreb=length(treb[k])
    A=vcat(pos[:,treb[k]],ones(Int,1,ltreb))
    b=vcat(V[:,1],Int(1))
    lambuda[k]=LinearSolve(A,b)
    del=[0]
    for i=1:ltreb
        bi=pos[:,treb[k][i]]
        loc=bfind(U,lU,bi,n)
        if loc!=0
           del=[del loc]
        end
    end
    del=del[2:end]
    ldel=length(del)
    Ub=[i for i=1:lU]
    for i=1:ldel
        loc=lbfind(Ub,lU,del[i])
        deleteat!(Ub,loc)
        lU=lU-1
    end
    U=U[:,Ub]
    V=V[:,2:end]
end
return inset[:,2:end],treb[1:k],lambuda[1:k],k
end

function lbfind(A,l,a)
    if l==0
        return 0
    end
    low=Int(1)
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        if A[mid]==a
           return mid
       elseif A[mid]<a
           low=mid+1
        else
           high=mid-1
        end
    end
    return 0
end

function comp(a,b,n)
    i=Int(1)
    while i<=n
          if a[i]<b[i]
             return -1
          elseif a[i]>b[i]
             return 1
          else
             i+=1
          end
    end
    if i==n+1
       return 0
    end
end

function bfind(A,l,a,n)
    if l==0
        return 0
    end
    low=Int(1)
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        order=comp(A[:,mid],a,n)
        if order==0
           return mid
        elseif order<0
           low=mid+1
        else
           high=mid-1
        end
    end
    return 0
end

function tpower(v)
    k=0
    while iseven(v)
        k+=1
        v=Int(v/2)
    end
    return v,k
end

function MedSeq(p,q)
    w=gcd(p,q)
    u=Int(p/w)
    v=Int(q/w)
    if iseven(u)
        if v==u/2
            A=[1;0;2]
        elseif v<Int(u/2)
            a=[Int(u/2);0;u]
            b=MedSeq(Int(u/2),v)
            A=[a b]
        else
            a=[Int(u/2);0;u]
            b=MedSeq(Int(u/2),v-Int(u/2)).+Int(u/2)
            A=[a b]
        end
    else
        if iseven(v)
            r,k=tpower(v)
            seq=[v-r*2^(k-i) for i=0:k]
            A=zeros(Int,3,k+1)
            for i=1:k
                A[:,i]=[seq[i+1];seq[i];v]
            end
            if v==u-r
                A[:,k+1]=[v;v-r;u]
            elseif v<u-r
                A[:,k+1]=[Int((v-r+u)/2);v-r;u]
                b=MedSeq(Int((u+r-v)/2),r).+(v-r)
                A=[A b]
            else
                A[:,k+1]=[Int((v-r+u)/2);v-r;u]
                b=MedSeq(Int((u+r-v)/2),Int((v+r-u)/2)).+Int((v+u-r)/2)
                A=[A b]
            end
        else
            A=u.-MedSeq(u,u-v)
        end
    end
    if w>1
        return w*A
    else
        return A
    end
end

function MedSet(inner,trellis,lambuda,pos,n)
    l=length(trellis)
    lamb=lambuda
    trel=trellis
    mid=inner
    M=zeros(Rational,n,3,1)
    for j=1:l-1
        p,k=findmin([denominator(lamb[i]) for i=1:l+1-j])
        q=numerator(lamb[k])
        end1=pos[:,trel[k]]
        treb=[i for i=1:l+1-j]
        deleteat!(treb,k)
        trel=trel[treb]
        lamb=lamb[treb]*(p//(p-q))
        mid=p//(p-q)*mid-q//(p-q)*end1
        A=MedSeq(p,q)
        lA=size(A,2)
        if lA==1
            new=[mid+A[1]//p*(end1-mid) mid+A[2]//p*(end1-mid) mid+A[3]//p*(end1-mid)]
            M=cat(M,new,dims=3)
        else
            for s=1:lA
                new=[mid+A[1,s]//p*(end1-mid) mid+A[2,s]//p*(end1-mid) mid+A[3,s]//p*(end1-mid)]
                M=cat(M,new,dims=3)
            end
        end
    end
    return M[:,:,2:end]
end

function LinearSolve(A,b)
    A=Rational{Int}.(A)
    m=length(b)
    n=size(A,2)
    x=zeros(Rational,n,1)
    y=zeros(Rational,n,1)
    F=lu(A)
    b=b[F.p]
    if n<m
       A=F.L[1:n,:]+F.U-Matrix(I,n,n)
    else
       A=F.L+F.U-Matrix(I,n,n)
    end
    for i=1:n
        alpha=0
        for k=1:i
            alpha=alpha+A[i,k]*y[k]
        end
        y[i]=b[i]-alpha
    end
    for i=n:(-1):1
        alpha=0
        for k=(i+1):n
            alpha=alpha+A[i,k]*x[k]
        end
        x[i]=(y[i]-alpha)//A[i,i]
    end
    return x
end

#function MedSet2(inner,trellis,lambuda,pos,n)
#    l=length(trellis)
#    lamb=lambuda
#    trel=trellis
#    mid=inner
#    M=zeros(Rational,n,3,1)
#    i=1
#    while i<=l&&isodd(denominator(lamb[i]))
#          i+=1
#    end
#    if i<=l
#       tp=zeros(Int,1,l)
#       for k=1:l
#           v,tp[k]=tpower(denominator(lamb[k]))
#       end
#       v,k=findmax([tp[i] for i=1:l])
#       p=denominator(lamb[k])
#       q=numerator(lamb[k])
#       end1=pos[:,trel[k]]
#       treb=[i for i=1:l]
#       deleteat!(treb,k)
#       l-=1
#       trel=trel[treb]
#       lamb=lamb[treb]*(p//(p-q))
#       mid=p//(p-q)*mid-q//(p-q)*end1
#       if q<=p/2
#          new=[inner mid 2*inner-mid]
#          A=MedSeq(p,2*q)
#          M=cat(M,new,dims=3)
#          lA=size(A,2)
#          if lA==1
#              new=[mid+A[1]//p*(end1-mid) mid+A[2]//p*(end1-mid) mid+A[3]//p*(end1-mid)]
#              M=cat(M,new,dims=3)
#          else
#              for s=1:lA
#                  new=[mid+A[1,s]//p*(end1-mid) mid+A[2,s]//p*(end1-mid) mid+A[3,s]//p*(end1-mid)]
#                  M=cat(M,new,dims=3)
#              end
#          end
#       else
#          new=[inner 2*inner-end1 end1]
#          M=cat(M,new,dims=3)
#          A=MedSeq(p,2*(p-q))
#          lA=size(A,2)
#          if lA==1
#              new=[end1+A[1]//p*(mid-end1) end1+A[2]//p*(mid-end1) end1+A[3]//p*(mid-end1)]
#              M=cat(M,new,dims=3)
#          else
#              for s=1:lA
#                  new=[end1+A[1,s]//p*(mid-end1) end1+A[2,s]//p*(mid-end1) end1+A[3,s]//p*(mid-end1)]
#                  M=cat(M,new,dims=3)
#              end
#          end
#       end
#    end
#    e=1
#    while e>0
#        even=[k for k=1:l]
#        even=even[[iseven(numerator(lamb[k])) for k=1:l]]
#        if length(even)>0
#            p,k=findmin(denominator.(lamb[even]))
#            q=numerator(lamb[even[k]])
#            end1=pos[:,trel[even[k]]]
#            treb=[i for i=1:l]
#            deleteat!(treb,even[k])
#            l-=1
#            trel=trel[treb]
#            lamb=lamb[treb]*(p//(p-q))
#            mid=p//(p-q)*mid-q//(p-q)*end1
#            A=MedSeq(p,q)
#            lA=size(A,2)
#            if lA==1
#                new=[mid+A[1]//p*(end1-mid) mid+A[2]//p*(end1-mid) mid+A[3]//p*(end1-mid)]
#                M=cat(M,new,dims=3)
#            else
#                for s=1:lA
#                    new=[mid+A[1,s]//p*(end1-mid) mid+A[2,s]//p*(end1-mid) mid+A[3,s]//p*(end1-mid)]
#                    M=cat(M,new,dims=3)
#                end
#            end
#        else
#           break
#        end
#    end
#    if l>1
#        q=numerator(lamb[1])*denominator(lamb[2])
#        p=q+numerator(lamb[2])*denominator(lamb[1])
#        mid1=mid+lamb[2]*(pos[:,trel[1]]-pos[:,trel[2]])
#        mid2=mid+lamb[1]*(pos[:,trel[2]]-pos[:,trel[1]])
#        A=MedSeq(p,q)
#        lA=size(A,2)
#        if lA==1
#            new=[mid2+A[1]//p*(mid1-mid2) mid2+A[2]//p*(mid1-mid2) mid2+A[3]//p*(mid1-mid2)]
#            M=cat(M,new,dims=3)
#        else
#            for s=1:lA
#                new=[mid2+A[1,s]//p*(mid1-mid2) mid2+A[2,s]//p*(mid1-mid2) mid2+A[3,s]//p*(mid1-mid2)]
#                M=cat(M,new,dims=3)
#            end
#        end
#        treb1=[i for i=1:l]
#        treb2=[i for i=1:l]
#        deleteat!(treb1,2)
#        deleteat!(treb2,1)
#        trel1=trel[treb1]
#        trel2=trel[treb2]
#        lamb1=lamb[2:end]
#        lamb1[1]+=lamb[1]
#        lamb2=lamb1
#        M1=MedSet2(mid1,trel1,lamb1,pos,n)
#        M2=MedSet2(mid2,trel2,lamb2,pos,n)
#        M=cat(M,M1,dims=3)
#        M=cat(M,M2,dims=3)
#    end
#    return M[:,:,2:end]
#end

#function soncsocp2(n,supp,coe)
#lsupp=size(supp,2)
#pos,neg,label=posneg(n,supp,coe)
#if size(pos,2)==n+1
#    l=size(neg,2)
#    A=vcat(pos,ones(Int,1,n+1))
#    B=zeros(Rational,n,3,1)
#    trellis=[k for k=1:n+1]
#    for i=1:l
#        b=vcat(neg[:,i],Int(1))
#        lambuda=LinearSolve(A,b)
#        if i∈label
#           new=MedSet(neg[:,i],trellis,lambuda,pos,n)
#        else
#           new=MedSet2(neg[:,i],trellis,lambuda,pos,n)
#        end
#        B=cat(B,new,dims=3)
#    end
#    B=B[:,:,2:end]
#else
#    pos=sortslices(pos,dims=2)
#    inset,trellis,lambuda,l=simplexcover(pos,neg,n)
#    B=zeros(Rational,n,3,1)
#    for i=1:l
#        if i∈label
#           new=MedSet(inset[:,i],trellis[i],lambuda[i],pos,n)
#        else
#           new=MedSet2(inset[:,i],trellis[i],lambuda[i],pos,n)
#        end
#        B=cat(B,new,dims=3)
#    end
#    B=B[:,:,2:end]
#end
#num=size(B,3)
#model=Model(with_optimizer(Mosek.Optimizer, QUIET=false))
#socp=Array{Any}(undef, num)
#supp1=[pos B[:,1,1:end]]
#for i=1:num
#    socp[i]=@variable(model, [1:3])
#    @constraint(model, socp[i] in RotatedSecondOrderCone())
#end
#supp1=sortslices(supp1,dims=2)
#supp1=unique(supp1,dims=2)
#lsupp1=size(supp1,2)
##cons=Array{Any}(undef, lsupp1)
#cons.=AffExpr(0)
#for i=1:num
#    Locb=bfind(supp1,lsupp1,B[:,2,i],n)
#    cons[Locb]+=2*socp[i][1]
#    Locb=bfind(supp1,lsupp1,B[:,3,i],n)
#    cons[Locb]+=socp[i][2]
#    Locb=bfind(supp1,lsupp1,B[:,1,i],n)
#    cons[Locb]-=2*socp[i][3]
#end
#bc=zeros(1,lsupp1)
#for i=1:lsupp
#    Locb=bfind(supp1,lsupp1,supp[:,i],n)
#    bc[Locb]=coe[i]
#end
#@constraint(model, cons[2:end].==bc[2:end])
#@variable(model, lower)
#@constraint(model, cons[1]+lower==bc[1])
#@objective(model, Max, lower)
#optimize!(model)
#status=termination_status(model)
#if status == MOI.OPTIMAL
#   objv = objective_value(model)
#   println("optimum = $objv")
#   sol=zeros(3,num)
#   for i=1:num
#       sol[:,i]=value.(socp[i])
#   end
#else
#   objv = objective_value(model)
#   println("$status")
#   println("optimum = $objv")
#   sol=zeros(3,num)
#   for i=1:num
#       sol[:,i]=value.(socp[i])
#   end
#end
#return objv,sol
#end
