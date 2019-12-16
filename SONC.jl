using JuMP
using Mosek
using MosekTools

n=2
ms=MSession()
mat"x = sym('x',[1 $n]);
poly=1+x(1)^4*x(2)^4+x(1)^4+x(2)^4-x(1)*x(2)^2-x(1)^2*x(2);
% poly=1+3*x(1)^2*x(2)^6+2*x(1)^6*x(2)^2+6*x(1)^2*x(2)^2-x(1)*x(2)^2-2*x(1)^2*x(2)-3*x(1)^3*x(2)^3;
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
supp=convert(Array{UInt8},supp)

@time begin
opt,sol=soncsocp(n,supp,coe)
end

function soncsocp(n,supp,coe)
lsupp=size(supp,2)
pos,neg,coe=posneg!(n,supp,coe)
pos=sortslices(pos,dims=2)
inset,trellis,lambuda,l=simplexcover(pos,neg,n)
B=MedSet(inset[:,1],trellis[1],lambuda[1],pos,n)
for i=2:l
    new=MedSet(inset[:,i],trellis[i],lambuda[i],pos,n)
    B=cat(B,new,dims=1)
end
num=length(B)
model=Model(with_optimizer(Mosek.Optimizer, QUIET=false))
socp=Array{Any}(undef, num)
supp1=pos
for i=1:num
    supp1=[supp1 B[i][1]]
    socp[i]=@variable(model, [1:3])
    @constraint(model, socp[i] in RotatedSecondOrderCone())
end
supp1=sortslices(supp1,dims=2)
supp1=unique(supp1,dims=2)
lsupp1=size(supp1,2)
cons=Array{Any}(undef, lsupp1)
cons.=AffExpr(0)
for i=1:num
    Locb=bfind(supp1,lsupp1,B[i][2],n)
    cons[Locb]+=2*socp[i][1]
    Locb=bfind(supp1,lsupp1,B[i][3],n)
    cons[Locb]+=socp[i][2]
    Locb=bfind(supp1,lsupp1,B[i][1],n)
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

function posneg!(n,supp,coe)
lo=size(supp,2)
posb=[0]
negb=[0]
for i=1:lo
    bi=supp[:,i]
    if sum(Int[iseven(bi[j]) for j=1:n])==n&&coe[i]>0
       posb=[posb i]
    else
       negb=[negb i]
       if coe[i]>0
           coe[i]=-coe[i]
       end
    end
end
posb=posb[2:end]
negb=negb[2:end]
return supp[:,posb],supp[:,negb],coe
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
        if bary[i]>0
           treb[k]=[treb[k] i]
        end
    end
    treb[k]=treb[k][2:end]
    ltreb=length(treb[k])
    lambuda[k]=zeros(Rational,1,ltreb)
    del=[0]
    for i=1:ltreb
        lambuda[k][i]=rationalize(bary[treb[k][i]],tol=0.0001)
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
            if bary[i]>0
               treb[k]=[treb[k] i]
            end
        end
        treb[k]=treb[k][2:end]
        ltreb=length(treb[k])
        lambuda[k]=zeros(Rational,1,ltreb)
        del=[0]
        for i=1:ltreb
            lambuda[k][i]=rationalize(bary[treb[k][i]],tol=0.0001)
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
            if bary[i]>0
               treb[k]=[treb[k] i]
            end
        end
        treb[k]=treb[k][2:end]
        ltreb=length(treb[k])
        lambuda[k]=zeros(Rational,1,ltreb)
        del=[0]
        for i=1:ltreb
            lambuda[k][i]=rationalize(bary[treb[k][i]],tol=0.0001)
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
        if bary[i]>0
           treb[k]=[treb[k] i]
        end
    end
    treb[k]=treb[k][2:end]
    ltreb=length(treb[k])
    lambuda[k]=zeros(Rational,1,ltreb)
    del=[0]
    for i=1:ltreb
        lambuda[k][i]=rationalize(bary[treb[k][i]],tol=0.0001)
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
    l=length(trellis[1])
    lamb=lambuda[1]
    trel=trellis[1]
    mid=inset[:,1]
    M=[[zeros(Rational,n,1) for i=1:3]]
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
            new=[[mid+A[i]//p*(end1-mid) for i=1:3]]
            M=cat(M,new,dims=3)
        else
            for s=1:lA
                new=[[mid+A[i,s]//p*(end1-mid) for i=1:3]]
                M=cat(M,new,dims=3)
            end
        end
    end
    return M[2:end]
end
