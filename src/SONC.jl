function soncsocp(f,x;ntype="general",itype="Int64",ltype="Int64")
n=length(x)
mon=monomials(f)
coe=coefficients(f)
lsupp=length(mon)
supp=zeros(Int8,n,lsupp)
for i=1:lsupp
    for j=1:n
        supp[j,i]=degree(mon[i],x[j])
    end
end
pos,neg,coe=posneg!(n,supp,coe)
pos=sortslices(pos,dims=2)
if size(pos,2)==n+1
    l=size(neg,2)
    B=zeros(Rational,n,3,1)
    trellis=[k for k=1:n+1]
    if ntype=="standard"
        d=pos[n,2]
        lambuda=zeros(Rational,n+1,1)
        for i=1:l
            lambuda[1]=(d-sum(neg[:,i]))//d
            for k=2:n+1
                lambuda[k]=neg[n+2-k,i]//d
            end
            new=MedSet(neg[:,i],trellis,lambuda,pos,n,itype=itype)
            B=cat(B,new,dims=3)
        end
    else
        A=vcat(pos,ones(Int,1,n+1))
        for i=1:l
            b=vcat(neg[:,i],Int(1))
            lambuda=LinearSolve(A,b,ltype=ltype)
            new=MedSet(neg[:,i],trellis,lambuda,pos,n,itype=itype)
            B=cat(B,new,dims=3)
        end
    end
    B=B[:,:,2:end]
else
    inset,trellis,lambuda,l=simplexcover(pos,neg,n,ltype=ltype)
    B=zeros(Rational,n,3,1)
    for i=1:l
        new=MedSet(inset[:,i],trellis[i],lambuda[i],pos,n,itype=itype)
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

function simplexcover(pos,neg,n;ltype="Int64")
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
    lambuda[k]=LinearSolve(A,b,ltype=ltype)
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
        inset=[inset V[:,1]]
        loc=bfind(pos,lpos,U[:,1],n)
        bary=simsel(V[:,1],pos,loc)
        k=k+1
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
        lambuda[k]=LinearSolve(A,b,ltype=ltype)
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
        inset=[inset V[:,1]]
        loc=bfind(pos,lpos,U[:,1],n)
        bary=simsel(V[:,1],pos,loc)
        k=k+1
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
        lambuda[k]=LinearSolve(A,b,ltype=ltype)
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

function MedSeq(p,q;itype="Int64")
    w=gcd(p,q)
    u=Int(p/w)
    v=Int(q/w)
    if iseven(u)
        if v==u/2
            A=[1;0;2]
        elseif v<Int(u/2)
            a=[Int(u/2);0;u]
            b=MedSeq(Int(u/2),v,itype=itype)
            A=[a b]
        else
            a=[Int(u/2);0;u]
            b=MedSeq(Int(u/2),v-Int(u/2),itype=itype).+Int(u/2)
            A=[a b]
        end
    else
        if iseven(v)
            r,k=tpower(v)
            seq=[v-r*2^(k-i) for i=0:k]
            if itype=="Int64"
                A=zeros(Int64,3,k+1)
            elseif itype=="Int32"
                A=zeros(Int32,3,k+1)
            else
                A=zeros(Int128,3,k+1)
            end
            for i=1:k
                A[:,i]=[seq[i+1];seq[i];v]
            end
            if v==u-r
                A[:,k+1]=[v;v-r;u]
            elseif v<u-r
                A[:,k+1]=[Int((v-r+u)/2);v-r;u]
                b=MedSeq(Int((u+r-v)/2),r,itype=itype).+(v-r)
                A=[A b]
            else
                A[:,k+1]=[Int((v-r+u)/2);v-r;u]
                b=MedSeq(Int((u+r-v)/2),Int((v+r-u)/2),itype=itype).+Int((v+u-r)/2)
                A=[A b]
            end
        else
            A=u.-MedSeq(u,u-v,itype=itype)
        end
    end
    if w>1
        return w*A
    else
        return A
    end
end

function MedSet(inner,trellis,lambuda,pos,n;itype="Int64")
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
        A=MedSeq(p,q,itype=itype)
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

function LinearSolve(A,b;ltype="Int64")
    if ltype=="Int64"
        A=Rational{Int}.(A)
    elseif ltype=="Int32"
        A=Rational{Int32}.(A)
    else
        A=Rational{Int128}.(A)
    end
    n=size(A,2)
    x=zeros(Rational,n,1)
    y=zeros(Rational,n,1)
    F=lu(A)
    b=b[F.p]
    A=F.L[1:n,:]+F.U-Matrix(I,n,n)
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
