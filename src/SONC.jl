function soncsocp(f, x; ntype="general", itype="Int64", ltype="Int64", QUIET=false, exact=false, tol=1e-5)
    n=length(x)
    mon=monomials(f)
    coe=coefficients(f)
    lsupp=length(mon)
    supp=zeros(Int, n, lsupp)
    for i=1:lsupp, j=1:n
        supp[j,i]=MultivariatePolynomials.degree(mon[i], x[j])
    end
    pos,neg,coe=posneg!(n, supp, coe)
    pos=sortslices(pos, dims=2)
    if itype=="Int32"
        B=Matrix{Rational{Int32}}[]
    else
        B=Matrix{Rational}[]
    end
    if size(pos,2)==n+1
        lneg=size(neg, 2)
        trellis=[k for k=1:n+1]
        if ntype=="standard"
            d=pos[n,2]
            lambuda=zeros(Rational, n+1)
            for i=1:lneg
                lambuda[1]=(d-sum(neg[:,i]))//d
                for k=2:n+1
                    lambuda[k]=neg[n+2-k,i]//d
                end
                new=MedSet(neg[:,i], trellis, lambuda, pos, n, itype=itype)
                append!(B, new)
            end
        else
            A=vcat(pos, ones(Int, 1, n+1))
            for i=1:lneg
                b=vcat(neg[:,i], 1)
                lambuda=LinearSolve(A, b, ltype=ltype)
                new=MedSet(neg[:,i], trellis, lambuda, pos, n, itype=itype)
                append!(B, new)
            end
        end
    else
        inset,trellis,lambuda,numCirc=simplexcover(pos, neg, n, ltype=ltype, itype=itype)
        for i=1:numCirc
            new=MedSet(inset[i], trellis[i], lambuda[i], pos, n, itype=itype)
            append!(B, new)
        end
    end
    num=length(B)
    model=Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    socp=Vector{Vector{VariableRef}}(undef, num)
    for i=1:num
        socp[i]=@variable(model, [1:3])
        @constraint(model, socp[i] in RotatedSecondOrderCone())
    end
    tsupp=[B[i][:, 3] for i=1:num]
    for i=1:size(pos, 2)
        push!(tsupp, pos[:, i])
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp=length(tsupp)
    cons=[AffExpr(0) for i=1:ltsupp]
    Locb=zeros(UInt32, 3, num)
    for i=1:num
        Locb[1,i]=bfind(tsupp, ltsupp, B[i][:,1])
        cons[Locb[1,i]]+=2*socp[i][1]
        Locb[2,i]=bfind(tsupp, ltsupp, B[i][:,2])
        cons[Locb[2,i]]+=socp[i][2]
        Locb[3,i]=bfind(tsupp, ltsupp, B[i][:,3])
        cons[Locb[3,i]]-=2*socp[i][3]
    end
    if exact==false
        bc=zeros(ltsupp)
    else
        bc=zeros(Rational{Int64}, ltsupp)
    end
    for i=1:lsupp
        locb=bfind(tsupp, ltsupp, supp[:,i])
        bc[locb]=coe[i]
    end
    if exact==false
        @constraint(model, cons[2:end].==bc[2:end])
        @variable(model, lower)
        @constraint(model, cons[1]+lower==bc[1])
        @objective(model, Max, lower)
        optimize!(model)
        objv = objective_value(model)
        println("optimum = $objv")
        sol=nothing
    else
        @constraint(model, cons.==bc)
        optimize!(model)
        objv=nothing
        if primal_status(model)==MOI.FEASIBLE_POINT
            sol=zeros(Rational{BigInt}, 3, num)
            for i=1:num
                sol[:,i]=rationalize.(BigInt, value.(socp[i]), tol=tol)
            end
            η=zeros(UInt8, ltsupp)
            val=zeros(Rational{BigInt}, ltsupp)
            for i=1:num
                val[Locb[1,i]]+=2*sol[1,i]
                η[Locb[1,i]]+=1
                val[Locb[2,i]]+=sol[2,i]
                η[Locb[2,i]]+=1
                val[Locb[3,i]]-=2*sol[3,i]
                η[Locb[3,i]]+=1
            end
            for i=1:num
                sol[1,i]-=(val[Locb[1,i]]-bc[Locb[1,i]])/(2*η[Locb[1,i]])
                sol[2,i]-=(val[Locb[2,i]]-bc[Locb[2,i]])/η[Locb[2,i]]
                sol[3,i]+=(val[Locb[3,i]]-bc[Locb[3,i]])/(2*η[Locb[3,i]])
            end
        else
            @warn("Exact decompositon failed!")
            sol=nothing
        end
    end
    return objv,sol,B
end

function posneg!(n, supp, coe)
    lo=size(supp, 2)
    posb=UInt16[]
    negb=UInt16[]
    for i=1:lo
        bi=supp[:,i]
        if sum(Int[iseven(bi[j]) for j=1:n])==n&&coe[i]>0
           posb=push!(posb, i)
        else
           negb=push!(negb, i)
           if coe[i]>0
              coe[i]=-coe[i]
           end
        end
    end
    return supp[:,posb],supp[:,negb],coe
end

function simsel(inner, pos, loc)
    lpos=size(pos, 2)
    model=Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
    @variable(model, x[1:lpos]>=0)
    @constraint(model, pos*x.==inner)
    @constraint(model, sum(x)==1)
    @objective(model, Max, x[loc])
    optimize!(model)
    return value.(x)
end

function simplexcover(pos, neg, n; ltype="Int64", itype="Int64")
    U=pos
    lpos=size(pos, 2)
    V=neg
    k=0
    lneg=size(neg, 2)
    treb=Vector{Vector{UInt16}}(undef, lpos+lneg)
    if itype=="Int32"
        lambuda=Vector{Vector{Rational{Int32}}}(undef, lpos+lneg)
    else
        lambuda=Vector{Vector{Rational}}(undef, lpos+lneg)
    end
    inset=Vector{Int}[]
    while size(U, 2)>0&&size(V, 2)>0
        lU=size(U, 2)
        k+=1
        push!(inset, V[:,1])
        loc=bfind(pos, lpos, U[:,1])
        bary=simsel(V[:,1], pos, loc)
        treb[k]=UInt16[]
        for i=1:lpos
            if bary[i]>0.000001
               push!(treb[k], i)
            end
        end
        ltreb=length(treb[k])
        A=vcat(pos[:,treb[k]], ones(Int, 1, ltreb))
        b=vcat(V[:,1], 1)
        lambuda[k]=LinearSolve(A, b, ltype=ltype)
        del=UInt16[]
        for i=1:ltreb
            bi=pos[:, treb[k][i]]
            loc=bfind(U, lU, bi)
            if loc!=0
               push!(del, loc)
            end
        end
        ldel=length(del)
        Ub=[i for i=1:lU]
        for i=1:ldel
            loc=bfind(Ub, lU, del[i])
            deleteat!(Ub, loc)
            lU=lU-1
        end
        U=U[:, Ub]
        V=V[:, 2:end]
    end
    if size(V, 2)>0
        while size(V, 2)>0
            if size(U, 2)==0
                U=pos
            end
            lU=size(U, 2)
            push!(inset, V[:,1])
            loc=bfind(pos, lpos, U[:,1])
            bary=simsel(V[:,1], pos, loc)
            k+=1
            treb[k]=UInt16[]
            for i=1:lpos
                if bary[i]>0.000001
                   push!(treb[k], i)
                end
            end
            ltreb=length(treb[k])
            A=vcat(pos[:,treb[k]], ones(Int, 1, ltreb))
            b=vcat(V[:,1], 1)
            lambuda[k]=LinearSolve(A, b, ltype=ltype)
            del=UInt16[]
            for i=1:ltreb
                bi=pos[:, treb[k][i]]
                loc=bfind(U, lU, bi)
                if loc!=0
                   push!(del, loc)
                end
            end
            ldel=length(del)
            Ub=[i for i=1:lU]
            for i=1:ldel
                loc=bfind(Ub, lU, del[i])
                deleteat!(Ub, loc)
                lU=lU-1
            end
            U=U[:, Ub]
            V=V[:, 2:end]
        end
    else
        while size(U, 2)>0
            if size(V, 2)==0
                V=neg
            end
            lU=size(U, 2)
            push!(inset, V[:,1])
            loc=bfind(pos, lpos, U[:,1])
            bary=simsel(V[:,1], pos, loc)
            k+=1
            treb[k]=UInt16[]
            for i=1:lpos
                if bary[i]>0.000001
                    push!(treb[k], i)
                end
            end
            ltreb=length(treb[k])
            A=vcat(pos[:, treb[k]], ones(Int, 1, ltreb))
            b=vcat(V[:,1], 1)
            lambuda[k]=LinearSolve(A, b, ltype=ltype)
            del=UInt16[]
            for i=1:ltreb
                bi=pos[:,treb[k][i]]
                loc=bfind(U, lU, bi)
                if loc!=0
                   push!(del, loc)
                end
            end
            ldel=length(del)
            Ub=[i for i=1:lU]
            for i=1:ldel
                loc=bfind(Ub, lU, del[i])
                deleteat!(Ub, loc)
                lU=lU-1
            end
            U=U[:, Ub]
            V=V[:, 2:end]
        end
    end
    return inset,treb[1:k],lambuda[1:k],k
end

function bfind(A, l, a)
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        if ndims(A)==1
            b=A[mid]
        else
            b=A[:, mid]
        end
        if a==b
           return mid
        elseif a>b
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

function MedSeq(p, q; itype="Int64")
    if itype=="Int128"
        w=Int128(gcd(p, q))
        u=Int128(p/w)
        v=Int128(q/w)
    elseif itype=="Int32"
        w=Int32(gcd(p, q))
        u=Int32(p/w)
        v=Int32(q/w)
    else
        w=Int(gcd(p, q))
        u=Int(p/w)
        v=Int(q/w)
    end
    if iseven(u)
        if v==u/2
            if itype=="Int128"
                A=Int128[1;0;2]
            elseif itype=="Int32"
                A=Int32[1;0;2]
            else
                A=[1;0;2]
            end
        else
            a=[Int(u/2);0;u]
            if v<Int(u/2)
                b=MedSeq(Int(u/2), v, itype=itype)
            else
                b=MedSeq(Int(u/2), v-Int(u/2), itype=itype).+Int(u/2)
            end
            if itype=="Int128"
                A=Int128[a b]
            elseif itype=="Int32"
                A=Int32[a b]
            else
                A=[a b]
            end
        end
    else
        if iseven(v)
            r,k=tpower(v)
            seq=[v-r*2^(k-i) for i=0:k]
            if itype=="Int64"
                A=zeros(Int64, 3, k+1)
            elseif itype=="Int32"
                A=zeros(Int32, 3, k+1)
            else
                A=zeros(Int128, 3, k+1)
            end
            for i=1:k
                A[:,i]=[seq[i+1];seq[i];v]
            end
            if v==u-r
                A[:,k+1]=[v;v-r;u]
            else
                A[:,k+1]=[Int((v-r+u)/2);v-r;u]
                if v<u-r
                    b=MedSeq(Int((u+r-v)/2), r, itype=itype).+(v-r)
                else
                    b=MedSeq(Int((u+r-v)/2), Int((v+r-u)/2), itype=itype).+Int((v+u-r)/2)
                end
                if itype=="Int128"
                    A=Int128[A b]
                elseif itype=="Int32"
                    A=Int32[A b]
                else
                    A=[A b]
                end
            end
        else
            A=u.-MedSeq(u, u-v, itype=itype)
        end
    end
    if w>1
        return w*A
    else
        return A
    end
end

function MedSet(inner, trellis, lambuda, pos, n; itype="Int64")
    l=length(trellis)
    lamb=lambuda
    trel=trellis
    mid=inner
    if itype=="Int32"
        M=Matrix{Rational{Int32}}[]
    else
        M=Matrix{Rational}[]
    end
    for j=1:l-1
        p,k=findmin([denominator(lamb[i]) for i=1:l+1-j])
        q=numerator(lamb[k])
        end1=pos[:, trel[k]]
        treb=[i for i=1:l+1-j]
        deleteat!(treb, k)
        trel=trel[treb]
        lamb=lamb[treb]*(p//(p-q))
        mid=p//(p-q)*mid-q//(p-q)*end1
        A=MedSeq(p, q, itype=itype)
        lA=size(A, 2)
        if lA==1
            new=[mid+A[2]//p*(end1-mid) mid+A[3]//p*(end1-mid) mid+A[1]//p*(end1-mid)]
            push!(M, new)
        else
            for s=1:lA
                new=[mid+A[2,s]//p*(end1-mid) mid+A[3,s]//p*(end1-mid) mid+A[1,s]//p*(end1-mid)]
                push!(M, new)
            end
        end
    end
    return M
end

function LinearSolve(A, b; ltype="Int64")
    n=size(A, 2)
    if ltype=="Int32"
        A=Rational{Int32}.(A)
        b=Rational{Int32}.(b)
        x=zeros(Rational{Int32}, n)
        y=zeros(Rational{Int32}, n)
    else
        if ltype=="Int128"
            A=Rational{Int128}.(A)
            b=Rational{Int128}.(b)
        else
            A=Rational.(A)
            b=Rational.(b)
        end
        x=zeros(Rational, n)
        y=zeros(Rational, n)
    end
    F=lu(A)
    b=b[F.p]
    A=F.L[1:n,:]+F.U-Matrix(I, n, n)
    for i=1:n
        alpha=0
        for k=1:i
            alpha+=A[i,k]*y[k]
        end
        y[i]=b[i]-alpha
    end
    for i=n:(-1):1
        alpha=0
        for k=(i+1):n
            alpha+=A[i,k]*x[k]
        end
        x[i]=(y[i]-alpha)//A[i,i]
    end
    return x
end
