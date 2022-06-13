function soncsocp(f, x; method="DPCP", alg="GPT", solver="Mosek", ntype="general", itype=Int, QUIET=false)
    n = length(x)
    mon = monomials(f)
    coe = coefficients(f)
    lsupp = length(mon)
    supp = zeros(UInt16, n, lsupp)
    for i = 1:lsupp, j = 1:n
        supp[j,i] = MultivariatePolynomials.degree(mon[i], x[j])
    end
    pos,neg,coe = posneg(n, supp, coe)
    if method == "DPCP"
        itype = Float64
    end
    if ntype != "general"
        inset = neg
        λ = Vector{Vector{itype}}(undef, length(neg))
        trellis = Vector{Vector{UInt16}}(undef, length(neg))
        for i = 1:length(neg)
            trellis[i] = pos
            if ntype == "standard"
                if method == "SOCP"
                    λ[i] = [supp[:,neg[i]];supp[1,1]-sum(supp[:,neg[i]])]
                else
                    λ[i] = [supp[:,neg[i]];supp[1,1]-sum(supp[:,neg[i]])]//supp[1,1]
                end
            else
                b = vcat(supp[:,neg[i]], UInt16(1))
                A = vcat(supp[:,pos], ones(UInt16, 1, length(pos)))
                lamb = LinearSolve(A, b, itype=itype)
                if method == "SOCP"
                    λ[i] = lamb .* lcm(denominator.(lamb))
                else
                    λ[i] = lamb
                end
            end
        end
    else
        inset,trellis,λ = simplexcover(supp, pos, neg, method=method, itype=itype)
    end
    if solver == "Mosek"
        model = Model(optimizer_with_attributes(Mosek.Optimizer))
    else
        model = Model(optimizer_with_attributes(ECOS.Optimizer))
    end
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    circ = Vector{Vector{VariableRef}}(undef, length(λ))
    cons = [AffExpr(0) for i=1:lsupp]
    for i = 1:length(λ)
        m = length(trellis[i])
        if method == "SOCP"
            if alg == "GPT"
                conf = GreedyPowertwo(λ[i])
            else
                conf = GreedyCommone(λ[i])
            end
            circ[i] = @variable(model, [1:m+length(conf)])
            for j = 1:length(conf)
                @constraint(model, [2circ[i][conf[j][1]],circ[i][conf[j][2]],2circ[i][conf[j][3]]] in RotatedSecondOrderCone())
            end
            sl = sum(λ[i])
            factor = prod((λ[i][j]/sl)^(λ[i][j]/sl) for j=1:length(λ[i]))
            cons[inset[i]] -= circ[i][m+1]/factor
        else
            circ[i] = @variable(model, [1:3(m-1)])
            for j = 1:m-2
                s = 1 - λ[i][end]
                λ[i] = λ[i][1:end-1]/s
                @constraint(model, [circ[i][2m+j-1];circ[i][m-j+1];circ[i][m+j]] in MOI.DualPowerCone(s))
                @constraint(model, circ[i][2m+j-1]==circ[i][m+j+1])
            end
            @constraint(model, [circ[i][1];circ[i][2];circ[i][2m-1]] in MOI.DualPowerCone(λ[i][1]))
            cons[inset[i]] -= circ[i][m+1]
        end
        cons[trellis[i]] += circ[i][1:m]
    end
    @variable(model, lower)
    cons[end] += lower
    @constraint(model, cons.==coe)
    @objective(model, Max, lower)
    optimize!(model)
    status = termination_status(model)
    objv = objective_value(model)
    if status != MOI.OPTIMAL
       println("termination status: $status")
    end
    println("optimum = $objv")
    return objv
end

function posneg(n, supp, coe)
    posb = UInt16[]
    negb = UInt16[]
    for i = 1:size(supp,2)
        if coe[i] > 0 && all(j->iseven(supp[j,i]), 1:n)
           posb = push!(posb, i)
        else
           negb = push!(negb, i)
           if coe[i] > 0
              coe[i] = -coe[i]
           end
        end
    end
    return posb,negb,coe
end

function simsel(inner, pos, loc)
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
    @variable(model, x[1:size(pos,2)]>=0)
    @constraint(model, pos*x.==inner)
    @constraint(model, sum(x)==1)
    @objective(model, Max, x[loc])
    optimize!(model)
    return value.(x)
end

function simplexcover(supp, pos, neg; method="SOCP", itype=Int)
    treb = Vector{UInt16}[]
    λ = Vector{itype}[]
    inset = UInt16[]
    U = pos
    V = neg
    while length(U) > 0 && length(V) > 0
        U,V = NewCirc!(supp, pos, U, V, inset, treb, λ, method=method, itype=itype)
    end
    if length(V) > 0
        while length(V) > 0
            if length(U) == 0
                U = pos
            end
            U,V = NewCirc!(supp, pos, U, V, inset, treb, λ, method=method, itype=itype)
        end
    else
        while length(U) > 0
            if length(V) == 0
                V = neg
            end
            U,V = NewCirc!(supp, pos, U, V, inset, treb, λ, method=method, itype=itype)
        end
    end
    return inset,treb,λ
end

function NewCirc!(supp, pos, U, V, inset, treb, λ; method="SOCP", itype=Int)
    push!(inset, V[1])
    loc = bfind(pos, length(pos), U[1])
    bary = simsel(supp[:,V[1]], supp[:,pos], loc)
    push!(treb, pos[[bary[i]>=1e-6 for i=1:length(pos)]])
    if method == "SOCP"
        A = vcat(supp[:,treb[end]], ones(UInt16, 1, length(treb[end])))
        b = vcat(supp[:,V[1]], UInt16(1))
        lamb = LinearSolve(A, b, itype=itype)
        push!(λ, lamb .* lcm(denominator.(lamb)))
    else
        push!(λ, bary[[bary[i]>=1e-6 for i=1:length(pos)]])
    end
    U = setdiff(U, treb[end])
    V = V[2:end]
    return U,V
end

function bfind(A, l, a)
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        b = A[mid]
        if a == b
           return mid
        elseif a > b
           low = mid + 1
        else
           high = mid - 1
        end
    end
    return 0
end

function LinearSolve(A, b; itype=Int)
    n = size(A, 2)
    A = Rational{itype}.(A)
    b = Rational{itype}.(b)
    x = zeros(Rational{itype}, n)
    y = zeros(Rational{itype}, n)
    F = lu(A)
    b = b[F.p]
    A = F.L[1:n,:] + F.U - Matrix(I, n, n)
    for i = 1:n
        alpha = 0
        for k = 1:i
            alpha += A[i,k]*y[k]
        end
        y[i] = b[i] - alpha
    end
    for i = n:(-1):1
        alpha = 0
        for k = i+1:n
            alpha += A[i,k]*x[k]
        end
        x[i] = (y[i]-alpha)//A[i,i]
    end
    return x
end
