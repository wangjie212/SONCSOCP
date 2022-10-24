using DynamicPolynomials
include("D:\\Programs\\SONCSOCP\\SONCSOCP\\src\\SONCSOCP.jl")
using .SONCSOCP

@polyvar x0 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31 x32 x33 x34 x35 x36 x37 x38 x39
x = tuple(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9)
f = include("D:\\Programs\\SONCSOCP\\benchmarks\\general_10_60_300_231.txt")
@time begin
opt = soncsocp(f, x, method="SOCP", alg="GCO", ltype=Int128, itype=Int, solver="Mosek", QUIET=true)
end
println(opt)
@time begin
opt = soncsocp(f, x, method="SOCP", alg="GPT", ltype=Int128, itype=Int, solver="Mosek", QUIET=true)
end
println(opt)
