# SONCSOCP
SONCSOCP is an uncontrained sparse polynomial optimization tool based on SONC decompositions.

To use SONCSOCP, run
```Julia
pkg> add https://github.com/wangjie212/SONCSOCP
 ```

## Dependencies
- Julia
- MOSEK

SONCSOCP has been tested on WINDOW 10, Julia 1.2.0, and MOSEK 8.1.
## Usage
SONCSOCP computes a lower bound for the unconstrained polynomial optimization problem:
```
Inf{f(x): x\in R^n}
```
where f is a polynomial with variables x1,...,xn.

Taking f=1+x1^4+x2^4-x1\*x2 as an example, to compute a SONC lower bound for f, run
```Julia
julia> using SONCSOCP
julia> using TypedPolynomials
julia> @polyvar x[1:2]
f=1+x[1]^4+x[2]^4-x[1]*x[2]
julia> opt,sol,basis=soncsocp(f, x)
```
By default, SONCSOCP uses Int64 as the integer type. You can modify it by setting ltype="Int128" (for solving linear equations) or itype="Int128" (for computing rational mediated sets) if an overflow occurs.
```Julia
julia> opt,data,basis=soncsocp(f, x, ltype="Int128", itype="Int128")
```

If the Newton polytope of f is a scaled standard simplex, you can set ntype="standard" so that the implementation will be more efficient.
```Julia
julia> opt,data,basis=soncsocp(f, x, ntype="standard")
```

## References
[1] [A second order cone characterization for sums of nonnegative circuits](https://arxiv.org/abs/1906.06179)  

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): jwang@laas.fr
