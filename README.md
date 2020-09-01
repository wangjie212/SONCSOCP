# SONCSOCP
SONCSOCP is an uncontrained sparse polynomial optimization tool based on SONC decompositions.

To use SONCSOCP, run
```Julia
pkg> add https://github.com/wangjie212/SONCSOCP
 ```

## Dependencies
- Julia
- MOSEK 8.1

SONCSOCP has been tested on WINDOW 10, Julia 1.2.0, and MOSEK 8.1.
## Usage
SONCSOCP computes a lower bound for the unconstrained polynomial optimization problem:
```
Inf{f(x): x\in R^n}
```
where f is a polynomial with variables x1,...,xn.

Taking f=x1^4+x2^4-x1\*x2 as an example, to compute a SONC lower bound for f, run
```Julia
julia> using TSSOS
julia> using TypedPolynomials
julia> @polyvar x[1:2]
f=x[1]^4+x[2]^4-x[1]*x[2]
julia> opt,sol=soncsocp(f,x)
```
By default, SONCSOCP uses Int64 as the integer type. You can change it by setting ltype="Int128" (for solving linear equations) or itype="Int128" (for computing mediated sets) if an overflow occurs.
```Julia
julia> opt,data=soncsocp(f,x,ltype="Int128",itype="Int128")
```

If the Newton polytope of f is a scaled standard simplex, you can set ntype="standard" so that the implementation will be more efficient.
```Julia
julia> opt,data=soncsocp(f,x,ntype="standard")
```

## References
[1] [A second order cone characterization for sums of nonnegative circuits](https://arxiv.org/abs/1906.06179)  

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): jwang@laas.fr
