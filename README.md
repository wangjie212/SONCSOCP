# SONCSOCP
SONCSOCP is an uncontrained sparse polynomial optimization tool based on SONC decompositions.

To use SONCSOCP, run
```Julia
pkg> add https://github.com/wangjie212/SONCSOCP
 ```

## Dependencies
- JuMP
- MOSEK
- ECOS

SONCSOCP has been tested on WINDOW 10, Julia 1.6.2, and MOSEK 9.0.
## Usage
SONCSOCP computes a lower bound for the unconstrained polynomial optimization problem:
$$\text{Inf}\ \lbrace f(\mathbf{x}): \mathbf{x}\in\mathbb{R}^n\rbrace,$$
where $f$ is a polynomial in variables $x_1,\ldots,x_n$.

Taking $f=1+x_1^4+x_2^4-x_1x_2$ as an example, to compute a SONC lower bound for $f$, run
```Julia
using SONCSOCP
using DynamicPolynomials
@polyvar x[1:2]
f = 1 + x[1]^4 + x[2]^4 - x[1]*x[2]
opt = soncsocp(f, x)
```

If the Newton polytope of $f$ is a scaled standard simplex (resp. simplex), then you may set ntype="standard" (resp. ntype="simplex") so that the computation will be more efficient.
```Julia
opt = soncsocp(f, x, ntype="standard")
```

## References
[1] [A second order cone characterization for sums of nonnegative circuits](https://arxiv.org/abs/1906.06179)  
[2] [SONC Optimization and Exact Nonnegativity Certificates via Second-Order Cone Programming](https://arxiv.org/abs/2012.07903)  

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn
