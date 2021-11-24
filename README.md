# InterRelax
InterRelax is a Julia package of solving polynomial optimization problem on the nonnegative orthant:
```
f* := inf_{x in R^n} {f(x) : x >= 0, gi(x) >= 0, hj(x) = 0}.
```
We provide a hierarchy of semidefinite relaxations based on Putinar-Valilescu's Positivstellensatz for solving this problem. It has the following advantages:
- The maxmal matrix size of each semidefinite relaxation can be arbitrarily chosen. Accordingly, each semidefinite relaxation of relative small maxmal matrix size can be solved efficiently by using interior point methods (Mosek, SDPT3, ...).
- We guarantee the convergence of the sequence of values returned by this hierarchy to ```f*``` under some mild condition. Moreover, the rate of convergence is at least ```O(eps^-c)```.


# Required softwares
InterRelax has been implemented on a desktop compute with the following softwares:
- Ubuntu 18.04.4
- Julia 1.3.1

Before installing InterRelax, you should install [TSSOS](https://github.com/wangjie212/TSSOS) and [PolyPowerModels](https://github.com/tweisser/PolyPowerModels) with the following commands:
```ruby
Pkg> add https://github.com/wangjie212/TSSOS
Pkg> add https://github.com/tweisser/PolyPowerModels.git
```

# Installation
- To use InterRelax in Julia, run
```ruby
Pkg> add https://github.com/maihoanganh/InterRelax.git
```

# Usage
The following examples briefly guide to use InterRelax:

## Polynomial optimization on the nonnegative orthant
Consider the following optimization polynomial problem:
```ruby
using DynamicPolynomials

@polyvar x[1:2] # nonnegative variables

f=x[1]^2+0.5*x[1]*x[2]-0.25*x[2]^2+0.75*x[1]-0.3*x[2] # the objective polynomial to minimize

g=[1.0-sum(x.^2)] # the inequality constraints
h=[(x[1]-1.0)*x[2]] # the equality constraints

k=1 # relaxation order
s=3 # sparsity order

using InterRelax

# get information from the input data f,gi,hj
n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=InterRelax.get_info(x,f,g,h,sparse=false);

# get an approximate optimal value and an approximate optimal solution of the polynomial optimization problem
    opt_val,opt_sol=InterRelax.RelaxDense(n,m,l,
                                          lmon_g,supp_g,coe_g, # information of the inequality constraints
                                          lmon_h,supp_h,coe_h, # information of the equality constraints
                                          lmon_f,supp_f,coe_f, # information of the objective polynomial
                                          dg,dh,k,s,
                                          solver="Mosek", # solver for the semidefinite program
                                          comp_opt_sol=true) # to get an approximate optimal solution
```

See other examples from .ipynb files in the [link](https://github.com/maihoanganh/InterRelax/tree/main/examples).


# References
For more details, please refer to:

**N. H. A. Mai, J.-B. Lasserre, V. Magron and K-C Toh. A hierarchy of interrupted relaxations forpolynomial optimization on the nonnegative orthant. 2020. Forthcoming.**

The following codes are to run the paper's benchmarks:
```ruby
using InterRelax

InterRelax.test_AMGM() #Table 1
InterRelax.test_dense_POP_arbcons() #Tables 3 and 11
InterRelax.test_CS_POP_arbcons() #Tables 4 and 12
InterRelax.test_PMSV() #Table 5
InterRelax.test_compute_stability_number_of_graph_random() #Table 6
InterRelax.test_compute_stability_number_of_graph() #Table 7
InterRelax.test_deciding_copositivity() #Table 8
InterRelax.test_deciding_nonegativity() #Table 9
InterRelax.test_dense_POP_binary_constr_random() #Table 10

```
