# EMD.jl
Empirical Mode Decomposition in Julia

> [!CAUTION]
> This is a refactored fork of the orginal project [jarrison/EMD.jl](https://github.com/jarrison/EMD.jl)

# Usage
Install this way.
```julia
using Pkg
Pkg.add(url="https://github.com/shayandavoodii/EMD.jl")
```

Set-up a signal to be decomposed using the EMD.
```julia
julia> using EMD, Plots
julia> t = LinRange(0,1.0,10^3)
julia> s = 3*sin.(2π*8*t) + sin.(2π*4*t)
```

Apply the decomposition with default values.
```julia
julia> imfs = sEMD(s)
julia> plot(imfs)
```

Typically we want fewer IMFs, the number of sifts will greatly affect the orthogonality of the resulting IMFs.
```julia
julia> imfs = sEMD(s,maximfs=5, nsifts=2)
julia> plot(imfs)
```
# References
## sfEMD
Order-Statistics envelopes for extracting IMFs.  
(1) Li, H.; Hu, Y.; Li, F.; Meng, G. Succinct and Fast Empirical Mode Decomposition. Mech. Syst. Signal Process. 2017, 85 (August 2015), 879–895.

## bEMD
B-Spline midpoint interpolation method.  
(2) Li, H.; Wang, C.; Zhao, D. An Improved EMD and Its Applications to Find the Basis Functions of EMI Signals. 2015, 2015.
