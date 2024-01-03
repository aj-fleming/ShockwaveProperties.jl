# ShockwaveProperties.jl

[![Build Status](https://github.com/aj-fleming/ShockwaveProperties.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/aj-fleming/ShockwaveProperties.jl/actions/workflows/CI.yml?query=branch%3Amaster)

---

This package provides routines for computing property changes across shock waves in perfect gases. The ultimate goal is to support a solver that can be used to analyze complex interactions of shock waves in dynamic systems.

## 

The `PrimitiveState` and `ConservedState` types allow for conversion between specifying gas properties as $\rho, \vec M, T$ and $\rho, \rho \vec v,\rho E$.

```julia
free_stream = PrimitiveState(1.225u"kg/m^3", [2.0, 0.0], 300u"K")
u = ConservedState(free_stream, DRY_AIR)
```


