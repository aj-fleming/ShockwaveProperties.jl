"""
**Equation 4.8** from Anderson&Anderson

Computes the density change across a shock wave where free stream mach number is ``M``
and the outward (away from body) normal is ``n̂``
"""
function shock_density_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    Mn2 = (M * n̂[1])^2
    return ((gas.γ + 1) * Mn2) / ((gas.γ - 1) * Mn2 + 2)
end

"""
**Equation 4.9** from Anderson&Anderson

Computes the pressure ratio across a shock wave where the free stream mach number is ``M``
and the outward (away from body) normal is ``n̂``.
"""
function shock_pressure_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    Mn2 = (M * n̂[1])^2
    return 1 + (2 * gas.γ) / (gas.γ + 1) * (Mn2 - 1)
end

"""
**Equation 4.10** from Anderson&Anderson

Computes the normal Mach number ratio across a shock wave where the free stream mach number is ``M``
and the outward (away from body) normal is ``n̂``.
"""
function shock_normal_mach_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    MnL2 = (M * n̂[1])^2
    Mn2sqr = (MnL2 + (2 / (gas.γ - 1))) / (2 * gas.γ / (gas.γ - 1) * MnL2 - 1)
    Mn2 = sqrt(Mn2sqr)
    return Mn2
end


"""
Computes the normal velocity ratio across a shock wave where the free stream mach number is ``M``
and the outward (away from body) normal is ``n̂``.

Derived from speed of sound proportional to the square root of temperature.

Useful for computing momentum ratio.
"""
function shock_normal_velocity_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    return sqrt(1 / shock_temperature_ratio(M, n̂; gas=gas)) * shock_normal_mach_ratio(M, n̂; gas=gas)
end

"""
Computes the tangential mach number ratio across a shock wave, where the free stream mach number is ``M``
and the outward (away from body) normal is ``n̂``.

Derived from speed of sound proportional to the square root of temperature.
"""
function shock_tangent_mach_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    return sqrt(shock_temperature_ratio(M, n̂; gas=gas))
end

"""
**Equation 4.11** from Anderson&Anderson

Computes the temperature across a shock wave where the free stream mach number is ``M``
and the outward (away from body) normal is ``n̂``.
"""
function shock_temperature_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    return shock_pressure_ratio(M, n̂; gas=gas) / shock_density_ratio(M, n̂; gas=gas)
end

"""
Computes the 
"""