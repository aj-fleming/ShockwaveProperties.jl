using LinearAlgebra

"""
**Equation 4.8** from Anderson&Anderson

Computes the density change across a shock wave. 
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.
"""
@inline function shock_density_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    Mn2 = (M ⋅ n̂)^2
    return ((gas.γ + 1) * Mn2) / ((gas.γ - 1) * Mn2 + 2)
end

"""
**Equation 4.9** from Anderson&Anderson

Computes the pressure ratio across a shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.
"""
@inline function shock_pressure_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    Mn2 = (M ⋅ n̂)^2
    return 1 + (2 * gas.γ) / (gas.γ + 1) * (Mn2 - 1)
end

"""
**Equation 4.10** from Anderson&Anderson

Computes the normal Mach number ratio across a shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.
"""
function shock_normal_mach_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    MnL2 = (M ⋅ n̂)^2
    Mn2sqr = (MnL2 + (2 / (gas.γ - 1))) / (2 * gas.γ / (gas.γ - 1) * MnL2 - 1)
    Mn2 = sqrt(Mn2sqr)
    return Mn2
end

"""
Computes the normal velocity ratio across a shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.

Derived from speed of sound proportional to the square root of temperature.
Useful for computing momentum ratio.
"""
@inline function shock_normal_velocity_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    return sqrt(1 / shock_temperature_ratio(M, n̂; gas=gas)) * shock_normal_mach_ratio(M, n̂; gas=gas)
end

"""
Computes the tangential mach number ratio across a shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.

Derived from speed of sound proportional to the square root of temperature.
"""
function shock_tangent_mach_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    return sqrt(shock_temperature_ratio(M, n̂; gas=gas))
end

@inline function shock_normal_momentum_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    return shock_density_ratio(M, n̂; gas=gas) * shock_normal_velocity_ratio(M, n̂; gas=gas)
end

"""
**Equation 4.11** from Anderson&Anderson

Computes the temperature ratio across a shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.
"""
function shock_temperature_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    return shock_pressure_ratio(M, n̂; gas=gas) / shock_density_ratio(M, n̂; gas=gas)
end

"""
Computes the conservation properties behind a shockwave.
The outward (away from body) normal is ``n̂``.
"""
function conserved_state_behind(state_L::ConservedState, n̂, t̂; gas::CaloricallyPerfectGas=DRY_AIR)
    @assert t̂ ⋅ n̂ == 0. "tangent and normal vectors should be normal."
    M_L = state_L.ρv / (state_L.ρ * speed_of_sound(gas, state_L))
    ρv_n_L = (state_L.ρv ⋅ n̂) * shock_normal_momentum_ratio(M_L, n̂; gas=gas)
    ρv_t_L = (state_L.ρv ⋅ n̂) * shock_density_ratio(M_L, n̂; gas=gas)
    
end