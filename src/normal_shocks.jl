using LinearAlgebra

"""
**Equation 4.8** from Anderson&Anderson

Computes the density change across a shock wave. 
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.
"""
function shock_density_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    Mn2 = (M ⋅ n̂)^2
    return ((gas.γ + 1) * Mn2) / ((gas.γ - 1) * Mn2 + 2)
end

"""
**Equation 4.9** from Anderson&Anderson

Computes the pressure ratio across a shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.
"""
function shock_pressure_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    Mn2 = (M ⋅ n̂)^2
    return 1 + (2 * gas.γ) / (gas.γ + 1) * (Mn2 - 1)
end

"""
**Equation 4.10** from Anderson&Anderson

Computes the normal Mach number ratio across a shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.
"""
function shock_normal_mach_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    Mnsqr = (M ⋅ n̂)^2
    Mn2sqr = (Mnsqr + (2 / (gas.γ - 1))) / ((2 * gas.γ / (gas.γ - 1)) * Mnsqr - 1)
    mach_ratio = sqrt(Mn2sqr/Mnsqr)
    return mach_ratio
end

"""
Computes the normal velocity ratio across a shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.

Derived from speed of sound proportional to the square root of temperature.
Useful for computing momentum ratio.
"""
@inline function shock_normal_velocity_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    return sqrt(shock_temperature_ratio(M, n̂; gas=gas)) * shock_normal_mach_ratio(M, n̂; gas=gas)
end

"""
Computes the tangential mach number ratio across a shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.

Derived from speed of sound proportional to the square root of temperature.
"""
@inline function shock_tangent_mach_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    return sqrt(1 / shock_temperature_ratio(M, n̂; gas=gas))
end

@inline function shock_normal_momentum_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    return shock_density_ratio(M, n̂; gas=gas) * shock_normal_velocity_ratio(M, n̂; gas=gas)
end

"""
**Equation 4.11** from Anderson&Anderson

Computes the temperature ratio across a shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.
"""
@inline function shock_temperature_ratio(M, n̂; gas::CaloricallyPerfectGas=DRY_AIR)
    return shock_pressure_ratio(M, n̂; gas=gas) / shock_density_ratio(M, n̂; gas=gas)
end

"""
    state_behind(state_L, n̂, t̂; gas::CaloricallyPerfectGas)
Computes the gas state behind a shockwave.
The outward (away from body) normal to the shockwave is ``n̂`` and the tangent to the shockwave is ``t̂``.
"""
function state_behind(state_L::ConservedState, n̂, t̂; gas::CaloricallyPerfectGas=DRY_AIR)
    @assert t̂ ⋅ n̂ == 0.0 "tangent and normal vectors should be normal."
    M_L = state_L.ρv / (state_L.ρ * speed_of_sound(state_L; gas=gas))
    # density change
    ρ_R = state_L.ρ * shock_density_ratio(M_L, n̂; gas=gas)
    # momentum change
    ρv_n_R = (state_L.ρv ⋅ n̂) * shock_normal_momentum_ratio(M_L, n̂; gas=gas)
    ρv_t_R = (state_L.ρv ⋅ t̂) * shock_density_ratio(M_L, n̂; gas=gas)
    ρv_R = ρv_n_R * n̂ + ρv_t_R * t̂
    # total internal energy change
    ρe_L = internal_energy_density(state_L)
    # e = cT
    ρe_R = ρe_L * shock_density_ratio(M_L, n̂; gas=gas) * shock_temperature_ratio(M_L, n̂; gas=gas)
    ρE_R = ρe_R + (ρv_R ⋅ ρv_R) / (2 * ρ_R)
    return ConservedState(ρ_R, ρv_R, ρE_R)
end

function state_behind(state_L::PrimitiveState, n̂, t̂; gas::CaloricallyPerfectGas=DRY_AIR)
    @assert t̂ ⋅ n̂ == 0.0 "tangent and normal vectors should be normal."
    M_n_L = state_L.M ⋅ n̂
    M_t_L = state_L.M ⋅ t̂
    M_n_R = M_n_L * shock_normal_mach_ratio(state_L.M, n̂; gas=gas)
    M_t_R = M_t_L * shock_tangent_mach_ratio(state_L.M, n̂; gas=gas)
    M_R = M_n_R * n̂ + M_t_R * t̂
    ρ_R = state_L.ρ * shock_density_ratio(state_L.M, n̂; gas=gas)
    T_R = state_L.T * shock_temperature_ratio(state_L.M, n̂; gas=gas)
    return PrimitiveState(ρ_R, M_R, T_R)
end

