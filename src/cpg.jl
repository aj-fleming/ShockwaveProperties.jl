using Unitful
using LinearAlgebra

const _units_cvcp = u"J/kg/K"
const _units_ℳ = u"kg/mol"
const _units_ρ = u"kg/m^3"
const _units_v = u"m/s"
const _units_T = u"K"

const _units_ρv = _units_ρ * _units_v
const _units_int_e = _units_cvcp * _units_T
const _units_ρE = _units_ρ * _units_int_e

const _dimension_ρ = Unitful.dimension(_units_ρ)
const _dimension_v = Unitful.dimension(_units_v)
const _dimension_int_e = Unitful.dimension(_units_int_e)
const _dimension_ρE = Unitful.dimension(_units_ρE)
const _dimension_ρv = Unitful.dimension(_units_ρv)

"""
Provides the properties of a calorically perfect gas (or mixture of gases). 
"""
struct CaloricallyPerfectGas
    c_p::Quantity{Float64,Unitful.dimension(_units_cvcp),typeof(_units_cvcp)}
    c_v::Quantity{Float64,Unitful.dimension(_units_cvcp),typeof(_units_cvcp)}
    ℳ::Quantity{Float64,Unitful.dimension(_units_ℳ),typeof(_units_ℳ)}

    γ::Float64
    R::Quantity{Float64,Unitful.dimension(_units_cvcp),typeof(_units_cvcp)}
end

function CaloricallyPerfectGas(c_p::Float64, c_v::Float64, ℳ::Float64)
    return CaloricallyPerfectGas(
        Quantity(c_p, _units_cvcp),
        Quantity(c_v, _units_cvcp),
        Quantity(ℳ, _units_ℳ),
        c_p / c_v,
        Quantity(c_p - c_v, _units_cvcp)
    )
end

enthalpy(gas::CaloricallyPerfectGas, T::Float64) = gas.c_p * Quantity(T, u"K")
enthalpy(gas::CaloricallyPerfectGas, T::Quantity{Float64,Unitful.𝚯,Units}) where {Units} = gas.c_p * T

internal_energy(gas::CaloricallyPerfectGas, T::Float64) = gas.c_v * Quantity(T, u"K")
internal_energy(gas::CaloricallyPerfectGas, T::Quantity{Float64,Unitful.𝚯,Units}) where {Units} = gas.c_v * T


function speed_of_sound(gas::CaloricallyPerfectGas, T::Float64)
    return sqrt(gas.γ * gas.R * Quantity(T, u"K"))
end

"""
Compute the speed of sound in an ideal gas at a temperature ``T``. We assume 
that the gas is a non-dispersive medium.
"""
function speed_of_sound(gas::CaloricallyPerfectGas, T::Quantity{Float64,Unitful.𝚯,Units}) where {Units}
    return sqrt(gas.γ * gas.R * T)
end

"""
Compute the pressure in a calorically perfect gas from its density and temperature.
"""
function pressure(gas::CaloricallyPerfectGas, ρ::Float64, T::Float64)
    # calculate ρe then 
    # calculate p from calorically perfect gas relations
    return (gas.γ - 1) * Quantity(ρ, _units_ρ) * internal_energy(gas, T)
end

function pressure(
    gas::CaloricallyPerfectGas,
    ρ::Quantity{Float64,_dimension_ρ,Units1},
    T::Quantity{Float64,Unitful.𝚯,Units2}) where {Units1,Units2}
    return (gas.γ - 1) * ρ * internal_energy(gas, T)
end

"""
Compute the pressure in a calorically perfect gas from its internal energy density.
"""
pressure(gas::CaloricallyPerfectGas, ρe::Float64) = (gas.γ - 1) * Quantity(ρe, _units_ρE)
pressure(gas::CaloricallyPerfectGas, ρe::Quantity{Float64,_dimension_ρE,Units}) where {Units} = (gas.γ - 1) * ρe

"""
The conserved quantities in the Euler equations.
"""
struct ConservedState
    ρ::Quantity{Float64,_dimension_ρ,_units_ρ}
    ρv::Vector{Quantity{Float64,_dimension_ρv,_units_ρv}}
    ρE::Quantity{Float64,_dimension_ρE,_units_ρE}
end

"""
Compute the internal energy volume density (ρe) from conserved state quantities.
"""
internal_energy_density(state::ConservedState) = state.ρE - (state.ρv ⋅ state.ρv) / (2 * state.ρ)

"""
Properties, that are easier to reason about than those in a `ConservedState`, 
that also completely determine the state of a calorically perfect gas.
"""
struct PrimitiveState
    ρ::Quantity{Float64,_dimension_ρ,_units_ρ}
    M::Vector{Float64}
    T::Quantity{Float64,Unitful.𝚯,_units_T}
end

"""
Compute the pressure at a given state in a gas.
"""
pressure(gas::CaloricallyPerfectGas, state::ConservedState) = pressure(gas, internal_energy_density(state))
pressure(gas::CaloricallyPerfectGas, state::PrimitiveState) = pressure(gas, state.ρ, state.T)

"""
Compute the temperature at a given state in a gas.
"""
temperature(gas::CaloricallyPerfectGas, state::ConservedState) = internal_energy_density(state) / gas.c_v
temperature(::CaloricallyPerfectGas, state::PrimitiveState) = state.T

"""
Compute the density at a given state in a gas.
"""
density(::CaloricallyPerfectGas, state::Union{ConservedState,PrimitiveState}) = state.ρ

"""
Compute the speed of sound in a gas at a given state. We assume that the gas is a non-dispersive medium.
"""
function speed_of_sound(gas::CaloricallyPerfectGas, state::Union{ConservedState,PrimitiveState})
    return speed_of_sound(gas, temperature(gas, state))
end