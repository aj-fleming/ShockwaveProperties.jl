using Unitful
using Unitful: Temperature, Density, Velocity, @derived_dimension
using Unitful: 𝐋, 𝐓, 𝐌, 𝚯, 𝐍
using LinearAlgebra

@derived_dimension HeatCapacity 𝐋^2 * 𝐓^-2 * 𝚯^-1
@derived_dimension MolarMass 𝐌 * 𝐍^-1
@derived_dimension MomentumDensity 𝐌 * 𝐋^-2 * 𝐓^-1
@derived_dimension SpecificEnergy 𝐋^2 * 𝐓^-2
@derived_dimension EnergyDensity 𝐌 * 𝐋^-1 * 𝐓^-2

# default units for certain quantities.
const _units_cvcp = u"J/kg/K"
const _units_ℳ = u"kg/mol"
const _units_ρ = u"kg/m^3"
const _units_v = u"m/s"
const _units_T = u"K"

const _units_ρv = _units_ρ * _units_v
const _units_int_e = _units_cvcp * _units_T
const _units_ρE = _units_ρ * _units_int_e

"""
    CaloricallyPerfectGas
Provides the properties of a calorically perfect gas (or mixture of gases).
 - ``c_p``: Heat capacity at constant pressure
 - ``c_v``: Heat capacity at constant volume
 - ``ℳ``: Molar mass
 - ``γ``: Heat capacity ratio
 - ``R``: Gas constant
"""
struct CaloricallyPerfectGas{U1<:HeatCapacity,U2<:MolarMass}
    c_p::U1
    c_v::U1
    ℳ::U2

    γ::Float64
    R::U1
end

function CaloricallyPerfectGas(c_p::Float64, c_v::Float64, ℳ::Float64)
    q_cp = Quantity(c_p, _units_cvcp)
    q_cv = Quantity(c_v, _units_cvcp)
    q_ℳ = Quantity(ℳ, _units_ℳ)
    q_R = q_cp - q_cv
    return CaloricallyPerfectGas{typeof(q_cp),typeof(q_ℳ)}(q_cp, q_cv, q_ℳ, c_p / c_v, q_R)
end

"""
    enthalpy(gas::CaloricallyPerfectGas, T)
Computes the enthalpy of a calorically perfect gas.
"""
enthalpy(gas::CaloricallyPerfectGas, T::Float64) = gas.c_p * Quantity(T, _units_T)
enthalpy(gas::CaloricallyPerfectGas, T::Temperature) = gas.c_p * T

"""
    internal_energy(gas::CaloricallyPerfectGas, T)
Computes the internal energy of a calorically perfect gas.
"""
internal_energy(gas::CaloricallyPerfectGas, T::Float64) = gas.c_v * Quantity(T, _units_T)
internal_energy(gas::CaloricallyPerfectGas, T::Temperature) = gas.c_v * T

"""
    speed_of_sound(gas::CaloricallyPerfectGas, T)
Computes the speed of sound in an ideal gas at a temperature ``T``. We assume 
that the gas is a non-dispersive medium.
"""
speed_of_sound(gas::CaloricallyPerfectGas, T::Float64) = sqrt(gas.γ * gas.R * Quantity(T, _units_T))
speed_of_sound(gas::CaloricallyPerfectGas, T::Temperature) = sqrt(gas.γ * gas.R * T)


"""
    pressure(gas::CaloricallyPerfectGas, ρ, T)
Compute the pressure in a calorically perfect gas from its density and temperature.
"""
function pressure(gas::CaloricallyPerfectGas, ρ::Float64, T::Float64)
    return (gas.γ - 1) * Quantity(ρ, _units_ρ) * internal_energy(gas, T)
end

function pressure(gas::CaloricallyPerfectGas, ρ::Density, T::Temperature)
    return (gas.γ - 1) * ρ * internal_energy(gas, T)
end

"""
    pressure(gas::CaloricallyPerfectGas, ρE)
Compute the pressure in a calorically perfect gas from its internal energy density.
"""
pressure(gas::CaloricallyPerfectGas, ρe::Float64) = (gas.γ - 1) * Quantity(ρe, _units_ρE)
pressure(gas::CaloricallyPerfectGas, ρe::EnergyDensity) = (gas.γ - 1) * ρe

"""
Properties, that are easier to reason about than those in a `ConservedState`. 
These completely determine the state of a calorically perfect gas.

 - ``ρ``: Density of the gas.
 - ``M``: The mach number, represented as a vector quantity.
 - ``T``: The absolute temperature of the gas.
"""
struct PrimitiveState{U1<:Density,U2<:Temperature}
    ρ::U1
    M::Vector{Float64}
    T::U2
end

"""
    PrimitiveState(ρ::Float64, M::Vector{Float64}, T::Float64)
Construct a PrimitiveState and assign the default units.
"""
function PrimitiveState(ρ::Float64, M::Vector{Float64}, T::Float64)
    return PrimitiveState(Quantity(ρ, _units_ρ), M, Quantity(T, _units_T))
end

"""
    ConservedState{<:Density, <:MomentumDensity, <:EnergyDensity}
The conserved quantities in the Euler equations.

 - ``ρ``: Density of the gas.
 - ``ρv``: Momentum density of the gas, represented as a vector quantity.
 - ``ρE``: The **sum** of the internal energy density of the gas and kinetic energy density of the moving gas.
"""
struct ConservedState{U1<:Density,U2<:MomentumDensity,U3<:EnergyDensity}
    ρ::U1
    ρv::Vector{U2}
    ρE::U3
end

"""
    PrimitiveState(state::ConservedState, gas::CaloricallyPerfectGas)
Compute the primitive state quantities from the conserved state quantities and 
the thermodynamic properies of the gas.
"""
function PrimitiveState(state::ConservedState, gas::CaloricallyPerfectGas)
    a = speed_of_sound(gas, state)
    return PrimitiveState(state.ρ, state.ρv / (state.ρ * a), temperature(gas, state))
end

"""
    ConservedState(ρ::Float64, ρv::Vector{Float64}, ρE::Float64)
Construct a ConservedState and assign the default units.
"""
function ConservedState(ρ::Float64, ρv::Vector{Float64}, ρE::Float64)
    return ConservedState(Quantity(ρ, _units_ρ), Quantity.(ρv, _units_ρv), Quantity(ρE, _units_ρE))
end

"""
    ConservedState(state::PrimitiveState, gas::CaloricallyPerfectGas)
Compute the conserved state quantities from the primitive state quantities and
the thermodynamic properties of the gas.
"""
function ConservedState(state::PrimitiveState, gas::CaloricallyPerfectGas)
    v = state.M * speed_of_sound(gas, state)
    e = gas.c_v * state.T
    return ConservedState(state.ρ, state.ρ * v, state.ρ * (e + v ⋅ v / 2))
end

"""
    internal_energy_density(state::ConservedState)
Compute the internal energy volume density (``ρe``) from conserved state quantities.
"""
internal_energy_density(state::ConservedState) = state.ρE - (state.ρv ⋅ state.ρv) / (2 * state.ρ)

"""
    internal_energy(gas::CaloricallyPerfectGas, state)
Compute the internal energy (``e``) of a gas at a given state.
"""
internal_energy(::CaloricallyPerfectGas, state::ConservedState) = internal_energy_density(state) / state.ρ
internal_energy(gas::CaloricallyPerfectGas, state::PrimitiveState) = gas.c_v * state.T 

"""
    pressure(gas::CaloricallyPerfectGas, state)
Compute the pressure at a given state in a gas.
"""
pressure(gas::CaloricallyPerfectGas, state::ConservedState) = pressure(gas, internal_energy_density(state))
pressure(gas::CaloricallyPerfectGas, state::PrimitiveState) = pressure(gas, state.ρ, state.T)

"""
    density(gas::CaloricallyPerfectGas, state)
Compute the density at a given state in a gas.
"""
density(::CaloricallyPerfectGas, state::Union{ConservedState,PrimitiveState}) = state.ρ

"""
    temperature(gas::CaloricallyPerfectGas, state)
Compute the temperature at a given state in a gas.
"""
function temperature(gas::CaloricallyPerfectGas, state::ConservedState)
    return internal_energy(gas, state) / gas.c_v
end

temperature(::CaloricallyPerfectGas, state::PrimitiveState) = state.T



"""
    speed_of_sound(gas::CaloricallyPerfectGas, state::Union{ConservedState, PrimitiveState})
Compute the speed of sound in a gas at a given state. We assume that the gas is a non-dispersive medium.
"""
function speed_of_sound(gas::CaloricallyPerfectGas, state)
    return speed_of_sound(gas, temperature(gas, state))
end