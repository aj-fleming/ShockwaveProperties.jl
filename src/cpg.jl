using Unitful

const _units_cvcp = u"J/kg/K"
const _units_ℳ = u"kg/mol"
const _units_ρ = u"kg/m^3"
const _dimension_ρ = Unitful.dimension(_units_ρ)

"""
Provides the properties of a calorically perfect gas (or mixture of gases). 
"""
struct CaloricallyPerfectGas{T}
    c_p::Quantity{T,Unitful.dimension(_units_cvcp),typeof(_units_cvcp)}
    c_v::Quantity{T,Unitful.dimension(_units_cvcp),typeof(_units_cvcp)}
    ℳ::Quantity{T,Unitful.dimension(_units_ℳ),typeof(_units_ℳ)}

    γ::T
    R::Quantity{T,Unitful.dimension(_units_cvcp),typeof(_units_cvcp)}
end

function CaloricallyPerfectGas(c_p::T, c_v::T, ℳ::T) where {T}
    return CaloricallyPerfectGas{T}(
        Quantity(c_p, _units_cvcp),
        Quantity(c_v, _units_cvcp),
        Quantity(ℳ, _units_ℳ),
        c_p / c_v,
        Quantity(c_p - c_v, _units_cvcp)
    )
end

const DRY_AIR = CaloricallyPerfectGas(1.0049, 0.7178, 0.0289647)

enthalpy(gas::CaloricallyPerfectGas, T::Real) = gas.c_p * Quantity(T, u"K")
enthalpy(gas::CaloricallyPerfectGas, T::Quantity{T1,Unitful.𝚯,Units}) where {T1,Units} = gas.c_p * T

int_energy(gas::CaloricallyPerfectGas, T::Real) = gas.c_v * Quantity(T, u"K")
int_energy(gas::CaloricallyPerfectGas, T::Quantity{T1,Unitful.𝚯,Units}) where {T1,Units} = gas.c_v * T

function speed_of_sound(gas::CaloricallyPerfectGas, T::Real)
    return sqrt(gas.γ * gas.R * Quantity(T, u"K"))
end

function speed_of_sound(gas::CaloricallyPerfectGas, T::Quantity{T1,Unitful.𝚯,Units}) where {T1,Units}
    return sqrt(gas.γ * gas.R * T)
end

function pressure(gas::CaloricallyPerfectGas, ρ::Real, T::Real)
    # calculate ρe then 
    # calculate p from calorically perfect gas relations
    return (gas.γ - 1) * Quantity(ρ, _units_ρ) * int_energy(gas, T)
end

function pressure(gas::CaloricallyPerfectGas,
    ρ::Quanity{U,_dimension_ρ,Units1},
    T::Quantity{U,Unitful.𝚯,Units2}) where {U,Units1,Units2}
    return (gas.γ-1) * ρ * T
end
