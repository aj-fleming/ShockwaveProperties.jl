module ShockwaveProperties

# constants
export DRY_AIR
# submodules
export BilligShockParametrization
# gas properties
export ConservedState, PrimitiveState
export density, temperature, pressure, speed_of_sound
export enthalpy, internal_energy, internal_energy_density

# shock wave jumps
export shock_density_ratio, shock_pressure_ratio, shock_temperature_ratio
export state_behind

include("cpg.jl")
include("billig.jl")
include("normal_shocks.jl")

const DRY_AIR = CaloricallyPerfectGas(1.0049, 0.7178, 0.0289647)

end
