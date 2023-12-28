module ShockwaveProperties

export DRY_AIR
export BilligShockParametrization
export shock_density_ratio, shock_pressure_ratio, shock_temperature_ratio
export conserved_state_behind

include("cpg.jl")
include("billig.jl")
include("normal_shocks.jl")

const DRY_AIR = CaloricallyPerfectGas(1.0049, 0.7178, 0.0289647)

end
