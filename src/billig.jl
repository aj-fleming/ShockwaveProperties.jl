"""
Implementation of the shockwave parametrization about cylindrical blunt bodies developed by Billig in
*Shock-wave shapes around spherical-and cylindrical-nosed bodies.*
"""
module BilligShockParametrization
using LinearAlgebra

export shock_front, shock_normal, shock_tangent

const BILLIG_CONSTANTS = (A=0.386, B=4.67, C=1.386, D=1.8)

# compute standoff distance relative to a disk of radius 1 in the plane
@inline function standoff_distance(M_inf)
    return BILLIG_CONSTANTS.A * exp(BILLIG_CONSTANTS.B / M_inf^2)
end

# compute critical radius relative to a disk of radius 1 in the plane
@inline function critical_radius(M_inf)
    return BILLIG_CONSTANTS.C * exp(BILLIG_CONSTANTS.D / (M_inf - 1)^0.75)
end

# compute sine and cosine of the critical angle 
# formed between the x-axis and the tangent to the shock front at x=0
@inline function critical_shock_angle(M_inf)
    sinθ = 1 / M_inf
    cosθ = sqrt(1 - sinθ^2)
    return (sinθ, cosθ)
end

"""
    shock_front(α, M_inf, R_b)
Maps a parameter ``α`` to the shock wave generated about a body with radius ``R_b`` 
by free-stream flow to the right with mach number ``M_infty``.

The shock wave parametrization developed by *Billig* parametrizes the shock wave on the y-coordinate.
"""
function shock_front(α, M_inf, R_b)
    y̅ = α / R_b
    Δ = standoff_distance(M_inf)
    R̅c = critical_radius(M_inf)
    sinθ, cosθ = critical_shock_angle(M_inf)
    cotθ = (cosθ / sinθ)
    tanθ = (sinθ / cosθ)
    return [R_b * (-1 - Δ + R̅c * cotθ^2 * (sqrt(1 + (y̅ * tanθ / R̅c)^2) - 1)), α]
end

# we'll manually take derivatives here... Zygote doesn't like 2nd derivatives.
# we can also optimize out cot^2 * tan^2 ;)

"""
    shock_normal(α, M_inf, R_b)
Maps a parameter ``α`` to the outward-facing normal to a shock wave generated about a body with radius ``R_b``
by a free-stream flow to the right with nach number ``M_infty``

The shock wave parametrization developed by *Billig* parametrizes the shock wave on the y-coordinate.
"""
function shock_normal(α, M_inf, R_b)
    R_c = R_b * critical_radius(M_inf)
    sinθ, cosθ = critical_shock_angle(M_inf)
    tanθ = (sinθ / cosθ)
    dxdα = α / (R_c * sqrt(1 + (α * tanθ / R_c)^2))
    n = [-1, dxdα]
    return n ./ norm(n)
end

"""
    shock_normal(α, M_inf, R_b)
Maps a parameter ``α`` to the tangent to a shock wave generated about a body with radius ``R_b``
by a free-stream flow to the right with nach number ``M_infty``

The shock wave parametrization developed by *Billig* parametrizes the shock wave on the y-coordinate.
"""
function shock_tangent(α, M_inf, R_b)
    n = shock_normal(α, M_inf, R_b)
    return [n[2], -n[1]]
end

end # module