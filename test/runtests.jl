using LinearAlgebra
using ShockwaveProperties
using Test

@testset verbose=true "ShockwaveProperties.jl" begin
    @testset "BilligShockParametrization" begin
        using .BilligShockParametrization
        @testset "Shockwave Normals" begin
            y = -10:0.1:10
            for R_b ∈ 1.0:0.5:4.0, M ∈ 2.0:1.0:5.0
                vals = mapreduce(hcat, y) do α
                    shock_normal(α, M, R_b) ⋅ shock_tangent(α, M, R_b)
                end
                @test all(≈(0), vals)
            end
        end
    end
end