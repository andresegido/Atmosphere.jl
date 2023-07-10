module Atmosphere

# ISA constants
export P0, T0, ρ0, μ0, γ, Rg, cs, h_levels, a_levels, T_levels, P_levels, ρ_levels
# ISA functions
export level, T, P, ρ, μ, sound
# unit conversion functions
export meters, feet, mps, knots, celsius, kelvin

include("isa.jl")
include("units.jl")

end
