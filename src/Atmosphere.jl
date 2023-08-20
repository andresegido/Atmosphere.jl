module Atmosphere
import IdealGas as gas
# ISA constants
export P0, T0, ρ0, μ0, γ, Rg, cs, a0
# ISA functions
export level, temperature, pressure, density, viscosity, sound
export T, P, ρ, μ, a
# unit conversion functions
export ft2m, m2ft, knots2ms, ms2knots, kel2cel, cel2kel
# FlightLevel struct
export FlightLevel

include("isa.jl")
include("units.jl")
include("flight_levels.jl")

end
