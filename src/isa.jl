# Atmosphere constants
"Gravity constant: `g = 9.80665 m/s2`"
const g = 9.80665
"ISA pressure at sea level: `P0 = 101325 Pa`"
const P0 = 101325
"ISA temperature at sea level: `T0 = 288.15 K (=15.15 C)`"
const T0 = 288.15
"ISA density at sea level: `ρ0 = 1.225 kg/m3`"
const ρ0 = 1.225
"ISA viscosity at sea level: `μ0 = 1.7894e-5 Pa·s`"
const μ0 = 1.7894e-5
"Air's heat capacity ratio: `γ = 1.4`"
const γ = 1.4
"Air's specific gas constant: `Rg = P0 / (ρ0 * T0)`"
const Rg = P0 / (ρ0 * T0)
"Sutherland's constant for air at sea level: `cs = 110 K`"
const cs = 110

const h_lim = [-610, 84852]
const h_levels = [0, 11, 20, 32, 47, 51, 71, 84.852] * 1e3
const grad_levels = [-6.5, 0, 1, 2.8, 0, -2.8, -2] * 1e-3

const T_levels = Array{Real}(undef, size(grad_levels, 1))
const P_levels = Array{Real}(undef, size(grad_levels, 1))
const ρ_levels = Array{Real}(undef, size(grad_levels, 1))
T_levels[1] = T0
P_levels[1] = P0
ρ_levels[1] = ρ0

T_gas(P::Real, ρ::Real) = gas.T(P, ρ, Rg)
P_gas(ρ::Real, T::Real) = gas.P(ρ, T, Rg)
ρ_gas(P::Real, T::Real) = gas.T(P, T, Rg)
μ_gas(T::Real) = gas.μ(μ0, T, T0, cs)
sound_gas(T::Real) = gas.sound(γ, Rg, T)

"""
    level(h::Real)
Return an integer (i) indicating the atmosphere layer given a specific
altitude in meters (h).

This function is used internally by other functions in this module
to determine its base-level characteristincs and compute from there.

The correspondance between the integer (i) and the layers is:
- i = 1: Troposphere
- i = 2: Tropopause
- i = 3: Low Stratosphere
- i = 4: High Stratosphere
- i = 5: Stratopause
- i = 6: Low Mesosphere
- i = 7: High Mesosphere
"""
function level(h::Real)
    h < h_lim[1] && throw(ArgumentError("Height $h set lower than ISA minimum $hMin"))
    h > h_lim[end] && throw(ArgumentError("Height $h set higher than ISA maximum $hMax"))
    return h ≤ 0 ? 1 : findfirst(x -> x ≥ h, h_levels) - 1
end

"""
    temperature(h::Real, ΔT::Real=0)
Return temperature in Kelvin given an altitude in meters (h) based on ISA model.
Temperature deviations from ISA may be introduced as an optional ΔT.

Optional alias: `T(h::Real, ΔT::Real=0)`
"""
temperature(h::Real, ΔT::Real=0) = (i=level(h); T_levels[i] + grad_levels[i] * (h-h_levels[i]) + ΔT)
const T = temperature
# Loop required to initialize T_levels correctly
for i=eachindex(T_levels)[2:end]
    T_levels[i] = T(h_levels[i])
end

"""
    pressure(h::Real)
Return pressure in Pascals given an altitude in meters (h) based on ISA model.

Optional alias: `P(h::Real)`
"""
function pressure(h::Real)
    i = level(h)
    if grad_levels[i] == 0
        return P_levels[i] * exp(-g*(h-h_levels[i])/(Rg*T_levels[i]))
    else
        return P_levels[i] * (T(h)/T_levels[i]) ^ (-g/(Rg*grad_levels[i]))
    end
end
const P = pressure
# Loop required to initialize P_levels correctly
for i=eachindex(P_levels)[2:end]
    P_levels[i] = P(h_levels[i])
end

"""
    density(h::Real, ΔT::Real=0)
Return density in meters per kilogram given an altitude in meters (h) based on ISA model.
Temperature deviations from ISA may be introduced as an optional ΔT.

Optional alias: `ρ(h::Real, ΔT::Real=0)`
"""
function density(h::Real, ΔT::Real=0)
    if ΔT == 0
        i = level(h)
        if grad_levels[i] == 0
            return ρ_levels[i] * exp(-g*(h-h_levels[i])/(Rg*T_levels[i]))
        else
            return ρ_levels[i] * (T(h)/T_levels[i]) ^ (-g/(Rg*grad_levels[i])-1)
        end
    else
        return ρ_gas(P(h), T(h, ΔT))
    end
end
const ρ = density
# Loop required to initialize ρ_levels correctly
for i=eachindex(ρ_levels)[2:end]
    ρ_levels[i] = ρ(h_levels[i])
end

"""
    viscosity(h::Real, ΔT::Real=0)
Return viscosity in Pascals·second given an altitude in meters (h) based on ISA model.
Temperature deviations from ISA may be introduced as an optional ΔT.

Optional alias: `μ(h::Real, ΔT::Real=0)`
"""
viscosity(h::Real, ΔT::Real=0) = μ_gas(T(h, ΔT))
const μ = viscosity
"""
    sound(h::Real, ΔT::Real=0)
Return speed of sound in meters per second (sound) given an altitude in meters (h) based on ISA model.
Temperature deviations from ISA may be introduced as an optional ΔT.

Optional alias: `a(h::Real, ΔT::Real=0)`
"""
sound(h::Real, ΔT::Real=0) = sound_gas(T(h, ΔT))
const a = sound

"ISA speed of sound at sea level"
const a0 = a(0)