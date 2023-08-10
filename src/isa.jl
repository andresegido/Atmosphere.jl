# Atmosphere constants
const g = 9.80665
const P0 = 101325
const T0 = 288.15
const ρ0 = 1.225
const μ0 = 1.7894e-5
const γ = 1.4
const Rg = P0 / (ρ0 * T0)
const cs = 110

const h_levels = [0, 11, 20, 32, 47, 51, 71, 84.852] * 1e3
const a_levels = [-6.5, 0, 1, 2.8, 0, -2.8, -2] * 1e-3

const T_levels = Array{Real}(undef, size(a_levels, 1))
const P_levels = Array{Real}(undef, size(a_levels, 1))
const ρ_levels = Array{Real}(undef, size(a_levels, 1))
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
    hMin, hMax = -610, 84852
    h < hMin && throw(ArgumentError("Height $h set lower than ISA minimum $hMin"))
    h > hMax && throw(ArgumentError("Height $h set higher than ISA maximum $hMax"))
    return h ≤ 0 ? 1 : findfirst(x -> x ≥ h, h_levels) - 1
end

T(h::Real, ΔT::Real=0) = (i=level(h); T_levels[i] + a_levels[i] * (h-h_levels[i]) + ΔT)
# Loop required to initialize T_levels correctly
for i=eachindex(T_levels)[2:end]
    T_levels[i] = T(h_levels[i])
end

function P(h::Real)
    i = level(h)
    if a_levels[i] == 0
        return P_levels[i] * exp(-g*(h-h_levels[i])/(Rg*T_levels[i]))
    else
        return P_levels[i] * (T(h)/T_levels[i]) ^ (-g/(Rg*a_levels[i]))
    end
end
# Loop required to initialize P_levels correctly
for i=eachindex(P_levels)[2:end]
    P_levels[i] = P(h_levels[i])
end

function ρ(h::Real, ΔT::Real=0)
    if ΔT ≈ 0
        i = level(h)
        if a_levels[i] == 0
            return ρ_levels[i] * exp(-g*(h-h_levels[i])/(Rg*T_levels[i]))
        else
            return ρ_levels[i] * (T(h)/T_levels[i]) ^ (-g/(Rg*a_levels[i])-1)
        end
    else
        return ρ_gas(P(h), T(h, ΔT))
    end
end
# Loop required to initialize ρ_levels correctly
for i=eachindex(ρ_levels)[2:end]
    ρ_levels[i] = ρ(h_levels[i])
end

μ(h::Real, ΔT::Real=0) = μ_gas(T(h, ΔT))
sound(h::Real, ΔT::Real=0) = sound_gas(T(h, ΔT))